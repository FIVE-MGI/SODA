// GeoMx DSP ROI layers (prepare.py geomx -> geomx.py; an expression layer of kind "rois").
// A row of such a layer is one ROI segment (AOI): a hand-drawn region with ONE count vector,
// a few dozen per study. viewer.js treats it as any layer (gene files, groupings, plots,
// volcano); this file holds what differs:
//   - the words (ROIs and counts instead of bins and UMI),
//   - the drawing: each ROI's own outline filled with the colour of the genes on show (a disc
//     per segment where an ROI has several), with a hover readout; no density grid, as ROIs
//     do not tile the slide,
//   - the notes that say why the UMAP view and Circle -> UMAP's own UMAP are not offered.
// Only unsigned typed arrays here: the de-identification check reads the signed names as IDs.
/* exported SodaRoi */
(function () {
  const WORDS = { bin: "ROI", bins: "ROIs", Bins: "ROIs", gene: "gene", genes: "genes", Gene: "Gene", Genes: "Genes" };
  const toWords = (t) => String(t)
    .replace("Bins are neighbours, not independent samples: the p-values describe this section only.",
             "ROIs of one slide or animal are not independent samples: the p-values describe these ROIs only.")
    .replace(/of the section's ([\d,]+) variable genes/, "of the $1 targets with a count")
    .replace(/^Rest of the section$/, "All the other ROIs").replace(/ of the section \(/, " (")
    .replace(/\bbins\b/g, "ROIs").replace(/\bbin\b/g, "ROI")
    .replace(/\bBins\b/g, "ROIs").replace(/\bBin\b/g, "ROI").replace(/\bUMI\b/g, "counts");
  const esc = (x) => String(x).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
  const GREY = [70, 76, 90];
  const fmt = (v) => (!v ? "0" : Number.isInteger(v) ? v.toLocaleString("en-US") : v >= 100 ? v.toFixed(0) : v.toFixed(1));
  const key2 = (x, y) => Math.round(x * 100) * 1e8 + Math.round(y * 100);
  let rowMap = { pos: null, map: null };

  // gene files hold x, y, value of the rows with a count: a row is found again by its position
  function byRow(g, pos) {
    if (g.roiRow && g.roiRowPos === pos) return g.roiRow;
    if (rowMap.pos !== pos) {
      const m = new Map();
      for (let r = 0; r < pos.length / 2; r++) m.set(key2(pos[2 * r], pos[2 * r + 1]), r);
      rowMap = { pos, map: m };
    }
    const out = new Float32Array(pos.length / 2);
    for (let i = 0; i < g.n; i++) {
      const r = rowMap.map.get(key2(g.positions[2 * i], g.positions[2 * i + 1]));
      if (r !== undefined) out[r] = g.values[i];
    }
    g.roiRow = out; g.roiRowPos = pos;
    return out;
  }

  // What is drawn for the genes on show.
  //   An ROI with one segment: its own outline, filled (the shape that was measured).
  //   Segments that share an ROI, or a layer whose outlines are not loaded: one disc per segment.
  // A ramp palette: the mean of the genes' positions in their windows (left handle = bottom of
  // the ramp, right handle = top). "Gene colours": each gene's colour weighted by that position,
  // added up. Every segment is drawn: low is the bottom of the ramp (or grey), never missing,
  // because an ROI that was measured and is low is a result.
  function rowsOf(a) {
    const on = a.genes.filter((g) => g.visible !== false);
    // the layer's positions; without them (the live tile server has no positions file in memory
    // until asked) the first gene's own rows stand in, and segments are drawn as discs
    const pos = a.pos || (on[0] && on[0].positions);
    if (!pos || !on.length) return null;
    const N = pos.length / 2, vals = on.map((g) => byRow(g, pos));
    const ramp = a.palette !== "genes" ? (a.palettes[a.palette] || a.palettes.viridis) : null;
    const names = a.rowName || (() => "");
    const rows = [];
    for (let r = 0; r < N; r++) {
      const x = pos[2 * r], y = pos[2 * r + 1];
      if (!(x === x)) continue;
      let sum = 0, top = 0, cr = 0, cg = 0, cb = 0;
      const lines = [];
      on.forEach((g, i) => {
        const lo = g.min || 0, t = Math.max(0, Math.min(1, (vals[i][r] - lo) / Math.max(g.max - lo, 1e-6)));
        sum += t; if (t > top) top = t;
        const c = a.hexToRgb255(g.hex);
        cr += c[0] * t; cg += c[1] * t; cb += c[2] * t;
        lines.push(`${g.name}: ${fmt(vals[i][r])} ${a.unit}`);
      });
      let color;
      if (ramp) {
        const k = Math.min(255, Math.round((sum / on.length) * 255)) * 3;
        color = [ramp[k], ramp[k + 1], ramp[k + 2], 240];
      } else {
        const m = Math.max(cr, cg, cb, 1), s = top > 0 ? Math.min(1, 255 / m) : 0;
        color = top > 0 ? [Math.min(255, cr * s * (0.35 + 0.65 * top) + GREY[0] * (1 - top)),
                           Math.min(255, cg * s * (0.35 + 0.65 * top) + GREY[1] * (1 - top)),
                           Math.min(255, cb * s * (0.35 + 0.65 * top) + GREY[2] * (1 - top)), 240] : [...GREY, 200];
      }
      const nm = names(r);
      rows.push({ row: r, position: [x, y], color, lines, properties: { name: (nm ? nm + "\n" : "") + lines.join("\n") } });
    }
    return { rows, sig: on.map((g) => `${g.name}:${g.hex}:${g.min}:${g.max}`).join("|") + `#${a.palette}` };
  }

  // region id -> the outline's rings, for ROIs that hold exactly one segment
  let shapes = { src: null, ids: null, map: null };
  function shapesOf(a) {
    if (!a.outlines || !a.regionIds) return null;
    if (shapes.src === a.outlines && shapes.ids === a.regionIds) return shapes.map;
    const n = new Map(), map = new Map();
    for (const id of a.regionIds) if (id) n.set(id, (n.get(id) || 0) + 1);
    for (const f of a.outlines.features || []) {
      const id = f.properties && f.properties.__region_id, g = f.geometry;
      if (id && n.get(id) === 1 && g && g.type === "Polygon") map.set(id, g.coordinates);
    }
    shapes = { src: a.outlines, ids: a.regionIds, map };
    return map;
  }

  function geneLayers(a) {
    const R = rowsOf(a);
    if (!R) return [];
    const map = a.pos ? shapesOf(a) : null, filled = [], discs = [];
    for (const d of R.rows) {
      const poly = map && map.get(a.regionIds[d.row]);
      if (poly) filled.push({ ...d, polygon: poly }); else discs.push(d);
    }
    const out = [];
    if (filled.length) out.push(new a.PolygonLayer({
      id: `roi-fill-${a.dsId}-${a.exprKey}`, data: filled,
      getPolygon: (d) => d.polygon, getFillColor: (d) => d.color,
      stroked: false, filled: true, opacity: a.opacity, pickable: true, parameters: a.NORMAL,
      updateTriggers: { getFillColor: [R.sig] },
    }));
    if (discs.length) out.push(new a.ScatterplotLayer({
      id: `roi-discs-${a.dsId}-${a.exprKey}`, data: discs,
      getPosition: (d) => d.position, getFillColor: (d) => d.color,
      radiusUnits: "common", getRadius: a.radius, radiusMinPixels: 3,
      stroked: true, getLineColor: [255, 255, 255, 200], lineWidthUnits: "pixels", getLineWidth: 1,
      filled: true, opacity: a.opacity, pickable: true, parameters: a.NORMAL,
      updateTriggers: { getFillColor: [R.sig] },
    }));
    return out;
  }

  // hovering an ROI outline (the overlay lies on top of the fills): its name, and the counts of
  // the genes on show for each of its segments
  function tip(a, regionId, name) {
    const R = rowsOf(a);
    if (!R || !a.regionIds || !a.pos) return name;
    const mine = R.rows.filter((d) => a.regionIds[d.row] === regionId);
    if (!mine.length) return name;
    return [name].concat(mine.length === 1 ? mine[0].lines
      : mine.map((d, i) => `segment ${i + 1}: ${d.lines.join(", ")}`)).join("\n");
  }

  // the colour key under the picture
  function keyHtml(a) {
    const on = a.genes.filter((g) => g.visible !== false);
    if (!on.length) return "";
    const win = (g) => `${fmt(g.min || 0)}–${fmt(g.max)}`;
    if (a.palette === "genes") {
      return `<div class="genes">${on.map((g) => `<span class="g"><i style="background:${g.hex}"></i>${esc(g.name)}</span>`).join("")}</div>` +
        `<div class="note">One filled outline per ROI (a disc per segment where an ROI has several). Grey = at or below the left handle of the gene's slider, full colour = at the right handle ` +
        `(${on.map((g) => `${esc(g.name)} ${win(g)}`).join("; ")} ${a.unit}).` + (on.length > 1 ? " Where genes overlap the colours add." : "") +
        ` Counts are not normalised here: a larger or deeper-sequenced ROI has more of everything.</div>`;
    }
    return `<div class="ramp"><span>low</span><i style="background:${a.rampCss}"></i><span>high</span></div>` +
      `<div class="note">One filled outline per ROI (a disc per segment where an ROI has several): ` +
      (on.length > 1 ? `mean of ${on.length} genes, each on its own slider window (${on.map((g) => `${esc(g.name)} ${win(g)}`).join("; ")} ${a.unit})`
                     : `${esc(on[0].name)}, ${win(on[0])} ${a.unit} (the slider's handles)`) +
      `. Counts are not normalised here: a larger or deeper-sequenced ROI has more of everything. The plots and the volcano use counts per 10,000.</div>`;
  }

  const umapNote = (e, n) =>
    `${(e.nBins || 0).toLocaleString()} ROI segments, each one measurement. ` +
    `The groups are the columns of the study's annotation sheet${n ? ` (${n})` : ""}. ` +
    `No UMAP or Leiden clusters: a few dozen hand-drawn regions are not a sample of the tissue, so an embedding of them would mostly redraw the annotation.`;

  // Circle -> UMAP panel on an ROI layer: drawing shapes stays (a shape is a group of ROIs for the
  // volcano and the plots); the panel's own UMAP run and its settings are hidden
  function syncLasso(on, e, hvg) {
    const $ = (id) => document.getElementById(id);
    for (const id of ["lassoRun", "lassoSaveCsv"]) if ($(id)) $(id).hidden = on;
    for (const id of ["lassoMode", "lassoSeed", "lassoRes"]) { const el = $(id) && $(id).closest("label"); if (el) el.hidden = on; }
    // layout.js wraps the heading's words in a button: the title is that button's first text
    const h = document.querySelector("#panelLasso h2 .sec-name") || document.querySelector("#panelLasso h2");
    if (h && h.firstChild && h.firstChild.nodeType === 3) {
      if (!h.dataset.orig) h.dataset.orig = h.firstChild.textContent;
      h.firstChild.textContent = on ? "Groups drawn on the slide " : h.dataset.orig;
    }
    if (on && $("lassoNote")) $("lassoNote").textContent =
      `Draw a circle or a loop around ROIs to make a group of them; two groups (drawn, or from the annotation) can be compared in the ` +
      `Volcano plot and Heatmap and violin sections. No UMAP is run on ROIs. The volcano tests ` +
      `${hvg && hvg.nGenes ? hvg.nGenes.toLocaleString() : "all"} targets (every target with a count), on counts per 10,000` +
      (hvg && hvg.rounded ? `; the multi-probe targets' geometric means are rounded to whole counts for it` : "") +
      (hvg && hvg.clipped ? `; ${hvg.clipped.toLocaleString()} values above 65,535 are capped there` : "") +
      `. With this few ROIs, and ROIs of one slide or animal not being independent, read the fold changes; the p-values are rough.`;
  }

  window.SodaRoi = { WORDS, toWords, geneLayers, tip, keyHtml, umapNote, syncLasso };
})();
