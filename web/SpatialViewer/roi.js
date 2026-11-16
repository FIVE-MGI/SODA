// GeoMx DSP ROI layers (prepare.py geomx -> geomx.py; an expression layer of kind "rois").
// A row of such a layer is one ROI segment (AOI): a hand-drawn region with ONE count vector,
// a few dozen per study. viewer.js treats it as any layer (gene files, groupings, plots,
// volcano); this file holds what differs:
//   - the words (ROIs and counts instead of bins and UMI),
//   - the drawing: each ROI's own outline filled with the colour of the genes on show (a disc
//     per segment where an ROI has several), with a hover readout; no density grid, as ROIs
//     do not tile the slide,
//   - the notes that say why the UMAP view and Circle -> UMAP's own UMAP are not offered,
//   - coloured by a grouping (an annotation column): each ROI's outline filled with its group's colour,
//   - Q3 normalisation: geomx.py writes GeomxTools' Q3 normFactors per row into the layer's "q3"
//     entry at export; one choice (Q3, the default, or counts per 10,000) shared by the gene colours,
//     the plots (plots.js) and the volcano (volcano.js). An older export has no factors: nothing changes.
// Only unsigned typed arrays here: the de-identification check reads the signed names as IDs.
/* exported SodaRoi */
(function () {
  // ---- Q3 normalisation --------------------------------------------------------------------------
  // value = count / normFactor (GeomxTools normalize(norm_method = "quant")): the browser multiplies
  // by 1 / normFactor; a row without a factor (Q3 = 0) gets 0 = left out of the plots and volcano
  let norm = "q3", rerender = null, scaleMemo = { f: null, s: null };
  const q3Of = (e) => (e && e.kind === "rois" && e.q3 && Array.isArray(e.q3.normFactors) && e.q3.normFactors.length === e.nBins ? e.q3 : null);
  function scale(e) {
    const q = norm === "q3" ? q3Of(e) : null;
    if (!q) return null;
    if (scaleMemo.f !== q.normFactors) scaleMemo = { f: q.normFactors, s: Float32Array.from(q.normFactors, (x) => (x > 0 ? 1 / x : 0)) };
    return scaleMemo.s;
  }
  function setNorm(v) {
    norm = v === "cp10k" ? "cp10k" : "q3";
    for (const s of document.querySelectorAll("select[data-roinorm]")) s.value = norm;
    if (rerender) rerender();
  }
  // the same menu in the plots panel and the volcano block; hidden until a layer has factors
  const normRow = (id) =>
    `<label class="row slim" id="${id}Row" hidden><span title="Q3: each ROI's counts divided by its normalisation factor (its 75th-percentile count over the targets above LOQ in at least 10% of the ROIs, over the geometric mean of all ROIs' Q3), as GeomxTools and the GeoMx software do by default. Counts per 10,000: each ROI's counts over its total. Also sets the colours on the slide (Q3 or raw counts).">Normalise</span>` +
    `<select id="${id}" data-roinorm="1"><option value="q3">Q3 (GeoMx standard)</option><option value="cp10k">counts per 10,000</option></select></label>`;
  function bindNorm(id) { const el = document.getElementById(id); if (el) { el.value = norm; el.onchange = () => setNorm(el.value); } }
  const normShown = (id, e) => { const r = document.getElementById(`${id}Row`); if (r) r.hidden = !q3Of(e); };
  const normLabel = (s) => (s ? "Q3-normalised counts" : "counts per 10,000");

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
    // Q3 on: the colours and the readout are Q3-normalised counts (count-like, so the slider windows still read in counts)
    const s = scale(a.expr) && scale(a.expr).length === N ? scale(a.expr) : null;
    const rows = [];
    for (let r = 0; r < N; r++) {
      const x = pos[2 * r], y = pos[2 * r + 1];
      if (!(x === x)) continue;
      let sum = 0, top = 0, cr = 0, cg = 0, cb = 0;
      const lines = [];
      on.forEach((g, i) => {
        const q = s && s[r] > 0, v = q ? vals[i][r] * s[r] : vals[i][r];
        const lo = g.min || 0, t = Math.max(0, Math.min(1, (v - lo) / Math.max(g.max - lo, 1e-6)));
        sum += t; if (t > top) top = t;
        const c = a.hexToRgb255(g.hex);
        cr += c[0] * t; cg += c[1] * t; cb += c[2] * t;
        lines.push(`${g.name}: ${fmt(v)} ${q ? "Q3-normalised counts" : s ? `${a.unit} (no Q3 factor: Q3 = 0)` : a.unit}`);
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
    return { rows, sig: on.map((g) => `${g.name}:${g.hex}:${g.min}:${g.max}`).join("|") + `#${a.palette}#${s ? "q3" : "raw"}` };
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
    if (a.render) rerender = a.render;
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

  // Coloured by a grouping (an annotation column): an ROI holding one segment is its own outline
  // filled with its group's colour; segments that share an ROI stay discs. Colours, the lit / faded
  // state and the hidden groups all come from the viewer's cloud C (rows of positions.bin, as
  // region.bin). hooks.discs draws a cloud as the viewer does; hooks.register(id, {n, rows}) lets the
  // viewer's hover and click find a fill's row as they find a disc's (layer ids hold "ucloud").
  function groupLayers(a, C, hooks) {
    if (a.render) rerender = a.render;
    const map = shapesOf(a);
    if (!map || !map.size || !C.n || !C.rows) return [hooks.discs(C)];
    const pos = C.data.attributes.getPosition.value, col = C.data.attributes.getFillColor.value;
    const fill = [], rest = [];
    for (let j = 0; j < C.n; j++) {
      const poly = map.get(a.regionIds[C.rows[j]]);
      if (poly) fill.push({ j, polygon: poly, color: [col[4 * j], col[4 * j + 1], col[4 * j + 2], col[4 * j + 3]] }); else rest.push(j);
    }
    const out = [];
    if (fill.length) {
      const id = `roi-ucloud-fill-${a.dsId}-${a.exprKey}`;
      hooks.register(id, { n: fill.length, rows: Uint32Array.from(fill, (f) => C.rows[f.j]) });
      out.push(new a.PolygonLayer({
        id, data: fill, getPolygon: (d) => d.polygon, getFillColor: (d) => d.color,
        // the fill lies over the outline overlay: a thin edge of its own keeps neighbouring ROIs apart
        stroked: true, getLineColor: [255, 255, 255, 150], lineWidthUnits: "pixels", getLineWidth: 1,
        filled: true, pickable: true, parameters: a.NORMAL, updateTriggers: { getFillColor: [C.key] },
      }));
    }
    if (rest.length) {
      const n = rest.length, p = new Float32Array(2 * n), c = new Uint8Array(4 * n);
      rest.forEach((j, k) => { p[2 * k] = pos[2 * j]; p[2 * k + 1] = pos[2 * j + 1]; for (let t = 0; t < 4; t++) c[4 * k + t] = col[4 * j + t]; });
      out.push(hooks.discs({ key: `${C.key}|roi-discs`, n, rows: Uint32Array.from(rest, (j) => C.rows[j]),
        data: { length: n, attributes: { getPosition: { value: p, size: 2 }, getFillColor: { value: c, size: 4, normalized: true } } } }));
    }
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
    const s = scale(a.expr), unit = s ? "Q3-normalised counts" : a.unit;
    const how = s
      ? ` Q3-normalised counts: each ROI's counts divided by its Q3 factor (GeomxTools' Q3 normalisation; "Normalise" in Heatmap and violin or Volcano). The plots and the volcano use the same.`
      : ` Counts are not normalised here: a larger or deeper-sequenced ROI has more of everything. The plots and the volcano use counts per 10,000.`;
    if (a.palette === "genes") {
      return `<div class="genes">${on.map((g) => `<span class="g"><i style="background:${g.hex}"></i>${esc(g.name)}</span>`).join("")}</div>` +
        `<div class="note">One filled outline per ROI (a disc per segment where an ROI has several). Grey = at or below the left handle of the gene's slider, full colour = at the right handle ` +
        `(${on.map((g) => `${esc(g.name)} ${win(g)}`).join("; ")} ${unit}).` + (on.length > 1 ? " Where genes overlap the colours add." : "") + how + `</div>`;
    }
    return `<div class="ramp"><span>low</span><i style="background:${a.rampCss}"></i><span>high</span></div>` +
      `<div class="note">One filled outline per ROI (a disc per segment where an ROI has several): ` +
      (on.length > 1 ? `mean of ${on.length} genes, each on its own slider window (${on.map((g) => `${esc(g.name)} ${win(g)}`).join("; ")} ${unit})`
                     : `${esc(on[0].name)}, ${win(on[0])} ${unit} (the slider's handles)`) + "." + how + `</div>`;
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
      `${hvg && hvg.nGenes ? hvg.nGenes.toLocaleString() : "all"} targets (every target with a count), on ${q3Of(e) ? "Q3-normalised counts or counts per 10,000 (its Normalise menu)" : "counts per 10,000"}` +
      (hvg && hvg.rounded ? `; the multi-probe targets' geometric means are rounded to whole counts for it` : "") +
      (hvg && hvg.clipped ? `; ${hvg.clipped.toLocaleString()} values above 65,535 are capped there` : "") +
      `. With this few ROIs, and ROIs of one slide or animal not being independent, read the fold changes; the p-values are rough.`;
  }

  window.SodaRoi = { WORDS, toWords, geneLayers, groupLayers, tip, keyHtml, umapNote, syncLasso,
    // Q3 normalisation, shared by the gene colours, plots.js and volcano.js
    q3: q3Of, scale, setNorm, normNow: () => norm, normRow, bindNorm, normShown, normLabel };
})();
