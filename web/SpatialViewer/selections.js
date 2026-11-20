// Export shapes: the Draw a region shapes as GeoJSON out and in, and a CSV of the bins (cells) inside the shapes
// with the instrument's own barcodes / cell ids.
//
// The shapes live in viewer.js's state (S.lasso.shapes). This file reads them, adds new ones and asks viewer.js to
// assign the bins again (api.changed = lassoShapesChanged); the volcano, the plots and Run UMAP then see a new shape
// exactly as they see a drawn one.
//
// Membership. viewer.js gives each bin to ONE shape, the one drawn last (L.assign). The exports here use the full
// membership instead: a bin is inside a shape when it is inside its outline (Circle, Freehand, an imported polygon)
// or in its list of bins (a selection), whatever was drawn after it. An imported polygon with holes, or a combined
// shape from a file or link made before 0.10.0 (Combine / Overlap / Subtract were removed then), keeps how it was made
// (s.def), so on another layer its bins are found again from the outlines (see saveDef below).
//
// GeoJSON. Coordinates are in the export's full-resolution image pixels (x right, y down, origin top-left), as every
// other file of the export. An outline is a Polygon; a list of bins is a MultiPoint of the bins' centres. Import takes
// Polygon / MultiPolygon features (holes kept), and the MultiPoints written here. A QuPath export lines up only when
// QuPath opened the same full-resolution image the export was made from (same crop, same flip, same scale).
//
// Barcodes. An export holds them by default (prepare.py export; --no-ids leaves them out, `prepare.py ids` adds them to
// an export made without): expr/<key>/ids.txt.gz, one id per positions.bin row, named in the manifest's layer entry
// (`ids`). They are the instrument's barcodes / cell ids, checked for sample IDs before writing (layerids.py): a
// layer whose ids failed the check has none.
/* exported SodaShapes */
(function () {
  "use strict";
  let api = null, sig = "";
  const picked = new Set();                     // ids of the ticked shapes
  const idCache = new Map();                    // layer key -> Promise<string[]>
  const $ = (id) => document.getElementById(id);
  const MAX_NAME = 60, MAX_IMPORT = 100;
  const EMPTY = new Uint32Array(0);
  const esc = (x) => String(x).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
  const W = () => (api && api.words ? api.words() : { bin: "bin", bins: "bins" });
  const cut = (t) => (t.length > MAX_NAME ? t.slice(0, MAX_NAME - 1) + "…" : t);
  const slug = (s) => String(s || "").normalize("NFKD").replace(/[^A-Za-z0-9.-]+/g, "_").replace(/^[_.]+|[_.]+$/g, "").slice(0, 40) || "shapes";
  const q = (x) => (/[",\n\r]/.test(String(x)) ? `"${String(x).replace(/"/g, '""')}"` : String(x));
  const L = () => { const S = api && api.state(); return S && S.lasso; };
  const msg = (t, err) => { const el = $("shMsg"); if (!el) return; el.textContent = t || ""; el.style.color = err ? "#ff8f8f" : ""; };

  // ---- the section, in the Analyze tab under Draw a region -------------------------------------------------
  function buildDom() {
    const sec = document.createElement("section");
    sec.className = "panel"; sec.id = "panelShapes"; sec.dataset.tab = "analyze"; sec.hidden = true;
    sec.innerHTML = `<h2>Export shapes</h2>
      <div id="shNote" class="meta cs-note"></div>
      <div id="shPick" class="cs-pick" role="group" aria-label="Shapes to export"></div>
      <div class="btnrow">
        <button id="shGeo" type="button" disabled title="The shapes (the ticked ones, else all) as GeoJSON in image pixels: for QuPath, Python or a later visit">Save GeoJSON</button>
        <button id="shGeoBtn" type="button" title="Add shapes from a GeoJSON file of polygons in this image's full-resolution pixels (SODA, QuPath)">Load GeoJSON</button>
        <input id="shGeoFile" type="file" accept=".geojson,.json,application/geo+json,application/json" hidden />
        <button id="shIds" type="button" disabled title="One row per bin inside the shapes (the ticked ones, else all): its barcode or cell id, position and shape">Save barcodes (CSV)</button>
      </div>
      <div id="shMsg" class="meta cs-msg" role="status"></div>`;
    const side = $("sidebar"), lasso = $("panelLasso");
    if (window.SodaLayout) { side.appendChild(sec); window.SodaLayout.adopt(sec); }
    else if (lasso) lasso.after(sec); else side.appendChild(sec);
    // shown exactly when Draw a region is
    if (lasso) {
      const f = () => { if (sec.hidden !== lasso.hidden) sec.hidden = lasso.hidden; };
      new MutationObserver(f).observe(lasso, { attributes: true, attributeFilter: ["hidden"] });
      f();
    }
    $("shGeo").onclick = saveGeoJson;
    $("shGeoBtn").onclick = () => $("shGeoFile").click();
    $("shGeoFile").onchange = (e) => { const f = e.target.files[0]; e.target.value = ""; if (f) loadGeoJson(f); };
    $("shIds").onclick = saveIds;
    $("shPick").addEventListener("change", (e) => {
      const id = Number(e.target && e.target.dataset ? e.target.dataset.s : NaN);
      if (!Number.isInteger(id)) return;
      if (e.target.checked) picked.add(id); else picked.delete(id);
      e.target.closest(".cs-chip").classList.toggle("on", e.target.checked);
      syncButtons();
    });
  }

  // ---- membership: the bins inside a shape, whatever was drawn after it --------------------------------------
  function bbox(ring) {
    let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity;
    for (const [x, y] of ring) { if (x < x0) x0 = x; if (y < y0) y0 = y; if (x > x1) x1 = x; if (y > y1) y1 = y; }
    return [x0, y0, x1, y1];
  }
  // polygons: [[outer ring, hole, ...], ...]; inside one outer ring and none of its holes
  function polysMask(polys, pos) {
    const n = pos.length / 2, m = new Uint8Array(n);
    for (const rings of polys) {
      if (!rings.length || rings[0].length < 3) continue;
      const [x0, y0, x1, y1] = bbox(rings[0]);
      for (let i = 0; i < n; i++) {
        if (m[i]) continue;
        const x = pos[2 * i], y = pos[2 * i + 1];
        if (!(x >= x0 && x <= x1 && y >= y0 && y <= y1)) continue;
        if (!api.inPoly(x, y, rings[0])) continue;
        let hole = false;
        for (let k = 1; k < rings.length && !hole; k++) hole = rings[k].length >= 3 && api.inPoly(x, y, rings[k]);
        if (!hole) m[i] = 1;
      }
    }
    return m;
  }
  function rowsMask(rows, n) { const m = new Uint8Array(n); for (const r of rows) if (r < n) m[r] = 1; return m; }
  const maskRows = (m) => { let k = 0; for (let i = 0; i < m.length; i++) if (m[i]) k++; const out = new Uint32Array(k); k = 0; for (let i = 0; i < m.length; i++) if (m[i]) out[k++] = i; return out; };
  // how a shape is made, independent of the layer where possible:
  //   { pts } an outline | { geo } polygons with holes | { rows, layer } bins of one layer | { op, parts } an operation
  //   (op: a combined shape from a file or link made before 0.10.0)
  // the mask of a definition on this layer's bins; null when part of it is bins of another layer
  function maskOfDef(d, pos, key) {
    const n = pos.length / 2;
    if (d.pts) return polysMask([[d.pts]], pos);
    if (d.geo) return polysMask(d.geo, pos);
    if (d.op) {
      const ms = d.parts.map((p) => maskOfDef(p, pos, key));
      if (ms.some((x) => !x)) return null;
      return combineMasks(d.op, ms, n);
    }
    if (d.rows) return d.layer === key ? rowsMask(d.rows, n) : null;
    return null;
  }
  function combineMasks(op, ms, n) {
    const out = new Uint8Array(n);
    for (let i = 0; i < n; i++) {
      if (op === "union") { for (const m of ms) if (m[i]) { out[i] = 1; break; } }
      else if (op === "inter") { let all = 1; for (const m of ms) if (!m[i]) { all = 0; break; } out[i] = all; }
      else { if (!ms[0][i]) continue; let any = 0; for (let k = 1; k < ms.length; k++) if (ms[k][i]) { any = 1; break; } out[i] = any ? 0 : 1; }
    }
    return out;
  }
  // the bins of one shape on the layer on show (its own list, or its outline)
  function maskOfShape(s, pos) {
    const n = pos.length / 2;
    if (s.rows) return rowsMask(s.rows, n);
    return s.pts && s.pts.length >= 3 ? polysMask([[s.pts]], pos) : new Uint8Array(n);
  }

  // viewer.js (lassoAssign), before the bins are given to the shapes: a list of bins belongs to one layer. A shape
  // made here is found again on another layer from its outlines; a plain selection keeps its bins for its own layer
  // and is empty on the others.
  function remap(shapes, pos, key) {
    if (!pos || !key) return;
    for (const s of shapes) {
      if (!s.rows) continue;
      if (!s.layer) { s.layer = key; s.home = { layer: key, rows: s.rows }; }
      if (!s.home) s.home = { layer: s.layer, rows: s.rows };
      if (s.layer === key) continue;
      let rows = null;
      if (s.def) { const m = maskOfDef(s.def, pos, key); if (m) rows = maskRows(m); }
      if (!rows) rows = s.home.layer === key ? s.home.rows : EMPTY;
      if (s.home.layer === PENDING) s.home = { layer: key, rows };      // restored from its recipe: this is its first layer
      s.rows = rows; s.layer = key;
      s.away = rows === EMPTY ? s.home.layer : null;
    }
  }

  // ---- how a shape was made, in Save JSON and in a link to the view (viewer.js sessionState / sessionApply) ----------
  // A shape made here from outlines (Circle, Freehand, imported polygons) is saved with its recipe as `def`:
  //   { points: [[x, y], ...] } | { geo: [polygon: [ring: [[x, y], ...], ...], ...] } | { op: "union"|"inter"|"sub", parts: [def, ...] }
  // in image pixels (2 decimals). Save JSON keeps its list of bins too (`rows`, as before); a link carries only the
  // recipe, which is far shorter, and the bins are found again from it when the link is opened. A shape with a plain
  // selection in it (bins of one layer, no outline) has no recipe and is saved as its list of bins only.
  const PENDING = "\u0000def";                  // layer of a shape restored from its recipe, before its bins are found
  const DEF_DEPTH = 6, DEF_LEAVES = 64, DEF_POINTS = 20000;
  const r2 = (v) => Math.round(v * 100) / 100;
  function saveDef(s) {
    if (!s || !s.def) return null;
    let leaves = 0;
    const pack = (d, depth) => {
      if (!d || depth > DEF_DEPTH) return null;
      if (d.pts) { leaves++; return { points: d.pts.map((p) => [r2(p[0]), r2(p[1])]) }; }
      if (d.geo) { leaves++; return { geo: d.geo.map((rings) => rings.map((ring) => ring.map((p) => [r2(p[0]), r2(p[1])]))) }; }
      if (d.op && OPS[d.op] && Array.isArray(d.parts)) {
        const parts = d.parts.map((x) => pack(x, depth + 1));
        return parts.every(Boolean) ? { op: d.op, parts } : null;
      }
      return null;                              // bins of one layer: no recipe
    };
    const out = pack(s.def, 0);
    return out && leaves <= DEF_LEAVES ? out : null;
  }
  // a recipe read from a file or a link (untrusted): checked, or null
  function readDef(j) {
    let leaves = 0, points = 0;
    const pt = (p) => Array.isArray(p) && p.length >= 2 && Number.isFinite(p[0]) && Number.isFinite(p[1]) && Math.abs(p[0]) < 1e9 && Math.abs(p[1]) < 1e9;
    const ring = (r) => { if (!Array.isArray(r) || r.length < 3 || !r.every(pt)) return null; points += r.length; return r.map((p) => [Number(p[0]), Number(p[1])]); };
    const read = (d, depth) => {
      if (!d || typeof d !== "object" || depth > DEF_DEPTH || ++leaves > DEF_LEAVES * 2) return null;
      if (Array.isArray(d.points)) { const r = ring(d.points); return r ? { pts: r } : null; }
      if (Array.isArray(d.geo)) {
        const geo = d.geo.slice(0, DEF_LEAVES).map((rings) => (Array.isArray(rings) && rings.length ? rings.slice(0, 64).map(ring) : null));
        return geo.length && geo.every((rs) => rs && rs[0]) ? { geo: geo.map((rs) => rs.filter(Boolean)) } : null;
      }
      if (typeof d.op === "string" && OPS[d.op] && Array.isArray(d.parts) && d.parts.length >= 2 && d.parts.length <= DEF_LEAVES) {
        const parts = d.parts.map((x) => read(x, depth + 1));
        return parts.every(Boolean) ? { op: d.op, parts } : null;
      }
      return null;
    };
    const d = read(j, 0);
    return d && points <= DEF_POINTS ? d : null;
  }
  // viewer.js sessionApply, right after it rebuilt the shapes from `saved` (the file's or the link's shapes, same order)
  function restore(shapes, saved) {
    if (!Array.isArray(shapes) || !Array.isArray(saved)) return;
    shapes.forEach((sh, i) => {
      const j = saved[i], d = j && j.def ? readDef(j.def) : null;
      if (d) {
        sh.def = d;
        if (!sh.rows) { sh.rows = EMPTY; sh.layer = PENDING; sh.home = null; }    // a link: the bins come from the recipe
      }
      if (sh.rows && !(sh.pts && sh.pts.length >= 3)) sh.kind = "rows";
    });
  }

  // the operations a recipe may hold (combined shapes made before 0.10.0: read, never made here any more)
  const OPS = { union: true, inter: true, sub: true };
  function tickedShapes() { const l = L(); return l ? l.shapes.filter((s) => picked.has(s.id)) : []; }
  // ---- what the section shows ---------------------------------------------------------------------------------
  function refresh() {
    const sec = $("panelShapes"), l = L();
    if (!api || !sec || !l) return;
    const S = api.state(), e = api.layer();
    const key = [S.dsId, S.exprKey, l.shapes.map((s) => `${s.id}:${s.bins}:${s.name}:${s.hex}`).join(","), e && e.ids ? 1 : 0].join("|");
    if (key === sig) { syncButtons(); return; }
    sig = key;
    const live = new Set(l.shapes.map((s) => s.id));
    for (const id of [...picked]) if (!live.has(id)) picked.delete(id);
    const w = W();
    const away = l.shapes.filter((s) => s.away);
    $("shNote").textContent = !l.shapes.length
      ? "Draw shapes in Draw a region (or load a GeoJSON file) to save them here."
      : `Tick shapes to save only those (none ticked: all of them).` +
        (away.length ? ` ${away.map((s) => `“${s.name}”`).join(", ")} ${away.length === 1 ? "was" : "were"} made on another layer (${[...new Set(away.map((s) => s.away))].join(", ")}): ${away.length === 1 ? "its" : "their"} ${w.bins} are not on this one.` : "");
    $("shPick").innerHTML = l.shapes.map((s) => {
      const on = picked.has(s.id);
      return `<label class="cs-chip${on ? " on" : ""}" title="${esc(s.name)}"><input type="checkbox" data-s="${s.id}"${on ? " checked" : ""} />` +
        `<i style="background:${esc(s.hex)}"></i><span>${esc(s.name)}</span></label>`;
    }).join("");
    $("shPick").hidden = !l.shapes.length;
    const t = $("shIds");
    t.title = e && e.ids
      ? `One row per ${w.bin} inside the shapes (the ticked ones, else all): its ${e.ids.column === "cell_id" ? "cell id" : "barcode"}, position and shape`
      : "This export has no barcode list for this layer (made with --no-ids, before lists were written by default, or its ids failed the sample ID check): run prepare.py ids on it";
    syncButtons();
  }
  function syncButtons() {
    const l = L(), any = !!(l && l.shapes.length && l.pos);
    $("shGeo").disabled = !any;
    $("shIds").disabled = !any;
  }

  // ---- GeoJSON out ------------------------------------------------------------------------------------------------
  const r3 = (v) => Math.round(v * 1000) / 1000;
  const hexRgb = (h) => { const m = /^#?([0-9a-f]{6})$/i.exec(String(h || "")); const v = m ? parseInt(m[1], 16) : 0x888888; return [(v >> 16) & 255, (v >> 8) & 255, v & 255]; };
  const closed = (ring) => { const r = ring.map(([x, y]) => [r3(x), r3(y)]); if (r.length && (r[0][0] !== r[r.length - 1][0] || r[0][1] !== r[r.length - 1][1])) r.push(r[0].slice()); return r; };
  function scope() { const t = tickedShapes(), l = L(); return t.length ? t : l ? l.shapes.slice() : []; }
  function geoJson() {
    const S = api.state(), l = L(), um = api.um(), key = S.exprKey;
    const features = [], skipped = [];
    for (const s of scope()) {
      const props = { objectType: "annotation", name: s.name, shape_type: s.kind === "circle" ? "circle" : s.kind === "free" ? "freehand" : s.def && s.def.geo ? "polygon" : "bins",
                      microns_per_pixel: um || null, dataset: S.dsId, color: hexRgb(s.hex) };
      let geometry;
      if (s.def && s.def.geo) geometry = { type: "MultiPolygon", coordinates: s.def.geo.map((rings) => rings.map(closed)) };
      else if (!s.rows && s.pts && s.pts.length >= 3) geometry = { type: "Polygon", coordinates: [closed(s.pts)] };
      else {
        if (!s.rows || !s.rows.length || !l.pos) { skipped.push(s.name); continue; }
        // a list of bins: their centres, exactly as positions.bin holds them (so they are found again)
        const pts = new Array(s.rows.length);
        for (let k = 0; k < s.rows.length; k++) { const r = s.rows[k]; pts[k] = [l.pos[2 * r], l.pos[2 * r + 1]]; }
        geometry = { type: "MultiPoint", coordinates: pts };
        props.layer = key; props.n_bins = s.rows.length;
        if (s.def && s.def.op) props.made_from = s.name;
      }
      features.push({ type: "Feature", geometry, properties: props });
    }
    return { gj: { type: "FeatureCollection",
      soda: { dataset: S.dsId, layer: key, microns_per_pixel: um || null,
              coordinates: "full-resolution image pixels of this export: x to the right, y down, origin at the top-left corner" },
      features }, skipped };
  }
  function saveGeoJson() {
    const S = api.state();
    const { gj, skipped } = geoJson();
    if (!gj.features.length) { msg("No shape to save.", true); return; }
    const name = `${slug(S.dsId)}_${slug(S.exprKey)}_shapes.geojson`;
    api.download(name, JSON.stringify(gj), "application/geo+json");
    msg(`Saved ${name}: ${gj.features.length} ${gj.features.length === 1 ? "shape" : "shapes"} in image pixels` +
        (api.um() ? ` (${api.um()} µm per pixel)` : "") + "." +
        (skipped.length ? ` Left out, no ${W().bins} on this layer: ${skipped.join(", ")}.` : "") +
        (api.touch ? " On a phone or tablet it goes to Downloads (iPhone and iPad: the Files app)." : ""));
  }

  // ---- GeoJSON in -------------------------------------------------------------------------------------------------
  const isXY = (p) => Array.isArray(p) && p.length >= 2 && Number.isFinite(+p[0]) && Number.isFinite(+p[1]);
  const ringOf = (r) => {
    if (!Array.isArray(r) || !r.every(isXY)) return null;
    const out = r.map((p) => [+p[0], +p[1]]);
    if (out.length > 1 && out[0][0] === out[out.length - 1][0] && out[0][1] === out[out.length - 1][1]) out.pop();
    return out.length >= 3 ? out : null;
  };
  const polyOf = (c) => { if (!Array.isArray(c)) return null; const rings = c.map(ringOf); return rings[0] ? rings.filter(Boolean) : null; };
  const colorOf = (c) => {
    if (typeof c === "string" && /^#[0-9a-f]{6}$/i.test(c)) return c;
    if (Array.isArray(c) && c.length >= 3 && c.every((v) => Number.isFinite(+v))) return "#" + c.slice(0, 3).map((v) => Math.max(0, Math.min(255, Math.round(+v))).toString(16).padStart(2, "0")).join("");
    if (Number.isInteger(c)) return "#" + ((c >>> 0) & 0xffffff).toString(16).padStart(6, "0");      // QuPath's packed RGB
    return null;
  };
  // positions.bin rows by place, for MultiPoints: nearest bin within 0.05 px
  function finder(pos) {
    const grid = new Map(), n = pos.length / 2;
    for (let i = 0; i < n; i++) { const k = `${Math.floor(pos[2 * i])},${Math.floor(pos[2 * i + 1])}`; const a = grid.get(k); if (a) a.push(i); else grid.set(k, [i]); }
    return (x, y) => {
      let best = -1, bd = 0.05 * 0.05;
      const fx = Math.floor(x), fy = Math.floor(y);
      for (let dx = -1; dx <= 1; dx++) for (let dy = -1; dy <= 1; dy++) {
        for (const i of grid.get(`${fx + dx},${fy + dy}`) || []) {
          const d = (pos[2 * i] - x) ** 2 + (pos[2 * i + 1] - y) ** 2;
          if (d <= bd) { bd = d; best = i; }
        }
      }
      return best;
    };
  }
  async function loadGeoJson(file) {
    let j;
    try { j = JSON.parse(await file.text()); } catch (e) { msg("Not a GeoJSON file (it does not read as JSON).", true); return; }
    return importGeoJson(j);
  }
  async function importGeoJson(j) {
    const feats = !j || typeof j !== "object" ? [] : j.type === "FeatureCollection" && Array.isArray(j.features) ? j.features
      : j.type === "Feature" ? [j] : j.type ? [{ type: "Feature", geometry: j, properties: {} }] : [];
    if (!feats.length) { msg("No features in that file: it needs a GeoJSON FeatureCollection of Polygon or MultiPolygon features.", true); return; }
    // a list of bins belongs to its layer: open that layer first, as Load JSON does
    const S0 = api.state(), want = feats.find((f) => f && f.geometry && f.geometry.type === "MultiPoint" && f.properties && f.properties.layer);
    if (want && want.properties.layer !== S0.exprKey && api.layers().includes(want.properties.layer)) {
      if (!(await api.setLayer(want.properties.layer))) return;
    }
    const l = L(), S = api.state();
    if (!l || !l.pos) { msg("The shapes are not ready yet: open a layer that Draw a region works on.", true); return; }
    const pos = l.pos, find = hasMultiPoint(feats) ? finder(pos) : null;
    const img = api.image(), Wd = (img && img.width) || 0, Ht = (img && img.height) || 0;
    const add = [], skip = { detections: 0, other: 0, unmatched: 0 };
    let outside = 0, total = 0, otherDs = null, otherUm = null;
    for (const f of feats) {
      if (add.length >= MAX_IMPORT) break;
      const g = f && f.geometry, p = (f && f.properties) || {};
      if (!g) { skip.other++; continue; }
      const ot = String(p.objectType || p.object_type || "").toLowerCase();
      if (ot === "detection" || ot === "cell" || ot === "tile") { skip.detections++; continue; }
      if (p.dataset && p.dataset !== S.dsId) otherDs = p.dataset;
      const um = api.um();
      if (Number.isFinite(+p.microns_per_pixel) && +p.microns_per_pixel > 0 && um && Math.abs(+p.microns_per_pixel / um - 1) > 0.01) otherUm = +p.microns_per_pixel;
      const cls = p.classification && typeof p.classification === "object" ? p.classification : null;
      const nm = cut(String(p.name || (cls && cls.name) || `Imported ${l.next + add.length}`).replace(/\s+/g, " ").trim());
      const hex = colorOf(p.color) || (cls && colorOf(cls.color)) || colorOf(p.colorRGB) || null;
      let polys = null;
      if (g.type === "Polygon") { const P = polyOf(g.coordinates); if (P) polys = [P]; }
      else if (g.type === "MultiPolygon" && Array.isArray(g.coordinates)) { polys = g.coordinates.map(polyOf).filter(Boolean); if (!polys.length) polys = null; }
      else if (g.type === "MultiPoint" && Array.isArray(g.coordinates)) {
        const seen = new Set();
        for (const c of g.coordinates) {
          if (!isXY(c)) continue;
          total++;
          const r = find(+c[0], +c[1]);
          if (r < 0) { skip.unmatched++; continue; }
          seen.add(r);
        }
        if (!seen.size) { skip.other++; continue; }
        const rows = Uint32Array.from([...seen].sort((a, b) => a - b));
        add.push({ name: nm, hex, kind: "rows", pts: [], rows, layer: S.exprKey, home: { layer: S.exprKey, rows } });
        continue;
      } else { skip.other++; continue; }
      if (!polys) { skip.other++; continue; }
      for (const rings of polys) for (const [x, y] of rings[0]) { total++; if (Wd && Ht && (x < -0.02 * Wd || y < -0.02 * Ht || x > 1.02 * Wd || y > 1.02 * Ht)) outside++; }
      const st = String(p.shape_type || "").toLowerCase();
      if (polys.length === 1 && polys[0].length === 1) add.push({ name: nm, hex, kind: st === "circle" ? "circle" : "free", pts: polys[0][0] });
      else {
        // holes or several parts: a list of bins, made from the polygons (and found again from them on other layers)
        const def = { geo: polys }, rows = maskRows(polysMask(polys, pos));
        add.push({ name: nm, hex, kind: "rows", pts: [], rows, layer: S.exprKey, home: { layer: S.exprKey, rows }, def });
      }
    }
    if (!add.length) {
      msg(`Nothing to add: no Polygon or MultiPolygon features` + (skip.detections ? ` (${skip.detections} detections / cells left out: export annotations from QuPath)` : "") + ".", true);
      return;
    }
    const first = l.shapes.length;
    for (const s of add) {
      const id = l.next++;
      l.shapes.push({ id, bins: 0, shared: 0, lost: 0, ...s, hex: s.hex || api.colors[(id - 1) % api.colors.length] });
    }
    picked.clear();
    api.changed();
    const got = l.shapes.slice(first), empty = got.filter((s) => !s.bins && !s.lost).map((s) => `“${s.name}”`);
    const notes = [];
    if (feats.length > add.length + skip.detections + skip.other) notes.push(`Only the first ${MAX_IMPORT} shapes were read.`);
    if (skip.detections) notes.push(`${skip.detections} detections / cells left out (only annotations are read).`);
    if (skip.other) notes.push(`${skip.other} features that are not polygons left out.`);
    if (skip.unmatched) notes.push(`${skip.unmatched.toLocaleString()} points matched no ${W().bin} of this layer.`);
    if (outside) notes.push(`${outside.toLocaleString()} of ${total.toLocaleString()} points lie outside the image (${Wd} × ${Ht} px): the file may be in another image's pixels (a crop, a flip or another resolution).`);
    if (otherUm) notes.push(`The file says ${otherUm} µm per pixel, this dataset ${api.um()}: the coordinates are probably on another scale.`);
    if (otherDs) notes.push(`Drawn on ${otherDs}, not ${S.dsId}.`);
    if (empty.length) notes.push(`No ${W().bins} inside ${empty.join(", ")}.`);
    msg(`Added ${got.length} ${got.length === 1 ? "shape" : "shapes"}: ${got.reduce((a, s) => a + s.bins + s.lost, 0).toLocaleString()} ${W().bins} inside. ${notes.join(" ")}`.trim(),
        !!(outside || otherUm || empty.length === got.length));
  }
  function hasMultiPoint(feats) { return feats.some((f) => f && f.geometry && f.geometry.type === "MultiPoint"); }

  // ---- barcodes / cell ids -----------------------------------------------------------------------------------------
  async function gunzipText(buf) {
    const u8 = new Uint8Array(buf);
    if (u8.length >= 2 && u8[0] === 0x1f && u8[1] === 0x8b) {
      if (typeof DecompressionStream === "undefined") throw new Error("this browser cannot unpack the barcode list; a current Chrome, Edge, Firefox or Safari can");
      const st = new Blob([u8]).stream().pipeThrough(new DecompressionStream("gzip"));
      return new Response(st).text();
    }
    return new TextDecoder().decode(u8);           // a host that already unpacked it
  }
  function layerIds(e, n) {
    const S = api.state(), k = `${S.dsId}|${e.key}`;
    if (!idCache.has(k)) {
      idCache.set(k, (async () => {
        const txt = await gunzipText(await api.buf(`expr/${e.key}/${e.ids.file || "ids.txt.gz"}`));
        const ids = txt.split("\n");
        if (ids.length && ids[ids.length - 1] === "") ids.pop();
        return ids;
      })());
      idCache.get(k).catch(() => idCache.delete(k));
    }
    return idCache.get(k).then((ids) => {
      if (ids.length !== n) throw new Error(`the barcode list has ${ids.length.toLocaleString()} rows, the layer ${n.toLocaleString()}`);
      return ids;
    });
  }
  async function idsCsv() {
    const l = L(), S = api.state(), e = api.layer();
    if (!l || !l.pos || !e) throw new Error("the shapes are not ready yet");
    if (!e.ids) return null;
    const pos = l.pos, n = pos.length / 2, ids = await layerIds(e, n), col = e.ids.column === "cell_id" ? "cell_id" : "barcode";
    const parts = [`${col},bin_row,x_px,y_px,shape\n`];
    let rowsOut = 0;
    for (const s of scope()) {
      const m = maskOfShape(s, pos), nm = q(s.name), lines = [];
      for (let i = 0; i < n; i++) {
        if (!m[i]) continue;
        lines.push(`${q(ids[i])},${i},${pos[2 * i].toFixed(2)},${pos[2 * i + 1].toFixed(2)},${nm}`);
        rowsOut++;
      }
      if (lines.length) parts.push(lines.join("\n") + "\n");
    }
    return { parts, rows: rowsOut, col, layer: S.exprKey };
  }
  async function saveIds() {
    const e = api.layer(), w = W();
    if (!e || !e.ids) {
      msg(`This export has no barcode list for this layer, so the ${w.bins} cannot be saved with their barcodes or cell ids. ` +
          `It was exported with --no-ids, by a SODA from before the lists were written by default, or its ids failed the sample ID check (the export log says which). Add the list with prepare.py ids.`, true);
      return;
    }
    msg("Reading the barcode list…");
    try {
      const r = await idsCsv(), S = api.state();
      if (!r.rows) { msg(`No ${w.bins} inside the shapes.`, true); return; }
      const name = `${slug(S.dsId)}_${slug(S.exprKey)}_shape_${r.col === "cell_id" ? "cell_ids" : "barcodes"}.csv`;
      api.download(name, new Blob(r.parts, { type: "text/csv" }), "text/csv");
      msg(`Saved ${name}: ${r.rows.toLocaleString()} rows, one per ${w.bin} and shape (a ${w.bin} inside two shapes is listed for both).` +
          (api.touch ? " On a phone or tablet it goes to Downloads (iPhone and iPad: the Files app)." : ""));
    } catch (err) { msg(`The barcodes could not be read: ${err && err.message ? err.message : err}`, true); }
  }

  window.SodaShapes = {
    init(a) { api = a; buildDom(); refresh(); },
    refresh() { if (api) refresh(); },
    remap, saveDef, restore,
    // for tests from the console
    _pick(ids) { picked.clear(); ids.forEach((i) => picked.add(i)); sig = ""; refresh(); },
    _geojson: () => geoJson().gj,
    _import: async (obj) => { await importGeoJson(obj); return $("shMsg").textContent; },
    _ids: async () => { const r = await idsCsv(); return r ? r.parts.join("") : null; },
    _members: (id) => { const l = L(), s = l && l.shapes.find((x) => x.id === id); return s && l.pos ? Array.from(maskRows(maskOfShape(s, l.pos))) : null; },
  };
})();
