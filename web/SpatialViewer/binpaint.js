// Colour by a value: the bins (cells, squares, ROIs) of the layer on show coloured by one number each,
// on the viridis ramp the gene density uses, with its legend on the map and in Save image.
//
//   A value per bin   any per-bin number the export already holds: the binqc/ measures named by the
//                     layer's binQc block (total counts, genes detected, % mitochondrial; on a GeoMx
//                     layer also raw reads, saturation, area, nuclei, negative-probe background, LOQ,
//                     Q3 ...) and, on a GeoMx layer, each AOI's Q3 normalisation factor (manifest q3).
//
// The section is in the View tab (layout.js), under the genes. While it colours the bins, the genes'
// own density is not drawn (one colour per bin, as Colour by cluster does); a colour by cluster or
// cell type in the Clusters tab wins over it, and the section says so. The tissue view draws its own
// dots here (layer()); the UMAP, Side by side and the cell outlines take the same colours per row
// through viewer.js's gene mix (mix()).
//
// The state travels in Save JSON and in "Copy link to this view" as the view block's `paint` entry
// (viewlink.js extend()): absent in older files and links, which then open with nothing coloured.
//
// Rows are positions.bin rows throughout. Only unsigned typed arrays here: the de-identification
// check reads the signed names as sample IDs.
/* exported SodaPaint */
(function () {
  "use strict";
  let api = null;
  const $ = (id) => document.getElementById(id);
  const esc = (x) => String(x).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
  const MODES = ["off", "value"];
  const ALPHA = 235;                               // as the gene density squares

  // P: what is on show. cols: the values this layer can be coloured by. cur: the loaded values.
  let P = blank();
  function blank() {
    return { ds: null, key: null, seq: 0, mode: "off", cols: [], col: null, scale: "lin", lo: null, hi: null,
             auto: true, cur: null, pos: null, colors: null, ver: 0, busy: null };
  }
  const W = () => (api && api.words ? api.words() : { bin: "bin", bins: "bins", Bins: "Bins" });
  const say = (t) => t.replace(/\bbins\b/g, W().bins).replace(/\bbin\b/g, W().bin);
  const fmt = (v) => {
    if (v === null || v === undefined || !Number.isFinite(v)) return "–";
    const a = Math.abs(v);
    if (Number.isInteger(v)) return v.toLocaleString();
    return a >= 1000 ? Math.round(v).toLocaleString() : a >= 100 ? v.toFixed(0) : a >= 10 ? v.toFixed(1) : a >= 1 ? v.toFixed(2) : v.toPrecision(2);
  };

  // ---- the section ---------------------------------------------------------------------------------
  function buildDom() {
    const sec = document.createElement("section");
    sec.className = "panel"; sec.id = "panelPaint"; sec.dataset.tab = "view"; sec.hidden = true;
    sec.innerHTML = `<h2>Colour by a value</h2>
      <label class="row slim"><span>Colour by</span><select id="paintMode" title="What colours the bins: nothing (the genes as picked above) or one value per bin">
        <option value="off">Nothing (genes as picked)</option><option value="value">A value per bin</option></select></label>
      <div id="paintValueBox" hidden>
        <label class="row slim"><span>Value</span><select id="paintCol"></select></label>
      </div>
      <div id="paintRangeBox" class="paint-range" hidden>
        <label class="row slim"><span>Scale</span><select id="paintScale"><option value="lin">Linear</option><option value="log">Log (log10 of 1 + value)</option></select></label>
        <div class="paint-lohi">
          <label><span>From</span><input id="paintLo" type="number" step="any" inputmode="decimal" /></label>
          <label><span>to</span><input id="paintHi" type="number" step="any" inputmode="decimal" /></label>
          <button id="paintAuto" type="button" title="From the 1st to the 99th percentile of the values">Auto</button>
        </div>
      </div>
      <div id="paintMsg" class="meta paint-msg" role="status"></div>`;
    const side = $("sidebar");
    if (window.SodaLayout) { side.appendChild(sec); window.SodaLayout.adopt(sec); }
    else { const after = $("panelGenes"); if (after) after.after(sec); else side.appendChild(sec); }
    $("paintMode").onchange = (e) => setMode(e.target.value);
    $("paintCol").onchange = (e) => setCol(e.target.value);
    $("paintScale").onchange = (e) => { P.scale = e.target.value === "log" ? "log" : "lin"; recolor(); };
    const num = (el) => { const v = parseFloat(el.value); return Number.isFinite(v) ? v : null; };
    $("paintLo").onchange = () => { P.lo = num($("paintLo")); P.auto = false; recolor(); };
    $("paintHi").onchange = () => { P.hi = num($("paintHi")); P.auto = false; recolor(); };
    $("paintAuto").onclick = () => { P.auto = true; autoRange(); recolor(); };
  }
  const msg = (t, err) => { const el = $("paintMsg"); if (!el) return; el.textContent = t || ""; el.classList.toggle("err", !!err); };

  // ---- the values a layer offers ---------------------------------------------------------------------
  function columnsOf(e, key, src) {
    const cols = [];
    if (!e || !src || !src.buf) return cols;
    const q = e.binQc;
    if (q && Array.isArray(q.metrics)) for (const m of q.metrics) {
      if (!m || typeof m.file !== "string" || !/^[\w.-]+\.bin$/.test(m.file)) continue;
      const dir = /^[\w.-]+$/.test(q.dir || "") ? q.dir : "binqc";
      cols.push({ id: `qc:${m.key}`, label: m.label || m.key, unit: m.unit || "", scale: m.scale === "log" ? "log" : "lin",
        load: async () => {
          const b = await src.buf(`expr/${key}/${dir}/${m.file}`);
          const v = m.dtype === "uint16" ? new Uint16Array(b) : m.dtype === "uint32" ? new Uint32Array(b) : new Float32Array(b);
          return Float32Array.from(v);
        } });
    }
    const nf = e.q3 && Array.isArray(e.q3.normFactors) ? e.q3.normFactors : null;
    if (nf && nf.length === e.nBins) cols.push({ id: "q3:normFactor", label: "Q3 normalisation factor", unit: "", scale: "lin",
      load: async () => Float32Array.from(nf, (x) => (typeof x === "number" ? x : NaN)) });
    return cols;
  }

  // a new layer (or dataset): what it offers, and the state carried over when it offers the same
  function onLayer(dsId, key, e, src) {
    const keep = P.ds === dsId ? { mode: P.mode, col: P.col, scale: P.scale, lo: P.lo, hi: P.hi, auto: P.auto } : null;
    P = { ...blank(), seq: P.seq + 1, ds: dsId, key };
    P.cols = columnsOf(e, key, src);
    const sel = $("paintCol");
    sel.innerHTML = P.cols.map((c) => `<option value="${esc(c.id)}">${esc(say(c.label))}${c.unit && !/^(genes|targets|nuclei)$/.test(c.unit) ? ` (${esc(c.unit)})` : ""}</option>`).join("");
    $("panelPaint").hidden = !P.cols.length;
    if (keep && keep.mode !== "off" && P.cols.some((c) => c.id === keep.col)) {
      Object.assign(P, { col: keep.col, scale: keep.scale, lo: keep.lo, hi: keep.hi, auto: keep.auto });
      setMode(keep.mode, true);
    } else { $("paintMode").value = "off"; sync(); }
  }

  async function setMode(m, quiet) {
    P.mode = MODES.includes(m) && P.cols.length ? m : "off";
    $("paintMode").value = P.mode;
    if (P.mode === "value") await setCol(P.col && P.cols.some((c) => c.id === P.col) ? P.col : P.cols[0].id, quiet);
    else { P.cur = null; P.colors = null; P.ver++; sync(); api.render(); }
  }
  async function setCol(id, keepRange) {
    const c = P.cols.find((x) => x.id === id);
    if (!c) return;
    const seq = P.seq;
    P.col = id; $("paintCol").value = id;
    if (!keepRange) { P.scale = c.scale; P.auto = true; P.lo = P.hi = null; }
    msg(`Loading ${say(c.label).toLowerCase()}…`);
    try {
      const [v, pos] = await Promise.all([c.load(), ensurePos()]);
      if (seq !== P.seq || P.col !== id) return;                       // another layer or value meanwhile
      if (v.length * 2 !== pos.length) throw new Error(`${v.length} values for ${pos.length / 2} ${W().bins}`);
      P.cur = { id, label: c.label, unit: c.unit, v };
      if (P.auto || P.lo === null || P.hi === null) autoRange();
      msg("");
    } catch (err) {
      if (seq !== P.seq) return;
      P.cur = null; msg(`Could not load ${say(c.label).toLowerCase()}: ${err && err.message ? err.message : err}`, true);
    }
    recolor();
  }
  function ensurePos() {
    if (P.pos) return Promise.resolve(P.pos);
    const seq = P.seq;
    if (!P.busy) P.busy = api.positions().then((p) => { if (seq === P.seq) P.pos = p; P.busy = null; return p; }, (e) => { P.busy = null; throw e; });
    return P.busy;
  }

  // ---- colours -------------------------------------------------------------------------------------
  const tf = (x) => (P.scale === "log" ? Math.log10(1 + Math.max(0, x)) : x);
  function autoRange() {
    const v = P.cur && P.cur.v; if (!v) return;
    const s = Array.from(v).filter(Number.isFinite).sort((a, b) => a - b);
    if (!s.length) { P.lo = 0; P.hi = 1; return; }
    const q = (f) => s[Math.min(s.length - 1, Math.max(0, Math.round((s.length - 1) * f)))];
    let lo = q(0.01), hi = q(0.99);
    if (!(hi > lo)) { lo = s[0]; hi = s[s.length - 1]; }
    if (!(hi > lo)) hi = lo + 1;
    P.lo = lo; P.hi = hi;
  }
  // t per row (0..1, -1 = no value) and the colours
  function recolor() {
    P.colors = null; P.t = null;
    if (P.mode === "value" && P.cur && P.pos) {
      const v = P.cur.v, n = v.length, lut = api.lut(), t = new Float32Array(n), rgb = new Uint8Array(3 * n), col = new Uint8Array(4 * n);
      const a = tf(P.lo), b = tf(P.hi), inv = 1 / Math.max(b - a, 1e-12);
      for (let r = 0; r < n; r++) {
        const x = v[r];
        if (!Number.isFinite(x)) { t[r] = -1; continue; }
        const f = Math.min(1, Math.max(0, (tf(x) - a) * inv)), k = Math.round(f * 255) * 3;
        t[r] = f; rgb[3 * r] = lut[k]; rgb[3 * r + 1] = lut[k + 1]; rgb[3 * r + 2] = lut[k + 2];
        col[4 * r] = lut[k]; col[4 * r + 1] = lut[k + 1]; col[4 * r + 2] = lut[k + 2]; col[4 * r + 3] = ALPHA;
      }
      P.t = t; P.rgb = rgb; P.colors = col;
    }
    P.ver++;
    sync();
    api.render();
  }
  function sync() {
    const on = P.mode !== "off";
    $("paintValueBox").hidden = P.mode !== "value";
    $("paintRangeBox").hidden = !(P.mode === "value" && P.cur);
    $("paintScale").value = P.scale;
    $("paintLo").value = P.lo === null ? "" : String(+P.lo.toPrecision(6));
    $("paintHi").value = P.hi === null ? "" : String(+P.hi.toPrecision(6));
    if (on && api.overridden()) msg(say("The bins are coloured by the Clusters tab's choice: set Colour bins by back to genes there to see this colouring."));
    else if ($("paintMsg").textContent.startsWith("The bins are coloured")) msg("");
  }

  // ---- what viewer.js draws --------------------------------------------------------------------------
  const active = () => !!(P.mode !== "off" && P.colors && P.pos && P.key === api.layerKey());
  // the tissue view: one dot per row, the size of the gene dots
  function layer(radius, opacity) {
    if (!active()) return null;
    const n = P.pos.length / 2;
    return new api.ScatterplotLayer({
      id: `paint-${P.ds}-${P.key}`,
      data: { length: n, attributes: { getPosition: { value: P.pos, size: 2 }, getFillColor: { value: P.colors, size: 4, normalized: true } } },
      radiusUnits: "common", getRadius: radius, radiusMinPixels: 0.6, stroked: false, filled: true,
      opacity, pickable: true, parameters: api.NORMAL,
      updateTriggers: { getFillColor: [P.ver] },
    });
  }
  // the UMAP, Side by side and the cell outlines: { key, t, rgb } per row, as viewer.js's gene mix
  function mix(nRows) {
    if (!active() || P.t.length !== nRows) return null;
    return { key: `paint${P.ver}|${P.key}`, t: P.t, rgb: P.rgb };
  }
  function valueAt(r) {
    return P.cur && r >= 0 && r < P.cur.v.length ? P.cur.v[r] : NaN;
  }
  const unitTxt = () => (P.cur && P.cur.unit && P.cur.unit !== "genes" ? ` ${P.cur.unit}` : "");
  function tip(r) {
    if (!active()) return null;
    const x = valueAt(r);
    return `${say(P.cur.label)}: ${Number.isFinite(x) ? fmt(x) + unitTxt() : "no value"}`;
  }
  // the key on the map
  function keyHtml() {
    if (!active()) return "";
    const log = P.scale === "log";
    return `<div class="ramp"><span>${esc(fmt(P.lo))}</span><i style="background:${api.rampCss()}"></i><span>${esc(fmt(P.hi))}</span></div>` +
      `<div class="note">${esc(say(P.cur.label))} per ${esc(W().bin)}${esc(unitTxt() ? ` (${unitTxt().trim()})` : "")}${log ? ", log scale" : ""}. ` +
      `Below ${esc(fmt(P.lo))} or above ${esc(fmt(P.hi))}: the end colours.</div>`;
  }
  // Save image's legend
  function legend() {
    if (!active()) return [];
    const lut = api.lut(), stops = [];
    for (let i = 0; i <= 6; i++) { const k = Math.round((i / 6) * 255) * 3; stops.push([lut[k], lut[k + 1], lut[k + 2]]); }
    return [{ title: `${say(P.cur.label)}${unitTxt() ? ` (${unitTxt().trim()})` : ""}`, items: [{ label: `${fmt(P.lo)} to ${fmt(P.hi)}${P.scale === "log" ? " (log scale)" : ""}`, stops }] }];
  }
  function summary() {
    if (P.mode === "off" || !P.cur) return $("panelPaint") && !$("panelPaint").hidden ? "off" : "";
    return say(P.cur.label).toLowerCase();
  }

  // ---- Save JSON and the link (viewlink.js) ------------------------------------------------------------
  function capture() {
    if (P.mode === "off" || !P.col) return null;
    const o = { m: P.mode, col: P.col, sc: P.scale };
    if (!P.auto && Number.isFinite(P.lo) && Number.isFinite(P.hi)) { o.lo = +P.lo.toPrecision(6); o.hi = +P.hi.toPrecision(6); }
    return o;
  }
  function clean(x, h) {
    if (!h.isObj(x)) return null;
    const o = {};
    h.put(o, "m", h.oneOf(x.m, MODES)); h.put(o, "col", h.str(x.col, 80)); h.put(o, "sc", h.oneOf(x.sc, ["lin", "log"]));
    const lo = h.num(x.lo, -1e15, 1e15), hi = h.num(x.hi, -1e15, 1e15);
    if (lo !== null && hi !== null && hi > lo) { o.lo = lo; o.hi = hi; }
    return o.m ? o : null;
  }
  async function apply(x) {
    const notes = [];
    if (!x || x.m === "off") { if (P.mode !== "off") await setMode("off"); return notes; }
    if (!P.cols.some((c) => c.id === x.col)) { notes.push("This dataset does not have the value the bins were coloured by: not coloured."); await setMode("off"); return notes; }
    P.col = x.col;
    P.scale = x.sc === "log" ? "log" : "lin";
    if (Number.isFinite(x.lo) && Number.isFinite(x.hi)) { P.lo = x.lo; P.hi = x.hi; P.auto = false; } else { P.auto = true; P.lo = P.hi = null; }
    P.mode = x.m; $("paintMode").value = P.mode;
    await setCol(x.col, true);
    if (P.auto) { autoRange(); recolor(); }
    return notes.map(say);
  }

  window.SodaPaint = {
    init(a) {
      api = a;
      buildDom();
      if (window.SodaViewLink && window.SodaViewLink.extend) window.SodaViewLink.extend({ key: "paint", capture, clean, apply });
    },
    onLayer, active, layer, mix, tip, keyHtml, legend, summary, valueAt,
    // tests
    _state: () => P,
  };
})();
