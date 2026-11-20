// Colour by a value: the bins (cells, squares, ROIs) of the layer on show coloured by one number each,
// or by where two genes are expressed, with a legend on the map and in Save image.
//
//   A value per bin   any per-bin number the export already holds: the binqc/ measures named by the
//                     layer's binQc block (total counts, genes detected, % mitochondrial; on a GeoMx
//                     layer also raw reads, saturation, area, nuclei, negative-probe background, LOQ,
//                     Q3 ...) and, on a GeoMx layer, each AOI's Q3 normalisation factor (manifest q3).
//                     Also the number columns of a file loaded under Your own results (ownres.js, viewer.js
//                     extraCols). On the viridis ramp the gene density uses, from a From to a to value.
//   Two genes         where gene A, gene B or both are expressed: a bin counts as expressing a gene when
//                     its value is above that gene's threshold (default: above 0, i.e. detected; on
//                     marker squares, above the marker's median square). Shown as three classes (A only,
//                     B only, both) or as a blend: A's colour plus B's, each from its threshold (at least
//                     a quarter bright) to its 99th percentile. Bins expressing neither are not drawn.
//                     The reader picks the colours ("Colours"): magenta / green, the default, which
//                     people with red-green colour blindness tell apart too (A magenta, B green, both
//                     light grey; the blend adds up to white), or the older red / green (both yellow).
//                     On the page's dark background each class colour has a contrast of at least 6:1.
//   Gene set score    a list of genes typed or pasted, or a cell type's marker genes (the dataset's own
//                     marker list, then SODA's generic one). Per bin, the MEAN over the listed genes
//                     found in the layer of the gene's normalised value:
//                       sequencing bins   log1p(10,000 x count / the bin's total count), the log1p CP10k
//                                         normalisation Draw a region uses; the total is binqc's total
//                                         counts (all genes of the layer). Without a total: log1p(count).
//                       GeoMx ROIs        log1p(count / the AOI's Q3 normalisation factor), the Q3
//                                         normalisation the ROI plots use; an AOI without one: no score.
//                       marker squares    the marker's intensity scaled from its 1st (0) to its 99th
//                                         percentile square (1).
//                     Not a Seurat module score: no control genes are subtracted, so a score reads as
//                     "how much of these genes, on the usual log scale". Genes not in the layer are
//                     left out and named. On the same viridis ramp and range controls as a value.
//
// Genes always come through the source object (viewer.js S.src.gene: loose files, a packed file or a
// tile server alike); their values are put on positions.bin rows by position, as the UMAP does.
//
// The section is in the View tab (layout.js), under More options. While it colours the bins, the genes'
// own density is not drawn (one colour per bin, as Colour by cluster does); a colour by cluster or
// cell type in the Clusters section (Analyze tab) wins over it, and the section says so. The tissue view draws its own
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
  const MODES = ["off", "value", "pair", "set"];
  const MAX_SET = 200;                             // genes in one set
  const SET_PARALLEL = 4;                          // gene files fetched at once for a set
  const PAIR_VIEWS = ["classes", "blend"];
  const ALPHA = 235;                               // as the gene density squares
  // two genes: the colours of A only, B only and both, and the blend's (a bin's colour = sa x blendA + sb x blendB
  // + min(sa, sb) x blendX, each s from 0 to 1). The reader's choice; a saved view without one gets magenta / green.
  // Contrast on the page's background (#0e1117): magenta 6.2:1, green 8.5:1, grey 15.9:1. Told apart under a
  // simulated deuteranopia (Machado 2009): magenta / green ΔE 84 (red / green: 14).
  const PAIR_PALS = ["mg", "rg"];
  const PAL = {
    mg: { a: [240, 70, 240], b: [40, 200, 90], both: [235, 235, 235], blendA: [255, 0, 255], blendB: [0, 255, 0], blendX: [0, 0, 0],
          nameA: "magenta", nameB: "green", nameBoth: "white", label: "magenta/green" },
    rg: { a: [255, 51, 85], b: [51, 221, 102], both: [255, 221, 51], blendA: [255, 0, 0], blendB: [0, 221, 0], blendX: [0, 0, 51],
          nameA: "red", nameB: "green", nameBoth: "yellow", label: "red/green" },
  };
  let pairPal = "mg";                              // kept across layers and datasets; a link or Save JSON sets it
  const pal = () => PAL[pairPal] || PAL.mg;
  const blendRgb = (p, sa, sb) => [0, 1, 2].map((i) => Math.min(255, Math.round(sa * p.blendA[i] + sb * p.blendB[i] + Math.min(sa, sb) * p.blendX[i])));
  const BLEND_FLOOR = 0.25;                        // a bin just above a threshold is a quarter bright, not black
  const GENE_CACHE = 4;                            // gene rows kept per layer

  // P: what is on show. cols: the values this layer can be coloured by. cur: the loaded values.
  let P = blank();
  function blank() {
    return { ds: null, key: null, seq: 0, mode: "off", cols: [], col: null, scale: "lin", lo: null, hi: null,
             auto: true, cur: null, pos: null, rowOf: null, colors: null, t: null, rgb: null, ver: 0, busy: null,
             genes: new Map(), pair: { a: null, b: null, ta: null, tb: null, view: "classes" }, pd: null,
             set: { genes: [], type: "" }, lists: null };
  }
  const W = () => (api && api.words ? api.words() : { bin: "bin", bins: "bins", Bins: "Bins", gene: "gene", Gene: "Gene", genes: "genes" });
  const say = (t) => t.replace(/\bbins\b/g, W().bins).replace(/\bbin\b/g, W().bin);
  const fmt = (v) => {
    if (v === null || v === undefined || !Number.isFinite(v)) return "–";
    const a = Math.abs(v);
    if (Number.isInteger(v)) return v.toLocaleString();
    return a >= 1000 ? Math.round(v).toLocaleString() : a >= 100 ? v.toFixed(0) : a >= 10 ? v.toFixed(1) : a >= 1 ? v.toFixed(2) : v.toPrecision(2);
  };
  const rgbCss = (c) => `rgb(${c[0]},${c[1]},${c[2]})`;

  // ---- the section ---------------------------------------------------------------------------------
  function buildDom() {
    const sec = document.createElement("section");
    sec.className = "panel"; sec.id = "panelPaint"; sec.dataset.tab = "view"; sec.hidden = true;
    sec.innerHTML = `<h2>Colour by a value</h2>
      <label class="row slim"><span>Colour by</span><select id="paintMode" title="What colours the bins: nothing (the genes as picked above), one value per bin, or where two genes are expressed">
        <option value="off">Nothing (genes as picked)</option><option value="value">A value per bin</option><option value="pair">Two genes</option><option value="set">Gene set score</option></select></label>
      <div id="paintValueBox" hidden>
        <label class="row slim"><span>Value</span><select id="paintCol"></select></label>
      </div>
      <div id="paintPairBox" class="paint-pair" hidden>
        <datalist id="paintGeneList"></datalist>
        <div class="paint-gene"><label><span id="paintLblA">Gene A</span><input id="paintGeneA" type="text" list="paintGeneList" autocomplete="off" spellcheck="false" placeholder="name, then Enter" /></label>
          <label><span>above</span><input id="paintThrA" type="number" step="any" inputmode="decimal" title="A bin expresses it when its value is above this" /></label></div>
        <div class="paint-gene"><label><span id="paintLblB">Gene B</span><input id="paintGeneB" type="text" list="paintGeneList" autocomplete="off" spellcheck="false" placeholder="name, then Enter" /></label>
          <label><span>above</span><input id="paintThrB" type="number" step="any" inputmode="decimal" title="A bin expresses it when its value is above this" /></label></div>
        <label class="row slim"><span>Show</span><select id="paintPairView"><option value="classes">Both / one only</option><option value="blend">Blend (magenta + green = white)</option></select></label>
        <label class="row slim"><span>Colours</span><select id="paintPairPal" title="Magenta / green: told apart with red-green colour blindness too (both: light grey). Red / green: the older colours (both: yellow)."><option value="mg">magenta/green</option><option value="rg">red/green</option></select></label>
      </div>
      <div id="paintSetBox" class="paint-set" hidden>
        <label class="row slim"><span>Cell type</span><select id="paintSetType" title="Fill the list with a cell type's marker genes and score it"><option value="">Your own list</option></select></label>
        <textarea id="paintSetGenes" rows="3" spellcheck="false" autocomplete="off" aria-label="Genes of the set" placeholder="Gene names, separated by spaces, commas or new lines"></textarea>
        <div class="btnrow"><button id="paintSetRun" type="button" title="Score every bin by this list (Ctrl+Enter in the box)">Score</button></div>
      </div>
      <div id="paintRangeBox" class="paint-range" hidden>
        <label class="row slim" id="paintScaleRow"><span>Scale</span><select id="paintScale"><option value="lin">Linear</option><option value="log">Log (log10 of 1 + value)</option></select></label>
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
    for (const w of ["A", "B"]) {
      const inp = $(`paintGene${w}`);
      inp.oninput = () => fillList(inp.value);
      inp.onfocus = () => fillList(inp.value);
      inp.onchange = () => setPairGene(w.toLowerCase(), inp.value);
      inp.onkeydown = (e) => { if (e.key === "Enter") { e.preventDefault(); setPairGene(w.toLowerCase(), inp.value); } };
      $(`paintThr${w}`).onchange = (e) => { P.pair[`t${w.toLowerCase()}`] = num(e.target); pairCompute(); };
    }
    $("paintSetRun").onclick = () => setScore(parseList($("paintSetGenes").value), $("paintSetType").value);
    $("paintSetGenes").onkeydown = (e) => { if (e.key === "Enter" && (e.ctrlKey || e.metaKey)) { e.preventDefault(); $("paintSetRun").click(); } };
    $("paintSetGenes").oninput = () => { if ($("paintSetType").value) $("paintSetType").value = ""; };
    $("paintSetType").onchange = (e) => {
      const t = (P.lists || []).find((x) => x.id === e.target.value);
      if (t) { $("paintSetGenes").value = t.genes.join(", "); setScore(t.genes, t.id); }
    };
    $("paintPairView").onchange = (e) => { P.pair.view = PAIR_VIEWS.includes(e.target.value) ? e.target.value : "classes"; pairCompute(); };
    $("paintPairPal").onchange = (e) => { pairPal = PAIR_PALS.includes(e.target.value) ? e.target.value : "mg"; pairCompute(); };
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
    if (api && api.extraCols) for (const c of api.extraCols(key) || []) cols.push(c);     // your own results' numbers (ownres.js)
    return cols;
  }
  const hasGenes = () => !!(api && api.geneIndex().length && P.key);

  // a new layer (or dataset): what it offers, and the state carried over when it offers the same
  function onLayer(dsId, key, e, src) {
    const keep = P.ds === dsId && key ? { mode: P.mode, col: P.col, scale: P.scale, lo: P.lo, hi: P.hi, auto: P.auto, pair: { ...P.pair }, set: { ...P.set } } : null;
    P = { ...blank(), seq: P.seq + 1, ds: dsId, key };
    P.cols = columnsOf(e, key, src);
    const sel = $("paintCol");
    sel.innerHTML = P.cols.map((c) => `<option value="${esc(c.id)}">${esc(say(c.label))}${c.unit && !/^(genes|targets|nuclei)$/.test(c.unit) ? ` (${esc(c.unit)})` : ""}</option>`).join("");
    const opt = (v) => $("paintMode").querySelector(`option[value=${v}]`);
    opt("value").hidden = opt("value").disabled = !P.cols.length;
    opt("pair").hidden = opt("pair").disabled = !hasGenes();
    opt("pair").textContent = `Two ${W().genes}`;
    opt("set").hidden = opt("set").disabled = !hasGenes();
    opt("set").textContent = `${W().Gene} set score`;
    $("paintSetGenes").placeholder = `${W().Gene} names, separated by spaces, commas or new lines`;
    $("paintSetType").innerHTML = `<option value="">Your own list</option>`;
    $("paintLblA").textContent = `${W().Gene} A`; $("paintLblB").textContent = `${W().Gene} B`;
    $("panelPaint").hidden = !P.cols.length && !hasGenes();
    $("paintGeneA").value = $("paintGeneB").value = ""; $("paintThrA").value = $("paintThrB").value = "";
    if (!keep || keep.mode !== "set") $("paintSetGenes").value = "";
    msg("");
    if (keep && keep.mode === "value" && P.cols.some((c) => c.id === keep.col)) {
      Object.assign(P, { col: keep.col, scale: keep.scale, lo: keep.lo, hi: keep.hi, auto: keep.auto });
      setMode("value", true);
    } else if (keep && keep.mode === "pair" && hasGenes() && keep.pair.a && keep.pair.b && lookup(keep.pair.a) && lookup(keep.pair.b)) {
      P.pair = { ...keep.pair, ta: null, tb: null };                 // thresholds are per layer: back to the default
      setMode("pair", true);
    } else if (keep && keep.mode === "set" && hasGenes() && keep.set.genes.length) {
      P.set = { genes: keep.set.genes, type: "" };           // a cell type's list is per dataset: kept as a plain list
      setMode("set", true);
    } else { $("paintMode").value = "off"; sync(); }
  }

  // keep: put back a state (a layer switch, a link) rather than start the mode afresh
  async function setMode(m, keep) {
    P.mode = MODES.includes(m) && !$("paintMode").querySelector(`option[value=${m}]`).disabled ? m : "off";
    $("paintMode").value = P.mode;
    if (P.mode === "value") await setCol(keep && P.cols.some((c) => c.id === P.col) ? P.col : P.cols[0].id, !!keep);
    else if (P.mode === "pair") await pairLoad();
    else if (P.mode === "set") {
      loadLists();
      if (P.set.genes.length) { $("paintSetGenes").value = P.set.genes.join(", "); await setScore(P.set.genes, P.set.type, !!keep); }
      else { P.cur = null; clearOut(); P.ver++; sync(); api.render(); msg(`Type or paste ${W().genes}, or pick a cell type, then Score.`); }
    }
    else { P.cur = null; clearOut(); P.ver++; sync(); api.render(); }
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
      if (seq !== P.seq || P.col !== id || P.mode !== "value") return;  // another layer, value or mode meanwhile
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

  // ---- genes: through the source object, onto rows ----------------------------------------------------
  // the dataset's own name for a name typed or read from a link: exact first, then ignoring case
  function lookup(name) {
    const n = String(name || "").trim(); if (!n) return null;
    const idx = api.geneIndex(), lo = n.toLowerCase();
    const hit = idx.find((g) => g.n === n) || idx.find((g) => g.n.toLowerCase() === lo);
    return hit ? hit.n : null;
  }
  function fillList(q) {
    const dl = $("paintGeneList"); q = String(q || "").trim().toLowerCase();
    const idx = api.geneIndex(), out = [];
    for (const g of idx) { if (!q || g.n.toLowerCase().startsWith(q)) out.push(g.n); if (out.length >= 30) break; }
    if (q && out.length < 30) for (const g of idx) { const n = g.n.toLowerCase(); if (!n.startsWith(q) && n.includes(q)) out.push(g.n); if (out.length >= 30) break; }
    dl.innerHTML = out.map((n) => `<option value="${esc(n)}"></option>`).join("");
  }
  const ukey = (x, y) => Math.round(x * 100) * 1e8 + Math.round(y * 100);      // as viewer.js's UMAP rows
  function rowOf() {
    if (!P.rowOf) { const m = new Map(), p = P.pos; for (let r = 0; r < p.length / 2; r++) m.set(ukey(p[2 * r], p[2 * r + 1]), r); P.rowOf = m; }
    return P.rowOf;
  }
  // one gene's value on every row (0 where its file has no entry), not kept (a set's many genes)
  async function geneRowsOnce(name) {
    const hit = P.genes.get(name);
    if (hit) return hit;
    const [xyv] = await Promise.all([api.gene(name), ensurePos()]);
    const m = rowOf(), out = new Float32Array(P.pos.length / 2);
    for (let i = 0; i + 2 < xyv.length; i += 3) { const r = m.get(ukey(xyv[i], xyv[i + 1])); if (r !== undefined) out[r] = xyv[i + 2]; }
    return { v: out };
  }
  // one gene's value on every row (0 where its file has no entry): fetched once, a few kept
  function geneRows(name) {
    const hit = P.genes.get(name);
    if (hit) { P.genes.delete(name); P.genes.set(name, hit); return hit; }
    const seq = P.seq;
    const p = (async () => {
      const [xyv] = await Promise.all([api.gene(name), ensurePos()]);
      if (seq !== P.seq) throw new Error("the layer changed");
      const m = rowOf(), out = new Float32Array(P.pos.length / 2);
      let miss = 0;
      for (let i = 0; i + 2 < xyv.length; i += 3) { const r = m.get(ukey(xyv[i], xyv[i + 1])); if (r === undefined) miss++; else out[r] = xyv[i + 2]; }
      return { v: out, miss };
    })();
    P.genes.set(name, p);
    p.catch(() => { if (P.genes.get(name) === p) P.genes.delete(name); });
    while (P.genes.size > GENE_CACHE) P.genes.delete(P.genes.keys().next().value);
    return p;
  }

  // ---- two genes -------------------------------------------------------------------------------------
  async function setPairGene(w, typed) {
    const n = lookup(typed);
    const inp = $(`paintGene${w.toUpperCase()}`);
    if (!n) { if (String(typed || "").trim()) msg(`No ${W().gene} named "${String(typed).trim().slice(0, 60)}" in this layer.`, true); return; }
    inp.value = n;
    if (P.pair[w] === n && P.pd) return;
    P.pair[w] = n; P.pair[`t${w}`] = null;
    if (P.mode !== "pair") { P.mode = "pair"; $("paintMode").value = "pair"; }
    await pairLoad();
  }
  async function pairLoad() {
    const { a, b } = P.pair;
    sync();
    if (!a || !b) { clearOut(); P.pd = null; P.ver++; api.render(); msg(`Pick two ${W().genes}: type a name and press Enter.`); return; }
    const seq = P.seq;
    msg(`Loading ${a} and ${b}…`);
    try {
      const [ga, gb] = await Promise.all([geneRows(a), geneRows(b)]);
      if (seq !== P.seq || P.pair.a !== a || P.pair.b !== b || P.mode !== "pair") return;
      P.pd = { a, b, va: ga.v, vb: gb.v, miss: ga.miss + gb.miss };
      msg(P.pd.miss ? `${P.pd.miss} values were not on a ${W().bin} of this layer and are left out.` : "");
    } catch (err) {
      if (seq !== P.seq) return;
      P.pd = null; msg(`Could not load ${a} / ${b}: ${err && err.message ? err.message : err}`, true);
    }
    pairCompute();
  }
  const median = (v) => { const s = Array.from(v).sort((x, y) => x - y); return s.length ? s[Math.floor((s.length - 1) / 2)] : 0; };
  // the threshold of a gene: as set, else 0 (counts) or the median square (markers)
  const thrOf = (w) => { const t = P.pair[`t${w}`]; return Number.isFinite(t) ? t : api.isProt() && P.pd ? median(P.pd[`v${w}`]) : 0; };
  function pairCompute() {
    clearOut();
    const D = P.pd;
    if (P.mode === "pair" && D && P.pos) {
      const n = D.va.length, ta = thrOf("a"), tb = thrOf("b"), view = P.pair.view;
      const t = new Float32Array(n).fill(-1), rgb = new Uint8Array(3 * n), col = new Uint8Array(4 * n), cls = new Uint8Array(n);
      const counts = { a: 0, b: 0, both: 0, none: 0 };
      // the blend's top: the 99th percentile of the values above the threshold
      const top = (v, th) => { const s = []; for (let r = 0; r < n; r++) if (v[r] > th) s.push(v[r]); s.sort((x, y) => x - y); return s.length ? Math.max(s[Math.floor((s.length - 1) * 0.99)], th + 1e-6) : th + 1; };
      const hiA = view === "blend" ? top(D.va, ta) : 0, hiB = view === "blend" ? top(D.vb, tb) : 0;
      const sOf = (v, th, hi) => (v > th ? BLEND_FLOOR + (1 - BLEND_FLOOR) * Math.min(1, (v - th) / Math.max(hi - th, 1e-6)) : 0);
      const p = pal();
      for (let r = 0; r < n; r++) {
        const ea = D.va[r] > ta, eb = D.vb[r] > tb;
        const k = ea && eb ? "both" : ea ? "a" : eb ? "b" : "none";
        counts[k]++;
        if (k === "none") continue;
        cls[r] = k === "both" ? 3 : k === "a" ? 1 : 2;
        let c;
        if (view === "blend") {
          const sa = sOf(D.va[r], ta, hiA), sb = sOf(D.vb[r], tb, hiB);
          c = blendRgb(p, sa, sb);
          t[r] = Math.max(sa, sb);
        } else { c = p[k]; t[r] = k === "both" ? 1 : 0.6; }
        rgb[3 * r] = c[0]; rgb[3 * r + 1] = c[1]; rgb[3 * r + 2] = c[2];
        col[4 * r] = c[0]; col[4 * r + 1] = c[1]; col[4 * r + 2] = c[2]; col[4 * r + 3] = ALPHA;
      }
      Object.assign(D, { ta, tb, view, counts, cls, hiA, hiB, pal: pairPal });
      P.t = t; P.rgb = rgb; P.colors = col;
    }
    P.ver++;
    sync();
    api.render();
  }

  // ---- gene set score ------------------------------------------------------------------------------
  const parseList = (t) => String(t || "").split(/[\s,;]+/).map((x) => x.trim()).filter(Boolean);
  // the cell types to pick from: the dataset's own marker list (celltypes.json), then SODA's generic one
  async function loadLists() {
    if (P.lists) return;
    const seq = P.seq, out = [], e = api.expr(), src = api.src();
    const add = (types, group, pre) => {
      for (const t of Array.isArray(types) ? types : []) {
        if (!t || !t.name || !Array.isArray(t.genes) || !t.genes.length) continue;
        out.push({ id: `${pre}:${t.name}`, group, name: String(t.name), genes: t.genes.map(String).slice(0, MAX_SET) });
      }
    };
    if (e && e.umap && e.umap.celltypes && src && src.json) {
      try {
        const j = await src.json(`expr/${P.key}/umap/celltypes.json`);
        add(j && j.types, (j && j.markerSet && j.markerSet.name) || "This dataset's marker list", "own");
      } catch (err) { /* none of its own */ }
    }
    try {
      const r = await fetch(api.markerFile(), { cache: "no-cache" });
      const j = r.ok ? await r.json() : null;
      add(j && j.types, (j && j.name) || "Generic marker list", "soda");
    } catch (err) { /* the page's lists are optional */ }
    if (seq !== P.seq) return;
    P.lists = out;
    const groups = [...new Set(out.map((t) => t.group))];
    $("paintSetType").innerHTML = `<option value="">Your own list</option>` + groups.map((g) => `<optgroup label="${esc(g)}">` +
      out.filter((t) => t.group === g).map((t) => `<option value="${esc(t.id)}">${esc(t.name)} (${t.genes.length})</option>`).join("") + "</optgroup>").join("");
    $("paintSetType").value = P.set.type && out.some((t) => t.id === P.set.type) ? P.set.type : "";
  }
  // which normalisation this layer's genes get (see the top of the file)
  function setKind() {
    const e = api.expr();
    if (api.isProt()) return "prot";
    if (api.isRois()) return e && e.q3 && Array.isArray(e.q3.normFactors) && e.q3.normFactors.length === P.pos.length / 2 ? "q3" : "log1p";
    return P.cols.some((c) => c.id === "qc:total") ? "cp10k" : "log1p";
  }
  const KIND_TXT = {
    cp10k: "mean of log1p(10,000 × count / bin total)",
    q3: "mean of log1p(count / Q3 factor)",
    log1p: "mean of log1p(count)",
    prot: "mean of each marker scaled 1st–99th percentile",
  };
  async function setScore(list, type, keepRange) {
    const seq = P.seq, want = [], missing = [];
    for (const x of list) { const n = lookup(x); if (n) { if (!want.includes(n)) want.push(n); } else if (!missing.includes(x)) missing.push(x); }
    if (want.length > MAX_SET) want.length = MAX_SET;
    if (P.mode !== "set") { P.mode = "set"; $("paintMode").value = "set"; }
    P.set = { genes: want, type: type || "" };
    const gone = () => (missing.length ? `Not in this layer, left out: ${missing.slice(0, 8).join(", ")}${missing.length > 8 ? ` and ${missing.length - 8} more` : ""}.` : "");
    if (!want.length) {
      P.cur = null; recolor();
      msg(list.length ? `None of these are ${W().genes} of this layer: ${missing.slice(0, 8).join(", ")}${missing.length > 8 ? " …" : ""}.` : `Type or paste ${W().genes}, or pick a cell type, then Score.`, !!list.length);
      return;
    }
    try {
      await ensurePos();
      if (seq !== P.seq) return;
      const kind = setKind(), n = P.pos.length / 2, sum = new Float64Array(n);
      let div = null;
      if (kind === "cp10k") div = await P.cols.find((c) => c.id === "qc:total").load();
      else if (kind === "q3") div = Float32Array.from(api.expr().q3.normFactors, (x) => (typeof x === "number" && x > 0 ? x : NaN));
      if (div && div.length !== n) throw new Error(`${div.length} totals for ${n} ${W().bins}`);
      let next = 0, done = 0;
      const one = async () => {
        while (next < want.length) {
          const { v } = await geneRowsOnce(want[next++]);
          if (seq !== P.seq) return;
          if (kind === "cp10k") { for (let r = 0; r < n; r++) if (v[r] > 0 && div[r] > 0) sum[r] += Math.log1p((1e4 * v[r]) / div[r]); }
          else if (kind === "q3") { for (let r = 0; r < n; r++) if (v[r] > 0) sum[r] += Math.log1p(v[r] / div[r]); }
          else if (kind === "log1p") { for (let r = 0; r < n; r++) if (v[r] > 0) sum[r] += Math.log1p(v[r]); }
          else {
            const q = Float32Array.from(v).sort(), lo = q[Math.round((n - 1) * 0.01)], hi = q[Math.round((n - 1) * 0.99)], inv = 1 / Math.max(hi - lo, 1e-6);
            for (let r = 0; r < n; r++) sum[r] += Math.min(1, Math.max(0, (v[r] - lo) * inv));
          }
          done++;
          if (want.length > 4 && P.set.genes === want) msg(`Scoring: ${done} of ${want.length} ${W().genes}…`);
        }
      };
      msg(`Scoring ${want.length} ${want.length === 1 ? W().gene : W().genes}…`);
      await Promise.all(Array.from({ length: Math.min(SET_PARALLEL, want.length) }, one));
      if (seq !== P.seq || P.mode !== "set" || P.set.genes !== want) return;   // another layer, mode or list meanwhile
      const score = new Float32Array(n);
      for (let r = 0; r < n; r++) score[r] = sum[r] / want.length;          // a q3 row without a factor stays NaN: no score
      P.cur = { id: "set", label: want.length <= 3 ? `Score of ${want.join(", ")}` : `Score of ${want.length} ${W().genes}`, unit: "", v: score, kind, genes: want };
      P.scale = "lin";
      if (!keepRange || !Number.isFinite(P.lo) || !Number.isFinite(P.hi) || !(P.hi > P.lo)) { P.auto = true; autoRange(); }
      msg(gone());
    } catch (err) {
      if (seq !== P.seq) return;
      P.cur = null; msg(`Could not score the set: ${err && err.message ? err.message : err}`, true);
    }
    recolor();
  }

  // ---- a value per bin ----------------------------------------------------------------------------------
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
  const clearOut = () => { P.colors = null; P.t = null; P.rgb = null; };
  // t per row (0..1, -1 = no value) and the colours
  function recolor() {
    clearOut();
    if ((P.mode === "value" || P.mode === "set") && P.cur && P.pos) {
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
    $("paintPairBox").hidden = P.mode !== "pair";
    $("paintRangeBox").hidden = !((P.mode === "value" || P.mode === "set") && P.cur);
    $("paintSetBox").hidden = P.mode !== "set";
    $("paintScaleRow").hidden = P.mode === "set";
    $("paintScale").value = P.scale;
    $("paintLo").value = P.lo === null ? "" : String(+P.lo.toPrecision(6));
    $("paintHi").value = P.hi === null ? "" : String(+P.hi.toPrecision(6));
    if (P.mode === "pair") {
      if (P.pair.a && document.activeElement !== $("paintGeneA")) $("paintGeneA").value = P.pair.a;
      if (P.pair.b && document.activeElement !== $("paintGeneB")) $("paintGeneB").value = P.pair.b;
      $("paintPairView").value = P.pair.view;
      const p = pal(), bo = $("paintPairView").querySelector("option[value=blend]");
      bo.textContent = `Blend (${p.nameA} + ${p.nameB} = ${p.nameBoth})`;
      $("paintPairPal").value = pairPal;
      const D = P.pd;
      for (const w of ["a", "b"]) {
        const el = $(`paintThr${w.toUpperCase()}`);
        el.placeholder = D ? fmt(thrOf(w)) : "";
        if (document.activeElement !== el) el.value = Number.isFinite(P.pair[`t${w}`]) ? String(P.pair[`t${w}`]) : "";
      }
    }
    if (on && api.overridden()) msg(say("The bins are coloured by the Clusters section's choice (Analyze tab): set Colour bins by back to genes there to see this colouring."));
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
  // the UMAP, Side by side and the cell outlines: { key, t, rgb } per row, as viewer.js's gene mix; alpha: the
  // outlines fill every coloured cell at the dots' one alpha (the colour says how much, not how see-through)
  function mix(nRows) {
    if (!active() || P.t.length !== nRows) return null;
    return { key: `paint${P.ver}|${P.key}`, t: P.t, rgb: P.rgb, alpha: ALPHA };
  }
  function valueAt(r) {
    return (P.mode === "value" || P.mode === "set") && P.cur && r >= 0 && r < P.cur.v.length ? P.cur.v[r] : NaN;
  }
  const unitTxt = () => (P.cur && P.cur.unit && P.cur.unit !== "genes" ? ` ${P.cur.unit}` : "");
  const CLASS_TXT = ["neither", "A only", "B only", "both"];
  function tip(r) {
    if (!active()) return null;
    if (P.mode === "pair") {
      const D = P.pd; if (!D || r < 0 || r >= D.va.length) return null;
      const k = D.cls[r];
      return `${D.a} ${fmt(D.va[r])} · ${D.b} ${fmt(D.vb[r])}: ${k === 3 ? "both" : k === 1 ? `${D.a} only` : k === 2 ? `${D.b} only` : CLASS_TXT[0]}`;
    }
    const x = valueAt(r);
    if (P.mode === "set") return `${W().Gene} set score: ${Number.isFinite(x) ? x.toFixed(3) : "no score"}`;
    return `${say(P.cur.label)}: ${Number.isFinite(x) ? fmt(x) + unitTxt() : "no value"}`;
  }
  const pct = (k, n) => (n ? ` (${((100 * k) / n).toFixed(1)}%)` : "");
  // the key on the map
  function keyHtml() {
    if (!active()) return "";
    if (P.mode === "pair") {
      const D = P.pd, n = D.va.length, c = D.counts, thr = `${W().Gene}s count above ${fmt(D.ta)} (${esc(D.a)}) and ${fmt(D.tb)} (${esc(D.b)}).`;
      const p = PAL[D.pal] || PAL.mg, up = (t) => t.charAt(0).toUpperCase() + t.slice(1);
      if (D.view === "blend") {
        return `<div class="genes"><span class="g"><i style="background:linear-gradient(90deg,#000,${rgbCss(p.blendA)})"></i>${esc(D.a)}</span>` +
          `<span class="g"><i style="background:linear-gradient(90deg,#000,${rgbCss(p.blendB)})"></i>${esc(D.b)}</span>` +
          `<span class="g"><i style="background:${rgbCss(blendRgb(p, 1, 1))}"></i>both high</span></div>` +
          `<div class="note">${up(p.nameA)} for ${esc(D.a)} plus ${p.nameB} for ${esc(D.b)}: ${p.nameBoth} where both are high. ${thr} ${c.both.toLocaleString()} ${esc(W().bins)} have both. ${esc(W().Bins)} with neither are not drawn.</div>`;
      }
      const row = (rgb, label, k) => `<span class="g"><i style="background:${rgbCss(rgb)}"></i>${esc(label)} · ${k.toLocaleString()}${pct(k, n)}</span>`;
      return `<div class="genes">${row(p.both, `${D.a} and ${D.b}`, c.both)}${row(p.a, `${D.a} only`, c.a)}${row(p.b, `${D.b} only`, c.b)}</div>` +
        `<div class="note">${esc(W().Bins)} per class, of ${n.toLocaleString()}. ${thr} ${c.none.toLocaleString()} with neither are not drawn.</div>`;
    }
    if (P.mode === "set") {
      const g = P.cur.genes;
      return `<div class="ramp"><span>${esc(fmt(P.lo))}</span><i style="background:${api.rampCss()}"></i><span>${esc(fmt(P.hi))}</span></div>` +
        `<div class="note">${esc(W().Gene)} set score per ${esc(W().bin)}: ${esc(say(KIND_TXT[P.cur.kind]))}, over ${g.length} ${esc(g.length === 1 ? W().gene : W().genes)} ` +
        `(${esc(g.slice(0, 6).join(", "))}${g.length > 6 ? ", …" : ""}).</div>`;
    }
    const log = P.scale === "log";
    return `<div class="ramp"><span>${esc(fmt(P.lo))}</span><i style="background:${api.rampCss()}"></i><span>${esc(fmt(P.hi))}</span></div>` +
      `<div class="note">${esc(say(P.cur.label))} per ${esc(W().bin)}${esc(unitTxt() ? ` (${unitTxt().trim()})` : "")}${log ? ", log scale" : ""}. ` +
      `Below ${esc(fmt(P.lo))} or above ${esc(fmt(P.hi))}: the end colours.</div>`;
  }
  // Save image's legend
  function legend() {
    if (!active()) return [];
    if (P.mode === "pair") {
      const D = P.pd, c = D.counts, n = D.va.length, title = `${D.a} and ${D.b} (above ${fmt(D.ta)} / ${fmt(D.tb)})`, p = PAL[D.pal] || PAL.mg;
      if (D.view === "blend") return [{ title, items: [{ label: D.a, stops: [[0, 0, 0], p.blendA] }, { label: D.b, stops: [[0, 0, 0], p.blendB] }, { label: `both high (${c.both.toLocaleString()} ${W().bins})`, rgb: blendRgb(p, 1, 1) }] }];
      return [{ title, items: [{ label: `both · ${c.both.toLocaleString()}${pct(c.both, n)}`, rgb: p.both },
        { label: `${D.a} only · ${c.a.toLocaleString()}${pct(c.a, n)}`, rgb: p.a }, { label: `${D.b} only · ${c.b.toLocaleString()}${pct(c.b, n)}`, rgb: p.b }] }];
    }
    const lut = api.lut(), stops = [];
    for (let i = 0; i <= 6; i++) { const k = Math.round((i / 6) * 255) * 3; stops.push([lut[k], lut[k + 1], lut[k + 2]]); }
    if (P.mode === "set") {
      const g = P.cur.genes;
      return [{ title: `${W().Gene} set score: ${g.slice(0, 5).join(", ")}${g.length > 5 ? ` and ${g.length - 5} more` : ""}`,
        items: [{ label: `${fmt(P.lo)} to ${fmt(P.hi)}, ${say(KIND_TXT[P.cur.kind])}`, stops }] }];
    }
    return [{ title: `${say(P.cur.label)}${unitTxt() ? ` (${unitTxt().trim()})` : ""}`, items: [{ label: `${fmt(P.lo)} to ${fmt(P.hi)}${P.scale === "log" ? " (log scale)" : ""}`, stops }] }];
  }
  function summary() {
    if (P.mode === "pair" && P.pd) return `${P.pd.a} + ${P.pd.b}`;
    if (P.mode === "set" && P.cur) return `${W().gene} set (${P.cur.genes.length})`;
    if (P.mode === "off" || !P.cur) return $("panelPaint") && !$("panelPaint").hidden ? "off" : "";
    return say(P.cur.label).toLowerCase();
  }

  // ---- Save JSON and the link (viewlink.js) ------------------------------------------------------------
  function capture() {
    if (P.mode === "pair") {
      if (!P.pair.a || !P.pair.b) return null;
      const o = { m: "pair", a: P.pair.a, b: P.pair.b, pv: P.pair.view, pc: pairPal };
      if (Number.isFinite(P.pair.ta)) o.ta = P.pair.ta;
      if (Number.isFinite(P.pair.tb)) o.tb = P.pair.tb;
      return o;
    }
    if (P.mode === "set") {
      if (!P.set.genes.length) return null;
      const o = { m: "set", g: P.set.genes.slice(0, MAX_SET) };
      if (P.set.type) o.ct = P.set.type;
      if (!P.auto && Number.isFinite(P.lo) && Number.isFinite(P.hi)) { o.lo = +P.lo.toPrecision(6); o.hi = +P.hi.toPrecision(6); }
      return o;
    }
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
    h.put(o, "a", h.str(x.a, 80)); h.put(o, "b", h.str(x.b, 80)); h.put(o, "pv", h.oneOf(x.pv, PAIR_VIEWS)); h.put(o, "pc", h.oneOf(x.pc, PAIR_PALS));
    h.put(o, "ta", h.num(x.ta, -1e15, 1e15)); h.put(o, "tb", h.num(x.tb, -1e15, 1e15));
    h.put(o, "g", h.list(x.g, MAX_SET, (g) => h.str(g, 60))); h.put(o, "ct", h.str(x.ct, 120));
    return o.m ? o : null;
  }
  async function apply(x) {
    const notes = [];
    if (!x || x.m === "off") { if (P.mode !== "off") await setMode("off"); return notes; }
    if (x.m === "pair") {
      const a = lookup(x.a), b = lookup(x.b);
      if (!a || !b || !hasGenes()) { notes.push(`Not in this layer: ${[!a ? x.a : null, !b ? x.b : null].filter(Boolean).join(", ") || "its genes"}: not coloured by two genes.`); await setMode("off"); return notes; }
      P.pair = { a, b, ta: Number.isFinite(x.ta) ? x.ta : null, tb: Number.isFinite(x.tb) ? x.tb : null, view: PAIR_VIEWS.includes(x.pv) ? x.pv : "classes" };
      pairPal = PAIR_PALS.includes(x.pc) ? x.pc : "mg";             // a link or file from before the choice: the default
      P.mode = "pair"; $("paintMode").value = "pair";
      await pairLoad();
      return notes;
    }
    if (x.m === "set") {
      if (!hasGenes() || !Array.isArray(x.g) || !x.g.length) { notes.push(`This layer has no ${W().genes}: not coloured by a set.`); await setMode("off"); return notes; }
      P.mode = "set"; $("paintMode").value = "set";
      if (Number.isFinite(x.lo) && Number.isFinite(x.hi)) { P.lo = x.lo; P.hi = x.hi; P.auto = false; } else { P.auto = true; P.lo = P.hi = null; }
      $("paintSetGenes").value = x.g.join(", ");
      await loadLists();
      const type = x.ct && (P.lists || []).some((t) => t.id === x.ct) ? x.ct : "";
      $("paintSetType").value = type;
      await setScore(x.g, type, !P.auto);
      const out = x.g.filter((g) => !lookup(g));
      if (out.length) notes.push(`Not in this layer, left out of the ${W().gene} set: ${out.slice(0, 6).join(", ")}${out.length > 6 ? ` and ${out.length - 6} more` : ""}.`);
      return notes;
    }
    if (!P.cols.some((c) => c.id === x.col)) { notes.push("This dataset does not have the value the bins were coloured by: not coloured."); await setMode("off"); return notes; }
    P.col = x.col;
    P.scale = x.sc === "log" ? "log" : "lin";
    if (Number.isFinite(x.lo) && Number.isFinite(x.hi)) { P.lo = x.lo; P.hi = x.hi; P.auto = false; } else { P.auto = true; P.lo = P.hi = null; }
    P.mode = "value"; $("paintMode").value = P.mode;
    await setCol(x.col, true);
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
    _pal: () => ({ id: pairPal, ...PAL[pairPal] }),
  };
})();
