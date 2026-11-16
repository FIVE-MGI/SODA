// Volcano plot: differential expression between two groups of bins, drawn as plain SVG.
// Groups: a drawn shape, a cluster of the last Circle -> UMAP run, a cluster of the section
// UMAP, the rest of the section, the rest of the run's bins, or the rest of the drawn shapes.
// The test runs in de_worker.js on the hvg/ chunks the Circle -> UMAP tool already fetches
// (plus, if asked, per-gene files for named genes). viewer.js hands over what it owns
// through SodaVolcano.init(api) and calls SodaVolcano.refresh() when the shapes or the run
// change; everything else lives here, so the viewer itself barely changes.
(function () {
  const NONE16 = 0xffff;
  const DEFAULTS = { minUmi: 10, cap: 500000, capTouch: 30000, seed: 42, minPct: 5, lfc: 1, padj: 0.05, labels: 8, maxExtra: 50 };
  const COL = { up: "#ff5a5f", down: "#4ea1ff", ns: "#6b7385", ring: "#ffffff" };
  let api = null, worker = null, running = false, R = null, hoverJ = -1, pickedJ = -1;
  const V = { lfc: DEFAULTS.lfc, padj: DEFAULTS.padj, yScale: "auto" };
  const $ = (id) => document.getElementById(id);
  const esc = (x) => String(x).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));

  // ---- DOM: a block in the Circle -> UMAP panel, a floating panel over the map -----------
  function buildDom() {
    const block = document.createElement("div");
    block.id = "volcanoBlock"; block.className = "vc-block"; block.hidden = true;
    block.innerHTML = `
      <h3 class="sub-h">Volcano plot (differential expression)</h3>
      <label class="row slim"><span>Group A</span><select id="vcA"></select></label>
      <label class="row slim"><span>vs B</span><select id="vcB"></select></label>
      <label class="row slim"><span title="Wilcoxon rank-sum (tie-corrected, as scanpy method='wilcoxon', tie_correct=True), or Welch's t on log1p CP10k (scanpy method='t-test')">Test</span>
        <select id="vcTest"><option value="wilcoxon">Wilcoxon rank-sum</option><option value="welch">Welch t-test</option></select></label>
      <label class="row slim"><span title="Test a gene only if it is counted in at least this share of the bins of A or of B">Min. % bins</span>
        <input id="vcMinPct" type="text" inputmode="decimal" value="${DEFAULTS.minPct}" /></label>
      <label class="row slim"><span title="Genes outside the section's variable genes, by name (comma or space separated, up to ${DEFAULTS.maxExtra}); each is one more download">Extra genes</span>
        <input id="vcExtra" type="text" placeholder="e.g. GAPDH, ACTB" autocomplete="off" /></label>
      <div class="btnrow"><button id="vcRun">Volcano</button><button id="vcCancel" hidden>Cancel</button><button id="vcShow" hidden title="Show the last volcano again">Show</button></div>
      <div id="vcMsg" class="meta vc-msg"></div>`;
    const anchor = $("lassoSaveJson") && $("lassoSaveJson").parentElement;
    if (anchor && anchor.parentElement) anchor.parentElement.insertBefore(block, anchor); else $("panelLasso").appendChild(block);
    if (window.SodaLayout) window.SodaLayout.adopt(block);       // layout.js: its own section in the Analyze tab

    const panel = document.createElement("div");
    panel.id = "volcanoPanel"; panel.className = "vc-panel"; panel.hidden = true;
    panel.innerHTML = `
      <div class="vc-head"><b id="vcTitle">Volcano</b>
        <span class="vc-btns"><button id="vcCsv" title="Every gene: log2FC, p, padj, % expressing, mean CP10k (no bins, no barcodes)">Save CSV</button><button id="vcClose" title="Hide">×</button></span></div>
      <div id="vcWarn" class="lp-warn" hidden></div>
      <div class="vc-ctl">
        <label title="Fold-change threshold, both sides">|log2FC| ≥ <input id="vcLfc" type="text" inputmode="decimal" value="${DEFAULTS.lfc}" /></label>
        <label title="Benjamini-Hochberg adjusted p threshold">padj ≤ <input id="vcPadj" type="text" inputmode="decimal" value="${DEFAULTS.padj}" /></label>
        <label title="A few genes with enormous -log10 padj squeeze the rest: square-root axis">y <select id="vcY"><option value="auto">auto</option><option value="linear">linear</option><option value="sqrt">√</option></select></label>
        <span id="vcCounts" class="vc-counts"></span>
      </div>
      <div id="vcHost" class="vc-host"><svg id="vcSvg" class="vc-svg"></svg><div id="vcTip" class="vc-tip" hidden></div></div>
      <div id="vcFoot" class="lp-foot"></div>`;
    document.querySelector(".main").appendChild(panel);

    $("vcRun").onclick = runVolcano;
    $("vcCancel").onclick = () => { cancel(); msg("Cancelled."); };
    $("vcShow").onclick = () => { if (R) { $("volcanoPanel").hidden = false; draw(); } };
    $("vcA").onchange = () => fillB();
    // the section's clustering can change in the Clusters panel without telling the lasso panel
    for (const id of ["vcA", "vcB"]) $(id).addEventListener("pointerdown", refresh);
    $("vcClose").onclick = () => { $("volcanoPanel").hidden = true; syncButtons(); };
    $("vcCsv").onclick = saveCsv;
    const num = (id, key, lo, hi) => { $(id).onchange = () => { const v = Number($(id).value); if (Number.isFinite(v) && v >= lo && v <= hi) V[key] = v; else $(id).value = V[key]; draw(); }; };
    num("vcLfc", "lfc", 0, 20); num("vcPadj", "padj", 1e-300, 1);
    $("vcY").onchange = () => { V.yScale = $("vcY").value; draw(); };
    const svg = $("vcSvg");
    svg.addEventListener("mousemove", (ev) => { const j = hitGene(ev); if (j !== hoverJ) { hoverJ = j; tip(ev); highlight(); } else if (j >= 0) tip(ev); });
    svg.addEventListener("mouseleave", () => { hoverJ = -1; tip(null); highlight(); });
    svg.addEventListener("click", (ev) => {
      const j = hitGene(ev); if (j < 0 || !R) return;
      pickedJ = j; hoverJ = j; highlight(); tip(ev);      // a tap has no hover: the numbers show on click too
      api.addGene(R.genes[j]);                     // onto the tissue, with the viewer's own colours
    });
    if (window.ResizeObserver) new ResizeObserver(() => { if (!$("volcanoPanel").hidden) draw(); }).observe($("vcHost"));
  }

  // the layer's own words (roi.js: ROI segments instead of bins), when the viewer hands a wording function over
  const say = (t) => (api && api.say ? api.say(t) : t);
  const msg = (t, err) => { const el = $("vcMsg"); if (!el) return; el.textContent = say(t); el.style.color = err ? "#ff8f8f" : ""; };
  function syncButtons() {
    if (!$("vcRun")) return;
    $("vcRun").disabled = running || !$("vcA").value || !$("vcB").value;
    $("vcCancel").hidden = !running;
    $("vcShow").hidden = !R || !$("volcanoPanel").hidden;
  }

  // ---- the groups on offer ------------------------------------------------------------
  // value strings: s:<shape index>, c:<run cluster>, k:<section cluster>, rest, runrest, shapesrest
  function groupsNow() {
    const S = api.state(), L = S.lasso, out = { shapes: [], run: [], section: [] };
    if (!L) return out;
    L.shapes.forEach((s, i) => { if (s.bins > 0) out.shapes.push({ v: `s:${i}`, t: `${s.name} (${s.bins.toLocaleString()} bins)` }); });
    const Rr = L.result;
    if (Rr && Rr.cluster && Rr.info.mode === "bins" && Rr.layer === S.exprKey) {
      const n = Rr.info.nClusters || 0, sizes = Rr.info.clusterSizes || [];
      for (let c = 0; c < n; c++) out.run.push({ v: `c:${c}`, t: `${(L.clNames && L.clNames[c]) || `Cluster ${c}`} of the UMAP run${sizes[c] ? ` (${sizes[c].toLocaleString()} bins)` : ""}` });
    }
    const sec = api.sectionClusters();
    if (sec) for (let c = 0; c < sec.n; c++) out.section.push({ v: `k:${c}`, t: `${sec.names[c] || `Cluster ${c}`} of the section (${sec.label || `res ${sec.res}`}${sec.sizes && sec.sizes[c] ? `, ${sec.sizes[c].toLocaleString()} bins` : ""})` });
    return out;
  }
  function fillSelect(sel, groups, extra, keep) {
    const old = keep !== undefined ? keep : sel.value;
    sel.innerHTML = "";
    const add = (parent, o) => { const e = document.createElement("option"); e.value = o.v; e.textContent = say(o.t); parent.appendChild(e); };
    for (const o of extra) add(sel, o);
    for (const [k, label] of [["shapes", "Drawn shapes"], ["run", "Clusters of the Draw a region run"], ["section", "Clusters of the whole section"]]) {
      const list = groups[k].filter((o) => o.v !== $("vcA").value || sel.id === "vcA");
      if (!list.length) continue;
      const g = document.createElement("optgroup"); g.label = label; list.forEach((o) => add(g, o)); sel.appendChild(g);
    }
    if ([...sel.options].some((o) => o.value === old)) sel.value = old;
  }
  function fillB() {
    const G = groupsNow(), a = $("vcA").value, sel = $("vcB"), had = sel.value;
    const extra = [{ v: "rest", t: "Rest of the section" }];
    if (a.startsWith("c:")) extra.unshift({ v: "runrest", t: "Rest of the UMAP run's bins" });
    if (G.shapes.length > 1 && a.startsWith("s:")) extra.push({ v: "shapesrest", t: "The other drawn shapes, together" });
    fillSelect(sel, G, extra, had);
    if (!had || ![...sel.options].some((o) => o.value === had)) {
      // a second shape when there is one, else the rest of the section
      const other = G.shapes.find((o) => o.v !== a);
      sel.value = a.startsWith("s:") && other ? other.v : a.startsWith("c:") ? "runrest" : "rest";
    }
    syncButtons();
  }

  function refresh() {
    if (!api || !$("volcanoBlock")) return;
    const S = api.state(), L = S.lasso;
    const ok = !!(L && L.hvg && !L.prot && L.pos);
    $("volcanoBlock").hidden = !ok;
    if (!ok) { if (!$("volcanoPanel").hidden) $("volcanoPanel").hidden = true; return; }
    // rebuild the menus only when what they offer changed (a rebuild closes an open menu)
    const G = groupsNow(), key = JSON.stringify(G);
    if (key !== lastGroups) { lastGroups = key; fillSelect($("vcA"), G, []); fillB(); }
    if (R && !$("volcanoPanel").hidden && R.sig !== sigNow(R.a, R.b)) { $("vcWarn").textContent = "The groups have changed since this volcano was computed: run it again."; $("vcWarn").hidden = false; }
    syncButtons();
  }
  let lastGroups = "";

  // what a result depends on: layer, shapes, run, section clustering
  function sigNow(a, b) {
    const S = api.state(), L = S.lasso, sec = api.sectionClusters();
    return [S.exprKey, L.shapes.map((s) => `${s.id}:${s.bins}`).join(","), L.result ? L.result.id : "", sec ? sec.res : "", a, b].join("|");
  }

  // ---- group codes per row: 1 = A, 2 = B (A wins a bin in both) ----------------------------
  async function groupMap(a, b) {
    const S = api.state(), L = S.lasso, n = L.pos.length / 2, g = new Uint8Array(n);
    let sec = null;
    const needSec = a.startsWith("k:") || b.startsWith("k:");
    if (needSec) { sec = await api.loadSectionClusters(); if (!sec) throw new Error("The section clusters could not be loaded."); }
    const Rr = L.result;
    const test = (v) => {
      if (v.startsWith("s:")) { const s = Number(v.slice(2)); return (r) => L.assign[r] === s; }
      if (v.startsWith("k:")) { const c = Number(v.slice(2)); return (r) => sec.cl[r] === c; }
      if (v.startsWith("c:")) {
        const c = Number(v.slice(2)), set = new Uint8Array(n);
        for (let j = 0; j < Rr.rows.length; j++) if (Rr.cluster[j] === c) set[Rr.rows[j]] = 1;
        return (r) => set[r] === 1;
      }
      if (v === "runrest") { const set = new Uint8Array(n); for (const r of Rr.rows) set[r] = 1; return (r) => set[r] === 1; }
      if (v === "shapesrest") return (r) => L.assign[r] !== NONE16;
      if (v === "rest") return () => true;
      throw new Error(`unknown group ${v}`);
    };
    const ta = test(a), tb = test(b);
    let nA = 0, nB = 0, both = 0;
    for (let r = 0; r < n; r++) {
      const inA = ta(r), inB = tb(r);
      if (inA) { g[r] = 1; nA++; if (inB && b !== "rest" && b !== "runrest" && b !== "shapesrest") both++; }
      else if (inB) { g[r] = 2; nB++; }
    }
    return { g, nA, nB, both };
  }
  const labelOf = (sel, v) => { const o = [...sel.options].find((x) => x.value === v); return o ? o.textContent.replace(/\s*\([^)]*(?:bins|ROIs)\)$/, "").replace(/\s*\(res [^)]*\)$/, "") : v; };

  function cancel() {
    if (worker && running) { worker.terminate(); worker = null; }
    running = false; syncButtons();
  }

  async function runVolcano() {
    if (running) return;
    const S = api.state(), L = S.lasso, a = $("vcA").value, b = $("vcB").value;
    if (!a || !b) { msg("Pick two groups.", true); return; }
    const minPct = Number($("vcMinPct").value);
    if (!Number.isFinite(minPct) || minPct < 0 || minPct > 50) { msg("Min. % bins must be a number from 0 to 50.", true); return; }
    const test = $("vcTest").value, nameA = labelOf($("vcA"), a), nameB = labelOf($("vcB"), b), sig = sigNow(a, b);
    running = true; syncButtons();
    const t0 = performance.now(), T = {};
    try {
      msg("Finding the bins…");
      const { g, nA, nB, both } = await groupMap(a, b);
      if (!running) return;
      if (nA < 3 || nB < 3) throw new Error(`Each group needs bins: A has ${nA}, B has ${nB}.`);
      // the chunks under either group
      const cp = L.hvg.chunkPx, keys = new Set();
      for (let r = 0; r < g.length; r++) if (g[r]) keys.add(api.chunkOf(L.pos[2 * r], L.pos[2 * r + 1], cp));
      const total = L.hvg.chunks.length, mb = L.hvg.chunks.filter((c) => keys.has(`${c.tx}_${c.ty}`)).reduce((s, c) => s + c.bytes, 0) / 1e6;
      if (api.touch && mb > 40 && !confirm(`This comparison reads ${keys.size} of ${total} data chunks (${mb.toFixed(0)} MB). Go on?`)) { msg("Not run."); return; }
      msg(`Fetching ${keys.size} data chunk(s) (${mb.toFixed(1)} MB; cached ones are not fetched again)…`);
      let { bufs, missing } = await api.fetchChunks([...keys]);
      for (let t = 0; t < 2 && missing.length; t++) {                    // a busy server can drop a few: those again
        const more = await api.fetchChunks(missing); bufs = bufs.concat(more.bufs); missing = more.missing;
      }
      if (!running) return;
      if (missing.length && !confirm(`${missing.length} of ${keys.size} data chunk(s) could not be loaded, so part of a group has no data. Run anyway on what did load?`)) { msg("Not run: data missing.", true); return; }
      T.fetch = performance.now() - t0;
      const genes = await api.hvgGenes();
      if (!genes) throw new Error("The list of variable genes (hvg/genes.json) could not be loaded.");
      // extra genes by name, outside the variable genes
      const extraNames = [], unknown = [], already = [];
      const have = new Set(genes.map((x) => x.toUpperCase()));
      const index = new Map((S.geneIndex || []).map((x) => [x.n.toUpperCase(), x.n]));
      for (const w of $("vcExtra").value.split(/[\s,;]+/).filter(Boolean)) {
        const u = w.toUpperCase();
        if (have.has(u)) { already.push(w); continue; }
        const nm = index.get(u);
        if (!nm) { unknown.push(w); continue; }
        if (!extraNames.includes(nm)) extraNames.push(nm);
      }
      if (extraNames.length > DEFAULTS.maxExtra) throw new Error(`Up to ${DEFAULTS.maxExtra} extra genes at a time (${extraNames.length} given).`);
      const extra = [], failed = [];
      for (const [i, nm] of extraNames.entries()) {
        msg(`Fetching extra gene ${i + 1} of ${extraNames.length}: ${nm}…`);
        try { extra.push({ name: nm, xyv: await api.geneXyv(nm) }); } catch (e) { failed.push(nm); }
        if (!running) return;
      }
      T.extra = performance.now() - t0 - T.fetch;
      msg("Testing…");
      if (!worker) worker = new Worker("de_worker.js?v=20261111a");
      const cap = api.touch ? DEFAULTS.capTouch : DEFAULTS.cap;
      const res = await new Promise((resolve, reject) => {
        worker.onmessage = (ev) => {
          const d = ev.data;
          if (d.type === "progress") { msg(`${d.step}… ${Math.round(d.frac * 100)}%`); return; }
          if (d.type === "error") reject(new Error(d.message)); else resolve(d);
        };
        worker.onerror = (e) => reject(new Error(e.message || "worker failed"));
        worker.postMessage({ chunks: bufs, group: g, nGenes: genes.length, extra, positions: extra.length ? L.pos : null,
                             params: { minUmi: DEFAULTS.minUmi, cap, seed: DEFAULTS.seed, minPct: minPct / 100, test } });
      });
      if (!running) return;
      T.worker = performance.now() - t0 - T.fetch - T.extra;
      R = { ...res, genes: genes.concat(extra.map((e) => e.name)), a, b, nameA, nameB, sig, test, layer: S.exprKey, dataset: S.dsId,
            nHvg: genes.length, both, missingChunks: missing.length, givenGroups: !!((api.sectionClusters() || {}).user), secs: (performance.now() - t0) / 1000, ms: { ...T, ...res.info.ms } };
      window.SodaVolcano.result = R;                // for tests from the console
      pickedJ = -1; hoverJ = -1; tip(null);
      const I = R.info;
      msg([`${I.nTested.toLocaleString()} of ${I.nGenes.toLocaleString()} genes tested in ${R.secs.toFixed(1)} s ` +
             `(A ${I.nA.toLocaleString()} bins, B ${I.nB.toLocaleString()} bins).`,
           I.lowA + I.lowB ? `Left out ${(I.lowA + I.lowB).toLocaleString()} bins with < ${I.minUmi} UMI.` : "",
           I.sampledFromA ? `A: ${I.nA.toLocaleString()} of ${I.sampledFromA.toLocaleString()} bins (seeded sample).` : "",
           I.sampledFromB ? `B: ${I.nB.toLocaleString()} of ${I.sampledFromB.toLocaleString()} bins (seeded sample).` : "",
           both ? `${both.toLocaleString()} bins in both groups count for A only.` : "",
           already.length ? `Already among the variable genes: ${already.join(", ")}.` : "",
           unknown.length ? `Not in this layer: ${unknown.join(", ")}.` : "",
           failed.length ? `Could not load: ${failed.join(", ")}.` : "",
           (I.extra || []).filter((e) => e.unplaced).map((e) => `${e.name}: ${e.unplaced} points not matched to a bin.`).join(" "),
          ].filter(Boolean).join("\n"));
      $("volcanoPanel").hidden = false;
      draw();
    } catch (e) {
      msg(e.message || String(e), true);
    } finally { running = false; syncButtons(); }
  }

  // ---- drawing ----------------------------------------------------------------------
  let geom = null;                                  // screen positions of the drawn genes, for hit tests
  const fmtP = (l10) => (l10 > -4 ? Math.pow(10, l10).toPrecision(3) : `${Math.pow(10, l10 - Math.floor(l10)).toFixed(2)}e${Math.floor(l10)}`);
  const cls = (j) => (R.log10padj[j] <= Math.log10(V.padj) && Math.abs(R.lfc[j]) >= V.lfc ? (R.lfc[j] > 0 ? "up" : "down") : "ns");

  function draw() {
    if (!R || $("volcanoPanel").hidden) return;
    const host = $("vcHost"), W = Math.max(240, host.clientWidth), H = Math.max(200, host.clientHeight);
    const svg = $("vcSvg");
    svg.setAttribute("viewBox", `0 0 ${W} ${H}`); svg.setAttribute("width", W); svg.setAttribute("height", H);
    const m = { l: 46, r: 14, t: 12, b: 34 }, pw = W - m.l - m.r, ph = H - m.t - m.b;
    const js = []; for (let j = 0; j < R.genes.length; j++) if (R.tested[j]) js.push(j);
    // x: symmetric, to the 99.5th percentile of |log2FC| (never under 1.5 x the threshold); beyond it, triangles at the edge
    const ab = js.map((j) => Math.abs(R.lfc[j])).sort((p, q) => p - q);
    const xMax = Math.max(V.lfc * 1.5, 1, (ab[Math.floor(0.995 * (ab.length - 1))] || 1) * 1.1);
    const yv = (j) => -R.log10padj[j];
    let yTop = 1; for (const j of js) yTop = Math.max(yTop, yv(j));
    const sqrtY = V.yScale === "sqrt" || (V.yScale === "auto" && yTop > 50);
    const fy = sqrtY ? Math.sqrt : (v) => v, yMax = fy(yTop) * 1.05;
    const X = (v) => m.l + ((Math.max(-xMax, Math.min(xMax, v)) + xMax) / (2 * xMax)) * pw;
    const Y = (v) => m.t + ph - (fy(Math.max(0, v)) / yMax) * ph;
    const parts = [];
    // axes and ticks
    const tick = (lo, hi, n) => { const step = niceStep((hi - lo) / n), out = []; for (let v = Math.ceil(lo / step) * step; v <= hi + 1e-9; v += step) out.push(+v.toFixed(10)); return out; };
    parts.push(`<rect x="${m.l}" y="${m.t}" width="${pw}" height="${ph}" class="vc-frame"/>`);
    for (const v of tick(-xMax, xMax, 6)) parts.push(`<line class="vc-grid" x1="${X(v)}" x2="${X(v)}" y1="${m.t}" y2="${m.t + ph}"/><text class="vc-tick" x="${X(v)}" y="${m.t + ph + 14}" text-anchor="middle">${v}</text>`);
    const yTicks = sqrtY ? sqrtTicks(yTop) : tick(0, yTop, 5);
    for (const v of yTicks) parts.push(`<line class="vc-grid" x1="${m.l}" x2="${m.l + pw}" y1="${Y(v)}" y2="${Y(v)}"/><text class="vc-tick" x="${m.l - 5}" y="${Y(v) + 4}" text-anchor="end">${v}</text>`);
    parts.push(`<text class="vc-axis" x="${m.l + pw / 2}" y="${H - 4}" text-anchor="middle">log2 fold change  (← higher in B · higher in A →)</text>`);
    parts.push(`<text class="vc-axis" transform="translate(12 ${m.t + ph / 2}) rotate(-90)" text-anchor="middle">−log10 adjusted p${sqrtY ? " (√ scale)" : ""}</text>`);
    // thresholds
    const yThr = -Math.log10(V.padj);
    if (yThr <= yTop) parts.push(`<line class="vc-thr" x1="${m.l}" x2="${m.l + pw}" y1="${Y(yThr)}" y2="${Y(yThr)}"/>`);
    for (const s of [-1, 1]) if (V.lfc > 0 && V.lfc < xMax) parts.push(`<line class="vc-thr" x1="${X(s * V.lfc)}" x2="${X(s * V.lfc)}" y1="${m.t}" y2="${m.t + ph}"/>`);
    // points: grey first, then the significant ones on top
    let nUp = 0, nDown = 0;
    geom = [];
    const order = js.slice().sort((p, q) => (cls(p) === "ns") - (cls(q) === "ns") || yv(p) - yv(q)).reverse();
    for (const j of order) {
      const c = cls(j), x = X(R.lfc[j]), y = Y(yv(j)), out = Math.abs(R.lfc[j]) > xMax;
      if (c === "up") nUp++; else if (c === "down") nDown++;
      geom.push([j, x, y]);
      parts.push(out ? `<path class="vc-pt ${c}" d="M${x} ${y - 4}l3.5 7h-7z"/>` : `<circle class="vc-pt ${c}" cx="${x.toFixed(1)}" cy="${y.toFixed(1)}" r="${c === "ns" ? 2.2 : 3}"/>`);
    }
    // labels: the most significant genes past both thresholds, per side, nudged apart
    for (const side of ["up", "down"]) {
      const top = js.filter((j) => cls(j) === side).sort((p, q) => R.log10padj[p] - R.log10padj[q] || Math.abs(R.lfc[q]) - Math.abs(R.lfc[p])).slice(0, DEFAULTS.labels);
      const placed = [];
      for (const j of top) {
        let x = X(R.lfc[j]), y = Y(yv(j)) - 6;
        const anchor = x > m.l + pw - 60 ? "end" : x < m.l + 60 ? "start" : side === "up" ? "start" : "end";
        x += anchor === "start" ? 5 : -5;
        for (let k = 0; k < 20 && placed.some((p) => Math.abs(p[1] - y) < 11 && Math.abs(p[0] - x) < 70); k++) y += 11;
        y = Math.max(m.t + 9, Math.min(m.t + ph - 2, y));
        placed.push([x, y]);
        parts.push(`<text class="vc-lab ${side}" x="${x.toFixed(1)}" y="${y.toFixed(1)}" text-anchor="${anchor}">${esc(R.genes[j])}</text>`);
      }
    }
    // nothing past the thresholds: still name the few most significant, in grey
    if (!nUp && !nDown) {
      const placed = [];
      for (const j of js.slice().sort((p, q) => R.log10padj[p] - R.log10padj[q]).slice(0, 5)) {
        const x = X(R.lfc[j]), anchor = x > m.l + pw / 2 ? "end" : "start";
        let y = Y(yv(j)) - 6;
        for (let k = 0; k < 20 && placed.some((p) => Math.abs(p - y) < 11); k++) y += 11;
        placed.push(y);
        parts.push(`<text class="vc-lab ns" x="${(x + (anchor === "start" ? 5 : -5)).toFixed(1)}" y="${Math.max(m.t + 9, y).toFixed(1)}" text-anchor="${anchor}">${esc(R.genes[j])}</text>`);
      }
    }
    parts.push(`<g id="vcHi"></g>`);
    svg.innerHTML = parts.join("");
    highlight();
    $("vcTitle").textContent = `${R.nameA}  vs  ${R.nameB}`;
    $("vcCounts").innerHTML = `<span style="color:${COL.up}">↑ A ${nUp}</span> · <span style="color:${COL.down}">↑ B ${nDown}</span> · tested ${R.info.nTested} of ${R.info.nGenes}`;
    const I = R.info, warns = [];
    if (R.sig !== sigNow(R.a, R.b)) warns.push("The groups have changed since this volcano was computed: run it again.");
    if (Math.min(I.nA, I.nB) < 30) warns.push(`Few bins (${Math.min(I.nA, I.nB)} in the smaller group): treat the p-values as rough.`);
    // groups that came with the data (an annotation column, a lab's own clustering) were not found from these counts here
    if (R.a.startsWith("c:") || R.b.startsWith("c:") || ((R.a.startsWith("k:") || R.b.startsWith("k:")) && !R.givenGroups))
      warns.push("Clusters were found from these same genes, so testing them against the rest is circular: p-values are inflated. Use the ranking, not the p-values.");
    $("vcWarn").innerHTML = warns.map((w) => esc(say(w))).join("<br>"); $("vcWarn").hidden = !warns.length;
    $("vcFoot").textContent = say(
      `${R.test === "welch" ? "Welch t-test" : "Wilcoxon rank-sum (tie-corrected)"} on log1p CP10k (library = all genes), bins with ≥ ${I.minUmi} UMI; ` +
      `A ${I.nA.toLocaleString()} bins, B ${I.nB.toLocaleString()} bins. Tested: genes counted in ≥ ${+(I.minPct * 100).toFixed(2)}% of the bins of A or B, ` +
      `of the section's ${R.nHvg.toLocaleString()} variable genes${I.nExtra ? ` plus ${I.nExtra} named` : ""}; Benjamini-Hochberg over the ${I.nTested.toLocaleString()} tested. ` +
      `log2FC = log2((mean CP10k A + 1) / (mean CP10k B + 1)). ` +
      (Math.max(I.nA, I.nB) > 5000 ? "With this many bins nearly every difference is 'significant': bins are neighbours, not independent samples, so read the fold change and % expressing, and confirm on other sections. " : "Bins are neighbours, not independent samples: the p-values describe this section only. ") +
      `Hover a dot for its numbers; click it to show the gene on the tissue.`);
  }
  function niceStep(raw) { const p = Math.pow(10, Math.floor(Math.log10(raw || 1))), f = raw / p; return (f < 1.5 ? 1 : f < 3 ? 2 : f < 7 ? 5 : 10) * p; }
  function sqrtTicks(top) { const out = [0]; for (const v of [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000]) if (v <= top) out.push(v); return out.length > 7 ? out.filter((_, i) => i % 2 === 0) : out; }

  function hitGene(ev) {
    if (!geom || !R) return -1;
    const r = $("vcSvg").getBoundingClientRect(), x = ev.clientX - r.left, y = ev.clientY - r.top;
    let best = -1, bd = 7 * 7;
    for (const [j, gx, gy] of geom) { const d = (gx - x) ** 2 + (gy - y) ** 2; if (d < bd) { bd = d; best = j; } }
    return best;
  }
  function highlight() {
    const g = $("vcHi"); if (!g || !geom) return;
    const at = new Map(geom.map((p) => [p[0], p]));
    g.innerHTML = [pickedJ, hoverJ].filter((j) => j >= 0 && at.has(j)).map((j) => `<circle class="vc-ring" cx="${at.get(j)[1]}" cy="${at.get(j)[2]}" r="6"/>`).join("");
  }
  function tip(ev) {
    const el = $("vcTip");
    if (!ev || hoverJ < 0 || !R) { el.hidden = true; return; }
    const j = hoverJ, r = $("vcHost").getBoundingClientRect();
    el.innerHTML = `<b>${esc(R.genes[j])}</b>${j >= R.nHvg ? " <small>(named gene)</small>" : ""}<br>` + say(
      `log2FC ${R.lfc[j] >= 0 ? "+" : ""}${R.lfc[j].toFixed(2)} <small>(scanpy-style ${Number.isFinite(R.lfcScanpy[j]) ? R.lfcScanpy[j].toFixed(2) : "–"})</small><br>` +
      `padj ${fmtP(R.log10padj[j])} · p ${fmtP(R.log10p[j])}<br>` +
      `in A ${(R.pctA[j] * 100).toFixed(1)}% of bins, mean ${R.cpA[j].toFixed(1)} CP10k<br>` +
      `in B ${(R.pctB[j] * 100).toFixed(1)}% of bins, mean ${R.cpB[j].toFixed(1)} CP10k<br><small>click: show on the tissue</small>`);
    el.hidden = false;
    const x = ev.clientX - r.left, y = ev.clientY - r.top;
    el.style.left = `${Math.min(x + 12, r.width - el.offsetWidth - 4)}px`;
    el.style.top = `${y + 12 + el.offsetHeight > r.height ? y - el.offsetHeight - 8 : y + 12}px`;
  }

  // ---- CSV: one row per gene; p-values as text so 1e-400 survives --------------------------
  function saveCsv() {
    if (!R) return;
    const q = (t) => `"${String(t).replace(/"/g, '""')}"`, f = (x, d) => (Number.isFinite(x) ? (+x).toFixed(d) : "");
    const pTxt = (l10) => (l10 > -300 ? Math.pow(10, l10).toExponential(4) : `${Math.pow(10, l10 - Math.floor(l10)).toFixed(4)}e${Math.floor(l10)}`);
    const I = R.info, lines = ["dataset,layer,group_a,group_b,test,bins_a,bins_b,gene,source,tested,log2fc,log2fc_scanpy,statistic,pval,padj,neg_log10_padj,pct_a,pct_b,mean_cp10k_a,mean_cp10k_b"];
    const ord = R.genes.map((_, j) => j).sort((p, q2) => R.tested[q2] - R.tested[p] || R.log10padj[p] - R.log10padj[q2] || Math.abs(R.lfc[q2]) - Math.abs(R.lfc[p]));
    for (const j of ord) {
      const t = R.tested[j] === 1;
      lines.push([q(R.dataset), R.layer, q(R.nameA), q(R.nameB), R.test, I.nA, I.nB, q(R.genes[j]), j < R.nHvg ? "hvg" : "named", t ? 1 : 0,
        f(R.lfc[j], 4), f(R.lfcScanpy[j], 4), t ? f(R.stat[j], 4) : "", t ? pTxt(R.log10p[j]) : "", t ? pTxt(R.log10padj[j]) : "", t ? f(-R.log10padj[j], 3) : "",
        f(R.pctA[j], 4), f(R.pctB[j], 4), f(R.cpA[j], 3), f(R.cpB[j], 3)].join(","));
    }
    const slug = (s) => String(s).replace(/[^A-Za-z0-9]+/g, "_").replace(/^_|_$/g, "").slice(0, 40);
    api.download(`${R.dataset}_${R.layer}_volcano_${slug(R.nameA)}_vs_${slug(R.nameB)}.csv`, lines.join("\n") + "\n", "text/csv");
  }

  window.SodaVolcano = {
    init(a) { api = a; buildDom(); refresh(); },
    refresh,
    // a new layer or dataset: the old result belongs to other bins
    reset() { cancel(); R = null; lastGroups = ""; window.SodaVolcano.result = null; if ($("volcanoPanel")) $("volcanoPanel").hidden = true; refresh(); },
    result: null,
  };
})();
