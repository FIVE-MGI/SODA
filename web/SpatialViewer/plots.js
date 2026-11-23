// Heatmap, dot plot and violin plots of genes across groups of bins (or spots, or cells).
// Groups: the clusters of any section clustering in umap/index.json, the drawn shapes, the
// clusters of the last Circle -> UMAP run, the regions of region.bin, and for single-cell
// datasets any categorical annotation. Genes: the top markers of each group (where the
// grouping has marker lists) plus the names typed. The per-gene files are fetched once and
// kept; the numbers come from plots_worker.js; the pictures are plain SVG, so Save SVG is the
// picture itself and Save PNG is that SVG drawn on a canvas.
// viewer.js hands over what it owns through SodaPlots.init(api) and calls SodaPlots.refresh()
// when shapes, runs or the layer change; everything else lives here.
// Only unsigned typed arrays here: the de-identification check reads the signed names as IDs.
(function () {
  const NONE = 0xffff;
  const D = { minUmi: 10, topN: 3, maxTop: 20, maxGenes: 150, maxGenesTouch: 40, softMB: 150, softMBTouch: 30, hardMB: 600, hardMBTouch: 80,
              cacheMB: 400, cacheMBTouch: 60, fetchers: 4, maxGroups: 200, fewBins: 10, kdeGrid: 64, kdeMin: 5, points: 300, pointsTouch: 120, pointsAll: 6000,
              maxViolin: 4, chunkMB: 80, chunkMBTouch: 15, zClip: 2.5, seed: 42, worker: "plots_worker.js?v=20261116q3" };
  const THEMES = {
    dark: { bg: "#161b24", ink: "#e6e9ef", muted: "#9aa3b2", line: "#3a4352", grid: "#232a37", warn: "#f6d68a", hi: "#ffffff",
            seq: ["#1b2230", "#213a57", "#2c5b88", "#3f83bd", "#73afe6", "#cfe6ff"], div: ["#4ea1ff", "#39404d", "#ff5a5f"] },
    light: { bg: "#ffffff", ink: "#1a1d23", muted: "#5b6472", line: "#b9c0cc", grid: "#e9ecf1", warn: "#a15c00", hi: "#000000",
             seq: ["#f4f8fc", "#d3e3f3", "#a3c7e6", "#659fd3", "#2f71b3", "#0c3f7d"], div: ["#2166ac", "#f1f1f1", "#b2182b"] },
  };
  const FONT = "system-ui, -apple-system, 'Segoe UI', Roboto, sans-serif";
  let api = null, worker = null, running = false, runId = 0, R = null, lastSvg = "", lastGroups = "", uid = 0;
  const V = { view: "heat", scale: "mean", orderG: "listed", orderK: "listed", light: false, exprOnly: false, points: true, vGenes: [], sel: null };
  let cache = blankCache("");
  function blankCache(key) { return { key, genes: new Map(), bytes: 0, totals: null, chunkTot: null, chunkDone: new Set(), pos: null, secIndex: null, clusters: new Map(), obs: new Map() }; }
  const $ = (id) => document.getElementById(id);
  const esc = (x) => String(x).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
  const clamp = (x, a, b) => Math.max(a, Math.min(b, x));
  const nf = (x) => Number(x).toLocaleString("en-US");

  // ---- what is open ---------------------------------------------------------------------------
  const curExpr = (S) => ((S.m && S.m.expr) || []).find((e) => e.key === S.exprKey) || null;
  function modeNow() {
    const S = api.state();
    if (!S || !S.m || !S.src) return null;
    if (S.mode === "rna") return S.rna && S.src.buf ? "rna" : null;
    const e = curExpr(S);
    return e && e.kind !== "protein" && e.modality !== "protein" && S.src.gene ? "spatial" : null;
  }
  // cells of a single-cell dataset hold stored values, bins hold counts
  const words = () => (modeNow() === "rna" ? { bin: "cell", bins: "cells", above: "a value above 0", indep: "Cells of one sample are not independent replicates" }
    : (curExpr(api.state()) || {}).kind === "rois"      // GeoMx ROI segments (roi.js)
    ? { bin: "ROI", bins: "ROIs", above: "a count above 0", indep: "ROIs of one slide or animal are not independent replicates" }
    : (curExpr(api.state()) || {}).kind === "cells"     // segmented cells (Xenium, CosMx, MERSCOPE)
    ? { bin: "cell", bins: "cells", above: "a count above 0", indep: "Neighbouring cells are not independent samples" }
    : { bin: "bin", bins: "bins", above: "a count above 0", indep: "Neighbouring bins are not independent samples" });
  const cacheKey = () => { const S = api.state(); return `${S.dsId}|${S.mode === "rna" ? "rna" : S.exprKey}`; };
  // a clustering is known by its "id" when it has one (a lab's own), else by its Leiden resolution
  const clKey = (c) => (c.id != null ? c.id : c.resolution);
  const clLabel = (c) => c.name || `Leiden ${c.resolution}`;
  const rgbCss = (a) => `rgb(${a[0]},${a[1]},${a[2]})`;

  // ---- DOM: a panel in the sidebar, a floating panel over the map ------------------------------
  function buildDom() {
    const sec = document.createElement("section");
    sec.id = "panelPlots"; sec.className = "panel pl-side"; sec.hidden = true;
    sec.dataset.tab = "analyze";                                 // layout.js: which tab this section belongs to
    sec.innerHTML = `
      <h2>Heatmap and violin</h2>
      <label class="row slim"><span title="Which groups to compare: every group is one column of the heatmap and one violin">Group by</span><select id="plGroup"></select></label>
      ${window.SodaRoi ? window.SodaRoi.normRow("plNorm") : ""}
      <label class="row slim" id="plTopRow"><span title="The first N marker genes of each group, from the grouping's own marker lists (0 = none)">Markers / group</span>
        <input id="plTop" type="text" inputmode="numeric" value="${D.topN}" /></label>
      <label class="row slim"><span title="More genes by name, comma or space separated">Genes</span>
        <input id="plGenes" type="text" placeholder="e.g. GAPDH, ACTB" autocomplete="off" /></label>
      <div class="btnrow"><button id="plRun">Draw plots</button><button id="plShown" title="Add the genes now shown on the picture to the list">+ shown genes</button>
        <button id="plCancel" hidden>Cancel</button></div>
      <div id="plMsg" class="meta pl-msg"></div>`;
    const before = $("panelOverlays") || null, side = $("sidebar") || document.body;
    if (before && before.parentElement === side) side.insertBefore(sec, before); else side.appendChild(sec);
    if (window.SodaLayout) window.SodaLayout.adopt(sec);

    const panel = document.createElement("div");
    panel.id = "plotsPanel"; panel.className = "pl-panel"; panel.hidden = true;
    panel.innerHTML = `
      <div class="pl-head"><b id="plTitle">Plots</b>
        <span class="pl-btns"><button id="plPng" title="The picture as shown, as a PNG image (2x)">Save PNG</button><button id="plSvg" title="The picture as shown, as an SVG (vector) file">Save SVG</button><button id="plCsv" title="The numbers plotted, one row per gene and group, with the number of bins of every group">Save CSV</button><button id="plClose" title="Hide">×</button></span></div>
      <div id="plWarn" class="lp-warn" hidden></div>
      <div class="pl-ctl">
        <span class="seg pl-tabs" id="plTabs"><button data-v="heat" class="on">Heatmap</button><button data-v="dot">Dot plot</button><button data-v="violin">Violin</button></span>
        <label class="pl-m" title="Colour: the plain mean, or each gene rescaled across the groups shown">Colour <select id="plScale">
          <option value="mean">mean expression</option><option value="z">z-score per gene</option><option value="mm">0–1 per gene</option></select></label>
        <label class="pl-m" title="Order of the genes (rows)">Genes <select id="plOrderG"><option value="listed">as listed</option><option value="cluster">clustered</option><option value="az">A–Z</option></select></label>
        <label title="Order of the groups">Groups <select id="plOrderK"><option value="listed">as listed</option><option value="size">by size</option><option value="cluster">clustered</option></select></label>
        <label class="pl-v" title="Leave out the bins (or cells) where the gene was not counted: the violin then shows only the expressing ones"><input id="plExpr" type="checkbox" /> expressing only</label>
        <label class="pl-v" id="lblPlPts" title="A seeded sample of the bins behind each violin"><input id="plPts" type="checkbox" checked /> points</label>
        <label title="White background, for print"><input id="plLight" type="checkbox" /> light</label>
      </div>
      <div class="pl-ctl pl-v" id="plVRow"><span id="plVChips" class="pl-chips"></span><select id="plVAdd"></select></div>
      <div id="plSel" class="pl-sel"></div>
      <div id="plHost" class="pl-host"><div id="plSvgBox" class="pl-svgbox"></div><div id="plTip" class="pl-tip" hidden></div></div>
      <div id="plFoot" class="lp-foot"></div>`;
    (document.querySelector(".main") || document.body).appendChild(panel);

    $("plRun").onclick = runOrShow;
    $("plCancel").onclick = () => { cancel(); msg("Cancelled."); };
    $("plShown").onclick = () => {
      const S = api.state(), shown = (modeNow() === "rna" ? (S.rna && S.rna.genes) || [] : S.genes || []).map((g) => g.name);
      const have = new Set($("plGenes").value.split(/[\s,;]+/).filter(Boolean).map((x) => x.toUpperCase()));
      const add = shown.filter((n) => !have.has(n.toUpperCase()));
      if (!shown.length) { msg("No gene is shown on the picture yet: search one above first, or type names here."); return; }
      $("plGenes").value = [$("plGenes").value.trim().replace(/[,;\s]+$/, ""), ...add].filter(Boolean).join(", ");
    };
    $("plGroup").addEventListener("pointerdown", refresh);       // clusterings can change without telling this panel
    $("plGroup").onchange = () => { syncTopRow(); syncButtons(); };
    for (const id of ["plTop", "plGenes"]) $(id).addEventListener("input", syncButtons);
    if (window.SodaRoi && $("plNorm")) { window.SodaRoi.bindNorm("plNorm"); $("plNorm").addEventListener("change", () => { if (R && !$("plotsPanel").hidden) syncWarn(); }); }
    $("plGenes").addEventListener("keydown", (e) => { if (e.key === "Enter") runOrShow(); });
    $("plClose").onclick = () => { $("plotsPanel").hidden = true; syncButtons(); };
    $("plTabs").onclick = (e) => { const v = e.target.dataset && e.target.dataset.v; if (v) setView(v); };
    $("plScale").onchange = () => { V.scale = $("plScale").value; draw(); };
    $("plOrderG").onchange = () => { V.orderG = $("plOrderG").value; draw(); };
    $("plOrderK").onchange = () => { V.orderK = $("plOrderK").value; draw(); };
    $("plExpr").onchange = () => { V.exprOnly = $("plExpr").checked; draw(); };
    $("plPts").onchange = () => { V.points = $("plPts").checked; draw(); };
    $("plLight").onchange = () => { V.light = $("plLight").checked; draw(); };
    $("plVAdd").onchange = () => { const j = Number($("plVAdd").value); $("plVAdd").value = ""; if (Number.isInteger(j) && j >= 0) addViolinGene(j); };
    $("plVChips").onclick = (e) => { const j = e.target.dataset && e.target.dataset.rm; if (j !== undefined) { V.vGenes = V.vGenes.filter((x) => x !== Number(j)); draw(); } };
    $("plSvg").onclick = saveSvg; $("plPng").onclick = savePng; $("plCsv").onclick = saveCsv;
    $("plSel").onclick = (e) => { const a = e.target.dataset && e.target.dataset.act; if (a) { e.preventDefault(); selAction(a); } };
    const box = $("plSvgBox");
    box.addEventListener("mousemove", (ev) => tip(ev));
    box.addEventListener("mouseleave", () => tip(null));
    box.addEventListener("click", (ev) => {
      const t = hit(ev); if (!t || !R) return;
      V.sel = t; tip(ev); syncSel(); markSel();
    });
    // the matrix follows the width alone, so the selection bar growing by a line after a click does not rebuild the picture under the pointer
    if (window.ResizeObserver) { let last = ""; new ResizeObserver(() => { const h = $("plHost"), k = V.view === "violin" ? `${h.clientWidth}x${h.clientHeight}` : `${h.clientWidth}`; if (k !== last && h.clientWidth && !$("plotsPanel").hidden) { last = k; draw(); } }).observe($("plHost")); }
  }

  const msg = (t, err) => { const el = $("plMsg"); if (!el) return; el.textContent = t; el.style.color = err ? "#ff8f8f" : ""; };
  // what plots are made from: the grouping (and its content), the markers per group, the genes typed
  const inputsNow = () => [sigOf($("plGroup").value), $("plTop").disabled ? "0" : $("plTop").value.trim(), $("plGenes").value.trim(), modeNow()].join("|");
  const canShow = () => !!R && !running && $("plotsPanel").hidden && R.inputs === inputsNow();
  // one button: plots of exactly these inputs are there but closed, so it shows them again; anything changed, it draws new ones
  function runOrShow() { if (canShow()) { $("plotsPanel").hidden = false; syncButtons(); draw(); } else runPlot(); }
  function syncButtons() {
    if (!$("plRun")) return;
    const show = canShow();
    $("plRun").disabled = running || (!show && !$("plGroup").value);
    $("plRun").textContent = show ? "Show plots" : "Draw plots";
    $("plRun").title = show ? "Show the last plots again (nothing has changed since they were drawn)" : "Draw the genes across the groups";
    $("plCancel").hidden = !running;
  }
  function syncTopRow() {
    const v = $("plGroup").value, has = v.startsWith("sec:") || v === "run";
    $("plTop").disabled = !has;
    $("plTopRow").title = has ? "" : "This grouping has no marker lists: type the genes";
    syncButtons();
  }

  // ---- the groupings on offer -------------------------------------------------------------------
  function groupingsNow() {
    const S = api.state(), out = [], mode = modeNow();
    if (mode === "rna") {
      for (const o of S.rna.obs || []) if (o.type === "category") out.push({ v: `obs:${o.name}`, t: `${o.name} · ${o.categories.length} groups`, og: "Annotations of the cells" });
      return out;
    }
    if (mode !== "spatial") return out;
    if (S.src.buf && S.src.json) for (const c of api.clusterings() || []) out.push({ v: `sec:${clKey(c)}`, t: `${clLabel(c)} · ${c.n} ${c.source === "user" ? "groups" : "clusters"}`, og: "Clusters of the whole section" });
    const L = S.lasso;
    if (L && L.assign && L.pos) {
      const live = L.shapes.filter((s) => s.bins > 0).length;
      if (live) out.push({ v: "shapes", t: `The drawn shapes (${live})`, og: "Drawn on the slide" });
      const Rr = L.result;
      if (Rr && Rr.cluster && Rr.info.mode === "bins" && Rr.layer === S.exprKey && Rr.info.nClusters > 0) out.push({ v: "run", t: `Clusters of the Draw a region run (${Rr.info.nClusters})`, og: "Drawn on the slide" });
    }
    if (S.region && S.region.meta && S.region.ids && S.region.meta.length) out.push({ v: "reg", t: `Regions (${S.region.meta.length})`, og: "Regions of the overlay" });
    return out;
  }
  function refresh() {
    if (!api || !$("panelPlots")) return;
    const list = groupingsNow(), key = JSON.stringify(list) + "|" + (modeNow() || "");
    $("panelPlots").hidden = !list.length;
    if (window.SodaRoi) window.SodaRoi.normShown("plNorm", modeNow() === "spatial" ? curExpr(api.state()) : null);   // GeoMx: Q3 or per 10,000
    if (!list.length && !$("plotsPanel").hidden) $("plotsPanel").hidden = true;
    if (key !== lastGroups) {
      lastGroups = key;
      const sel = $("plGroup"), old = sel.value;
      sel.innerHTML = "";
      let og = null;
      for (const o of list) {
        if (!og || og.label !== o.og) { og = document.createElement("optgroup"); og.label = o.og; sel.appendChild(og); }
        const e = document.createElement("option"); e.value = o.v; e.textContent = o.t; og.appendChild(e);
      }
      if (list.some((o) => o.v === old)) sel.value = old;
      else if (list.length) {
        // the section's default clustering when there is one
        const S = api.state(), e = modeNow() === "spatial" ? curExpr(S) : null, cls = (modeNow() === "spatial" && api.clusterings()) || [];
        const def = cls.length ? cls[Math.min(cls.length - 1, (e && e.umap && e.umap.default) || 0)] : null;
        const want = def ? `sec:${clKey(def)}` : list[0].v;
        sel.value = list.some((o) => o.v === want) ? want : list[0].v;
      }
      syncTopRow();
    }
    if (R && !$("plotsPanel").hidden) syncWarn();
    syncButtons();
  }

  // what a result depends on: dataset, layer, the grouping's content, the QC filter
  function sigOf(gv) {
    const S = api.state(), L = S.lasso, qc = api.qc && api.qc();
    let part = "";
    if (gv === "shapes") part = L ? L.shapes.map((s) => `${s.id}:${s.bins}`).join(",") : "";
    else if (gv === "run") part = L && L.result ? `${L.result.id}` : "";
    const q3 = modeNow() === "spatial" && window.SodaRoi && window.SodaRoi.scale(curExpr(S)) ? "q3" : "";
    return [cacheKey(), gv, part, modeNow() === "spatial" && qc ? JSON.stringify(qc.state()) : "", q3].join("|");
  }
  function staleWarn() {
    const w = [];
    if (R && R.sig !== sigOf(R.grouping.id)) w.push(R.grouping.id === "shapes" || R.grouping.id === "run" ? "The shapes, the run or the QC filter have changed since these plots were computed: draw them again."
      : window.SodaRoi && R.mode === "spatial" && window.SodaRoi.q3(curExpr(api.state())) ? "The QC filter or the normalisation has changed since these plots were computed: draw them again."
      : "The QC filter has changed since these plots were computed: draw them again.");
    return w;
  }

  // ---- groups: a code per row -------------------------------------------------------------------
  async function positions() {
    const S = api.state(), L = S.lasso;
    if (cache.pos) return cache.pos;
    const p = L && L.pos ? L.pos : S.tissue && S.tissue.positions && S.tissue.key === S.exprKey ? S.tissue.positions : await S.src.positions(S.exprKey);
    cache.pos = p;
    return p;
  }
  async function sectionIndex() {
    const S = api.state();
    if (S.umap && S.umap.key === S.exprKey && S.umap.index) return S.umap.index;
    if (!cache.secIndex) cache.secIndex = await S.src.json(`expr/${S.exprKey}/umap/index.json`);
    return cache.secIndex;
  }
  async function buildGroups(gv) {
    const S = api.state(), L = S.lasso;
    if (gv.startsWith("obs:")) {
      const o = (S.rna.obs || []).find((x) => x.name === gv.slice(4));
      if (!o) throw new Error("That annotation is not in this dataset.");
      let codes = cache.obs.get(o.name);
      if (!codes) { codes = new Uint16Array(await S.src.buf(`obs/${o.file}`)); cache.obs.set(o.name, codes); }
      if (codes.length !== S.rna.n) throw new Error(`The annotation has ${codes.length} cells, the dataset ${S.rna.n}.`);
      return { id: gv, kind: "obs", label: o.name, K: o.categories.length, N: S.rna.n, group: codes, names: o.categories.map(String),
               colors: o.categories.map((_, i) => (o.colors && o.colors[i]) || rgbCss(api.catRgb(i))), markers: null };
    }
    const pos = await positions(), N = pos.length / 2;
    if (gv.startsWith("sec:")) {
      // the viewer's own description first (your own results' groupings are only there)
      const own0 = (api.clusterings() || []).find((x) => String(clKey(x)) === gv.slice(4));
      const idx = own0 ? null : await sectionIndex(), c = own0 || (idx.clusterings || []).find((x) => String(clKey(x)) === gv.slice(4));
      if (!c) throw new Error("That clustering is not in this layer.");
      const ck = clKey(c);
      let cl = api.clusterCodes ? await api.clusterCodes(c) : null;
      if (!cl) cl = S.umap && S.umap.key === S.exprKey && S.umap.clusters ? S.umap.clusters.get(ck) : null;
      if (!cl) cl = cache.clusters.get(ck);
      if (!cl) { cl = new Uint16Array(await S.src.buf(`expr/${S.exprKey}/umap/${c.file || `clusters_r${c.resolution}.bin`}`)); cache.clusters.set(ck, cl); }
      if (cl.length !== N) throw new Error(`The clusters file has ${cl.length} rows, the layer ${N}.`);
      const own = (S.umap && S.umap.names && S.umap.names[ck]) || {};
      const names = [], colors = [];
      for (let k = 0; k < c.n; k++) {
        names.push(own[k] || (c.categories && c.categories[k] != null ? String(c.categories[k]) : `Cluster ${k}`));
        colors.push((Array.isArray(c.colors) && c.colors[k]) || rgbCss(api.clusterRgb(k)));
      }
      const markers = Array.isArray(c.markers) && c.markers.length ? c.markers.map((l) => (l || []).map((m) => m && m.gene).filter((g) => typeof g === "string")) : null;
      return { id: gv, kind: "sec", label: clLabel(c), K: c.n, N, group: cl, names, colors, markers };
    }
    if (gv === "shapes") {
      if (!L || !L.assign || L.assign.length !== N) throw new Error("The shapes are not ready yet.");
      return { id: gv, kind: "shapes", label: "Drawn shapes", K: L.shapes.length, N, group: L.assign.slice(), names: L.shapes.map((s) => s.name), colors: L.shapes.map((s) => s.hex), markers: null };
    }
    if (gv === "run") {
      const Rr = L && L.result;
      if (!Rr || !Rr.cluster) throw new Error("There is no Draw a region run with clusters.");
      const K = Rr.info.nClusters, g = new Uint16Array(N).fill(NONE);
      for (let j = 0; j < Rr.rows.length; j++) { const c = Rr.cluster[j]; if (c !== NONE && c < K) g[Rr.rows[j]] = c; }
      const names = [], colors = [];
      for (let k = 0; k < K; k++) { names.push((L.clNames && L.clNames[k]) || `Cluster ${k}`); colors.push(rgbCss(api.clusterRgb(k))); }
      const mk = Rr.info.markers, markers = Array.isArray(mk) && mk.length === K ? mk.map((l) => (l || []).map((m) => m && m.gene).filter((x) => typeof x === "string")) : null;
      return { id: gv, kind: "run", label: "Circle → UMAP clusters", K, N, group: g, names, colors, markers };
    }
    if (gv === "reg") {
      const meta = S.region.meta, ids = S.region.ids;
      if (!meta || !ids) throw new Error("This layer has no regions.");
      const at = new Map(meta.map((r, i) => [r.id, i])), g = new Uint16Array(N).fill(NONE);
      for (let r = 0; r < Math.min(N, ids.length); r++) { const i = at.get(ids[r]); if (i !== undefined) g[r] = i; }
      return { id: gv, kind: "reg", label: "Regions", K: meta.length, N, group: g, names: meta.map((r) => String(r.name)), colors: meta.map((_, i) => rgbCss(api.catRgb(i))), markers: null };
    }
    throw new Error(`unknown grouping ${gv}`);
  }

  // ---- genes ------------------------------------------------------------------------------------
  const geneIndex = () => { const S = api.state(); return (modeNow() === "rna" ? S.rna.index : S.geneIndex) || []; };
  function geneList(G, topN, typed) {
    const index = new Map(geneIndex().map((x) => [x.n.toUpperCase(), x])), genes = [], seen = new Set(), unknown = [], from = {};
    const add = (name, src) => {
      const e = index.get(String(name).toUpperCase());
      if (!e) { if (!unknown.includes(name)) unknown.push(name); return; }
      if (seen.has(e.n)) return;
      seen.add(e.n); genes.push(e.n); from[e.n] = src;
    };
    if (G.markers && topN > 0) G.markers.forEach((list, k) => list.slice(0, topN).forEach((g) => add(g, G.names[k] || `group ${k}`)));
    for (const w of typed.split(/[\s,;]+/).filter(Boolean)) add(w, "typed");
    return { genes, unknown, from, index };
  }
  // an upper bound on a gene file's size, from the gene list alone
  function estBytes(e, N) { return modeNow() === "rna" ? 4 + 8 * Math.ceil((e.p || 0) * N) : 12 * Math.min(e.t || 0, N); }

  function evict(needed) {
    const cap = (api.touch ? D.cacheMBTouch : D.cacheMB) * 1e6;
    if (cache.bytes <= cap) return;
    const old = [...cache.genes.entries()].filter(([n]) => !needed.has(n)).sort((a, b) => a[1].at - b[1].at);
    for (const [n, e] of old) { if (cache.bytes <= cap) break; cache.genes.delete(n); cache.bytes -= e.bytes; }
  }
  // names -> { bufs: Map name -> ArrayBuffer | null (no counts), failed, fetched: {files, bytes}, cached }
  async function fetchGenes(names, index, dead, onProg) {
    const S = api.state(), rna = modeNow() === "rna", key = S.exprKey, bufs = new Map(), failed = [], need = [];
    let cached = 0, files = 0, bytes = 0;
    for (const n of names) {
      const e = index.get(n.toUpperCase()), c = cache.genes.get(n);
      if (c) { c.at = Date.now(); bufs.set(n, c.buf); cached++; }
      else if (!rna && e && !(e.t > 0)) bufs.set(n, null);          // listed, but never counted: no file
      else need.push(n);
    }
    let next = 0;
    const one = async () => {
      while (next < need.length && !dead()) {
        const n = need[next++];
        try {
          const buf = rna ? await S.src.rnaGene(n) : (await S.src.gene(key, n)).buffer;
          if (dead()) return;
          cache.genes.set(n, { buf, bytes: buf.byteLength, at: Date.now() }); cache.bytes += buf.byteLength;
          bufs.set(n, buf); files++; bytes += buf.byteLength;
        } catch (err) { failed.push(n); }
        onProg(files + failed.length, need.length, bytes);
      }
    };
    await Promise.all(Array.from({ length: Math.min(D.fetchers, need.length) }, one));
    evict(new Set(names));
    return { bufs, failed, fetched: { files, bytes }, cached };
  }

  // ---- per-bin totals: binqc/total.bin, else the hvg/ chunks' library sizes, else none ---------
  async function loadTotals(group, dead) {
    const S = api.state(), e = curExpr(S), N = group.length, force = window.SodaPlots.force || "";
    const m = e && e.binQc && Array.isArray(e.binQc.metrics) ? e.binQc.metrics.find((x) => x.key === "total") : null;
    if (m && S.src.buf && force !== "chunks" && force !== "raw") {
      if (cache.totals) return { totals: cache.totals, from: "binqc", bytes: 0, note: "" };
      try {
        const buf = await S.src.buf(`expr/${S.exprKey}/${e.binQc.dir || "binqc"}/${m.file}`);
        const t = m.dtype === "uint32" ? Float32Array.from(new Uint32Array(buf)) : m.dtype === "uint16" ? Float32Array.from(new Uint16Array(buf)) : new Float32Array(buf);
        if (t.length === N) { cache.totals = t; return { totals: t, from: "binqc", bytes: buf.byteLength, note: "" }; }
      } catch (err) { /* fall through */ }
    }
    const L = S.lasso;
    if (L && L.hvg && L.pos && L.pos.length === 2 * N && !L.prot && api.fetchChunks && force !== "raw") {
      const cp = L.hvg.chunkPx, keys = new Set();
      for (let r = 0; r < N; r++) if (group[r] !== NONE) keys.add(api.chunkOf(L.pos[2 * r], L.pos[2 * r + 1], cp));
      const todo = [...keys].filter((k) => !cache.chunkDone.has(k));
      const mb = L.hvg.chunks.filter((c) => todo.includes(`${c.tx}_${c.ty}`)).reduce((s, c) => s + c.bytes, 0) / 1e6, lim = api.touch ? D.chunkMBTouch : D.chunkMB;
      if (mb <= lim) {
        if (todo.length) msg(`This export has no per-${words().bin} totals file: reading them from ${todo.length} data chunk(s) (${mb.toFixed(1)} MB)…`);
        const got = todo.length ? await api.fetchChunks(todo) : { bufs: [], missing: [] };
        if (dead()) return null;
        if (!got.missing.length) {
          if (!cache.chunkTot) cache.chunkTot = new Float32Array(N).fill(NaN);
          for (const buf of got.bufs) {
            const n = new DataView(buf).getUint32(0, true), row = new Uint32Array(buf, 4, n), lib = new Float32Array(buf, 4 + 4 * n, n);
            for (let i = 0; i < n; i++) if (row[i] < N) cache.chunkTot[row[i]] = lib[i];
          }
          todo.forEach((k) => cache.chunkDone.add(k));
          return { totals: cache.chunkTot, from: "chunks", bytes: mb * 1e6, note: `Per-${words().bin} totals read from the Draw a region data chunks (this export has no binqc/ totals).` };
        }
        return { totals: null, from: "none", bytes: 0, note: `${got.missing.length} data chunk(s) with the per-${words().bin} totals could not be loaded: showing log1p of the raw counts, NOT depth-normalised.` };
      }
      return { totals: null, from: "none", bytes: 0, note: `This export has no per-${words().bin} totals file and reading them from the data chunks would take ${mb.toFixed(0)} MB: showing log1p of the raw counts, NOT depth-normalised.` };
    }
    return { totals: null, from: "none", bytes: 0, note: `This export has no per-${words().bin} totals (older export): showing log1p of the raw counts, NOT depth-normalised. Re-export to get counts per 10k.` };
  }

  // ---- the worker ---------------------------------------------------------------------------------
  function cancel() {
    runId++;
    if (worker) { worker.terminate(); worker = null; }
    running = false; syncButtons();
  }
  function callWorker(payload, onProg) {
    if (!worker) worker = new Worker(D.worker);
    return new Promise((resolve, reject) => {
      worker.onmessage = (ev) => {
        const d = ev.data;
        if (d.type === "progress") { onProg(d.step, d.frac); return; }
        if (d.type === "error") reject(new Error(d.message)); else resolve(d);
      };
      worker.onerror = (e) => reject(new Error(e.message || "the plot worker failed"));
      worker.postMessage(payload);
    });
  }
  const unitOf = (norm, S) => (norm === "cp10k" ? "log1p(counts per 10k)" : norm === "q3" ? "log1p(Q3-normalised counts)" : norm === "raw" ? "log1p(UMI count), not depth-normalised" : (S.m.valueLabel || "expression (as provided)"));

  async function runPlot() {
    if (running || !api) return;
    const S = api.state(), mode = modeNow(), gv = $("plGroup").value, W = words(), inputs = inputsNow();
    if (!mode || !gv) { msg("Pick what to group by.", true); return; }
    const topN = $("plTop").disabled ? 0 : Number($("plTop").value);
    if (!Number.isInteger(topN) || topN < 0 || topN > D.maxTop) { msg(`Markers per group must be a whole number from 0 to ${D.maxTop}.`, true); return; }
    const my = ++runId, dead = () => my !== runId;
    running = true; syncButtons();
    const t0 = performance.now(), T = {};
    try {
      if (cache.key !== cacheKey()) cache = blankCache(cacheKey());
      msg("Finding the groups…");
      const sig = sigOf(gv), G = await buildGroups(gv);
      if (dead()) return;
      if (G.K > D.maxGroups) throw new Error(`${G.K} groups: more than ${D.maxGroups} cannot be drawn.`);
      if (G.K < 1) throw new Error("There is no group to plot.");
      const gl = geneList(G, topN, $("plGenes").value);
      if (!gl.genes.length) throw new Error(gl.unknown.length ? `Not in this ${mode === "rna" ? "dataset" : "layer"}: ${gl.unknown.join(", ")}.`
        : G.markers ? "No genes: set the markers per group above 0, or type gene names." : "This grouping has no marker lists: type gene names (or add the genes shown on the picture).");
      const maxG = api.touch ? D.maxGenesTouch : D.maxGenes;
      if (gl.genes.length > maxG) throw new Error(`${gl.genes.length} genes: at most ${maxG} at a time${api.touch ? " on a phone or tablet" : ""}. Lower the markers per group or shorten the list.`);
      // per-bin totals and the QC filter
      let totals = null, totalsFrom = "", totNote = "", totBytes = 0, norm = "asis", keep = null, qcState = null, qcTxt = "", scale = null;
      if (mode === "spatial") {
        const t = await loadTotals(G.group, dead);
        if (dead() || !t) return;
        totals = t.totals; totalsFrom = t.from; totNote = t.note; totBytes = t.bytes; norm = totals ? "cp10k" : "raw";
        // GeoMx ROI layer with Q3 factors and Q3 chosen (roi.js): count x 1 / normFactor per row
        const sc = window.SodaRoi && window.SodaRoi.scale(curExpr(S));
        if (sc && sc.length === G.N) { scale = sc; norm = "q3"; }
        const qc = api.qc && api.qc(), k = qc ? qc.keepMask() : null;
        if (k && k.length === G.N) { keep = k.slice(); qcState = qc.state(); qcTxt = qc.describe(); }
      }
      T.groups = performance.now() - t0;
      // what it will cost, before fetching anything
      let est = 0, toFetch = 0;
      for (const n of gl.genes) if (!cache.genes.has(n)) { const b = estBytes(gl.index.get(n.toUpperCase()), G.N); if (b > 0) { est += b; toFetch++; } }
      const soft = api.touch ? D.softMBTouch : D.softMB, hard = api.touch ? D.hardMBTouch : D.hardMB;
      if (est / 1e6 > hard) throw new Error(`These ${toFetch} genes need up to ${(est / 1e6).toFixed(0)} MB of downloads (limit ${hard} MB${api.touch ? " on a phone or tablet" : ""}). Use fewer genes or fewer markers per group.`);
      if (est / 1e6 > soft && !confirm(`These ${toFetch} genes need up to ${(est / 1e6).toFixed(0)} MB of downloads (one file per gene, each covering the whole ${mode === "rna" ? "dataset" : "section"}). Go on?`)) { msg("Not run."); return; }
      msg(toFetch ? `Fetching ${toFetch} gene file(s), up to ${(est / 1e6).toFixed(1)} MB…` : "Averaging…");
      const got = await fetchGenes(gl.genes, gl.index, dead, (i, n, b) => msg(`Fetching gene ${i} of ${n} (${(b / 1e6).toFixed(1)} MB so far)…`));
      if (dead()) return;
      T.fetch = performance.now() - t0 - T.groups;
      const genes = gl.genes.filter((n) => got.bufs.has(n));
      if (!genes.length) throw new Error(`None of the genes could be loaded (${got.failed.join(", ")}).`);
      const input = { layout: mode === "rna" ? "indexed" : "xyv", nRows: G.N, positions: mode === "rna" ? null : await positions(), group: G.group, nGroups: G.K, totals, keep, scale,
                      base: { norm, minUmi: D.minUmi, asisLog: mode === "rna" && /^log1p/i.test(S.m.valueLabel || ""), kdeGrid: D.kdeGrid, kdeMin: D.kdeMin,
                              maxPoints: Math.max(20, Math.min(api.touch ? D.pointsTouch : D.points, Math.floor(D.pointsAll / G.K))), seed: D.seed } };
      const res = await callWorker({ ...input, base: undefined, genes: genes.map((n) => ({ name: n, buf: got.bufs.get(n) })), params: { ...input.base, detail: false } },
        (step, f) => { if (!dead()) msg(`${step}… ${Math.round(f * 100)}%`); });
      if (dead()) return;
      T.worker = performance.now() - t0 - T.groups - T.fetch;
      const cols = []; for (let k = 0; k < G.K; k++) if (res.n[k] > 0) cols.push(k);
      if (!cols.length) throw new Error(`No ${W.bins} are left in any group` + (keep ? " after the QC filter" : "") + (norm === "cp10k" ? ` (${W.bins} under ${D.minUmi} UMI are left out).` : norm === "q3" ? ` (${W.bins} without a Q3 factor are left out).` : "."));
      const sum = (a) => a.reduce((s, x) => s + x, 0);
      R = { ...res, sig, inputs, input, dataset: S.dsId, layer: mode === "rna" ? "cells" : S.exprKey, mode, grouping: { id: G.id, kind: G.kind, label: G.label }, K: G.K, names: G.names, colors: G.colors, cols,
            genes, from: gl.from, noCounts: genes.filter((n) => got.bufs.get(n) === null), unit: unitOf(norm, S), norm, totalsFrom, qcState, qcTxt, words: W,
            detail: new Map(), nLow: sum(Array.from(res.low)), nQc: sum(Array.from(res.qcOut)), empty: G.names.filter((_, k) => res.found[k] > 0 && res.n[k] === 0),
            bytes: { genes: got.fetched.bytes, files: got.fetched.files, cached: got.cached, totals: totBytes, estimate: est }, secs: (performance.now() - t0) / 1000, ms: { ...T, workerInner: res.info.ms } };
      window.SodaPlots.result = R;
      V.sel = null; V.vGenes = [genes.find((n) => gl.from[n] === "typed" && !R.noCounts.includes(n)) || genes[0]].map((n) => genes.indexOf(n));
      const unplaced = res.info.genes.filter((g) => g.unplaced).map((g) => `${g.name}: ${nf(g.unplaced)} points not matched to a ${W.bin}`);
      msg([`${genes.length} gene${genes.length === 1 ? "" : "s"} × ${cols.length} group${cols.length === 1 ? "" : "s"} in ${R.secs.toFixed(1)} s ` +
             `(${got.fetched.files} file${got.fetched.files === 1 ? "" : "s"}, ${(got.fetched.bytes / 1e6).toFixed(1)} MB fetched${got.cached ? `; ${got.cached} already loaded` : ""}).`,
           totNote,
           R.nLow ? `Left out ${nf(R.nLow)} ${W.bins} with < ${D.minUmi} UMI${norm === "q3" ? " or no Q3 factor (Q3 = 0)" : ""}.` : "",
           R.nQc ? `QC filter: left out ${nf(R.nQc)} ${W.bins} (${qcTxt}).` : "",
           R.empty.length ? `No ${W.bins} left in: ${R.empty.join(", ")}.` : "",
           gl.unknown.length ? `Not in this ${mode === "rna" ? "dataset" : "layer"}: ${gl.unknown.join(", ")}.` : "",
           got.failed.length ? `Could not load: ${got.failed.join(", ")}.` : "",
           R.noCounts.length ? `No counts in this layer: ${R.noCounts.join(", ")}.` : "",
           unplaced.join("; "),
          ].filter(Boolean).join("\n"));
      $("plotsPanel").hidden = false;
      if (V.view === "violin") await ensureDetail(my); else draw();
    } catch (e) {
      if (!dead()) msg(e.message || String(e), true);
    } finally { if (!dead()) { running = false; syncButtons(); } }
  }

  // the violin's numbers for the genes on show (quartiles, density, points): a second, small run
  async function ensureDetail(owner) {
    if (!R) return;
    const need = V.vGenes.filter((j) => !R.detail.has(j));
    if (!need.length) { draw(); return; }
    if (running && owner === undefined) return;
    const my = owner !== undefined ? owner : ++runId, dead = () => my !== runId, Rr = R;
    running = true; syncButtons();
    try {
      const names = need.map((j) => Rr.genes[j]);
      const got = await fetchGenes(names, new Map(geneIndex().map((x) => [x.n.toUpperCase(), x])), dead, () => {});
      if (dead() || R !== Rr) return;
      const I = Rr.input;
      const res = await callWorker({ layout: I.layout, nRows: I.nRows, positions: I.positions, group: I.group, nGroups: I.nGroups, totals: I.totals, keep: I.keep, scale: I.scale,
                                     genes: names.map((n) => ({ name: n, buf: got.bufs.get(n) || null })), params: { ...I.base, detail: true } }, () => {});
      if (dead() || R !== Rr) return;
      need.forEach((j, i) => Rr.detail.set(j, res.detail[i]));
      Rr.kde = res.info.kde; Rr.msDetail = res.info.ms;
      draw();
    } catch (e) {
      if (!dead()) msg(e.message || String(e), true);
    } finally { if (!dead()) { running = false; syncButtons(); } }
  }
  function addViolinGene(j) {
    if (!R || V.vGenes.includes(j)) return;
    V.vGenes = V.vGenes.concat(j).slice(-D.maxViolin);
    if (V.view !== "violin") setView("violin"); else ensureDetail();
  }
  function setView(v) {
    V.view = v;
    $("plTabs").querySelectorAll("button").forEach((b) => b.classList.toggle("on", b.dataset.v === v));
    if (v === "violin") ensureDetail(); else draw();
  }

  // ---- ordering and scaling -----------------------------------------------------------------------
  // average-linkage clustering on 1 - Pearson r; returns the leaves left to right
  function clusterOrder(vecs) {
    const n = vecs.length;
    if (n < 3) return vecs.map((_, i) => i);
    const z = vecs.map((v) => {
      const f = Array.from(v, (x) => (Number.isFinite(x) ? x : 0)), m = f.reduce((s, x) => s + x, 0) / f.length;
      let q = 0; for (const x of f) q += (x - m) * (x - m);
      const sd = Math.sqrt(q); return f.map((x) => (sd > 0 ? (x - m) / sd : 0));
    });
    const d = Array.from({ length: n }, () => new Float64Array(n));
    for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) { let r = 0; for (let k = 0; k < z[i].length; k++) r += z[i][k] * z[j][k]; d[i][j] = d[j][i] = 1 - r; }
    const live = new Set(vecs.map((_, i) => i)), leaves = vecs.map((_, i) => [i]), size = new Float64Array(n).fill(1);
    while (live.size > 1) {
      let bi = -1, bj = -1, bd = Infinity;
      for (const i of live) for (const j of live) if (j > i && d[i][j] < bd) { bd = d[i][j]; bi = i; bj = j; }
      for (const k of live) if (k !== bi && k !== bj) d[bi][k] = d[k][bi] = (size[bi] * d[bi][k] + size[bj] * d[bj][k]) / (size[bi] + size[bj]);
      leaves[bi] = leaves[bi].concat(leaves[bj]); size[bi] += size[bj]; live.delete(bj);
    }
    return leaves[[...live][0]];
  }
  function orderGroups() {
    const cols = R.cols.slice(), G = R.genes.length;
    if (V.orderK === "size") cols.sort((a, b) => R.n[b] - R.n[a] || a - b);
    else if (V.orderK === "cluster" && cols.length > 2 && G > 1) {
      const sc = scaledRows(R.cols, "z");
      return clusterOrder(R.cols.map((_, c) => R.genes.map((__, j) => sc[j][c]))).map((i) => R.cols[i]);
    }
    return cols;
  }
  function orderGenes(cols) {
    const idx = R.genes.map((_, j) => j);
    if (V.orderG === "az") idx.sort((a, b) => R.genes[a].localeCompare(R.genes[b], undefined, { numeric: true, sensitivity: "base" }));
    else if (V.orderG === "cluster" && idx.length > 2 && cols.length > 1) return clusterOrder(idx.map((j) => cols.map((k) => R.mean[j * R.K + k])));
    return idx;
  }
  // per gene, across the groups given: "z" = (mean - mean of the group means) / their sd (n - 1);
  // "mm" = (mean - min) / (max - min); a gene that does not vary gives 0 everywhere
  function scaledRows(cols, how) {
    return R.genes.map((_, j) => {
      const v = cols.map((k) => R.mean[j * R.K + k]);
      if (how === "mean") return v;
      if (how === "mm") { const lo = Math.min(...v), hi = Math.max(...v); return v.map((x) => (hi > lo ? (x - lo) / (hi - lo) : 0)); }
      const m = v.reduce((s, x) => s + x, 0) / v.length, sd = v.length > 1 ? Math.sqrt(v.reduce((s, x) => s + (x - m) * (x - m), 0) / (v.length - 1)) : 0;
      return v.map((x) => (sd > 0 ? (x - m) / sd : 0));
    });
  }
  const hex3 = (h) => [1, 3, 5].map((i) => parseInt(h.slice(i, i + 2), 16));
  function ramp(stops, t) {
    t = clamp(Number.isFinite(t) ? t : 0, 0, 1);
    const x = t * (stops.length - 1), i = Math.min(stops.length - 2, Math.floor(x)), f = x - i, a = hex3(stops[i]), b = hex3(stops[i + 1]);
    return `rgb(${[0, 1, 2].map((c) => Math.round(a[c] + (b[c] - a[c]) * f)).join(",")})`;
  }
  function niceStep(raw) { const p = Math.pow(10, Math.floor(Math.log10(raw || 1))), f = raw / p; return (f < 1.5 ? 1 : f < 3 ? 2 : f < 7 ? 5 : 10) * p; }
  let mctx = null;
  function txtW(s, px, weight) {
    if (!mctx) mctx = document.createElement("canvas").getContext("2d");
    mctx.font = `${weight || 400} ${px}px ${FONT}`;
    return mctx.measureText(String(s)).width;
  }
  const short = (s, n) => (String(s).length > n ? String(s).slice(0, n - 1) + "…" : String(s));
  function wrap(text, maxChars) {
    const out = []; let line = "";
    for (const w of text.split(" ")) { if (line && (line + " " + w).length > maxChars) { out.push(line); line = w; } else line = line ? line + " " + w : w; }
    if (line) out.push(line);
    return out;
  }

  // ---- drawing --------------------------------------------------------------------------------------
  function draw() {
    if (!R || $("plotsPanel").hidden) return;
    const panel = $("plotsPanel");
    panel.classList.toggle("pl-is-violin", V.view === "violin");
    panel.classList.toggle("pl-light", V.light);
    $("plTitle").textContent = `${V.view === "violin" ? "Violin" : V.view === "dot" ? "Dot plot" : "Heatmap"} · ${R.grouping.label.replace("Circle → UMAP", "Draw a region")}`;   // (the label is also in the saved files: unchanged there)
    const svg = V.view === "violin" ? drawViolin() : drawMatrix(V.view === "dot");
    lastSvg = svg;
    $("plSvgBox").innerHTML = svg;
    syncViolinRow(); syncSel(); markSel();
    const W = R.words;
    tip(null); syncWarn();
    $("plFoot").textContent =
      (R.norm === "cp10k" ? `Value of a ${W.bin} = log1p(count × 10,000 / its total over all genes); ${W.bins} under ${D.minUmi} UMI left out${R.nLow ? ` (${nf(R.nLow)})` : ""}. `
        : R.norm === "q3" ? `Value of a ${W.bin} = log1p(count / its Q3 normalisation factor), as GeomxTools' Q3 normalisation: factor = the ${W.bin}'s 75th-percentile count ` +
          `(over the targets above LOQ in ≥ 10% of the ${W.bins}, plus the negative probes) / the geometric mean of all ${W.bins}' Q3${R.nLow ? `; ${nf(R.nLow)} ${W.bins} left out (no factor, or under ${D.minUmi} counts)` : ""}. `
        : R.norm === "raw" ? `Value of a ${W.bin} = log1p(count). ` : `Value of a cell = ${R.unit}, as stored in the dataset. `) +
      (R.nQc ? `QC filter on (${R.qcTxt}): ${nf(R.nQc)} ${W.bins} left out. ` : R.mode === "spatial" ? "No QC filter. " : "") +
      (V.view === "violin"
        ? `Violin: Gaussian kernel density, Scott's bandwidth (h = sd × n^−1/5, sd with n − 1), drawn only between the smallest and largest value; every violin has the same greatest width. ` +
          `Box = quartiles, line = median (linear interpolation). Groups with fewer than ${D.kdeMin} ${W.bins} (or all values equal) show points only. Points: a seeded sample of at most ${R.input.base.maxPoints} per group. % = ${W.bins} with ${W.above}. `
        : `Mean = over all ${W.bins} of the group, zeros included. % expressing = ${W.bins} with ${W.above}. n = ${W.bins} of the group after the filters. ` +
          (V.scale === "z" ? `z-score: per gene, across the ${R.cols.length} groups shown (sd with n − 1), colour clipped at ±${D.zClip}. ` : V.scale === "mm" ? "0–1: per gene, (mean − lowest group) / (highest − lowest group). " : "")) +
      `${W.indep}: read these as a description of this ${R.mode === "rna" ? "dataset" : "section"}. Click a cell, a gene or a group for its numbers and for Show on slide.`;
  }

  function syncWarn() {
    if (!R) return;
    const W = R.words, warns = staleWarn(), few = R.cols.filter((k) => R.n[k] < D.fewBins).map((k) => `${R.names[k]} (${R.n[k]})`);
    if (few.length) warns.push(`Fewer than ${D.fewBins} ${W.bins}: ${few.join(", ")}. Their means and violins are rough.`);
    if (R.norm === "raw") warns.push(`No per-${W.bin} totals in this export: values are log1p of the raw counts, not depth-normalised, so a group of deeper ${W.bins} reads higher for every gene.`);
    if (R.cols.length > 30) warns.push(`${R.cols.length} groups: scroll sideways to see them all.`);
    $("plWarn").innerHTML = warns.map(esc).join("<br>"); $("plWarn").hidden = !warns.length;
  }
  function theme() { return V.light ? THEMES.light : THEMES.dark; }
  const svgOpen = (w, h, T) => `<svg xmlns="http://www.w3.org/2000/svg" width="${w}" height="${h}" viewBox="0 0 ${w} ${h}" font-family="${FONT}"><rect width="${w}" height="${h}" fill="${T.bg}"/>`;
  const text = (x, y, s, o = {}) => `<text x="${(+x).toFixed(1)}" y="${(+y).toFixed(1)}" font-size="${o.size || 11}" fill="${o.fill}"${o.anchor ? ` text-anchor="${o.anchor}"` : ""}${o.weight ? ` font-weight="${o.weight}"` : ""}${o.rot ? ` transform="rotate(${o.rot} ${(+x).toFixed(1)} ${(+y).toFixed(1)})"` : ""}${o.attr || ""}>${esc(s)}</text>`;
  const nLabel = (k) => `n = ${nf(R.n[k])}`;

  // rows = genes, columns = groups; a filled cell (heatmap) or a dot (area = % expressing)
  function drawMatrix(dot) {
    const T = theme(), host = $("plHost"), Wh = Math.max(260, host.clientWidth - 2);
    const cols = orderGroups(), rows = orderGenes(cols), K = cols.length, G = rows.length;
    const sc = scaledRows(cols, V.scale);
    // the colour scale's ends
    let lo = 0, hi = 1, stops = T.seq, barTitle = `Mean expression · ${R.unit}`;
    if (V.scale === "mean") { hi = 0; for (const r of sc) for (const x of r) if (x > hi) hi = x; if (!(hi > 0)) hi = 1; }
    else if (V.scale === "z") { let m = 0; for (const r of sc) for (const x of r) m = Math.max(m, Math.abs(x)); hi = clamp(m || 1, 1, D.zClip); lo = -hi; stops = T.div; barTitle = `z-score per gene, across the ${K} groups shown${m > D.zClip ? ` (clipped at ±${D.zClip})` : ""}`; }
    else barTitle = "0–1 per gene (lowest to highest group mean)";
    const col = (x) => ramp(stops, (x - lo) / (hi - lo));
    let fMax = 0; for (const j of rows) for (const k of cols) fMax = Math.max(fMax, R.frac[j * R.K + k] || 0);
    const dotMax = clamp(Math.ceil(fMax * 10 - 1e-9) / 10, 0.1, 1);
    // geometry
    const lab = cols.map((k) => short(R.names[k], 28));
    const geneLab = rows.map((j) => short(R.genes[j] + (R.noCounts.includes(R.genes[j]) ? " (no counts)" : ""), 30));
    const left = Math.max(...geneLab.map((s) => txtW(s, 11))) + 16;
    const lab2 = Math.max(...cols.map((k, i) => Math.max(txtW(lab[i], 11), txtW(nLabel(k), 10) + 12)));
    const cw = clamp(Math.floor((Wh - left - Math.min(Math.cos(50 * Math.PI / 180) * lab2 + 8, 90)) / K), dot ? 18 : 16, 56), ch = dot ? clamp(cw, 18, 26) : G <= 12 ? 22 : 16;
    // two lines (name, n) need room between neighbours; narrow columns get "name · n = ..." on one steeper line
    const oneLine = cw * Math.sin(50 * Math.PI / 180) < 24, deg = oneLine ? 62 : 50, ang = deg * Math.PI / 180;
    const labW = oneLine ? Math.max(...cols.map((k, i) => txtW(`${lab[i]} · ${nLabel(k)}`, 11))) : lab2;
    const top = Math.sin(ang) * labW + (oneLine ? 18 : 30), over = Math.cos(ang) * labW + 8;
    const right = Math.max(14, over - cw / 2);
    const legendW = 190, legendDotW = dot ? 200 : 0, stack = left + legendW + legendDotW + 20 > left + K * cw + right && dot;
    const legH = (stack ? 112 : 62);
    const W = Math.max(left + K * cw + right, left + legendW + (stack ? 0 : legendDotW) + 16, 300), H = top + G * ch + legH + 34;
    const p = [svgOpen(W, H, T)], gid = `plg${++uid}`;
    p.push(`<defs><linearGradient id="${gid}" x1="0" x2="1" y1="0" y2="0">${Array.from({ length: 11 }, (_, i) => `<stop offset="${i * 10}%" stop-color="${ramp(stops, i / 10)}"/>`).join("")}</linearGradient></defs>`);
    // column labels: name and n, slanted
    cols.forEach((k, c) => {
      const x = left + c * cw + cw / 2, y = top - 8;
      const nCol = R.n[k] < D.fewBins ? T.warn : T.muted;
      p.push(`<g transform="translate(${x.toFixed(1)} ${y}) rotate(-${deg})" data-k="${k}" style="cursor:pointer">` + (oneLine
        ? `<text x="0" y="3" font-size="11" fill="${T.ink}">${esc(lab[c])}<tspan fill="${nCol}" font-size="10"> · ${esc(nLabel(k))}</tspan></text></g>`
        : `<text x="0" y="0" font-size="11" fill="${T.ink}">${esc(lab[c])}</text><text x="12" y="12" font-size="10" fill="${nCol}">${esc(nLabel(k))}</text></g>`));
    });
    // rows
    rows.forEach((j, r) => {
      const y = top + r * ch, none = R.noCounts.includes(R.genes[j]);
      p.push(text(left - 8, y + ch / 2 + 4, geneLab[r], { fill: none ? T.muted : T.ink, anchor: "end", attr: ` data-g="${j}" style="cursor:pointer"` }));
      cols.forEach((k, c) => {
        const x = left + c * cw, o = j * R.K + k, v = sc[j][c];
        if (dot) {
          const f = R.frac[o] || 0, rad = Math.sqrt(Math.min(1, f / dotMax)) * (Math.min(cw, ch) / 2 - 1.5);
          p.push(`<rect x="${x}" y="${y}" width="${cw}" height="${ch}" fill="${T.bg}" data-g="${j}" data-k="${k}"/>`);
          if (f > 0) p.push(`<circle cx="${x + cw / 2}" cy="${y + ch / 2}" r="${Math.max(0.8, rad).toFixed(2)}" fill="${col(v)}" stroke="${T.line}" stroke-width="0.6" pointer-events="none"/>`);
        } else p.push(`<rect x="${x + 0.5}" y="${y + 0.5}" width="${cw - 1}" height="${ch - 1}" fill="${none ? T.bg : col(v)}"${none ? ` stroke="${T.grid}"` : ""} data-g="${j}" data-k="${k}"/>`);
      });
    });
    if (dot) p.push(`<rect x="${left}" y="${top}" width="${K * cw}" height="${G * ch}" fill="none" stroke="${T.grid}"/>`);
    p.push(`<g id="plMark"></g>`);
    // legends
    const ly = top + G * ch + 22, bw = 150;
    p.push(text(left, ly, barTitle, { fill: T.ink, size: 10.5 }));
    p.push(`<rect x="${left}" y="${ly + 7}" width="${bw}" height="10" fill="url(#${gid})" stroke="${T.line}" stroke-width="0.5"/>`);
    const fmt = (x) => (Math.abs(x) >= 10 ? x.toFixed(0) : x.toFixed(Math.abs(x) < 1 && x !== 0 ? 2 : 1));
    for (const [f, v] of [[0, lo], [0.5, (lo + hi) / 2], [1, hi]]) p.push(text(left + f * bw, ly + 30, V.scale === "mm" ? (f === 0 ? "0 (lowest)" : f === 1 ? "1 (highest)" : "0.5") : fmt(v), { fill: T.muted, size: 10, anchor: f === 0 ? "start" : f === 1 ? "end" : "middle" }));
    if (dot) {
      const dx = stack ? left : left + legendW + 24, dy = stack ? ly + 48 : ly, rm = Math.min(cw, ch) / 2 - 1.5;
      p.push(text(dx, dy, `% of ${R.words.bins} expressing (dot area)`, { fill: T.ink, size: 10.5 }));
      const steps = [0.25, 0.5, 0.75, 1].map((s) => s * dotMax);
      steps.forEach((f, i) => {
        const cx = dx + 14 + i * 46, rad = Math.sqrt(f / dotMax) * rm;
        p.push(`<circle cx="${cx}" cy="${dy + 16}" r="${rad.toFixed(2)}" fill="${T.muted}"/>`);
        p.push(text(cx, dy + 16 + rm + 12, `${Math.round(f * 100)}%`, { fill: T.muted, size: 10, anchor: "middle" }));
      });
    }
    const cap = `${R.dataset} · ${R.layer} · ${R.grouping.label} · n = ${R.words.bins} per group` + (R.nQc ? " · QC filter on" : "");
    p.push(text(8, H - 8, short(cap, Math.floor((W - 16) / 5.2)), { fill: T.muted, size: 10 }));
    p.push("</svg>");
    return p.join("");
  }

  // one panel per gene; a violin per group
  function drawViolin() {
    const T = theme(), host = $("plHost"), Wh = Math.max(260, host.clientWidth - 2), Hh = Math.max(220, host.clientHeight - 4);
    const cols = orderGroups(), K = cols.length, genes = V.vGenes.filter((j) => R.detail.has(j)), pop = V.exprOnly ? "expr" : "all", W_ = R.words;
    const lab = cols.map((k) => short(R.names[k], 26)), ang = 50 * Math.PI / 180;
    const labW = Math.max(...cols.map((k, i) => Math.max(txtW(lab[i], 11), txtW(nLabel(k), 10))));
    const left = 48 + Math.max(0, Math.cos(ang) * labW - 60), right = 14, bottom = Math.sin(ang) * labW + 34;
    const bw = clamp(Math.floor((Wh - left - right) / K), 34, 96), W = Math.max(left + K * bw + right, 330);
    const cap = wrap(`% above each violin = ${W_.bins} with ${W_.above} (of n). Box: quartiles; line: median. Density: Gaussian kernel, Scott's bandwidth, cut at the data's range, equal greatest width.` + (R.nQc ? " QC filter on." : ""), Math.max(40, Math.floor((W - 16) / 5.2)));
    const capH = cap.length * 12 + 10;
    const headH = 22, annH = V.exprOnly ? 26 : 15;
    const ph = clamp(Math.floor((Hh - bottom - capH - genes.length * (headH + annH + 10)) / Math.max(1, genes.length)), 120, 340);
    const H = genes.length * (headH + annH + ph + 10) + bottom + capH;
    const p = [svgOpen(W, Math.max(H, 120), T)];
    if (!genes.length) { p.push(text(left, 40, "Pick a gene above (or click a gene in the heatmap).", { fill: T.muted, size: 12 })); p.push("</svg>"); return p.join(""); }
    let y0 = 0;
    genes.forEach((j, gi) => {
      const det = R.detail.get(j), name = R.genes[j];
      let top = 0; for (const k of cols) { const d = det[k][pop]; if (d.n && d.max > top) top = d.max; }
      if (!(top > 0)) top = 1;
      const step = niceStep(top / 5), yMax = Math.ceil(top / step - 1e-9) * step;
      const py = y0 + headH + annH, Y = (v) => py + ph - (clamp(v, 0, yMax) / yMax) * ph;
      p.push(text(left, y0 + 15, name, { fill: T.ink, size: 12.5, weight: 600, attr: ` data-g="${j}" style="cursor:pointer"` }));
      p.push(text(left + txtW(name, 12.5, 600) + 10, y0 + 15, V.exprOnly ? `expressing ${W_.bins} only (${W_.above})` : `all ${W_.bins}, zeros included`, { fill: T.muted, size: 10.5 }));
      for (let v = 0; v <= yMax + 1e-9; v += step) {
        p.push(`<line x1="${left}" x2="${left + K * bw}" y1="${Y(v).toFixed(1)}" y2="${Y(v).toFixed(1)}" stroke="${T.grid}"/>`);
        p.push(text(left - 6, Y(v) + 3.5, +v.toFixed(6), { fill: T.muted, size: 10, anchor: "end" }));
      }
      p.push(`<text transform="translate(12 ${(py + ph / 2).toFixed(1)}) rotate(-90)" text-anchor="middle" font-size="10.5" fill="${T.muted}">${esc(R.unit)}</text>`);
      cols.forEach((k, c) => {
        const cx = left + c * bw + bw / 2, d = det[k][pop], f = R.frac[j * R.K + k], colr = R.colors[k] || T.muted, half = bw * 0.42;
        p.push(`<rect x="${left + c * bw}" y="${py}" width="${bw}" height="${ph}" fill="${T.bg}" fill-opacity="0" data-g="${j}" data-k="${k}"/>`);
        p.push(text(cx, y0 + headH + 10, `${(f * 100).toFixed(f < 0.1 ? 1 : 0)}%`, { fill: T.muted, size: 9.5, anchor: "middle" }));
        if (V.exprOnly) p.push(text(cx, y0 + headH + 21, `n=${nf(d.n)}`, { fill: d.n < D.fewBins ? T.warn : T.muted, size: 9, anchor: "middle" }));
        if (!d.n) return;
        let widthAt = () => half * 0.6;
        if (d.dens) {
          let dm = 0; for (const x of d.dens) if (x > dm) dm = x;
          const Gp = d.dens.length, pts = [];
          for (let i = 0; i < Gp; i++) pts.push([(d.dens[i] / dm) * half, Y(d.lo + ((d.hi - d.lo) * i) / (Gp - 1))]);
          const path = pts.map((q, i) => `${i ? "L" : "M"}${(cx + q[0]).toFixed(1)} ${q[1].toFixed(1)}`).join("") + pts.reverse().map((q) => `L${(cx - q[0]).toFixed(1)} ${q[1].toFixed(1)}`).join("") + "Z";
          p.push(`<path d="${path}" fill="${colr}" fill-opacity="0.5" stroke="${colr}" stroke-width="1.2" pointer-events="none"/>`);
          widthAt = (v) => { const t = clamp(((v - d.lo) / (d.hi - d.lo)) * (Gp - 1), 0, Gp - 1), i = Math.min(Gp - 2, Math.floor(t)); return Math.max(1.5, ((d.dens[i] + (d.dens[i + 1] - d.dens[i]) * (t - i)) / dm) * half * 0.9); };
        }
        if (V.points && d.pts) {
          const q = [];
          for (let i = 0; i < d.pts.length; i++) {
            const u = Math.sin((i + 1) * 12.9898 + (k + 1) * 78.233 + (j + 1) * 37.719) * 43758.5453, rnd = u - Math.floor(u);
            q.push(`<circle cx="${(cx + (rnd * 2 - 1) * widthAt(d.pts[i])).toFixed(1)}" cy="${Y(d.pts[i]).toFixed(1)}" r="1.3"/>`);
          }
          p.push(`<g fill="${T.ink}" fill-opacity="0.45" pointer-events="none">${q.join("")}</g>`);
        }
        // quartile box and median
        p.push(`<rect x="${(cx - 2.5).toFixed(1)}" y="${Y(d.q3).toFixed(1)}" width="5" height="${Math.max(0.8, Y(d.q1) - Y(d.q3)).toFixed(1)}" fill="${T.ink}" fill-opacity="0.85" pointer-events="none"/>`);
        p.push(`<line x1="${(cx - 7).toFixed(1)}" x2="${(cx + 7).toFixed(1)}" y1="${Y(d.median).toFixed(1)}" y2="${Y(d.median).toFixed(1)}" stroke="${T.hi}" stroke-width="2" pointer-events="none"/>`);
        p.push(`<line x1="${(cx - 7).toFixed(1)}" x2="${(cx + 7).toFixed(1)}" y1="${Y(d.median).toFixed(1)}" y2="${Y(d.median).toFixed(1)}" stroke="${T.bg}" stroke-width="0.6" pointer-events="none"/>`);
      });
      p.push(`<line x1="${left}" x2="${left + K * bw}" y1="${py + ph}" y2="${py + ph}" stroke="${T.line}"/>`);
      y0 = py + ph + 10;
      if (gi === genes.length - 1) cols.forEach((k, c) => {
        const x = left + c * bw + bw / 2 + 4, y = y0 + 6;
        p.push(`<g transform="translate(${x.toFixed(1)} ${y}) rotate(-50)" data-k="${k}" style="cursor:pointer">` +
          `<text x="0" y="0" text-anchor="end" font-size="11" fill="${T.ink}">${esc(lab[c])}</text>` +
          `<text x="-8" y="12" text-anchor="end" font-size="10" fill="${R.n[k] < D.fewBins ? T.warn : T.muted}">${esc(nLabel(k))}</text></g>`);
      });
    });
    p.push(`<g id="plMark"></g>`);
    cap.forEach((l, i) => p.push(text(8, Math.max(H, 120) - capH + 14 + i * 12, l, { fill: T.muted, size: 10 })));
    p.push("</svg>");
    return p.join("");
  }

  function syncViolinRow() {
    if (!R) return;
    $("plVChips").innerHTML = V.vGenes.map((j) => `<span class="pl-chip">${esc(R.genes[j])}<button data-rm="${j}" title="Remove">×</button></span>`).join("");
    const opts = R.genes.map((n, j) => (V.vGenes.includes(j) ? "" : `<option value="${j}">${esc(n)}</option>`)).join("");
    $("plVAdd").innerHTML = `<option value="">${V.vGenes.length >= D.maxViolin ? "Replace the first with…" : "Add a gene…"}</option>${opts}`;
  }

  // ---- hover, selection, Show on slide ------------------------------------------------------------
  function hit(ev) {
    const el = ev.target && ev.target.closest ? ev.target.closest("[data-g],[data-k]") : null;
    if (!el || !R) return null;
    const g = el.getAttribute("data-g"), k = el.getAttribute("data-k");
    return { g: g === null ? null : Number(g), k: k === null ? null : Number(k) };
  }
  function cellTxt(g, k) {
    const o = g * R.K + k, W = R.words, lines = [`<b>${esc(R.genes[g])}</b> in <b>${esc(R.names[k])}</b>`];
    lines.push(`mean ${R.mean[o].toFixed(3)} <small>${esc(R.unit)}</small>`);
    if (R.norm !== "asis" || R.input.base.asisLog) lines.push(`<small>mean on the linear scale ${R.meanLin[o].toFixed(2)} ${R.norm === "raw" ? "UMI" : R.norm === "q3" ? "Q3-normalised counts" : "per 10k"}</small>`);
    lines.push(`${(R.frac[o] * 100).toFixed(1)}% of ${W.bins} expressing (${nf(R.nExpr[o])})`);
    lines.push(`n = ${nf(R.n[k])} ${W.bins}`);
    const det = R.detail.get(g);
    if (V.view === "violin" && det) {
      const d = det[k][V.exprOnly ? "expr" : "all"];
      if (d.n) lines.push(`median ${d.median.toFixed(3)} · quartiles ${d.q1.toFixed(3)} – ${d.q3.toFixed(3)}` + (V.exprOnly ? ` <small>(of the ${nf(d.n)} expressing)</small>` : ""),
        d.dens ? `<small>bandwidth ${d.h.toPrecision(3)}</small>` : `<small>too few ${W.bins} (or one value) for a density: points only</small>`);
    } else if (V.scale !== "mean" && V.view !== "violin") {
      const cols = orderGroups(), c = cols.indexOf(k);
      if (c >= 0) lines.push(`<small>${V.scale === "z" ? "z-score" : "0–1"} across groups: ${scaledRows(cols, V.scale)[g][c].toFixed(2)}</small>`);
    }
    return lines.join("<br>");
  }
  function tip(ev) {
    const el = $("plTip"), t = ev ? hit(ev) : null;
    if (!t || !R) { el.hidden = true; return; }
    el.innerHTML = t.g !== null && t.k !== null ? cellTxt(t.g, t.k)
      : t.g !== null ? `<b>${esc(R.genes[t.g])}</b>${R.from[R.genes[t.g]] && R.from[R.genes[t.g]] !== "typed" ? `<br><small>marker of ${esc(R.from[R.genes[t.g]])}</small>` : ""}<br><small>click: violin, Show on slide</small>`
      : `<b>${esc(R.names[t.k])}</b><br>n = ${nf(R.n[t.k])} ${R.words.bins}` + (R.qcOut[t.k] ? `<br><small>${nf(R.qcOut[t.k])} left out by the QC filter</small>` : "") + (R.low[t.k] ? `<br><small>${nf(R.low[t.k])} under ${D.minUmi} UMI left out</small>` : "") + `<br><small>click: Show on slide</small>`;
    el.hidden = false;
    const r = $("plHost").getBoundingClientRect(), x = ev.clientX - r.left, y = ev.clientY - r.top;
    el.style.left = `${Math.max(2, Math.min(x + 12, r.width - el.offsetWidth - 4))}px`;
    el.style.top = `${y + 14 + el.offsetHeight > r.height ? Math.max(2, y - el.offsetHeight - 8) : y + 14}px`;
  }
  function syncSel() {
    const el = $("plSel"), s = V.sel;
    if (!R || !s) { el.innerHTML = `<small>Click a cell, a gene name or a group name.</small>`; return; }
    const parts = [];
    if (s.g !== null && s.k !== null) parts.push(cellTxt(s.g, s.k).split("<br>").slice(0, 4).join(" · "));
    else if (s.g !== null) parts.push(`<b>${esc(R.genes[s.g])}</b>`);
    else parts.push(`<b>${esc(R.names[s.k])}</b> · n = ${nf(R.n[s.k])} ${R.words.bins}`);
    const links = [];
    if (s.g !== null) links.push(`<a href="#" data-act="gene">Show ${esc(R.genes[s.g])} on ${R.mode === "rna" ? "the UMAP" : "the slide"}</a>`);
    if (s.g !== null && !(V.view === "violin" && V.vGenes.includes(s.g))) links.push(`<a href="#" data-act="violin">Violin of ${esc(R.genes[s.g])}</a>`);
    if (s.k !== null) links.push(`<a href="#" data-act="group">${R.mode === "rna" ? `Colour the UMAP by ${esc(R.grouping.label)}` : `Show ${esc(short(R.names[s.k], 24))} on the slide`}</a>`);
    el.innerHTML = `${parts.join("")}<span class="pl-links">${links.join("")}</span>`;
  }
  function markSel() {
    const g = document.querySelector("#plSvgBox #plMark"), s = V.sel;
    if (!g) return;
    g.innerHTML = "";
    if (!s || s.g === null || s.k === null) return;
    const el = document.querySelector(`#plSvgBox rect[data-g="${s.g}"][data-k="${s.k}"]`);
    if (!el) return;
    const b = el.getBBox();
    g.innerHTML = `<rect x="${b.x}" y="${b.y}" width="${b.width}" height="${b.height}" fill="none" stroke="${theme().hi}" stroke-width="1.5" pointer-events="none"/>`;
  }
  async function selAction(act) {
    const s = V.sel; if (!R || !s) return;
    if (R.sig.split("|").slice(0, 2).join("|") !== cacheKey()) { msg("These plots belong to another layer or dataset: draw them again.", true); return; }
    if (act === "violin" && s.g !== null) { addViolinGene(s.g); return; }
    if (act === "gene" && s.g !== null) { if (R.mode === "rna") api.addRnaGene(R.genes[s.g]); else api.addGene(R.genes[s.g]); }
    if (act === "group" && s.k !== null) {
      if (R.mode === "rna") { api.rnaColorBy(R.grouping.label); }
      else {
        const N = R.eff.length, mask = new Uint8Array(N); let n = 0;
        for (let r = 0; r < N; r++) if (R.eff[r] === s.k) { mask[r] = 1; n++; }
        const shown = await api.lightRows(mask, n);
        msg(shown ? `${R.names[s.k]}: ${nf(n)} ${R.words.bins} lit on the slide (Esc or "Clear selection" puts it back).` : `Could not show ${R.names[s.k]} on the slide.`);
      }
    }
    if (api.touch) { $("plotsPanel").hidden = true; syncButtons(); }     // a phone's panel covers the slide
  }

  // ---- saving ---------------------------------------------------------------------------------------
  const slug = (s) => String(s).replace(/[^A-Za-z0-9]+/g, "_").replace(/^_|_$/g, "").slice(0, 40);
  const fileBase = () => `${R.dataset}_${R.layer}_${V.view === "violin" ? (V.exprOnly ? "violin_expressing" : "violin") : V.view === "dot" ? "dotplot" : "heatmap"}_${slug(R.grouping.label)}`;
  // without the interaction-only attributes, so the file is a plain picture
  const cleanSvg = () => lastSvg.replace(/ data-[gk]="\d+"/g, "").replace(/ style="cursor:pointer"/g, "").replace(/<g id="plMark"><\/g>/, "");
  function saveSvg() { if (R && lastSvg) api.download(`${fileBase()}.svg`, `<?xml version="1.0" encoding="UTF-8"?>\n${cleanSvg()}\n`, "image/svg+xml"); }
  function savePng() {
    if (!R || !lastSvg) return;
    const svg = cleanSvg(), m = /width="([\d.]+)" height="([\d.]+)"/.exec(svg), w = Number(m[1]), h = Number(m[2]);
    const scale = Math.min(2, 8000 / Math.max(w, h)), img = new Image(), name = `${fileBase()}.png`;
    img.onload = () => {
      const c = document.createElement("canvas"); c.width = Math.round(w * scale); c.height = Math.round(h * scale);
      c.getContext("2d").drawImage(img, 0, 0, c.width, c.height);
      c.toBlob((b) => { if (b) api.download(name, b, "image/png"); else msg("The PNG could not be made: save the SVG instead.", true); }, "image/png");
    };
    img.onerror = () => msg("The PNG could not be made: save the SVG instead.", true);
    img.src = "data:image/svg+xml;charset=utf-8," + encodeURIComponent(svg);
  }
  // the numbers plotted, one row per gene and group, in the order drawn
  function saveCsv() {
    if (!R) return;
    const q = (t) => `"${String(t).replace(/"/g, '""')}"`, f = (x, d) => (Number.isFinite(x) ? (+x).toFixed(d) : "");
    const cols = orderGroups(), rows = V.view === "violin" ? V.vGenes.filter((j) => R.detail.has(j)) : orderGenes(cols), lines = [];
    const common = (j, k) => [q(R.dataset), R.layer, q(R.grouping.label), q(R.names[k]), R.n[k], R.qcOut[k], R.low[k], q(R.genes[j]), q(R.unit)];
    const head = "dataset,layer,grouping,group,n_bins,bins_left_out_qc,bins_left_out_low_umi,gene,value_unit";
    if (V.view === "violin") {
      lines.push(`${head},population,n_plotted,n_expressing,frac_expressing,mean,q1,median,q3,min,max,sd,kde_bandwidth_scott,kde_from,kde_to,kde_density`);
      for (const j of rows) for (const k of cols) for (const pop of V.exprOnly ? ["expr"] : ["all"]) {
        const d = R.detail.get(j)[k][pop], o = j * R.K + k;
        lines.push([...common(j, k), pop === "expr" ? "expressing bins only" : "all bins", d.n, R.nExpr[o], f(R.frac[o], 6), f(d.mean, 6), f(d.q1, 6), f(d.median, 6), f(d.q3, 6), f(d.min, 6), f(d.max, 6),
          f(d.sd, 6), d.dens ? f(d.h, 6) : "", d.dens ? f(d.lo, 6) : "", d.dens ? f(d.hi, 6) : "", d.dens ? q(Array.from(d.dens, (x) => x.toPrecision(5)).join(" ")) : ""].join(","));
      }
    } else {
      const sc = V.scale === "mean" ? null : scaledRows(cols, V.scale);
      lines.push(`${head},mean,frac_expressing,n_expressing,mean_linear,colour_scale,colour_value,marker_of`);
      for (const j of rows) cols.forEach((k, c) => {
        const o = j * R.K + k;
        lines.push([...common(j, k), f(R.mean[o], 6), f(R.frac[o], 6), R.nExpr[o], f(R.meanLin[o], 4), V.scale === "mean" ? "mean" : V.scale === "z" ? "z-score per gene" : "0-1 per gene",
          f(sc ? sc[j][c] : R.mean[o], 6), q(R.from[R.genes[j]] === "typed" ? "" : R.from[R.genes[j]] || "")].join(","));
      });
    }
    const note = `# ${R.norm === "cp10k" ? `value = log1p(count x 10000 / total); bins under ${D.minUmi} UMI left out`
      : R.norm === "q3" ? "value = log1p(count / Q3 normFactor) (GeomxTools Q3 normalisation, factors from the export's q3 entry); rows without a factor left out" : R.norm === "raw" ? "value = log1p(count); NOT depth-normalised (no per-bin totals in the export)" : `value = ${R.unit}`}` +
      `; QC filter: ${R.qcTxt || "none"}; mean over all bins of the group (zeros included); frac_expressing = bins with count > 0 / n_bins`;
    api.download(`${fileBase()}.csv`, note + "\n" + lines.join("\n") + "\n", "text/csv");
  }

  window.SodaPlots = {
    init(a) { api = a; buildDom(); refresh(); },
    refresh,
    // a new layer or dataset: the old result belongs to other bins
    reset() { if (!api) return; cancel(); R = null; lastGroups = ""; lastSvg = ""; V.sel = null; V.vGenes = []; window.SodaPlots.result = null; cache = blankCache("");
              if ($("plotsPanel")) { $("plotsPanel").hidden = true; msg(""); } refresh(); },
    result: null, force: "",
    // for tests from the console
    _view: () => ({ ...V }), _svg: () => lastSvg, _draw: runPlot, _setView: setView, _cache: () => cache,
  };
})();
