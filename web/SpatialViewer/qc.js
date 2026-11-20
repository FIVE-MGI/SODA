// Per-bin QC filter: total counts, genes detected and % mitochondrial for every bin of a
// sequencing layer (export.py -> binqc.py: expr/<key>/binqc/*.bin, named by the manifest
// layer's `binQc` block). For each metric a histogram with a min and a max slider; a bin
// outside any threshold is "removed": greyed on the slide and left out of Circle -> UMAP
// (the worker gets a keep mask, one byte per row of positions.bin).
//
// OFF BY DEFAULT. Every bin is kept until a slider is moved or a number typed, so no earlier
// result changes by opening a newer export. A slider at either end means "no bound" (not
// "the data's min / max"), and a filter that removes no bin sends no mask at all, so the
// worker runs exactly the unfiltered path.
//
// Older exports have no `binQc` block (or a source without .buf): the
// panel stays hidden and keepMask() is null.
//
// The thresholds travel with Circle -> UMAP's Save JSON (`qcFilter`) and in the page URL
// (?qc=<layer>:total=20-,pctMt=-25), so a result can be repeated.
//
// Only unsigned typed arrays here: the de-identification check reads the signed names as
// sample IDs.
/* exported SodaBinQc */
(function () {
  const STEPS = 1000;                      // slider positions; 0 and STEPS = no bound
  const NBARS = 72;                        // histogram bars
  const GREY = [88, 88, 88, 215];          // removed bins on the slide
  const NONE_TXT = "any";

  function create(ctx) {
    const $ = ctx.$, esc = ctx.esc;
    const say = (t) => (ctx.words ? ctx.words(t) : t);   // "bin" -> "cell" on a cells layer (viewer.js)
    // the page's ?qc= value, applied once to the layer it names (on the dataset in ?ds=)
    const params = new URLSearchParams(location.search);
    let pending = parseParam(params.get("qc"), params.get("ds"));
    const saved = new Map();               // "<ds>|<layer>" -> thresholds, kept while switching layers
    let greyPref = true;                   // "Grey out removed bins", kept across layers
    let Q = blank();

    function blank() {
      return { ds: null, key: null, meta: null, n: 0, metrics: [], th: {}, keep: null, removed: 0,
               perMetric: {}, grey: greyPref, pos: null, layerData: null, ver: 0, raf: 0 };
    }

    // ---- loading ------------------------------------------------------------------
    function reset() {
      Q = blank();
      $("panelBinQc").hidden = true;
    }

    async function onLayer(dsId, key, e, src) {
      if (Q.ds && Q.key) saved.set(`${Q.ds}|${Q.key}`, { ...Q.th });
      Q = blank();
      Q.ds = dsId; Q.key = key;
      $("bqStale").hidden = true;
      const meta = e && e.kind !== "protein" ? e.binQc : null;
      if (!meta || !Array.isArray(meta.metrics) || !meta.metrics.length || !src || !src.buf) {
        $("panelBinQc").hidden = true; syncUrl(); ctx.onChange(); return;
      }
      try {
        const dir = meta.dir || "binqc";
        const bufs = await Promise.all(meta.metrics.map((m) => src.buf(`expr/${key}/${dir}/${m.file}`)));
        if (Q.key !== key || Q.ds !== dsId) return;        // the layer changed meanwhile
        Q.metrics = meta.metrics.map((m, i) => {
          const v = m.dtype === "uint16" ? new Uint16Array(bufs[i]) : m.dtype === "uint32" ? new Uint32Array(bufs[i]) : new Float32Array(bufs[i]);
          return { ...m, v, log: m.scale === "log" };
        });
        Q.n = Q.metrics[0].v.length;
        if (Q.metrics.some((m) => m.v.length !== Q.n)) throw new Error("QC arrays differ in length");
        Q.metrics.forEach((m) => { if (!Number.isFinite(m.min) || !Number.isFinite(m.max)) { const r = range(m.v); m.min = r[0]; m.max = r[1]; } });
        Q.meta = meta;
      } catch (err) {
        Q.meta = null; Q.metrics = [];
        $("panelBinQc").hidden = true; syncUrl(); ctx.onChange(); return;
      }
      let th = saved.get(`${dsId}|${key}`) || {};
      if (pending && pending.ds === dsId && pending.key === key) { th = pending.th; pending = null; }
      Q.th = cleanTh(th);
      build();
      recompute(true);
    }

    function range(v) {
      let lo = Infinity, hi = -Infinity;
      for (let i = 0; i < v.length; i++) { const x = v[i]; if (x < lo) lo = x; if (x > hi) hi = x; }
      return [lo === Infinity ? 0 : lo, hi === -Infinity ? 0 : hi];
    }
    // only known metrics, finite numbers; {min: null, max: null} entries dropped
    function cleanTh(th) {
      const out = {};
      for (const m of Q.metrics) {
        const t = th && th[m.key];
        if (!t) continue;
        const lo = Number.isFinite(t.min) ? t.min : null, hi = Number.isFinite(t.max) ? t.max : null;
        if (lo !== null || hi !== null) out[m.key] = { min: lo, max: hi };
      }
      return out;
    }

    // ---- axis: slider position <-> value ------------------------------------------------
    const tf = (m, x) => (m.log ? Math.log10(1 + Math.max(0, x)) : x);
    const itf = (m, t) => (m.log ? Math.pow(10, t) - 1 : t);
    function posOf(m, x) {
      const a = tf(m, m.min), b = tf(m, m.max);
      return b > a ? Math.round(((tf(m, x) - a) / (b - a)) * STEPS) : 0;
    }
    function valOf(m, p) {
      const a = tf(m, m.min), b = tf(m, m.max), x = itf(m, a + (p / STEPS) * (b - a));
      // counts are whole numbers; % keeps two decimals
      return m.key === "pctMt" || m.dtype === "float32" && !m.log ? Math.round(x * 100) / 100 : Math.round(x);
    }
    const fmt = (m, x) => (x === null ? NONE_TXT : m.key === "pctMt" ? `${(+x).toFixed(x >= 10 ? 1 : 2)}%` : Math.round(x).toLocaleString());

    // ---- the panel ------------------------------------------------------------------
    function build() {
      const box = $("bqMetrics");
      box.innerHTML = "";
      for (const m of Q.metrics) {
        const el = document.createElement("div");
        el.className = "bq-metric";
        el.innerHTML =
          `<div class="bq-head"><span class="bq-name" title="${esc(tipOf(m))}">${esc(m.label)}</span><span class="bq-cut"></span></div>` +
          `<canvas class="bq-hist" title="${say("Bins")} per value (square-root height${m.log ? ", log axis" : ""}); grey = removed"></canvas>` +
          `<div class="bq-axis"><span>${esc(fmt(m, m.min))}</span><span>${esc(fmt(m, m.max))}</span></div>` +
          `<label class="bq-sl"><span>Min</span><input type="range" class="bq-lo" min="0" max="${STEPS}" step="1" value="0">` +
          `<input type="text" inputmode="decimal" class="bq-lon" placeholder="${NONE_TXT}"></label>` +
          `<label class="bq-sl"><span>Max</span><input type="range" class="bq-hi" min="0" max="${STEPS}" step="1" value="${STEPS}">` +
          `<input type="text" inputmode="decimal" class="bq-hin" placeholder="${NONE_TXT}"></label>`;
        m.el = el;
        const lo = el.querySelector(".bq-lo"), hi = el.querySelector(".bq-hi");
        const lon = el.querySelector(".bq-lon"), hin = el.querySelector(".bq-hin");
        const setT = (which, val) => {
          const t = { ...(Q.th[m.key] || { min: null, max: null }) };
          t[which] = val;
          if (t.min !== null && t.max !== null && t.min > t.max) { if (which === "min") t.max = t.min; else t.min = t.max; }
          if (t.min === null && t.max === null) delete Q.th[m.key]; else Q.th[m.key] = t;
          syncMetric(m); schedule();
        };
        lo.oninput = () => { const p = +lo.value; if (p >= +hi.value && +hi.value < STEPS) lo.value = hi.value; setT("min", +lo.value <= 0 ? null : valOf(m, +lo.value)); };
        hi.oninput = () => { const p = +hi.value; if (p <= +lo.value && +lo.value > 0) hi.value = lo.value; setT("max", +hi.value >= STEPS ? null : valOf(m, +hi.value)); };
        const typed = (inp, which) => () => {
          const s = inp.value.trim().replace(/%$/, "").replace(/,/g, "");
          if (s === "" || /^any$/i.test(s)) { setT(which, null); return; }
          const x = Number(s);
          if (!Number.isFinite(x) || x < 0) { syncMetric(m); return; }
          setT(which, x);
        };
        lon.onchange = typed(lon, "min"); hin.onchange = typed(hin, "max");
        [lon, hin].forEach((i) => { i.onkeydown = (ev) => { ev.stopPropagation(); if (ev.key === "Enter") i.blur(); }; });
        box.appendChild(el);
        hist(m);
        syncMetric(m);
      }
      $("bqGrey").checked = Q.grey;
      $("panelBinQc").hidden = false;
      // the canvases get their size only once the panel is shown
      requestAnimationFrame(() => Q.metrics.forEach(drawHist));
    }
    function tipOf(m) {
      const base = m.key === "total" ? "UMI over all genes of the bin (the library size the UMAP normalises by)"
        : m.key === "genes" ? "Genes with at least one UMI in the bin"
        : m.key === "pctMt" ? `UMI of mitochondrial genes (${(Q.meta.mtGenes || []).length} genes named MT-…) as a share of the bin's total; 0 for an empty bin`
        : m.label;
      const b = say(base);
      return `${b.charAt(0).toUpperCase() + b.slice(1)}. Range ${fmt(m, m.min)} – ${fmt(m, m.max)}, median ${fmt(m, m.median)}.`;
    }

    // per metric: the bar counts on its own axis (computed once per layer)
    function hist(m) {
      const a = tf(m, m.min), b = tf(m, m.max), w = b > a ? (b - a) / NBARS : 1, c = new Uint32Array(NBARS);
      for (let i = 0; i < m.v.length; i++) c[Math.min(NBARS - 1, Math.max(0, Math.floor((tf(m, m.v[i]) - a) / w)))]++;
      m.bars = c; m.barW = w; m.a = a;
    }
    function drawHist(m) {
      const cv = m.el && m.el.querySelector(".bq-hist");
      if (!cv || !m.bars) return;
      const W = cv.clientWidth || 260, H = cv.clientHeight || 56, dpr = window.devicePixelRatio || 1;
      cv.width = Math.round(W * dpr); cv.height = Math.round(H * dpr);
      const g = cv.getContext("2d"); g.setTransform(dpr, 0, 0, dpr, 0, 0); g.clearRect(0, 0, W, H);
      let top = 0; for (const x of m.bars) top = Math.max(top, Math.sqrt(x));
      const t = Q.th[m.key] || {}, lo = t.min != null ? tf(m, t.min) : -Infinity, hi = t.max != null ? tf(m, t.max) : Infinity;
      const bw = W / NBARS;
      for (let k = 0; k < NBARS; k++) {
        if (!m.bars[k]) continue;
        const h = Math.max(1, (Math.sqrt(m.bars[k]) / (top || 1)) * (H - 2));
        const x0 = m.a + k * m.barW, x1 = x0 + m.barW;
        // a bar is drawn as kept if any of it lies inside the thresholds
        const kept = x1 > lo && x0 <= hi;
        g.fillStyle = kept ? "#4da3ff" : "#4b5263";
        g.fillRect(k * bw + 0.5, H - h, Math.max(1, bw - 1), h);
      }
      g.strokeStyle = "#ffd84d"; g.lineWidth = 1.5;
      const xOf = (v) => ((tf(m, v) - m.a) / (m.barW * NBARS)) * W;
      for (const v of [t.min, t.max]) if (v != null) { const x = Math.max(1, Math.min(W - 1, xOf(v))); g.beginPath(); g.moveTo(x, 0); g.lineTo(x, H); g.stroke(); }
    }
    function syncMetric(m) {
      const el = m.el; if (!el) return;
      const t = Q.th[m.key] || { min: null, max: null };
      el.querySelector(".bq-lo").value = t.min === null ? 0 : Math.max(1, Math.min(STEPS - 1, posOf(m, t.min)));
      el.querySelector(".bq-hi").value = t.max === null ? STEPS : Math.max(1, Math.min(STEPS - 1, posOf(m, t.max)));
      const lon = el.querySelector(".bq-lon"), hin = el.querySelector(".bq-hin");
      if (document.activeElement !== lon) lon.value = t.min === null ? "" : String(t.min);
      if (document.activeElement !== hin) hin.value = t.max === null ? "" : String(t.max);
      drawHist(m);
    }

    // ---- the keep mask ------------------------------------------------------------------
    function schedule() {
      if (Q.raf) return;
      Q.raf = requestAnimationFrame(() => { Q.raf = 0; recompute(false); });
    }
    function recompute(quiet) {
      const n = Q.n, keep = new Uint8Array(n).fill(1), per = {};
      for (const m of Q.metrics) {
        const t = Q.th[m.key]; per[m.key] = 0;
        if (!t) continue;
        const lo = t.min === null ? -Infinity : t.min, hi = t.max === null ? Infinity : t.max, v = m.v;
        let r = 0;
        for (let i = 0; i < n; i++) { const x = v[i]; if (x < lo || x > hi) { r++; keep[i] = 0; } }
        per[m.key] = r;
      }
      let removed = 0; for (let i = 0; i < n; i++) removed += 1 - keep[i];
      Q.keep = removed ? keep : null; Q.removed = removed; Q.perMetric = per; Q.layerData = null; Q.ver++;
      summary();
      syncUrl();
      if (Q.removed && Q.grey) ensurePositions();
      ctx.requestRender();
      if (!quiet || removed) ctx.onChange();
    }
    function summary() {
      const n = Q.n, r = Q.removed, pct = (x) => (n ? ((100 * x) / n).toFixed(1) : "0");
      Q.metrics.forEach((m) => {
        const t = Q.th[m.key], cut = m.el && m.el.querySelector(".bq-cut");
        if (!cut) return;
        cut.textContent = t ? `${t.min !== null ? "≥ " + fmt(m, t.min) : ""}${t.min !== null && t.max !== null ? " · " : ""}${t.max !== null ? "≤ " + fmt(m, t.max) : ""}` +
          ` · removes ${Q.perMetric[m.key].toLocaleString()}` : "no filter";
        cut.classList.toggle("on", !!t);
      });
      $("bqCount").textContent = r ? `(${pct(r)}% removed)` : "";
      $("bqSummary").innerHTML = say(r
        ? `Kept <b>${(n - r).toLocaleString()}</b> of ${n.toLocaleString()} bins · removed <b>${r.toLocaleString()}</b> (${pct(r)}%). ` +
          `Removed bins are grey on the slide and left out of Draw a region the next time it runs.`
        : Object.keys(Q.th).length
          ? `Kept all ${n.toLocaleString()} bins: no bin is outside these thresholds.`
          : `Off: all ${n.toLocaleString()} bins are kept. Move a slider (or type a number) to remove bins from Draw a region.`);
      $("bqReset").disabled = !Object.keys(Q.th).length;
    }

    // ---- the grey layer ------------------------------------------------------------------
    async function ensurePositions() {
      if (Q.pos) return;
      const key = Q.key;
      try {
        const p = await ctx.positions();
        if (Q.key === key && p && p.length === 2 * Q.n) { Q.pos = p; Q.layerData = null; ctx.requestRender(); }
      } catch (e) { /* no grey layer; the filter still works */ }
    }
    function layer(radius, opacity) {
      if (!Q.keep || !Q.grey || !Q.pos) return null;
      if (!Q.layerData) {
        const xy = new Float32Array(2 * Q.removed);
        let k = 0;
        for (let i = 0; i < Q.n; i++) if (!Q.keep[i]) { xy[2 * k] = Q.pos[2 * i]; xy[2 * k + 1] = Q.pos[2 * i + 1]; k++; }
        Q.layerData = { length: Q.removed, attributes: { getPosition: { value: xy, size: 2 } } };
      }
      return new ctx.ScatterplotLayer({
        id: `binqc-${Q.ds}-${Q.key}-${Q.ver}`, data: Q.layerData,
        radiusUnits: "common", getRadius: radius, radiusMinPixels: 0.8, stroked: false, filled: true,
        opacity: opacity == null ? 1 : opacity,      // cell outlines (outlines.js) take over from the dots when zoomed in
        getFillColor: GREY, parameters: ctx.NORMAL,
      });
    }

    // ---- state: Save JSON, the URL --------------------------------------------------------
    // what a result was computed with; null = no filter (or no QC values in this export)
    function state() {
      if (!Q.meta || !Object.keys(Q.th).length) return null;
      return { layer: Q.key, thresholds: JSON.parse(JSON.stringify(Q.th)), nRows: Q.n, kept: Q.n - Q.removed, removed: Q.removed,
               rule: "a bin is kept when every metric is >= its min and <= its max (null = no bound)" };
    }
    // null / missing = no filter. Returns a note when the thresholds cannot be applied.
    function setState(o) {
      const th = o && o.thresholds && typeof o.thresholds === "object" ? o.thresholds : {};
      if (!Q.meta) {
        return Object.keys(th).length ? "These shapes were run with a QC filter, but this export has no per-bin QC values: running without it." : "";
      }
      Q.th = cleanTh(th);
      Q.metrics.forEach(syncMetric);
      recompute(false);
      return o && o.layer && o.layer !== Q.key ? `The QC thresholds were set on the ${o.layer} layer; applied to ${Q.key}.` : "";
    }
    function syncUrl() {
      const u = new URL(location.href), s = paramOf();
      if (s) u.searchParams.set("qc", s); else u.searchParams.delete("qc");
      if (u.href !== location.href) history.replaceState(null, "", u);
    }
    function paramOf() {
      if (!Q.meta) return pending && pending.ds === Q.ds ? pending.raw : "";   // not applied yet: keep it in the link
      const parts = Q.metrics.filter((m) => Q.th[m.key]).map((m) => { const t = Q.th[m.key]; return `${m.key}=${t.min === null ? "" : t.min}-${t.max === null ? "" : t.max}`; });
      return parts.length ? `${Q.key}:${parts.join(",")}` : "";
    }
    function parseParam(s, ds) {
      const m = s && /^([\w.-]+):(.+)$/.exec(s);
      if (!m) return null;
      const th = {};
      for (const p of m[2].split(",")) {
        const q = /^(\w+)=([0-9.]*)-([0-9.]*)$/.exec(p.trim());
        if (!q) continue;
        th[q[1]] = { min: q[2] === "" ? null : Number(q[2]), max: q[3] === "" ? null : Number(q[3]) };
      }
      return { ds, key: m[1], th, raw: s };
    }

    function init() {
      $("bqReset").onclick = () => { Q.th = {}; Q.metrics.forEach(syncMetric); recompute(false); };
      $("bqGrey").onchange = (e) => { Q.grey = greyPref = e.target.checked; if (Q.grey && Q.removed) ensurePositions(); ctx.requestRender(); };
      window.addEventListener("resize", () => Q.metrics.forEach(drawHist));
      window.addEventListener("soda:shown", () => Q.metrics.forEach(drawHist));   // layout.js: the section was folded or in another tab when drawn
    }
    init();

    return {
      reset, onLayer, layer, state, setState,
      // the mask for the worker: one byte per row, 1 = keep; null when nothing is removed
      keepMask: () => (Q.keep && Q.key === ctx.layerKey() ? Q.keep : null),
      // the same mask while "Grey out removed bins" is on (what outlines.js greys), and a number that changes with it
      greyMask: () => (Q.keep && Q.grey && Q.key === ctx.layerKey() ? Q.keep : null),
      ver: () => Q.ver,
      removed: () => (Q.key === ctx.layerKey() ? Q.removed : 0),
      available: () => !!Q.meta,
      describe: () => {
        if (!Q.removed) return "";
        const t = Q.metrics.filter((m) => Q.th[m.key]).map((m) => { const x = Q.th[m.key];
          return `${m.label.toLowerCase()} ${x.min !== null ? "≥ " + fmt(m, x.min) : ""}${x.min !== null && x.max !== null ? ", " : ""}${x.max !== null ? "≤ " + fmt(m, x.max) : ""}`; });
        return t.join("; ");
      },
      _debug: () => Q,
    };
  }

  window.SodaBinQc = { create };
})();
