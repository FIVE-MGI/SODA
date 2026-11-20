// Circle -> UMAP, off the main thread. The viewer hands over the hvg/ chunks under the
// drawn shapes plus a row -> shape map; this returns a 2-D embedding of the selected bins
// (or of one summed profile per shape). Steps, all seeded so a rerun is identical:
//   gather -> drop low-UMI bins -> cap (stratified sample) -> log1p CP10k ->
//   the section's variable genes counted in enough of the bins -> scale ->
//   randomized PCA (PCs by bin count) ->
//   umap-js started from PC1 / PC2 -> Leiden on the UMAP's fuzzy neighbour graph, the
//   resolution chosen by stability unless given -> marker genes per cluster.
// Marker squares (protein grid, msg.protein) take a shorter road through runProtein: the
// mat/ chunks' raw means -> x with the whole slide's background and p99.5 (asinh; older
// exports: x' with the 1st..99.5th percentile window) -> z-score within the selection ->
// PCA -> the same UMAP, Leiden, markers and suggested names.
// Speed (the same results as before, bit for bit on one browser): umapfit.js runs umap-js's fit over typed
// arrays; the Leiden runs of the stability sweep go to helper workers (leiden_worker.js) when the browser
// has cores to spare, and run while the layout's epochs do (fitUmapEpochs); the last run's layout is kept, so the same shapes again with other clustering
// settings only re-cluster. Fast mode (params.sampling, above params.fastAbove bins): the UMAP and Leiden
// of a seeded, shape-stratified sample of params.fastSample bins, every other bin placed by its 15 nearest
// sample bins in PCA space (layout, below); markers and the result cover every bin.
/* global UMAP, leiden, SODA_LIBS, SodaUmapFit */
importScripts("libs.js?v=20261119s3");
importScripts("leiden.js?v=20261102a");
importScripts("umapfit.js?v=20261117b");
const LEIDEN_WORKER = "leiden_worker.js?v=20261117a";
const now = () => performance.now();
// umap-js (pinned in libs.js): the copy in vendor/ when it is there and matches its SHA-384, else the CDN's,
// as before. importScripts has no integrity check of its own, so fetch() checks the hash first; the script
// then comes from the browser's cache. Loaded before the first job.
const umapReady = (async () => {
  const L = SODA_LIBS.umap;
  for (const [url, opts] of [[L.local, {}], [L.cdn, { mode: "cors", credentials: "omit" }]]) {
    try {
      // vendor/ first asked without the hash: when it is not there, WebKit reports a failed check of the hash
      // ("access control checks") as an error rather than the plain 404. No answer at all is not "not there":
      // the hash-checked GET below may still come from the browser's cache.
      if (url === L.local) {
        let h = null;
        try { h = await fetch(url, { method: "HEAD" }); } catch (e) { /* offline */ }
        if (h && !h.ok) continue;
      }
      const r = await fetch(url, { ...opts, integrity: L.sri });
      if (!r.ok) continue;
      await r.arrayBuffer();
      importScripts(url);
      return url === L.local ? "vendor" : "cdn";
    } catch (e) { /* the next one */ }
  }
  importScripts(L.cdn);                        // the last try: its error says why (offline, blocked)
  return "cdn";
})();
umapReady.catch(() => {});                     // reported with the first job

const NO_SHAPE = 0xffff;                       // rowShape value for a bin in no shape

function mulberry32(seed) {
  let s = seed >>> 0;
  return () => {
    s = (s + 0x6d2b79f5) | 0;
    let t = Math.imul(s ^ (s >>> 15), 1 | s);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}
function gauss(rand) {                         // Box-Muller
  let u = 0; while (u === 0) u = rand();
  return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * rand());
}
const progress = (step, frac) => postMessage({ type: "progress", step, frac });

// ---- gather the selected bins out of the chunks ----------------------------------
function gather(chunks, rowShape) {
  const rows = [], shape = [], lib = [], ptr = [0], gi = [], ct = [];
  for (const buf of chunks) {
    const dv = new DataView(buf);
    const n = dv.getUint32(0, true);
    let o = 4;
    const row = new Uint32Array(buf.slice(o, o + 4 * n)); o += 4 * n;
    const libs = new Float32Array(buf.slice(o, o + 4 * n)); o += 4 * n;
    const ip = new Uint32Array(buf.slice(o, o + 4 * (n + 1))); o += 4 * (n + 1);
    const nnz = ip[n];
    const g = new Uint16Array(buf.slice(o, o + 2 * nnz)); o += 2 * nnz;
    const c = new Uint16Array(buf.slice(o, o + 2 * nnz));
    for (let i = 0; i < n; i++) {
      const s = rowShape[row[i]];
      if (s === NO_SHAPE) continue;
      rows.push(row[i]); shape.push(s); lib.push(libs[i]);
      for (let k = ip[i]; k < ip[i + 1]; k++) { gi.push(g[k]); ct.push(c[k]); }
      ptr.push(gi.length);
    }
  }
  return { rows, shape, lib, ptr, gi, ct };
}

// keep[] of bin indices -> a new compact CSR
function subset(d, keep) {
  const out = { rows: [], shape: [], lib: [], ptr: [0], gi: [], ct: [] };
  for (const i of keep) {
    out.rows.push(d.rows[i]); out.shape.push(d.shape[i]); out.lib.push(d.lib[i]);
    for (let k = d.ptr[i]; k < d.ptr[i + 1]; k++) { out.gi.push(d.gi[k]); out.ct.push(d.ct[k]); }
    out.ptr.push(out.gi.length);
  }
  return out;
}

// Proportional per shape, but never fewer than min(n_s, floor) so a small shape is not
// drowned by a big one. Seeded partial Fisher-Yates inside each shape.
function stratifiedSample(shapeOf, nShapes, cap, rand) {
  const by = Array.from({ length: nShapes }, () => []);
  shapeOf.forEach((s, i) => by[s].push(i));
  const N = shapeOf.length, floor = Math.min(500, Math.floor(cap / Math.max(1, nShapes)));
  let quota = by.map((b) => Math.min(b.length, Math.max(Math.min(b.length, floor), Math.floor(cap * b.length / N))));
  let over = quota.reduce((a, b) => a + b, 0) - cap;
  while (over > 0) {                              // trim the largest quotas until it fits
    const j = quota.indexOf(Math.max(...quota));
    const t = Math.min(over, quota[j] - Math.min(by[j].length, floor) || 1);
    quota[j] -= t; over -= t;
  }
  const keep = [];
  by.forEach((b, s) => {
    const a = b.slice();
    for (let i = 0; i < quota[s]; i++) {
      const j = i + Math.floor(rand() * (a.length - i));
      [a[i], a[j]] = [a[j], a[i]];
      keep.push(a[i]);
    }
  });
  return keep.sort((x, y) => x - y);
}

// ---- clustering ----------------------------------------------------------------
// The fuzzy simplicial set (scanpy's "connectivities"), as umap-js keeps it (umapfit.js graph: entries in
// umap-js's Map order). Symmetrise to CSR: w_ij = (a_ij + a_ji) / 2 (the fuzzy union is already
// symmetric, so this only guards), zeros and self loops dropped.
function fuzzyGraphCSR(G, N) {
  const R = G.rows, C = G.cols, V = G.vals, m = G.nnz, deg = new Uint32Array(N + 1);
  for (let e = 0; e < m; e++) if (V[e] > 0 && R[e] !== C[e]) { deg[R[e] + 1]++; deg[C[e] + 1]++; }
  for (let i = 0; i < N; i++) deg[i + 1] += deg[i];
  const fill = deg.slice(0, N), col = new Uint32Array(deg[N]), val = new Float64Array(deg[N]);
  for (let e = 0; e < m; e++) {
    const v = V[e], r = R[e], c = C[e];
    if (!(v > 0) || r === c) continue;
    col[fill[r]] = c; val[fill[r]++] = v / 2;
    col[fill[c]] = r; val[fill[c]++] = v / 2;
  }
  // merge the duplicate (r, c) pairs inside each row
  const indptr = new Uint32Array(N + 1), at = new Float64Array(N).fill(-1);
  const indices = new Uint32Array(deg[N]), weights = new Float64Array(deg[N]);
  let n = 0;
  for (let i = 0; i < N; i++) {
    const start = n;
    for (let p = deg[i]; p < deg[i + 1]; p++) {
      const j = col[p];
      if (at[j] >= start) weights[at[j]] += val[p];
      else { at[j] = n; indices[n] = j; weights[n++] = val[p]; }
    }
    indptr[i + 1] = n;
  }
  return { indptr, indices: indices.slice(0, n), weights: weights.slice(0, n) };
}

// PCs and neighbours by the number of bins: a few hundred bins carry only a few real axes,
// and at 30 PCs the rest is noise the neighbour graph would follow.
function sizeRule(N) {
  return N < 1000 ? { pcs: 15, k: 10 } : N <= 5000 ? { pcs: 20, k: 15 } : { pcs: 30, k: 15 };
}

// umap-js starts from uniform random positions; start from PC1 / PC2 instead, as umap-learn's
// init="pca" in sectionumap.py does: each axis min-max scaled to [0, 10] plus a 1e-4 jitter so
// identical bins do not sit on one point. The optimiser holds umap.embedding by reference, so
// the start is written into it in place after initializeFit (a new array would be ignored).
// umap-js's fit itself runs as umapfit.js (the same numbers over typed arrays); the returned fit holds
// the embedding flat (point i at embedding[2i], [2i + 1]) and the fuzzy graph. T: stage times (ms).
function fitUmap(scores, N, k, minDist, rand, T) {
  const fit = fitUmapStart(scores, N, k, minDist, rand, T), nEpochs = fit.nEpochs;
  const t = now();
  for (let e = 0; e < nEpochs; e++) {
    fit.step(e);
    if (e % 10 === 0) progress("UMAP", 0.2 + 0.8 * e / nEpochs);
  }
  if (T) T.layout = now() - t;
  return fit;
}
// the neighbours, the graph and the start (T gets knn, graph, layoutPrep); the epochs are the caller's
function fitUmapStart(scores, N, k, minDist, rand, T) {
  const D = scores[0].length, X = new Float64Array(N * D);
  for (let i = 0; i < N; i++) { const r = scores[i]; for (let d = 0; d < D; d++) X[i * D + d] = r[d]; }
  const [a, b] = SodaUmapFit.abParams(UMAP, minDist);
  const fit = SodaUmapFit.init(X, N, D, { nNeighbors: Math.max(2, Math.min(k, N - 1)), random: rand, a, b, timing: T });
  const nEpochs = fit.nEpochs, E = fit.embedding;
  for (let a = 0; a < 2; a++) {
    let lo = Infinity, hi = -Infinity;
    for (let i = 0; i < N; i++) { const v = scores[i][a] || 0; if (v < lo) lo = v; if (v > hi) hi = v; }
    const s = hi > lo ? 10 / (hi - lo) : 0;
    for (let i = 0; i < N; i++) E[2 * i + a] = ((scores[i][a] || 0) - lo) * s + 1e-4 * gauss(rand);
  }
  return fit;
}

// The Leiden sweep needs only the graph, not the layout, so a run with clustering starts it on the helper
// workers before the layout's epochs and lays out meanwhile: the same epochs, edge for edge in the same
// order (umapfit.js edges / endEpoch), handing this thread back every ~25 ms while the sweep is pending
// so the helpers' results are taken and their next runs given out. The numbers are those of fitUmap;
// only the waiting overlaps. side: the sweep's promise (or null: fitUmap's loop, no hand-backs).
const sideWork = { on: false, pool: false, step: null, frac: 0 };
const yieldTurn = (() => {
  let ch = null; const q = [];
  return () => new Promise((r) => {
    if (typeof MessageChannel === "undefined") { setTimeout(r, 0); return; }
    if (!ch) { ch = new MessageChannel(); ch.port1.onmessage = () => q.shift()(); if (ch.port1.unref) ch.port1.unref(); }
    q.push(r); ch.port2.postMessage(0);
  });
})();
async function fitUmapEpochs(fit, T, side) {
  const nEpochs = fit.nEpochs, nE = fit.nEdges, SLICE = 4096;
  // without helpers (a small run, few cores) the sweep has already run, here, before the epochs: no hand-backs
  let pending = !!side && sideWork.pool, last = now();
  if (side) side.then(() => (pending = false), () => (pending = false));
  const t = now();
  for (let e = 0; e < nEpochs; e++) {
    if (!pending) fit.step(e);
    else {
      for (let i = 0; i < nE; i += SLICE) {
        fit.edges(e, i, Math.min(nE, i + SLICE));
        if (pending && now() - last > 25) { await yieldTurn(); last = now(); }
      }
      fit.endEpoch(e);
    }
    if (e % 10 === 0) progress(pending && sideWork.step ? `UMAP, and alongside it ${sideWork.step.toLowerCase()}` : "UMAP", 0.2 + 0.8 * e / nEpochs);
  }
  if (T) T.layout = now() - t;
  return fit;
}

// "auto" = choose the resolution by stability; a number = that resolution (its stability is
// still measured); 0 = no clusters
function wantedResolution(msg, P) {
  const r = msg.resolution !== undefined ? msg.resolution : P.resolution !== undefined ? P.resolution : "auto";
  return r === "auto" ? "auto" : Number(r) || 0;
}

// ---- how stable is a clustering ------------------------------------------------------
// Leiden at each resolution on the run's fixed UMAP graph, with several seeds on the whole
// graph and on seeded 80% induced subgraphs (bins dropped, graph not rebuilt, so this is
// optimistic: "stability", not significance). stability = median adjusted Rand index over
// every pair of a whole-graph run with another run (subgraph runs compared on their own bins).
// The clusters shown are the medoid whole-graph run (highest mean ARI to the others), so one
// unlucky seed neither decides the score nor the picture; per cluster, the median best
// Jaccard overlap of the medoid's cluster with any cluster of each other run. "auto" takes,
// among resolutions with > 1 cluster and stability >= minStable, the most stable one; a finer
// resolution wins a tie (within tieStable), so of several equally reproducible answers the more
// detailed is shown. (Taking simply the finest stable one over-split planted populations in
// the synthetic tests: ARI 0.65-0.79 against 1.0 for the most stable.) If none qualifies, the coarsest resolution is shown, flagged, with the closest
// candidate named: on a homogeneous patch the honest answer is one population, not the
// least-bad dozen. The sweep stops after two resolutions in a row below stopBelow: on the
// test lassos stability never climbed back once it had fallen that far, and the finest
// resolutions cost the most. Leiden runs 2 iterations, as scanpy and sectionumap.py do. Same
// graph and seed give the same answer.
const CLR_EPS = 0.1, LOW_SIGNAL_X = 0.2, CLR_MIN_CHANNELS = 10;   // = protgrid.CLR_EPS, LOW_SIGNAL_X, CLR_MIN_CHANNELS
const SWEEP = { res: [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2], frac: 0.8, minStable: 0.8, tieStable: 0.02, stopBelow: 0.6, iters: 2 };
const sweepRuns = (N) => (N <= 20000 ? { seeds: 10, subs: 10 } : { seeds: 5, subs: 5 });

// contingency of two labelings (labels 0..n-1) as a dense table
function contingency(a, b) {
  let na = 0, nb = 0;
  for (let i = 0; i < a.length; i++) { if (a[i] >= na) na = a[i] + 1; if (b[i] >= nb) nb = b[i] + 1; }
  const t = new Uint32Array(na * nb), ra = new Float64Array(na), rb = new Float64Array(nb);
  for (let i = 0; i < a.length; i++) { t[a[i] * nb + b[i]]++; ra[a[i]]++; rb[b[i]]++; }
  return { t, ra, rb, na, nb };
}

function ari(a, b) {
  const { t, ra, rb } = contingency(a, b), c2 = (m) => (m * (m - 1)) / 2;
  let sij = 0, sa = 0, sb = 0;
  for (let i = 0; i < t.length; i++) if (t[i] > 1) sij += c2(t[i]);
  ra.forEach((v) => (sa += c2(v))); rb.forEach((v) => (sb += c2(v)));
  const exp = (sa * sb) / c2(a.length), max = (sa + sb) / 2;
  return max === exp ? 1 : (sij - exp) / (max - exp);
}

// for each cluster c of ref: the best Jaccard |c & d| / |c | d| over the clusters d of lab
function bestJaccard(ref, lab, nCl, into) {
  const { t, ra, rb, na, nb } = contingency(ref, lab);
  for (let c = 0; c < Math.min(na, nCl); c++) {
    if (!ra[c]) continue;
    let best = 0;
    for (let d = 0; d < nb; d++) { const n = t[c * nb + d]; if (n) best = Math.max(best, n / (ra[c] + rb[d] - n)); }
    into[c].push(best);
  }
}

// the graph restricted to the bins in keep (ascending), renumbered 0..keep.length-1 (= leiden_worker.js)
function inducedGraph(g, keep, N) {
  const NONE = 0xFFFFFFFF, at = new Uint32Array(N).fill(NONE); keep.forEach((v, j) => (at[v] = j));
  const indptr = new Uint32Array(keep.length + 1), indices = new Uint32Array(g.indices.length), weights = new Float64Array(g.indices.length);
  let n = 0;
  for (let j = 0; j < keep.length; j++) {
    const v = keep[j];
    for (let p = g.indptr[v]; p < g.indptr[v + 1]; p++) { const u = at[g.indices[p]]; if (u !== NONE) { indices[n] = u; weights[n++] = g.weights[p]; } }
    indptr[j + 1] = n;
  }
  return { indptr, indices: indices.slice(0, n), weights: weights.slice(0, n) };
}

const median = (xs) => { if (!xs.length) return null; const s = xs.slice().sort((a, b) => a - b), m = s.length >> 1; return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2; };

// The runs of one resolution are independent: with cores to spare they go to helper workers
// (leiden_worker.js, the same leiden.js), each holding the graph and the subsamples. Which worker runs
// which job does not matter: every run is the same function of (graph, resolution, seed). Up to 8
// helpers, one core left for the page (and one more for the layout when the sweep runs alongside it,
// fitUmapEpochs); none when that leaves fewer than 2 or the browser cannot start a worker from a worker,
// and any failure falls back to running them here. Below POOL_MIN bins none either: the whole sweep of a
// 1,000-bin run takes ~1.4 s here, and starting 8 helpers cost more than they saved on a CPU slowed 4x
// (17.7 s against 13.7 s, Leiden alongside the layout).
const POOL_MIN = 2000;
function leidenPool(g, N, keeps) {
  const cores = (self.navigator && navigator.hardwareConcurrency) || 1, n = Math.min(8, cores - (sideWork.on ? 2 : 1));
  if (n < 2 || N < POOL_MIN || typeof Worker === "undefined") return null;
  try {
    const ws = [];
    for (let i = 0; i < n; i++) {
      const w = new Worker(LEIDEN_WORKER);
      w.postMessage({ type: "graph", n: N, indptr: g.indptr, indices: g.indices, weights: g.weights, keeps });
      ws.push(w);
    }
    return ws;
  } catch (e) { return null; }
}
function onPool(ws, jobs) {
  return new Promise((resolve, reject) => {
    const out = new Array(jobs.length);
    let next = 0, done = 0;
    const give = (w) => {
      if (next >= jobs.length) return;
      const id = next++;
      w.onmessage = (ev) => {
        if (ev.data.error) { reject(new Error(ev.data.error)); return; }
        out[id] = ev.data.lab;
        if (++done === jobs.length) resolve(out); else give(w);
      };
      w.onerror = (e) => { e.preventDefault && e.preventDefault(); reject(new Error(e.message || "Leiden helper failed")); };
      w.postMessage({ type: "run", ...jobs[id] });
    };
    ws.forEach(give);
  });
}

// clusters at the wanted resolution ("auto" or a number) with their stability; fills info
// fullSweep (fast mode): every resolution, no early stop. A sample's coarse resolutions can be unstable (its
// small populations are thinner) while finer ones are not: on a 50,000-bin region the 10,000-bin sample's
// rows at 0.1 and 0.2 scored 0.35 and 0.50 (the stop) and 0.3-0.4 scored 0.90; the sweep of a sample is cheap.
async function clusterGraph(g, N, resolution, seed, info, fullSweep) {
  const list = resolution === "auto" ? SWEEP.res : [resolution], R = sweepRuns(N), r3 = (x) => Math.round(x * 1000) / 1000;
  const rows = [];
  // the subsampled runs' bins: pick restarts from the same seed at every resolution, so they are the
  // same bins (and induced graphs) each time
  const pick = mulberry32((seed ^ 0x9e3779b9) >>> 0), keeps = [];
  for (let s = 0; s < R.subs; s++) {
    const keep = []; for (let i = 0; i < N; i++) if (pick() < SWEEP.frac) keep.push(i);
    keeps.push(Uint32Array.from(keep));
  }
  const subG = [];
  const runHere = (j) => {
    if (j.sub < 0) return leiden(N, g.indptr, g.indices, g.weights, { resolution: j.res, seed: j.seed, maxIter: SWEEP.iters });
    const sub = subG[j.sub] || (subG[j.sub] = inducedGraph(g, keeps[j.sub], N));
    return leiden(keeps[j.sub].length, sub.indptr, sub.indices, sub.weights, { resolution: j.res, seed: j.seed, maxIter: SWEEP.iters });
  };
  let pool = leidenPool(g, N, keeps);
  if (sideWork.on) sideWork.pool = !!pool;
  info.leidenWorkers = pool ? pool.length : 0;
  try {
  for (const [li, res] of list.entries()) {
    const last2 = rows.slice(-2);
    if (!fullSweep && last2.length === 2 && last2.every((r) => r.stability < SWEEP.stopBelow)) break;
    const step = `Clustering, resolution ${res}`;
    if (sideWork.on && pool) { sideWork.step = step; sideWork.frac = li / list.length; } else progress(step, li / list.length);
    const jobs = [];
    for (let s = 0; s < R.seeds; s++) jobs.push({ res, seed: seed + s, sub: -1, iters: SWEEP.iters });
    for (let s = 0; s < R.subs; s++) jobs.push({ res, seed: seed + s, sub: s, iters: SWEEP.iters });
    let labs = null;
    if (pool) {
      try { labs = await onPool(pool, jobs); }
      catch (e) { pool.forEach((w) => w.terminate()); pool = null; info.leidenWorkers = 0; }
    }
    if (!labs) labs = jobs.map(runHere);
    const full = labs.slice(0, R.seeds), subs = [];
    for (let s = 0; s < R.subs; s++) subs.push({ keep: keeps[s], lab: labs[R.seeds + s] });
    // every whole-graph run against every other run; sc[a] = run a's scores
    const sc = full.map(() => []), scores = [];
    for (let a = 0; a < full.length; a++) {
      for (let b = a + 1; b < full.length; b++) { const x = ari(full[a], full[b]); sc[a].push(x); sc[b].push(x); scores.push(x); }
      for (const s of subs) { const x = ari(Uint16Array.from(s.keep, (i) => full[a][i]), s.lab); sc[a].push(x); scores.push(x); }
    }
    const mean = (xs) => xs.reduce((p, q) => p + q, 0) / (xs.length || 1);
    let m = 0; for (let a = 1; a < full.length; a++) if (mean(sc[a]) > mean(sc[m])) m = a;
    const ref = full[m];
    let nCl = 0; for (let i = 0; i < N; i++) if (ref[i] + 1 > nCl) nCl = ref[i] + 1;
    const match = Array.from({ length: nCl }, () => []);
    full.forEach((lab, a) => { if (a !== m) bestJaccard(ref, lab, nCl, match); });
    subs.forEach((s) => bestJaccard(Uint16Array.from(s.keep, (i) => ref[i]), s.lab, nCl, match));
    rows.push({ res, nCl, ref, stability: median(scores), perCluster: match.map(median) });
  }
  } finally { if (pool) pool.forEach((w) => w.terminate()); }
  const multi = rows.filter((r) => r.nCl > 1), ok = multi.filter((r) => r.stability >= SWEEP.minStable);
  const maxStable = ok.reduce((m, r) => Math.max(m, r.stability), -1);
  const top = ok.filter((r) => r.stability >= maxStable - SWEEP.tieStable).pop();   // rows run coarse -> fine
  const best = resolution !== "auto" ? rows[0] : top || rows[0];
  const closest = multi.length ? multi.reduce((a, b) => (b.stability > a.stability ? b : a)) : null;
  const sizes = new Array(best.nCl).fill(0); best.ref.forEach((c) => sizes[c]++);
  info.nClusters = best.nCl; info.clusterSizes = sizes; info.resolution = best.res;
  info.autoResolution = resolution === "auto"; info.stability = r3(best.stability);
  info.stable = best.nCl > 1 ? best.stability >= SWEEP.minStable : false;
  info.closest = resolution === "auto" && !ok.length && closest
    ? { resolution: closest.res, nClusters: closest.nCl, stability: r3(closest.stability) } : null;
  info.clusterStability = best.perCluster.map((x) => (x === null ? null : r3(x)));
  info.sweep = rows.map((r) => ({ resolution: r.res, nClusters: r.nCl, stability: r3(r.stability) }));
  info.stabilityRuns = { seeds: R.seeds, subsamples: R.subs, fraction: SWEEP.frac, minStable: SWEEP.minStable, tieStable: SWEEP.tieStable,
                         stopBelow: fullSweep ? null : SWEEP.stopBelow, leidenIterations: SWEEP.iters, shown: "medoid run" };
  return best.ref;
}

// Top 10 genes per cluster by log2 fold change (scanpy's convention: log2 of expm1(mean
// log1p CP10k) in the cluster over the rest of the selection), among the genes used for the
// layout that are detected in >= 25% of the cluster's bins (lfc > 0 only, so a cluster
// can list fewer than 10).
function clusterMarkers(d, val, genes, colOf, cl, nCl, names) {
  const N = cl.length, K = genes.length;
  const sumIn = new Float64Array(nCl * K), nzIn = new Uint32Array(nCl * K), tot = new Float64Array(K), size = new Uint32Array(nCl);
  for (let i = 0; i < N; i++) {
    const c = cl[i]; size[c]++;
    for (let k = d.ptr[i]; k < d.ptr[i + 1]; k++) {
      const j = colOf[d.gi[k]]; if (j < 0 || !(val[k] > 0)) continue;
      sumIn[c * K + j] += val[k]; nzIn[c * K + j]++; tot[j] += val[k];
    }
  }
  const r3 = (x) => Math.round(x * 1000) / 1000;
  const out = [];
  for (let c = 0; c < nCl; c++) {
    const nIn = size[c], nOut = N - nIn, rows = [];
    for (let j = 0; j < K; j++) {
      const pct = nzIn[c * K + j] / nIn;
      if (pct < 0.25) continue;
      const mIn = sumIn[c * K + j] / nIn, mOut = nOut ? (tot[j] - sumIn[c * K + j]) / nOut : 0;
      const lfc = Math.log2((Math.expm1(mIn) + 1e-9) / (Math.expm1(Math.max(0, mOut)) + 1e-9));
      if (lfc > 0) rows.push([j, lfc, pct]);              // up in the cluster only
    }
    rows.sort((a, b) => b[1] - a[1] || a[0] - b[0]);
    out.push(rows.slice(0, 10).map(([j, lfc, pct]) => ({ gene: names ? names[genes[j]] : genes[j], gi: genes[j], lfc: r3(lfc), pct: r3(pct) })));
  }
  return out;
}

// ---- suggested cell types from marker lists ----------------------------------------
// The same rule as celltypes.py (keep the two identical), on the run's bins:
//   a marker gene is missing as "absent" when it is not among the shipped genes (matched
//   case-insensitively), as "rare" when it is counted in < 0.5% of the bins (a single UMI of
//   a gene almost never seen would z-score so high it alone decided the bin; sd 0 is rare too);
//   x = log1p CP10k of each remaining gene; z = (x - mean) / sd over the bins (population sd);
//   score[bin, type] = mean z over the type's genes; a type with < 2 genes left is not scored.
//   A bin's type = its best score (first in list order on ties) if >= 0.25, else none.
//   Per cluster: pct[t] = share of its bins (none included) of type t; top / second by pct;
//   high: top >= 50% and >= 2x second; medium: >= 30%; low ("maybe"): >= 15%; else "Unclear".
const CT = { minGenes: 2, minDetect: 0.005, binMin: 0.25, high: 0.5, medium: 0.3, low: 0.15 };
function cellTypes(d, val, cl, nCl, set, names) {
  const N = cl.length, first = new Map();
  names.forEach((g, j) => { const k = String(g).toUpperCase(); if (!first.has(k)) first.set(k, j); });
  const types = set.types.map((t) => {
    const cols = [], genes = [], missing = [];
    for (const g of t.genes || []) {
      const j = first.get(String(g).toUpperCase());
      if (j === undefined) missing.push({ gene: g, reason: "absent" }); else { cols.push(j); genes.push(names[j]); }
    }
    return { name: t.name, cols, genes, missing };
  });
  // the marker genes as dense columns over the bins (two-pass mean / sd, as numpy)
  const want = new Map(), seen = new Map();
  for (const t of types) for (const j of t.cols) if (!want.has(j)) { want.set(j, new Float64Array(N)); seen.set(j, 0); }
  for (let i = 0; i < N; i++)
    for (let k = d.ptr[i]; k < d.ptr[i + 1]; k++) {
      const x = want.get(d.gi[k]); if (!x) continue;
      x[i] = val[k]; if (d.ct[k] > 0) seen.set(d.gi[k], seen.get(d.gi[k]) + 1);
    }
  const z = new Map();
  for (const [j, x] of want) {
    if (seen.get(j) < CT.minDetect * N) continue;            // rare
    let m = 0; for (let i = 0; i < N; i++) m += x[i];
    m /= N;
    let v = 0; for (let i = 0; i < N; i++) v += (x[i] - m) ** 2;
    const sd = Math.sqrt(v / N);
    if (sd > 0) { for (let i = 0; i < N; i++) x[i] = (x[i] - m) / sd; z.set(j, x); }
  }
  const T = types.length, score = new Float64Array(N * T), scored = types.map(() => false);
  types.forEach((t, ti) => {
    const cols = t.cols.filter((j) => z.has(j));
    t.used = t.genes.filter((_, k) => z.has(t.cols[k]));
    t.missing = t.missing.concat(t.genes.filter((_, k) => !z.has(t.cols[k])).map((g) => ({ gene: g, reason: "rare" })));
    if (cols.length < CT.minGenes) return;
    scored[ti] = true;
    for (const j of cols) { const x = z.get(j); for (let i = 0; i < N; i++) score[i * T + ti] += x[i]; }
    for (let i = 0; i < N; i++) score[i * T + ti] /= cols.length;
  });
  return ctSummarise(score, scored, types, cl, nCl, CT.binMin, set,
    (t) => `only ${t.used.length} of its ${t.used.length + t.missing.length} genes are measurable here (needs ${CT.minGenes})`);
}

// each bin's type (best score >= binMin), then per cluster the share of each type and the mean score
function ctSummarise(score, scored, types, cl, nCl, binMin, set, why) {
  const N = cl.length, T = types.length;
  const NONE = -1, cnt = new Float64Array(nCl * T), sum = new Float64Array(nCl * T), size = new Float64Array(nCl);
  for (let i = 0; i < N; i++) {
    const c = cl[i]; size[c]++;
    let best = NONE, bs = -Infinity;
    for (let t = 0; t < T; t++) {
      if (!scored[t]) continue;
      const s = score[i * T + t]; sum[c * T + t] += s;
      if (s > bs) { bs = s; best = t; }
    }
    if (best !== NONE && bs >= binMin) cnt[c * T + best]++;
  }
  const r3 = (x) => Math.round(x * 1000) / 1000;
  const clusters = [];
  for (let c = 0; c < nCl; c++) {
    const nb = size[c], pct = Array.from({ length: T }, (_, t) => cnt[c * T + t] / Math.max(nb, 1));
    const order = pct.map((_, t) => t).sort((a, b) => pct[b] - pct[a] || a - b);    // stable: list order on ties
    const top = order[0], sec = T > 1 ? order[1] : null, p1 = pct[top], p2 = sec !== null ? pct[sec] : 0;
    const confidence = p1 >= CT.high && p1 >= 2 * p2 ? "high" : p1 >= CT.medium ? "medium" : p1 >= CT.low ? "low" : "none";
    const nm = types[top].name;
    clusters.push({ cluster: c, suggest: confidence === "low" ? `maybe ${nm}` : confidence === "none" ? "Unclear" : nm,
                    top: nm, pct: r3(p1), second: sec !== null ? types[sec].name : null, pct2: r3(p2), confidence,
                    meanScore: types.map((_, t) => (scored[t] && nb ? r3(sum[c * T + t] / nb) : null)) });
  }
  return { markerSet: { name: set.name || null, version: set.version != null ? set.version : null },
           types: types.map((t, ti) => ({ name: t.name, genes: t.used, missing: t.missing, scored: scored[ti],
             ...(scored[ti] ? {} : { reason: why(t) }) })),
           clusters };
}

// ---- marker squares (protein grid) -------------------------------------------------
// The same rule as celltypes.write_celltypes_protein (keep the two identical), on the run's
// squares: z = (x - mean) / sd per channel (population sd; sd 0 = "flat"), mean and sd those of the
// whole slide when the export carries them (sec: mat/index.json mean / sd), else of the run's squares -
// so an all-immune circle still reads as immune, which against its own squares it cannot;
// score = mean z over the type's channels found (>= the type's "min", else the list's
// "minChannels", else 1); bin threshold 0.5; the cluster summary as for genes. Channel names
// are matched after normName, each list name also under its aliases (set.aliases).
// A cluster with more than half its squares flagged everything-bright (flags bit 1) is
// suggested as ARTEFACT.
const PROTEIN_CT = { minChannels: 1, binMin: 0.5, artefactShare: 0.5 };
const ARTEFACT = "Artefact? (fold / edge)";
// = celltypes.norm_name: a leading index, trailing fluorophore / FFPE tags, then case,
// spaces, "-", "_" and "." ignored
const NN_LEAD = /^\s*\d+\s*[_\-.:)]\s*/;
const NN_TAG = /(?:[\s_\-]+|\s*\()(?:opal\s*\d+|atto\s*\d+|af\s*\d+|alexa(?:\s*fluor)?\s*\d+|cy\s*\d(?:\.\d)?|fitc|tritc|pe|apc|ffpe|405|488|514|520|532|546|555|568|570|594|620|647|650|660|680|690|700|750|780|790)\)?\s*$/i;
function normName(s) {
  s = String(s).replace(NN_LEAD, "");
  for (;;) { const t = s.replace(NN_TAG, ""); if (t === s || !t.trim()) break; s = t; }
  return s.replace(/[\s_\-.]/g, "").toLowerCase();
}
const typeMin = (set, t) => Number(t.min != null ? t.min : set.minChannels != null ? set.minChannels : PROTEIN_CT.minChannels);
function cellTypesProtein(X, N, K, cl, nCl, set, names, flags, sec) {
  const first = new Map(), aliases = set.aliases || {};
  names.forEach((g, j) => { const k = normName(g); if (!first.has(k)) first.set(k, j); });
  const z = new Map();
  const zOf = (j) => {
    if (z.has(j)) return z.get(j);
    let m, sd;
    if (sec) { m = sec.mean[j]; sd = sec.sd[j]; }
    else {
      m = 0; for (let i = 0; i < N; i++) m += X[i * K + j];
      m /= N;
      let v = 0; for (let i = 0; i < N; i++) v += (X[i * K + j] - m) ** 2;
      sd = Math.sqrt(v / N);
    }
    let col = null;
    if (sd > 0) { col = new Float64Array(N); for (let i = 0; i < N; i++) col[i] = (X[i * K + j] - m) / sd; }
    z.set(j, col); return col;
  };
  const T = set.types.length, score = new Float64Array(N * T), scored = new Array(T).fill(false);
  const types = set.types.map((t, ti) => {
    const used = [], missing = [], cols = [], need = typeMin(set, t);
    for (const g of t.genes || []) {
      const k = [normName(g), ...(aliases[g] || []).map(normName)].find((x) => first.has(x));
      const j = k === undefined ? undefined : first.get(k);
      if (j === undefined || cols.includes(j)) missing.push({ gene: g, reason: "absent" });
      else if (!zOf(j)) missing.push({ gene: names[j], reason: "flat" });
      else { used.push(names[j]); cols.push(j); }
    }
    if (cols.length >= need) {
      scored[ti] = true;
      for (const j of cols) { const x = z.get(j); for (let i = 0; i < N; i++) score[i * T + ti] += x[i]; }
      for (let i = 0; i < N; i++) score[i * T + ti] /= cols.length;
    }
    return { name: t.name, used, missing, need };
  });
  const out = ctSummarise(score, scored, types, cl, nCl, PROTEIN_CT.binMin, set,
    (t) => `${t.used.length} of ${t.used.length + t.missing.length} channels measurable here, needs ${t.need}`);
  out.types.forEach((t, ti) => { t.min = types[ti].need; });
  if (flags) {                                 // flags[i]: the qc_flags byte of the run's square i
    const br = new Float64Array(nCl), size = new Float64Array(nCl);
    for (let i = 0; i < N; i++) { size[cl[i]]++; if (flags[i] & 1) br[cl[i]]++; }
    for (const c of out.clusters) {
      const share = size[c.cluster] ? br[c.cluster] / size[c.cluster] : 0;
      c.artefact = Math.round(share * 1000) / 1000;
      if (share > PROTEIN_CT.artefactShare) Object.assign(c, { suggest: ARTEFACT, confidence: "low", artefactRule: true });
    }
  }
  out.modality = "protein";
  out.scale = sec ? "slide" : "run";
  return out;
}

// Top 10 channels per cluster by the difference of the mean z (cluster - rest), lfc > 0 only;
// pct = share of the cluster's squares whose x' is above the channel's median over the whole
// section (sectionumap.protein_markers, on the run's squares). Field names as for genes.
function proteinMarkers(Z, XP, N, K, med, cl, nCl, names) {
  const sz = new Float64Array(nCl * K), ab = new Float64Array(nCl * K), tot = new Float64Array(K), size = new Float64Array(nCl);
  for (let i = 0; i < N; i++) {
    const c = cl[i]; size[c]++;
    for (let j = 0; j < K; j++) { const v = Z[i * K + j]; sz[c * K + j] += v; tot[j] += v; if (XP[i * K + j] > med[j]) ab[c * K + j]++; }
  }
  const r3 = (x) => Math.round(x * 1000) / 1000;
  const out = [];
  for (let c = 0; c < nCl; c++) {
    const rest = Math.max(N - size[c], 1), rows = [];
    for (let j = 0; j < K; j++) {
      const diff = sz[c * K + j] / size[c] - (tot[j] - sz[c * K + j]) / rest;
      if (diff > 0) rows.push([j, diff, ab[c * K + j] / size[c]]);
    }
    rows.sort((a, b) => b[1] - a[1] || a[0] - b[0]);
    out.push(rows.slice(0, 10).map(([j, lfc, pct]) => ({ gene: names[j], gi: j, lfc: r3(lfc), pct: r3(pct) })));
  }
  return out;
}

// The same per cluster, but against the whole slide: the cluster's mean z with the slide's mean and
// sd (sec), lfc > 0 only. What a cluster has that the rest of the slide lacks - a marker high in
// every square drawn (CD45 over an immune circle) leads here and is absent from proteinMarkers.
function proteinMarkersVsSlide(XP, N, K, med, sec, cl, nCl, names) {
  const sz = new Float64Array(nCl * K), ab = new Float64Array(nCl * K), size = new Float64Array(nCl);
  for (let i = 0; i < N; i++) {
    const c = cl[i]; size[c]++;
    for (let j = 0; j < K; j++) {
      const x = XP[i * K + j];
      if (sec.sd[j] > 0) sz[c * K + j] += (x - sec.mean[j]) / sec.sd[j];
      if (x > med[j]) ab[c * K + j]++;
    }
  }
  const r3 = (x) => Math.round(x * 1000) / 1000;
  const out = [];
  for (let c = 0; c < nCl; c++) {
    const rows = [];
    for (let j = 0; j < K; j++) { const z = sz[c * K + j] / Math.max(size[c], 1); if (z > 0) rows.push([j, z, ab[c * K + j] / Math.max(size[c], 1)]); }
    rows.sort((a, b) => b[1] - a[1] || a[0] - b[0]);
    out.push(rows.slice(0, 10).map(([j, lfc, pct]) => ({ gene: names[j], gi: j, lfc: r3(lfc), pct: r3(pct) })));
  }
  return out;
}

async function runProtein(msg) {
  const P = msg.params, rand = mulberry32(P.seed), pr = msg.protein, K = pr.channels.length;
  progress("Reading squares", 0);
  // gather: row, shape and where the raw means of every selected square sit
  const sel = [];
  for (const buf of msg.chunks) {
    const n = new DataView(buf).getUint32(0, true);
    const row = new Uint32Array(buf.slice(4, 4 + 4 * n));
    const q = new Uint16Array(buf.slice(4 + 4 * n, 4 + 4 * n + 2 * n * K));
    for (let i = 0; i < n; i++) {
      const s = msg.rowShape[row[i]];
      if (s !== NO_SHAPE) sel.push({ row: row[i], shape: s, q, at: i * K });
    }
  }
  sel.sort((a, b) => a.row - b.row);           // row order before anything seeded (see run)
  const nShapes = msg.nShapes, found = new Array(nShapes).fill(0), dropped = new Array(nShapes).fill(0);
  sel.forEach((r) => found[r.shape]++);
  // x' with the whole section's window, so a circle reads the slide as the section UMAP did
  // = protgrid.normalise_bg_asinh ("normalize": "bg-asinh-v1"); older exports: the p1..p99.5 window.
  // uOf = the linear version, 0 at background, 1 at p99.5 (the CLR features below start from it).
  const cof = pr.cofactor || 0.2, asinh1 = Math.asinh(1 / cof), bgAsinh = pr.normalize === "bg-asinh-v1";
  const uOf = bgAsinh
    ? (r, j) => { const hi = pr.hi[j]; return hi > 0 ? Math.max(r.q[r.at + j] / pr.scale - pr.bg[j], 0) / hi : 0; }
    : (r, j) => { const lo = pr.p1[j], rg = pr.p99_5[j] - lo; return rg > 0 ? Math.min(1, Math.max(0, (r.q[r.at + j] / pr.scale - lo) / rg)) : 0; };
  const xpOf = bgAsinh ? (r, j) => Math.asinh(uOf(r, j) / cof) / asinh1 : uOf;
  // squares with no marker at x' >= LOW_SIGNAL_X are rims and gaps: left out of the dots and
  // clusters, as the section UMAP leaves them out (protgrid.LOW_SIGNAL_X)
  let sel2 = sel;
  if (P.mode === "bins") {
    sel2 = sel.filter((r) => { for (let j = 0; j < K; j++) if (xpOf(r, j) >= LOW_SIGNAL_X) return true; dropped[r.shape]++; return false; });
  }
  sel.length = 0; for (const r of sel2) sel.push(r);
  let keep = sel.map((_, i) => i), sampledFrom = 0;
  if (P.mode === "bins" && sel.length > P.cap) {
    sampledFrom = sel.length;
    keep = stratifiedSample(sel.map((r) => r.shape), nShapes, P.cap, rand);
  }
  const kept = new Array(nShapes).fill(0);
  keep.forEach((i) => kept[sel[i].shape]++);
  let rows, shape, XP, UL, N;
  if (P.mode === "shapes") {                   // one mean profile per non-empty shape
    const live = kept.map((n, s) => (n ? s : -1)).filter((s) => s >= 0), at = new Map(live.map((s, i) => [s, i]));
    N = live.length; XP = new Float64Array(N * K); UL = new Float64Array(N * K);
    for (const i of keep) { const r = sel[i], a = at.get(r.shape); for (let j = 0; j < K; j++) { XP[a * K + j] += xpOf(r, j); UL[a * K + j] += uOf(r, j); } }
    live.forEach((s, a) => { for (let j = 0; j < K; j++) { XP[a * K + j] /= kept[s]; UL[a * K + j] /= kept[s]; } });
    rows = live.slice(); shape = live.slice();
  } else {
    N = keep.length; XP = new Float64Array(N * K); UL = new Float64Array(N * K); rows = new Array(N); shape = new Array(N);
    keep.forEach((i, a) => { const r = sel[i]; rows[a] = r.row; shape[a] = r.shape; for (let j = 0; j < K; j++) { XP[a * K + j] = xpOf(r, j); UL[a * K + j] = uOf(r, j); } });
  }
  // the layout's input: centred log-ratio of u, log(u + CLR_EPS) minus the square's mean
  // (protgrid.clr_features): which markers are up relative to each other, not how bright the
  // square is overall (on x' the first PC was overall brightness and clusters were its bands)
  // (fewer than CLR_MIN_CHANNELS channels: one bright marker moves the square's mean too much; x' as before)
  const useClr = K >= CLR_MIN_CHANNELS;
  if (!useClr) UL = XP;
  for (let i = 0; useClr && i < N; i++) {
    let m = 0; for (let j = 0; j < K; j++) { const l = Math.log(UL[i * K + j] + CLR_EPS); UL[i * K + j] = l; m += l; }
    m /= K; for (let j = 0; j < K; j++) UL[i * K + j] -= m;
  }
  const info = { timings: {}, found, dropped, kept, sampledFrom, nIn: N, modality: "protein", features: useClr ? "clr" : "x", clrEps: CLR_EPS, lowSignalX: LOW_SIGNAL_X };
  if (P.mode === "bins" && N < P.minPoints) throw new Error(`UMAP needs at least ${P.minPoints} squares; the shapes hold ${N}.`);
  if (P.mode === "shapes" && N < 3) throw new Error(`"One dot per shape" needs at least 3 non-empty shapes; there are ${N}.`);

  // z-score per channel within the selection (population sd), clipped at 10 as the section's
  progress("Scaling", 0.05);
  const Z = new Float64Array(N * K), Zc = new Float64Array(N * K);
  let varying = 0;
  for (let j = 0; j < K; j++) {
    let m = 0; for (let i = 0; i < N; i++) m += XP[i * K + j];
    m /= N;
    let v = 0; for (let i = 0; i < N; i++) v += (XP[i * K + j] - m) ** 2;
    const sd = Math.sqrt(v / N);
    if (sd > 0) varying++;
    for (let i = 0; i < N; i++) Z[i * K + j] = sd > 0 ? (XP[i * K + j] - m) / sd : 0;
    // the layout's columns: z of the CLR features, clipped at 10 as the section's
    let mf = 0; for (let i = 0; i < N; i++) mf += UL[i * K + j];
    mf /= N;
    let vf = 0; for (let i = 0; i < N; i++) vf += (UL[i * K + j] - mf) ** 2;
    const sf = Math.sqrt(vf / N);
    for (let i = 0; i < N; i++) Zc[i * K + j] = sf > 0 && sd > 0 ? Math.min(10, (UL[i * K + j] - mf) / sf) : 0;
  }
  if (varying < 2) throw new Error("Fewer than 2 markers vary inside the shapes: nothing to lay out.");
  // dense PCA through the sparse routine: every entry present, columns already centred
  progress("PCA", 0.15);
  const ptr = new Uint32Array(N + 1), col = new Uint32Array(N * K);
  for (let i = 0; i < N; i++) { ptr[i + 1] = (i + 1) * K; for (let j = 0; j < K; j++) col[i * K + j] = j; }
  const nPcs = Math.max(2, Math.min(P.nPcsProtein, N - 1, K - 1));
  const pca = randomizedPCA(N, K, ptr, col, Zc, new Float64Array(K), nPcs, rand);

  let emb = new Float32Array(N * 2), method = "umap", umap = null;
  if (P.mode === "shapes" && N < P.minShapesForUmap) {
    method = "pca";
    for (let i = 0; i < N; i++) { emb[2 * i] = pca.scores[i][0]; emb[2 * i + 1] = pca.scores[i][1] || 0; }
  } else {
    umap = fitUmap(pca.scores, N, P.nNeighbors, P.minDist, rand, info.timings);
    const Y = umap.embedding;
    for (let i = 0; i < 2 * N; i++) emb[i] = Y[i];
  }
  info.nGenesUsed = K; info.nPcs = nPcs; info.varFrac = pca.varFrac.slice(0, 2); info.method = method;
  info.nNeighbors = umap ? Math.max(2, Math.min(P.nNeighbors, N - 1)) : null; info.init = umap ? "pca" : null;
  const out = { emb, rows: Uint32Array.from(rows), shape: Uint16Array.from(shape), info };

  const resolution = wantedResolution(msg, P);
  if (P.mode === "bins" && umap && resolution) {
    const cluster = await clusterGraph(fuzzyGraphCSR(umap.graph, N), N, resolution, P.seed, info);
    const nCl = info.nClusters;
    // one cluster has no "rest of the selection" to compare with (vs the slide still works)
    info.markers = nCl > 1 ? proteinMarkers(Z, XP, N, K, pr.median, cluster, nCl, pr.channels) : [[]];
    // the slide's mean / sd per channel: from mat/index.json, or worked out by the page for older exports
    const sec = Array.isArray(pr.mean) && Array.isArray(pr.sd) && pr.mean.length === K && pr.sd.length === K ? { mean: pr.mean, sd: pr.sd } : null;
    if (sec) info.markersVsSlide = proteinMarkersVsSlide(XP, N, K, pr.median, sec, cluster, nCl, pr.channels);
    const ms = msg.markerSet;
    if (ms && Array.isArray(ms.types) && ms.types.length) {
      progress("Suggesting cell types", 1);
      const fl = msg.flags ? Uint8Array.from(rows, (r) => msg.flags[r] || 0) : null;
      info.celltypes = cellTypesProtein(XP, N, K, cluster, nCl, ms, pr.channels, fl, sec);
    }
    out.cluster = cluster;
  }
  return out;
}

// ---- dense helpers for the randomized PCA (column-major: cols[j] is Float64Array(N)) ----
function orthonormalize(cols) {
  for (let j = 0; j < cols.length; j++) {
    const v = cols[j];
    for (let pass = 0; pass < 2; pass++)         // twice is enough for numerical orthogonality
      for (let i = 0; i < j; i++) {
        const u = cols[i]; let d = 0;
        for (let r = 0; r < v.length; r++) d += u[r] * v[r];
        for (let r = 0; r < v.length; r++) v[r] -= d * u[r];
      }
    let nrm = 0; for (let r = 0; r < v.length; r++) nrm += v[r] * v[r];
    nrm = Math.sqrt(nrm) || 1;
    for (let r = 0; r < v.length; r++) v[r] /= nrm;
  }
  return cols;
}

// symmetric eigen (cyclic Jacobi), small l x l
function jacobiEigen(A, l) {
  const a = A.map((r) => r.slice()), V = Array.from({ length: l }, (_, i) => { const r = new Float64Array(l); r[i] = 1; return r; });
  for (let sweep = 0; sweep < 100; sweep++) {
    let off = 0;
    for (let p = 0; p < l; p++) for (let q = p + 1; q < l; q++) off += a[p][q] * a[p][q];
    if (off < 1e-20) break;
    for (let p = 0; p < l; p++) for (let q = p + 1; q < l; q++) {
      if (Math.abs(a[p][q]) < 1e-300) continue;
      const th = (a[q][q] - a[p][p]) / (2 * a[p][q]);
      const t = Math.sign(th || 1) / (Math.abs(th) + Math.sqrt(th * th + 1));
      const c = 1 / Math.sqrt(t * t + 1), s = t * c;
      for (let k = 0; k < l; k++) {
        const akp = a[k][p], akq = a[k][q];
        a[k][p] = c * akp - s * akq; a[k][q] = s * akp + c * akq;
      }
      for (let k = 0; k < l; k++) {
        const apk = a[p][k], aqk = a[q][k];
        a[p][k] = c * apk - s * aqk; a[q][k] = s * apk + c * aqk;
      }
      for (let k = 0; k < l; k++) {
        const vkp = V[k][p], vkq = V[k][q];
        V[k][p] = c * vkp - s * vkq; V[k][q] = s * vkp + c * vkq;
      }
    }
  }
  const vals = a.map((r, i) => r[i]);
  const order = vals.map((_, i) => i).sort((x, y) => vals[y] - vals[x]);
  return { vals: order.map((i) => vals[i]), vecs: order.map((i) => V.map((r) => r[i])) };   // vecs[j] = j-th eigenvector
}

// Z = W - 1 c^T, W sparse CSR (ptr, col, val) N x K. Returns top-k PC scores (N x k).
function randomizedPCA(N, K, ptr, col, val, c, k, rand, wantLoadings) {
  const l = Math.min(k + 10, N, K);
  const omega = Array.from({ length: l }, () => { const v = new Float64Array(K); for (let g = 0; g < K; g++) v[g] = gauss(rand); return v; });
  const ZM = (M) => M.map((m) => {                       // Z m, m in R^K
    let cm = 0; for (let g = 0; g < K; g++) cm += c[g] * m[g];
    const out = new Float64Array(N);
    for (let i = 0; i < N; i++) { let s = 0; for (let p = ptr[i]; p < ptr[i + 1]; p++) s += val[p] * m[col[p]]; out[i] = s - cm; }
    return out;
  });
  const ZtY = (Y) => Y.map((y) => {                      // Z^T y, y in R^N
    let sy = 0; for (let i = 0; i < N; i++) sy += y[i];
    const out = new Float64Array(K);
    for (let i = 0; i < N; i++) { const yi = y[i]; if (!yi) continue; for (let p = ptr[i]; p < ptr[i + 1]; p++) out[col[p]] += val[p] * yi; }
    for (let g = 0; g < K; g++) out[g] -= c[g] * sy;
    return out;
  });
  let Q = orthonormalize(ZM(omega));
  for (let it = 0; it < 3; it++) Q = orthonormalize(ZM(orthonormalize(ZtY(Q))));
  const Bt = ZtY(Q);                                     // columns of B^T: l vectors in R^K
  const G = Array.from({ length: l }, (_, a) => { const r = new Float64Array(l); for (let b = 0; b < l; b++) { let s = 0; const x = Bt[a], y = Bt[b]; for (let g = 0; g < K; g++) s += x[g] * y[g]; r[b] = s; } return r; });
  const { vals, vecs } = jacobiEigen(G, l);
  const kk = Math.min(k, l);
  const scores = Array.from({ length: N }, () => new Array(kk));
  // wantLoadings (fast mode): v_j = sgn B^T u_j / s_j, so that a row's score_j = (its row of Z) . v_j,
  // for placing bins that were not in the PCA
  const loadings = wantLoadings ? [] : null;
  for (let j = 0; j < kk; j++) {
    const u = vecs[j];
    // B B^T u = s^2 u, and Z v = Q B v = Q u s  ->  score_j = Q u_j * s_j
    // Fix the sign so the largest-magnitude loading is positive: reruns do not mirror.
    const sc = new Float64Array(N);
    for (let a = 0; a < l; a++) { const q = Q[a], w = u[a]; for (let i = 0; i < N; i++) sc[i] += q[i] * w; }
    let big = 0; for (let i = 0; i < N; i++) if (Math.abs(sc[i]) > Math.abs(big)) big = sc[i];
    const sgn = big < 0 ? -1 : 1, sv = Math.sqrt(Math.max(vals[j], 0));
    for (let i = 0; i < N; i++) scores[i][j] = sgn * sc[i] * sv;
    if (loadings) {
      const v = new Float64Array(K);
      if (sv > 0) { for (let a = 0; a < l; a++) { const bt = Bt[a], w = u[a]; for (let g = 0; g < K; g++) v[g] += bt[g] * w; } for (let g = 0; g < K; g++) v[g] *= sgn / sv; }
      loadings.push(v);
    }
  }
  const total = vals.reduce((a, b) => a + Math.max(0, b), 0);
  return { scores, loadings, varFrac: vals.slice(0, kk).map((v) => Math.max(0, v) / (total || 1)) };
}

// ---- fast mode: place the bins outside the sample ----------------------------------------
// A small random-projection forest over the sample (the same kind of split as umap-js's trees), then for
// each other bin a best-first walk of the sample's own neighbour graph from the leaves it lands in
// (pynndescent's query, epsilon 0.1): its k nearest sample bins in PCA space. Seeded on its own stream.
function rpTree(Y, n, D, leafSize, rand) {
  const hyper = [], offs = [], kids = [], leaves = [];
  const build = (idx) => {
    if (idx.length <= leafSize) { leaves.push(idx); return -leaves.length; }      // leaf l -> -(l + 1)
    const m = idx.length;
    let li = Math.floor(rand() * m), ri = Math.floor(rand() * m);
    if (li === ri) ri = (ri + 1) % m;
    const lo = idx[li] * D, ro = idx[ri] * D, hv = new Float64Array(D);
    let off = 0;
    for (let i = 0; i < D; i++) { hv[i] = Y[lo + i] - Y[ro + i]; off -= (hv[i] * (Y[lo + i] + Y[ro + i])) / 2; }
    const L = [], R = [];
    for (let i = 0; i < m; i++) {
      let mg = off; const o = idx[i] * D;
      for (let d = 0; d < D; d++) mg += hv[d] * Y[o + d];
      if (mg > 0 || (mg === 0 && rand() < 0.5)) L.push(idx[i]); else R.push(idx[i]);
    }
    if (!L.length || !R.length) { leaves.push(idx); return -leaves.length; }      // no split: one leaf
    const me = hyper.length;
    hyper.push(hv); offs.push(off); kids.push([0, 0]);
    kids[me][0] = build(Uint32Array.from(L)); kids[me][1] = build(Uint32Array.from(R));
    return me;
  };
  const all = new Uint32Array(n); for (let i = 0; i < n; i++) all[i] = i;
  const root = build(all);
  return (Q, qo) => {
    let node = root;
    while (node >= 0) {
      let mg = offs[node]; const hv = hyper[node];
      for (let d = 0; d < D; d++) mg += hv[d] * Q[qo + d];
      node = kids[node][mg >= 0 ? 0 : 1];
    }
    return leaves[-node - 1];
  };
}

// k nearest of each query row of Q (nQ x D) among the rows of Y (nS x D); graph: Y's own neighbour
// lists (gk per row). Returns squared distances, ascending per row.
function nearestInSample(Y, nS, Q, nQ, D, graph, gk, k, rand) {
  const trees = []; for (let t = 0; t < 4; t++) trees.push(rpTree(Y, nS, D, 30, rand));
  const outI = new Float64Array(nQ * k), outD = new Float64Array(nQ * k);   // -1: fewer than k found (never in practice)
  const stamp = new Uint32Array(nS);
  const bi = new Float64Array(k), bd = new Float64Array(k);     // the best k: a max-heap
  let ci = new Uint32Array(1024), cd = new Float64Array(1024), nc = 0;  // candidates to expand: a min-heap
  const d2 = (qo, p) => { const o = p * D; let r = 0; for (let i = 0; i < D; i++) { const x = Q[qo + i] - Y[o + i]; r += x * x; } return r; };
  const bestPush = (d, p) => {                                  // replace the root, sift down
    if (d >= bd[0]) return false;
    let i = 0;
    for (;;) {
      const l = 2 * i + 1, r = l + 1;
      let s = i, sd = d;
      if (l < k && bd[l] > sd) { s = l; sd = bd[l]; }
      if (r < k && bd[r] > sd) { s = r; }
      if (s === i) break;
      bd[i] = bd[s]; bi[i] = bi[s]; i = s;
    }
    bd[i] = d; bi[i] = p;
    return true;
  };
  const candPush = (d, p) => {
    if (nc === ci.length) { const a = new Uint32Array(2 * nc), b = new Float64Array(2 * nc); a.set(ci); b.set(cd); ci = a; cd = b; }
    let i = nc++;
    while (i > 0) { const up = (i - 1) >> 1; if (cd[up] <= d) break; cd[i] = cd[up]; ci[i] = ci[up]; i = up; }
    cd[i] = d; ci[i] = p;
  };
  const candPop = () => {
    const top = ci[0], last = --nc, ld = cd[last], lp = ci[last];
    let i = 0;
    for (;;) {
      const l = 2 * i + 1, r = l + 1;
      if (l >= nc) break;
      const c = r < nc && cd[r] < cd[l] ? r : l;
      if (cd[c] >= ld) break;
      cd[i] = cd[c]; ci[i] = ci[c]; i = c;
    }
    cd[i] = ld; ci[i] = lp;
    return top;
  };
  const eps2 = 1.1 * 1.1;
  for (let q = 0, s = 1; q < nQ; q++, s++) {
    const qo = q * D;
    bd.fill(Infinity); bi.fill(-1); nc = 0;
    for (const tree of trees) {
      for (const p of tree(Q, qo)) {
        if (stamp[p] === s) continue;
        stamp[p] = s;
        const d = d2(qo, p);
        bestPush(d, p); candPush(d, p);
      }
    }
    while (nc) {
      const dTop = cd[0];
      if (dTop > bd[0] * eps2) break;
      const p = candPop(), o = p * gk;
      for (let j = 0; j < gk; j++) {
        const u = graph[o + j];
        if (u < 0 || stamp[u] === s) continue;
        stamp[u] = s;
        const d = d2(qo, u);
        if (bestPush(d, u)) candPush(d, u);
      }
    }
    // ascending
    const ord = Array.from({ length: k }, (_, j) => j).sort((a, b) => bd[a] - bd[b] || bi[a] - bi[b]);
    for (let j = 0; j < k; j++) { outI[q * k + j] = bi[ord[j]]; outD[q * k + j] = bd[ord[j]]; }
  }
  return { idx: outI, d2: outD };
}

// weights of a placed bin's k sample neighbours, as umap-js's transform starts them: exp(-d / sigma)
// with sigma found as for the fuzzy graph (rho 0: transform's local connectivity), normalised to sum 1
function placeWeights(nb, nQ, k) {
  const w = new Float64Array(nQ * k), target = Math.log(k) / Math.log(2);
  for (let q = 0; q < nQ; q++) {
    const o = q * k, dist = (j) => Math.sqrt(nb.d2[o + j]);
    let lo = 0, hi = Infinity, mid = 1;
    for (let it = 0; it < 64; it++) {
      let ps = 0; for (let j = 1; j < k; j++) { const d = dist(j); ps += d > 0 ? Math.exp(-(d / mid)) : 1; }
      if (Math.abs(ps - target) < 1e-5) break;
      if (ps > target) { hi = mid; mid = (lo + hi) / 2; } else { lo = mid; if (hi === Infinity) mid *= 2; else mid = (lo + hi) / 2; }
    }
    let tot = 0;
    for (let j = 0; j < k; j++) { const d = dist(j), v = d <= 0 ? 1 : Math.exp(-(d / mid)); w[o + j] = v; tot += v; }
    for (let j = 0; j < k; j++) w[o + j] = tot > 0 ? w[o + j] / tot : 1 / k;
  }
  return w;
}

// ---- the last run's layout, so the same shapes again with other clustering settings only re-cluster ----
let lastLayout = null;
function hashWords(h, buf, off, len) {
  const n4 = len >> 2, u32 = (off & 3) === 0 ? new Uint32Array(buf, off, n4) : new Uint32Array(buf.slice(off, off + 4 * n4));
  for (let i = 0; i < n4; i++) h = Math.imul(h ^ u32[i], 0x01000193) >>> 0;
  const u8 = new Uint8Array(buf, off + 4 * n4, len - 4 * n4);
  for (let i = 0; i < u8.length; i++) h = Math.imul(h ^ u8[i], 0x01000193) >>> 0;
  return h;
}
function layoutKey(msg) {
  const P = { ...msg.params }; delete P.resolution;
  const view = (a) => (a && a.buffer ? hashWords(0x811c9dc5, a.buffer, a.byteOffset, a.byteLength) : 0);
  let sum = 0, xor = 0;                       // the chunks in any order
  for (const b of msg.chunks) { const h = hashWords(0x811c9dc5, b, 0, b.byteLength); sum = (sum + h) >>> 0; xor = (xor ^ h) >>> 0; }
  return [msg.dataKey || "", JSON.stringify(P), msg.nGenes, msg.nShapes, view(msg.rowShape), view(msg.qcKeep), msg.chunks.length, sum, xor].join("|");
}
const clone = (o) => JSON.parse(JSON.stringify(o));

// log1p CP10k, gene choice, scaling, PCA and UMAP of the bins in d: the exact path (and fast mode's sample)
// during(fit, N): called once the graph is built, before the epochs; returns the side work's promise or null
async function exactCore(d, msg, rand, T, wantLoadings, during) {
  const P = msg.params, N = d.rows.length;
  let t = now();
  // log1p CP10k, then gene stats over the selection (zeros included) and in how many bins
  // each gene is counted at all
  progress("Normalising", 0.05);
  const G = msg.nGenes, val = new Float64Array(d.ct.length);
  const sum = new Float64Array(G), sq = new Float64Array(G), det = new Uint32Array(G);
  for (let i = 0; i < N; i++) {
    const f = 1e4 / d.lib[i];
    for (let k = d.ptr[i]; k < d.ptr[i + 1]; k++) {
      const v = Math.log1p(d.ct[k] * f); val[k] = v; sum[d.gi[k]] += v; sq[d.gi[k]] += v * v;
      if (d.ct[k] > 0) det[d.gi[k]]++;
    }
  }
  const mu = new Float64Array(G), vr = new Float64Array(G);
  for (let g = 0; g < G; g++) { mu[g] = sum[g] / N; vr[g] = Math.max(0, sq[g] / N - mu[g] * mu[g]); }

  // Genes for the layout: the section's variable genes (the hvg/ chunks), minus any counted in
  // fewer than max(minDetect, min(minDetectFrac * N, minDetectCap)) of the bins. At a few
  // hundred bins such genes top any dispersion ranking and, scaled to sd 1 below, a handful of
  // counts would steer the PCs; the cap (30 bins, the smallest cluster worth reporting) keeps
  // the markers of a small population in a big selection. No re-ranking inside the selection:
  // on a 30,000-bin test it gave 2 stable clusters that matched the section clusters at ARI
  // 0.02, against 4 at ARI 0.29 with the section's genes; and every selection now lays out on
  // the same genes.
  progress("Choosing genes", 0.1);
  const minDet = Math.max(P.minDetect !== undefined ? P.minDetect : 5,
    Math.min(Math.ceil((P.minDetectFrac !== undefined ? P.minDetectFrac : 0.01) * N), P.minDetectCap !== undefined ? P.minDetectCap : 30));
  const genes = []; for (let g = 0; g < G; g++) if (det[g] >= minDet && vr[g] > 0) genes.push(g);
  const K = genes.length;
  if (K < 3) throw new Error("Too few genes vary inside the shapes to compute a layout.");
  const colOf = new Array(G).fill(-1); genes.forEach((g, j) => (colOf[g] = j));

  // scale to sd 1, clip at 10; Z = W - 1 c^T keeps it sparse
  const ddof = Math.max(1, N - 1);
  const sd = genes.map((g) => Math.sqrt(vr[g] * N / ddof) || 1);
  const c = new Float64Array(K); genes.forEach((g, j) => (c[j] = mu[g] / sd[j]));
  let nnzW = 0;
  for (let k = 0; k < d.gi.length; k++) if (colOf[d.gi[k]] >= 0) nnzW++;
  const wptr = new Uint32Array(N + 1), wcol = new Uint32Array(nnzW), wval = new Float64Array(nnzW);
  let nw = 0;
  for (let i = 0; i < N; i++) {
    for (let k = d.ptr[i]; k < d.ptr[i + 1]; k++) {
      const j = colOf[d.gi[k]]; if (j < 0) continue;
      wcol[nw] = j; wval[nw++] = Math.min(10, (val[k] - mu[d.gi[k]]) / sd[j]) + c[j];
    }
    wptr[i + 1] = nw;
  }
  T.normalise = now() - t; t = now();
  progress("PCA", 0.15);
  const size = sizeRule(N);
  const nPcs = Math.max(2, Math.min(size.pcs, P.nPcs, N - 1, K - 1));
  const pca = randomizedPCA(N, K, wptr, wcol, wval, c, nPcs, rand, wantLoadings);
  T.pca = now() - t;

  let emb = new Float32Array(N * 2), method = "umap", umap = null;
  if (P.mode === "shapes" && N < P.minShapesForUmap) {
    method = "pca";
    for (let i = 0; i < N; i++) { emb[2 * i] = pca.scores[i][0]; emb[2 * i + 1] = pca.scores[i][1] || 0; }
  } else {
    umap = fitUmapStart(pca.scores, N, size.k, P.minDist, rand, T);
    await fitUmapEpochs(umap, T, during ? during(umap, N) : null);
    const Y = umap.embedding;
    for (let i = 0; i < 2 * N; i++) emb[i] = Y[i];
  }
  return { N, val, mu, sd, c, genes, colOf, K, minDet, nPcs, size, pca, umap, emb, method };
}

// gather, the QC filter, the UMI floor, the cap; then the exact path, or in fast mode the sample's
// exact path and every other bin placed by its nearest sample bins. startSweep(csr, n, fast): the Leiden
// sweep on the layout's graph, started before its epochs (run below), or null
async function layout(msg, startSweep) {
  const P = msg.params, rand = mulberry32(P.seed), T = {};
  let t = now();
  progress("Reading bins", 0);
  let d = gather(msg.chunks, msg.rowShape);
  // Chunks arrive in whatever order the network finishes them. Put the bins in row order
  // before anything seeded touches them, or the same shapes and seed give a new layout.
  d = subset(d, d.rows.map((_, i) => i).sort((a, b) => d.rows[a] - d.rows[b]));
  const nShapes = msg.nShapes;
  const found = new Array(nShapes).fill(0), dropped = new Array(nShapes).fill(0);
  // the page's QC filter (qc.js): one byte per row, 0 = removed; absent when nothing is removed
  const qk = msg.qcKeep && msg.qcKeep.length ? msg.qcKeep : null, qcOut = new Array(nShapes).fill(0);
  const ok = [];
  for (let i = 0; i < d.rows.length; i++) {
    found[d.shape[i]]++;
    if (qk && !qk[d.rows[i]]) { qcOut[d.shape[i]]++; continue; }
    if (d.lib[i] >= P.minUmi) ok.push(i); else dropped[d.shape[i]]++;
  }
  let sampledFrom = 0;
  let keep = ok;
  if (P.mode === "bins" && ok.length > P.cap) {
    sampledFrom = ok.length;
    const s = stratifiedSample(ok.map((i) => d.shape[i]), nShapes, P.cap, rand);
    keep = s.map((j) => ok[j]);
  }
  d = subset(d, keep);
  const kept = new Array(nShapes).fill(0);
  d.shape.forEach((s) => kept[s]++);
  const info = { timings: T, found, dropped, kept, sampledFrom, nIn: d.rows.length };
  if (qk) info.qcRemoved = qcOut;

  // one summed profile per shape
  if (P.mode === "shapes") {
    const live = kept.map((n, s) => (n ? s : -1)).filter((s) => s >= 0);
    const acc = live.map(() => new Map()), libs = live.map(() => 0), pos = new Map(live.map((s, i) => [s, i]));
    for (let i = 0; i < d.rows.length; i++) {
      const j = pos.get(d.shape[i]); libs[j] += d.lib[i];
      for (let k = d.ptr[i]; k < d.ptr[i + 1]; k++) acc[j].set(d.gi[k], (acc[j].get(d.gi[k]) || 0) + d.ct[k]);
    }
    const agg = { rows: live.slice(), shape: live.slice(), lib: libs, ptr: [0], gi: [], ct: [] };
    acc.forEach((m) => { [...m.keys()].sort((a, b) => a - b).forEach((g) => { agg.gi.push(g); agg.ct.push(m.get(g)); }); agg.ptr.push(agg.gi.length); });
    d = agg;
  }
  const N = d.rows.length;
  if (P.mode === "bins" && N < P.minPoints) throw new Error(`UMAP needs at least ${P.minPoints} bins with ≥ ${P.minUmi} UMI${qk ? " that pass the QC filter" : ""}; the shapes hold ${N}.`);
  if (P.mode === "shapes" && N < 3) throw new Error(`"One dot per shape" needs at least 3 non-empty shapes; there are ${N}.`);
  T.read = now() - t;

  // fast mode: "fast", or "auto" above fastAbove bins; never when the sample would be all of them
  const nSample = P.fastSample || 10000;
  const fast = P.mode === "bins" && N > nSample && (P.sampling === "fast" || (P.sampling !== "exact" && P.fastAbove > 0 && N > P.fastAbove));
  const sweep = { csr: null, promise: null };
  const during = startSweep ? (fit, n) => { sweep.csr = fuzzyGraphCSR(fit.graph, n); return (sweep.promise = startSweep(sweep.csr, n, fast)); } : null;
  if (!fast) {
    const core = await exactCore(d, msg, rand, T, false, during);
    info.nGenesUsed = core.K; info.minDetect = core.minDet; info.nPcs = core.nPcs;
    info.nNeighbors = core.umap ? Math.max(2, Math.min(core.size.k, N - 1)) : null; info.init = core.umap ? "pca" : null;
    info.varFrac = core.pca.varFrac.slice(0, 2); info.method = core.method;
    return { d, val: core.val, genes: core.genes, colOf: core.colOf, emb: core.emb, info, umap: !!core.umap,
             csr: core.umap ? sweep.csr || fuzzyGraphCSR(core.umap.graph, N) : null, nG: N, fast: null, sweep: sweep.promise };
  }

  // ---- fast mode ----
  const sIdx = stratifiedSample(d.shape, nShapes, nSample, rand), nS = sIdx.length;
  const ds = subset(d, sIdx);
  const core = await exactCore(ds, msg, rand, T, true, during);
  t = now();
  progress("Placing the other bins", 0.97);
  // every bin: log1p CP10k (for the markers), and its row of Z (the sample's scaling) on the sample's PCs
  const D = core.nPcs, Lg = core.pca.loadings, val = new Float64Array(d.ct.length), proj = new Float64Array(N * D);
  const cL = new Float64Array(D);
  for (let j = 0; j < D; j++) { let s = 0; for (let g = 0; g < core.K; g++) s += core.c[g] * Lg[j][g]; cL[j] = s; }
  for (let i = 0; i < N; i++) {
    const f = 1e4 / d.lib[i], o = i * D;
    for (let j = 0; j < D; j++) proj[o + j] = -cL[j];
    for (let k = d.ptr[i]; k < d.ptr[i + 1]; k++) {
      const v = Math.log1p(d.ct[k] * f); val[k] = v;
      const g = core.colOf[d.gi[k]]; if (g < 0) continue;
      const w = Math.min(10, (v - core.mu[d.gi[k]]) / core.sd[g]) + core.c[g];
      for (let j = 0; j < D; j++) proj[o + j] += w * Lg[j][g];
    }
  }
  const isSample = new Uint8Array(N);
  sIdx.forEach((i) => (isSample[i] = 1));
  const Y = new Float64Array(nS * D);
  for (let j = 0; j < nS; j++) Y.set(proj.subarray(sIdx[j] * D, sIdx[j] * D + D), j * D);
  const others = []; for (let i = 0; i < N; i++) if (!isSample[i]) others.push(i);
  const nQ = others.length, Q = new Float64Array(nQ * D);
  others.forEach((i, q) => Q.set(proj.subarray(i * D, i * D + D), q * D));
  const k = FAST.k, knn = core.umap.knn;
  const nb = nearestInSample(Y, nS, Q, nQ, D, knn.idx, knn.k, k, mulberry32((P.seed ^ 0x2545f491) >>> 0));
  const w = placeWeights(nb, nQ, k);
  const emb = new Float32Array(2 * N), E = core.emb;
  for (let j = 0; j < nS; j++) { emb[2 * sIdx[j]] = E[2 * j]; emb[2 * sIdx[j] + 1] = E[2 * j + 1]; }
  others.forEach((i, q) => {
    let x = 0, y = 0;
    for (let j = 0; j < k; j++) { const s = nb.idx[q * k + j], ww = w[q * k + j]; if (s >= 0) { x += ww * E[2 * s]; y += ww * E[2 * s + 1]; } }
    emb[2 * i] = x; emb[2 * i + 1] = y;
  });
  T.place = now() - t;
  info.nGenesUsed = core.K; info.minDetect = core.minDet; info.nPcs = core.nPcs;
  info.nNeighbors = Math.max(2, Math.min(core.size.k, nS - 1)); info.init = "pca";
  info.varFrac = core.pca.varFrac.slice(0, 2); info.method = "umap";
  info.fast = { sample: nS, of: N, k, isSample, placedBy: `the ${k} nearest sample bins in PCA space (weighted mean of their UMAP positions; their clusters by weighted vote)`,
                above: P.fastAbove || null, chosen: P.sampling === "fast" ? "asked" : "auto" };
  return { d, val, genes: core.genes, colOf: core.colOf, emb, info, umap: true, csr: sweep.csr || fuzzyGraphCSR(core.umap.graph, nS), nG: nS,
           fast: { sIdx, others, nb, w, k }, sweep: sweep.promise };
}
const FAST = { k: 15 };

// Leiden (bins mode only, seeded with the run's seed), markers and suggested names; in fast mode the
// sample's clusters, every other bin taking the weighted vote of its nearest sample bins
async function clusterStage(Lay, msg, info) {
  const P = msg.params, resolution = wantedResolution(msg, P);
  if (!(P.mode === "bins" && Lay.umap && resolution)) return null;
  let t = now();
  // the sweep started alongside the layout (run below) is this run's only: a re-cluster of a kept layout sweeps again
  const pre = Lay.sweep; Lay.sweep = null;
  let cluster;
  if (pre) {
    const r = await pre;
    cluster = r.cluster; Object.assign(info, r.info); info.timings.leiden = r.ms; info.timings.leidenAlongside = 1;
  } else cluster = await clusterGraph(Lay.csr, Lay.nG, resolution, P.seed, info, !!Lay.fast);
  const nCl = info.nClusters;
  if (Lay.fast) {
    const F = Lay.fast, N = Lay.d.rows.length, all = new Uint16Array(N), votes = new Float64Array(nCl);
    F.sIdx.forEach((i, j) => (all[i] = cluster[j]));
    F.others.forEach((i, q) => {
      votes.fill(0);
      for (let j = 0; j < F.k; j++) { const s = F.nb.idx[q * F.k + j]; if (s >= 0) votes[cluster[s]] += F.w[q * F.k + j]; }
      let b = 0; for (let c = 1; c < nCl; c++) if (votes[c] > votes[b]) b = c;
      all[i] = b;
    });
    info.fast.sampleClusterSizes = info.clusterSizes;
    const sizes = new Array(nCl).fill(0); all.forEach((c) => sizes[c]++);
    info.clusterSizes = sizes;
    cluster = all;
  }
  if (pre) info.timings.leidenWait = now() - t; else info.timings.leiden = now() - t;
  t = now();
  info.markers = nCl > 1 ? clusterMarkers(Lay.d, Lay.val, Lay.genes, Lay.colOf, cluster, nCl, msg.geneNames || null) : [[]];
  // a suggested name per cluster, when the page found a marker set and the gene names
  const ms = msg.markerSet;
  if (ms && Array.isArray(ms.types) && ms.types.length && msg.geneNames) {
    progress("Suggesting cell types", 1);
    info.celltypes = cellTypes(Lay.d, Lay.val, cluster, nCl, ms, msg.geneNames);
  }
  info.timings.markers = now() - t;
  return cluster;
}

async function run(msg) {
  const key = layoutKey(msg);
  let Lay, info;
  if (lastLayout && lastLayout.key === key) {
    Lay = lastLayout.lay;
    info = clone({ ...lastLayout.info, fast: null });
    if (lastLayout.info.fast) info.fast = { ...clone({ ...lastLayout.info.fast, isSample: null }), isSample: lastLayout.info.fast.isSample.slice() };
    info.timings = { reusedLayout: 1 };
    progress("Same shapes and settings: re-clustering the last layout", 0.9);
  } else {
    lastLayout = null;                          // a new layout: drop the old one first (memory)
    // with clustering on, the sweep runs on the helper workers while this thread lays out (fitUmapEpochs);
    // it fills its own info, merged in clusterStage, so the kept layout's info holds no clusters
    const P = msg.params, resolution = wantedResolution(msg, P);
    const startSweep = P.mode === "bins" && resolution ? (csr, n, fast) => {
      const t = now(), sw = {};
      sideWork.on = true; sideWork.pool = false; sideWork.step = null;
      return clusterGraph(csr, n, resolution, P.seed, sw, fast).then((cluster) => ({ cluster, info: sw, ms: now() - t }));
    } : null;
    try { Lay = await layout(msg, startSweep); }
    finally { sideWork.on = false; }
    if (Lay.sweep && sideWork.step) progress(sideWork.step, sideWork.frac);
    info = Lay.info;
    const keepInfo = { ...info, fast: null, timings: null };
    lastLayout = { key, lay: Lay, info: { ...clone(keepInfo), fast: info.fast ? { ...clone({ ...info.fast, isSample: null }), isSample: info.fast.isSample.slice() } : null } };
  }
  const out = { emb: Lay.emb.slice(), rows: Uint32Array.from(Lay.d.rows), shape: Uint16Array.from(Lay.d.shape), info };
  const cluster = await clusterStage(Lay, msg, info);
  if (cluster) out.cluster = cluster;
  return out;
}

onmessage = async (ev) => {
  try {
    await umapReady;
    const t0 = now();
    const r = ev.data.protein ? await runProtein(ev.data) : await run(ev.data);
    if (r.info && r.info.timings) r.info.timings.total = now() - t0;
    const tr = [r.emb.buffer, r.rows.buffer, r.shape.buffer];
    if (r.cluster) tr.push(r.cluster.buffer);
    postMessage({ type: "done", ...r }, tr);
  } catch (e) {
    postMessage({ type: "error", message: e.message || String(e) });
  }
};
