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
/* global UMAP, leiden */
importScripts("https://cdn.jsdelivr.net/npm/umap-js@1.4.0/lib/umap-js.min.js");
importScripts("leiden.js?v=20261102a");

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
// umap-js keeps the fuzzy simplicial set (scanpy's "connectivities") as umap.graph, a
// SparseMatrix over a Map of {value, row, col}. Symmetrise to CSR: w_ij = (a_ij + a_ji) / 2
// (the fuzzy union is already symmetric, so this only guards), zeros and self loops dropped.
function fuzzyGraphCSR(umap, N) {
  const deg = new Uint32Array(N + 1);
  umap.graph.forEach((v, r, c) => { if (v > 0 && r !== c) { deg[r + 1]++; deg[c + 1]++; } });
  for (let i = 0; i < N; i++) deg[i + 1] += deg[i];
  const fill = deg.slice(0, N), col = new Uint32Array(deg[N]), val = new Float64Array(deg[N]);
  umap.graph.forEach((v, r, c) => {
    if (!(v > 0) || r === c) return;
    col[fill[r]] = c; val[fill[r]++] = v / 2;
    col[fill[c]] = r; val[fill[c]++] = v / 2;
  });
  // merge the duplicate (r, c) pairs inside each row
  const indptr = new Uint32Array(N + 1), at = new Float64Array(N).fill(-1);
  const indices = [], weights = [];
  for (let i = 0; i < N; i++) {
    const start = indices.length;
    for (let p = deg[i]; p < deg[i + 1]; p++) {
      const j = col[p];
      if (at[j] >= start) weights[at[j]] += val[p];
      else { at[j] = indices.length; indices.push(j); weights.push(val[p]); }
    }
    indptr[i + 1] = indices.length;
  }
  return { indptr, indices: Uint32Array.from(indices), weights: Float64Array.from(weights) };
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
function fitUmap(scores, N, k, minDist, rand) {
  const umap = new UMAP.UMAP({ nComponents: 2, nNeighbors: Math.max(2, Math.min(k, N - 1)), minDist, random: rand });
  const nEpochs = umap.initializeFit(scores), E = umap.embedding;
  for (let a = 0; a < 2; a++) {
    let lo = Infinity, hi = -Infinity;
    for (let i = 0; i < N; i++) { const v = scores[i][a] || 0; if (v < lo) lo = v; if (v > hi) hi = v; }
    const s = hi > lo ? 10 / (hi - lo) : 0;
    for (let i = 0; i < N; i++) E[i][a] = ((scores[i][a] || 0) - lo) * s + 1e-4 * gauss(rand);
  }
  for (let e = 0; e < nEpochs; e++) {
    umap.step();
    if (e % 10 === 0) progress("UMAP", 0.2 + 0.8 * e / nEpochs);
  }
  return umap;
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
// Jaccard overlap of the medoid's cluster with any cluster of each other run. "auto" takes
// the finest resolution with stability >= minStable and > 1 cluster, i.e. the most detail
// that reproduces (the coarsest would nearly always win, since 2-3 clusters are trivially
// stable). If none qualifies, the coarsest resolution is shown, flagged, with the closest
// candidate named: on a homogeneous patch the honest answer is one population, not the
// least-bad dozen. The sweep stops after two resolutions in a row below stopBelow: on the
// test lassos stability never climbed back once it had fallen that far, and the finest
// resolutions cost the most. Leiden runs 2 iterations, as scanpy and sectionumap.py do. Same
// graph and seed give the same answer.
const SWEEP = { res: [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2], frac: 0.8, minStable: 0.8, stopBelow: 0.6, iters: 2 };
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

// the graph restricted to the bins in keep (ascending), renumbered 0..keep.length-1
function inducedGraph(g, keep, N) {
  const at = new Int32Array(N).fill(-1); keep.forEach((v, j) => (at[v] = j));
  const indptr = new Uint32Array(keep.length + 1), indices = [], weights = [];
  keep.forEach((v, j) => {
    for (let p = g.indptr[v]; p < g.indptr[v + 1]; p++) { const u = at[g.indices[p]]; if (u >= 0) { indices.push(u); weights.push(g.weights[p]); } }
    indptr[j + 1] = indices.length;
  });
  return { indptr, indices: Uint32Array.from(indices), weights: Float64Array.from(weights) };
}

const median = (xs) => { if (!xs.length) return null; const s = xs.slice().sort((a, b) => a - b), m = s.length >> 1; return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2; };

// clusters at the wanted resolution ("auto" or a number) with their stability; fills info
function clusterGraph(g, N, resolution, seed, info) {
  const list = resolution === "auto" ? SWEEP.res : [resolution], R = sweepRuns(N), r3 = (x) => Math.round(x * 1000) / 1000;
  const rows = [];
  for (const [li, res] of list.entries()) {
    const last2 = rows.slice(-2);
    if (last2.length === 2 && last2.every((r) => r.stability < SWEEP.stopBelow)) break;
    progress(`Clustering, resolution ${res}`, li / list.length);
    const full = [], subs = [];
    for (let s = 0; s < R.seeds; s++) full.push(leiden(N, g.indptr, g.indices, g.weights, { resolution: res, seed: seed + s, maxIter: SWEEP.iters }));
    const pick = mulberry32((seed ^ 0x9e3779b9) >>> 0);
    for (let s = 0; s < R.subs; s++) {
      const keep = []; for (let i = 0; i < N; i++) if (pick() < SWEEP.frac) keep.push(i);
      const sub = inducedGraph(g, keep, N);
      subs.push({ keep, lab: leiden(keep.length, sub.indptr, sub.indices, sub.weights, { resolution: res, seed: seed + s, maxIter: SWEEP.iters }) });
    }
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
  const multi = rows.filter((r) => r.nCl > 1), ok = multi.filter((r) => r.stability >= SWEEP.minStable);
  const best = resolution !== "auto" ? rows[0] : ok.length ? ok[ok.length - 1] : rows[0];
  const closest = multi.length ? multi.reduce((a, b) => (b.stability > a.stability ? b : a)) : null;
  const sizes = new Array(best.nCl).fill(0); best.ref.forEach((c) => sizes[c]++);
  info.nClusters = best.nCl; info.clusterSizes = sizes; info.resolution = best.res;
  info.autoResolution = resolution === "auto"; info.stability = r3(best.stability);
  info.stable = best.nCl > 1 ? best.stability >= SWEEP.minStable : false;
  info.closest = resolution === "auto" && !ok.length && closest
    ? { resolution: closest.res, nClusters: closest.nCl, stability: r3(closest.stability) } : null;
  info.clusterStability = best.perCluster.map((x) => (x === null ? null : r3(x)));
  info.sweep = rows.map((r) => ({ resolution: r.res, nClusters: r.nCl, stability: r3(r.stability) }));
  info.stabilityRuns = { seeds: R.seeds, subsamples: R.subs, fraction: SWEEP.frac, minStable: SWEEP.minStable,
                         stopBelow: SWEEP.stopBelow, leidenIterations: SWEEP.iters, shown: "medoid run" };
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
// squares: z = (x - mean) / sd per channel over the squares (population sd; sd 0 = "flat");
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
function cellTypesProtein(X, N, K, cl, nCl, set, names, flags) {
  const first = new Map(), aliases = set.aliases || {};
  names.forEach((g, j) => { const k = normName(g); if (!first.has(k)) first.set(k, j); });
  const z = new Map();
  const zOf = (j) => {
    if (z.has(j)) return z.get(j);
    let m = 0; for (let i = 0; i < N; i++) m += X[i * K + j];
    m /= N;
    let v = 0; for (let i = 0; i < N; i++) v += (X[i * K + j] - m) ** 2;
    const sd = Math.sqrt(v / N);
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

function runProtein(msg) {
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
  let keep = sel.map((_, i) => i), sampledFrom = 0;
  if (P.mode === "bins" && sel.length > P.cap) {
    sampledFrom = sel.length;
    keep = stratifiedSample(sel.map((r) => r.shape), nShapes, P.cap, rand);
  }
  const kept = new Array(nShapes).fill(0);
  keep.forEach((i) => kept[sel[i].shape]++);
  // x' with the whole section's window, so a circle reads the slide as the section UMAP did
  // = protgrid.normalise_bg_asinh ("normalize": "bg-asinh-v1"); older exports: the p1..p99.5 window
  const asinh1 = Math.asinh(1 / (pr.cofactor || 0.2));
  const xpOf = pr.normalize === "bg-asinh-v1"
    ? (r, j) => {
      const hi = pr.hi[j], v = r.q[r.at + j] / pr.scale;
      return hi > 0 ? Math.asinh(Math.max(v - pr.bg[j], 0) / hi / pr.cofactor) / asinh1 : 0;
    }
    : (r, j) => {
      const lo = pr.p1[j], rg = pr.p99_5[j] - lo, v = r.q[r.at + j] / pr.scale;
      return rg > 0 ? Math.min(1, Math.max(0, (v - lo) / rg)) : 0;
    };
  let rows, shape, XP, N;
  if (P.mode === "shapes") {                   // one mean profile per non-empty shape
    const live = kept.map((n, s) => (n ? s : -1)).filter((s) => s >= 0), at = new Map(live.map((s, i) => [s, i]));
    N = live.length; XP = new Float64Array(N * K);
    for (const i of keep) { const r = sel[i], a = at.get(r.shape); for (let j = 0; j < K; j++) XP[a * K + j] += xpOf(r, j); }
    live.forEach((s, a) => { for (let j = 0; j < K; j++) XP[a * K + j] /= kept[s]; });
    rows = live.slice(); shape = live.slice();
  } else {
    N = keep.length; XP = new Float64Array(N * K); rows = new Array(N); shape = new Array(N);
    keep.forEach((i, a) => { const r = sel[i]; rows[a] = r.row; shape[a] = r.shape; for (let j = 0; j < K; j++) XP[a * K + j] = xpOf(r, j); });
  }
  const info = { found, dropped, kept, sampledFrom, nIn: N, modality: "protein" };
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
    for (let i = 0; i < N; i++) { const zz = sd > 0 ? (XP[i * K + j] - m) / sd : 0; Z[i * K + j] = zz; Zc[i * K + j] = Math.min(10, zz); }
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
    umap = fitUmap(pca.scores, N, P.nNeighbors, P.minDist, rand);
    const Y = umap.getEmbedding();
    for (let i = 0; i < N; i++) { emb[2 * i] = Y[i][0]; emb[2 * i + 1] = Y[i][1]; }
  }
  info.nGenesUsed = K; info.nPcs = nPcs; info.varFrac = pca.varFrac.slice(0, 2); info.method = method;
  info.nNeighbors = umap ? Math.max(2, Math.min(P.nNeighbors, N - 1)) : null; info.init = umap ? "pca" : null;
  const out = { emb, rows: Uint32Array.from(rows), shape: Uint16Array.from(shape), info };

  const resolution = wantedResolution(msg, P);
  if (P.mode === "bins" && umap && resolution) {
    const cluster = clusterGraph(fuzzyGraphCSR(umap, N), N, resolution, P.seed, info);
    const nCl = info.nClusters;
    info.markers = proteinMarkers(Z, XP, N, K, pr.median, cluster, nCl, pr.channels);
    const ms = msg.markerSet;
    if (ms && Array.isArray(ms.types) && ms.types.length) {
      progress("Suggesting cell types", 1);
      const fl = msg.flags ? Uint8Array.from(rows, (r) => msg.flags[r] || 0) : null;
      info.celltypes = cellTypesProtein(XP, N, K, cluster, nCl, ms, pr.channels, fl);
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
function randomizedPCA(N, K, ptr, col, val, c, k, rand) {
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
  for (let j = 0; j < kk; j++) {
    const u = vecs[j];
    // B B^T u = s^2 u, and Z v = Q B v = Q u s  ->  score_j = Q u_j * s_j
    // Fix the sign so the largest-magnitude loading is positive: reruns do not mirror.
    const sc = new Float64Array(N);
    for (let a = 0; a < l; a++) { const q = Q[a], w = u[a]; for (let i = 0; i < N; i++) sc[i] += q[i] * w; }
    let big = 0; for (let i = 0; i < N; i++) if (Math.abs(sc[i]) > Math.abs(big)) big = sc[i];
    const sgn = big < 0 ? -1 : 1;
    for (let i = 0; i < N; i++) scores[i][j] = sgn * sc[i];
  }
  const total = vals.reduce((a, b) => a + Math.max(0, b), 0);
  return { scores, varFrac: vals.slice(0, kk).map((v) => Math.max(0, v) / (total || 1)) };
}

function run(msg) {
  const P = msg.params, rand = mulberry32(P.seed);
  progress("Reading bins", 0);
  let d = gather(msg.chunks, msg.rowShape);
  // Chunks arrive in whatever order the network finishes them. Put the bins in row order
  // before anything seeded touches them, or the same shapes and seed give a new layout.
  d = subset(d, d.rows.map((_, i) => i).sort((a, b) => d.rows[a] - d.rows[b]));
  const nShapes = msg.nShapes;
  const found = new Array(nShapes).fill(0), dropped = new Array(nShapes).fill(0);
  const ok = [];
  for (let i = 0; i < d.rows.length; i++) {
    found[d.shape[i]]++;
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
  const info = { found, dropped, kept, sampledFrom, nIn: d.rows.length };

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
  if (P.mode === "bins" && N < P.minPoints) throw new Error(`UMAP needs at least ${P.minPoints} bins with ≥ ${P.minUmi} UMI; the shapes hold ${N}.`);
  if (P.mode === "shapes" && N < 3) throw new Error(`"One dot per shape" needs at least 3 non-empty shapes; there are ${N}.`);

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
  const wptr = new Uint32Array(N + 1), wcol = [], wval = [];
  for (let i = 0; i < N; i++) {
    for (let k = d.ptr[i]; k < d.ptr[i + 1]; k++) {
      const j = colOf[d.gi[k]]; if (j < 0) continue;
      wcol.push(j); wval.push(Math.min(10, (val[k] - mu[d.gi[k]]) / sd[j]) + c[j]);
    }
    wptr[i + 1] = wcol.length;
  }
  progress("PCA", 0.15);
  const size = sizeRule(N);
  const nPcs = Math.max(2, Math.min(size.pcs, P.nPcs, N - 1, K - 1));
  const pca = randomizedPCA(N, K, wptr, Uint32Array.from(wcol), Float64Array.from(wval), c, nPcs, rand);

  let emb = new Float32Array(N * 2), method = "umap", umap = null;
  if (P.mode === "shapes" && N < P.minShapesForUmap) {
    method = "pca";
    for (let i = 0; i < N; i++) { emb[2 * i] = pca.scores[i][0]; emb[2 * i + 1] = pca.scores[i][1] || 0; }
  } else {
    umap = fitUmap(pca.scores, N, size.k, P.minDist, rand);
    const Y = umap.getEmbedding();
    for (let i = 0; i < N; i++) { emb[2 * i] = Y[i][0]; emb[2 * i + 1] = Y[i][1]; }
  }
  info.nGenesUsed = K; info.minDetect = minDet; info.nPcs = nPcs;
  info.nNeighbors = umap ? Math.max(2, Math.min(size.k, N - 1)) : null; info.init = umap ? "pca" : null;
  info.varFrac = pca.varFrac.slice(0, 2); info.method = method;
  const out = { emb, rows: Uint32Array.from(d.rows), shape: Uint16Array.from(d.shape), info };

  // Leiden on the UMAP's own neighbour graph (bins mode only), seeded with the run's seed
  const resolution = wantedResolution(msg, P);
  if (P.mode === "bins" && umap && resolution) {
    const cluster = clusterGraph(fuzzyGraphCSR(umap, N), N, resolution, P.seed, info);
    const nCl = info.nClusters;
    info.markers = clusterMarkers(d, val, genes, colOf, cluster, nCl, msg.geneNames || null);
    // a suggested name per cluster, when the page found a marker set and the gene names
    const ms = msg.markerSet;
    if (ms && Array.isArray(ms.types) && ms.types.length && msg.geneNames) {
      progress("Suggesting cell types", 1);
      info.celltypes = cellTypes(d, val, cluster, nCl, ms, msg.geneNames);
    }
    out.cluster = cluster;
  }
  return out;
}

onmessage = (ev) => {
  try {
    const r = ev.data.protein ? runProtein(ev.data) : run(ev.data);
    const tr = [r.emb.buffer, r.rows.buffer, r.shape.buffer];
    if (r.cluster) tr.push(r.cluster.buffer);
    postMessage({ type: "done", ...r }, tr);
  } catch (e) {
    postMessage({ type: "error", message: e.message || String(e) });
  }
};
