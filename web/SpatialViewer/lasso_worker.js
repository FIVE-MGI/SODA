// Circle -> UMAP, off the main thread. The viewer hands over the hvg/ chunks under the
// drawn shapes plus a row -> shape map; this returns a 2-D embedding of the selected bins
// (or of one summed profile per shape). Steps, all seeded so a rerun is identical:
//   gather -> drop low-UMI bins -> cap (stratified sample) -> log1p CP10k ->
//   top variable genes within the selection -> scale -> randomized PCA -> umap-js
//   -> Leiden on the UMAP's fuzzy neighbour graph -> marker genes per cluster.
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
  // each bin's type, then per cluster the share of each type and the mean score
  const NONE = -1, cnt = new Float64Array(nCl * T), sum = new Float64Array(nCl * T), size = new Float64Array(nCl);
  for (let i = 0; i < N; i++) {
    const c = cl[i]; size[c]++;
    let best = NONE, bs = -Infinity;
    for (let t = 0; t < T; t++) {
      if (!scored[t]) continue;
      const s = score[i * T + t]; sum[c * T + t] += s;
      if (s > bs) { bs = s; best = t; }
    }
    if (best !== NONE && bs >= CT.binMin) cnt[c * T + best]++;
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
             ...(scored[ti] ? {} : { reason: `only ${t.used.length} of its ${t.used.length + t.missing.length} genes are measurable here (needs ${CT.minGenes})` }) })),
           clusters };
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

  // log1p CP10k, then gene stats over the selection (zeros included)
  progress("Normalising", 0.05);
  const G = msg.nGenes, val = new Float64Array(d.ct.length);
  const sum = new Float64Array(G), sq = new Float64Array(G);
  for (let i = 0; i < N; i++) {
    const f = 1e4 / d.lib[i];
    for (let k = d.ptr[i]; k < d.ptr[i + 1]; k++) {
      const v = Math.log1p(d.ct[k] * f); val[k] = v; sum[d.gi[k]] += v; sq[d.gi[k]] += v * v;
    }
  }
  const mu = new Float64Array(G), vr = new Float64Array(G);
  for (let g = 0; g < G; g++) { mu[g] = sum[g] / N; vr[g] = Math.max(0, sq[g] / N - mu[g] * mu[g]); }

  // variable genes within the selection: dispersion z-scored inside 20 mean bins
  progress("Choosing genes", 0.1);
  const cand = []; for (let g = 0; g < G; g++) if (mu[g] > 0 && vr[g] > 0) cand.push(g);
  let genes = cand;
  if (cand.length > P.nHvg) {
    const lm = cand.map((g) => Math.log1p(mu[g])), disp = cand.map((g) => Math.log(vr[g] / mu[g]));
    const lo = Math.min(...lm), hi = Math.max(...lm), w = (hi - lo) / 20 || 1;
    const bin = lm.map((x) => Math.min(19, Math.floor((x - lo) / w)));
    const z = new Float64Array(cand.length);
    for (let b = 0; b < 20; b++) {
      const ix = []; bin.forEach((v, i) => { if (v === b) ix.push(i); });
      if (!ix.length) continue;
      const m = ix.reduce((a, i) => a + disp[i], 0) / ix.length;
      const s = Math.sqrt(ix.reduce((a, i) => a + (disp[i] - m) ** 2, 0) / ix.length) || 1;
      ix.forEach((i) => (z[i] = (disp[i] - m) / s));
    }
    genes = cand.map((g, i) => [g, z[i]]).sort((a, b) => b[1] - a[1] || a[0] - b[0]).slice(0, P.nHvg).map((x) => x[0]).sort((a, b) => a - b);
  }
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
  const nPcs = Math.max(2, Math.min(P.nPcs, N - 1, K - 1));
  const pca = randomizedPCA(N, K, wptr, Uint32Array.from(wcol), Float64Array.from(wval), c, nPcs, rand);

  let emb = new Float32Array(N * 2), method = "umap", umap = null;
  if (P.mode === "shapes" && N < P.minShapesForUmap) {
    method = "pca";
    for (let i = 0; i < N; i++) { emb[2 * i] = pca.scores[i][0]; emb[2 * i + 1] = pca.scores[i][1] || 0; }
  } else {
    umap = new UMAP.UMAP({ nComponents: 2, nNeighbors: Math.max(2, Math.min(P.nNeighbors, N - 1)), minDist: P.minDist, random: rand });
    const nEpochs = umap.initializeFit(pca.scores);
    for (let e = 0; e < nEpochs; e++) {
      umap.step();
      if (e % 10 === 0) progress("UMAP", 0.2 + 0.8 * e / nEpochs);
    }
    const Y = umap.getEmbedding();
    for (let i = 0; i < N; i++) { emb[2 * i] = Y[i][0]; emb[2 * i + 1] = Y[i][1]; }
  }
  info.nGenesUsed = K; info.nPcs = nPcs; info.varFrac = pca.varFrac.slice(0, 2); info.method = method;
  const out = { emb, rows: Uint32Array.from(d.rows), shape: Uint16Array.from(d.shape), info };

  // Leiden on the UMAP's own neighbour graph (bins mode only), seeded with the run's seed
  const resolution = Number(msg.resolution !== undefined ? msg.resolution : P.resolution !== undefined ? P.resolution : 1.0);
  if (P.mode === "bins" && umap && resolution) {
    progress("Clustering", 1);
    const g = fuzzyGraphCSR(umap, N);
    const cluster = leiden(N, g.indptr, g.indices, g.weights, { resolution, seed: P.seed });
    let nCl = 0; for (let i = 0; i < N; i++) if (cluster[i] + 1 > nCl) nCl = cluster[i] + 1;
    const sizes = new Array(nCl).fill(0); cluster.forEach((c) => sizes[c]++);
    info.nClusters = nCl; info.clusterSizes = sizes; info.resolution = resolution;
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
    const r = run(ev.data);
    const tr = [r.emb.buffer, r.rows.buffer, r.shape.buffer];
    if (r.cluster) tr.push(r.cluster.buffer);
    postMessage({ type: "done", ...r }, tr);
  } catch (e) {
    postMessage({ type: "error", message: e.message || String(e) });
  }
};
