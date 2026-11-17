// umap-js 1.4.0's fit (UMAP.initializeFit + step), rewritten over typed arrays for Draw a region.
// Same algorithm, same arithmetic, same order of operations and the same draws from the same random
// function, so the same input and seed give the same embedding and graph as umap-js, bit for bit
// (checked against umap-js itself: tests/smoke_test.py, the Draw-a-region fixture). What changed is only
// how the numbers are held: flat typed arrays instead of arrays of arrays and a Map keyed
// by "row:col" strings, no transform search graph (Draw a region never calls transform), and in the
// nearest-neighbour search a distance is left unfinished once it is certain to lose (see dist below).
// Plain script: importScripts("umapfit.js") in a worker sets self.SodaUmapFit.
//
//   const f = SodaUmapFit.init(X, N, D, { nNeighbors, random, a, b })   X: Float64Array N*D, row-major
//   f.nEpochs; f.embedding (Float64Array 2N, start written in place as umap-js's embedding is)
//   for (let e = 0; e < f.nEpochs; e++) f.step(e);
//   f.graph = { n, rows, cols, vals }: the fuzzy simplicial set in umap-js's Map order (umap.graph.forEach)
//   f.knn = { k, idx, dist }: the nearest neighbours, each row sorted by distance (umap.knnIndices / knnDistances)
// a and b: umap-js's own fit of the curve for minDist (abParams below asks umap-js for them).
(function (root) {
  "use strict";
  const SMOOTH_K_TOLERANCE = 1e-5, MIN_K_DIST_SCALE = 1e-3;

  // ---- heaps (umap-js heap.ts): one row of `size` slots per point, flat ----------------------------
  function makeHeap(n, size) {
    // idx in a Float64Array: -1 marks an empty slot, as in umap-js (no signed integer arrays in the page's code: deid.py)
    return { size, idx: new Float64Array(n * size).fill(-1), w: new Float64Array(n * size).fill(Infinity), flag: new Uint8Array(n * size) };
  }
  // heap.heapPush: 0 when the weight is not below the row's largest or the index is already there
  function heapPush(H, row, weight, index, flag) {
    const size = H.size, o = row * size, I = H.idx, W = H.w;
    if (weight >= W[o]) return 0;
    for (let i = 0; i < size; i++) if (index === I[o + i]) return 0;
    return unchecked(H, o, weight, index, flag);
  }
  function unchecked(H, o, weight, index, flag) {
    const size = H.size, I = H.idx, W = H.w, F = H.flag;
    if (weight >= W[o]) return 0;
    W[o] = weight; I[o] = index; F[o] = flag;
    let i = 0, iSwap = 0;
    for (;;) {
      const ic1 = 2 * i + 1, ic2 = ic1 + 1;
      if (ic1 >= size) break;
      else if (ic2 >= size) { if (W[o + ic1] > weight) iSwap = ic1; else break; }
      else if (W[o + ic1] >= W[o + ic2]) { if (weight < W[o + ic1]) iSwap = ic1; else break; }
      else { if (weight < W[o + ic2]) iSwap = ic2; else break; }
      W[o + i] = W[o + iSwap]; I[o + i] = I[o + iSwap]; F[o + i] = F[o + iSwap];
      i = iSwap;
    }
    W[o + i] = weight; I[o + i] = index; F[o + i] = flag;
    return 1;
  }
  // heap.deheapSort: each row ascending by weight
  function deheapSort(H, n) {
    const size = H.size, I = H.idx, W = H.w;
    for (let r = 0; r < n; r++) {
      const o = r * size;
      for (let j = 0; j < size - 1; j++) {
        const last = size - j - 1;
        let t = I[o]; I[o] = I[o + last]; I[o + last] = t;
        t = W[o]; W[o] = W[o + last]; W[o + last] = t;
        // siftDown(dist, ind, ceiling = last, 0)
        let elt = 0;
        while (elt * 2 + 1 < last) {
          const l = elt * 2 + 1, rr = l + 1;
          let swap = elt;
          if (W[o + swap] < W[o + l]) swap = l;
          if (rr < last && W[o + swap] < W[o + rr]) swap = rr;
          if (swap === elt) break;
          let u = W[o + elt]; W[o + elt] = W[o + swap]; W[o + swap] = u;
          u = I[o + elt]; I[o + elt] = I[o + swap]; I[o + swap] = u;
          elt = swap;
        }
      }
    }
  }

  // ---- random projection forest (umap-js tree.ts): the leaves only, in umap-js's order ---------------
  function makeLeaves(X, N, D, leafSize, nTrees, random, leaves) {
    const hv = new Float64Array(D);
    const split = (idx) => {
      if (idx.length <= leafSize) { leaves.push(idx); return; }
      const n = idx.length;
      const leftIndex = Math.floor(random() * n);
      let rightIndex = Math.floor(random() * n);
      rightIndex += leftIndex === rightIndex ? 1 : 0;
      rightIndex = rightIndex % n;
      const lo = idx[leftIndex] * D, ro = idx[rightIndex] * D;
      let offset = 0;
      for (let i = 0; i < D; i++) {
        hv[i] = X[lo + i] - X[ro + i];
        offset -= (hv[i] * (X[lo + i] + X[ro + i])) / 2.0;
      }
      const side = new Uint8Array(n);
      let nLeft = 0, nRight = 0;
      for (let i = 0; i < n; i++) {
        let margin = offset;
        const xo = idx[i] * D;
        for (let d = 0; d < D; d++) margin += hv[d] * X[xo + d];
        if (margin === 0) { side[i] = Math.floor(random() * 2); if (side[i] === 0) nLeft++; else nRight++; }
        else if (margin > 0) { side[i] = 0; nLeft++; }
        else { side[i] = 1; nRight++; }
      }
      const L = new Uint32Array(nLeft), R = new Uint32Array(nRight);
      nLeft = 0; nRight = 0;
      for (let i = 0; i < n; i++) { if (side[i] === 0) L[nLeft++] = idx[i]; else R[nRight++] = idx[i]; }
      split(L); split(R);
    };
    for (let t = 0; t < nTrees; t++) { const all = new Uint32Array(N); for (let i = 0; i < N; i++) all[i] = i; split(all); }
  }

  // ---- nearest neighbours: umap-js nn_descent.ts with its defaults ---------------------------------
  // Euclidean as umap-js's: the sum of (x - y)^2 over the dimensions, in order, then sqrt. (x - y) * (x - y)
  // is the same double as Math.pow(x - y, 2): checked on every engine the tests run (powSquareSame).
  // distBound(p, q, T): the distance when it is below T; otherwise any number >= T. Every distance here
  // only ever feeds two heap pushes that both drop it when it is >= both rows' largest kept distance, so
  // once the partial sum alone makes the final distance >= T the rest of the sum cannot change any
  // outcome and is skipped. The test is safe against rounding: the partial sum only grows (adding
  // squares), sqrt rounds monotonically, and T*T is given a relative margin of 1e-9 (far above 2^-52).
  function nnDescent(X, N, D, k, leaves, random, nIters) {
    const maxCandidates = 50, delta = 0.001, rho = 0.5;
    const dist = (p, q) => {
      const a = p * D, b = q * D;
      let r = 0;
      for (let i = 0; i < D; i++) { const x = X[a + i] - X[b + i]; r += x * x; }
      return Math.sqrt(r);
    };
    const distBound = (p, q, T) => {
      const a = p * D, b = q * D, T2 = T * T * (1 + 1e-9);
      let r = 0;
      for (let i = 0; i < D; i++) {
        const x = X[a + i] - X[b + i]; r += x * x;
        if (r > T2) return T;
      }
      return Math.sqrt(r);
    };
    const G = makeHeap(N, k), W = G.w;
    // a pair goes into both rows' heaps (umap-js pushes p then q, the same here)
    const both = (p, q) => {
      if (p === q) { const d = dist(p, q); return heapPush(G, p, d, q, 1) + heapPush(G, q, d, p, 1); }
      const wp = W[p * k], wq = W[q * k], T = wp > wq ? wp : wq;
      const d = distBound(p, q, T);
      if (d >= T) return 0;
      return heapPush(G, p, d, q, 1) + heapPush(G, q, d, p, 1);
    };
    // heap.rejectionSample: nNeighbors distinct random points per row
    const sample = new Uint32Array(k);
    for (let i = 0; i < N; i++) {
      for (let s = 0; s < k; s++) {
        let j = 0;
        for (;;) {
          j = Math.floor(random() * N);
          let broken = false;
          for (let t = 0; t < s; t++) if (j === sample[t]) { broken = true; break; }
          if (!broken) break;
        }
        sample[s] = j;
      }
      for (let s = 0; s < k; s++) both(i, sample[s]);
    }
    for (const leaf of leaves) {
      for (let i = 0; i < leaf.length; i++) for (let j = i + 1; j < leaf.length; j++) both(leaf[i], leaf[j]);
    }
    const C = makeHeap(N, maxCandidates), CI = C.idx, CF = C.flag;
    for (let n = 0; n < nIters; n++) {
      // heap.buildCandidates
      CI.fill(-1); C.w.fill(Infinity); CF.fill(0);
      for (let i = 0; i < N; i++) {
        for (let j = 0; j < k; j++) {
          const o = i * k + j;
          if (G.idx[o] < 0) continue;
          const idx = G.idx[o], isn = G.flag[o], d = random();
          heapPush(C, i, d, idx, isn);
          heapPush(C, idx, d, i, isn);
          G.flag[o] = 0;
        }
      }
      let c = 0;
      for (let i = 0; i < N; i++) {
        const o = i * maxCandidates;
        for (let j = 0; j < maxCandidates; j++) {
          const p = CI[o + j];
          if (p < 0 || random() < rho) continue;
          const cj = CF[o + j];
          for (let kk = 0; kk < maxCandidates; kk++) {
            const q = CI[o + kk];
            if (q < 0 || (!cj && !CF[o + kk])) continue;
            c += both(p, q);
          }
        }
      }
      if (c <= delta * k * N) break;
    }
    deheapSort(G, N);
    return { k, idx: G.idx, dist: G.w };
  }

  // ---- the fuzzy simplicial set (UMAP.fuzzySimplicialSet with setOpMixRatio 1, localConnectivity 1) --
  function smoothKnnDist(knn, N) {
    const k = knn.k, Dd = knn.dist, target = (Math.log(k) / Math.log(2)) * 1.0;
    const sigmas = new Float64Array(N), rhos = new Float64Array(N);
    const rowMean = (i) => { const o = i * k; let s = Dd[o]; for (let j = 1; j < k; j++) s += Dd[o + j]; return s / k; };
    let meanAll = null;                       // utils.mean(distances.map(utils.mean)), when a row needs it
    for (let i = 0; i < N; i++) {
      const o = i * k;
      let lo = 0.0, hi = Infinity, mid = 1.0;
      // the first non-zero distance (rows are sorted; localConnectivity 1: index 1, no interpolation)
      let first = -1;
      for (let j = 0; j < k; j++) if (Dd[o + j] > 0.0) { first = j; break; }
      if (first >= 0) rhos[i] = Dd[o + first];
      for (let n = 0; n < 64; n++) {
        let psum = 0.0;
        for (let j = 1; j < k; j++) {
          const d = Dd[o + j] - rhos[i];
          if (d > 0) psum += Math.exp(-(d / mid)); else psum += 1.0;
        }
        if (Math.abs(psum - target) < SMOOTH_K_TOLERANCE) break;
        if (psum > target) { hi = mid; mid = (lo + hi) / 2.0; }
        else { lo = mid; if (hi === Infinity) mid *= 2; else mid = (lo + hi) / 2.0; }
      }
      sigmas[i] = mid;
      if (rhos[i] > 0.0) {
        const m = rowMean(i);
        if (sigmas[i] < MIN_K_DIST_SCALE * m) sigmas[i] = MIN_K_DIST_SCALE * m;
      } else {
        if (meanAll === null) { let s = rowMean(0); for (let r = 1; r < N; r++) s += rowMean(r); meanAll = s / N; }
        if (sigmas[i] < MIN_K_DIST_SCALE * meanAll) sigmas[i] = MIN_K_DIST_SCALE * meanAll;
      }
    }
    return { sigmas, rhos };
  }

  // The union a + b - a*b of the membership strengths and their transpose, as umap-js builds it through
  // SparseMatrix (a Map in insertion order): first every (i, knn[i][j]) in row order (a neighbour of -1
  // writes (0, 0, 0), as umap-js's zero-filled arrays do), then every transposed (j, i) not already there,
  // in the same order. Entries with value 0 stay (umap-js keeps them; the users skip them).
  function fuzzyGraph(knn, N) {
    const k = knn.k, KI = knn.idx, KD = knn.dist, { sigmas, rhos } = smoothKnnDist(knn, N);
    const m = N * k, val = new Float64Array(m);
    for (let i = 0; i < N; i++) {
      for (let j = 0; j < k; j++) {
        const o = i * k + j, c = KI[o];
        if (c === -1) continue;
        if (c === i) val[o] = 0.0;
        else if (KD[o] - rhos[i] <= 0.0) val[o] = 1.0;
        else val[o] = Math.exp(-((KD[o] - rhos[i]) / sigmas[i]));
      }
    }
    // S(r, c): the strength of c in r's list (0 when absent); (0, 0) also gathers the -1 slots
    const S = (r, c) => {
      const o = r * k;
      for (let j = 0; j < k; j++) if (KI[o + j] === c) return val[o + j];
      return 0;
    };
    const has = (r, c) => {
      const o = r * k;
      for (let j = 0; j < k; j++) if (KI[o + j] === c) return true;
      return r === 0 && c === 0 && hasMinus1;
    };
    let hasMinus1 = false;
    for (let o = 0; o < m; o++) if (KI[o] === -1) { hasMinus1 = true; break; }
    const rows = new Uint32Array(2 * m), cols = new Uint32Array(2 * m), vals = new Float64Array(2 * m);
    let e = 0, seen00 = false;
    const put = (r, c) => {
      const a = S(r, c), b = S(c, r);
      rows[e] = r; cols[e] = c; vals[e] = (a + b) - a * b; e++;
    };
    for (let i = 0; i < N; i++) {
      for (let j = 0; j < k; j++) {
        let r = i, c = KI[i * k + j];
        if (c === -1) { r = 0; c = 0; }
        if (r === 0 && c === 0) { if (seen00) continue; seen00 = true; }
        put(r, c);
      }
    }
    const nS = e;
    for (let t = 0; t < nS; t++) {
      const r = cols[t], c = rows[t];
      if (has(r, c)) continue;
      put(r, c);
    }
    return { n: N, nnz: e, rows: rows.subarray(0, e), cols: cols.subarray(0, e), vals: vals.subarray(0, e),
             forEach(fn) { for (let t = 0; t < this.nnz; t++) fn(this.vals[t], this.rows[t], this.cols[t]); } };
  }

  // ---- UMAP.getNEpochs, initializeSimplicialSetEmbedding, makeEpochsPerSample -----------------------
  const nEpochsFor = (N) => (N <= 2500 ? 500 : N <= 5000 ? 400 : N <= 7500 ? 300 : 200);

  function init(X, N, D, opt) {
    const nNeighbors = opt.nNeighbors, random = opt.random;
    if (N <= nNeighbors) throw new Error(`Not enough data points (${N}) to create nNeighbors: ${nNeighbors}.  Add more data points or adjust the configuration.`);
    const timing = opt.timing || {}, now = () => (typeof performance !== "undefined" ? performance.now() : Date.now());
    let t = now();
    // UMAP.nearestNeighbors
    const round = (n) => (n === 0.5 ? 0 : Math.round(n));
    const nTrees = 5 + Math.floor(round(Math.pow(N, 0.5) / 20.0));
    const nIters = Math.max(5, Math.floor(Math.round(Math.log(N) / Math.log(2))));
    const leaves = [];
    makeLeaves(X, N, D, Math.max(10, nNeighbors), nTrees, random, leaves);
    const knn = nnDescent(X, N, D, nNeighbors, leaves, random, nIters);
    timing.knn = now() - t; t = now();
    const graph = fuzzyGraph(knn, N);
    timing.graph = now() - t; t = now();
    const nEpochs = nEpochsFor(N);
    let graphMax = 0;
    for (let e = 0; e < graph.nnz; e++) if (graphMax < graph.vals[e]) graphMax = graph.vals[e];
    // umap-js's random start (overwritten by the caller's PCA start, but its draws count)
    const E = new Float64Array(2 * N);
    for (let i = 0; i < 2 * N; i++) E[i] = random() * 20 + -10;
    // the edges kept for the layout, ordered by (row, col) as SparseMatrix.getAll; tail = row, head = col
    const cut = graphMax / nEpochs, cnt = new Uint32Array(N + 1);
    let nE = 0;
    for (let e = 0; e < graph.nnz; e++) { const v = graph.vals[e]; if (v < cut ? 0 : v) { cnt[graph.rows[e] + 1]++; nE++; } }
    for (let i = 0; i < N; i++) cnt[i + 1] += cnt[i];
    const fill = cnt.slice(0, N), ecol = new Uint32Array(nE), ew = new Float64Array(nE);
    for (let e = 0; e < graph.nnz; e++) {
      const v = graph.vals[e];
      if (!(v < cut ? 0 : v)) continue;
      const p = fill[graph.rows[e]]++; ecol[p] = graph.cols[e]; ew[p] = v;
    }
    const head = new Uint32Array(nE), tail = new Uint32Array(nE), weights = new Float64Array(nE);
    for (let i = 0; i < N; i++) {
      const a = cnt[i], b = cnt[i + 1];
      const ord = Array.from({ length: b - a }, (_, x) => a + x).sort((x, y) => ecol[x] - ecol[y]);
      for (let x = 0; x < ord.length; x++) { tail[a + x] = i; head[a + x] = ecol[ord[x]]; weights[a + x] = ew[ord[x]]; }
    }
    let max = 0;
    for (let e = 0; e < nE; e++) max = weights[e] > max ? weights[e] : max;
    const eps = new Float64Array(nE).fill(-1.0);
    for (let e = 0; e < nE; e++) { const ns = (weights[e] / max) * nEpochs; if (ns > 0) eps[e] = nEpochs / ns; }
    // prepareForOptimizationLoop (negativeSampleRate 5, learningRate 1, repulsionStrength 1)
    const epns = new Float64Array(nE), eonn = new Float64Array(nE), eon = Float64Array.from(eps);
    for (let e = 0; e < nE; e++) { epns[e] = eps[e] / 5; eonn[e] = epns[e]; }
    timing.layoutPrep = now() - t;
    const a = opt.a, b = opt.b, gamma = 1.0, initialAlpha = 1.0, clipValue = 4.0, nV = N;
    let alpha = initialAlpha;
    const clip = (x) => (x > clipValue ? clipValue : x < -clipValue ? -clipValue : x);
    // UMAP.optimizeLayoutStep(n), the embedding held flat: point j at E[2j], E[2j + 1]
    function step(n) {
      for (let i = 0; i < nE; i++) {
        if (eon[i] > n) continue;
        const j2 = 2 * head[i], k2 = 2 * tail[i];
        let x = E[j2] - E[k2], y = E[j2 + 1] - E[k2 + 1];
        const distSquared = 0.0 + x * x + y * y;
        let gradCoeff = 0;
        if (distSquared > 0) {
          gradCoeff = -2.0 * a * b * Math.pow(distSquared, b - 1.0);
          gradCoeff /= a * Math.pow(distSquared, b) + 1.0;
        }
        let gradD = clip(gradCoeff * x);
        E[j2] += gradD * alpha; E[k2] += -gradD * alpha;
        gradD = clip(gradCoeff * y);
        E[j2 + 1] += gradD * alpha; E[k2 + 1] += -gradD * alpha;
        eon[i] += eps[i];
        const nNeg = Math.floor((n - eonn[i]) / epns[i]);
        for (let p = 0; p < nNeg; p++) {
          const k1 = Math.floor(random() * nV), o2 = 2 * k1;
          x = E[j2] - E[o2]; y = E[j2 + 1] - E[o2 + 1];
          const ds = 0.0 + x * x + y * y;
          let gc = 0.0;
          if (ds > 0.0) {
            gc = 2.0 * gamma * b;
            gc /= (0.001 + ds) * (a * Math.pow(ds, b) + 1);
          } else if (j2 === o2) continue;
          if (gc > 0.0) {
            E[j2] += clip(gc * x) * alpha;
            E[j2 + 1] += clip(gc * (E[j2 + 1] - E[o2 + 1])) * alpha;
          } else { E[j2] += 4.0 * alpha; E[j2 + 1] += 4.0 * alpha; }
        }
        eonn[i] += nNeg * epns[i];
      }
      alpha = initialAlpha * (1.0 - n / nEpochs);
    }
    return { nEpochs, embedding: E, graph, knn, step, nEdges: nE };
  }

  // a and b of UMAP's curve for this minDist, from umap-js itself (its findABParams, a Levenberg-Marquardt
  // fit, is not exported: a UMAP with an empty graph runs it in initializeOptimization)
  const abMemo = new Map();
  function abParams(UMAPlib, minDist) {
    if (!abMemo.has(minDist)) {
      const u = new UMAPlib.UMAP({ minDist, nEpochs: 1 });
      u.graph = { nRows: 0, nCols: 0 };
      u.initializeOptimization();
      abMemo.set(minDist, [u.optimizationState.a, u.optimizationState.b]);
    }
    return abMemo.get(minDist);
  }

  // (x - y) * (x - y) and Math.pow(x - y, 2) give the same double on this engine (umap-js uses Math.pow)
  function powSquareSame(rand, n) {
    for (let i = 0; i < n; i++) {
      const x = (rand() - 0.5) * Math.pow(10, Math.floor(rand() * 12) - 6);
      if (x * x !== Math.pow(x, 2)) return false;
    }
    return true;
  }

  root.SodaUmapFit = { init, abParams, powSquareSame, nEpochsFor, _test: { makeLeaves, nnDescent, fuzzyGraph, smoothKnnDist } };
})(typeof self !== "undefined" ? self : globalThis);
