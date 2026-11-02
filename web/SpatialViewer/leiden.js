// Leiden community detection (Traag, Waltman & van Eck 2019, "From Louvain to Leiden"),
// no dependencies. Plain script: importScripts("leiden.js") in a worker sets self.leiden,
// require("./leiden.js") in Node returns { leiden, modularity }.
//
//   leiden(n, indptr, indices, weights, { resolution = 1, seed = 42, maxIter = -1, theta = 0.01, stats })
//     -> Uint16Array labels, cluster 0 the largest
//
// The graph is undirected and weighted, given as a symmetric CSR (each edge stored in both
// rows). Self loops and non-positive weights are ignored. Objective: RB-configuration
// modularity  sum_ij (A_ij - resolution * k_i k_j / 2m) [c_i == c_j],  the same one
// leidenalg's RBConfigurationVertexPartition and scanpy's tl.leiden optimise.
// One iteration = fast local moving -> refinement inside each community (singletons merged
// at random among the well-connected options, so every cluster is connected) ->
// aggregate the refined clusters, start from the unrefined partition, repeat until every
// community is one node. Iterations repeat until one changes nothing (maxIter < 0) or
// maxIter is reached. All randomness comes from a seeded mulberry32: same input and seed
// give the same labels.
(function (root) {
  "use strict";

  function mulberry32(seed) {
    let s = seed >>> 0;
    return () => {
      s = (s + 0x6d2b79f5) | 0;
      let t = Math.imul(s ^ (s >>> 15), 1 | s);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }

  function permutation(n, rand) {             // seeded Fisher-Yates
    const p = new Uint32Array(n);
    for (let i = 0; i < n; i++) p[i] = i;
    for (let i = n - 1; i > 0; i--) {
      const j = Math.floor(rand() * (i + 1));
      const t = p[i]; p[i] = p[j]; p[j] = t;
    }
    return p;
  }

  // Input CSR -> { n, ptr, adj, w, k } with self loops / zero weights dropped; k = strength.
  function baseGraph(n, indptr, indices, weights) {
    const ptr = new Uint32Array(n + 1);
    for (let v = 0; v < n; v++) {
      let c = 0;
      for (let p = indptr[v]; p < indptr[v + 1]; p++) if (indices[p] !== v && weights[p] > 0) c++;
      ptr[v + 1] = ptr[v] + c;
    }
    const adj = new Uint32Array(ptr[n]), w = new Float64Array(ptr[n]), k = new Float64Array(n);
    for (let v = 0, o = 0; v < n; v++) {
      for (let p = indptr[v]; p < indptr[v + 1]; p++) {
        const u = indices[p], x = weights[p];
        if (u === v || !(x > 0)) continue;
        adj[o] = u; w[o] = x; k[v] += x; o++;
      }
    }
    return { n, ptr, adj, w, k };
  }

  // Shared scratch space, sized once for the base graph (aggregates are never larger).
  function scratch(n) {
    return { nbW: new Float64Array(n), seen: new Uint8Array(n), touched: new Uint32Array(n),
             cand: new Uint32Array(n + 1), prob: new Float64Array(n + 1) };
  }

  // Fast local moving. comm[] holds ids < g.n; K[] the community strengths (updated in place).
  // Returns whether any node changed community.
  function moveNodes(g, comm, K, r, rand, sc) {
    const { n, ptr, adj, w, k } = g, { nbW, seen, touched } = sc;
    const size = new Uint32Array(n);
    for (let v = 0; v < n; v++) size[comm[v]]++;
    const empty = new Uint32Array(n); let nEmpty = 0;
    for (let c = n - 1; c >= 0; c--) if (!size[c]) empty[nEmpty++] = c;
    const queue = permutation(n, rand), inQ = new Uint8Array(n).fill(1);
    let head = 0, left = n, moved = false;
    while (left > 0) {
      const v = queue[head]; head = head + 1 === n ? 0 : head + 1; left--; inQ[v] = 0;
      const cv = comm[v], kv = k[v];
      let nt = 0;
      for (let p = ptr[v]; p < ptr[v + 1]; p++) {
        const c = comm[adj[p]];
        if (!seen[c]) { seen[c] = 1; nbW[c] = 0; touched[nt++] = c; }
        nbW[c] += w[p];
      }
      K[cv] -= kv;                                 // take v out, then put it where it gains most
      let best = cv, bestGain = (seen[cv] ? nbW[cv] : 0) - r * kv * K[cv];
      for (let t = 0; t < nt; t++) {
        const c = touched[t];
        if (c === cv) continue;
        const gain = nbW[c] - r * kv * K[c];
        if (gain > bestGain) { bestGain = gain; best = c; }
      }
      if (bestGain < 0 && size[cv] > 1) { best = empty[--nEmpty]; bestGain = 0; }   // alone is better
      K[best] += kv;
      for (let t = 0; t < nt; t++) seen[touched[t]] = 0;
      if (best === cv) continue;
      comm[v] = best; moved = true;
      if (--size[cv] === 0) empty[nEmpty++] = cv;
      size[best]++;
      for (let p = ptr[v]; p < ptr[v + 1]; p++) {   // neighbours outside the new community get revisited
        const u = adj[p];
        if (!inQ[u] && comm[u] !== best) {
          inQ[u] = 1; queue[(head + left) % n] = u; left++;
        }
      }
    }
    return moved;
  }

  // Refinement: inside each community of comm, start from singletons and merge each
  // still-singleton, well-connected node into a well-connected refined cluster of the same
  // community, chosen at random with probability ~ exp(gain / theta) among gains >= 0.
  function refine(g, comm, K, r, theta, rand, sc) {
    const { n, ptr, adj, w, k } = g, { nbW, seen, touched, cand, prob } = sc;
    const ref = new Uint32Array(n), refK = Float64Array.from(k), refSize = new Uint32Array(n).fill(1);
    const ext = new Float64Array(n);               // weight from v to the rest of its community
    for (let v = 0; v < n; v++) {
      ref[v] = v;
      for (let p = ptr[v]; p < ptr[v + 1]; p++) if (comm[adj[p]] === comm[v]) ext[v] += w[p];
    }
    const refExt = Float64Array.from(ext);         // weight from refined cluster to rest of its community
    const order = permutation(n, rand);
    for (let o = 0; o < n; o++) {
      const v = order[o], cv = ref[v];
      if (refSize[cv] !== 1) continue;             // only nodes still on their own
      const S = comm[v], KS = K[S], kv = k[v];
      if (ext[v] < r * kv * (KS - kv)) continue;   // v not well connected to its community
      refK[cv] = 0; refSize[cv] = 0;
      let nt = 0;
      for (let p = ptr[v]; p < ptr[v + 1]; p++) {
        const u = adj[p];
        if (comm[u] !== S) continue;
        const c = ref[u];
        if (!seen[c]) { seen[c] = 1; nbW[c] = 0; touched[nt++] = c; }
        nbW[c] += w[p];
      }
      let nc = 0, maxGain = 0;
      cand[nc] = cv; prob[nc++] = 0;               // staying alone: gain 0
      for (let t = 0; t < nt; t++) {
        const c = touched[t];
        if (refExt[c] < r * refK[c] * (KS - refK[c])) continue;   // cluster not well connected
        const gain = nbW[c] - r * kv * refK[c];
        if (gain < 0) continue;
        cand[nc] = c; prob[nc++] = gain;
        if (gain > maxGain) maxGain = gain;
      }
      let best = cv;
      if (nc > 1) {
        let tot = 0;
        for (let i = 0; i < nc; i++) { prob[i] = Math.exp((prob[i] - maxGain) / theta); tot += prob[i]; }
        let x = rand() * tot, i = 0;
        while (i < nc - 1 && x >= prob[i]) { x -= prob[i]; i++; }
        best = cand[i];
      }
      ref[v] = best; refK[best] += kv; refSize[best]++;
      if (best !== cv) refExt[best] += ext[v] - 2 * nbW[best];
      for (let t = 0; t < nt; t++) seen[touched[t]] = 0;
    }
    return ref;
  }

  // Collapse each cluster of part into one node (internal weight dropped: moves never need it;
  // the strength k keeps it). Returns the new graph and node -> new node.
  function aggregate(g, part) {
    const { n, ptr, adj, w, k } = g;
    const NONE = 0xffffffff, id = new Uint32Array(n).fill(NONE), map = new Uint32Array(n);
    let m = 0;
    for (let v = 0; v < n; v++) { const c = part[v]; if (id[c] === NONE) id[c] = m++; map[v] = id[c]; }
    const start = new Uint32Array(m + 1);          // members of each new node, by counting sort
    for (let v = 0; v < n; v++) start[map[v] + 1]++;
    for (let c = 0; c < m; c++) start[c + 1] += start[c];
    const fill = start.slice(0, m), mem = new Uint32Array(n);
    for (let v = 0; v < n; v++) mem[fill[map[v]]++] = v;
    const nk = new Float64Array(m), nptr = new Uint32Array(m + 1);
    const acc = new Float64Array(m), seen = new Uint8Array(m), touched = new Uint32Array(m);
    const nadj = new Uint32Array(ptr[n]), nw = new Float64Array(ptr[n]);   // never more edges than g
    let o = 0;
    for (let c = 0; c < m; c++) {
      let nt = 0;
      for (let q = start[c]; q < start[c + 1]; q++) {
        const v = mem[q]; nk[c] += k[v];
        for (let p = ptr[v]; p < ptr[v + 1]; p++) {
          const d = map[adj[p]];
          if (d === c) continue;
          if (!seen[d]) { seen[d] = 1; acc[d] = 0; touched[nt++] = d; }
          acc[d] += w[p];
        }
      }
      for (let t = 0; t < nt; t++) { const d = touched[t]; nadj[o] = d; nw[o++] = acc[d]; seen[d] = 0; }
      nptr[c + 1] = o;
    }
    return { g: { n: m, ptr: nptr, adj: nadj.slice(0, o), w: nw.slice(0, o), k: nk }, map };
  }

  // ids -> 0..c-1 in order of first appearance; returns the count
  function compact(a) {
    const id = new Map();
    for (let i = 0; i < a.length; i++) {
      let x = id.get(a[i]);
      if (x === undefined) { x = id.size; id.set(a[i], x); }
      a[i] = x;
    }
    return id.size;
  }

  // One Leiden iteration starting from labels (changed in place). Returns whether any node moved.
  function iterate(base, labels, r, theta, rand, sc, stats) {
    let g = base, P = labels.slice();
    compact(P);
    const nodeOf = new Uint32Array(base.n);
    for (let i = 0; i < base.n; i++) nodeOf[i] = i;
    let changed = false;
    for (let level = 0; ; level++) {
      const K = new Float64Array(g.n);
      for (let v = 0; v < g.n; v++) K[P[v]] += g.k[v];
      if (moveNodes(g, P, K, r, rand, sc)) changed = true;
      const nComm = compact(P);
      stats.levels++;
      if (nComm === g.n) break;                    // every community is one node: done
      K.fill(0);                                   // ids were compacted: recount strengths
      for (let v = 0; v < g.n; v++) K[P[v]] += g.k[v];
      const ref = refine(g, P, K, r, theta, rand, sc);
      const agg = aggregate(g, ref);
      if (agg.g.n === g.n) break;                  // refinement merged nothing (rare); stop here
      const Pn = new Uint32Array(agg.g.n);
      for (let v = 0; v < g.n; v++) Pn[agg.map[v]] = P[v];
      for (let i = 0; i < base.n; i++) nodeOf[i] = agg.map[nodeOf[i]];
      g = agg.g; P = Pn;
    }
    for (let i = 0; i < base.n; i++) labels[i] = P[nodeOf[i]];
    return changed;
  }

  // Split any cluster that is not connected in the base graph into its components.
  // Leiden already guarantees connected clusters; this is a cheap safety net.
  function splitDisconnected(base, labels) {
    const { n, ptr, adj } = base, out = new Uint32Array(n).fill(0xffffffff), stack = new Uint32Array(n);
    const firstComp = new Map();
    let next = 0, splits = 0;
    for (let s = 0; s < n; s++) {
      if (out[s] !== 0xffffffff) continue;
      const lab = labels[s], id = next++;
      if (firstComp.has(lab)) splits++; else firstComp.set(lab, id);
      let top = 0; stack[top++] = s; out[s] = id;
      while (top) {
        const v = stack[--top];
        for (let p = ptr[v]; p < ptr[v + 1]; p++) {
          const u = adj[p];
          if (out[u] === 0xffffffff && labels[u] === lab) { out[u] = id; stack[top++] = u; }
        }
      }
    }
    labels.set(out);
    return splits;
  }

  function leiden(n, indptr, indices, weights, opts) {
    const o = opts || {};
    const resolution = o.resolution == null ? 1 : +o.resolution;
    const seed = o.seed == null ? 42 : o.seed;
    const maxIter = o.maxIter == null ? -1 : o.maxIter;
    const theta = o.theta == null ? 0.01 : o.theta;
    const stats = { iterations: 0, levels: 0, splits: 0, quality: 0 };
    const base = baseGraph(n, indptr, indices, weights);
    let twoM = 0; for (let v = 0; v < n; v++) twoM += base.k[v];
    const labels = new Uint32Array(n);
    for (let v = 0; v < n; v++) labels[v] = v;
    if (twoM > 0 && n > 0) {
      const r = resolution / twoM, rand = mulberry32(seed), sc = scratch(n);
      const limit = maxIter < 0 ? 100 : maxIter;
      for (let it = 0; it < limit; it++) {
        stats.iterations++;
        if (!iterate(base, labels, r, theta, rand, sc, stats)) break;
      }
    }
    stats.splits = splitDisconnected(base, labels);
    // renumber: 0 = largest, ties by first node
    const nC = compact(labels), size = new Uint32Array(nC);
    for (let v = 0; v < n; v++) size[labels[v]]++;
    if (nC > 0xffff) throw new Error(`Leiden found ${nC} clusters; at most 65535 fit in the labels.`);
    const order = Array.from({ length: nC }, (_, c) => c).sort((a, b) => size[b] - size[a] || a - b);
    const rank = new Uint32Array(nC); order.forEach((c, i) => (rank[c] = i));
    const out = new Uint16Array(n);
    for (let v = 0; v < n; v++) out[v] = rank[labels[v]];
    stats.quality = modularity(n, base.ptr, base.adj, base.w, out, resolution);
    if (o.stats) Object.assign(o.stats, stats);
    return out;
  }

  // Q = (1/2m) sum_ij (A_ij - resolution k_i k_j / 2m) [c_i == c_j]  (igraph's modularity)
  function modularity(n, indptr, indices, weights, labels, resolution) {
    const res = resolution == null ? 1 : resolution;
    let nC = 0; for (let v = 0; v < n; v++) if (labels[v] + 1 > nC) nC = labels[v] + 1;
    const K = new Float64Array(nC);
    let twoM = 0, inW = 0;
    for (let v = 0; v < n; v++) {
      for (let p = indptr[v]; p < indptr[v + 1]; p++) {
        const x = weights[p];
        K[labels[v]] += x; twoM += x;
        if (labels[indices[p]] === labels[v]) inW += x;
      }
    }
    if (!twoM) return 0;
    let kk = 0; for (let c = 0; c < nC; c++) kk += K[c] * K[c];
    return (inW - res * kk / twoM) / twoM;
  }

  const api = { leiden, modularity };
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  if (root) { root.leiden = leiden; root.leidenModularity = modularity; }
})(typeof self !== "undefined" ? self : typeof globalThis !== "undefined" ? globalThis : this);
