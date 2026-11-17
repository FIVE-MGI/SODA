// A helper of lasso_worker.js: Leiden runs of Draw a region's stability sweep, several at once on several
// cores. Holds the run's graph and its subsampled bin lists (sent once), then answers
//   { type: "run", res, seed, sub, iters }  ->  { lab }  (sub -1 = the whole graph, else that subsample)
// with the same leiden.js and the same induced graphs as lasso_worker.js would make itself, so a run gives
// the same labels here as there.
/* global leiden */
importScripts("leiden.js?v=20261102a");
let G = null, keeps = [], subG = [];

// = lasso_worker.js inducedGraph
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

onmessage = (ev) => {
  const m = ev.data;
  try {
    if (m.type === "graph") { G = m; keeps = m.keeps || []; subG = []; return; }
    if (m.type !== "run" || !G) return;
    let lab;
    if (m.sub < 0) lab = leiden(G.n, G.indptr, G.indices, G.weights, { resolution: m.res, seed: m.seed, maxIter: m.iters });
    else {
      const sub = subG[m.sub] || (subG[m.sub] = inducedGraph(G, keeps[m.sub], G.n));
      lab = leiden(keeps[m.sub].length, sub.indptr, sub.indices, sub.weights, { resolution: m.res, seed: m.seed, maxIter: m.iters });
    }
    postMessage({ lab }, [lab.buffer]);
  } catch (e) {
    postMessage({ error: e.message || String(e) });
  }
};
