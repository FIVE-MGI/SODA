// Heatmap / dot plot / violin numbers (plots.js), off the main thread.
// The page hands over a group code per row of positions.bin (0xFFFF = in no group), the
// per-row totals (binqc/total.bin, or the hvg/ chunks' library sizes) when the export has
// them, the QC keep mask when a filter is on, and the per-gene files of the genes asked for.
//   bins kept: in a group, kept by the QC filter, and (with totals) total >= minUmi
//   value of a bin for a gene:
//     "cp10k": log1p(count x 10,000 / total)   (the volcano worker's normalisation)
//     "raw":   log1p(count)                    (no totals in the export: not depth-normalised)
//     "asis":  the stored value                (single-cell datasets: normalised at export)
//   per gene x group: mean over ALL kept bins (zeros included), fraction of bins with a
//   count > 0, mean on the linear scale; with detail: quartiles (linear interpolation, numpy's
//   default), a Gaussian kernel density (Scott's rule: h = sd x n^(-1/5), sd with n - 1; cut at
//   the data's min and max) and a seeded sample of points, each for all bins and for the
//   expressing bins only.
// A gene file has no row index: it is in row order, so one merge walk against positions.bin
// places it; if the walk fails, an exact (x, y) lookup does, and the misses are counted.
// The core is plain functions (PLOTS) so Node can run the same code for the scanpy check.
// Only unsigned typed arrays here: the de-identification check reads the signed names as IDs.
/* eslint-env worker */

const PLOTS = (() => {
  const NONE = 0xffff, NOROW = 0xffffffff;
  const KDE_EXACT = 2048, KDE_BINS = 2048;

  function mulberry32(seed) {
    let s = seed >>> 0;
    return () => {
      s = (s + 0x6d2b79f5) | 0;
      let t = Math.imul(s ^ (s >>> 15), 1 | s);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }

  // ---- order statistics on "z zeros + the sorted non-zeros" without building the array ----
  // nz ascending; zi = how many of nz are below 0 (0 for counts); element i of the whole:
  const elem = (nz, z, zi, i) => (i < zi ? nz[i] : i < zi + z ? 0 : nz[i - z]);
  // numpy.quantile(method="linear"): position p x (n - 1), linear between the neighbours
  function quantile(nz, z, zi, p) {
    const n = z + nz.length;
    if (!n) return NaN;
    const pos = p * (n - 1), lo = Math.floor(pos), hi = Math.min(n - 1, lo + 1);
    const a = elem(nz, z, zi, lo), b = elem(nz, z, zi, hi);
    return a + (b - a) * (pos - lo);
  }
  function countBelowZero(nz) { let lo = 0, hi = nz.length; while (lo < hi) { const m = (lo + hi) >> 1; if (nz[m] < 0) lo = m + 1; else hi = m; } return lo; }

  // ---- Gaussian kernel density on a grid of G points from lo to hi -------------------------
  // The points are nz (ascending) plus z zeros. Exact sum for up to KDE_EXACT non-zeros; above
  // that the non-zeros are first spread on KDE_BINS bins (linear binning: each value shared
  // between its two nearest bin centres), which changes the curve by far less than a pixel.
  function kde(nz, z, lo, hi, h, G) {
    const n = z + nz.length, out = new Float32Array(G);
    if (!(h > 0) || !(hi > lo) || n < 2) return null;
    let xs = nz, ws = null;
    if (nz.length > KDE_EXACT) {
      const a = nz[0], b = nz[nz.length - 1], step = (b - a) / (KDE_BINS - 1) || 1;
      xs = new Float64Array(KDE_BINS); ws = new Float64Array(KDE_BINS);
      for (let i = 0; i < KDE_BINS; i++) xs[i] = a + i * step;
      for (let i = 0; i < nz.length; i++) {
        const t = (nz[i] - a) / step, k = Math.min(KDE_BINS - 2, Math.floor(t)), f = t - k;
        ws[k] += 1 - f; ws[k + 1] += f;
      }
    }
    const c = 1 / (n * h * Math.sqrt(2 * Math.PI)), cut = 8 * h;
    let start = 0;
    for (let g = 0; g < G; g++) {
      const x = lo + ((hi - lo) * g) / (G - 1);
      let s = 0;
      if (z) { const u = x / h; s += z * Math.exp(-0.5 * u * u); }
      while (start < xs.length && xs[start] < x - cut) start++;      // xs ascending, x ascending
      for (let i = start; i < xs.length; i++) {
        const d = xs[i] - x;
        if (d > cut) break;
        const u = d / h, e = Math.exp(-0.5 * u * u);
        s += ws ? ws[i] * e : e;
      }
      out[g] = s * c;
    }
    return out;
  }

  // seeded sample of k of the n values (z zeros + nz), without replacement (Floyd's algorithm)
  function samplePoints(nz, z, zi, k, rand) {
    const n = z + nz.length, out = new Float32Array(Math.min(k, n));
    if (n <= k) { for (let i = 0; i < n; i++) out[i] = elem(nz, z, zi, i); return out; }
    const seen = new Set();
    let q = 0;
    for (let j = n - k; j < n; j++) {
      const t = Math.floor(rand() * (j + 1)), pick = seen.has(t) ? j : t;
      seen.add(pick); out[q++] = elem(nz, z, zi, pick);
    }
    return out;
  }

  // summary of one population: nz ascending (Float64Array), z zeros among n = z + nz.length
  function describe(nz, z, sum, sumSq, P, rand) {
    const n = z + nz.length;
    const o = { n, mean: n ? sum / n : NaN, q1: NaN, median: NaN, q3: NaN, min: NaN, max: NaN, sd: NaN, h: NaN, lo: NaN, hi: NaN, dens: null, pts: new Float32Array(0) };
    if (!n) return o;
    const zi = countBelowZero(nz);
    o.q1 = quantile(nz, z, zi, 0.25); o.median = quantile(nz, z, zi, 0.5); o.q3 = quantile(nz, z, zi, 0.75);
    o.min = elem(nz, z, zi, 0); o.max = elem(nz, z, zi, n - 1);
    o.lo = o.min; o.hi = o.max;
    if (n > 1) {
      o.sd = Math.sqrt(Math.max(0, (sumSq - n * o.mean * o.mean) / (n - 1)));
      o.h = o.sd * Math.pow(n, -0.2);                               // Scott's rule
      if (n >= P.kdeMin) o.dens = kde(nz, z, o.lo, o.hi, o.h, P.kdeGrid);
    }
    o.pts = samplePoints(nz, z, zi, P.maxPoints, rand);
    return o;
  }

  // ---- rows of a gene file's points ---------------------------------------------------------
  function makePlacer(pos) {
    const nRows = pos.length / 2;
    let hash = null;
    const keyOf = (x, y) => `${x},${y}`;
    return (xyv) => {
      const n = xyv.length / 3, rows = new Uint32Array(n);
      let r = 0, ok = true;
      for (let k = 0; k < n && ok; k++) {                // a merge walk: both files are in row order
        const x = xyv[3 * k], y = xyv[3 * k + 1];
        while (r < nRows && (pos[2 * r] !== x || pos[2 * r + 1] !== y)) r++;
        if (r >= nRows) ok = false; else rows[k] = r++;
      }
      let unplaced = 0;
      if (!ok) {                                          // fall back to an exact position lookup
        if (!hash) { hash = new Map(); for (let i = 0; i < nRows; i++) { const kk = keyOf(pos[2 * i], pos[2 * i + 1]); if (!hash.has(kk)) hash.set(kk, i); } }
        for (let k = 0; k < n; k++) { const i = hash.get(keyOf(xyv[3 * k], xyv[3 * k + 1])); if (i === undefined) { rows[k] = NOROW; unplaced++; } else rows[k] = i; }
      }
      return { rows, unplaced, matched: ok ? "row order" : "position" };
    };
  }

  // ---- the run ------------------------------------------------------------------------------
  // msg: { layout: "xyv" | "indexed", nRows, positions (xyv only), group: Uint16Array per row,
  //        nGroups, totals: Float32Array | null, keep: Uint8Array | null,
  //        genes: [{ name, buf: ArrayBuffer | null }]   (null = the gene has no counts: all zeros)
  //        params: { norm: "cp10k" | "raw" | "asis", minUmi, detail, kdeGrid, kdeMin, maxPoints, seed } }
  function run(msg, progress = () => {}) {
    const t0 = Date.now(), P = Object.assign({ norm: "cp10k", minUmi: 10, detail: false, kdeGrid: 64, kdeMin: 5, maxPoints: 300, seed: 42 }, msg.params || {});
    const N = msg.nRows, K = msg.nGroups, grp = msg.group, tot = msg.totals, keep = msg.keep;
    if (P.norm === "cp10k" && !tot) throw new Error("cp10k normalisation needs the per-bin totals");
    if (grp.length !== N) throw new Error(`group map has ${grp.length} rows, the layer ${N}`);
    progress("Finding the bins", 0);
    // the group each row counts for, after the QC filter and the minimum total
    const eff = new Uint16Array(N).fill(NONE);
    const found = new Uint32Array(K), qcOut = new Uint32Array(K), low = new Uint32Array(K), n = new Uint32Array(K);
    for (let r = 0; r < N; r++) {
      const g = grp[r];
      if (g === NONE || g >= K) continue;
      found[g]++;
      if (keep && !keep[r]) { qcOut[g]++; continue; }
      if (P.norm === "cp10k" && !(tot[r] >= P.minUmi)) { low[g]++; continue; }
      eff[r] = g; n[g]++;
    }
    const Gn = msg.genes.length;
    const mean = new Float64Array(Gn * K), meanLin = new Float64Array(Gn * K), frac = new Float64Array(Gn * K), nExpr = new Uint32Array(Gn * K);
    const geneInfo = [], detail = P.detail ? [] : null;
    const place = msg.layout === "xyv" ? makePlacer(msg.positions) : null;
    const sum = new Float64Array(K), sumSq = new Float64Array(K), sumLin = new Float64Array(K), cnt = new Uint32Array(K);

    for (let j = 0; j < Gn; j++) {
      const ge = msg.genes[j];
      progress("Averaging", (j + 0.5) / Gn);
      sum.fill(0); sumSq.fill(0); sumLin.fill(0); cnt.fill(0);
      let rows = null, vals = null, np = 0, info = { name: ge.name, points: 0, unplaced: 0, matched: "none" };
      if (ge.buf) {
        if (msg.layout === "xyv") {
          const xyv = new Float32Array(ge.buf), pl = place(xyv);
          np = xyv.length / 3; rows = pl.rows; info = { name: ge.name, points: np, unplaced: pl.unplaced, matched: pl.matched };
          vals = { at: (k) => xyv[3 * k + 2] };
        } else {                                          // uint32 k, uint32[k] cell index, float32[k] value
          np = new Uint32Array(ge.buf, 0, 1)[0];
          rows = new Uint32Array(ge.buf, 4, np);
          const v = new Float32Array(ge.buf, 4 + 4 * np, np);
          vals = { at: (k) => v[k] };
          info = { name: ge.name, points: np, unplaced: 0, matched: "index" };
        }
      }
      // pass 1: the value of each point, sums per group
      const v = new Float64Array(np), gOf = new Uint16Array(np);
      for (let k = 0; k < np; k++) {
        const r = rows[k];
        gOf[k] = NONE;
        if (r === NOROW || r >= N) continue;
        const g = eff[r]; if (g === NONE) continue;
        const x = vals.at(k);
        if (P.norm === "asis" ? x === 0 || x !== x : !(x > 0)) continue;
        let val, lin;
        if (P.norm === "cp10k") { lin = (x * 1e4) / tot[r]; val = Math.log1p(lin); }
        else if (P.norm === "raw") { lin = x; val = Math.log1p(x); }
        else { val = x; lin = P.asisLog ? Math.expm1(x) : x; }
        v[k] = val; gOf[k] = g;
        sum[g] += val; sumSq[g] += val * val; sumLin[g] += lin; cnt[g]++;
      }
      for (let g = 0; g < K; g++) {
        const o = j * K + g;
        mean[o] = n[g] ? sum[g] / n[g] : NaN; meanLin[o] = n[g] ? sumLin[g] / n[g] : NaN;
        frac[o] = n[g] ? cnt[g] / n[g] : NaN; nExpr[o] = cnt[g];
      }
      geneInfo.push(info);
      if (!P.detail) continue;
      // pass 2: each group's non-zero values, sorted
      const off = new Uint32Array(K + 1);
      for (let g = 0; g < K; g++) off[g + 1] = off[g] + cnt[g];
      const fill = off.slice(0, K), all = new Float64Array(off[K]);
      for (let k = 0; k < np; k++) if (gOf[k] !== NONE) all[fill[gOf[k]]++] = v[k];
      const per = [];
      for (let g = 0; g < K; g++) {
        const nz = all.subarray(off[g], off[g + 1]).sort();
        const rand = mulberry32((P.seed * 1000003 + j * 8191 + g) >>> 0);
        per.push({ all: describe(nz, n[g] - nz.length, sum[g], sumSq[g], P, rand), expr: describe(nz, 0, sum[g], sumSq[g], P, rand) });
      }
      detail.push(per);
    }
    const info = { norm: P.norm, minUmi: P.norm === "cp10k" ? P.minUmi : null, nGroups: K, nGenes: Gn, genes: geneInfo,
                   kde: P.detail ? { rule: "scott", grid: P.kdeGrid, minPoints: P.kdeMin, exactUpTo: KDE_EXACT, bins: KDE_BINS } : null,
                   maxPoints: P.maxPoints, seed: P.seed, ms: Date.now() - t0 };
    return { mean, meanLin, frac, nExpr, n, found, qcOut, low, eff, detail, info };
  }

  return { run, quantile, kde, describe, samplePoints, makePlacer, NONE };
})();

if (typeof module !== "undefined" && module.exports) module.exports = PLOTS;      // Node: the scanpy check
else if (typeof self !== "undefined") {
  self.onmessage = (ev) => {
    try {
      const r = PLOTS.run(ev.data, (step, frac) => postMessage({ type: "progress", step, frac }));
      postMessage({ type: "done", ...r }, [r.mean, r.meanLin, r.frac, r.nExpr, r.n, r.found, r.qcOut, r.low, r.eff].map((a) => a.buffer));
    } catch (e) {
      postMessage({ type: "error", message: e.message || String(e) });
    }
  };
}
