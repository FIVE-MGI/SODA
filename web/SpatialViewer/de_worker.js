// Differential expression for the volcano plot (volcano.js), off the main thread.
// The page hands over the hvg/ chunks under the two groups and a group code per
// positions.bin row (0 = neither, 1 = A, 2 = B); optionally per-gene files for genes
// outside the section's variable genes. Steps:
//   gather (bins with >= minUmi UMI) -> row order -> seeded cap per group ->
//   log1p CP10k per non-zero (library = all genes, from the chunk) ->
//   per gene: Wilcoxon rank-sum with tie correction (only the non-zeros are sorted: the zeros
//   are one tie block) or Welch's t on log1p CP10k -> log p (never underflows) ->
//   Benjamini-Hochberg in log space over the genes tested.
// log2 fold change (plotted) = log2((mean CP10k A + 1) / (mean CP10k B + 1));
// log2fcScanpy = scanpy's log2((expm1(mean log1p A) + 1e-9) / (expm1(mean log1p B) + 1e-9)).
// The core is plain functions (DE) so Node can run the same code for the scanpy check.
/* eslint-env worker */

const DE = (() => {
  const LN10 = Math.LN10;

  function mulberry32(seed) {
    let s = seed >>> 0;
    return () => {
      s = (s + 0x6d2b79f5) | 0;
      let t = Math.imul(s ^ (s >>> 15), 1 | s);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }

  // ---- special functions in log space ---------------------------------------------
  // log erfc(x), x >= 0 (Numerical Recipes erfcc: fractional error < 1.2e-7 for every x,
  // written as a log so it never underflows)
  function logErfc(x) {
    const t = 1 / (1 + 0.5 * x);
    return Math.log(t) - x * x - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 +
      t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277))))))));
  }
  // log of the two-sided normal p for a z score
  const logP2Normal = (z) => Math.min(0, Math.log(2) + Math.log(0.5) + logErfc(Math.abs(z) / Math.SQRT2));

  // Lanczos log-gamma (g = 7, n = 9)
  const LG = [0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313,
    -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];
  function lgamma(x) {
    if (x < 0.5) return Math.log(Math.PI / Math.abs(Math.sin(Math.PI * x))) - lgamma(1 - x);
    x -= 1;
    let a = LG[0];
    const t = x + 7.5;
    for (let i = 1; i < 9; i++) a += LG[i] / (x + i);
    return 0.5 * Math.log(2 * Math.PI) + (x + 0.5) * Math.log(t) - t + Math.log(a);
  }
  // continued fraction for the incomplete beta (Numerical Recipes betacf)
  function betacf(a, b, x) {
    const FPMIN = 1e-300;
    let qab = a + b, qap = a + 1, qam = a - 1, c = 1, d = 1 - (qab * x) / qap;
    if (Math.abs(d) < FPMIN) d = FPMIN;
    d = 1 / d;
    let h = d;
    for (let m = 1; m <= 300; m++) {
      const m2 = 2 * m;
      let aa = (m * (b - m) * x) / ((qam + m2) * (a + m2));
      d = 1 + aa * d; if (Math.abs(d) < FPMIN) d = FPMIN;
      c = 1 + aa / c; if (Math.abs(c) < FPMIN) c = FPMIN;
      d = 1 / d; h *= d * c;
      aa = (-(a + m) * (qab + m) * x) / ((a + m2) * (qap + m2));
      d = 1 + aa * d; if (Math.abs(d) < FPMIN) d = FPMIN;
      c = 1 + aa / c; if (Math.abs(c) < FPMIN) c = FPMIN;
      d = 1 / d;
      const del = d * c;
      h *= del;
      if (Math.abs(del - 1) < 3e-14) break;
    }
    return h;
  }
  // log of the regularized incomplete beta I_x(a, b)
  function logBetaInc(a, b, x) {
    if (x <= 0) return -Infinity;
    if (x >= 1) return 0;
    const lbt = lgamma(a + b) - lgamma(a) - lgamma(b) + a * Math.log(x) + b * Math.log1p(-x);
    if (x < (a + 1) / (a + b + 2)) return lbt + Math.log(betacf(a, b, x)) - Math.log(a);
    const rest = Math.exp(lbt) * betacf(b, a, 1 - x) / b;            // I = 1 - rest, not small here
    return Math.log(Math.max(1 - rest, 1e-300));
  }
  // log of the two-sided Student t p: P(|T| > |t|) = I_{df / (df + t^2)}(df / 2, 1 / 2)
  function logP2Student(t, df) {
    if (!(df > 0) || !Number.isFinite(t)) return 0;
    if (df > 1e7) return logP2Normal(t);
    return Math.min(0, logBetaInc(df / 2, 0.5, df / (df + t * t)));
  }

  // ---- reading the chunks -----------------------------------------------------------
  // views, not copies: every array of the hvg chunk layout starts on a multiple of its size
  function parseChunk(buf) {
    const n = new DataView(buf).getUint32(0, true);
    let o = 4;
    const row = new Uint32Array(buf, o, n); o += 4 * n;
    const lib = new Float32Array(buf, o, n); o += 4 * n;
    const ip = new Uint32Array(buf, o, n + 1); o += 4 * (n + 1);
    const nnz = ip[n];
    const gi = new Uint16Array(buf, o, nnz); o += 2 * nnz;
    const ct = new Uint16Array(buf, o, nnz);
    return { n, row, lib, ip, gi, ct };
  }

  // seeded choice of k of the indices in idx (partial Fisher-Yates), returned ascending
  function sample(idx, k, rand) {
    const a = idx.slice();
    for (let i = 0; i < k; i++) { const j = i + Math.floor(rand() * (a.length - i)); const t = a[i]; a[i] = a[j]; a[j] = t; }
    return a.slice(0, k).sort((x, y) => x - y);
  }

  // ---- one gene -----------------------------------------------------------------------
  // a, b: the non-zero log1p values of the gene in A and B (sorted ascending, in place);
  // nA, nB: group sizes (zeros included). Returns the rank-sum z with tie correction.
  function wilcoxonZ(a, b, nA, nB) {
    const N = nA + nB, z0A = nA - a.length, z0B = nB - b.length, z0 = z0A + z0B;
    let RA = z0A * (z0 + 1) / 2, tie = z0 > 1 ? z0 * z0 * z0 - z0 : 0, pos = z0, i = 0, j = 0;
    while (i < a.length || j < b.length) {
      const v = i < a.length && (j >= b.length || a[i] <= b[j]) ? a[i] : b[j];
      let ta = 0, tb = 0;
      while (i < a.length && a[i] === v) { ta++; i++; }
      while (j < b.length && b[j] === v) { tb++; j++; }
      const t = ta + tb, avg = pos + (t + 1) / 2;
      RA += ta * avg;
      if (t > 1) tie += t * t * t - t;
      pos += t;
    }
    const v = (nA * nB / 12) * ((N + 1) - tie / (N * (N - 1)));
    return v > 0 ? (RA - nA * (N + 1) / 2) / Math.sqrt(v) : 0;
  }

  // Welch's t on log1p values (zeros included) from the sums
  function welch(sA, qA, nA, sB, qB, nB) {
    const mA = sA / nA, mB = sB / nB;
    const vA = nA > 1 ? Math.max(0, (qA - nA * mA * mA) / (nA - 1)) : 0, vB = nB > 1 ? Math.max(0, (qB - nB * mB * mB) / (nB - 1)) : 0;
    const se2 = vA / nA + vB / nB;
    if (!(se2 > 0)) return { t: 0, df: 1 };
    const df = (se2 * se2) / ((vA / nA) ** 2 / Math.max(nA - 1, 1) + (vB / nB) ** 2 / Math.max(nB - 1, 1));
    return { t: (mA - mB) / Math.sqrt(se2), df };
  }

  // Benjamini-Hochberg on log p (natural log), over the indices in idx
  function bhLog(logp, idx) {
    const m = idx.length, out = new Float64Array(logp.length).fill(0);
    const ord = idx.slice().sort((x, y) => logp[x] - logp[y] || x - y);
    let run = 0, lm = Math.log(m);
    for (let k = m; k >= 1; k--) {
      const g = ord[k - 1], v = Math.min(0, logp[g] + lm - Math.log(k));
      run = k === m ? v : Math.min(run, v);
      out[g] = run;
    }
    return out;
  }

  // ---- the run ------------------------------------------------------------------------
  // msg: { chunks: [ArrayBuffer], group: Uint8Array per row, nGenes, geneNames,
  //        extra: [{ name, xyv: Float32Array }], positions: Float32Array (only with extra),
  //        params: { minUmi, cap, seed, minPct, test: "wilcoxon" | "welch" } }
  function run(msg, progress = () => {}) {
    const P = msg.params, t0 = Date.now(), times = {};
    const minUmi = P.minUmi != null ? P.minUmi : 10, cap = P.cap || Infinity, minPct = P.minPct != null ? P.minPct : 0.05;
    progress("Reading bins", 0);
    const parsed = msg.chunks.map(parseChunk);
    // refs to every bin of A or B: (row, chunk, index in chunk)
    const refRow = [], refC = [], refI = [], found = [0, 0, 0], low = [0, 0, 0];
    parsed.forEach((c, ci) => {
      for (let i = 0; i < c.n; i++) {
        const g = msg.group[c.row[i]];
        if (g !== 1 && g !== 2) continue;
        found[g]++;
        if (c.lib[i] < minUmi) { low[g]++; continue; }
        refRow.push(c.row[i]); refC.push(ci); refI.push(i);
      }
    });
    // row order before anything seeded: chunks arrive in network order
    const ord = refRow.map((_, k) => k).sort((x, y) => refRow[x] - refRow[y]);
    const byG = [[], [], []];
    for (const k of ord) byG[msg.group[refRow[k]]].push(k);
    const rand = mulberry32(P.seed != null ? P.seed : 42), sampledFrom = [0, 0, 0];
    for (const g of [1, 2]) if (byG[g].length > cap) { sampledFrom[g] = byG[g].length; byG[g] = sample(byG[g], cap, rand); }
    const nA = byG[1].length, nB = byG[2].length;
    const info = { foundA: found[1], foundB: found[2], lowA: low[1], lowB: low[2], nA, nB,
                   sampledFromA: sampledFrom[1], sampledFromB: sampledFrom[2], minUmi, cap: Number.isFinite(cap) ? cap : null, minPct, test: P.test };
    if (nA < 3 || nB < 3) throw new Error(`Each group needs at least 3 bins with >= ${minUmi} UMI; A has ${nA}, B has ${nB}.`);
    times.gather = Date.now() - t0;

    // per gene: the non-zero log1p values of each group, laid out gene by gene
    progress("Normalising", 0.1);
    const G = msg.nGenes, E = (msg.extra || []).length, GT = G + E;
    const cnt = [null, new Uint32Array(G + 1), new Uint32Array(G + 1)];
    for (const g of [1, 2]) for (const k of byG[g]) {
      const c = parsed[refC[k]], i = refI[k];
      for (let p = c.ip[i]; p < c.ip[i + 1]; p++) if (c.ct[p] > 0) cnt[g][c.gi[p] + 1]++;
    }
    const off = [null, cnt[1], cnt[2]];
    for (const g of [1, 2]) for (let j = 0; j < G; j++) off[g][j + 1] += off[g][j];
    const vals = [null, new Float64Array(off[1][G]), new Float64Array(off[2][G])];
    const sLog = [null, new Float64Array(GT), new Float64Array(GT)], qLog = [null, new Float64Array(GT), new Float64Array(GT)];
    const sCp = [null, new Float64Array(GT), new Float64Array(GT)], nz = [null, new Uint32Array(GT), new Uint32Array(GT)];
    for (const g of [1, 2]) {
      const fill = off[g].slice(0, G), V = vals[g];
      for (const k of byG[g]) {
        const c = parsed[refC[k]], i = refI[k], f = 1e4 / c.lib[i];
        for (let p = c.ip[i]; p < c.ip[i + 1]; p++) {
          const x = c.ct[p]; if (!(x > 0)) continue;
          const j = c.gi[p], cp = x * f, v = Math.log1p(cp);
          V[fill[j]++] = v; sLog[g][j] += v; qLog[g][j] += v * v; sCp[g][j] += cp; nz[g][j]++;
        }
      }
    }
    times.normalise = Date.now() - t0 - times.gather;

    // genes outside the chunks: their gene/ files, matched to rows
    const extraVals = [], extraInfo = [];
    if (E) {
      progress("Placing extra genes", 0.3);
      const pos = msg.positions, nRows = pos.length / 2;
      // kept bin of each row (A or B) and its library size
      const NONE = 0xFFFFFFFF, keptOf = new Uint32Array(nRows).fill(NONE), libOf = new Float64Array(nA + nB), grpOf = new Uint8Array(nA + nB);
      let q = 0;
      for (const g of [1, 2]) for (const k of byG[g]) { keptOf[refRow[k]] = q; libOf[q] = parsed[refC[k]].lib[refI[k]]; grpOf[q] = g; q++; }
      let hash = null;
      const keyOf = (x, y) => `${x},${y}`;
      msg.extra.forEach((ex, e) => {
        const xyv = ex.xyv, n = xyv.length / 3, rows = new Uint32Array(n);
        let r = 0, ok = true;
        for (let k = 0; k < n && ok; k++) {                  // a merge walk: both files are in row order
          const x = xyv[3 * k], y = xyv[3 * k + 1];
          while (r < nRows && (pos[2 * r] !== x || pos[2 * r + 1] !== y)) r++;
          if (r >= nRows) ok = false; else rows[k] = r++;
        }
        let unplaced = 0;
        if (!ok) {                                         // fall back to an exact position lookup
          if (!hash) { hash = new Map(); for (let i = 0; i < nRows; i++) { const kk = keyOf(pos[2 * i], pos[2 * i + 1]); if (!hash.has(kk)) hash.set(kk, i); } }
          for (let k = 0; k < n; k++) { const i = hash.get(keyOf(xyv[3 * k], xyv[3 * k + 1])); if (i === undefined) { rows[k] = NONE; unplaced++; } else rows[k] = i; }
        }
        const a = [], b = [], j = G + e;
        for (let k = 0; k < n; k++) {
          if (rows[k] === NONE) continue;
          const kb = keptOf[rows[k]]; if (kb === NONE) continue;
          const x = xyv[3 * k + 2]; if (!(x > 0)) continue;
          const cp = x * (1e4 / libOf[kb]), v = Math.log1p(cp), g = grpOf[kb];   // rounded as the chunk genes are: ties stay ties
          (g === 1 ? a : b).push(v); sLog[g][j] += v; qLog[g][j] += v * v; sCp[g][j] += cp; nz[g][j]++;
        }
        extraVals.push([Float64Array.from(a), Float64Array.from(b)]);
        extraInfo.push({ name: ex.name, points: n, unplaced, matched: ok ? "row order" : "position" });
      });
    }

    // the tests
    progress("Testing", 0.4);
    const wil = P.test !== "welch";
    const stat = new Float64Array(GT), logp = new Float64Array(GT), tested = new Uint8Array(GT), df = wil ? null : new Float64Array(GT);
    const idx = [];
    for (let j = 0; j < GT; j++) {
      const pa = nz[1][j] / nA, pb = nz[2][j] / nB;
      if (Math.max(pa, pb) < minPct || (nz[1][j] + nz[2][j]) === 0) { logp[j] = 0; continue; }
      tested[j] = 1; idx.push(j);
      if (wil) {
        let a, b;
        if (j < G) { a = vals[1].subarray(off[1][j], off[1][j + 1]); b = vals[2].subarray(off[2][j], off[2][j + 1]); }
        else [a, b] = extraVals[j - G];
        a.sort(); b.sort();
        stat[j] = wilcoxonZ(a, b, nA, nB);
        logp[j] = logP2Normal(stat[j]);
      } else {
        const w = welch(sLog[1][j], qLog[1][j], nA, sLog[2][j], qLog[2][j], nB);
        stat[j] = w.t; df[j] = w.df; logp[j] = logP2Student(w.t, w.df);
      }
      if ((j & 127) === 0) progress("Testing", 0.4 + 0.55 * j / GT);
    }
    const logpadj = bhLog(logp, idx);
    times.test = Date.now() - t0 - times.gather - times.normalise;

    const lfc = new Float64Array(GT), lfcScanpy = new Float64Array(GT), pctA = new Float64Array(GT), pctB = new Float64Array(GT);
    const cpA = new Float64Array(GT), cpB = new Float64Array(GT), log10p = new Float64Array(GT), log10padj = new Float64Array(GT);
    for (let j = 0; j < GT; j++) {
      cpA[j] = sCp[1][j] / nA; cpB[j] = sCp[2][j] / nB;
      lfc[j] = Math.log2((cpA[j] + 1) / (cpB[j] + 1));
      lfcScanpy[j] = Math.log2((Math.expm1(sLog[1][j] / nA) + 1e-9) / (Math.expm1(sLog[2][j] / nB) + 1e-9));
      pctA[j] = nz[1][j] / nA; pctB[j] = nz[2][j] / nB;
      log10p[j] = logp[j] / LN10; log10padj[j] = logpadj[j] / LN10;
    }
    info.nTested = idx.length; info.nGenes = GT; info.nExtra = E; info.extra = extraInfo;
    times.total = Date.now() - t0;
    info.ms = times;
    // the bins actually compared (rows of positions.bin, never barcodes): to light them or recheck them
    const rowsA = Uint32Array.from(byG[1], (k) => refRow[k]), rowsB = Uint32Array.from(byG[2], (k) => refRow[k]);
    return { lfc, lfcScanpy, stat, df, log10p, log10padj, pctA, pctB, cpA, cpB, tested, rowsA, rowsB, info };
  }

  return { run, wilcoxonZ, welch, logP2Normal, logP2Student, bhLog, logErfc, lgamma };
})();

if (typeof module !== "undefined" && module.exports) module.exports = DE;      // Node: the scanpy check
else if (typeof onmessage !== "undefined" || typeof self !== "undefined") {
  self.onmessage = (ev) => {
    try {
      const r = DE.run(ev.data, (step, frac) => postMessage({ type: "progress", step, frac }));
      const tr = [r.lfc, r.lfcScanpy, r.stat, r.log10p, r.log10padj, r.pctA, r.pctB, r.cpA, r.cpB, r.tested, r.rowsA, r.rowsB].map((a) => a.buffer);
      if (r.df) tr.push(r.df.buffer);
      postMessage({ type: "done", ...r }, tr);
    } catch (e) {
      postMessage({ type: "error", message: e.message || String(e) });
    }
  };
}
