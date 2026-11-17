#!/usr/bin/env node
// Draw a region test (Chromium, Playwright): the UMAP + Leiden worker and its fast mode.
//   worker  lasso_worker.js on a seeded synthetic region (2,400 bins in 2 shapes, 4 populations, 60 genes), built
//           here, so the check does not hang on any export:
//           - the exact path gives the stored layout and clusters (tests/fixtures/drawregion_exact.json, written by
//             the worker as it was before its speed work, commit b7eb968): the same embedding to float noise
//             (max |difference| under 1e-4 of the layout's extent) and the same clusters;
//           - a second worker gives the same bits; the same region again with only the resolution changed
//             re-clusters the kept layout (timings.reusedLayout) and keeps its embedding;
//           - fast mode (its threshold lowered to 1,000 bins, a 1,000-bin sample): every bin comes back, placed and
//             clustered, the sample marked; seeded (a second worker gives the same bits); its clusters agree with
//             the exact path's (ARI).
//   page    demo_visiumhd from the smoke test's export (skipped without it): a region over the whole section with
//           fast mode's threshold lowered from the console (window.__spatialViewer.lassoFast): the panel, the message
//           and Save image's legend say it is a sample; Save JSON records the mode (sampling, result.sampling,
//           result.fast with the sample's rows) and Save CSV marks the sampled bins; Run on every bin re-runs on
//           every bin; a link to the view keeps the choice, and opens with it; the saved JSON loads back to the
//           same result.
//
//   node tests/drawregion_test.js [--export demo/export] [--page web/SpatialViewer] [--shots <dir>]
//   node tests/drawregion_test.js --write-fixture --page <a folder with b7eb968's lasso_worker.js, leiden.js, libs.js>
//
// Playwright as for tests/browsers_test.js (NODE_PATH or --playwright; PLAYWRIGHT_BROWSERS_PATH). Without it, or
// without the umap-js script from the CDN (no internet), the test says SKIPPED and exits 0.
"use strict";
const http = require("http");
const fs = require("fs");
const path = require("path");
const { sendFile, preflight } = require("./static_files");

const arg = (k, d) => { const i = process.argv.indexOf(k); return i > 0 ? process.argv[i + 1] : d; };
const REPO = path.resolve(__dirname, "..");
const EXPORT = path.resolve(arg("--export", path.join(REPO, "demo", "export")));
const PAGE = path.resolve(arg("--page", path.join(REPO, "web", "SpatialViewer")));
const FIXTURE = path.join(__dirname, "fixtures", "drawregion_exact.json");
const WRITE = process.argv.includes("--write-fixture");
const SHOTS = arg("--shots", "");

let pw = null;
for (const base of [arg("--playwright", ""), ...(process.env.NODE_PATH || "").split(path.delimiter), path.join(REPO, "node_modules")].filter(Boolean)) {
  try { pw = require(path.join(base, "playwright")); break; } catch (e) { /* next */ }
}
if (!pw) { try { pw = require("playwright"); } catch (e) { /* none */ } }
if (!pw) { console.log("SKIPPED: playwright is not installed (see tests/browsers_test.js)"); process.exit(0); }

let fails = 0, oks = 0;
function check(ok, what, detail) {
  if (ok) oks++; else fails++;
  console.log(`  ${ok ? "ok  " : "FAIL"} ${what}${!ok && detail !== undefined ? ` -- ${typeof detail === "string" ? detail : JSON.stringify(detail).slice(0, 400)}` : ""}`);
  return ok;
}

// ---- the synthetic region: one hvg/ chunk in lasso_worker.js gather()'s layout ---------------------------
function mulberry32(seed) {
  let s = seed >>> 0;
  return () => { s = (s + 0x6d2b79f5) | 0; let t = Math.imul(s ^ (s >>> 15), 1 | s); t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t; return ((t ^ (t >>> 14)) >>> 0) / 4294967296; };
}
function poisson(rand, lam) { const L = Math.exp(-lam); let k = 0, p = 1; do { k++; p *= rand(); } while (p > L); return k - 1; }
function syntheticRegion() {
  const rand = mulberry32(20261117), G = 60, nRows = 3000;
  const rows = [], libs = [], ip = [0], gi = [], ct = [];
  const rowShape = new Uint16Array(nRows).fill(0xffff);
  for (let r = 0; r < nRows; r++) {
    if (r % 5 === 0) continue;                                       // a bin outside the shapes (not in the chunk)
    const pop = Math.floor(rand() * 4), depth = 0.5 + rand();
    let lib = 0;
    for (let g = 0; g < G; g++) {
      const lam = depth * (Math.floor(g / 10) === pop ? 3 : g >= 40 ? 1 : 0.3);
      const c = poisson(rand, lam);
      if (c) { gi.push(g); ct.push(c); lib += c; }
    }
    rows.push(r); libs.push(lib * 1.5 + 10); ip.push(gi.length);
    rowShape[r] = r < nRows / 2 ? 0 : 1;
  }
  const n = rows.length, nnz = gi.length;
  const buf = new ArrayBuffer(4 + 4 * n + 4 * n + 4 * (n + 1) + 2 * nnz + 2 * nnz);
  const dv = new DataView(buf); dv.setUint32(0, n, true);
  let o = 4;
  new Uint32Array(buf, o, n).set(rows); o += 4 * n;
  new Float32Array(buf, o, n).set(libs); o += 4 * n;
  new Uint32Array(buf, o, n + 1).set(ip); o += 4 * (n + 1);
  const g16 = new Uint16Array(nnz); g16.set(gi); new Uint8Array(buf, o, 2 * nnz).set(new Uint8Array(g16.buffer)); o += 2 * nnz;
  const c16 = new Uint16Array(nnz); c16.set(ct); new Uint8Array(buf, o, 2 * nnz).set(new Uint8Array(c16.buffer));
  return { chunk: Array.from(new Uint8Array(buf)), rowShape: Array.from(rowShape), nGenes: G, n };
}
// the page's parameters (viewer.js LASSO) for a run of bins
const PARAMS = { mode: "bins", seed: 42, minUmi: 10, cap: 50000, minPoints: 50, nPcs: 30, nPcsProtein: 15, nNeighbors: 15, minDist: 0.3,
  minShapesForUmap: 15, minDetect: 5, minDetectFrac: 0.01, minDetectCap: 30 };

function ari(a, b) {
  const c2 = (x) => (x * (x - 1)) / 2, t = new Map(), ra = new Map(), rb = new Map();
  for (let i = 0; i < a.length; i++) { const k = a[i] + "," + b[i]; t.set(k, (t.get(k) || 0) + 1); ra.set(a[i], (ra.get(a[i]) || 0) + 1); rb.set(b[i], (rb.get(b[i]) || 0) + 1); }
  let sij = 0, sa = 0, sb = 0; t.forEach((v) => (sij += c2(v))); ra.forEach((v) => (sa += c2(v))); rb.forEach((v) => (sb += c2(v)));
  const exp = (sa * sb) / c2(a.length), max = (sa + sb) / 2;
  return max === exp ? 1 : (sij - exp) / (max - exp);
}

// ---- servers: the page folder (plus a blank test page), and the export ------------------------------------
const TYPES = { ".html": "text/html", ".js": "text/javascript", ".css": "text/css", ".json": "application/json", ".png": "image/png", ".jpg": "image/jpeg", ".bin": "application/octet-stream", ".geojson": "application/json" };
const inside = (root, rel) => { const f = path.join(root, rel); return f.startsWith(root + path.sep) ? f : null; };
function serve(route) {
  return new Promise((resolve) => {
    const srv = http.createServer((req, res) => {
      if (req.method === "OPTIONS") { preflight(req, res); return; }
      const u = decodeURIComponent(req.url.split("?")[0]);
      if (u === "/__drawregion_test.html") { res.writeHead(200, { "Content-Type": "text/html" }); res.end("<!doctype html><title>Draw a region test</title>"); return; }
      const file = route(u);
      if (!file || !fs.existsSync(file) || fs.statSync(file).isDirectory()) { res.writeHead(404, { "Access-Control-Allow-Origin": "*" }); res.end("not found"); return; }
      sendFile(req, res, file, TYPES[path.extname(file)]);
    });
    srv.listen(0, "127.0.0.1", () => resolve(srv));
  });
}

// one run of the worker in the page: a fresh worker unless `keep` (the same one, for the kept layout)
const RUN_IN_PAGE = async ([url, msg, keep]) => {
  const w = keep && window.__w ? window.__w : new Worker(url);
  if (keep) window.__w = w;
  const chunk = new Uint8Array(msg.chunk).buffer;
  const m = { chunks: [chunk], rowShape: Uint16Array.from(msg.rowShape), nShapes: 2, nGenes: msg.nGenes, params: msg.params,
              resolution: msg.params.resolution, geneNames: Array.from({ length: msg.nGenes }, (_, g) => `G${g}`), markerSet: null,
              protein: null, flags: null, qcKeep: null, dataKey: "drawregion_test" };
  const t0 = performance.now();
  const r = await new Promise((resolve) => {
    w.onmessage = (ev) => { if (ev.data.type !== "progress") resolve(ev.data); };
    w.onerror = (e) => resolve({ type: "error", message: e.message || "worker error" });
    w.postMessage(m);
  });
  if (!keep) w.terminate();
  if (r.type === "error") return { error: r.message };
  const I = r.info;
  return { secs: (performance.now() - t0) / 1000, rows: Array.from(r.rows), emb: Array.from(r.emb), cluster: r.cluster ? Array.from(r.cluster) : null,
           info: { nIn: I.nIn, nClusters: I.nClusters, resolution: I.resolution, nGenesUsed: I.nGenesUsed, nPcs: I.nPcs, timings: I.timings || null,
                   fast: I.fast ? { sample: I.fast.sample, of: I.fast.of, k: I.fast.k, nSampled: I.fast.isSample.reduce((a, b) => a + b, 0), isSample: Array.from(I.fast.isSample) } : null } };
};
const sameBits = (a, b) => !!a && !!b && a.rows.length === b.rows.length && a.emb.every((v, i) => v === b.emb[i]) && a.rows.every((v, i) => v === b.rows[i]) &&
  (!a.cluster || a.cluster.every((v, i) => v === b.cluster[i]));

async function workerChecks(browser, pageUrl, workerUrl) {
  const page = await browser.newPage();
  const errs = [];
  page.on("pageerror", (e) => errs.push(String(e)));
  await page.goto(`${pageUrl}/__drawregion_test.html`);
  const R = syntheticRegion();
  const run = (params, keep) => page.evaluate(RUN_IN_PAGE, [workerUrl, { ...R, params: { ...PARAMS, ...params } }, !!keep]);
  const exact = await run({ resolution: "auto", sampling: "exact" });
  if (exact.error && /umap|importScripts|NetworkError|Failed to fetch|load/i.test(exact.error)) {
    console.log(`SKIPPED: the worker could not load umap-js (no internet?): ${exact.error}`); await page.close(); return "skip";
  }
  if (!check(!exact.error && exact.rows.length === R.n, `exact path ran on ${R.n} synthetic bins in ${exact.secs ? exact.secs.toFixed(1) : "?"} s`, exact.error || exact.rows.length)) { await page.close(); return; }
  if (WRITE) {
    const r6 = (x) => Math.round(x * 1e6) / 1e6;
    fs.mkdirSync(path.dirname(FIXTURE), { recursive: true });
    fs.writeFileSync(FIXTURE, JSON.stringify({
      about: "tests/drawregion_test.js: lasso_worker.js on its synthetic region (seed 42, resolution auto), written by the worker before its speed work (commit b7eb968), in Chromium",
      written: new Date().toISOString().slice(0, 10), browser: browser.version(), n: R.n, nClusters: exact.info.nClusters, resolution: exact.info.resolution,
      nGenesUsed: exact.info.nGenesUsed, nPcs: exact.info.nPcs, rows: exact.rows, cluster: exact.cluster, emb: exact.emb.map(r6) }) + "\n");
    console.log(`  wrote ${FIXTURE}: ${R.n} bins, ${exact.info.nClusters} clusters at resolution ${exact.info.resolution}`);
    await page.close(); return;
  }
  const F = JSON.parse(fs.readFileSync(FIXTURE, "utf8"));
  let lo = Infinity, hi = -Infinity, maxd = 0;
  for (const v of F.emb) { lo = Math.min(lo, v); hi = Math.max(hi, v); }
  const rowsSame = F.rows.length === exact.rows.length && F.rows.every((v, i) => v === exact.rows[i]);
  if (rowsSame) for (let i = 0; i < F.emb.length; i++) maxd = Math.max(maxd, Math.abs(F.emb[i] - exact.emb[i]));
  check(rowsSame && maxd <= 1e-4 * (hi - lo), `exact path = the stored fixture (before the speed work): max |embedding difference| ${maxd.toExponential(2)} of an extent ${(hi - lo).toFixed(2)}`, { rowsSame, maxd });
  const clEq = rowsSame && exact.cluster && F.cluster.every((v, i) => v === exact.cluster[i]);
  check(clEq && exact.info.nClusters === F.nClusters && exact.info.resolution === F.resolution,
    `exact path: the fixture's clusters (${F.nClusters} at resolution ${F.resolution}), bin for bin`, { got: [exact.info.nClusters, exact.info.resolution], ari: exact.cluster ? ari(F.cluster, exact.cluster) : null });
  const again = await run({ resolution: "auto", sampling: "exact" });
  check(sameBits(exact, again), "exact path: a second worker gives the same bits");
  // the kept layout: the same region with another resolution re-clusters only
  const k1 = await run({ resolution: "auto", sampling: "exact" }, true);
  const k2 = await run({ resolution: 0.5, sampling: "exact" }, true);
  check(!k2.error && k2.info.timings && k2.info.timings.reusedLayout === 1 && k1.emb.every((v, i) => v === k2.emb[i]) && k2.cluster && k2.info.resolution === 0.5,
    `the same region at resolution 0.5: the layout is reused (${k1.secs.toFixed(1)} s, then ${k2.secs.toFixed(1)} s) and the embedding is unchanged`, k2.error || k2.info.timings);
  const k3 = await run({ resolution: "auto", sampling: "exact", seed: 7 }, true);
  check(!k3.error && !(k3.info.timings && k3.info.timings.reusedLayout) && !sameBits(k1, k3), "another seed: a new layout, not the kept one", k3.info && k3.info.timings);
  // fast mode, thresholds lowered so the synthetic region is above them
  const fp = { resolution: "auto", sampling: "auto", fastAbove: 1000, fastSample: 1000 };
  const fast = await run(fp);
  const F2 = fast.info && fast.info.fast;
  check(!fast.error && F2 && F2.sample === 1000 && F2.of === R.n && F2.nSampled === 1000 && fast.rows.length === R.n && fast.cluster && fast.cluster.length === R.n &&
        fast.cluster.every((c) => c < fast.info.nClusters) && fast.emb.every(Number.isFinite),
    `fast mode above the threshold: a 1,000-bin sample, all ${R.n} bins placed and clustered (${fast.info && fast.info.nClusters} clusters, ${fast.secs ? fast.secs.toFixed(1) : "?"} s)`, fast.error || F2);
  const fast2 = await run(fp);
  check(sameBits(fast, fast2) && fast2.info.fast.isSample.every((v, i) => v === F2.isSample[i]), "fast mode: a second worker gives the same bits and the same sample");
  if (fast.cluster && exact.cluster) {
    const a = ari(exact.cluster, fast.cluster);
    check(a >= 0.8, `fast mode's clusters agree with the exact path's: ARI ${a.toFixed(3)} (at least 0.8)`);
  }
  const ex = await run({ ...fp, sampling: "exact" });
  check(!ex.error && !ex.info.fast && sameBits(exact, ex), "sampling \"exact\" above the threshold: every bin, the exact path's bits");
  const below = await run({ resolution: "auto", sampling: "auto", fastAbove: 20000, fastSample: 10000 });
  check(!below.error && !below.info.fast && sameBits(exact, below), "below the threshold (20,000): the exact path, whatever sampling says");
  check(!errs.length, "no page errors (worker test)", errs);
  await page.close();
}

// ---- the page: fast mode as the user meets it -------------------------------------------------------------
async function pageChecks(browser, pageUrl, dataUrl) {
  const id = "demo_visiumhd";
  if (!fs.existsSync(path.join(EXPORT, id, "manifest.json"))) { console.log(`  page checks SKIPPED: no ${id} in ${EXPORT} (run tests/smoke_test.py first, or pass --export)`); return; }
  const ctx = await browser.newContext({ acceptDownloads: true, viewport: { width: 1400, height: 900 } });
  const page = await ctx.newPage();
  const errs = [];
  page.on("pageerror", (e) => errs.push(String(e)));
  const cons = [];
  page.on("console", (m) => { if (m.type() === "error" || m.type() === "warning") cons.push(m.text().slice(0, 200)); });
  page.on("requestfailed", (r) => cons.push(`request failed: ${r.url().slice(0, 120)} ${r.failure() ? r.failure().errorText : ""}`));
  page.on("dialog", (d) => d.dismiss().catch(() => {}));
  const ev = (f, a) => page.evaluate(f, a);
  // a function, never a string: the page's CSP has no 'unsafe-eval'
  const until = async (f, ms, arg) => { try { await page.waitForFunction(f, arg === undefined ? null : arg, { timeout: ms, polling: 250 }); return true; } catch (e) { cons.push(`wait: ${String(e.message).split("\n")[0]}`); return false; } };
  await page.goto(`${pageUrl}/index.html?static=${dataUrl}&ds=${id}`);
  // Draw a region ready: positions and the hvg/ index loaded (a run started before that is dropped when they land)
  const ready = await until(() => { const S = window.__spatialViewer && window.__spatialViewer.state, L = S && S.lasso;
    return S && S.m && S.exprKey && L && L.pos && L.pos.length > 0 && L.hvg && !document.getElementById("panelLasso").hidden; }, 90000);
  if (!check(ready, `${id} opened, Draw a region ready`)) { await ctx.close(); return; }
  const nPos = await ev(() => window.__spatialViewer.state.lasso.pos.length / 2);
  // a rectangle around every bin, so the region is the whole section
  const box = await ev(() => { const p = window.__spatialViewer.state.lasso.pos; let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity;
    for (let i = 0; i < p.length; i += 2) { x0 = Math.min(x0, p[i]); x1 = Math.max(x1, p[i]); y0 = Math.min(y0, p[i + 1]); y1 = Math.max(y1, p[i + 1]); }
    return [[x0 - 5, y0 - 5], [x1 + 5, y0 - 5], [x1 + 5, y1 + 5], [x0 - 5, y1 + 5]]; });
  // thresholds well under the section's bins
  const above = Math.max(100, Math.floor(nPos / 4)), sample = Math.max(60, Math.floor(nPos / 8));
  await ev(([a, s]) => window.__spatialViewer.lassoFast(a, s), [above, sample]);
  const runDone = (prev) => until((p) => { const L = window.__spatialViewer.state.lasso; return L.result && !L.running && L.result.id !== p; }, 240000, prev || 0);
  await ev((pts) => window.__spatialViewer.lassoLoad({ shapes: [{ name: "Whole section", kind: "free", points: pts }], seed: 42, resolution: "auto", mode: "bins" }), box);
  if (!check(await runDone(0), `Run UMAP on the whole section (${nPos} bins; fast mode above ${above}, a ${sample}-bin sample)`,
    { msg: await ev(() => document.getElementById("lassoMsg").textContent), console: cons.filter((c) => /^wait:|rror/.test(c)).slice(-4) })) { await ctx.close(); return; }
  const A = await ev(() => { const L = window.__spatialViewer.state.lasso, R = L.result;
    return { id: R.id, n: R.rows.length, nIn: R.info.nIn, fast: R.info.fast ? { sample: R.info.fast.sample, of: R.info.fast.of } : null, sampling: L.sampling,
             emb: Array.from(R.emb), cluster: R.cluster ? Array.from(R.cluster) : null, nCl: R.info.nClusters,
             title: document.getElementById("lpTitle").textContent, warn: document.getElementById("lpWarn").textContent, warnHidden: document.getElementById("lpWarn").hidden,
             button: !!document.getElementById("lpExact"), msg: document.getElementById("lassoMsg").textContent }; });
  if (!check(A.fast && A.fast.sample === sample && A.fast.of === A.nIn && A.n === A.nIn, `fast mode: the UMAP of a ${sample}-bin sample, all ${A.nIn} bins in the result`, A.fast)) { await ctx.close(); return; }
  check(/fast mode/.test(A.title), `the panel's title says fast mode: "${A.title}"`);
  check(!A.warnHidden && /sample; other bins placed by nearest neighbours/.test(A.warn) && A.button, `the panel says it is a sample, with a Run on every bin button: "${A.warn.slice(0, 120)}"`);
  check(/Fast mode: UMAP/.test(A.msg) && /use every bin/.test(A.msg), "the run's message says fast mode, and that markers, Save CSV and volcano plots use every bin");
  if (SHOTS) { fs.mkdirSync(SHOTS, { recursive: true }); await page.screenshot({ path: path.join(SHOTS, "drawregion_fast.png") }); }
  const leg = await ev(() => window.__spatialViewer.snapLegend());
  const dr = leg.find((x) => /^Drawn (region|shapes)/.test(x.title));
  check(dr && dr.items.some((it) => /sample;/.test(it.label)) && dr.items.some((it) => /placed by nearest neighbours/.test(it.label)),
    "Save image's legend says the layout is of a sample", leg.map((x) => x.title));
  // Save JSON / Save CSV: what was written
  await ev(() => { window.__blobs = []; const orig = URL.createObjectURL.bind(URL); URL.createObjectURL = (b) => { window.__blobs.push(b); return orig(b); }; });
  const saved = async (btn) => { await ev(() => (window.__blobs = [])); await ev((b) => document.getElementById(b).click(), btn); await until(() => window.__blobs.length > 0, 20000); return ev(() => window.__blobs[0].text()); };
  const J = JSON.parse(await saved("lassoSaveJson"));
  const jf = J.result && J.result.fast;
  check(J.sampling === "auto" && J.result && J.result.sampling === "fast" && jf && jf.sample === sample && jf.of === A.nIn && Array.isArray(jf.sampleRows) && jf.sampleRows.length === sample &&
        J.params.fastAbove === above && J.params.fastSample === sample,
    `Save JSON records fast mode: sampling "${J.sampling}", result.sampling "${J.result && J.result.sampling}", the ${jf ? jf.sampleRows.length : 0} sampled rows, the thresholds`, { sampling: J.sampling, result: J.result && J.result.sampling, params: J.params });
  const csv = (await saved("lassoSaveCsv")).trim().split("\n");
  const col = csv[0].split(",").indexOf("in_umap_sample");
  const ones = col < 0 ? -1 : csv.slice(1).filter((l) => l.split(",").slice(-1)[0] === "1").length;
  check(col > 0 && csv.length - 1 === A.nIn && ones === sample, `Save CSV: every bin (${csv.length - 1}), in_umap_sample = 1 on the ${ones} sampled ones`, csv[0]);
  // Run on every bin
  await ev(() => document.getElementById("lpExact").click());
  check(await runDone(A.id), "Run on every bin finished");
  const B = await ev(() => { const L = window.__spatialViewer.state.lasso, R = L.result;
    return { id: R.id, n: R.rows.length, fast: !!R.info.fast, sampling: L.sampling, title: document.getElementById("lpTitle").textContent, button: !!document.getElementById("lpExact") }; });
  check(!B.fast && B.n === A.nIn && B.sampling === "exact" && !/fast mode/.test(B.title) && !B.button, `Run on every bin: the exact path on all ${B.n} bins, no sample note`, B);
  const JB = JSON.parse(await saved("lassoSaveJson"));
  check(JB.sampling === "exact" && JB.result.sampling === "exact" && !JB.result.fast, `Save JSON after it: sampling "${JB.sampling}", result.sampling "${JB.result.sampling}"`);
  // the link to this view keeps the choice, and a page opened from it has it
  const link = await ev(async () => { const r = await window.SodaViewLink.make(); const p = r.url.split("#v=")[1]; const d = await window.SodaViewLink.decode(p); return { url: r.url, state: d.ok ? d.state : null }; });
  check(link.state && link.state.sampling === "exact", `the link to this view carries sampling "exact" (${link.url.length} characters)`, link.state && Object.keys(link.state));
  const p2 = await ctx.newPage();
  p2.on("pageerror", (e) => errs.push(String(e)));
  await p2.goto(link.url);
  const opened = await p2.waitForFunction(() => { const S = window.__spatialViewer && window.__spatialViewer.state; return S && S.lasso && S.lasso.shapes.length === 1 && S.lasso.sampling === "exact"; }, null, { timeout: 60000, polling: 250 }).then(() => true, () => false);
  check(opened, "the link opens with the shape and sampling \"exact\" (for the next Run UMAP)", await p2.evaluate(() => { const L = window.__spatialViewer && window.__spatialViewer.state.lasso; return L ? { shapes: L.shapes.length, sampling: L.sampling } : null; }));
  await p2.close();
  // the saved fast-mode JSON, loaded back: the same result
  await ev((j) => window.__spatialViewer.lassoLoad(j), J);
  check(await runDone(B.id), "the fast-mode Save JSON, loaded back, runs");
  const C = await ev(() => { const R = window.__spatialViewer.state.lasso.result; return { fast: R.info.fast ? R.info.fast.sample : null, emb: Array.from(R.emb), cluster: R.cluster ? Array.from(R.cluster) : null }; });
  check(C.fast === sample && C.emb.length === A.emb.length && C.emb.every((v, i) => v === A.emb[i]) && (!A.cluster || C.cluster.every((v, i) => v === A.cluster[i])),
    "and gives the same fast-mode result, bit for bit");
  check(!errs.length, "no page errors (page test)", errs);
  await ctx.close();
}

(async () => {
  const pageSrv = await serve((u) => inside(PAGE, u === "/" ? "index.html" : u));
  const dataSrv = await serve((u) => inside(EXPORT, u));
  const pageUrl = `http://127.0.0.1:${pageSrv.address().port}`, dataUrl = `http://127.0.0.1:${dataSrv.address().port}`;
  let browser;
  try { browser = await pw.chromium.launch(); } catch (e) { console.log(`SKIPPED: Chromium could not be launched: ${e.message.split("\n")[0]}`); process.exit(0); }
  const wsrc = fs.readFileSync(path.join(PAGE, "lasso_worker.js"), "utf8");
  console.log(`Draw a region test: Chromium ${browser.version()}, page ${PAGE}`);
  let skip = false;
  try {
    skip = (await workerChecks(browser, pageUrl, `lasso_worker.js?v=test${wsrc.length}`)) === "skip";
    if (!WRITE && !skip) await pageChecks(browser, pageUrl, dataUrl);
  } catch (e) {
    check(false, "the test ran to the end", e.stack || String(e));
  }
  await browser.close(); pageSrv.close(); dataSrv.close();
  if (skip) process.exit(0);
  console.log(`Draw a region test: ${oks} passed, ${fails} failed`);
  process.exit(fails ? 1 : 0);
})();
