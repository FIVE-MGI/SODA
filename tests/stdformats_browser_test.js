#!/usr/bin/env node
// Standard and instrument formats in the browser: opens the public exports made from a SpatialData store (xenium_pub,
// web/SpatialViewer/stdformats.py) and from the CosMx, MERSCOPE, Stereo-seq and GeoMx outputs (cosmx_1, merscope_2,
// stereoseq_1, geomx_1) in Playwright's Chromium (and each other browser of --browsers, e.g. firefox) and checks that
// each one is an ordinary SODA dataset to the page: the image tiles are drawn, the expression layer loads with every
// gene of its genes.json and a gene from the search shows on the slide, the side panel matches the image (a photo
// shows "Image", a fluorescence image or a density background its channels), the file's own groupings are there,
// cells have their outlines, transcripts are the molecules layer, and the page logs no error and no failed request.
// The exports of spatial AnnData and OME-Zarr files have no public example: tests/smoke_test.py's synthetic ones went
// with the demo data. cosmx_1, stereoseq_1 and geomx_1 are restricted (tests/test_data.js): without
// SODA_TEST_PRIVATE_DATA or a local copy they are SKIPPED, each with a line saying so.
//
//   node tests/stdformats_browser_test.js                  (the public exports: SODA_TEST_DATA or --data-root, see tests/test_data.js)
//   node tests/stdformats_browser_test.js --data-root <url or folder> --page <dir> --shots <dir> --ports 8868,8869 --gpu
//        --browsers chromium,firefox   (default chromium; a browser other than Chromium that is not installed is
//                                       SKIPPED on its own, with a "SKIPPED <browser>: ..." line)
//
// Playwright is not a dependency of this repo. Install it anywhere (npm i playwright; npx playwright install
// chromium) and point NODE_PATH or --playwright at its node_modules. Without it, or without any of the exports, the
// test says SKIPPED and exits 0. Chromium draws with software rendering unless --gpu is given.
"use strict";
const unfold = (page, ids) => page.evaluate((ids) => { for (const id of ids) { const s = document.getElementById(id); if (s && s.classList.contains("closed")) document.getElementById(`secBtn-${id}`).click(); } }, ids);   // every section starts closed
const http = require("http");
const fs = require("fs");
const path = require("path");
const { sendFile, preflight } = require("./static_files");

const arg = (k, d) => { const i = process.argv.indexOf(k); return i > 0 ? process.argv[i + 1] : d; };
const REPO = path.resolve(__dirname, "..");
const DR = require("./test_data");                             // the public exports: SODA_TEST_DATA / --data-root
const PAGE = path.resolve(arg("--page", path.join(REPO, "web", "SpatialViewer")));
const SHOTS = arg("--shots", "");
const PORTS = arg("--ports", "0,0").split(",").map(Number);      // data server, page server; 0 = any free port
const GPU = process.argv.includes("--gpu");
const BROWSERS = arg("--browsers", "chromium").split(",").map((s) => s.trim()).filter(Boolean);
// Firefox's console lines that are not the page's faults (see tests/browsers_test.js): the vendor/ tag's SHA-384
// checked against a 404 when there is no vendor/, and deck.gl's WebAssembly probe refused by the page's CSP
const FIREFOX_QUIET = [/integrity attribute match the content of the subresource at .[^ ]*\/vendor\//,
  /Content-Security-Policy: .*blocked WebAssembly .*Missing 'wasm-unsafe-eval'/];

// id, what it is, the panel it must show (photo: "Image" without a channel list), its layer, what else to check
// (channels: true = the side panel lists the manifest's channel names, in order)
const DATASETS = [
  { id: "xenium_pub", what: "SpatialData (Xenium, zarr v3)", photo: false, layer: "cells", channels: true, outlines: true, molecules: true },
  { id: "cosmx_1", what: "CosMx SMI flat files", photo: false, layer: "cells", channels: true, outlines: true, molecules: true },
  { id: "merscope_2", what: "MERSCOPE, no images (density background)", photo: false, layer: "cells", channels: true, outlines: true, molecules: true },
  { id: "stereoseq_1", what: "Stereo-seq GEF + cell bin GEF, H&E", photo: true, layer: "bin020" },
  { id: "geomx_1", what: "GeoMx DSP (DCC, PKC, annotations, two scans)", photo: false, layer: "rois", user: "Tissue" },
];

let pw = null;
for (const base of [arg("--playwright", ""), ...(process.env.NODE_PATH || "").split(path.delimiter), path.join(REPO, "node_modules")].filter(Boolean)) {
  try { pw = require(path.join(base, "playwright")); break; } catch (e) { /* next */ }
}
if (!pw) { try { pw = require("playwright"); } catch (e) { /* none */ } }
if (!pw) { console.log("SKIPPED: playwright is not installed (see the top of this file)"); process.exit(0); }
let present = [];                                               // the datasets the data root has (found before the run)

// files the page asks for that no export has to have (as in tests/compat_test.js): a lab's published gates, and the
// page's markers/<file> named by a protein grid layer exported with a markers list the page does not ship
// vendor/: deck.gl from the page's own site when `prepare.py vendor` put it there, else the CDN (libs.js)
const OPTIONAL = [/\/gates\.json$/, /\/markers\/[\w.-]+\.json$/, /\/vendor\//];
const TYPES = { ".html": "text/html", ".js": "text/javascript", ".css": "text/css", ".json": "application/json", ".png": "image/png", ".jpg": "image/jpeg", ".bin": "application/octet-stream" };
// root null: the data (tests/test_data.js: a folder, or passed on to a URL)
function serve(root, port) {
  return new Promise((resolve) => {
    const srv = http.createServer((req, res) => {
      const rel = decodeURIComponent(req.url.split("?")[0]);
      if (!root) { if (req.method === "OPTIONS") { preflight(req, res); return; } DR.serve(req, res, rel.slice(1), TYPES[path.extname(rel)]); return; }
      if (rel === "/config.js") { res.writeHead(200, { "Content-Type": "text/javascript" }); res.end("window.SpatialViewerConfig = { tokenEndpoint: \"\" };"); return; }
      const file = path.join(root, rel === "/" ? "index.html" : rel);
      if (req.method === "OPTIONS") { preflight(req, res); return; }
      if (!file.startsWith(root) || !fs.existsSync(file) || fs.statSync(file).isDirectory()) { res.writeHead(404, { "Access-Control-Allow-Origin": "*" }); res.end("not found"); return; }
      sendFile(req, res, file, TYPES[path.extname(file)]);
    });
    srv.listen(port || 0, "127.0.0.1", () => resolve(srv));
  });
}

const results = [];
let TAG = "";                                    // "firefox: " before each check after Chromium's
const check = (cond, what, got) => { what = TAG + what; results.push({ ok: !!cond, what }); console.log(`  ${cond ? "ok  " : "FAIL"} ${what}${cond || got === undefined ? "" : "  got " + JSON.stringify(got)}`); return !!cond; };

// every dataset in one browser (B: its name)
async function run(B, browser, DATA, URL0) {
  // no service worker (an old one of the removed offline use must not answer requests here)
  const ctx = await browser.newContext({ viewport: { width: 1366, height: 768 }, serviceWorkers: "block" });
  const errors = [];
  for (const D of present) {
    console.log(`${TAG}${D.id} (${D.what})`);
    const page = await ctx.newPage();
    // a failed request is reported by its URL below; the browser's own "Failed to load resource" line repeats it
    page.on("console", (m) => { if (m.type() === "error" && !/^Failed to load resource/.test(m.text()) && !(B === "firefox" && FIREFOX_QUIET.some((re) => re.test(m.text())))) errors.push(`${D.id}: console: ${m.text()}`); });
    page.on("pageerror", (e) => errors.push(`${D.id}: pageerror: ${e.message}`));
    page.on("response", (r) => { if (r.status() >= 400 && !/favicon/.test(r.url()) && !OPTIONAL.some((re) => re.test(r.url().split("?")[0]))) errors.push(`${D.id}: http ${r.status()}: ${r.url()}`); });
    await page.goto(`${URL0}?static=${encodeURIComponent(DATA)}&ds=${D.id}`);
    try {
      await page.waitForFunction((id) => { const sv = window.__spatialViewer; return !!(sv && sv.state && sv.state.dsId === id && sv.state.m) && document.getElementById("status").hidden; }, D.id, { timeout: 120000 });
      await page.waitForFunction(() => { const S = window.__spatialViewer.state; return !(S.m.expr || []).length || (S.exprKey && S.geneIndex.length > 0); }, null, { timeout: 120000 });
    } catch (e) {
      check(false, `${D.id}: the page opens the dataset`, await page.evaluate(() => document.getElementById("status").textContent));
      await page.close(); continue;
    }
    await page.waitForTimeout(1500);
    const st = await page.evaluate(() => {
      const sv = window.__spatialViewer, S = sv.state, e = (S.m.expr || []).find((x) => x.key === S.exprKey) || {};
      let lit = -1;
      if (sv.deck && sv.deck.device) {
        const gl = sv.deck.device.gl; sv.deck.redraw("test");
        const w = gl.drawingBufferWidth, h = gl.drawingBufferHeight, px = new Uint8Array(w * h * 4);
        gl.bindFramebuffer(gl.FRAMEBUFFER, null); gl.readPixels(0, 0, w, h, gl.RGBA, gl.UNSIGNED_BYTE, px);
        let n = 0; for (let i = 0; i < px.length; i += 4) if (px[i] > 8 || px[i + 1] > 8 || px[i + 2] > 8) n++;
        lit = n / (w * h);
      }
      return { lit, exprKey: S.exprKey || null, genes: (S.geneIndex || []).length, chanHidden: document.getElementById("chanList").hidden,
               chanTitle: document.getElementById("chanTitle").textContent, channels: (S.channels || []).map((c) => c.name),
               user: ((e.umap || {}).clusterings || []).filter((c) => c.source === "user").map((c) => c.name),
               outlines: !!(sv.outlines && sv.outlines.available()), molecules: !!S.m.molecules };
    });
    check(st.lit > 0.02, `${D.id}: tiles drawn (${(100 * st.lit).toFixed(1)} % of the map lit)`, st.lit);
    check(st.chanHidden === D.photo, `${D.id}: the side panel shows ${D.photo ? "Image (a photo)" : "the channels"}`, { hidden: st.chanHidden, title: st.chanTitle });
    if (D.channels) check(JSON.stringify(st.channels) === JSON.stringify(D.channels), `${D.id}: channel names from the file`, st.channels);
    if (D.layer) {
      check(st.exprKey === D.layer && st.genes === D.nGenes, `${D.id}: layer ${D.layer} with its ${D.nGenes} genes (genes.json)`, { exprKey: st.exprKey, genes: st.genes });
      const gene = await page.evaluate(() => [...window.__spatialViewer.state.geneIndex].sort((a, b) => b.t - a.t)[0].n);
      await unfold(page, ["panelGenes"]); await page.fill("#geneSearch", gene); await page.keyboard.press("Enter");
      let shown = false;
      try {
        await page.waitForFunction(() => window.__spatialViewer.state.genes.length === 1 && window.__spatialViewer.state.genes[0].n > 0, null, { timeout: 60000 });
        shown = true;
      } catch (e) { /* reported below */ }
      check(shown, `${D.id}: a gene from the search shows on the slide (${gene})`);
    }
    if (D.user) check(st.user.includes(D.user), `${D.id}: the file's own '${D.user}' column is a clustering of the layer`, st.user);
    if (D.outlines) check(st.outlines, `${D.id}: the cells have outlines (Cell shapes)`);
    if (D.molecules) check(st.molecules, `${D.id}: the transcripts are the molecules layer`);
    // Chromium's in <shots>/ as before, another browser's in <shots>/<browser>/
    const dir = B === "chromium" ? SHOTS : path.join(SHOTS, B);
    if (SHOTS) { fs.mkdirSync(dir, { recursive: true }); await page.screenshot({ path: path.join(dir, `${D.id}.png`) }); }
    await page.close();
  }
  check(errors.length === 0, "no console error, page error or failed request", errors.slice(0, 5));
  await ctx.close();
}

(async () => {
  for (const d of DATASETS) {
    if (!(await DR.has(d.id))) { console.log(DR.skipLine(d.id, d.what)); continue; }
    if (d.layer) d.nGenes = ((await DR.json(`${d.id}/expr/${d.layer}/genes.json`)) || []).length;
    // the channel names the export read from the file (the page must list them as they are, in order): a real export
    // names them as its reader found them, which can change with the reader's version (xenium_pub: "0", later "DAPI")
    if (d.channels === true) d.channels = (((await DR.manifest(d.id)).image || {}).channels || []).map((c) => c.name);
    present.push(d);
  }
  if (!present.length) { console.log(`SKIPPED: none of ${DATASETS.map((d) => d.id).join(", ")} can be read from ${DR.ROOT}`); process.exit(0); }
  console.log(`data: ${DR.describe()}`);
  const data = await serve(null, PORTS[0]), pageSrv = await serve(PAGE, PORTS[1]);
  const DATA = `http://127.0.0.1:${data.address().port}`, URL0 = `http://127.0.0.1:${pageSrv.address().port}/index.html`;
  const ran = [];
  for (const B of BROWSERS) {
    if (!pw[B]) { console.log(`SKIPPED ${B}: not a Playwright browser`); continue; }
    let browser;
    try {
      browser = await pw[B].launch(B === "chromium" ? { args: GPU ? ["--use-angle=d3d11", "--ignore-gpu-blocklist", "--enable-gpu"] : ["--use-angle=swiftshader", "--enable-unsafe-swiftshader", "--ignore-gpu-blocklist"] } : {});
    } catch (e) {
      if (B === "chromium") throw e;
      console.log(`SKIPPED ${B}: cannot launch (${String(e.message || e).split("\n")[0]})`); continue;
    }
    TAG = B === "chromium" ? "" : `${B}: `;
    console.log(`== ${B} ${browser.version()}`);
    try { await run(B, browser, DATA, URL0); } finally { await browser.close(); }
    ran.push(B);
  }
  data.close(); pageSrv.close();
  const bad = results.filter((r) => !r.ok);
  console.log(bad.length ? `\nFAILED: ${bad.length} of ${results.length} checks (${ran.join(", ")})` : `\nPASSED: ${results.length} checks (${ran.join(", ")})`);
  process.exit(bad.length ? 1 : 0);
})().catch((e) => { console.error(e); process.exit(2); });
