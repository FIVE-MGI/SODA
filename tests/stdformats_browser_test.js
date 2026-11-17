#!/usr/bin/env node
// Standard formats in the browser: opens the exports made from the synthetic spatial AnnData, OME-Zarr and
// SpatialData files (examples/make_stdformats_demo.py, web/SpatialViewer/stdformats.py) in Playwright's Chromium
// and checks that each one is an ordinary SODA dataset to the page: the image tiles are drawn, the expression
// layer loads and a gene from the search shows on the slide, the side panel matches the image (a photo shows
// "Image", a fluorescence image or a density background its channels), the file's own clusterings are there,
// the SpatialData cells have their outlines, and the page logs no error and no failed request.
//
//   node tests/stdformats_browser_test.js                  (needs demo/export: run tests/smoke_test.py first)
//   node tests/stdformats_browser_test.js --export <dir> --page <dir> --shots <dir> --ports 8868,8869 --gpu
//
// Playwright is not a dependency of this repo. Install it anywhere (npm i playwright; npx playwright install
// chromium) and point NODE_PATH or --playwright at its node_modules. Without it, or without the exports, the
// test says SKIPPED and exits 0. It draws with software rendering unless --gpu is given.
"use strict";
const http = require("http");
const fs = require("fs");
const path = require("path");
const { sendFile, preflight } = require("./static_files");

const arg = (k, d) => { const i = process.argv.indexOf(k); return i > 0 ? process.argv[i + 1] : d; };
const REPO = path.resolve(__dirname, "..");
const EXPORT = path.resolve(arg("--export", path.join(REPO, "demo", "export")));
const PAGE = path.resolve(arg("--page", path.join(REPO, "web", "SpatialViewer")));
const SHOTS = arg("--shots", "");
const PORTS = arg("--ports", "0,0").split(",").map(Number);      // data server, page server; 0 = any free port
const GPU = process.argv.includes("--gpu");

// id, what it is, the panel it must show (photo: "Image" without a channel list), its layer, what else to check
const DATASETS = [
  { id: "demo_anndata_visium", what: "spatial AnnData (.h5ad, Visium)", photo: true, layer: "spots", user: "cluster" },
  { id: "demo_anndata_cells", what: "AnnData-Zarr cells, no image", photo: false, layer: "cells" },
  { id: "demo_omezarr", what: "OME-Zarr image", photo: false, layer: null, channels: ["DAPI", "CD45", "PanCK", "CD3"] },
  { id: "demo_sdata_xenium", what: "SpatialData (Xenium-like)", photo: false, layer: "cells", outlines: true, molecules: true, user: "cell_type" },
  { id: "demo_sdata_visiumhd", what: "SpatialData (Visium HD-like)", photo: true, layer: "016um" },
  { id: "demo_omezarr_v3", what: "OME-Zarr image, OME-NGFF 0.5 (zarr v3, sharded)", photo: false, layer: null, channels: ["DAPI", "CD45", "PanCK", "CD3"] },
];

let pw = null;
for (const base of [arg("--playwright", ""), ...(process.env.NODE_PATH || "").split(path.delimiter), path.join(REPO, "node_modules")].filter(Boolean)) {
  try { pw = require(path.join(base, "playwright")); break; } catch (e) { /* next */ }
}
if (!pw) { try { pw = require("playwright"); } catch (e) { /* none */ } }
if (!pw) { console.log("SKIPPED: playwright is not installed (see the top of this file)"); process.exit(0); }
const present = DATASETS.filter((d) => fs.existsSync(path.join(EXPORT, d.id, "manifest.json")));
if (!present.length) { console.log(`SKIPPED: no standard-format export at ${EXPORT} (run tests/smoke_test.py first, or pass --export)`); process.exit(0); }

// files the page asks for that no export has to have (as in tests/compat_test.js): a lab's published gates, and the
// page's markers/<file> named by a protein grid layer exported with a markers list the page does not ship
// vendor/: deck.gl from the page's own site when `prepare.py vendor` put it there, else the CDN (libs.js)
const OPTIONAL = [/\/gates\.json$/, /\/markers\/[\w.-]+\.json$/, /\/vendor\//];
const TYPES = { ".html": "text/html", ".js": "text/javascript", ".css": "text/css", ".json": "application/json", ".png": "image/png", ".jpg": "image/jpeg", ".bin": "application/octet-stream" };
function serve(root, port) {
  return new Promise((resolve) => {
    const srv = http.createServer((req, res) => {
      const rel = decodeURIComponent(req.url.split("?")[0]);
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
const check = (cond, what, got) => { results.push({ ok: !!cond, what }); console.log(`  ${cond ? "ok  " : "FAIL"} ${what}${cond || got === undefined ? "" : "  got " + JSON.stringify(got)}`); return !!cond; };

(async () => {
  const data = await serve(EXPORT, PORTS[0]), pageSrv = await serve(PAGE, PORTS[1]);
  const DATA = `http://127.0.0.1:${data.address().port}`, URL0 = `http://127.0.0.1:${pageSrv.address().port}/index.html`;
  const browser = await pw.chromium.launch({ args: GPU ? ["--use-angle=d3d11", "--ignore-gpu-blocklist", "--enable-gpu"] : ["--use-angle=swiftshader", "--enable-unsafe-swiftshader", "--ignore-gpu-blocklist"] });
  // the service worker (sw.js) is tested by tests/offline_test.js; here the page as it loads
  const ctx = await browser.newContext({ viewport: { width: 1366, height: 768 }, serviceWorkers: "block" });
  const errors = [];
  for (const D of present) {
    console.log(`${D.id} (${D.what})`);
    const page = await ctx.newPage();
    // a failed request is reported by its URL below; the browser's own "Failed to load resource" line repeats it
    page.on("console", (m) => { if (m.type() === "error" && !/^Failed to load resource/.test(m.text())) errors.push(`${D.id}: console: ${m.text()}`); });
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
      check(st.exprKey === D.layer && st.genes === 20, `${D.id}: layer ${D.layer} with its 20 genes`, { exprKey: st.exprKey, genes: st.genes });
      const gene = await page.evaluate(() => [...window.__spatialViewer.state.geneIndex].sort((a, b) => b.t - a.t)[0].n);
      await page.fill("#geneSearch", gene); await page.keyboard.press("Enter");
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
    if (SHOTS) { fs.mkdirSync(SHOTS, { recursive: true }); await page.screenshot({ path: path.join(SHOTS, `${D.id}.png`) }); }
    await page.close();
  }
  check(errors.length === 0, "no console error, page error or failed request", errors.slice(0, 5));
  await browser.close(); data.close(); pageSrv.close();
  const bad = results.filter((r) => !r.ok);
  console.log(bad.length ? `\nFAILED: ${bad.length} of ${results.length} checks` : `\nPASSED: ${results.length} checks`);
  process.exit(bad.length ? 1 : 0);
})().catch((e) => { console.error(e); process.exit(2); });
