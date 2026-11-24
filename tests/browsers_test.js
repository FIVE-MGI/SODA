#!/usr/bin/env node
// Browsers test: the page on the public example exports (tests/test_data.js) in Chromium, WebKit (Safari's engine), Microsoft Edge and Firefox,
// plus an iPhone-sized WebKit and Chromium with touch, and a phone-sized Firefox (the size only: Playwright's
// Firefox has no isMobile / hasTouch emulation, so it is a narrow desktop window, not a touch phone). The same
// steps in each browser:
//   probe     what the browser offers (WebGL2 and its limits, createImageBitmap, OffscreenCanvas, module
//             workers, BigInt64Array, CSS dvh / :has / backdrop-filter, clipboard, the largest 2D canvas
//             that still encodes a PNG); written to <shots>/<browser>/probe.json, a few must hold.
//   datasets  one public dataset of each kind (Visium HD visiumhd_pub, multiplex phenocycler_pub, Xenium
//             xenium_pub, GeoMx geomx_1, Stereo-seq stereoseq_1, CosMx cosmx_1, MERSCOPE merscope_2) opened with ?static=: tiles fetched, decoded and
//             drawn, one gene added from the search (from a gene pack: read with a Range request, 206), no
//             console error, no page error, no failed request. The first view (fit, before the gene) is compared with Chromium's pixel by pixel
//             on a 64 x 40 grid of the map (same JPEG / PNG tiles, same view: only rounding may differ).
//   walk      visiumhd_pub: Tissue / UMAP / Side by side; Draw a region with the mouse (two circles)
//             and Run UMAP (a worker that loads umap-js); Save JSON and Save CSV (file downloads);
//             volcano (its worker) and its PNG (an SVG drawn on a canvas); heatmap / violin and the
//             plot's PNG; the QC filter; Save image at 1x, 2x and 4x (the
//             PNG's size read from the file); Copy link to this view, and that link opened in a new tab
//             (same dataset, layer, genes, mode and camera).
//   compare   Compare side by side (compare.js) on ?data=<folder>/ with a datasets.json: the Compare button in
//             the top bar beside Fit, a second dataset opened in a linked window (both draw, both join the group,
//             they split the tile cache: ~512 MB each), Link views moves the other window to the same place in µm
//             at the same scale (A to B and B to A), a gene added in one follows into the other (or the other
//             says it does not have it), and closing the second window leaves the first alone with the whole cache.
//   links     ?data=<folder>/&ds=<id> (no datasets.json), and ?data=<folder>/ with a datasets.json
//             carrying tags: the picker lists every dataset, a tag filter narrows the list, a row opens.
//   outlines  a cells layer with segmentation boundaries (xenium_pub), Cell shapes on Outlines, zoomed in
//             until the outlines take over from the dots: tiles fetched and drawn, no error; compared with
//             Chromium by pixels.
//   phone     375 x 812, touch (WebKit as iPhone X, Chromium as a phone): tiles drawn, no sideways
//             scrolling, the tab rail there, Save image 1x; every Save image size offered only when the whole
//             picture (map, name and legend) fits the phone's limit, and the largest one saved at that size;
//             WebKit compared with Chromium by pixels. firefox-phone: the same 375 x 812 window without touch
//             (clicked, not tapped); the touch-only checks (the page sees a touch screen, the phone's Save image
//             limits) are left out there and the page is checked to use the desktop limits instead. The Compare button
//             is not shown on a phone (no linked windows there).
//
//   node tests/browsers_test.js                         (the public exports: SODA_TEST_DATA or --data-root, see tests/test_data.js)
//   node tests/browsers_test.js --only webkit,msedge --shots <dir> --data-root <url or folder> --ports 8892,8893 (default: free ports)
//        --ds visiumhd_pub,xenium_pub (only these datasets)  --skip probe,datasets,walk,compare,links,outlines
//   browsers: chromium, webkit, msedge (Edge; skipped when it is not installed), firefox (skipped when Playwright's
//             Firefox is not installed), webkit-phone, chromium-phone, firefox-phone
//
// Playwright is not a dependency of this repo (see tests/phone_browser_test.js): install it anywhere and
// point NODE_PATH or --playwright at its node_modules (PLAYWRIGHT_BROWSERS_PATH if its browsers are
// elsewhere). Without it, or without the public exports, or without deck.gl from the CDN (no internet), the
// test says SKIPPED and exits 0; a browser that cannot be launched is SKIPPED on its own (a line
// "SKIPPED <browser>: cannot launch (...)"), so a machine without Edge or Firefox still passes.
//
// Firefox reports two things in its console that Chromium and WebKit keep quiet about, and that are not
// faults of the page; only these two, and only in Firefox, are not counted as console errors:
//   - the vendor/ deck.gl tag when there is no vendor/ folder: Firefox checks the 404 page against the
//     tag's SHA-384 and says the hashes do not match (Chromium says "404", which is allowed for vendor/
//     already); libs.js then loads the CDN copy, as in every browser.
//   - deck.gl's bundled Long.js probing WebAssembly when it loads: the page's Content-Security-Policy has
//     no 'wasm-unsafe-eval' (on purpose: tests/smoke_test.py checks it is no wider), Firefox logs the
//     refusal, deck.gl catches it and uses plain JavaScript. Chromium refuses it the same way, silently.
//
// Playwright's WebKit is the WebKit engine, not Safari: on Windows it does not use Apple's image decoders,
// Metal, or an iPhone's memory limits (iOS caps a canvas at about 16.7 million pixels; the page keeps a
// phone's Save image under that). Passing here proves the page's code paths work in WebKit, not on a
// given iPhone or Mac.
"use strict";
const http = require("http");
const fs = require("fs");
const os = require("os");
const path = require("path");

const arg = (k, d) => { const i = process.argv.indexOf(k); return i > 0 ? process.argv[i + 1] : d; };
const REPO = path.resolve(__dirname, "..");
const D = require("./test_data");                                    // the public exports: SODA_TEST_DATA / --data-root
const PAGE = path.resolve(arg("--page", path.join(REPO, "web", "SpatialViewer")));
const ONLY = arg("--only", "");
const SHOTS = arg("--shots", "");
const ALL = process.argv.includes("--all");
const DS = arg("--ds", "");                                          // e.g. visiumhd_pub,xenium_pub: only these
const SKIP = arg("--skip", "").split(",").filter(Boolean);           // steps to leave out: probe, datasets, walk, compare, links, hardening, outlines
const PORTS = arg("--ports", "0,0").split(",").map(Number);         // data server, page server; 0 = any free port
const OPEN_MS = 120000;                                             // the public exports are large, and may come over the internet

let pw = null;
for (const base of [arg("--playwright", ""), ...(process.env.NODE_PATH || "").split(path.delimiter), path.join(REPO, "node_modules")].filter(Boolean)) {
  try { pw = require(path.join(base, "playwright")); break; } catch (e) { /* next */ }
}
if (!pw) { try { pw = require("playwright"); } catch (e) { /* none */ } }
if (!pw) { console.log("SKIPPED: playwright is not installed (see the top of this file)"); process.exit(0); }

// ---- the datasets: one public export of each kind (each kind has one, so --all is the same list) ----------
// GeoMx, Stereo-seq and CosMx are restricted (tests/test_data.js): without SODA_TEST_PRIVATE_DATA or a local copy, SKIPPED with a line
const allIds = D.PUBLIC.slice();
const KINDS = [["Visium HD", "visiumhd_pub"], ["multiplex protein", "phenocycler_pub"], ["Xenium", "xenium_pub"], ["GeoMx", "geomx_1"],
  ["Stereo-seq", "stereoseq_1"], ["CosMx", "cosmx_1"], ["MERSCOPE", "merscope_2"]];
let picked = (ALL ? allIds.map((id) => [id, id]) : KINDS)
  .filter(([, id]) => !DS || DS.split(",").includes(id));

// ---- two small static servers: the export (with CORS, as blob storage) and the page ----------------
// /plain/...  the export folder as it is;  /cat/...  the same, plus a datasets.json with tags (made here)
const catDir = fs.mkdtempSync(path.join(os.tmpdir(), "soda-browsers-"));
const CAT_IDS = picked.map((x) => x[1]).filter((id) => D.HOSTED.includes(id)).slice(0, 4);   // all on the hosted copy
const tagOf = (id) => ({ platform: /visium/.test(id) ? "Visium HD" : /phenocycler|multiplex/.test(id) ? "multiplex" : "other", tissue: "public" });
fs.writeFileSync(path.join(catDir, "datasets.json"), JSON.stringify({ soda: 1, title: "Browsers test", default: CAT_IDS[0],
  datasets: CAT_IDS.map((id) => ({ id, label: `Public ${id}`, group: "Public", tags: { ...tagOf(id), slide: `S${id.length}` } })) }));
const TYPES = { ".html": "text/html", ".js": "text/javascript", ".css": "text/css", ".json": "application/json", ".png": "image/png", ".jpg": "image/jpeg", ".bin": "application/octet-stream", ".geojson": "application/json" };
const inside = (root, rel) => { const f = path.join(root, rel); return f.startsWith(root + path.sep) || f === root ? f : null; };
// served as blob storage serves them (tests/static_files.js): CORS, no caching, HTTP Range (206) for gene packs
const { sendFile, preflight } = require("./static_files");
// a route gives a local file, or { rel } for a file of the data root (tests/test_data.js: a folder or a URL)
function handler(route) {
  return (req, res) => {
    if (req.method === "OPTIONS") { preflight(req, res); return; }
    const file = route(decodeURIComponent(req.url.split("?")[0]));
    if (file && file.rel) { D.serve(req, res, file.rel, TYPES[path.extname(file.rel)]); return; }
    if (!file || !fs.existsSync(file) || fs.statSync(file).isDirectory()) { res.writeHead(404, { "Access-Control-Allow-Origin": "*" }); res.end("not found"); return; }
    sendFile(req, res, file, TYPES[path.extname(file)]);
  };
}
const dataRoute = (u) => {
  const m = u.match(/^\/(plain|cat)(\/.*)$/);
  if (!m) return null;
  if (m[2] === "/datasets.json") return m[1] === "cat" ? path.join(catDir, "datasets.json") : null;   // /plain/: no datasets.json
  return { rel: m[2].slice(1) };
};
const pageRoute = (u) => inside(PAGE, u === "/" ? "index.html" : u);
function listen(h, port) {
  return new Promise((resolve, reject) => {
    const srv = http.createServer(h);
    srv.once("error", (e) => { if (e.code !== "EADDRINUSE" || !port) return reject(e); console.log(`  (port ${port} is taken: using a free one)`); listen(h, 0).then(resolve, reject); });
    srv.listen(port || 0, "127.0.0.1", () => resolve(srv));
  });
}

// ---- what the page holds -------------------------------------------------------------------------------
// lit: share of map pixels that are not black; thumb: the map averaged on a 64 x 40 grid (RGB), for the diff
const look = (page) => page.evaluate(() => {
  const sv = window.__spatialViewer, S = sv && sv.state, st = document.getElementById("status"), out = {};
  out.status = st && !st.hidden ? st.textContent : "";
  out.err = !!(st && !st.hidden && st.classList.contains("err"));
  out.deck = !!window.deck;
  if (!S) return out;
  out.id = S.m && S.m.id;
  out.cached = sv.diag ? sv.diag().cached : 0;
  out.decoder = sv.diag ? sv.diag().decoder : null;
  out.exprKey = S.exprKey || null;
  out.geneIndex = (S.geneIndex || []).length;
  out.genes = (S.genes || []).map((g) => ({ name: g.name, n: g.n }));
  out.view = S.umap ? S.umap.view : null;
  out.lit = -1;
  if (sv.deck && sv.deck.device && !sv.deck.device.gl.isContextLost()) {
    const gl = sv.deck.device.gl; sv.deck.redraw("test");
    const w = gl.drawingBufferWidth, h = gl.drawingBufferHeight, px = new Uint8Array(w * h * 4);
    gl.bindFramebuffer(gl.FRAMEBUFFER, null); gl.readPixels(0, 0, w, h, gl.RGBA, gl.UNSIGNED_BYTE, px);
    let lit = 0; for (let i = 0; i < px.length; i += 4) if (px[i] > 8 || px[i + 1] > 8 || px[i + 2] > 8) lit++;
    out.lit = lit / (w * h);
    const GX = 64, GY = 40, sum = new Float64Array(GX * GY * 3), cnt = new Float64Array(GX * GY);
    for (let y = 0; y < h; y++) {
      const gy = Math.min(GY - 1, Math.floor(((h - 1 - y) / h) * GY));         // readPixels starts at the bottom
      for (let x = 0; x < w; x++) {
        const g = gy * GX + Math.min(GX - 1, Math.floor((x / w) * GX)), i = 4 * (y * w + x);
        sum[3 * g] += px[i]; sum[3 * g + 1] += px[i + 1]; sum[3 * g + 2] += px[i + 2]; cnt[g]++;
      }
    }
    out.thumb = Array.from(sum, (v, i) => Math.round(v / Math.max(1, cnt[Math.floor(i / 3)])));
  }
  return out;
});
// every image tile on show has arrived and is drawn (as snapshot.js waitTiles: a few checks in a row)
const tilesDone = (page, ms = 60000) => page.evaluate(async (ms) => {
  const deck = window.__spatialViewer.deck, t0 = performance.now();
  let calm = 0;
  while (performance.now() - t0 < ms) {
    await new Promise((r) => setTimeout(r, 150));
    const tl = (deck.layerManager ? deck.layerManager.getLayers() : []).filter((l) => l.state && l.state.tileset);
    const all = tl.every((l) => l.isLoaded && (l.state.tileset.selectedTiles || []).every((t) => t.isLoaded));
    calm = all ? calm + 1 : 0;
    if (calm >= 4) return Math.round(performance.now() - t0);
  }
  return -1;
}, ms);
// look() again until the map has drawn something (a tile can be decoded a frame before it is drawn)
async function lookDrawn(page, min, ms = 15000) {
  const t0 = Date.now();
  let r = await look(page);
  while (!(r.lit > min) && Date.now() - t0 < ms) { await page.waitForTimeout(300); r = await look(page); }
  return r;
}
const until = (page, fn, a, ms = OPEN_MS) => page.waitForFunction(fn, a, { timeout: ms, polling: 200 }).then(() => true, () => false);
const settled = (page) => until(page, () => { const st = document.getElementById("status"); return !st || st.hidden || !/Loading|Opening/.test(st.textContent); }, null, 60000);
const pickFromSearch = (page, box, results, name) => page.evaluate(([box, results, name]) => {
  const el = document.getElementById(box);
  el.value = name; el.dispatchEvent(new Event("input", { bubbles: true }));
  const hit = document.querySelector(`#${results} > div`);
  if (!hit) return null;
  const shown = hit.querySelector("span") ? hit.querySelector("span").textContent : "";
  hit.dispatchEvent(new MouseEvent("mousedown", { bubbles: true, cancelable: true }));
  return shown;
}, [box, results, name]);
const diffThumb = (a, b) => { if (!a || !b || a.length !== b.length) return null; let s = 0; for (let i = 0; i < a.length; i++) s += Math.abs(a[i] - b[i]); return s / a.length; };
const pngSize = (buf) => (buf && buf.length > 24 && buf.readUInt32BE(12) === 0x49484452 ? { w: buf.readUInt32BE(16), h: buf.readUInt32BE(20) } : null);
// the layer's gene with the most counts that has data, and its points: from its own file, or from its slice of a
// gene pack (genepack.py: genes.json gives each gene's byte count `b`), as tests/compat_test.js reads them
const geneOf = async (id, key) => {
  const L = ((await D.manifest(id)).expr || []).find((e) => e.key === key) || {};
  const genes = ((await D.json(`${id}/expr/${key}/genes.json`)) || []).slice().sort((a, b) => b.t - a.t);
  for (const g of L.genePack ? genes : genes.slice(0, 20)) {               // loose files: asked for one at a time
    const b = L.genePack ? (g.b > 0 ? g.b : -1) : await D.size(`${id}/expr/${key}/gene/${encodeURIComponent(g.n)}.bin`);
    if (b > 0) return { name: g.n, pts: b / 12, packed: !!L.genePack };
  }
  return null;
};

// A form control as the browser drew it: its background (the most common luminance) must be dark, as on the
// rest of the page, and its text must stand out from it (WCAG contrast ratio of the pixel furthest from the
// background). Catches a browser that draws a <select> natively, ignoring the page's colours.
async function controlContrast(page, sel) {
  const el = await page.$(sel);
  if (!el || !(await el.isVisible())) return null;
  const b64 = (await el.screenshot()).toString("base64");
  return page.evaluate(async (b64) => {
    const bmp = await createImageBitmap(await (await fetch(`data:image/png;base64,${b64}`)).blob());
    const c = document.createElement("canvas"); c.width = bmp.width; c.height = bmp.height;
    const g = c.getContext("2d"); g.drawImage(bmp, 0, 0);
    const d = g.getImageData(0, 0, c.width, c.height).data, lin = (v) => { v /= 255; return v <= 0.04045 ? v / 12.92 : Math.pow((v + 0.055) / 1.055, 2.4); };
    const hist = new Map(), L = [];
    for (let i = 0; i < d.length; i += 4) { const l = 0.2126 * lin(d[i]) + 0.7152 * lin(d[i + 1]) + 0.0722 * lin(d[i + 2]); L.push(l); const k = Math.round(l * 200); hist.set(k, (hist.get(k) || 0) + 1); }
    const bg = [...hist.entries()].sort((a, b) => b[1] - a[1])[0][0] / 200;
    let fg = bg; for (const l of L) if (Math.abs(l - bg) > Math.abs(fg - bg)) fg = l;
    return { bg: +bg.toFixed(3), fg: +fg.toFixed(3), ratio: +((Math.max(bg, fg) + 0.05) / (Math.min(bg, fg) + 0.05)).toFixed(2) };
  }, b64);
}

// ---- results ------------------------------------------------------------------------------------------
let pass = 0, fail = 0;
const failures = [], table = {};                  // table[feature][browser] = "works" | "BROKEN" | "skipped"
const mark = (feature, eng, cond) => { const row = (table[feature] = table[feature] || {}); if (row[eng] !== "BROKEN") row[eng] = cond ? "works" : "BROKEN"; };
const info = (feature, eng, yes) => { (table[feature] = table[feature] || {})[eng] = yes ? "yes" : "no"; };   // offered or not: not a check
const thumbs = {};                                // thumbs[`${group}|${key}`][engine]

// files the page asks for that an export need not have (see tests/compat_test.js)
// vendor/: deck.gl and umap-js are asked for there first and come from the CDN when a copy has no vendor/ (libs.js)
const OPTIONAL = [/\/gates\.json$/, /\/markers\/[\w.-]+\.json$/, /\/about\.html$/, /\/vendor\//];
// Firefox's console lines that are not the page's faults (see the top of this file): the vendor/ tag's hash
// checked against a 404, and deck.gl's WebAssembly probe refused by the page's Content-Security-Policy
const FIREFOX_QUIET = [/integrity attribute match the content of the subresource at .[^ ]*\/vendor\//,
  /Content-Security-Policy: .*blocked WebAssembly .*Missing 'wasm-unsafe-eval'/];

async function openPage(browser, eng, url, label, feature, before) {
  // serviceWorkers "block": an old service worker of the removed offline use must not answer requests out of
  // page.route's sight
  const ctx = await browser.newContext({ ...(eng.device || { viewport: { width: 1280, height: 800 } }), acceptDownloads: true, serviceWorkers: "block" });
  if (eng.clipboard) await ctx.grantPermissions(["clipboard-read", "clipboard-write"]).catch(() => {});
  const page = await ctx.newPage();
  const log = { console: [], pageErrors: [], bad: [], ok: [] };
  page.on("console", (m) => { if (m.type() === "error") log.console.push({ text: m.text(), url: (m.location() || {}).url || "" }); });
  page.on("pageerror", (e) => log.pageErrors.push(String(e.message || e).split("\n")[0]));
  page.on("requestfailed", (r) => { const f = r.failure() ? r.failure().errorText : "failed"; if (!/ERR_ABORTED|cancelled|NS_BINDING_ABORTED/i.test(f)) log.bad.push({ url: r.url(), status: f }); });
  page.on("response", (r) => (r.status() >= 400 ? log.bad : log.ok).push({ url: r.url(), status: r.status(), range: (r.request().headers() || {}).range || "" }));
  page.on("dialog", (d) => d.accept().catch(() => {}));         // a confirm() (e.g. a large volcano on touch) is answered yes
  await page.route("**/config.js*", (rt) => rt.fulfill({ status: 200, contentType: "text/javascript", body: "window.SpatialViewerConfig = {};" }));
  await page.addInitScript(() => {          // how many WebGL contexts the page makes (a browser keeps about 16; the oldest is lost)
    const real = HTMLCanvasElement.prototype.getContext;
    HTMLCanvasElement.prototype.getContext = function (type, ...a) {
      const had = this.__glMade, c = real.call(this, type, ...a);
      if (c && /^(webgl|experimental-webgl)/.test(type) && !had) { this.__glMade = true; window.__glContexts = (window.__glContexts || 0) + 1; }
      return c;
    };
  });
  if (before) await before(page);                                 // routes of a failure check, set before the page loads
  await page.goto(url, { waitUntil: "load" });
  const ok = (cond, what, feat) => {
    if (cond) pass++; else { fail++; failures.push(`${eng.name} | ${label} | ${what}`); }
    if (feat || feature) mark(feat || feature, eng.name, cond);
    console.log(`${cond ? "  ok  " : "  FAIL"} ${eng.name} | ${label} | ${what}`);
    return cond;
  };
  const shot = async (name) => {
    if (!SHOTS) return;
    const dir = path.join(SHOTS, eng.name); fs.mkdirSync(dir, { recursive: true });
    await page.screenshot({ path: path.join(dir, `${name.replace(/[^\w]+/g, "_")}.png`) }).catch(() => {});
  };
  return { ctx, page, log, ok, shot, eng };
}
// A connection to this test's own servers that the machine dropped (refused / reset at the TCP level, seen
// now and then on a busy Windows machine with other servers running) says nothing about the browser: up to
// two per page are noted, not failed. Any HTTP error, a third drop, or a drop elsewhere still fails.
const DROPPED = /ERR_CONNECTION_(FAILED|RESET|REFUSED|ABORTED)|ERR_EMPTY_RESPONSE|Could not connect|connection was (reset|refused)/i;
let LOCAL = [];
function checkClean(t, feature, allowed = []) {
  const drop = (b) => DROPPED.test(b.status) && LOCAL.some((o) => b.url.startsWith(o));
  const dropped = t.log.bad.filter(drop);
  const flaky = dropped.length && dropped.length <= 2 ? new Set(dropped.map((b) => b.url)) : new Set();
  if (flaky.size) console.log(`        note: ${dropped.length} connection(s) to the test servers dropped, not counted: ${[...flaky].map((u) => u.replace(/^https?:\/\/[^/]+/, "")).join(", ")}`);
  const allow = (u) => flaky.has(u) || allowed.concat(OPTIONAL).some((re) => re.test(u));
  const bad = t.log.bad.filter((b) => !allow(b.url));
  const cons = t.log.console.filter((c) => !(/Failed to load resource|status of 404/.test(c.text) && (!c.url || allow(c.url)))
    && !(t.eng.type === pw.firefox && FIREFOX_QUIET.some((re) => re.test(c.text))));
  t.ok(!t.log.pageErrors.length, `no page error${t.log.pageErrors.length ? ": " + t.log.pageErrors.slice(0, 3).join(" / ") : ""}`, feature);
  t.ok(!cons.length, `no console error${cons.length ? ": " + cons.slice(0, 3).map((c) => `${c.text} ${c.url}`).join(" / ") : ""}`, feature);
  t.ok(!bad.length, `no failed request${bad.length ? ": " + bad.slice(0, 5).map((b) => `${b.status} ${b.url.replace(/^https?:\/\/[^/]+/, "")}`).join(", ") : ""}`, feature);
}
// a click that starts a file download; returns { name, buf } or null
async function download(t, click, ms = 120000) {
  const wait = t.page.waitForEvent("download", { timeout: ms }).catch(() => null);
  await click();
  const d = await wait;
  if (!d) return null;
  const p = await d.path().catch(() => null);
  return { name: d.suggestedFilename(), buf: p && fs.existsSync(p) ? fs.readFileSync(p) : null };
}

// ---- probe: what this browser offers ----------------------------------------------------------------
async function probe(browser, eng, urls) {
  const t = await openPage(browser, eng, `${urls.page}/index.html`, "probe", null);
  try {
    const p = await t.page.evaluate(async () => {
      const out = {};
      const c = document.createElement("canvas"), gl = c.getContext("webgl2");
      out.webgl2 = !!gl;
      if (gl) {
        out.maxTexture = gl.getParameter(gl.MAX_TEXTURE_SIZE);
        out.maxRenderbuffer = gl.getParameter(gl.MAX_RENDERBUFFER_SIZE);
        out.maxViewport = Array.from(gl.getParameter(gl.MAX_VIEWPORT_DIMS) || []);
        out.extensions = gl.getSupportedExtensions();
        const dbg = gl.getExtension("WEBGL_debug_renderer_info");
        out.renderer = dbg ? gl.getParameter(dbg.UNMASKED_RENDERER_WEBGL) : gl.getParameter(gl.RENDERER);
      }
      out.createImageBitmap = typeof createImageBitmap === "function";
      if (out.createImageBitmap) {
        const blob = await new Promise((r) => { const k = document.createElement("canvas"); k.width = k.height = 4; k.getContext("2d").fillRect(0, 0, 2, 2); k.toBlob(r, "image/png"); });
        try { await createImageBitmap(blob, { premultiplyAlpha: "none", colorSpaceConversion: "none" }); out.createImageBitmapOptions = true; } catch (e) { out.createImageBitmapOptions = String(e.message || e); }
      }
      out.offscreenCanvas = typeof OffscreenCanvas === "function";
      out.offscreenWebgl2 = false;
      try { out.offscreenWebgl2 = !!new OffscreenCanvas(4, 4).getContext("webgl2"); } catch (e) { /* no */ }
      out.worker = typeof Worker === "function";
      out.moduleWorker = await new Promise((res) => {
        try {
          const u = URL.createObjectURL(new Blob(["export const x = 1; postMessage(x);"], { type: "text/javascript" }));
          const w = new Worker(u, { type: "module" }); const tm = setTimeout(() => { w.terminate(); res(false); }, 3000);
          w.onmessage = () => { clearTimeout(tm); w.terminate(); res(true); }; w.onerror = () => { clearTimeout(tm); res(false); };
        } catch (e) { res(false); }
      });
      out.bigInt64Array = typeof BigInt64Array === "function";
      out.arrayAt = typeof [].at === "function" && typeof new Float32Array(1).at === "function";
      out.structuredClone = typeof structuredClone === "function";
      out.roundRect = typeof CanvasRenderingContext2D.prototype.roundRect === "function";
      out.css = { dvh: CSS.supports("height", "100dvh"), has: CSS.supports("selector(:has(a))"), backdropFilter: CSS.supports("backdrop-filter", "blur(2px)"),
        webkitBackdropFilter: CSS.supports("-webkit-backdrop-filter", "blur(2px)"), inset: CSS.supports("inset", "0"), gap: CSS.supports("gap", "1px"), aspectRatio: CSS.supports("aspect-ratio", "1") };
      out.clipboardWrite = !!(navigator.clipboard && navigator.clipboard.writeText);
      out.anchorDownload = "download" in document.createElement("a");
      out.compressionStream = typeof CompressionStream === "function";
      // the largest 2D canvas that still encodes a PNG (Save image draws on one)
      out.canvas = {};
      for (const [w, h] of [[8192, 8192], [16384, 4096], [11000, 7000], [16384, 16384]]) {
        const k = document.createElement("canvas"); k.width = w; k.height = h;
        const ctx = k.getContext("2d"); let res = "no context";
        if (ctx) {
          ctx.fillStyle = "#f00"; ctx.fillRect(w - 2, h - 2, 2, 2);
          const px = ctx.getImageData(w - 1, h - 1, 1, 1).data[0];
          const b = await new Promise((r) => { try { k.toBlob(r, "image/png"); } catch (e) { r(null); } });
          res = px === 255 && b && b.size > 0 ? "ok" : `drawn ${px}, blob ${b ? b.size : "null"}`;
        }
        out.canvas[`${w}x${h}`] = res; k.width = k.height = 0;
      }
      out.userAgent = navigator.userAgent;
      return out;
    });
    if (SHOTS) { fs.mkdirSync(path.join(SHOTS, eng.name), { recursive: true }); fs.writeFileSync(path.join(SHOTS, eng.name, "probe.json"), JSON.stringify(p, null, 1)); }
    console.log(`  ${eng.name}: WebGL2 ${p.webgl2} (texture ${p.maxTexture}, renderbuffer ${p.maxRenderbuffer}, ${p.renderer}); createImageBitmap ${p.createImageBitmap} (options ${p.createImageBitmapOptions}); ` +
      `OffscreenCanvas ${p.offscreenCanvas}; module worker ${p.moduleWorker}; BigInt64Array ${p.bigInt64Array}; .at ${p.arrayAt}; CSS ${JSON.stringify(p.css)}; clipboard ${p.clipboardWrite}; canvas ${JSON.stringify(p.canvas)}`);
    t.ok(p.webgl2 && p.maxTexture >= 4096, `WebGL2 with textures of 4096 px or more (${p.maxTexture})`, "WebGL2");
    t.ok(p.worker && p.arrayAt, "Web Workers and Array / TypedArray .at() (plots_worker.js)", "Workers, .at()");
    t.ok(p.anchorDownload, "<a download> (Save JSON / CSV / PNG)", "file download via <a download>");
    t.ok(p.canvas["8192x8192"] === "ok", `a 2D canvas of 8192 x 8192 encodes a PNG (${p.canvas["8192x8192"]})`, "2D canvas 8192 x 8192 to PNG");
    mark("createImageBitmap with options", eng.name, p.createImageBitmapOptions === true);
    info("CompressionStream (shorter view links; without it the link is plain JSON)", eng.name, p.compressionStream);
    for (const [k, v] of Object.entries(p.css)) if (k !== "webkitBackdropFilter") mark(`CSS ${k}`, eng.name, v);
    info("OffscreenCanvas (not used by the page)", eng.name, p.offscreenCanvas);
    info("module workers (not used by the page)", eng.name, p.moduleWorker);
    return p;
  } catch (e) { t.ok(false, `probe crashed: ${String(e.message).split("\n")[0]}`, "WebGL2"); return null; }
  finally { await t.ctx.close(); }
}

// ---- one dataset: open, tiles, a gene ----------------------------------------------------------------
async function checkDataset(browser, eng, urls, kind, id) {
  const m = await D.manifest(id);
  const feature = `open ${kind}`;
  const t = await openPage(browser, eng, `${urls.page}/index.html?static=${urls.data}/plain&ds=${id}`, `${id}`, feature);
  try {
    const opened = await until(t.page, (id) => { const sv = window.__spatialViewer; return sv && sv.state && sv.state.m && sv.state.m.id === id; }, id);
    let r = await look(t.page);
    if (!t.ok(r.deck && opened, `opened (${m.kind}${r.status ? `; status "${r.status}"` : ""})`)) return;
    const hasImage = !!(m.image && m.image.channels && m.image.channels.length);
    if (hasImage) {
      const tiles = await until(t.page, () => window.__spatialViewer.diag().cached > 0);
      await settled(t.page);
      const ms = await tilesDone(t.page);
      r = await lookDrawn(t.page, 0.005);
      t.ok(tiles && ms >= 0, `image tiles fetched, decoded and all on show (${r.cached} kept, decoder: ${r.decoder}; ${ms >= 0 ? `settled in ${ms} ms` : "NOT all loaded after 60 s"})`, `tiles decode and draw`);
    } else { await settled(t.page); r = await lookDrawn(t.page, 0.005); }
    t.ok(r.lit > 0.005, `slide drawn (lit ${r.lit.toFixed(3)})`, hasImage ? "tiles decode and draw" : null);
    // the fit view, before anything is added: compared with Chromium
    const k = `${eng.group}|${id}`; (thumbs[k] = thumbs[k] || {})[eng.name] = r.thumb;
    await t.shot(`${id}_1_open`);
    if ((m.expr || []).length) {
      await until(t.page, () => { const S = window.__spatialViewer.state; return S.exprKey && S.geneIndex && S.geneIndex.length; });
      await settled(t.page);
      r = await look(t.page);
      const g = r.exprKey && (await geneOf(id, r.exprKey));
      if (g) {
        const shown = await pickFromSearch(t.page, "geneSearch", "geneResults", g.name);
        await until(t.page, () => window.__spatialViewer.state.genes.length > 0, null, 60000);
        await settled(t.page);
        await t.page.waitForTimeout(500);
        r = await look(t.page);
        const got = r.genes[0];
        t.ok(shown === g.name && got && got.name === g.name && got.n === g.pts, `gene ${g.name} added from the search: ${got ? got.n : "none"} points (file ${g.pts})`, "add a gene");
        if (g.packed) {
          const ranged = t.log.ok.filter((q) => /\/genepack\/[^/]+\.pack$/.test(q.url.split("?")[0]));
          t.ok(ranged.length >= 1 && ranged.every((q) => q.status === 206 && /^bytes=\d+-\d+$/.test(q.range)),
            `packed gene read with HTTP Range: ${ranged.map((q) => `${q.status} ${q.range}`).join(", ") || "no request into a pack"}`, "packed gene (HTTP Range, 206)");
        }
        await t.shot(`${id}_2_gene`);
      }
    }
    t.ok(!r.err, `no error line (${r.status || "none"})`);
    checkClean(t, feature);
  } catch (e) { t.ok(false, `crashed: ${String(e.message).split("\n")[0]}`); }
  finally { await t.ctx.close(); }
}

// ---- the walk through the features, on visiumhd_pub ---------------------------------------------------
async function walk(browser, eng, urls) {
  const id = "visiumhd_pub";
  const t = await openPage(browser, eng, `${urls.page}/index.html?static=${urls.data}/plain&ds=${id}`, "walk", null);
  const P = t.page, ev = (f, a) => P.evaluate(f, a);
  const step = async (feature, fn) => {
    try { await fn(); } catch (e) { t.ok(false, `${feature} crashed: ${String(e.message).split("\n")[0]}`, feature); }
    const gl = await ev(() => ({ lost: window.__spatialViewer.deck.device.gl.isContextLost(), contexts: window.__glContexts || null }));
    t.ok(!gl.lost, `after ${feature}: the map's graphics context is still there${gl.contexts ? ` (${gl.contexts} WebGL contexts made so far)` : ""}`, "map keeps its graphics context");
  };
  try {
    await until(P, () => { const sv = window.__spatialViewer; const S = sv && sv.state; return S && S.m && S.exprKey && S.lasso && S.lasso.pos && S.geneIndex && S.geneIndex.length; }, null, 60000);
    await until(P, () => window.__spatialViewer.diag().cached > 0);
    await settled(P);
    const key = await ev(() => window.__spatialViewer.state.exprKey);
    const g = await geneOf(id, key);
    await pickFromSearch(P, "geneSearch", "geneResults", g.name);
    await until(P, () => window.__spatialViewer.state.genes.length > 0, null, 60000);

    // the page's dark form controls stay dark and readable (Layer drop-down, gene search box)
    await step("form controls readable (dark <select>)", async () => {
      // every section starts closed: open the genes' section to see its drop-downs
      await ev(() => { const s = document.getElementById("panelGenes"); if (s && s.classList.contains("closed")) document.getElementById("secBtn-panelGenes").click(); });
      await P.waitForTimeout(200);
      for (const sel of ["#exprSelect", "#geneBlend"]) {
        const c = await controlContrast(P, sel);
        t.ok(c && c.bg < 0.2 && c.ratio >= 4, `${sel} drawn dark and readable: ${JSON.stringify(c)} (background luminance < 0.2, contrast >= 4)`, "form controls readable (dark <select>)");
      }
      await t.shot("walk_controls");
    });

    // Tissue / UMAP / Side by side
    await step("Tissue / UMAP / Side by side", async () => {
      for (const v of ["umap", "sbs", "tissue"]) {
        await ev((v) => window.__spatialViewer.umapSetView(v), v);
        // asked once: a view asked for during a flight is kept (U.wantView) and applied when the flight ends, so wait
        // for the flight, and the one it queued, to end on this view (a slow runner takes a while to load the UMAP)
        const ok = await until(P, (v) => { const U = window.__spatialViewer.state.umap;
          return U && U.view === v && !U.anim && !U.wantView; }, v, 60000);
        await settled(P); await P.waitForTimeout(600);
        const r = await look(P);
        const why = ok ? null : await ev(() => { const U = window.__spatialViewer.state.umap, s = document.getElementById("status");
          return { view: U && U.view, anim: U && U.anim, status: s && !s.hidden ? s.textContent.slice(0, 120) : "" }; });
        t.ok(ok && r.lit > 0.002, `view ${v} (lit ${r.lit.toFixed(3)})${why ? " " + JSON.stringify(why) : ""}`, "Tissue / UMAP / Side by side");
        await t.shot(`walk_view_${v}`);
      }
    });

    // Draw a region: two circles with the mouse, then Run UMAP (worker + umap-js from the CDN)
    await step("Draw a region (mouse) + Run UMAP (worker)", async () => {
      await ev(() => window.__spatialViewer.fitView && window.__spatialViewer.fitView());
      await P.waitForTimeout(400);
      const box = await ev(() => { const r = document.getElementById("viewport").getBoundingClientRect(); return { x: r.left, y: r.top, w: r.width, h: r.height }; });
      // centred on two bins (a third and two thirds down the bin list), so each circle is over tissue whatever its shape
      const at = await ev(() => { const sv = window.__spatialViewer, pos = sv.state.lasso.pos, n = pos.length / 2, vp = sv.deck.getViewports()[0];
        return [0.37, 0.63].map((f) => { const r = Math.floor(n * f), p = vp.project([pos[2 * r], pos[2 * r + 1]]); return [p[0], p[1]]; }); });
      for (const [px, py] of at) {
        await ev(() => { if (window.__spatialViewer.state.lasso.tool !== "circle") document.getElementById("lassoCircle").click(); });
        // small enough that the two circles stay under Draw a region's fast-mode size (LASSO.fastAbove, 20,000 bins) on
        // visiumhd_pub's 8 um layer (about 15,000 bins): the walk runs the exact path, as it did on the demo data
        // (tests/drawregion_test.js runs fast mode)
        const cx = box.x + px, cy = box.y + py, rad = Math.min(box.w, box.h) * 0.025;
        await P.mouse.move(cx, cy); await P.mouse.down();
        for (let i = 1; i <= 6; i++) await P.mouse.move(cx + (rad * i) / 6, cy);
        await P.mouse.up();
        await P.waitForTimeout(250);
      }
      await ev(() => { const L = window.__spatialViewer.state.lasso; if (L.tool) document.dispatchEvent(new KeyboardEvent("keydown", { key: "Escape" })); });
      const shapes = await ev(() => window.__spatialViewer.state.lasso.shapes.map((s) => ({ kind: s.kind, bins: s.bins })));
      t.ok(shapes.length === 2 && shapes.every((s) => s.kind === "circle" && s.bins > 20), `two circles drawn with the mouse: ${JSON.stringify(shapes)}`, "Draw a region (mouse)");
      await t.shot("walk_circles");
      const t0 = Date.now();
      await ev(() => document.getElementById("lassoRun").click());
      const done = await until(P, () => { const L = window.__spatialViewer.state.lasso; return L.result && !L.running; }, null, 600000);   // a real section: Firefox on software rendering takes minutes
      const R = await ev(() => { const L = window.__spatialViewer.state.lasso, R = L.result; return R ? { n: R.info && R.info.n, emb: R.emb ? R.emb.length : R.xy ? R.xy.length : null, k: R.info && R.info.nClusters, msg: document.getElementById("lassoMsg").textContent } : { msg: document.getElementById("lassoMsg").textContent }; });
      t.ok(done && R.emb > 0, `Run UMAP finished in ${((Date.now() - t0) / 1000).toFixed(1)} s: ${JSON.stringify(R).slice(0, 200)}`, "Run UMAP (worker, umap-js)");
      await P.waitForTimeout(600);
      await t.shot("walk_lasso_umap");
      // Save JSON / Save CSV: <a download> of a Blob
      const j = await download(t, () => ev(() => document.getElementById("lassoSaveJson").click()), 20000);
      let parsed = null; try { parsed = j && j.buf ? JSON.parse(j.buf.toString("utf-8")) : null; } catch (e) { /* not JSON */ }
      t.ok(parsed && Array.isArray(parsed.shapes) && parsed.shapes.length === 2, `Save JSON downloaded ${j ? j.name : "nothing"} (${parsed ? parsed.shapes.length : 0} shapes)`, "Save JSON / CSV (download)");
      const c = await download(t, () => ev(() => document.getElementById("lassoSaveCsv").click()), 20000);
      const lines = c && c.buf ? c.buf.toString("utf-8").trim().split("\n").length : 0;
      t.ok(lines > 20, `Save CSV downloaded ${c ? c.name : "nothing"} (${lines} lines)`, "Save JSON / CSV (download)");
    });

    // volcano: two groups, its worker, the plot and its PNG (an SVG drawn on a canvas)
    await step("volcano (worker) + PNG", async () => {
      const groups = await ev(() => {
        const a = document.getElementById("vcA"); a.dispatchEvent(new Event("pointerdown"));
        const opts = [...a.options].filter((o) => !o.disabled && o.value).map((o) => o.value);
        if (opts.length < 2) return opts;
        a.value = opts[0]; a.dispatchEvent(new Event("change"));
        const b = document.getElementById("vcB"); const ob = [...b.options].filter((o) => !o.disabled && o.value && o.value !== opts[0]).map((o) => o.value);
        if (ob.length) { b.value = ob[0]; b.dispatchEvent(new Event("change")); }
        return [a.value, b.value];
      });
      if (!t.ok(groups.length === 2 && groups[0] !== groups[1], `volcano groups: ${groups.join(" vs ")}`, "volcano (worker)")) return;
      await ev(() => { window.SodaVolcano.result = null; document.getElementById("vcRun").click(); });
      const ok = await until(P, () => window.SodaVolcano.result || /could not|needs|failed|error/i.test(document.getElementById("vcMsg").textContent), null, 180000);
      const R = await ev(() => { const R = window.SodaVolcano.result; return { ok: !!R, genes: R ? R.genes.length : 0, tested: R ? R.info.nTested : 0, msg: document.getElementById("vcMsg").textContent }; });
      t.ok(ok && R.ok && R.tested > 10, `volcano: ${R.tested} of ${R.genes} genes tested ("${R.msg.slice(0, 120)}")`, "volcano (worker)");
      await P.waitForTimeout(500);
      await t.shot("walk_volcano");
      const png = await download(t, () => ev(() => document.getElementById("vcPng").click()), 20000);
      const z = png && pngSize(png.buf);
      t.ok(!!z && z.w > 200, `volcano Save PNG: ${png ? png.name : "nothing"} ${z ? `${z.w} x ${z.h}` : ""}`, "volcano / plot PNG (SVG on a canvas)");
      await ev(() => { const c = document.getElementById("vcClose"); if (c) c.click(); });
    });

    // heatmap / violin
    await step("heatmap / violin", async () => {
      const grp = await ev(() => { const s = document.getElementById("plGroup"); s.dispatchEvent(new Event("pointerdown")); const o = [...s.options].find((o) => !o.disabled && o.value); if (o) { s.value = o.value; s.dispatchEvent(new Event("change")); } return s.value; });
      await ev(() => { window.SodaPlots.result = null; document.getElementById("plTop").value = 3; document.getElementById("plGenes").value = ""; document.getElementById("plRun").click(); });
      const ok = await until(P, () => window.SodaPlots.result && document.querySelector("#plSvgBox svg") && !document.getElementById("plRun").disabled, null, 180000);
      const R = await ev(() => { const R = window.SodaPlots.result, s = document.querySelector("#plSvgBox svg"); return { genes: R ? R.genes.length : 0, groups: R ? R.cols.length : 0, rects: s ? s.querySelectorAll("rect").length : 0, msg: document.getElementById("plMsg").textContent }; });
      t.ok(ok && R.genes > 0 && R.rects > 0, `heatmap of ${R.genes} genes x ${R.groups} groups by ${grp} (${R.rects} rects)`, "heatmap / violin");
      await t.shot("walk_heatmap");
      const vio = await ev(async () => {
        const b = [...document.getElementById("plTabs").querySelectorAll("button")].find((x) => /violin/i.test(x.dataset.v || x.textContent));
        if (!b) return "no violin tab"; b.click();
        const t0 = Date.now(); while (Date.now() - t0 < 60000 && (document.getElementById("plRun").disabled || !document.querySelector("#plSvgBox svg path"))) await new Promise((r) => setTimeout(r, 100));
        return document.querySelectorAll("#plSvgBox svg path").length;
      });
      t.ok(typeof vio === "number" && vio > 0, `violin drawn (${vio} paths)`, "heatmap / violin");
      await t.shot("walk_violin");
      const btn = await ev(() => { const b = [...document.querySelectorAll("#plotsPanel button")].find((x) => /png/i.test(x.textContent)); return b ? b.id || b.textContent : null; });
      if (btn) {
        const png = await download(t, () => ev(() => [...document.querySelectorAll("#plotsPanel button")].find((x) => /png/i.test(x.textContent)).click()), 20000);
        const z = png && pngSize(png.buf);
        t.ok(!!z && z.w > 100, `plot Save PNG: ${png ? png.name : "nothing"} ${z ? `${z.w} x ${z.h}` : ""}`, "volcano / plot PNG (SVG on a canvas)");
      }
    });

    // the QC filter: total counts >= the median removes about half the bins
    await step("QC filter", async () => {
      const r = await ev(() => {
        const q = window.__spatialViewer.binQc; if (!q || !q.available()) return { none: true };
        const Q = q._debug(), m = Q.metrics[0], v = Array.from(m.v).sort((a, b) => a - b), med = v[Math.floor(v.length / 2)];
        let want = 0; for (const x of m.v) if (x < med) want++;
        q.setState({ thresholds: { [m.key]: { min: med, max: null } } });
        return { key: m.key, med, want, removed: q.removed(), url: new URL(location.href).searchParams.get("qc") };
      });
      t.ok(!r.none && r.removed === r.want && r.removed > 0 && /=/.test(r.url || ""), `QC ${r.key} >= ${r.med}: ${r.removed} removed (expected ${r.want}); ?qc=${r.url}`, "QC filter");
      await P.waitForTimeout(500);
      await t.shot("walk_qc");
      await ev(() => window.__spatialViewer.binQc.setState(null));
    });

    // Save image 1x / 2x / 4x: the PNG's own size
    await step("Save image PNG 1x / 2x / 4x", async () => {
      // Headless Chromium draws with SwiftShader (a software GPU): a 4x picture there can exhaust its memory and
      // the map's context is lost (seen about one run in two). 4x is checked on the browsers with a real GPU.
      // With visiumhd_pub's 516,880 bins a 2x picture lost it too (19 Nov 2026): on SwiftShader only 1x is saved.
      const soft = await ev(() => { const gl = window.__spatialViewer.deck.device.gl, d = gl.getExtension("WEBGL_debug_renderer_info"); return /SwiftShader|llvmpipe|softpipe|Software/i.test(d ? gl.getParameter(d.UNMASKED_RENDERER_WEBGL) : gl.getParameter(gl.RENDERER)); });
      if (soft) console.log(`        ${eng.name}: software renderer (SwiftShader): Save image 2x and 4x not run here`);
      // the map's size at rest, read once before any save: after a 2x save deck puts the canvas back to its own size
      // on its next frame, which can come late under load (WebKit), so reading it again later can see the 2x canvas
      const base = await ev(() => { const c = document.getElementById("viewport").querySelector("canvas"); return { w: c.width, h: c.height, cw: c.clientWidth }; });
      for (const s of soft ? [1] : [1, 2, 4]) {
        const want = await ev(({ s, base }) => { const sel = document.getElementById("snapScale"); const o = sel && [...sel.options].find((o) => o.value === String(s)); return { w: Math.round(base.w * s), h: Math.round(base.h * s), k: (base.w * s) / base.cw, off: !o || o.disabled, txt: o ? o.textContent : "" }; }, { s, base });
        if (want.off) { t.ok(true, `Save image ${s}x not offered here: ${want.txt}`); continue; }
        let said = "";
        const png = await download(t, async () => { said = await ev((s) => window.SodaSnapshot._save({ scale: s, bg: "black" }), s); }, 120000);
        const z = png && pngSize(png.buf);
        // the PNG: the map, a 30 css px line with the name above it, the legend beside it (k = PNG px per css px)
        t.ok(!!z && z.w >= want.w && Math.abs(z.h - want.h) <= Math.round(30 * want.k) + 2 && !/could not/.test(said),
          `Save image ${s}x: ${png ? png.name : "nothing"} ${z ? `${z.w} x ${z.h}` : ""} (map ${want.w} x ${want.h}) "${said.slice(0, 120)}"`, "Save image PNG 1x / 2x / 4x");
        if (s === 2 && z && SHOTS) { fs.mkdirSync(path.join(SHOTS, eng.name), { recursive: true }); fs.writeFileSync(path.join(SHOTS, eng.name, "walk_saveimage_2x.png"), png.buf); }
        const lost = await ev(() => window.__spatialViewer.deck.device.gl.isContextLost());
        t.ok(!lost, `after Save image ${s}x the map is still drawn (graphics context not lost)`, "Save image PNG 1x / 2x / 4x");
      }
    });

    // Copy link: the clipboard API, or the textarea fallback; the message says which worked
    // Copy link to this view (viewlink.js): the clipboard (or the box it shows when the browser will not copy),
    // then the link opened in a new tab of the same browser: the same dataset, layer, genes, mode and camera
    await step("Copy link to this view", async () => {
      await ev(() => { const sv = window.__spatialViewer, S = sv.state, i = S.m.image;
        sv.deck.props.onViewStateChange({ viewState: { ...S.viewState, target: [i.width * 0.4, i.height * 0.55, 0], zoom: S.viewState.zoom + 1 }, viewId: "ortho" }); });
      await settled(P); await P.waitForTimeout(400);
      const want = await ev(() => { const S = window.__spatialViewer.state; return { ds: S.dsId, layer: S.exprKey, genes: S.genes.map((g) => g.name), mode: S.umap ? S.umap.view : null,
        x: S.viewState.target[0], y: S.viewState.target[1], zoom: S.viewState.zoom }; });
      await ev(() => document.getElementById("btnCopyLink").click());
      await until(P, () => /link to this view|could not/i.test(document.getElementById("saveMsg").textContent), null, 15000);
      const said = await ev(() => document.getElementById("saveMsg").textContent);
      let clip = null;
      if (eng.clipboard) clip = await ev(() => navigator.clipboard.readText()).catch(() => null);
      const box = await ev(() => { const b = document.getElementById("saveLinkBox"); return b && !b.hidden ? b.value : null; });
      const link = clip || box || (await ev(async () => (await window.SodaViewLink.make()).url));
      t.ok(/^Link to this view copied/.test(said) && /#v=1[dj]\./.test(link) && (clip === null || clip === link),
        `Copy link to this view: "${said.slice(0, 70)}" (${clip !== null ? "clipboard" : box ? "box" : "made"}: ${link.length} characters)`, "Copy link to this view (clipboard)");
      const t2 = await openPage(browser, eng, link, "view link opened", "view link restores the view");
      try {
        const done = await until(t2.page, (ds) => { const sv = window.__spatialViewer, S = sv && sv.state; return S && S.dsId === ds && S.m && window.SodaViewLink && !window.SodaViewLink.pending() && S.genes.length > 0; }, want.ds, 60000);
        await settled(t2.page); await t2.page.waitForTimeout(800);
        const got = await t2.page.evaluate(() => { const S = window.__spatialViewer.state; return { ds: S.dsId, layer: S.exprKey, genes: S.genes.map((g) => g.name), mode: S.umap ? S.umap.view : null,
          x: S.viewState.target[0], y: S.viewState.target[1], zoom: S.viewState.zoom, notes: window.SodaViewLink.notes() }; });
        const camOk = Math.abs(got.zoom - want.zoom) < 0.05 && Math.hypot(got.x - want.x, got.y - want.y) < 4 * Math.pow(2, -want.zoom) + 2;
        t2.ok(done && got.ds === want.ds && got.layer === want.layer && JSON.stringify(got.genes) === JSON.stringify(want.genes) && got.mode === want.mode && camOk,
          `restored: ${JSON.stringify({ ...got, x: Math.round(got.x), y: Math.round(got.y), zoom: +got.zoom.toFixed(3) })} (wanted ${JSON.stringify({ ...want, x: Math.round(want.x), y: Math.round(want.y), zoom: +want.zoom.toFixed(3) })})`);
        const r2 = await lookDrawn(t2.page, 0.01);
        t2.ok(r2.lit > 0.01 && !r2.err, `map drawn (lit ${r2.lit.toFixed(3)}${r2.status ? `; status "${r2.status}"` : ""})`);
        await t2.shot("viewlink_opened");
        checkClean(t2, "view link restores the view");
      } finally { await t2.ctx.close(); }
    });
    const r = await look(P);
    t.ok(!r.err, `no error line at the end (${r.status || "none"})`);
    checkClean(t, "walk: no console / page error, no failed request");
  } catch (e) { t.ok(false, `walk crashed: ${String(e.message).split("\n")[0]}`); }
  finally { await t.ctx.close(); }
}

// ---- cell outlines (outlines.js): a cells layer with segmentation boundaries, zoomed in, Outlines on ------
const hasShapes = (id) => D.exists(`${id}/expr/cells/shapes/cell/index.json`);
let OUTLINE_DS = null;                                                        // xenium_pub first (found before the run)
async function cellOutlines(browser, eng, urls) {
  const id = OUTLINE_DS, feature = "cell outlines (Outlines on, zoomed in)";
  const t = await openPage(browser, eng, `${urls.page}/index.html?static=${urls.data}/plain&ds=${id}`, `${id} outlines`, feature);
  const P = t.page, ev = (f, a) => P.evaluate(f, a);
  const idle = () => P.waitForFunction(() => { const s = window.__spatialViewer.outlines.stats(); if (s.busy) { window.__idle = 0; return false; } window.__idle = (window.__idle || 0) + 1; return window.__idle >= 4; }, null, { timeout: 60000, polling: 100 }).then(() => true, () => false);
  try {
    const ready = await until(P, (id) => { const sv = window.__spatialViewer, S = sv && sv.state; return S && S.dsId === id && S.m && sv.outlines && S.exprKey && S.geneIndex && S.geneIndex.length; }, id, 60000);
    if (!t.ok(ready, `opened ${id}`)) return;
    await settled(P); await tilesDone(P);
    // the cells layer, Outlines picked as a reader would
    const lay = await ev(() => {
      const S = window.__spatialViewer.state, sel = document.getElementById("exprSelect");
      if (S.exprKey !== "cells" && [...sel.options].some((o) => o.value === "cells")) { sel.value = "cells"; sel.dispatchEvent(new Event("change")); }
      return S.exprKey;
    });
    await until(P, () => { const S = window.__spatialViewer.state; return S.exprKey === "cells" && window.__spatialViewer.outlines.available(); }, null, 30000);
    const shown = await ev(() => { const s = document.getElementById("cellShapes"); s.value = "outlines"; s.dispatchEvent(new Event("change")); return !document.getElementById("rowCellShapes").hidden; });
    t.ok(shown, `the Cell shapes control shows on the cells layer (was on ${lay})`);
    // zoomed in so the median cell is 40 px across, at the image's centre
    await ev(() => { const sv = window.__spatialViewer, S = sv.state, e = S.m.expr.find((x) => x.key === S.exprKey), i = S.m.image;
      window.__idle = 0; sv.deck.props.onViewStateChange({ viewState: { ...S.viewState, target: [i.width / 2, i.height / 2, 0], zoom: Math.log2(40 / e.spotDiameterPx) }, viewId: "ortho" }); });
    const quiet = await idle();
    await tilesDone(P);
    const r = await lookDrawn(P, 0.01);
    const s = await ev(() => ({ ...window.__spatialViewer.outlines.stats(),
      layers: window.__spatialViewer.deck.layerManager.getLayers().filter((l) => /^ol-line-/.test(l.id) && l.props.visible !== false).length }));
    t.ok(quiet && s.mode === "outlines" && s.alpha === 1 && s.dot === 0 && s.layers > 0 && s.shown > 50 && s.sets[0].tiles > 0 && !s.sets.some((x) => x.bad),
      `outlines drawn: ${s.shown} shown, ${s.sets.map((x) => `${x.key}: ${x.tiles} tiles, ${x.polygons} polygons, ${x.bad} bad`).join("; ")} (alpha ${s.alpha}, dots ${s.dot}, ${s.layers} outline layers drawn, mode ${s.mode}, ${quiet ? "settled" : "STILL LOADING after 60 s"}${s.note ? `; "${s.note}"` : ""})`);
    t.ok(r.lit > 0.01 && !r.err, `map drawn (lit ${r.lit.toFixed(3)}${r.status ? `; status "${r.status}"` : ""})`);
    const k = `${eng.group}|${id} outlines`; (thumbs[k] = thumbs[k] || {})[eng.name] = r.thumb;
    await t.shot(`${id}_outlines`);
    checkClean(t, feature);
  } catch (e) { t.ok(false, `crashed: ${String(e.message).split("\n")[0]}`); }
  finally { await t.ctx.close(); }
}

// ---- when things fail: the page says so in plain words and stays usable ------------------------------------
// (1) the layer's gene list (genes.json) will not load: the layer still opens, with a message; (2) Draw a region's
// worker fails as it starts: the run stops with a message, never a "Computing…" that stays; (3) a Save JSON file and
// a link made by an older page, with parts since cut (Measure's meas, the flight speed, cell sets): they open with
// those parts ignored, and a cell set's lowest member gets its own name back.
async function hardening(browser, eng, urls) {
  const id = "visiumhd_pub", url = `${urls.page}/index.html?static=${urls.data}/plain&ds=${id}`;
  const errs = (t) => t.log.pageErrors.filter((e) => !/the worker fails/.test(e));
  // (1)
  let feature = "a layer whose gene list will not load: message, page usable";
  let t = await openPage(browser, eng, url, "genes.json refused", feature,
    (page) => page.route(/\/expr\/[^/]+\/genes\.json/, (rt) => rt.fulfill({ status: 500, body: "refused by the test" })));
  try {
    const P = t.page;
    const ok = await until(P, () => { const S = window.__spatialViewer && window.__spatialViewer.state, st = document.getElementById("status");
      return S && S.m && S.exprKey && st && !st.hidden && /could not be loaded/.test(st.textContent); }, null, 90000);
    const r = await P.evaluate(() => { const S = window.__spatialViewer.state;
      return { status: document.getElementById("status").textContent.slice(0, 200), genes: S.geneIndex.length };
    });
    t.ok(ok && r.genes === 0 && /list of this layer could not be loaded/.test(r.status), `the page says the gene list could not be loaded ("${r.status}")`, feature);
    const lasso = await until(P, () => { const S = window.__spatialViewer.state; return !document.getElementById("panelLasso").hidden && !!(S.lasso && S.lasso.pos); }, null, 60000);
    t.ok(lasso, "the rest of the layer still opened (Draw a region has its bins)", feature);
    const e = errs(t);
    t.ok(!e.length, `no page error${e.length ? ": " + e.slice(0, 3).join(" / ") : ""}`, feature);
    await t.shot("hardening_genes_refused");
  } catch (e) { t.ok(false, `crashed: ${String(e.message).split("\n")[0]}`, feature); }
  finally { await t.ctx.close(); }
  // (2) and (3), on one page
  feature = "Draw a region worker that fails: message, no stuck Computing";
  t = await openPage(browser, eng, url, "worker fails", feature,
    (page) => page.route("**/lasso_worker.js*", (rt) => rt.fulfill({ status: 200, contentType: "text/javascript", body: "throw new Error('the worker fails (test)');" })));
  try {
    const P = t.page, ev = (f, a) => P.evaluate(f, a);
    const ready = await until(P, () => { const S = window.__spatialViewer && window.__spatialViewer.state; return S && S.m && S.exprKey && S.lasso && S.lasso.pos && S.lasso.hvg; }, null, 90000);
    if (!t.ok(ready, "opened, with Draw a region ready", feature)) return;
    await settled(P);
    // two squares round bins a third and two thirds down the list, then Run (lassoLoad runs it, as Load JSON does)
    await ev(() => { const sv = window.__spatialViewer, pos = sv.state.lasso.pos, n = pos.length / 2;
      const sq = (f) => { const r = Math.floor(n * f), x = pos[2 * r], y = pos[2 * r + 1], d = 300; return [[x - d, y - d], [x + d, y - d], [x + d, y + d], [x - d, y + d]]; };
      sv.lassoLoad({ shapes: [{ name: "A", kind: "free", points: sq(0.37) }, { name: "B", kind: "free", points: sq(0.63) }], seed: 42, resolution: "auto", mode: "bins" }); });
    const done = await until(P, () => { const L = window.__spatialViewer.state.lasso, m = document.getElementById("lassoMsg").textContent;
      return L && !L.running && /UMAP stopped/.test(m); }, null, 90000);
    const r = await ev(() => ({ msg: document.getElementById("lassoMsg").textContent.slice(0, 200), running: window.__spatialViewer.state.lasso.running, cancel: !document.getElementById("lassoCancel").hidden }));
    t.ok(done && !r.running && !r.cancel && !/Computing/.test(r.msg), `the run stops with a message ("${r.msg}")`, feature);
    await t.shot("hardening_worker_fails");
    // (3) an older link, then an older Save JSON
    const feat3 = "an older link and Save JSON (speed, Measure, cell sets) open";
    // the section's UMAP and clusters loaded (as a reader does with the UMAP button), so a clustering is on hand
    await ev(() => window.__spatialViewer.umapSetView("umap"));
    await until(P, () => { const U = window.__spatialViewer.state.umap; return U && U.view === "umap" && !U.anim; }, null, 60000);
    const res = await ev(() => { const U = window.__spatialViewer.state.umap; return U && U.res != null ? U.res : null; });
    const st = await ev((res) => { const S = window.__spatialViewer.state;
      return { dataset: S.dsId, layer: S.exprKey,
        view: { v: window.SodaViewLink.VIEW_V, mode: "tissue", tab: "analyze", meas: [{ k: "ruler", pts: [[10, 10], [200, 200]] }], fly: 8000 },
        ...(res != null ? { sectionClusters: { layer: S.exprKey, resolution: res, names: { [res]: { 0: "Old set", 2: "Two" } }, hidden: [],
                                                sets: { [res]: [{ ids: [0, 1], prev: "Zero" }] } } } : {}) }; }, res);
    await ev(async (st) => { location.hash = "v=" + (await window.SodaViewLink.encode(st)); }, st);
    const applied = await until(P, () => !window.SodaViewLink.pending() && window.SodaViewLink.notes() !== null, null, 60000);
    const a = await ev((res) => { const U = window.__spatialViewer.state.umap, n = U && res != null ? (U.names[res] || {}) : null;
      return { notes: window.SodaViewLink.notes(), n0: n && n[0], n2: n && n[2] }; }, res);
    t.ok(applied && (res == null || (a.n0 === "Zero" && a.n2 === "Two")),
      `an older link opens; its cell set is ignored and cluster 0 has its own name back (${JSON.stringify(a)}${res == null ? "; no clustering here" : ""})`, feat3);
    await ev((res) => { const sv = window.__spatialViewer, S = sv.state;
      sv.lassoLoad({ dataset: S.dsId, layer: S.exprKey, seed: 42, shapes: [], view: { meas: [{ k: "area", pts: [[0, 0], [9, 0], [9, 9]] }] },
        ...(res != null ? { sectionClusters: { layer: S.exprKey, resolution: res, names: { [res]: { 0: "File set" } }, hidden: [],
                                                sets: { [res]: [{ ids: [0, 3], prev: null }] } } } : {}) }); }, res);
    await P.waitForTimeout(1500);
    const b = await ev((res) => { const U = window.__spatialViewer.state.umap; return { n0: U && res != null ? (U.names[res] || {})[0] : null }; }, res);
    t.ok(res == null || b.n0 === undefined, `an older Save JSON opens; its cell set is ignored (cluster 0 named: ${JSON.stringify(b.n0)})`, feat3);
    const e = errs(t);
    t.ok(!e.length, `no page error${e.length ? ": " + e.slice(0, 3).join(" / ") : ""}`, feat3);
  } catch (e) { t.ok(false, `crashed: ${String(e.message).split("\n")[0]}`, feature); }
  finally { await t.ctx.close(); }
}

// ---- ?data= links: one dataset without datasets.json; the list with tags --------------------------------
async function links(browser, eng, urls) {
  const id = CAT_IDS[1] || CAT_IDS[0];
  let t = await openPage(browser, eng, `${urls.page}/index.html?data=${urls.data}/plain/&ds=${id}`, "?data= link", "?data= link");
  try {
    const opened = await until(t.page, (id) => { const sv = window.__spatialViewer; return sv && sv.state && sv.state.m && sv.state.m.id === id && sv.diag().cached > 0; }, id);
    await settled(t.page);
    const r = await lookDrawn(t.page, 0.005);
    t.ok(opened && r.lit > 0.005 && !r.err, `?data=<folder>/&ds=${id} opens with no datasets.json (lit ${r.lit.toFixed(3)})`);
    await t.shot("data_link");
    checkClean(t, "?data= link", [/\/datasets\.json$/]);
  } catch (e) { t.ok(false, `crashed: ${String(e.message).split("\n")[0]}`); }
  finally { await t.ctx.close(); }

  t = await openPage(browser, eng, `${urls.page}/index.html?data=${urls.data}/cat/`, "picker with tags", "dataset picker + tag filter");
  try {
    const P = t.page;
    const listed = await until(P, (n) => { const sv = window.__spatialViewer; return sv && sv.dsPicker && sv.dsPicker.count() === n && sv.state.m; }, CAT_IDS.length);
    t.ok(listed, `datasets.json lists ${CAT_IDS.length}: ${(await ev2(P, () => window.__spatialViewer.state.datasets.map((d) => d.id))).join(", ")}`);
    // open the list, pick a value in the first tag filter, count the rows, click one
    const btn = await P.$(".dsp-btn");
    if (btn) await btn.click(); else await ev2(P, () => window.__spatialViewer.dsPicker.open && window.__spatialViewer.dsPicker.open());
    await P.waitForTimeout(300);
    const f = await ev2(P, () => {
      const sels = [...document.querySelectorAll(".dsp-fil select.dsp-fs")].filter((s) => s.offsetParent);
      const s = sels.find((x) => [...x.options].some((o) => /multiplex/i.test(o.textContent))) || sels[0];
      if (!s) return { filters: 0 };
      const o = [...s.options].find((o) => /multiplex/i.test(o.textContent)) || [...s.options].find((o) => o.value);
      s.value = o.value; s.dispatchEvent(new Event("change", { bubbles: true }));
      const rows = [...document.querySelectorAll(".dsp-o[role=option]")].filter((r) => r.offsetParent);
      return { filters: sels.length, value: o.textContent, rows: rows.length, first: rows[0] ? rows[0].textContent : "" };
    });
    const nMulti = CAT_IDS.filter((x) => tagOf(x).platform === "multiplex").length;
    t.ok(f.filters >= 1 && f.rows === nMulti, `tag filter "${f.value}" leaves ${f.rows} rows (expected ${nMulti}; ${f.filters} filters)`);
    await t.shot("picker_filtered");
    const row = await P.$(".dsp-o[role=option]:visible");
    if (row) await row.click();
    const want = CAT_IDS.find((x) => tagOf(x).platform === "multiplex");
    const opened = await until(P, (want) => window.__spatialViewer.state.m && window.__spatialViewer.state.m.id === want && window.__spatialViewer.diag().cached > 0, want);
    t.ok(opened, `clicking the row opens ${want} (now ${await ev2(P, () => window.__spatialViewer.state.m && window.__spatialViewer.state.m.id)})`);
    await settled(P);
    checkClean(t, "dataset picker + tag filter");
  } catch (e) { t.ok(false, `crashed: ${String(e.message).split("\n")[0]}`); }
  finally { await t.ctx.close(); }
}
const ev2 = (P, f, a) => P.evaluate(f, a);

// ---- Compare side by side (compare.js): a second dataset in a linked window ------------------------------
// where a window shows in the frame the group shares: µm from its own offset, and log2(screen px per µm)
const sharedOf = (P) => P.evaluate(() => {
  const S = window.__spatialViewer.state, C = window.SodaCompare.state(), um = S.m.pixelSizeUm > 0 ? S.m.pixelSizeUm : 1, vs = S.viewState;
  return { c: [vs.target[0] * um - C.off.o[0], vs.target[1] * um - C.off.o[1]], s: vs.zoom - Math.log2(um) };
});
const sameShared = (a, b) => Math.abs(a.s - b.s) < 1e-3 && Math.hypot(a.c[0] - b.c[0], a.c[1] - b.c[1]) < 0.5 + 2 * Math.pow(2, -a.s);
async function compare(browser, eng, urls) {
  const F = "Compare side by side (linked window)";
  const A_ID = "visiumhd_pub", B_ID = CAT_IDS.find((x) => x !== A_ID) || A_ID;
  if (!CAT_IDS.includes(A_ID)) { console.log("  (visiumhd_pub is not in the catalog: Compare not checked)"); return; }
  const t = await openPage(browser, eng, `${urls.page}/index.html?data=${urls.data}/cat/&ds=${A_ID}`, "compare", F);
  const P = t.page, ev = (f, a) => P.evaluate(f, a);
  let b = null;
  const bLog = { pageErrors: [], console: [], bad: [] };
  try {
    // the second window is a page of the same context: config.js answered there too
    await t.ctx.route("**/config.js*", (rt) => rt.fulfill({ status: 200, contentType: "text/javascript", body: "window.SpatialViewerConfig = {};" }));
    const ready = (pg, id) => until(pg, (id) => { const sv = window.__spatialViewer, S = sv && sv.state; return S && S.m && S.m.id === id && S.exprKey && S.geneIndex && S.geneIndex.length && sv.diag().cached > 0 && window.SodaCompare.state().ready; }, id);
    t.ok(await ready(P, A_ID), `${A_ID} opens from the catalog`);
    await settled(P);
    // the Compare button: in the top bar, just before Fit; it brings up the Share tab's Compare section
    const btn = await ev(() => { const b = document.getElementById("btnCompare"), s = b && getComputedStyle(b);
      return b ? { inBar: b.parentElement && b.parentElement.id === "topBar", next: b.nextElementSibling ? b.nextElementSibling.id : null, shown: s.display !== "none" && !!b.offsetParent, text: b.textContent } : null; });
    t.ok(btn && btn.inBar && btn.next === "btnFit" && btn.shown && btn.text === "Compare", `Compare button in the top bar before Fit: ${JSON.stringify(btn)}`);
    await P.click("#btnCompare");
    await P.waitForTimeout(300);
    const sec = await ev(() => { const L = window.SodaLayout.state(), q = (id) => document.getElementById(id), c = q("panelCompare");
      const before = (a, b) => !!(a && b && (a.compareDocumentPosition(b) & Node.DOCUMENT_POSITION_FOLLOWING));
      return { tab: L.tab, open: L.open.includes("panelCompare"), shown: !!(c && c.offsetParent), order: before(q("panelSave"), c) && before(c, q("panelDataset")),
               canLink: window.SodaCompare.state().canLink, peers: document.getElementById("cmpPeers").textContent }; });
    t.ok(sec.tab === "share" && sec.open && sec.shown && sec.order && sec.canLink, `the button opens Compare side by side in the Share tab, under Save and load: ${JSON.stringify(sec)}`);
    const mem0 = await ev(() => window.__spatialViewer.diag());
    t.ok(mem0.memWindows === 1 && mem0.cacheMaxMB === 1024, `alone: the whole tile cache (${mem0.cacheMaxMB} MB, ${mem0.memWindows} window)`);
    // open the other dataset beside it (as the dataset list in the section does)
    const popup = t.ctx.waitForEvent("page", { timeout: 30000 });
    await ev((id) => window.SodaCompare.open(id), B_ID);
    b = await popup;
    b.on("pageerror", (e) => bLog.pageErrors.push(String(e.message || e).split("\n")[0]));
    b.on("console", (m) => { if (m.type() === "error") bLog.console.push({ text: m.text(), url: (m.location() || {}).url || "" }); });
    b.on("response", (r) => { if (r.status() >= 400) bLog.bad.push({ url: r.url(), status: r.status() }); });
    const bReady = await ready(b, B_ID);
    const joined = await until(P, () => { const C = window.SodaCompare.state(); return C.peers.length === 1 && C.peers[0].ready && C.linked; }, null, 60000);
    const bJoined = await until(b, () => { const C = window.SodaCompare.state(); return C.peers.length === 1 && C.linked; }, null, 30000);
    t.ok(bReady && joined && bJoined, `${B_ID} opens in a second window, both join one group and are linked (${await ev(() => document.getElementById("cmpPeers").textContent)})`);
    await settled(P); await settled(b);
    const rA = await lookDrawn(P, 0.01), rB = await lookDrawn(b, 0.01);
    t.ok(rA.lit > 0.01 && !rA.err && rB.lit > 0.01 && !rB.err, `both windows draw (lit ${rA.lit.toFixed(3)} and ${rB.lit.toFixed(3)})`);
    await t.shot("compare_a");
    if (SHOTS) await b.screenshot({ path: path.join(SHOTS, eng.name, "compare_b.png") }).catch(() => {});
    // memory: the two windows split the ~1 GB tile cache (and Draw a region's kept chunks)
    const mA = await ev(() => window.__spatialViewer.diag()), mB = await b.evaluate(() => window.__spatialViewer.diag());
    t.ok(mA.memWindows === 2 && mB.memWindows === 2 && mA.cacheMaxMB === 512 && mB.cacheMaxMB === 512 && mA.cachedMB <= 512 && mB.cachedMB <= 512,
      `two windows split the tile cache: ${mA.cacheMaxMB} + ${mB.cacheMaxMB} MB (holding ${mA.cachedMB} + ${mB.cachedMB} MB)`);
    // linked: the centres lined up when the second window came; pan and zoom A, B goes to the same place in µm
    const s0A = await sharedOf(P), s0B = await sharedOf(b);
    t.ok(sameShared(s0A, s0B), `linked from the start: the same place and scale (${JSON.stringify(s0A)} / ${JSON.stringify(s0B)})`);
    await ev(() => { const sv = window.__spatialViewer, S = sv.state, i = S.m.image;
      sv.deck.props.onViewStateChange({ viewState: { ...S.viewState, target: [i.width * 0.4, i.height * 0.55, 0], zoom: S.viewState.zoom + 1 }, viewId: "ortho" }); });
    const wantA = await sharedOf(P);
    const followB = await until(b, (w) => { const S = window.__spatialViewer.state, C = window.SodaCompare.state(), um = S.m.pixelSizeUm > 0 ? S.m.pixelSizeUm : 1, vs = S.viewState;
      const c = [vs.target[0] * um - C.off.o[0], vs.target[1] * um - C.off.o[1]], s = vs.zoom - Math.log2(um);
      return Math.abs(s - w.s) < 1e-3 && Math.hypot(c[0] - w.c[0], c[1] - w.c[1]) < 0.5 + 2 * Math.pow(2, -w.s); }, wantA, 15000);
    t.ok(followB && sameShared(wantA, await sharedOf(b)), `pan + zoom in the first window: the second follows (${JSON.stringify(await sharedOf(b))})`);
    // and the other way round
    await b.evaluate(() => { const sv = window.__spatialViewer, S = sv.state;
      sv.deck.props.onViewStateChange({ viewState: { ...S.viewState, target: [S.viewState.target[0] * 1.1, S.viewState.target[1] * 0.9, 0], zoom: S.viewState.zoom - 0.5 }, viewId: "ortho" }); });
    const wantB = await sharedOf(b);
    const followA = await until(P, (w) => { const S = window.__spatialViewer.state, C = window.SodaCompare.state(), um = S.m.pixelSizeUm > 0 ? S.m.pixelSizeUm : 1, vs = S.viewState;
      const c = [vs.target[0] * um - C.off.o[0], vs.target[1] * um - C.off.o[1]], s = vs.zoom - Math.log2(um);
      return Math.abs(s - w.s) < 1e-3 && Math.hypot(c[0] - w.c[0], c[1] - w.c[1]) < 0.5 + 2 * Math.pow(2, -w.s); }, wantB, 15000);
    t.ok(followA && sameShared(wantB, await sharedOf(P)), `pan + zoom in the second window: the first follows (${JSON.stringify(await sharedOf(P))})`);
    // a gene picked in the first window: it follows when the other dataset has it, else that window names it as missing
    const namesB = new Set((await b.evaluate(() => window.__spatialViewer.state.geneIndex.map((g) => g.n.toLowerCase()))));
    const keyA = await ev(() => window.__spatialViewer.state.exprKey);
    const gA = await geneOf(A_ID, keyA);
    const common = namesB.has(gA.name.toLowerCase()) ? gA.name
      : await ev((nb) => { const s = new Set(nb); const g = window.__spatialViewer.state.geneIndex.find((x) => s.has(x.n.toLowerCase()) && !(x.b <= 0)); return g ? g.n : null; }, [...namesB]);
    const gene = common || gA.name;
    await pickFromSearch(P, "geneSearch", "geneResults", gene);
    await until(P, () => window.__spatialViewer.state.genes.length > 0, null, 60000);
    if (common) {
      const got = await until(b, (n) => window.__spatialViewer.state.genes.some((g) => g.name.toLowerCase() === n.toLowerCase()), gene, 60000);
      t.ok(got, `gene ${gene} picked in the first window is drawn in the second (${(await b.evaluate(() => window.__spatialViewer.state.genes.map((g) => g.name))).join(", ") || "none"})`);
    } else {
      const said = await until(b, (n) => document.getElementById("cmpMsg").textContent.includes(n), gene, 30000);
      t.ok(said, `gene ${gene} is not in ${B_ID}: the second window says so ("${await b.evaluate(() => document.getElementById("cmpMsg").textContent)}")`);
    }
    await t.shot("compare_gene");
    // close the second window: the first is alone again, with the whole cache. Its "bye" (pagehide) says so at once
    // where the browser sends it; otherwise its heartbeat stops and the first window drops it after 7 s
    await b.close();
    const alone = await until(P, () => window.SodaCompare.state().peers.length === 0 && window.__spatialViewer.diag().memWindows === 1, null, 30000);
    const mE = await ev(() => ({ ...window.__spatialViewer.diag(), msg: document.getElementById("cmpMsg").textContent, btn: document.getElementById("btnCompare").textContent }));
    t.ok(alone && /was closed|stopped answering/.test(mE.msg) && mE.cacheMaxMB === 1024 && mE.btn === "Compare", `closing the second window: the first is alone (${mE.msg || "no message"}), the whole cache back (${mE.cacheMaxMB} MB)`);
    // the second window's own errors (the first's are checked below)
    const allow = (u) => OPTIONAL.some((re) => re.test(String(u).split("?")[0]));
    const bCons = bLog.console.filter((c) => !(/Failed to load resource|status of 404/.test(c.text) && (!c.url || allow(c.url)))
      && !(eng.type === pw.firefox && FIREFOX_QUIET.some((re) => re.test(c.text))));
    const bBad = bLog.bad.filter((x) => !allow(x.url));
    t.ok(!bLog.pageErrors.length && !bCons.length && !bBad.length, `second window: no page error, console error or failed request${bLog.pageErrors.length + bCons.length + bBad.length
      ? ": " + [...bLog.pageErrors, ...bCons.map((c) => `${c.text} ${c.url}`), ...bBad.map((x) => `${x.status} ${x.url}`)].slice(0, 3).join(" / ") : ""}`);
    checkClean(t, F);
  } catch (e) { t.ok(false, `crashed: ${String(e.message).split("\n")[0]}`); }
  finally { if (b && !b.isClosed()) await b.close().catch(() => {}); await t.ctx.close(); }
}

// ---- phone layout ----------------------------------------------------------------------------------------
async function phone(browser, eng, urls) {
  const id = "visiumhd_pub";
  const t = await openPage(browser, eng, `${urls.page}/index.html?static=${urls.data}/plain&ds=${id}`, "phone 375 x 812",
    eng.noTouch ? "phone-sized window 375 x 812, no touch" : "phone layout 375 x 812, touch");
  try {
    const P = t.page;
    await until(P, (id) => { const sv = window.__spatialViewer; return sv && sv.state && sv.state.m && sv.state.m.id === id; }, id);
    await until(P, () => window.__spatialViewer.diag().cached > 0);
    await settled(P); await tilesDone(P);
    const r = await lookDrawn(P, 0.05);
    const L = await ev2(P, () => ({ sw: document.documentElement.scrollWidth, iw: innerWidth, ih: innerHeight, touch: window.__spatialViewer.diag().touch,
      rail: (() => { const e = document.querySelector("[role=tablist]"); if (!e) return null; const b = e.getBoundingClientRect(); return { w: b.width, h: b.height, vis: b.width > 0 && b.height > 0 && b.bottom <= innerHeight + 1 }; })(),
      vp: (() => { const b = document.getElementById("viewport").getBoundingClientRect(); return { w: b.width, h: b.height }; })() }));
    t.ok(r.lit > 0.05, `slide drawn (lit ${r.lit.toFixed(3)}, decoder ${r.decoder})`);
    if (eng.noTouch) t.ok(L.touch === false, "no touch screen here (Playwright's Firefox cannot emulate one): the page does not claim one");
    else t.ok(L.touch === true, "the page sees a touch screen");
    t.ok(L.sw <= L.iw + 1, `no sideways scrolling (page ${L.sw} px wide, window ${L.iw})`);
    t.ok(L.rail && L.rail.vis, `tab rail on screen ${JSON.stringify(L.rail)}`);
    t.ok(L.vp.w >= 300 && L.vp.h >= 250, `map ${Math.round(L.vp.w)} x ${Math.round(L.vp.h)} px`);
    // Compare: no button in a phone's top bar; on a touch screen the section opens a plain new tab (no linked windows)
    const cmp = await ev2(P, () => { const b = document.getElementById("btnCompare"); return { display: b ? getComputedStyle(b).display : "missing", canLink: window.SodaCompare ? window.SodaCompare.state().canLink : null }; });
    t.ok(cmp.display === "none" && (eng.noTouch || cmp.canLink === false), `no Compare button on a phone (${JSON.stringify(cmp)})`);
    const k = `${eng.group}|phone`; (thumbs[k] = thumbs[k] || {})[eng.name] = r.thumb;
    await t.shot("phone_open");
    // a tab opens its sheet (tapped)
    const tab = await P.$("[role=tab]:not([aria-selected=true])");
    if (tab) { if (eng.noTouch) await tab.click(); else await tab.tap(); await P.waitForTimeout(500); }
    const sheet = await ev2(P, () => { const p = [...document.querySelectorAll("[role=tabpanel]")].find((x) => !x.hidden && x.getBoundingClientRect().height > 40); return p ? Math.round(p.getBoundingClientRect().height) : 0; });
    t.ok(sheet > 0, `tapping a tab shows its panel (${sheet} px high)`);
    // every section starts closed: open the first one on show that holds a drop-down
    await ev2(P, () => { const s = [...document.querySelectorAll("[role=tabpanel]:not([hidden]) section.panel.closed:not([hidden])")].find((x) => x.querySelector("select"));
      if (s) document.getElementById(`secBtn-${s.id}`).click(); });
    await P.waitForTimeout(300);
    const sid =await ev2(P, () => { const s = [...document.querySelectorAll("select")].find((x) => x.id && x.getClientRects().length && x.getBoundingClientRect().bottom < innerHeight - 80 && x.getBoundingClientRect().top > 0); return s ? s.id : null; });
    const cc = sid ? await controlContrast(P, `#${sid}`) : null;
    t.ok(cc && cc.bg < 0.2 && cc.ratio >= 4, `a drop-down (#${sid}) drawn dark and readable: ${JSON.stringify(cc)}`, "form controls readable (dark <select>)");
    await t.shot("phone_tab");
    let said = "";
    const png = await download(t, async () => { said = await ev2(P, () => window.SodaSnapshot._save({ scale: 1 })); }, 60000);
    const z = png && pngSize(png.buf);
    t.ok(!!z && z.w > 300, `Save image 1x on the phone: ${png ? png.name : "nothing"} ${z ? `${z.w} x ${z.h}` : ""} "${said.slice(0, 100)}"`, "Save image on the phone");
    const sel = await ev2(P, () => [...document.getElementById("snapScale").options].map((o) => `${o.textContent}${o.disabled ? " [off]" : ""}`));
    console.log(`        phone Save image sizes: ${sel.join(" | ")}`);
    // the phone's limits count the whole picture (map + name + legend): a size is offered only when all of it fits
    const Z = await ev2(P, () => window.SodaSnapshot._sizes ? window.SodaSnapshot._sizes() : null);
    // without touch (firefox-phone) the page uses the desktop limits, as on any narrow desktop window
    if (eng.noTouch) t.ok(Z && Z.touch === false, "no touch screen: Save image uses the desktop limits (the phone's limits are checked in WebKit and Chromium)");
    else if (t.ok(Z && Z.touch, "Save image knows it is on a touch screen", "Save image on the phone: whole picture within the limit")) {
      const fits = (z) => z.W * z.H <= Z.limits.areaTouch && Math.max(z.W, z.H) <= Z.limits.dimTouch;
      const offered = new Set(await ev2(P, () => [...document.getElementById("snapScale").options].filter((o) => !o.disabled).map((o) => Number(o.value))));
      t.ok(Z.sizes.every((z) => (z.W > z.w || z.H > z.h) && fits(z) === offered.has(z.s) && !z.why === fits(z)),
        `each size offered only when the whole picture fits (${Z.sizes.map((z) => `${z.s}x ${z.W} x ${z.H} = ${(z.W * z.H / 1e6).toFixed(1)} M px${offered.has(z.s) ? "" : " refused"}`).join("; ")}; limit ${Z.limits.areaTouch / 1e6} M px)`,
        "Save image on the phone: whole picture within the limit");
      // the largest size offered, saved: the file is the size counted (and within the limit). Not on a software
      // GPU (SwiftShader), where a large picture can exhaust its memory.
      const soft = await ev2(P, () => { const gl = window.__spatialViewer.deck.device.gl, d = gl.getExtension("WEBGL_debug_renderer_info"); return /SwiftShader|llvmpipe|softpipe|Software/i.test(d ? gl.getParameter(d.UNMASKED_RENDERER_WEBGL) : gl.getParameter(gl.RENDERER)); });
      const big = Z.sizes.filter((z) => offered.has(z.s)).pop();
      if (big && big.s > 1 && !soft) {
        let said2 = "";
        const f = await download(t, async () => { said2 = await ev2(P, (s) => window.SodaSnapshot._save({ scale: s }), big.s); }, 120000);
        const zz = f && pngSize(f.buf);
        t.ok(!!zz && Math.abs(zz.w - big.W) <= 2 && Math.abs(zz.h - big.H) <= 2 && zz.w * zz.h <= Z.limits.areaTouch,
          `Save image ${big.s}x on the phone: ${zz ? `${zz.w} x ${zz.h}` : "nothing"} (counted ${big.W} x ${big.H}) "${said2.slice(0, 80)}"`, "Save image on the phone: whole picture within the limit");
      }
      // a taller screen, where the map alone at some size fits the limit but the map with its name does not:
      // that size must be refused, with the phone's message
      let edge = null;
      for (let hgt = 560; hgt <= 1200 && !edge; hgt += 4) {
        await P.setViewportSize({ width: 375, height: hgt });
        await P.waitForTimeout(150);
        const Y = await ev2(P, () => window.SodaSnapshot._sizes());
        const z = Y.sizes.find((q) => q.w * q.h <= Y.limits.areaTouch && Math.max(q.w, q.h) <= Y.limits.dimTouch && q.W * q.H > Y.limits.areaTouch);
        if (z) edge = { hgt, z };
      }
      if (t.ok(!!edge, `a screen height where the map fits the phone's limit but the whole picture does not${edge ? ` (375 x ${edge.hgt}: ${edge.z.s}x map ${edge.z.w} x ${edge.z.h}, whole ${edge.z.W} x ${edge.z.H})` : ""}`, "Save image on the phone: whole picture within the limit")) {
        const o = await ev2(P, (s) => { window.dispatchEvent(new Event("resize")); const x = [...document.getElementById("snapScale").options].find((q) => Number(q.value) === s); return x ? { off: x.disabled, txt: x.textContent } : null; }, edge.z.s);
        t.ok(o && o.off && /too large for a phone or tablet's memory/.test(o.txt), `${edge.z.s}x refused before drawing: "${o ? o.txt : "no option"}"`, "Save image on the phone: whole picture within the limit");
      }
      await P.setViewportSize({ width: 375, height: 812 });
    }
    checkClean(t, eng.noTouch ? "phone-sized window 375 x 812, no touch" : "phone layout 375 x 812, touch");
  } catch (e) { t.ok(false, `crashed: ${String(e.message).split("\n")[0]}`); }
  finally { await t.ctx.close(); }
}

(async () => {
  await D.need(["visiumhd_pub"].concat(CAT_IDS));
  const can = [];
  for (const p of picked) if (await D.has(p[1])) can.push(p); else console.log(D.skipLine(p[1], `dataset ${p[0]}`));
  picked = can;
  for (const id of ["xenium_pub"].concat(allIds)) if ((await D.has(id)) && (await hasShapes(id))) { OUTLINE_DS = id; break; }
  console.log(`data: ${D.describe()}`);
  const data = await listen(handler(dataRoute), PORTS[0]), pageSrv = await listen(handler(pageRoute), PORTS[1]);
  const urls = { data: `http://127.0.0.1:${data.address().port}`, page: `http://127.0.0.1:${pageSrv.address().port}` };
  LOCAL = [urls.data + "/", urls.page + "/"];
  const phoneDev = { ...pw.devices["iPhone X"] };                    // 375 x 812, touch, 3x
  const ENGINES = [
    { name: "chromium", type: pw.chromium, group: "desktop", clipboard: true },
    { name: "webkit", type: pw.webkit, group: "desktop" },
    { name: "msedge", type: pw.chromium, launch: { channel: "msedge" }, group: "desktop", clipboard: true, optional: true },
    { name: "firefox", type: pw.firefox, group: "desktop", clipboard: true, optional: true },
    { name: "chromium-phone", type: pw.chromium, group: "phone", device: { ...phoneDev, userAgent: pw.devices["Pixel 7"].userAgent }, phoneOnly: true },
    { name: "webkit-phone", type: pw.webkit, group: "phone", device: phoneDev, phoneOnly: true },
    // Playwright's Firefox has no isMobile / hasTouch: a 375 x 812 window, nothing more (see phone())
    { name: "firefox-phone", type: pw.firefox, group: "phone", device: { viewport: { width: 375, height: 812 } }, phoneOnly: true, noTouch: true, optional: true },
  ].filter((e) => (ONLY ? ONLY.split(",").includes(e.name) : !e.onlyIfAsked));
  console.log(`datasets: ${picked.map(([k, id]) => `${k} = ${id}`).join(", ")}`);
  let deckOk = null, ran = [];
  const skipped = [];
  try {
    for (const eng of ENGINES) {
      let browser;
      try { browser = await eng.type.launch({ headless: true, ...(eng.launch || {}) }); }
      catch (e) { const why = String(e.message).split("\n")[0]; console.log(`SKIPPED ${eng.name}: cannot launch (${why})`); skipped.push(`${eng.name}: ${why.slice(0, 120)}`); continue; }
      console.log(`\n== ${eng.name} ${browser.version()}`);
      try {
        if (deckOk === null) {                         // deck.gl comes from the CDN: without it nothing can be checked
          const ctx = await browser.newContext({ serviceWorkers: "block" }), page = await ctx.newPage();
          await page.route("**/config.js*", (rt) => rt.fulfill({ status: 200, contentType: "text/javascript", body: "window.SpatialViewerConfig = {};" }));
          await page.goto(`${urls.page}/index.html`, { waitUntil: "load" }).catch(() => {});
          deckOk = await page.evaluate(() => !!window.deck).catch(() => false);
          await ctx.close();
          if (!deckOk) { console.log("SKIPPED: deck.gl did not load from the CDN (no internet?)"); break; }
        }
        ran.push(eng.name);
        if (eng.phoneOnly) { await phone(browser, eng, urls); continue; }
        if (!SKIP.includes("probe")) await probe(browser, eng, urls);
        if (!SKIP.includes("datasets")) for (const [kind, id] of picked) await checkDataset(browser, eng, urls, kind, id);
        if (!SKIP.includes("walk")) await walk(browser, eng, urls);
        if (!SKIP.includes("compare")) await compare(browser, eng, urls);
        if (!SKIP.includes("links")) await links(browser, eng, urls);
        if (!SKIP.includes("hardening")) await hardening(browser, eng, urls);
        if (!SKIP.includes("outlines")) { if (OUTLINE_DS) await cellOutlines(browser, eng, urls); else console.log("  (no export with cell outlines: outlines not checked)"); }
      } finally { await browser.close(); }
    }
  } finally {
    data.close(); pageSrv.close();
    fs.rmSync(catDir, { recursive: true, force: true });
  }
  if (deckOk === false) process.exit(0);
  // the same views, compared with Chromium (desktop) / chromium-phone (phone): mean difference per channel, 0-255
  const diffs = [];
  for (const [k, by] of Object.entries(thumbs)) {
    const ref = by[k.startsWith("phone|") ? "chromium-phone" : "chromium"];
    for (const [e, th] of Object.entries(by)) {
      if (e === "chromium" || e === "chromium-phone" || !ref) continue;
      const d = diffThumb(ref, th);
      const okd = d !== null && d <= 8;
      if (okd) pass++; else { fail++; failures.push(`${e} | ${k.split("|")[1]} | map differs from Chromium's by ${d === null ? "?" : d.toFixed(2)} (allowed 8)`); }
      mark("same picture as Chromium (pixel diff)", e, okd);
      diffs.push(`${e} ${k.split("|")[1]} ${d === null ? "?" : d.toFixed(2)}`);
    }
  }
  if (diffs.length) console.log(`\npixel diff against Chromium (mean |difference| on a 64 x 40 grid, 0-255; allowed 8):\n  ${diffs.join("\n  ")}`);
  // the feature x browser table
  if (ran.length) {
    const w = Math.max(...Object.keys(table).map((k) => k.length)) + 2;
    console.log(`\n${"feature".padEnd(w)}${ran.map((e) => e.padEnd(16)).join("")}`);
    for (const [f, by] of Object.entries(table)) console.log(`${f.padEnd(w)}${ran.map((e) => String(by[e] || "-").padEnd(16)).join("")}`);
    if (SHOTS) fs.writeFileSync(path.join(SHOTS, "table.json"), JSON.stringify({ browsers: ran, skipped, table, diffs }, null, 1));
  }
  if (failures.length) console.log("\nFailed:\n" + failures.map((f) => "  " + f).join("\n"));
  console.log(`\nbrowsers (${ran.join(", ") || "none"}${skipped.length ? `; skipped ${skipped.map((s) => s.split(":")[0]).join(", ")}` : ""}): ${pass} passed, ${fail} failed`);
  process.exit(fail ? 1 : 0);
})();
