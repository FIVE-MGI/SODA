#!/usr/bin/env node
// Offline test: the page keeps itself (sw.js, offline.js) and opens with no network, in Chromium and WebKit.
// The same steps in each browser, on a page server and a data server of its own (the data on another
// origin, as blob storage):
//   online    the demo dataset opens (?static=); the service worker installs and takes over; its caches
//             hold index.html, the ?v= files, the workers, the marker lists, deck.gl and umap-js; nothing
//             of the data is kept. The CDN bytes the browser received are saved (temp folder) and their
//             SHA-384 checked against libs.js.
//   keep      Save tab > Offline use > Keep this dataset: a size is shown first, then the files are kept;
//             a view link is made.
//   offline   context.setOffline(true) (WebKit: both servers refuse connections, see run()), reload: the page and the dataset open, tiles drawn; Draw a region
//             (two circles) and Run UMAP (umap-js in the worker); a marker list reads; the view link opens
//             in a second window; a dataset that was not kept says so.
//   update    back online, the page served with every ?v= string changed (a new publish): reload, the new
//             version's cache replaces the old one (only one soda-app cache of this folder), and the page
//             opens offline again from the new one.
//   remove    the kept copy removed: offline, the dataset no longer opens.
//   vendor    (a new context) the page server also serves vendor/ with the saved files: deck.gl and umap-js
//             come from there, not the CDN (no request to cdn.jsdelivr.net), and the service worker keeps
//             the vendor/ copies; offline, Run UMAP works from them. Then vendor.py --from <saved> checks and
//             copies them into a temp folder (exit 0), and --check on an empty one says missing (exit 1).
//
//   node tests/offline_test.js [--export demo/export] [--shots <dir>] [--only chromium,webkit] [--python <exe>]
//
// Playwright as for tests/browsers_test.js (NODE_PATH / --playwright, PLAYWRIGHT_BROWSERS_PATH). Without it,
// without the demo export, or without the CDN (no internet: the first visit needs it), SKIPPED, exit 0.
// The existing browser tests block service workers (serviceWorkers: "block") so that their page.route
// rules see every request; this one allows them and routes nothing.
"use strict";
const http = require("http");
const fs = require("fs");
const os = require("os");
const path = require("path");
const crypto = require("crypto");
const { spawnSync } = require("child_process");

const arg = (k, d) => { const i = process.argv.indexOf(k); return i > 0 ? process.argv[i + 1] : d; };
const REPO = path.resolve(__dirname, "..");
const EXPORT = path.resolve(arg("--export", path.join(REPO, "demo", "export")));
const PAGE = path.resolve(arg("--page", path.join(REPO, "web", "SpatialViewer")));
const SHOTS = arg("--shots", "");
const ONLY = arg("--only", "chromium,webkit").split(",").filter(Boolean);
const PYTHON = arg("--python", "");
const DS = "demo_visiumhd", OTHER = "demo_multiplex";

let pw = null;
for (const base of [arg("--playwright", ""), ...(process.env.NODE_PATH || "").split(path.delimiter), path.join(REPO, "node_modules")].filter(Boolean)) {
  try { pw = require(path.join(base, "playwright")); break; } catch (e) { /* next */ }
}
if (!pw) { try { pw = require("playwright"); } catch (e) { /* none */ } }
if (!pw) { console.log("SKIPPED: playwright is not installed (see tests/browsers_test.js)"); process.exit(0); }
if (!fs.existsSync(path.join(EXPORT, DS, "manifest.json"))) {
  console.log(`SKIPPED: no demo export at ${EXPORT} (run tests/smoke_test.py first, or pass --export)`); process.exit(0);
}

// the pinned libraries, from libs.js (as vendor.py reads them)
const LIBS = JSON.parse(/\/\* LIBS-JSON \*\/([\s\S]*?)\/\* END-LIBS-JSON \*\//.exec(fs.readFileSync(path.join(PAGE, "libs.js"), "utf-8"))[1]);
const CDN = (L) => `https://cdn.jsdelivr.net/npm/${L.path}`;
const sri = (buf) => "sha384-" + crypto.createHash("sha384").update(buf).digest("base64");
const TMP = fs.mkdtempSync(path.join(os.tmpdir(), "soda-offline-"));
const SAVED = path.join(TMP, "saved");               // the CDN bytes the test browser received, laid out as vendor/

// ---- servers: the data (CORS, as blob storage) and the page --------------------------------------------
const TYPES = { ".html": "text/html; charset=utf-8", ".js": "text/javascript", ".css": "text/css", ".json": "application/json", ".png": "image/png", ".jpg": "image/jpeg", ".bin": "application/octet-stream" };
const { sendFile, preflight } = require("./static_files");
const inside = (root, rel) => { const f = path.join(root, rel); return f.startsWith(root + path.sep) ? f : null; };
const catDir = path.join(TMP, "cat");
fs.mkdirSync(catDir);
fs.writeFileSync(path.join(catDir, "datasets.json"), JSON.stringify({ soda: 1, title: "Offline test", default: DS,
  datasets: [DS, OTHER].filter((id) => fs.existsSync(path.join(EXPORT, id, "manifest.json"))).map((id) => ({ id, label: `Demo ${id}` })) }));
function dataHandler(req, res) {
  if (NET.down) { req.socket.destroy(); return; }
  if (req.method === "OPTIONS") { preflight(req, res); return; }
  const u = decodeURIComponent(req.url.split("?")[0]);
  const m = /^\/plain(\/.*)$/.exec(u);
  const file = m ? (m[1] === "/datasets.json" ? path.join(catDir, "datasets.json") : inside(EXPORT, m[1])) : null;
  if (!file || !fs.existsSync(file) || fs.statSync(file).isDirectory()) { res.writeHead(404, { "Access-Control-Allow-Origin": "*" }); res.end("not found"); return; }
  sendFile(req, res, file, TYPES[path.extname(file)]);
}
const PG = { variant: "", vendor: false, hits: [] };
const NET = { down: false };              // WebKit's "offline": both servers refuse every connection
function pageHandler(req, res) {
  if (NET.down) { req.socket.destroy(); return; }
  const u = decodeURIComponent(req.url.split("?")[0]);
  PG.hits.push(u);
  const send = (code, type, body) => { res.writeHead(code, { "Content-Type": type, "Cache-Control": "no-store" }); res.end(req.method === "HEAD" ? "" : body); };
  if (u === "/config.js") { send(200, "text/javascript", "window.SpatialViewerConfig = {};"); return; }
  if (u === "/" || u === "/index.html") {
    let html = fs.readFileSync(path.join(PAGE, "index.html"), "utf-8");
    // a new publish: every ?v= string of the page changed (the files themselves are the same)
    if (PG.variant) html = html.replace(/(<(?:script|link)\b[^>]*?\s(?:src|href)="(?!https?:|vendor\/)[^"?]+\?v=[\w.-]+)"/g, `$1${PG.variant}"`);
    send(200, "text/html; charset=utf-8", html); return;
  }
  if (u.startsWith("/vendor/")) {
    const f = PG.vendor ? inside(SAVED, u.slice("/vendor/".length)) : null;
    if (!f || !fs.existsSync(f)) { send(404, "text/plain", "not found"); return; }
    sendFile(req, res, f, "text/javascript"); return;
  }
  const f = inside(PAGE, u.slice(1));
  if (!f || !fs.existsSync(f) || fs.statSync(f).isDirectory() || !/\.(js|css|json|html|png|svg|ico)$/.test(f)) { send(404, "text/plain", "not found"); return; }
  sendFile(req, res, f, TYPES[path.extname(f)]);
}
const listen = (h) => new Promise((resolve) => { const s = http.createServer(h); s.listen(0, "127.0.0.1", () => resolve(s)); });

// ---- results -------------------------------------------------------------------------------------------
let pass = 0, fail = 0;
const failures = [], notes = [], table = {};
const mark = (feature, eng, cond) => { const row = (table[feature] = table[feature] || {}); if (row[eng] !== "BROKEN") row[eng] = cond ? "works" : "BROKEN"; };

async function run(browserType, eng, urls) {
  let browser;
  try { browser = await browserType.launch(); } catch (e) { console.log(`  SKIPPED ${eng}: cannot launch (${String(e.message).split("\n")[0]})`); return "skipped"; }
  const ok = (cond, what, feat) => {
    if (cond) pass++; else { fail++; failures.push(`${eng} | ${what}`); }
    if (feat) mark(feat, eng, cond);
    console.log(`${cond ? "  ok  " : "  FAIL"} ${eng} | ${what}`);
    return cond;
  };
  const newCtx = async () => {
    const ctx = await browser.newContext({ viewport: { width: 1280, height: 800 }, serviceWorkers: "allow" });
    const log = { console: [], pageErrors: [], cdn: [], net: [] };
    ctx.on("request", (r) => { if (/cdn\.jsdelivr\.net/.test(r.url())) log.cdn.push(r.url()); log.net.push(r.url()); });
    return { ctx, log };
  };
  const watch = (page, log) => {
    page.on("console", (m) => { if (m.type() === "error") log.console.push({ text: m.text(), url: (m.location() || {}).url || "" }); });
    page.on("pageerror", (e) => log.pageErrors.push(String(e.message || e).split("\n")[0]));
    return page;
  };
  const shot = async (page, name) => {
    if (!SHOTS) return;
    const dir = path.join(SHOTS, eng); fs.mkdirSync(dir, { recursive: true });
    await page.screenshot({ path: path.join(dir, `${name}.png`) }).catch(() => {});
  };
  let untilErr = "";
  const until = (page, fn, a, ms = 45000) => page.waitForFunction(fn, a, { timeout: ms, polling: 200 }).then(() => true, (e) => { untilErr = String(e.message).split(String.fromCharCode(10))[0].slice(0, 200); return false; });
  const opened = (page, id = DS, ms = 45000) => until(page, (id) => { const sv = window.__spatialViewer; return sv && sv.state && sv.state.m && sv.state.m.id === id && sv.state.dsId === id; }, id, ms);
  const where = (page) => page.evaluate(() => { const sv = window.__spatialViewer, s = document.getElementById("status");
    return { m: sv && sv.state && sv.state.m ? sv.state.m.id : null, ds: sv && sv.state ? sv.state.dsId : null, status: s && !s.hidden ? s.textContent.slice(0, 160) : "" }; }).catch((e) => ({ err: String(e.message).slice(0, 160) }));
  const statusText = (page) => page.evaluate(() => { const s = document.getElementById("status"); return s && !s.hidden ? s.textContent : ""; });
  const lit = (page) => page.evaluate(() => {
    const sv = window.__spatialViewer;
    if (!sv || !sv.deck || !sv.deck.device) return -1;
    const gl = sv.deck.device.gl; sv.deck.redraw("test");
    const w = gl.drawingBufferWidth, h = gl.drawingBufferHeight, px = new Uint8Array(w * h * 4);
    gl.bindFramebuffer(gl.FRAMEBUFFER, null); gl.readPixels(0, 0, w, h, gl.RGBA, gl.UNSIGNED_BYTE, px);
    let n = 0; for (let i = 0; i < px.length; i += 4) if (px[i] > 8 || px[i + 1] > 8 || px[i + 2] > 8) n++;
    return n / (w * h);
  });
  const drawn = async (page, ms = 20000) => { const t0 = Date.now(); let v = await lit(page); while (!(v > 0.01) && Date.now() - t0 < ms) { await page.waitForTimeout(400); v = await lit(page); } return v; };
  const offState = (page) => page.evaluate(() => (window.SodaOffline ? window.SodaOffline.state() : null));
  const swReady = (page, ms = 90000) => page.evaluate((ms) => Promise.race([window.SodaOffline.ready, new Promise((r) => setTimeout(() => r("timeout"), ms))]), ms);
  const cacheInfo = (page) => page.evaluate(async () => {
    const scope = new URL("./", location.href).href, out = { app: [], lib: [], data: [], files: {} };
    for (const n of await caches.keys()) {
      const kind = n.startsWith(`soda-app:${scope}|`) ? "app" : n.startsWith(`soda-lib:${scope}|`) ? "lib" : n.startsWith(`soda-data:${scope}|`) ? "data" : null;
      if (!kind) continue;
      out[kind].push(n);
      out.files[n] = (await (await caches.open(n)).keys()).map((r) => r.url.replace(location.origin, ""));
    }
    return out;
  });
  // two circles on the map and Run UMAP (as tests/browsers_test.js does)
  async function lassoUmap(page) {
    const ev = (f, a) => page.evaluate(f, a);
    await ev(() => window.__spatialViewer.fitView && window.__spatialViewer.fitView());
    await page.waitForTimeout(500);
    const ready = await until(page, () => !document.getElementById("panelLasso").hidden && window.__spatialViewer.state.lasso && window.__spatialViewer.state.lasso.hvg, null, 30000);
    if (!ready) return { error: "Draw a region is not offered" };
    const box = await ev(() => { const r = document.getElementById("viewport").getBoundingClientRect(); return { x: r.left, y: r.top, w: r.width, h: r.height }; });
    for (const [fx, fy] of [[0.42, 0.45], [0.6, 0.55]]) {
      await ev(() => { if (window.__spatialViewer.state.lasso.tool !== "circle") document.getElementById("lassoCircle").click(); });
      const cx = box.x + box.w * fx, cy = box.y + box.h * fy, rad = Math.min(box.w, box.h) * 0.09;
      await page.mouse.move(cx, cy); await page.mouse.down();
      for (let i = 1; i <= 6; i++) await page.mouse.move(cx + (rad * i) / 6, cy);
      await page.mouse.up();
      await page.waitForTimeout(250);
    }
    await ev(() => { if (window.__spatialViewer.state.lasso.tool) document.dispatchEvent(new KeyboardEvent("keydown", { key: "Escape" })); });
    const shapes = await ev(() => window.__spatialViewer.state.lasso.shapes.map((s) => s.bins));
    await ev(() => document.getElementById("lassoRun").click());
    const done = await until(page, () => { const L = window.__spatialViewer.state.lasso; return L.result && !L.running; }, null, 240000);
    return ev((a) => {
      const L = window.__spatialViewer.state.lasso, R = L.result;
      return { shapes: a.shapes, done: a.done, emb: R && R.emb ? R.emb.length : 0, k: R && R.info ? R.info.nClusters : null, method: R && R.info ? R.info.method : null,
               msg: document.getElementById("lassoMsg").textContent.slice(0, 160) };
    }, { shapes, done });
  }
  const showOffline = (page) => page.evaluate(() => { if (window.SodaLayout) window.SodaLayout.show("panelOffline"); });
  const errorsOf = (log, allow) => [...log.pageErrors.map((e) => `page: ${e}`),
    ...log.console.filter((c) => !/\/vendor\//.test(c.url) && !/\/vendor\//.test(c.text) && !(allow && allow.test(`${c.text} ${c.url}`))).map((c) => `console: ${c.text.slice(0, 160)} ${c.url}`)];

  const pageUrl = `${urls.page}/?static=${encodeURIComponent(urls.data + "/plain")}&ds=${DS}`;
  const rootUrl = `${urls.data}/plain/${DS}/`;
  let swOk = false;
  // offline: Chromium's context.setOffline. WebKit's (Playwright) fails every navigation before the service
  // worker sees it ("WebKit encountered an internal error"), so there both servers refuse connections
  // instead: the network error a page meets when its site is out of reach (navigator.onLine stays true,
  // and the CDN stays reachable, which the kept libraries do not need)
  const EMULATE = eng !== "webkit";
  const offline = async (ctx, on) => { if (EMULATE) await ctx.setOffline(on); NET.down = !!on && !EMULATE; };
  if (!EMULATE) notes.push(`${eng}: "offline" = the page and data servers refuse connections (Playwright's WebKit setOffline fails navigations before the service worker)`);
  PG.variant = ""; PG.vendor = false;

  // ---- online: first visit, the service worker installs ----------------------------------------------
  const A = await newCtx();
  const page = watch(await A.ctx.newPage(), A.log);
  const bodies = {};
  page.on("response", async (r) => {
    const L = Object.values(LIBS).find((x) => CDN(x) === r.url());
    if (L && r.ok()) { try { bodies[L.path] = await r.body(); } catch (e) { /* a worker's: from the cache later */ } }
  });
  await page.goto(pageUrl, { waitUntil: "load" });
  ok(await opened(page), `online: ${DS} opens: ${JSON.stringify(await where(page))} ${untilErr}`, "first visit online");
  ok(await drawn(page) > 0.01, "online: tiles drawn", "first visit online");
  const st0 = await swReady(page);
  const s0 = await offState(page);
  swOk = ok(st0 === "ready", `the service worker is installed and in charge: ${JSON.stringify(s0)}`, "service worker installs");
  if (!swOk) {
    notes.push(`${eng}: the service worker did not install (${JSON.stringify(s0)}); the offline steps are not run in this browser`);
    await shot(page, "online_no_sw");
    await browser.close();
    return "no service worker";
  }
  const c0 = await cacheInfo(page);
  const appFiles = c0.app.length === 1 ? c0.files[c0.app[0]] : [];
  ok(c0.app.length === 1 && c0.app[0].endsWith(`|${s0.app}`), `one app cache, named after the page's ?v= hash ${s0.app}: ${c0.app.join(", ")}`, "service worker installs");
  const want = ["/index.html", "/viewer.js?v=", "/lasso_worker.js?v=", "/leiden.js?v=", "/libs.js?v=", "/offline.js?v=", "/de_worker.js?v=", "/plots_worker.js?v=", "/styles.css?v=", "/markers/rna_generic.json", "/markers/protein_generic.json", "/config.js?v="];
  const missing = want.filter((w) => !appFiles.some((f) => f.startsWith(w)));
  ok(!missing.length, `the app cache holds the page, its ?v= files, the workers and the marker lists (${appFiles.length} files); missing: ${missing.join(", ") || "none"}`, "service worker installs");
  const htmlHash = await page.evaluate(async (name) => { const r = await (await caches.open(name)).match(new URL("index.html", location.href).href); const t = r ? await r.text() : ""; return window.SodaOffline.hashRefs(window.SodaOffline.htmlRefs(t)); }, c0.app[0]);
  ok(htmlHash === s0.app, `the kept index.html is the version the cache is named after (${htmlHash})`, "service worker installs");
  const libFiles = c0.lib.length === 1 ? c0.files[c0.lib[0]] : [];
  ok(Object.values(LIBS).every((L) => libFiles.includes(CDN(L))), `deck.gl and umap-js kept from the CDN (no vendor/ here): ${libFiles.join(", ")}`, "libraries kept");
  ok(!c0.data.length, "no data kept without asking", "data kept only on request");
  // the CDN bytes this browser received: their SHA-384 is the pinned one; saved (temp folder) for the vendor step
  for (const L of Object.values(LIBS)) {
    let buf = bodies[L.path];
    if (!buf) buf = Buffer.from(await page.evaluate(async (u) => { const r = await (await caches.match(u)); return r ? Array.from(new Uint8Array(await r.arrayBuffer())) : []; }, CDN(L)));
    ok(buf.length === L.bytes && sri(buf) === L.sri, `${L.name} ${L.version} from the CDN: ${buf.length} bytes, ${sri(buf)} (libs.js: ${L.sri})`, "pinned SHA-384");
    if (buf.length) { const f = path.join(SAVED, L.path); fs.mkdirSync(path.dirname(f), { recursive: true }); fs.writeFileSync(f, buf); }
  }
  ok(!errorsOf(A.log).length, `online: no error: ${errorsOf(A.log).slice(0, 3).join(" | ") || "none"}`, "first visit online");
  await shot(page, "1_online");

  // ---- keep the dataset ---------------------------------------------------------------------------------
  await showOffline(page);
  ok(await page.isVisible("#offKeep"), "Save tab > Offline use offers Keep this dataset", "keep a dataset");
  await page.click("#offKeep");
  const asked = await until(page, () => !document.getElementById("offGo").hidden, null, 30000);
  const est = await page.textContent("#offMsg");
  ok(asked && /About [\d.]+ (kB|MB) .*in \d+ files/.test(est), `a size is shown before anything is kept: "${est}"`, "keep a dataset");
  await shot(page, "2_keep_estimate");
  await page.click("#offGo");
  const kept = await until(page, (root) => { const s = window.SodaOffline.state(); return !s.busy && s.kept.some((k) => k.root === root); }, rootUrl, 120000);
  const s1 = await offState(page);
  const k1 = s1.kept.find((k) => k.root === rootUrl) || {};
  ok(kept && k1.files > 10 && k1.bytes > 0, `kept: ${JSON.stringify(k1)}; "${await page.textContent("#offMsg")}"`, "keep a dataset");
  const c1 = await cacheInfo(page);
  const dataFiles = c1.data.length ? c1.files[c1.data[0]] : [];
  const hvgN = dataFiles.filter((f) => /\/hvg\/\d+_\d+\.bin$/.test(f)).length;
  ok(dataFiles.some((f) => f.endsWith(`${DS}/manifest.json`)) && hvgN > 0 && dataFiles.some((f) => /\/tiles\//.test(f)),
     `the kept copy holds the manifest, ${dataFiles.filter((f) => /\/tiles\//.test(f)).length} tiles looked at and ${hvgN} hvg chunks (Draw a region)`, "keep a dataset");
  ok(!dataFiles.some((f) => /\?(?!soda-offline-meta)/.test(f)), "no address with a query is kept", "keep a dataset");
  // zoom in a little (the mouse wheel), then a link to this view
  const vb = await page.evaluate(() => { const r = document.getElementById("viewport").getBoundingClientRect(); return { x: r.left + r.width / 2, y: r.top + r.height / 2 }; });
  await page.mouse.move(vb.x, vb.y); await page.mouse.wheel(0, -240); await page.waitForTimeout(800);
  const view = await page.evaluate(async () => { const m = await window.SodaViewLink.make(); return { url: m.url, zoom: window.__spatialViewer.state.viewState.zoom }; });
  await shot(page, "3_kept");

  // ---- offline -------------------------------------------------------------------------------------------
  await offline(A.ctx, true);
  A.log.console.length = 0; A.log.pageErrors.length = 0;
  await page.reload({ waitUntil: "load" });
  const onl = await page.evaluate(() => navigator.onLine);
  ok(await opened(page), `offline (navigator.onLine ${onl}): the page and ${DS} open`, "offline: page + dataset open");
  ok(await drawn(page) > 0.01, "offline: tiles drawn", "offline: page + dataset open");
  ok(await page.evaluate(() => !!window.deck && !!navigator.serviceWorker.controller), "offline: deck.gl loaded, the service worker answers", "offline: page + dataset open");
  const mk = await page.evaluate(async () => { const r = await fetch("markers/rna_generic.json", { cache: "no-cache" }); const j = r.ok ? await r.json() : null; return j && j.types ? j.types.length : 0; });
  ok(mk > 0, `offline: a marker list reads (${mk} types)`, "offline: marker lists");
  const u1 = await lassoUmap(page);
  ok(u1.done && u1.emb > 0, `offline: Draw a region + Run UMAP (umap-js in the worker): ${JSON.stringify(u1)}`, "offline: Draw a region + Run UMAP");
  await shot(page, "4_offline_umap");
  ok(!errorsOf(A.log, /net::ERR_INTERNET_DISCONNECTED|FetchEvent|Failed to load resource|network error|Load failed|offline/i).length,
     `offline: no page error: ${errorsOf(A.log, /net::ERR_INTERNET_DISCONNECTED|FetchEvent|Failed to load resource|network error|Load failed|offline/i).slice(0, 3).join(" | ") || "none"}`, "offline: page + dataset open");
  // the view link in a second window (as Compare and a pasted link open one)
  const p2 = watch(await A.ctx.newPage(), A.log);
  await p2.goto(view.url, { waitUntil: "load" }).catch(() => {});
  const v2 = await opened(p2) && await until(p2, (z) => { const v = window.__spatialViewer.state.viewState; return v && Math.abs(v.zoom - z) < 0.05; }, view.zoom, 20000);
  ok(v2, `offline: the view link opens in a second window, same dataset and zoom (${view.zoom.toFixed(2)})`, "offline: view link / second window");
  await shot(p2, "5_offline_view_link");
  await p2.close();
  // a dataset that was not kept: says so
  if (fs.existsSync(path.join(EXPORT, OTHER, "manifest.json"))) {
    const p3 = watch(await A.ctx.newPage(), { console: [], pageErrors: [] });
    await p3.goto(`${urls.page}/?static=${encodeURIComponent(urls.data + "/plain")}&ds=${OTHER}`, { waitUntil: "load" }).catch(() => {});
    const said = await until(p3, () => /not kept for offline use/.test((document.getElementById("status") || {}).textContent || ""), null, 30000);
    ok(said, `offline: ${OTHER} (not kept) says so: "${await statusText(p3)}"`, "offline: not-kept dataset explained");
    await p3.close();
  }

  // ---- a new publish: every ?v= changed ---------------------------------------------------------------
  await offline(A.ctx, false);
  PG.variant = "u2";
  await page.reload({ waitUntil: "load" });
  await opened(page);
  const st2 = await swReady(page);
  const s2 = await offState(page);
  const c2 = await cacheInfo(page);
  const files2 = c2.app.length === 1 ? c2.files[c2.app[0]] : [];
  ok(st2 === "ready" && s2.app !== s0.app && c2.app.length === 1 && c2.app[0].endsWith(`|${s2.app}`) && !c2.app.includes(c0.app[0]),
     `a newer ?v= replaces the cache: ${c0.app[0].split("|").pop()} -> ${c2.app.join(", ").split("|").pop()} (${c2.app.length} app cache: ${c2.app.map((n) => `${n.split("|").pop()} ${c2.files[n].length} files`).join(", ")})`, "new ?v= replaces the cache");
  ok(files2.some((f) => /\/viewer\.js\?v=[\w.-]+u2$/.test(f)) && !files2.some((f) => /\/viewer\.js\?v=[\w.-]+[^2]$/.test(f)),
     "the new cache holds the new version's files only", "new ?v= replaces the cache");
  ok(c2.lib.length === 1 && c2.lib[0] === c0.lib[0] && c2.data.length === 1, "the libraries and the kept dataset stay", "new ?v= replaces the cache");
  await offline(A.ctx, true);
  await page.reload({ waitUntil: "load" });
  const new2 = await page.evaluate(() => !!document.querySelector('script[src$="u2"]'));
  ok(await opened(page) && new2 && await drawn(page) > 0.01, "offline after the update: the new version opens, with the dataset", "new ?v= replaces the cache");
  await shot(page, "6_offline_after_update");

  // ---- remove ----------------------------------------------------------------------------------------------
  await offline(A.ctx, false);
  await showOffline(page);
  await until(page, () => !document.getElementById("offDrop").hidden, null, 10000);
  await page.click("#offDrop");
  const gone = await until(page, () => !window.SodaOffline.state().kept.length, null, 10000);
  ok(gone && !(await cacheInfo(page)).data.length, "Remove offline copy: the kept data is gone", "remove a kept dataset");
  await offline(A.ctx, true);
  await page.reload({ waitUntil: "load" });
  const refused = await until(page, () => /not kept for offline use/.test((document.getElementById("status") || {}).textContent || ""), null, 30000);
  ok(refused, `offline after Remove: the page opens, the dataset says it is not kept: "${await statusText(page)}"`, "remove a kept dataset");
  await offline(A.ctx, false);
  await A.ctx.close();

  // ---- vendor/: the libraries from the page's own site --------------------------------------------------
  PG.variant = ""; PG.vendor = true;
  const haveSaved = Object.values(LIBS).every((L) => fs.existsSync(path.join(SAVED, L.path)));
  if (!haveSaved) { ok(false, "vendor: the CDN bytes were not saved, so vendor/ cannot be tested"); await browser.close(); return "done"; }
  const B = await newCtx();
  const pv = watch(await B.ctx.newPage(), B.log);
  await pv.goto(pageUrl, { waitUntil: "load" });
  ok(await opened(pv) && await drawn(pv) > 0.01, "vendor: the page opens with deck.gl from vendor/", "vendor/ (no CDN)");
  ok(await pv.evaluate(() => window.SODA_LIBS_FROM && window.SODA_LIBS_FROM.deck === "local"), "vendor: libs.js says deck.gl came from vendor/", "vendor/ (no CDN)");
  const uv = await lassoUmap(pv);
  ok(uv.done && uv.emb > 0, `vendor: Run UMAP with umap-js from vendor/: ${JSON.stringify(uv)}`, "vendor/ (no CDN)");
  ok(!B.log.cdn.length, `vendor: nothing asked of the CDN: ${B.log.cdn.join(", ") || "none"}`, "vendor/ (no CDN)");
  const stv = await swReady(pv);
  const cv = await cacheInfo(pv);
  const libv = cv.lib.length === 1 ? cv.files[cv.lib[0]] : [];
  ok(stv === "ready" && Object.values(LIBS).every((L) => libv.includes(`/vendor/${L.path}`)), `vendor: the service worker kept the vendor/ copies: ${libv.join(", ")}`, "vendor/ (no CDN)");
  ok(!errorsOf(B.log).length, `vendor: no error: ${errorsOf(B.log).slice(0, 3).join(" | ") || "none"}`, "vendor/ (no CDN)");
  // keep the dataset here too (a new context keeps nothing of the first), for Run UMAP offline below
  await showOffline(pv);
  await pv.click("#offKeep");
  await until(pv, () => !document.getElementById("offGo").hidden, null, 30000);
  await pv.click("#offGo");
  await until(pv, (root) => { const s = window.SodaOffline.state(); return !s.busy && s.kept.some((k) => k.root === root); }, rootUrl, 120000);
  await offline(B.ctx, true);
  await pv.reload({ waitUntil: "load" });
  const ov = await until(pv, () => !!window.deck && window.SODA_LIBS_FROM && window.SODA_LIBS_FROM.deck === "local", null, 30000);
  ok(ov, "vendor, offline: the page loads deck.gl from the kept vendor/ copy", "vendor/ (no CDN)");
  const uvo = ov && await opened(pv) ? await lassoUmap(pv) : { error: "not open" };
  ok(uvo.done && uvo.emb > 0, `vendor, offline: Run UMAP with the kept vendor/ umap-js: ${JSON.stringify(uvo)}`, "vendor/ (no CDN)");
  await shot(pv, "7_vendor_offline");
  await offline(B.ctx, false);
  await B.ctx.close();
  // a vendor/ file that is not the pinned one is refused by the browser: the CDN copy is used
  const bad = path.join(SAVED, LIBS.deck.path), good = fs.readFileSync(bad);
  fs.writeFileSync(bad, Buffer.concat([good, Buffer.from("\n// changed\n")]));
  const C = await newCtx();
  const pc = watch(await C.ctx.newPage(), C.log);
  await pc.goto(pageUrl, { waitUntil: "load" });
  const fromC = await until(pc, () => !!window.deck, null, 30000) && await pc.evaluate(() => window.SODA_LIBS_FROM.deck);
  ok(fromC === "cdn" && C.log.cdn.some((u) => u === CDN(LIBS.deck)), `vendor: a changed vendor/deck.gl is refused (SHA-384) and the CDN copy runs (${fromC})`, "vendor/ (no CDN)");
  fs.writeFileSync(bad, good);
  await C.ctx.close();
  await browser.close();
  return "done";
}

(async () => {
  const data = await listen(dataHandler), pg = await listen(pageHandler);
  const urls = { data: `http://127.0.0.1:${data.address().port}`, page: `http://127.0.0.1:${pg.address().port}` };
  console.log(`offline test: page ${urls.page}, data ${urls.data}, export ${EXPORT}`);
  const done = {};
  // the first visit needs the CDN: without internet, skip
  try {
    const b = await pw.chromium.launch(), p = await b.newPage();
    const r = await p.request.get(CDN(LIBS.umap), { timeout: 20000 }).catch(() => null);
    await b.close();
    if (!r || !r.ok()) { console.log("SKIPPED: the CDN (cdn.jsdelivr.net) cannot be reached: no internet?"); process.exit(0); }
  } catch (e) { console.log(`SKIPPED: ${String(e.message).split("\n")[0]}`); process.exit(0); }
  for (const name of ONLY) {
    const bt = pw[name];
    if (!bt) { console.log(`  SKIPPED ${name}: not a Playwright browser`); continue; }
    console.log(`-- ${name}`);
    NET.down = false;
    try { done[name] = await run(bt, name, urls); }
    catch (e) { fail++; failures.push(`${name} | crashed: ${String(e.stack || e).split("\n").slice(0, 3).join(" ")}`); console.log(`  FAIL ${name} crashed: ${e.stack || e}`); done[name] = "crashed"; }
  }
  // vendor.py: copies the saved (checked) files, refuses nothing it should keep; --check on an empty folder fails
  if (PYTHON && Object.values(LIBS).every((L) => fs.existsSync(path.join(SAVED, L.path)))) {
    const dest = path.join(TMP, "site"), empty = path.join(TMP, "empty");
    fs.mkdirSync(dest); fs.mkdirSync(empty);
    const vp = path.join(PAGE, "vendor.py");
    const a = spawnSync(PYTHON, [vp, "--from", SAVED, "--dest", dest], { encoding: "utf-8" });
    const okFiles = Object.values(LIBS).every((L) => { const f = path.join(dest, "vendor", L.path); return fs.existsSync(f) && sri(fs.readFileSync(f)) === L.sri; });
    const r1 = a.status === 0 && okFiles; if (r1) pass++; else { fail++; failures.push(`vendor.py --from: exit ${a.status}: ${a.stdout} ${a.stderr}`); }
    console.log(`${r1 ? "  ok  " : "  FAIL"} vendor.py --from <saved> --dest <temp>: exit ${a.status}, files in place with the pinned SHA-384`);
    const c = spawnSync(PYTHON, [vp, "--check", "--dest", dest], { encoding: "utf-8" });
    const e = spawnSync(PYTHON, [vp, "--check", "--dest", empty], { encoding: "utf-8" });
    const r2 = c.status === 0 && e.status === 1 && /MISSING/.test(e.stdout); if (r2) pass++; else { fail++; failures.push(`vendor.py --check: ${c.status} / ${e.status}: ${e.stdout}`); }
    console.log(`${r2 ? "  ok  " : "  FAIL"} vendor.py --check: complete folder exit ${c.status}, empty folder exit ${e.status} (MISSING)`);
    fs.writeFileSync(path.join(TMP, "bad.js"), "x");
    const badDir = path.join(TMP, "badsrc"); for (const L of Object.values(LIBS)) { const f = path.join(badDir, L.path); fs.mkdirSync(path.dirname(f), { recursive: true }); fs.writeFileSync(f, "not the library"); }
    const d = spawnSync(PYTHON, [vp, "--from", badDir, "--dest", empty], { encoding: "utf-8" });
    const r3 = d.status === 1 && /REFUSED/.test(d.stdout) && !fs.existsSync(path.join(empty, "vendor", LIBS.deck.path)); if (r3) pass++; else { fail++; failures.push(`vendor.py wrong file: ${d.status}: ${d.stdout}`); }
    console.log(`${r3 ? "  ok  " : "  FAIL"} vendor.py --from <wrong files>: refused, nothing written (exit ${d.status})`);
  }
  NET.down = false;
  data.close(); pg.close();
  try { fs.rmSync(TMP, { recursive: true, force: true }); } catch (e) { /* temp */ }

  const engs = ONLY.filter((n) => done[n]);
  console.log("\nfeature table:");
  for (const [f, row] of Object.entries(table)) console.log(`  ${f.padEnd(40)} ${engs.map((e) => `${e}: ${row[e] || "-"}`).join("   ")}`);
  for (const n of notes) console.log(`  note: ${n}`);
  if (fail) { console.log(`\nOFFLINE TEST FAILED: ${fail} of ${pass + fail} checks`); failures.forEach((f) => console.log("  - " + f)); process.exit(1); }
  console.log(`\nOFFLINE TEST PASSED: ${pass} checks (${engs.map((e) => `${e} ${done[e]}`).join(", ")})`);
})();
