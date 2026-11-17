#!/usr/bin/env node
// Phone browser test: opens the viewer on the demo export in Playwright's WebKit (iPhone
// profile) and Chromium (phone emulation) and checks the image tile path, including what
// happens when a step fails: createImageBitmap refusing or missing, tiles refused, dropped
// or not pictures, no WebGL2, a lost graphics context. Also checks the ?diag=1 panel and
// that it never prints a token or password.
//
//   node tests/phone_browser_test.js                      (needs demo/export: run tests/smoke_test.py first)
//   node tests/phone_browser_test.js --export <dir> --page <dir> --only webkit --shots <dir> --ports 8864,8865
//
// Playwright is not a dependency of this repo. Install it anywhere (npm i playwright; npx
// playwright install webkit chromium) and point NODE_PATH or --playwright at its node_modules.
// Without it the test says so and exits 0.
//
// Playwright's WebKit is the WebKit engine, not iOS Safari: on Windows and Linux it does not
// use Apple's image decoders, Metal, or an iPhone's memory limits. It proves the fallbacks
// and messages work in a WebKit; it cannot prove what a given iPhone does.
"use strict";
const http = require("http");
const fs = require("fs");
const path = require("path");

const arg = (k, d) => { const i = process.argv.indexOf(k); return i > 0 ? process.argv[i + 1] : d; };
const REPO = path.resolve(__dirname, "..");
const EXPORT = path.resolve(arg("--export", path.join(REPO, "demo", "export")));
const PAGE = path.resolve(arg("--page", path.join(REPO, "web", "SpatialViewer")));
const ONLY = arg("--only", "");
const SHOTS = arg("--shots", "");
const SCEN = arg("--scenario", "");
const PORTS = arg("--ports", "0,0").split(",").map(Number);      // data server, page server; 0 = any free port

let pw = null;
for (const base of [arg("--playwright", ""), ...(process.env.NODE_PATH || "").split(path.delimiter), path.join(REPO, "node_modules")].filter(Boolean)) {
  try { pw = require(path.join(base, "playwright")); break; } catch (e) { /* next */ }
}
if (!pw) { try { pw = require("playwright"); } catch (e) { /* none */ } }
if (!pw) { console.log("SKIPPED: playwright is not installed (see the top of this file)"); process.exit(0); }
if (!fs.existsSync(path.join(EXPORT, "demo_visiumhd", "manifest.json"))) {
  console.log(`SKIPPED: no demo export at ${EXPORT} (run tests/smoke_test.py first, or pass --export)`); process.exit(0);
}

const TYPES = { ".html": "text/html", ".js": "text/javascript", ".css": "text/css", ".json": "application/json", ".png": "image/png", ".jpg": "image/jpeg", ".bin": "application/octet-stream" };
const { sendFile, preflight } = require("./static_files");
function serve(root, port) {
  return new Promise((resolve) => {
    const srv = http.createServer((req, res) => {
      if (req.method === "OPTIONS") { preflight(req, res); return; }
      const rel = decodeURIComponent(req.url.split("?")[0]);
      const file = path.join(root, rel === "/" ? "index.html" : rel);
      if (!file.startsWith(root) || !fs.existsSync(file) || fs.statSync(file).isDirectory()) { res.writeHead(404, { "Access-Control-Allow-Origin": "*" }); res.end("not found"); return; }
      sendFile(req, res, file, TYPES[path.extname(file)]);            // with Range: gene packs are read by range
    });
    srv.listen(port || 0, "127.0.0.1", () => resolve(srv));
  });
}

const SECRET_SIG = "SIGNATURE-THAT-MUST-NOT-BE-SHOWN", TEST_PW = "local-test-" + Math.random().toString(36).slice(2, 10);
// what the page shows and what the map holds
const look = async (page) => page.evaluate(() => {
  const st = document.getElementById("status"), sv = window.__spatialViewer, out = {};
  out.status = st && !st.hidden ? st.textContent : "";
  out.err = !!(st && !st.hidden && st.classList.contains("err"));
  out.lit = -1;
  if (sv && sv.deck && sv.deck.device && !sv.deck.device.gl.isContextLost()) {
    const gl = sv.deck.device.gl; sv.deck.redraw("test");
    const w = gl.drawingBufferWidth, h = gl.drawingBufferHeight, px = new Uint8Array(w * h * 4);
    gl.bindFramebuffer(gl.FRAMEBUFFER, null); gl.readPixels(0, 0, w, h, gl.RGBA, gl.UNSIGNED_BYTE, px);
    let lit = 0; for (let i = 0; i < px.length; i += 4) if (px[i] > 8 || px[i + 1] > 8 || px[i + 2] > 8) lit++;
    out.lit = lit / (w * h);
  }
  const d = sv && sv.diag ? sv.diag() : null;
  out.decoder = d ? d.decoder : null; out.cacheMax = d ? d.cacheMax : null; out.touch = d ? d.touch : null;
  out.errors = window.__sodaDiag ? window.__sodaDiag.errors.map((e) => `${e.kind}: ${e.msg}`) : null;
  const dt = document.getElementById("diagText"); out.diag = dt ? dt.textContent : null;
  return out;
});

// each scenario: set up the fault, open the page, then say what must be true
const SCENARIOS = [
  { name: "normal", ds: "demo_visiumhd", check: (r, ok) => { ok(r.lit > 0.05, `image drawn (lit ${r.lit.toFixed(3)})`); ok(!r.err, `no error line (${r.status})`); ok(r.decoder === "createImageBitmap", `decoder ${r.decoder}`); } },
  { name: "normal-jpeg", ds: "demo_multiplex", check: (r, ok) => { ok(r.lit > 0.02, `image drawn (lit ${r.lit.toFixed(3)})`); ok(!r.err, `no error line (${r.status})`); } },
  { name: "createImageBitmap refuses the options", ds: "demo_visiumhd",
    init: () => { const real = window.createImageBitmap.bind(window); window.createImageBitmap = (src, o) => (o && o.premultiplyAlpha ? Promise.reject(new DOMException("options refused (test)", "InvalidStateError")) : real(src)); },
    check: (r, ok) => { ok(r.lit > 0.05, `image drawn through the fallback (lit ${r.lit.toFixed(3)})`); ok(r.decoder === "createImageBitmap without options", `decoder ${r.decoder}`); ok(!r.err, `no error line (${r.status})`); } },
  { name: "createImageBitmap always refuses", ds: "demo_multiplex",
    init: () => { window.createImageBitmap = () => Promise.reject(new DOMException("refused (test)", "InvalidStateError")); },
    check: (r, ok) => { ok(r.lit > 0.02, `image drawn from img elements (lit ${r.lit.toFixed(3)})`); ok(r.decoder === "an img element", `decoder ${r.decoder}`); } },
  { name: "createImageBitmap missing", ds: "demo_visiumhd",
    init: () => { delete window.createImageBitmap; window.createImageBitmap = undefined; },
    check: (r, ok) => { ok(r.lit > 0.05, `image drawn from img elements (lit ${r.lit.toFixed(3)})`); ok(r.decoder === "an img element", `decoder ${r.decoder}`); } },
  { name: "?decode=img", ds: "demo_visiumhd", extra: "&decode=img",
    check: (r, ok) => { ok(r.lit > 0.05, `image drawn (lit ${r.lit.toFixed(3)})`); ok(r.decoder === "an img element", `decoder ${r.decoder}`); } },
  { name: "tiles refused (403)", ds: "demo_visiumhd",
    route: async (page) => page.route("**/tiles/**", (rt) => rt.fulfill({ status: 403, headers: { "Access-Control-Allow-Origin": "*" }, body: "denied" })),
    check: (r, ok) => { ok(r.err && /refused the image tiles \(403/.test(r.status), `status line says so: "${r.status}"`); ok(r.lit < 0.01, "map is black (nothing to draw)"); } },
  { name: "tiles are not pictures", ds: "demo_visiumhd",
    route: async (page) => page.route("**/tiles/**", (rt) => rt.fulfill({ status: 200, headers: { "Access-Control-Allow-Origin": "*", "Content-Type": "text/html" }, body: "<html>sign in</html>" })),
    check: (r, ok) => { ok(r.err && /could not decode an image tile/.test(r.status), `status line says so: "${r.status}"`); } },
  { name: "connection drops, then returns", ds: "demo_visiumhd", wait: 14000,
    route: async (page) => { const t0 = Date.now(); await page.route("**/tiles/**", (rt) => (Date.now() - t0 < 5000 ? rt.abort("connectionreset") : rt.continue())); },
    mid: 3500, midCheck: (r, ok) => { ok(r.err && /could not be downloaded/.test(r.status), `while down, the status line says so: "${r.status}"`); },
    check: (r, ok) => { ok(r.lit > 0.05, `image drawn once the connection is back, without a reload (lit ${r.lit.toFixed(3)})`); ok(!/could not/.test(r.status), `error line cleared ("${r.status}")`); } },
  { name: "no WebGL2", ds: "demo_visiumhd",
    init: () => { const real = HTMLCanvasElement.prototype.getContext; HTMLCanvasElement.prototype.getContext = function (t, ...a) { return /^webgl/.test(t) ? null : real.call(this, t, ...a); }; },
    check: (r, ok) => { ok(r.err && /WebGL2/.test(r.status), `status line says so: "${r.status}"`); } },
  { name: "graphics context lost", ds: "demo_visiumhd",
    after: async (page) => page.evaluate(() => { const gl = window.__spatialViewer.deck.device.gl; gl.getExtension("WEBGL_lose_context").loseContext(); return new Promise((r) => setTimeout(r, 800)); }),
    check: (r, ok) => { ok(r.err && /graphics memory/.test(r.status), `status line says so: "${r.status}"`); } },
  { name: "diag panel", ds: "demo_multiplex", extra: "&diag=1", wait: 12000,
    check: (r, ok) => {
      ok(!!r.diag, "panel is on the page");
      const d = r.diag || "";
      for (const want of ["RESULT", "User agent", "WebGL2", "Largest texture", "createImageBitmap", "Manifest: HTTP 200", "Tile: HTTP 200", "JPEG 1024 x 1024, 8-bit grey",
        "createImageBitmap with options", "img element: ok", "texStorage2D + texSubImage2D", "Map pixels that are not black", "Decoded tiles kept", "Errors caught"]) ok(d.includes(want), `panel has "${want}"`);
      ok(/RESULT: these checks found no problem/.test(d), `result line: ${(d.match(/RESULT: [^\n]*/) || [""])[0]}`);
      ok(!/\?/.test(d.replace(/\? \(none yet\)/g, "")), "no address in the panel keeps its ?query");
    } },
  { name: "diag panel behind the password (token flow)", ds: "demo_visiumhd", extra: "&diag=1", noStatic: true, wait: 12000,
    route: async (page, urls) => {
      await page.route("**/config.js*", (rt) => rt.fulfill({ status: 200, contentType: "text/javascript",
        body: `window.SpatialViewerConfig = { tokenEndpoint: "/api/token", tileCacheSize: 400, datasets: { demo_visiumhd: { label: "demo" } } };` }));
      await page.route("**/api/token", (rt) => {
        const req = rt.request();
        const good = req.method() === "POST" && JSON.parse(req.postData() || "{}").password === TEST_PW;
        return good ? rt.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify({ base: urls.data + "/", sas: `sv=2024-11-04&sp=r&sig=${SECRET_SIG}`, expires: Date.now() + 3600e3 }) })
                    : rt.fulfill({ status: 401, contentType: "application/json", body: `{"error":"no session"}` });
      });
    },
    during: async (page) => { await page.waitForSelector("#login:not([hidden])", { timeout: 15000 }); await page.fill("#loginPw", TEST_PW); await page.press("#loginPw", "Enter"); },
    check: (r, ok) => {
      const d = r.diag || "";
      ok(r.lit > 0.05, `image drawn after signing in (lit ${r.lit.toFixed(3)})`);
      ok(/Manifest: HTTP 200/.test(d) && /Tile: HTTP 200/.test(d), "panel fetched the manifest and a tile with the pass");
      ok(/signed in through the token service; pass valid for \d+ more min/.test(d), "panel says how long the pass lasts");
      ok(!d.includes(SECRET_SIG) && !/sig=/.test(d), "panel does not show the signature");
      ok(!d.includes(TEST_PW), "panel does not show the password");
      ok(!(r.errors || []).join("\n").includes(SECRET_SIG), "error list does not hold the signature");
      ok(r.touch !== true || r.cacheMax <= 64, `on touch, config tileCacheSize 400 does not lift the phone's limit (${r.cacheMax})`);
    } },
].filter((s) => !SCEN || s.name.includes(SCEN));

(async () => {
  const data = await serve(EXPORT, PORTS[0]), pageSrv = await serve(PAGE, PORTS[1]);
  const urls = { data: `http://127.0.0.1:${data.address().port}`, page: `http://127.0.0.1:${pageSrv.address().port}` };
  const engines = [
    { name: "webkit-iphone", type: pw.webkit, device: pw.devices["iPhone 15"] },
    { name: "chromium-phone", type: pw.chromium, device: pw.devices["Pixel 7"] },
  ].filter((e) => !ONLY || e.name.startsWith(ONLY));
  let pass = 0, fail = 0;
  if (SHOTS) fs.mkdirSync(SHOTS, { recursive: true });
  for (const eng of engines) {
    let browser;
    try { browser = await eng.type.launch(); } catch (e) { console.log(`SKIPPED ${eng.name}: cannot launch (${String(e.message).split("\n")[0]})`); continue; }
    for (const sc of SCENARIOS) {
      const ctx = await browser.newContext({ ...eng.device });
      const page = await ctx.newPage();
      const ok = (cond, what) => { if (cond) pass++; else fail++; console.log(`${cond ? "  ok  " : "  FAIL"} ${eng.name} | ${sc.name} | ${what}`); };
      try {
        if (sc.init) await page.addInitScript(sc.init);
        if (sc.route) await sc.route(page, urls);
        await page.goto(`${urls.page}/index.html?${sc.noStatic ? "" : `static=${urls.data}&`}ds=${sc.ds}${sc.extra || ""}`, { waitUntil: "load" });
        if (sc.during) await sc.during(page);
        if (sc.mid) { await page.waitForTimeout(sc.mid); sc.midCheck(await look(page), ok); }
        await page.waitForTimeout((sc.wait || 6000) - (sc.mid || 0));
        if (sc.after) await sc.after(page);
        const r = await look(page);
        sc.check(r, ok);
        if (SHOTS) {
          await page.screenshot({ path: path.join(SHOTS, `${eng.name}_${sc.name.replace(/[^\w]+/g, "_")}.png`) });
          if (r.diag) fs.writeFileSync(path.join(SHOTS, `${eng.name}_${sc.name.replace(/[^\w]+/g, "_")}.txt`), r.diag);
        }
      } catch (e) { ok(false, `scenario crashed: ${String(e.message).split("\n")[0]}`); }
      await ctx.close();
    }
    await browser.close();
  }
  data.close(); pageSrv.close();
  console.log(`\n${pass} passed, ${fail} failed`);
  process.exit(fail ? 1 : 0);
})();
