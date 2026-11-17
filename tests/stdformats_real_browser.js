#!/usr/bin/env node
// Real public standard-format datasets in the browser (tests/stdformats_real_test.py runs this; not part of the
// smoke test). Opens each export in Playwright's Chromium and checks that it is drawn: the tiles, the scale bar
// (µm when the manifest has a pixel size, image pixels otherwise, and its length matches the pixel size), the
// expression layer with a gene from the search on the slide, the cell outlines when the layer has them, the
// molecules layer when there is one. Screenshots: <id>.png (the whole section), <id>_gene.png (a gene),
// <id>_outlines.png (zoomed in until the median cell is 40 screen pixels, outlines drawn).
//
//   node tests/stdformats_real_browser.js --export <dir> --ids a,b,c --shots <dir> [--gpu]
//
// Prints one JSON line per dataset ("RESULT {...}") for the Python side, then PASSED / FAILED.
"use strict";
const http = require("http");
const fs = require("fs");
const path = require("path");
const { sendFile, preflight } = require("./static_files");

const arg = (k, d) => { const i = process.argv.indexOf(k); return i > 0 ? process.argv[i + 1] : d; };
const REPO = path.resolve(__dirname, "..");
const EXPORT = path.resolve(arg("--export", ""));
const PAGE = path.resolve(arg("--page", path.join(REPO, "web", "SpatialViewer")));
const SHOTS = arg("--shots", "");
const GPU = process.argv.includes("--gpu");
const IDS = arg("--ids", "").split(",").filter(Boolean);

let pw = null;
for (const base of [arg("--playwright", ""), ...(process.env.NODE_PATH || "").split(path.delimiter), path.join(REPO, "node_modules")].filter(Boolean)) {
  try { pw = require(path.join(base, "playwright")); break; } catch (e) { /* next */ }
}
if (!pw) { try { pw = require("playwright"); } catch (e) { /* none */ } }
if (!pw) { console.log("SKIPPED: playwright is not installed"); process.exit(0); }

const OPTIONAL = [/\/gates\.json$/, /\/markers\/[\w.-]+\.json$/, /\/vendor\//];
const TYPES = { ".html": "text/html", ".js": "text/javascript", ".css": "text/css", ".json": "application/json", ".png": "image/png", ".jpg": "image/jpeg", ".bin": "application/octet-stream" };
function serve(root) {
  return new Promise((resolve) => {
    const srv = http.createServer((req, res) => {
      const rel = decodeURIComponent(req.url.split("?")[0]);
      if (rel === "/config.js") { res.writeHead(200, { "Content-Type": "text/javascript" }); res.end("window.SpatialViewerConfig = { tokenEndpoint: \"\" };"); return; }
      const file = path.join(root, rel === "/" ? "index.html" : rel);
      if (req.method === "OPTIONS") { preflight(req, res); return; }
      if (!file.startsWith(root) || !fs.existsSync(file) || fs.statSync(file).isDirectory()) { res.writeHead(404, { "Access-Control-Allow-Origin": "*" }); res.end("not found"); return; }
      sendFile(req, res, file, TYPES[path.extname(file)]);
    });
    srv.listen(0, "127.0.0.1", () => resolve(srv));
  });
}

const results = [];
const check = (cond, what, got) => { results.push({ ok: !!cond, what }); console.log(`  ${cond ? "ok  " : "FAIL"} ${what}${cond || got === undefined ? "" : "  got " + JSON.stringify(got)}`); return !!cond; };

(async () => {
  const data = await serve(EXPORT), pageSrv = await serve(PAGE);
  const DATA = `http://127.0.0.1:${data.address().port}`, URL0 = `http://127.0.0.1:${pageSrv.address().port}/index.html`;
  const browser = await pw.chromium.launch({ args: GPU ? ["--use-angle=d3d11", "--ignore-gpu-blocklist", "--enable-gpu"] : ["--use-angle=swiftshader", "--enable-unsafe-swiftshader", "--ignore-gpu-blocklist"] });
  const ctx = await browser.newContext({ viewport: { width: 1366, height: 768 }, serviceWorkers: "block" });
  const shot = async (page, name) => { if (SHOTS) { fs.mkdirSync(SHOTS, { recursive: true }); await page.screenshot({ path: path.join(SHOTS, name) }); } };
  for (const id of IDS) {
    console.log(id);
    const errors = [];
    const page = await ctx.newPage();
    page.on("console", (m) => { if (m.type() === "error" && !/^Failed to load resource/.test(m.text())) errors.push(`console: ${m.text()}`); });
    page.on("pageerror", (e) => errors.push(`pageerror: ${e.message}`));
    page.on("response", (r) => { if (r.status() >= 400 && !/favicon/.test(r.url()) && !OPTIONAL.some((re) => re.test(r.url().split("?")[0]))) errors.push(`http ${r.status()}: ${r.url()}`); });
    const t0 = Date.now();
    await page.goto(`${URL0}?static=${encodeURIComponent(DATA)}&ds=${id}`);
    try {
      await page.waitForFunction((id) => { const sv = window.__spatialViewer; return !!(sv && sv.state && sv.state.dsId === id && sv.state.m) && document.getElementById("status").hidden; }, id, { timeout: 180000 });
      await page.waitForFunction(() => { const S = window.__spatialViewer.state; return !(S.m.expr || []).length || (S.exprKey && S.geneIndex.length > 0); }, null, { timeout: 180000 });
    } catch (e) {
      check(false, `${id}: the page opens the dataset`, await page.evaluate(() => document.getElementById("status").textContent));
      await page.close(); continue;
    }
    const openS = (Date.now() - t0) / 1000;
    await page.waitForTimeout(2500);
    const lit = () => page.evaluate(() => {
      const sv = window.__spatialViewer; if (!sv.deck || !sv.deck.device) return -1;
      const gl = sv.deck.device.gl; sv.deck.redraw("test");
      const w = gl.drawingBufferWidth, h = gl.drawingBufferHeight, px = new Uint8Array(w * h * 4);
      gl.bindFramebuffer(gl.FRAMEBUFFER, null); gl.readPixels(0, 0, w, h, gl.RGBA, gl.UNSIGNED_BYTE, px);
      let n = 0; for (let i = 0; i < px.length; i += 4) if (px[i] > 8 || px[i + 1] > 8 || px[i + 2] > 8) n++;
      return n / (w * h);
    });
    const st = await page.evaluate(() => {
      const sv = window.__spatialViewer, S = sv.state, e = (S.m.expr || []).find((x) => x.key === S.exprKey) || {};
      const sb = document.getElementById("scalebar"), bar = sb.querySelector(".bar"), lbl = sb.querySelector(".lbl").textContent;
      return { exprKey: S.exprKey || null, layers: (S.m.expr || []).map((x) => x.key), genes: (S.geneIndex || []).length,
               pixelSizeUm: S.m.pixelSizeUm || null, zoom: S.viewState.zoom, scaleLabel: lbl, scaleBarPx: parseFloat(bar.style.width),
               scaleHidden: sb.hidden, chanHidden: document.getElementById("chanList").hidden, channels: (S.channels || []).map((c) => c.name),
               user: ((e.umap || {}).clusterings || []).filter((c) => c.source === "user").map((c) => c.name),
               clusterings: ((e.umap || {}).clusterings || []).length,
               outlines: !!(sv.outlines && sv.outlines.available()), molecules: !!S.m.molecules,
               image: [S.m.image.width, S.m.image.height], spot: e.spotDiameterPx || null, nBins: e.nBins || null };
    });
    st.openSeconds = openS;
    st.lit = await lit();
    check(st.lit > 0.02, `${id}: tiles drawn (${(100 * st.lit).toFixed(1)} % of the map lit), opened in ${openS.toFixed(1)} s`, st.lit);
    // the scale bar: its label's length / its width in screen pixels = µm (or image px) per screen pixel
    const m = /^([\d.]+) (µm|mm|px)$/.exec(st.scaleLabel || "");
    if (check(!st.scaleHidden && m, `${id}: a scale bar (${st.scaleLabel})`, st.scaleLabel)) {
      const len = parseFloat(m[1]) * (m[2] === "mm" ? 1000 : 1);
      const perScreen = Math.pow(2, -st.zoom) * (st.pixelSizeUm || 1);
      check((m[2] === "px") === !st.pixelSizeUm && Math.abs(len / st.scaleBarPx - perScreen) / perScreen < 0.01,
        `${id}: the scale bar is ${st.pixelSizeUm ? "in µm and matches the pixel size " + st.pixelSizeUm.toFixed(4) : "in image pixels (no pixel size)"}`,
        { label: st.scaleLabel, barPx: st.scaleBarPx, perScreen });
    }
    await shot(page, `${id}.png`);
    if (st.exprKey) {
      check(st.genes > 0, `${id}: layer ${st.exprKey} with ${st.genes} genes`, st);
      const gene = await page.evaluate(() => [...window.__spatialViewer.state.geneIndex].sort((a, b) => b.t - a.t)[0].n);
      // an image's marker squares keep the channel list first: their panel starts folded; open it as a person would
      if (await page.evaluate(() => document.getElementById("panelGenes").classList.contains("closed"))) await page.click("#panelGenes h2");
      await page.fill("#geneSearch", gene); await page.keyboard.press("Enter");
      let shown = false;
      try {
        await page.waitForFunction(() => window.__spatialViewer.state.genes.length === 1 && window.__spatialViewer.state.genes[0].n > 0, null, { timeout: 60000 });
        shown = true;
      } catch (e) { /* reported below */ }
      st.gene = gene;
      check(shown, `${id}: a gene from the search shows on the slide (${gene})`);
      await page.waitForTimeout(1500);
      await shot(page, `${id}_gene.png`);
    }
    if (st.outlines) {
      // zoom to the densest part: the median row position, until the median cell is 40 screen pixels across
      const v = await page.evaluate(async () => {
        const S = window.__spatialViewer.state, e = S.m.expr.find((x) => x.key === S.exprKey);
        const pos = (S.umap && S.umap.pos) || await S.src.positions(S.exprKey), n = pos.length / 2;   // the UMAP's copy is fetched lazily
        const xs = [], ys = []; for (let r = 0; r < n; r++) if (Number.isFinite(pos[2 * r])) { xs.push(pos[2 * r]); ys.push(pos[2 * r + 1]); }
        xs.sort((a, b) => a - b); ys.sort((a, b) => a - b);
        return [xs[xs.length >> 1], ys[ys.length >> 1], Math.log2(40 / e.spotDiameterPx)];
      });
      await page.evaluate(([x, y, zoom]) => { const sv = window.__spatialViewer; sv.deck.props.onViewStateChange({ viewState: { ...sv.state.viewState, target: [x, y, 0], zoom }, viewId: "ortho" }); }, v);
      await page.evaluate(() => { window.__idle = 0; });
      await page.waitForFunction(() => { const s = window.__spatialViewer.outlines.stats(); if (s.busy) { window.__idle = 0; return false; } window.__idle = (window.__idle || 0) + 1; return window.__idle >= 4; }, null, { timeout: 60000, polling: 100 });
      await page.waitForTimeout(800);
      const os = await page.evaluate(() => window.__spatialViewer.outlines.stats());
      st.outlineStats = os;
      check((os.shown || 0) > 0 && (os.tiles || 0) > 0, `${id}: outlines drawn when zoomed in (${os.shown} cells in ${os.tiles} tiles)`, os);
      await shot(page, `${id}_outlines.png`);
    }
    check(errors.length === 0, `${id}: no console error, page error or failed request`, errors.slice(0, 5));
    console.log("RESULT " + JSON.stringify({ id, ...st, errors: errors.slice(0, 5) }));
    await page.close();
  }
  await browser.close(); data.close(); pageSrv.close();
  const bad = results.filter((r) => !r.ok);
  console.log(bad.length ? `\nFAILED: ${bad.length} of ${results.length} checks` : `\nPASSED: ${results.length} checks`);
  process.exit(bad.length ? 1 : 0);
})().catch((e) => { console.error(e); process.exit(2); });
