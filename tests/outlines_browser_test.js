#!/usr/bin/env node
// Cell outlines in the browser: opens the synthetic demo exports in Playwright's Chromium and
// checks what outlines.js draws for a cells layer with segmentation boundaries: nothing fetched
// until zoomed in, outlines taking over from the dots, the colouring (gene, cluster, selection of
// a cluster, quality filter), hover and click answering as the dot does, nucleus outlines, the
// Dots choice, Fit, a dataset without outlines, and that the page logs no error and no
// Content-Security-Policy violation that it does not also log without outlines.
//
//   node tests/outlines_browser_test.js                   (needs demo/export: run tests/smoke_test.py first)
//   node tests/outlines_browser_test.js --export <dir> --page <dir> --shots <dir> --ports 8866,8867 --gpu
//
// Playwright is not a dependency of this repo. Install it anywhere (npm i playwright; npx
// playwright install chromium) and point NODE_PATH or --playwright at its node_modules.
// Without it the test says so and exits 0. It draws with software rendering unless --gpu is given.
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
const FROM_PX = 24, FULL_PX = 30;

let pw = null;
for (const base of [arg("--playwright", ""), ...(process.env.NODE_PATH || "").split(path.delimiter), path.join(REPO, "node_modules")].filter(Boolean)) {
  try { pw = require(path.join(base, "playwright")); break; } catch (e) { /* next */ }
}
if (!pw) { try { pw = require("playwright"); } catch (e) { /* none */ } }
if (!pw) { console.log("SKIPPED: playwright is not installed (see the top of this file)"); process.exit(0); }
if (!fs.existsSync(path.join(EXPORT, "demo_xenium", "expr", "cells", "shapes", "cell", "index.json"))) {
  console.log(`SKIPPED: no demo export with outlines at ${EXPORT} (run tests/smoke_test.py first, or pass --export)`); process.exit(0);
}

const TYPES = { ".html": "text/html", ".js": "text/javascript", ".css": "text/css", ".json": "application/json", ".png": "image/png", ".jpg": "image/jpeg", ".bin": "application/octet-stream" };
const asked = [];
function serve(root, port, log) {
  return new Promise((resolve) => {
    const srv = http.createServer((req, res) => {
      const rel = decodeURIComponent(req.url.split("?")[0]);
      if (log) asked.push(rel);
      // config.js is each deployment's own file: the test brings its own. The demo sections are small (a cell is
      // ~16 px across with the whole section in view), so the hand-over is moved up to 24 - 30 px for the test.
      if (rel === "/config.js") { res.writeHead(200, { "Content-Type": "text/javascript" }); res.end(`window.SpatialViewerConfig = { tokenEndpoint: "", outlineFromPx: ${FROM_PX}, outlineFullPx: ${FULL_PX} };`); return; }
      const file = path.join(root, rel === "/" ? "index.html" : rel);
      if (req.method === "OPTIONS") { preflight(req, res); return; }
      if (!file.startsWith(root) || !fs.existsSync(file) || fs.statSync(file).isDirectory()) { res.writeHead(404, { "Access-Control-Allow-Origin": "*" }); res.end("not found"); return; }
      sendFile(req, res, file, TYPES[path.extname(file)]);            // with Range: genes are read out of the gene packs by range
    });
    srv.listen(port || 0, "127.0.0.1", () => resolve(srv));
  });
}

const results = [];
const check = (cond, what, got) => { results.push({ ok: !!cond, what }); console.log(`  ${cond ? "ok  " : "FAIL"} ${what}${cond || got === undefined ? "" : "  got " + JSON.stringify(got)}`); return !!cond; };

(async () => {
  const data = await serve(EXPORT, PORTS[0], true), pageSrv = await serve(PAGE, PORTS[1], false);
  const DATA = `http://127.0.0.1:${data.address().port}`, URL0 = `http://127.0.0.1:${pageSrv.address().port}/index.html`;
  const browser = await pw.chromium.launch({ args: GPU ? ["--use-angle=d3d11", "--ignore-gpu-blocklist", "--enable-gpu"] : ["--use-angle=swiftshader", "--enable-unsafe-swiftshader", "--ignore-gpu-blocklist"] });
  const ctx = await browser.newContext({ viewport: { width: 1366, height: 768 } });
  const errors = [];
  const open = async (ds) => {
    const page = await ctx.newPage();
    page.on("console", (m) => { if (m.type() === "error") errors.push(`${ds}: console: ${m.text()}`); });
    page.on("pageerror", (e) => errors.push(`${ds}: pageerror: ${e.message}`));
    page.on("response", (r) => { if (r.status() >= 400 && !/favicon/.test(r.url())) errors.push(`${ds}: http ${r.status()}: ${r.url()}`); });
    await page.addInitScript(() => { window.__csp = []; document.addEventListener("securitypolicyviolation", (e) => window.__csp.push(`${e.violatedDirective} ${e.blockedURI}`)); });
    await page.goto(`${URL0}?static=${encodeURIComponent(DATA)}&ds=${ds}`);
    await page.waitForFunction((id) => { const sv = window.__spatialViewer; return !!(sv && sv.state && sv.state.dsId === id && sv.state.m && sv.outlines) && document.getElementById("status").hidden; }, ds, { timeout: 120000 });
    await page.waitForFunction(() => { const S = window.__spatialViewer.state; return !(S.m.expr || []).length || (S.exprKey && S.geneIndex.length > 0); }, null, { timeout: 120000 });
    await page.waitForTimeout(800);
    return page;
  };
  const view = (page, x, y, zoom) => page.evaluate(([x, y, zoom]) => {
    const sv = window.__spatialViewer; sv.deck.props.onViewStateChange({ viewState: { ...sv.state.viewState, target: [x, y, 0], zoom }, viewId: "ortho" });
  }, [x, y, zoom]);
  const settled = async (page) => {
    await page.evaluate(() => { window.__idle = 0; });
    await page.waitForFunction(() => { const s = window.__spatialViewer.outlines.stats(); if (s.busy) { window.__idle = 0; return false; } window.__idle = (window.__idle || 0) + 1; return window.__idle >= 4; }, null, { timeout: 60000, polling: 100 });
    await page.waitForTimeout(300);
  };
  const st = (page) => page.evaluate(() => window.__spatialViewer.outlines.stats());
  // the zoom at which the layer's median cell is `px` screen pixels across, at the image's centre
  const zoomIn = async (page, px = 40) => {
    const v = await page.evaluate((px) => { const S = window.__spatialViewer.state, e = S.m.expr.find((x) => x.key === S.exprKey), i = S.m.image; return [i.width / 2, i.height / 2, Math.log2(px / e.spotDiameterPx)]; }, px);
    await view(page, v[0], v[1], v[2]); await settled(page);
  };
  const shot = async (page, name) => { if (SHOTS) { fs.mkdirSync(SHOTS, { recursive: true }); await page.screenshot({ path: path.join(SHOTS, name) }); } };
  // rows in view with the screen position of their centre
  const cells = (page, n) => page.evaluate((n) => {
    const sv = window.__spatialViewer, O = sv.outlines._debug(), vp = sv.deck.getViewports()[0], pos = sv.state.umap.pos, out = [];
    for (const T of O.sets[0].tiles.values()) for (let i = 0; i < T.ringRow.length; i += 3) {
      const r = T.ringRow[i], p = vp.project([pos[2 * r], pos[2 * r + 1]]);
      if (p[0] > 20 && p[1] > 20 && p[0] < vp.width - 20 && p[1] < vp.height - 20) out.push([r, p[0], p[1]]);
    }
    return out.slice(0, n);
  }, n);

  // ---- a dataset with no cells layer: the baseline for console and CSP messages
  console.log("demo_xenium_pq (no cells layer)");
  let page = await open("demo_xenium_pq");
  check(await page.evaluate(() => document.getElementById("rowCellShapes").hidden && !window.__spatialViewer.outlines.available()), "no cells layer: the Cell shapes control is hidden");
  const cspBase = [...new Set(await page.evaluate(() => window.__csp))].sort();
  await page.close();

  // ---- the Xenium fixture: cells + nuclei, clusters
  console.log("demo_xenium (cell and nucleus boundaries)");
  asked.length = 0;
  page = await open("demo_xenium");
  check(!asked.some((u) => u.includes("/shapes/")), "the first, zoomed-out view asks for nothing under shapes/", asked.filter((u) => u.includes("/shapes/")));
  check(await page.evaluate(() => !document.getElementById("rowCellShapes").hidden && document.getElementById("cellShapes").value === "outlines"), "the Cell shapes control shows, on Outlines");
  let s = await st(page);
  check(s.alpha === 0 && s.dot === 1 && s.layers === 0, "whole section in view: dots", { alpha: s.alpha, dot: s.dot });
  await zoomIn(page, 27); await page.waitForTimeout(300);
  s = await st(page);
  check(s.alpha > 0.1 && s.alpha < 0.9 && Math.abs(s.alpha + s.dot - 1) < 1e-6, "between the two sizes: outlines and dots cross-faded", { alpha: s.alpha, dot: s.dot });
  await zoomIn(page);
  s = await st(page);
  check(s.alpha === 1 && s.dot === 0 && s.shown > 100, "zoomed in: outlines drawn, dots handed over", { alpha: s.alpha, dot: s.dot, shown: s.shown });
  check(asked.filter((u) => u.includes("/shapes/cell/")).length >= 2 && !asked.some((u) => u.includes("/shapes/nucleus/")), "only the cell set's index and tiles were fetched");
  await shot(page, "outlines_xenium.png");
  const gene = await page.evaluate(() => [...window.__spatialViewer.state.geneIndex].sort((a, b) => b.t - a.t)[0].n);
  await page.fill("#geneSearch", gene); await page.keyboard.press("Enter");
  await page.waitForFunction(() => window.__spatialViewer.state.genes.length === 1 && window.__spatialViewer.state.genes[0].n > 0, null, { timeout: 60000 });
  await settled(page);
  let c = await page.evaluate(() => { const sv = window.__spatialViewer, col = sv.outlines._debug().col, N = sv.state.umap.pos.length / 2; let f = 0; for (let r = 0; r < N; r++) if (col.fill[4 * r + 3] > 0) f++; return { filled: f, withGene: sv.state.genes[0].n }; });
  check(c.filled === c.withGene && c.filled > 0, "colour by gene: the cells holding the gene are filled, the others are not", c);
  await shot(page, "outlines_xenium_gene.png");
  const hasClusters = await page.evaluate(() => { const S = window.__spatialViewer.state, e = S.m.expr.find((x) => x.key === S.exprKey); return !!(e.umap && (e.umap.clusterings || []).length); });
  if (hasClusters) {
    await page.click("#tabClusters"); await page.selectOption("#colorBy", "cluster");
    await page.waitForFunction(() => { const U = window.__spatialViewer.state.umap; return U && U.cl && U.colorBy === "cluster"; }, null, { timeout: 60000 });
    await settled(page);
    c = await page.evaluate(() => { const sv = window.__spatialViewer, col = sv.outlines._debug().col, U = sv.state.umap, N = U.pos.length / 2; let f = 0, k = 0; for (let r = 0; r < N; r++) { if (col.fill[4 * r + 3] > 0) f++; if (U.cl[r] !== 0xffff) k++; } return { filled: f, clustered: k }; });
    check(c.filled === c.clustered && c.filled > 0, "colour by cluster: every clustered cell is filled", c);
    const pts = await cells(page, 60);
    const got = await page.evaluate((pts) => pts.map(([r, x, y]) => { const sv = window.__spatialViewer, info = sv.deck.pickObject({ x, y, radius: 0 }); return info ? sv.outlines.rowAt(info) : -1; }), pts);
    const own = got.filter((r, i) => r === pts[i][0]).length;
    check(own >= 0.95 * pts.length && pts.length >= 30, "picking at a cell's centre gives that cell's row", { own, of: pts.length });
    const box = await page.evaluate(() => { const r = document.getElementById("viewport").getBoundingClientRect(); return [r.left, r.top]; });
    const t = pts[Math.floor(pts.length / 2)], want = await page.evaluate((r) => window.__spatialViewer.state.umap.cl[r], t[0]);
    await page.mouse.move(box[0] + t[1], box[1] + t[2]); await page.waitForTimeout(500);
    const tip = await page.evaluate(() => { const e = document.querySelector(".deck-tooltip"); return e && e.style.display !== "none" ? e.textContent : ""; });
    check(/Cluster|\w/.test(tip) && tip.length > 0, "hovering an outline shows the cell's tooltip", tip);
    await page.mouse.click(box[0] + t[1], box[1] + t[2]); await page.waitForTimeout(700);
    check((await page.evaluate(() => window.__spatialViewer.state.umap.pick)) === want, "clicking an outline highlights its cluster, as its dot does", want);
    c = await page.evaluate(() => { const sv = window.__spatialViewer, col = sv.outlines._debug().col, U = sv.state.umap, N = U.pos.length / 2; let f = 0, k = 0; for (let r = 0; r < N; r++) { if (col.fill[4 * r + 3] > 0) f++; if (U.cl[r] === U.pick) k++; } return { lit: f, inCluster: k, over: col.over }; });
    check(c.lit === c.inCluster && c.over, "a highlighted cluster: only its cells are filled, above the dimmed slide", c);
    await shot(page, "outlines_xenium_cluster_picked.png");
    await page.evaluate(() => { const b = document.getElementById("selClear"); if (!b.disabled) b.click(); });
  } else console.log("  (no clusters in this export: cluster checks skipped)");
  if (await page.evaluate(() => !!(window.__spatialViewer.binQc && window.__spatialViewer.binQc.available()))) {
    const q = await page.evaluate(() => { const sv = window.__spatialViewer; sv.binQc.setState({ thresholds: { total: { min: 45, max: null } } }); return sv.binQc.removed(); });
    await page.waitForTimeout(600);
    c = await page.evaluate(() => { const sv = window.__spatialViewer, col = sv.outlines._debug().col, keep = sv.binQc.keepMask(); let g = 0; if (keep) for (let r = 0; r < keep.length; r++) if (!keep[r] && col.fill[4 * r] === 88 && col.fill[4 * r + 3] === 200) g++; return g; });
    check(q > 0 && c === q, "quality filter: the removed cells are grey", { removed: q, grey: c });
    await page.evaluate(() => window.__spatialViewer.binQc.setState(null));
  }
  await page.evaluate(() => window.__spatialViewer.outlines.setNuclei(true)); await settled(page);
  s = await st(page);
  check(s.sets.length === 2 && s.sets[1].polygons > 100 && asked.some((u) => u.includes("/shapes/nucleus/")), "Nucleus outlines: the second set is fetched and drawn", s.sets[1]);
  await shot(page, "outlines_xenium_nuclei.png");
  await page.evaluate(() => window.__spatialViewer.outlines.setNuclei(false));
  await page.click("#tabView"); await page.selectOption("#cellShapes", "dots"); await page.waitForTimeout(400);
  s = await st(page);
  check(s.alpha === 0 && s.dot === 1 && s.layers === 0, "Cell shapes: Dots brings the dots back", { alpha: s.alpha, dot: s.dot });
  await page.selectOption("#cellShapes", "outlines"); await settled(page);
  check((await st(page)).alpha === 1, "Cell shapes: Outlines again");
  await page.click("#btnFit"); await page.waitForTimeout(600);
  s = await st(page);
  check(s.alpha === 0 && s.dot === 1 && s.layers === 0, "Fit: back to dots, no outline layer left on screen", { alpha: s.alpha, dot: s.dot, layers: s.layers });
  const csp = [...new Set(await page.evaluate(() => window.__csp))].sort();
  check(JSON.stringify(csp) === JSON.stringify(cspBase), "no Content-Security-Policy violation beyond those of a page without outlines", { withOutlines: csp, without: cspBase });
  await page.close();

  // ---- the other platforms' fixtures: outlines appear
  for (const ds of ["demo_cosmx", "demo_merscope", "demo_stereoseq"]) {
    if (!fs.existsSync(path.join(EXPORT, ds, "manifest.json"))) { console.log(`${ds}: not in this export, skipped`); continue; }
    console.log(ds);
    page = await open(ds);
    if ((await page.evaluate(() => window.__spatialViewer.state.exprKey)) !== "cells") {
      await page.selectOption("#exprSelect", "cells");
      await page.waitForFunction(() => window.__spatialViewer.state.exprKey === "cells" && window.__spatialViewer.outlines.available(), null, { timeout: 60000 });
      await page.waitForTimeout(800);
    }
    await zoomIn(page);
    s = await st(page);
    check(s.alpha === 1 && s.shown > 100, `${ds}: outlines drawn when zoomed in`, { alpha: s.alpha, shown: s.shown, note: s.note });
    await shot(page, `outlines_${ds}.png`);
    await page.close();
  }
  // ---- a layer of bins: no outlines, nothing asked for
  if (fs.existsSync(path.join(EXPORT, "demo_visiumhd", "manifest.json"))) {
    console.log("demo_visiumhd (bins)");
    asked.length = 0;
    page = await open("demo_visiumhd");
    await page.evaluate(() => { const sv = window.__spatialViewer, i = sv.state.m.image; sv.deck.props.onViewStateChange({ viewState: { ...sv.state.viewState, target: [i.width / 2, i.height / 2, 0], zoom: 1 }, viewId: "ortho" }); });
    await page.waitForTimeout(800);
    check(await page.evaluate(() => document.getElementById("rowCellShapes").hidden) && !asked.some((u) => u.includes("/shapes/")), "bins: no Cell shapes control, nothing asked for under shapes/");
    await page.close();
  }
  check(errors.length === 0, "no console error, page error or failed request", errors.slice(0, 5));

  await browser.close(); data.close(); pageSrv.close();
  const bad = results.filter((r) => !r.ok);
  console.log(bad.length ? `\nFAILED: ${bad.length} of ${results.length} checks` : `\nPASSED: ${results.length} checks`);
  process.exit(bad.length ? 1 : 0);
})().catch((e) => { console.error(e); process.exit(2); });
