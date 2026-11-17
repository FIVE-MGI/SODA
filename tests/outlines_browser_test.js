#!/usr/bin/env node
// Cell outlines in the browser: opens the synthetic demo exports in Playwright's Chromium and
// checks what outlines.js draws for a cells layer with segmentation boundaries: nothing fetched
// until zoomed in, outlines taking over from the dots, the colouring (gene, cluster, selection of
// a cluster, quality filter), hover and click answering as the dot does, nucleus outlines, the
// Dots choice, Fit, a dataset without outlines, and that the page logs no error and no
// Content-Security-Policy violation that it does not also log without outlines.
// Then, in each browser of --browsers (chromium, webkit, firefox), what draws a subset of the cells, on demo_xenium_nt: Colour by
// a value (one value, two genes, a gene set score), the region click-filter, Draw a region (and a picked shape), the
// Dots choice bringing the dots back, Save image (its PNG and legend) and a link to the view, checked by pixel where it
// matters: a selected cell's polygon is coloured inside, also where its dot would not reach.
//
//   node tests/outlines_browser_test.js                   (needs demo/export: run tests/smoke_test.py first)
//   node tests/outlines_browser_test.js --export <dir> --page <dir> --shots <dir> --ports 8866,8867 --gpu --browsers chromium,webkit
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
// Firefox's console lines that are not the page's faults (see tests/browsers_test.js): the vendor/ tag's SHA-384
// checked against a 404 when there is no vendor/, and deck.gl's WebAssembly probe refused by the page's CSP
const FIREFOX_QUIET = [/integrity attribute match the content of the subresource at .[^ ]*\/vendor\//,
  /Content-Security-Policy: .*blocked WebAssembly .*Missing 'wasm-unsafe-eval'/];
// the browsers the checks of what draws a subset of the cells run in (Draw a region, the region click-filter, Colour by
// a value, Save image, a link to the view); the checks before them run in Chromium only
const BROWSERS = arg("--browsers", "chromium").split(",").map((s) => s.trim()).filter(Boolean);
const GPU_ARGS = ["--use-angle=d3d11", "--ignore-gpu-blocklist", "--enable-gpu"];
const SOFT_ARGS = ["--use-angle=swiftshader", "--enable-unsafe-swiftshader", "--ignore-gpu-blocklist"];
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

// ---- in the page: pixels of the map, and places inside the cells to read them at ------------------------
// The map's pixels at points given in CSS pixels of the viewport, read from deck's canvas (it keeps its last frame),
// over black as the page shows it (the canvas itself is see-through where little is drawn).
function pagePixels(pts) {
  const sv = window.__spatialViewer, cv = sv.deck.getCanvas ? sv.deck.getCanvas() : sv.deck.canvas, k = cv.width / cv.clientWidth;
  const c = document.createElement("canvas"); c.width = cv.width; c.height = cv.height;
  const g = c.getContext("2d"); g.fillStyle = "#000"; g.fillRect(0, 0, c.width, c.height); g.drawImage(cv, 0, 0);
  return pts.map(([x, y]) => Array.from(g.getImageData(Math.round(x * k), Math.round(y * k), 1, 1).data));
}
// Cells drawn as outlines whose row passes `want` (a list of rows, or null = any): for each, its centre on screen and a
// point inside its polygon away from the centre (toward its farthest vertex, at least 3 screen px from the outline);
// beyond: that point lies outside the cell's own dot and every other cell's dot (a dot cannot be what colours it).
function pageCellSamples(arg) {
  const [want, n] = arg;
  const sv = window.__spatialViewer, O = sv.outlines._debug(), vp = sv.deck.getViewports()[0], S = sv.state, pos = S.umap.pos;
  const e = S.m.expr.find((x) => x.key === S.exprKey), rad = Math.max(0.5, (e.spotDiameterPx / 2) * S.pointScale), scale = Math.pow(2, S.viewState.zoom);
  const ok = want ? new Set(want) : null, out = [];
  const box = document.getElementById("viewport").getBoundingClientRect();
  const segD = (px, py, ax, ay, bx, by) => { const dx = bx - ax, dy = by - ay, L = dx * dx + dy * dy, t = L ? Math.max(0, Math.min(1, ((px - ax) * dx + (py - ay) * dy) / L)) : 0; return Math.hypot(px - ax - t * dx, py - ay - t * dy); };
  for (const T of O.sets[0].tiles.values()) {
    const P = T.fill.attributes.getPolygon.value, st = T.fill.startIndices;
    for (let p = 0; p < T.n; p++) {
      const r = T.ringRow[p];
      if (ok && !ok.has(r)) continue;
      const cx = pos[2 * r], cy = pos[2 * r + 1], c = vp.project([cx, cy]);
      if (c[0] < 0.15 * box.width || c[0] > 0.85 * box.width || c[1] < 0.15 * box.height || c[1] > 0.8 * box.height) continue;
      const ring = []; for (let v = st[p]; v < st[p + 1]; v++) ring.push([P[2 * v], P[2 * v + 1]]);
      let far = ring[0], fd = 0;
      for (const q of ring) { const d = Math.hypot(q[0] - cx, q[1] - cy); if (d > fd) { fd = d; far = q; } }
      let hit = null;
      for (let f = 0.85; f >= 0.45 && !hit; f -= 0.05) {
        const x = cx + f * (far[0] - cx), y = cy + f * (far[1] - cy);
        let inside = false;
        for (let i = 0, j = ring.length - 1; i < ring.length; j = i++) { const [xi, yi] = ring[i], [xj, yj] = ring[j]; if ((yi > y) !== (yj > y) && x < xi + ((xj - xi) * (y - yi)) / (yj - yi)) inside = !inside; }
        if (!inside) continue;
        let edge = Infinity; for (let i = 0; i + 1 < ring.length; i++) edge = Math.min(edge, segD(x, y, ring[i][0], ring[i][1], ring[i + 1][0], ring[i + 1][1]));
        if (edge * scale < 3 || Math.hypot(x - cx, y - cy) * scale < 4) continue;
        hit = [x, y];
      }
      if (!hit) continue;
      let beyond = Math.hypot(hit[0] - cx, hit[1] - cy) > 1.15 * rad;
      if (beyond) for (let q = 0; q < pos.length / 2; q++) if (q !== r && Math.hypot(pos[2 * q] - hit[0], pos[2 * q + 1] - hit[1]) < 1.15 * rad) { beyond = false; break; }
      const s = vp.project(hit);
      out.push({ r, c: [c[0], c[1]], s: [s[0], s[1]], beyond });
      if (out.length >= n) return out;
    }
  }
  return out;
}
const reddish = (p) => p[0] >= 40 && p[0] > p[1] + 25 && p[0] > p[2] + 12;          // Draw a region's first colour, #e6194b
const share = (list, f) => (list.length ? list.filter(f).length / list.length : 0);

// Draw a region, the region click-filter, Colour by a value, Save image and a link to the view, with the cells drawn as
// their outlines, in one browser (chromium, webkit or firefox): the selected cells' polygons are filled, not only dotted.
async function selectionChecks(B, browser, DATA, URL0) {
  const DS = "demo_xenium_nt";
  if (!fs.existsSync(path.join(EXPORT, DS, "expr", "cells", "shapes", "cell", "index.json"))) { console.log(`${B}: no ${DS} with outlines, selection checks skipped`); return; }
  console.log(`${B}: ${DS} (what draws a subset of the cells)`);
  const ctx = await browser.newContext({ viewport: { width: 1366, height: 768 }, serviceWorkers: "block", acceptDownloads: true });
  const errors = [];
  const watch = (page) => {
    page.on("console", (m) => { if (m.type() === "error" && !/\/vendor\//.test((m.location() || {}).url || "") && !(B === "firefox" && FIREFOX_QUIET.some((re) => re.test(m.text())))) errors.push(`console: ${m.text()}`); });
    page.on("pageerror", (e) => errors.push(`pageerror: ${e.message}`));
    page.on("response", (r) => { if (r.status() >= 400 && !/favicon|\/vendor\//.test(r.url())) errors.push(`http ${r.status()}: ${r.url()}`); });
  };
  const ready = async (page) => {
    await page.waitForFunction((id) => { const sv = window.__spatialViewer; return !!(sv && sv.state && sv.state.dsId === id && sv.state.m && sv.outlines) && document.getElementById("status").hidden; }, DS, { timeout: 120000 });
    await page.waitForFunction(() => { const S = window.__spatialViewer.state; return S.exprKey && S.geneIndex.length > 0 && S.lasso && S.lasso.pos && S.umap && !(window.SodaViewLink && window.SodaViewLink.pending()); }, null, { timeout: 120000 });
    await page.waitForTimeout(600);
  };
  const settled = async (page) => {
    await page.evaluate(() => { window.__idle = 0; });
    await page.waitForFunction(() => {
      const sv = window.__spatialViewer, s = sv.outlines.stats(), O = sv.outlines._debug();
      // busy: tiles on their way; outlines drawn before the rows' positions are in: not yet coloured
      if (s.busy || (s.alpha > 0 && !(sv.state.umap && sv.state.umap.pos && O.col && O.col.pos))) { window.__idle = 0; return false; }
      window.__idle = (window.__idle || 0) + 1; return window.__idle >= 4;
    }, null, { timeout: 60000, polling: 100 });
    await page.evaluate(() => new Promise((r) => requestAnimationFrame(() => requestAnimationFrame(r))));
    await page.waitForTimeout(250);
  };
  const shotB = async (page, name) => { if (SHOTS) { fs.mkdirSync(SHOTS, { recursive: true }); await page.screenshot({ path: path.join(SHOTS, `${B}_${name}`) }); } };
  const st = (page) => page.evaluate(() => window.__spatialViewer.outlines.stats());
  const px = (page, pts) => page.evaluate(pagePixels, pts);
  const samples = (page, want, n = 40) => page.evaluate(pageCellSamples, [want, n]);
  const layersNow = (page) => page.evaluate(() => window.__spatialViewer.deck.layerManager.getLayers().filter((l) => !l.isComposite || /^lasso-bins/.test(l.id)).map((l) => ({ id: l.id, op: l.props.opacity, visible: l.props.visible })));
  // a control in a folded section or another tab: its value set and its own handler run, as a person's pick would
  const setVal = (page, id, v, ev = "change") => page.evaluate(([id, v, ev]) => { const el = document.getElementById(id); el.value = v; el.dispatchEvent(new Event(ev, { bubbles: true })); }, [id, v, ev]);
  const colOf = (page) => page.evaluate(() => { const c = window.__spatialViewer.outlines._debug().col; return c ? { fill: Array.from(c.fill), line: Array.from(c.line), opacity: c.opacity } : null; });

  let page = await ctx.newPage(); watch(page);
  await page.goto(`${URL0}?static=${encodeURIComponent(DATA)}&ds=${DS}`);
  await ready(page);
  // the image off: the map is black where nothing is drawn, so a coloured pixel is a cell's colour
  await page.evaluate(() => { if (window.__spatialViewer.state.imgOn) document.getElementById("btnImage").click(); });
  const zv = await page.evaluate(() => { const S = window.__spatialViewer.state, e = S.m.expr.find((x) => x.key === S.exprKey), i = S.m.image; return [i.width / 2, i.height / 2, Math.log2(40 / e.spotDiameterPx)]; });
  await page.evaluate(([x, y, z]) => { const sv = window.__spatialViewer; sv.deck.props.onViewStateChange({ viewState: { ...sv.state.viewState, target: [x, y, 0], zoom: z }, viewId: "ortho" }); }, zv);
  // dots at half size: a point inside a cell's polygon but outside its dot then tells the two apart
  await setVal(page, "pointScale", "0.5", "input");
  await settled(page);
  let s = await st(page);
  check(s.alpha === 1 && s.dot === 0, `${B}: zoomed in, the cells are outlines`, { alpha: s.alpha, dot: s.dot });

  // ---- Colour by a value (binpaint.js): the outlines take its colours, at its dots' one alpha
  const N = await page.evaluate(() => window.__spatialViewer.state.umap.pos.length / 2);
  const paintSame = async (what) => {
    await page.waitForFunction(() => { const P = window.SodaPaint._state(); return P.colors && window.SodaPaint.active(); }, null, { timeout: 60000 });
    await settled(page);
    const r = await page.evaluate(() => {
      const P = window.SodaPaint._state(), f = window.__spatialViewer.outlines._debug().col.fill, pc = P.colors;
      let same = 0, coloured = 0, bad = 0;
      for (let r = 0; r < pc.length / 4; r++) {
        const k = 4 * r;
        if (pc[k + 3]) { coloured++; if (f[k] === pc[k] && f[k + 1] === pc[k + 1] && f[k + 2] === pc[k + 2] && f[k + 3] === pc[k + 3]) same++; }
        else if (f[k + 3]) bad++;
      }
      return { coloured, same, bad };
    });
    check(r.coloured > 0 && r.same === r.coloured && r.bad === 0, `${B}: ${what}: every coloured cell's outline is filled in its dot's colour and alpha, no other`, r);
  };
  if (await page.evaluate(() => !document.querySelector("#paintMode option[value=value]").disabled)) {
    await setVal(page, "paintMode", "value");
    await paintSame("Colour by a value");
    await shotB(page, "select_paint_value.png");
    if (B === "chromium") {
      const t = (await samples(page, null, 10))[0];
      const box = await page.evaluate(() => { const r = document.getElementById("viewport").getBoundingClientRect(); return [r.left, r.top]; });
      if (t) {
        await page.mouse.move(box[0] + t.c[0], box[1] + t.c[1]); await page.waitForTimeout(500);
        const tip = await page.evaluate(() => { const e = document.querySelector(".deck-tooltip"); return e && e.style.display !== "none" ? e.textContent : ""; });
        const want = await page.evaluate((r) => window.SodaPaint.tip(r), t.r);
        check(!!want && tip.includes(want.split("\n")[0]), `${B}: hovering an outline coloured by a value tells its value, as its dot does`, { tip, want });
        await page.mouse.move(box[0] + 5, box[1] + 5);
      }
    }
  } else console.log(`  ${B}: no per-cell value in this export: Colour by a value skipped`);
  const genes = await page.evaluate(() => [...window.__spatialViewer.state.geneIndex].sort((a, b) => b.t - a.t).slice(0, 2).map((g) => g.n));
  if (genes.length === 2) {
    await setVal(page, "paintMode", "pair");
    for (const [w, g] of [["A", genes[0]], ["B", genes[1]]]) await setVal(page, `paintGene${w}`, g);
    await page.waitForFunction(() => { const P = window.SodaPaint._state(); return P.pd && P.colors; }, null, { timeout: 60000 });
    await paintSame("Two genes (co-expression)");
    await shotB(page, "select_paint_pair.png");
    await setVal(page, "paintMode", "set");
    await setVal(page, "paintSetGenes", genes.join(", "), "input"); await page.evaluate(() => document.getElementById("paintSetRun").click());
    await page.waitForFunction(() => { const P = window.SodaPaint._state(); return P.mode === "set" && P.cur && P.colors; }, null, { timeout: 60000 });
    await paintSame("Gene set score");
  }
  await setVal(page, "paintMode", "off");
  await settled(page);

  // ---- the region click-filter: a region picked with a click cuts the gene density to its bins, and the outlines with it
  const gene = genes[0];
  await page.evaluate(() => { const b = document.getElementById("tabView"); if (b) b.click(); });
  await page.fill("#geneSearch", gene); await page.keyboard.press("Enter");
  await page.waitForFunction(() => window.__spatialViewer.state.genes.length === 1 && window.__spatialViewer.state.genes[0].n > 0, null, { timeout: 60000 });
  await settled(page);
  // two regions, left and right of the image's middle, as region.bin and a region overlay would give them
  await page.evaluate(() => {
    const sv = window.__spatialViewer, S = sv.state, pos = S.umap.pos, n = pos.length / 2, mid = S.m.image.width / 2, h = S.m.image.height;
    const ids = new Uint16Array(n), lookup = new Map();
    for (let r = 0; r < n; r++) { ids[r] = pos[2 * r] < mid ? 1 : 2; lookup.set(`${Math.round(pos[2 * r] * 100)},${Math.round(pos[2 * r + 1] * 100)}`, ids[r]); }
    S.region.meta = [{ id: 1, name: "Left half" }, { id: 2, name: "Right half" }]; S.region.ids = ids; S.region.lookup = lookup;
    const f = (id, x0, x1) => ({ type: "Feature", properties: { name: id === 1 ? "Left half" : "Right half", __region_id: id }, geometry: { type: "Polygon", coordinates: [[[x0, 0], [x1, 0], [x1, h], [x0, h], [x0, 0]]] } });
    S.overlays.push({ key: "test_regions", file: "test_regions.geojson", label: "Test regions", on: true, data: { type: "FeatureCollection", features: [f(1, 0, mid), f(2, mid, 2 * mid)] } });
    sv.render();
  });
  const regionFill = () => page.evaluate(() => {
    const sv = window.__spatialViewer, S = sv.state, f = sv.outlines._debug().col.fill, g = S.genes[0], ids = S.region.ids, pos = S.umap.pos;
    const has = new Set(); for (let i = 0; i < g.n; i++) has.add(`${Math.round(g.positions[2 * i] * 100)},${Math.round(g.positions[2 * i + 1] * 100)}`);
    const o = { sel: S.region.selected, inGene: 0, inFilled: 0, outGene: 0, outFilled: 0 };
    for (let r = 0; r < ids.length; r++) {
      if (!has.has(`${Math.round(pos[2 * r] * 100)},${Math.round(pos[2 * r + 1] * 100)}`)) continue;
      const lit = f[4 * r + 3] > 0;
      if (!S.region.selected || ids[r] === S.region.selected) { o.inGene++; if (lit) o.inFilled++; } else { o.outGene++; if (lit) o.outFilled++; }
    }
    return o;
  });
  await page.evaluate(() => window.__spatialViewer.deck.props.onClick({ object: { properties: { __region_id: 1, name: "Left half" } }, x: 5, y: 5 }, {}));
  await settled(page);
  let rf = await regionFill();
  check(rf.sel === 1 && rf.inGene > 0 && rf.inFilled === rf.inGene && rf.outGene > 0 && rf.outFilled === 0,
    `${B}: a region picked with a click: only its cells with the gene are filled (as the density map is cut to it)`, rf);
  await shotB(page, "select_region_filter.png");
  // something lit puts the outlines over the region outlines: a click on a cell still picks the region under it
  const right = (await samples(page, await page.evaluate(() => { const ids = window.__spatialViewer.state.region.ids, o = []; for (let r = 0; r < ids.length; r++) if (ids[r] === 2) o.push(r); return o; }), 5))[0];
  if (right) {
    await page.evaluate(() => { const U = window.__spatialViewer.state.umap, n = U.pos.length / 2; U.sel = new Uint8Array(n); for (let r = 0; r < n; r += 2) U.sel[r] = 1; U.selN = Math.ceil(n / 2); U.ver++; window.__spatialViewer.render(); });
    await settled(page);
    const over = await page.evaluate(() => window.__spatialViewer.outlines._debug().col.over);
    const box = await page.evaluate(() => { const r = document.getElementById("viewport").getBoundingClientRect(); return [r.left, r.top]; });
    // the pointer first, then the click once the hover's frame is drawn: deck counts a press as a click only when it is
    // short, and with software rendering a frame drawn in between can make it too long
    await page.mouse.move(box[0] + right.c[0], box[1] + right.c[1]); await page.waitForTimeout(800);
    await page.mouse.click(box[0] + right.c[0], box[1] + right.c[1]);
    // deck holds a click back a moment (it might be the first of a double-click)
    await page.waitForFunction(() => window.__spatialViewer.state.region.selected === 2, null, { timeout: 8000 }).catch(() => {});
    check(over && (await page.evaluate(() => window.__spatialViewer.state.region.selected)) === 2, `${B}: with the outlines drawn over the regions (something lit), clicking a cell picks the region under it`, { over });
    await page.evaluate(() => { const U = window.__spatialViewer.state.umap; U.sel = null; U.selN = 0; U.ver++; window.__spatialViewer.render(); });
  }
  await page.evaluate(() => window.__spatialViewer.deck.props.onClick({ x: 5, y: 5 }, {}));
  await settled(page);
  rf = await regionFill();
  check(!rf.sel && rf.inFilled === rf.inGene && rf.inGene > 0, `${B}: the region cleared: every cell with the gene is filled again`, rf);
  // leave no gene and no regions for what follows
  await page.evaluate(() => { const sv = window.__spatialViewer, S = sv.state; S.overlays = S.overlays.filter((o) => o.key !== "test_regions"); S.region.meta = null; S.region.ids = null; S.region.lookup = null; });
  await page.evaluate(() => { const b = document.querySelector("#geneList .chan .x"); if (b) b.click(); window.__spatialViewer.render(); });
  await settled(page);
  check(await page.evaluate(() => window.__spatialViewer.state.genes.length === 0), `${B}: the gene removed again`);

  // ---- Draw a region: a circle over the middle of the map, with the mouse
  const vbox = await page.evaluate(() => { const r = document.getElementById("viewport").getBoundingClientRect(); return { x: r.left, y: r.top, w: r.width, h: r.height }; });
  await page.evaluate(() => { if (window.__spatialViewer.state.lasso.tool !== "circle") document.getElementById("lassoCircle").click(); });
  const cx = vbox.x + vbox.w / 2, cy = vbox.y + vbox.h / 2, rad = Math.min(vbox.w, vbox.h) * 0.22;
  await page.mouse.move(cx, cy); await page.mouse.down();
  for (let i = 1; i <= 8; i++) await page.mouse.move(cx + (rad * i) / 8, cy);
  await page.mouse.up();
  await page.evaluate(() => { if (window.__spatialViewer.state.lasso.tool) document.getElementById("lassoCircle").click(); });
  await page.waitForFunction(() => { const L = window.__spatialViewer.state.lasso; return L.shapes.length === 1 && L.bin && L.bin.n > 0; }, null, { timeout: 30000 });
  await settled(page);
  const inShape = await page.evaluate(() => { const a = window.__spatialViewer.state.lasso.assign, o = []; for (let r = 0; r < a.length; r++) if (a[r] !== 0xffff) o.push(r); return o; });
  let col = await colOf(page);
  const filledIn = inShape.filter((r) => col.fill[4 * r + 3] > 0).length;
  let filledOut = 0; for (let r = 0; r < N; r++) if (col.fill[4 * r + 3] > 0 && !inShape.includes(r)) filledOut++;
  const hlRgb = await page.evaluate(() => { const h = window.__spatialViewer.state.lasso.shapes[0].hex; return [1, 3, 5].map((i) => parseInt(h.slice(i, i + 2), 16)); });
  const wantA = inShape.length ? Math.round(255 * (80 / 255)) : 0;
  check(inShape.length > 5 && filledIn === inShape.length && filledOut === 0 && inShape.every((r) => col.fill[4 * r] === hlRgb[0] && col.fill[4 * r + 3] === wantA),
    `${B}: Draw a region: the outlines of the cells inside the shape are filled in its colour (alpha as its dots), no other`, { inShape: inShape.length, filledIn, filledOut });
  let L = await layersNow(page);
  check(!L.some((l) => /^lasso-bins/.test(l.id)) && L.some((l) => /^ol-fill-/.test(l.id)), `${B}: Draw a region: the dots have handed over to the outlines (no dot layer for the shape)`, L.filter((l) => /^(lasso|ol-)/.test(l.id)).map((l) => l.id).slice(0, 8));
  const sm = await samples(page, inShape, 40);
  let pin = await px(page, sm.map((q) => q.s)), pc = await px(page, sm.map((q) => q.c));
  const beyond = sm.filter((q) => q.beyond).length;
  const pinBeyond = pin.filter((p, i) => sm[i].beyond);
  check(sm.length >= 5 && beyond >= 5 && share(pin, reddish) >= 0.9 && share(pinBeyond, reddish) >= 0.9 && share(pc, reddish) >= 0.9,
    `${B}: pixels: a selected cell's polygon is coloured inside, away from its centre, not only at the centre (${sm.length} cells, ${beyond} points outside every dot)`,
    { cells: sm.length, beyond, inside: share(pin, reddish), insideBeyond: share(pinBeyond, reddish), centre: share(pc, reddish), sample: pin.slice(0, 4) });
  await shotB(page, "select_drawn_region.png");
  // picking the shape lights its cells and dims the rest of the shapes' cells, as with the dots
  await page.evaluate(() => document.querySelector("#lassoShapes .lshape").click());
  await settled(page);
  col = await colOf(page);
  check(inShape.every((r) => col.fill[4 * r + 3] === 235 && col.line[4 * r + 3] === 255), `${B}: the shape picked: its cells' outlines filled at the dots' lit alpha (235), their line in its colour`,
    inShape.slice(0, 3).map((r) => [col.fill[4 * r + 3], col.line[4 * r + 3]]));
  const pinLit = await px(page, sm.map((q) => q.s));
  check(share(pinLit, (p) => reddish(p) && p[0] >= 150) >= 0.9, `${B}: pixels: the picked shape's cells are bright inside`, pinLit.slice(0, 4));
  await shotB(page, "select_drawn_region_picked.png");
  await page.evaluate(() => document.querySelector("#lassoShapes .lshape").click());
  await settled(page);

  // ---- the dots come back when the outlines are off (Cell shapes: Dots)
  await page.evaluate(() => { const b = document.getElementById("tabView"); if (b) b.click(); });
  await setVal(page, "cellShapes", "dots"); await settled(page);
  s = await st(page); L = await layersNow(page);
  const dotsLayer = L.find((l) => /^lasso-bins/.test(l.id));
  pc = await px(page, sm.map((q) => q.c));
  const pBeyond = await px(page, sm.filter((q) => q.beyond).map((q) => q.s));
  check(s.alpha === 0 && !L.some((l) => /^ol-/.test(l.id)) && dotsLayer && dotsLayer.op === 1 && share(pc, reddish) >= 0.9 && pBeyond.length >= 5 && share(pBeyond, reddish) <= 0.1,
    `${B}: Cell shapes: Dots: no outline layer, the shape's dots are back (centres coloured, points inside the polygons but outside the dots not)`,
    { alpha: s.alpha, dots: dotsLayer || null, centre: share(pc, reddish), beyondRed: pBeyond.length ? share(pBeyond, reddish) : null });
  await shotB(page, "select_dots.png");
  await setVal(page, "cellShapes", "outlines"); await settled(page);
  check((await st(page)).alpha === 1, `${B}: Cell shapes: Outlines again`);

  // ---- Save image: the PNG shows the outlines as on screen, and its legend names them
  const leg = await page.evaluate(() => window.__spatialViewer.snapLegend());
  const cs = leg.find((x) => x.title === "Cell shapes");
  check(cs && cs.items.length >= 1 && cs.items[0].outline, `${B}: Save image's legend: a Cell shapes entry while the cells are outlines`, leg.map((x) => x.title));
  await page.evaluate(() => {
    window.__blobs = []; const orig = URL.createObjectURL.bind(URL);
    URL.createObjectURL = (b) => { window.__blobs.push(b); return orig(b); };
    const sc = document.getElementById("snapScale"); sc.dispatchEvent(new Event("focus")); sc.value = "1";
    for (const [id, on] of [["snapLegend", true], ["snapName", true], ["snapBar", false]]) document.getElementById(id).checked = on;
    document.getElementById("snapBg").value = "black";
  });
  await page.evaluate(() => document.getElementById("snapSave").click());
  await page.waitForFunction(() => window.__blobs.some((b) => b.type === "image/png") && !document.getElementById("snapSave").disabled, null, { timeout: 60000 });
  const png = await page.evaluate(async (pts) => {
    const b = window.__blobs.filter((x) => x.type === "image/png").pop(), url = URL.createObjectURL(b);
    const img = new Image(); await new Promise((res, rej) => { img.onload = res; img.onerror = rej; img.src = url; });
    const sv = window.__spatialViewer, cv = sv.deck.getCanvas ? sv.deck.getCanvas() : sv.deck.canvas, k = cv.width / cv.clientWidth, top = img.height - cv.height;
    const c = document.createElement("canvas"); c.width = img.width; c.height = img.height;
    const g = c.getContext("2d"); g.drawImage(img, 0, 0);
    return { w: img.width, h: img.height, mapW: cv.width, mapH: cv.height, px: pts.map(([x, y]) => Array.from(g.getImageData(Math.round(x * k), top + Math.round(y * k), 1, 1).data)) };
  }, sm.map((q) => q.s));
  check(png.w > png.mapW && png.h > png.mapH && share(png.px, reddish) >= 0.9,
    `${B}: Save image: the PNG (${png.w} x ${png.h}, map ${png.mapW} x ${png.mapH}, with legend and name) shows the selected cells' polygons filled`, { share: share(png.px, reddish), sample: png.px.slice(0, 4) });

  // ---- a link to the view: the camera, the shape and the outlines come back; and the link's Dots choice wins
  const link = await page.evaluate(async () => (await window.SodaViewLink.make()).url);
  const vs0 = await page.evaluate(() => { const v = window.__spatialViewer.state.viewState; return [v.target[0], v.target[1], v.zoom]; });
  let p2 = await ctx.newPage(); watch(p2);
  await p2.goto(link); await ready(p2); await settled(p2);
  const back = await p2.evaluate(() => { const sv = window.__spatialViewer, S = sv.state, v = S.viewState, L = S.lasso; return { vs: [v.target[0], v.target[1], v.zoom], shapes: L.shapes.length, bins: L.bin ? L.bin.n : 0, imgOn: S.imgOn, notes: window.SodaViewLink.notes() }; });
  s = await st(p2);
  const pLink = await px(p2, sm.map((q) => q.s));
  check(Math.abs(back.vs[0] - vs0[0]) < 1 && Math.abs(back.vs[1] - vs0[1]) < 1 && Math.abs(back.vs[2] - vs0[2]) < 0.01 && back.shapes === 1 && back.bins === inShape.length && s.alpha === 1 && share(pLink, reddish) >= 0.9,
    `${B}: a link to the view: same camera, the shape and its cells' filled outlines come back`, { back, alpha: s.alpha, share: share(pLink, reddish) });
  await shotB(p2, "select_link_back.png");
  await p2.close();
  await page.evaluate(() => window.__spatialViewer.outlines.setMode("dots"));
  const linkDots = await page.evaluate(async () => (await window.SodaViewLink.make()).url);
  await page.evaluate(() => window.__spatialViewer.outlines.setMode("outlines"));
  p2 = await ctx.newPage(); watch(p2);
  await p2.goto(linkDots); await ready(p2); await settled(p2);
  s = await st(p2);
  check(s.mode === "dots" && s.alpha === 0 && (await p2.evaluate(() => document.getElementById("cellShapes").value)) === "dots",
    `${B}: a link made with Cell shapes on Dots opens on Dots (whatever this browser had chosen)`, { mode: s.mode, alpha: s.alpha });
  await p2.close();

  check(errors.length === 0, `${B}: no console error, page error or failed request in the selection checks`, errors.slice(0, 5));
  await ctx.close();
}

(async () => {
  const data = await serve(EXPORT, PORTS[0], true), pageSrv = await serve(PAGE, PORTS[1], false);
  const DATA = `http://127.0.0.1:${data.address().port}`, URL0 = `http://127.0.0.1:${pageSrv.address().port}/index.html`;
  const browser = await pw.chromium.launch({ args: GPU ? GPU_ARGS : SOFT_ARGS });
  // serviceWorkers "block": the page registers a service worker (sw.js, offline use) that would answer some requests
  // itself, out of page.route's sight; tests/offline_test.js tests it
  const ctx = await browser.newContext({ viewport: { width: 1366, height: 768 }, serviceWorkers: "block" });
  const errors = [];
  const open = async (ds) => {
    const page = await ctx.newPage();
    // vendor/: deck.gl and umap-js are asked for there first and come from the CDN when there is no vendor/ (libs.js)
    page.on("console", (m) => { if (m.type() === "error" && !/\/vendor\//.test((m.location() || {}).url || "")) errors.push(`${ds}: console: ${m.text()}`); });
    page.on("pageerror", (e) => errors.push(`${ds}: pageerror: ${e.message}`));
    page.on("response", (r) => { if (r.status() >= 400 && !/favicon|\/vendor\//.test(r.url())) errors.push(`${ds}: http ${r.status()}: ${r.url()}`); });
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
    // Cell sets: two clusters merged into one set: the second one's cells take the set's colour, as their dots do
    await page.evaluate(() => { const b = document.getElementById("csShowCl"); if (b && !b.closest("[hidden]")) b.click(); });
    const csReady = await page.waitForFunction(() => document.querySelectorAll("#csPick input[data-c]").length >= 2, null, { timeout: 30000 }).then(() => true, () => false);
    if (csReady) {
      const before = await page.evaluate(() => Array.from(window.__spatialViewer.state.umap.cl));
      const codes = await page.evaluate(() => [...document.querySelectorAll("#csPick input[data-c]")].slice(0, 2).map((i) => Number(i.dataset.c)));
      await page.evaluate((c) => window.SodaSets._merge(c, "Outline test set"), codes);
      await settled(page);
      c = await page.evaluate(([before, codes]) => {
        const sv = window.__spatialViewer, col = sv.outlines._debug().col, U = sv.state.umap, f = col.fill;
        const rgb = (r) => `${f[4 * r]},${f[4 * r + 1]},${f[4 * r + 2]},${f[4 * r + 3]}`;
        const a = [], b = [];
        for (let r = 0; r < before.length; r++) { if (before[r] === codes[0]) a.push(rgb(r)); else if (before[r] === codes[1]) b.push(rgb(r)); }
        const ua = [...new Set(a)], ub = [...new Set(b)];
        return { set: U.cl[before.indexOf(codes[1])], head: codes[0], first: ua, second: ub, n: [a.length, b.length] };
      }, [before, codes]);
      check(c.first.length === 1 && c.second.length === 1 && c.first[0] === c.second[0] && c.set === c.head && c.n[1] > 0,
        "Cell sets: clusters merged into a set fill their cells' outlines in the set's one colour", c);
      await shot(page, "outlines_xenium_cell_set.png");
    } else check(false, "Cell sets: the section lists the clusters to merge");
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
  await browser.close();

  // ---- what draws a subset of the cells, with outlines: in each browser asked for
  for (const name of BROWSERS) {
    if (!pw[name]) { console.log(`${name}: not in this Playwright, skipped`); continue; }
    let b = null;
    try { b = await pw[name].launch(name === "chromium" ? { args: GPU ? GPU_ARGS : SOFT_ARGS } : {}); }
    catch (e) { console.log(`SKIPPED ${name}: not installed (${String(e.message || e).split("\n")[0]})`); continue; }
    try { await selectionChecks(name, b, DATA, URL0); }
    catch (e) { check(false, `${name}: the selection checks ran to the end`, String(e && e.stack || e).slice(0, 400)); }
    await b.close();
  }

  data.close(); pageSrv.close();
  const bad = results.filter((r) => !r.ok);
  console.log(bad.length ? `\nFAILED: ${bad.length} of ${results.length} checks` : `\nPASSED: ${results.length} checks`);
  process.exit(bad.length ? 1 : 0);
})().catch((e) => { console.error(e); process.exit(2); });
