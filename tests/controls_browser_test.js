#!/usr/bin/env node
// One control for each job, in the browser, on the public merscope_2 (cells; tests/test_data.js), in each browser of --browsers (chromium,
// webkit, firefox): Show on slide and Clear selection once each (the chip above the picture, which also names a shape
// picked in Draw a region, with or without a section UMAP), no "Image off" button
// (Image opacity at 0 is off: the image is not drawn), one run button each for the volcano and the plots that shows the
// last result again while nothing has changed (and runs anew once something has), and old links and Save JSON files
// made while the image was off opening with it off.
//
//   node tests/controls_browser_test.js                   (the public exports: SODA_TEST_DATA or --data-root, see tests/test_data.js)
//   node tests/controls_browser_test.js --data-root <url or folder> --page <dir> --shots <dir> --gpu --browsers chromium,webkit,firefox
//
// Playwright is not a dependency of this repo. Install it anywhere (npm i playwright; npx playwright install chromium)
// and point NODE_PATH or --playwright at its node_modules. Without it the test says so and exits 0; a browser that is
// not installed is skipped with a line saying so. Chromium draws with software rendering unless --gpu is given.
"use strict";
const http = require("http");
const fs = require("fs");
const path = require("path");
const { sendFile, preflight } = require("./static_files");

const arg = (k, d) => { const i = process.argv.indexOf(k); return i > 0 ? process.argv[i + 1] : d; };
const REPO = path.resolve(__dirname, "..");
const D = require("./test_data");                              // the public exports: SODA_TEST_DATA / --data-root
const PAGE = path.resolve(arg("--page", path.join(REPO, "web", "SpatialViewer")));
const SHOTS = arg("--shots", "");
const GPU = process.argv.includes("--gpu");
const BROWSERS = arg("--browsers", "chromium").split(",").map((s) => s.trim()).filter(Boolean);
const GPU_ARGS = ["--use-angle=d3d11", "--ignore-gpu-blocklist", "--enable-gpu"];
const SOFT_ARGS = ["--use-angle=swiftshader", "--enable-unsafe-swiftshader", "--ignore-gpu-blocklist"];
// Firefox's console lines that are not the page's faults (as in tests/browsers_test.js): the vendor/ tag's SHA-384
// checked against a 404 when there is no vendor/, and deck.gl's WebAssembly probe refused by the page's CSP
const FIREFOX_QUIET = [/integrity attribute match the content of the subresource at .[^ ]*\/vendor\//,
  /Content-Security-Policy: .*blocked WebAssembly .*Missing 'wasm-unsafe-eval'/];
const DS = "merscope_2";                                       // a cells layer, small enough for three browsers

let pw = null;
for (const base of [arg("--playwright", ""), ...(process.env.NODE_PATH || "").split(path.delimiter), path.join(REPO, "node_modules")].filter(Boolean)) {
  try { pw = require(path.join(base, "playwright")); break; } catch (e) { /* next */ }
}
if (!pw) { try { pw = require("playwright"); } catch (e) { /* none */ } }
if (!pw) { console.log("SKIPPED: playwright is not installed (see the top of this file)"); process.exit(0); }

const TYPES = { ".html": "text/html", ".js": "text/javascript", ".css": "text/css", ".json": "application/json", ".png": "image/png", ".jpg": "image/jpeg", ".bin": "application/octet-stream" };
// root null: the data (tests/test_data.js: a folder, or passed on to a URL)
function serve(root) {
  return new Promise((resolve) => {
    const srv = http.createServer((req, res) => {
      const rel = decodeURIComponent(req.url.split("?")[0]);
      if (!root) { if (req.method === "OPTIONS") { preflight(req, res); return; } D.serve(req, res, rel.slice(1), TYPES[path.extname(rel)]); return; }
      // config.js is each deployment's own file: the test brings its own
      if (rel === "/config.js") { res.writeHead(200, { "Content-Type": "text/javascript" }); res.end(`window.SpatialViewerConfig = { tokenEndpoint: "" };`); return; }
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

async function run(B, browser, DATA, URL0) {
  console.log(`${B}: ${DS}`);
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
      if (s.busy || (s.alpha > 0 && !(sv.state.umap && sv.state.umap.pos && O.col && O.col.pos))) { window.__idle = 0; return false; }
      window.__idle = (window.__idle || 0) + 1; return window.__idle >= 4;
    }, null, { timeout: 60000, polling: 100 });
    await page.evaluate(() => new Promise((r) => requestAnimationFrame(() => requestAnimationFrame(r))));
    await page.waitForTimeout(300);
  };
  const shot = async (page, name) => { if (SHOTS) { fs.mkdirSync(SHOTS, { recursive: true }); await page.screenshot({ path: path.join(SHOTS, `${B}_${name}`) }); } };
  const setVal = (page, id, v, ev = "change") => page.evaluate(([id, v, ev]) => { const el = document.getElementById(id); el.value = v; el.dispatchEvent(new Event(ev, { bubbles: true })); }, [id, v, ev]);
  const opacity = async (page, v) => { await setVal(page, "imgOpacity", String(v), "input"); await setVal(page, "imgOpacity", String(v), "change"); };
  const zoomTo = async (page, cellPx) => {
    // centred on the middle row's cell (a real section is not centred in its image: its middle can be bare slide)
    const zv = await page.evaluate((cp) => { const S = window.__spatialViewer.state, e = S.m.expr.find((x) => x.key === S.exprKey), p = S.lasso.pos, r = Math.floor(p.length / 4); return [p[2 * r], p[2 * r + 1], Math.log2(cp / e.spotDiameterPx)]; }, cellPx);
    await page.evaluate(([x, y, z]) => { const sv = window.__spatialViewer; sv.deck.props.onViewStateChange({ viewState: { ...sv.state.viewState, target: [x, y, 0], zoom: z }, viewId: "ortho" }); }, zv);
    await settled(page);
  };

  let page = await ctx.newPage(); watch(page);
  await page.goto(`${URL0}?static=${encodeURIComponent(DATA)}&ds=${DS}`);
  await ready(page);

  // ---- one control for each job
  const ctl = await page.evaluate(() => {
    const btn = [...document.querySelectorAll("button")], n = (t) => btn.filter((b) => b.textContent.trim() === t).length;
    return { gone: ["btnImage", "selShow", "selClear", "vcShow", "plShow", "lpSel", "lpSelShow", "lpSelClear"].filter((id) => document.getElementById(id)),
             showOnSlide: n("Show on slide"), clear: n("Clear selection"), imageBtn: n("Image off") + n("Image on"), show: n("Show"),
             chip: [document.getElementById("selChipShow"), document.getElementById("selChipClear")].every((b) => b && b.closest("#selChip")),
             opacity: !!document.getElementById("imgOpacity"), runs: !!document.getElementById("vcRun") && !!document.getElementById("plRun") };
  });
  check(!ctl.gone.length && ctl.showOnSlide === 1 && ctl.clear === 1 && ctl.chip && ctl.imageBtn === 0 && ctl.show === 0 && ctl.opacity && ctl.runs,
    `${B}: one of each: Show on slide and Clear selection once (the chip above the picture), no Image off/on button, no separate Show button`, ctl);
  // Image opacity: 0 is off (not drawn), above 0 on again
  const imgState = () => page.evaluate(() => { const sv = window.__spatialViewer, S = sv.state; return { on: S.imgOn, op: S.imgOpacity, slider: +document.getElementById("imgOpacity").value,
    layers: sv.deck.layerManager.getLayers().filter((l) => /^ch-/.test(l.id)).map((l) => +l.props.opacity.toFixed(2)) }; });
  // deck takes new layers on its next frame: wait for it (a loaded machine can take a while)
  const layersSettle = (want) => page.waitForFunction((w) => { const n = window.__spatialViewer.deck.layerManager.getLayers().filter((l) => /^ch-/.test(l.id)).length; return w ? n > 0 : n === 0; }, want, { timeout: 15000, polling: 100 }).catch(() => {});
  await opacity(page, 0); await layersSettle(false);
  let im = await imgState();
  check(im.on === false && im.op === 0 && im.layers.length === 0, `${B}: Image opacity at 0: the image is off (no stain layer drawn)`, im);
  await opacity(page, 0.5); await layersSettle(true); await page.waitForTimeout(300);
  im = await imgState();
  check(im.on === true && im.op === 0.5 && im.layers.length > 0 && im.layers.includes(0.5) && im.layers.every((o) => o > 0), `${B}: Image opacity back to 0.5: the image is on again at 0.5`, im);
  await opacity(page, 0);

  // a shape drawn on the slide: a group for the volcano and the plots
  await zoomTo(page, 12);
  const vbox = await page.evaluate(() => { const r = document.getElementById("viewport").getBoundingClientRect(); return { x: r.left, y: r.top, w: r.width, h: r.height }; });
  await page.evaluate(() => { if (window.__spatialViewer.state.lasso.tool !== "circle") document.getElementById("lassoCircle").click(); });
  const cx = vbox.x + vbox.w / 2, cy = vbox.y + vbox.h / 2, rad = Math.min(vbox.w, vbox.h) * 0.22;
  await page.mouse.move(cx, cy); await page.mouse.down();
  for (let i = 1; i <= 8; i++) await page.mouse.move(cx + (rad * i) / 8, cy);
  await page.mouse.up();
  await page.evaluate(() => { if (window.__spatialViewer.state.lasso.tool) document.getElementById("lassoCircle").click(); });
  await page.mouse.move(vbox.x + 3, vbox.y + 3);
  await page.waitForFunction(() => { const L = window.__spatialViewer.state.lasso; return L.shapes.length === 1 && L.bin && L.bin.n > 0; }, null, { timeout: 30000 });

  // ---- a drawn shape picked: named in the chip above the picture, with its one Show on slide and Clear selection
  // (Draw a region's panel had its own copy of the two; now the chip is the one place, also without a section UMAP)
  const chip = () => page.evaluate(() => {
    const c = document.getElementById("selChip"), up = (el) => !!el && !el.closest("[hidden]") && el.offsetParent !== null;
    const vis = [...document.querySelectorAll("button")].filter(up), n = (t) => vis.filter((b) => b.textContent.trim() === t).length;
    return { shown: up(c), txt: document.getElementById("selChipTxt").textContent, show: n("Show on slide"), clear: n("Clear selection"),
             seg: up(document.querySelector("#viewSwitch .seg")), pick: window.__spatialViewer.state.lasso.selShape };
  });
  for (const noUmap of [false, true]) {
    const what = noUmap ? " (a layer without a section UMAP)" : "";
    if (noUmap) await page.evaluate(() => { const S = window.__spatialViewer.state, e = S.m.expr.find((x) => x.key === S.exprKey); window.__umapWas = e.umap; e.umap = null; });
    await page.evaluate(() => document.querySelector("#lassoShapes .lshape").click());
    await page.waitForTimeout(400);
    let ch = await chip();
    check(ch.shown && ch.pick === 0 && /Shape 1|Drawn region/.test(ch.txt) && ch.show === 1 && ch.clear === 1 && (!noUmap || !ch.seg),
      `${B}: a drawn shape picked: the chip above the picture names it, with one Show on slide and one Clear selection${what}`, ch);
    // away from the shape first, so framing it has somewhere to move to
    await page.evaluate(() => { const sv = window.__spatialViewer, v = sv.state.viewState, i = sv.state.m.image; sv.deck.props.onViewStateChange({ viewState: { ...v, target: [i.width * 0.1, i.height * 0.1, 0] }, viewId: "ortho" }); });
    await page.waitForTimeout(300);
    const vs0 = await page.evaluate(() => JSON.stringify(window.__spatialViewer.state.viewState.target));
    await page.evaluate(() => document.getElementById("selChipShow").click());
    await page.waitForFunction((v) => JSON.stringify(window.__spatialViewer.state.viewState.target) !== v, vs0, { timeout: 10000 }).catch(() => {});
    const moved = await page.evaluate((v) => JSON.stringify(window.__spatialViewer.state.viewState.target) !== v, vs0);
    await page.evaluate(() => document.getElementById("selChipClear").click());
    await page.waitForTimeout(400);
    ch = await chip();
    check(moved && ch.pick === null && !ch.shown, `${B}: the chip's Show on slide frames the shape, its Clear selection unpicks it and the chip goes${what}`, { moved, ch });
    if (noUmap) await page.evaluate(() => { const S = window.__spatialViewer.state; S.m.expr.find((x) => x.key === S.exprKey).umap = window.__umapWas; });
  }

  // ---- the volcano's and the plots' one button
  const one = async (run, panel, close, label, wait) => {
    await page.evaluate((id) => document.getElementById(id).click(), run);
    const ok = await page.waitForFunction(wait, null, { timeout: 180000 }).then(() => true, () => false);
    await page.evaluate((id) => document.getElementById(id).click(), close);
    const t1 = await page.evaluate((id) => document.getElementById(id).textContent, run);
    await page.evaluate(([w]) => { window.__res = window[w].result; }, [label === "volcano" ? "SodaVolcano" : "SodaPlots"]);
    await page.evaluate((id) => document.getElementById(id).click(), run);
    await page.waitForTimeout(400);
    const again = await page.evaluate(([w, p]) => ({ same: window[w].result === window.__res, shown: !document.getElementById(p).hidden }), [label === "volcano" ? "SodaVolcano" : "SodaPlots", panel]);
    const t2 = await page.evaluate((id) => document.getElementById(id).textContent, run);
    return { ok, t1, again, t2 };
  };
  await page.evaluate(() => { const b = document.getElementById("tabAnalyze"); if (b) b.click(); window.SodaVolcano.refresh(); });
  const vg = await page.evaluate(() => {
    const a = document.getElementById("vcA"); a.dispatchEvent(new Event("pointerdown"));
    if (![...a.options].some((o) => o.value === "s:0")) return null;
    a.value = "s:0"; a.dispatchEvent(new Event("change"));
    const b = document.getElementById("vcB"); b.value = "rest"; b.dispatchEvent(new Event("change"));
    return [a.value, b.value];
  });
  if (vg) {
    const v = await one("vcRun", "volcanoPanel", "vcClose", "volcano", () => window.SodaVolcano.result && !document.getElementById("vcRun").disabled);
    check(v.ok && v.t1 === "Show volcano" && v.again.same && v.again.shown && v.t2 === "Volcano",
      `${B}: the volcano's one button: closed, it reads Show volcano and brings the same volcano back without running again`, v);
    await page.evaluate(() => document.getElementById("vcClose").click());
    await setVal(page, "vcTest", (await page.evaluate(() => document.getElementById("vcTest").value)) === "welch" ? "wilcoxon" : "welch");
    check((await page.evaluate(() => document.getElementById("vcRun").textContent)) === "Volcano", `${B}: the volcano's settings changed: the button runs a new one (reads Volcano)`);
  } else check(false, `${B}: the volcano offers the drawn shape`);
  await page.evaluate(() => window.SodaPlots.refresh());
  const pg = await page.evaluate(() => { const s = document.getElementById("plGroup"); s.dispatchEvent(new Event("pointerdown")); if (![...s.options].some((o) => o.value === "shapes")) return null; s.value = "shapes"; s.dispatchEvent(new Event("change")); document.getElementById("plGenes").value = ""; return s.value; });
  const g2 = await page.evaluate(() => [...window.__spatialViewer.state.geneIndex].sort((a, b) => b.t - a.t).slice(0, 3).map((g) => g.n).join(", "));
  if (pg) {
    await setVal(page, "plGenes", g2, "input");
    const p = await one("plRun", "plotsPanel", "plClose", "plots", () => window.SodaPlots.result && document.querySelector("#plSvgBox svg") && !document.getElementById("plRun").disabled);
    check(p.ok && p.t1 === "Show plots" && p.again.same && p.again.shown && p.t2 === "Draw plots",
      `${B}: the plots' one button: closed, it reads Show plots and brings the same plots back without drawing again`, p);
    await page.evaluate(() => document.getElementById("plClose").click());
    await setVal(page, "plGenes", g2.split(", ").slice(0, 2).join(", "), "input");
    check((await page.evaluate(() => document.getElementById("plRun").textContent)) === "Draw plots", `${B}: the plots' genes changed: the button draws new ones (reads Draw plots)`);
  } else check(false, `${B}: the plots offer the drawn shapes as a grouping`);

  // ---- old links and Save JSON files: the image off stays off
  await opacity(page, 0.6);
  const mk = await page.evaluate(async () => {
    const url = (await window.SodaViewLink.make()).url, d = await window.SodaViewLink.decode(url.split("#v=")[1]);
    const v = (d.state && d.state.view) || d.state || d;
    const find = (o) => { if (!o || typeof o !== "object") return null; if (o.img && typeof o.img === "object") return o; for (const x of Object.values(o)) { const f = find(x); if (f) return f; } return null; };
    const host = find(d.state), out = { now: JSON.stringify(host && host.img) };
    const enc = async (img) => { host.img = img; return url.split("#v=")[0] + "#v=" + (await window.SodaViewLink.encode(d.state)); };
    out.off = await enc({ on: false, op: 0.6 });     // made while "Image off" was pressed: opacity kept, image off
    out.on = await enc({ on: true, op: 0.6 });
    out.bare = await enc({ on: false });
    return out;
  });
  check(mk.now === JSON.stringify({ on: true, op: 0.6 }), `${B}: a link holds the image as before ({on, op})`, mk.now);
  for (const [u, want, what] of [[mk.off, { on: false, op: 0 }, "an old link with the image off (Image off pressed, opacity 0.6) opens with it off (Image opacity 0)"],
                                 [mk.on, { on: true, op: 0.6 }, "a link with the image on at 0.6 opens on at 0.6"],
                                 [mk.bare, { on: false, op: 0 }, "a link with only on: false opens with the image off"]]) {
    const p2 = await ctx.newPage(); watch(p2);
    await p2.goto(u); await ready(p2);
    await p2.waitForFunction((on) => (window.__spatialViewer.deck.layerManager.getLayers().filter((l) => /^ch-/.test(l.id)).length > 0) === on, want.on, { timeout: 15000, polling: 100 }).catch(() => {});
    const got = await p2.evaluate(() => { const S = window.__spatialViewer.state; return { on: S.imgOn, op: S.imgOpacity, slider: +document.getElementById("imgOpacity").value, layers: window.__spatialViewer.deck.layerManager.getLayers().filter((l) => /^ch-/.test(l.id)).length }; });
    check(got.on === want.on && got.op === want.op && got.slider === want.op && (got.layers > 0) === want.on, `${B}: ${what}`, got);
    await p2.close();
  }
  // Save JSON (Draw a region's file holds the same view): one saved while Image off was pressed loads with the image off
  const sj = await page.evaluate(async () => {
    const sv = window.__spatialViewer;
    const view = window.SodaViewLink.capture(); view.img = { on: false, op: 0.6 };
    await sv.lassoLoad({ dataset: sv.state.dsId, layer: sv.state.exprKey, shapes: [], view });
    await new Promise((r) => setTimeout(r, 800));
    return { on: sv.state.imgOn, op: sv.state.imgOpacity, slider: +document.getElementById("imgOpacity").value };
  });
  check(sj.on === false && sj.op === 0 && sj.slider === 0, `${B}: a Save JSON file saved with the image off loads with it off`, sj);

  check(errors.length === 0, `${B}: no console error, page error or failed request`, errors.slice(0, 5));
  await ctx.close();
}

(async () => {
  await D.need([DS]);
  console.log(`data: ${D.describe()}`);
  const data = await serve(null), pageSrv = await serve(PAGE);
  const DATA = `http://127.0.0.1:${data.address().port}`, URL0 = `http://127.0.0.1:${pageSrv.address().port}/index.html`;
  for (const name of BROWSERS) {
    if (!pw[name]) { console.log(`SKIPPED: ${name} is not in this Playwright`); continue; }
    let b = null;
    try { b = await pw[name].launch(name === "chromium" ? { args: GPU ? GPU_ARGS : SOFT_ARGS } : {}); }
    catch (e) { console.log(`SKIPPED: ${name} is not installed (${String(e.message || e).split("\n")[0]})`); continue; }
    try { await run(name, b, DATA, URL0); }
    catch (e) { check(false, `${name}: the checks ran to the end`, String(e && e.stack || e).slice(0, 600)); }
    await b.close();
  }
  data.close(); pageSrv.close();
  const bad = results.filter((r) => !r.ok);
  console.log(bad.length ? `\nFAILED: ${bad.length} of ${results.length} checks` : `\nPASSED: ${results.length} checks`);
  process.exit(bad.length ? 1 : 0);
})().catch((e) => { console.error(e); process.exit(2); });
