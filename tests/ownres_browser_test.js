#!/usr/bin/env node
// Your own results (web/SpatialViewer/ownres.js) in the browser: a reader loads a table of results they already have
// (UMAP, clusterings, cell types, numbers) from a file, with no re-export. On the demo's Visium HD export with barcode
// lists (export_user/demo_visiumhd, 16 um layer; tests/smoke_test.py makes it with `prepare.py ids`), in each browser of
// --browsers (chromium, webkit):
//   a Seurat-style CSV (write.csv: an empty first header cell, a sample prefix on the barcodes, a tenth of the bins
//   dropped, seurat_clusters, cell_type with a colour column, a number, UMAP_1 / UMAP_2): roles guessed, the match rate,
//   colours on the tissue and on the file's UMAP checked by pixel, Colour by a value with the file's number;
//   a scanpy-style TSV, gzipped (obs.to_csv: "-1-0" endings, leiden as text, X_umap1 / X_umap2), with an uns-style
//   colour list pasted in the role table; a partial match (Use anyway); a file whose roles are guessed wrong and put
//   right in the table; DE (the volcano) on a loaded clustering; cell sets on it; Save JSON with the table in it and
//   Load JSON on a fresh page; a link to the view (a fingerprint: it asks for the file, then the colouring comes back);
//   drag and drop; junk files (each a clear message, no page error); a layer without a barcode list (how to add it, and
//   a "# soda row order" file) and a layer without SODA's UMAP; and the section at 375 px phone width.
//
//   node tests/ownres_browser_test.js                         (needs demo/: run tests/smoke_test.py first)
//   node tests/ownres_browser_test.js --root <demo dir> --page <dir> --shots <dir> --browsers chromium,webkit --gpu
//
// Playwright is not a dependency of this repo (see tests/outlines_browser_test.js): without it the test says so and
// exits 0. The fixtures are made here from the export's own ids.txt.gz and positions.bin; nothing is downloaded.
"use strict";
const http = require("http");
const fs = require("fs");
const path = require("path");
const zlib = require("zlib");
const { sendFile, preflight } = require("./static_files");

const arg = (k, d) => { const i = process.argv.indexOf(k); return i > 0 ? process.argv[i + 1] : d; };
const REPO = path.resolve(__dirname, "..");
const ROOT = path.resolve(arg("--root", path.join(REPO, "demo")));
const PAGE = path.resolve(arg("--page", path.join(REPO, "web", "SpatialViewer")));
const SHOTS = arg("--shots", "");
const GPU = process.argv.includes("--gpu");
const BROWSERS = arg("--browsers", "chromium").split(",").map((s) => s.trim()).filter(Boolean);
const GPU_ARGS = ["--use-angle=d3d11", "--ignore-gpu-blocklist", "--enable-gpu"];
const SOFT_ARGS = ["--use-angle=swiftshader", "--enable-unsafe-swiftshader", "--ignore-gpu-blocklist"];
const DS = "demo_visiumhd", KEY = "016um";

let pw = null;
for (const base of [arg("--playwright", ""), ...(process.env.NODE_PATH || "").split(path.delimiter), path.join(REPO, "node_modules")].filter(Boolean)) {
  try { pw = require(path.join(base, "playwright")); break; } catch (e) { /* next */ }
}
if (!pw) { try { pw = require("playwright"); } catch (e) { /* none */ } }
if (!pw) { console.log("SKIPPED: playwright is not installed (see the top of tests/outlines_browser_test.js)"); process.exit(0); }
const EXP = path.join(ROOT, "export_user", DS, "expr", KEY);
if (!fs.existsSync(path.join(EXP, "ids.txt.gz"))) { console.log(`SKIPPED: no ${DS} export with barcode lists under ${ROOT}/export_user (run tests/smoke_test.py first, or pass --root)`); process.exit(0); }

// ---- the fixtures, from the export's own ids and positions --------------------------------------------------------
const f32 = (p) => { const b = fs.readFileSync(p); return new Float32Array(b.buffer.slice(b.byteOffset, b.byteOffset + b.length)); };
const IDS = zlib.gunzipSync(fs.readFileSync(path.join(EXP, "ids.txt.gz"))).toString("utf-8").split("\n").filter(Boolean);
const POS = f32(path.join(EXP, "positions.bin"));
const N = IDS.length;
if (POS.length !== 2 * N) { console.log(`FAIL: ${KEY} ids ${N} vs positions ${POS.length / 2}`); process.exit(1); }
// three bands across the tissue, by x: a known group per bin, far enough from the band edges to read by pixel
const xs = Array.from({ length: N }, (_, r) => POS[2 * r]).sort((a, b) => a - b);
const T1 = xs[Math.floor(N / 3)], T2 = xs[Math.floor((2 * N) / 3)], X0 = xs[0], X1 = xs[N - 1];
let Y0 = Infinity, Y1 = -Infinity; for (let r = 0; r < N; r++) { Y0 = Math.min(Y0, POS[2 * r + 1]); Y1 = Math.max(Y1, POS[2 * r + 1]); }
const band = (r) => (POS[2 * r] < T1 ? 0 : POS[2 * r] < T2 ? 1 : 2);
const NAMES = ["Left band", "Middle band", "Right band"], COLORS = ["#ff2020", "#20e020", "#3050ff"];
const jit = (r, k) => (((r * 7919 + k * 104729) % 1000) / 1000 - 0.5) * 2;
const umapOf = (r) => [band(r) * 10 + jit(r, 1), ((POS[2 * r + 1] - Y0) / Math.max(1, Y1 - Y0)) * 8 + jit(r, 2) * 0.3];
const scoreOf = (r) => ((POS[2 * r] - X0) / Math.max(1, X1 - X0)) * 100;
const q = (s) => `"${s}"`;
const FIX = {};
// Seurat: write.csv(cbind(obj@meta.data, Embeddings(obj, "umap")), "seurat.csv") of a merged object (slice1_ prefix),
// a tenth of the bins dropped in QC
{
  const L = [`"","orig.ident","nCount_Spatial","seurat_clusters","cell_type","cell_type_color","UMAP_1","UMAP_2"`];
  for (let r = 0; r < N; r++) {
    if (r % 10 === 3) continue;
    const g = band(r), [ux, uy] = umapOf(r);
    L.push([q(`slice1_${IDS[r]}`), q("SeuratProject"), (scoreOf(r)).toFixed(3), q(String(g)), q(NAMES[g]), q(COLORS[g]), ux.toFixed(5), uy.toFixed(5)].join(","));
  }
  FIX.seurat = { name: "seurat_results.csv", mimeType: "text/csv", buffer: Buffer.from(L.join("\n") + "\n", "utf-8"), rows: L.length - 1 };
}
// scanpy: obs.to_csv(sep="\t") with X_umap as two columns, gzipped; concatenated objects end barcodes in "-1-0";
// leiden as text, numbered differently from the bands
{
  const L = ["\tleiden\ttotal_counts\tX_umap1\tX_umap2"];
  for (let r = 0; r < N; r++) { const g = band(r), [ux, uy] = umapOf(r); L.push([`${IDS[r]}-0`, String((g + 1) % 3), (scoreOf(r) * 3).toFixed(2), ux.toFixed(4), uy.toFixed(4)].join("\t")); }
  FIX.scanpy = { name: "scanpy_obs.tsv.gz", mimeType: "application/gzip", buffer: zlib.gzipSync(Buffer.from(L.join("\n") + "\n", "utf-8")), rows: N };
}
// partial: 65 % of the rows carry ids of no bin
{
  const L = ["barcode,cluster"];
  for (let r = 0; r < N; r++) L.push(`${r % 20 < 13 ? "zz_" + IDS[r].replace(/^s_/, "q_") : IDS[r]},${band(r)}`);
  FIX.partial = { name: "partial.csv", mimeType: "text/csv", buffer: Buffer.from(L.join("\n") + "\n", "utf-8"), rows: N, matched: L.length - 1 - Array.from({ length: N }, (_, r) => r % 20 < 13).filter(Boolean).length };
}
// guessed wrong: the UMAP under names nobody guesses, a score of whole numbers 0-4 (looks like groups)
{
  const L = ["cellname,dim_a,dim_b,level"];
  for (let r = 0; r < N; r++) { const [ux, uy] = umapOf(r); L.push(`${IDS[r]},${ux.toFixed(4)},${uy.toFixed(4)},${Math.min(4, Math.floor(scoreOf(r) / 20))}`); }
  FIX.wrong = { name: "guess_me.csv", mimeType: "text/csv", buffer: Buffer.from(L.join("\n") + "\n", "utf-8") };
}
// junk: each must give a clear message and no page error
const JUNK = [
  { name: "picture.png", mimeType: "image/png", buffer: Buffer.concat([Buffer.from([0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a]), Buffer.alloc(200)]), want: /not a text table/ },
  { name: "empty.csv", mimeType: "text/csv", buffer: Buffer.from(""), want: /empty/ },
  { name: "one_column.csv", mimeType: "text/csv", buffer: Buffer.from("barcode\nA\nB\n"), want: /Only one column/ },
  { name: "analysis.h5ad", mimeType: "application/octet-stream", buffer: Buffer.concat([Buffer.from([0x89, 0x48, 0x44, 0x46, 0x0d, 0x0a, 0x1a, 0x0a]), Buffer.alloc(100)]), want: /HDF5/ },
  { name: "broken.csv.gz", mimeType: "application/gzip", buffer: Buffer.concat([Buffer.from([0x1f, 0x8b, 0x08, 0x00]), Buffer.from("this is not deflate data at all")]), want: /Could not unpack/ },
  { name: "other_layer.csv", mimeType: "text/csv", buffer: Buffer.from("barcode,cluster\n" + IDS.slice(0, 300).map((b, i) => `${b.replace("_016um_", "_002um_")}x,${i % 3}`).join("\n")), want: /No column of this file matches/ },
  { name: "notes.txt", mimeType: "text/plain", buffer: Buffer.from("Some notes about the run.\nNothing tabular here at all\n"), want: /Only one column|different number of fields/ },
  { name: "workbook.xlsx", mimeType: "application/octet-stream", buffer: Buffer.concat([Buffer.from([0x50, 0x4b, 0x03, 0x04]), Buffer.alloc(64)]), want: /zip or Excel/ },
];

// ---- the servers: the demo root as blob storage would answer it; /noids/ = export_user without its id lists -----------
const TYPES = { ".html": "text/html", ".js": "text/javascript", ".css": "text/css", ".json": "application/json", ".png": "image/png", ".jpg": "image/jpeg", ".bin": "application/octet-stream", ".gz": "application/gzip" };
function serve(root, data) {
  return new Promise((resolve) => {
    const srv = http.createServer((req, res) => {
      let rel = decodeURIComponent(req.url.split("?")[0]);
      if (req.method === "OPTIONS") { preflight(req, res); return; }
      if (!data && rel === "/config.js") { res.writeHead(200, { "Content-Type": "text/javascript" }); res.end("window.SpatialViewerConfig = { tokenEndpoint: \"\" };"); return; }
      let strip = false;
      if (data && rel.startsWith("/noids/")) { strip = true; rel = "/export_user/" + rel.slice("/noids/".length); }
      const file = path.join(root, rel === "/" ? "index.html" : rel);
      if (!file.startsWith(root) || !fs.existsSync(file) || fs.statSync(file).isDirectory()) { res.writeHead(404, { "Access-Control-Allow-Origin": "*" }); res.end("not found"); return; }
      if (strip && /manifest\.json$/.test(file)) {
        const m = JSON.parse(fs.readFileSync(file, "utf-8"));
        for (const L of m.expr || []) delete L.ids;
        res.writeHead(200, { "Content-Type": "application/json", "Access-Control-Allow-Origin": "*", "Cache-Control": "no-store" }); res.end(JSON.stringify(m)); return;
      }
      sendFile(req, res, file, TYPES[path.extname(file)]);
    });
    srv.listen(0, "127.0.0.1", () => resolve(srv));
  });
}

const results = [];
const check = (cond, what, got) => { results.push({ ok: !!cond, what }); console.log(`  ${cond ? "ok  " : "FAIL"} ${what}${cond || got === undefined ? "" : "  got " + JSON.stringify(got).slice(0, 600)}`); return !!cond; };

// ---- in the page: pixels at bins, on the tissue or in the UMAP ----------------------------------------------------
// rows -> the map's pixel under each row's dot (tissue: positions; umap: the UMAP's places), over black
function pagePixelsAt(arg) {
  const [rows, where] = arg;
  const sv = window.__spatialViewer, S = sv.state, U = S.umap, cv = sv.deck.getCanvas ? sv.deck.getCanvas() : sv.deck.canvas, k = cv.width / cv.clientWidth;
  const vp = sv.deck.getViewports().find((v) => v.id === "ortho") || sv.deck.getViewports()[0];
  const src = where === "umap" ? U.uxy : U.pos;
  const c = document.createElement("canvas"); c.width = cv.width; c.height = cv.height;
  const g = c.getContext("2d"); g.fillStyle = "#000"; g.fillRect(0, 0, c.width, c.height); g.drawImage(cv, 0, 0);
  const out = [];
  for (const r of rows) {
    if (!(src[2 * r] === src[2 * r]) || !(src[2 * r + 1] === src[2 * r + 1])) continue;      // not in this UMAP
    const p = vp.project([src[2 * r], src[2 * r + 1]]);
    if (!(p[0] >= 2 && p[1] >= 2 && p[0] < cv.clientWidth - 2 && p[1] < cv.clientHeight - 2)) continue;
    out.push(Array.from(g.getImageData(Math.round((vp.x + p[0]) * k), Math.round((vp.y + p[1]) * k), 1, 1).data));
  }
  return out;
}
const HUE = [(p) => p[0] > p[1] + 40 && p[0] > p[2] + 40, (p) => p[1] > p[0] + 40 && p[1] > p[2] + 40, (p) => p[2] > p[0] + 30 && p[2] > p[1] + 20];
const share = (list, f) => (list.length ? list.filter(f).length / list.length : 0);
// rows well inside each band (and away from the tissue's edge), a few per band
function bandRows(per) {
  const out = [[], [], []], pitch = (X1 - X0) / Math.sqrt(N) * 2.5;
  for (let r = 0; r < N; r++) {
    const x = POS[2 * r], g = band(r), lo = g === 0 ? X0 : g === 1 ? T1 : T2, hi = g === 0 ? T1 : g === 1 ? T2 : X1;
    if (x - lo < pitch || hi - x < pitch || POS[2 * r + 1] - Y0 < pitch || Y1 - POS[2 * r + 1] < pitch) continue;
    if (out[g].length < per && r % 7 === 0 && r % 10 !== 3) out[g].push(r);      // r % 10 = 3: not in the Seurat file
  }
  return out;
}

async function run(B, browser, URL0, DATA) {
  console.log(`${B}:`);
  const ctx = await browser.newContext({ viewport: { width: 1366, height: 768 }, serviceWorkers: "block", acceptDownloads: true });
  const errors = [];
  const watch = (page, tag) => {
    page.on("console", (m) => { if (m.type() === "error" && !/\/vendor\//.test((m.location() || {}).url || "")) errors.push(`${tag}: console: ${m.text()}`); });
    page.on("pageerror", (e) => errors.push(`${tag}: pageerror: ${e.message}`));
    page.on("response", (r) => { if (r.status() >= 400 && !/favicon|\/vendor\//.test(r.url())) errors.push(`${tag}: http ${r.status()}: ${r.url()}`); });
  };
  const shot = async (page, name) => { if (SHOTS) { fs.mkdirSync(SHOTS, { recursive: true }); await page.screenshot({ path: path.join(SHOTS, `${B}_${name}`) }); } };
  const setVal = (page, id, v) => page.evaluate(([id, v]) => { const el = document.getElementById(id); el.value = v; el.dispatchEvent(new Event("change", { bubbles: true })); }, [id, v]);
  const settle = async (page) => {
    await page.waitForFunction(() => { const S = window.__spatialViewer.state; return document.getElementById("status").hidden && !S.ownBusy && !(S.umap && (S.umap.anim || S.umap.busy)); }, null, { timeout: 60000 });
    await page.evaluate(() => new Promise((r) => requestAnimationFrame(() => requestAnimationFrame(r))));
    await page.waitForTimeout(400);
  };
  const open = async (data, ds, key, tag, page) => {
    page = page || await ctx.newPage(); watch(page, tag);
    await page.goto(`${URL0}?static=${encodeURIComponent(data)}&ds=${ds}`);
    await page.waitForFunction((id) => { const sv = window.__spatialViewer; return !!(sv && sv.state && sv.state.dsId === id && sv.state.m && sv.own) && document.getElementById("status").hidden; }, ds, { timeout: 120000 });
    await page.waitForFunction(() => { const S = window.__spatialViewer.state; return S.exprKey && S.umap && !(window.SodaViewLink && window.SodaViewLink.pending()); }, null, { timeout: 120000 });
    if (key && (await page.evaluate(() => window.__spatialViewer.state.exprKey)) !== key) {
      await setVal(page, "exprSelect", key);
      await page.waitForFunction((k) => { const S = window.__spatialViewer.state; return S.exprKey === k && S.umap && document.getElementById("status").hidden; }, key, { timeout: 60000 });
    }
    await settle(page);
    await page.evaluate(() => { if (window.__spatialViewer.state.imgOn) document.getElementById("btnImage").click(); });
    return page;
  };
  const load = async (page, f) => {
    await page.evaluate(() => { document.getElementById("ownMsg").textContent = ""; });
    await page.setInputFiles("#ownFile", { name: f.name, mimeType: f.mimeType, buffer: f.buffer });
    await page.waitForFunction(() => { const t = document.getElementById("ownMsg").textContent; return !window.SodaOwn._state().busy && t.length > 0 && !/^Reading/.test(t); }, null, { timeout: 60000 });
    await settle(page);
  };
  const ownMsg = (page) => page.evaluate(() => ({ text: document.getElementById("ownMsg").textContent, err: document.getElementById("ownMsg").classList.contains("err"), warn: document.getElementById("ownMsg").classList.contains("warn") }));
  const um = (page) => page.evaluate(() => { const U = window.__spatialViewer.state.umap; return { key: U.key, res: U.clRes, colorBy: U.colorBy, view: U.view, nEmb: U.nEmb, embOwn: !!U.embOwn, n: U.cl ? U.cl.length : 0, opts: [...document.getElementById("clusterRes").options].map((o) => o.value) }; });
  const roles = (page) => page.evaluate(() => { const s = window.SodaOwn._state(); return s.D ? Object.fromEntries(s.D.T.header.map((h, j) => [h, s.D.roles[j]])) : null; });
  const px = (page, rows, where) => page.evaluate(pagePixelsAt, [rows, where]);
  const BR = bandRows(30);
  const bandsBy = async (page, where, hues) => {
    const out = [];
    for (let g = 0; g < 3; g++) { const p = await px(page, BR[g], where); out.push({ n: p.length, share: share(p, hues[g]), sample: p.slice(0, 2) }); }
    return out;
  };

  // ---- 1. a Seurat CSV: roles, match rate, the clustering and its colours on the tissue, the UMAP, a number --------
  let page = await open(`${DATA}/export_user`, DS, KEY, "seurat");
  check(await page.evaluate(() => !document.getElementById("panelOwn").hidden), `${B}: the Your own results section is there on a spatial layer`);
  await load(page, FIX.seurat);
  let m = await ownMsg(page), r = await roles(page), u = await um(page);
  const want = { "": "id", "orig.ident": "skip", nCount_Spatial: "num", seurat_clusters: "cat", cell_type: "cat", cell_type_color: "color:4", UMAP_1: "x", UMAP_2: "y" };
  check(JSON.stringify(r) === JSON.stringify(want), `${B}: Seurat CSV: the roles guessed from names and values`, r);
  const matched = FIX.seurat.rows;
  check(m.text.includes(`${matched.toLocaleString()} of ${matched.toLocaleString()} rows matched`) && m.text.includes("sample prefix") && !m.err,
    `${B}: Seurat CSV: the match rate (every row, through the "slice1_" prefix), said`, m.text.slice(0, 300));
  check(!IDS.slice(0, 2000).some((b) => m.text.includes(b)) && !m.text.includes("slice1_"), `${B}: no id or prefix is printed`);
  check(u.res === "own:seurat_clusters" && u.colorBy === "cluster" && u.opts[0] === "own:seurat_clusters" && u.opts.includes("own:cell_type") && u.embOwn && u.nEmb === matched,
    `${B}: the file's groupings are listed first and shown; its UMAP is the one in use`, u);
  await setVal(page, "clusterRes", "own:cell_type"); await settle(page);
  let b = await bandsBy(page, "tissue", HUE);
  check(b.every((x) => x.n >= 10 && x.share >= 0.85), `${B}: cell_type coloured on the tissue in the file's own colours (red / green / blue bands, by pixel)`, b);
  await shot(page, "own_seurat_tissue.png");
  await page.evaluate(() => window.__spatialViewer.umapSetView("umap"));
  await page.waitForFunction(() => { const U = window.__spatialViewer.state.umap; return U.view === "umap" && !U.anim; }, null, { timeout: 30000 }); await settle(page);
  b = await bandsBy(page, "umap", HUE);
  const spread = await page.evaluate((rows) => { const U = window.__spatialViewer.state.umap, xs = rows.map((g) => g.reduce((a, r) => a + U.uxy[2 * r], 0) / g.length); return xs; }, BR);
  check(b.every((x) => x.n >= 10 && x.share >= 0.85) && spread[0] < spread[1] && spread[1] < spread[2],
    `${B}: the file's UMAP (three bands left to right) coloured by cell_type, by pixel`, { b, spread });
  await shot(page, "own_seurat_umap.png");
  // SODA's UMAP instead, and back
  await setVal(page, "ownEmb", "soda");
  await page.waitForFunction(() => { const U = window.__spatialViewer.state.umap; return !U.embOwn && U.uxy && !U.busy; }, null, { timeout: 30000 }); await settle(page);
  u = await um(page);
  check(!u.embOwn && u.view === "umap" && u.res === "own:cell_type", `${B}: UMAP shown: SODA's, with the file's grouping still on it`, u);
  await setVal(page, "ownEmb", "file");
  await page.waitForFunction(() => window.__spatialViewer.state.umap.embOwn, null, { timeout: 30000 }); await settle(page);
  await page.evaluate(() => window.__spatialViewer.umapSetView("tissue"));
  await page.waitForFunction(() => { const U = window.__spatialViewer.state.umap; return U.view === "tissue" && !U.anim; }, null, { timeout: 30000 }); await settle(page);
  // a number of the file in Colour by a value: viridis from the left (dark) to the right (yellow)
  await setVal(page, "colorBy", "gene"); await settle(page);
  const opts = await page.evaluate(() => [...document.getElementById("paintCol").options].map((o) => [o.value, o.textContent]));
  check(opts.some(([v, t]) => v === "own:nCount_Spatial" && /your file/.test(t)), `${B}: the file's number is offered in Colour by a value`, opts);
  await setVal(page, "paintMode", "value"); await page.waitForTimeout(100);
  await setVal(page, "paintCol", "own:nCount_Spatial"); await settle(page);
  const dark = (p) => p[2] > p[1] && p[0] < 120, yellow = (p) => p[1] > 120 && p[1] > p[2] && p[1] > p[0];      // viridis: purple .. green-yellow
  const pl = await px(page, BR[0], "tissue"), pr = await px(page, BR[2], "tissue");
  check(share(pl, dark) >= 0.85 && share(pr, yellow) >= 0.5, `${B}: Colour by a value with the file's number: dark on the left, yellow on the right (by pixel)`, { left: pl.slice(0, 3), right: pr.slice(0, 3), l: share(pl, dark), r: share(pr, yellow) });
  await shot(page, "own_value.png");
  await setVal(page, "paintMode", "off");
  await setVal(page, "colorBy", "cluster"); await settle(page);

  // ---- 2. DE (volcano) on a loaded clustering; cell sets on it -----------------------------------------------------
  await setVal(page, "clusterRes", "own:cell_type"); await settle(page);
  await page.evaluate(() => window.SodaVolcano.refresh());
  const groups = await page.evaluate(() => {
    const a = document.getElementById("vcA"), o = [...a.options].map((x) => [x.value, x.textContent]);
    a.value = "k:0"; a.dispatchEvent(new Event("change"));
    const bb = document.getElementById("vcB"); bb.value = "k:2"; bb.dispatchEvent(new Event("change"));
    return { o, a: a.value, b: bb.value };
  });
  check(groups.o.some(([v, t]) => v === "k:0" && t.startsWith("Left band")) && groups.a === "k:0" && groups.b === "k:2", `${B}: the volcano offers the loaded clustering's groups by name`, groups);
  await page.evaluate(() => { window.SodaVolcano.result = null; document.getElementById("vcRun").click(); });
  const vok = await page.waitForFunction(() => window.SodaVolcano.result || /could not|needs|failed|error/i.test(document.getElementById("vcMsg").textContent), null, { timeout: 180000 }).then(() => true, () => false);
  const V = await page.evaluate(() => { const R = window.SodaVolcano.result, I = R && R.info; return R ? { tested: I.nTested, a: (I.sampledFromA || I.nA) + I.lowA, b: (I.sampledFromB || I.nB) + I.lowB, names: [R.nameA, R.nameB] } : { msg: document.getElementById("vcMsg").textContent }; });
  const sizes = [0, 0, 0]; for (let rr = 0; rr < N; rr++) if (rr % 10 !== 3) sizes[band(rr)]++;
  check(vok && V.tested > 5 && V.a === sizes[0] && V.b === sizes[2], `${B}: DE (volcano) on the loaded clustering: Left band vs Right band, with their sizes`, { V, sizes });
  await page.evaluate(() => { const p = document.getElementById("volcanoPanel"); if (p) p.hidden = true; });
  const before = await page.evaluate(() => Array.from(window.__spatialViewer.state.umap.cl.slice(0, 4000)));
  await page.evaluate(() => window.SodaSets._merge([0, 1], "Left and middle"));
  await settle(page);
  const cs = await page.evaluate((bf) => { const U = window.__spatialViewer.state.umap; let ok = 0, n = 0; for (let r = 0; r < bf.length; r++) if (bf[r] === 1) { n++; if (U.cl[r] === 0) ok++; } return { ok, n, sets: U.sets["own:cell_type"] }; }, before);
  check(cs.n > 0 && cs.ok === cs.n && cs.sets && cs.sets.length === 1, `${B}: cell sets: two of the file's groups merged into one`, cs);
  await page.evaluate(() => window.SodaSets._unmerge(0)); await settle(page);

  // ---- 3. a link to the view: a fingerprint, the file asked for, the colouring back --------------------------------
  const fp = await page.evaluate(() => window.SodaOwn.current(window.__spatialViewer.state.dsId, window.__spatialViewer.state.exprKey).fp);
  const link = await page.evaluate(async () => (await window.SodaViewLink.make()).url);
  check(link.length < 2400 && !IDS.slice(0, 50).some((x) => link.includes(x)), `${B}: the link stays short and holds no data (${link.length} characters)`);
  const p2 = await ctx.newPage(); watch(p2, "link");
  await p2.goto(link);
  await p2.waitForFunction(() => { const sv = window.__spatialViewer; return sv && sv.state && sv.state.m && sv.own && !(window.SodaViewLink && window.SodaViewLink.pending()) && document.getElementById("status").hidden; }, null, { timeout: 120000 });
  await settle(p2);
  const notes = await p2.evaluate(() => window.SodaViewLink.notes());
  const pend = await p2.evaluate(() => ({ p: window.SodaOwn._state().pending, msg: document.getElementById("ownMsg").textContent, key: window.__spatialViewer.state.exprKey }));
  check(pend.key === KEY && pend.p && pend.p.fp === fp && pend.msg.includes(fp) && notes.some((t) => t.includes(fp)), `${B}: a link with own results: it names the file's fingerprint and asks for it`, { notes, pend });
  await load(p2, FIX.scanpy);                                  // another file: the link still waits for its own
  check(!!(await p2.evaluate(() => window.SodaOwn._state().pending)), `${B}: another file does not count as the link's`);
  await load(p2, FIX.seurat);
  await p2.waitForFunction(() => window.__spatialViewer.state.umap.clRes === "own:cell_type", null, { timeout: 20000 }).catch(() => null);
  u = await um(p2);
  check(u.res === "own:cell_type" && u.colorBy === "cluster" && !(await p2.evaluate(() => window.SodaOwn._state().pending)), `${B}: the same file loaded: the link's grouping and colouring come back`, u);
  await p2.close();

  // ---- 4. Save JSON with the table in it; Load JSON on a fresh page -------------------------------------------------
  const saveTxt = await page.evaluate(() => { const b = document.getElementById("ownSave"); b.checked = true; b.dispatchEvent(new Event("change")); return document.getElementById("ownSaveTxt").textContent; });
  check(/about \d+(\.\d)? (kB|MB)/.test(saveTxt), `${B}: Save JSON shows the size of the table before it is put in`, saveTxt);
  await page.evaluate(() => { window.__blobs = []; const orig = URL.createObjectURL.bind(URL); URL.createObjectURL = (x) => { window.__blobs.push(x); return orig(x); }; document.getElementById("lassoSaveJson").click(); });
  await page.waitForFunction(() => window.__blobs.some((x) => /json/.test(x.type)), null, { timeout: 20000 });
  const saved = await page.evaluate(async () => (await window.__blobs.filter((x) => /json/.test(x.type)).pop().text()));
  const sj = JSON.parse(saved);
  const ob = sj.ownResults || {};
  check(ob.rows && ob.rows.length === matched && ob.groups.length === 2 && ob.values.length === 1 && ob.emb && ob.layer === KEY && !IDS.slice(0, 2000).some((x) => saved.includes(x)) && !saved.includes("seurat_results"),
    `${B}: Save JSON holds the matched table per bin row (no ids, no file name)`, { keys: Object.keys(ob), rows: ob.rows && ob.rows.length });
  const p3 = await open(`${DATA}/export_user`, DS, null, "loadjson");
  await p3.setInputFiles("#lassoLoad", { name: "view.json", mimeType: "application/json", buffer: Buffer.from(saved, "utf-8") });
  await p3.waitForFunction(() => { const S = window.__spatialViewer.state; return S.exprKey === "016um" && S.umap.clRes === "own:cell_type" && S.umap.cl; }, null, { timeout: 60000 }).catch(() => null);
  await settle(p3);
  await p3.evaluate(() => { if (window.__spatialViewer.state.imgOn) document.getElementById("btnImage").click(); });
  await settle(p3);
  u = await um(p3);
  b = await bandsBy(p3, "tissue", HUE);
  check(u.res === "own:cell_type" && u.embOwn && b.every((x) => x.n >= 10 && x.share >= 0.85), `${B}: Load JSON on a fresh page: the file's groups, colours and UMAP are back (by pixel)`, { u, b });
  await shot(p3, "own_loadjson.png");
  await p3.close();
  await page.close();

  // ---- 5. scanpy TSV.gz with an uns-style colour list; a partial match; roles put right ---------------------------
  page = await open(`${DATA}/export_user`, DS, KEY, "scanpy");
  await load(page, FIX.scanpy);
  m = await ownMsg(page); r = await roles(page); u = await um(page);
  check(r && r[""] === "id" && r.leiden === "cat" && r.total_counts === "num" && r.X_umap1 === "x" && r.X_umap2 === "y" && m.text.includes(`${N.toLocaleString()} of ${N.toLocaleString()} rows matched`) && /ending like "-1"/.test(m.text),
    `${B}: scanpy TSV.gz: unpacked, "-1-0" endings matched, roles guessed`, { r, m: m.text.slice(0, 240) });
  // leiden is (band + 1) % 3: "0" = right band, "1" = left, "2" = middle; an uns list in group order
  await page.fill(".own-pal[data-j='1']", "['#3050ff', '#ff2020', '#20e020']");
  await page.press(".own-pal[data-j='1']", "Tab"); await settle(page);
  b = await bandsBy(page, "tissue", HUE);
  u = await um(page);
  check(u.res === "own:leiden" && b.every((x) => x.share >= 0.85), `${B}: an uns-style colour list pasted for leiden colours its groups (by pixel)`, { b, u, cols: await page.evaluate(() => { const S = window.__spatialViewer.state, R = window.SodaOwn.current(S.dsId, S.exprKey); return R && R.clusterings.map((c) => c.colors); }) });
  await shot(page, "own_scanpy_uns.png");
  await load(page, FIX.partial);
  m = await ownMsg(page);
  const pm = FIX.partial.matched;
  const usable = await page.evaluate(() => !document.getElementById("ownUse").hidden);
  check(m.warn && m.text.includes(`${pm.toLocaleString()} of ${N.toLocaleString()} rows matched`) && m.text.includes("not listed") && usable,
    `${B}: a partial match (${Math.round((100 * pm) / N)}%): the rate, a warning, and Use anyway`, m.text.slice(0, 300));
  await page.click("#ownUse").catch(async () => page.evaluate(() => document.getElementById("ownUse").click()));
  await page.waitForFunction(() => window.__spatialViewer.state.umap.clRes === "own:cluster", null, { timeout: 20000 }).catch(() => null);
  u = await um(page);
  check(u.res === "own:cluster", `${B}: Use anyway puts the matched rows to use`, { u, s: await page.evaluate(() => { const s = window.SodaOwn._state(); return { applied: s.D && s.D.applied, force: s.D && s.D.force, store: s.store }; }) });
  await load(page, FIX.wrong);
  r = await roles(page);
  check(r.cellname === "id" && r.dim_a === "num" && r.level === "cat", `${B}: a file whose names say little: the guesses to correct`, r);
  await page.selectOption("#ownRoles select[data-j='1']", "x");
  await page.selectOption("#ownRoles select[data-j='2']", "y");
  await page.selectOption("#ownRoles select[data-j='3']", "num");
  await settle(page);
  const st = await page.evaluate(() => { const S = window.__spatialViewer.state, R = window.SodaOwn.current(S.dsId, S.exprKey); return { emb: !!(R && R.emb), nums: R ? R.nums.map((x) => x.id) : [], cl: R ? R.clusterings.length : -1, paint: [...document.getElementById("paintCol").options].map((o) => o.value) }; });
  check(st.emb && st.nums.includes("own:level") && st.cl === 0 && st.paint.includes("own:level"), `${B}: roles corrected in the table: the UMAP from dim_a / dim_b, level a value`, st);
  await page.evaluate(() => window.__spatialViewer.umapSetView("umap"));
  await page.waitForFunction(() => { const U = window.__spatialViewer.state.umap; return U.view === "umap" && !U.anim && U.embOwn; }, null, { timeout: 30000 }).catch(() => null);
  check((await um(page)).embOwn, `${B}: and that UMAP is shown`);
  await shot(page, "own_corrected_umap.png");
  await page.evaluate(() => window.__spatialViewer.umapSetView("tissue"));
  await settle(page);

  // ---- 6. drag and drop onto the page --------------------------------------------------------------------------------
  await page.evaluate(() => document.getElementById("ownClear").click()); await settle(page);
  const dropped = await page.evaluate(async (txt) => {
    const dt = new DataTransfer(); dt.items.add(new File([txt], "dropped.csv", { type: "text/csv" }));
    const fire = (type) => document.body.dispatchEvent(new DragEvent(type, { dataTransfer: dt, bubbles: true, cancelable: true }));
    fire("dragenter"); fire("dragover");
    const shown = !document.getElementById("ownDrop").hidden;
    fire("drop");
    return { shown, hidden: document.getElementById("ownDrop").hidden };
  }, FIX.seurat.buffer.toString("utf-8"));
  await page.waitForFunction(() => window.__spatialViewer.state.umap.clRes === "own:seurat_clusters", null, { timeout: 30000 }).catch(() => null);
  check(dropped.shown && dropped.hidden && (await um(page)).res === "own:seurat_clusters", `${B}: a file dropped on the page: the drop target shows, the file is loaded`, dropped);

  // ---- 7. junk files: a clear message each, nothing applied, no page error ----------------------------------------
  for (const j of JUNK) {
    await load(page, j);
    m = await ownMsg(page);
    check(m.err && j.want.test(m.text), `${B}: junk ${j.name}: "${m.text.slice(0, 110)}"`, m);
  }
  u = await um(page);
  check(u.res === "own:seurat_clusters", `${B}: after the junk, the results in use are still the last good file's`, u);
  await page.close();

  // ---- 8. no barcode list: how to add one; a "# soda row order" file; a layer without SODA's UMAP ------------------
  page = await open(`${DATA}/noids`, DS, KEY, "noids");
  await load(page, FIX.seurat);
  m = await ownMsg(page);
  check(m.err && m.text.includes("prepare.py ids --id demo_visiumhd") && m.text.includes("soda row order"), `${B}: no id list: says how to add one and what else works`, m.text.slice(0, 260));
  const order = Buffer.from(["# soda row order", "cell_type,score"].concat(Array.from({ length: N }, (_, rr) => `${NAMES[band(rr)]},${scoreOf(rr).toFixed(2)}`)).join("\n") + "\n");
  await load(page, { name: "in_order.csv", mimeType: "text/csv", buffer: order });
  m = await ownMsg(page); u = await um(page);
  check(m.warn && /row order, not by barcode/.test(m.text) && u.res === "own:cell_type", `${B}: a file in SODA's row order (said so in its first line): used, with a warning`, m.text.slice(0, 200));
  const short = Buffer.from(["# soda row order", "cell_type"].concat(Array.from({ length: N - 5 }, (_, rr) => NAMES[band(rr)])).join("\n") + "\n");
  await load(page, { name: "short.csv", mimeType: "text/csv", buffer: short });
  m = await ownMsg(page);
  check(m.err, `${B}: row order with a row count that is not the layer's: refused`, m.text.slice(0, 200));
  await page.close();
  if (fs.existsSync(path.join(ROOT, "export_loose", DS, "manifest.json"))) {
    page = await open(`${DATA}/export_loose`, DS, KEY, "noumap");
    const hadUmap = await page.evaluate(() => !document.getElementById("panelUmap").hidden);
    const lpos = f32(path.join(ROOT, "export_loose", DS, "expr", KEY, "positions.bin"));
    const same = lpos.length === POS.length && lpos.every((v, i) => v === POS[i]);
    const rowFile = Buffer.from(["soda_row,cell_type,cell_type_color,UMAP_1,UMAP_2"].concat(Array.from({ length: N }, (_, rr) => { const [ux, uy] = umapOf(rr); return `${rr},${NAMES[band(rr)]},${COLORS[band(rr)]},${ux},${uy}`; })).join("\n") + "\n");
    await load(page, { name: "rows.csv", mimeType: "text/csv", buffer: rowFile });
    m = await ownMsg(page); u = await um(page);
    check(same && !hadUmap && m.warn && /row number/.test(m.text) && u.res === "own:cell_type" && u.embOwn, `${B}: a layer without SODA's UMAP or ids: a soda_row file brings the Clusters panel, groups and a UMAP`, { hadUmap, m: m.text.slice(0, 160), u });
    await settle(page);
    b = await bandsBy(page, "tissue", HUE);
    check(b.every((x) => x.share >= 0.85), `${B}: ... and its groups colour the tissue (by pixel)`, b);
    await shot(page, "own_noumap.png");
    await page.close();
  }
  check(errors.length === 0, `${B}: no console error, page error or failed request`, errors.slice(0, 6));
  await ctx.close();

  // ---- 9. 375 px wide (a phone held upright) ---------------------------------------------------------------------
  const pctx = await browser.newContext({ viewport: { width: 375, height: 812 }, serviceWorkers: "block", hasTouch: true });
  const pp = await pctx.newPage(); const perr = []; pp.on("pageerror", (e) => perr.push(e.message));
  await pp.goto(`${URL0}?static=${encodeURIComponent(`${DATA}/export_user`)}&ds=${DS}`);
  await pp.waitForFunction(() => { const sv = window.__spatialViewer; return sv && sv.state && sv.state.m && sv.own && document.getElementById("status").hidden; }, null, { timeout: 120000 });
  await pp.evaluate(([k]) => { const el = document.getElementById("exprSelect"); el.value = k; el.dispatchEvent(new Event("change")); }, [KEY]);
  await pp.waitForFunction((k) => window.__spatialViewer.state.exprKey === k && document.getElementById("status").hidden, KEY, { timeout: 60000 });
  await pp.setInputFiles("#ownFile", { name: FIX.seurat.name, mimeType: FIX.seurat.mimeType, buffer: FIX.seurat.buffer });
  await pp.waitForFunction(() => /rows matched/.test(document.getElementById("ownMsg").textContent), null, { timeout: 60000 });
  await pp.evaluate(() => window.SodaLayout && window.SodaLayout.show("panelOwn"));
  await pp.waitForTimeout(700);
  const fit = await pp.evaluate(() => {
    const W = window.innerWidth, sec = document.getElementById("panelOwn"), r = sec.getBoundingClientRect();
    const over = [...sec.querySelectorAll("select, input, button, .own-col")].filter((e) => e.offsetParent).map((e) => e.getBoundingClientRect()).filter((b) => b.right > W + 1 || b.left < -1).length;
    return { W, doc: document.documentElement.scrollWidth, right: r.right, over, roles: sec.querySelectorAll(".own-role select").length };
  });
  check(fit.doc <= fit.W + 1 && fit.right <= fit.W + 1 && fit.over === 0 && fit.roles === 8 && !perr.length, `${B}: 375 px wide: the section and its role table fit the screen`, { fit, perr });
  await shot(pp, "own_phone.png");
  await pctx.close();
}

(async () => {
  const data = await serve(ROOT, true), pageSrv = await serve(PAGE, false);
  const DATA = `http://127.0.0.1:${data.address().port}`, URL0 = `http://127.0.0.1:${pageSrv.address().port}/index.html`;
  for (const name of BROWSERS) {
    if (!pw[name]) { console.log(`${name}: not in this Playwright, skipped`); continue; }
    let b = null;
    try { b = await pw[name].launch(name === "chromium" ? { args: GPU ? GPU_ARGS : SOFT_ARGS } : {}); }
    catch (e) { console.log(`${name}: not installed (${String(e.message || e).split("\n")[0]}), skipped`); continue; }
    try { await run(name, b, URL0, DATA); }
    catch (e) { check(false, `${name}: the checks ran to the end`, String(e && e.stack || e).slice(0, 600)); }
    await b.close();
  }
  data.close(); pageSrv.close();
  const bad = results.filter((x) => !x.ok);
  console.log(bad.length ? `\nFAILED: ${bad.length} of ${results.length} checks` : `\nPASSED: ${results.length} checks`);
  process.exit(bad.length ? 1 : 0);
})().catch((e) => { console.error(e); process.exit(2); });
