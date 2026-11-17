#!/usr/bin/env node
// Old exports test: the CURRENT page must still open folders exported by older versions of
// prepare.py, because those folders are already on people's storage.
//
// tests/fixtures/old_exports/<version>/<id>/ holds small synthetic exports written by old code
// (tests/fixtures/old_exports/make_old_exports.py rebuilds them from git history):
//   v1-20261002  first public version: tiles, expr (no region.bin / hvg / umap), overlays, an RNA-seq set
//   v2-20261104  regions, hvg/, umap/ (no clustering names), first protein grid (mat/ with p1 / p99_5)
//   v2-20261111  protein grid with bg-asinh-v1 + mean/sd, qc.json; before binqc/ and user clusters
// The page (web/SpatialViewer, served as it is) is opened in headless Chromium against each one:
//   static   ?static=<folder>&ds=<id>      every dataset: no console error, no failed request, the
//            dataset in the list, the channels of the manifest, image tiles fetched, decoded and
//            drawn, the layer's genes listed, one gene added through the gene search (file fetched,
//            decoded, the right number of points), clusters / UMAP and region ids when the export
//            has them, every overlay switched on and read; RNA-seq: UMAP, colour-by, one gene.
//   data     ?data=<folder>/&ds=<id>       a shared link to an old folder with no datasets.json
//   catalog  ?data=<folder>/               with --python: datasets.json written by the current
//            catalog.py over the old manifests; the list shows every dataset and each one opens.
//
// With --packed <export root> it also opens CURRENT exports, whose gene files are packed
// (web/SpatialViewer/genepack.py), with the same checks; there a gene must arrive as one 206 answer
// to a Range request into a pack. On top of that, for those exports:
//   sas      ?data=<root>/?sv=..&sig=..&ds=<id>  a data link with a read pass: the pass rides on the
//            Range request too, and the gene is right
//   norange  ?static=<the same folder on a server that ignores Range>  the whole shard comes back
//            (200) and the gene is cut out of it, right; with config packWholeMaxMB below the
//            shard size the gene is refused with a message that names Range, and nothing else breaks
// Its own line, "packed exports (...): N passed, M failed", comes before the old exports line.
//
//   node tests/compat_test.js [--fixtures <dir>] [--only v1-20261002] [--python <python>]
//                             [--packed <export root>] [--packed-ids demo_visiumhd,demo_xenium]
//                             [--ports 8890,8891] [--shots <dir>] [--playwright <node_modules>]
//
// Playwright is not a dependency of this repo (see tests/phone_browser_test.js): install it
// anywhere and point NODE_PATH or --playwright at its node_modules (PLAYWRIGHT_BROWSERS_PATH if
// its browsers are elsewhere). Without it, or without the fixtures, the test says SKIPPED and
// exits 0. deck.gl comes from the CDN, as on the real page: no internet = SKIPPED too.
"use strict";
const http = require("http");
const fs = require("fs");
const os = require("os");
const path = require("path");
const { spawnSync } = require("child_process");

const arg = (k, d) => { const i = process.argv.indexOf(k); return i > 0 ? process.argv[i + 1] : d; };
const REPO = path.resolve(__dirname, "..");
const FIX = path.resolve(arg("--fixtures", path.join(REPO, "tests", "fixtures", "old_exports")));
const PAGE = path.resolve(arg("--page", path.join(REPO, "web", "SpatialViewer")));
const ONLY = arg("--only", "");
const PYTHON = arg("--python", "");
const SHOTS = arg("--shots", "");
const PORTS = arg("--ports", "8890,8891").split(",").map(Number);    // data server, page server
const PACKED = arg("--packed", "") ? path.resolve(arg("--packed", "")) : "";
const PACKED_IDS = arg("--packed-ids", "");
const OPEN_MS = 30000;

let pw = null;
for (const base of [arg("--playwright", ""), ...(process.env.NODE_PATH || "").split(path.delimiter), path.join(REPO, "node_modules")].filter(Boolean)) {
  try { pw = require(path.join(base, "playwright")); break; } catch (e) { /* next */ }
}
if (!pw) { try { pw = require("playwright"); } catch (e) { /* none */ } }
if (!pw) { console.log("SKIPPED: playwright is not installed (see the top of this file)"); process.exit(0); }

// <version>/<id>/manifest.json (any other folder, e.g. __pycache__ of the generator, is not a version)
const datasetsOf = (v) => fs.readdirSync(path.join(FIX, v)).filter((id) => fs.existsSync(path.join(FIX, v, id, "manifest.json"))).sort();
const versions = fs.existsSync(FIX) ? fs.readdirSync(FIX).filter((v) => fs.statSync(path.join(FIX, v)).isDirectory())
  .filter((v) => (!ONLY || ONLY.split(",").includes(v)) && datasetsOf(v).length).sort() : [];
// current (packed) exports to open: those named, else every dataset in the folder with a genePack
const packedIds = !PACKED || !fs.existsSync(PACKED) ? [] : fs.readdirSync(PACKED).filter((id) => fs.existsSync(path.join(PACKED, id, "manifest.json")))
  .filter((id) => (PACKED_IDS ? PACKED_IDS.split(",").includes(id) : true)).sort();
if (!versions.length) {
  console.log(`SKIPPED: no old exports in ${FIX} (rebuild them with tests/fixtures/old_exports/make_old_exports.py)`); process.exit(0);
}
const readJson = (...p) => JSON.parse(fs.readFileSync(path.join(...p), "utf-8").replace(/\bNaN\b/g, "null"));

// ---- datasets.json from the current catalog.py, written next to copies of the old manifests ----
let catDir = null;
const catFail = [];
if (PYTHON) {
  catDir = fs.mkdtempSync(path.join(os.tmpdir(), "soda-compat-"));
  for (const v of versions) for (const id of datasetsOf(v)) {
    fs.mkdirSync(path.join(catDir, v, id), { recursive: true });
    fs.copyFileSync(path.join(FIX, v, id, "manifest.json"), path.join(catDir, v, id, "manifest.json"));
  }
  for (const v of versions) {
    const p = spawnSync(PYTHON, [path.join(PAGE, "catalog.py"), path.join(catDir, v)], { encoding: "utf-8" });
    if (p.status !== 0) { catFail.push(`${v} | catalog | catalog.py on the old manifests exited ${p.status}: ${(p.stderr || p.stdout || "").trim().slice(-300)}`); console.log(`  FAIL ${catFail[catFail.length - 1]}`); }
  }
}

// ---- two small static servers: the exports (with CORS and Range, as blob storage) and the page ----
const { sendFile, preflight } = require("./static_files");
const TYPES = { ".html": "text/html", ".js": "text/javascript", ".css": "text/css", ".json": "application/json", ".png": "image/png", ".jpg": "image/jpeg", ".bin": "application/octet-stream", ".geojson": "application/json" };
const seen = [];                  // every data request: url, Range header, status (the packed checks read it)
function handler(route) {
  return (req, res) => {
    if (req.method === "OPTIONS") { preflight(req, res); return; }
    const u = decodeURIComponent(req.url.split("?")[0]);
    const file = route(u);
    if (!file || !fs.existsSync(file) || fs.statSync(file).isDirectory()) { res.writeHead(404, { "Access-Control-Allow-Origin": "*" }); res.end("not found"); return; }
    const noRange = /^\/norange\//.test(u);              // a host that ignores Range
    res.on("finish", () => seen.push({ url: req.url, range: req.headers.range || "", status: res.statusCode }));
    sendFile(req, res, file, TYPES[path.extname(file)], { range: !noRange });
  };
}
const inside = (root, rel) => { const f = path.join(root, rel); return f.startsWith(root + path.sep) || f === root ? f : null; };
// /plain/<version>/...  the old folder as it is (no datasets.json)
// /cat/<version>/...    the same, plus the datasets.json catalog.py wrote
const dataRoute = (u) => {
  const pk = u.match(/^\/(packed|norange)(\/.*)$/);           // the current export root (--packed)
  if (pk) return PACKED ? inside(PACKED, pk[2]) : null;
  const m = u.match(/^\/(plain|cat)(\/.*)$/);
  if (!m) return null;
  if (m[1] === "cat" && catDir && /\/datasets\.json$/.test(m[2])) return inside(catDir, m[2]);
  return inside(FIX, m[2]);
};
const pageRoute = (u) => inside(PAGE, u === "/" ? "index.html" : u);
function listen(h, port) {
  return new Promise((resolve, reject) => {
    const srv = http.createServer(h);
    srv.once("error", (e) => {
      if (e.code !== "EADDRINUSE" || !port) return reject(e);
      console.log(`  (port ${port} is taken: using a free one)`);
      listen(h, 0).then(resolve, reject);
    });
    srv.listen(port || 0, "127.0.0.1", () => resolve(srv));
  });
}

// what the page holds, read from its test hooks (window.__spatialViewer)
const look = (page) => page.evaluate(() => {
  const sv = window.__spatialViewer, S = sv && sv.state, st = document.getElementById("status"), out = {};
  out.status = st && !st.hidden ? st.textContent : "";
  out.err = !!(st && !st.hidden && st.classList.contains("err"));
  out.deck = !!window.deck;
  if (!S) return out;
  out.id = S.m && S.m.id;
  out.picker = sv.dsPicker ? { n: sv.dsPicker.count(), value: sv.dsPicker.value() } : null;
  out.datasets = (S.datasets || []).map((d) => d.id);
  out.channels = (S.channels || []).map((c) => c.name);
  out.cached = sv.diag ? sv.diag().cached : 0;
  out.exprKey = S.exprKey || null;
  out.geneIndex = (S.geneIndex || []).length;
  out.genes = (S.genes || []).map((g) => ({ name: g.name, n: g.n }));
  out.regionIds = S.region && S.region.ids ? S.region.ids.length : null;
  out.overlays = (S.overlays || []).map((o) => ({ file: o.file, on: o.on, nFeat: o.data && o.data.features ? o.data.features.length : null }));
  const U = S.umap;
  out.umap = U && U.uxy ? { nEmb: U.nEmb, rows: U.uxy.length / 2, cl: U.cl ? U.cl.length : null, view: U.view } : null;
  const R = S.rna;
  out.rna = R ? { n: R.n, index: (R.index || []).length, colorBy: R.colorBy, cat: !!R.cat, genes: R.genes.map((g) => ({ name: g.name, k: g.k })) } : null;
  out.lit = -1;
  if (sv.deck && sv.deck.device && !sv.deck.device.gl.isContextLost()) {
    const gl = sv.deck.device.gl; sv.deck.redraw("test");
    const w = gl.drawingBufferWidth, h = gl.drawingBufferHeight, px = new Uint8Array(w * h * 4);
    gl.bindFramebuffer(gl.FRAMEBUFFER, null); gl.readPixels(0, 0, w, h, gl.RGBA, gl.UNSIGNED_BYTE, px);
    let lit = 0; for (let i = 0; i < px.length; i += 4) if (px[i] > 8 || px[i + 1] > 8 || px[i + 2] > 8) lit++;
    out.lit = lit / (w * h);
  }
  return out;
});
const until = (page, fn, arg, ms = OPEN_MS) => page.waitForFunction(fn, arg, { timeout: ms, polling: 200 }).then(() => true, () => false);
const settled = (page) => until(page, () => { const st = document.getElementById("status"); return !st || st.hidden || !/Loading|Opening/.test(st.textContent); }, null, 15000);
// type into a search box and click the first hit, as a reader would (works on a folded panel too)
const pickFromSearch = (page, box, results, name) => page.evaluate(([box, results, name]) => {
  const el = document.getElementById(box);
  el.value = name; el.dispatchEvent(new Event("input", { bubbles: true }));
  const hit = document.querySelector(`#${results} > div`);
  if (!hit) return null;
  const shown = hit.querySelector("span") ? hit.querySelector("span").textContent : "";
  hit.dispatchEvent(new MouseEvent("mousedown", { bubbles: true, cancelable: true }));
  return shown;
}, [box, results, name]);

const fileSize = (...p) => { const f = path.join(...p); return fs.existsSync(f) ? fs.statSync(f).size : -1; };
const encName = (n) => encodeURIComponent(n);
// a gene's bytes in a packed group: the slice of its shard (genepack.py), or -1 / null when none
const packOf = (pack, ...dir) => pack ? { pack, dir: path.join(...dir) } : null;
function packSlice(pk, e) {
  if (!pk || !e || !(e.b > 0)) return null;
  const fd = fs.openSync(path.join(pk.dir, pk.pack.dir || "genepack", pk.pack.shards[e.s].file), "r");
  try { const b = Buffer.alloc(e.b); fs.readSync(fd, b, 0, e.b, e.o); return b; } finally { fs.closeSync(fd); }
}
// the layer's first gene that has a file (genes with no counts have none), by total; its points
function geneToAdd(dir, L) {
  const genes = readJson(dir, "expr", L.key, "genes.json");
  const bytes = (g) => (L.genePack ? (g.b > 0 ? g.b : -1) : fileSize(dir, "expr", L.key, "gene", `${encName(g.n)}.bin`));
  const withFile = genes.filter((g) => bytes(g) > 0);
  withFile.sort((a, b) => b.t - a.t);
  return { genes, gene: withFile[0], pts: withFile[0] ? bytes(withFile[0]) / 12 : -1 };
}

let pass = 0, fail = catFail.length;
const failures = [...catFail];
const packedTally = { pass: 0, fail: 0 };

async function openPage(browser, url, version, label, config) {
  const ctx = await browser.newContext({ viewport: { width: 1280, height: 800 } });
  const page = await ctx.newPage();
  const log = { console: [], pageErrors: [], bad: [], ok: [] };
  page.on("console", (m) => { if (m.type() === "error") log.console.push({ text: m.text(), url: (m.location() || {}).url || "" }); });
  page.on("pageerror", (e) => log.pageErrors.push(String(e.message || e).split("\n")[0]));
  page.on("requestfailed", (r) => log.bad.push({ url: r.url(), status: r.failure() ? r.failure().errorText : "failed" }));
  page.on("response", (r) => (r.status() >= 400 ? log.bad : log.ok).push({ url: r.url(), status: r.status() }));
  // the repo has no config.js (only config.example.js): an empty one, as on a fresh copy
  await page.route("**/config.js*", (rt) => rt.fulfill({ status: 200, contentType: "text/javascript", body: `window.SpatialViewerConfig = ${JSON.stringify(config || {})};` }));
  await page.goto(url, { waitUntil: "load" });
  const cur = version === "current";
  const ok = (cond, what) => {
    if (cur) { if (cond) packedTally.pass++; else packedTally.fail++; } else if (cond) pass++; else fail++;
    if (!cond) failures.push(`${version} | ${label} | ${what}`);
    console.log(`${cond ? "  ok  " : "  FAIL"} ${version} | ${label} | ${what}`);
  };
  return { ctx, page, log, ok };
}
// no console error, no page error, no failed request except the ones `allowed` names
// Files the page asks for that no export is required to have, current ones included (the browser
// still logs their 404 as a console error): <dataset>/gates.json, a lab's published gates
// (viewer.js gateFetchPublished, "absent = none, quietly"), and the page's own markers/<file> named by
// a grid layer's celltypes markerSet (the fixtures were exported with examples/demo_protein_markers.json,
// which the page does not ship; no list = no suggested names, as for a current export made that way).
const OPTIONAL = [/\/gates\.json$/, /\/markers\/[\w.-]+\.json$/];
function checkClean(t, allowed) {
  // a data link's read pass (?sv=...&sig=...) follows every file name: match the path without it too
  const allow = (u) => allowed.concat(OPTIONAL).some((re) => re.test(u) || re.test(String(u).split("?")[0]));
  const bad = t.log.bad.filter((b) => !allow(b.url));
  const cons = t.log.console.filter((c) => !(/Failed to load resource/.test(c.text) && c.url && allow(c.url)));
  t.ok(!t.log.pageErrors.length, `no page error${t.log.pageErrors.length ? ": " + t.log.pageErrors.slice(0, 3).join(" / ") : ""}`);
  t.ok(!cons.length, `no console error${cons.length ? ": " + cons.slice(0, 3).map((c) => `${c.text} ${c.url}`).join(" / ") : ""}`);
  t.ok(!bad.length, `no failed request${bad.length ? ": " + bad.slice(0, 5).map((b) => `${b.status} ${b.url.replace(/^https?:\/\/[^/]+/, "")}`).join(", ") : ""}`);
}
const fetched = (t, re, status = 200) => t.log.ok.some((r) => r.status === status && re.test(r.url));
async function shot(t, version, label) {
  if (!SHOTS) return;
  fs.mkdirSync(SHOTS, { recursive: true });
  await t.page.screenshot({ path: path.join(SHOTS, `${version}_${label.replace(/[^\w]+/g, "_")}.png`) });
}

// the full set of checks on one dataset, opened with ?static= (an old fixture, or with version
// "current" a dataset of the --packed export root)
async function checkDataset(browser, urls, version, id) {
  const cur = version === "current";
  const dir = cur ? path.join(PACKED, id) : path.join(FIX, version, id), m = readJson(dir, "manifest.json");
  const label = `static ${id}`;
  const t = await openPage(browser, `${urls.page}/index.html?static=${urls.data}/${cur ? "packed" : `plain/${version}`}&ds=${id}`, version, label);
  try {
    const opened = await until(t.page, (id) => { const sv = window.__spatialViewer; return sv && sv.state && sv.state.m && sv.state.m.id === id; }, id);
    let r = await look(t.page);
    if (!r.deck) { t.ok(false, "deck.gl did not load (no internet?)"); return; }
    t.ok(opened, `manifest read (kind ${m.kind}${r.status ? `; status "${r.status}"` : ""})`);
    if (!opened) return;
    t.ok(r.picker && r.picker.n >= 1 && r.picker.value === id && r.datasets.includes(id), `dataset list shows ${id} (${r.picker ? r.picker.n : 0} listed, picked ${r.picker && r.picker.value})`);
    if (m.kind === "rnaseq") {
      await until(t.page, () => { const R = window.__spatialViewer.state.rna; return R && R.index && R.index.length; });
      await settled(t.page);
      const genes = readJson(dir, "genes.json");
      r = await look(t.page);
      t.ok(r.rna && r.rna.n === m.nCells, `UMAP read: ${r.rna && r.rna.n} cells (manifest ${m.nCells})`);
      t.ok(r.rna && r.rna.index === genes.length, `genes listed: ${r.rna && r.rna.index} (genes.json ${genes.length})`);
      const cat = (m.obs || []).find((o) => o.type === "category");
      t.ok(!cat || (r.rna.colorBy === cat.name && r.rna.cat), `coloured by "${r.rna && r.rna.colorBy}" (obs/${cat ? cat.file : "-"} read)`);
      const pk = packOf(m.genePack, dir);
      const g = genes[0].n, k = (pk ? packSlice(pk, genes[0]) : fs.readFileSync(path.join(dir, "gene", `${encName(g)}.bin`))).readUInt32LE(0);
      const shown = await pickFromSearch(t.page, "rnaSearch", "rnaResults", g);
      await until(t.page, () => window.__spatialViewer.state.rna.genes.length > 0, null, 15000);
      r = await look(t.page);
      const got = r.rna.genes[0];
      t.ok(shown === g && got && got.name === g && got.k === k, `gene ${g} added from the search: ${got ? got.k : "none"} cells (file ${k})`);
      t.ok(fetched(t, /\/umap\.bin$/) && (pk ? fetched(t, /\/genepack\/[0-9a-f]+-\d+\.pack$/, 206) : fetched(t, new RegExp(`/gene/${g}\\.bin$`))),
        pk ? "umap.bin fetched, the gene by a Range request into the pack (206)" : "umap.bin and the gene file fetched");
      await shot(t, version, label);
      t.ok(r.lit > 0.002, `UMAP drawn (lit ${r.lit.toFixed(4)})`);
      t.ok(!r.err, `no error line (${r.status || "none"})`);
    } else {
      const want = m.image.channels.map((c) => c.name).filter((n) => !/^\s*(open|blank|empty|unused|none)([\s_-]|$)/i.test(n));
      t.ok(r.channels.length === want.length && want.every((n) => r.channels.includes(n)), `channels listed: ${r.channels.length} (${r.channels.slice(0, 4).join(", ")}${r.channels.length > 4 ? ", ..." : ""}; manifest ${want.length})`);
      const tiles = await until(t.page, () => window.__spatialViewer.diag().cached > 0);
      await settled(t.page);
      const fmt = m.tileFormat || "png";
      t.ok(tiles && fetched(t, new RegExp(`/tiles/\\d+/\\d+/\\d+/\\d+\\.${fmt}$`)), `image tiles fetched and decoded (.${fmt})`);
      const ex = m.expr || [];
      if (ex.length) {
        await until(t.page, () => { const S = window.__spatialViewer.state; return S.exprKey && S.geneIndex && S.geneIndex.length; });
        await settled(t.page);
        r = await look(t.page);
        const L = ex.find((e) => e.key === r.exprKey);
        t.ok(!!L, `expression layer opened: ${r.exprKey} (manifest: ${ex.map((e) => e.key).join(", ")})`);
        if (L) {
          const { genes, gene, pts } = geneToAdd(dir, L);
          t.ok(r.geneIndex === genes.length, `${L.kind === "protein" ? "markers" : "genes"} listed: ${r.geneIndex} (genes.json ${genes.length})`);
          const shown = await pickFromSearch(t.page, "geneSearch", "geneResults", gene.n);
          await until(t.page, () => window.__spatialViewer.state.genes.length > 0, null, 15000);
          await settled(t.page);
          r = await look(t.page);
          const got = r.genes[0];
          t.ok(shown === gene.n && got && got.name === gene.n && got.n === pts, `gene ${gene.n} added from the search: ${got ? got.n : "none"} points (file ${pts})`);
          if (L.genePack) t.ok(fetched(t, new RegExp(`/expr/${L.key}/genepack/[0-9a-f]+-\\d+\\.pack$`), 206), `expr/${L.key}: the gene read by a Range request into the pack (206)`);
          else t.ok(fetched(t, new RegExp(`/expr/${L.key}/gene/${encName(gene.n)}\\.bin$`)), `expr/${L.key}/gene/${gene.n}.bin fetched`);
          if (L.regions) t.ok(r.regionIds === fileSize(dir, "expr", L.key, "positions.bin") / 8, `region.bin read: ${r.regionIds} rows`);
          if (L.umap) {
            const idx = readJson(dir, "expr", L.key, "umap", "index.json");
            const rows = fileSize(dir, "expr", L.key, "positions.bin") / 8;
            await t.page.evaluate(() => window.__spatialViewer.umapSetView("umap"));
            await until(t.page, () => { const U = window.__spatialViewer.state.umap; return U && U.cl && U.view === "umap" && !U.anim; }, null, 15000);   // the flight from the tissue has landed
            await settled(t.page);
            r = await look(t.page);
            t.ok(r.umap && r.umap.nEmb === idx.nEmbedded && r.umap.cl === rows && r.umap.view === "umap",
              `UMAP and clusters read: ${r.umap ? `${r.umap.nEmb} embedded, ${r.umap.cl} cluster rows, view ${r.umap.view}` : "none"} (index.json ${idx.nEmbedded}, ${rows} rows)`);
            t.ok(r.lit > 0.002, `UMAP drawn (lit ${r.lit.toFixed(4)})`);
            await shot(t, version, label + " umap");
            await t.page.evaluate(() => window.__spatialViewer.umapSetView("tissue"));
            await settled(t.page);
          }
        }
      }
      // every overlay switched on, as a reader would tick it
      if ((m.overlays || []).length) {
        await t.page.evaluate(() => {
          for (const cb of document.querySelectorAll("#overlayList input[type=checkbox]")) if (!cb.checked) { cb.checked = true; cb.dispatchEvent(new Event("change")); }
        });
        await until(t.page, () => window.__spatialViewer.state.overlays.every((o) => o.data), null, 15000);
        r = await look(t.page);
        t.ok(r.overlays.length === m.overlays.length && r.overlays.every((o) => o.on && o.nFeat > 0),
          `overlays listed and read: ${r.overlays.map((o) => `${o.file} (${o.nFeat} features)`).join(", ")}`);
      }
      await settled(t.page);
      await shot(t, version, label);
      r = await look(t.page);
      t.ok(r.lit > 0.01, `slide drawn (lit ${r.lit.toFixed(3)})`);
      t.ok(!r.err, `no error line (${r.status || "none"})`);
    }
    checkClean(t, []);
  } catch (e) { t.ok(false, `crashed: ${String(e.message).split("\n")[0]}`); }
  finally { await t.ctx.close(); }
}

// a current export with packed gene files: a ?data= link with a read pass, and a host that ignores Range
const SAS = "sv=2099-01-01&sr=c&sp=rl&sig=FAKE-SIG-FOR-THE-TEST";
async function addLayerGene(t, dir, m) {
  await until(t.page, () => { const S = window.__spatialViewer.state; return S.exprKey && S.geneIndex && S.geneIndex.length; });
  await settled(t.page);
  let r = await look(t.page);
  const L = (m.expr || []).find((e) => e.key === r.exprKey);
  if (!L) { t.ok(false, `expression layer opened: ${r.exprKey}`); return null; }
  const { gene, pts } = geneToAdd(dir, L);
  await pickFromSearch(t.page, "geneSearch", "geneResults", gene.n);
  await until(t.page, () => { const S = window.__spatialViewer.state, st = document.getElementById("status");
    return S.genes.length > 0 || (st && !st.hidden && st.classList.contains("err")); }, null, 15000);
  await settled(t.page);
  r = await look(t.page);
  return { L, gene, pts, r, got: r.genes[0] };
}
async function checkPacked(browser, urls, id) {
  const dir = path.join(PACKED, id), m = readJson(dir, "manifest.json");
  const packedLayers = (m.expr || []).filter((L) => L.genePack).map((L) => L.key);
  const V = "current";
  if (!(m.expr || []).length && !m.genePack && !(m.molecules && m.molecules.genePack)) return;
  // 1. a ?data= link whose files carry a read pass (SAS): the pass goes on the Range request too
  {
    const n0 = seen.length;
    const link = encodeURIComponent(`${urls.data}/packed/?${SAS}`);
    const t = await openPage(browser, `${urls.page}/index.html?data=${link}&ds=${id}`, V, `sas link ${id}`);
    try {
      const opened = await until(t.page, (id) => { const sv = window.__spatialViewer; return sv && sv.state && sv.state.m && sv.state.m.id === id; }, id);
      t.ok(opened, "opens from a ?data= link with a read pass");
      if (m.kind === "rnaseq") {
        await until(t.page, () => { const R = window.__spatialViewer.state.rna; return R && R.index && R.index.length; });
        const genes = readJson(dir, "genes.json"), e = genes.find((g) => g.b > 0);
        const k = packSlice(packOf(m.genePack, dir), e).readUInt32LE(0);
        await pickFromSearch(t.page, "rnaSearch", "rnaResults", e.n);
        await until(t.page, () => window.__spatialViewer.state.rna.genes.length > 0, null, 15000);
        const r = await look(t.page), got = r.rna.genes[0];
        t.ok(got && got.name === e.n && got.k === k, `RNA-seq gene ${e.n} from the pack: ${got ? got.k : "none"} cells (pack ${k})`);
      } else if ((m.expr || []).length) {
        const a = await addLayerGene(t, dir, m);
        if (a) {
          t.ok(a.got && a.got.name === a.gene.n && a.got.n === a.pts, `gene ${a.gene.n} added: ${a.got ? a.got.n : "none"} points (pack ${a.pts})`);
          if (a.L.genePack) {
            const hits = seen.slice(n0).filter((q) => q.url.includes(`/expr/${a.L.key}/genepack/`));
            t.ok(hits.length > 0 && hits.every((q) => q.status === 206 && /^bytes=\d+-\d+$/.test(q.range) && q.url.includes("sig=FAKE-SIG")),
              `the pack read with Range and the read pass: ${hits.map((q) => `${q.status} ${q.range}${q.url.includes("sig=") ? " +pass" : " NO PASS"}`).join(", ") || "no request"}`);
          }
        }
      }
      // molecules: the most abundant gene's dots, from the molecules pack
      if (m.molecules && m.molecules.genePack) {
        const idx = readJson(dir, m.molecules.dir || "molecules", "index.json");
        const g = idx.genes.filter((x) => x.b > 0).sort((a, b) => b.c - a.c)[0];
        await until(t.page, () => !!window.__spatialViewer.state.mol.index);
        await pickFromSearch(t.page, "molSearch", "molResults", g.n);
        await until(t.page, () => window.__spatialViewer.state.mol.genes.length > 0, null, 15000);
        const got = await t.page.evaluate(() => window.__spatialViewer.state.mol.genes.map((x) => ({ name: x.name, n: x.n })));
        t.ok(got.length && got[0].name === g.n && got[0].n === g.c && g.b === 4 * g.c, `molecules of ${g.n}: ${got.length ? got[0].n : "none"} (index ${g.c})`);
        t.ok(seen.slice(n0).some((q) => q.url.includes("/molecules/genepack/") && q.status === 206 && q.url.includes("sig=FAKE-SIG")), "molecules read with Range and the read pass");
      }
      const r = await look(t.page);
      t.ok(!r.err, `no error line (${r.status || "none"})`);
      checkClean(t, [/\/datasets\.json(\?|$)/, /\/about\.html(\?|$)/]);
    } catch (e) { t.ok(false, `crashed: ${String(e.message).split("\n")[0]}`); }
    finally { await t.ctx.close(); }
  }
  if (!packedLayers.length) return;
  // 2. a host that ignores Range: the whole shard comes back, the gene is cut out of it
  {
    const n0 = seen.length;
    const t = await openPage(browser, `${urls.page}/index.html?static=${urls.data}/norange&ds=${id}`, V, `no-Range host ${id}`);
    try {
      await until(t.page, (id) => { const sv = window.__spatialViewer; return sv && sv.state && sv.state.m && sv.state.m.id === id; }, id);
      const a = await addLayerGene(t, dir, m);
      if (a) {
        t.ok(a.got && a.got.name === a.gene.n && a.got.n === a.pts, `gene ${a.gene.n} cut out of the whole shard: ${a.got ? a.got.n : "none"} points (pack ${a.pts}${a.r.err ? `; "${a.r.status}"` : ""})`);
        const hits = seen.slice(n0).filter((q) => q.url.includes("/genepack/"));
        t.ok(hits.length > 0 && hits.every((q) => q.status === 200 && q.range), `the host answered the Range request with the whole shard: ${hits.map((q) => q.status).join(", ")}`);
      }
      checkClean(t, []);
    } catch (e) { t.ok(false, `crashed: ${String(e.message).split("\n")[0]}`); }
    finally { await t.ctx.close(); }
  }
  // 3. the same host with a shard larger than the page will download for one gene: a clear refusal
  {
    const t = await openPage(browser, `${urls.page}/index.html?static=${urls.data}/norange&ds=${id}`, V, `no-Range host, big shard ${id}`, { packWholeMaxMB: 1e-6 });
    try {
      await until(t.page, (id) => { const sv = window.__spatialViewer; return sv && sv.state && sv.state.m && sv.state.m.id === id; }, id);
      const a = await addLayerGene(t, dir, m);
      if (a) {
        t.ok(!a.r.genes.length && a.r.err && /Range/.test(a.r.status) && /--gene-files loose/.test(a.r.status), `gene refused with a message: "${a.r.status}"`);
        const r = await look(t.page);
        t.ok(r.lit > 0.01, `slide still drawn (lit ${r.lit.toFixed(3)})`);
      }
      checkClean(t, [/\/genepack\//]);                // the cancelled shard download
    } catch (e) { t.ok(false, `crashed: ${String(e.message).split("\n")[0]}`); }
    finally { await t.ctx.close(); }
  }
}

// a shared ?data= link to an old folder: no datasets.json there, so the page falls back to ?ds=
async function checkDataLink(browser, urls, version, id) {
  const m = readJson(FIX, version, id, "manifest.json");
  const t = await openPage(browser, `${urls.page}/index.html?data=${urls.data}/plain/${version}/&ds=${id}`, version, `data link ${id}`);
  try {
    const opened = await until(t.page, (id) => { const sv = window.__spatialViewer; return sv && sv.state && sv.state.m && sv.state.m.id === id; }, id);
    t.ok(opened, "opens with no datasets.json in the folder");
    if (m.kind === "rnaseq") await until(t.page, () => { const R = window.__spatialViewer.state.rna; return R && R.index; });
    else await until(t.page, () => window.__spatialViewer.diag().cached > 0);
    await settled(t.page);
    const r = await look(t.page);
    t.ok(r.picker && r.picker.value === id, `dataset list shows ${id}`);
    t.ok(!r.err, `no error line (${r.status || "none"})`);
    checkClean(t, [/\/datasets\.json$/, /\/about\.html$/]);
  } catch (e) { t.ok(false, `crashed: ${String(e.message).split("\n")[0]}`); }
  finally { await t.ctx.close(); }
}

// the old folder with a datasets.json from the current catalog.py: the list, then each dataset
async function checkCatalog(browser, urls, version) {
  const ids = datasetsOf(version);
  const t = await openPage(browser, `${urls.page}/index.html?data=${urls.data}/cat/${version}/`, version, "catalog");
  try {
    const listed = await until(t.page, (n) => { const sv = window.__spatialViewer; return sv && sv.dsPicker && sv.dsPicker.count() === n && sv.state.m; }, ids.length);
    let r = await look(t.page);
    t.ok(listed && ids.every((id) => r.datasets.includes(id)), `datasets.json lists all ${ids.length}: ${r.datasets.join(", ")}`);
    for (const id of ids) {
      const m = readJson(FIX, version, id, "manifest.json");
      await t.page.evaluate((id) => window.__spatialViewer.dsPicker.pick(id), id);
      const opened = await until(t.page, ([id, rna]) => {
        const S = window.__spatialViewer.state;
        return S.m && S.m.id === id && (rna ? S.rna && S.rna.index : window.__spatialViewer.diag().cached > 0);
      }, [id, m.kind === "rnaseq"]);
      await settled(t.page);
      r = await look(t.page);
      t.ok(opened && !r.err, `${id} picked from the list and opened${r.err ? `: "${r.status}"` : ""}`);
    }
    checkClean(t, [/\/about\.html$/]);
  } catch (e) { t.ok(false, `crashed: ${String(e.message).split("\n")[0]}`); }
  finally { await t.ctx.close(); }
}

(async () => {
  const data = await listen(handler(dataRoute), PORTS[0]), pageSrv = await listen(handler(pageRoute), PORTS[1]);
  const urls = { data: `http://127.0.0.1:${data.address().port}`, page: `http://127.0.0.1:${pageSrv.address().port}` };
  let browser;
  try { browser = await pw.chromium.launch(); }
  catch (e) { console.log(`SKIPPED: cannot launch Chromium (${String(e.message).split("\n")[0]})`); data.close(); pageSrv.close(); process.exit(0); }
  // deck.gl is loaded from the CDN by the page: without it nothing can be checked
  {
    const ctx = await browser.newContext(), page = await ctx.newPage();
    await page.route("**/config.js*", (rt) => rt.fulfill({ status: 200, contentType: "text/javascript", body: "window.SpatialViewerConfig = {};" }));
    await page.goto(`${urls.page}/index.html`, { waitUntil: "load" }).catch(() => {});
    const deck = await page.evaluate(() => !!window.deck).catch(() => false);
    await ctx.close();
    if (!deck) { console.log("SKIPPED: deck.gl did not load from the CDN (no internet?)"); await browser.close(); data.close(); pageSrv.close(); process.exit(0); }
  }
  try {
    if (packedIds.length) {
      console.log(`current (packed) exports in ${PACKED}: ${packedIds.join(", ")}`);
      for (const id of packedIds) { await checkDataset(browser, urls, "current", id); await checkPacked(browser, urls, id); }
    }
    for (const v of versions) {
      const ids = datasetsOf(v);
      console.log(`${v}: ${ids.join(", ")}`);
      for (const id of ids) await checkDataset(browser, urls, v, id);
      for (const id of ids) await checkDataLink(browser, urls, v, id);
      if (catDir) await checkCatalog(browser, urls, v);
    }
  } finally {
    await browser.close();
    data.close(); pageSrv.close();
    if (catDir) fs.rmSync(catDir, { recursive: true, force: true });     // this run's own temporary copies
  }
  if (failures.length) console.log("\nFailed:\n" + failures.map((f) => "  " + f).join("\n"));
  if (packedIds.length) console.log(`\npacked exports (${packedIds.join(", ")}): ${packedTally.pass} passed, ${packedTally.fail} failed`);
  console.log(`\nold exports (${versions.join(", ")}): ${pass} passed, ${fail} failed`);
  process.exit(fail || packedTally.fail ? 1 : 0);
})();
