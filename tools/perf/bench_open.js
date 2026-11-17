"use strict";
// How long SODA takes to open a dataset, show a gene, and run Draw a region, in headless Chromium (Playwright).
//
//   node tools/perf/bench_open.js --data-root <export root> --ds <id> [--label <name>] [options]
//
// It serves the page (web/SpatialViewer, or --page-root) and the export (--data-root) from two local
// HTTP servers in this process, the way blob storage answers: Range requests (206), ETag and
// Last-Modified, open CORS, and by default no Cache-Control header (Azure Blob's default; --cache-control
// sets one). Every response's bytes are counted on the server side too.
//
// Mode "open" (default): for each repetition, a fresh browser profile, then
//   cold   the first visit with an empty cache
//   warm   a second visit in the same browser (HTTP cache and service worker from the first visit)
// and for each visit records, in the page's clock (ms from navigation start):
//   first_tile_ms     the first frame in which an image tile is on screen
//   view_ms           the whole first view drawn: every image tile on show loaded (as the tests' tilesDone,
//                     4 calm checks 150 ms apart; the time recorded is the first of them) and no
//                     "Loading/Opening" status
//   bytes_open        bytes received until the network went quiet after the view (no request in flight for
//                     2 s): everything opening the dataset fetched, before any gene is picked
//   bytes_view        bytes the page received (Chrome's encodedDataLength over the network, from the page and
//                     its workers, CDN included; cache hits count 0) from navigation to view_ms
//   gene_ms           typing a gene and pressing Enter (dispatched in the page) to the first frame with that
//                     gene drawn; --genes picks which: top (most counts), median, or names
//   bytes_gene        bytes received for that gene
//   js_heap_mb        the page's used JS heap after the view settled (CDP Runtime.getHeapUsage)
//   browser_private_mb / browser_wset_mb   Windows: private bytes / working set summed over every process of
//                     the browser (browser, GPU, renderer, utility), 2 s after the view settled
// --net home applies Chrome DevTools network emulation to the page: 50 Mbit/s down, 10 Mbit/s up, 20 ms
// latency (--down/--up/--latency to change). It throttles what the page and its dedicated workers fetch.
//
// Mode "region": opens the dataset, switches to --layer, waits until Draw a region is ready, then hands the
// circles of --regions (JSON: {regions:[{target, points:[[x,y]...]}]}, image pixels) one at a time to
// window.__spatialViewer.lassoLoad, which starts the same run as the "Run UMAP" button, and times
//   region_ms   "Fetching N data chunk(s)" shown (the run started) to the UMAP panel drawn (two frames after)
// with Leiden off and/or on (--leiden off,on).
//
// Output: --out <dir>/<label>_<mode>_<stamp>.json (everything) and rows appended to --csv in the long format
//   date,category,dataset,platform,step,condition,metric,value,unit,rep,cpu,notes
//
// Needs Playwright: NODE_PATH=<folder holding node_modules/playwright>, PLAYWRIGHT_BROWSERS_PATH as installed.
// --gpu (default on) asks Chromium for the real GPU (ANGLE D3D11 on Windows); --no-gpu uses SwiftShader.
// CPU: start it under the affinity you want (Windows: start /affinity FF0000 ...); Chromium inherits it.
const http = require("http");
const fs = require("fs");
const path = require("path");
const crypto = require("crypto");
const { execFileSync } = require("child_process");
const { chromium } = require("playwright");

const argv = process.argv.slice(2);
const arg = (k, d) => { const i = argv.indexOf(k); return i >= 0 ? argv[i + 1] : d; };
const flag = (k) => argv.includes(k);
const REPO = path.resolve(__dirname, "..", "..");
const PAGE_ROOT = path.resolve(arg("--page-root", path.join(REPO, "web", "SpatialViewer")));
const DATA_ROOT = arg("--data-root") && path.resolve(arg("--data-root"));
const DS = arg("--ds");
const LABEL = arg("--label", DS);
const PLATFORM = arg("--platform", "");
const MODE = arg("--mode", "open");
const REPS = +arg("--reps", 3);
const NET = arg("--net", "none");
const DOWN = +arg("--down", 50), UP = +arg("--up", 10), LAT = +arg("--latency", 20);
const VISITS = arg("--visits", "cold,warm").split(",");
const GENES = arg("--genes", "top").split(",").filter(Boolean);
const CACHE_CONTROL = arg("--cache-control", "");
const OUT = path.resolve(arg("--out", "."));
const CSV = arg("--csv");
const COND_EXTRA = arg("--cond", "");
const NOTES = arg("--notes", "");
const GPU = !flag("--no-gpu");
const VIEW_W = +arg("--width", 1600), VIEW_H = +arg("--height", 1000);
const LAYER = arg("--layer");
const REGIONS = arg("--regions");
const TARGETS = (arg("--targets", "") || "").split(",").filter(Boolean).map(Number);
const LEIDEN = arg("--leiden", "off").split(",");
const TIMEOUT = +arg("--timeout-s", 180) * 1000;
if (!DATA_ROOT || !DS) { console.error("need --data-root and --ds"); process.exit(2); }
fs.mkdirSync(OUT, { recursive: true });
const sleep = (ms) => new Promise((r) => setTimeout(r, ms));
const log = (...a) => console.error(`[${new Date().toISOString().slice(11, 19)}] ${LABEL}:`, ...a);

// ---------------------------------------------------------------- servers (blob-storage-like)
const TYPES = { ".html": "text/html", ".js": "text/javascript", ".css": "text/css", ".json": "application/json",
  ".png": "image/png", ".jpg": "image/jpeg", ".svg": "image/svg+xml", ".webp": "image/webp", ".ico": "image/x-icon" };
function makeServer(root) {
  const st = { bytes: 0, requests: 0, s206: 0, s304: 0 };
  const srv = http.createServer((req, res) => {
    const cors = { "Access-Control-Allow-Origin": "*", "Access-Control-Expose-Headers": "Content-Range, Content-Length, Accept-Ranges, ETag" };
    if (req.method === "OPTIONS") {
      res.writeHead(204, { ...cors, "Access-Control-Allow-Methods": "GET, HEAD, OPTIONS",
        "Access-Control-Allow-Headers": req.headers["access-control-request-headers"] || "Range", "Access-Control-Max-Age": "3600" });
      return res.end();
    }
    st.requests++;
    let rel; try { rel = decodeURIComponent(new URL(req.url, "http://x").pathname); } catch { res.writeHead(400); return res.end(); }
    const file = path.join(root, rel);
    if (!file.startsWith(root)) { res.writeHead(403, cors); return res.end(); }
    let s; try { s = fs.statSync(file); if (s.isDirectory()) throw 0; } catch { res.writeHead(404, cors); return res.end(); }
    const etag = `"${s.size.toString(16)}-${Math.floor(s.mtimeMs).toString(16)}"`;
    const head = { ...cors, "Content-Type": TYPES[path.extname(file).toLowerCase()] || "application/octet-stream",
      "Accept-Ranges": "bytes", ETag: etag, "Last-Modified": new Date(s.mtimeMs).toUTCString(),
      ...(CACHE_CONTROL ? { "Cache-Control": CACHE_CONTROL } : {}) };
    if (req.headers["if-none-match"] === etag) { st.s304++; res.writeHead(304, head); return res.end(); }
    const m = /^bytes=(\d*)-(\d*)$/.exec((req.headers.range || "").trim());
    let start = 0, end = s.size - 1, code = 200;
    if (m && (m[1] || m[2])) {
      if (m[1]) { start = +m[1]; end = Math.min(m[2] ? +m[2] : s.size - 1, s.size - 1); } else { start = Math.max(0, s.size - +m[2]); }
      if (start >= s.size || start > end) { res.writeHead(416, { ...head, "Content-Range": `bytes */${s.size}` }); return res.end(); }
      code = 206; head["Content-Range"] = `bytes ${start}-${end}/${s.size}`; st.s206++;
    }
    head["Content-Length"] = String(end - start + 1);
    res.writeHead(code, head);
    if (req.method === "HEAD") return res.end();
    st.bytes += end - start + 1;
    fs.createReadStream(file, { start, end }).pipe(res);
  });
  return new Promise((r) => srv.listen(0, "127.0.0.1", () => r({ srv, st, port: srv.address().port })));
}

// ---------------------------------------------------------------- browser processes' memory (Windows)
function browserMem(profile) {
  if (process.platform !== "win32") return null;
  const ps = `$all = Get-CimInstance Win32_Process | Select-Object ProcessId,ParentProcessId,CommandLine;
    $root = $all | Where-Object { $_.CommandLine -and $_.CommandLine.Contains('${profile.replace(/'/g, "''")}') -and -not $_.CommandLine.Contains('--type=') } | Select-Object -First 1;
    if (-not $root) { '{}'; exit }
    $ids = New-Object System.Collections.Generic.HashSet[int]; [void]$ids.Add([int]$root.ProcessId);
    do { $n = $ids.Count; foreach ($p in $all) { if ($ids.Contains([int]$p.ParentProcessId)) { [void]$ids.Add([int]$p.ProcessId) } } } while ($ids.Count -ne $n)
    $pr = Get-Process -Id @($ids) -ErrorAction SilentlyContinue;
    @{ n = @($pr).Count; private_mb = [math]::Round((($pr | Measure-Object PrivateMemorySize64 -Sum).Sum)/1MB,1);
       wset_mb = [math]::Round((($pr | Measure-Object WorkingSet64 -Sum).Sum)/1MB,1) } | ConvertTo-Json -Compress`;
  try { return JSON.parse(execFileSync("powershell.exe", ["-NoProfile", "-Command", ps], { encoding: "utf8", timeout: 60000 }).trim() || "{}"); }
  catch (e) { return { error: String(e).slice(0, 200) }; }
}
function cpuName() {
  try { return execFileSync("powershell.exe", ["-NoProfile", "-Command", "(Get-CimInstance Win32_Processor).Name"], { encoding: "utf8" }).trim(); }
  catch { return process.arch; }
}
function affinityText() {
  if (process.platform !== "win32") return "";
  try {
    const m = +execFileSync("powershell.exe", ["-NoProfile", "-Command", `(Get-Process -Id ${process.pid}).ProcessorAffinity`], { encoding: "utf8" }).trim();
    const ids = []; for (let i = 0; i < 64; i++) if (Math.floor(m / 2 ** i) % 2) ids.push(i);
    return `${ids.length} of ${require("os").cpus().length} logical CPUs (${ids[0]}-${ids[ids.length - 1]})`;
  } catch { return ""; }
}

// ---------------------------------------------------------------- in-page probes
// installed before the page's own scripts: watches every animation frame for the first image tile drawn
const INIT = () => {
  const P = (window.__perf = { firstTile: null, frames: 0 });
  const tick = () => {
    P.frames++;
    try {
      const sv = window.__spatialViewer, lm = sv && sv.deck && sv.deck.layerManager;
      if (P.firstTile == null && lm) {
        // a sublayer of a TileLayer that carries a decoded tile image (renderSubLayers only builds one then)
        if (lm.getLayers().some((l) => l.props && l.props.image && l.parent && l.parent.state && l.parent.state.tileset)) P.firstTile = performance.now();
      }
    } catch (e) { P.err = String(e); }
    if (P.firstTile == null) requestAnimationFrame(tick);
  };
  requestAnimationFrame(tick);
};
// resolves when the first view is complete; returns the page time of the first of 4 calm checks
const VIEW_DONE = async (ms) => {
  const t0 = performance.now();
  let calm = 0, first = null;
  while (performance.now() - t0 < ms) {
    await new Promise((r) => setTimeout(r, 150));
    const sv = window.__spatialViewer, S = sv && sv.state, st = document.getElementById("status");
    let ok = !!(S && S.m && sv.deck && sv.deck.layerManager) && (!st || st.hidden || !/Loading|Opening/.test(st.textContent));
    if (ok && S.m.expr && S.m.expr.length) ok = !!(S.exprKey && S.geneIndex && S.geneIndex.length);
    if (ok) {
      const tl = sv.deck.layerManager.getLayers().filter((l) => l.state && l.state.tileset);
      ok = tl.every((l) => l.isLoaded && (l.state.tileset.selectedTiles || []).every((t) => t.isLoaded));
      if (S.m.image && !tl.length) ok = false;
    }
    if (ok) { if (calm === 0) first = performance.now(); calm++; } else calm = 0;
    if (calm >= 4) return { view: first, firstTile: window.__perf ? window.__perf.firstTile : null, err: st && st.classList.contains("err") ? st.textContent : null };
  }
  // what was still missing when the time ran out
  const sv = window.__spatialViewer, S = sv && sv.state, st = document.getElementById("status");
  const tl = sv && sv.deck && sv.deck.layerManager ? sv.deck.layerManager.getLayers().filter((l) => l.state && l.state.tileset) : [];
  return { view: null, firstTile: window.__perf ? window.__perf.firstTile : null, timeout: true,
           status: st ? `${st.hidden ? "hidden" : "shown"}: ${st.textContent}` : null,
           missing: { manifest: !!(S && S.m), exprKey: S ? S.exprKey || null : null, geneIndex: S && S.geneIndex ? S.geneIndex.length : null,
                      tileLayers: tl.length, tileLayersLoaded: tl.filter((l) => l.isLoaded).length,
                      tilesSelected: tl.reduce((a, l) => a + (l.state.tileset.selectedTiles || []).length, 0),
                      tilesSelectedLoaded: tl.reduce((a, l) => a + (l.state.tileset.selectedTiles || []).filter((t) => t.isLoaded).length, 0) } };
};
const PICK_GENES = (which) => {
  const S = window.__spatialViewer.state, idx = [...(S.geneIndex || [])].filter((g) => g && g.n && g.t > 0).sort((a, b) => b.t - a.t);
  return which.map((w) => {
    if (w === "top") return idx[0] ? { want: w, name: idx[0].n, total: idx[0].t } : null;
    if (w === "median") { const g = idx[Math.floor(idx.length / 2)]; return g ? { want: w, name: g.n, total: g.t } : null; }
    const g = idx.find((x) => x.n === w); return g ? { want: w, name: g.n, total: g.t } : { want: w, name: w };
  }).filter(Boolean);
};
const ADD_GENE = async (name, ms) => {
  const S = window.__spatialViewer.state, el = document.getElementById("geneSearch");
  const frame = () => new Promise((r) => requestAnimationFrame(() => r()));
  const t0 = performance.now();
  el.focus(); el.value = name; el.dispatchEvent(new Event("input", { bubbles: true }));
  el.dispatchEvent(new KeyboardEvent("keydown", { key: "Enter", code: "Enter", bubbles: true }));
  while (performance.now() - t0 < ms) {
    const g = S.genes.find((x) => x.name === name);
    if (g) { await frame(); await frame(); return { ms: performance.now() - t0, n: g.n, nGenes: S.genes.length }; }
    const st = document.getElementById("status");
    if (st && !st.hidden && st.classList.contains("err")) return { ms: null, err: st.textContent };
    await new Promise((r) => setTimeout(r, 10));
  }
  return { ms: null, err: "timeout" };
};

// ---------------------------------------------------------------- one visit
async function visit(ctx, url, kind, profile, servers, rec) {
  const page = await ctx.newPage();
  const cdp = await ctx.newCDPSession(page);
  const net = { bytes: 0, n: 0, cached: 0, fromSW: 0, byHost: {}, inflight: new Set(), last: Date.now() };
  const reqHost = new Map();
  cdp.on("Network.requestWillBeSent", (e) => { net.inflight.add(e.requestId); net.last = Date.now(); try { reqHost.set(e.requestId, new URL(e.request.url).host); } catch { /* data: */ } });
  cdp.on("Network.responseReceived", (e) => { if (e.response.fromDiskCache || e.response.fromPrefetchCache) net.cached++; if (e.response.fromServiceWorker) net.fromSW++; });
  cdp.on("Network.loadingFailed", (e) => { net.inflight.delete(e.requestId); net.last = Date.now(); });
  cdp.on("Network.loadingFinished", (e) => {
    net.inflight.delete(e.requestId); net.last = Date.now();
    net.bytes += e.encodedDataLength; net.n++;
    const h = reqHost.get(e.requestId) || "?"; net.byHost[h] = (net.byHost[h] || 0) + e.encodedDataLength;
  });
  await cdp.send("Network.enable");
  if (NET !== "none") await cdp.send("Network.emulateNetworkConditions", { offline: false, latency: LAT,
    downloadThroughput: (DOWN * 1e6) / 8, uploadThroughput: (UP * 1e6) / 8 });
  const errors = [];
  page.on("pageerror", (e) => errors.push(String(e).slice(0, 300)));
  page.on("dialog", (d) => d.dismiss().catch(() => {}));
  const srv0 = servers.map((s) => s.st.bytes);
  const out = { visit: kind, url };
  await page.goto(url, { waitUntil: "commit", timeout: TIMEOUT });
  const v = await page.evaluate(VIEW_DONE, TIMEOUT).catch((e) => ({ view: null, err: String(e) }));
  out.first_tile_ms = v.firstTile == null ? null : Math.round(v.firstTile);
  out.view_ms = v.view == null ? null : Math.round(v.view);
  out.bytes_view = net.bytes; out.requests_view = net.n; out.cache_hits_view = net.cached; out.sw_view = net.fromSW;
  out.bytes_by_host_view = { ...net.byHost };
  out.server_bytes_view = servers.map((s, i) => s.st.bytes - srv0[i]);
  out.view_err = v.err || (v.timeout ? `timeout after ${TIMEOUT / 1000} s (status ${v.status}; ${JSON.stringify(v.missing)})` : null);
  if (out.view_ms == null) {
    out.genes = []; out.errors = errors;
    await page.screenshot({ path: path.join(OUT, `${LABEL}_${kind}_viewtimeout_${Date.now()}.png`) }).catch(() => {});
    await page.close(); return out;
  }
  // everything the open started (positions, QC values, UMAP ...) has arrived: no request in flight for 2 s
  const tq = Date.now();
  while (Date.now() - tq < 120000 && (net.inflight.size || Date.now() - net.last < 2000)) await sleep(100);
  out.bytes_open = net.bytes; out.requests_open = net.n;
  out.idle_ms = await page.evaluate(() => Math.round(performance.now())) - 2000;
  await sleep(2000);
  const h = await cdp.send("Runtime.getHeapUsage").catch(() => null);
  out.js_heap_mb = h ? +(h.usedSize / 1e6).toFixed(1) : null;
  const mem = browserMem(profile);
  out.browser_private_mb = mem && mem.private_mb; out.browser_wset_mb = mem && mem.wset_mb; out.browser_procs = mem && mem.n;
  out.info = await page.evaluate(() => { const S = window.__spatialViewer.state; return { exprKey: S.exprKey || null, layers: (S.m.expr || []).map((e) => e.key),
    image: S.m.image ? [S.m.image.width, S.m.image.height, S.m.image.levels] : null, nGenes: (S.geneIndex || []).length, frames: window.__perf.frames }; });
  // where the time went: when the page's HTML, deck.gl and the dataset's manifest had arrived (Resource Timing)
  out.milestones = await page.evaluate(() => {
    const n = performance.getEntriesByType("navigation")[0], res = performance.getEntriesByType("resource");
    const end = (re) => { const r = res.find((x) => re.test(x.name)); return r ? Math.round(r.responseEnd) : null; };
    return { html_ms: n ? Math.round(n.responseEnd) : null, dcl_ms: n ? Math.round(n.domContentLoadedEventEnd) : null,
             deckgl_ms: end(/deck\.gl/), manifest_ms: end(/manifest\.json/), first_tile_fetched_ms: end(/\/tiles\//) };
  });
  out.genes = [];
  if (out.info.nGenes && GENES.length) {
    const picks = rec.genes || (rec.genes = await page.evaluate(PICK_GENES, GENES));
    for (const g of picks) {
      const b0 = net.bytes;
      const r = await page.evaluate(([n, ms]) => window.__ADD_GENE(n, ms), [g.name, TIMEOUT]).catch((e) => ({ ms: null, err: String(e) }));
      out.genes.push({ ...g, gene_ms: r.ms == null ? null : Math.round(r.ms), points: r.n, bytes_gene: net.bytes - b0, err: r.err || null });
      await sleep(500);
    }
  }
  out.errors = errors;
  await page.close();
  return out;
}

async function launch(profile) {
  const args = GPU ? ["--enable-gpu", "--use-angle=d3d11", "--ignore-gpu-blocklist"] : [];
  const ctx = await chromium.launchPersistentContext(profile, { headless: true, args, viewport: { width: VIEW_W, height: VIEW_H } });
  await ctx.addInitScript(INIT);
  await ctx.addInitScript(`window.__ADD_GENE = ${ADD_GENE.toString()};`);
  return ctx;
}

// ---------------------------------------------------------------- region (Draw a region -> Run UMAP)
async function regionRuns(ctx, url, profile, rec) {
  const page = await ctx.newPage();
  page.on("dialog", (d) => d.dismiss().catch(() => {}));
  await page.goto(url, { waitUntil: "commit", timeout: TIMEOUT });
  const v = await page.evaluate(VIEW_DONE, TIMEOUT);
  if (v.view == null) throw new Error("dataset did not open: " + JSON.stringify(v));
  if (LAYER) {
    await page.evaluate((k) => { const s = document.getElementById("exprSelect"); if (s.value !== k) { s.value = k; s.dispatchEvent(new Event("change")); } }, LAYER);
  }
  await page.waitForFunction((k) => { const S = window.__spatialViewer.state, L = S.lasso;
    return (!k || S.exprKey === k) && L && L.pos && L.pos.length > 0 && L.hvg && !document.getElementById("panelLasso").hidden; }, LAYER, { timeout: TIMEOUT, polling: 250 });
  await sleep(4000);
  const regions = JSON.parse(fs.readFileSync(REGIONS, "utf8")).regions.filter((r) => !TARGETS.length || TARGETS.includes(r.target));
  const out = [];
  for (const leiden of LEIDEN) {
    for (const region of regions) {
      const job = { tool: "SODA circle -> UMAP", version: 1, dataset: DS, layer: LAYER, seed: 42, mode: "bins",
        resolution: leiden === "on" ? "auto" : 0, shapes: [{ name: `circle ${region.target}`, color: "#e6194b", kind: "circle", points: region.points }], qcFilter: null };
      const peak = { private_mb: 0 };
      let sampling = true;
      const sampler = (async () => { while (sampling) { const m = browserMem(profile); if (m && m.private_mb) peak.private_mb = Math.max(peak.private_mb, m.private_mb); await sleep(1000); } })();
      const r = await page.evaluate(async ([job, ms]) => {
        const V = window.__spatialViewer, S = V.state, msgEl = document.getElementById("lassoMsg"), panel = document.getElementById("lassoPanel");
        const msgs = [];
        const mo = new MutationObserver(() => msgs.push([performance.now(), msgEl.textContent.slice(0, 160)]));
        mo.observe(msgEl, { childList: true, characterData: true, subtree: true });
        const frame = () => new Promise((r) => requestAnimationFrame(() => r()));
        const prev = S.lasso.result;
        const t_load = performance.now();
        await V.lassoLoad(job);
        const t_loaded = performance.now();
        let err = null;
        while (true) {
          const L = S.lasso;
          if (L.result && L.result !== prev && !panel.hidden) break;
          if (!L.running && msgEl.style.color && performance.now() - t_loaded > 500) { err = msgEl.textContent; break; }
          if (performance.now() - t_load > ms) { err = "timeout"; break; }
          await new Promise((r) => setTimeout(r, 50));
        }
        await frame(); await frame();
        const t_shown = performance.now();
        mo.disconnect();
        const R = S.lasso.result, click = msgs.find((m) => /^Fetching/.test(m[1]));
        return { err, t_click: click ? click[0] : null, t_load, t_shown, secs: R ? R.secs : null,
                 nIn: R ? R.info.nIn : null, nClusters: R ? R.info.nClusters || null : null, bins: S.lasso.shapes.map((s) => s.bins)[0],
                 heap_mb: performance.memory ? performance.memory.usedJSHeapSize / 1e6 : null };
      }, [job, TIMEOUT]);
      sampling = false; await sampler;
      const row = { target: region.target, leiden, ...r, region_ms: r.t_click == null ? null : Math.round(r.t_shown - r.t_click),
                    load_to_shown_ms: Math.round(r.t_shown - r.t_load), browser_private_peak_mb: peak.private_mb };
      log(`region ${region.target} leiden ${leiden}: ${row.region_ms} ms (nIn ${row.nIn}, app ${r.secs} s)`);
      out.push(row);
      await sleep(1500);
    }
  }
  await page.close();
  return out;
}

// ---------------------------------------------------------------- main
(async () => {
  const pageSrv = await makeServer(PAGE_ROOT), dataSrv = await makeServer(DATA_ROOT);
  const url = `http://127.0.0.1:${pageSrv.port}/index.html?static=http://127.0.0.1:${dataSrv.port}&ds=${encodeURIComponent(DS)}`;
  const cpu = `${cpuName()}; ${affinityText()}`;
  const netCond = NET === "none" ? "no throttling (localhost)" : `CDP ${DOWN}/${UP} Mbit/s, ${LAT} ms`;
  const cond = [netCond, GPU ? "GPU" : "SwiftShader", COND_EXTRA].filter(Boolean).join("; ");
  const stamp = new Date().toISOString().replace(/[:.]/g, "-").slice(0, 19);
  const res = { label: LABEL, ds: DS, mode: MODE, url: url.replace(/127\.0\.0\.1:\d+/g, "127.0.0.1:<port>"), dataRoot: DATA_ROOT, pageRoot: PAGE_ROOT,
    net: NET === "none" ? null : { down_mbit: DOWN, up_mbit: UP, latency_ms: LAT }, cacheControl: CACHE_CONTROL || null,
    gpu: GPU, viewport: [VIEW_W, VIEW_H], cpu, started: new Date().toISOString(), runs: [] };
  const rec = {};
  const rows = [];
  const add = (step, visit, metric, value, unit, rep, extra = "") => rows.push([res.started, MODE === "region" ? "region" : "open", LABEL, PLATFORM, step,
    [visit, cond].filter(Boolean).join("; "), metric, value, unit, rep, cpu, [NOTES, extra].filter(Boolean).join("; ")]);
  for (let rep = 1; rep <= REPS; rep++) {
    const profile = path.join(OUT, "profiles", `${LABEL}_${MODE}_${stamp}_r${rep}_${crypto.randomBytes(3).toString("hex")}`);
    const ctx = await launch(profile);
    res.browser = ctx.browser() ? ctx.browser().version() : null;
    try {
      if (MODE === "region") {
        const runs = await regionRuns(ctx, url, profile, rec);
        res.runs.push({ rep, runs });
        for (const r of runs) {
          const st = `region ${r.target} bins, Leiden ${r.leiden}`;
          add(st, "", "region_ms", r.region_ms, "ms", rep, `nIn ${r.nIn}; bins in shape ${r.bins}; app-reported ${r.secs} s` + (r.err ? `; ERROR ${r.err}` : ""));
          add(st, "", "browser_private_peak_mb", r.browser_private_peak_mb, "MB", rep);
        }
      } else {
        for (const kind of VISITS) {
          const v = await visit(ctx, url, kind, profile, [pageSrv, dataSrv], rec);
          log(`${kind} r${rep}: first tile ${v.first_tile_ms} ms, view ${v.view_ms} ms, ${(v.bytes_view / 1e6).toFixed(2)} MB (all of the open ${((v.bytes_open || 0) / 1e6).toFixed(2)} MB)` +
              (v.genes || []).map((g) => `; gene ${g.want} ${g.gene_ms} ms ${(g.bytes_gene / 1e6).toFixed(2)} MB`).join("") +
              `; heap ${v.js_heap_mb} MB, browser ${v.browser_private_mb} MB` + (v.view_err ? `; ${v.view_err}` : "") +
              (v.milestones ? `; ${JSON.stringify(v.milestones)}` : ""));
          res.runs.push({ rep, ...v });
          add("open", kind, "first_tile_ms", v.first_tile_ms, "ms", rep, v.view_err || "");
          add("open", kind, "view_ms", v.view_ms, "ms", rep);
          for (const k of ["deckgl_ms", "manifest_ms"]) if (v.milestones && v.milestones[k] != null) add("open", kind, k, v.milestones[k], "ms", rep);
          add("open", kind, "bytes_view", v.bytes_view, "B", rep, `${v.requests_view} responses; ${v.cache_hits_view} from disk cache`);
          if (v.bytes_open != null) add("open", kind, "bytes_open", v.bytes_open, "B", rep, `${v.requests_open} responses; network quiet at ${v.idle_ms} ms`);
          if (v.js_heap_mb != null) add("open", kind, "js_heap_mb", v.js_heap_mb, "MB", rep);
          if (v.browser_private_mb != null) add("open", kind, "browser_private_mb", v.browser_private_mb, "MB", rep, `${v.browser_procs} processes`);
          if (v.browser_wset_mb != null) add("open", kind, "browser_wset_mb", v.browser_wset_mb, "MB", rep);
          for (const g of v.genes || []) {
            add(`gene (${g.want})`, kind, "gene_ms", g.gene_ms, "ms", rep, `${g.points} points` + (g.err ? `; ERROR ${g.err}` : ""));
            add(`gene (${g.want})`, kind, "bytes_gene", g.bytes_gene, "B", rep);
          }
        }
      }
    } catch (e) { log("FAILED", e); res.runs.push({ rep, error: String(e).slice(0, 500) }); add(MODE, "", "error", 1, "", rep, String(e).slice(0, 200)); }
    await ctx.close().catch(() => {});
    await sleep(1000);
  }
  res.genes = rec.genes || null;
  res.server = { page: pageSrv.st, data: dataSrv.st };
  res.finished = new Date().toISOString();
  const file = path.join(OUT, `${LABEL}_${MODE}${NET === "none" ? "" : "_" + NET}_${stamp}.json`);
  fs.writeFileSync(file, JSON.stringify(res, null, 1));
  if (CSV) {
    const q = (x) => (x == null ? "" : /[",\n]/.test(String(x)) ? `"${String(x).replace(/"/g, '""')}"` : String(x));
    if (!fs.existsSync(CSV)) fs.writeFileSync(CSV, "date,category,dataset,platform,step,condition,metric,value,unit,rep,cpu,notes\n");
    fs.appendFileSync(CSV, rows.map((r) => r.map(q).join(",")).join("\n") + "\n");
  }
  log("wrote", file);
  pageSrv.srv.close(); dataSrv.srv.close();
  process.exit(0);
})().catch((e) => { console.error(e); process.exit(1); });
