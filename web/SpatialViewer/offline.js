// Offline use: the service worker (sw.js) and the "Offline use" section of the Save tab.
//
// After one visit online the page opens and works with no network: sw.js keeps index.html, every file
// it names with a ?v= string, the workers, the marker lists, deck.gl and umap-js. That needs a page
// served over https or from localhost; a page opened as a file (file://) cannot register a service
// worker, so it is not kept (the section says so).
//
// Data is NOT kept on its own: an export can be gigabytes, and a private site's data sits behind a
// password. A dataset is kept only when the reader presses "Keep this dataset for offline use", after
// seeing the size: what the page has fetched from that dataset so far (the image tiles looked at, genes,
// clusters, the UMAP) plus the files Draw a region needs (each layer's hvg/ or mat/ chunks), up to
// offlineMaxMB (config.js, default 500 MB). It can be updated (adds what was looked at since) and
// removed, one dataset or all. Never kept: a dataset of a password-protected site (tokenEndpoint), one
// from an encrypted link file (linkPath), one whose link carries a read pass (?data= with a query),
// one from a SODA server (?api=, server.py). Nothing decrypted with a password is written to disk.
//
// config.js:  offline: false   turns all of it off (and removes what this page kept)
//             offline: true    turns it on even on a password-protected site (the page only; data
//                              of such a site is still never kept)
//             offlineMaxMB: n  the most one kept dataset may hold (default 500)
// Default: on, except on a password-protected site (tokenEndpoint), where nothing opens without the
// password service anyway.
//
// Caches (named after this page's folder, so another copy of SODA on the same site keeps its own):
// soda-app:<folder>|<hash of the ?v= strings>, soda-lib:<folder>|<libraries>, soda-data:<folder>|<dataset
// folder>. To clear them all: offline: false in config.js, or the browser's "clear site data".
/* exported SodaOffline */
(function () {
  "use strict";
  const cfg = window.SpatialViewerConfig || {};
  const $ = (id) => document.getElementById(id);
  const esc = (x) => String(x).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
  const SCOPE = new URL("./", location.href).href;
  const PFX = { app: `soda-app:${SCOPE}|`, lib: `soda-lib:${SCOPE}|`, data: `soda-data:${SCOPE}|` };
  const META = "?soda-offline-meta";
  const MAX_MB = Number(cfg.offlineMaxMB) > 0 ? Number(cfg.offlineMaxMB) : 500;
  const WANT = cfg.offline !== undefined ? !!cfg.offline : !cfg.tokenEndpoint;
  const CAN = typeof navigator !== "undefined" && "serviceWorker" in navigator && window.isSecureContext &&
              /^https?:$/.test(location.protocol) && typeof caches !== "undefined";
  // state: "starting" | "registering" | "ready" | "failed" | "off" | "unsupported"
  const O = { state: "starting", why: "", app: "", busy: null, ds: null, kept: [] };
  let api = null, readyResolve = null;
  const ready = new Promise((r) => { readyResolve = r; });
  const sleep = (ms) => new Promise((r) => setTimeout(r, ms));
  const bare = (u) => String(u).split("#")[0].split("?")[0];
  const mb = (b) => (b >= 1e9 ? `${(b / 1e9).toFixed(2)} GB` : b >= 1e6 ? `${(b / 1e6).toFixed(1)} MB` : `${Math.max(1, Math.round(b / 1e3))} kB`);

  // ---- the files this page fetched (address without its query, size when the browser gives it) ----------
  // A PerformanceObserver sees every request, also past the timeline's 250-entry buffer. Kept in memory
  // only, for "Keep this dataset": a part after "?" (a read pass, a token) is never kept.
  const seen = new Map();
  function note(e) {
    const u = bare(e.name);
    if (!/^https?:/.test(u) || seen.size > 200000) return;
    const b = e.encodedBodySize || 0;
    if (!seen.has(u) || b) seen.set(u, b);
  }
  try { new PerformanceObserver((l) => l.getEntries().forEach(note)).observe({ type: "resource", buffered: true }); }
  catch (e) { try { performance.getEntriesByType("resource").forEach(note); } catch (err) { /* no timing */ } }

  // ---- the page's ?v= files and their hash: the service worker's version (the same rule as sw.js) ------
  function appRefs(list) {
    const out = new Set();
    for (const r of list) {
      if (!r || /^[a-z][\w+.-]*:|^\/\//i.test(r) || /^vendor\//.test(r) || !/[?&]v=/.test(r)) continue;
      out.add(r);
    }
    return [...out].sort();
  }
  const pageRefs = () => appRefs([...document.querySelectorAll("script[src],link[href]")].map((e) => e.getAttribute("src") || e.getAttribute("href")));
  const htmlRefs = (html) => appRefs([...String(html).matchAll(/<(?:script|link)\b[^>]*?\s(?:src|href)="([^"]+)"/gi)].map((m) => m[1]));
  function hashRefs(list) {                        // FNV-1a, 32 bits, as hex
    let h = 0x811c9dc5;
    const s = list.join("\n");
    for (let i = 0; i < s.length; i++) { h ^= s.charCodeAt(i); h = Math.imul(h, 0x01000193) >>> 0; }
    return ("0000000" + h.toString(16)).slice(-8);
  }

  // ---- the service worker ---------------------------------------------------------------------------------
  const ours = (r) => { const w = r && (r.active || r.waiting || r.installing); return !!w && r.scope === SCOPE && /\/sw\.js$/.test(new URL(w.scriptURL).pathname); };
  async function start() {
    if (!CAN) {
      O.state = "unsupported";
      O.why = location.protocol === "file:" ? "file" : !window.isSecureContext ? "insecure" : "browser";
      return done();
    }
    if (!WANT) { O.state = "off"; await turnOff(); return done(); }
    O.app = hashRefs(pageRefs());
    O.state = "registering"; paint();
    try {
      const reg = await navigator.serviceWorker.register(`sw.js?app=${O.app}`, { updateViaCache: "none" });
      // this version installed (its files kept) and activated (older versions' caches deleted)
      const mine = () => reg.active && reg.active.state === "activated" && reg.active.scriptURL.indexOf(`app=${O.app}`) > 0;
      // wait for this version to be installed (its files kept) and to take over: at most two minutes
      for (let i = 0; i < 480 && !(mine() && navigator.serviceWorker.controller); i++) {
        if (!reg.installing && !reg.waiting && !(reg.active && reg.active.state === "activating") && i > 40 && !mine()) break;   // nothing coming: the install failed
        await sleep(250);
      }
      O.state = mine() ? "ready" : "failed";
      if (!mine()) O.why = reg.active ? "an earlier version stays kept" : "the browser did not keep it";
    } catch (e) {
      // offline, some browsers (WebKit) fail register() because they cannot fetch sw.js to compare it, while
      // this same version, already kept, answers the page: that is ready, not failed
      let reg = null;
      try { reg = await navigator.serviceWorker.getRegistration(SCOPE); } catch (err) { /* none */ }
      if (reg && reg.active && reg.active.scriptURL.indexOf(`app=${O.app}`) > 0 && navigator.serviceWorker.controller) O.state = "ready";
      else { O.state = "failed"; O.why = e && e.message ? e.message.replace(/\?[^\s]*/g, "") : String(e); }
    }
    return done();
  }
  function done() { paint(); readyResolve(O.state); return O.state; }
  // offline: false in config.js: no service worker, and what an earlier visit kept goes
  async function turnOff() {
    try { for (const r of await navigator.serviceWorker.getRegistrations()) if (ours(r)) await r.unregister(); } catch (e) { /* none */ }
    try { for (const n of await caches.keys()) if (Object.values(PFX).some((p) => n.startsWith(p))) await caches.delete(n); } catch (e) { /* none */ }
  }
  async function tellWorker() {
    try {
      const reg = await navigator.serviceWorker.getRegistration(SCOPE);
      const msg = { type: "soda-kept-changed" };
      if (reg && reg.active) reg.active.postMessage(msg);
      if (navigator.serviceWorker.controller && (!reg || navigator.serviceWorker.controller !== reg.active)) navigator.serviceWorker.controller.postMessage(msg);
    } catch (e) { /* no worker: nothing reads the caches */ }
  }

  // ---- kept datasets ----------------------------------------------------------------------------------------
  async function readKept() {
    const out = [];
    if (typeof caches === "undefined") return out;
    for (const name of await caches.keys()) {
      if (!name.startsWith(PFX.data)) continue;
      const root = name.slice(PFX.data.length);
      let meta = null;
      try { const r = await (await caches.open(name)).match(root + META); meta = r ? await r.json() : null; } catch (e) { /* none */ }
      out.push({ name, root, meta: meta || {} });
    }
    O.kept = out;
    return out;
  }
  // what keeping the open dataset takes: what was fetched from it, plus Draw a region's chunks
  async function plan(ds) {
    const root = ds.root, m = ds.manifest || {}, files = new Map();
    const add = (u, b) => { u = bare(u); if (!files.has(u) || b) files.set(u, Number(b) || 0); };
    add(root + "manifest.json");
    for (const [u, b] of seen) if (u.startsWith(root)) add(u, b);
    for (const e of Array.isArray(m.expr) ? m.expr : []) {
      if (!e || typeof e.key !== "string") continue;
      const dirs = [];
      if (e.hvg) dirs.push("hvg");
      if (e.kind === "protein" && e.mat) dirs.push("mat");
      for (const dir of dirs) {
        const base = `${root}expr/${e.key}/${dir}/`;
        try {
          const r = await fetch(base + "index.json", { cache: "no-store" });
          if (!r.ok) continue;
          const idx = await r.json();
          add(base + "index.json");
          if (dir === "hvg") add(base + "genes.json");
          for (const c of Array.isArray(idx.chunks) ? idx.chunks : []) if (Number.isInteger(c.tx) && Number.isInteger(c.ty)) add(`${base}${c.tx}_${c.ty}.bin`, c.bytes);
        } catch (err) { /* this layer cannot be drawn on offline */ }
      }
      if (dirs.length) add(`${root}expr/${e.key}/positions.bin`);
    }
    // the folder's list of datasets and its description (a ?data= or ?static= folder), when it has them
    const extras = ds.folder ? ["datasets.json", "about.html"].map((f) => ds.folder + f) : [];
    extras.forEach((u) => add(u));
    // sizes the browser did not give (another site's files): ask with HEAD; a file that is not there is left out
    let unknown = 0;
    await pool([...files].filter(([, b]) => !b).map(([u]) => u), 6, async (u) => {
      if (O.busy && O.busy.stop) return;
      try {
        const r = await fetch(u, { method: "HEAD", cache: "no-store" });
        if (r.status === 404 || r.status === 410) { files.delete(u); return; }
        const n = Number(r.headers.get("Content-Length")) || 0;
        if (n) files.set(u, n); else unknown++;
      } catch (e) { unknown++; }
    });
    let bytes = 0;
    for (const b of files.values()) bytes += b;
    return { files, bytes, unknown, extras: extras.filter((u) => files.has(u)) };
  }
  async function pool(items, n, fn) {
    let next = 0;
    await Promise.all(Array.from({ length: Math.min(n, items.length) }, async () => { while (next < items.length) await fn(items[next++]); }));
  }
  async function keep(ds, p) {
    const name = PFX.data + ds.root, had = await caches.has(name), cache = await caches.open(name);
    const old = had ? ((await readKept()).find((k) => k.name === name) || {}).meta || {} : {};
    const sizes = old.sizes && typeof old.sizes === "object" ? { ...old.sizes } : {};
    const cap = MAX_MB * 1e6, urls = [...p.files.keys()];
    let got = 0, failed = 0, full = false;
    const total = () => Object.values(sizes).reduce((a, b) => a + b, 0);
    await pool(urls, 4, async (u) => {
      if (O.busy.stop || full) return;
      try {
        const r = await fetch(u, { cache: "no-store" });
        if (r.status !== 200) { failed++; return; }
        const buf = await r.arrayBuffer();
        if (total() - (sizes[u] || 0) + buf.byteLength > cap) { full = true; return; }
        await cache.put(u, new Response(buf, { headers: { "Content-Type": r.headers.get("Content-Type") || "application/octet-stream" } }));
        sizes[u] = buf.byteLength; got++;
        msg(`Keeping… ${got} of ${urls.length} files, ${mb(total())}`);
      } catch (e) {
        if (e && e.name === "QuotaExceededError") full = true; else failed++;
      }
    });
    if (O.busy.stop && !had) { await caches.delete(name); await tellWorker(); return { stopped: true }; }
    const meta = { v: 1, id: ds.id, label: ds.label, root: ds.root, extras: p.extras, files: Object.keys(sizes).length,
                   bytes: total(), failed, full, kept: new Date().toISOString(), sizes };
    await cache.put(ds.root + META, new Response(JSON.stringify(meta), { headers: { "Content-Type": "application/json" } }));
    await tellWorker();
    return { got, failed, full, bytes: meta.bytes, stopped: !!O.busy.stop };
  }
  async function remove(name) {
    await caches.delete(name);
    await tellWorker();
  }

  // ---- the section ------------------------------------------------------------------------------------------
  function msg(t, err) { const m = $("offMsg"); if (m) { m.textContent = t || ""; m.classList.toggle("err", !!err); } }
  function build() {
    if ($("panelOffline")) return;
    const sec = document.createElement("section");
    sec.className = "panel"; sec.id = "panelOffline"; sec.dataset.tab = "save";
    sec.innerHTML =
      `<h2>Offline use</h2>` +
      `<div id="offState" class="meta"></div>` +
      `<div id="offDsBox" hidden>` +
        `<div class="sub-h">This dataset</div>` +
        `<div id="offDs" class="meta"></div>` +
        `<div class="btnrow">` +
          `<button id="offKeep" type="button" title="Keep what you have looked at in this dataset, and what Draw a region needs, on this computer">Keep this dataset for offline use</button>` +
          `<button id="offGo" type="button" hidden>Keep it</button>` +
          `<button id="offStop" type="button" hidden>Cancel</button>` +
          `<button id="offDrop" type="button" hidden>Remove offline copy</button>` +
        `</div>` +
        `<div id="offMsg" class="meta" role="status"></div>` +
      `</div>` +
      `<div id="offListBox" hidden>` +
        `<div class="sub-h">Kept on this computer</div>` +
        `<div id="offList" class="list"></div>` +
        `<div class="btnrow"><button id="offDropAll" type="button">Remove every offline copy</button></div>` +
      `</div>`;
    (document.getElementById("sidebar") || document.body).appendChild(sec);
    if (window.SodaLayout) window.SodaLayout.adopt(sec);
    $("offKeep").onclick = askKeep;
    $("offGo").onclick = goKeep;
    $("offStop").onclick = () => { if (O.busy) O.busy.stop = "user"; if (O.busy && O.busy.stage === "confirm") { O.busy = null; msg(""); paint(); } };
    $("offDrop").onclick = async () => { const ds = O.ds; if (!ds || !ds.root) return; await remove(PFX.data + ds.root); msg("The offline copy of this dataset was removed."); paint(); };
    $("offDropAll").onclick = async () => {
      for (const k of await readKept()) await caches.delete(k.name);
      await tellWorker(); msg("Every offline copy on this computer was removed."); paint();
    };
    $("offList").onclick = async (e) => {
      const b = e.target.closest("button[data-name]");
      if (!b) return;
      await remove(b.dataset.name); paint();
    };
    window.addEventListener("online", () => paint());
    window.addEventListener("offline", () => paint());
    window.addEventListener("soda:opened", () => { if (!O.busy) msg(""); paint(); });
  }
  const NOT_KEPT = {
    password: "This dataset is not kept for offline use: it is password-protected, and nothing from a password-protected source is written to this computer.",
    link: "This dataset is not kept for offline use: it is opened through a password-protected link.",
    pass: "This dataset is not kept for offline use: its link carries a read pass, which is not written to this computer.",
    server: "This dataset is not kept for offline use: it comes from a SODA server (server.py), not from exported files.",
  };
  function stateText() {
    const offline = navigator.onLine === false;
    switch (O.state) {
      case "ready": return offline
        ? "You are offline. The page runs from the copy this browser kept; only datasets kept for offline use open."
        : "This page is kept by this browser (with deck.gl, umap-js and the marker lists), so it opens with no network after this visit. Data is kept only when you ask, one dataset at a time.";
      case "registering": return "Getting the page ready for offline use…";
      case "failed": return `This version of the page could not be kept for offline use (${O.why}). It works as usual while online.`;
      case "off": return cfg.offline === false ? "Offline use is turned off on this site." : "Offline use is off on this password-protected site: its data needs the password service.";
      case "unsupported": return O.why === "file"
        ? "Offline use needs the page served from a web address (https, or http://localhost). A page opened as a file cannot keep itself."
        : O.why === "insecure" ? "Offline use needs https (or http://localhost)." : "This browser cannot keep pages for offline use.";
      default: return "";
    }
  }
  async function paint() {
    if (!$("panelOffline")) return;
    $("offState").textContent = stateText();
    const ds = O.ds = api && api.dataset ? api.dataset() : null;
    const kept = O.state === "ready" || O.state === "failed" ? await readKept() : [];
    const box = $("offDsBox");
    box.hidden = !ds || !(O.state === "ready" || O.state === "failed");
    if (!box.hidden) {
      const mine = ds.root ? kept.find((k) => k.root === ds.root) : null;
      const busy = !!O.busy;
      let t;
      if (!ds.root) t = NOT_KEPT[ds.why] || "This dataset cannot be kept for offline use.";
      else if (mine) t = `Kept for offline use: ${mb(mine.meta.bytes || 0)} in ${mine.meta.files || 0} files, ${dateOf(mine.meta.kept)}${mine.meta.full ? " (stopped at the size limit)" : ""}. ` +
                         "Update adds what you have looked at since.";
      else t = `Not kept. Keeping it lets it open with no network: what you have looked at in it so far (image tiles, genes, clusters) ` +
               `and what Draw a region needs, up to ${MAX_MB} MB. Zoom and pan over what you will need first.`;
      $("offDs").textContent = t;
      $("offKeep").hidden = !ds.root || (busy && O.busy.stage !== "idle");
      $("offKeep").textContent = mine ? "Update offline copy" : "Keep this dataset for offline use";
      $("offKeep").disabled = O.state !== "ready" || navigator.onLine === false;
      $("offDrop").hidden = !mine || busy;
      $("offGo").hidden = !(busy && O.busy.stage === "confirm");
      $("offStop").hidden = !busy;
    }
    const lb = $("offListBox");
    lb.hidden = !kept.length;
    $("offList").innerHTML = kept.map((k) => `<div class="row slim off-row"><span>${esc(k.meta.label || k.meta.id || k.root)} · ${mb(k.meta.bytes || 0)}` +
      `${ds && ds.root === k.root ? " · open now" : ""}</span> <button type="button" data-name="${esc(k.name)}" title="Remove this offline copy">Remove</button></div>`).join("");
  }
  const dateOf = (iso) => { const d = new Date(iso); return isNaN(d) ? "" : `on ${d.toLocaleDateString(undefined, { year: "numeric", month: "short", day: "numeric" })}`; };

  async function askKeep() {
    const ds = O.ds;
    if (!ds || !ds.root || O.busy) return;
    O.busy = { stage: "plan", stop: null };
    paint();
    msg("Working out the size…");
    try {
      const p = await plan(ds);
      if (O.busy.stop) { O.busy = null; msg(""); paint(); return; }
      if (p.bytes > MAX_MB * 1e6) {
        O.busy = null; paint();
        msg(`That is about ${mb(p.bytes)}, more than the ${MAX_MB} MB this site keeps for one dataset (offlineMaxMB in config.js). ` +
            "Reload the page, look at less of it (fewer genes, less zooming in), and keep it then.", true);
        return;
      }
      let room = "";
      try {
        const est = navigator.storage && navigator.storage.estimate ? await navigator.storage.estimate() : null;
        if (est && est.quota && est.quota - (est.usage || 0) < p.bytes * 1.2) {
          O.busy = null; paint();
          msg(`This browser has room for about ${mb(Math.max(0, est.quota - (est.usage || 0)))} more, not enough for ${mb(p.bytes)}.`, true);
          return;
        }
        if (est && est.quota) room = ` (this browser has room for ${mb(est.quota - (est.usage || 0))})`;
      } catch (e) { /* no estimate */ }
      O.busy = { stage: "confirm", stop: null, ds, plan: p };
      paint();
      msg(`About ${mb(p.bytes)}${p.unknown ? ` plus ${p.unknown} files of unknown size` : ""} in ${p.files.size} files${room}. Keep them on this computer?`);
    } catch (e) {
      O.busy = null; paint();
      msg(`Could not work out the size: ${e.message}`, true);
    }
  }
  async function goKeep() {
    const b = O.busy;
    if (!b || b.stage !== "confirm") return;
    O.ds = api && api.dataset ? api.dataset() : null;
    // the same dataset is still open (dataset() returns a new object on each call: compare its folder)
    if (!b || b.stage !== "confirm" || !O.ds || b.ds.root !== O.ds.root) {
      O.busy = null; msg("Another dataset was opened: press Keep again."); paint(); return;
    }
    b.stage = "keep";
    paint();
    try {
      const r = await keep(b.ds, b.plan);
      O.busy = null;
      if (r.stopped) msg("Stopped: nothing more was kept.");
      else msg(`Kept: ${mb(r.bytes)}.${r.failed ? ` ${r.failed} files could not be fetched and are missing offline.` : ""}${r.full ? ` Stopped at the size limit (${MAX_MB} MB) or the browser's storage limit: part of it is missing offline.` : ""}`, r.failed > 0 || r.full);
    } catch (e) {
      O.busy = null;
      msg(`Could not keep it: ${e.message}`, true);
    }
    paint();
  }

  window.SodaOffline = {
    // viewer.js: what it knows about the open dataset. dataset() -> null, or { id, label, manifest, and either
    // root (the dataset's folder address, no query) and folder (the ?data= / ?static= folder above it, or null),
    // or why ("password" | "link" | "pass" | "server") it is not kept }
    init(a) { api = a; build(); paint(); },
    summary() {
      if (O.state === "ready") return navigator.onLine === false ? "offline" : O.kept.length ? `${O.kept.length} kept` : "ready";
      return O.state === "registering" ? "getting ready" : O.state === "failed" ? "not kept" : "off";
    },
    ready,                                         // resolves with the state once the service worker is settled
    state: () => ({ state: O.state, why: O.why, app: O.app, scope: SCOPE, busy: O.busy ? O.busy.stage : null, seen: seen.size,
                    kept: O.kept.map((k) => ({ root: k.root, id: k.meta.id, bytes: k.meta.bytes, files: k.meta.files })) }),
    htmlRefs, hashRefs, pageRefs,                  // tests: the same hash as sw.js
  };
  // after the page has loaded, so the service worker's own downloads never compete with the first picture
  const go = () => { build(); start(); };
  if (document.readyState === "complete") setTimeout(go, 0); else window.addEventListener("load", () => setTimeout(go, 0));
})();
