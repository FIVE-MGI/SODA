// SODA's service worker (registered by offline.js): after one visit online, the page opens and works with
// no network.
//
// What it keeps, in caches named after this page's folder (its scope), so another copy of SODA on the
// same site never touches these:
//   soda-app:<scope>|<app>    index.html, every file it names with a ?v= string, what those name in turn
//                             (the workers, leiden.js), and the marker lists. <app> is a hash of the page's
//                             ?v= strings (offline.js puts it on this script's address), so a new publish
//                             installs a new cache and the old one is deleted when it takes over.
//   soda-lib:<scope>|<libs>   deck.gl and umap-js (libs.js): the vendor/ copy when there, else the CDN's,
//                             each kept only when it matches its pinned SHA-384.
//   soda-data:<scope>|<root>  a dataset kept with "Keep this dataset for offline use" (offline.js). Only the
//                             page writes these, only on that click, and never for a password-protected
//                             dataset; this file only reads them, when the network fails.
// Anything else (every other dataset, tiles, the token service, other sites) is not touched: the browser
// fetches it as it would without a service worker.
//
// How it answers:
//   the page (a navigation to index.html)   the network first; the kept copy when the network fails
//   the app's files (?v=, marker lists)     the network first; the kept copy when the network fails
//   deck.gl, umap-js                        the kept copy first (pinned by version and hash, never change)
//   a kept dataset's files                  the network first; the kept copy when the network fails
//                                           (a Range request is cut out of the kept whole file)
// So while online the page behaves as before (a file edited without a new ?v= still arrives), and the
// kept copies are used only offline.
/* global SODA_LIBS */
"use strict";
importScripts("libs.js");

const SCOPE = self.registration.scope;
const APP = new URL(self.location.href).searchParams.get("app") || "0";
const APP_CACHE = `soda-app:${SCOPE}|${APP}`;
const LIB_CACHE = `soda-lib:${SCOPE}|${Object.keys(SODA_LIBS).map((k) => SODA_LIBS[k].path).join("+")}`;
const DATA_PREFIX = `soda-data:${SCOPE}|`;
const META = "?soda-offline-meta";                 // a kept dataset's own record, at <root>?soda-offline-meta
const PAGE = new URL("index.html", SCOPE).href;
const SCOPE_PATH = new URL(SCOPE).pathname;
// the marker lists in markers/ (tests/smoke_test.py checks this list against the folder)
const MARKERS = ["akoya_nerve.json", "brain.json", "nerve.json", "protein_generic.json", "rna_generic.json"];
const LIBS = Object.keys(SODA_LIBS).map((k) => ({ ...SODA_LIBS[k], localUrl: new URL(SODA_LIBS[k].local, SCOPE).href }));

// ---- the page's ?v= files (the same rule and hash as offline.js) ------------------------------------
function appRefs(html) {
  const out = new Set();
  for (const m of String(html).matchAll(/<(?:script|link)\b[^>]*?\s(?:src|href)="([^"]+)"/gi)) {
    const r = m[1];
    if (/^[a-z][\w+.-]*:|^\/\//i.test(r) || /^vendor\//.test(r) || !/[?&]v=/.test(r)) continue;
    out.add(r);
  }
  return [...out].sort();
}
function hashRefs(list) {                          // FNV-1a, 32 bits, as hex
  let h = 0x811c9dc5;
  const s = list.join("\n");
  for (let i = 0; i < s.length; i++) { h ^= s.charCodeAt(i); h = Math.imul(h, 0x01000193) >>> 0; }
  return ("0000000" + h.toString(16)).slice(-8);
}

// a file of this site, as it is now (not an old copy from the browser's HTTP cache), and only a real answer:
// not an error page, not a sign-in page the server redirected to. Read whole at once (a response left
// unread holds its connection, and the browser has only a few per site), kept as a new response.
async function getOk(url) {
  const r = await fetch(url, { cache: "no-cache", credentials: "same-origin", redirect: "follow" });
  if (!r.ok || r.redirected || r.type !== "basic") { try { await r.body.cancel(); } catch (e) { /* none */ } throw new Error(`${r.status} ${url}`); }
  const buf = await r.arrayBuffer();
  return { buf, response: () => new Response(buf, { headers: { "Content-Type": r.headers.get("Content-Type") || "application/octet-stream" } }) };
}
const inScope = (u) => u.startsWith(SCOPE) && new URL(u).origin === self.location.origin;

// ---- install: keep this version's files ---------------------------------------------------------------
async function install() {
  const cache = await caches.open(APP_CACHE);
  const page = await getOk(PAGE);
  const html = new TextDecoder().decode(page.buf);
  if (hashRefs(appRefs(html)) !== APP) console.warn(`SODA offline: index.html changed while the page was open (${APP})`);
  const want = new Map();                          // address -> needed for the page to start at all
  for (const r of appRefs(html)) want.set(new URL(r, PAGE).href, !/^config\.js(\?|$)/.test(r));
  const files = [[PAGE, page]];
  const seen = new Set(want.keys());
  let queue = [...want.keys()];
  while (queue.length) {
    const got = await Promise.all(queue.map((u) => getOk(u).then((r) => [u, r], (e) => { if (want.get(u)) throw e; return null; })));
    queue = [];
    for (const g of got) {
      if (!g) continue;
      files.push(g);
      if (!/\.js(\?|$)/.test(g[0])) continue;
      // what a script names in turn: workers and the scripts they import ("lasso_worker.js?v=...")
      const text = new TextDecoder().decode(g[1].buf);
      for (const m of text.matchAll(/["'`]([\w-]+(?:\/[\w-]+)*\.js\?v=[\w.-]+)["'`]/g)) {
        const u = new URL(m[1], g[0]).href;
        if (inScope(u) && !seen.has(u)) { seen.add(u); queue.push(u); }
      }
    }
  }
  for (const f of MARKERS) {
    const u = new URL(`markers/${f}`, SCOPE).href;
    try { files.push([u, await getOk(u)]); } catch (e) { /* a site without that list */ }
  }
  await Promise.all(files.map(([u, r]) => cache.put(u, r.response())));
  await keepLibs();
}
// deck.gl and umap-js: the vendor/ copy, else the CDN's; either only when its SHA-384 is the pinned one
async function keepLibs() {
  const cache = await caches.open(LIB_CACHE);
  await Promise.all(LIBS.map(async (L) => {
    if ((await cache.match(L.localUrl)) || (await cache.match(L.cdn))) return;
    try { await cache.put(L.localUrl, await fetch(L.localUrl, { integrity: L.sri, cache: "no-cache", credentials: "same-origin" }).then(whole)); return; }
    catch (e) { /* no vendor/ copy: the CDN's */ }
    try { await cache.put(L.cdn, await fetch(L.cdn, { integrity: L.sri, mode: "cors", credentials: "omit" }).then(whole)); }
    catch (e) { console.warn(`SODA offline: ${L.name} could not be kept (${e.message}); it loads from the network`); }
  }));
}
// a library's answer, read whole and checked (fetch's integrity): kept as a response of its own
async function whole(r) {
  if (!r.ok || r.redirected) { try { await r.body.cancel(); } catch (e) { /* none */ } throw new Error(`${r.status}`); }
  const buf = await r.arrayBuffer();
  return new Response(buf, { headers: { "Content-Type": r.headers.get("Content-Type") || "text/javascript" } });
}

self.addEventListener("install", (e) => {
  e.waitUntil(install().then(() => self.skipWaiting()));
});
// a new version takes over: this scope's older app and library caches go; kept datasets stay
self.addEventListener("activate", (e) => {
  e.waitUntil((async () => {
    for (const name of await caches.keys()) {
      const app = name.startsWith(`soda-app:${SCOPE}|`), lib = name.startsWith(`soda-lib:${SCOPE}|`);
      if ((app && name !== APP_CACHE) || (lib && name !== LIB_CACHE)) await caches.delete(name);
    }
    await self.clients.claim();
  })());
});

// ---- kept datasets: read from their caches (the page writes them) ----------------------------------------
let KEPT = null;                                   // [{ name, root, extras }] once read
async function readKept() {
  const out = [];
  for (const name of await caches.keys()) {
    if (!name.startsWith(DATA_PREFIX)) continue;
    const root = name.slice(DATA_PREFIX.length);
    let extras = [];
    try {
      const m = await (await caches.open(name)).match(root + META);
      const j = m ? await m.json() : null;
      if (j && Array.isArray(j.extras)) extras = j.extras.filter((x) => typeof x === "string");
    } catch (e) { /* no record: the root alone */ }
    out.push({ name, root, extras });
  }
  KEPT = out;
  return out;
}
let keptReady = readKept().catch(() => (KEPT = []));
self.addEventListener("message", (e) => {
  if (e.data && e.data.type === "soda-kept-changed") e.waitUntil(keptReady = readKept().catch(() => (KEPT = [])));
});
const bare = (u) => u.split("#")[0].split("?")[0];
function keptFor(list, href) {
  const u = bare(href);
  return list.find((k) => u.startsWith(k.root) || k.extras.includes(u)) || null;
}
// a Range request answered from the whole file that was kept
async function cut(resp, range) {
  const m = /^bytes=(\d*)-(\d*)$/.exec(range.trim());
  const buf = await resp.arrayBuffer(), size = buf.byteLength;
  let a, b;
  if (m && m[1]) { a = +m[1]; b = m[2] ? Math.min(+m[2], size - 1) : size - 1; }
  else if (m && m[2]) { a = Math.max(0, size - +m[2]); b = size - 1; }
  if (!m || a > b || a >= size) return new Response("", { status: 416, headers: { "Content-Range": `bytes */${size}` } });
  return new Response(buf.slice(a, b + 1), { status: 206, headers: {
    "Content-Type": resp.headers.get("Content-Type") || "application/octet-stream",
    "Content-Range": `bytes ${a}-${b}/${size}`, "Content-Length": String(b - a + 1) } });
}
async function fromKept(req, k) {
  try { return await fetch(req); }
  catch (err) {
    const hit = await caches.match(bare(req.url), { cacheName: k.name });   // match, not open: open would make a removed copy again
    if (!hit) return Response.error();
    const range = req.headers.get("range");
    return range ? cut(hit, range) : hit;
  }
}

// ---- answers --------------------------------------------------------------------------------------------
async function fromApp(e) {
  const req = e.request;
  try {
    const r = await fetch(req);
    if (r.ok && r.type === "basic" && !r.redirected) {
      // keep what this version uses: a file with a ?v= once, a marker list or config.js whenever it changes
      const u = req.url, copy = r.clone();
      const fixed = /[?&]v=/.test(u) && !/\/config\.js\?/.test(u);
      e.waitUntil((async () => {
        // a newer version took over (and deleted this cache) while this worker still answered: keep nothing
        if (!(await caches.has(APP_CACHE))) return;
        const cache = await caches.open(APP_CACHE);
        if (!fixed || !(await cache.match(u))) await cache.put(u, new Response(await copy.arrayBuffer(), { headers: { "Content-Type": copy.headers.get("Content-Type") || "application/octet-stream" } }));
      })().catch(() => {}));
    }
    return r;
  } catch (err) {
    const hit = await caches.match(req.url, { cacheName: APP_CACHE });
    return hit || Response.error();
  }
}
async function fromPage(e) {
  try {
    const r = await fetch(e.request);
    if (r.ok && r.type === "basic" && !r.redirected && /text\/html/i.test(r.headers.get("Content-Type") || "")) {
      // the same version's page as it is now (its ?v= strings unchanged): that is the one to open offline
      const copy = r.clone();
      e.waitUntil(copy.text().then(async (t) => {
        if (hashRefs(appRefs(t)) === APP && await caches.has(APP_CACHE)) await (await caches.open(APP_CACHE)).put(PAGE, new Response(t, { headers: { "Content-Type": "text/html; charset=utf-8" } }));
      }).catch(() => {}));
    }
    return r;
  } catch (err) {
    const hit = await caches.match(PAGE, { cacheName: APP_CACHE });
    return hit || Response.error();
  }
}
async function fromLibs(e) {
  const req = e.request;
  const hit = await caches.match(req.url, { cacheName: LIB_CACHE });
  if (hit) return hit;
  let r;
  // not kept and out of reach (offline, or no vendor/ copy kept): an error status, not a rejected answer,
  // which WebKit reports as a page error; the page then falls back as it does online (libs.js, the worker)
  try { r = await fetch(req); } catch (err) { return new Response("", { status: 504, statusText: "Offline" }); }
  // kept now if it was not at install (the CDN was out of reach then): only a copy the browser checked
  // against the pinned hash (the page's script tags carry it), never an opaque one
  const L = LIBS.find((x) => x.localUrl === req.url || x.cdn === req.url);
  if (L && r.ok && (r.type === "basic" || r.type === "cors") && req.integrity === L.sri) {
    const copy = r.clone();
    e.waitUntil(whole(copy).then(async (w) => { if (await caches.has(LIB_CACHE)) await (await caches.open(LIB_CACHE)).put(req.url, w); }).catch(() => {}));
  }
  return r;
}

self.addEventListener("fetch", (e) => {
  const req = e.request;
  if (req.method !== "GET") return;
  const url = new URL(req.url);
  const same = url.origin === self.location.origin;
  if (req.mode === "navigate") {
    if (same && (url.pathname === SCOPE_PATH || url.pathname === SCOPE_PATH + "index.html")) e.respondWith(fromPage(e));
    return;
  }
  const href = bare(req.url);
  if (LIBS.some((L) => L.localUrl === href || L.cdn === href)) { e.respondWith(fromLibs(e)); return; }
  if (same && inScope(href)) {
    const rel = href.slice(SCOPE.length);
    if (/^sw\.js$/.test(rel)) return;
    if (url.searchParams.has("v") || /^markers\/[\w.-]+\.json$/.test(rel)) { e.respondWith(fromApp(e)); return; }
  }
  if (KEPT) {
    const k = keptFor(KEPT, req.url);
    if (k) e.respondWith(fromKept(req, k));
    return;
  }
  // just started and the list of kept datasets is still being read: wait for it (a few ms, once)
  if (req.cache === "only-if-cached" && req.mode !== "same-origin") return;
  e.respondWith(keptReady.then(() => { const k = keptFor(KEPT || [], req.url); return k ? fromKept(req, k) : fetch(req); }));
});
