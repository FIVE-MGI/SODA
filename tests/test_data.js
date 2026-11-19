"use strict";
// The data the browser tests open: the public example exports (REPRODUCE.md), one place for every test.
//
//   SODA_TEST_DATA=<url or folder>          or  --data-root <url or folder>          (the flag wins)
//       default https://fivestore.blob.core.windows.net/soda/examples/ : the hosted public exports, anonymous read.
//       It holds HOSTED only: xenium_pub, visiumhd_pub, phenocycler_pub, merscope_2.
//   SODA_TEST_PRIVATE_DATA=<url or folder>  or  --private-data-root <url or folder>  (optional)
//       RESTRICTED: cosmx_1, stereoseq_1, geomx_1. Public data whose terms do not allow passing it on (stereoseq_1:
//       "do not redistribute"), so it is kept in a private container (soda-test/examples/, a URL with a read SAS
//       query: ...?sv=...&sig=...) or a local folder. Without it, those ids come from SODA_TEST_DATA when that is a
//       folder or another URL (a local copy of all seven), and not at all from the default hosted copy: the cases
//       that need them are SKIPPED with a one-line message (skipLine), never failed.
//
//   a URL     a copy of the exports (blob storage, or `python web/SpatialViewer/static_serve.py <root> <port>`)
//   a folder  a local copy: <folder>/<id>/manifest.json
//
// A SAS query is a secret: it is added to each request and never printed (describe(), errors and SKIPPED lines show
// the address without its query). Never write one into a file of the repo.
//
// Each test keeps its own small local data server (so it can still refuse, drop or alter answers, as before);
// serve() answers a request from the root that holds that id: read from the folder, or passed on to the URL with
// the Range header (a 206 stays a 206), with open CORS and no caching, as blob storage answers (tests/static_files.js).
//
// Only the seven public ids are ever served or read, whatever else a folder holds: any other path gets 404. A test
// that needs a dataset calls need(): when its manifest cannot be read (no internet, the copy not uploaded yet, a
// folder without it) the test prints a line starting SKIPPED and exits 0, as it does without Playwright or without
// deck.gl from the CDN. CI checks the hosted URL first, so there a missing copy fails instead.
const fs = require("fs");
const path = require("path");
const http = require("http");
const https = require("https");
const { sendFile, CORS } = require("./static_files");

const DEFAULT_ROOT = "https://fivestore.blob.core.windows.net/soda/examples/";
const HOSTED = ["xenium_pub", "visiumhd_pub", "phenocycler_pub", "merscope_2"];
const RESTRICTED = ["cosmx_1", "stereoseq_1", "geomx_1"];
const PUBLIC = HOSTED.concat(RESTRICTED);

const flag = (k) => { const i = process.argv.indexOf(k); return i > 0 ? process.argv[i + 1] : ""; };
// a root: { url, query } for a URL (query = a SAS or "", kept apart), { dir } for a folder
function rootOf(given) {
  if (!given) return null;
  if (/^https?:\/\//i.test(given)) {
    const u = new URL(given), query = u.search;
    u.search = ""; u.hash = "";
    return { url: u.href.replace(/\/*$/, "/"), query };
  }
  return { dir: path.resolve(given) };
}
const MAIN = rootOf(flag("--data-root") || process.env.SODA_TEST_DATA || DEFAULT_ROOT);
const PRIVATE = rootOf(flag("--private-data-root") || process.env.SODA_TEST_PRIVATE_DATA || "");
const ON_HOSTED = !!MAIN.url && MAIN.url === DEFAULT_ROOT;
const show = (R) => (R ? (R.url ? R.url + (R.query ? "?<SAS hidden>" : "") : R.dir) : "(none)");
const ROOT = show(MAIN), IS_URL = !!MAIN.url;

// the root that holds a path "<id>/...", or null when the path may not be read; datasets.json: the main root's
function rootFor(r) {
  const top = r.split("/")[0];
  if (r === "datasets.json" || HOSTED.includes(top)) return MAIN;
  if (!RESTRICTED.includes(top)) return null;
  if (PRIVATE) return PRIVATE;
  return ON_HOSTED ? null : MAIN;          // never asked of the hosted public copy, which does not have them
}
function resolve(rel) {
  const r = String(rel || "").replace(/^\/+/, "");
  if (!r || r.split("/").some((p) => p === ".." || p === ".")) return null;
  const R = rootFor(r);
  return R ? { r, R } : null;
}
const allowed = (rel) => (resolve(rel) ? String(rel).replace(/^\/+/, "") : null);
const urlOf = ({ r, R }) => R.url + r.split("/").map(encodeURIComponent).join("/") + R.query;
const safeUrl = ({ r, R }) => R.url + r;                                  // for messages: no SAS
const fileOf = ({ r, R }) => { const f = path.join(R.dir, ...r.split("/")); return f.startsWith(R.dir + path.sep) ? f : null; };

const agent = { "http:": new http.Agent({ keepAlive: true, maxSockets: 32 }), "https:": new https.Agent({ keepAlive: true, maxSockets: 32 }) };
// one GET (or HEAD) to a URL root: resolves to the response, unread
function upstream(x, method = "GET", headers = {}) {
  return new Promise((resolveP, reject) => {
    const u = new URL(urlOf(x)), mod = u.protocol === "https:" ? https : http;
    const q = mod.request(u, { method, headers: { "User-Agent": "soda-tests", ...headers }, agent: agent[u.protocol], timeout: 120000 }, resolveP);
    q.on("timeout", () => q.destroy(new Error(`timed out: ${safeUrl(x)}`)));
    q.on("error", (e) => reject(new Error(`${e.code || "error"} for ${safeUrl(x)}`)));      // the message could carry the SAS
    q.end();
  });
}

// the whole file as a Buffer, or null when it is not there (404) or may not be read
async function bytes(rel) {
  const x = resolve(rel);
  if (!x) return null;
  if (x.R.dir) { const f = fileOf(x); return f && fs.existsSync(f) && fs.statSync(f).isFile() ? fs.readFileSync(f) : null; }
  const res = await upstream(x);
  const chunks = [];
  for await (const c of res) chunks.push(c);
  if (res.statusCode === 404) return null;
  if (res.statusCode !== 200) throw new Error(`HTTP ${res.statusCode} for ${safeUrl(x)}`);
  return Buffer.concat(chunks);
}
// n bytes from byte o of a file (a Range request to a URL), or null when it is not there
async function range(rel, o, n) {
  const x = resolve(rel);
  if (!x) return null;
  if (x.R.dir) {
    const f = fileOf(x);
    if (!f || !fs.existsSync(f)) return null;
    const fd = fs.openSync(f, "r");
    try { const b = Buffer.alloc(n); fs.readSync(fd, b, 0, n, o); return b; } finally { fs.closeSync(fd); }
  }
  const res = await upstream(x, "GET", { Range: `bytes=${o}-${o + n - 1}` });
  const chunks = [];
  for await (const c of res) chunks.push(c);
  if (res.statusCode === 404) return null;
  const b = Buffer.concat(chunks);
  if (res.statusCode === 206) return b;
  if (res.statusCode === 200) return b.subarray(o, o + n);                  // a host that ignores Range
  throw new Error(`HTTP ${res.statusCode} for ${safeUrl(x)}`);
}
const memo = new Map();
// a JSON file (NaN read as null, as the page does); null when it is not there
async function json(rel) {
  if (!memo.has(rel)) memo.set(rel, bytes(rel).then((b) => (b ? JSON.parse(b.toString("utf-8").replace(/\bNaN\b/g, "null")) : null)));
  return memo.get(rel);
}
async function head(x) { const res = await upstream(x, "HEAD"); res.resume(); return res; }
async function exists(rel) {
  const x = resolve(rel);
  if (!x) return false;
  if (x.R.dir) { const f = fileOf(x); return !!f && fs.existsSync(f) && fs.statSync(f).isFile(); }
  return (await head(x)).statusCode === 200;
}
// a file's size in bytes, or -1
async function size(rel) {
  const x = resolve(rel);
  if (!x) return -1;
  if (x.R.dir) { const f = fileOf(x); return f && fs.existsSync(f) ? fs.statSync(f).size : -1; }
  const res = await head(x);
  return res.statusCode === 200 ? Number(res.headers["content-length"]) : -1;
}

// answer one request for `rel` (a path inside the data, e.g. "visiumhd_pub/manifest.json") as blob storage would;
// opts as sendFile's ({ range: false }: a host that ignores Range)
function serve(req, res, rel, type, opts = {}) {
  const x = resolve(rel);
  const notFound = () => { if (!res.headersSent) res.writeHead(404, { "Access-Control-Allow-Origin": "*" }); res.end("not found"); };
  if (!x) return notFound();
  if (x.R.dir) {
    const f = fileOf(x);
    if (!f || !fs.existsSync(f) || fs.statSync(f).isDirectory()) return notFound();
    return sendFile(req, res, f, type, opts);
  }
  const h = {};
  if (req.headers.range && opts.range !== false) h.Range = req.headers.range;
  upstream(x, req.method === "HEAD" ? "HEAD" : "GET", h).then((up) => {
    if (up.statusCode >= 400) { up.resume(); if (up.statusCode === 404) return notFound(); res.writeHead(up.statusCode, { "Access-Control-Allow-Origin": "*" }); return res.end(); }
    const hd = { ...CORS, "Content-Type": type || up.headers["content-type"] || "application/octet-stream", "Cache-Control": "no-store" };
    for (const k of ["content-length", "content-range", "content-encoding"]) if (up.headers[k]) hd[k.replace(/(^|-)\w/g, (m) => m.toUpperCase())] = up.headers[k];
    if (opts.range !== false) hd["Accept-Ranges"] = "bytes";
    res.writeHead(up.statusCode, hd);
    up.pipe(res);
    res.on("close", () => up.destroy());
  }, (e) => { if (!res.headersSent) res.writeHead(502, { "Access-Control-Allow-Origin": "*" }); res.end(String(e.message || e)); });
}

// can this dataset be read? (a restricted id needs SODA_TEST_PRIVATE_DATA, or a data root other than the hosted one)
async function has(id) {
  if (!resolve(`${id}/manifest.json`)) return false;
  try { return !!(await json(`${id}/manifest.json`)); } catch (e) { return false; }
}
// the line a test prints for a case it leaves out because has(id) is false
function skipLine(id, what) {
  const x = resolve(`${id}/manifest.json`);
  const why = !x ? (RESTRICTED.includes(id) ? "is restricted (its terms do not allow public hosting) and SODA_TEST_PRIVATE_DATA is not set"
    : "is not one of the public example ids") : `cannot be read from ${x.R.dir || x.R.url}`;
  return `SKIPPED ${what}: ${id} ${why}`;
}
// every id here must be readable, or the test says SKIPPED and exits 0
async function need(ids, what = "") {
  for (const id of ids) {
    let m = null, why = "";
    try { m = await json(`${id}/manifest.json`); } catch (e) { why = String(e.message || e).split("\n")[0]; }
    if (!m) {
      const x = resolve(`${id}/manifest.json`);
      console.log(`SKIPPED${what ? ` ${what}` : ""}: cannot read ${id}/manifest.json from ${x ? x.R.dir || x.R.url : "(no root for it)"} (${why || "not there"}); ` +
        "set SODA_TEST_DATA / --data-root (and SODA_TEST_PRIVATE_DATA for cosmx_1, stereoseq_1, geomx_1) to a copy of the public exports");
      process.exit(0);
    }
  }
}
const manifest = (id) => json(`${id}/manifest.json`);
const describe = () => `${show(MAIN)}${PRIVATE ? `; restricted ids from ${show(PRIVATE)}` : ""} (SODA_TEST_DATA${PRIVATE ? " / SODA_TEST_PRIVATE_DATA" : ""})`;

module.exports = { ROOT, IS_URL, PUBLIC, HOSTED, RESTRICTED, DEFAULT_ROOT, allowed, bytes, range, json, exists, size, serve, need, has, skipLine, manifest, describe };
