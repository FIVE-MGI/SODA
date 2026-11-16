"use strict";
// One static file as blob storage would answer it, for the test servers (compat_test.js,
// phone_browser_test.js): open CORS, no caching, and HTTP Range requests (206 Partial Content),
// which the viewer uses to read one gene out of an export's gene packs (web/SpatialViewer/genepack.py).
//   sendFile(req, res, file, type, { range: false })   a host that ignores Range: always the whole file
const fs = require("fs");

const CORS = {
  "Access-Control-Allow-Origin": "*",
  "Access-Control-Expose-Headers": "Content-Range, Content-Length, Accept-Ranges",
};

// a CORS preflight (a browser that does not count "Range: bytes=a-b" as a simple header asks first)
function preflight(req, res) {
  res.writeHead(204, { ...CORS, "Access-Control-Allow-Methods": "GET, HEAD, OPTIONS",
    "Access-Control-Allow-Headers": req.headers["access-control-request-headers"] || "Range", "Access-Control-Max-Age": "3600" });
  res.end();
}

function sendFile(req, res, file, type, opts = {}) {
  const size = fs.statSync(file).size;
  const head = { ...CORS, "Content-Type": type || "application/octet-stream", "Cache-Control": "no-store" };
  const m = opts.range === false ? null : /^bytes=(\d*)-(\d*)$/.exec((req.headers.range || "").trim());
  if (m && (m[1] || m[2])) {
    let start, end;
    if (m[1]) { start = +m[1]; end = Math.min(m[2] ? +m[2] : size - 1, size - 1); }
    else { start = Math.max(0, size - +m[2]); end = size - 1; }
    if (start >= size || start > end) { res.writeHead(416, { ...head, "Content-Range": `bytes */${size}` }); res.end(); return; }
    res.writeHead(206, { ...head, "Accept-Ranges": "bytes", "Content-Range": `bytes ${start}-${end}/${size}`, "Content-Length": String(end - start + 1) });
    if (req.method === "HEAD") { res.end(); return; }
    fs.createReadStream(file, { start, end }).pipe(res);
    return;
  }
  res.writeHead(200, { ...head, ...(opts.range === false ? {} : { "Accept-Ranges": "bytes" }), "Content-Length": String(size) });
  if (req.method === "HEAD") { res.end(); return; }
  fs.createReadStream(file).pipe(res);
}

module.exports = { sendFile, preflight, CORS };
