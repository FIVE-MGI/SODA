// The two libraries SODA does not ship in this repository, pinned by version and SHA-384:
//   deck.gl   draws the map (index.html, the page)
//   umap-js   the UMAP of drawn shapes (lasso_worker.js, a worker)
// Each is read from vendor/<path> next to index.html when that file is there (python prepare.py vendor
// puts it there, checked against the hash below), otherwise from https://cdn.jsdelivr.net/npm/<path>.
// The page's own <script> for vendor/ and the CDN copy written below both carry the SHA-384
// (Subresource Integrity), so a file that differs by one byte is not run; the worker checks the same
// hash with fetch() before importScripts (which has no integrity attribute of its own).
//
// This file is read in three places: by the page (after the vendor/ tag in index.html), by
// lasso_worker.js (importScripts), and by vendor.py, which reads the JSON between the two
// LIBS markers. Changing a version: change it here, in index.html's vendor/ tag (src and integrity)
// and bump the ?v= of this file; tests/smoke_test.py checks that the copies agree.
(function (g) {
  "use strict";
  var LIBS = /* LIBS-JSON */ {
    "deck": { "name": "deck.gl", "version": "9.1.12", "path": "deck.gl@9.1.12/dist.min.js", "bytes": 1472815,
              "sri": "sha384-5HL6K250f0iIGF3yxNktkzN7iXNkWaduwyG3fie9Izlsj535dRJ8qPJ9BbrP28wV" },
    "umap": { "name": "umap-js", "version": "1.4.0", "path": "umap-js@1.4.0/lib/umap-js.min.js", "bytes": 83349,
              "sri": "sha384-5bjAtWKrlgIl+nR7Xjv+eEp7QlOEeyXjMil3ll3fIHu+cusk8DfbJLrsxm3Ftpfo" }
  } /* END-LIBS-JSON */;
  var CDN = "https://cdn.jsdelivr.net/npm/";
  for (var k in LIBS) { LIBS[k].local = "vendor/" + LIBS[k].path; LIBS[k].cdn = CDN + LIBS[k].path; }
  g.SODA_LIBS = LIBS;
  if (typeof document === "undefined" || !document.currentScript) return;      // a worker or the service worker
  // The page: the <script> just above this one asked for vendor/<path>. When it did not run (no vendor/
  // folder, a 404, or a file that failed its hash), deck.gl comes from the CDN instead, as it always did.
  // document.write keeps the order: the CDN copy runs before config.js and the rest, which need window.deck.
  var d = LIBS.deck;
  g.SODA_LIBS_FROM = { deck: g.deck ? "local" : "cdn" };
  if (!g.deck) document.write('<script src="' + d.cdn + '" integrity="' + d.sri + '" crossorigin="anonymous"><\/script>');
})(typeof self !== "undefined" ? self : this);
