// Offline use ("Keep this dataset", offline.js) was removed in SODA 0.10.0. A reader who used it still has
// the old service worker and its caches; the browser fetches this file when it checks that worker for an
// update, so this one takes its place, deletes the caches it kept (soda-app:, soda-lib:, soda-data: for this
// page's folder) and unregisters itself. It answers no requests: the page works as if there never was one.
"use strict";
const SCOPE = self.registration.scope;
const OURS = ["soda-app:", "soda-lib:", "soda-data:"].map((p) => p + SCOPE + "|");
self.addEventListener("install", () => self.skipWaiting());
self.addEventListener("activate", (e) => {
  e.waitUntil((async () => {
    try { for (const n of await caches.keys()) if (OURS.some((p) => n.startsWith(p))) await caches.delete(n); } catch (err) { /* none */ }
    try { await self.registration.unregister(); } catch (err) { /* none */ }
  })());
});
