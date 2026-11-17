#!/usr/bin/env node
// The SODA desktop app's pages in headless Chromium (run by tests/app_test.py --browser).
//
//   node tests/app_browser_test.js --viewer <viewer url from `SODA view`> --datasets id1,id2 [--shots <dir>]
//   node tests/app_browser_test.js --home <http://127.0.0.1:PORT/app/> [--shots <dir>]
//
// viewer: the page opens the export's datasets.json through the app's own server, lists every
//         dataset, opens the first, fetches and draws image tiles; no console error, no page error.
// home:   the home page renders its three sections and its status line, with no console error.
//         (Nothing is clicked: Browse... would open a real folder picker on the desktop.)
// Playwright is not a dependency of this repo: point NODE_PATH at a node_modules that has it
// (PLAYWRIGHT_BROWSERS_PATH if its browsers are elsewhere). Without it, or without internet for
// deck.gl, the test says SKIPPED and exits 0.
"use strict";
const fs = require("fs");
const path = require("path");

const arg = (k, d) => { const i = process.argv.indexOf(k); return i > 0 ? process.argv[i + 1] : d; };
const VIEWER = arg("--viewer", ""), HOME = arg("--home", ""), SHOTS = arg("--shots", "");
const WANT = arg("--datasets", "").split(",").filter(Boolean);

let pw = null;
for (const base of [...(process.env.NODE_PATH || "").split(path.delimiter)].filter(Boolean)) {
  try { pw = require(path.join(base, "playwright")); break; } catch (e) { /* next */ }
}
if (!pw) { try { pw = require("playwright"); } catch (e) { /* none */ } }
if (!pw) { console.log("SKIPPED: playwright is not installed"); process.exit(0); }

let fails = 0;
const ok = (cond, what) => { console.log(`${cond ? "  ok  " : "  FAIL"} ${what}`); if (!cond) fails++; };
const until = (page, fn, a, ms = 30000) => page.waitForFunction(fn, a, { timeout: ms, polling: 200 }).then(() => true, () => false);

(async () => {
  let browser;
  try { browser = await pw.chromium.launch(); }
  catch (e) { console.log(`SKIPPED: cannot launch Chromium (${String(e.message).split("\n")[0]})`); process.exit(0); }
  const ctx = await browser.newContext({ viewport: { width: 1400, height: 900 } });
  const page = await ctx.newPage();
  const errors = [], failed = [];
  // gates.json is optional (a lab's published gates); the page asks for it and the browser logs the 404,
  // as tests/compat_test.js allows too
  const optional = (u) => /\/gates\.json(\?|$)/.test(u || "");
  page.on("console", (m) => {
    if (m.type() !== "error") return;
    const u = (m.location() || {}).url || "";
    if (/Failed to load resource/.test(m.text()) && optional(u)) return;
    errors.push(`${m.text()} ${u}`);
  });
  page.on("pageerror", (e) => errors.push(String(e.message || e)));
  page.on("requestfailed", (r) => { if (!optional(r.url())) failed.push(r.url()); });
  page.on("response", (r) => { if (r.status() >= 400 && !optional(r.url())) failed.push(`${r.status()} ${r.url()}`); });
  try {
    if (VIEWER) {
      await page.goto(VIEWER, { waitUntil: "load" });
      if (!(await page.evaluate(() => !!window.deck))) { console.log("SKIPPED: deck.gl did not load from the CDN (no internet?)"); process.exit(0); }
      const opened = await until(page, () => { const sv = window.__spatialViewer; return sv && sv.state && sv.state.m && sv.state.m.id; });
      const r = await page.evaluate(() => {
        const sv = window.__spatialViewer;
        const sel = document.querySelector("select#dataset, select#ds, select[id*=ataset]");
        return { id: sv && sv.state.m && sv.state.m.id, list: (sv && sv.state.datasets || []).map((d) => d.id),
                 picker: sel ? sel.options.length : null };
      });
      ok(opened, `the viewer opened a dataset (${r.id})`);
      ok(WANT.every((w) => r.list.includes(w)), `every dataset listed from datasets.json: ${r.list.join(", ")}`);
      const tiles = await until(page, () => window.__spatialViewer.diag().cached > 0);
      ok(tiles, "image tiles fetched and decoded");
      await page.waitForTimeout(1500);
      if (SHOTS) await page.screenshot({ path: path.join(SHOTS, "viewer.png") });
    }
    if (HOME) {
      await page.goto(HOME, { waitUntil: "load" });
      await page.waitForTimeout(1000);
      const r = await page.evaluate(() => ({ h2: [...document.querySelectorAll("h2")].map((h) => h.textContent),
                                             key: document.getElementById("keyinfo").textContent }));
      ok(r.h2.length === 3, `home page sections: ${r.h2.join(" / ")}`);
      ok(!!r.key, "home page status line filled from /app/api/status");
      if (SHOTS) await page.screenshot({ path: path.join(SHOTS, "home.png"), fullPage: true });
    }
    ok(!errors.length, `no console or page error${errors.length ? ": " + errors.slice(0, 3).join(" / ") : ""}`);
    ok(!failed.length, `no failed request${failed.length ? ": " + failed.slice(0, 3).join(" / ") : ""}`);
  } finally {
    await browser.close();
  }
  console.log(fails ? `${fails} failed` : "all passed");
  process.exit(fails ? 1 : 0);
})();
