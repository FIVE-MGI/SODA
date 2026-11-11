#!/usr/bin/env node
// Unit checks for the pure parts of web/SpatialViewer/diag.js (the ?diag=1 panel): the
// function that strips tokens and passwords from anything shown, the user-agent reader and
// the tile sniffer. Needs only node.   node tests/diag_test.js [<export root>]
"use strict";
const fs = require("fs");
const path = require("path");
const { clean, parseUA, imageInfo } = require(path.join(__dirname, "..", "web", "SpatialViewer", "diag.js"));

let fail = 0, n = 0;
const ok = (cond, what) => { n++; if (!cond) { fail++; console.log(`  FAIL  ${what}`); } };
const eq = (got, want, what) => ok(got === want, `${what}: got ${JSON.stringify(got)}, want ${JSON.stringify(want)}`);

// ---- clean(): nothing after "?" survives, wherever the address sits in the text
const SAS = "sv=2024-11-04&st=2026-11-11T00%3A00%3A00Z&se=2026-11-11T01%3A00%3A00Z&sr=c&sp=r&spr=https&sig=AbCdEf123%2B456%3D";
eq(clean(`https://acct.blob.core.windows.net/box/ds1/tiles/0/0/0/0.png?${SAS}`), "https://acct.blob.core.windows.net/box/ds1/tiles/0/0/0/0.png", "signed tile address");
eq(clean(`403 https://acct.blob.core.windows.net/box/ds1/manifest.json?${SAS}`), "403 https://acct.blob.core.windows.net/box/ds1/manifest.json", "error text with an address");
eq(clean(`could not load (https://x.test/a.js?v=1) twice: "https://x.test/b?${SAS}"`), `could not load (https://x.test/a.js) twice: "https://x.test/b"`, "addresses inside brackets and quotes");
eq(clean(`/api/token?password=hunter2&x=1 failed`), "/api/token failed", "relative address");
for (const leak of [`retry with ${SAS}`, `sig=AbCdEf123%2B456%3D`, `token=abc.def.ghi`, `password=hunter2`, `Authorization: Basic dmlld2VyOmh1bnRlcjI=`, `Bearer eyJhbGciOiJIUzI1NiJ9.e30.abc`]) {
  const out = clean(leak);
  ok(!/AbCdEf123|hunter2|dmlld2Vy|eyJhbGci|abc\.def/.test(out), `secret survives clean(): ${JSON.stringify(leak)} -> ${JSON.stringify(out)}`);
}
eq(clean(null), "", "null");
eq(clean("HTTP 200, 761,923 bytes, image/png"), "HTTP 200, 761,923 bytes, image/png", "plain text is left alone");

// ---- parseUA()
let u = parseUA("Mozilla/5.0 (iPhone; CPU iPhone OS 17_5_1 like Mac OS X) AppleWebKit/605.1.15 (KHTML, like Gecko) Version/17.5 Mobile/15E148 Safari/604.1", 5);
eq(`${u.device}|${u.os}|${u.browser}|${u.ios}`, "iPhone|iOS 17.5.1|Safari 17.5|true", "iPhone Safari");
u = parseUA("Mozilla/5.0 (iPhone; CPU iPhone OS 18_6 like Mac OS X) AppleWebKit/605.1.15 (KHTML, like Gecko) Version/26.0 Mobile/15E148 Safari/604.1", 5);
eq(`${u.os}|${u.browser}`, "iOS 18.6|Safari 26.0", "iOS 26 (frozen OS number, true Safari number)");
u = parseUA("Mozilla/5.0 (iPhone; CPU iPhone OS 16_7 like Mac OS X) AppleWebKit/605.1.15 (KHTML, like Gecko) CriOS/120.0.6099.119 Mobile/15E148 Safari/604.1", 5);
ok(u.ios && /^CriOS 120/.test(u.browser) && /WebKit/.test(u.browser), `Chrome on iOS is read as WebKit: ${u.browser}`);
u = parseUA("Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/605.1.15 (KHTML, like Gecko) Version/17.4 Safari/605.1.15", 5);
ok(u.ios && /iPad/.test(u.device), `iPad asking for desktop sites: ${u.device}`);
u = parseUA("Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/605.1.15 (KHTML, like Gecko) Version/17.4 Safari/605.1.15", 0);
eq(`${u.device}|${u.browser}|${u.ios}`, "Mac|Safari 17.4|false", "Mac Safari");
u = parseUA("Mozilla/5.0 (Linux; Android 14; Pixel 7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/141.0.7390.37 Mobile Safari/537.36", 5);
eq(`${u.device}|${u.browser}|${u.ios}`, "Android|Chrome 141.0.7390.37|false", "Android Chrome");

// ---- imageInfo(): on made-up headers, and on real demo tiles when an export is at hand
const png = new Uint8Array(33); png.set([0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a, 0, 0, 0, 13, 0x49, 0x48, 0x44, 0x52, 0, 0, 4, 0, 0, 0, 4, 0, 8, 0, 0, 0, 0]);
eq(imageInfo(png), "PNG 1024 x 1024, 8-bit grey", "grey PNG header");
eq(imageInfo(new TextEncoder().encode("<html><body>Sign in</body></html>")), "NOT an image: text (an error page?)", "an error page");
eq(imageInfo(new Uint8Array([1, 2, 3, 4])), "NOT a PNG or JPEG", "junk");
const root = process.argv[2] || path.join(__dirname, "..", "demo", "export");
const tiles = [["demo_visiumhd/tiles/0/0/0/0.png", /^PNG \d+ x \d+, 8-bit grey$/], ["demo_multiplex/tiles/0/0/0/0.jpg", /^JPEG \d+ x \d+, 8-bit grey$/]];
for (const [rel, want] of tiles) {
  const p = path.join(root, rel);
  if (!fs.existsSync(p)) { console.log(`  (no ${rel} under ${root}: skipped)`); continue; }
  const got = imageInfo(new Uint8Array(fs.readFileSync(p).subarray(0, 4096)));
  ok(want.test(got), `${rel}: read as ${JSON.stringify(got)}`);
}

console.log(fail ? `diag_test: ${fail} of ${n} checks FAILED` : `diag_test: ${n} checks passed`);
process.exit(fail ? 1 : 0);
