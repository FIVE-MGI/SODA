#!/usr/bin/env node
// Unit checks for the pure half of web/SpatialViewer/viewlink.js (the link that restores a view):
// encode / decode round trips in both forms, the length budget and what is left out when a state
// is too big, and that malformed, truncated, oversize or hostile payloads fail safely.
// Needs only node (18 or newer: CompressionStream).   node tests/viewlink_test.js
"use strict";
const path = require("path");
const VL = require(path.join(__dirname, "..", "web", "SpatialViewer", "viewlink.js"));

let fail = 0, n = 0;
const ok = (cond, what) => { n++; if (!cond) { fail++; console.log(`  FAIL  ${what}`); } };
const eq = (got, want, what) => ok(JSON.stringify(got) === JSON.stringify(want), `${what}: got ${JSON.stringify(got)}, want ${JSON.stringify(want)}`);

// a typical state, as viewer.js's sessionState(false) + forLink gives it
const typical = () => ({
  dataset: "demo_visiumhd", layer: "008um",
  sectionClusters: { layer: "008um", resolution: 0.6, names: { "0.6": { 0: "Fascicle", 3: "Epineurium" } }, hidden: [5] },
  qcFilter: { layer: "008um", thresholds: { total: { min: 20, max: null }, pctMt: { min: null, max: 25 } } },
  view: {
    v: 1, mode: "sbs", cam: [5120.5, 3890.2, -1.25], scam: [5000, 3900, -2], win: [1366, 722], tab: "clusters",
    img: { on: true, op: 0.8 }, genes: [{ n: "MARKER_FASC1", c: "#ff3355", r: [0, 12] }, { n: "MARKER_EPI2", c: "#33ddff", r: [1, 7.5], on: false }],
    blend: "heatmap", palette: "viridis", gop: 0.9, ps: 1.5, tb: true, colorBy: "cluster", pick: 3, ov: ["regions"],
  },
});
const circle = (cx, cy, r, k) => Array.from({ length: k }, (_, i) => [cx + r * Math.cos((2 * Math.PI * i) / k), cy + r * Math.sin((2 * Math.PI * i) / k)]);

(async () => {
  // ---- round trip, both forms ---------------------------------------------------------------------
  const st = typical();
  const packed = await VL.encode(st), plain = await VL.encode(st, { compress: false });
  ok(/^1d\.[A-Za-z0-9_-]+$/.test(packed), `the compressed form is "1d." + base64url: ${packed.slice(0, 20)}`);
  ok(/^1j\.[A-Za-z0-9_-]+$/.test(plain), `the plain form is "1j." + base64url: ${plain.slice(0, 20)}`);
  ok(packed.length < plain.length, `compression shortens a typical state (${packed.length} < ${plain.length})`);
  for (const [name, p] of [["compressed", packed], ["plain", plain]]) {
    const d = await VL.decode(p);
    ok(d.ok, `${name}: decodes`);
    eq(d.state, VL.clean(st), `${name}: round trip gives the cleaned state back`);
    eq(d.state.view.genes, st.view.genes, `${name}: genes with colours and ranges survive`);
    eq(d.state.view.cam, st.view.cam, `${name}: the camera survives`);
    eq(d.state.qcFilter.thresholds, st.qcFilter.thresholds, `${name}: QC thresholds survive (null = no bound)`);
    eq(d.state.sectionClusters, st.sectionClusters, `${name}: the clustering, its names and hidden clusters survive`);
  }
  eq(VL.clean(VL.clean(st)), VL.clean(st), "clean() is idempotent");
  // cell sets (clusters merged into named sets) travel with the clustering, as in Save JSON
  const withSets = typical();
  withSets.sectionClusters.sets = { "0.6": [{ ids: [1, 4, 7], prev: "Old name" }, { ids: [2, 3], prev: null }] };
  eq((await VL.decode(await VL.encode(VL.forLink(withSets)))).state.sectionClusters.sets, withSets.sectionClusters.sets, "cell sets survive the link");
  eq(VL.clean({ dataset: "d", sectionClusters: { layer: "l", resolution: 1, names: {}, hidden: [], sets: JSON.parse('{"__proto__": [{"ids": [1, 2]}], "1": [{"ids": [5]}, {"ids": [3, "x", 9], "prev": "<b>n</b>"}, "junk"]}') } }).sectionClusters.sets,
     { 1: [{ ids: [3, 9], prev: "bn/b" }] }, "cell sets: unsafe keys, one-member sets and non-numbers are dropped, text is stripped");
  // non-ASCII text (µ, accents, CJK) survives both forms
  const uni = { dataset: "d", sectionClusters: { layer: "l", resolution: 1, names: { 1: { 0: "Müller glia 視細胞 – β" } }, hidden: [] } };
  for (const c of [true, false]) eq((await VL.decode(await VL.encode(uni, { compress: c }))).state.sectionClusters.names["1"]["0"], "Müller glia 視細胞 – β", `non-ASCII names survive (compress ${c})`);

  // ---- the budget: what is left out when the state is too big -----------------------------------------
  const base = "https://example.org/SODA/index.html?ds=demo_visiumhd".length;
  let r = await VL.fit(st, base, VL.BUDGET);
  ok(r.tier === 0 && !r.left.length && !r.over, `a typical state fits whole (payload ${r.payload.length} characters)`);
  ok(base + 3 + r.payload.length <= VL.BUDGET, `a typical link is within the ${VL.BUDGET} character budget (${base + 3 + r.payload.length})`);
  // a few circles still fit whole
  const few = { ...typical(), shapes: [0, 1, 2].map((i) => ({ name: `Shape ${i + 1}`, color: "#e6194b", kind: "circle", points: circle(4000 + 900 * i, 3000, 350, 64).map((p) => [Math.round(p[0] * 100) / 100, Math.round(p[1] * 100) / 100]) })) };
  r = await VL.fit(few, base, VL.BUDGET);
  ok(!r.over && r.tier <= 1, `three circles fit in a link (tier ${r.tier}, payload ${r.payload.length})`);
  // long freehand outlines: thinned first
  const rnd = (() => { let s = 7; return () => { s = (s * 16807) % 2147483647; return s / 2147483647; }; })();
  const wobble = (k) => circle(5000, 4000, 1500, k).map((p) => [p[0] + 40 * rnd(), p[1] + 40 * rnd()]);
  const free = { ...typical(), shapes: [0, 1].map((i) => ({ name: `Region ${i}`, color: "#3cb44b", kind: "free", points: wobble(900) })) };
  r = await VL.fit(free, base, VL.BUDGET);
  ok(r.tier === 1 && !r.over && r.left.some((t) => /fine detail/.test(t)), `long freehand outlines are thinned to fit (tier ${r.tier}, left: ${r.left})`);
  let back = (await VL.decode(r.payload)).state;
  ok(back.shapes.length === 2 && back.shapes.every((s) => s.points.length === 40), "thinned outlines keep 40 points each");
  // many shapes, and a shape made from a selection (a list of bins): the shapes are left out
  const many = { ...typical(), shapes: Array.from({ length: 30 }, (_, i) => ({ name: `Region ${i}`, color: "#4363d8", kind: "free", points: wobble(600) }))
    .concat([{ name: "Selection", color: "#f58231", kind: "sel", points: [], rows: Array.from({ length: 60000 }, (_, i) => i * 3) }]) };
  r = await VL.fit(many, base, VL.BUDGET);
  ok(r.tier >= 2 && !r.over && r.left.includes("the drawn shapes"), `too many shapes are left out, and said so (tier ${r.tier}, left: ${r.left})`);
  back = (await VL.decode(r.payload)).state;
  ok(!back.shapes && back.view.genes.length === 2 && back.sectionClusters.names["0.6"]["0"] === "Fascicle", "without the shapes the rest of the view is intact");
  ok(base + 3 + r.payload.length <= VL.BUDGET, `the oversize state still gives a link within the budget (${base + 3 + r.payload.length})`);
  // hundreds of named clusters: the names go last
  const named = typical();
  named.sectionClusters.names = { "0.6": Object.fromEntries(Array.from({ length: 400 }, (_, i) => [i, `Cell type number ${i} ${(i * 7919).toString(36)} ${(i * 104729).toString(36)}`])) };
  named.shapes = many.shapes;
  r = await VL.fit(named, base, VL.BUDGET);
  ok(r.tier === 3 && r.left.includes("the names given to clusters") && !r.over, `too many cluster names are left out last (tier ${r.tier})`);
  // the worst case of what is never left out: every list full
  const worst = { dataset: "demo_dataset_with_a_long_name_0123456789", layer: "square_016um",
    sectionClusters: { layer: "square_016um", resolution: "seurat_clusters", names: {}, hidden: Array.from({ length: 40 }, (_, i) => i) },
    qcFilter: { layer: "square_016um", thresholds: { total: { min: 20, max: 90000 }, genes: { min: 10, max: 9000 }, pctMt: { min: 0.5, max: 25 } } },
    view: { v: 1, mode: "sbs", cam: [123456.7, 98765.4, -5.123], scam: [123456.7, 98765.4, -5.123], win: [2560, 1300], tab: "analyze", img: { on: true, op: 0.66 },
      ch: Array.from({ length: 40 }, (_, i) => ({ n: `Marker-${i} (Cy${i})`, on: i < 12, c: `#${(0x123456 + i * 0x070503).toString(16).slice(0, 6)}`, w: [i, 200 + (i % 50)] })),
      genes: Array.from({ length: 16 }, (_, i) => ({ n: `GENE${i}LONGNAME`, c: `#${(0xa1b2c3 + i * 0x010203).toString(16).slice(0, 6)}`, r: [i * 0.5, 10.5 + i] })),
      mol: { genes: Array.from({ length: 24 }, (_, i) => ({ n: `MOLGENE${i}`, c: `#${(0x334455 + i * 0x050709).toString(16).slice(0, 6)}` })), size: 1.5, op: 0.8 },
      blend: "normal", palette: "genes", gop: 0.85, ps: 2.2, tb: true, colorBy: "celltype", pick: 12, ov: ["regions", "cells", "rois"], qcOutline: true } };
  r = await VL.fit(worst, base, VL.BUDGET);
  ok(r.tier === 0 && !r.left.length, "the worst-case view has no shapes or names to leave out, and nothing is claimed to be");
  ok(r.payload.length < 1.25 * VL.BUDGET, `the worst-case view (40 changed stains, 16 genes, 24 molecule genes) stays near the budget (${r.payload.length})`);
  console.log(`  lengths: typical payload ${packed.length} (plain ${plain.length}); worst-case view payload ${r.payload.length} (${r.over ? "over" : "within"} the ${VL.BUDGET} budget with a ${base}-character address)`);
  ok(r.payload.length < VL.MAX_PAYLOAD, `the worst-case view is far below what the decoder accepts (${r.payload.length} < ${VL.MAX_PAYLOAD})`);
  ok((await VL.decode(r.payload)).ok, "the worst-case view decodes");

  // ---- malformed, truncated, oversize: a message, never an exception ------------------------------------
  const bad = ["", "garbage", "1d.", "1d.!!!!", "1x.AAAA", "0d.AAAA", "d.AAAA", "1d", "1j.e30x$", "1j." + Buffer.from("not json").toString("base64url"),
    "1j." + Buffer.from("[1,2,3]").toString("base64url"), "1j." + Buffer.from("null").toString("base64url"), "1j." + Buffer.from('"text"').toString("base64url"),
    "1j." + Buffer.from([0xff, 0xfe, 0x80]).toString("base64url"), "1d." + Buffer.from("this is not deflate data at all").toString("base64url")];
  for (const p of bad) {
    let d = null, threw = false;
    try { d = await VL.decode(p); } catch (e) { threw = true; }
    ok(!threw && d && d.ok === false && typeof d.message === "string" && d.message.length > 10, `malformed payload fails with a message: ${JSON.stringify(p.slice(0, 24))}`);
  }
  for (const v of [null, undefined, 42, {}, []]) { const d = await VL.decode(v); ok(d.ok === false, `a payload that is not text fails safely: ${JSON.stringify(v)}`); }
  let truncOk = 0, truncBad = 0;
  for (let cut = 1; cut < packed.length; cut += Math.max(1, Math.floor(packed.length / 60))) {
    let d = null, threw = false;
    try { d = await VL.decode(packed.slice(0, cut)); } catch (e) { threw = true; }
    if (threw || !d || (d.ok && JSON.stringify(d.state) !== JSON.stringify(VL.clean(st)) && !d.state)) truncBad++; else truncOk++;
    ok(!threw && d && (d.ok === false || typeof d.state === "object"), `a link cut at ${cut} of ${packed.length} characters does not throw`);
    if (d && d.ok) ok(false, `a link cut at ${cut} characters was read as a whole view`);
  }
  console.log(`  truncated links: ${truncOk} cuts tried, ${truncBad} failed unsafely`);
  let d = await VL.decode("9d.AAAA");
  ok(!d.ok && d.why === "newer" && /newer SODA/.test(d.message), "a newer format version says so");
  d = await VL.decode("1j." + "A".repeat(VL.MAX_PAYLOAD));
  ok(!d.ok && d.why === "toobig", "a payload over the length limit is not read");
  // a small payload that unpacks to something huge (a "zip bomb"): stopped at the limit
  const zeros = Buffer.alloc(40 * 1024 * 1024, 0x20);
  const cs = new CompressionStream("deflate-raw"), w = cs.writable.getWriter();
  w.write(zeros); w.close();
  const bomb = Buffer.from(await new Response(cs.readable).arrayBuffer());
  ok(bomb.length < 60000, `the test bomb is small (${bomb.length} bytes for 40 MB)`);
  const t0 = Date.now();
  d = await VL.decode("1d." + bomb.toString("base64url").slice(0, VL.MAX_PAYLOAD - 10));
  ok(!d.ok, `a payload that unpacks past ${VL.MAX_JSON} bytes is refused (${Date.now() - t0} ms)`);

  // ---- hostile or sloppy content: only known keys of the right type come through -----------------------
  const evil = JSON.parse(JSON.stringify(typical()));
  evil.view.genes = [{ n: '<img src=x onerror="alert(1)">', c: "red; background:url(//evil.test/x)", r: ["a", "b"] }, { n: "OK1", c: "#ABCDEF", r: [2, 1] }, 17, null,
    { n: "x".repeat(500), c: "#12345", r: [0, 5] }];
  evil.view.cam = [NaN, 5, 1]; evil.view.scam = ["1", "2", "3"]; evil.view.ucam = [1, 2, 1e9];
  evil.view.mode = "javascript:alert(1)"; evil.view.tab = "__proto__"; evil.view.onload = "alert(1)"; evil.view.url = "https://evil.test/data/";
  evil.view.ch = [{ n: "DAPI", on: 1, c: "#00ff00", w: [300, -5] }, { n: "CD8", on: "yes", w: [10, 200] }];
  evil.view.rna = { colorBy: "leiden<script>", genes: Array.from({ length: 500 }, (_, i) => `G${i}`), bin: 3, ps: 99 };
  evil.sectionClusters.names = JSON.parse('{"__proto__": {"polluted": "yes"}, "constructor": {"1": "x"}, "0.6": {"1": "<b>bold</b>", "x": "bad key", "2": 5}}');
  evil.sectionClusters.hidden = [1, "2", -3, 1e9, 4.5];
  evil.qcFilter.thresholds = JSON.parse('{"__proto__": {"min": 1}, "total": {"min": "20", "max": 50}, "bad key!": {"min": 1}, "genes": {"min": 5, "max": null}}');
  evil.shapes = [{ name: "<script>alert(1)</script>", color: "url(javascript:alert(1))", kind: "free", points: [[1, 2], [3, "x"], [5, 6], [7, 8], [9, 10]] },
    { name: "two points", points: [[1, 2], [3, 4]] }, { rows: [1, 2, -1, 2.5, "7", 9] }];
  evil.dataset = "demo<script>"; evil.base = "https://evil.test/"; evil.tokenEndpoint = "https://evil.test/token"; evil.__proto__ = { polluted: true };
  const c = VL.clean(evil), flat = JSON.stringify(c);
  ok(!/[<>]/.test(flat), "no angle bracket survives clean()");
  ok(!/evil\.test|javascript:|onerror|onload/.test(flat.replace(/"n":"img src=x onerror=alert\(1\)"/, "")), "no address or handler survives clean() outside a name that is only looked up");
  ok(c.view.genes.length === 3 && c.view.genes[0].c === undefined && c.view.genes[0].r === undefined, "a gene with a bad colour and range keeps only its (stripped) name");
  ok(c.view.genes[1].c === "#abcdef" && c.view.genes[1].r === undefined, "a colour is lower-cased; a range with max < min is dropped");
  ok(c.view.genes[2].n.length === 80, "a name is cut to 80 characters; entries that are not genes are dropped");
  ok(c.view.cam === undefined && c.view.scam === undefined && c.view.ucam[2] === 40, "cameras: NaN and text dropped, zoom clamped");
  ok(c.view.mode === undefined && c.view.tab === undefined && c.view.onload === undefined && c.view.url === undefined, "unknown modes, tabs and keys are dropped");
  eq(c.view.ch, [{ n: "DAPI", on: true, c: "#00ff00" }, { n: "CD8", on: false, w: [10, 200] }], "stains: a bad window is dropped, on is strictly true/false");
  ok(c.view.rna.genes.length === 24 && c.view.rna.bin === undefined && c.view.rna.ps === 6 && c.view.rna.colorBy === "leidenscript", "single-cell view: gene list capped, bad values dropped or clamped");
  eq(Object.keys(c.sectionClusters.names), ["0.6"], "cluster names: __proto__ and constructor keys are dropped");
  eq(c.sectionClusters.names["0.6"], { 1: "bbold/b" }, "cluster names: text stripped, non-numeric keys and non-text values dropped");
  eq(c.sectionClusters.hidden, [1], "hidden clusters: whole numbers in range only");
  eq(c.qcFilter.thresholds, { total: { min: null, max: 50 }, genes: { min: 5, max: null } }, "QC thresholds: known-shaped keys, numbers only");
  ok(c.shapes.length === 2 && c.shapes[0].name === "scriptalert(1)/script" && c.shapes[0].color === undefined && c.shapes[0].points.length === 4, "shapes: text stripped, bad colour and bad points dropped, too-short outlines dropped");
  eq(c.shapes[1].rows, [1, 2, 9], "shape rows: non-negative whole numbers only");
  ok(c.base === undefined && c.tokenEndpoint === undefined && ({}).polluted === undefined && c.polluted === undefined, "no address key is read, and Object.prototype is untouched");
  ok(c.dataset === "demoscript", "the dataset id is only compared with the dataset that is open");
  eq(VL.clean({ dataset: "d" }), { dataset: "d" }, "a state with nothing but a dataset is fine (every field is optional)");
  eq(VL.cleanView({}), { v: 1 }, "an empty view is fine");
  eq(VL.cleanView({ v: 7, future: { a: 1 }, mode: "tissue" }), { v: 7, mode: "tissue" }, "a view from a later version keeps what this version knows");
  ok(VL.clean(null) === null && VL.clean([1]) === null && VL.cleanView("x") === null, "clean() of a non-object is null");
  // forLink(): defaults are not written
  eq(VL.forLink({ tool: "SODA circle -> UMAP", version: 1, dataset: "d", layer: "l", seed: 42, mode: "bins", resolution: "auto", shapes: [], sectionClusters: null, qcFilter: null, view: { v: 1 } }),
     { dataset: "d", layer: "l", view: { v: 1 } }, "forLink() leaves the defaults out");
  eq(VL.forLink({ dataset: "d", seed: 7, mode: "shapes", resolution: 0.8, qcFilter: { layer: "l", thresholds: { total: { min: 5, max: null } }, nRows: 9, kept: 8, removed: 1, rule: "..." } }),
     { dataset: "d", seed: 7, mode: "shapes", resolution: 0.8, qcFilter: { layer: "l", thresholds: { total: { min: 5, max: null } } } }, "forLink() keeps what differs, and only the thresholds of the QC filter");

  // ---- combined shapes (selections.js): a link carries the recipe (def), not the bins ----------------------------
  const ring = (cx, cy, rr) => circle(cx, cy, rr, 32).map((p) => [Math.round(p[0] * 100) / 100, Math.round(p[1] * 100) / 100]);
  const def = { op: "sub", parts: [{ op: "union", parts: [{ points: ring(4000, 3000, 400) }, { points: ring(4500, 3000, 400) }] },
    { geo: [[ring(4200, 3000, 150), ring(4200, 3000, 50)]] }] };
  const comb = { ...typical(), shapes: [{ name: "A + B − C", color: "#911eb4", kind: "rows", points: [], def }] };
  r = await VL.fit(comb, base, VL.BUDGET);
  back = (await VL.decode(r.payload)).state;
  ok(r.tier === 0 && back.shapes.length === 1, `a combined shape fits a link as its recipe (tier ${r.tier}, payload ${r.payload.length})`);
  eq(back.shapes[0].def, def, "the recipe survives the round trip (operations, outlines, polygon with a hole)");
  ok(back.shapes[0].rows === undefined && back.shapes[0].kind === undefined, "no bins in the link; the kind \"rows\" is not a link kind (selections.js sets it again)");
  const longDef = { ...typical(), shapes: [{ name: "Big", color: "#911eb4", points: [], def: { op: "inter", parts: [{ points: wobble(900) }, { points: wobble(900) }] } }] };
  r = await VL.fit(longDef, base, VL.BUDGET);
  back = (await VL.decode(r.payload)).state;
  ok(r.tier === 1 && !r.over && r.left.some((t) => /fine detail/.test(t)) && back.shapes[0].def.parts.every((p) => p.points.length === 40),
     `a combined shape of long outlines is thinned like a drawn one (tier ${r.tier}, left: ${r.left})`);
  let deep = { points: ring(0, 0, 10) };
  for (let i = 0; i < 9; i++) deep = { op: "union", parts: [deep, { points: ring(0, 0, 5) }] };
  const badDefs = VL.clean({ dataset: "d", shapes: [
    { name: "deep", def: deep },
    { name: "xor", def: { op: "xor", parts: [{ points: ring(0, 0, 9) }, { points: ring(0, 0, 8) }] } },
    { name: "one part", def: { op: "union", parts: [{ points: ring(0, 0, 9) }] } },
    { name: "huge", def: { op: "union", parts: Array.from({ length: 6 }, () => ({ points: wobble(4000) })) } },
    { name: "rows leaf", def: { op: "union", parts: [{ rows: [1, 2, 3] }, { points: ring(0, 0, 9) }] } },
    { name: "proto", def: JSON.parse(`{"op":"union","parts":[{"points":${JSON.stringify(ring(0, 0, 9))}},{"__proto__":{"points":${JSON.stringify(ring(0, 0, 9))}}}]}`) },
    { name: "ok", def: { op: "inter", parts: [{ points: ring(0, 0, 9) }, { points: [[0, 0], [1, "x"], [5, 0], [5, 5], [0, 5]] }] } },
  ] });
  eq(badDefs.shapes.map((s) => s.name), ["ok"], "a recipe too deep, with an unknown operation, one part, too many points, bins or a borrowed key is dropped (with its shape)");
  ok(badDefs.shapes[0].def.parts[1].points.length === 4 && ({}).points === undefined, "a recipe's bad points are dropped like an outline's; Object.prototype is untouched");

  console.log(fail ? `viewlink_test: ${fail} of ${n} checks FAILED` : `viewlink_test: ${n} checks passed`);
  process.exit(fail ? 1 : 0);
})().catch((e) => { console.log(`viewlink_test: crashed: ${e && e.stack || e}`); process.exit(1); });
