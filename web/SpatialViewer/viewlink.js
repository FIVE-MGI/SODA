// A link that restores the view: "Copy link to this view" in the Share tab.
//
// The state is the one Save JSON writes (viewer.js: sessionState), without the UMAP run's
// parameters and result, plus a `view` block made here: the camera, the mode, stains, genes,
// molecules, outlines, colour-by, the open tab. One serialiser, two outputs: the file and the link.
// Load JSON and a link both go back in through viewer.js's sessionApply, which calls apply() here
// for the `view` block.
//
// The link:   <page>?ds=<dataset>[&static=...][&qc=...]#v=<payload>
//   payload   "1d." + base64url(deflate-raw(JSON))     where the browser has CompressionStream
//             "1j." + base64url(UTF-8 JSON)            where it has not
// The first character is the format version. Everything after "#" stays in the browser: it is
// never sent to a server. ?ds= and ?qc= are as before, so the link also opens on an older page.
//
// What comes out of a link is data only. It is parsed as JSON and passed through clean(), which
// keeps known keys of the right type and range and nothing else. Names (genes, stains, outlines,
// clusterings) are only looked up in the open dataset's own lists; no address is ever read from
// a link, so a link cannot make the page load from anywhere else or run anything.
//
// Only unsigned typed arrays here (the de-identification check reads the signed names as sample IDs).
/* exported SodaViewLink */
(function (root) {
  "use strict";
  const FORMAT = 1;               // the payload's first character
  const VIEW_V = 1;               // view.v
  const BUDGET = 2000;            // characters of the whole link that every browser, chat app and mail client carries
  const MAX_PAYLOAD = 16000;      // longer payloads are not read at all
  const MAX_JSON = 512 * 1024;    // bytes of JSON a payload may unpack to
  const THIN_POINTS = 40;         // a freehand outline in a link that is too long keeps this many points
  const TABS = ["view", "analyze", "share"];
  const OLD_TAB = { clusters: "analyze", save: "share" };   // a link made before 0.10.0 named four tabs
  // more of the view, from other scripts (binpaint.js ...): { key, capture() -> value or null,
  // clean(value, helpers) -> value or null, apply(value or null) -> Promise of notes }. Each is one more
  // key of the view block: a file or link without it opens as before, and apply(null) puts it back to off.
  const EXT = [];

  // ---- base64url ------------------------------------------------------------------------------
  function toB64u(u8) {
    let s = "";
    for (let i = 0; i < u8.length; i += 0x8000) s += String.fromCharCode.apply(null, u8.subarray(i, i + 0x8000));
    return btoa(s).replace(/\+/g, "-").replace(/\//g, "_").replace(/=+$/, "");
  }
  function fromB64u(s) {
    if (!/^[A-Za-z0-9_-]*$/.test(s)) throw new Error("not base64url");
    const bin = atob(s.replace(/-/g, "+").replace(/_/g, "/") + "===".slice((s.length + 3) % 4));
    const u8 = new Uint8Array(bin.length);
    for (let i = 0; i < bin.length; i++) u8[i] = bin.charCodeAt(i);
    return u8;
  }

  // ---- deflate, where the browser has it ---------------------------------------------------------
  const canDeflate = () => typeof CompressionStream === "function" && typeof DecompressionStream === "function";
  async function pump(stream, u8, max) {
    const w = stream.writable.getWriter();
    // the write only finishes once the other side is read: wait for it at the end, and let a refusal surface there
    const sent = w.write(u8).then(() => w.close());
    sent.catch(() => {});
    const r = stream.readable.getReader(), parts = [];
    let n = 0;
    for (;;) {
      const { done, value } = await r.read();
      if (done) break;
      n += value.length;
      if (max && n > max) { try { await r.cancel(); } catch (e) { /* stopping is all that matters */ } throw new Error("too big"); }
      parts.push(value);
    }
    await sent;
    const out = new Uint8Array(n);
    let at = 0;
    for (const p of parts) { out.set(p, at); at += p.length; }
    return out;
  }

  // ---- encode / decode -----------------------------------------------------------------------------
  // opts.compress === false: the "1j." form whatever the browser can do (tests, and the fallback)
  async function encode(obj, opts) {
    const raw = new TextEncoder().encode(JSON.stringify(obj));
    const plain = `${FORMAT}j.${toB64u(raw)}`;
    if ((opts && opts.compress === false) || !canDeflate()) return plain;
    try {
      const packed = `${FORMAT}d.${toB64u(await pump(new CompressionStream("deflate-raw"), raw))}`;
      return packed.length < plain.length ? packed : plain;
    } catch (e) { return plain; }
  }
  // -> { ok: true, state } or { ok: false, why, message }. Never throws.
  async function decode(payload) {
    const bad = (why, message) => ({ ok: false, why, message });
    const DAMAGED = "This link's view could not be read (damaged or cut short): showing the default view.";
    try {
      if (typeof payload !== "string" || !payload) return bad("empty", DAMAGED);
      if (payload.length > MAX_PAYLOAD) return bad("toobig", "This link is too long to be a SODA view: showing the default view.");
      const m = /^(\d{1,3})([a-z])\.(.*)$/.exec(payload);
      if (!m) return bad("damaged", DAMAGED);
      if (Number(m[1]) > FORMAT) return bad("newer", "This link was made by a newer SODA than this page: showing the default view.");
      if (Number(m[1]) < 1) return bad("damaged", DAMAGED);
      let bytes = fromB64u(m[3]);
      if (m[2] === "d") {
        if (!canDeflate()) return bad("nodeflate", "This browser cannot unpack this link's view: showing the default view.");
        bytes = await pump(new DecompressionStream("deflate-raw"), bytes, MAX_JSON);
      } else if (m[2] !== "j") return bad("damaged", DAMAGED);
      if (bytes.length > MAX_JSON) return bad("toobig", "This link is too long to be a SODA view: showing the default view.");
      const j = JSON.parse(new TextDecoder("utf-8", { fatal: true }).decode(bytes));
      const state = clean(j);
      if (!state) return bad("damaged", DAMAGED);
      return { ok: true, state };
    } catch (e) { return bad("damaged", DAMAGED); }
  }

  // ---- clean(): known keys, right types, sane ranges; everything else is dropped ---------------------
  const isObj = (x) => !!x && typeof x === "object" && !Array.isArray(x);
  const UNSAFE_KEY = /^(__proto__|prototype|constructor)$/;
  // text that came from outside: no markup characters, no control characters, a length limit
  function str(x, max) {
    if (typeof x !== "string") return null;
    // eslint-disable-next-line no-control-regex
    const t = x.replace(/[\u0000-\u001f\u007f<>"'`]/g, "").trim().slice(0, max || 80);
    return t || null;
  }
  const num = (x, lo, hi) => (typeof x === "number" && Number.isFinite(x) ? Math.min(hi, Math.max(lo, x)) : null);
  const int = (x, lo, hi) => (Number.isInteger(x) && x >= lo && x <= hi ? x : null);
  const hex = (x) => (typeof x === "string" && /^#[0-9a-fA-F]{6}$/.test(x) ? x.toLowerCase() : null);
  const oneOf = (x, list) => (list.includes(x) ? x : null);
  const bool = (x) => (x === true || x === 1 ? true : x === false || x === 0 ? false : null);
  const list = (x, max, f) => (Array.isArray(x) ? x.slice(0, max).map(f).filter((v) => v !== null && v !== undefined) : null);
  const put = (o, k, v) => { if (v !== null && v !== undefined) o[k] = v; };
  function cam(x) {
    if (!Array.isArray(x) || x.length < 3) return null;
    const a = num(x[0], -1e9, 1e9), b = num(x[1], -1e9, 1e9), z = num(x[2], -40, 40);
    return a === null || b === null || z === null ? null : [a, b, z];
  }
  function range(x) {
    if (!Array.isArray(x) || x.length < 2) return null;
    const a = num(x[0], -1e12, 1e12), b = num(x[1], -1e12, 1e12);
    return a === null || b === null || b < a ? null : [a, b];
  }
  function gene(x) {
    if (typeof x === "string") x = { n: x };
    if (!isObj(x)) return null;
    const n = str(x.n, 80);
    if (!n) return null;
    const o = { n };
    put(o, "c", hex(x.c)); put(o, "r", range(x.r)); if (bool(x.on) === false) o.on = false;
    return o;
  }
  function cleanView(v) {
    if (!isObj(v)) return null;
    const o = { v: int(v.v, 1, 999) || VIEW_V };
    put(o, "mode", oneOf(v.mode, ["tissue", "sbs", "umap", "rna"]));
    put(o, "cam", cam(v.cam)); put(o, "ucam", cam(v.ucam)); put(o, "scam", cam(v.scam));
    if (Array.isArray(v.win) && v.win.length >= 2) { const w = num(v.win[0], 1, 20000), h = num(v.win[1], 1, 20000); if (w && h) o.win = [w, h]; }
    if (isObj(v.img)) { const g = {}; put(g, "on", bool(v.img.on)); put(g, "op", num(v.img.op, 0, 1)); o.img = g; }
    put(o, "ch", list(v.ch, 300, (c) => {
      if (!isObj(c)) return null;
      const n = str(c.n, 120); if (!n) return null;
      const e = { n, on: bool(c.on) === true };
      put(e, "c", hex(c.c));
      if (Array.isArray(c.w) && c.w.length >= 2) { const a = num(c.w[0], 0, 255), b = num(c.w[1], 0, 255); if (a !== null && b !== null && b >= a) e.w = [Math.round(a), Math.round(b)]; }
      // gamma and colour map (absent: 1 and the stain's own colour, as before they existed)
      put(e, "g", num(c.g, 0.2, 5)); put(e, "m", oneOf(c.m, ["solid", "grey", "viridis", "magma"]));
      return e;
    }));
    put(o, "genes", list(v.genes, 24, gene));
    put(o, "blend", oneOf(v.blend, ["normal", "heatmap"]));
    put(o, "palette", oneOf(v.palette, ["genes", "viridis"]));
    put(o, "gop", num(v.gop, 0, 1)); put(o, "ps", num(v.ps, 0.25, 4)); put(o, "tb", bool(v.tb));
    if (isObj(v.mol)) {
      const g = {};
      put(g, "genes", list(v.mol.genes, 24, gene)); put(g, "size", num(v.mol.size, 0.5, 4)); put(g, "op", num(v.mol.op, 0.1, 1));
      o.mol = g;
    }
    put(o, "ov", list(v.ov, 50, (s) => str(s, 120)));
    put(o, "colorBy", oneOf(v.colorBy, ["gene", "cluster", "celltype"]));
    put(o, "pick", int(v.pick, 0, 65534));
    put(o, "tab", oneOf(OLD_TAB[v.tab] || v.tab, TABS));
    for (const x of EXT) { if (o[x.key] !== undefined || UNSAFE_KEY.test(x.key)) continue; try { put(o, x.key, x.clean(v[x.key], HELP_FNS)); } catch (e) { /* that part is dropped */ } }
    if (isObj(v.rna)) {
      const r = {}, x = v.rna;
      put(r, "colorBy", str(x.colorBy, 80)); put(r, "genes", list(x.genes, 24, gene));
      put(r, "view", oneOf(x.view, ["normal", "density"])); put(r, "palette", oneOf(x.palette, ["genes", "viridis"]));
      put(r, "bin", oneOf(x.bin, [1, 2, 4, 8])); put(r, "ps", num(x.ps, 0.5, 6)); put(r, "labels", bool(x.labels)); put(r, "cam", cam(x.cam));
      o.rna = r;
    }
    return o;
  }
  const HELP_FNS = { isObj, str, num, int, hex, oneOf, bool, list, put };
  function cleanShape(s) {
    if (!isObj(s)) return null;
    const pts = list(s.points, 5000, (p) => {
      if (!Array.isArray(p) || p.length < 2) return null;
      const a = num(p[0], -1e9, 1e9), b = num(p[1], -1e9, 1e9);
      return a === null || b === null ? null : [a, b];
    }) || [];
    const o = { points: pts };
    put(o, "name", str(s.name, 60)); put(o, "color", hex(s.color)); put(o, "kind", oneOf(s.kind, ["circle", "free", "sel"]));
    if (Array.isArray(s.rows)) o.rows = s.rows.slice(0, 500000).filter((r) => Number.isInteger(r) && r >= 0 && r < 4294967295);
    put(o, "def", cleanDef(s.def));
    return o.rows || pts.length >= 3 || o.def ? o : null;
  }
  // How a combined shape was made (selections.js saveDef): outlines and set operations on them. A link carries this
  // instead of the shape's bins; selections.js checks it again as it reads it.
  const DEF_DEPTH = 6, DEF_NODES = 128, DEF_POINTS = 20000;
  function cleanDef(d) {
    let nodes = 0, points = 0;
    const pt = (p) => {
      if (!Array.isArray(p) || p.length < 2) return null;
      const a = num(p[0], -1e9, 1e9), b = num(p[1], -1e9, 1e9);
      return a === null || b === null ? null : [a, b];
    };
    const ring = (r) => {
      const l = list(r, 5000, pt);
      if (!l || l.length < 3 || (points += l.length) > DEF_POINTS) return null;
      return l;
    };
    const walk = (x, depth) => {
      if (!isObj(x) || depth > DEF_DEPTH || ++nodes > DEF_NODES) return null;
      if (Array.isArray(x.points)) { const r = ring(x.points); return r ? { points: r } : null; }
      if (Array.isArray(x.geo)) {
        const geo = x.geo.slice(0, 64).map((rings) => (Array.isArray(rings) ? rings.slice(0, 64).map(ring) : [null]));
        return geo.length && geo.every((rs) => rs.length && rs.every(Boolean)) ? { geo } : null;
      }
      const op = oneOf(x.op, ["union", "inter", "sub"]);
      if (op && Array.isArray(x.parts) && x.parts.length >= 2 && x.parts.length <= 64) {
        const parts = x.parts.map((y) => walk(y, depth + 1));
        return parts.every(Boolean) ? { op, parts } : null;
      }
      return null;
    };
    return d === undefined ? null : walk(d, 0);
  }
  function thinDef(d) {
    if (d.points) return { points: thin(d.points, THIN_POINTS) };
    if (d.geo) return { geo: d.geo.map((rings) => rings.map((r) => thin(r, THIN_POINTS))) };
    return { op: d.op, parts: d.parts.map(thinDef) };
  }
  const defPoints = (d) => (d.points ? d.points.length : d.geo ? d.geo.flat().reduce((n, r) => Math.max(n, r.length), 0) : Math.max(0, ...d.parts.map(defPoints)));
  // The whole state of a link (not of a file: Load JSON reads its file as it always did, and only
  // the `view` block of a file goes through cleanView).
  function clean(j) {
    if (!isObj(j)) return null;
    const o = {};
    put(o, "dataset", str(j.dataset, 200)); put(o, "layer", str(j.layer, 120));
    put(o, "seed", int(j.seed, 0, 2147483647)); put(o, "mode", oneOf(j.mode, ["bins", "shapes"]));
    if (j.resolution === "auto") o.resolution = "auto"; else put(o, "resolution", num(j.resolution, 0, 100));
    put(o, "sampling", oneOf(j.sampling, ["auto", "exact", "fast"]));     // Draw a region: every bin, or a sample above the threshold
    put(o, "shapes", list(j.shapes, 64, cleanShape));
    if (isObj(j.sectionClusters)) {
      const c = j.sectionClusters, s = {};
      put(s, "layer", str(c.layer, 120));
      if (typeof c.resolution === "number") put(s, "resolution", num(c.resolution, 0, 1000)); else put(s, "resolution", str(c.resolution, 80));
      s.names = {};
      if (isObj(c.names)) for (const rk of Object.keys(c.names).slice(0, 40)) {
        if (UNSAFE_KEY.test(rk) || rk.length > 80 || !isObj(c.names[rk])) continue;
        const per = {};
        for (const ck of Object.keys(c.names[rk]).slice(0, 2000)) { const t = /^\d{1,5}$/.test(ck) ? str(c.names[rk][ck], 60) : null; if (t) per[ck] = t; }
        if (Object.keys(per).length) s.names[rk] = per;
      }
      s.hidden = list(c.hidden, 2000, (x) => int(x, 0, 65534)) || [];
      // cell sets (cellsets.js): clusters merged into one set, per clustering. viewer.js cleans them again as it reads them.
      if (isObj(c.sets)) {
        const sets = {};
        for (const rk of Object.keys(c.sets).slice(0, 40)) {
          if (UNSAFE_KEY.test(rk) || rk.length > 80 || !Array.isArray(c.sets[rk])) continue;
          const l = c.sets[rk].slice(0, 500).map((x) => {
            const ids = isObj(x) ? list(x.ids, 2000, (k) => int(k, 0, 65534)) : null;
            return ids && ids.length >= 2 ? { ids, prev: str(x.prev, 60) } : null;
          }).filter(Boolean);
          if (l.length) sets[rk] = l;
        }
        if (Object.keys(sets).length) s.sets = sets;
      }
      if (s.layer && s.resolution !== undefined) o.sectionClusters = s;
    }
    if (isObj(j.qcFilter) && isObj(j.qcFilter.thresholds)) {
      const th = {};
      for (const k of Object.keys(j.qcFilter.thresholds).slice(0, 20)) {
        const t = j.qcFilter.thresholds[k];
        if (UNSAFE_KEY.test(k) || !/^\w{1,40}$/.test(k) || !isObj(t)) continue;
        const lo = num(t.min, 0, 1e12), hi = num(t.max, 0, 1e12);
        if (lo !== null || hi !== null) th[k] = { min: lo, max: hi };
      }
      o.qcFilter = { thresholds: th };
      put(o.qcFilter, "layer", str(j.qcFilter.layer, 120));
    }
    put(o, "view", cleanView(j.view));
    return o;
  }

  // ---- from the full state to what a link carries ----------------------------------------------------
  const r2 = (x) => Math.round(x * 100) / 100, r3 = (x) => Math.round(x * 1000) / 1000, r4 = (x) => Math.round(x * 10000) / 10000;
  // sessionState(false) minus what only a file needs and what is at its default
  function forLink(st) {
    const o = { dataset: st.dataset };
    if (st.layer) o.layer = st.layer;
    if (Number.isInteger(st.seed) && st.seed !== 42) o.seed = st.seed;
    if (st.mode === "shapes") o.mode = st.mode;
    if (st.resolution !== "auto" && Number.isFinite(st.resolution)) o.resolution = st.resolution;
    if (st.sampling === "exact" || st.sampling === "fast") o.sampling = st.sampling;
    if (Array.isArray(st.shapes) && st.shapes.length) o.shapes = st.shapes;
    const sc = st.sectionClusters;
    if (sc) o.sectionClusters = { layer: sc.layer, resolution: sc.resolution, names: sc.names || {}, hidden: sc.hidden || [], ...(sc.sets ? { sets: sc.sets } : {}) };
    if (st.qcFilter && st.qcFilter.thresholds) o.qcFilter = { layer: st.qcFilter.layer, thresholds: st.qcFilter.thresholds };
    if (st.view) o.view = st.view;
    return o;
  }
  function thin(pts, max) {
    if (pts.length <= max) return pts.map((p) => [Math.round(p[0]), Math.round(p[1])]);
    const out = [];
    for (let i = 0; i < max; i++) { const p = pts[Math.floor((i * pts.length) / max)]; out.push([Math.round(p[0]), Math.round(p[1])]); }
    return out;
  }
  // The tiers a link goes down when it is too long: each returns [state, what was left out]
  const TIERS = [
    (s) => [s, []],
    (s) => {
      const all = s.shapes || [], drawn = all.filter((x) => !x.rows);
      const left = [];
      if (drawn.length < all.length) left.push(all.length - drawn.length === 1 ? "a shape made from a selection" : `${all.length - drawn.length} shapes made from selections`);
      if (drawn.some((x) => (x.points || []).length > THIN_POINTS || (x.def && defPoints(x.def) > THIN_POINTS))) left.push("the fine detail of the drawn outlines");
      // a combined shape (no rows in a link: its recipe, def) is thinned like a drawn outline
      return [{ ...s, shapes: drawn.map((x) => ({ ...x, points: thin(x.points || [], THIN_POINTS), ...(x.def ? { def: thinDef(x.def) } : {}) })) }, left];
    },
    (s) => { const o = { ...s }; delete o.shapes; return [o, (s.shapes || []).length ? ["the drawn shapes"] : []]; },
    (s) => {
      const o = { ...s }; delete o.shapes;
      const named = !!(o.sectionClusters && Object.keys(o.sectionClusters.names || {}).length);
      if (named) o.sectionClusters = { ...o.sectionClusters, names: {} };
      return [o, named ? [...((s.shapes || []).length ? ["the drawn shapes"] : []), "the names given to clusters"] : []];
    },
  ];
  // -> { payload, left: [what was left out], tier, over: true when even the last tier is over the budget }
  async function fit(state, baseLength, budget) {
    const room = (budget || BUDGET) - baseLength - "#v=".length;
    let last = null;
    for (let t = 0; t < TIERS.length; t++) {
      const [s, left] = TIERS[t](state);
      if (t > 0 && !left.length) continue;                 // nothing (more) to leave out at this tier
      const payload = await encode(s);
      last = { payload, left, tier: t, over: payload.length > room };
      if (!last.over) return last;
    }
    return last;                                           // everything that can go has gone: a long link, and it says so
  }

  // ====================================================================================================
  // The page half: capture the view, put it back, the button. (Above: no DOM, also run by node in
  // tests/viewlink_test.js.)
  // ====================================================================================================
  let H = null;                                    // what viewer.js hands over (init)
  const $ = (id) => document.getElementById(id);
  // image px to 0.01 and zoom to 0.0001: under a tenth of a screen pixel anywhere in the window, at any zoom SODA allows
  const camOf = (vs) => (vs && Array.isArray(vs.target) && Number.isFinite(vs.zoom) ? [r2(vs.target[0]), r2(vs.target[1]), r4(vs.zoom)] : null);
  const ci = (a, b) => String(a).toLowerCase() === String(b).toLowerCase();
  // the dataset's own entry for a name from outside: exact first, then ignoring case
  const lookup = (items, name, key) => items.find((x) => x[key] === name) || items.find((x) => ci(x[key], name)) || null;

  // ---- capture ----------------------------------------------------------------------------------------
  function capture() {
    if (!H) return null;
    const S = H.state();
    if (!S.m) return null;
    const v = { v: VIEW_V }, el = $("viewport");
    if (el && el.clientWidth && el.clientHeight) v.win = [el.clientWidth, el.clientHeight];
    const lay = window.SodaLayout && window.SodaLayout.state ? window.SodaLayout.state() : null;
    if (lay && TABS.includes(lay.tab)) v.tab = lay.tab;
    if (S.mode === "rna") {
      const R = S.rna; v.mode = "rna";
      if (!R) return v;
      const r = { view: R.view, palette: R.palette, bin: R.bin, ps: R.pointSize, labels: !!R.showLabels };
      if (R.colorBy) r.colorBy = R.colorBy;
      r.genes = R.genes.map((g) => ({ n: g.name, c: g.hex, r: [r3(g.min), r3(g.max)], ...(g.visible === false ? { on: false } : {}) }));
      if (camOf(R.vs)) r.cam = [r4(R.vs.target[0]), r4(R.vs.target[1]), r4(R.vs.zoom)];   // UMAP units are small: four decimals
      v.rna = r;
      return v;
    }
    const U = S.umap, mode = U && ["tissue", "sbs", "umap"].includes(U.view) ? U.view : "tissue";
    v.mode = mode;
    const tissue = camOf(mode === "umap" ? U.tissueVs : S.viewState);
    if (tissue) v.cam = tissue;
    if (mode === "umap" && camOf(S.viewState)) v.ucam = camOf(S.viewState);
    if (mode === "sbs" && camOf(U.sideVs)) v.scam = camOf(U.sideVs);
    v.img = { on: !!S.imgOn, op: r3(S.imgOpacity) };
    // stains: those on, and those whose colour or window was changed (an H&E image has no list)
    if (!$("chanList").hidden) {
      v.ch = [];
      for (const ch of S.channels) {
        const d = S.m.image.channels[ch.index] || {}, w = d.window || [0, 255];
        const col = String(ch.hex).toLowerCase() !== String(d.color || "").toLowerCase(), win = ch.lo !== Math.round(w[0]) || ch.hi !== Math.round(w[1]);
        const gm = ch.gamma && ch.gamma !== 1 ? r3(ch.gamma) : null, cm = ch.cmap && ch.cmap !== "solid" ? ch.cmap : null;
        if (!ch.visible && !col && !win && !gm && !cm) continue;
        v.ch.push({ n: ch.name, on: !!ch.visible, ...(col ? { c: ch.hex } : {}), ...(win ? { w: [ch.lo, ch.hi] } : {}), ...(gm ? { g: gm } : {}), ...(cm ? { m: cm } : {}) });
      }
    }
    if (S.exprKey) {
      v.genes = S.genes.map((g) => ({ n: g.name, c: g.hex, r: [r3(g.min), r3(g.max)], ...(g.visible === false ? { on: false } : {}) }));
      v.blend = S.geneBlend; v.palette = S.genePalette; v.gop = r3(S.geneOpacity); v.ps = r3(S.pointScale);
      if (S.tissue && S.tissue.on) v.tb = true;
      if (U && U.key) { v.colorBy = U.colorBy; if (U.pick !== null && U.pick !== undefined && U.cl) v.pick = U.pick; }
    }
    if (S.mol && S.mol.index) v.mol = { genes: S.mol.genes.map((g) => ({ n: g.name, c: g.hex, ...(g.visible ? {} : { on: false }) })), size: r3(S.mol.size), op: r3(S.mol.opacity) };
    if (S.overlays && S.overlays.length) v.ov = S.overlays.filter((o) => o.on).map((o) => o.key);
    for (const x of EXT) { try { const c = x.capture(); if (c !== null && c !== undefined) v[x.key] = c; } catch (e) { /* left out */ } }
    return v;
  }

  // ---- apply ------------------------------------------------------------------------------------------
  // set a control and run its own handler, so the page's state and what the control shows stay one thing
  async function fire(id, value, prop) {
    const el = $(id);
    if (!el || el.disabled) return false;
    if (prop === "checked") { if (el.checked === value) return true; el.checked = value; }
    else {
      if (el.tagName === "SELECT" && ![...el.options].some((o) => o.value === String(value) && !o.disabled && !o.hidden)) return false;
      if (el.value === String(value)) return true;
      el.value = String(value);
    }
    const f = el.oninput || el.onchange;
    if (f) await f.call(el, { target: el });
    return true;
  }
  // the camera for this window: the same area in view as when the link was made, whatever the window's size
  function vsFrom(c, win, fit) {
    const el = $("viewport");
    let zoom = c[2];
    if (win && el && el.clientWidth && el.clientHeight) zoom += Math.log2(Math.min(el.clientWidth / win[0], el.clientHeight / win[1]));
    const f = fit || { minZoom: zoom - 2, maxZoom: zoom + 10 };
    return { target: [c[0], c[1], 0], zoom, minZoom: Math.min(f.minZoom, zoom), maxZoom: Math.max(f.maxZoom, zoom) };
  }
  // -> notes for the user (names that are not in this dataset, and the like)
  async function apply(view, whole) {
    const notes = [];
    view = cleanView(view);
    if (!H || !view) return notes;
    const S = H.state(), seq = H.seq(), stale = () => seq !== H.seq();
    if (!S.m) return notes;
    if (S.mode === "rna") { await applyRna(S, view, notes, stale); return notes; }
    if (view.mode === "rna") return notes;                   // a single-cell view on a slide: nothing to put back
    const missing = [];

    // the image and its stains
    if (view.img) {
      // Image opacity is the image's one control (0 = off). A link made while the "Image off" button was there says
      // on: false and keeps the opacity it had: the image is off, as it was then.
      const op = view.img.on === false ? 0 : view.img.op !== undefined ? view.img.op : view.img.on === true && !S.imgOn ? 1 : undefined;
      if (op !== undefined) { await fire("imgOpacity", op); S.imgOpacity = op; S.imgOn = op > 0; }      // the slider's steps are coarser than the value
    }
    if (view.ch && !$("chanList").hidden) {
      const want = new Map();
      for (const c of view.ch) { const ch = lookup(S.channels, c.n, "name"); if (ch) want.set(ch, c); else missing.push(c.n); }
      for (const ch of S.channels) {
        const c = want.get(ch), d = S.m.image.channels[ch.index] || {}, w = d.window || [0, 255];
        ch.visible = !!(c && c.on);
        ch.hex = (c && c.c) || d.color || ch.hex; ch.color = H.hexToRgb01(ch.hex);
        ch.lo = c && c.w ? c.w[0] : Math.round(w[0]); ch.hi = c && c.w ? c.w[1] : Math.round(w[1]);
        ch.gamma = (c && c.g) || 1; ch.cmap = (c && c.m) || "solid";
      }
      H.buildChannelUI();
    }

    // genes of the layer on show (sessionApply has already switched to the link's layer)
    if (S.exprKey) {
      if (view.blend && !$("rowGeneBlend").hidden) await fire("geneBlend", view.blend);
      if (view.palette) await fire("genePalette", view.palette);
      if (view.gop !== undefined) { await fire("geneOpacity", view.gop); S.geneOpacity = view.gop; }
      if (view.ps !== undefined) { await fire("pointScale", view.ps); S.pointScale = view.ps; }
      if (view.tb !== undefined) await fire("tissueBins", !!view.tb, "checked");
      if (stale()) return notes;
      if (view.genes) {
        S.genes = []; S.heat = null;
        for (const g of view.genes) {
          const hit = lookup(S.geneIndex, g.n, "n");
          if (!hit) { missing.push(g.n); continue; }
          await H.addGene(hit.n);
          if (stale()) return notes;
          const got = S.genes.find((x) => x.name === hit.n);
          if (!got) continue;
          if (g.c) got.hex = g.c;
          if (g.r) { got.min = g.r[0]; got.max = Math.max(g.r[1], g.r[0] + 1e-6); }
          if (g.on === false) got.visible = false;
          H.recolorGene(got);
        }
        S.heat = null;
        H.buildGeneUI();
      }
    }
    if (stale()) return notes;

    // molecules
    if (view.mol && S.mol && S.mol.index) {
      if (view.mol.size !== undefined) { await fire("molSize", view.mol.size); S.mol.size = view.mol.size; }
      if (view.mol.op !== undefined) { await fire("molOpacity", view.mol.op); S.mol.opacity = view.mol.op; }
      if (view.mol.genes) {
        S.mol.genes = [];
        for (const g of view.mol.genes) {
          const hit = lookup(S.mol.index.genes, g.n, "n");
          if (!hit) { missing.push(g.n); continue; }
          await H.addMolGene(hit.n);
          if (stale()) return notes;
          const got = S.mol.genes.find((x) => x.name === hit.n);
          if (!got) continue;
          if (g.c) { got.hex = g.c; got.rgb = H.hexToRgb255(g.c); }
          if (g.on === false) got.visible = false;
        }
        H.buildMolUI();
      }
    }

    // outlines: exactly the ones in the view
    if (view.ov && S.overlays.length) {
      const boxes = $("overlayList").querySelectorAll("input[type=checkbox]");
      for (let i = 0; i < S.overlays.length; i++) {
        const on = view.ov.some((k) => k === S.overlays[i].key), box = boxes[i];
        if (!box || box.checked === on) continue;
        box.checked = on;
        // (an outline that will not load unticks itself and says so on the page: noted here too)
        try { await box.onchange({ target: box }); if (on && !box.checked) throw new Error("not loaded"); }
        catch (e) { notes.push(`The outline ${S.overlays[i].label} could not be loaded.`); }
        if (stale()) return notes;
      }
    }

    // more of the view from other scripts (extend): each one, also when the view does not have it (back to off)
    for (const x of EXT) {
      try { const n = await x.apply(view[x.key] === undefined ? null : view[x.key]); if (n && n.length) notes.push(...n); }
      catch (e) { notes.push("Part of this view could not be put back."); }
      if (stale()) return notes;
    }

    // the section's clustering, colour-by, the highlighted cluster, the mode
    const U = S.umap, sc = whole && whole.sectionClusters, hasUmap = !!(U && U.key && !$("viewSwitch").hidden);
    let mode = view.mode || "tissue";
    if (U && U.key) {
      const resOpt = sc && sc.layer === U.key ? [...$("clusterRes").options].find((o) => o.value === String(sc.resolution)) : null;
      const needCl = (view.colorBy && view.colorBy !== "gene") || mode !== "tissue" || view.pick !== undefined;
      try {
        if (resOpt && (needCl || U.cl)) { $("clusterRes").value = resOpt.value; await $("clusterRes").onchange({ target: $("clusterRes") }); }
        if (stale()) return notes;
        if (view.colorBy && !(await fire("colorBy", view.colorBy)) && view.colorBy !== "gene") notes.push("This dataset cannot be coloured the way the link was: coloured by genes.");
        if (stale()) return notes;
        if (sc && sc.layer === U.key && U.cl && String(U.clRes) === String(sc.resolution)) {
          const n = H.clusterN();
          U.hidden = new Set((sc.hidden || []).filter((c) => c < n));
          if (view.pick !== undefined && view.pick < n && !U.hidden.has(view.pick)) U.pick = view.pick;
          H.umapDirty();
        }
        const noEmb = [...$("viewSwitch").querySelectorAll(".seg")].some((e) => e.hidden);
        if (mode !== "tissue" && (!hasUmap || noEmb)) { notes.push("This dataset has no UMAP: showing the tissue."); mode = "tissue"; }
        if (mode !== "tissue") { await H.umapLoad(); if (!U.cl && H.hasClusterings()) await H.umapSetRes(U.res); }
        if (stale() || S.umap !== U) return notes;
        if (mode !== "tissue" && !(U.uxy && U.nEmb)) mode = "tissue";
      } catch (e) { notes.push("The UMAP of this layer could not be loaded: showing the tissue."); mode = "tissue"; }
      if (U.anim) U.anim = null;
      if (U.view !== mode) { U.view = mode; H.umapApplyViews(); }
    } else if (mode !== "tissue") { notes.push("This dataset has no UMAP: showing the tissue."); mode = "tissue"; }

    // the cameras, last: everything above may have re-framed the picture
    const fitT = H.imageFitVs(), tissueVs = view.cam ? vsFrom(view.cam, view.win, fitT) : fitT;
    if (mode === "umap" && U) {
      U.tissueVs = tissueVs;
      H.setCamera(view.ucam ? vsFrom(view.ucam, view.win, H.umapFitVs()) : H.umapFitVs());
    } else {
      if (mode === "sbs" && U) {
        const [w2, h2] = H.umapPaneSize("uside");
        U.sideVs = view.scam ? vsFrom(view.scam, view.win, H.umapFitVs(w2, h2)) : H.umapFitVs(w2, h2);
      }
      H.setCamera(tissueVs);
    }
    H.syncUmapUI();
    H.render();
    if (missing.length) notes.push(`Not in this dataset: ${missing.slice(0, 6).join(", ")}${missing.length > 6 ? ` and ${missing.length - 6} more` : ""}.`);
    return notes;
  }
  async function applyRna(S, view, notes, stale) {
    const R = S.rna, r = view.rna;
    if (!R || !r) return;
    const missing = [];
    if (r.colorBy && r.colorBy !== R.colorBy) { const o = lookup(R.obs, r.colorBy, "name"); if (o) await H.setColorBy(o.name); }
    if (stale() || S.rna !== R) return;
    if (r.palette) await fire("rnaPalette", r.palette);
    if (r.view) await fire("rnaView", r.view);
    if (r.bin) await fire("rnaBin", r.bin);
    if (r.ps !== undefined) { await fire("rnaPointSize", r.ps); R.pointSize = r.ps; }
    if (r.labels !== undefined) await fire("rnaLabels", !!r.labels, "checked");
    if (r.genes) {
      R.genes = []; R.heat = null; R.grid = "";
      for (const g of r.genes) {
        const hit = lookup(R.index, g.n, "n");
        if (!hit) { missing.push(g.n); continue; }
        await H.addRnaGene(hit.n);
        if (stale() || S.rna !== R) return;
        const got = R.genes.find((x) => x.name === hit.n);
        if (!got) continue;
        if (g.c) got.hex = g.c;
        if (g.r) { got.min = g.r[0]; got.max = Math.max(g.r[1], g.r[0] + 1e-6); }
        if (g.on === false) got.visible = false;
        H.recolorRnaGene(got);
      }
      R.heat = null; R.grid = "";
      H.buildRnaGeneUI();
    }
    H.renderRna();                                 // lays the panels out and fits them: the camera goes on top of that
    if (r.cam && R.vs) { H.setCamera(vsFrom(r.cam, view.win, R.vs)); H.renderRna(); }
    if (missing.length) notes.push(`Not in this dataset: ${missing.slice(0, 6).join(", ")}${missing.length > 6 ? ` and ${missing.length - 6} more` : ""}.`);
  }

  // ---- the link ---------------------------------------------------------------------------------------
  // what else of the address a link keeps: where the data is (?data= a shared folder, ?static= a local test
  // folder) and the quality filter. ?diag= and ?decode= are for looking at a problem, not for sharing.
  const KEEP_PARAMS = ["data", "static", "qc"];
  function baseUrl(dsId) {
    const now = new URLSearchParams(location.search), q = new URLSearchParams();
    q.set("ds", dsId);
    for (const k of KEEP_PARAMS) if (now.get(k)) q.set(k, now.get(k));
    return `${location.origin}${location.pathname}?${q.toString()}`;
  }
  // -> { url, length, left, local, over } for the view on show
  async function make() {
    const S = H.state();
    const base = baseUrl(S.dsId);
    const r = await fit(forLink(H.sessionState(false)), base.length, BUDGET);
    const url = `${base}#v=${r.payload}`;
    return { url, length: url.length, left: r.left, over: r.over, tier: r.tier, local: !!new URLSearchParams(location.search).get("static") };
  }
  function describe(r) {
    let t = `Link to this view copied (${r.length.toLocaleString()} characters).`;
    if (r.left.length) t += ` Left out, too big for a link: ${r.left.join(", ")}. Save JSON keeps everything.`;
    if (r.over) t += " It is still a long link: some chat apps may cut it short.";
    if (r.local) t += " It points at a test folder on this computer, so it only opens here.";
    return t;
  }
  async function copy() {
    const msg = $("saveMsg");
    if (!H || !H.state().m) { if (msg) msg.textContent = "Open a dataset first."; return null; }
    let r;
    try { r = await make(); } catch (e) { if (msg) msg.textContent = "Could not make a link to this view."; return null; }
    let ok = false;
    try { await navigator.clipboard.writeText(r.url); ok = true; }
    catch (e) {
      // no clipboard permission (an http page, an older browser): the old way
      try {
        const ta = document.createElement("textarea");
        ta.value = r.url; ta.setAttribute("readonly", ""); ta.style.cssText = "position:fixed;left:0;top:0;opacity:0;font-size:16px";
        document.body.appendChild(ta); ta.select(); ok = document.execCommand("copy"); ta.remove();
      } catch (e2) { ok = false; }
    }
    let box = $("saveLinkBox");
    if (!ok && msg) {
      // the address bar does not hold the view: show the link so it can be copied by hand
      if (!box) {
        box = document.createElement("input");
        box.id = "saveLinkBox"; box.type = "text"; box.readOnly = true; box.className = "save-link-box";
        box.setAttribute("aria-label", "Link to this view");
        msg.after(box);
      }
      box.value = r.url; box.hidden = false; box.focus(); box.select();
    } else if (box) box.hidden = true;
    if (msg) msg.textContent = ok ? describe(r) : describe(r).replace(/^Link to this view copied/, "Could not copy by itself: select the link in the box below and copy it");
    return r;
  }

  // ---- a link being opened ----------------------------------------------------------------------------
  const hashPayload = () => { const m = /(?:^#|&)v=([^&]*)/.exec(location.hash || ""); return m ? m[1] : null; };
  let pending = typeof location !== "undefined" ? hashPayload() : null;      // waits for its dataset to open (after the password, if any)
  let lastNotes = null;
  function say(text, ms) {
    if (!H) return;
    H.status(text);
    setTimeout(() => { const st = $("status"); if (st && st.textContent === text) H.status(null); }, ms || 7000);
    const msg = $("saveMsg"); if (msg) msg.textContent = text;
  }
  function dropHash() {
    pending = null;
    try { history.replaceState(null, "", location.pathname + location.search); } catch (e) { /* the address keeps its tail: harmless */ }
  }
  let busy = false;
  async function opened(id) {
    if (!pending || !H || busy) return;
    const S = H.state();
    if (!S.m || S.dsId !== id) return;             // the open failed or was overtaken: the view waits for the next one
    busy = true;
    const raw = pending;
    try {
      const d = await decode(raw);
      if (!d.ok) { dropHash(); lastNotes = [d.message]; say(d.message); return; }
      if (d.state.dataset && d.state.dataset !== S.dsId) {
        dropHash(); lastNotes = ["This link's view belongs to another dataset: showing the default view."]; say(lastNotes[0]); return;
      }
      const lay = window.SodaLayout, tok = lay ? lay.opening(true) : 0;
      let notes = [];
      try { notes = (await H.sessionApply(d.state, false)) || []; }
      catch (e) { notes = ["Part of this link's view could not be put back."]; }
      finally { if (lay) lay.opening(false, tok); }
      if (H.state().dsId !== id) return;           // another dataset was picked while the view went back
      dropHash();
      const tab = d.state.view && d.state.view.tab;
      if (tab && lay && lay.state().tabs.includes(tab)) lay.select(tab);
      lastNotes = notes;
      say(["View restored from the link.", ...notes].join(" "), notes.length ? 12000 : 4000);
    } finally { busy = false; }
  }

  // The Share tab's "Copy link" button (layout.js builds it, with no click of its own) becomes "Copy link to this
  // view": its words, its click and the sentence about it in the section's "?" note.
  const HELP = "Copy link to this view copies an address that opens this dataset as it looks now: position and zoom, layer, genes, stains, colours, clusters, the quality filter and the open tab. Drawn shapes too big for a link are left out, and it says so. Save JSON keeps the same view in a file.";
  function takeButton() {
    const b = $("btnCopyLink");
    if (!b) return;
    b.textContent = "Copy link to this view";
    b.title = "Copy an address that opens this dataset as it looks now: position and zoom, layer, genes, stains, colours, clusters and the quality filter";
    b.onclick = () => { copy(); };
    const p = document.querySelector("#secHelp-panelSave p");
    if (p && p.dataset.text && !p.dataset.text.includes(HELP)) {
      const old = /\s*Copy link copies this page's address[^.]*\.?/;
      p.dataset.text = old.test(p.dataset.text) ? p.dataset.text.replace(old, " " + HELP) : `${p.dataset.text} ${HELP}`;
      p.textContent = p.dataset.text;
    }
  }

  const api = {
    FORMAT, VIEW_V, BUDGET, MAX_PAYLOAD, MAX_JSON,
    encode, decode, clean, cleanView, forLink, fit,
    // viewer.js
    init(hooks) {
      H = hooks;
      takeButton();
      window.addEventListener("soda:opened", (e) => { opened(e.detail && e.detail.id); });
      // a link pasted into a page that is already open
      window.addEventListener("hashchange", () => { const p = hashPayload(); if (p) { pending = p; opened(H.state().dsId); } });
    },
    capture, apply, make, copy,
    extend(x) { if (x && typeof x.key === "string" && /^\w{1,20}$/.test(x.key) && !EXT.some((e) => e.key === x.key)) EXT.push(x); },
    // tests: is a link still waiting, and what was said when the last one was applied
    pending: () => !!pending, notes: () => lastNotes,
  };
  root.SodaViewLink = api;
  if (typeof module !== "undefined" && module.exports) module.exports = api;
})(typeof window !== "undefined" ? window : globalThis);
