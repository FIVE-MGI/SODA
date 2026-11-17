// Cell outlines: the segmentation boundaries of a cells layer (Xenium, MERSCOPE, CosMx, Stereo-seq
// cell bins), drawn instead of the dots once the cells are large enough on screen.
//
// The export holds them as spatial tiles of quantised polygons (cellshapes.py, format
// "poly-u16-tiles", docs/DATA_FORMAT.md): expr/<key>/shapes/<set>/index.json and <tx>_<ty>.bin,
// set = "cell" or "nucleus". A tile lists, for each polygon, the row of positions.bin it belongs
// to, so an outline is coloured, lit, greyed and picked by row exactly as its dot is.
//
//   - Level of detail: nothing is fetched or drawn while the median cell is under ~5 px across;
//     between 5 and 9 px the outlines fade in and the dots fade out; above, outlines only.
//   - Only the tiles under the view (grown by the reach of the largest cell) are fetched; a tile's
//     polygons are tessellated once, when it arrives. Too many polygons in view (a phone, a very
//     large window): the view stays on dots and the panel says so.
//   - Colour never rebuilds geometry: the fill and line colours are accessors that look a row up
//     in one RGBA array per layer, re-read when the colouring's key changes.
//   - A cell listed without an outline (count 0 in the tile) keeps its dot.
// viewer.js owns the colouring (genes, clusters, selection, quality filter) and hands it over
// through ctx.colors(); this file owns the tiles, the layers and the two controls.
// Only unsigned typed arrays here: the de-identification check reads the signed names as IDs.
/* exported SodaOutlines */
(function () {
  "use strict";
  const FORMAT = "poly-u16-tiles", VERSION = 1;
  const KEY = "soda.cellShapes";
  const LINE = [226, 232, 244], DIM = [150, 158, 176], GREY = [88, 88, 88], NUCLEUS = [150, 205, 255];

  function create(ctx) {
    const $ = ctx.$, cfg = ctx.cfg || {}, D = window.deck;
    // the median cell's width on screen, in pixels, at which outlines start and are complete
    const T0 = cfg.outlineFromPx || 5, T1 = Math.max((cfg.outlineFullPx || 9), T0 + 0.5);
    const BUDGET = cfg.outlineBudget || (ctx.TOUCH ? 20000 : 120000);        // polygons in view
    const KEEP_VERTS = ctx.TOUCH ? 4e5 : 3e6;                                // decoded vertices kept for a quick return
    // deck builds a tile's triangles and line segments the first time it is drawn (some tens of ms for
    // a full tile): tiles that arrive together go in over a few frames, about this many polygons a frame
    const NEW_PER_FRAME = 700;
    const PARALLEL = 6;
    let mode = "outlines", nuclei = false;
    try { if (localStorage.getItem(KEY) === "dots") mode = "dots"; } catch (e) { /* private mode */ }
    let O = blank();

    function blank() {
      return { ds: null, key: null, e: null, src: null, sets: [], seq: 0, alpha: 0, dot: 1, sig: "", note: "",
               col: null, byLayer: new Map(), drawn: new Set(), waiting: 0, stats: { fetched: 0, bytes: 0, tiles: 0, decodeMs: 0, firstMs: null, t0: 0 } };
    }
    const setOf = (key) => O.sets.find((s) => s.key === key) || null;
    const primary = () => O.sets[0] || null;
    const available = () => !!(O.sets.length && O.src && O.src.json && O.src.buf);

    function reset() { O = blank(); syncUI(); }

    // Called when the expression layer changes (and so on every dataset open).
    function onLayer(dsId, key, e, src) {
      O = blank();
      O.ds = dsId; O.key = key; O.e = e; O.src = src;
      const sh = e && e.kind === "cells" ? e.shapes : null;
      if (sh && sh.format === FORMAT && sh.version === VERSION && Array.isArray(sh.sets) && src && src.json && src.buf) {
        O.sets = sh.sets.map((m) => ({ key: m.key, meta: m, index: null, busy: null, failed: false, chunks: new Map(),
                                       tiles: new Map(), loading: new Set(), bad: new Set() }));
      }
      if (!setOf("nucleus") || primary().key === "nucleus") nuclei = false;
      syncUI();
    }

    function syncUI() {
      const row = $("rowCellShapes"), nrow = $("rowCellNuclei"), note = $("cellShapesNote");
      if (!row) return;
      const on = available();
      row.hidden = !on;
      $("cellShapes").value = mode;
      const hasNuc = on && !!setOf("nucleus") && primary().key !== "nucleus";
      if (nrow) { nrow.hidden = !hasNuc || mode !== "outlines"; $("cellNuclei").checked = nuclei; }
      if (note) { note.hidden = !on || !O.note; note.textContent = O.note; }
    }
    function say(t) {
      if (O.note === t) return;
      O.note = t;
      const note = $("cellShapesNote");
      if (note) { note.hidden = !available() || !t; note.textContent = t; }
    }

    // ---- tiles ---------------------------------------------------------------------------------
    function loadIndex(S) {
      if (S.index || S.busy || S.failed) return;
      const seq = O.seq, own = O;
      S.busy = O.src.json(`expr/${O.key}/shapes/${S.key}/index.json`).then((ix) => {
        if (own !== O || seq !== O.seq) return;
        if (!ix || ix.format !== FORMAT || ix.version !== VERSION || !Array.isArray(ix.chunks)) throw new Error("not a poly-u16-tiles index");
        for (const c of ix.chunks) S.chunks.set(`${c.tx}_${c.ty}`, c);
        S.index = ix; S.busy = null;
        ctx.requestRender();
      }).catch((err) => {
        if (own !== O) return;
        S.failed = true; S.busy = null;
        ctx.status(`Cell outlines could not be loaded (${err.message}): the cells stay dots.`, true);
        ctx.requestRender();
      });
    }

    // One tile file -> closed rings in one Float32Array, ready for deck's binary polygon / path layers.
    function decode(buf, step) {
      const dv = new DataView(buf), n = dv.getUint32(0, true), nv = dv.getUint32(4, true);
      if (buf.byteLength !== 16 + 6 * n + 4 * nv) throw new Error(`${buf.byteLength} bytes for ${n} polygons, ${nv} vertices`);
      const x0 = dv.getFloat32(8, true), y0 = dv.getFloat32(12, true);
      const rows = new Uint32Array(buf, 16, n), cnt = new Uint16Array(buf, 16 + 4 * n, n), q = new Uint16Array(buf, 16 + 6 * n, 2 * nv);
      let np = 0;
      for (let i = 0; i < n; i++) if (cnt[i]) np++;
      const ringRow = new Uint32Array(np), start = new Uint32Array(np + 1), pos = new Float32Array(2 * (nv + np));
      const dotRows = new Uint32Array(n - np);
      let p = 0, v = 0, w = 0, d = 0;
      for (let i = 0; i < n; i++) {
        const c = cnt[i];
        if (!c) { dotRows[d++] = rows[i]; continue; }
        ringRow[p] = rows[i]; start[p] = w;
        const w0 = w;
        for (let k = 0; k < c; k++, v++, w++) { pos[2 * w] = x0 + q[2 * v] * step; pos[2 * w + 1] = y0 + q[2 * v + 1] * step; }
        pos[2 * w] = pos[2 * w0]; pos[2 * w + 1] = pos[2 * w0 + 1]; w++;       // closed: the first vertex again
        p++;
      }
      start[np] = w;
      const attributes = { getPolygon: { value: pos, size: 2 } };
      return { n: np, nVert: nv, ringRow, dotRows,
               fill: { length: np, startIndices: start, attributes },
               path: { length: np, startIndices: start, attributes: { getPath: { value: pos, size: 2 } } } };
    }

    function fetchTile(S, name) {
      const seq = O.seq, own = O, c = S.chunks.get(name);
      S.loading.add(name);
      if (!O.stats.t0) O.stats.t0 = performance.now();
      O.src.buf(`expr/${O.key}/shapes/${S.key}/${name}.bin`).then((buf) => {
        if (own !== O || seq !== O.seq) return;
        const t0 = performance.now();
        const T = decode(buf, S.index.step);
        T.name = name; T.used = performance.now();
        O.stats.decodeMs += T.used - t0; O.stats.fetched++; O.stats.bytes += buf.byteLength;
        if (c && T.n + T.dotRows.length !== c.n) throw new Error(`tile ${name}: ${T.n + T.dotRows.length} entries, the index says ${c.n}`);
        S.tiles.set(name, T); S.loading.delete(name);
        if (O.stats.firstMs === null) O.stats.firstMs = performance.now() - O.stats.t0;
        ctx.requestRender();
      }).catch((err) => {
        if (own !== O) return;
        S.loading.delete(name); S.bad.add(name);
        if (S.bad.size === 1) ctx.status(`Some cell outlines could not be loaded (${err.message}): those cells stay dots.`, true);
        ctx.requestRender();
      });
    }

    // the tiles of a set under the view: [names nearest the centre first, polygons in them]
    function tilesInView(S, box) {
      const ix = S.index, cp = ix.chunkPx, m = ix.overhangPx || 0;
      const tx0 = Math.floor((box[0] - m) / cp), tx1 = Math.floor((box[2] + m) / cp);
      const ty0 = Math.floor((box[1] - m) / cp), ty1 = Math.floor((box[3] + m) / cp);
      const out = [];
      let n = 0;
      if ((tx1 - tx0 + 1) * (ty1 - ty0 + 1) > 4 * S.chunks.size + 64) {            // a view far larger than the section
        for (const c of S.chunks.values()) if (c.tx >= tx0 && c.tx <= tx1 && c.ty >= ty0 && c.ty <= ty1) { out.push(c); n += c.n; }
      } else {
        for (let ty = ty0; ty <= ty1; ty++) for (let tx = tx0; tx <= tx1; tx++) {
          const c = S.chunks.get(`${tx}_${ty}`);
          if (c) { out.push(c); n += c.n; }
        }
      }
      const cx = (box[0] + box[2]) / 2 / cp - 0.5, cy = (box[1] + box[3]) / 2 / cp - 0.5;
      out.sort((a, b) => ((a.tx - cx) ** 2 + (a.ty - cy) ** 2) - ((b.tx - cx) ** 2 + (b.ty - cy) ** 2));
      return [out.map((c) => `${c.tx}_${c.ty}`), n];
    }

    // ---- what the view asks for -----------------------------------------------------------------
    // { alpha, box, cellPx } for the view; alpha 0 = dots only
    function plan(vs) {
      if (mode !== "outlines" || !available() || !vs || !Number.isFinite(vs.zoom) || !ctx.tissueView()) return null;
      const e = O.e, scale = Math.pow(2, vs.zoom);
      const cellPx = (e.spotDiameterPx || 8) * scale;
      if (cellPx < T0 - 1) return { alpha: 0, cellPx, early: false };
      const t = Math.max(0, Math.min(1, (cellPx - T0) / (T1 - T0)));
      const [w, h] = ctx.viewSize();
      const hw = (w || 1000) / scale / 2, hh = (h || 800) / scale / 2;
      return { alpha: t * t * (3 - 2 * t), cellPx, early: cellPx < T0,
               box: [vs.target[0] - hw, vs.target[1] - hh, vs.target[0] + hw, vs.target[1] + hh] };
    }

    // ---- colours: one RGBA per row for the fill, one for the line ---------------------------------
    // ctx.colors() may also hand over
    //   only  a mask per row (1 = keep): the gene colouring is drawn only where it is 1 (a region picked
    //         with a click: the density map is cut to that region's bins the same way)
    //   hl    { key, rgba } one RGBA per row (alpha 0 = none): a highlight drawn over the colouring, as the
    //         dots of Draw a region are drawn over the gene dots. Mixed in here ("over", straight alpha),
    //         so a highlight is a colour update of the layers there are, never another layer of polygons.
    //         The gene opacity then goes into the fill's alpha (the highlight is not faded by it, as its
    //         dots are not), and the layer's opacity is only the fade between dots and outlines.
    //   mix.alpha  one fill alpha for every coloured row, instead of one from its strength (Colour by a
    //         value: its dots are all drawn at that alpha, the colour alone says how much)
    function colours(N) {
      const C = ctx.colors();
      const hl = C && C.hl && C.hl.rgba && C.hl.rgba.length === 4 * N ? C.hl : null;
      const op = C ? C.opacity : 1;
      const key = C ? `${C.key}|${N}${hl ? `|hl${hl.key}|${op}` : ""}` : `none|${N}`;
      if (O.col && O.col.key === key) {
        // the opacity and dot size sliders change no colour: they are read every frame, not built in
        if (C) { O.col.radius = C.radius; if (!hl) O.col.opacity = op; }
        return O.col;
      }
      const fill = new Uint8Array(4 * N), line = new Uint8Array(4 * N);
      const mix = C && C.mix, cat = C && C.cat, lit = C && C.lit, keep = C && C.keep, sel = (C && C.sel) || [255, 214, 64];
      const only = C && C.only && C.only.length === N ? C.only : null, mixA = mix && mix.alpha > 0 ? mix.alpha : 0;
      const H = hl ? hl.rgba : null;
      let anyFill = false;
      for (let r = 0; r < N; r++) {
        const k = 4 * r;
        let c = null, fa = 0, lc = LINE, la = 150;
        const kept = !(keep && !keep[r]);
        if (!kept) { c = GREY; fa = 200; lc = GREY; la = 235; }                    // removed by the quality filter
        else {
          const isLit = lit ? !!lit(r) : false;
          let t = 1;
          if (cat) {
            const v = cat(r);
            if (v === null) { la = 0; }                                             // a hidden cluster: not drawn
            else if (v) c = v;
          } else if (mix && mix.t[r] >= 0 && (!only || only[r])) { c = [mix.rgb[3 * r], mix.rgb[3 * r + 1], mix.rgb[3 * r + 2]]; t = mix.t[r]; }
          if (lit && !isLit) { c = null; if (la) { lc = DIM; la = 55; } }            // something else is lit: a faint line
          else if (c) { fa = lit ? 235 : cat ? 175 : mixA || Math.round(70 + 170 * t); lc = [c[0] + (255 - c[0]) * 0.35, c[1] + (255 - c[1]) * 0.35, c[2] + (255 - c[2]) * 0.35]; la = 255; }
          else if (isLit) { c = sel; fa = 200; lc = sel; la = 255; }
        }
        if (!c || !fa) { c = null; fa = 0; }
        if (H) {
          if (fa) fa = Math.round(fa * op);
          const ha = H[k + 3];
          // the highlight over the colouring; a cell removed by the quality filter stays grey (its grey dot is drawn last)
          if (ha && kept) {
            const h = [H[k], H[k + 1], H[k + 2]], A = ha / 255, B = (fa / 255) * (1 - A), out = A + B;
            c = c ? [0, 1, 2].map((i) => Math.round((h[i] * A + c[i] * B) / out)) : h;
            fa = Math.round(255 * out);
            // a highlight that shows (not one dimmed because another shape is picked) also draws the line
            if (ha >= 64 && la) { lc = [h[0] + (255 - h[0]) * 0.35, h[1] + (255 - h[1]) * 0.35, h[2] + (255 - h[2]) * 0.35]; la = 255; }
          }
        }
        if (c && fa) { fill[k] = c[0]; fill[k + 1] = c[1]; fill[k + 2] = c[2]; fill[k + 3] = fa; anyFill = true; }
        line[k] = lc[0]; line[k + 1] = lc[1]; line[k + 2] = lc[2]; line[k + 3] = la;
      }
      O.col = { key, fill, line, anyFill, over: !!lit, pos: C ? C.pos : null, opacity: H ? 1 : op, radius: C ? C.radius : 4,
                pick: !!(C && C.pick) };
      return O.col;
    }

    // ---- the layers of one frame ------------------------------------------------------------------
    // Returns { layers, over } (over: draw them above the slide's dimming, because something is lit),
    // or null when the cells are dots. Sets what dotAlpha() answers until the next frame.
    function frame(vs) {
      const P = plan(vs), S = primary();
      O.byLayer.clear();
      const off = (note) => { O.alpha = 0; O.dot = 1; O.sig = sigOf(P, []); O.drawn = new Set(); O.waiting = 0; say(note || ""); return null; };
      if (!P) return off("");
      if (S.failed) return off("The outlines could not be loaded: cells are drawn as dots.");
      if (!P.box) return off("Zoom in to see the cell outlines.");
      loadIndex(S);
      const N = nuclei ? setOf("nucleus") : null;
      if (N) loadIndex(N);
      if (!S.index) return off("");
      const [names, nPoly] = tilesInView(S, P.box);
      if (nPoly > BUDGET) {
        return off(`Too many cells in view for outlines (${nPoly.toLocaleString()}; up to ${BUDGET.toLocaleString()} ${ctx.TOUCH ? "on this device" : "at once"}): zoom in.`);
      }
      const want = [[S, names]];
      if (N && N.index) want.push([N, tilesInView(N, P.box)[0]]);
      let inFlight = 0, missing = 0, failed = 0;
      for (const [X] of want) inFlight += X.loading.size;
      for (const [X, list] of want) for (const name of list) {
        if (X.tiles.has(name)) continue;
        if (X.bad.has(name)) { if (X === S) failed++; continue; }
        if (X === S) missing++;
        if (!X.loading.has(name) && inFlight < PARALLEL) { fetchTile(X, name); inFlight++; }
      }
      O.sig = sigOf(P, names);
      if (P.early || P.alpha <= 0) { O.alpha = 0; O.dot = 1; O.drawn = new Set(); O.waiting = missing; say("Zoom in to see the cell outlines."); return null; }
      const a = P.alpha;
      O.alpha = a;
      // while tiles are still on their way (or failed) the dots stay half visible, so no cell goes missing
      O.dot = missing || failed ? Math.max(1 - a, 0.5) : 1 - a;
      const col = colours(S.index.nRows);
      const now = performance.now(), layers = [], drawn = new Set();
      let shown = 0, fresh = 0, held = 0;
      for (const name of names) {
        const T = S.tiles.get(name);
        if (!T) continue;
        if (!O.drawn.has(name)) {                    // new this frame: only so many at once
          if (fresh && fresh + T.n > NEW_PER_FRAME) { held++; continue; }
          fresh += T.n;
        }
        drawn.add(name);
        T.used = now; shown += T.n;
        const base = `${O.ds}-${O.key}-${S.key}-${name}`;
        const rr = T.ringRow, f = col.fill, l = col.line;
        if (T.n && (col.anyFill || col.pick)) {
          const id = `ol-fill-${base}`;
          O.byLayer.set(id, T.ringRow);
          layers.push(new D.SolidPolygonLayer({
            id, data: T.fill, positionFormat: "XY", _normalize: false, filled: true, extruded: false,
            getFillColor: (_, { index, target }) => { const k = rr[index] * 4; target[0] = f[k]; target[1] = f[k + 1]; target[2] = f[k + 2]; target[3] = f[k + 3]; return target; },
            opacity: a * col.opacity, pickable: true, autoHighlight: true, highlightColor: [255, 255, 255, 70],
            parameters: ctx.NORMAL, updateTriggers: { getFillColor: [col.key] },
          }));
        }
        if (T.n) layers.push(new D.PathLayer({
          id: `ol-line-${base}`, data: T.path, positionFormat: "XY", _pathType: "loop",
          getColor: (_, { index, target }) => { const k = rr[index] * 4; target[0] = l[k]; target[1] = l[k + 1]; target[2] = l[k + 2]; target[3] = l[k + 3]; return target; },
          widthUnits: "pixels", getWidth: 1.25, widthMinPixels: 1, jointRounded: false, capRounded: false,
          opacity: a, pickable: false, parameters: ctx.NORMAL, updateTriggers: { getColor: [col.key] },
        }));
        // cells of this tile that have no outline keep a dot, in the colour their outline would have
        if (T.dotRows.length && col.pos) {
          const dr = T.dotRows, pos = col.pos, id = `ol-dot-${base}`;
          O.byLayer.set(id, T.dotRows);
          layers.push(new D.ScatterplotLayer({
            id, data: { length: dr.length }, getPosition: (_, { index, target }) => { target[0] = pos[2 * dr[index]]; target[1] = pos[2 * dr[index] + 1]; target[2] = 0; return target; },
            getFillColor: (_, { index, target }) => { const k = dr[index] * 4; target[0] = l[k]; target[1] = l[k + 1]; target[2] = l[k + 2]; target[3] = f[k + 3] || 110; return target; },
            radiusUnits: "common", getRadius: col.radius, radiusMinPixels: 0.6, stroked: false, filled: true, opacity: a,
            pickable: col.pick, parameters: ctx.NORMAL, updateTriggers: { getFillColor: [col.key], getPosition: [name] },
          }));
        }
      }
      if (N && N.index) for (const name of want[1][1]) {
        const T = N.tiles.get(name);
        if (!T || !T.n) continue;
        T.used = now;
        layers.push(new D.PathLayer({
          id: `ol-nuc-${O.ds}-${O.key}-${name}`, data: T.path, positionFormat: "XY", _pathType: "loop",
          getColor: [NUCLEUS[0], NUCLEUS[1], NUCLEUS[2], 210], widthUnits: "pixels", getWidth: 1, widthMinPixels: 1,
          opacity: a * (col.over ? 0.45 : 1), pickable: false, parameters: ctx.NORMAL,
        }));
      }
      evict(want);
      O.drawn = drawn; O.waiting = missing + held;
      if (held) { O.dot = Math.max(O.dot, 0.5); ctx.requestRender(); }      // the rest of them next frame
      O.stats.tiles = names.length; O.stats.shown = shown;
      say(failed ? "Some outlines could not be loaded: those cells are drawn as dots."
        : missing ? "Loading cell outlines…" : "");
      return { layers, over: col.over };
    }
    const sigOf = (P, names) => (P ? `${P.alpha.toFixed(2)}|${P.early ? 1 : 0}|${names.join(",")}` : "off");

    // decoded tiles out of view are kept for a quick return, up to a number of vertices
    function evict(want) {
      for (const [X, list] of want) {
        let total = 0;
        for (const T of X.tiles.values()) total += T.nVert;
        if (total <= KEEP_VERTS) continue;
        const inView = new Set(list);
        const old = [...X.tiles.values()].filter((T) => !inView.has(T.name)).sort((a, b) => a.used - b.used);
        for (const T of old) { if (total <= KEEP_VERTS) break; X.tiles.delete(T.name); total -= T.nVert; }
      }
    }

    // has the view moved enough to need another frame (another tile, another step of the fade)?
    function viewChanged(vs) {
      if (!available() || mode !== "outlines") return false;
      const P = plan(vs), S = primary();
      let names = [];
      if (P && P.box && S && S.index) { const got = tilesInView(S, P.box); names = got[1] > BUDGET ? [] : got[0]; }
      return sigOf(P, names) !== O.sig;
    }

    // the row of positions.bin under a picked outline (or the dot of a cell without one); -1 = none
    function rowAt(info) {
      if (!info || !info.layer || info.index < 0) return -1;
      const rows = O.byLayer.get(info.layer.id);
      return rows && info.index < rows.length ? rows[info.index] : -1;
    }

    function init() {
      const sel = $("cellShapes"), nb = $("cellNuclei");
      if (sel) sel.onchange = () => {
        mode = sel.value === "dots" ? "dots" : "outlines";
        try { localStorage.setItem(KEY, mode); } catch (e) { /* private mode */ }
        syncUI(); ctx.requestRender();
      };
      if (nb) nb.onchange = () => { nuclei = nb.checked; ctx.requestRender(); };
      syncUI();
      // a link to the view and Save JSON carry the choice (viewlink.js extend()): { m: "outlines" | "dots", n: 1 = nuclei }.
      // Only for a layer with outlines. A link or file without it leaves the viewer's own choice as it is; putting a
      // link's choice back does not change what this browser remembers (that is the menu's alone).
      if (window.SodaViewLink && window.SodaViewLink.extend) window.SodaViewLink.extend({
        key: "cs",
        capture: () => (available() ? { m: mode, ...(nuclei && setOf("nucleus") ? { n: 1 } : {}) } : null),
        clean: (x, h) => {
          if (!h.isObj(x)) return null;
          const o = {};
          h.put(o, "m", h.oneOf(x.m, ["outlines", "dots"])); h.put(o, "n", h.bool(x.n));
          return o.m || o.n !== undefined ? o : null;
        },
        apply: async (v) => {
          if (!v) return [];
          if (v.m) mode = v.m;
          if (v.n !== undefined) nuclei = !!v.n && !!setOf("nucleus");
          syncUI(); ctx.requestRender();
          return [];
        },
      });
    }
    init();

    return {
      reset, onLayer, frame, viewChanged, rowAt, available,
      // how much of the dots to draw this frame: 1 = as always, 0 = the outlines have taken over
      dotAlpha: () => O.dot,
      mode: () => mode,
      // Save image's legend: what the outlines on show are (none while the cells are dots)
      legend: () => {
        if (!(O.alpha > 0) || !primary()) return [];
        const items = [{ label: (primary().meta && primary().meta.label) || "Cell outlines", outline: LINE }];
        const N = nuclei ? setOf("nucleus") : null;
        if (N && N.tiles.size) items.push({ label: (N.meta && N.meta.label) || "Nucleus outlines", outline: NUCLEUS });
        return items;
      },
      setMode: (m) => { mode = m === "dots" ? "dots" : "outlines"; syncUI(); ctx.requestRender(); },
      setNuclei: (on) => { nuclei = !!on && !!setOf("nucleus"); syncUI(); ctx.requestRender(); },
      // for tests and the console: what is loaded and what it cost
      stats: () => ({ ...O.stats, alpha: O.alpha, dot: O.dot, note: O.note, mode, nuclei, layers: O.byLayer.size,
        busy: O.waiting > 0 || O.sets.some((s) => !!s.busy || s.loading.size > 0),
        sets: O.sets.map((s) => ({ key: s.key, index: !!s.index, tiles: s.tiles.size, loading: s.loading.size, bad: s.bad.size,
          vertices: [...s.tiles.values()].reduce((n, T) => n + T.nVert, 0), polygons: [...s.tiles.values()].reduce((n, T) => n + T.n, 0) })) }),
      _debug: () => O,
    };
  }

  window.SodaOutlines = { create, FORMAT, VERSION };
})();
