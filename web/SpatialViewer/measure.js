// Measure: a ruler (the distance between two points) and an area (a shape drawn round something),
// in µm and µm² from the export's microns per pixel (manifest pixelSizeUm), drawn on the map.
//
//   Ruler   click two points, or press and drag; a finger drags or taps twice.
//   Area    click the corners and close the shape with a double-click, Enter, Finish, or a click on
//           its first corner; or press and drag round the shape (it closes itself on release).
//           The area is the shoelace formula on the corners (a shape crossing itself counts its
//           loops with their signs, so draw round, not in figures of eight); the perimeter is
//           every side, the closing one included.
//   Keys    the map pad takes the focus when a tool is turned on from the keyboard: the arrow keys
//           move the picture, + and - zoom, Enter or Space puts a point at the centre (the cross),
//           Enter again without moving closes an area, Backspace takes the last point back, Escape
//           drops the shape being drawn, a second Escape turns the tool off.
//
// World units on the slide are full-resolution image pixels (as the scale bar and the HUD read
// them), so a length is pixels x pixelSizeUm and an area pixels² x pixelSizeUm². An export without a
// pixel size gives lengths in pixels and says so.
//
// The shapes are deck.gl layers on top of the tissue view (viewer.js render()), so Save image has
// them, and its legend lists each measurement. They travel in Save JSON and "Copy link to this view"
// as the view block's `meas` entry (viewlink.js extend()): absent in older files and links. They
// belong to the dataset: opening another one clears them.
/* exported SodaMeasure */
(function () {
  "use strict";
  let api = null;
  const $ = (id) => document.getElementById(id);
  const esc = (x) => String(x).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
  const MAX_ITEMS = 50, MAX_PTS = 2000;
  const RGB = [255, 214, 10];                      // the lines and corners
  const FILL = [255, 214, 10, 46];
  const DRAG_PX = 6;                               // a press that moves this far is a drag, not a click
  const CLOSE_PX = 12;                             // a click this close to the first corner closes an area
  const STROKE_PX = 3;
  const SAME_PX = 3;                               // a click this close to the last point is the same point                             // a drawn stroke keeps a point every few screen pixels
  const PAN_PX = 60;                               // an arrow key moves the picture this far

  // tool: null | "line" | "area"; items: finished shapes; cur: the shape being drawn
  const M = { tool: null, items: [], cur: null, hover: null, drag: null, ver: 0, next: 1 };

  // ---- lengths and areas ------------------------------------------------------------------------------
  const um = () => (api ? api.um() : 0);
  function lengthPx(p, closed) {
    let s = 0;
    for (let i = 1; i < p.length; i++) s += Math.hypot(p[i][0] - p[i - 1][0], p[i][1] - p[i - 1][1]);
    if (closed && p.length > 2) s += Math.hypot(p[0][0] - p[p.length - 1][0], p[0][1] - p[p.length - 1][1]);
    return s;
  }
  function areaPx(p) {
    let s = 0;
    for (let i = 0, j = p.length - 1; i < p.length; j = i++) s += p[j][0] * p[i][1] - p[i][0] * p[j][1];
    return Math.abs(s) / 2;
  }
  const num = (v, d) => v.toLocaleString(undefined, { minimumFractionDigits: d, maximumFractionDigits: d });
  function fmtLen(px) {
    const u = um();
    if (!u) return `${num(px, px < 100 ? 1 : 0)} px`;
    const v = px * u;
    return v >= 1000 ? `${num(v, 0)} µm (${num(v / 1000, 2)} mm)` : `${num(v, v < 10 ? 2 : 1)} µm`;
  }
  function fmtArea(px2) {
    const u = um();
    if (!u) return `${num(px2, 0)} px²`;
    const v = px2 * u * u;
    return v >= 1e6 ? `${num(v, 0)} µm² (${num(v / 1e6, 3)} mm²)` : `${num(v, v < 100 ? 1 : 0)} µm²`;
  }
  // the numbers of one shape, for the list, the labels, the legend and the tests
  function stats(it) {
    if (it.k === "line") { const l = lengthPx(it.p, false); return { lengthPx: l, um: um() ? l * um() : null, text: fmtLen(l) }; }
    const a = areaPx(it.p), per = lengthPx(it.p, true);
    return { areaPx: a, perimeterPx: per, um2: um() ? a * um() * um() : null, umPer: um() ? per * um() : null, text: fmtArea(a), per: fmtLen(per) };
  }
  const nameOf = (it) => `${it.k === "line" ? "Line" : "Area"} ${it.n}`;

  // ---- the section ---------------------------------------------------------------------------------
  function buildDom() {
    const sec = document.createElement("section");
    sec.className = "panel"; sec.id = "panelMeasure"; sec.dataset.tab = "view"; sec.hidden = true;
    sec.innerHTML = `<h2>Measure</h2>
      <div class="btnrow meas-tools">
        <button id="measLine" type="button" aria-pressed="false" title="Measure a distance: click two points, or drag">Ruler</button>
        <button id="measArea" type="button" aria-pressed="false" title="Measure an area: click its corners, or drag round it">Area</button>
        <button id="measFinish" type="button" hidden title="Close the shape being drawn (Enter)">Finish</button>
        <button id="measClear" type="button" disabled title="Remove every measurement">Clear all</button>
      </div>
      <div id="measHint" class="meta meas-hint"></div>
      <ol id="measList" class="meas-list"></ol>`;
    const side = $("sidebar");
    if (window.SodaLayout) { side.appendChild(sec); window.SodaLayout.adopt(sec); }
    else side.appendChild(sec);
    // a tool turned on from the keyboard (detail 0: Enter or Space, not a click) hands the map pad the focus
    $("measLine").onclick = (e) => { setTool(M.tool === "line" ? null : "line"); if (M.tool && e.detail === 0) $("measurePad").focus(); };
    $("measArea").onclick = (e) => { setTool(M.tool === "area" ? null : "area"); if (M.tool && e.detail === 0) $("measurePad").focus(); };
    $("measFinish").onclick = () => finish();
    $("measClear").onclick = () => { M.items = []; M.cur = null; M.next = 1; changed(); };
    $("measList").onclick = (e) => {
      const b = e.target.closest("button[data-i]"); if (!b) return;
      M.items.splice(+b.dataset.i, 1); changed();
      const left = $("measList").querySelectorAll("button[data-i]");
      if (left.length) left[Math.min(+b.dataset.i, left.length - 1)].focus(); else $("measClear").focus();
    };

    const pad = document.createElement("div");
    pad.id = "measurePad"; pad.className = "draw-pad meas-pad"; pad.hidden = true; pad.tabIndex = 0;
    pad.setAttribute("role", "application");
    pad.setAttribute("aria-label", "Measuring pad: arrow keys move the picture, plus and minus zoom, Enter or Space puts a point at the centre, Enter again without moving closes an area, Backspace takes the last point back, Escape stops");
    const after = $("selPad") || $("drawPad");
    if (after) after.after(pad); else $("viewport").after(pad);
    padEvents(pad);
  }
  const hint = (t) => { $("measHint").textContent = t; };
  function syncUi() {
    for (const [id, t] of [["measLine", "line"], ["measArea", "area"]]) {
      const b = $(id); b.classList.toggle("on", M.tool === t); b.setAttribute("aria-pressed", String(M.tool === t));
    }
    $("measurePad").hidden = !M.tool;
    $("measFinish").hidden = !(M.tool === "area" && M.cur && M.cur.p.length >= 3);
    $("measClear").disabled = !M.items.length;
    const noUm = um() ? "" : " This dataset does not say how many µm a pixel is: lengths are in image pixels.";
    hint(M.tool === "line" ? `Click two points, or press and drag. Esc stops.${noUm}`
      : M.tool === "area" ? `Click the corners, then double-click, press Enter or click the first corner to close it; or press and drag round it. Esc stops.${noUm}`
      : M.items.length ? noUm.trim() : `Ruler measures a distance, Area the inside of a shape you draw.${noUm}`);
    $("measList").innerHTML = M.items.map((it, i) => {
      const s = stats(it);
      return `<li><span class="nm">${esc(nameOf(it))}</span><span class="v">${esc(s.text)}${it.k === "area" ? `<small> · perimeter ${esc(s.per)}</small>` : ""}</span>` +
        `<button type="button" class="x" data-i="${i}" aria-label="Remove ${esc(nameOf(it))}" title="Remove">✕</button></li>`;
    }).join("");
  }
  function changed() { M.ver++; syncUi(); api.render(); }

  function setTool(t) {
    t = t === "line" || t === "area" ? t : null;
    if (t && !api.ok()) t = null;                  // the tissue view only
    if (t) api.othersOff();                        // one drawing tool at a time
    if (t !== M.tool) { M.cur = null; M.hover = null; }
    M.tool = t;
    changed();
  }

  // ---- drawing -------------------------------------------------------------------------------------
  const vp = () => api.viewport();
  function local(ev) { const r = $("measurePad").getBoundingClientRect(); return [ev.clientX - r.left, ev.clientY - r.top]; }
  const world = (xy) => vp().unproject(xy).slice(0, 2);
  const screen = (w) => vp().project(w);
  function add(kind, p) {
    if (M.items.length >= MAX_ITEMS) { hint(`At most ${MAX_ITEMS} measurements: remove one first.`); return; }
    M.items.push({ k: kind, p: p.slice(0, MAX_PTS), n: M.next++ });
  }
  // a click (or tap, or Enter at the centre) at world point w
  function point(w) {
    if (!M.tool) return;
    if (!M.cur) M.cur = { k: M.tool, p: [] };
    const p = M.cur.p;
    if (M.tool === "area" && p.length >= 3) {
      const a = screen(p[0]), b = screen(w);
      if (Math.hypot(a[0] - b[0], a[1] - b[1]) <= CLOSE_PX) { finish(); return; }
    }
    if (p.length) { const a = screen(p[p.length - 1]), b = screen(w); if (Math.hypot(a[0] - b[0], a[1] - b[1]) < SAME_PX) return; }   // the same spot twice (a double-click)
    p.push(w);
    if (M.tool === "line" && p.length === 2) { add("line", p); M.cur = null; }
    changed();
  }
  function finish() {
    const c = M.cur;
    if (c && c.k === "area" && c.p.length >= 3) { add("area", c.p); M.cur = null; changed(); return true; }
    if (c && c.k === "area") hint("An area needs at least three corners.");
    return false;
  }
  function undo() { if (M.cur && M.cur.p.length) { M.cur.p.pop(); if (!M.cur.p.length) M.cur = null; changed(); } }
  function escape() { if (M.cur) { M.cur = null; changed(); } else setTool(null); }

  function padEvents(pad) {
    let D = null;                                  // a press: where it started, whether it moved, the stroke
    if (api.pinch) api.pinch(pad, () => { D = null; M.drag = null; api.render(); });   // two fingers move the picture
    pad.addEventListener("pointerdown", (ev) => {
      if (!M.tool || (ev.pointerType === "mouse" && ev.button !== 0)) return;
      pad.setPointerCapture(ev.pointerId);
      const s = local(ev);
      D = { s, w: world(s), moved: false, stroke: [world(s)], last: s };
    });
    pad.addEventListener("pointermove", (ev) => {
      const s = local(ev);
      if (!D) {                                    // the mouse over the pad: the rubber band to the next point
        if (M.cur && ev.pointerType === "mouse") { M.hover = world(s); M.ver++; api.render(); }
        return;
      }
      if (!D.moved && Math.hypot(s[0] - D.s[0], s[1] - D.s[1]) >= DRAG_PX) D.moved = true;
      if (!D.moved) return;
      if (M.cur && M.cur.p.length) return;         // clicking corners: a drag in between does nothing
      if (M.tool === "line") D.stroke = [D.w, world(s)];
      else if (Math.hypot(s[0] - D.last[0], s[1] - D.last[1]) >= STROKE_PX) { D.stroke.push(world(s)); D.last = s; }
      M.drag = { k: M.tool, p: D.stroke };
      M.ver++; api.render();
    });
    const up = (ev) => {
      const d = D; D = null; M.drag = null;
      if (!d || !M.tool) return;
      if (ev.type === "pointercancel") { api.render(); return; }
      if (!d.moved) { point(world(local(ev))); return; }
      if (M.cur && M.cur.p.length) return;
      if (M.tool === "line" && d.stroke.length === 2) add("line", d.stroke);
      else if (M.tool === "area" && d.stroke.length >= 3) add("area", d.stroke);
      changed();
    };
    pad.addEventListener("pointerup", up);
    pad.addEventListener("pointercancel", up);
    pad.addEventListener("pointerleave", () => { if (M.hover) { M.hover = null; M.ver++; api.render(); } });
    // a double-click closes an area (its second click fell on the first's point, so it put none down)
    pad.addEventListener("dblclick", (ev) => { ev.preventDefault(); if (M.tool === "area") finish(); });
    pad.addEventListener("wheel", (ev) => { ev.preventDefault(); zoomBy(ev.deltaY < 0 ? 0.25 : -0.25); }, { passive: false });
    pad.addEventListener("keydown", (ev) => {
      const k = ev.key;
      if (k === "ArrowLeft" || k === "ArrowRight" || k === "ArrowUp" || k === "ArrowDown") {
        const s = Math.pow(2, api.viewState().zoom), d = (ev.shiftKey ? 4 : 1) * PAN_PX / s, vs = api.viewState(), t = vs.target;
        const dx = k === "ArrowLeft" ? -d : k === "ArrowRight" ? d : 0, dy = k === "ArrowUp" ? -d : k === "ArrowDown" ? d : 0;
        api.camera({ ...vs, target: [t[0] + dx, t[1] + dy, 0] });
      } else if (k === "+" || k === "=") zoomBy(0.5);
      else if (k === "-" || k === "_") zoomBy(-0.5);
      else if (k === "Enter" || k === " ") {
        const c = M.cur, t = api.viewState().target;
        if (k === "Enter" && M.tool === "area" && c && c.p.length >= 3 && dist(screen(c.p[c.p.length - 1]), screen([t[0], t[1]])) < SAME_PX) finish();
        else centrePoint();
      }
      else if (k === "Backspace" || k === "Delete") undo();
      else if (k === "Escape") escape();
      else return;
      ev.preventDefault(); ev.stopPropagation();
    });
  }
  function zoomBy(dz) {
    const vs = api.viewState();
    const lo = Number.isFinite(vs.minZoom) ? vs.minZoom : -Infinity, hi = Number.isFinite(vs.maxZoom) ? vs.maxZoom : Infinity;
    api.camera({ ...vs, zoom: Math.max(lo, Math.min(hi, vs.zoom + dz)) });
  }
  // Enter / Space on the pad: a point at the centre cross; Enter on the first corner again closes the area
  function centrePoint() { const t = api.viewState().target; point([t[0], t[1]]); }
  // Escape and Enter while a tool is on but the focus is elsewhere on the page (not in a box)
  function keys(e) {
    if (!M.tool || e.defaultPrevented) return;
    if (e.target && (/^(INPUT|SELECT|TEXTAREA)$/.test(e.target.tagName) || e.target.id === "measurePad")) return;
    if (e.key === "Escape") { escape(); e.preventDefault(); }
    else if (e.key === "Enter" && M.tool === "area" && M.cur && M.cur.p.length >= 3 && !(e.target && e.target.closest && e.target.closest("button"))) { finish(); e.preventDefault(); }
    else if (e.key === "Backspace" && M.cur && !(e.target && e.target.closest && e.target.closest("button"))) { undo(); e.preventDefault(); }
  }

  // ---- what viewer.js draws ---------------------------------------------------------------------------
  function layers() {
    if (!api || !api.ok() || (!M.items.length && !M.cur && !M.tool)) return [];
    const D = window.deck, out = [], lines = [], polys = [], dots = [], labels = [];
    const u = M.ver;
    for (const it of M.items) {
      const s = stats(it);
      if (it.k === "line") { lines.push(it.p); dots.push(...it.p); labels.push({ at: mid(it.p[0], it.p[1]), text: s.text }); }
      else { polys.push(it.p); lines.push([...it.p, it.p[0]]); labels.push({ at: centroid(it.p), text: s.text }); }
    }
    const c = M.drag && M.drag.p.length >= 2 ? M.drag : M.cur;
    if (c && c.p.length) {
      const p = M.hover && c === M.cur ? [...c.p, M.hover] : c.p;
      if (c === M.cur) dots.push(...c.p);
      if (p.length >= 2) lines.push(c.k === "area" && p.length >= 3 ? [...p, p[0]] : p);
      if (c.k === "area" && p.length >= 3) { polys.push(p); labels.push({ at: centroid(p), text: fmtArea(areaPx(p)) }); }
      else if (c.k === "line" && p.length === 2) labels.push({ at: mid(p[0], p[1]), text: fmtLen(lengthPx(p, false)) });
    }
    if (polys.length) out.push(new D.PolygonLayer({ id: "measure-fill", data: polys, getPolygon: (d) => d, filled: true, stroked: false,
      getFillColor: FILL, parameters: api.NORMAL, pickable: false, updateTriggers: { getPolygon: [u] } }));
    if (lines.length) {
      out.push(new D.PathLayer({ id: "measure-under", data: lines, getPath: (d) => d, getColor: [0, 0, 0, 200], widthUnits: "pixels", getWidth: 4.5,
        parameters: api.NORMAL, pickable: false, updateTriggers: { getPath: [u] } }));
      out.push(new D.PathLayer({ id: "measure-line", data: lines, getPath: (d) => d, getColor: [...RGB, 255], widthUnits: "pixels", getWidth: 2,
        parameters: api.NORMAL, pickable: false, updateTriggers: { getPath: [u] } }));
    }
    if (dots.length) out.push(new D.ScatterplotLayer({ id: "measure-dots", data: dots, getPosition: (d) => d, radiusUnits: "pixels", getRadius: 3.5,
      getFillColor: [...RGB, 255], stroked: true, lineWidthUnits: "pixels", getLineWidth: 1, getLineColor: [0, 0, 0, 255],
      parameters: api.NORMAL, pickable: false, updateTriggers: { getPosition: [u] } }));
    if (labels.length) out.push(new D.TextLayer({ id: "measure-labels", data: labels, getPosition: (d) => d.at, getText: (d) => d.text,
      sizeUnits: "pixels", getSize: 13, characterSet: "auto", fontFamily: "system-ui, sans-serif", fontWeight: 700, getColor: [255, 236, 140, 255],
      background: true, getBackgroundColor: [8, 11, 18, 220], backgroundPadding: [5, 2, 5, 2], getPixelOffset: [0, -14],
      parameters: api.NORMAL, pickable: false, updateTriggers: { getText: [u], getPosition: [u] } }));
    return out;
  }
  const dist = (a, b) => Math.hypot(a[0] - b[0], a[1] - b[1]);
  const mid = (a, b) => [(a[0] + b[0]) / 2, (a[1] + b[1]) / 2];
  // the label's place: the polygon's own centroid (area-weighted), the mean of its corners if degenerate
  function centroid(p) {
    let a = 0, x = 0, y = 0;
    for (let i = 0, j = p.length - 1; i < p.length; j = i++) { const f = p[j][0] * p[i][1] - p[i][0] * p[j][1]; a += f; x += (p[j][0] + p[i][0]) * f; y += (p[j][1] + p[i][1]) * f; }
    if (Math.abs(a) < 1e-9) return [p.reduce((s, q) => s + q[0], 0) / p.length, p.reduce((s, q) => s + q[1], 0) / p.length];
    return [x / (3 * a), y / (3 * a)];
  }
  function legend() {
    if (!api || !api.ok() || !M.items.length) return [];
    return [{ title: um() ? "Measurements" : "Measurements (image pixels)", items: M.items.map((it) => { const s = stats(it); return { label: `${nameOf(it)}: ${s.text}${it.k === "area" ? `, perimeter ${s.per}` : ""}`, outline: RGB }; }) }];
  }
  const summary = () => (M.items.length ? `${M.items.length} drawn` : "");

  // ---- Save JSON and the link (viewlink.js) ------------------------------------------------------------
  const r1 = (v) => Math.round(v * 10) / 10;
  function capture() {
    if (!M.items.length) return null;
    return { s: M.items.map((it) => ({ k: it.k === "area" ? "a" : "l", p: it.p.map((q) => [r1(q[0]), r1(q[1])]) })) };
  }
  function clean(x, h) {
    if (!h.isObj(x)) return null;
    const pt = (q) => { if (!Array.isArray(q) || q.length !== 2) return null; const a = h.num(q[0], -1e9, 1e9), b = h.num(q[1], -1e9, 1e9); return a === null || b === null ? null : [a, b]; };
    const s = h.list(x.s, MAX_ITEMS, (it) => {
      if (!h.isObj(it)) return null;
      const k = h.oneOf(it.k, ["l", "a"]), p = h.list(it.p, MAX_PTS, pt);
      if (!k || !p) return null;
      if (k === "l") return p.length >= 2 ? { k, p: p.slice(0, 2) } : null;
      return p.length >= 3 ? { k, p } : null;
    });
    return s && s.length ? { s } : null;
  }
  function apply(x) {
    M.items = []; M.cur = null; M.next = 1;
    if (x && Array.isArray(x.s)) for (const it of x.s) add(it.k === "a" ? "area" : "line", it.p.map((q) => [q[0], q[1]]));
    changed();
    return [];
  }

  window.SodaMeasure = {
    init(a) {
      api = a;
      buildDom();
      window.addEventListener("keydown", keys);
      if (window.SodaViewLink && window.SodaViewLink.extend) window.SodaViewLink.extend({ key: "meas", capture, clean, apply });
      syncUi();
    },
    // a dataset opened (or the view changed): measurements belong to one dataset's pixels
    reset() { M.items = []; M.cur = null; M.hover = null; M.drag = null; M.next = 1; $("panelMeasure").hidden = true; setTool(null); },
    // the section shows once a slide is drawn (viewer.js render(); never for an RNA-seq set)
    show(on) { const sec = $("panelMeasure"); if (sec.hidden === !on) return; sec.hidden = !on; if (!on && M.tool) setTool(null); else syncUi(); },
    setTool, finish, layers, legend, summary,
    // "Done drawing" on a touch screen: close an area being drawn, then put the tool away
    done() { if (M.tool === "area") finish(); setTool(null); },
    active: () => !!M.tool,
    // tests
    _state: () => M, _stats: (i) => stats(M.items[i]),
  };
})();
