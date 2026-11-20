// Save image: the map as a PNG file, as it is on screen (Tissue, UMAP or Side by side, or a single-cell
// dataset's UMAP panels), with a legend, the scale bar and the dataset's name around it.
//
// 1x is the map's own pixels: deck keeps its last frame (preserveDrawingBuffer), copied as it is.
// 2x and 4x draw the same view again with two or four times the pixels (deck's useDevicePixels) and the
// image tiles one or two levels finer (TileLayer zoomOffset, through viewer.js), wait until every tile
// on show has arrived, copy the frame and put the screen back. Dots, lines and labels are sized in
// screen pixels, so they grow with the picture: the file looks like the screen, only sharper.
// The legend, the name and the scale bar are drawn here on a 2D canvas, from what viewer.js says is
// on show (SodaSnapshot.init(api)); the picture itself is not touched.
//
// The tiles come from fetch() and are decoded to bitmaps, so a picture from another site (a ?data=
// link, a storage account) cannot taint the canvas; if a browser still refuses to hand the pixels
// back, the message says so. Nothing from the address goes into the picture or the file name: the
// name is the dataset's label as the picker shows it, never the link, its token or a password.
/* exported SodaSnapshot */
(function () {
  "use strict";
  let api = null, busy = false;
  const $ = (id) => document.getElementById(id);
  const FONT = 'system-ui, -apple-system, "Segoe UI", Roboto, sans-serif';
  // the largest picture drawn: WebGL and canvas limits, and memory (a phone reloads a tab that takes too much)
  const LIMITS = { dim: 16384, area: 50e6, dimTouch: 8192, areaTouch: 16e6, waitMs: 60000 };
  const SCALES = [1, 2, 4];
  const MODE_NAME = { tissue: "Tissue", umap: "UMAP", sbs: "Side by side" };
  const MODE_FILE = { tissue: "tissue", umap: "umap", sbs: "sidebyside" };
  const THEME = {
    black: { bg: "#000000", fg: "#e6e9ef", muted: "#9aa2b1", edge: "rgba(255,255,255,0.28)" },
    white: { bg: "#ffffff", fg: "#1a1f2b", muted: "#555d6b", edge: "rgba(0,0,0,0.35)" },
    transparent: { bg: null, fg: "#1a1f2b", muted: "#555d6b", edge: "rgba(0,0,0,0.35)" },
  };

  // ---- the controls: a section in the Share tab, and a button beside Fit ------------------------------
  function buildDom() {
    const sec = document.createElement("section");
    sec.className = "panel"; sec.id = "panelImage"; sec.dataset.tab = "share";
    sec.innerHTML = `<h2>Save image</h2>
      <label class="row slim"><span>Size</span><select id="snapScale" title="How many pixels: 1x is the map as on screen; 2x and 4x draw it again, sharper, with finer image tiles"></select></label>
      <label class="row slim"><span>Background</span><select id="snapBg">
        <option value="black">Black</option><option value="white">White</option>
        <option value="transparent">Transparent (dark text)</option></select></label>
      <label class="row slim check"><input id="snapLegend" type="checkbox" checked /><span>Legend (colours and names)</span></label>
      <label class="row slim check"><input id="snapBar" type="checkbox" checked /><span>Scale bar (tissue)</span></label>
      <label class="row slim check"><input id="snapName" type="checkbox" checked /><span>Dataset name and view</span></label>
      <div class="btnrow"><button id="snapSave" type="button" title="Save the map as it is on screen as a PNG file">Save PNG</button></div>
      <div id="snapMsg" class="meta snap-msg" role="status"></div>`;
    const side = $("sidebar");
    if (window.SodaLayout) { side.appendChild(sec); window.SodaLayout.adopt(sec); }
    else { const help = $("panelHelp"); if (help) help.before(sec); else side.appendChild(sec); }

    // the same, one click away: a camera beside Fit in the bar over the map (the options are the Share tab's)
    const quick = document.createElement("button");
    quick.type = "button"; quick.id = "snapQuick"; quick.className = "snap-quick";
    quick.title = "Save the map as a PNG image (size and options: Share tab)";
    quick.setAttribute("aria-label", "Save the map as a PNG image");
    quick.innerHTML = '<svg viewBox="0 0 24 24" aria-hidden="true"><path d="M4 8h3l2-3h6l2 3h3v11H4z"/><circle cx="12" cy="13" r="3.5"/></svg>';
    const fit = $("btnFit");
    if (fit && fit.parentElement) fit.after(quick);

    $("snapSave").onclick = save;
    quick.onclick = save;
    $("snapScale").addEventListener("pointerdown", fillScales);
    $("snapScale").addEventListener("focus", fillScales);
    window.addEventListener("resize", () => { if (!busy) fillScales(); });
    window.addEventListener("soda:shown", () => { if (!busy) fillScales(); });
    fillScales();
  }
  const msg = (t, err) => { const el = $("snapMsg"); if (!el) return; el.textContent = t; el.style.color = err ? "#ff8f8f" : ""; if (t && err && window.SodaLayout) window.SodaLayout.flag("panelImage", "warn"); };
  function sync() {
    for (const id of ["snapSave", "snapQuick"]) if ($(id)) $(id).disabled = busy;
    if ($("snapSave")) $("snapSave").textContent = busy ? "Saving…" : "Save PNG";
  }

  // ---- sizes ------------------------------------------------------------------------------------------
  const mapCanvas = () => { const d = api && api.deck(); const c = d && (d.getCanvas ? d.getCanvas() : d.canvas); return c || ($("viewport") && $("viewport").querySelector("canvas")); };
  function glLimit() {
    try {
      const d = api.deck(), gl = d && d.device && d.device.gl;
      if (!gl) return Infinity;
      const vp = gl.getParameter(gl.MAX_VIEWPORT_DIMS);
      return Math.min(gl.getParameter(gl.MAX_RENDERBUFFER_SIZE) || Infinity, vp ? Math.min(vp[0], vp[1]) : Infinity);
    } catch (e) { return Infinity; }
  }
  // the pixels of the map at scale s, and why not, when it is too large
  // On a phone or tablet the limits count the whole picture (the map, the name above it and the legend
  // beside it, with the options as set), so a picture over them is refused before anything is drawn.
  // the map's device pixel ratio at rest (deck's useDevicePixels), not read from the canvas: just after a 2x or 4x
  // save, deck puts the canvas back to its own size on its next frame, and until then it still has the larger size
  function restDpr() {
    const d = api && api.deck(), u = d ? d.props.useDevicePixels : true;
    return typeof u === "number" ? u : u === false ? 1 : (window.devicePixelRatio || 1);
  }
  function sizeAt(s) {
    const c = mapCanvas();
    if (!c || !c.clientWidth) return null;
    const r = restDpr(), w = Math.round(Math.floor(c.clientWidth * r) * s), h = Math.round(Math.floor(c.clientHeight * r) * s);
    let W = w, H = h;
    if (api.touch) {
      const opts = { legend: !$("snapLegend") || $("snapLegend").checked, name: !$("snapName") || $("snapName").checked };
      const F = frameAround(w / Math.max(1, c.clientWidth), h, api.mode(), opts);
      W = w + F.legW; H = h + F.headH;
    }
    const dim = Math.min(api.touch ? LIMITS.dimTouch : LIMITS.dim, glLimit()), area = api.touch ? LIMITS.areaTouch : LIMITS.area;
    const why = Math.max(w, h) > dim || Math.max(W, H) > (api.touch ? LIMITS.dimTouch : Infinity) ? `wider than ${dim.toLocaleString()} px, this browser's limit` : W * H > area
      ? (api.touch ? "too large for a phone or tablet's memory" : `over ${Math.round(area / 1e6)} million pixels`) : "";
    return { w, h, W, H, why };
  }
  function fillScales() {
    const sel = $("snapScale"); if (!sel || !api) return;
    const keep = Number(sel.value) || (api.touch ? 1 : 2), opts = [];
    for (const s of SCALES) {
      const z = sizeAt(s);
      const txt = z ? `${s}x · ${z.w.toLocaleString()} × ${z.h.toLocaleString()} px${s === 1 ? " (as on screen)" : ""}${z.why ? ` (too large: ${z.why})` : ""}` : `${s}x`;
      opts.push({ s, txt, off: !!(z && z.why) });
    }
    const key = opts.map((o) => `${o.txt}|${o.off}`).join(";");
    if (sel.dataset.key === key) return;
    sel.dataset.key = key; sel.innerHTML = "";
    for (const o of opts) { const e = document.createElement("option"); e.value = o.s; e.textContent = o.txt; e.disabled = o.off; sel.appendChild(e); }
    const ok = opts.filter((o) => !o.off).map((o) => o.s);
    sel.value = String(ok.includes(keep) ? keep : ok.filter((s) => s <= keep).pop() || 1);
  }

  // ---- drawing it larger: wait until every image tile on show has arrived ----------------------------
  const frame = () => new Promise((res) => { let done = false; requestAnimationFrame(() => { done = true; res(); }); setTimeout(() => { if (!done) res(); }, 120); });
  async function frames(n) { for (let i = 0; i < n; i++) await frame(); }
  const tileLayers = (deck) => (deck.layerManager ? deck.layerManager.getLayers() : []).filter((l) => l.state && l.state.tileset);
  async function waitTiles(deck, onProg) {
    const t0 = performance.now();
    let calm = 0, n = 0, done = 0;
    while (performance.now() - t0 < LIMITS.waitMs) {
      await frame();
      const tl = tileLayers(deck);
      n = 0; done = 0;
      for (const l of tl) { const sel = l.state.tileset.selectedTiles || []; n += sel.length; for (const t of sel) if (t.isLoaded) done++; }
      const all = done === n && tl.every((l) => l.isLoaded);
      calm = all ? calm + 1 : 0;
      if (calm >= 4) return { ok: true, n, done };          // a few frames in a row: no new tiles picked meanwhile
      if (onProg && n) onProg(done, n);
    }
    return { ok: false, n, done };
  }

  // ---- the file -----------------------------------------------------------------------------------------
  const pad2 = (v) => String(v).padStart(2, "0");
  function stamp() { const d = new Date(); return `${d.getFullYear()}${pad2(d.getMonth() + 1)}${pad2(d.getDate())}-${pad2(d.getHours())}${pad2(d.getMinutes())}`; }
  const slug = (s) => String(s || "").normalize("NFKD").replace(/[^A-Za-z0-9.-]+/g, "_").replace(/^[_.]+|[_.]+$/g, "").slice(0, 60) || "dataset";

  async function save() {
    if (busy || !api) return;
    const deck = api.deck(), canvas = mapCanvas();
    if (!deck || !canvas || !api.ready()) { msg("Nothing to save yet: open a dataset first (or wait for the flight between the tissue and the UMAP to end).", true); return; }
    fillScales();
    const s = Number($("snapScale").value) || 1, z = sizeAt(s);
    if (!z) { msg("The map has no size right now (is it hidden?): nothing to save.", true); return; }
    if (z.why) { msg(`${s}x is ${z.why}: pick a smaller size.`, true); return; }
    const opts = { bg: $("snapBg").value, legend: $("snapLegend").checked, bar: $("snapBar").checked, name: $("snapName").checked };
    const mode = api.mode(), notes = [];
    busy = true; sync();
    let restore = null, shot = null;
    try {
      const dpr0 = restDpr();
      if (s > 1) {
        const was = deck.props.useDevicePixels;
        restore = () => { restore = null; try { deck.setProps({ useDevicePixels: was }); api.setTileBoost(0); } catch (e) { /* the page goes on as it was */ } };
        msg(`Drawing the picture at ${s}x…`);
        deck.setProps({ useDevicePixels: dpr0 * s });
        api.setTileBoost(Math.round(Math.log2(s)));         // s times finer image tiles than the screen asks for
        await frames(2);
        const w = await waitTiles(deck, (d, n) => msg(`Drawing the picture at ${s}x: image tiles ${d} of ${n}…`));
        if (!w.ok) notes.push(`${w.n - w.done} of ${w.n} image tiles had not arrived after ${LIMITS.waitMs / 1000} s, so parts of the image may be coarser than the rest.`);
        if (Math.abs(canvas.width - z.w) > 2) notes.push(`The browser drew ${canvas.width.toLocaleString()} px wide instead of ${z.w.toLocaleString()}.`);
      }
      deck.redraw(true);                                     // draws now, into the buffer copied next
      shot = copyFrame(deck, canvas);
      if (restore) restore();
      if (shot.small) notes.push(`This browser limits the picture to ${shot.w.toLocaleString()} × ${shot.h.toLocaleString()} px.`);
      const out = compose(shot, mode, opts);
      const name = `${slug(api.label())}_${MODE_FILE[mode] || mode}_${stamp()}.png`;
      const blob = await toPng(out);
      api.download(name, blob, "image/png");
      msg(`Saved ${name} (${out.width.toLocaleString()} × ${out.height.toLocaleString()} px).` + (notes.length ? " " + notes.join(" ") : "") +
          (api.touch ? " On a phone or tablet it goes to Downloads (iPhone and iPad: the Files app)." : ""), notes.length > 0);
    } catch (e) {
      msg(`The image could not be saved: ${e && e.message ? e.message : e}`, true);
    } finally {
      if (restore) restore();
      busy = false; sync();
      if (api.render) api.render();
    }
  }
  // the map's frame on a canvas of its own, at the pixels the GPU actually drew
  function copyFrame(deck, canvas) {
    const gl = deck.device && deck.device.gl;
    const w = gl && gl.drawingBufferWidth ? gl.drawingBufferWidth : canvas.width, h = gl && gl.drawingBufferHeight ? gl.drawingBufferHeight : canvas.height;
    const c = document.createElement("canvas");
    c.width = w; c.height = h;
    const ctx = c.getContext("2d");
    if (!ctx) throw new Error(`the browser could not make a ${w} × ${h} px picture (memory): pick a smaller size.`);
    ctx.drawImage(canvas, 0, 0, w, h);
    return { canvas: c, w, h, k: w / Math.max(1, canvas.clientWidth), small: w < canvas.width - 2, rect: canvas.getBoundingClientRect() };
  }
  function toPng(c) {
    return new Promise((res, rej) => {
      try {
        c.toBlob((b) => (b ? res(b) : rej(new Error("the browser could not encode the PNG (most often: not enough memory for this size). Pick a smaller size."))), "image/png");
      } catch (e) {
        rej(new Error(e && e.name === "SecurityError"
          ? "the browser will not let the page read the picture back, because part of it came from another site that does not allow that. Use your system's screenshot tool instead."
          : (e && e.message) || String(e)));
      }
    });
  }

  // ---- the picture with its legend, name and scale bar -------------------------------------------------
  // what goes around a map of ph pixels high, at k picture pixels per screen pixel: the name and the view
  // above it (headH high), the legend beside it (legW wide)
  function frameAround(k, ph, mode, opts) {
    const font = (px, bold) => `${bold ? "600 " : ""}${Math.round(px * k)}px ${FONT}`;
    const meas = document.createElement("canvas").getContext("2d");
    const pad = Math.round(10 * k);
    const head = opts.name ? [api.label(), api.layer(), MODE_NAME[mode] || mode].filter(Boolean).join("  ·  ") : "";
    const headH = head ? Math.round(30 * k) : 0;
    const L = opts.legend && meas ? legendLayout(meas, font, k, ph - 2 * pad) : null;
    const legW = L && L.w ? Math.ceil(L.w + 2 * pad) : 0;               // a margin on both sides
    return { font, pad, head, headH, L, legW };
  }
  function compose(shot, mode, opts) {
    const T = THEME[opts.bg] || THEME.black, k = shot.k, pw = shot.w, ph = shot.h;
    const { font, pad, head, headH, L, legW } = frameAround(k, ph, mode, opts);
    const out = document.createElement("canvas");
    out.width = pw + legW; out.height = ph + headH;
    const ctx = out.getContext("2d");
    if (!ctx) throw new Error("the browser could not make a picture this large (memory): pick a smaller size.");
    if (T.bg) { ctx.fillStyle = T.bg; ctx.fillRect(0, 0, out.width, out.height); }
    ctx.drawImage(shot.canvas, 0, headH);
    ctx.textBaseline = "middle";
    if (head) {
      ctx.font = font(14, true); ctx.fillStyle = T.fg;
      ctx.fillText(fitText(ctx, head, out.width - 2 * pad), pad, headH / 2);
    }
    drawCaptions(ctx, shot, headH, k, font);
    if (opts.bar) drawScaleBar(ctx, shot, headH, k, font);
    if (L && L.w) drawLegend(ctx, L, pw + pad, headH + pad, k, font, T);
    return out;
  }
  function fitText(ctx, t, max) {
    if (ctx.measureText(t).width <= max) return t;
    let lo = 0, hi = t.length;
    while (lo < hi) { const m = (lo + hi + 1) >> 1; if (ctx.measureText(t.slice(0, m) + "…").width <= max) lo = m; else hi = m - 1; }
    return t.slice(0, lo) + "…";
  }
  const rgbCss = (c, a) => `rgba(${c[0]},${c[1]},${c[2]},${a == null ? 1 : a})`;
  function plate(ctx, x, y, w, h, r) {
    ctx.fillStyle = "rgba(0,0,0,0.6)";
    ctx.beginPath();
    if (ctx.roundRect) ctx.roundRect(x, y, w, h, r); else ctx.rect(x, y, w, h);
    ctx.fill();
  }

  // the captions over the map (Side by side: "Tissue", "UMAP · whole section"; a single-cell dataset: each
  // panel's title), where they are on screen
  function drawCaptions(ctx, shot, top, k, font) {
    const els = [];
    const caps = $("sbsCaps"), labs = $("umapLabels");
    if (caps && !caps.hidden) els.push(...caps.querySelectorAll("span"));
    if (labs && !labs.hidden) {
      // the panel grid lines first
      ctx.strokeStyle = "#1c2230"; ctx.lineWidth = Math.max(1, k);
      for (const c of labs.querySelectorAll(".umap-cell")) {
        const r = c.getBoundingClientRect(), x = (r.left - shot.rect.left) * k, y = top + (r.top - shot.rect.top) * k;
        ctx.beginPath(); ctx.moveTo(x + r.width * k, y); ctx.lineTo(x + r.width * k, y + r.height * k); ctx.lineTo(x, y + r.height * k); ctx.stroke();
      }
      els.push(...labs.querySelectorAll(".umap-cell > .t"));
    }
    if (caps && !caps.hidden) {                   // the line between the two halves
      const vps = api.deck().getViewports(), side = vps.find((v) => v.id === "uside");
      if (side) {
        ctx.strokeStyle = "#3a4356"; ctx.lineWidth = Math.max(1, k); ctx.beginPath();
        if (side.x > 0) { ctx.moveTo(side.x * k, top); ctx.lineTo(side.x * k, top + shot.h); } else { ctx.moveTo(0, top + side.y * k); ctx.lineTo(shot.w, top + side.y * k); }
        ctx.stroke();
      }
    }
    for (const el of els) {
      const t = el.textContent.replace(/\s+/g, " ").trim();
      if (!t || !el.getClientRects().length) continue;
      const r = el.getBoundingClientRect(), x = (r.left - shot.rect.left) * k, y = top + (r.top - shot.rect.top) * k;
      ctx.font = font(12.5, true);
      const w = ctx.measureText(t).width + 12 * k, h = 20 * k;
      plate(ctx, x, y, w, h, 3 * k);
      ctx.fillStyle = "#e6e9ef"; ctx.fillText(t, x + 6 * k, y + h / 2);
    }
  }

  // the scale bar, at the bottom right of the slide (the left half when side by side), as on screen
  function drawScaleBar(ctx, shot, top, k, font) {
    const sc = api.scale(); if (!sc || !Number.isFinite(sc.zoom)) return;
    const vp = api.deck().getViewports().find((v) => v.id === "ortho"); if (!vp) return;
    const perCss = (sc.um || 1) / Math.pow(2, sc.zoom);              // µm (or image px) per screen px
    const nice = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000];
    const len = nice.find((n) => n / perCss >= 80) || nice[nice.length - 1];
    const bar = (len / perCss) * k, lbl = sc.um ? (len >= 1000 ? `${len / 1000} mm` : `${len} µm`) : `${len} px`;
    ctx.font = font(12, true);
    const tw = ctx.measureText(lbl).width, w = Math.max(bar, tw) + 16 * k, h = 30 * k;
    const x = (vp.x + vp.width) * k - 14 * k - w, y = top + (vp.y + vp.height) * k - 12 * k - h;
    plate(ctx, x, y, w, h, 4 * k);
    ctx.fillStyle = "#ffffff";
    ctx.fillRect(x + (w - bar) / 2, y + h - 10 * k, bar, 4 * k);
    ctx.textAlign = "center"; ctx.fillText(lbl, x + w / 2, y + 10 * k); ctx.textAlign = "left";
  }

  // the legend in columns as tall as the picture: section titles, then a swatch and a name per row
  function legendLayout(meas, font, k, maxH) {
    let secs = [];
    try { secs = api.legend() || []; } catch (e) { secs = []; }
    const MAXROWS = 300, rowH = 19 * k, titleH = 21 * k, gap = 8 * k, sw = 14 * k, grad = 34 * k, maxText = 240 * k;
    const ents = [];
    let rows = 0;
    for (const sct of secs) {
      if (!sct.items || !sct.items.length) continue;
      if (ents.length) ents.push({ kind: "gap", h: gap });
      ents.push({ kind: "title", text: sct.title, h: titleH });
      const items = sct.items.slice(0, Math.max(0, MAXROWS - rows));
      for (const it of items) ents.push({ kind: "item", it, h: rowH });
      rows += items.length;
      if (items.length < sct.items.length) { ents.push({ kind: "more", text: `+ ${sct.items.length - items.length} more`, h: rowH }); break; }
    }
    if (!ents.length) return null;
    const limit = Math.max(maxH, 6 * rowH + titleH), cols = [];
    let col = { ents: [], h: 0, w: 0 };
    for (const e of ents) {
      if (col.h + e.h > limit && col.ents.length) { cols.push(col); col = { ents: [], h: 0, w: 0 }; if (e.kind === "gap") continue; }
      col.ents.push(e); col.h += e.h;
    }
    cols.push(col);
    for (const c of cols) for (const e of c.ents) {
      if (e.kind === "title") { meas.font = font(12, true); e.w = Math.min(maxText + grad, meas.measureText(e.text).width); }
      else if (e.kind === "item") { meas.font = font(12); e.sw = e.it.stops ? grad : sw; e.w = e.sw + 6 * k + Math.min(maxText, meas.measureText(e.it.label).width); }
      else if (e.kind === "more") { meas.font = font(12); e.w = meas.measureText(e.text).width; }
      else e.w = 0;
      c.w = Math.max(c.w, e.w);
    }
    const colGap = 14 * k;
    return { cols, w: cols.reduce((a, c) => a + c.w, 0) + colGap * (cols.length - 1), colGap, maxText, sw };
  }
  function drawLegend(ctx, L, x0, y0, k, font, T) {
    let x = x0;
    for (const c of L.cols) {
      let y = y0;
      for (const e of c.ents) {
        const mid = y + e.h / 2;
        if (e.kind === "title") {
          ctx.font = font(12, true); ctx.fillStyle = T.fg; ctx.fillText(fitText(ctx, e.text, c.w), x, mid + 1 * k);
        } else if (e.kind === "item") {
          const it = e.it;
          if (it.stops) {
            const g = ctx.createLinearGradient(x, 0, x + e.sw, 0);
            it.stops.forEach((s, i) => g.addColorStop(it.stops.length > 1 ? i / (it.stops.length - 1) : 0, rgbCss(s)));
            ctx.fillStyle = g; ctx.fillRect(x, mid - 5 * k, e.sw, 10 * k);
            ctx.strokeStyle = T.edge; ctx.lineWidth = Math.max(1, 0.75 * k); ctx.strokeRect(x, mid - 5 * k, e.sw, 10 * k);
          } else if (it.outline) {
            ctx.strokeStyle = rgbCss(it.outline); ctx.lineWidth = 2 * k; ctx.strokeRect(x + k, mid - L.sw / 2 + k, L.sw - 2 * k, L.sw - 2 * k);
          } else if (it.rgb) {
            ctx.fillStyle = rgbCss(it.rgb); ctx.fillRect(x, mid - L.sw / 2, L.sw, L.sw);
            ctx.strokeStyle = T.edge; ctx.lineWidth = Math.max(1, 0.75 * k); ctx.strokeRect(x, mid - L.sw / 2, L.sw, L.sw);
          }
          ctx.font = font(12); ctx.fillStyle = T.fg;
          ctx.fillText(fitText(ctx, it.label, L.maxText), x + e.sw + 6 * k, mid + 1 * k);
        } else if (e.kind === "more") {
          ctx.font = font(12); ctx.fillStyle = T.muted; ctx.fillText(e.text, x, mid + 1 * k);
        }
        y += e.h;
      }
      x += c.w + L.colGap;
    }
  }

  window.SodaSnapshot = {
    init(a) { api = a; buildDom(); },
    save,
    // for tests from the console: each size with the map's pixels (w, h), the whole picture's on a phone or
    // tablet (W, H; the map's elsewhere) and why it is not offered; the limits
    _sizes: () => ({ limits: { ...LIMITS }, touch: !!(api && api.touch), sizes: SCALES.map((s) => ({ s, ...(sizeAt(s) || {}) })) }),
    // for tests from the console: set the options and save, resolving when it is done
    _save: async (o = {}) => {
      if (o.scale) { fillScales(); $("snapScale").value = String(o.scale); }
      if (o.bg) $("snapBg").value = o.bg;
      for (const [key, id] of [["legend", "snapLegend"], ["bar", "snapBar"], ["name", "snapName"]]) if (key in o) $(id).checked = !!o[key];
      await save();
      return $("snapMsg").textContent;
    },
  };
})();
