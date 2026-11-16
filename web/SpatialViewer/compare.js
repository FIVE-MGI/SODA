// Compare side by side: a second SODA window on another dataset, linked to this one.
//
// "Open beside" (Save tab, and the Compare button in the top bar) opens another SODA window, sized
// to the right half of the screen, on the dataset picked. The windows opened from one another form
// a GROUP (up to MAX_WINDOWS). While "Link views" is on, the window you move leads and the others
// follow to the same place in real-world units:
//
//     shared centre (µm) = view centre (image px) x pixel size (µm/px) - this window's offset (µm)
//     shared scale       = zoom - log2(pixel size)          = log2(screen px per µm)
//
// so two datasets with different pixel sizes and image sizes show the same physical place at the
// same physical scale. The offset is each window's own: "Line up centres" sets it to the middle of
// the picture (what a new window starts with), "Link from here" to whatever the window shows now,
// so two sections that sit differently on their slides can be lined up by eye and then locked.
// SODA's view has no rotation or mirroring, so none is linked. Genes, the colour-by and the
// Tissue / Side by side / UMAP mode follow where the other dataset has them (each can be unticked).
//
// How the windows talk: a BroadcastChannel named after the group. The group's id is 128 random bits
// handed to the new window in window.name (never in the address), so groups cannot hear each other
// and nothing outside this origin can listen. A window only sends what ITS user changed: what it
// applies from another window becomes its "last seen" state, so nothing is ever sent back (no loops).
//
// Private data: a new window asks its group for the access token before it would show the password
// box (viewer.js: ensureToken). A member that holds a valid token answers, but only to a member of
// its group on the same page and the same token address. The token (read-only, short-lived) is
// handed over, never the password, and it is kept in memory only.
//
// Phones and tablets (a finger for a pointer), and browsers without BroadcastChannel: "Open in a new
// tab", no linking.
/* exported SodaCompare */
(function () {
  "use strict";
  const $ = (id) => document.getElementById(id);
  const MAX_WINDOWS = 4;              // in one group: each is a whole dataset in memory
  const BEAT_MS = 2000, GONE_MS = 7000;
  const TOKEN_WAIT_MS = 1200, TOKEN_QUIET_MS = 400;
  const POLL_MS = 200;
  const NAME_RE = /^soda-cmp:([0-9a-f]{32}):([0-9a-f]{8})$/;
  const KEY = "soda.compare.v1";
  const TOUCH = !!(window.matchMedia && matchMedia("(pointer: coarse)").matches);
  const CAN_LINK = typeof BroadcastChannel === "function" && !TOUCH;

  const hexRand = (bytes) => { const a = new Uint8Array(bytes); crypto.getRandomValues(a); return [...a].map((b) => b.toString(16).padStart(2, "0")).join(""); };
  const fin = (x) => typeof x === "number" && Number.isFinite(x);
  const me = hexRand(8);              // this page load
  let H = null;                       // what viewer.js hands over (init)
  let group = null, claim = "root", chan = null;
  let linked = false;
  let off = { ds: null, o: [0, 0], own: false };   // this window's offset (µm) for the dataset it was set on; own: set by "Link from here"
  const follow = { genes: true, colorBy: true, mode: true };
  const peers = new Map();            // id -> { ds, label, ready, spatial, um, d, seen }
  const expect = new Set();           // windows this one opened that have not opened their dataset yet
  const stats = { sent: {}, got: {} };
  let picker = null, note = "", missing = [];
  let openedDs = null;                // the dataset that has finished opening here (its gene list and layers are in)
  const isReady = () => { const S = H && H.state(); return !!(S && S.m && openedDs === S.dsId); };

  // ---- remembered for this tab (a reload keeps the group and the offset) ---------------------------
  function load() {
    const m = NAME_RE.exec(window.name || "");
    if (m) { group = m[1]; claim = m[2]; }
    let rec = null;
    try { rec = JSON.parse(sessionStorage.getItem(KEY) || "null"); } catch (e) { rec = null; }
    // a window opened by window.open starts with a COPY of its opener's storage: that record is the opener's
    if (!rec || typeof rec !== "object" || rec.claim !== claim || (group && rec.group !== group)) return;
    if (!group && /^[0-9a-f]{32}$/.test(rec.group || "")) group = rec.group;
    linked = !!rec.linked;
    if (rec.off && Array.isArray(rec.off.o) && fin(rec.off.o[0]) && fin(rec.off.o[1])) off = { ds: rec.off.ds || null, o: [rec.off.o[0], rec.off.o[1]], own: !!rec.off.own };
    if (rec.follow) for (const k of Object.keys(follow)) if (typeof rec.follow[k] === "boolean") follow[k] = rec.follow[k];
  }
  function save() {
    try { sessionStorage.setItem(KEY, JSON.stringify({ group, claim, linked, off, follow })); } catch (e) { /* storage blocked: a reload forgets the link */ }
  }

  // ---- the channel -------------------------------------------------------------------------------
  function join() {
    if (chan || !group || !CAN_LINK) return;
    chan = new BroadcastChannel(`soda-compare:${group}`);
    chan.onmessage = (e) => { try { onMessage(e.data); } catch (err) { /* a message this version does not understand */ } };
    post("hello", { info: info() });
  }
  function post(t, data) {
    if (!chan) return;
    stats.sent[t] = (stats.sent[t] || 0) + 1;
    chan.postMessage({ t, from: me, ...(data || {}) });
  }
  function info() {
    const S = H && H.state();
    if (!S || !S.m) return { ds: S ? S.dsId : null, ready: false, claim };
    const d = (S.datasets || []).find((x) => x.id === S.dsId);
    return { ds: S.dsId, label: String((d && d.label) || S.dsId).slice(0, 120), ready: isReady(), claim, spatial: S.mode !== "rna",
             um: S.m.pixelSizeUm > 0 ? S.m.pixelSizeUm : null, d: shiftNow() };
  }
  function seePeer(id, inf) {
    if (!inf || typeof inf !== "object") return null;
    const was = peers.get(id);
    const p = { ds: typeof inf.ds === "string" ? inf.ds : null, label: typeof inf.label === "string" ? inf.label.slice(0, 120) : "another dataset",
                ready: !!inf.ready, spatial: !!inf.spatial, um: fin(inf.um) ? inf.um : null, claim: typeof inf.claim === "string" ? inf.claim : "",
                d: Array.isArray(inf.d) && fin(inf.d[0]) && fin(inf.d[1]) ? [inf.d[0], inf.d[1]] : [0, 0], seen: Date.now() };
    peers.set(id, p);
    if (!was || was.ready !== p.ready || was.ds !== p.ds || was.d[0] !== p.d[0] || was.d[1] !== p.d[1]) paint();
    return p;
  }

  // ---- the camera in the shared frame ---------------------------------------------------------------
  const umOf = (S) => (S.m && S.m.pixelSizeUm > 0 ? S.m.pixelSizeUm : 1);      // no pixel size: one pixel counts as 1 µm
  function spatialNow() {
    const S = H && H.state();
    return !!(S && S.m && S.mode !== "rna" && S.umap && S.umap.view !== "umap" && !S.umap.anim && S.viewState && fin(S.viewState.zoom));
  }
  const centreUm = (S) => [S.m.image.width * umOf(S) / 2, S.m.image.height * umOf(S) / 2];
  // the offset belongs to one dataset: another dataset in this window starts from its own centre
  function offNow() {
    const S = H.state();
    if (off.ds !== S.dsId) off = { ds: S.dsId, o: centreUm(S), own: false };
    return off.o;
  }
  // how far this window's offset is from its picture's centre (µm): [0, 0] = centres lined up
  function shiftNow() {
    const S = H && H.state();
    if (!S || !S.m || S.mode === "rna" || !S.m.image) return [0, 0];
    const o = offNow(), c = centreUm(S);
    return [o[0] - c[0], o[1] - c[1]];
  }
  function sharedNow() {
    const S = H.state(), um = umOf(S), vs = S.viewState, o = offNow();
    return { c: [vs.target[0] * um - o[0], vs.target[1] * um - o[1]], s: vs.zoom - Math.log2(um) };
  }
  let lastCam = null;                 // the camera this window last sent or was given
  function applyShared(m) {
    if (!m || !Array.isArray(m.c) || !fin(m.c[0]) || !fin(m.c[1]) || !fin(m.s) || !spatialNow()) return false;
    const S = H.state(), um = umOf(S), o = offNow(), fit = H.imageFitVs();
    const zoom = Math.max(-30, Math.min(30, m.s + Math.log2(um)));
    const vs = { target: [(m.c[0] + o[0]) / um, (m.c[1] + o[1]) / um, 0], zoom, minZoom: Math.min(fit.minZoom, zoom), maxZoom: Math.max(fit.maxZoom, zoom) };
    H.setCamera(vs);
    lastCam = [vs.target[0], vs.target[1], vs.zoom];
    return true;
  }
  // once per frame: has this window's own user moved the picture?
  function tick() {
    requestAnimationFrame(tick);
    if (!H || !linked || !peers.size || !spatialNow()) { lastCam = null; return; }
    const vs = H.state().viewState, c = [vs.target[0], vs.target[1], vs.zoom];
    const first = !lastCam;
    if (!first && c[0] === lastCam[0] && c[1] === lastCam[1] && c[2] === lastCam[2]) return;
    lastCam = c;
    if (!first) post("view", sharedNow());             // the first look after linking only takes note
  }

  // ---- genes, colour-by, mode: changes made in this window ---------------------------------------------
  let applying = 0;                   // > 0 while a change from another window is being put in
  let ctxKey = null, lastGenes = new Map(), lastColorBy = null, lastMode = null;
  function genesNow() {
    const S = H.state();
    if (S.mode === "rna") return S.rna ? S.rna.genes : [];
    return S.genes || [];
  }
  function snap() {
    const S = H.state();
    ctxKey = `${S.dsId}|${S.mode}|${S.exprKey}`;
    lastGenes = new Map(genesNow().map((g) => [g.name, g.hex]));
    lastColorBy = S.mode !== "rna" && S.umap ? S.umap.colorBy : null;
    lastMode = S.mode !== "rna" && S.umap ? S.umap.view : null;
  }
  function poll() {
    if (!H) return;
    const S = H.state();
    if (!S.m) return;
    // another dataset or layer, or a change being applied: take note, send nothing
    if (applying || ctxKey !== `${S.dsId}|${S.mode}|${S.exprKey}` || !linked || !peers.size) { snap(); return; }
    const now = new Map(genesNow().map((g) => [g.name, g.hex]));
    for (const [name, hex] of now) {
      if (!lastGenes.has(name)) post("gene", { op: "add", name, hex });
      else if (lastGenes.get(name) !== hex) post("gene", { op: "col", name, hex });
    }
    for (const name of lastGenes.keys()) if (!now.has(name)) post("gene", { op: "rm", name });
    if (S.mode !== "rna" && S.umap) {
      if (S.umap.colorBy !== lastColorBy) post("colorBy", { v: S.umap.colorBy });
      if (S.umap.view !== lastMode) post("mode", { v: S.umap.view });
    }
    snap();
  }
  // everything this window shows that the others can follow: when linking turns on from here
  function pushAll() {
    if (!isReady()) return;
    const S = H.state();
    for (const g of genesNow()) post("gene", { op: "add", name: g.name, hex: g.hex });
    if (S.mode !== "rna" && S.umap) { post("colorBy", { v: S.umap.colorBy }); post("mode", { v: S.umap.view }); }
  }

  // ---- changes from another window ---------------------------------------------------------------------
  let queue = Promise.resolve();
  const later = (f) => { queue = queue.then(f).catch(() => {}); };
  async function during(f) { applying++; try { await f(); } finally { applying--; snap(); } }
  const sameName = (a, b) => String(a).toLowerCase() === String(b).toLowerCase();
  function miss(name, gone) {
    missing = missing.filter((n) => !sameName(n, name));
    if (!gone) missing.push(name);
    missing = missing.slice(-8);
    paint();
  }
  async function onGene(m) {
    if (!linked || !follow.genes || typeof m.name !== "string" || !m.name || m.name.length > 120) return;
    const S = H.state();
    if (!isReady()) return;                        // still opening: it asks for everything once it has opened
    const rna = S.mode === "rna", R = S.rna;
    if (rna && !R) return;
    const index = rna ? R.index : S.geneIndex, list = () => (rna ? R.genes : S.genes);
    const hex = typeof m.hex === "string" && /^#[0-9a-fA-F]{6}$/.test(m.hex) ? m.hex.toLowerCase() : null;
    const cur = list().find((g) => g.name === m.name) || list().find((g) => sameName(g.name, m.name));
    const redraw = () => { if (rna) { R.heat = null; H.buildRnaGeneUI(); H.renderRna(); } else { S.heat = null; H.buildGeneUI(); H.render(); } };
    const colour = (g) => { if (!hex || g.hex === hex) return; g.hex = hex; if (rna) H.recolorRnaGene(g); else H.recolorGene(g); };
    if (m.op === "rm") {
      miss(m.name, true);
      if (!cur) return;
      await during(() => { if (rna) { R.genes = R.genes.filter((g) => g !== cur); R.grid = ""; } else S.genes = S.genes.filter((g) => g !== cur); redraw(); });
      return;
    }
    if (cur) { await during(() => { colour(cur); redraw(); }); return; }
    if (m.op !== "add") return;
    // only a name this dataset has itself, and its own spelling of it
    const hit = (index || []).find((g) => g.n === m.name) || (index || []).find((g) => sameName(g.n, m.name));
    if (!hit) { miss(m.name, false); return; }
    await during(async () => {
      await (rna ? H.addRnaGene(hit.n) : H.addGene(hit.n));
      const got = list().find((g) => g.name === hit.n);
      if (got) { colour(got); redraw(); }
    });
  }
  async function onColorBy(m) {
    if (!linked || !follow.colorBy || !["gene", "cluster", "celltype"].includes(m.v)) return;
    const S = H.state(), el = $("colorBy");
    if (!isReady() || S.mode === "rna" || !S.umap || S.umap.colorBy === m.v || $("panelUmap").hidden) return;
    const opt = [...el.options].find((o) => o.value === m.v);
    if (!opt || opt.disabled || opt.hidden) { say(`This dataset cannot be coloured by ${m.v === "cluster" ? "clusters" : "suggested cell types"}: it stays as it is.`); return; }
    await during(async () => { el.value = m.v; await el.onchange({ target: el }); });
  }
  async function onMode(m) {
    if (!linked || !follow.mode || !["tissue", "sbs", "umap"].includes(m.v)) return;
    const S = H.state();
    if (!isReady() || S.mode === "rna" || !S.umap || S.umap.view === m.v) return;
    const sw = $("viewSwitch"), none = sw.hidden || [...sw.querySelectorAll(".seg")].some((e) => e.hidden);
    if (none) { if (m.v !== "tissue") say("This dataset has no UMAP: it stays on the tissue."); return; }
    await during(() => H.umapSetView(m.v));
  }

  // ---- messages -----------------------------------------------------------------------------------------
  function onMessage(m) {
    if (!m || typeof m !== "object" || typeof m.t !== "string" || typeof m.from !== "string" || m.from === me) return;
    stats.got[m.t] = (stats.got[m.t] || 0) + 1;
    const known = peers.get(m.from);
    if (known) known.seen = Date.now();
    switch (m.t) {
      case "hello":
        seePeer(m.from, m.info);
        post("here", { info: info(), linked, to: m.from });
        break;
      case "here":
        seePeer(m.from, m.info);
        // a window that has just come (or come back) takes the group's link state
        if (m.to === me && typeof m.linked === "boolean" && linked !== m.linked) { linked = m.linked; save(); paint(); if (linked) askSync(); }
        break;
      case "beat": {
        const p = seePeer(m.from, m.info);
        if (!known) post("here", { info: info(), linked, to: m.from });
        // the window this one opened has its dataset: line the two up and link them
        if (p && p.ready && expect.has(p.claim)) { expect.delete(p.claim); firstLink(); }
        break;
      }
      case "bye": {
        const p = peers.get(m.from);
        if (peers.delete(m.from)) { say(`${p && p.ready ? p.label : "The other window"} was closed.`); paint(); }
        break;
      }
      case "link":
        linked = !!m.on; save();
        if (linked && m.view) applyShared(m.view);
        paint();
        break;
      case "view":
        if (linked) applyShared(m);
        break;
      case "align":
        setOwnOffset(true);
        linked = true; save();
        if (fin(m.s)) applyShared({ c: [0, 0], s: m.s });       // stay where this window is, at the leader's scale
        paint();
        break;
      case "centre":
        setOwnOffset(false);
        linked = true; save();
        applyShared(m.view);
        paint();
        break;
      case "sync":
        if (linked && isReady()) { if (spatialNow()) post("view", sharedNow()); pushAll(); }
        break;
      case "gene": later(() => onGene(m)); break;
      case "colorBy": later(() => onColorBy(m)); break;
      case "mode": later(() => onMode(m)); break;
      case "token-ask": giveToken(m); break;
      case "token-give": takeToken(m); break;
      default: break;
    }
  }
  function askSync() { if (isReady()) post("sync"); }
  // from: true = what this window shows now ("Link from here"); false = the middle of its picture
  function setOwnOffset(here) {
    const S = H && H.state();
    if (!S || !S.m || S.mode === "rna" || !S.m.image) return;
    if (here && spatialNow()) { const um = umOf(S), vs = S.viewState; off = { ds: S.dsId, o: [vs.target[0] * um, vs.target[1] * um], own: true }; }
    else off = { ds: S.dsId, o: centreUm(S), own: false };
    lastCam = null;
    post("beat", { info: info() });                // the others show how the windows are shifted: tell them now
  }

  // ---- the three actions ------------------------------------------------------------------------------
  function setLinked(on) {
    linked = on; save();
    post("link", { on, view: on && spatialNow() ? sharedNow() : null });
    if (on) pushAll();
    lastCam = null;
    paint();
  }
  function linkFromHere() {
    if (!peers.size) return;
    const s = spatialNow() ? H.state().viewState.zoom - Math.log2(umOf(H.state())) : null;
    setOwnOffset(true);
    linked = true; save();
    post("align", { s });
    pushAll();
    say("Linked from here: the windows now move together from what each one shows.");
    paint();
  }
  function lineUpCentres(quiet) {
    if (!peers.size) return;
    setOwnOffset(false);
    linked = true; save();
    post("centre", { view: spatialNow() ? sharedNow() : null });
    pushAll();
    if (!quiet) say("Centres lined up: the middle of each picture is the same place, at the same scale.");
    paint();
  }
  function firstLink() {
    lineUpCentres(true);
    say("The two windows are linked: move one and the other follows. Link from here re-aligns them.");
  }

  // ---- opening the other window --------------------------------------------------------------------------
  function urlFor(id) {
    const u = new URL(location.href);
    u.hash = ""; u.searchParams.set("ds", id);
    for (const k of ["qc", "diag", "decode"]) u.searchParams.delete(k);      // the quality filter is this dataset's own
    return u.href;
  }
  function openBeside(id) {
    const S = H && H.state(), d = S ? (S.datasets || []).find((x) => x.id === id) : null, label = (d && d.label) || id;
    $("cmpBlocked").hidden = true;
    if (!CAN_LINK) {                                             // a phone, a tablet, an old browser: just another tab
      const w = window.open(urlFor(id), "_blank", "noopener");
      if (!w) blocked(urlFor(id), "_blank", label);
      return;
    }
    if (peers.size + 1 >= MAX_WINDOWS) { say(`Up to ${MAX_WINDOWS} windows can be compared at once: close one first.`); return; }
    if (!group) { group = hexRand(16); save(); }
    join();
    const nonce = hexRand(4), name = `soda-cmp:${group}:${nonce}`, url = urlFor(id);
    expect.add(nonce);
    const sw = screen.availWidth || 1600, sh = screen.availHeight || 900, sl = screen.availLeft || 0, st = screen.availTop || 0;
    const w = Math.max(520, Math.floor(sw / 2));
    let win = null;
    try { win = window.open(url, name, `popup=yes,width=${w},height=${sh},left=${sl + sw - w},top=${st}`); } catch (e) { win = null; }
    if (!win) { blocked(url, name, label); return; }
    say(`Opening ${label} in a new window. The two link up when it has loaded.`);
  }
  // the browser stopped the window: a plain link, which no blocker stops
  function blocked(url, name, label) {
    const a = $("cmpBlockedLink");
    a.href = url; a.target = name; a.rel = name === "_blank" ? "noopener" : "opener";
    a.textContent = `Open ${label} in a new tab`;
    $("cmpBlocked").hidden = false;
    a.focus();
  }

  // ---- private data: the access token, from a window of this group ------------------------------------------
  const tokenWait = new Map();        // k -> resolve
  function askToken(url) {
    if (!chan || typeof url !== "string" || !url) return Promise.resolve(null);
    return new Promise((resolve) => {
      const k = hexRand(8);
      let done = false;
      const end = (t) => { if (done) return; done = true; tokenWait.delete(k); resolve(t); };
      tokenWait.set(k, end);
      post("token-ask", { k, url, path: location.pathname });
      setTimeout(() => { if (!peers.size) end(null); }, TOKEN_QUIET_MS);       // nobody answered the hello: alone in the group
      setTimeout(() => end(null), TOKEN_WAIT_MS);
    });
  }
  function giveToken(m) {
    if (!H || !H.tokenUrl || !peers.has(m.from) || typeof m.k !== "string") return;
    if (m.url !== H.tokenUrl || m.path !== location.pathname) return;          // another SODA site on this origin: not ours to give
    const t = H.token();
    if (t) post("token-give", { to: m.from, k: m.k, token: t });
  }
  function takeToken(m) {
    const end = m.to === me && typeof m.k === "string" ? tokenWait.get(m.k) : null;
    const t = m.token;
    if (!end || !t || typeof t.sas !== "string" || typeof t.base !== "string" || !fin(t.expires) || t.expires < Date.now() + 60000) return;
    end({ sas: t.sas, base: t.base, expires: t.expires });
  }

  // ---- the section and the top-bar button -------------------------------------------------------------------
  function build() {
    const sec = document.createElement("section");
    sec.className = "panel"; sec.id = "panelCompare"; sec.dataset.tab = "save";
    sec.innerHTML =
      `<h2>Compare side by side</h2>` +
      `<div class="row cmp-pick"><span id="cmpOpenLabel">${CAN_LINK ? "Open beside" : "Open in a new tab"}</span><div id="cmpPicker"></div></div>` +
      `<div id="cmpBlocked" class="meta cmp-blocked" hidden>Your browser stopped the new window. <a id="cmpBlockedLink" href="#"></a></div>` +
      `<div id="cmpLinkBox"${CAN_LINK ? "" : " hidden"}>` +
        `<div id="cmpPeers" class="meta cmp-peers"></div>` +
        `<label class="row slim check"><input id="cmpLink" type="checkbox" disabled /><span>Link views</span></label>` +
        `<div class="btnrow">` +
          `<button id="cmpAlign" type="button" disabled title="Keep what each window shows now as the same place, and move them together from there">Link from here</button>` +
          `<button id="cmpCentre" type="button" disabled title="The middle of each picture is the same place, at the same scale">Line up centres</button></div>` +
        `<div id="cmpState" class="meta cmp-state"></div>` +
        `<div class="sub-h">Also follow</div>` +
        `<label class="row slim check"><input id="cmpGenes" type="checkbox" /><span>Genes picked</span></label>` +
        `<label class="row slim check"><input id="cmpColorBy" type="checkbox" /><span>Colour by</span></label>` +
        `<label class="row slim check"><input id="cmpMode" type="checkbox" /><span>Tissue, side by side or UMAP</span></label>` +
      `</div>` +
      `<div id="cmpNoLink" class="meta"${CAN_LINK ? " hidden" : ""}>The other dataset opens in a new tab. The two are not linked here: linked windows need a computer.</div>` +
      `<div id="cmpMsg" class="meta cmp-msg" role="status"></div>`;
    (document.getElementById("sidebar") || document.body).appendChild(sec);
    if (window.SodaLayout) window.SodaLayout.adopt(sec);
    // under Save and load and Save image (the layout puts a section it does not know last)
    const after = [$("panelImage"), $("panelSave")].find((x) => x && x.parentElement === sec.parentElement);
    if (after) after.after(sec);

    $("cmpLink").onchange = (e) => setLinked(e.target.checked);
    $("cmpAlign").onclick = linkFromHere;
    $("cmpCentre").onclick = () => lineUpCentres(false);
    for (const [id, k] of [["cmpGenes", "genes"], ["cmpColorBy", "colorBy"], ["cmpMode", "mode"]]) {
      $(id).checked = follow[k];
      $(id).onchange = (e) => { follow[k] = e.target.checked; save(); };
    }
    if (window.SodaDatasetPicker) {
      picker = window.SodaDatasetPicker.create({ host: $("cmpPicker"), labelledBy: "cmpOpenLabel", onPick: (id) => { resetPicker(); openBeside(id); } });
      resetPicker();
    }
    // the top bar: one button that brings this section up; it says so while windows are linked
    const bar = $("topBar"), fit = $("btnFit");
    if (bar) {
      const b = document.createElement("button");
      b.id = "btnCompare"; b.type = "button"; b.className = "cmp-top"; b.textContent = "Compare";
      b.title = "Open another dataset beside this one, and link the two views";
      b.onclick = () => {
        if (window.SodaLayout) window.SodaLayout.show("panelCompare");
        const pb = $("cmpPicker") && $("cmpPicker").querySelector(".dsp-btn");
        if (pb) requestAnimationFrame(() => pb.focus());
      };
      if (fit && fit.parentElement === bar) fit.before(b); else bar.appendChild(b);
    }
    paint();
  }
  function resetPicker() {
    if (!picker) return;
    picker.setValue(null);
    const v = $("cmpPicker").querySelector(".dsp-val");
    if (v) v.textContent = "Pick a dataset…";
  }
  let sayTimer = 0;
  function say(t) {
    note = t || "";
    clearTimeout(sayTimer);
    if (note) sayTimer = setTimeout(() => { note = ""; paint(); }, 12000);
    paint();
  }
  const set = (el, t) => { if (el && el.textContent !== t) el.textContent = t; };
  const um0 = (x) => `${Math.round(Math.abs(x)).toLocaleString()} µm`;
  function stateText() {
    const S = H && H.state();
    if (!peers.size) return "";
    if (!linked) return "Not linked: each window moves on its own.";
    if (!S || !S.m) return "";
    if (S.mode === "rna") return "Linked. A single-cell dataset has no slide, so only the genes follow here.";
    const parts = [];
    const others = [...peers.values()].filter((p) => p.ready && p.spatial);
    if (!others.length) parts.push("Linked. The other window has no slide to move with.");
    else if (S.umap && S.umap.view === "umap") parts.push("Linked. The UMAP view has no scale, so the windows move together again on the tissue.");
    else {
      const d = shiftNow();
      const rel = others.length === 1 ? [d[0] - others[0].d[0], d[1] - others[0].d[1]] : d;
      if (Math.abs(rel[0]) < 0.5 && Math.abs(rel[1]) < 0.5) parts.push("Linked at the same scale, centres lined up.");
      else parts.push(`Linked at the same scale. This window looks ${um0(rel[0])} further ${rel[0] >= 0 ? "right" : "left"} and ${um0(rel[1])} further ${rel[1] >= 0 ? "down" : "up"} ` +
        `than ${others.length === 1 ? "the other" : "the centre-to-centre place"}.`);
    }
    if (!(S.m.pixelSizeUm > 0)) parts.push("This dataset has no pixel size: one pixel is taken as 1 µm.");
    return parts.join(" ");
  }
  function summary() {
    if (!CAN_LINK) return "";
    if (!peers.size) return "";
    return `${linked ? "linked" : "not linked"} · ${peers.size + 1} windows`;
  }
  function paint() {
    if (!$("panelCompare")) return;
    const n = peers.size, ready = [...peers.values()].filter((p) => p.ready);
    set($("cmpPeers"), !n ? "No other window yet. Pick a dataset to open it beside this one."
      : `Beside this one: ${ready.length ? ready.map((p) => p.label).join(", ") : "a window that is still loading"}.`);
    const link = $("cmpLink");
    if (link.disabled !== !n) link.disabled = !n;
    if (link.checked !== (linked && n > 0)) link.checked = linked && n > 0;
    for (const id of ["cmpAlign", "cmpCentre"]) if ($(id).disabled !== !n) $(id).disabled = !n;
    set($("cmpState"), stateText());
    set($("cmpMsg"), [note, missing.length && linked ? `Not in this dataset: ${missing.join(", ")}.` : ""].filter(Boolean).join(" "));
    const b = $("btnCompare");
    if (b) {
      set(b, n && linked ? `Linked · ${n + 1}` : "Compare");
      b.classList.toggle("on", n > 0 && linked);
    }
  }

  // ---- start ------------------------------------------------------------------------------------------------
  load();
  join();                             // before viewer.js opens a dataset: a private one may ask this group for its token
  build();

  window.SodaCompare = {
    // viewer.js: its functions, once
    init(hooks) {
      H = hooks;
      window.addEventListener("soda:opened", () => {
        const S = H.state();
        openedDs = S.m ? S.dsId : null;
        if (picker && S.datasets) { picker.setItems(S.datasets, S.dataGroups || (window.SpatialViewerConfig || {}).groups || {}); resetPicker(); }
        missing = []; lastCam = null; snap();
        if (chan) { post("beat", { info: info() }); if (linked && peers.size) askSync(); }
        paint();
      });
      if (!CAN_LINK) return;
      requestAnimationFrame(tick);
      setInterval(poll, POLL_MS);
      setInterval(() => {
        if (!chan) return;
        post("beat", { info: info() });
        const now = Date.now();
        let gone = false;
        for (const [id, p] of peers) if (now - p.seen > GONE_MS) { peers.delete(id); gone = true; }
        if (gone) { say("A linked window stopped answering."); }
        paint();
      }, BEAT_MS);
      window.addEventListener("pagehide", () => post("bye"));
    },
    askToken,
    summary,
    // tests and the console
    state: () => ({ me, group, claim, linked, off: { ...off, o: [...off.o] }, follow: { ...follow }, canLink: CAN_LINK,
                    peers: [...peers.entries()].map(([id, p]) => ({ id, ...p })), expect: [...expect] }),
    stats: () => JSON.parse(JSON.stringify(stats)),
    resetStats: () => { stats.sent = {}; stats.got = {}; },
    open: openBeside, setLinked, linkFromHere, lineUpCentres,
    MAX_WINDOWS,
  };
})();
