// SpatialViewer - a thin viewing client over server.py in this folder.
// Images never come to the browser whole: each visible channel is a deck.gl
// TileLayer pulling 8-bit grey tiles, tinted + windowed on the GPU and
// composited with MAX blending like Xenium Explorer / Viv. Transcripts are per-gene
// point clouds fetched as raw float32 and drawn with ScatterplotLayer.
(function () {
  const cfg = window.SpatialViewerConfig || {};
  const params = new URLSearchParams(location.search);
  const API = (params.get("api") || cfg.apiBase || "").replace(/\/$/, "");
  const TILE_FMT = params.get("fmt") || cfg.tileFormat || "png";
  const GENE_COLORS = cfg.geneColors || ["#ff3355", "#33ddff", "#ffdd33", "#66ff66"];
  const NUCLEAR_CHANNEL = /^\s*(dapi|hoechst)\b/i;   // pinned to the top of the list
  const HIDE_CHANNEL = cfg.hideChannels ? new RegExp(cfg.hideChannels, "i")
                                        : /^\s*(open|blank|empty|unused|none)([\s_-]|$)/i;
  const $ = (id) => document.getElementById(id);

  if (!window.deck) {
    $("status").textContent = "deck.gl failed to load (no internet?)";
    $("status").hidden = false; $("status").classList.add("err");
    return;
  }
  const { Deck, OrthographicView, TileLayer, BitmapLayer, ScatterplotLayer, GeoJsonLayer, TextLayer, ColumnLayer, PolygonLayer } = deck;

  // ---- channel tile layer: grey texture -> window -> color, on the GPU ------
  const chanModule = {
    name: "chan",
    vs: "uniform chanUniforms { float lo; float hi; vec3 color; } chan;",
    fs: "uniform chanUniforms { float lo; float hi; vec3 color; } chan;",
    uniformTypes: { lo: "f32", hi: "f32", color: "vec3<f32>" },
  };
  class ChannelBitmapLayer extends BitmapLayer {
    getShaders() {
      const s = super.getShaders();
      return {
        ...s,
        modules: [...(s.modules || []), chanModule],
        inject: {
          ...(s.inject || {}),
          "fs:DECKGL_FILTER_COLOR": `
            float g = color.r;
            float v = clamp((g - chan.lo) / max(chan.hi - chan.lo, 1e-4), 0.0, 1.0);
            color = vec4(chan.color * v * layer.opacity, v * layer.opacity);`,
        },
      };
    }
    draw(opts) {
      const { lo, hi, color } = this.props;
      this.setShaderModuleProps({ chan: { lo, hi, color } });
      super.draw(opts);
    }
  }
  ChannelBitmapLayer.layerName = "ChannelBitmapLayer";
  ChannelBitmapLayer.defaultProps = {
    ...BitmapLayer.defaultProps,
    lo: { type: "number", value: 0 },
    hi: { type: "number", value: 1 },
    color: { type: "array", value: [1, 1, 1] },
  };
  // Additive sums whatever overlaps. That is why a coarse tile drawn under its own
  // children flashes bright mid-zoom, and it is also why two markers on the same spot
  // blow out. MAX takes the brighter of the two instead: overlapping copies of the same
  // tile are then a no-op, so the parent can stay on screen while children load - sharp
  // refinement and no flash at once.
  const MAXBLEND = {
    blend: true, depthTest: false,
    blendColorOperation: "max", blendColorSrcFactor: "one", blendColorDstFactor: "one",
    blendAlphaOperation: "max", blendAlphaSrcFactor: "one", blendAlphaDstFactor: "one",
  };
  const NORMAL = {
    blend: true, depthTest: false,
    blendColorOperation: "add", blendColorSrcFactor: "src-alpha", blendColorDstFactor: "one-minus-src-alpha",
    blendAlphaOperation: "add", blendAlphaSrcFactor: "one", blendAlphaDstFactor: "one-minus-src-alpha",
  };

  // ---- state ------------------------------------------------------------------
  const S = {
    datasets: [], dsId: null, m: null, src: null, mode: "image", rna: null,
    channels: [],          // {index,name,color:[r,g,b] 0-1, hex, lo, hi, visible}
    imgOpacity: 1, imgOn: true, geneOpacity: 1,
    exprKey: null, geneIndex: [], genes: [],   // genes: {name,hex,color,max,p99,n,positions,values,colors}
    geneBlend: "heatmap", pointScale: 1, genePalette: "genes", geneBin: 1, heat: null,
    tissue: { on: false, positions: null, key: null },
    // marker squares: the export's slide QC (qc.json) and per-square flags (qc_flags.bin)
    qc: { key: null, info: null, flags: null, positions: null, outline: false },
    overlays: [],          // {key,label,on,data}
    // Tier A region linking: meta from the manifest, ids aligned to positions.bin,
    // and a quantized-xy -> region-id lookup so a gene point finds its region.
    region: { meta: null, ids: null, lookup: null, selected: null, overlayKey: null },
    mol: molBlank(),       // molecules (Xenium ...): the gene list and the genes shown as dots
    lasso: null,           // circle -> UMAP: shapes, the bins under them, the last result
    umap: null,            // whole-section UMAP + Leiden: tissue <-> UMAP switch, clusters, selection
    viewState: null,
  };
  // A finger for a pointer: a phone or tablet. It has a fraction of a desktop's memory and
  // single-core speed, and iOS reloads a tab that outgrows ~1-1.5 GB, so a few limits are
  // tighter on it. Decided once, at load.
  const TOUCH = !!(window.matchMedia && matchMedia("(pointer: coarse)").matches);
  // A phone screen is 3x its CSS size: drawing every one of those pixels is 9x the work for
  // detail nobody can see at arm's length. 2x is as sharp to the eye.
  const DEVICE_PX = TOUCH ? Math.min(window.devicePixelRatio || 1, 2) : true;
  const tileCache = new Map();   // "ds/c/z/x/y" -> ImageBitmap (LRU)
  // On touch the decoded tiles kept for a quick return get ~256 MB (4 bytes a pixel), not 1500 tiles.
  // config.js's tileCacheSize is a desktop number (the example ships 400 = 1.6 GB of 1024 px
  // tiles), so on touch it can only lower that limit, never raise it.
  const tileCacheMax = () => {
    if (!(TOUCH && S.m)) return cfg.tileCacheSize || 1500;
    const touchMax = Math.max(32, Math.floor(256e6 / (4 * S.m.image.tileSize * S.m.image.tileSize)));
    return Math.min(cfg.tileCacheSize || touchMax, touchMax);
  };
  // On touch each stain's tile layer keeps at most this many tiles (GPU textures) between them,
  // rather than deck's five screens' worth per stain
  const TOUCH_TILE_BUDGET = 96;
  const TOUCH_MANY_CHANNELS = 10;      // more stains than this on at once: say a phone may reload
  // How many genes may be shown at once, so no dataset can exhaust the browser tab.
  // UMAP panels are light (one small panel per gene); a spatial transcript layer draws every
  // bin of that gene, so it is capped more tightly.
  const MAX_RNA_GENES = Math.max(1, Math.min(cfg.rnaseqMaxGenes || 12, 24));
  const MAX_SPATIAL_GENES = Math.max(1, Math.min(cfg.maxGenes || 8, 16));
  let deckgl = null;

  // ---- helpers ----------------------------------------------------------------
  const hexToRgb01 = (h) => [1, 3, 5].map((i) => parseInt(h.slice(i, i + 2), 16) / 255);
  const hexToRgb255 = (h) => [1, 3, 5].map((i) => parseInt(h.slice(i, i + 2), 16));
  // A problem that stays true (no WebGL2, graphics memory taken back, tiles that fail) is
  // kept in stickyErr: the many status(null) calls that end a "Loading…" fall back to it
  // instead of wiping the one line that says why the map is black.
  let stickyErr = null;
  function status(msg, err) {
    const el = $("status");
    if (!msg && stickyErr) { msg = stickyErr; err = true; }
    if (!msg) { el.hidden = true; return; }
    el.textContent = msg; el.hidden = false; el.classList.toggle("err", !!err);
  }
  function stickyStatus(msg, ms) {
    stickyErr = msg; status(msg, true);
    if (ms) setTimeout(() => { if (stickyErr === msg) { stickyErr = null; if ($("status").textContent === msg) status(null); } }, ms);
  }
  const noQuery = (t) => String(t).replace(/\?[^\s)"']*/g, "");   // never show a token: drop "?..." from anything printed
  function apiUrl(p) { return `${API}/api/${p}`; }

  // ---- data sources -------------------------------------------------------------
  // "static": a folder of files exported by prepare.py export, hosted anywhere.
  //           Normally a PRIVATE blob container: `qs` appends the access token
  //           from /api/token to every request (see below).
  // "api":    the tile server (server.py), password protected.
  function staticSource(base, qs) {
    base = base.replace(/\/?$/, "/");
    const u = (p) => base + p + (qs ? qs() : "");
    return {
      kind: "static",
      manifest: () => getJSON(u("manifest.json"), true),
      tileUrl: (c, z, x, y, fmt) => u(`tiles/${c}/${z}/${x}/${y}.${fmt}`),
      genes: (key) => getJSON(u(`expr/${key}/genes.json`), true),
      // The export names a gene file encodeURIComponent(name) + ".bin", so for a name with a
      // character that gets encoded ("TRBC1/2", "KIR Family": GeoMx panels have them) the file
      // name itself holds a "%", and the URL must encode it again. Names of letters and digits
      // are unchanged by either encoding. The singly encoded URL is tried second, for a host
      // whose upload decoded the names.
      gene: async (key, name) => {
        const e1 = encodeURIComponent(name), e2 = encodeURIComponent(e1);
        if (e2 === e1) return getF32(u(`expr/${key}/gene/${e1}.bin`), true);
        try { return await getF32(u(`expr/${key}/gene/${e2}.bin`), true); }
        catch (err) { return getF32(u(`expr/${key}/gene/${e1}.bin`), true); }
      },
      positions: (key) => getF32(u(`expr/${key}/positions.bin`), true),
      region: (key) => getU16(u(`expr/${key}/region.bin`), true),
      overlay: (file) => getJSON(u(`overlays/${file}`), true),
      json: (p) => getJSON(u(p), true),
      buf: async (p) => { const r = await dataFetch(u(p)); if (!r.ok) throw new Error(`${r.status} ${p}`); return r.arrayBuffer(); },
      mol: async (p) => { const r = await dataFetch(u(`molecules/${p}`)); if (!r.ok) throw new Error(`${r.status} molecules/${p}`); return r.arrayBuffer(); },
    };
  }
  function apiSource(dsId) {
    return {
      kind: "api",
      manifest: () => getJSON(apiUrl(`datasets/${dsId}`)),
      tileUrl: (c, z, x, y, fmt) => apiUrl(`datasets/${dsId}/tiles/${c}/${z}/${x}/${y}.${fmt}`),
      genes: (key) => getJSON(apiUrl(`datasets/${dsId}/expr/${key}/genes`)),
      gene: (key, name) => getF32(apiUrl(`datasets/${dsId}/expr/${key}/gene/${encodeURIComponent(name)}`)),
      positions: (key) => getF32(apiUrl(`datasets/${dsId}/expr/${key}/positions`)),
      region: (key) => getU16(apiUrl(`datasets/${dsId}/expr/${key}/region`)),
      overlay: (key) => getJSON(apiUrl(`datasets/${dsId}/overlays/${key}`)),
      mol: async (p) => { const r = await apiFetch(apiUrl(`datasets/${dsId}/molecules/${p}`)); if (!r.ok) throw new Error(`${r.status} molecules/${p}`); return r.arrayBuffer(); },
    };
  }
  // ---- one dataset at a time, one layer at a time ---------------------------------
  // Picking a second dataset while the first is still loading used to let the first one's
  // late answers land on the second. openSeq numbers the opens and layerSeq the layer
  // switches: code that awaits keeps the number it started under and stops when it is no
  // longer the current one. On top of that every request of a dataset goes through its
  // source, and opening another dataset retires the source: its fetches are aborted and
  // whatever was still in flight never settles (resolve or reject), so no loader of the old
  // dataset can run on into the new one or report the abort as an error. Anything waiting
  // on a retired source is simply dropped and garbage-collected.
  let openSeq = 0, layerSeq = 0, openCtl = null, srcSignal = null;
  function liveSource(src, seq, signal) {
    const live = () => seq === openSeq;
    const dropped = () => new Promise(() => {});           // a fresh one each time: nothing holds on to it
    const hold = (p) => Promise.resolve(p).then((v) => (live() ? v : dropped()),
                                                (e) => { if (live()) throw e; return dropped(); });
    for (const k of Object.keys(src)) {
      const f = src[k];
      if (typeof f !== "function") continue;
      // dataFetch / apiFetch pick srcSignal up as they start, so every method gets the signal.
      // A method that answers at once (tileUrl builds a string) is passed through as it is.
      src[k] = (...a) => {
        srcSignal = signal;
        try { const r = f(...a); return r && typeof r.then === "function" ? hold(r) : r; } finally { srcSignal = null; }
      };
    }
    return src;
  }
  // ---- private container: short-lived read-only access token --------------------
  // The exported datasets live in a private blob container, so a URL on its own
  // fetches nothing. /api/token checks the shared password on the server and
  // returns a read-only SAS that expires after TOKEN_HOURS (default 1, at most 24; a preset
  // BLOB_SAS keeps its own expiry; see deploy/azure/api/token/index.js); the storage key never
  // reaches the browser. A session cookie set by that call lets an open tab renew
  // silently, so nobody has to retype the password in the middle of a session.
  const TOKEN_URL = cfg.tokenEndpoint || "";
  // ---- ?data=<folder link>: someone else's export folder, no config.js needed ----------------
  // The page reads <folder>/datasets.json (written by catalog.py) and lists every dataset in
  // it; each one is read from <folder>/<id>/. A query string on the link (a read-only
  // pass) is added to every file. Off on a password-protected site unless config.js sets
  // allowDataLinks: a link built by anyone could otherwise point a signed-in reader's page
  // at a folder they control.
  // A copy of SODA on someone's own site can name its own folder in config.js (`data:`); the
  // page then opens it with no ?data= on the link. That one is the site owner's choice, so it
  // is allowed even where links from outside are not. A ?data= on the link still wins.
  const DATA_LINKS = cfg.allowDataLinks !== undefined ? !!cfg.allowDataLinks : !TOKEN_URL;
  const DATA = (() => {
    const fromLink = (params.get("data") || "").trim();
    const raw = fromLink || (typeof cfg.data === "string" ? cfg.data.trim() : "");
    if (!raw) return null;
    if (fromLink && !DATA_LINKS) return { error: "this site does not open data links (it is password-protected)" };
    let u;
    try { u = new URL(raw, location.href); } catch (e) { return { error: `"${raw}" is not a web address` }; }
    if (!/^https?:$/.test(u.protocol)) return { error: "the data link must start with https://" };
    const qs = u.search;
    u.search = ""; u.hash = "";
    return { root: u.href.replace(/\/?$/, "/"), qs, fromLink: !!fromLink };
  })();
  // Settings a data folder may carry in its datasets.json. The folder can be anyone's, so only
  // these, each checked for its type: nothing in it can name a script, a server or a password.
  const strList = (v) => (Array.isArray(v) ? v.filter((x) => typeof x === "string").slice(0, 1000) : undefined);
  const plainObj = (v) => (v && typeof v === "object" && !Array.isArray(v) ? v : undefined);
  function folderSettings(d) {
    const s = {};
    if (strList(d.hideChannels)) s.hideChannels = strList(d.hideChannels);
    if (strList(d.defaultChannels)) s.defaultChannels = strList(d.defaultChannels);
    if (plainObj(d.labelPositions)) s.labelPositions = plainObj(d.labelPositions);
    return s;
  }
  // `tags` ({tissue: "nerve", ...}) for the picker's filters: short strings only, at most 20.
  const folderTags = (v) => (plainObj(v) ? Object.fromEntries(Object.entries(v)
    .filter(([k, x]) => typeof x === "string" && k.length <= 40 && x.length <= 80).slice(0, 20)) : undefined);
  // Per-dataset settings: from the data folder for a ?data= dataset, else config.js.
  const dsSettings = (id) => {
    const d = (S.datasets || []).find((x) => x.id === id);
    return (d && d.settings) || (cfg.datasets && cfg.datasets[id]) || {};
  };
  let token = null, tokenWait = null;
  const tokenQuery = () => (token ? "?" + token.sas : "");
  async function requestToken(password) {
    const opts = password == null
      ? { method: "GET", cache: "no-store", credentials: "same-origin" }
      : { method: "POST", cache: "no-store", credentials: "same-origin",
          headers: { "Content-Type": "application/json" }, body: JSON.stringify({ password }) };
    const r = await fetch(TOKEN_URL, opts);
    if (r.status === 401) return null;                     // no session yet, or wrong password
    if (r.status === 429) throw new Error("too many password attempts - wait a few minutes");
    if (!r.ok) throw new Error(`${r.status} from ${TOKEN_URL}`);
    return r.json();
  }
  async function ensureToken(force) {
    if (!force && token && Date.now() < token.expires - 60000) return token;
    if (tokenWait) return tokenWait;                       // one request, however many tiles wait on it
    tokenWait = (async () => {
      try {
        let t = await requestToken(null);                  // cookie still valid? then no prompt
        for (let tries = 0; !t && tries < 6; tries++) {
          const pw = await askPassword(tries > 0, "This viewer is password protected.");
          t = await requestToken(pw);
        }
        if (!t) throw new Error("password not accepted");
        token = t;
        return t;
      } finally { tokenWait = null; }
    })();
    return tokenWait;
  }
  // A URL built with a token that expired while it was queued: sign it again.
  const reSign = (url) => (token ? url.split("?")[0] + tokenQuery() : url);
  async function dataFetch(url, opts = {}) {
    if (srcSignal && !opts.signal) opts = { ...opts, signal: srcSignal };   // see liveSource
    let r = await fetch(url, opts);
    if ((r.status === 403 || r.status === 401) && TOKEN_URL && token) {
      await ensureToken(true);
      r = await fetch(reSign(url), opts);
    }
    return r;
  }
  // Renew a little before expiry so a long-open tab never sees a failed tile.
  if (TOKEN_URL) setInterval(() => {
    if (token && Date.now() > token.expires - 5 * 60000) ensureToken(true).catch(() => {});
  }, 60000);

  async function resolveLink(path) {
    // "<id>.link.txt" holds the real URL,
    // either plain or encrypted with the shared password ("enc:v1:...", see links.py).
    // Only used for datasets configured with `linkPath`: the no-server password option
    // ("Encrypted link files" in deploy/azure/README.md). Production uses /api/token instead.
    if (!/\.link\.txt$/i.test(path)) return path;
    const r = await fetch(path, { cache: "no-cache" });
    if (!r.ok) throw new Error(`${r.status} ${path}`);
    const line = (await r.text()).trim().split(/\r?\n/)[0].trim();
    if (!line.startsWith("enc:v1:")) return line;
    return decryptLink(line);
  }
  // ---- encrypted links: PBKDF2-SHA256 (200k) -> AES-GCM, all in the browser ----
  let linkPassword = null;
  try { linkPassword = sessionStorage.getItem("spatialviewer.pw"); } catch (e) { /* ignore */ }
  const b64 = (str) => Uint8Array.from(atob(str), (ch) => ch.charCodeAt(0));
  async function tryDecrypt(payload, password) {
    const [, , salt, iv, ct] = payload.split(":");
    const base = await crypto.subtle.importKey("raw", new TextEncoder().encode(password), "PBKDF2", false, ["deriveKey"]);
    const key = await crypto.subtle.deriveKey({ name: "PBKDF2", salt: b64(salt), iterations: 200000, hash: "SHA-256" },
      base, { name: "AES-GCM", length: 256 }, false, ["decrypt"]);
    const plain = await crypto.subtle.decrypt({ name: "AES-GCM", iv: b64(iv) }, key, b64(ct));
    return new TextDecoder().decode(plain);
  }
  async function decryptLink(payload) {
    if (!window.isSecureContext || !crypto.subtle) throw new Error("this page needs https (or localhost) to decrypt links");
    for (let tries = 0; tries < 6; tries++) {
      if (linkPassword) {
        try { return await tryDecrypt(payload, linkPassword); } catch (e) { linkPassword = null; }
      }
      const pw = await askPassword(tries > 0, "This dataset is password protected.");
      linkPassword = pw;
      try { sessionStorage.setItem("spatialviewer.pw", pw); } catch (e) { /* ignore */ }
    }
    throw new Error("password not accepted");
  }

  // ---- password (HTTP Basic) ----------------------------------------------------
  // Same-origin: the browser already asked once when the page loaded and sends
  // the credentials itself. Cross-origin (?api=...): browsers never show the
  // login box for fetch(), so we ask here and send the header ourselves.
  let authHeader = null;
  try { const s = sessionStorage.getItem("spatialviewer.auth"); if (s) authHeader = s; } catch (e) { /* ignore */ }
  // One inline login box shared by every request that hits a 401; concurrent
  // requests (tiles) all wait on the same submission.
  let loginWait = null, loginResolve = null;
  // The form must never actually submit (that would put the password in the URL).
  $("login").addEventListener("submit", (e) => {
    e.preventDefault();
    const pw = $("loginPw").value;
    if (!pw) return;
    authHeader = "Basic " + btoa("viewer:" + pw);
    try { sessionStorage.setItem("spatialviewer.auth", authHeader); } catch (err) { /* ignore */ }
    $("login").hidden = true;
    const r = loginResolve; loginWait = null; loginResolve = null;
    if (r) r(pw);
  });
  function askPassword(wrong, text) {
    if (loginWait) return loginWait;
    const form = $("login"), msg = $("loginMsg"), pw = $("loginPw");
    msg.textContent = wrong ? "Wrong password - try again." : (text || `This server (${API || location.origin}) needs a password.`);
    msg.classList.toggle("err", !!wrong);
    form.hidden = false; pw.value = ""; pw.focus();
    loginWait = new Promise((resolve) => { loginResolve = resolve; });
    return loginWait;
  }
  async function apiFetch(url, opts = {}) {
    if (srcSignal && !opts.signal) opts = { ...opts, signal: srcSignal };   // see liveSource
    let r = null;
    for (let tries = 0; tries < 6; tries++) {
      const sent = authHeader;
      r = await fetch(url, { ...opts, headers: { ...(opts.headers || {}), ...(sent ? { Authorization: sent } : {}) } });
      if (r.status !== 401) return r;
      if (authHeader !== sent) continue;           // a password arrived while this was in flight: just retry
      if (sent) {                                  // this request carried the current password, so it is wrong
        authHeader = null;
        try { sessionStorage.removeItem("spatialviewer.auth"); } catch (e) { /* ignore */ }
      }
      await askPassword(!!sent);
    }
    return r;
  }
  async function getJSON(url, plain) {
    const r = plain ? await dataFetch(url) : await apiFetch(url);
    if (!r.ok) throw new Error(`${r.status} ${url}`);
    const t = await r.text();
    try { return JSON.parse(t); } catch (e) {
      // Python's json writes NaN, which JSON.parse refuses (umap/index.json "bounds" of
      // exports made before that was fixed): read those as null rather than lose the file
      if (/\bNaN\b/.test(t)) return JSON.parse(t.replace(/\bNaN\b/g, "null"));
      throw e;
    }
  }
  async function getF32(url, plain) {
    const r = plain ? await dataFetch(url) : await apiFetch(url);
    if (!r.ok) throw new Error(`${r.status} ${url}`);
    return new Float32Array(await r.arrayBuffer());
  }
  async function getU16(url, plain) {
    const r = plain ? await dataFetch(url) : await apiFetch(url);
    if (!r.ok) throw new Error(`${r.status} ${url}`);
    return new Uint16Array(await r.arrayBuffer());
  }
  function percentile(sortedCopy, q) {
    if (!sortedCopy.length) return 1;
    return sortedCopy[Math.min(sortedCopy.length - 1, Math.floor(q * sortedCopy.length))];
  }

  // ---- layers -----------------------------------------------------------------
  function tileUrl(c, z, x, y) {
    const fmt = S.src.kind === "static" ? (S.m.tileFormat || "png") : TILE_FMT;
    return S.src.tileUrl(c, z, x, y, fmt);
  }
  // ---- tiles: decode with a fallback, and say so when one fails ----------------------
  // A tile is decoded with createImageBitmap, as before. Where that is missing or refuses
  // the tile, the same bytes are tried without the options and then through an <img>
  // element, and the first way that works is kept for the rest of the visit (so only the
  // first tile pays for the failed attempt). ?decode=plain or ?decode=img starts at that
  // step, to compare them on a phone.
  const TILE_DECODERS = [
    ["createImageBitmap", (b) => createImageBitmap(b, { premultiplyAlpha: "none", colorSpaceConversion: "none" })],
    ["createImageBitmap without options", (b) => createImageBitmap(b)],
    ["an img element", (b) => new Promise((resolve, reject) => {
      const img = new Image(), u = URL.createObjectURL(b);
      img.onload = () => { URL.revokeObjectURL(u); resolve(img); };
      img.onerror = () => { URL.revokeObjectURL(u); reject(new Error("the picture could not be decoded")); };
      img.src = u;
    })],
  ];
  let tileDecoder = typeof createImageBitmap !== "function" || params.get("decode") === "img" ? 2
    : params.get("decode") === "plain" ? 1 : 0;
  async function decodeTile(blob) {
    let first = null;
    for (let i = tileDecoder; i < TILE_DECODERS.length; i++) {
      try {
        const img = await TILE_DECODERS[i][1](blob);
        if (!img || !img.width) throw new Error("the picture came back empty");
        if (i > tileDecoder) {
          diagNote("decode", `${TILE_DECODERS[tileDecoder][0]} failed (${first && first.message}); using ${TILE_DECODERS[i][0]}`);
          tileDecoder = i;
        }
        return img;
      } catch (e) { first = first || e; }
    }
    throw Object.assign(new Error(`this browser could not decode an image tile (${first && first.message})`), { tileDecode: true });
  }
  // Failures used to be dropped without a word: a phone showed a black map and no reason.
  // The first one of a dataset goes to the status line in plain words; every one is kept
  // for the ?diag=1 panel (diag.js), which never prints the part of an address after "?".
  const diagNote = (kind, msg) => { try { window.__sodaDiag.note(kind, msg); } catch (e) { /* diag.js not loaded */ } };
  let tileProblems = 0, tileMsg = null;
  function tileProblem(what, where) {
    diagNote("tile", `${what}${where ? " - " + String(where).split("?")[0] : ""}`);
    tilesFailed = true; retryTilesSoon();
    if (tileProblems++) return;
    const again = !!tileRetry;                 // no retry left: the line stays until the page is reloaded or another dataset opens
    tileMsg = `Part of the image could not be shown: ${noQuery(what)}. ` +
      (again ? "Trying again…" : "Reload the page; if it stays, add &diag=1 to the address for details.");
    stickyStatus(tileMsg, again ? 30000 : 0);
  }
  // A tile that failed for a real reason stays blank in deck until it is asked for again, so
  // one dropped connection (a phone changing networks, or coming back from the lock screen)
  // left holes for the rest of the visit. Ask again: a few times with a growing wait, and
  // whenever the page comes back to the front or the network returns. tileEpoch is in the
  // tile layers' updateTriggers, so raising it makes deck reload its tiles (the ones that
  // did arrive come straight from tileCache).
  let tileEpoch = 0, tileRetry = null, tileRetries = 0, tilesFailed = false;
  function retryTilesSoon(ms) {
    if (tileRetry || !tilesFailed) return;
    if (ms == null) { if (tileRetries >= 4) return; ms = 2000 * Math.pow(3, tileRetries++); }   // 2 s, 6 s, 18 s, 54 s
    tileRetry = setTimeout(() => {
      tileRetry = null; tilesFailed = false; tileProblems = 0; tileEpoch++;
      if (stickyErr && stickyErr === tileMsg) { stickyErr = null; if ($("status").textContent === tileMsg) status(null); }
      render();
    }, ms);
  }
  const retryTilesNow = () => { if (!tilesFailed) return; clearTimeout(tileRetry); tileRetry = null; tileRetries = 0; retryTilesSoon(1000); };
  document.addEventListener("visibilitychange", () => { if (!document.hidden) retryTilesNow(); });
  window.addEventListener("online", retryTilesNow);
  // The map itself cannot draw: no WebGL2, a layer deck could not build, or the browser
  // took the graphics memory back (phones do when memory runs short; luma.gl cannot rebuild
  // its textures on a restored context, so the way back is a reload).
  let drawProblems = 0, noGraphics = false;
  function drawProblem(err, layer) {
    const msg = noQuery((err && err.message) || err).slice(0, 160);
    console.error(err);
    diagNote("draw", `${layer && layer.id ? layer.id + ": " : ""}${msg}`);
    if (drawProblems++ || noGraphics) return;        // without graphics the line already says why
    stickyStatus(`The map could not be drawn (${msg}). Reload the page; if it stays, add &diag=1 to the address for details.`, 20000);
  }
  function checkWebGL2() {
    let gl2 = null;
    try { gl2 = document.createElement("canvas").getContext("webgl2"); } catch (e) { /* treated as missing */ }
    if (gl2) { const lose = gl2.getExtension("WEBGL_lose_context"); if (lose) lose.loseContext(); return; }
    noGraphics = true;
    diagNote("draw", "no WebGL2 context");
    stickyStatus("This browser cannot draw the map: WebGL2 (3D graphics) is off or missing. On an iPhone or iPad: update to iOS 15 or newer, " +
      "and if Lockdown Mode is on, turn it off for this site (Safari: Website Settings for this page).");
  }
  function watchContext(device) {
    const canvas = device && device.canvasContext && device.canvasContext.canvas;
    if (!canvas || !canvas.addEventListener) return;
    canvas.addEventListener("webglcontextlost", () => {
      noGraphics = true;
      diagNote("draw", "the graphics context was lost (webglcontextlost)");
      tileCache.clear();
      stickyStatus("The browser took back the map's graphics memory (a phone does this when memory runs short). " +
        "Reload the page to draw it again; fewer stains on at once keeps it from happening.");
    });
  }
  function fetchTile(c, backendZ, x, y, signal) {
    const key = `${S.dsId}/${c}/${backendZ}/${x}/${y}`;
    const hit = tileCache.get(key);
    if (hit) { tileCache.delete(key); tileCache.set(key, hit); return Promise.resolve(hit); }
    const url = tileUrl(c, backendZ, x, y);
    const req = S.src.kind === "static" ? dataFetch(url, { signal }) : apiFetch(url, { signal });
    return req
      .then((r) => {
        if (r.ok && r.status !== 204) return r.blob();
        if (!r.ok) tileProblem(r.status === 401 || r.status === 403 ? `the data store refused the image tiles (${r.status}; the sign-in may have run out)`
          : r.status === 404 ? "image tiles are missing on the server (404)" : `the server answered ${r.status} for image tiles`, url);
        return null;
      })
      .then((b) => (b ? decodeTile(b) : null))
      .then((bmp) => {
        if (bmp) {
          tileCache.set(key, bmp);
          // never close() evicted bitmaps: deck may still be drawing them; GC reclaims them
          while (tileCache.size > tileCacheMax()) tileCache.delete(tileCache.keys().next().value);
        }
        return bmp;
      })
      .catch((e) => {
        // an abort is deck moving on, not a failure (the caller re-throws it to deck)
        if (!(e && (e.name === "AbortError" || (signal && signal.aborted))))
          tileProblem(e && e.tileDecode ? e.message : `image tiles could not be downloaded (${(e && e.message) || e}; the network, or a content blocker)`, url);
        throw e;
      });
  }

  function channelLayer(ch, op = S.imgOpacity, nOn = 1) {
    const img = S.m.image;
    const maxLevel = img.levels - 1;
    return new TileLayer({
      id: `ch-${S.dsId}-${ch.index}`,
      tileSize: img.tileSize,
      minZoom: -maxLevel,
      maxZoom: 0,
      extent: [0, 0, img.width, img.height],
      maxRequests: 12,
      // touch: the stains share one budget, but each keeps enough for a phone screen and its parents
      ...(TOUCH ? { maxCacheSize: Math.max(8, Math.floor(TOUCH_TILE_BUDGET / nOn)) } : {}),
      // Safe under MAX: an overlapping parent tile cannot brighten anything, so a coarse
      // tile can show immediately while its children load.
      refinementStrategy: "best-available",
      // renderSubLayers closes over ch.lo / ch.hi / ch.hex and the global image
      // settings. Without them here deck keeps the sublayers it already built, so
      // moving a window slider changed nothing on screen until you panned.
      updateTriggers: {
        getTileData: [S.dsId, tileEpoch],
        renderSubLayers: [ch.lo, ch.hi, ch.hex, op],
      },
      getTileData: ({ index, signal }) => {
        const z = maxLevel + index.z;
        if (z < 0 || z > maxLevel) return Promise.resolve(null);
        return fetchTile(ch.index, z, index.x, index.y, signal).catch((e) => {
          // An abort is "we moved on", not "this tile is empty". Swallowing it as
          // null marked the tile loaded-with-no-content, so deck stopped drawing the
          // parent tile in its place and never asked for it again: holes that only a
          // pan would fill, and detail that never sharpened past the coarse level.
          // Re-throw so deck discards the tile and re-requests it.
          if (e && (e.name === "AbortError" || signal?.aborted)) throw e;
          return null;                                  // a real error: leave it blank
        });
      },
      onTileError: () => {},
      // Force a draw once the visible tiles are in. On the big slides (8 levels, first
      // open, cold cache) the sublayers were built and loaded but nothing was painted
      // until something else forced a redraw - the dataset opened black until a manual
      // Fit. Costs one redraw per viewport load.
      onViewportLoad: () => { if (deckgl) deckgl.redraw(true); },
      renderSubLayers: (props) => {
        if (!props.data) return null;
        const bb = props.tile.boundingBox;
        return new ChannelBitmapLayer(props, {
          data: null,
          image: props.data,
          bounds: [bb[0][0], bb[1][1], bb[1][0], bb[0][1]],
          lo: ch.lo / 255, hi: ch.hi / 255, color: ch.color,
          opacity: op,
          parameters: MAXBLEND,
          // Always interpolated. This was a checkbox, but BitmapLayer declares
          // textureParameters as `ignore: true`, so the sampler is fixed when the
          // texture is built and nothing revisits it - the control could only ever
          // be made to work by rebuilding every tile layer to change one sampler.
          textureParameters: { minFilter: "linear", magFilter: "linear" },
        });
      },
    });
  }

  // ---- transcript heat map, the way FIVE Atlas draws one --------------------------
  // One square per bin, colored by the COMBINED density of every gene on show,
  // through a single scientific ramp. Two differences from coloring each gene on
  // its own ramp, which is what this used to do:
  //   - genes add up instead of hiding one another. With a layer per gene the last
  //     one drawn simply covered the others, so a second gene could make the first
  //     disappear rather than reading as more signal.
  //   - a bin holding transcripts never falls off the bottom of the ramp. The
  //     window's low end lifts to 6% rather than 0, so "a little here" stays
  //     visible and only genuinely empty bins are blank.
  // Squares, not discs, and sized to the bin pitch: the bins tile the slide, so
  // circles left gaps that read as structure that is not there.
  function geneHeat() {
    const e = currentExpr();
    const pitch = (e && e.spotDiameterPx) ? e.spotDiameterPx : 8;
    const cell = pitch * Math.max(1, S.geneBin);
    const ink = S.genePalette === "genes";
    const on = S.genes.filter((g) => g.visible !== false);
    // Tier A: when a region is selected, the density collapses to its bins only.
    const sel = S.region.selected || 0;
    const inSel = (g, i) => !sel || regionAt(g.positions[2 * i], g.positions[2 * i + 1]) === sel;
    // protein squares are windowed by both handles (like a stain), so the top one counts too
    const prot = isProt();
    const key = on.map((g) => `${g.name}:${g.hex}:${g.min}` + (prot ? `:${g.max}` : "")).join("|") +
                `#${S.genePalette}#${S.geneBin}#${S.exprKey}#${cell}#r${sel}`;
    if (S.heat && S.heat.key === key) return S.heat;

    const idx = new Map(), KX = 1 << 21;
    const binOf = (g, i) => Math.floor(g.positions[2 * i + 1] / cell) * KX + Math.floor(g.positions[2 * i] / cell);
    for (const g of on) for (let i = 0; i < g.n; i++) {
      if (!inSel(g, i)) continue;
      const k = binOf(g, i);
      if (!idx.has(k)) idx.set(k, idx.size);
    }
    const n = idx.size;
    const positions = new Float32Array(n * 3), colors = new Uint8Array(n * 4);
    for (const [k, j] of idx) {
      const iy = Math.floor(k / KX), ix = k - iy * KX;
      positions[j * 3] = (ix + 0.5) * cell; positions[j * 3 + 1] = (iy + 0.5) * cell;
    }

    // DENSITY, not presence. Sum the counts falling in each square and scale by how
    // busy the busiest squares are. Normalising each spot against the gene's own
    // per-spot window instead was the bug that made every bin the same color: at
    // 2 um almost every occupied bin holds exactly one UMI, so a 0-1 window sent all
    // of them to the top of the scale. A density has to count over an area, and the
    // scale has to come from the spread of those counts.
    const norm = on.map((g) => {
      const sum = new Float32Array(n);
      const lo = g.min || 0;
      if (prot) {
        // One value per square already (a mean intensity, not a count), so no summing and
        // no percentile: the square's place in the marker's min..max display window.
        const inv = 1 / Math.max(g.max - lo, 1e-6);
        for (let i = 0; i < g.n; i++) {
          if (!inSel(g, i)) continue;
          const v = g.values[i];
          if (v > lo) sum[idx.get(binOf(g, i))] = Math.min(1, (v - lo) * inv);
        }
        return sum;
      }
      for (let i = 0; i < g.n; i++) {
        if (!inSel(g, i)) continue;
        const v = g.values[i];
        if (v >= lo) sum[idx.get(binOf(g, i))] += v;       // the window gates spots, not bins
      }
      // top of scale = a high percentile of the occupied squares, so one hot square
      // cannot flatten everything else, and empty squares are excluded
      const nz = Array.from(sum).filter((v) => v > 0).sort((x, y) => x - y);
      const hi = nz.length ? Math.max(nz[Math.floor((nz.length - 1) * 0.99)], 1e-6) : 1;
      for (let j = 0; j < n; j++) sum[j] = Math.min(1, sum[j] / hi);
      return sum;
    });

    if (ink) {
      const paper = new Float32Array(n * 3).fill(1), miss = new Float32Array(n).fill(1);
      // Marker squares sit over a dark fluorescence image, not white paper: a faint square
      // must read as faint red, not as white haze. So the hue mixes the inks relative to the
      // strongest marker in the square, and how strong it is goes into the opacity alone.
      let dmax = null;
      if (prot) { dmax = new Float32Array(n); for (const d of norm) for (let j = 0; j < n; j++) if (d[j] > dmax[j]) dmax[j] = d[j]; }
      on.forEach((g, gi) => {
        const d = norm[gi], [tr, tg, tb] = hexToRgb255(g.hex);
        for (let j = 0; j < n; j++) {
          if (d[j] <= 0) continue;
          // a square with anything in it is at least a fifth of full ink
          const sEff = 0.2 + 0.8 * d[j];
          const sInk = prot ? d[j] / dmax[j] : sEff;
          paper[j * 3] *= 1 - sInk * (1 - tr / 255);
          paper[j * 3 + 1] *= 1 - sInk * (1 - tg / 255);
          paper[j * 3 + 2] *= 1 - sInk * (1 - tb / 255);
          miss[j] *= 1 - (prot ? d[j] : sEff);
        }
      });
      for (let j = 0; j < n; j++) {
        const c = j * 4;
        colors[c] = Math.round(paper[j * 3] * 255); colors[c + 1] = Math.round(paper[j * 3 + 1] * 255);
        colors[c + 2] = Math.round(paper[j * 3 + 2] * 255);
        const a = Math.max(0, 1 - miss[j]);
        colors[c + 3] = Math.round((prot ? a : Math.sqrt(a)) * 255);
      }
    } else {
      const lut = PALETTES[S.genePalette] || PALETTES.viridis;
      // The AVERAGE of the genes on show, not the sum. Summing meant two genes at
      // half density each already hit the top of the ramp, so adding a second gene
      // flooded the map and every dense area looked identical. An average keeps the
      // ramp meaning the same thing however many genes are on.
      const inv = 1 / Math.max(1, norm.length);
      for (let j = 0; j < n; j++) {
        let total = 0, any = false;
        for (const d of norm) { total += d[j]; if (d[j] > 0) any = true; }
        const t = Math.min(255, Math.round((0.06 + 0.94 * Math.min(1, total * inv)) * 255)) * 3;
        const c = j * 4;
        colors[c] = lut[t]; colors[c + 1] = lut[t + 1]; colors[c + 2] = lut[t + 2];
        colors[c + 3] = any ? 235 : 0;
      }
    }
    S.heat = { key, n, positions, colors, cell, index: Array.from({ length: n }, (_, j) => j) };
    return S.heat;
  }


  // With one gene the picture explains itself. With three it does not: nothing on
  // screen says whether a dark square means "dense" or "two genes overlapping", or
  // what the ramp is averaged over. This says so in a sentence.
  function drawDensityKey() {
    const el = $("densityKey");
    const on = S.genes.filter((g) => g.visible !== false);
    const umap = inUmap();
    if (S.mode === "image" && gatePreviewOn() && !umapByCat()) {
      const G = S.gate, c = G.cut.get(G.cur);
      el.innerHTML = `<div class="genes"><span class="g"><i style="background:transparent;outline:2px solid rgb(${gateLineRgb(G.cur).join(",")});outline-offset:-2px"></i>${escHtml(G.cur)} above ${gateFmt(c.v)}</span></div>` +
        `<div class="note">${umap ? "Lit dots: squares above the line." : "Outlined: squares above the line. Other tissue squares are darkened; the stain is underneath."}` +
        ` ${(gateFracPos(G.vals.get(G.cur), c.v) * 100).toFixed(1)}% of squares.</div>`;
      el.hidden = false;
      return;
    }
    if (isRois() && S.mode === "image" && on.length && !umapByCat() && !umap) { el.innerHTML = window.SodaRoi.keyHtml(roiCtx()); el.hidden = false; return; }
    if (S.mode !== "image" || (S.geneBlend !== "heatmap" && !umap) || !on.length || umapByCat()) { el.hidden = true; return; }
    const e = currentExpr();
    // in the UMAP every dot is one bin, so the key talks about bins rather than squares
    const um = umap ? wd().bin : (e && e.kind === "cells") ? `${e.cellDiameterUm ? `${e.cellDiameterUm.toFixed(0)} µm` : "cell-sized"} square`
      : `${(e && e.binSizeUm) ? `${e.binSizeUm.toLocaleString()} µm` : "each"} square`;
    if (isProt()) {
      const none = umap ? " Grey dots: below the left handle of every marker shown." : "";
      const names = on.map((g) => g.name).join(", ");
      el.innerHTML = S.genePalette === "genes"
        ? `<div class="genes">${on.map((g) => `<span class="g"><i style="background:${g.hex}"></i>${g.name}</span>`).join("")}</div>` +
          `<div class="note">Mean marker intensity per ${um}: ${umap ? "grey" : "clear"} = at or below the left handle, full colour = at the right handle.` +
          (on.length > 1 ? " Where markers overlap the colours mix." : "") + `</div>`
        : `<div class="ramp"><span>low</span><i style="background:${rampCss(S.genePalette)}"></i><span>high</span></div>` +
          `<div class="note">` + (on.length > 1 ? `Mean intensity per ${um}, averaged over ${on.length} markers (${names})` : `Mean ${names} intensity per ${um}`) +
          `.${none}</div>`;
      el.hidden = false;
      return;
    }
    const none = umap ? " Grey dots: none of these genes." : "";
    if (S.genePalette === "genes") {
      el.innerHTML =
        `<div class="genes">${on.map((g) => `<span class="g"><i style="background:${inkCss(g.hex)}"></i>${g.name}</span>`).join("")}</div>` +
        `<div class="note">pale = few transcripts, saturated = many, per ${um}.` +
        (on.length > 1 ? " Where genes overlap the colors mix and darken." : "") + none + `</div>`;
    } else {
      el.innerHTML =
        `<div class="ramp"><span>low</span><i style="background:${rampCss(S.genePalette)}"></i><span>high</span></div>` +
        `<div class="note">` +
        (on.length > 1
          ? `transcripts per ${um}, averaged over ${on.length} genes (${on.map((g) => g.name).join(", ")})`
          : `${on[0].name} transcripts per ${um}`) + `.${none}</div>`;
    }
    el.hidden = false;
  }

  function geneHeatLayer() {
    const h = geneHeat();
    if (!h.n) return null;
    // A 2 um bin is under a tenth of a screen pixel when the whole slide is in
    // view, so the squares vanished entirely at the fit zoom. ScatterplotLayer
    // has radiusMinPixels for exactly this; ColumnLayer has no equivalent, so
    // hold a floor in screen pixels by hand. zoom is log2(pixels per unit).
    const ppu = Math.pow(2, (S.viewState && S.viewState.zoom) || 0);
    const side = Math.max(h.cell * S.pointScale, HEAT_MIN_PX / ppu);
    return new ColumnLayer({
      id: `geneheat-${S.dsId}-${S.exprKey}`,
      // ColumnLayer does not take the binary-attribute form ScatterplotLayer does;
      // fed that way it draws nothing at all, silently. An index array with plain
      // accessors reading out of the same typed arrays costs no extra memory.
      data: h.index,
      getPosition: (i) => [h.positions[i * 3], h.positions[i * 3 + 1], 0],
      getFillColor: (i) => [h.colors[i * 4], h.colors[i * 4 + 1], h.colors[i * 4 + 2], h.colors[i * 4 + 3]],
      diskResolution: 4, angle: 45,            // four sides, corner-up: an axis-aligned square
      radius: side / Math.SQRT2,               // circumradius of a square of side `side`
      radiusUnits: "common", extruded: false, filled: true, stroked: false,
      opacity: S.geneOpacity, parameters: NORMAL,
      updateTriggers: { getFillColor: [h.key], getPosition: [h.key] },
    });
  }

  function geneLayer(g) {
    const e = currentExpr();
    const r = Math.max(0.5, (e && e.spotDiameterPx ? e.spotDiameterPx / 2 : 4) * S.pointScale);
    return new ScatterplotLayer({
      id: `gene-${S.dsId}-${S.exprKey}-${g.name}`,
      data: { length: g.n, attributes: { getPosition: { value: g.positions, size: 2 }, getFillColor: { value: g.colors, size: 4, normalized: true } } },
      dataComparator: () => false,
      radiusUnits: "common", getRadius: r, radiusMinPixels: 0.6,
      stroked: false, filled: true,
      opacity: S.geneOpacity,
      parameters: NORMAL,
      updateTriggers: { getFillColor: [g.colorKey] },
    });
  }

  function tissueLayer() {
    const t = S.tissue;
    if (!t.on || !t.positions) return null;
    const e = currentExpr();
    const r = Math.max(0.5, (e && e.spotDiameterPx ? e.spotDiameterPx / 2 : 4) * S.pointScale * 0.9);
    return new ScatterplotLayer({
      id: `tissue-${S.dsId}-${t.key}`,
      data: { length: t.positions.length / 2, attributes: { getPosition: { value: t.positions, size: 2 } } },
      radiusUnits: "common", getRadius: r, radiusMinPixels: 0.5,
      getFillColor: [120, 130, 150, 90], stroked: false,
      parameters: NORMAL,
    });
  }

  // ---- molecules: every decoded transcript of a picked gene (Xenium ...) --------------------
  // molecules/index.json lists the genes and their counts; molecules/gene/<NAME>.bin holds
  // uint16 qx, qy per molecule, quantised over the image (x = qx * scale[0], y = qy * scale[1]).
  // A file is in random order, so any prefix is a uniform sample of the gene: with the whole
  // section in view only the first k are drawn (k from a budget of dots on screen, in powers of
  // two so zooming redraws only when the sample doubles), zoomed in all of them. One
  // ScatterplotLayer per gene, one colour each, drawn over the cells / bins.
  const MOL_COLORS = cfg.moleculeColors || ["#ff3355", "#33ddff", "#ffdd33", "#66ff66", "#ff9933", "#cc88ff",
    "#ff66cc", "#ffffff", "#3399ff", "#bbff33", "#ff5555", "#00e0b0"];
  const MAX_MOL_GENES = Math.max(1, Math.min(cfg.maxMoleculeGenes || 12, 24));
  // dots in view at once over all genes shown: a desktop GPU draws a few million smoothly
  const MOL_BUDGET = cfg.moleculeBudget || (TOUCH ? 1.5e6 : 6e6);
  // ... and no more than about one dot per 4 screen pixels: zoomed out, a dense gene then
  // reads as density rather than a solid fill
  const MOL_PER_PX = cfg.moleculesPerPixel || 0.25;
  const MOL_RADIUS_UM = 0.8;                       // drawn dot; a transcript is a diffraction-limited spot
  function molBlank() { return { index: null, genes: [], size: 1, opacity: 0.9, sig: "", loading: new Set() }; }
  function molReset() {
    S.mol = molBlank();
    $("panelMolecules").hidden = true; $("molList").innerHTML = ""; $("molNote").textContent = "";
    $("molSize").value = 1; $("molOpacity").value = S.mol.opacity;
  }
  async function molOnDataset() {
    const meta = S.m && S.m.molecules;
    if (!meta || !S.src || !S.src.mol) return;
    try {
      const idx = JSON.parse(new TextDecoder().decode(await S.src.mol("index.json")));
      idx.genes.sort((a, b) => b.c - a.c || a.n.localeCompare(b.n));
      S.mol.index = idx;
      $("panelMolecules").hidden = false;
      $("molMeta").textContent = `${idx.n.toLocaleString()} molecules of ${idx.nGenes.toLocaleString()} genes` +
        (idx.qvMin != null ? `, quality ≥ ${idx.qvMin}` : "") + ". Controls left out.";
      buildMolUI();
    } catch (e) { status(`Could not load the molecule list: ${e.message}`, true); }
  }
  function showMolResults(q) {
    const box = $("molResults"), idx = S.mol.index;
    if (!idx) { box.hidden = true; return; }
    q = q.trim().toLowerCase();
    const hits = [];
    const exact = q ? idx.genes.find((g) => g.n.toLowerCase() === q) : null;
    if (exact) hits.push(exact);
    for (const g of idx.genes) {                   // nothing typed: the most abundant genes
      if (g !== exact && (!q || g.n.toLowerCase().startsWith(q))) hits.push(g);
      if (hits.length >= 40) break;
    }
    if (q && hits.length < 40) for (const g of idx.genes) {
      const n = g.n.toLowerCase();
      if (!n.startsWith(q) && n.includes(q)) hits.push(g);
      if (hits.length >= 40) break;
    }
    box.innerHTML = "";
    hits.forEach((g, i) => {
      const d = document.createElement("div");
      const on = S.mol.genes.some((x) => x.name === g.n);
      d.innerHTML = `<span>${escHtml(g.n)}</span><small>${on ? "shown · " : ""}${g.c.toLocaleString()}</small>`;
      if (i === 0) d.classList.add("hot");
      d.onmousedown = (e) => { e.preventDefault(); addMolGene(g.n); };
      box.appendChild(d);
    });
    box.hidden = hits.length === 0;
  }
  async function addMolGene(name) {
    const M = S.mol, idx = M.index;
    if (!idx || M.loading.has(name) || M.genes.some((g) => g.name === name)) return;
    $("molResults").hidden = true; $("molSearch").value = "";
    const ent = idx.genes.find((g) => g.n === name);
    if (!ent || !ent.c) { status(`${name}: no molecules in this section`, true); setTimeout(() => status(null), 3500); return; }
    if (M.genes.length >= MAX_MOL_GENES) {
      status(`Up to ${MAX_MOL_GENES} genes' molecules at once. Remove one to add ${name}.`, true);
      setTimeout(() => status(null), 4000);
      return;
    }
    status(`Loading ${ent.c.toLocaleString()} ${name} molecules…`);
    M.loading.add(name);
    try {
      const q = new Uint16Array(await S.src.mol(`gene/${encodeURIComponent(name)}.bin`));
      const n = q.length >> 1, [sx, sy] = idx.scale, pos = new Float32Array(2 * n);
      for (let i = 0; i < n; i++) { pos[2 * i] = q[2 * i] * sx; pos[2 * i + 1] = q[2 * i + 1] * sy; }
      const used = new Set(M.genes.map((g) => g.hex));
      const hex = MOL_COLORS.find((c) => !used.has(c)) || MOL_COLORS[M.genes.length % MOL_COLORS.length];
      M.genes.push({ name, n, pos, hex, rgb: hexToRgb255(hex), visible: true, data: null });
      buildMolUI(); render();
    } catch (e) { status(`Could not load ${name} molecules: ${e.message}`, true); return; }
    finally { M.loading.delete(name); }
    status(null);
  }
  function buildMolUI() {
    const M = S.mol, list = $("molList");
    $("molCount").textContent = `(${M.genes.length} / ${MAX_MOL_GENES})`;
    list.innerHTML = "";
    for (const g of M.genes) {
      const row = document.createElement("div");
      row.className = "chan" + (g.visible ? "" : " off");
      row.innerHTML = `
        <input class="cb" type="checkbox" ${g.visible ? "checked" : ""} title="${g.n.toLocaleString()} molecules" />
        <input class="col" type="color" value="${g.hex}" />
        <span class="name" title="${escHtml(g.name)}: ${g.n.toLocaleString()} molecules">${escHtml(g.name)} <small class="muted">${g.n.toLocaleString()}</small></span>
        <button class="rm x" title="remove">×</button>`;
      row.querySelector(".x").onclick = () => { M.genes = M.genes.filter((x) => x !== g); buildMolUI(); render(); };
      row.querySelector(".col").oninput = (e) => { g.hex = e.target.value; g.rgb = hexToRgb255(g.hex); requestRender(); };
      const cb = row.querySelector(".cb");
      cb.onchange = () => { g.visible = cb.checked; row.classList.toggle("off", !cb.checked); render(); };
      row.querySelector(".name").onclick = () => { cb.checked = !cb.checked; cb.onchange(); };
      list.appendChild(row);
    }
  }
  // share of the image in view (the slide pane): how thinly a prefix may be drawn
  function molVisibleFrac() {
    const vs = S.viewState, img = S.m && S.m.image;
    if (!vs || !img || !Number.isFinite(vs.zoom)) return 1;
    const [w, h] = umapSbs() ? umapPaneSize("ortho") : [$("viewport").clientWidth, $("viewport").clientHeight];
    const s = Math.pow(2, vs.zoom), hw = (w || 1000) / s / 2, hh = (h || 1000) / s / 2;
    const x0 = Math.max(0, vs.target[0] - hw), x1 = Math.min(img.width, vs.target[0] + hw);
    const y0 = Math.max(0, vs.target[1] - hh), y1 = Math.min(img.height, vs.target[1] + hh);
    return Math.min(1, Math.max(1e-9, (Math.max(0, x1 - x0) * Math.max(0, y1 - y0)) / (img.width * img.height)));
  }
  function molDrawCounts() {
    const on = S.mol.genes.filter((g) => g.visible);
    const [w, h] = umapSbs() ? umapPaneSize("ortho") : [$("viewport").clientWidth, $("viewport").clientHeight];
    const budget = Math.min(MOL_BUDGET, Math.max(2e5, MOL_PER_PX * (w || 1000) * (h || 1000)));
    const want = budget / Math.max(1, on.length) / molVisibleFrac();
    return on.map((g) => [g, want >= g.n ? g.n : Math.min(g.n, Math.pow(2, Math.floor(Math.log2(Math.max(4096, want)))))]);
  }
  function molLayers() {
    const M = S.mol;
    if (!M || !M.genes.length || inUmap()) { if (M) M.sig = ""; return []; }
    const um = S.m.pixelSizeUm || 0.2125;
    const r = (MOL_RADIUS_UM / um) * M.size;
    const counts = molDrawCounts();
    M.sig = counts.map(([g, k]) => `${g.name}:${k}`).join("|");
    const thin = counts.filter(([g, k]) => k < g.n);
    $("molNote").textContent = thin.length
      ? `Zoomed out: drawing a random ${thin.length === counts.length && counts.length > 1 ? "sample of each gene" : "sample"} (` +
        thin.map(([g, k]) => `${g.name} 1 in ${Math.round(g.n / k).toLocaleString()}`).join(", ") + "). Zoom in to see every molecule."
      : (counts.length ? "Every molecule of the genes shown is drawn." : "");
    return counts.map(([g, k]) => {
      // the same data object while k holds: deck then keeps the uploaded buffer
      if (!g.data || g.data.length !== k) g.data = { length: k, attributes: { getPosition: { value: g.pos.subarray(0, 2 * k), size: 2 } } };
      return new ScatterplotLayer({
        id: `mol-${S.dsId}-${g.name}`, data: g.data,
        getFillColor: [g.rgb[0], g.rgb[1], g.rgb[2], 255], updateTriggers: { getFillColor: [g.hex] },
        radiusUnits: "common", getRadius: r, radiusMinPixels: 1 * M.size, radiusMaxPixels: 8 * M.size,
        stroked: false, filled: true, opacity: M.opacity, parameters: NORMAL,
      });
    });
  }
  // a view change redraws the molecules only when a gene's drawn sample changes size
  function molViewChanged() {
    const M = S.mol;
    if (!M || !M.genes.length || inUmap()) return;
    const sig = molDrawCounts().map(([g, k]) => `${g.name}:${k}`).join("|");
    if (sig !== M.sig) requestRender();
  }

  function overlayLayer(o) {
    if (!o.on || !o.data) return null;
    return new GeoJsonLayer({
      id: `ov-${S.dsId}-${o.key}`, data: o.data,
      stroked: true, filled: true, getFillColor: [255, 255, 255, 18],
      getLineColor: (f) => {
        const c = f.properties && f.properties.classification && f.properties.classification.color;
        return Array.isArray(c) ? [c[0], c[1], c[2], 220] : [255, 220, 80, 220];
      },
      lineWidthUnits: "pixels", getLineWidth: 1.5, pointRadiusMinPixels: 2,
      pickable: true, parameters: NORMAL,
    });
  }

  // what roi.js needs to draw an ROI layer's genes and their key
  function roiCtx() {
    const e = currentExpr(), R = S.region, names = R.meta && R.ids ? new Map(R.meta.map((r) => [r.id, r.name])) : null;
    return { genes: S.genes, pos: (S.lasso && S.lasso.pos) || (S.tissue.key === S.exprKey && S.tissue.positions) || null,
      palette: S.genePalette, palettes: PALETTES, rampCss: rampCss(S.genePalette), hexToRgb255, unit: countUnit(),
      radius: Math.max(0.5, (e && e.spotDiameterPx ? e.spotDiameterPx / 2 : 4) * S.pointScale), opacity: S.geneOpacity,
      rowName: names ? (r) => names.get(R.ids[r]) || "" : null, dsId: S.dsId, exprKey: S.exprKey, ScatterplotLayer, PolygonLayer, NORMAL,
      // the ROI outlines (the layer's region overlay, once loaded) and each row's ROI
      regionIds: R.ids, outlines: ((S.overlays || []).find((o) => o.data && (o.file === R.overlayKey || o.key === R.overlayKey)) || {}).data || null };
  }

  let renderPending = false;  // asked to draw while the pane had no size (hidden tab)
  // Dragging a slider fires `input` continuously, and render() rebuilds layers - on a
  // multi-channel slide that is a tile layer and all its per-tile sublayers, per pixel of
  // travel. Coalesce to at most one render per animation frame: the work still happens,
  // just once per frame the screen can actually show.
  let renderQueued = false;
  function requestRender() {
    if (renderQueued) return;
    renderQueued = true;
    requestAnimationFrame(() => { renderQueued = false; render(); });
  }

  function render() {
    if (!deckgl || !S.m) return;
    const el = $("viewport");
    if (!el.clientWidth || !el.clientHeight) { renderPending = true; return; }   // deck can't build a viewport at 0 px
    renderPending = false;
    if (S.mode === "rna") { renderRna(); return; }
    if (S.umap && S.umap.anim) return;          // the tissue <-> UMAP flight draws its own frames
    if (!S.viewState || !Number.isFinite(S.viewState.zoom) || fitPending) fitView();
    if (inUmap()) { renderUmap(); return; }
    const layers = [];
    const nOn = S.imgOn ? S.channels.filter((ch) => ch.visible).length : 0;
    manyChannelsNote(nOn);
    if (S.imgOn) for (const ch of S.channels) if (ch.visible) layers.push(channelLayer(ch, S.imgOpacity, nOn));
    const tl = tissueLayer(); if (tl) layers.push(tl);
    // Colour by cluster replaces the gene density: one colour per bin, not two on top of each other.
    const gPrev = gatePreviewOn();
    if (!umapByCat() && !gPrev) {
      if (isRois()) {
        layers.push(...window.SodaRoi.geneLayers(roiCtx()));
      } else if (S.geneBlend === "heatmap") {
        const hl = geneHeatLayer();
        if (hl) layers.push(hl);
      } else {
        for (const g of S.genes) if (g.visible !== false) layers.push(geneLayer(g));
      }
    }
    layers.push(...molLayers());
    drawDensityKey();
    updateRegionReadout();
    for (const o of S.overlays) { const l = overlayLayer(o); if (l) layers.push(l); }
    // the gating preview replaces the coloured squares on the slide while it is on
    if (gPrev) layers.push(gatePreviewLayer()); else layers.push(...umapTissueLayers());
    const ql = qcLayer(); if (ql) layers.push(ql);
    layers.push(...lassoLayers());
    const bq = binQcLayer(); if (bq) layers.push(bq);
    layers.push(...umapSideLayers());
    deckgl.setProps({ layers });
  }

  // touch: every stain on is a set of full-size textures. Past a handful a phone can run out
  // of memory and reload the page, so say so once each time the count goes over (not blocked:
  // a newer phone may well cope).
  let manyWarned = false;
  function manyChannelsNote(n) {
    if (!TOUCH) return;
    if (n <= TOUCH_MANY_CHANNELS) { manyWarned = false; return; }
    if (manyWarned) return;
    manyWarned = true;
    const m = `${n} stains on at once: a phone may run out of memory and reload the page. Turning some off keeps it smooth.`;
    status(m);
    setTimeout(() => { if ($("status").textContent === m) status(null); }, 7000);
  }

  // ---- gene colors -------------------------------------------------------------
  function recolorGene(g) {
    const n = g.n, out = new Uint8Array(n * 4);
    const [r, gg, b] = hexToRgb255(g.hex);
    const lo = g.min || 0, inv = 1 / Math.max(g.max - lo, 1e-6);
    const mode = S.geneBlend;
    for (let i = 0; i < n; i++) {
      const raw = g.values[i];
      const k = i * 4;
      if (raw < lo) { out[k + 3] = 0; continue; }                       // below min display: hidden
      const v = Math.min(1, (raw - lo) * inv);
      if (mode === "heatmap") {
        // Viridis, opaque - the same ramp the RNA-seq panels use, so a gene reads the
        // same way in both views. Low is dark blue, high is yellow, and unlike the glow
        // the brightness of a bin is not confounded with how many bins overlap it.
        const t = Math.min(255, Math.round(v * 255)) * 3;
        out[k] = VIRIDIS[t]; out[k + 1] = VIRIDIS[t + 1]; out[k + 2] = VIRIDIS[t + 2]; out[k + 3] = 255;
      } else {
        out[k] = r; out[k + 1] = gg; out[k + 2] = b; out[k + 3] = 40 + 215 * v;
      }
    }
    g.colors = out;
    g.colorKey = `${g.hex}|${lo}|${g.max}|${S.geneBlend}`;
    // Only the low gate decides which spots the density counts, so only a change to
    // it need throw the bins away. Clearing on every window change meant dragging
    // the upper handle recomputed 300k+ squares, and a percentile sort, per pixel of
    // drag - to redraw exactly the same picture.
    // (Protein squares are windowed by both handles, so for them the top one counts too.)
    const hk = g.prot ? `${g.hex}|${lo}|${g.max}` : `${g.hex}|${lo}`;
    if (g.heatKey !== hk) { g.heatKey = hk; S.heat = null; }
  }

  // ---- view -----------------------------------------------------------------------
  let fitPending = false;   // the pane was 0 px when we fitted (hidden tab): refit on resize
  // still showing the fitted view (nobody has panned or zoomed since): a resize, a rotated
  // phone or the controls sheet opening refits it rather than leaving it framed for the old size
  let autoFit = false;
  function fitView() {
    if (!S.m) return;
    autoFit = true;
    if (S.mode === "rna") { fitRna(); return; }
    if (inUmap()) { S.viewState = umapFitVs(); deckgl.setProps({ viewState: deckVs() }); updateHud(); return; }
    if (umapSbs()) {
      // side by side: the tissue fills the left half, the section UMAP the right
      const [w, h] = umapPaneSize("ortho"), [w2, h2] = umapPaneSize("uside");
      S.viewState = imageFitVs(w, h); S.umap.sideVs = umapFitVs(w2, h2);
      fitPending = !w;
      deckgl.setProps({ viewState: deckVs() }); updateHud(); return;
    }
    const img = S.m.image;
    const el = $("viewport");
    const w = el.clientWidth || 0, h = el.clientHeight || 0;
    let zoom = Math.log2(Math.min(w / img.width, h / img.height)) - 0.03;
    if (!Number.isFinite(zoom)) { zoom = Math.log2(1000 / Math.max(img.width, img.height)); fitPending = true; }
    else fitPending = false;
    S.viewState = { target: [img.width / 2, img.height / 2, 0], zoom, minZoom: zoom - 2, maxZoom: 6 };
    deckgl.setProps({ viewState: deckVs() });
    updateHud();
  }
  function updateHud(xy) {
    if (S.mode === "rna") return;
    if (!S.m) return;                               // between two datasets (or one that failed to open): nothing to measure
    const vs = S.viewState; if (!vs) return;
    // UMAP axes have no unit: no coordinates, no scale bar
    $("scalebar").hidden = inUmap();
    if (inUmap()) { $("hudXY").textContent = "UMAP of the whole section"; $("hudZoom").textContent = ""; return; }
    if (!xy && $("hudXY").textContent === "UMAP of the whole section") $("hudXY").textContent = "–";
    $("hudZoom").textContent = `zoom ${vs.zoom.toFixed(2)}  (${(Math.pow(2, vs.zoom)).toFixed(3)} px/px)`;
    if (xy) {
      const um = S.m.pixelSizeUm;
      $("hudXY").textContent = `x ${xy[0].toFixed(0)}  y ${xy[1].toFixed(0)} px` +
        (um ? `   (${(xy[0] * um).toFixed(1)}, ${(xy[1] * um).toFixed(1)} µm)` : "");
    }
    // scale bar
    const sb = $("scalebar");
    const um = S.m.pixelSizeUm;
    const unitsPerPx = 1 / Math.pow(2, vs.zoom);              // image px per screen px
    const worldPerScreen = um ? unitsPerPx * um : unitsPerPx; // µm (or px) per screen px
    const nice = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000];
    let len = nice.find((n) => n / worldPerScreen >= 70) || nice[nice.length - 1];
    sb.querySelector(".bar").style.width = `${len / worldPerScreen}px`;
    sb.querySelector(".lbl").textContent = um ? (len >= 1000 ? `${len / 1000} mm` : `${len} µm`) : `${len} px`;
  }

  function tooltipFor(info) {
    if (info.layer && info.layer.id.includes("ucloud") && info.index >= 0) return umapTooltip(umapRowAt(info.layer.id, info.index));
    if (!info.object || !info.object.properties) return null;
    const p = info.object.properties;
    const name = p.name || (p.classification && p.classification.name) || p.label || p.id;
    // an ROI outline of a GeoMx layer: its name plus the counts of the genes on show (roi.js)
    if (name && p.__region_id && isRois() && !umapByCat()) return { text: window.SodaRoi.tip(roiCtx(), p.__region_id, String(name)) };
    return name ? { text: String(name) } : null;
  }
  // touch: the hover readout (coordinates, the linked bin in the other pane) and the tooltip,
  // for the spot tapped. The tip stays until the next tap or until the view moves.
  function tapInspect(info) {
    const tip = $("tapTip");
    if (S.mode !== "rna" && info.coordinate && !(info.viewport && info.viewport.id === "uside")) updateHud(info.coordinate);
    if (S.mode !== "rna") umapHover(info);
    const t = tooltipFor(info);
    tip.hidden = !t;
    if (!t) return;
    tip.textContent = t.text;
    // centred over the finger, but kept inside the map at the sides and the top
    const W = $("viewport").clientWidth, half = tip.offsetWidth / 2;
    tip.style.left = `${Math.max(half + 6, Math.min(W - half - 6, info.x))}px`; tip.style.top = `${info.y}px`;
    tip.classList.toggle("below", info.y < tip.offsetHeight + 24);
  }

  const ORTHO_CTL = { inertia: 200, scrollZoom: { speed: 0.01, smooth: false }, doubleClickZoom: true, keyboard: false };
  function imageViews() {
    return [new OrthographicView({ id: "ortho", flipY: true, controller: ORTHO_CTL })];
  }
  // Image-mode view states, keyed by view: "ortho" is the slide (or the full-screen UMAP),
  // "uside" the section UMAP in the right half when the two are side by side.
  function deckVs() {
    const out = { ortho: S.viewState };
    if (umapSbs() && S.umap.sideVs) out.uside = S.umap.sideVs;
    return out;
  }
  function initDeck() {
    deckgl = new Deck({
      parent: $("viewport"),
      views: imageViews(),
      // UMAP mode: "umap-base" is drawn in every panel (except under the color-by panel's own
      // colors); "umap-<i>-..." layers only in panel i. Image-mode layers are unaffected.
      layerFilter: ({ layer, viewport }) => {
        const id = layer.id;
        // side by side: "uside-..." layers only in the right half, everything else only in the left
        if (viewport.id === "uside") return id.startsWith("uside-");
        if (id.startsWith("uside-")) return false;
        if (!id.startsWith("umap-")) return true;
        const m = id.match(/^umap-(\d+)-/);
        if (m) return viewport.id === `umap-${m[1]}`;
        return !(viewport.id === "umap-0" && S.rna && S.rna.obsData);
      },
      viewState: { target: [0, 0, 0], zoom: 0 },
      layers: [],
      useDevicePixels: DEVICE_PX,
      // keep the last frame in the buffer so screenshots / "save image" see it
      deviceProps: { webgl: { preserveDrawingBuffer: true } },
      style: { background: "#000" },
      getCursor: ({ isDragging }) => (isDragging ? "grabbing" : "crosshair"),
      onViewStateChange: ({ viewState, viewId }) => {
        autoFit = false;
        $("tapTip").hidden = true;                   // it pointed at a spot that has just moved
        if (S.mode === "rna" && S.rna) { S.rna.vs = viewState; deckgl.setProps({ viewState: rnaViewStates() }); return; }
        if (viewId === "uside") { S.umap.sideVs = viewState; deckgl.setProps({ viewState: deckVs() }); return; }
        const zoomed = !S.viewState || Math.abs(S.viewState.zoom - viewState.zoom) > 1e-6;
        S.viewState = viewState; deckgl.setProps({ viewState: deckVs() }); updateHud();
        // the heat map's squares have a screen-space floor, so their size in world
        // units changes with the zoom
        if (zoomed && ((S.geneBlend === "heatmap" && S.genes.length) || gatePreviewOn())) requestRender();
        molViewChanged();
      },
      onHover: (info) => {
        if (info.coordinate && !(info.viewport && info.viewport.id === "uside")) updateHud(info.coordinate);
        umapHover(info);
      },
      // Tier A: click a region polygon to filter expression to it; click elsewhere clears.
      onClick: (info, ev) => {
        // a finger has no hover: a tap shows what hovering the same spot would
        if (ev && (ev.pointerType === "touch" || (ev.srcEvent && ev.srcEvent.pointerType === "touch"))) tapInspect(info);
        if (S.mode === "rna") return;
        // a dot coloured by cluster: clicking it highlights its cluster, in either view
        if (info && info.layer && info.layer.id.includes("ucloud")) {
          const row = umapRowAt(info.layer.id, info.index), U = S.umap;
          if (row >= 0 && U.cl && U.colorBy === "cluster" && U.cl[row] !== UMAP_NONE) umapPick(U.cl[row]);
          else if (row >= 0 && umapByType() && U.ct[row] !== CT_OUT) umapPickType(U.ct[row]);
          return;
        }
        if (inUmap() || (info && info.viewport && info.viewport.id === "uside")) return;
        if (info && info.layer && info.layer.id === "lasso-shapes") {
          if (info.object && info.object.i >= 0) lassoPickShape(info.object.i);
          return;
        }
        const p = info && info.object && info.object.properties;
        const rid = p && p.__region_id;
        if (rid) selectRegion(rid);
        else if (S.region.selected) selectRegion(null);
      },
      getTooltip: tooltipFor,
      onError: drawProblem,
      onDeviceInitialized: watchContext,
    });
    checkWebGL2();
    window.addEventListener("keydown", (e) => {
      if (e.target && /^(INPUT|SELECT|TEXTAREA)$/.test(e.target.tagName)) return;   // typing a name, not fitting
      if (e.key === "r" || e.key === "R") fitView();
    });
    window.__spatialViewer = { deck: deckgl, state: S, render, fitView,       // for debugging from the console
      // for diag.js (?diag=1): fetch the way the viewer does (its token or password applied), and how tiles are being handled
      diag: () => ({
        kind: S.src && S.src.kind, tokenFlow: !!(TOKEN_URL && token), tokenLeftMin: token ? Math.round((token.expires - Date.now()) / 60000) : null,
        fetch: (url, o) => (S.src && S.src.kind === "api" ? apiFetch(url, o) : dataFetch(url, o)),
        tileUrl, manifestUrl: !S.src ? null : S.src.kind === "api" ? apiUrl(`datasets/${S.dsId}`)
          : S.src.tileUrl(0, 0, 0, 0, "png").replace(/tiles\/0\/0\/0\/0\.png/, "manifest.json"),
        decoder: TILE_DECODERS[tileDecoder][0], decoders: TILE_DECODERS, cached: tileCache.size, cacheMax: tileCacheMax(), touch: TOUCH,
      }) };
  }

  // ---- dual-handle range slider (QuPath "min display / max display") ------------------
  // Two native range inputs stacked on one track; the band between the handles is drawn
  // in the layer's color; an optional histogram sits behind the track. `fmt` turns the
  // 0..max slider units into the label the user should read (raw intensity, counts...).
  // Fluorescence intensities pile up near zero: on a linear 0-255 axis every useful
  // window lands in the first few percent of travel, so the handles are unusable
  // exactly where the detail is. The slider therefore works in position units and
  // maps them through a power curve, stretching the dark end. The histogram is drawn
  // on the same axis, so a peak sits directly under the handle that windows it.
  const SLIDER_GAMMA = cfg.sliderGamma || 2.5;
  const HIST_GAMMA = cfg.histGamma || 0.35;        // count -> bar height; 1 = linear, lower = more shape
  const HIST_CLIP = cfg.histClip || 0.985;         // normalise to this percentile, not the max
  const EXPR_MIN_UM = cfg.exprMinBinUm || 8;        // smallest bin layer offered
  const HEAT_MIN_PX = cfg.heatMinPixels || 1.5;    // smallest a heat-map square may draw
  const HIST_BAR = cfg.histBarWidth || 3;          // px, on the 256-wide histogram canvas
  const HIST_GAP = cfg.histBarGap || 1;
  const POS_MAX = 1000;
  const posToVal = (p, min, span) => min + span * Math.pow(Math.min(1, Math.max(0, p / POS_MAX)), SLIDER_GAMMA);
  const valToPos = (v, min, span) => POS_MAX * Math.pow(Math.min(1, Math.max(0, (v - min) / span)), 1 / SLIDER_GAMMA);

  function rangeSlider(host, o) {
    host.classList.add("rs");
    const span = o.max - o.min || 1;
    const toVal = (p) => Math.round(posToVal(+p, o.min, span));
    const toPos = (v) => Math.round(valToPos(v, o.min, span));
    host.innerHTML = `<canvas class="rs-hist"></canvas><div class="rs-track"></div><div class="rs-band"></div>
      <input class="rs-lo" type="range" min="0" max="${POS_MAX}" step="1" value="${toPos(o.lo)}" title="min display" />
      <input class="rs-hi" type="range" min="0" max="${POS_MAX}" step="1" value="${toPos(o.hi)}" title="max display" />
      <span class="rs-lbl rs-lbl-lo"></span><span class="rs-lbl rs-lbl-hi"></span>`;
    const lo = host.querySelector(".rs-lo"), hi = host.querySelector(".rs-hi"), band = host.querySelector(".rs-band");
    const llo = host.querySelector(".rs-lbl-lo"), lhi = host.querySelector(".rs-lbl-hi");
    const hist = host.querySelector(".rs-hist");

    const paintHist = () => {
      const cols = hist.__cols;
      if (!cols) return;
      const W = hist.width, H = hist.height;
      const g = hist.getContext("2d");
      g.clearRect(0, 0, W, H);
      const a = (+lo.value / POS_MAX) * W, b = (+hi.value / POS_MAX) * W;
      // Drawn as separated bars rather than a filled curve: 1-px columns butted
      // together read as one smooth shape, which is what made a histogram with
      // plenty of structure in it still look like a flat line.
      const step = HIST_BAR + HIST_GAP;
      for (let x0 = 0; x0 < W; x0 += step) {
        let v = 0;
        for (let x = x0; x < Math.min(W, x0 + step); x++) if (cols[x] > v) v = cols[x];
        const bh = Math.round(v * (H - 1));
        if (bh <= 0) continue;
        const mid = x0 + step / 2;
        g.fillStyle = (mid >= a && mid <= b) ? o.color() : "rgba(255,255,255,0.20)";
        g.fillRect(x0, H - bh, HIST_BAR, bh);
      }
    };
    const paint = () => {
      const a = +lo.value / POS_MAX, b = +hi.value / POS_MAX;
      band.style.left = `${a * 100}%`; band.style.width = `${Math.max(0, b - a) * 100}%`;
      band.style.background = o.color();
      llo.textContent = o.fmt(toVal(lo.value)); lhi.textContent = o.fmt(toVal(hi.value));
      paintHist();
    };
    // keep at least one raw unit between the handles, in value space not position space
    lo.oninput = () => { if (toVal(lo.value) >= toVal(hi.value)) lo.value = Math.max(0, +hi.value - 1); paint(); o.onChange(toVal(lo.value), toVal(hi.value)); };
    hi.oninput = () => { if (toVal(hi.value) <= toVal(lo.value)) hi.value = Math.min(POS_MAX, +lo.value + 1); paint(); o.onChange(toVal(lo.value), toVal(hi.value)); };
    // handle on top = the one nearer the pointer, so both are reachable when they overlap
    host.onpointermove = (e) => {
      const r = host.getBoundingClientRect(); const f = (e.clientX - r.left) / r.width;
      const a = +lo.value / POS_MAX, b = +hi.value / POS_MAX;
      const nearLo = Math.abs(f - a) < Math.abs(f - b);
      lo.style.zIndex = nearLo ? 4 : 3; hi.style.zIndex = nearLo ? 3 : 4;
    };
    paint();
    return {
      set: (a, b) => { lo.value = toPos(a); hi.value = toPos(b); paint(); },
      repaint: paint, hist, toPos,
    };
  }
  // 64-bin histogram of one 8-bit channel from its coarsest tile (background pixels = 0 dropped)
  async function drawHistogram(canvas, ch, slider) {
    try {
      // The coarsest level, not level 0. Level 0 is full resolution, so this used to
      // pull one full-res tile per channel on open - 33 of them on the big Akoya
      // slides - which saturated the connection pool and starved the tiles actually
      // being drawn. Worse, tile (0,0) at full res is the top-left corner, which is
      // usually blank glass, so the histogram described background rather than tissue.
      // The coarsest level is one small tile covering the whole slide.
      // z 0 is the coarsest level here: getTileData maps deck's index.z through
      // `maxLevel + index.z`, so `levels - 1` is FULL resolution. Asking for it by
      // that name fetched tile (0,0) of the largest level - the top-left corner of
      // the slide, which is blank glass - so most channels measured nothing and drew
      // a flat line. One small tile covering the whole slide is what was wanted.
      const bmp = await fetchTile(ch.index, 0, 0, 0);
      if (!bmp) return;
      // Sample the tile at its own resolution. Squeezing a 1024-px tile into a
      // 256x128 thumbnail averaged neighbouring pixels together before counting
      // them, which is precisely what destroys a distribution: ~32k blurred samples
      // with the tails averaged off, so every channel drew the same smooth hump.
      const w = Math.min(bmp.width || 256, 1024), h = Math.min(bmp.height || 256, 1024);
      const oc = document.createElement("canvas"); oc.width = w; oc.height = h;
      const cx = oc.getContext("2d", { willReadFrequently: true });
      cx.imageSmoothingEnabled = false; cx.drawImage(bmp, 0, 0, w, h);
      const d = cx.getImageData(0, 0, w, h).data;
      const raw = new Float32Array(256);                      // one bin per 8-bit value
      for (let i = 0; i < d.length; i += 4) { const v = d[i]; if (v > 0) raw[v] += 1; }
      oc.width = oc.height = 0;      // WebKit counts canvas memory until GC and caps the total: hand it back now

      // Spread each value bin across the columns it occupies on the stretched axis and
      // divide by that width, so the bars show density rather than piling every dark
      // value into one spike. Without this the curve would misrepresent the low end
      // exactly where the axis gives it the most room.
      const W = 256, cols = new Float32Array(W);
      const colOf = (v) => (valToPos(v, 0, 255) / POS_MAX) * W;
      for (let v = 1; v < 256; v++) {
        if (!raw[v]) continue;
        const c0 = Math.max(0, Math.floor(colOf(v)));
        const c1 = Math.min(W, Math.max(c0 + 1, Math.ceil(colOf(v + 1))));
        const share = raw[v] / (c1 - c0);
        for (let c = c0; c < c1; c++) cols[c] += share;
      }
      // Counts run over orders of magnitude - a background bin can hold 1e5 pixels where
      // real signal holds 1e2. log1p squashed that to within a factor of two, so every
      // channel drew the same lazy hump. A gentler power compression keeps the shape,
      // and normalising to a high percentile rather than the maximum stops one
      // background spike from flattening everything else against the floor.
      for (let x = 0; x < W; x++) cols[x] = Math.pow(cols[x], HIST_GAMMA);
      const sorted = Array.from(cols).filter((v) => v > 0).sort((a, b) => a - b);
      let top = sorted.length ? sorted[Math.floor((sorted.length - 1) * HIST_CLIP)] : 0;
      if (!(top > 0)) for (const v of cols) if (v > top) top = v;
      if (top > 0) for (let x = 0; x < W; x++) cols[x] = Math.min(1, cols[x] / top);

      canvas.width = W; canvas.height = 34;
      canvas.__cols = cols;
      if (slider) slider.repaint();                            // draws it, tinted by the window
    } catch (e) { /* no histogram is fine */ }
  }
  function rawLabel(ch, v255) {
    // slider units are 0..255 display; show the intensity in the file's own units
    const bits = (S.m.image.channels[ch.index] || {}).bits || [0, 255];
    const raw = bits[0] + (v255 / 255) * (bits[1] - bits[0]);
    return bits[1] > 255 ? String(Math.round(raw)) : String(Math.round(v255));
  }

  // ---- UI: channels ---------------------------------------------------------------
  function syncImageButton() {
    const btn = $("btnImage"); if (!btn) return;
    btn.textContent = S.imgOn ? "Image off" : "Image on";
    btn.title = S.imgOn ? "Hide the image, keeping your channel settings" : "Show the image again";
    $("imgOpacity").disabled = !S.imgOn;
  }

  function filterChannelUI() {
    const q = ($("chanSearch").value || "").trim().toLowerCase();
    let shown = 0;
    for (const row of $("chanList").children) {
      const hit = !q || (row.dataset.name || "").includes(q);
      row.hidden = !hit;
      if (hit) shown += 1;
    }
    $("chanCount").textContent = q ? `(${shown} / ${S.channels.length})` : `(${S.channels.length})`;
  }
  function buildChannelUI() {
    const list = $("chanList"); list.innerHTML = "";
    $("chanCount").textContent = `(${S.channels.length})`;
    for (const ch of S.channels) {
      const row = document.createElement("div");
      row.className = "chan" + (ch.visible ? "" : " off");
      row.innerHTML = `
        <input class="cb" type="checkbox" ${ch.visible ? "checked" : ""} />
        <input class="col" type="color" value="${ch.hex}" />
        <span class="name" title="${ch.name}">${ch.name}</span>
        <button class="x reset" title="reset to the default display range">reset</button>
        <div class="range"></div>`;
      row.dataset.name = ch.name.toLowerCase();
      const cb = row.querySelector(".cb"), col = row.querySelector(".col");
      const slider = rangeSlider(row.querySelector(".range"), {
        min: 0, max: 255, step: 1, lo: ch.lo, hi: ch.hi,
        color: () => ch.hex, fmt: (v) => rawLabel(ch, v),
        onChange: (a, b) => { ch.lo = a; ch.hi = b; requestRender(); },
      });
      drawHistogram(slider.hist, ch, slider);
      cb.onchange = () => { ch.visible = cb.checked; row.classList.toggle("off", !ch.visible); render(); };
      col.oninput = () => { ch.hex = col.value; ch.color = hexToRgb01(col.value); slider.repaint(); requestRender(); };
      row.querySelector(".reset").onclick = () => {
        const w = S.m.image.channels[ch.index].window; ch.lo = Math.round(w[0]); ch.hi = Math.round(w[1]);
        slider.set(ch.lo, ch.hi); render();
      };
      row.querySelector(".name").onclick = () => { cb.checked = !cb.checked; cb.onchange(); };
      list.appendChild(row);
    }
    filterChannelUI();
  }

  // ---- UI: genes ------------------------------------------------------------------
  function currentExpr() {
    return (S.m.expr || []).find((e) => e.key === S.exprKey) || null;
  }
  // ---- protein grid layers (Akoya): the same machinery, other words --------------------
  // An image dataset can carry an expression layer of kind "protein": the slide cut into
  // squares, each holding the mean intensity of every marker. Searching, colouring, the
  // section UMAP, clusters and suggested names all run on it exactly as on Visium bins;
  // only the words change ("squares", "markers", intensity instead of UMI).
  const isProt = () => { const e = S.m && S.mode !== "rna" ? currentExpr() : null; return !!(e && e.kind === "protein"); };
  const WORDS_RNA = { bin: "bin", bins: "bins", Bins: "Bins", gene: "gene", genes: "genes", Gene: "Gene", Genes: "Genes" };
  const WORDS_PROT = { bin: "square", bins: "squares", Bins: "Squares", gene: "marker", genes: "markers", Gene: "Marker", Genes: "Markers" };
  // Cells of a segmented platform (Xenium "cells" layer): the bin machinery, said as cells, and
  // counts are transcripts rather than UMI.
  const isCells = () => { const e = S.m && S.mode !== "rna" ? currentExpr() : null; return !!(e && e.kind === "cells"); };
  const WORDS_CELL = { bin: "cell", bins: "cells", Bins: "Cells", gene: "gene", genes: "genes", Gene: "Gene", Genes: "Genes" };
  // GeoMx ROI segments (geomx.py, a layer of kind "rois"): a few dozen hand-drawn regions, one count
  // vector each. Said as ROIs and counts; drawn as one disc per segment; no UMAP (roi.js).
  const isRois = () => { const e = S.m && S.mode !== "rna" ? currentExpr() : null; return !!(e && e.kind === "rois" && window.SodaRoi); };
  const toRoiWords = (t) => window.SodaRoi.toWords(t);
  const wd = () => (isProt() ? WORDS_PROT : isCells() ? WORDS_CELL : isRois() ? window.SodaRoi.WORDS : WORDS_RNA);
  const countUnit = () => (isCells() ? "transcripts" : isRois() ? "counts" : "UMI");
  const toCellWords = (t) => t.replace(/\bbins\b/g, "cells").replace(/\bbin\b/g, "cell")
    .replace(/\bBins\b/g, "Cells").replace(/\bBin\b/g, "Cell").replace(/\bUMI\b/g, "transcripts");
  // a marker's value in a square: mean intensity on the image's own 0-255 scale
  const protVal = (v) => (v ? (v >= 10 ? v.toFixed(0) : v.toFixed(1)) : "0");
  // Static labels in index.html are the sequencing wording. For a protein layer they are
  // swapped; the originals are kept so a Visium dataset opened afterwards reads exactly as before.
  const toProtWords = (t) => t.replace(/\bbins\b/g, "squares").replace(/\bbin\b/g, "square")
    .replace(/\bgenes\b/g, "markers").replace(/\bgene\b/g, "marker").replace(/\bUMI\b/g, "intensity");
  const WORDING = [
    ["#geneTitle", "textContent", "Markers per square"],
    ["#geneSearch", "placeholder", "Search marker… (Enter adds top hit)"],
    ["#genePalette option[value=genes]", "textContent", "Marker colors (overlap darkens)"],
    ["#lblGeneOpacity", "textContent", "Square opacity"],
    ["#lblPointScale", "textContent", "Square size ×"],
    ["#lblTissueBins", "textContent", "Show all tissue squares"],
    ["#lblColorBy", "textContent", "Colour squares by"],
    ["#colorBy option[value=gene]", "textContent", "Markers (picked in View)"],
    ["#subSel", "textContent", null],
    ["#rowRes span", "title", null], ["#clusterRes", "title", null],
    ["#viewSwitch button[data-v=tissue]", "title", null], ["#viewSwitch button[data-v=sbs]", "title", null],
    ["#viewSwitch button[data-v=umap]", "title", null], ["#flyPlay", "title", null],
    ["#lassoMode option[value=bins]", "textContent", null],
    ["#lassoSaveCsv", "title", null], ["#selToShape", "title", null],
  ];
  const wordingOrig = new Map();
  function syncWording() {
    const prot = isProt(), cells = isCells(), rois = isRois();
    for (const [sel, prop, txt] of WORDING) {
      const el = document.querySelector(sel); if (!el) continue;
      const k = sel + "|" + prop;
      if (!wordingOrig.has(k)) wordingOrig.set(k, el[prop]);
      const orig = wordingOrig.get(k);
      el[prop] = prot ? (txt == null ? toProtWords(orig) : txt) : cells ? toCellWords(orig) : rois ? toRoiWords(orig) : orig;
    }
    if (cells) $("geneTitle").textContent = "Genes per cell";
    if (rois) $("geneTitle").textContent = "Genes per ROI";
    // squares are always drawn as squares: the round-spot "Normal" view has no meaning here
    // (ROI segments are always one disc each: neither view applies)
    $("rowGeneBlend").hidden = prot || rois;
    $("lassoOff").hidden = !prot || !!(currentExpr() || {}).mat;
    if (window.SodaLayout) window.SodaLayout.wording(prot ? toProtWords : cells ? toCellWords : rois ? toRoiWords : null);
  }
  async function loadGeneIndex() {
    S.geneIndex = [];
    if (!S.exprKey) return;
    status(`Loading ${wd().gene} list…`);
    const seq = layerSeq;
    try {
      const list = await S.src.genes(S.exprKey);
      if (seq !== layerSeq) return;                // another layer was picked meanwhile: not this list
      // a panel of ~30 markers reads best alphabetically; a transcriptome by total count
      if (isProt()) list.sort((a, b) => a.n.localeCompare(b.n, undefined, { numeric: true, sensitivity: "base" }));
      else list.sort((a, b) => b.t - a.t);
      S.geneIndex = list;
    } finally { if (seq === layerSeq) status(null); }
  }
  function showGeneResults(q) {
    const box = $("geneResults");
    q = q.trim().toLowerCase();
    const prot = isProt();
    // a marker panel is short: with nothing typed, list all of it (alphabetical)
    if (!q && !prot) { box.hidden = true; return; }
    const hits = [];
    const exact = S.geneIndex.find((g) => g.n.toLowerCase() === q);
    if (exact) hits.push(exact);
    for (const g of S.geneIndex) {
      const n = g.n.toLowerCase();
      if (g !== exact && n.startsWith(q)) hits.push(g);
      if (hits.length >= 40) break;
    }
    if (hits.length < 40) for (const g of S.geneIndex) {
      const n = g.n.toLowerCase();
      if (!n.startsWith(q) && n.includes(q)) hits.push(g);
      if (hits.length >= 40) break;
    }
    box.innerHTML = "";
    hits.forEach((g, i) => {
      const d = document.createElement("div");
      const on = S.genes.some((x) => x.name === g.n);
      d.innerHTML = prot ? `<span>${g.n}</span><small>${on ? "shown" : ""}</small>`
        : `<span>${g.n}</span><small>${Math.round(g.t).toLocaleString()} ${countUnit()}</small>`;
      if (i === 0) d.classList.add("hot");
      d.onmousedown = (e) => { e.preventDefault(); addGene(g.n); };
      box.appendChild(d);
    });
    box.hidden = hits.length === 0;
  }
  async function addGene(name) {
    if (!S.exprKey || S.genes.some((g) => g.name === name)) return;
    $("geneResults").hidden = true; $("geneSearch").value = "";
    if (S.genes.length >= MAX_SPATIAL_GENES) {
      status(`Up to ${MAX_SPATIAL_GENES} ${wd().genes} at once. Remove one to add ${name}.`, true);
      setTimeout(() => status(null), 4000);
      return;
    }
    status(`Loading ${name}…`);
    try {
      const xyv = await S.src.gene(S.exprKey, name);
      const n = xyv.length / 3;
      const positions = new Float32Array(n * 2), values = new Float32Array(n);
      for (let i = 0; i < n; i++) { positions[2 * i] = xyv[3 * i]; positions[2 * i + 1] = xyv[3 * i + 1]; values[i] = xyv[3 * i + 2]; }
      const sorted = Float32Array.from(values).sort();
      const p99 = Math.max(1, percentile(sorted, 0.99));
      const g = { name, hex: GENE_COLORS[S.genes.length % GENE_COLORS.length], n, positions, values,
                  min: 0, max: p99, p99, vmax: sorted[sorted.length - 1] || 1, visible: true };
      if (isProt()) {
        // Every tissue square has a value, and most of it is background. Start the window
        // like a stain's: squares at or below the typical (median) square stay clear, the
        // brightest 1% are full colour. The sliders move it as for an image channel.
        g.prot = true;
        g.min = Math.floor(percentile(sorted, 0.5));
        g.max = Math.max(g.min + 1, Math.ceil(p99));
      }
      recolorGene(g);
      S.genes.push(g);
      buildGeneUI(); render();
    } catch (e) { status(`Could not load ${name}: ${e.message}`, true); return; }
    status(null);
  }
  function buildGeneUI() {
    $("geneCount").textContent = `(${S.genes.length} / ${MAX_SPATIAL_GENES})`;
    const list = $("geneList"); list.innerHTML = "";
    for (const g of S.genes) {
      const row = document.createElement("div");
      row.className = "chan";
      const top = g.prot ? Math.max(2, Math.ceil(g.vmax)) : Math.max(2, Math.ceil(Math.min(g.vmax, g.p99 * 3)));
      row.innerHTML = `
        <input class="cb" type="checkbox" ${g.visible === false ? "" : "checked"} />
        <input class="col" type="color" value="${g.hex}" />
        <span class="name" title="${g.name}">${g.name}</span>
        <button class="rm x" title="remove">×</button>
        <div class="range"></div>`;
      row.querySelector(".x").onclick = () => { S.genes = S.genes.filter((x) => x !== g); buildGeneUI(); render(); };
      const slider = rangeSlider(row.querySelector(".range"), {
        min: 0, max: top, step: top > 50 ? 1 : 0.5, lo: g.min, hi: g.max,
        color: () => g.hex, fmt: (v) => (g.prot ? `${v}` : `${v} ${countUnit()}`),
        onChange: (a, b) => { g.min = a; g.max = b; recolorGene(g); requestRender(); },
      });
      if (g.prot) row.querySelector(".range").title = "Mean intensity per square (0–255, as the image). " +
        "Squares at or below the left handle stay clear; at the right handle and above, full colour.";
      row.querySelector(".col").oninput = (e) => { g.hex = e.target.value; slider.repaint(); recolorGene(g); requestRender(); };
      const cb = row.querySelector(".cb");
      cb.title = g.prot ? `${g.n.toLocaleString()} squares` : `${g.n.toLocaleString()} ${wd().bins} > 0`;
      cb.onchange = () => {
        g.visible = cb.checked;
        row.classList.toggle("off", !cb.checked);
        S.heat = null;                  // the density view sums only what is on
        render();
      };
      row.classList.toggle("off", g.visible === false);
      list.appendChild(row);
    }
  }
  // Resolves true when the layer is fully loaded, false when another layer (or dataset) was
  // picked while it loaded: the newer call does the work, this one stops where it is.
  async function setExprKey(key) {
    const seq = ++layerSeq, stale = () => seq !== layerSeq;
    S.exprKey = key || null;
    S.genes = []; S.tissue.positions = null; S.heat = null; S.tissue.key = key;
    // One square per bin, always: the layer already sets the resolution, and a
    // second multiplier over it only invited the picker and the squares to disagree.
    S.geneBin = 1;
    syncWording();
    if (isProt()) { S.geneBlend = "heatmap"; $("geneBlend").value = "heatmap"; $("rowGenePalette").hidden = false; }
    buildGeneUI();
    await loadGeneIndex(); if (stale()) return false;
    S.region.selected = null;
    await loadRegions(); if (stale()) return false;
    await loadQc(); if (stale()) return false;
    if (binQc) { await binQc.onLayer(S.dsId, S.exprKey, currentExpr(), S.src); if (stale()) return false; }
    if (S.tissue.on) { await loadTissue(); if (stale()) return false; }
    await lassoOnLayer(); if (stale()) return false;
    await umapOnLayer(); if (stale()) return false;
    render();
    await gateOnLayer(); if (stale()) return false;
    syncUmapUI(); render();
    return true;
  }
  async function loadTissue() {
    if (!S.exprKey) return;
    status(`Loading ${wd().bin} positions…`);
    const seq = layerSeq;
    try { const pos = await S.src.positions(S.exprKey); if (seq === layerSeq) S.tissue.positions = pos; }
    finally { if (seq === layerSeq) status(null); }
  }

  // ---- slide QC (marker squares): protgrid.py qc.json + qc_flags.bin -------------------
  // Channels left out or flagged (weak, saturated), squares flagged "everything bright"
  // (folds, bubbles, edges) or with many saturated pixels, the seam check. Older exports
  // have no qc: the box stays hidden.
  const QC_STATUS = { weak: "weak", saturated: "saturated", flat: "flat: not in the UMAP", nuclear: "nuclear: tissue mask only" };
  async function loadQc() {
    const Q = S.qc, e = currentExpr();
    Q.key = null; Q.info = null; Q.flags = null; Q.positions = null;
    if (e && e.qc && e.qc.file && S.src && S.src.buf) {
      const seq = layerSeq;
      try {
        const [info, buf] = await Promise.all([S.src.json(`expr/${S.exprKey}/${e.qc.file}`),
          S.src.buf(`expr/${S.exprKey}/${e.qc.flags || "qc_flags.bin"}`)]);
        if (seq !== layerSeq) return;
        Q.info = info; Q.flags = new Uint8Array(buf); Q.key = S.exprKey;
      } catch (err) { if (seq !== layerSeq) return; Q.info = null; Q.flags = null; }
    }
    buildQcBox();
  }
  function buildQcBox() {
    const Q = S.qc, box = $("qcBox");
    box.hidden = !Q.info;
    if (!Q.info) { Q.outline = false; $("qcOutline").checked = false; return; }
    const I = Q.info, esc = escHtml;
    const ch = I.channels || [], weak = ch.filter((c) => c.weak && c.inUmap), sat = ch.filter((c) => c.saturated);
    const nUmap = ch.filter((c) => c.inUmap).length, drop = I.dropped || [];
    const pct = (n) => (I.squares ? `${(100 * n / I.squares).toFixed(2)}%` : "");
    const br = I.everythingBright || {}, ss = I.saturatedSquares || {}, sm = I.seams;
    const lines = [
      `<b>${nUmap}</b> markers in the UMAP; ${drop.length} left out` +
        (drop.length ? ` (<span title="${esc(drop.map((d) => `${d.name}: ${d.reason}`).join("\n"))}">${esc(drop.map((d) => d.name).join(", "))}</span>)` : "") + ".",
      weak.length ? `<span class="qc-warn" title="${esc((I.rules && I.rules.weak) || "")}">Weak stain</span>: ${esc(weak.map((c) => c.name).join(", "))}. Kept in the UMAP; check them by eye.` : "No weak stains.",
      sat.length ? `<span class="qc-warn" title="${esc((I.rules && I.rules.saturated) || "")}">Saturated</span>: ${esc(sat.map((c) => `${c.name} ${c.saturatedPct.toFixed(1)}%`).join(", "))} of tissue pixels at the detector's limit.` : "No saturated stains.",
      `<span title="${esc(br.rule || "")}">Everything bright</span> (folds, bubbles, edges): ${(br.n || 0).toLocaleString()} squares (${pct(br.n || 0)}).` +
        ` <span title="${esc(ss.rule || "")}">Many saturated pixels</span>: ${(ss.n || 0).toLocaleString()} (${pct(ss.n || 0)}).`,
      sm ? (sm.flagged ? `<span class="qc-warn" title="${esc(sm.rule || "")}">Possible tile seams</span>: ` +
                         ["columns", "rows"].filter((k) => sm[k] && sm[k].flagged).map((k) =>
                           `<span title="${esc(sm[k].channels.join(", "))}">${sm[k].channels.length} of ${sm[k].channelsChecked} markers ` +
                           `(${esc(sm[k].channels.slice(0, 4).join(", "))}${sm[k].channels.length > 4 ? " …" : ""})</span> step every ` +
                           `${sm[k].spacingPx} px (~${sm[k].spacingUm} µm) ${k === "columns" ? "across" : "down"} the slide`).join("; ") +
                         (sm.correction ? "." : `. Clusters can follow such steps: if a cluster pattern repeats at this spacing, compare it with the image before trusting it.` +
                                                (sm.gridProblem ? ` Not corrected: ${esc(sm.gridProblem)}.` :
                                                 sm.correctionSkipped ? ` Not corrected: ${esc(sm.correctionSkipped)}.` : ""))
                       : `<span title="${esc(sm.rule || "")}">No tile seams found.</span>`) : "",
      sm && sm.correction ? (() => {
        const co = sm.correction, nc = co.notCorrectable || [];
        return `<span class="qc-ok" title="${esc(co.rule || "")}">Seams corrected</span> in ${co.corrected.length} of ` +
          `${co.corrected.length + nc.length} markers with steps (step at tile edges ${co.excessBeforeMedian} → ` +
          `${co.excessAfterMedian} points, median). The numbers are corrected; the image is shown as recorded.` +
          (nc.length ? ` <span class="qc-warn" title="${esc(nc.map((x) => `${x.name}: ${x.reason}`).join("\n"))}">Not correctable</span>: ` +
                       `${esc(nc.map((x) => x.name).join(", "))} (kept in the UMAP; their clusters may still follow the tiles).` +
                       (nc.filter((x) => /shading shared/.test(x.reason || "")).length
                         ? ` For ${nc.filter((x) => /shading shared/.test(x.reason || "")).length} of them the seam is shading ` +
                           `shared by every field (one side of every field brighter), which needs a flat-field correction (not built yet).` : "") : "");
      })() : "",
      `<i>${esc(I.note || "")}</i>`,
    ].filter(Boolean);
    $("qcSummary").innerHTML = lines.map((l) => `<div>${l}</div>`).join("");
    $("rowQcOutline").hidden = !((br.n || 0) + (ss.n || 0));
    $("qcChannels").innerHTML = `<table class="qc-tab"><tr><th>Marker</th><th title="Off-tissue level (0–255)">Bg</th>` +
      `<th title="(99.5th percentile − background) / background">Signal</th><th title="Share of tissue pixels at the detector's limit">Sat.</th>` +
      `<th title="Tile seams: corrected (and how), or not correctable">Seams</th><th></th></tr>` +
      ch.map((c) => `<tr class="${c.status === "ok" ? "" : "qc-" + c.status}"><td>${esc(c.name)}</td><td>${c.background.toFixed(1)}</td>` +
        `<td>${c.signalToBackground.toFixed(0)}</td><td>${c.saturatedPct.toFixed(2)}%</td>` +
        `<td class="${c.seams === "not correctable" ? "qc-warnc" : ""}">${esc(c.seams || "")}</td><td>${esc(QC_STATUS[c.status] || "")}</td></tr>`).join("") +
      `</table>`;
  }
  async function qcSetOutline(on) {
    const Q = S.qc;
    Q.outline = !!on && !!Q.flags;
    if (Q.outline && !Q.positions) {
      const pos = S.tissue.positions && S.tissue.key === S.exprKey ? S.tissue.positions : await S.src.positions(S.exprKey);
      const idx = [];
      for (let i = 0; i < Q.flags.length; i++) if (Q.flags[i]) idx.push(i);
      const xy = new Float32Array(idx.length * 2), col = new Uint8Array(idx.length * 4);
      idx.forEach((i, k) => {
        xy[2 * k] = pos[2 * i]; xy[2 * k + 1] = pos[2 * i + 1];
        col.set(Q.flags[i] & 1 ? [255, 170, 40, 255] : [230, 60, 230, 255], 4 * k);   // orange: everything bright, magenta: saturated
      });
      Q.positions = { xy, col, n: idx.length };
    }
    render();
  }
  function qcLayer() {
    const Q = S.qc;
    if (!Q.outline || !Q.positions || Q.key !== S.exprKey) return null;
    const e = currentExpr(), r = Math.max(0.5, e && e.spotDiameterPx ? e.spotDiameterPx / 2 : 4);
    return new ScatterplotLayer({
      id: `qc-${S.dsId}-${Q.key}`,
      data: { length: Q.positions.n, attributes: { getPosition: { value: Q.positions.xy, size: 2 },
                                                   getLineColor: { value: Q.positions.col, size: 4, normalized: true } } },
      radiusUnits: "common", getRadius: r, filled: false, stroked: true,
      lineWidthUnits: "pixels", getLineWidth: 1.5, parameters: NORMAL,
    });
  }

  // ---- Tier A: region linking --------------------------------------------------------
  // region.bin is aligned to positions.bin, so build one quantized-xy -> region-id map
  // on load. Gene points carry the same bin-centre coordinates (a subset), so a point
  // looks up its region by the same key. Quantizing to 1/100 px: bins are spaced several
  // px apart, and the float32 values come from the same source as positions.bin, so the
  // rounding matches exactly.
  const binKey = (x, y) => `${Math.round(x * 100)},${Math.round(y * 100)}`;
  const regionAt = (x, y) => (S.region.lookup ? (S.region.lookup.get(binKey(x, y)) || 0) : 0);

  async function loadRegions() {
    const e = currentExpr();
    S.region.meta = (e && e.regions) || null;
    S.region.ids = null; S.region.lookup = null;
    S.region.overlayKey = (e && e.regionOverlay) || null;
    if (!S.region.meta) return;
    const seq = layerSeq;
    try {
      const [pos, ids] = await Promise.all([S.src.positions(S.exprKey), S.src.region(S.exprKey)]);
      if (seq !== layerSeq) return;
      const lookup = new Map();
      const n = Math.min(ids.length, pos.length / 2);
      for (let i = 0; i < n; i++) if (ids[i]) lookup.set(binKey(pos[2 * i], pos[2 * i + 1]), ids[i]);
      S.region.ids = ids; S.region.lookup = lookup;
    } catch (e2) { if (seq === layerSeq) S.region.meta = null; }     // region.bin absent -> feature simply stays off
  }

  function selectRegion(rid) {
    S.region.selected = rid || null;
    S.heat = null;                              // density is now over the selected region only
    render();
  }

  function updateRegionReadout() {
    let el = $("regionReadout");
    if (!el) {
      el = document.createElement("div");
      el.id = "regionReadout";
      el.style.cssText = "position:fixed;left:12px;bottom:44px;z-index:30;max-width:320px;" +
        "background:rgba(10,12,16,.85);color:#e8eef5;font:12px/1.45 system-ui,sans-serif;" +
        "padding:8px 10px;border:1px solid #2a3340;border-radius:6px;";
      document.body.appendChild(el);
    }
    const sel = S.region.selected;
    if (S.mode === "rna" || inUmap() || !sel || !S.region.meta) { el.hidden = true; el.dataset.sig = ""; return; }
    const on = S.genes.filter((g) => g.visible !== false);
    // recompute only when the selection or gene set changes, not on every pan/zoom frame
    const sig = `${sel}|${S.exprKey}|${on.map((g) => g.name + ":" + g.hex).join(",")}`;
    if (!el.hidden && el.dataset.sig === sig) return;
    el.dataset.sig = sig;
    const meta = S.region.meta.find((r) => r.id === sel);
    const rows = on.map((g) => {
      let sum = 0, k = 0;
      for (let i = 0; i < g.n; i++)
        if (regionAt(g.positions[2 * i], g.positions[2 * i + 1]) === sel) { sum += g.values[i]; k++; }
      return `<div style="display:flex;justify-content:space-between;gap:12px">` +
        `<span><i style="display:inline-block;width:8px;height:8px;background:${g.hex};margin-right:5px"></i>${g.name}</span>` +
        `<span>${k.toLocaleString()} bins · ${(k ? sum / k : 0).toFixed(2)} mean</span></div>`;
    }).join("");
    el.innerHTML =
      `<div style="display:flex;justify-content:space-between;gap:12px;margin-bottom:4px">` +
      `<b>${meta ? meta.name : "region " + sel}</b>` +
      `<a href="#" id="regionClear" style="color:#8fb6ff;text-decoration:none">clear ✕</a></div>` +
      (rows || `<div style="opacity:.7">add a gene to see its expression in this region</div>`);
    el.hidden = false;
    const clr = $("regionClear");
    if (clr) clr.onclick = (ev) => { ev.preventDefault(); selectRegion(null); };
  }

  // ---- Circle -> UMAP ----------------------------------------------------------------
  // Draw circles or freehand shapes over a Visium section. The bins inside go through a
  // UMAP computed in a worker (lasso_worker.js) from the export's hvg/ chunks - only the
  // chunks under the shapes are fetched - and land in a side panel coloured by shape.
  // Marker squares (protein grid) work the same way from their mat/ chunks (the embedded
  // channels' mean intensities; protgrid.write_mat_chunks): L.prot, and the worker's runProtein.
  // Brushing the panel lights the bins on the slide; picking a shape lights its points.
  // A bin inside two shapes belongs to the one drawn last, never to both.
  const LASSO_COLORS = ["#e6194b", "#3cb44b", "#4363d8", "#f58231", "#911eb4", "#42d4f4",
                        "#f032e6", "#bfef45", "#469990", "#9a6324", "#fabed4", "#dcbeff"];
  // nPcs caps the worker's PCs-by-bin-count rule (15 / 20 / 30 PCs, 10 / 15 neighbours below
  // 1,000 / up to 5,000 / above); of the section's variable genes, those counted in fewer than
  // max(minDetect, min(minDetectFrac x bins, minDetectCap)) bins are left out.
  const LASSO = { minPoints: 50, fewPoints: 300, smallShape: 20, softCap: 20000, cap: 50000, minUmi: 10,
                  nPcs: 30, nPcsProtein: 15, nNeighbors: 15, minDist: 0.3, minShapesForUmap: 15, umapJs: "1.4.0",
                  minDetect: 5, minDetectFrac: 0.01, minDetectCap: 30 };
  // The run is single-threaded: 50,000 bins took 76 s on a desktop core and would take 2-5 min on
  // a phone, long enough to be frozen by a locked screen. On touch the sample stops at the soft cap.
  const LASSO_CAP = TOUCH ? LASSO.softCap : LASSO.cap;
  const lassoBlank = () => ({ tool: null, shapes: [], next: 1, drawing: null, pos: null, hvg: null,
    chunks: new Map(), assign: null, bin: null, result: null, brush: null, selShape: null, hover: -1,
    worker: null, running: false, lpVs: null,
    // Leiden clusters of the last result (absent when the worker skipped clustering)
    selCluster: null, lpColor: "cluster", clNames: {} });
  let lpDeck = null;
  // per-bin QC filter (qc.js): removed bins are greyed on the slide and left out of the run
  let binQc = null;
  function binQcInit() {
    if (!window.SodaBinQc) return;
    binQc = window.SodaBinQc.create({
      $, esc: escHtml, ScatterplotLayer, NORMAL, requestRender, layerKey: () => S.exprKey,
      positions: () => (S.lasso && S.lasso.pos ? Promise.resolve(S.lasso.pos)
        : S.tissue.positions && S.tissue.key === S.exprKey ? Promise.resolve(S.tissue.positions) : S.src.positions(S.exprKey)),
      onChange: () => {
        const L = S.lasso; if (!L) return;
        buildLassoUI();
        $("bqStale").hidden = !(L.result && JSON.stringify(L.result.qc || null) !== JSON.stringify(binQc.state()));
      },
    });
    if (window.__spatialViewer) window.__spatialViewer.binQc = binQc;   // for testing from the console
  }
  function binQcLayer() {
    if (!binQc) return null;
    const e = currentExpr();
    return binQc.layer(Math.max(0.5, (e && e.spotDiameterPx ? e.spotDiameterPx / 2 : 4) * S.pointScale));
  }
  // bins inside the shapes that the QC filter removes
  function binQcInShapes() {
    const L = S.lasso, k = binQc && binQc.keepMask();
    if (!k || !L.assign || k.length !== L.assign.length) return 0;
    let n = 0; for (let i = 0; i < k.length; i++) if (!k[i] && L.assign[i] !== LASSO_NONE) n++;
    return n;
  }
  const LASSO_NONE = 0xffff;

  function lassoReset() {
    lassoCancel();                                 // a run of the dataset being left: its Cancel button, the screen lock
    if (S.lasso && S.lasso.worker) S.lasso.worker.terminate();
    S.lasso = lassoBlank();
    $("panelLasso").hidden = true; $("lassoPanel").hidden = true; $("drawPad").hidden = true;
    if (window.SodaVolcano) window.SodaVolcano.reset();   // a volcano of the other dataset's bins
    buildLassoUI();
  }

  // Called whenever the expression layer changes: the shapes stay (they are in image
  // pixels), the bins under them are found again in the new layer.
  async function lassoOnLayer() {
    const L = S.lasso, e = currentExpr();
    lassoCancel(); L.result = null; L.brush = null; L.selShape = null; L.selCluster = null; L.chunks = new Map(); L.hvgGenes = null;
    $("lassoPanel").hidden = true;
    if (window.SodaVolcano) window.SodaVolcano.reset();
    if (window.SodaPlots) window.SodaPlots.reset();
    const prot = !!(e && e.kind === "protein" && e.mat);
    const ok = S.src && S.src.kind === "static" && e && (e.hvg || prot) && S.mode !== "rna";
    $("panelLasso").hidden = !ok;
    if (!ok) { lassoSetTool(null); L.pos = null; L.hvg = null; L.prot = false; return; }
    const seq = layerSeq;
    try {
      const [hvg, pos] = await Promise.all([S.src.json(`expr/${S.exprKey}/${prot ? "mat" : "hvg"}/index.json`),
        S.tissue.positions && S.tissue.key === S.exprKey ? S.tissue.positions : S.src.positions(S.exprKey)]);
      if (seq !== layerSeq || S.lasso !== L) return;        // the layer changed meanwhile
      // a protein index names its channels: they stand in for the gene count and names
      if (prot) { hvg.nGenes = hvg.nChannels; L.hvgGenes = hvg.channels; }
      L.hvg = hvg; L.pos = pos; L.prot = prot;
    } catch (err) { if (seq === layerSeq && S.lasso === L) $("panelLasso").hidden = true; return; }
    $("lassoNote").textContent = prot
      ? `Uses the ${e.label} layer: the ${L.hvg.nChannels} markers of the section UMAP, each scaled as the whole ` +
        (L.hvg.normalize === "bg-asinh-v1" ? `slide was (background taken off, divided by its 99.5th percentile, asinh), `
                                           : `section was (1st to 99.5th percentile), `) +
        `then compared within the shapes for the UMAP. Each shape is also compared with the rest of the slide.`
      : `Uses the ${e.label} layer. The UMAP can only use the ${L.hvg.nGenes.toLocaleString()} ` +
        `most variable genes of the whole section; inside the shapes it keeps those counted in at least ` +
        `${LASSO.minDetect}–${LASSO.minDetectCap} of the ${wd().bins} (1% of them, within those limits).`;
    if (window.SodaRoi) window.SodaRoi.syncLasso(isRois(), e, L.hvg);     // ROI layers: shapes yes, the panel's own UMAP no
    lassoAssign();
  }

  function lassoSetTool(t) {
    const L = S.lasso;
    if (t && inUmap()) t = null;                 // shapes are drawn on the tissue, not on the UMAP
    if (t) umapSetTool(null);                    // one drawing tool at a time
    L.tool = t && L.tool !== t ? t : null;
    $("lassoCircle").classList.toggle("on", L.tool === "circle");
    $("lassoFree").classList.toggle("on", L.tool === "free");
    $("drawPad").hidden = !L.tool;
  }

  const lassoRgb = (s) => hexToRgb255(s.hex);
  function circlePts(c, r) {
    const pts = [];
    for (let k = 0; k < 72; k++) { const a = (k / 72) * 2 * Math.PI; pts.push([c[0] + r * Math.cos(a), c[1] + r * Math.sin(a)]); }
    return pts;
  }
  function inPoly(x, y, pts) {
    let inside = false;
    for (let i = 0, j = pts.length - 1; i < pts.length; j = i++) {
      const [xi, yi] = pts[i], [xj, yj] = pts[j];
      if ((yi > y) !== (yj > y) && x < xi + ((xj - xi) * (y - yi)) / (yj - yi)) inside = !inside;
    }
    return inside;
  }

  // shape index per positions.bin row (LASSO_NONE = none); later shapes win overlaps.
  // Unsigned on purpose: the deid gate reads the signed typed-array names as sample IDs.
  function lassoAssign() {
    const L = S.lasso;
    if (!L.pos) return;
    const pos = L.pos, n = pos.length / 2, a = new Uint16Array(n).fill(LASSO_NONE);
    L.shapes.forEach((s) => { s.bins = 0; s.shared = 0; s.lost = 0; });
    L.shapes.forEach((s, si) => {
      const take = (i) => {
        if (a[i] !== LASSO_NONE) { const o = L.shapes[a[i]]; o.bins--; o.lost++; s.shared++; }
        a[i] = si; s.bins++;
      };
      // a shape made from a selection (Tissue / UMAP "Use selection as a shape") is a list of
      // bins, not an outline: it can be any set of bins, scattered over the section
      if (s.rows) { for (const i of s.rows) if (i < n) take(i); return; }
      let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity;
      for (const [x, y] of s.pts) { x0 = Math.min(x0, x); y0 = Math.min(y0, y); x1 = Math.max(x1, x); y1 = Math.max(y1, y); }
      for (let i = 0; i < n; i++) {
        const x = pos[2 * i], y = pos[2 * i + 1];
        if (x < x0 || x > x1 || y < y0 || y > y1 || !inPoly(x, y, s.pts)) continue;
        take(i);
      }
    });
    L.assign = a;
    lassoRefresh();
    buildLassoUI();
  }

  // Whatever changes the shapes makes any result stale: drop it rather than show a
  // UMAP that no longer matches what is drawn.
  function lassoShapesChanged() {
    const L = S.lasso;
    lassoCancel(); L.result = null; L.brush = null; L.hover = -1; L.selCluster = null;
    if (L.selShape !== null && L.selShape >= L.shapes.length) L.selShape = null;
    $("lassoPanel").hidden = true;
    lassoAssign();
  }

  // ---- two fingers on a drawing pad move the picture under it -------------------------------
  // While a drawing tool is on, its pad sits over the map, so the map never sees those fingers.
  // One finger draws; a second ends that stroke (dropStroke) and the pair pans and pinch-zooms
  // the pane they are on instead. Listens in the capture phase so that, while two fingers are
  // down, the pad's own drawing handlers never see them.
  function padPinch(pad, dropStroke) {
    const pts = new Map();
    let P = null;                                    // the pinch: pane, its start view, the world point under the fingers
    const local = (ev) => { const r = $("viewport").getBoundingClientRect(); return [ev.clientX - r.left, ev.clientY - r.top]; };
    const centre = () => {
      const [a, b] = [...pts.values()];
      return [(a[0] + b[0]) / 2, (a[1] + b[1]) / 2, Math.hypot(a[0] - b[0], a[1] - b[1])];
    };
    pad.addEventListener("pointerdown", (ev) => {
      if (ev.pointerType !== "touch") return;
      pts.set(ev.pointerId, local(ev));
      if (pts.size < 2) return;
      ev.stopImmediatePropagation();
      if (P || pts.size > 2) return;
      dropStroke();
      const [mx, my, d] = centre(), vps = deckgl.getViewports();
      const vp = vps.find((v) => mx >= v.x && mx < v.x + v.width && my >= v.y && my < v.y + v.height) || vps[0];
      P = { vp, vs0: vp.id === "uside" ? S.umap.sideVs : S.viewState, w: vp.unproject([mx - vp.x, my - vp.y]), d0: Math.max(d, 1) };
    }, true);
    pad.addEventListener("pointermove", (ev) => {
      if (!pts.has(ev.pointerId)) return;
      pts.set(ev.pointerId, local(ev));
      if (!P) return;
      ev.stopImmediatePropagation();
      if (pts.size < 2) return;
      const [mx, my, d] = centre(), { vp, vs0, w, d0 } = P;
      const lo = Number.isFinite(vs0.minZoom) ? vs0.minZoom : -Infinity, hi = Number.isFinite(vs0.maxZoom) ? vs0.maxZoom : Infinity;
      const zoom = Math.max(lo, Math.min(hi, vs0.zoom + Math.log2(d / d0))), s = Math.pow(2, zoom);
      // keep the world point that started under the fingers under them (every pane here is flipY)
      const target = [w[0] - (mx - vp.x - vp.width / 2) / s, w[1] - (my - vp.y - vp.height / 2) / s, 0];
      deckgl.props.onViewStateChange({ viewState: { ...vs0, target, zoom }, viewId: vp.id });
    }, true);
    const up = (ev) => {
      if (!pts.delete(ev.pointerId)) return;
      if (P) ev.stopImmediatePropagation();
      if (!pts.size) P = null;                         // the last finger lifted: the next touch draws again
    };
    pad.addEventListener("pointerup", up, true);
    pad.addEventListener("pointercancel", up, true);
  }

  function lassoPadEvents() {
    const pad = $("drawPad");
    padPinch(pad, () => { S.lasso.drawing = null; requestRender(); });
    const world = (ev) => {
      const r = $("viewport").getBoundingClientRect();
      return deckgl.getViewports()[0].unproject([ev.clientX - r.left, ev.clientY - r.top]).slice(0, 2);
    };
    pad.addEventListener("pointerdown", (ev) => {
      const L = S.lasso;
      if (!L.tool || ev.button !== 0) return;
      pad.setPointerCapture(ev.pointerId);
      const p = world(ev);
      L.drawing = { kind: L.tool, c: p, pts: [p], sx: ev.clientX, sy: ev.clientY, lx: ev.clientX, ly: ev.clientY };
    });
    pad.addEventListener("pointermove", (ev) => {
      const D = S.lasso.drawing; if (!D) return;
      const p = world(ev);
      if (D.kind === "circle") D.pts = circlePts(D.c, Math.hypot(p[0] - D.c[0], p[1] - D.c[1]));
      else if (Math.hypot(ev.clientX - D.lx, ev.clientY - D.ly) >= 3) { D.pts.push(p); D.lx = ev.clientX; D.ly = ev.clientY; }
      requestRender();
    });
    const finish = (ev) => {
      const L = S.lasso, D = L.drawing; if (!D) return;
      L.drawing = null;
      const big = Math.hypot(ev.clientX - D.sx, ev.clientY - D.sy) >= 4;
      if ((D.kind === "circle" && big) || (D.kind === "free" && D.pts.length >= 3)) {
        const id = L.next++;
        L.shapes.push({ id, name: `Shape ${id}`, kind: D.kind, pts: D.pts, hex: LASSO_COLORS[(id - 1) % LASSO_COLORS.length],
                        bins: 0, shared: 0, lost: 0 });
        lassoShapesChanged();
      }
      requestRender();
    };
    pad.addEventListener("pointerup", finish);
    pad.addEventListener("pointercancel", finish);
  }

  // Slide-side colours: each bin in its shape's colour; when a brush or a picked shape
  // is active, everything outside it dims.
  function lassoRefresh() {
    const L = S.lasso;
    if (!L.assign || !L.pos) { L.bin = null; lassoSyncSel(); requestRender(); renderLassoPanel(); return; }
    const a = L.assign, pos = L.pos, R = L.result;
    let n = 0; for (let i = 0; i < a.length; i++) if (a[i] !== LASSO_NONE) n++;
    const positions = new Float32Array(2 * n), colors = new Uint8Array(4 * n), rgb = L.shapes.map(lassoRgb);
    const lit = lassoLitRows();
    // coloured by cluster in the panel: the same colours on the slide, so the two read as one
    const byCl = lassoByCluster() && R.info.mode === "bins";
    let k = 0;
    for (let i = 0; i < a.length; i++) {
      if (a[i] === LASSO_NONE) continue;
      positions[2 * k] = pos[2 * i]; positions[2 * k + 1] = pos[2 * i + 1];
      let c = rgb[a[i]];
      if (byCl) { const j = R.idxOfRow.get(i); c = j === undefined ? [120, 120, 120] : clusterRgb(R.cluster[j]); }
      colors[4 * k] = c[0]; colors[4 * k + 1] = c[1]; colors[4 * k + 2] = c[2];
      // see-through by default so the H&E under a shape stays readable
      colors[4 * k + 3] = lit ? (lit(i, a[i]) ? 235 : 12) : (R ? 110 : 80);
      k++;
    }
    L.bin = { n, data: { length: n, attributes: { getPosition: { value: positions, size: 2 }, getFillColor: { value: colors, size: 4, normalized: true } } },
              key: `${n}|${Date.now()}` };
    umapDirty();                                   // a picked shape / brush also lights the main UMAP
    if (S.umap) syncUmapUI();                      // the chip names the pick
    lassoSyncSel();
    requestRender();
    renderLassoPanel();
  }
  // the panel's own line under the plot: what is picked, Show on slide, Clear
  function lassoSyncSel() {
    const t = lassoPickTxt();
    $("lpSel").hidden = !t;
    $("lpSelTxt").textContent = t.replace(/^Drawn region: /, "");
    $("lpHint").hidden = !!t || !(S.lasso && S.lasso.result) || !!S.lasso.hinted;
  }

  // (row, shape) -> lit? for the current brush / picked shape; null when nothing is picked
  function lassoLitRows() {
    const L = S.lasso, R = L.result;
    if (L.brush && R) {
      if (R.info.mode === "shapes") {
        const on = new Set(); L.brush.forEach((b, j) => { if (b) on.add(R.shape[j]); });
        return (row, s) => on.has(s);
      }
      return (row) => { const j = R.idxOfRow.get(row); return j !== undefined && L.brush[j] === 1; };
    }
    if (L.selShape !== null) return (row, s) => s === L.selShape;
    if (L.selCluster !== null && R && R.cluster) {
      return (row) => { const j = R.idxOfRow.get(row); return j !== undefined && R.cluster[j] === L.selCluster; };
    }
    return null;
  }
  function lassoPickTxt() {
    const L = S.lasso, R = L && L.result;
    if (!L) return "";
    if (L.brush) { let k = 0; for (const b of L.brush) k += b; return `Drawn region: ${k.toLocaleString()} point${k === 1 ? "" : "s"} selected`; }
    if (L.selShape !== null && L.shapes[L.selShape]) return `Drawn region: ${L.shapes[L.selShape].name}`;
    if (L.selCluster !== null && R) return `Drawn region: ${L.clNames[L.selCluster] || `cluster ${L.selCluster}`}`;
    return "";
  }
  const lassoByCluster = () => { const R = S.lasso && S.lasso.result; return !!(R && R.cluster && S.lasso.lpColor === "cluster"); };

  function lassoLayers() {
    const L = S.lasso;
    if (!L || S.mode === "rna" || inUmap() || (!L.shapes.length && !L.drawing)) return [];
    const e = currentExpr();
    const r = Math.max(0.5, (e && e.spotDiameterPx ? e.spotDiameterPx / 2 : 4) * S.pointScale);
    const out = [];
    if (L.bin && L.bin.n) out.push(new ScatterplotLayer({
      id: `lasso-bins-${S.exprKey}`, data: L.bin.data,
      radiusUnits: "common", getRadius: r, radiusMinPixels: 0.8, stroked: false, filled: true,
      parameters: NORMAL, updateTriggers: { getFillColor: [L.bin.key] },
    }));
    const R = L.result;
    if (R && L.hover >= 0 && R.info.mode === "bins" && L.pos) {
      const row = R.rows[L.hover];
      out.push(new ScatterplotLayer({
        id: "lasso-hover", data: [[L.pos[2 * row], L.pos[2 * row + 1]]], getPosition: (d) => d,
        radiusUnits: "common", getRadius: r * 2.2, radiusMinPixels: 7, filled: false, stroked: true,
        lineWidthUnits: "pixels", getLineWidth: 2.5, getLineColor: [255, 255, 255, 255], parameters: NORMAL,
      }));
    }
    const polys = L.shapes.map((s, i) => ({ s, i })).filter((d) => d.s.pts.length >= 3);   // selection shapes have no outline
    if (L.drawing) polys.push({ s: { pts: L.drawing.pts, hex: "#ffffff" }, i: -1 });
    out.push(new PolygonLayer({
      id: "lasso-shapes", data: polys, getPolygon: (d) => d.s.pts,
      getLineColor: (d) => [...lassoRgb(d.s), 255],
      getFillColor: (d) => [...lassoRgb(d.s), d.i >= 0 && d.i === L.selShape ? 45 : 0],
      lineWidthUnits: "pixels", getLineWidth: (d) => (d.i >= 0 && d.i === L.selShape ? 3.5 : 2),
      stroked: true, filled: true, pickable: true, parameters: NORMAL,
      updateTriggers: { getFillColor: [L.selShape, L.shapes.length], getLineWidth: [L.selShape], getLineColor: [L.shapes.length] },
    }));
    return out;
  }

  function lassoPickShape(i) {
    const L = S.lasso;
    L.selShape = i === null || L.selShape === i ? null : i;
    L.brush = null; L.selCluster = null;
    lassoRefresh(); buildLassoUI();
  }

  function buildLassoUI() {
    const L = S.lasso; if (!L) return;
    $("lassoCount").textContent = L.shapes.length ? `(${L.shapes.length})` : "";
    const list = $("lassoShapes"); list.innerHTML = "";
    const R = L.result;
    L.shapes.forEach((s, i) => {
      const row = document.createElement("div");
      row.className = "lshape" + (L.selShape === i ? " sel" : "");
      const notes = [];
      if (!s.bins) notes.push(`<small class="bad">0 bins — empty, ignored</small>`);
      else {
        let t = `${s.bins.toLocaleString()} bins`;
        if (s.shared) t += ` · ${s.shared.toLocaleString()} shared with earlier shapes, counted here`;
        if (s.lost) t += ` · ${s.lost.toLocaleString()} given to later shapes`;
        notes.push(`<small>${t}</small>`);
        if (R && R.info.found[i] !== undefined)
          notes.push(`<small>${R.info.kept[i].toLocaleString()} in UMAP` +
            (L.prot ? "" : ` · ${R.info.dropped[i].toLocaleString()} under ${LASSO.minUmi} UMI`) +
            (R.info.qcRemoved ? ` · ${R.info.qcRemoved[i].toLocaleString()} removed by QC` : "") + `</small>`);
        const n = R ? R.info.kept[i] : s.bins;
        if (n > 0 && n < LASSO.smallShape) notes.push(`<small class="warn">too few bins for this shape to mean much</small>`);
        const rr = L.prot && L.roi && L.roi.sig === lassoSig() && L.roi.key === S.exprKey ? L.roi.rows.find((x) => x.shape === i) : null;
        if (rr) notes.push(`<small class="roi">${roiLine(rr)}</small>`);
      }
      if (s.missing) notes.push(`<small class="bad">data missing for part of this shape</small>`);
      const nt = notes.join("");
      row.innerHTML = `<i style="background:${s.hex}"></i><div><span class="nm"></span>${L.prot ? toProtWords(nt) : isCells() ? toCellWords(nt) : isRois() ? toRoiWords(nt) : nt}</div>` +
        `<span class="lbtns"><button class="ren" title="rename / label">✎</button><button class="rm" title="remove">×</button></span>`;
      const nm = row.querySelector(".nm");
      nm.textContent = s.name;                     // typed by the user: never as HTML
      // A label is what the shape is ("fascicle", "epineurium"): it follows the shape into
      // the legend, the CSV and the saved JSON. Renaming does not change the UMAP.
      // A button, not a double-click: the first click of a double-click rebuilds this list.
      row.querySelector(".ren").onclick = (ev) => {
        ev.stopPropagation();
        const inp = document.createElement("input");
        inp.type = "text"; inp.value = s.name; inp.className = "nm-edit";
        let fin = false;                           // Escape, then the blur it causes, must not commit
        const done = (keep) => {
          if (fin) return; fin = true;
          if (keep && inp.value.trim()) s.name = inp.value.trim().slice(0, 60);
          buildLassoUI(); renderLassoPanel();
        };
        inp.onkeydown = (e) => { e.stopPropagation(); if (e.key === "Enter") done(true); if (e.key === "Escape") done(false); };
        inp.onblur = () => done(true);
        inp.onclick = (e) => e.stopPropagation();
        nm.replaceWith(inp); inp.focus(); inp.select();
      };
      row.querySelector(".rm").onclick = (ev) => {
        ev.stopPropagation();
        L.shapes.splice(i, 1);
        if (L.selShape === i) L.selShape = null; else if (L.selShape > i) L.selShape--;
        lassoShapesChanged();
      };
      row.onclick = () => lassoPickShape(i);
      list.appendChild(row);
    });
    const total = L.shapes.reduce((a, s) => a + s.bins, 0);
    if (!L.running && !R) {
      let m = L.shapes.length ? `${total.toLocaleString()} bins inside the shapes.` : "Pick Circle or Freehand, then draw on the slide.";
      const qcOut = L.shapes.length ? binQcInShapes() : 0;
      if (qcOut) m += ` The QC filter leaves out ${qcOut.toLocaleString()} of them.`;
      if (total > LASSO_CAP) m += ` Over ${LASSO_CAP.toLocaleString()}: a seeded sample of ${LASSO_CAP.toLocaleString()} will be used` +
        (TOUCH ? ` (a phone's limit; a computer takes up to ${LASSO.cap.toLocaleString()}).` : ".");
      const secs = Math.round(Math.min(total, LASSO_CAP) * 0.00095);
      if (TOUCH && total > 5000) m += ` On a phone this takes about ${secs}-${3 * secs} s: keep this page open until it finishes.`;
      else if (total > LASSO.softCap) m += ` This will take about ${secs} s.`;
      $("lassoMsg").textContent = L.prot ? toProtWords(m) : isCells() ? toCellWords(m) : isRois() ? toRoiWords(m) : m;
    }
    $("lassoRun").disabled = L.running || !L.shapes.some((s) => s.bins > 0);
    $("lassoCompare").hidden = !L.prot;
    $("lassoCompare").disabled = L.running || !L.shapes.some((s) => s.bins > 0);
    $("lassoSaveCsv").disabled = !R;
    buildRoiTable();
    if (window.SodaVolcano) window.SodaVolcano.refresh();
    if (window.SodaPlots) window.SodaPlots.refresh();
  }

  function lassoCancel() {
    const L = S.lasso;
    if (L && L.worker && L.running) { L.worker.terminate(); L.worker = null; }
    if (L) L.running = false;
    $("lassoCancel").hidden = true;
    lassoWake(false);
  }

  // Keep the screen on while a UMAP runs: a phone that locks freezes the page, and the run with
  // it. Refused (battery saver, plain http) or unsupported: the run simply goes on without it.
  let wakeLock = null, wakeWant = false;
  async function lassoWake(on) {
    if (!TOUCH) return;                                        // a desktop's screen settings are its own
    wakeWant = on;
    try {
      if (on && !wakeLock && navigator.wakeLock) {
        const w = await navigator.wakeLock.request("screen");
        w.addEventListener("release", () => { if (wakeLock === w) wakeLock = null; });
        if (wakeWant) wakeLock = w; else w.release();          // the run ended while we waited
      } else if (!on && wakeLock) { const w = wakeLock; wakeLock = null; await w.release(); }
    } catch (e) { /* not allowed here */ }
  }

  const lassoChunkOf = (x, y, cp) => `${Math.floor(x / cp)}_${Math.floor(y / cp)}`;

  // how far a run's clusters can be trusted (stability = median ARI over reseeded and 80%
  // subsampled Leiden runs, lasso_worker.js clusterGraph), with the sweep when there was one
  function lassoStabilityNote(I) {
    if (I.stability == null) return "";
    const sw = I.sweep && I.sweep.length > 1
      ? ` Sweep (resolution → clusters, stability): ${I.sweep.map((r) => `${r.resolution} → ${r.nClusters} (${r.stability.toFixed(2)})`).join(", ")}.` : "";
    if (I.stable) return `These clusters reproduce across reseeded and subsampled runs.${sw}`;
    if (I.nClusters < 2) return `No split of these bins into clusters reproduces (stability ≥ 0.8), so they read as one population.` +
      (I.closest ? ` Closest: ${I.closest.nClusters} clusters at ${I.closest.resolution} (${I.closest.stability.toFixed(2)}); type that resolution to see them.` : "") + sw;
    if (I.closest) return `No resolution gave clusters that reproduce (stability ≥ 0.8), so this shows the coarsest. ` +
      `Closest: ${I.closest.nClusters} clusters at ${I.closest.resolution} (${I.closest.stability.toFixed(2)}); type that resolution to see them.${sw}`;
    if (I.nClusters < 2) return `One cluster: no split of these bins reproduces.${sw}`;
    return `These clusters do not reproduce (stability below 0.8): treat them as a sketch.${sw}`;
  }

  // The marker lists behind the suggested cluster names (config markerSet, default
  // markers/nerve.json next to this page). Fetched once; missing or broken -> no suggestions.
  let markerSetP = null;
  function loadMarkerSet() {
    if (!markerSetP) markerSetP = (async () => {
      try {
        const r = await fetch(cfg.markerSet || "markers/rna_generic.json", { cache: "no-cache" });
        const j = r.ok ? await r.json() : null;
        S.markerSet = j && Array.isArray(j.types) && j.types.length ? j : null;
      } catch (e) { S.markerSet = null; }
      return S.markerSet;
    })();
    return markerSetP;
  }

  let protMarkerSetP = null;
  function loadProteinMarkerSet() {
    const e = currentExpr() || {}, ci = (S.umap && S.umap.ctInfo) || (e.umap && e.umap.celltypes), own = ci && ci.markerSet && ci.markerSet.file;
    const file = own && /^markers\/[\w.-]+\.json$/.test(own) ? own : "markers/protein_generic.json";
    if (!protMarkerSetP || protMarkerSetP.file !== file) {
      protMarkerSetP = (async () => {
        try {
          const r = await fetch(file, { cache: "no-cache" });
          const j = r.ok ? await r.json() : null;
          return j && j.modality === "protein" && Array.isArray(j.types) && j.types.length ? j : null;
        } catch (e) { return null; }
      })();
      protMarkerSetP.file = file;
    }
    return protMarkerSetP;
  }

  // ---- each shape against the rest of the slide (marker squares) ---------------------
  // A circle's own UMAP can only tell its squares apart from each other: a marker that is high
  // in every square drawn (CD45 over an immune-rich circle) differs nowhere inside it. So each
  // shape is also set against the rest of the slide, marker by marker, from the slide-wide
  // mean and sd that mat/index.json carries (older exports: worked out once from every chunk).
  const lassoSig = () => S.lasso.shapes.map((s) => `${s.id}:${s.bins}`).join(",");

  // the chunks with these keys, from the cache or the server (6 at a time)
  async function lassoFetch(keys) {
    const L = S.lasso, have = new Set(L.hvg.chunks.map((c) => `${c.tx}_${c.ty}`)), bufs = [], missing = [];
    let next = 0;
    const one = async () => {
      while (next < keys.length) {
        const k = keys[next++];
        if (L.chunks.has(k)) { bufs.push(L.chunks.get(k)); continue; }
        if (!have.has(k)) { missing.push(k); continue; }
        try { const b = await S.src.buf(`expr/${S.exprKey}/${L.prot ? "mat" : "hvg"}/${k}.bin`); L.chunks.set(k, b); bufs.push(b); }
        catch (e) { missing.push(k); }
      }
    };
    await Promise.all(Array.from({ length: 6 }, one));
    return { bufs, missing };
  }

  // x of a raw square mean v, channel j: the slide's own scaling (= lasso_worker.js xpOf)
  function protXp(H) {
    if (H.normalize === "bg-asinh-v1") {
      const cof = H.cofactor || 0.2, a1 = Math.asinh(1 / cof);
      return (v, j) => (H.hi[j] > 0 ? Math.asinh(Math.max(v - H.bg[j], 0) / H.hi[j] / cof) / a1 : 0);
    }
    return (v, j) => { const lo = H.p1[j], rg = H.p99_5[j] - lo; return rg > 0 ? Math.min(1, Math.max(0, (v - lo) / rg)) : 0; };
  }

  // each square of mat chunks: fn(row, values (Uint16, / scale), offset)
  function protEach(bufs, K, fn) {
    for (const buf of bufs) {
      const n = new DataView(buf).getUint32(0, true);
      const row = new Uint32Array(buf.slice(4, 4 + 4 * n)), q = new Uint16Array(buf.slice(4 + 4 * n, 4 + 4 * n + 2 * n * K));
      for (let i = 0; i < n; i++) fn(row[i], q, i * K);
    }
  }

  // The slide's mean and sd of x, and mean raw intensity, per marker. { error } when an older
  // export would need more than a phone should download, or some chunks would not load.
  async function lassoSlideStats(say) {
    const L = S.lasso, H = L.hvg, K = H.nChannels;
    if (Array.isArray(H.mean) && Array.isArray(H.sd) && Array.isArray(H.rawMean))
      return { mean: H.mean, sd: H.sd, rawMean: H.rawMean, n: H.nSquares || H.chunks.reduce((a, c) => a + c.n, 0), from: "index" };
    if (L.slideStats && L.slideStats.key === S.exprKey) return L.slideStats;
    const mb = H.chunks.reduce((a, c) => a + c.bytes, 0) / 1e6;
    if (TOUCH && mb > 30) return { error: "This export has no slide-wide averages, and working them out needs more than a phone should download: compare on a computer, or re-export." };
    say(`This export has no slide-wide averages yet: reading every square once (${mb.toFixed(0)} MB)…`);
    const all = H.chunks.map((c) => `${c.tx}_${c.ty}`);
    let { missing } = await lassoFetch(all);
    // a busy server can drop a few of hundreds of requests: those again, twice (the rest are cached)
    for (let t = 0; t < 2 && missing.length; t++) missing = (await lassoFetch(missing)).missing;
    if (missing.length) return { error: `${missing.length} of ${all.length} chunks would not load, so the slide-wide averages are unknown: try again.` };
    const bufs = all.map((k) => L.chunks.get(k));
    const xp = protXp(H), sc = H.scale, sx = new Float64Array(K), sxx = new Float64Array(K), sv = new Float64Array(K);
    let n = 0;
    protEach(bufs, K, (r, q, at) => {
      n++;
      for (let j = 0; j < K; j++) { const v = q[at + j] / sc, x = xp(v, j); sx[j] += x; sxx[j] += x * x; sv[j] += v; }
    });
    const mean = Array.from(sx, (a) => a / n);
    L.slideStats = { key: S.exprKey, mean, sd: mean.map((m, j) => Math.sqrt(Math.max(sxx[j] / n - m * m, 0))),
                     rawMean: Array.from(sv, (a) => a / n), n, from: "chunks" };
    return L.slideStats;
  }

  // per shape (and all shapes together): mean x and raw intensity inside, the same for the rest of
  // the slide, z = (inside - rest) / the slide's sd, and the ratio of raw means. Every square of
  // the shape counts (no sampling); a square in two shapes belongs to the one drawn last, as in the UMAP.
  function lassoRoiCompare(bufs, st) {
    const L = S.lasso, H = L.hvg, K = H.nChannels, xp = protXp(H), sc = H.scale, nS = L.shapes.length;
    const n = new Float64Array(nS + 1), sx = new Float64Array((nS + 1) * K), sv = new Float64Array((nS + 1) * K);
    protEach(bufs, K, (r, q, at) => {
      const s = L.assign[r]; if (s === LASSO_NONE || s >= nS) return;
      n[s]++; n[nS]++;
      for (let j = 0; j < K; j++) { const v = q[at + j] / sc, x = xp(v, j); sx[s * K + j] += x; sv[s * K + j] += v; sx[nS * K + j] += x; sv[nS * K + j] += v; }
    });
    const live = L.shapes.map((_, i) => i).filter((i) => n[i] > 0);
    const rows = [...live, ...(live.length > 1 ? [nS] : [])].map((i) => {
      const k = n[i], rest = Math.max(st.n - k, 1), x = [], raw = [], xRest = [], rawRest = [], z = [], fold = [];
      for (let j = 0; j < K; j++) {
        x.push(sx[i * K + j] / k); raw.push(sv[i * K + j] / k);
        xRest.push((st.mean[j] * st.n - sx[i * K + j]) / rest); rawRest.push((st.rawMean[j] * st.n - sv[i * K + j]) / rest);
        z.push(st.sd[j] > 0 ? (x[j] - xRest[j]) / st.sd[j] : 0);
        fold.push(rawRest[j] > 0 ? raw[j] / rawRest[j] : null);
      }
      return { shape: i === nS ? -1 : i, n: k, x, raw, xRest, rawRest, z, fold };
    });
    L.roi = { sig: lassoSig(), key: S.exprKey, channels: H.channels, nSlide: st.n, from: st.from, rows };
    if (L.roiPick === undefined || !rows.some((r) => r.shape === L.roiPick)) L.roiPick = rows.length ? rows[0].shape : null;
  }

  // "Compare with the slide" on its own: the chunks under the shapes and the table, no UMAP
  async function lassoCompareOnly() {
    const L = S.lasso;
    if (L.running || !L.prot || !L.hvg || !L.assign) return;
    const msg = (t, err) => { $("lassoMsg").textContent = t; $("lassoMsg").style.color = err ? "#ff8f8f" : ""; };
    if (!L.shapes.some((s) => s.bins > 0)) { msg("Draw a shape over tissue first.", true); return; }
    const keys = new Set();
    for (let i = 0; i < L.assign.length; i++) if (L.assign[i] !== LASSO_NONE) keys.add(lassoChunkOf(L.pos[2 * i], L.pos[2 * i + 1], L.hvg.chunkPx));
    L.running = true; buildLassoUI();
    try {
      msg(`Fetching ${keys.size} data chunk(s)…`);
      const { bufs, missing } = await lassoFetch([...keys]);
      const st = await lassoSlideStats(msg);
      if (st.error) { msg(st.error, true); return; }
      lassoRoiCompare(bufs, st);
      msg(missing.length ? `Compared with the rest of the slide; ${missing.length} chunk(s) could not be loaded, so part of a shape is missing.` : "Compared with the rest of the slide.");
    } finally { L.running = false; buildLassoUI(); }
  }

  // raised / lowered only when the raw ratio and z (background off, asinh) agree on the direction
  const roiAgree = (r, j) => r.fold[j] == null || (r.z[j] > 0) === (r.fold[j] > 1);

  function roiTxt(r, j) {
    const f = r.fold[j];
    return f == null ? `z ${r.z[j] >= 0 ? "+" : ""}${r.z[j].toFixed(1)}` : f >= 10 ? `×${f.toFixed(0)}` : `×${f.toFixed(1)}`;
  }

  // one line per shape in the list: the markers that set it apart most
  function roiLine(r) {
    const ch = S.lasso.roi.channels, ord = ch.map((_, j) => j).sort((a, b) => r.z[b] - r.z[a]);
    const up = ord.filter((j) => r.z[j] >= 0.5 && roiAgree(r, j)).slice(0, 4), down = ord.filter((j) => r.z[j] <= -0.5 && roiAgree(r, j)).reverse().slice(0, 2);
    if (!up.length && !down.length) return "vs rest of slide: no marker differs by half a standard deviation";
    return "vs rest of slide:" + (up.length ? " ↑ " + up.map((j) => `${escHtml(ch[j])} ${roiTxt(r, j)}`).join(", ") : "") +
      (down.length ? `${up.length ? " ·" : ""} ↓ ` + down.map((j) => `${escHtml(ch[j])} ${roiTxt(r, j)}`).join(", ") : "");
  }

  // the full table for the picked shape, every marker, most raised first
  function buildRoiTable() {
    const L = S.lasso, box = $("lassoRoi");
    const R = L.roi && L.roi.sig === lassoSig() && L.roi.key === S.exprKey ? L.roi : null;
    box.hidden = !L.prot || !R || !R.rows.length;
    if (box.hidden) return;
    const nameOf = (i) => (i === -1 ? "All shapes together" : L.shapes[i] ? L.shapes[i].name : `Shape ${i + 1}`);
    const sel = $("lassoRoiPick");
    sel.innerHTML = "";
    for (const r of R.rows) { const o = document.createElement("option"); o.value = r.shape; o.textContent = `${nameOf(r.shape)} (${r.n.toLocaleString()} squares)`; sel.appendChild(o); }
    sel.value = L.roiPick;
    const r = R.rows.find((x) => x.shape === L.roiPick) || R.rows[0];
    const ord = R.channels.map((_, j) => j).sort((a, b) => r.z[b] - r.z[a]);
    $("lassoRoiTable").innerHTML = `<table><thead><tr><th>Marker</th><th title="Mean intensity of the shape's squares (0-255, as the image)">Inside</th>` +
      `<th title="Mean intensity of every other tissue square on the slide">Rest</th><th title="Inside / rest">Ratio</th>` +
      `<th title="(inside - rest) in the slide's standard deviations, on the scaled values the UMAP uses">z</th></tr></thead><tbody>` +
      ord.map((j) => `<tr class="${Math.abs(r.z[j]) < 0.25 || !roiAgree(r, j) ? "dim" : r.z[j] > 0 ? "up" : "down"}"><td>${escHtml(R.channels[j])}</td>` +
        `<td>${protVal(r.raw[j])}</td><td>${protVal(r.rawRest[j])}</td><td>${r.fold[j] == null ? "–" : r.fold[j].toFixed(2)}</td>` +
        `<td>${r.z[j] >= 0 ? "+" : ""}${r.z[j].toFixed(2)}</td></tr>`).join("") + `</tbody></table>`;
    $("lassoRoiNote").textContent = `${r.n.toLocaleString()} squares inside, ${(R.nSlide - r.n).toLocaleString()} in the rest of the slide` +
      (R.from === "chunks" ? " (slide averages worked out in this browser: the export predates them). " : ". ") +
      `Ratio = mean intensity inside / rest. z = the difference on the scale the UMAP uses (${L.hvg.normalize === "bg-asinh-v1" ? "background off, asinh" : "1st–99.5th percentile window"}), ` +
      "in the slide's standard deviations. Grey: under 0.25, or ratio and z disagree on the direction.";
  }

  function lassoSaveRoiCsv() {
    const L = S.lasso, R = L.roi; if (!R) return;
    const q = (t) => `"${String(t).replace(/"/g, '""')}"`, f = (x, d) => (x == null ? "" : (+x).toFixed(d));
    const lines = ["dataset,layer,shape,squares_inside,squares_rest,marker,mean_inside,mean_rest,ratio,z_vs_rest,scaled_inside,scaled_rest"];
    for (const r of R.rows) {
      const nm = r.shape === -1 ? "All shapes together" : L.shapes[r.shape].name;
      R.channels.forEach((c, j) => lines.push([S.dsId, S.exprKey, q(nm), r.n, R.nSlide - r.n, q(c), f(r.raw[j], 3), f(r.rawRest[j], 3),
        f(r.fold[j], 4), f(r.z[j], 4), f(r.x[j], 5), f(r.xRest[j], 5)].join(",")));
    }
    lassoDownload(`${S.dsId}_${S.exprKey}_shapes_vs_slide.csv`, lines.join("\n") + "\n", "text/csv");
  }

  async function lassoRun() {
    const L = S.lasso;
    if (L.running || !L.hvg || !L.assign) return;
    const qcKeep = binQc && !L.prot ? binQc.keepMask() : null, qcState = binQc && !L.prot ? binQc.state() : null, qcTxt = qcKeep ? binQc.describe() : "", qcNote = L.qcNote || "";
    const msg = (t, err) => { $("lassoMsg").textContent = L.prot ? toProtWords(t) : isCells() ? toCellWords(t) : t; $("lassoMsg").style.color = err ? "#ff8f8f" : ""; };
    const mode = $("lassoMode").value;
    const seed = Number($("lassoSeed").value);
    if (!Number.isInteger(seed) || seed < 0) { msg("Seed must be a whole number ≥ 0.", true); return; }
    // "auto" (or empty) = the worker picks the finest resolution whose clusters reproduce
    const resTxt = $("lassoRes").value.trim(), resolution = resTxt === "" || /^auto$/i.test(resTxt) ? "auto" : Number(resTxt);
    if (resolution !== "auto" && (!Number.isFinite(resolution) || resolution < 0 || resolution > 10)) { msg("Resolution must be auto or a number from 0 to 10 (0 = no clusters).", true); return; }
    const clustering = resolution === "auto" || resolution > 0;
    const live = L.shapes.filter((s) => s.bins > 0);
    const total = live.reduce((a, s) => a + s.bins, 0);
    if (!live.length) { msg("Every shape is empty — draw over tissue.", true); return; }
    if (mode === "bins" && total < LASSO.minPoints) { msg(`UMAP needs at least ${LASSO.minPoints} bins; the shapes hold ${total}.`, true); return; }
    if (mode === "shapes" && live.length < 3) { msg(`"One per shape" needs at least 3 non-empty shapes; there are ${live.length}.`, true); return; }

    // the chunks under the shapes, and which shapes each one feeds
    const cp = L.hvg.chunkPx;
    const need = new Map();
    for (let i = 0; i < L.assign.length; i++) {
      const s = L.assign[i]; if (s === LASSO_NONE) continue;
      const k = lassoChunkOf(L.pos[2 * i], L.pos[2 * i + 1], cp);
      if (!need.has(k)) need.set(k, new Set());
      need.get(k).add(s);
    }
    L.shapes.forEach((s) => (s.missing = false));
    L.running = true; $("lassoCancel").hidden = false; buildLassoUI();
    lassoWake(true);
    msg(`Fetching ${need.size} data chunk(s)…`);
    const { bufs, missing } = await lassoFetch([...need.keys()]);
    if (!L.running) return;                         // cancelled while fetching
    if (missing.length) {
      missing.forEach((k) => need.get(k).forEach((s) => (L.shapes[s].missing = true)));
      L.running = false; buildLassoUI();
      if (!confirm(`${missing.length} of ${need.size} data chunk(s) could not be loaded, so part of a shape has no data. Run anyway on what did load?`)) {
        msg("Not run: data missing for part of a shape.", true); $("lassoCancel").hidden = true; lassoWake(false); return;
      }
      L.running = true;
    }
    // gene names for the marker lists (column order of the hvg chunks); markers fall back to
    // column numbers if the file is missing
    if (!L.hvgGenes) { try { L.hvgGenes = await S.src.json(`expr/${S.exprKey}/hvg/genes.json`); } catch (e) { L.hvgGenes = null; } }
    const markerSet = mode === "bins" && clustering ? await (L.prot ? loadProteinMarkerSet() : loadMarkerSet()) : null;
    if (!L.running) return;                         // cancelled meanwhile
    if (!L.worker) L.worker = new Worker("lasso_worker.js?v=20261111q");
    const sig = L.shapes.map((s) => `${s.id}:${s.bins}`).join(",");
    const params = { mode, seed, minUmi: LASSO.minUmi, cap: LASSO_CAP, minPoints: LASSO.minPoints,
                     nPcs: LASSO.nPcs, nPcsProtein: LASSO.nPcsProtein, nNeighbors: LASSO.nNeighbors, minDist: LASSO.minDist, minShapesForUmap: LASSO.minShapesForUmap,
                     minDetect: LASSO.minDetect, minDetectFrac: LASSO.minDetectFrac, minDetectCap: LASSO.minDetectCap, resolution };
    const t0 = performance.now();
    L.worker.onmessage = (ev) => {
      const d = ev.data;
      if (d.type === "progress") { msg(`${d.step}… ${Math.round(d.frac * 100)}%`); return; }
      L.running = false; $("lassoCancel").hidden = true; lassoWake(false);
      if (d.type === "error") { msg(d.message, true); buildLassoUI(); return; }
      if (sig !== L.shapes.map((s) => `${s.id}:${s.bins}`).join(",")) { buildLassoUI(); return; }   // shapes moved on
      const idxOfRow = new Map(); d.rows.forEach((r, j) => idxOfRow.set(r, j));
      d.info.mode = mode;
      // Leiden clusters ride along when the worker made them; an older worker, "one dot per
      // shape" or resolution 0 leave them out, and the panel simply colours by shape
      const cluster = d.cluster && d.cluster.length === d.rows.length ? d.cluster : null;
      L.result = { emb: d.emb, rows: d.rows, shape: d.shape, cluster, info: d.info, idxOfRow, seed, params, layer: S.exprKey,
                   resolution: cluster ? (d.info.resolution != null ? d.info.resolution : resolution) : null,
                   secs: (performance.now() - t0) / 1000, id: Date.now(), qc: qcState };
      if (cluster && !d.info.nClusters) { let m = 0; for (const c of cluster) if (c !== LASSO_NONE && c + 1 > m) m = c + 1; d.info.nClusters = m; }
      L.brush = null; L.selShape = null; L.hover = -1; L.lpVs = null; L.selCluster = null; L.clNames = {};
      const I = d.info, drop = I.dropped.reduce((a, b) => a + b, 0);
      msg([`${I.method === "pca" ? "PCA" : "UMAP"} of ${I.nIn.toLocaleString()} ${mode === "bins" ? "bins" : "bins in " + d.rows.length + " shapes"} ` +
           `in ${L.result.secs.toFixed(1)} s · ${I.nGenesUsed} genes · ${I.nPcs} PCs${I.nNeighbors ? ` · ${I.nNeighbors} neighbours` : ""} · seed ${seed}`,
           drop ? (L.prot ? `Left out ${drop.toLocaleString()} low-signal squares (no marker above ${Math.round((I.lowSignalX || 0.2) * 100)}% of its range: rims and gaps).`
                          : `Dropped ${drop.toLocaleString()} bins with < ${LASSO.minUmi} UMI.`) : "",
           I.sampledFrom ? `Showing ${I.nIn.toLocaleString()} of ${I.sampledFrom.toLocaleString()} bins (seeded sample, stratified by shape).` : "",
           cluster && I.nClusters < 2 ? `Leiden found one cluster (resolution ${L.result.resolution}).` :
           cluster ? `${I.nClusters} Leiden clusters at resolution ${L.result.resolution}${I.autoResolution ? " (chosen by stability)" : ""}` +
             (I.stability != null ? `, stability ${I.stability.toFixed(2)}.` : ".") :
             mode === "bins" && clustering ? "No clusters came back (this worker version does not cluster yet)." : "",
           cluster ? lassoStabilityNote(I) : "",
           L.prot && L.roiError ? `Not compared with the rest of the slide: ${L.roiError}` : "",
           I.qcRemoved ? `QC filter: left out ${I.qcRemoved.reduce((a, b) => a + b, 0).toLocaleString()} bins (${qcTxt}).` : "",
           qcNote,
          ].filter(Boolean).join("\n"));
      $("lassoPanel").hidden = false;
      lassoRefresh(); buildLassoUI();
      requestAnimationFrame(() => { lpDeck && lpDeck.setProps({ width: "100%", height: "100%" }); fitLassoPanel(); lassoFrameShapes(); });
    };
    const H = L.hvg;
    let st = null;
    if (L.prot) {
      st = await lassoSlideStats(msg);
      if (!L.running) return;                       // cancelled meanwhile
      if (st.error) { L.roiError = st.error; st = null; }
      else { L.roiError = null; lassoRoiCompare(bufs, st); buildLassoUI(); }
    }
    msg("Computing…");
    const protein = L.prot ? { channels: H.channels, p1: H.p1, p99_5: H.p99_5, median: H.median, scale: H.scale,
                               normalize: H.normalize || null, bg: H.bg, hi: H.hi, cofactor: H.cofactor,
                               mean: st ? st.mean : null, sd: st ? st.sd : null } : null;
    const flags = L.prot && S.qc.flags && S.qc.key === S.exprKey ? S.qc.flags : null;   // for the "Artefact?" name
    L.worker.postMessage({ chunks: bufs, rowShape: L.assign, nShapes: L.shapes.length, nGenes: L.hvg.nGenes, params, resolution,
                           geneNames: L.hvgGenes || null, markerSet, protein, flags, qcKeep });
    L.qcNote = null; $("bqStale").hidden = true;
  }

  // ---- the UMAP panel: its own small deck --------------------------------------------
  function lpEnsureDeck() {
    if (lpDeck) return lpDeck;
    lpDeck = new Deck({
      parent: $("lpHost"),
      views: [new OrthographicView({ id: "lp", flipY: false,
        controller: { dragPan: false, dragRotate: false, doubleClickZoom: false, keyboard: false, scrollZoom: { speed: 0.01, smooth: false } } })],
      viewState: { target: [0, 0, 0], zoom: 0 },
      layers: [], useDevicePixels: DEVICE_PX,
      deviceProps: { webgl: { preserveDrawingBuffer: true } },
      getCursor: () => "crosshair",
      onViewStateChange: ({ viewState }) => { S.lasso.lpVs = viewState; lpDeck.setProps({ viewState }); },
      onHover: (info) => {
        const L = S.lasso;
        const i = info && info.layer && info.layer.id.startsWith("lp-pts") ? info.index : -1;
        if (i !== L.hover) { L.hover = i; requestRender(); renderLassoPanel(); }
      },
    });
    // drag = select (a box, or a loop with the Lasso tool; Shift or "+" adds to what is
    // selected); a click picks the point's shape or cluster; a click on nothing clears
    const host = $("lpHost"), box = $("lpBox"), svg = $("lpSvg");
    let start = null, path = null, add = false;
    const local = (ev) => { const r = host.getBoundingClientRect(); return [ev.clientX - r.left, ev.clientY - r.top]; };
    const free = () => S.lasso.lpTool === "free";
    const clearDraw = () => { box.hidden = true; svg.innerHTML = ""; };
    host.addEventListener("pointerdown", (ev) => {
      if (ev.button !== 0 || !S.lasso.result || start) return;
      start = local(ev); path = [start]; add = ev.shiftKey || !!S.lasso.lpAdd;
      host.setPointerCapture(ev.pointerId);
    });
    // a second finger (a pinch, which the panel's own deck zooms) or the browser taking the touch over
    host.addEventListener("pointercancel", () => { start = null; path = null; clearDraw(); });
    host.addEventListener("pointermove", (ev) => {
      if (!start) return;
      const p = local(ev);
      if (free()) {
        const q = path[path.length - 1];
        if (Math.hypot(p[0] - q[0], p[1] - q[1]) >= 3) path.push(p);
        svg.innerHTML = `<polygon points="${path.map((x) => x.join(",")).join(" ")}" />`;
        return;
      }
      Object.assign(box.style, { left: `${Math.min(p[0], start[0])}px`, top: `${Math.min(p[1], start[1])}px`,
        width: `${Math.abs(p[0] - start[0])}px`, height: `${Math.abs(p[1] - start[1])}px` });
      box.hidden = false;
    });
    host.addEventListener("pointerup", (ev) => {
      if (!start) return;
      const p = local(ev), s = start, L = S.lasso, R = L.result, loop = free() ? path : null;
      start = null; path = null; clearDraw();
      if (!R) return;
      let far = 0; for (const q of loop || [p]) far = Math.max(far, Math.hypot(q[0] - s[0], q[1] - s[1]));
      if (far < 4) {
        const hit = lpDeck.pickObject({ x: p[0], y: p[1], radius: 6 });
        if (hit && hit.layer && hit.layer.id.startsWith("lp-pts")) {
          if (lassoByCluster()) lassoPickCluster(R.cluster[hit.index]); else lassoPickShape(R.shape[hit.index]);
        }
        else { L.brush = null; L.selShape = null; L.selCluster = null; lassoRefresh(); buildLassoUI(); }
        return;
      }
      if (loop && loop.length < 3) return;
      const vp = lpDeck.getViewports()[0];
      if (!vp) return;                              // not drawn yet
      const pts = loop || [s, p];
      let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity;
      for (const [x, y] of pts) { x0 = Math.min(x0, x); y0 = Math.min(y0, y); x1 = Math.max(x1, x); y1 = Math.max(y1, y); }
      const n = R.emb.length / 2, keep = add && L.brush && L.brush.length === n;
      const mask = keep ? L.brush.slice() : new Uint8Array(n);
      for (let j = 0; j < n; j++) {
        const [sx, sy] = vp.project([R.emb[2 * j], R.emb[2 * j + 1]]);
        if (!(sx >= x0 && sx <= x1 && sy >= y0 && sy <= y1)) continue;
        if (!loop || inPoly(sx, sy, loop)) mask[j] = 1;
      }
      let k = 0; for (let j = 0; j < n; j++) k += mask[j];
      L.brush = k ? mask : null; L.selShape = null; L.selCluster = null; L.hinted = true;
      $("lpFoot").dataset.brush = k ? `${k.toLocaleString()} point(s) selected — lit on the slide.` : "";
      lassoRefresh(); buildLassoUI();
    });
    return lpDeck;
  }

  // The panel sits over the right of the slide: frame the shapes in what is left of it,
  // so both sides of the link are on screen at once.
  function lassoFrameShapes() {
    const L = S.lasso; if (!L.shapes.length || !S.viewState || inUmap()) return;
    let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity;
    const grow = (x, y) => { x0 = Math.min(x0, x); y0 = Math.min(y0, y); x1 = Math.max(x1, x); y1 = Math.max(y1, y); };
    for (const s of L.shapes) {
      for (const [x, y] of s.pts) grow(x, y);
      if (s.rows && L.pos) for (const i of s.rows) grow(L.pos[2 * i], L.pos[2 * i + 1]);
    }
    if (!(x1 > x0)) return;
    const [W, H] = umapPaneSize("ortho"), lp = $("lassoPanel");
    // side by side, the panel sits over the UMAP pane, not over the slide; otherwise it covers
    // the right of the slide, or on a phone its lower part
    const open = !lp.hidden && !umapSbs(), below = open && lp.offsetWidth > W * 0.8;
    const free = Math.max(200, W - (open && !below ? lp.offsetWidth + 20 : 0)), freeH = Math.max(150, H - (below ? lp.offsetHeight + 50 : 0));
    const zoom = Math.min(1, Math.log2(Math.min(free / ((x1 - x0) * 1.3 || 1), freeH / ((y1 - y0) * 1.3 || 1))));
    const scale = Math.pow(2, zoom);
    S.viewState = { ...S.viewState, target: [(x0 + x1) / 2 + (W - free) / 2 / scale, (y0 + y1) / 2 + (H - freeH) / 2 / scale, 0], zoom };
    autoFit = false;
    deckgl.setProps({ viewState: deckVs() }); updateHud(); requestRender();
  }

  function fitLassoPanel() {
    const R = S.lasso.result; if (!R) return;
    lpEnsureDeck();
    const host = $("lpHost"), w = host.clientWidth || 300, h = host.clientHeight || 300;
    let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity;
    for (let j = 0; j < R.emb.length; j += 2) {
      x0 = Math.min(x0, R.emb[j]); x1 = Math.max(x1, R.emb[j]); y0 = Math.min(y0, R.emb[j + 1]); y1 = Math.max(y1, R.emb[j + 1]);
    }
    const dx = Math.max(x1 - x0, 1e-6), dy = Math.max(y1 - y0, 1e-6);
    S.lasso.lpVs = { target: [(x0 + x1) / 2, (y0 + y1) / 2, 0], zoom: Math.log2(Math.min(w / (dx * 1.15), h / (dy * 1.15))) };
    lpDeck.setProps({ viewState: S.lasso.lpVs });
    renderLassoPanel();
  }

  function renderLassoPanel() {
    const L = S.lasso, R = L && L.result;
    if (!R || $("lassoPanel").hidden) return;
    lpEnsureDeck();
    const n = R.emb.length / 2, colors = new Uint8Array(4 * n), rgb = L.shapes.map(lassoRgb);
    const byCl = lassoByCluster();
    const pick = L.brush ? (j) => L.brush[j] === 1 : L.selShape !== null ? (j) => R.shape[j] === L.selShape
      : L.selCluster !== null && R.cluster ? (j) => R.cluster[j] === L.selCluster : null;
    for (let j = 0; j < n; j++) {
      const c = byCl ? clusterRgb(R.cluster[j]) : rgb[R.shape[j]] || [200, 200, 200];
      colors[4 * j] = c[0]; colors[4 * j + 1] = c[1]; colors[4 * j + 2] = c[2];
      colors[4 * j + 3] = pick ? (pick(j) ? 255 : 30) : 220;
    }
    const big = R.info.mode === "shapes";
    const layers = [new ScatterplotLayer({
      id: `lp-pts-${R.id}`, data: { length: n, attributes: { getPosition: { value: R.emb, size: 2 }, getFillColor: { value: colors, size: 4, normalized: true } } },
      radiusUnits: "pixels", getRadius: big ? 7 : n > 20000 ? 1.5 : 2.2, stroked: big, getLineColor: [255, 255, 255, 200],
      lineWidthUnits: "pixels", getLineWidth: 1, pickable: true, parameters: NORMAL,
      updateTriggers: { getFillColor: [Date.now()] },
    })];
    const seg = $("lpColorBy");
    seg.hidden = !R.cluster;
    seg.querySelectorAll("button").forEach((b) => b.classList.toggle("on", b.dataset.c === (byCl ? "cluster" : "shape")));
    if (L.hover >= 0 && L.hover < n) layers.push(new ScatterplotLayer({
      id: "lp-hover", data: [[R.emb[2 * L.hover], R.emb[2 * L.hover + 1]]], getPosition: (d) => d,
      radiusUnits: "pixels", getRadius: big ? 11 : 6, filled: false, stroked: true, lineWidthUnits: "pixels", getLineWidth: 2,
      getLineColor: [255, 255, 255, 255], parameters: NORMAL,
    }));
    lpDeck.setProps({ layers, viewState: L.lpVs || { target: [0, 0, 0], zoom: 0 } });
    const I = R.info;
    $("lpTitle").textContent = `${I.method === "pca" ? "PCA (PC1 vs PC2)" : "UMAP"} of drawn shapes · ${R.info.mode === "bins" ? `one dot per ${wd().bin}` : "one dot per shape"}`;
    const warns = [];
    if (I.method === "pca") warns.push(`Fewer than ${LASSO.minShapesForUmap} shapes: a UMAP of so few points is meaningless, so this is PC1 vs PC2 (${(I.varFrac[0] * 100).toFixed(0)}% / ${((I.varFrac[1] || 0) * 100).toFixed(0)}% of variance).`);
    else if (I.mode === "bins" && I.nIn < LASSO.fewPoints) warns.push(`Few points (${I.nIn}): treat this layout as a sketch.`);
    if (R.cluster && I.stability != null && !I.stable) warns.push(lassoStabilityNote({ ...I, sweep: null }));   // numbers only: safe as HTML
    if (I.sampledFrom) warns.push(`Showing ${I.nIn.toLocaleString()} of ${I.sampledFrom.toLocaleString()} bins (seeded sample, stratified by shape).`);
    $("lpWarn").innerHTML = warns.join("<br>"); $("lpWarn").hidden = !warns.length;
    const esc = (x) => String(x).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
    if (byCl) {
      // one row per cluster with its top markers; click lights its bins on the slide
      const sizes = I.clusterSizes || (() => { const z = new Array(I.nClusters).fill(0); for (const c of R.cluster) if (c < z.length) z[c]++; return z; })();
      const ct = I.celltypes && Array.isArray(I.celltypes.clusters) ? I.celltypes : null;
      const sug = (c) => (ct && ct.clusters[c]) || null;
      const vsSlide = !!I.markersVsSlide && L.mkVs !== "shapes";
      $("lpMkVs").hidden = !I.markersVsSlide;
      $("lpMkVs").querySelectorAll("button").forEach((b) => b.classList.toggle("on", b.dataset.v === (vsSlide ? "slide" : "shapes")));
      clusterLegend($("lpLegend"), {
        n: I.nClusters, sizes, markers: (vsSlide ? I.markersVsSlide : I.markers) || [], picked: L.selCluster,
        markerTitle: L.prot ? (vsSlide
          ? "Top markers against the whole slide (the cluster's mean, in the slide's standard deviations; share of its squares above the marker's typical level):\n"
          : "Top markers against the other squares drawn (how much higher than in the rest of the shapes, in their standard deviations; share of its squares above the marker's typical level):\n") : null,
        name: (c) => L.clNames[c] || `Cluster ${c}`,
        onPick: (c) => lassoPickCluster(c),
        onRename: (c, t) => { if (t) L.clNames[c] = t; renderLassoPanel(); },
        suggest: ct ? sug : null,
        markerInfo: () => ({ ...ctMarkerInfo(ct, L.prot
          ? `Scored on the ${I.nIn.toLocaleString()} squares of this run, on the ${L.hvg ? L.hvg.nChannels : ""} markers of the section UMAP, ` +
            (ct && ct.scale === "slide" ? "each against the whole slide's average and spread." : "each against the other squares drawn (this export has no slide-wide averages).")
          : `Scored on the ${I.nIn.toLocaleString()} bins of this run, using only the section's ` +
            `${L.hvg ? L.hvg.nGenes.toLocaleString() : ""} most variable genes: marker genes outside them are listed as not found.`) }),
        onAccept: (c, t) => { L.clNames[c] = t; renderLassoPanel(); },
        acceptAll: (count) => {
          const todo = [];
          for (let c = 0; c < I.nClusters; c++) { const s = sug(c); if (s && s.confidence === "high" && !L.clNames[c]) todo.push([c, s.suggest]); }
          if (!count) { todo.forEach(([c, t]) => (L.clNames[c] = t)); renderLassoPanel(); }
          return todo.length;
        },
      });
    } else {
      $("lpMkVs").hidden = true;
      $("lpLegend").innerHTML = L.shapes.map((s, i) => I.kept[i] ? `<div class="leg"><i style="background:${s.hex}"></i><span>${esc(s.name)}</span><small>${I.kept[i].toLocaleString()} ${wd().bins}</small></div>` : "").join("");
    }
    $("lpFoot").textContent = ($("lpFoot").dataset.brush && L.brush ? $("lpFoot").dataset.brush + " " : "") +
      (L.prot ? `${I.nGenesUsed} markers · intensity ${L.hvg.normalize === "bg-asinh-v1" ? "with the slide's background off, scaled to its 99.5th percentile, asinh" : "scaled to the section's 1st–99.5th percentile"}, z-scored here · ${I.nPcs} PCs · `
              : `${I.nGenesUsed} variable genes (of the section's ${L.hvg ? L.hvg.nGenes : "?"}) · log1p CP10k · ${I.nPcs} PCs · `) +
      `${I.method === "umap" ? `umap-js ${LASSO.umapJs}, ${I.nNeighbors || LASSO.nNeighbors} neighbours, min_dist ${LASSO.minDist}${I.init ? `, ${I.init.toUpperCase()} start` : ""}, ` : ""}seed ${R.seed}` +
      `${R.cluster ? ` · Leiden resolution ${R.resolution}${I.autoResolution ? " (auto)" : ""}, ${I.nClusters} clusters` +
        `${I.stability != null ? `, stability ${I.stability.toFixed(2)}` : ""}` : ""}. ` +
      `Drag ${L.lpTool === "free" ? "a loop" : "a box"} to select (Shift or + adds) · click a point to pick its ${byCl ? "cluster" : "shape"} · ` +
      `Show on slide frames the picked squares · scroll to zoom.`;
  }

  function lassoPickCluster(c) {
    const L = S.lasso; L.hinted = true;
    L.selCluster = c === null || c === undefined || L.selCluster === c ? null : c;
    L.brush = null; L.selShape = null;
    lassoRefresh(); buildLassoUI();
  }

  // ---- save / load ------------------------------------------------------------------
  function lassoDownload(name, text, type) {
    const a = document.createElement("a");
    a.href = URL.createObjectURL(new Blob([text], { type }));
    a.download = name; document.body.appendChild(a); a.click(); a.remove();
    setTimeout(() => URL.revokeObjectURL(a.href), 1000);
  }
  function lassoSaveJson() {
    const L = S.lasso, R = L.result;
    const out = {
      tool: "SODA circle -> UMAP", version: 1, umapJs: LASSO.umapJs, dataset: S.dsId, layer: S.exprKey,
      seed: Number($("lassoSeed").value), mode: $("lassoMode").value,
      params: L.prot
        ? { modality: "protein", channels: L.hvg ? L.hvg.channels : null, nPcs: LASSO.nPcsProtein, nNeighbors: LASSO.nNeighbors,
            minDist: LASSO.minDist, cap: R ? R.params.cap : LASSO_CAP, normalisation: L.hvg ? L.hvg.normalize : null, scaling: "z-score within the shapes, clip 10",
            overlap: "drawn last wins" }
        : { minUmi: LASSO.minUmi, nPcs: LASSO.nPcs, nNeighbors: LASSO.nNeighbors, minDist: LASSO.minDist,
            minDetect: LASSO.minDetect, minDetectFrac: LASSO.minDetectFrac, minDetectCap: LASSO.minDetectCap,
            cap: R ? R.params.cap : LASSO_CAP, sectionHvgs: L.hvg ? L.hvg.nGenes : null, normalisation: "log1p CP10k", overlap: "drawn last wins" },
      shapes: L.shapes.map((s) => ({ name: s.name, color: s.hex, kind: s.kind,
                                      points: s.pts.map((p) => [Math.round(p[0] * 100) / 100, Math.round(p[1] * 100) / 100]),
                                      ...(s.rows ? { rows: Array.from(s.rows) } : {}) })),
      resolution: /^\s*(auto)?\s*$/i.test($("lassoRes").value) ? "auto" : Number($("lassoRes").value) || 0,
      result: R ? { method: R.info.method, points: R.info.nIn, genesUsed: R.info.nGenesUsed, minDetect: R.info.minDetect || null,
                    nPcs: R.info.nPcs, nNeighbors: R.info.nNeighbors || null, init: R.info.init || null, seconds: +R.secs.toFixed(2),
                    ...(R.cluster ? { resolution: R.resolution, nClusters: R.info.nClusters, clusterSizes: R.info.clusterSizes || null,
                                      autoResolution: !!R.info.autoResolution, stability: R.info.stability != null ? R.info.stability : null,
                                      clusterStability: R.info.clusterStability || null, sweep: R.info.sweep || null,
                                      stabilityRuns: R.info.stabilityRuns || null,
                                      clusterNames: Object.keys(L.clNames).length ? { ...L.clNames } : null,
                                      // suggested names from the marker lists, with the genes found / missing
                                      cellTypes: R.info.celltypes || null,
                                      markers: R.info.markers || null,
                                      markersVsSlide: R.info.markersVsSlide || null,
                                      // bin row (positions.bin order, never a barcode) and its cluster
                                      rows: Array.from(R.rows), cluster: Array.from(R.cluster) } : {}) } : null,
      // the whole-section clustering, if one is shown: its resolution and any names given to clusters
      sectionClusters: umapSaveState(),
      // the per-bin QC thresholds the result was computed with (null = no filter)
      qcFilter: R ? R.qc || null : binQc ? binQc.state() : null,
      saved: new Date().toISOString(),
    };
    lassoDownload(`${S.dsId}_${S.exprKey}_shapes.json`, JSON.stringify(out, null, 1), "application/json");
  }
  function lassoSaveCsv() {
    const L = S.lasso, R = L.result; if (!R) return;
    const q = (x) => `"${String(x).replace(/"/g, '""')}"`;
    const lines = [];
    if (R.info.mode === "bins") {
      const cl = R.cluster;
      lines.push("bin_row,x_px,y_px,shape,umap_x,umap_y" + (cl ? ",cluster,cluster_name,leiden_resolution" : ""));
      for (let j = 0; j < R.rows.length; j++) {
        const r = R.rows[j];
        lines.push(`${r},${L.pos[2 * r].toFixed(2)},${L.pos[2 * r + 1].toFixed(2)},${q(L.shapes[R.shape[j]].name)},${R.emb[2 * j]},${R.emb[2 * j + 1]}` +
          (cl ? `,${cl[j] === LASSO_NONE ? "" : cl[j]},${cl[j] === LASSO_NONE ? "" : q(L.clNames[cl[j]] || `Cluster ${cl[j]}`)},${R.resolution}` : ""));
      }
    } else {
      lines.push("shape,bins,x,y");
      for (let j = 0; j < R.rows.length; j++) lines.push(`${q(L.shapes[R.shape[j]].name)},${R.info.kept[R.shape[j]]},${R.emb[2 * j]},${R.emb[2 * j + 1]}`);
    }
    lassoDownload(`${S.dsId}_${S.exprKey}_${R.info.method}.csv`, lines.join("\n") + "\n", "text/csv");
  }
  async function lassoLoadJson(file) {
    let j;
    try { j = JSON.parse(await file.text()); } catch (e) { $("lassoMsg").textContent = "Not a JSON file."; return; }
    if (!Array.isArray(j.shapes)) { $("lassoMsg").textContent = "No shapes in that file."; return; }
    if (j.dataset && j.dataset !== S.dsId && !confirm(`These shapes were drawn on ${j.dataset}, not ${S.dsId}. Load them anyway?`)) return;
    if (j.layer && j.layer !== S.exprKey && [...$("exprSelect").options].some((o) => o.value === j.layer)) {
      $("exprSelect").value = j.layer;
      if (!(await setExprKey(j.layer))) return;    // another layer or dataset was picked while it loaded
    }
    const L = S.lasso;
    L.shapes = j.shapes.map((s, i) => ({ id: i + 1, name: s.name || `Shape ${i + 1}`, kind: s.kind || "free", pts: s.points || [],
                                         ...(Array.isArray(s.rows) ? { rows: Uint32Array.from(s.rows) } : {}),
                                         hex: s.color || LASSO_COLORS[i % LASSO_COLORS.length], bins: 0, shared: 0, lost: 0 }));
    L.next = L.shapes.length + 1; L.selShape = null;
    if (Number.isInteger(j.seed)) $("lassoSeed").value = j.seed;
    if (Number.isFinite(j.resolution) || j.resolution === "auto") $("lassoRes").value = j.resolution;
    umapLoadState(j.sectionClusters);
    L.qcNote = binQc ? binQc.setState(j.qcFilter || null) : null;   // a file without qcFilter = no filter
    if (j.mode === "bins" || j.mode === "shapes") $("lassoMode").value = j.mode;
    lassoShapesChanged();
    lassoRun();
  }

  function lassoInit() {
    S.lasso = lassoBlank();
    lassoPadEvents();
    $("lassoCircle").onclick = () => lassoSetTool("circle");
    $("lassoFree").onclick = () => lassoSetTool("free");
    $("lassoClear").onclick = () => { S.lasso.shapes = []; S.lasso.next = 1; S.lasso.selShape = null; lassoShapesChanged(); };
    $("lassoRun").onclick = lassoRun;
    $("lassoCompare").onclick = lassoCompareOnly;
    $("lassoRoiPick").onchange = (e) => { S.lasso.roiPick = Number(e.target.value); buildRoiTable(); };
    $("lassoRoiCsv").onclick = lassoSaveRoiCsv;
    $("lpMkVs").onclick = (e) => { const v = e.target.dataset && e.target.dataset.v; if (!v) return; S.lasso.mkVs = v; renderLassoPanel(); };
    $("lassoCancel").onclick = () => { lassoCancel(); $("lassoMsg").textContent = "Cancelled."; buildLassoUI(); };
    $("lassoMode").onchange = () => { if (S.lasso.result) lassoRun(); };
    $("lassoSaveJson").onclick = lassoSaveJson;
    $("lassoSaveCsv").onclick = lassoSaveCsv;
    $("lassoLoadBtn").onclick = () => $("lassoLoad").click();
    $("lassoLoad").onchange = (e) => { const f = e.target.files[0]; e.target.value = ""; if (f) lassoLoadJson(f); };
    $("lpFit").onclick = fitLassoPanel;
    $("lpColorBy").onclick = (e) => {
      const c = e.target.dataset && e.target.dataset.c; if (!c) return;
      S.lasso.lpColor = c; S.lasso.selCluster = null; lassoRefresh();
    };
    $("lpClose").onclick = () => { $("lassoPanel").hidden = true; lassoClearPick(); };
    $("lpSelShow").onclick = () => { const L = S.lasso, g = lassoLitRows(); if (g) showRowsOnSlide((r) => L.assign[r] !== LASSO_NONE && g(r, L.assign[r]), L.pos); };
    $("lpSelClear").onclick = () => lassoClearPick();
    $("lpTool").onclick = (e) => {
      const t = e.target.dataset && e.target.dataset.t; if (!t) return;
      S.lasso.lpTool = t; $("lpTool").querySelectorAll("button").forEach((b) => b.classList.toggle("on", b.dataset.t === t));
    };
    $("lpAdd").onclick = () => { S.lasso.lpAdd = !S.lasso.lpAdd; $("lpAdd").classList.toggle("on", S.lasso.lpAdd); };
    // the panel's plot follows its box: side by side, the phone drawer, a resized window
    if (window.ResizeObserver) {
      let last = "";
      new ResizeObserver(() => {
        const host = $("lpHost"), k = `${host.clientWidth}x${host.clientHeight}`;
        if (k === last || !host.clientWidth || !host.clientHeight) return;
        last = k;
        if (lpDeck && S.lasso && S.lasso.result) { lpDeck.setProps({ width: "100%", height: "100%" }); fitLassoPanel(); }
      }).observe($("lpHost"));
    }
    window.addEventListener("keydown", (e) => {
      if (e.key !== "Escape" || !S.lasso) return;
      if (S.lasso.tool) lassoSetTool(null);
      else if (S.lasso.brush || S.lasso.selShape !== null || S.lasso.selCluster !== null) {
        S.lasso.brush = null; S.lasso.selShape = null; S.lasso.selCluster = null; lassoRefresh(); buildLassoUI();
      }
    });
    window.addEventListener("resize", () => { if (S.lasso && S.lasso.result) renderLassoPanel(); });
    // for testing from the console: restore shapes saved with Save JSON, without a file
    window.__spatialViewer.lassoLoad = (obj) => lassoLoadJson({ text: async () => JSON.stringify(obj) });
    // volcano plots (volcano.js, de_worker.js): handed only what they need from here
    if (window.SodaVolcano) window.SodaVolcano.init({
      state: () => S, touch: TOUCH, fetchChunks: lassoFetch, chunkOf: lassoChunkOf, addGene, download: lassoDownload,
      say: (t) => (isRois() ? toRoiWords(t) : t),
      geneXyv: (name) => S.src.gene(S.exprKey, name),
      hvgGenes: async () => { const L = S.lasso; if (!L.hvgGenes) { try { L.hvgGenes = await S.src.json(`expr/${S.exprKey}/hvg/genes.json`); } catch (e) { return null; } } return L.hvgGenes; },
      // the section's clustering on show (Leiden or the lab's own), with the names given in the Clusters panel
      sectionClusters: () => {
        if (!umapAvail() || !S.umap || S.umap.res == null) return null;
        const c = umapClusterings().find((x) => clKey(x) === S.umap.res);
        if (!c) return null;
        const names = { ...(c.categories || []) }, own = (S.umap.names && S.umap.names[clKey(c)]) || {};
        Object.assign(names, own);
        return { n: c.n, res: clKey(c), label: clLabel(c), sizes: c.sizes || null, names, user: clUser(c) };
      },
      loadSectionClusters: async () => { await umapSetRes(S.umap.res); return S.umap.cl ? { cl: S.umap.cl, res: S.umap.clRes } : null; },
    });
    // heatmap / dot plot / violin (plots.js, plots_worker.js): handed only what they need from here
    if (window.SodaPlots) window.SodaPlots.init({
      state: () => S, touch: TOUCH, download: lassoDownload, addGene, addRnaGene, rnaColorBy: setColorBy, fetchChunks: lassoFetch, chunkOf: lassoChunkOf,
      qc: () => binQc, clusterings: umapClusterings, clusterRgb, catRgb: (i) => liftForDark(catColor(i)), lightRows: umapLightRows,
    });
    buildLassoUI();
  }

  // ---- Tissue <-> UMAP: the whole section's own embedding ------------------------------
  // export.py embeds every bin of the 8 / 16 um layers that has >= 10 UMI and clusters them
  // with Leiden (expr/<key>/umap/). The main view shows the bins where they sit on the slide
  // (Tissue), the same bins laid out by that UMAP (UMAP), or both at once (Side by side).
  // Every view draws rows of positions.bin, so a gene, a highlighted cluster or a selection
  // is one set of bins and shows everywhere.
  // The UMAP is scaled into the tissue's own bounding box, so one camera serves both and
  // each dot can fly from its place on the slide to its place in the UMAP.
  const UMAP_NONE = 0xffff;
  const UMAP_GREY = [92, 99, 115], UMAP_LIT_GREY = [232, 236, 244], UMAP_SEL = [255, 214, 64];
  // how long the flight takes: the "Speed" menu beside Play (slow ones are for talks and
  // for following one region by eye); remembered per browser
  const FLY_MS = 1300, FLY_KEY = "sv.flyMs";
  const flyMs = () => Number($("flySpeed") && $("flySpeed").value) || FLY_MS;
  const umapBlank = () => ({ key: null, index: null, pos: null, uxy: null, box: null, nEmb: 0, radius: 1,
    clusters: new Map(), cl: null, clRes: null, res: null, centroids: new Map(), names: {}, hidden: new Set(),
    pick: null, sel: null, selN: 0, selEmb: 0, view: "tissue", colorBy: "gene", tool: null, hover: -1, grid: null,
    tissueVs: null, umapVs: null, sideVs: null, anim: null, lastFly: null, rowOf: null, mix: null, ver: 0,
    cloud: new Map(), drawn: new Map(), busy: null,
    // suggested cell types: celltypes.json (per clustering), celltype.bin (per bin), a picked type
    ctInfo: null, ctTried: false, ct: null, ctCounts: null, ctPick: null,
    // gating: a highlighted phenotype (rule index, or GATE_NONE for Unassigned)
    phPick: null });
  const inUmap = () => !!(S.umap && S.umap.view === "umap" && S.mode !== "rna");
  const umapSbs = () => !!(S.umap && S.umap.view === "sbs" && S.mode !== "rna");
  const umapMeta = () => { const e = S.m && S.mode !== "rna" ? currentExpr() : null; return e && e.umap ? e.umap : null; };
  const umapAvail = () => !!(umapMeta() && S.src && S.src.json && S.src.buf);
  const umapByCluster = () => !!(S.umap && S.umap.colorBy === "cluster" && S.umap.cl);
  const umapClusterings = () => {
    const U = S.umap, m = umapMeta();
    return (U && U.index && U.index.clusterings) || (m && m.clusterings) || [];
  };
  // a clustering is known by its "id" (the lab's own, userclusters.py) or its Leiden resolution
  const clKey = (c) => (c ? (c.id != null ? c.id : c.resolution) : null);
  const clUser = (c) => !!(c && c.source === "user");
  const clLabel = (c) => (c ? c.name || `Leiden ${c.resolution}` : "");
  const clFromValue = (v) => { const c = umapClusterings().find((x) => String(clKey(x)) === String(v)); return c ? clKey(c) : null; };
  const umapClusterInfo = () => umapClusterings().find((c) => clKey(c) === S.umap.clRes) || null;
  // a name given in the legend wins, then the lab's own category name, then "Cluster N"
  const umapName = (c) => {
    const own = (S.umap.names[S.umap.clRes] || {})[c];
    if (own) return own;
    const info = umapClusterInfo();
    return info && info.categories && info.categories[c] != null ? String(info.categories[c]) : `Cluster ${c}`;
  };
  // tab20, lifted for the dark background (see liftForDark); fixed per cluster id, so cluster 3
  // is the same colour in every view, the legend and the circle panel
  const clusterLut = [];
  const clusterRgb = (c) => (c === UMAP_NONE || c == null ? [120, 120, 120] : clusterLut[c] || (clusterLut[c] = liftForDark(catColor(c))));
  // the section's clustering on show: the lab's own colours when it brought them (index.json "colors")
  let secLut = { key: null, lut: null };
  const sectionRgb = (c) => {
    if (c === UMAP_NONE || c == null) return [120, 120, 120];
    const U = S.umap, key = `${U.key}|${U.clRes}`;
    if (secLut.key !== key) {
      const info = umapClusterInfo(), cols = info && Array.isArray(info.colors) ? info.colors : null;
      secLut = { key, lut: cols ? cols.map((h) => liftForDark(hexToRgb255(h))) : null };
    }
    return (secLut.lut && secLut.lut[c]) || clusterRgb(c);
  };
  function umapDirty() { if (S.umap) S.umap.ver++; }

  // ---- the flight: one ScatterplotLayer, two positions per dot, mixed on the GPU --------
  // instancePositions is where the bin sits on the slide, instanceFlyTargets.xy where it sits
  // in the UMAP (z = 1 if it is in the UMAP at all). fly.t moves every dot at once, so an
  // animation frame changes one uniform and the camera - no buffer is rebuilt. Bins that are
  // not in the UMAP (< 10 UMI) stay put and fade out.
  const flyModule = { name: "fly", vs: "uniform flyUniforms { float t; } fly;", uniformTypes: { t: "f32" } };
  class FlyScatterLayer extends ScatterplotLayer {
    getShaders() {
      const s = super.getShaders();
      const cut = s.vs.indexOf("void main(void) {");
      const head = s.vs.slice(0, cut).replace("in vec3 instancePositions;", "in vec3 instancePositions;\nin vec3 instanceFlyTargets;");
      const body = s.vs.slice(cut + "void main(void) {".length)
        .replace(/instancePositions\b/g, "flyPos")
        .replace(/vFillColor = vec4\(instanceFillColors\.rgb, instanceFillColors\.a \* layer\.opacity\);/,
                 "$&\nvFillColor.a *= mix(1.0, instanceFlyTargets.z, fly.t);");
      const vs = head + "void main(void) {\nvec3 flyPos = vec3(mix(instancePositions.xy, instanceFlyTargets.xy, fly.t * instanceFlyTargets.z), instancePositions.z);" + body;
      return { ...s, vs, modules: [...(s.modules || []), flyModule] };
    }
    initializeState() {
      super.initializeState();
      this.getAttributeManager().addInstanced({ instanceFlyTargets: { size: 3, accessor: "getFlyTarget" } });
    }
    draw(opts) {
      this.setShaderModuleProps({ fly: { t: this.props.flyT } });
      super.draw(opts);
    }
  }
  FlyScatterLayer.layerName = "FlyScatterLayer";
  FlyScatterLayer.defaultProps = {
    ...ScatterplotLayer.defaultProps,
    flyT: { type: "number", value: 0 },
    getFlyTarget: { type: "accessor", value: [0, 0, 0] },
  };

  function umapReset() {
    if (S.umap) umapSetTool(null);
    S.umap = umapBlank();
    umapApplyViews();
    syncUmapUI();
  }

  // the bins' slide positions: whichever copy is already in memory for this layer
  function umapPositions() {
    const e = currentExpr(), n2 = e ? e.nBins * 2 : -1;
    if (S.tissue.positions && S.tissue.key === S.exprKey && S.tissue.positions.length === n2) return Promise.resolve(S.tissue.positions);
    if (S.lasso && S.lasso.pos && S.lasso.pos.length === n2) return Promise.resolve(S.lasso.pos);
    return S.src.positions(S.exprKey);
  }

  // Fetched only when someone asks for the UMAP or for clusters, never on open.
  function umapLoad() {
    const U = S.umap;
    if (U.uxy) return Promise.resolve();
    if (U.busy) return U.busy;
    U.busy = (async () => {
      status("Loading the UMAP…");
      try {
        const [index, buf, pos] = await Promise.all([
          S.src.json(`expr/${U.key}/umap/index.json`),
          S.src.buf(`expr/${U.key}/umap/coords.bin`),
          U.pos ? U.pos : umapPositions()]);
        const xy = new Float32Array(buf), N = pos.length / 2;
        if (xy.length !== 2 * N) throw new Error(`coords.bin has ${xy.length / 2} rows, positions.bin ${N}`);
        let nEmb = 0;
        let tx0 = Infinity, ty0 = Infinity, tx1 = -Infinity, ty1 = -Infinity;
        const ex = new Float32Array(N), ey = new Float32Array(N);
        for (let i = 0; i < N; i++) {
          const px = pos[2 * i], py = pos[2 * i + 1];
          if (px < tx0) tx0 = px; if (px > tx1) tx1 = px; if (py < ty0) ty0 = py; if (py > ty1) ty1 = py;
          const x = xy[2 * i], y = xy[2 * i + 1];
          if (!(x === x) || !(y === y)) continue;                  // NaN: under 10 UMI, not embedded
          ex[nEmb] = x; ey[nEmb] = y; nEmb++;
        }
        // an ROI layer (geomx.py) has groupings and no embedding: only the tissue view exists
        if (!nEmb && index.embedding && index.embedding.source === "none") {
          U.index = index; U.pos = pos; U.uxy = new Float32Array(2 * N).fill(NaN); U.nEmb = 0; U.radius = 1;
          U.box = [tx0, ty0, Math.max(tx1, tx0 + 1), Math.max(ty1, ty0 + 1)];
          const gl = index.clusterings || [];
          if (U.res == null && gl.length) U.res = clKey(gl[Math.min(gl.length - 1, index.default || 0)]);
          return;
        }
        if (!nEmb) throw new Error(`no ${wd().bin} is embedded`);
        // frame on the middle 99%: a few stray bins far from the cloud would otherwise shrink
        // it to a corner of the picture (they are still drawn, just outside the first frame)
        const sx = ex.subarray(0, nEmb).sort(), sy = ey.subarray(0, nEmb).sort();
        const q = (a, f) => a[Math.min(nEmb - 1, Math.max(0, Math.round(f * (nEmb - 1))))];
        const x0 = q(sx, 0.005), x1 = q(sx, 0.995), y0 = q(sy, 0.005), y1 = q(sy, 0.995);
        // into the tissue's bounding box, aspect kept: the camera hardly moves in the flight
        const bw = Math.max(x1 - x0, 1e-6), bh = Math.max(y1 - y0, 1e-6);
        const tw = Math.max(tx1 - tx0, 1), th = Math.max(ty1 - ty0, 1);
        const sc = Math.min(tw / bw, th / bh) * 0.96;
        const cx = (tx0 + tx1) / 2, cy = (ty0 + ty1) / 2, mx = (x0 + x1) / 2, my = (y0 + y1) / 2;
        const uxy = new Float32Array(2 * N);
        for (let i = 0; i < N; i++) {
          const x = xy[2 * i], y = xy[2 * i + 1];
          if (!(x === x) || !(y === y)) { uxy[2 * i] = NaN; uxy[2 * i + 1] = NaN; continue; }
          uxy[2 * i] = cx + (x - mx) * sc;
          uxy[2 * i + 1] = cy - (y - my) * sc;                    // UMAP y up = up on screen
        }
        U.index = index; U.pos = pos; U.uxy = uxy; U.nEmb = nEmb;
        U.box = [cx - (bw * sc) / 2, cy - (bh * sc) / 2, cx + (bw * sc) / 2, cy + (bh * sc) / 2];
        // about a third of the mean spacing: dense areas still read as dense, sparse dots stay dots
        U.radius = Math.sqrt((bw * sc) * (bh * sc) / nEmb) * 0.35;
        const list = index.clusterings || [];
        if (U.res == null && list.length) U.res = clKey(list[Math.min(list.length - 1, index.default || 0)]);
      } finally { status(null); U.busy = null; }
    })();
    return U.busy;
  }

  async function umapSetRes(res) {
    const U = S.umap;
    await umapLoad();
    if (S.umap !== U) return;
    const list = umapClusterings();
    const c = list.find((x) => clKey(x) === res) || list[Math.min(list.length - 1, (U.index && U.index.default) || 0)];
    if (!c) throw new Error("this layer has no clusterings");
    const ck = clKey(c);
    let arr = U.clusters.get(ck);
    if (!arr) {
      status("Loading clusters…");
      try { arr = new Uint16Array(await S.src.buf(`expr/${U.key}/umap/${c.file || `clusters_r${c.resolution}.bin`}`)); }
      finally { status(null); }
      if (S.umap !== U) return;
      if (arr.length !== U.pos.length / 2) throw new Error(`clusters file has ${arr.length} rows, positions.bin ${U.pos.length / 2}`);
      U.clusters.set(ck, arr);
    }
    if (U.clRes !== ck) { U.pick = null; U.hidden = new Set(); }
    U.res = ck; U.cl = arr; U.clRes = ck;
    await umapLoadCtInfo();                         // suggested names, if the export has them
    umapDirty();
  }

  // Called when the expression layer changes: everything here belongs to one layer's bins.
  async function umapOnLayer() {
    const old = S.umap, was = old.view;
    umapSetTool(null);
    const U = S.umap = umapBlank();
    U.colorBy = old.colorBy === "phenotype" ? "gene" : old.colorBy; U.tissueVs = was === "umap" ? old.tissueVs : S.viewState;
    if (umapAvail()) {
      U.key = S.exprKey;
      const m = umapMeta(), list = m.clusterings || [];
      if (list.length) U.res = clKey(list[Math.min(list.length - 1, m.default || 0)]);
      if (U.colorBy === "celltype" && !umapCtAvail()) U.colorBy = "gene";
      if (was !== "tissue" || U.colorBy !== "gene") {
        try {
          await umapLoad();
          if (S.umap === U && list.length) await umapSetRes(U.res);
          if (S.umap === U && U.colorBy === "celltype") await umapLoadTypes();
        } catch (e) { if (S.umap === U) { status(`Could not load the UMAP for this layer: ${e.message}`, true); setTimeout(() => status(null), 5000); } }
        if (S.umap !== U) return;
        if (was !== "tissue" && U.uxy && U.nEmb) U.view = was;     // stay in UMAP / side by side
      }
    }
    if (was === "umap") S.viewState = U.view === "umap" ? null : U.tissueVs;   // null: render() refits
    umapApplyViews();
    if (U.view === "sbs") fitView();
    syncUmapUI();
  }

  // side by side on a screen taller than wide (a phone held upright): the tissue on top and
  // the UMAP below, rather than two slivers next to each other
  const sbsStacked = () => { const el = $("viewport"); return el.clientHeight > el.clientWidth; };
  // the size of one pane: the whole picture, or a half of it side by side. The UMAP half gives
  // its lower part to the circle -> UMAP panel when that is open; stacked, the panel covers
  // the whole UMAP pane instead, as a quarter of a phone is too small to read either in.
  function umapPaneSize(id) {
    const el = $("viewport"), W = el.clientWidth || 800, H = el.clientHeight || 600;
    if (!umapSbs()) return [W, H];
    if (sbsStacked()) return [W, H / 2];
    if (id === "uside") return [W / 2, $("lassoPanel").hidden ? H : H / 2];
    return [W / 2, H];
  }
  function umapFitVs(w, h) {
    if (w == null) [w, h] = umapPaneSize("ortho");
    const b = S.umap.box;
    const zoom = Math.log2(Math.min(w / ((b[2] - b[0]) * 1.12), h / ((b[3] - b[1]) * 1.18)));
    return { target: [(b[0] + b[2]) / 2, (b[1] + b[3]) / 2, 0], zoom, minZoom: zoom - 2, maxZoom: zoom + 10 };
  }
  function imageFitVs(w, h) {
    if (w == null) [w, h] = umapPaneSize("ortho");
    const img = S.m.image;
    const zoom = Math.log2(Math.min(w / img.width, h / img.height)) - 0.03;
    return { target: [img.width / 2, img.height / 2, 0], zoom, minZoom: zoom - 2, maxZoom: 6 };
  }
  function umapApplyViews() {
    const main = $("viewport").parentElement;
    if (!deckgl || S.mode === "rna") { main.classList.remove("sbs", "stack"); return; }
    const sbs = umapSbs(), st = sbs && sbsStacked(), lp = !$("lassoPanel").hidden;
    main.classList.toggle("sbs", sbs);
    main.classList.toggle("stack", st);
    main.classList.toggle("lp-open", lp);
    const views = st
      ? [new OrthographicView({ id: "ortho", x: 0, y: 0, width: "100%", height: "50%", flipY: true, controller: ORTHO_CTL }),
         new OrthographicView({ id: "uside", x: 0, y: "50%", width: "100%", height: "50%", flipY: true, controller: ORTHO_CTL })]
      : sbs
      ? [new OrthographicView({ id: "ortho", x: 0, y: 0, width: "50%", height: "100%", flipY: true, controller: ORTHO_CTL }),
         new OrthographicView({ id: "uside", x: "50%", y: 0, width: "50%", height: lp ? "50%" : "100%",
                                flipY: true, controller: ORTHO_CTL })]
      : imageViews();
    deckgl.setProps(S.viewState ? { views, viewState: deckVs() } : { views });
  }

  // force: the person asked to watch the flight (Play), so it runs even with reduced motion
  async function umapSetView(v, force) {
    const U = S.umap;
    if (!U || U.anim || v === U.view || (v !== "tissue" && (!U.key || !umapAvail() || umapMeta().embedding === "none"))) return;
    if (v !== "tissue") {
      try {
        await umapLoad();
        // the clusters come along so hovering a dot can say which cluster it is in
        if (!U.cl && umapClusterings().length) await umapSetRes(U.res);
      } catch (e) { status(`Could not load the UMAP: ${e.message}`, true); setTimeout(() => status(null), 5000); return; }
      if (S.umap !== U) return;
    }
    const from = U.view;
    if (from === "umap") U.umapVs = S.viewState; else U.tissueVs = S.viewState;
    umapSetTool(null);
    if (v === "umap") lassoSetTool(null);
    U.view = v;
    umapApplyViews();
    syncUmapUI();
    if (v === "sbs") {
      // static: both panes framed afresh (the tissue keeps where you were if it fits)
      const [w2, h2] = umapPaneSize("uside");
      // from the tissue: the same place, one step out, as the pane is now half as wide
      S.viewState = from === "tissue" && S.viewState ? { ...S.viewState, zoom: S.viewState.zoom - 1 } : imageFitVs();
      U.sideVs = U.sideVs || umapFitVs(w2, h2);
      deckgl.setProps({ viewState: deckVs() }); updateHud(); render();
      return;
    }
    if (v === "umap") umapFly(true, S.viewState || imageFitVs(), U.umapVs || umapFitVs(), force);
    else if (from === "umap") umapFly(false, S.viewState, U.tissueVs || imageFitVs());
    else { S.viewState = U.tissueVs || imageFitVs(); deckgl.setProps({ viewState: deckVs() }); updateHud(); render(); }
  }

  // "Play transition": always tissue -> UMAP, ending in the UMAP view.
  async function umapPlay() {
    const U = S.umap;
    if (!U || U.anim) return;
    if (U.view !== "umap") { umapSetView("umap", true); return; }
    umapFly(true, U.tissueVs || imageFitVs(), S.viewState || umapFitVs(), true);
  }

  // The flight, eased (1.3 s, or slower from the Speed menu): dots move from the slide to the UMAP (or back), the H&E fades,
  // the camera moves between the two framings. Instant for reduced motion (unless asked for
  // with Play) or a hidden tab.
  function umapFly(toUmap, vsFrom, vsTo, force) {
    const U = S.umap, vs0 = vsFrom || S.viewState || vsTo;
    const done = () => {
      if (U.anim !== "fly") return;                               // already finished
      U.anim = null;
      if (S.umap !== U) return;
      S.viewState = vsTo; deckgl.setProps({ viewState: deckVs() }); updateHud(); syncUmapUI(); render();
    };
    U.anim = "fly";
    const reduce = window.matchMedia && window.matchMedia("(prefers-reduced-motion: reduce)").matches;
    if ((reduce && !force) || document.hidden) { done(); return; }
    const F = umapCloud("fly");
    if (!F.n) { done(); return; }
    $("densityKey").hidden = true;
    const t0 = performance.now(), ms = flyMs();
    let frames = 0;
    const step = (now) => {
      if (S.umap !== U || U.anim !== "fly") return;               // the layer changed mid-flight
      const f = Math.min(1, (now - t0) / ms);
      const ez = f < 0.5 ? 4 * f * f * f : 1 - Math.pow(-2 * f + 2, 3) / 2;
      umapFlyFrame(F, toUmap ? ez : 1 - ez, vs0, vsTo, ez);
      frames++;
      if (f < 1) requestAnimationFrame(step);
      else { U.lastFly = { frames, ms: now - t0, fps: (frames * 1000) / Math.max(1, now - t0), bins: F.n }; done(); }
    };
    requestAnimationFrame(step);
    // a tab sent to the background stops animation frames: never leave the view mid-flight
    setTimeout(() => { if (U.anim === "fly" && S.umap === U) done(); }, ms + 1500);
  }
  // one frame: a = 0 slide layout .. 1 UMAP layout; e = camera progress from vs0 to vs1
  function umapFlyFrame(F, a, vs0, vs1, e) {
    const U = S.umap, ex = currentExpr();
    const rT = Math.max(0.5, (ex && ex.spotDiameterPx ? ex.spotDiameterPx / 2 : 4) * S.pointScale), rU = U.radius * S.pointScale;
    const vs = { target: [0, 1].map((i) => vs0.target[i] + (vs1.target[i] - vs0.target[i]) * e).concat(0),
                 zoom: vs0.zoom + (vs1.zoom - vs0.zoom) * e, minZoom: -40, maxZoom: 40 };
    S.viewState = vs;
    const layers = [];
    const op = S.imgOpacity * Math.max(0, 1 - a * 1.6);           // the H&E is gone a little before the dots land
    if (S.imgOn && op > 0.01) {
      const nOn = S.channels.filter((ch) => ch.visible).length;
      for (const ch of S.channels) if (ch.visible) layers.push(channelLayer(ch, op, nOn));
    }
    layers.push(new FlyScatterLayer({
      id: `ufly-${U.key}`, data: F.data, flyT: a,
      radiusUnits: "common", getRadius: 1, radiusScale: rT + (rU - rT) * a, radiusMinPixels: 0.8,
      stroked: false, filled: true, parameters: NORMAL,
    }));
    deckgl.setProps({ viewState: deckVs(), layers });
  }

  // ---- what is lit: a selection, a picked cluster, or (in a UMAP) the circle tool's pick
  function umapLit(where) {
    const U = S.umap;
    if (!U || !U.pos) return null;
    if (U.sel) return (r) => U.sel[r] === 1;
    if (U.pick !== null && U.cl) return (r) => U.cl[r] === U.pick;
    if (umapCtLit()) return (r) => U.ct[r] === U.ctPick;
    if (umapPhLit()) { const calls = S.gate.calls, k = U.phPick; return (r) => calls[r] === k; }
    // gating preview: the squares above the cut-off light up in the UMAP (on the slide the
    // preview draws its own outlines)
    if (where !== "tissue" && gatePreviewOn()) { const v = S.gate.vals.get(S.gate.cur), t = S.gate.cut.get(S.gate.cur).v; return (r) => v[r] > t; }
    // a shape, brushed points or a cluster picked in Circle -> UMAP light up in the UMAP too;
    // on the slide the circle tool draws its own highlight, and the section's cluster colours
    // under it dim (same palette: undimmed, panel cluster 0 and section cluster 0 read as one)
    const L = S.lasso;
    if (L && L.assign && L.assign.length * 2 === U.pos.length && (where !== "tissue" || umapByCat())) {
      const f = lassoLitRows();
      if (f) return (r) => L.assign[r] !== LASSO_NONE && f(r, L.assign[r]);
    }
    return null;
  }

  // ---- genes per bin, for colouring the UMAP's dots -------------------------------------
  // Gene files hold only the bins where the gene is found, at the same coordinates as
  // positions.bin, so a bin's row is found from its coordinates (as region linking does).
  const ukey = (x, y) => Math.round(x * 100) * 1e8 + Math.round(y * 100);
  function umapRowOf() {
    const U = S.umap;
    if (!U.rowOf) {
      const m = new Map(), p = U.pos;
      for (let r = 0; r < p.length / 2; r++) m.set(ukey(p[2 * r], p[2 * r + 1]), r);
      U.rowOf = m;
    }
    return U.rowOf;
  }
  function geneByRow(g) {
    const U = S.umap;
    if (g.byRow && g.byRowKey === U.key) return g.byRow;
    const m = umapRowOf(), out = new Float32Array(U.pos.length / 2);      // 0 = not found in that bin
    for (let i = 0; i < g.n; i++) {
      const r = m.get(ukey(g.positions[2 * i], g.positions[2 * i + 1]));
      if (r !== undefined) out[r] = g.values[i];
    }
    g.byRow = out; g.byRowKey = U.key;
    return out;
  }
  // the bin under a point on the slide (hover in the side-by-side tissue pane): a hash of
  // bin centres one pitch wide, nearest centre within the 3x3 cells around the point
  function umapBinAt(x, y) {
    const U = S.umap, e = currentExpr();
    if (!U.pos) return -1;
    const cell = (e && e.spotDiameterPx) || 8;
    if (!U.grid) {
      const g = new Map(), p = U.pos;
      for (let r = 0; r < p.length / 2; r++) {
        const k = Math.floor(p[2 * r] / cell) * 1e6 + Math.floor(p[2 * r + 1] / cell);
        const a = g.get(k); if (a) a.push(r); else g.set(k, [r]);
      }
      U.grid = g;
    }
    const cx = Math.floor(x / cell), cy = Math.floor(y / cell);
    let best = -1, bd = (cell * 0.75) ** 2;
    for (let i = -1; i <= 1; i++) for (let j = -1; j <= 1; j++) {
      const a = U.grid.get((cx + i) * 1e6 + cy + j); if (!a) continue;
      for (const r of a) { const d = (U.pos[2 * r] - x) ** 2 + (U.pos[2 * r + 1] - y) ** 2; if (d < bd) { bd = d; best = r; } }
    }
    return best;
  }
  // One colour per bin from the genes on show, mixed the way the slide's density map mixes
  // them (inks that multiply, or the average on one ramp). The window's top is the gene's
  // max display value; t = how strong (-1 = none of these genes in this bin).
  function umapGeneMix() {
    const U = S.umap, on = S.genes.filter((g) => g.visible !== false);
    if (!on.length) return null;
    const key = on.map((g) => `${g.name}:${g.hex}:${g.min}:${g.max}`).join("|") + `#${S.genePalette}#${U.key}`;
    if (U.mix && U.mix.key === key) return U.mix;
    const N = U.pos.length / 2, t = new Float32Array(N).fill(-1), rgb = new Uint8Array(3 * N);
    if (S.genePalette === "genes") {
      const paper = new Float32Array(3 * N).fill(1), miss = new Float32Array(N).fill(1);
      // marker squares: each one's place in its marker's min..max display window, as on the
      // slide; the hue mixes the inks relative to the strongest marker, strength is the opacity
      const prot = isProt(), strength = prot ? on.map((g) => {
        const br = geneByRow(g), lo = g.min || 0, inv = 1 / Math.max(g.max - lo, 1e-6), s = new Float32Array(N);
        for (let r = 0; r < N; r++) { const v = br[r]; if (v > lo) s[r] = Math.min(1, (v - lo) * inv); }
        return s;
      }) : null;
      const smax = prot ? new Float32Array(N) : null;
      if (prot) for (const s of strength) for (let r = 0; r < N; r++) if (s[r] > smax[r]) smax[r] = s[r];
      for (const [gi, g] of on.entries()) {
        const br = geneByRow(g), [tr, tg, tb] = hexToRgb255(g.hex), lo = g.min || 0, hi = Math.max(g.max, 1e-6);
        for (let r = 0; r < N; r++) {
          const v = br[r];
          if (prot) {
            const s = strength[gi][r]; if (!(s > 0)) continue;
            const k = s / smax[r];
            paper[3 * r] *= 1 - k * (1 - tr / 255); paper[3 * r + 1] *= 1 - k * (1 - tg / 255); paper[3 * r + 2] *= 1 - k * (1 - tb / 255);
            miss[r] *= 1 - s;
            continue;
          }
          if (!(v > 0) || v < lo) continue;
          const sEff = 0.2 + 0.8 * Math.min(1, v / hi);
          paper[3 * r] *= 1 - sEff * (1 - tr / 255); paper[3 * r + 1] *= 1 - sEff * (1 - tg / 255); paper[3 * r + 2] *= 1 - sEff * (1 - tb / 255);
          miss[r] *= 1 - sEff;
        }
      }
      for (let r = 0; r < N; r++) {
        if (miss[r] >= 1) continue;
        const a = Math.sqrt(1 - miss[r]); t[r] = a;
        for (let c = 0; c < 3; c++) rgb[3 * r + c] = Math.round(UMAP_GREY[c] + (paper[3 * r + c] * 255 - UMAP_GREY[c]) * a);
      }
    } else {
      const lut = PALETTES[S.genePalette] || VIRIDIS, sum = new Float32Array(N), any = new Uint8Array(N);
      for (const g of on) {
        const br = geneByRow(g), lo = g.min || 0, hi = Math.max(g.max, 1e-6);
        if (g.prot) {
          const inv = 1 / Math.max(g.max - lo, 1e-6);
          for (let r = 0; r < N; r++) { const v = br[r]; if (v > lo) { sum[r] += Math.min(1, (v - lo) * inv); any[r] = 1; } }
          continue;
        }
        for (let r = 0; r < N; r++) { const v = br[r]; if (v > 0 && v >= lo) { sum[r] += Math.min(1, v / hi); any[r] = 1; } }
      }
      for (let r = 0; r < N; r++) {
        if (!any[r]) continue;
        const f = Math.min(1, sum[r] / on.length), k = Math.min(255, Math.round((0.06 + 0.94 * f) * 255)) * 3;
        t[r] = f; rgb[3 * r] = lut[k]; rgb[3 * r + 1] = lut[k + 1]; rgb[3 * r + 2] = lut[k + 2];
      }
    }
    U.mix = { key, t, rgb };
    return U.mix;
  }

  // The dots to draw: which rows, in what order (lit and strongly expressing last, so on
  // top), where and in what colour. Binary attributes throughout - 345k bins at 8 um.
  //   where = "umap":   every embedded bin, at its UMAP place
  //   where = "tissue": colour by cluster -> every clustered bin over the H&E;
  //                     colour by genes   -> only the lit bins (the density map stays below)
  //   where = "fly":    every bin, slide place + UMAP place (for FlyScatterLayer)
  function umapCloud(where) {
    const U = S.umap, N = U.pos.length / 2, lit = umapLit(where), byCl = umapByCluster(), byCt = !byCl && umapByType();
    const byPh = !byCl && !byCt && umapByPh(), calls = byPh ? S.gate.calls : null;
    const mix = byCl || byCt || byPh ? null : umapGeneMix();
    const key = `${where}|${U.key}|${U.ver}|${byCl ? U.clRes : byCt ? "ct" : byPh ? "ph" + S.gate.ver : "g"}|${mix ? mix.key : ""}|${[...U.hidden].join(",")}|${lit ? 1 : 0}`;
    const hit = U.cloud.get(where);
    if (hit && hit.key === key) return hit;
    if (where === "tissue" && !byCl && !byCt && !byPh && !lit) { const c = { key, n: 0 }; U.cloud.set(where, c); return c; }
    // counting sort on a 0..511 bucket: expression strength, +256 when lit
    const bucket = new Uint16Array(N), keep = new Uint8Array(N), start = new Uint32Array(513);
    let n = 0;
    for (let r = 0; r < N; r++) {
      if (where === "umap" && !(U.uxy[2 * r] === U.uxy[2 * r])) continue;
      let b = 0;
      if (byCl) { const c = U.cl[r]; if (U.hidden.has(c) || (c === UMAP_NONE && where !== "fly")) continue; }
      else if (byCt) { if (U.ct[r] === CT_OUT && where !== "fly") continue; }
      // phenotypes: Unassigned squares are left off the slide, so the stains show through them
      else if (byPh) { if (where === "tissue" && calls[r] === GATE_NONE && !(lit && lit(r))) continue; }
      else if (mix && mix.t[r] >= 0) b = 1 + Math.min(254, Math.round(mix.t[r] * 254));
      const l = lit ? lit(r) : false;
      if (where === "tissue" && !byCl && !byCt && !byPh && !l) continue;
      if (l) b += 256;
      keep[r] = 1; bucket[r] = b; start[b + 1]++; n++;
    }
    for (let b = 0; b < 512; b++) start[b + 1] += start[b];
    const rows = new Uint32Array(n);
    for (let r = 0; r < N; r++) if (keep[r]) rows[start[bucket[r]]++] = r;
    const src = where === "umap" ? U.uxy : U.pos, pos = new Float32Array(2 * n), col = new Uint8Array(4 * n);
    const tgt = where === "fly" ? new Float32Array(3 * n) : null;
    // over the H&E a little see-through, as the circle tool's bins are
    const base = where === "tissue" ? 185 : 225, dim = where === "tissue" ? 14 : 24;
    for (let j = 0; j < n; j++) {
      const r = rows[j], l = bucket[r] >= 256;
      pos[2 * j] = src[2 * r]; pos[2 * j + 1] = src[2 * r + 1];
      if (tgt) {
        const ux = U.uxy[2 * r], emb = ux === ux;
        tgt[3 * j] = emb ? ux : pos[2 * j]; tgt[3 * j + 1] = emb ? U.uxy[2 * r + 1] : pos[2 * j + 1]; tgt[3 * j + 2] = emb ? 1 : 0;
      }
      let c;
      if (lit && !l) c = UMAP_GREY;               // a faded cluster colour still read as the cluster in a dense UMAP
      else if (byCl) c = sectionRgb(U.cl[r]);
      else if (byCt) c = ctRgb(U.ct[r]);
      else if (byPh) c = gateRuleRgb(calls[r]);
      else if (mix && mix.t[r] >= 0) c = [mix.rgb[3 * r], mix.rgb[3 * r + 1], mix.rgb[3 * r + 2]];
      else c = l ? (where === "tissue" ? UMAP_SEL : UMAP_LIT_GREY) : UMAP_GREY;
      col[4 * j] = c[0]; col[4 * j + 1] = c[1]; col[4 * j + 2] = c[2];
      col[4 * j + 3] = lit ? (l ? 245 : dim) : base;
    }
    const attributes = { getPosition: { value: pos, size: 2 }, getFillColor: { value: col, size: 4, normalized: true } };
    if (tgt) attributes.getFlyTarget = { value: tgt, size: 3 };
    const C = { key, n, rows, data: { length: n, attributes } };
    U.cloud.set(where, C);
    return C;
  }

  function umapCloudLayer(C, where, prefix = "") {
    const U = S.umap, e = currentExpr();
    const id = `${prefix}ucloud-${where}-${U.key}`;
    U.drawn.set(id, C);
    const r = where === "umap" ? U.radius * S.pointScale
      : Math.max(0.5, (e && e.spotDiameterPx ? e.spotDiameterPx / 2 : 4) * S.pointScale);
    return new ScatterplotLayer({
      id, data: C.data,
      radiusUnits: "common", getRadius: r, radiusMinPixels: where === "umap" ? 0.8 : 0.6,
      stroked: false, filled: true,
      // hover says which cluster a dot is in; on the slide only when the dots are the clusters
      pickable: where === "umap" || umapByCat(),
      // marker squares: "Square opacity" fades the cluster colours too, so the stains show through
      opacity: where === "tissue" && isProt() ? S.geneOpacity : 1,
      parameters: NORMAL,
      updateTriggers: { getFillColor: [C.key], getPosition: [C.key] },
    });
  }

  // a label per cluster, on the dot nearest the cluster's mean (a mean alone can land in
  // the gap between two lobes of one cluster)
  function umapCentroids() {
    const U = S.umap;
    if (!U.cl || !U.uxy) return [];
    const hit = U.centroids.get(U.clRes); if (hit) return hit;
    const N = U.cl.length;
    let K = 0; for (let r = 0; r < N; r++) if (U.cl[r] !== UMAP_NONE && U.cl[r] + 1 > K) K = U.cl[r] + 1;
    const sx = new Float64Array(K), sy = new Float64Array(K), cnt = new Float64Array(K);
    for (let r = 0; r < N; r++) {
      const c = U.cl[r], x = U.uxy[2 * r]; if (c === UMAP_NONE || !(x === x)) continue;
      sx[c] += x; sy[c] += U.uxy[2 * r + 1]; cnt[c]++;
    }
    const best = new Float64Array(K).fill(Infinity), bx = new Float64Array(K), by = new Float64Array(K);
    for (let r = 0; r < N; r++) {
      const c = U.cl[r], x = U.uxy[2 * r], y = U.uxy[2 * r + 1]; if (c === UMAP_NONE || !(x === x)) continue;
      const d = (x - sx[c] / cnt[c]) ** 2 + (y - sy[c] / cnt[c]) ** 2;
      if (d < best[c]) { best[c] = d; bx[c] = x; by[c] = y; }
    }
    const out = [];
    for (let c = 0; c < K; c++) if (cnt[c]) out.push({ c, x: bx[c], y: by[c] });
    U.centroids.set(U.clRes, out);
    return out;
  }
  function umapLabelLayer(prefix = "") {
    const U = S.umap;
    const data = umapCentroids().filter((d) => !U.hidden.has(d.c));
    return new TextLayer({
      id: `${prefix}ulabels-${U.key}-${U.clRes}`, data,
      getPosition: (d) => [d.x, d.y],
      getText: (d) => ((U.names[U.clRes] || {})[d.c]) || (clUser(umapClusterInfo()) ? umapName(d.c) : String(d.c)),
      getColor: (d) => [...labelColor(sectionRgb(d.c)), U.pick === null || U.pick === d.c ? 255 : 90],
      sizeUnits: "pixels", getSize: 13, characterSet: "auto", fontFamily: "system-ui, sans-serif", fontWeight: 700,
      background: true, getBackgroundColor: (d) => [8, 11, 18, U.pick === null || U.pick === d.c ? 215 : 70], backgroundPadding: [5, 2, 5, 2],
      parameters: NORMAL, pickable: false,
      updateTriggers: { getText: [U.ver], getColor: [U.ver], getBackgroundColor: [U.ver] },
    });
  }
  // a white ring round the hovered bin: side by side it marks the same bin in both panes
  function umapRing(id, xy) {
    return new ScatterplotLayer({
      id, data: [xy], getPosition: (d) => d, radiusUnits: "pixels", getRadius: 7, filled: false, stroked: true,
      lineWidthUnits: "pixels", getLineWidth: 2, getLineColor: [255, 255, 255, 255], parameters: NORMAL,
    });
  }

  function renderUmap() {
    const U = S.umap;
    if (!U.uxy) return;
    const layers = [umapCloudLayer(umapCloud("umap"), "umap")];
    if (umapByCluster()) layers.push(umapLabelLayer());
    deckgl.setProps({ layers });
    drawDensityKey();
    updateRegionReadout();
  }

  // Side by side: the section UMAP in the right pane, plus the hover rings in both panes.
  function umapSideLayers() {
    const U = S.umap;
    if (!umapSbs() || !U.uxy || U.key !== S.exprKey) return [];
    const out = [umapCloudLayer(umapCloud("umap"), "umap", "uside-")];
    if (umapByCluster()) out.push(umapLabelLayer("uside-"));
    const r = U.hover;
    if (r >= 0 && r < U.pos.length / 2) {
      out.push(umapRing("uhover-tissue", [U.pos[2 * r], U.pos[2 * r + 1]]));
      if (U.uxy[2 * r] === U.uxy[2 * r]) out.push(umapRing("uside-hover", [U.uxy[2 * r], U.uxy[2 * r + 1]]));
    }
    return out;
  }

  // On the slide: clusters as coloured bins over the H&E; a selection or picked cluster
  // dims everything else (the slide included) and redraws the lit bins on top.
  function umapTissueLayers() {
    const U = S.umap;
    if (!U || !U.pos || inUmap() || U.key !== S.exprKey) return [];
    const lit = umapLit("tissue"), byCl = umapByCat();
    if (!byCl && !lit) return [];
    const out = [];
    if (lit) {
      const img = S.m.image;
      out.push(new PolygonLayer({
        id: "udim", data: [[[0, 0], [img.width, 0], [img.width, img.height], [0, img.height]]], getPolygon: (d) => d,
        getFillColor: [0, 0, 0, 150], stroked: false, filled: true, pickable: false, parameters: NORMAL,
      }));
    }
    const C = umapCloud("tissue");
    if (C.n) out.push(umapCloudLayer(C, "tissue"));
    return out;
  }

  function umapRowAt(layerId, i) {
    const C = S.umap && S.umap.drawn.get(layerId);
    return C && C.rows && i >= 0 && i < C.n ? C.rows[i] : -1;
  }
  function umapHover(info) {
    const U = S.umap;
    if (!U || !umapSbs() || !U.pos) return;
    let row = -1;
    if (info && info.layer && info.layer.id.includes("ucloud") && info.index >= 0) row = umapRowAt(info.layer.id, info.index);
    else if (info && info.viewport && info.viewport.id === "ortho" && info.coordinate) row = umapBinAt(info.coordinate[0], info.coordinate[1]);
    if (row !== U.hover) { U.hover = row; requestRender(); }
  }
  function umapTooltip(r) {
    if (r < 0) return null;
    const U = S.umap, parts = [], prot = isProt();
    if (U.cl) parts.push(U.cl[r] === UMAP_NONE ? (prot ? "No cluster (low signal: not in the UMAP)" : "No cluster (too few counts)") : umapName(U.cl[r]));
    if (umapByType() && U.ct[r] !== CT_OUT) parts.push(`Suggested type: ${ctName(U.ct[r])}`);
    if (gateOn() && S.gate.calls && r < S.gate.N) parts.push(`Phenotype: ${gatePhName(S.gate.calls[r])}`);
    if (!umapByCat()) {
      for (const g of S.genes) {
        if (g.visible === false) continue;
        const v = geneByRow(g)[r];
        parts.push(prot ? `${g.name}: ${protVal(v)} mean intensity` : `${g.name}: ${v ? (Number.isInteger(v) ? v : v.toFixed(2)) : 0} UMI`);
      }
    }
    if (prot && parts.length) {
      const e = currentExpr();
      parts.push(`(one ${e && e.binSizeUm ? `${e.binSizeUm} µm ` : ""}square)`);
    }
    return parts.length ? { text: parts.join("\n") } : null;
  }

  function umapPick(c) {
    const U = S.umap;
    U.pick = c === null || c === undefined || U.pick === c ? null : c;
    U.sel = null; U.selN = 0; U.ctPick = null; U.phPick = null;
    umapDirty(); syncUmapUI(); render();
  }
  function umapClearSel() {
    const U = S.umap;
    if (U) { U.sel = null; U.selN = 0; U.pick = null; U.ctPick = null; U.phPick = null; }
    lassoClearPick();
    if (U) { umapDirty(); syncUmapUI(); }
    render();
  }
  // a brush, picked shape or picked cluster of Circle -> UMAP
  function lassoClearPick() {
    const L = S.lasso;
    if (!L || (!L.brush && L.selShape === null && L.selCluster === null)) return false;
    L.brush = null; L.selShape = null; L.selCluster = null; lassoRefresh(); buildLassoUI();
    return true;
  }

  // ---- "Show on slide": from a UMAP pick back to the tissue ------------------------------------
  // Frames the lit squares on the slide. From the UMAP view the dots fly back to the tissue,
  // ending on them; side by side and on the tissue the tissue pane moves. Tiny picks are
  // framed with some room (at least 30 squares across) so their neighbourhood shows too.
  function slideBoxVs(x0, y0, x1, y1) {
    const e = currentExpr(), sq = (e && e.spotDiameterPx) || 8, minW = sq * 30;
    if (x1 - x0 < minW) { const m = (x0 + x1) / 2; x0 = m - minW / 2; x1 = m + minW / 2; }
    if (y1 - y0 < minW) { const m = (y0 + y1) / 2; y0 = m - minW / 2; y1 = m + minW / 2; }
    const [W, H] = umapPaneSize("ortho"), lp = $("lassoPanel");
    // the Circle -> UMAP panel covers the right of the slide (on a phone its lower part):
    // frame the squares in what is left
    const open = !lp.hidden && !umapSbs(), below = open && lp.offsetWidth > W * 0.8;
    const free = Math.max(200, W - (open && !below ? lp.offsetWidth + 20 : 0)), freeH = Math.max(150, H - (below ? lp.offsetHeight + 50 : 0));
    const zoom = Math.min(2, Math.log2(Math.min(free / ((x1 - x0) * 1.2 || 1), freeH / ((y1 - y0) * 1.2 || 1))));
    const scale = Math.pow(2, zoom), base = (S.umap && S.umap.view === "umap" ? S.umap.tissueVs : S.viewState) || imageFitVs();
    return { ...base, target: [(x0 + x1) / 2 + (W - free) / 2 / scale, (y0 + y1) / 2 + (H - freeH) / 2 / scale, 0], zoom };
  }
  // test(row) -> lit?; pos = the rows' slide positions. Returns how many were framed.
  function showRowsOnSlide(test, pos) {
    if (!test || !pos) return 0;
    let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity, n = 0;
    for (let r = 0; r < pos.length / 2; r++) {
      if (!test(r)) continue;
      const x = pos[2 * r], y = pos[2 * r + 1];
      if (!(x === x)) continue;
      n++; if (x < x0) x0 = x; if (x > x1) x1 = x; if (y < y0) y0 = y; if (y > y1) y1 = y;
    }
    if (!n) return 0;
    const vs = slideBoxVs(x0, y0, x1, y1), U = S.umap;
    autoFit = false;
    if (U && U.view === "umap") { U.tissueVs = vs; umapSetView("tissue"); return n; }
    S.viewState = vs;
    deckgl.setProps({ viewState: deckVs() }); updateHud(); requestRender();
    return n;
  }
  // any set of rows (mask: 1 = lit), e.g. a group of the heatmap: lit as a selection where the layer has a UMAP, framed either way
  async function umapLightRows(mask, n) {
    const U = S.umap;
    if (!umapAvail()) return showRowsOnSlide((r) => mask[r] === 1, S.lasso && S.lasso.pos);
    await umapLoad();
    if (S.umap !== U || !U.pos || mask.length * 2 !== U.pos.length) return 0;
    let emb = 0; for (let r = 0; r < mask.length; r++) if (mask[r] && U.uxy[2 * r] === U.uxy[2 * r]) emb++;
    U.sel = mask; U.selN = n; U.selEmb = emb; U.pick = null; U.ctPick = null; U.phPick = null;
    umapDirty(); syncUmapUI(); render();
    return showLitOnSlide();
  }
  // whatever is lit: the section selection / cluster / type / phenotype, else the circle tool's pick
  function showLitOnSlide() {
    const U = S.umap, L = S.lasso;
    const f = U && U.pos ? umapLit("tissue") : null;
    if (f) return showRowsOnSlide(f, U.pos);
    const g = L && L.assign && lassoLitRows();
    if (g) return showRowsOnSlide((r) => L.assign[r] !== LASSO_NONE && g(r, L.assign[r]), L.pos);
    return 0;
  }

  // ---- box / freehand selection, in any view ---------------------------------------------
  function umapSetTool(t) {
    const U = S.umap; if (!U) return;
    U.tool = t && U.tool !== t ? t : null;
    if (U.tool && S.lasso && S.lasso.tool) lassoSetTool(null);
    $("selBox").classList.toggle("on", U.tool === "box");
    $("selFree").classList.toggle("on", U.tool === "free");
    $("selPad").hidden = !U.tool; $("selSvg").innerHTML = "";
  }
  function umapSelEvents() {
    const pad = $("selPad"), svg = $("selSvg");
    let D = null;
    const local = (ev) => { const r = pad.getBoundingClientRect(); return [ev.clientX - r.left, ev.clientY - r.top]; };
    const boxPts = (a, b) => [[a[0], a[1]], [b[0], a[1]], [b[0], b[1]], [a[0], b[1]]];
    const outline = () => (D.kind === "box" ? (D.pts.length > 1 ? boxPts(D.pts[0], D.pts[1]) : []) : D.pts);
    const draw = () => { svg.innerHTML = D ? `<polygon points="${outline().map((p) => p.join(",")).join(" ")}" />` : ""; };
    padPinch(pad, () => { D = null; draw(); });
    pad.addEventListener("pointerdown", (ev) => {
      const U = S.umap; if (!U.tool || ev.button !== 0) return;
      pad.setPointerCapture(ev.pointerId);
      const p = local(ev);
      // side by side: the pane the drag starts in decides which layout is being selected
      const vp = deckgl.getViewports().find((v) => p[0] >= v.x && p[0] < v.x + v.width && p[1] >= v.y && p[1] < v.y + v.height)
        || deckgl.getViewports()[0];
      // no Shift on a phone: the "Add to the selection" box does the same
      D = { kind: U.tool, pts: [p], add: ev.shiftKey || $("selAdd").checked, vp }; draw();
    });
    pad.addEventListener("pointermove", (ev) => {
      if (!D) return;
      const p = local(ev);
      if (D.kind === "box") D.pts[1] = p;
      else { const q = D.pts[D.pts.length - 1]; if (Math.hypot(p[0] - q[0], p[1] - q[1]) >= 3) D.pts.push(p); }
      draw();
    });
    const finish = () => {
      if (!D) return;
      const pts = outline(), add = D.add, vp = D.vp; D = null; svg.innerHTML = "";
      if (pts.length < 3) return;
      let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity;
      for (const [x, y] of pts) { x0 = Math.min(x0, x); y0 = Math.min(y0, y); x1 = Math.max(x1, x); y1 = Math.max(y1, y); }
      if (x1 - x0 < 4 && y1 - y0 < 4) return;                    // a click, not a selection
      const space = vp.id === "uside" || inUmap() ? "umap" : "tissue";
      umapSelect(pts.map((p) => vp.unproject([p[0] - vp.x, p[1] - vp.y]).slice(0, 2)), add, space);
      umapSetTool(null);                                          // one selection per click of the button
    };
    pad.addEventListener("pointerup", finish);
    pad.addEventListener("pointercancel", finish);
  }
  async function umapStartTool(t) {
    const U = S.umap;
    if (!U.pos) {
      try { U.pos = await umapPositions(); } catch (e) { status(`Could not load ${wd().bin} positions: ${e.message}`, true); return; }
    }
    umapSetTool(t);
  }
  function umapSelect(poly, add, space) {
    const U = S.umap, src = space === "umap" ? U.uxy : U.pos;
    if (!src) return;
    const N = src.length / 2, sel = add && U.sel ? U.sel : new Uint8Array(N);
    let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity;
    for (const [x, y] of poly) { x0 = Math.min(x0, x); y0 = Math.min(y0, y); x1 = Math.max(x1, x); y1 = Math.max(y1, y); }
    const hid = umapByCluster() && U.hidden.size ? U.hidden : null;   // hidden clusters cannot be picked up
    for (let r = 0; r < N; r++) {
      const x = src[2 * r], y = src[2 * r + 1];
      if (!(x >= x0 && x <= x1 && y >= y0 && y <= y1)) continue;    // NaN fails this too
      if (hid && hid.has(U.cl[r])) continue;
      if (inPoly(x, y, poly)) sel[r] = 1;
    }
    let n = 0, emb = 0;
    for (let r = 0; r < N; r++) if (sel[r]) { n++; if (U.uxy && U.uxy[2 * r] === U.uxy[2 * r]) emb++; }
    U.sel = n ? sel : null; U.selN = n; U.selEmb = emb; U.pick = null; U.ctPick = null; U.phPick = null;
    umapDirty(); syncUmapUI(); render();
  }
  // The selection becomes one more Circle -> UMAP shape: a list of bins rather than an outline.
  function umapSelToShape() {
    const U = S.umap, L = S.lasso;
    if (!U.sel || !L || !L.pos || L.pos.length !== U.sel.length * 2) return;
    const rows = [];
    for (let r = 0; r < U.sel.length; r++) if (U.sel[r]) rows.push(r);
    const id = L.next++, name = `Selection ${id}`;
    L.shapes.push({ id, name, kind: "rows", pts: [], rows: Uint32Array.from(rows),
                    hex: LASSO_COLORS[(id - 1) % LASSO_COLORS.length], bins: 0, shared: 0, lost: 0 });
    // the shape carries these bins from here on; keeping the selection lit as well would
    // stack two highlights on the same bins
    U.sel = null; U.selN = 0; umapDirty(); syncUmapUI();
    lassoShapesChanged();
    status(`Added "${name}" (${rows.length.toLocaleString()} bins) to Draw a region`);
    setTimeout(() => status(null), 3500);
  }

  // ---- legend: one row per cluster ----------------------------------------------------------
  // Colour, name (✎ renames it: "Schwann cells"), size and the top marker genes. Click a row
  // to highlight that cluster; the box hides it. Shared by the main view and the circle panel.
  const markerName = (m) => {
    if (m == null) return "";
    if (typeof m === "string") return m;
    const g = m.gene;
    const names = S.lasso && S.lasso.hvgGenes;
    return typeof g === "number" && names ? names[g] || String(g) : String(g);
  };
  // With a marker set (o.suggest), each row also carries a muted line: the suggested cell
  // type, a confidence badge and "✓ Accept" (which names the cluster). A name given with ✎
  // always wins; the suggestion stays visible under it.
  const escHtml = (x) => String(x).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
  const CT_BADGE = { high: "high", medium: "medium", low: "low", none: "unclear" };
  const pctTxt = (p) => `${Math.round((p || 0) * 100)}%`;
  function clusterLegend(host, o) {
    const esc = escHtml, W = wd(), prot = isProt();
    host.innerHTML = "";
    if (o.suggest) {
      const nHigh = o.acceptAll ? o.acceptAll(true) : 0;
      const head = document.createElement("div");
      head.className = "ct-head";
      head.innerHTML = `<span>Suggested names from <button class="linkish ct-info" title="The ${W.genes} behind each suggestion, and how it is made">marker lists</button></span>` +
        (o.acceptAll ? `<button class="ct-all" ${nHigh ? "" : "disabled"} title="Name every high-confidence cluster after its suggestion (clusters you already named keep their name)">Accept all high${nHigh ? ` (${nHigh})` : ""}</button>` : "");
      head.querySelector(".ct-info").onclick = (ev) => { ev.stopPropagation(); ctPopover(ev.currentTarget, o.markerInfo()); };
      const all = head.querySelector(".ct-all");
      if (all) all.onclick = (ev) => { ev.stopPropagation(); o.acceptAll(false); };
      host.appendChild(head);
    }
    for (let c = 0; c < o.n; c++) {
      const off = !!(o.hidden && o.hidden.has(c));
      const mk = o.markers[c] || [];
      const top = mk.slice(0, 3).map(markerName).filter(Boolean);
      const [r, g, b] = (o.rgb || clusterRgb)(c);
      const row = document.createElement("div");
      row.className = "cl-row" + (o.picked === c ? " sel" : "") + (off ? " off" : "");
      row.innerHTML = (o.onHide ? `<input type="checkbox" class="cb" ${off ? "" : "checked"} title="Show or hide this cluster" />` : `<span></span>`) +
        `<i style="background:rgb(${r},${g},${b})"></i>` +
        `<div class="cl-txt"><span class="nm"></span><small>${(o.sizes[c] || 0).toLocaleString()} ${W.bins}` +
        `${top.length ? ` · ${esc(top.join(", "))}` : ""}</small></div>` +
        `<button class="ren" title="Rename this cluster">✎</button>`;
      const nm = row.querySelector(".nm");
      nm.textContent = o.name(c);                                 // typed by the user: never as HTML
      row.title = mk.length
        ? (o.markerTitle ? o.markerTitle
           : prot ? `Top markers (how much higher than in the rest of the section, in standard deviations; ` +
                  `share of the cluster's squares above that marker's typical level):\n`
                : `Top marker genes (log fold change vs the rest, share of bins with the gene):\n`) +
          mk.slice(0, 10).map((m) => `${markerName(m)}  ${m.lfc != null ? (+m.lfc).toFixed(2) : ""}  ${m.pct != null ? Math.round(m.pct * 100) + "%" : ""}`).join("\n") +
          `\nClick to highlight this cluster in every view.`
        : "Click to highlight this cluster in every view.";
      const sg = o.suggest && o.suggest(c);
      if (sg) {
        const line = document.createElement("div");
        line.className = "cl-sug";
        const conf = sg.confidence in CT_BADGE ? sg.confidence : "none";
        const named = o.name(c) === sg.suggest;
        const mix = conf !== "none" && sg.second && sg.pct2 >= 0.1 ? ` · ${sg.second} ${pctTxt(sg.pct2)}` : "";
        line.innerHTML = (sg.artefactRule ? `<b class="ct-badge low qc-warn" title="Most of this cluster's squares are flagged everything bright">check</b>`
                                          : `<b class="ct-badge ${conf}">${CT_BADGE[conf]}</b>`) +
          `<span class="txt">${sg.artefactRule ? `${esc(sg.suggest)} · ${pctTxt(sg.artefact)} of ${W.bins} flagged everything bright`
                               : conf === "none" ? `Unclear · no list fits most ${W.bins} (best: ${esc(sg.top)} ${pctTxt(sg.pct)})`
                                                 : `${esc(sg.suggest)} · ${pctTxt(sg.pct)} of ${W.bins}${esc(mix)}`}</span>` +
          (conf === "none" ? "" : named ? `<span class="acc-ok" title="This cluster carries the suggested name">✓</span>`
                                        : `<button class="acc" title="Name this cluster “${esc(sg.suggest)}”">✓ Accept</button>`);
        line.title = `Suggested: ${sg.suggest} (${CT_BADGE[conf]} confidence)\n` +
          `${pctTxt(sg.pct)} of this cluster's ${W.bins} score best for ${sg.top}` +
          (sg.second ? `, ${pctTxt(sg.pct2)} for ${sg.second}` : "") + `.\nClick “marker lists” above for the ${W.genes} and the rule.`;
        const acc = line.querySelector(".acc");
        if (acc) acc.onclick = (ev) => { ev.stopPropagation(); o.onAccept(c, sg.suggest); };
        row.appendChild(line);                                    // under the name, the full width of the row
      }
      const cb = row.querySelector(".cb");
      if (cb) { cb.onclick = (e) => e.stopPropagation(); cb.onchange = () => o.onHide(c, !cb.checked); }
      row.querySelector(".ren").onclick = (ev) => {
        ev.stopPropagation();
        const inp = document.createElement("input");
        inp.type = "text"; inp.value = o.name(c); inp.className = "nm-edit";
        let fin = false;                                          // Escape, then the blur it causes, must not commit
        const done = (keep) => {
          if (fin) return; fin = true;
          o.onRename(c, keep && inp.value.trim() ? inp.value.trim().slice(0, 60) : null);
        };
        inp.onkeydown = (e) => { e.stopPropagation(); if (e.key === "Enter") done(true); if (e.key === "Escape") done(false); };
        inp.onblur = () => done(true);
        inp.onclick = (e) => e.stopPropagation();
        nm.replaceWith(inp); inp.focus(); inp.select();
      };
      row.onclick = () => o.onPick(c);
      host.appendChild(row);
    }
  }

  // ---- "Marker lists": which genes each suggestion rests on, and the rule --------------
  // One floating box, opened from either legend; Esc, × or a click elsewhere closes it.
  const CT_RULE = "Each bin gets a score for each cell type: how far above their usual level that type's genes are in the bin " +
    "(log-normalised counts, as z-scores), averaged over the type's genes. A bin counts as the type with the best score " +
    "if that score is at least 0.25, otherwise as none. A cluster is named after the type most of its bins count as: " +
    "high = at least half of its bins and twice the next type; medium = at least 30%; low (shown as “maybe …”) = at least 15%; " +
    "anything less is Unclear. Genes found in under 0.5% of the bins are left out, and a type needs at least 2 genes left to be scored.";
  // the same rule on marker squares (celltypes.json "modality": "protein"): intensities, not counts
  const CT_BIN_MIN_PROT = 0.5;
  const CT_RULE_PROT = "Each square gets a score for each cell type: how far above their usual level that type's markers are in the square " +
    "(mean intensity with the slide's background taken off, divided by the marker's brightest 0.5% of tissue squares and put through asinh, " +
    "as z-scores), averaged over the type's markers. " +
    `A square counts as the type with the best score if that score is at least ${CT_BIN_MIN_PROT}, otherwise as none. ` +
    "A cluster is named after the type most of its squares count as: high = at least half of its squares and twice the next type; " +
    "medium = at least 30%; low (shown as “maybe …”) = at least 15%; anything less is Unclear. " +
    "A type is scored when enough of its markers are on this panel (the number is next to each type below); names match ignoring " +
    "case, number prefixes, dye tags and the list's aliases; a marker that is nearly the same in every square is left out. " +
    "A cluster with more than half of its squares flagged “everything bright” (folds, bubbles, edges) is suggested as “Artefact? (fold / edge)”. " +
    "Squares hold several cells, so the names describe neighbourhoods. (Exports made before November 2026 scale each marker from its " +
    "dimmest 1% to its brightest 0.5% of tissue squares instead.)";
  const ctProt = (ci) => !!((ci && ci.modality === "protein") || isProt());
  function ctPopover(anchor, info) {
    let el = $("ctPop");
    if (!el) {
      el = document.createElement("div");
      el.id = "ctPop"; el.className = "ct-pop"; el.hidden = true;
      document.body.appendChild(el);
      document.addEventListener("pointerdown", (e) => { if (!el.hidden && !el.contains(e.target)) el.hidden = true; }, true);
      window.addEventListener("keydown", (e) => { if (e.key === "Escape") el.hidden = true; });
    }
    if (!info) return;
    const esc = escHtml, set = info.markerSet || {}, prot = !!info.prot;
    el.innerHTML =
      `<div class="ct-pop-h"><b>Marker lists${set.name ? `: ${esc(set.name)}` : ""}${set.version != null ? ` <small>v${esc(set.version)}</small>` : ""}</b>` +
      `<button class="ct-x" title="Close">×</button></div>` +
      (info.scope ? `<p class="meta">${esc(info.scope)}</p>` : "") +
      `<p>${esc(prot ? CT_RULE_PROT : CT_RULE)}</p>` +
      (info.note ? `<p class="meta">${esc(info.note)}</p>` : "") +
      `<div class="ct-types">` + (info.types || []).map((t) => {
        // missing genes: "absent" = no such gene in the data used, "rare" = in under 0.5% of the bins
        const miss = (t.missing || []).map((m) => (typeof m === "string" ? { gene: m, reason: "absent" } : m));
        // protein: "flat" = a channel with no contrast over the tissue, left out like a rare gene
        const flat = prot ? miss.filter((m) => m.reason === "flat").map((m) => m.gene) : [];
        const absent = miss.filter((m) => m.reason !== "rare" && !(prot && m.reason === "flat")).map((m) => m.gene), rare = miss.filter((m) => m.reason === "rare").map((m) => m.gene);
        const need = t.min != null ? t.min : 1;
        const why = prot ? (flat.length ? `too few of its markers are on this panel or vary here (needs ${need})`
                                        : need > 1 ? `fewer than ${need} of its markers are on this panel` : "none of its markers are on this panel")
          : rare.length && !absent.length ? "too rare in this data" : absent.length && !rare.length ? "its genes are not in the data" : "too few of its genes are measurable here";
        return `<div class="ct-type${t.scored ? "" : " off"}"><b>${esc(t.name)}</b>${prot && t.min != null ? ` <small class="meta" title="Markers this type needs on the panel">needs ${t.min}</small>` : ""}${t.scored ? "" : ` <small title="${esc(t.reason || "")}">not scored: ${why}</small>`}` +
          `<div>${t.genes.length ? `<span class="found" title="${prot ? "Markers" : "Genes"} used for the score">${esc(t.genes.join(", "))}</span>` : `<span class="meta">no usable ${prot ? "markers" : "genes"}</span>`}` +
          (absent.length && prot ? ` <span class="miss" title="No channel of this name on this panel">· not on this panel: ${esc(absent.join(", "))}</span>` : "") +
          (flat.length ? ` <span class="miss" title="Nearly the same in every tissue square: left out, as it cannot tell squares apart">· flat here: ${esc(flat.join(", "))}</span>` : "") +
          (absent.length && !prot ? ` <span class="miss" title="No gene of this name in the data used here">· not in the data: ${esc(absent.join(", "))}</span>` : "") +
          (rare.length ? ` <span class="miss" title="Counted in under 0.5% of the bins: left out, as one stray count would decide a bin">· too rare: ${esc(rare.join(", "))}</span>` : "") +
          `</div></div>`;
      }).join("") +
      `</div>`;
    el.querySelector(".ct-x").onclick = () => { el.hidden = true; };
    el.hidden = false;
    const r = anchor.getBoundingClientRect(), w = Math.min(420, window.innerWidth - 24);
    el.style.width = `${w}px`;
    el.style.left = `${Math.max(12, Math.min(r.left, window.innerWidth - w - 12))}px`;
    const h = el.offsetHeight, below = r.bottom + 6;
    el.style.top = `${below + h < window.innerHeight - 8 ? below : Math.max(8, r.top - h - 6)}px`;
  }

  // ---- whole section: suggested cell types (celltypes.json, celltype.bin) ---------------
  const CT_NONE = 254, CT_OUT = 255;                // celltype.bin: embedded but no type / not embedded
  const CT_COLORS = ["#3fa7ff", "#9fd3ff", "#ff9f1c", "#5cc95c", "#ff5c5c", "#ff9ecf", "#b388ff", "#ffd84d",
                     "#c8a27a", "#2ec4b6", "#fff3b0", "#e056fd"];
  const ctLut = [];
  const ctRgb = (t) => (t >= CT_NONE ? [120, 120, 120]
    : ctLut[t] || (ctLut[t] = liftForDark(t < CT_COLORS.length ? hexToRgb255(CT_COLORS[t]) : hslRgb((t * 137.508 + 40) % 360, 0.6, 0.6))));
  const umapCtAvail = () => { const m = umapMeta(); return !!(m && m.celltypes); };
  const umapByType = () => !!(S.umap && S.umap.colorBy === "celltype" && S.umap.ct);
  const umapByCat = () => umapByCluster() || umapByType() || umapByPh();      // one colour per bin replaces the genes
  const umapSuggest = (c) => { const U = S.umap, a = U.ctInfo && U.ctInfo.byRes.get(U.clRes); return (a && a[c]) || null; };

  // the suggestions per clustering: loaded with the clusters, quietly absent for older exports
  async function umapLoadCtInfo() {
    const U = S.umap;
    if (U.ctInfo || U.ctTried || !umapCtAvail() || !U.key) return;
    U.ctTried = true;
    try {
      const j = await S.src.json(`expr/${U.key}/umap/celltypes.json`);
      if (S.umap !== U || !j || !Array.isArray(j.types)) return;
      const byRes = new Map();
      for (const [k, list] of Object.entries(j.clusterings || {})) {
        const a = []; for (const s of list || []) a[s.cluster] = s;
        // "1" in the file is resolution 1.0 here; the lab's own clusterings are keyed by their id ("user1")
        byRes.set(/^-?[\d.]+$/.test(k) ? Number(k) : k, a);
      }
      U.ctInfo = { ...j, byRes };
      // only for its note in "marker lists"; a protein set names its own file (markers/akoya_*.json)
      if (ctProt(j) && j.markerSet && j.markerSet.file) loadMarkerNote(j.markerSet.file);
      else loadMarkerSet();
    } catch (e) { U.ctInfo = null; }
  }
  // each bin's suggested type, for Colour bins by "Suggested cell type"
  async function umapLoadTypes() {
    const U = S.umap;
    await umapLoad();
    if (S.umap !== U) return;
    await umapLoadCtInfo();
    if (S.umap !== U) return;
    if (!U.ctInfo) throw new Error("this layer has no suggested cell types");
    if (!U.ct) {
      status("Loading cell types…");
      let arr;
      try { arr = new Uint8Array(await S.src.buf(`expr/${U.key}/umap/celltype.bin`)); }
      finally { status(null); }
      if (S.umap !== U) return;
      if (arr.length !== U.pos.length / 2) throw new Error(`celltype.bin has ${arr.length} rows, positions.bin ${U.pos.length / 2}`);
      const counts = new Array(256).fill(0);
      for (let r = 0; r < arr.length; r++) counts[arr[r]]++;
      U.ct = arr; U.ctCounts = counts;
    }
    umapDirty();
  }
  // A protein marker set (its own file, not the circle tool's RNA set): fetched once, for its note.
  const markerNotes = new Map();
  function loadMarkerNote(file) {
    if (!/^markers\/[\w.-]+\.json$/.test(file) || markerNotes.has(file)) return;
    markerNotes.set(file, null);
    fetch(file, { cache: "no-cache" }).then((r) => (r.ok ? r.json() : null))
      .then((j) => { if (j && j.name) markerNotes.set(file, j); }).catch(() => {});
  }
  const ctMarkerInfo = (ci, scope) => {
    if (!ci) return ci;
    const own = ci.markerSet && ci.markerSet.file && markerNotes.get(ci.markerSet.file);
    const ms = own && own.name === ci.markerSet.name ? own
      : S.markerSet && ci.markerSet && S.markerSet.name === ci.markerSet.name ? S.markerSet : null;
    return { markerSet: ci.markerSet, types: ci.types, scope, prot: ctProt(ci), note: ms ? ms.note : null };
  };
  const ctScope = () => (isProt() ? "Scored on every marker of the panel, over all squares in the section UMAP."
    : "Scored on every gene, over all bins in the section UMAP.");
  function umapPickType(t) {
    const U = S.umap;
    U.ctPick = t === null || t === undefined || U.ctPick === t ? null : t;
    U.pick = null; U.sel = null; U.selN = 0; U.phPick = null;
    umapDirty(); syncUmapUI(); render();
  }
  const umapCtLit = () => !!(S.umap && S.umap.ctPick !== null && umapByType());
  const ctName = (t) => (t === CT_NONE ? "No type" : (S.umap.ctInfo && S.umap.ctInfo.types[t] && S.umap.ctInfo.types[t].name) || `Type ${t}`);

  // Colour bins by "Suggested cell type": one row per type in the marker set
  function typeLegend(host) {
    const U = S.umap, ci = U.ctInfo, esc = escHtml, W = wd(), prot = ctProt(ci);
    host.innerHTML = "";
    const head = document.createElement("div");
    head.className = "ct-head";
    head.innerHTML = `<span>Each ${W.bin}'s best-scoring type from the <button class="linkish ct-info" title="The ${W.genes} behind each type, and how a ${W.bin} is called">marker lists</button></span>`;
    head.querySelector(".ct-info").onclick = (ev) => { ev.stopPropagation(); ctPopover(ev.currentTarget, ctMarkerInfo(ci, ctScope())); };
    host.appendChild(head);
    const rows = ci.types.map((t, i) => ({ t: i, name: t.name, scored: t.scored, n: U.ctCounts[i] }));
    rows.sort((a, b) => b.scored - a.scored || b.n - a.n || a.t - b.t);
    rows.push({ t: CT_NONE, name: "No type", scored: true, n: U.ctCounts[CT_NONE], none: true });
    for (const d of rows) {
      const [r, g, b] = ctRgb(d.t);
      const row = document.createElement("div");
      row.className = "cl-row" + (U.ctPick === d.t ? " sel" : "") + (d.scored ? "" : " off");
      row.innerHTML = `<span></span><i style="background:${d.scored ? `rgb(${r},${g},${b})` : "transparent"}"></i>` +
        `<div class="cl-txt"><span class="nm">${esc(d.name)}</span><small>${!d.scored ? (prot ? "not scored: its markers are missing from this panel or flat here" : "not scored: too few of its genes in the data")
          : d.none ? `${d.n.toLocaleString()} ${W.bins} · no type scores ${prot ? CT_BIN_MIN_PROT : "0.25"} or more` : `${d.n.toLocaleString()} ${W.bins}`}</small></div><span></span>`;
      if (d.scored) { row.title = `Click to highlight these ${W.bins} in every view.`; row.onclick = () => umapPickType(d.t); }
      host.appendChild(row);
    }
  }

  // Colour squares by "Phenotype": one row per rule, in rule order, then Unassigned
  function phenoLegend(host) {
    const U = S.umap, G = S.gate, W = wd(), esc = escHtml;
    host.innerHTML = "";
    const head = document.createElement("div");
    head.className = "ct-head";
    head.innerHTML = `<span>Your cell types: the first matching rule per ${W.bin} (Combine markers…)</span>`;
    host.appendChild(head);
    const ids = G.rules.map((_, i) => i).slice(0, GATE_NONE).concat([GATE_NONE]);
    for (const k of ids) {
      const bad = k !== GATE_NONE && gateRuleProblems(G.rules[k]).length > 0;
      const [r, g, b] = gateRuleRgb(k), n = G.counts ? G.counts[k] : 0;
      const row = document.createElement("div");
      row.className = "cl-row" + (U.phPick === k ? " sel" : "") + (bad ? " off" : "");
      row.innerHTML = `<span></span><i style="background:rgb(${r},${g},${b})"></i>` +
        `<div class="cl-txt"><span class="nm">${esc(gatePhName(k))}</span><small>${bad ? "matches nothing yet: see its rule"
          : `${n.toLocaleString()} ${W.bins} · ${G.N ? ((100 * n) / G.N).toFixed(1) : 0}%`}</small></div><span></span>`;
      if (!bad) { row.title = `Click to highlight these ${W.bins} in every view.`; row.onclick = () => umapPickPh(k); }
      host.appendChild(row);
    }
  }

  function buildUmapLegend() {
    const U = S.umap, host = $("clusterLegend"), info = umapClusterInfo();
    if (U.colorBy === "celltype" && U.ct && U.ctInfo) { typeLegend(host); return; }
    if (U.colorBy === "phenotype" && umapByPh()) { phenoLegend(host); return; }
    if (U.colorBy !== "cluster" || !U.cl || !info) { host.innerHTML = ""; return; }
    let sizes = info.sizes;
    if (!sizes) { sizes = []; for (const c of U.cl) if (c !== UMAP_NONE) sizes[c] = (sizes[c] || 0) + 1; }
    const n = info.n || sizes.length;
    const setName = (c, t) => { (U.names[U.clRes] = U.names[U.clRes] || {})[c] = t; };
    const hasSug = !!(U.ctInfo && U.ctInfo.byRes.get(U.clRes));
    clusterLegend(host, {
      n, sizes, markers: info.markers || [], picked: U.pick, hidden: U.hidden, name: umapName, rgb: sectionRgb,
      onPick: (c) => umapPick(c),
      onRename: (c, t) => {
        if (t) setName(c, t);
        umapDirty(); syncUmapUI(); render();
      },
      suggest: hasSug ? umapSuggest : null,
      markerInfo: () => ctMarkerInfo(U.ctInfo, ctScope()),
      onAccept: (c, t) => { setName(c, t); umapDirty(); syncUmapUI(); render(); },
      // count only = how many high clusters still have no name of their own
      acceptAll: (count) => {
        const named = U.names[U.clRes] || {};
        const todo = [];
        for (let c = 0; c < n; c++) { const s = umapSuggest(c); if (s && s.confidence === "high" && !named[c]) todo.push([c, s.suggest]); }
        if (count) return todo.length;
        todo.forEach(([c, t]) => setName(c, t));
        umapDirty(); syncUmapUI(); render();
        return todo.length;
      },
      onHide: (c, hide) => {
        if (hide) U.hidden.add(c); else U.hidden.delete(c);
        if (hide && U.pick === c) U.pick = null;
        umapDirty(); syncUmapUI(); render();
      },
    });
  }

  function syncUmapUI() {
    const U = S.umap, ok = !!U && umapAvail();
    $("viewSwitch").hidden = !ok; $("panelUmap").hidden = !ok;
    const um = inUmap();
    const why = "Shapes are drawn on the tissue: switch the view to Tissue or Side by side to draw one.";
    $("lassoCircle").disabled = um; $("lassoFree").disabled = um;
    $("lassoCircle").title = um ? why : "Drag from the centre outwards (Esc to stop)";
    $("lassoFree").title = um ? why : `Hold and draw around the ${wd().bins} (Esc to stop)`;
    $("sbsCaps").hidden = !ok || !umapSbs();
    // groupings without an embedding (GeoMx ROIs): the Tissue / UMAP switch has nothing to switch to
    const noEmb = ok && umapMeta().embedding === "none";
    for (const el of $("viewSwitch").querySelectorAll(".seg, .fly-play, .fly-speed")) el.hidden = noEmb;
    if (!ok) { syncGateColourOpt(); gateCountsRefresh(); return; }
    $("viewSwitch").querySelectorAll("button[data-v]").forEach((b) => b.classList.toggle("on", b.dataset.v === U.view));
    $("flyPlay").disabled = !!U.anim;
    const list = umapClusterings();
    $("colorBy").value = U.colorBy;
    $("colorBy").options[1].disabled = !list.length;
    $("colorBy").options[1].textContent = noEmb ? "Groups (annotation)" : list.some(clUser) ? "Clusters" : "Clusters (Leiden)";
    const ctOpt = $("colorBy").options[2], ctOk = umapCtAvail();      // older exports: no such option
    ctOpt.hidden = !ctOk; ctOpt.disabled = !ctOk;
    syncGateColourOpt();
    const resSel = $("clusterRes");
    const opts = list.map((c) => `<option value="${escHtml(clKey(c))}">${escHtml(clLabel(c))} · ${c.n} ${clUser(c) ? "groups" : "clusters"}</option>`).join("");
    if (resSel.dataset.opts !== opts) { resSel.innerHTML = opts; resSel.dataset.opts = opts; }
    if (U.res != null) resSel.value = String(U.res);
    $("rowRes").hidden = U.colorBy !== "cluster";
    $("clusterTools").hidden = U.colorBy !== "cluster" || !U.cl;
    $("clShowAll").disabled = !U.hidden.size;
    buildUmapLegend();
    const m = umapMeta(), e = currentExpr(), W = wd(), prot = isProt();
    const nEmb = U.nEmb || m.nEmbedded;
    $("sbsCapUmap").textContent = `UMAP · whole section${nEmb ? ` · ${nEmb.toLocaleString()} ${W.bins}` : ""}`;
    $("umapNote").textContent = noEmb && window.SodaRoi ? window.SodaRoi.umapNote(e, list.length) : prot
      ? (nEmb ? (nEmb >= e.nBins ? `All ${e.nBins.toLocaleString()} tissue squares` : `${nEmb.toLocaleString()} of the ${e.nBins.toLocaleString()} tissue squares`) +
                ` (${e.binSizeUm ? `${e.binSizeUm} µm` : "each"} across) are in the UMAP, placed by how alike their markers are. ` : "") +
        `Switch between Tissue, Side by side and UMAP in the bar above the picture; markers, clusters and selections stay.`
      : m.embedding === "user" ? (nEmb ? `${nEmb.toLocaleString()} of the ${e.nBins.toLocaleString()} ${W.bins} are in the UMAP, which came with the data (the lab's own analysis). ` : "") +
        `Switch between Tissue, Side by side and UMAP in the bar above the picture; genes, clusters and selections stay.`
      : (nEmb ? `${nEmb.toLocaleString()} of the ${e.nBins.toLocaleString()} ${W.bins} (those with at least 10 ${countUnit()}${m.minHvgDetected ? ` and ${m.minHvgDetected} variable genes` : ""}) are in the UMAP. ` : "") +
        `Switch between Tissue, Side by side and UMAP in the bar above the picture; genes, clusters and selections stay.`;
    const ctLit = umapCtLit(), phLit = umapPhLit(), has = !!(U.sel || U.pick !== null || ctLit || phLit);
    const lpTxt = !has ? lassoPickTxt() : "";
    $("selClear").disabled = !has && !lpTxt;
    $("selShow").hidden = !has && !lpTxt;
    $("selToShape").hidden = !U.sel || $("panelLasso").hidden;
    const selTxt = U.sel ? `${U.selN.toLocaleString()} ${U.selN === 1 ? W.bin : W.bins} selected` : U.pick !== null ? `${umapName(U.pick)} highlighted`
      : ctLit ? `${ctName(U.ctPick)} highlighted` : phLit ? `${gatePhName(U.phPick)} highlighted` : "";
    $("selMsg").textContent = U.sel
      ? `${selTxt}, lit in every view` + (U.uxy && U.selEmb < U.selN
        ? ` (${U.selEmb.toLocaleString()} of them are in the UMAP; the others ${prot ? "are low signal (rims, gaps)" : "have too few counts"}).` : ".")
      : U.pick !== null || ctLit || phLit ? `${selTxt} in every view.`
      : noEmb ? `Draw a box or a loop around ${W.bins} on the slide to select them.`
      : `Draw a box or a loop around dots, on the tissue or the UMAP, to light the same ${W.bins} in both.`;
    $("selChip").hidden = !has && !lpTxt;
    $("selChipTxt").textContent = selTxt || lpTxt;
    gateCountsRefresh();
  }

  // Saved with Circle -> UMAP's JSON: which whole-section clustering, and the names given to clusters.
  function umapSaveState() {
    const U = S.umap;
    if (!U || !U.key || U.clRes == null) return null;
    return { layer: U.key, resolution: U.clRes, names: U.names, hidden: [...U.hidden] };
  }
  function umapLoadState(o) {
    const U = S.umap;
    if (!o || !U || o.layer !== U.key || !o.names || typeof o.names !== "object") return;
    U.names = o.names; umapDirty(); syncUmapUI(); render();
  }

  function umapInit() {
    S.umap = umapBlank();
    umapSelEvents();
    $("viewSwitch").querySelectorAll("button[data-v]").forEach((b) => { b.onclick = () => umapSetView(b.dataset.v); });
    $("flyPlay").onclick = umapPlay;
    try { const v = localStorage.getItem(FLY_KEY); if (v && $("flySpeed").querySelector(`option[value="${v}"]`)) $("flySpeed").value = v; } catch (e) {}
    $("flySpeed").onchange = (e) => { try { localStorage.setItem(FLY_KEY, e.target.value); } catch (err) {} };
    $("colorBy").onchange = async (e) => {
      const U = S.umap;
      U.colorBy = e.target.value;
      if (U.colorBy === "cluster" && !U.cl) {
        try { await umapSetRes(U.res); }
        catch (err) { status(`Could not load clusters: ${err.message}`, true); setTimeout(() => status(null), 5000); U.colorBy = "gene"; }
      }
      if (U.colorBy === "phenotype") { await gateSetColour(true); return; }
      if (U.colorBy === "celltype" && !U.ct) {
        try { await umapLoadTypes(); }
        catch (err) { status(`Could not load cell types: ${err.message}`, true); setTimeout(() => status(null), 5000); U.colorBy = "gene"; }
      }
      umapDirty(); syncUmapUI(); render();
    };
    $("clusterRes").onchange = async (e) => {
      try { await umapSetRes(clFromValue(e.target.value)); }
      catch (err) { status(`Could not load clusters: ${err.message}`, true); setTimeout(() => status(null), 5000); }
      syncUmapUI(); render();
    };
    $("clShowAll").onclick = () => { S.umap.hidden = new Set(); umapDirty(); syncUmapUI(); render(); };
    $("selBox").onclick = () => umapStartTool("box");
    $("selFree").onclick = () => umapStartTool("free");
    $("selClear").onclick = umapClearSel;
    $("selChipClear").onclick = umapClearSel;
    $("selChipShow").onclick = showLitOnSlide;
    $("selShow").onclick = showLitOnSlide;
    $("selToShape").onclick = umapSelToShape;
    window.addEventListener("keydown", (e) => {
      if (e.key !== "Escape" || !S.umap) return;
      if (e.target && /^(INPUT|TEXTAREA)$/.test(e.target.tagName)) return;
      if (S.umap.tool) umapSetTool(null);
      else if (S.umap.sel || S.umap.pick !== null || S.umap.ctPick !== null || S.umap.phPick != null) umapClearSel();
    });
    // side by side, the circle -> UMAP panel takes the lower half of the UMAP pane: re-lay
    // the panes whenever it opens or closes
    new MutationObserver(() => {
      if (!umapSbs()) return;
      umapApplyViews();
      const [w2, h2] = umapPaneSize("uside");
      if (S.umap.uxy) S.umap.sideVs = umapFitVs(w2, h2);
      deckgl.setProps({ viewState: deckVs() }); render();
    }).observe($("lassoPanel"), { attributes: true, attributeFilter: ["hidden"] });
    // for testing from the console
    Object.assign(window.__spatialViewer, {
      umap: () => S.umap, umapSetView, umapSelect,
      // draw one still frame of the flight (a: 0 = tissue .. 1 = UMAP); flyPreview(null) ends it
      flyPreview: (a) => {
        const U = S.umap;
        if (a == null) { U.anim = null; S.viewState = U.view === "umap" ? (U.umapVs || umapFitVs()) : (U.tissueVs || imageFitVs()); deckgl.setProps({ viewState: deckVs() }); render(); return; }
        U.anim = "preview";
        const vs0 = U.tissueVs || imageFitVs(), vs1 = U.umapVs || umapFitVs();
        umapFlyFrame(umapCloud("fly"), a, vs0, vs1, a);
      },
    });
  }

  // ---- Counting tools for marker squares: % positive, counts per cluster, combined rules -----
  // For an expression layer with one value per square per marker (protein grid). Each marker
  // gets a cut-off on its squares' mean intensity (0-255, the image's own units): set by hand
  // on a histogram, or from Otsu's suggestion with one click - never silently. Rules then
  // name the squares: the FIRST rule (top to bottom) whose conditions all hold, else
  // "Unassigned". Calls live per positions.bin row, so they colour, highlight and select
  // through the same machinery as clusters, in every view.
  // Only shown for protein layers for now; nothing here runs for a sequencing layer.
  const GATE_VERSION = "20261102gate3";
  const GATE_NONE = 255;                                   // calls: no rule matched
  const GATE_COLORS = ["#ff7f0e", "#e6194b", "#3cb44b", "#4c8dff", "#f032e6", "#42d4f4", "#ffe119", "#b07cff",
                       "#bfef45", "#fabed4", "#2ec4b6", "#dcbeff", "#c8a27a", "#aaffc3", "#ff9e9e", "#9fd3ff"];
  const GATE_GREY = [120, 120, 120];
  // A starting point drawn from common immune / vascular / nerve markers. Only rules whose
  // markers are all on the panel are offered, and the panel is labelled "example".
  // Order matters (first match wins): the more specific rule sits above the general one.
  const GATE_EXAMPLES = [
    ["Treg", ["CD45+", "CD4+", "FoxP3+"]],
    ["Cytotoxic T cell", ["CD45+", "CD8+"]],
    ["Helper T cell", ["CD45+", "CD4+"]],
    ["B cell", ["CD45+", "CD19+"]],
    ["Plasma cell", ["CD138+", "CD38+"]],
    ["Macrophage", ["CD45+", "CD68+"]],
    ["Dendritic cell", ["CD45+", "CD11c+", "HLA-DR+"]],
    ["Lymphatic endothelium", ["CD45-", "CD31+", "Podoplanin+"]],
    ["Endothelial", ["CD45-", "CD31+"]],
    ["Axons (NF-L)", ["NF-L+"]],
    ["Glia / Schwann (GFAP)", ["GFAP+"]],
  ];
  const GATE_GAP_DIP = 1.5, GATE_GAP_ETA = 0.7, GATE_POS_HI = 0.6, GATE_POS_LO = 0.005;
  const gateBlank = () => ({ key: null, pos: null, N: 0, rowIdx: null, markers: [], vals: new Map(), info: new Map(),
    cut: new Map(), rules: [], cur: null, preview: false, calls: null, counts: null, ver: 0,
    source: null, busy: null, masks: new Map(), tableKey: "", table: null, posT: null, clT: null });
  const gateOn = () => !!(S.gate && S.gate.key && S.gate.key === S.exprKey && isProt());
  const gateFmt = (v) => (v == null ? "–" : v >= 100 ? v.toFixed(0) : v.toFixed(1));
  const gateRuleRgb = (i) => (i === GATE_NONE || !S.gate.rules[i] ? GATE_GREY : hexToRgb255(S.gate.rules[i].hex));
  const gatePhName = (i) => (i === GATE_NONE ? "Unassigned" : (S.gate.rules[i] && S.gate.rules[i].name) || `Rule ${i + 1}`);
  // the outlines belong to the % positive tool, while it is open on a marker with a line
  const gatePreviewOn = () => !!(gateOn() && S.gateUi && S.gateUi.open && S.gateUi.mode === "pos" &&
    S.gate.preview && S.gate.cur && S.gate.cut.has(S.gate.cur) && S.gate.vals.has(S.gate.cur));
  const umapByPh = () => !!(S.umap && S.umap.colorBy === "phenotype" && gateOn() && S.gate.calls);
  const umapPhLit = () => !!(S.umap && S.umap.phPick != null && gateOn() && S.gate.calls);
  function gateMsg(t) { $("gateMsg").textContent = t || ""; }

  // Called when the expression layer changes (after the UMAP state is rebuilt).
  async function gateOnLayer() {
    S.gate = gateBlank();
    const prot = isProt();
    $("rowGateOpen").hidden = !prot;
    $("gateOpenCl").hidden = !umapClusterings().length;
    $("gatePub").hidden = true;
    $("panelGate").hidden = !prot || !S.gateUi.open;
    syncGateColourOpt();
    if (!prot) return;
    const G = S.gate, e = currentExpr();
    G.key = S.exprKey;
    G.markers = S.geneIndex.map((g) => g.n);
    $("gateCaveat").textContent = `Squares are ${e.binSizeUm ? `${e.binSizeUm} µm` : "small"} and can mix cells: counts are approximate.`;
    if (S.gateUi.open) await gateOpenPanel(S.gateUi.mode);
    // a lab can publish its cut-offs with the data: <dataset>/gates.json (absent = none, quietly)
    const j = await gateFetchPublished();
    if (S.gate !== G || !j) return;
    if (j.layer && j.layer !== G.key) return;
    await gateEnsurePos();
    if (S.gate !== G) return;
    if (await gateApply(j, "this dataset's gates.json (published with the data)")) $("gatePub").hidden = false;
  }
  // fetched once per dataset; only static exports have a dataset folder to look in
  function gateFetchPublished() {
    if (!S.src || !S.src.json) return Promise.resolve(null);
    if (S.gatePub && S.gatePub.src === S.src) return S.gatePub.p;   // per source: a retired source never answers
    const p = S.src.json("gates.json").catch(() => null);
    S.gatePub = { ds: S.dsId, src: S.src, p };
    return p;
  }

  // positions + one value per row for a marker. Gene files of a grid layer hold every square
  // in positions.bin order; if one ever does not, rows are matched by coordinates.
  async function gateEnsurePos() {
    const G = S.gate;
    if (G.pos) return G.pos;
    if (!G.busy) G.busy = (async () => {
      status(`Loading ${wd().bin} positions…`);
      try {
        const pos = S.umap && S.umap.pos && S.umap.key === S.exprKey ? S.umap.pos : await umapPositions();
        if (S.gate !== G) return null;
        G.pos = pos; G.N = pos.length / 2;
        G.rowIdx = Array.from({ length: G.N }, (_, i) => i);
        // the phenotype colours are drawn by the UMAP machinery, which needs the same rows
        if (S.umap && !S.umap.pos) S.umap.pos = pos;
        if (S.umap && !S.umap.key) S.umap.key = S.exprKey;
        return pos;
      } finally { status(null); G.busy = null; }
    })();
    return G.busy;
  }
  async function gateValues(name) {
    const G = S.gate;
    if (G.vals.has(name)) return G.vals.get(name);
    await gateEnsurePos();
    const shown = S.genes.find((g) => g.name === name);
    let xyv = null, n, px;
    if (shown) { n = shown.n; px = shown.positions; }
    else {
      status(`Loading ${name}…`);
      try { xyv = await S.src.gene(S.exprKey, name); } finally { status(null); }
      n = xyv.length / 3;
    }
    if (S.gate !== G) return null;
    const N = G.N, P = G.pos, out = new Float32Array(N);
    const X = (i) => (xyv ? xyv[3 * i] : px[2 * i]), Y = (i) => (xyv ? xyv[3 * i + 1] : px[2 * i + 1]);
    const V = (i) => (xyv ? xyv[3 * i + 2] : shown.values[i]);
    let same = n === N;
    for (let i = 0; same && i < n; i++) if (X(i) !== P[2 * i] || Y(i) !== P[2 * i + 1]) same = false;
    if (same) for (let i = 0; i < N; i++) out[i] = V(i);
    else {
      const m = new Map();
      for (let r = 0; r < N; r++) m.set(ukey(P[2 * r], P[2 * r + 1]), r);
      for (let i = 0; i < n; i++) { const r = m.get(ukey(X(i), Y(i))); if (r !== undefined) out[r] = V(i); }
    }
    G.vals.set(name, out);
    G.info.set(name, gateStats(out));
    return out;
  }

  // Otsu's threshold over the tissue squares (256 bins up to the 99.9th percentile, the very
  // brightest clipped into the top bin), plus two plain checks on whether it means anything:
  //   eta  = between-class / total variance at that threshold (1 = two perfectly separate groups)
  //   dip  = is there a valley? the highest bar above the threshold vs the lowest bar between
  //          the dim peak and it (64 smoothed bars). A tail that just fades out has dip ~ 1.
  function gateStats(v) {
    const N = v.length, s = Float32Array.from(v).sort();
    const q = (f) => s[Math.min(N - 1, Math.max(0, Math.floor(f * (N - 1))))];
    const top = Math.max(q(0.999), 1e-3), NB = 256, h = new Float64Array(NB);
    for (let i = 0; i < N; i++) h[Math.min(NB - 1, Math.max(0, Math.floor((v[i] / top) * NB)))]++;
    let muT = 0; for (let b = 0; b < NB; b++) muT += h[b] * ((b + 0.5) / NB) * top;
    muT /= N;
    let sT = 0; for (let b = 0; b < NB; b++) sT += h[b] * (((b + 0.5) / NB) * top - muT) ** 2;
    sT /= N;
    let w = 0, mu = 0, best = -1, kb = 0;
    for (let b = 0; b < NB - 1; b++) {
      const p = h[b] / N; w += p; mu += p * ((b + 0.5) / NB) * top;
      if (w <= 0 || w >= 1) continue;
      const sb = (muT * w - mu) ** 2 / (w * (1 - w));
      if (sb > best) { best = sb; kb = b; }
    }
    const otsu = ((kb + 1) / NB) * top, eta = sT > 0 ? best / sT : 0;
    // valley check on 64 smoothed bars; the top bar is left out (it holds the clipped brightest)
    const NS = 64, h2 = new Float64Array(NS);
    for (let i = 0; i < N; i++) h2[Math.min(NS - 1, Math.max(0, Math.floor((v[i] / top) * NS)))]++;
    const sm = new Float64Array(NS);
    for (let b = 0; b < NS; b++) sm[b] = (h2[Math.max(0, b - 1)] + h2[b] + h2[Math.min(NS - 1, b + 1)]) / 3;
    const kt = Math.min(NS - 2, Math.floor((otsu / top) * NS));
    let L = 0; for (let b = 0; b <= kt; b++) if (sm[b] > sm[L]) L = b;
    let R = kt + 1; for (let b = kt + 1; b < NS - 1; b++) if (sm[b] > sm[R]) R = b;
    let V = Infinity; for (let b = L; b <= R; b++) V = Math.min(V, sm[b]);
    const dip = sm[R] / Math.max(V, 1e-9);
    let pos = 0; for (let i = 0; i < N; i++) if (v[i] > otsu) pos++;
    // display histogram: 96 bars from 0 to a little past the bright tail (or the cut-off)
    const xmax = Math.min(255, Math.max(q(0.995) * 1.25, otsu * 1.5, 1));
    return { otsu, eta, dip, fracOtsu: pos / N, xmax, median: q(0.5), p99: q(0.99) };
  }
  function gateHistBars(G, name) {
    const I = G.info.get(name), v = G.vals.get(name), NB = 96;
    if (I.bars && I.barsMax === I.xmax) return I.bars;
    const h = new Float64Array(NB);
    for (let i = 0; i < v.length; i++) h[Math.min(NB - 1, Math.max(0, Math.floor((v[i] / I.xmax) * NB)))]++;
    I.bars = h; I.barsMax = I.xmax;
    return h;
  }
  function gateWarnings(I) {
    const out = [];
    if (I.dip < GATE_GAP_DIP) out.push("No clear gap: the squares fade from dim to bright without a dip, so Otsu only splits the tail. Check by eye on the tissue.");
    else if (I.eta < GATE_GAP_ETA) out.push(`Weak split (separation ${I.eta.toFixed(2)}; 1 = two clean groups). Check by eye on the tissue.`);
    if (I.fracOtsu > GATE_POS_HI) out.push("Over 60% of the squares would be positive: check by eye.");
    if (I.fracOtsu < GATE_POS_LO) out.push("Under 0.5% of the squares would be positive: check by eye.");
    return out;
  }
  const gateFracPos = (v, c) => { let k = 0; for (let i = 0; i < v.length; i++) if (v[i] > c) k++; return v.length ? k / v.length : 0; };

  // ---- calls: first matching rule per square ----------------------------------------------
  function gateRuleProblems(rule) {
    const G = S.gate, out = [];
    if (!rule.conds.length) out.push("no conditions yet");
    const off = rule.conds.filter((c) => !G.markers.includes(c.m)).map((c) => c.m);
    if (off.length) out.push(`not on this panel: ${off.join(", ")}`);
    const nocut = rule.conds.filter((c) => G.markers.includes(c.m) && !G.cut.has(c.m)).map((c) => c.m);
    if (nocut.length) out.push(`no cut-off yet for ${nocut.join(", ")}`);
    return out;
  }
  function gateRecompute() {
    const G = S.gate;
    if (!G.pos) return;
    const N = G.N, calls = new Uint8Array(N).fill(GATE_NONE);
    const act = G.rules.slice(0, GATE_NONE).map((rule) => (gateRuleProblems(rule).length || rule.conds.some((c) => !G.vals.has(c.m))
      ? null : rule.conds.map((c) => ({ v: G.vals.get(c.m), t: G.cut.get(c.m).v, pos: c.s > 0 }))));
    if (act.some(Boolean)) {
      for (let r = 0; r < N; r++) {
        for (let i = 0; i < act.length; i++) {
          const a = act[i]; if (!a) continue;
          let ok = true;
          for (let j = 0; j < a.length; j++) { if ((a[j].v[r] > a[j].t) !== a[j].pos) { ok = false; break; } }
          if (ok) { calls[r] = i; break; }
        }
      }
    }
    const counts = new Array(256).fill(0);
    for (let r = 0; r < N; r++) counts[calls[r]]++;
    G.calls = calls; G.counts = counts; G.ver++;
    if (S.umap && S.umap.phPick != null && S.umap.phPick !== GATE_NONE && !G.rules[S.umap.phPick]) S.umap.phPick = null;
    umapDirty();
  }
  // everything that follows a change of cut-off or rules: calls, the panel, the legend, the picture
  function gateChanged() {
    gateRecompute();
    buildGateUI();
    syncUmapUI(); render();
  }
  function gateSetCut(name, v, how) {
    const G = S.gate;
    if (v == null || !Number.isFinite(v)) G.cut.delete(name);
    else G.cut.set(name, { v: Math.max(0, Math.min(255, Math.round(v * 10) / 10)), how });
  }

  // ---- the panel: one tool at a time ---------------------------------------------------------
  //   pos      % positive for one marker: a line on its histogram, % above it per region
  //   clusters counts per region of the current clusters, under their names
  //   rules    optional: name cell types from several markers' lines, count them per region
  const GATE_TITLES = { pos: "% positive for a marker", clusters: "Count clusters", rules: "Combine markers" };
  function buildGateUI() {
    if (!gateOn() || !S.gateUi.open) return;
    const mode = S.gateUi.mode, G = S.gate;
    $("panelGate").hidden = false;
    $("gateTitle").textContent = GATE_TITLES[mode];
    $("gateModePos").hidden = mode !== "pos"; $("gateModeCl").hidden = mode !== "clusters"; $("gateModeRules").hidden = mode !== "rules";
    $("gateSource").hidden = !G.source; $("gateSource").textContent = G.source ? `Cut-offs loaded from ${G.source}.` : "";
    if (mode === "pos") { buildGateMarkerSelect(); buildGateEdit(); }
    if (mode === "rules") buildGateRules();
    gateCountsRefresh(true);
    syncGateColourOpt();
  }
  async function gateOpenPanel(mode) {
    const G = S.gate;
    S.gateUi.open = !!mode;
    if (mode) S.gateUi.mode = mode;
    gateMsg("");
    if (!mode || !gateOn()) { $("panelGate").hidden = true; render(); return; }
    try { await gateEnsurePos(); } catch (e) { gateMsg(`Could not load the ${wd().bins}: ${e.message}`); return; }
    if (mode === "clusters" && !S.umap.cl && umapClusterings().length) {
      try { await umapSetRes(S.umap.res); } catch (e) { gateMsg(`Could not load the clusters: ${e.message}`); }
    }
    if (S.gate !== G) return;
    if (!G.calls) gateRecompute();
    buildGateUI(); syncUmapUI();
    $("panelGate").scrollIntoView({ block: "nearest" });
    render();
  }
  // markers to offer: the panel minus the nuclear counterstain; the status is in the menu
  const gateListMarkers = () => S.gate.markers.filter((m) => !NUCLEAR_CHANNEL.test(m));
  function buildGateMarkerSelect() {
    const G = S.gate, sel = $("gateMarker"), esc = escHtml;
    sel.innerHTML = `<option value="">Choose a marker…</option>` + gateListMarkers().map((m) => {
      const c = G.cut.get(m), v = G.vals.get(m);
      const st = c && v ? ` · ${(gateFracPos(v, c.v) * 100).toFixed(1)}% positive` : c ? " · set" : "";
      return `<option value="${esc(m)}">${esc(m)}${st}</option>`;
    }).join("");
    sel.value = G.cur || "";
  }
  async function gateSelect(m) {
    const G = S.gate;
    G.cur = m || null;
    if (G.cur) {
      try { await gateValues(m); } catch (e) { status(`Could not load ${m}: ${e.message}`, true); setTimeout(() => status(null), 4000); G.cur = null; }
      if (S.gate !== G) return;
      G.preview = true;                                       // outlines on by default while setting
      if (G.cur) gateShowStain(G.cur);                        // the stain is what the line is judged against
    }
    buildGateMarkerSelect(); buildGateEdit(); gateCountsRefresh(true); render();
  }
  function gateShowStain(m) {
    const ch = S.channels.find((c) => c.name.toLowerCase() === m.toLowerCase());
    if (!ch || ch.visible) return;
    ch.visible = true;
    const rows = $("chanList").querySelectorAll(".chan"), i = S.channels.indexOf(ch);
    if (rows[i]) { rows[i].querySelector(".cb").checked = true; rows[i].classList.remove("off"); }
  }
  function gateSentence() {
    const G = S.gate, m = G.cur, c = G.cut.get(m), v = G.vals.get(m), I = G.info.get(m);
    const el = $("gateSentence");
    el.textContent = c ? `${(gateFracPos(v, c.v) * 100).toFixed(1)}% of tissue squares are above the line` +
                         (G.preview ? " — they're outlined on the image." : ".")
                       : "Drag on the chart to place the line where real staining starts, or press Suggest.";
    el.title = (c ? `Line at ${gateFmt(c.v)} (mean intensity per square, 0–255 as the image; ${c.how === "otsu" ? "the suggestion" : "set by hand"}). ` : "") +
      `Suggestion (Otsu): ${gateFmt(I.otsu)}, ${(I.fracOtsu * 100).toFixed(1)}% above. Typical square ${gateFmt(I.median)}; ` +
      `brightest 1% above ${gateFmt(I.p99)}. Separation ${I.eta.toFixed(2)} (1 = two clean groups).` +
      (gateWarnings(I).length ? " " + gateWarnings(I).join(" ") : "");
  }
  function buildGateEdit() {
    const G = S.gate, m = G.cur, edit = $("gateEdit");
    edit.hidden = !m || !G.vals.has(m);
    if (edit.hidden) return;
    const I = G.info.get(m), c = G.cut.get(m);
    if (c && c.v > I.xmax * 0.95) I.xmax = Math.min(255, c.v * 1.3);
    // the one warning worth a line of its own: no valley, or a split that leaves almost all / none positive
    const w = I.dip < GATE_GAP_DIP || I.eta < GATE_GAP_ETA || I.fracOtsu > GATE_POS_HI || I.fracOtsu < GATE_POS_LO;
    $("gateNoGap").hidden = !w;
    $("gateNoGap").title = gateWarnings(I).join(" ");
    $("gatePreview").checked = !!G.preview;
    $("gateClearCut").hidden = !c;
    $("gateAxisHi").textContent = `bright (${gateFmt(I.xmax)})`;
    gateSentence();
    drawGateHist();
  }
  function drawGateHist() {
    const G = S.gate, m = G.cur, cv = $("gateHist");
    if (!m || !G.vals.has(m) || $("gateEdit").hidden) return;
    const I = G.info.get(m), c = G.cut.get(m), bars = gateHistBars(G, m);
    const dpr = window.devicePixelRatio || 1, W = Math.max(100, cv.clientWidth || 300), H = cv.clientHeight || 96;
    cv.width = Math.round(W * dpr); cv.height = Math.round(H * dpr);
    const g = cv.getContext("2d");
    g.setTransform(dpr, 0, 0, dpr, 0, 0); g.clearRect(0, 0, W, H);
    // bar height on a log scale: the dim background holds most squares by far
    let top = 0; for (const b of bars) top = Math.max(top, Math.log1p(b));
    const NB = bars.length, bw = W / NB, cx = c ? (c.v / I.xmax) * W : -1;
    for (let b = 0; b < NB; b++) {
      if (!bars[b]) continue;
      const bh = Math.max(1, (Math.log1p(bars[b]) / top) * (H - 14));
      g.fillStyle = c && (b + 0.5) * bw > cx ? `rgba(${gateLineRgb(m).join(",")},0.9)` : "rgba(160,170,190,0.45)";
      g.fillRect(b * bw, H - bh, Math.max(1, bw - 0.6), bh);
    }
    if (c) {
      g.strokeStyle = "#fff"; g.lineWidth = 2;
      g.beginPath(); g.moveTo(cx, 0); g.lineTo(cx, H); g.stroke();
      g.fillStyle = "#fff"; g.font = "11px system-ui, sans-serif";
      const t = gateFmt(c.v), tw = g.measureText(t).width;
      g.fillText(t, cx + 4 + tw > W ? cx - tw - 4 : cx + 4, 11);
    }
  }
  function gateHistEvents() {
    const cv = $("gateHist");
    let drag = false;
    const at = (ev) => {
      const G = S.gate, I = G.info.get(G.cur); if (!I) return null;
      const r = cv.getBoundingClientRect();
      return Math.max(0, Math.min(1, (ev.clientX - r.left) / r.width)) * I.xmax;
    };
    const live = (v, how) => {
      const G = S.gate; if (v == null) return;
      gateSetCut(G.cur, v, how || "manual");
      $("gateClearCut").hidden = false;
      gateSentence(); drawGateHist(); umapDirty(); requestRender();
    };
    cv.addEventListener("pointerdown", (ev) => {
      if (ev.button !== 0 || !S.gate || !S.gate.cur) return;
      drag = true; cv.setPointerCapture(ev.pointerId); live(at(ev));
    });
    cv.addEventListener("pointermove", (ev) => { if (drag) live(at(ev)); });
    const end = () => { if (!drag) return; drag = false; gateChanged(); };
    cv.addEventListener("pointerup", end);
    cv.addEventListener("pointercancel", end);
    // arrow keys nudge the line by half a percent of the axis
    cv.addEventListener("keydown", (ev) => {
      const G = S.gate, I = G && G.info.get(G.cur); if (!I || !/^Arrow(Left|Right)$/.test(ev.key)) return;
      ev.preventDefault(); ev.stopPropagation();
      const c = G.cut.get(G.cur), step = Math.max(0.1, I.xmax / 200) * (ev.key === "ArrowLeft" ? -1 : 1);
      live((c ? c.v : I.otsu) + step); gateChanged();
    });
  }

  // ---- Combine markers: cell types as sentences ------------------------------------------------
  // "[Macrophage] = CD45 + and CD68 +". The marker menu lists only markers that have a line.
  const gateCutMarkers = () => S.gate.markers.filter((m) => S.gate.cut.has(m));
  function buildGateRules() {
    const G = S.gate, host = $("gateRules"), esc = escHtml, cutM = gateCutMarkers();
    host.innerHTML = "";
    $("gateAddRule").disabled = !cutM.length;
    $("gateAddRule").title = cutM.length ? "Add a cell type at the bottom" : "Set a line for at least one marker first (← Back to markers)";
    const ex = gateExampleRules();
    $("gateExamples").disabled = !ex.length;
    $("gateExamples").title = ex.length ? `Adds ${ex.map((r) => r[0]).join(", ")} (examples: edit them for your panel)`
      : "No example fits the markers you have set yet (examples use CD45, CD4, CD8, CD68, CD31, NF-L …)";
    if (!G.rules.length) {
      host.innerHTML = `<div class="meta gate-note">${cutM.length ? `No cell types yet. Markers with a line: ${esc(cutM.join(", "))}.` : "Set a line for at least one marker first (← Back to markers)."}</div>`;
      return;
    }
    G.rules.forEach((rule, i) => {
      const probs = gateRuleProblems(rule);
      const row = document.createElement("div");
      row.className = "gr-row" + (probs.length && rule.conds.length ? " off" : "");
      const n = G.counts && !probs.length ? G.counts[i] : null;
      const conds = rule.conds.map((c, j) => `<span class="cond">${esc(c.m)} <button class="sg" data-j="${j}" title="Click to switch between + (above its line) and − (below)">${c.s > 0 ? "+" : "−"}</button>` +
        `<button class="rmc" data-j="${j}" title="Remove ${esc(c.m)}">×</button></span>`).join(`<span class="and">and</span>`);
      const free = cutM.filter((m) => !rule.conds.some((c) => c.m === m));
      row.innerHTML = `<input type="color" value="${rule.hex}" title="Colour" />` +
        `<div class="gr-txt"><input type="text" class="nm" maxlength="60" /><span class="eq">=</span>${conds}` +
        (free.length ? `<select class="addc" title="Add a marker"><option value="">${rule.conds.length ? "+ and…" : "+ marker…"}</option>${free.map((m) => `<option value="${esc(m)}">${esc(m)}</option>`).join("")}</select>` : "") +
        (n != null ? ` <span class="gr-n" title="Squares called ${esc(rule.name)}">· ${n.toLocaleString()}</span>` : "") +
        (probs.length && rule.conds.length ? `<span class="gr-skip">skipped: ${esc(probs.join("; "))}</span>` : "") + `</div>` +
        `<div class="gr-btns"><button class="up" title="Move up (checked earlier)" ${i ? "" : "disabled"}>↑</button>` +
        `<button class="dn" title="Move down" ${i < G.rules.length - 1 ? "" : "disabled"}>↓</button><button class="del" title="Delete">×</button></div>`;
      const nm = row.querySelector(".nm");
      nm.value = rule.name;
      nm.onchange = () => { rule.name = nm.value.trim().slice(0, 60) || `Cell type ${i + 1}`; G.ver++; umapDirty(); syncUmapUI(); render(); };
      nm.onkeydown = (e) => { e.stopPropagation(); if (e.key === "Enter") nm.blur(); };
      row.querySelector("input[type=color]").oninput = (e) => { rule.hex = e.target.value; G.ver++; umapDirty(); syncUmapUI(); requestRender(); };
      row.querySelector(".up").onclick = () => { G.rules.splice(i - 1, 0, G.rules.splice(i, 1)[0]); gateRulesEdited(); };
      row.querySelector(".dn").onclick = () => { G.rules.splice(i + 1, 0, G.rules.splice(i, 1)[0]); gateRulesEdited(); };
      row.querySelector(".del").onclick = () => { G.rules.splice(i, 1); gateRulesEdited(); };
      row.querySelectorAll(".sg").forEach((b) => { b.onclick = () => { const c = rule.conds[+b.dataset.j]; c.s = -c.s; gateRulesEdited(); }; });
      row.querySelectorAll(".rmc").forEach((b) => { b.onclick = () => { rule.conds.splice(+b.dataset.j, 1); gateRulesEdited(); }; });
      const add = row.querySelector(".addc");
      if (add) add.onchange = (e) => { const m = e.target.value; if (m) rule.conds.push({ m, s: 1 }); gateRulesEdited(); };
      host.appendChild(row);
    });
  }
  // the examples whose markers all have a line (and are not in the list already)
  function gateExampleRules() {
    const G = S.gate, byLow = new Map(G.markers.map((m) => [m.toLowerCase(), m]));
    const have = new Set(G.rules.map((r) => r.name.toLowerCase()));
    const out = [];
    for (const [name, conds] of GATE_EXAMPLES) {
      const cs = conds.map((c) => ({ m: byLow.get(c.slice(0, -1).toLowerCase()), s: c.endsWith("-") ? -1 : 1 }));
      if (cs.every((c) => c.m && G.cut.has(c.m)) && !have.has(name.toLowerCase())) out.push([name, cs]);
    }
    return out;
  }
  async function gateRulesEdited() {
    const G = S.gate;
    if (S.umap) S.umap.phPick = null;                  // rule numbers may have moved
    // a rule whose markers have lines needs their values (loaded on demand)
    const need = new Set();
    for (const r of G.rules) for (const c of r.conds) if (G.cut.has(c.m) && !G.vals.has(c.m)) need.add(c.m);
    for (const m of need) { try { await gateValues(m); } catch (e) { /* stays skipped */ } }
    if (S.gate !== G) return;
    gateChanged();
  }

  // ---- regions: whole tissue, Tier A regions, overlays switched on, the selection -------------
  function gateOverlayMask(o) {
    const G = S.gate, k = `${o.key || o.file}|${G.N}`;
    if (G.masks.has(k)) return G.masks.get(k);
    const groups = new Map();                          // region name -> Uint8Array over rows
    const feats = (o.data && (o.data.features || (o.data.type === "Feature" ? [o.data] : []))) || [];
    const P = G.pos, N = G.N;
    for (const f of feats) {
      const geo = f && f.geometry; if (!geo) continue;
      const polys = geo.type === "Polygon" ? [geo.coordinates] : geo.type === "MultiPolygon" ? geo.coordinates : [];
      if (!polys.length) continue;
      const p = f.properties || {};
      const name = String(p.name || (p.classification && p.classification.name) || p.label || p.id || o.label || "Region");
      let mk = groups.get(name); if (!mk) { mk = new Uint8Array(N); groups.set(name, mk); }
      for (const rings of polys) {
        let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity;
        for (const [x, y] of rings[0]) { x0 = Math.min(x0, x); y0 = Math.min(y0, y); x1 = Math.max(x1, x); y1 = Math.max(y1, y); }
        for (let r = 0; r < N; r++) {
          const x = P[2 * r], y = P[2 * r + 1];
          if (mk[r] || x < x0 || x > x1 || y < y0 || y > y1) continue;
          let inside = false;                            // even-odd over the rings: holes come out
          for (const ring of rings) if (inPoly(x, y, ring)) inside = !inside;
          if (inside) mk[r] = 1;
        }
      }
    }
    const out = [...groups].map(([name, mk]) => ({ name: `${o.label}: ${name}`, mask: mk }));
    G.masks.set(k, out);
    return out;
  }
  function gateGroups() {
    const G = S.gate, U = S.umap, out = [{ name: "Whole tissue", test: null }];
    const R = S.region;
    if (R.meta && R.ids && R.ids.length === G.N) {
      for (const m of R.meta) { const id = m.id; out.push({ name: m.name || `Region ${id}`, test: (r) => R.ids[r] === id }); }
    }
    for (const o of S.overlays) {
      if (!o.on || !o.data) continue;
      if (R.meta && R.ids && (o.file === R.overlayKey || o.key === R.overlayKey)) continue;   // already counted by region.bin
      for (const g of gateOverlayMask(o)) { const mk = g.mask; out.push({ name: g.name, test: (r) => mk[r] === 1 }); }
    }
    if (U && U.sel && U.sel.length === G.N) { const sel = U.sel; out.push({ name: `Selection (${U.selN.toLocaleString()} ${wd().bins})`, test: (r) => sel[r] === 1, sel: true }); }
    return out;
  }
  // counts per group x phenotype; rows = rules in order, then Unassigned
  function gateTable() {
    const G = S.gate;
    if (!G.calls) return null;
    const key = `${G.ver}|${gateGroupKey()}`;
    if (G.table && G.tableKey === key) return G.table;
    const groups = gateGroups(), K = G.rules.length, calls = G.calls;
    for (const g of groups) {
      const c = new Float64Array(K + 1); let tot = 0;
      for (let r = 0; r < G.N; r++) {
        if (g.test && !g.test(r)) continue;
        const k = calls[r]; c[k === GATE_NONE ? K : k]++; tot++;
      }
      g.counts = c; g.total = tot;
    }
    G.table = { groups, K }; G.tableKey = key;
    return G.table;
  }
  // what the region columns depend on: the selection, overlays switched on, region.bin
  function gateGroupKey() {
    const U = S.umap;
    return `${U && U.sel ? `${U.selN}:${U.selEmb}` : ""}|${S.overlays.filter((o) => o.on && o.data).map((o) => o.key || o.file).join(",")}|${S.region.ids ? 1 : 0}`;
  }
  // squares above one marker's line, per region
  function gatePosCounts(m) {
    const G = S.gate, c = G.cut.get(m), v = G.vals.get(m);
    if (!c || !v || !G.pos) return null;
    const groups = gateGroups(), t = c.v;
    for (const g of groups) {
      let pos = 0, tot = 0;
      for (let r = 0; r < G.N; r++) { if (g.test && !g.test(r)) continue; tot++; if (v[r] > t) pos++; }
      g.pos = pos; g.total = tot;
    }
    return { m, cut: c, groups };
  }
  function gatePosTable() {
    const G = S.gate, m = G.cur, c = m && G.cut.get(m);
    if (!c) return null;
    const key = `${m}|${c.v}|${gateGroupKey()}`;
    if (!(G.posT && G.posT.key === key)) G.posT = { key, ...gatePosCounts(m) };
    return G.posT;
  }
  // squares per cluster (the clustering on show), per region; the last column is "not in a cluster"
  function gateClTable() {
    const U = S.umap, G = S.gate;
    if (!U || !U.cl || !G.pos || U.cl.length !== G.N) return null;
    const key = `${U.clRes}|${gateGroupKey()}`;
    if (G.clT && G.clT.key === key) return G.clT;
    let K = 0; for (let r = 0; r < G.N; r++) if (U.cl[r] !== UMAP_NONE && U.cl[r] + 1 > K) K = U.cl[r] + 1;
    const groups = gateGroups(), cl = U.cl;
    for (const g of groups) {
      const cnt = new Float64Array(K + 1); let tot = 0;
      for (let r = 0; r < G.N; r++) { if (g.test && !g.test(r)) continue; const c = cl[r]; cnt[c === UMAP_NONE ? K : c]++; tot++; }
      g.counts = cnt; g.total = tot;
    }
    G.clT = { key, K, res: U.clRes, groups };
    return G.clT;
  }
  const gateColHead = (groups) => `<tr><th></th>${groups.map((g) => `<th title="${escHtml(g.name)}">${escHtml(g.name.length > 22 ? g.name.slice(0, 21) + "…" : g.name)}</th>`).join("")}</tr>`;
  const gatePct = (a, b) => (b ? ((100 * a) / b).toFixed(1) : "0.0");
  const gateMoreCols = (n, what = "column") => (n === 1 ? ` A Box or Freehand selection (or a region overlay) adds a ${what}.` : "");
  // redraws the table of whichever tool is open (cheap when nothing it counts has changed)
  function gateCountsRefresh(force) {
    const G = S.gate;
    if (!G || !gateOn() || !S.gateUi || !S.gateUi.open) return;
    if (force) { G.tableKey = ""; }
    const mode = S.gateUi.mode, esc = escHtml, U = S.umap, W = wd();
    if (mode === "pos") {
      const host = $("gatePosTable"), T = gatePosTable();
      $("gatePosCsv").disabled = !gateCutMarkers().length;
      if (!T) { host.innerHTML = ""; return; }
      host.innerHTML = `<table><thead><tr><th></th><th>Above the line</th><th>Tissue ${W.bins}</th><th>%</th></tr></thead><tbody>` +
        T.groups.map((g) => `<tr><td class="lbl" title="${esc(g.name)}">${esc(g.name)}</td><td>${g.pos.toLocaleString()}</td><td>${g.total.toLocaleString()}</td><td>${gatePct(g.pos, g.total)}%</td></tr>`).join("") +
        `</tbody></table><div class="meta gate-note">${esc(T.m)} above ${gateFmt(T.cut.v)}.${gateMoreCols(T.groups.length, "row")}` +
        (gateCutMarkers().length > 1 ? ` The CSV holds every marker with a line: ${esc(gateCutMarkers().join(", "))}.` : "") + `</div>`;
      return;
    }
    if (mode === "clusters") {
      const host = $("gateClTable"), list = umapClusterings(), rs = $("gateClRes");
      const opts = list.map((c) => `<option value="${escHtml(clKey(c))}">${escHtml(clLabel(c))} · ${c.n} ${clUser(c) ? "groups" : "clusters"}</option>`).join("");
      if (rs.dataset.opts !== opts) { rs.innerHTML = opts; rs.dataset.opts = opts; }
      if (U && U.clRes != null) rs.value = String(U.clRes);
      $("gateClColour").checked = !!(U && U.colorBy === "cluster");
      const T = gateClTable();
      $("gateClCsv").disabled = !T;
      if (!T) { host.innerHTML = `<div class="meta gate-note">${list.length ? "Loading the clusters…" : "This layer has no clusters."}</div>`; return; }
      const rows = [];
      for (let k = 0; k <= T.K; k++) {
        const tot = T.groups.reduce((a, g) => a + g.counts[k], 0);
        if (k === T.K && !tot) continue;
        const none = k === T.K, [r, gg, b] = none ? GATE_GREY : sectionRgb(k);
        rows.push(`<tr class="ph${!none && U.pick === k ? " sel" : ""}" data-k="${none ? -1 : k}"><td title="${none ? "Squares left out of the UMAP" : "Click to highlight in every view"}"><i style="background:rgb(${r},${gg},${b})"></i>${esc(none ? "Not in a cluster" : umapName(k))}</td>` +
          T.groups.map((g) => `<td>${g.counts[k].toLocaleString()}<small>${gatePct(g.counts[k], g.total)}%</small></td>`).join("") + `</tr>`);
      }
      rows.push(`<tr class="tot"><td>All ${W.bins}</td>${T.groups.map((g) => `<td>${g.total.toLocaleString()}</td>`).join("")}</tr>`);
      host.innerHTML = `<table><thead>${gateColHead(T.groups)}</thead><tbody>${rows.join("")}</tbody></table>` +
        `<div class="meta gate-note">Names as in Clusters &amp; selection (rename with ✎ there). Click a cluster to highlight it.${gateMoreCols(T.groups.length)}</div>`;
      host.querySelectorAll("tr.ph").forEach((tr) => { const k = +tr.dataset.k; if (k >= 0) tr.onclick = () => umapPick(k); });
      return;
    }
    // combined rules
    const host = $("gateCounts");
    const T = G.rules.length ? gateTable() : null;
    $("gateCsvCounts").disabled = !T; $("gateCsvSquares").disabled = !T;
    if (!T) { host.innerHTML = ""; return; }
    const rows = [];
    for (let k = 0; k <= T.K; k++) {
      const id = k === T.K ? GATE_NONE : k, [r, gg, b] = gateRuleRgb(id);
      rows.push(`<tr class="ph${U && U.phPick === id ? " sel" : ""}" data-k="${id}"><td title="Click to highlight in every view"><i style="background:rgb(${r},${gg},${b})"></i>${esc(gatePhName(id))}</td>` +
        T.groups.map((g) => `<td>${g.counts[k].toLocaleString()}<small>${gatePct(g.counts[k], g.total)}%</small></td>`).join("") + `</tr>`);
    }
    rows.push(`<tr class="tot"><td>All ${W.bins}</td>${T.groups.map((g) => `<td>${g.total.toLocaleString()}</td>`).join("")}</tr>`);
    host.innerHTML = `<table><thead>${gateColHead(T.groups)}</thead><tbody>${rows.join("")}</tbody></table>` +
      `<div class="meta gate-note">Click a cell type to highlight it.${gateMoreCols(T.groups.length)}</div>`;
    host.querySelectorAll("tr.ph").forEach((tr) => { tr.onclick = () => umapPickPh(+tr.dataset.k); });
  }
  function umapPickPh(k) {
    const U = S.umap; if (!U) return;
    U.phPick = k === null || k === undefined || U.phPick === k ? null : k;
    U.pick = null; U.sel = null; U.selN = 0; U.ctPick = null;
    umapDirty(); syncUmapUI(); render();
  }

  // ---- exports: counts, per-square calls, the gates themselves ---------------------------------
  const csvQ = (x) => `"${String(x).replace(/"/g, '""')}"`;
  function gateCountsCsv() {
    const G = S.gate, T = gateTable(); if (!T) return "";
    const lines = ["dataset,layer,region,cell_type,count,total_squares,percent"];
    for (const g of T.groups) for (let k = 0; k <= T.K; k++) {
      const id = k === T.K ? GATE_NONE : k;
      lines.push([csvQ(S.dsId), csvQ(G.key), csvQ(g.sel ? "Selection" : g.name), csvQ(gatePhName(id)), g.counts[k], g.total,
                  g.total ? ((100 * g.counts[k]) / g.total).toFixed(3) : "0.000"].join(","));
    }
    return lines.join("\n") + "\n";
  }
  // % positive: one row per marker with a line x region
  function gatePosCsv() {
    const G = S.gate, ms = gateCutMarkers().filter((m) => G.vals.has(m));
    if (!ms.length) return "";
    const lines = ["dataset,layer,marker,cutoff,how_set,region,squares_above,total_squares,percent"];
    for (const m of ms) {
      const T = gatePosCounts(m);
      for (const g of T.groups) lines.push([csvQ(S.dsId), csvQ(G.key), csvQ(m), T.cut.v, T.cut.how === "otsu" ? "otsu" : "manual",
        csvQ(g.sel ? "Selection" : g.name), g.pos, g.total, g.total ? ((100 * g.pos) / g.total).toFixed(3) : "0.000"].join(","));
    }
    return lines.join("\n") + "\n";
  }
  function gateClCsv() {
    const G = S.gate, T = gateClTable(); if (!T) return "";
    const lines = ["dataset,layer,leiden_resolution,region,cluster,cluster_name,count,total_squares,percent"];
    for (const g of T.groups) for (let k = 0; k <= T.K; k++) {
      if (k === T.K && !g.counts[k]) continue;
      lines.push([csvQ(S.dsId), csvQ(G.key), T.res, csvQ(g.sel ? "Selection" : g.name), k === T.K ? "" : k,
        csvQ(k === T.K ? "Not in a cluster" : umapName(k)), g.counts[k], g.total, g.total ? ((100 * g.counts[k]) / g.total).toFixed(3) : "0.000"].join(","));
    }
    return lines.join("\n") + "\n";
  }
  function gateSquaresCsv() {
    const G = S.gate; if (!G.calls) return "";
    const ms = G.markers.filter((m) => G.cut.has(m) && G.vals.has(m));
    const cols = ms.map((m) => ({ v: G.vals.get(m), t: G.cut.get(m).v }));
    const names = G.rules.map((r, i) => csvQ(gatePhName(i))), none = csvQ("Unassigned");
    const out = ["row,x_px,y_px,cell_type" + ms.map((m) => "," + csvQ(`${m}_pos`)).join("")];
    const P = G.pos;
    for (let r = 0; r < G.N; r++) {
      let s = `${r},${P[2 * r].toFixed(1)},${P[2 * r + 1].toFixed(1)},${G.calls[r] === GATE_NONE ? none : names[G.calls[r]]}`;
      for (const c of cols) s += c.v[r] > c.t ? ",1" : ",0";
      out.push(s);
    }
    return out.join("\n") + "\n";
  }
  function gateJson() {
    const G = S.gate, e = currentExpr();
    const cutoffs = {};
    for (const m of G.markers) {
      const c = G.cut.get(m); if (!c) continue;
      const I = G.info.get(m);
      cutoffs[m] = { value: c.v, how: c.how === "otsu" ? "otsu" : "manual", ...(I ? { otsu: Math.round(I.otsu * 10) / 10 } : {}) };
    }
    return {
      tool: "SODA gating", viewerVersion: GATE_VERSION, dataset: S.dsId, layer: G.key,
      unit: "mean intensity per square, 0-255 as the image (8-bit); positive = above the cut-off",
      squareUm: e && e.binSizeUm ? e.binSizeUm : null,
      perSlide: "Cut-offs belong to the slide they were set on; check them again before using them on another slide.",
      cutoffs,
      rules: G.rules.map((r) => ({ name: r.name, color: r.hex, conditions: r.conds.map((c) => `${c.m}${c.s > 0 ? "+" : "-"}`) })),
      firstMatchWins: true,
      date: new Date().toISOString(),
    };
  }
  // Apply a gates JSON (a file, or the dataset's own gates.json). Nothing is guessed: cut-offs
  // for markers not on this panel are listed and skipped; rules keep their order and flags.
  async function gateApply(j, source) {
    const G = S.gate;
    if (!j || typeof j !== "object" || (!j.cutoffs && !Array.isArray(j.rules))) { gateMsg("That file holds no cut-offs or rules.", true); return false; }
    const notes = [];
    if (j.dataset && j.dataset !== S.dsId) notes.push(`These cut-offs were set on ${j.dataset}, not this slide: check each one before trusting the calls.`);
    if (j.layer && j.layer !== G.key) notes.push(`They were set on layer ${j.layer}.`);
    const byLow = new Map(G.markers.map((m) => [m.toLowerCase(), m]));
    const skipped = [];
    G.cut = new Map();
    for (const [m0, c] of Object.entries(j.cutoffs || {})) {
      const m = byLow.get(String(m0).toLowerCase()), v = c && typeof c === "object" ? Number(c.value) : Number(c);
      if (!m || !Number.isFinite(v)) { skipped.push(m0); continue; }
      gateSetCut(m, v, c && c.how === "otsu" ? "otsu" : "manual");
    }
    if (skipped.length) notes.push(`Not on this panel, skipped: ${skipped.join(", ")}.`);
    if (Array.isArray(j.rules)) {
      G.rules = j.rules.slice(0, GATE_NONE).map((r, i) => ({
        name: String((r && r.name) || `Rule ${i + 1}`).slice(0, 60),
        hex: /^#[0-9a-f]{6}$/i.test(r && r.color) ? r.color : GATE_COLORS[i % GATE_COLORS.length],
        conds: (Array.isArray(r && r.conditions) ? r.conditions : []).map((s) => {
          const t = String(s).trim().replace(/−/g, "-"), sg = t.endsWith("-") ? -1 : 1, nm = t.replace(/[+-]$/, "");
          return { m: byLow.get(nm.toLowerCase()) || nm, s: sg };
        }).filter((c) => c.m),
      }));
    }
    G.source = source;
    try {
      for (const m of G.cut.keys()) { await gateValues(m); if (S.gate !== G) return false; }
    } catch (e) { notes.push(`Could not load every marker: ${e.message}`); }
    if (S.umap) S.umap.phPick = null;
    gateRecompute();
    buildGateUI();
    gateMsg(notes.join(" "));
    syncUmapUI(); render();
    return true;
  }
  async function gateLoadFile(file) {
    let j;
    try { j = JSON.parse(await file.text()); } catch (e) { gateMsg("Not a JSON file.", true); return; }
    await gateEnsurePos();
    await gateApply(j, `the file ${file.name || "you chose"}`);
  }

  // ---- drawing: positives preview on the tissue; the colour-by option ---------------------------
  // The squares above the cut-off keep the image as it is and get a yellow outline; every other
  // tissue square is darkened, so stained structures inside the outlines are easy to check.
  // The outline must stand out from the marker's own stain colour: magenta, or yellow over a red stain.
  function gateLineRgb(m) {
    const ch = S.channels.find((c) => c.name.toLowerCase() === String(m || "").toLowerCase());
    const [r, g, b] = ch ? ch.color : [0, 0, 0];
    return r > 0.6 && g < 0.5 ? [255, 214, 64] : [255, 72, 220];
  }
  function gatePreviewLayer() {
    const G = S.gate, e = currentExpr(), v = G.vals.get(G.cur), t = G.cut.get(G.cur).v, P = G.pos;
    const ppu = Math.pow(2, (S.viewState && S.viewState.zoom) || 0);
    const side = Math.max(((e && e.spotDiameterPx) || 8) * S.pointScale, HEAT_MIN_PX / ppu);
    const CLEAR = [0, 0, 0, 0], DIM = [0, 0, 0, 170], LINE = [...gateLineRgb(G.cur), 255];
    return new ColumnLayer({
      id: `gateprev-${S.dsId}-${G.key}`, data: G.rowIdx,
      getPosition: (r) => [P[2 * r], P[2 * r + 1], 0],
      getFillColor: (r) => (v[r] > t ? CLEAR : DIM),
      getLineColor: (r) => (v[r] > t ? LINE : CLEAR),
      diskResolution: 4, angle: 45, radius: side / Math.SQRT2, radiusUnits: "common",
      extruded: false, filled: true, stroked: true, lineWidthUnits: "pixels", getLineWidth: 1.25,
      parameters: NORMAL, pickable: false,
      updateTriggers: { getFillColor: [G.cur, t], getLineColor: [G.cur, t, LINE.join()] },
    });
  }
  function syncGateColourOpt() {
    const opt = $("colorBy").querySelector("option[value=phenotype]"), ok = gateOn() && S.gate.rules.length > 0;
    opt.hidden = !ok; opt.disabled = !ok;
    $("gateColour").checked = !!(S.umap && S.umap.colorBy === "phenotype");
    $("gateClColour").checked = !!(S.umap && S.umap.colorBy === "cluster");
  }
  async function gateSetColour(on) {
    const U = S.umap; if (!U) return;
    if (on) { await gateEnsurePos(); if (!S.gate.calls) gateRecompute(); }
    U.colorBy = on ? "phenotype" : "gene";
    umapDirty(); syncUmapUI(); render();
  }

  function gateInit() {
    S.gate = gateBlank();
    S.gateUi = { open: false, mode: "pos" };
    gateHistEvents();
    const toggle = (mode) => gateOpenPanel(S.gateUi.open && S.gateUi.mode === mode ? null : mode);
    $("gateOpenPos").onclick = () => toggle("pos");
    $("gateOpenCl").onclick = () => toggle("clusters");
    $("gateClose").onclick = () => gateOpenPanel(null);
    $("gateCombine").onclick = () => gateOpenPanel("rules");
    $("gateBackPos").onclick = () => gateOpenPanel("pos");
    $("gateMarker").onchange = (e) => gateSelect(e.target.value);
    $("gateSuggest").onclick = () => {
      const G = S.gate, I = G.info.get(G.cur); if (!I) return;
      gateSetCut(G.cur, I.otsu, "otsu"); gateChanged();
    };
    $("gateClearCut").onclick = () => {
      const G = S.gate; if (!G.cur) return;
      gateSetCut(G.cur, null); gateChanged();
    };
    $("gatePreview").onchange = (e) => {
      const G = S.gate; G.preview = e.target.checked;
      if (G.preview && G.cur) { gateShowStain(G.cur); if (S.umap && S.umap.view === "umap") umapSetView("sbs"); }
      gateSentence(); umapDirty(); render();
    };
    $("gateMore").onclick = () => { $("gateMoreBox").hidden = !$("gateMoreBox").hidden; };
    $("gatePosCsv").onclick = () => { const t = gatePosCsv(); if (t) lassoDownload(`${S.dsId}_${S.gate.key}_percent_positive.csv`, t, "text/csv"); };
    $("gateClRes").onchange = async (e) => {
      try { await umapSetRes(clFromValue(e.target.value)); } catch (err) { gateMsg(`Could not load the clusters: ${err.message}`); }
      syncUmapUI(); render();
    };
    $("gateClColour").onchange = async (e) => {
      const U = S.umap;
      if (e.target.checked && !U.cl) { try { await umapSetRes(U.res); } catch (err) { gateMsg(err.message); return; } }
      U.colorBy = e.target.checked ? "cluster" : "gene";
      umapDirty(); syncUmapUI(); render();
    };
    $("gateClCsv").onclick = () => { const t = gateClCsv(); if (t) lassoDownload(`${S.dsId}_${S.gate.key}_cluster_counts_r${S.umap.clRes}.csv`, t, "text/csv"); };
    $("gateAddRule").onclick = () => {
      const G = S.gate;
      G.rules.push({ name: `Cell type ${G.rules.length + 1}`, hex: GATE_COLORS[G.rules.length % GATE_COLORS.length], conds: [] });
      gateRulesEdited().then(() => { const nm = $("gateRules").querySelectorAll("input.nm"); if (nm.length) { nm[nm.length - 1].focus(); nm[nm.length - 1].select(); } });
    };
    $("gateExamples").onclick = () => {
      const G = S.gate, ex = gateExampleRules(); if (!ex.length) return;
      for (const [name, conds] of ex) G.rules.push({ name, hex: GATE_COLORS[G.rules.length % GATE_COLORS.length], conds });
      gateRulesEdited().then(() => gateMsg("Example rules added: edit them for your panel."));
    };
    $("gateColour").onchange = (e) => gateSetColour(e.target.checked);
    $("gateCsvCounts").onclick = () => { const t = gateCountsCsv(); if (t) lassoDownload(`${S.dsId}_${S.gate.key}_celltype_counts.csv`, t, "text/csv"); };
    $("gateCsvSquares").onclick = () => { const t = gateSquaresCsv(); if (t) lassoDownload(`${S.dsId}_${S.gate.key}_celltype_squares.csv`, t, "text/csv"); };
    $("gateSave").onclick = () => lassoDownload(`${S.dsId}_${S.gate.key}_gates.json`, JSON.stringify(gateJson(), null, 1), "application/json");
    $("gateLoadBtn").onclick = () => $("gateLoad").click();
    $("gateLoad").onchange = (e) => { const f = e.target.files[0]; e.target.value = ""; if (f) gateLoadFile(f); };
    window.addEventListener("resize", () => { if (S.gate && S.gate.cur) drawGateHist(); });
    window.addEventListener("soda:shown", () => { if (S.gate && S.gate.cur) drawGateHist(); });   // layout.js: a folded section or another tab opened
    // for testing from the console: the texts the buttons would download, and loading without a file
    Object.assign(window.__spatialViewer, {
      gate: () => S.gate, gateUi: () => S.gateUi, gateSelect, gateOpenPanel,
      gateSetCut: (m, v, how) => { gateSetCut(m, v, how || "manual"); gateChanged(); },
      gatePosCsv, gateClCsv, gateCountsCsv, gateSquaresCsv, gateJson, gateApply: (j) => gateApply(j, "a test object"),
      gateSetColour, umapPickPh, gateRulesEdited,
    });
  }

  // ---- UI: overlays -----------------------------------------------------------------
  function buildOverlayUI() {
    const list = $("overlayList"); list.innerHTML = "";
    $("panelOverlays").hidden = S.overlays.length === 0;
    for (const o of S.overlays) {
      const row = document.createElement("label");
      row.className = "row check";
      row.innerHTML = `<input type="checkbox" ${o.on ? "checked" : ""}/><span>${o.label}</span>`;
      row.querySelector("input").onchange = async (e) => {
        o.on = e.target.checked;
        if (o.on && !o.data) {
          status(`Loading ${o.label}…`);
          try { o.data = await S.src.overlay(S.src.kind === "static" ? o.file : o.key); } finally { status(null); }
        }
        gateCountsRefresh();
        render();
      };
      list.appendChild(row);
    }
  }

  // ---- dataset ------------------------------------------------------------------------

  // ---- RNA-seq: one linked UMAP panel per gene --------------------------------------
  // All panels share one WebGL context (one Deck, several views) and one pan/zoom state,
  // and the number of gene panels is capped so large datasets cannot exhaust memory.
  const CAT_COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#bcbd22",
    "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5", "#c49c94", "#f7b6d2", "#dbdb8d", "#9edae5",
    "#7f7f7f", "#393b79"];
  const catColor = (i) => (i < CAT_COLORS.length ? hexToRgb255(CAT_COLORS[i]) : hslRgb((i * 137.508) % 360, 0.65, 0.55));
  function hslRgb(h, s, l) {
    const a = s * Math.min(l, 1 - l), f = (n) => { const k = (n + h / 30) % 12; return Math.round(255 * (l - a * Math.max(-1, Math.min(k - 3, 9 - k, 1)))); };
    return [f(0), f(8), f(4)];
  }
  // tab20 was drawn for white paper. On this background #1f77b4 (blue), #8c564b
  // (brown) and #393b79 (navy) sink into it, and a cluster you cannot see is one
  // whose label you cannot read either. Lift anything under the floor toward its
  // own lighter self, keeping the hue and most of the saturation so the legend
  // swatch, the points and the label stay recognisably the same color.
  const CAT_MIN_Y = cfg.catMinLuminance || 0.21;
  const LABEL_MIN_Y = cfg.labelMinLuminance || 0.5;
  const _lin = (c) => { c /= 255; return c <= 0.04045 ? c / 12.92 : Math.pow((c + 0.055) / 1.055, 2.4); };
  const _lum = ([r, g, b]) => 0.2126 * _lin(r) + 0.7152 * _lin(g) + 0.0722 * _lin(b);
  function liftTo(rgb, target) {
    // Relative luminance, not HSL lightness: by lightness pure orange scores 0.53 and
    // would get "brightened" into a duller, darker orange. Anything already bright
    // enough is left exactly alone; the rest is mixed toward white by the smallest
    // amount that clears the floor, which keeps the hue.
    if (_lum(rgb) >= target) return rgb;
    let lo = 0, hi = 0.95;
    for (let i = 0; i < 14; i++) {
      const t = (lo + hi) / 2;
      if (_lum(rgb.map((v) => v + (255 - v) * t)) >= target) hi = t; else lo = t;
    }
    return rgb.map((v) => Math.round(v + (255 - v) * hi));
  }
  const liftForDark = (rgb) => liftTo(rgb, CAT_MIN_Y);
  // Label text sits over its own cluster's points, so clearing the background is not
  // enough - it has to clear the color it is sitting on. A higher floor than the
  // points themselves, and the halo below does the rest.
  const labelColor = (rgb) => liftTo(rgb, LABEL_MIN_Y);


  // Viridis, 16 anchors interpolated to 256. FIVE Atlas offers five maps; one is
  // enough here, and viridis is the one to keep - perceptually uniform and readable
  // to color-blind viewers, which turbo and the rest are not.
  const RAMP = (stops) => {
    const lut = new Uint8Array(256 * 3);
    for (let i = 0; i < 256; i++) {
      const t = (i / 255) * (stops.length - 1), a = Math.floor(t), b = Math.min(stops.length - 1, a + 1), f = t - a;
      for (let c = 0; c < 3; c++) lut[i * 3 + c] = Math.round(stops[a][c] * (1 - f) + stops[b][c] * f);
    }
    return lut;
  };
  const PALETTES = {
    viridis: RAMP([[68,1,84],[72,26,108],[71,47,125],[65,68,135],[57,86,140],[49,104,142],[42,120,142],
      [35,136,142],[31,152,139],[34,168,132],[53,183,121],[84,197,104],[122,209,81],[165,219,54],[210,226,27],[253,231,37]]),
  };
  const VIRIDIS = PALETTES.viridis;
  const VIRIDIS_CSS = "linear-gradient(90deg, #440154, #3b528b, #21918c, #5ec962, #fde725)";
  // A gene in ink mode spans from its lightest printable tint - the 0.2 floor every
  // occupied square gets - to the full color. Show that span rather than one chip,
  // so the legend reads on the same scale as the map.
  const inkCss = (hex) => {
    const [r, g, b] = hexToRgb255(hex);
    const pale = [r, g, b].map((v) => Math.round((1 - 0.2 * (1 - v / 255)) * 255));
    return `linear-gradient(90deg, rgb(${pale.join(",")}), ${hex})`;
  };
  const rampCss = (name) => {
    const lut = PALETTES[name] || PALETTES.viridis;
    const stops = [];
    for (let i = 0; i <= 6; i++) {
      const t = Math.round((i / 6) * 255) * 3;
      stops.push(`rgb(${lut[t]},${lut[t + 1]},${lut[t + 2]})`);
    }
    return `linear-gradient(90deg, ${stops.join(", ")})`;
  };

  // The only clustering columns we publish. First match wins; a dataset carries one.
  const RNA_COLOR_BY = cfg.rnaColorBy || ["leiden"];

  // Put a name on every cluster, outside its own cloud of points and out of anyone
  // else's. Placement is solved once, in UMAP coordinates rather than screen pixels,
  // so a label keeps its bearing relative to its cluster at every zoom instead of
  // swimming around as you scroll.
  //
  // For each cluster: the centroid, then a ring of eight candidate positions just
  // beyond the cluster's own median radius. Each candidate is scored on how many
  // cells - of any cluster - sit under it, plus a penalty for crowding a label that
  // is already placed. Big clusters are placed first, so the crowded small ones
  // inherit the leftover room rather than the other way round.
  function placeClusterLabels(R, o, codes, lut) {
    const K = o.categories.length, n = R.n, pos = R.pos;
    const b = S.m.bounds, spanX = Math.max(b[2] - b[0], 1e-6), spanY = Math.max(b[3] - b[1], 1e-6);
    const span = Math.max(spanX, spanY);

    // occupancy grid over the whole embedding: one O(1) density lookup per candidate
    const GN = 160, gx = (x) => Math.min(GN - 1, Math.max(0, ((x - b[0]) / spanX * GN) | 0));
    const gy = (y) => Math.min(GN - 1, Math.max(0, ((y - b[1]) / spanY * GN) | 0));
    const grid = new Float32Array(GN * GN);
    for (let i = 0; i < n; i++) grid[gy(pos[i * 2 + 1]) * GN + gx(pos[i * 2])] += 1;

    let mx0 = 0, my0 = 0;
    for (let i = 0; i < n; i++) { mx0 += pos[i * 2]; my0 += pos[i * 2 + 1]; }
    mx0 /= Math.max(n, 1); my0 /= Math.max(n, 1);   // centre of mass of the whole plot

    const sum = new Float64Array(K * 2), cnt = new Float64Array(K);
    for (let i = 0; i < n; i++) {
      const k = codes[i]; if (k >= K) continue;
      sum[k * 2] += pos[i * 2]; sum[k * 2 + 1] += pos[i * 2 + 1]; cnt[k] += 1;
    }
    // median distance to the centroid: robust to the stragglers every cluster has
    const rad = new Float64Array(K), bucket = [];
    for (let k = 0; k < K; k++) bucket.push([]);
    const stride = Math.max(1, Math.floor(n / 60000));        // sampled: exact radius is not needed
    for (let i = 0; i < n; i += stride) {
      const k = codes[i]; if (k >= K || !cnt[k]) continue;
      const dx = pos[i * 2] - sum[k * 2] / cnt[k], dy = pos[i * 2 + 1] - sum[k * 2 + 1] / cnt[k];
      bucket[k].push(Math.hypot(dx, dy));
    }
    for (let k = 0; k < K; k++) {
      const a = bucket[k]; if (!a.length) { rad[k] = span * 0.02; continue; }
      a.sort((x, y) => x - y); rad[k] = Math.max(a[(a.length * 0.6) | 0], span * 0.01);
    }

    // east first, so a label lands to the right of its cluster whenever that is free
    const DIRS = [[1, 0], [0.71, 0.71], [0.71, -0.71], [0, 1], [0, -1],
                  [-0.71, 0.71], [-0.71, -0.71], [-1, 0]];
    const probe = span * 0.035;                               // radius the label text occupies
    // The text is drawn at a fixed pixel size, so its footprint in UMAP units depends
    // on the zoom. Guessing it as a fraction of `span` was out by half on a narrow
    // panel, which is why two labels could overlap while the solver thought they were
    // clear. Derive it from the same fit the renderer uses instead.
    const L0 = rnaLayout();
    const ppu = Math.max(1e-9, Math.min((L0.W / L0.cols) / spanX, (L0.H / L0.rows) / spanY) * 0.88);
    const labelBox = (name, x, y, anchor) => {
      const hw = (3.7 * name.length + 7) / ppu, hh = 9 / ppu;
      return { x: anchor === "start" ? x + hw : anchor === "end" ? x - hw : x, y, hw, hh };
    };
    const density = (x, y) => {
      let t = 0;
      const r = Math.max(1, Math.round(probe / spanX * GN)), cx = gx(x), cy = gy(y);
      for (let j = -r; j <= r; j++) for (let i2 = -r; i2 <= r; i2++) {
        const a = cx + i2, c = cy + j;
        if (a >= 0 && a < GN && c >= 0 && c < GN) t += grid[c * GN + a];
      }
      return t / ((2 * r + 1) * (2 * r + 1));
    };

    // Positions placed by hand in content/label_editor.html win outright. No solver
    // gets this as right as someone looking at the plot, and once a figure is settled
    // it should not move because the heuristic was retuned.
    const fixed = dsSettings(S.dsId).labelPositions || (cfg.labelPositions && cfg.labelPositions[S.dsId]) || {};
    const order = Array.from({ length: K }, (_, k) => k).filter((k) => cnt[k] > 0)
      .sort((a, c) => cnt[c] - cnt[a]);
    const placed = [], out = [];
    const STEPS = [0.15, 0.3, 0.45, 0.65, 0.9, 1.2, 1.6];
    // the fit leaves a small margin round the data, so a label pushed past the
    // bounds is a label nobody can read - keep every one inside the frame
    const mx = spanX * 0.04, my = spanY * 0.04;
    const clampX = (v) => Math.min(b[2] - mx, Math.max(b[0] + mx, v));
    const clampY = (v) => Math.min(b[3] - my, Math.max(b[1] + my, v));
    for (const k of order) {
      const name = o.categories[k];
      const hand = fixed[name];
      if (hand && hand.length === 2 && Number.isFinite(hand[0]) && Number.isFinite(hand[1])) {
        // the editor stores the centre of the label, so it anchors dead centre
        const box = labelBox(name, hand[0], hand[1], "middle");
        placed.push(box);
        out.push({
          name, x: hand[0], y: hand[1], count: cnt[k],
          color: [...labelColor(lut[k] || [200, 200, 200]), 255],
          anchor: "middle", baseline: "center",
        });
        continue;
      }
      const cx = sum[k * 2] / cnt[k], cy = sum[k * 2 + 1] / cnt[k];
      // Scale both scoring terms to 0..1 before weighing them against each other.
      // density() counts cells per grid square, which on 175k cells runs into the
      // hundreds, while the distance term was a fraction of the plot width - about
      // three. So distance counted for nothing and every label went wherever the
      // plot was emptiest, which is a long way from the cluster it names.
      const core = Math.max(density(cx, cy), 1);
      const reach = Math.max(rad[k] * 2, 1e-6);
      // Which way is "out" for this cluster: away from the middle of the plot. A
      // label on an outlying cluster belongs on its far side, in the empty space,
      // rather than tucked back toward the crowd - that is where there is room, and
      // it is the side a reader's eye is already on. Fixed east said nothing about
      // the plot and put Perineurial's label between it and everything else.
      const ox0 = cx - mx0, oy0 = cy - my0, om = Math.hypot(ox0, oy0) || 1;
      const outX = ox0 / om, outY = oy0 / om;
      let best = null;
      // Every direction at every distance is a candidate, not just the first gap
      // along each ray. With one candidate per direction a label boxed in by its
      // neighbours had nowhere to go and simply sat on top of one of them; now it
      // can step further out, and the distance term keeps it as close as it can be.
      DIRS.forEach(([ux, uy], di) => {
        for (const t of STEPS) {
          const d = rad[k] * t + probe * 0.5;
          const x = clampX(cx + ux * d), y = clampY(cy + uy * d);
          const anchor = ux > 0.3 ? "start" : ux < -0.3 ? "end" : "middle";
          const box = labelBox(o.categories[k], x, y, anchor);
          // clamping the anchor is not enough: the text grows away from it, so a
          // label anchored just inside the frame can still run off the edge
          const over = Math.max(0, (box.x + box.hw) - (b[2] - mx)) + Math.max(0, (b[0] + mx) - (box.x - box.hw))
                     + Math.max(0, (box.y + box.hh) - (b[3] - my)) + Math.max(0, (b[1] + my) - (box.y - box.hh));
          // Sit on the edge of the cloud: far enough out that the label is not lost
          // in the thick of the points, no further. Both terms are now fractions of
          // the same size, so the trade between them is the one intended.
          // Emptiness decides, the rest only breaks ties. These were close enough in
          // magnitude that "do not point back toward the middle" could outvote "do
          // not sit on top of the points", which is how a label ended up over its own
          // cluster with clear space just the other side of it.
          const inward = 1 - (ux * outX + uy * outY);   // 0 straight out, 2 straight in
          let sc = (density(x, y) / core) * 6 + (d / reach) * 3 + inward * 0.3
                 + (over / span) * 150;
          for (const q of placed) {
            const ox = (box.hw + q.hw) - Math.abs(box.x - q.x);
            const oy = (box.hh + q.hh) - Math.abs(box.y - q.y);
            if (ox > 0 && oy > 0) sc += Math.min(ox / (box.hw + q.hw), oy / (box.hh + q.hh)) * 200;
          }
          if (!best || sc < best.sc) best = { sc, x, y, ux, uy, box };
        }
      });
      placed.push(best.box);
      out.push({
        name, x: best.x, y: best.y, count: cnt[k],
        color: [...labelColor(lut[k] || [200, 200, 200]), 255],
        anchor: best.ux > 0.3 ? "start" : best.ux < -0.3 ? "end" : "middle",
        // world +y is up here (flipY is false), so a label above its cluster sits on
        // its own bottom edge
        baseline: best.uy > 0.3 ? "bottom" : best.uy < -0.3 ? "top" : "center",
      });
    }
    return out;
  }

  function leaveRna() {
    if (S.mode !== "rna") return;
    S.mode = "image"; S.rna = null; S.viewState = null;
    deckgl.setProps({ views: imageViews(), viewState: { target: [0, 0, 0], zoom: 0 }, layers: [] });
    $("umapLabels").hidden = true; $("umapLabels").innerHTML = "";
    $("densityKey").hidden = true;
    $("panelRna").hidden = true; $("panelChannels").hidden = false;
    $("hud").hidden = false; $("scalebar").hidden = false;
  }

  async function openRna() {
    S.mode = "rna";
    const m = S.m;
    for (const id of ["panelChannels", "panelGenes", "panelOverlays", "btnAllOn", "btnAllOff", "btnImage", "hud", "scalebar", "densityKey"]) $(id).hidden = true;
    $("panelRna").hidden = false; $("umapLabels").hidden = false;
    $("dsInfo").textContent = `${m.nCells.toLocaleString()} cells · ${m.nGenes.toLocaleString()} genes · ${m.valueLabel}`;
    status("Loading UMAP…");
    let pos, index;
    try {
      [pos, index] = await Promise.all([S.src.buf("umap.bin").then((b) => new Float32Array(b)),
        S.src.buf("genes.json").then((b) => JSON.parse(new TextDecoder().decode(b).replace(/:(NaN|-?Infinity)\b/g, ":null")))]);
    } catch (e) { status(`Cannot load UMAP: ${e.message}`, true); return; }
    const n = pos.length / 2;
    S.rna = { pos, n, index, genes: [], obs: m.obs || [], colorBy: null, obsData: null, grid: "", vs: null,
              labels: null, labelsKey: null, cat: null, showLabels: $("rnaLabels").checked,
              view: "normal", palette: "viridis", bin: 1, heat: null,
              pointSize: n > 200000 ? 1 : n > 50000 ? 1.5 : 2.5,
              baseData: { length: n, attributes: { getPosition: { value: pos, size: 2 } } } };
    $("rnaPointSize").value = S.rna.pointSize;
    $("rnaView").value = S.rna.view; $("rnaPalette").value = S.rna.palette;
    $("rnaBin").value = String(S.rna.bin); $("rowRnaBin").hidden = S.rna.view !== "density";
    // One clustering per dataset, so there is nothing to choose between and no picker.
    // RNA_COLOR_BY is tried in order; the fallback only matters if a future dataset
    // arrives without either column.
    const cats = S.rna.obs.filter((o) => o.type === "category" && !/sample|doublet/i.test(o.name));
    const pick = RNA_COLOR_BY.map((want) => cats.find((o) => o.name === want)).find(Boolean)
      || cats.find((o) => /cluster|leiden|louvain/i.test(o.name)) || cats[0] || null;
    await setColorBy(pick ? pick.name : "");
    buildRnaGeneUI();
    status(null);
    if (window.SodaPlots) window.SodaPlots.reset();
  }

  async function setColorBy(name) {
    const R = S.rna; if (!R) return;
    R.colorBy = name || null; R.obsData = null; R.labels = null; R.labelsKey = null; R.cat = null;
    const o = R.obs.find((x) => x.name === name);
    const leg = $("rnaLegend"); leg.innerHTML = "";
    if (o) {
      const buf = await S.src.buf(`obs/${o.file}`);
      const colors = new Uint8Array(R.n * 4);
      if (o.type === "category") {
        const codes = new Uint16Array(buf);
        const lut = o.categories.map((_, i) => liftForDark(o.colors && o.colors[i] ? hexToRgb255(o.colors[i]) : catColor(i)));
        for (let i = 0; i < R.n; i++) { const c = lut[codes[i]] || [128, 128, 128]; const k = i * 4; colors[k] = c[0]; colors[k + 1] = c[1]; colors[k + 2] = c[2]; colors[k + 3] = 255; }
        leg.innerHTML = o.categories.map((c, i) => `<div class="leg"><i style="background:rgb(${lut[i].join(",")})"></i><span>${c}</span><small>${(o.counts[i] || 0).toLocaleString()}</small></div>`).join("");
        // Placement depends on how big a label is in UMAP units, which depends on
        // the fitted zoom - and at this point the panel has not been laid out yet
        // (clientWidth is still 0 on a first open). Keep the inputs and let
        // renderRna solve it once it knows the real geometry.
        R.cat = { o, codes, lut };
      } else {
        const vals = new Float32Array(buf), lo = o.min, span = Math.max(o.max - o.min, 1e-6);
        for (let i = 0; i < R.n; i++) {
          const t = Math.max(0, Math.min(255, Math.round(((vals[i] - lo) / span) * 255))), k = i * 4;
          colors[k] = VIRIDIS[t * 3]; colors[k + 1] = VIRIDIS[t * 3 + 1]; colors[k + 2] = VIRIDIS[t * 3 + 2]; colors[k + 3] = Number.isFinite(vals[i]) ? 255 : 0;
        }
        leg.innerHTML = `<div class="cbar"><span>${o.min.toFixed(2)}</span><i style="background:${VIRIDIS_CSS}"></i><span>${o.max.toFixed(2)}</span></div>`;
      }
      R.obsData = { length: R.n, attributes: { getPosition: { value: R.pos, size: 2 }, getFillColor: { value: colors, size: 4, normalized: true } } };
    }
    renderRna();
  }

  function showRnaResults(q) {
    const box = $("rnaResults"), R = S.rna;
    q = q.trim().toLowerCase();
    if (!q || !R) { box.hidden = true; return; }
    const hits = [];
    const exact = R.index.find((g) => g.n.toLowerCase() === q);          // an exact name always comes first
    if (exact) hits.push(exact);
    for (const g of R.index) if (g !== exact && g.n.toLowerCase().startsWith(q)) { hits.push(g); if (hits.length >= 40) break; }
    if (hits.length < 40) for (const g of R.index) { const n = g.n.toLowerCase(); if (!n.startsWith(q) && n.includes(q)) { hits.push(g); if (hits.length >= 40) break; } }
    box.innerHTML = "";
    hits.forEach((g, i) => {
      const d = document.createElement("div");
      d.innerHTML = `<span>${g.n}</span><small>${(g.p * 100).toFixed(1)}% of cells</small>`;
      if (i === 0) d.classList.add("hot");
      d.onmousedown = (e) => { e.preventDefault(); addRnaGene(g.n); };
      box.appendChild(d);
    });
    box.hidden = hits.length === 0;
  }

  async function addRnaGene(name) {
    const R = S.rna; if (!R) return;
    $("rnaResults").hidden = true; $("rnaSearch").value = "";
    if (R.genes.some((g) => g.name === name)) return;
    if (R.genes.length >= MAX_RNA_GENES) {
      status(`Up to ${MAX_RNA_GENES} genes at once. Remove one to add ${name}.`, true);
      setTimeout(() => status(null), 4000);
      return;
    }
    status(`Loading ${name}…`);
    let buf;
    try { buf = await S.src.buf(`gene/${encodeURIComponent(name)}.bin`); }
    catch (e) { status(`${name}: not expressed in this dataset`, true); setTimeout(() => status(null), 4000); return; }
    const k = new Uint32Array(buf, 0, 1)[0];
    const idx = new Uint32Array(buf, 4, k), val = new Float32Array(buf, 4 + 4 * k, k);
    const order = Array.from({ length: k }, (_, i) => i).sort((a, b) => val[a] - val[b]);   // brightest drawn last
    const pos = new Float32Array(k * 2), v = new Float32Array(k);
    order.forEach((o, j) => { pos[2 * j] = R.pos[2 * idx[o]]; pos[2 * j + 1] = R.pos[2 * idx[o] + 1]; v[j] = val[o]; });
    // default display range = 5th..99th percentile of the expressing cells (quantile cutoffs, as Seurat FeaturePlot)
    const q = (f) => (k ? v[Math.min(k - 1, Math.floor(f * k))] : 0);
    const vmax = k ? v[k - 1] : 1, p99 = q(0.99) || vmax;
    let p05 = q(0.05);
    if (p05 >= p99) p05 = 0;
    const g = { name, k, pos, val: v, min: +p05.toFixed(2), max: +p99.toFixed(2), vmax: Math.max(vmax, 0.01),
                frac: k / R.n, hex: GENE_COLORS[R.genes.length % GENE_COLORS.length], visible: true };
    recolorRnaGene(g);
    if (!S.rna || S.rna !== R) return;                 // dataset switched while loading
    R.genes.push(g); R.heat = null; R.grid = "";
    buildRnaGeneUI();
    renderRna();
    status(null);
  }

  function recolorRnaGene(g) {
    const R = S.rna;
    const colors = new Uint8Array(g.k * 4), lo = g.min, span = Math.max(g.max - lo, 1e-6);
    const useInk = R && R.palette === "genes";
    const lut = (R && PALETTES[R.palette]) || VIRIDIS;
    const [hr, hg, hb] = hexToRgb255(g.hex || "#33ddff");
    for (let i = 0; i < g.k; i++) {
      const x = g.val[i], k = i * 4;
      if (x < lo || x <= 0) { colors[k + 3] = 0; continue; }
      const f = Math.min(1, (x - lo) / span);
      if (useInk) {
        // the gene's own color, faint where it is barely expressed
        colors[k] = hr; colors[k + 1] = hg; colors[k + 2] = hb; colors[k + 3] = 40 + Math.round(215 * f);
      } else {
        const t = Math.min(255, Math.round(f * 255)) * 3;
        colors[k] = lut[t]; colors[k + 1] = lut[t + 1]; colors[k + 2] = lut[t + 2]; colors[k + 3] = 255;
      }
    }
    g.data = { length: g.k, attributes: { getPosition: { value: g.pos, size: 2 }, getFillColor: { value: colors, size: 4, normalized: true } } };
  }

  function buildRnaGeneUI() {
    const R = S.rna; if (!R) return;
    $("rnaCount").textContent = `(${R.genes.length} / ${MAX_RNA_GENES})`;
    const list = $("rnaGeneList"); list.innerHTML = "";
    R.genes.forEach((g) => {
      const row = document.createElement("div");
      row.className = "chan";
      const sw = R.palette === "genes" ? inkCss(g.hex) : rampCss(R.palette);
      row.innerHTML = `<input class="cb" type="checkbox" ${g.visible === false ? "" : "checked"} /><span class="col sw" style="background:${sw}"></span>
        <span class="name" title="${g.name}">${g.name}</span><button class="rm x" title="remove">×</button><div class="range"></div>`;
      row.querySelector(".x").onclick = () => { R.genes = R.genes.filter((x) => x !== g); buildRnaGeneUI(); renderRna(); };
      const step = g.vmax > 20 ? 0.5 : 0.05;
      rangeSlider(row.querySelector(".range"), {
        min: 0, max: +g.vmax.toFixed(2), step, lo: g.min, hi: Math.min(g.max, +g.vmax.toFixed(2)),
        color: () => "#5ec962", fmt: (x) => x.toFixed(2),
        onChange: (a, b) => {
          g.min = a; g.max = b; recolorRnaGene(g);
          const hk = `${g.hex}|${a}`;
          if (g.heatKey !== hk) { g.heatKey = hk; R.heat = null; }
          renderRna();
        },
      });
      const cb = row.querySelector(".cb");
      cb.title = `${g.k.toLocaleString()} cells expressing`;
      cb.onchange = () => {
        g.visible = cb.checked;
        row.classList.toggle("off", !cb.checked);
        // a hidden gene loses its panel in Normal view, so the grid is re-fitted
        R.heat = null; R.grid = "";
        renderRna();
      };
      row.classList.toggle("off", g.visible === false);
      list.appendChild(row);
    });
  }

  // The UMAP's density view: the same square-bin model as the slide, in embedding
  // space. Density collapses the small multiples into ONE panel on purpose - the
  // point of it is that genes combine, either as inks that multiply or as one ramp
  // over their total, and neither says anything if each gene keeps its own panel.
  function rnaHeat() {
    const R = S.rna;
    const b = S.m.bounds, span = Math.max(b[2] - b[0], b[3] - b[1], 1e-6);
    const cell = (span / 200) * Math.max(1, R.bin);
    const ink = R.palette === "genes";
    const on = R.genes.filter((g) => g.visible !== false);
    const key = on.map((g) => `${g.name}:${g.hex}:${g.min}`).join("|") + `#${R.palette}#${R.bin}`;
    if (R.heat && R.heat.key === key) return R.heat;

    const idx = new Map(), KX = 1 << 21;
    const binOf = (x, y) => (Math.floor((y - b[1]) / cell)) * KX + Math.floor((x - b[0]) / cell);
    for (const g of on) for (let i = 0; i < g.k; i++) {
      const k = binOf(g.pos[2 * i], g.pos[2 * i + 1]);
      if (!idx.has(k)) idx.set(k, idx.size);
    }
    const n = idx.size;
    const positions = new Float32Array(n * 3), colors = new Uint8Array(n * 4);
    for (const [k, j] of idx) {
      const iy = Math.floor(k / KX), ix = k - iy * KX;
      positions[j * 3] = b[0] + (ix + 0.5) * cell; positions[j * 3 + 1] = b[1] + (iy + 0.5) * cell;
    }
    // expression summed over the cells in each square, scaled by how much the busiest
    // squares hold - the same correction the slide view needed
    const norm = on.map((g) => {
      const sum = new Float32Array(n), lo = g.min;
      for (let i = 0; i < g.k; i++) {
        const x = g.val[i];
        if (x >= lo && x > 0) sum[idx.get(binOf(g.pos[2 * i], g.pos[2 * i + 1]))] += x;
      }
      const nz = Array.from(sum).filter((v) => v > 0).sort((x, y) => x - y);
      const hi = nz.length ? Math.max(nz[Math.floor((nz.length - 1) * 0.99)], 1e-6) : 1;
      for (let j = 0; j < n; j++) sum[j] = Math.min(1, sum[j] / hi);
      return sum;
    });

    if (ink) {
      const paper = new Float32Array(n * 3).fill(1), miss = new Float32Array(n).fill(1);
      on.forEach((g, gi) => {
        const d = norm[gi], [tr, tg, tb] = hexToRgb255(g.hex);
        for (let j = 0; j < n; j++) {
          if (d[j] <= 0) continue;
          const sEff = 0.2 + 0.8 * d[j];
          paper[j * 3] *= 1 - sEff * (1 - tr / 255);
          paper[j * 3 + 1] *= 1 - sEff * (1 - tg / 255);
          paper[j * 3 + 2] *= 1 - sEff * (1 - tb / 255);
          miss[j] *= 1 - sEff;
        }
      });
      for (let j = 0; j < n; j++) {
        const c = j * 4;
        colors[c] = Math.round(paper[j * 3] * 255); colors[c + 1] = Math.round(paper[j * 3 + 1] * 255);
        colors[c + 2] = Math.round(paper[j * 3 + 2] * 255);
        colors[c + 3] = Math.round(Math.sqrt(Math.max(0, 1 - miss[j])) * 255);
      }
    } else {
      const lut = PALETTES[R.palette] || VIRIDIS;
      const inv = 1 / Math.max(1, norm.length);      // average, not sum - see geneHeat
      for (let j = 0; j < n; j++) {
        let total = 0, any = false;
        for (const d of norm) { total += d[j]; if (d[j] > 0) any = true; }
        const t = Math.min(255, Math.round((0.06 + 0.94 * Math.min(1, total * inv)) * 255)) * 3;
        const c = j * 4;
        colors[c] = lut[t]; colors[c + 1] = lut[t + 1]; colors[c + 2] = lut[t + 2];
        colors[c + 3] = any ? 235 : 0;
      }
    }
    R.heat = { key, n, positions, colors, cell, index: Array.from({ length: n }, (_, j) => j) };
    return R.heat;
  }


  // A gene switched off keeps its row and its window, but loses its panel and
  // stops contributing to the combined density.
  const rnaShownGenes = () => (S.rna ? S.rna.genes.filter((g) => g.visible !== false) : []);

  function rnaLayout() {
    const R = S.rna, el = $("viewport");
    const shown = rnaShownGenes();
    const dens = R.view === "density" && shown.length > 0;
    const n = dens ? 2 : 1 + shown.length;
    const W = el.clientWidth, H = el.clientHeight;
    let best = { cols: 1, rows: n, size: 0 };
    for (let c = 1; c <= n; c++) { const r = Math.ceil(n / c), size = Math.min(W / c, H / r); if (size > best.size) best = { cols: c, rows: r, size }; }
    return { n, W, H, ...best };
  }

  function fitRna() {
    const R = S.rna; if (!R) return;
    R.vs = null;                                              // renderRna refits when there is no view state
    renderRna();
  }

  function rnaViewStates() {
    const R = S.rna, out = {};
    for (let i = 0; i <= R.genes.length; i++) out[`umap-${i}`] = R.vs;
    return out;
  }

  function renderRna() {
    const R = S.rna; if (!R || !deckgl) return;
    const L = rnaLayout();
    if (!L.W || !L.H) { renderPending = true; return; }
    const gridKey = `${L.cols}x${L.rows}@${L.W}x${L.H}`;   // L.n already follows what is shown
    if (!R.vs || R.grid !== gridKey) {                        // first draw, or panel size changed: refit
      const b = S.m.bounds, pw = L.W / L.cols, ph = L.H / L.rows;
      const zoom = Math.log2(Math.min(pw / Math.max(b[2] - b[0], 1e-6), ph / Math.max(b[3] - b[1], 1e-6)) * 0.88);
      R.vs = { target: [(b[0] + b[2]) / 2, (b[1] + b[3]) / 2, 0], zoom, minZoom: zoom - 3, maxZoom: zoom + 8 };
    }
    R.grid = gridKey;
    // Re-solved whenever the fit changes - a resize changes how much of the UMAP a
    // label covers, and so which placement is actually free.
    const key = `${R.colorBy}|${gridKey}`;
    if (R.cat && R.labelsKey !== key) {
      R.labels = placeClusterLabels(R, R.cat.o, R.cat.codes, R.cat.lut);
      R.labelsKey = key;
    }
    const ctl = { scrollZoom: { speed: 0.01, smooth: false }, inertia: false, doubleClickZoom: true, keyboard: false };
    const views = [], labels = [];
    const shown = rnaShownGenes();
    const dens = R.view === "density" && shown.length > 0;
    for (let i = 0; i < L.n; i++) {
      const c = i % L.cols, r = Math.floor(i / L.cols);
      const box = { x: `${(c * 100) / L.cols}%`, y: `${(r * 100) / L.rows}%`, width: `${100 / L.cols}%`, height: `${100 / L.rows}%` };
      views.push(new OrthographicView({ id: `umap-${i}`, ...box, flipY: false, controller: ctl }));
      const style = `left:${box.x};top:${box.y};width:${box.width};height:${box.height}`;
      if (i === 0) {
        labels.push(`<div class="umap-cell" style="${style}"><div class="t">${R.colorBy ? `Color: ${R.colorBy}` : "All cells"}</div></div>`);
      } else if (dens) {
        const names = shown.map((g) => g.name).join(" + ");
        // say what the colors mean, same as the slide view's key
        const key = R.palette === "genes"
          ? `<div class="genes">${shown.map((g) => `<span class="g"><i style="background:${inkCss(g.hex)}"></i>${g.name}</span>`).join("")}</div>` +
            `<div class="note">pale = few cells expressing, saturated = many` +
            (shown.length > 1 ? ". Where genes overlap the colors mix and darken." : "") + `</div>`
          : `<div class="ramp"><span>low</span><i style="background:${rampCss(R.palette)}"></i><span>high</span></div>` +
            `<div class="note">expression per square` +
            (shown.length > 1 ? `, averaged over ${shown.length} genes` : ``) + `</div>`;
        labels.push(`<div class="umap-cell" style="${style}"><div class="t">${names} <small>combined density</small></div>` +
          `<div class="dkey inpanel">${key}</div></div>`);
      } else {
        const g = shown[i - 1];
        const sw = R.palette === "genes" ? inkCss(g.hex) : rampCss(R.palette);
        labels.push(`<div class="umap-cell" style="${style}"><div class="t">${g.name} <small>${(g.frac * 100).toFixed(1)}% expressing</small></div>` +
          `<div class="cbar"><span>${g.min.toFixed(2)}</span><i style="background:${sw}"></i><span>${g.max.toFixed(2)}</span></div></div>`);
      }
    }
    $("umapLabels").innerHTML = labels.join("");
    const common = { radiusUnits: "pixels", getRadius: R.pointSize, radiusMinPixels: 0.5, stroked: false, filled: true, parameters: NORMAL, antialiasing: true };
    const layers = [new ScatterplotLayer({ id: "umap-base", data: R.baseData, getFillColor: [80, 86, 100, 255], ...common,
      updateTriggers: { getRadius: [R.pointSize] } })];
    if (R.obsData) layers.push(new ScatterplotLayer({ id: `umap-0-obs-${R.colorBy}`, data: R.obsData, ...common }));
    if (R.obsData && R.showLabels && R.labels && R.labels.length) {
      layers.push(new TextLayer({
        id: `umap-0-labels-${R.colorBy}`, data: R.labels,
        getPosition: (d) => [d.x, d.y], getText: (d) => d.name, getColor: (d) => d.color,
        getTextAnchor: (d) => d.anchor, getAlignmentBaseline: (d) => d.baseline,
        sizeUnits: "pixels", getSize: 14, sizeMinPixels: 12, sizeMaxPixels: 15,
        characterSet: "auto", fontFamily: "system-ui, sans-serif", fontWeight: 700,
        // An SDF outline was tried here and made the glyphs thinner without reading as
        // a halo at this size. A near-opaque plate with generous padding does the job:
        // it is the plate, not the letterform, that lifts the label off the cloud.
        background: true, getBackgroundColor: [8, 11, 18, 232],
        backgroundPadding: [7, 4, 7, 4],
        getBorderColor: (d) => [...d.color.slice(0, 3), 110], getBorderWidth: 1,
        parameters: NORMAL, pickable: false,
        updateTriggers: { getPosition: [R.colorBy], getText: [R.colorBy] },
      }));
    }
    if (dens) {
      const h = rnaHeat();
      if (h.n) layers.push(new ColumnLayer({
        id: `umap-1-heat-${R.palette}`, data: h.index,
        getPosition: (j) => [h.positions[j * 3], h.positions[j * 3 + 1], 0],
        getFillColor: (j) => [h.colors[j * 4], h.colors[j * 4 + 1], h.colors[j * 4 + 2], h.colors[j * 4 + 3]],
        diskResolution: 4, angle: 45, radius: (h.cell * 1.02) / Math.SQRT2,
        radiusUnits: "common", extruded: false, filled: true, stroked: false,
        parameters: NORMAL, pickable: false,
        updateTriggers: { getFillColor: [h.key], getPosition: [h.key] },
      }));
    } else {
      shown.forEach((g, i) => layers.push(new ScatterplotLayer({ id: `umap-${i + 1}-g-${g.name}`, data: g.data, ...common })));
    }
    deckgl.setProps({ views, viewState: rnaViewStates(), layers });
  }

  let dsPicker = null;                             // dspicker.js: the searchable dataset list
  // layout.js (the tabs) is told when a dataset starts and stops opening: every panel shows or hides on the way,
  // and the tab on show is only re-chosen once they have settled. An open that a later one overtakes never
  // returns (see liveSource), so each start gets a number and only the latest one's end counts.
  async function openDataset(id) {
    const lay = window.SodaLayout, tok = lay ? lay.opening(true) : 0;
    try { return await openDatasetNow(id); } finally { if (lay) lay.opening(false, tok); }
  }
  async function openDatasetNow(id) {
    if (dsPicker) dsPicker.setValue(id);           // whoever asked (the picker, ?ds=, the console): the picker shows it
    // This open is now the only one that counts (see liveSource): the one before it stops at
    // its next await, its downloads are cancelled, and so is a layer switch still under way.
    const seq = ++openSeq, stale = () => seq !== openSeq;
    layerSeq++;
    if (openCtl) openCtl.abort();
    openCtl = new AbortController();
    status("Opening dataset…");
    S.dsId = id;
    tileCache.clear();
    tileProblems = 0; drawProblems = 0; tileRetries = 0; tilesFailed = false;
    if (stickyErr && stickyErr === tileMsg) stickyErr = null;
    S.channels = []; S.overlays = []; S.m = null;
    deckgl.setProps({ layers: [] });           // drop the previous dataset's layers right away
    S.genes = []; S.tissue.positions = null; S.heat = null;
    S.region = { meta: null, ids: null, lookup: null, selected: null, overlayKey: null };
    lassoReset();
    molReset();
    if (binQc) binQc.reset();
    umapReset();
    S.gate = gateBlank(); $("panelGate").hidden = true; $("rowGateOpen").hidden = true; $("gatePub").hidden = true;
    leaveRna();
    try {
      const d = S.datasets.find((x) => x.id === id);
      const staticRoot = params.get("static");   // ?static=http://host/folder -> <folder>/<id>/ (test exports before upload)
      let src;
      if (staticRoot) src = staticSource(`${staticRoot.replace(/\/$/, "")}/${id}/`);
      else if (TOKEN_URL && !params.get("api") && !(d && (d.linkPath || d.base))) {   // ?api= wins: a local tile server
        const t = await ensureToken();           // asks for the password once per session
        src = staticSource(`${t.base}${id}/`, tokenQuery);
      }
      else if (d && d.linkPath) src = staticSource(await resolveLink(d.linkPath));
      else if (d && d.base) src = staticSource(d.base, d.qs ? () => d.qs : undefined);   // d.qs: a ?data= link's read pass
      else src = apiSource(id);
      if (stale()) return;                       // another dataset was picked while the password box was up
      S.src = liveSource(src, seq, openCtl.signal);
      S.m = await S.src.manifest();
    } catch (e) { if (!stale()) status(`Cannot open ${id}: ${e.message}`, true); return; }
    if (stale()) return;
    if (S.m.kind === "rnaseq") {
      await openRna();
      if (stale()) return;
      const u = new URL(location.href); u.searchParams.set("ds", id); history.replaceState(null, "", u);
      return;
    }
    const img = S.m.image;
    const dsCfg = dsSettings(S.dsId);
    const wanted = (dsCfg.defaultChannels || cfg.defaultChannels || []).map((s) => String(s).toLowerCase());
    // "open-atto550", "open AF488", "open-cy5": detector slots with no antibody in
    // them. They carry autofluorescence and nothing a viewer should reason about,
    // and on the Akoya panels they are the three biggest files in the set.
    // Panel order out of the scanner is acquisition order, which is nothing to the
    // reader: finding CD8 among thirty markers meant reading the whole list. Sort
    // them, with the nuclear counterstain pinned at the top because it is the one
    // channel people leave on as a reference for all the others. `index` keeps
    // pointing at the right tiles, so only the row order moves.
    // Left alone on an RGB photo, where the channels are not markers at all.
    // Per-dataset exclusions, by exact name, on top of the blanket `open-*` rule.
    // Some panels carry an antibody that was only ever there to test the stain; it is
    // real data but not part of the panel, and leaving it in the list invites someone
    // to read it as one more marker.
    const dropped = new Set((dsCfg.hideChannels || []).map((x) => String(x).trim().toLowerCase()));
    const named = img.channels.filter((c) =>
      !HIDE_CHANNEL.test(c.name) && !dropped.has(c.name.trim().toLowerCase()));
    const photo = named.length === 3 && named.every((c) => /^(red|green|blue)$/i.test(c.name));
    if (!photo && named.length > 4) {
      const nuclear = (c) => (NUCLEAR_CHANNEL.test(c.name) ? 0 : 1);
      named.sort((a, b) => nuclear(a) - nuclear(b) ||
        a.name.localeCompare(b.name, undefined, { numeric: true, sensitivity: "base" }));
    }
    S.channels = named.map((c, i) => {
      const on = wanted.length ? wanted.some((w) => c.name.toLowerCase().includes(w)) : i < 3;
      return { index: c.index, name: c.name, hex: c.color, color: hexToRgb01(c.color),
               lo: Math.round(c.window[0]), hi: Math.round(c.window[1]), visible: on };
    });
    if (!S.channels.some((c) => c.visible)) S.channels.slice(0, 3).forEach((c) => (c.visible = true));
    buildChannelUI();
    const um = S.m.pixelSizeUm;
    $("dsInfo").textContent = `${img.width.toLocaleString()} × ${img.height.toLocaleString()} px · ${img.levels} levels · ${img.dtype}` +
      (um ? ` · ${um.toFixed(4)} µm/px` : "");

    const ex = S.m.expr || [];
    // Panels follow the platform: sequencing datasets (Visium) get transcripts and just an
    // image opacity control for the H&E; imaging datasets (Akoya, Xenium stains, TIFF) get
    // the per-channel marker list and no transcript panel.
    // An image dataset can also carry protein squares (Akoya grid): it keeps its marker
    // list and gets the expression panel too, worded for markers (see syncWording).
    // Xenium, CosMx and MERSCOPE are imaging platforms, and a GeoMx scan is an image too: their
    // stains stay a marker list even with a cells / ROI layer.
    const isVisium = S.m.kind === "visium" || (!["xenium", "cosmx", "merscope", "geomx"].includes(S.m.kind) && ex.some((e) => e.kind !== "protein"));
    const hasProt = !isVisium && ex.length > 0;
    $("panelGenes").hidden = !isVisium && !hasProt;
    $("chanList").hidden = isVisium;
    $("chanCount").hidden = isVisium;
    $("btnAllOn").hidden = isVisium;
    $("btnAllOff").hidden = isVisium;
    // Only where there is something underneath to see. On a Visium slide the image
    // is an H&E the transcripts are drawn over, so hiding it is useful; on a marker
    // panel the image IS the data, and the button only ever blanked the view.
    $("btnImage").hidden = !isVisium;
    if (!isVisium) S.imgOn = true;
    $("chanTitle").textContent = isVisium ? "Image" : "Markers / stains";
    if (isVisium) S.channels.forEach((c) => (c.visible = true));
    if (!ex.length) syncWording();             // no layer to switch to: back to the default words
    // Only the resolutions worth offering. The 2 um grid stays in the exports but is
    // not listed: at that size a bin holds a single transcript nearly everywhere, so
    // a density map of it is a presence map. The counts that give the view meaning
    // only appear once a bin is 8 um or larger.
    const offered = ex.filter((e) => !e.binSizeUm || e.binSizeUm >= EXPR_MIN_UM);
    const shownEx = offered.length ? offered : ex;
    const sel = $("exprSelect"); sel.innerHTML = "";
    for (const e of shownEx) {
      const o = document.createElement("option"); o.value = e.key;
      o.textContent = e.kind === "protein"
        ? `${e.label} · ${e.nBins.toLocaleString()} squares · ${e.nGenes.toLocaleString()} markers`
        : `${e.label} · ${e.nBins.toLocaleString()} ${e.kind === "cells" ? "cells" : e.kind === "rois" ? "ROI segments" : "bins"} · ${e.nGenes.toLocaleString()} genes`;
      sel.appendChild(o);
    }
    const isPhoto = img.channels.length === 3 && img.channels.every((c) => /^(red|green|blue)$/i.test(c.name));
    S.geneBlend = "heatmap";        // density: glow conflates intensity with overlap
    S.genePalette = "genes"; S.geneBin = 1; S.heat = null;
    $("geneBlend").value = S.geneBlend;
    $("genePalette").value = S.genePalette;
    // marker squares start a little see-through, so the stains under them always show
    S.imgOn = true; S.imgOpacity = 1; S.geneOpacity = hasProt ? 0.85 : 1;
    $("imgOpacity").value = 1; $("geneOpacity").value = S.geneOpacity;
    syncImageButton();
    S.overlays = (S.m.overlays || []).map((o) => ({ ...o, on: false, data: null }));
    buildOverlayUI();
    fitView();
    render();
    status(null);
    await molOnDataset();
    if (stale()) return;
    if (shownEx.length) { sel.value = shownEx[0].key; await setExprKey(shownEx[0].key); }
    if (stale()) return;
    // Tier A: turn the region overlay on so it's clickable without hunting for the toggle.
    if (S.region.overlayKey) {
      const ov = S.overlays.find((o) => o.file === S.region.overlayKey || o.key === S.region.overlayKey);
      if (ov && !ov.on) {
        ov.on = true;
        try { ov.data = await S.src.overlay(S.src.kind === "static" ? ov.file : ov.key); } catch (e) { /* leave off */ }
        if (stale()) return;
        buildOverlayUI(); render();
      }
    }
    const u = new URL(location.href); u.searchParams.set("ds", id); history.replaceState(null, "", u);
  }

  async function init() {
    sheetEvents();                                  // first, so the controls open even if loading fails
    initDeck();
    binQcInit();
    lassoInit();
    umapInit();
    gateInit();
    touchInit();
    const cfgSets = cfg.datasets || {};
    S.datasets = Object.keys(cfgSets).map((id) => ({ id, label: cfgSets[id].label || id, linkPath: cfgSets[id].linkPath,
      group: cfgSets[id].group, base: cfgSets[id].base, tags: cfgSets[id].tags }));
    // ?static=<root>&ds=<id> can open an exported dataset that is not in config.js yet
    if (params.get("static") && params.get("ds") && !S.datasets.some((d) => d.id === params.get("ds")))
      S.datasets.push({ id: params.get("ds"), label: `${params.get("ds")} (local export)` });
    let dataTitle = "", dataByline = "", dataPaper = "", dataDefault = "";
    if (DATA && DATA.error) { status(`Cannot open the data link: ${DATA.error}.`, true); return; }
    if (DATA) {                                     // ?data=: that folder's datasets replace config.js's
      const at = (id) => ({ base: `${DATA.root}${encodeURIComponent(id)}/`, qs: DATA.qs });
      const text = (v, max) => (typeof v === "string" ? v.trim().slice(0, max) : "");
      try {
        const cat = await getJSON(`${DATA.root}datasets.json${DATA.qs}`, true);
        S.datasets = (Array.isArray(cat.datasets) ? cat.datasets : []).filter((d) => d && d.id).map((d) => ({
          id: String(d.id), label: String(d.label || d.id), group: d.group ? String(d.group) : undefined,
          tags: folderTags(d.tags), settings: folderSettings(d), ...at(String(d.id)) }));
        S.dataGroups = {};
        for (const d of S.datasets) if (d.group) S.dataGroups[d.group] = d.group;   // first-seen order
        dataTitle = text(cat.title, 120);
        dataByline = text(cat.byline, 300);
        const paper = text(cat.paper, 2000);
        dataPaper = /^https:\/\//i.test(paper) ? paper : "";         // a link that opens a page, nothing else
        dataDefault = text(cat.default, 200);
        if (!S.datasets.length) { status(`The folder's datasets.json lists no datasets: ${DATA.root}`, true); return; }
      } catch (e) {
        const ds = params.get("ds");
        // the link may name one dataset's own folder (the one holding manifest.json), not the folder above
        const one = /^404 /.test(e.message) && !ds
          ? await getJSON(`${DATA.root}manifest.json${DATA.qs}`, true).catch(() => null) : null;
        if (one) {
          const id = String(one.id || "dataset");
          S.datasets = [{ id, label: String(one.label || id), base: DATA.root, qs: DATA.qs }];
        }
        else if (ds && /^404 /.test(e.message)) S.datasets = [{ id: ds, label: ds, ...at(ds) }];   // no list, one named set
        else {
          status(/^404 /.test(e.message)
            ? `No datasets.json in ${DATA.root} — run catalog.py on the export folder (the notebook does this) and upload it, or check the link.`
            : /^40[13] /.test(e.message)
            ? `${DATA.root} is not readable (${e.message.slice(0, 3)}). The folder must allow public reading, or the link must carry a read pass.`
            : `Could not reach ${DATA.root}. Check the link. If it is right, the storage must allow this page to read it (a CORS rule, see docs/HOSTING.md).`, true);
          return;
        }
      }
    }
    const useServer = !DATA && !!(params.get("api") || cfg.apiBase || cfg.useLocalServer);
    $("apiInfo").textContent = DATA ? `data: ${new URL(DATA.root).host}`
      : useServer ? `server: ${API || location.origin}` : "data: static links (blob)";
    if (useServer) {
      try {
        const fromServer = await getJSON(apiUrl("datasets"));
        for (const d of fromServer) if (!S.datasets.some((x) => x.id === d.id)) S.datasets.push(d);
      } catch (e) {
        if (!S.datasets.length) {
          status(/^401/.test(e.message) ? "Password not accepted." :
            `No tile server at ${API || location.origin} — run run.ps1 or add ?api=http://host:8760`, true);
          return;
        }
      }
    }
    // Grouped by what was measured. The channel count used to be appended here, but
    // it counted raw channels including the `open-*` slots the viewer hides, so it
    // disagreed with both the label and the stain list; the labels carry it now.
    // The list itself is dspicker.js: search, group headings (cfg.groups in their order, then
    // "Other"; an ungrouped ?ds= local export or server listing goes in flat when there are no
    // groups), arrow keys that only move the highlight. Opening happens on Enter or a click.
    // A ?data= folder brings its own groups (S.dataGroups, from its datasets.json).
    if (!window.SodaDatasetPicker || !$("dsPicker")) { status("This page is out of date: reload it (Ctrl+F5).", true); return; }
    dsPicker = window.SodaDatasetPicker.create({ host: $("dsPicker"), labelledBy: "dsLabel", onPick: (id) => openDataset(id) });
    dsPicker.setItems(S.datasets, S.dataGroups || cfg.groups || {});
    Object.assign(window.__spatialViewer, { dsPicker, openDataset });   // for testing from the console
    $("btnFit").onclick = fitView;
    // One button for both directions: if anything is on it turns everything off,
    // otherwise it turns everything on. Rebuilding the whole channel list would refetch
    // every histogram, so just flip the checkboxes that are already there.
    // Two buttons rather than one toggle: "turn everything on" was otherwise two clicks
    // (off, then on) and you could not tell which way the toggle would go.
    const setAllChannels = (on) => {
      S.channels.forEach((c) => (c.visible = on));
      const rows = $("chanList").querySelectorAll(".chan");
      S.channels.forEach((c, i) => {
        const row = rows[i]; if (!row) return;
        row.querySelector(".cb").checked = on;
        row.classList.toggle("off", !on);
      });
      render();
    };
    $("btnAllOn").onclick = () => setAllChannels(true);
    $("btnAllOff").onclick = () => setAllChannels(false);
    $("imgOpacity").oninput = (e) => { S.imgOpacity = +e.target.value; requestRender(); };
    $("btnImage").onclick = () => { S.imgOn = !S.imgOn; syncImageButton(); render(); };
    $("geneOpacity").oninput = (e) => { S.geneOpacity = +e.target.value; requestRender(); };
    $("exprSelect").onchange = (e) => setExprKey(e.target.value);
    $("geneSearch").oninput = (e) => showGeneResults(e.target.value);
    $("geneSearch").onkeydown = (e) => {
      if (e.key === "Enter" || e.key === "Return" || e.keyCode === 13) {
        e.preventDefault();
        const hot = $("geneResults").querySelector(".hot span") || $("geneResults").querySelector("div span");
        if (hot) addGene(hot.textContent);
        else if (e.target.value.trim()) { status(`No ${wd().gene} named "${e.target.value.trim()}" in this layer`, true); setTimeout(() => status(null), 3500); }
      }
      if (e.key === "Escape") $("geneResults").hidden = true;
    };
    // a marker panel: clicking into the box lists every marker, nothing to remember or type
    $("geneSearch").onfocus = (e) => { if (isProt()) showGeneResults(e.target.value); };
    $("geneSearch").onblur = () => setTimeout(() => ($("geneResults").hidden = true), 150);
    const syncGeneMode = () => {
      const heat = S.geneBlend === "heatmap";
      $("rowGenePalette").hidden = !heat;
    };
    $("geneBlend").onchange = (e) => {
      S.geneBlend = e.target.value; S.heat = null;
      S.genes.forEach(recolorGene); syncGeneMode(); buildGeneUI(); render();
    };
    $("genePalette").onchange = (e) => { S.genePalette = e.target.value; S.heat = null; render(); };
    syncGeneMode();
    $("pointScale").oninput = (e) => { S.pointScale = +e.target.value; requestRender(); };
    $("molSearch").oninput = (e) => showMolResults(e.target.value);
    $("molSearch").onfocus = (e) => showMolResults(e.target.value);
    $("molSearch").onblur = () => setTimeout(() => ($("molResults").hidden = true), 150);
    $("molSearch").onkeydown = (e) => {
      if (e.key === "Enter" || e.keyCode === 13) {
        e.preventDefault();
        const hot = $("molResults").querySelector(".hot span") || $("molResults").querySelector("div span");
        if (hot) addMolGene(hot.textContent);
        else if (e.target.value.trim()) { status(`No gene named "${e.target.value.trim()}" among the molecules`, true); setTimeout(() => status(null), 3500); }
      }
      if (e.key === "Escape") $("molResults").hidden = true;
    };
    $("molSize").oninput = (e) => { S.mol.size = +e.target.value; requestRender(); };
    $("molOpacity").oninput = (e) => { S.mol.opacity = +e.target.value; requestRender(); };
    $("tissueBins").onchange = async (e) => { S.tissue.on = e.target.checked; if (S.tissue.on && !S.tissue.positions) await loadTissue(); render(); };
    $("qcOutline").onchange = (e) => qcSetOutline(e.target.checked);
    // The description is held with the data, not in this page. Anything written
    // into index.html is readable by anyone who opens the page - including the
    // parts of the project that are not public yet - so the body is fetched from
    // the private container, behind the same password and the same expiring
    // credential the imaging data uses.
    let aboutLoaded = false;
    // about.html from a folder named in the link: rebuild it from an allow-list of tags,
    // so no script, handler, style or frame from someone else's file reaches this page
    function safeAbout(html) {
      const KEEP = new Set(["P", "H1", "H2", "H3", "H4", "UL", "OL", "LI", "B", "STRONG", "I", "EM", "BR", "CODE", "A", "SMALL"]);
      const src = new DOMParser().parseFromString(html, "text/html").body;
      const copy = (from, to) => {
        for (const n of from.childNodes) {
          if (n.nodeType === 3) { to.appendChild(document.createTextNode(n.nodeValue)); continue; }
          if (n.nodeType !== 1) continue;
          if (!KEEP.has(n.tagName)) { if (!/^(SCRIPT|STYLE|IFRAME|OBJECT|TEMPLATE)$/.test(n.tagName)) copy(n, to); continue; }
          const el = document.createElement(n.tagName);
          if (n.tagName === "A") {
            const href = n.getAttribute("href") || "";
            if (/^https?:\/\//i.test(href)) { el.href = href; el.target = "_blank"; el.rel = "noopener noreferrer"; }
          }
          copy(n, el); to.appendChild(el);
        }
      };
      const frag = document.createDocumentFragment();
      copy(src, frag);
      return frag;
    }
    async function loadAbout() {
      if (aboutLoaded) return;
      const body = $("aboutBody");
      body.innerHTML = '<p class="pending">Loading…</p>';
      try {
        const staticRoot = params.get("static");
        let url, trusted = false;
        if (DATA) url = `${DATA.root}about.html${DATA.qs}`;
        else if (staticRoot) url = `${staticRoot.replace(/\/$/, "")}/about.html`;
        else if (TOKEN_URL) { const t = await ensureToken(); url = `${t.base}about.html${tokenQuery()}`; trusted = true; }
        else throw new Error("no source configured");
        const r = await dataFetch(url, { cache: "no-cache" });
        if (r.status === 404) {
          body.innerHTML = '<p class="pending">No description has been published yet.</p>';
          return;                                   // not an error: it just is not written
        }
        if (!r.ok) throw new Error(`${r.status}`);
        // our own container's copy goes in as written; a folder named in the link could be
        // anyone's, so its text keeps only plain formatting and https links
        const html = await r.text();
        if (trusted) body.innerHTML = html; else body.replaceChildren(safeAbout(html));
        aboutLoaded = true;
      } catch (e) {
        body.innerHTML = `<p class="pending">Could not load the description: ${e.message}</p>`;
      }
    }
    // Everything naming this particular project lives in config, so the viewer
    // itself carries no dataset, no institution and no byline. A different study
    // points it at different data and supplies its own three lines.
    const brand = cfg.brand || {};
    if (brand.title) { $("brandTitle").textContent = brand.title; document.title = brand.title; }
    $("brandLogos").innerHTML = (brand.logos || []).map((l) =>
      `<img class="${l.className || ""}" src="${l.src}" alt="${l.alt || ""}" />`).join("");
    $("colophon").textContent = brand.byline || "";
    if (dataTitle) { $("brandTitle").textContent = dataTitle; document.title = `${dataTitle} · SODA`; }
    if (dataByline) $("colophon").textContent = dataByline;

    // someone else's folder (?data= on the link) never shows this site's own paper
    const paper = $("lnkPaper"), paperUrl = DATA && DATA.fromLink ? dataPaper : dataPaper || (cfg.paperUrl || "").trim();
    if (paperUrl) {
      paper.href = paperUrl;
      paper.title = paperUrl;
    } else {
      paper.classList.add("pending");
      paper.removeAttribute("href");                 // not a link until there is one
      paper.title = "Added once the paper is published";
      paper.onclick = (e) => e.preventDefault();
    }
    const about = (open) => { $("aboutBack").hidden = !open; if (open) loadAbout(); };
    $("btnAbout").onclick = () => about(true);
    $("aboutClose").onclick = () => about(false);
    // clicking the backdrop closes; clicking inside the dialog must not
    $("aboutBack").onclick = (e) => { if (e.target === $("aboutBack")) about(false); };
    document.addEventListener("keydown", (e) => { if (e.key === "Escape") about(false); });
    $("chanSearch").oninput = () => filterChannelUI();
    const syncRnaMode = () => { $("rowRnaBin").hidden = !(S.rna && S.rna.view === "density"); };
    $("rnaView").onchange = (e) => {
      if (!S.rna) return;
      S.rna.view = e.target.value; S.rna.heat = null; S.rna.grid = "";
      syncRnaMode(); renderRna();
    };
    $("rnaPalette").onchange = (e) => {
      if (!S.rna) return;
      S.rna.palette = e.target.value; S.rna.heat = null;
      S.rna.genes.forEach(recolorRnaGene); buildRnaGeneUI(); renderRna();
    };
    $("rnaBin").onchange = (e) => { if (!S.rna) return; S.rna.bin = +e.target.value; S.rna.heat = null; renderRna(); };
    $("rnaLabels").onchange = (e) => { if (S.rna) { S.rna.showLabels = e.target.checked; renderRna(); } };
    $("rnaSearch").oninput = (e) => showRnaResults(e.target.value);
    $("rnaSearch").onkeydown = (e) => {
      if (e.key === "Enter") {
        e.preventDefault();
        const hot = $("rnaResults").querySelector(".hot span") || $("rnaResults").querySelector("div span");
        if (hot) addRnaGene(hot.textContent);
        else if (e.target.value.trim()) { status(`No gene named "${e.target.value.trim()}" in this dataset`, true); setTimeout(() => status(null), 3500); }
      }
      if (e.key === "Escape") $("rnaResults").hidden = true;
    };
    $("rnaSearch").onblur = () => setTimeout(() => ($("rnaResults").hidden = true), 150);
    $("rnaPointSize").oninput = (e) => { if (S.rna) { S.rna.pointSize = +e.target.value; renderRna(); } };
    $("rnaClear").onclick = () => {
      if (!S.rna) return;
      S.rna.genes = []; S.rna.heat = null; S.rna.grid = "";
      buildRnaGeneUI(); renderRna();
    };
    window.addEventListener("resize", () => updateHud());
    let wasStacked = null;
    if (window.ResizeObserver) new ResizeObserver(() => {
      const el = $("viewport");
      if (!el.clientWidth || !el.clientHeight) return;
      // side by side flips between left | right and top / bottom as the screen turns: re-lay
      // the panes and frame both afresh, as the old framing no longer fits either
      const st = sbsStacked(), flipped = wasStacked !== null && st !== wasStacked && umapSbs();
      wasStacked = st;
      if (flipped) umapApplyViews();
      if (fitPending || flipped || (autoFit && S.mode !== "rna")) fitView();   // the RNA panels refit themselves
      if (renderPending) render(); else updateHud();
    }).observe($("viewport"));

    if (!S.datasets.length) { status("No datasets: add entries to config.js (static links) or set apiBase", true); return; }
    const want = params.get("ds") || dataDefault || cfg.defaultDataset;
    const first = S.datasets.find((d) => d.id === want) || S.datasets[0];
    await openDataset(first.id);
  }

  // ---- phones: the controls sheet under the map --------------------------------------------
  // Tap the handle to open or close it; drag it to any height and it settles at the nearest of
  // a strip, half the screen or nearly all of it. The handle only shows on narrow screens
  // (styles.css), so on a desktop none of this ever fires.
  function sheetEvents() {
    const app = document.querySelector(".app"), side = $("sidebar"), h = $("sheetHandle");
    const set = (s) => {
      app.dataset.sheet = s; side.style.height = "";
      h.setAttribute("aria-expanded", String(s !== "peek"));
      if (s === "peek") side.scrollTop = 0;
    };
    let drag = null, dragged = false;
    h.addEventListener("pointerdown", (ev) => {
      if (ev.button !== 0) return;
      h.setPointerCapture(ev.pointerId);
      drag = { y: ev.clientY, h0: side.offsetHeight, moved: false };
    });
    h.addEventListener("pointermove", (ev) => {
      if (!drag) return;
      const dy = drag.y - ev.clientY;
      if (Math.abs(dy) > 6) drag.moved = true;
      if (drag.moved) side.style.height = `${Math.max(h.offsetHeight, Math.min(window.innerHeight - 60, drag.h0 + dy))}px`;
    });
    const end = () => {
      const D = drag; drag = null;
      if (!D || !D.moved) return;
      dragged = true;                                          // the click that follows is not a tap
      const H = window.innerHeight, now = side.offsetHeight;
      const stops = { peek: h.offsetHeight, half: H * 0.45, full: H * 0.88 };
      set(Object.keys(stops).reduce((a, b) => (Math.abs(stops[b] - now) < Math.abs(stops[a] - now) ? b : a)));
    };
    h.addEventListener("pointerup", end);
    h.addEventListener("pointercancel", end);
    h.addEventListener("click", () => {
      if (dragged) { dragged = false; return; }
      set(app.dataset.sheet === "peek" ? "half" : "peek");
    });
  }

  // ---- touch screens: what a keyboard does on a desktop ---------------------------------------
  // "Done drawing" stands in for Esc while a drawing pad is up (the button is touch-only in
  // styles.css). A pad going up with the controls sheet at full height also drops the sheet
  // to its strip, as the map is then too small to draw on.
  function touchInit() {
    const done = $("drawDone"), pads = [$("drawPad"), $("selPad")], handle = $("sheetHandle");
    const sync = () => {
      const was = done.hidden;
      done.hidden = pads.every((p) => p.hidden);
      if (was && !done.hidden && document.querySelector(".app").dataset.sheet === "full" &&
          getComputedStyle(handle).display !== "none") handle.click();
    };
    const mo = new MutationObserver(sync);
    for (const p of pads) mo.observe(p, { attributes: true, attributeFilter: ["hidden"] });
    done.onclick = () => { lassoSetTool(null); umapSetTool(null); };
    // A phone freezes a page it hides (another app, the screen locked anyway): the screen lock is
    // dropped and a running UMAP waits. Coming back, take the lock again and say why it took longer.
    let hiddenRun = false;
    document.addEventListener("visibilitychange", () => {
      const L = S.lasso; if (!TOUCH || !L) return;
      if (document.hidden) { hiddenRun = L.running; return; }
      if (!hiddenRun || !L.running) { hiddenRun = false; return; }
      hiddenRun = false;
      lassoWake(true);
      const m = "The UMAP was paused while this page was hidden; it carries on now (its time includes the pause).";
      status(m);
      setTimeout(() => { if ($("status").textContent === m) status(null); }, 7000);
    });
  }

  init();
})();
