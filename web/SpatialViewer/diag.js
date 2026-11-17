// SODA diagnostic panel. Add &diag=1 to the page address (after signing in as usual) and a
// panel lists what this browser can do and where the image path stops: the browser and
// its graphics, one manifest and one tile fetched the way the viewer fetches them, each
// way of decoding that tile, whether it reaches a texture and draws, what the map holds,
// and the first errors caught. Made for phones, where there is no console to open.
//
// It never prints a password, a token or a signed address: clean() drops everything after
// "?" from anything shown or copied, and the fetches go through the viewer's own functions,
// so the panel never sees the password at all.
//
// Loaded before viewer.js so errors during start-up are caught too. Without ?diag=1 it only
// keeps a short list of errors (window.__sodaDiag.errors) and draws nothing.
(function () {
  "use strict";
  // ---- pure helpers (also used by tests/diag_test.js under node) -------------------
  function clean(t) {
    return String(t == null ? "" : t)
      .replace(/\?[^\s"'<>)\]]*/g, "")                                        // everything after "?" in an address
      .replace(/\b(sig|se|sp|sv|sr|st|spr|skoid|sktid|skt|ske|sks|skv|token|sas|password|pw)=[^\s&"'<>)\]]*/gi, "$1=[hidden]")
      .replace(/\b(Basic|Bearer)\s+[A-Za-z0-9+/=._-]{6,}/g, "$1 [hidden]");
  }
  function parseUA(ua, touchPoints) {
    ua = String(ua || "");
    const os = ua.match(/(iPhone|iPad|iPod)[^)]*? OS (\d+)[_.](\d+)(?:[_.](\d+))?/);
    const ver = ua.match(/Version\/([\d.]+)/);
    const other = ua.match(/(CriOS|FxiOS|EdgiOS|OPiOS|DuckDuckGo)\/([\d.]+)/);
    const out = { device: "", os: "", browser: "", ios: false };
    if (os) {
      out.ios = true; out.device = os[1];
      out.os = `iOS ${os[2]}.${os[3]}${os[4] ? "." + os[4] : ""}`;
      out.browser = other ? `${other[1]} ${other[2]} (WebKit underneath, as every iOS browser)` : ver ? `Safari ${ver[1]}` : "an in-app browser (WebKit)";
    } else if (/Macintosh/.test(ua) && touchPoints > 1) {
      out.ios = true; out.device = "iPad (desktop-style address)"; out.os = "iPadOS, version not given";
      out.browser = ver ? `Safari ${ver[1]}` : "WebKit";
    } else {
      const m = ua.match(/(Edg|Chrome|Firefox)\/([\d.]+)/) || (ver && /Safari/.test(ua) ? ["", "Safari", ver[1]] : null);
      out.device = /Android/.test(ua) ? "Android" : /Windows/.test(ua) ? "Windows" : /Macintosh/.test(ua) ? "Mac" : /Linux/.test(ua) ? "Linux" : "unknown";
      out.browser = m ? `${m[1] === "Edg" ? "Edge" : m[1]} ${m[2]}` : "unknown";
    }
    return out;
  }
  // what the first bytes of a tile say it is
  function imageInfo(bytes) {
    const b = bytes;
    if (b.length > 24 && b[0] === 0x89 && b[1] === 0x50 && b[2] === 0x4e && b[3] === 0x47) {
      const w = (b[16] << 24 | b[17] << 16 | b[18] << 8 | b[19]) >>> 0, h = (b[20] << 24 | b[21] << 16 | b[22] << 8 | b[23]) >>> 0;
      const kinds = { 0: "grey", 2: "RGB", 3: "palette", 4: "grey + alpha", 6: "RGB + alpha" };
      return `PNG ${w} x ${h}, ${b[24]}-bit ${kinds[b[25]] || "type " + b[25]}`;
    }
    if (b.length > 3 && b[0] === 0xff && b[1] === 0xd8) {
      for (let i = 2; i + 9 < b.length;) {
        if (b[i] !== 0xff) break;
        const m = b[i + 1], len = b[i + 2] << 8 | b[i + 3];
        if (m >= 0xc0 && m <= 0xc2) {
          const n = b[i + 9];
          return `JPEG ${b[i + 7] << 8 | b[i + 8]} x ${b[i + 5] << 8 | b[i + 6]}, ${b[i + 4]}-bit ${n === 1 ? "grey" : n + " channels"}${m === 0xc2 ? ", progressive" : ""}`;
        }
        i += 2 + len;
      }
      return "JPEG";
    }
    if (b.length > 1 && (b[0] === 0x3c || b[0] === 0x7b)) return "NOT an image: text (an error page?)";
    return "NOT a PNG or JPEG";
  }
  if (typeof module === "object" && module.exports) { module.exports = { clean, parseUA, imageInfo }; return; }

  // ---- always on: keep the first errors ------------------------------------------------
  const T0 = Date.now(), MAX_ERR = 40;
  const D = (window.__sodaDiag = { errors: [], dropped: 0, clean });
  D.note = (kind, msg) => {
    if (D.errors.length >= MAX_ERR) { D.dropped++; return; }
    D.errors.push({ t: ((Date.now() - T0) / 1000).toFixed(1), kind, msg: clean(msg).slice(0, 300) });
    if (D.onNote) D.onNote();
  };
  window.addEventListener("error", (e) => {
    const el = e.target;
    // vendor/deck.gl not there: libs.js takes the CDN's copy instead, so that is not a problem in itself
    if (el && el.dataset && el.dataset.sodaLib) return;
    if (el && el !== window && (el.src || el.href)) D.note("load", `could not load ${el.src || el.href}`);
    else D.note("script", `${e.message || "error"}${e.filename ? ` (${String(e.filename).split("/").pop()}:${e.lineno})` : ""}`);
  }, true);
  window.addEventListener("unhandledrejection", (e) => {
    const r = e.reason;
    D.note("promise", r && (r.message ? `${r.name || "Error"}: ${r.message}` : String(r)));
  });
  if (new URLSearchParams(location.search).get("diag") !== "1") return;

  // ---- the panel ---------------------------------------------------------------------
  const css = `
#diagPanel { position: fixed; left: 0; right: 0; bottom: 0; z-index: 10000; max-height: 62vh; display: flex; flex-direction: column;
  background: #10141c; color: #e8ecf2; border-top: 2px solid #4da3ff; font: 16px/1.45 -apple-system, "Segoe UI", Roboto, sans-serif;
  padding-bottom: env(safe-area-inset-bottom); -webkit-text-size-adjust: 100%; }
#diagPanel.small #diagText { display: none; }
#diagBar { display: flex; gap: 8px; align-items: center; padding: 8px 10px; flex-wrap: wrap; }
#diagBar b { flex: 1 1 auto; font-size: 17px; }
#diagBar button { font: inherit; font-size: 16px; padding: 8px 14px; border-radius: 8px; border: 1px solid #4da3ff; background: #17324f; color: #fff; }
#diagText { overflow: auto; -webkit-overflow-scrolling: touch; padding: 2px 12px 14px; margin: 0; white-space: pre-wrap; word-break: break-word;
  -webkit-user-select: text; user-select: text; font: 15px/1.5 ui-monospace, Menlo, Consolas, monospace; }
#diagText .bad { color: #ff8f8f; font-weight: 700; }
#diagText .good { color: #8fe3a0; }
#diagText .head { color: #7fc0ff; font-weight: 700; }`;
  let panel, textEl, titleEl, lines = [], plain = "";
  function build() {
    if (panel) return;
    const st = document.createElement("style"); st.textContent = css; document.head.appendChild(st);
    panel = document.createElement("div"); panel.id = "diagPanel"; panel.className = "small";
    panel.innerHTML = `<div id="diagBar"><b id="diagTitle">Diagnostic: starting…</b>` +
      `<button id="diagCopy" type="button">Copy</button><button id="diagShare" type="button" hidden>Share</button>` +
      `<button id="diagRun" type="button">Run again</button><button id="diagHide" type="button">Show</button></div><pre id="diagText"></pre>`;
    document.body.appendChild(panel);
    textEl = panel.querySelector("#diagText"); titleEl = panel.querySelector("#diagTitle");
    panel.querySelector("#diagHide").onclick = () => setSmall(!panel.classList.contains("small"));
    panel.querySelector("#diagRun").onclick = () => run();
    panel.querySelector("#diagCopy").onclick = copy;
    if (navigator.share) { const s = panel.querySelector("#diagShare"); s.hidden = false; s.onclick = () => navigator.share({ title: "SODA diagnostic", text: plain }).catch(() => {}); }
  }
  function setSmall(on) { panel.classList.toggle("small", on); panel.querySelector("#diagHide").textContent = on ? "Show" : "Hide"; }
  async function copy() {
    const btn = panel.querySelector("#diagCopy");
    let ok = false;
    try { await navigator.clipboard.writeText(plain); ok = true; } catch (e) {
      try {                                             // older Safari, or a page that is not https
        const ta = document.createElement("textarea"); ta.value = plain; ta.setAttribute("readonly", "");
        ta.style.cssText = "position:fixed;left:0;top:0;opacity:0;font-size:16px"; document.body.appendChild(ta);
        ta.select(); ta.setSelectionRange(0, plain.length); ok = document.execCommand("copy"); ta.remove();
      } catch (e2) { ok = false; }
    }
    btn.textContent = ok ? "Copied" : "Select the text and copy it";
    setTimeout(() => { btn.textContent = "Copy"; }, 2500);
  }
  // one report line: label, value, and whether it is fine ("good"), a problem ("bad") or plain
  function add(label, value, mark) { lines.push({ label, value: clean(value), mark: mark || "" }); paint(); }
  function head(t) { lines.push({ head: t }); paint(); }
  function paint() {
    if (!textEl) return;
    textEl.textContent = "";
    const out = [];
    for (const l of lines) {
      const row = document.createElement("div");
      if (l.head) { row.className = "head"; row.textContent = `\n${l.head}`; out.push("", `== ${l.head} ==`); }
      else {
        const tag = l.mark === "bad" ? "PROBLEM  " : "";
        row.textContent = `${tag}${l.label}: ${l.value}`; if (l.mark) row.className = l.mark;
        out.push(`${tag}${l.label}: ${l.value}`);
      }
      textEl.appendChild(row);
    }
    plain = out.join("\n").trim();
  }
  const sleep = (ms) => new Promise((r) => setTimeout(r, ms));
  const errText = (e) => clean(e && e.name ? `${e.name}: ${e.message}` : String(e));
  const timed = async (fn) => { const t = performance.now(); const v = await fn(); return [v, Math.round(performance.now() - t)]; };

  // ---- graphics: a scratch WebGL2 context, the same calls the map makes -------------------
  function glInfo() {
    const c = document.createElement("canvas"); c.width = c.height = 64;
    let gl = null, why = "";
    try { gl = c.getContext("webgl2", { preserveDrawingBuffer: true }); } catch (e) { why = errText(e); }
    if (!gl) {
      let gl1 = null; try { gl1 = c.getContext("webgl"); } catch (e) { /* none */ }
      add("WebGL2", `NOT available${why ? " (" + why + ")" : ""}; WebGL1 ${gl1 ? "is" : "is not either"}. The map needs WebGL2: iOS 15 or newer, and Lockdown Mode off for this site`, "bad");
      return null;
    }
    const dbg = gl.getExtension("WEBGL_debug_renderer_info");
    add("WebGL2", `${gl.getParameter(gl.VERSION)} · ${gl.getParameter(gl.SHADING_LANGUAGE_VERSION)}`, "good");
    add("Graphics chip", `${dbg ? gl.getParameter(dbg.UNMASKED_RENDERER_WEBGL) : gl.getParameter(gl.RENDERER)} (${dbg ? gl.getParameter(dbg.UNMASKED_VENDOR_WEBGL) : gl.getParameter(gl.VENDOR)})`);
    const mts = gl.getParameter(gl.MAX_TEXTURE_SIZE);
    add("Largest texture", `${mts} px (tiles are 1024 or smaller); largest drawing surface ${gl.getParameter(gl.MAX_RENDERBUFFER_SIZE)} px`, mts >= 1024 ? "good" : "bad");
    return gl;
  }
  // upload a decoded tile and read it back, then draw it through a shader with MAX blending
  // (what ChannelBitmapLayer does) and read the canvas
  function textureTest(gl, src, how) {
    const w = src.naturalWidth || src.width, h = src.naturalHeight || src.height;
    const tex = gl.createTexture();
    gl.bindTexture(gl.TEXTURE_2D, tex);
    gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);
    while (gl.getError()) { /* clear */ }
    if (how === "storage") { gl.texStorage2D(gl.TEXTURE_2D, 1, gl.RGBA8, w, h); gl.texSubImage2D(gl.TEXTURE_2D, 0, 0, 0, gl.RGBA, gl.UNSIGNED_BYTE, src); }
    else gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, gl.RGBA, gl.UNSIGNED_BYTE, src);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR); gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE); gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    const upErr = gl.getError();
    // read the texture itself
    const fb = gl.createFramebuffer(); gl.bindFramebuffer(gl.FRAMEBUFFER, fb);
    gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, tex, 0);
    const px = new Uint8Array(w * h * 4);
    let texMax = -1, texMean = 0;
    if (gl.checkFramebufferStatus(gl.FRAMEBUFFER) === gl.FRAMEBUFFER_COMPLETE) {
      gl.readPixels(0, 0, w, h, gl.RGBA, gl.UNSIGNED_BYTE, px);
      let s = 0; texMax = 0; for (let i = 0; i < px.length; i += 4) { if (px[i] > texMax) texMax = px[i]; s += px[i]; }
      texMean = s / (w * h);
    }
    gl.bindFramebuffer(gl.FRAMEBUFFER, null); gl.deleteFramebuffer(fb);
    // draw it: grey -> colour in the fragment shader, MAX blending
    const sh = (type, code) => { const s = gl.createShader(type); gl.shaderSource(s, code); gl.compileShader(s); if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) throw new Error("shader: " + gl.getShaderInfoLog(s)); return s; };
    const prog = gl.createProgram();
    gl.attachShader(prog, sh(gl.VERTEX_SHADER, "#version 300 es\nout vec2 uv; void main(){ vec2 p = vec2(float(gl_VertexID & 1), float(gl_VertexID >> 1)); uv = p; gl_Position = vec4(p * 2.0 - 1.0, 0.0, 1.0); }"));
    gl.attachShader(prog, sh(gl.FRAGMENT_SHADER, "#version 300 es\nprecision highp float; uniform sampler2D t; in vec2 uv; out vec4 o; void main(){ float g = texture(t, uv).r; o = vec4(vec3(g), g); }"));
    gl.linkProgram(prog);
    if (!gl.getProgramParameter(prog, gl.LINK_STATUS)) throw new Error("program: " + gl.getProgramInfoLog(prog));
    gl.canvas.width = gl.canvas.height = 64; gl.viewport(0, 0, 64, 64);
    gl.clearColor(0, 0, 0, 0); gl.clear(gl.COLOR_BUFFER_BIT);
    gl.useProgram(prog); gl.activeTexture(gl.TEXTURE0); gl.bindTexture(gl.TEXTURE_2D, tex); gl.uniform1i(gl.getUniformLocation(prog, "t"), 0);
    gl.enable(gl.BLEND); gl.blendEquation(gl.MAX); gl.blendFunc(gl.ONE, gl.ONE);
    gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);
    const out = new Uint8Array(64 * 64 * 4); gl.readPixels(0, 0, 64, 64, gl.RGBA, gl.UNSIGNED_BYTE, out);
    let drawMax = 0; for (let i = 0; i < out.length; i += 4) if (out[i] > drawMax) drawMax = out[i];
    const drawErr = gl.getError();
    gl.disable(gl.BLEND); gl.deleteProgram(prog); gl.deleteTexture(tex);
    return { upErr, texMax, texMean, drawMax, drawErr };
  }

  // ---- the report ------------------------------------------------------------------------
  let running = false;
  async function run() {
    if (running) return;
    running = true; lines = []; build();
    titleEl.textContent = "Diagnostic: running…";
    const verdicts = [];
    try {
      const now = new Date();
      head("SODA diagnostic");
      add("When", `${now.toISOString()} (page open ${((Date.now() - T0) / 1000).toFixed(0)} s)`);
      add("Page", location.origin + location.pathname);
      const vOf = (name) => { const s = document.querySelector(`script[src*="${name}"],link[href*="${name}"]`); const m = s && (s.src || s.href || "").match(/[?&]v=([\w.-]+)/); return m ? m[1] : s ? "no version" : "not on the page"; };
      lines.push({ label: "Files", value: `viewer.js ${vOf("viewer.js")} · diag.js ${vOf("diag.js")} · styles.css ${vOf("styles.css")}`, mark: "" }); paint();

      head("Browser");
      const ua = parseUA(navigator.userAgent, navigator.maxTouchPoints || 0);
      add("Device", `${ua.device}${ua.os ? " · " + ua.os : ""} · ${ua.browser}`);
      if (ua.ios) add("Note", "from iOS 26 on, the iOS number in the user agent is frozen; the Safari number is the true one");
      add("User agent", navigator.userAgent);
      const vp = document.getElementById("viewport");
      add("Screen", `${innerWidth} x ${innerHeight} CSS px, ${window.devicePixelRatio || 1}x; touch ${matchMedia("(pointer: coarse)").matches ? "yes" : "no"}`);
      const vpOk = vp && vp.clientWidth > 0 && vp.clientHeight > 0;
      add("Map area", vp ? `${vp.clientWidth} x ${vp.clientHeight} CSS px` : "missing", vpOk ? "good" : "bad");
      if (!vpOk) verdicts.push("the map area has no size on this screen, so nothing can be drawn in it");
      add("Online / cookies / https", `${navigator.onLine ? "yes" : "NO"} / ${navigator.cookieEnabled ? "yes" : "NO"} / ${window.isSecureContext ? "yes" : "no"}`);

      head("Graphics");
      const from = window.SODA_LIBS_FROM && window.SODA_LIBS_FROM.deck === "local" ? "this site's vendor/" : "cdn.jsdelivr.net";
      add("deck.gl", window.deck ? `loaded, version ${window.deck.VERSION || "?"}, from ${from}` : "NOT loaded (the drawing library did not arrive: offline, or blocked)", window.deck ? "good" : "bad");
      if (!window.deck) verdicts.push("deck.gl (the drawing library, from this site's vendor/ or cdn.jsdelivr.net) did not load: a content blocker or the network stopped it, and there is no copy kept for offline use");
      const off = window.SodaOffline ? window.SodaOffline.state() : null;
      if (off) add("Offline use", `${off.state}${off.why ? " (" + off.why + ")" : ""}; service worker ${navigator.serviceWorker && navigator.serviceWorker.controller ? "in charge of this page" : "not in charge"}; ${off.kept.length} dataset(s) kept`);
      const gl = glInfo();
      if (!gl) verdicts.push("this browser has no WebGL2, which the map needs");
      const hasCIB = typeof createImageBitmap === "function";
      add("createImageBitmap", hasCIB ? "exists" : "MISSING (the viewer decodes through an img element instead)", hasCIB ? "good" : "");

      // the viewer and its dataset: wait for the sign-in and the manifest
      const sv = () => window.__spatialViewer;
      if (!sv() || !sv().diag) {
        add("Viewer", sv() ? "an older viewer.js without the diagnostic hooks: only the lines above apply" : "did not start (see errors below)", "bad");
        if (!sv()) verdicts.push("the viewer did not start");
      } else {
        titleEl.textContent = "Diagnostic: waiting for a dataset (sign in if asked)…";
        for (let i = 0; i < 240 && !sv().state.m; i++) {
          const st = document.getElementById("status");
          if (st && !st.hidden && st.classList.contains("err") && i > 4) break;
          await sleep(500);
        }
        titleEl.textContent = "Diagnostic: running…";
        const S = sv().state;
        head("Data");
        if (!S.m) {
          const st = document.getElementById("status");
          add("Dataset", `did not open. The page says: ${st && !st.hidden ? st.textContent : "(nothing)"}`, "bad");
          verdicts.push("the dataset did not open (its manifest could not be read): see the Data lines");
        }
        const info = sv().diag();
        add("Source", `${info.kind || "none"}${info.tokenFlow ? `, signed in through the token service; pass valid for ${info.tokenLeftMin == null ? "? (none yet)" : info.tokenLeftMin + " more min"}` : ""}`);
        if (info.manifestUrl) {
          try {
            const [r, ms] = await timed(async () => { const r = await info.fetch(info.manifestUrl); return { status: r.status, type: r.headers.get("content-type"), size: (await r.blob()).size }; });
            add("Manifest", `HTTP ${r.status}, ${r.size.toLocaleString()} bytes, ${r.type || "no type"}, ${ms} ms — ${info.manifestUrl}`, r.status === 200 ? "good" : "bad");
            if (r.status !== 200) verdicts.push(`the manifest answered HTTP ${r.status}`);
          } catch (e) { add("Manifest", `could not be fetched: ${errText(e)} — ${info.manifestUrl}`, "bad"); verdicts.push("the manifest could not be downloaded (network, CORS or a content blocker)"); }
        }
        if (S.m && !S.m.image) { add("Image", "this dataset has no image (single-cell data): open an image dataset and press Run again"); verdicts.push("no image in this dataset to check: open an image dataset and run again"); }
        if (S.m && S.m.image && S.channels && S.channels.length) {
          await sleep(3500);                               // let the map ask for its first tiles
          const img = S.m.image;
          add("Image", `${img.width} x ${img.height} px, ${img.levels} levels, ${img.tileSize} px ${S.m.tileFormat || "png"} tiles, ${img.channels.length} channels`);
          const ch = S.channels.find((c) => c.visible) || S.channels[0];
          const url = info.tileUrl(ch.index, 0, 0, 0);
          let blob = null;
          try {
            const ac = new AbortController();
            const [r, ms] = await timed(async () => { const r = await info.fetch(url, { signal: ac.signal }); return { status: r.status, type: r.headers.get("content-type"), blob: await r.blob() }; });
            const kind = imageInfo(new Uint8Array(await r.blob.slice(0, 4096).arrayBuffer()));
            const ok = r.status === 200 && !/^NOT/.test(kind);
            add("Tile", `HTTP ${r.status}, ${r.blob.size.toLocaleString()} bytes, ${r.type || "no type"}, ${kind}, ${ms} ms — ${url}`, ok ? "good" : "bad");
            if (ok) blob = r.blob; else verdicts.push(r.status !== 200 ? `image tiles answer HTTP ${r.status}` : "what arrives as a tile is not a picture");
          } catch (e) { add("Tile", `could not be fetched: ${errText(e)} — ${url}`, "bad"); verdicts.push("image tiles could not be downloaded although the manifest could"); }

          if (blob) {
            head("Decoding that tile");
            const ways = [
              ["createImageBitmap with options (the viewer's first choice)", hasCIB && (() => createImageBitmap(blob, { premultiplyAlpha: "none", colorSpaceConversion: "none" }))],
              ["createImageBitmap, no options", hasCIB && (() => createImageBitmap(blob))],
              ["img element", () => info.decoders[2][1](blob)],
            ];
            const decoded = [];
            for (const [name, fn] of ways) {
              if (!fn) { add(name, "not available"); continue; }
              try { const [im, ms] = await timed(fn); add(name, `ok, ${im.naturalWidth || im.width} x ${im.naturalHeight || im.height}, ${ms} ms`, "good"); decoded.push([name, im]); }
              catch (e) { add(name, `FAILED: ${errText(e)}`, "bad"); }
            }
            if (!decoded.length) verdicts.push("this browser cannot decode the tile in any of the three ways");
            else if (decoded.length < ways.filter((w) => w[1]).length) verdicts.push(`decoding fails with ${ways.filter((w) => w[1] && !decoded.some((d) => d[0] === w[0])).map((w) => w[0]).join(" and ")}; the viewer falls back to a way that works`);
            add("The viewer is using", `${info.decoder} (adding &decode=img to the address makes it use the img element; &decode=plain, no options)`);

            if (gl && decoded.length) {
              head("Texture and draw (scratch canvas)");
              // what the tile holds according to a plain 2D canvas, to tell a black tile from a failed upload
              let ref = -1;
              try {
                const im = decoded[decoded.length - 1][1], c2 = document.createElement("canvas");
                c2.width = Math.min(256, im.naturalWidth || im.width); c2.height = Math.min(256, im.naturalHeight || im.height);
                const g2 = c2.getContext("2d"); g2.drawImage(im, 0, 0, c2.width, c2.height);
                const d2 = g2.getImageData(0, 0, c2.width, c2.height).data; ref = 0;
                for (let i = 0; i < d2.length; i += 4) if (d2[i] > ref) ref = d2[i];
                c2.width = c2.height = 0;
                add("The tile on a 2D canvas (for comparison)", `brightest ${ref}${ref === 0 ? ": this tile is black, so the lines below cannot tell much. Turn on another stain and run again" : ""}`);
              } catch (e) { add("The tile on a 2D canvas (for comparison)", `FAILED: ${errText(e)}`); }
              const good = [], bad = [];
              for (const [name, im] of decoded) for (const how of ["image", "storage"]) {
                const label = `${name.split(" (")[0]} > ${how === "image" ? "texImage2D" : "texStorage2D + texSubImage2D"}`;
                try {
                  const t = textureTest(gl, im, how);
                  const ok = !t.upErr && !t.drawErr && t.texMax > 0 && t.drawMax > 0;
                  (ok ? good : bad).push(label);
                  add(label, `upload ${t.upErr ? "GL error " + t.upErr : "ok"}; read back brightest ${t.texMax}, mean ${t.texMean.toFixed(1)}; drawn with MAX blending brightest ${t.drawMax}${t.drawErr ? ", GL error " + t.drawErr : ""}`, ok ? "good" : ref === 0 ? "" : "bad");
                } catch (e) { bad.push(label); add(label, `FAILED: ${errText(e)}`, "bad"); }
              }
              if (ref !== 0 && !good.length) verdicts.push("a decoded tile does not reach a texture or does not draw (every way reads back black or fails)");
              else if (ref !== 0 && bad.length) verdicts.push(`the tile reads back black or fails this way: ${bad.join("; ")}` +
                (good.some((g) => /^img element/.test(g)) && bad.some((b) => /^createImageBitmap/.test(b)) ? " (adding &decode=img to the address makes the viewer use the img element, which works here)" : ""));
            }
          }

          head("The map itself");
          try {
            const deckgl = sv().deck, dgl = deckgl.device && deckgl.device.gl;
            if (!dgl) { add("Map context", "deck.gl has no graphics context (it could not start)", "bad"); verdicts.push("deck.gl could not start its graphics context"); }
            else {
              const lost = dgl.isContextLost();
              add("Map context", `${dgl.drawingBufferWidth} x ${dgl.drawingBufferHeight} px${lost ? ", LOST (the browser took the graphics memory back)" : ""}`, lost ? "bad" : "good");
              if (lost) verdicts.push("the map's graphics context is lost (memory): reload, with fewer stains on");
              const layers = (deckgl.layerManager ? deckgl.layerManager.getLayers() : []).filter((l) => /^ch-/.test(l.id) && l.state && l.state.tileset);
              let asked = 0, loaded = 0, pic = 0;
              for (const l of layers) for (const t of l.state.tileset.tiles) { asked++; if (t.isLoaded) loaded++; if (t.content) pic++; }
              const on = S.channels.filter((c) => c.visible).length;
              add("Stains on", `${on} of ${S.channels.length}; image ${S.imgOn ? "on" : "OFF (Image opacity at 0)"}, opacity ${S.imgOpacity}`, on && S.imgOn ? "" : "bad");
              add("Tiles", `${asked} asked, ${loaded} finished, ${pic} with a picture`, asked && pic === asked ? "good" : pic ? "" : "bad");
              if (!lost) {
                deckgl.redraw("diag");
                const w = dgl.drawingBufferWidth, h = dgl.drawingBufferHeight, px = new Uint8Array(w * h * 4);
                dgl.bindFramebuffer(dgl.FRAMEBUFFER, null);
                dgl.readPixels(0, 0, w, h, dgl.RGBA, dgl.UNSIGNED_BYTE, px);
                let lit = 0; for (let i = 0; i < px.length; i += 4) if (px[i] > 8 || px[i + 1] > 8 || px[i + 2] > 8) lit++;
                const pct = 100 * lit / (w * h);
                add("Map pixels that are not black", `${pct.toFixed(1)}%`, pct > 0.5 ? "good" : "bad");
                if (pct <= 0.5 && S.imgOn && on) verdicts.push(pic ? "tiles arrive and decode, but the map stays black: the drawing step fails on this device" : asked ? "the map asked for tiles but none produced a picture" : "the map has not asked for any tiles");
              }
            }
          } catch (e) { add("Map check", `FAILED: ${errText(e)}`, "bad"); }
          head("Memory");
          const mb = info.cached * img.tileSize * img.tileSize * 4 / 1e6;
          add("Decoded tiles kept", `${info.cached} of at most ${info.cacheMax} (about ${mb.toFixed(0)} MB if full size); touch limits ${info.touch ? "on" : "off"}`);
          add("Device hints", `memory ${navigator.deviceMemory ? navigator.deviceMemory + " GB" : "not told (Safari never says)"}; ${navigator.hardwareConcurrency || "?"} cores` +
            (performance.memory ? `; script memory ${(performance.memory.usedJSHeapSize / 1e6).toFixed(0)} of ${(performance.memory.jsHeapSizeLimit / 1e6).toFixed(0)} MB` : ""));
        }
      }
      head(`Errors caught (${D.errors.length}${D.dropped ? " + " + D.dropped + " more" : ""})`);
      if (!D.errors.length) add("None", "no errors since the page opened", "good");
      D.errors.slice(0, 12).forEach((e) => add(`${e.t} s ${e.kind}`, e.msg, "bad"));
      const st = document.getElementById("status");
      add("Status line on the page", st && !st.hidden ? st.textContent : "(empty)");
    } catch (e) { add("The diagnostic itself failed", errText(e), "bad"); }
    lines.unshift({ label: "RESULT", value: clean(verdicts.length ? verdicts.join("; ") : "these checks found no problem: the tile downloads, decodes, uploads and draws"), mark: verdicts.length ? "bad" : "good" });
    paint();
    titleEl.textContent = verdicts.length ? "Diagnostic: problem found" : "Diagnostic: no problem found";
    setSmall(false);
    running = false;
  }
  D.run = run;
  const start = () => { build(); setTimeout(run, 300); };
  if (document.readyState === "complete") start(); else window.addEventListener("load", start);
})();
