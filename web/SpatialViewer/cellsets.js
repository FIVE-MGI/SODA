// Cell sets: merge clusters of the clustering on show into named sets, undo, and export one row per bin.
//
// The sets live in viewer.js's state (S.umap.sets, per clustering key) next to the names given to clusters
// (S.umap.names); viewer.js applies them where the clusters are read, so every view, legend, group picker and
// count sees the sets (viewer.js: clSetsApply, clWithSets). A set takes the number of its lowest member: its
// colour is that cluster's, its name is that cluster's name entry (the member's own name is kept in `prev`
// and comes back on Un-merge). This file only builds the section, changes the sets and writes the CSV.
//
// The CSV holds positions.bin rows (bin_row) and their pixel positions, as the other exports of the page do:
// never a barcode or a sample ID.
/* exported SodaSets */
(function () {
  "use strict";
  let api = null, sig = "", picked = new Set(), undo = [], undoFor = null;
  const $ = (id) => document.getElementById(id);
  const MAX_NAME = 60, MAX_UNDO = 30;
  const esc = (x) => String(x).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
  const rgbCss = (c) => `rgb(${c[0]},${c[1]},${c[2]})`;
  const cut = (t) => (t.length > MAX_NAME ? t.slice(0, MAX_NAME - 1) + "…" : t);
  const W = () => (api && api.words ? api.words() : { bin: "bin", bins: "bins" });

  // ---- the section, in the Clusters tab under Clusters & selection -------------------------------------
  function buildDom() {
    const sec = document.createElement("section");
    sec.className = "panel"; sec.id = "panelSets"; sec.dataset.tab = "clusters"; sec.hidden = true;
    sec.innerHTML = `<h2>Cell sets</h2>
      <div id="csNote" class="meta cs-note"></div>
      <div class="btnrow" id="csShowRow" hidden><button id="csShowCl" type="button" title="Colour the bins by the clustering, so its clusters can be merged">Show the clusters</button></div>
      <div id="csPick" class="cs-pick" role="group" aria-label="Clusters to merge"></div>
      <label class="row slim"><span>Set name</span><input id="csName" type="text" maxlength="${MAX_NAME}" autocomplete="off" placeholder="e.g. T cells" title="The name of the new set (empty: the clusters' names joined with +)" /></label>
      <div class="btnrow">
        <button id="csMerge" type="button" disabled title="Merge the ticked clusters into one set">Merge</button>
        <button id="csUndo" type="button" disabled title="Undo the last merge, un-merge or rename made here">Undo</button>
        <button id="csCsv" type="button" title="One row per bin: its row number, position, cluster and set (no barcodes)">Export CSV</button>
      </div>
      <div id="csSets" class="cs-sets"></div>
      <div id="csMsg" class="meta cs-msg" role="status"></div>`;
    const side = $("sidebar");
    if (window.SodaLayout) { side.appendChild(sec); window.SodaLayout.adopt(sec); }
    else { const after = $("panelUmap"); if (after) after.after(sec); else side.appendChild(sec); }

    $("csMerge").onclick = merge;
    $("csUndo").onclick = undoLast;
    $("csCsv").onclick = exportCsv;
    $("csShowCl").onclick = () => api.showClusters();
    $("csName").addEventListener("keydown", (e) => { if (e.key === "Enter" && !$("csMerge").disabled) { e.preventDefault(); merge(); } });
    $("csPick").addEventListener("change", (e) => {
      const c = Number(e.target && e.target.dataset ? e.target.dataset.c : NaN);
      if (!Number.isInteger(c)) return;
      if (e.target.checked) picked.add(c); else picked.delete(c);
      e.target.closest(".cs-chip").classList.toggle("on", e.target.checked);
      syncButtons();
    });
  }
  const msg = (t, err) => { const el = $("csMsg"); if (!el) return; el.textContent = t || ""; el.style.color = err ? "#ff8f8f" : ""; };

  // ---- what is on show ----------------------------------------------------------------------------------
  // the clustering on show, with its codes loaded: { U, ck, c (its description, sets applied) }
  function now() {
    const S = api.state(), U = S && S.umap;
    if (!U || !api.available() || U.clRes == null || !U.cl) return null;
    const c = api.clustering();
    return c ? { U, ck: U.clRes, c } : null;
  }
  const setsOf = (cur) => cur.U.sets[cur.ck] || [];
  const setOf = (cur, head) => setsOf(cur).find((s) => s.ids[0] === head) || null;
  // bins per code of the codes on show (sets applied), counted once per array
  const sizeCache = new WeakMap();
  function sizes(cl) {
    let z = sizeCache.get(cl);
    if (!z) { z = new Map(); for (let r = 0; r < cl.length; r++) { const c = cl[r]; if (c !== api.NONE) z.set(c, (z.get(c) || 0) + 1); } sizeCache.set(cl, z); }
    return z;
  }
  // a cluster's own name: what it was called before it went into a set
  function ownName(cur, c) {
    const s = setOf(cur, c);
    return s ? (s.prev != null ? s.prev : api.baseName(c)) : api.name(c);
  }
  let uid = 0;
  const uids = new WeakMap();
  const idOf = (o) => { if (!uids.has(o)) uids.set(o, ++uid); return uids.get(o); };

  // ---- drawing the section: only when what it shows has changed -------------------------------------------
  function refresh() {
    const sec = $("panelSets"); if (!api || !sec) return;
    const ok = api.available();
    if (sec.hidden === ok) sec.hidden = !ok;
    if (!ok) return;
    const cur = now(), U = cur && cur.U;
    if (!cur || undoFor !== `${idOf(U)}|${cur.ck}`) { undo = []; undoFor = cur ? `${idOf(U)}|${cur.ck}` : null; }
    const key = cur ? [idOf(U), cur.ck, idOf(U.cl), U.setsVer, JSON.stringify(U.names[cur.ck] || {})].join("|") : "none";
    if (key === sig) { syncButtons(); return; }
    sig = key;
    const w = W();
    $("csShowRow").hidden = !!cur;
    $("csPick").hidden = !cur;
    if (!cur) {
      picked.clear();
      $("csNote").textContent = `Merge clusters into named sets: first colour the ${w.bins} by clusters (Clusters & selection, or the button below).`;
      $("csPick").innerHTML = ""; $("csSets").innerHTML = "";
      syncButtons();
      return;
    }
    const z = sizes(U.cl), codes = [...z.keys()].sort((a, b) => a - b), sets = setsOf(cur);
    for (const c of [...picked]) if (!z.has(c)) picked.delete(c);
    $("csNote").textContent = `${api.label(cur.c)}: ${codes.length} ${codes.length === 1 ? "group" : "groups"}` +
      (sets.length ? `, ${sets.length} of them ${sets.length === 1 ? "a set" : "sets"} of merged clusters` : "") +
      `. Tick two or more and press Merge.`;
    // one tick box per cluster or set, in the legend's colours
    $("csPick").innerHTML = codes.map((c) => {
      const on = picked.has(c), set = setOf(cur, c), nm = api.name(c);
      return `<label class="cs-chip${on ? " on" : ""}${set ? " set" : ""}" title="${esc(nm)}: ${z.get(c).toLocaleString()} ${esc(w.bins)}${set ? ` (a set of ${set.ids.length} clusters)` : ""}">` +
        `<input type="checkbox" data-c="${c}"${on ? " checked" : ""} /><i style="background:${rgbCss(api.rgb(c))}"></i><span>${esc(nm)}</span></label>`;
    }).join("");
    // the sets: a name to edit, the clusters in it, Un-merge
    const host = $("csSets");
    host.innerHTML = "";
    for (const s of sets) {
      const h = s.ids[0], row = document.createElement("div");
      row.className = "cs-set"; row.dataset.c = h;
      row.innerHTML = `<i style="background:${rgbCss(api.rgb(h))}"></i><input type="text" class="cs-setname" maxlength="${MAX_NAME}" aria-label="Name of this set" title="Rename this set (Enter)" />` +
        `<button type="button" class="cs-unmerge">Un-merge</button><small></small>`;
      const inp = row.querySelector("input"), name = api.name(h);
      inp.value = name;
      row.querySelector("small").textContent = `${s.ids.map((k) => ownName(cur, k)).join(" + ")} · ${(z.get(h) || 0).toLocaleString()} ${w.bins}`;
      row.querySelector("button").setAttribute("aria-label", `Un-merge ${name}`);
      row.querySelector("button").onclick = () => unmerge(h);
      inp.addEventListener("keydown", (e) => {
        if (e.key === "Enter") { e.preventDefault(); inp.blur(); }
        else if (e.key === "Escape") { e.stopPropagation(); inp.value = api.name(h); inp.blur(); }
      });
      inp.addEventListener("change", () => rename(h, inp.value));
      host.appendChild(row);
    }
    syncButtons();
  }
  function syncButtons() {
    const cur = api && now();
    const n = cur ? [...picked].length : 0;
    $("csMerge").disabled = n < 2;
    $("csMerge").textContent = n >= 2 ? `Merge ${n}` : "Merge";
    $("csUndo").disabled = !undo.length;
    $("csUndo").title = undo.length ? `Undo: ${undo[undo.length - 1].what}` : "Undo the last merge, un-merge or rename made here";
    $("csCsv").disabled = !cur;
  }

  // ---- changes: each one can be undone -----------------------------------------------------------------------
  // what an undo puts back: the clustering's sets, the names of the clusters it touched, the hidden ones
  function remember(cur, codes, what) {
    const own = cur.U.names[cur.ck] || {}, names = {};
    for (const k of codes) names[k] = Object.prototype.hasOwnProperty.call(own, k) ? own[k] : undefined;
    undo.push({ ck: cur.ck, what, sets: JSON.parse(JSON.stringify(setsOf(cur))), names, hidden: [...cur.U.hidden], pick: cur.U.pick });
    if (undo.length > MAX_UNDO) undo.shift();
  }
  function namesOf(cur) { return cur.U.names[cur.ck] || (cur.U.names[cur.ck] = {}); }
  function setName(names, k, v) { if (v == null) delete names[k]; else names[k] = v; }

  function merge() {
    const cur = now(); if (!cur) return;
    const z = sizes(cur.U.cl), pick = [...picked].filter((c) => z.has(c)).sort((a, b) => a - b);
    if (pick.length < 2) { msg("Tick two or more clusters to merge.", true); return; }
    const sets = setsOf(cur), olds = pick.map((c) => setOf(cur, c)).filter(Boolean);
    const typed = $("csName").value.replace(/\s+/g, " ").trim();
    const name = cut(typed || (olds.length === 1 ? api.name(olds[0].ids[0]) : pick.map((c) => api.name(c)).join(" + ")));
    const ids = [...new Set(pick.flatMap((c) => { const s = setOf(cur, c); return s ? s.ids : [c]; }))].sort((a, b) => a - b);
    remember(cur, ids, `merge into “${name}”`);
    const names = namesOf(cur);
    for (const s of olds) setName(names, s.ids[0], s.prev);        // a set merged again: its head is a cluster again first
    const head = ids[0], prev = Object.prototype.hasOwnProperty.call(names, head) ? names[head] : null;
    names[head] = name;
    const keep = sets.filter((s) => !olds.includes(s));
    keep.push({ ids, prev });
    keep.sort((a, b) => a.ids[0] - b.ids[0]);
    cur.U.sets[cur.ck] = keep;
    for (const k of ids) cur.U.hidden.delete(k);
    if (cur.U.pick != null && ids.includes(cur.U.pick)) cur.U.pick = head;
    picked.clear(); $("csName").value = "";
    const n = pick.reduce((a, c) => a + z.get(c), 0);
    api.changed(cur.ck);
    msg(`Merged ${pick.length} into “${name}”: ${ids.length} clusters, ${n.toLocaleString()} ${W().bins}.`);
  }
  function unmerge(head) {
    const cur = now(); if (!cur) return;
    const s = setOf(cur, head); if (!s) return;
    const name = api.name(head);
    remember(cur, s.ids, `un-merge “${name}”`);
    setName(namesOf(cur), head, s.prev);
    const keep = setsOf(cur).filter((x) => x !== s);
    if (keep.length) cur.U.sets[cur.ck] = keep; else delete cur.U.sets[cur.ck];
    api.changed(cur.ck);
    msg(`“${name}” is ${s.ids.length} clusters again.`);
  }
  function rename(head, t) {
    const cur = now(); if (!cur) return;
    const name = cut(String(t || "").replace(/\s+/g, " ").trim()), old = api.name(head);
    if (!name || name === old) { refreshNow(); return; }
    remember(cur, [head], `rename “${old}”`);
    namesOf(cur)[head] = name;
    api.changed(cur.ck);
    msg(`Renamed “${old}” to “${name}”.`);
  }
  function undoLast() {
    const cur = now(), u = undo.pop();
    if (!cur || !u || u.ck !== cur.ck) { undo = []; syncButtons(); return; }
    if (u.sets.length) cur.U.sets[cur.ck] = u.sets; else delete cur.U.sets[cur.ck];
    const names = namesOf(cur);
    for (const [k, v] of Object.entries(u.names)) setName(names, k, v);
    cur.U.hidden = new Set(u.hidden); cur.U.pick = u.pick;
    api.changed(cur.ck);
    msg(`Undone: ${u.what}.`);
  }
  function refreshNow() { sig = ""; refresh(); }

  // ---- the CSV: one row per bin of the layer -------------------------------------------------------------------
  const q = (x) => (/[",\n\r]/.test(x) ? `"${String(x).replace(/"/g, '""')}"` : String(x));
  const slug = (s) => String(s || "").normalize("NFKD").replace(/[^A-Za-z0-9.-]+/g, "_").replace(/^[_.]+|[_.]+$/g, "").slice(0, 40) || "clusters";
  function csvParts(cur) {
    const U = cur.U, raw = U.clRaw.get(cur.ck), cl = U.cl, pos = U.pos, N = cl.length, NONE = api.NONE;
    if (!raw || raw.length !== N) throw new Error("the clusters are not loaded");
    const own = new Map(), set = new Map();
    const ownQ = (c) => { if (!own.has(c)) own.set(c, q(ownName(cur, c))); return own.get(c); };
    const setQ = (c) => { if (!set.has(c)) set.set(c, q(api.name(c))); return set.get(c); };
    const parts = ["bin_row,x_px,y_px,cluster,cluster_name,set_name\n"];
    const STEP = 50000;
    for (let r0 = 0; r0 < N; r0 += STEP) {
      const lines = [];
      for (let r = r0; r < Math.min(N, r0 + STEP); r++) {
        const xy = pos && pos.length >= 2 * N ? `${pos[2 * r].toFixed(2)},${pos[2 * r + 1].toFixed(2)}` : ",";
        const c = raw[r];
        lines.push(c === NONE ? `${r},${xy},,,` : `${r},${xy},${c},${ownQ(c)},${setQ(cl[r])}`);
      }
      parts.push(lines.join("\n") + "\n");
    }
    return parts;
  }
  function exportCsv() {
    const cur = now();
    if (!cur) { msg("Nothing to export: show a clustering first.", true); return; }
    try {
      const parts = csvParts(cur), name = `${slug(api.dsId())}_${slug(api.layer())}_cell_sets_${slug(api.label(cur.c))}.csv`;
      api.download(name, new Blob(parts, { type: "text/csv" }), "text/csv");
      const nSets = setsOf(cur).length;
      msg(`Saved ${name}: ${cur.U.cl.length.toLocaleString()} ${W().bins}, ${nSets ? `${nSets} ${nSets === 1 ? "set" : "sets"}` : "no sets yet (each cluster is its own set)"}.` +
          (api.touch ? " On a phone or tablet it goes to Downloads (iPhone and iPad: the Files app)." : ""));
    } catch (e) { msg(`The CSV could not be made: ${e && e.message ? e.message : e}`, true); }
  }

  // sets read from a saved JSON file: kept only when well formed (whole numbers, two or more, each cluster once)
  function clean(o) {
    const out = {};
    if (!o || typeof o !== "object" || Array.isArray(o)) return out;
    for (const [ck, l] of Object.entries(o)) {
      if (!Array.isArray(l)) continue;
      const used = new Set(), list = [];
      for (const x of l) {
        const ids = x && Array.isArray(x.ids) ? [...new Set(x.ids.filter((k) => Number.isInteger(k) && k >= 0 && k < 0xffff))].sort((a, b) => a - b) : [];
        if (ids.length < 2 || ids.some((k) => used.has(k))) continue;
        ids.forEach((k) => used.add(k));
        list.push({ ids, prev: typeof x.prev === "string" ? x.prev.slice(0, MAX_NAME) : null });
      }
      if (list.length) out[ck] = list;
    }
    return out;
  }

  window.SodaSets = {
    init(a) { api = a; buildDom(); refresh(); },
    refresh() { if (api) refresh(); },
    clean,
    // for tests from the console
    _pick(codes) { picked = new Set(codes); sig = ""; refresh(); },
    _merge(codes, name) { picked = new Set(codes); $("csName").value = name || ""; merge(); return $("csMsg").textContent; },
    _unmerge: (h) => { unmerge(h); return $("csMsg").textContent; },
    _rename: (h, t) => { rename(h, t); return $("csMsg").textContent; },
    _undo: () => { undoLast(); return $("csMsg").textContent; },
    _csv: () => { const cur = now(); return cur ? csvParts(cur).join("") : null; },
  };
})();
