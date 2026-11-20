// Your own results: a table of per-bin (per-cell, per-spot) results that a reader already has, loaded into the
// page from a file, with no re-export. A Seurat or scanpy UMAP, clusterings, cell types, samples, scores, QC.
//
//   The file   CSV or TSV (comma, tab or semicolon), optionally gzipped (.gz), one row per bin / cell / spot:
//              a barcode or cell id column and any of
//                two embedding columns   UMAP_1 / UMAP_2, umap_1, UMAP1, X_umap1, tSNE_1, PC_1 ...
//                categorical columns     seurat_clusters, leiden, louvain, cell_type, celltype, sample ...
//                                        (text with 2-200 values; whole numbers with up to 30, or up to 200 when
//                                        the name says clustering)
//                number columns          scores, QC (any other numeric column)
//                colour columns          <column>_color / <column>.colour holding #rrggbb per row (Seurat-style),
//                                        or a bare "color" next to exactly one other column
//              Roles are guessed from the names and the values, and can be corrected in a table of drop-downs.
//              A colour list in the order of the groups (scanpy's uns["<column>_colors"]) can be pasted there.
//
//   Matching   rows are matched to the layer's bins by barcode / cell id, through the export's id list (on by default)
//              (expr/<key>/ids.txt.gz, docs/DATA_FORMAT.md; `prepare.py ids` adds it to an export). Exact first,
//              then the Seurat / scanpy variants userclusters.py also accepts: a different or missing "-1"-style
//              ending ("-1", "_1", "-1-0", ".1") and a sample prefix ("S1_", "slice1_", "S1:"). A merged object
//              (several samples on the same bins) is refused until one sample is kept (starts with / ends with).
//              Without an id list: a soda_row (or bin_row, as SODA's own CSVs write) column of positions.bin
//              row numbers, or a first line "# soda row order" and exactly one row per bin, in the export's order;
//              both said with a warning. The match rate is reported; unmatched rows only as a count, and no id
//              is ever shown, logged or saved.
//
//   Where it goes  the UMAP: the section UMAP (Analyze tab, Clusters), in place of SODA's or as the only one, the
//              reader picks which; clusterings: the Clusters section's clusterings (tissue and UMAP colours, legend,
//              names, cell sets, volcano / DE and plots groups); numbers: View tab, Colour by a value.
//
//   Privacy    the file is read here and stays in this browser: nothing is uploaded. Save JSON can carry the
//              matched table (per positions.bin row, no ids, no file name; its size is shown first). A link to the
//              view cannot carry data: it carries a fingerprint of the file, and asks for the same file.
//
// viewer.js hands over what this needs (init) and reads the result back (current): clusterings in the format
// of umap/index.json entries (source "user", categories, colours) with their codes per positions.bin row
// (uint16, 0xffff = none), the embedding as positions.bin rows (NaN = not in it), numbers per row (NaN = none).
// Only unsigned typed arrays here: the de-identification check reads the signed names as sample IDs.
/* exported SodaOwn */
(function () {
  "use strict";
  let api = null;
  const $ = (id) => document.getElementById(id);
  const esc = (x) => String(x).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
  const NONE = 0xffff, NOROW = 0xffffffff;
  const MAX_FILE = 1024 * 1024 * 1024;           // bytes of the file as picked
  const MAX_TEXT = 512 * 1024 * 1024;            // bytes once unpacked
  const MAX_CATEGORIES = 200, MAX_INT_CATEGORIES = 30;
  const MIN_MATCH = 0.5;                         // below this share of the file's rows matched: not used until asked
  const SAMPLE_IDS = 400;                        // values of a column tried against the id list to find the id column
  const NA = new Set(["", "na", "nan", "none", "null", "<na>", "n/a", "#n/a"]);
  const SUFFIX = /(?:[-_.]\d{1,3})+$/;           // -1, _1, -1-0, -1_2, .1 (userclusters.py _SUFFIX)
  const SEPS = /[_:#|]/g;                        // where a sample prefix ends (userclusters.py _SEPS)
  const HEX = /^#?([0-9a-fA-F]{6})(?:[0-9a-fA-F]{2})?$/;
  const COLOR_COL = /^(.*?)[._]colou?rs?$/i;
  const BARE_COLOR = /^colou?rs?$/i;
  const CLUSTER_NAME = /cluster|leiden|louvain|seurat|ident|snn_res|_res\.|^res\.|celltype|cell_type|annot|type|label|group|sample|class|niche|domain|region|state/i;
  const ID_NAME = /^(|barcode|barcodes|cell|cells|cell_?id|cellid|cell_?barcode|spot|spot_?id|bin|bin_?id|obs_?names|index|_index|unnamed: ?0|rownames|row\.names)$/i;
  const ROW_NAME = /^(soda_row|bin_row)$/i;
  const ROW_ORDER = /^#\s*soda[\s_-]*row[\s_-]*order\b/i;
  // embedding columns: kind + axis. Spatial coordinates are not an embedding.
  const EMB = /^(?:x_)?(umap|tsne|t_sne|pca|pc|phate|mde|densmap|draw_graph_fa|fa|diffmap|harmony)[\s._-]*([0-9]|x|y)$/i;
  const EMB_RANK = { umap: 0, densmap: 0, tsne: 1, t_sne: 1, phate: 2, mde: 2, draw_graph_fa: 2, fa: 2, diffmap: 3, harmony: 4, pca: 5, pc: 5 };
  const ROLES = [["skip", "Not used"], ["id", "Barcode / cell id"], ["row", "Row number (SODA's order)"], ["x", "UMAP x"], ["y", "UMAP y"],
                 ["cat", "Groups (clustering, cell type …)"], ["num", "Value (a number)"]];

  // the results that are in use, per dataset and layer; the file read last (its table, for the role table)
  const store = new Map();          // `${ds}|${key}` -> R
  let D = null;                     // the draft: { ds, key, name, fp, T, roles, pal, prefix, suffix, m, R, applied, ... }
  let pending = null;               // a link's own block whose file is not loaded yet
  let busy = false;
  const idCache = new Map();        // `${ds}|${key}` -> Promise of the layer's ids

  const W = () => (api && api.words ? api.words() : { bin: "bin", bins: "bins", Bins: "Bins" });
  const n2 = (x) => Number(x).toLocaleString();
  const pct = (a, b) => (b ? `${((100 * a) / b).toFixed(a === b ? 0 : 1)}%` : "0%");
  const msg = (t, kind) => { const el = $("ownMsg"); if (!el) return; el.textContent = t || ""; el.classList.toggle("err", kind === "err"); el.classList.toggle("warn", kind === "warn"); };

  // ---- the section -----------------------------------------------------------------------------------
  function buildDom() {
    const sec = document.createElement("section");
    sec.className = "panel"; sec.id = "panelOwn"; sec.dataset.tab = "analyze"; sec.hidden = true;
    sec.innerHTML = `<h2>Your own results</h2>
      <div class="btnrow"><button id="ownLoadBtn" type="button" title="A CSV or TSV (also .gz) with a barcode or cell id column and your UMAP, clusters, cell types or scores">Load your own results…</button>
        <input id="ownFile" type="file" accept=".csv,.tsv,.txt,.gz,.csv.gz,.tsv.gz,text/csv,text/tab-separated-values,text/plain,application/gzip" hidden /></div>
      <div id="ownHint" class="meta own-hint">A CSV or TSV (also gzipped) with one row per bin: a barcode or cell id column, then your UMAP, clusters, cell types or scores. Or drop the file anywhere on the page. It is read in this browser and not uploaded.</div>
      <div id="ownMsg" class="meta own-msg" role="status"></div>
      <div id="ownBox" class="own-box" hidden>
        <div id="ownRoles" class="own-roles"></div>
        <div id="ownKeep" class="own-keep" hidden>
          <div class="meta">Several samples land on the same bins. Keep the rows of this section only:</div>
          <label class="row slim"><span>Ids starting with</span><input id="ownPrefix" type="text" autocomplete="off" spellcheck="false" /></label>
          <label class="row slim"><span>Ids ending with</span><input id="ownSuffix" type="text" autocomplete="off" spellcheck="false" /></label>
        </div>
        <label class="row slim" id="ownEmbRow" hidden><span>UMAP shown</span><select id="ownEmb" title="Which UMAP the Tissue / UMAP switch and Side by side show"><option value="file">From your file</option><option value="soda">SODA's (made at export)</option></select></label>
        <label class="row slim check" id="ownSaveRow" hidden><input id="ownSave" type="checkbox" /><span id="ownSaveTxt">Put the table in Save JSON</span></label>
        <div class="btnrow"><button id="ownUse" type="button" hidden>Use anyway</button><button id="ownClear" type="button" title="Forget the loaded results (the file itself is not touched)">Remove</button></div>
      </div>`;
    const side = $("sidebar");
    if (window.SodaLayout) { side.appendChild(sec); window.SodaLayout.adopt(sec); }
    else { const after = $("panelUmap"); if (after) after.after(sec); else side.appendChild(sec); }
    $("ownLoadBtn").onclick = () => $("ownFile").click();
    $("ownFile").onchange = (e) => { const f = e.target.files && e.target.files[0]; e.target.value = ""; if (f) loadFile(f); };
    $("ownClear").onclick = clear;
    $("ownUse").onclick = () => { if (D && D.R) { D.force = true; applyDraft(true); } };
    $("ownEmb").onchange = (e) => { const R = cur(); if (!R) return; R.useEmb = e.target.value !== "soda"; api.changed({ emb: true }); render(); };
    const keep = () => { if (!D) return; D.prefix = $("ownPrefix").value.trim(); D.suffix = $("ownSuffix").value.trim(); rematch(); };
    $("ownPrefix").onchange = keep; $("ownSuffix").onchange = keep;
    // drop a file anywhere on the page
    const drop = document.createElement("div");
    drop.id = "ownDrop"; drop.className = "own-drop"; drop.hidden = true;
    drop.innerHTML = `<div>Drop to load your own results<small>CSV or TSV, also .gz · read in this browser, not uploaded</small></div>`;
    document.body.appendChild(drop);
    const files = (e) => !!(e.dataTransfer && Array.from(e.dataTransfer.types || []).includes("Files"));
    let depth = 0;
    document.addEventListener("dragenter", (e) => { if (!files(e)) return; depth++; drop.hidden = false; });
    document.addEventListener("dragover", (e) => { if (!files(e)) return; e.preventDefault(); e.dataTransfer.dropEffect = "copy"; drop.hidden = false; });
    document.addEventListener("dragleave", (e) => { if (!files(e)) return; depth = Math.max(0, depth - 1); if (!depth) drop.hidden = true; });
    document.addEventListener("drop", (e) => {
      if (!files(e)) return;
      e.preventDefault(); depth = 0; drop.hidden = true;
      const f = e.dataTransfer.files && e.dataTransfer.files[0];
      if (f) loadFile(f, true);
    });
  }

  // ---- reading the file ---------------------------------------------------------------------------------
  async function gunzip(u8) {
    if (typeof DecompressionStream === "undefined") throw new Error("this browser cannot unpack .gz files; unzip it first, or use a current Chrome, Edge, Firefox or Safari");
    const rd = new Blob([u8]).stream().pipeThrough(new DecompressionStream("gzip")).getReader(), parts = [];
    let n = 0;
    for (;;) {
      const { done, value } = await rd.read();
      if (done) break;
      n += value.length;
      if (n > MAX_TEXT) { try { await rd.cancel(); } catch (e) { /* stopping is all that matters */ } throw new Error(`it unpacks to more than ${MAX_TEXT / 1048576} MB`); }
      parts.push(value);
    }
    const out = new Uint8Array(n);
    let at = 0;
    for (const p of parts) { out.set(p, at); at += p.length; }
    return out;
  }
  // what a file that is not a text table is, in words, or "" when it may be one
  function notTable(u8, name) {
    const h = (i) => u8[i];
    if (u8.length >= 4 && h(0) === 0x50 && h(1) === 0x4b && h(2) === 3 && h(3) === 4) return "This is a zip or Excel file, not a text table: save the sheet as CSV (comma or tab separated) and load that.";
    if (u8.length >= 4 && h(0) === 0x89 && h(1) === 0x48 && h(2) === 0x44 && h(3) === 0x46) return "This is an HDF5 file (.h5ad, .h5): the browser cannot read it. Write the columns to a CSV (see Advanced, \"Load your own results\"), or attach it at export with prepare.py export --user-h5ad.";
    if (u8.length >= 2 && ((h(0) === 0x58 && h(1) === 0x0a) || (h(0) === 0x52 && h(1) === 0x44 && h(2) === 0x58))) return "This is an R object (.rds / .RData): write the columns to a CSV with write.csv (see Advanced, \"Load your own results\") and load that.";
    if (/\.(h5ad|h5|loom|rds|rdata|xlsx|xls|zarr|parquet|feather)$/i.test(name || "")) return `A ${name.split(".").pop()} file cannot be read here: write the columns to a CSV or TSV (see Advanced, "Load your own results").`;
    const k = Math.min(u8.length, 65536);
    let ctl = 0;
    for (let i = 0; i < k; i++) { const b = u8[i]; if (b === 0) return "This is not a text table (it holds binary data). Load a CSV or TSV, optionally gzipped."; if (b < 9 || (b > 13 && b < 32)) ctl++; }
    if (ctl > k / 100) return "This is not a text table (it holds binary data). Load a CSV or TSV, optionally gzipped.";
    return "";
  }
  // a fingerprint of the unpacked bytes: two 32-bit FNV-style hashes and the length (not a security measure;
  // it tells "the same file" from "another file" when a link asks for one)
  function fingerprint(u8) {
    let a = 0x811c9dc5, b = 0x9747b28c;
    for (let i = 0; i < u8.length; i++) { const x = u8[i]; a = Math.imul(a ^ x, 0x01000193); b = Math.imul(b ^ x, 0x5bd1e995); b ^= b >>> 15; }
    const hx = (v) => (v >>> 0).toString(16).padStart(8, "0");
    return (hx(a) + hx(b ^ u8.length)).slice(0, 12);
  }

  // CSV / TSV into columns of strings. Leading "#" lines are comments ("# soda row order" is read).
  function parseTable(text) {
    if (text.charCodeAt(0) === 0xfeff) text = text.slice(1);
    let at = 0, rowOrder = false;
    while (at < text.length && text[at] === "#") {
      const nl = text.indexOf("\n", at), line = text.slice(at, nl < 0 ? text.length : nl);
      if (ROW_ORDER.test(line)) rowOrder = true;
      at = nl < 0 ? text.length : nl + 1;
    }
    const nl0 = text.indexOf("\n", at), first = text.slice(at, nl0 < 0 ? text.length : nl0);
    const count = (ch) => { let n = 0, q = false; for (const c of first) { if (c === '"') q = !q; else if (!q && c === ch) n++; } return n; };
    const tabs = count("\t"), commas = count(","), semis = count(";");
    const delim = tabs && tabs >= commas ? "\t" : semis > commas ? ";" : commas ? "," : tabs ? "\t" : ",";
    const rows = [];
    if (text.indexOf('"', at) < 0) {               // no quotes: plain splits
      const lines = text.slice(at).split(/\r?\n/);
      for (const ln of lines) if (ln.length && ln[0] !== "#") rows.push(ln.split(delim));
    } else {
      let row = [], field = "", q = false, i = at;
      const n = text.length, endRow = () => { row.push(field); field = ""; if (row.length > 1 || row[0] !== "") rows.push(row); row = []; };
      while (i < n) {
        const c = text[i];
        if (q) {
          if (c === '"') { if (text[i + 1] === '"') { field += '"'; i += 2; continue; } q = false; i++; continue; }
          const j = text.indexOf('"', i);
          if (j < 0) { field += text.slice(i); i = n; continue; }
          field += text.slice(i, j); i = j; continue;
        }
        if (c === '"') { q = true; i++; continue; }
        if (c === delim) { row.push(field); field = ""; i++; continue; }
        if (c === "\n" || c === "\r") { endRow(); if (c === "\r" && text[i + 1] === "\n") i++; i++; continue; }
        field += c; i++;
      }
      if (field !== "" || row.length) endRow();
    }
    if (!rows.length) return { error: "The file is empty." };
    let header = rows[0].map((h) => String(h).trim());
    const body = rows.slice(1);
    if (!body.length) return { error: "The file has a header but no rows." };
    // R's write.table without col.names = NA: the header has one name fewer than the rows (the row names)
    if (body[0].length === header.length + 1) header = ["", ...header];
    const nc = header.length;
    if (nc < 2 && !rowOrder) return { error: `Only one column was found${delim === "," ? "" : ` (separated by ${delim === "\t" ? "tabs" : "semicolons"})`}: the file needs a barcode or cell id column and at least one more.` };
    let ragged = 0;
    const cols = header.map(() => new Array(body.length));
    for (let r = 0; r < body.length; r++) {
      const f = body[r];
      if (f.length !== nc) ragged++;
      for (let j = 0; j < nc; j++) { const v = f[j]; cols[j][r] = v === undefined ? "" : v.trim(); }
    }
    if (ragged > body.length * 0.1) return { error: `${n2(ragged)} of ${n2(body.length)} rows have a different number of fields than the header: is this a CSV or TSV?` };
    // names: unique, never empty except the row-name column
    const seen = new Map();
    header = header.map((h, j) => {
      let s = h || (j === 0 ? "" : `column ${j + 1}`);
      if (seen.has(s)) { let k = 2; while (seen.has(`${s} (${k})`)) k++; s = `${s} (${k})`; }
      seen.set(s, 1); return s;
    });
    return { header, cols, nRows: body.length, delim, rowOrder, ragged };
  }

  const isNA = (s) => NA.has(s.toLowerCase());
  function num(s, comma) {
    if (isNA(s)) return NaN;
    const x = Number(comma ? s.replace(",", ".") : s);
    return Number.isFinite(x) ? x : NaN;
  }
  // what a column holds (all rows): numbers, whole numbers, colours, how many distinct values (capped)
  function stats(v, comma) {
    let na = 0, nn = 0, ints = 0, hex = 0;
    const seen = new Set();
    for (let i = 0; i < v.length; i++) {
      const s = v[i];
      if (isNA(s)) { na++; continue; }
      const x = num(s, comma);
      if (x === x) { nn++; if (Number.isInteger(x)) ints++; }
      else if (HEX.test(s)) hex++;
      if (seen.size <= MAX_CATEGORIES + 1) seen.add(s);
    }
    const on = v.length - na;
    return { on, numeric: on > 0 && nn >= on * 0.98, ints: on > 0 && ints === nn && nn >= on * 0.98, hex: on > 0 && hex >= on * 0.9, distinct: seen.size };
  }

  // ---- the layer's ids, and matching -------------------------------------------------------------------
  async function textOf(buf) {
    const u8 = new Uint8Array(buf);
    return new TextDecoder().decode(u8.length >= 2 && u8[0] === 0x1f && u8[1] === 0x8b ? await gunzip(u8) : u8);   // a host may have unpacked it
  }
  function layerIds(ctx) {
    const e = ctx.e;
    if (!e || !e.ids) return Promise.resolve(null);
    const k = `${ctx.ds}|${ctx.key}`;
    if (!idCache.has(k)) {
      idCache.set(k, (async () => {
        const ids = (await textOf(await ctx.src.buf(`expr/${ctx.key}/${e.ids.file || "ids.txt.gz"}`))).split("\n");
        if (ids.length && ids[ids.length - 1] === "") ids.pop();
        if (ids.length !== ctx.nBins) throw new Error(`the export's id list has ${n2(ids.length)} rows, the layer ${n2(ctx.nBins)}`);
        return ids;
      })());
      idCache.get(k).catch(() => idCache.delete(k));
    }
    return idCache.get(k);
  }
  // a matcher over the layer's ids: (id) -> [row, how, tag] (userclusters.py match_barcodes, without printing)
  //   how: 1 exact, 2 a different "-1"-style ending, 3 a sample prefix, 4 both; tag: what told the sample apart
  function matcher(ids) {
    const exact = new Map(), core = new Map(), clash = new Set();
    for (let i = 0; i < ids.length; i++) {
      const b = ids[i]; exact.set(b, i);
      const c = b.replace(SUFFIX, "");
      if (core.has(c)) clash.add(c); else core.set(c, i);
    }
    for (const c of clash) core.delete(c);           // a stripped form two ids share is no evidence
    return (u) => {
      let j = exact.get(u);
      if (j !== undefined) return [j, 1, ""];
      const cu = u.replace(SUFFIX, "");
      j = core.get(cu);
      if (j !== undefined) return [j, 2, u.slice(cu.length)];
      SEPS.lastIndex = 0;
      let m, k = 0;
      while ((m = SEPS.exec(u)) && k++ < 3) {
        const rest = u.slice(m.index + 1);
        j = exact.get(rest);
        if (j !== undefined) return [j, 3, u.slice(0, m.index + 1)];
        const cr = rest.replace(SUFFIX, "");
        j = core.get(cr);
        if (j !== undefined) return [j, 4, u.slice(0, m.index + 1) + "…" + rest.slice(cr.length)];
      }
      return [-1, 0, ""];
    };
  }

  // ---- roles: guessed from names and values ------------------------------------------------------------
  function guessRoles(T, ids, nBins) {
    const nc = T.header.length, comma = T.delim === ";";
    const st = T.cols.map((c) => stats(c, comma));
    const roles = new Array(nc).fill(null), notes = new Array(nc).fill("");
    // the id column: the one whose values match the layer's ids best (a sample of each), else by name
    let idj = -1;
    if (ids) {
      const M = matcher(ids);
      let best = 0;
      for (let j = 0; j < nc; j++) {
        if (st[j].numeric && st[j].distinct < 3) continue;
        const v = T.cols[j], step = Math.max(1, Math.floor(v.length / SAMPLE_IDS));
        let hit = 0;
        for (let i = 0; i < v.length; i += step) if (v[i] && M(v[i])[0] >= 0) hit++;
        if (hit > best || (hit === best && hit > 0 && ID_NAME.test(T.header[j]) && !ID_NAME.test(T.header[idj]))) { best = hit; idj = j; }
      }
    }
    const rowj = T.header.findIndex((h, j) => ROW_NAME.test(h) && st[j].ints);
    if (idj < 0 && rowj < 0 && !T.rowOrder) {            // "# soda row order": no id column at all
      idj = T.header.findIndex((h, j) => ID_NAME.test(h) && !st[j].numeric);
      if (idj < 0 && !st[0].numeric) idj = 0;
    }
    if (idj >= 0) roles[idj] = "id";
    if (rowj >= 0 && idj < 0) roles[rowj] = "row";
    else if (rowj >= 0) { roles[rowj] = "skip"; notes[rowj] = "row numbers: not needed, the ids match"; }
    // colour columns: <col>_color, <col>.colour ..., or a bare "color" beside one other column
    for (let j = 0; j < nc; j++) {
      if (roles[j]) continue;
      const m = COLOR_COL.exec(T.header[j]);
      let t = m ? T.header.indexOf(m[1]) : -1;
      if (t < 0 && BARE_COLOR.test(T.header[j])) { const others = T.header.map((_, k) => k).filter((k) => k !== j && roles[k] !== "id" && roles[k] !== "row"); if (others.length === 1) t = others[0]; }
      if (t >= 0 && t !== j && st[j].hex) roles[j] = `color:${t}`;
      else if (st[j].hex) { roles[j] = "skip"; notes[j] = "colours, but of which column? Pick it"; }
    }
    // an embedding: the best-ranked kind with two numeric axes
    const emb = new Map();
    for (let j = 0; j < nc; j++) {
      if (roles[j] || !st[j].numeric) continue;
      const m = EMB.exec(T.header[j]); if (!m) continue;
      const k = m[1].toLowerCase(), a = m[2].toLowerCase();
      if (!emb.has(k)) emb.set(k, []);
      emb.get(k).push({ j, a: a === "x" ? 0 : a === "y" ? 1 : Number(a) });
    }
    const kinds = [...emb.entries()].filter(([, l]) => l.length >= 2).sort((p, q) => (EMB_RANK[p[0]] ?? 9) - (EMB_RANK[q[0]] ?? 9));
    if (kinds.length) { const l = kinds[0][1].sort((p, q) => p.a - q.a); roles[l[0].j] = "x"; roles[l[1].j] = "y"; }
    // the rest: groups or values
    for (let j = 0; j < nc; j++) {
      if (roles[j]) continue;
      const s = st[j], named = CLUSTER_NAME.test(T.header[j]);
      if (s.distinct < 2) { roles[j] = "skip"; notes[j] = s.on ? "one value only" : "empty"; continue; }
      if (s.numeric) {
        if (s.ints && (s.distinct <= MAX_INT_CATEGORIES || (named && s.distinct <= MAX_CATEGORIES)) && !EMB.test(T.header[j])) roles[j] = "cat";
        else roles[j] = "num";
        continue;
      }
      if (s.distinct <= MAX_CATEGORIES) roles[j] = "cat";
      else { roles[j] = "skip"; notes[j] = `more than ${MAX_CATEGORIES} different values: an id or free text?`; }
    }
    return { roles, notes, st };
  }

  // file row -> positions.bin row, by the id column, a row-number column or the row order
  function matchRows(Dr, ids) {
    const T = Dr.T, nR = T.nRows, nB = Dr.nBins, rowOf = new Uint32Array(nR).fill(NOROW);
    const idj = Dr.roles.indexOf("id"), rowj = Dr.roles.indexOf("row");
    const m = { mode: "none", nFile: nR, considered: nR, matched: 0, other: 0, dup: 0, bad: 0, how: [0, 0, 0, 0, 0], merged: false };
    let tags = null;
    if (idj >= 0 && ids) {
      m.mode = "id";
      const M = matcher(ids), v = T.cols[idj], pre = Dr.prefix || "", suf = Dr.suffix || "";
      tags = new Array(nR);
      for (let i = 0; i < nR; i++) {
        let u = v[i];
        if (pre) { if (!u.startsWith(pre)) { m.other++; continue; } u = u.slice(pre.length); }
        if (suf) { if (!u.endsWith(suf)) { m.other++; continue; } u = u.slice(0, u.length - suf.length); }
        if (!u) continue;
        const [r, how, tag] = M(u);
        if (r >= 0) { rowOf[i] = r; m.how[how]++; tags[i] = tag; }
      }
      m.considered = nR - m.other;
    } else if (rowj >= 0) {
      m.mode = "row";
      const v = T.cols[rowj];
      for (let i = 0; i < nR; i++) { const x = Number(v[i]); if (Number.isInteger(x) && x >= 0 && x < nB) rowOf[i] = x; else m.bad++; }
    } else if (T.rowOrder && nR === nB) {
      m.mode = "order";
      for (let i = 0; i < nR; i++) rowOf[i] = i;
    } else if (T.rowOrder && idj < 0) m.mode = "badorder";
    else if (idj >= 0 && !ids) m.mode = "noids";
    // one bin claimed by several rows: the first is kept; with different sample tags the file is a merged object
    const owner = new Uint32Array(nB).fill(NOROW);
    let clashTags = 0;
    for (let i = 0; i < nR; i++) {
      const r = rowOf[i]; if (r === NOROW) continue;
      if (owner[r] === NOROW) { owner[r] = i; continue; }
      m.dup++;
      if (tags && tags[owner[r]] !== tags[i]) clashTags++;
      rowOf[i] = NOROW;
    }
    if (clashTags > 0 && clashTags >= m.dup * 0.5) m.merged = true;
    for (let i = 0; i < nR; i++) if (rowOf[i] !== NOROW) m.matched++;
    m.rate = m.considered ? m.matched / m.considered : 0;
    return { rowOf, m };
  }

  // ---- the result: what the viewer uses -------------------------------------------------------------------
  const safeId = (s, used) => {
    let b = String(s || "column").replace(/[^\w.-]+/g, "_").replace(/^_+|_+$/g, "").slice(0, 60) || "column", id = `own:${b}`, k = 2;
    while (used.has(id)) id = `own:${b}_${k++}`;
    used.add(id); return id;
  };
  // numeric-aware order, then code points (pandas' order for text categories)
  const catOrder = (a, b) => {
    const da = /^-?\d+$/.test(a), db = /^-?\d+$/.test(b);
    if (da && db) return Number(a) - Number(b);
    if (da !== db) return da ? -1 : 1;
    return a < b ? -1 : a > b ? 1 : 0;
  };
  function parseHexList(t) {
    return String(t || "").split(/[\s,;]+/).map((s) => s.replace(/^['"[]+|['"\]]+$/g, "")).filter(Boolean).map((s) => { const h = HEX.exec(s); return h ? `#${h[1].toLowerCase()}` : null; });
  }
  function build(Dr, rowOf, m) {
    const T = Dr.T, N = Dr.nBins, comma = T.delim === ";", used = new Set(), notes = [];
    const R = { ds: Dr.ds, key: Dr.key, fp: Dr.fp, nFile: T.nRows, matched: m.matched, m, emb: null, nEmb: 0, embName: "",
                useEmb: true, clusterings: [], nums: [], saved: false, ver: (Dr.R ? Dr.R.ver : 0) + 1 };
    const xj = Dr.roles.indexOf("x"), yj = Dr.roles.indexOf("y");
    if (xj >= 0 && yj >= 0) {
      const emb = new Float32Array(2 * N).fill(NaN), vx = T.cols[xj], vy = T.cols[yj];
      let k = 0;
      for (let i = 0; i < T.nRows; i++) {
        const r = rowOf[i]; if (r === NOROW) continue;
        const x = num(vx[i], comma), y = num(vy[i], comma);
        if (x === x && y === y) { emb[2 * r] = x; emb[2 * r + 1] = y; k++; }
      }
      if (k >= 3) { R.emb = emb; R.nEmb = k; R.embName = `${T.header[xj]} / ${T.header[yj]}`; }
      else notes.push(`The UMAP columns have fewer than 3 matched rows with numbers: no UMAP from the file.`);
    } else if (xj >= 0 || yj >= 0) notes.push("A UMAP needs both an x and a y column.");
    for (let j = 0; j < T.header.length; j++) {
      const role = Dr.roles[j];
      if (role === "cat") {
        const v = T.cols[j], numeric = Dr.st[j].numeric, seen = new Map();
        const key = (s) => (numeric ? String(num(s, comma)) : s);
        for (let i = 0; i < T.nRows; i++) { if (rowOf[i] === NOROW || isNA(v[i])) continue; const s = key(v[i]); seen.set(s, (seen.get(s) || 0) + 1); }
        const cats = [...seen.keys()].sort(catOrder);
        if (cats.length < 2) { notes.push(`${T.header[j]}: ${cats.length} group among the matched rows: left out.`); continue; }
        if (cats.length > MAX_CATEGORIES) { notes.push(`${T.header[j]}: ${cats.length} groups (more than ${MAX_CATEGORIES}): left out.`); continue; }
        const pos = new Map(cats.map((c, i) => [c, i])), codes = new Uint16Array(N).fill(NONE), sizes = new Array(cats.length).fill(0);
        for (let i = 0; i < T.nRows; i++) {
          const r = rowOf[i]; if (r === NOROW || isNA(v[i])) continue;
          const c = pos.get(key(v[i])); codes[r] = c; sizes[c]++;
        }
        // colours: a colour column (first colour per group), or a list pasted in the groups' order
        let colors = null;
        const cj = Dr.roles.indexOf(`color:${j}`), list = Dr.pal[j] ? parseHexList(Dr.pal[j]) : null;
        if (list && list.length) {
          if (list.length >= cats.length && list.slice(0, cats.length).every(Boolean)) colors = list.slice(0, cats.length);
          else notes.push(`${T.header[j]}: the colour list has ${list.filter(Boolean).length} colours for ${cats.length} groups: the usual colours are used.`);
        } else if (cj >= 0) {
          const first = new Map(), cv = T.cols[cj];
          for (let i = 0; i < T.nRows; i++) {
            if (rowOf[i] === NOROW || isNA(v[i])) continue;
            const c = key(v[i]); if (first.has(c)) continue;
            const h = HEX.exec(cv[i]); if (h) first.set(c, `#${h[1].toLowerCase()}`);
          }
          if (cats.every((c) => first.has(c))) colors = cats.map((c) => first.get(c));
          else notes.push(`${T.header[j]}: colours for ${first.size} of ${cats.length} groups: the usual colours are used.`);
        }
        R.clusterings.push({ id: safeId(T.header[j] || "groups", used), name: `${T.header[j] || "groups"} (your file)`, col: T.header[j], source: "user", own: true,
                             n: cats.length, sizes, categories: cats, ...(colors ? { colors } : {}), codes });
      } else if (role === "num") {
        const v = T.cols[j], out = new Float32Array(N).fill(NaN);
        let k = 0;
        for (let i = 0; i < T.nRows; i++) { const r = rowOf[i]; if (r === NOROW) continue; const x = num(v[i], comma); if (x === x) { out[r] = x; k++; } }
        if (!k) { notes.push(`${T.header[j]}: no numbers among the matched rows: left out.`); continue; }
        R.nums.push({ id: safeId(T.header[j] || "value", used), name: T.header[j] || "value", v: out, n: k });
      }
    }
    R.meta = { embedding: R.emb ? "own" : "none", clusterings: [], nEmbedded: R.nEmb, own: true };
    R.notes = notes;
    return R;
  }
  const usable = (R) => !!(R && (R.emb || R.clusterings.length || R.nums.length));

  // ---- loading, matching, applying ---------------------------------------------------------------------
  async function loadFile(file, dropped) {
    const ctx = api.ctx();
    if (window.SodaLayout && $("panelOwn") && !$("panelOwn").hidden) window.SodaLayout.show("panelOwn");
    if (!ctx) { api.status("Open a spatial dataset first: your own results are matched to its bins.", true); setTimeout(() => api.status(null), 5000); return; }
    if (busy) { msg("Still reading the last file…", "warn"); return; }
    if (/\.json$/i.test(file.name || "") && dropped) { msg("A .json file is a saved view: load it with Load JSON (Analyze tab, Draw a region). Your own results come as a CSV or TSV.", "err"); return; }
    if (file.size > MAX_FILE) { msg(`This file is ${(file.size / 1073741824).toFixed(1)} GB: too big to read in a browser. Keep only the columns you need (and gzip it).`, "err"); return; }
    busy = true;
    msg(`Reading ${file.name || "the file"}…`);
    try {
      let u8 = new Uint8Array(await file.arrayBuffer());
      if (u8.length >= 2 && u8[0] === 0x1f && u8[1] === 0x8b) {
        try { u8 = await gunzip(u8); } catch (e) { throw new Error(`Could not unpack this .gz file: ${e && e.message ? e.message : e}`); }
      }
      if (u8.length > MAX_TEXT) throw new Error(`This file is more than ${MAX_TEXT / 1048576} MB of text: too big to read in a browser. Keep only the columns you need.`);
      const bad = notTable(u8, file.name);
      if (bad) throw new Error(bad);
      const fp = fingerprint(u8);
      const T = parseTable(new TextDecoder().decode(u8));
      u8 = null;
      if (T.error) throw new Error(`${T.error} Expected: a CSV or TSV with a header row.`);
      let ids = null, idErr = "";
      try { ids = await layerIds(ctx); } catch (e) { idErr = e && e.message ? e.message : String(e); }
      const now = api.ctx();
      if (!now || now.ds !== ctx.ds || now.key !== ctx.key) { msg("The dataset or layer changed while the file was read: load it again.", "warn"); return; }
      const g = guessRoles(T, ids, ctx.nBins);
      D = { ds: ctx.ds, key: ctx.key, label: ctx.label, nBins: ctx.nBins, hasIds: !!ids, idErr, name: file.name || "your file", fp, T,
            roles: g.roles, notes: g.notes, st: g.st, pal: {}, prefix: "", suffix: "", force: false, R: null, applied: false };
      await rematch(ids, true);
    } catch (e) {
      render();                                     // the results in use (if any) stay as they were
      const inUse = D && D.applied && api.ctx() && D.ds === api.ctx().ds && D.key === api.ctx().key ? ` The results in use are still those of ${D.name}.` : "";
      msg(`${file.name ? `${file.name}: ` : ""}${e && e.message ? e.message : String(e)}${inUse}`, "err");
    } finally { busy = false; }
  }
  // match again (the id column or the sample filter changed) and build; first: the file was just read
  async function rematch(ids, first) {
    if (!D) return;
    if (ids === undefined) { try { ids = await layerIds(api.ctx() || D); } catch (e) { ids = null; } }
    const { rowOf, m } = matchRows(D, ids);
    D.rowOf = rowOf; D.m = m;
    rebuild(first);
  }
  function rebuild(first) {
    if (!D) return;
    D.R = D.m.matched ? build(D, D.rowOf, D.m) : null;
    const ok = D.R && usable(D.R) && !D.m.merged && (D.force || D.m.rate >= MIN_MATCH || D.m.mode === "row" || D.m.mode === "order");
    if (ok && (first || D.applied || D.force)) applyDraft(first);
    else {
      if (D.applied && !ok) { store.delete(`${D.ds}|${D.key}`); D.applied = false; api.changed({}); }
      render();
    }
  }
  function applyDraft(first) {
    if (!D || !D.R) return;
    const R = D.R, k = `${D.ds}|${D.key}`, prev = store.get(k);
    if (prev) { R.useEmb = prev.useEmb; R.saveOn = prev.saveOn; }
    store.set(k, R);
    D.applied = true;
    // show what came: the first clustering (when nothing of the file was on show), else the one on show stays
    const onShow = api.clusterKey(), wasOwn = typeof onShow === "string" && onShow.startsWith("own:");
    const show = R.clusterings.length && (first || (wasOwn && !R.clusterings.some((c) => c.id === onShow))) ? R.clusterings[0].id : null;
    Promise.resolve(api.changed({ show, emb: true })).then(() => {
      if (pending && pending.fp === R.fp) applyPending(R);
      render();
    });
    render();
  }
  function clear() {
    if (D) { store.delete(`${D.ds}|${D.key}`); D = null; }
    const ctx = api.ctx();
    if (ctx) store.delete(`${ctx.ds}|${ctx.key}`);
    msg("");
    api.changed({});
    render();
  }

  // ---- the table of roles ------------------------------------------------------------------------------
  function roleOptions(j) {
    const T = D.T, cur = D.roles[j], st = D.st[j];
    const opts = ROLES.filter(([v]) => v !== "row" || st.ints || cur === "row").map(([v, t]) => `<option value="${v}"${cur === v ? " selected" : ""}>${esc(t)}</option>`);
    const targets = D.roles.map((r, k) => (r === "cat" && k !== j ? k : -1)).filter((k) => k >= 0);
    if (st.hex || (cur || "").startsWith("color:")) for (const k of targets) opts.push(`<option value="color:${k}"${cur === `color:${k}` ? " selected" : ""}>Colours of ${esc(T.header[k] || "(first column)")}</option>`);
    return opts.join("");
  }
  function matchText() {
    const m = D.m, w = W(), lay = D.label || D.key;
    if (m.mode === "noids") return { kind: "err", t: `This export has no id list for the ${lay} layer, so the rows cannot be matched by barcode or cell id. ` +
      `Ask whoever made the export to add it: prepare.py ids --id ${D.ds} --out <cache folder> --dest <export folder> (an export holds one by default, unless it was made with --no-ids or the list failed its sample-ID check). ` +
      `Or give the file a soda_row column (each row's ${w.bin} number in this export, as SODA's own CSVs write it), or a first line "# soda row order" ` +
      `and exactly ${n2(D.nBins)} rows in the export's order.${D.idErr ? ` (${D.idErr})` : ""}` };
    if (m.mode === "badorder") return { kind: "err", t: `The file says it is in SODA's row order ("# soda row order") but has ${n2(m.nFile)} rows; the ${lay} layer has ${n2(D.nBins)} ${w.bins}. Row order only works with exactly one row per ${w.bin}: not used.` };
    if (m.mode === "none" || (m.mode === "id" && !m.matched && !m.other)) return { kind: "err", t: D.hasIds ? `No column of this file matches the ${lay} layer's barcodes or cell ids. Pick the id column below, or check that the file belongs to this dataset and layer.`
      : `Nothing to match the rows with: pick the barcode / cell id column, or a row-number column, below.` };
    let t = `${n2(m.matched)} of ${n2(m.considered)} rows matched ${w.bins} of the ${lay} layer (${pct(m.matched, m.considered)})`;
    if (m.mode === "id") {
      const alt = m.how[2] + m.how[3] + m.how[4];
      t += alt ? `: ${n2(m.how[1])} by exact id, ${n2(alt)} after setting aside an ending like "-1" or a sample prefix.` : ", by exact id.";
    } else t += ".";
    if (m.other) t += ` ${n2(m.other)} rows of other samples set aside.`;
    if (m.considered > m.matched) t += ` ${n2(m.considered - m.matched)} rows matched no ${w.bin} (not listed).`;
    if (m.dup) t += ` ${n2(m.dup)} rows matched a ${w.bin} another row already had: the first one is used.`;
    if (m.bad) t += ` ${n2(m.bad)} rows have no valid row number.`;
    t += ` ${n2(m.matched)} of the layer's ${n2(D.nBins)} ${w.bins} have results.`;
    let kind = "";
    if (m.mode === "row") { t = `Matched by row number (the soda_row column), not by barcode: right only if the numbers come from this export's own ${w.bins}. ` + t; kind = "warn"; }
    if (m.mode === "order") { t = `Matched by row order, not by barcode: right only if the file was written from this export's ${w.bins}, in their order. ` + t; kind = "warn"; }
    if (m.merged) { t = `This file holds several samples: ${n2(m.dup)} rows land on ${w.bins} that rows with another prefix or ending already took. Keep one sample's rows below.`; kind = "err"; }
    else if (m.rate < MIN_MATCH && m.mode === "id") { t += ` That is under half of the file: is it the right dataset and layer?${D.applied ? "" : " Use anyway puts the matched rows to use."}`; kind = "warn"; }
    return { kind, t };
  }
  function render() {
    const box = $("ownBox"); if (!box) return;
    const ctx = api.ctx(), here = ctx && D && D.ds === ctx.ds && D.key === ctx.key;
    const R = cur();
    if (!here) {
      box.hidden = true;
      const other = ctx ? [...store.values()].find((x) => x.ds === ctx.ds && x.key !== ctx.key) : null;
      if (!busy) msg(other ? `Your results are for the ${other.key} layer: pick that layer to see them.` : pending ? pendingText() : "", other ? "" : pending ? "warn" : "");
      if (R) renderSaved(R);
      return;
    }
    box.hidden = false;
    const T = D.T, mt = matchText(), notes = [...(D.R && D.R.notes ? D.R.notes : [])];
    if (D.applied && D.R) {
      const parts = [];
      if (D.R.emb) parts.push(`a UMAP (${D.R.embName}, ${n2(D.R.nEmb)} ${W().bins})`);
      if (D.R.clusterings.length) parts.push(`${D.R.clusterings.length} ${D.R.clusterings.length === 1 ? "grouping" : "groupings"} (Analyze tab, Clusters: Colour bins by, Clustering)`);
      if (D.R.nums.length) parts.push(`${D.R.nums.length} ${D.R.nums.length === 1 ? "value" : "values"} (View tab, More options: Colour by a value)`);
      notes.unshift(`In use: ${parts.join(", ")}.`);
    } else if (D.R && !usable(D.R)) notes.unshift("No column is in use: pick a role for at least one column.");
    msg(`${D.name}: ${T.nRows.toLocaleString()} rows, ${T.header.length} columns. ${mt.t}${notes.length ? " " + notes.join(" ") : ""}`, mt.kind);
    const host = $("ownRoles");
    host.innerHTML = `<div class="own-role own-head"><span>Column</span><span>Use as</span></div>` + T.header.map((h, j) => {
      const note = D.notes[j] && D.roles[j] === "skip" ? `<small>${esc(D.notes[j])}</small>` : "";
      const pal = D.roles[j] === "cat" ? `<input class="own-pal" data-j="${j}" type="text" spellcheck="false" autocomplete="off" placeholder="Colours in group order (optional): #1f77b4, #ff7f0e …" value="${esc(D.pal[j] || "")}" />` : "";
      return `<div class="own-role"><span class="own-col" title="${esc(h || "(row names)")}">${esc(h || "(row names)")}${note}</span>` +
             `<select data-j="${j}" aria-label="Use of the column ${esc(h || "row names")}">${roleOptions(j)}</select>${pal}</div>`;
    }).join("");
    for (const sel of host.querySelectorAll("select[data-j]")) sel.onchange = (e) => setRole(Number(e.target.dataset.j), e.target.value);
    for (const inp of host.querySelectorAll("input.own-pal")) inp.onchange = (e) => { D.pal[Number(e.target.dataset.j)] = e.target.value; rebuild(); };
    $("ownKeep").hidden = !(D.m && (D.m.merged || D.prefix || D.suffix));
    $("ownPrefix").value = D.prefix; $("ownSuffix").value = D.suffix;
    $("ownUse").hidden = !(D.R && usable(D.R) && !D.applied && !D.m.merged);
    renderSaved(D.applied ? R : null);
  }
  // the parts that hold for results in use, also those put back from Save JSON (no table to change)
  function renderSaved(R) {
    const sodaEmb = api.sodaUmap();
    $("ownEmbRow").hidden = !(R && R.emb && sodaEmb);
    if (R) $("ownEmb").value = R.useEmb ? "file" : "soda";
    $("ownSaveRow").hidden = !R;
    if (R) {
      if (R.saveOn === undefined) R.saveOn = false;
      R.saveOn = $("ownSave").checked = R.saveOn || false;
      $("ownSave").onchange = () => { R.saveOn = $("ownSave").checked; };
      $("ownSaveTxt").textContent = `Put the table in Save JSON (adds about ${sizeText(saveSize(R))}: the matched rows' results, no ids)`;
    }
    if (R && R.saved && !(D && D.applied && D.R === R)) {
      $("ownBox").hidden = false; $("ownRoles").innerHTML = ""; $("ownKeep").hidden = true; $("ownUse").hidden = true;
      msg(`Your own results, put back from a saved file: ${[R.emb ? "a UMAP" : "", R.clusterings.length ? `${R.clusterings.length} grouping(s)` : "", R.nums.length ? `${R.nums.length} value(s)` : ""].filter(Boolean).join(", ")} for ${n2(R.matched)} ${W().bins}. Load the file again to change the columns' roles.`);
    }
  }
  function setRole(j, v) {
    if (!D) return;
    const was = D.roles[j];
    // one id column, one row-number column, one x, one y: picking it here takes it from elsewhere
    if (["id", "row", "x", "y"].includes(v)) D.roles = D.roles.map((r, k) => (k !== j && r === v ? "skip" : r));
    D.roles[j] = v;
    if (was === "cat" && v !== "cat") D.roles = D.roles.map((r) => (r === `color:${j}` ? "skip" : r));
    if (["id", "row"].includes(v) || ["id", "row"].includes(was)) { rematch(); return; }
    rebuild();
  }
  const sizeText = (b) => (b < 1048576 ? `${Math.max(1, Math.round(b / 1024))} kB` : `${(b / 1048576).toFixed(1)} MB`);

  // ---- what the viewer reads -----------------------------------------------------------------------------
  function cur() { const c = api && api.ctx(); return c ? store.get(`${c.ds}|${c.key}`) || null : null; }
  function current(ds, key) { return store.get(`${ds}|${key}`) || null; }
  // Colour by a value (binpaint.js): the file's numbers
  function paintCols(ds, key) {
    const R = store.get(`${ds}|${key}`);
    return R ? R.nums.map((x) => ({ id: x.id, label: `${x.name} (your file)`, unit: "", scale: "lin", load: async () => x.v })) : [];
  }

  // ---- Save JSON: the matched table, per positions.bin row (no ids, no file name) --------------------------
  const r6 = (x) => +x.toPrecision(6);
  function saveBlock(R) {
    R = R || cur();
    if (!R || !R.saveOn) return null;
    if (R.block) return { ...R.block, useEmb: R.useEmb };
    const N = R.clusterings.length ? R.clusterings[0].codes.length : R.emb ? R.emb.length / 2 : R.nums.length ? R.nums[0].v.length : 0;
    const has = new Uint8Array(N);
    if (R.emb) for (let r = 0; r < N; r++) if (R.emb[2 * r] === R.emb[2 * r]) has[r] = 1;
    for (const c of R.clusterings) for (let r = 0; r < N; r++) if (c.codes[r] !== NONE) has[r] = 1;
    for (const x of R.nums) for (let r = 0; r < N; r++) if (x.v[r] === x.v[r]) has[r] = 1;
    const rows = [];
    for (let r = 0; r < N; r++) if (has[r]) rows.push(r);
    const b = { v: 1, layer: R.key, nBins: N, fp: R.fp, nFile: R.nFile, matched: R.matched, useEmb: R.useEmb, rows };
    if (R.emb) b.emb = { name: R.embName, x: rows.map((r) => (R.emb[2 * r] === R.emb[2 * r] ? r6(R.emb[2 * r]) : null)), y: rows.map((r) => (R.emb[2 * r + 1] === R.emb[2 * r + 1] ? r6(R.emb[2 * r + 1]) : null)) };
    b.groups = R.clusterings.map((c) => ({ id: c.id, name: c.col, categories: c.categories, ...(c.colors ? { colors: c.colors } : {}), codes: rows.map((r) => (c.codes[r] === NONE ? null : c.codes[r])) }));
    b.values = R.nums.map((x) => ({ id: x.id, name: x.name, v: rows.map((r) => (x.v[r] === x.v[r] ? r6(x.v[r]) : null)) }));
    R.block = b;
    return { ...b, useEmb: R.useEmb };
  }
  function saveSize(R) {
    if (R.blockBytes === undefined) { const on = R.saveOn; R.saveOn = true; const b = saveBlock(R); R.saveOn = on; R.blockBytes = b ? JSON.stringify(b).length : 0; }
    return R.blockBytes;
  }
  // Load JSON: put a saved table back (viewer.js sessionApply, after the layer switch). -> notes
  async function restore(b) {
    const notes = [], ctx = api.ctx();
    const isObj = (x) => !!x && typeof x === "object" && !Array.isArray(x);
    if (!isObj(b) || !ctx) return notes;
    const bad = () => { notes.push("The file's own results could not be read: left out."); return notes; };
    if (b.v !== 1 || !Array.isArray(b.rows)) return bad();
    if (b.layer !== ctx.key) { notes.push(`The file's own results are for the ${String(b.layer).slice(0, 40)} layer, not this one: left out.`); return notes; }
    const N = ctx.nBins, rows = b.rows, n = rows.length;
    if (b.nBins !== N || !rows.every((r, i) => Number.isInteger(r) && r >= 0 && r < N && (i === 0 || r > rows[i - 1]))) return bad();
    const str = (x, m) => String(x == null ? "" : x).replace(/[\u0000-\u001f<>"'`]/g, "").slice(0, m || 120);
    const used = new Set();
    const R = { ds: ctx.ds, key: ctx.key, fp: /^[0-9a-f]{6,16}$/.test(b.fp) ? b.fp : "", nFile: Number.isInteger(b.nFile) ? b.nFile : n, matched: Number.isInteger(b.matched) ? b.matched : n,
                m: null, emb: null, nEmb: 0, embName: "", useEmb: b.useEmb !== false, clusterings: [], nums: [], saved: true, saveOn: true, ver: 1, notes: [] };
    const fin = (x) => typeof x === "number" && Number.isFinite(x);
    if (isObj(b.emb) && Array.isArray(b.emb.x) && Array.isArray(b.emb.y) && b.emb.x.length === n && b.emb.y.length === n) {
      const emb = new Float32Array(2 * N).fill(NaN);
      let k = 0;
      for (let i = 0; i < n; i++) if (fin(b.emb.x[i]) && fin(b.emb.y[i])) { emb[2 * rows[i]] = b.emb.x[i]; emb[2 * rows[i] + 1] = b.emb.y[i]; k++; }
      if (k >= 3) { R.emb = emb; R.nEmb = k; R.embName = str(b.emb.name, 120); }
    }
    for (const g of Array.isArray(b.groups) ? b.groups.slice(0, 100) : []) {
      if (!isObj(g) || !Array.isArray(g.categories) || !Array.isArray(g.codes) || g.codes.length !== n) return bad();
      const cats = g.categories.slice(0, MAX_CATEGORIES).map((c) => str(c, 120));
      if (cats.length < 2) continue;
      const codes = new Uint16Array(N).fill(NONE), sizes = new Array(cats.length).fill(0);
      for (let i = 0; i < n; i++) { const c = g.codes[i]; if (Number.isInteger(c) && c >= 0 && c < cats.length) { codes[rows[i]] = c; sizes[c]++; } }
      const colors = Array.isArray(g.colors) && g.colors.length === cats.length && g.colors.every((h) => typeof h === "string" && /^#[0-9a-f]{6}$/i.test(h)) ? g.colors.map((h) => h.toLowerCase()) : null;
      const col = str(g.name, 120) || "groups";
      R.clusterings.push({ id: safeId(col, used), name: `${col} (your file)`, col, source: "user", own: true, n: cats.length, sizes, categories: cats, ...(colors ? { colors } : {}), codes });
    }
    for (const x of Array.isArray(b.values) ? b.values.slice(0, 200) : []) {
      if (!isObj(x) || !Array.isArray(x.v) || x.v.length !== n) return bad();
      const v = new Float32Array(N).fill(NaN);
      let k = 0;
      for (let i = 0; i < n; i++) if (fin(x.v[i])) { v[rows[i]] = x.v[i]; k++; }
      const name = str(x.name, 120) || "value";
      if (k) R.nums.push({ id: safeId(name, used), name, v, n: k });
    }
    if (!usable(R)) return bad();
    R.meta = { embedding: R.emb ? "own" : "none", clusterings: [], nEmbedded: R.nEmb, own: true };
    store.set(`${ctx.ds}|${ctx.key}`, R);
    if (D && D.ds === ctx.ds && D.key === ctx.key) D = null;
    pending = null;
    await api.changed({ show: null, emb: true });
    render();
    return notes;
  }

  // ---- a link to the view: the file's fingerprint, and what of it was on show ------------------------------
  function capture() {
    const R = cur(); if (!R || !R.fp) return null;
    const o = { fp: R.fp, n: R.nFile, m: R.matched };
    const ck = api.clusterKey();
    if (typeof ck === "string" && R.clusterings.some((c) => c.id === ck)) { o.res = ck; o.cb = api.colorBy() === "cluster" ? 1 : 0; }
    const pc = $("paintMode") && $("paintMode").value === "value" && $("paintCol") ? $("paintCol").value : "";
    if (pc && R.nums.some((x) => x.id === pc)) o.pc = pc;
    if (R.emb && !R.useEmb) o.e = 0;
    return o;
  }
  function clean(x, h) {
    if (!h.isObj(x) || typeof x.fp !== "string" || !/^[0-9a-f]{6,16}$/.test(x.fp)) return null;
    const o = { fp: x.fp };
    h.put(o, "n", h.int(x.n, 0, 1e9)); h.put(o, "m", h.int(x.m, 0, 1e9));
    h.put(o, "res", h.str(x.res, 80)); h.put(o, "pc", h.str(x.pc, 80)); h.put(o, "cb", h.int(x.cb, 0, 1)); h.put(o, "e", h.int(x.e, 0, 1));
    return o;
  }
  const pendingText = () => `This view was made with your own results from a file (fingerprint ${pending.fp}${pending.n ? `, ${n2(pending.n)} rows` : ""}). ` +
    "Load the same file here to see them; it stays in this browser.";
  async function apply(x) {
    if (!x) { pending = null; render(); return []; }       // a view without own results leaves loaded ones as they are
    const R = cur();
    if (R && R.fp === x.fp) { pending = null; await applyPending(R, x); return []; }
    pending = x;
    render();
    return [`${pendingText()} (Analyze tab, More options: Your own results.)`];
  }
  async function applyPending(R, x) {
    x = x || pending; pending = null;
    if (!x) return;
    if (x.e === 0 && R.emb && R.useEmb) { R.useEmb = false; await api.changed({ emb: true }); }
    if (x.res && R.clusterings.some((c) => c.id === x.res)) await api.showClustering(x.res, x.cb !== 0);
    if (x.pc && R.nums.some((v) => v.id === x.pc) && $("paintMode") && $("paintCol")) {
      const mode = $("paintMode");
      mode.value = "value"; mode.dispatchEvent(new Event("change"));
      await new Promise((r) => setTimeout(r, 0));
      const col = $("paintCol"); col.value = x.pc; col.dispatchEvent(new Event("change"));
    }
    if (R === cur()) msg("This is the file the view was made with: its colouring is back.");
  }

  // ---- the viewer's calls ----------------------------------------------------------------------------
  function onLayer() {
    const ctx = api && api.ctx();
    const sec = $("panelOwn");
    if (sec) sec.hidden = !ctx;
    render();
  }
  function reset() {
    store.clear(); D = null; pending = null; idCache.clear();
    if ($("panelOwn")) $("panelOwn").hidden = true;
    msg("");
  }

  window.SodaOwn = {
    init(a) {
      api = a;
      buildDom();
      if (window.SodaViewLink && window.SodaViewLink.extend) window.SodaViewLink.extend({ key: "own", capture, clean, apply });
    },
    onLayer, reset, current, paintCols, saveBlock: () => saveBlock(), restore, load: loadFile,
    // tests
    _state: () => ({ D, pending, busy, store: [...store.keys()] }),
    _parse: parseTable, _match: (ids) => matcher(ids),
  };
})();
