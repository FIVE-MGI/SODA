// The control layout: three tabs on a rail (View, Analyze, Share), sections that fold, "More options", "?" notes, a top
// bar over the map, and on phones a bottom bar with a sheet.
//
// Nothing here makes or replaces a control. index.html holds every control with its id, exactly
// as viewer.js and the other scripts expect; each section.panel names its tab in data-tab. This
// file MOVES those nodes: a section into its tab's box, its heading into a button, the rest of it
// into a body that folds, a few rows into a "More options" box, the long notes behind a "?".
// It runs before viewer.js (it needs only the markup); panels that other scripts build later
// (volcano.js, plots.js) are handed over with SodaLayout.adopt().
//
// What decides which sections show is unchanged: the `hidden` attribute, set by viewer.js, qc.js,
// volcano.js, plots.js and roi.js. A tab with no visible section hides; the tab on show is re-chosen
// when a dataset has finished opening (viewer.js: SodaLayout.opening).
//
// Without this file the page is the plain sidebar: the rail and the top bar stay hidden.
/* exported SodaLayout */
(function () {
  "use strict";
  const $ = (id) => document.getElementById(id);
  const app = document.querySelector(".app"), side = $("sidebar"), rail = $("tabRail"), tabList = $("tabList"),
        panelsBox = $("tabPanels"), topBar = $("topBar");
  if (!app || !side || !rail || !tabList || !panelsBox || !topBar) return;

  const TABS = ["view", "analyze", "share"];
  const tabBtn = {}, tabBox = {};
  for (const t of TABS) {
    tabBtn[t] = tabList.querySelector(`[role=tab][data-tab=${t}]`);
    tabBox[t] = $(tabBtn[t] ? tabBtn[t].getAttribute("aria-controls") : "");
    if (!tabBtn[t] || !tabBox[t]) return;
  }

  // ---- what each section is: its tab comes from data-tab, the rest from here -------------------
  //   more   ids of controls whose row goes behind "More options" (moreAt: where the button sits)
  //   notes  ids of existing notes that go behind the "?"
  //   help   the "?" text, in the sequencing words (bins, genes); viewer.js's wording function
  //          turns them into squares / cells / ROIs and markers, as it does for the labels
  const SECTIONS = {
    panelChannels: {
      help: "The picture under everything else. Image opacity: how much of it shows; at 0 the image is off (the stain settings stay). Stains: tick one to show it, click its colour to change it, and drag the two handles under its histogram to set the darkest and the brightest value shown. The γ button opens two more settings for that stain: Gamma (below 1 brightens faint signal, above 1 darkens it; 1 changes nothing, a double-click puts it back) and Colours (its own colour, or grey, viridis or magma from dark to bright, best for one stain on its own). reset puts all of these back. Type in the filter to shorten a long list. An H&E image has no stains to pick: Image opacity is its one control." },
    panelGenes: { more: ["geneOpacity", "pointScale", "tissueBins", "cellNuclei"],
      help: "Pick the layer, then search a gene and press Enter (or click a hit) to draw it on the slide. Each gene on show gets a colour and a slider for the range of values drawn. View: Normal draws each bin as it is, Density adds up neighbours so busy areas stand out. Cell shapes (segmented cells that came with their boundaries): Outlines draws every cell's own outline once you zoom in, Dots keeps the dots. More options: opacity, point size, an outline of every bin in the tissue, and the nucleus outlines." },
    panelPaint: {
      help: "Colours every bin by one number, or by two genes, instead of by the genes picked above. A value per bin: pick the value (total counts, genes detected, % mitochondrial and whatever else the export measured per bin, such as a GeoMx ROI's Q3 normalisation factor). The colours run along viridis from the From value (dark) to the to value (yellow); Auto sets them to the 1st and 99th percentile, and Log spreads out values that pile up at the low end. Hover or tap a bin to read its value. Two genes: type a gene in A and in B and press Enter; a bin counts as having a gene when its value is above that gene's number (empty: above 0, that is found at all; for markers, above the median square). Both / one only colours each bin by which of the two it has (yellow both, red A only, green B only) and counts them; Blend adds red for A to green for B, brighter for more, so bins high in both turn yellow. Bins with neither are not drawn. Gene set score: type or paste a list of genes (spaces, commas or new lines between them) and press Score or Ctrl+Enter, or pick a cell type to use its marker genes (the dataset's own list first, then SODA's generic one). Each bin's score is the mean, over the listed genes the layer has, of its normalised value: log1p(10,000 × count / the bin's total count) for sequencing bins and cells, log1p(count / Q3 factor) for GeoMx ROIs, and for markers the marker scaled from its 1st to its 99th percentile square. It is a plain mean, not a module score with control genes. Genes the layer does not have are left out and named under the list. The colours and From / to work as for a value. The legend on the map and in Save image says what the colours mean. Colour bins by cluster (Clusters, Analyze tab) wins over this while it is on." },
    panelMolecules: { more: ["molSize", "molOpacity"],
      help: "Every decoded molecule of the genes picked here, one dot each, one colour per gene. Zoomed out only a share of the dots is drawn; zoom in to see them all. More options: dot size and opacity." },
    panelOverlays: { help: "Outlines that came with the dataset (regions, cells, ROIs), drawn over the picture. Tick one to show it; hover an outline to read its name." },
    panelRna: { more: ["rnaBin", "rnaPointSize"], moreAt: "end",
      help: "Each dot is a cell, placed by the UMAP that came with the data and coloured by its cluster. Add a gene to get one small UMAP per gene, coloured by how much of it each cell has. View: Normal draws the cells, Density a smoothed map. More options: bin size (Density) and point size." },
    panelUmap: { notes: ["umapNote"],
      help: "Colour the bins by the genes picked in View, by cluster, or by suggested cell type. Click a cluster in the list to light it in every view; the pencil renames it and the tick box hides it. Box and Freehand select bins by drawing around them, on the tissue or on the UMAP. What is selected or highlighted is named above the picture, beside Tissue / Side by side / UMAP, with Show on slide (frame it on the tissue) and Clear selection (or press Esc)." },
    panelOwn: {
      help: "Load results you already have for these bins, with no new export: a CSV or TSV (also .gz) with one row per bin, a barcode or cell id column, and any of a UMAP (two columns such as UMAP_1, UMAP_2), groupings (seurat_clusters, leiden, cell_type, sample …), numbers (scores, QC) and a colour per group (a column named like cell_type_color with #rrggbb). Press Load your own results or drop the file anywhere on the page. The file is read in this browser and never uploaded. Each column's use is guessed from its name and values; change it in the table (a colour list in the order of the groups, as scanpy keeps it, can be pasted under a grouping). The rows are matched to the bins by id, also with a -1 style ending or a sample prefix that differs; the share matched is shown, and unmatched rows are counted, never listed. The groupings then appear under Clustering in the Clusters section (colours on the tissue and the UMAP, names, volcano and plots), the numbers in View, Colour by a value, and the UMAP in the Tissue / UMAP switch (UMAP shown picks yours or SODA's). Matching needs the export's barcode list (prepare.py ids adds it). Save JSON can keep the matched table (its size is shown first); a link to the view keeps only a fingerprint of the file and asks for the same file. How to write the file from Seurat or scanpy: docs/ADVANCED.md." },
    panelBinQc: { more: ["bqGrey"], moreAt: "end", notes: ["bqNote"],
      help: "One histogram per measure (counts, genes detected, % mitochondrial). Move a slider or type a number to set a minimum or a maximum: bins outside are removed." },
    panelLasso: { more: ["lassoMode", "lassoSeed", "lassoRes"], notes: ["lassoNote"],
      help: "Pick Circle or Freehand and draw one or more shapes on the slide (Esc stops drawing). Run UMAP makes a new UMAP and clusters from only the bins inside the shapes. More options: Dots (one per bin, or one per shape with its bins added up), Seed (the same seed gives the same picture), Resolution (auto, or a number: higher gives more, smaller clusters; 0 = no clusters)." },
    panelShapes: {
      help: "Save GeoJSON writes the Draw a region shapes (the ticked ones, or all) in the image's full-resolution pixels, with the pixel size in µm: an outline as a polygon, a selection as the centres of its bins. Load GeoJSON adds polygons from such a file; one from QuPath lines up only if QuPath opened the same full-resolution image (same crop, flip and scale). Save barcodes writes one row per bin inside the shapes with its barcode or cell id, its position and its shape; it needs the export's barcode lists, which prepare.py export writes by default (not with --no-ids)." },
    panelVolcano: { more: ["vcTest", "vcMinPct", "vcExtra"],
      help: "Compares two groups of bins gene by gene: a drawn shape, a cluster, or the rest of the section. Genes to the right are higher in A, to the left higher in B; click a dot to draw that gene on the slide. × closes the plot; while nothing has changed the button reads Show volcano and brings it back without running again. More options: the test (Wilcoxon rank-sum or Welch's t), the smallest share of bins a gene must be counted in to be tested, and extra genes to test by name." },
    panelPlots: { more: ["plTop"], moreAt: { before: "plRun" },
      help: "Draws genes across groups (clusters, drawn shapes, regions or annotations) as a heatmap, a dot plot or violins. Pick the grouping, type gene names or let it take each group's top marker genes (how many: More options), then Draw plots. + shown genes adds the genes now on the picture. × closes the plots; while nothing has changed the button reads Show plots and brings them back without drawing again." },
    panelSave: { help: "Save JSON writes a file with the drawn shapes, the layer, the UMAP settings, the clustering on show with your cluster names, and the quality filter. Load JSON reads such a file and runs it again. Copy link to this view copies an address that opens this dataset as it looks now: position and zoom, layer, genes, stains, colours, clusters, the quality filter and the open tab. Drawn shapes too big for a link are left out, and it says so. Save JSON keeps the same view in a file." },
    panelImage: { more: ["snapBar", "snapName"], moreAt: "end",
      help: "Saves the map as it is on screen (Tissue, UMAP or Side by side) as a PNG file, with a legend of the colours on show, the scale bar and the dataset's name. 1x is the screen's own pixels; 2x and 4x draw the same view again with more pixels and finer image tiles, which takes a moment while the tiles load. Panels over the map (volcano, plots) are not in the picture: they have their own Save buttons. The camera beside Fit saves with these settings." },
    panelDataset: { help: "The size of the image, its pixel size, and how many bins and genes each layer holds." },
    panelHelp: {},
  };
  const ORDER = {
    view: ["panelChannels", "panelGenes", "panelMolecules", "panelOverlays", "panelRna"],
    analyze: ["panelUmap", "panelLasso", "panelShapes", "panelVolcano", "panelPlots"],
    share: ["panelImage", "panelSave", "panelDataset", "panelHelp"],
  };
  // the sections few people need: behind one "More options" at the end of their tab, in this order
  const MORE = {
    view: ["panelPaint"],
    analyze: ["panelOwn", "panelBinQc"],
  };

  // ---- remembered per browser: the tab, which sections are open, which "More options" ----------
  const KEY = "soda.layout.v1";
  const store = { tab: null, open: {}, more: {}, hidden: false };
  try {
    const j = JSON.parse(localStorage.getItem(KEY) || "null");
    if (j && typeof j === "object") {
      const tab = { clusters: "analyze", save: "share" }[j.tab] || j.tab;      // the four tabs before 0.10.0
      if (TABS.includes(tab)) store.tab = tab;
      if (j.open && typeof j.open === "object") Object.assign(store.open, j.open);
      if (j.more && typeof j.more === "object") Object.assign(store.more, j.more);
      store.hidden = !!j.hidden;
    }
  } catch (e) { /* private mode, storage blocked: the page works without it */ }
  const remember = () => { try { localStorage.setItem(KEY, JSON.stringify(store)); } catch (e) { /* as above */ } };

  // ---- small helpers ----------------------------------------------------------------------------
  const CHEV = '<svg class="chev" viewBox="0 0 12 12" aria-hidden="true"><path d="M4 2l4 4-4 4"/></svg>';
  const el = (tag, cls, id) => { const e = document.createElement(tag); if (cls) e.className = cls; if (id) e.id = id; return e; };
  const txt = (id) => { const e = $(id); return e ? e.textContent.replace(/\s+/g, " ").trim() : ""; };
  const up = (id) => { const e = $(id); return !!e && !e.hidden; };                 // there, and not hidden itself
  const hide = (e, v) => { if (e && e.hidden !== v) e.hidden = v; };            // only when it changes: every write wakes the observer
  const plural = (n, w) => `${n} ${w}${n === 1 ? "" : "s"}`;
  const phone = window.matchMedia ? window.matchMedia("(max-width: 700px)") : { matches: false };
  const sections = (tab) => [...tabBox[tab].querySelectorAll(":scope > section.panel, :scope > .tab-more > .tab-more-box > section.panel")];
  const tabOf = (sec) => { const b = sec.closest(".tabpanel"); return b ? TABS.find((t) => tabBox[t] === b) : null; };
  // the section to fold for an element: the outermost one
  const sectionOf = (node) => { let s = node && node.closest ? node.closest("section.panel") : null; while (s && s.parentElement && s.parentElement.closest("section.panel")) s = s.parentElement.closest("section.panel"); return s; };
  let words = null;                                 // viewer.js's wording function for the layer on show
  const say = (t) => (words ? words(t) : t);
  // Keep `dst` hidden exactly when any of `srcs` is. A block that used to sit inside a panel was hidden
  // with that panel; in a section of its own it has to follow the panel as well as itself, or it would show a
  // previous dataset's content (viewer.js only resets what the panel's hiding used to cover).
  function mirror(dst, ...srcs) {
    srcs = srcs.filter(Boolean);
    const f = () => hide(dst, srcs.some((x) => x.hidden));
    const mo = new MutationObserver(f);
    srcs.forEach((x) => mo.observe(x, { attributes: true, attributeFilter: ["hidden"] }));
    f();
  }

  // ---- build: the body, the wrapper sections, the top bar ---------------------------------------
  let current = "view", body = null;
  function build() {
    body = el("div", "side-body", "sideBody");
    side.appendChild(body);
    const brand = side.querySelector("header.brand"), mark = side.querySelector("footer.soda-mark");
    if (brand) body.appendChild(brand);
    body.appendChild(panelsBox);
    if (mark) panelsBox.appendChild(mark);          // under the sections: pinned to the bottom edge on a desktop, at the end of the sheet on a phone
    panelsBox.hidden = false;

    // the first panel held the picker, Fit and the stain buttons: it becomes "Dataset details"
    const dsRow = side.querySelector(".row.ds-row"), fit = $("btnFit"), vs = $("viewSwitch");
    const stainBtns = $("btnAllOn") ? $("btnAllOn").parentElement : null;
    if (dsRow) topBar.appendChild(dsRow);
    if (vs) {
      topBar.appendChild(vs);
      // Play is rarely used: behind a small button after Tissue / Side by side / UMAP
      const play = $("flyPlay");
      if (play && play.parentElement === vs) {
        const more = el("button", "fly-more", "flyMore"), box = el("span", "fly-more-box", "flyPlayBox");
        more.type = "button"; more.textContent = "…"; more.title = "Play the flight between the tissue and the UMAP";
        more.setAttribute("aria-label", "Play the flight"); more.setAttribute("aria-expanded", "false"); more.setAttribute("aria-controls", "flyPlayBox");
        box.hidden = true;
        play.before(more); more.after(box);
        box.appendChild(play);
        more.onclick = () => { box.hidden = !box.hidden; more.setAttribute("aria-expanded", String(!box.hidden)); if (!box.hidden && !play.hidden) play.focus(); };
        mirror(more, play);
      }
    }
    if (fit) topBar.appendChild(fit);
    topBar.hidden = false;
    const chan = $("panelChannels");
    if (stainBtns && chan && stainBtns.classList.contains("btnrow")) {
      stainBtns.classList.add("stain-btns"); chan.insertBefore(stainBtns, $("imgOpacity") ? $("imgOpacity").closest("label") : null);
      // an H&E has no stain buttons: no empty row above Image opacity
      const syncRow = () => hide(stainBtns, [...stainBtns.children].every((b) => b.hidden));
      new MutationObserver(syncRow).observe(stainBtns, { subtree: true, attributes: true, attributeFilter: ["hidden"] });
      syncRow();
    }
    const info = $("apiInfo"), dataset = $("panelDataset");
    if (dataset) {
      const ul = el("ul", "ds-layers", "dsLayers");
      dataset.appendChild(ul);
      if (info) dataset.appendChild(info);
    }

    // About and help: the logos and the two links join the help text; the title stays on top
    const help = $("panelHelp");
    if (help && brand) {
      const about = el("div", "brand-about");
      for (const id of ["brandLogos", "btnAbout", "lnkPaper"]) if ($(id)) about.appendChild($(id));
      const h = help.querySelector("h2");
      help.insertBefore(about, h ? h.nextSibling : help.firstChild);
    }

    // Save and load: Save JSON and Load JSON leave the Draw a region panel; "Copy link" is wired by viewlink.js
    // (Copy link to this view)
    const lasso = $("panelLasso");
    const save = newSection("panelSave", "share", "Save and load");
    const jsonRow = el("div", "btnrow", "saveJsonRow"), linkRow = el("div", "btnrow"), msg = el("div", "meta save-msg", "saveMsg");
    for (const id of ["lassoSaveJson", "lassoLoadBtn", "lassoLoad"]) if ($(id)) jsonRow.appendChild($(id));
    const copy = el("button", "", "btnCopyLink");
    copy.type = "button"; copy.textContent = "Copy link";
    linkRow.appendChild(copy);
    msg.setAttribute("role", "status");
    save.append(jsonRow, linkRow, msg);
    side.appendChild(save);
    if (lasso) mirror(jsonRow, lasso);

    // every section into its tab, in the tab's order; then its heading and body
    for (const tab of TABS) {
      const want = ORDER[tab].map($).filter(Boolean), extra = MORE[tab] || [];
      const rest = [...side.querySelectorAll(`section.panel[data-tab=${tab}]`)].filter((s) => !want.includes(s) && !extra.includes(s.id));
      for (const s of [...want, ...rest]) tabBox[tab].appendChild(s);
      if (extra.length) buildTabMore(tab);
      const none = el("div", "tab-empty"); none.textContent = "Nothing here for this dataset."; none.hidden = true;
      tabBox[tab].appendChild(none);
    }
    for (const tab of TABS) sections(tab).forEach(wrap);

    app.classList.add("tabbed");
    rail.hidden = false;
  }
  // a tab's "More options": one button at the end of the tab, the sections of MORE under it (folded at first)
  const tabMore = {};                               // tab -> { wrap, btn, box, sum, dot }
  function buildTabMore(tab) {
    const wrapEl = el("div", "tab-more", `tabMore-${tab}`), btn = el("button", "more-btn tab-more-btn", `tabMoreBtn-${tab}`),
          box = el("div", "tab-more-box", `tabMoreBox-${tab}`), sum = el("span", "tab-more-sum"), dot = el("i", "sec-dot");
    btn.type = "button"; btn.innerHTML = CHEV + "<span>More options</span>";
    btn.setAttribute("aria-controls", box.id);
    dot.hidden = true; dot.title = "Something changed here"; dot.setAttribute("role", "img"); dot.setAttribute("aria-label", "changed");
    btn.append(sum, dot);
    wrapEl.append(btn, box);
    tabBox[tab].appendChild(wrapEl);
    for (const s of MORE[tab].map($).filter(Boolean)) box.appendChild(s);
    tabMore[tab] = { wrap: wrapEl, btn, box, sum, dot };
    setTabMore(tab, !!store.more[`tab:${tab}`], false);
    btn.onclick = () => setTabMore(tab, box.hidden, true);
  }
  function setTabMore(tab, open, byUser) {
    const m = tabMore[tab]; if (!m) return;
    m.box.hidden = !open; m.btn.setAttribute("aria-expanded", String(open)); m.btn.classList.toggle("open", open);
    if (byUser) { store.more[`tab:${tab}`] = open; remember(); }
    if (open) { seen(); announceShown(); }
    schedule();
  }
  // the fold shows only when one of its sections does; folded, it names them
  function paintTabMore() {
    for (const tab of Object.keys(tabMore)) {
      const m = tabMore[tab], live = [...m.box.children].filter((s) => s.matches("section.panel") && !s.hidden);
      hide(m.wrap, !live.length);
      const t = m.box.hidden ? live.map((s) => { const n = s.querySelector(":scope > .sec-head .sec-name"); return n ? n.textContent.trim() : ""; }).filter(Boolean).join(", ") : "";
      if (m.sum.textContent !== t) m.sum.textContent = t;
    }
  }
  const behindMore = (sec) => !!sec && !!sec.closest(".tab-more-box[hidden]");
  function newSection(id, tab, title) {
    const s = el("section", "panel", id), h = el("h2");
    s.dataset.tab = tab; h.textContent = title;
    s.appendChild(h);
    return s;
  }

  // ---- a section: heading button (title, summary, dot), "?", body, "More options" ---------------
  function wrap(sec) {
    if (sec.dataset.lay) return;
    sec.dataset.lay = "1";
    const id = sec.id, cfg = SECTIONS[id] || {};
    let h2 = sec.querySelector(":scope > h2");
    if (!h2) { h2 = el("h2"); h2.textContent = id; }
    h2.classList.remove("lay-only");
    const head = el("div", "sec-head"), btn = el("button", "sec-toggle", `secBtn-${id}`), name = el("span", "sec-name"),
          sum = el("span", "sec-sum", `secSum-${id}`), dot = el("i", "sec-dot"), bodyEl = el("div", "sec-body", `secBody-${id}`);
    btn.type = "button"; btn.setAttribute("aria-controls", bodyEl.id); btn.setAttribute("aria-expanded", "true");
    btn.innerHTML = CHEV;
    while (h2.firstChild) name.appendChild(h2.firstChild);      // the title's own nodes (#geneTitle, #chanCount ...) keep their ids
    dot.hidden = true; dot.title = "Something changed here";
    dot.setAttribute("role", "img"); dot.setAttribute("aria-label", "changed");      // read out with the title while it shows
    btn.append(name, sum, dot);
    h2.appendChild(btn);
    head.appendChild(h2);
    while (sec.firstChild) bodyEl.appendChild(sec.firstChild);
    sec.append(head, bodyEl);
    btn.onclick = () => setOpen(sec, sec.classList.contains("closed"), true);

    if (cfg.help || (cfg.notes && cfg.notes.length)) {
      const q = el("button", "sec-q", `secQ-${id}`), help = el("div", "sec-help", `secHelp-${id}`);
      q.type = "button"; q.textContent = "?"; q.title = "What this section does";
      q.setAttribute("aria-label", "What this section does"); q.setAttribute("aria-expanded", "false"); q.setAttribute("aria-controls", help.id);
      help.hidden = true; help.setAttribute("role", "note");
      if (cfg.help) { const p = el("p"); p.dataset.text = cfg.help; p.textContent = say(cfg.help); help.appendChild(p); }
      for (const n of cfg.notes || []) if ($(n)) help.appendChild($(n));
      q.onclick = () => { help.hidden = !help.hidden; q.setAttribute("aria-expanded", String(!help.hidden)); };
      q.addEventListener("keydown", (e) => { if (e.key === "Escape" && !help.hidden) { help.hidden = true; q.setAttribute("aria-expanded", "false"); } });
      head.appendChild(q);
      sec.insertBefore(help, bodyEl);
    }
    if (cfg.more) buildMore(sec, cfg, bodyEl);
    setOpen(sec, isOpenWanted(id), false);
  }
  const rowOf = (id) => { const e = $(id); return e ? (e.closest("label") || e) : null; };
  function buildMore(sec, cfg, bodyEl) {
    const rows = cfg.more.map(rowOf).filter((r) => r && sec.contains(r));
    if (!rows.length) return;
    const btn = el("button", "more-btn", `secMoreBtn-${sec.id}`), box = el("div", "more-box", `secMore-${sec.id}`);
    btn.type = "button"; btn.innerHTML = CHEV + "<span>More options</span>";
    btn.setAttribute("aria-controls", box.id);
    const at = cfg.moreAt;
    if (at === "end") bodyEl.appendChild(btn);
    else if (at && at.before && $(at.before)) { const r = $(at.before).closest(".btnrow") || $(at.before); r.parentNode.insertBefore(btn, r); }
    else rows[0].parentNode.insertBefore(btn, rows[0]);
    btn.after(box);
    rows.forEach((r) => box.appendChild(r));
    const set = (open) => { box.hidden = !open; btn.setAttribute("aria-expanded", String(open)); btn.classList.toggle("open", open); };
    set(!!store.more[sec.id]);
    btn.onclick = () => { const open = box.hidden; set(open); store.more[sec.id] = open; remember(); if (open) announceShown(); };
    // nothing behind the button (roi.js hides Dots, Seed and Resolution on ROI layers): no button
    const sync = () => { const none = rows.every((r) => r.hidden); if (btn.hidden !== none) btn.hidden = none; box.classList.toggle("none", none); };
    const mo = new MutationObserver(sync);
    rows.forEach((r) => mo.observe(r, { attributes: true, attributeFilter: ["hidden"] }));
    sync();
  }

  // ---- open / closed -----------------------------------------------------------------------------
  // which section is open when nothing is remembered: the one in use for this kind of dataset
  function layerKind() {
    try {
      const S = window.__spatialViewer && window.__spatialViewer.state;
      if (!S || !S.m) return null;
      if (S.mode === "rna" || S.m.kind === "rnaseq") return "rna";
      const e = (S.m.expr || []).find((x) => x.key === S.exprKey);
      return e ? (e.modality === "protein" ? "protein" : e.kind || "bins") : "none";   // a lab's own cells of a protein image: protein
    } catch (e) { return null; }
  }
  function defaultOpen(id) {
    if (id === "panelGenes") return up("panelGenes") && layerKind() !== "protein";
    if (id === "panelChannels") return !(up("panelGenes") && layerKind() !== "protein");
    return ["panelRna", "panelUmap", "panelLasso", "panelSave"].includes(id);
  }
  const isOpenWanted = (id) => (id in store.open ? !!store.open[id] : defaultOpen(id));
  function setOpen(sec, open, byUser) {
    sec.classList.toggle("closed", !open);
    const btn = $(`secBtn-${sec.id}`);
    if (btn) btn.setAttribute("aria-expanded", String(open));
    if (byUser) { store.open[sec.id] = open; remember(); }
    if (open) { seen(); announceShown(); }
    schedule();
  }
  // When a dataset has opened: the remembered state, else the default. A tab's only section is always open
  // (a lone folded header helps nobody), and a tab whose sections are all folded by default, because the usual
  // open one is not there for this dataset, opens its first; sections the user folded themselves stay folded.
  function applyOpen() {
    for (const tab of TABS) {
      const secs = sections(tab);
      secs.forEach((s) => setOpen(s, isOpenWanted(s.id), false));
      const live = secs.filter((s) => !s.hidden);
      if (!live.length || !live.every((s) => s.classList.contains("closed"))) continue;
      if (live.length === 1 || !live.some((s) => store.open[s.id] === false)) setOpen(live.find((s) => defaultOpen(s.id)) || live[0], true, false);
    }
  }
  // histograms drawn while their box was folded or in another tab need drawing again (qc.js, viewer.js listen)
  let shownRaf = 0;
  function announceShown() {
    if (shownRaf) return;
    shownRaf = requestAnimationFrame(() => { shownRaf = 0; window.dispatchEvent(new CustomEvent("soda:shown")); });
  }

  // ---- tabs -----------------------------------------------------------------------------------
  let wanted = store.tab || "view";                 // the tab the user last chose
  const scrollAt = {};
  const available = (tab) => sections(tab).some((s) => !s.hidden);
  const sheet = () => app.dataset.sheet || "peek";
  const bodyHidden = () => (phone.matches ? sheet() === "peek" : app.classList.contains("side-hidden"));
  function select(tab, byUser) {
    if (!TABS.includes(tab)) return;
    if (tab !== current) scrollAt[current] = panelsBox.scrollTop;
    current = tab;
    for (const t of TABS) {
      const on = t === tab;
      tabBtn[t].setAttribute("aria-selected", String(on));
      tabBtn[t].tabIndex = on ? 0 : -1;
      tabBox[t].hidden = !on;
    }
    panelsBox.scrollTop = scrollAt[tab] || 0;
    if (byUser) { wanted = tab; store.tab = tab; remember(); }
    syncSheetAria();
    seen(); paintFlags(); announceShown(); schedule();
  }
  function setSideHidden(hide, byUser) {
    app.classList.toggle("side-hidden", hide);
    const b = $("railHide");
    if (b) { b.setAttribute("aria-expanded", String(!hide)); b.title = hide ? "Show the controls" : "Hide the controls (more room for the picture)"; const s = b.querySelector(".sr-only"); if (s) s.textContent = hide ? "Show the controls" : "Hide the controls"; }
    if (byUser) { store.hidden = hide; remember(); }
    if (!hide) { seen(); announceShown(); }
    paintFlags();
  }
  // phones: the sheet's three heights are viewer.js's (data-sheet on .app); a tab opens it, the tab on show closes it
  function setSheet(s) {
    app.dataset.sheet = s; side.style.height = "";
    const h = $("sheetHandle"); if (h) h.setAttribute("aria-expanded", String(s !== "peek"));
  }
  function syncSheetAria() {
    for (const t of TABS) {
      if (phone.matches && t === current) tabBtn[t].setAttribute("aria-expanded", String(sheet() !== "peek"));
      else tabBtn[t].removeAttribute("aria-expanded");
    }
  }
  function tabEvents() {
    for (const t of TABS) {
      const d = tabBtn[t].querySelector(".tab-dot");
      if (d) { d.setAttribute("role", "img"); d.setAttribute("aria-label", "something changed"); }
      tabBtn[t].addEventListener("click", () => {
        if (phone.matches) {
          if (sheet() === "peek") { select(t, true); setSheet("half"); }
          else if (t === current) setSheet("peek");
          else select(t, true);
        } else {
          if (app.classList.contains("side-hidden")) setSideHidden(false, true);
          select(t, true);
        }
      });
    }
    tabList.addEventListener("keydown", (e) => {
      const live = TABS.filter((t) => !tabBtn[t].hidden), i = live.indexOf(e.target.dataset ? e.target.dataset.tab : "");
      if (i < 0) return;
      let j = -1;
      if (e.key === "ArrowDown" || e.key === "ArrowRight") j = (i + 1) % live.length;
      else if (e.key === "ArrowUp" || e.key === "ArrowLeft") j = (i - 1 + live.length) % live.length;
      else if (e.key === "Home") j = 0;
      else if (e.key === "End") j = live.length - 1;
      if (j < 0) return;
      e.preventDefault();
      select(live[j], true);
      tabBtn[live[j]].focus();
    });
    const hideBtn = $("railHide");
    if (hideBtn) hideBtn.onclick = () => setSideHidden(!app.classList.contains("side-hidden"), true);
    new MutationObserver(() => { syncSheetAria(); if (sheet() !== "peek") { seen(); announceShown(); } paintFlags(); })
      .observe(app, { attributes: true, attributeFilter: ["data-sheet"] });
    const onMq = () => { syncSheetAria(); paintFlags(); announceShown(); };
    if (phone.addEventListener) phone.addEventListener("change", onMq); else if (phone.addListener) phone.addListener(onMq);
  }
  // Which tabs exist for this dataset, and which one shows. `opened`: a dataset has just finished opening, so
  // the tab the user last chose comes back if it is there; otherwise the tab on show only changes when it has emptied.
  function settle(opened) {
    for (const t of TABS) { const ok = available(t); if (tabBtn[t].hidden === ok) tabBtn[t].hidden = !ok; }
    const first = TABS.find(available) || "view";
    const next = opened && available(wanted) ? wanted : available(current) ? current : available(wanted) ? wanted : first;
    if (next !== current) select(next, false);
    for (const t of TABS) hide(tabBox[t].querySelector(":scope > .tab-empty"), available(t));
  }

  // ---- signals: a dot on a tab and on a section header when something changed out of sight -------
  const flags = new Map();                          // section id -> "info" | "warn"
  const onShow = (sec) => !!sec && !sec.hidden && tabOf(sec) === current && !sec.classList.contains("closed") && !behindMore(sec) && !bodyHidden();
  function flag(id, kind) {
    const sec = $(id);
    if (!sec || sec.hidden || opening > 0 || onShow(sec)) return;
    if (flags.get(id) !== "warn") flags.set(id, kind || "info");
    paintFlags();
  }
  function seen() {
    let any = false;
    for (const id of [...flags.keys()]) if (onShow($(id)) || !$(id) || $(id).hidden) { flags.delete(id); any = true; }
    if (any) paintFlags();
  }
  // a state that lasts (not a one-off change): marked for as long as it holds and its section is out of sight
  const standing = (id) => (id === "panelLasso" && up("bqStale") && up("panelBinQc") ? "warn" : null);   // the result on show was made before the quality filter changed
  const toolOn = (tab) => (tab === "analyze" ? ["lassoCircle", "lassoFree", "selBox", "selFree"] : []).some((id) => $(id) && $(id).classList.contains("on"));
  function paintFlags() {
    for (const tab of TABS) {
      let kind = null;
      let behind = null;                            // the strongest mark among the sections behind a folded More options
      for (const s of sections(tab)) {
        const k = flags.get(s.id) || (!s.hidden && !onShow(s) ? standing(s.id) : null), d = s.querySelector(":scope > .sec-head .sec-dot");
        if (d) { hide(d, !k); d.classList.toggle("warn", k === "warn"); }
        if (k && kind !== "warn") kind = k;
        if (k && behindMore(s) && behind !== "warn") behind = k;
      }
      if (tabMore[tab]) { const d = tabMore[tab].dot; hide(d, !behind); d.classList.toggle("warn", behind === "warn"); }
      const away = tab !== current || bodyHidden();
      const d = tabBtn[tab].querySelector(".tab-dot"), show = away && (!!kind || toolOn(tab));
      if (d) {
        hide(d, !show); d.classList.toggle("warn", kind === "warn");
        const tip = kind ? "Something changed in this tab" : "A drawing tool is on in this tab";
        if (d.title !== tip) d.title = tip;
      }
      tabBtn[tab].classList.toggle("has-dot", show);
    }
  }
  function watchSignals() {
    // a run ends (its Cancel button goes away): the result or the error is in that section
    watchRun("lassoCancel", "panelLasso");
    // the quality thresholds changed since the last run: the warning is beside the sliders (Analyze tab, More options); the
    // button that applies them is Run UMAP, so Draw a region says it too and carries a dot while it holds (standing)
    const stale = $("bqStale"), run = $("lassoRun");
    if (stale && run) {
      const note = el("div", "meta bq-stale", "lassoStale");
      note.textContent = "The quality filter changed since this run: press Run UMAP to apply it.";
      (run.closest(".btnrow") || run).after(note);
      mirror(note, stale, $("panelBinQc"));
    }
    // a gene added from somewhere else (a volcano dot, a plot): it is now in the View tab's list
    for (const [list, secId] of [["geneList", "panelGenes"], ["rnaGeneList", "panelRna"]]) {
      const box = $(list); if (!box) continue;
      let n = box.children.length;
      new MutationObserver(() => { const now = box.children.length; if (now > n) flag(secId); n = now; }).observe(box, { childList: true });
    }
    // the user's own action sends them to another tab: show it
    const load = $("lassoLoad");
    if (load) load.addEventListener("change", () => show("panelLasso"));      // a file was chosen: the run and its messages are in Draw a region
  }
  function watchRun(cancelId, secId) {
    const c = $(cancelId); if (!c) return;
    let was = c.hidden;
    new MutationObserver(() => { if (c.hidden && !was) flag(secId); was = c.hidden; schedule(); }).observe(c, { attributes: true, attributeFilter: ["hidden"] });
  }
  function watchMsg(msgId, secId) {
    const m = $(msgId); if (!m) return;
    new MutationObserver(() => { if (m.textContent.trim() && m.style.color) flag(secId, "warn"); })
      .observe(m, { childList: true, characterData: true, subtree: true, attributes: true, attributeFilter: ["style"] });
  }
  // after a dataset has opened: what is already on that a first look at the View tab would miss
  function flagOpenState() {
    if (txt("bqCount")) flag("panelBinQc", "warn");
  }

  // bring a section in front of the user: its tab, unfolded, in view (and on a phone the sheet up)
  function show(id) {
    const sec = sectionOf($(id));
    if (!sec || sec.hidden) return false;
    const tab = tabOf(sec);
    if (!tab) return false;
    if (tabBtn[tab].hidden) tabBtn[tab].hidden = false;
    if (phone.matches) { if (sheet() === "peek") setSheet("half"); }
    else if (app.classList.contains("side-hidden")) setSideHidden(false, false);
    select(tab, true);
    if (behindMore(sec)) setTabMore(tab, true, false);
    setOpen(sec, true, false);
    requestAnimationFrame(() => { const t = $(id); if (t && t.scrollIntoView) t.scrollIntoView({ block: "nearest" }); });
    return true;
  }

  // ---- summaries on the folded headers ---------------------------------------------------------------
  const names = (listId, none) => {
    const rows = $(listId) ? [...$(listId).querySelectorAll(".chan .name")] : [];
    const n = rows.map((r) => (r.firstChild && r.firstChild.nodeType === 3 ? r.firstChild.textContent : r.textContent).trim()).filter(Boolean);
    return !n.length ? none : n.length <= 2 ? n.join(", ") : `${n[0]}, ${n[1]} +${n.length - 2}`;
  };
  const SUMMARY = {
    panelChannels() {
      const list = $("chanList");
      if (list && !list.hidden) return list.children.length ? `${plural(list.querySelectorAll(".chan:not(.off)").length, "stain")} on` : "";
      const o = $("imgOpacity");                    // Image opacity at 0: the image is off
      return o && +o.value <= 0 ? "image off" : "image on";
    },
    panelGenes: () => names("geneList", "none picked"),
    panelPaint: () => (window.SodaPaint ? window.SodaPaint.summary() : ""),
    panelMolecules: () => names("molList", "none picked"),
    panelOverlays() { const n = $("overlayList") ? $("overlayList").querySelectorAll("input[type=checkbox]:checked").length : 0; return n ? `${n} on` : "off"; },
    panelRna: () => names("rnaGeneList", "clusters, no gene picked"),
    panelUmap() {
      const by = $("colorBy"); if (!by) return "";
      const sel = up("selChip") ? txt("selChipTxt") : "";
      let t;
      if (by.value === "cluster") { const o = $("clusterRes") && $("clusterRes").selectedOptions[0]; t = o ? o.textContent.split(" · ")[0] : "clusters"; }
      else if (by.value === "celltype") t = "suggested cell types";
      else t = `by ${(by.options[0] ? by.options[0].textContent : "genes").split(" ")[0].toLowerCase()}`;
      return sel ? `${t} · ${sel}` : t;
    },
    panelBinQc: () => txt("bqCount").replace(/[()]/g, "") || "off",
    panelLasso() {
      const n = $("lassoShapes") ? $("lassoShapes").querySelectorAll(".lshape").length : 0;
      return (n ? plural(n, "shape") : "no shapes") + (up("lassoCancel") ? " · running…" : toolOn("analyze") ? " · drawing…" : "");
    },
    panelShapes() { const n = $("shPick") ? $("shPick").querySelectorAll("input:checked").length : 0; return n ? `${n} ticked` : ""; },
    panelVolcano() {
      if (up("vcCancel")) return "running…";
      const t = txt("vcTitle"), has = !!(window.SodaVolcano && window.SodaVolcano.result) || up("volcanoPanel");
      return has && t && t !== "Volcano" ? t : "not run";
    },
    panelPlots() {
      if (up("plCancel")) return "running…";
      const t = txt("plTitle"), has = !!(window.SodaPlots && window.SodaPlots.result) || up("plotsPanel");
      return has && t && t !== "Plots" ? t : "not drawn";
    },
    panelDataset() { const p = txt("dsInfo").split(" · "); return p.length > 3 ? `${p[0]} · ${p[p.length - 1]}` : p.slice(0, 2).join(" · "); },   // size and pixel size
  };
  function paintSummaries() {
    for (const tab of TABS) for (const s of sections(tab)) {
      const out = $(`secSum-${s.id}`), f = SUMMARY[s.id];
      if (!out || !f) continue;
      let t = ""; try { t = f() || ""; } catch (e) { t = ""; }
      if (out.textContent !== t) out.textContent = t;
      out.classList.toggle("warn", s.id === "panelBinQc" && t !== "off");
    }
  }
  // Dataset details: the layers and their counts (the Layer menu's own lines), the stains, the molecules
  function paintDetails() {
    const ul = $("dsLayers"); if (!ul) return;
    const lines = [];
    const sel = $("exprSelect");
    if (sel && up("panelGenes")) for (const o of sel.options) lines.push(o.textContent);
    const cl = $("chanList");
    if (cl && up("panelChannels") && !cl.hidden && cl.children.length) lines.push(plural(cl.children.length, "stain"));
    if (up("panelMolecules") && txt("molMeta")) lines.push(txt("molMeta"));
    const key = lines.join("\n");
    if (ul.dataset.key === key) return;
    ul.dataset.key = key; ul.textContent = "";
    for (const t of lines) { const li = el("li"); li.textContent = t; ul.appendChild(li); }
  }

  // ---- keeping up with the page -----------------------------------------------------------------
  let timer = 0, opening = 0, openTok = 0, firstOpen = false;
  function schedule() { if (!timer) timer = setTimeout(refresh, 60); }
  function refresh() {
    timer = 0;
    paintSummaries(); paintDetails(); paintTabMore();
    if (opening === 0 && firstOpen) settle(false);
    const stains = $("chanList"), chan = $("panelChannels");
    if (stains && chan) chan.classList.toggle("no-stains", stains.hidden);
    seen(); paintFlags();
  }
  function watch() {
    const mo = new MutationObserver(schedule);
    for (const root of [side, topBar]) mo.observe(root, { subtree: true, childList: true, characterData: true, attributes: true, attributeFilter: ["hidden", "class", "disabled"] });
    for (const ev of ["change", "input", "click"]) document.addEventListener(ev, schedule, true);
  }

  // ---- what the other scripts call ---------------------------------------------------------------
  window.SodaLayout = {
    // viewer.js: a dataset starts opening (true; returns that open's number) / has finished (false, its number).
    // An open overtaken by a later one never finishes, so only the latest number's end counts.
    opening(on, tok) {
      if (on) { opening = 1; return ++openTok; }
      if (tok !== openTok) return;
      opening = 0;
      firstOpen = true;
      if (layerKind() === "rna") this.wording(null);          // the labels of a single-cell dataset are not reworded
      applyOpen(); settle(true); flags.clear(); flagOpenState(); refresh();
    },
    // viewer.js (syncWording): the function that turns "bins / genes" into this layer's words, or null
    wording(fn) {
      words = typeof fn === "function" ? fn : null;
      for (const p of side.querySelectorAll("[data-text]")) { const t = say(p.dataset.text); if (p.textContent !== t) p.textContent = t; }
      schedule();
    },
    // volcano.js, plots.js, snapshot.js, selections.js: a block or a section built after this file ran
    adopt(node) {
      if (!node) return;
      if (node.id === "volcanoBlock") {
        let sec = $("panelVolcano");
        if (!sec) {
          sec = newSection("panelVolcano", "analyze", "Volcano plot");
          const after = $("panelLasso");
          if (after && after.parentElement === tabBox.analyze) after.after(sec); else tabBox.analyze.insertBefore(sec, tabBox.analyze.firstChild);
        }
        sec.appendChild(node);
        mirror(sec, node, $("panelLasso"));
        wrap(sec);
        watchRun("vcCancel", "panelVolcano"); watchMsg("vcMsg", "panelVolcano");
      } else if (node.matches && node.matches("section.panel")) {
        const tab = TABS.includes(node.dataset.tab) ? node.dataset.tab : "view";
        const more = tabMore[tab] && MORE[tab].includes(node.id);
        const box = more ? tabMore[tab].box : tabBox[tab], order = more ? MORE[tab] : ORDER[tab], at = order.indexOf(node.id);
        const next = at < 0 ? null : order.slice(at + 1).map($).find((s) => s && s.parentElement === box);
        const end = more ? null : tabMore[tab] ? tabMore[tab].wrap : tabBox[tab].querySelector(":scope > .tab-empty");
        box.insertBefore(node, next || end);
        wrap(node);
        if (node.id === "panelPlots") { watchRun("plCancel", "panelPlots"); watchMsg("plMsg", "panelPlots"); }
      }
      schedule();
    },
    show, select: (t) => select(t, true), flag,
    state: () => ({ tab: current, wanted, tabs: TABS.filter((t) => !tabBtn[t].hidden), open: TABS.flatMap((t) => sections(t)).filter((s) => !s.hidden && !s.classList.contains("closed")).map((s) => s.id), flags: [...flags.keys()] }),
  };

  build();
  tabEvents();
  watchSignals();
  watch();
  setSideHidden(!!store.hidden, false);
  select(available(wanted) ? wanted : "view", false);
  refresh();
  // no dataset ever opens (deck.gl did not load, no datasets listed): settle on what is there
  setTimeout(() => { if (!firstOpen && opening === 0) { firstOpen = true; settle(false); } }, 4000);
})();
