// Copy to config.js and edit. config.js is gitignored: it names your datasets,
// your institution and your storage, none of which belong in this repository.
//
// The viewer reads nothing else. Everything specific to a study lives here.
window.SpatialViewerConfig = {

  // ---- datasets -----------------------------------------------------------------
  // One entry per folder in the storage container. The key is the folder name.
  //   group         which heading it appears under (see `groups` below)
  //   label         what the picker shows; the tissue and its condition, nothing else
  //   hideChannels  exact channel names to leave out of this dataset, e.g. stain
  //                 tests that are not panel markers. Matched case-insensitively.
  //   tags          optional short labels, e.g. { tissue: "nerve", condition: "injured",
  //                 slide: "S12" } (any names). The picker then gets one drop-down filter
  //                 per name, its search also finds the values, and a name that has the
  //                 same value everywhere gets no filter. Text only, a name up to 40
  //                 characters, a value up to 80. No tags on any dataset: no filters.
  datasets: {
    // example_visium:    { group: "transcriptomics", label: "Normal tissue",
    //                      tags: { tissue: "nerve", condition: "normal", slide: "S1" } },
    // example_proteomics: { group: "proteomics", label: "Condition A",
    //                       hideChannels: ["CD11c FFPE"] },
    // example_rnaseq:    { group: "scrnaseq", label: "All cells" },
  },

  // Headings for the picker, in the order they should appear. A dataset whose
  // `group` is missing from here falls under "Other"; with no groups at all the
  // picker is a flat list.
  groups: {
    transcriptomics: "Spatial transcriptomics",
    proteomics: "Spatial proteomics",
    scrnaseq: "Single-cell RNA-seq",
  },

  defaultDataset: "",

  // ---- whose instance this is ----------------------------------------------------
  // The viewer has no institution of its own. Logo paths are relative to this folder.
  brand: {
    title: "SODA",
    logos: [
      // { src: "logos/institution.png", alt: "Institution", className: "logo-wide" },
    ],
    byline: "",
  },

  // Shown as "Link to paper". Empty renders it greyed and inert.
  paperUrl: "",

  // ---- your own copy of SODA -------------------------------------------------------
  // The link to your export folder (the one holding datasets.json, written by catalog.py).
  // Set this and leave `datasets` above empty: the page then lists every dataset in that
  // folder, with the titles, groups and settings in its datasets.json.
  // data: "https://yourstorage.example.org/soda-export/",

  // Whether this page also opens other people's folders from a link
  // (index.html?data=<their folder>). On by default, off when tokenEndpoint is set.
  // allowDataLinks: true,

  // ---- data access ----------------------------------------------------------------
  // The managed function that checks the shared password and returns a short-lived,
  // read-only SAS for the private container. Empty disables the password entirely,
  // which only makes sense when the data is public or served locally.
  tokenEndpoint: "/api/token",

  // A live tile server (server.py) instead of, or as well as, static exports.
  // Empty means none; ?api=http://host at runtime also works.
  apiBase: "",
  useLocalServer: false,

  // ---- UMAP cluster labels --------------------------------------------------------
  // Hand-placed label positions in UMAP coordinates, keyed by dataset then cluster.
  // Anything listed wins over the automatic placement; anything missing is placed by
  // the solver, so a partial list is fine. Produced by tools/label_editor.html.
  labelPositions: {},

  // Marker-gene lists the circle -> UMAP tool uses to suggest cluster names (relative to
  // this folder). The export takes its own copy via prepare.py export --markers.
  markerSet: "markers/nerve.json",

  // Clustering columns to color the UMAP by, first match wins. The picker was removed
  // because each dataset here carries exactly one; add names as needed.
  rnaColorBy: ["leiden"],

  // ---- display --------------------------------------------------------------------
  // Most genes at once, so one dataset cannot exhaust the tab.
  rnaseqMaxGenes: 12,   // UMAP panels, one per gene
  maxGenes: 8,          // spatial transcript layers

  // Gene packs (exports since 16 Nov 2026) are read one gene at a time with HTTP Range
  // requests. A host that ignores Range sends the whole pack instead: up to this many MB the
  // page keeps it and cuts the genes out of it (slower), above it the gene is refused with a
  // message. Export with --gene-files loose for such a host.
  // packWholeMaxMB: 160,

  // Smallest bin layer offered. Below roughly 8 um a bin holds a single transcript
  // nearly everywhere, which makes a density map of it a presence map.
  exprMinBinUm: 8,

  // Channels never shown, whatever dataset they are in: unassigned detector slots.
  // A JavaScript regex source string; omit to use the built-in default.
  // hideChannels: "^\\s*(open|blank|empty|unused)([\\s_-]|$)",

  tileFormat: "png",    // "jpg" is ~4x smaller per tile (server mode only)
  tileCacheSize: 400,   // decoded tiles held per page, all channels together

  // Window-slider axis. Fluorescence piles up near zero, so a linear 0-255 slider
  // crams every useful window into the first few percent of travel. Positions map
  // through v = max * (p ** gamma); 1 is linear.
  sliderGamma: 3.2,

  // Histogram shape: count -> bar height, then normalised to this percentile rather
  // than the maximum so one background spike cannot flatten the rest.
  histGamma: 0.35,
  histClip: 0.985,

  // Colors cycled through as genes are added.
  geneColors: ["#ff3355", "#33ddff", "#ffdd33", "#66ff66", "#ff88ff", "#ff9933", "#99aaff", "#ffffff"],

  // Single molecules (Xenium): up to this many genes at once, their colors, and how many dots
  // are drawn when zoomed out (at most moleculeBudget, and about moleculesPerPixel per screen
  // pixel; the rest appear as you zoom in).
  // maxMoleculeGenes: 12,
  // moleculeColors: ["#ff3355", "#33ddff", "#ffdd33", "#66ff66", "#ff9933", "#cc88ff"],
  // moleculeBudget: 6e6,
  // moleculesPerPixel: 0.25,

  // Cell outlines (segmented cells with boundaries): drawn from the zoom at which the median cell is
  // outlineFromPx screen pixels across, alone from outlineFullPx (dots below, both in between), and
  // only while no more than outlineBudget outlines are in view (20,000 on a phone or tablet).
  // outlineFromPx: 5,
  // outlineFullPx: 9,
  // outlineBudget: 120000,

  // Floors for categorical colors, as relative luminance, against a dark background.
  // tab20 was drawn for white paper and several of its entries vanish otherwise.
  catMinLuminance: 0.21,    // the points
  labelMinLuminance: 0.5,   // cluster label text, which sits over its own points
};
