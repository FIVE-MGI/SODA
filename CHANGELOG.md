# Changelog

What is in each version of SODA, in plain words. The details are in the [docs](README.md#documentation).

## 0.9.0 (17 November 2026)

The first version for release, frozen for the bioRxiv preprint. It is what the paper describes.

### Input: what SODA reads

- Multiplexed immunofluorescence: OME-TIFF, QPTIFF, other multichannel TIFFs, and folders with
  one TIFF per marker (Akoya / PhenoCycler).
- H&E and other brightfield images.
- Visium and Visium HD (Space Ranger output). Visium HD builds the 8 µm and 16 µm layers by
  default; `--bins all` adds the 2 µm layer.
- Xenium: stains, cells, cell and nucleus outlines, and every transcript. A zipped Xenium
  output is read as it is, and an output without stain images opens on a transcript-density
  background.
- CosMx SMI and MERSCOPE: stains, cells, outlines and transcripts.
- Stereo-seq (STOmics): bin layers from a GEF or GEM file, and the cells of the cell bin file.
- GeoMx DSP: scans, ROI outlines and ROI groups, Q3-normalised by default.
- Single-cell RNA-seq from `.h5ad` or 10x output.
- Standard formats: spatial AnnData (`.h5ad` or AnnData-Zarr), OME-Zarr images (OME-NGFF 0.4
  and 0.5) and SpatialData stores, in zarr v2 or v3, without installing anndata or
  spatialdata.
- SODA looks at a folder and says what it holds (`prepare.py detect` / `auto`). When it is not
  sure it says so, and never guesses.

### Preparing and sharing

- One notebook (`notebooks/SODA_pipeline.ipynb`) from instrument folders to one link: two
  settings, a check of the inputs, a memory warning, a local preview, the upload command, and
  a test that the storage lets the page read it.
- `prepare.py batch`: a whole folder of datasets in one command, with a log of any that fail.
- One-command install (`install.sh`, `install.ps1`) from a lock file of tested package versions.
- Share a `?data=` link to your own storage folder, or copy SODA on GitHub as your own site.
  `datasets.json` sets names, groups, tags and which stains show.
- Gene data packed into a few large files per layer, read one gene at a time with HTTP Range
  requests (far fewer files to upload). `--gene-files loose` for hosts without Range.
- A password-protected site on Azure, with a small token function.
- Offline use: the two browser libraries from your own site (`prepare.py vendor`), and the
  page keeps itself after one visit.
- A desktop app for Windows and Mac (no Python needed), with a Windows installer.

### In the viewer

- Four tabs (View, Clusters, Analyze, Save), on desktop, tablet and phone.
- Any number of fluorescence channels, each with its own colour, range, gamma and colour map.
- Up to eight genes at once on the tissue, as points or a density map; two-gene co-expression
  (magenta / green by default, red / green as an option); gene-set scores.
- Tissue | Side by side | UMAP: a whole-section UMAP with Leiden clusters at three
  resolutions, linked selection between the tissue and the UMAP.
- Suggested cluster names from editable marker lists; cell sets (merge clusters into named
  sets).
- Your own UMAP and clusters (Seurat, scanpy), at export or loaded into the page from a CSV /
  TSV with no new export.
- Draw a region: a UMAP and clusters of just the bins inside the shapes you draw, computed in
  the browser, with a fast mode for regions of more than 20,000 bins. Shapes can be
  combined, saved and loaded as GeoJSON, and the run's clusters merged into named sets.
- Per-bin quality filter; slide QC and "% positive" for protein images.
- Volcano plot, heatmap, dot plot and violin plots of genes across groups.
- Cell and nucleus outlines from the instrument's segmentation or a lab's own file.
- Save image (PNG at 1x, 2x, 4x with legend and scale bar); a link that restores the whole
  view; Compare, to open a second dataset in a linked window.

### Privacy

- Sample IDs are replaced through a private key, and `deid.py check` gates every upload.
- Barcode and cell id lists are now in an export by default, after the sample ID check;
  `--no-ids` leaves them out.

### Testing

- `tests/smoke_test.py` runs every synthetic demo end to end and checks every output file.
  It runs on Linux on every push.
- Browser tests in Chromium, WebKit (Safari's engine), Firefox and Edge, including phone-sized
  windows, with Playwright.
- Performance scripts (`tools/perf/`) and measured times in `docs/PERFORMANCE.md`.

### Deprecated

- Encrypted link files (`links.py`, `linkPath` in `config.js`) still work but will be removed
  in a later version. Use the password site or a link with a read pass.

### Not in this version

- A SpatialData table whose region is a labels element (a segmentation mask) is refused with
  the reason: save the cells as shapes instead.
- Cell segmentation of protein images: they are measured on a grid of 16 µm squares.
- Opening `.h5ad`, OME-Zarr or SpatialData files in the browser without converting them first.
