# Changelog

What is in each version of SODA, in plain words. The details are in the [docs](README.md#documentation).

## 0.10.0 (22 November 2026): simpler

Fewer ways to do the same thing, so there is less to learn and less that can break.

### Changed

- The notebook, the desktop app and `prepare.py batch` now build their commands with the same
  code, so a setting works the same way in all three.
- Exports make JPEG image tiles by default everywhere (it was PNG for `prepare.py export` on its
  own). `--fmt png` still gives lossless tiles.
- Two places that skipped data without a word now print a warning: an `.h5ad` column of a type
  SODA cannot read, and an unreadable Xenium `gene_category`.

- The page has three tabs instead of four: **View**, **Analyze** (which now also holds Clusters
  and selection) and **Share** (Save image, Save and load, Copy link, dataset details, help).
  The tools few people need are under **More options** at the end of their tab: Colour by a
  value and Measure (View); Your own results, Cell sets and the Quality filter (Analyze). Play
  (the flight between the tissue and the UMAP) is behind the **…** beside Tissue / Side by
  side / UMAP. A link or browser that remembered the Clusters or Save tab opens Analyze or
  Share.
- The docs start with the SODA app: the README is one page (see the examples, the app in four
  steps, the notebook for more control), and [GETTING_STARTED.md](docs/GETTING_STARTED.md)
  has the app as Part A before the notebook. What the README held beyond that (features,
  installing from the code, the command line, requirements, browsers and phones, the
  repository) is in [docs/REFERENCE.md](docs/REFERENCE.md).

### Removed

- Encrypted link files (`links.py`, `linkPath` in `config.js`), deprecated in 0.9.0. Use the
  password site or a link with a read pass.
- Reading a zipped dataset folder as it is. Unzip it first. (Xenium's own `transcripts.zarr.zip`
  and `cells.zarr.zip` inside a folder are still read.)
- `--seam-correct`. Tile seams are still found and listed in the slide QC.
- `--gene-files loose / both` and `--shard-mb`. Exports always use gene packs, which need a host
  that answers HTTP Range requests (blob storage, S3, GitHub Pages, nginx). The page still opens
  older exports with one file per gene.
- The older inputs layout with a folder per kind (`inputs/visium/<id>/`). Put each dataset's
  folder straight in the inputs folder; a folder holding several datasets is still split into
  one dataset each. The notebook now shows the same kind names as `prepare.py` (`tiff`, `akoya`,
  `rnaseq`); the old names still work in `OPTIONS`.
- The live tile server (`server.py`, `run.ps1`), which showed the prepare cache on a lab network
  without an export. With it go `?api=` and `?fmt=` in the address, `apiBase`, `useLocalServer`
  and `tileFormat` in `config.js`, and FastAPI, Uvicorn and the 17 packages only they needed. To
  look at a dataset on your own computer, export it and use the SODA app's **Open viewer** or
  [HOSTING.md A](docs/HOSTING.md#a-on-your-own-computer).

Removed from the page (old links and Save JSON files still open; the parts they held for these
tools are ignored):

- Counting tools for protein images: % positive for a marker, Count clusters, Combine markers,
  saved cut-offs and a dataset's `gates.json`.
- Compare side by side (a second, linked window).
- Combine / Overlap / Subtract in the shapes section, now called Export shapes (Save and Load
  GeoJSON, Save barcodes). A combined shape in an older link or file still opens.
- Merge clusters of the run (Analyze tab). Cell sets, for the whole section's clusters, stays.
- Compare with slide and its table (protein Draw a region). Run UMAP still sets the clusters'
  markers against the whole slide.
- The Slide QC box. The export still writes `qc.json` with the same checks; the page still uses
  `qc_flags.bin` for the "Artefact?" cluster flag.
- Offline use ("Keep this dataset" and the page keeping itself). For readers who used it,
  `sw.js` now only deletes what it kept and removes itself on their next visit.

## 0.9.0 (19 November 2026)

The first version for release, frozen for the bioRxiv preprint. It is what the paper describes.
New to SODA? Start with [docs/GETTING_STARTED.md](docs/GETTING_STARTED.md).

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
- A lab's own cell segmentation of a protein image (`--cells`: a label mask from Cellpose, Mesmer,
  StarDist or QuPath, or QuPath GeoJSON) becomes a cells layer with outlines, marker means per
  cell, UMAP and clusters, beside the grid of squares.
- Single-cell RNA-seq from `.h5ad` or 10x output. A file without a UMAP is shown with its t-SNE
  (or another 2-D layout, else its PCA), and SODA says which.
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
- One-command install (`install.sh`, `install.ps1`) from a lock file of tested package versions,
  into its own environment, `soda-env`, which VS Code and JupyterLab list as the kernel **SODA**.
- Share a `?data=` link to your own storage folder, or copy SODA on GitHub as your own site.
  `datasets.json` sets names, groups, tags and which stains show.
- `site.js` names the data folder a copy of SODA opens when the link names none (by default the
  four public example datasets). A link's `?data=` or `config.js` wins.
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

- No demo data: the synthetic demo datasets and their generators are removed; the examples are
  four public datasets at https://fivestore.blob.core.windows.net/soda/examples/.
- `tests/smoke_test.py` checks the tools and the public example exports. It runs on Linux on
  every push.
- Browser tests in Chromium, WebKit (Safari's engine), Firefox and Edge, including phone-sized
  windows, with Playwright.
- Performance scripts (`tools/perf/`) and measured times in `docs/PERFORMANCE.md`.

### Fixed before release

- Float images (32-bit) are scaled from 0 to their 99.99th percentile into 16 bits. Before,
  they were cast without scaling and looked posterized.
- Save image: a save right after a 2x save could come out twice too large. The picture's size
  now comes from the map's resting pixel ratio.
- Offline use: the service worker no longer re-creates an app cache that a publish deleted,
  and keeps the page only while its own worker is in charge.
- Whole-section UMAP: umap-learn 0.6 is asked for its compatibility layout, the same embedding
  as before, so every export no longer prints a FutureWarning.

### Deprecated

- Encrypted link files (`links.py`, `linkPath` in `config.js`) still work but will be removed
  in a later version. Use the password site or a link with a read pass.

### Not in this version

- A SpatialData table whose region is a labels element (a segmentation mask) is refused with
  the reason: save the cells as shapes instead.
- Cell segmentation of protein images: SODA reads a lab's own segmentation (`--cells`) but does
  not segment; without one, protein images are measured on a grid of 16 µm squares.
- Opening `.h5ad`, OME-Zarr or SpatialData files in the browser without converting them first.
