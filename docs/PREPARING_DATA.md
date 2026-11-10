# Preparing your data

Getting a dataset into the viewer takes two steps:

1. **`prepare.py <kind>`** reads the instrument output once and writes a **cache**:
   a multi-resolution copy of the image and a gene-by-bin index of the counts.
   This is the slow step, and the cache stays on your computer.
2. **`prepare.py export`** turns the cache into the **static file tree** the viewer reads
   (image tiles, per-gene files, the whole-section UMAP, and so on). This is what you upload.

Single-cell RNA-seq is the exception: `prepare.py rnaseq` writes the static tree directly.

All commands below are run from the repository folder, with the virtual environment
active (see the [README](../README.md#try-it-in-5-minutes)). Forward slashes in paths
work on Windows too. Each command is one line; paste it as is.

> **Before your first real dataset**, decide where the de-identification key lives and
> set `SPATIALVIZ_KEY` (see [DEIDENTIFICATION.md](DEIDENTIFICATION.md)). Otherwise the key
> file is created in whatever folder you run the command from.

## Options shared by every image kind

| Option | Meaning |
|---|---|
| `--id` | dataset id: the folder name in the cache and the export, and what appears in URLs. Use letters, digits, `_` and `-`. Do not put sample or patient identifiers in it. |
| `--out` | the cache folder (created if needed). Many datasets can share one. |
| `--label` | display name stored in the manifest. `config.js` can override it. |
| `--pixel-size` | µm per full-resolution pixel, for the scale bar. Read from the file when it carries one (OME-TIFF, Xenium, Space Ranger); give it otherwise. |
| `--geojson "Label=path.geojson"` | add an overlay (regions, annotations). Repeatable. Coordinates must be in full-resolution pixels of the final image. See [Region overlays](#region-overlays). |
| `--tile` | tile size in pixels, default 1024. Leave it. |
| `--rgb-split` | `akoya` and `xenium` only: keep RGB files as three channels instead of one. `tiff` and `visium` always split RGB. |
| `--flip h\|v\|hv` | `akoya`, `tiff`, `visium`: mirror the image (and, for Visium, the bins with it). |
| `--crop x0,y0,x1,y1` | `akoya`, `tiff`, `visium`: keep only that rectangle, in full-resolution pixels. To split one slide into several samples, run once per rectangle with a different `--id`. |

The output is `<out>/<id>/` with `manifest.json`, `image.zarr/` and, for Visium,
`expr/`. The command prints `[done] ...manifest.json` when it finishes.

---

## Visium HD

You need Space Ranger's `outs/` folder and **the same full-resolution microscope image
you gave `spaceranger count --image`** (H&E or immunofluorescence). The bin positions in
`outs` are in that image's pixel grid.

The demo command (see the [README](../README.md#try-it-in-5-minutes)):

```bash
python web/SpatialViewer/prepare.py visium --id demo_visiumhd --image demo/raw/visiumhd/image.tif --outs demo/raw/visiumhd/outs --out demo/cache
```

For your data:

```bash
python web/SpatialViewer/prepare.py visium --id mouse_brain_hd --label "Mouse brain" --image path/to/microscope_image.tif --outs path/to/spaceranger/outs --out cache --bins 008um,016um
```

- It reads every `outs/binned_outputs/square_<N>um/` that has
  `filtered_feature_bc_matrix.h5` and `spatial/tissue_positions.parquet`. Each becomes
  one expression layer named `002um`, `008um`, `016um` ...
- `--bins 008um,016um` keeps only those layers. **Recommended:** the 2 µm layer is by far
  the largest, the viewer does not offer layers below 8 µm by default (`exprMinBinUm` in
  `config.js`), and at 2 µm nearly every bin holds a single transcript.
- Only features of type "Gene Expression" are kept.
- An RGB (H&E) image becomes three channels, Red/Green/Blue, shown as a colour photo.
  A multichannel IF image keeps its channels.

## Visium (classic, 55 µm spots)

Same command. With no `binned_outputs/` folder it reads `outs/filtered_feature_bc_matrix.h5`
and `outs/spatial/tissue_positions.parquet` (or `tissue_positions.csv` /
`tissue_positions_list.csv`) as one layer called `spots`.

```bash
python web/SpatialViewer/prepare.py visium --id lymph_node --image path/to/fullres_image.tif --outs path/to/spaceranger/outs --out cache
```

Give the full-resolution image, not `tissue_hires_image.png` (that one is downscaled; if
you have nothing else, see `--image-scale` under [Common problems](#common-problems)).

## Xenium

```bash
python web/SpatialViewer/prepare.py xenium --src path/to/xenium_output_folder --id my_xenium --out cache
```

- Reads the morphology stains from `<src>/morphology_focus/*.ome.tif` (the layout written
  by recent Xenium Onboard Analysis versions). Channel names and pixel size come from the
  OME metadata.
- If `<src>/Custom_annotation.geojson` exists and you gave no `--geojson`, it is added as
  an overlay called "Annotations". It must be in image pixels (see
  [Region overlays](#region-overlays)).
- **Transcripts and cell boundaries are not imported.** For Xenium the viewer currently
  shows the morphology images and overlays only.
- `--crop` and `--flip` are not available for `xenium`.

## Multiplexed immunofluorescence

### One multichannel file (OME-TIFF, QPTIFF, ImageJ or other TIFF stack)

Look inside first. This prints the axes, size, data type, pyramid levels and channel names
without building anything:

```bash
python web/SpatialViewer/prepare.py inspect --src path/to/slide.ome.tif
```

Then:

```bash
python web/SpatialViewer/prepare.py tiff --id tonsil_panel --src path/to/slide.ome.tif --out cache
```

The demo uses the same command:

```bash
python web/SpatialViewer/prepare.py tiff   --id demo_multiplex --src demo/raw/multiplex.ome.tif --out demo/cache
```

- Channel names come from the file. If they are missing or wrong, set them in order with
  `--channels "DAPI,CD3,CD8,CD20,PanCK,Ki67"`.
- Any non-spatial axis (C, I, Z, Q) is taken as the channel axis. 8-bit and 16-bit images
  are kept at their depth; other types are stored as 16-bit.
- Channels whose names start with `open`, `blank`, `empty` or `unused` are hidden by the
  viewer (`hideChannels` in `config.js` changes this); the nuclear stain (DAPI / Hoechst)
  is pinned to the top of the list.

### A folder with one TIFF per marker (Akoya-style export)

```bash
python web/SpatialViewer/prepare.py akoya --id tonsil_panel --src path/to/marker_folder --out cache --pixel-size 0.5
```

Files named like `NN_Marker_<folder name>_gray.tif` are read in order of `NN`, and the
marker name is taken from between. `_gray` files are preferred; `_red`/`_green`/`_blue`/`_rgb`
renderings are used only when a marker has no grey file. All files must have the same
width and height.

On `export`, a multi-channel fluorescence image (TIFF or Akoya, with a pixel size) also gets a
grid layer: the mean of each marker per 16 µm square (`--grid-um`; 8 gives about four times the
squares and files), with a UMAP and clusters. The means are taken from full-resolution pixels
(the coarser pyramid levels keep the brightest pixel for display and would overstate markers).
Only tissue squares are kept. The tissue mask: DAPI (or Hoechst) square means, averaged over
100 µm (7 x 7 squares, `--tissue-smooth-um`), Otsu threshold on log1p, enclosed holes of up to
4 squares filled. `--tissue-smooth-um 0` thresholds the raw square means instead, which keeps
only nucleus-dense squares and drops axon and collagen stretches. The rule and window are
recorded in the layer's `grid` block in `manifest.json` and in `umap/index.json`.

Everything after that is automatic, so a slide needs no hand-tuning:

- **Channel cleanup.** Channels named open / blank / empty / unused / none are dropped. Two
  channels with the same marker name (ignoring case, `01_` prefixes and dye tags such as
  `-488` or `FFPE`) are reduced to the one with the stronger signal. A channel that is flat
  over the tissue is kept in the files but left out of the UMAP; so is DAPI.
- **Normalisation, per channel and per slide** (pieces on one slide were stained together):
  the background (median of the off-tissue squares) is taken off, the result is divided by
  the channel's 99.5th percentile over the tissue (folds and other everything-bright squares
  left out of that percentile), then put through asinh (cofactor 0.2), so
  a few very bright squares cannot dominate. Then z-scores, PCA and UMAP as for genes.
- **Slide QC** (`expr/<layer>/qc.json`, shown in the viewer under the layer): weak stains
  (signal under 15x background: kept, flagged), saturated stains (more than 1% of tissue
  pixels at the detector's limit: kept, flagged), squares where nearly every marker is bright
  (folds, bubbles, bright edges; a cluster made mostly of them is suggested as
  "Artefact? (fold / edge)"), squares with many saturated pixels, and a check for tile seams.
  UMAPs are per slide: they are not comparable across slides.
- **Tile seams.** A stitched slide can show steps where the camera fields meet. On our s48
  (fields of 960 x 720 px) there are three parts: a dark line along the last pixel column / row
  of every field; shading shared by every field (one side of every field brighter than the
  other, so the step goes the same way at most edges); and brightness / background differences
  between fields. The seams are always detected and reported. With `--seam-correct`, when the
  field grid is found, the border line is measured and left out, the
  shared shading is measured and reported, and every marker with steps gets one brightness
  factor (and, if it helps, a background offset) per field, fitted from the pixels on both sides
  of each field edge. It is kept only if it at least halves the step at field edges without
  raising it inside the fields; the squares are then recomputed from corrected values. Shared
  shading is **not** corrected (that needs a flat-field correction, not built yet): markers
  dominated by it are listed as not correctable, with the reason, and stay in the UMAP. The
  image itself is shown as recorded. The correction is off by default: on s48 it removed the
  steps at field edges in 11 of 24 markers but made no measurable difference to the squares,
  because most of the seam there is the shared shading. Details: `web/SpatialViewer/seams.py`.
- **Large slides.** Above 200,000 tissue squares the UMAP is fitted on a seeded random
  sample of 200,000 and the other squares are placed into it (`--umap-sample N`; 0 = fit on
  every square, slower).

## H&E or other brightfield image on its own

```bash
python web/SpatialViewer/prepare.py tiff --id he_section --src path/to/he_image.tif --out cache --pixel-size 0.25
```

An RGB image becomes three channels shown together as a photo.

## Single-cell / single-nucleus RNA-seq

The UMAP and clusters must already exist: the viewer shows your analysis, it does not
redo it. Output goes straight to the export folder (`--dest`); there is no cache and no
`export` step.

From an AnnData file (scanpy; or Seurat converted with SeuratDisk or sceasy):

```bash
python web/SpatialViewer/prepare.py rnaseq --id pbmc_cells --label "All cells" --h5ad path/to/cells.h5ad --dest export
```

- The embedding is read from `obsm["X_umap"]` (change with `--embedding X_tsne` etc.).
- Gene names come from the first of the `var` columns `gene_symbols`, `feature_name`,
  `gene_name`, `gene_names`, `features`, `name`, else the `var` index.
- `obs` columns with 2-200 distinct values become colour-by categories; numeric columns
  are kept as numbers. The viewer colours the UMAP by the first column listed in
  `rnaColorBy` in `config.js`, else the first column whose name contains
  `cluster`, `leiden` or `louvain`.
- Integer counts are normalised to log1p(counts per 10,000); anything else is shown as
  provided.
- Files over 2 GB are streamed instead of loaded: it uses `raw/X` when present, keeps
  only categorical `obs` columns plus common QC columns, and needs free disk space of about
  8 bytes per non-zero value for temporary files (put them elsewhere with `--tmp path`).

From a 10x matrix folder plus CSV files:

```bash
python web/SpatialViewer/prepare.py rnaseq --id pbmc_cells --mtx path/to/filtered_feature_bc_matrix --umap path/to/umap.csv --obs path/to/clusters.csv --dest export
```

- `--mtx`: folder with `matrix.mtx`, `features.tsv` (or `genes.tsv`) and `barcodes.tsv`,
  gzipped or not.
- `--umap`: CSV with a header row, then `barcode,UMAP1,UMAP2`. Cell Ranger's
  `analysis/umap/gene_expression_2_components/projection.csv` works. Cells without a row
  are left out.
- `--obs`: CSV with a header, then `barcode,<column>,...`. Repeatable. A column called
  `cluster` or `clusters` is renamed after its folder, so Cell Ranger's
  `analysis/clustering/gene_expression_graphclust/clusters.csv` becomes `graphclust`.

The dataset id, label, column names and values pass through de-identification. If the id
contained something that looks like a sample ID, the output folder gets the replaced name;
the last log line shows the real path.

---

## Export

```bash
python web/SpatialViewer/prepare.py export --id demo_visiumhd  --out demo/cache --dest demo/export --markers examples/demo_markers.json
python web/SpatialViewer/prepare.py export --id demo_multiplex --out demo/cache --dest demo/export --fmt jpg --markers examples/demo_protein_markers.json
```

| Option | Meaning |
|---|---|
| `--id` | the dataset to export |
| `--out` | the cache folder `prepare.py` wrote into |
| `--dest` | where to write `<dest>/<id>/` |
| `--fmt png\|jpg` | tile format. `png` (default) is lossless. `jpg` is several times smaller and fine for viewing; it slightly shifts very dark pixels. Use `png` when exact intensities matter. |
| `--workers` | parallel tile writers, default 8 |

What it writes, per dataset (layout in [DATA_FORMAT.md](DATA_FORMAT.md)):

1. Image tiles for every channel and zoom level, then a check that one tile per channel
   matches the cache (`[export] verified N channel(s)`).
2. Per expression layer: the gene list, bin positions, one file per gene, and, if a region
   overlay exists, a region id per bin.
3. Per layer with bins of 8 µm or more (or classic spots): the data for the circle → UMAP
   tool (`hvg/`) and the whole-section UMAP with Leiden clusters (`umap/`). See below.
4. Overlays and the manifest, without your local file paths.

**Re-running** skips tiles and gene files that already exist, so an interrupted export
can simply be started again. The flip side: after you re-run `prepare.py` on a dataset
(new crop, new image), **delete `<dest>/<id>/` before exporting again**, or old tiles stay.

## Clustering and the whole-section UMAP

During export, every expression layer with bins of 8 µm or more (or classic spots) gets a
UMAP of all its bins and Leiden clusters. In the viewer this drives the
**Tissue | Side by side | UMAP** switch, colouring by cluster, and the marker gene lists.

What it computes, in order (scanpy's defaults, written out in `sectionumap.py`):

1. keep bins with at least 10 UMI;
2. normalise to log1p(counts per 10,000);
3. pick the 2,000 most variable genes over the whole section (dispersion, z-scored within
   20 bins of mean expression);
4. scale each gene, clip at 10, PCA to 30 components (exact, deterministic);
5. UMAP with umap-learn: 15 neighbours, `min_dist` 0.3, seed 42;
6. Leiden clustering (leidenalg, RB configuration) on UMAP's neighbour graph at
   resolutions 0.3, 0.6 and 1.0 (0.6 is shown first), 2 iterations, seed 42.
   Cluster 0 is the largest;
7. markers: per cluster, the top 10 genes by log2 fold change of mean expression against
   all other bins, among genes seen in at least 25 % of the cluster's bins.

Everything is seeded, so the same cache gives the same files. The UMAP runs on one core
because seeding umap-learn makes it single-threaded; on large layers this is the slowest
part of the export.

**It needs the optional packages:**

```bash
pip install -r web/SpatialViewer/requirements-clustering.txt
```

Without them the export logs `[umap] <layer>: skipped - <package> is not installed` and
carries on; everything else, including the circle → UMAP data, is still written.

**Turning it off.** Add switches to `prepare.py export`:

| Switch | Skips |
|---|---|
| `--no-umap` | the whole-section UMAP and Leiden clusters (and so the suggested names) |
| `--no-celltypes` | only the suggested cell-type names |
| `--no-hvg` | the gene chunks the circle → UMAP tool reads |

**Recomputing.** An existing `expr/<layer>/umap/` made with the same parameters is reused.
Delete that folder to compute it again.

## Suggested cluster names (marker lists)

A cluster number says nothing about what the cluster is. After clustering, the export
scores every bin against lists of known marker genes and suggests a name for each cluster,
with a confidence you can see, accept or overrule.

By default the export uses general lists that work for any tissue:
`web/SpatialViewer/markers/rna_generic.json` (genes, for Visium) and
`web/SpatialViewer/markers/protein_generic.json` (protein panels). They name
neighbourhoods rather than single cells ("T-cell rich", "Vascular", "Neural / glial"), because
a bin or square holds several cells. They were written from widely used textbook markers, not
copied from a marker database. The protein list also has an alias table (CD45 = PTPRC,
Pan-Cytokeratin = panCK, NF-L = NEFL, ...), and channel names are matched ignoring case, `01_`
prefixes and dye tags; a protein type needs at least 2 of its markers on the panel (CD45,
lymphatic and Ki67 types need 1), otherwise it is listed as absent.

Tissue-specific lists can be passed instead: `markers/nerve.json` has marker genes for twelve
cell types found in and around peripheral nerve (myelinating and non-myelinating Schwann cells,
perineurial cells, fibroblasts, endothelial cells, pericytes / smooth muscle, macrophages,
mast cells, T / NK cells, B / plasma cells, adipocytes, skeletal muscle), and
`markers/akoya_nerve.json` the protein version. **Edit them for your tissue**, or write your
own file in the same shape and pass it with `prepare.py export --markers my_markers.json`.
Gene names are matched without regard to case, so one file serves human and mouse. In the
browser the circle tool uses the list the export used for protein layers, and for genes the
file named by `markerSet` in `config.js` (default `markers/rna_generic.json`).

How a name is suggested:

1. For each type, a bin's score is the mean, over the type's genes, of the gene's
   log1p(counts per 10,000) expressed as a z-score across the bins.
2. A gene detected in fewer than 0.5 % of bins is left out as too rare to score (one stray
   count would otherwise decide the bin). A type needs at least 2 measurable genes;
   otherwise it is listed as "not scored" with the reason.
3. A bin's type is its best-scoring type if that score is at least 0.25, otherwise none.
4. A cluster's suggestion is the type most of its bins have. Confidence comes from that
   share: **high** (at least 50 %, and at least twice the runner-up), **medium** (at least
   30 %), **low** (at least 15 %, shown as "maybe ..."), otherwise **unclear**. The runner-up
   and its share are shown too, because a bin of 8–16 µm usually holds parts of several
   cells and many clusters are mixtures.

In the viewer each cluster row shows the suggestion with **✓ Accept**; **Accept all high**
takes every high-confidence one; a name typed with ✎ always wins. "Colour bins by →
Suggested cell type" paints every bin by its best type. The "marker lists" link shows which
genes were found, which were too rare and which are absent from the data.

These are suggestions from marker genes, not a cell-type annotation: check them against the
cluster's own marker genes and where it sits in the tissue. The same rule runs in the
browser on circle → UMAP clusters, but there only the section's 2,000 most variable genes
are available, so some markers are reported as absent.

## The circle → UMAP tool

On an exported Visium / Visium HD dataset, you can draw circles or freehand shapes on the
tissue and get a UMAP of just the bins inside them, computed in the browser. It uses the
`hvg/` files written at export.

- Available for layers with bins of 8 µm or more (and classic spots). The 2 µm layer is not
  offered: almost every bin there holds one transcript, too sparse to embed.
- Only the section's 2,000 most variable genes are available to it (chosen at export over
  the whole section). Inside the shapes it keeps the 1,000 of those that vary most.
- Bins with fewer than 10 UMI are dropped. "One per bin" needs at least 50 bins.
  "One per shape (summed)" needs at least 3 shapes and shows PC1 vs PC2 instead of a UMAP
  when there are fewer than 15.
- At most **50,000 bins** are embedded. Above that it takes a seeded sample, stratified by
  shape, and says so. The page estimates about one second per 1,000 bins.
- A bin inside two shapes belongs to the one drawn last.
- Leiden clusters of the result at a resolution you choose (0 = none).
- **Marker squares (Akoya and other multi-channel images)** use the `mat/` files instead, and
  add a comparison with the rest of the slide, because a UMAP of only the shapes cannot show
  what the shapes have in common (a marker high in every square drawn, such as CD45 over an
  immune-rich area, differs nowhere inside them):
  - **Each shape vs the rest of the slide.** Under each shape, the markers most raised and
    lowered compared with every other tissue square of the slide; below the shapes, the full
    table for one shape (or all shapes together): mean intensity inside and in the rest, their
    ratio, and z (the difference on the scale the UMAP uses, in the slide's standard
    deviations). Every square of the shape counts, also above the 50,000 sample. **Compare
    with slide** does only this, without a UMAP; **Save table (CSV)** writes every shape and
    marker.
  - **Cluster markers against the slide or against the other drawn squares.** A switch in the
    UMAP panel. Against the slide (the default) lists what each cluster has more of than the
    slide as a whole; against the other drawn squares lists what tells the clusters apart.
  - **Suggested names** are scored against the whole slide, so an all-immune circle is named
    immune.
- "Save JSON" stores the shapes, layer, seed and every parameter so a result can be redone;
  "Save CSV" saves bin row, shape and UMAP coordinates (no barcodes).
- It works with exported data only, not with the live tile server.

## Region overlays

`--geojson "Cortex regions=path/to/regions.geojson"` adds a GeoJSON overlay, for example
regions drawn in QuPath. Rules:

- Coordinates in full-resolution pixels of the image you gave `prepare.py`. If your tool
  writes micrometres, divide by the pixel size first.
- `--crop` and `--flip` are **not** applied to overlays. If you crop or flip, the overlay
  must already be in the cropped / flipped frame.
- The first overlay with polygons becomes the region set: at export every bin is assigned
  to the smallest polygon containing it, and the viewer can then filter genes by region.
  A polygon named `hemi`, or covering 85 % or more of the overlay's extent, is treated as a
  whole-section outline, not a region. Region names come from QuPath's
  `classification.name`, or `name` / `label` / `region` / `id` properties.

---

## How big and how long

We have not benchmarked preparation or export on a reference machine, so there are no
timing figures here. What can be said from the formats:

- **Prepare** time is dominated by reading the source image once. Large uncompressed
  multi-gigabyte TIFFs read fastest from a local SSD.
- **Tiles:** about 1.33 × channels × (width / 1024) × (height / 1024) files. For example,
  a 40,000 × 30,000 pixel slide with 30 channels gives roughly 48,000 tiles. Mostly-black
  fluorescence tiles are small; JPEG is smaller than PNG.
- **Expression**, per layer: 8 bytes per bin for positions, 12 bytes per non-zero
  (bin, gene) pair across the gene files, one file per expressed gene (often 15,000-20,000
  files).
- **`hvg/`**: about 12 bytes per bin plus 4 bytes per non-zero value among the 2,000 genes,
  in one file per 500 µm square of tissue.
- **`umap/`**: 8 bytes per bin, plus 2 bytes per bin for each of the three clusterings.
- **UMAP time** grows with the number of bins with ≥ 10 UMI and runs on one core.
  Expect it to dominate the export of a large Visium HD 8 µm layer.

Run the [checklist in HOSTING.md](HOSTING.md#before-you-upload-checklist) to count files
and bytes before you upload.

## Common problems

**Bins do not line up with the tissue.** The positions in Space Ranger's output are in
the pixel grid of the image given to `spaceranger count --image`. Check, in this order:

| What you see | Cause | Fix |
|---|---|---|
| The pattern of bins is right but too large or too small, anchored at the top-left | your image has a different resolution than the one Space Ranger used | `--image-scale <your image width ÷ Space Ranger's image width>`, e.g. `0.5` for an image downsampled by 2 |
| The bins are a mirror image of the tissue | the image was flipped relative to the one Space Ranger used | `--flip h` (left-right) or `--flip v` (top-bottom); both image and bins are mirrored |
| You only want part of the slide, or the image has a margin Space Ranger did not have | | `--crop x0,y0,x1,y1`; bins move with the crop, bins outside it are dropped |
| Close but not matching, or rotated | a different scan or a rotated export of the slide | use the exact file Space Ranger was given. Rotation is not supported. |

**Windows: `FileNotFoundError` while the UMAP starts.** umap-learn's compiler writes cache
files inside the Python environment. If that environment sits in a deeply nested folder the
full path can pass Windows' 260-character limit. Keep the virtual environment at a short
path such as `<repo>\.venv`, or enable long paths in Windows.

**`[umap] ...: skipped - umap is not installed`.** Install
`requirements-clustering.txt`, delete nothing, and run the export again: existing tiles and
gene files are skipped, so it only adds the missing part.

**`[export] *** N channel(s) do not match the pyramid`.** Some tiles on disk differ from the
cache (an earlier export of an older cache, or a disk problem). Delete
`<dest>/<id>/tiles/` (at least for the channels listed) and export again.

**`[deid] REFUSING to publish ...: sample IDs still present`.** Something about to be written
looks like a sample ID (see [DEIDENTIFICATION.md](DEIDENTIFICATION.md)). It can be a gene
name that happens to match the pattern. Read the listed values; rename the source if it
really is an identifier.

**`[visium] WARNING: no count layers found`.** `--outs` must point at the `outs` folder
itself, which contains `binned_outputs/` (HD) or `filtered_feature_bc_matrix.h5` (classic).

**`does not match` (akoya).** The marker files are not all the same size; they must come
from the same export.

**Wrong or missing channel names.** Use `prepare.py inspect`, then `--channels` (`tiff`).

**No scale bar.** Give `--pixel-size`.

**Running out of memory** on Visium HD: leave out the 2 µm layer with `--bins 008um,016um`.
