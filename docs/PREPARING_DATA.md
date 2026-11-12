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

`--src` is the Xenium Onboard Analysis output folder (the one holding `morphology_focus/`;
if you unpacked an archive with a `data/` folder inside, either folder works). What is read,
when present:

| File | Becomes |
|---|---|
| `morphology_focus/*.ome.tif` | the stain images (DAPI and the cell-boundary stains). Channel names come from the OME metadata. |
| `experiment.xenium` | the pixel size (µm per pixel), used for every position below |
| `cells.parquet` or `cells.csv.gz` | cell centroids (`x_centroid`, `y_centroid`, µm) and the median cell area |
| `cell_feature_matrix.h5` or `cell_feature_matrix/` (`matrix.mtx.gz`, `features.tsv.gz`, `barcodes.tsv.gz`) | counts per cell and gene |
| `transcripts.parquet`, else `transcripts.zarr.zip` | every decoded transcript, for the molecules layer |

- **Cells** become an expression layer called "Cells" (key `cells`), in the same format as
  Visium bins: each cell is drawn as a dot at its centroid, sized to the median cell area,
  and everything that works on bins works on it (gene maps, the whole-section UMAP and
  Leiden clusters, suggested names, the QC metrics, circle → UMAP and the volcano panel).
  The viewer says "cells" and "transcripts" instead of "bins" and "UMI". Only features of
  type "Gene Expression" are kept: negative-control probes and codewords, blank, unassigned
  and deprecated codewords are left out. Both files (cells table and matrix) are needed; if
  either is missing the layer is skipped with a message and the rest is still prepared.
- **Molecules**: every transcript with quality value ≥ 20 (`--qv-min`; the same cut Xenium
  uses for its count matrix), control codewords left out, is written to one small file per
  gene (4 bytes per molecule; format in [DATA_FORMAT.md](DATA_FORMAT.md#molecules-molecules)).
  In the viewer, the **Molecules** panel lets you pick genes; each picked gene's molecules
  are drawn as coloured dots on the image (up to 12 genes at once), and only the files of
  the picked genes are downloaded. Zoomed out, a random sample is drawn and the panel says
  so; zoomed in, every molecule is drawn. If there is no transcripts file the molecules are
  skipped with a message. Reading the zarr of a full slide (about 200 million transcripts)
  takes a while, see [How big and how long](#how-big-and-how-long).
- `--no-cells` / `--no-transcripts` skip either part, e.g. to re-prepare the image quickly.
- **Cell and nucleus outlines** (`cell_boundaries`, `nucleus_boundaries`) are not imported
  yet; the cells are drawn as dots of the median cell size.
- If `<src>/Custom_annotation.geojson` exists and you gave no `--geojson`, it is added as
  an overlay called "Annotations". It must be in image pixels (see
  [Region overlays](#region-overlays)).
- `--crop` and `--flip` are not available for `xenium`.
- Positions need no alignment: Xenium's µm coordinates start at the top-left corner of the
  morphology image, so pixel = µm ÷ pixel size.

For scale: a whole coronal mouse brain section (morphology_focus 54,002 × 40,815 px, 4 stains;
155,754 cells; 480 genes; 209 million transcripts) took 43 min to prepare on a loaded
workstation reading from a network share (41 min of it the image pyramid, 7 s the cells, 1.5 to
6 min the transcripts) and 7 min to export (tiles 2.4 min, section UMAP and Leiden 3.2 min). The
export is 3.85 GB in 12,679 files: tiles 2.78 GB, molecules 0.67 GB (480 files, largest 37 MB,
median 0.7 MB), cell gene files 0.30 GB, `hvg/` 0.10 GB. Opening it loads about 3.7 MB (11
files) before any gene is picked.

The writers in `expr.py` (`build_cell_layer`) and `molecules.py` (`MoleculeWriter`) take
plain arrays (centroids in µm, a cells x genes count matrix, gene ids and µm positions of
molecules), so a CosMx or MERSCOPE reader can reuse them; `xenium.py` is the Xenium reader.

## Stereo-seq (STOmics SAW)

```bash
python web/SpatialViewer/prepare.py stereoseq --id my_stereoseq --image path/to/SN_HE_regist.tif --gef path/to/SN.raw.gef --cellbin path/to/SN.cellbin.gef --out cache
```

Stereo-seq counts are per **DNB** (a capture spot; 0.5 µm apart on every chip so far), and
people read them as square **bins** of N × N DNBs: bin20 = 10 µm, bin50 = 25 µm,
bin100 = 50 µm, bin200 = 100 µm. `prepare.py stereoseq` makes those bin layers (like the
Visium HD bin layers) and, from the cell bin file, a cells layer. The files it reads are SAW
outputs (`<SN>` is the chip serial number; names as in STOmics' demo data, where they sit
flat in `outs/`; other SAW versions put the same files in subfolders, so give each path):

| Option | SAW file | Becomes |
|---|---|---|
| `--image` (required) | `<SN>_HE_regist.tif`, or `<SN>_ssDNA_regist.tif` / `<SN>_DAPI_regist.tif` | the image. It must be the **registered** image: one pixel per DNB, the size of the chip. RGB (H&E) is shown as a colour photo, a grey ssDNA / DAPI image as one white channel. |
| `--gef` | `<SN>.tissue.gef`, `<SN>.raw.gef` or `<SN>.gef` (bin GEF, HDF5) | the bin layers |
| `--gem` | `<SN>.tissue.gem.gz` or `<SN>.gem.gz` (text: `geneID geneName x y MIDCount ExonCount`) | the bin layers when there is no `--gef` (if both are given the GEF is read) |
| `--cellbin` | `<SN>.cellbin.gef` (or `<SN>.adjusted.cellbin.gef`) | the cells layer |
| `--tissue-mask` | `<SN>_<stain>_tissue_cut.tif` (same size as the image, non-zero = tissue) | which DNBs count as tissue, instead of the outline inside the GEF |

Not needed: `<SN>.report.tar.gz`, `visualization.tar.gz`, `*.barcodeToPos.h5`, the image
`.tar.gz` package, FASTQ files. Each part is optional except the image: without `--gef` /
`--gem` there are no bin layers, without `--cellbin` no cells layer. A file that is not found,
or is not the kind its option expects, is skipped with a message (a bin GEF and a cell bin GEF
given the wrong way round are swapped). A bin size the file cannot give stops the command
before the image is read.

**Positions.** The registered image has one pixel per DNB, with pixel (0, 0) at chip DNB
(0, 0), x to the right and y down. x and y in a file count DNBs from the file's own offset:
STOmics' writer stores x − minX and y − minY and keeps the offset in the file, as `minX` /
`minY` on a bin GEF's `expression` dataset, `offsetX` / `offsetY` on a cell bin GEF and
`#OffsetX` / `#OffsetY` in a GEM header. A whole-chip file has offset 0; a file cut to a
region has that region's corner. `prepare.py` adds the offset back (chip DNB = x + offset), so
a DNB sits at image pixel (chip x + 0.5, chip y + 0.5) and a bin's centre at its first DNB
+ N / 2, with no flip. The pixel size is the files' `resolution` (500 nm) = **0.5 µm per
pixel**; a GEM does not state it, so 0.5 µm is assumed (`--dnb-um` to change). `--offset X,Y`
is only for an image that was itself cropped (the chip DNB at its first pixel), and
`--image-scale` for an image that is not at one pixel per DNB. Records or cells that fall
outside the image are left out and counted in a warning: that means the image is not the
registered one. This was checked on whole-chip files and on files cut to a region, written by
STOmics' own library (gefpy 1.3.3).

**Bin sizes.** `--bins 20,50,100` is the default:

- **bin20 (10 µm)** is about a cell wide and close to Visium HD's 8 µm bins; it is the finest
  size the viewer offers (8 µm and up) and the largest layer.
- **bin50 (25 µm)** is the usual size for clustering Stereo-seq sections.
- **bin100 (50 µm)** is close to a classic Visium spot (55 µm): few bins, many counts each,
  for whole-section overviews.

Add `200` for 100 µm bins, or pass `--bins 50,100` to leave out the largest layer (about
40 % of the gene files and most of the clustering time). Bins are counted from the file's
origin (bin i covers file x from i × N to i × N + N − 1), as in SAW and Stereopy, so a bin
here is the same bin there. If the GEF already stores a wanted size (`geneExp/bin50`, as
`tissue.gef` and `<SN>.gef` do) it is used as it is; a size that is not stored is summed from
the largest stored size that divides it (bin1 in a `raw.gef`).

**Tissue.** `<SN>.raw.gef` holds the whole chip, including the DNBs beside the tissue (8 % of
the counts in the mouse brain used for testing); `<SN>.tissue.gef` and `<SN>.tissue.gem.gz`
hold only the DNBs under SAW's tissue mask. A whole-chip GEF is cut to the tissue here, using
`--tissue-mask` if given, else the tissue outline stored in the GEF (`contour/tissueContour`;
holes in the outline are left out). On the test chip that cut of `raw.gef` gave exactly the
records of `tissue.gem.gz`. A file with neither mask nor outline is kept whole, with a
message, and so is a GEF that is itself cut to a region (`minX` / `minY` not 0: taken as a
tissue GEF); `--tissue contour` applies the outline anyway and `--tissue none` keeps a whole
chip on purpose. The cut is per DNB when counts are summed from bin1; a stored coarser level
is kept bin by bin (a bin stays if any of its DNBs is tissue).

**Genes.** The gene name is shown (the gene id where there is no name). When several ids
share one name, the first keeps the name and the others are shown as `Name_<geneID>`, so no
counts are merged (a mouse annotation has about twenty such names).

**Cells.** Each cell of `cellBin/cell` becomes a dot at its centre, sized to the median cell
area; its counts are `cellBin/cellExp`. Everything that works on bins works on the cells layer
(see [Xenium](#xenium)). Cell outlines (`cellBin/cellBorder`) are not imported yet, and SAW's
own `clusterID` per cell is not read.

**Your own clusters.** The cache keeps a name per row so that
[your own UMAP and clusters](#your-own-umap-and-clusters-seurat-scanpy) can be attached with
`--user-layer bin050` (or `bin020`, `bin100`, `cells`): for a bin, `x_y` of its first DNB as
the file stores it (e.g. `5300_9500`); for a cell, the cell `id` of the cell bin GEF. Stereopy
names bins and cells `x * 2^32 + y` instead (one integer); prepare with `--barcodes stereopy`
to keep those names, so the `obs_names` of a Stereopy / AnnData object match as they are.

**Memory and time.** The bin1 table is never held whole: it is read gene by gene in chunks
of about 8 million records, each chunk is summed to (gene, bin) pairs and appended to a
temporary file per bin size inside the cache, and the layers are then written one at a time.
Peak memory is about 13 bytes per non-zero (gene, bin) pair of the largest layer. For a whole
mouse brain on a 1 cm chip (23,520 × 23,520 px image; 270 million bin1 records of which 248
million lie under the tissue; about 27,600 genes; about 130,000 cells), prepare with the
default bins took 5.3 min on a loaded workstation (2.3 min the image pyramid,
1.6 min reading the GEF, the rest writing the layers) with a peak of 2.7 GB.
Reading the same counts from `tissue.gem.gz` took 2.6 min for the text, with a peak of
3.9 GB. The export (53 min, 6.5 GB in 111,395 files) is described under
[How big and how long](#how-big-and-how-long). `--crop` and `--flip` are not available for
`stereoseq`.

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
  a few very bright squares cannot dominate. That value is what is shown, gated and ranked
  as cluster markers.
- **What the UMAP and clusters use: the mix of markers, not the brightness.** Every marker of
  a square rises and falls with its overall brightness (cell density, section thickness,
  the tissue edge), so on the values above the clusters became brightness bands. The UMAP
  instead gets the centred log-ratio of each square's background-free, percentile-scaled
  values (log(u + 0.1) minus the square's mean): which markers are up relative to the
  others. Panels with fewer than 10 embedded markers keep the asinh values (with so few, one
  bright marker moves the mean as much as brightness does). Then z-scores, PCA and UMAP as
  for genes. Squares with no marker above about 10%
  of its range (rims, gaps, empty squares: a quarter of one test slide) are left out of the
  UMAP and the clusters, and counted in the slide QC as low signal; Circle → UMAP leaves
  them out the same way.
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
   overlay exists, a region id per bin. Layers of 8 µm or more (and classic spots) also get
   per-bin QC values (`binqc/`: total counts, genes detected, % mitochondrial) for the
   viewer's Bin QC filter.
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

## Your own UMAP and clusters (Seurat, scanpy)

If you have already clustered the same bins or spots (in Seurat or scanpy, with cell-type
names), the export can show **your** results on a layer: by name in the Clustering menu, with
your category names in the legend and on hover, and your colours. SODA's Leiden clusterings
stay next to yours (yours first) unless you pass `--user-replace`. If you give clusters but
no embedding, they are drawn on SODA's UMAP.

What you give, for **one layer** (`--user-layer`, e.g. `008um`, `016um` or `spots`; the
barcodes must be that layer's):

| Option | File |
|---|---|
| `--user-umap umap.csv` | `barcode,UMAP_1,UMAP_2`: the first two numeric columns are used |
| `--user-obs clusters.csv` | `barcode,<column>,...`, repeatable. Text columns and integer columns with at most 30 values become clusterings; an integer column named like a clustering (`seurat_clusters`, `leiden`, `SCT_snn_res.0.8`, ...) or listed in `--user-columns` may have up to 200. Other numeric columns are skipped. A column `<name>_color` with `#rrggbb` values gives the colours of `<name>` |
| `--user-h5ad analysis.h5ad` | instead of the CSVs: obs columns, `obsm/X_umap` (or `--user-embedding <key>`) and the `uns/<column>_colors` palettes |
| `--user-columns a,b` | only these columns (default: every usable one) |
| `--user-barcode-prefix S1_` / `--user-barcode-suffix _2` | a merged object holding several sections: keep this section's barcodes |
| `--user-replace` | leave SODA's Leiden clusterings out |

```bash
python web/SpatialViewer/prepare.py export --id sample1 --out <cache> --dest <export> \
    --user-layer 008um --user-umap my_umap.csv --user-obs my_clusters.csv
```

On an export that already exists, `prepare.py userclusters` does only this step (same
options, plus `--id`, `--out`, `--dest`), in seconds to minutes instead of a whole export.

**From Seurat** (R; a Visium HD object loaded at one bin size):

```r
emb <- Embeddings(obj, "umap"); write.csv(emb, "my_umap.csv")
meta <- obj@meta.data[, c("seurat_clusters", "cell_type")]; write.csv(meta, "my_clusters.csv")
# colours (optional): meta$cell_type_color <- my_palette[as.character(meta$cell_type)]
```

**From scanpy** (Python):

```python
import pandas as pd
pd.DataFrame(adata.obsm["X_umap"][:, :2], index=adata.obs_names, columns=["UMAP_1", "UMAP_2"]).to_csv("my_umap.csv")
adata.obs[["leiden", "cell_type"]].to_csv("my_clusters.csv")      # or pass the .h5ad itself with --user-h5ad
```

**How rows are matched.** By barcode: exactly first, then the usual variants: a different or
missing `-1` style suffix (`-1`, `_1`, `-1-0`, `.1`) and a sample prefix (`S1_`, `slice1_`,
`S1:`). The export prints how many of your barcodes were found and by which rule.

- Under 50 % found: the export **stops** and shows examples of your barcodes and of the
  layer's (the usual cause is the wrong `--user-layer`, e.g. 16 µm clusters given for the 8 µm
  layer; it says so when it can tell).
- Some not found: a **warning**; the rest is imported.
- Bins you did not cluster (filtered out in QC, or not in your file) show as "not assigned",
  as SODA's own clusters show bins with too few counts.
- A merged object (two sections whose barcodes become the same once the sample tag is
  removed) stops with the tags it found: pick one with `--user-barcode-prefix` / `-suffix`.

Marker genes per group are computed by SODA from the counts, the same way as for its own
clusters, so the legend's marker lists, suggested names, the volcano plot's group pickers and
the QC panel work with your groups too.

**De-identification.** Column names and category values pass through `deid.py` as for
single-cell data: a column that names samples (`sample`, `orig.ident`, `donor`, ...) has every
value replaced by a number; barcodes and file names are never written.

**Caches made before this feature** do not keep barcodes; they are then read again from the
Space Ranger `outs` folder the cache was prepared from, which must still be reachable.
Re-running `prepare.py visium` stores them in the cache.

**Taking it back out.** Export again without the `--user-*` options: SODA's own UMAP and
clusters (kept as `umap/index_soda.json` and `coords_soda.bin`) are put back.

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
- **Bin QC filter.** On exports that carry `binqc/` (see [DATA_FORMAT.md](DATA_FORMAT.md)), a
  panel above Circle → UMAP shows a histogram of total counts, genes detected and % mitochondrial
  per bin, each with a min and a max slider (or type the number). It is **off by default**:
  nothing is removed until you move a slider, so results do not change by themselves. Bins
  outside a threshold are greyed on the slide and left out of the next run (PCA, UMAP,
  clusters, markers); the per-shape counts say how many each shape lost. The thresholds go into
  Save JSON (`qcFilter`) and the page link (`?qc=`), and Load JSON restores them (a file
  without them loads with the filter off). The whole-section UMAP made at export is not filtered.
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
- **Molecules** (Xenium): 4 bytes per transcript kept, one file per gene. A mouse brain
  section with a 480-gene panel (209 million transcripts, 167 million of them quality ≥ 20)
  gave 667 MB in 480 files, the largest gene 37 MB; the viewer downloads only the genes
  someone picks. Reading its `transcripts.zarr.zip` took about 6 minutes.
- **UMAP time** grows with the number of bins with ≥ 10 UMI and runs on one core.
  Expect it to dominate the export of a large Visium HD 8 µm layer.
- **Stereo-seq:** a whole mouse brain on a 1 cm chip (about 320 million counts under the tissue,
  27,600 genes, 131,000 cells) with the default bins exported in 53 min on one core, with a
  peak of 8.5 GB of memory: 6.5 GB in 111,395 files. bin20 (523,000 bins, 206 million non-zero
  values) is 2.6 GB and took 12 min to cluster, bin50 (84,000 bins) 1.9 GB and 1.5 min, bin100
  (21,000 bins) 1.2 GB and half a minute, the cells 0.6 GB and 2 min; the H&E tiles are 0.28 GB
  as JPEG. Each layer writes one file per gene that has counts (about 27,400 files; 4 to 10 min
  per layer onto a network drive). Opening the dataset loads 11 MB in 22 files before a gene is
  picked, most of it the bin positions and QC values of bin20; a widely expressed gene drawn on
  bin20 is 3 to 6 MB. `--bins 50,100` leaves out the largest layer.

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
