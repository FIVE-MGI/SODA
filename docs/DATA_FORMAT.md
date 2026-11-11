# Data format

This is the layout that `prepare.py export` (images, spatial expression) and
`prepare.py rnaseq` (single-cell UMAPs) write, and that the viewer reads. It is derived
from `export.py`, `hvgchunks.py`, `binqc.py`, `sectionumap.py`, `regions.py`, `rnaseq.py`,
`expr.py`, `molecules.py`, `xenium.py`, `viewer.js` and `qc.js`. You only need it if you want to write the files with your own tools or
read them from another program.

General rules:

- One folder per dataset, named by its id: `<root>/<id>/`. The viewer fetches
  `<root>/<id>/manifest.json` first and everything else relative to that folder.
- All binary files are **little-endian**, with no header unless one is described below.
  "float32" is IEEE 754 single precision.
- **Coordinates** of bins, cells, molecules, regions and overlays are in **full-resolution image pixels**
  of the image the pyramid was built from (after any `--crop` and `--flip`): x to the
  right, y down, origin at the top-left corner of the image.
- A **row** means a position in `expr/<key>/positions.bin`. `region.bin`, `binqc/`, `hvg/`
  and `umap/` all refer to bins by row. No file contains a barcode.
- File names made from gene names use JavaScript's `encodeURIComponent` (Python:
  `urllib.parse.quote(name, safe="-_.!~*'()")`). Known issue: for a name that actually
  changes under encoding (a space, `/`, `|`, `+` ...), the file on disk contains a literal
  `%`, while a web server decodes the `%` in the request, so that gene is not found.
  Names made of letters, digits and `- _ . ! ~ * ' ( )` (nearly all gene symbols) are
  unaffected.

```
<id>/
  manifest.json
  tiles/<c>/<z>/<x>/<y>.png              (or .jpg)
  overlays/<NN>_<label>.geojson
  expr/<key>/genes.json
  expr/<key>/positions.bin
  expr/<key>/gene/<NAME>.bin
  expr/<key>/region.bin                  only if an overlay has polygons
  expr/<key>/binqc/{total.bin, genes.bin, pct_mt.bin}   per-bin QC values (RNA layers, exports from 11 Nov 2026)
  expr/<key>/hvg/{genes.json, index.json, <tx>_<ty>.bin}
  expr/<key>/mat/{index.json, <tx>_<ty>.bin}   protein grid layers only (instead of hvg/)
  expr/<key>/umap/{index.json, coords.bin, clusters_r<res>.bin}
  molecules/index.json                   imaging platforms (Xenium): every transcript, per gene
  molecules/gene/<NAME>.bin
<root>/about.html                        optional, shown by "Read more about this project"
```

---

## manifest.json (imaging and spatial expression)

```json
{
 "id": "demo_visiumhd",
 "label": "demo_visiumhd",
 "kind": "visium",
 "pixelSizeUm": 0.2737,
 "created": "2026-11-02T10:00:00",
 "expr": [ ... ],
 "overlays": [ {"key": "00_Regions", "label": "Regions", "file": "00_Regions.geojson"} ],
 "image": {
   "width": 20000, "height": 18000, "levels": 6, "tileSize": 1024, "dtype": "uint8",
   "levelShapes": [[18000, 20000], [9000, 10000], ...],
   "channels": [
     {"index": 0, "name": "Red", "color": "#ff0000", "bits": [0, 255], "window": [0.0, 255.0], "dataMax": 255}
   ]
 },
 "static": true,
 "tileFormat": "png"
}
```

(The numbers are illustrative.)

| Field | Meaning |
|---|---|
| `id` | dataset id = folder name |
| `label` | display name (the `--label` given to `prepare.py`, else the id) |
| `kind` | `akoya`, `tiff`, `xenium` or `visium` (the `prepare.py` subcommand used) |
| `pixelSizeUm` | µm per full-resolution pixel, for the scale bar; `null` if unknown |
| `created` | when the cache was prepared, local time, ISO 8601 |
| `image.width`, `image.height` | full-resolution size in pixels |
| `image.levels` | number of pyramid levels; the coarsest fits in one tile |
| `image.tileSize` | tile edge in pixels (default 1024) |
| `image.dtype` | `uint8` or `uint16`: the type of the source image (other types are stored as uint16) |
| `image.levelShapes` | `[height, width]` per level, **index 0 = full resolution**, each level half the previous (rounded up) |
| `image.channels[].index` | channel number used in tile paths |
| `image.channels[].name` | marker / stain name; `Red`, `Green`, `Blue` for an RGB photo |
| `image.channels[].color` | default display color `#rrggbb` |
| `image.channels[].bits` | `[lo, hi]` in raw units: how raw values were mapped to the 8-bit tiles (below) |
| `image.channels[].window` | default display window `[lo, hi]` in 8-bit tile units |
| `image.channels[].dataMax` | largest raw value seen at a mid-resolution level |
| `crop` | present when `--crop` was used: `{"box": [x0, y0, x1, y1], "ofFullSize": [w, h]}` |
| `flip` | present when `--flip` was used: `h`, `v` or `hv` |
| `expr` | spatial expression layers (Visium bins, Xenium cells, protein grid squares; empty list otherwise), see below |
| `molecules` | only when single molecules were imported (Xenium): `{"dir": "molecules", "format": "xy-uint16", "n", "nGenes", "bytes", "qvMin", "platform"}`, see [Molecules](#molecules-molecules) |
| `overlays` | GeoJSON overlays: `key` (file stem), `label` (from `--geojson "Label=path"`), `file` |
| `static` | always `true` in an export |
| `tileFormat` | `png` or `jpg`: extension of every tile |

The prepare cache's own manifest also has a `source` block with local file paths;
`export` removes it.

### `expr[]` entries

| Field | Meaning |
|---|---|
| `key` | layer id and folder name: `002um`, `008um`, `016um` ... for Visium HD bins, `spots` for classic Visium, `cells` for Xenium cells, `grid016um` for a protein grid |
| `label` | e.g. `8 µm bins`, `Visium spots (55 µm)`, `Cells` |
| `kind` | absent for sequencing bins; `cells` for segmented cells (Xenium), `protein` for a protein grid |
| `nBins` | barcodes (cells) in the count matrix (including any without a position) |
| `nGenes` | genes of feature type "Gene Expression" |
| `binSizeUm` | bin edge in µm; `null` for classic spots and cells |
| `spotDiameterPx` | spot / bin size in full-resolution pixels; for cells the diameter of a circle of the median cell area (`cellDiameterUm`, 10 µm when the area is unknown) |
| `cellDiameterUm` | cells only: that diameter in µm |
| `platform` | cells only: `xenium` |
| `micronsPerPixel` | from Space Ranger's `scalefactors_json.json`, divided by `--image-scale` |
| `imageScale` | the `--image-scale` used |
| `regions` | only with a region overlay: `[{"id": 1, "name": "Cortex", "nBins": 5123}, ...]`, ids as in `region.bin`, regions with no bins left out |
| `regionOverlay` | only with a region overlay: the overlay `file` the ids refer to |
| `binQc` | only if `binqc/` was written (RNA layers of 8 µm or more, and classic spots; exports from 11 November 2026): the files and ranges of the per-bin QC values, see below |
| `hvg` | only if `hvg/` was written: `{"nGenes", "chunkPx", "chunks", "bytes"}` |
| `mat` | protein grid layers (`"kind": "protein"`) only, if `mat/` was written: `{"nChannels", "chunkPx", "chunks", "bytes"}` |
| `umap` | only if `umap/` was written: `{"nEmbedded", "clusterings": [{"resolution", "n"}], "default"}` |

## Tiles: `tiles/<c>/<z>/<x>/<y>.<tileFormat>`

- `c`: channel index.
- `z`: zoom, `0` = coarsest level, `levels - 1` = full resolution. Pyramid level
  `L = levels - 1 - z`, and `levelShapes[L]` is that level's size.
- `x`, `y`: tile column and row at that level. Tile `(x, y)` covers level pixels
  `[x*T, (x+1)*T) × [y*T, (y+1)*T)` with `T = tileSize`, which is full-resolution pixels
  `[x*T*2^L, ...)`.
- Every tile is a **single-channel 8-bit greyscale** image, exactly `T × T`. Tiles at the
  right and bottom edges are padded with black (0). Tiles are written for every
  `(c, z, x, y)` inside the level; nothing outside it.
- 8-bit value = `clip((raw - bits[0]) * 255 / (bits[1] - bits[0]), 0, 255)`. For uint16
  images `bits[1]` is the 99.99th percentile of non-zero pixels, so the brightest
  0.01 % saturate. For uint8 images `bits` is `[0, 255]` (unchanged).
- Coarser levels are made by 2×2 **maximum** for fluorescence (so small bright objects
  survive zooming out) and 2×2 **mean** for an RGB photo.
- PNG is written with zlib level 1; JPEG at quality 88.
- Colour, window and blending are applied in the browser on the GPU; tiles carry no colour.

Number of tile files ≈ channels × Σ over levels of ⌈width_L / T⌉ × ⌈height_L / T⌉,
roughly 1.33 × channels × (width / 1024) × (height / 1024).

## Overlays: `overlays/<NN>_<label>.geojson`

GeoJSON `FeatureCollection`s, copied from `--geojson "Label=path"` (and, for Xenium,
`Custom_annotation.geojson`). `NN` is the order given; the label is reduced to
`[A-Za-z0-9_.-]`. Coordinates must be in full-resolution image pixels of the final
image. `prepare.py` copies them unchanged: it does **not** apply `--crop` or `--flip`
to overlays.

The **region overlay** is the first overlay that contains Polygon / MultiPolygon features.
In the export each of its polygon features gets either

- `properties.__region_id`: integer ≥ 1, the id used in `region.bin`, or
- `properties.__region_outline: true`: a whole-section outline, not a region. A feature
  is an outline if its name is `hemi` or its area is ≥ 85 % of the overlay's bounding box.

A feature's name is `properties.classification.name` (QuPath style), else
`properties.name`, `label`, `region` or `id`. A bin belongs to the **smallest** polygon
that contains its centre.

## Spatial expression: `expr/<key>/`

### `genes.json`

```json
[{"n": "Gapdh", "t": 123456.0}, ...]
```

Every gene of the layer in matrix order: `n` name, `t` total UMI count over all bins.
Genes with no counts are listed but have no `gene/` file (the request returns 404).

### `positions.bin`

float32 pairs `x, y`, one per bin that has a position, in count-matrix barcode order.
Number of rows = file size / 8. This order defines **row** numbers.

### `gene/<NAME>.bin`

float32 triples `x, y, value` for each bin where the gene's count is > 0 (and the bin has a
position), ordered by bin. `value` is the raw UMI count. Number of points = file size / 12.

### `region.bin`

uint16 per row (same length and order as `positions.bin`): region id, `0` = in no region.
Ids match `regions[]` in the manifest and `__region_id` in the region overlay.

### `binqc/`: per-bin QC values (total counts, genes detected, % mitochondrial)

Written by `binqc.py` for every sequencing layer with bins of 8 µm or more (and classic
spots), from the same count matrix as the gene files. One value per **row**, same length and
order as `positions.bin`, no header:

| File | Type | Content |
|---|---|---|
| `total.bin` | float32 | UMI over all genes of the layer (equal to the `hvg/` chunks' library size) |
| `genes.bin` | uint16 (uint32 if the layer has ≥ 65,536 genes; see `dtype`) | genes with a count > 0 |
| `pct_mt.bin` | float32 | 100 × UMI of the mitochondrial genes / total; 0 for a bin with no UMI. Only written when the layer has genes whose name starts with `MT-` (ignoring case, so mouse `mt-` too) |

The manifest layer names them:

```json
"binQc": {"dir": "binqc", "nRows": 144509,
          "metrics": [
            {"key": "total", "label": "Total counts", "file": "total.bin", "dtype": "float32", "unit": "UMI",
             "scale": "log", "min": 0.0, "max": 740.0, "median": 8.0, "q1": 1.0, "q3": 35.0},
            {"key": "genes", "label": "Genes detected", "file": "genes.bin", "dtype": "uint16", ...},
            {"key": "pctMt", "label": "% mitochondrial", "file": "pct_mt.bin", "dtype": "float32", "scale": "linear", ...}],
          "mtPattern": "^MT- (ignoring case)", "mtGenes": ["MT-ND1", "MT-ND2", ...]}
```

`scale` is the histogram axis the viewer uses (`log` = log10(1 + x)); `min` … `q3` are over all
rows. The viewer's **Bin QC filter** panel (`qc.js`) draws a histogram per metric with a min and
a max threshold. It is off until a threshold is set; a bin outside any threshold is greyed on
the slide and left out of Circle → UMAP (PCA, UMAP, Leiden, markers). The thresholds are stored
in Circle → UMAP's saved JSON as

```json
"qcFilter": {"layer": "008um", "thresholds": {"total": {"min": 13, "max": null}, "pctMt": {"min": null, "max": 25}},
             "nRows": 144509, "kept": 61318, "removed": 83191, "rule": "..."}
```

(`null` = no bound; `qcFilter: null` = no filter) and in the page URL as
`?qc=<layer>:<metric>=<min>-<max>,...`, e.g. `?qc=008um:total=20-,pctMt=-25`. Exports made
before 11 November 2026 have no `binQc` block; the panel is then hidden and nothing is filtered.

### `hvg/`: data for the circle → UMAP tool

Written for layers with bins of 8 µm or more (and classic spots), when at least 50 bins
have ≥ 10 UMI.

- `genes.json`: `["Gene", ...]`, the highly variable genes; their order is the column
  number used in the chunks. Up to 2,000, chosen once over the whole section
  (log1p CP10k, dispersion z-scored within 20 mean bins, on bins with ≥ 10 UMI).
- `index.json`:

  ```json
  {"chunkPx": 1827, "chunkUm": 500.0, "nGenes": 2000, "minUmi": 10,
   "selection": "seurat dispersion, log1p CP10k, whole section",
   "chunks": [{"tx": 0, "ty": 3, "n": 4021, "bytes": 512345}, ...]}
  ```

  `chunkPx` = round(500 µm / micronsPerPixel), or 1024 when the pixel size is unknown.
- `<tx>_<ty>.bin`: all bins (every bin with a position, low-UMI ones included) whose
  position satisfies `floor(x / chunkPx) = tx` and `floor(y / chunkPx) = ty`, computed on
  the float32 values of `positions.bin`. `tx`/`ty` can be negative for bins just off
  the image edge (the file name then has a minus sign). Layout, with `n` bins and `nnz`
  stored values:

  | Offset (bytes) | Type | Count | Content |
  |---|---|---|---|
  | 0 | uint32 | 1 | `n` |
  | 4 | uint32 | n | row in `positions.bin` |
  | 4 + 4n | float32 | n | library size: UMI over **all** genes, so normalisation is honest |
  | 4 + 8n | uint32 | n + 1 | CSR row pointer into the next two arrays; `nnz` = last value |
  | 8 + 12n | uint16 | nnz | gene column (index into `hvg/genes.json`) |
  | 8 + 12n + 2nnz | uint16 | nnz | UMI count, clipped at 65535 |

  File size = 8 + 12n + 4nnz.

### Cells as a layer (`kind: "cells"`, Xenium)

A cell-segmented platform's cells are written in exactly the bin format above, one "bin" per
cell (`expr.build_cell_layer`, called by `xenium.py`; a CosMx or MERSCOPE adapter can call it
with centroids in µm, the pixel size, a cells x genes sparse matrix and the gene names):

- `positions.bin`: the cell centroid (`x_centroid`, `y_centroid` of `cells.csv.gz` /
  `cells.parquet`, in µm) divided by the pixel size (`pixel_size` in `experiment.xenium`), in
  the order of the count matrix's barcodes. Xenium's µm frame has its origin at the top-left
  corner of the morphology image, so no other transform is applied.
- `genes.json`, `gene/<NAME>.bin`: features of type "Gene Expression" only; control probes
  and codewords (`NegControlProbe`, `NegControlCodeword`, `BLANK`, `UnassignedCodeword`,
  `DeprecatedCodeword`, `Intergenic`, `antisense`, genomic controls) are left out. Values are
  transcript counts from `cell_feature_matrix` (quality ≥ 20 molecules assigned to the cell).
- `binqc/`, `hvg/`, `umap/` and `region.bin` as for bins. Targeted panels rarely contain `MT-`
  genes, so `pct_mt.bin` is then not written.

The viewer words such a layer as cells ("cells", "transcripts") instead of bins / UMI.

## Molecules: `molecules/`

Written by `prepare.py xenium` (`molecules.py`, platform-neutral; `xenium.py` reads
`transcripts.parquet` when present, else `transcripts.zarr.zip` level 0) and copied unchanged
by `export`. One file per gene, fetched only when the gene is picked in the Molecules panel.

- Kept: transcripts with quality value (`qv` / `quality_score`) ≥ 20 (`qvMin`; Xenium's own
  cut for its cell-feature matrix), of genes only: the control codewords listed above (and
  any the zarr's `gene_category` does not mark `is_gene`) are left out. Molecules outside the
  image are dropped. All of them are counted in `index.json` `excluded` / `lowQv`.
- `gene/<NAME>.bin`: **uint16 pairs `qx, qy`**, one per molecule, little-endian, no header.
  Number of molecules = file size / 4. Full-resolution pixel `x = qx * scale[0]`,
  `y = qy * scale[1]`, with `scale = [width / 65535, height / 65535]` of the image (`qx =
  round(x / scale[0])`), so a step is ~0.8 px (0.18 µm) on a 54,000 px wide Xenium image;
  float32 x, y would be twice the size. Molecules are in a **seeded random order**: any prefix
  is a uniform sample of the gene, which the viewer uses to draw fewer dots when zoomed out.
  Genes with no molecules have no file.
- `index.json`:

  ```json
  {"format": "xy-uint16", "n": 166624389, "nGenes": 480,
   "bounds": [0, 0, 54002, 40815], "scale": [0.824, 0.623],
   "units": "full-resolution image pixels", "pixelSizeUm": 0.2125,
   "order": "seeded random: any prefix is a uniform sample of the gene",
   "excluded": {"controls": 76991, "outsideImage": 0, "unknownGeneId": 0},
   "platform": "xenium", "source": "transcripts.zarr", "qvMin": 20.0, "read": 208946850, "lowQv": 42245470,
   "genes": [{"n": "Gad1", "c": 1234567}, ...]}
  ```

  (numbers from a 480-gene mouse brain section). `genes` is in the source's gene order; `c` = molecules in the file.

The viewer draws a picked gene's molecules as small dots in one colour per gene, up to 12 genes
at once. With the whole section in view it draws only the first *k* of each file (*k* a power of
two; over all genes shown at most about one dot per 4 screen pixels, at least 200,000 and at most
6 million, 1.5 million on a phone; `config.js` `moleculeBudget` / `moleculesPerPixel` change
this), and the whole file once zoomed in; the panel says when it is showing a sample.

### `mat/`: the same, for marker squares (protein grid layers)

A protein grid layer (`protgrid.py`, Akoya and other multi-channel images) has no counts,
so instead of `hvg/` it carries every embedded channel square by square, cut into the
same 500 µm chunks:

- `index.json`:

  ```json
  {"chunkPx": 987, "chunkUm": 500.0, "scale": 100, "nChannels": 26,
   "channels": ["CD104", "CD45", ...], "normalize": "bg-asinh-v1",
   "bg": [...], "hi": [...], "cofactor": 0.2, "median": [...],
   "nSquares": 1084213, "mean": [...], "sd": [...], "rawMean": [...],
   "chunks": [{"tx": 0, "ty": 3, "n": 970, "bytes": 54324}, ...]}
  ```

  `channels` = the layer's embedded channels (manifest `channels`), in column order.
  `bg` / `hi` = each channel's whole-slide background (median of the off-tissue squares) and
  99.5th percentile above it over the tissue (leaving out the squares flagged everything
  bright, so a fold cannot set it), so the circle tool scales a selection exactly
  as the section did: x = asinh(max(v − bg, 0) / hi / cofactor) / asinh(1 / cofactor) for
  markers and gating, and for the layout the centred log-ratio of u = max(v − bg, 0) / hi,
  log(u + 0.1) minus the square's mean (protgrid.clr_features); squares with no channel at
  x ≥ 0.2 are left out (low signal).
  `median` = each channel's median x over the section (the `pct` of the cluster markers).
  `mean` / `sd` = each channel's mean and (population) standard deviation of x, and `rawMean`
  its mean raw value v, over all `nSquares` squares of the slide, from the stored (rounded)
  values. The circle tool compares each shape with the rest of the slide from these, and
  reads its clusters' markers and suggested names against the slide. Exports made before
  10 November 2026 lack them; the viewer then works them out once from every chunk (a
  computer only; a phone skips it).
  Exports made before November 2026 have no `normalize`, and `p1` / `p99_5` instead of
  `bg` / `hi` (x′ = clip((v − p1) / (p99.5 − p1), 0, 1)); the viewer reads both.
- `<tx>_<ty>.bin`: every square of the chunk (chunk of a square as for `hvg/`):

  | Offset (bytes) | Type | Count | Content |
  |---|---|---|---|
  | 0 | uint32 | 1 | `n` |
  | 4 | uint32 | n | row in `positions.bin` |
  | 4 + 4n | uint16 | n × nChannels | square-major mean intensities, `round(v × scale)`, clipped at 65535 |

  File size = 4 + 4n + 2n·nChannels.

### `qc.json` and `qc_flags.bin`: slide QC (protein grid layers)

Written by `protgrid.py` next to `positions.bin`; the manifest layer's `qc` block names
them and repeats the headline numbers (`everythingBright`, `saturatedSquares`, `weak`,
`saturated`, `seamsFlagged`).

- `qc.json`:

  ```json
  {"channels": [{"name": "CD90", "inUmap": true, "status": "saturated", "background": 1.09,
                 "p99_5AboveBg": 253.7, "signalToBackground": 232.0, "saturatedPct": 2.84,
                 "weak": false, "saturated": true}, ...],
   "dropped": [{"name": "open-cy5", "reason": "open / blank channel: not written"}, ...],
   "squares": 677828,
   "everythingBright": {"n": 757, "share": 0.00112, "rule": "..."},
   "saturatedSquares": {"n": 1203, "share": 0.00177, "rule": "..."},
   "seams": {"columns": {"channelsChecked": 27, "channelsWithSteps": 9, "spacingPx": 960, "spacingUm": 486,
                         "channels": ["CD45", "CD90", ...], "repeats": 24, "flagged": true},
             "rows": {...}, "flagged": true, "rule": "..."},
   "rules": {"weak": "...", "saturated": "...", "background": "...", "limit": [255]},
   "note": "UMAPs, clusters and names are per slide: ...",
   "flagsFile": "qc_flags.bin", "flagBits": {"everythingBright": 1, "saturated": 2}}
  ```

  `status` is one of `ok`, `weak`, `saturated`, `flat` (not in the UMAP), `nuclear`
  (tissue mask only). Values are in the 8-bit units of the gene files. `seams`: per axis, the
  channels whose full-resolution column (row) means step sharply at one repeated spacing - a
  camera tile edge - and that spacing; `flagged` when at least two channels share it. `grid`
  ({x0, y0, w, h} in px) when the edges sit in the same place in both halves of the slide;
  `correction` when the square means were corrected for the seams: `corrected` / `notCorrectable`
  markers, the median excess step at tile edges before / after (percentage points), and
  `file` = `qc_seams.json`.
- `qc_seams.json` (only when a correction was made): `grid`, `tiles` ({rows, cols, row0, col0}:
  tile (r, c) of the grid is entry (r - row0) * cols + (c - col0)), per marker with steps a record
  (`method` = `gain`, `gain+offset` or `none`, `excessBefore` / `excessAfter`, `controlBefore` /
  `controlAfter`, `borderLine` (the factors of the pixels next to each edge), `shading` (per
  orientation the share of edges whose step goes up, and the median step), the reason when not
  corrected), and `factors`: per corrected marker `g` and `c`
  per tile (corrected = raw / g - c, in raw image units) and the tiles that were `fitted`. The
  gene files and `mat/` hold the corrected means; the image tiles are not changed.
- `qc_flags.bin`: one uint8 per row of `positions.bin`. Bit 1: everything bright (more
  than 80% of the UMAP channels above their 95th percentile over the tissue: folds,
  bubbles, bright edges). Bit 2: more than half of the square's pixels at the detector's
  limit in some written channel.

### `umap/`: whole-section UMAP and Leiden clusters

Written under the same conditions as `hvg/`, when the clustering packages are installed
(see [PREPARING_DATA.md](PREPARING_DATA.md#clustering-and-the-whole-section-umap)).

- `coords.bin`: float32 pairs `x, y` (UMAP coordinates) per row; `NaN, NaN` for bins
  that were not embedded: RNA bins with < 10 UMI or with fewer than 3 of the variable genes
  detected; protein squares at background in every embedded marker. Such bins have no
  profile to place (hundreds of them are one identical point).
- `clusters_r<res>.bin`, one per resolution (`clusters_r0.3.bin`, `clusters_r0.6.bin`,
  `clusters_r1.0.bin`): uint16 cluster per row, `0` = the largest cluster,
  `65535` (0xFFFF) = not embedded.
- `index.json`:

  | Field | Meaning |
  |---|---|
  | `nBins` | rows in `positions.bin` |
  | `nEmbedded` | bins that were embedded (RNA: ≥ `minUmi` UMI and ≥ `minHvgDetected` variable genes) |
  | `minUmi` | 10 (RNA) |
  | `minHvgDetected` | 3 (RNA; absent in older exports, where every bin ≥ `minUmi` was embedded) |
  | `bounds` | `[xmin, ymin, xmax, ymax]` of the embedding (never NaN: the file is strict JSON) |
  | `params` | every parameter of the pipeline (HVGs, PCs, neighbours, min_dist, seed, normalisation, clustering, markers) |
  | `versions` | versions of umap-learn, leidenalg, igraph and numpy used |
  | `clusterings` | per resolution: `resolution`, `file`, `n` clusters, `sizes` (bins per cluster), `markers` |
  | `default` | index into `clusterings` shown first (the 0.6 resolution) |

  `markers[k]` lists up to 10 genes for cluster `k` as
  `{"gene", "lfc", "pct"}`: log2 fold change of mean CP10k expression (cluster vs rest)
  and the fraction of the cluster's bins expressing the gene, among genes expressed in
  ≥ 25 % of the cluster. Only genes with `lfc > 0` are listed.

## Single-cell / single-nucleus RNA-seq datasets

`prepare.py rnaseq` writes a different, smaller layout straight into the export folder:

```
<id>/manifest.json
<id>/umap.bin
<id>/obs/<NN>.bin
<id>/genes.json
<id>/gene/<NAME>.bin
```

`manifest.json`:

```json
{"id": "...", "label": "...", "kind": "rnaseq", "static": true,
 "nCells": 52000, "nGenes": 24000, "valueLabel": "log1p(counts per 10k)",
 "bounds": [xmin, ymin, xmax, ymax],
 "obs": [
   {"name": "leiden", "type": "category", "file": "03.bin",
    "categories": ["0", "1", ...], "counts": [9000, 7000, ...], "colors": ["#1f77b4", ...]},
   {"name": "total_counts", "type": "numeric", "file": "05.bin", "min": 800.0, "max": 25000.0}
 ]}
```

| File | Content |
|---|---|
| `umap.bin` | float32 pairs `x, y` per cell. Cell order = the input's order (in `--mtx` mode, cells without a UMAP row are dropped first). |
| `obs/<NN>.bin`, `type: "category"` | uint16 per cell: index into `categories`. Columns with fewer than 2 or more than 200 distinct values are not exported. `colors` is present only if the input carried a palette (`uns/<column>_colors`, streaming mode). |
| `obs/<NN>.bin`, `type: "numeric"` | float32 per cell. `min`/`max` are the 1st and 99th percentiles. Integer columns with ≤ 30 distinct values are treated as categories. |
| `genes.json` | `[{"n": name, "p": fraction of cells expressing, "m": maximum value}]`, sorted by `p`, most widely expressed first. |
| `gene/<NAME>.bin` | uint32 `k`, then uint32[k] cell indices, then float32[k] values, for the cells where the gene is non-zero. |

`valueLabel` is `log1p(counts per 10k)` when the input held integer counts (the export
normalises them), otherwise `expression (as provided)`. Genes expressed in no cell are
skipped; in streaming mode (`.h5ad` files over 2 GB) genes in fewer than 3 cells are
skipped too.

Column names, category values, the label and the id pass through de-identification
before they are written (see [DEIDENTIFICATION.md](DEIDENTIFICATION.md)). Cell barcodes
are never written.

## `about.html` (optional)

An HTML fragment at the root of the data (next to the dataset folders, e.g. the top of the
blob container). The page shows it in the "Read more about this project" box, fetched
with the same access token as the data. It is inserted into the page as HTML, so only put
your own trusted content there.

## The prepare cache (not for upload)

For reference: `prepare.py` writes `<cache>/<id>/manifest.json`, `image.zarr/` (an
OME-NGFF 0.4 multiscale image with `omero` channel metadata, readable by napari and
other OME-Zarr tools), `expr/<key>.zarr/` (a gene-major sparse count matrix plus bin
positions) and `overlays/`. `server.py` and `prepare.py export` read it. It holds the
full-resolution raw image and local file paths, so it should stay on your machine.
