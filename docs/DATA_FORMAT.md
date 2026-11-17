# Data format

This is the layout that `prepare.py export` (images, spatial expression) and
`prepare.py rnaseq` (single-cell UMAPs) write, and that the viewer reads. It is derived
from `export.py`, `hvgchunks.py`, `binqc.py`, `sectionumap.py`, `regions.py`, `rnaseq.py`,
`expr.py`, `molecules.py`, `cellshapes.py`, `xenium.py`, `cosmx.py`, `merscope.py`, `stdformats.py`, `viewer.js`, `outlines.js` and `qc.js`. You only need it if you want to write the files with your own tools or
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
  and `umap/` all refer to bins by row. No file contains a barcode, except the optional
  `ids.txt.gz` written only when asked for (`prepare.py export --ids`, see below). It is written
  only when every id fully matches its platform's instrument id format (an allow-list) and none
  looks like a sample ID; otherwise the export stops, writes no list, and says which layer and how
  many ids failed without printing any of them.
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
  expr/<key>/genepack/<hash>-<i>.pack    every gene's points, packed (exports from 16 Nov 2026; see "Gene packs")
  expr/<key>/gene/<NAME>.bin             the same, one file per gene (older exports, or --gene-files loose / both)
  expr/<key>/region.bin                  only if an overlay has polygons
  expr/<key>/binqc/{total.bin, genes.bin, pct_mt.bin}   per-bin QC values (RNA layers, exports from 11 Nov 2026)
  expr/<key>/hvg/{genes.json, index.json, <tx>_<ty>.bin}
  expr/<key>/mat/{index.json, <tx>_<ty>.bin}   protein grid layers only (instead of hvg/)
  expr/<key>/umap/{index.json, coords.bin, clusters_r<res>.bin}
  expr/<key>/shapes/<set>/{index.json, <tx>_<ty>.bin}   cells layers with segmentation boundaries: cell (and nucleus) outlines
  molecules/index.json                   imaging platforms (Xenium): every transcript, per gene
  molecules/genepack/<hash>-<i>.pack     packed like expr/ (molecules/gene/<NAME>.bin in older or loose exports)
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
| `kind` | `akoya`, `tiff`, `xenium`, `cosmx`, `merscope`, `stereoseq`, `geomx`, `visium`, `anndata`, `omezarr` or `spatialdata` (the `prepare.py` subcommand used) |
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
| `ids` | only with `prepare.py export --ids` or `prepare.py ids`: `{"file": "ids.txt.gz", "n": rows, "bytes": file size, "column": "barcode" \| "cell_id"}`, see `ids.txt.gz` below |
| `shapes` | cells layers only, if the platform's cell outlines were imported: `{"version": 1, "format": "poly-u16-tiles", "dir": "shapes", "chunkPx", "sets": [{"key", "label", "nPolygons", "nRowsWithShape", "nVertices", "chunks", "bytes", "source"}]}`, see [Cell outlines](#cell-outlines-shapesset-cells-layers) |
| `umap` | only if `umap/` was written: `{"nEmbedded", "clusterings": [{"resolution", "name", "source", "n"}], "default"}`; with the lab's own results also `id` per clustering and `"embedding": "user"` (see umap/ below) |

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
Genes with no counts are listed but have no points (no `gene/` file, no place in the pack).
In a packed layer each gene with points also has `s`, `o`, `b`: its shard, byte offset and
length in the layer's pack (`{"n": "Gapdh", "t": 123456.0, "s": 0, "o": 4128, "b": 2436}`; see
"Gene packs" below).

### `positions.bin`

float32 pairs `x, y`, one per bin that has a position, in count-matrix barcode order.
Number of rows = file size / 8. This order defines **row** numbers.

### `ids.txt.gz` (optional): barcode / cell id per row

Only when the export was made with `prepare.py export --ids`, or `prepare.py ids --id <id> --out <cache>
--dest <export root>` was run on it later (it reads the barcodes from the prepare cache and adds the files and
the manifest entries; nothing else changes). gzip of UTF-8 text, one id per line, one line per row of
`positions.bin` and in its order, no header. `column` in the manifest entry says what they are: `barcode`
(Visium spots, Visium HD and Stereo-seq bins) or `cell_id` (cells layers: Xenium `cell_id`, CosMx
`<cell_ID>_<fov>`, MERSCOPE `EntityID`, Stereo-seq cell id). GeoMx ROI layers and protein grid layers get none.
The viewer reads it for "Save barcodes (CSV)" in Combine and export shapes; without it that button says how
to add one. A host that unpacks `.gz` files itself is fine: the viewer reads the text either way.

Before anything is written, every layer's ids pass two gates (layerids.py), and one refused layer means no
list for any layer of the dataset:

1. **Allow-list.** Each id must match, in full, the instrument id format of its platform, as the readers put it
   in the cache:

   | Platform | Accepted ids |
   |---|---|
   | Visium HD bins | `s_<3 digits>um_<5 digits>_<5 digits>-1` |
   | Visium spots | 16 of `ACGT`, `-` and a 1-3 digit GEM well (`-1`; aggr gives `-2` ...) |
   | Xenium cells | 8 letters `a`-`p`, `-`, 1-3 digits (`aaabbbcc-1`), or an integer (Xenium before Ranger 1.3) |
   | CosMx cells | `c_<slide>_<fov>_<cell>` (AtoMx; all integers) or `<cell_ID>_<fov>` (made from the matrix) |
   | MERSCOPE cells | an integer `EntityID` |
   | Stereo-seq | bins `<x>_<y>` or an integer (`--barcodes stereopy`); cells an integer |

   A CosMx id keeps its slide number: it is the instrument's slide index, an integer, and the form AtoMx and
   Seurat use. Anything else, such as a merge prefix (`S1_`), a sample or run name glued to a barcode, or a
   GeoMx Sample_ID, is refused, never stripped or rewritten. Other platforms get no list.
2. **Sample-ID patterns.** `deid.find_ids` over the whole list, as for every other file.

A refusal prints the layer, the platform, the format and the count of failing ids, never an id. `deid.py
check` also reads `ids.txt.gz` (any `.gz` that unpacks to text) and reports a hit there as file and line
only.

Size, measured on the smoke-test demo data: 8 µm bins 13,291 rows 35.1 kB (2.64 bytes per row), 16 µm
bins 3,388 rows 8.5 kB (2.50 B/row), Xenium cells 838 rows 2.3 kB (2.69 B/row). As plain text a Visium HD
barcode is 22 bytes per row and a Xenium cell id 11, so gzip takes off 85-90 %; for comparison
`positions.bin` is 8 bytes per row. Grid barcodes compress this well because neighbouring rows differ in a
digit or two; ids with no order (random cell ids) cost more, up to their plain-text size divided by 2-3.

### `gene/<NAME>.bin` (and each gene's bytes in the pack)

float32 triples `x, y, value` for each bin where the gene's count is > 0 (and the bin has a
position), ordered by bin. `value` is the raw UMI count. Number of points = file size / 12.
A packed layer holds exactly these bytes per gene, end to end (next section).

### Gene packs: `genepack/<hash>-<i>.pack`

Since 16 November 2026 an export writes each group of per-gene files (an expression layer,
the molecules, an RNA-seq set) as a few **packs** instead of one file per gene: a Visium HD
layer of ~16,000 genes becomes 1-23 files of about 64 MB instead of ~16,000 small ones
(written by `genepack.py`; `prepare.py export --gene-files loose` writes the old layout,
`--gene-files both` the two side by side).

- A shard holds the payloads of the group's genes **end to end, byte for byte what
  `gene/<NAME>.bin` would hold**, in listing order. A gene is never split across shards; a new
  shard starts when the next gene would take it past the size (`--shard-mb`, default 64), so a
  gene larger than that has a shard of its own.
- `<hash>` (16 hex digits) is a hash over every name and payload of the group: shards of one
  export run share it, and a listing can never point into shards of another run.
- Each listing entry with points carries `s` (shard), `o` (byte offset) and `b` (bytes):
  `genes.json` of a layer or an RNA-seq set, the `genes` of `molecules/index.json`.
- The shard list is in the manifest, on the expression layer, on the `molecules` block, or at
  the top of an RNA-seq manifest:

```json
"genePack": {"version": 1, "dir": "genepack", "unit": 12, "genes": 16203, "bytes": 1458201936,
             "shards": [{"file": "3f9c0a1b2d4e5f60-0.pack", "bytes": 67108860}, ...]}
```

`unit` is the record size (12 for expression, 4 for molecules; absent for RNA-seq, whose
records have a header). The viewer reads one gene with **one HTTP Range request**,
`Range: bytes=o-(o+b-1)`, into `<dir>/<file>` of shard `s`; it requires a `206` answer of `b`
bytes (and, when the host exposes `Content-Range`, the range asked for out of a shard of the
manifest's size). A host that ignores Range answers `200` with the whole shard: the viewer then
keeps the shard and cuts genes out of it if it is at most `packWholeMaxMB` (config.js, default
160 MB), and otherwise says the host does not support Range requests. Export with
`--gene-files loose` for such a host.

`python genepack.py check <export>/<id>` checks a dataset's packs against its listings (every
byte accounted for); `python genepack.py repack <export>/<id> [--delete-loose]` packs the gene
files of an export made before the packs in place (no cache needed), verifies every gene
against its file, and only then (with `--delete-loose`) deletes the files.

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
rows. The viewer's **Quality filter** section (`qc.js`, Clusters tab) draws a histogram per metric with a min and
a max threshold. It is off until a threshold is set; a bin outside any threshold is greyed on
the slide and left out of Draw a region (PCA, UMAP, Leiden, markers). The thresholds are stored
in the JSON that Save JSON writes as

```json
"qcFilter": {"layer": "008um", "thresholds": {"total": {"min": 13, "max": null}, "pctMt": {"min": null, "max": 25}},
             "nRows": 144509, "kept": 61318, "removed": 83191, "rule": "..."}
```

(`null` = no bound; `qcFilter: null` = no filter) and in the page URL as
`?qc=<layer>:<metric>=<min>-<max>,...`, e.g. `?qc=008um:total=20-,pctMt=-25`. Exports made
before 11 November 2026 have no `binQc` block; the panel is then hidden and nothing is filtered.

### The saved view: Save JSON's `view` block and the link after `#`

Not part of an export: this is what the page itself writes. **Save JSON** and **Copy link to
this view** come from one function (`sessionState` in `viewer.js`); the link is the file's
state without the UMAP run's parameters and result. Both carry a `view` block (made and read
by `viewlink.js`):

```json
"view": {"v": 1, "mode": "sbs", "cam": [5120.5, 3890.25, -1.25], "scam": [5000, 3900, -2], "win": [1366, 722],
         "tab": "clusters", "img": {"on": true, "op": 0.8},
         "ch": [{"n": "DAPI", "on": true, "c": "#0000ff", "w": [10, 180]}],
         "genes": [{"n": "MARKER_FASC1", "c": "#ff3355", "r": [0, 12]}, {"n": "MARKER_EPI2", "c": "#33ddff", "r": [1, 7.5], "on": false}],
         "blend": "heatmap", "palette": "viridis", "gop": 0.9, "ps": 1.5, "tb": true,
         "mol": {"genes": [{"n": "GENE_A", "c": "#ffffff"}], "size": 1.5, "op": 0.8},
         "ov": ["regions"], "colorBy": "cluster", "pick": 3, "qcOutline": true}
```

- `mode`: `tissue`, `sbs` (side by side), `umap`, or `rna` (a single-cell dataset, whose settings
  are in a `rna` block: `colorBy`, `genes`, `view`, `palette`, `bin`, `ps`, `labels`, `cam`).
- `cam` = centre x, y in image pixels and deck.gl zoom (log2 of screen px per image px) of the
  tissue; `ucam` / `scam` the same for the full UMAP view and the UMAP half of side by side;
  `win` = the map's size in CSS px when saved. On opening, the zoom is corrected by
  log2(min(w'/w, h'/h)) so the area that was in view is in view again in a window of another size.
- `ch` lists the stains that are on or whose colour (`c`) or display range (`w`, 0–255) was
  changed; every other stain is off with its default. `genes[].r` is the range drawn, `on:
  false` a gene that is listed but switched off. `ov` = the keys of the outlines switched on.
- The layer, drawn shapes, the clustering on show with cluster names, cell sets and hidden
  clusters (`sectionClusters`) and the quality thresholds (`qcFilter`) are where Save JSON always had
  them, one level up. A file without `view` (saved before 16 November 2026) loads as it did.
- A shape made by Combine / Overlap / Subtract (selections.js) is a list of bins (`rows`, as
  before) and, when it was made from outlines only (Circle, Freehand, imported polygons), also
  carries its recipe `def`: `{"points": [[x, y], ...]}`, `{"geo": [[ring, hole, ...], ...]}` or
  `{"op": "union" | "inter" | "sub", "parts": [def, ...]}`, in image pixels. Save JSON writes
  both; the link writes only `def`, and the bins are found again from it on opening. A shape with
  a plain selection inside it has no `def` and stays a list of bins.

The link is `<page>?ds=<dataset>[&data=…|&static=…][&qc=…]#v=<payload>`. The payload is
`1d.` + base64url(deflate-raw(UTF-8 JSON)), or `1j.` + base64url(UTF-8 JSON) from a browser
without `CompressionStream`; the leading digit is the format version. The JSON is the state
above with defaults left out (`seed` 42, `mode` "bins", `resolution` "auto", no shapes). The
page aims at 2,000 characters for the whole link: outlines are first thinned to 40 points
(the outlines in a combined shape's `def` too, and shapes made from a selection dropped), then
shapes are left out, then cluster names;
what was left out is said under the button. A reader accepts payloads up to 16,000 characters
that unpack to at most 512 kB, keeps only known keys of the right type and range, looks names
up in the open dataset's own lists (it never fetches by a name from a link) and reads no
address from the payload; a payload it cannot read leaves the default view and a message.

### What the heatmap, dot plot and violin plots read

The **Heatmap and violin** section (heatmap, dot plot and violins; `plots.js`, numbers from `plots_worker.js`) computes
everything in the browser from files already described here; nothing extra is exported.

| Needed for | File |
|---|---|
| each gene shown | the gene's bytes: one Range request into the pack, or `gene/<NAME>.bin` in older exports (one request per gene, whole layer; placed on rows by walking `positions.bin` in order, or by exact position if the file is not in row order) |
| counts per 10k | `binqc/total.bin`. Without it (exports before 11 November 2026): the library sizes in the `hvg/` chunks under the groups; without those: log1p of the raw counts, labelled "not depth-normalised" |
| groups | `umap/index.json` and `umap/<file>` of the chosen clustering (labelled by its `name` if it has one, else "Leiden <resolution>"); `region.bin` and the layer's `regions[]`; drawn shapes and the clusters of a Draw a region run come from the page itself |
| default genes | `clusterings[].markers` (the first N per cluster) |
| single-cell datasets | `obs/<NN>.bin` of a `category` column, each gene's bytes (pack or `gene/<NAME>.bin`; values used as stored) |

Definitions: a bin counts for a group if it is in it, kept by the QC filter (when one is set)
and has ≥ 10 UMI (only when totals are available). Value of a bin = log1p(count × 10,000 /
total). **Mean** = over all kept bins of the group, zeros included; **fraction expressing** =
bins with a count > 0 / kept bins; per-gene scaling (optional) = z-score across the groups
shown (sd with n − 1, colour clipped at ±2.5) or (mean − min) / (max − min). **Violin** =
Gaussian kernel density with Scott's bandwidth (sd × n^(−1/5)), drawn between the smallest and
largest value only, quartiles by linear interpolation; "expressing only" leaves out the bins
with a count of 0. The saved CSV holds the numbers plotted, with the number of bins per group
and how many were left out by the QC filter and the UMI minimum.

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

### Cell outlines: `shapes/<set>/` (cells layers)

Written by `cellshapes.py` for a cells layer whose platform came with segmentation boundaries
(or for which a lab gave its own, `prepare.py outlines`): the outline of every cell, and of
every nucleus when the platform has them, as polygons in spatial tiles. The viewer draws them
instead of the dots once the cells are large enough on screen (`outlines.js`) and fetches only
the tiles under the view. Format name `poly-u16-tiles`, version 1. A layer has them when its
`expr[]` entry has

```json
"shapes": {"version": 1, "format": "poly-u16-tiles", "dir": "shapes", "chunkPx": 2353,
           "sets": [{"key": "cell", "label": "Cell boundaries", "nPolygons": 155754, "nRowsWithShape": 155754,
                     "nVertices": 3505410, "chunks": 296, "bytes": 14960900, "source": "cell_boundaries.csv.gz"},
                    {"key": "nucleus", "label": "Nucleus boundaries", "nPolygons": 156299, "nRowsWithShape": 155560,
                     "nVertices": 3215919, "chunks": 294, "bytes": 13806174, "source": "nucleus_boundaries.csv.gz"}]}
```

(numbers from a 155,754-cell mouse brain section: 96 bytes per cell for 22.5 vertices). The
first set is the one drawn for the cells; a set with key `nucleus` is the optional nucleus
outlines. An export without `shapes` is drawn as dots, as before.

- **Tiles** are the grid `hvg/` uses: `chunkPx` = round(500 µm / µm per pixel) (1024 when the
  pixel size is unknown), and a polygon is in the tile of **its cell's position**:
  `tx = floor(x / chunkPx)`, `ty = floor(y / chunkPx)` on the float32 values of
  `positions.bin`. A polygon is never cut or repeated at a tile edge; `overhangPx` in
  `index.json` is how far a vertex can reach outside its tile, so a viewer fetches the tiles of
  its view grown by that much.
- `shapes/<set>/<tx>_<ty>.bin`, with `n` polygons and `nv` vertices:

  | Offset (bytes) | Type | Count | Content |
  |---|---|---|---|
  | 0 | uint32 | 1 | `n` |
  | 4 | uint32 | 1 | `nv` |
  | 8 | float32 | 2 | `x0`, `y0`: the origin of this tile's vertices (image pixels) |
  | 16 | uint32 | n | row in `positions.bin` |
  | 16 + 4n | uint16 | n | vertices of each polygon; `0` = the cell has no outline |
  | 16 + 6n | uint16 | 2 nv | `qx, qy` per vertex: x = `x0 + qx * step`, y = `y0 + qy * step` |

  File size = 16 + 6n + 4nv. `step` (in `index.json`) is 1/16 pixel, or the next power of two
  under which the largest tile still fits in 16 bits. Entries are ordered by row.
- A **row** may appear more than once: a cell in several parts (a MultiPolygon, a cell with
  two nuclei) is several polygons with the same row. In the **first** set every row of
  `positions.bin` appears at least once; one with `count` 0 has no outline and keeps its dot.
  Further sets list only cells that have a polygon. No cell id is written.
- **Rings** are open (the first vertex is not repeated), **clockwise on screen** (x right, y
  down; the shoelace sum Σ(xᵢ yᵢ₊₁ − xᵢ₊₁ yᵢ) is positive), outer rings only (holes are left out),
  3 to 64 vertices.
- **What is changed on the way in** (all counted in `index.json` `stats`): repeated points and
  the closing vertex are dropped; vertices are removed where that moves the outline by no more
  than 0.1 µm (Douglas-Peucker; `--outline-tolerance-um`, 0 = keep all; 1.5 mask pixels for
  outlines traced from a label image), more for a polygon that would otherwise have over 64;
  parts under 2 % of their cell's area, polygons with fewer than 3 points or no area and
  polygons wider than two tiles are dropped; polygons of cells that are not in the count
  matrix, or have no position, are left out.
- `shapes/<set>/index.json`:

  ```json
  {"format": "poly-u16-tiles", "version": 1, "set": "cell", "label": "Cell boundaries",
   "units": "full-resolution image pixels", "pixelSizeUm": 0.2125, "chunkPx": 2353, "chunkUm": 500.0,
   "step": 0.0625, "overhangPx": 72.062,
   "rings": "open (first vertex not repeated), clockwise on screen, outer rings only",
   "nRows": 155754, "nRowsWithShape": 155754, "complete": true,
   "nPolygons": 155754, "nVertices": 3505410, "bytes": 14960900, "source": "cell_boundaries.csv.gz",
   "stats": {"polygonsRead": 155754, "polygonsNotInMatrix": 0, "polygonsOfCellsWithoutPosition": 0,
             "dropped": {"degenerate": 0, "tooLarge": 0, "smallParts": 0, "holes": 0},
             "rowsInSeveralParts": 0, "selfIntersecting": 10,
             "centreInside": 0.99976, "centreWithin1um": 0.99998, "centreWithinOneDiameter": 1.0,
             "verticesRead": 3890798, "verticesPerPolygon": [7, 23, 24], "meanVertices": 22.51,
             "areaUm2": [49.513, 74.825, 114.831],
             "simplify": {"toleranceUm": 0.1, "tolerancePx": 0.4706, "maxVertices": 64, "toleranceRaised": 0}},
   "chunks": [{"tx": 0, "ty": 5, "n": 80, "nv": 1798, "bytes": 7688}, ...]}
  ```

  `nRows` is the number of rows of `positions.bin` the set was made for (the export leaves a
  set out when it differs). `complete`: every row is listed (the first set). `centreInside`:
  the share of cells whose own position lies inside (one of) its own polygon(s), and
  `centreWithin1um` inside or within 1 µm; this is the importer's check that the outlines are
  in the cells' frame: under 90 % it warns, under 50 % (or when under half of the polygons
  belong to a cell of the matrix) it writes nothing. `selfIntersecting` polygons are kept as
  they are. `source` names the kind of file, never a path. Readers add what they know
  (`zPlane`, `zPlanes`, `otherPlane` for MERSCOPE; `traced` for label images).
- Where the polygons come from, and the transform to image pixels (the same as the centroids'):

  | Platform | Source | To image pixels |
  |---|---|---|
  | Xenium | `cell_boundaries` / `nucleus_boundaries` `.parquet`, `.csv.gz` (`cell_id, vertex_x, vertex_y, label_id`, µm); else `cells.zarr.zip` (`polygon_sets/<i>/{cell_index, num_vertices, vertices}`, or `polygon_vertices` of older runs) | µm ÷ pixel size |
  | MERSCOPE | `cell_boundaries.parquet` or a VPT run's `*_micron_space.parquet` (`EntityID, ZIndex, Geometry` as WKB, µm); else `cell_boundaries/feature_data_<fov>.hdf5` | the manifest `frame`; the z-plane of the imported mosaic (else the middle one), and for a cell with no polygon there its nearest plane |
  | CosMx | `*polygons.csv` (`fov, cellID, x_global_px, y_global_px`) when the export has one; else traced from `CellLabels/CellLabels_F<fov>.tif`, exactly along the pixel edges of every label, then simplified | the manifest `frame`; a traced FOV is placed as its image is in the mosaic |
  | Stereo-seq | `cellBin/cellBorder` of the cell bin GEF: up to 32 (x, y) offsets from the cell's x, y in DNB | (cell + offset + file offset + 0.5 − image offset) × image pixels per DNB. The border points are the outermost DNBs of the cell, so the outline runs through their centres |
  | any | your own file (`prepare.py outlines --boundaries`): GeoJSON polygons, a vertex table, or a parquet with WKB geometry | `--boundary-units` px, um or frame |

The viewer draws dots while the median cell (`spotDiameterPx`) is under 5 screen pixels across,
outlines alone from 9, both cross-faded in between, and stays on dots when more than 120,000
outlines (20,000 on a touch device) would be in view (`config.js` `outlineFromPx`,
`outlineFullPx`, `outlineBudget`). Colour never rebuilds a polygon: whatever colours or picks out
cells (genes, clusters, cell sets, Colour by a value, Draw a region's shapes, a region picked with a
click, the quality filter) is one RGBA per row that the outline layers read. A link to the view
carries the Outlines / Dots choice as the view block's `cs` entry (`{"m": "outlines" | "dots", "n": 1}`
with nucleus outlines on); a link without it leaves the viewer's own choice.

In the prepare cache the same files are in `<cache>/<id>/shapes/<layer key>/<set>/`.

### ROI segments as a layer (`kind: "rois"`, GeoMx DSP)

`prepare.py geomx` (`geomx.py`) writes a GeoMx DSP study as one layer with key `rois`, in the
bin format above with one "bin" per ROI segment (AOI). The dataset's `kind` is `geomx`.
What differs from a layer of bins or cells:

- `expr[]` entry: `"kind": "rois"`, `"platform": "geomx"`, `"unit": "counts"`, `nBins` =
  segments, `nRois` = ROIs, `binSizeUm` null, `spotDiameterPx` = 0.6 × the diameter of a circle of
  the median segment area, `"layout": "scan"` or `"grid"` (no scan image: positions are a
  drawn layout).
- `positions.bin`: the centre of the segment's mask from the scan's OME-XML (else of the ROI
  polygon), moved inside the ROI polygon when it falls outside, plus the scan's x offset when
  several scans share the canvas. No two rows share a position. Rows are in the scans'
  left-to-right order, then annotation order.
- `genes.json`, `gene/<NAME>.bin`: target values = DCC probe counts collapsed with the PKC as
  GeomxTools `aggregateCounts` (one probe: its count; several: the geometric mean, 0.5 added
  to all of them first when any is 0). Values of multi-probe targets are not whole numbers.
  Targets whose PKC `CodeClass` starts with `Negative` are not written.
- `overlays/rois.geojson`: one polygon per ROI (`properties.name`, `roi`, `slide`, `aois`),
  the first overlay and the layer's `regionOverlay`; `region.bin` links each segment to its ROI.
- `binqc/`: the usual `total.bin`, `genes.bin` (and `pct_mt.bin`) plus, each float32 per row
  and listed in the `binQc.metrics` block like the others: `raw_reads.bin` (`rawReads`),
  `pct_aligned.bin` (`pctAligned`, 100 × Aligned ÷ Raw), `saturation.bin` (`saturation`,
  100 × (1 − sum of the DCC counts ÷ Aligned)), `area.bin` (`area`, µm²), `nuclei.bin`
  (`nuclei`), `neg_geomean.bin` (`negGeoMean`), `loq.bin` (`loq`), `genes_above_loq.bin`
  (`genesAboveLoq`), `q3.bin` (`q3`), `q3_over_neg.bin` (`q3OverNeg`). A metric that is missing
  for some segment (no nuclei column) is left out. With several panels the background and LOQ
  metrics are those of the panel with the most targets.
- `hvg/`: the same chunk format, but `genes.json` lists **every target with a count** (no
  variable-gene selection, no minimum number of rows); `index.json` has `"selection"`,
  `"rounded"` (values that were not whole numbers and were rounded for the uint16 counts) and
  `"clipped"` (values above 65,535); the manifest's `hvg` entry has `"allTargets": true`.
- `umap/`: **no embedding.** `index.json` has `"embedding": {"source": "none", "why": ...}`,
  `"nEmbedded": 0`, and `coords.bin` is all NaN; the manifest's `umap.embedding` is `"none"`.
  `clusterings` holds the annotation columns, each exactly like a lab's own clustering
  (`id` `user<N>`, `name`, `"source": "user"`, `file` `clusters_user<N>.bin`, `categories`,
  `sizes`, `nAssigned`, `markers`), de-identified the same way; slide and scan names are
  "Slide 1", "Scan 1" ... unless prepared with `--keep-slide-names`. No `clusters_r*.bin`, no
  `celltypes.json`.
- manifest `geomx` block: `nAoi`, `nRoi`, `nSlides`, `layout`, `panels` (name, version,
  targets, probes, negativeProbes), `collapse` and `loq` (the rules, as text),
  `normalisation`, `dcc` (`used`, `inFolder`, `sharePkcMapped`), `scans` (order, `offsetPx`,
  `size`, ROIs, AOIs, `unplaced`, `maskAreaOverAnnotation`: the median ratio of the masks' areas to the annotation's, 1.0 when the file and pixel size are right), `ntc` (the No Template Control's total count per plate),
  `background` (median background and LOQ per panel).

The viewer words such a layer as ROIs and counts, fills each ROI's outline with the colour of
the genes on show (one disc per segment where an ROI has several; no density grid), hides the Tissue / UMAP switch and Draw a region's own UMAP run, and leaves the rest
(legend, groups, plots, volcano, QC filter) as for any layer (`roi.js`).

In the prepare cache, `geomx/rows.json` keeps the annotation's raw text and the QC values per
row for the export; like the rest of the cache it is never uploaded.

### CosMx and MERSCOPE datasets (`kind: "cosmx"`, `"merscope"`)

Written by `prepare.py cosmx` / `prepare.py merscope` (`cosmx.py`, `merscope.py`, with
`flatfiles.py`). They use the same files as a Xenium dataset: tiles, a `cells` layer and
`molecules/`. What differs:

- `manifest.json` has a `frame` entry, the map from the platform's own coordinates to
  full-resolution image pixels, so that other data in those coordinates (cell outlines, your
  own annotations) can be placed (`cellshapes.py` uses it for the cell outlines): `{"from": "<what the coordinates are>", "affine": [[a, b, c],
  [d, e, f]]}` with `x_px = a*x + b*y + c`, `y_px = d*x + e*y + f`.
  - CosMx: from `x_global_px, y_global_px` (instrument pixels). `a` and `e` are +1 or -1
    (y is turned over when global y points up) and `c`, `f` move the mosaic's corner to
    (0, 0); with a density background they are scaled by `instrumentPixelUm` / `pixelSizeUm`.
    `instrumentPixelUm` is the µm per instrument pixel used.
  - MERSCOPE: from `global_x, global_y` (µm): the first two rows of
    `micron_to_mosaic_pixel_transform.csv`. `zPlane` is the mosaic z-plane imported.
- `image.background` is present only when the dataset had no images: `"transcript density"`
  or `"cell density"`. The single channel ("Transcript density (no image)") is then the
  number of transcripts (cells) per pixel of `pixelSizeUm` (1 µm unless `--background-um`),
  Gaussian-smoothed (σ 1.5 µm; 5 µm for cells) and scaled so that the 99.9th percentile of
  the non-empty pixels is 255. No protein grid layer is built on such an image.
- CosMx image channels are either the raw morphology stains (uint16) or "CellComposite
  red / green / blue", the vendor's 8-bit rendering, stored as three channels.
- The `cells` layer has `platform: "cosmx"` or `"merscope"`. `positions.bin` is the centroid
  (`CenterX/Y_global_px`; `center_x`, `center_y`) put through `frame`, in count-matrix order.
  Left out of `genes.json` and `gene/`: CosMx `NegPrb*`, `Negative*`, `SystemControl*`,
  `FalseCode*` and the `cell_ID` 0 rows; MERSCOPE `Blank-*`. `cellDiameterUm` comes from the
  median `Area` (CosMx) or the median √(bounding-box width × height) (MERSCOPE).
- `molecules/index.json` has `"platform": "cosmx"` / `"merscope"`, `"source": "tx_file"` /
  `"detected_transcripts"` and `"read"` (rows of the table); there is no `qvMin` / `lowQv`
  because these tables carry no quality value. `excluded.controls` counts the control codes
  above. `genes` is sorted by name (the names found in the table).

### Stereo-seq layers (`kind: "stereoseq"`)

`prepare.py stereoseq` (`stereoseq.py`) writes ordinary layers in the formats above: one bin
layer per bin size, summed from the DNB counts of a bin GEF or a GEM, and a `cells` layer from
the cell bin GEF. The manifest `kind` is `stereoseq`; the viewer shows it like a Visium
dataset (the image under a transcripts panel). What is specific to it:

| Where | Field | Meaning |
|---|---|---|
| `expr[]` | `key`, `label` | `bin020`, `bin050`, `bin100` ... (`bin` + DNBs per side, three digits), labelled `bin20 (10 µm)` ...; `cells`, labelled `Cells (cell bin)` |
| `expr[]` | `binSizeUm` | DNBs per side × µm per DNB (10, 25, 50 ... at 0.5 µm per DNB); `null` for cells |
| `expr[]` | `binDnb` | bin layers: DNBs per side of a bin |
| `expr[]` | `platform` | `stereoseq`, on the bin layers and on cells |
| `expr[]` | `midCount` | the sum of the layer's counts (MID is Stereo-seq's name for UMI) |
| `expr[]` | `medianAreaUm2`, `fileOffsetDnb` | cells: the median cell area; the cell bin GEF's `offsetX`, `offsetY` |
| manifest | `stereoseq` | how the layers were made: `dnbUm` (µm per DNB), `imageScale` (image pixels per DNB), `offsetDnb` (`--offset`), `source` (`gef` or `gem`), `gefVersion` or `gemFormat`, `levelsStored` (the bin sizes in the file), `fileOffsetDnb` (the file's `minX`, `minY` or `#OffsetX`, `#OffsetY`), `recordsRead`, `midRead`, `recordsInTissue`, `midInTissue`, `tissueCut` (`GEF contour (n rings)`, `mask image` or `none (source kept whole)`), `nGenes`, and `recordsOffImage` when records were left out |

- `positions.bin` of a bin layer holds the bin centres, (first chip DNB of the bin + half a
  bin − image offset) × image pixels per DNB, where chip DNB = the x, y stored in the file +
  the file's offset. Only bins with at least one count under the tissue are rows, in row-major
  order of the bin grid. A cell's position is its centre + half a DNB.
- `genes.json` holds gene names (the gene id where a gene has no name). A name that several
  gene ids share is written `Name_<geneID>` for all but the first, so each id keeps its own file.
- In the cache, the `barcodes` of a bin layer are `x_y` of the bin's first DNB as the file
  stores it, and of the cells layer the cell `id` (with `--barcodes stereopy`, Stereopy's
  `x * 2^32 + y` for both). They are not exported.
- The cells layer's outlines (`shapes/cell/`) come from `cellBin/cellBorder`, see
  [Cell outlines](#cell-outlines-shapesset-cells-layers).

### Standard-format datasets (`kind: "anndata"`, `"omezarr"`, `"spatialdata"`)

Written by `prepare.py anndata` / `omezarr` / `spatialdata` (`stdformats.py`, reading with
`annstore.py` and `ngff.py`) from a spatial AnnData, an OME-Zarr image or a SpatialData store
(see [Standard formats](PREPARING_DATA.md#standard-formats-spatial-anndata-ome-zarr-spatialdata)).
The export is an ordinary one: tiles, `expr/<key>/` layers, and for SpatialData `molecules/`
and `shapes/<set>/`. What differs:

| Where | Field | Meaning |
|---|---|---|
| manifest | `imaging` | `true` when the image is a fluorescence image or a density background (not a photo): the side panel shows its channels, as for Xenium. Absent for an H&E / RGB photo |
| manifest | `image.background` | `"cell density"` when the file had no image: one grey channel (`#b4b4b4`) drawn from the positions |
| manifest | `anndata` | `source` (`h5ad` or `zarr`), `counts` (`layers/counts`, `raw/X` or `X`), `countsAreIntegers`, `positions` (`obsm <key>`), `libraries`, `pixelSize` (where it came from), `rows`, `rowsTotal` |
| manifest | `omezarr` | `ngffVersion`, `channels`, `zPlanes`, `zPlane`, `levels`, `dtype` (of the file) |
| manifest | `spatialdata` | `coordinateSystem` (`global` or `chosen`), `tables`, `imageFrom` (`images` or `density background`), `transforms`, and `layers[]`: `key`, `rows`, `counts`, `countsAreIntegers`, `region` (`circles` or `polygons`) |
| `expr[]` | `platform` | `visium` for a classic Visium AnnData (its barcodes are Visium barcodes), `anndata` for other AnnData layers, `spatialdata` for SpatialData tables |
| `molecules` | `platform`, `source` | `spatialdata`, `points` |

- No file, library, table, region, element or coordinate-system name is written to the export:
  the blocks above hold counts and fixed words only.
- Layer keys: `spots`, `cells`, or a Visium HD bin size (`016um`) as for Space Ranger output.
- Only `platform: "visium"` layers can get per-row id lists (`export --ids`, Visium spot
  barcodes); the `anndata` and `spatialdata` platforms are not on the allow-list, so their ids
  are refused there.
- The file's own `obs` clusterings and `X_umap` are kept in the cache (`userobs/<key>/obs.csv`,
  `umap.csv`, and the manifest's `userObs` list) and attached at export as your own clusters
  (`source: "user"` in `umap/index.json`). `userObs` and `userobs/` stay in the cache: the
  export's manifest does not have them.

## Molecules: `molecules/`

Written by `prepare.py xenium` (`molecules.py`, platform-neutral; `xenium.py` reads
`transcripts.parquet` when present, else `transcripts.zarr.zip` level 0). The cache holds one
file per gene (`gene/<NAME>.bin`, which `server.py` serves as it is); `export` packs them in
`index.json` order (`molecules/genepack/`, each gene's `s`, `o`, `b` on its `index.json` entry and
`genePack` on the manifest's `molecules` block, as described under "Gene packs") or, with
`--gene-files loose`, copies the files unchanged. A gene is fetched only when it is picked in
the Molecules panel.

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
  Genes with no molecules have no file (no place in the pack).
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

  Each SODA clustering also carries `"name": "Leiden 0.6"` and `"source": "soda"` (exports
  from 11 November 2026; older ones have neither, and the viewer names them by resolution).

#### The lab's own UMAP and clusterings (`userclusters.py`)

An export made with `--user-layer` and `--user-umap` / `--user-obs` / `--user-h5ad`
([PREPARING_DATA.md](PREPARING_DATA.md#your-own-umap-and-clusters-seurat-scanpy)) writes into
the same folder:

- `clusters_user<N>.bin`, one per imported clustering (N = 1, 2, ...): uint16 per row, the
  index into that clustering's `categories`, `65535` = not assigned (a row the user's file did
  not have, or had as NA). Not tied to `coords.bin`: a row can be assigned and not embedded,
  or the reverse.
- `coords.bin`: the user's embedding when one was given (NaN for rows without one), else
  SODA's, unchanged.
- `index_soda.json`, `coords_soda.bin`: SODA's own `index.json` and `coords.bin`, kept so an
  export without user files restores them and so the suggested cell types are scored on SODA's
  embedded bins. The viewer does not read them. With `--user-replace` SODA's
  `clusters_r*.bin` stay on disk but are not listed.

`index.json` then has, in addition:

| Field | Meaning |
|---|---|
| `clusterings[]` with `"source": "user"` | listed first; `id` (`"user1"` ...), `name` (the column name, de-identified), `file`, `n`, `categories` (the group names, in the user's order: a categorical's own order, else numbers first, then text), `colors` (`#rrggbb` per category, only if the user gave a palette), `sizes` (rows per category), `nAssigned`, `markers` (as above, over the assigned rows with counts). No `resolution` |
| `clusterings[]` with `"source": "soda"` | SODA's own, as above; left out with `--user-replace` |
| `default` | index of the clustering shown first: the first user one |
| `embedding` | `{"source": "user", "name": "<obsm key or the two column names>"}` or `{"source": "soda"}` |
| `nEmbedded`, `bounds` | of `coords.bin` as written (the user's embedding when there is one) |
| `userImport` | `layer`, `kind` (`csv` / `h5ad`), `barcodes`: one report per input file (`inFile`, `considered`, `matched`, `rate`, `offImage` = matched but outside a `--crop`, `methods` = how many matched exactly / by suffix / by prefix / by both, `otherSample` = left out by `--user-barcode-prefix` / `-suffix`), `sodaClusterings` (`kept`, `replaced` or `none`), `markers` (the rule) |

The manifest's `umap` entry lists every clustering as `{"id"?, "resolution"?, "name", "source",
"n"}` and has `"embedding": "user"` when the coordinates are the user's. In `celltypes.json`
(next section) a user clustering's suggestions are keyed by its `id` instead of the
resolution, and cover only its assigned rows that SODA embedded.

The viewer keys a clustering by `id` when present, else by `resolution`; it lists them by
`name`, names groups after `categories` (a name given with ✎ wins) and colours them with
`colors` when present.

## Single-cell / single-nucleus RNA-seq datasets

`prepare.py rnaseq` writes a different, smaller layout straight into the export folder:

```
<id>/manifest.json
<id>/umap.bin
<id>/obs/<NN>.bin
<id>/genes.json
<id>/genepack/<hash>-<i>.pack          every gene's record (below) end to end ("Gene packs"); s, o, b in genes.json
<id>/gene/<NAME>.bin                   the same, one file per gene (older exports, or --gene-files loose)
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
| `genes.json` | `[{"n": name, "p": fraction of cells expressing, "m": maximum value, "s", "o", "b"}]`, sorted by `p`, most widely expressed first; `s`, `o`, `b` = the gene's place in the pack (manifest `genePack`). |
| a gene's record (in the pack, or `gene/<NAME>.bin`) | uint32 `k`, then uint32[k] cell indices, then float32[k] values, for the cells where the gene is non-zero. |

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
positions and, from 11 November 2026, the barcode of each bin, by which `userclusters.py` matches a lab's own clusters and `cellshapes.py` a platform's or a lab's own cell outlines), `shapes/` (the cell outlines, already in their export format) and `overlays/`. `server.py` and `prepare.py export` read it. It holds the
full-resolution raw image and local file paths, so it should stay on your machine.
