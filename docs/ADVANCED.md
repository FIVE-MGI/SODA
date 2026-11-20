# The advanced way

The [quick start in the README](../README.md#quick-start) and
[SHARE_YOUR_DATA.md](SHARE_YOUR_DATA.md) are enough to share your data. This page is for
everything past that: your own website, a password, your own names and colours, the
notebook's settings, your own clusters, marker lists, the code, a batch run over many datasets, and offline use.
Each section is short
and links to the full reference.

1. [Your own copy and website](#1-your-own-copy-and-website)
2. [Private data with a password](#2-private-data-with-a-password)
3. [Your own labels, groups, colours and stains](#3-your-own-labels-groups-colours-and-stains)
4. [Every notebook setting](#4-every-notebook-setting)
5. [Your own clusters and UMAP](#5-your-own-clusters-and-umap)
6. [Marker lists for suggested names](#6-marker-lists-for-suggested-names)
7. [Changing the code](#7-changing-the-code)
8. [Many datasets at once, without the notebook](#8-many-datasets-at-once-without-the-notebook)
9. [Offline use](#9-offline-use)

---

## 1. Your own copy and website

Copy SODA on GitHub, turn on its free website and point it at your data folder with one
line in `config.js` (`data: "<your export folder's link>"`). The steps are in
[SHARE_YOUR_DATA.md, "Make it your own"](SHARE_YOUR_DATA.md#2-make-it-your-own), including
how to pull in SODA's later fixes.

Other places the page can live:

- your own computer, for a look before anything is uploaded:
  [HOSTING.md A](HOSTING.md#a-on-your-own-computer) (the notebook's Section 6 does this for you);
- any web server or static host: [HOSTING.md B](HOSTING.md#b-any-static-web-host-public-data);
- a live tile server on your lab network: [HOSTING.md E](HOSTING.md#e-optional-live-tile-server-serverpy).

Please put `soda` in your site's name, for example `yourlab-soda`.

## 2. Private data with a password

| Way | How it works | Steps |
|---|---|---|
| **Password site on Azure** (what runs in production) | The data stays in a private blob container. A small function, `/api/token`, checks one shared password and gives the browser a read-only pass that expires after an hour by default (`TOKEN_HOURS`, at most 24). The storage key stays on the server. In `config.js`: `tokenEndpoint: "/api/token"`. | [HOSTING.md C](HOSTING.md#c-azure-private-data-with-a-password), then [deploy/azure/README.md](../deploy/azure/README.md) |
| **A link with a read pass** | The folder stays private; the link carries a read-only pass (an Azure SAS, an S3 pre-signed link). In the notebook: `READ_PASS`. Anyone holding the link can read the data until the pass expires, so treat the link like a password. | [SHARE_YOUR_DATA.md Step 5](SHARE_YOUR_DATA.md#step-5-make-the-link) |
| **AWS or another cloud** | Your own token endpoint that behaves like the Azure one. A design, not built or tested. | [HOSTING.md D](HOSTING.md#d-aws-s3--cloudfront-or-another-backend) |

Keep in mind:

- On a password site (`tokenEndpoint` set) the page does not open `?data=` links from
  outside, unless `config.js` sets `allowDataLinks: true`.
- `config.js` is sent to every visitor before the password is asked. Put no name in it that
  identifies a person.
- Before any upload, go through the [checklist in HOSTING.md](HOSTING.md#before-you-upload-checklist).

## 3. Your own labels, groups, colours and stains

There are two places for these settings:

- **`datasets.json`**, next to your data (written by `catalog.py` or the notebook). It works
  on the shared SODA page and on your own copy. Edit it in any text editor and upload only
  that file again.
- **`web/SpatialViewer/config.js`**, in your own copy only. Copy it from `config.example.js`,
  which explains every option. If `config.js` has `data:`, the dataset list, names and
  groups come from that folder's `datasets.json`, so set them there.

| What | In `datasets.json` | In `config.js` |
|---|---|---|
| title at the top of the page | `title` | `brand.title` |
| line under the title | `byline` | `brand.byline` |
| "Link to paper" | `paper` | `paperUrl` |
| logos | (none) | `brand.logos`; make a folder `web/SpatialViewer/logos/` for the files (git ignores it: add them with `git add -f`) |
| dataset opened first | `default` | `defaultDataset` |
| a dataset's name in the list | `label` | `datasets.<id>.label` |
| the heading it sits under | `group`: the heading text; headings appear in the order they first occur | `datasets.<id>.group`: a key of `groups`, whose order sets the order of the headings |
| stains left out | `hideChannels` (exact names, per dataset) | `datasets.<id>.hideChannels`; the top-level `hideChannels` is one pattern for every dataset (default: open, blank, empty, unused, none) |
| stains on when a dataset opens | `defaultChannels` (part of a name is enough) | `datasets.<id>.defaultChannels`, or `defaultChannels` for every dataset |
| UMAP cluster label positions | `labelPositions` | `labelPositions`, keyed by dataset id; made with `tools/label_editor.html` |
| gene colours | (none) | `geneColors` |
| molecule colours | (none) | `moleculeColors` |
| cluster colours | come with your own clusters (see [section 5](#5-your-own-clusters-and-umap)) | `catMinLuminance` and `labelMinLuminance` set how dark a colour may be |

The full `datasets.json` example is in
[SHARE_YOUR_DATA.md](SHARE_YOUR_DATA.md#what-you-can-change-in-datasetsjson). Never change a
dataset's `id`: it is its folder's name. An `about.html` next to `datasets.json` is shown
under "About".

Stain colours and display windows are picked in the page (View tab) and are not saved in
either file. The other `config.js` options (how many genes at once, the window slider, the
tile cache, how many molecules are drawn) are explained one by one in
[`config.example.js`](../web/SpatialViewer/config.example.js).

## 4. Every notebook setting

All in [`notebooks/SODA_pipeline.ipynb`](../notebooks/SODA_pipeline.ipynb) (the code behind its cells is
[`notebooks/soda_notebook.py`](../notebooks/soda_notebook.py)). Only the first
two need filling in; the rest have defaults, in the "Advanced settings (optional)" cell of
Section 0. Run that cell even if you change nothing.

| Setting | Default | What it does |
|---|---|---|
| `WORK` | `D:\SODA_work` | one folder on a big disk; `inputs/`, `cache/` and `export/` go inside. Plan for 2-3 times the size of your raw files |
| `CLOUD_URL` | empty | your cloud folder's address, without any key or `?sv=...` pass. Needed before Section 7 |
| `KEY_FILE` | the file named by the `SPATIALVIZ_KEY` environment variable, else `sample_id_key.csv` in `private/soda/` in your home folder | the de-identification key. Keep it outside `WORK` and the repository; never upload it. Section 5 checks for sample IDs only when it exists |
| `ALWAYS_CHECK_IDS` | `False` | `True` runs the sample ID check even without a key file |
| `SODA_PAGE` | `"https://soda.mgifive.org/"` | the address of the SODA page your link opens: the shared SODA page, or your own copy's address. Section 7 builds the link from it |
| `TITLE` | empty | the name of the whole collection at the top of the page (passed to `catalog.py`); empty keeps the title already in `datasets.json` |
| `READ_PASS` | empty | private folder only: a read-only pass (Azure SAS, starts with `sv=`), added to the link. Never a write pass |
| `TILE_FORMAT` | `"jpg"` | image tiles as JPEG (several times smaller) or `"png"` (lossless, for exact intensities) |
| `VISIUM_BINS` | `"008um,016um"` | the Visium HD bin sizes to keep. Add `002um` for the 2 µm layer (very large) |
| `WORKERS` | `8` | how many tile writers run at once during the export |
| `MARKERS` | `{"visium": "", "proteomics": ""}` | marker lists for suggested cluster names (see [section 6](#6-marker-lists-for-suggested-names)); empty uses the general list. `"visium"` is the gene list (Visium, Xenium, CosMx, MERSCOPE, Stereo-seq, GeoMx), `"proteomics"` the protein list (multichannel protein images) |
| `INCLUDE_IDS` | `True` | the barcode / cell id lists (`expr/<layer>/ids.txt.gz`) in the export, checked for sample IDs; `False` passes `--no-ids` and writes none (see [below](#barcode-and-cell-id-lists)) |
| `REBUILD` | `False` | `True` prepares datasets again even when `cache/<id>/` has them; also delete `export/<id>/`, or old tiles stay |

Two more, outside Section 0:

- **`OPTIONS`** (Section 2): settings per dataset, such as `label`, `pixel_size`, `crop`,
  `flip`, `channels` or `user_layer`. Every key:
  [PREPARING_DATA.md, Settings per dataset](PREPARING_DATA.md#settings-per-dataset-options). For the
  standard formats (spatial AnnData, OME-Zarr, SpatialData) there are also `library`,
  `matrix`, `spatial_key`, `image_res`, `layer_kind`, `coords_um`, `series`, `table`,
  `region`, `coordinate_system`, `image_element`, `points` and `no_user_clusters`, the
  `prepare.py` options of the same names
  ([PREPARING_DATA.md](PREPARING_DATA.md#standard-formats-spatial-anndata-ome-zarr-spatialdata)).
- **`IS_THIS_DATA_PUBLIC`** (Section 7): `"yes"` or `"no"`. The cell stops until you answer.

The same steps without the notebook, one command each:
[PREPARING_DATA.md](PREPARING_DATA.md). All of them for a whole folder in one command:
[section 8](#8-many-datasets-at-once-without-the-notebook).

### Barcode and cell id lists

**Included by default**, as Loupe Browser, Xenium Explorer and Vitessce keep them: every export
holds each layer's barcodes or cell ids in `expr/<layer>/ids.txt.gz` (one per bin or cell, in
row order), which the viewer's *Save barcodes (CSV)* reads. This holds for `prepare.py export`,
`prepare.py batch`, the notebook and the desktop app.

- **They are checked for sample IDs.** Before a list is written, `deid.py`'s sample-ID patterns run
  over every id; for Visium / Visium HD, Xenium, CosMx, MERSCOPE and Stereo-seq each id must also
  match the instrument's own id format in full (so a merge prefix such as `S1_`, or a sample name
  glued to a barcode, fails). A layer that fails gets no list, a `***` line in the export log says
  why (the layer, the reason and how many ids, never an id), and the export goes on. AnnData and
  SpatialData ids (obs names) have no format check: only the sample-ID patterns guard them.
- **You are responsible for what you publish.** The checks know only the ID patterns in
  `deid.py` (edit `_ID` there for your lab, [DEIDENTIFICATION.md](DEIDENTIFICATION.md)) and the
  formats above. Read your ids before you share a public link.
- **To turn them off:** `prepare.py export ... --no-ids`, `prepare.py batch ... --no-ids` (also
  `SODA prepare ... --no-ids`), `INCLUDE_IDS = False` in the notebook's advanced settings, or untick
  *Include barcode / cell id lists* under More options in the desktop app. An export made without
  them can get them later with `prepare.py ids --id <id> --out <cache> --dest <export folder>`.
- `--ids` is still accepted (old commands work). It asks for the lists explicitly: then a failed
  check stops the export and nothing is written, instead of leaving one list out.
- To hold another platform's ids to a format, add a regular expression that matches the whole id
  to `FORMATS` in `web/SpatialViewer/layerids.py` (and its name to `platform_of` if the cache
  calls it something else).

## 5. Your own clusters and UMAP

If you already clustered a Visium HD or Visium dataset in Seurat or scanpy, SODA can show
your clusters, with your names and colours, next to its own (or instead of them). In the
notebook, put `my_umap.csv` and `my_clusters.csv` (or one `my_analysis.h5ad`) next to
`outs/` and set `user_layer` for that dataset in `OPTIONS`. On the command line, use the
`--user-*` options of `prepare.py export`, or `prepare.py userclusters` on an export that
already exists. A column `<name>_color` with `#rrggbb` values (or the `.h5ad` palettes) sets
the colours.

File formats, how barcodes are matched and how to write the files from R or Python:
[PREPARING_DATA.md, "Your own UMAP and clusters"](PREPARING_DATA.md#your-own-umap-and-clusters-seurat-scanpy).

### Load your own results in the viewer (no new export)

A reader can also bring results to a page that is already shared, with no new export:
**Clusters tab, Your own results, Load your own results...** (or drop the file anywhere on
the page). The file is read in the browser and never uploaded.

**The file.** A CSV or TSV, gzipped or not (`.csv`, `.tsv`, `.txt`, `.gz`), one row per bin
or cell, with:

- a **barcode / cell id** column (or R row names, as `write.csv` writes them);
- any of: two **UMAP** columns (also tSNE, PCA, PHATE ...: `UMAP_1`/`UMAP_2`, `umap_1`,
  `X_umap1`, `tSNE_1`, `PC_1` ...), **groupings** (`seurat_clusters`, `leiden`, `louvain`,
  `cell_type`, `sample` ..., text or whole numbers), **numbers** (scores, QC), and a
  **colour per group**: a column named `<grouping>_color` (or `_colour`, `.colors`) holding
  `#rrggbb`.

Each column's use is guessed from its name and values. The section lists every column with a
drop-down to change it (Not used, Barcode / cell id, Row number, UMAP x, UMAP y, Groups, Value,
Colours of a grouping). Under a grouping you can paste a colour list in group order, as scanpy
keeps it in `uns["<grouping>_colors"]`. Groups are ordered numerically, then alphabetically.

**Matching.** Rows are matched to the bins through the export's barcode list
(`expr/<layer>/ids.txt.gz`; see [DATA_FORMAT.md](DATA_FORMAT.md)). Matching is by exact id, then
with a `-1` style ending set aside (Seurat, scanpy `concat`), then with a sample prefix set aside
(`sampleA_` from Seurat `merge`). The section shows how many rows matched, for example "9,812 of
10,000 rows matched", and counts unmatched rows; it never lists them or any id. Below half
matched, nothing is applied until you press **Use anyway**. When a merged object holds several
samples, type the start or ending of the ids to keep ("Ids starting with", "Ids ending with"). An export without a barcode list says so, with the
command that adds one: `prepare.py ids --id <id> --out <cache> --dest <export root>`.

Without a barcode list there are two other ways. Each needs its own clear statement:

- a column `soda_row` (or `bin_row`) holding each bin's row number in the export (the order of
  `positions.bin`); or
- a first line `# soda row order`, then exactly one row per bin, in the export's order. A file
  with any other number of rows is refused.

Both show a warning, because nothing checks that the rows really belong to those bins.

**What you get.** The groupings appear first under **Clustering**, named `<column> (your
file)`. The first one is shown with its colours on the tissue and the UMAP, and works with
names, cell sets and the volcano plot, as SODA's own clusterings do. Numbers
appear in **View, Colour by a value**. If the file has a UMAP, **UMAP shown** picks between
yours and SODA's. A layer with no SODA UMAP uses yours. The results belong to one layer; each
layer can have its own.

**Saving and links.** **Save JSON** can keep the matched table (no ids and no file name, one
entry per matched bin). Its size is shown before you tick it, and Load JSON brings it back. A
**view link** never holds the table: it keeps a short fingerprint of the file and the grouping
or value on show. Whoever opens the link is asked to load the same file, and a different file is
reported, and the link's view is not applied until the right file is loaded.

**From Seurat (R).** Writes the cell metadata you choose, the UMAP and a colour per cell type:

```r
library(Seurat)
obj <- readRDS("my_object.rds")                     # your analysed object
df  <- obj@meta.data[, c("seurat_clusters", "nCount_Spatial")]   # any columns you want
df$cell_type <- Idents(obj)                         # or a metadata column of names
umap <- Embeddings(obj, reduction = "umap")[, 1:2]
colnames(umap) <- c("UMAP_1", "UMAP_2")
df <- cbind(df, umap[rownames(df), ])
# optional: your colours, one per cell type
pal <- setNames(scales::hue_pal()(nlevels(df$cell_type)), levels(df$cell_type))
df$cell_type_color <- unname(pal[as.character(df$cell_type)])
write.csv(df, gzfile("my_results.csv.gz"))          # row names = barcodes
```

**From scanpy (Python).** Writes `obs` columns, the UMAP and the colours scanpy keeps in `uns`:

```python
import scanpy as sc

adata = sc.read_h5ad("my_analysis.h5ad")            # your analysed object
df = adata.obs[["leiden", "total_counts"]].copy()   # any columns you want
df["X_umap1"] = adata.obsm["X_umap"][:, 0]
df["X_umap2"] = adata.obsm["X_umap"][:, 1]
# optional: your colours, one per group
if "leiden_colors" in adata.uns:
    pal = dict(zip(adata.obs["leiden"].cat.categories, adata.uns["leiden_colors"]))
    df["leiden_color"] = adata.obs["leiden"].map(pal).astype(str)
df.to_csv("my_results.tsv.gz", sep="	")            # index = barcodes; gzipped from the name
```

An `.h5ad` or `.rds` can't be read in the browser; write the CSV as above, or use
`--user-h5ad` when you export. Limits: 1 GB per file (512 MB unpacked) and 200 groups per
grouping. A column of whole numbers with more than 30 different values is taken as a value unless its
name says it is a grouping (`seurat_clusters`, `leiden` ...).

## 6. Marker lists for suggested names

Each cluster gets a suggested name from a list of marker genes (or proteins). The lists are
plain JSON files in [`web/SpatialViewer/markers/`](../web/SpatialViewer/markers/):

| File | For |
|---|---|
| `rna_generic.json` | genes, any tissue (the default for Visium and the other sequencing data) |
| `protein_generic.json` | protein panels, any tissue (the default for fluorescence images), with a table of marker-name aliases |
| `nerve.json` | genes, peripheral nerve |
| `brain.json` | genes, brain (mouse or human): excitatory and inhibitory neurons, astrocytes, oligodendrocytes, OPCs, microglia, border macrophages, vascular, meningeal, ependymal and choroid plexus cells |
| `akoya_nerve.json` | protein panel, peripheral nerve |

To write your own, copy one and edit it: a `name`, and a list of `types`, each with a
`name` and its `genes`. Gene names are matched without regard to case. Then use it:

- notebook: `MARKERS` in the advanced settings;
- command line: `prepare.py export --markers my_markers.json`;
- in the browser, for Draw a region on gene data: `markerSet` in `config.js`. Without that
  line the page uses `markers/rna_generic.json`; `config.example.js` sets
  `markers/nerve.json`, so change or remove it.

How a name is scored, and the confidence shown:
[PREPARING_DATA.md, "Suggested cluster names"](PREPARING_DATA.md#suggested-cluster-names-marker-lists).

## 7. Changing the code

In your own copy you can change anything. Start with the table in
[SHARE_YOUR_DATA.md, Step 4](SHARE_YOUR_DATA.md#step-4-change-what-you-want). The page is
plain HTML, CSS and JavaScript with no build step: edit a file, reload the page.

The page, in `web/SpatialViewer/`:

| File | What it does |
|---|---|
| `index.html` | every control, with its id; loads the scripts and styles below |
| `config.js` | your settings (copied from `config.example.js`) |
| `viewer.js` | the main viewer: opens datasets, draws image tiles, genes, molecules and clusters with deck.gl, the password and `?data=` links, the UMAP view |
| `layout.js` | the four tabs, folding sections, "More options", the "?" notes, the phone layout. It moves controls; it makes none |
| `dspicker.js` | the searchable dataset list |
| `qc.js` | the quality filter: per-bin histograms and thresholds |
| `roi.js` | what differs for GeoMx ROI layers: words, drawing, notes |
| `plots.js`, `plots_worker.js` | heatmap, dot plot and violin plots; the worker does the numbers |
| `volcano.js`, `de_worker.js` | the volcano plot; the worker runs the differential expression test |
| `lasso_worker.js` | Draw a region (circle → UMAP): PCA, UMAP, Leiden and markers in the browser; fast mode above 20,000 bins |
| `umapfit.js` | umap-js's fit (nearest neighbours, fuzzy graph, layout) over typed arrays, the same numbers as umap-js; used by `lasso_worker.js` |
| `leiden.js`, `leiden_worker.js` | Leiden clustering, used by `lasso_worker.js`; the worker runs the stability sweep's Leiden runs on spare cores |
| `diag.js` | the `&diag=1` diagnostic panel |
| `styles.css`, `layout.css`, `dspicker.css`, `plots.css`, `volcano.css` | how the page looks; each `.css` belongs to the `.js` of the same name, `styles.css` to the rest |

The Python tools that write the data are in the same folder. `prepare.py` is the entry
point (one subcommand per input type, plus `export`, `userclusters`, `inspect`, `detect`, `auto` and
`batch`, which is `batch.py`); each
other `.py` file starts with a note on what it does. The standard formats are four of them:
`stdformats.py` (the `anndata`, `omezarr` and `spatialdata` subcommands), `annstore.py` (reads an
AnnData's groups from `.h5ad` or zarr with h5py / zarr, without the anndata package), `ngff.py`
(OME-NGFF images and coordinate transformations) and `zarrv3.py` (a small reader for zarr v3
stores: `zarr.json`, sharding and the common codecs, since SODA's zarr stays at version 2). What
they write, file by file:
[DATA_FORMAT.md](DATA_FORMAT.md). The password function is
`deploy/azure/api/token/index.js`.

When you change the page:

- `index.html` loads each file with `?v=...` after its name (a date and a letter). Change
  that part when you change the file, so browsers fetch the new copy instead of an old one.
- Use only unsigned typed arrays in the page code. `deid.py check` reads the names of the
  signed ones as sample IDs and stops the upload.
- `index.html` allows scripts only from the page's own site and `cdn.jsdelivr.net` (its
  Content-Security-Policy line). A library from anywhere else needs that line changed. The two
  outside libraries are pinned by version and SHA-384 in `libs.js` ([Offline use](#9-offline-use)).
- The service worker (`sw.js`) keeps every file `index.html` loads with a `?v=` string, and every
  `"name.js?v=..."` a script names (the workers). Name a new worker the same way, with its
  `?v=`, or it will be missing offline.
- Test before you share: `python tests/smoke_test.py` checks the tools and the public example
  exports ([REPRODUCE.md, section 6](../REPRODUCE.md#6-the-smoke-test); and `tests/diag_test.js`
  when Node is installed). Then open your own data.

## 8. Many datasets at once, without the notebook

`prepare.py batch` does what the notebook's Sections 2 to 5 do, for a whole folder, in one
command. Use it on a server, in a scheduled job, or for a folder of many datasets.

```bash
python web/SpatialViewer/prepare.py batch --src D:/SODA_work/inputs --out D:/SODA_work/cache --dest D:/SODA_work/export
```

- `--src` holds one folder per dataset, as the notebook's `inputs/`, and the folder's name
  becomes the dataset id. Each folder is checked as `prepare.py detect` does. Only a
  **certain** one is run. Any other is skipped, and the log says why and what to add.
- `--options options.json` is the notebook's `OPTIONS` as JSON, with the same keys (and
  `markers`, a marker list for that dataset only). Relative paths are read from the options
  file's folder.
- A dataset can also be a single `.h5ad` file or `.zarr` folder (spatial AnnData, OME-Zarr,
  SpatialData) directly in `--src`: its id is the name without `.h5ad`, `.zarr` or `.ome.zarr`.
  A choice the file leaves open (which library, table or coordinate system) makes it
  skipped, with the option to add (`library`, `table`, `region`, `coordinate_system`,
  `image_element`) for `options.json`. A zipped Xenium output (`.zip`) works the same way: its
  id is the name without `.zip`, and it is read as it is (the stain images are copied out into
  `<out>/_unzipped/<id>/` while it is prepared, then removed; see
  [PREPARING_DATA.md, Xenium](PREPARING_DATA.md#xenium)). A zip next to a dataset folder of the
  same name (`run1.zip` beside `run1/`) is prepared too, as a second dataset with the id
  `run1_zip`, and the log warns that the two may be the same data (keep one); a zip beside its
  folder that cannot be read as it is (any other kind) is left out with a warning. A zip that
  fails the zip-slip / zip-bomb checks is skipped with the reason. A Xenium transcripts store on
  its own, `<name>_transcripts.zarr.zip`, is a dataset too (id `<name>_transcripts`): its
  molecules on a transcript-density background, no images and no cells.
- Names starting with `_` or `.` are left out (a place for what is not a dataset), and the log
  says how many `_` names there were. The notebook's Section 2 and the desktop app find the
  datasets with this same code (`batch.find_datasets`).

  ```json
  {"nerve_section1": {"label": "Normal nerve"},
   "two_runs": {"kind": "merscope", "cells_dir": "D:/SODA_work/inputs/two_runs/vpt-run-a"}}
  ```

- Each dataset is prepared (`prepare.py auto`) and exported (`prepare.py export`) in its own
  process, so one that fails or runs out of memory does not stop the others.
- Run the same command again after a fix. A dataset already in the export folder is skipped,
  and one already in the cache is not prepared again. `--rebuild` does them again (delete
  `export/<id>/` first, or old tiles stay). `--only id1,id2` runs just those. `--dry-run`
  shows what would run and writes nothing.
- At the end it writes `datasets.json` (`catalog.py`; `--title` names the collection) and runs
  the sample ID check on the export and on its own logs. With a key file (`--key`,
  `SPATIALVIZ_KEY`, or the notebook's default place) the check is a gate, as in the notebook.
  Without one it only warns.
- The notebook's export settings: `--fmt` (default `jpg`), `--bins`, `--workers`,
  `--markers-genes`, `--markers-proteins`.

**The log** goes to `logs/` next to the export folder (`--logs` to change it; never inside
the export, as it holds local paths): `batch_<date>_<time>.md`, a table to read, and the same
as `.json` and `.csv`. One row per dataset, then the last lines of each failure's error and
what to add for each skipped one. Each step's full output is in `<id>_prepare.log` and
`<id>_export.log`. Anything that looks like a sample ID is written as `[sample ID]`, and a
folder whose name looks like one is skipped (rename it).

| id | kind | status | time | peak memory | export size | note |
|---|---|---|---|---|---|---|
| bx_broken | xenium (molecules) | failed | 12 s | 0.08 GB |  | prepare failed (exit 1) |
| bx_cosmx | cosmx (molecules) | ok | 31 s | 0.12 GB | 17 files, 0.2 MB |  |
| bx_two_runs | merscope (cells+molecules) | skipped |  |  |  | not sure: ... a choice is open: 2 folders hold cell_by_gene + cell_metadata (segmentation runs); which one? |

**Exit code:** 0 when nothing failed (skipped datasets are listed, not errors), 1 when a
dataset or `catalog.py` failed, 2 when there was nothing to do or an option was wrong, 3 when
the sample ID check is a gate and found something.

**Memory.** `--jobs` (default 1) runs that many datasets at once, and each needs its own
memory: small datasets under 1 GB (a CosMx slide's counts 0.76 GB), whole-section data such
as a big Visium HD section or a whole Stereo-seq chip about 8.5 GB (the largest measured, see
the [README](../README.md#requirements)). Keep `--jobs 1` for whole sections unless the computer
has the number of jobs times 8.5 GB free. The batch prints the free memory at the start and
warns. With more than one job the BLAS threads are split between them (`--threads` sets them).

## 9. Offline use

SODA can open with no network: on a train, in a lab without internet, or in the desktop app.
Two parts make that work, and you can use either on its own.

### The two libraries from your own site: `prepare.py vendor`

The page needs two libraries that are not in this repository: deck.gl (it draws the map) and
umap-js (Draw a region's UMAP). By default the browser loads them from `cdn.jsdelivr.net`. To
serve them from your own site instead, run once:

```bash
python web/SpatialViewer/prepare.py vendor          # or: python web/SpatialViewer/vendor.py
```

It downloads the two pinned files, checks each one's size and SHA-384, and writes them to
`web/SpatialViewer/vendor/`, next to `index.html`:

```text
vendor/deck.gl@9.1.12/dist.min.js
vendor/umap-js@1.4.0/lib/umap-js.min.js
```

The page and the UMAP worker look in `vendor/` first and use the CDN only when a file is not
there. The browser checks the SHA-384 again before it runs a file (a changed or damaged file is
refused and the CDN copy is used). Other options:

- `--check` only checks what is in `vendor/` and downloads nothing (exit code 1 if anything is
  missing or wrong).
- `--from <folder>` copies the files from a folder with the same layout, for example another
  computer's `vendor/` on a USB stick, for a machine with no internet at all.
- `--dest <folder>` writes `vendor/` into another copy of the page (the folder that holds
  `index.html`), for example your website's folder.

`vendor/` is git-ignored: the files are not part of the repository. Copy the folder along with the
page when you publish it (or run the command where the site is built). The versions and hashes
are in `web/SpatialViewer/libs.js`; change them there only, and in `index.html`'s `vendor/`
line (`tests/smoke_test.py` checks that they agree). Opened as a file (`file://`), the page
cannot use `vendor/` (a browser does not check a file's hash there) and uses the CDN.

### The page keeps itself: after one visit, it opens with no network

When the page is served over `https://` or from `localhost` (`static_serve.py`, `server.py`,
the desktop app), the browser keeps a copy of it after the first visit (a service worker,
`sw.js`): `index.html`, every file it loads, the workers, the marker lists, deck.gl and
umap-js (from `vendor/`, or the CDN's copy with the same hash). Online, the page still loads
everything from the network as before; the kept copy is used only when the network fails.
When you publish a new version (the `?v=` strings in `index.html` change), the next visit keeps
the new version and deletes the old copy.

**Data is not kept unless you ask.** An export can be gigabytes, and a private site's data sits
behind a password, so nothing of a dataset is kept on its own. To use a dataset offline, open it,
look at what you will need (zoom in where you will want detail, add the genes), then
**Save tab > Offline use > Keep this dataset for offline use**. It shows the size first, then
keeps on this computer what the page has fetched from that dataset so far plus the files Draw
a region needs (each layer's `hvg/` chunks), at most 500 MB per dataset (`offlineMaxMB` in
`config.js`). **Update offline copy** adds what you have looked at since; **Remove offline
copy** and **Remove every offline copy** take it off again. Offline, a dataset that was not
kept says so instead of opening.

Never kept: a dataset of a password-protected site (`tokenEndpoint`), one whose `?data=` link
carries a read pass, and one from `server.py`.

Settings in `config.js`:

| | |
|---|---|
| `offline: false` | no service worker on this site, and what earlier visits kept is deleted |
| `offline: true` | on also on a password-protected site (the page only; its data is still never kept) |
| `offlineMaxMB: 500` | the most one kept dataset may hold |
| `offlineKeep: false` | no "Keep this dataset for offline use" (the page still keeps itself): for a page whose data is already on the reader's computer. The desktop app serves this |

By default offline use is on, except on a password-protected site (`tokenEndpoint`), where no
dataset opens without the password service anyway.

Limits:

- **A page opened as a file (`file://`) cannot keep itself**: browsers allow a service worker
  only on `https://` and `localhost`. Serve the folder instead
  (`python web/SpatialViewer/static_serve.py`, or any web server).
- The browser may delete kept copies when the disk is nearly full, and a private window keeps
  nothing after it closes.
- The copies belong to one browser on one computer, and to one address: a page opened from
  another address (another port on `localhost`, too) keeps its own.

**Clearing it:** Save tab > Offline use > Remove every offline copy removes the kept datasets. To
remove everything, the page as well, clear the site's data in the browser (Chrome and Edge:
the padlock or site icon left of the address > Site settings > Delete data; Safari: Settings >
Privacy > Manage Website Data; Firefox: the padlock > Clear cookies and site data), or set
`offline: false` in `config.js`, which removes it on the next visit. The `&diag=1` panel shows
whether the service worker is in charge and how many datasets are kept.
