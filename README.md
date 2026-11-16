# SODA: Spatial Omics Data Annex

SODA shows very large spatial biology datasets in a web browser:
multiplexed immunofluorescence, H&E, Visium and Visium HD spatial transcriptomics,
Xenium (stains, cells and single transcripts), and single-cell RNA-seq UMAPs. Each slide
is converted once into a pyramid of small image tiles plus compact per-gene files. For the
image, the browser fetches only the tiles on screen, so a 30-channel whole slide opens in
seconds and pans without downloading the rest. Genes are different: picking a gene downloads
its whole file (its value in every bin or cell of the section, or every molecule of it).
The output is plain static files: you can look at it on
your own computer, put it on any web host for public data, or keep it in a private cloud
container behind a password. No database and no always-on server are needed, but the page
needs an internet connection: it loads its libraries from a public CDN.

<!-- SCREENSHOT PLACEHOLDER: docs/img/overview.png
     Whole-slide multiplex IF view with the channel panel open. Lead adds the image. -->

## Show your own data in 4 steps

No coding needed: you paste a few commands and fill in two lines.

1. **Install** with one command: `bash install.sh` (macOS / Linux) or
   `powershell -ExecutionPolicy Bypass -File install.ps1` (Windows). It ends with **SODA is
   installed**. Then activate it and run `pip install jupyter` ([details](#try-it-in-5-minutes)).
2. **Convert.** Open `notebooks/SODA_pipeline.ipynb` and fill in your work folder (`WORK`) and
   your cloud folder's link (`CLOUD_URL`). Run Sections 0 and 1, put your instrument files in
   the `inputs/` folders Section 1 makes (one subfolder per dataset), then **Run All**.
3. **Upload.** Answer Section 7's one question (is this data public?), then upload the
   `export/` folder, as it is, to your own storage: Azure, AWS, Google Cloud or a web server.
   Box, Google Drive and Dropbox do not work.
4. **Share** the link the notebook prints: `<the SODA page>/index.html?data=<your folder>`.

Every step in detail, with the one storage setting (CORS) the link needs:
[docs/SHARE_YOUR_DATA.md](docs/SHARE_YOUR_DATA.md). Everything beyond that: [docs/ADVANCED.md](docs/ADVANCED.md).

## Features

- **Four tabs, one job each.** The controls sit beside the picture in four tabs: **View**
  (image and stains, genes, molecules, outlines), **Clusters** (colour by clusters or cell
  types, select, quality filter, slide QC and counting tools for protein images),
  **Analyze** (draw a region, volcano plot, heatmap and violin) and **Save** (save / load,
  a link, dataset details, about and help). A tab shows only when the dataset has something
  for it. Sections fold, and a folded one says its state on its header ("3 stains on",
  "Leiden 0.6", "12% removed"); rarely used settings are under **More options** and each
  section's explanation is behind its **?**. The dataset picker and the Tissue / Side by
  side / UMAP switch are in a bar over the picture.
- **Multiplexed fluorescence.** Any number of channels, each with its own colour and
  display window (slider with histogram), blended on the GPU. Channel search, all on/off,
  nuclear stain pinned to the top, empty detector slots hidden.
- **H&E and other brightfield images** shown as a colour photo, alone or under
  expression data.
- **Spatial transcriptomics (Visium, Visium HD).** Search genes and show up to eight at
  once on the tissue, as points or as a density map, at 8 µm, 16 µm or coarser bins.
- **Xenium.** The morphology stains, the segmented cells as an expression layer (gene maps,
  UMAP and Leiden clusters, suggested names, QC, Draw a region and the volcano panel work on
  cells exactly as on bins), and the individual transcripts: pick genes in the Molecules
  panel and each one's molecules are drawn as coloured dots on the stains, several genes at
  once, fetched one gene at a time. Cell outlines are not drawn yet: each cell is a dot at
  its centre.
- **GeoMx DSP.** A study of hand-drawn regions (ROIs) rather than cells or bins: the scans
  (several slides side by side), the ROI outlines, and one disc per ROI segment coloured by a
  gene. The columns of the annotation sheet are the groups: legend, heatmap / dot plot /
  violin and the volcano compare groups of ROIs, and the QC panel has the GeoMx values (reads,
  saturation, nuclei, area, negative-probe background, LOQ). What works for ROI data is the
  comparing of groups; a UMAP, clusters, cell-type suggestions and density maps are not made
  for a few dozen regions, and the viewer says so instead of showing them. Counts are shown as
  counts per 10,000 in the plots (Q3 normalisation is not built in).
- **CosMx SMI and MERSCOPE.** The same three parts as Xenium, from the vendors' flat files:
  the per-FOV CosMx images stitched into one mosaic (raw stains or the rendered composite),
  or one z-plane of the MERSCOPE mosaic stains; the cells as an expression layer; and
  every transcript in the Molecules panel. Control probes and blanks are left out, and a
  dataset that comes without images still opens on a transcript-density background.
- **Stereo-seq (STOmics).** The registered H&E or ssDNA image with square-bin layers summed
  from the DNB counts of a GEF or GEM file (bin20 / bin50 / bin100 = 10 / 25 / 50 µm, cut to
  the tissue) and the segmented cells of the cell bin file as a layer of their own; gene maps,
  UMAP and clusters, QC, Draw a region and the plots work on both.
- **Tissue | Side by side | UMAP.** Each section's bins are embedded once into a UMAP
  with Leiden clusters (three resolutions, marker genes per cluster). Switch the main view
  between the tissue and the UMAP, or show both linked side by side: select bins in one
  (box or freehand) and the same bins light up in the other. **Show on slide** takes any
  selection, picked cluster or cell type back to the tissue and frames it there.
- **Suggested cluster names.** Each cluster gets a suggested cell type from editable
  marker lists, with a confidence and the runner-up shown; accept it, accept all
  high-confidence ones, or type your own. General lists for any tissue (genes and protein
  panels, with common marker-name aliases) are used by default; peripheral nerve lists are
  included too.
- **Draw a region (circle → UMAP).** Draw circles or freehand shapes on the tissue and get a UMAP, with
  optional Leiden clusters, of just the bins inside, computed in the browser. Drag a box or
  a loop over dots in its UMAP (Shift or **+** adds more) to light those bins on the slide,
  and **Show on slide** to frame them. Save the shapes and every parameter to redo it later.
- **Quality filter (per-bin QC).** Histograms of total counts, genes detected and % mitochondrial per bin,
  each with min / max thresholds. Off until you move a slider; removed bins are greyed on the
  slide and left out of Draw a region, and the thresholds are saved with the shapes and in
  the page link.
- **Save and share.** **Copy link** (Save tab) copies the page's address. It opens the same
  dataset (`ds=`) with the quality filter, if one is set (`qc=`), and keeps a `data=` or
  `static=` already in the address. Zoom, genes, stains, the Tissue / UMAP view, clusters
  and drawn shapes are not in it. To pass on a Draw a region result, use **Save JSON**: the
  file holds the shapes, every setting, the clusters and their names and the quality filter,
  and the other person loads it with **Load JSON** on the same dataset. Zoom, genes and
  stains are in neither.
- **Heatmap, dot plot and violin plots.** Compare genes across groups: the section's clusters
  at any resolution, the shapes you drew, the clusters of a Draw a region run, the regions of an overlay,
  or (single-cell datasets) any annotation column. Genes are each group's top markers plus any
  you type. A bin's value is log1p(counts per 10k); the heatmap shows the group mean (as it
  is, or z-scored / scaled 0–1 per gene), the dot plot adds the share of bins expressing, and
  violins show the spread with median and quartiles, with an "expressing only" switch because
  most bins hold zero. Every group shows its number of bins and the QC filter applies. Click a
  gene or a group to **show it on the slide**; save PNG, SVG, or the plotted numbers as CSV.
  Computed in your browser from the per-gene files (a typical cluster heatmap of a Visium HD
  section reads roughly 10–30 MB).
- **Multiplexed protein images (Akoya / PhenoCycler and other multi-channel TIFFs).** The
  slide is cut into 16 µm squares holding each marker's mean intensity (from full-resolution
  pixels), and the squares get everything above: UMAP, Leiden clusters, suggested names from
  protein marker lists, and Draw a region. Each marker is normalised per slide (background
  off, scaled to its 99.5th percentile, asinh), empty and duplicate channels are dropped
  automatically, and a slide QC box lists weak and saturated stains and outlines folds and
  bright edges. Set a cut-off per marker (by eye or Otsu) to read the % of tissue squares
  positive for it, or count squares per cluster. No cell segmentation (yet).
- **Region overlays.** GeoJSON regions (for example from QuPath) drawn over the image;
  click a region to restrict the gene density to its bins.
- **Single-cell RNA-seq.** UMAP from AnnData (`.h5ad`) or 10x output, coloured by cluster
  with automatically placed labels, plus up to twelve gene expression panels.
- **Host it yourself.** Static files work on any web host. Optional password protection
  for private data on Azure, or a small live tile server for a lab network.
- **De-identification tools.** Sample IDs are replaced by numbers through a private key,
  and a check command gates uploads. They know only our lab's ID formats: edit `_ID` in
  `deid.py` for yours first (see [DEIDENTIFICATION.md](docs/DEIDENTIFICATION.md)).

![The synthetic demo section in Side by side mode: tissue on the left, its UMAP on the right, both coloured by Leiden cluster](docs/img/demo_side_by_side.jpg)

*The synthetic demo dataset (below) in "Side by side" mode, coloured by Leiden cluster.*

<!-- SCREENSHOT PLACEHOLDER: docs/img/circle_umap.png
     Draw a region (circle -> UMAP): two circles on the tissue and the resulting UMAP panel. -->

## Try it in 5 minutes

This makes a small synthetic dataset (no real tissue or patient data), prepares it and
opens it. You need Python 3.12 or 3.13 (64-bit), Git, and an internet connection.

**1. Get the code and install it with one command.**

macOS / Linux:

```bash
git clone <this repository's URL> spatial-viewer
cd spatial-viewer
bash install.sh
source .venv/bin/activate
```

Windows (PowerShell or Command Prompt):

```powershell
git clone <this repository's URL> spatial-viewer
cd spatial-viewer
powershell -ExecutionPolicy Bypass -File install.ps1
.venv\Scripts\Activate.ps1
```

The install script finds Python 3.12 or 3.13, makes a virtual environment `.venv` in the
repository, installs the exact, tested package versions from
`web/SpatialViewer/requirements-lock.txt`, checks that they all import, and ends with
**SODA is installed**. It needs no administrator rights, takes a few minutes the first time
and is safe to run again. If no suitable Python is found, it says so and what to install.
(`-ExecutionPolicy Bypass` only lets this one script run; it changes no setting. If
PowerShell refuses the activation line, use Command Prompt and run
`.venv\Scripts\activate.bat`.) Once the environment is active, `python` means the
environment's Python on every system.

> **By hand instead** (or where the script stops: Intel Macs, Linux on ARM, Python 3.10 or
> 3.11): make and activate a venv yourself, then
> `pip install -r web/SpatialViewer/requirements.txt -r web/SpatialViewer/requirements-clustering.txt`.
> That takes the newest versions the ranges allow, not the tested set; with Python 3.12 or
> 3.13, `pip install -r web/SpatialViewer/requirements-lock.txt` gives the tested set.

**2. Make the demo data, prepare and export it.** Same commands on every system:

```bash
python examples/make_demo.py --out demo/raw
python web/SpatialViewer/prepare.py visium --id demo_visiumhd --image demo/raw/visiumhd/image.tif --outs demo/raw/visiumhd/outs --out demo/cache
python web/SpatialViewer/prepare.py tiff   --id demo_multiplex --src demo/raw/multiplex.ome.tif --out demo/cache
python web/SpatialViewer/prepare.py export --id demo_visiumhd  --out demo/cache --dest demo/export --markers examples/demo_markers.json
python web/SpatialViewer/prepare.py export --id demo_multiplex --out demo/cache --dest demo/export --fmt jpg --markers examples/demo_protein_markers.json
```

`make_demo.py` writes a Space Ranger-shaped Visium HD folder (`demo/raw/visiumhd/`: an
H&E-like image and `outs/` with 8 µm and 16 µm bins) and a small 6-channel fluorescence
image (`demo/raw/multiplex.ome.tif`). Everything goes under `demo/`, which Git ignores.

**3. Start two local web servers**, each in its own terminal window (activate the
environment in each). One serves the data, one serves the page:

```bash
python web/SpatialViewer/static_serve.py demo/export 8765
```

```bash
python -m http.server 8766 --directory web/SpatialViewer
```

**4. Open** this address in Chrome:

```
http://localhost:8766/index.html?static=http://localhost:8765&ds=demo_visiumhd
```

Search a gene (for example `MARKER_FASC1`), then try the **Tissue | Side by side | UMAP**
switch. Change `ds=demo_visiumhd` to `ds=demo_multiplex` for the fluorescence image.
Press Ctrl+C in both terminals to stop.

To check that everything above works without opening a browser:

```bash
python tests/smoke_test.py
```

It runs the same steps and then checks every output file.

## Use your own data

1. **[Prepare your data](docs/PREPARING_DATA.md):** one command per input type
   (Visium HD, Visium, Xenium, CosMx, MERSCOPE, Stereo-seq, GeoMx DSP, OME-TIFF / QPTIFF, per-marker TIFF folders, H&E,
   `.h5ad` or 10x RNA-seq), then `prepare.py export`.
2. **[De-identify and check](docs/DEIDENTIFICATION.md)** before anything leaves your
   computer.
3. **[Share it](docs/SHARE_YOUR_DATA.md), one of two ways:**
   - **Share a link.** Upload the export folder to your own storage and add its address
     to the SODA page's link (`index.html?data=<your folder>`). Names, groups and which
     stains show are set in one text file next to your data. No code, no website of your own.
   - **Make it your own.** Copy SODA on GitHub ("Use this template") and get your own
     SODA website that you can change however you like, still reading your data from your
     storage.

   For everything else (your own computer, any web server, privately on Azure with a
   password) see [Hosting](docs/HOSTING.md) and the [Azure steps](deploy/azure/README.md).

The page reads `web/SpatialViewer/config.js` for your dataset list, institution name,
logos and password endpoint. Copy `config.example.js` to `config.js`; every option is
explained in that file. `config.js` is listed in `.gitignore`.

**Rather not set it up yourself?** We can prepare and host your datasets for you. Please
email wbuchser@wustl.edu with any questions about hosting.

## Documentation

| | |
|---|---|
| [docs/PREPARING_DATA.md](docs/PREPARING_DATA.md) | commands and options per input type, export, clustering, the circle tool, common problems |
| [docs/SHARE_YOUR_DATA.md](docs/SHARE_YOUR_DATA.md) | the two easy ways: a link to your data folder, or your own copy of SODA |
| [docs/ADVANCED.md](docs/ADVANCED.md) | the advanced way in one place: your own site, passwords, labels and colours, every notebook setting, your own clusters, marker lists, which code file does what |
| [docs/HOSTING.md](docs/HOSTING.md) | hosting options from a laptop to a password-protected cloud site; pre-upload checklist |
| [deploy/azure/README.md](deploy/azure/README.md) | Azure Static Web App + private blob container + password, step by step |
| [docs/DEIDENTIFICATION.md](docs/DEIDENTIFICATION.md) | `deid.py`, the private key file, the upload gate, and its limits |
| [docs/DATA_FORMAT.md](docs/DATA_FORMAT.md) | the exported file tree, byte by byte, for people writing their own tools |

## Requirements

- **Python 3.12 or 3.13 (64-bit)** for preparing data. Install everything with one command:
  `bash install.sh` (macOS / Linux) or `powershell -ExecutionPolicy Bypass -File install.ps1`
  (Windows), see [Try it in 5 minutes](#try-it-in-5-minutes). The script installs
  `web/SpatialViewer/requirements-lock.txt`, which pins every package to one version. That
  set is tested end to end on Windows (Python 3.13), and the automated test installs and
  runs it on Ubuntu (Python 3.12) on every push. Every version in it is also published for
  macOS 12 or newer on Apple Silicon and for Linux x86_64 with glibc 2.28 or newer, but we
  have not run it on a Mac. Intel Macs, Linux on ARM and Python 3.14 are not covered by the
  lock file yet.
- Packages, as ranges (the human-edited lists the lock file is made from):
  `web/SpatialViewer/requirements.txt` (reading images and counts, export, tile server) and
  `web/SpatialViewer/requirements-clustering.txt` (optional: whole-section UMAP and Leiden
  clusters with scikit-learn, umap-learn, igraph, leidenalg). Installing these by hand needs
  Python 3.10 or newer (NumPy 2, SciPy 1.14) and gets the newest versions they allow.
- **Memory (RAM):** no minimum has been benchmarked. Measured peaks: the demo above (its
  Visium HD part, on Windows) 0.6 GB; a CosMx slide's counts (26 million transcripts, without its images) 0.76 GB; a
  whole mouse brain Stereo-seq chip 2.7 GB to prepare (3.9 GB from `tissue.gem.gz`) and
  **8.5 GB to export, the largest measured so far**. For whole-section datasets, have more
  free memory than that. Visium HD with its 2 µm layer has not been measured and is
  expected to need the most; if it runs out, leave that layer out with `--bins 008um,016um`
  (see [How big and how long](docs/PREPARING_DATA.md#how-big-and-how-long)).
- Disk space: the prepare cache holds a compressed full-resolution copy of each image, and
  the export a tiled copy. Time has not been benchmarked on a reference machine; measured
  examples are in [How big and how long](docs/PREPARING_DATA.md#how-big-and-how-long).
- Viewing needs no installation, only a browser and an internet connection, also for a copy
  on your own computer: the page loads its drawing library (deck.gl), and Draw a region its
  UMAP library (umap-js), from `cdn.jsdelivr.net`.

## Browser support

The viewer needs **WebGL2**. It is developed and tested in **Google Chrome** on desktop.
Other desktop browsers (Edge, Firefox, Safari) have not been tested; they may work if they
support WebGL2. On an iPhone every browser runs on Safari's engine (WebKit); what was tested
there is below. Password-encrypted dataset links need the page to be served over
`https://` or from `localhost`.

**Phones and tablets** open the same page and link; the layout follows the screen.

- On a phone held upright the map fills the screen and the four tabs are a bar along the
  bottom: tap a tab to open its controls in a sheet over the lower half, tap it again (or
  the sheet's handle) to close it, drag the handle for more or less of it. The dataset
  picker and the view buttons float at the top of the map. Side by side puts the tissue
  above the UMAP. Tablets and phones on their side keep the tabs and the panel at the left.
- Touch: pinch to zoom, drag to pan, tap a spot to read it (what hovering shows on a
  computer). While drawing a shape, one finger draws and two fingers move the picture.
  "Add to the selection" and "Done drawing" stand in for Shift and Esc.
- A phone has far less memory and a slower single core than a desktop, so on a touch
  screen the viewer draws at most 2× the screen's CSS pixels, keeps fewer image tiles,
  notes when more than 10 stains are on (a phone can run out of memory and reload the
  page), and Draw a region samples at most 20,000 bins (50,000 on a computer) and keeps
  the screen on while it runs. Heatmaps and violins take up to 40 genes at a time (150 on a
  computer) and ask before downloading more than 30 MB of gene files. A shape file saved on a phone records that limit; rerun on
  a computer, a selection of more than 20,000 bins is sampled up to 50,000 instead.

Touch handling and the image path were tested in Chrome's phone emulation and on an
iPhone (so on WebKit, Safari's engine), opening a password-protected Azure deployment. The
tabbed layout (the bottom bar and its sheet) has so far been tested in Chromium's phone
emulation only, not on a real phone or on WebKit.

**If the image does not appear** (most often on a phone), the page now says why in the
status line: tiles the data store refused, tiles that did not download or decode, no
WebGL2 (iOS before 15, or Lockdown Mode), or graphics memory the browser took back. Tiles
that failed are asked for again by themselves, and again when the page comes back to the
front. For the rest, add `&diag=1` to the address after signing in as usual
(`.../index.html?ds=<id>&diag=1`): a panel lists the browser and iOS version, WebGL2 and
its limits, one manifest and one tile fetched the way the viewer fetches them, each way of
decoding that tile, whether it reaches a texture and draws, what the map holds, and the
first errors caught. **Copy** puts the report on the clipboard. It never shows the
password or the access token: everything after `?` in an address is dropped. `&decode=img`
(or `&decode=plain`) makes the viewer decode tiles through an `<img>` element (or without
createImageBitmap's options), to compare on a device where the default fails.
`tests/phone_browser_test.js` runs these cases in Playwright's WebKit and Chromium; that
WebKit is the engine, not iOS Safari, so a real iPhone is still the only proof for one.

## Repository layout

```
install.sh / install.ps1  one-command install (macOS, Linux / Windows) into .venv
web/SpatialViewer/    the page (index.html, viewer.js, ...) and the Python tools (prepare.py, export.py, ...)
web/SpatialViewer/requirements-lock.txt  every package at one tested version (made from requirements*.txt)
web/SpatialViewer/markers/  marker-gene lists for suggested cluster names (edit for your tissue)
deploy/azure/         Azure Static Web Apps config and the /api/token function
docs/                 documentation
examples/make_demo.py synthetic demo data
examples/make_xenium_demo.py  synthetic Xenium output folder (cells, transcripts, stains)
examples/make_geomx_demo.py   synthetic GeoMx DSP study (DCC, PKC, annotation sheets, scans with ROIs)
examples/make_imaging_demo.py  synthetic CosMx and MERSCOPE folders (FOV images / mosaics, cells, transcripts)
examples/make_stereoseq_demo.py  synthetic Stereo-seq output (bin GEF, cell bin GEF, GEM, registered image)
tests/smoke_test.py   end-to-end test
tests/phone_browser_test.js  the image path in phone browsers, with each failure injected (needs Playwright)
tests/diag_test.js    the ?diag=1 panel never prints tokens (node only; the smoke test runs it)
tools/label_editor.html  place UMAP cluster labels by hand (writes labelPositions for config.js)
tools/check_install.py   checks that the Python packages import, with the right versions
```

## Citing

Manuscript in preparation. Until it is published, please cite this repository by its URL.
<!-- CITATION PLACEHOLDER: replace with the paper reference and DOI once available. -->

## License

©2026 Washington University. Non-commercial, non-clinical research use; see
[LICENSE](LICENSE). For a commercial license, contact Washington University's Office
of Technology Management at otm@wustl.edu.
