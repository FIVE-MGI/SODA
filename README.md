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
- **Save and share.** **Copy link to this view** (Save tab) copies an address that opens the
  dataset as it looks now: position and zoom, the Tissue / Side by side / UMAP view, the
  layer, the genes on show with their colours and ranges, stains on or off with their colours
  and display ranges, molecule genes, outlines, the colour-by, the clustering with your
  cluster names and hidden clusters, the quality filter and the open tab. Drawn shapes are in
  it when they fit: a link is kept to about 2,000 characters, so long outlines are simplified
  or left out, and the message under the button says so. The view travels in the part of the
  address after `#`, which the browser never sends to a server; `ds=`, `qc=` and a `data=` or
  `static=` already in the address stay as they are. Not in a link: a box or freehand
  selection, marker cut-offs and cell-type rules, volcano and plot settings, and the result
  of a Draw a region run (the shapes come back; press Run UMAP). **Save JSON** writes the same
  view to a file together with the shapes in full, every setting and the run's result, and
  **Load JSON** on the same dataset puts it back and runs it again.
- **Compare side by side.** **Compare** (top bar, or the Save tab) opens another dataset in a
  second window beside the first. While **Link views** is on, moving or zooming one window
  moves the other to the same place at the same scale in micrometres, whatever the two
  datasets' pixel sizes and image sizes; **Link from here** keeps what each window shows now
  as the same place (to line up two sections that sit differently on their slides), **Line
  up centres** goes back to the middle of each picture. Genes picked, the colour-by and the
  Tissue / UMAP view follow where the other dataset has them, and the window says so where it
  has not. Up to four windows can be linked. Views are not rotated or mirrored, and the
  windows are separate: there is no split view inside one page. On a private site the second
  window does not ask for the password again. On a phone or tablet the other dataset opens
  in a new tab, not linked.
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
opens it. You need Python 3.10 or newer, Git, and an internet connection.

**1. Get the code and make a virtual environment.**

macOS / Linux:

```bash
git clone <this repository's URL> spatial-viewer
cd spatial-viewer
python3 -m venv .venv
source .venv/bin/activate
```

Windows (PowerShell):

```powershell
git clone <this repository's URL> spatial-viewer
cd spatial-viewer
py -m venv .venv
.venv\Scripts\Activate.ps1
```

(If PowerShell refuses to run the activation script, use Command Prompt instead and run
`.venv\Scripts\activate.bat`.) Once the environment is active, `python` means the
environment's Python on every system.

**2. Install, make the demo data, prepare and export it.** Same commands on every system:

```bash
pip install -r web/SpatialViewer/requirements.txt -r web/SpatialViewer/requirements-clustering.txt
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
| [docs/HOSTING.md](docs/HOSTING.md) | hosting options from a laptop to a password-protected cloud site; pre-upload checklist |
| [deploy/azure/README.md](deploy/azure/README.md) | Azure Static Web App + private blob container + password, step by step |
| [docs/DEIDENTIFICATION.md](docs/DEIDENTIFICATION.md) | `deid.py`, the private key file, the upload gate, and its limits |
| [docs/DATA_FORMAT.md](docs/DATA_FORMAT.md) | the exported file tree, byte by byte, for people writing their own tools |

## Requirements

- **Python 3.10 or newer** for preparing data (the pinned NumPy 2 and SciPy 1.14 need it).
  The automated test runs on Python 3.12 on Ubuntu. The tools are developed on Windows;
  macOS has not been tested by us.
- Packages: `web/SpatialViewer/requirements.txt` (reading images and counts, export, tile
  server). `web/SpatialViewer/requirements-clustering.txt` is optional and adds the
  whole-section UMAP and Leiden clusters (scikit-learn, umap-learn, igraph, leidenalg).
- Disk space: the prepare cache holds a compressed full-resolution copy of each image, and
  the export a tiled copy. Memory and time have not been benchmarked; large Visium HD
  datasets need the most (see [How big and how long](docs/PREPARING_DATA.md#how-big-and-how-long)).
- Viewing needs no installation, only a browser and an internet connection, also for a copy
  on your own computer: the page loads its drawing library (deck.gl), and Draw a region its
  UMAP library (umap-js), from `cdn.jsdelivr.net`.

## Browser support

The viewer needs **WebGL2**. It is developed and tested in **Google Chrome** on desktop.
Other desktop browsers (Edge, Firefox, Safari) have not been tested; they may work if they
support WebGL2. On an iPhone every browser runs on Safari's engine (WebKit); what was tested
there is below. Password-encrypted dataset links need the page to be served over
`https://` or from `localhost`.

**Copy link to this view** and **Compare** were checked with automated tests in Chromium
and in Playwright's WebKit build (the engine behind Safari), at desktop size and at
390 × 844; they have not been tried in Safari itself, in Firefox, or on a real phone.
A link made in one browser opens in the other.

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
web/SpatialViewer/    the page (index.html, viewer.js, ...) and the Python tools (prepare.py, export.py, ...)
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
tests/viewlink_test.js  the link that restores a view: round trips, the length budget, damaged and hostile links (node only; the smoke test runs it)
tools/label_editor.html  place UMAP cluster labels by hand (writes labelPositions for config.js)
```

## Citing

Manuscript in preparation. Until it is published, please cite this repository by its URL.
<!-- CITATION PLACEHOLDER: replace with the paper reference and DOI once available. -->

## License

©2026 Washington University. Non-commercial, non-clinical research use; see
[LICENSE](LICENSE). For a commercial license, contact Washington University's Office
of Technology Management at otm@wustl.edu.
