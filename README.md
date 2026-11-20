# SODA: Spatial Omics Data Annex

> **Never done this before? [Start here](docs/GETTING_STARTED.md):** every click and command,
> from downloading SODA to a link that shows your data, for people new to coding and the cloud.
>
> **Rather not use code at all?** [Download the SODA app](https://github.com/FIVE-MGI/SODA/releases/latest)
> (Windows or Apple Silicon Mac): pick your data folder, click Prepare, and it writes the folder to upload.

SODA turns what a spatial biology instrument writes into a link that opens the whole slide
in a web browser. You run it once on your own computer, from a Python notebook; it writes
plain files that you upload to your own storage, and anyone with the link can look, with
nothing to install.

**Licence:** GNU Affero General Public License, version 3 or later (AGPL-3.0-or-later), the same licence as SiDoLa-NS. See [Licence](#licence).

<!-- SCREENSHOT PLACEHOLDER: docs/img/overview.png
     Whole-slide multiplex IF view with the channel panel open. Lead adds the image. -->

## Quick start

From instrument folder to shared link. Looking needs only a browser; preparing your own data
needs Python 3.12 or 3.13 (64-bit), Git and an internet connection. Nothing leaves your
computer until step 6.

1. **Look at the examples.** Four public datasets (Xenium, Visium HD, PhenoCycler and
   MERSCOPE), prepared with SODA, open in Chrome with nothing to install:

   ```
   https://soda.mgifive.org/index.html?data=https://fivestore.blob.core.windows.net/soda/examples/
   ```

   Pick a dataset, search a gene, then try **Tissue | Side by side | UMAP**. Where each one
   comes from: [REPRODUCE.md, section 9](REPRODUCE.md#9-the-public-datasets).
2. **Install** (once, a few minutes). In a terminal:

   ```bash
   git clone https://github.com/FIVE-MGI/SODA.git soda
   cd soda
   bash install.sh                                        # macOS / Linux
   powershell -ExecutionPolicy Bypass -File install.ps1   # Windows
   ```

   It ends with **SODA is installed**. Then turn the environment on and add Jupyter:

   ```bash
   source soda-env/bin/activate          # macOS / Linux
   soda-env\Scripts\activate             # Windows
   pip install jupyter
   ```

3. **Put your data in a work folder.** Pick a work folder on a big disk (here `D:\SODA_work`;
   on a Mac, for example `~/SODA_work`) and put each dataset in its own folder under its
   `inputs/` (any name: SODA works out what each one is). To practise on public data first,
   [REPRODUCE.md, section 9](REPRODUCE.md#9-the-public-datasets) downloads the datasets the
   examples were made from.
4. **Run the notebook.** `jupyter notebook notebooks/SODA_pipeline.ipynb`. In the
   *Settings* cell set `WORK` to the same work folder, leave `CLOUD_URL` empty, then
   **Run > Run All Cells**. The run stops at Section 7 with a question: that is expected.
   Measured times per kind of dataset: [docs/PERFORMANCE.md](docs/PERFORMANCE.md).
5. **Open the preview.** Section 6 starts a small local server and the cell after it stops it
   again, so after Run All click into Section 6's server cell and press **Shift+Enter**. Open
   the first link it prints in Chrome. Search a gene, then try
   **Tissue | Side by side | UMAP**. Run the next cell to stop the server.
6. **Share it.** Set your storage folder's link (`CLOUD_URL`) and run the notebook again.
   Section 5 checks the export for sample IDs when you have a key file. In Section 7 answer
   one question (is this data public?), upload the `export/` folder as it is to your own
   storage (Azure, AWS, Google Cloud or a web server; not Box, Google Drive or Dropbox;
   never `cache/` or `inputs/`), allow the SODA page to read it once (CORS), and share the
   link the notebook prints:
   `https://soda.mgifive.org/index.html?data=<your folder>`.

Every step in detail: [docs/SHARE_YOUR_DATA.md](docs/SHARE_YOUR_DATA.md). The same from the
command line, without the notebook: [Use your own data](#use-your-own-data). A whole
folder of datasets in one command:
[ADVANCED.md, section 8](docs/ADVANCED.md#8-many-datasets-at-once-without-the-notebook).

## Documentation

| Guide | What it covers |
|---|---|
| [docs/PREPARING_DATA.md](docs/PREPARING_DATA.md) | every input type, its command and options; export, clustering, Draw a region; common problems |
| [docs/SHARE_YOUR_DATA.md](docs/SHARE_YOUR_DATA.md) | the two easy ways to share: a link to your data folder, or your own copy of SODA |
| [docs/HOSTING.md](docs/HOSTING.md) | hosting from a laptop to a password-protected cloud site; what the host must do; pre-upload checklist |
| [deploy/azure/README.md](deploy/azure/README.md) | Azure, step by step: private container, password, Static Web App |
| [docs/ADVANCED.md](docs/ADVANCED.md) | your own site, passwords, names and colours, every notebook setting, your own clusters, marker lists, the code, batch runs |
| [docs/DEIDENTIFICATION.md](docs/DEIDENTIFICATION.md) | `deid.py`, the private key file, the sample ID check before upload, and its limits |
| [docs/DATA_FORMAT.md](docs/DATA_FORMAT.md) | the exported files, byte by byte, for people writing their own tools |
| [docs/DESKTOP_APP.md](docs/DESKTOP_APP.md) | the SODA desktop app for Windows and Mac (no Python needed) |
| [docs/PERFORMANCE.md](docs/PERFORMANCE.md) | measured times and memory, the machine they came from, and how to measure your own |
| [REPRODUCE.md](REPRODUCE.md) | how to reproduce the numbers and figures in the paper |
| [CHANGELOG.md](CHANGELOG.md) | what is in each version |
| [CONTRIBUTING.md](CONTRIBUTING.md) | reporting a problem, asking a question, sending a change |

## Troubleshooting

The common failures, and where the fix is written down.

- **The page says it "could not reach" your folder, or genes do not load.** Your storage
  does not yet let the SODA page read it. Set the CORS rule once, and allow the `Range`
  header: [SHARE_YOUR_DATA.md, Step 4](docs/SHARE_YOUR_DATA.md#step-4-let-the-soda-page-read-it-one-setting-once).
  The notebook's last cell tests this and prints the fix for your storage. A host that does
  not answer Range requests at all: [HOSTING.md, What the host must do](docs/HOSTING.md#what-the-host-must-do).
- **A Box, Google Drive, Dropbox or OneDrive link does not work.** Their share links open a
  preview page, not the files, and other sites may not read them. Use Azure, AWS, Google
  Cloud or a web server: [SHARE_YOUR_DATA.md, Step 3](docs/SHARE_YOUR_DATA.md#step-3-upload-the-whole-export-folder).
- **A `?data=` link opens nothing on a password-protected site.** A site with a password
  (`tokenEndpoint` in `config.js`) does not open `?data=` links from outside, unless
  `config.js` sets `allowDataLinks: true`: [ADVANCED.md, section 2](docs/ADVANCED.md#2-private-data-with-a-password).
- **The install stops, or says no suitable Python was found.** The install script needs
  Python 3.12 or 3.13, 64-bit. On Intel Macs, Linux on ARM, or with Python 3.10 or 3.11, use
  the manual install under [Try it in 5 minutes](#try-it-in-5-minutes). Python 3.14 is not
  covered yet ([Requirements](#requirements)). `python tools/check_install.py` says which
  package is missing or has the wrong version.
- **The image does not appear (often on a phone).** The status line says why; add `&diag=1`
  to the address for a full report: [Browser support](#browser-support).
- **Bins do not line up, out of memory, `[deid] REFUSING to publish`, and other messages
  while preparing:** [PREPARING_DATA.md, Common problems](docs/PREPARING_DATA.md#common-problems).

## How it works

Each slide is converted once into a pyramid of small image tiles plus compact per-gene data.
For the image, the browser fetches only the tiles on screen, so a 30-channel whole slide
opens in seconds and pans without downloading the rest. Genes are different: picking a gene
downloads all of that gene (its value in every bin or cell of the section, or every molecule
of it), as one HTTP Range request into a pack that holds a layer's genes in a few large
files. The output is plain static files: you can look at it on your own computer, put it on
any web host for public data, or keep it in a private cloud container behind a password. No
database and no always-on server are needed. By default the page loads its two libraries
(deck.gl and umap-js) from a public CDN; `prepare.py vendor` puts them on your own site
instead.

SODA reads multiplexed immunofluorescence (OME-TIFF, QPTIFF, Akoya), H&E, Visium and Visium
HD, Xenium, CosMx, MERSCOPE, Stereo-seq, GeoMx DSP, single-cell RNA-seq, and the standard
formats: spatial AnnData, OME-Zarr and SpatialData.

## Features

- **Three tabs, one job each.** The controls sit beside the picture in three tabs: **View**
  (image and stains, genes, molecules, outlines), **Analyze** (clusters and selection, draw a
  region, volcano plot, heatmap and violin) and **Share** (save image, save / load, a link,
  dataset details, about and help). A tab shows only when the dataset has something for it.
  The tools few people need are under **More options** at the end of their tab (View: Colour
  by a value, Measure; Analyze: Your own results, Cell sets, Quality filter). Sections fold,
  and a folded one says its state on its header ("3 stains on", "Leiden 0.6", "12% removed");
  rarely used settings are under a section's own **More options** and each section's
  explanation is behind its **?**. The dataset picker and the Tissue / Side by
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
  once, fetched one gene at a time. Zoomed in, each cell is drawn with its own outline (the
  instrument's cell and nucleus boundaries, also for CosMx, MERSCOPE and Stereo-seq cell bins,
  or a lab's own segmentation file), filled with the colour of the gene or cluster on show;
  zoomed out, a dot at its centre. A zipped Xenium output is read as it is (no unzipping), and
  an output without its stain images opens on a transcript-density background.
- **GeoMx DSP.** A study of hand-drawn regions (ROIs) rather than cells or bins: the scans
  (several slides side by side), and each ROI's outline filled with the colour of a gene or
  of its group (a disc per segment where an ROI has several). The columns of the annotation
  sheet are the groups: legend, heatmap / dot plot / violin and the volcano compare groups of
  ROIs, and the QC panel has the GeoMx values (reads, saturation, nuclei, area, negative-probe
  background, LOQ). What works for ROI data is the comparing of groups; a UMAP, clusters,
  cell-type suggestions and density maps are not made for a few dozen regions, and the viewer
  says so instead of showing them. Values are **Q3-normalised** by default, as GeomxTools and
  the GeoMx software do (factors worked out at export; "Normalise" in the plots and volcano
  switches to counts per 10,000; older exports stay on counts per 10,000). Untested: several
  segments per ROI on real data; and the statistics treat ROIs as independent samples, which
  ROIs of one slide or animal are not.
- **CosMx SMI and MERSCOPE.** The same three parts as Xenium, from the vendors' flat files:
  the per-FOV CosMx images stitched into one mosaic (raw stains or the rendered composite),
  or one z-plane of the MERSCOPE mosaic stains; the cells as an expression layer; and
  every transcript in the Molecules panel. Control probes and blanks are left out, and a
  dataset that comes without images still opens on a transcript-density background.
- **Stereo-seq (STOmics).** The registered H&E or ssDNA image with square-bin layers summed
  from the DNB counts of a GEF or GEM file (bin20 / bin50 / bin100 = 10 / 25 / 50 µm, cut to
  the tissue) and the segmented cells of the cell bin file as a layer of their own; gene maps,
  UMAP and clusters, QC, Draw a region and the plots work on both.
- **Standard formats.** A spatial AnnData (`.h5ad` or AnnData-Zarr, the Squidpy layout), an
  OME-Zarr image (NGFF 0.4 or 0.5) or a SpatialData store (tables on circles or polygons, images,
  points, coordinate transformations), in zarr v2 or v3, is read without installing anndata or
  spatialdata, and becomes an ordinary SODA dataset: the image, the expression layer, cell
  outlines, molecules, and the file's own clusters, colours and UMAP. What is not read (a zarr v3
  codec outside the common ones, labels regions, plates, transformations other than scale /
  translation / affine) is refused with the reason. The page reads the converted export; opening these files in the browser
  directly is a possible later step.
- **Tissue | Side by side | UMAP.** Each section's bins are embedded once into a UMAP
  with Leiden clusters (three resolutions, marker genes per cluster). Switch the main view
  between the tissue and the UMAP, or show both linked side by side: select bins in one
  (box or freehand) and the same bins light up in the other. What is selected is named above
  the picture with **Show on slide**, which takes any selection, picked cluster or cell type
  back to the tissue and frames it there, and **Clear selection** (or Esc).
- **Suggested cluster names.** Each cluster gets a suggested cell type from editable
  marker lists, with a confidence and the runner-up shown; accept it, accept all
  high-confidence ones, or type your own. General lists for any tissue (genes and protein
  panels, with common marker-name aliases) are used by default; peripheral nerve lists are
  included too.
- **Draw a region (circle → UMAP).** Draw circles or freehand shapes on the tissue and get a UMAP, with
  optional Leiden clusters, of just the bins inside, computed in the browser. Drag a box or
  a loop over dots in its UMAP (Shift or **+** adds more) to light those bins on the slide,
  and **Show on slide** (above the picture) to frame them. Save the shapes and every parameter to redo it later.
- **Quality filter (per-bin QC).** Histograms of total counts, genes detected and % mitochondrial per bin,
  each with min / max thresholds. Off until you move a slider; removed bins are greyed on the
  slide and left out of Draw a region, and the thresholds are saved with the shapes and in
  the page link.
- **Save and share.** **Copy link to this view** (Share tab) copies an address that opens the
  dataset as it looks now: position and zoom, the Tissue / Side by side / UMAP view, the
  layer, the genes on show with their colours and ranges, stains on or off with their colours
  and display ranges, molecule genes, outlines, the colour-by, the clustering with your
  cluster names, cell sets and hidden clusters, the quality filter and the open tab. Drawn shapes are in
  it when they fit: a link is kept to about 2,000 characters, so long outlines are simplified
  or left out, and the message under the button says so. The view travels in the part of the
  address after `#`, which the browser never sends to a server; `ds=`, `qc=` and a `data=` or
  `static=` already in the address stay as they are. Not in a link: a box or freehand
  selection, volcano and plot settings, and the result
  of a Draw a region run (the shapes come back; press Run UMAP). **Save JSON** writes the same
  view to a file together with the shapes in full, every setting and the run's result, and
  **Load JSON** on the same dataset puts it back and runs it again.
- **Heatmap, dot plot and violin plots.** Compare genes across groups: the section's clusters
  at any resolution, the shapes you drew, the clusters of a Draw a region run, the regions of an overlay,
  or (single-cell datasets) any annotation column. Genes are each group's top markers plus any
  you type. A bin's value is log1p(counts per 10k); the heatmap shows the group mean (as it
  is, or z-scored / scaled 0–1 per gene), the dot plot adds the share of bins expressing, and
  violins show the spread with median and quartiles, with an "expressing only" switch because
  most bins hold zero. Every group shows its number of bins and the QC filter applies. Click a
  gene or a group to **show it on the slide**; save PNG, SVG, or the plotted numbers as CSV.
  Computed in your browser from the per-gene data (a typical cluster heatmap of a Visium HD
  section reads roughly 10–30 MB).
- **Multiplexed protein images (Akoya / PhenoCycler and other multi-channel TIFFs).** The
  slide is cut into 16 µm squares holding each marker's mean intensity (from full-resolution
  pixels), and the squares get everything above: UMAP, Leiden clusters, suggested names from
  protein marker lists, and Draw a region. Each marker is normalised per slide (background
  off, scaled to its 99.5th percentile, asinh), empty and duplicate channels are dropped
  automatically, and the export's slide QC (`qc.json`) lists weak and saturated stains, folds
  and bright edges, and possible tile seams. No cell segmentation (yet).
- **Save as image.** **Save image** (Share tab, or the camera beside Fit) saves the map as a
  PNG exactly as shown in Tissue, UMAP or Side by side, with a legend, the scale bar and the
  dataset's name, at 1x, 2x or 4x the screen's pixels (finer image tiles at 2x / 4x), on a
  black, white or transparent background. The volcano plot saves as PNG or SVG, like the plots.
- **Region overlays.** GeoJSON regions (for example from QuPath) drawn over the image;
  click a region to restrict the gene density to its bins.
- **Single-cell RNA-seq.** UMAP from AnnData (`.h5ad`) or 10x output, coloured by cluster
  with automatically placed labels, plus up to twelve gene expression panels.
- **Host it yourself.** Static files work on any web host. Optional password protection
  for private data on Azure, or a small live tile server for a lab network.
- **De-identification tools.** Sample IDs are replaced by numbers through a private key,
  and a check command gates uploads. They know only our lab's ID formats: edit `_ID` in
  `deid.py` for yours first (see [DEIDENTIFICATION.md](docs/DEIDENTIFICATION.md)). Barcode and
  cell id lists are in an export by default and are checked for sample IDs; `--no-ids` leaves them out.

<!-- SCREENSHOT PLACEHOLDER: docs/img/side_by_side.png
     A public example in Side by side mode: tissue on the left, its UMAP on the right, both
     coloured by Leiden cluster. Lead adds the image. -->

<!-- SCREENSHOT PLACEHOLDER: docs/img/circle_umap.png
     Draw a region (circle -> UMAP): two circles on the tissue and the resulting UMAP panel. -->

## Try it in 5 minutes

This installs SODA and opens the public examples from the page in your own copy. You need
Python 3.12 or 3.13 (64-bit), Git, and an internet connection.

**1. Get the code and install it with one command.**

macOS / Linux:

```bash
git clone https://github.com/FIVE-MGI/SODA.git soda
cd soda
bash install.sh
source soda-env/bin/activate
```

Windows (PowerShell or Command Prompt):

```powershell
git clone https://github.com/FIVE-MGI/SODA.git soda
cd soda
powershell -ExecutionPolicy Bypass -File install.ps1
soda-env\Scripts\Activate.ps1
```

The install script finds Python 3.12 or 3.13, makes a virtual environment `soda-env` in the
repository, installs the exact, tested package versions from
`web/SpatialViewer/requirements-lock.txt`, checks that they all import, and ends with
**SODA is installed**. It needs no administrator rights, takes a few minutes the first time
and is safe to run again. If no suitable Python is found, it says so and what to install.
(`-ExecutionPolicy Bypass` only lets this one script run; it changes no setting. If
PowerShell refuses the activation line, use Command Prompt and run
`soda-env\Scripts\activate.bat`.) Once the environment is active, `python` means the
environment's Python on every system.

> **By hand instead** (or where the script stops: Intel Macs, Linux on ARM, Python 3.10 or
> 3.11): make and activate a venv yourself, then
> `pip install -r web/SpatialViewer/requirements.txt -r web/SpatialViewer/requirements-clustering.txt`.
> That takes the newest versions the ranges allow, not the tested set; with Python 3.12 or
> 3.13, `pip install -r web/SpatialViewer/requirements-lock.txt` gives the tested set.

**2. Serve the page** from your copy, in a terminal with the environment active:

```bash
python -m http.server 8766 --directory web/SpatialViewer
```

**3. Open** this address in Chrome:

```
http://localhost:8766/index.html?data=https://fivestore.blob.core.windows.net/soda/examples/
```

Pick a dataset, search a gene, then try the **Tissue | Side by side | UMAP** switch. Press
Ctrl+C in the terminal to stop. To open your own export the same way, see
[HOSTING.md, On your own computer](docs/HOSTING.md#a-on-your-own-computer).

To test the code: `python tests/smoke_test.py`. What it runs and what it needs is in
[REPRODUCE.md, section 6](REPRODUCE.md#6-the-smoke-test).

## Use your own data

1. **[Prepare your data](docs/PREPARING_DATA.md):** point SODA at a dataset's folder and it works
   out what it is (Visium HD, Visium, Xenium, CosMx, MERSCOPE, Stereo-seq, GeoMx DSP, OME-TIFF / QPTIFF,
   per-marker TIFF folders, H&E, `.h5ad` or 10x RNA-seq), and the standard formats: spatial AnnData
   (`.h5ad` or AnnData-Zarr, as scanpy / Squidpy write it), OME-Zarr images and SpatialData stores:

   ```bash
   python web/SpatialViewer/prepare.py detect path/to/folder                                # what is in there? (writes nothing)
   python web/SpatialViewer/prepare.py auto --src path/to/folder --id my_dataset --out cache   # prepare it
   python web/SpatialViewer/prepare.py export --id my_dataset --out cache --dest export
   ```

   When it is not certain it does not guess: it lists what it found and the `--kind ...` to add.
   Every kind also has its own command with all its options, and
   `notebooks/SODA_pipeline.ipynb` does the same for a folder of datasets.

   **Barcode and cell id lists are included by default.** The export holds each layer's barcodes or
   cell ids (`expr/<layer>/ids.txt.gz`), as Loupe Browser, Xenium Explorer and Vitessce keep them, so
   people can save the bins inside a shape they draw with their barcodes. Every list is checked for
   sample IDs first (and, for the instruments SODA knows, against the instrument's id format); a list
   that fails is left out with a message, never printing the ids. The check catches the usual
   mistakes: **you, the data owner, are responsible for what you publish.** To leave every list out,
   add `--no-ids` to `prepare.py export` (or `prepare.py batch`, `SODA prepare`), set
   `INCLUDE_IDS = False` in the notebook, or untick "Include barcode / cell id lists" in the desktop
   app ([ADVANCED.md](docs/ADVANCED.md#barcode-and-cell-id-lists)).
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
- **Memory (RAM):** no minimum has been benchmarked. Measured peaks: a CosMx slide's counts (26
  million transcripts, without its images) 0.76 GB; a whole mouse brain Stereo-seq chip
  2.7 GB to prepare (3.9 GB from `tissue.gem.gz`) and **8.5 GB to export, the largest measured so far**. For whole-section datasets, have more
  free memory than that. Visium HD is prepared with its 8 µm and 16 µm layers by default;
  its 2 µm layer (`--bins all`) has not been measured and is expected to need the most
  (see [How big and how long](docs/PREPARING_DATA.md#how-big-and-how-long)). The notebook
  compares your free memory with these numbers before the Prepare step and warns.
- Disk space: the prepare cache holds a compressed full-resolution copy of each image, and
  the export a tiled copy. Measured times for one dataset of each kind, with the machine and
  conditions, are in [docs/PERFORMANCE.md](docs/PERFORMANCE.md); sizes are in
  [How big and how long](docs/PREPARING_DATA.md#how-big-and-how-long).
- Viewing needs no installation, only a browser. The page loads its drawing library (deck.gl),
  and Draw a region its UMAP library (umap-js), from `cdn.jsdelivr.net`, or from your own copy
  after `python web/SpatialViewer/prepare.py vendor`.

## Browser support

The viewer needs **WebGL2**. It is developed and tested in **Google Chrome** on desktop.
Other desktop browsers (Edge, Firefox, Safari) have not been tested by hand; they may work if they
support WebGL2. Automated tests (`tests/browsers_test.js`, `tests/offline_test.js`) run the page in
Playwright's Chromium, WebKit, Edge and Firefox builds, Firefox also in a phone-sized window (the size
only: Playwright cannot emulate a touch screen in Firefox). On an iPhone every browser runs on
Safari's engine (WebKit); what was tested there is below.

**Copy link to this view** and **Compare** were checked with automated tests in Chromium
and in Playwright's WebKit build (the engine behind Safari), at desktop size and at
390 × 844; they have not been tried in Safari itself, in Firefox, or on a real phone.
A link made in one browser opens in the other. Copy link to this view (not Compare) was
also checked in Playwright's Firefox, at desktop size.

**Phones and tablets** open the same page and link; the layout follows the screen.

- On a phone held upright the map fills the screen and the three tabs are a bar along the
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
install.sh / install.ps1  one-command install (macOS, Linux / Windows) into soda-env
web/SpatialViewer/    the page (index.html, viewer.js, ...) and the Python tools (prepare.py, export.py, ...)
web/SpatialViewer/requirements-lock.txt  every package at one tested version (made from requirements*.txt)
web/SpatialViewer/markers/  marker-gene lists for suggested cluster names (edit for your tissue)
web/SpatialViewer/stdformats.py  spatial AnnData, OME-Zarr and SpatialData into a SODA dataset (with annstore.py, ngff.py)
web/SpatialViewer/detect.py  says what kind of dataset a folder holds (prepare.py detect / auto, the notebook)
deploy/azure/         Azure Static Web Apps config and the /api/token function
docs/                 documentation
tests/smoke_test.py   end-to-end test
tests/stdformats_browser_test.js  the standard-format exports open in Chromium (needs Playwright; the smoke test runs it)
tests/phone_browser_test.js  the image path in phone browsers, with each failure injected (needs Playwright)
tests/diag_test.js    the ?diag=1 panel never prints tokens (node only; the smoke test runs it)
tests/viewlink_test.js  the link that restores a view: round trips, the length budget, damaged and hostile links (node only; the smoke test runs it)
tools/label_editor.html  place UMAP cluster labels by hand (writes labelPositions for config.js)
tools/check_install.py   checks that the Python packages import, with the right versions
tools/perf/           the timing scripts behind docs/PERFORMANCE.md
REPRODUCE.md          how to reproduce the paper's numbers and figures
CITATION.cff          how to cite SODA (GitHub's "Cite this repository")
CHANGELOG.md          what is in each version
CONTRIBUTING.md       reporting a problem, asking a question, sending a change
.github/              the automated test (smoke.yml), app builds, Pages, issue templates
```

## Citing

Manuscript in preparation. Until it is published, please cite this repository: GitHub's
**Cite this repository** button reads [CITATION.cff](CITATION.cff) (version 0.9.0; no DOI yet).
<!-- CITATION PLACEHOLDER: replace with the paper reference and DOI once available. -->

## Licence

SODA is free software: you can redistribute it and/or modify it under the terms of the
GNU Affero General Public License as published by the Free Software Foundation, either
version 3 of the License, or (at your option) any later version. See [LICENSE](LICENSE).
Copyright (C) 2026 Washington University.

In short: you may use, change and share SODA, including for commercial work. If you share
a changed version, or run a changed version as a website that other people use, you must
make your changed source code available under the same licence. The licence covers SODA's
code, including the viewer page copied into each export; it does not cover your data.
