# SODA: Spatial Omics Data Annex

## [Download the SODA app for Windows or Mac](https://github.com/FIVE-MGI/SODA/releases/latest)

SODA turns what a spatial biology instrument writes into a website that opens the whole slide
in a web browser. You prepare your data once on your own computer with the SODA app; you upload
the result to your own cloud storage; your own free copy of SODA on GitHub shows it. Anyone
with the link can look, with nothing to install.

**Licence:** GNU Affero General Public License, version 3 or later (AGPL-3.0-or-later), the same licence as SiDoLa-NS. See [Licence](#licence).

<!-- SCREENSHOT PLACEHOLDER: docs/img/overview.png
     Whole-slide multiplex IF view with the channel panel open. Lead adds the image. -->

## Quick start

**See it first.** Four public datasets (Xenium, Visium HD, PhenoCycler and MERSCOPE), prepared
with SODA, open in Chrome with nothing to install: **https://sidolans01.mgifive.org/SODA-examples/**

Pick a dataset, search a gene, then try **Tissue | Side by side | UMAP**. Where each one comes
from: [REPRODUCE.md, section 9](REPRODUCE.md#9-the-public-datasets).

**Your own data, no code.** Nothing leaves your computer until step 3.

1. **Download the app** from the [latest release](https://github.com/FIVE-MGI/SODA/releases/latest):
   `SODA-<version>-setup.exe` (Windows 10 or 11) or `SODA-<version>-macos-arm64.dmg` (Apple
   Silicon Mac, macOS 14 or newer). The app is not signed yet: the first time, Windows says
   "Windows protected your PC" (click **More info**, then **Run anyway**); on a Mac, follow the
   README inside the .dmg.
2. **Prepare.** Put each dataset in its own folder inside one data folder (what each instrument
   needs: [PREPARING_DATA.md](docs/PREPARING_DATA.md)). Open SODA, pick the data folder and an
   empty output folder, click **Prepare**. A green banner says the export is clean.
3. **Upload** the `export` folder (inside the output folder) to your own cloud storage (Azure,
   AWS or Google Cloud; not Box, Google Drive or Dropbox), readable by anyone, with a CORS rule
   that lets your site read it.
4. **Make your site.** On this page click **Use this template** to make your own copy of SODA,
   turn on GitHub Pages, and paste your storage folder's address into `web/SpatialViewer/config.js`.
   Your site is `https://<your GitHub name>.github.io/<your copy>/SODA/`.

Every click for steps 1-2: [GETTING_STARTED.md](docs/GETTING_STARTED.md); for steps 3-4:
[SHARE_YOUR_DATA.md](docs/SHARE_YOUR_DATA.md). More about the app: [DESKTOP_APP.md](docs/DESKTOP_APP.md).

**More control: the notebook.** `notebooks/SODA_pipeline.ipynb` runs the same code as the app,
with every setting ([GETTING_STARTED.md, Part 1](docs/GETTING_STARTED.md#part-1-or-the-notebook-install-once-about-15-minutes)).
Installing from the code, the command line and a whole folder of datasets in one command:
[REFERENCE.md](docs/REFERENCE.md).

## Documentation

| Guide | What it covers |
|---|---|
| [docs/GETTING_STARTED.md](docs/GETTING_STARTED.md) | the first time, step by step: preparing your data with the app or the notebook |
| [docs/DESKTOP_APP.md](docs/DESKTOP_APP.md) | the SODA desktop app for Windows and Mac (no Python needed) |
| [docs/REFERENCE.md](docs/REFERENCE.md) | what the viewer does, installing from the code, the command line, requirements, browsers and phones, the repository |
| [docs/PREPARING_DATA.md](docs/PREPARING_DATA.md) | every input type, its command and options; export, clustering, Draw a region; common problems |
| [docs/SHARE_YOUR_DATA.md](docs/SHARE_YOUR_DATA.md) | your own SODA website: upload to your cloud storage, copy SODA on GitHub, connect your data |
| [docs/HOSTING.md](docs/HOSTING.md) | your own computer or web server; what the storage or host must do; pre-upload checklist |
| [deploy/azure/README.md](deploy/azure/README.md) | Azure, step by step: private container, password, Static Web App |
| [docs/ADVANCED.md](docs/ADVANCED.md) | other places for the data and the page, passwords, names and colours, every notebook setting, your own clusters, marker lists, the code, batch runs |
| [docs/DEIDENTIFICATION.md](docs/DEIDENTIFICATION.md) | `deid.py`, the private key file, the sample ID check before upload, and its limits |
| [docs/DATA_FORMAT.md](docs/DATA_FORMAT.md) | the exported files, byte by byte, for people writing their own tools |
| [docs/PERFORMANCE.md](docs/PERFORMANCE.md) | measured times and memory, the machine they came from, and how to measure your own |
| [REPRODUCE.md](REPRODUCE.md) | how to reproduce the numbers and figures in the paper |
| [CHANGELOG.md](CHANGELOG.md) | what is in each version |
| [CONTRIBUTING.md](CONTRIBUTING.md) | reporting a problem, asking a question, sending a change |

## Troubleshooting

The common failures, and where the fix is written down.

- **The page says "No data connected yet".** Your copy's `web/SpatialViewer/config.js` does
  not name your storage folder yet: [SHARE_YOUR_DATA.md, Step 7](docs/SHARE_YOUR_DATA.md#step-7-connect-your-data).
- **The page says it "could not reach" your folder, or genes do not load.** Your storage
  does not yet let your site read it. Set the CORS rule once, and allow the `Range`
  header: [SHARE_YOUR_DATA.md, Step 4](docs/SHARE_YOUR_DATA.md#step-4-let-your-site-read-your-storage-cors-once).
  The notebook's last cell tests this and prints the fix for your storage. A host that does
  not answer Range requests at all: [HOSTING.md, What the host must do](docs/HOSTING.md#what-the-host-must-do).
- **A Box, Google Drive, Dropbox or OneDrive link does not work.** Their share links open a
  preview page, not the files, and other sites may not read them. Use Azure, AWS or Google
  Cloud: [SHARE_YOUR_DATA.md, Step 2](docs/SHARE_YOUR_DATA.md#step-2-make-your-cloud-storage-once).
- **A `?data=` link opens nothing on a password-protected site.** A site with a password
  (`tokenEndpoint` in `config.js`) does not open `?data=` links from outside, unless
  `config.js` sets `allowDataLinks: true`: [ADVANCED.md, section 2](docs/ADVANCED.md#2-private-data-with-a-password).
- **The install stops, or says no suitable Python was found.** The install script needs
  Python 3.12 or 3.13, 64-bit. On Intel Macs, Linux on ARM, or with Python 3.10 or 3.11, use
  the manual install in [REFERENCE.md, Try it in 5 minutes](docs/REFERENCE.md#try-it-in-5-minutes). Python 3.14 is not
  covered yet ([Requirements](docs/REFERENCE.md#requirements)). The SODA app needs no Python. `python tools/check_install.py` says which
  package is missing or has the wrong version.
- **The image does not appear (often on a phone).** The status line says why; add `&diag=1`
  to the address for a full report: [Browser support](docs/REFERENCE.md#browser-support).
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

## Citing

Manuscript in preparation. Until it is published, please cite this repository: GitHub's
**Cite this repository** button reads [CITATION.cff](CITATION.cff) (version 0.10.0; no DOI yet).
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
