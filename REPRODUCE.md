# Reproducing the paper's numbers and figures

This page says how to rebuild what the SODA paper reports, with the version of SODA it used
(0.9.0, see [CHANGELOG.md](CHANGELOG.md)). Every command below is in this repository, and
every dataset in the paper is public: you need the repository and an internet connection.

Three kinds of result, three kinds of promise:

| Kind | Example | What you should get |
|---|---|---|
| The tests | the smoke test and the browser tests, on the public example exports | the same checks passing, with the same package versions |
| Times and memory | Supplementary Table 2, the tables in [docs/PERFORMANCE.md](docs/PERFORMANCE.md) | the same order of magnitude; they depend on the machine and on the data |
| Results on real datasets | Supplementary Table 2, Figures 2 and 3 | the same commands on the same public datasets (sections [9](#9-the-public-datasets) to [11](#11-the-figures)) |

## 1. Get the version the paper used

```bash
git clone https://github.com/FIVE-MGI/SODA.git soda
cd soda
git checkout v0.9.0
```

<!-- RELEASE: create the tag v0.9.0 on the release commit before the paper posts. -->

## 2. The environment

**Python 3.12 or 3.13, 64-bit.** The install script makes `.venv` in the repository and
installs every package at one exact version from the lock file:

```bash
bash install.sh                                        # macOS / Linux
powershell -ExecutionPolicy Bypass -File install.ps1   # Windows
```

- To use one particular Python: `bash install.sh --python /path/to/python3.13`, or
  `install.ps1 -Python C:\path\to\python.exe`.
- By hand, the same set: make and activate a venv with Python 3.12 or 3.13, then
  `pip install -r web/SpatialViewer/requirements-lock.txt`.
- Check it: `python tools/check_install.py` imports every package, prints its version and
  checks the version rules the code depends on (zarr 2.x, tifffile older than 2025, Python
  3.12 or 3.13). Exit code 0 means everything is fine.

Where it was run:

| System | Python | How |
|---|---|---|
| Windows 11 x64 | 3.13 | the lock file was made here (pip freeze of a clean 3.13.9 venv) and the smoke test passes |
| Ubuntu (GitHub Actions) | 3.12 | `.github/workflows/smoke.yml` installs with `install.sh` and runs the smoke test on every push |
| macOS 12+ (Apple Silicon), Linux x86_64 | 3.12, 3.13 | every pinned version has a ready-made wheel (checked with a resolver); not run by us on a Mac |

Not covered by the lock file: Python 3.11 or older, Python 3.14, Intel Macs, Linux on ARM.
There the manual install (`requirements.txt` + `requirements-clustering.txt`) gets the newest
versions the ranges allow, which may give different numbers.

**Node.js 18 or newer** for the JavaScript tests (the automated test uses Node 22).

## 3. Where the versions are pinned

| What | Where |
|---|---|
| Every Python package, one exact version | `web/SpatialViewer/requirements-lock.txt` (do not edit by hand; its header says how to refresh it) |
| The ranges the lock file is made from | `web/SpatialViewer/requirements.txt` (core), `web/SpatialViewer/requirements-clustering.txt` (UMAP and Leiden) |
| The two browser libraries (deck.gl 9.1.12, umap-js 1.4.0), with SHA-384 hashes | `web/SpatialViewer/libs.js` (also in `index.html`'s `vendor/` tag and used by `vendor.py`) |
| Python and Node of the automated test | `.github/workflows/smoke.yml` (Python 3.12, Node 22) |
| The desktop app builds | `.github/workflows/build-windows.yml`, `build-macos.yml` (Python 3.13, PyInstaller 6.21.0) |
| The analysis settings (genes, PCs, neighbours, resolutions, seed) | constants at the top of `web/SpatialViewer/sectionumap.py` |

Each export records what made its clustering: `expr/<layer>/umap/index.json` holds the
settings (`params`, including the seed) and the versions of umap-learn, leidenalg, igraph and
numpy that ran ([DATA_FORMAT.md](docs/DATA_FORMAT.md#umap-whole-section-umap-and-leiden-clusters)).

## 4. Seeds

Everything random is seeded. The default seed is **42**.

| Where | Seed |
|---|---|
| Whole-section UMAP and Leiden at export (`sectionumap.py`: umap-learn, leidenalg, the sample for large layers); the protein squares use the same code | `SEED = 42` |
| The shuffled order each gene's molecules are written in (`molecules.py`) | `SEED = 42` |
| Draw a region in the browser (`lasso_worker.js`, `leiden.js`) | 42; Save JSON stores the seed with the shapes, so Load JSON reruns the same |
| Volcano plot, heatmap / violin sampling (`volcano.js`, `de_worker.js`, `plots.js`, `plots_worker.js`) | 42 |

**Threads do not change the result.** `prepare.py` caps the maths libraries' threads (at
most 8; `--threads <n>` or `SODA_THREADS=<n>` to change it). The thread count changes only the
last bits of some floating-point values; the exported UMAP coordinates and clusters came out
byte-identical with 1, 8 and 24 threads on the datasets compared. The smoke test checks this
(see `web/SpatialViewer/threadcap.py`).

## 5. The public example exports

Four of the datasets of section [9](#9-the-public-datasets) (`xenium_pub`, `visiumhd_pub`,
`phenocycler_pub`, `merscope_2`), prepared and exported with SODA, are online in one folder
with its `datasets.json`:

```
https://fivestore.blob.core.windows.net/soda/examples/
```

Open them on the shared page, or from your own copy as in the README's
[Try it in 5 minutes](README.md#try-it-in-5-minutes), with
`index.html?data=https://fivestore.blob.core.windows.net/soda/examples/`. The tests open these
exports (sections [6](#6-the-smoke-test) and [7](#7-the-browser-tests)).

The other three (`cosmx_1`, `stereoseq_1`, `geomx_1`) are not redistributed, for licence
reasons (see the table in section 9): download them from their sources with the commands in
section [9](#9-the-public-datasets).

## 6. The smoke test

```bash
python tests/smoke_test.py
```

From the repository folder, with the environment active. It checks the tools and the public
example exports of section [5](#5-the-public-example-exports), read by default from
https://fivestore.blob.core.windows.net/soda/examples/: the tiles, the gene data, the sample
ID check, and much more (the full list is at the top of `tests/smoke_test.py`). To use a
local copy of the exports instead, pass `--exports <folder>` (or set `SODA_TEST_DATA` to the
folder). With `--raw-data <folder>` (or `SODA_RAW_DATA`), the raw files laid out as in section
[9](#9-the-public-datasets), it also runs the pipeline (prepare, export, ID check) on them;
without it that part prints SKIPPED. It ends with
**SMOKE TEST PASSED: N checks**.

- `--root <dir>` sets the folder it writes to (default: `soda-smoke` in the temp folder);
  `--reuse` keeps earlier outputs; `--checks-only` runs no commands and only checks the
  outputs already there.
- The log is `smoke_test.log`, the summary `smoke_test_summary.json`, in that folder.
- Exit code 0: every check passed. 1: a check failed. 2: a command failed.
- Without the clustering packages the UMAP checks are skipped, and it says so.
- Without Node.js the JavaScript tests are skipped. Without Playwright the browser tests print
  **SKIPPED** and still pass: for a complete run, install Playwright (below) and make sure no
  line says SKIPPED.

## 7. The browser tests

Playwright is not a dependency of this repository. Install it once, in the repository folder
(it goes into `node_modules/`, which Git ignores), with its browsers:

```bash
npm install --no-save playwright@1.56
npx playwright install chromium webkit firefox
```

The performance numbers used Playwright 1.56. The tests also find Playwright through
`NODE_PATH` or `--playwright <folder holding node_modules>`, and its browsers through
`PLAYWRIGHT_BROWSERS_PATH`. Microsoft Edge is tested only where it is installed; a browser
that cannot be started is reported as SKIPPED on its own.

The browser tests open the public example exports of section
[5](#5-the-public-example-exports), at https://fivestore.blob.core.windows.net/soda/examples/;
to use a local copy, pass `--data-root <folder>` (or set `SODA_TEST_DATA`). The three
restricted datasets (`cosmx_1`, `geomx_1`, `stereoseq_1`) are tested only when
`SODA_TEST_PRIVATE_DATA` (or `--private-data-root`) points at a copy of them; otherwise
those cases print SKIPPED. The smoke test already runs most of
them. To run one by itself:

| Command | What it checks |
|---|---|
| `node tests/browsers_test.js` | the page on every example dataset in Chromium, WebKit (Safari's engine), Firefox and Edge, plus phone-sized windows; pixel comparison with Chromium |
| `node tests/phone_browser_test.js` | the image path in phone browsers (WebKit as an iPhone, Chromium), with each failure injected. Not run by the smoke test |
| `node tests/offline_test.js` | the page and a kept dataset open with no network |
| `node tests/outlines_browser_test.js` | cell outlines, and what draws a subset of the cells with them |
| `node tests/ownres_browser_test.js` | loading your own results (CSV / TSV) into the page |
| `node tests/stdformats_browser_test.js` | the standard-format exports open |
| `node tests/compat_test.js` | the page opens exports made by older versions |
| `node tests/viewlink_test.js`, `node tests/diag_test.js` | the view link and the `?diag=1` panel (Node only, no browser) |
| `python tests/app_test.py` | the desktop app, from source (`--browser` also opens its pages in Chromium). Not run by the smoke test |

Each file's header lists its options (`--only`, `--shots <dir>`, `--export <dir>` ...).
Playwright's WebKit is the WebKit engine, not Safari on an iPhone or a Mac; see
[Browser support](README.md#browser-support) for what was tested on real devices.

## 8. Times and memory

The numbers, the machine and the conditions are in [docs/PERFORMANCE.md](docs/PERFORMANCE.md).
In short: an Intel Core i7-13700KF using **only its 8 efficiency cores**, 32 GB RAM, Windows 11,
Python 3.13, headless Chromium (Playwright) on the real GPU, data on network shares, the page
served by a local server that answers like blob storage, and Chrome's network emulation for
the "home connection" (50 Mbit/s down, 10 up, 20 ms). Other jobs sometimes ran on the same
cores, so treat the times as an order of magnitude. The times in the paper were measured on
commit 57acd2f (18 November 2026); section [10](#10-supplementary-table-2-times-and-memory)
has the command for every row.

The scripts are in `tools/perf/`. `bench_prepare.py` needs `psutil`
(`pip install psutil`; it is not in the lock file). `bench_open.js` needs Playwright, as above.

```bash
# preparing or exporting: time, peak memory, size
python tools/perf/bench_prepare.py --label my_ds --step export --csv perf.csv --logdir perf_logs \
    --measure export/my_ds -- python web/SpatialViewer/prepare.py export --id my_ds --out cache --dest export

# opening a dataset, cold and warm, unthrottled and on a home connection
node tools/perf/bench_open.js --data-root export --ds my_ds --reps 3 --genes top,median --csv perf.csv --out perf_logs
node tools/perf/bench_open.js --data-root export --ds my_ds --reps 3 --net home --csv perf.csv --out perf_logs

# the tables
python tools/perf/summarize.py perf.csv
```

- **The same cores on Windows:** start the command under the same affinity, for example
  `cmd /c "start /affinity FF0000 /wait /b python tools/perf/bench_prepare.py ..."` (logical
  CPUs 16-23, the efficiency cores of that CPU; another CPU needs its own mask; on Linux,
  `taskset -c 16-23 <command>`). Child processes, Chromium included, inherit it. The affinity
  the command had is recorded in the output.
- **Draw a region times:** `bench_open.js --mode region` with `--layer`, `--regions <circles.json>`
  and `--leiden off,on` (see the top of `tools/perf/bench_open.js`).
- Each script's header lists every option and what each number means.

## 9. The public datasets

Every real dataset in the paper is a public download. None of them is in this repository.

| Id | What | Source, licence | Download |
|---|---|---|---|
| `xenium_pub` | 10x Genomics Xenium, human breast cancer, replicate 1 (Janesick et al. 2023), as the SpatialData sandbox store (zarr v3, spatialdata 0.7.1) | SpatialData sandbox (s3.embl.de); 10x Genomics' terms for the original data | 3,616,589,814 B |
| `visiumhd_pub` | 10x Genomics Visium HD, human colorectal cancer (FFPE), Space Ranger 4.0.1 | 10x Genomics, CC BY 4.0 | 15,912,158,814 B (binned outputs) + 11,098,424,743 B (full-resolution H&E) |
| `phenocycler_pub` | HuBMAP HBM274.QXJB.759, human lung, PhenoCycler-Fusion QPTIFF, 46 channels, 36,480 × 37,440 px | HuBMAP, doi:10.35079/HBM274.QXJB.759, CC BY 4.0 | 14,154,986,909 B |
| `cosmx_1` | NanoString CosMx SMI, non-small-cell lung cancer (FFPE), Lung9_Rep1 (2021 release), 20 FOVs | NanoString public data share; no dataset licence stated (Bruker Spatial Biology terms of use, research use only) | 2.95 GB |
| `merscope_2` | Vizgen MERSCOPE, human fetal colon, sample AAN4, one whole region without images | Zenodo 19450420, CC BY 4.0 | 275,238,975 B |
| `stereoseq_1` | STOmics Stereo-seq demo chip C04042E3, mouse whole brain (fresh frozen, H&E) | STOmics demo data (https://www.stomics.tech/col1317); no open licence stated, the site's terms apply (do not redistribute) | 4,497,528,479 B (3 files) |
| `geomx_1` | NanoString GeoMx DSP Spatial Organ Atlas, mouse embryo E11 (2 scans, 84 AOIs, mouse whole transcriptome) | NanoString Spatial Organ Atlas; no licence stated (Bruker Spatial Biology terms of use, research use only) | 3,525,766,712 B (2 files) |
| `merscope_crop` (Figure 2) | Vizgen MERSCOPE small data sample: a cropped region with its images (500 genes plus 50 blank control probes, 1,179 cells) | Zenodo 7852005, CC BY 4.0 | 129,368,282 B |

The commands below put each dataset in its own folder under `data/`. `curl -C -` resumes a
download that stopped; compare each file's size with the table afterwards. They run in bash
and in the Windows command prompt (Windows 10 and later ship `curl` and `tar`).

**xenium_pub.** Download, then unzip; this gives `data/xenium_pub/data.zarr`:

```bash
curl -fL -C - --create-dirs -o data/xenium_pub/xenium_rep1_io_spatialdata_0.7.1.zip https://s3.embl.de/spatialdata/spatialdata-sandbox/xenium_rep1_io_spatialdata_0.7.1.zip
```

```bash
python -m zipfile -e data/xenium_pub/xenium_rep1_io_spatialdata_0.7.1.zip data/xenium_pub
```

**visiumhd_pub.** The binned outputs and the full-resolution H&E:

```bash
curl -fL -C - --create-dirs -o data/visiumhd_pub/Visium_HD_Human_Colon_Cancer_binned_outputs.tar.gz https://cf.10xgenomics.com/samples/spatial-exp/4.0.1/Visium_HD_Human_Colon_Cancer/Visium_HD_Human_Colon_Cancer_binned_outputs.tar.gz
```

```bash
curl -fL -C - --create-dirs -o data/visiumhd_pub/Visium_HD_Human_Colon_Cancer_tissue_image.btf https://cf.10xgenomics.com/samples/spatial-exp/4.0.1/Visium_HD_Human_Colon_Cancer/Visium_HD_Human_Colon_Cancer_tissue_image.btf
```

Extract only the 8 µm and 16 µm bins (the default layers); this gives
`data/visiumhd_pub/binned_outputs/square_008um` and `square_016um`:

```bash
tar -xzf data/visiumhd_pub/Visium_HD_Human_Colon_Cancer_binned_outputs.tar.gz -C data/visiumhd_pub binned_outputs/square_008um/filtered_feature_bc_matrix.h5 binned_outputs/square_008um/spatial binned_outputs/square_016um/filtered_feature_bc_matrix.h5 binned_outputs/square_016um/spatial
```

**phenocycler_pub.** The HuBMAP server refuses curl's default User-Agent; `-A` sends a client
name instead (any honest name works):

```bash
curl -fL -C - --create-dirs -A "SODA-data-fetch/0.9" -o data/phenocycler_pub/D438-RLL-16B4-17_Scan1.qptiff https://assets.hubmapconsortium.org/be9e06cea900e49292b5668a6f3e1789/raw/images/D438-RLL-16B4-17_Scan1.qptiff
```

**cosmx_1.** From NanoString's public bucket: the flat files, the three per-FOV image folders
(FOVs F001 to F020), and the one raw morphology image the timed folder also held (SODA uses the
CellComposite images unless every FOV has a raw one, so it does not change the export):

```bash
curl -fL --create-dirs -o "data/cosmx_1/#1" "https://nanostring-public-share.s3.us-west-2.amazonaws.com/SMI/Lung9_Rep1/{Lung9_Rep1_exprMat_file.csv,Lung9_Rep1_metadata_file.csv,Lung9_Rep1_fov_positions_file.csv,Lung9_Rep1_tx_file.csv,SMI-ReadMe.html}"
```

```bash
curl -fL --create-dirs -o "data/cosmx_1/CellComposite/CellComposite_F#1.jpg" "https://nanostring-public-share.s3.us-west-2.amazonaws.com/SMI/Lung9_Rep1/CellComposite/CellComposite_F[001-020].jpg"
```

```bash
curl -fL --create-dirs -o "data/cosmx_1/CellLabels/CellLabels_F#1.tif" "https://nanostring-public-share.s3.us-west-2.amazonaws.com/SMI/Lung9_Rep1/CellLabels/CellLabels_F[001-020].tif"
```

```bash
curl -fL --create-dirs -o "data/cosmx_1/CompartmentLabels/CompartmentLabels_F#1.tif" "https://nanostring-public-share.s3.us-west-2.amazonaws.com/SMI/Lung9_Rep1/CompartmentLabels/CompartmentLabels_F[001-020].tif"
```

```bash
curl -fL --create-dirs -o data/cosmx_1/ResegmentationImages/20210706_105344_S1_C902_P99_N99_F001_Z001.TIF https://nanostring-public-share.s3.us-west-2.amazonaws.com/SMI/Lung9_Rep1/ResegmentationImages/20210706_105344_S1_C902_P99_N99_F001_Z001.TIF
```

**merscope_2.** Download, then unzip; this gives `data/merscope_2/AAN4`:

```bash
curl -fL -C - --create-dirs -o data/merscope_2/AAN4.zip https://zenodo.org/records/19450420/files/AAN4.zip
```

```bash
python -m zipfile -e data/merscope_2/AAN4.zip data/merscope_2
```

**stereoseq_1.** The registered H&E, the bin GEF and the cell bin GEF:

```bash
curl -fL --create-dirs -o "data/stereoseq_1/#1" "https://demo.stomicsdb.tech/C04042E3_Mouse_Whole_Brain_Stereo-seq_FF_V1.3_HE/outs/{C04042E3_HE_regist.tif,C04042E3.raw.gef,C04042E3.cellbin.gef}"
```

**geomx_1.** The workflow-and-counts archive and the scan images, from the Spatial Organ Atlas
bucket:

```bash
curl -fL --create-dirs -o "data/geomx_1/#1" "https://external-soa-downloads-p-1.s3.us-west-2.amazonaws.com/{mu_dev_E11_workflow_and_counts_files.tar.gz,mu_dev_E11_image_files.tar.gz}"
```

```bash
tar -xzf data/geomx_1/mu_dev_E11_workflow_and_counts_files.tar.gz -C data/geomx_1
```

```bash
tar -xzf data/geomx_1/mu_dev_E11_image_files.tar.gz -C data/geomx_1
```

The timed run read the DCC files and the probe kit unzipped (`prepare.py geomx` also takes the
two .zip files as they are):

```bash
python -m zipfile -e data/geomx_1/mu_dev_E11_workflow_and_count_files/workflow/dcc/E9-E11-E15_DCC.zip data/geomx_1/mu_dev_E11_workflow_and_count_files/workflow/dcc/E9-E11-E15_DCC_unzipped
```

```bash
python -m zipfile -e data/geomx_1/mu_dev_E11_workflow_and_count_files/workflow/pkc/Mm_R_NGS_WTA_v1.0.zip data/geomx_1/mu_dev_E11_workflow_and_count_files/workflow/pkc/Mm_R_NGS_WTA_v1.0_unzipped
```

**merscope_crop (Figure 2).** Download, then unzip into a folder named `merscope_crop`. That
folder should hold `detected_transcripts.csv`, `images/` and `vpt-run/` directly; if they land
in a subfolder, move them up:

```bash
curl -fL -C - --create-dirs -o data/vizgen_sample.zip https://zenodo.org/records/7852005/files/vizgen_sample.zip
```

```bash
python -m zipfile -e data/vizgen_sample.zip data/merscope_crop
```

## 10. Supplementary Table 2: times and memory

Every timed command runs on the 8 efficiency cores only (logical CPUs 16-23 of the
i7-13700KF; section [8](#8-times-and-memory)), from the repository folder, with the environment
active. Run them with the machine otherwise quiet. Each prepare or export step runs under
`bench_prepare.py`, and the step itself is the part after `--`. The first one in full, on
Windows:

```bash
start "" /affinity FF0000 /wait /b python tools/perf/bench_prepare.py --label xenium_pub --step prepare --csv perf.csv --logdir perf_logs --input data/xenium_pub/data.zarr --measure cache/xenium_pub -- python web/SpatialViewer/prepare.py spatialdata --src data/xenium_pub/data.zarr --id xenium_pub --out cache --pixel-size 0.2125 --label "Breast cancer, Xenium (10x public data)"
```

Every other step is the same line with its own `--label`, `--step` (`prepare` or `export`),
`--input` (the files or folders the step reads; only recorded, for the input size) and
`--measure` (`cache/<id>` for a prepare, `export/<id>` for an export).

**Prepare, one real dataset per kind.** The `--label` text is the name the page shows:

```bash
python web/SpatialViewer/prepare.py spatialdata --src data/xenium_pub/data.zarr --id xenium_pub --out cache --pixel-size 0.2125 --label "Breast cancer, Xenium (10x public data)"
```

```bash
python web/SpatialViewer/prepare.py visium --id visiumhd_pub --image data/visiumhd_pub/Visium_HD_Human_Colon_Cancer_tissue_image.btf --outs data/visiumhd_pub --out cache --label "Colorectal cancer, Visium HD (10x public data)"
```

```bash
python web/SpatialViewer/prepare.py tiff --src data/phenocycler_pub/D438-RLL-16B4-17_Scan1.qptiff --id phenocycler_pub --out cache --label "Lung, PhenoCycler (HuBMAP public data)"
```

```bash
python web/SpatialViewer/prepare.py cosmx --src data/cosmx_1 --id cosmx_1 --out cache
```

```bash
python web/SpatialViewer/prepare.py merscope --src data/merscope_2/AAN4 --id merscope_2 --out cache
```

```bash
python web/SpatialViewer/prepare.py stereoseq --image data/stereoseq_1/C04042E3_HE_regist.tif --gef data/stereoseq_1/C04042E3.raw.gef --cellbin data/stereoseq_1/C04042E3.cellbin.gef --id stereoseq_1 --out cache
```

```bash
python web/SpatialViewer/prepare.py geomx --dcc data/geomx_1/mu_dev_E11_workflow_and_count_files/workflow/dcc/E9-E11-E15_DCC_unzipped --pkc data/geomx_1/mu_dev_E11_workflow_and_count_files/workflow/pkc/Mm_R_NGS_WTA_v1.0_unzipped/Mm_R_NGS_WTA_v1.0.pkc --annotation data/geomx_1/mu_dev_E11_workflow_and_count_files/workflow/readout_package/E11E15readoutpackage/SZ021_20211213T1749_LabWorksheet.txt --annotation data/geomx_1/mu_dev_E11_workflow_and_count_files/count/Export1_InitialDataset.xlsx --image data/geomx_1/mu_dev_E11_image_files/mu_dev_E11_004.ome.tiff --image data/geomx_1/mu_dev_E11_image_files/mu_dev_E11_005.ome.tiff --id geomx_1 --out cache
```

**Export.** The same command for each of the seven ids (`xenium_pub`, `visiumhd_pub`,
`phenocycler_pub`, `cosmx_1`, `merscope_2`, `stereoseq_1`, `geomx_1`), after its prepare:

```bash
python web/SpatialViewer/prepare.py export --id <id> --out cache --dest export --fmt jpg
```

**Opening.** For each of the seven exports above, once unthrottled and once on the home
connection. Chromium inherits the affinity:

```bash
start "" /affinity FF0000 /wait /b node tools/perf/bench_open.js --data-root export --ds <id> --net none --reps 3 --genes top,median --csv open.csv --out open_logs
```

```bash
start "" /affinity FF0000 /wait /b node tools/perf/bench_open.js --data-root export --ds <id> --net home --reps 3 --genes top,median --csv open.csv --out open_logs
```

**Draw a region → Run UMAP**, on the 8 µm bins of `visiumhd_pub`. The regions are three circles
around one centre (the median position of the bins with at least 10 UMI) that hold 1,000,
10,000 and 50,000 such bins, each a 72-point polygon as the page's circle tool draws it, in a
JSON file `{"regions": [{"target": 1000, "points": [[x, y], ...]}, ...]}` (image pixels).
The script that chose them is not in this repository.
<!-- RELEASE: add the circle-choosing script (select_regions.mjs) to tools/perf/, or publish the
     regions file of visiumhd_pub, and give its command here. -->

```bash
start "" /affinity FF0000 /wait /b node tools/perf/bench_open.js --mode region --data-root export --ds visiumhd_pub --layer 008um --regions regions_visiumhd_pub.json --leiden off,on --reps 3 --timeout-s 1800 --csv region.csv --out region_logs
```

**The tables:**

```bash
python tools/perf/summarize.py perf.csv open.csv region.csv
```

The paper's runs read the data from, and wrote the cache and exports to, a network share; a
local disk may be faster.

## 11. The figures

**Figure 2** is the notebook `notebooks/SODA_pipeline.ipynb` run on `merscope_crop` (commit
1ba11a5, on the efficiency cores). Change only these settings, in a copy of the notebook:

| Setting | Value |
|---|---|
| `WORK` | a folder of your own (its `inputs` folder holds `merscope_crop`) |
| `OPTIONS` | `merscope_crop`: label "MERSCOPE region (Vizgen sample data)" |
| `TITLE` | "Public example data" |
| `IS_THIS_DATA_PUBLIC` | "yes" |
| `ALWAYS_CHECK_IDS` | `True` (so Section 5 prints the sample ID check) |
| `CLOUD_URL`, `SODA_PAGE` | two local servers (`static_serve.py`) for the export and for `web/SpatialViewer`, so nothing is uploaded |

The panels are screenshots, in headless Chromium, of the executed notebook, of the export
folder and of the page opening the export with the link the notebook printed.

**Figure 3** is the page on `xenium_pub`, `visiumhd_pub` and `phenocycler_pub`, exported as in
section 10, or the copies of section [5](#5-the-public-example-exports). To use your own
exports, serve them and the page as in
[HOSTING.md, On your own computer](docs/HOSTING.md#a-on-your-own-computer):

```bash
python web/SpatialViewer/static_serve.py export 8765
```

```bash
python -m http.server 8766 --directory web/SpatialViewer
```

and open `http://localhost:8766/index.html?static=http://localhost:8765&ds=xenium_pub` (or
`ds=visiumhd_pub`, `ds=phenocycler_pub`).

<!-- RELEASE: give the view link of each Figure 3 panel (Save view link in the page). -->

## If your result differs

- Check the versions first: `python tools/check_install.py`, and the `versions` block in
  `expr/<layer>/umap/index.json` of both exports.
- A different Python, operating system or package version may move UMAP coordinates
  slightly and so renumber clusters.
- Open an issue with the bug report template (see [CONTRIBUTING.md](CONTRIBUTING.md)).
