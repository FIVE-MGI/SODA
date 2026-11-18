# Reproducing the paper's numbers and figures

This page says how to rebuild what the SODA paper reports, with the version of SODA it used
(0.9.0, see [CHANGELOG.md](CHANGELOG.md)). Every command below is in this repository.

Three kinds of result, three kinds of promise:

| Kind | Example | What you should get |
|---|---|---|
| Synthetic data and the files made from it | the demo dataset, its UMAP and clusters, the smoke test | the same files and the same checks passing, with the same package versions |
| Times and memory | the tables in [docs/PERFORMANCE.md](docs/PERFORMANCE.md) | the same order of magnitude; they depend on the machine and on the data |
| Results on real datasets | timings of real Visium HD, Xenium, CosMx ... sections | the same commands, on the same public or shared dataset (see [Real datasets](#real-datasets)) |

## 1. Get the version the paper used

```bash
git clone <this repository's URL> soda
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
| Synthetic demo data | `--seed` of each generator: `make_demo.py` 0, `make_xenium_demo.py` 7, `make_geomx_demo.py`, `make_imaging_demo.py`, `make_stereoseq_demo.py`, `make_stdformats_demo.py` 11 |

**Threads do not change the result.** `prepare.py` caps the maths libraries' threads (at
most 8; `--threads <n>` or `SODA_THREADS=<n>` to change it). The thread count changes only the
last bits of some floating-point values; the exported UMAP coordinates and clusters came out
byte-identical with 1, 8 and 24 threads on the datasets compared. The smoke test checks this
(see `web/SpatialViewer/threadcap.py`).

## 5. The synthetic demo data

Nothing in it comes from a tissue, a patient or an instrument. Each generator draws every
pixel and count from its seed, so the same `--seed` gives the same files.

```bash
python examples/make_demo.py --out demo/raw
```

writes a Visium HD-shaped folder (8 µm and 16 µm bins, an H&E-like image) and a 9-channel
fluorescence image. The README's [Try it in 5 minutes](README.md#try-it-in-5-minutes) prepares,
exports and opens it. The other generators each take `--out` and `--seed`:

| Generator | Makes |
|---|---|
| `examples/make_xenium_demo.py` | a Xenium output folder (cells, transcripts, stains) |
| `examples/make_geomx_demo.py` | a GeoMx DSP study (DCC, PKC, annotation sheets, scans with ROIs) |
| `examples/make_imaging_demo.py` | CosMx and MERSCOPE folders |
| `examples/make_stereoseq_demo.py` | Stereo-seq output (bin GEF, cell bin GEF, GEM, registered image) |
| `examples/make_stdformats_demo.py` | spatial AnnData (`.h5ad`, AnnData-Zarr), OME-Zarr and SpatialData stores |

The figure in the README (`docs/img/demo_side_by_side.jpg`) is this demo in Side by side
mode, coloured by Leiden cluster.

## 6. The smoke test

```bash
python tests/smoke_test.py
```

From the repository folder, with the environment active. It makes all the synthetic data
above, prepares and exports it with the README's commands, and checks every file: the tiles,
the gene data, that the clusters match the regions the demo was drawn from, the sample ID
check, and much more (the full list is at the top of `tests/smoke_test.py`). It ends with
**SMOKE TEST PASSED: N checks**.

- `--root <dir>` writes somewhere other than `demo/`; `--reuse` keeps earlier outputs;
  `--checks-only` runs no commands and only checks the outputs already there.
- The log is `demo/smoke_test.log`, the summary `demo/smoke_test_summary.json`.
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

Run the smoke test first: it makes `demo/export`, which the browser tests open. The smoke test
already runs most of them. To run one by itself:

| Command | What it checks |
|---|---|
| `node tests/browsers_test.js` | the page on every demo dataset in Chromium, WebKit (Safari's engine), Firefox and Edge, plus phone-sized windows; pixel comparison with Chromium |
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
cores, so treat the times as an order of magnitude.

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
  CPUs 16-23, the efficiency cores of that CPU; another CPU needs its own mask). The affinity
  the command had is recorded in the output.
- **Draw a region times:** `bench_open.js --mode region` with `--layer`, `--regions <circles.json>`
  and `--leiden off,on` (see the top of `tools/perf/bench_open.js`).
- Each script's header lists every option and what each number means.

## Real datasets

The performance tables use one real dataset per kind (listed in
[docs/PERFORMANCE.md](docs/PERFORMANCE.md#preparing-one-real-dataset-per-kind)). They are not
in this repository. To repeat those rows, prepare the same dataset with
`python web/SpatialViewer/prepare.py auto --src <folder> --id <id> --out cache` and export it
with `--fmt jpg` (as in PERFORMANCE.md), under `bench_prepare.py` as above.

<!-- RELEASE: name each real dataset and where to get it (public download or on request),
     and map each figure panel of the paper to its dataset and view link. -->

## If your result differs

- Check the versions first: `python tools/check_install.py`, and the `versions` block in
  `expr/<layer>/umap/index.json` of both exports.
- A different Python, operating system or package version may move UMAP coordinates
  slightly and so renumber clusters; the smoke test checks that clusters still follow the
  drawn regions, not exact coordinates.
- Open an issue with the bug report template (see [CONTRIBUTING.md](CONTRIBUTING.md)).
