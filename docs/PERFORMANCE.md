# Performance: how long preparing and opening take

These are measurements from one workstation, not guarantees. They show the order of magnitude to
expect and how to measure your own setup with the same scripts. The full tables, every
repetition and the exact conditions are in the benchmark report kept with the raw results
(outside the repository).

## The machine

- CPU: Intel Core i7-13700KF. **All runs used only its 8 efficiency cores** (logical CPUs
  16-23, `start /affinity FF0000` on Windows); the performance cores were not used.
- 32 GB RAM, Windows 11, Python 3.13 (NumPy 2.5), Chromium 141 (Playwright 1.56, headless, real GPU
  through ANGLE/Direct3D 11: NVIDIA RTX 3070 Ti).
- Raw data and the prepare cache on network shares (SMB). Exports were written to the network
  share and, for comparison, to a local SSD.
- The page was served from a local HTTP server that answers like blob storage (Range
  requests, ETag/Last-Modified, no Cache-Control header). Network throttling, where used, is
  Chrome DevTools emulation of a home connection: 50 Mbit/s down, 10 Mbit/s up, 20 ms latency.

Other jobs were sometimes running on the same cores: about 2 cores of unrelated scripts,
virus scanning while files were written, and another export. A quiet machine with 8
performance cores should be faster. Measured in November 2026.

## Preparing (one real dataset per kind)

| Dataset | Prepare | Export | Peak memory | Export size |
|---|---|---|---|---|
| README demo, end to end (make + prepare + export, both datasets) | 2.6 min in total | | 0.6 GB | 80 MB |
| Visium HD section (345,179 bins of 8 µm + 89,525 of 16 µm, 18,132 genes) | 1.8 min | 8.2 min | 2.0 GB | 0.98 GB, 1,792 files |
| Xenium mouse brain section (54,002 × 40,815 px, 155,754 cells, 209 M transcripts) | 16 min | 5.6 min | 1.7 GB | 1.75 GB, 12,032 files |
| Stereo-seq whole mouse brain chip (bins 20/50/100 + 130,588 cells) | 5.2 min | 18 min | 8.0 GB | 6.5 GB, 3,541 files |
| CosMx slide (20 FOVs, 91,972 cells, 26 M transcripts) | 3.8 min | 4.0 min | 0.9 GB | 0.48 GB, 1,750 files |
| Akoya PhenoCycler section (34 markers, 14,080 × 9,600 px) | 4.5 min | 3.0 min | 0.8 GB | 0.23 GB, 6,650 files |
| GeoMx DSP study (2 scans, 84 AOIs, mouse whole transcriptome) | 4.1 min | 0.7 min | 0.4 GB | 0.15 GB, 2,922 files |
| MERSCOPE region (65,029 cells, 2.35 M transcripts, no images) | 0.6 min | 1.9 min | 1.6 GB | 36 MB, 210 files |
| Classic Visium section (10x V1_Adult_Mouse_Brain: 2,702 spots of 55 µm, 32,285 genes, 11,291 × 11,757 px H&E) | 0.9 min | 1.0 min | 0.7 GB | 0.29 GB, 657 files |
| SpatialData Xenium store, zarr v3 (sandbox Xenium rep1: 167,780 cells, 35,416 × 25,778 px, 34 M transcripts) | 3.2 min | 6.4 min | 1.1 GB | 0.45 GB, 2,343 files |

- Exports used `--fmt jpg` with UMAP and Leiden clustering on (the defaults otherwise).
- Most of a large export is the whole-section UMAP and clustering. For the Visium HD 8 µm layer
  that is 5.7 min of the 8.2; for Stereo-seq bin20 it is 10.6 min of the 18.
- Most of the Xenium prepare is the image pyramid: 13.9 min of the 16.
- Classic Visium: `prepare.py visium` on Space Ranger 1.1.0's outputs as 10x hands them out
  (`V1_Adult_Mouse_Brain_filtered_feature_bc_matrix.h5`, `spatial/`, the 0.4 GB full-resolution
  TIFF), conditions as above: input and cache on the network share (`input=D: (network)
  cache=D: (network)`), export to it (`cache=D: (network) dest=D: (network)`), the 8 efficiency
  cores. Prepare peak 0.52 GB, export 0.73 GB.
- SpatialData Xenium (zarr v3, read by `zarrv3.py`; `--pixel-size 0.2125`): the same conditions.
  The image pyramid is 2.0 min of the 3.2; the export's UMAP and clustering 4.0 min of the 6.4.
  The export time is the first full export (UMAP on) before the molecules were added; adding them
  was a second, incremental export of 0.8 min. Peak 0.86 GB prepare, 1.13 GB export.
- Writing the export to the local SSD instead of the network share made no consistent difference.

## Opening a dataset

Times run from navigation start, as the median of 3 runs. *Cold* is a fresh browser profile.
`bench_open.js` now opens about:blank (untimed, then 1 s) before the cold visit, so cold means an empty cache, not a browser still starting; the cold numbers below predate that and include browser start-up.
*Warm* is a second visit in the same browser.

| | Unthrottled, cold | Home connection, cold | Warm |
|---|---|---|---|
| First image tile on screen | 2.3-3.9 s | 3.5-5.3 s | 1.4-1.8 s |
| Whole first view drawn | 2.4-4.0 s | 3.9-5.9 s | 1.5-2.1 s |
| Downloaded for the open (deck.gl included) | 2.1-13 MB | same | none |
| A gene drawn, timed from pressing Enter | 0.04-0.97 s | 0.08-2.0 s | 0.04-0.66 s |

- Some datasets in detail:
  - Visium HD: first tile 2.3 s, 9.3 MB, its most widely expressed gene (2.2 MB) in 0.45 s.
  - Xenium: 3.0 s and 6.4 MB.
  - Stereo-seq: 2.9 s and 13 MB, its most widely expressed gene (6.3 MB) in 0.97 s.
- Most of the cold time is the page starting up: the dataset is only asked for once the page
  is loaded (1.9-2.9 s cold, 1.3 s warm). Everything specific to the dataset takes 0.2-1.4 s after that.
- The browser held 0.4-0.7 GB after a first open (all its processes, GPU process included). The
  page's JavaScript heap was 7-12 MB.

## Draw a region → Run UMAP

Visium HD 8 µm layer, in the browser's web worker, timed from the run starting to the UMAP drawn.
These were measured on 17 November 2026 while other test runs shared the same 8 efficiency cores,
so they are slower than an idle machine; the page before the speed-up (b7eb968) and now were run
in turn, twice each (both runs shown):

| Bins in the region | Leiden | Before | Now, every bin | Now, default |
|---|---|---|---|---|
| 1,000 | off | 2.5 / 3.3 s | 2.6 / 3.0 s | the same |
| 1,000 | on | 3.9 / 4.4 s | 3.2 / 2.7 s | the same |
| 10,000 | off | 16.1 / 19.7 s | 16.7 / 14.9 s | the same |
| 10,000 | on | 36.2 / 32.7 s | 19.3 / 20.1 s | the same |
| 50,000 | off | 120 / 118 s | 90 / 104 s | 12.1 / 14.1 s (fast mode) |
| 50,000 | on | 228 / 196 s | 93 / 142 s | 16.7 / 33.0 s (fast mode) |

Every bin gives the same layout and clusters as before, bit for bit in the same browser. Above
20,000 bins the default is fast mode (a 10,000-bin sample; see PREPARING_DATA.md). On every bin,
the layout's epochs are most of the time (63 of the 93 s at 50,000 bins, Leiden on), then the
nearest neighbours (11 s) and Leiden (15 s, run alongside the layout). Changing only the
resolution of a 50,000-bin run re-clustered the kept layout in 14-24 s.

On a CPU about 4 times slower (the page capped at a quarter of one core), with Leiden off: 1,000
bins 8.3 s before and 7.3 s now, 10,000 bins 53 s and 35 s, 50,000 bins 401 s and 184 s (40 s in
fast mode). With Leiden on: 10,000 bins 112 s and 102 s, 50,000 bins 593 s and 497 s (118 s in
fast mode).

## What was not measured

- Visium HD's 2 µm layer.
- A whole Akoya slide (one section was cropped out of it).
- Real hosting over the internet. Only deck.gl came from its CDN, and the data came from a local server.
- Browsers other than Chromium, and phones.

## Measure your own

The scripts are in `tools/perf/`:

| Script | What it does |
|---|---|
| `bench_prepare.py` | runs any `prepare.py` command (or anything else) after `--`, and records wall time, peak memory of the whole process tree, CPU seconds, input and output size and file count; appends rows to a CSV |
| `bench_open.js` | serves the page and an export like blob storage, opens the dataset in headless Chromium (Playwright) cold and warm, optionally throttled, and records the time to the first tile, to the whole first view, to a gene after typing it, the bytes downloaded and the browser's memory. `--mode region` times Draw a region → Run UMAP for given circles |
| `summarize.py` | turns the CSV rows into Markdown tables (medians and ranges over repetitions) |

```bash
# preparing (Windows: prefix with  start "" /affinity FF0000 /wait /b  to pick cores)
python tools/perf/bench_prepare.py --label my_ds --step export --csv perf.csv --logdir perf_logs \
    --measure export/my_ds -- python web/SpatialViewer/prepare.py export --id my_ds --out cache --dest export

# opening (needs Playwright: NODE_PATH pointing at its node_modules)
node tools/perf/bench_open.js --data-root export --ds my_ds --reps 3 --genes top,median --csv perf.csv --out perf_logs
node tools/perf/bench_open.js --data-root export --ds my_ds --reps 3 --net home --csv perf.csv --out perf_logs

python tools/perf/summarize.py perf.csv
```

`bench_prepare.py` needs `psutil` in the Python that runs it (the command after `--` can use
another Python). Each script's header lists every option and what each number means.
