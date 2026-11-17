#!/usr/bin/env python
r"""End-to-end smoke test: synthetic demo data -> prepare -> export -> check every file.

    python tests/smoke_test.py                 (from the repo root; writes demo/)
    python tests/smoke_test.py --root <dir>    (anywhere else)
    python -m pytest tests/smoke_test.py       (same test, if pytest is installed)

It runs the demo commands from the README, in this order, with the Python that runs
this file:

    python examples/make_demo.py --out demo/raw
    python web/SpatialViewer/prepare.py visium --id demo_visiumhd --image demo/raw/visiumhd/image.tif --outs demo/raw/visiumhd/outs --out demo/cache
    python web/SpatialViewer/prepare.py tiff   --id demo_multiplex --src demo/raw/multiplex.ome.tif --out demo/cache
    python web/SpatialViewer/prepare.py export --id demo_visiumhd  --out demo/cache --dest demo/export --markers examples/demo_markers.json
    python web/SpatialViewer/prepare.py export --id demo_multiplex --out demo/cache --dest demo/export --fmt jpg --markers examples/demo_protein_markers.json
    python web/SpatialViewer/prepare.py export --id demo_multiplex --out demo/cache --dest demo/export_sample --fmt jpg --markers examples/demo_protein_markers.json --umap-sample 600
    python web/SpatialViewer/prepare.py tiff   --id demo_tiled --src demo/raw/multiplex_tiled.ome.tif --out demo/cache
    python web/SpatialViewer/prepare.py export --id demo_tiled --out demo/cache --dest demo/export --fmt jpg --markers examples/demo_protein_markers.json --no-umap --seam-correct
    python web/SpatialViewer/prepare.py export --id demo_visiumhd --out demo/cache --dest demo/export_user --markers examples/demo_markers.json --user-layer 008um --user-obs demo/user/clusters.csv --user-umap demo/user/umap.csv

and then checks the static tree the viewer reads: manifest.json, a tile for every
channel / level / position, expr/<key>/positions.bin, genes.json and the per-gene
files, the hvg/ chunks (layout in web/SpatialViewer/hvgchunks.py), the umap/ folder
(web/SpatialViewer/sectionumap.py) including that the clusters are spatially coherent
and match the regions the demo was drawn from, and that `deid.py check` passes.
For the multiplex image it checks the protein grid layer (web/SpatialViewer/protgrid.py,
expr/grid016um/) the same way: files, row counts, values, umap/ and celltypes.json, and
that its clusters are coherent and follow the drawn regions; and the slide QC (qc.json,
qc_flags.bin) on the demo's built-in problems: an empty "open" channel and a weaker
duplicate (both not written), a saturated channel (flagged), a fold (its squares flagged
"everything bright"). The second multiplex export fits the UMAP on a 600-square sample and
places the rest (sectionumap.UMAP_SAMPLE path); its clusters must agree with the first.
demo_tiled is the same image with per-tile steps in two markers: the export must find the
tile grid, correct those two markers (seams.py) and recover the factors that were put in.
A float32 OME-TIFF (~0..66 with fractions; NaN, inf and negatives; a constant channel) must be scaled into the
16-bit cache (0 .. its 99.99th percentile -> 0 .. 65535, order kept, far more than 67 grey levels) and a uint16
one must come through exactly as before (check_float_image).
It also checks the built-in marker lists and that the channel-name cleanup is the same
in celltypes.py and lasso_worker.js.
Stereo-seq (examples/make_stereoseq_demo.py, web/SpatialViewer/stereoseq.py): a synthetic chip as a
whole-chip bin GEF with a tissue contour + a cell bin GEF, as a tissue GEM, as an older multi-level
tissue GEF and as a shuffled whole-chip GEM with a mask image (on the grey ssDNA image). Checked: the tissue cut, bin totals,
positions, per-gene files and barcodes of the bin20 / bin50 layers against the drawn DNB counts, the
cells layer against the cell table, that all four routes give the same layers, and the image under
bins and cells.
The last export brings "the lab's own" results for the 8 µm layer (userclusters.py): a
Seurat-style clusters CSV (seurat_clusters, cell_type with a colour column, a numeric column)
whose barcodes end in "-1-0" instead of "-1" and lack a tenth of the bins, and a UMAP CSV.
Checked: the match rate and rule, index.json / manifest names and sources, the files, labels
and colours, markers and suggested names per user group; then a partial match (warning, import
goes ahead), no match and a merged two-sample file (both fail loudly, nothing changed).
The synthetic GeoMx DSP study (examples/make_geomx_demo.py: DCC files, two PKC files, an
annotation workbook / CSV / lab worksheet, two scans with ROI polygons and masks) is prepared
twice, with its scans and without (a drawn grid): one row per ROI segment, target values =
the generator's own probe collapse, negative probes out of the genes and behind the LOQ,
positions inside the ROI outlines, the annotation columns as groupings with no UMAP, the GeoMx
QC values, the Q3 normalisation factors (= the generator's plain-loop GeomxTools Q3), and no slide
name or Sample_ID in the export (check_geomx).
Standard formats (examples/make_stdformats_demo.py, web/SpatialViewer/stdformats.py, annstore.py, ngff.py): a spatial
AnnData .h5ad in the Squidpy layout (classic Visium: the hires image, scale factors, obs clusters with colours, X_umap),
an AnnData-Zarr of cells without an image, an OME-Zarr image (NGFF 0.4), and two SpatialData stores (Xenium-like:
circles, cell / nucleus polygons, transcripts, a second coordinate system; Visium HD-like: equal squares under an image
with a scale). Checked (check_stdformats): positions through scale factors and NGFF transformations, counts, genes, the
image pixels, outlines, molecules, the file's own clusterings, no made-up sample ID in any message or export file, the
cases that are refused (zarr v3, a plate, two libraries, a labels region, mapAxis), and the pages in Chromium and
Firefox (tests/stdformats_browser_test.js).
Xenium inputs (web/SpatialViewer/xenium.py, zipinput.py): the Xenium demo without its morphology_focus/ (cells and
molecules on a transcript-density background whose pixels are two 0.5 µm morphology pixels, origin at Xenium's 0, 0),
zipped as a whole output folder (deflated, two folders down) and as Xenium's own bundle (stored, files at the top),
and its cells-only variant zipped and read with every compressed member copied out (SODA_ZIP_MEMORY_MB=0). The zips
must give the folder's cache, leave no copies behind and not change the zip; zips with a '..' or an absolute member
name, a member that holds more than it declares, more than SODA_ZIP_MAX_GB declared, or too little free space must be
refused with the reason (check_xenium_inputs).
Cell outlines (web/SpatialViewer/cellshapes.py): the Xenium fixture carries cell and nucleus boundaries as
csv.gz, parquet, cells.zarr.zip (both layouts) and "a lab's own" GeoJSON + vertex table, CosMx CellLabels
images and a polygons file, MERSCOPE an older cell_boundaries HDF5 and a WKB parquet with z-planes, Stereo-seq
its cellBorder. Checked (check_outlines): the tile format byte for byte, every outline against the polygon
it came from (traced ones against their label mask, by IoU), rows (the files are not in matrix order), a
cell's centre inside its own outline, cells in parts, holes, crumbs, missing outlines, simplification, and
that outlines in the wrong units or of another run are refused and change nothing, and --no-outlines. The
browser half of that is tests/outlines_browser_test.js (Playwright, run here when node is installed): the outlines
in Chromium, and in Chromium, WebKit and Firefox what draws a subset of the cells with them (Draw a region, the region
click-filter, Colour by a value, cell sets, Save image and its legend, a link to the view), checked by pixel.
One control for each job (check_layout; in Chromium, WebKit and Firefox by tests/controls_browser_test.js):
Show on slide and Clear selection once (the chip above the picture), Image opacity at 0 for the old Image off
button, one run-or-show button for the volcano and the plots; old links and Save JSON files that had the image
off open with it off.
Last, the page itself: every local script and
stylesheet index.html names exists and carries a ?v= string, diag.js (the ?diag=1 panel)
loads before viewer.js, viewer.js still has its three tile decoders in order and reports
tile and graphics failures, and, when node is installed, tests/diag_test.js passes (the
panel's token-stripping). The browser half is tests/phone_browser_test.js (Playwright); tests/browsers_test.js
runs the page on this export in Chromium, WebKit (Safari's engine), Edge, Firefox, an iPhone-sized WebKit and Chromium,
and a phone-sized Firefox (the window size only: Playwright's Firefox cannot emulate touch) (Playwright). A browser
that is not installed here (Edge, Playwright's Firefox) is SKIPPED by the test with a line saying so.
tests/drawregion_test.js (Playwright, Chromium): Draw a region's worker on a seeded synthetic region gives the
layout and clusters stored in tests/fixtures/drawregion_exact.json (written by the worker before its speed work),
re-clusters a kept layout, and its fast mode (a sample laid out, the other bins placed by their nearest sample
bins) is seeded and agrees with the exact path; on this export, fast mode is labelled (panel, message, Save image
legend), recorded by Save JSON and Save CSV, undone by Run on every bin, and kept by a link to the view; merged
clusters saved with a fast-mode run come back only on a fast-mode run of the same sample.
Also: prepare.py caps the BLAS / OpenMP threads before numpy loads (threadcap.py: default,
SODA_THREADS, --threads, a variable already set), the PCA runs on one thread and gives the same
numbers as on the cap's threads; and the dataset picker (dspicker.js) is wired into the page with its
roles, and every .js of the page parses when node is installed. The tabbed layout (layout.js): the
tab rail, its four tabs and tab panels with their roles and ids, every section.panel named into one of
the four tabs (those in index.html and those the scripts build), no id twice, the controls layout.js
moves by id all there, and the hooks the other scripts need.
Your own results (web/SpatialViewer/ownres.js): a reader's CSV / TSV (also .gz) with a UMAP, groupings and
numbers, matched to the bins by barcode in the browser; wired into the page (check_ownres), and in Chromium,
WebKit and Firefox by tests/ownres_browser_test.js on export_user/demo_visiumhd: Seurat and scanpy files, partial
matches, roles put right, colours by pixel, the volcano on a loaded grouping, Save JSON, the view link, junk
files and a 375 px phone.
Automatic detection (web/SpatialViewer/detect.py): every synthetic folder above must be named with the right
kind and variant and verdict "certain"; `prepare.py auto` on five of them (and once through --kind / --image)
must write the same cache as the explicit subcommand, manifest and files; and folders built here from stub
files (empty, two platforms mixed, Space Ranger output with no / two / one right-sized image, a 10x matrix
folder with and without spatial/, a zip, missing stains, a missing path ...) must get the planned verdict and
wording, with `auto` stopping (exit 3 / 4 / 5) and writing nothing. The notebook must be valid, without
outputs, compile, and find its datasets through detect.py.
Batch export (web/SpatialViewer/batch.py, `prepare.py batch`) over a folder of dataset folders built from the
synthetic ones (batch/inputs): three that detection is sure of, one that only the options file settles ("kind"),
one that is certain but cannot be read, one with an open choice and one whose name looks like a sample ID. The
good ones must export (two at a time) and be in datasets.json, the broken one must fail with its last lines
logged without stopping the others, the other two must be skipped with what to add, the exit code must be 1, the
log (.md / .json / .csv) must hold time, peak memory and size per dataset and pass `deid.py check`, a second run
must skip what is done, and a dry run must write nothing.
Sharing (viewlink.js: the link that
restores a view; compare.js: linked windows): both load before viewer.js with ?v= strings, the ids they
build are there and clash with none of the page's, viewer.js hands them its functions and has one
serialiser for Save JSON and the link, nothing in them needs an inline handler or eval (the page's
Content-Security-Policy), the access token is kept out of storage and addresses, and, when node is
installed, tests/viewlink_test.js passes (encode / decode round trips, the length budget, oversize,
truncated, malformed and hostile links).
Offline use (libs.js, vendor.py, sw.js, offline.js): deck.gl and umap-js pinned by version and SHA-384 in one
place (libs.js) and the same in index.html's vendor/ tag, the worker and vendor.py; vendor/ git-ignored; the
Content-Security-Policy no wider (no eval of any kind); every worker named with a ?v= string (the service
worker finds them that way); sw.js's marker list = markers/; and tests/offline_test.js (Playwright,
Chromium, WebKit and Firefox): the page and a kept dataset open with no network, Run UMAP works offline, a new ?v=
replaces the cache, the libraries from vendor/ instead of the CDN, and vendor.py.

Without the optional clustering packages (requirements-clustering.txt) the export
skips umap/ and so do these checks, with a message saying so. The earlier outputs
of these two demo datasets are removed first so nothing stale is checked; pass
--reuse to keep them. Exit code 0 = all checks passed, 1 = a check failed,
2 = a command failed. Needs only the standard library and numpy to run the checks.
"""
from __future__ import annotations

import argparse
import gzip
import importlib.util
import json
import math
import os
import re
import shutil
import struct
import subprocess
import sys
import time
import zipfile
from pathlib import Path
from urllib.parse import quote

import numpy as np

REPO = Path(__file__).resolve().parents[1]
APP = REPO / "web" / "SpatialViewer"
VIS, MX, TL = "demo_visiumhd", "demo_multiplex", "demo_tiled"
# the lab's own names for the demo's drawn regions (user clusters test), and colours for some
USER_NAMES = {"fascicle": "Schwann cells", "perineurium": "Perineurial cells", "vessel": "Endothelial cells",
              "stroma": "Fibroblasts", "adipose": "Adipocytes", "edge": "Edge"}
USER_COLORS = {"Schwann cells": "#e41a1c", "Perineurial cells": "#377eb8", "Endothelial cells": "#4daf4a",
               "Fibroblasts": "#984ea3", "Adipocytes": "#ff7f00", "Edge": "#a65628"}
# synthetic Xenium folders (examples/make_xenium_demo.py): everything / transcripts.parquet without
# cell files / cells (parquet + h5) without transcripts
XE, XP, XN = "demo_xenium", "demo_xenium_pq", "demo_xenium_nt"
XEN_DS = (XE, XP, XN)
# Xenium inputs (xenium.py without stain images, zipinput.py): the demo folder without morphology_focus/, the demo
# zipped as a whole output folder (deflated, two folders down) and as Xenium's own bundle (stored, files at the top)
XNI, XZ, XZS = "demo_xenium_noimg", "demo_xenium_zip", "demo_xenium_bundle"
XIN_DS = (XNI, XZ, XZS)
# the lab's own names for the Xenium demo's three cell types (user clusters on cells), with colours
XE_USER_NAMES = ["Left-band cells", "Middle-band cells", "Right-band cells"]
XE_USER_COLORS = ["#1b9e77", "#d95f02", "#7570b3"]
# synthetic GeoMx DSP study (examples/make_geomx_demo.py): with its two scans / without images (drawn grid)
GX, GXG = "demo_geomx", "demo_geomx_grid"
GEOMX_DS = (GX, GXG)
# a small synthetic RNA-seq set written by rnaseq.export_rnaseq in-process (gene packs of that kind)
RNA = "demo_rnaseq"
# synthetic Stereo-seq folder (examples/make_stereoseq_demo.py): raw GEF + cell bin GEF / tissue GEM /
# older tissue GEF with stored bin10 + bin50 / shuffled whole-chip GEM cut with the mask image
ST, STG, STL, STM = "demo_stereoseq", "demo_stereoseq_gem", "demo_stereoseq_lv", "demo_stereoseq_mask"
ST_DS = (ST, STG, STL, STM)
# float32 / uint16 OME-TIFF fixtures (write_float_fixtures): the cast into the 16-bit cache
FL, U16, FL_SHAPE = "demo_float32", "demo_uint16", (192, 320)
NO_CLUSTER = 0xFFFF
CLUSTER_MODULES = ("umap", "igraph", "leidenalg")
# a barcode or a local path must never reach the export
BARCODE = re.compile(r"s_\d{3}um_\d{5}_\d{5}|[ACGT]{16}-1")
LOCAL_PATH = re.compile(r"(?<![A-Za-z])[A-Za-z]:[\\/]|\\\\[A-Za-z0-9]|/home/|/Users/|/mnt/|/tmp/")


class Failures(Exception):
    pass


class Checker:
    def __init__(self):
        self.failed: list[str] = []
        self.n = 0

    def ok(self, cond, what: str) -> bool:
        self.n += 1
        if not cond:
            self.failed.append(what)
            print(f"  FAIL  {what}", flush=True)
        return bool(cond)

    def eq(self, got, want, what: str) -> bool:
        return self.ok(got == want, f"{what}: got {got!r}, want {want!r}")


# ------------------------------------------------------------------ running commands
def demo_commands(root: Path) -> list[list[str]]:
    """The README's demo commands. Paths stay relative when root is <repo>/demo."""
    try:
        r = root.resolve().relative_to(REPO.resolve()).as_posix()
    except ValueError:
        r = root.resolve().as_posix()
    raw, cache, export = f"{r}/raw", f"{r}/cache", f"{r}/export"
    prep = "web/SpatialViewer/prepare.py"
    return [
        ["examples/make_demo.py", "--out", raw],
        [prep, "visium", "--id", VIS, "--image", f"{raw}/visiumhd/image.tif", "--outs", f"{raw}/visiumhd/outs",
         "--out", cache],
        [prep, "tiff", "--id", MX, "--src", f"{raw}/multiplex.ome.tif", "--out", cache],
        [prep, "export", "--id", VIS, "--out", cache, "--dest", export, "--markers", "examples/demo_markers.json"],
        [prep, "export", "--id", MX, "--out", cache, "--dest", export, "--fmt", "jpg",
         "--markers", "examples/demo_protein_markers.json"],
        [prep, "export", "--id", MX, "--out", cache, "--dest", f"{r}/export_sample", "--fmt", "jpg",
         "--markers", "examples/demo_protein_markers.json", "--umap-sample", str(SAMPLE)],
        [prep, "tiff", "--id", TL, "--src", f"{raw}/multiplex_tiled.ome.tif", "--out", cache],
        [prep, "export", "--id", TL, "--out", cache, "--dest", export, "--fmt", "jpg",
         "--markers", "examples/demo_protein_markers.json", "--no-umap", "--seam-correct"],
        [prep, "export", "--id", VIS, "--out", cache, "--dest", f"{r}/export_user", "--markers", "examples/demo_markers.json",
         "--user-layer", "008um", "--user-obs", f"{r}/user/clusters.csv", "--user-umap", f"{r}/user/umap.csv", "--no-ids"],
        # the old layout, one file per gene, for hosts without Range requests (check_packs compares it with the packs)
        [prep, "export", "--id", VIS, "--out", cache, "--dest", f"{r}/export_loose", "--markers", "examples/demo_markers.json",
         "--no-umap", "--no-hvg", "--gene-files", "loose", "--no-ids"],
        # the barcode lists (layerids.py): written by every export by default; added here to export_user's Visium HD
        # (exported with --no-ids above); asked for explicitly with --ids (the old flag) on export_user's Xenium below
        [prep, "ids", "--id", VIS, "--out", cache, "--dest", f"{r}/export_user"],
        ["examples/make_xenium_demo.py", "--out", raw],
        [prep, "xenium", "--id", XE, "--src", f"{raw}/xenium", "--out", cache, "--label", "Synthetic Xenium"],
        [prep, "export", "--id", XE, "--out", cache, "--dest", export, "--markers", "examples/demo_xenium_markers.json"],
        [prep, "export", "--id", XE, "--out", cache, "--dest", f"{r}/export_user", "--markers", "examples/demo_xenium_markers.json",
         "--user-layer", "cells", "--user-obs", f"{r}/user/xenium_cells.csv", "--user-obs", f"{r}/user/xenium_groups.csv", "--ids"],
        # packs AND files, in tiny shards so the cells layer and the molecules each need several
        [prep, "export", "--id", XE, "--out", cache, "--dest", f"{r}/export_both", "--no-umap", "--grid-um", "0",
         "--gene-files", "both", "--shard-mb", "0.05"],
        [prep, "xenium", "--id", XP, "--src", f"{raw}/xenium_pq", "--out", cache],
        [prep, "export", "--id", XP, "--out", cache, "--dest", export, "--grid-um", "0"],
        [prep, "xenium", "--id", XN, "--src", f"{raw}/xenium_nt", "--out", cache],
        [prep, "export", "--id", XN, "--out", cache, "--dest", export, "--no-umap"],
        # no stain images: a transcript-density background; zips read as they are (write_xenium_inputs made them)
        [prep, "xenium", "--id", XNI, "--src", f"{raw}/xenium_noimg", "--out", cache],
        [prep, "export", "--id", XNI, "--out", cache, "--dest", export, "--no-umap"],
        [prep, "xenium", "--id", XZ, "--src", f"{raw}/xenium_output.zip", "--out", cache, "--label", "Synthetic Xenium"],
        [prep, "export", "--id", XZ, "--out", cache, "--dest", export, "--no-umap", "--no-hvg"],
        [prep, "xenium", "--id", XZS, "--src", f"{raw}/xenium_bundle.zip", "--out", cache, "--label", "Synthetic Xenium"],
        ["examples/make_geomx_demo.py", "--out", raw],
        [prep, "geomx", "--id", GX, "--dcc", f"{raw}/geomx/dcc", "--pkc", f"{raw}/geomx/pkc",
         "--annotation", f"{raw}/geomx/annotations.xlsx", "--annotation", f"{raw}/geomx/segment_properties.xlsx",
         "--image", f"{raw}/geomx/demoA.ome.tiff", "--image", f"{raw}/geomx/demoB.ome.tiff", "--out", cache,
         "--label", "Synthetic GeoMx"],
        [prep, "export", "--id", GX, "--out", cache, "--dest", export],
        [prep, "geomx", "--id", GXG, "--dcc", f"{raw}/geomx/dcc", "--pkc", f"{raw}/geomx/pkc",
         "--annotation", f"{raw}/geomx/Demo_LabWorksheet.txt", "--annotation", f"{raw}/geomx/annotations.csv", "--out", cache],
        [prep, "export", "--id", GXG, "--out", cache, "--dest", export],
    ]


def stereoseq_commands(root: Path) -> list[list[str]]:
    """The Stereo-seq demo commands (kept apart from demo_commands so the lists stay easy to merge)."""
    try:
        r = root.resolve().relative_to(REPO.resolve()).as_posix()
    except ValueError:
        r = root.resolve().as_posix()
    raw, cache, export = f"{r}/raw/stereoseq", f"{r}/cache", f"{r}/export"
    prep = "web/SpatialViewer/prepare.py"
    img = ["--image", f"{raw}/SYNCHIP_HE_regist.tif"]
    return [
        ["examples/make_stereoseq_demo.py", "--out", f"{r}/raw"],
        [prep, "stereoseq", "--id", ST, *img, "--gef", f"{raw}/SYNCHIP.raw.gef", "--cellbin", f"{raw}/SYNCHIP.cellbin.gef",
         "--bins", "20,50", "--out", cache, "--label", "Synthetic Stereo-seq"],
        [prep, "export", "--id", ST, "--out", cache, "--dest", export, "--markers", "examples/demo_xenium_markers.json"],
        [prep, "stereoseq", "--id", STG, *img, "--gem", f"{raw}/SYNCHIP.tissue.gem.gz", "--bins", "20,50", "--out", cache],
        [prep, "export", "--id", STG, "--out", cache, "--dest", export, "--no-umap"],
        [prep, "stereoseq", "--id", STL, *img, "--gef", f"{raw}/SYNCHIP.tissue.gef", "--cellbin", f"{raw}/SYNCHIP.cut.cellbin.gef",
         "--bins", "20,50", "--barcodes", "stereopy", "--out", cache],
        [prep, "export", "--id", STL, "--out", cache, "--dest", export, "--no-umap"],
        [prep, "stereoseq", "--id", STM, "--image", f"{raw}/SYNCHIP_ssDNA_regist.tif", "--gem", f"{raw}/SYNCHIP.raw.shuffled.gem",
         "--tissue-mask", f"{raw}/SYNCHIP_HE_tissue_cut.tif", "--bins", "50", "--out", cache],
        [prep, "export", "--id", STM, "--out", cache, "--dest", export, "--no-umap"],
    ]


def write_float_fixtures(root: Path):
    """Two small OME-TIFFs for the image cast (pyramid._float_scale): float32 as a PhenoChart / inForm component
    export holds it (~0..66 with fractions; a channel with negatives, NaN and +-inf; a constant channel), and uint16
    (must come through exactly as before)."""
    import tifffile
    d = root / "floatimg"            # not under raw/: detect.py reads every dataset there
    d.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(11)
    H, W = FL_SHAPE
    ramp = rng.uniform(0, 66, (H, W)).astype(np.float32)
    ramp[:4] = 0                                                   # a background strip
    odd = rng.uniform(-3, 66, (H, W)).astype(np.float32)
    odd[0, :8], odd[1, :8], odd[2, :8] = np.nan, np.inf, -np.inf
    flat = np.full((H, W), 5.0, np.float32)
    kw = dict(ome=True, photometric="minisblack")
    tifffile.imwrite(str(d / "float32.ome.tif"), np.stack([ramp, odd, flat]),
                     metadata={"axes": "CYX", "Channel": {"Name": ["Ramp", "Odd", "Flat"]}}, **kw)
    u = rng.integers(0, 4000, (2, H, W)).astype(np.uint16)
    u[:, :4] = 0
    tifffile.imwrite(str(d / "uint16.ome.tif"), u, metadata={"axes": "CYX", "Channel": {"Name": ["A", "B"]}}, **kw)


def float_commands(root: Path) -> list[list[str]]:
    """prepare + export of the float32 / uint16 fixtures (write_float_fixtures), into export_float/."""
    try:
        r = root.resolve().relative_to(REPO.resolve()).as_posix()
    except ValueError:
        r = root.resolve().as_posix()
    prep = "web/SpatialViewer/prepare.py"
    out = []
    for ds, f in ((FL, "float32"), (U16, "uint16")):
        out += [[prep, "tiff", "--id", ds, "--src", f"{r}/floatimg/{f}.ome.tif", "--out", f"{r}/cache", "--pixel-size", "0.5"],
                [prep, "export", "--id", ds, "--out", f"{r}/cache", "--dest", f"{r}/export_float", "--grid-um", "0", "--no-umap"]]
    return out


def clean_stereoseq(root: Path):
    for p in [root / "raw" / "stereoseq"] + [root / sub / ds for sub in ("cache", "export") for ds in ST_DS]:
        if p.exists():
            print(f"  removing earlier output {p}", flush=True)
            shutil.rmtree(p)


def write_user_files(root: Path):
    """The lab's own results for the demo's 8 um layer, as Seurat writes them: write.csv of the
    metadata (empty first header cell), barcodes "-1-0" instead of "-1", a tenth of the bins
    dropped in QC; a cell_type column with a colour column; a numeric column; a UMAP CSV."""
    rows = [ln.split(",") for ln in (root / "raw" / "visiumhd" / "truth_square_008um.csv")
            .read_text(encoding="utf-8").split("\n")[1:] if ln]
    regions = sorted({r for _, r in rows})
    out = root / "user"
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(3)
    cl, um = ['"","seurat_clusters","cell_type","cell_type_color","nCount_Spatial"'], ["barcode,UMAP_1,UMAP_2"]
    for i, (bc, reg) in enumerate(rows):
        if i % 10 == 0:
            continue
        k = regions.index(reg)
        name = USER_NAMES[reg]
        cl.append(f'"{bc}-0",{k},"{name}","{USER_COLORS[name]}",{int(rng.integers(10, 900))}')
        um.append(f"{bc},{k * 4 + rng.normal():.5f},{(k % 3) * 4 + rng.normal():.5f}")
    (out / "clusters.csv").write_text("\n".join(cl) + "\n", encoding="utf-8")
    (out / "umap.csv").write_text("\n".join(um) + "\n", encoding="utf-8")
    # partial: 40 % of the barcodes unknown; none: all unknown; merged: two samples on the same bins
    part = [cl[0]] + [ln.replace("s_008um", "zz_008um") if j % 5 < 2 else ln for j, ln in enumerate(cl[1:])]
    (out / "partial.csv").write_text("\n".join(part) + "\n", encoding="utf-8")
    none = [cl[0]] + [ln.replace("s_008um", "s_016um") for ln in cl[1:]]
    (out / "nomatch.csv").write_text("\n".join(none) + "\n", encoding="utf-8")
    merged = ["barcode,cluster"] + [f"{t}_{bc},{k % 5}" for t in ("A", "B") for k, (bc, _) in enumerate(rows[:2000])]
    (out / "merged.csv").write_text("\n".join(merged) + "\n", encoding="utf-8")


def write_xenium_user_files(root: Path):
    """The lab's own clusters for the Xenium demo's cells, as scanpy writes them from a spatialdata-io
    table (adata.obs.to_csv): the row index ("0", "1", ...) first under an empty header, then the
    Xenium cell id in a cell_id column, a leiden column, a cell_type column with a colour column.
    Every ninth cell is missing (dropped in QC). And Xenium Explorer's cell groups CSV for all cells:
    cell_id,group,color (a bare colour column)."""
    raw = root / "raw"
    ids = gzip.open(raw / "xenium" / "cell_feature_matrix" / "barcodes.tsv.gz", "rt").read().split()
    ctype = np.load(raw / "truth_xenium.npz")["ctype"]
    out = root / "user"
    out.mkdir(parents=True, exist_ok=True)
    lines = [",cell_id,leiden,cell_type,cell_type_color"]
    for i, (cid, t) in enumerate(zip(ids, ctype)):
        if i % 9 == 4:
            continue
        lines.append(f"{i},{cid},{(int(t) + 1) % 3},{XE_USER_NAMES[t]},{XE_USER_COLORS[t]}")
    (out / "xenium_cells.csv").write_text("\n".join(lines) + "\n", encoding="utf-8")
    grp = ["cell_id,group,color"] + [f"{cid},Band {'ABC'[t]},{XE_USER_COLORS[t]}" for cid, t in zip(ids, ctype)]
    (out / "xenium_groups.csv").write_text("\n".join(grp) + "\n", encoding="utf-8")


def run(args: list[str], log_path: Path) -> float:
    print(f"$ python {' '.join(args)}", flush=True)
    t0 = time.time()
    with open(log_path, "a", encoding="utf-8") as lg:
        lg.write(f"\n$ python {' '.join(args)}\n")
        lg.flush()
        # the child must print UTF-8 (its logs have µ): on Windows a piped child defaults to the
        # ANSI code page, which this side would then read as broken characters
        p = subprocess.run([sys.executable, *args], cwd=str(REPO), stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace",
                           env={**os.environ, "PYTHONIOENCODING": "utf-8"})
        lg.write(p.stdout)
    dt = time.time() - t0
    tail = "\n".join(p.stdout.rstrip().splitlines()[-6:])
    print("    " + tail.replace("\n", "\n    "), flush=True)
    if p.returncode != 0:
        print(p.stdout[-6000:], flush=True)
        raise SystemExit(f"command failed (exit {p.returncode}) after {dt:.0f}s: python {' '.join(args)}\n"
                         f"full log: {log_path}")
    print(f"    ({dt:.1f}s)", flush=True)
    return dt


def print_skipped(out: str):
    """The browser tests' SKIPPED lines, shown here: a browser that is not installed on this machine (Playwright's
    Firefox or Edge) is skipped by the test itself, which still passes, and says so with a line starting SKIPPED."""
    for ln in out.splitlines():
        if ln.strip().startswith("SKIPPED"):
            print("    " + ln.strip())


def clean(root: Path):
    """Remove only what this test itself writes, so nothing stale gets checked."""
    for p in [root / "raw", root / "user", root / "export_sample" / MX, root / "export_user" / VIS,
              root / "export_user" / XE, root / "export_partial" / VIS, root / "export_loose" / VIS,
              root / "export_both" / XE, root / "export" / RNA, root / "export_rnaseq_loose", root / "export_repack",
              root / "export_tamper", root / "xenin", root / "cache" / "_unzipped", root / "export_float", root / "floatimg",
              root / "cache" / FL, root / "cache" / U16] + [root / sub / ds for sub in ("cache", "export")
                                                            for ds in (VIS, MX, TL) + XEN_DS + XIN_DS + GEOMX_DS]:
        if p.exists():
            print(f"  removing earlier output {p}", flush=True)
            shutil.rmtree(p)


# ------------------------------------------------------------------ small helpers
def read_json(p: Path):
    return json.loads(p.read_text(encoding="utf-8"))


def enc(name: str) -> str:
    return quote(name, safe="-_.!~*'()")              # = JavaScript encodeURIComponent


def f32(p: Path, cols: int) -> np.ndarray:
    return np.fromfile(str(p), dtype="<f4").reshape(-1, cols)


# ---- gene files: packed (genepack.py, the default) or loose (gene/<NAME>.bin) -------------------
def app_module(name: str):
    sys.path.insert(0, str(APP))
    try:
        return importlib.import_module(name)
    finally:
        sys.path.pop(0)


_GROUPS: dict = {}


def gene_group(d: Path) -> tuple[dict | None, dict]:
    """An exported gene-group folder (expr/<key>, molecules, an RNA-seq set): its genePack (None =
    loose files) and {name: listing entry with s, o, b}."""
    gp = app_module("genepack")
    d = Path(d)
    for ds in (d.parent.parent, d.parent, d):
        mf = ds / "manifest.json"
        if not mf.is_file():
            continue
        stamp = (str(d.resolve()), mf.stat().st_mtime_ns)
        if stamp in _GROUPS:
            return _GROUPS[stamp]
        for g in gp.groups(read_json(mf)):
            if (ds / g["dir"] if g["dir"] else ds).resolve() == d.resolve():
                idx = {}
                if g["pack"]:
                    for e in gp.read_listing(ds, g):
                        if "s" in e:
                            idx.setdefault(e["n"], e)
                _GROUPS[stamp] = (g["pack"], idx)
                return _GROUPS[stamp]
    return None, {}


def gene_bytes(d: Path, name: str) -> bytes | None:
    """One gene's payload in an exported group: its slice of the pack, or its loose file. None = none."""
    pack, idx = gene_group(d)
    if pack:
        e = idx.get(name)
        if e is None:
            return None
        with open(Path(d) / pack["dir"] / pack["shards"][e["s"]]["file"], "rb") as f:
            f.seek(e["o"])
            return f.read(e["b"])
    p = Path(d) / "gene" / f"{enc(name)}.bin"
    return p.read_bytes() if p.exists() else None


def gene_f32(d: Path, name: str, cols: int = 3) -> np.ndarray | None:
    b = gene_bytes(d, name)
    return None if b is None else np.frombuffer(b, dtype="<f4").reshape(-1, cols).copy()


def gene_count(d: Path) -> int:
    """Genes with a payload in the group (pack entries, or loose files)."""
    pack, idx = gene_group(d)
    return len(idx) if pack else len(list((Path(d) / "gene").glob("*.bin")))


def ari(a: np.ndarray, b: np.ndarray) -> float:
    """Adjusted Rand index of two labelings (numpy only)."""
    _, ai = np.unique(a, return_inverse=True)
    _, bi = np.unique(b, return_inverse=True)
    ct = np.zeros((ai.max() + 1, bi.max() + 1), np.int64)
    np.add.at(ct, (ai, bi), 1)
    c2 = lambda x: x * (x - 1) / 2.0
    s = c2(ct).sum()
    sa, sb = c2(ct.sum(1)).sum(), c2(ct.sum(0)).sum()
    exp = sa * sb / c2(len(a))
    mx = (sa + sb) / 2
    return float((s - exp) / (mx - exp)) if mx != exp else 1.0


def grid_pairs(xy: np.ndarray, step: float) -> np.ndarray:
    """Pairs of rows that are 4-neighbours on the bin grid."""
    ij = np.rint(xy / step - 0.5).astype(np.int64)
    key = {(int(i), int(j)): r for r, (i, j) in enumerate(ij)}
    pairs = []
    for r, (i, j) in enumerate(ij):
        for di, dj in ((1, 0), (0, 1)):
            q = key.get((int(i) + di, int(j) + dj))
            if q is not None:
                pairs.append((r, q))
    return np.array(pairs, np.int64).reshape(-1, 2)


# ------------------------------------------------------------------ checks
def check_tiles(c: Checker, out: Path, m: dict, ext: str):
    img = m["image"]
    T, nl = img["tileSize"], img["levels"]
    magic = b"\x89PNG" if ext == "png" else b"\xff\xd8\xff"
    c.eq(len(img["levelShapes"]), nl, f"{m['id']}: levelShapes length")
    c.ok(max(img["levelShapes"][-1]) <= T, f"{m['id']}: coarsest level fits one tile")
    want = set()
    for L, (H, W) in enumerate(img["levelShapes"]):
        z = nl - 1 - L
        for ch in range(len(img["channels"])):
            for tx in range(math.ceil(W / T)):
                for ty in range(math.ceil(H / T)):
                    want.add(out / "tiles" / str(ch) / str(z) / str(tx) / f"{ty}.{ext}")
    have = set((out / "tiles").rglob(f"*.{ext}"))
    missing = want - have
    c.ok(not missing, f"{m['id']}: {len(missing)} tile(s) missing, e.g. {sorted(missing)[:3]}")
    c.ok(not (have - want), f"{m['id']}: unexpected tiles {sorted(have - want)[:3]}")
    c.ok(not list((out / "tiles").rglob("*.part")), f"{m['id']}: leftover .part files")
    bad = [p for p in sorted(want & have) if p.read_bytes()[:len(magic)] != magic]
    c.ok(not bad, f"{m['id']}: tiles that are not {ext}: {bad[:3]}")
    print(f"  tiles: {len(have)} {ext} files, {nl} levels x {len(img['channels'])} channels")


def check_manifest_common(c: Checker, m: dict, ds: str, ext: str):
    c.eq(m.get("id"), ds, "manifest id")
    c.eq(m.get("static"), True, f"{ds}: manifest static")
    c.eq(m.get("tileFormat"), ext, f"{ds}: manifest tileFormat")
    c.ok("source" not in m, f"{ds}: manifest still has 'source' (local paths)")
    img = m.get("image", {})
    for k in ("width", "height", "levels", "tileSize", "dtype", "levelShapes", "channels"):
        c.ok(k in img, f"{ds}: manifest image.{k} missing")
    for ch in img.get("channels", []):
        for k in ("index", "name", "color", "bits", "window"):
            c.ok(k in ch, f"{ds}: channel {ch.get('index')} has no {k}")


def check_visium(c: Checker, root: Path, clustering: bool) -> dict:
    out = root / "export" / VIS
    m = read_json(out / "manifest.json")
    check_manifest_common(c, m, VIS, "png")
    c.eq(m.get("kind"), "visium", "visium manifest kind")
    c.eq((m["image"]["width"], m["image"]["height"]), (4000, 4000), "visium image size")
    c.eq([ch["name"] for ch in m["image"]["channels"]], ["Red", "Green", "Blue"], "visium channel names")
    c.ok(abs((m.get("pixelSizeUm") or 0) - 0.32) < 1e-6, f"visium pixelSizeUm {m.get('pixelSizeUm')}")
    check_tiles(c, out, m, "png")
    layers = {L["key"]: L for L in m.get("expr", [])}
    c.eq(sorted(layers), ["008um", "016um"], "visium expr layers")
    summary = {}
    for key, L in sorted(layers.items()):
        print(f"  expr {key}:")
        summary[key] = check_layer(c, root, out, L, clustering)
    return summary


def check_layer(c: Checker, root: Path, out: Path, L: dict, clustering: bool) -> dict:
    key = L["key"]
    d = out / "expr" / key
    bin_um = float(key[:3])
    for k in ("label", "nBins", "nGenes", "binSizeUm", "spotDiameterPx", "micronsPerPixel"):
        c.ok(k in L, f"{key}: manifest layer has no {k}")
    c.eq(L.get("binSizeUm"), bin_um, f"{key}: binSizeUm")
    c.ok(abs(L.get("spotDiameterPx", 0) - bin_um / 0.32) < 1e-3, f"{key}: spotDiameterPx {L.get('spotDiameterPx')}")
    c.ok(abs(L.get("micronsPerPixel", 0) - 0.32) < 1e-6, f"{key}: micronsPerPixel {L.get('micronsPerPixel')}")

    # positions.bin: one float32 (x, y) per bin, all inside the image
    pos = f32(d / "positions.bin", 2)
    n = len(pos)
    c.eq((d / "positions.bin").stat().st_size, 8 * L["nBins"], f"{key}: positions.bin bytes (8 x nBins)")
    c.ok(np.isfinite(pos).all(), f"{key}: positions.bin has NaN")
    c.ok(((pos >= 0) & (pos <= 4000)).all(), f"{key}: positions outside the image")
    truth = [ln.split(",")[1] for ln in
             (root / "raw" / "visiumhd" / f"truth_square_{key}.csv").read_text(encoding="utf-8").split("\n")[1:] if ln]
    c.eq(len(truth), n, f"{key}: rows in positions.bin vs bins the demo wrote")

    # genes.json + each gene's [x, y, v] float32 (pack or gene/<NAME>.bin), values summing to the gene total
    genes = read_json(d / "genes.json")
    c.eq(len(genes), L["nGenes"], f"{key}: genes.json length")
    c.eq(len(genes), 300, f"{key}: demo gene count")
    row_of = {p.tobytes(): i for i, p in enumerate(pos)}
    lib = np.zeros(n, np.float64)
    dense = {}
    n_files, bad = 0, []
    for g in genes:
        if g["t"] <= 0:
            continue
        a = gene_f32(d, g["n"])
        if a is None:
            bad.append(f"missing {g['n']}")
            continue
        n_files += 1
        if abs(float(a[:, 2].sum()) - g["t"]) > 1e-3 * max(1.0, g["t"]):
            bad.append(f"{g['n']}: sum {a[:, 2].sum()} != total {g['t']}")
        rows = np.array([row_of.get(r.tobytes(), -1) for r in a[:, :2]])
        if (rows < 0).any():
            bad.append(f"{g['n']}: points not on positions.bin")
            continue
        np.add.at(lib, rows, a[:, 2])
        dense[g["n"]] = (rows, a[:, 2])
    c.ok(not bad, f"{key}: gene files: {bad[:5]}")
    c.eq(gene_count(d), n_files, f"{key}: genes with a payload (pack entries / gene files)")
    print(f"    positions.bin {n:,} rows, {n_files} genes, median {np.median(lib):.0f} UMI per bin")

    # a marker gene sits in the region it was drawn for
    fasc = np.array([t == "fascicle" for t in truth])
    if "MARKER_FASC1" in dense:
        rows, v = dense["MARKER_FASC1"]
        share = float(v[fasc[rows]].sum() / max(v.sum(), 1))
        c.ok(share > 0.8, f"{key}: only {share:.0%} of MARKER_FASC1 counts lie in fascicle bins")

    hv = check_hvg(c, d, L, pos, lib, dense)
    check_binqc(c, d, L, lib, dense, [g["n"] for g in genes])
    res = {"bins": n, "medianUmi": float(np.median(lib))}
    if clustering:
        res.update(check_umap(c, d, L, pos, lib, truth, bin_um / 0.32))
    else:
        c.ok("umap" not in L, f"{key}: umap in manifest although clustering packages are missing")
        print("    SKIP umap checks: " + ", ".join(CLUSTER_MODULES) + " not all installed "
              "(pip install -r web/SpatialViewer/requirements-clustering.txt)")
    res["hvgGenes"] = hv
    return res


def check_binqc(c: Checker, d: Path, L: dict, lib, dense, gene_names) -> None:
    """binqc/ (web/SpatialViewer/binqc.py): per row total counts, genes detected, % MT,
    recomputed here from the gene files."""
    key, q = L["key"], L.get("binQc")
    if not c.ok(q and (d / "binqc").is_dir(), f"{key}: manifest binQc / binqc/ missing"):
        return
    n = len(lib)
    c.eq(q.get("nRows"), n, f"{key}: binQc.nRows")
    mets = {m["key"]: m for m in q.get("metrics", [])}
    for k in ("total", "genes"):
        c.ok(k in mets, f"{key}: binQc has no {k} metric")
    tot = np.fromfile(d / "binqc" / mets["total"]["file"], dtype="<f4")
    c.eq(len(tot), n, f"{key}: binqc total.bin rows")
    c.ok(np.allclose(tot, lib, rtol=1e-5, atol=1e-3), f"{key}: binqc total != sum of the gene files per bin")
    ng = np.zeros(n, np.int64)
    for rows, v in dense.values():
        np.add.at(ng, rows[v > 0], 1)
    gdt = "<u2" if mets["genes"]["dtype"] == "uint16" else "<u4"
    gq = np.fromfile(d / "binqc" / mets["genes"]["file"], dtype=gdt)
    c.ok(len(gq) == n and (gq == ng).all(), f"{key}: binqc genes detected != genes with a count per bin")
    for k, a in (("total", tot), ("genes", gq)):
        c.ok(abs(mets[k]["min"] - float(a.min())) < 1e-3 and abs(mets[k]["max"] - float(a.max())) < 1e-3,
             f"{key}: binQc {k} min/max vs the file")
    mt = [g for g in gene_names if g.upper().startswith("MT-")]
    if mt:
        pm = np.fromfile(d / "binqc" / mets["pctMt"]["file"], dtype="<f4")
        mts = np.zeros(n)
        for g in mt:
            if g in dense:
                np.add.at(mts, dense[g][0], dense[g][1])
        want = np.where(lib > 0, 100 * mts / np.maximum(lib, 1e-12), 0)
        c.ok(np.allclose(pm, want, atol=1e-3), f"{key}: binqc pct_mt != MT- share per bin")
    else:
        c.ok("pctMt" not in mets and not (d / "binqc" / "pct_mt.bin").exists(),
             f"{key}: pct_mt written although the layer has no MT- genes")
    print(f"    binqc: {len(mets)} metrics, median {np.median(tot):.0f} UMI, {np.median(gq):.0f} genes"
          + ("" if mt else ", no MT- genes (no % MT)"))


def check_hvg(c: Checker, d: Path, L: dict, pos, lib, dense) -> int:
    key = L["key"]
    hd = d / "hvg"
    if not c.ok((hd / "index.json").exists(), f"{key}: hvg/index.json missing"):
        return 0
    idx = read_json(hd / "index.json")
    hg = read_json(hd / "genes.json")
    for k in ("chunkPx", "chunkUm", "nGenes", "minUmi", "chunks"):
        c.ok(k in idx, f"{key}: hvg/index.json has no {k}")
    c.eq(idx["nGenes"], len(hg), f"{key}: hvg nGenes vs genes.json")
    c.ok(0 < len(hg) <= 65535, f"{key}: hvg gene count {len(hg)}")
    c.eq(idx["chunkPx"], int(round(idx["chunkUm"] / L["micronsPerPixel"])), f"{key}: hvg chunkPx")
    h = L.get("hvg") or {}
    c.eq(h.get("nGenes"), len(hg), f"{key}: manifest hvg.nGenes")
    c.eq(h.get("chunks"), len(idx["chunks"]), f"{key}: manifest hvg.chunks")
    c.eq(h.get("bytes"), sum(ch["bytes"] for ch in idx["chunks"]), f"{key}: manifest hvg.bytes")
    c.eq(sorted(p.name for p in hd.glob("*.bin")), sorted(f"{ch['tx']}_{ch['ty']}.bin" for ch in idx["chunks"]),
         f"{key}: hvg chunk files vs index")
    px = idx["chunkPx"]
    seen = np.zeros(len(pos), np.int64)
    bad = []
    for ch in idx["chunks"]:
        buf = (hd / f"{ch['tx']}_{ch['ty']}.bin").read_bytes()
        if len(buf) != ch["bytes"]:
            bad.append(f"{ch['tx']}_{ch['ty']}: {len(buf)} bytes, index says {ch['bytes']}")
            continue
        nb = struct.unpack_from("<I", buf, 0)[0]
        o = 4
        rows = np.frombuffer(buf, "<u4", nb, o); o += 4 * nb
        ls = np.frombuffer(buf, "<f4", nb, o); o += 4 * nb
        ip = np.frombuffer(buf, "<u4", nb + 1, o); o += 4 * (nb + 1)
        nnz = int(ip[-1]) if nb else 0
        gi = np.frombuffer(buf, "<u2", nnz, o); o += 2 * nnz
        cnt = np.frombuffer(buf, "<u2", nnz, o); o += 2 * nnz
        if nb != ch["n"] or o != len(buf) or ip[0] != 0 or (np.diff(ip.astype(np.int64)) < 0).any():
            bad.append(f"{ch['tx']}_{ch['ty']}: layout (n {nb}/{ch['n']}, {o}/{len(buf)} bytes)")
            continue
        if (rows >= len(pos)).any() or (gi >= len(hg)).any():
            bad.append(f"{ch['tx']}_{ch['ty']}: row or gene index out of range")
            continue
        seen[rows] += 1
        tx = np.floor(pos[rows, 0].astype(np.float64) / px)
        ty = np.floor(pos[rows, 1].astype(np.float64) / px)
        if not ((tx == ch["tx"]) & (ty == ch["ty"])).all():
            bad.append(f"{ch['tx']}_{ch['ty']}: holds bins of another chunk")
        if not np.allclose(ls, lib[rows], rtol=1e-5, atol=1e-3):
            bad.append(f"{ch['tx']}_{ch['ty']}: libSize differs from the gene files")
        # spot-check counts against the per-gene files
        for k in range(0, nb, max(1, nb // 25)):
            for j in range(ip[k], ip[k + 1]):
                r_g = dense.get(hg[gi[j]])
                if r_g is None:
                    bad.append(f"hvg gene {hg[gi[j]]} has no gene file"); break
                hit = r_g[1][r_g[0] == rows[k]]
                if len(hit) != 1 or min(round(float(hit[0])), 65535) != int(cnt[j]):
                    bad.append(f"{ch['tx']}_{ch['ty']}: count of {hg[gi[j]]} at row {rows[k]} differs"); break
    c.ok(not bad, f"{key}: hvg chunks: {bad[:5]}")
    c.ok((seen == 1).all(), f"{key}: hvg chunks cover {int((seen == 1).sum())}/{len(pos)} rows exactly once")
    print(f"    hvg: {len(hg)} genes, {len(idx['chunks'])} chunks of {px} px")
    return len(hg)


def check_umap(c: Checker, d: Path, L: dict, pos, lib, truth, step_px) -> dict:
    key = L["key"]
    ud = d / "umap"
    if not c.ok((ud / "index.json").exists(), f"{key}: umap/index.json missing"):
        return {}
    idx = read_json(ud / "index.json")
    n = len(pos)
    for k in ("nBins", "nEmbedded", "minUmi", "bounds", "params", "versions", "clusterings", "default"):
        c.ok(k in idx, f"{key}: umap/index.json has no {k}")
    c.eq(idx["nBins"], n, f"{key}: umap nBins vs positions.bin rows")
    co = f32(ud / "coords.bin", 2)
    c.eq(len(co), n, f"{key}: coords.bin rows vs positions.bin rows")
    emb = ~np.isnan(co[:, 0])
    c.ok((np.isnan(co[:, 0]) == np.isnan(co[:, 1])).all(), f"{key}: coords.bin half-NaN rows")
    c.eq(int(emb.sum()), idx["nEmbedded"], f"{key}: embedded rows vs nEmbedded")
    # embedded = >= minUmi UMI, less the few bins with < minHvgDetected variable genes detected
    c.ok((emb <= (lib >= idx["minUmi"])).all(), f"{key}: an embedded bin has < {idx['minUmi']} UMI")
    out_few = int(((lib >= idx["minUmi"]) & ~emb).sum())
    c.ok(out_few <= 0.02 * max(1, int((lib >= idx["minUmi"]).sum())),
         f"{key}: {out_few} bins >= {idx['minUmi']} UMI left out for too few variable genes (expected a handful)")
    raw = (ud / "index.json").read_text(encoding="utf-8")
    c.ok("NaN" not in raw and "Infinity" not in raw, f"{key}: umap/index.json has NaN / Infinity (a browser cannot parse it)")
    b = idx["bounds"]
    c.ok((co[emb, 0] >= b[0] - 1e-3).all() and (co[emb, 0] <= b[2] + 1e-3).all()
         and (co[emb, 1] >= b[1] - 1e-3).all() and (co[emb, 1] <= b[3] + 1e-3).all(), f"{key}: coords outside bounds")
    u = L.get("umap") or {}
    c.eq(u.get("nEmbedded"), idx["nEmbedded"], f"{key}: manifest umap.nEmbedded")
    c.eq(u.get("default"), idx["default"], f"{key}: manifest umap.default")
    c.eq([(x["resolution"], x["n"]) for x in u.get("clusterings", [])],
         [(x["resolution"], x["n"]) for x in idx["clusterings"]], f"{key}: manifest umap.clusterings")
    c.ok(0 <= idx["default"] < len(idx["clusterings"]), f"{key}: umap default index")
    files = {"index.json", "coords.bin"} | {x["file"] for x in idx["clusterings"]}
    if (ud / "celltypes.json").exists():
        files |= {"celltypes.json", "celltype.bin"}
    c.eq(sorted(p.name for p in ud.iterdir()), sorted(files), f"{key}: files in umap/")

    pairs = grid_pairs(pos, step_px)
    rng = np.random.default_rng(0)
    out = {"nEmbedded": int(emb.sum()), "clusterings": []}
    truth = np.array(truth)
    for cl in idx["clusterings"]:
        p = ud / cl["file"]
        c.eq(p.stat().st_size, 2 * n, f"{key}: {cl['file']} bytes (uint16 x rows)")
        lab = np.fromfile(str(p), dtype="<u2")
        c.ok(((lab == NO_CLUSTER) == ~emb).all(), f"{key}: {cl['file']}: 0xFFFF rows != NaN coords rows")
        le = lab[emb].astype(np.int64)
        c.ok(le.max() < cl["n"], f"{key}: {cl['file']}: label >= n")
        sizes = np.bincount(le, minlength=cl["n"])
        c.eq(sizes.tolist(), cl["sizes"], f"{key}: {cl['file']}: sizes")
        c.ok((np.diff(sizes) <= 0).all(), f"{key}: {cl['file']}: cluster 0 is not the largest")
        c.eq(len(cl["markers"]), cl["n"], f"{key}: {cl['file']}: one marker list per cluster")
        # spatial coherence: neighbouring bins share a cluster far more often than chance
        pe = pairs[emb[pairs[:, 0]] & emb[pairs[:, 1]]]
        same = float((lab[pe[:, 0]] == lab[pe[:, 1]]).mean())
        sh = lab.copy()
        sh[emb] = rng.permutation(sh[emb])
        chance = float((sh[pe[:, 0]] == sh[pe[:, 1]]).mean())
        a = ari(le, truth[emb])
        c.ok(same > chance + 0.25, f"{key}: r={cl['resolution']}: neighbours agree {same:.2f} vs shuffled {chance:.2f}")
        print(f"    leiden r={cl['resolution']:g}: {cl['n']} clusters, neighbour agreement {same:.2f} "
              f"(shuffled {chance:.2f}), ARI vs drawn regions {a:.2f}")
        out["clusterings"].append({"res": cl["resolution"], "n": cl["n"], "same": same, "chance": chance, "ari": a})
    best = max((x["ari"] for x in out["clusterings"]), default=0)
    c.ok(best > 0.3, f"{key}: clusters do not recover the drawn regions (best ARI {best:.2f})")
    check_celltypes(c, ud, L, idx, emb, truth)
    return out


def check_celltypes(c: Checker, ud: Path, L: dict, idx: dict, emb, truth):
    """Suggested names from examples/demo_markers.json: each type is named after the region
    whose marker genes it lists, so a cluster drawn mostly from one region should be
    suggested that region's name."""
    key = L["key"]
    if not c.ok((ud / "celltypes.json").exists(), f"{key}: umap/celltypes.json missing"):
        return
    ct = read_json(ud / "celltypes.json")
    names = [t["name"] for t in ct["types"]]
    c.ok(all(t["scored"] for t in ct["types"]), f"{key}: a demo marker type was not scored")
    c.eq(((L.get("umap") or {}).get("celltypes") or {}).get("nTypes"), len(names), f"{key}: manifest umap.celltypes.nTypes")
    tb = np.fromfile(str(ud / "celltype.bin"), dtype=np.uint8)
    c.eq(len(tb), len(emb), f"{key}: celltype.bin rows vs positions.bin rows")
    c.ok(((tb == 255) == ~emb).all(), f"{key}: celltype.bin 255 rows != not-embedded rows")
    c.ok(((tb[emb] < len(names)) | (tb[emb] == 254)).all(), f"{key}: celltype.bin values out of range")
    typed = tb[emb] < len(names)
    agree = float(np.mean([names[v].lower() == r for v, r in zip(tb[emb][typed], truth[emb][typed])]))
    c.ok(agree > 0.8, f"{key}: bins' suggested type matches the drawn region for {agree:.0%} of typed bins")
    cl = idx["clusterings"][idx["default"]]
    lab = np.fromfile(str(ud / cl["file"]), dtype="<u2")[emb].astype(np.int64)
    sug = ct["clusterings"][f"{cl['resolution']:g}"]
    c.eq(len(sug), cl["n"], f"{key}: one suggestion per cluster")
    hits = total = 0
    for s in sug:
        regions, counts = np.unique(truth[emb][lab == s["cluster"]], return_counts=True)
        major = regions[counts.argmax()]
        if s["confidence"] in ("high", "medium"):
            total += 1
            hits += s["top"].lower() == major
    c.ok(total > 0 and hits == total, f"{key}: confident suggestions naming the cluster's main region: {hits}/{total}")
    print(f"    suggested names: {hits}/{total} confident suggestions name the cluster's main region; "
          f"{agree:.0%} of typed bins match their region")


def check_user(c: Checker, root: Path, clustering: bool, log_path: Path) -> dict:
    """The lab's own clusters / UMAP on the 8 um layer (userclusters.py); see the module docstring."""
    key = "008um"
    out = root / "export_user" / VIS
    m = read_json(out / "manifest.json")
    L = next((x for x in m.get("expr", []) if x["key"] == key), {})
    if not clustering:
        print("    SKIP user cluster checks: clustering packages not installed")
        return {}
    ud = out / "expr" / key / "umap"
    idx = read_json(ud / "index.json")
    raw = (ud / "index.json").read_text(encoding="utf-8")
    pos = f32(out / "expr" / key / "positions.bin", 2)
    n = len(pos)
    truth_rows = [ln.split(",") for ln in (root / "raw" / "visiumhd" / f"truth_square_{key}.csv")
                  .read_text(encoding="utf-8").split("\n")[1:] if ln]
    truth = np.array([r for _, r in truth_rows])
    in_file = np.arange(n) % 10 != 0                     # write_user_files dropped every tenth bin
    rep = idx.get("userImport", {}).get("barcodes", [])
    c.eq(len(rep), 2, "user: one barcode report per input file")
    by = {}
    for r in rep:
        by["umap" if r["methods"]["exact"] else "obs"] = r
    ro, ru = by.get("obs", {}), by.get("umap", {})
    c.eq((ro.get("inFile"), ro.get("matched"), ro.get("rate")), (int(in_file.sum()), int(in_file.sum()), 1.0),
         "user: clusters CSV barcodes matched")
    c.eq(ro.get("methods", {}).get("suffix"), int(in_file.sum()), "user: clusters CSV matched by the -1-0 suffix rule")
    c.eq((ru.get("matched"), ru.get("methods", {}).get("exact")), (int(in_file.sum()),) * 2, "user: UMAP CSV matched exactly")
    cls = idx["clusterings"]
    c.eq([(x.get("name"), x.get("source")) for x in cls],
         [("seurat_clusters", "user"), ("cell_type", "user"), ("Leiden 0.3", "soda"), ("Leiden 0.6", "soda"),
          ("Leiden 1", "soda")], "user: clusterings by name and source (yours first, nCount_Spatial skipped)")
    c.eq(idx["default"], 0, "user: default is the first user clustering")
    c.eq(idx.get("embedding", {}).get("source"), "user", "user: embedding source")
    c.ok("NaN" not in raw and "Infinity" not in raw, "user: index.json has NaN / Infinity")
    c.ok(not BARCODE.search(raw), "user: index.json contains a barcode")
    u = L.get("umap") or {}
    c.eq([(x.get("id"), x.get("resolution"), x.get("name"), x.get("source"), x.get("n")) for x in u.get("clusterings", [])],
         [(x.get("id"), x.get("resolution"), x.get("name"), x.get("source"), x.get("n")) for x in cls],
         "user: manifest umap.clusterings")
    c.eq((u.get("embedding"), u.get("default"), u.get("nEmbedded")), ("user", 0, idx["nEmbedded"]), "user: manifest umap entry")
    files = {"index.json", "coords.bin", "index_soda.json", "coords_soda.bin", "celltypes.json", "celltype.bin"} \
        | {x["file"] for x in cls}
    c.eq(sorted(p.name for p in ud.iterdir()), sorted(files), "user: files in umap/")
    base = read_json(root / "export" / VIS / "expr" / key / "umap" / "index.json")
    c.eq([(x["resolution"], x["n"]) for x in read_json(ud / "index_soda.json")["clusterings"]],
         [(x["resolution"], x["n"]) for x in base["clusterings"]], "user: index_soda.json = SODA's own clusterings")
    # coordinates: the CSV's, for exactly the bins in the file
    co = f32(ud / "coords.bin", 2)
    um = {ln.split(",")[0]: (float(ln.split(",")[1]), float(ln.split(",")[2]))
          for ln in (root / "user" / "umap.csv").read_text(encoding="utf-8").split("\n")[1:] if ln}
    want = np.full((n, 2), np.nan, np.float32)
    for i, (bc, _) in enumerate(truth_rows):
        if bc in um:
            want[i] = um[bc]
    c.ok(((co == want) | (np.isnan(co) & np.isnan(want))).all(), "user: coords.bin = the UMAP CSV by barcode")
    c.eq(idx["nEmbedded"], int(in_file.sum()), "user: nEmbedded")
    # cell_type: names, colours, labels per bin
    ct = next(x for x in cls if x["name"] == "cell_type")
    names = [USER_NAMES[r] for r in truth]
    c.eq(sorted(ct["categories"]), sorted(set(USER_NAMES.values())), "user: cell_type categories (the user's names)")
    c.eq(ct.get("colors"), [USER_COLORS[x] for x in ct["categories"]], "user: cell_type colours from the colour column")
    lab = np.fromfile(str(ud / ct["file"]), dtype="<u2")
    c.eq(len(lab), n, "user: clusters_user file rows")
    c.ok(((lab == NO_CLUSTER) == ~in_file).all(), "user: 0xFFFF exactly for the bins not in the file")
    got = np.array([ct["categories"][v] if v != NO_CLUSTER else "" for v in lab])
    agree = float((got[in_file] == np.array(names)[in_file]).mean())
    c.eq(agree, 1.0, "user: cell_type label per bin = the file's")
    c.eq(ct["sizes"], [int((got == x).sum()) for x in ct["categories"]], "user: cell_type sizes")
    c.eq(ct.get("nAssigned"), int(in_file.sum()), "user: cell_type nAssigned")
    sc = next(x for x in cls if x["name"] == "seurat_clusters")
    c.eq(sc["categories"], [str(k) for k in range(len(set(truth)))], "user: seurat_clusters categories 0..5 in order")
    # markers and suggested names per user group
    c.eq(len(ct["markers"]), ct["n"], "user: one marker list per group")
    fi = ct["categories"].index("Schwann cells")
    c.ok(any(mk["gene"].startswith("MARKER_FASC") for mk in ct["markers"][fi][:3]),
         f"user: Schwann cells' top markers {[mk['gene'] for mk in ct['markers'][fi][:3]]} lack MARKER_FASC*")
    cj = read_json(ud / "celltypes.json")
    c.eq(sorted(cj["clusterings"]), sorted([x["id"] for x in cls if x.get("source") == "user"] + ["0.3", "0.6", "1"]),
         "user: celltypes.json keyed by user id and resolution")
    sg = cj["clusterings"].get(ct["id"], [])
    c.eq(len(sg), ct["n"], "user: one suggestion per user group")
    c.ok(sg and sg[fi]["top"] == "Fascicle" and sg[fi]["confidence"] in ("high", "medium"),
         f"user: Schwann cells suggested as Fascicle ({sg[fi] if sg else None})")
    res = {"matchRate": ro.get("rate"), "method": ro.get("methods"), "groups": ct["n"], "agree": agree}

    # partial / none / merged, with `prepare.py userclusters` on a copy
    part = root / "export_partial" / VIS
    if part.exists():
        shutil.rmtree(part)
    shutil.copytree(out, part)
    pidx = part / "expr" / key / "umap" / "index.json"

    def uc(csv_name, *extra):
        args = ["web/SpatialViewer/prepare.py", "userclusters", "--id", VIS, "--out", str(root / "cache"),
                "--dest", str(part.parent), "--user-layer", key, "--user-obs", str(root / "user" / csv_name), *extra]
        p = subprocess.run([sys.executable, *args], cwd=str(REPO), stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                           text=True, encoding="utf-8", errors="replace")
        with open(log_path, "a", encoding="utf-8") as lg:
            lg.write(f"\n$ python {' '.join(args)}\n{p.stdout}")
        return p.returncode, p.stdout

    before = pidx.read_bytes()
    rc, txt = uc("nomatch.csv")
    c.ok(rc != 0 and "match this layer" in txt and "016" in txt, f"user: no match must fail loudly and name the layer (exit {rc})")
    c.ok(pidx.read_bytes() == before, "user: a failed import left index.json unchanged")
    rc, txt = uc("merged.csv")
    c.ok(rc != 0 and "sample tags" in txt, f"user: a merged two-sample file must fail and ask for a prefix (exit {rc})")
    rc, txt = uc("merged.csv", "--user-barcode-prefix", "A_")
    pi = read_json(pidx)
    c.ok(rc == 0 and pi["userImport"]["barcodes"][0].get("otherSample") == 2000,
         f"user: with --user-barcode-prefix the merged file imports one sample (exit {rc})")
    rc, txt = uc("partial.csv")
    pi = read_json(pidx)
    pr = pi["userImport"]["barcodes"][0] if pi.get("userImport") else {}
    exp = int(in_file.sum()) - sum(1 for j in range(int(in_file.sum())) if j % 5 < 2)
    c.ok(rc == 0 and "WARNING" in txt, f"user: a partial match imports with a warning (exit {rc})")
    c.eq((pr.get("matched"), round(pr.get("rate", 0), 2)), (exp, round(exp / int(in_file.sum()), 2)), "user: partial match report")
    c.eq([x.get("source") for x in pi["clusterings"]].count("soda"), 3, "user: SODA's clusterings kept next to the partial import")
    c.eq(pi.get("embedding", {}).get("source"), "soda", "user: clusters without a UMAP stay on SODA's UMAP")
    soda_co = np.fromfile(str(part / "expr" / key / "umap" / "coords_soda.bin"), dtype="<f4")
    c.ok(np.array_equal(np.fromfile(str(part / "expr" / key / "umap" / "coords.bin"), dtype="<f4"), soda_co, equal_nan=True),
         "user: coords.bin is SODA's when no UMAP is given")
    res["partialRate"] = pr.get("rate")
    print(f"    user clusters: {ro.get('matched')}/{ro.get('inFile')} barcodes matched ({ro.get('methods')}), "
          f"cell_type {ct['n']} groups, labels agree {agree:.0%}; partial {pr.get('rate')}")
    return res


def check_ids(c: Checker, root: Path) -> dict:
    """The barcode / cell id lists (layerids.py): in every default export (one per layer with barcodes, named in
    the manifest, the same as `prepare.py ids` writes), none with --no-ids; `prepare.py ids` on the Visium HD
    export_user copy (exported with --no-ids) and `export --ids` (the old flag, now asking explicitly) on the
    Xenium one write expr/<key>/ids.txt.gz with one id per positions.bin row in that order (= the cache's
    barcodes / the matrix's cell ids)."""
    out = {}
    for ds in (VIS, XE):
        m = read_json(root / "export" / ds / "manifest.json")
        rna = [L for L in m.get("expr", []) if L.get("kind") not in ("protein", "rois")]
        ok = bool(rna)
        for L in rna:
            f, e = root / "export" / ds / "expr" / L["key"] / "ids.txt.gz", L.get("ids") or {}
            pos = root / "export" / ds / "expr" / L["key"] / "positions.bin"
            n = len(gzip.decompress(f.read_bytes()).decode("utf-8").split("\n")) - 1 if f.exists() else -1
            ok = ok and e.get("file") == "ids.txt.gz" and n == e.get("n") == pos.stat().st_size // 8 and e.get("bytes") == f.stat().st_size
        c.ok(ok, f"{ds}: a default export has an id list for every layer ({[L['key'] for L in rna]}), one id per positions.bin row")
    # the same ids as prepare.py ids writes, byte for byte
    for key in ("008um", "016um"):
        a, b = root / "export" / VIS / "expr" / key / "ids.txt.gz", root / "export_user" / VIS / "expr" / key / "ids.txt.gz"
        c.ok(a.exists() and b.exists() and a.read_bytes() == b.read_bytes(), f"{VIS} {key}: the default export's id list = prepare.py ids's")
    lo = root / "export_loose" / VIS
    c.ok(not any("ids" in L for L in read_json(lo / "manifest.json").get("expr", [])) and not list(lo.rglob("ids.txt*")),
         f"{VIS} --no-ids: no id list and no manifest ids entry")
    # Xenium: the cell ids of the matrix, in its (= positions) order
    xo = root / "export_user" / XE
    L = next((x for x in read_json(xo / "manifest.json").get("expr", []) if x["key"] == "cells"), {})
    e = L.get("ids") or {}
    c.eq((e.get("file"), e.get("column")), ("ids.txt.gz", "cell_id"), f"{XE} --ids: manifest ids entry")
    f = xo / "expr" / "cells" / "ids.txt.gz"
    if c.ok(f.exists(), f"{XE} --ids: expr/cells/ids.txt.gz written"):
        got = gzip.decompress(f.read_bytes()).decode("utf-8").split("\n")[:-1]
        want = gzip.open(root / "raw" / "xenium" / "cell_feature_matrix" / "barcodes.tsv.gz", "rt").read().split()
        c.eq(got, want, f"{XE} --ids: one cell id per positions.bin row, in matrix order")
        c.eq((e.get("n"), e.get("bytes")), (len(f32(xo / "expr" / "cells" / "positions.bin", 2)), f.stat().st_size), f"{XE} --ids: n and bytes")
        out[f"{XE}/cells"] = {"rows": len(got), "bytes": f.stat().st_size, "bytesPerRow": round(f.stat().st_size / max(len(got), 1), 2)}
    # Visium HD: prepare.py ids on an existing export; the cache's barcodes of the rows with a position
    vo = root / "export_user" / VIS
    for L in read_json(vo / "manifest.json").get("expr", []):
        key, e = L["key"], L.get("ids") or {}
        f = vo / "expr" / key / "ids.txt.gz"
        if not c.ok(f.exists() and e.get("column") == "barcode", f"{VIS} ids: {key}/ids.txt.gz and its manifest entry"):
            continue
        got = gzip.decompress(f.read_bytes()).decode("utf-8").split("\n")[:-1]
        pos = f32(vo / "expr" / key / "positions.bin", 2)
        c.eq((len(got), e.get("n"), e.get("bytes")), (len(pos), len(pos), f.stat().st_size), f"{VIS} ids {key}: rows, n, bytes")
        c.eq(len(set(got)), len(got), f"{VIS} ids {key}: every barcode once")
        c.ok(all(BARCODE.match(b) for b in got[:50]), f"{VIS} ids {key}: Space Ranger bin barcodes ({got[:2]})")
        # a bin barcode names its grid row and column: they must follow the positions (same order)
        rc = np.array([[int(b.split("_")[2]), int(b.split("_")[3].split("-")[0])] for b in got], np.float64)
        ok = len(rc) > 3
        if ok:
            for k, col in ((0, 1), (1, 0)):          # x follows the column, y the row
                A = np.column_stack([rc[:, col], rc[:, 1 - col], np.ones(len(rc))])
                coef, res, *_ = np.linalg.lstsq(A, pos[:, k].astype(np.float64), rcond=None)
                ok = ok and float(np.abs(A @ coef - pos[:, k]).max()) < 1.0
        c.ok(ok, f"{VIS} ids {key}: the barcodes' grid row / column follow positions.bin (same order)")
        try:
            import zarr
            zg = zarr.open_group(str(root / "cache" / VIS / "expr" / f"{key}.zarr"), mode="r")
            xy = np.asarray(zg["xy"])
            want = [b.decode() for b, on in zip(zg["barcodes"][()], ~np.isnan(xy[:, 0])) if on]
            c.eq(got, want, f"{VIS} ids {key}: = the cache's barcodes of the bins with a position")
        except ImportError:
            print("    SKIP ids vs cache barcodes (zarr not installed)")
        out[f"{VIS}/{key}"] = {"rows": len(got), "bytes": f.stat().st_size, "bytesPerRow": round(f.stat().st_size / max(len(got), 1), 2)}
    print("    id lists: " + ", ".join(f"{k} {v['rows']:,} rows {v['bytes'] / 1e3:.1f} kB ({v['bytesPerRow']} B/row)" for k, v in out.items()))
    check_ids_gates(c, root)
    return out


# a made-up sample ID of the shape deid.py knows (never a real one)
FAKE_SAMPLE = "INT9907"   # deid:ok (fake)


def check_ids_gates(c: Checker, root: Path):
    """The id lists are published, so: (1) `deid.py check` reads ids.txt.gz and reports a sample ID planted in it
    by file and line only; (2) `prepare.py ids` refuses a layer whose ids do not fit the platform's format
    (a merge prefix, a sample name glued on), writes nothing for any layer, and never prints an id."""
    work = root / "ids_gates"
    if work.exists():
        shutil.rmtree(work)
    src = root / "export_user" / VIS / "expr" / "008um" / "ids.txt.gz"
    if not c.ok(src.exists(), "ids gates: the Visium HD id list to start from"):
        return
    # (1) a sample ID planted on line 7 of a copy
    lines = gzip.decompress(src.read_bytes()).decode("utf-8").split("\n")
    planted = FAKE_SAMPLE + "_" + lines[6]
    lines[6] = planted
    bad = work / "deid" / VIS / "expr" / "008um" / "ids.txt.gz"
    bad.parent.mkdir(parents=True)
    bad.write_bytes(gzip.compress("\n".join(lines).encode("utf-8")))
    p = subprocess.run([sys.executable, "web/SpatialViewer/deid.py", "check", str(work / "deid")], cwd=str(REPO),
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace")
    c.eq(p.returncode, 1, "deid.py check finds a sample ID planted inside an ids.txt.gz (exit code)")
    c.ok("ids.txt.gz:7" in p.stdout and "[deid] 1 sample ID(s) found" in p.stdout,
         f"deid.py check names the .gz file and line 7: {p.stdout.strip().splitlines()[:2]}")
    c.ok(FAKE_SAMPLE not in p.stdout and planted not in p.stdout, "deid.py check does not print the ID it found in a .gz file")
    # (2) a cache whose 016um barcodes carry a merge prefix (4) and a sample ID (1): prepare.py ids refuses
    try:
        import zarr
    except ImportError:
        print("    SKIP the ids refusal check (zarr not installed)")
        return
    cache, exp = work / "cache" / VIS, work / "export" / VIS
    (cache / "expr").mkdir(parents=True)
    shutil.copyfile(root / "cache" / VIS / "manifest.json", cache / "manifest.json")
    keys = [L["key"] for L in read_json(cache / "manifest.json").get("expr", [])]
    for k in keys:
        shutil.copytree(root / "cache" / VIS / "expr" / f"{k}.zarr", cache / "expr" / f"{k}.zarr")
        (exp / "expr" / k).mkdir(parents=True)
        shutil.copyfile(root / "export" / VIS / "expr" / k / "positions.bin", exp / "expr" / k / "positions.bin")
    shutil.copyfile(root / "export" / VIS / "manifest.json", exp / "manifest.json")
    g = zarr.open_group(str(cache / "expr" / "016um.zarr"), mode="r+")
    bcs = [b.decode() for b in g["barcodes"][()]]
    xy = np.asarray(g["xy"])
    on = np.flatnonzero(~np.isnan(xy[:, 0]))
    for j, r in enumerate(on[:5]):
        bcs[r] = (f"{FAKE_SAMPLE}_" if j == 0 else "S1_") + bcs[r]
    g.create_dataset("barcodes", data=np.array([b.encode() for b in bcs], dtype="S"), chunks=(1 << 20,), overwrite=True)
    before = (exp / "manifest.json").read_bytes()
    p = subprocess.run([sys.executable, "web/SpatialViewer/prepare.py", "ids", "--id", VIS, "--out", str(work / "cache"),
                        "--dest", str(work / "export")], cwd=str(REPO), stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                       text=True, encoding="utf-8", errors="replace")
    c.ok(p.returncode != 0, f"prepare.py ids refuses ids that do not fit the Visium HD format (exit {p.returncode})")
    c.ok("REFUSING to write 016um/ids.txt.gz: 5 of" in p.stdout and "visium_hd id format" in p.stdout,
         f"the refusal names the layer, the count and the format: {[ln for ln in p.stdout.splitlines() if 'REFUSING' in ln][:1]}")
    c.ok(FAKE_SAMPLE not in p.stdout and "S1_" not in p.stdout and not any(b in p.stdout for b in bcs[:50]),
         "the refusal prints no id (neither the bad ones nor the good ones)")
    c.ok(not list(exp.rglob("ids.txt*")) and (exp / "manifest.json").read_bytes() == before,
         "a refusal writes nothing: no id list for any layer (the good 008um included), manifest unchanged")
    # (3) the same cache exported by default (no flag): 016um's list is left out with a *** line that names the layer,
    # the reason and the count but no id; 008um keeps its list; the export goes on (exit 0)
    full = work / "cache_full" / VIS
    shutil.copytree(root / "cache" / VIS, full)
    shutil.rmtree(full / "expr" / "016um.zarr")
    shutil.copytree(cache / "expr" / "016um.zarr", full / "expr" / "016um.zarr")
    p = subprocess.run([sys.executable, "web/SpatialViewer/prepare.py", "export", "--id", VIS, "--out", str(work / "cache_full"),
                        "--dest", str(work / "export_default"), "--no-umap", "--no-hvg"], cwd=str(REPO), stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace")
    ed = work / "export_default" / VIS
    skip = [ln for ln in p.stdout.splitlines() if "LEFT OUT" in ln]
    c.ok(p.returncode == 0 and len(skip) == 1 and "*** 016um" in skip[0] and "5 of" in skip[0] and "visium_hd id format" in skip[0],
         f"a default export leaves out the failing layer's list with a *** line naming it, the reason and the count (exit {p.returncode}): {skip[:1]}")
    c.ok(FAKE_SAMPLE not in p.stdout and "S1_" not in p.stdout and not any(b in p.stdout for b in bcs[:50]),
         "the default export's message prints no id")
    m2 = read_json(ed / "manifest.json") if (ed / "manifest.json").exists() else {}
    got = {L["key"]: bool(L.get("ids")) for L in m2.get("expr", [])}
    c.ok(got.get("008um") is True and got.get("016um") is False and (ed / "expr" / "008um" / "ids.txt.gz").exists()
         and not (ed / "expr" / "016um" / "ids.txt.gz").exists(),
         f"a default export goes on: 008um keeps its list, 016um has none ({got})")
    # (4) a platform without an id format (AnnData / SpatialData obs names): any list passes unless it holds a sample ID
    # or a line break; the reason never names an id
    code = ("import sys, json; sys.path.insert(0, 'web/SpatialViewer'); import layerids as L; "
            "ok = ['cell_' + str(i) for i in range(50)] + ['AAACCCAAGG-1', 'my cell #3']; "
            f"print(json.dumps([L.problem(ok, None), L.problem(ok + ['{FAKE_SAMPLE}_cell'], None), L.problem(ok + ['a' + chr(10) + 'b'], 'anndata'), "
            "L.platform_of('anndata', {'key': 'cells'}), L.platform_of('spatialdata', {'key': 'cells'})]))")
    p = subprocess.run([sys.executable, "-c", code], cwd=str(REPO), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                       encoding="utf-8", errors="replace")
    try:
        r = json.loads(p.stdout.strip().splitlines()[-1])
    except (ValueError, IndexError):
        r = None
    c.ok(r is not None and r[0] is None and r[1] and r[1][0] == "sample" and r[2] and r[2][0] == "shape" and r[3] is None and r[4] is None
         and FAKE_SAMPLE not in p.stdout,
         f"ids of a platform with no format (AnnData, SpatialData): pass unless a sample ID or a line break is in them ({r})")
    shutil.rmtree(work)


def check_multiplex(c: Checker, root: Path, clustering: bool) -> dict:
    out = root / "export" / MX
    m = read_json(out / "manifest.json")
    check_manifest_common(c, m, MX, "jpg")
    c.eq(m.get("kind"), "tiff", "multiplex manifest kind")
    c.eq((m["image"]["width"], m["image"]["height"]), (2000, 2000), "multiplex image size")
    c.eq([ch["name"] for ch in m["image"]["channels"]], MX_ALL, "multiplex channel names (from the OME-XML)")
    c.ok(abs((m.get("pixelSizeUm") or 0) - 0.5) < 1e-6, f"multiplex pixelSizeUm {m.get('pixelSizeUm')}")
    check_tiles(c, out, m, "jpg")
    layers = m.get("expr") or []
    c.eq([L.get("key") for L in layers], [GRID_KEY], "multiplex expr layers (the protein grid only)")
    if not layers:
        return {}
    res = check_grid(c, out, layers[0], clustering)
    res["qc"] = check_qc(c, root, out, layers[0])
    if clustering:
        res["sample"] = check_sample(c, root, layers[0])
    return res


GRID_UM, GRID_KEY, MX_MPP = 16.0, "grid016um", 0.5
MX_ALL = ["DAPI", "Synth-Membrane", "Synth-Axon", "Synth-Vessel", "Synth-Immune", "Synth-Autofluor",
          "Synth-Saturated", "Synth-Axon-488", "open-cy5"]
# written to the grid layer: open-cy5 is dropped by name, Synth-Axon-488 as a weaker duplicate of Synth-Axon
MX_CHANNELS = [x for x in MX_ALL if x not in ("Synth-Axon-488", "open-cy5")]
SAMPLE = 600                                                    # --umap-sample of the second multiplex export


def demo_regions(pos_px: np.ndarray) -> np.ndarray:
    """Main drawn region of each grid square (5 x 5 samples), from examples/make_demo.py's layout."""
    spec = importlib.util.spec_from_file_location("make_demo", REPO / "examples" / "make_demo.py")
    md = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(md)
    scale = md.FIELD_UM / (md.MX_PX * md.MX_MPP)               # layout um per image um
    off = (np.arange(5) + 0.5) / 5 - 0.5
    dx, dy = np.meshgrid(off, off)
    sx = (pos_px[:, 0:1] * MX_MPP + dx.ravel()[None, :] * GRID_UM) * scale
    sy = (pos_px[:, 1:2] * MX_MPP + dy.ravel()[None, :] * GRID_UM) * scale
    lab = md.region_at(sx, sy)
    major = np.array([np.bincount(r, minlength=len(md.REGIONS)).argmax() for r in lab])
    return np.array([md.REGIONS[k] for k in major])


def check_grid(c: Checker, out: Path, L: dict, clustering: bool) -> dict:
    key = L["key"]
    d = out / "expr" / key
    step = GRID_UM / MX_MPP                                     # 32 px squares
    print(f"  expr {key} (protein grid):")
    for k in ("label", "kind", "binSizeUm", "nBins", "nGenes", "micronsPerPixel", "spotDiameterPx",
              "imageScale", "channels"):
        c.ok(k in L, f"{key}: manifest layer has no {k}")
    c.eq(L.get("kind"), "protein", f"{key}: kind")
    c.eq(L.get("binSizeUm"), GRID_UM, f"{key}: binSizeUm")
    c.ok(abs(L.get("spotDiameterPx", 0) - step) < 1e-6, f"{key}: spotDiameterPx {L.get('spotDiameterPx')}")
    c.ok(abs(L.get("micronsPerPixel", 0) - MX_MPP) < 1e-9, f"{key}: micronsPerPixel {L.get('micronsPerPixel')}")
    c.ok("DAPI" not in L.get("channels", []) and len(L.get("channels", [])) >= 4,
         f"{key}: embedded channels {L.get('channels')} (nuclear stain must be left out)")
    c.ok(not (d / "hvg").exists(), f"{key}: protein layer must not have hvg/")

    pos = f32(d / "positions.bin", 2)
    n = len(pos)
    c.eq(n, L["nBins"], f"{key}: positions.bin rows vs nBins")
    c.ok(np.isfinite(pos).all() and ((pos >= 0) & (pos <= 2000)).all(), f"{key}: positions NaN / outside the image")
    g = pos / step - 0.5
    c.ok(np.abs(g - np.rint(g)).max() < 1e-4, f"{key}: positions are not square centres of the {step:g} px grid")
    c.ok(len({p.tobytes() for p in pos}) == n, f"{key}: duplicate positions")
    total = round(2000 / step) ** 2
    c.ok(total * 0.08 < n < total, f"{key}: {n} tissue squares of {total} (mask kept too few / everything)")

    genes = read_json(d / "genes.json")
    c.eq([x["n"] for x in genes], MX_CHANNELS, f"{key}: genes.json = the channels, in order")
    c.eq(len(genes), L["nGenes"], f"{key}: genes.json length vs nGenes")
    V = {}
    bad = []
    for x in genes:
        a = gene_f32(d, x["n"])
        if a is None:
            bad.append(f"missing {x['n']}"); continue
        if len(a) != n or not np.array_equal(a[:, :2], pos):
            bad.append(f"{x['n']}: rows / xy differ from positions.bin"); continue
        if not (np.isfinite(a[:, 2]).all() and (a[:, 2] >= 0).all()):
            bad.append(f"{x['n']}: values NaN or < 0"); continue
        if abs(float(a[:, 2].mean()) - x["t"]) > 1e-3 * max(1.0, x["t"]):
            bad.append(f"{x['n']}: mean {a[:, 2].mean()} != t {x['t']}")
        V[x["n"]] = a[:, 2].astype(np.float64)
    c.ok(not bad, f"{key}: gene files: {bad[:5]}")
    c.eq(gene_count(d), len(genes), f"{key}: genes with a payload (pack entries / gene files)")

    # mat/: the embedded channels square by square, for circle -> UMAP (protgrid.write_mat_chunks)
    mi = read_json(d / "mat" / "index.json") if (d / "mat" / "index.json").exists() else None
    c.ok(mi is not None and isinstance(L.get("mat"), dict), f"{key}: no mat/index.json or manifest mat entry")
    if mi:
        ch = mi.get("channels", [])
        c.eq(ch, L.get("channels"), f"{key}: mat channels vs the embedded channels")
        c.eq(len(mi.get("chunks", [])), L["mat"].get("chunks"), f"{key}: mat chunk count vs manifest")
        K, sc, cp = len(ch), mi.get("scale"), mi.get("chunkPx")
        seen = np.zeros(n, np.int64)
        bad = []
        for e in mi.get("chunks", []):
            p = d / "mat" / f"{e['tx']}_{e['ty']}.bin"
            b = p.read_bytes() if p.exists() else b""
            if len(b) != e["bytes"]:
                bad.append(f"{p.name}: {len(b)} bytes, index says {e['bytes']}"); continue
            m = int(np.frombuffer(b[:4], "<u4")[0])
            rows = np.frombuffer(b[4:4 + 4 * m], "<u4")
            q = np.frombuffer(b[4 + 4 * m:], "<u2").reshape(m, K).astype(np.float64) / sc
            seen[rows] += 1
            if (np.floor(pos[rows] / cp) != [e["tx"], e["ty"]]).any():
                bad.append(f"{p.name}: a square outside its chunk")
            for j, nm in enumerate(ch):
                if nm in V and np.abs(q[:, j] - V[nm][rows]).max() > 0.5 / sc + 1e-6:
                    bad.append(f"{p.name}: {nm} values differ from the gene's [x, y, v]"); break
        c.ok(not bad, f"{key}: mat chunks: {bad[:5]}")
        c.ok((seen == 1).all(), f"{key}: mat chunks must hold every square exactly once")
        c.ok(all(len(mi.get(k, [])) == K for k in ("bg", "hi", "median")), f"{key}: mat bg / hi / median per channel")
        # slide-wide mean / sd of x and mean raw value: what a shape is compared with (viewer lassoRoiCompare)
        if c.ok(all(len(mi.get(k, [])) == K for k in ("mean", "sd", "rawMean")) and mi.get("nSquares") == n,
                f"{key}: mat mean / sd / rawMean per channel and nSquares"):
            vq = np.zeros((n, K))
            for e in mi.get("chunks", []):
                b = (d / "mat" / f"{e['tx']}_{e['ty']}.bin").read_bytes()
                m = int(np.frombuffer(b[:4], "<u4")[0])
                vq[np.frombuffer(b[4:4 + 4 * m], "<u4")] = np.frombuffer(b[4 + 4 * m:], "<u2").reshape(m, K) / sc
            cof, bgv, hiv = mi["cofactor"], np.array(mi["bg"]), np.array(mi["hi"])
            xq = np.where(hiv > 0, np.arcsinh(np.maximum(vq - bgv, 0) / np.where(hiv > 0, hiv, 1) / cof) / np.arcsinh(1 / cof), 0)
            c.ok(np.allclose(mi["mean"], xq.mean(0), atol=2e-6) and np.allclose(mi["sd"], xq.std(0), atol=2e-6)
                 and np.allclose(mi["rawMean"], vq.mean(0), atol=2e-6), f"{key}: mat mean / sd / rawMean differ from the chunks")
        c.eq(mi.get("normalize"), "bg-asinh-v1", f"{key}: mat normalize")
        c.ok(all(h > 0 for h in mi.get("hi", [])), f"{key}: mat hi must be > 0 for embedded channels")
    gr = L.get("grid") or {}
    # default tissue rule: Otsu after a 100 µm (7 x 7 squares of 16 µm = 112 µm) moving average of the DAPI means
    c.eq((gr.get("tissueSmoothUm"), gr.get("tissueSmoothSquares")), (112.0, 7), f"{key}: grid tissue smoothing")
    c.eq(gr.get("level"), 0, f"{key}: square means from full resolution (level 0)")
    c.ok("moving average" in gr.get("tissueRule", ""), f"{key}: grid tissueRule {gr.get('tissueRule')!r}")
    if "DAPI" in V and "tissueThreshold" in gr and gr.get("tissueSmoothSquares", 1) <= 1:
        # raw rule only: with smoothing the threshold applies to the averaged means, not to each square's
        below = int((V["DAPI"] <= gr["tissueThreshold"]).sum())
        c.ok(below <= gr.get("holesFilled", 0), f"{key}: {below} kept squares have DAPI <= the tissue threshold "
                                                 f"(only the {gr.get('holesFilled', 0)} filled holes may)")

    truth = demo_regions(pos)
    if "Synth-Axon" in V:
        fa = truth == "fascicle"
        # the demo's fold (every channel bright) is not tissue to compare with: left out
        fold = (np.fromfile(str(d / "qc_flags.bin"), dtype=np.uint8) & 1) > 0 if (d / "qc_flags.bin").exists() \
            else np.zeros(n, bool)
        ratio = V["Synth-Axon"][fa & ~fold].mean() / max(V["Synth-Axon"][~fa & ~fold].mean(), 1e-9)
        c.ok(fa.any() and ratio > 3, f"{key}: Synth-Axon only {ratio:.1f}x higher in fascicle squares")
    print(f"    positions.bin {n:,} squares, {len(genes)} gene files; regions: "
          + ", ".join(f"{r} {k}" for r, k in zip(*np.unique(truth, return_counts=True))))
    if not clustering:
        c.ok("umap" not in L, f"{key}: umap in manifest although clustering packages are missing")
        print("    SKIP umap checks: clustering packages not installed")
        return {"squares": n}

    ud = d / "umap"
    if not c.ok((ud / "index.json").exists(), f"{key}: umap/index.json missing"):
        return {"squares": n}
    idx = read_json(ud / "index.json")
    c.eq(idx.get("nBins"), n, f"{key}: umap nBins vs squares")
    raw = (ud / "index.json").read_text(encoding="utf-8")
    c.ok("NaN" not in raw and "Infinity" not in raw, f"{key}: umap/index.json has NaN / Infinity (a browser cannot parse it)")
    c.eq(idx.get("modality"), "protein", f"{key}: umap modality")
    c.eq(idx["params"].get("channels"), L["channels"], f"{key}: umap params.channels vs manifest channels")
    co = f32(ud / "coords.bin", 2)
    # every square is embedded except those at background in every embedded channel
    c.ok(len(co) == n and (np.isnan(co[:, 0]) == np.isnan(co[:, 1])).all(), f"{key}: coords.bin rows / half-NaN rows")
    pe = np.isfinite(co[:, 0])
    c.eq(int(pe.sum()), idx.get("nEmbedded"), f"{key}: embedded squares vs nEmbedded")
    c.ok(idx.get("nEmbedded", 0) >= 0.5 * n, f"{key}: only {idx.get('nEmbedded')} of {n} squares embedded")
    for cl in idx["clusterings"]:
        lab = np.fromfile(str(ud / cl["file"]), dtype="<u2")
        c.ok(((lab == NO_CLUSTER) == ~pe).all(), f"{key}: {cl['file']}: 0xFFFF rows != not-embedded squares")
    u = L.get("umap") or {}
    c.eq([(x["resolution"], x["n"]) for x in u.get("clusterings", [])],
         [(x["resolution"], x["n"]) for x in idx["clusterings"]], f"{key}: manifest umap.clusterings")
    c.eq(sorted(p.name for p in ud.iterdir()),
         sorted({"index.json", "coords.bin", "celltypes.json", "celltype.bin"} | {x["file"] for x in idx["clusterings"]}),
         f"{key}: files in umap/")
    pairs = grid_pairs(pos, step)
    rng = np.random.default_rng(0)
    res = {"squares": n, "clusterings": []}
    for cl in idx["clusterings"]:
        lab = np.fromfile(str(ud / cl["file"]), dtype="<u2").astype(np.int64)
        c.eq(len(lab), n, f"{key}: {cl['file']} rows")
        c.ok((lab[pe] < cl["n"]).all() and (lab[~pe] == NO_CLUSTER).all(),
             f"{key}: {cl['file']}: label >= n, or 0xFFFF on an embedded square")
        sizes = np.bincount(lab[pe], minlength=cl["n"])
        c.eq(sizes.tolist(), cl["sizes"], f"{key}: {cl['file']}: sizes")
        c.ok((np.diff(sizes) <= 0).all(), f"{key}: {cl['file']}: cluster 0 is not the largest")
        c.eq(len(cl["markers"]), cl["n"], f"{key}: {cl['file']}: one marker list per cluster")
        c.ok(all(m["gene"] in V for ms in cl["markers"] for m in ms), f"{key}: {cl['file']}: marker not a channel")
        pp = pairs[pe[pairs[:, 0]] & pe[pairs[:, 1]]]
        same = float((lab[pp[:, 0]] == lab[pp[:, 1]]).mean())
        sh = lab.copy()
        sh[pe] = rng.permutation(sh[pe])
        chance = float((sh[pp[:, 0]] == sh[pp[:, 1]]).mean())
        a = ari(lab[pe], np.asarray(truth)[pe])
        # ~1,300 squares in 13-18 clusters at r=1: chance is tiny, so also accept 3x chance
        c.ok(same > chance + 0.25 or same > 3 * chance,
             f"{key}: r={cl['resolution']}: neighbours agree {same:.2f} vs shuffled {chance:.2f}")
        print(f"    leiden r={cl['resolution']:g}: {cl['n']} clusters, neighbour agreement {same:.2f} "
              f"(shuffled {chance:.2f}), ARI vs drawn regions {a:.2f}")
        res["clusterings"].append({"res": cl["resolution"], "n": cl["n"], "same": same, "chance": chance, "ari": a})
    best = max((x["ari"] for x in res["clusterings"]), default=0)
    c.ok(best > 0.3, f"{key}: clusters do not recover the drawn regions (best ARI {best:.2f})")

    # suggested names (protein rule, examples/demo_protein_markers.json)
    ct = read_json(ud / "celltypes.json")
    names = [t["name"] for t in ct["types"]]
    c.eq(ct.get("modality"), "protein", f"{key}: celltypes modality")
    c.ok("protein" in ct.get("rule", "").lower(), f"{key}: celltypes rule does not say protein")
    c.ok(all(t["scored"] for t in ct["types"]), f"{key}: a demo protein type was not scored")
    c.eq((u.get("celltypes") or {}).get("nTypes"), len(names), f"{key}: manifest umap.celltypes.nTypes")
    tb = np.fromfile(str(ud / "celltype.bin"), dtype=np.uint8)
    c.eq(len(tb), n, f"{key}: celltype.bin rows")
    c.ok(((tb < len(names)) | (tb == 254) | ((tb == 255) & ~pe)).all(), f"{key}: celltype.bin values out of range (255 only for squares not embedded)")
    typed = tb < len(names)
    agree = float(np.mean([names[v].lower() == r for v, r in zip(tb[typed], truth[typed])])) if typed.any() else 0.0
    c.ok(agree > 0.6, f"{key}: squares' suggested type matches the drawn region for {agree:.0%} of typed squares")
    cl = idx["clusterings"][idx["default"]]
    lab = np.fromfile(str(ud / cl["file"]), dtype="<u2").astype(np.int64)
    sug = ct["clusterings"][f"{cl['resolution']:g}"]
    c.eq(len(sug), cl["n"], f"{key}: one suggestion per cluster")
    # Vessels and the perineurium rim are about one 20 um square wide, so their squares are
    # mixed: a confident name counts when its region is >= 20% of the cluster's squares.
    # Squares on glass ("background") are left out of that share: the tissue mask averages the
    # DAPI means over 100 um (5 x 5 squares) first, so it keeps a 1-2 square fringe of glass
    # around the tissue - on this 1 mm demo ~140 squares - and a cluster of fringe squares has
    # no drawn region to be right about (the typed-square check above still counts them).
    hits = total = 0
    for s in sug:
        t = truth[lab == s["cluster"]]
        t = t[t != "background"]
        share = float(np.mean(t == s["top"].lower())) if len(t) else 0.0
        if s["confidence"] in ("high", "medium"):
            total += 1
            hits += share >= 0.2
    c.ok(total > 0 and hits == total, f"{key}: confident suggestions whose region is >= 20% of the cluster's "
                                      f"tissue squares: {hits}/{total}")
    print(f"    suggested names: {hits}/{total} confident suggestions name a region that is >= 20% of the cluster's tissue squares; "
          f"{agree:.0%} of typed squares match their region")
    return res


def check_qc(c: Checker, root: Path, out: Path, L: dict) -> dict:
    """Slide QC (protgrid.py qc.json / qc_flags.bin) on the demo's built-in problems."""
    key = L["key"]
    d = out / "expr" / key
    print(f"  expr {key} slide QC:")
    qc = read_json(d / "qc.json") if (d / "qc.json").exists() else None
    if not c.ok(qc is not None and isinstance(L.get("qc"), dict), f"{key}: no qc.json or manifest qc entry"):
        return {}
    n = L["nBins"]
    dropped = {x["name"]: x["reason"] for x in L.get("dropped", [])}
    c.ok("open" in dropped.get("open-cy5", ""), f"{key}: open-cy5 not dropped by name: {dropped.get('open-cy5')!r}")
    c.ok("duplicate" in dropped.get("Synth-Axon-488", ""),
         f"{key}: Synth-Axon-488 not dropped as a duplicate: {dropped.get('Synth-Axon-488')!r}")
    c.ok("Synth-Axon" not in dropped, f"{key}: the stronger Synth-Axon must be kept")
    ch = {x["name"]: x for x in qc.get("channels", [])}
    c.eq(list(ch), MX_CHANNELS, f"{key}: qc.json channels = the written channels")
    sat = ch.get("Synth-Saturated", {})
    c.ok(sat.get("saturated") and sat.get("saturatedPct", 0) > 1, f"{key}: Synth-Saturated not flagged saturated: {sat}")
    c.ok(sat.get("inUmap"), f"{key}: a saturated channel stays in the UMAP")
    others = [nm for nm, x in ch.items() if x.get("saturated") and nm != "Synth-Saturated"]
    c.ok(not others, f"{key}: channels wrongly flagged saturated: {others}")
    c.eq(L["qc"].get("saturated"), ["Synth-Saturated"], f"{key}: manifest qc.saturated")
    c.eq(qc.get("rules", {}).get("limit"), [65535], f"{key}: detector limit of the uint16 demo")
    c.ok(all(x["background"] >= 0 and x["p99_5AboveBg"] >= 0 for x in ch.values()), f"{key}: qc background / p99.5 < 0")
    fl = np.fromfile(str(d / "qc_flags.bin"), dtype=np.uint8)
    c.eq(len(fl), n, f"{key}: qc_flags.bin rows")
    c.eq(int(((fl & 1) > 0).sum()), qc["everythingBright"]["n"], f"{key}: everything-bright count vs qc_flags.bin")
    fold = read_json(root / "raw" / "truth_multiplex_fold.json")
    pos = f32(d / "positions.bin", 2)
    s = fold["sizePx"]
    inside = ((pos[:, 0] > fold["x0"]) & (pos[:, 0] < fold["x0"] + s) & (pos[:, 1] > fold["y0"]) & (pos[:, 1] < fold["y0"] + s))
    core = ((pos[:, 0] > fold["x0"] + 32) & (pos[:, 0] < fold["x0"] + s - 32)
            & (pos[:, 1] > fold["y0"] + 32) & (pos[:, 1] < fold["y0"] + s - 32))
    br = (fl & 1) > 0
    hit = float(br[core].mean()) if core.any() else 0.0
    out_rate = float(br[~inside].mean())
    c.ok(core.sum() >= 20 and hit > 0.8, f"{key}: fold squares flagged everything bright: {hit:.0%} of {int(core.sum())}")
    c.ok(out_rate < 0.02, f"{key}: {out_rate:.1%} of squares outside the fold flagged everything bright")
    c.ok(isinstance(qc.get("seams"), dict) and "flagged" in qc["seams"], f"{key}: qc.json seams check missing")
    c.ok(not qc["seams"]["flagged"], f"{key}: the demo has no camera tiles, but seams were flagged: {qc['seams']}")
    c.ok("correction" not in qc["seams"] and not (d / "qc_seams.json").exists(), f"{key}: untiled demo was 'corrected'")
    c.ok("not comparable across slides" in qc.get("note", ""), f"{key}: qc.json note")
    res = {"foldFlagged": round(hit, 3), "outsideFlagged": round(out_rate, 4), "saturatedPct": sat.get("saturatedPct")}
    ud = d / "umap"
    if (ud / "celltypes.json").exists():
        ct = read_json(ud / "celltypes.json")
        idx = read_json(ud / "index.json")
        art = []
        for cl in idx["clusterings"]:
            lab = np.fromfile(str(ud / cl["file"]), dtype="<u2").astype(np.int64)
            for sg in ct["clusterings"][f"{cl['resolution']:g}"]:
                m = lab == sg["cluster"]
                share = float(br[m].mean()) if m.any() else 0.0
                if abs(share - sg.get("artefact", -1)) > 1e-3:
                    art.append(f"r{cl['resolution']} c{sg['cluster']}: artefact {sg.get('artefact')} vs {share:.3f}")
                if (share > 0.5) != bool(sg.get("artefactRule")):
                    art.append(f"r{cl['resolution']} c{sg['cluster']}: artefactRule {sg.get('artefactRule')} at share {share:.2f}")
        c.ok(not art, f"{key}: artefact rule: {art[:4]}")
        named = [(cl, sg["cluster"]) for cl, v in ct["clusterings"].items() for sg in v if sg.get("artefactRule")]
        c.ok(named, f"{key}: no cluster suggested as an artefact (the fold should make one)")
        res["artefactClusters"] = len(named)
    print(f"    dropped {sorted(dropped)}; saturated {L['qc'].get('saturated')}; fold squares flagged {hit:.0%}, "
          f"outside {out_rate:.1%}; artefact clusters {res.get('artefactClusters')}")
    return res


def check_tiled(c: Checker, root: Path) -> dict:
    """demo_tiled: the grid found, the two stepped markers corrected, the factors recovered."""
    out = root / "export" / TL
    m = read_json(out / "manifest.json")
    L = (m.get("expr") or [{}])[0]
    d = out / "expr" / L.get("key", GRID_KEY)
    truth = read_json(root / "raw" / "truth_multiplex_tiles.json")
    qc = read_json(d / "qc.json")
    sm = qc.get("seams") or {}
    gr = sm.get("grid") or {}
    c.ok(sm.get("flagged") and gr, f"demo_tiled: tile grid not found: {sm.get('gridProblem')} {sm.get('columns')} {sm.get('rows')}")
    if not gr:
        return {}
    w, h = truth["w"], truth["h"]
    c.eq((gr["w"], gr["h"]), (w, h), "demo_tiled: tile size")
    off = lambda a, b, p: min(abs(a - b) % p, p - abs(a - b) % p)
    c.ok(off(gr["x0"], truth["x0"], w) <= 2 and off(gr["y0"], truth["y0"], h) <= 2,
         f"demo_tiled: tile offset {gr['x0']}, {gr['y0']} vs {truth['x0']}, {truth['y0']}")
    cor = sm.get("correction") or {}
    c.eq(sorted(cor.get("corrected", [])), sorted(truth["channels"]), "demo_tiled: corrected markers")
    fx = read_json(d / "qc_seams.json")
    chans = {x["name"]: x for x in fx["channels"]}
    res = {}
    for nm, t in truth["channels"].items():
        x = chans.get(nm, {})
        c.ok(x.get("excessAfter", 99) <= 0.5 * x.get("excessBefore", 0),
             f"demo_tiled: {nm} excess step {x.get('excessBefore')} -> {x.get('excessAfter')}")
        f = fx["factors"].get(nm)
        if not c.ok(f is not None, f"demo_tiled: no factors for {nm}"):
            continue
        tl = fx["tiles"]
        g = np.array(f["g"]).reshape(tl["rows"], tl["cols"])
        fit = np.zeros(tl["rows"] * tl["cols"], bool)
        fit[f["fitted"]] = True                                 # tiles with tissue on their edges
        fit = fit.reshape(tl["rows"], tl["cols"])
        gt = np.array(t["gain"])
        # the same tiles in both: tile (row, col) = index + row0 / col0
        rr = np.arange(tl["rows"]) + tl["row0"] - truth["row0"]
        cc = np.arange(tl["cols"]) + tl["col0"] - truth["col0"]
        okr, okc = (rr >= 0) & (rr < gt.shape[0]), (cc >= 0) & (cc < gt.shape[1])
        g, fit, gt = g[np.ix_(okr, okc)], fit[np.ix_(okr, okc)], gt[np.ix_(rr[okr], cc[okc])]
        fitted = fit
        e = np.log(g[fitted]) - np.log(gt[fitted])
        e -= np.median(e)                                       # only ratios between tiles are defined
        res[nm] = {"tiles": int(fitted.sum()), "medianErrPct": round(float(np.median(np.abs(e)) * 100), 2),
                   "p90ErrPct": round(float(np.percentile(np.abs(e), 90) * 100), 2), "method": x.get("method")}
        if not t["offset"] or not np.any(t["offset"]):
            # a factor only: the factors themselves must come back (with an offset too, a smooth
            # marker cannot tell the two apart - seams.py - and only the image test below applies)
            c.ok(fitted.sum() >= 20 and res[nm]["medianErrPct"] <= 6,
                 f"demo_tiled: {nm} factors not recovered: {res[nm]}")
        # the image test: the corrected squares must match the same image without tiles
        # (demo_multiplex, same squares) far better than the tile differences put in
        pa = f32(root / "export" / MX / "expr" / GRID_KEY / "positions.bin", 2)
        pb = f32(d / "positions.bin", 2)
        if not c.ok(len(pa) == len(pb) and np.array_equal(pa, pb), "demo_tiled: squares differ from demo_multiplex"):
            continue
        plain = gene_f32(root / "export" / MX / "expr" / GRID_KEY, nm)[:, 2]
        fixed = gene_f32(d, nm)[:, 2]
        sel = (plain > np.median(plain)) & (fixed > 0)
        lr = np.log(fixed[sel] / plain[sel])
        rsd = float(1.4826 * np.median(np.abs(lr - np.median(lr))))
        tcol = ((pb[sel, 0] - truth["x0"]) // w).astype(int) - truth["col0"]
        trow = ((pb[sel, 1] - truth["y0"]) // h).astype(int) - truth["row0"]
        lg = np.log(np.array(t["gain"])[trow, tcol])
        inj = float(1.4826 * np.median(np.abs(lg - np.median(lg))))
        res[nm].update(residualSdPct=round(rsd * 100, 2), injectedSdPct=round(inj * 100, 2))
        # edges pin each tile only to its neighbours, so small errors add up across the slide: on
        # this demo the factors come back within ~5-6%, and the corrected squares keep less than
        # 55% of the tile variation put in (measured 0.19-0.52 while building, 2026-11-05)
        c.ok(rsd <= 0.55 * inj, f"demo_tiled: {nm} corrected squares still differ from the untiled image by "
                               f"{rsd:.1%} (tile factors put in: {inj:.1%})")
    others = [x["name"] for x in fx["channels"] if x["name"] not in truth["channels"] and x["method"] != "none"]
    c.ok(not others, f"demo_tiled: markers without steps were corrected: {others}")
    print(f"    grid {gr['w']} x {gr['h']} px at {gr['x0']}, {gr['y0']}; corrected {cor.get('corrected')}; "
          f"recovered factors: {res}")
    return res


def _monotone(src: np.ndarray, out: np.ndarray) -> bool:
    """out never goes down where src goes up: the order of the values is kept (ties allowed)."""
    o = np.argsort(src.ravel(), kind="stable")
    return bool(np.all(np.diff(out.ravel()[o].astype(np.int64)) >= 0))


def check_float_image(c: Checker, root: Path) -> dict:
    """The float32 fixture is scaled into the 16-bit cache (0 .. its 99.99th percentile -> 0 .. 65535,
    pyramid._float_scale) instead of being cast as-is (~66 grey levels); uint16 comes through unchanged."""
    try:
        import tifffile
        import zarr
        from PIL import Image  # noqa: F401  (_img_tile)
    except ImportError as e:
        print(f"    SKIP float / uint16 image checks ({e.name} not installed)")
        return {}
    res = {}
    raw = root / "floatimg"
    # ---- float32
    src = tifffile.imread(str(raw / "float32.ome.tif"))
    m = read_json(root / "export_float" / FL / "manifest.json")
    c.eq(m["image"]["dtype"], "uint16", f"{FL}: cache dtype")
    z0 = np.asarray(zarr.open_group(str(root / "cache" / FL / "image.zarr"), mode="r")["0"][:])
    c.eq(z0.shape, src.shape, f"{FL}: level 0 shape")
    ramp, cr = src[0], z0[0]
    pos = ramp[ramp > 0].astype(np.float64)
    top = float(np.percentile(pos, 99.99))
    want = np.rint(np.clip(ramp.astype(np.float64) * (65535.0 / top), 0, 65535)).astype(np.uint16)
    # the scale is taken from every pixel here (the fixture is far below pyramid.FLOAT_SAMPLE_SIDE)
    c.ok(np.array_equal(cr, want), f"{FL}: Ramp not 0 .. p99.99 ({top:.4f}) -> 0 .. 65535: "
                                   f"{int((cr != want).sum())} pixels differ, max {int(cr.max())}")
    t0 = _img_tile(root / "export_float" / FL, m, 0)
    lv_cache, lv_tile = len(np.unique(cr)), len(np.unique(t0))
    c.ok(lv_cache > 10000, f"{FL}: Ramp has {lv_cache} distinct levels in the cache (a plain cast keeps 67)")
    c.ok(lv_tile >= 200, f"{FL}: Ramp's exported tile has {lv_tile} distinct levels (a plain cast keeps <= 67)")
    c.ok(_monotone(ramp, cr), f"{FL}: Ramp's order of values not kept in the cache")
    c.ok(_monotone(ramp, t0), f"{FL}: Ramp's order of values not kept in the exported tile")
    odd, co = src[1], z0[1]
    bad = ~np.isfinite(odd) | (odd <= 0)
    c.ok(int(co[bad].max()) == 0, f"{FL}: NaN / +-inf / negative pixels not 0 in the cache (max {int(co[bad].max())})")
    c.ok(_monotone(odd[~bad], co[~bad]), f"{FL}: Odd's order of values not kept")
    c.ok(len(np.unique(co[~bad])) > 10000, f"{FL}: Odd has {len(np.unique(co[~bad]))} distinct levels")
    c.ok(len(np.unique(z0[2])) == 1, f"{FL}: the constant channel is not one value in the cache: {np.unique(z0[2])[:5]}")
    res[FL] = {"levelsCache": lv_cache, "levelsTile": lv_tile, "p9999": round(top, 4)}
    # ---- uint16: stored as it is, `bits` / `window` and the tiles by the same rule as before
    src = tifffile.imread(str(raw / "uint16.ome.tif"))
    m = read_json(root / "export_float" / U16 / "manifest.json")
    z0 = np.asarray(zarr.open_group(str(root / "cache" / U16 / "image.zarr"), mode="r")["0"][:])
    c.ok(z0.dtype == np.uint16 and np.array_equal(z0, src), f"{U16}: level 0 is not the input as it is")
    for k, ch in enumerate(m["image"]["channels"]):
        nz = src[k][src[k] > 0]
        lo, hi = float(np.percentile(nz, 1.0)), float(np.percentile(nz, 99.8))
        t = float(max(float(np.percentile(nz, 99.99)), hi, 1))
        c.ok(ch["bits"] == [0, t] and ch["window"] == [round(lo * 255.0 / t, 2), round(hi * 255.0 / t, 2)]
             and ch["dataMax"] == int(src[k].max()), f"{U16}: channel {k} bits / window / dataMax {ch}")
        tile = _img_tile(root / "export_float" / U16, m, k)
        exp = np.clip((src[k].astype(np.float32) - 0.0) * (255.0 / (t - 0.0)), 0, 255).astype(np.uint8)
        c.ok(np.array_equal(tile, exp), f"{U16}: channel {k} tile differs from the 16-bit rule in "
                                        f"{int((tile != exp).sum())} pixels")
    print(f"    float32 Ramp: {lv_cache} levels in the cache, {lv_tile} in the tile (p99.99 {top:.3f}); uint16 unchanged")
    return res


def check_sample(c: Checker, root: Path, L: dict) -> dict:
    """The second multiplex export: UMAP fitted on SAMPLE squares, the rest placed."""
    key = L["key"]
    a = root / "export" / MX / "expr" / key / "umap"
    b = root / "export_sample" / MX / "expr" / key / "umap"
    print(f"  expr {key} UMAP on a {SAMPLE}-square sample:")
    if not c.ok((b / "index.json").exists(), f"{key}: export_sample has no umap/index.json"):
        return {}
    ia, ib = read_json(a / "index.json"), read_json(b / "index.json")
    n = ia["nBins"]
    c.eq((ib.get("nBins"), ib.get("nEmbedded"), ib.get("fitOn")), (n, ia["nEmbedded"], SAMPLE), f"{key}: sample nBins / nEmbedded / fitOn")
    c.ok("fitted on" in ib["params"].get("fit", "") and ia["params"].get("fit") == "all rows",
         f"{key}: params.fit {ib['params'].get('fit')!r} / {ia['params'].get('fit')!r}")
    co = f32(b / "coords.bin", 2)
    c.ok(len(co) == n and int(np.isfinite(co).all(1).sum()) == ia["nEmbedded"], f"{key}: sample coords.bin rows / NaN")
    out = {}
    for ca, cb in zip(ia["clusterings"], ib["clusterings"]):
        la = np.fromfile(str(a / ca["file"]), dtype="<u2").astype(np.int64)
        lb = np.fromfile(str(b / cb["file"]), dtype="<u2").astype(np.int64)
        m = lb != NO_CLUSTER                         # squares at background everywhere: not embedded
        sizes = np.bincount(lb[m], minlength=cb["n"])
        c.ok(len(lb) == n and (lb[m] < cb["n"]).all() and sizes.tolist() == cb["sizes"] and (np.diff(sizes) <= 0).all()
             and (m == (la != NO_CLUSTER)).all(), f"{key}: sample {cb['file']}: labels / sizes / order")
        out[ca["resolution"]] = round(ari(la[m], lb[m]), 3)
    best = max(out.values(), default=0)
    c.ok(best > 0.4, f"{key}: clusters from the sample disagree with the full fit (ARI {out})")
    print(f"    fitOn {ib.get('fitOn')} of {n}; ARI vs fitting on all: {out}")
    return {"ari": out}


def _mol_files(c: Checker, out: Path, ds: str) -> tuple[dict, dict]:
    """molecules/index.json and every gene file of an export as packed uint32 (qx | qy << 16)."""
    md = out / "molecules"
    if not c.ok((md / "index.json").exists(), f"{ds}: molecules/index.json missing"):
        return {}, {}
    idx = read_json(md / "index.json")
    packed, bad = {}, []
    for g in idx["genes"]:
        b = gene_bytes(md, g["n"])
        if not g["c"]:
            if b is not None:
                bad.append(f"{g['n']}: a payload although c = 0")
            continue
        if b is None or len(b) != 4 * g["c"]:
            bad.append(f"{g['n']}: {len(b) if b is not None else 'no'} bytes, want {4 * g['c']}")
            continue
        packed[g["n"]] = np.frombuffer(b, dtype="<u4").copy()
    c.ok(not bad, f"{ds}: molecule files: {bad[:4]}")
    c.eq(gene_count(md), len(packed), f"{ds}: molecule genes with a payload (pack entries / gene files)")
    return idx, packed


def check_xenium(c: Checker, root: Path, clustering: bool) -> dict:
    """The synthetic Xenium folders: cells as an expression layer (positions = centroid um / pixel
    size in matrix order, controls left out, totals = matrix sums, binqc / hvg / umap as for bins),
    molecules (quality >= 20, controls and off-image ones left out, counts per gene, exact quantised
    positions, random order), alignment of cells and molecules with the DAPI nuclei, and the
    parquet / no-cells and no-transcripts variants."""
    raw = root / "raw"
    tr = read_json(raw / "truth_xenium.json")
    tz = np.load(raw / "truth_xenium.npz")
    px, genes, W, H = tr["pixelSizeUm"], tr["genes"], tr["width"], tr["height"]
    out = root / "export" / XE
    m = read_json(out / "manifest.json")
    check_manifest_common(c, m, XE, "png")
    c.eq(m.get("kind"), "xenium", "xenium manifest kind")
    c.eq((m["image"]["width"], m["image"]["height"]), (W, H), "xenium image size")
    c.eq([ch["name"] for ch in m["image"]["channels"]], ["DAPI", "Boundary"], "xenium channel names (from the OME block)")
    c.ok(abs((m.get("pixelSizeUm") or 0) - px) < 1e-9, f"xenium pixelSizeUm {m.get('pixelSizeUm')}")
    check_tiles(c, out, m, "png")
    layers = {L["key"]: L for L in m.get("expr", [])}
    c.eq(sorted(layers), ["cells"], "xenium expr layers (a cells layer, no protein grid beside it)")
    L = layers.get("cells", {})
    d = out / "expr" / "cells"
    c.eq(L.get("kind"), "cells", "cells layer kind")
    c.eq(L.get("nBins"), tr["nCells"], "cells layer nBins = cells in the matrix")
    c.eq(L.get("nGenes"), len(genes), "cells layer nGenes (controls left out)")
    c.eq(L.get("binSizeUm"), None, "cells layer binSizeUm")
    c.ok(abs((L.get("micronsPerPixel") or 0) - px) < 1e-9, "cells layer micronsPerPixel")
    c.ok(10 / px < (L.get("spotDiameterPx") or 0) < 20 / px, f"cells layer spotDiameterPx {L.get('spotDiameterPx')}")
    pos = f32(d / "positions.bin", 2)
    want = (tz["xy_um"] / px).astype(np.float32)
    c.ok(pos.shape == want.shape and np.allclose(pos, want, atol=1e-3),
         "cells positions.bin = centroid um / pixel size, in matrix order (cells.csv rows are in another order)")
    gj = read_json(d / "genes.json")
    c.eq([g["n"] for g in gj], genes, "cells genes.json: the genes in matrix order, control features left out")
    cnt = tz["counts"].astype(np.float64)
    c.ok(np.allclose([g["t"] for g in gj], cnt.sum(0)), "cells genes.json totals = count-matrix column sums")
    lib = cnt.sum(1)
    dense, bad = {}, []
    for j, g in enumerate(genes):
        rows = np.flatnonzero(cnt[:, j] > 0)
        a = gene_f32(d, g)
        a = np.zeros((0, 3), np.float32) if a is None else a
        if len(a) != len(rows) or not np.allclose(a[:, :2], pos[rows]) or not np.allclose(a[:, 2], cnt[rows, j]):
            bad.append(g)
        dense[g] = (rows, cnt[rows, j])
    c.ok(not bad, f"cells gene files differ from the count matrix: {bad[:4]}")
    c.ok(not any(gene_bytes(d, x) is not None for x in tr["controls"]), "cells: a control feature has a gene file")
    # the cache keeps each row's cell id (user clusters are matched to cells by it)
    try:
        import gzip
        import zarr
        want_bc = gzip.open(raw / "xenium" / "cell_feature_matrix" / "barcodes.tsv.gz", "rt").read().split()
        zc = zarr.open_group(str(root / "cache" / XE / "expr" / "cells.zarr"), mode="r")
        got_bc = [b.decode() for b in zc["barcodes"][()]] if "barcodes" in zc else None
        c.eq(got_bc, want_bc, "cells cache barcodes = cell ids in matrix (= positions) order")
    except ImportError:
        print("    SKIP cells barcodes check (zarr not installed)")
    check_hvg(c, d, L, pos, lib, dense)
    check_binqc(c, d, L, lib, dense, genes)
    res = {"cells": int(len(pos))}
    if clustering:
        truth = np.array(["typea", "typeb", "typec"])[tz["ctype"]].tolist()
        res["umap"] = check_umap(c, d, L, pos, lib, truth, 20 / px)
    else:
        print("    SKIP umap checks for the cells layer (clustering packages not installed)")

    # cells sit on DAPI nuclei: the full-resolution DAPI tile is bright at every centroid
    try:
        from PIL import Image
        T, z = m["image"]["tileSize"], m["image"]["levels"] - 1
        tiles, vals = {}, []
        for x, y in pos:
            k = (int(x // T), int(y // T))
            if k not in tiles:
                tiles[k] = np.asarray(Image.open(out / "tiles" / "0" / str(z) / str(k[0]) / f"{k[1]}.png"))
            vals.append(int(tiles[k][int(y) - k[1] * T, int(x) - k[0] * T]))
        vals = np.array(vals)
        bg = float(np.median(np.concatenate([t.ravel() for t in tiles.values()])))
        c.ok(np.median(vals) > max(3 * bg, 60) and (vals > 2 * bg).mean() > 0.97,
             f"cells on DAPI nuclei: median {np.median(vals):.0f} at the centroids vs {bg:.0f} overall")
        print(f"    DAPI at the {len(vals)} centroids: median {np.median(vals):.0f} (8-bit), image median {bg:.0f}")
    except ImportError:
        print("    SKIP DAPI alignment check (Pillow not installed)")

    # molecules
    idx, packed = _mol_files(c, out, XE)
    if idx:
        c.eq(idx.get("format"), "xy-uint16", "molecules format")
        c.eq({g["n"]: g["c"] for g in idx["genes"]}, tr["moleculesPerGene"],
             "molecules per gene = quality >= 20, on the image, controls left out")
        c.eq(idx.get("excluded", {}).get("controls"), tr["controlMolecules"], "molecules: control codewords counted out")
        c.eq(idx.get("excluded", {}).get("outsideImage"), tr["outsideImage"], "molecules: off-image ones counted out")
        c.eq(idx.get("lowQv"), tr["lowQuality"], "molecules: below quality 20")
        c.eq(idx.get("read"), tr["read"], "molecules: transcripts read")
        mb = m.get("molecules") or {}
        c.eq((mb.get("n"), mb.get("nGenes"), mb.get("qvMin")), (sum(tr["moleculesPerGene"].values()), len(genes), 20.0),
             "manifest molecules block")
        sx, sy = idx["scale"]
        c.ok(abs(sx - W / 65535) < 1e-9 and abs(sy - H / 65535) < 1e-9, "molecules scale = image size / 65535")
        mg, mxy = tz["mol_gene"], tz["mol_xy_um"].astype(np.float64) / px
        bad, near = [], []
        from_cells = {g: set(dense[g][0].tolist()) for g in genes}
        for j, g in enumerate(genes):
            sel = mg == j
            qx = np.minimum(np.rint(mxy[sel, 0] / sx), 65535).astype(np.uint32)
            qy = np.minimum(np.rint(mxy[sel, 1] / sy), 65535).astype(np.uint32)
            got = packed.get(g, np.zeros(0, np.uint32))
            if not np.array_equal(np.sort(got), np.sort(qx | (qy << 16))):
                bad.append(g)
                continue
            # random order: the first half is a fair sample (its mean x is the gene's mean x)
            gx = (got & 0xFFFF).astype(np.float64) * sx
            if len(gx) > 200 and abs(gx[: len(gx) // 2].mean() - gx.mean()) > 0.05 * W:
                bad.append(f"{g} (not shuffled)")
            # molecule dots fall on cells that show the gene in the cells layer
            mxy_g = np.column_stack([gx, (got >> 16).astype(np.float64) * sy])
            dd = ((mxy_g[:, None, :] - pos[None, :, :].astype(np.float64)) ** 2).sum(-1)
            nn = dd.argmin(1)
            ok = (np.sqrt(dd[np.arange(len(nn)), nn]) <= 5.5 / px) & np.array([r in from_cells[g] for r in nn])
            near.append(float(ok.mean()))
        c.ok(not bad, f"molecule files differ from the quantised truth: {bad[:4]}")
        c.ok(bool(near) and min(near) > 0.8, f"molecules on cells expressing the gene: worst gene {min(near or [0]):.0%}")
        res["molecules"] = {"n": int(sum(len(v) for v in packed.values())), "onExpressingCell": round(min(near or [0]), 3)}
        print(f"    molecules: {res['molecules']['n']:,} in {len(packed)} files; per gene {min(near or [0]):.0%}-"
              f"{max(near or [0]):.0%} lie within 5.5 um of a cell whose gene file has that gene")

    # transcripts.parquet, no cell files: the same molecules, no expression layer
    out_p = root / "export" / XP
    mp = read_json(out_p / "manifest.json")
    c.eq(mp.get("expr"), [], f"{XP}: no expression layer without cell files (and --grid-um 0)")
    ip, pp = _mol_files(c, out_p, XP)
    if ip and idx:
        c.eq({g["n"]: g["c"] for g in ip["genes"]}, tr["moleculesPerGene"], f"{XP}: molecules per gene from the parquet")
        empty = np.zeros(0, np.uint32)
        c.ok(all(np.array_equal(np.sort(pp.get(g, empty)), np.sort(packed.get(g, empty))) for g in genes),
             f"{XP}: parquet molecules = zarr molecules")
    # cells.parquet + cell_feature_matrix.h5, no transcripts: the same cells layer, no molecules
    out_n = root / "export" / XN
    mn = read_json(out_n / "manifest.json")
    c.ok("molecules" not in mn and not (out_n / "molecules").exists(), f"{XN}: no molecules without transcripts")
    c.eq([L["key"] for L in mn.get("expr", [])], ["cells"], f"{XN}: cells layer from cells.parquet + .h5")
    # (cells.csv.gz holds 4 decimals of a µm, cells.parquet the full value: equal to well under a pixel)
    c.ok(np.allclose(f32(out_n / "expr" / "cells" / "positions.bin", 2), pos, atol=1e-3),
         f"{XN}: positions.bin = the csv / mtx import")
    c.ok(all((gene_f32(out_n / "expr" / "cells", g) is None) == (gene_f32(d, g) is None) and
             (gene_f32(d, g) is None or np.allclose(gene_f32(out_n / "expr" / "cells", g), gene_f32(d, g), atol=1e-3))
             for g in genes), f"{XN}: gene files = the csv / mtx import")
    return res


def _zip_folder(src: Path, dst: Path, top: str, method: int) -> None:
    with zipfile.ZipFile(dst, "w", method) as z:
        for q in sorted(src.rglob("*")):
            if q.is_file():
                z.write(q, (top + "/" if top else "") + q.relative_to(src).as_posix())


def _sha(p: Path) -> str:
    import hashlib
    return hashlib.sha256(p.read_bytes()).hexdigest()


def write_xenium_inputs(root: Path):
    """The Xenium input fixtures (run right after make_xenium_demo.py): raw/xenium_noimg (the demo without
    morphology_focus/), raw/xenium_output.zip (the demo folder as xenium_output/outs/ in a deflated zip),
    raw/xenium_bundle.zip (the demo's files at the top of a stored zip); and under xenin/: the cells-only variant
    zipped, and zips that must be refused: a '..' member, an absolute member, a member holding more than it declares."""
    raw, d = root / "raw", root / "xenin"
    noimg = raw / "xenium_noimg"
    if noimg.exists():
        shutil.rmtree(noimg)
    shutil.copytree(raw / "xenium", noimg, ignore=shutil.ignore_patterns("morphology_focus"))
    _zip_folder(raw / "xenium", raw / "xenium_output.zip", "xenium_output/outs", zipfile.ZIP_DEFLATED)
    _zip_folder(raw / "xenium", raw / "xenium_bundle.zip", "", zipfile.ZIP_STORED)
    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True)
    _zip_folder(raw / "xenium_nt", d / "cells_only.zip", "run", zipfile.ZIP_DEFLATED)
    for name, extra in (("bad_slip.zip", "run/../../outside.txt"), ("bad_absolute.zip", "/outside/evil.txt")):
        _zip_folder(raw / "xenium", d / name, "run", zipfile.ZIP_DEFLATED)
        with zipfile.ZipFile(d / name, "a") as z:
            z.writestr(extra, "x")
    # cells.csv.gz declared at half its size, in the file list and in its local header
    bad = d / "bad_size.zip"
    _zip_folder(raw / "xenium", bad, "run", zipfile.ZIP_DEFLATED)
    b = bytearray(bad.read_bytes())
    with zipfile.ZipFile(bad) as z:
        info = next(i for i in z.infolist() if i.filename.endswith("cells.csv.gz"))
    cd = b.rfind(info.filename.encode()) - 46
    assert b[cd:cd + 4] == b"PK\x01\x02" and b[info.header_offset:info.header_offset + 4] == b"PK\x03\x04"
    struct.pack_into("<I", b, cd + 24, info.file_size // 2)
    struct.pack_into("<I", b, info.header_offset + 22, info.file_size // 2)
    bad.write_bytes(bytes(b))
    (d / "zip_sha.json").write_text(json.dumps({n: _sha(raw / n) for n in ("xenium_output.zip", "xenium_bundle.zip")}),
                                    encoding="utf-8")


# (name, zip under xenin/ or raw/, extra environment): prepare.py xenium runs into xenin/cache
XIN_RUNS = {"cells_only": ("xenin/cells_only.zip", {"SODA_ZIP_MEMORY_MB": "0"}), "cells_only_folder": ("raw/xenium_nt", {}),
            "slip": ("xenin/bad_slip.zip", {}), "absolute": ("xenin/bad_absolute.zip", {}), "size": ("xenin/bad_size.zip", {}),
            "oversize": ("raw/xenium_output.zip", {"SODA_ZIP_MAX_GB": "0.001"}),
            "disk": ("raw/xenium_output.zip", {"SODA_UNZIP_RESERVE_GB": "1e9"})}


def run_xenium_inputs(root: Path, log_path: Path) -> dict:
    """prepare.py xenium on the cells-only zip with every compressed member copied out, and on the zips that must be
    refused; each one's exit code and output is kept in xenin/runs.json for check_xenium_inputs."""
    d = root / "xenin"
    for sub in ("cache", "cache_folder"):
        if (d / sub).exists():
            shutil.rmtree(d / sub)
    res = {}
    for name, (src, env) in XIN_RUNS.items():
        ds = XN if name.startswith("cells_only") else f"x_{name}"
        args = ["web/SpatialViewer/prepare.py", "xenium", "--id", ds, "--src", str(root / src),
                "--out", str(d / ("cache_folder" if name == "cells_only_folder" else "cache"))]
        print(f"$ {' '.join(f'{k}={v}' for k, v in env.items())} python {' '.join(args)}", flush=True)
        p = subprocess.run([sys.executable, *args], cwd=str(REPO), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                           encoding="utf-8", errors="replace", env={**os.environ, "PYTHONIOENCODING": "utf-8", **env})
        with open(log_path, "a", encoding="utf-8") as lg:
            lg.write(f"\n$ {' '.join(f'{k}={v}' for k, v in env.items())} python {' '.join(args)}\n{p.stdout}")
        print("    " + (p.stdout.strip().splitlines() or ["(no output)"])[-1] + f"  (exit {p.returncode})", flush=True)
        res[name] = {"rc": p.returncode, "out": p.stdout[-4000:]}
    (d / "runs.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    return res


def check_xenium_inputs(c: Checker, root: Path) -> dict:
    """No stain images: cells and molecules on a transcript-density background in Xenium's own frame; zips read as
    they are give the folder's cache; unsafe or oversized zips are refused, and nothing is left behind."""
    raw, d = root / "raw", root / "xenin"
    tr = read_json(raw / "truth_xenium.json")
    tz = np.load(raw / "truth_xenium.npz")
    px, genes = tr["pixelSizeUm"], tr["genes"]
    res = {}

    # ---------------- no stain images
    out = root / "export" / XNI
    m = read_json(out / "manifest.json")
    check_manifest_common(c, m, XNI, "png")
    bg = 2 * px                                                    # about 1 µm: two 0.5 µm morphology pixels
    c.ok(abs((m.get("pixelSizeUm") or 0) - bg) < 1e-9, f"{XNI}: background pixel = 2 morphology pixels ({m.get('pixelSizeUm')})")
    # extent: the transcripts' 250 µm chunks (an upper bound) + 10 µm, from Xenium's 0, 0
    mxy = tz["mol_xy_um"].astype(np.float64)
    want_wh = tuple(int(math.ceil(((math.floor(mxy[:, k].max() / 250.0) + 1) * 250.0 + 10.0) / bg)) for k in (0, 1))
    c.eq((m["image"]["width"], m["image"]["height"]), want_wh, f"{XNI}: background size = transcript chunks + 10 um, origin 0, 0")
    cm = read_json(root / "cache" / XNI / "manifest.json")
    c.eq((cm.get("frame") or {}).get("xeniumPixelSizeUm"), px, f"{XNI}: the cache's frame keeps experiment.xenium's pixel size")
    layers = {L["key"]: L for L in m.get("expr", [])}
    c.eq(sorted(layers), ["cells"], f"{XNI}: a cells layer and no protein grid on the background")
    pos = f32(out / "expr" / "cells" / "positions.bin", 2)
    c.ok(pos.shape == tz["xy_um"].shape and np.allclose(pos, tz["xy_um"] / bg, atol=1e-3),
         f"{XNI}: cells positions.bin = centroid um / background pixel size, matrix order")
    c.ok(abs((layers.get("cells", {}).get("micronsPerPixel") or 0) - bg) < 1e-9, f"{XNI}: cells layer micronsPerPixel")
    ok_genes = [g for g in genes if (gene_f32(out / "expr" / "cells", g) is None) == (gene_f32(root / "export" / XE / "expr" / "cells", g) is None)]
    c.eq(len(ok_genes), len(genes), f"{XNI}: the same genes with cells as {XE}")
    c.ok(bool((out / "expr" / "cells" / "shapes").exists()) == bool((root / "export" / XE / "expr" / "cells" / "shapes").exists()),
         f"{XNI}: cell outlines as for {XE}")
    _img_background(c, root, XNI, m, mxy / bg)
    idx, packed = _mol_files(c, out, XNI)
    if idx:
        c.eq({g["n"]: g["c"] for g in idx["genes"]}, tr["moleculesPerGene"], f"{XNI}: molecules per gene (quality >= 20, controls out)")
        c.eq(idx.get("excluded", {}).get("outsideImage"), tr["outsideImage"], f"{XNI}: the two molecules at x < 0 left out")
        W, H = m["image"]["width"], m["image"]["height"]
        sx, sy = idx["scale"]
        c.ok(abs(sx - W / 65535) < 1e-9 and abs(sy - H / 65535) < 1e-9, f"{XNI}: molecules scale = background size / 65535")
        mg, bad = tz["mol_gene"], []
        for j, g in enumerate(genes):
            q = mxy[mg == j] / bg
            want = np.minimum(np.rint(q[:, 0] / sx), 65535).astype(np.uint32) | (np.minimum(np.rint(q[:, 1] / sy), 65535).astype(np.uint32) << 16)
            if not np.array_equal(np.sort(packed.get(g, np.zeros(0, np.uint32))), np.sort(want)):
                bad.append(g)
        c.ok(not bad, f"{XNI}: molecule positions = um / background pixel, quantised (differ: {bad[:4]})")
    res["noImage"] = {"size": [m["image"]["width"], m["image"]["height"]], "pixelSizeUm": m.get("pixelSizeUm")}

    # ---------------- zips: the folder's cache, nothing left behind, the zip unchanged
    for ds, what in ((XZ, "deflated zip of the output folder, two folders down"), (XZS, "stored bundle zip, files at the top")):
        _same_cache(c, root / "cache" / ds, root / "cache" / XE, f"{ds} ({what}) = the folder's cache", source_differs=True)
    c.ok((root / "export" / XZ / "manifest.json").exists(), f"{XZ}: exported")
    c.ok(not (root / "cache" / "_unzipped").exists(), "zips: the files copied out are removed after a prepare (cache/_unzipped gone)")
    sha = read_json(d / "zip_sha.json")
    c.ok(all(_sha(raw / n) == h for n, h in sha.items()), "zips: the zips themselves are unchanged")
    c.ok(not (raw / "xenium_output").exists() and not (raw / "xenium_bundle").exists(), "zips: nothing written next to the zips")
    log_txt = (root / "smoke_test.log").read_text(encoding="utf-8", errors="replace") if (root / "smoke_test.log").exists() else ""
    c.ok("transcripts.zarr.zip: read in place (stored in the zip" in log_txt,
         "zips: a stored transcripts.zarr.zip is read in place (no copy)")
    c.ok("which is in outs/" in log_txt and "which is at the top of the zip" in log_txt, "zips: the dataset folder found at any depth")

    runs = read_json(d / "runs.json") if (d / "runs.json").exists() else {}
    r = runs.get("cells_only", {})
    c.ok(r.get("rc") == 0 and "cell_feature_matrix.h5" in r.get("out", "") and "cells.parquet" in r.get("out", "")
         and "copying" in r.get("out", ""), "zip, SODA_ZIP_MEMORY_MB=0: the .h5 and the parquet copied out, not read in memory")
    _same_cache(c, d / "cache" / XN, d / "cache_folder" / XN, "cells-only zip (every compressed member copied out) = the folder's cache",
                source_differs=True)
    for name, say in (("slip", ["refused", "'..'"]), ("absolute", ["refused", "absolute path"]),
                      ("size", ["does not match the file list", "refused"]), ("oversize", ["SODA_ZIP_MAX_GB", "zip bomb"]),
                      ("disk", ["not enough free space", "SODA_UNZIP_RESERVE_GB"])):
        r = runs.get(name, {})
        last = (r.get("out", "").strip().splitlines() or [""])[-1]
        c.ok(r.get("rc") == 1 and all(w in last for w in say), f"zip {name}: refused with the reason (exit {r.get('rc')}: {last[:200]})")
        c.ok(not (d / "cache" / f"x_{name}" / "manifest.json").exists(), f"zip {name}: no dataset written")
        if name != "size":                                     # the size lie is only seen while reading: the image is built by then
            c.ok(not (d / "cache" / f"x_{name}").exists() and not (d / "cache" / "_unzipped" / f"x_{name}").exists(),
                 f"zip {name}: refused before anything was written")
    c.ok(not (root / "outside.txt").exists() and not (d / "outside.txt").exists() and not (root.parent / "outside.txt").exists(),
         "zip slip: nothing written outside")
    res["refusals"] = {k: v.get("rc") for k, v in runs.items()}
    print(f"    {XNI}: background {res['noImage']['size']} px of {res['noImage']['pixelSizeUm']} um; zips = folder; refusals "
          + ", ".join(f"{k} exit {v}" for k, v in res["refusals"].items()))
    return res


def check_xenium_user(c: Checker, root: Path, clustering: bool) -> dict:
    """The lab's own clusters attached to the Xenium demo's cells (write_xenium_user_files): matched
    by cell id through the cell_id column (and a Xenium Explorer groups CSV), no UMAP given so they sit on SODA's UMAP next to SODA's
    clusterings, labels / colours / sizes per cell, markers and suggested names per group, and no
    cell id anywhere in the export's text files."""
    key = "cells"
    out = root / "export_user" / XE
    if not clustering:
        print("    SKIP Xenium user cluster checks: clustering packages not installed")
        return {}
    m = read_json(out / "manifest.json")
    L = next((x for x in m.get("expr", []) if x["key"] == key), {})
    ud = out / "expr" / key / "umap"
    if not c.ok((ud / "index.json").exists(), f"{XE} user: umap/index.json missing"):
        return {}
    idx = read_json(ud / "index.json")
    raw = root / "raw"
    ids = gzip.open(raw / "xenium" / "cell_feature_matrix" / "barcodes.tsv.gz", "rt").read().split()
    ctype = np.load(raw / "truth_xenium.npz")["ctype"]
    n = len(f32(out / "expr" / key / "positions.bin", 2))
    c.eq(n, len(ids), f"{XE} user: positions.bin rows = cells")
    in_file = np.arange(n) % 9 != 4                       # write_xenium_user_files left out every ninth cell
    rep = idx.get("userImport", {}).get("barcodes", [])
    c.eq(len(rep), 2, f"{XE} user: one barcode report per file")
    r0, r1 = (rep + [{}, {}])[:2]
    c.eq((r0.get("inFile"), r0.get("matched"), r0.get("rate"), (r0.get("methods") or {}).get("exact")),
         (int(in_file.sum()), int(in_file.sum()), 1.0, int(in_file.sum())), f"{XE} user: cell ids matched exactly")
    c.eq((r1.get("inFile"), r1.get("matched"), r1.get("rate")), (n, n, 1.0), f"{XE} user: Xenium Explorer groups CSV matched")
    base = read_json(root / "export" / XE / "expr" / key / "umap" / "index.json")
    cls = idx.get("clusterings", [])
    c.eq([(x.get("name"), x.get("source")) for x in cls],
         [("leiden", "user"), ("cell_type", "user"), ("group", "user")]
         + [(f"Leiden {x['resolution']:g}", "soda") for x in base["clusterings"]],
         f"{XE} user: clusterings (yours first, the row-index column skipped, SODA's kept)")
    c.eq(idx.get("embedding", {}).get("source"), "soda", f"{XE} user: no UMAP given -> SODA's embedding")
    c.ok(np.array_equal(np.fromfile(str(ud / "coords.bin"), dtype="<f4"),
                        np.fromfile(str(root / "export" / XE / "expr" / key / "umap" / "coords.bin"), dtype="<f4"), equal_nan=True),
         f"{XE} user: coords.bin = SODA's whole-section UMAP")
    u = L.get("umap") or {}
    c.eq([(x.get("id"), x.get("name"), x.get("source")) for x in u.get("clusterings", [])],
         [(x.get("id"), x.get("name"), x.get("source")) for x in cls], f"{XE} user: manifest umap.clusterings")
    ct = next((x for x in cls if x.get("name") == "cell_type"), {})
    c.eq(ct.get("categories"), XE_USER_NAMES, f"{XE} user: cell_type categories")
    c.eq(ct.get("colors"), XE_USER_COLORS, f"{XE} user: cell_type colours from the colour column")
    lab = np.fromfile(str(ud / ct.get("file", "missing")), dtype="<u2") if ct else np.zeros(0, np.uint16)
    c.eq(len(lab), n, f"{XE} user: clusters_user file rows")
    agree = 0.0
    if len(lab) == n:
        c.ok(((lab == NO_CLUSTER) == ~in_file).all(), f"{XE} user: 0xFFFF exactly for the cells not in the file")
        agree = float((lab[in_file] == ctype[in_file]).mean())
        c.eq(agree, 1.0, f"{XE} user: cell_type per cell = the file's")
        c.eq(ct.get("sizes"), [int(((lab == k) & in_file).sum()) for k in range(3)], f"{XE} user: cell_type sizes")
    ld = next((x for x in cls if x.get("name") == "leiden"), {})
    c.eq(ld.get("categories"), ["0", "1", "2"], f"{XE} user: leiden categories")
    gp = next((x for x in cls if x.get("name") == "group"), {})
    c.eq((gp.get("categories"), gp.get("colors"), gp.get("nAssigned")), (["Band A", "Band B", "Band C"], XE_USER_COLORS, n),
         f"{XE} user: Xenium Explorer groups, their colours from the bare color column")
    # markers and suggested names per group: the left band is TypeA (SynA1-4), and so on
    mk = ct.get("markers") or [[]]
    c.ok(len(mk) == 3 and all(any(x["gene"].startswith(f"Syn{t}") for x in mk[i][:3]) for i, t in enumerate("ABC")),
         f"{XE} user: top markers per group {[[x['gene'] for x in g[:3]] for g in mk]}")
    sg = read_json(ud / "celltypes.json").get("clusterings", {}).get(ct.get("id"), [])
    c.eq([(s.get("top"), s.get("confidence")) for s in sg], [("TypeA", "high"), ("TypeB", "high"), ("TypeC", "high")],
         f"{XE} user: suggested names per group")
    # no cell id in any text file of this export (they are the user's row keys; only ids.txt.gz holds them)
    ids_set = set(ids)
    tok = re.compile(r"[a-p]{8}-1")
    leaked = [p.name for p in out.rglob("*") if p.is_file() and p.suffix in (".json", ".geojson", ".txt")
              and any(t in ids_set for t in tok.findall(p.read_text(encoding="utf-8", errors="ignore")))]
    c.ok(not leaked, f"{XE} user: cell ids in the export's text files: {leaked[:3]}")
    print(f"    Xenium user clusters: {r0.get('matched')}/{r0.get('inFile')} cell ids matched, "
          f"cell_type {ct.get('n')} groups, labels agree {agree:.0%}")
    return {"matchRate": r0.get("rate"), "groups": ct.get("n"), "agree": agree}


def check_geomx(c: Checker, root: Path) -> dict:
    """The synthetic GeoMx study (examples/make_geomx_demo.py): a "rois" layer with one row per AOI,
    target values = the probe collapse worked out by plain loops in the generator, negative probes
    kept out of the genes and used for the background / LOQ, positions = mask centroids inside the
    ROI polygons of two scans put side by side, the annotation columns as groupings with no
    embedding, GeoMx QC values in binqc/, every target in hvg/, nothing of the slide names or
    Sample_IDs in the export; and the same study without images on the drawn grid."""
    tr = read_json(root / "raw" / "truth_geomx.json")
    px, genes = tr["pixelSizeUm"], tr["genes"]
    (wa, ha), (wb, hb) = tr["scans"]["demoA"], tr["scans"]["demoB"]
    off = {"demoA": 0, "demoB": wa + tr["gap"]}
    aois = [a for a in tr["aois"] if a["scan"] in off]
    out = root / "export" / GX
    m = read_json(out / "manifest.json")
    check_manifest_common(c, m, GX, "png")
    c.eq(m.get("kind"), "geomx", "geomx manifest kind")
    c.eq((m["image"]["width"], m["image"]["height"]), (wa + tr["gap"] + wb, max(ha, hb)), "geomx image = the two scans side by side")
    c.eq([ch["name"] for ch in m["image"]["channels"]], ["SYTO 13 (FITC/525nm)", "Alexa 532 (Cy3/568nm)"],
         "geomx channel names (Fluor + Name of the OME block)")
    c.ok(abs((m.get("pixelSizeUm") or 0) - px) < 1e-9, f"geomx pixelSizeUm {m.get('pixelSizeUm')}")
    check_tiles(c, out, m, "png")
    gi = m.get("geomx") or {}
    c.eq((gi.get("nAoi"), gi.get("nRoi"), gi.get("nSlides"), gi.get("layout")), (len(aois), 12, 2, "scan"), "geomx block: AOIs, ROIs, slides, layout")
    c.eq([(p["targets"], p["probes"], p["negativeProbes"]) for p in gi.get("panels", [])], [(4, 9, 3), (18, 32, 6)], "geomx block: the two panels")
    c.eq((gi.get("dcc") or {}).get("used"), len(aois), "geomx block: DCC files used (the third scan's and the NTC's are not)")
    c.eq((gi.get("dcc") or {}).get("inFolder"), len(tr["aois"]) + 1, "geomx block: DCC files in the folder")
    c.eq([x.get("total") for x in gi.get("ntc", [])], [tr["ntcTotal"]], "geomx block: the No Template Control's total")
    c.eq([s["offsetPx"][0] for s in gi.get("scans", [])], [0, off["demoB"]], "geomx block: scan offsets")
    layers = {L["key"]: L for L in m.get("expr", [])}
    c.eq(sorted(layers), ["rois"], "geomx expr layers (one rois layer, no protein grid of the stains)")
    L = layers.get("rois", {})
    d = out / "expr" / "rois"
    c.eq((L.get("kind"), L.get("nBins"), L.get("nGenes"), L.get("binSizeUm")), ("rois", len(aois), len(genes), None), "rois layer kind / rows / targets")
    med_d = 0.6 * 2 * math.sqrt(float(np.median([a["area"] for a in aois])) / math.pi) / px
    c.ok(abs((L.get("spotDiameterPx") or 0) - med_d) < 0.05, f"rois spotDiameterPx {L.get('spotDiameterPx')} = 0.6 x the diameter of a circle of the median AOI area ({med_d:.2f})")

    # ---- positions: mask centroids, inside the ROI outlines; one outline per ROI ---------------
    pos = f32(d / "positions.bin", 2)
    c.eq(len(pos), len(aois), "rois positions.bin rows")
    gj = read_json(out / "overlays" / "rois.geojson")
    c.eq(len(gj["features"]), 12, "rois.geojson: one outline per ROI (the two segments of one ROI share theirs)")
    c.eq((m.get("overlays") or [{}])[0].get("file"), "rois.geojson", "the ROI outlines are the first overlay")
    c.eq(L.get("regionOverlay"), "rois.geojson", "the ROI outlines are the layer's region overlay")
    rid = np.fromfile(d / "region.bin", dtype="<u2")
    polys = {f["properties"].get("__region_id"): np.asarray(f["geometry"]["coordinates"][0], float) for f in gj["features"]}

    def inside(x, y, P):
        n, hit = len(P) - 1, False
        for i in range(n):
            (x1, y1), (x2, y2) = P[i], P[i + 1]
            if (y1 > y) != (y2 > y) and x < x1 + (x2 - x1) * (y - y1) / (y2 - y1):
                hit = not hit
        return hit
    c.ok(len(rid) == len(pos) and all(r in polys and inside(float(x), float(y), polys[r]) for (x, y), r in zip(pos, rid)),
         "every AOI position lies inside the outline region.bin links it to")
    want = np.array([[a["mask"]["cx"] + off[a["scan"]], a["mask"]["cy"]] for a in aois])
    moved = [i for i, a in enumerate(aois) if (a["scan"], a["roi"]) == ("demoA", "005")]       # the "C": centroid outside
    same = [i for i in range(len(aois)) if i not in moved]
    c.ok(np.allclose(pos[same], want[same], atol=0.02), "AOI positions = centroid of the AOI's own mask (+ the scan's offset), in annotation order")
    c.ok(len(moved) == 1 and not np.allclose(pos[moved], want[moved], atol=1.0),
         "an AOI whose mask centroid falls outside its ROI is moved inside it")
    seg = [i for i, a in enumerate(aois) if (a["scan"], a["roi"]) == ("demoA", "003")]
    c.ok(len(seg) == 2 and rid[seg[0]] == rid[seg[1]] and abs(float(pos[seg[0], 0] - pos[seg[1], 0])) > 50,
         "two segments of one ROI: two rows, one outline, two positions")
    c.eq(len({(round(float(x), 2), round(float(y), 2)) for x, y in pos}), len(pos), "no two AOIs share a position")
    c.eq(sorted(r["nBins"] for r in L.get("regions", [])), [1] * 11 + [2], "region table: 11 ROIs of one AOI, one of two")

    # ---- counts: targets as collapsed from the probes ------------------------------------------
    gl = read_json(d / "genes.json")
    c.eq([g["n"] for g in gl], genes, "rois genes.json: the targets in PKC order, a name shared by two panels kept apart")
    c.ok(not any(n in [g["n"] for g in gl] for n in tr["negativeTargets"]), "negative-probe targets are not genes")
    V = np.array([a["values"] for a in aois], np.float64).T                     # genes x AOIs
    c.ok(np.allclose([g["t"] for g in gl], V.sum(1), rtol=1e-5, atol=1e-3), "genes.json totals = sums of the collapsed values")
    worst, files = 0.0, 0
    for k, name in enumerate(genes):
        a = gene_f32(d, name)
        if V[k].sum() == 0:
            c.ok(a is None, f"target {name} was never counted: no gene file")
            continue
        files += 1
        if not c.ok(a is not None, f"gene file of {name}"):
            continue
        got = np.zeros(len(aois))
        for x, y, v in a:
            j = int(np.argmin(np.abs(pos[:, 0] - x) + np.abs(pos[:, 1] - y)))
            got[j] = v
        worst = max(worst, float(np.abs(got - V[k]).max() / max(1.0, V[k].max())))
    c.ok(worst < 1e-5, f"gene file values = geometric mean of the probes with the + 0.5 rule (worst relative difference {worst:.2e})")
    k2 = genes.index("Multi2")
    c.ok(np.allclose(V[k2], 0.5) and gene_f32(d, "Multi2")[:, 2].tolist() == [0.5] * len(aois),
         "a multi-probe target with every probe at 0 reads 0.5 in every AOI (as GeomxTools)")
    c.ok(gene_bytes(d, "TRB1/2") is not None and gene_bytes(d, "H2-T23") is not None, "gene files of 'TRB1/2' and 'H2-T23'")

    # ---- QC values per AOI -----------------------------------------------------------------------
    q = L.get("binQc") or {}
    met = {x["key"]: x for x in q.get("metrics", [])}
    for key in ("total", "genes", "rawReads", "pctAligned", "saturation", "area", "nuclei", "negGeoMean", "loq", "genesAboveLoq", "q3", "q3OverNeg"):
        c.ok(key in met and (d / "binqc" / met[key]["file"]).exists(), f"binqc metric {key}")
    qv = lambda key: np.fromfile(d / "binqc" / met[key]["file"], dtype="<f4") if key in met else np.zeros(0)
    c.ok(np.allclose(qv("total"), V.sum(0), rtol=1e-5), "binqc total = sum of the target values (negative probes not in it)")
    c.ok(np.allclose(qv("rawReads"), [a["raw"] for a in aois]), "binqc rawReads = Raw of the DCC")
    c.ok(np.allclose(qv("saturation"), [100 * (1 - a["dedup"] / a["aligned"]) for a in aois], atol=1e-3),
         "binqc saturation = 100 x (1 - deduplicated / aligned)")
    c.ok(np.allclose(qv("negGeoMean"), [a["negGeoMean"] for a in aois], rtol=1e-5), "binqc negGeoMean = geometric mean of the main panel's negative probes")
    c.ok(np.allclose(qv("loq"), [a["loq"] for a in aois], rtol=1e-5), "binqc loq = max(2, NegGeoMean x NegGeoSD^2)")
    c.ok(np.allclose(qv("area"), [a["area"] for a in aois], rtol=1e-5) and np.allclose(qv("nuclei"), [a["nuclei"] for a in aois]),
         "binqc area and nuclei from the annotation")
    c.ok(np.allclose(qv("q3"), np.percentile(V, 75, axis=0), rtol=1e-5), "binqc q3 = third quartile of the AOI's target values")
    q3e, q3t = L.get("q3") or {}, (tr.get("q3") or {}).get("scan") or {}
    c.ok(len(q3e.get("normFactors", [])) == len(aois) == len(q3t.get("normFactors", []))
         and np.allclose(q3e["normFactors"], q3t["normFactors"], rtol=1e-5),
         "rois q3 normFactors = GeomxTools Q3 (75th pct over targets above LOQ in >= 10% of AOIs + negative probes) / ngeoMean, as the generator's plain loops")
    c.eq((q3e.get("genesKept"), q3e.get("negativeTargets"), q3e.get("default")), (len(q3t.get("targets", [])), 2, True),
         "rois q3 entry: targets kept by the LOQ filter, two negative-probe targets, Q3 the default")

    # ---- groupings from the annotation, no embedding, nothing computed for cells / bins ----------
    u = L.get("umap") or {}
    ix = read_json(d / "umap" / "index.json")
    c.eq((u.get("embedding"), u.get("nEmbedded"), (ix.get("embedding") or {}).get("source")), ("none", 0, "none"), "rois layer: groupings without an embedding")
    c.ok(bool(np.isnan(f32(d / "umap" / "coords.bin", 2)).all()), "umap/coords.bin is all NaN (no UMAP of ROIs)")
    c.ok(not list((d / "umap").glob("clusters_r*.bin")) and not (d / "umap" / "celltypes.json").exists() and "celltypes" not in u,
         "no Leiden clusters and no suggested cell types for an ROI layer")
    names = [x["name"] for x in ix.get("clusterings", [])]
    c.eq(names, ["tissue", "genotype", "Region", "segment", "Slide"],
         "groupings: the study's columns, segment, slide (scan = the same groups, panel = one value, the display name = an id: left out)")
    c.ok(all(x.get("source") == "user" for x in ix.get("clusterings", [])), "groupings carry source 'user'")
    for col in ("tissue", "genotype", "segment"):
        cl = next((x for x in ix.get("clusterings", []) if x["name"] == col), None)
        if not c.ok(cl is not None, f"grouping {col}"):
            continue
        lab = np.fromfile(d / "umap" / cl["file"], dtype="<u2")
        got = [cl["categories"][v] for v in lab]
        c.eq(got, [a[col] for a in aois], f"grouping {col}: each AOI in its annotation group")
        c.eq(cl.get("sizes"), [sum(a[col] == k for a in aois) for k in cl["categories"]], f"grouping {col}: group sizes")
    sl = next((x for x in ix.get("clusterings", []) if x["name"] == "Slide"), {})
    c.eq((sl.get("categories"), sl.get("sizes")), (["Slide 1", "Slide 2"], [7, 6]), "slide names replaced by Slide 1, Slide 2")
    ti = next((x for x in ix.get("clusterings", []) if x["name"] == "tissue"), {})
    tops = {cat: (mk[0]["gene"] if mk else None) for cat, mk in zip(ti.get("categories", []), ti.get("markers", []))}
    c.eq(tops, {"cortex": "Gene1", "stroma": "Gene3", "tumour": "Gene2"}, "top marker of each tissue group = the gene the demo raised there")
    hv = read_json(d / "hvg" / "index.json")
    hg = read_json(d / "hvg" / "genes.json")
    c.eq((L.get("hvg") or {}).get("allTargets"), True, "rois hvg entry says all targets")
    c.eq(hg, [g for k, g in enumerate(genes) if V[k].sum() > 0], "hvg/genes.json = every target with a count (no variable-gene selection)")
    c.eq(sum(ch["n"] for ch in hv.get("chunks", [])), len(aois), "hvg chunks hold every AOI (no minimum number of rows)")
    c.ok(all((d / "hvg" / f"{ch['tx']}_{ch['ty']}.bin").stat().st_size == ch["bytes"] for ch in hv.get("chunks", [])), "hvg chunk files and sizes")

    # ---- nothing of the slide names, scan names or Sample_IDs leaves the cache --------------------
    texts = [p for p in out.rglob("*") if p.is_file() and p.suffix in (".json", ".geojson")]
    blob = "\n".join(p.read_text(encoding="utf-8", errors="ignore") for p in texts)
    for bad, what in ((tr["fakeId"][-4:], "the slide names' sample ID"), ("block", "a slide name"), ("demoA", "a scan name"),
                      ("DSP-", "a Sample_ID"), ("1001000000001", "the plate number")):
        c.ok(bad not in blob, f"geomx export text contains {what} ({bad!r})")

    # ---- the same study without images: a drawn grid ---------------------------------------------
    og = root / "export" / GXG
    mg = read_json(og / "manifest.json")
    check_manifest_common(c, mg, GXG, "png")
    Lg = (mg.get("expr") or [{}])[0]
    c.eq((mg.get("kind"), mg.get("pixelSizeUm"), (mg.get("geomx") or {}).get("layout"), Lg.get("nBins")),
         ("geomx", None, "grid", len(tr["aois"])), "grid dataset: kind, no pixel size, layout, every AOI of the annotation")
    c.eq([ch["name"] for ch in mg["image"]["channels"]], ["ROI grid (no scan image)"], "grid dataset: one drawn channel that says what it is")
    check_tiles(c, og, mg, "png")
    pg = f32(og / "expr" / "rois" / "positions.bin", 2)
    c.ok(len(pg) == len(tr["aois"]) and len({(float(x), float(y)) for x, y in pg}) == len(pg)
         and bool((pg[:, 0] > 0).all() and (pg[:, 0] < mg["image"]["width"]).all() and (pg[:, 1] > 0).all() and (pg[:, 1] < mg["image"]["height"]).all()),
         "grid dataset: one cell per AOI, all on the canvas")
    c.eq(len(read_json(og / "overlays" / "rois.geojson")["features"]), len(tr["aois"]), "grid dataset: one circle per AOI")
    ig = read_json(og / "expr" / "rois" / "umap" / "index.json")
    c.eq([(x["name"], x["sizes"]) for x in ig.get("clusterings", [])],
         [("tissue", [5, 4, 6]), ("genotype", [7, 8]), ("segment", [13, 1, 1]), ("Slide", [7, 6, 2])],
         "grid dataset (lab worksheet + CSV): groupings and sizes")
    tg = read_json(og / "expr" / "rois" / "genes.json")
    Vg = np.array([a["values"] for a in tr["aois"]], np.float64).T
    c.ok([g["n"] for g in tg] == genes and np.allclose([g["t"] for g in tg], Vg.sum(1), rtol=1e-5, atol=1e-3),
         "grid dataset: the same targets and totals over all 15 AOIs")
    fg, ft = (Lg.get("q3") or {}).get("normFactors", []), ((tr.get("q3") or {}).get("grid") or {}).get("normFactors", [])
    c.ok(len(fg) == len(ft) == len(tr["aois"]) and np.allclose(fg, ft, rtol=1e-5), "grid dataset: Q3 normFactors over its 15 AOIs, as the generator's plain loops")
    blob = "\n".join(p.read_text(encoding="utf-8", errors="ignore") for p in og.rglob("*") if p.is_file() and p.suffix in (".json", ".geojson"))
    c.ok(tr["fakeId"][-4:] not in blob and "block" not in blob and "DSP-" not in blob, "grid dataset: no slide name or Sample_ID in the export text")
    return {"aois": len(aois), "targets": len(genes), "geneFiles": files, "groupings": names, "gridAois": int(len(pg))}


# ------------------------------------------------------------------ CosMx / MERSCOPE (examples/make_imaging_demo.py)
# synthetic CosMx folders: older names + CellComposite + CellLabels / newer names, .csv.gz and raw
# morphology TIFFs / tx parquet only; synthetic MERSCOPE folders: mosaics + old tables / no images,
# cell tables one folder down, re-saved gzipped transcripts / mosaics + transcripts.parquet only
CX, CXN, CXT = "demo_cosmx", "demo_cosmx_new", "demo_cosmx_tx"
MS, MSN, MSP = "demo_merscope", "demo_merscope_noimg", "demo_merscope_pq"
IMG_DS = (CX, CXN, CXT, MS, MSN, MSP)


def imaging_commands(root: Path) -> list[list[str]]:
    """The CosMx and MERSCOPE fixtures and their prepare / export commands (run after demo_commands)."""
    try:
        r = root.resolve().relative_to(REPO.resolve()).as_posix()
    except ValueError:
        r = root.resolve().as_posix()
    raw, cache, export = f"{r}/raw", f"{r}/cache", f"{r}/export"
    prep, mk = "web/SpatialViewer/prepare.py", "examples/demo_xenium_markers.json"
    return [
        ["examples/make_imaging_demo.py", "--out", raw],
        [prep, "cosmx", "--id", CX, "--src", f"{raw}/cosmx", "--out", cache, "--pixel-size", "0.5", "--label", "Synthetic CosMx"],
        [prep, "export", "--id", CX, "--out", cache, "--dest", export, "--markers", mk],
        [prep, "cosmx", "--id", CXN, "--src", f"{raw}/cosmx_new", "--out", cache, "--channels", "Membrane,DAPI"],
        [prep, "export", "--id", CXN, "--out", cache, "--dest", export, "--no-umap"],
        [prep, "cosmx", "--id", CXT, "--src", f"{raw}/cosmx_tx", "--out", cache, "--pixel-size", "0.5"],
        [prep, "export", "--id", CXT, "--out", cache, "--dest", export],
        [prep, "merscope", "--id", MS, "--src", f"{raw}/merscope", "--out", cache, "--label", "Synthetic MERSCOPE"],
        [prep, "export", "--id", MS, "--out", cache, "--dest", export, "--markers", mk],
        [prep, "merscope", "--id", MSN, "--src", f"{raw}/merscope_noimg", "--out", cache],
        [prep, "export", "--id", MSN, "--out", cache, "--dest", export, "--no-umap"],
        [prep, "merscope", "--id", MSP, "--src", f"{raw}/merscope_pq", "--out", cache],
        [prep, "export", "--id", MSP, "--out", cache, "--dest", export, "--grid-um", "0"],
    ]


def clean_imaging(root: Path):
    for p in [root / sub / ds for sub in ("cache", "export") for ds in IMG_DS]:
        if p.exists():
            print(f"  removing earlier output {p}", flush=True)
            shutil.rmtree(p)


def _img_cells(c: Checker, root: Path, ds: str, L: dict, want_px, counts, genes, controls, barcodes, px: float,
               platform: str, clustering: bool, ctype=None) -> tuple:
    """A cells layer of an imaging platform against the fixture's truth: positions (image pixels,
    matrix order), genes and totals, gene files, controls left out, cell ids kept in the cache."""
    out = root / "export" / ds
    d = out / "expr" / "cells"
    c.eq((L.get("kind"), L.get("platform")), ("cells", platform), f"{ds}: cells layer kind / platform")
    c.eq(L.get("nBins"), len(want_px), f"{ds}: cells layer nBins = cells of the count matrix")
    c.eq(L.get("nGenes"), len(genes), f"{ds}: cells layer nGenes (controls left out)")
    c.eq(L.get("binSizeUm"), None, f"{ds}: cells layer binSizeUm")
    c.ok(abs((L.get("micronsPerPixel") or 0) - px) < 1e-9, f"{ds}: cells layer micronsPerPixel {L.get('micronsPerPixel')}")
    c.ok(14 < (L.get("cellDiameterUm") or 0) < 18, f"{ds}: cells layer cellDiameterUm {L.get('cellDiameterUm')} (cells are 14-18 um across)")
    pos = f32(d / "positions.bin", 2)
    c.ok(pos.shape == want_px.shape and np.allclose(pos, want_px, atol=2e-2),
         f"{ds}: cells positions.bin = centroids in image pixels, in matrix order (the cell table is in another order)")
    gj = read_json(d / "genes.json")
    c.eq([g["n"] for g in gj], genes, f"{ds}: cells genes.json: the genes in matrix order, controls left out")
    cnt = counts.astype(np.float64)
    c.ok(np.allclose([g["t"] for g in gj], cnt.sum(0)), f"{ds}: cells genes.json totals = count-matrix column sums")
    dense, bad = {}, []
    for j, g in enumerate(genes):
        rows = np.flatnonzero(cnt[:, j] > 0)
        a = gene_f32(d, g)
        a = np.zeros((0, 3), np.float32) if a is None else a
        if len(a) != len(rows) or not np.allclose(a[:, :2], pos[rows]) or not np.allclose(a[:, 2], cnt[rows, j]):
            bad.append(g)
        dense[g] = (rows, cnt[rows, j])
    c.ok(not bad, f"{ds}: cells gene files differ from the count matrix: {bad[:4]}")
    c.ok(not any(gene_bytes(d, x) is not None for x in controls), f"{ds}: a control code has a gene file")
    try:
        import zarr
        zc = zarr.open_group(str(root / "cache" / ds / "expr" / "cells.zarr"), mode="r")
        got = [b.decode() for b in zc["barcodes"][()]] if "barcodes" in zc else None
        c.eq(got, list(barcodes), f"{ds}: cells cache barcodes = cell ids in matrix (= positions) order")
    except ImportError:
        print("    SKIP cells barcodes check (zarr not installed)")
    lib = cnt.sum(1)
    res = {"cells": int(len(pos))}
    if ctype is not None:
        check_hvg(c, d, L, pos, lib, dense)
        check_binqc(c, d, L, lib, dense, genes)
        if clustering:
            res["umap"] = check_umap(c, d, L, pos, lib, np.array(["typea", "typeb", "typec"])[ctype].tolist(), 20 / px)
        else:
            print("    SKIP umap checks for the cells layer (clustering packages not installed)")
    return pos, dense, res


def _img_tile(out: Path, m: dict, ch: int):
    """Channel `ch` of an export at full resolution (the fixtures are a few tiles), or None without Pillow."""
    try:
        from PIL import Image
    except ImportError:
        return None
    T, z, ext = m["image"]["tileSize"], m["image"]["levels"] - 1, m.get("tileFormat", "png")
    W, H = m["image"]["width"], m["image"]["height"]
    a = np.zeros((H, W), np.uint8)
    for tx in range(math.ceil(W / T)):
        for ty in range(math.ceil(H / T)):
            t = np.asarray(Image.open(out / "tiles" / str(ch) / str(z) / str(tx) / f"{ty}.{ext}"))
            a[ty * T:(ty + 1) * T, tx * T:(tx + 1) * T] = t[:min(T, H - ty * T), :min(T, W - tx * T)]
    return a


def _img_on_nuclei(c: Checker, ds: str, img, pos, what: str):
    if img is None:
        print("    SKIP stain alignment check (Pillow not installed)")
        return
    vals = img[np.clip(pos[:, 1].astype(int), 0, img.shape[0] - 1), np.clip(pos[:, 0].astype(int), 0, img.shape[1] - 1)].astype(float)
    bg = float(np.median(img))
    c.ok(np.median(vals) > max(3 * bg, 60) and (vals > 2 * bg).mean() > 0.97,
         f"{ds}: cells on {what} nuclei: median {np.median(vals):.0f} at the centroids vs {bg:.0f} overall")
    print(f"    {ds}: {what} at the {len(vals)} centroids: median {np.median(vals):.0f} (8-bit), image median {bg:.0f}")


def _img_molecules(c: Checker, root: Path, ds: str, m: dict, genes, per_gene: dict, n_ctrl: int, read: int, platform: str,
                   mol_px=None, mol_gene=None, pos=None, dense=None, px: float = 0.5) -> dict:
    """molecules/ of an imaging platform: counts per gene, controls counted out, quantised
    positions, random order, and (with cells) molecules on cells that have the gene."""
    out = root / "export" / ds
    idx, packed = _mol_files(c, out, ds)
    if not idx:
        return {}
    W, H = m["image"]["width"], m["image"]["height"]
    c.eq(idx.get("format"), "xy-uint16", f"{ds}: molecules format")
    c.eq({g["n"]: g["c"] for g in idx["genes"]}, per_gene, f"{ds}: molecules per gene (controls left out)")
    c.eq(idx.get("excluded", {}).get("controls"), n_ctrl, f"{ds}: molecules: control codes counted out")
    c.eq(idx.get("excluded", {}).get("outsideImage"), 0, f"{ds}: molecules: none outside the image")
    c.eq((idx.get("read"), idx.get("platform")), (read, platform), f"{ds}: molecules: rows read / platform")
    c.ok("qvMin" not in idx and "lowQv" not in idx, f"{ds}: molecules index carries a Xenium quality cut")
    mb = m.get("molecules") or {}
    c.eq((mb.get("n"), mb.get("nGenes"), mb.get("platform")), (sum(per_gene.values()), len(genes), platform), f"{ds}: manifest molecules block")
    sx, sy = idx["scale"]
    c.ok(abs(sx - W / 65535) < 1e-9 and abs(sy - H / 65535) < 1e-9, f"{ds}: molecules scale = image size / 65535")
    if mol_px is None:
        return {"n": int(sum(len(v) for v in packed.values()))}
    bad, near = [], []
    for j, g in enumerate(genes):
        sel = mol_gene == j
        got = packed.get(g, np.zeros(0, np.uint32))
        gxy = np.column_stack([(got & 0xFFFF).astype(np.float64) * sx, (got >> 16).astype(np.float64) * sy])
        if len(got) != int(sel.sum()):
            bad.append(g)
            continue
        # every molecule of the truth has a written one within 1.5 steps of the 16-bit grid (the
        # tables hold 4-6 decimals, so a value can round to the neighbouring step) and vice versa
        want = mol_px[sel]
        dd = ((want[:, None, :] - gxy[None, :, :]) ** 2).sum(-1)
        tol = (1.5 * max(sx, sy)) ** 2
        if dd.min(1).max(initial=0) > tol or dd.min(0).max(initial=0) > tol:
            bad.append(g)
            continue
        fx = gxy[:, 0]
        if len(fx) > 200 and abs(fx[: len(fx) // 2].mean() - fx.mean()) > 0.05 * W:
            bad.append(f"{g} (not shuffled)")
        if pos is not None:
            dd = ((gxy[:, None, :] - pos[None, :, :].astype(np.float64)) ** 2).sum(-1)
            nn = dd.argmin(1)
            has = set(dense[g][0].tolist())
            ok = (np.sqrt(dd[np.arange(len(nn)), nn]) <= 5.5 / px) & np.array([r in has for r in nn])
            near.append(float(ok.mean()))
    c.ok(not bad, f"{ds}: molecule files differ from the truth: {bad[:4]}")
    res = {"n": int(sum(len(v) for v in packed.values()))}
    if pos is not None:
        c.ok(bool(near) and min(near) > 0.8, f"{ds}: molecules on cells that have the gene: worst gene {min(near or [0]):.0%}")
        res["onExpressingCell"] = round(min(near or [0]), 3)
        print(f"    {ds}: {res['n']:,} molecules in {len(packed)} files; per gene {min(near or [0]):.0%}-{max(near or [0]):.0%} "
              f"lie within 5.5 um of a cell whose gene file has that gene")
    return res


def _img_same_molecules(c: Checker, root: Path, ds: str, ref: str, genes, what: str):
    """Two imports of the same molecules (another container / column generation) give the same files."""
    ia, pa = _mol_files(c, root / "export" / ds, ds)
    ib, pb = _mol_files(c, root / "export" / ref, ref)
    if not (ia and ib):
        return
    c.eq({g["n"]: g["c"] for g in ia["genes"]}, {g["n"]: g["c"] for g in ib["genes"]}, f"{ds}: molecules per gene = {ref} ({what})")
    c.eq(ia["bounds"], ib["bounds"], f"{ds}: molecules bounds = {ref}")
    e = np.zeros(0, np.uint32)
    far = [g for g in genes if len(pa.get(g, e)) != len(pb.get(g, e)) or
           np.abs(np.sort(pa.get(g, e) & 0xFFFF).astype(np.int64) - np.sort(pb.get(g, e) & 0xFFFF).astype(np.int64)).max(initial=0) > 1]
    c.ok(not far, f"{ds}: molecule x positions differ from {ref}: {far[:4]}")


def _img_background(c: Checker, root: Path, ds: str, m: dict, mol_px):
    """No images: a one-channel transcript-density background, labelled as such, bright where the
    molecules are and empty in the margin; no protein grid is built on it."""
    img = m.get("image", {})
    c.eq(img.get("background"), "transcript density", f"{ds}: manifest image.background")
    c.eq([ch["name"] for ch in img.get("channels", [])], ["Transcript density (no image)"], f"{ds}: background channel name")
    c.ok(abs((m.get("pixelSizeUm") or 0) - 1.0) < 1e-9, f"{ds}: background pixel size {m.get('pixelSizeUm')} (1 um)")
    c.ok(not [L for L in m.get("expr", []) if L.get("kind") == "protein"], f"{ds}: a protein grid was built on the density background")
    a = _img_tile(root / "export" / ds, m, 0)
    if a is None:
        return
    at = a[np.clip(mol_px[:, 1].astype(int), 0, a.shape[0] - 1), np.clip(mol_px[:, 0].astype(int), 0, a.shape[1] - 1)]
    c.ok(float(np.median(at)) > 40 and int(a[:5].max()) == 0 and int(a[:, :5].max()) == 0,
         f"{ds}: density background: median {np.median(at):.0f} at the molecules, margin max {max(int(a[:5].max()), int(a[:, :5].max()))}")


def check_imaging(c: Checker, root: Path, clustering: bool, log_path: Path) -> dict:
    """The synthetic CosMx and MERSCOPE folders (examples/make_imaging_demo.py).

    CosMx: four FOV images stitched by their positions with y_local turned over (global y up),
    2 / 3 px gaps kept black, cell_ID 0 rows and NegPrb / SystemControl codes left out, ids
    <cell_ID>_<fov>; the newer names (csv.gz, cell column, Negative*, mm positions, raw
    morphology TIFFs) give the same cells and molecules; transcripts alone (parquet) give a
    density background and the same molecules. MERSCOPE: the middle z-plane, positions through
    the micron -> mosaic transform, Blank-* left out, 39-digit cell ids kept as text; no images
    + tables one folder down + a re-saved transcript table give the same cells and molecules on
    a density background; mosaics + transcripts.parquet alone give the molecules and no cells."""
    raw = root / "raw"
    tr = read_json(raw / "truth_imaging.json")
    tz = np.load(raw / "truth_imaging.npz")
    px, genes, n = tr["pixelSizeUm"], tr["genes"], tr["nCells"]
    G = len(genes)
    mol_gene = tz["mol_gene"]
    res = {}

    # ---------------- CosMx, older names: CellComposite mosaic, cells, molecules
    t = tr["cosmx"]
    W, H = t["width"], t["height"]
    out = root / "export" / CX
    m = read_json(out / "manifest.json")
    check_manifest_common(c, m, CX, "png")
    c.eq(m.get("kind"), "cosmx", "cosmx manifest kind")
    c.eq((m["image"]["width"], m["image"]["height"]), (W, H), "cosmx mosaic size (2 x 2 FOVs and the gaps between them)")
    c.eq([ch["name"] for ch in m["image"]["channels"]], ["CellComposite red", "CellComposite green", "CellComposite blue"],
         "cosmx channel names (CellComposite)")
    c.ok(abs((m.get("pixelSizeUm") or 0) - px) < 1e-9, f"cosmx pixelSizeUm {m.get('pixelSizeUm')}")
    aff = np.array((m.get("frame") or {}).get("affine", [[0, 0, 0], [0, 0, 0]]), float)
    c.ok(aff.shape == (2, 3) and aff[0, 0] == 1 and aff[1, 1] == -1 and aff[0, 1] == 0 and aff[1, 0] == 0,
         f"cosmx manifest frame: x kept, global y turned over (got {aff.tolist()})")
    g_xy = tz["cosmx_mol_global"].astype(np.float64)
    back = np.column_stack([g_xy[:, 0] * aff[0, 0] + aff[0, 2], g_xy[:, 1] * aff[1, 1] + aff[1, 2]])
    c.ok(np.allclose(back, tz["cosmx_mol_px"], atol=1e-3), "cosmx manifest frame maps x/y_global_px to the mosaic pixels")
    check_tiles(c, out, m, "png")
    layers = {L["key"]: L for L in m.get("expr", [])}
    c.eq(sorted(layers), ["cells"], "cosmx expr layers (a cells layer, no protein grid beside it)")
    # the count matrix lists the cells FOV by FOV (after each FOV's cell_ID 0 row): that is the layer's order
    order = np.concatenate([np.flatnonzero(tz["cosmx_fov"] == f) for f in range(1, 5)])
    c_fov, c_id = tz["cosmx_fov"][order], tz["cosmx_cell_id"][order]
    c_px, c_cnt = tz["cosmx_cells_px"][order].astype(np.float32), tz["cosmx_counts"][order]
    bc = [f"{i}_{f}" for i, f in zip(c_id.tolist(), c_fov.tolist())]
    c.eq(len(set(bc)), n, "cosmx fixture: <cell_ID>_<fov> ids are unique (cell_ID restarts in every FOV)")
    pos, dense, res[CX] = _img_cells(c, root, CX, layers.get("cells", {}), c_px, c_cnt,
                                     genes, t["controls"], bc, px, "cosmx", clustering, ctype=tz["ctype"][order])
    blue = _img_tile(out, m, 2)
    _img_on_nuclei(c, CX, blue, pos, "CellComposite blue")
    if blue is not None:
        fw, fh = t["fovSize"]
        c.ok(int(blue[fh:H - fh].max()) == 0 and int(blue[:, fw:W - fw].max()) == 0,
             "cosmx mosaic: the gap between the FOV images is not black (a FOV is misplaced)")
        for f in range(1, 5):                              # every FOV carries cells that sit on its own nuclei
            sel = c_fov == f
            v = blue[pos[sel, 1].astype(int), pos[sel, 0].astype(int)]
            c.ok(float(np.median(v)) > 100, f"cosmx FOV {f}: nuclei under its {int(sel.sum())} centroids (median {np.median(v):.0f})")
    keep = mol_gene < G
    res[CX]["molecules"] = _img_molecules(c, root, CX, m, genes, t["moleculesPerGene"], t["controlMolecules"], t["read"], "cosmx",
                                          tz["cosmx_mol_px"][keep], mol_gene[keep], pos, dense, px)

    # ---------------- CosMx, newer names: csv.gz, cell column, Negative*, mm positions, raw morphology
    mn = read_json(root / "export" / CXN / "manifest.json")
    c.eq((mn["image"]["width"], mn["image"]["height"], mn["image"]["dtype"]), (W, H, "uint16"), f"{CXN}: mosaic of the raw morphology TIFFs")
    c.eq([ch["name"] for ch in mn["image"]["channels"]], ["Membrane", "DAPI"], f"{CXN}: channel names from --channels")
    c.ok(abs((mn.get("pixelSizeUm") or 0) - px) < 1e-6, f"{CXN}: pixel size from the mm and px columns of the FOV positions ({mn.get('pixelSizeUm')})")
    ln = {L["key"]: L for L in mn.get("expr", [])}
    c.eq(sorted(ln), ["cells"], f"{CXN}: expr layers")
    bcn = [f"c_1_{f}_{i}" for i, f in zip(c_id.tolist(), c_fov.tolist())]
    posn, _d, _r = _img_cells(c, root, CXN, ln.get("cells", {}), c_px, c_cnt, genes,
                              t["controlsNew"], bcn, px, "cosmx", clustering)
    _img_on_nuclei(c, CXN, _img_tile(root / "export" / CXN, mn, 1), posn, "DAPI (raw morphology)")
    _img_molecules(c, root, CXN, mn, genes, t["moleculesPerGene"], t["controlMolecules"], t["read"], "cosmx")
    _img_same_molecules(c, root, CXN, CX, genes, "csv.gz with the newer column set")

    # ---------------- CosMx, transcripts only (parquet): density background, the same molecules
    mt = read_json(root / "export" / CXT / "manifest.json")
    c.eq(mt.get("expr"), [], f"{CXT}: no expression layer without cell files (and no grid on the background)")
    lo, hi = g_xy.min(0), g_xy.max(0)
    bg_px = np.column_stack([(g_xy[:, 0] - lo[0]) * px + 10.0, (hi[1] - g_xy[:, 1]) * px + 10.0])    # 1 um pixels, 10 um margin
    c.eq((mt["image"]["width"], mt["image"]["height"]),
         (math.ceil((hi[0] - lo[0]) * px + 20), math.ceil((hi[1] - lo[1]) * px + 20)), f"{CXT}: background size = extent + margins at 1 um")
    _img_background(c, root, CXT, mt, bg_px[keep])
    _img_molecules(c, root, CXT, mt, genes, t["moleculesPerGene"], t["controlMolecules"], t["read"], "cosmx", bg_px[keep], mol_gene[keep], px=1.0)

    # ---------------- MERSCOPE: mosaics (middle z), cells, molecules
    t = tr["merscope"]
    W, H = t["width"], t["height"]
    out = root / "export" / MS
    m = read_json(out / "manifest.json")
    check_manifest_common(c, m, MS, "png")
    c.eq(m.get("kind"), "merscope", "merscope manifest kind")
    c.eq((m["image"]["width"], m["image"]["height"]), (W, H), "merscope image size")
    c.eq([ch["name"] for ch in m["image"]["channels"]], ["DAPI", "PolyT"], "merscope channel names (DAPI first)")
    c.ok(abs((m.get("pixelSizeUm") or 0) - px) < 1e-9, f"merscope pixelSizeUm from the transform ({m.get('pixelSizeUm')})")
    c.eq((m.get("frame") or {}).get("zPlane"), 1, "merscope: the middle z-plane of 0, 1, 2")
    check_tiles(c, out, m, "png")
    layers = {L["key"]: L for L in m.get("expr", [])}
    c.eq(sorted(layers), ["cells"], "merscope expr layers (a cells layer, no protein grid beside it)")
    org = np.array(t["origin"], np.float64)
    cells_px = ((tz["merscope_cells_um"] - org) / px).astype(np.float32)
    pos, dense, res[MS] = _img_cells(c, root, MS, layers.get("cells", {}), cells_px, tz["merscope_counts"], genes, t["controls"],
                                     t["idsOld"], px, "merscope", clustering, ctype=tz["ctype"])
    _img_on_nuclei(c, MS, _img_tile(out, m, 0), pos, "DAPI")
    try:
        import zarr
        z0 = zarr.open_group(str(root / "cache" / MS / "image.zarr"), mode="r")["0"]
        v = np.asarray(z0[0])[pos[:, 1].astype(int), pos[:, 0].astype(int)]
        c.ok(float(np.median(v)) > 4000, f"merscope: the in-focus plane was imported (DAPI {np.median(v):.0f} at the centroids; the other planes give ~1800)")
    except ImportError:
        print("    SKIP z-plane brightness check (zarr not installed)")
    mg = tz["merscope_mol_gene"]
    mpx = (tz["merscope_mol_um"] - org) / px
    keep = mg < G
    res[MS]["molecules"] = _img_molecules(c, root, MS, m, genes, t["moleculesPerGene"], t["controlMolecules"], t["read"], "merscope",
                                          mpx[keep], mg[keep], pos, dense, px)

    # ---------------- MERSCOPE, no images: tables one folder down, EntityID, re-saved gz transcripts
    mn = read_json(root / "export" / MSN / "manifest.json")
    ln = {L["key"]: L for L in mn.get("expr", [])}
    c.eq(sorted(ln), ["cells"], f"{MSN}: expr layers (cell tables found one folder down; no grid on the background)")
    um_all = np.vstack([tz["merscope_mol_um"], tz["merscope_cells_um"]])
    lo, hi = um_all.min(0), um_all.max(0)
    c.eq((mn["image"]["width"], mn["image"]["height"]), (math.ceil(hi[0] - lo[0] + 20), math.ceil(hi[1] - lo[1] + 20)),
         f"{MSN}: background size = extent + margins at 1 um")
    bgc = (tz["merscope_cells_um"] - lo + 10.0).astype(np.float32)
    _p, _d, _r = _img_cells(c, root, MSN, ln.get("cells", {}), bgc, tz["merscope_counts"], genes, t["controls"], t["idsNew"], 1.0,
                            "merscope", clustering)
    bgm = tz["merscope_mol_um"] - lo + 10.0
    _img_background(c, root, MSN, mn, bgm[keep])
    _img_molecules(c, root, MSN, mn, genes, t["moleculesPerGene"], t["controlMolecules"], t["read"], "merscope", bgm[keep], mg[keep], px=1.0)

    # ---------------- MERSCOPE, mosaics + transcripts.parquet only
    mp = read_json(root / "export" / MSP / "manifest.json")
    c.eq(mp.get("expr"), [], f"{MSP}: no expression layer without cell tables (and --grid-um 0)")
    _img_molecules(c, root, MSP, mp, genes, t["moleculesPerGene"], t["controlMolecules"], t["read"], "merscope")
    _img_same_molecules(c, root, MSP, MS, genes, "parquet")

    # ---------------- publishable: no local path, no cell id, no leftover sample id
    export = root / "export"
    texts = [p for ds in IMG_DS for p in (export / ds).rglob("*") if p.is_file() and p.suffix in (".json", ".geojson", ".txt")]
    lp = [p for p in texts if LOCAL_PATH.search(p.read_text(encoding="utf-8", errors="ignore"))]
    c.ok(not lp, f"imaging exports: text files contain local paths: {lp[:3]}")
    ids = [p for p in texts if any(s in p.read_text(encoding="utf-8", errors="ignore") for s in (t["idsOld"][0], t["idsNew"][0], "c_1_1_1"))]
    c.ok(not ids, f"imaging exports: text files contain cell ids: {ids[:3]}")
    args = ["web/SpatialViewer/deid.py", "check"] + [str(export / ds) for ds in IMG_DS]
    p = subprocess.run([sys.executable, *args], cwd=str(REPO), stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                       text=True, encoding="utf-8", errors="replace")
    with open(log_path, "a", encoding="utf-8") as lg:
        lg.write(p.stdout)
    c.eq(p.returncode, 0, "deid.py check on the CosMx / MERSCOPE exports (exit code)")
    return res


def _st_truth(tz, n: int, W: int, G: int, whole_chip: bool = False):
    """Bins of n x n DNBs from the drawn records: (occupied bin ids in row order, centres px, bins x genes counts)."""
    sel = np.ones(len(tz["x"]), bool) if whole_chip else tz["inside"]
    nbx = -(-W // n)
    b = (tz["y"][sel] // n) * nbx + tz["x"][sel] // n
    cells = np.unique(b)
    dense = np.zeros((len(cells), G))
    np.add.at(dense, (np.searchsorted(cells, b), tz["gene"][sel]), tz["count"][sel])
    pos = np.column_stack([(cells % nbx + 0.5) * n, (cells // nbx + 0.5) * n]).astype(np.float32)
    return cells, pos, dense


def _st_layer_files(c: Checker, d: Path, what: str, names: list[str], pos_want, dense_want):
    """positions.bin, genes.json and every gene file of one exported layer against a dense truth matrix."""
    pos = f32(d / "positions.bin", 2)
    c.ok(pos.shape == pos_want.shape and np.array_equal(pos, pos_want), f"{what}: positions.bin = the drawn rows' centres")
    gj = read_json(d / "genes.json")
    c.eq([g["n"] for g in gj], names, f"{what}: genes.json names")
    c.ok(len(gj) == dense_want.shape[1] and np.allclose([g["t"] for g in gj], dense_want.sum(0)),
         f"{what}: genes.json totals = the drawn counts per gene")
    dense, bad = {}, []
    for j, g in enumerate(names):
        rows = np.flatnonzero(dense_want[:, j] > 0)
        a = gene_f32(d, g)
        a = np.zeros((0, 3), np.float32) if a is None else a
        if len(a) != len(rows) or not np.array_equal(a[:, :2], pos_want[rows]) or not np.array_equal(a[:, 2], dense_want[rows, j]):
            bad.append(g)
        dense[g] = (rows, dense_want[rows, j])
    c.ok(not bad, f"{what}: gene files differ from the drawn counts: {bad[:4]}")
    c.eq(gene_count(d), int((dense_want.sum(0) > 0).sum()), f"{what}: genes with a payload (pack entries / gene files)")
    return pos, dense


def _tile_values(out: Path, m: dict, ch: int, pos) -> np.ndarray:
    """8-bit value of channel `ch` of the full-resolution tiles at each (x, y) pixel position."""
    from PIL import Image
    T, z = m["image"]["tileSize"], m["image"]["levels"] - 1
    tiles, vals = {}, []
    for x, y in pos:
        k = (int(x // T), int(y // T))
        if k not in tiles:
            tiles[k] = np.asarray(Image.open(out / "tiles" / str(ch) / str(z) / str(k[0]) / f"{k[1]}.png"))
        vals.append(int(tiles[k][int(y) - k[1] * T, int(x) - k[0] * T]))
    return np.array(vals)


def check_stereoseq(c: Checker, root: Path, clustering: bool, log_path: Path) -> dict:
    """The synthetic Stereo-seq chip (examples/make_stereoseq_demo.py) through its four routes."""
    raw = root / "raw"
    tr = read_json(raw / "truth_stereoseq.json")
    tz = np.load(raw / "truth_stereoseq.npz")
    W, H, um, labels = tr["width"], tr["height"], tr["dnbUm"], tr["labels"]
    G = len(labels)
    out = root / "export" / ST
    m = read_json(out / "manifest.json")
    check_manifest_common(c, m, ST, "png")
    c.eq(m.get("kind"), "stereoseq", "stereoseq manifest kind")
    c.eq((m["image"]["width"], m["image"]["height"]), (W, H), "stereoseq image size")
    c.eq([ch["name"] for ch in m["image"]["channels"]], ["Red", "Green", "Blue"], "stereoseq channel names")
    c.ok(abs((m.get("pixelSizeUm") or 0) - um) < 1e-9, f"stereoseq pixelSizeUm {m.get('pixelSizeUm')} (resolution 500 nm, 1 px per DNB)")
    check_tiles(c, out, m, "png")
    sb = m.get("stereoseq") or {}
    c.eq((sb.get("source"), sb.get("gefVersion"), sb.get("levelsStored")), ("gef", 4, [1]), "stereoseq block: source / version / levels")
    c.eq((sb.get("recordsRead"), sb.get("midRead")), (tr["records"], tr["mid"]), "stereoseq block: whole-chip records / MID read")
    c.eq((sb.get("recordsInTissue"), sb.get("midInTissue")), (tr["recordsInTissue"], tr["midInTissue"]),
         "stereoseq block: records / MID inside the tissue contour (hole cut out, island kept, background left out)")
    c.ok(tr["midInTissue"] < tr["mid"], "the demo chip has background counts off the tissue")
    c.ok(str(sb.get("tissueCut", "")).startswith("GEF contour"), f"stereoseq block: tissueCut {sb.get('tissueCut')!r}")
    layers = {L["key"]: L for L in m.get("expr", [])}
    c.eq(sorted(layers), ["bin020", "bin050", "cells"], "stereoseq expr layers")
    try:
        import zarr
    except ImportError:
        zarr = None
        print("    SKIP barcode checks (zarr not installed)")
    res = {}
    for n in (20, 50):
        key = f"bin{n:03d}"
        L = layers.get(key, {})
        d = out / "expr" / key
        print(f"  expr {key}:")
        cells, pos_w, dense_w = _st_truth(tz, n, W, G)
        c.eq(L.get("nBins"), len(cells), f"{key}: nBins = bins holding a tissue DNB count")
        c.eq((L.get("nGenes"), L.get("binSizeUm"), L.get("binDnb")), (G, n * um, n), f"{key}: nGenes / binSizeUm / binDnb")
        c.ok(abs((L.get("spotDiameterPx") or 0) - n) < 1e-9 and abs((L.get("micronsPerPixel") or 0) - um) < 1e-9,
             f"{key}: spotDiameterPx {L.get('spotDiameterPx')} / micronsPerPixel {L.get('micronsPerPixel')}")
        c.eq(L.get("midCount"), tr["midInTissue"], f"{key}: midCount = MID inside the tissue")
        c.ok("kind" not in L, f"{key}: a bin layer has no kind")
        pos, dense = _st_layer_files(c, d, key, labels, pos_w, dense_w)
        lib = dense_w.sum(1)
        c.ok(abs(lib.sum() - tr["midInTissue"]) < 1e-6, f"{key}: sum of the layer = MID inside the tissue")
        if zarr is not None:
            zc = zarr.open_group(str(root / "cache" / ST / "expr" / f"{key}.zarr"), mode="r")
            nbx = -(-W // n)
            want_bc = [f"{(b % nbx) * n}_{(b // nbx) * n}" for b in cells.tolist()]
            c.eq([b.decode() for b in zc["barcodes"][()]], want_bc, f"{key}: cache barcodes = x_y of each bin's first DNB")
        check_hvg(c, d, L, pos, lib, dense)
        check_binqc(c, d, L, lib, dense, labels)
        res[key] = {"bins": int(len(pos)), "medianMid": float(np.median(lib))}
        truth = np.array(tr["types"])[np.searchsorted(np.array(tr["bandEdges"]), pos[:, 0], side="right")].tolist()
        if clustering:
            res[key]["umap"] = check_umap(c, d, L, pos, lib, truth, n)
        else:
            print(f"    SKIP umap checks for {key} (clustering packages not installed)")
        if n == 20:
            try:                                  # bins lie on pink tissue, not on white glass (green is low on tissue)
                gv = _tile_values(out, m, 1, pos)
                c.ok((gv < 200).mean() > 0.9, f"{key}: {(gv < 200).mean():.0%} of bin centres are on tissue in the image")
                # orientation: the small island (top right in the image) has bins on it only as imported
                ix, iy, ia = tr["island"][:3]
                isl = (np.abs(pos[:, 0] - ix) <= ia) & (np.abs(pos[:, 1] - iy) <= ia)
                on = float((gv[isl] < 200).mean())
                fy = float((_tile_values(out, m, 1, np.column_stack([pos[isl, 0], H - pos[isl, 1]])) < 200).mean())
                fx = float((_tile_values(out, m, 1, np.column_stack([W - pos[isl, 0], pos[isl, 1]])) < 200).mean())
                c.ok(isl.sum() >= 10 and on > 0.6 and fy < 0.1 and fx < 0.1,
                     f"{key}: island bins on tissue {on:.0%} (y flipped {fy:.0%}, x flipped {fx:.0%}): x right, y down, no flip")
                print(f"    image under the bin centres: {(gv < 200).mean():.0%} tissue; the {int(isl.sum())} island bins {on:.0%} "
                      f"(y flipped {fy:.0%}, x flipped {fx:.0%})")
            except ImportError:
                print("    SKIP image alignment check (Pillow not installed)")

    # cells
    print("  expr cells:")
    L = layers.get("cells", {})
    d = out / "expr" / "cells"
    cxy, ccnt, carea = tz["cell_xy"], tz["cell_counts"].astype(np.float64), tz["cell_area"]
    c.eq((L.get("kind"), L.get("nBins"), L.get("nGenes"), L.get("binSizeUm")), ("cells", tr["nCells"], G, None),
         "cells layer kind / nBins / nGenes / binSizeUm")
    c.eq(L.get("midCount"), tr["cellMid"], "cells layer midCount = MID of the cell table")
    want_d = 2 * math.sqrt(float(np.median(carea)) / math.pi)
    c.ok(abs((L.get("spotDiameterPx") or 0) - want_d) < 0.02, f"cells spotDiameterPx {L.get('spotDiameterPx')} (median cell area)")
    pos, dense = _st_layer_files(c, d, "cells", labels, (cxy + 0.5).astype(np.float32), ccnt)
    lib = ccnt.sum(1)
    if zarr is not None:
        zc = zarr.open_group(str(root / "cache" / ST / "expr" / "cells.zarr"), mode="r")
        c.eq([b.decode() for b in zc["barcodes"][()]], [str(i) for i in range(tr["nCells"])], "cells cache barcodes = cell ids")
    check_hvg(c, d, L, pos, lib, dense)
    check_binqc(c, d, L, lib, dense, labels)
    res["cells"] = {"cells": int(len(pos))}
    if clustering:
        res["cells"]["umap"] = check_umap(c, d, L, pos, lib, np.array(tr["types"])[tz["cell_type"]].tolist(), 40)
    else:
        print("    SKIP umap checks for the cells layer (clustering packages not installed)")
    try:                                          # every cell centre is a dark nucleus in the H&E
        rv = _tile_values(out, m, 0, pos)
        c.ok(np.median(rv) < 110 and (rv < 150).mean() > 0.97, f"cells on nuclei: red {np.median(rv):.0f} at the centres")
        print(f"    red channel at the {len(rv)} cell centres: median {np.median(rv):.0f} (tissue ~226, glass ~236)")
    except ImportError:
        print("    SKIP nucleus alignment check (Pillow not installed)")

    # the other three routes give the same layers
    def same(ds, key, names_a, names_b, ref_ds=ST):
        a, b = root / "export" / ref_ds / "expr" / key, root / "export" / ds / "expr" / key
        ok = (b / "positions.bin").exists() and (a / "positions.bin").read_bytes() == (b / "positions.bin").read_bytes()
        for ga, gb in zip(names_a, names_b):
            ok = ok and gene_bytes(a, ga) == gene_bytes(b, gb)
        return ok

    mg = read_json(root / "export" / STG / "manifest.json")
    c.eq([L["key"] for L in mg.get("expr", [])], ["bin020", "bin050"], f"{STG}: layers from the tissue GEM (no cells)")
    c.eq(((mg.get("stereoseq") or {}).get("source"), (mg.get("stereoseq") or {}).get("midInTissue")), ("gem", tr["midInTissue"]),
         f"{STG}: stereoseq block (GEM source, already cut to the tissue)")
    c.eq((mg.get("stereoseq") or {}).get("fileOffsetDnb"), tr["cut"], f"{STG}: the GEM's #OffsetX / #OffsetY are read as its offset")
    c.ok(abs((mg.get("pixelSizeUm") or 0) - um) < 1e-9, f"{STG}: pixelSizeUm (0.5 um per DNB assumed for a GEM)")
    for key in ("bin020", "bin050"):
        c.ok(same(STG, key, labels, labels), f"{STG}: {key} from the tissue GEM = from the raw GEF cut by its contour")
        c.eq([g["n"] for g in read_json(root / "export" / STG / "expr" / key / "genes.json")], labels, f"{STG}: {key} gene names")
    old = tr["oldLayoutNames"]
    ml = read_json(root / "export" / STL / "manifest.json")
    c.eq(((ml.get("stereoseq") or {}).get("gefVersion"), (ml.get("stereoseq") or {}).get("levelsStored")), (2, [10, 50]),
         f"{STL}: stereoseq block (older GEF, stored bin10 + bin50)")
    cut = tr["cut"]
    c.eq((ml.get("stereoseq") or {}).get("fileOffsetDnb"), cut, f"{STL}: the GEF's minX / minY (it is cut to a region) are read as its offset")
    c.eq([L["key"] for L in ml.get("expr", [])], ["bin020", "bin050", "cells"], f"{STL}: layers")
    c.ok(same(STL, "cells", labels, labels), f"{STL}: cells of the cut cell bin GEF (x stored minus offsetX) land where the whole-chip file's do")
    for key in ("bin020", "bin050"):
        c.ok(same(STL, key, labels, old), f"{STL}: {key} from stored levels ({'bin10 summed' if key == 'bin020' else 'bin50 as it is'}) = from bin1")
        c.eq([g["n"] for g in read_json(root / "export" / STL / "expr" / key / "genes.json")], old, f"{STL}: {key} gene names (one `gene` field)")
    if zarr is not None:
        cells, _, _ = _st_truth(tz, 20, W, G)
        nbx = -(-W // 20)
        zc = zarr.open_group(str(root / "cache" / STL / "expr" / "bin020.zarr"), mode="r")
        c.eq([b.decode() for b in zc["barcodes"][()]],
             [str((((b % nbx) * 20 - cut[0]) << 32) | ((b // nbx) * 20 - cut[1])) for b in cells.tolist()],
             f"{STL}: --barcodes stereopy = x * 2^32 + y, x and y as the file stores them")
        zc = zarr.open_group(str(root / "cache" / STL / "expr" / "cells.zarr"), mode="r")
        c.eq([b.decode() for b in zc["barcodes"][()]], [str(((int(x) - cut[0]) << 32) | (int(y) - cut[1])) for x, y in cxy],
             f"{STL}: cells --barcodes stereopy = x * 2^32 + y of the centre, as stored")
    mm = read_json(root / "export" / STM / "manifest.json")
    c.eq([L["key"] for L in mm.get("expr", [])], ["bin050"], f"{STM}: layers (and no protein grid of the stain)")
    c.eq([(ch["name"], ch["color"]) for ch in mm["image"]["channels"]], [("ssDNA", "#ffffff")], f"{STM}: grey ssDNA image = one white channel")
    c.eq(mm["image"]["dtype"], "uint16", f"{STM}: ssDNA image dtype")
    check_tiles(c, root / "export" / STM, mm, "png")
    bm = mm.get("stereoseq") or {}
    c.eq((bm.get("midRead"), bm.get("midInTissue")), (tr["mid"], tr["midInTissue"]), f"{STM}: whole-chip GEM cut with the mask image")
    c.ok(str(bm.get("tissueCut", "")).startswith("mask image"), f"{STM}: tissueCut {bm.get('tissueCut')!r}")
    c.ok(same(STM, "bin050", old, old, ref_ds=STL), f"{STM}: bin050 from the shuffled 4-column GEM = the other routes (by gene name)")
    c.eq(sorted(g["n"] for g in read_json(root / "export" / STM / "expr" / "bin050" / "genes.json")), sorted(old), f"{STM}: gene names")

    # the reader's parts on their own
    try:
        sys.path.insert(0, str(APP))
        import stereoseq
        mask_w = np.unpackbits(tz["mask"])[: W * H].reshape(H, W).astype(bool)
        tm = stereoseq.TissueMask.from_contour(stereoseq.read_contour(raw / "stereoseq" / "SYNCHIP.raw.gef"), W, H)
        got = np.unpackbits(tm.bits, axis=1)[:, :W].astype(bool)
        c.ok(np.array_equal(got, mask_w), f"tissue contour -> mask: {int((got != mask_w).sum())} DNBs differ from the drawn mask")
        c.eq(tm.area, tr["tissueDnb"], "tissue mask area")
        pad = np.zeros((-(-H // 50) * 50, -(-W // 50) * 50), bool)
        pad[:H, :W] = mask_w
        c.ok(np.array_equal(tm.blocks(50), pad.reshape(pad.shape[0] // 50, 50, -1, 50).any(axis=(1, 3))), "tissue mask blocks(50)")
        sh = np.zeros((1200, 1350), bool)                 # blocks counted from DNB (70, 30)
        sh[:H - 30, :W - 70] = mask_w[30:, 70:]
        c.ok(np.array_equal(tm.blocks(50, 70, 30), sh.reshape(24, 50, 27, 50).any(axis=(1, 3))), "tissue mask blocks(50) from another origin")
        xs, ys = np.array([-1, 0, W - 1, W, 700]), np.array([0, -1, H - 1, 5, 610])
        c.eq(tm.contains(xs, ys).tolist(), [False, False, bool(mask_w[H - 1, W - 1]), False, True], "tissue mask contains(): off-chip = not tissue")
        for name, whole, want_ordered in (("SYNCHIP.tissue.gem.gz", False, True), ("SYNCHIP.raw.shuffled.gem", True, False)):
            hdr = stereoseq.gem_header(raw / "stereoseq" / name)
            genes = {"ids": [], "names": [], "index": {}}
            acc = stereoseq.BinAccumulator([50], {50: 1}, W, H, root / "cache" / "_st_acc_test")
            try:
                nb = 0
                for g, x, y, v in stereoseq.gem_chunks(raw / "stereoseq" / name, hdr, genes, block_bytes=200_000):
                    acc.add(1, g, x + hdr["offset"][0], y + hdr["offset"][1], v)        # chip DNBs, as the drawn truth
                    nb += 1
                ordered = acc.levels[50].ordered
                indptr, idx, val = acc.finish(50, len(genes["ids"]))
            finally:
                acc.cleanup()
            cells_w, _, dw = _st_truth(tz, 50, W, G, whole_chip=whole)
            order = [tr["geneIds"].index(i) if i in tr["geneIds"] else old.index(i) for i in genes["ids"]]
            got = np.zeros((len(cells_w), G))
            for j, gj in enumerate(order):
                a, b = int(indptr[j]), int(indptr[j + 1])
                got[np.searchsorted(cells_w, idx[a:b]), gj] = val[a:b]
            c.ok(nb > 5 and ordered == want_ordered, f"{name} in {nb} small batches: arrives gene by gene = {ordered}, want {want_ordered}")
            c.ok(np.array_equal(got, dw), f"{name} in small batches ({'sorted at the end' if not want_ordered else 'concatenated'}) = the drawn counts")
        c.eq(stereoseq.gene_labels(["E1", "E2", "E3", "E4"], ["Abc", "ABC", "", "Abc"]), ["Abc", "ABC_E2", "E3", "Abc_E4"],
             "gene_labels: names shared by several ids (or differing only by case) stay distinct")
        c.eq(stereoseq.plan_levels([20, 50, 100], [1, 10, 50]), {20: 10, 50: 50, 100: 50}, "plan_levels: largest stored level that divides")
        # the mouse gene names H2-T23 ... are not sample IDs (the demo has one); other IDs still are
        import deid
        c.eq([deid.find_ids(s) for s in ("H2-T23", "h2-t10", "Cd4, H2-T24;")], [[], [], []], "deid: mouse MHC gene names H2-Tnn pass")
        c.eq([deid.find_ids(s) for s in ("T23", "H2-INT95", "H2_T23")], [["T23"], ["INT95"], ["T23"]],        # deid:ok
             "deid: sample IDs are still found")
        c.ok("H2-T23" in labels, "the demo has a gene named H2-T23")
        try:
            stereoseq.plan_levels([20], [50])
            c.ok(False, "plan_levels: bin20 from a file that only holds bin50 must stop")
        except SystemExit:
            c.ok(True, "")
    except ImportError as e:
        print(f"    SKIP reader unit checks ({e})")

    # nothing identifying in the exports: no chip serial number, no local path, deid clean
    texts = [p for ds in ST_DS for p in (root / "export" / ds).rglob("*") if p.is_file() and p.suffix in (".json", ".geojson", ".txt")]
    c.ok(not [p for p in texts if tr["sn"] in p.read_text(encoding="utf-8", errors="ignore")], "stereoseq exports name the chip serial number")
    c.ok(not [p for p in texts if LOCAL_PATH.search(p.read_text(encoding="utf-8", errors="ignore"))], "stereoseq exports contain local paths")
    args = ["web/SpatialViewer/deid.py", "check"] + [str(root / "export" / ds) for ds in ST_DS]
    p = subprocess.run([sys.executable, *args], cwd=str(REPO), stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                       text=True, encoding="utf-8", errors="replace")
    with open(log_path, "a", encoding="utf-8") as lg:
        lg.write(p.stdout)
    c.eq(p.returncode, 0, "deid.py check on the stereoseq exports (exit code)")
    return res



# ------------------------------------------------------------------ standard formats (examples/make_stdformats_demo.py, stdformats.py)
# a spatial AnnData .h5ad (Squidpy layout, classic Visium) / an AnnData-Zarr of cells without an image / an OME-Zarr
# image / a SpatialData store shaped like spatialdata-io's Xenium / one shaped like a Visium HD 16 µm layer
SA_AV, SA_AC, SA_OZ, SA_SX, SA_SH = "demo_anndata_visium", "demo_anndata_cells", "demo_omezarr", "demo_sdata_xenium", "demo_sdata_visiumhd"
# zarr v3 (zarrv3.py): the OME-Zarr image as OME-NGFF 0.5 (sharded, transpose, big-endian, v2 chunk keys) and the
# cells AnnData-Zarr again
SA_OZ3, SA_AC3 = "demo_omezarr_v3", "demo_anndata_cells_v3"
STD_DS = (SA_AV, SA_AC, SA_OZ, SA_SX, SA_SH, SA_OZ3, SA_AC3)
STD_FAKE = "INT" + "9961"           # the generator's made-up sample ID (library name, a column's values; deid:ok)
STD_CASES = "stdformats_cases"


def stdformats_commands(root: Path) -> list[list[str]]:
    try:
        r = root.resolve().relative_to(REPO.resolve()).as_posix()
    except ValueError:
        r = root.resolve().as_posix()
    raw, cache, export = f"{r}/raw", f"{r}/cache", f"{r}/export"
    prep, mk = "web/SpatialViewer/prepare.py", "examples/demo_xenium_markers.json"
    return [
        ["examples/make_stdformats_demo.py", "--out", raw, "--cases", f"{r}/{STD_CASES}"],
        [prep, "anndata", "--id", SA_AV, "--src", f"{raw}/anndata_visium.h5ad", "--out", cache, "--label", "Synthetic spatial AnnData"],
        [prep, "export", "--id", SA_AV, "--out", cache, "--dest", export, "--markers", mk],
        [prep, "anndata", "--id", SA_AC, "--src", f"{raw}/anndata_cells.zarr", "--out", cache, "--coords-um", "1"],
        # no UMAP: the file's leiden column has no embedding to sit on -> a message, the export goes on
        [prep, "export", "--id", SA_AC, "--out", cache, "--dest", export, "--no-umap"],
        [prep, "omezarr", "--id", SA_OZ, "--src", f"{raw}/image.ome.zarr", "--out", cache],
        [prep, "export", "--id", SA_OZ, "--out", cache, "--dest", export, "--markers", "examples/demo_protein_markers.json"],
        [prep, "spatialdata", "--id", SA_SX, "--src", f"{raw}/sdata_xenium.zarr", "--out", cache, "--pixel-size", "0.5"],
        [prep, "export", "--id", SA_SX, "--out", cache, "--dest", export, "--markers", mk],
        [prep, "spatialdata", "--id", SA_SH, "--src", f"{raw}/sdata_visiumhd.zarr", "--out", cache, "--pixel-size", "2"],
        [prep, "export", "--id", SA_SH, "--out", cache, "--dest", export, "--no-umap"],
        [prep, "omezarr", "--id", SA_OZ3, "--src", f"{raw}/image_v3.ome.zarr", "--out", cache],
        [prep, "export", "--id", SA_OZ3, "--out", cache, "--dest", export, "--no-umap"],
        [prep, "anndata", "--id", SA_AC3, "--src", f"{raw}/anndata_cells_v3.zarr", "--out", cache, "--coords-um", "1"],
        [prep, "export", "--id", SA_AC3, "--out", cache, "--dest", export, "--no-umap"],
    ]


def clean_stdformats(root: Path):
    for p in [root / STD_CASES] + [root / sub / ds for sub in ("cache", "export") for ds in STD_DS]:
        if p.exists():
            print(f"  removing earlier output {p}", flush=True)
            shutil.rmtree(p)


def _std_layer(c: Checker, root: Path, ds: str, key: str, want_px: np.ndarray, counts: np.ndarray, genes: list,
               what: str) -> np.ndarray:
    """An expression layer of a standard-format export against the fixture: positions (image pixels, table
    order), genes, totals and every gene file."""
    d = root / "export" / ds / "expr" / key
    pos = f32(d / "positions.bin", 2)
    c.ok(pos.shape == want_px.shape and np.allclose(pos, want_px, atol=2e-2),
         f"{ds}: {key} positions.bin = {what}, in the file's row order (max off "
         f"{float(np.abs(pos - want_px).max()) if pos.shape == want_px.shape else 'shape ' + str(pos.shape)})")
    gj = read_json(d / "genes.json")
    c.eq([g["n"] for g in gj], list(genes), f"{ds}: {key} genes.json (features that are not genes and controls left out)")
    c.ok(np.allclose([g["t"] for g in gj], counts.sum(0)), f"{ds}: {key} genes.json totals = the matrix's column sums")
    bad = []
    for j, g in enumerate(genes):
        rows = np.flatnonzero(counts[:, j] > 0)
        a = gene_f32(d, g)
        a = np.zeros((0, 3), np.float32) if a is None else a
        if len(a) != len(rows) or not np.allclose(a[:, :2], pos[rows]) or not np.allclose(a[:, 2], counts[rows, j]):
            bad.append(g)
    c.ok(not bad, f"{ds}: {key} gene files = the matrix (differ: {bad[:4]})")
    return pos


def _check_zarr_v3(c: Checker, root: Path, tr, tj) -> dict:
    """The zarr v3 fixtures (zarrv3.py): the OME-NGFF 0.5 image = the v2 image, the AnnData-Zarr v3 = the v2 one."""
    raw, export, cache = root / "raw", root / "export", root / "cache"
    out = {}
    for ds in (SA_OZ3, SA_AC3):
        if not c.ok((export / ds / "manifest.json").exists(), f"{ds}: exported"):
            return out
    # ---- OME-NGFF 0.5 image (zarr v3: level 0 sharded + blosc, level 1 transpose + big-endian + gzip, "." keys)
    m3, m2 = read_json(export / SA_OZ3 / "manifest.json"), read_json(export / SA_OZ / "manifest.json")
    c.eq([(ch["name"], ch["color"]) for ch in m3["image"]["channels"]], [(ch["name"], ch["color"]) for ch in m2["image"]["channels"]],
         f"{SA_OZ3}: channel names and colours from ome.omero (NGFF 0.5 under attributes.ome)")
    c.eq((m3.get("pixelSizeUm"), m3["image"]["width"], m3["image"]["height"]), (0.5, m2["image"]["width"], m2["image"]["height"]),
         f"{SA_OZ3}: pixel size and size = the zarr v2 image's")
    try:
        import zarr
        img = zarr.open_group(str(cache / SA_OZ3 / "image.zarr"), mode="r")["0"][:]
        c.ok(np.array_equal(img, tr["oz_image"]), f"{SA_OZ3}: level 0 of the pyramid = the image, pixel for pixel (read from shards)")
    except ImportError:
        print("    SKIP the OME-Zarr v3 pixel check (zarr not installed)")
    ng, zv = app_module("ngff"), app_module("zarrv3")
    for lev in ("1", "2"):
        a = ng.open_array(raw / "image_v3.ome.zarr" / lev)
        b = tr["oz_image"][:, ::2 ** int(lev), ::2 ** int(lev)]
        c.ok(isinstance(a, zv.Array) and a.shape == b.shape and np.array_equal(a[...], b) and np.array_equal(a[1, 3:17, 5:40], b[1, 3:17, 5:40]),
             f"zarrv3: image_v3 level {lev} ({read_json(raw / 'image_v3.ome.zarr' / lev / 'zarr.json')['codecs'][0]['name']} ...) = the image, whole and a window")
    # ---- AnnData-Zarr v3 = the v2 one
    for f in ("positions.bin", "genes.json"):
        a, b = export / SA_AC3 / "expr" / "cells" / f, export / SA_AC / "expr" / "cells" / f
        c.ok(a.exists() and b.exists() and a.read_bytes() == b.read_bytes(), f"{SA_AC3}: {f} = the zarr v2 AnnData's, byte for byte")
    return out


def check_stdformats(c: Checker, root: Path, clustering: bool, log_path: Path) -> dict:
    """The standard formats (stdformats.py, annstore.py, ngff.py) against the generator's truth: positions through
    scale factors / NGFF transformations, counts, genes, images, outlines, molecules, the file's own clusterings,
    what the export must not hold, the cases that are refused, and the page in a browser."""
    raw, export, cache = root / "raw", root / "export", root / "cache"
    tr = np.load(raw / "truth_stdformats.npz")
    tj = read_json(raw / "truth_stdformats.json")
    genes = tj["genes"]
    out = {}
    for ds in STD_DS:
        if not c.ok((export / ds / "manifest.json").exists(), f"{ds}: exported"):
            return out
    # ---- spatial AnnData, classic Visium (uns spatial hires image, scale factors)
    m = read_json(export / SA_AV / "manifest.json")
    L = (m.get("expr") or [{}])[0]
    um = 55.0 / tj["spotPx"] / tj["hires"]
    c.ok(m.get("kind") == "anndata" and abs((m.get("pixelSizeUm") or 0) - um) < 1e-9,
         f"{SA_AV}: kind anndata, pixel size = 55 µm / spot_diameter_fullres / tissue_hires_scalef ({m.get('pixelSizeUm')})")
    c.eq([ch["name"] for ch in m["image"]["channels"]], ["Red", "Green", "Blue"], f"{SA_AV}: the hires image as an RGB photo")
    c.eq((m["image"]["width"], m["image"]["height"]), (1000, 600), f"{SA_AV}: image size = the hires image")
    c.ok("source" not in m and "userObs" not in m and not m.get("imaging"), f"{SA_AV}: no local paths / cache tables / imaging flag in the export manifest")
    c.eq((L.get("key"), L.get("platform"), L.get("nBins"), L.get("nGenes"), L.get("spotDiameterPx")),
         ("spots", "visium", 240, 20, tj["spotPx"] * tj["hires"]), f"{SA_AV}: spots layer (key, platform, rows, genes, dot)")
    _std_layer(c, root, SA_AV, "spots", tr["av_xy"] * tj["hires"], tr["av_counts"], genes, "obsm spatial x tissue_hires_scalef")
    c.ok((cache / SA_AV / "userobs" / "spots" / "obs.csv").exists() and not (export / SA_AV / "userobs").exists(),
         f"{SA_AV}: the obs tables stay in the cache")
    if clustering:
        idx = read_json(export / SA_AV / "expr" / "spots" / "umap" / "index.json")
        cls = {x.get("name"): x for x in idx.get("clusterings", []) if x.get("source") == "user"}
        cl = cls.get("cluster", {})
        c.ok(cl.get("categories") and sorted(cl["categories"]) == ["Region A", "Region B", "Region C"]
             and dict(zip(cl["categories"], cl.get("colors") or [])) == {"Region A": "#1f77b4", "Region B": "#ff7f0e", "Region C": "#2ca02c"},
             f"{SA_AV}: obs 'cluster' became your own clustering with uns cluster_colors ({cl.get('categories')}, {cl.get('colors')})")
        c.eq(idx.get("embedding", {}).get("source"), "user", f"{SA_AV}: obsm X_umap became your own UMAP")
        lab = np.fromfile(str(export / SA_AV / "expr" / "spots" / "umap" / cl.get("file", "missing")), dtype="<u2") if cl else np.zeros(0)
        if len(lab) == len(tr["av_types"]):
            c.ok(ari(lab, tr["av_types"]) > 0.99, f"{SA_AV}: your own clustering = the cluster column, row by row")
        else:
            c.ok(False, f"{SA_AV}: clusters_user rows {len(lab)}")
        smp = cls.get("sample", {})
        c.ok(len(smp.get("categories") or []) == 2 and not any(STD_FAKE in str(x) for x in smp.get("categories") or []),
             f"{SA_AV}: the sample column (made-up sample IDs) kept, with its values de-identified ({smp.get('categories')})")
        out["anndataUserClusters"] = sorted(cls)
    # ---- AnnData-Zarr of cells, no image: density background in µm (--coords-um 1)
    m = read_json(export / SA_AC / "manifest.json")
    L = (m.get("expr") or [{}])[0]
    c.ok(m["image"].get("background") == "cell density" and m.get("pixelSizeUm") == 1.0,
         f"{SA_AC}: a cell-density background at 1 µm per pixel ({m['image'].get('background')}, {m.get('pixelSizeUm')})")
    c.eq((L.get("key"), L.get("kind"), L.get("nBins")), ("cells", "cells", 300), f"{SA_AC}: cells layer")
    pos = f32(export / SA_AC / "expr" / "cells" / "positions.bin", 2)
    xy = tr["ac_xy"]
    c.ok(pos.shape == xy.shape and np.allclose(pos - pos.min(0), xy - xy.min(0), atol=2e-3),
         f"{SA_AC}: positions = obsm spatial in µm, shifted onto the background")
    lt = log_path.read_text(encoding="utf-8", errors="replace") if log_path.exists() else ""
    c.ok("the file's own clusterings / UMAP were not attached" in lt and not (export / SA_AC / "expr" / "cells" / "umap").exists(),
         f"{SA_AC}: export --no-umap: the leiden column has no embedding to sit on, a message and the export goes on")
    # ---- OME-Zarr image
    m = read_json(export / SA_OZ / "manifest.json")
    c.eq([(ch["name"], ch["color"]) for ch in m["image"]["channels"]],
         [("DAPI", "#0000ff"), ("CD45", "#ff0000"), ("PanCK", "#00ff00"), ("CD3", "#ffff00")],
         f"{SA_OZ}: channel names and colours from omero")
    c.eq(m.get("pixelSizeUm"), 0.5, f"{SA_OZ}: pixel size from the level-0 scale (micrometer axes)")
    try:
        import zarr
        img = zarr.open_group(str(cache / SA_OZ / "image.zarr"), mode="r")["0"][:]
        c.ok(np.array_equal(img, tr["oz_image"]), f"{SA_OZ}: level 0 of the pyramid = the OME-Zarr's level 0, pixel for pixel")
    except ImportError:
        print("    SKIP the OME-Zarr pixel check (zarr not installed)")
    c.ok(any(x.get("kind") == "protein" for x in m.get("expr", [])), f"{SA_OZ}: the protein grid layer, as for an OME-TIFF")
    # ---- SpatialData, Xenium-like: circles, outlines, transcripts, a second coordinate system
    m = read_json(export / SA_SX / "manifest.json")
    L = next((x for x in m.get("expr", []) if x.get("key") == "cells"), {})
    c.ok(m.get("imaging") is True and m.get("pixelSizeUm") == 0.5 and (m.get("spatialdata") or {}).get("coordinateSystem") == "global",
         f"{SA_SX}: fluorescence image (imaging), --pixel-size, the global coordinate system of two")
    c.eq((L.get("kind"), L.get("platform"), L.get("nBins"), L.get("nGenes")), ("cells", "spatialdata", 80, 20), f"{SA_SX}: cells layer")
    c.ok(abs((L.get("spotDiameterPx") or 0) - float(np.median(2 * tr["sx_rad_px"]))) < 1e-3,
         f"{SA_SX}: dot = the median circle diameter ({L.get('spotDiameterPx')})")
    pos = _std_layer(c, root, SA_SX, "cells", tr["sx_xy_px"], tr["sx_counts"], genes, "circle centres x scale (µm -> pixels)")
    sets = [s.get("key") for s in ((L.get("shapes") or {}).get("sets") or [])]
    c.eq(sets, ["cell", "nucleus"], f"{SA_SX}: cell outlines (cell_boundaries, scale) and nucleus outlines (nucleus_boundaries, affine)")
    cs = app_module("cellshapes")
    for key in sets:
        try:
            _i, rows, start, vx = cs.read_set(export / SA_SX / "expr" / "cells" / "shapes" / key)
            def centroid(r):                         # area centroid of a ring (the vertex mean moves when one is simplified away)
                x, y = r[:, 0], r[:, 1]
                x1, y1 = np.roll(x, -1), np.roll(y, -1)
                cr = x * y1 - x1 * y
                a = cr.sum() / 2.0
                return np.array([((x + x1) * cr).sum() / (6 * a), ((y + y1) * cr).sum() / (6 * a)])
            cen = np.array([centroid(vx[start[k]:start[k + 1]]) for k in range(len(rows)) if start[k + 1] > start[k]])
            rr = np.array([rows[k] for k in range(len(rows)) if start[k + 1] > start[k]])
            off = np.abs(cen - pos[rr]).max() if len(rr) else 99
            c.ok(len(rr) == 80 and off < 0.5, f"{SA_SX}: {key} outlines: 80 polygons, each around its own cell (max off {off:.3f} px)")
        except (OSError, ValueError, KeyError) as e:
            c.ok(False, f"{SA_SX}: {key} outlines readable ({e})")
    mol = m.get("molecules") or {}
    c.eq((mol.get("n"), mol.get("nGenes")), (int(tr["sx_molecules"][0]), 20),
         f"{SA_SX}: molecules = the transcripts with qv >= 20, controls left out (two parquet parts, dictionary and plain)")
    if clustering:
        idx = read_json(export / SA_SX / "expr" / "cells" / "umap" / "index.json")
        ct = next((x for x in idx.get("clusterings", []) if x.get("name") == "cell_type" and x.get("source") == "user"), {})
        c.ok(sorted(ct.get("categories") or []) == ["Type A", "Type B", "Type C"], f"{SA_SX}: the table's cell_type as your own clustering")
    # ---- SpatialData, Visium HD-like: equal squares, an image with a scale of 4
    m = read_json(export / SA_SH / "manifest.json")
    L = (m.get("expr") or [{}])[0]
    c.eq((L.get("key"), L.get("binSizeUm"), L.get("spotDiameterPx"), L.get("nBins")), ("016um", 16.0, 8.0, 140),
         f"{SA_SH}: equal squares = a 16 µm bin layer (side through the inverse image transformation)")
    c.ok(not m.get("imaging"), f"{SA_SH}: an RGB image is a photo")
    _std_layer(c, root, SA_SH, "016um", tr["sh_xy_px"], tr["sh_counts"], genes, "square centres / the image's scale")
    out.update(_check_zarr_v3(c, root, tr, tj))

    # ---- what is refused, and that no message names the made-up sample ID
    def cli(*argv):
        p = subprocess.run([sys.executable, str(APP / "prepare.py"), *argv], cwd=str(REPO), stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace",
                           env={**os.environ, "PYTHONIOENCODING": "utf-8"})
        with open(log_path, "a", encoding="utf-8") as lg:
            lg.write(f"\n$ prepare.py {' '.join(argv)}\n{p.stdout}")
        return p.returncode, p.stdout
    cases, never = root / STD_CASES, root / "std_never"
    for kind, name, say, extra in (("omezarr", "v3_badcodec.ome.zarr", ["zarr v3", "'zfp'", "zarr_write_format = 2"], []),
                                   ("spatialdata", "labels_3d.zarr", ["labels element", "not read yet"], []),
                                   ("omezarr", "plate.ome.zarr", ["plate", "one field"], []),
                                   ("anndata", "two_libraries.h5ad", ["2 image libraries", "--library"], []),
                                   ("anndata", "two_libraries.h5ad", ["not one of the 2 libraries"], ["--library", "zz"]),
                                   ("spatialdata", "labels_region.zarr", ["labels element", "not read yet"], []),
                                   ("spatialdata", "mapaxis.zarr", ["mapAxis", "identity, scale"], [])):
        rc, txt = cli(kind, "--src", str(cases / name), "--id", "x", "--out", str(never), *extra)
        c.ok(rc != 0 and all(w in txt for w in say) and "Traceback" not in txt and STD_FAKE not in txt,
             f"prepare.py {kind} {name} {' '.join(extra)}: refused with {say}, no traceback, no sample ID (exit {rc}: {txt[-240:]!r})")
    c.ok(not (never / "x" / "manifest.json").exists(), "standard formats: nothing prepared from a refused file")
    rc, txt = cli("anndata", "--src", str(cases / "two_libraries.h5ad"), "--id", "lib2", "--out", str(never), "--library", "2")
    m2 = read_json(never / "lib2" / "manifest.json") if (never / "lib2" / "manifest.json").exists() else {}
    c.ok(rc == 0 and ((m2.get("expr") or [{}])[0]).get("nBins") == 120 and STD_FAKE not in txt,
         f"prepare.py anndata --library 2: that section's 120 rows, no sample ID in the log (exit {rc})")
    # the cell ids of a SpatialData table: no id format to hold them to, so they are written when no sample ID is in them
    sx = export / SA_SX / "expr"
    lists = sorted(sx.rglob("ids.txt.gz"))
    before = {f: f.read_bytes() for f in lists}
    rc, txt = cli("ids", "--id", SA_SX, "--out", str(cache), "--dest", str(export))
    c.ok(lists and rc == 0 and "REFUSING" not in txt and all(f.read_bytes() == b for f, b in before.items())
         and not any(STD_FAKE in gzip.decompress(f.read_bytes()).decode("utf-8") for f in lists),
         f"SpatialData cell ids: in the default export (no id format, no sample ID in them), prepare.py ids writes the same (exit {rc}, {[f.parent.name for f in lists]})")
    lt = log_path.read_text(encoding="utf-8", errors="replace") if log_path.exists() else ""
    nl = chr(10)
    ours = [blk for blk in lt.split(nl + "$ ") if any(k in blk.split(nl, 1)[0] for k in STD_DS + (STD_CASES, "make_stdformats_demo"))]
    c.ok(len(ours) >= 10 and not [blk.split(nl, 1)[0] for blk in ours if STD_FAKE in blk],
         f"standard formats: no prepare / export message holds the made-up sample ID (library name, column values; {len(ours)} commands read)")
    texts = [p for ds in STD_DS for p in (export / ds).rglob("*") if p.is_file() and p.suffix in (".json", ".geojson", ".txt")]
    c.ok(not [p for p in texts if STD_FAKE in p.read_text(encoding="utf-8", errors="ignore")],
         "standard formats: no export text file holds the made-up sample ID")

    # ---- the page, in Chromium and Firefox (SKIPPED without node / Playwright; Firefox SKIPPED when not installed)
    node = shutil.which("node")
    if not node:
        print("    node not found: tests/stdformats_browser_test.js skipped")
    else:
        p = subprocess.run([node, str(REPO / "tests" / "stdformats_browser_test.js"), "--export", str(export),
                            "--shots", str(root / "stdformats_shots"), "--browsers", "chromium,firefox"], cwd=str(REPO),
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace")
        with open(log_path, "a", encoding="utf-8") as lg:
            lg.write(p.stdout)
        print_skipped(p.stdout)
        print("    " + (p.stdout.strip().splitlines() or ["(no output)"])[-1])
        c.eq(p.returncode, 0, "node tests/stdformats_browser_test.js: the standard-format exports open in Chromium and Firefox (exit code)")
        out["browser"] = (p.stdout.strip().splitlines() or [""])[-1]
    return out


# ------------------------------------------------------------------ automatic detection (detect.py)
# every synthetic folder the generators above write: the kind and variant detect.py must name, with verdict "certain"
DETECT_WANT = {
    "visiumhd": ("visium", "visium-hd"), "multiplex.ome.tif": ("tiff", "ome-tiff"), "multiplex_tiled.ome.tif": ("tiff", "ome-tiff"),
    "xenium": ("xenium", "cells+molecules"), "xenium_pq": ("xenium", "molecules"), "xenium_nt": ("xenium", "cells"),
    "xenium_noimg": ("xenium", "cells+molecules"), "xenium_output.zip": ("xenium", "cells+molecules"),
    "xenium_bundle.zip": ("xenium", "cells+molecules"),
    "geomx": ("geomx", "scans"), "cosmx": ("cosmx", "images+cells+molecules"), "cosmx_new": ("cosmx", "images+cells+molecules"),
    "cosmx_tx": ("cosmx", "molecules"), "merscope": ("merscope", "images+cells+molecules"),
    "merscope_noimg": ("merscope", "cells+molecules"), "merscope_pq": ("merscope", "images+molecules"),
    "stereoseq": ("stereoseq", "bins+cells"),
    "anndata_visium.h5ad": ("anndata", "h5ad"), "anndata_cells.zarr": ("anndata", "zarr"), "image.ome.zarr": ("omezarr", "zarr"),
    "sdata_xenium.zarr": ("spatialdata", "zarr"), "sdata_visiumhd.zarr": ("spatialdata", "zarr"),
    "image_v3.ome.zarr": ("omezarr", "zarr-v3"), "anndata_cells_v3.zarr": ("anndata", "zarr-v3"),
}
# `prepare.py auto` runs whose cache must equal the explicit subcommand's: (id, source under raw/, options passed through)
LONE_STORE, XST = "demo_transcripts.zarr.zip", "demo_xenium_store"   # a Xenium transcripts store on its own
AUTO_SAME = ((VIS, "visiumhd", []), (XE, "xenium", ["--label", "Synthetic Xenium"]), (MX, "multiplex.ome.tif", []),
             (XZ, "xenium_output.zip", ["--label", "Synthetic Xenium"]),
             (CXT, "cosmx_tx", ["--pixel-size", "0.5"]), (MSN, "merscope_noimg", []),
             (SA_AV, "anndata_visium.h5ad", ["--label", "Synthetic spatial AnnData"]), (SA_OZ, "image.ome.zarr", []),
             (SA_SX, "sdata_xenium.zarr", ["--pixel-size", "0.5"]))


def _touch(p: Path, data: bytes = b"x") -> Path:
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_bytes(data)
    return p


def make_detect_cases(root: Path) -> Path:
    """Small folders detection must NOT call certain (and a few it must), built from stub files:
    detection reads names, sizes and headers, so most files only need to exist."""
    raw, neg = root / "raw", root / "detect" / "cases"
    if neg.exists():
        shutil.rmtree(neg)
    (neg / "empty").mkdir(parents=True)

    def xenium(d: Path):
        _touch(d / "experiment.xenium", b"{}")
        _touch(d / "morphology_focus" / "ch0000_dapi.ome.tif")
        _touch(d / "transcripts.parquet")

    def spaceranger(d: Path):
        for b in ("square_008um", "square_016um"):
            _touch(d / "outs" / "binned_outputs" / b / "filtered_feature_bc_matrix.h5")
            _touch(d / "outs" / "binned_outputs" / b / "spatial" / "tissue_positions.parquet")

    def matrix(d: Path):
        for n in ("matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz"):
            _touch(d / "filtered_feature_bc_matrix" / n)

    xenium(neg / "two_platforms")                                  # Xenium and MERSCOPE files in one folder
    for n in ("detected_transcripts.csv", "cell_by_gene.csv", "cell_metadata.csv"):
        _touch(neg / "two_platforms" / n, b"a,b\n")
    spaceranger(neg / "visium_no_image")
    shutil.copytree(raw / "visiumhd", neg / "visium_two_images", ignore=shutil.ignore_patterns("truth_*"))
    shutil.copyfile(raw / "visiumhd" / "image.tif", neg / "visium_two_images" / "image_again.tif")
    matrix(neg / "tenx_cells")
    _touch(neg / "tenx_cells" / "umap.csv", b"barcode,UMAP1,UMAP2\n")
    _touch(neg / "tenx_cells" / "clusters.csv", b"barcode,cluster\n")
    _touch(neg / "tenx_cells" / "notes.csv", b"gene,comment\n")
    matrix(neg / "tenx_no_umap")
    matrix(neg / "tenx_spatial")
    _touch(neg / "tenx_spatial" / "spatial" / "tissue_positions.csv")
    _touch(neg / "tenx_spatial" / "spatial" / "scalefactors_json.json", b"{}")
    (neg / "zipped").mkdir()
    with zipfile.ZipFile(neg / "zipped" / "run1.zip", "w") as z:
        for n in ("run1/experiment.xenium", "run1/morphology_focus/ch0000_dapi.ome.tif", "run1/cells.parquet",
                  "run1/cell_feature_matrix.h5", "run1/transcripts.parquet"):
            z.writestr(n, "x")
    _touch(neg / "xenium_no_stains" / "transcripts.zarr.zip")
    _touch(neg / "xenium_old_layout" / "morphology.ome.tif")          # Xenium 1.x: one image file, no morphology_focus/
    _touch(neg / "xenium_old_layout" / "transcripts.zarr.zip")
    (neg / "zipped_merscope").mkdir()
    with zipfile.ZipFile(neg / "zipped_merscope" / "region0.zip", "w") as z:
        for n in ("region_0/detected_transcripts.csv", "region_0/cell_by_gene.csv", "region_0/cell_metadata.csv"):
            z.writestr(n, "a,b\n")
    (neg / "zip_and_folder").mkdir()                                 # run1.zip with its unzipped copy run1/ beside it
    shutil.copyfile(neg / "zipped" / "run1.zip", neg / "zip_and_folder" / "run1.zip")
    xenium(neg / "zip_and_folder" / "run1")
    store = raw / "xenium" / "transcripts.zarr.zip"                 # Xenium transcripts stores on their own
    (neg / "xenium_lone_store").mkdir()
    shutil.copyfile(store, neg / "xenium_lone_store" / LONE_STORE)
    _touch(neg / "xenium_store_bad" / "bad_transcripts.zarr.zip", b"not a zip")
    xenium(neg / "xenium_extra_store")                               # one inside a Xenium folder: part of the folder
    shutil.copyfile(store, neg / "xenium_extra_store" / "copy_transcripts.zarr.zip")
    (neg / "zip_and_folder_merscope").mkdir()                        # a zip that must be unzipped, beside its folder
    shutil.copyfile(neg / "zipped_merscope" / "region0.zip", neg / "zip_and_folder_merscope" / "region0.zip")
    for n in ("detected_transcripts.csv", "cell_by_gene.csv", "cell_metadata.csv"):
        _touch(neg / "zip_and_folder_merscope" / "region0" / n, b"a,b\n")
    (neg / "zipped_slip").mkdir()                                     # a Xenium zip with a name that climbs out (zip slip)
    shutil.copyfile(neg / "zipped" / "run1.zip", neg / "zipped_slip" / "slip.zip")
    with zipfile.ZipFile(neg / "zipped_slip" / "slip.zip", "a") as z:
        z.writestr("../outside.txt", "x")
    (neg / "two_zips").mkdir()
    for n in ("a.zip", "b.zip"):
        shutil.copyfile(neg / "zipped" / "run1.zip", neg / "two_zips" / n)
    xenium(neg / "xenium_extra_tiffs")
    for n in ("GeneA.tif", "GeneB.tif"):
        _touch(neg / "xenium_extra_tiffs" / n)
    _touch(neg / "stereoseq_no_image" / "CHIP9.tissue.gef")
    _touch(neg / "merscope_two_runs" / "detected_transcripts.csv", b"a,b\n")
    for run_ in ("vpt-run-a", "vpt-run-b"):
        for n in ("cell_by_gene.csv", "cell_metadata.csv"):
            _touch(neg / "merscope_two_runs" / run_ / n, b"a,b\n")
    _touch(neg / "not_data" / "notes.docx")
    shutil.copytree(raw / "geomx", neg / "geomx_two_versions", ignore=shutil.ignore_patterns("*.ome.tiff"))
    shutil.copyfile(raw / "geomx" / "pkc" / "Demo_Panel_v1.0.pkc", neg / "geomx_two_versions" / "pkc" / "Demo_Panel_v1.1_copy.pkc")
    try:                                                           # cases that need a real header
        import tifffile
        grey = np.zeros((50, 100), np.uint8)
        for d, names in (("grey_unnamed", ("DAPI.tif", "CD45.tif", "PanCK.tif")),
                         ("grey_markers", ("01_DAPI_s_gray.tif", "02_CD45_s_gray.tif", "03_PanCK_s_gray.tif"))):
            (neg / d).mkdir()
            for n in names:
                tifffile.imwrite(str(neg / d / n), grey)
        # Space Ranger recorded a 200 x 100 image (hires 100 x 50 at scale 0.5): of two TIFFs only one is that size
        d = neg / "visium_size_check"
        spaceranger(d)
        sp = d / "outs" / "binned_outputs" / "square_008um" / "spatial"
        (sp / "scalefactors_json.json").write_text(json.dumps({"tissue_hires_scalef": 0.5}), encoding="utf-8")
        (sp / "tissue_hires_image.png").write_bytes(b"\x89PNG\r\n\x1a\n" + struct.pack(">I4sII", 13, b"IHDR", 100, 50) + b"\x08\x02\x00\x00\x00")
        tifffile.imwrite(str(d / "section.tif"), np.zeros((100, 200, 3), np.uint8), photometric="rgb")
        tifffile.imwrite(str(d / "section_small_copy.tif"), np.zeros((50, 100, 3), np.uint8), photometric="rgb")
    except ImportError:
        pass
    try:
        import h5py
        for n, key in (("with_umap.h5ad", "X_umap"), ("no_umap.h5ad", "X_pca")):
            (neg / "h5ad").mkdir(exist_ok=True)
            with h5py.File(str(neg / "h5ad" / n), "w") as f:
                f.create_dataset("X", data=np.zeros((3, 2), np.float32))
                f.create_group("obs")
                f.create_group("obsm").create_dataset(key, data=np.zeros((3, 2), np.float32))
    except ImportError:
        pass
    return neg


def detect_commands(root: Path) -> list[list[str]]:
    """`prepare.py auto` on five of the synthetic folders (into detect/cache, same ids as the explicit runs), and
    once through --kind / --image on the folder with two candidate images."""
    try:
        r = root.resolve().relative_to(REPO.resolve()).as_posix()
    except ValueError:
        r = root.resolve().as_posix()
    prep = "web/SpatialViewer/prepare.py"
    cmds = [[prep, "auto", "--id", ds, "--src", f"{r}/raw/{src}", "--out", f"{r}/detect/cache", *extra] for ds, src, extra in AUTO_SAME]
    two = f"{r}/detect/cases/visium_two_images"
    cmds.append([prep, "auto", "--id", VIS, "--src", two, "--out", f"{r}/detect/cache_kind", "--kind", "visium",
                 "--image", f"{two}/image.tif"])
    cmds.append([prep, "auto", "--id", XST, "--src", f"{r}/detect/cases/xenium_lone_store/{LONE_STORE}", "--out",
                 f"{r}/detect/cache"])
    return cmds


def clean_detect(root: Path):
    # never/: what an earlier run's `prepare.py auto` must not have written (left over, it would fail this run too)
    for p in (root / "detect" / "cache", root / "detect" / "cache_kind", root / "detect" / "never"):
        if p.exists():
            print(f"  removing earlier output {p}", flush=True)
            shutil.rmtree(p)


def _same_cache(c: Checker, a: Path, b: Path, what: str, source_differs: bool = False):
    """Two cache folders of one dataset hold the same thing: manifests equal (creation time and the spelling of
    source paths aside), the same files, and the same values in every array. Caches of two ids (a.name, b.name)
    may differ by the id alone."""
    if not c.ok((a / "manifest.json").exists() and (b / "manifest.json").exists(), f"{what}: both caches have a manifest"):
        return

    def norm(v):
        if isinstance(v, dict):
            return {k: norm(x) for k, x in v.items()}
        if isinstance(v, list):
            return [norm(x) for x in v]
        if isinstance(v, str) and v and ("/" in v or "\\" in v):
            q = Path(v) if Path(v).is_absolute() else REPO / v
            if q.exists():
                return os.path.normcase(str(q.resolve()))
        return v
    ma, mb = read_json(a / "manifest.json"), read_json(b / "manifest.json")
    if a.name != b.name and mb.get("id") == b.name:
        mb["id"] = a.name
    for m in (ma, mb):
        m.pop("created", None)
        if source_differs:
            m.pop("source", None)
    c.ok(norm(ma) == norm(mb), f"{what}: the manifest is the explicit subcommand's (kind {ma.get('kind')!r} / {mb.get('kind')!r})")
    fa = sorted(q.relative_to(a).as_posix() for q in a.rglob("*") if q.is_file() and q.name != "manifest.json")
    fb = sorted(q.relative_to(b).as_posix() for q in b.rglob("*") if q.is_file() and q.name != "manifest.json")
    if not c.ok(fa == fb and len(fa) > 3, f"{what}: the same {len(fa)} cache files ({len(fb)} from the explicit subcommand)"):
        return
    # Blosc does not pack a chunk into the same bytes on every run, so array chunks are compared by their
    # values; everything else (array metadata, attributes, molecule files, overlays) byte for byte
    arrays = sorted(q.parent.relative_to(a).as_posix() for q in a.rglob(".zarray"))
    chunks = {n for n in fa for d in arrays if n.startswith(d + "/") and not n.rsplit("/", 1)[1].startswith(".z")}
    diff = [n for n in fa if n not in chunks
            and (a / n).read_bytes() != (b / n).read_bytes().replace(f'"{b.name}"'.encode(), f'"{a.name}"'.encode())]
    c.ok(not diff, f"{what}: array metadata, attributes and the other files are byte-identical to the explicit subcommand's (differ: {diff[:5]})")
    try:
        import zarr
    except ImportError:
        print(f"    ({what}: zarr not importable here, array values not compared)")
        return
    bad = [d for d in arrays if not np.array_equal(zarr.open(str(a / d), mode="r")[:], zarr.open(str(b / d), mode="r")[:])]
    c.ok(len(arrays) > 0 and not bad, f"{what}: all {len(arrays)} arrays (image levels, counts, positions, barcodes) hold the explicit subcommand's values (differ: {bad[:5]})")


def check_detect(c: Checker, root: Path, log_path: Path) -> dict:
    """detect.py on every synthetic folder, on folders it must not call certain, and `prepare.py auto`."""
    sys.path.insert(0, str(APP))
    try:
        import detect
    finally:
        sys.path.pop(0)
    raw = root / "raw"
    out = {"certain": 0, "seconds": 0.0}
    for name, (kind, variant) in DETECT_WANT.items():
        r = detect.detect(raw / name)
        out["seconds"] += r["scan"]["seconds"]
        if c.ok(r["verdict"] == "certain" and r["kind"] == kind and r["variant"] == variant,
                f"detect raw/{name}: certain {kind} / {variant} (got {r['verdict']} {r['kind']} / {r['variant']}: {r['message']})"):
            out["certain"] += 1
        c.ok(len(r["candidates"]) == 1 and r["candidates"][0]["argv"][:1] == [kind] and r["candidates"][0]["evidence"]
             and not r["candidates"][0]["missing"], f"detect raw/{name}: one candidate with evidence and the {kind} command")
    argv = {n: detect.detect(raw / n)["candidates"][0]["args"] for n in ("visiumhd", "stereoseq", "geomx", "merscope_noimg")}
    c.ok(Path(argv["visiumhd"].get("image", "")) == raw / "visiumhd" / "image.tif" and Path(argv["visiumhd"].get("outs", "")) == raw / "visiumhd" / "outs",
         f"detect visiumhd: --image is the TIFF next to outs/, --outs the Space Ranger folder: {argv['visiumhd']}")
    st = {k: Path(v).name for k, v in argv["stereoseq"].items()}
    c.eq(st, {"image": "SYNCHIP_HE_regist.tif", "gef": "SYNCHIP.tissue.gef", "cellbin": "SYNCHIP.cellbin.gef",
              "tissue_mask": "SYNCHIP_HE_tissue_cut.tif"}, "detect stereoseq: H&E before ssDNA, tissue.gef before raw.gef, the GEF before a GEM")
    c.ok(Path(argv["geomx"].get("dcc", "")).name == "dcc" and len(argv["geomx"].get("pkc", [])) == 2 and len(argv["geomx"].get("image", [])) == 2
         and {Path(a).name for a in argv["geomx"].get("annotation", [])} >= {"annotations.xlsx", "Demo_LabWorksheet.txt"},
         f"detect geomx: dcc folder, both PKC files, the sheets with Sample_ID, both scans: {argv['geomx']}")
    c.eq(list(argv["merscope_noimg"]), ["src"], "detect merscope_noimg: one segmentation run one folder down needs no --cells-dir")
    try:                                                 # the layer list is the reader's own
        sys.path.insert(0, str(APP))
        from expr import find_visium_layers
        want = [layer["key"] for layer in find_visium_layers(raw / "visiumhd" / "outs")]
        c.eq(detect.detect(raw / "visiumhd")["candidates"][0]["info"]["layers"], want, "detect visiumhd: layers = expr.find_visium_layers")
    except ImportError:
        print("    (expr.py not importable here: layer list not compared)")
    finally:
        sys.path.pop(0)
    r = detect.detect(raw)                               # the folder of all of them: each reported, none merged
    got = sorted(Path(d["path"]).name for d in r["datasets"])
    c.ok(r["verdict"] == "ambiguous" and r["reason"] == "several-datasets" and got == sorted(DETECT_WANT),
         f"detect raw/: {len(DETECT_WANT)} datasets reported one by one (got {r['reason']}: {got})")
    c.ok(all(d["verdict"] == "certain" and (d["kind"], d["variant"]) == DETECT_WANT[Path(d["path"]).name] for d in r["datasets"]
             if Path(d["path"]).name in DETECT_WANT), "detect raw/: every dataset of the folder is certain, with its own kind")

    # ---- folders that must not be called certain
    neg = root / "detect" / "cases"
    if not neg.exists():
        make_detect_cases(root)

    def case(name, verdict, reason, want=None, say=(), **kw):
        r = detect.detect(neg / name, **kw)
        text = (r["message"] + " | " + " | ".join(r["suggestions"]) + " | " + detect.format_report(r)).lower()
        c.ok(r["verdict"] == verdict and r["reason"] == reason and (want is None or r["kind"] == want)
             and all(w.lower() in text for w in say),
             f"detect {name}{' ' + str(kw) if kw else ''}: {verdict} / {reason}" + (f" / {want}" if want else "")
             + (f", saying {list(say)}" if say else "") + f" (got {r['verdict']} / {r['reason']} / {r['kind']}: {r['message'][:160]})")
        return r
    case("empty", "unknown", "empty", say=["empty"])
    case("no_such_folder", "unknown", "not-found", say=["does not exist"])
    case("not_data", "unknown", "nothing-recognised", say=["--kind"])
    case("two_platforms", "ambiguous", "two-kinds", say=["--kind xenium", "--kind merscope"])
    case("two_platforms", "certain", "one-kind", "xenium", kind="xenium")
    case("two_platforms", "ambiguous", "kind-mismatch", hint="visium")
    case("two_platforms", "incomplete", "kind-not-found", "visium", say=["outs"], kind="visium")
    case("visium_no_image", "incomplete", "missing-files", "visium", say=["full-resolution image", "--image"])
    case("visium_no_image", "certain", "one-kind", "visium", overrides={"image": str(raw / "visiumhd" / "image.tif")})
    r = case("visium_two_images", "ambiguous", "open-choice", "visium", say=["--kind visium --image", "image.tif", "image_again.tif", "mb"])
    opts = [o for ch in (r["candidates"][0]["choices"] if r["candidates"] else []) for o in ch["options"]]
    c.ok(len(opts) == 2 and all(o["size"] > 1000 for o in opts), "detect visium_two_images: both candidate images listed with their sizes")
    case("tenx_no_umap", "incomplete", "missing-files", "rnaseq", say=["umap"])
    r = case("tenx_cells", "certain", "one-kind", "rnaseq")
    a = r["candidates"][0]["args"] if r["candidates"] else {}
    c.ok(Path(a.get("mtx", "")).name == "filtered_feature_bc_matrix" and Path(a.get("umap", "")).name == "umap.csv"
         and [Path(o).name for o in a.get("obs", [])] == ["clusters.csv"],
         f"detect tenx_cells: --mtx, --umap and the barcode CSV as --obs (not the CSV without a barcode column): {a}")
    case("tenx_spatial", "incomplete", "missing-files", "visium", say=["filtered_feature_bc_matrix.h5", "full-resolution image"])
    r = case("zipped", "certain", "one-kind", "xenium", say=["read from the zip as it is", "run1.zip"])
    c.ok(r["candidates"] and r["candidates"][0]["argv"][:2] == ["xenium", "--src"]
         and Path(r["candidates"][0]["args"].get("src", "")).name == "run1.zip", "detect zipped: a Xenium zip is --src itself")
    case("zipped/run1.zip", "certain", "one-kind", "xenium", say=["zipped", "never next to the zip"])
    case("zipped_merscope", "incomplete", "zipped", "merscope", say=["unzip", "region0.zip"])
    r = case("zip_and_folder", "certain", "one-kind", "xenium", say=["run1.zip", "not prepared on its own"])
    c.ok(r["candidates"] and Path(r["candidates"][0]["args"].get("src", "")).name == "run1",
         "detect zip_and_folder: the unzipped folder is prepared, the zip beside it noted")
    sys.path.insert(0, str(APP))
    try:
        import batch
    finally:
        sys.path.pop(0)
    rows, warn = batch.find_datasets(neg / "zip_and_folder", {}, None)
    c.ok([(r["id"], r["verdict"]) for r in rows] == [("run1", "certain"), ("run1_zip", "certain")]
         and Path(rows[1]["src"]).name == "run1.zip" and any("may be the same data" in w and "run1_zip" in w for w in warn),
         f"batch zip_and_folder: run1/ and run1.zip beside it are two datasets, run1 and run1_zip, with a warning that "
         f"they may be the same data: {[r['id'] for r in rows]} {warn}")
    rows, warn = batch.find_datasets(neg / "zip_and_folder_merscope", {}, None)
    c.ok([r["id"] for r in rows] == ["region0"] and any("region0.zip is not prepared" in w for w in warn),
         f"batch zip_and_folder_merscope: a zip that must be unzipped, beside its folder, is left out with a warning: "
         f"{[r['id'] for r in rows]} {warn}")
    rows, warn = batch.find_datasets(neg / "xenium_lone_store", {}, None)
    c.ok([(r["id"], r["kind"], r["variant"], r["verdict"]) for r in rows] == [("demo_transcripts", "xenium", "molecules", "certain")],
         f"batch: a loose *_transcripts.zarr.zip is a Xenium dataset of its own: {[(r['id'], r['kind'], r['verdict']) for r in rows]}")
    for where in (f"xenium_lone_store/{LONE_STORE}", "xenium_lone_store"):
        r = case(where, "certain", "one-kind", "xenium", say=["density background", "0.2125", "no cells layer"])
        c.ok(r["variant"] == "molecules" and r["candidates"] and r["candidates"][0]["argv"][:2] == ["xenium", "--src"]
             and Path(r["candidates"][0]["args"].get("src", "")).name == LONE_STORE,
             f"detect {where}: a transcripts store on its own is --src itself, molecules only")
    case("xenium_store_bad/bad_transcripts.zarr.zip", "incomplete", "missing-files", "xenium", say=["not a readable zip"])
    r = case("xenium_extra_store", "certain", "one-kind", "xenium", say=["copy_transcripts.zarr.zip", "not prepared on their own"])
    c.ok(r["candidates"] and Path(r["candidates"][0]["args"].get("src", "")).name == "xenium_extra_store",
         "detect xenium_extra_store: a *_transcripts.zarr.zip inside a Xenium folder stays part of it (noted)")
    case("zipped_slip/slip.zip", "incomplete", "zipped", "xenium", say=["not read as it is", "'..'", "refused"])
    r = case("two_zips", "ambiguous", "several-datasets", say=["a.zip", "b.zip"])
    c.ok(len(r["datasets"]) == 2 and all(x["verdict"] == "certain" for x in r["datasets"]),
         "detect two_zips: two zipped Xenium outputs are two datasets, each certain")
    r = case("xenium_no_stains", "certain", "one-kind", "xenium", say=["no stain images", "density background", "0.2125"])
    c.eq(r["variant"], "molecules", "detect xenium_no_stains: variant")
    case("xenium_old_layout", "incomplete", "missing-files", "xenium", say=["morphology_focus", "--kind tiff"])
    r = case("xenium_extra_tiffs", "certain", "one-kind", "xenium", say=["2 other tiff", "not used"])
    c.ok(r["variant"] == "molecules" and any("cells layer" in n for n in (r["candidates"][0]["notes"] if r["candidates"] else [])),
         "detect xenium_extra_tiffs: stains + transcripts without cells is certain, and says there is no cells layer")
    case("stereoseq_no_image", "incomplete", "missing-files", "stereoseq", say=["regist", "--image"])
    case("merscope_two_runs", "ambiguous", "open-choice", "merscope", say=["--cells-dir", "vpt-run-a", "vpt-run-b"])
    case("geomx_two_versions", "ambiguous", "open-choice", "geomx", say=["demo panel", "--kind geomx --pkc", "Demo_Panel_v1.1_copy.pkc"])
    if (neg / "grey_unnamed").exists():
        case("grey_unnamed", "ambiguous", "two-kinds", say=["--kind akoya", "--src"])
        case("grey_unnamed", "certain", "one-kind", "akoya", kind="akoya")
        case("grey_markers", "certain", "one-kind", "akoya")
        r = case("grey_unnamed/DAPI.tif", "certain", "one-kind", "tiff")
        c.eq(r["variant"], "single-channel", "detect a lone grey TIFF: variant")
        c.eq(detect.detect(raw / "visiumhd" / "image.tif")["variant"], "he", "detect a lone RGB TIFF: variant he")
        r = case("visium_size_check", "certain", "one-kind", "visium", say=["space ranger recorded"])
        c.eq(Path((r["candidates"][0]["args"] if r["candidates"] else {}).get("image", "")).name, "section.tif",
             "detect visium_size_check: of two TIFFs, the one whose size Space Ranger recorded is the image")
    else:
        print("    (tifffile not importable here: the cases that need a TIFF header are skipped)")
    if (neg / "h5ad").exists():
        case("h5ad/with_umap.h5ad", "certain", "one-kind", "rnaseq")
        case("h5ad/no_umap.h5ad", "incomplete", "missing-files", "rnaseq", say=["x_umap", "x_pca", "--embedding"])
        case("h5ad/no_umap.h5ad", "certain", "one-kind", "rnaseq", overrides={"embedding": "X_pca"})
        case("h5ad", "ambiguous", "several-datasets")
    std = root / STD_CASES                               # standard formats (stdformats.py) that are refused or asked about
    if std.exists():
        r = detect.detect(std / "v3_badcodec.ome.zarr")
        c.ok(r["verdict"] == "incomplete" and r["kind"] == "omezarr" and "zarr v3" in r["message"] and "'zfp'" in r["message"],
             f"detect a zarr v3 OME-Zarr with a codec SODA does not read: incomplete, names it (got {r['verdict']} {r['kind']}: {r['message'][:120]})")
        r = detect.detect(std / "plate.ome.zarr")
        c.ok(r["verdict"] == "incomplete" and "plate" in r["message"], f"detect an HCS plate: incomplete ({r['verdict']}: {r['message'][:120]})")
        for name in ("labels_region.zarr", "labels_3d.zarr"):             # zarr v2 / zarr v3: refused as prepare refuses it
            r = detect.detect(std / name)
            c.ok(r["verdict"] == "incomplete" and r["kind"] == "spatialdata" and "labels element" in r["message"]
                 and "not read yet" in r["message"],
                 f"detect a SpatialData table on a labels element ({name}): incomplete, says labels are not read "
                 f"({r['verdict']} {r['kind']}: {r['message'][:120]})")
        r = detect.detect(std / "two_libraries.h5ad")
        text = " | ".join(r["suggestions"]) + " | " + detect.format_report(r)
        c.ok(r["verdict"] == "ambiguous" and r["reason"] == "open-choice" and r["kind"] == "anndata" and "--library 2" in text
             and STD_FAKE not in text, f"detect two libraries in one AnnData: a choice by number, no library name ({r['verdict']}: {r['message'][:120]})")
        r = detect.detect(std / "two_libraries.h5ad", overrides={"library": "2"})
        c.eq(r["verdict"], "certain", "detect two libraries with library 2 given: certain")
        c.ok(any('"library": r"2"' in ln for ln in detect.option_lines(detect.detect(std / "two_libraries.h5ad"), "two_libraries")),
             "detect two libraries: the notebook's option line names library")
        r = detect.detect(std / "two_libraries.h5ad", kind="rnaseq")
        c.eq((r["kind"], r["verdict"]), ("rnaseq", "certain"), "detect a spatial .h5ad with --kind rnaseq: the single-cell reading")
    if (neg / "h5ad").exists():
        r = detect.detect(neg / "h5ad" / "with_umap.h5ad", kind="anndata")
        c.ok(r["verdict"] == "incomplete" and "obsm['spatial']" in r["message"],
             f"detect an .h5ad without positions with --kind anndata: incomplete, says obsm['spatial'] ({r['message'][:120]})")
    out["cases"] = c.n

    # ---- the command line
    def cli(*argv):
        p = subprocess.run([sys.executable, str(APP / "prepare.py"), *argv], cwd=str(REPO), stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace",
                           env={**os.environ, "PYTHONIOENCODING": "utf-8"})
        return p.returncode, p.stdout
    rc, txt = cli("detect", str(raw / "xenium"), "--json")
    try:
        j = json.loads(txt)
    except ValueError:
        j = {}
    c.ok(rc == 0 and j.get("verdict") == "certain" and j.get("kind") == "xenium" and j.get("candidates", [{}])[0].get("argv", [])[:2] == ["xenium", "--src"],
         f"prepare.py detect --json: a report a script can read (exit {rc})")
    rc, txt = cli("detect", str(raw / "visiumhd"))
    c.ok(rc == 0 and "CERTAIN" in txt and "prepare.py auto --src" in txt and "prepare.py visium --image" in txt,
         "prepare.py detect: plain words, the auto command and the command it would run")
    before = sorted(q.name for q in (root / "detect").iterdir())
    for name, code, say, extra in (("visium_two_images", 3, ["--kind visium --image", "ambiguous", "Nothing was written"], []),
                                   ("visium_no_image", 4, ["MISSING", "--image"], []),
                                   ("empty", 5, ["empty"], []),
                                   ("two_platforms", 3, ["--kind xenium", "--kind merscope"], []),
                                   ("two_platforms", 4, ["Visium", "outs"], ["--kind", "visium"]),
                                   ("zipped_merscope", 4, ["unzip", "region0.zip"], []),
                                   ("zipped_slip/slip.zip", 4, ["'..'", "refused"], [])):
        rc, txt = cli("auto", "--src", str(neg / name), "--id", "never_written", "--out", str(root / "detect" / "never"), *extra)
        c.ok(rc == code and all(w in txt for w in say),
             f"prepare.py auto {name} {' '.join(extra)}: stops with exit {code} and says {say} (exit {rc}: {txt[-300:]!r})")
    c.ok(sorted(q.name for q in (root / "detect").iterdir()) == before and not (root / "detect" / "never").exists(),
         "prepare.py auto: nothing is written when it stops")
    rc, txt = cli("auto", "--src", str(neg / "empty"), "--id", "x", "--out", str(root / "detect" / "never"), "--kind", "nonsense")
    c.ok(rc != 0 and "unknown kind" in txt, "prepare.py auto --kind nonsense: refused with the list of kinds")
    rc, txt = cli("xenium", "--src", str(raw / "xenium"), "--id", "x", "--out", str(root / "detect" / "never"), "--not-an-option")
    c.ok(rc == 2 and "unrecognized arguments" in txt and not (root / "detect" / "never").exists(),
         "prepare.py <kind>: an unknown option is still refused (only auto passes options through)")

    # ---- auto = the explicit subcommand
    log_txt = log_path.read_text(encoding="utf-8", errors="replace") if log_path.exists() else ""
    for ds, src, _extra in AUTO_SAME:
        _same_cache(c, root / "detect" / "cache" / ds, root / "cache" / ds, f"auto {src} -> {ds}")
    _same_cache(c, root / "detect" / "cache_kind" / VIS, root / "cache" / VIS, "auto --kind visium --image (two candidate images)",
                source_differs=True)
    if "prepare.py auto" in log_txt:
        c.ok("[auto] running: python prepare.py xenium --src" in log_txt and "(detected)" in log_txt and "(--kind visium)" in log_txt,
             "prepare.py auto: the log names what was detected and the command that ran")
    out["autoSame"] = len(AUTO_SAME) + 1
    ms = root / "detect" / "cache" / XST / "manifest.json"
    if c.ok(ms.exists(), f"prepare.py auto on a lone {LONE_STORE}: prepared ({ms})"):
        m, ref = read_json(ms), read_json(root / "cache" / XNI / "manifest.json") if (root / "cache" / XNI / "manifest.json").exists() else {}
        mol, rmol = m.get("molecules") or {}, ref.get("molecules") or {}
        c.ok(m.get("kind") == "xenium" and (m.get("image") or {}).get("background") == "transcript density"
             and not m.get("expr") and abs((m.get("pixelSizeUm") or 0) - 5 * 0.2125) < 1e-9
             and (m.get("frame") or {}).get("xeniumPixelSizeUm") == 0.2125,
             f"lone transcripts store: Xenium, a transcript-density background of 5 x 0.2125 um pixels (no experiment.xenium), "
             f"no cells layer (kind {m.get('kind')}, background {(m.get('image') or {}).get('background')}, {m.get('pixelSizeUm')} um)")
        c.ok(mol.get("n") and mol.get("n") == rmol.get("n") and mol.get("nGenes") == rmol.get("nGenes"),
             f"lone transcripts store: the same molecules as the Xenium folder without images ({mol.get('n')} / {rmol.get('n')}, "
             f"{mol.get('nGenes')} genes)")
    print(f"    {out['certain']} of {len(DETECT_WANT)} synthetic folders certain, listing {out['seconds']:.2f}s; "
          f"auto = explicit for {out['autoSame']} datasets")
    return out


# ------------------------------------------------------------------ batch export (batch.py, prepare.py batch)
# a folder of dataset folders: three demo folders copied as they are, one that only the options file settles,
# one that detection calls certain but prepare cannot read, one with an open choice, one whose name looks like
# a sample ID, and a stray file that is not a dataset
BATCH_GOOD = {"bx_xenium": "xenium_pq", "bx_geomx": "geomx", "bx_cosmx": "cosmx_tx"}
BATCH_MIXED = "bx_mixed"            # xenium_pq + MERSCOPE tables: two kinds; the options file says "kind": "xenium"
BATCH_BROKEN = "bx_broken"          # a Xenium folder whose stain file is not a TIFF: certain, then prepare fails
BATCH_AMBIG = "bx_two_runs"         # MERSCOPE with two segmentation runs: skipped, with the --cells-dir to add
BATCH_IDLIKE = "INT" + "9961_scan"  # a folder name that looks like a sample ID (made up; deid:ok): skipped, never written
BATCH_STD = "bx_sdata"              # sdata_xenium.zarr copied as bx_sdata.zarr: the id drops .zarr, pixel_size from the options
BATCH_STORE = "bx_store_transcripts"  # bx_store_transcripts.zarr.zip: a Xenium transcripts store on its own (molecules)
BATCH_OK = sorted(list(BATCH_GOOD) + [BATCH_MIXED, BATCH_STD, BATCH_STORE])


def make_batch_inputs(root: Path) -> tuple[Path, Path]:
    """root/batch/inputs (the dataset folders) and root/batch/options.json, from the synthetic raw/ folders."""
    b = root / "batch"
    if b.exists():
        shutil.rmtree(b)
    inp = b / "inputs"
    for name, raw in BATCH_GOOD.items():
        shutil.copytree(root / "raw" / raw, inp / name)
    shutil.copytree(root / "raw" / "xenium_pq", inp / BATCH_MIXED)
    shutil.copytree(root / "raw" / "sdata_xenium.zarr", inp / f"{BATCH_STD}.zarr")
    for n in ("detected_transcripts.csv", "cell_by_gene.csv", "cell_metadata.csv"):
        _touch(inp / BATCH_MIXED / n, b"a,b\n")
    for d in (BATCH_BROKEN, BATCH_IDLIKE):
        _touch(inp / d / "experiment.xenium", b"{}")
        _touch(inp / d / "morphology_focus" / "ch0000_dapi.ome.tif", b"not a TIFF")
        _touch(inp / d / "transcripts.parquet", b"not parquet")
    _touch(inp / BATCH_AMBIG / "detected_transcripts.csv", b"a,b\n")
    for run_ in ("vpt-run-a", "vpt-run-b"):
        for n in ("cell_by_gene.csv", "cell_metadata.csv"):
            _touch(inp / BATCH_AMBIG / run_ / n, b"a,b\n")
    _touch(inp / "notes.txt", b"not a dataset\n")
    shutil.copyfile(root / "raw" / "xenium" / "transcripts.zarr.zip", inp / f"{BATCH_STORE}.zarr.zip")
    shutil.copytree(root / "raw" / "xenium_pq", inp / "_parked_copy")   # a name starting with _: left out
    opts = b / "options.json"
    opts.write_text(json.dumps({"_comment": "the notebook's OPTIONS as JSON",
                                "bx_cosmx": {"pixel_size": 0.5, "label": "Batch CosMx"},
                                BATCH_MIXED: {"kind": "xenium"}, BATCH_STD: {"pixel_size": 0.5}}, indent=1), encoding="utf-8")
    return inp, opts


def batch_args(root: Path, *extra: str) -> list[str]:
    b = root / "batch"
    return ["web/SpatialViewer/prepare.py", "batch", "--src", str(b / "inputs"), "--out", str(b / "cache"),
            "--dest", str(b / "export"), "--options", str(b / "options.json"), *extra]


def run_batch(root: Path, log_path: Path) -> dict:
    """The batch three times: the first run (two at a time), the same again (resume), a dry run of one dataset
    with --rebuild. A failing dataset makes the exit code 1, so these do not go through run()."""
    make_batch_inputs(root)
    runs = {}
    for name, extra in (("first", ["--jobs", "2", "--markers-genes", "examples/demo_xenium_markers.json",
                                   "--title", "Batch test"]),
                        ("again", []),
                        ("dry", ["--only", "bx_cosmx", "--rebuild", "--dry-run", "--no-ids"])):
        args = batch_args(root, *extra)
        before = {p.name for p in (root / "batch" / "logs").glob("batch_*.json")} if (root / "batch" / "logs").exists() else set()
        print(f"$ python {' '.join(args)}", flush=True)
        t0 = time.time()
        p = subprocess.run([sys.executable, *args], cwd=str(REPO), stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                           text=True, encoding="utf-8", errors="replace", env={**os.environ, "PYTHONIOENCODING": "utf-8"})
        dt = time.time() - t0
        with open(log_path, "a", encoding="utf-8") as lg:
            lg.write(f"\n$ python {' '.join(args)}\n{p.stdout}")
        print("    " + "\n    ".join(p.stdout.rstrip().splitlines()[-8:]), flush=True)
        print(f"    (exit {p.returncode}, {dt:.1f}s)", flush=True)
        new = sorted({q.name for q in (root / "batch" / "logs").glob("batch_*.json")} - before) if (root / "batch" / "logs").exists() else []
        runs[name] = {"exit": p.returncode, "seconds": round(dt, 1), "stdout": p.stdout,
                      "report": str(root / "batch" / "logs" / new[-1]) if new else None, "newReports": len(new)}
    (root / "batch" / "runs.json").write_text(json.dumps(runs, indent=1), encoding="utf-8")
    return runs


def check_batch(c: Checker, root: Path, log_path: Path) -> dict:
    """prepare.py batch: the good datasets export, the broken one fails and is logged without stopping the others,
    the unsure one and the ID-like one are skipped with what to do, the exit code says so, datasets.json lists the
    good ones, the log is a table (md, json, csv) with no sample ID, and a second run skips what is done."""
    b = root / "batch"
    if not c.ok((b / "runs.json").exists(), "batch: the batch runs happened (batch/runs.json)"):
        return {}
    runs = read_json(b / "runs.json")
    first, again, dry = runs["first"], runs["again"], runs["dry"]
    c.eq(first["exit"], 1, "batch first run: exit code 1 (one dataset failed)")
    if not c.ok(first["report"] and Path(first["report"]).exists(), "batch first run: a JSON report in batch/logs"):
        return {}
    rep = read_json(Path(first["report"]))
    rows = {r["id"]: r for r in rep["datasets"]}
    idlike_shown = [k for k in rows if "[sample ID]" in k]
    c.eq(sorted(k for k in rows if k not in idlike_shown), sorted(BATCH_OK + [BATCH_BROKEN, BATCH_AMBIG]),
         "batch: one row per dataset folder, none for the stray notes.txt")
    c.ok(len(idlike_shown) == 1 and rows[idlike_shown[0]]["status"] == "skipped"
         and "sample ID" in " ".join(rows[idlike_shown[0]]["problems"]),
         f"batch: the folder whose name looks like a sample ID is skipped and shown as [sample ID] ({idlike_shown})")
    for ds in BATCH_OK:
        r = rows.get(ds, {})
        c.ok(r.get("status") == "ok" and (b / "export" / ds / "manifest.json").exists(),
             f"batch: {ds} ok and exported (got {r.get('status')}: {r.get('reason')})")
        c.ok((r.get("exportFiles") or 0) > 3 and (r.get("exportBytes") or 0) > 1000 and (r.get("seconds") or 0) > 0,
             f"batch: {ds} has its time and export size in the log ({r.get('seconds')} s, {r.get('exportFiles')} files)")
        if sys.platform == "win32" or sys.platform.startswith("linux"):
            c.ok((r.get("peakMemoryBytes") or 0) > 20 * 2**20,
                 f"batch: {ds} peak memory measured ({r.get('peakMemoryBytes')} bytes)")
    br = rows.get(BATCH_BROKEN, {})
    c.ok(br.get("status") == "failed" and "prepare failed" in br.get("reason", "") and len(br.get("errorTail") or []) >= 2,
         f"batch: {BATCH_BROKEN} failed in prepare, with the last lines of its error ({br.get('reason')})")
    c.ok(not (b / "export" / BATCH_BROKEN).exists() and not (b / "cache" / BATCH_BROKEN / "manifest.json").exists(),
         f"batch: nothing finished was written for {BATCH_BROKEN}")
    c.ok((b / "logs" / f"{BATCH_BROKEN}_prepare.log").exists(), f"batch: {BATCH_BROKEN}_prepare.log kept")
    am = rows.get(BATCH_AMBIG, {})
    hint = "\n".join(am.get("hint") or [])
    c.ok(am.get("status") == "skipped" and "not sure" in am.get("reason", "") and f'"{BATCH_AMBIG}": {{"kind": "merscope", "cells_dir": ' in hint
         and "--kind merscope --cells-dir" in hint and "vpt-run-a" in hint and "vpt-run-b" in hint,
         f"batch: {BATCH_AMBIG} skipped, with the options-file entry and the flag that settle it ({am.get('reason')}; {hint[:300]!r})")
    c.ok(not (b / "export" / BATCH_AMBIG).exists() and not (b / "cache" / BATCH_AMBIG).exists(),
         f"batch: nothing written for {BATCH_AMBIG}")
    st = read_json(b / "export" / BATCH_STORE / "manifest.json") if (b / "export" / BATCH_STORE / "manifest.json").exists() else {}
    c.ok(st.get("kind") == "xenium" and (st.get("image") or {}).get("background") == "transcript density" and st.get("molecules"),
         f"batch: {BATCH_STORE}.zarr.zip prepared on its own, molecules on a transcript-density background")
    c.ok("1 name(s) under the inputs folder start with _" in first["stdout"] and "_parked_copy" not in rows,
         "batch: a folder whose name starts with _ is left out, and the log says so")
    mixed = read_json(b / "export" / BATCH_MIXED / "manifest.json") if (b / "export" / BATCH_MIXED / "manifest.json").exists() else {}
    c.eq(mixed.get("kind"), "xenium", f"batch: {BATCH_MIXED} (two kinds) prepared as the options file's kind")
    cx = read_json(b / "export" / "bx_cosmx" / "manifest.json") if (b / "export" / "bx_cosmx" / "manifest.json").exists() else {}
    cxlog = (b / "logs" / "bx_cosmx_prepare.log").read_text(encoding="utf-8") if (b / "logs" / "bx_cosmx_prepare.log").exists() else ""
    c.ok(cx.get("label") == "Batch CosMx" and re.search(r"\[auto\] running: python prepare\.py cosmx .*--pixel-size 0\.5", cxlog),
         f"batch: the options file's label and pixel_size reach bx_cosmx's prepare ({cx.get('label')!r})")
    cat = read_json(b / "export" / "datasets.json") if (b / "export" / "datasets.json").exists() else {}
    c.eq(sorted(e["id"] for e in cat.get("datasets", [])), BATCH_OK, "batch: datasets.json lists the datasets that exported")
    c.eq(cat.get("title"), "Batch test", "batch: --title reaches datasets.json")
    c.ok(rep.get("deid", {}).get("clean") is True and rep.get("catalog", {}).get("exit") == 0,
         f"batch: the report has catalog.py's result and a clean sample ID check ({rep.get('deid', {}).get('text')})")
    stem = Path(first["report"]).with_suffix("")
    md = stem.with_suffix(".md").read_text(encoding="utf-8") if stem.with_suffix(".md").exists() else ""
    c.ok("| id | kind | status | time | peak memory | export size | note |" in md and f"| {BATCH_BROKEN} | xenium" in md
         and "## Failed" in md and "## Skipped" in md and "```" in md, "batch: the .md log has the table, the failure's last lines and the skipped ones")
    try:
        import csv as _csv
        with open(stem.with_suffix(".csv"), newline="", encoding="utf-8") as f:
            crow = {r["id"]: r for r in _csv.DictReader(f)}
        c.ok({k: v["status"] for k, v in crow.items()} == {k: v["status"] for k, v in rows.items()}
             and crow[BATCH_BROKEN]["error_tail"] and crow[BATCH_OK[0]]["peak_memory_mb"] != "",
             "batch: the .csv log has the same rows and statuses, the error tail and the peak memory")
    except (OSError, KeyError) as e:
        c.ok(False, f"batch: the .csv log could not be read ({e})")
    out = first["stdout"]
    c.ok("1 failed" in out and "[batch] log:" in out and str(b / "export") in out and "exit code 1" in out,
         "batch: the closing summary gives the counts, the paths and the exit code")
    # nothing in the export or the logs looks like a sample ID; the ID-like folder name is nowhere in them
    args = ["web/SpatialViewer/deid.py", "check", str(b / "export"), str(b / "logs")]
    print(f"$ python {' '.join(args)}", flush=True)
    p = subprocess.run([sys.executable, *args], cwd=str(REPO), stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                       text=True, encoding="utf-8", errors="replace")
    with open(log_path, "a", encoding="utf-8") as lg:
        lg.write(p.stdout)
    print("    " + (p.stdout.strip().splitlines() or ["(no output)"])[-1])
    c.eq(p.returncode, 0, "deid.py check on the batch export and its logs (exit code)")
    leak = [q for q in (b / "logs").rglob("*") if q.is_file() and BATCH_IDLIKE in q.read_text(encoding="utf-8", errors="ignore")]
    c.ok(not leak and BATCH_IDLIKE not in first["stdout"], f"batch: the ID-like folder name is in no log file and not printed ({leak[:2]})")
    # the second run: the finished ones are skipped, the broken one fails again on its own
    c.eq(again["exit"], 1, "batch second run: exit code 1 again (the broken dataset still fails)")
    rep2 = read_json(Path(again["report"])) if again.get("report") else {"datasets": []}
    rows2 = {r["id"]: r for r in rep2["datasets"]}
    c.ok(all(rows2.get(ds, {}).get("status") == "skipped" and rows2[ds]["reason"].startswith("already exported")
             for ds in BATCH_OK), "batch second run: every exported dataset skipped as already exported")
    c.eq(rows2.get(BATCH_BROKEN, {}).get("status"), "failed", f"batch second run: {BATCH_BROKEN} tried again and failed again")
    c.ok(not any("prepare" in s or "export" in s for ds in BATCH_OK for s in [rows2.get(ds, {}).get("steps", {})]),
         "batch second run: no prepare or export ran for the finished datasets")
    c.ok(again["seconds"] < first["seconds"], f"batch second run: faster ({again['seconds']}s vs {first['seconds']}s)")
    # the dry run
    c.ok(dry["exit"] == 0 and "would run" in dry["stdout"] and "prepare.py auto --src" in dry["stdout"]
         and dry["newReports"] == 0, f"batch --dry-run --only --rebuild: says what would run, writes no report (exit {dry['exit']})")
    ex = [ln for ln in dry["stdout"].splitlines() if "prepare.py export" in ln]
    c.ok(ex and all("--no-ids" in ln for ln in ex), f"batch --no-ids: passed on to prepare.py export ({ex[:1]})")
    withids = sorted(ds for ds in BATCH_OK if list((b / "export" / ds / "expr").rglob("ids.txt.gz")))
    c.ok(BATCH_STD in withids, f"batch: the exports hold their id lists by default (SpatialData cells; with lists: {withids})")
    print(f"    batch: {len(BATCH_OK)} exported, 1 failed, 2 skipped; first run {first['seconds']}s, again {again['seconds']}s")
    return {"exported": len(BATCH_OK), "firstSec": first["seconds"], "againSec": again["seconds"]}


def check_notebook_finder(c: Checker, root: Path):
    """The notebook's Section 2 cell, run on the batch inputs (as inputs/) and on detect/cases/zip_and_folder: the
    same datasets as batch.find_datasets, in the notebook's words, and the zip beside its folder noted."""
    import contextlib
    import io
    nb = json.loads((REPO / "notebooks" / "SODA_pipeline.ipynb").read_text(encoding="utf-8"))
    cell = next(("".join(x["source"]) for x in nb["cells"] if x["cell_type"] == "code" and "DATASETS = find_datasets()" in "".join(x["source"])), "")
    sys.path.insert(0, str(APP))
    try:
        import batch
    finally:
        sys.path.pop(0)
    for inputs, opts_file in ((root / "batch" / "inputs", root / "batch" / "options.json"),
                              (root / "detect" / "cases" / "zip_and_folder", None)):
        if not c.ok(cell and inputs.exists(), f"notebook Section 2 on {inputs.name}: the cell and the folder are there"):
            continue
        opts = batch.load_options(opts_file)[0] if opts_file else {}
        g = {"__name__": "nb_section2", "Path": Path, "sys": sys, "importlib": importlib, "SV": APP, "INPUTS": inputs,
             "VISIUM_BINS": "008um,016um"}
        g["_OPTS"] = opts                       # the options file stands in for the cell's OPTIONS
        buf = io.StringIO()
        try:
            with contextlib.redirect_stdout(buf):
                exec(compile(cell.replace("OPTIONS = {", "OPTIONS = _OPTS or {", 1), "notebook Section 2", "exec"), g)
        except Exception as e:  # noqa: BLE001
            c.ok(False, f"notebook Section 2 on {inputs.name}: runs ({type(e).__name__}: {e})")
            continue
        rows, _ = batch.find_datasets(inputs, opts, "008um,016um")
        nbrows = g.get("DATASETS") or []
        c.eq([(d["id"], d["report"]["verdict"], bool(d["problems"])) for d in nbrows],
             [(r["id"], r["verdict"], bool(r["problems"])) for r in rows],
             f"notebook Section 2 on {inputs.name}: the same datasets, verdicts and ready / not ready as batch.find_datasets")
        out = buf.getvalue()
        if inputs.name == "zip_and_folder":
            c.ok("Note:" in out and "may be the same data" in out and "run1_zip" in out and "READY" in out,
                 f"notebook Section 2: run1.zip beside run1/ is a second dataset, run1_zip, with a note ({out[-300:]!r})")
        else:
            c.ok("OPTIONS" in out and "the options file" not in out and BATCH_STORE in out,
                 "notebook Section 2 on the batch inputs: batch's messages in the notebook's words (OPTIONS), the lone store found")


def check_notebook(c: Checker):
    """notebooks/SODA_pipeline.ipynb: valid, no saved outputs, code that compiles, and its dataset finding is batch.py's."""
    p = REPO / "notebooks" / "SODA_pipeline.ipynb"
    try:
        nb = json.loads(p.read_text(encoding="utf-8"))
    except ValueError as e:
        c.ok(False, f"notebook: valid JSON ({e})")
        return
    code = [x for x in nb["cells"] if x["cell_type"] == "code"]
    c.ok(len(code) >= 10 and all(not x.get("outputs") and x.get("execution_count") is None for x in code),
         "notebook: no saved outputs or execution counts")
    bad = []
    for x in code:
        try:
            compile("".join(x["source"]), x.get("id", "cell"), "exec")
        except SyntaxError as e:
            bad.append(f"{x.get('id')}: {e}")
    c.ok(not bad, f"notebook: every code cell compiles {bad}")
    src = "\n".join("".join(x["source"]) for x in code)
    c.ok("import detect" in src and "import batch" in src and "batch.find_datasets(" in src and "OPTIONS = {" in src
         and "def check_dataset" not in src and "is_kind_folder" not in src,
         "notebook: Section 2 finds the datasets with batch.find_datasets (no copy of its own) and has the OPTIONS dictionary")
    c.ok("cell_by_gene" not in src and 'glob("*.qptiff' not in src and "morphology_focus" not in src and "_regist" not in src,
         "notebook: no per-kind file finding of its own (detect.py has it)")
    text = "\n".join("".join(x["source"]) for x in nb["cells"] if x["cell_type"] == "markdown")
    c.ok("own folder under `inputs/`" in text and "How sure" not in src.split("HOW_SURE", 1)[0] and "`kind`" in text,
         "notebook: Section 1 says one folder per dataset, Section 2 documents the kind override")


# ------------------------------------------------------------------ cell outlines (web/SpatialViewer/cellshapes.py)
# The fixtures above carry the cells' outlines in every form a platform writes them; two more exports of
# demo_xenium_nt take them from cells.zarr.zip and from "a lab's own" files (prepare.py outlines).
OL_ZARR, OL_OWN, OL_NONE, OL_NONE_CACHE = "export_zarr", "export_own", "export_none", "cache_none"


def outline_commands(root: Path) -> list[list[str]]:
    try:
        r = root.resolve().relative_to(REPO.resolve()).as_posix()
    except ValueError:
        r = root.resolve().as_posix()
    raw, cache, own = f"{r}/raw", f"{r}/cache", f"{r}/outline_files"
    prep, small = "web/SpatialViewer/prepare.py", ["--no-umap", "--no-hvg", "--grid-um", "0"]
    return [
        [prep, "outlines", "--id", XN, "--out", cache, "--src", f"{own}/xenium_zarr"],
        [prep, "export", "--id", XN, "--out", cache, "--dest", f"{r}/{OL_ZARR}", *small],
        [prep, "export", "--id", XN, "--out", cache, "--dest", f"{r}/{OL_OWN}", *small],
        [prep, "outlines", "--id", XN, "--out", cache, "--dest", f"{r}/{OL_OWN}", "--boundaries", f"{own}/xenium_own/cells.geojson",
         "--nuclei", f"{own}/xenium_own/nuclei.csv", "--boundary-units", "um"],
        # the same folder prepared with --no-outlines: the export must be what it was before outlines existed
        [prep, "xenium", "--id", XN, "--src", f"{raw}/xenium_nt", "--out", f"{r}/{OL_NONE_CACHE}", "--no-outlines"],
        [prep, "export", "--id", XN, "--out", f"{r}/{OL_NONE_CACHE}", "--dest", f"{r}/{OL_NONE}", *small],
    ]


def clean_outlines(root: Path):
    for p in [root / sub / XN for sub in (OL_ZARR, OL_OWN, OL_NONE, OL_NONE_CACHE)]:
        if p.exists():
            print(f"  removing earlier output {p}", flush=True)
            shutil.rmtree(p)


def _ol_area(r: np.ndarray) -> float:
    """Shoelace area; positive = clockwise on screen (x right, y down)."""
    x, y = r[:, 0] - r[0, 0], r[:, 1] - r[0, 1]
    return 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))


def _ol_inside(r: np.ndarray, x, y) -> np.ndarray:
    """Points inside a ring (even-odd)."""
    x, y = np.atleast_1d(np.asarray(x, float)), np.atleast_1d(np.asarray(y, float))
    x1, y1 = r[:, 0][:, None], r[:, 1][:, None]
    x2, y2 = np.roll(r[:, 0], -1)[:, None], np.roll(r[:, 1], -1)[:, None]
    with np.errstate(divide="ignore", invalid="ignore"):
        hit = ((y1 > y) != (y2 > y)) & (x < x1 + (y - y1) * (x2 - x1) / (y2 - y1))
    return hit.sum(0) % 2 == 1


def _ol_dist(p: np.ndarray, r: np.ndarray) -> np.ndarray:
    """Distance from each point of p to the outline of ring r."""
    a, b = r[None, :, :], np.roll(r, -1, axis=0)[None, :, :]
    q = p[:, None, :]
    ab = b - a
    L2 = (ab ** 2).sum(-1)
    t = np.clip(((q - a) * ab).sum(-1) / np.where(L2 > 0, L2, 1), 0, 1)
    return np.sqrt(((q - (a + t[..., None] * ab)) ** 2).sum(-1)).min(1)


def _ol_same(stored: np.ndarray, truth: np.ndarray, tol: float = 0.0) -> bool:
    """A stored ring against the ring it was made from: every stored vertex is one of the truth's
    (to the vertex step), and no truth vertex lies farther than `tol` px from the stored outline;
    tol = 0: the very same vertices."""
    d = np.hypot(stored[:, None, 0] - truth[None, :, 0], stored[:, None, 1] - truth[None, :, 1]).min(1)
    if d.max() > 0.05 or (tol == 0 and len(stored) != len(truth)):
        return False
    return bool(_ol_dist(truth, stored).max() <= tol + 0.05)


def _ol_read(c: Checker, out: Path, what: str, key: str = "cells") -> dict:
    """Every outline set of a layer, decoded and checked against the format (docs/DATA_FORMAT.md):
    {set: {"index", "rings": {row: [ring, ...]}, "empty": rows listed without an outline, "pos"}}."""
    m = read_json(out / "manifest.json")
    L = next((x for x in m.get("expr", []) if x.get("key") == key), {})
    blk = L.get("shapes")
    if not c.ok(isinstance(blk, dict) and blk.get("sets"), f"{what}: the {key} layer has no shapes block in the manifest"):
        return {}
    c.eq((blk.get("version"), blk.get("format"), blk.get("dir")), (1, "poly-u16-tiles", "shapes"), f"{what}: shapes version / format / dir")
    pos = f32(out / "expr" / key / "positions.bin", 2)
    N = len(pos)
    if L.get("hvg"):
        c.eq(blk.get("chunkPx"), L["hvg"]["chunkPx"], f"{what}: outlines use the same tile grid as hvg/")
    res = {}
    for si, s in enumerate(blk["sets"]):
        tag = f"{what} {s.get('key')}"
        d = out / "expr" / key / "shapes" / s["key"]
        if not c.ok((d / "index.json").exists(), f"{tag}: index.json missing"):
            continue
        idx = read_json(d / "index.json")
        cp, step = idx["chunkPx"], idx["step"]
        c.eq((idx.get("format"), idx.get("version"), idx.get("set"), idx.get("nRows"), cp),
             ("poly-u16-tiles", 1, s["key"], N, blk.get("chunkPx")), f"{tag}: index format / version / set / nRows / chunkPx")
        c.ok(step >= 1 / 16 and math.log2(step) == int(math.log2(step)), f"{tag}: vertex step {step} is a power of two, 1/16 px at the finest")
        c.eq(sorted(p.name for p in d.glob("*.bin")), sorted(f"{k['tx']}_{k['ty']}.bin" for k in idx["chunks"]),
             f"{tag}: the tile files are those index.json lists")
        rings, empty, bad, nP, nV, total, over = {}, set(), [], 0, 0, 0, 0.0
        for k in idx["chunks"]:
            buf = (d / f"{k['tx']}_{k['ty']}.bin").read_bytes()
            n, nv = (int(v) for v in np.frombuffer(buf, "<u4", 2))
            if len(buf) != 16 + 6 * n + 4 * nv or (n, nv, len(buf)) != (k["n"], k["nv"], k["bytes"]):
                bad.append(f"{k['tx']}_{k['ty']}: {len(buf)} bytes for n {n}, nv {nv} (index: {k})")
                continue
            origin = np.frombuffer(buf, "<f4", 2, 8).astype(np.float64)
            rows = np.frombuffer(buf, "<u4", n, 16).astype(np.int64)
            cnt = np.frombuffer(buf, "<u2", n, 16 + 4 * n).astype(np.int64)
            xy = origin + np.frombuffer(buf, "<u2", 2 * nv, 16 + 6 * n).reshape(-1, 2).astype(np.float64) * step
            if int(cnt.sum()) != nv or (n and int(rows.max()) >= N):
                bad.append(f"{k['tx']}_{k['ty']}: vertex counts or rows out of range")
                continue
            t = np.floor(pos[rows].astype(np.float64) / cp).astype(np.int64)
            if not ((t[:, 0] == k["tx"]) & (t[:, 1] == k["ty"])).all():
                bad.append(f"{k['tx']}_{k['ty']}: a polygon is not in the tile of its cell's position")
            if nv:
                over = max(over, k["tx"] * cp - xy[:, 0].min(), k["ty"] * cp - xy[:, 1].min(),
                           xy[:, 0].max() - (k["tx"] + 1) * cp, xy[:, 1].max() - (k["ty"] + 1) * cp)
            o = 0
            for r, k_ in zip(rows.tolist(), cnt.tolist()):
                if k_:
                    rings.setdefault(r, []).append(xy[o:o + k_])
                else:
                    empty.add(r)
                o += k_
            nP, nV, total = nP + int((cnt > 0).sum()), nV + nv, total + len(buf)
        c.ok(not bad, f"{tag}: tile files: {bad[:3]}")
        c.eq((nP, nV, total, len(rings)), (idx.get("nPolygons"), idx.get("nVertices"), idx.get("bytes"), idx.get("nRowsWithShape")),
             f"{tag}: polygons / vertices / bytes / cells with an outline as index.json says")
        c.eq((s.get("nPolygons"), s.get("nVertices"), s.get("bytes"), s.get("chunks"), s.get("nRowsWithShape")),
             (nP, nV, total, len(idx["chunks"]), len(rings)), f"{tag}: the manifest's numbers")
        c.ok(over <= idx.get("overhangPx", -1) + 1e-6, f"{tag}: a vertex reaches {over:.2f} px outside its tile, index.json says {idx.get('overhangPx')}")
        if si == 0:
            c.ok(idx.get("complete") is True and (set(rings) | empty) == set(range(N)) and not (set(rings) & empty),
                 f"{tag}: the first set lists every row (with or without an outline)")
        else:
            c.ok(not empty and idx.get("complete") is False, f"{tag}: a further set lists only cells with an outline")
        flat = [r for rs in rings.values() for r in rs]
        c.ok(all(3 <= len(r) <= 64 for r in flat), f"{tag}: every ring has 3 to 64 vertices")
        c.ok(all(_ol_area(r) > 0 for r in flat), f"{tag}: every ring is clockwise on screen")
        c.ok(all((r[0] != r[-1]).any() for r in flat), f"{tag}: rings are stored open (the first vertex is not repeated)")
        txt = (d / "index.json").read_text(encoding="utf-8")
        c.ok(not BARCODE.search(txt) and not LOCAL_PATH.search(txt), f"{tag}: index.json holds a barcode or a local path")
        res[s["key"]] = {"index": idx, "rings": rings, "empty": empty, "pos": pos}
    return res


def _ol_centres(c: Checker, S: dict, what: str, want: float = 1.0):
    """Each cell's own centre lies inside (one of) its own outline(s)."""
    pos, rings = S["pos"], S["rings"]
    ok = sum(any(_ol_inside(r, pos[i, 0], pos[i, 1])[0] for r in rs) for i, rs in rings.items())
    c.ok(ok >= want * len(rings), f"{what}: {ok} of {len(rings)} cell centres lie inside their own outline")
    c.ok(abs(S["index"]["stats"]["centreInside"] - ok / max(1, len(rings))) < 1e-4, f"{what}: index.json centreInside is the share counted here")


def _ol_xenium(c: Checker, out: Path, what: str, tr: dict, tz, own: bool = False, not_in_matrix: int = 0) -> dict:
    """The Xenium fixture's outlines in an export against the generator's truth."""
    px, sp = tr["pixelSizeUm"], tr["outlines"]
    R = _ol_read(c, out, what)
    if "cell" not in R or "nucleus" not in R:
        c.ok(False, f"{what}: expected a cell and a nucleus set, got {sorted(R)}")
        return {}
    C, Nu = R["cell"], R["nucleus"]
    n = len(C["pos"])
    cell, nuc, dense = tz["cell_poly_um"] / px, tz["nucleus_poly_um"] / px, tz["dense_poly_um"] / px
    none = set() if own else set(sp["noBoundary"])
    c.eq(C["empty"], none, f"{what}: the cells without a boundary are listed without an outline (they stay dots)")
    multi = {sp["ownMulti"]} if own else set()
    bad = [i for i in range(n) if i not in none and i not in multi and (own or i != sp["dense"])
           and not (len(C["rings"].get(i, [])) == 1 and _ol_same(C["rings"][i][0], cell[i]))]
    c.ok(not bad, f"{what}: cell outlines differ from the fixture's polygons (matrix order; the file is in another order): rows {bad[:5]}")
    st = C["index"]["stats"]
    if own:
        rs = C["rings"].get(sp["ownMulti"], [])
        c.ok(len(rs) == 2 and sorted(len(r) for r in rs) == [8, 12], f"{what}: a cell in two parts keeps both ({[len(r) for r in rs]})")
        c.eq((st["rowsInSeveralParts"], st["dropped"]["holes"], st["dropped"]["smallParts"]), (1, 1, 1),
             f"{what}: one cell in parts, one hole left out, one crumb (1 % of its cell) dropped")
        c.ok(len(C["rings"].get(sp["ownCrumb"], [])) == 1 and len(C["rings"].get(sp["ownHole"], [])) == 1, f"{what}: the crumb and the hole are gone")
    else:
        d = C["rings"].get(sp["dense"], [np.zeros((0, 2))])[0]
        c.ok(8 <= len(d) < sp["denseN"] and _ol_same(d, dense, tol=0.1 / px),
             f"{what}: the {sp['denseN']}-vertex outline is simplified ({len(d)} left) and no vertex moved by more than 0.1 um")
        c.eq(st["polygonsNotInMatrix"], not_in_matrix, f"{what}: polygons of cells that are not in the matrix are counted and left out")
    _ol_centres(c, C, f"{what} cell")
    area = tz["cell_area_um2"]
    ratio = np.array([_ol_area(C["rings"][i][0]) * px * px / area[i] for i in range(n) if i in C["rings"] and i not in multi
                      and (own or i != sp["dense"])])
    c.ok(0.94 < float(np.median(ratio)) < 0.97 and ratio.min() > 0.93 and ratio.max() < 0.98,
         f"{what}: polygon area / cell_area {np.median(ratio):.3f} ({ratio.min():.3f}-{ratio.max():.3f}); a 12-gon holds 95.5 % of its circle")
    twins = {} if own else {i: tz["twin_poly_um"][k] / px for k, i in enumerate(sp["binucleate"])}
    no_nuc = set() if own else set(sp["noNucleus"])
    c.eq(set(range(n)) - set(Nu["rings"]), no_nuc, f"{what}: cells without a nucleus have no entry in the nucleus set")
    bad = [i for i in range(n) if i not in no_nuc and i not in twins
           and not (len(Nu["rings"].get(i, [])) == 1 and _ol_same(Nu["rings"][i][0], nuc[i]))]
    c.ok(not bad, f"{what}: nucleus outlines differ from the fixture's: rows {bad[:5]}")
    for i, tw in twins.items():
        rs = Nu["rings"].get(i, [])
        c.ok(len(rs) == 2 and all(any(_ol_same(r, t) for t in tw) for r in rs), f"{what}: a cell with two nuclei keeps both (row {i})")
    return {"cells": len(C["rings"]), "nuclei": sum(len(v) for v in Nu["rings"].values()), "bytesPerCell": round(C["index"]["bytes"] / n, 1)}


def check_outlines(c: Checker, root: Path, log_path: Path) -> dict:
    """Cell outlines of every fixture (cellshapes.py): the tile format, each outline against the
    polygon or label mask it came from, rows, the special cases, and what is refused."""
    raw, export = root / "raw", root / "export"
    res = {}
    tr, tz = read_json(raw / "truth_xenium.json"), np.load(raw / "truth_xenium.npz")
    px = tr["pixelSizeUm"]
    print("  xenium: cell_boundaries / nucleus_boundaries .csv.gz, .parquet, cells.zarr.zip, a lab's own files")
    res["csv"] = _ol_xenium(c, export / XE, f"{XE} (csv.gz)", tr, tz, not_in_matrix=tr["outlines"]["notInMatrix"])
    _ol_xenium(c, export / XN, f"{XN} (parquet)", tr, tz, not_in_matrix=tr["outlines"]["notInMatrix"])
    _ol_xenium(c, root / OL_ZARR / XN, f"{XN} (cells.zarr.zip)", tr, tz)
    _ol_xenium(c, root / OL_OWN / XN, f"{XN} (own GeoJSON + vertex table)", tr, tz, own=True)
    for ds in (XE, XN):
        ix = read_json(export / ds / "expr" / "cells" / "shapes" / "cell" / "index.json")
        c.eq(ix.get("source"), "cell_boundaries.csv.gz" if ds == XE else "cell_boundaries.parquet", f"{ds}: outlines source")
    c.eq(read_json(root / OL_ZARR / XN / "expr" / "cells" / "shapes" / "cell" / "index.json").get("source"), "cells.zarr.zip", "zarr export: source")
    try:
        import gzip
        ids = gzip.open(raw / "xenium" / "cell_feature_matrix" / "barcodes.tsv.gz", "rt").read().split()
        blob = b"".join(p.read_bytes() for p in (export / XE / "expr" / "cells" / "shapes").rglob("*") if p.is_file())
        c.ok(not any(i.encode() in blob for i in ids), f"{XE}: a cell id is written in the outline files")
    except OSError as e:
        c.ok(False, f"{XE}: cannot read the outline files ({e})")
    c.ok("shapes" not in json.dumps(read_json(export / XP / "manifest.json")), f"{XP}: no cells, so no outlines")
    try:                                               # the older cells.zarr.zip layout, read in-process
        sys.path.insert(0, str(APP))
        import cellshapes
        got = cellshapes.read_xenium_zarr(root / "outline_files" / "xenium_zarr_v1" / "cells.zarr.zip")
        ids1, st1, xy1 = got["cell"]
        want_ids = [b for i, b in enumerate(ids) if i not in tr["outlines"]["noBoundary"]]
        keep = np.flatnonzero(np.diff(st1) > 0)
        c.eq([ids1[i] for i in keep], want_ids, "cells.zarr.zip, older layout: cell ids from (prefix, suffix), cells without a polygon skipped")
        k0 = int(keep[0])
        r0 = xy1[st1[k0]:st1[k0 + 1]]
        c.ok(len(r0) == 13 and np.allclose(r0[:12], tz["cell_poly_um"][0], atol=1e-3), "cells.zarr.zip, older layout: polygon_vertices read as x, y pairs")
        c.eq(len(got["nucleus"][0]), len(ids), "cells.zarr.zip, older layout: one nucleus slot per cell")
        t = cellshapes.trace_labels(np.array([[0, 1, 1, 0], [0, 1, 0, 2], [3, 0, 2, 2]], np.uint8))
        c.eq((t[0].tolist(), np.diff(t[1]).tolist(), t[3]), ([1, 2, 3], [6, 6, 4], 0), "trace_labels: three labels, corners only, no holes")
    except ImportError as e:
        print(f"    SKIP in-process reader checks ({e})")

    # ---- CosMx: traced from the CellLabels images / the polygons file
    ti, tzi = read_json(raw / "truth_imaging.json"), np.load(raw / "truth_imaging.npz")
    order = np.concatenate([np.flatnonzero(tzi["cosmx_fov"] == f) for f in range(1, 5)])
    print("  cosmx: traced CellLabels, polygons file")
    R = _ol_read(c, export / CX, CX)
    if "cell" in R:
        S, lab = R["cell"], tzi["cosmx_labels"]
        fw, fh = ti["cosmx"]["fovSize"]
        gx, gy = ti["cosmx"]["width"] - 2 * fw, ti["cosmx"]["height"] - 2 * fh
        c.eq(S["index"].get("source"), "CellLabels images (traced)", f"{CX}: outlines source")
        c.ok(not S["empty"] and S["index"]["stats"].get("traced") is True, f"{CX}: every cell has a traced outline")
        iou = []
        for row, ti_ in enumerate(order.tolist()):
            f, cid = int(tzi["cosmx_fov"][ti_]), int(tzi["cosmx_cell_id"][ti_])
            x0, y0 = ((f - 1) % 2) * (fw + gx), (1 - (f - 1) // 2) * (fh + gy)
            mask = lab[y0:y0 + fh, x0:x0 + fw] == cid
            ys, xs = np.nonzero(mask)
            a0, a1, b0, b1 = xs.min() - 2, xs.max() + 3, ys.min() - 2, ys.max() + 3
            qx, qy = np.meshgrid(np.arange(a0, a1) + 0.5 + x0, np.arange(b0, b1) + 0.5 + y0)
            got = np.zeros(qx.shape, bool)
            for r in S["rings"].get(row, []):
                got |= _ol_inside(r, qx.ravel(), qy.ravel()).reshape(qx.shape)
            want = np.zeros(qx.shape, bool)
            want[ys - b0, xs - a0] = True
            iou.append((got & want).sum() / max(1, (got | want).sum()))
        iou = np.array(iou)
        c.ok(float(np.median(iou)) > 0.93 and float(iou.min()) > 0.85,
             f"{CX}: traced outlines against the label masks: IoU median {np.median(iou):.3f}, worst {iou.min():.3f}")
        _ol_centres(c, S, CX)
        mv = float(np.mean([len(r) for rs in S["rings"].values() for r in rs]))
        c.ok(6 <= mv <= 24, f"{CX}: {mv:.1f} vertices a traced outline (a 28 px disk; the pixel staircase is simplified away)")
        res["cosmxTraced"] = {"iouMedian": round(float(np.median(iou)), 4), "iouMin": round(float(iou.min()), 4), "meanVertices": round(mv, 1)}
        print(f"    traced outlines: IoU with the label mask median {np.median(iou):.3f}, worst {iou.min():.3f}; {mv:.1f} vertices each")
    R = _ol_read(c, export / CXN, CXN)
    if "cell" in R:
        S, want = R["cell"], tzi["cosmx_poly_px"][order]
        c.eq(S["index"].get("source"), "polygons file", f"{CXN}: the polygons file is used when the export has one")
        bad = [i for i in range(len(want)) if not (len(S["rings"].get(i, [])) == 1 and _ol_same(S["rings"][i][0], want[i]))]
        c.ok(not bad, f"{CXN}: outlines differ from the polygons file put through the frame (y turned over): rows {bad[:5]}")
        c.ok(all(_ol_area(r) > 0 for rs in S["rings"].values() for r in rs), f"{CXN}: rings turned clockwise after the y flip")
    c.ok("shapes" not in json.dumps(read_json(export / CXT / "manifest.json")), f"{CXT}: no cells, so no outlines")

    # ---- MERSCOPE: cell_boundaries/*.hdf5 (with images) / *_micron_space.parquet (WKB, no images)
    print("  merscope: cell_boundaries hdf5, micron_space parquet (WKB)")
    t = ti["merscope"]
    org = np.array(t["origin"], np.float64)
    pu = tzi["merscope_poly_um"]
    cen = tzi["merscope_cells_um"]
    um_all = np.vstack([tzi["merscope_mol_um"], cen])
    for ds, to_px, src in ((MS, lambda u: (u - org) / ti["pixelSizeUm"], "cell_boundaries/*.hdf5"),
                           (MSN, lambda u: u - um_all.min(0) + 10.0, "micron_space.parquet")):
        R = _ol_read(c, export / ds, ds)
        if "cell" not in R:
            continue
        S = R["cell"]
        st = S["index"]["stats"]
        c.eq(S["index"].get("source"), src, f"{ds}: outlines source")
        c.eq(S["empty"], {4}, f"{ds}: the cell without a boundary stays a dot")
        c.eq((st.get("zPlane"), st.get("otherPlane")), (1, 1), f"{ds}: the middle z-plane, and the nearest one for the cell that has none there")
        other = to_px(cen[3] + 0.8 * (pu[3] - cen[3]))
        c.ok(len(S["rings"].get(3, [])) == 1 and _ol_same(S["rings"][3][0], other), f"{ds}: that cell's outline is the one of its nearest plane")
        skip = {3, 4} | ({0} if ds == MSN else set())
        bad = [i for i in range(len(pu)) if i not in skip and not (len(S["rings"].get(i, [])) == 1 and _ol_same(S["rings"][i][0], to_px(pu[i])))]
        c.ok(not bad, f"{ds}: outlines differ from the fixture's (the middle plane, through the micron -> pixel transform): rows {bad[:5]}")
        if ds == MSN:
            rs = S["rings"].get(0, [])
            c.ok(len(rs) == 2 and any(_ol_same(r, to_px(pu[0])) for r in rs), f"{ds}: a MultiPolygon cell keeps both parts")
            c.eq((st["rowsInSeveralParts"], st["dropped"]["holes"], st["dropped"]["smallParts"]), (1, 1, 1),
                 f"{ds}: one cell in parts, one hole left out, one crumb dropped")
        _ol_centres(c, S, ds)
    c.ok("shapes" not in json.dumps(read_json(export / MSP / "manifest.json")), f"{MSP}: no cells, so no outlines")

    # ---- Stereo-seq: cellBin/cellBorder, whole chip and cut (file offset)
    print("  stereoseq: cellBin/cellBorder")
    tzs = np.load(raw / "truth_stereoseq.npz")
    R = _ol_read(c, export / ST, ST)
    if "cell" in R:
        S, cxy = R["cell"], tzs["cell_xy"].astype(np.float64)
        ang = np.linspace(0, 2 * np.pi, 16, endpoint=False)
        rr = read_json(raw / "truth_stereoseq.json")["cellRadiusDnb"]
        off = np.column_stack([np.rint(rr * np.cos(ang)), np.rint(rr * np.sin(ang))])
        bad = [i for i in range(len(cxy)) if not (len(S["rings"].get(i, [])) == 1 and _ol_same(S["rings"][i][0], cxy[i] + 0.5 + off, tol=0.2))]
        c.ok(not bad and not S["empty"], f"{ST}: outlines = cell x, y + the border offsets + half a DNB: rows {bad[:5]}")
        c.eq(S["index"].get("source"), "cellBin/cellBorder", f"{ST}: outlines source")
        _ol_centres(c, S, ST)
        ratio = np.array([_ol_area(S["rings"][i][0]) for i in range(len(cxy))]) / tzs["cell_area"]
        c.ok(0.8 < float(np.median(ratio)) < 1.0, f"{ST}: outline area / the cell's DNB count {np.median(ratio):.3f} (the border runs through the outermost DNBs)")
        a, b = export / ST / "expr" / "cells" / "shapes" / "cell", export / STL / "expr" / "cells" / "shapes" / "cell"
        c.ok(b.is_dir() and all((b / p.name).exists() and (b / p.name).read_bytes() == p.read_bytes() for p in a.glob("*.bin")),
             f"{STL}: a cut cell bin GEF (offsetX) gives the same outline tiles")
    c.ok("shapes" not in json.dumps(read_json(export / STG / "manifest.json")), f"{STG}: no cell bin GEF, so no outlines")

    # ---- what is refused (and then changes nothing); a dataset without outlines exports as it did before
    print("  refusals: wrong units, another run's ids, a folder without boundaries")
    def ol(*extra):
        args = ["web/SpatialViewer/prepare.py", "outlines", "--id", XN, "--out", str(root / "cache"), *extra]
        print(f"$ python {' '.join(args)}", flush=True)
        p = subprocess.run([sys.executable, *args], cwd=str(REPO), stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                           text=True, encoding="utf-8", errors="replace")
        with open(log_path, "a", encoding="utf-8") as lg:
            lg.write(p.stdout)
        return p
    def cached():
        return next((L for L in read_json(root / "cache" / XN / "manifest.json")["expr"] if L["key"] == "cells"), {}).get("shapes")
    def files():
        d = root / "cache" / XN / "shapes" / "cells"
        return {p.relative_to(d).as_posix(): p.read_bytes() for p in sorted(d.rglob("*")) if p.is_file()}
    gj = str(root / "outline_files" / "xenium_own" / "cells.geojson")
    had, held = cached(), files()
    c.ok(bool(had) and had["sets"][0]["source"] == "GeoJSON polygons (own file)" and len(held) > 4, "the cache holds the own-file outlines before the refusals")
    p = ol("--boundaries", gj, "--boundary-units", "px")
    c.ok(p.returncode == 0 and "do not sit on the cells" in p.stdout and "nothing written" in p.stdout,
         "outlines in the wrong units (um read as pixels) are refused: the centres are not inside them")
    p = ol("--boundaries", gj, "--boundary-units", "um", "--boundary-id", "other_run")
    c.ok(p.returncode == 0 and "not this dataset's outlines" in p.stdout and "nothing written" in p.stdout,
         "outlines whose ids are not this dataset's cells are refused")
    p = ol("--src", str(raw / "xenium_pq"))
    c.ok(p.returncode == 0 and "no cell outlines" in p.stdout and "nothing written" in p.stdout, "a folder without boundary files: a message, nothing written")
    c.ok(cached() == had and files() == held, "a refused import changes nothing: the outlines the dataset already had are kept, byte for byte")
    none = root / OL_NONE / XN
    c.ok((none / "manifest.json").exists() and "shapes" not in (none / "manifest.json").read_text(encoding="utf-8")
         and not (none / "expr" / "cells" / "shapes").exists() and not (root / OL_NONE_CACHE / XN / "shapes").exists(),
         "--no-outlines: no shapes in the cache, in the manifest or in the export")
    same = [f for f in ("positions.bin", "genes.json") if (none / "expr" / "cells" / f).read_bytes() != (export / XN / "expr" / "cells" / f).read_bytes()]
    c.ok(not same, f"--no-outlines: the cells layer itself is the same as with outlines ({same})")
    args = ["web/SpatialViewer/deid.py", "check"] + [str(p_) for p_ in (root / OL_ZARR / XN, root / OL_OWN / XN, none)] + \
           [str(export / ds / "expr" / "cells" / "shapes") for ds in (CX, CXN, MS, MSN, ST, STL)]
    p = subprocess.run([sys.executable, *args], cwd=str(REPO), stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                       text=True, encoding="utf-8", errors="replace")
    c.eq(p.returncode, 0, "deid.py check on the outline exports (exit code)")

    # ---- the page
    html = (APP / "index.html").read_text(encoding="utf-8")
    js = (APP / "outlines.js").read_text(encoding="utf-8") if (APP / "outlines.js").exists() else ""
    c.ok(re.search(r'<script src="outlines\.js\?v=[\w.-]+"></script>', html) is not None and html.find("outlines.js") < html.find("viewer.js?"),
         "index.html loads outlines.js (with a ?v= string) before viewer.js")
    for i in ("rowCellShapes", "cellShapes", "rowCellNuclei", "cellNuclei"):
        c.ok(f'id="{i}"' in html, f"index.html: the cell shapes control #{i} is missing")
    c.ok("window.SodaOutlines" in js and "poly-u16-tiles" in js, "outlines.js: SodaOutlines and the format name")
    c.ok(not re.search(r"\bInt(8|16|32)Array\b", js), "outlines.js uses a signed typed array (deid.py check reads those names as sample IDs)")
    vj = (APP / "viewer.js").read_text(encoding="utf-8")
    c.ok("SodaOutlines" in vj, "viewer.js: the outlines hooks are missing")
    # what draws a subset of the cells hands its colours to the outlines (a colour update, never more polygons)
    for needle, f, what in (("hl = lassoHl()", vj, "viewer.js: Draw a region's colours go to the outlines"),
                            ("regionOnlyRows(", vj, "viewer.js: the region click-filter cuts the outlines' gene colouring"),
                            ("outlines.legend()", vj, "viewer.js: Save image's legend names the cell outlines"),
                            ('key: "cs"', js, "outlines.js: a link to the view carries Cell shapes")):
        c.ok(needle in f, f"{what} ({needle!r} not found)")
    # in the browser: tests/outlines_browser_test.js (SKIPPED without Playwright; WebKit / Firefox SKIPPED when not installed)
    node = shutil.which("node")
    if not node:
        print("    node not found: tests/outlines_browser_test.js skipped")
        return res
    p = subprocess.run([node, str(REPO / "tests" / "outlines_browser_test.js"), "--export", str(export),
                        "--shots", str(root / "outlines_shots"), "--browsers", "chromium,webkit,firefox"], cwd=str(REPO),
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace")
    with open(log_path, "a", encoding="utf-8") as lg:
        lg.write(p.stdout)
    print_skipped(p.stdout)
    print("    " + (p.stdout.strip().splitlines() or ["(no output)"])[-1])
    c.eq(p.returncode, 0, "node tests/outlines_browser_test.js: outlines, and what draws a subset of the cells with them, in Chromium, WebKit and Firefox (exit code)")
    # tests/controls_browser_test.js: one control for each job, old links and Save JSON with the image off
    p = subprocess.run([node, str(REPO / "tests" / "controls_browser_test.js"), "--export", str(export),
                        "--shots", str(root / "controls_shots"), "--browsers", "chromium,webkit,firefox"], cwd=str(REPO),
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace")
    with open(log_path, "a", encoding="utf-8") as lg:
        lg.write(p.stdout)
    print_skipped(p.stdout)
    print("    " + (p.stdout.strip().splitlines() or ["(no output)"])[-1])
    c.eq(p.returncode, 0, "node tests/controls_browser_test.js: one control for each job, old links with the image off, "
                          "in Chromium, WebKit and Firefox (exit code)")
    return res


def check_marker_lists(c: Checker):
    """The built-in lists load and are usable, and the channel-name cleanup is the same in Python and JS."""
    sys.path.insert(0, str(APP))
    try:
        import celltypes
        import seams  # noqa: F401  (used below)
    finally:
        sys.path.pop(0)
    for name, prot in (("protein_generic.json", True), ("rna_generic.json", False), ("akoya_nerve.json", True),
                       ("nerve.json", False), ("brain.json", False)):
        ms = read_json(APP / "markers" / name)
        ok = (isinstance(ms.get("types"), list) and ms["types"]
              and all(t.get("name") and t.get("genes") for t in ms["types"])
              and (ms.get("modality") == "protein") == prot)
        c.ok(ok, f"markers/{name}: types / modality")
    pg = read_json(APP / "markers" / "protein_generic.json")
    c.eq(pg.get("minChannels"), 2, "protein_generic.json minChannels")
    c.eq(sorted(t["name"] for t in pg["types"] if t.get("min") == 1), ["Immune rich", "Lymphatic", "Proliferating"],
         "protein_generic.json one-marker types")
    al = pg.get("aliases", {})
    c.ok(all(k in {g for t in pg["types"] for g in t["genes"]} for k in al), "protein_generic.json: alias of a name not in any type")
    cases = {"GFAP-488": "gfap", "CD11c FFPE": "cd11c", "01_CD45": "cd45", "CD45 (Opal520)": "cd45", "CD3_AF647": "cd3",
             "NF-L": "nfl", "PD-1": "pd1", "Ki-67": "ki67", "CD138": "cd138", "CD8 Cy5.5": "cd8", "12 - CD4": "cd4",
             "Pan-CytokeratinGFAP": "pancytokeratingfap", "HLA-DR": "hladr"}
    bad = {k: celltypes.norm_name(k) for k, v in cases.items() if celltypes.norm_name(k) != v}
    c.ok(not bad, f"celltypes.norm_name: {bad}")
    m = celltypes.match_channels([{"name": "T", "genes": ["CD3", "CD4"]}], ["01_CD3e", "CD4-488"], {"CD3": ["CD3e"]})
    c.eq(m[0]["cols"], [0, 1], "celltypes.match_channels with aliases and tags")
    # the JS copy (lasso_worker.js normName) must use the same tag list
    js = (APP / "lasso_worker.js").read_text(encoding="utf-8")
    py = (APP / "celltypes.py").read_text(encoding="utf-8")
    tag = lambda txt: re.search(r"fitc\|tritc\|pe\|apc\|ffpe\|([0-9|]+)\)", txt)
    tj, tp = tag(js), tag(py)
    c.ok(tj and tp and tj.group(1) == tp.group(1), "lasso_worker.js normName tags differ from celltypes.norm_name")
    print(f"    marker lists ok; norm_name {len(cases) - len(bad)}/{len(cases)} cases")

    # the seam check on profiles with a known answer: camera tiles 480 px wide, 360 px high
    import seams
    rng = np.random.default_rng(0)
    W, H = 9600, 7200
    smooth = lambda n: np.convolve(rng.normal(0, 1, n + 400), np.ones(200) / 200, mode="valid")[:n] * 20 + 60
    tile = lambda n, sp: (1 + 0.15 * rng.normal(0, 1, n // sp + 1))[np.arange(n) // sp]   # each tile its own gain
    prof = lambda n, sp: smooth(n) * (tile(n, sp) if sp else 1) + rng.normal(0, 0.3, n)
    six = lambda cw, rw: (cw, rw, cw, cw, rw, rw)                 # whole slide and halves alike
    plain = [six(prof(W, 0), prof(H, 0)) for _ in range(4)]
    tiled = [six(prof(W, 480), prof(H, 360)) for _ in range(3)] + plain[:1]
    nm = ["A", "B", "C", "D"]
    sp, sn = seams.find_grid(tiled, nm, 0.5), seams.find_grid(plain, nm, 0.5)
    c.ok(sp["flagged"] and sp["columns"].get("spacingPx") in (476, 480, 484) and sp["rows"].get("spacingPx") in (356, 360, 364),
         f"protgrid.seam_check missed synthetic 480 x 360 px tiles: {sp}")
    c.ok(not sn["flagged"], f"protgrid.seam_check flagged a slide without tiles: {sn}")
    print(f"    seam check: tiles found at {sp['columns'].get('spacingPx')} x {sp['rows'].get('spacingPx')} px; "
          f"none on the plain profiles")


def check_viewer_page(c: Checker, root: Path, log_path: Path):
    """The page's files hang together, and the tile path keeps its fallbacks and its error reports."""
    html = (APP / "index.html").read_text(encoding="utf-8")
    refs = re.findall(r'<(?:script|link)[^>]+?(?:src|href)="([^"]+)"', html)
    # vendor/ holds the libraries only after `prepare.py vendor` (git-ignored, version in the path): check_offline
    local = [r for r in refs if not re.match(r"https?:", r) and not r.startswith("vendor/")]
    for r in local:
        name, _, query = r.partition("?")
        if name != "config.js":                      # config.js is each deployment's own file (git-ignored)
            c.ok((APP / name).exists(), f"index.html names {name}, which is not in web/SpatialViewer")
        c.ok(re.fullmatch(r"v=[\w.-]+", query) is not None, f"index.html: {name} has no ?v= string")
    names = [r.partition("?")[0] for r in refs]
    c.ok("diag.js" in names and "viewer.js" in names and names.index("diag.js") < names.index("viewer.js"),
         "index.html must load diag.js before viewer.js (it catches start-up errors)")
    js = (APP / "viewer.js").read_text(encoding="utf-8")
    order = [js.find(k) for k in ('["createImageBitmap", (b) => createImageBitmap(b, {',
                                  '["createImageBitmap without options"', '["an img element"')]
    c.ok(all(i > 0 for i in order) and order == sorted(order), "viewer.js: the three tile decoders are missing or out of order")
    for needle, what in (("onError: drawProblem", "deck errors go to the status line"),
                         ('addEventListener("webglcontextlost"', "a lost graphics context is reported"),
                         ("getTileData: [S.dsId, tileEpoch]", "failed tiles are asked for again"),
                         ("onTileError: () => {}", "tile errors are handled in getTileData (reported by fetchTile)"),
                         ("Math.min(cfg.tileCacheSize || touchMax, touchMax)", "config.js cannot lift the touch tile-cache limit")):
        c.ok(needle in js, f"viewer.js: {what} ({needle!r} not found)")
    dj = (APP / "diag.js").read_text(encoding="utf-8")
    c.ok('get("diag") !== "1"' in dj, "diag.js: the panel must stay off without ?diag=1")
    node = shutil.which("node")
    if not node:
        print("    node not found: tests/diag_test.js skipped")
        return
    p = subprocess.run([node, str(REPO / "tests" / "diag_test.js"), str(root / "export")], cwd=str(REPO),
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace")
    with open(log_path, "a", encoding="utf-8") as lg:
        lg.write(p.stdout)
    print("    " + (p.stdout.strip().splitlines()[-1] if p.stdout.strip() else "(no output)"))
    c.eq(p.returncode, 0, "node tests/diag_test.js (exit code)")
    for f in ("viewer.js", "diag.js", "lasso_worker.js", "umapfit.js", "leiden_worker.js"):
        q = subprocess.run([node, "--check", str(APP / f)], cwd=str(REPO), stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        c.eq(q.returncode, 0, f"node --check {f}")
    # old exports still open (SKIPPED without Playwright); and the current, packed ones (--packed): a gene by
    # Range, with a read pass on a ?data= link, and from a host that ignores Range
    packed_ids = [ds for ds in (VIS, MX, XE, RNA) if (root / "export" / ds / "manifest.json").exists()]
    p = subprocess.run([node, str(REPO / "tests" / "compat_test.js"), "--python", sys.executable,
                        "--packed", str(root / "export"), "--packed-ids", ",".join(packed_ids)], cwd=str(REPO),
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace")
    with open(log_path, "a", encoding="utf-8") as lg:
        lg.write(p.stdout)
    lines = p.stdout.strip().splitlines() or ["(no output)"]
    for ln in [x for x in lines if x.startswith("packed exports (")] + [lines[-1]]:
        print("    " + ln)
    c.eq(p.returncode, 0, "node tests/compat_test.js: the page opens exports made by older versions and packed current ones (exit code)")
    # the same page in Chromium, WebKit (Safari's engine), Edge, Firefox and phone-sized WebKit, Chromium and Firefox
    # (Firefox: the size only, no touch), on this run's demo export (SKIPPED without Playwright; a browser that is
    # not installed, such as Edge or Playwright's Firefox, is SKIPPED on its own and the rest still pass)
    p = subprocess.run([node, str(REPO / "tests" / "browsers_test.js"), "--export", str(root / "export"),
                        "--shots", str(root / "browsers_shots")], cwd=str(REPO),
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace")
    with open(log_path, "a", encoding="utf-8") as lg:
        lg.write(p.stdout)
    print_skipped(p.stdout)
    print("    " + (p.stdout.strip().splitlines() or ["(no output)"])[-1])
    c.eq(p.returncode, 0, "node tests/browsers_test.js: the page works in Chromium, WebKit, Edge and Firefox (exit code)")
    # Draw a region: the exact path against its fixture, fast mode in the worker and on this export (SKIPPED
    # without Playwright, or without umap-js from the CDN)
    p = subprocess.run([node, str(REPO / "tests" / "drawregion_test.js"), "--export", str(root / "export"),
                        "--shots", str(root / "browsers_shots")], cwd=str(REPO),
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace")
    with open(log_path, "a", encoding="utf-8") as lg:
        lg.write(p.stdout)
    print("    " + (p.stdout.strip().splitlines() or ["(no output)"])[-1])
    c.eq(p.returncode, 0, "node tests/drawregion_test.js: Draw a region's exact path = its fixture; fast mode labelled, saved and undone (exit code)")


THREAD_PROBE = r"""
import json, os, sys
sys.path.insert(0, sys.argv[1])
import threadcap
early = "numpy" in sys.modules                 # threadcap itself must not load numpy
import prepare                                 # caps first, then imports numpy (through pyramid, expr)
import numpy as np
from scipy import linalg
out = {"numpyBeforeCap": early, "cap": threadcap.cap(), "usable": threadcap.usable_cpus(),
       "env": {v: os.environ.get(v) for v in threadcap.VARS}, "describe": threadcap.describe()}
try:
    from threadpoolctl import threadpool_info
    out["pools"] = [p["num_threads"] for p in threadpool_info()]
except ImportError:
    out["pools"] = None
print("PROBE " + json.dumps(out))
"""


def check_threadcap(c: Checker, log_path: Path):
    """prepare.py caps the BLAS / OpenMP threads before numpy loads; the user's number wins;
    the cap does not change what the PCA computes."""
    import os
    sys.path.insert(0, str(APP))
    try:
        import sectionumap
        import threadcap
    finally:
        sys.path.pop(0)
    base = {k: v for k, v in os.environ.items() if k not in threadcap.VARS and k != threadcap.ENV}

    def probe(env: dict, *argv: str) -> dict | None:
        p = subprocess.run([sys.executable, "-c", THREAD_PROBE, str(APP), *argv], env={**base, **env},
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace")
        line = next((ln for ln in p.stdout.splitlines() if ln.startswith("PROBE ")), None)
        if p.returncode != 0 or line is None:
            print(p.stdout[-2000:], flush=True)
            return None
        return json.loads(line[6:])

    names = threadcap.VARS
    d = probe({})
    if c.ok(d is not None, "threadcap: probe ran (import prepare with no thread variables set)"):
        want = min(threadcap.DEFAULT_MAX, d["usable"])
        c.ok(not d["numpyBeforeCap"], "threadcap: importing threadcap does not load numpy")
        c.eq(d["cap"], want, "threadcap: default cap = min(8, usable CPUs)")
        c.ok(all(d["env"][v] == str(want) for v in names), f"threadcap: every thread variable set to {want}: {d['env']}")
        if d["pools"] is not None:
            c.ok(d["pools"] and max(d["pools"]) <= want, f"threadcap: loaded BLAS pools within the cap {want}: {d['pools']}")
        print(f"    default: {d['describe']}, pools {d['pools']}")
    d = probe({threadcap.ENV: "3"})
    if c.ok(d is not None, "threadcap: probe ran (SODA_THREADS=3)"):
        c.ok(d["cap"] == 3 and all(d["env"][v] == "3" for v in names), f"threadcap: SODA_THREADS=3 sets every variable to 3: {d['env']}")
        if d["pools"] is not None:
            c.ok(max(d["pools"]) <= 3, f"threadcap: SODA_THREADS=3 caps the loaded pools: {d['pools']}")
    d = probe({threadcap.ENV: "3", "OMP_NUM_THREADS": "6"}, "export", "--threads", "2")
    if c.ok(d is not None, "threadcap: probe ran (--threads 2)"):
        c.ok(d["cap"] == 2 and all(d["env"][v] == "2" for v in names), f"threadcap: --threads 2 wins over the environment: {d['env']}")
    d = probe({"OMP_NUM_THREADS": "5"})
    if c.ok(d is not None, "threadcap: probe ran (OMP_NUM_THREADS=5 already set)"):
        c.ok(d["env"]["OMP_NUM_THREADS"] == "5" and d["env"]["OPENBLAS_NUM_THREADS"] == str(d["cap"]),
             f"threadcap: a variable the user set is kept, the others get the cap: {d['env']}")
    log_txt = log_path.read_text(encoding="utf-8", errors="ignore") if log_path.exists() else ""
    if "prepare.py export" in log_txt:             # no export in the log: --checks-only on a tree made earlier
        c.ok(bool(re.search(r"\[export\] \d+ BLAS/OpenMP threads", log_txt)),
             "threadcap: the export log says how many BLAS/OpenMP threads it used")

    # The PCA runs on sectionumap.PCA_THREADS thread(s) whatever the cap, and gives the same
    # numbers as on the cap's threads (a seeded sparse matrix)
    try:
        from threadpoolctl import threadpool_info
    except ImportError:
        print("    threadpoolctl not installed: PCA thread checks skipped")
        return
    from scipy import sparse
    rng = np.random.default_rng(7)
    X = sparse.random(4000, 400, density=0.05, random_state=rng, data_rvs=lambda n: rng.poisson(2.0, n) + 1.0, format="csr")
    X.data = np.log1p(X.data).astype(np.float32)
    seen = []
    plain = sectionumap._pca_uncapped

    def spy(Xh, n_pcs):
        seen.append([p["num_threads"] for p in threadpool_info()])
        return plain(Xh, n_pcs)

    before = [p["num_threads"] for p in threadpool_info()]   # this process loaded numpy first: its own count
    sectionumap._pca_uncapped = spy
    try:
        a = sectionumap._pca(X, 30)                   # as the export calls it
    finally:
        sectionumap._pca_uncapped = plain
    c.ok(bool(seen) and seen[0] and max(seen[0]) <= sectionumap.PCA_THREADS,
         f"threadcap: the PCA runs on {sectionumap.PCA_THREADS} BLAS thread(s): pools inside it {seen[:1]}")
    after = [p["num_threads"] for p in threadpool_info()]
    c.eq(after, before, "threadcap: the pools are as they were once the PCA is done")
    b = plain(X, 30)                                  # on all of this process's threads
    diff = float(np.abs(a - b).max())
    c.ok(a.shape == (4000, 30) and diff <= 1e-4 * float(np.abs(a).max()),
         f"threadcap: PCA on {sectionumap.PCA_THREADS} thread(s) and on {max(before)} agree (max abs difference {diff:.3g})")
    print(f"    PCA on {sectionumap.PCA_THREADS} thread(s) vs {max(before)}: max abs difference {diff:.3g}"
          + (" (bit-identical)" if np.array_equal(a, b) else ""))


def check_picker(c: Checker):
    """The dataset picker is wired into the page (check_viewer_page covers the files and ?v= strings),
    carries its ARIA roles, and every .js of the page parses."""
    html = (APP / "index.html").read_text(encoding="utf-8")
    names = [r.partition("?")[0] for r in re.findall(r'<(?:script|link)[^>]+?(?:src|href)="([^"]+)"', html)]
    c.ok("dspicker.js" in names and "dspicker.css" in names, "index.html loads dspicker.js and dspicker.css")
    c.ok("dspicker.js" in names and "viewer.js" in names and names.index("dspicker.js") < names.index("viewer.js"),
         "dspicker.js loads before viewer.js")
    c.ok('id="dsPicker"' in html and 'id="dsLabel"' in html and 'id="dsSelect"' not in html,
         "index.html: the dataset picker's host and label, no <select id=dsSelect>")
    js = (APP / "dspicker.js").read_text(encoding="utf-8")
    for attr in ('role="combobox"', 'role="listbox"', '"role", "option"', '"role", "group"', 'aria-haspopup="listbox"',
                 "aria-activedescendant", "aria-selected", "aria-expanded", "aria-controls"):
        c.ok(attr in js, f"dspicker.js sets {attr}")
    viewer = (APP / "viewer.js").read_text(encoding="utf-8")
    c.ok("SodaDatasetPicker.create" in viewer and "openSeq" in viewer and "AbortController" in viewer,
         "viewer.js: picker created; dataset opens are numbered and abort the one before")
    node = shutil.which("node")
    if not node:
        print("    node not found: JavaScript syntax check skipped")
        return
    for f in sorted(APP.glob("*.js")):
        if f.name == "config.js":
            continue
        p = subprocess.run([node, "--check", str(f)], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        c.ok(p.returncode == 0, f"node --check {f.name}" + ("" if p.returncode == 0 else f": {p.stdout.strip()[:200]}"))


def check_layout(c: Checker):
    """The tabbed layout (layout.js, layout.css): the markup it needs is in index.html, every section
    belongs to a tab, no id is used twice, and the hooks in the other scripts are in place."""
    from html.parser import HTMLParser

    html = (APP / "index.html").read_text(encoding="utf-8")
    js = {f: (APP / f).read_text(encoding="utf-8") for f in ("layout.js", "viewer.js", "volcano.js", "plots.js", "qc.js", "roi.js", "snapshot.js", "cellsets.js", "binpaint.js", "selections.js", "measure.js")}
    names = [r.partition("?")[0] for r in re.findall(r'<(?:script|link)[^>]+?(?:src|href)="([^"]+)"', html)]
    c.ok("layout.js" in names and "layout.css" in names, "index.html loads layout.js and layout.css")
    c.ok("layout.js" in names and "viewer.js" in names and names.index("layout.js") < names.index("viewer.js"),
         "layout.js loads before viewer.js (and before the scripts that hand it their panels)")

    class Tags(HTMLParser):
        def __init__(self):
            super().__init__()
            self.tags = []

        def handle_starttag(self, tag, attrs):
            self.tags.append((tag, dict(attrs)))

    t = Tags()
    t.feed(html)
    ids = [a["id"] for _, a in t.tags if a.get("id")]
    dup = sorted({i for i in ids if ids.count(i) > 1})
    c.ok(not dup, f"index.html: an id is used twice: {dup}")
    by_id = {a["id"]: (tag, a) for tag, a in t.tags if a.get("id")}
    tabs = {"view": "tabView", "clusters": "tabClusters", "analyze": "tabAnalyze", "save": "tabSave"}
    c.ok(by_id.get("tabList", ("", {}))[1].get("role") == "tablist" and "tabRail" in by_id and "topBar" in by_id,
         "index.html: the tab rail (#tabRail, #tabList role=tablist) and the top bar (#topBar)")
    for key, tid in tabs.items():
        tag, a = by_id.get(tid, ("", {}))
        c.ok(tag == "button" and a.get("role") == "tab" and a.get("data-tab") == key and a.get("aria-selected") in ("true", "false"),
             f"index.html: #{tid} is a button with role=tab, data-tab={key} and aria-selected")
        box = by_id.get(a.get("aria-controls", ""), ("", {}))[1]
        c.ok(box.get("role") == "tabpanel" and box.get("aria-labelledby") == tid,
             f"index.html: #{tid} controls a role=tabpanel box that is labelled by it")
    c.eq(sum(1 for _, a in t.tags if a.get("role") == "tab" and a.get("aria-selected") == "true"), 1,
         "index.html: one tab is selected to begin with")
    # every section names one of the four tabs: in the page, and in the scripts that build one
    secs = [(a.get("id", "(no id)"), a.get("data-tab")) for tag, a in t.tags
            if tag == "section" and "panel" in (a.get("class") or "").split()]
    c.ok(len(secs) >= 10, f"index.html: the sections are there ({len(secs)} section.panel)")
    for sid, tab in secs:
        c.ok(tab in tabs and sid != "(no id)", f"index.html: section #{sid} has an id and names its tab (data-tab={tab!r})")
    made = re.findall(r'newSection\("(\w+)", "(\w+)"', js["layout.js"])
    c.ok(len(made) >= 4 and all(tab in tabs for _, tab in made), f"layout.js: the sections it makes name a tab: {made}")
    c.ok('sec.dataset.tab = "analyze"' in js["plots.js"] and "SodaLayout.adopt(sec)" in js["plots.js"],
         "plots.js: its section names its tab and is handed to the layout")
    c.ok("SodaLayout.adopt(block)" in js["volcano.js"], "volcano.js: its block is handed to the layout")
    c.ok('sec.dataset.tab = "save"' in js["snapshot.js"] and "SodaLayout.adopt(sec)" in js["snapshot.js"],
         "snapshot.js: its Save image section names its tab and is handed to the layout")
    c.ok("sec.dataset.tab = cfg.tab" in js["cellsets.js"] and "SodaLayout.adopt(sec)" in js["cellsets.js"]
         and 'id: "panelSets", tab: "clusters"' in js["cellsets.js"] and 'id: "panelRunSets", tab: "analyze"' in js["cellsets.js"],
         "cellsets.js: its Cell sets (Clusters) and Merge clusters of the run (Analyze) sections name their tab and are handed to the layout")
    c.ok('sec.dataset.tab = "view"' in js["binpaint.js"] and "SodaLayout.adopt(sec)" in js["binpaint.js"],
         "binpaint.js: its Colour by a value section names its tab and is handed to the layout")
    c.ok("binpaint.js" in names and not re.search(r"\bInt(8|16|32)Array\b", js["binpaint.js"]),
         "index.html loads binpaint.js, which uses only unsigned typed arrays (deid.py check reads the signed names as sample IDs)")
    c.ok('sec.dataset.tab = "view"' in js["measure.js"] and "SodaLayout.adopt(sec)" in js["measure.js"] and "measure.js" in names
         and names.index("measure.js") < names.index("viewer.js") and not re.search(r"\bInt(8|16|32)Array\b", js["measure.js"]),
         "measure.js: its Measure section names its tab and is handed to the layout, loaded before viewer.js, unsigned arrays only")
    c.ok('sec.dataset.tab = "analyze"' in js["selections.js"] and "SodaLayout.adopt(sec)" in js["selections.js"],
         "selections.js: its Combine and export shapes section names its tab and is handed to the layout")
    # ids the layout makes must be new; ids it moves must exist (in the page or in the script that builds them)
    new_ids = set(i for i, _ in made) | set(re.findall(r'el\("\w+", "[^"]*", "(\w+)"\)', js["layout.js"]))
    clash = sorted(new_ids & set(ids))
    c.ok(not clash, f"layout.js makes ids that index.html already has: {clash}")
    built = set(re.findall(r'id="(\w+)"', js["volcano.js"] + js["plots.js"] + js["snapshot.js"] + js["cellsets.js"] + js["binpaint.js"] + js["selections.js"] + js["measure.js"])) | set(ids) | new_ids
    block = js["layout.js"][js["layout.js"].index("const SECTIONS = {"):js["layout.js"].index("const ORDER = {")]
    moved = set(re.findall(r'"(\w+)"', " ".join(re.findall(r'(?:more|notes): \[([^\]]*)\]', block))))
    c.ok(len(moved) >= 15 and not (moved - built), f"layout.js moves controls that do not exist: {sorted(moved - built)}")
    order = js["layout.js"][js["layout.js"].index("const ORDER = {"):js["layout.js"].index("const NESTED")]
    ordered = set(re.findall(r'"(panel\w+)"', order))
    known = {sid for sid, _ in secs} | {i for i, _ in made} | {"panelPlots", "panelImage", "panelSets", "panelRunSets", "panelPaint", "panelShapes", "panelMeasure",
                                                               "panelOffline", "panelOwn"}
    c.ok(not (ordered - known), f"layout.js orders sections that do not exist: {sorted(ordered - known)}")
    for needle, f, what in (("lay.opening(false, tok)", "viewer.js", "tells the layout when a dataset starts and stops opening"),
                            ("SodaLayout.wording", "viewer.js", "hands its wording function to the layout"),
                            ('addEventListener("soda:shown"', "viewer.js", "redraws the gate histogram when its section is shown"),
                            ('addEventListener("soda:shown"', "qc.js", "redraws the QC histograms when their section is shown"),
                            ("#panelLasso h2 .sec-name", "roi.js", "finds the section title inside the layout's heading button"),
                            ('"role", "note"', "layout.js", "the ? notes carry role=note"),
                            ("aria-expanded", "layout.js", "section headers and More options carry aria-expanded"),
                            ("localStorage", "layout.js", "remembers the tab and the open sections")):
        c.ok(needle in js[f], f"{f}: {what} ({needle!r} not found)")
    # the page must work with storage blocked: every access inside a try
    lines = [ln for ln in js["layout.js"].splitlines() if "localStorage." in ln]
    c.ok(lines and all("try" in ln or "JSON.parse(localStorage" in ln for ln in lines),
         f"layout.js: a localStorage access outside try / catch: {[ln.strip()[:60] for ln in lines]}")
    # one control for each job (the duplicates are gone; old links that name them still open, viewlink.js)
    built_all = html + js["volcano.js"] + js["plots.js"]
    for gone in ("btnImage", "selShow", "selClear", "vcShow", "plShow", "lpSel", "lpSelShow", "lpSelClear"):
        c.ok(f'id="{gone}"' not in built_all and f'$("{gone}")' not in "".join(js.values()),
             f"#{gone} is back: Show on slide, Clear selection, the image's opacity and the run buttons are each one control")
    for label, where in (("Show on slide", "selChipShow"), ("Clear selection", "selChipClear")):
        c.eq(len(re.findall(rf">{label}</button>", built_all)), 1, f"index.html: one {label} button")
        c.ok(re.search(rf'id="{where}"[^>]*>{label}</button>', html) is not None, f"index.html: {label} is #{where}, in the chip above the picture")
    c.ok(not re.search(r">Show</button>", built_all), "index.html / volcano.js / plots.js: no separate Show button (the run button shows the last result)")
    c.ok("Show volcano" in js["volcano.js"] and "Show plots" in js["plots.js"], "volcano.js / plots.js: the run button reads Show … while nothing has changed")
    vl = (APP / "viewlink.js").read_text(encoding="utf-8")
    c.ok("view.img.on === false ? 0" in vl, "viewlink.js: a link made with the image off opens with Image opacity at 0")
    css = (APP / "layout.css").read_text(encoding="utf-8")
    c.ok("@media (max-width: 700px)" in css and ".app.tabbed" in css and ":focus-visible" in css,
         "layout.css: the phone layout, the .tabbed scope and a visible focus ring")


def check_ownres(c: Checker, root: Path, log_path: Path):
    """Your own results (ownres.js, ownres.css): a reader's CSV / TSV (also .gz) of results they already have,
    matched to the bins in the browser with no new export. Wired into the page (its section in the Clusters tab,
    the hooks in viewer.js and binpaint.js, the link extension), no signed typed arrays, and, when node is
    installed, tests/ownres_browser_test.js in Chromium, WebKit and Firefox (SKIPPED when not installed) on
    export_user/demo_visiumhd (barcode lists from
    `prepare.py ids`, check_ids): Seurat and scanpy files, .gz, partial matches, roles put right, colours by
    pixel, the volcano and cell sets on a loaded grouping, Save JSON and the view link, junk files, 375 px."""
    html = (APP / "index.html").read_text(encoding="utf-8")
    names = [r.partition("?")[0] for r in re.findall(r'<(?:script|link)[^>]+?(?:src|href)="([^"]+)"', html)]
    js = {f: (APP / f).read_text(encoding="utf-8") for f in ("ownres.js", "viewer.js", "binpaint.js", "layout.js")}
    c.ok("ownres.js" in names and "ownres.css" in names and "binpaint.js" in names
         and names.index("layout.js") < names.index("ownres.js") < names.index("viewer.js"),
         "index.html loads ownres.css and ownres.js (after layout.js, before viewer.js, which starts it)")
    o = js["ownres.js"]
    c.ok('sec.dataset.tab = "clusters"' in o and "SodaLayout.adopt(sec)" in o,
         "ownres.js: its Your own results section names its tab and is handed to the layout")
    c.ok(not re.search(r"Int(8|16|32)Array", o), "ownres.js uses only unsigned typed arrays (deid.py check reads the signed names as sample IDs)")
    c.ok('key: "own"' in o and "SodaViewLink.extend" in o, "ownres.js: a link to the view carries the file's fingerprint")
    c.ok('DecompressionStream("gzip")' in o, "ownres.js: reads gzipped files")
    for needle, f, what in (("ownInit()", "viewer.js", "starts Your own results"),
                            ("ownWith(base)", "viewer.js", "lists the file's groupings with the clusterings"),
                            ("c.own ? c.codes", "viewer.js", "reads a loaded grouping's codes from memory"),
                            ("own.saveBlock()", "viewer.js", "Save JSON can keep the matched table"),
                            ("own.restore(j.ownResults)", "viewer.js", "Load JSON brings it back"),
                            ("extraCols:", "viewer.js", "hands the file's numbers to Colour by a value"),
                            ("api.extraCols", "binpaint.js", "offers the file's numbers"),
                            ("panelOwn", "layout.js", "has the section's help and its place in the Clusters tab")):
        c.ok(needle in js[f], f"{f}: {what} ({needle!r} not found)")
    node = shutil.which("node")
    if not node:
        print("    node not found: tests/ownres_browser_test.js skipped")
        return
    p = subprocess.run([node, str(REPO / "tests" / "ownres_browser_test.js"), "--root", str(root),
                        "--shots", str(root / "ownres_shots"), "--browsers", "chromium,webkit,firefox"], cwd=str(REPO),
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace")
    with open(log_path, "a", encoding="utf-8") as lg:
        lg.write(p.stdout)
    print_skipped(p.stdout)
    print("    " + (p.stdout.strip().splitlines() or ["(no output)"])[-1])
    c.eq(p.returncode, 0, "node tests/ownres_browser_test.js: your own results loaded in the page, in Chromium, WebKit and Firefox (exit code)")


def check_sharing(c: Checker, log_path: Path):
    """The link that restores a view (viewlink.js) and linked windows (compare.js, compare.css): wired into
    the page, their ids, the hooks in viewer.js and layout.js, the page's Content-Security-Policy, and the
    node unit test of the link's encoding."""
    html = (APP / "index.html").read_text(encoding="utf-8")
    names = [r.partition("?")[0] for r in re.findall(r'<(?:script|link)[^>]+?(?:src|href)="([^"]+)"', html)]
    for f in ("viewlink.js", "compare.js", "compare.css"):
        c.ok(f in names and (APP / f).exists(), f"index.html loads {f}")
    for f in ("viewlink.js", "compare.js"):
        c.ok(f in names and names.index("layout.js") < names.index(f) < names.index("viewer.js")
             and names.index("dspicker.js") < names.index(f),
             f"{f} loads after layout.js and dspicker.js and before viewer.js (which hands it its functions)")
    js = {f: (APP / f).read_text(encoding="utf-8") for f in ("viewlink.js", "compare.js", "viewer.js", "layout.js")}
    page_ids = set(re.findall(r'\bid="(\w+)"', html))
    made = set(re.findall(r'\bid="(\w+)"', js["compare.js"])) | set(re.findall(r'\.id = "(\w+)"', js["compare.js"] + js["viewlink.js"]))
    want = {"panelCompare", "cmpOpenLabel", "cmpPicker", "cmpBlocked", "cmpBlockedLink", "cmpPeers", "cmpLink", "cmpAlign",
            "cmpCentre", "cmpState", "cmpGenes", "cmpColorBy", "cmpMode", "cmpMsg", "btnCompare", "saveLinkBox"}
    c.ok(want <= made, f"compare.js / viewlink.js build their controls: missing {sorted(want - made)}")
    c.ok(not (made & page_ids), f"ids built by compare.js / viewlink.js that index.html already has: {sorted(made & page_ids)}")
    lay_ids = set(re.findall(r'el\("\w+", "[^"]*", "(\w+)"\)', js["layout.js"]))
    c.ok({"btnCopyLink", "saveMsg"} <= lay_ids and not (made & lay_ids),
         "layout.js still builds the Copy link button and its message line (viewlink.js takes the button over)")
    c.ok('$("btnCopyLink")' in js["viewlink.js"] and "Copy link to this view" in js["viewlink.js"],
         "viewlink.js re-labels and re-wires the Copy link button")
    c.ok('sec.dataset.tab = "save"' in js["compare.js"] and "SodaLayout.adopt(sec)" in js["compare.js"],
         "compare.js: its section names its tab and is handed to the layout")
    c.ok("panelCompare: { help:" in js["layout.js"] and "panelCompare: () =>" in js["layout.js"],
         "layout.js: the Compare section's ? note and its summary")
    for needle, what in (("function sessionState(full)", "one serialiser for Save JSON and the link"),
                         ("JSON.stringify(sessionState(true), null, 1)", "Save JSON writes sessionState(true)"),
                         ("async function sessionApply(j, run)", "Load JSON and a link go back in through one function"),
                         ("window.SodaViewLink.capture()", "the view block is captured into the saved state"),
                         ("window.SodaViewLink.apply(j.view, j)", "the view block is applied from a file and from a link"),
                         ("function setCamera(viewState, viewId)", "one function moves the camera (mouse, link, linked window)"),
                         ("window.SodaViewLink.init(share)", "viewlink.js is handed viewer.js's functions"),
                         ("window.SodaCompare.init(share)", "compare.js is handed viewer.js's functions"),
                         ('new CustomEvent("soda:opened"', "the scripts are told when a dataset has opened"),
                         ("window.SodaCompare.askToken(", "a private site's new window asks its group for the token before the password box")):
        c.ok(needle in js["viewer.js"], f"viewer.js: {what} ({needle!r} not found)")
    # the page's Content-Security-Policy (no inline script, no handler attributes, no eval): nothing in the new files needs more
    for f in ("viewlink.js", "compare.js"):
        t = js[f]
        c.ok(not re.search(r"""<[a-z][^>`]*\son[a-z]+\s*=""", t) and "javascript:" not in t
             and not re.search(r"""\beval\(|new Function\(|setTimeout\(\s*["'`]""", t),
             f"{f}: no inline handler attribute, no eval, no javascript: address")
        c.ok(not re.search(r"\b(Int8Array|Int16Array|Int32Array|BigInt64Array)\b", t),  # deid:ok (type names)
             f"{f}: only unsigned typed arrays (the de-identification check reads the signed names as sample IDs)")
        lines = [ln for ln in t.splitlines() if re.search(r"\b(local|session)Storage\.", ln)]
        c.ok(all("try" in ln for ln in lines), f"{f}: a storage access outside try / catch: {[ln.strip()[:60] for ln in lines if 'try' not in ln]}")
    # the access token: memory only. What compare.js remembers for a tab is exactly this, and no address is built from a token
    saved = re.search(r"sessionStorage\.setItem\(KEY, JSON\.stringify\(\{([^}]*)\}\)\)", js["compare.js"])
    c.ok(bool(saved) and [k.strip() for k in saved.group(1).split(",")] == ["group", "claim", "linked", "off", "follow"],
         "compare.js remembers the group, the link state and the offset for a tab: no token, no password")
    c.ok(js["compare.js"].count("sessionStorage.setItem") == 1 and "localStorage" not in js["compare.js"],
         "compare.js writes to storage in one place only")
    c.ok("m.url !== H.tokenUrl || m.path !== location.pathname" in js["compare.js"] and "peers.has(m.from)" in js["compare.js"],
         "compare.js gives the token only to a member of its group on the same page and token address")
    c.ok("window.name" in js["compare.js"] and re.findall(r"""searchParams\.set\(\s*["'](\w+)""", js["compare.js"]) == ["ds"],
         "compare.js passes the group in window.name; the only thing it sets in an address is ?ds=")
    c.ok('KEEP_PARAMS = ["data", "static", "api", "fmt", "qc"]' in js["viewlink.js"] and "location.hash" in js["viewlink.js"]
         and "fetch(" not in js["viewlink.js"] and "fetch(" not in js["compare.js"],
         "viewlink.js: the view travels after #, the link keeps ?data= / ?static= / ?qc=, and neither file fetches anything itself")
    css = (APP / "compare.css").read_text(encoding="utf-8")
    c.ok("#btnCompare" in css and "@media (pointer: coarse)" in css, "compare.css: the top-bar button, hidden on touch screens")
    node = shutil.which("node")
    if not node:
        print("    node not found: tests/viewlink_test.js skipped")
        return
    p = subprocess.run([node, str(REPO / "tests" / "viewlink_test.js")], cwd=str(REPO),
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace")
    with open(log_path, "a", encoding="utf-8") as lg:
        lg.write(p.stdout)
    for ln in p.stdout.strip().splitlines()[-3:]:
        print("    " + ln.strip())
    c.eq(p.returncode, 0, "node tests/viewlink_test.js (exit code)")
    c.ok(bool(re.search(r"viewlink_test: \d+ checks passed", p.stdout)), "tests/viewlink_test.js: every check passed")


def check_link_files(c: Checker):
    """Encrypted link files (links.py, linkPath): marked deprecated in the docs, and still working."""
    url = "https://example.org/data/demo_visiumhd/"
    try:
        sys.path.insert(0, str(APP))
        import links as _links
        payload = _links.encrypt("a long test password", url)
        c.eq(_links.decrypt("a long test password", payload), url, "links.py: encrypt / decrypt still round-trip")
    except ImportError as e:
        print(f"    links.py round trip skipped ({e})")
    finally:
        sys.path.pop(0)
    for doc, text in (("docs/HOSTING.md", "Encrypted link files are deprecated"), ("docs/ADVANCED.md", "(no server; **deprecated**)"),
                      ("README.md", "deprecated: they will be removed"), ("deploy/azure/README.md", "> **Deprecated.**")):
        c.ok(text in (REPO / doc).read_text(encoding="utf-8"), f"{doc}: encrypted link files marked deprecated ({text!r})")


def check_offline(c: Checker, root: Path, log_path: Path):
    """Offline use: the pinned libraries agree everywhere, nothing widens the page's Content-Security-Policy, the
    service worker can find every file it must keep, and tests/offline_test.js passes in Chromium, WebKit and Firefox
    (Firefox SKIPPED when Playwright's Firefox is not installed)."""
    import tempfile
    html = (APP / "index.html").read_text(encoding="utf-8")
    js = {f: (APP / f).read_text(encoding="utf-8") for f in ("libs.js", "sw.js", "offline.js", "lasso_worker.js", "viewer.js", "layout.js", "diag.js")}
    m = re.search(r"/\* LIBS-JSON \*/(.*?)/\* END-LIBS-JSON \*/", js["libs.js"], re.S)
    libs = json.loads(m.group(1)) if m else {}
    c.ok(set(libs) == {"deck", "umap"} and all(re.fullmatch(r"sha384-[A-Za-z0-9+/]{64}", L.get("sri", "")) and L.get("bytes", 0) > 0
                                               and L["path"].startswith(f"{L['name']}@{L['version']}/") for L in libs.values()),
         f"libs.js: deck.gl and umap-js, each with a version, a path, a size and a SHA-384: {sorted(libs)}")
    tags = re.findall(r"<script\b[^>]*>", html)
    srcs = [re.search(r'src="([^"]+)"', s).group(1) for s in tags if 'src="' in s]
    deck = libs.get("deck", {})
    vtag = next((s for s in tags if 'src="vendor/' in s), "")
    c.ok(f'src="vendor/{deck.get("path")}"' in vtag and f'integrity="{deck.get("sri")}"' in vtag and "data-soda-lib" in vtag,
         "index.html: deck.gl from vendor/<path> with libs.js's SHA-384 as its integrity, marked data-soda-lib")
    names = [s.split("?")[0] for s in srcs]
    vname = f"vendor/{deck.get('path')}"
    c.ok(all(n in names for n in ("libs.js", vname, "diag.js", "config.js", "viewer.js"))
         and names.index("libs.js") == names.index(vname) + 1
         and names.index("diag.js") < names.index("libs.js") < names.index("config.js") < names.index("viewer.js"),
         "index.html: libs.js right after the vendor/ tag (it writes the CDN copy when that did not run), before config.js and viewer.js")
    c.ok("offline.js" in names and names.index("layout.js") < names.index("offline.js") < names.index("viewer.js"),
         "index.html: offline.js loads after layout.js and before viewer.js")
    c.ok(not any("cdn.jsdelivr.net" in s for s in srcs), "index.html names no CDN script itself (libs.js writes the fallback)")
    c.ok("cdn.jsdelivr.net/npm/" in js["libs.js"] and "integrity=" in js["libs.js"] and 'crossorigin="anonymous"' in js["libs.js"],
         "libs.js: the CDN fallback carries the same SHA-384 (integrity, crossorigin)")
    csp = re.search(r'http-equiv="Content-Security-Policy" content="([^"]+)"', html)
    csp = csp.group(1) if csp else ""
    c.ok("script-src 'self' https://cdn.jsdelivr.net;" in csp and "worker-src 'self' blob:;" in csp and "unsafe-eval" not in csp,
         f"index.html: the Content-Security-Policy is no wider than before (no 'unsafe-eval', no 'wasm-unsafe-eval'): {csp[:120]}")
    for f in ("libs.js", "sw.js", "offline.js"):
        t = js[f]
        c.ok(not re.search(r"""\beval\(|new Function\(|setTimeout\(\s*["'`]""", t) and "javascript:" not in t
             and not re.search(r"\b(Int8Array|Int16Array|Int32Array|BigInt64Array)\b", t),  # deid:ok (type names)
             f"{f}: no eval, no javascript: address, only unsigned typed arrays")
        lines = [ln for ln in t.splitlines() if re.search(r"\b(local|session)Storage\.", ln)]
        c.ok(not lines, f"{f}: no local / session storage (offline copies live in the Cache API only): {lines[:2]}")
    w = js["lasso_worker.js"]
    c.ok("https://" not in "".join(re.findall(r"importScripts\([^)]*\)", w)) and "SODA_LIBS.umap" in w and "integrity: L.sri" in w,
         "lasso_worker.js: umap-js from libs.js (vendor/ first, the CDN after), its SHA-384 checked with fetch() before importScripts")
    v_index = re.search(r'libs\.js\?v=([\w.-]+)"', html)
    v_worker = re.search(r'importScripts\("libs\.js\?v=([\w.-]+)"\)', w)
    c.ok(bool(v_index and v_worker and v_index.group(1) == v_worker.group(1)),
         "lasso_worker.js imports the same libs.js ?v= as index.html loads")
    # its helpers: umapfit.js (umap-js's fit over typed arrays) and leiden_worker.js (the stability sweep's runs)
    c.ok(re.search(r'importScripts\("umapfit\.js\?v=[\w.-]+"\)', w) is not None and re.search(r'"leiden_worker\.js\?v=[\w.-]+"', w) is not None
         and (APP / "umapfit.js").exists() and (APP / "leiden_worker.js").exists(),
         "lasso_worker.js names umapfit.js and leiden_worker.js with ?v= strings, and both are in the page folder")
    lw = (APP / "leiden_worker.js").read_text(encoding="utf-8")
    v_l1 = re.search(r'importScripts\("leiden\.js\?v=([\w.-]+)"\)', w)
    v_l2 = re.search(r'importScripts\("leiden\.js\?v=([\w.-]+)"\)', lw)
    c.ok(bool(v_l1 and v_l2 and v_l1.group(1) == v_l2.group(1)), "lasso_worker.js and leiden_worker.js import the same leiden.js ?v=")
    # the service worker keeps what a script names as "<name>.js?v=...": every worker and import must be named that way
    allj = "".join((APP / f).read_text(encoding="utf-8") for f in sorted(os.listdir(APP))
                   if f.endswith(".js") and f not in ("config.js", "sw.js"))
    wk = re.findall(r"""["'`]([\w-]+_worker\.js)(\?v=[\w.-]+)?["'`]""", allj)
    c.ok(len(wk) >= 3 and all(v for _, v in wk), f"every worker is named with a ?v= string (sw.js keeps it): {sorted(set(n + v for n, v in wk))}")
    imp = re.findall(r"""importScripts\(\s*["'`]([^"'`]+)["'`]""", "".join((APP / f).read_text(encoding="utf-8") for f in sorted(os.listdir(APP)) if f.endswith("_worker.js")))
    c.ok(imp and all(re.fullmatch(r"[\w-]+\.js\?v=[\w.-]+", i) for i in imp),
         f"every importScripts of the page's workers names a local script with a ?v= string: {imp}")
    mk = re.search(r"const MARKERS = \[([^\]]*)\]", js["sw.js"])
    listed = sorted(re.findall(r'"([\w.-]+\.json)"', mk.group(1))) if mk else []
    c.eq(listed, sorted(p.name for p in (APP / "markers").glob("*.json")), "sw.js MARKERS = the files in markers/")
    body = lambda t: re.sub(r"\s+", " ", (re.search(r"function hashRefs\(list\) \{.*?return [^\n]*\n", t, re.S) or re.search("x^", "")).group(0)) \
        if re.search(r"function hashRefs\(list\) \{", t) else None
    c.ok(body(js["sw.js"]) is not None and body(js["sw.js"]) == body(js["offline.js"]), "sw.js and offline.js hash the page's ?v= strings the same way")
    rule = r"""if (!r || /^[a-z][\w+.-]*:|^\/\//i.test(r) || /^vendor\//.test(r) || !/[?&]v=/.test(r)) continue;"""
    c.ok(rule.replace("!r || ", "") in js["sw.js"] and rule in js["offline.js"], "sw.js and offline.js pick the page's ?v= files by the same rule")
    for needle, f, what in (("soda-data:${SCOPE}|", "sw.js", "kept datasets in caches of their own, per page folder"),
                            ("name.startsWith(`soda-app:${SCOPE}|`)", "sw.js", "on activate, only this folder's own app caches are deleted"),
                            ('if (req.method !== "GET") return;', "sw.js", "only GET requests are answered"),
                            ("e.respondWith(fromKept(req, k))", "sw.js", "data: only a kept dataset's files are answered"),
                            ("!cfg.tokenEndpoint", "offline.js", "off by default on a password-protected site"),
                            ('why: "password"', "viewer.js", "a password-protected dataset is never kept"),
                            ("window.SodaOffline.init(", "viewer.js", "viewer.js hands offline.js the open dataset"),
                            ('sec.dataset.tab = "save"', "offline.js", "its section names its tab"),
                            ("SodaLayout.adopt(sec)", "offline.js", "its section is handed to the layout"),
                            ("panelOffline: { help:", "layout.js", "the Offline use section's ? note"),
                            ("el.dataset.sodaLib", "diag.js", "the vendor/ probe is not reported as an error")):
        c.ok(needle in js[f], f"{f}: {what} ({needle!r} not found)")
    gi = (REPO / ".gitignore").read_text(encoding="utf-8")
    c.ok("web/SpatialViewer/vendor/" in gi.splitlines(), ".gitignore: web/SpatialViewer/vendor/ (the library files are not committed)")
    prep = (APP / "prepare.py").read_text(encoding="utf-8")
    c.ok('raw[:1] == ["vendor"]' in prep and "vendor.main(raw[1:])" in prep, "prepare.py vendor runs vendor.py")
    # vendor.py: reads libs.js; --check on a folder without vendor/ says MISSING and exits 1 (downloads nothing)
    with tempfile.TemporaryDirectory() as td:
        p = subprocess.run([sys.executable, str(APP / "vendor.py"), "--check", "--dest", td], cwd=str(REPO),
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace")
        c.ok(p.returncode == 1 and p.stdout.count("MISSING") == len(libs) and not (Path(td) / "vendor").exists(),
             f"vendor.py --check on an empty folder: MISSING for each library, exit 1, nothing written ({p.returncode})")
    node = shutil.which("node")
    if not node:
        print("    node not found: tests/offline_test.js skipped")
        return
    p = subprocess.run([node, str(REPO / "tests" / "offline_test.js"), "--export", str(root / "export"),
                        "--shots", str(root / "offline_shots"), "--python", sys.executable], cwd=str(REPO),
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace")
    with open(log_path, "a", encoding="utf-8") as lg:
        lg.write(p.stdout)
    lines = p.stdout.strip().splitlines() or ["(no output)"]
    tab = lines.index("feature table:") if "feature table:" in lines else len(lines) - 1
    print_skipped("\n".join(lines[:tab]))
    for ln in lines[tab:]:
        print("    " + ln)
    c.eq(p.returncode, 0, "node tests/offline_test.js: offline in Chromium, WebKit and Firefox (exit code)")


def check_publishable(c: Checker, root: Path, log_path: Path):
    export = root / "export"
    texts = [p for ds in (VIS, MX, TL) + XEN_DS + XIN_DS + GEOMX_DS + STD_DS for p in (export / ds).rglob("*")
             if p.is_file() and p.suffix in (".json", ".geojson", ".txt")]
    texts += [p for p in (root / "export_sample" / MX).rglob("*") if p.is_file() and p.suffix == ".json"]
    texts += [p for ds in (VIS, XE) for p in (root / "export_user" / ds).rglob("*")
              if p.is_file() and p.suffix in (".json", ".geojson")]
    bc = [p for p in texts if BARCODE.search(p.read_text(encoding="utf-8", errors="ignore"))]
    c.ok(not bc, f"export text files contain barcodes: {bc[:3]}")
    lp = [p for p in texts if LOCAL_PATH.search(p.read_text(encoding="utf-8", errors="ignore"))]
    c.ok(not lp, f"export text files contain local paths: {lp[:3]}")
    args = (["web/SpatialViewer/deid.py", "check", str(export / VIS), str(export / MX), str(export / TL)] +
            [str(export / ds) for ds in XEN_DS + GEOMX_DS] + [str(export / ds) for ds in STD_DS + XIN_DS if (export / ds).exists()])
    if (root / "export_sample" / MX).exists():
        args.append(str(root / "export_sample" / MX))
    for extra in (root / "export_user" / VIS, root / "export_user" / XE, root / "export_partial" / VIS, root / "export" / RNA,
                  root / "export_loose" / VIS, root / "export_both" / XE):
        if extra.exists():
            args.append(str(extra))
    print(f"$ python {' '.join(args)}", flush=True)
    p = subprocess.run([sys.executable, *args], cwd=str(REPO), stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                       text=True, encoding="utf-8", errors="replace")
    with open(log_path, "a", encoding="utf-8") as lg:
        lg.write(p.stdout)
    print("    " + p.stdout.strip().splitlines()[-1] if p.stdout.strip() else "    (no output)")
    c.eq(p.returncode, 0, "deid.py check on the export (exit code)")


# ------------------------------------------------------------------ gene packs
def _http(url: str, method: str = "GET", headers: dict | None = None):
    """(status, headers, body) of one request; HTTP errors are answers too."""
    import urllib.error
    import urllib.request
    req = urllib.request.Request(url, method=method, headers=headers or {})
    try:
        with urllib.request.urlopen(req, timeout=30) as r:
            return r.status, r.headers, r.read()
    except urllib.error.HTTPError as e:
        return e.code, e.headers, e.read()


def write_rnaseq_demo(root: Path, log_path: Path) -> dict:
    """A small synthetic RNA-seq set, packed (root/export/demo_rnaseq) and loose
    (root/export_rnaseq_loose/demo_rnaseq), by rnaseq.export_rnaseq in this process."""
    from scipy import sparse
    rnaseq = app_module("rnaseq")
    rng = np.random.default_rng(7)
    n_cells, n_genes = 400, 60
    dense = rng.poisson(rng.gamma(0.6, 2.0, size=(1, n_genes)), size=(n_cells, n_genes)).astype(np.float32)
    dense[:, 5] = 0                                   # a gene nobody expresses: listed nowhere, no payload
    dense[0, 0] = 40                                  # integer counts above 20: normalised to log1p CP10k
    genes = [f"Rg{i:03d}" for i in range(n_genes)]
    genes[7] = "TRB1/2"                               # a name that is encoded in a loose file name
    emb = rng.normal(size=(n_cells, 2)).astype(np.float32)
    obs = {"cluster": np.array([f"c{k % 4}" for k in range(n_cells)], dtype=object)}
    X = sparse.csr_matrix(dense)
    with open(log_path, "a", encoding="utf-8") as lg:
        say = lambda *a: lg.write(" ".join(map(str, a)) + "\n")
        rnaseq.export_rnaseq(X, genes, emb, obs, root / "export", RNA, "Synthetic RNA-seq", log=say)
        rnaseq.export_rnaseq(X, genes, emb, obs, root / "export_rnaseq_loose", RNA, "Synthetic RNA-seq", log=say,
                             gene_files="loose")
    return {"cells": n_cells, "genes": n_genes, "expressed": int((dense.sum(0) > 0).sum())}


def check_packs(c: Checker, root: Path, log_path: Path, rna_demo: dict) -> dict:
    """Gene packs (genepack.py), the default since they replaced one file per gene:
    - every current export is packed, with no gene/ folder, and every pack agrees with its listing
      byte for byte (genepack.check_group with the content hash re-read);
    - --gene-files loose writes the old layout, and each of its files equals the gene's slice of the
      pack; --gene-files both writes the two side by side, equal again, in shards of 0.05 MB: several
      shards, a gene never split, no shard over the size unless it holds a single gene;
    - static_serve.py and server.py answer a Range request with 206, the asked bytes (equal to the
      loose file) and Content-Range, also with a read pass (?sv=..&sig=..) on the URL; a plain GET
      gets the whole shard, a preflight is allowed, a range past the end is refused (416);
    - genepack.py repack turns the loose export into the same packs the export writes, and deletes
      the files only after checking them; an RNA-seq set packs like the rest;
    - deid.py check on a pack with an extra byte or an extra listing field reports it."""
    gp = app_module("genepack")
    deid = app_module("deid")
    export = root / "export"
    res = {"rnaseqDemo": rna_demo}
    # ---- every current export: packed, consistent ----------------------------------------------
    # (the cell-outline exports too: their shapes/ tiles sit beside the packs, unpacked and untouched)
    ids = [(export, ds) for ds in (VIS, MX, TL, RNA) + XEN_DS + XIN_DS + GEOMX_DS + ST_DS + IMG_DS + STD_DS
           if (export / ds / "manifest.json").exists()]
    ids += [(root / r, ds) for r in (OL_ZARR, OL_OWN, OL_NONE) for ds in XEN_DS if (root / r / ds / "manifest.json").exists()]
    n_groups = n_shards = n_genes = n_shapes = 0
    for base, ds in ids:
        m = read_json(base / ds / "manifest.json")
        for L in m.get("expr", []):
            if L.get("shapes"):
                n_shapes += 1
                sd = base / ds / "expr" / L["key"] / "shapes"
                c.ok(sd.is_dir() and not any(sd.rglob("*.pack")), f"{base.name}/{ds} {L['key']}: shapes/ present, nothing packed in it")
        for g in gp.groups(m):
            gdir = base / ds / g["dir"] if g["dir"] else base / ds
            c.ok(bool(g["pack"]), f"{base.name}/{ds} {g['what']}: not packed by default")
            c.ok(not (gdir / "gene").exists(), f"{base.name}/{ds} {g['what']}: a gene/ folder in a packed export")
            probs = gp.check_group(base / ds, g, verify_hash=True)
            c.ok(not probs, f"{base.name}/{ds} {g['what']}: pack vs listing: {probs[:3]}")
            if g["pack"]:
                n_groups += 1
                n_shards += len(g["pack"]["shards"])
                n_genes += g["pack"]["genes"]
    res["current"] = {"datasets": len(ids), "groups": n_groups, "shards": n_shards, "genes": n_genes, "layersWithShapes": n_shapes}
    print(f"    {len(ids)} exports: {n_groups} packed groups, {n_genes:,} genes in {n_shards} shards; {n_shapes} layer(s) with outlines")

    # ---- prepare.py batch: its exports are packed and consistent too; --gene-files / --shard-mb pass through
    bx = root / "batch" / "export"
    bids = sorted(p.parent.name for p in bx.glob("*/manifest.json")) if bx.exists() else []
    nb_groups = 0
    for ds in bids:
        m = read_json(bx / ds / "manifest.json")
        for g in gp.groups(m):
            nb_groups += 1
            c.ok(bool(g["pack"]) and not gp.check_group(bx / ds, g, verify_hash=True), f"batch export {ds} {g['what']}: packed and consistent")
    bm = app_module("batch")
    ba = bm.build_parser().parse_args(["--src", "s", "--out", "o", "--dest", "d", "--gene-files", "loose", "--shard-mb", "32"])
    d0 = {"id": "x", "kind": "xenium", "variant": "", "opts": {}, "user": None, "src": "s", "forceKind": False}
    ex, rn = bm.export_argv(d0, ba), bm.prepare_argv({**d0, "kind": "rnaseq"}, ba)
    c.ok(ex[ex.index("--gene-files") + 1] == "loose" and ex[ex.index("--shard-mb") + 1] == "32.0"
         and "--gene-files" in rn and "--gene-files" not in bm.prepare_argv(d0, ba)
         and "--gene-files" not in bm.export_argv(d0, bm.build_parser().parse_args(["--src", "s", "--out", "o", "--dest", "d"])),
         f"batch --gene-files / --shard-mb: on every export and the RNA-seq step, nowhere else, absent by default ({ex})")
    res["batch"] = {"datasets": bids, "groups": nb_groups}
    print(f"    batch: {len(bids)} exports, {nb_groups} packed groups; --gene-files passes through")

    # ---- --gene-files loose: the old layout, equal to the packs byte for byte ------------------
    lo, pk = root / "export_loose" / VIS, export / VIS
    ml, mp = read_json(lo / "manifest.json"), read_json(pk / "manifest.json")
    c.ok(not any(L.get("genePack") for L in ml.get("expr", [])) and not (lo / "expr" / "008um" / "genepack").exists(),
         "export --gene-files loose: no genePack, no genepack/")
    n_eq = n_files = 0
    bad = []
    for L in mp.get("expr", []):
        key = L["key"]
        gl, gpk = read_json(lo / "expr" / key / "genes.json"), read_json(pk / "expr" / key / "genes.json")
        c.eq([(g["n"], g["t"]) for g in gl], [(g["n"], g["t"]) for g in gpk], f"loose vs packed {key}: genes.json names and totals")
        c.ok(not any(set(g) - {"n", "t"} for g in gl), f"loose {key}: genes.json entries carry only n, t")
        files = {f.name for f in (lo / "expr" / key / "gene").glob("*.bin")}
        n_files += len(files)
        c.eq(len(files), L["genePack"]["genes"], f"loose {key}: one file per packed gene")
        for g in gpk:
            a, b = gene_bytes(lo / "expr" / key, g["n"]), gene_bytes(pk / "expr" / key, g["n"])
            if a != b:
                bad.append(f"{key}/{g['n']}")
            elif a is not None:
                n_eq += 1
    c.ok(not bad and n_eq, f"loose files = pack slices byte for byte ({n_eq} genes; differ: {bad[:4]})")
    res["loose"] = {"files": n_files, "equal": n_eq}

    # ---- --gene-files both, tiny shards ----------------------------------------------------------
    bo = root / "export_both" / XE
    mb = read_json(bo / "manifest.json")
    multi = 0
    for g in gp.groups(mb):
        gdir = bo / g["dir"] if g["dir"] else bo
        if not c.ok(bool(g["pack"]) and (gdir / "gene").is_dir(), f"both {g['what']}: pack and gene/ files"):
            continue
        c.ok(not gp.check_group(bo, g, verify_hash=True), f"both {g['what']}: pack vs listing")
        sh = g["pack"]["shards"]
        multi += len(sh) > 1
        per = {}
        for e in gp.read_listing(bo, g):
            if "s" in e:
                per.setdefault(e["s"], []).append(e["b"])
        lim = 0.05 * (1 << 20)
        c.ok(all(x["bytes"] <= lim or len(per.get(i, [])) == 1 for i, x in enumerate(sh)),
             f"both {g['what']}: a shard over 0.05 MB holds more than one gene")
        c.ok(all(sum(per[i][:-1]) + 0 < lim for i in per if len(per[i]) > 1),
             f"both {g['what']}: a shard was filled past its size before the last gene")
        diff = [e["n"] for e in gp.read_listing(bo, g) if "s" in e and
                gp.gene_bytes(bo, g, e) != (gdir / "gene" / f"{enc(e['n'])}.bin").read_bytes()]
        c.ok(not diff, f"both {g['what']}: pack slices = files ({diff[:3]})")
        # the same genes as the default export (64 MB shards): sharding changes where, not what
        ref = export / XE / g["dir"]
        diff = [e["n"] for e in gp.read_listing(bo, g) if gene_bytes(ref, e["n"]) != gene_bytes(gdir, e["n"])]
        c.ok(not diff, f"both {g['what']}: same genes as {XE} in one shard ({diff[:3]})")
    c.ok(multi >= 2, f"both: several shards in the cells layer and the molecules ({multi} groups with more than one)")

    # ---- HTTP Range: static_serve.py ---------------------------------------------------------------
    import threading
    ss = app_module("static_serve")
    srv = ss.serve(root, 0)
    th = threading.Thread(target=srv.serve_forever, daemon=True)
    th.start()
    base = f"http://127.0.0.1:{srv.server_address[1]}"
    try:
        res["staticServe"] = _range_checks(c, base, root, "static_serve.py")
    finally:
        srv.shutdown()
        srv.server_close()

    # ---- HTTP Range: server.py (FastAPI StaticFiles at /, the repo's web/ in real use) -------------
    if all(importlib.util.find_spec(mod) for mod in ("fastapi", "uvicorn", "zarr")):
        import socket
        with socket.socket() as so:
            so.bind(("127.0.0.1", 0))
            port = so.getsockname()[1]
        proc = subprocess.Popen([sys.executable, str(APP / "server.py"), "--cache", str(root / "cache"), "--web", str(root),
                                 "--no-password", "--host", "127.0.0.1", "--port", str(port)], cwd=str(APP),
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace")
        try:
            up = False
            for _ in range(120):
                try:
                    with socket.create_connection(("127.0.0.1", port), timeout=1):
                        up = True
                        break
                except OSError:
                    time.sleep(0.25)
            if c.ok(up, "server.py started"):
                res["serverPy"] = _range_checks(c, f"http://127.0.0.1:{port}", root, "server.py")
        finally:
            proc.terminate()
            try:
                out, _ = proc.communicate(timeout=10)
            except subprocess.TimeoutExpired:
                proc.kill()
                out, _ = proc.communicate()
            with open(log_path, "a", encoding="utf-8") as lg:
                lg.write(out or "")
    else:
        print("    SKIP server.py Range check (fastapi / uvicorn not installed)")

    # ---- repack: the loose export -> the export's own packs --------------------------------------
    rp = root / "export_repack" / VIS
    shutil.copytree(lo, rp)
    p = subprocess.run([sys.executable, str(APP / "genepack.py"), "repack", str(rp), "--delete-loose"], cwd=str(REPO),
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace")
    with open(log_path, "a", encoding="utf-8") as lg:
        lg.write(p.stdout)
    c.eq(p.returncode, 0, "genepack.py repack (exit code)")
    mr = read_json(rp / "manifest.json")
    for L in mr.get("expr", []):
        ref = next(x for x in mp["expr"] if x["key"] == L["key"])
        c.eq(L.get("genePack"), ref["genePack"], f"repack {L['key']}: the same shards (names, sizes) as the export's own packs")
        c.eq(read_json(rp / "expr" / L["key"] / "genes.json"), read_json(pk / "expr" / L["key"] / "genes.json"),
             f"repack {L['key']}: the same genes.json as the export's")
        c.ok(not (rp / "expr" / L["key"] / "gene").exists(), f"repack {L['key']}: --delete-loose removed gene/")
    c.ok(not gp.check_dataset(rp, verify_hash=True), "repack: packs consistent")

    # ---- RNA-seq set: packed and loose agree; repack of the loose one ------------------------------
    rs, rl = export / RNA, root / "export_rnaseq_loose" / RNA
    mr = read_json(rs / "manifest.json")
    gj = read_json(rs / "genes.json")
    c.ok(bool(mr.get("genePack")) and not (rs / "gene").exists() and (rl / "gene").is_dir() and not read_json(rl / "manifest.json").get("genePack"),
         "rnaseq: packed by default (genePack, no gene/), loose with gene_files='loose'")
    c.eq(sum(1 for e in gj if "s" in e), len(list((rl / "gene").glob("*.bin"))), "rnaseq: one pack entry per loose gene file")
    diff = [e["n"] for e in gj if gene_bytes(rs, e["n"]) != gene_bytes(rl, e["n"])]
    c.ok(not diff, f"rnaseq: pack slices = loose files ({diff[:3]})")
    c.ok((rl / "gene" / "TRB1%2F2.bin").exists() and gene_bytes(rs, "TRB1/2") == (rl / "gene" / "TRB1%2F2.bin").read_bytes(),
         "rnaseq: 'TRB1/2' in the pack = its encoded loose file")
    c.ok(not gp.check_dataset(rs, verify_hash=True), "rnaseq: pack consistent")
    k = np.frombuffer(gene_bytes(rs, gj[0]["n"])[:4], "<u4")[0]
    c.ok(abs(k / mr["nCells"] - gj[0]["p"]) < 1e-3, f"rnaseq: {gj[0]['n']} k = {k} matches p = {gj[0]['p']}")
    rr = root / "export_repack" / RNA
    shutil.copytree(rl, rr)
    gp.repack(rr, log=lambda *a: None)
    diff = [e["n"] for e in read_json(rr / "genes.json") if gene_bytes(rr, e["n"]) != gene_bytes(rl, e["n"])]
    c.ok(read_json(rr / "manifest.json").get("genePack") and not diff, f"rnaseq repack: slices = files ({diff[:3]})")

    # ---- deid.py check covers the packs ---------------------------------------------------------
    c.eq(deid.check_packs([export / RNA / "manifest.json"]), [], "deid check_packs: a good pack is clean")
    tp = root / "export_tamper" / RNA
    shutil.copytree(rs, tp)
    sh0 = tp / mr["genePack"]["dir"] / mr["genePack"]["shards"][0]["file"]
    with open(sh0, "ab") as f:
        f.write(b"\0" * 12)                          # bytes no listing entry accounts for
    probs = deid.check_paths([tp])
    c.ok(any("bytes" in x and "pack" in x for x in probs), f"deid check: an extra 12 bytes in a shard is reported ({probs[:2]})")
    with open(sh0, "r+b") as f:
        f.truncate(mr["genePack"]["shards"][0]["bytes"])
    gl = read_json(tp / "genes.json")
    gl[0]["sample"] = "a sample label"                # a field a listing must not carry
    (tp / "genes.json").write_text(json.dumps(gl), encoding="utf-8")
    probs = deid.check_paths([tp])
    c.ok(any("fields" in x and "sample" in x for x in probs), f"deid check: an extra listing field is reported ({probs[:2]})")
    return res


def _range_checks(c: Checker, base: str, root: Path, what: str) -> dict:
    """Range requests into the demo_visiumhd packs, answered by the server at `base` serving `root`."""
    pk, lo = root / "export" / VIS, root / "export_loose" / VIS
    m = read_json(pk / "manifest.json")
    n = 0
    bad = []
    for L in m["expr"]:
        key, pack = L["key"], L["genePack"]
        gl = [e for e in read_json(pk / "expr" / key / "genes.json") if "s" in e]
        pick = {gl[0]["n"]: gl[0], gl[-1]["n"]: gl[-1], gl[len(gl) // 2]["n"]: gl[len(gl) // 2]}
        big = max(gl, key=lambda e: e["b"])
        pick[big["n"]] = big
        for e in pick.values():
            sh = pack["shards"][e["s"]]
            url = f"{base}/export/{VIS}/expr/{key}/{pack['dir']}/{sh['file']}?sv=2099-01-01&sig=FAKE%2Fsig%3D"
            st, hd, body = _http(url, headers={"Range": f"bytes={e['o']}-{e['o'] + e['b'] - 1}", "Origin": "http://example.org"})
            want = (lo / "expr" / key / "gene" / f"{enc(e['n'])}.bin").read_bytes()
            cr = (hd.get("Content-Range") or "").strip()
            if st != 206 or body != want or cr != f"bytes {e['o']}-{e['o'] + e['b'] - 1}/{sh['bytes']}":
                bad.append(f"{key}/{e['n']}: {st}, {len(body)} bytes (file {len(want)}), Content-Range {cr!r}")
            n += 1
        expose = (hd.get("Access-Control-Expose-Headers") or "").lower()
        c.ok("content-range" in expose and hd.get("Access-Control-Allow-Origin") in ("*", "http://example.org"),
             f"{what}: CORS on a ranged answer lets a page read Content-Range ({expose!r})")
    c.ok(not bad, f"{what}: Range requests with a read pass give 206 and the gene's file bytes: {bad[:3]}")
    L = m["expr"][0]
    sh = L["genePack"]["shards"][0]
    url = f"{base}/export/{VIS}/expr/{L['key']}/{L['genePack']['dir']}/{sh['file']}"
    st, hd, body = _http(url)
    c.ok(st == 200 and len(body) == sh["bytes"], f"{what}: a plain GET gets the whole shard ({st}, {len(body)} bytes)")
    st, hd, _ = _http(url, headers={"Range": f"bytes={sh['bytes'] + 10}-{sh['bytes'] + 20}"})
    c.eq(st, 416, f"{what}: a range past the end is refused")
    st, hd, _ = _http(url, method="OPTIONS", headers={"Origin": "http://example.org", "Access-Control-Request-Method": "GET",
                                                      "Access-Control-Request-Headers": "range"})
    c.ok(st in (200, 204) and "range" in (hd.get("Access-Control-Allow-Headers") or "").lower()
         or hd.get("Access-Control-Allow-Headers") == "*", f"{what}: a CORS preflight for Range is allowed ({st})")
    print(f"    {what}: {n} genes fetched by Range = their loose files")
    return {"genes": n}


# ------------------------------------------------------------------ main
def main(argv=None) -> int:
    for stream in (sys.stdout, sys.stderr):      # a redirected Windows console (cp1252) must not stop the test
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(errors="replace")
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=str(REPO / "demo"), help="where raw/, cache/, export/ go (default <repo>/demo)")
    ap.add_argument("--reuse", action="store_true", help="keep earlier outputs (export skips existing files)")
    ap.add_argument("--checks-only", action="store_true", help="do not run the commands, only check the outputs")
    args = ap.parse_args(argv)
    for k, v in (("NODE_PATH", r"D:\SpatialVizBench\iphone\pw\node_modules"),
                 ("PLAYWRIGHT_BROWSERS_PATH", r"D:\SpatialVizBench\iphone\browsers")):
        if not os.environ.get(k) and Path(v).is_dir():   # unset: use this machine's bench Playwright, so no browser test skips
            os.environ[k] = v
            print(f"{k} was not set: using {v}")
    root = Path(args.root)
    root.mkdir(parents=True, exist_ok=True)
    log_path = root / "smoke_test.log"
    clustering = all(importlib.util.find_spec(mod) is not None for mod in CLUSTER_MODULES)
    print(f"python {sys.version.split()[0]} at {sys.executable}")
    print(f"clustering packages: {'yes' if clustering else 'NO - umap checks will be skipped'}")
    timings = {}
    t_all = time.time()
    if not args.checks_only:
        if not args.reuse:
            clean(root)
            clean_stereoseq(root)
        log_path.write_text("", encoding="utf-8")
        for cmd in demo_commands(root) + stereoseq_commands(root):
            label = Path(cmd[0]).stem if cmd[1].startswith("--") else f"{cmd[1]} {cmd[3]}"
            label += " (sample)" if "--umap-sample" in cmd else ""
            label += " (user clusters)" if "--user-layer" in cmd else ""
            label += f" ({cmd[cmd.index('--gene-files') + 1]} gene files)" if "--gene-files" in cmd else ""
            timings[label] = run(cmd, log_path)
            if cmd[0] == "examples/make_demo.py":
                write_user_files(root)
            if cmd[0] == "examples/make_xenium_demo.py":
                write_xenium_user_files(root)
                write_xenium_inputs(root)
        write_float_fixtures(root)
        for cmd in float_commands(root):                 # float32 / uint16 OME-TIFF into the 16-bit cache
            timings[f"{cmd[1]} {cmd[3]}"] = run(cmd, log_path)
        if not args.reuse:
            clean_imaging(root)
        for cmd in imaging_commands(root):               # CosMx / MERSCOPE fixtures
            timings[Path(cmd[0]).stem if cmd[1].startswith("--") else f"{cmd[1]} {cmd[3]}"] = run(cmd, log_path)
        if not args.reuse:
            clean_stdformats(root)
        for cmd in stdformats_commands(root):            # OME-Zarr, spatial AnnData, SpatialData
            timings[Path(cmd[0]).stem if cmd[1].startswith("--") else f"{cmd[1]} {cmd[3]}"] = run(cmd, log_path)
        if not args.reuse:
            clean_outlines(root)
        for i, cmd in enumerate(outline_commands(root)):  # outlines from cells.zarr.zip and from a lab's own files
            timings[f"{cmd[1]} {cmd[3]} (outlines {i + 1})"] = run(cmd, log_path)
        clean_detect(root)                               # prepare.py auto: always fresh, never skipped by --reuse
        make_detect_cases(root)
        for i, cmd in enumerate(detect_commands(root)):
            timings[f"auto {cmd[3]}" + (" (--kind)" if "--kind" in cmd else "")] = run(cmd, log_path)
        t_x = time.time()                                # zips read with every member copied out, and refused zips
        run_xenium_inputs(root, log_path)
        timings["xenium inputs (zip runs)"] = time.time() - t_x
        t_b = time.time()                                # prepare.py batch: always fresh, never skipped by --reuse
        run_batch(root, log_path)
        timings["batch (3 runs)"] = time.time() - t_b
        t1 = time.time()
        rna_demo = write_rnaseq_demo(root, log_path)
        timings["rnaseq demo"] = time.time() - t1
    else:
        rna_demo = {}

    c = Checker()
    t0 = time.time()
    print("checks: demo_visiumhd")
    summary = check_visium(c, root, clustering)
    print("checks: demo_multiplex")
    summary[GRID_KEY] = check_multiplex(c, root, clustering)
    print("checks: demo_tiled (tile seams)")
    summary["tiles"] = check_tiled(c, root)
    print("checks: float32 / uint16 OME-TIFF into the 16-bit cache (pyramid._float_scale)")
    summary["floatImage"] = check_float_image(c, root)
    print("checks: the lab's own clusters / UMAP (userclusters.py)")
    summary["userClusters"] = check_user(c, root, clustering, log_path)
    print("checks: demo_xenium (cells layer + molecules), demo_xenium_pq, demo_xenium_nt")
    summary["xenium"] = check_xenium(c, root, clustering)
    print("checks: Xenium inputs (no stain images, zips, refused zips)")
    summary["xeniumInputs"] = check_xenium_inputs(c, root)
    print("checks: the lab's own clusters on the demo_xenium cells (userclusters.py)")
    summary["xeniumUserClusters"] = check_xenium_user(c, root, clustering)
    print("checks: the optional barcode / cell id lists (layerids.py)")
    summary["ids"] = check_ids(c, root)
    print("checks: demo_geomx (GeoMx ROI layer, two scans), demo_geomx_grid (no image)")
    summary["geomx"] = check_geomx(c, root)
    print("checks: demo_cosmx*, demo_merscope* (FOV mosaic / mosaics, cells layer, molecules, density background)")
    summary["imaging"] = check_imaging(c, root, clustering, log_path)
    print("checks: standard formats (spatial AnnData, OME-Zarr, SpatialData)")
    summary["stdformats"] = check_stdformats(c, root, clustering, log_path)
    print("checks: demo_stereoseq (bin layers + cells), _gem, _lv, _mask")
    summary["stereoseq"] = check_stereoseq(c, root, clustering, log_path)
    print("checks: automatic detection (detect.py, prepare.py detect / auto)")
    summary["detect"] = check_detect(c, root, log_path)
    print("checks: batch export (batch.py, prepare.py batch)")
    summary["batch"] = check_batch(c, root, log_path)
    print("checks: the notebook")
    check_notebook(c)
    check_notebook_finder(c, root)
    print("checks: cell outlines (every platform's boundaries, traced label masks, a lab's own file)")
    summary["outlines"] = check_outlines(c, root, log_path)
    print("checks: gene packs (genepack.py): every export, loose / both, Range over HTTP, repack, deid")
    summary["genePacks"] = check_packs(c, root, log_path, rna_demo)
    print("checks: marker lists and name cleanup")
    check_marker_lists(c)
    print("checks: the page (index.html, viewer.js tile path, diag.js)")
    check_viewer_page(c, root, log_path)
    print("checks: BLAS / OpenMP thread cap")
    check_threadcap(c, log_path)
    print("checks: dataset picker, JavaScript syntax")
    check_picker(c)
    print("checks: tabbed layout (tabs, sections, ids)")
    check_layout(c)
    print("checks: sharing (the link that restores a view, linked windows)")
    check_sharing(c, log_path)
    print("checks: your own results (ownres.js: a reader's CSV / TSV matched in the browser)")
    check_ownres(c, root, log_path)
    print("checks: encrypted link files (links.py, linkPath): deprecated, still working")
    check_link_files(c)
    print("checks: offline use (libs.js, vendor.py, sw.js, offline.js)")
    check_offline(c, root, log_path)
    print("checks: publishable")
    check_publishable(c, root, log_path)
    timings["checks"] = time.time() - t0

    print("\ntimings: " + ", ".join(f"{k} {v:.0f}s" for k, v in timings.items())
          + f"; total {time.time() - t_all:.0f}s")
    (root / "smoke_test_summary.json").write_text(json.dumps(
        {"passed": not c.failed, "checks": c.n, "failed": c.failed, "clustering": clustering,
         "timingsSec": {k: round(v, 1) for k, v in timings.items()}, "layers": summary}, indent=1),
        encoding="utf-8")
    if c.failed:
        print(f"\nSMOKE TEST FAILED: {len(c.failed)} of {c.n} checks")
        for f in c.failed:
            print("  - " + f)
        return 1
    print(f"\nSMOKE TEST PASSED: {c.n} checks" + ("" if clustering else " (umap checks skipped)"))
    print(f"export: {(root / 'export').resolve()}")
    return 0


def test_smoke():
    """pytest entry point."""
    assert main([]) == 0


if __name__ == "__main__":
    sys.exit(main())
