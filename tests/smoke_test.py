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
Last, the page itself: every local script and
stylesheet index.html names exists and carries a ?v= string, diag.js (the ?diag=1 panel)
loads before viewer.js, viewer.js still has its three tile decoders in order and reports
tile and graphics failures, and, when node is installed, tests/diag_test.js passes (the
panel's token-stripping). The browser half is tests/phone_browser_test.js (Playwright).

Without the optional clustering packages (requirements-clustering.txt) the export
skips umap/ and so do these checks, with a message saying so. The earlier outputs
of these two demo datasets are removed first so nothing stale is checked; pass
--reuse to keep them. Exit code 0 = all checks passed, 1 = a check failed,
2 = a command failed. Needs only the standard library and numpy to run the checks.
"""
from __future__ import annotations

import argparse
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
# synthetic Stereo-seq folder (examples/make_stereoseq_demo.py): raw GEF + cell bin GEF / tissue GEM /
# older tissue GEF with stored bin10 + bin50 / shuffled whole-chip GEM cut with the mask image
ST, STG, STL, STM = "demo_stereoseq", "demo_stereoseq_gem", "demo_stereoseq_lv", "demo_stereoseq_mask"
ST_DS = (ST, STG, STL, STM)
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
         "--user-layer", "008um", "--user-obs", f"{r}/user/clusters.csv", "--user-umap", f"{r}/user/umap.csv"],
        ["examples/make_xenium_demo.py", "--out", raw],
        [prep, "xenium", "--id", XE, "--src", f"{raw}/xenium", "--out", cache, "--label", "Synthetic Xenium"],
        [prep, "export", "--id", XE, "--out", cache, "--dest", export, "--markers", "examples/demo_xenium_markers.json"],
        [prep, "xenium", "--id", XP, "--src", f"{raw}/xenium_pq", "--out", cache],
        [prep, "export", "--id", XP, "--out", cache, "--dest", export, "--grid-um", "0"],
        [prep, "xenium", "--id", XN, "--src", f"{raw}/xenium_nt", "--out", cache],
        [prep, "export", "--id", XN, "--out", cache, "--dest", export, "--no-umap"],
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


def clean(root: Path):
    """Remove only what this test itself writes, so nothing stale gets checked."""
    for p in [root / "raw", root / "user", root / "export_sample" / MX, root / "export_user" / VIS,
              root / "export_partial" / VIS] + [root / sub / ds for sub in ("cache", "export") for ds in (VIS, MX, TL) + XEN_DS]:
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

    # genes.json + gene/<NAME>.bin: [x, y, v] float32, values summing to the gene total
    genes = read_json(d / "genes.json")
    c.eq(len(genes), L["nGenes"], f"{key}: genes.json length")
    c.eq(len(genes), 300, f"{key}: demo gene count")
    row_of = {p.tobytes(): i for i, p in enumerate(pos)}
    lib = np.zeros(n, np.float64)
    dense = {}
    n_files, bad = 0, []
    for g in genes:
        p = d / "gene" / f"{enc(g['n'])}.bin"
        if g["t"] <= 0:
            continue
        if not p.exists():
            bad.append(f"missing {p.name}")
            continue
        n_files += 1
        a = f32(p, 3)
        if abs(float(a[:, 2].sum()) - g["t"]) > 1e-3 * max(1.0, g["t"]):
            bad.append(f"{g['n']}: sum {a[:, 2].sum()} != total {g['t']}")
        rows = np.array([row_of.get(r.tobytes(), -1) for r in a[:, :2]])
        if (rows < 0).any():
            bad.append(f"{g['n']}: points not on positions.bin")
            continue
        np.add.at(lib, rows, a[:, 2])
        dense[g["n"]] = (rows, a[:, 2])
    c.ok(not bad, f"{key}: gene files: {bad[:5]}")
    c.eq(len(list((d / "gene").glob("*.bin"))), n_files, f"{key}: gene/*.bin count")
    print(f"    positions.bin {n:,} rows, {n_files} gene files, median {np.median(lib):.0f} UMI per bin")

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
        p = d / "gene" / f"{enc(x['n'])}.bin"
        if not p.exists():
            bad.append(f"missing {p.name}"); continue
        a = f32(p, 3)
        if len(a) != n or not np.array_equal(a[:, :2], pos):
            bad.append(f"{x['n']}: rows / xy differ from positions.bin"); continue
        if not (np.isfinite(a[:, 2]).all() and (a[:, 2] >= 0).all()):
            bad.append(f"{x['n']}: values NaN or < 0"); continue
        if abs(float(a[:, 2].mean()) - x["t"]) > 1e-3 * max(1.0, x["t"]):
            bad.append(f"{x['n']}: mean {a[:, 2].mean()} != t {x['t']}")
        V[x["n"]] = a[:, 2].astype(np.float64)
    c.ok(not bad, f"{key}: gene files: {bad[:5]}")
    c.eq(len(list((d / "gene").glob("*.bin"))), len(genes), f"{key}: gene/*.bin count")

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
                    bad.append(f"{p.name}: {nm} values differ from gene/{enc(nm)}.bin"); break
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
        plain = f32(root / "export" / MX / "expr" / GRID_KEY / "gene" / f"{enc(nm)}.bin", 3)[:, 2]
        fixed = f32(d / "gene" / f"{enc(nm)}.bin", 3)[:, 2]
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
        f = md / "gene" / f"{enc(g['n'])}.bin"
        if not g["c"]:
            if f.exists():
                bad.append(f"{g['n']}: file although c = 0")
            continue
        if not f.exists() or f.stat().st_size != 4 * g["c"]:
            bad.append(f"{g['n']}: {f.stat().st_size if f.exists() else 'no'} bytes, want {4 * g['c']}")
            continue
        packed[g["n"]] = np.fromfile(str(f), dtype="<u4")
    c.ok(not bad, f"{ds}: molecule files: {bad[:4]}")
    c.eq(len(list((md / "gene").glob("*.bin"))), len(packed), f"{ds}: molecules/gene/*.bin count")
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
        f = d / "gene" / f"{enc(g)}.bin"
        rows = np.flatnonzero(cnt[:, j] > 0)
        a = f32(f, 3) if f.exists() else np.zeros((0, 3), np.float32)
        if len(a) != len(rows) or not np.allclose(a[:, :2], pos[rows]) or not np.allclose(a[:, 2], cnt[rows, j]):
            bad.append(g)
        dense[g] = (rows, cnt[rows, j])
    c.ok(not bad, f"cells gene files differ from the count matrix: {bad[:4]}")
    c.ok(not any((d / "gene" / f"{enc(x)}.bin").exists() for x in tr["controls"]), "cells: a control feature has a gene file")
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
    c.ok(all(np.allclose(f32(out_n / "expr" / "cells" / "gene" / f"{enc(g)}.bin", 3), f32(d / "gene" / f"{enc(g)}.bin", 3),
                         atol=1e-3) for g in genes), f"{XN}: gene files = the csv / mtx import")
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
        f = d / "gene" / f"{enc(g)}.bin"
        a = f32(f, 3) if f.exists() else np.zeros((0, 3), np.float32)
        if len(a) != len(rows) or not np.array_equal(a[:, :2], pos_want[rows]) or not np.array_equal(a[:, 2], dense_want[rows, j]):
            bad.append(g)
        dense[g] = (rows, dense_want[rows, j])
    c.ok(not bad, f"{what}: gene files differ from the drawn counts: {bad[:4]}")
    c.eq(len(list((d / "gene").glob("*.bin"))), int((dense_want.sum(0) > 0).sum()), f"{what}: gene/*.bin count")
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
            fa, fb = a / "gene" / f"{enc(ga)}.bin", b / "gene" / f"{enc(gb)}.bin"
            ok = ok and fa.exists() == fb.exists() and (not fa.exists() or fa.read_bytes() == fb.read_bytes())
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



def check_marker_lists(c: Checker):
    """The built-in lists load and are usable, and the channel-name cleanup is the same in Python and JS."""
    sys.path.insert(0, str(APP))
    try:
        import celltypes
        import seams  # noqa: F401  (used below)
    finally:
        sys.path.pop(0)
    for name, prot in (("protein_generic.json", True), ("rna_generic.json", False), ("akoya_nerve.json", True),
                       ("nerve.json", False)):
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
    local = [r for r in refs if not re.match(r"https?:", r)]
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
    for f in ("viewer.js", "diag.js"):
        q = subprocess.run([node, "--check", str(APP / f)], cwd=str(REPO), stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        c.eq(q.returncode, 0, f"node --check {f}")


def check_publishable(c: Checker, root: Path, log_path: Path):
    export = root / "export"
    texts = [p for ds in (VIS, MX, TL) + XEN_DS for p in (export / ds).rglob("*")
             if p.is_file() and p.suffix in (".json", ".geojson", ".txt")]
    texts += [p for p in (root / "export_sample" / MX).rglob("*") if p.is_file() and p.suffix == ".json"]
    texts += [p for p in (root / "export_user" / VIS).rglob("*") if p.is_file() and p.suffix in (".json", ".geojson")]
    bc = [p for p in texts if BARCODE.search(p.read_text(encoding="utf-8", errors="ignore"))]
    c.ok(not bc, f"export text files contain barcodes: {bc[:3]}")
    lp = [p for p in texts if LOCAL_PATH.search(p.read_text(encoding="utf-8", errors="ignore"))]
    c.ok(not lp, f"export text files contain local paths: {lp[:3]}")
    args = (["web/SpatialViewer/deid.py", "check", str(export / VIS), str(export / MX), str(export / TL)] +
            [str(export / ds) for ds in XEN_DS])
    if (root / "export_sample" / MX).exists():
        args.append(str(root / "export_sample" / MX))
    for extra in (root / "export_user" / VIS, root / "export_partial" / VIS):
        if extra.exists():
            args.append(str(extra))
    print(f"$ python {' '.join(args)}", flush=True)
    p = subprocess.run([sys.executable, *args], cwd=str(REPO), stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                       text=True, encoding="utf-8", errors="replace")
    with open(log_path, "a", encoding="utf-8") as lg:
        lg.write(p.stdout)
    print("    " + p.stdout.strip().splitlines()[-1] if p.stdout.strip() else "    (no output)")
    c.eq(p.returncode, 0, "deid.py check on the export (exit code)")


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
            timings[label] = run(cmd, log_path)
            if cmd[0] == "examples/make_demo.py":
                write_user_files(root)

    c = Checker()
    t0 = time.time()
    print("checks: demo_visiumhd")
    summary = check_visium(c, root, clustering)
    print("checks: demo_multiplex")
    summary[GRID_KEY] = check_multiplex(c, root, clustering)
    print("checks: demo_tiled (tile seams)")
    summary["tiles"] = check_tiled(c, root)
    print("checks: the lab's own clusters / UMAP (userclusters.py)")
    summary["userClusters"] = check_user(c, root, clustering, log_path)
    print("checks: demo_xenium (cells layer + molecules), demo_xenium_pq, demo_xenium_nt")
    summary["xenium"] = check_xenium(c, root, clustering)
    print("checks: demo_stereoseq (bin layers + cells), _gem, _lv, _mask")
    summary["stereoseq"] = check_stereoseq(c, root, clustering, log_path)
    print("checks: marker lists and name cleanup")
    check_marker_lists(c)
    print("checks: the page (index.html, viewer.js tile path, diag.js)")
    check_viewer_page(c, root, log_path)
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
