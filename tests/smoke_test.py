#!/usr/bin/env python
r"""Smoke test: the page's own files, the public example exports, and (with their raw files) the pipeline that makes them.

    python tests/smoke_test.py                                     (from the repo root)
    python tests/smoke_test.py --raw-data <dir> [--root <dir>]     (also prepare + export from the raw files)
    python tests/smoke_test.py --exports <folder>                  (check exports made earlier, run nothing)
    python -m pytest tests/smoke_test.py                           (same test, if pytest is installed)

The synthetic demo data and its generators are gone (19 Nov 2026); the test runs on the public datasets of
REPRODUCE.md section 9 instead. Its three parts:

1. Always (no data needed): the built-in marker lists and the channel-name cleanup (celltypes.py = lasso_worker.js),
   the seam check on profiles with a known answer, the page's files (index.html, viewer.js tile path, diag.js, every
   .js parses), the BLAS / OpenMP thread cap (threadcap.py) and the PCA on one thread, the dataset picker, the tabbed
   layout, the cell-outline controls, sharing (viewlink.js), your own results (ownres.js), the pinned libraries
   (libs.js, vendor.py, sw.js removing the old offline use) and the notebook.

2. The pipeline, only with the raw files: --raw-data <dir> or SODA_RAW_DATA=<dir>, laid out as REPRODUCE.md section 9
   downloads them (<dir>/xenium_pub/data.zarr, <dir>/visiumhd_pub/..., one folder per dataset id). For each public
   dataset whose input is there (or those of --ids), REPRODUCE.md's `prepare.py <kind> ...` and
   `prepare.py export --id <id> --out <root>/cache --dest <root>/export --fmt jpg` are run, then the export is checked
   (part 3) and passed through `deid.py check`. A dataset whose raw files are not there is SKIPPED with a line saying
   so. Without a raw-data folder this part is SKIPPED with one line: nothing is prepared, nothing is passed.
   The real datasets take from minutes (merscope_2) to hours (visiumhd_pub); --root (default: soda-smoke in the
   system temp folder) needs room for their caches and exports.

3. The exports: every public export in <root>/export after part 2, or in --exports <folder> (a local copy, read
   only; SODA_TEST_DATA when that is a folder): manifest fields, every image tile there and of its format, each
   layer's positions.bin (8 bytes x nBins, finite), genes.json (nGenes), the gene packs against their listings
   (genepack.check_dataset) and a sample of genes (values summing to the gene's total, points on positions.bin),
   umap/ (rows, clusterings, sizes), hvg/ chunks, binqc/ and the barcode lists' row counts, the molecules pack, no
   barcode or local path in a text file, and `deid.py check`. Without exports this part is SKIPPED with one line.
   Only the seven public ids are ever read (xenium_pub, visiumhd_pub, phenocycler_pub, merscope_2, cosmx_1,
   stereoseq_1, geomx_1), whatever else a folder holds.

The browser tests (Playwright, when node is installed; each SKIPPED with a line saying so without it) run on the
shared test data, tests/test_data.js: SODA_TEST_DATA (default the hosted public exports,
https://fivestore.blob.core.windows.net/soda/examples/) and SODA_TEST_PRIVATE_DATA for the restricted three
(cosmx_1, stereoseq_1, geomx_1). When part 2 made all four hosted exports, they run on those instead
(SODA_TEST_DATA=<root>/export). They are: tests/diag_test.js, compat_test.js (old exports and the current packed
ones), browsers_test.js, drawregion_test.js, stdformats_browser_test.js, outlines_browser_test.js,
controls_browser_test.js, ownres_browser_test.js and viewlink_test.js.

What the synthetic data checked and nothing checks now (no public dataset has a known answer): the export against
the generator's ground truth (counts per bin, clusters against drawn regions, suggested cell-type names), the
Stereo-seq / GeoMx / CosMx / MERSCOPE / Xenium readers' input variants (GEM, older GEF, mask, zips, no image,
parquet / h5 only), the float32 image cast, tile-seam correction, the lab's own clusters and segmentation (--user-obs,
--cells), standard formats other than SpatialData (spatial AnnData, OME-Zarr), cell outlines from every file form,
automatic detection and batch export, the notebook's dataset finder, and gene packs against loose files.

Without the optional clustering packages (requirements-clustering.txt) an export has no umap/ and its checks say so.
Exit code 0 = all checks passed, 1 = a check failed, 2 = a command failed.
"""
from __future__ import annotations

import argparse
import gzip
import importlib
import importlib.util
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from urllib.parse import quote

import numpy as np

REPO = Path(__file__).resolve().parents[1]
APP = REPO / "web" / "SpatialViewer"
# the public example datasets (REPRODUCE.md section 9); tests/test_data.js has the same lists
HOSTED = ("xenium_pub", "visiumhd_pub", "phenocycler_pub", "merscope_2")       # on the public blob
RESTRICTED = ("cosmx_1", "stereoseq_1", "geomx_1")                            # public, but not to be passed on
PUBLIC = HOSTED + RESTRICTED
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



# ------------------------------------------------------------------ small helpers
def read_json(p: Path):
    return json.loads(p.read_text(encoding="utf-8"))



def enc(name: str) -> str:
    return quote(name, safe="-_.!~*'()")              # = JavaScript encodeURIComponent



def f32(p: Path, cols: int) -> np.ndarray:
    return np.fromfile(str(p), dtype="<f4").reshape(-1, cols)



def app_module(name: str):
    sys.path.insert(0, str(APP))
    try:
        return importlib.import_module(name)
    finally:
        sys.path.pop(0)



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



# ------------------------------------------------------------------ the pipeline on the raw files (REPRODUCE.md)
GX_DIR = "mu_dev_E11_workflow_and_count_files"


def raw_commands(raw: Path, cache: Path) -> dict[str, tuple[Path, list[str]]]:
    """Each public id: (the input that must be there, REPRODUCE.md's prepare command)."""
    P = "web/SpatialViewer/prepare.py"
    d = lambda *p: raw.joinpath(*p)
    s = lambda p: str(p)
    gx = d("geomx_1", GX_DIR)
    return {
        "xenium_pub": (d("xenium_pub", "data.zarr"),
                       [P, "spatialdata", "--src", s(d("xenium_pub", "data.zarr")), "--id", "xenium_pub", "--out", s(cache),
                        "--pixel-size", "0.2125", "--label", "Breast cancer, Xenium (10x public data)"]),
        "visiumhd_pub": (d("visiumhd_pub", "Visium_HD_Human_Colon_Cancer_tissue_image.btf"),
                         [P, "visium", "--id", "visiumhd_pub", "--image", s(d("visiumhd_pub", "Visium_HD_Human_Colon_Cancer_tissue_image.btf")),
                          "--outs", s(d("visiumhd_pub")), "--out", s(cache), "--label", "Colorectal cancer, Visium HD (10x public data)"]),
        "phenocycler_pub": (d("phenocycler_pub", "D438-RLL-16B4-17_Scan1.qptiff"),
                            [P, "tiff", "--src", s(d("phenocycler_pub", "D438-RLL-16B4-17_Scan1.qptiff")), "--id", "phenocycler_pub",
                             "--out", s(cache), "--label", "Lung, PhenoCycler (HuBMAP public data)"]),
        "merscope_2": (d("merscope_2", "AAN4"), [P, "merscope", "--src", s(d("merscope_2", "AAN4")), "--id", "merscope_2", "--out", s(cache)]),
        "cosmx_1": (d("cosmx_1", "Lung9_Rep1_exprMat_file.csv"), [P, "cosmx", "--src", s(d("cosmx_1")), "--id", "cosmx_1", "--out", s(cache)]),
        "stereoseq_1": (d("stereoseq_1", "C04042E3.raw.gef"),
                        [P, "stereoseq", "--image", s(d("stereoseq_1", "C04042E3_HE_regist.tif")), "--gef", s(d("stereoseq_1", "C04042E3.raw.gef")),
                         "--cellbin", s(d("stereoseq_1", "C04042E3.cellbin.gef")), "--id", "stereoseq_1", "--out", s(cache)]),
        "geomx_1": (gx / "workflow" / "dcc" / "E9-E11-E15_DCC_unzipped",
                    [P, "geomx", "--dcc", s(gx / "workflow" / "dcc" / "E9-E11-E15_DCC_unzipped"),
                     "--pkc", s(gx / "workflow" / "pkc" / "Mm_R_NGS_WTA_v1.0_unzipped" / "Mm_R_NGS_WTA_v1.0.pkc"),
                     "--annotation", s(gx / "workflow" / "readout_package" / "E11E15readoutpackage" / "SZ021_20211213T1749_LabWorksheet.txt"),
                     "--annotation", s(gx / "count" / "Export1_InitialDataset.xlsx"),
                     "--image", s(d("geomx_1", "mu_dev_E11_image_files", "mu_dev_E11_004.ome.tiff")),
                     "--image", s(d("geomx_1", "mu_dev_E11_image_files", "mu_dev_E11_005.ome.tiff")),
                     "--id", "geomx_1", "--out", s(cache)]),
    }


def run_pipeline(raw: Path, root: Path, ids: list[str], log_path: Path, timings: dict) -> list[str]:
    """prepare + export for each id whose input is under raw/; the ids exported."""
    cache, export = root / "cache", root / "export"
    cmds = raw_commands(raw, cache)
    done = []
    for ds in ids:
        need, prep = cmds[ds]
        if not need.exists():
            print(f"  SKIPPED {ds}: no {need.relative_to(raw).as_posix()} under the raw-data folder")
            continue
        for p in (cache / ds, export / ds):            # what this test wrote last time: nothing stale is checked
            if p.exists():
                print(f"  removing earlier output {p}", flush=True)
                shutil.rmtree(p)
        timings[f"prepare {ds}"] = run(prep, log_path)
        timings[f"export {ds}"] = run(["web/SpatialViewer/prepare.py", "export", "--id", ds, "--out", str(cache),
                                       "--dest", str(export), "--fmt", "jpg"], log_path)
        done.append(ds)
    return done


# ------------------------------------------------------------------ an export, without a known answer
def _sample(n: int, k: int, seed: int = 0) -> np.ndarray:
    return np.unique(np.r_[0, n - 1, np.random.default_rng(seed).integers(0, n, k)]) if n else np.zeros(0, np.int64)


def check_layer_files(c: Checker, out: Path, m: dict, L: dict, clustering: bool) -> dict:
    """One expression layer: positions, genes, a sample of genes from the packs, umap/, hvg/, binqc/, ids."""
    gp = app_module("genepack")
    key, ds = L["key"], m["id"]
    d = out / "expr" / key
    what = f"{ds} {key}"
    for k in ("label", "nBins", "nGenes", "binSizeUm", "micronsPerPixel", "spotDiameterPx"):
        c.ok(k in L, f"{what}: manifest layer has no {k}")
    n = int(L.get("nBins", -1))
    c.eq((d / "positions.bin").stat().st_size if (d / "positions.bin").exists() else -1, 8 * n, f"{what}: positions.bin bytes (8 x nBins)")
    pos = f32(d / "positions.bin", 2)
    c.ok(np.isfinite(pos).all(), f"{what}: positions.bin has NaN or inf")
    W, H = m["image"]["width"], m["image"]["height"]
    inside = float(((pos[:, 0] >= 0) & (pos[:, 0] <= W) & (pos[:, 1] >= 0) & (pos[:, 1] <= H)).mean()) if len(pos) else 0.0
    c.ok(inside > 0.99, f"{what}: {inside:.1%} of positions inside the {W} x {H} image")
    genes = read_json(d / "genes.json")
    c.eq(len(genes), L.get("nGenes"), f"{what}: genes.json length vs nGenes")
    # the packs: structure against the listing (genepack.check_dataset covers every group of the dataset), then
    # a sample of genes, largest included: values sum to the gene's total, every point on a row of positions.bin
    g = next((x for x in gp.groups(m) if x["dir"] == f"expr/{key}"), None)
    c.ok(g is not None and g["pack"], f"{what}: a packed gene group")
    row_of = {p.tobytes(): i for i, p in enumerate(pos)}
    listing = gp.read_listing(out, g) if g else []
    have = [e for e in listing if e.get("b", 0) > 0 or "s" in e]
    pick = {e["n"]: e for e in [have[i] for i in _sample(len(have), 12)]} if have else {}
    if have:
        big = max(have, key=lambda e: e.get("b", 0))
        pick[big["n"]] = big
    bad = []
    for e in pick.values():
        b = gp.gene_bytes(out, g, e)
        if b is None or len(b) % 12:
            bad.append(f"{e['n']}: {None if b is None else len(b)} bytes")
            continue
        a = np.frombuffer(b, dtype="<f4").reshape(-1, 3)
        t = e.get("t") if L.get("kind") != "protein" else None    # a protein grid's t is not a sum of counts
        if t is not None and abs(float(a[:, 2].astype(np.float64).sum()) - t) > 1e-3 * max(1.0, abs(t)):
            bad.append(f"{e['n']}: values sum {a[:, 2].sum():.4g}, total {t}")
        if any(r.tobytes() not in row_of for r in a[:, :2]):
            bad.append(f"{e['n']}: points not on positions.bin")
    c.ok(not bad, f"{what}: genes read from the packs ({len(pick)} checked): {bad[:4]}")
    res = {"rows": n, "genes": len(genes), "genesChecked": len(pick)}
    # umap/ (sectionumap.py; GeoMx: its own embedding)
    ud = d / "umap"
    if L.get("umap"):
        idx = read_json(ud / "index.json")
        c.eq(idx.get("nBins", n), n, f"{what}: umap nBins vs positions.bin rows")
        co = f32(ud / "coords.bin", 2) if (ud / "coords.bin").exists() else None
        if c.ok(co is not None and len(co) == n, f"{what}: umap/coords.bin, one row per bin"):
            emb = ~np.isnan(co[:, 0])
            c.eq(int(emb.sum()), idx.get("nEmbedded"), f"{what}: embedded rows vs nEmbedded")
            for cl in idx.get("clusterings", []):
                p = ud / cl["file"]
                lab = np.fromfile(str(p), dtype="<u2") if p.exists() else np.zeros(0, np.uint16)
                c.eq(len(lab), n, f"{what}: {cl['file']} rows")
                if len(lab) == n and cl.get("source") == "user":
                    # a grouping from the data's own annotations: every row it names, embedded or not
                    got = lab[lab != NO_CLUSTER].astype(np.int64)
                    c.ok(len(got) == cl.get("nAssigned", len(got)) and (got < cl["n"]).all(),
                         f"{what}: {cl['file']}: {len(got)} rows assigned (nAssigned {cl.get('nAssigned')}), labels below n {cl['n']}")
                elif len(lab) == n:
                    c.ok(((lab == NO_CLUSTER) == ~emb).all(), f"{what}: {cl['file']}: 0xFFFF rows = rows not embedded")
                    if "sizes" in cl:
                        c.eq(np.bincount(lab[emb].astype(np.int64), minlength=cl["n"]).tolist(), cl["sizes"], f"{what}: {cl['file']}: sizes")
        raw = (ud / "index.json").read_text(encoding="utf-8")
        c.ok("NaN" not in raw and "Infinity" not in raw, f"{what}: umap/index.json has no NaN / Infinity (a browser cannot parse it)")
        res["nEmbedded"] = idx.get("nEmbedded")
    elif not clustering:
        print(f"    {what}: no umap/ (clustering packages not installed)")
    # hvg/: the chunks the index names, each of its size
    if L.get("hvg"):
        hi = read_json(d / "hvg" / "index.json")
        miss = [f"{ch['tx']}_{ch['ty']}" for ch in hi.get("chunks", [])
                if not (d / "hvg" / f"{ch['tx']}_{ch['ty']}.bin").exists() or (d / "hvg" / f"{ch['tx']}_{ch['ty']}.bin").stat().st_size != ch["bytes"]]
        c.ok(hi.get("chunks") and not miss, f"{what}: hvg chunks there with their sizes ({len(hi.get('chunks', []))}; wrong: {miss[:3]})")
    # binqc/: one value per row
    q = L.get("binQc")
    if q:
        for mt in q.get("metrics", []):
            f = d / "binqc" / mt["file"]
            size = {"uint16": 2, "uint8": 1}.get(mt.get("dtype"), 4)
            c.eq(f.stat().st_size if f.exists() else -1, size * n, f"{what}: binqc/{mt['file']} bytes ({size} x rows)")
    # the barcode / cell id lists (layerids.py): one line per row
    if L.get("ids"):
        ids = [x for x in gzip.decompress((d / "ids.txt.gz").read_bytes()).decode("utf-8").split("\n") if x]
        c.eq(len(ids), n, f"{what}: ids.txt.gz lines vs rows")
    return res


def check_export(c: Checker, out: Path, clustering: bool) -> dict:
    """One public export: the manifest, the tiles, each layer, the molecules, the packs, no barcode or local path."""
    gp = app_module("genepack")
    m = read_json(out / "manifest.json")
    ds = out.name
    ext = m.get("tileFormat")
    c.ok(ext in ("png", "jpg"), f"{ds}: tileFormat {ext!r}")
    check_manifest_common(c, m, ds, ext)
    check_tiles(c, out, m, ext)
    probs = gp.check_dataset(out)
    c.ok(not probs, f"{ds}: gene packs against their listings (genepack.check_dataset): {probs[:3]}")
    res = {"kind": m.get("kind"), "layers": {}}
    for L in m.get("expr", []):
        print(f"  expr {L['key']}:")
        res["layers"][L["key"]] = check_layer_files(c, out, m, L, clustering)
        print(f"    {res['layers'][L['key']]}")
    mol = m.get("molecules")
    if mol:
        md = out / (mol.get("dir") or "molecules")
        idx = read_json(md / "index.json") if (md / "index.json").exists() else {}
        c.ok(bool(idx.get("genes")), f"{ds}: molecules/index.json lists genes")
        tot = sum(x.get("c", 0) for x in idx.get("genes", []))
        c.ok(tot == mol.get("n", tot), f"{ds}: molecules: {tot} in index.json, manifest n {mol.get('n')}")
        res["molecules"] = tot
    texts = [p for p in out.rglob("*") if p.is_file() and p.suffix in (".json", ".geojson", ".txt")]
    bc = [p.relative_to(out).as_posix() for p in texts if BARCODE.search(p.read_text(encoding="utf-8", errors="ignore"))]
    c.ok(not bc, f"{ds}: export text files contain barcodes: {bc[:3]}")
    lp = [p.relative_to(out).as_posix() for p in texts if LOCAL_PATH.search(p.read_text(encoding="utf-8", errors="ignore"))]
    c.ok(not lp, f"{ds}: export text files contain local paths: {lp[:3]}")
    return res


def check_exports(c: Checker, folder: Path, ids: list[str], clustering: bool, log_path: Path) -> dict:
    out = {}
    for ds in ids:
        print(f"checks: {ds} ({folder / ds})")
        out[ds] = check_export(c, folder / ds, clustering)
    args = ["web/SpatialViewer/deid.py", "check", *[str(folder / ds) for ds in ids]]
    print(f"$ python {' '.join(args)}", flush=True)
    p = subprocess.run([sys.executable, *args], cwd=str(REPO), stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                       text=True, encoding="utf-8", errors="replace")
    with open(log_path, "a", encoding="utf-8") as lg:
        lg.write(p.stdout)
    print("    " + (p.stdout.strip().splitlines() or ["(no output)"])[-1])
    c.eq(p.returncode, 0, "deid.py check on the exports (exit code)")
    return out


# ------------------------------------------------------------------ the browser tests (tests/test_data.js data)
def node_test(c: Checker, log_path: Path, script: str, args: list[str], what: str, tail: int = 1) -> str:
    """Run one node test; its SKIPPED lines and last line(s) shown; exit code 0 wanted. '' when node is missing."""
    node = shutil.which("node")
    if not node:
        print(f"    node not found: tests/{script} skipped")
        return ""
    p = subprocess.run([node, str(REPO / "tests" / script), *args], cwd=str(REPO),
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace")
    with open(log_path, "a", encoding="utf-8") as lg:
        lg.write(f"\n$ node tests/{script} {' '.join(args)}\n{p.stdout}")
    print_skipped(p.stdout)
    for ln in (p.stdout.strip().splitlines() or ["(no output)"])[-tail:]:
        print("    " + ln)
    c.eq(p.returncode, 0, f"node tests/{script}: {what} (exit code)")
    return p.stdout


def check_format_pages(c: Checker, root: Path, log_path: Path):
    """Cell outlines in the page (outlines.js) and, in the browser, the instrument / standard-format exports, the
    outlines and what draws a subset of the cells with them, and one control for each job."""
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
    node_test(c, log_path, "stdformats_browser_test.js", ["--shots", str(root / "stdformats_shots"), "--browsers", "chromium,firefox"],
              "the SpatialData and instrument exports open in Chromium and Firefox")
    node_test(c, log_path, "outlines_browser_test.js", ["--shots", str(root / "outlines_shots"), "--browsers", "chromium,webkit,firefox"],
              "outlines, and what draws a subset of the cells with them, in Chromium, WebKit and Firefox")
    node_test(c, log_path, "controls_browser_test.js", ["--shots", str(root / "controls_shots"), "--browsers", "chromium,webkit,firefox"],
              "one control for each job, old links with the image off, in Chromium, WebKit and Firefox")


# ------------------------------------------------------------------ the page and the code (no data)
NB_HELPER = REPO / "notebooks" / "soda_notebook.py"



def notebook_helper():
    """notebooks/soda_notebook.py, the notebook's code, imported from notebooks/ as the notebook's first cell does."""
    folder = str(NB_HELPER.parent)
    sys.path.insert(0, folder)
    try:
        import soda_notebook
    finally:
        sys.path.remove(folder)
    return soda_notebook



def check_notebook(c: Checker):
    """notebooks/SODA_pipeline.ipynb and its code, notebooks/soda_notebook.py: valid, no saved outputs, code that
    compiles, the helper found from notebooks/ and from the repository folder, and its dataset finding is batch.py's."""
    import ast
    import inspect
    p = REPO / "notebooks" / "SODA_pipeline.ipynb"
    try:
        nb = json.loads(p.read_text(encoding="utf-8"))
    except ValueError as e:
        c.ok(False, f"notebook: valid JSON ({e})")
        return
    c.ok(nb.get("nbformat") == 4 and isinstance(nb.get("cells"), list), "notebook: nbformat 4")
    code = [x for x in nb["cells"] if x["cell_type"] == "code"]
    c.ok(len(code) >= 10 and all(not x.get("outputs") and x.get("execution_count") is None for x in code),
         "notebook: no saved outputs or execution counts")
    bad = []
    for x in code:
        try:
            compile("".join(x["source"]), x.get("id", "cell"), "exec")
        except SyntaxError as e:
            bad.append(f"{x.get('id')}: {e}")
    try:
        helper = NB_HELPER.read_text(encoding="utf-8")
        compile(helper, str(NB_HELPER), "exec")
    except (OSError, SyntaxError) as e:
        bad.append(f"soda_notebook.py: {e}")
        helper = ""
    c.ok(not bad, f"notebook: every code cell and soda_notebook.py compile {bad}")
    src = "\n".join("".join(x["source"]) for x in code)
    # the first cell finds the helper whether Jupyter's working folder is notebooks/ or the repository folder
    first = "".join(code[0]["source"]) if code else ""
    found = []
    for cwd in (REPO / "notebooks", REPO):
        r = subprocess.run([sys.executable, "-c", first.replace("soda.setup()", "print(soda.__file__)")], cwd=str(cwd),
                           capture_output=True, text=True, encoding="utf-8", errors="replace")
        found.append(bool(r.returncode == 0 and r.stdout.strip() and Path(r.stdout.strip()).resolve() == NB_HELPER.resolve()))
    c.ok("import soda_notebook as soda" in first and "soda.setup()" in first and all(found),
         f"notebook: the first cell imports soda_notebook.py from notebooks/ and from the repository folder ({found})")
    # Section 2: batch.find_datasets (no copy of its own), detect.option_lines, and the OPTIONS dictionary in the notebook
    c.ok("import detect" in helper and "import batch" in helper and "batch.find_datasets(" in helper
         and "detect.option_lines(" in helper and "OPTIONS = {" in src and "soda.check_inputs(OPTIONS)" in src
         and not any(s in src + helper for s in ("def check_dataset", "is_kind_folder")),
         "notebook: Section 2 finds the datasets with batch.find_datasets (no copy of its own) and has the OPTIONS dictionary")
    c.ok(not any(s in src + helper for s in ("cell_by_gene", 'glob("*.qptiff', "morphology_focus", "_regist")),
         "notebook: no per-kind file finding of its own, in the cells or soda_notebook.py (detect.py has it)")
    # the same scripts as the command lines, and the printed lines (and safety gates) the docs and the paper show
    kept = ["prepare.py", "catalog.py", "deid.py", "static_serve.py", "[deid] clean", "OK to upload.",
            "'READY' if not d['problems'] else 'FIX  '", "THE LINK TO SHARE:", "OK: a page on {origin} can read your folder.",
            'Answer the question first: set IS_THIS_DATA_PUBLIC = "yes" or "no"', "Run Section 5 (the de-identification check) first.",
            "Section 5 found sample IDs.", "removed the '?...' part of CLOUD_URL",
            "anyone holding it can read the data", "keep the key file outside the work folder and the repository"]
    gone = [s for s in kept if s not in helper]
    c.ok(not gone, f"notebook: soda_notebook.py runs the same scripts and keeps the printed lines and gates ({gone})")
    c.ok('IS_THIS_DATA_PUBLIC = ""' in src and "soda.share(IS_THIS_DATA_PUBLIC)" in src,
         "notebook: Section 7 asks is this data public? (empty until answered)")
    # the advanced settings cell names every setting of soda.advanced(), with its default
    try:
        soda = notebook_helper()
        adv = next("".join(x["source"]) for x in code if "soda.advanced(" in "".join(x["source"]))
        call = next(n for n in ast.walk(ast.parse(adv)) if isinstance(n, ast.Call) and getattr(n.func, "attr", "") == "advanced")
        params = inspect.signature(soda.advanced).parameters
        defaults = {k: (v.default if v.default is not None else {"visium": "", "proteomics": ""}) for k, v in params.items()}
        given = {}
        for k in call.keywords:
            try:
                given[k.arg] = ast.literal_eval(k.value)
            except ValueError:
                given[k.arg] = defaults.get(k.arg) if ast.unparse(k.value) == "soda.DEFAULT_KEY_FILE" else ast.unparse(k.value)
        c.eq(given, defaults, "notebook: the advanced settings cell lists every setting with soda.advanced()'s default")
    except Exception as e:  # noqa: BLE001
        c.ok(False, f"notebook: the advanced settings cell ({type(e).__name__}: {e})")
    text = "\n".join("".join(x["source"]) for x in nb["cells"] if x["cell_type"] == "markdown")
    heads = [ln for ln in text.splitlines() if ln.startswith("## ")]
    c.ok([h.split(".", 1)[0] for h in heads] == [f"## {i}" for i in range(8)],
         f"notebook: Sections 0 to 7 in order ({heads})")
    c.ok("own folder under `inputs/`" in text and "How sure" not in helper.split("HOW_SURE", 1)[0] and "`kind`" in text,
         "notebook: Section 1 says one folder per dataset, Section 2 documents the kind override")



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
    # vendor/ holds the libraries only after `prepare.py vendor` (git-ignored, version in the path): check_libs
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
    # site.js: the tracked default data folder (a site published straight from the repository has no config.js)
    sj = (APP / "site.js").read_text(encoding="utf-8")
    m = re.search(r'\bdata:\s*"([^"]*)"', sj)
    c.ok("window.SodaSite" in sj and m is not None and (m.group(1) == "" or m.group(1).startswith("https://")),
         f"site.js: window.SodaSite with data empty or an https:// link ({m.group(1) if m else 'no data:'})")
    # site.js is public (in git and on GitHub Pages): a read pass (SAS) in its link would publish the pass
    c.ok(m is None or not re.search(r"[?&](sig|sv|se|sp|skoid)=", m.group(1)),
         "site.js: the data link has no read pass (sig=/sv=/se=/sp=); only a public-read folder goes there")
    c.ok("site.js" in names and "config.js" in names and names.index("site.js") < names.index("config.js") < names.index("viewer.js"),
         "index.html: site.js loads before config.js and viewer.js")
    c.ok(all(k in js for k in ("window.SodaSite", 'params.get("static")', "TOKEN_URL || (cfg.datasets", '|| SITE_DATA;')),
         "viewer.js: site.js's folder is used only when no link, ?static=, config.js data/datasets or password names data")
    c.ok("No datasets: add entries to config.js" not in js and "open the SODA link you were given" in js,
         "viewer.js: the page with no data says in plain words how to open a link")
    dj = (APP / "diag.js").read_text(encoding="utf-8")
    c.ok('get("diag") !== "1"' in dj, "diag.js: the panel must stay off without ?diag=1")
    node = shutil.which("node")
    if not node:
        print("    node not found: tests/diag_test.js skipped")
        return
    p = subprocess.run([node, str(REPO / "tests" / "diag_test.js")], cwd=str(REPO),
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace")
    with open(log_path, "a", encoding="utf-8") as lg:
        lg.write(p.stdout)
    print("    " + (p.stdout.strip().splitlines()[-1] if p.stdout.strip() else "(no output)"))
    c.eq(p.returncode, 0, "node tests/diag_test.js (exit code)")
    for f in ("viewer.js", "diag.js", "lasso_worker.js", "umapfit.js", "leiden_worker.js"):
        q = subprocess.run([node, "--check", str(APP / f)], cwd=str(REPO), stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        c.eq(q.returncode, 0, f"node --check {f}")
    # old exports still open (SKIPPED without Playwright); and the current, packed public ones: a gene by Range,
    # with a read pass on a ?data= link, and from a host that ignores Range
    out = node_test(c, log_path, "compat_test.js", ["--python", sys.executable],
                    "the page opens exports made by older versions and packed current ones")
    for ln in [x for x in out.splitlines() if x.startswith("packed exports (")]:
        print("    " + ln)
    # the same page in Chromium, WebKit (Safari's engine), Edge, Firefox and phone-sized WebKit, Chromium and Firefox
    # (Firefox: the size only, no touch), on the public exports (SKIPPED without Playwright; a browser that is not
    # installed, such as Edge or Playwright's Firefox, is SKIPPED on its own and the rest still pass)
    node_test(c, log_path, "browsers_test.js", ["--shots", str(root / "browsers_shots")],
              "the page works in Chromium, WebKit, Edge and Firefox")
    # Draw a region: the exact path against its fixture, fast mode in the worker and on a public export (SKIPPED
    # without Playwright, or without umap-js from the CDN)
    node_test(c, log_path, "drawregion_test.js", ["--shots", str(root / "browsers_shots")],
              "Draw a region's exact path = its fixture; fast mode labelled, saved and undone")



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
    js = {f: (APP / f).read_text(encoding="utf-8") for f in ("layout.js", "viewer.js", "volcano.js", "plots.js", "qc.js", "roi.js", "snapshot.js", "cellsets.js", "binpaint.js", "selections.js")}
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
    tabs = {"view": "tabView", "analyze": "tabAnalyze", "share": "tabShare"}
    c.ok("tabClusters" not in by_id and "tabSave" not in by_id, "index.html: three tabs (Clusters joined Analyze, Save became Share in 0.10.0)")
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
    # every section names one of the three tabs: in the page, and in the scripts that build one
    secs = [(a.get("id", "(no id)"), a.get("data-tab")) for tag, a in t.tags
            if tag == "section" and "panel" in (a.get("class") or "").split()]
    c.ok(len(secs) >= 10, f"index.html: the sections are there ({len(secs)} section.panel)")
    for sid, tab in secs:
        c.ok(tab in tabs and sid != "(no id)", f"index.html: section #{sid} has an id and names its tab (data-tab={tab!r})")
    made = re.findall(r'newSection\("(\w+)", "(\w+)"', js["layout.js"])
    c.ok({i for i, _ in made} == {"panelSave", "panelVolcano"} and all(tab in tabs for _, tab in made),
         f"layout.js: the sections it makes (Save and load, Volcano plot) name a tab: {made}")
    c.ok('sec.dataset.tab = "analyze"' in js["plots.js"] and "SodaLayout.adopt(sec)" in js["plots.js"],
         "plots.js: its section names its tab and is handed to the layout")
    c.ok("SodaLayout.adopt(block)" in js["volcano.js"], "volcano.js: its block is handed to the layout")
    c.ok('sec.dataset.tab = "share"' in js["snapshot.js"] and "SodaLayout.adopt(sec)" in js["snapshot.js"],
         "snapshot.js: its Save image section names its tab and is handed to the layout")
    c.ok("sec.dataset.tab = cfg.tab" in js["cellsets.js"] and "SodaLayout.adopt(sec)" in js["cellsets.js"]
         and 'id: "panelSets", tab: "analyze"' in js["cellsets.js"] and "panelRunSets" not in js["cellsets.js"],
         "cellsets.js: its Cell sets section (Analyze tab) names its tab and is handed to the layout")
    c.ok('sec.dataset.tab = "view"' in js["binpaint.js"] and "SodaLayout.adopt(sec)" in js["binpaint.js"],
         "binpaint.js: its Colour by a value section names its tab and is handed to the layout")
    c.ok("binpaint.js" in names and not re.search(r"\bInt(8|16|32)Array\b", js["binpaint.js"]),
         "index.html loads binpaint.js, which uses only unsigned typed arrays (deid.py check reads the signed names as sample IDs)")
    c.ok("measure.js" not in names and not (APP / "measure.js").exists() and "panelMeasure" not in js["layout.js"],
         "Measure (measure.js) is gone (cut in the final release): no script, no section")
    c.ok('sec.dataset.tab = "analyze"' in js["selections.js"] and "SodaLayout.adopt(sec)" in js["selections.js"],
         "selections.js: its Export shapes section names its tab and is handed to the layout")
    # ids the layout makes must be new; ids it moves must exist (in the page or in the script that builds them)
    new_ids = set(i for i, _ in made) | set(re.findall(r'el\("\w+", "[^"]*", "(\w+)"\)', js["layout.js"]))
    clash = sorted(new_ids & set(ids))
    c.ok(not clash, f"layout.js makes ids that index.html already has: {clash}")
    built = set(re.findall(r'id="(\w+)"', js["volcano.js"] + js["plots.js"] + js["snapshot.js"] + js["cellsets.js"] + js["binpaint.js"] + js["selections.js"])) | set(ids) | new_ids
    block = js["layout.js"][js["layout.js"].index("const SECTIONS = {"):js["layout.js"].index("const ORDER = {")]
    moved = set(re.findall(r'"(\w+)"', " ".join(re.findall(r'(?:more|notes): \[([^\]]*)\]', block))))
    c.ok(len(moved) >= 15 and not (moved - built), f"layout.js moves controls that do not exist: {sorted(moved - built)}")
    order = js["layout.js"][js["layout.js"].index("const ORDER = {"):js["layout.js"].index("const KEY = ")]
    ordered = set(re.findall(r'"(panel\w+)"', order))
    more = order[order.index("const MORE = {"):]
    c.ok(set(re.findall(r'"(panel\w+)"', more)) == {"panelPaint", "panelOwn", "panelSets", "panelBinQc"},
         "layout.js: Colour by a value, Your own results, Cell sets and the quality filter are behind their tab's More options")
    c.ok("function buildTabMore(tab)" in js["layout.js"] and "if (behindMore(sec)) setTabMore(tab, true, false);" in js["layout.js"],
         "layout.js: a tab's More options is built, and show() opens it for a section behind it")
    known = {sid for sid, _ in secs} | {i for i, _ in made} | {"panelPlots", "panelImage", "panelSets", "panelPaint", "panelShapes",
                                                               "panelOwn"}
    c.ok(not (ordered - known), f"layout.js orders sections that do not exist: {sorted(ordered - known)}")
    for needle, f, what in (("lay.opening(false, tok)", "viewer.js", "tells the layout when a dataset starts and stops opening"),
                            ("SodaLayout.wording", "viewer.js", "hands its wording function to the layout"),
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
    matched to the bins in the browser with no new export. Wired into the page (its section in the Analyze tab,
    the hooks in viewer.js and binpaint.js, the link extension), no signed typed arrays, and, when node is
    installed, tests/ownres_browser_test.js in Chromium, WebKit and Firefox (SKIPPED when not installed) on the
    public visiumhd_pub (its barcode lists): Seurat and scanpy files, .gz, partial matches, roles put right, colours by
    pixel, the volcano and cell sets on a loaded grouping, Save JSON and the view link, junk files, 375 px."""
    html = (APP / "index.html").read_text(encoding="utf-8")
    names = [r.partition("?")[0] for r in re.findall(r'<(?:script|link)[^>]+?(?:src|href)="([^"]+)"', html)]
    js = {f: (APP / f).read_text(encoding="utf-8") for f in ("ownres.js", "viewer.js", "binpaint.js", "layout.js")}
    c.ok("ownres.js" in names and "ownres.css" in names and "binpaint.js" in names
         and names.index("layout.js") < names.index("ownres.js") < names.index("viewer.js"),
         "index.html loads ownres.css and ownres.js (after layout.js, before viewer.js, which starts it)")
    o = js["ownres.js"]
    c.ok('sec.dataset.tab = "analyze"' in o and "SodaLayout.adopt(sec)" in o,
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
                            ("panelOwn", "layout.js", "has the section's help and its place in the Analyze tab")):
        c.ok(needle in js[f], f"{f}: {what} ({needle!r} not found)")
    node_test(c, log_path, "ownres_browser_test.js", ["--shots", str(root / "ownres_shots"), "--browsers", "chromium,webkit,firefox"],
              "your own results loaded in the page, in Chromium, WebKit and Firefox")



def check_sharing(c: Checker, log_path: Path):
    """The link that restores a view (viewlink.js): wired into the page, its ids, the hooks in viewer.js and
    layout.js, the page's Content-Security-Policy, and the node unit test of the link's encoding."""
    html = (APP / "index.html").read_text(encoding="utf-8")
    names = [r.partition("?")[0] for r in re.findall(r'<(?:script|link)[^>]+?(?:src|href)="([^"]+)"', html)]
    c.ok("viewlink.js" in names and (APP / "viewlink.js").exists(), "index.html loads viewlink.js")
    c.ok("viewlink.js" in names and names.index("layout.js") < names.index("viewlink.js") < names.index("viewer.js")
         and names.index("dspicker.js") < names.index("viewlink.js"),
         "viewlink.js loads after layout.js and dspicker.js and before viewer.js (which hands it its functions)")
    c.ok(not any(n in names for n in ("compare.js", "compare.css", "offline.js")) and not (APP / "compare.js").exists()
         and not (APP / "offline.js").exists(), "index.html loads no compare.js / compare.css / offline.js (removed in 0.10.0)")
    js = {f: (APP / f).read_text(encoding="utf-8") for f in ("viewlink.js", "viewer.js", "layout.js")}
    page_ids = set(re.findall(r'\bid="(\w+)"', html))
    made = set(re.findall(r'\.id = "(\w+)"', js["viewlink.js"]))
    c.ok("saveLinkBox" in made, f"viewlink.js builds the box for a link the browser would not copy: {sorted(made)}")
    c.ok(not (made & page_ids), f"ids built by viewlink.js that index.html already has: {sorted(made & page_ids)}")
    c.ok(".save-link-box" in (APP / "layout.css").read_text(encoding="utf-8"), "layout.css styles the link box")
    lay_ids = set(re.findall(r'el\("\w+", "[^"]*", "(\w+)"\)', js["layout.js"]))
    c.ok({"btnCopyLink", "saveMsg"} <= lay_ids and not (made & lay_ids),
         "layout.js still builds the Copy link button and its message line (viewlink.js takes the button over)")
    c.ok('$("btnCopyLink")' in js["viewlink.js"] and "Copy link to this view" in js["viewlink.js"],
         "viewlink.js re-labels and re-wires the Copy link button")
    for needle, what in (("function sessionState(full)", "one serialiser for Save JSON and the link"),
                         ("JSON.stringify(sessionState(true), null, 1)", "Save JSON writes sessionState(true)"),
                         ("async function sessionApply(j, run)", "Load JSON and a link go back in through one function"),
                         ("window.SodaViewLink.capture()", "the view block is captured into the saved state"),
                         ("window.SodaViewLink.apply(j.view, j)", "the view block is applied from a file and from a link"),
                         ("function setCamera(viewState, viewId)", "one function moves the camera (mouse, link)"),
                         ("window.SodaViewLink.init(share)", "viewlink.js is handed viewer.js's functions"),
                         ('new CustomEvent("soda:opened"', "the scripts are told when a dataset has opened")):
        c.ok(needle in js["viewer.js"], f"viewer.js: {what} ({needle!r} not found)")
    # the page's Content-Security-Policy (no inline script, no handler attributes, no eval): nothing in the file needs more
    t = js["viewlink.js"]
    c.ok(not re.search(r"""<[a-z][^>`]*\son[a-z]+\s*=""", t) and "javascript:" not in t
         and not re.search(r"""\beval\(|new Function\(|setTimeout\(\s*["'`]""", t),
         "viewlink.js: no inline handler attribute, no eval, no javascript: address")
    c.ok(not re.search(r"\b(Int8Array|Int16Array|Int32Array|BigInt64Array)\b", t),  # deid:ok (type names)
         "viewlink.js: only unsigned typed arrays (the de-identification check reads the signed names as sample IDs)")
    lines = [ln for ln in t.splitlines() if re.search(r"\b(local|session)Storage\.", ln)]
    c.ok(all("try" in ln for ln in lines), f"viewlink.js: a storage access outside try / catch: {[ln.strip()[:60] for ln in lines if 'try' not in ln]}")
    c.ok('KEEP_PARAMS = ["data", "static", "qc"]' in t and "location.hash" in t and "fetch(" not in t,
         "viewlink.js: the view travels after #, the link keeps ?data= / ?static= / ?qc=, and it fetches nothing itself")
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



def check_libs(c: Checker):
    """The pinned libraries agree everywhere, nothing widens the page's Content-Security-Policy, and sw.js only
    removes the service worker and caches of the old offline use (removed in 0.10.0)."""
    import tempfile
    html = (APP / "index.html").read_text(encoding="utf-8")
    js = {f: (APP / f).read_text(encoding="utf-8") for f in ("libs.js", "sw.js", "lasso_worker.js", "viewer.js", "layout.js", "diag.js")}
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
    c.ok(not any("cdn.jsdelivr.net" in s for s in srcs), "index.html names no CDN script itself (libs.js writes the fallback)")
    c.ok("cdn.jsdelivr.net/npm/" in js["libs.js"] and "integrity=" in js["libs.js"] and 'crossorigin="anonymous"' in js["libs.js"],
         "libs.js: the CDN fallback carries the same SHA-384 (integrity, crossorigin)")
    csp = re.search(r'http-equiv="Content-Security-Policy" content="([^"]+)"', html)
    csp = csp.group(1) if csp else ""
    c.ok("script-src 'self' https://cdn.jsdelivr.net;" in csp and "worker-src 'self' blob:;" in csp and "unsafe-eval" not in csp,
         f"index.html: the Content-Security-Policy is no wider than before (no 'unsafe-eval', no 'wasm-unsafe-eval'): {csp[:120]}")
    for f in ("libs.js", "sw.js"):
        t = js[f]
        c.ok(not re.search(r"""\beval\(|new Function\(|setTimeout\(\s*["'`]""", t) and "javascript:" not in t
             and not re.search(r"\b(Int8Array|Int16Array|Int32Array|BigInt64Array)\b", t),  # deid:ok (type names)
             f"{f}: no eval, no javascript: address, only unsigned typed arrays")
        lines = [ln for ln in t.splitlines() if re.search(r"\b(local|session)Storage\.", ln)]
        c.ok(not lines, f"{f}: no local / session storage: {lines[:2]}")
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
    # every worker and import is named "<name>.js?v=...", so a new publish is never mixed with a cached old file
    allj = "".join((APP / f).read_text(encoding="utf-8") for f in sorted(os.listdir(APP))
                   if f.endswith(".js") and f not in ("config.js", "sw.js"))
    wk = re.findall(r"""["'`]([\w-]+_worker\.js)(\?v=[\w.-]+)?["'`]""", allj)
    c.ok(len(wk) >= 3 and all(v for _, v in wk), f"every worker is named with a ?v= string: {sorted(set(n + v for n, v in wk))}")
    imp = re.findall(r"""importScripts\(\s*["'`]([^"'`]+)["'`]""", "".join((APP / f).read_text(encoding="utf-8") for f in sorted(os.listdir(APP)) if f.endswith("_worker.js")))
    c.ok(imp and all(re.fullmatch(r"[\w-]+\.js\?v=[\w.-]+", i) for i in imp),
         f"every importScripts of the page's workers names a local script with a ?v= string: {imp}")
    sw = js["sw.js"]
    c.ok("self.registration.unregister()" in sw and "caches.delete(" in sw and '"fetch"' not in sw and "importScripts" not in sw
         and all(p in sw for p in ('"soda-app:"', '"soda-lib:"', '"soda-data:"')) and "+ SCOPE +" in sw,
         "sw.js: deletes only this folder's soda-app: / soda-lib: / soda-data: caches, unregisters itself, answers no request")
    c.ok("navigator.serviceWorker.register" not in js["viewer.js"] and "r.update()" in js["viewer.js"],
         "viewer.js registers no service worker; it only asks an old one to update (to sw.js, which removes it)")
    c.ok("el.dataset.sodaLib" in js["diag.js"], "diag.js: the vendor/ probe is not reported as an error")
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


# ------------------------------------------------------------------ main
def main(argv=None) -> int:
    for stream in (sys.stdout, sys.stderr):      # a redirected Windows console (cp1252) must not stop the test
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(errors="replace")
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--raw-data", default=os.environ.get("SODA_RAW_DATA", ""),
                    help="the public datasets' raw files, one folder per id as REPRODUCE.md section 9 lays them out "
                         "(default SODA_RAW_DATA); without it the pipeline part is SKIPPED")
    ap.add_argument("--ids", default="", help="comma-separated public ids to prepare / check (default: all seven)")
    ap.add_argument("--root", default=str(Path(tempfile.gettempdir()) / "soda-smoke"),
                    help="where cache/, export/, the log and the browser screenshots go (default: soda-smoke in the temp folder)")
    ap.add_argument("--exports", default="",
                    help="check the public exports in this folder (read only) instead of making them; default "
                         "SODA_TEST_DATA when that is a folder")
    args = ap.parse_args(argv)
    for k, v in (("NODE_PATH", r"D:\SpatialVizBench\iphone\pw\node_modules"),
                 ("PLAYWRIGHT_BROWSERS_PATH", r"D:\SpatialVizBench\iphone\browsers")):
        if not os.environ.get(k) and Path(v).is_dir():   # unset: use this machine's bench Playwright, so no browser test skips
            os.environ[k] = v
            print(f"{k} was not set: using {v}")
    ids = [x for x in args.ids.split(",") if x] or list(PUBLIC)
    other = [x for x in ids if x not in PUBLIC]
    if other:
        print(f"--ids: {', '.join(other)} not public example ids ({', '.join(PUBLIC)})")
        return 2
    root = Path(args.root)
    root.mkdir(parents=True, exist_ok=True)
    log_path = root / "smoke_test.log"
    log_path.write_text("", encoding="utf-8")
    clustering = all(importlib.util.find_spec(mod) is not None for mod in CLUSTER_MODULES)
    print(f"python {sys.version.split()[0]} at {sys.executable}")
    print(f"clustering packages: {'yes' if clustering else 'NO - umap checks will be skipped'}")
    timings = {}
    t_all = time.time()

    # part 2: the pipeline, only with the raw files
    made = []
    if args.raw_data:
        raw = Path(args.raw_data)
        if not raw.is_dir():
            print(f"--raw-data / SODA_RAW_DATA: {raw} is not a folder")
            return 2
        print(f"pipeline: REPRODUCE.md's prepare + export on the raw files in {raw}, into {root}")
        made = run_pipeline(raw, root, ids, log_path, timings)
    else:
        print("SKIPPED pipeline (prepare + export): no raw-data folder (--raw-data or SODA_RAW_DATA, laid out as REPRODUCE.md section 9)")

    c = Checker()
    t0 = time.time()
    summary = {}
    # part 3: the exports just made, or a local copy made earlier
    if made:
        summary["exports"] = check_exports(c, root / "export", made, clustering, log_path)
        if all(ds in made for ds in HOSTED):        # the browser tests on this run's exports
            os.environ["SODA_TEST_DATA"] = str((root / "export").resolve())
    else:
        folder = args.exports or (os.environ.get("SODA_TEST_DATA", "") if not re.match(r"https?://", os.environ.get("SODA_TEST_DATA", ""), re.I) else "")
        have = [ds for ds in ids if folder and (Path(folder) / ds / "manifest.json").is_file()]
        if have:
            summary["exports"] = check_exports(c, Path(folder), have, clustering, log_path)
        else:
            print("SKIPPED export checks: no exports made here and no local export folder (--exports, or SODA_TEST_DATA "
                  "set to a folder) holding the public ids")
    print(f"browser tests' data: SODA_TEST_DATA={os.environ.get('SODA_TEST_DATA') or '(unset: the hosted public exports)'}"
          + ("; SODA_TEST_PRIVATE_DATA is set" if os.environ.get("SODA_TEST_PRIVATE_DATA") else ""))

    # part 1 (and the browser tests)
    print("checks: the notebook")
    check_notebook(c)
    print("checks: marker lists and name cleanup")
    check_marker_lists(c)
    print("checks: the page (index.html, viewer.js tile path, diag.js; old and packed exports, the browsers, Draw a region)")
    check_viewer_page(c, root, log_path)
    print("checks: cell outlines in the page; the instrument and SpatialData exports, outlines and controls in the browsers")
    check_format_pages(c, root, log_path)
    print("checks: BLAS / OpenMP thread cap")
    check_threadcap(c, log_path)
    print("checks: dataset picker, JavaScript syntax")
    check_picker(c)
    print("checks: tabbed layout (tabs, sections, ids)")
    check_layout(c)
    print("checks: sharing (the link that restores a view)")
    check_sharing(c, log_path)
    print("checks: your own results (ownres.js: a reader's CSV / TSV matched in the browser)")
    check_ownres(c, root, log_path)
    print("checks: pinned libraries (libs.js, vendor.py) and sw.js removing the old offline use")
    check_libs(c)
    timings["checks"] = time.time() - t0

    print("\ntimings: " + ", ".join(f"{k} {v:.0f}s" for k, v in timings.items())
          + f"; total {time.time() - t_all:.0f}s")
    (root / "smoke_test_summary.json").write_text(json.dumps(
        {"passed": not c.failed, "checks": c.n, "failed": c.failed, "clustering": clustering, "pipeline": made,
         "timingsSec": {k: round(v, 1) for k, v in timings.items()}, "exports": summary.get("exports", {})}, indent=1),
        encoding="utf-8")
    if c.failed:
        print(f"\nSMOKE TEST FAILED: {len(c.failed)} of {c.n} checks")
        for f in c.failed:
            print("  - " + f)
        return 1
    print(f"\nSMOKE TEST PASSED: {c.n} checks" + ("" if clustering else " (umap checks skipped)")
          + ("" if made else "; pipeline SKIPPED (no raw data)"))
    print(f"log: {log_path.resolve()}")
    return 0


def test_smoke():
    """pytest entry point."""
    assert main([]) == 0


if __name__ == "__main__":
    sys.exit(main())
