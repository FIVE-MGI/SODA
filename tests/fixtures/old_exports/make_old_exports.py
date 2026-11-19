#!/usr/bin/env python
r"""Rebuild the old-export fixtures that tests/compat_test.js opens with the current page.

Each fixture is a small SYNTHETIC export written by an OLD version of the code, taken read-only
from git history (git archive into a scratch folder; no commit is changed):

    v1-20261002  FIVE-MGI/spatial-viewer 09b12fa (2 Oct 2026, the first public version):
                 tiles + expr (genes.json, positions.bin, gene/) + overlays, no region.bin,
                 hvg/, umap/ or protein grid; plus a single-cell RNA-seq dataset (rnaseq.py).
    v2-20261104  spatial-viewer-v2 a5d0e2d (4 Nov 2026): adds region.bin / regions[], hvg/, umap/
                 (clusterings without name/source), celltypes, and the first protein grid layer
                 (grid020um, mat/ with p1 / p99_5, no normalize).
    v2-20261111  spatial-viewer-v2 9812dca (11 Nov 2026, before binqc/, the lab's own clusters and
                 Xenium cells as a layer): protein grid with bg-asinh-v1 and mean/sd, qc.json,
                 qc_flags.bin.

Sizes are kept small with the old code's own options: a crop of the demo images, 256 px tiles,
JPEG tiles (PNG ones of the noisy synthetic stains were ~2 MB per version) and only the 16 um
bin layer. The demo generator of each version makes the raw data
(v1 had none, so it uses the earliest one in spatial-viewer-v2, 891ca86). Only synthetic data
from that old commit's examples/make_demo.py and the made-up inputs below are used.

The synthetic demo generators (examples/make_*demo.py, examples/demo_*markers.json) are no longer in the
current tree: the tests use the public example exports now. The fixtures under v1-*/ and v2-*/ stay as
they were committed, and this script can still rebuild them, because it takes each generator and its
marker lists from the old commit named in VERSIONS (`git archive`), never from the working tree.

    python tests/fixtures/old_exports/make_old_exports.py --work <scratch folder outside the repo> \
        --v1-repo <clone of FIVE-MGI/spatial-viewer> [--install]

--install copies each <work>/<version>/export into tests/fixtures/old_exports/<version>/.
Run it with the Python environment of web/SpatialViewer/requirements.txt (+ clustering).
"""
from __future__ import annotations

import argparse
import io
import json
import shutil
import subprocess
import sys
import tarfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]

# version -> (repo: "v2" = this repository, "v1" = --v1-repo), commit, commit of the demo generator
VERSIONS = {
    "v1-20261002": ("v1", "09b12fa", "891ca86"),
    "v2-20261104": ("v2", "a5d0e2d", "a5d0e2d"),
    "v2-20261111": ("v2", "9812dca", "9812dca"),
}
VIS_CROP = (1744, 1744, 2256, 2256)       # 512 x 512 px of the 4000 x 4000 demo H&E
MX_CROP = (680, 680, 1320, 1320)          # 640 x 640 px of the 2000 x 2000 demo multiplex
TILE = 256


def run(cmd, log):
    log.write("$ " + " ".join(str(c) for c in cmd) + "\n"); log.flush()
    p = subprocess.run([str(c) for c in cmd], stdout=log, stderr=subprocess.STDOUT)
    if p.returncode:
        raise SystemExit(f"failed ({p.returncode}): {' '.join(str(c) for c in cmd)}")


def extract(repo: Path, commit: str, dest: Path):
    if (dest / "web" / "SpatialViewer" / "prepare.py").exists():
        return
    dest.mkdir(parents=True, exist_ok=True)
    tar = subprocess.run(["git", "-C", str(repo), "archive", commit], stdout=subprocess.PIPE, check=True).stdout
    with tarfile.open(fileobj=io.BytesIO(tar)) as t:
        t.extractall(dest, **({"filter": "data"} if hasattr(tarfile, "data_filter") else {}))


def regions_geojson(w: int, h: int) -> dict:
    """Three made-up regions in final image pixels: two halves and a small core (smallest wins)."""
    def poly(name, x0, y0, x1, y1):
        return {"type": "Feature", "properties": {"classification": {"name": name}},
                "geometry": {"type": "Polygon", "coordinates": [[[x0, y0], [x1, y0], [x1, y1], [x0, y1], [x0, y0]]]}}
    return {"type": "FeatureCollection", "features": [
        poly("Left", 0, 0, w // 2, h), poly("Right", w // 2, 0, w, h),
        poly("Core", w * 3 // 8, h * 3 // 8, w * 5 // 8, h * 5 // 8)]}


def tiny_rnaseq(folder: Path, n_cells=240, n_genes=40, seed=0):
    """A made-up 10x-style matrix folder + UMAP and obs CSVs (three groups of cells)."""
    rng = np.random.default_rng(seed)
    folder.mkdir(parents=True, exist_ok=True)
    grp = np.repeat(np.arange(3), n_cells // 3)
    lam = rng.gamma(2.0, 1.0, (3, n_genes))
    X = rng.poisson(lam[grp])                                  # cells x genes
    barcodes = [f"cell{i:04d}" for i in range(n_cells)]
    genes = [f"GENE{i + 1:03d}" for i in range(n_genes)]
    (folder / "barcodes.tsv").write_text("".join(b + "\n" for b in barcodes))
    (folder / "features.tsv").write_text("".join(f"DEMO{i + 1:06d}\t{g}\tGene Expression\n" for i, g in enumerate(genes)))
    nz = np.argwhere(X.T > 0)                                  # genes x cells, 1-based in the file
    lines = [f"{g + 1} {c + 1} {X[c, g]}" for g, c in nz]
    (folder / "matrix.mtx").write_text("%%MatrixMarket matrix coordinate integer general\n"
                                       f"{n_genes} {n_cells} {len(lines)}\n" + "\n".join(lines) + "\n")
    centre = np.array([[0, 0], [6, 0], [3, 5]], float)
    emb = centre[grp] + rng.normal(0, 1, (n_cells, 2))
    (folder / "umap.csv").write_text("barcode,UMAP1,UMAP2\n" + "".join(
        f"{b},{x:.4f},{y:.4f}\n" for b, (x, y) in zip(barcodes, emb)))
    names = np.array(["Group A", "Group B", "Group C"])
    score = rng.uniform(0, 1, n_cells)
    (folder / "obs.csv").write_text("barcode,group,score\n" + "".join(
        f"{b},{names[g]},{s:.4f}\n" for b, g, s in zip(barcodes, grp, score)))


def build(version: str, work: Path, v1_repo: Path | None, log) -> Path:
    which, commit, gen_commit = VERSIONS[version]
    repo = REPO if which == "v2" else v1_repo
    if repo is None or not (repo / ".git").exists():
        raise SystemExit(f"{version}: needs the FIVE-MGI/spatial-viewer clone (--v1-repo)")
    vdir = work / version
    src, gen = vdir / "src", vdir / "gen"
    extract(repo, commit, src)
    extract(REPO, gen_commit, gen)
    app = src / "web" / "SpatialViewer"
    raw, cache, export, inputs = vdir / "raw", vdir / "cache", vdir / "export", vdir / "inputs"
    for d in (cache, export):
        if d.exists():                         # an old export skips files already there: never mix two runs
            raise SystemExit(f"{d} exists: use a new --work folder (this script deletes nothing)")
    inputs.mkdir(parents=True, exist_ok=True)
    py = sys.executable
    if not (raw / "multiplex.ome.tif").exists():
        run([py, gen / "examples" / "make_demo.py", "--out", raw], log)
    vw, mw = VIS_CROP[2] - VIS_CROP[0], MX_CROP[2] - MX_CROP[0]
    (inputs / "vis_regions.geojson").write_text(json.dumps(regions_geojson(vw, vw)))
    (inputs / "mx_regions.geojson").write_text(json.dumps(regions_geojson(mw, mw)))
    crop = lambda c: ",".join(map(str, c))
    run([py, app / "prepare.py", "visium", "--id", "old_visiumhd", "--image", raw / "visiumhd" / "image.tif",
         "--outs", raw / "visiumhd" / "outs", "--out", cache, "--bins", "016um", "--crop", crop(VIS_CROP),
         "--tile", TILE, "--geojson", f"Regions={inputs / 'vis_regions.geojson'}"], log)
    run([py, app / "prepare.py", "tiff", "--id", "old_multiplex", "--src", raw / "multiplex.ome.tif",
         "--out", cache, "--crop", crop(MX_CROP), "--tile", TILE,
         "--geojson", f"Regions={inputs / 'mx_regions.geojson'}"], log)
    rna_markers = ["--markers", gen / "examples" / "demo_markers.json"] if which == "v2" else []
    prot = gen / "examples" / "demo_protein_markers.json"
    mx_markers = ["--markers", prot] if which == "v2" and prot.exists() else []
    run([py, app / "prepare.py", "export", "--id", "old_visiumhd", "--out", cache, "--dest", export, "--fmt", "jpg"]
        + rna_markers, log)
    run([py, app / "prepare.py", "export", "--id", "old_multiplex", "--out", cache, "--dest", export, "--fmt", "jpg"]
        + mx_markers, log)
    if which == "v1":
        tiny_rnaseq(inputs / "rnaseq")
        run([py, app / "prepare.py", "rnaseq", "--id", "old_rnaseq", "--label", "old_rnaseq", "--dest", export,
             "--mtx", inputs / "rnaseq", "--umap", inputs / "rnaseq" / "umap.csv",
             "--obs", inputs / "rnaseq" / "obs.csv"], log)
    return export


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work", required=True, help="scratch folder OUTSIDE the repo")
    ap.add_argument("--v1-repo", help="clone of FIVE-MGI/spatial-viewer (for v1-20261002)")
    ap.add_argument("--only", help="comma list of versions")
    ap.add_argument("--install", action="store_true", help="copy the exports into tests/fixtures/old_exports/")
    a = ap.parse_args(argv)
    work = Path(a.work).resolve()
    if REPO in work.parents or work == REPO:
        raise SystemExit("--work must be outside the repository")
    work.mkdir(parents=True, exist_ok=True)
    for version in (a.only.split(",") if a.only else VERSIONS):
        with open(work / f"{version}.log", "w", encoding="utf-8") as log:
            export = build(version, work, Path(a.v1_repo).resolve() if a.v1_repo else None, log)
        size = sum(f.stat().st_size for f in export.rglob("*") if f.is_file())
        print(f"{version}: {export}  {size / 1e6:.2f} MB")
        if a.install:
            dst = HERE / version
            if dst.exists():
                raise SystemExit(f"{dst} exists: move it away first (this script never deletes it)")
            shutil.copytree(export, dst)
            print(f"  installed -> {dst}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
