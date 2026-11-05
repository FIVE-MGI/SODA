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
It also checks the built-in marker lists and that the channel-name cleanup is the same
in celltypes.py and lasso_worker.js.

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
VIS, MX = "demo_visiumhd", "demo_multiplex"
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
    ]


def run(args: list[str], log_path: Path) -> float:
    print(f"$ python {' '.join(args)}", flush=True)
    t0 = time.time()
    with open(log_path, "a", encoding="utf-8") as lg:
        lg.write(f"\n$ python {' '.join(args)}\n")
        lg.flush()
        p = subprocess.run([sys.executable, *args], cwd=str(REPO), stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace")
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
    for p in [root / "raw", root / "export_sample" / MX] + [root / sub / ds for sub in ("cache", "export") for ds in (VIS, MX)]:
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
    res = {"bins": n, "medianUmi": float(np.median(lib))}
    if clustering:
        res.update(check_umap(c, d, L, pos, lib, truth, bin_um / 0.32))
    else:
        c.ok("umap" not in L, f"{key}: umap in manifest although clustering packages are missing")
        print("    SKIP umap checks: " + ", ".join(CLUSTER_MODULES) + " not all installed "
              "(pip install -r web/SpatialViewer/requirements-clustering.txt)")
    res["hvgGenes"] = hv
    return res


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
    c.ok(((lib >= idx["minUmi"]) == emb).all(),
         f"{key}: embedded rows are not exactly the bins with >= {idx['minUmi']} UMI")
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
    c.eq((idx.get("nBins"), idx.get("nEmbedded")), (n, n), f"{key}: umap nBins / nEmbedded vs squares")
    c.eq(idx.get("modality"), "protein", f"{key}: umap modality")
    c.eq(idx["params"].get("channels"), L["channels"], f"{key}: umap params.channels vs manifest channels")
    co = f32(ud / "coords.bin", 2)
    c.ok(len(co) == n and np.isfinite(co).all(), f"{key}: coords.bin rows / NaN (every square is embedded)")
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
        c.ok((lab < cl["n"]).all(), f"{key}: {cl['file']}: label >= n (or 0xFFFF)")
        sizes = np.bincount(lab, minlength=cl["n"])
        c.eq(sizes.tolist(), cl["sizes"], f"{key}: {cl['file']}: sizes")
        c.ok((np.diff(sizes) <= 0).all(), f"{key}: {cl['file']}: cluster 0 is not the largest")
        c.eq(len(cl["markers"]), cl["n"], f"{key}: {cl['file']}: one marker list per cluster")
        c.ok(all(m["gene"] in V for ms in cl["markers"] for m in ms), f"{key}: {cl['file']}: marker not a channel")
        same = float((lab[pairs[:, 0]] == lab[pairs[:, 1]]).mean())
        sh = rng.permutation(lab)
        chance = float((sh[pairs[:, 0]] == sh[pairs[:, 1]]).mean())
        a = ari(lab, truth)
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
    c.ok(((tb < len(names)) | (tb == 254)).all(), f"{key}: celltype.bin values out of range (255 = not embedded is impossible here)")
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
    c.eq((ib.get("nBins"), ib.get("nEmbedded"), ib.get("fitOn")), (n, n, SAMPLE), f"{key}: sample nBins / nEmbedded / fitOn")
    c.ok("fitted on" in ib["params"].get("fit", "") and ia["params"].get("fit") == "all rows",
         f"{key}: params.fit {ib['params'].get('fit')!r} / {ia['params'].get('fit')!r}")
    co = f32(b / "coords.bin", 2)
    c.ok(len(co) == n and np.isfinite(co).all(), f"{key}: sample coords.bin rows / NaN")
    out = {}
    for ca, cb in zip(ia["clusterings"], ib["clusterings"]):
        la = np.fromfile(str(a / ca["file"]), dtype="<u2").astype(np.int64)
        lb = np.fromfile(str(b / cb["file"]), dtype="<u2").astype(np.int64)
        sizes = np.bincount(lb, minlength=cb["n"])
        c.ok(len(lb) == n and (lb < cb["n"]).all() and sizes.tolist() == cb["sizes"] and (np.diff(sizes) <= 0).all(),
             f"{key}: sample {cb['file']}: labels / sizes / order")
        out[ca["resolution"]] = round(ari(la, lb), 3)
    best = max(out.values(), default=0)
    c.ok(best > 0.4, f"{key}: clusters from the sample disagree with the full fit (ARI {out})")
    print(f"    fitOn {ib.get('fitOn')} of {n}; ARI vs fitting on all: {out}")
    return {"ari": out}


def check_marker_lists(c: Checker):
    """The built-in lists load and are usable, and the channel-name cleanup is the same in Python and JS."""
    sys.path.insert(0, str(APP))
    try:
        import celltypes
        import protgrid  # noqa: F401  (used below)
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
    import protgrid
    rng = np.random.default_rng(0)
    W, H = 9600, 7200
    smooth = lambda n: np.convolve(rng.normal(0, 1, n + 400), np.ones(200) / 200, mode="valid")[:n] * 20 + 60
    tile = lambda n, sp: (1 + 0.15 * rng.normal(0, 1, n // sp + 1))[np.arange(n) // sp]   # each tile its own gain
    prof = lambda n, sp: smooth(n) * (tile(n, sp) if sp else 1) + rng.normal(0, 0.3, n)
    plain = [(prof(W, 0), prof(H, 0)) for _ in range(4)]
    tiled = [(prof(W, 480), prof(H, 360)) for _ in range(3)] + plain[:1]
    nm = ["A", "B", "C", "D"]
    sp, sn = protgrid.seam_check(tiled, nm, 0.5), protgrid.seam_check(plain, nm, 0.5)
    c.ok(sp["flagged"] and sp["columns"].get("spacingPx") in (476, 480, 484) and sp["rows"].get("spacingPx") in (356, 360, 364),
         f"protgrid.seam_check missed synthetic 480 x 360 px tiles: {sp}")
    c.ok(not sn["flagged"], f"protgrid.seam_check flagged a slide without tiles: {sn}")
    print(f"    seam check: tiles found at {sp['columns'].get('spacingPx')} x {sp['rows'].get('spacingPx')} px; "
          f"none on the plain profiles")


def check_publishable(c: Checker, root: Path, log_path: Path):
    export = root / "export"
    texts = [p for ds in (VIS, MX) for p in (export / ds).rglob("*")
             if p.is_file() and p.suffix in (".json", ".geojson", ".txt")]
    texts += [p for p in (root / "export_sample" / MX).rglob("*") if p.is_file() and p.suffix == ".json"]
    bc = [p for p in texts if BARCODE.search(p.read_text(encoding="utf-8", errors="ignore"))]
    c.ok(not bc, f"export text files contain barcodes: {bc[:3]}")
    lp = [p for p in texts if LOCAL_PATH.search(p.read_text(encoding="utf-8", errors="ignore"))]
    c.ok(not lp, f"export text files contain local paths: {lp[:3]}")
    args = ["web/SpatialViewer/deid.py", "check", str(export / VIS), str(export / MX)]
    if (root / "export_sample" / MX).exists():
        args.append(str(root / "export_sample" / MX))
    print(f"$ python {' '.join(args)}", flush=True)
    p = subprocess.run([sys.executable, *args], cwd=str(REPO), stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                       text=True, encoding="utf-8", errors="replace")
    with open(log_path, "a", encoding="utf-8") as lg:
        lg.write(p.stdout)
    print("    " + p.stdout.strip().splitlines()[-1] if p.stdout.strip() else "    (no output)")
    c.eq(p.returncode, 0, "deid.py check on the export (exit code)")


# ------------------------------------------------------------------ main
def main(argv=None) -> int:
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
        log_path.write_text("", encoding="utf-8")
        for cmd in demo_commands(root):
            label = Path(cmd[0]).stem if cmd[1].startswith("--") else f"{cmd[1]} {cmd[3]}"
            label += " (sample)" if "--umap-sample" in cmd else ""
            timings[label] = run(cmd, log_path)

    c = Checker()
    t0 = time.time()
    print("checks: demo_visiumhd")
    summary = check_visium(c, root, clustering)
    print("checks: demo_multiplex")
    summary[GRID_KEY] = check_multiplex(c, root, clustering)
    print("checks: marker lists and name cleanup")
    check_marker_lists(c)
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
