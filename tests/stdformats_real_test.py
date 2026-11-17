#!/usr/bin/env python3
r"""Standard formats on real public example data (not part of the smoke test).

Runs only when the downloads are there (default D:\SpatialVizBench\stdreal\data; --root for another
place, SODA_STDREAL in the environment too); otherwise prints SKIPPED and exits 0. The datasets are
the projects' own examples, each in its own folder under data\ (how they were fetched and their
sizes: downloads.md next to data\):

    sd_mibitof        SpatialData sandbox MIBI-TOF (zarr v3): a table on a labels element (refused: labels
                      are not read, so this checks the refusal)
    sd_merfish        SpatialData sandbox MERFISH (zarr v3): circles, points, an anatomical shapes element
    sd_visium         SpatialData sandbox Visium (zarr v3): classic 55 µm spots, no pixel size
    sd_visium_hd      SpatialData sandbox Visium HD 4.0.1 (zarr v3): cell + nucleus segmentations
    sd_xenium_rep1    SpatialData sandbox Xenium rep1 (zarr v3): 168k cells, 43 M transcripts (--large)
    squidpy_visium_full  Squidpy visium_hne_adata.h5ad (the V1_Adult_Mouse_Brain section)
    idr_6001240_v04 / _v05   IDR 6001240 as OME-NGFF 0.4 (zarr v2) and 0.5 (zarr v3, sharded)
    tenx_visium_mouse_brain  10x V1_Adult_Mouse_Brain, classic Visium (Space Ranger 1.1.0 outputs)

For each one that is present: detect (certain, the right kind), prepare, export, then checks
against the source files read directly (positions, counts, pixel size, spot size),
deid.py check on the export, and the page in Chromium (tests/stdformats_real_browser.js, when
node and Playwright are found: NODE_PATH / PLAYWRIGHT_BROWSERS_PATH as for the smoke test).

    python tests/stdformats_real_test.py [--root D:\SpatialVizBench\stdreal] [--work <dir>] [--only a,b] [--large]

Run it on the E-cores on this machine (start /affinity FF0000 /wait /b python ...): the children
inherit the affinity. The downloads are only read (h5py, zarr, pyarrow, zarrv3.py), never run.
Writes <work>\cache, <work>\export (default <root>\realtest) and <work>\shots."""
from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import struct
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
APP = REPO / "web" / "SpatialViewer"
sys.path.insert(0, str(APP))
DEFAULT_ROOT = Path(os.environ.get("SODA_STDREAL", r"D:\SpatialVizBench\stdreal"))

# id -> (download folder, path under it to give prepare / detect, detect kind, prepare subcommand + options, checks)
SD = "extracted/data.zarr"
DATASETS = {
    "mibitof_p8": ("sd_mibitof", SD, "spatialdata", ["spatialdata", "--region", "point8_labels"],
                   {"refused": ["labels element", "not read yet"]}),
    "sd_merfish": ("sd_merfish", SD, "spatialdata", ["spatialdata", "--coords-um", "1"],
                   {"sdata": ("table", "cells", "shapes"), "pixel": 1 / 0.27165}),
    "sd_visium": ("sd_visium", SD, "spatialdata", ["spatialdata", "--region", "ST8059048"],
                  {"sdata": ("table", "ST8059048", "shapes"), "visium55": True}),
    "sd_visium_hd": ("sd_visium_hd", SD, "spatialdata", ["spatialdata", "--coordinate-system", "Visium_HD_3prime_Mouse_Brain"],
                     {"sdata": ("cell_segmentations", "Visium_HD_3prime_Mouse_Brain_cell_segmentations", "shapes"),
                      "layers": ["cells", "nuclei"]}),
    "sd_xenium": ("sd_xenium_rep1", SD, "spatialdata", ["spatialdata", "--pixel-size", "0.2125"],
                  {"sdata": ("table", "cell_circles", "shapes"), "large": True, "molecules": True}),
    "squidpy_visium_full": ("squidpy_visium_full", "visium_hne_adata.h5ad", "anndata", ["anndata"], {}),
    "idr_v04": ("idr_6001240_v04", "download", "omezarr", ["omezarr"], {"pixel": 0.3604}),
    "idr_v05": ("idr_6001240_v05", "download", "omezarr", ["omezarr"], {"pixel": 0.3604, "same_as": "idr_v04"}),
    "tenx_visium": ("tenx_visium_mouse_brain", ".", "visium", None, {"tenx": True}),
}


class Checker:
    def __init__(self):
        self.n, self.failed = 0, []

    def ok(self, cond, what):
        self.n += 1
        print(f"  {'ok  ' if cond else 'FAIL'} {what}", flush=True)
        if not cond:
            self.failed.append(what)
        return bool(cond)


def prepare(args, log: Path) -> tuple[int, str, float]:
    t0 = time.time()
    p = subprocess.run([sys.executable, str(APP / "prepare.py"), *args], cwd=str(REPO), stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace",
                       env={**os.environ, "PYTHONIOENCODING": "utf-8"})
    with open(log, "a", encoding="utf-8") as fh:
        fh.write(f"\n$ prepare.py {' '.join(args)}\n{p.stdout}")
    return p.returncode, p.stdout, time.time() - t0


# ---------------------------------------------------------------- reading the sources
def _node(p: Path):
    import zarr
    import zarrv3
    return zarrv3.open_node(p) if (p / "zarr.json").exists() else zarr.open(str(p), mode="r")


def _attrs(n) -> dict:
    a = n.attrs
    return a.asdict() if hasattr(a, "asdict") else dict(a)


def _strings(a) -> list[str]:
    return [x.decode() if isinstance(x, bytes) else str(x) for x in np.asarray(a[...], dtype=object).ravel()]


def _obs(t, name):
    a = t["obs"][name]
    if hasattr(a, "keys") and "codes" in list(a.keys()):
        cats, codes = _strings(a["categories"]), np.asarray(a["codes"][...])
        return [cats[c] if c >= 0 else None for c in codes]
    return [x.decode() if isinstance(x, bytes) else x for x in np.asarray(a[...], dtype=object).ravel()]


def _wkb_centres(blobs) -> np.ndarray:
    """Point centres, or polygon (outer ring) area centroids, from WKB."""
    out = np.full((len(blobs), 2), np.nan)
    for i, b in enumerate(blobs):
        if b is None:
            continue
        e = "<" if b[0] == 1 else ">"
        typ = struct.unpack(e + "I", b[1:5])[0] & 0xFFFF
        if typ == 1:
            out[i] = struct.unpack(e + "dd", b[5:21])
        elif typ in (3, 6):
            off = 5
            if typ == 6:
                off += 4
                e = "<" if b[off] == 1 else ">"
                off += 5
            n = struct.unpack(e + "I", b[off + 4:off + 8])[0]
            xy = np.frombuffer(b[off + 8:off + 8 + 16 * n], dtype=e + "f8").reshape(n, 2)
            x, y = xy[:, 0], xy[:, 1]
            cr = x[:-1] * y[1:] - x[1:] * y[:-1]
            a = cr.sum() / 2
            out[i] = ((x[:-1] + x[1:]) * cr).sum() / (6 * a), ((y[:-1] + y[1:]) * cr).sum() / (6 * a)
    return out


def check_sdata(c: Checker, sid: str, src: Path, cache: Path, spec: tuple, opts: dict):
    """Rows, per-row totals, positions (one affine map of the source centres), dot size."""
    import pyarrow.parquet as pq
    import zarr
    tname, region, ekind = spec
    m = json.loads((cache / sid / "manifest.json").read_text(encoding="utf-8"))
    key = m["spatialdata"]["layers"][0]["key"]
    g = zarr.open_group(str(cache / sid / "expr" / f"{key}.zarr"), mode="r")
    bc = _strings(g["barcodes"])
    xy = np.asarray(g["xy"][...], np.float64)
    t = _node(src / "tables" / tname)
    sa = _attrs(t)
    sa = sa.get("spatialdata_attrs") or sa
    if "instance_key" not in sa:
        sa = _attrs(t["uns"]["spatialdata_attrs"])
    regs, inst = _obs(t, sa["region_key"]), _obs(t, sa["instance_key"])
    sel = [i for i, r in enumerate(regs) if r == region]
    c.ok(len(bc) == len(sel), f"{sid}: {len(bc):,} rows = the table rows annotating its region ({len(sel):,})")
    va = _attrs(t["var"])
    vn = _strings(t["var"][va.get("_index", "_index")])
    keep = set(g.attrs["genes"])
    cj = np.array([n in keep for n in vn])
    X = t["X"]
    if hasattr(X, "keys"):
        data, ind, ptr = (np.asarray(X[k][...]) for k in ("data", "indices", "indptr"))
        rows = np.repeat(np.arange(len(ptr) - 1), np.diff(ptr))
        tot_src = np.bincount(rows, weights=np.where(cj[ind], data.astype(np.float64), 0.0), minlength=len(ptr) - 1)[sel]
    else:
        tot_src = np.asarray(X[...], np.float64)[:, cj].sum(1)[sel]
    tot = np.bincount(np.asarray(g["indices"][...]), weights=np.asarray(g["data"][...], np.float64), minlength=len(bc))
    c.ok(np.allclose(tot, tot_src, rtol=1e-5, atol=1e-3),
         f"{sid}: every row's total = the table's X row sum over the genes kept ({int(cj.sum())} of {len(vn)}; controls left out)")
    tab = pq.read_table(src / "shapes" / region / "shapes.parquet")
    cols = tab.column_names
    ic = ((tab.schema.pandas_metadata or {}).get("index_columns") or [None])[0]
    sidx = [str(v) for v in tab.column(ic).to_pylist()] if isinstance(ic, str) else [str(i) for i in range(tab.num_rows)]
    cen = _wkb_centres(tab.column("geometry").to_pylist())
    pos = {k: i for i, k in enumerate(sidx)}
    want = np.array([cen[pos[str(inst[i])]] if str(inst[i]) in pos else [np.nan, np.nan] for i in sel])
    ok = np.isfinite(want).all(1) & np.isfinite(xy).all(1)
    A = np.c_[want[ok], np.ones(ok.sum())]
    coef = np.linalg.lstsq(A, xy[ok], rcond=None)[0]
    res = float(np.abs(A @ coef - xy[ok]).max())
    sc = float(np.sqrt(abs(np.linalg.det(coef[:2]))))
    c.ok(ok.sum() == len(sel) and res < 0.05, f"{sid}: {ok.sum():,} positions = one affine map of the source centres "
                                             f"(max residual {res:.4f} image px, scale {sc:.5f})")
    pu = m.get("pixelSizeUm")
    if "pixel" in opts:
        c.ok(pu and abs(pu - opts["pixel"]) / opts["pixel"] < 1e-3, f"{sid}: pixel size {pu} µm (want {opts['pixel']:.4f})")
    e = next(e for e in m["expr"] if e["key"] == key)
    if "radius" in cols:
        r = np.asarray(tab.column("radius").to_pylist(), np.float64)
        d = float(np.median(r) * 2 * sc)
        c.ok(abs(d - e["spotDiameterPx"]) / e["spotDiameterPx"] < 0.01, f"{sid}: dot {e['spotDiameterPx']:.2f} px = the circles' diameter ({d:.2f})")
        if opts.get("visium55"):
            c.ok(pu and abs(d * pu - 55.0) < 0.01 and "classic Visium" in str((m.get("spatialdata") or {}).get("pixelSize", "")),
                 f"{sid}: no pixel size in the store: 55 µm Visium spots assumed ({d * (pu or 0):.3f} µm a spot)")
    if opts.get("layers"):
        c.ok([x["key"] for x in m["expr"]] == opts["layers"], f"{sid}: layers {[x['key'] for x in m['expr']]}")
    if opts.get("molecules"):
        c.ok((m.get("molecules") or {}).get("n", 0) > 0, f"{sid}: the points are the molecules ({(m.get('molecules') or {}).get('n')})")


def check_tenx(c: Checker, d: Path, cache: Path, export: Path, work: Path):
    """Classic Visium through prepare.py visium: positions, pixel size, spot size, gene totals; the same section
    through the AnnData path (Squidpy's visium_hne_adata.h5ad) when it was prepared."""
    import h5py
    import zarr
    from scipy import sparse
    sp = d / "spatial"
    sf = json.loads((sp / "scalefactors_json.json").read_text())
    pos = {}
    tp = sp / "tissue_positions_list.csv" if (sp / "tissue_positions_list.csv").exists() else sp / "tissue_positions.csv"
    with open(tp, newline="") as fh:
        for row in csv.reader(fh):
            if row[0] == "barcode":
                continue
            pos[row[0]] = (float(row[5]), float(row[4]), int(row[1]))
    m = json.loads((export / "tenx_visium" / "manifest.json").read_text(encoding="utf-8"))
    L = next(e for e in m["expr"] if e["key"] == "spots")
    um = 55.0 / sf["spot_diameter_fullres"]
    c.ok(abs(m["pixelSizeUm"] - um) / um < 1e-6, f"tenx_visium: pixel size {m['pixelSizeUm']:.5f} µm = 55 µm / spot_diameter_fullres")
    c.ok(abs(L["spotDiameterPx"] - sf["spot_diameter_fullres"]) < 1e-3 and abs(L["spotDiameterPx"] * m["pixelSizeUm"] - 55) < 1e-6,
         f"tenx_visium: spot {L['spotDiameterPx']:.3f} px = spot_diameter_fullres = 55 µm")
    g = zarr.open_group(str(cache / "tenx_visium" / "expr" / "spots.zarr"), mode="r")
    bc = _strings(g["barcodes"])
    xy = np.asarray(g["xy"][...], np.float64)
    off = float(np.abs(xy - np.array([pos[b][:2] for b in bc])).max())
    c.ok(off < 1e-3, f"tenx_visium: every spot at tissue_positions pxl_col / pxl_row_in_fullres (max off {off:.1e} px)")
    nin = sum(1 for v in pos.values() if v[2] == 1)
    c.ok(len(bc) == nin and all(pos[b][2] == 1 for b in bc), f"tenx_visium: the {len(bc):,} in-tissue spots")
    h5 = next(d.glob("*filtered_feature_bc_matrix.h5"))
    with h5py.File(h5, "r") as h:
        mat = h["matrix"]
        names = [x.decode() for x in mat["features"]["name"][:]]
        X = sparse.csc_matrix((mat["data"][:], mat["indices"][:], mat["indptr"][:]), shape=tuple(mat["shape"][:]))
    col = np.asarray(X.sum(1)).ravel()
    gj = json.loads((export / "tenx_visium" / "expr" / "spots" / "genes.json").read_text(encoding="utf-8"))
    top = sorted(gj, key=lambda e: -e["t"])[:5]
    idx = {n: i for i, n in enumerate(names)}
    c.ok(all(abs(e["t"] - col[idx[e["n"]]]) < 0.5 for e in top), f"tenx_visium: gene totals = the .h5 sums ({[e['n'] for e in top]})")
    sq = cache / "squidpy_visium_full" / "expr" / "spots.zarr"
    if sq.exists():
        g2 = zarr.open_group(str(sq), mode="r")
        bc2 = _strings(g2["barcodes"])
        xy2 = np.asarray(g2["xy"][...], np.float64) / sf["tissue_hires_scalef"]
        common = [i for i, b in enumerate(bc2) if b in pos]
        dd = float(np.abs(xy2[common] - np.array([pos[bc2[i]][:2] for i in common])).max())
        # Squidpy's file comes from another Space Ranger run of this sample: positions differ by 0-2 px
        c.ok(len(common) == len(bc2) and dd < 2.5, f"tenx_visium: the AnnData path (Squidpy's .h5ad) puts the same spots "
                                                    f"within {dd:.2f} full-res px")
        m2 = json.loads((cache / "squidpy_visium_full" / "manifest.json").read_text(encoding="utf-8"))
        c.ok(abs(m2["pixelSizeUm"] * sf["tissue_hires_scalef"] - m["pixelSizeUm"]) / m["pixelSizeUm"] < 1e-4,
             "tenx_visium: the AnnData path's pixel size (hires) x tissue_hires_scalef = the Visium path's")


def check_idr(c: Checker, cache: Path, a: str, b: str, src_a: Path, src_b: Path):
    import zarr
    import zarrv3
    for lev in ("0", "1", "2"):
        A = zarr.open_array(str(src_a / lev), mode="r")[...]
        B = zarrv3.open_array(src_b / lev)[...]
        c.ok(A.shape == B.shape and A.dtype == B.dtype and np.array_equal(A, B),
             f"IDR 6001240 level {lev}: OME-NGFF 0.5 (zarr v3, sharded, read by zarrv3.py) = 0.4 (zarr v2), every pixel {A.shape}")
    ga, gb = zarr.open_group(str(cache / a / "image.zarr"), mode="r"), zarr.open_group(str(cache / b / "image.zarr"), mode="r")
    c.ok(sorted(ga.array_keys()) == sorted(gb.array_keys()) and all(np.array_equal(ga[k][...], gb[k][...]) for k in ga.array_keys()),
         f"{a} / {b}: the two prepared pyramids are equal")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--root", default=str(DEFAULT_ROOT))
    ap.add_argument("--work", default=None)
    ap.add_argument("--only", default="")
    ap.add_argument("--large", action="store_true", help="also the Xenium rep1 store (3.6 GB; about 4 min prepare + 7 min export)")
    ap.add_argument("--reuse", action="store_true", help="keep earlier cache / export (prepare and export are not run again)")
    a = ap.parse_args()
    root = Path(a.root)
    data = root / "data"
    if not data.is_dir():
        print(f"SKIPPED: no real data at {data} (see the top of this file)")
        return 0
    work = Path(a.work) if a.work else root / "realtest"
    cache, export, shots = work / "cache", work / "export", work / "shots"
    log = work / "log.txt"
    work.mkdir(parents=True, exist_ok=True)
    only = [x for x in a.only.split(",") if x]
    import detect
    c = Checker()
    done, times = [], {}
    order = sorted(DATASETS, key=lambda k: (k == "tenx_visium", k))      # the AnnData copy of the Visium section first
    for sid in order:
        folder, rel, kind, cmd, opts = DATASETS[sid]
        if only and sid not in only:
            continue
        if opts.get("large") and not a.large and sid not in only:
            print(f"{sid}: left out (--large)")
            continue
        d = data / folder
        if not (d / ".done").exists() and not (d / "extracted").exists() and not any(d.glob("*")):
            print(f"{sid}: not downloaded, skipped")
            continue
        src = d / rel
        if not src.exists():
            print(f"{sid}: {src} not there (not extracted?), skipped")
            continue
        print(sid, flush=True)
        r = detect.detect(d)
        if opts.get("refused"):
            c.ok(r["verdict"] == "incomplete" and r["kind"] == kind and all(w in r["message"] for w in opts["refused"]),
                 f"{sid}: detect: incomplete, {kind}, says {opts['refused']} (got {r['verdict']} {r['kind']}: {r['message'][:100]})")
        else:
            c.ok(r["verdict"] == "certain" and r["kind"] == kind, f"{sid}: detect: certain, {kind} (got {r['verdict']} {r['kind']}: {r['message'][:100]})")
        if opts.get("refused") or not a.reuse or not (export / sid / "manifest.json").exists():
            for p in (cache / sid, export / sid):
                if p.exists():
                    shutil.rmtree(p)
            if kind == "visium":
                img = next(d.glob("*_image.tif"), None) or next(d.glob("*.tif"))
                args = ["visium", "--image", str(img), "--outs", str(d)]
            else:
                args = [cmd[0], "--src", str(src), *cmd[1:]]
            rc, out, t = prepare(args + ["--id", sid, "--out", str(cache)], log)
            times[f"{sid} prepare"] = t
            if opts.get("refused"):                      # SODA does not read it: a clear refusal, nothing prepared
                c.ok(rc != 0 and all(w in out for w in opts["refused"]) and "Traceback" not in out
                     and not (cache / sid / "manifest.json").exists(),
                     f"{sid}: prepare refuses it with {opts['refused']}, no traceback, nothing prepared "
                     f"(exit {rc}: {out.strip().splitlines()[-1][:200] if out.strip() else ''})")
                continue
            if not c.ok(rc == 0 and "Traceback" not in out, f"{sid}: prepare ({t:.0f} s){'' if rc == 0 else ': ' + out[-300:]}"):
                continue
            rc, out, t = prepare(["export", "--id", sid, "--out", str(cache), "--dest", str(export), "--fmt", "jpg", "--no-umap"], log)
            times[f"{sid} export"] = t
            if not c.ok(rc == 0 and "Traceback" not in out, f"{sid}: export ({t:.0f} s){'' if rc == 0 else ': ' + out[-300:]}"):
                continue
        if "sdata" in opts:
            check_sdata(c, sid, src, cache, opts["sdata"], opts)
        elif opts.get("pixel"):
            m = json.loads((cache / sid / "manifest.json").read_text(encoding="utf-8"))
            c.ok(abs((m.get("pixelSizeUm") or 0) - opts["pixel"]) < 1e-3, f"{sid}: pixel size {m.get('pixelSizeUm')} µm from the level-0 scale")
        if opts.get("same_as") and (cache / opts["same_as"] / "manifest.json").exists():
            other = DATASETS[opts["same_as"]]
            check_idr(c, cache, opts["same_as"], sid, data / other[0] / other[1], src)
        if opts.get("tenx"):
            check_tenx(c, d, cache, export, work)
        p = subprocess.run([sys.executable, str(APP / "deid.py"), "check", str(export / sid)], cwd=str(REPO), stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace", env={**os.environ, "PYTHONIOENCODING": "utf-8"})
        c.ok(p.returncode == 0, f"{sid}: deid.py check on the export: {p.stdout.strip().splitlines()[-1] if p.stdout.strip() else p.returncode}")
        done.append(sid)
    node = shutil.which("node")
    if done and node:
        p = subprocess.run([node, str(REPO / "tests" / "stdformats_real_browser.js"), "--export", str(export), "--ids", ",".join(done),
                            "--shots", str(shots)], cwd=str(REPO), stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                           text=True, encoding="utf-8", errors="replace")
        print("\n".join(ln for ln in p.stdout.splitlines() if not ln.startswith("RESULT")))
        c.ok(p.returncode == 0 or "SKIPPED" in p.stdout, "the real exports open in Chromium (tests/stdformats_real_browser.js)")
    elif done:
        print("node not found: the browser part is skipped")
    for k, v in times.items():
        print(f"  {k}: {v:.0f} s")
    print(f"\n{'FAILED' if c.failed else 'PASSED'}: {c.n - len(c.failed)} of {c.n} checks ({len(done)} datasets)")
    for f in c.failed:
        print("  FAILED", f)
    return 1 if c.failed else 0


if __name__ == "__main__":
    sys.exit(main())
