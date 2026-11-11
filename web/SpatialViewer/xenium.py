"""Xenium Onboard Analysis output folder -> cells layer + molecules layer.

Reads, when present (each part is optional; a folder with stains only still works):

    experiment.xenium                       pixel_size (µm per pixel of the morphology images)
    cells.parquet | cells.csv.gz            cell_id, x_centroid, y_centroid (µm), cell_area
    cell_feature_matrix.h5 | cell_feature_matrix/ (barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz)
    transcripts.parquet | transcripts.zarr.zip | transcripts.zarr/

Xenium coordinates are µm in the frame of the morphology images, origin at the top-left
corner of the full-resolution image, so pixel = µm / pixel_size (the same transform Xenium
Explorer and FIVE Atlas use; no affine).

Cells become an ordinary expression layer (expr.build_cell_layer), so the whole-section
UMAP, Leiden clusters, suggested names, bin QC, circle -> UMAP and the volcano panel run on
them unchanged. Transcripts become the molecules layer (molecules.py): quality value >= 20
(Xenium's own "high quality" cut, used for its cell-feature matrix), control probes and
codewords left out.
"""
from __future__ import annotations

import ast
import gzip
import json
import time
from pathlib import Path

import numpy as np

QV_MIN = 20.0
# names of codewords / probes that are not genes, in every Xenium version seen so far
CONTROL_PREFIXES = ("negcontrolprobe", "negcontrolcodeword", "blank", "unassignedcodeword",
                    "deprecatedcodeword", "intergenic", "antisense", "genomiccontrol")
GENE_FEATURE_TYPES = ("gene expression",)


def is_control(name: str) -> bool:
    return str(name).strip().lower().startswith(CONTROL_PREFIXES)


def _first(folder: Path, *names):
    for n in names:
        p = folder / n
        if p.exists():
            return p
    return None


def find_inputs(folder: Path) -> dict:
    """The Xenium files present under `folder` (None for each one that is missing)."""
    folder = Path(folder)
    mtx_dir = folder / "cell_feature_matrix"
    return {
        "experiment": _first(folder, "experiment.xenium"),
        "cells": _first(folder, "cells.parquet", "cells.csv.gz", "cells.csv"),
        "matrix_h5": _first(folder, "cell_feature_matrix.h5"),
        "matrix_mtx": mtx_dir if (mtx_dir / "matrix.mtx.gz").exists() or (mtx_dir / "matrix.mtx").exists() else None,
        "transcripts": _first(folder, "transcripts.parquet", "transcripts.zarr.zip", "transcripts.zarr"),
        "cell_boundaries": _first(folder, "cell_boundaries.parquet", "cell_boundaries.csv.gz"),
    }


def pixel_size(folder: Path) -> float | None:
    p = Path(folder) / "experiment.xenium"
    if not p.exists():
        return None
    try:
        v = json.loads(p.read_text(encoding="utf-8")).get("pixel_size")
        return float(v) if v else None
    except (ValueError, TypeError):
        return None


# ---------------------------------------------------------------------------------- cells
def _table(path: Path, columns: list[str]) -> dict:
    """Columns of a CSV (optionally gzipped) or Parquet table as numpy / lists."""
    path = Path(path)
    if path.suffix == ".parquet":
        import pyarrow.parquet as pq
        have = set(pq.ParquetFile(str(path)).schema_arrow.names)
        t = pq.read_table(str(path), columns=[c for c in columns if c in have])
    else:
        import pyarrow.csv as pcsv
        t = pcsv.read_csv(str(path))
    out = {}
    for c in columns:
        if c in t.column_names:
            col = t.column(c)
            out[c] = col.to_pylist() if c == "cell_id" else col.to_numpy()
    return out


def _cell_ids(v) -> list[str]:
    return [x.decode() if isinstance(x, bytes) else str(x) for x in v]


def read_matrix(inp: dict):
    """(X cells x features CSR float32, barcodes, feature names, feature types)."""
    from scipy import io as sio, sparse
    if inp.get("matrix_h5"):
        import h5py
        with h5py.File(inp["matrix_h5"], "r") as f:
            m = f["matrix"]
            G, B = (int(v) for v in m["shape"][()])
            X = sparse.csc_matrix((m["data"][()].astype(np.float32), m["indices"][()], m["indptr"][()]),
                                  shape=(G, B)).T.tocsr()
            bc = _cell_ids(m["barcodes"][()])
            names = _cell_ids(m["features/name"][()])
            types = _cell_ids(m["features/feature_type"][()])
        return X, bc, names, types
    d = Path(inp["matrix_mtx"])

    def lines(stem):
        p = _first(d, stem + ".gz", stem)
        op = gzip.open if p.suffix == ".gz" else open
        with op(p, "rt", encoding="utf-8") as f:
            return [ln.rstrip("\n\r").split("\t") for ln in f if ln.strip()]
    bc = [r[0] for r in lines("barcodes.tsv")]
    feats = lines("features.tsv")
    names = [r[1] if len(r) > 1 else r[0] for r in feats]
    types = [r[2] if len(r) > 2 else "Gene Expression" for r in feats]
    M = sio.mmread(str(_first(d, "matrix.mtx.gz", "matrix.mtx")))     # features x cells
    X = sparse.csr_matrix(M.T, dtype=np.float32)
    return X, bc, names, types


def read_cells(inp: dict, log=print) -> dict | None:
    """Cell centroids (µm) in count-matrix order + the gene matrix, or None when either is missing."""
    if not inp.get("cells") or not (inp.get("matrix_h5") or inp.get("matrix_mtx")):
        missing = (([] if inp.get("cells") else ["cells.parquet / cells.csv.gz"]) +
                   ([] if (inp.get("matrix_h5") or inp.get("matrix_mtx")) else ["cell_feature_matrix"]))
        log(f"[xenium] no cells layer: {' and '.join(missing)} not found")
        return None
    t0 = time.time()
    X, bc, names, types = read_matrix(inp)
    keep = np.array([t.strip().lower() in GENE_FEATURE_TYPES and not is_control(n) for n, t in zip(names, types)])
    n_ctrl = int((~keep).sum())
    X = X[:, np.flatnonzero(keep)].tocsr()
    genes = [n for n, k in zip(names, keep) if k]
    tab = _table(inp["cells"], ["cell_id", "x_centroid", "y_centroid", "cell_area"])
    ids = _cell_ids(tab["cell_id"])
    row = {c: i for i, c in enumerate(ids)}
    xy = np.full((len(bc), 2), np.nan, np.float64)
    area = np.full(len(bc), np.nan, np.float64)
    hit = np.array([row.get(b, -1) for b in bc], np.int64)
    ok = hit >= 0
    xy[ok, 0] = np.asarray(tab["x_centroid"], np.float64)[hit[ok]]
    xy[ok, 1] = np.asarray(tab["y_centroid"], np.float64)[hit[ok]]
    if "cell_area" in tab:
        area[ok] = np.asarray(tab["cell_area"], np.float64)[hit[ok]]
    log(f"[xenium] cells: {len(bc):,} in the matrix, {len(ids):,} in {Path(inp['cells']).name}, "
        f"{int(ok.sum()):,} matched; {len(genes)} genes ({n_ctrl} control features left out), "
        f"nnz {X.nnz:,}  ({time.time() - t0:.0f}s)")
    if int((~ok).sum()):
        log(f"[xenium]   WARNING {int((~ok).sum())} matrix cells have no centroid: no position")
    med_area = float(np.nanmedian(area)) if np.isfinite(area).any() else None
    return {"X": X, "genes": genes, "xy_um": xy, "barcodes": list(bc), "nCellsTable": len(ids), "medianAreaUm2": med_area,
            "nControlFeatures": n_ctrl}


# ---------------------------------------------------------------------------- transcripts
def _attr_list(v):
    return v if isinstance(v, (list, dict)) else ast.literal_eval(v)


def _open_zarr(path: Path):
    import zarr
    p = Path(path)
    if p.is_dir():
        return zarr.open(str(p), mode="r"), None
    store = zarr.ZipStore(str(p), mode="r")
    return zarr.open(store, mode="r"), store


def zarr_genes(root) -> tuple[list[str], np.ndarray]:
    """Gene names indexed by gene_identity, and which of them are genes (not controls)."""
    names = [str(n) for n in _attr_list(root.attrs["gene_names"])]
    keep = np.array([not is_control(n) for n in names])
    if "gene_category" in root:                                   # (genes x flags), flag 0 = is_gene
        try:
            cols = _attr_list(root["gene_category"].attrs.get("column_names", []))
            cat = np.asarray(root["gene_category"][:])
            if cat.ndim == 2 and cat.shape[0] == len(names) and "is_gene" in cols:
                keep &= cat[:, cols.index("is_gene")].astype(bool)
        except Exception:
            pass
    return names, keep


def zarr_batches(root, qv_min: float, stats: dict):
    """(gene_identity, x µm, y µm) per level-0 grid chunk, quality value >= qv_min."""
    g0 = root["grids/0"]
    keys = list(g0.group_keys())
    for k in keys:
        c = g0[k]
        if "gene_identity" not in c:
            continue
        gid = np.asarray(c["gene_identity"][:]).reshape(-1)
        if not len(gid):
            continue
        loc = np.asarray(c["location"][:])
        q = np.asarray(c["quality_score"][:]).reshape(-1)
        ok = q >= qv_min
        if "valid" in c:
            ok &= np.asarray(c["valid"][:]).reshape(-1).astype(bool)
        stats["read"] += len(gid)
        stats["lowQv"] += int((q < qv_min).sum())
        yield gid[ok], loc[ok, 0], loc[ok, 1]


def parquet_batches(path: Path, names: list[str], qv_min: float, stats: dict):
    import pyarrow.parquet as pq
    pf = pq.ParquetFile(str(path))
    idx = {n: i for i, n in enumerate(names)}
    cols = ["feature_name", "x_location", "y_location", "qv"]
    for b in pf.iter_batches(batch_size=4_000_000, columns=cols):
        fn = b.column(0).to_pylist()
        x = b.column(1).to_numpy()
        y = b.column(2).to_numpy()
        q = b.column(3).to_numpy()
        fn = [f.decode() if isinstance(f, bytes) else str(f) for f in fn]
        gid = np.fromiter((idx.get(f, -1) for f in fn), np.int64, len(fn))
        ok = q >= qv_min
        stats["read"] += len(gid)
        stats["lowQv"] += int((~ok).sum())
        yield gid[ok], x[ok], y[ok]


def parquet_genes(path: Path) -> list[str]:
    import pyarrow.compute as pc
    import pyarrow.parquet as pq
    t = pq.read_table(str(path), columns=["feature_name"])
    u = pc.unique(t.column(0)).to_pylist()
    return sorted(f.decode() if isinstance(f, bytes) else str(f) for f in u)


def write_transcripts(inp: dict, ds_dir: Path, image_size, pixel_size_um: float, qv_min: float = QV_MIN,
                      log=print) -> dict | None:
    """molecules/ from transcripts.parquet or transcripts.zarr(.zip); returns the manifest block."""
    import molecules
    src = inp.get("transcripts")
    if not src:
        log("[xenium] no molecules layer: transcripts.parquet / transcripts.zarr.zip not found")
        return None
    t0 = time.time()
    stats = {"read": 0, "lowQv": 0}
    store = None
    try:
        if Path(src).suffix == ".parquet":
            names = parquet_genes(src)
            keep = np.array([not is_control(n) for n in names])
            batches = parquet_batches(src, names, qv_min, stats)
            kind = "transcripts.parquet"
        else:
            root, store = _open_zarr(src)
            names, keep = zarr_genes(root)
            batches = zarr_batches(root, qv_min, stats)
            kind = "transcripts.zarr"
        log(f"[xenium] transcripts from {kind}: {int(keep.sum())} genes, {int((~keep).sum())} control codewords "
            f"left out, quality >= {qv_min:g}")
        w = molecules.MoleculeWriter(names, image_size, pixel_size_um, keep=keep, log=log)
        for i, (g, x, y) in enumerate(batches):
            w.add(g, x, y)
            if (i + 1) % 200 == 0:
                log(f"[xenium]   {stats['read']:,} transcripts read  ({time.time() - t0:.0f}s)")
        block = w.finish(ds_dir, meta={"platform": "xenium", "source": kind, "qvMin": qv_min,
                                       "read": int(stats["read"]), "lowQv": int(stats["lowQv"])})
    finally:
        if store is not None:
            store.close()
    log(f"[xenium] transcripts: {stats['read']:,} read, {stats['lowQv']:,} below quality {qv_min:g}  "
        f"({time.time() - t0:.0f}s)")
    return block
