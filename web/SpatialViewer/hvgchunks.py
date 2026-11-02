"""Circle -> UMAP: a bin-major table of the most variable genes, cut into spatial chunks.

The per-gene export answers "where is gene X", one file per gene. A UMAP of the bins
under a drawn shape needs the opposite - many genes for a few bins - so this writes,
per expression layer, the top N highly variable genes bin by bin, in squares of
tissue. The viewer fetches only the squares its shapes touch.

    expr/<key>/hvg/genes.json        ["Gene", ...]   column order of every chunk
    expr/<key>/hvg/index.json        {chunkPx, chunkUm, nGenes, minUmi, chunks: [{tx, ty, n, bytes}]}
    expr/<key>/hvg/<tx>_<ty>.bin     little-endian, one chunk:
        u32 nBins
        u32 row[nBins]          index into positions.bin (never a barcode)
        f32 libSize[nBins]      UMI over ALL genes, so CP10k normalises honestly
        u32 indptr[nBins + 1]
        u16 gene[nnz]           column in hvg/genes.json
        u16 count[nnz]          UMI count, clipped at 65535

Chunk of a bin = floor(x / chunkPx), floor(y / chunkPx) on the positions.bin
float32 values, with chunkPx an integer so the browser reproduces it exactly.

HVGs are chosen once over the whole section (log1p CP10k, Seurat-style dispersion
z-scored within 20 mean bins), on bins with at least `min_umi` UMI.
"""
from __future__ import annotations

import json
import re

import numpy as np
from scipy import sparse

import deid

DEFAULT_N_HVG = 2000
DEFAULT_CHUNK_UM = 500.0
MIN_UMI = 10
# Space Ranger barcodes: Visium HD bins and classic spot barcodes. None may be written.
_BARCODE = re.compile(r"s_\d{3}um_\d{5}_\d{5}|[ACGT]{16}-1")


def select_hvg(X: sparse.csr_matrix, lib: np.ndarray, ok: np.ndarray, n_top: int) -> np.ndarray:
    """Column indices of the top `n_top` genes by mean-binned dispersion z-score."""
    Xo = X[ok]
    Xn = sparse.diags(1e4 / lib[ok]) @ Xo
    Xn.data = np.log1p(Xn.data)
    mu = np.asarray(Xn.mean(0)).ravel()
    var = np.asarray(Xn.multiply(Xn).mean(0)).ravel() - mu * mu
    disp = np.log(np.maximum(var, 1e-12) / np.maximum(mu, 1e-12))
    expressed = mu > 1e-3
    z = np.full(X.shape[1], -np.inf)
    if not expressed.any():
        return np.zeros(0, np.int64)
    lm = np.log1p(mu)
    edges = np.linspace(lm[expressed].min(), lm[expressed].max(), 21)[1:-1]
    bins = np.digitize(lm, edges)
    for b in np.unique(bins[expressed]):
        m = expressed & (bins == b)
        s = disp[m].std()
        z[m] = (disp[m] - disp[m].mean()) / (s if s > 0 else 1.0)
    order = np.argsort(-z, kind="stable")
    order = order[np.isfinite(z[order])][:n_top]
    return np.sort(order)


def write_hvg_chunks(er, d, mpp: float | None, n_top: int = DEFAULT_N_HVG,
                     chunk_um: float = DEFAULT_CHUNK_UM, log=print) -> dict | None:
    """Write expr/<key>/hvg/ for one layer (ExprReader `er`, output dir `d`).

    Returns the manifest entry, or None when there is nothing to write.
    """
    xy = er.xy
    on_image = ~np.isnan(xy[:, 0])
    # positions.bin keeps only bins with a position, in order: that is the row space
    row_of = np.cumsum(on_image) - 1
    G, B = len(er.genes), len(xy)
    gm = sparse.csr_matrix((np.asarray(er.g["data"], dtype=np.float32), np.asarray(er.g["indices"]),
                            er.indptr), shape=(G, B))
    X = gm.T.tocsr()                                   # bins x genes
    del gm
    X.sort_indices()
    lib = np.asarray(X.sum(1)).ravel().astype(np.float32)
    ok = on_image & (lib >= MIN_UMI)
    if ok.sum() < 50:
        log(f"[hvg] {er.meta.get('key')}: only {int(ok.sum())} bins with >= {MIN_UMI} UMI, skipped")
        return None
    cols = select_hvg(X, lib, ok, n_top)
    genes = [er.genes[i] for i in cols]
    deid.assert_clean(genes, "hvg/genes.json")

    Xs = X[:, cols].tocsr()
    Xs.sort_indices()
    chunk_px = int(round(chunk_um / mpp)) if mpp else 1024
    keep = np.where(on_image)[0]                       # every bin on the image, low-UMI included
    xy64 = xy[keep].astype(np.float64)
    tx = np.floor(xy64[:, 0] / chunk_px).astype(np.int64)
    ty = np.floor(xy64[:, 1] / chunk_px).astype(np.int64)
    # a few bins can sit just off the image edge: chunk names then carry a minus sign

    hd = d / "hvg"
    hd.mkdir(parents=True, exist_ok=True)
    order = np.lexsort((keep, ty, tx))
    keys = np.stack([tx[order], ty[order]], 1)
    split = np.flatnonzero(np.any(np.diff(keys, axis=0) != 0, axis=1)) + 1
    chunks, total = [], 0
    for grp in np.split(order, split):
        bins = keep[grp]
        sub = Xs[bins]
        counts = np.minimum(np.rint(sub.data), 65535).astype("<u2")
        buf = b"".join([
            np.uint32(len(bins)).astype("<u4").tobytes(),
            row_of[bins].astype("<u4").tobytes(),
            lib[bins].astype("<f4").tobytes(),
            sub.indptr.astype("<u4").tobytes(),
            sub.indices.astype("<u2").tobytes(),
            counts.tobytes(),
        ])
        cx, cy = int(tx[grp[0]]), int(ty[grp[0]])
        (hd / f"{cx}_{cy}.bin").write_bytes(buf)
        chunks.append({"tx": cx, "ty": cy, "n": int(len(bins)), "bytes": len(buf)})
        total += len(buf)
    index = {"chunkPx": chunk_px, "chunkUm": chunk_um, "nGenes": len(genes), "minUmi": MIN_UMI,
             "selection": "seurat dispersion, log1p CP10k, whole section", "chunks": chunks}
    deid.assert_clean(index, "hvg/index.json")
    for name, obj in (("genes.json", genes), ("index.json", index)):
        txt = json.dumps(obj, separators=(",", ":"))
        if _BARCODE.search(txt):
            raise SystemExit(f"[hvg] REFUSING to write {name}: it contains a barcode")
        (hd / name).write_text(txt, encoding="utf-8")
    log(f"[hvg] {er.meta.get('key')}: {len(genes)} genes, {len(chunks)} chunks of {chunk_um:g} um, "
        f"{total / 1e6:.1f} MB, nnz {Xs[keep].nnz:,}")
    return {"nGenes": len(genes), "chunkPx": chunk_px, "chunks": len(chunks), "bytes": total}
