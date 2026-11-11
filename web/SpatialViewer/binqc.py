"""Per-bin QC values for an expression layer: total counts, genes detected, % mitochondrial.

    expr/<key>/binqc/total.bin     float32 per row   UMI over all genes of the layer
    expr/<key>/binqc/genes.bin     uint16 per row    genes with a count > 0 (uint32 when the
                                                     layer has 65,536 genes or more)
    expr/<key>/binqc/pct_mt.bin    float32 per row   100 x UMI of MT- genes / total (0 when
                                                     total is 0); only when MT- genes exist

A row is a position in positions.bin (bins with a position, in count-matrix order), the
same row space as region.bin, hvg/ and umap/. No file carries a barcode. The manifest
layer gets a `binQc` block naming the files with each metric's range and quartiles; the
viewer draws a histogram per metric with min / max sliders and leaves the bins outside
them out of the circle -> UMAP tool. Exports made before this have no `binQc` block and
the panel stays hidden.

The total equals the hvg/ chunks' library size (both are the row sums of the same matrix),
so filtering on it and normalising by it agree.
"""
from __future__ import annotations

import re

import numpy as np

import deid

# Human mitochondrial genes are MT-ND1, MT-CO1, ...; mouse ones mt-Nd1, ... The match
# ignores case so either works. The layer's own gene list decides, never a fixed list.
MT_PATTERN = re.compile(r"^MT-", re.IGNORECASE)
DIR = "binqc"


def _stats(v: np.ndarray) -> dict:
    if not len(v):
        return {"min": 0.0, "max": 0.0, "median": 0.0, "q1": 0.0, "q3": 0.0}
    q1, med, q3 = np.percentile(v.astype(np.float64), [25, 50, 75])
    r = lambda x: float(round(float(x), 4))
    return {"min": r(v.min()), "max": r(v.max()), "median": r(med), "q1": r(q1), "q3": r(q3)}


def compute(er) -> dict:
    """Total counts, genes detected and % MT per row of positions.bin, from the cache's
    gene-major CSR (ExprReader `er`). No transpose: one bincount per metric."""
    xy = er.xy
    B = len(xy)
    on_image = ~np.isnan(xy[:, 0])
    indices = np.asarray(er.g["indices"])
    data = np.asarray(er.g["data"], dtype=np.float32)
    indptr = er.indptr
    total = np.bincount(indices, weights=data.astype(np.float64), minlength=B)
    ngenes = np.bincount(indices[data > 0], minlength=B)
    mt_idx = [i for i, g in enumerate(er.genes) if MT_PATTERN.match(str(g))]
    out = {"total": total[on_image].astype(np.float32),
           "genes": ngenes[on_image],
           "mtGenes": [er.genes[i] for i in mt_idx]}
    if mt_idx:
        mt = np.zeros(B, np.float64)
        for gi in mt_idx:
            a, b = int(indptr[gi]), int(indptr[gi + 1])
            if b > a:
                np.add.at(mt, indices[a:b], data[a:b])
        tot = total[on_image]
        mt = mt[on_image]
        out["pctMt"] = np.where(tot > 0, 100.0 * mt / np.maximum(tot, 1e-12), 0.0).astype(np.float32)
    return out


def write_bin_qc(er, d, log=print) -> dict | None:
    """Write expr/<key>/binqc/ for one layer; returns the manifest `binQc` block."""
    q = compute(er)
    n = len(q["total"])
    if not n:
        return None
    hd = d / DIR
    hd.mkdir(parents=True, exist_ok=True)
    gdt = "<u2" if len(er.genes) < 65536 else "<u4"
    metrics = []
    (hd / "total.bin").write_bytes(np.ascontiguousarray(q["total"], dtype="<f4").tobytes())
    metrics.append({"key": "total", "label": "Total counts", "file": "total.bin", "dtype": "float32",
                    "unit": "UMI", "scale": "log", **_stats(q["total"])})
    (hd / "genes.bin").write_bytes(np.ascontiguousarray(q["genes"], dtype=gdt).tobytes())
    metrics.append({"key": "genes", "label": "Genes detected", "file": "genes.bin",
                    "dtype": "uint16" if gdt == "<u2" else "uint32", "unit": "genes", "scale": "log",
                    **_stats(q["genes"])})
    if "pctMt" in q:
        (hd / "pct_mt.bin").write_bytes(np.ascontiguousarray(q["pctMt"], dtype="<f4").tobytes())
        metrics.append({"key": "pctMt", "label": "% mitochondrial", "file": "pct_mt.bin", "dtype": "float32",
                        "unit": "%", "scale": "linear", **_stats(q["pctMt"])})
    block = {"dir": DIR, "nRows": int(n), "metrics": metrics,
             "mtPattern": MT_PATTERN.pattern + " (ignoring case)", "mtGenes": q["mtGenes"]}
    deid.assert_clean(block, f"{DIR}/ manifest block")
    mt_txt = (f", {len(q['mtGenes'])} MT- genes, median {metrics[2]['median']:.2f}% MT"
              if "pctMt" in q else ", no MT- genes in this layer")
    log(f"[binqc] {er.meta.get('key')}: {n:,} bins, median {metrics[0]['median']:g} UMI, "
        f"{metrics[1]['median']:g} genes{mt_txt}")
    return block
