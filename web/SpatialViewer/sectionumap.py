"""Whole-section UMAP + Leiden clusters, made from the same bins the viewer draws.

The circle tool embeds the bins under a shape in the browser. This embeds every
bin of a layer once, at export, so the viewer can flip the main view between the
tissue and the UMAP and colour both by cluster.

    expr/<key>/umap/index.json            {nBins, nEmbedded, minUmi, bounds, params, versions,
                                           clusterings: [{resolution, file, n, sizes, markers}], default}
    expr/<key>/umap/coords.bin            float32 x,y per positions.bin row; NaN,NaN = not embedded
    expr/<key>/umap/clusters_r<res>.bin   uint16 per positions.bin row; 0xFFFF = not embedded

Pipeline (scanpy's defaults, written out): bins with >= MIN_UMI UMI -> log1p CP10k
-> top 2000 HVGs over the whole section (hvgchunks.select_hvg) -> scale, clip 10
-> PCA 30 -> umap-learn (15 neighbours, min_dist 0.3) -> Leiden (leidenalg,
RB configuration) on umap's fuzzy graph at 0.3 / 0.6 / 1.0. Cluster 0 is the
largest. Markers: top 10 genes per cluster by log2 fold change of the mean
expression (cluster vs rest) among genes seen in >= 25% of the cluster's bins.

Everything is seeded, so the same cache gives the same bytes. umap-learn, igraph
and leidenalg are imported lazily: without them the step is skipped and the rest
of the export carries on.
"""
from __future__ import annotations

import json
import sys
import time

import numpy as np
from scipy import linalg, sparse

import deid
from hvgchunks import _BARCODE, MIN_UMI, select_hvg

N_HVG = 2000
N_PCS = 30
N_NEIGHBORS = 15
MIN_DIST = 0.3
SEED = 42
CLIP = 10.0
RESOLUTIONS = (0.3, 0.6, 1.0)
DEFAULT_RES = 1                                         # index into RESOLUTIONS
N_MARKERS = 10
MIN_PCT = 0.25
NO_CLUSTER = 0xFFFF
# scanpy's default. Iterating to convergence (-1) took 15-35x longer on 16 um bins
# for under 1% more modularity and the same cluster count.
LEIDEN_ITERS = 2


def _params() -> dict:
    return {"nHvg": N_HVG, "nPcs": N_PCS, "nNeighbors": N_NEIGHBORS, "minDist": MIN_DIST, "seed": SEED,
            "normalize": "log1p CP10k", "scale": f"clip {CLIP:g}",
            "cluster": f"leiden RBConfiguration on umap fuzzy graph, {LEIDEN_ITERS} iterations",
            "markers": f"log2 fold change of mean expression vs rest, pct >= {MIN_PCT:g}, lfc > 0"}


def _pca(Xh: sparse.csr_matrix, n_pcs: int) -> np.ndarray:
    """PCA of the z-scored, clipped HVG matrix without ever holding it dense.

    After scaling and clipping, every zero of a gene column becomes the same value
    z0 = max(-mu/sd, -10), so the scaled matrix is Z = 1 z0' + S with S sparse
    (the clipped nonzeros minus z0). Its covariance is then S'S plus rank-one
    terms: exact, float64, a 2000 x 2000 eigenproblem, and deterministic - no
    randomised solver. Signs are fixed so each component's largest loading is positive.
    """
    n, p = Xh.shape
    Xh = Xh.astype(np.float64).tocsr()
    mu = np.asarray(Xh.mean(0)).ravel()
    sq = np.asarray(Xh.multiply(Xh).mean(0)).ravel()
    sd = np.sqrt(np.maximum(sq - mu * mu, 0) * n / max(n - 1, 1))   # ddof 1, as scanpy
    sd[sd == 0] = 1.0
    z0 = np.maximum(-mu / sd, -CLIP)
    S = Xh.copy()
    col = S.indices
    S.data = np.clip((S.data - mu[col]) / sd[col], -CLIP, CLIP) - z0[col]
    cs = np.asarray(S.sum(0)).ravel()
    m = z0 + cs / n                                    # column means of Z (clipping moves them off zero)
    StS = (S.T @ S).toarray()
    ZtZ = StS + np.outer(cs, z0) + np.outer(z0, cs) + n * np.outer(z0, z0)
    C = (ZtZ - n * np.outer(m, m)) / max(n - 1, 1)
    k = min(n_pcs, p)
    w, V = linalg.eigh(C, subset_by_index=[p - k, p - 1])
    V = V[:, np.argsort(-w, kind="stable")]
    V *= np.where(V[np.abs(V).argmax(0), np.arange(V.shape[1])] < 0, -1.0, 1.0)
    # (Z - 1 m') V = S V + 1 (z0 - m)' V
    return (np.asarray(S @ V) + (z0 - m) @ V).astype(np.float32)


def _import_umap():
    """umap-learn, without the TensorFlow its package __init__ drags in for ParametricUMAP.

    Loading TensorFlow costs seconds and brings its own thread pools, which fought
    OpenBLAS for the cores (a dense eigensolve ran five times slower after it).
    ParametricUMAP is not used.
    """
    block = "tensorflow" not in sys.modules
    if block:
        sys.modules["tensorflow"] = None               # import raises ImportError -> umap skips it
    try:
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ImportWarning)
            import umap
    finally:
        if block and sys.modules.get("tensorflow", 0) is None:
            del sys.modules["tensorflow"]
    return umap


def _graph(graph: sparse.spmatrix, ig):
    g = sparse.triu(sparse.csr_matrix(graph), k=1).tocoo()
    G = ig.Graph(n=graph.shape[0], edges=np.column_stack([g.row, g.col]).tolist(), directed=False)
    G.es["weight"] = g.data.astype(np.float64).tolist()
    return G


def _leiden(G, resolution: float, la) -> np.ndarray:
    part = la.find_partition(G, la.RBConfigurationVertexPartition, weights="weight",
                             resolution_parameter=resolution, n_iterations=LEIDEN_ITERS, seed=SEED)
    lab = np.asarray(part.membership, dtype=np.int64)
    # renumber: 0 = largest, ties broken by the old id so the order is stable
    sizes = np.bincount(lab)
    order = np.argsort(-sizes, kind="stable")
    remap = np.empty_like(order)
    remap[order] = np.arange(len(order))
    return remap[lab]


def _markers(Xn: sparse.csr_matrix, lab: np.ndarray, k: int, genes: list[str]) -> list[list[dict]]:
    """Top N_MARKERS per cluster: log2 fold change of mean CP10k (cluster vs rest), pct >= MIN_PCT."""
    n = len(lab)
    M = sparse.csr_matrix((np.ones(n, np.float64), (lab, np.arange(n))), shape=(k, n))
    nk = np.bincount(lab, minlength=k).astype(np.float64)
    tot = np.asarray(Xn.sum(0), dtype=np.float64).ravel()
    s = np.asarray((M @ Xn).todense(), dtype=np.float64)                # k x G sums of log1p CP10k
    nz = Xn.copy()
    nz.data = np.ones_like(nz.data)
    c = np.asarray((M @ nz).todense(), dtype=np.float64)               # k x G bins expressing
    out = []
    for i in range(k):
        rest = max(n - nk[i], 1.0)
        mi = s[i] / nk[i]
        mo = (tot - s[i]) / rest
        # the scanpy convention: fold change of the means back on the CP10k scale
        lfc = np.log2((np.expm1(mi) + 1e-9) / (np.expm1(mo) + 1e-9))
        pct = c[i] / nk[i]
        lfc[pct < MIN_PCT] = -np.inf
        top = np.argsort(-lfc, kind="stable")[:N_MARKERS]
        # only genes up in the cluster: a low-count cluster can have none, and then lists none
        out.append([{"gene": genes[j], "lfc": round(float(lfc[j]), 2), "pct": round(float(pct[j]), 2)}
                    for j in top if np.isfinite(lfc[j]) and lfc[j] > 0])
    return out


def write_section_umap(er, d, log=print, force: bool = False) -> dict | None:
    """Write expr/<key>/umap/ for one layer (ExprReader `er`, output dir `d`).

    Returns the manifest entry, or None when skipped. An index.json from an earlier
    run with the same params is reused unless `force`, like the rest of the export.
    """
    key = er.meta.get("key")
    ud = d / "umap"
    if not force and (ud / "index.json").exists():
        try:
            old = json.loads((ud / "index.json").read_text(encoding="utf-8"))
            if old.get("params") == _params() and all(
                    (ud / f).exists() for f in ["coords.bin"] + [c["file"] for c in old["clusterings"]]):
                log(f"[umap] {key}: reusing {ud / 'index.json'}")
                return _entry(old)
        except (ValueError, KeyError):
            pass
    try:
        import igraph as ig
        import leidenalg as la
        umap = _import_umap()
    except ImportError as e:
        log(f"[umap] {key}: skipped - {e.name} is not installed (pip install umap-learn igraph leidenalg)")
        return None

    t0 = time.time()
    xy = er.xy
    on_image = ~np.isnan(xy[:, 0])
    # positions.bin keeps only bins with a position, in order: that is the row space
    row_of = np.cumsum(on_image) - 1
    n_rows = int(on_image.sum())
    G, B = len(er.genes), len(xy)
    gm = sparse.csr_matrix((np.asarray(er.g["data"], dtype=np.float32), np.asarray(er.g["indices"]),
                            er.indptr), shape=(G, B))
    X = gm.T.tocsr()                                   # bins x genes
    del gm
    X.sort_indices()
    lib = np.asarray(X.sum(1)).ravel().astype(np.float32)
    ok = on_image & (lib >= MIN_UMI)
    if ok.sum() < 50:
        log(f"[umap] {key}: only {int(ok.sum())} bins with >= {MIN_UMI} UMI, skipped")
        return None
    cols = select_hvg(X, lib, ok, N_HVG)
    bins = np.flatnonzero(ok)
    Xn = (sparse.diags(1e4 / lib[bins]) @ X[bins]).tocsr().astype(np.float32)
    Xn.data = np.log1p(Xn.data)
    del X
    n = len(bins)
    log(f"[umap] {key}: {n:,} bins >= {MIN_UMI} UMI, {len(cols)} HVGs  {time.time() - t0:.0f}s")

    pcs = _pca(Xn[:, cols].tocsr(), N_PCS)
    log(f"[umap] {key}: PCA {pcs.shape[1]}  {time.time() - t0:.0f}s")
    # random_state makes umap-learn single-threaded: slower, but the same bytes every run
    reducer = umap.UMAP(n_neighbors=N_NEIGHBORS, min_dist=MIN_DIST, n_components=2, metric="euclidean",
                        random_state=SEED, transform_seed=SEED)
    emb = reducer.fit_transform(pcs).astype(np.float32)
    log(f"[umap] {key}: UMAP  {time.time() - t0:.0f}s")

    clusterings, labels = [], []
    G = _graph(reducer.graph_, ig)
    for res in RESOLUTIONS:
        lab = _leiden(G, res, la)
        k = int(lab.max()) + 1
        sizes = np.bincount(lab, minlength=k)
        clusterings.append({"resolution": res, "file": f"clusters_r{res:.1f}.bin", "n": k,
                            "sizes": [int(v) for v in sizes],
                            "markers": _markers(Xn, lab, k, er.genes)})
        labels.append(lab)
        log(f"[umap] {key}: leiden r={res:g} -> {k} clusters  {time.time() - t0:.0f}s")

    index = {"nBins": n_rows, "nEmbedded": n,
             "minUmi": MIN_UMI,
             "bounds": [round(float(v), 4) for v in (emb[:, 0].min(), emb[:, 1].min(),
                                                     emb[:, 0].max(), emb[:, 1].max())],
             "params": _params(),
             "versions": {"umap-learn": getattr(umap, "__version__", "?"),
                          "leidenalg": getattr(la, "__version__", getattr(la, "version", "?")),
                          "igraph": getattr(ig, "__version__", "?"), "numpy": np.__version__},
             "clusterings": clusterings, "default": DEFAULT_RES}
    deid.assert_clean(index, "umap/index.json")
    txt = json.dumps(index, separators=(",", ":"))
    if _BARCODE.search(txt):
        raise SystemExit("[umap] REFUSING to write index.json: it contains a barcode")

    ud.mkdir(parents=True, exist_ok=True)
    rows = row_of[bins]
    coords = np.full((n_rows, 2), np.nan, np.float32)
    coords[rows] = emb
    files = {"coords.bin": np.ascontiguousarray(coords, dtype="<f4").tobytes()}
    for c, lab in zip(clusterings, labels):
        cl = np.full(n_rows, NO_CLUSTER, np.uint16)
        cl[rows] = lab
        files[c["file"]] = cl.astype("<u2").tobytes()
    files["index.json"] = txt.encode("utf-8")           # last: its presence means the folder is complete
    total = 0
    for name, buf in files.items():
        tmp = ud / (name + ".part")
        tmp.write_bytes(buf)
        tmp.replace(ud / name)
        total += len(buf)
    log(f"[umap] {key}: {n:,}/{n_rows:,} bins embedded, clusters "
        + ", ".join(f"r{c['resolution']:g}={c['n']}" for c in clusterings)
        + f", {total / 1e6:.1f} MB  {time.time() - t0:.0f}s")
    return _entry(index)


def _entry(index: dict) -> dict:
    """The manifest's view of an index.json."""
    return {"nEmbedded": index["nEmbedded"],
            "clusterings": [{"resolution": c["resolution"], "n": c["n"]} for c in index["clusterings"]],
            "default": index["default"]}
