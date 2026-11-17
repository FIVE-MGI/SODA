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
-> PCA 30 -> umap-learn (15 neighbours, min_dist 0.3, started from PC1 / PC2 like the
circle tool's umap-js, which has no spectral start) -> Leiden (leidenalg,
RB configuration) on umap's fuzzy graph at 0.3 / 0.6 / 1.0. Cluster 0 is the
largest. Markers: top 10 genes per cluster by log2 fold change of the mean
expression (cluster vs rest) among genes seen in >= 25% of the cluster's bins.

Protein grid layers (protgrid.py) go through the same UMAP / Leiden / renumbering /
file writing (embed_and_cluster, write_umap_files) from a dense matrix of marker means
instead: write_feature_umap; their index.json says "modality": "protein". Above
UMAP_SAMPLE rows they are fitted on a seeded sample and the rest placed with umap's
transform (index.json "fitOn"; params "fit" says how).

Everything is seeded, so the same cache gives the same bytes. umap-learn, igraph
and leidenalg are imported lazily: without them the step is skipped and the rest
of the export carries on.
"""
from __future__ import annotations

import hashlib
import json
import sys
import time

import threadcap
threadcap.apply()                                      # before numpy, when this is the first module to load it

import numpy as np
from scipy import linalg, sparse

import deid
from hvgchunks import _BARCODE, MIN_UMI, select_hvg

N_HVG = 2000
N_PCS = 30
N_NEIGHBORS = 15
MIN_DIST = 0.3
INIT = "pca"                                            # umap-learn >= 0.5.4
SEED = 42
CLIP = 10.0
RESOLUTIONS = (0.3, 0.6, 1.0)
DEFAULT_RES = 1                                         # index into RESOLUTIONS
N_MARKERS = 10
MIN_PCT = 0.25
NO_CLUSTER = 0xFFFF
MIN_HVG_DETECTED = 3   # RNA bins with fewer variable genes detected are near-identical points: left out
# scanpy's default. Iterating to convergence (-1) took 15-35x longer on 16 um bins
# for under 1% more modularity and the same cluster count.
LEIDEN_ITERS = 2


def _params() -> dict:
    return {"nHvg": N_HVG, "nPcs": N_PCS, "nNeighbors": N_NEIGHBORS, "minDist": MIN_DIST, "init": INIT, "seed": SEED,
            "normalize": "log1p CP10k", "scale": f"clip {CLIP:g}",
            "cluster": f"leiden RBConfiguration on umap fuzzy graph, {LEIDEN_ITERS} iterations",
            "markers": f"log2 fold change of mean expression vs rest, pct >= {MIN_PCT:g}, lfc > 0"}


# The PCA's dense part is one 2000 x 2000 eigensolve: about a second on one thread. More threads
# save a fraction of that at best, and at worst cost minutes: with OpenBLAS 0.3.27 on 8 of 24
# CPUs it took 3-4 s on 1 or 4 threads, 4-17 s on 8 (more when the cores were busy), 42 s on 16
# and 114-126 s on 24. So it runs on one, whatever the process-wide cap is.
PCA_THREADS = 1


def _pca(Xh: sparse.csr_matrix, n_pcs: int) -> np.ndarray:
    """_pca_uncapped on PCA_THREADS BLAS / OpenMP thread(s) (threadcap.limit; without
    threadpoolctl the process-wide cap of threadcap.apply applies instead)."""
    with threadcap.limit(PCA_THREADS):
        return _pca_uncapped(Xh, n_pcs)


def _pca_uncapped(Xh: sparse.csr_matrix, n_pcs: int) -> np.ndarray:
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


def _renumber(lab: np.ndarray) -> np.ndarray:
    """0 = largest cluster, ties broken by the old id so the order is stable."""
    sizes = np.bincount(lab)
    order = np.argsort(-sizes, kind="stable")
    remap = np.empty_like(order)
    remap[order] = np.arange(len(order))
    return remap[lab]


def _leiden(G, resolution: float, la) -> np.ndarray:
    part = la.find_partition(G, la.RBConfigurationVertexPartition, weights="weight",
                             resolution_parameter=resolution, n_iterations=LEIDEN_ITERS, seed=SEED)
    return _renumber(np.asarray(part.membership, dtype=np.int64))


def sample_rows(n: int, size: int | None) -> np.ndarray | None:
    """Sorted, seeded random rows to fit the UMAP on, or None to fit on all `n`
    (size None = UMAP_SAMPLE; 0 or >= n = all)."""
    size = UMAP_SAMPLE if size is None else int(size)
    if size <= 0 or n <= size:
        return None
    return np.sort(np.random.default_rng(SEED).choice(n, size=size, replace=False))


def _sample_neighbours(reducer, X: np.ndarray, k: int) -> np.ndarray:
    """Indices (into the fitted sample) of the k nearest sampled rows of each row of X: umap's own
    search index when it has one (the one transform used), else an exact search."""
    idx = getattr(reducer, "_knn_search_index", None)
    if idx is not None:
        try:
            return np.asarray(idx.query(X, k=k, epsilon=0.12)[0], dtype=np.int64)
        except Exception:                               # older pynndescent: fall through
            pass
    from sklearn.neighbors import NearestNeighbors
    return NearestNeighbors(n_neighbors=k).fit(reducer._raw_data).kneighbors(X, return_distance=False)


def _majority(lab_s: np.ndarray, nn: np.ndarray) -> np.ndarray:
    """Per row, the most common label among its neighbours (lowest label on ties)."""
    k = int(lab_s.max()) + 1
    votes = np.zeros((len(nn), k), np.int32)
    rows = np.repeat(np.arange(len(nn)), nn.shape[1])
    np.add.at(votes, (rows, lab_s[nn].ravel()), 1)
    return votes.argmax(1).astype(np.int64)


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


def _clustering_libs(key, log):
    """(umap, igraph, leidenalg), or None with a message when one is not installed."""
    try:
        import igraph as ig
        import leidenalg as la
        umap = _import_umap()
    except ImportError as e:
        log(f"[umap] {key}: skipped - {e.name} is not installed (pip install umap-learn igraph leidenalg)")
        return None
    return umap, ig, la


def _reusable(ud, params: dict, key, log) -> dict | None:
    """The manifest entry of an earlier umap/ with the same params and all its files, else None."""
    if not (ud / "index.json").exists():
        return None
    try:
        old = json.loads((ud / "index.json").read_text(encoding="utf-8"))
        if old.get("params") == params and all(
                (ud / f).exists() for f in ["coords.bin"] + [c["file"] for c in old["clusterings"]]):
            log(f"[umap] {key}: reusing {ud / 'index.json'}")
            return _entry(old)
    except (ValueError, KeyError):
        pass
    return None


def embed_and_cluster(pcs: np.ndarray, markers_of, libs, key, log=print, t0: float | None = None,
                      sample: np.ndarray | None = None):
    """UMAP of `pcs` (rows x components) + Leiden at every resolution on umap's fuzzy graph.

    `markers_of(lab, k)` gives the marker lists of one clustering. Returns the 2-D
    embedding (float32) and, per resolution, the index.json entry and the labels.
    Shared by the RNA (write_section_umap) and protein (write_feature_umap) layers.

    `sample` (sorted rows, sample_rows): fit the UMAP and Leiden on those rows only; every
    other row is placed into the same map with umap's transform and takes the cluster most
    common among its SAMPLE_K nearest sampled rows (PCA space). Clusters are then renumbered
    by their size over all rows.
    """
    umap, ig, la = libs
    t0 = time.time() if t0 is None else t0
    # random_state makes umap-learn single-threaded: slower, but the same bytes every run
    reducer = umap.UMAP(n_neighbors=N_NEIGHBORS, min_dist=MIN_DIST, n_components=2, metric="euclidean",
                        init=INIT, random_state=SEED, transform_seed=SEED)
    rest = nn = None
    if sample is None:
        emb = reducer.fit_transform(pcs).astype(np.float32)
        log(f"[umap] {key}: UMAP  {time.time() - t0:.0f}s")
    else:
        emb = np.empty((len(pcs), 2), np.float32)
        emb[sample] = reducer.fit_transform(pcs[sample]).astype(np.float32)
        log(f"[umap] {key}: UMAP fitted on {len(sample):,} of {len(pcs):,} rows  {time.time() - t0:.0f}s")
        rest = np.setdiff1d(np.arange(len(pcs)), sample, assume_unique=True)
        emb[rest] = reducer.transform(pcs[rest]).astype(np.float32)
        nn = _sample_neighbours(reducer, pcs[rest], SAMPLE_K)
        # transform gives NaN for a row whose nearest sampled rows are all at distance 0 (e.g.
        # identical all-background squares): no graph edges to place it by. Those take the
        # median place of the same sampled neighbours that decide their cluster.
        bad = np.flatnonzero(~np.isfinite(emb[rest]).all(1))
        if len(bad):
            emb[rest[bad]] = np.median(emb[sample][nn[bad]], axis=1)
            log(f"[umap] {key}: {len(bad)} placed rows had no place from transform: put at their neighbours' median")
        log(f"[umap] {key}: {len(rest):,} rows placed with transform  {time.time() - t0:.0f}s")

    clusterings, labels = [], []
    G = _graph(reducer.graph_, ig)
    for res in RESOLUTIONS:
        lab = _leiden(G, res, la)
        if sample is not None:
            full = np.empty(len(pcs), np.int64)
            full[sample] = lab
            full[rest] = _majority(lab, nn)
            lab = _renumber(full)
        k = int(lab.max()) + 1
        sizes = np.bincount(lab, minlength=k)
        clusterings.append({"resolution": res, "name": f"Leiden {res:g}", "source": "soda",
                            "file": f"clusters_r{res:.1f}.bin", "n": k,
                            "sizes": [int(v) for v in sizes],
                            "markers": markers_of(lab, k)})
        labels.append(lab)
        log(f"[umap] {key}: leiden r={res:g} -> {k} clusters  {time.time() - t0:.0f}s")
    return emb, clusterings, labels


def write_umap_files(ud, head: dict, params: dict, emb: np.ndarray, clusterings: list, labels: list,
                     n_rows: int, rows: np.ndarray, libs, key, log=print, t0: float | None = None) -> dict:
    """index.json + coords.bin + clusters_r*.bin; returns the manifest entry.

    `head` = the first keys of index.json (nBins, nEmbedded, ...); `rows` = the
    positions.bin row of each row of `emb`.
    """
    umap, ig, la = libs
    t0 = time.time() if t0 is None else t0
    n = len(emb)
    fin = np.isfinite(emb).all(1)
    if not fin.all():                                  # should not happen (see embed_and_cluster)
        log(f"[umap] {key}: {int((~fin).sum())} rows without a place: written as not embedded")
        emb = emb.copy(); emb[~fin] = np.nan
        labels = [np.where(fin, lab, -1) for lab in labels]
        head = {**head, "nEmbedded": int(fin.sum())}
    if not fin.any():
        raise SystemExit(f"[umap] {key}: no row has a place in the UMAP")
    e = emb[fin]
    index = {**head,
             "bounds": [round(float(v), 4) for v in (e[:, 0].min(), e[:, 1].min(), e[:, 0].max(), e[:, 1].max())],
             "params": params,
             "versions": {"umap-learn": getattr(umap, "__version__", "?"),
                          "leidenalg": getattr(la, "__version__", getattr(la, "version", "?")),
                          "igraph": getattr(ig, "__version__", "?"), "numpy": np.__version__},
             "clusterings": clusterings, "default": DEFAULT_RES}
    deid.assert_clean(index, "umap/index.json")
    # allow_nan=False: Python would write NaN, which a browser's JSON.parse refuses
    txt = json.dumps(index, separators=(",", ":"), allow_nan=False)
    if _BARCODE.search(txt):
        raise SystemExit("[umap] REFUSING to write index.json: it contains a barcode")

    ud.mkdir(parents=True, exist_ok=True)
    coords = np.full((n_rows, 2), np.nan, np.float32)
    coords[rows] = emb
    files = {"coords.bin": np.ascontiguousarray(coords, dtype="<f4").tobytes()}
    for c, lab in zip(clusterings, labels):
        cl = np.full(n_rows, NO_CLUSTER, np.uint16)
        cl[rows] = np.where(lab >= 0, lab, NO_CLUSTER)
        files[c["file"]] = cl.astype("<u2").tobytes()
    files["index.json"] = txt.encode("utf-8")           # last: its presence means the folder is complete
    total = 0
    for name, buf in files.items():
        tmp = ud / (name + ".part")
        tmp.write_bytes(buf)
        tmp.replace(ud / name)
        total += len(buf)
    log(f"[umap] {key}: {n:,}/{n_rows:,} rows embedded, clusters "
        + ", ".join(f"r{c['resolution']:g}={c['n']}" for c in clusterings)
        + f", {total / 1e6:.1f} MB  {time.time() - t0:.0f}s")
    return _entry(index)


def write_section_umap(er, d, log=print, force: bool = False) -> dict | None:
    """Write expr/<key>/umap/ for one layer (ExprReader `er`, output dir `d`).

    Returns the manifest entry, or None when skipped. An index.json from an earlier
    run with the same params is reused unless `force`, like the rest of the export.
    """
    key = er.meta.get("key")
    ud = d / "umap"
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
    # a bin with (almost) no variable gene has no profile: such bins pile up on one point and
    # Leiden cuts that pile into clusters of nothing
    det = np.diff(X[:, cols].tocsr().indptr)
    few = ok & (det < MIN_HVG_DETECTED)
    ok &= ~few
    bins = np.flatnonzero(ok)
    # reuse an earlier umap/ only if it was made from these very bins and counts
    h = hashlib.sha256()
    for a in (bins.astype("<i8"), lib[bins].astype("<f4"), np.asarray(cols, "<i8"), det[bins].astype("<i8")):
        h.update(np.ascontiguousarray(a).tobytes())
    params = {**_params(), "minHvgDetected": MIN_HVG_DETECTED, "input": h.hexdigest()[:16]}
    if not force:
        old = _reusable(ud, params, key, log)
        if old:
            return old
    libs = _clustering_libs(key, log)
    if libs is None:
        return None
    Xn = (sparse.diags(1e4 / lib[bins]) @ X[bins]).tocsr().astype(np.float32)
    Xn.data = np.log1p(Xn.data)
    del X
    n = len(bins)
    log(f"[umap] {key}: {n:,} bins >= {MIN_UMI} UMI, {len(cols)} HVGs"
        + (f", {int(few.sum()):,} bins with < {MIN_HVG_DETECTED} of them detected left out" if few.any() else "")
        + f"  {time.time() - t0:.0f}s")

    pcs = _pca(Xn[:, cols].tocsr(), N_PCS)
    log(f"[umap] {key}: PCA {pcs.shape[1]}  {time.time() - t0:.0f}s")
    emb, clusterings, labels = embed_and_cluster(
        pcs, lambda lab, k: _markers(Xn, lab, k, er.genes), libs, key, log, t0)
    return write_umap_files(ud, {"nBins": n_rows, "nEmbedded": n, "minUmi": MIN_UMI,
                                 "minHvgDetected": MIN_HVG_DETECTED}, params,
                            emb, clusterings, labels, n_rows, row_of[bins], libs, key, log, t0)


# ---------------------------------------------------------------------------------------
# Protein layers (multiplexed immunofluorescence, protgrid.py): a dense matrix of marker
# means per grid square instead of counts. Same UMAP / Leiden / renumbering / files.
# ---------------------------------------------------------------------------------------
N_PCS_PROTEIN = 15
UMAP_SAMPLE = 200_000        # protein layers with more rows: fit on a seeded sample of this many, place the rest
SAMPLE_K = 15                # placed rows take the majority cluster of this many nearest sampled rows


def protein_params(channels: list[str], extra: dict | None = None, fit: str = "all rows") -> dict:
    return {"nPcs": f"min({N_PCS_PROTEIN}, channels - 1)", "nNeighbors": N_NEIGHBORS, "minDist": MIN_DIST,
            "init": INIT, "seed": SEED,
            "normalize": "per channel asinh(max(v - bg, 0) / hi / 0.2) / asinh(5); bg = off-tissue median, "
                         "hi = p99.5 of v - bg over the tissue squares (protgrid.normalise_bg_asinh)",
            "fit": fit,
            "scale": f"z-score per channel, clip {CLIP:g}",
            "cluster": f"leiden RBConfiguration on umap fuzzy graph, {LEIDEN_ITERS} iterations",
            "markers": (f"difference of the mean z-scored x' (cluster - rest, population sd), top {N_MARKERS} "
                        "with lfc > 0; pct = share of the cluster's squares whose x' is above that channel's "
                        "median over all squares"),
            "channels": list(channels), **(extra or {})}


def protein_markers(xp: np.ndarray, lab: np.ndarray, k: int, names: list[str], cols) -> list[list[dict]]:
    """Top N_MARKERS channels per cluster by mean z difference (cluster - rest); see protein_params.

    `xp` = rows x channels of x' in [0, 1]; only the columns in `cols` are ranked.
    Field names as for RNA: lfc = the mean z difference, pct = share above the median.
    """
    cols = np.asarray(cols, dtype=np.int64)
    X = xp[:, cols].astype(np.float64)
    sd = X.std(0)
    sd[sd == 0] = 1.0
    Z = (X - X.mean(0)) / sd
    above = (X > np.median(X, axis=0)).astype(np.float64)
    n = len(lab)
    nk = np.bincount(lab, minlength=k).astype(np.float64)
    M = sparse.csr_matrix((np.ones(n, np.float64), (lab, np.arange(n))), shape=(k, n))
    zs = np.asarray(M @ Z)                                   # k x channels sums of z
    ab = np.asarray(M @ above)
    tot = Z.sum(0)
    out = []
    for i in range(k):
        rest = max(n - nk[i], 1.0)
        diff = zs[i] / nk[i] - (tot - zs[i]) / rest
        pct = ab[i] / nk[i]
        top = np.argsort(-diff, kind="stable")[:N_MARKERS]
        out.append([{"gene": names[cols[j]], "lfc": round(float(diff[j]), 2), "pct": round(float(pct[j]), 2)}
                    for j in top if diff[j] > 0])
    return out


def write_feature_umap(d, xp: np.ndarray, names: list[str], embed_cols, marker_cols, key: str,
                       log=print, force: bool = False, extra_params: dict | None = None,
                       sample: int | None = None, features: np.ndarray | None = None,
                       keep: np.ndarray | None = None) -> dict | None:
    """expr/<key>/umap/ for a dense layer whose every positions.bin row is embedded.

    `xp` = rows x channels, already normalised (protgrid.normalise_bg_asinh: 0 at background,
    1 at p99.5). The embedding uses the columns `embed_cols` (z-scored, PCA to min(15, n - 1));
    the cluster markers rank `marker_cols`. `features` (rows x embedded channels): embed these
    instead of xp[:, embed_cols] (protgrid.clr_features); `keep`: the rows to embed (default:
    rows above background in some embedded channel). `sample`: above this many rows the UMAP is
    fitted on a seeded sample and the rest placed (None = UMAP_SAMPLE, 0 = always all rows).
    Returns the manifest entry, or None when skipped.
    """
    ud = d / "umap"
    embed_cols = [int(c) for c in embed_cols]
    # a square at or below background in every embedded channel has no profile: hundreds of
    # them are one identical point, which umap cannot place and Leiden splits into clusters of
    # nothing. They stay on the slide, unclustered (like RNA bins under MIN_UMI).
    if keep is None:
        keep = np.flatnonzero((xp[:, embed_cols] > 0).any(1)) if embed_cols else np.arange(0)
    keep = np.asarray(keep, dtype=np.int64)
    rows_fit = sample_rows(len(keep), sample)
    fit = ("all rows" if rows_fit is None else
           f"UMAP fitted on {len(rows_fit)} rows (seeded random sample), the rest placed with umap transform; "
           f"Leiden on the sample's graph, the rest take the most common cluster of their {SAMPLE_K} nearest "
           "sampled rows (PCA space); clusters renumbered by size over all rows")
    params = {**protein_params([names[c] for c in embed_cols], extra_params, fit)}
    params.setdefault("lowSignal", "at background in every embedded channel: not embedded")
    if not force:
        old = _reusable(ud, params, key, log)
        if old:
            return old
    n, ne = len(xp), len(keep)
    if ne < 50 or len(embed_cols) < 3:
        # two channels give one principal component, and umap's PCA start needs two (it stopped the export)
        log(f"[umap] {key}: {ne} rows x {len(embed_cols)} channels - too few to embed (3 channels or more), skipped")
        return None
    libs = _clustering_libs(key, log)
    if libs is None:
        return None
    t0 = time.time()
    n_pcs = min(N_PCS_PROTEIN, len(embed_cols) - 1)
    xk = xp[keep]
    fk = features[keep] if features is not None else xk[:, embed_cols]
    pcs = _pca(sparse.csr_matrix(np.ascontiguousarray(fk, dtype=np.float64)), n_pcs)
    log(f"[umap] {key}: {ne:,} of {n:,} rows ({n - ne:,} low signal left out), "
        f"{len(embed_cols)} channels, PCA {pcs.shape[1]}  {time.time() - t0:.0f}s")
    emb, clusterings, labels = embed_and_cluster(
        pcs, lambda lab, k: protein_markers(xk, lab, k, names, marker_cols), libs, key, log, t0, sample=rows_fit)
    head = {"nBins": n, "nEmbedded": ne, "modality": "protein"}
    if rows_fit is not None:
        head["fitOn"] = int(len(rows_fit))
    return write_umap_files(ud, head, params,
                            emb, clusterings, labels, n, keep, libs, key, log, t0)


def _entry(index: dict) -> dict:
    """The manifest's view of an index.json."""
    return {"nEmbedded": index["nEmbedded"],
            "clusterings": [{**{k: c[k] for k in ("id", "resolution", "name", "source") if c.get(k) is not None},
                             "n": c["n"]} for c in index["clusterings"]],
            "default": index["default"]}
