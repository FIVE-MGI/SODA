"""Marker expression on a grid: square means of a multiplexed fluorescence image as an expr layer.

Akoya / PhenoCycler (and any multi-channel fluorescence TIFF) has no count matrix, so
there is nothing for the gene search, the whole-section UMAP or the clusters to read.
This lays a square grid over the image, takes each marker channel's mean intensity per
square, and writes it as an ordinary expression layer of the image dataset, so all of
that machinery works unchanged (a lab's own cell segmentation becomes a cells layer beside it, by the same
rules: protcells.py; SODA segments nothing):

    expr/grid<NNN>um/positions.bin     float32 x,y  square centres in full-res image px (tissue squares only)
    expr/grid<NNN>um/genes.json        [{n, t}]     one per marker channel written; t = mean over the squares
    expr/grid<NNN>um/genepack/         the per-channel float32 x,y,v arrays packed (genepack.py; s/o/b in
                                       genes.json):
    expr/grid<NNN>um/gene/<NAME>.bin   float32 x,y,v  EVERY square, in positions.bin order (v may be 0;
                                                    the viewer treats v <= 0 as "nothing here", as for RNA)
    expr/grid<NNN>um/umap/...          as sectionumap.py (index.json, coords.bin, clusters_r*.bin)
                                       + celltypes.json / celltype.bin (celltypes.py, protein rule)
    expr/grid<NNN>um/mat/...           the embedded channels square by square, in spatial chunks, for the
                                       circle -> UMAP tool (write_mat_chunks). No hvg/: that is for counts.
    expr/grid<NNN>um/qc.json           per channel: background, p99.5 above it, signal-to-background,
                                       saturated %, status (ok / weak / saturated / flat / nuclear); the
                                       channels dropped and why; flagged squares; the seam check
    expr/grid<NNN>um/qc_flags.bin      uint8 per positions.bin row: bit 1 = everything bright (fold, bubble,
                                       edge), bit 2 = many saturated pixels

v is in the 8-bit units of the tiles: the raw mean for a uint8 pyramid, the raw mean
mapped through the channel's `bits` window for uint16 (not clipped at 255).

Squares: `grid_um` wide (default 16, as the Visium HD 16 µm bins), starting at the image's
top-left corner; a last partial row / column is kept when it is at least half a square. The
means come from FULL-RESOLUTION pixels (level 0, streamed one chunk row at a time, each chunk read
once per channel): the coarser pyramid levels are max-pooled for display and overstate markers
(square_means). The same pass counts pixels at the detector's limit per square and channel.

Tissue: the DAPI (or Hoechst) square means (without a nuclear stain, the mean over the written
channels) are first averaged over a window of `tissue_smooth_um` (default 100 µm = 7 x 7 squares
of 16 µm; moving average, edges by nearest-square replication); a square is kept when that value
is above an Otsu threshold over all squares, taken on log1p (see tissue_mask for why not linear).
Enclosed holes of <= HOLE_MAX squares are filled. tissue_smooth_um = 0 thresholds the raw square
means (the rule before 2026-11: it keeps only nucleus-dense squares, a lacy mask that misses axon /
collagen stretches).

Channels, in order (each reason goes to the manifest's "dropped" and qc.json):
  1. names matching the viewer's hidden-channel default (open / blank / empty / unused / none):
     not written;
  2. duplicates - the same name after celltypes.norm_name (case, "01_" prefixes, fluorophore /
     FFPE tags): the copy with the stronger signal-to-background is kept, the others not written;
  3. flat over the tissue (p99.5 - background < FLAT_RANGE units): written, not embedded;
  4. nuclear stains: written (they make the tissue mask), not embedded.
Weak (signal-to-background < WEAK_SNR) and saturated (> SAT_SHARE of in-tissue pixels at the
limit) channels stay in the embedding and are flagged in qc.json.

Normalization (per channel, per slide - pieces on one slide were stained together):
  bg = median over the off-tissue squares (background); hi = p99.5 of (v - bg) over the tissue
  squares, leaving out the "everything bright" ones (bright_squares: folds would set it);
  x = asinh(max(v - bg, 0) / hi / 0.2) / asinh(5)   (0 at background, 1 at p99.5, no upper clip).
Embedding (sectionumap.write_feature_umap): squares with an embedded channel at x >= LOW_SIGNAL_X;
CLR of u = max(v - bg, 0) / hi (clr_features) -> z-score (clip 10) -> PCA min(15, channels - 1) ->
UMAP (15 neighbours, min_dist 0.3, seed 42; above sectionumap.UMAP_SAMPLE squares fitted on a seeded
sample, the rest placed with transform) -> Leiden 0.3 / 0.6 / 1.0, cluster 0 = largest. Markers:
difference of the mean z of x (cluster - rest); pct = share of the cluster's squares whose x is
above the channel's median. Exports before 2026-11 used x' = clip((v - p1) / (p99.5 - p1), 0, 1)
on means from a max-pooled level (normalise, kept for reference).
"""
from __future__ import annotations

import hashlib
import json
import math
import re
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import quote

import numpy as np

import deid
import genepack
import seams

DEFAULT_GRID_UM = 16.0
NUCLEAR = re.compile(r"^\s*(dapi|hoechst)\b", re.I)                    # = viewer.js NUCLEAR_CHANNEL
HIDDEN = re.compile(r"^\s*(open|blank|empty|unused|none)([\s_-]|$)", re.I)   # = viewer.js HIDE_CHANNEL default
FLAT_RANGE = 2.0              # p99.5 - background over the tissue squares, 8-bit units: below this a channel is flat
P_LO, P_HI = 1.0, 99.5
HOLE_MAX = 4                  # tissue mask: enclosed background holes up to this many squares are filled
DEFAULT_TISSUE_SMOOTH_UM = 100.0   # tissue mask: moving-average window over the square means (0 = none)


def _encode(name: str) -> str:
    return quote(name, safe="-_.!~*'()")            # = export._encode = JavaScript encodeURIComponent


def layer_key(grid_um: float) -> str:
    g = float(grid_um)
    return f"grid{int(g):03d}um" if g.is_integer() else "grid" + f"{g:g}".replace(".", "p") + "um"


def is_photo(channels: list[dict]) -> bool:
    """An RGB photo (H&E / brightfield) has no markers to average."""
    names = [str(c.get("name", "")).strip().lower() for c in channels]
    return bool(names) and all(n.split(" ")[-1] in ("red", "green", "blue") for n in names)


def applies_to(m: dict) -> bool:
    """Image datasets only: not Visium (it has real counts), not RNA-seq, not a colour photo,
    not a density background standing in for missing images (cosmx.py / merscope.py). A lab's own cells of
    the image (protcells.py, modality "protein") do not count as expression: the grid is still made."""
    expr = [L for L in m.get("expr") or [] if L.get("modality") != "protein"]
    return (m.get("kind") not in ("visium", "rnaseq") and not expr and "image" in m
            and not m["image"].get("background") and not is_photo(m["image"].get("channels", [])))


def edges(n_px: int, step: float) -> np.ndarray:
    """Pixel edges of the squares along one axis (a last partial square is kept if >= half a square)."""
    n = int(math.floor(n_px / step + 0.5))
    return np.minimum(np.round(np.arange(n + 1) * step).astype(np.int64), n_px)


def channel_sums(a, c: int, ye: np.ndarray, xe: np.ndarray, top: int | None = None,
                 profiles: bool = False):
    """Sum of channel `c` of zarr array `a` (C, H, W) over every square, one chunk row at a time.

    With `top`, also the count of pixels >= top per square (saturation), returned as
    (sums, counts). With `profiles`, also the mean of every image column and row, and the same
    over the top / bottom half of the rows and the left / right half of the columns (seams.find_grid):
    (sums, counts, (colMean, rowMean, colTop, colBottom, rowLeft, rowRight))."""
    ny, nx = len(ye) - 1, len(xe) - 1
    S = np.zeros((ny, nx), np.float64)
    T = np.zeros((ny, nx), np.int64) if top is not None else None
    rows = a.chunks[1]
    h_used, w_used = int(ye[-1]), int(xe[-1])
    if profiles:
        colsum, coltop = np.zeros(w_used, np.float64), np.zeros(w_used, np.float64)
        rowmean, rowleft, rowright = (np.zeros(h_used, np.float64) for _ in range(3))
        hmid, wmid = h_used // 2, w_used // 2
    buf, b0, j = None, 0, 0
    for y0 in range(0, h_used, rows):
        y1 = min(y0 + rows, h_used)
        band = np.asarray(a[c, y0:y1, :w_used])
        if profiles:
            colsum += band.sum(0, dtype=np.float64)
            if y0 < hmid:
                coltop += band[:hmid - y0].sum(0, dtype=np.float64)
            rowmean[y0:y1] = band.mean(1, dtype=np.float64)
            rowleft[y0:y1] = band[:, :wmid].mean(1, dtype=np.float64)
            rowright[y0:y1] = band[:, wmid:].mean(1, dtype=np.float64)
        buf = band if buf is None or len(buf) == 0 else np.concatenate([buf, band])
        j1 = j
        while j1 < ny and ye[j1 + 1] <= y1:
            j1 += 1
        if j1 > j:
            blk = buf[ye[j] - b0: ye[j1] - b0]
            rs = np.add.reduceat(blk, ye[j:j1] - ye[j], axis=0, dtype=np.uint64)
            S[j:j1] = np.add.reduceat(rs, xe[:-1], axis=1)
            if T is not None:
                hit = np.add.reduceat((blk >= top).view(np.uint8), ye[j:j1] - ye[j], axis=0, dtype=np.int64)
                T[j:j1] = np.add.reduceat(hit, xe[:-1], axis=1)
            buf = buf[ye[j1] - b0:]
            b0, j = int(ye[j1]), j1
    if T is None:
        return S
    if profiles:
        return S, T, (colsum / max(h_used, 1), rowmean, coltop / max(hmid, 1),
                      (colsum - coltop) / max(h_used - hmid, 1), rowleft, rowright)
    return S, T


def otsu(v: np.ndarray, bins: int = 256) -> float:
    """Otsu threshold of the values (histogram of `bins` bins between min and max)."""
    v = np.asarray(v, np.float64)
    lo, hi = float(v.min()), float(v.max())
    if hi <= lo:
        return lo
    h, e = np.histogram(v, bins=bins, range=(lo, hi))
    p = h.astype(np.float64) / h.sum()
    mid = (e[:-1] + e[1:]) / 2
    w0 = np.cumsum(p)
    m0 = np.cumsum(p * mid)
    mt = m0[-1]
    w1 = 1 - w0
    with np.errstate(divide="ignore", invalid="ignore"):
        between = (mt * w0 - m0) ** 2 / (w0 * w1)
    between[~np.isfinite(between)] = -1
    k = int(np.argmax(between))
    return float(e[k + 1])


def tissue_mask(src: np.ndarray, ny: int, nx: int) -> tuple[np.ndarray, float, int, int]:
    """Squares on tissue: Otsu on log1p of the nuclear-stain means, then enclosed holes of
    <= HOLE_MAX squares filled. Returns (keep, threshold in intensity units, holes, squares filled).

    On the linear scale the near-black background (most of a slide) and the long tail of
    nucleus-dense squares make Otsu split the TISSUE in two: on a 30-marker nerve section it
    kept 16% of the squares, leaving every nucleus-poor stretch (axons, collagen) out. On
    log1p it keeps the tissue outline (43%). The hole fill only closes single squares or
    tiny gaps fully surrounded by tissue, never a tear or a lumen.
    """
    from scipy import ndimage

    t = otsu(np.log1p(np.maximum(src, 0)))
    keep = np.log1p(np.maximum(src, 0)) > t
    grid = keep.reshape(ny, nx)
    lab, k = ndimage.label(~grid)                      # background components (4-connected)
    sizes = np.bincount(lab.ravel(), minlength=k + 1)
    edge = np.unique(np.concatenate([lab[0], lab[-1], lab[:, 0], lab[:, -1]]))
    small = (sizes <= HOLE_MAX)
    small[0] = False
    small[edge] = False                                # touching the border: not enclosed
    fill = small[lab]
    grid = grid | fill
    return grid.ravel(), float(np.expm1(t)), int(small.sum()), int(fill.sum())


def smooth_squares(n_um: float, grid_um: float) -> int:
    """Moving-average window in squares for `n_um` µm: the nearest odd count (so it centres on the
    square), >= 1; 0 or 1 = no smoothing. 100 µm over 20 µm squares = 5."""
    if not n_um or n_um <= 0:
        return 1
    k = int(round(float(n_um) / float(grid_um)))
    if k % 2 == 0:
        k = k + 1 if float(n_um) / float(grid_um) >= k else k - 1
    return max(k, 1)


def tissue_source(src: np.ndarray, ny: int, nx: int, size: int) -> np.ndarray:
    """The values the tissue mask thresholds: the square means after a `size` x `size`-square moving
    average (scipy uniform_filter, mode "nearest": edge squares replicated), or the means as they are
    for size <= 1. Float64 in, float64 out (= make_tiles.py's area mask)."""
    from scipy import ndimage

    src = np.asarray(src, np.float64)
    if size <= 1:
        return src
    return ndimage.uniform_filter(src.reshape(ny, nx), size=size, mode="nearest").ravel()


def tissue_rule(size: int, grid_um: float) -> str:
    if size <= 1:
        return "Otsu on log1p of the means"
    return (f"Otsu on log1p of the means after a {size}x{size}-square ({size * grid_um:g} µm) "
            "moving average (edges: nearest square)")


def normalise(V: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """The rule before 2026-11 (old exports; the circle tool still reads it from them):
    x' = clip((v - p1) / (p99.5 - p1), 0, 1) per column; also returns p1 and p99.5."""
    V = np.asarray(V, np.float64)
    lo = np.percentile(V, P_LO, axis=0)
    hi = np.percentile(V, P_HI, axis=0)
    rng = hi - lo
    xp = np.where(rng > 0, (V - lo) / np.where(rng > 0, rng, 1.0), 0.0)
    return np.clip(xp, 0.0, 1.0), lo, hi


# ---- normalization (2026-11): background off, / p99.5, asinh -----------------------------
NORMALIZE = "bg-asinh-v1"     # mat/index.json "normalize": lasso_worker.js xpOf switches on it
COFACTOR = 0.2                # asinh cofactor, in units of the channel's p99.5 above background
# What the UMAP and Leiden see (not what is shown or how markers are ranked: those stay on x).
# Every marker of a square rises and falls with its overall brightness (cell density, section
# thickness, the tissue edge): on x the first PC was that brightness (48% of the variance, r =
# -0.999 with the square's total) and the clusters were brightness bands. The centred log-ratio
# of u = max(v - bg, 0) / hi, log(u + CLR_EPS) minus the square's mean, keeps which markers are
# up relative to each other and drops the overall level. On s48 that split the perineurium,
# endothelium, myeloid, T-cell, B-cell and Ki67 squares into their own clusters, where before
# they shared bright clusters (audit 2026-11-10; eps 0.02 and 0.2 did worse).
CLR_EPS = 0.1
# The log-ratio is taken against the square's mean over the embedded channels: with only a few
# channels one bright marker moves that mean as much as the brightness does (the 6-channel
# synthetic demo lost ARI 0.43 -> 0.34 with it). Smaller panels keep x as the UMAP's input.
CLR_MIN_CHANNELS = 10
# A square with no embedded marker at x >= LOW_SIGNAL_X (about 10% of its p99.5 above
# background) is a rim, gap or empty square: its log-ratios are noise, and on s48 such squares
# (26%) formed clusters with no marker. They stay on the slide but are not embedded or clustered.
LOW_SIGNAL_X = 0.2
BG_MIN_OFF = 1000             # off-tissue squares needed for the background; fewer -> 1st pct in tissue
WEAK_SNR = 15.0               # (p99.5 - bg) / bg below this: "weak" (kept, flagged)
SAT_SHARE = 0.01              # share of in-tissue pixels at the detector's limit above this: "saturated"
SAT_SQUARE = 0.5              # a square with more than this share of its pixels at the limit is flagged
                              # (0.10 flagged 25% of s48's squares: bright stains clip in every dense patch)
BRIGHT_PCT, BRIGHT_SHARE = 95.0, 0.8   # "everything bright": > 80% of the marker channels above their p95
FLAG_BRIGHT, FLAG_SAT = 1, 2  # qc_flags.bin bits


def background(V: np.ndarray, keep: np.ndarray) -> tuple[np.ndarray, str]:
    """Per column: the median over OFF-tissue squares, or, with fewer than BG_MIN_OFF of them
    (tissue fills the slide), the 1st percentile over the tissue squares."""
    V = np.asarray(V, np.float64)
    off = ~np.asarray(keep, bool)
    if off.sum() >= BG_MIN_OFF:
        return np.median(V[off], axis=0), "median of the off-tissue squares"
    return np.percentile(V[keep], P_LO, axis=0), "1st percentile of the tissue squares (too little background)"


def clr_features(Vt: np.ndarray, bg: np.ndarray, hi: np.ndarray, eps: float = CLR_EPS) -> np.ndarray:
    """The UMAP's input for protein squares: centred log-ratio of u = max(v - bg, 0) / hi.
    lasso_worker.js runProtein computes the same (keep the two identical)."""
    h = np.where(hi > 0, hi, 1.0)
    u = np.where(hi > 0, np.maximum(Vt.astype(np.float64) - bg, 0) / h, 0.0)
    L = np.log(u + eps)
    return (L - L.mean(1, keepdims=True)).astype(np.float32)


def normalise_bg_asinh(Vt: np.ndarray, bg: np.ndarray, hi: np.ndarray, cof: float = COFACTOR) -> np.ndarray:
    """x = asinh(u / cof) / asinh(1 / cof), u = max(v - bg, 0) / hi: 0 at background, 1 at the
    channel's p99.5, no upper clip (a fold at 3x p99.5 reads 1.47). hi <= 0 -> 0."""
    Vt = np.asarray(Vt, np.float64)
    ok = np.asarray(hi, np.float64) > 0
    u = np.maximum(Vt - bg, 0.0) / np.where(ok, hi, 1.0)
    return np.where(ok, np.arcsinh(u / cof) / np.arcsinh(1.0 / cof), 0.0)


def saturation_top(dtype, channels: list[dict]) -> int:
    """The pixel value that means "at the detector's limit", one per slide: 255 for uint8; for
    uint16 4,095 when no channel goes above it (a 12-bit camera), else 65,535."""
    if str(dtype) == "uint8":
        return 255
    dm = [float(ch["dataMax"]) for ch in channels if ch.get("dataMax") is not None]
    return 4095 if dm and max(dm) <= 4095 else 65535


def bright_squares(Vt: np.ndarray) -> np.ndarray:
    """"Everything bright" (folds, bubbles, bright edges): more than BRIGHT_SHARE of the columns
    above their own BRIGHT_PCT-th percentile. Ranks per column, so it does not depend on the
    normalization, and the normalization can leave these squares out of its p99.5."""
    Vt = np.asarray(Vt, np.float64)
    if Vt.shape[1] == 0:
        return np.zeros(len(Vt), bool)
    return (Vt > np.percentile(Vt, BRIGHT_PCT, axis=0)).mean(1) > BRIGHT_SHARE


def to_8bit(mean: np.ndarray, channel: dict, dtype) -> np.ndarray:
    """Raw means -> the 8-bit units of the tiles (the channel's `bits` window; not clipped at 255)."""
    b0, b1 = (float(x) for x in channel["bits"])
    if str(dtype) == "uint8" and b0 == 0 and b1 == 255:
        return mean
    return np.maximum((mean - b0) * (255.0 / max(b1 - b0, 1e-9)), 0.0)


def square_means(r, grid_um: float, px_um: float, log=print, workers: int = 4) -> dict:
    """Mean per square and channel (8-bit units) over the whole grid, from FULL-RESOLUTION pixels.

    The coarser pyramid levels of a fluorescence image are 2x2 MAX-pooled (pyramid.py, so the
    display keeps small bright structures), so a mean taken from them overstates every marker:
    on s48, 16 µm, DAPI by ~20% and sparse markers (Ki67) by ~2x. Level 0 is read once per
    channel, a chunk row at a time. Also counts, per square and channel, the pixels at the
    detector's limit (saturation_top), and every channel's column / row means (seams.find_grid).
    """
    a = r._arrays[0]
    step = grid_um / px_um
    if step < 2:
        raise ValueError(f"{grid_um} µm squares are {step:.1f} px at full resolution: too small")
    ye, xe = edges(r.H, step), edges(r.W, step)
    ny, nx = len(ye) - 1, len(xe) - 1
    tops = [saturation_top(a.dtype, r.channels)] * len(r.channels)
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=workers) as ex:
        res = list(ex.map(lambda c: channel_sums(a, c, ye, xe, top=tops[c], profiles=True), range(len(r.channels))))
    npx = np.outer(np.diff(ye), np.diff(xe)).astype(np.float64)
    V = np.empty((ny * nx, len(r.channels)), np.float32)
    V64 = np.empty((ny * nx, len(r.channels)), np.float64)   # for the tissue mask (not rounded to float32)
    SAT = np.empty((ny * nx, len(r.channels)), np.int32)     # pixels at the limit per square
    prof = []
    for c, (S, T, pr) in enumerate(res):
        prof.append(pr)
        mean = to_8bit((S / npx).ravel(), r.channels[c], a.dtype)
        V[:, c] = mean
        V64[:, c] = mean
        SAT[:, c] = T.ravel()
    log(f"[grid] full resolution ({px_um:.3f} µm/px, {r.W}x{r.H}), {step:.2f} px squares, {nx}x{ny} = {nx * ny:,} "
        f"squares x {len(r.channels)} channels  {time.time() - t0:.0f}s")
    return {"V": V, "V64": V64, "nx": nx, "ny": ny, "level": 0, "pxUm": float(px_um), "sat": SAT,
            "npx": npx.ravel(), "tops": tops, "profiles": prof, "ye": ye, "xe": xe}


MAT_SCALE = 100              # mat/ chunks store round(v * MAT_SCALE) as uint16 (0.01 of an 8-bit unit)
MAT_CHUNK_UM = 500.0         # = hvgchunks.DEFAULT_CHUNK_UM


def write_mat_chunks(d: Path, pos: np.ndarray, V: np.ndarray, names: list[str], bg, hi, mpp: float,
                     chunk_um: float = MAT_CHUNK_UM, log=print, cofactor: float = COFACTOR) -> dict:
    """expr/<key>/mat/: the embedded channels square by square, cut into squares of tissue.

    The circle -> UMAP tool needs every channel for the few squares under a shape - the
    gene files hold one channel for every square - so, like hvg/ for counts:

        mat/index.json      {chunkPx, chunkUm, scale, channels, normalize, bg, hi, cofactor, median,
                             nSquares, mean, sd, rawMean, chunks: [{tx, ty, n, bytes}]}
        mat/<tx>_<ty>.bin   little-endian: u32 nSquares, u32 row[nSquares] (index into positions.bin),
                            u16 value[nSquares * nChannels] (square-major, round(v * scale), clipped at 65535)

    `V` = squares x channels (positions.bin order), raw means; `bg` / `hi` each channel's whole-slide
    background and p99.5 above it, so a circle is scaled as the section UMAP was (normalise_bg_asinh,
    "normalize": NORMALIZE). `median` = each channel's median x over the section (the markers' pct).
    `mean` / `sd` (of x) and `rawMean` (of v) over every square of the slide, from the stored
    (rounded) values, so a shape can be compared with the rest of the slide and its clusters'
    markers and suggested names read against the slide, not only against the other squares drawn.
    Exports before 2026-11 carry p1 / p99_5 instead and no "normalize" (x' = clip((v - p1) / (p99.5 - p1), 0, 1)).
    """
    V = np.asarray(V, np.float64)
    bg, hi = np.asarray(bg, np.float64), np.asarray(hi, np.float64)
    xp = normalise_bg_asinh(V, bg, hi, cofactor)
    med = np.median(xp, axis=0)
    chunk_px = int(round(chunk_um / mpp))
    xy = np.asarray(pos, np.float64)
    tx = np.floor(xy[:, 0] / chunk_px).astype(np.int64)
    ty = np.floor(xy[:, 1] / chunk_px).astype(np.int64)
    rows = np.arange(len(xy))
    q = np.minimum(np.rint(np.maximum(V, 0.0) * MAT_SCALE), 65535).astype("<u2")
    vq = q.astype(np.float64) / MAT_SCALE                # what the browser reads back
    xq = normalise_bg_asinh(vq, bg, hi, cofactor)
    md = Path(d) / "mat"
    md.mkdir(parents=True, exist_ok=True)
    order = np.lexsort((rows, ty, tx))
    keys = np.stack([tx[order], ty[order]], 1)
    split = np.flatnonzero(np.any(np.diff(keys, axis=0) != 0, axis=1)) + 1
    chunks, total = [], 0
    for grp in np.split(order, split):
        buf = b"".join([np.uint32(len(grp)).astype("<u4").tobytes(), grp.astype("<u4").tobytes(),
                        np.ascontiguousarray(q[grp]).tobytes()])
        cx, cy = int(tx[grp[0]]), int(ty[grp[0]])
        _atomic(md / f"{cx}_{cy}.bin", buf)
        chunks.append({"tx": cx, "ty": cy, "n": int(len(grp)), "bytes": len(buf)})
        total += len(buf)
    r6 = lambda a: [round(float(x), 6) for x in a]
    index = {"chunkPx": chunk_px, "chunkUm": chunk_um, "scale": MAT_SCALE, "nChannels": len(names),
             "channels": list(names), "normalize": NORMALIZE, "bg": r6(bg), "hi": r6(hi), "cofactor": cofactor,
             "median": r6(med), "nSquares": int(len(xy)),
             "mean": r6(xq.mean(axis=0)), "sd": r6(xq.std(axis=0)), "rawMean": r6(vq.mean(axis=0)),
             "normalizeRule": "x = asinh(max(v - bg, 0) / hi / cofactor) / asinh(1 / cofactor), "
                              "bg and hi over the whole slide",
             "chunks": chunks}
    deid.assert_clean(index, "mat/index.json")
    _atomic(md / "index.json", json.dumps(index, separators=(",", ":")).encode("utf-8"))
    log(f"[grid] mat/: {len(names)} channels, {len(chunks)} chunks of {chunk_um:g} um, {total / 1e6:.1f} MB")
    return {"nChannels": len(names), "chunkPx": chunk_px, "chunks": len(chunks), "bytes": total}


def _atomic(p: Path, buf: bytes):
    tmp = p.with_name(p.name + ".part")
    tmp.write_bytes(buf)
    tmp.replace(p)


def write_grid_layer(r, m: dict, out: Path, grid_um: float = DEFAULT_GRID_UM, markers: Path | None = None,
                     umap: bool = True, celltypes: bool = True, log=print, force: bool = False,
                     tissue_smooth_um: float = DEFAULT_TISSUE_SMOOTH_UM, umap_sample: int | None = None,
                     stats: dict | None = None) -> dict | None:
    """Write expr/<key>/ for an image dataset (PyramidReader `r`, manifest `m`) into export folder `out`.

    `stats`: a dict to fill with {"background": {channel name: background}} for every channel (protcells.py
    takes the slide's background from it); nothing written depends on it.

    `tissue_smooth_um`: moving-average window (µm) over the square means before the tissue
    threshold; 0 = threshold the raw means (the old, lacy mask). `umap_sample`: fit the section
    UMAP on this many squares and place the rest (None = sectionumap.UMAP_SAMPLE, 0 = fit on all).
    Tile seams are found and reported in qc.json (seams.py); the square means are not corrected.

    Returns the manifest's expr entry, or None when the dataset has no pixel size.
    """
    import celltypes as celltypes_mod
    import sectionumap

    px_um = m.get("pixelSizeUm")
    if not px_um:
        log("[grid] skipped: the dataset has no pixel size (prepare with --pixel-size)")
        return None
    grid_um = float(grid_um)
    key = layer_key(grid_um)
    t_all = time.time()
    names = [str(ch["name"]) for ch in r.channels]

    # ---- 1. channels: open / blank by name ----------------------------------------------
    written = [j for j, n in enumerate(names) if not HIDDEN.match(n)]
    dropped = [{"name": names[j], "reason": "open / blank channel: not written"}
               for j in range(len(names)) if j not in written]
    if not written:
        log("[grid] skipped: every channel is an open / blank channel")
        return None
    nuc = [j for j in written if NUCLEAR.match(names[j])]

    g = square_means(r, grid_um, float(px_um), log=log)
    V, V64, nx, ny = g["V"], g["V64"], g["nx"], g["ny"]

    # ---- 2. tissue mask --------------------------------------------------------------------
    by = names[nuc[0]] if nuc else "mean of the channels"
    src = V64[:, nuc[0]] if nuc else V64[:, written].mean(1)
    size = smooth_squares(tissue_smooth_um, grid_um)
    rule = tissue_rule(size, grid_um)
    keep, thr, n_holes, n_filled = tissue_mask(tissue_source(src, ny, nx, size), ny, nx)
    if keep.sum() < 2:
        log("[grid] skipped: the tissue mask keeps no squares")
        return None
    rows = np.flatnonzero(keep)
    n = len(rows)
    sm_txt = f"{size}x{size}-square ({size * grid_um:g} µm) average of the {by} mean" if size > 1 else f"{by} mean"
    log(f"[grid] tissue mask: {sm_txt} > {thr:.2f} (Otsu on log1p over {len(src):,} squares), "
        f"+ {n_filled} squares in {n_holes} enclosed holes of <= {HOLE_MAX} squares: keeps {n:,} ({n / len(src):.0%})")

    # ---- 2b. tile seams: the camera grid (found and reported in qc.json; nothing is corrected) ---
    seam = seams.find_grid([g["profiles"][j] for j in written], [names[j] for j in written], float(px_um))
    if seam["flagged"]:
        log("[grid] possible tile seams (camera-tile edges) found: listed in the slide QC; the square means are not "
            "corrected" + (f" ({seam['gridProblem']})" if seam.get("gridProblem") else ""))
    V, V64 = g["V"], g["V64"]

    # ---- 3. per channel: background, p99.5 above it, signal, saturation --------------------
    # p99.5 leaves out the "everything bright" squares (bright_squares, on the non-nuclear written
    # channels): a fold covering more than 0.5% of the tissue would otherwise set every channel's
    # p99.5 and squash the real signal (the demo's 150 µm fold made fascicles read as perineurium)
    bg_all, bg_rule = background(V64, keep)
    if stats is not None:
        stats["background"] = {nm: float(b) for nm, b in zip(names, bg_all)}
        stats["backgroundRule"] = bg_rule
    Vt_all = V64[rows]
    bright = bright_squares(Vt_all[:, [j for j in written if not NUCLEAR.match(names[j])]])
    ref = ~bright if (~bright).sum() >= max(50, n // 2) else np.ones(n, bool)
    hi_all = np.percentile(np.maximum(Vt_all[ref] - bg_all, 0.0), P_HI, axis=0)
    snr_all = hi_all / np.maximum(bg_all, 0.5)
    npx_t = float(g["npx"][rows].sum())
    sat_all = g["sat"][rows].sum(0) / max(npx_t, 1.0)

    # ---- 4. duplicates (same name after the name cleanup): keep the stronger ----------------
    groups: dict[str, list[int]] = {}
    for j in written:
        groups.setdefault(celltypes_mod.norm_name(names[j]), []).append(j)
    for js in groups.values():
        if len(js) < 2:
            continue
        best = max(js, key=lambda j: (snr_all[j], -j))
        for j in js:
            if j != best:
                written.remove(j)
                dropped.append({"name": names[j], "reason": f"duplicate of channel {best + 1} ({names[best]}), "
                                                            f"which has the stronger signal: not written"})
    nuc = [j for j in written if NUCLEAR.match(names[j])]

    Vt = Vt_all[:, written]                            # tissue squares x written channels
    wnames = [names[j] for j in written]
    bg, hi, snr, sat = bg_all[written], hi_all[written], snr_all[written], sat_all[written]
    flat = hi < FLAT_RANGE
    embed, marker_cols = [], []
    for k, j in enumerate(written):
        if flat[k]:
            dropped.append({"name": names[j], "reason": f"flat over the tissue (p99.5 - background = {hi[k]:.2f}): "
                                                         "written, not embedded"})
            continue
        marker_cols.append(k)
        if NUCLEAR.match(names[j]):
            dropped.append({"name": names[j], "reason": "nuclear stain: written (and the tissue mask), not embedded"})
        else:
            embed.append(k)
    x = normalise_bg_asinh(Vt, bg, hi)
    log(f"[grid] {len(written)} channels written, {len(embed)} embedded: " + ", ".join(wnames[k] for k in embed))
    for dr in dropped:
        log(f"[grid]   left out: {dr['name']} - {dr['reason']}")

    # ---- 5. QC: weak / saturated channels, everything-bright and saturated squares, seams ----
    flags = np.zeros(n, np.uint8)
    flags[bright] |= FLAG_BRIGHT
    sat_sq = (g["sat"][rows][:, written] / g["npx"][rows][:, None]) > SAT_SQUARE
    flags[sat_sq.any(1)] |= FLAG_SAT
    seams_qc = seam
    qc_ch = []
    for k, nm in enumerate(wnames):
        status = ("flat" if flat[k] else "nuclear" if k in marker_cols and k not in embed
                  else "saturated" if sat[k] > SAT_SHARE else "weak" if snr[k] < WEAK_SNR else "ok")
        qc_ch.append({"name": nm, "inUmap": k in embed, "status": status,
                      "background": round(float(bg[k]), 3), "p99_5AboveBg": round(float(hi[k]), 3),
                      "signalToBackground": round(float(snr[k]), 2), "saturatedPct": round(float(sat[k]) * 100, 3),
                      "weak": bool(snr[k] < WEAK_SNR), "saturated": bool(sat[k] > SAT_SHARE)})
    n_bright, n_sat = int((flags & FLAG_BRIGHT).astype(bool).sum()), int((flags & FLAG_SAT).astype(bool).sum())
    qc = {"channels": qc_ch,
          "dropped": [{"name": dr["name"], "reason": dr["reason"]} for dr in dropped],
          "squares": n,
          "everythingBright": {"n": n_bright, "share": round(n_bright / n, 5),
                               "rule": f"more than {BRIGHT_SHARE:.0%} of the marker channels (all written but the "
                                       f"nuclear stain) above their {BRIGHT_PCT:g}th percentile over the tissue "
                                       "(folds, bubbles, bright edges); left out of each channel's p99.5"},
          "saturatedSquares": {"n": n_sat, "share": round(n_sat / n, 5),
                               "rule": f"more than {SAT_SQUARE:.0%} of the square's pixels at the detector's limit "
                                       "in any written channel"},
          "seams": seams_qc,
          "rules": {"weak": f"(p99.5 - background) / background < {WEAK_SNR:g}: kept in the UMAP, flagged",
                    "saturated": f"more than {SAT_SHARE:.0%} of the in-tissue pixels at the detector's limit "
                                 f"(255 for 8-bit images; 4,095 or 65,535 for 16-bit): kept, flagged",
                    "background": bg_rule, "limit": sorted({int(t) for t in g["tops"]}),
                    "p99_5": "99.5th percentile of (v - background) over the tissue squares, leaving out the "
                             "everything-bright squares"},
          "note": "UMAPs, clusters and names are per slide: their colours and positions are not comparable "
                  "across slides.",
          "flagsFile": "qc_flags.bin", "flagBits": {"everythingBright": FLAG_BRIGHT, "saturated": FLAG_SAT}}
    weak = [c["name"] for c in qc_ch if c["weak"] and c["inUmap"]]
    satd = [f"{c['name']} {c['saturatedPct']:.1f}%" for c in qc_ch if c["saturated"]]
    log(f"[grid] QC: weak {', '.join(weak) or 'none'}; saturated {', '.join(satd) or 'none'}; "
        f"everything-bright squares {n_bright:,}; squares with saturated pixels {n_sat:,}; "
        f"seams {'FLAGGED' if seams_qc['flagged'] else 'not flagged'}"
        + (f", corrected in {len(seams_qc['correction']['corrected'])} channels" if seams_qc.get("correction") else ""))

    # ---- files ---------------------------------------------------------------------------
    step_full = grid_um / float(px_um)
    iy, ix = np.divmod(rows, nx)
    pos = np.empty((n, 2), np.float32)
    pos[:, 0] = (ix + 0.5) * step_full
    pos[:, 1] = (iy + 0.5) * step_full
    d = Path(out) / "expr" / key
    d.mkdir(parents=True, exist_ok=True)
    _atomic(d / "positions.bin", np.ascontiguousarray(pos, dtype="<f4").tobytes())
    genes = [{"n": nm, "t": round(float(Vt[:, k].mean()), 4)} for k, nm in enumerate(wnames)]
    deid.assert_clean(genes, "grid genes.json")
    xyv = np.empty((n, 3), np.float32)
    xyv[:, :2] = pos
    pw = genepack.PackWriter(d, unit=12, log=log)
    where = {}
    for k, nm in enumerate(wnames):
        xyv[:, 2] = Vt[:, k]
        buf = np.ascontiguousarray(xyv, dtype="<f4")
        if pw is not None and n and nm not in where:
            where[nm] = pw.add(nm, buf)
    pack = pw.finish() if pw is not None else None
    if pack:
        genepack.note_loose(d, log=log)
    _atomic(d / "genes.json", json.dumps(genepack.with_index(genes, where), separators=(",", ":")).encode("utf-8"))
    log(f"[grid] {key}: {n:,} squares, {len(wnames)} channels -> {d}"
        + (f" ({len(pack['shards'])} pack shard(s))" if pack else ""))
    deid.assert_clean(qc, "grid qc.json")
    _atomic(d / "qc.json", json.dumps(qc, separators=(",", ":")).encode("utf-8"))
    _atomic(d / "qc_flags.bin", flags.tobytes())
    if (d / "qc_seams.json").exists():
        (d / "qc_seams.json").unlink()                 # from an older export that ran the seam correction
    mat = write_mat_chunks(d, pos, Vt[:, embed], [wnames[k] for k in embed], bg[embed], hi[embed],
                           float(px_um), log=log) if embed else None

    entry = {"key": key, "label": f"{grid_um:g} µm squares · marker means", "kind": "protein",
             "binSizeUm": grid_um, "nBins": n, "nGenes": len(wnames), "micronsPerPixel": float(px_um),
             "spotDiameterPx": step_full, "imageScale": 1, "channels": [wnames[k] for k in embed],
             "geneSummary": "mean intensity per square (8-bit units, full-resolution pixels)",
             # exact (not rounded) so the square edges can be rebuilt: protgrid.edges(H, binSizeUm / levelPixelUm)
             "grid": {"level": g["level"], "levelPixelUm": g["pxUm"], "nx": nx, "ny": ny,
                      "tissueBy": by, "tissueThreshold": round(thr, 3), "tissueRule": rule,
                      "tissueSmoothUm": size * grid_um if size > 1 else 0, "tissueSmoothSquares": size,
                      "holesFilled": n_filled, "squares": nx * ny},
             "normalize": NORMALIZE,
             "qc": {"file": "qc.json", "flags": "qc_flags.bin", "everythingBright": n_bright, "saturatedSquares": n_sat,
                    "lowSignalSquares": int((~(x[:, embed] >= LOW_SIGNAL_X).any(1)).sum()) if embed else n,
                    "lowSignalRule": f"no embedded marker at x >= {LOW_SIGNAL_X:g}",
                    "weak": weak, "saturated": [c["name"] for c in qc_ch if c["saturated"]],
                    "seamsFlagged": bool(seams_qc["flagged"]),
                    "seamsCorrected": (seams_qc.get("correction") or {}).get("corrected", [])},
             "dropped": dropped}
    if mat:
        entry["mat"] = mat
    if pack:
        entry["genePack"] = pack

    if umap:
        feats = clr_features(Vt[:, embed], bg[embed], hi[embed]) if len(embed) >= CLR_MIN_CHANNELS else None
        keep = np.flatnonzero((x[:, embed] >= LOW_SIGNAL_X).any(1)) if embed else np.arange(0)
        log(f"[grid] {n - len(keep):,} of {n:,} squares low signal (no embedded marker at x >= {LOW_SIGNAL_X:g}): "
            "not embedded or clustered")
        # reuse an earlier umap/ only if it was made from these very numbers
        hh = hashlib.sha256(np.ascontiguousarray(feats).tobytes() if feats is not None else b"")
        hh.update(np.ascontiguousarray(keep.astype("<i8")).tobytes())
        h = hh.hexdigest()[:16]
        # markers: the embedded channels (a nuclear stain only says how many nuclei, not which cells)
        u = sectionumap.write_feature_umap(d, x, wnames, embed, embed, key, log=log, force=force,
                                           sample=umap_sample, features=feats, keep=keep,
                                           extra_params={"gridUm": grid_um, "tissueRule": rule,
                                                         "tissueSmoothUm": size * grid_um if size > 1 else 0,
                                                         "normalize": NORMALIZE, "cofactor": COFACTOR,
                                                         "embedFeatures": (f"centred log-ratio of max(v - bg, 0) / hi + {CLR_EPS:g}"
                                                                           if feats is not None else
                                                                           f"x (fewer than {CLR_MIN_CHANNELS} channels)"),
                                                         "lowSignal": f"no embedded marker at x >= {LOW_SIGNAL_X:g}: not embedded",
                                                         "input": h})
        if u:
            entry["umap"] = u
            if celltypes:
                usable = ~flat
                c = celltypes_mod.write_celltypes_protein(d, x, wnames, usable, markers=markers, key=key, log=log,
                                                          flags=flags)
                if c:
                    entry["umap"] = {**u, "celltypes": c}
    deid.assert_clean(entry, "grid manifest entry")
    log(f"[grid] {key} done  {time.time() - t_all:.0f}s")
    return entry
