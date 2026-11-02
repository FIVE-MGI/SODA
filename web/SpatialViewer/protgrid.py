"""Marker expression on a grid: square means of a multiplexed fluorescence image as an expr layer.

Akoya / PhenoCycler (and any multi-channel fluorescence TIFF) has no count matrix, so
there is nothing for the gene search, the whole-section UMAP or the clusters to read.
This lays a square grid over the image, takes each marker channel's mean intensity per
square, and writes it as an ordinary expression layer of the image dataset, so all of
that machinery works unchanged (cell segmentation is a later, separate step):

    expr/grid<NNN>um/positions.bin     float32 x,y  square centres in full-res image px (tissue squares only)
    expr/grid<NNN>um/genes.json        [{n, t}]     one per marker channel written; t = mean over the squares
    expr/grid<NNN>um/gene/<NAME>.bin   float32 x,y,v  EVERY square, in positions.bin order (v may be 0;
                                                    the viewer treats v <= 0 as "nothing here", as for RNA)
    expr/grid<NNN>um/umap/...          as sectionumap.py (index.json, coords.bin, clusters_r*.bin)
                                       + celltypes.json / celltype.bin (celltypes.py, protein rule)
    No hvg/ (the circle -> UMAP tool is for counts).

v is in the 8-bit units of the tiles: the raw mean for a uint8 pyramid, the raw mean
mapped through the channel's `bits` window for uint16 (not clipped at 255).

Squares: `grid_um` wide (default 20), starting at the image's top-left corner; a last
partial row / column is kept when it is at least half a square. The means come from the
pyramid level whose pixel is closest to grid/8 and never coarser than grid/4, streamed one
chunk row at a time (each chunk is read once per channel; the full-res level is never loaded
unless it is that level).

Tissue: the DAPI (or Hoechst) square means (without a nuclear stain, the mean over the written
channels) are first averaged over a window of `tissue_smooth_um` (default 100 µm = 5 x 5 squares
of 20 µm; moving average, edges by nearest-square replication); a square is kept when that value
is above an Otsu threshold over all squares, taken on log1p (see tissue_mask for why not linear).
Enclosed holes of <= HOLE_MAX squares are filled. tissue_smooth_um = 0 thresholds the raw square
means (the rule before 2026-11: it keeps only nucleus-dense squares, a lacy mask that misses axon /
collagen stretches). The default is the validation study's "area mask"
(D:/SpatialVizValidation/scripts/make_tiles.py), square for square.

Channels: names matching the viewer's hidden-channel default (open / blank / empty / unused)
are dropped altogether. Written but left out of the embedding: nuclear stains, and channels
that are flat over the tissue (p99.5 - p1 of the square means < FLAT_RANGE units).
Embedding (sectionumap.write_feature_umap): x' = clip((v - p1) / (p99.5 - p1), 0, 1) per channel
over the tissue squares -> z-score -> PCA min(15, channels - 1) -> UMAP (15 neighbours, min_dist 0.3,
seed 42) -> Leiden 0.3 / 0.6 / 1.0, cluster 0 = largest. Markers: difference of the mean z of x'
(cluster - rest); pct = share of the cluster's squares whose x' is above the channel's median.
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

DEFAULT_GRID_UM = 20.0
NUCLEAR = re.compile(r"^\s*(dapi|hoechst)\b", re.I)                    # = viewer.js NUCLEAR_CHANNEL
HIDDEN = re.compile(r"^\s*(open|blank|empty|unused)([\s_-]|$)", re.I)   # = viewer.js HIDE_CHANNEL default
FLAT_RANGE = 2.0              # p99.5 - p1 of the tissue squares, 8-bit units: below this a channel is flat
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
    """Image datasets only: not Visium (it has real counts), not RNA-seq, not a colour photo."""
    return (m.get("kind") not in ("visium", "rnaseq") and not m.get("expr") and "image" in m
            and not is_photo(m["image"].get("channels", [])))


def pick_level(shapes, width: int, px_um: float, grid_um: float) -> tuple[int, float]:
    """(level, µm per pixel at that level): pixel closest to grid/8 (log scale), never coarser than grid/4."""
    best = None
    for L, (_h, w) in enumerate(shapes):
        p = px_um * width / w
        if p > grid_um / 4 * (1 + 1e-9):
            continue
        d = abs(math.log(p / (grid_um / 8)))
        if best is None or d < best[0] - 1e-12:
            best = (d, L, p)
    if best is None:                                  # even full res is coarser than grid/4
        return 0, px_um
    return best[1], best[2]


def edges(n_px: int, step: float) -> np.ndarray:
    """Pixel edges of the squares along one axis (a last partial square is kept if >= half a square)."""
    n = int(math.floor(n_px / step + 0.5))
    return np.minimum(np.round(np.arange(n + 1) * step).astype(np.int64), n_px)


def channel_sums(a, c: int, ye: np.ndarray, xe: np.ndarray) -> np.ndarray:
    """Sum of channel `c` of zarr array `a` (C, H, W) over every square, one chunk row at a time."""
    ny, nx = len(ye) - 1, len(xe) - 1
    S = np.zeros((ny, nx), np.float64)
    rows = a.chunks[1]
    h_used, w_used = int(ye[-1]), int(xe[-1])
    buf, b0, j = None, 0, 0
    for y0 in range(0, h_used, rows):
        y1 = min(y0 + rows, h_used)
        band = np.asarray(a[c, y0:y1, :w_used])
        buf = band if buf is None or len(buf) == 0 else np.concatenate([buf, band])
        j1 = j
        while j1 < ny and ye[j1 + 1] <= y1:
            j1 += 1
        if j1 > j:
            blk = buf[ye[j] - b0: ye[j1] - b0]
            rs = np.add.reduceat(blk, ye[j:j1] - ye[j], axis=0, dtype=np.uint64)
            S[j:j1] = np.add.reduceat(rs, xe[:-1], axis=1)
            buf = buf[ye[j1] - b0:]
            b0, j = int(ye[j1]), j1
    return S


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
    """x' = clip((v - p1) / (p99.5 - p1), 0, 1) per column; also returns p1 and p99.5."""
    V = np.asarray(V, np.float64)
    lo = np.percentile(V, P_LO, axis=0)
    hi = np.percentile(V, P_HI, axis=0)
    rng = hi - lo
    xp = np.where(rng > 0, (V - lo) / np.where(rng > 0, rng, 1.0), 0.0)
    return np.clip(xp, 0.0, 1.0), lo, hi


def square_means(r, grid_um: float, px_um: float, log=print, workers: int = 4) -> dict:
    """Mean per square and channel (8-bit units) over the whole grid, from the chosen level."""
    L, px_L = pick_level(r.shapes, r.W, px_um, grid_um)
    H_L, W_L = r.shapes[L]
    step_L = grid_um / px_L
    if step_L < 2:
        raise ValueError(f"{grid_um} µm squares are {step_L:.1f} px at full resolution: too small")
    if px_L > grid_um / 4:
        log(f"[grid] WARNING: full-res pixels ({px_um:.2f} µm) are coarser than grid/4")
    ye, xe = edges(H_L, step_L), edges(W_L, step_L)
    ny, nx = len(ye) - 1, len(xe) - 1
    a = r._arrays[L]
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=workers) as ex:
        sums = list(ex.map(lambda c: channel_sums(a, c, ye, xe), range(len(r.channels))))
    npx = np.outer(np.diff(ye), np.diff(xe)).astype(np.float64)
    V = np.empty((ny * nx, len(r.channels)), np.float32)
    V64 = np.empty((ny * nx, len(r.channels)), np.float64)   # for the tissue mask (not rounded to float32)
    for c, S in enumerate(sums):
        mean = (S / npx).ravel()
        b0, b1 = (float(x) for x in r.channels[c]["bits"])
        if not (str(a.dtype) == "uint8" and b0 == 0 and b1 == 255):
            mean = np.maximum((mean - b0) * (255.0 / max(b1 - b0, 1e-9)), 0.0)
        V[:, c] = mean
        V64[:, c] = mean
    log(f"[grid] level {L} ({px_L:.3f} µm/px, {W_L}x{H_L}), {step_L:.2f} px squares, {nx}x{ny} = {nx * ny:,} "
        f"squares x {len(r.channels)} channels  {time.time() - t0:.0f}s")
    return {"V": V, "V64": V64, "nx": nx, "ny": ny, "level": L, "pxUm": px_L}


def _atomic(p: Path, buf: bytes):
    tmp = p.with_name(p.name + ".part")
    tmp.write_bytes(buf)
    tmp.replace(p)


def write_grid_layer(r, m: dict, out: Path, grid_um: float = DEFAULT_GRID_UM, markers: Path | None = None,
                     umap: bool = True, celltypes: bool = True, log=print, force: bool = False,
                     tissue_smooth_um: float = DEFAULT_TISSUE_SMOOTH_UM) -> dict | None:
    """Write expr/<key>/ for an image dataset (PyramidReader `r`, manifest `m`) into export folder `out`.

    `tissue_smooth_um`: moving-average window (µm) over the square means before the tissue
    threshold; 0 = threshold the raw means (the old, lacy mask).

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
    g = square_means(r, grid_um, float(px_um), log=log)
    V, V64, nx, ny = g["V"], g["V64"], g["nx"], g["ny"]
    names = [str(ch["name"]) for ch in r.channels]

    # ---- channels written / embedded --------------------------------------------------
    written = [j for j, n in enumerate(names) if not HIDDEN.match(n)]
    dropped = [{"name": names[j], "reason": "open / blank channel: not written"}
               for j in range(len(names)) if j not in written]
    if not written:
        log("[grid] skipped: every channel is an open / blank channel")
        return None

    # ---- tissue mask ---------------------------------------------------------------------
    nuc = [j for j in written if NUCLEAR.match(names[j])]
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

    Vt = V[rows][:, written].astype(np.float64)       # tissue squares x written channels
    wnames = [names[j] for j in written]
    xp, p1, p995 = normalise(Vt)
    flat = (p995 - p1) < FLAT_RANGE
    embed, marker_cols = [], []
    for k, j in enumerate(written):
        if flat[k]:
            dropped.append({"name": names[j], "reason": f"flat over the tissue (p99.5 - p1 = {p995[k] - p1[k]:.2f}): "
                                                         "written, not embedded"})
            continue
        marker_cols.append(k)
        if NUCLEAR.match(names[j]):
            dropped.append({"name": names[j], "reason": "nuclear stain: written (and the tissue mask), not embedded"})
        else:
            embed.append(k)
    log(f"[grid] {len(written)} channels written, {len(embed)} embedded: " + ", ".join(wnames[k] for k in embed))
    for dr in dropped:
        log(f"[grid]   left out: {dr['name']} - {dr['reason']}")

    # ---- files ---------------------------------------------------------------------------
    step_full = grid_um / float(px_um)
    iy, ix = np.divmod(rows, nx)
    pos = np.empty((n, 2), np.float32)
    pos[:, 0] = (ix + 0.5) * step_full
    pos[:, 1] = (iy + 0.5) * step_full
    d = Path(out) / "expr" / key
    (d / "gene").mkdir(parents=True, exist_ok=True)
    _atomic(d / "positions.bin", np.ascontiguousarray(pos, dtype="<f4").tobytes())
    genes = [{"n": nm, "t": round(float(Vt[:, k].mean()), 4)} for k, nm in enumerate(wnames)]
    deid.assert_clean(genes, "grid genes.json")
    _atomic(d / "genes.json", json.dumps(genes, separators=(",", ":")).encode("utf-8"))
    xyv = np.empty((n, 3), np.float32)
    xyv[:, :2] = pos
    for k, nm in enumerate(wnames):
        xyv[:, 2] = Vt[:, k]
        _atomic(d / "gene" / f"{_encode(nm)}.bin", np.ascontiguousarray(xyv, dtype="<f4").tobytes())
    log(f"[grid] {key}: {n:,} squares, {len(wnames)} gene files -> {d}")

    entry = {"key": key, "label": f"{grid_um:g} µm squares · marker means", "kind": "protein",
             "binSizeUm": grid_um, "nBins": n, "nGenes": len(wnames), "micronsPerPixel": float(px_um),
             "spotDiameterPx": step_full, "imageScale": 1, "channels": [wnames[k] for k in embed],
             "geneSummary": "mean intensity per square (8-bit units)",
             # exact (not rounded) so the square edges can be rebuilt: protgrid.edges(H, binSizeUm / levelPixelUm)
             "grid": {"level": g["level"], "levelPixelUm": g["pxUm"], "nx": nx, "ny": ny,
                      "tissueBy": by, "tissueThreshold": round(thr, 3), "tissueRule": rule,
                      "tissueSmoothUm": size * grid_um if size > 1 else 0, "tissueSmoothSquares": size,
                      "holesFilled": n_filled, "squares": nx * ny},
             "dropped": dropped}

    if umap:
        # reuse an earlier umap/ only if it was made from these very numbers
        h = hashlib.sha256(np.ascontiguousarray(xp[:, embed]).tobytes()).hexdigest()[:16]
        u = sectionumap.write_feature_umap(d, xp, wnames, embed, marker_cols, key, log=log, force=force,
                                           extra_params={"gridUm": grid_um, "tissueRule": rule,
                                                         "tissueSmoothUm": size * grid_um if size > 1 else 0,
                                                         "input": h})
        if u:
            entry["umap"] = u
            if celltypes:
                usable = ~flat
                c = celltypes_mod.write_celltypes_protein(d, xp, wnames, usable, markers=markers, key=key, log=log)
                if c:
                    entry["umap"] = {**u, "celltypes": c}
    deid.assert_clean(entry, "grid manifest entry")
    log(f"[grid] {key} done  {time.time() - t_all:.0f}s")
    return entry
