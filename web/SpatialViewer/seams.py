"""Tile seams in a stitched multiplex image: find the camera-tile grid, and correct the square means.

A PhenoCycler / Akoya slide is photographed as a grid of camera fields ("tiles") and stitched.
Each tile can come out a little brighter or dimmer than its neighbours, separately for each
marker (each is imaged in its own cycle), so the image has sharp steps at tile edges that have
nothing to do with the tissue, and clusters can follow them. On s48 (16 µm squares) the column
means of CD45 jump 12.6% at a tile edge against 0.65% between other neighbouring columns, the
tiles are 960 x 720 px (486 x 365 µm), and circle -> UMAP clusters drawn there came out in bands.

1. Find the grid (find_grid), from the column / row means of every channel at full resolution
   (protgrid.channel_sums profiles): the sharpest jumps between neighbouring columns, the most
   common distance between them, and where they fall modulo that distance (the offset). A
   channel counts when the lattice phase + k * spacing is hit (+-3 px) at >= SEAM_COVERAGE of its
   positions; the slide has a grid along an axis when >= SEAM_CHANNELS channels agree on spacing
   (1%) and offset (3 px), and the offset found from each half of the slide (top / bottom for the
   columns, left / right for the rows) agrees too. Anatomy that repeats (fascicles) does not
   cover a whole lattice; a rescan or several pieces with different grids fail the halves test.

2. Read every channel with seams again at full resolution (channel_pieces): square sums split at
   the tile edges (so a square that straddles an edge can be corrected piece by piece), and thin
   strips (STRIP px) on both sides of every tile edge and of every mid-tile control line, tissue
   squares only, saturated pixels left out (clipped values do not scale).

3. Fit, per channel (fit_channel), one factor g and one offset c per tile, corrected = v / g - c:
   tissue is continuous across a tile edge, so the median log ratio of the pixels on the two
   sides of an edge between tiles a and b measures log g_b - log g_a; all edges together are one
   least-squares problem, each edge weighted by how consistent it is along its length, every tile
   with a prior of g within ~PRIOR_SD of 1, and g kept within GAIN_LIMITS (without the prior, small
   errors add up over ~1,200 tiles: factors drifted to 0.4-2.3x on s48). The offsets come the same way from the median differences after the gains.
   ("Each tile against its neighbourhood" was measured too and made the steps worse: tiles really
   differ in tissue.) Where a marker is smooth along an edge, a factor and an offset explain a step
   equally well, so the two are not separately identifiable there; either removes the step.

4. Keep a correction only if it works (choose): "excess step" = median |log ratio| across tile
   edges minus the same across mid-tile lines (tissue texture alone), in % points. Gain only and
   Shading shared by every tile (brighter at one side of every tile than at the other) puts a step
   of the same sign at every edge of one orientation; fitted as per-tile factors it would add up
   into a slope across the slide (on s48 96-98% of such a fit was that slope). So the median
   signed step per orientation is taken off the edge measurements before the fit, and reported:
   per-tile factors cannot remove that shading (a flat-field correction would; not built). Gain
   only and gain + offset are tried; the better one is kept if it at least halves the excess step
   and does not raise the mid-tile control by more than CONTROL_TOL points. Otherwise the channel
   stays uncorrected and is listed as not correctable, with the reason. The image is not changed.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage, sparse
from scipy.sparse.linalg import lsqr

JUMP_PCT = 99.7              # a "step" = a jump between neighbouring pixel columns above this percentile ...
MIN_JUMPS = 30               # ... but at least this many of the largest jumps
MIN_PX, MAX_PX = 128, 4000   # tile spacings looked for, in full-resolution pixels
REPEATS = 8                  # a channel shows steps when one spacing repeats this often between its jumps ...
COVERAGE = 0.5               # ... and the lattice phase + k * spacing is hit (+-3 px) at this share of its positions
CHANNELS = 2                 # the slide has a grid when this many channels agree (spacing 1%, offset 3 px)
PHASE_TOL = 3                # px
SEPARATION = 8               # px: of jumps closer than this only the largest counts
STRIP = 1                    # px on each side of a tile edge used to compare the two tiles: the
                             # pixels right at the edge are the most alike (on the demo's membrane
                             # marker the measured edge ratios follow the true ones r 0.93 at 1 px,
                             # 0.79 at 2, 0.59 at 4, 0.41 at 8) ...
BORDER, REF = 3, 3           # ... after taking off the border line: on s48 the last pixel column / row
                             # of every tile is dark (CD4 0.67 / 0.84 of the rest, CD45 0.93 / 0.95),
                             # the same in every tile. Each of the BORDER px next to an edge is
                             # compared with the REF px beyond it on the same side (the tile's own
                             # brightness cancels), over all edges; without this the edge pixels
                             # measure the line (steps of 20-55%, 97-99% one way), not the tiles.
OFFS = np.arange(-(BORDER + REF), BORDER + REF)   # pixel offsets from a line collected (0 = first px after it)
BORDER_MIN = 0.04            # a border factor within 4% of 1 is taken as no line (s48's lines: 5-33%)
OFFSET_GAIN = 0.8            # gain + offset is preferred over gain only when its excess step is this much smaller
SEG = 16                     # px: an edge is compared in pieces this long
MIN_SEGMENTS = 4             # pieces needed for one pair of tiles
PRIOR_SD = 0.15              # tiles are assumed within ~15% of each other (log factor sd) ...
PRIOR_SD_OFFSET = 2.0        # ... offsets: a wide prior (x the channel's median level at the edges):
                             # sparse markers sit near zero in tissue (s48 CD19 / CD38: median 0.31 of
                             # 255) and their seams are background offsets of ~0.2 between tiles,
                             # which per-channel scaling then magnifies; a prior of 0.1 x left most of
                             # the step (CD38 30 -> 25 points), 2 x removes it (-> 1.9)
GAIN_LIMITS = (0.7, 1.4)
MIN_EXCESS = 0.5             # % points: a smaller excess step is not worth correcting
KEEP = 0.5                   # the excess step must fall to at most this share
CONTROL_TOL = 0.3            # % points the mid-tile control may rise (a fall is no harm: an offset taken
                             # off near-zero values lowers it)
SMOOTH = 5                   # tiles: the share of a correction field that survives smoothing over this many
                             # tiles is reported (slow variation from tile to tile is allowed)
SHADING_SHARE = 0.7          # steps going the same way at this share of the edges of one orientation: shading
                             # shared by every tile (s48 after the border line: CD38 0.97, CD19 0.95, NF-L 0.80
                             # of vertical edges; CD68, PDL1 ~0.5)
OFFSET_ITERS = 3             # factor / offset fits alternated this many times


# ---- 1. the grid ----------------------------------------------------------------------------
def _jumps(prof) -> np.ndarray:
    """Positions p of the largest jumps (the step lies between pixel p - 1 and p), at most one per
    SEPARATION px (one boundary often spreads over a few neighbouring columns)."""
    g = np.abs(np.diff(np.asarray(prof, np.float64)))
    n = min(len(g), max(int(len(g) * (100 - JUMP_PCT) / 100), MIN_JUMPS))
    cand = np.argsort(-g, kind="stable")[:min(len(g), 20 * n)]
    taken = np.zeros(len(g) + 2 * SEPARATION, bool)
    out = []
    for i in cand:
        if not taken[i:i + 2 * SEPARATION + 1].any():
            out.append(i)
            taken[i + SEPARATION] = True
            if len(out) == n:
                break
    return np.sort(np.asarray(out, np.int64)) + 1


def _phase(pos: np.ndarray, sp: int) -> tuple[int, float]:
    """Offset of the lattice that the positions fall on most (+-PHASE_TOL), and that share."""
    h = np.bincount(pos % sp, minlength=sp).astype(np.float64)
    w = sum(np.roll(h, k) for k in range(-PHASE_TOL, PHASE_TOL + 1))
    ph0 = int(np.argmax(w))
    win = (ph0 + np.arange(-PHASE_TOL, PHASE_TOL + 1)) % sp       # the peak inside the best window
    ph = int(win[int(np.argmax(h[win]))])
    return ph, float(w[ph0] / max(len(pos), 1))


def _coverage(pos: np.ndarray, sp: int, ph: int, n: int) -> float:
    lat = np.arange(ph, n, sp)
    lat = lat[(lat > 0) & (lat < n)]
    if not len(lat):
        return 0.0
    d = np.abs(pos[None, :] - lat[:, None]).min(1)
    return float((d <= PHASE_TOL).mean())


def steps(prof) -> dict | None:
    """Sharp steps in one column (row) profile: the most common distance between any two of its
    largest jumps (to 4 px, refined to 1 px), the lattice offset, how many times that distance
    repeats and how much of the lattice is hit."""
    prof = np.asarray(prof, np.float64)
    n = len(prof)
    if n < 4 * MIN_PX:
        return None
    pos = _jumps(prof)
    dif = np.abs(pos[None, :] - pos[:, None])[np.triu_indices(len(pos), 1)]
    dif = dif[(dif >= MIN_PX) & (dif <= min(MAX_PX, n // 2))]
    if not len(dif):
        return None
    vals, cnt = np.unique(np.round(dif / 4).astype(np.int64) * 4, return_counts=True)
    # candidate spacings: the most common distances (to 4 px); each refined to the px whose lattice
    # catches the most jumps; the winner catches the most (a multiple of the spacing catches fewer)
    best = None
    for sp0 in vals[np.argsort(-cnt, kind="stable")[:6]]:
        for sp in range(max(MIN_PX, int(sp0) - 4), int(sp0) + 5):
            ph, share = _phase(pos, sp)
            if best is None or share > best[2] + 1e-9 or (abs(share - best[2]) <= 1e-9 and sp < best[0]):
                best = (sp, ph, share)
    sp, ph, _ = best
    reps = int(np.sum(np.abs(dif - sp) <= 2))
    return {"spacingPx": sp, "phasePx": ph, "repeats": reps, "coverage": round(_coverage(pos, sp, ph, n), 3),
            "jumps": int(len(pos))}


def _same(a: dict, b: dict) -> bool:
    if abs(a["spacingPx"] - b["spacingPx"]) > 0.01 * a["spacingPx"] + 1:
        return False
    d = abs(a["phasePx"] - b["phasePx"]) % a["spacingPx"]
    return min(d, a["spacingPx"] - d) <= PHASE_TOL


def find_grid(profiles: list, names: list[str], px_um: float) -> dict:
    """The seam check of the slide QC, and the tile grid if there is one.

    `profiles[k]` = (colMean, rowMean, colMeanTop, colMeanBottom, rowMeanLeft, rowMeanRight) of
    channel names[k] (protgrid.channel_sums). Per axis: the channels whose steps form a lattice,
    the largest group agreeing on spacing and offset, and whether each half of the slide puts the
    lattice in the same place."""
    out = {"rule": f"per channel, the sharpest jumps between neighbouring pixel columns / rows and the most "
                   f"common distance between them; a channel shows steps when that distance repeats >= {REPEATS} "
                   f"times and the lattice is hit at >= {COVERAGE:.0%} of its positions; flagged when >= "
                   f"{CHANNELS} channels agree on the spacing (1%) and offset ({PHASE_TOL} px)"}
    for axis, k, halves in (("columns", 0, (2, 3)), ("rows", 1, (4, 5))):
        hits = []
        for nm, pr in zip(names, profiles):
            st = steps(pr[k])
            if st and st["repeats"] >= REPEATS and st["coverage"] >= COVERAGE:
                hits.append((nm, st, pr))
        best = []
        for _nm, st, _pr in hits:
            grp = [h for h in hits if _same(st, h[1])]
            if len(grp) > len(best):
                best = grp
        res = {"channelsChecked": len(names), "channelsWithSteps": len(hits)}
        ok = len(best) >= CHANNELS
        if best:
            sp = int(np.median([h[1]["spacingPx"] for h in best]))
            ph = int(np.median([h[1]["phasePx"] for h in best]))
            # the halves: each half's jumps must sit on the same lattice
            agree = 0
            for _nm, _st, pr in best:
                hs = [_phase(_jumps(pr[i]), sp)[0] for i in halves]
                agree += all(min(abs(h - ph) % sp, sp - abs(h - ph) % sp) <= PHASE_TOL for h in hs)
            res.update(spacingPx=sp, spacingUm=round(sp * px_um), phasePx=ph, channels=[h[0] for h in best],
                       repeats=int(np.median([h[1]["repeats"] for h in best])),
                       coverage=round(float(np.median([h[1]["coverage"] for h in best])), 2),
                       halvesAgree=int(agree))
            res["consistent"] = bool(agree >= max(CHANNELS, len(best) // 2 + 1))
        res["flagged"] = bool(ok)
        out[axis] = res
    out["flagged"] = bool(out["columns"]["flagged"] or out["rows"]["flagged"])
    both = out["columns"]["flagged"] and out["rows"]["flagged"]
    consistent = both and out["columns"].get("consistent") and out["rows"].get("consistent")
    out["grid"] = ({"x0": out["columns"]["phasePx"], "y0": out["rows"]["phasePx"],
                    "w": out["columns"]["spacingPx"], "h": out["rows"]["spacingPx"]} if consistent else None)
    if out["flagged"] and not consistent:
        out["gridProblem"] = ("steps along one axis only" if not both else
                              "the two halves of the slide put the tile edges in different places "
                              "(several pieces or a rescan?)")
    return out


def edge_lines(n_px: int, sp: int, ph: int) -> tuple[np.ndarray, np.ndarray]:
    """Tile edges (phase + k * sp) and mid-tile control lines, far enough from the image edge for
    the pixels collected beside them (OFFS)."""
    lo, hi = -int(OFFS[0]), n_px - int(OFFS[-1]) - 1
    e = np.arange(ph % sp, n_px, sp)
    e = e[(e >= lo) & (e <= hi)]
    m = np.arange((ph + sp // 2) % sp, n_px, sp)
    m = m[(m >= lo) & (m <= hi)]
    return e, m


def tile_index(p: np.ndarray, sp: int, ph: int) -> np.ndarray:
    return np.floor((np.asarray(p) - ph) / sp).astype(np.int64)


# ---- 2. the second pass ---------------------------------------------------------------------
def channel_pieces(a, c: int, ye: np.ndarray, xe: np.ndarray, ye_r: np.ndarray, xe_r: np.ndarray, top: int,
                   vlines: np.ndarray, hlines: np.ndarray) -> dict:
    """One channel at full resolution: sums over the pieces of squares cut at the tile edges
    (edges ye_r x xe_r, which include the square edges ye x xe), and beside every vertical and
    horizontal line the pixels at each offset OFFS (sum and count of unsaturated pixels, per SEG px
    of the line)."""
    ny, nx = len(ye) - 1, len(xe) - 1
    nyr, nxr = len(ye_r) - 1, len(xe_r) - 1
    h_used, w_used = int(ye[-1]), int(xe[-1])
    S = np.zeros((nyr, nxr), np.float64)
    nv, nh = len(vlines), len(hlines)
    nby, nbx = (h_used + SEG - 1) // SEG, (w_used + SEG - 1) // SEG
    no = len(OFFS)
    VS, VN = np.zeros((nby, nv, no), np.float64), np.zeros((nby, nv, no), np.float64)
    HS, HN = np.zeros((nh, no, nbx), np.float64), np.zeros((nh, no, nbx), np.float64)
    row_bin = np.arange(h_used) // SEG
    xb = np.arange(0, w_used, SEG)
    vi = vlines[:, None] + OFFS[None, :]
    rows = a.chunks[1]
    buf, b0, j = None, 0, 0
    for y0 in range(0, h_used, rows):
        y1 = min(y0 + rows, h_used)
        band = np.asarray(a[c, y0:y1, :w_used])
        # strips beside the vertical lines, per square row
        if nv:
            sub = band[:, vi]                                        # rows x lines x offsets
            ok = sub < top
            np.add.at(VS, row_bin[y0:y1], np.where(ok, sub, 0).astype(np.float64))
            np.add.at(VN, row_bin[y0:y1], ok.astype(np.float64))
        # the rows at each offset from the horizontal lines, per SEG px
        for k, yl in enumerate(hlines):
            for oi, o in enumerate(OFFS):
                y = yl + o
                if y0 <= y < y1:
                    row = band[y - y0]
                    ok = row < top
                    HS[k, oi] += np.add.reduceat(np.where(ok, row, 0).astype(np.float64), xb)
                    HN[k, oi] += np.add.reduceat(ok.astype(np.float64), xb)
        # piece sums (as protgrid.channel_sums, on the finer edges)
        buf = band if buf is None or len(buf) == 0 else np.concatenate([buf, band])
        j1 = j
        while j1 < nyr and ye_r[j1 + 1] <= y1:
            j1 += 1
        if j1 > j:
            blk = buf[ye_r[j] - b0: ye_r[j1] - b0]
            rs = np.add.reduceat(blk, ye_r[j:j1] - ye_r[j], axis=0, dtype=np.uint64)
            S[j:j1] = np.add.reduceat(rs, xe_r[:-1], axis=1)
            buf = buf[ye_r[j1] - b0:]
            b0, j = int(ye_r[j1]), j1
    return {"S": S, "VS": VS, "VN": VN, "HS": HS, "HN": HN}


def fine_edges(e: np.ndarray, tile_edges: np.ndarray) -> np.ndarray:
    return np.unique(np.concatenate([e, tile_edges[(tile_edges > 0) & (tile_edges < e[-1])]]))


# ---- 3. the fit -------------------------------------------------------------------------------
def segments(P: dict, keep: np.ndarray, ye, xe, grid: dict, vlines, hlines, is_edge_v, is_edge_h):
    """Edge (and control) pieces with tissue on both sides: the mean at every offset OFFS, the two
    tiles, whether the piece is on a tile edge, and its orientation (0 = vertical line). A piece is
    SEG px of a line; it counts when the squares on both sides are tissue, it lies within one tile
    row (column), and at every offset at least half of its pixels are unsaturated."""
    w, h, x0, y0 = grid["w"], grid["h"], grid["x0"], grid["y0"]
    H, W = int(ye[-1]), int(xe[-1])
    ntc = int(tile_index(W - 1, w, x0) - tile_index(0, w, x0)) + 1
    tc0, tr0 = int(tile_index(0, w, x0)), int(tile_index(0, h, y0))

    def tid(tr, tc):
        return (np.asarray(tr) - tr0) * ntc + (np.asarray(tc) - tc0)

    def sq_col(x):
        return np.searchsorted(xe, x, side="right") - 1

    def sq_row(y):
        return np.searchsorted(ye, y, side="right") - 1

    V, A, B, E, O = [], [], [], [], []
    b0 = np.arange(0, H, SEG)
    b1 = np.minimum(b0 + SEG, H)                                    # pieces along a vertical line
    rs, blen = sq_row((b0 + b1) // 2), b1 - b0
    one_tr = tile_index(b0, h, y0) == tile_index(b1 - 1, h, y0)
    for v, x in enumerate(vlines):
        cl, cr = sq_col(x - 1), sq_col(x)
        ok = keep[rs, cl] & keep[rs, cr] & one_tr & (P["VN"][:, v, :] >= 0.5 * blen[:, None]).all(1)
        i_ = np.flatnonzero(ok)
        V.append(P["VS"][i_, v, :] / P["VN"][i_, v, :])
        trow = tile_index(b0[i_], h, y0)
        A.append(tid(trow, int(tile_index(x - 1, w, x0))))
        B.append(tid(trow, int(tile_index(x, w, x0))))
        E.append(np.full(len(i_), bool(is_edge_v[v])))
        O.append(np.zeros(len(i_), np.int8))
    c0 = np.arange(0, W, SEG)
    c1 = np.minimum(c0 + SEG, W)                                    # pieces along a horizontal line
    cs, clen = sq_col((c0 + c1) // 2), c1 - c0
    one_tc = tile_index(c0, w, x0) == tile_index(c1 - 1, w, x0)
    for k, y in enumerate(hlines):
        ra, rb = sq_row(y - 1), sq_row(y)
        ok = keep[ra, cs] & keep[rb, cs] & one_tc & (P["HN"][k] >= 0.5 * clen[None, :]).all(0)
        i_ = np.flatnonzero(ok)
        V.append((P["HS"][k][:, i_] / P["HN"][k][:, i_]).T)
        tcol = tile_index(c0[i_], w, x0)
        A.append(tid(int(tile_index(y - 1, h, y0)), tcol))
        B.append(tid(int(tile_index(y, h, y0)), tcol))
        E.append(np.full(len(i_), bool(is_edge_h[k])))
        O.append(np.ones(len(i_), np.int8))

    def cat(z, shape=(0,)):
        return np.concatenate(z) if z else np.zeros(shape)
    return (cat(V, (0, len(OFFS))), cat(A).astype(np.int64), cat(B).astype(np.int64), cat(E).astype(bool),
            cat(O).astype(np.int8))


def border_line(V: np.ndarray, E: np.ndarray, O: np.ndarray, eps: float) -> np.ndarray:
    """Factor per orientation and offset: the BORDER px on each side of a tile edge against the REF px
    beyond them on the same side (median over edge pieces); 1 for the reference px. (2, len(OFFS))."""
    F = np.ones((2, len(OFFS)))
    left_ref = (OFFS >= -(BORDER + REF)) & (OFFS < -BORDER)
    right_ref = (OFFS >= BORDER) & (OFFS < BORDER + REF)
    for o in (0, 1):
        sel = E & (O == o)
        if sel.sum() < 20:
            continue
        refL = V[sel][:, left_ref].mean(1)
        refR = V[sel][:, right_ref].mean(1)
        for k, d in enumerate(OFFS):
            if -BORDER <= d < 0:
                F[o, k] = np.exp(np.median(np.log((V[sel, k] + eps) / (refL + eps))))
            elif 0 <= d < BORDER:
                F[o, k] = np.exp(np.median(np.log((V[sel, k] + eps) / (refR + eps))))
    F[np.abs(np.log(F)) < BORDER_MIN] = 1.0
    return F


def _solve(pairs: dict, nt: int, scale: float = 1.0, prior_sd: float = PRIOR_SD) -> np.ndarray:
    """x_b - x_a = the median of each pair's values, weighted by 1 / its standard error (from the
    pair's own spread: 1.253 x 1.4826 x MAD / sqrt(n)); every x has a prior N(0, prior_sd); mean 0.
    Values, prior and x in units of `scale`."""
    if not pairs:
        return np.zeros(nt)
    rows, cols, vals, b = [], [], [], []
    for i, ((pa, pb), v) in enumerate(pairs.items()):
        v = np.asarray(v, np.float64) / scale
        med = float(np.median(v))
        mad = max(float(np.median(np.abs(v - med))) * 1.4826, 0.01)
        w = np.sqrt(len(v)) / (1.253 * mad)
        rows += [i, i]
        cols += [pb, pa]
        vals += [w, -w]
        b.append(w * med)
    k = len(pairs)
    rows += list(range(k, k + nt))
    cols += list(range(nt))
    vals += [1.0 / prior_sd] * nt
    b += [0.0] * nt
    A = sparse.csr_matrix((vals, (rows, cols)), shape=(k + nt, nt))
    x = lsqr(A, np.asarray(b), atol=1e-10, btol=1e-10)[0] * scale
    used = np.zeros(nt, bool)
    used[np.unique(np.array(list(pairs)).ravel())] = True
    x = x - x[used].mean()
    x[~used] = 0.0                                   # tiles with no usable edge: no change
    return x


def _pairs(a, b, v, edge) -> dict:
    out: dict = {}
    for pa, pb, val in zip(a[edge], b[edge], v[edge]):
        out.setdefault((int(pa), int(pb)), []).append(val)
    return {k: v for k, v in out.items() if len(v) >= MIN_SEGMENTS}


def _excess(L, R, A, B, E, eps) -> tuple[float, float]:
    """Median over tile pairs of |median log ratio across their edge| (texture cancels within a
    pair; a tile step does not), for tile edges and for mid-tile control lines, in %."""
    r = np.log((np.maximum(R, 0) + eps) / (np.maximum(L, 0) + eps))
    out = []
    for sel in (E, ~E):
        grp = _pairs(A, B, r, sel)
        out.append(float(np.median([abs(np.median(v)) for v in grp.values()]) * 100) if grp else 0.0)
    return out[0], out[1]


def _smooth(x: np.ndarray, used: np.ndarray, shape: tuple[int, int]) -> np.ndarray:
    """x averaged over SMOOTH x SMOOTH tiles (fitted tiles only)."""
    w = used.reshape(shape).astype(np.float64)
    return (ndimage.uniform_filter(x.reshape(shape) * w, SMOOTH, mode="constant") / np.maximum(
        ndimage.uniform_filter(w, SMOOTH, mode="constant"), 1e-9)).ravel()


def lowfreq_share(x: np.ndarray, used: np.ndarray, shape: tuple[int, int]) -> float:
    """sd of x smoothed over SMOOTH x SMOOTH tiles / sd of x, over the fitted tiles."""
    if used.sum() < 10 or np.std(x[used]) == 0:
        return 0.0
    return float(np.std(_smooth(x, used, shape)[used]) / np.std(x[used]))


def fit_channel(V, A, B, E, O, nt: int, shape: tuple[int, int] | None = None) -> dict:
    """Per-tile factors g and offsets c for one channel, and the before / after excess steps.
    `V` = pieces x OFFS means (segments); the border line (border_line) is divided out of the
    edge pieces first, then the STRIP px on each side of the line are compared."""
    eps = max(float(np.percentile(V, 5)) if V.size else 1.0, 0.5)
    F = border_line(V, E, O, eps)
    Vc = V / np.where(E[:, None], F[O], 1.0)                  # mid-tile lines have no border line
    li, ri = (OFFS >= -STRIP) & (OFFS < 0), (OFFS >= 0) & (OFFS < STRIP)
    L, R = Vc[:, li].mean(1), Vc[:, ri].mean(1)
    L0, R0 = V[:, li].mean(1), V[:, ri].mean(1)
    e_raw, _ = _excess(L0, R0, A, B, E, eps)
    e0, c0 = _excess(L, R, A, B, E, eps)
    res = {"excessBefore": round(e0 - c0, 2), "edgeBefore": round(e0, 2), "controlBefore": round(c0, 2),
           "edgeRawPixels": round(e_raw, 2), "segments": int(E.sum()),
           "borderLine": {ax: {"before": [round(float(v), 3) for v in F[o, (OFFS >= -BORDER) & (OFFS < 0)]],
                               "after": [round(float(v), 3) for v in F[o, (OFFS >= 0) & (OFFS < BORDER)]]}
                          for o, ax in ((0, "columns"), (1, "rows"))}}
    if not E.any() or not (~E).any():
        res.update(method="none", reason="too few tile edges with tissue on both sides")
        return res
    s = float(np.median(np.r_[L, R])) or 1.0
    # shading shared by every tile: the median signed step per orientation, taken off before the fit
    lr0 = np.log((R + eps) / (L + eps))
    shading = {}
    for o, ax in ((0, "columns"), (1, "rows")):
        pr = _pairs(A, B, lr0, E & (O == o))
        med = np.array([np.median(v) for v in pr.values()]) if pr else np.zeros(0)
        shading[ax] = {"stepsUp": round(float((med > 0).mean()), 3) if len(med) else None,
                       "medianStepPct": round(float(np.median(med) * 100), 2) if len(med) else None,
                       # taken off only when the steps really go one way (the median of steps that
                       # alternate up and down is just one of the two)
                       "_b": float(np.median(med)) if len(med) and max((med > 0).mean(), (med < 0).mean())
                       >= SHADING_SHARE else 0.0}
    res["shading"] = {ax: {k: v for k, v in d.items() if k != "_b"} for ax, d in shading.items()}
    bias = np.where(O == 0, shading["columns"]["_b"], shading["rows"]["_b"]) * E
    g_only = None
    g, c = np.ones(nt), np.zeros(nt)
    for it in range(OFFSET_ITERS):                # v = g (u + c): factors on v - g c, then offsets on v / g
        Lu = np.maximum(L - g[A] * c[A], 0.0)
        Ru = np.maximum(R - g[B] * c[B], 0.0)
        g = np.clip(np.exp(_solve(_pairs(A, B, np.log((Ru + eps) / (Lu + eps)) - bias, E), nt)), *GAIN_LIMITS)
        if g_only is None:
            g_only = g
        d = R / g[B] - L / g[A]
        dbias = {o: float(np.median(d[E & (O == o)])) if (E & (O == o)).any() and
                 shading["columns" if o == 0 else "rows"]["_b"] != 0.0 else 0.0 for o in (0, 1)}
        c = _solve(_pairs(A, B, d - np.where(O == 0, dbias[0], dbias[1]) * E, E), nt, scale=s,
                   prior_sd=PRIOR_SD_OFFSET)
    used = np.zeros(nt, bool)
    used[np.unique(np.r_[A[E], B[E]])] = True
    tries = {"gain": (g_only, np.zeros(nt)), "gain+offset": (g, c)}
    for name, (gg, cc) in tries.items():
        e1, c1 = _excess(L / gg[A] - cc[A], R / gg[B] - cc[B], A, B, E, eps)
        res[name] = {"excessAfter": round(e1 - c1, 2), "edgeAfter": round(e1, 2), "controlAfter": round(c1, 2)}
        if shape is not None:
            res[name]["slowShare"] = round(max(lowfreq_share(np.log(gg), used, shape),
                                               lowfreq_share(cc / s, used, shape) if np.any(cc) else 0.0), 3)
    res["_g"], res["_c"], res["_g_only"] = g, c, g_only
    return res


def choose(res: dict) -> dict:
    """Keep the better of gain / gain+offset if it works (module docstring, step 4)."""
    ex0 = res["excessBefore"]
    if res.get("method") == "none":
        return res
    if ex0 < MIN_EXCESS:
        res.update(method="none", reason=f"excess step {ex0:.2f} points: nothing to correct")
        return res
    ok = []
    for name in ("gain", "gain+offset"):
        t = res[name]
        why = []
        if t["excessAfter"] > KEEP * ex0:
            why.append(f"excess step only {ex0:.2f} -> {t['excessAfter']:.2f}")
        if t["controlAfter"] - res["controlBefore"] > CONTROL_TOL:
            why.append(f"mid-tile control rose {res['controlBefore']:.2f} -> {t['controlAfter']:.2f}")
        if why:
            t["rejected"] = "; ".join(why)
        else:
            ok.append(name)
    sh = res.get("shading") or {}
    up = [(ax, v["stepsUp"], v["medianStepPct"]) for ax, v in sh.items()
          if v.get("stepsUp") is not None and max(v["stepsUp"], 1 - v["stepsUp"]) >= SHADING_SHARE]
    if not ok:
        why = " / ".join(res[n]["rejected"] for n in ("gain", "gain+offset"))
        if up:
            why = ("shading shared by every tile (the step goes the same way at " +
                   ", ".join(f"{max(u, 1 - u):.0%} of the {ax} edges, median {m:+.0f}%" for ax, u, m in up) +
                   "), which per-tile factors cannot remove; " + why)
        res.update(method="none", reason="not correctable: " + why)
    else:
        res["method"] = ("gain+offset" if "gain+offset" in ok and ("gain" not in ok or res["gain+offset"]["excessAfter"]
                                                                   < OFFSET_GAIN * res["gain"]["excessAfter"]) else "gain")
        res["excessAfter"] = res[res["method"]]["excessAfter"]
    return res


# ---- 4. applying it to the squares ------------------------------------------------------------
def piece_maps(ye, xe, ye_r, xe_r, grid: dict):
    """For each piece row / column: its square row / column and tile row / column (pieces never cross)."""
    w, h, x0, y0 = grid["w"], grid["h"], grid["x0"], grid["y0"]
    sr = np.searchsorted(ye, ye_r[:-1], side="right") - 1
    sc = np.searchsorted(xe, xe_r[:-1], side="right") - 1
    tr = tile_index(ye_r[:-1], h, y0)
    tc = tile_index(xe_r[:-1], w, x0)
    ntc = int(tile_index(xe[-1] - 1, w, x0) - tile_index(0, w, x0)) + 1
    tid = (tr - tile_index(0, h, y0))[:, None] * ntc + (tc - tile_index(0, w, x0))[None, :]
    return sr, sc, tid


def corrected_sums(S: np.ndarray, ye_r, xe_r, sr, sc, tid, g, c, ny: int, nx: int) -> np.ndarray:
    """Square sums of v / g - c, from the piece sums."""
    npx = np.outer(np.diff(ye_r), np.diff(xe_r)).astype(np.float64)
    P = S / g[tid] - c[tid] * npx
    rs = np.add.reduceat(P, np.searchsorted(sr, np.arange(ny)), axis=0)
    return np.add.reduceat(rs, np.searchsorted(sc, np.arange(nx)), axis=1)
