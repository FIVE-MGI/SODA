"""Tile seams in a stitched multiplex image: find the camera-tile grid.

A PhenoCycler / Akoya slide is photographed as a grid of camera fields ("tiles") and stitched.
Each tile can come out a little brighter or dimmer than its neighbours, separately for each
marker (each is imaged in its own cycle), so the image has sharp steps at tile edges that have
nothing to do with the tissue, and clusters can follow them. On s48 (16 µm squares) the column
means of CD45 jump 12.6% at a tile edge against 0.65% between other neighbouring columns, the
tiles are 960 x 720 px (486 x 365 µm). (On s48 the clusters did not follow the tiles: the
dependence of cluster labels on the position within a tile was no higher than on a fake 1,000 x
750 px grid, over the whole slide and in 12 circle-sized regions; bands once seen in a circle ->
UMAP were ~1.5-2 tiles apart, not tile-aligned.)

1. Find the grid (find_grid), from the column / row means of every channel at full resolution
   (protgrid.channel_sums profiles): the sharpest jumps between neighbouring columns, the most
   common distance between them, and where they fall modulo that distance (the offset). A
   channel counts when the lattice phase + k * spacing is hit (+-3 px) at >= SEAM_COVERAGE of its
   positions; the slide has a grid along an axis when >= SEAM_CHANNELS channels agree on spacing
   (1%) and offset (3 px), and the offset found from each half of the slide (top / bottom for the
   columns, left / right for the rows) agrees too. Anatomy that repeats (fascicles) does not
   cover a whole lattice; a rescan or several pieces with different grids fail the halves test.

The seams are reported in the slide QC (qc.json "seams"); the image and the square means are not changed.
"""
from __future__ import annotations

import numpy as np

JUMP_PCT = 99.7              # a "step" = a jump between neighbouring pixel columns above this percentile ...
MIN_JUMPS = 30               # ... but at least this many of the largest jumps
MIN_PX, MAX_PX = 128, 4000   # tile spacings looked for, in full-resolution pixels
REPEATS = 8                  # a channel shows steps when one spacing repeats this often between its jumps ...
COVERAGE = 0.5               # ... and the lattice phase + k * spacing is hit (+-3 px) at this share of its positions
CHANNELS = 2                 # the slide has a grid when this many channels agree (spacing 1%, offset 3 px)
PHASE_TOL = 3                # px
SEPARATION = 8               # px: of jumps closer than this only the largest counts


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
