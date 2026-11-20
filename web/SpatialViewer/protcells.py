"""A lab's own cell segmentation of a protein image -> a "cells" layer. SODA segments nothing.

    prepare.py tiff | akoya | omezarr ... --cells <label mask .tif | cells.geojson>
                                          [--nuclei <label mask .tif | nuclei.geojson>] [--cells-units px|um]

What is read (--cells, and --nuclei for nucleus outlines only):
  label mask TIFF  integer labels, 0 = background, one number per cell: Cellpose (*_cp_masks.tif), Mesmer /
                   DeepCell, StarDist, QuPath's labelled-image export. It is the image's full-resolution size, or
                   a whole-number downscale of it (s = 2, 3, 4 ...: each mask pixel covers s x s image pixels, the
                   mask's size being the image's divided by s, rounded either way). Any other size is refused.
  GeoJSON          Polygon / MultiPolygon features in full-resolution pixels of the image (QuPath "Export objects
                   as GeoJSON"; --cells-units um: µm from the image's top-left corner), read with
                   cellshapes.read_geojson. Features whose objectType is annotation / tile / TMA core are skipped;
                   a feature without an id gets its number. QuPath's nucleusGeometry gives the nucleus outlines
                   when --nuclei is not given.
Both are in the pixels of the image file as it was segmented; --crop and --flip are applied to them as to the image.

Per cell, from the cache pyramid's full-resolution level (image.zarr level 0, one tile row at a time, each chunk
read once per channel, never the whole image): its pixels, area, centroid (the mean of its pixel centres, in
full-resolution image pixels) and the mean of every channel over its pixels, in the 8-bit units of the grid
layer's squares (protgrid.to_8bit: the raw mean through the channel's `bits` window, not clipped at 255). A pixel
belongs to a cell when its mask label is the cell's (a downscaled mask: the label of the mask pixel it falls in),
or for GeoJSON when its centre lies inside the cell's polygon (even-odd; a pixel inside the polygons of two
different cells belongs to neither). A cell with no full-resolution pixel in the image is left out.

Cache:
    expr/cells.zarr      as expr.py's layers: xy = centroids, genes = the image's channels, values = the means
                         (CSR, zeros not stored), barcodes = the mask label / the feature id, + area (pixels)
    shapes/cells/...     cellshapes.py: outlines traced from the mask (trace_labels, in windows of TRACE_TILE mask
                         pixels with a TRACE_HALO border; a cell larger than the border stays a dot) or the polygons
    manifest expr[]      {"key": "cells", "kind": "cells", "modality": "protein", "segmentation": {...}, ...}

Export (export_layer, called by export.py after the grid layer): the grid layer's rules (protgrid.py), on cells
instead of squares. Channels: open / blank names not written, duplicates reduced to the stronger, flat channels and
the nuclear stain written but not embedded. bg = the grid layer's background of that channel (the median of its
off-tissue squares: the slide's own background, which cells, all on tissue, cannot give), else the 1st percentile
over the cells; hi = p99.5 of v - bg over the cells, leaving out the everything-bright ones;
x = asinh(max(v - bg, 0) / hi / 0.2) / asinh(5). UMAP and Leiden on the centred log-ratio (10 or more embedded
channels, else x) of the cells with an embedded marker at x >= 0.2 (sectionumap.write_feature_umap), mat/ chunks
for Draw a region (protgrid.write_mat_chunks), suggested names by the protein rule (celltypes.py), and one gene
file per channel holding every cell in positions.bin order (v may be 0), as for the grid's squares.
"""
from __future__ import annotations

import hashlib
import math
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

import cellshapes

KEY = "cells"
LABEL = "Cells (own segmentation)"
TRACE_TILE = 2048                # mask pixels per side of a tracing window (its core)
TRACE_HALO = 256                 # mask pixels read around the core: a cell up to this size is traced whole
BAND_PX = 1 << 22                # pixels per label slice of the measuring pass
LUT_MAX = 1 << 26                # labels up to this: a lookup table; above it: a sorted search
SKIP_TYPES = {"annotation", "tile", "tma_core", "tmacore", "root"}   # QuPath objectType values that are no cell
MASK_EXT = (".tif", ".tiff")
GEOJSON_EXT = (".geojson", ".json")
UNITS = ("px", "um")


def add_args(p) -> None:
    """--cells / --nuclei / --cells-units of `prepare.py tiff | akoya | omezarr`."""
    p.add_argument("--cells", help="your own cell segmentation of this image (SODA segments nothing): a label mask TIFF "
                                   "(0 = background, one number per cell: Cellpose, Mesmer / DeepCell, StarDist, QuPath; "
                                   "the image's full-resolution size or a whole-number downscale of it) or a GeoJSON of "
                                   "cell polygons (QuPath) in full-resolution pixels. Adds a 'cells' layer: per cell its "
                                   "area, centroid, the mean of every channel, outlines; at export UMAP and clusters as "
                                   "for the grid squares (protcells.py)")
    p.add_argument("--nuclei", help="with --cells: nucleus outlines, a label mask or GeoJSON with the same ids as --cells "
                                    "(optional; QuPath's nucleusGeometry in the --cells GeoJSON is read without it)")
    p.add_argument("--cells-units", choices=UNITS, default="px",
                   help="coordinates of a GeoJSON --cells / --nuclei: px = full-resolution image pixels (default), "
                        "um = µm from the image's top-left corner")


# ------------------------------------------------------------------------------------- geometry
class Frame:
    """The cache image's pixels (after --crop and --flip) against the pixels of the image file as segmented."""

    def __init__(self, m: dict):
        img = m["image"]
        self.W, self.H = int(img["width"]), int(img["height"])
        crop = m.get("crop") or {}
        box = crop.get("box") or [0, 0, self.W, self.H]
        self.x0, self.y0 = int(box[0]), int(box[1])
        full = crop.get("ofFullSize") or [self.W, self.H]
        self.full_w, self.full_h = int(full[0]), int(full[1])
        flip = m.get("flip", "") or ""
        self.fh, self.fv = "h" in flip, "v" in flip

    def src_cols(self) -> np.ndarray:
        """Column of the original image under every column of the frame."""
        x = np.arange(self.W, dtype=np.int64)
        return (self.W - 1 - x if self.fh else x) + self.x0

    def src_rows(self, y0: int, y1: int) -> np.ndarray:
        y = np.arange(y0, y1, dtype=np.int64)
        return (self.H - 1 - y if self.fv else y) + self.y0

    def to_frame(self, xy: np.ndarray) -> np.ndarray:
        """Points (continuous, pixel corners at integers) of the original image -> the frame."""
        x, y = xy[:, 0] - self.x0, xy[:, 1] - self.y0
        if self.fh:
            x = self.W - x
        if self.fv:
            y = self.H - y
        return np.column_stack([x, y]).astype(np.float64)


def mask_scale(hw_mask, hw_image, what: str) -> int:
    """The whole-number downscale s of a mask against the image (1 = the same size), or SystemExit."""
    h, w = (int(v) for v in hw_mask)
    H, W = (int(v) for v in hw_image)
    if (h, w) == (H, W):
        return 1
    s = int(round(W / w)) if w else 0
    if s >= 2 and W // s <= w <= -(-W // s) and H // s <= h <= -(-H // s):
        return s
    hint = " (it is the image's size turned by 90 degrees: save it the other way round)" if (h, w) == (W, H) else ""
    near = int(round(W / w)) if w else 0
    sug = f"; a {near}x downscale would be {W // near} or {-(-W // near)} x {H // near} or {-(-H // near)} px" if near >= 2 else ""
    raise SystemExit(f"[cells] {what} is {w} x {h} px, the image {W} x {H} px: not the image's full-resolution size nor a "
                     f"whole-number downscale of it{hint}{sug}. Segment the full-resolution image (or a whole-number "
                     f"downscale of it, as Cellpose's --rescale is not), or give the mask of this image.")


# ------------------------------------------------------------------------------------- label masks
class Mask:
    """A 2-D label image in a TIFF, read a band of rows at a time (pyramid._page_reader)."""

    def __init__(self, path: Path, what: str = "the label mask"):
        import tifffile
        import pyramid
        self.path, self.what = Path(path), what
        try:
            self._tf = tifffile.TiffFile(str(path))
        except Exception as e:                                         # not a TIFF at all
            raise SystemExit(f"[cells] {what} {self.path.name} cannot be read as a TIFF ({type(e).__name__}: {e})")
        s = self._tf.series[0]
        dims = [(n, a) for n, a in zip(s.shape, s.axes) if n > 1 or a in "YX"]
        if "S" in s.axes and s.shape[s.axes.index("S")] > 1:
            raise SystemExit(f"[cells] {what} {self.path.name} is a colour image ({s.shape}), not labels: give the label "
                             f"mask itself (one integer per cell, 0 = background)")
        if len(dims) != 2:
            raise SystemExit(f"[cells] {what} {self.path.name} has shape {tuple(s.shape)} (axes {s.axes}): one 2-D label "
                             f"image is needed. A file holding a whole-cell and a nuclear mask (Mesmer) must be split: the "
                             f"cells for --cells, the nuclei for --nuclei")
        self.dtype = np.dtype(s.dtype)
        if self.dtype == np.bool_ or self.dtype.kind not in "uif":
            raise SystemExit(f"[cells] {what} {self.path.name} holds {self.dtype} values: a label mask holds one integer per "
                             f"cell (0 = background), not yes / no")
        page = s.levels[0].keyframe if hasattr(s.levels[0], "keyframe") else s.pages[0]
        self.h, self.w = int(page.imagelength), int(page.imagewidth)
        self._get, self._close = pyramid._page_reader(page)

    def rows(self, y0: int, y1: int) -> np.ndarray:
        """Mask rows [y0, y1) as integers (labels as they are; a float mask's whole numbers as int64)."""
        b = np.asarray(self._get(int(y0), int(y1))).reshape(int(y1) - int(y0), self.w)
        if self.dtype.kind == "f":
            if not np.all(np.isfinite(b)) or np.any(b != np.round(b)):
                raise SystemExit(f"[cells] {self.what} {self.path.name}: floating-point values that are not whole numbers: "
                                 f"not a label mask")
            b = b.astype(np.int64)
        if self.dtype.kind == "i" and b.size and int(b.min()) < 0:
            raise SystemExit(f"[cells] {self.what} {self.path.name}: negative labels; 0 = background, cells 1, 2, 3 ...")
        return b

    def close(self):
        try:
            self._close()
        finally:
            self._tf.close()


def _mask_window(F: Frame, s: int, M: Mask):
    """Mask rows and columns [r0, r1) x [c0, c1) under the frame."""
    cols = F.src_cols() // s
    rows = F.src_rows(0, F.H) // s
    r0, r1 = int(max(rows.min(), 0)), int(min(rows.max() + 1, M.h))
    c0, c1 = int(max(cols.min(), 0)), int(min(cols.max() + 1, M.w))
    return r0, r1, c0, c1


class Labels:
    """The labels present in a mask window -> rows 0 .. n-1 (ascending label)."""

    def __init__(self, M: Mask, win, log=print):
        r0, r1, c0, c1 = win
        step = max(1, BAND_PX // max(c1 - c0, 1))
        present, uniq = np.zeros(1, bool), []          # a table of the labels seen, until one is above LUT_MAX
        for y in range(r0, r1, step):
            b = M.rows(y, min(y + step, r1))[:, c0:c1]
            if not b.size:
                continue
            mx = int(b.max())
            if present is not None and mx <= LUT_MAX:
                if mx >= len(present):
                    present = np.concatenate([present, np.zeros(mx + 1 - len(present), bool)])
                present[b.ravel()] = True
                continue
            if present is not None:                    # large label numbers: the sorted labels from here on
                uniq.append(np.flatnonzero(present))
                present = None
            uniq.append(np.unique(b).astype(np.int64))
        lab = np.flatnonzero(present).astype(np.int64) if present is not None else np.unique(np.concatenate(uniq))
        self.lut = None
        self.labels = lab[lab > 0]
        self.n = len(self.labels)
        if self.labels.size and self.labels[-1] <= LUT_MAX:
            self.lut = np.full(int(self.labels[-1]) + 1, -1, np.int32)
            self.lut[self.labels] = np.arange(self.n, dtype=np.int32)

    def index(self, b: np.ndarray) -> np.ndarray:
        """Row of each label in `b` (-1 = background or not a label of the window)."""
        if self.lut is not None:
            bi = b.astype(np.int64, copy=False)
            out = np.full(b.shape, -1, np.int32)
            ok = bi < len(self.lut)
            out[ok] = self.lut[bi[ok]]
            return out
        p = np.minimum(np.searchsorted(self.labels, b), max(self.n - 1, 0))
        return np.where((self.n > 0) & (self.labels[p] == b), p, -1).astype(np.int32)


def mask_labeller(F: Frame, s: int, M: Mask, L: Labels):
    """labels_of(y0, y1) for the measuring pass: the cell row of every frame pixel of rows [y0, y1) (-1 = none)."""
    mc = F.src_cols() // s
    okc = (mc >= 0) & (mc < M.w)
    mcc = np.clip(mc, 0, M.w - 1)

    def labels_of(y0: int, y1: int) -> np.ndarray:
        mr = F.src_rows(y0, y1) // s
        okr = (mr >= 0) & (mr < M.h)
        out = np.full((y1 - y0, F.W), -1, np.int32)
        if not okr.any():
            return out
        lo, hi = int(mr[okr].min()), int(mr[okr].max()) + 1
        blk = M.rows(lo, hi)
        sub = blk[np.clip(mr, lo, hi - 1) - lo][:, mcc]
        idx = L.index(sub)
        idx[~okr] = -1
        idx[:, ~okc] = -1
        return idx
    return labels_of


def trace_mask(M: Mask, s: int, F: Frame, win, row_of=None, log=print):
    """Outlines of the labels of a mask window, traced in windows (cellshapes.trace_labels) and put in the frame.

    `row_of(labels) -> rows` (-1 = leave out) picks the labels wanted; None keeps every label and returns its
    number as the id. Returns (rows or ids per ring, start, xy frame pixels, info)."""
    r0, r1, c0, c1 = win
    T, Hh = TRACE_TILE, TRACE_HALO
    tol = cellshapes.TRACE_TOL_PX * s
    keys, starts, xys, nv, holes = [], [], [], 0, 0
    traced, cut = set(), set()
    for wy in range(r0, r1, T):
        for wx in range(c0, c1, T):
            cy1, cx1 = min(wy + T, r1), min(wx + T, c1)
            ry0, ry1 = max(0, wy - Hh), min(M.h, cy1 + Hh)
            rx0, rx1 = max(0, wx - Hh), min(M.w, cx1 + Hh)
            lb, st, xy, h = cellshapes.trace_labels(M.rows(ry0, ry1)[:, rx0:rx1])
            if not len(lb):
                continue
            holes += h
            # a label is whole in this window when none of its vertices lies on an inner edge of the window
            ww, wh = rx1 - rx0, ry1 - ry0
            edge = (((xy[:, 0] == 0) & (rx0 > 0)) | ((xy[:, 0] == ww) & (rx1 < M.w))
                    | ((xy[:, 1] == 0) & (ry0 > 0)) | ((xy[:, 1] == wh) & (ry1 < M.h)))
            ring_cut = np.add.reduceat(edge.astype(np.int64), st[:-1]) > 0
            # its anchor: the top-left corner of its first pixel in raster order (the least (y, x) of its vertices)
            ymin = np.minimum.reduceat(xy[:, 1], st[:-1])
            ring = np.repeat(np.arange(len(lb)), np.diff(st))
            xmin = np.minimum.reduceat(np.where(xy[:, 1] == ymin[ring], xy[:, 0], np.inf), st[:-1])
            u, inv = np.unique(lb, return_inverse=True)
            lab_cut = np.bincount(inv, weights=ring_cut, minlength=len(u)) > 0
            order = np.lexsort((xmin, ymin, inv))
            first = order[np.concatenate([[0], np.flatnonzero(np.diff(inv[order])) + 1])]
            ax, ay = xmin[first] + rx0, ymin[first] + ry0
            own = (ay >= wy) & (ay < cy1) & (ax >= wx) & (ax < cx1)
            cut.update(int(v) for v in u[own & lab_cut])
            take = (own & ~lab_cut)[inv]
            if row_of is not None:
                rws = row_of(lb)
                take &= rws >= 0
            if not take.any():
                continue
            traced.update(int(v) for v in np.unique(lb[take]))
            st2, xy2 = cellshapes._keep_rings(st, xy, take)
            pts = F.to_frame((xy2 + [rx0, ry0]) * s)
            st3, pts, _ = cellshapes.simplify(st2, pts, tol)
            keys.append(rws[take] if row_of is not None else lb[take])
            starts.append(st3[1:] + nv)
            xys.append(pts)
            nv += len(pts)
    big = len(cut - traced)
    if big:
        log(f"[cells] {big:,} labels of {M.path.name} reach more than {TRACE_HALO} mask pixels beyond a {TRACE_TILE}-pixel "
            f"tracing window: no outline (they stay dots)")
    start = np.concatenate([[0]] + starts) if starts else np.zeros(1, np.int64)
    key = np.concatenate(keys) if keys else np.zeros(0, np.int64)
    info = {"holes": holes, "traced": True, "traceTolerancePx": cellshapes.TRACE_TOL_PX, "maskScale": s,
            "tooLargeToTrace": big}
    return key, start, (np.concatenate(xys) if xys else np.zeros((0, 2))), info


# ------------------------------------------------------------------------------------- GeoJSON polygons
def has_nucleus_geometry(path: Path) -> bool:
    """Does the file mention QuPath's nucleusGeometry (a byte scan, nothing parsed)?"""
    tail, pat = b"", b'"nucleusGeometry"'
    with open(path, "rb") as f:
        while True:
            buf = f.read(1 << 22)
            if not buf:
                return False
            if pat in tail + buf:
                return True
            tail = buf[-len(pat):]


def read_polygons(path: Path, units: str, um: float, F: Frame, geometry: str = "geometry"):
    """GeoJSON polygons -> (ids per ring, start, xy in frame pixels, info)."""
    ids, start, xy, info = cellshapes.read_geojson(path, skip_types=SKIP_TYPES, number_missing=True, geometry=geometry)
    if units == "um":
        xy = xy / float(um)
    return ids, start, F.to_frame(xy) if len(xy) else xy, info


class Polygons:
    """Cells as polygons (frame pixels), rasterised a band of rows at a time: a pixel is in a ring when its centre
    is (even-odd, centres at x + 0.5, y + 0.5; a centre on a left or top edge is in)."""

    def __init__(self, start, xy, cell):
        start, xy = cellshapes.clean_rings(np.asarray(start, np.int64), np.asarray(xy, np.float64))
        cnt = np.diff(start)
        ok = cnt >= 3
        cell = np.asarray(cell, np.int64)[ok]
        start, xy = cellshapes._keep_rings(start, xy, ok)
        ring = np.repeat(np.arange(len(start) - 1), np.diff(start))
        nxt = cellshapes._next(start) if len(xy) else np.zeros(0, np.int64)
        x0, y0, x1, y1 = xy[:, 0], xy[:, 1], xy[nxt, 0], xy[nxt, 1]
        keep = y0 != y1                                                  # horizontal edges cross no row centre
        x0, y0, x1, y1, ring = x0[keep], y0[keep], x1[keep], y1[keep], ring[keep]
        lo, hi = np.minimum(y0, y1), np.maximum(y0, y1)
        o = np.argsort(lo, kind="stable")
        self.x0, self.y0, self.x1, self.y1 = x0[o], y0[o], x1[o], y1[o]
        self.lo, self.hi, self.ring = lo[o], hi[o], ring[o]
        self.cell = cell                                                  # cell row of each ring
        self.maxh = float((hi - lo).max()) if len(lo) else 0.0
        self.start, self.xy = start, xy
        self.overlap = 0                                                  # pixels inside two cells: left out

    def labels_of(self, y0: int, y1: int, W: int) -> np.ndarray:
        out = np.full((y1 - y0, W), -1, np.int32)
        a = np.searchsorted(self.lo, y0 - self.maxh - 1.0, "left")
        b = np.searchsorted(self.lo, y1, "right")
        sel = np.arange(a, b)
        sel = sel[self.hi[sel] > y0]
        if not len(sel):
            return out
        lo, hi = self.lo[sel], self.hi[sel]
        ja = np.maximum(np.ceil(lo - 0.5).astype(np.int64), y0)
        jb = np.minimum(np.ceil(hi - 0.5).astype(np.int64) - 1, y1 - 1)
        n = np.maximum(jb - ja + 1, 0)
        if not n.sum():
            return out
        e = np.repeat(sel, n)
        j = np.repeat(ja, n) + (np.arange(int(n.sum())) - np.repeat(np.cumsum(n) - n, n))
        yc = j + 0.5
        x = self.x0[e] + (yc - self.y0[e]) * (self.x1[e] - self.x0[e]) / (self.y1[e] - self.y0[e])
        rg = self.ring[e]
        o = np.lexsort((x, j, rg))
        x, j, rg = x[o], j[o], rg[o]
        # crossings pair up within each (ring, row): a closed ring crosses a row centre an even number of times
        grp = np.concatenate([[True], (np.diff(rg) != 0) | (np.diff(j) != 0)])
        gid = np.cumsum(grp) - 1
        odd = np.bincount(gid) % 2 == 1
        if odd.any():
            keep = ~odd[gid]
            x, j, rg = x[keep], j[keep], rg[keep]
        xa, xb, row, rr = x[0::2], x[1::2], j[0::2], rg[0::2]
        ca = np.clip(np.ceil(xa - 0.5).astype(np.int64), 0, W)
        cb = np.clip(np.ceil(xb - 0.5).astype(np.int64), 0, W)
        ok = ca < cb
        ca, cb, row, cid = ca[ok], cb[ok], row[ok] - y0, self.cell[rr[ok]].astype(np.float64) + 1
        R, W1 = y1 - y0, W + 1
        size = R * W1
        ia, ib = row * W1 + ca, row * W1 + cb

        def cum(w):
            d = np.bincount(ia, weights=w, minlength=size) - np.bincount(ib, weights=w, minlength=size)
            return np.cumsum(d.reshape(R, W1), axis=1)[:, :W]
        one = np.ones(len(ia))
        cnt, sm, sq = cum(one), cum(cid), cum(cid * cid)
        cnt = np.rint(cnt)
        single = (cnt == 1) | ((cnt > 1) & (np.abs(sm * sm - cnt * sq) < 0.5))
        lab = np.where(single, np.rint(sm / np.maximum(cnt, 1)), 0).astype(np.int64)
        self.overlap += int(((cnt > 0) & ~single).sum())
        out[:] = (lab - 1).astype(np.int32)
        return out


# ------------------------------------------------------------------------------------- measuring
def measure(r, labels_of, n: int, log=print, workers: int = 4) -> dict:
    """Area, centroid and per-channel sums of every cell row 0 .. n-1 from the cache's level 0 (8-bit units at the
    end). `labels_of(y0, y1)` -> (y1 - y0, W) int32 cell rows (-1 none)."""
    import protgrid
    a = r._arrays[0]
    C, H, W = len(r.channels), int(r.H), int(r.W)
    rows = int(a.chunks[1])
    step = max(1, BAND_PX // max(W, 1))
    A, SX, SY = np.zeros(n), np.zeros(n), np.zeros(n)
    S = np.zeros((C, n))
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=workers) as ex:
        for y0 in range(0, H, rows):
            y1 = min(y0 + rows, H)
            subs = []
            for s0 in range(y0, y1, step):
                s1 = min(s0 + step, y1)
                flat = labels_of(s0, s1).ravel()
                m = flat >= 0
                if not m.any():
                    continue
                fi = flat[m]
                pos = np.flatnonzero(m)
                A += np.bincount(fi, minlength=n)
                SX += np.bincount(fi, weights=pos % W + 0.5, minlength=n)
                SY += np.bincount(fi, weights=pos // W + (s0 + 0.5), minlength=n)
                subs.append((s0 - y0, s1 - y0, m, fi))
            if not subs:
                continue

            def chan(c, y0=y0, y1=y1, subs=subs):
                band = np.asarray(a[c, y0:y1, :W])
                out = np.zeros(n)
                for r0_, r1_, m, fi in subs:
                    out += np.bincount(fi, weights=band[r0_:r1_].ravel()[m], minlength=n)
                return out
            for c, s_ in enumerate(ex.map(chan, range(C))):
                S[c] += s_
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = S / np.where(A > 0, A, 1.0)
    V = np.empty((n, C), np.float32)
    for c in range(C):
        V[:, c] = protgrid.to_8bit(mean[c], r.channels[c], a.dtype)
    with np.errstate(invalid="ignore", divide="ignore"):
        cx, cy = SX / A, SY / A
    log(f"[cells] measured {n:,} cells x {C} channels on full-resolution pixels ({W} x {H})  {time.time() - t0:.0f}s")
    return {"area": A, "x": cx, "y": cy, "V": V}


# ------------------------------------------------------------------------------------- prepare
def precheck(args, manifest: dict, full_size, log=print) -> dict | None:
    """Before the pyramid is built: the files exist, are of a kind read here, and a mask's size fits the image.
    `full_size` = (width, height) of the image file as segmented (before --crop)."""
    cells, nuclei = getattr(args, "cells", None), getattr(args, "nuclei", None)
    if not cells and not nuclei:
        return None
    if nuclei and not cells:
        raise SystemExit("[cells] --nuclei needs --cells (the nucleus outlines belong to the cells of --cells)")
    if not manifest.get("pixelSizeUm"):
        raise SystemExit("[cells] --cells needs the image's pixel size (µm per pixel): give --pixel-size")
    plan = {"units": getattr(args, "cells_units", "px") or "px"}
    W, H = (int(v) for v in full_size)
    for k, f in (("cells", cells), ("nuclei", nuclei)):
        if not f:
            continue
        p = Path(f)
        if not p.is_file():
            raise SystemExit(f"[cells] --{k} {p}: not found")
        low = p.name.lower()
        if low.endswith(MASK_EXT):
            M = Mask(p, f"the --{k} mask")
            try:
                s = mask_scale((M.h, M.w), (H, W), f"the --{k} mask {p.name}")
            finally:
                M.close()
            plan[k] = {"path": p, "type": "mask", "scale": s}
            log(f"[cells] --{k}: label mask {M.w} x {M.h} px" + (f", {s}x downscale of the image" if s > 1 else
                                                                  ", the image's full-resolution size"))
        elif low.endswith(GEOJSON_EXT):
            plan[k] = {"path": p, "type": "geojson"}
            log(f"[cells] --{k}: GeoJSON polygons in {'µm' if plan['units'] == 'um' else 'full-resolution pixels'}")
        else:
            raise SystemExit(f"[cells] --{k} {p.name}: expected a label mask TIFF (.tif / .tiff) or a GeoJSON "
                             f"(.geojson / .json). Cellpose: save the masks as TIFF (--save_tif); QuPath: "
                             f"Export objects as GeoJSON, or a labelled image as TIFF")
    if plan["units"] == "um" and not any(plan[k]["type"] == "geojson" for k in ("cells", "nuclei") if k in plan):
        log("[cells] --cells-units applies to GeoJSON only: a mask has its own pixel grid")
    manifest.setdefault("source", {})["cells"] = str(cells)
    if nuclei:
        manifest["source"]["nuclei"] = str(nuclei)
    return plan


def _quartiles(a):
    return [round(float(v), 2) for v in np.percentile(a, [25, 50, 75])] if len(a) else []


def build(plan: dict, ds_dir: Path, manifest: dict, args=None, log=print) -> dict | None:
    """After the pyramid: measure the cells, write expr/cells.zarr and the outlines, add the manifest entry."""
    import expr as expr_mod
    from pyramid import PyramidReader
    from scipy import sparse
    t0 = time.time()
    ds_dir = Path(ds_dir)
    um = float(manifest["pixelSizeUm"])
    F = Frame(manifest)
    r = PyramidReader(ds_dir, manifest)
    cp = plan["cells"]
    seg = {"from": "label mask" if cp["type"] == "mask" else "GeoJSON polygons",
           "measured": "mean of every channel over the cell's full-resolution pixels (8-bit units, as the grid squares)",
           "segmentedBy": "the lab (SODA does not segment)"}
    sets = []
    if cp["type"] == "mask":
        s = cp["scale"]
        M = Mask(cp["path"], "the --cells mask")
        try:
            win = _mask_window(F, s, M)
            L = Labels(M, win, log=log)
            if L.n == 0:
                raise SystemExit(f"[cells] the --cells mask {cp['path'].name} holds no cell under the image (every pixel 0)")
            if L.n == 1:
                raise SystemExit(f"[cells] the --cells mask {cp['path'].name} holds one label only: a yes / no (tissue or "
                                 f"foreground) mask, not one number per cell")
            log(f"[cells] {L.n:,} labels in the mask (largest {int(L.labels[-1]):,})")
            meas = measure(r, mask_labeller(F, s, M, L), L.n, log=log)
            keep = meas["area"] > 0
            labels = L.labels[keep]
            rows_of_label = np.full(L.n, -1, np.int64)
            rows_of_label[keep] = np.arange(int(keep.sum()))

            def row_of(lb):
                return rows_of_label[L.index(lb)]
            key, start, xy, info = trace_mask(M, s, F, win, row_of=row_of, log=log)
            if len(key):
                sets.append(cellshapes.ShapeSet("cell", start, xy, rows=key, source="label mask (own file, traced)",
                                                tol_px_min=cellshapes.TRACE_TOL_PX * s, info=info))
        finally:
            M.close()
        barcodes = [str(int(v)) for v in labels]
        seg["maskScale"] = s
    else:
        ids, start, xy, info = read_polygons(cp["path"], plan["units"], um, F)
        if not len(ids):
            raise SystemExit(f"[cells] {cp['path'].name}: no Polygon / MultiPolygon features"
                             + (f" ({info.get('featuresSkipped', 0):,} skipped by objectType)" if info.get("featuresSkipped") else ""))
        uid, cell = np.unique(np.asarray(ids, dtype=object).astype(str), return_inverse=True)
        first = np.full(len(uid), len(ids), np.int64)
        np.minimum.at(first, cell, np.arange(len(ids)))
        o = np.argsort(first, kind="stable")                            # cells in the file's order
        rank = np.empty(len(uid), np.int64)
        rank[o] = np.arange(len(uid))
        cell = rank[cell]
        uid = uid[o]
        P = Polygons(start, xy, cell)
        lo, hi = xy.min(0) if len(xy) else (0, 0), xy.max(0) if len(xy) else (0, 0)
        meas = measure(r, lambda y0, y1: P.labels_of(y0, y1, F.W), len(uid), log=log)
        keep = meas["area"] > 0
        if keep.sum() < 0.5 * len(uid):
            raise SystemExit(f"[cells] only {int(keep.sum()):,} of {len(uid):,} cells of {cp['path'].name} cover a pixel of "
                             f"the image (the polygons span x {lo[0]:.0f} .. {hi[0]:.0f}, y {lo[1]:.0f} .. {hi[1]:.0f}; the "
                             f"image is {F.W} x {F.H} px): wrong units (--cells-units "
                             f"{'px' if plan['units'] == 'um' else 'um'}?), another image, or another "
                             f"crop. Nothing written.")
        if not keep.all():
            log(f"[cells] WARNING {int((~keep).sum()):,} of {len(uid):,} cells cover no pixel of the image and are left out")
        if P.overlap:
            log(f"[cells] {P.overlap:,} pixels lie inside the polygons of two cells: left out of both")
        rows_of_cell = np.full(len(uid), -1, np.int64)
        rows_of_cell[keep] = np.arange(int(keep.sum()))
        ring_rows = rows_of_cell[P.cell]
        ok = ring_rows >= 0
        if ok.any():
            st2, xy2 = cellshapes._keep_rings(P.start, P.xy, ok)
            sets.append(cellshapes.ShapeSet("cell", st2, xy2, rows=ring_rows[ok], source="GeoJSON polygons (own file)",
                                            info={"holes": info.get("holes", 0)}))
        barcodes = [str(v) for v in uid[keep]]
        seg.update({"featuresSkipped": info.get("featuresSkipped", 0), "featuresNumbered": info.get("featuresNumbered", 0),
                    "pixelsInTwoCells": P.overlap})
    n = len(barcodes)
    area = meas["area"][keep]
    xy_c = np.column_stack([meas["x"][keep], meas["y"][keep]]).astype(np.float32)
    V = meas["V"][keep]
    d_px = 2.0 * math.sqrt(float(np.median(area)) / math.pi)
    a_um = area * um * um
    seg.update({"cells": n, "areaUm2": _quartiles(a_um)})
    log(f"[cells] {n:,} cells: area {seg['areaUm2'][0]:g} / {seg['areaUm2'][1]:g} / {seg['areaUm2'][2]:g} µm² "
        f"(quartiles), diameter about {d_px * um:.1f} µm")
    if not 3.0 <= d_px * um <= 60.0:
        log(f"[cells] WARNING a typical cell is {d_px * um:.1f} µm across: check the pixel size and the mask's scale")
    names = [str(ch["name"]) for ch in r.channels]
    meta = {"key": KEY, "label": LABEL, "kind": "cells", "modality": "protein", "nBins": int(n), "nGenes": len(names),
            "binSizeUm": None, "cellDiameterUm": round(d_px * um, 2), "spotDiameterPx": d_px, "micronsPerPixel": um,
            "imageScale": 1.0, "segmentation": seg}
    csr = sparse.csr_matrix(V.T.astype(np.float32))
    csr.sort_indices()
    expr_mod._write_layer(ds_dir, xy_c, csr, names, meta, barcodes)
    import zarr
    g = zarr.open_group(str(ds_dir / "expr" / f"{KEY}.zarr"), mode="a")
    g.create_dataset("area", data=area.astype(np.float32), overwrite=True)
    # nuclei: by id (the cells' mask labels / feature ids), or QuPath's nucleusGeometry of the cells' own file
    npl = plan.get("nuclei")
    where = {b: i for i, b in enumerate(barcodes)}
    nuc = None
    if npl and npl["type"] == "mask":
        Mn = Mask(npl["path"], "the --nuclei mask")
        try:
            ids_n, st_n, xy_n, info_n = trace_mask(Mn, npl["scale"], F, _mask_window(F, npl["scale"], Mn), log=log)
        finally:
            Mn.close()
        nuc = ([str(int(v)) for v in ids_n], st_n, xy_n, info_n, "label mask (own file, traced)",
               cellshapes.TRACE_TOL_PX * npl["scale"])
    elif npl or (cp["type"] == "geojson" and has_nucleus_geometry(cp["path"])):
        src, geo = (npl["path"], "geometry") if npl else (cp["path"], "nucleusGeometry")
        ids_n, st_n, xy_n, info_n = read_polygons(src, plan["units"], um, F, geometry=geo)
        nuc = (list(ids_n), st_n, xy_n, info_n, "GeoJSON polygons (own file)" if npl else "QuPath nucleusGeometry (own file)", 0.0)
        if not npl:
            log(f"[cells] nucleus outlines from the nucleusGeometry of {cp['path'].name}")
    if nuc is not None:
        ids_n, st_n, xy_n, info_n, what, tmin = nuc
        rws = np.fromiter((where.get(i, -1) for i in ids_n), np.int64, len(ids_n))
        hit = rws >= 0
        if hit.sum() < len(rws):
            log(f"[cells] {int((~hit).sum()):,} of {len(rws):,} nucleus outlines have no cell of --cells with their id: left out")
        if hit.any():
            st_k, xy_k = cellshapes._keep_rings(np.asarray(st_n, np.int64), np.asarray(xy_n, np.float64), hit)
            sets.append(cellshapes.ShapeSet("nucleus", st_k, xy_k, rows=rws[hit], source=what, tol_px_min=tmin,
                                            info={k: v for k, v in info_n.items() if k != "tooLargeToTrace"}))
        else:
            log("[cells] no nucleus outline belongs to a cell of --cells (the ids differ): none written")
    cellshapes.attach(ds_dir, meta, sets, um, args=args, log=log)
    if not meta.get("shapes"):
        log("[cells] no outlines written: the cells are drawn as dots")
    manifest["expr"] = [L_ for L_ in manifest.get("expr", []) if L_.get("key") != KEY] + [meta]
    log(f"[cells] layer '{KEY}': {n:,} cells x {len(names)} channels  ({time.time() - t0:.0f}s)")
    return meta


# ------------------------------------------------------------------------------------- export
def is_layer(layer: dict) -> bool:
    return layer.get("kind") == "cells" and layer.get("modality") == "protein"


def export_layer(er, d: Path, layer: dict, px_um: float, grid_bg: dict | None = None, markers: Path | None = None,
                 umap: bool = True, celltypes: bool = True, umap_sample: int | None = None, force: bool = False,
                 log=print) -> tuple[dict, dict | None]:
    """expr/cells/: gene files, mat/, umap/ (+ suggested names) of a protein cells layer, by the grid's rules.
    `grid_bg` {channel name: background} from the grid layer of the same export (protgrid.write_grid_layer's
    `stats`), None without one. Returns (fields for the manifest entry, the umap entry or None)."""
    import celltypes as celltypes_mod
    import genepack
    import protgrid as pg
    import sectionumap
    import deid
    t0 = time.time()
    d = Path(d)
    key = layer["key"]
    pos = er.positions()
    n = len(pos)
    names = list(er.genes)
    if len(er.xy) != n:
        raise SystemExit(f"[cells] {key}: {len(er.xy) - n} cells without a position (the cache is not protcells.py's)")
    V = np.zeros((n, len(names)), np.float64)
    indices, data, indptr = np.asarray(er.g["indices"]), np.asarray(er.g["data"], np.float64), er.indptr
    for gi in range(len(names)):
        a, b = int(indptr[gi]), int(indptr[gi + 1])
        V[indices[a:b], gi] = data[a:b]
    written = [j for j, nm in enumerate(names) if not pg.HIDDEN.match(nm)]
    dropped = [{"name": names[j], "reason": "open / blank channel: not written"} for j in range(len(names)) if j not in written]
    if not written or n < 2:
        log(f"[cells] {key}: nothing to write ({n} cells, {len(written)} channels)")
        return {}, None
    if grid_bg and all(names[j] in grid_bg for j in written):
        bg_all = np.array([grid_bg.get(nm, 0.0) for nm in names], np.float64)
        bg_rule = "the grid layer's background of the channel (median of its off-tissue squares)"
    else:
        bg_all, _ = pg.background(V, np.ones(n, bool))
        bg_rule = "1st percentile of the cells (no grid layer in this export to give the slide's background)"
    nonnuc = [j for j in written if not pg.NUCLEAR.match(names[j])]
    bright = pg.bright_squares(V[:, nonnuc])
    ref = ~bright if (~bright).sum() >= max(50, n // 2) else np.ones(n, bool)
    hi_all = np.percentile(np.maximum(V[ref] - bg_all, 0.0), pg.P_HI, axis=0)
    snr_all = hi_all / np.maximum(bg_all, 0.5)
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
    Vt = V[:, written]
    wnames = [names[j] for j in written]
    bg, hi = bg_all[written], hi_all[written]
    flat = hi < pg.FLAT_RANGE
    embed = []
    for k, j in enumerate(written):
        if flat[k]:
            dropped.append({"name": names[j], "reason": f"flat over the cells (p99.5 - background = {hi[k]:.2f}): "
                                                         "written, not embedded"})
        elif pg.NUCLEAR.match(names[j]):
            dropped.append({"name": names[j], "reason": "nuclear stain: written, not embedded"})
        else:
            embed.append(k)
    x = pg.normalise_bg_asinh(Vt, bg, hi)
    log(f"[cells] {key}: {len(wnames)} channels written, {len(embed)} embedded; background: {bg_rule}")
    # gene files: every cell, in positions.bin order (as the grid's squares)
    xyv = np.empty((n, 3), np.float32)
    xyv[:, :2] = pos
    pw = genepack.PackWriter(d, unit=12, log=log)
    where = {}
    for k, nm in enumerate(wnames):
        xyv[:, 2] = Vt[:, k]
        buf = np.ascontiguousarray(xyv, dtype="<f4")
        if pw is not None and nm not in where:
            where[nm] = pw.add(nm, buf)
    pack = pw.finish() if pw is not None else None
    if pack:
        genepack.note_loose(d, log=log)
    genes = [{"n": nm, "t": round(float(Vt[:, k].mean()), 4)} for k, nm in enumerate(wnames)]
    deid.assert_clean(genes, "cells genes.json")
    genepack.atomic_json(d / "genes.json", genepack.with_index(genes, where))
    flags = np.zeros(n, np.uint8)
    flags[bright] |= pg.FLAG_BRIGHT
    mat = pg.write_mat_chunks(d, pos, Vt[:, embed], [wnames[k] for k in embed], bg[embed], hi[embed], float(px_um),
                              log=log) if embed else None
    entry = {"nGenes": len(wnames), "channels": [wnames[k] for k in embed], "normalize": pg.NORMALIZE,
             "geneSummary": "mean intensity per cell (8-bit units, full-resolution pixels)",
             "background": bg_rule, "dropped": dropped,
             "lowSignalCells": int((~(x[:, embed] >= pg.LOW_SIGNAL_X).any(1)).sum()) if embed else n,
             "everythingBrightCells": int(bright.sum())}
    if mat:
        entry["mat"] = mat
    if pack:
        entry["genePack"] = pack
    u = None
    if umap and embed:
        feats = pg.clr_features(Vt[:, embed], bg[embed], hi[embed]) if len(embed) >= pg.CLR_MIN_CHANNELS else None
        keep = np.flatnonzero((x[:, embed] >= pg.LOW_SIGNAL_X).any(1))
        log(f"[cells] {n - len(keep):,} of {n:,} cells low signal (no embedded marker at x >= {pg.LOW_SIGNAL_X:g}): "
            "not embedded or clustered")
        hh = hashlib.sha256(np.ascontiguousarray(feats).tobytes() if feats is not None else
                            np.ascontiguousarray(x[:, embed], dtype=np.float32).tobytes())
        hh.update(np.ascontiguousarray(keep.astype("<i8")).tobytes())
        u = sectionumap.write_feature_umap(d, x, wnames, embed, embed, key, log=log, force=force, sample=umap_sample,
                                           features=feats, keep=keep,
                                           extra_params={"unit": "cell (the lab's own segmentation)", "normalize": pg.NORMALIZE,
                                                         "cofactor": pg.COFACTOR, "background": bg_rule,
                                                         "embedFeatures": (f"centred log-ratio of max(v - bg, 0) / hi + {pg.CLR_EPS:g}"
                                                                           if feats is not None else
                                                                           f"x (fewer than {pg.CLR_MIN_CHANNELS} channels)"),
                                                         "lowSignal": f"no embedded marker at x >= {pg.LOW_SIGNAL_X:g}: not embedded",
                                                         "input": hh.hexdigest()[:16]})
        if u and celltypes:
            c = celltypes_mod.write_celltypes_protein(d, x, wnames, ~flat, markers=markers, key=key, log=log, flags=flags)
            if c:
                u = {**u, "celltypes": c}
    deid.assert_clean(entry, "cells manifest entry")
    log(f"[cells] {key} done  {time.time() - t0:.0f}s")
    return entry, u
