"""Cell outlines (segmentation boundaries) -> spatial tiles of quantised polygons.

Segmented platforms (Xenium, MERSCOPE, CosMx, Stereo-seq cell bins) come with the outline of
every cell, and often of every nucleus. This module LOADS what the instrument or its software
wrote (it segments nothing) and writes it so the viewer can fetch the outlines of the part of
the section in view:

    shapes/<layer>/<set>/index.json      in the prepare cache          set = "cell" | "nucleus"
    shapes/<layer>/<set>/<tx>_<ty>.bin
    expr/<layer>/shapes/<set>/...        the same files in the export (copied by export_shapes)

    <tx>_<ty>.bin, little-endian:
        u32 n                 polygons in this tile
        u32 nv                vertices in this tile
        f32 x0, y0            origin of the tile's vertices, full-resolution image pixels
        u32 row[n]            row of positions.bin (never a cell id); a row may repeat (a cell in parts)
        u16 count[n]          vertices of each polygon; 0 = the cell has no outline (it stays a dot)
        u16 qx, qy [nv]       vertex = (x0 + qx * step, y0 + qy * step), step in index.json

A polygon is in the tile of its cell's own position, floor(x / chunkPx), floor(y / chunkPx) on
the float32 values of positions.bin with chunkPx = round(500 um / um per pixel): the grid hvg/
uses. Polygons are never cut at a tile edge; index.json "overhangPx" says how far a vertex can
reach outside its tile. Rings are open (the first vertex is not repeated), clockwise on screen
(x right, y down), outer rings only. In the first set every row appears at least once.

Readers (each returns rings in the source's own coordinates; the platform glue below maps them
to image pixels with the same transform as the cell centroids):

    read_vertex_table     Xenium cell_boundaries / nucleus_boundaries .csv.gz | .parquet, a lab's own table
    read_wkb_table        MERSCOPE cell_boundaries.parquet / *_micron_space.parquet (WKB), GeoParquet
    read_geojson          a lab's own polygons
    read_xenium_zarr      Xenium cells.zarr.zip (polygon_sets, or polygon_vertices of older runs)
    read_merscope_hdf5    MERSCOPE cell_boundaries/feature_data_<fov>.hdf5 (older runs)
    read_cell_border      Stereo-seq cell bin GEF, cellBin/cellBorder
    trace_labels          CosMx CellLabels images: the outline of every label, exactly along its pixels

Needs numpy and what requirements.txt already lists (pyarrow, h5py, zarr, tifffile); WKB is
parsed here, label masks are traced here.
"""
from __future__ import annotations

import json
import re
import shutil
import struct
import time
from pathlib import Path

import numpy as np

import deid

DIR = "shapes"
FORMAT = "poly-u16-tiles"
VERSION = 1
CHUNK_UM = 500.0                 # the hvg/ chunk size, so both use one grid
TOLERANCE_UM = 0.1               # Douglas-Peucker tolerance; 0 keeps every vertex
MAX_VERTICES = 64                # a polygon over this is simplified further
STEP_MIN = 1.0 / 16              # finest vertex step, image pixels
Q = 65535
MIN_PART = 0.02                  # a part under this share of its cell's area is dropped
TRACE_TOL_PX = 1.5               # outlines traced from a label image: tolerance in mask pixels (removes the staircase)
MATCH_MIN = 0.5                  # under this share of polygons with a row, or of centres on their outline: not written
MATCH_WARN = 0.9
LABELS = {"cell": "Cell boundaries", "nucleus": "Nucleus boundaries"}


class ShapeSet:
    """One set of outlines for a layer: `start` (n + 1 offsets) into `xy` (vertices, image pixels),
    and for each polygon either the cache row of its cell (`rows`, -1 = none) or its cell id
    (`ids`, matched to the layer's barcodes). `tol_px_min`: a floor for the simplification
    tolerance in image pixels (traced masks: TRACE_TOL_PX mask pixels). `info`: reader notes for index.json."""

    def __init__(self, key, start, xy, rows=None, ids=None, source="", tol_px_min=0.0, info=None, label=None):
        self.key, self.label = key, label or LABELS.get(key, key)
        self.start = np.asarray(start, np.int64)
        self.xy = np.asarray(xy, np.float64).reshape(-1, 2)
        self.rows = None if rows is None else np.asarray(rows, np.int64)
        self.ids = ids
        self.source, self.tol_px_min, self.info = source, float(tol_px_min), dict(info or {})
        if len(self.start) < 1 or self.start[-1] != len(self.xy):
            raise ValueError(f"{key}: offsets do not cover the vertices")


# ------------------------------------------------------------------------- ring arithmetic
def _ring_of(start):
    return np.repeat(np.arange(len(start) - 1), np.diff(start))


def _sum(values, start):
    """Sum of `values` over each ring (empty rings give 0)."""
    c = np.concatenate([[0], np.cumsum(values)])
    return c[start[1:]] - c[start[:-1]]


def _next(start):
    """Index of each vertex's successor in its ring (rings must not be empty)."""
    nxt = np.arange(1, start[-1] + 1)
    nxt[start[1:] - 1] = start[:-1]
    return nxt


def _keep_vertices(start, xy, keep):
    n = _sum(keep.astype(np.int64), start)
    return np.concatenate([[0], np.cumsum(n)]), xy[keep]


def _keep_rings(start, xy, keep):
    cnt = np.diff(start)
    return np.concatenate([[0], np.cumsum(cnt[keep])]), xy[np.repeat(keep, cnt)]


def clean_rings(start, xy):
    """Drop vertices that are not finite or repeat the one before them in the ring (so also the
    closing vertex of a closed ring)."""
    cnt = np.diff(start)
    if not len(xy):
        return start, xy
    ring = _ring_of(start)
    prev = np.arange(-1, len(xy) - 1)
    prev[start[:-1][cnt > 0]] = start[1:][cnt > 0] - 1
    keep = np.isfinite(xy).all(1) & ~((xy == xy[prev]).all(1) & (cnt[ring] > 1))
    return _keep_vertices(start, xy, keep)


def signed_area(start, xy):
    """Shoelace area of each ring; positive = clockwise on screen (x right, y down)."""
    if not len(xy):
        return np.zeros(len(start) - 1)
    cnt = np.diff(start)
    first = np.repeat(start[:-1], cnt)
    x, y = xy[:, 0] - xy[first, 0], xy[:, 1] - xy[first, 1]
    nxt = _next(start)
    return 0.5 * _sum(x * y[nxt] - x[nxt] * y, start)


def _seg_dist(px, py, ax, ay, bx, by):
    """Distance from points to segments, element by element."""
    dx, dy = bx - ax, by - ay
    L2 = dx * dx + dy * dy
    t = np.clip(((px - ax) * dx + (py - ay) * dy) / np.where(L2 > 0, L2, 1.0), 0.0, 1.0)
    return np.hypot(px - (ax + t * dx), py - (ay + t * dy))


def _dp(P, keep, tol):
    """Douglas-Peucker on many open chains at once. `keep` marks the points every chain starts
    with (its two ends at least); a point farther than `tol` from the segment between its kept
    neighbours is kept, the farthest of each segment first, until none is left."""
    x, y = P[:, 0], P[:, 1]
    while True:
        kidx = np.flatnonzero(keep)
        seg = np.cumsum(keep) - 1
        a = kidx[seg]
        b = kidx[np.minimum(seg + 1, len(kidx) - 1)]
        d = _seg_dist(x, y, x[a], y[a], x[b], y[b])
        d[keep] = 0.0
        far = d > tol
        if not far.any():
            return keep
        mx = np.maximum.reduceat(d, kidx)
        cand = np.flatnonzero(far & (d >= mx[seg]))
        sc = seg[cand]
        keep[cand[np.r_[True, sc[1:] != sc[:-1]]]] = True


def _simplify_once(start, xy, idx, tol):
    """One Douglas-Peucker pass over the rings `idx`: (their vertex indices, which of them stay)."""
    c = np.diff(start)[idx]
    s = start[:-1][idx]
    cs = np.concatenate([[0], np.cumsum(c + 1)])              # chains: a ring's vertices + its first one again
    ring = np.repeat(np.arange(len(idx)), c + 1)
    local = np.arange(cs[-1]) - cs[:-1][ring]
    closing = local == c[ring]
    src = s[ring] + np.where(closing, 0, local)
    P = xy[src]
    # cut each ring at its first vertex and at the vertex farthest from it
    d0 = np.hypot(P[:, 0] - P[cs[:-1], 0][ring], P[:, 1] - P[cs[:-1], 1][ring])
    far = np.flatnonzero(d0 >= np.maximum.reduceat(d0, cs[:-1])[ring])
    far = far[np.r_[True, ring[far][1:] != ring[far][:-1]]]
    keep = np.zeros(len(P), bool)
    keep[cs[:-1]] = keep[cs[1:] - 1] = keep[far] = True
    keep = _dp(P, keep, tol)
    return src[~closing], keep[~closing]


def simplify(start, xy, tol, max_vertices=MAX_VERTICES):
    """Douglas-Peucker on closed rings, all at once: no vertex that is dropped lies farther than
    `tol` (image pixels) from the simplified outline. A ring left with more than `max_vertices`
    is done again from all its vertices with 1.5 x the tolerance, until it fits; a ring that
    would be left with fewer than 3 vertices (a cell smaller than the tolerance) is not
    simplified. tol = 0 only enforces `max_vertices`. Returns (start, xy, number of rings whose tolerance was raised)."""
    cnt = np.diff(start)
    keep_all = np.ones(len(xy), bool)
    todo = (cnt > 3) if tol > 0 else (cnt > max_vertices)
    t = float(tol) if tol > 0 else 0.05
    raised, first = 0, True
    while todo.any():
        idx = np.flatnonzero(todo)
        vsel, keep = _simplify_once(start, xy, idx, t)
        c = cnt[idx]
        left = _sum(keep.astype(np.int64), np.concatenate([[0], np.cumsum(c)]))
        keep |= np.repeat((left < 3) & (c <= max_vertices), c)   # a ring that would lose its area keeps all its vertices
        keep_all[vsel] = keep
        over = left > max_vertices
        if first and tol > 0:
            raised = int(over.sum())
        elif first:
            raised = len(idx)
        first = False
        todo[:] = False
        todo[idx[over]] = True
        t *= 1.5
    s2, xy2 = _keep_vertices(start, xy, keep_all)
    return s2, xy2, raised


def orient(start, xy, area):
    """Reverse the rings whose area is negative: every ring clockwise on screen."""
    rev = area < 0
    if not rev.any():
        return xy
    cnt = np.diff(start)
    ring = _ring_of(start)
    i = np.arange(len(xy))
    return xy[np.where(rev[ring], start[:-1][ring] + start[1:][ring] - 1 - i, i)]


def contains(start, xy, px, py):
    """Is point i inside ring i (even-odd rule)? px, py: one point per ring."""
    if not len(xy):
        return np.zeros(len(start) - 1, bool)
    ring = _ring_of(start)
    nxt = _next(start)
    x1, y1, x2, y2 = xy[:, 0], xy[:, 1], xy[nxt, 0], xy[nxt, 1]
    qx, qy = px[ring], py[ring]
    cross = (y1 > qy) != (y2 > qy)
    with np.errstate(divide="ignore", invalid="ignore"):
        xi = x1 + (qy - y1) * (x2 - x1) / (y2 - y1)
    return (_sum((cross & (qx < xi)).astype(np.int64), start) % 2) == 1


def distance(start, xy, px, py):
    """Distance from point i to the outline of ring i (rings must not be empty)."""
    ring = _ring_of(start)
    nxt = _next(start)
    d = _seg_dist(px[ring], py[ring], xy[:, 0], xy[:, 1], xy[nxt, 0], xy[nxt, 1])
    return np.minimum.reduceat(d, start[:-1])


def self_intersecting(start, xy, block=100_000):
    """Which rings cross themselves (two edges that are not neighbours properly intersect)."""
    cnt = np.diff(start)
    out = np.zeros(len(cnt), bool)
    for k in np.unique(cnt):
        if k < 4:
            continue
        idx = np.flatnonzero(cnt == k)
        for b0 in range(0, len(idx), block):
            sel = idx[b0:b0 + block]
            P = xy[(start[:-1][sel][:, None] + np.arange(k)[None, :])]            # (m, k, 2)
            A, B = P, np.roll(P, -1, axis=1)
            hit = np.zeros(len(sel), bool)

            def side(o, a, b):
                return (a[..., 0] - o[..., 0]) * (b[..., 1] - o[..., 1]) - (a[..., 1] - o[..., 1]) * (b[..., 0] - o[..., 0])
            for i in range(k - 2):
                j = np.arange(i + 2, k if i > 0 else k - 1)                       # edge k-1 is next to edge 0
                if not len(j):
                    continue
                a, b = A[:, i:i + 1], B[:, i:i + 1]
                c, d = A[:, j], B[:, j]
                hit |= ((side(a, b, c) * side(a, b, d) < 0) & (side(c, d, a) * side(c, d, b) < 0)).any(1)
            out[sel] = hit
    return out


# ------------------------------------------------------------------------- label images
def trace_labels(lab):
    """The outline of every label > 0 of a 2-D label image, exactly along its pixel edges.

    Returns (label per ring, start, xy): closed rings (not repeated) on the lattice of pixel
    corners, pixel (row r, column c) covering x in [c, c + 1], y in [r, r + 1]; clockwise on
    screen; corners only (the points along a straight run are left out). A label made of
    several pieces (touching by a corner counts as apart) gives one ring per piece; holes are
    not traced (their count is the fourth value returned).
    All labels are followed at once: every pixel edge between a label and something else is a
    directed edge with the label on its right, an edge's successor is the edge of the same
    label leaving its end point (the rightmost turn first), and the rings are the cycles of
    that map."""
    lab = np.asarray(lab)
    if lab.ndim != 2:
        lab = lab.reshape(lab.shape[-2], lab.shape[-1])
    H, W = lab.shape
    L = np.zeros((H + 2, W + 2), lab.dtype)
    L[1:-1, 1:-1] = lab
    C = L[1:-1, 1:-1]
    parts = []
    # the pixel's top, right, bottom and left edge, walked east, south, west, north
    for d, nb, dx0, dy0 in ((0, L[:-2, 1:-1], 0, 0), (1, L[1:-1, 2:], 1, 0), (2, L[2:, 1:-1], 1, 1), (3, L[1:-1, :-2], 0, 1)):
        r, c = np.nonzero((C != nb) & (C > 0))
        parts.append((C[r, c].astype(np.int64), c.astype(np.int64) + dx0, r.astype(np.int64) + dy0, np.full(len(r), d, np.int64)))
    le = np.concatenate([p[0] for p in parts])
    x0 = np.concatenate([p[1] for p in parts])
    y0 = np.concatenate([p[2] for p in parts])
    d = np.concatenate([p[3] for p in parts])
    E = len(le)
    if not E:
        return np.zeros(0, np.int64), np.zeros(1, np.int64), np.zeros((0, 2)), 0
    x1 = x0 + np.array([1, 0, -1, 0])[d]
    y1 = y0 + np.array([0, 1, 0, -1])[d]
    key = ((le * (H + 1) + y0) * (W + 1) + x0) * 4 + d
    order = np.argsort(key, kind="stable")
    skey = key[order]
    base = ((le * (H + 1) + y1) * (W + 1) + x1) * 4
    nxt = np.full(E, -1, np.int64)
    for turn in (1, 0, 3):                                    # right, straight on, left
        want = base + (d + turn) % 4
        pos = np.minimum(np.searchsorted(skey, want), E - 1)
        hit = (skey[pos] == want) & (nxt < 0)
        nxt[hit] = order[pos[hit]]
    if (nxt < 0).any():
        raise ValueError("label image: an outline does not close")
    # cycles: the smallest edge number of each, then every edge's place after it
    rep = np.arange(E)
    j = nxt.copy()
    while True:
        new = np.minimum(rep, rep[j])
        j = j[j]
        if np.array_equal(new, rep):
            break
        rep = new
    is_rep = rep == np.arange(E)
    dist = np.where(is_rep, 0, 1)
    tgt = np.where(is_rep, np.arange(E), nxt)
    while not is_rep[tgt].all():
        dist = dist + dist[tgt]
        tgt = tgt[tgt]
    length = np.bincount(rep, minlength=E)[rep]
    place = (length - dist) % length
    o = np.lexsort((place, rep))
    le, x0, y0, d, rep = le[o], x0[o], y0[o], d[o], rep[o]
    start = np.concatenate([[0], np.flatnonzero(np.diff(rep)) + 1, [E]])
    xy = np.column_stack([x0, y0]).astype(np.float64)
    # corners only: a vertex stays when the edge leaving it turns
    cnt = np.diff(start)
    prev = np.arange(-1, E - 1)
    prev[start[:-1]] = start[1:] - 1
    start, xy2 = _keep_vertices(start, xy, d != d[prev])
    labels = le[np.concatenate([[0], np.cumsum(cnt)])[:-1]]
    area = signed_area(start, xy2)
    outer = area > 0
    holes = int((~outer).sum())
    start, xy2 = _keep_rings(start, xy2, outer)
    return labels[outer], start, xy2, holes


# ------------------------------------------------------------------------------- readers
def _is_parquet(path) -> bool:
    return str(path).lower().endswith(".parquet")


def _pick(names, *cands):
    low = {n.lower(): n for n in names}
    for c in cands:
        if c.lower() in low:
            return low[c.lower()]
    return None


ID_COLUMNS = ("cell_id", "cell", "EntityID", "cellID", "cell_ID", "id", "label", "name")
X_COLUMNS = ("vertex_x", "x", "x_px", "x_location")
Y_COLUMNS = ("vertex_y", "y", "y_px", "y_location")
PART_COLUMNS = ("label_id", "part", "polygon", "ring")


def read_vertex_table(path, id_col=None, x_col=None, y_col=None, part_col=None):
    """A table with one row per vertex (.csv, .csv.gz, .parquet): an id column, x, y and
    optionally a column that tells the polygons of one cell apart (Xenium `label_id`).
    Rows of one polygon must follow each other. Returns (ids per polygon, start, xy)."""
    import pyarrow as pa
    import pyarrow.compute as pc
    path = Path(path)
    if _is_parquet(path):
        import pyarrow.parquet as pq
        names = list(pq.ParquetFile(str(path)).schema_arrow.names)
    else:
        import flatfiles
        names = flatfiles.header(path)
    idc = id_col or _pick(names, *ID_COLUMNS)
    xc, yc = x_col or _pick(names, *X_COLUMNS), y_col or _pick(names, *Y_COLUMNS)
    pc_ = part_col or _pick(names, *PART_COLUMNS)
    if not (idc and xc and yc) or any(c not in names for c in (idc, xc, yc)):
        raise ValueError(f"{path.name}: needs an id column ({', '.join(ID_COLUMNS[:4])} ...), vertex_x and vertex_y "
                         f"(found {names[:8]})")
    cols = [idc, xc, yc] + ([pc_] if pc_ else [])
    if _is_parquet(path):
        t = pq.read_table(str(path), columns=cols)
    else:
        import pyarrow.csv as pcsv
        t = pcsv.read_csv(str(path), read_options=pcsv.ReadOptions(column_names=names, skip_rows=1),
                          convert_options=pcsv.ConvertOptions(include_columns=cols, column_types={idc: pa.string()}))
    col = t.column(idc)
    if not pa.types.is_string(col.type) and not pa.types.is_large_string(col.type):
        col = pc.cast(col, pa.string())
    enc = pc.dictionary_encode(col).combine_chunks()
    code = enc.indices.to_numpy(zero_copy_only=False).astype(np.int64)
    names_ = enc.dictionary.to_pylist()
    x = t.column(xc).to_numpy().astype(np.float64)
    y = t.column(yc).to_numpy().astype(np.float64)
    brk = np.diff(code) != 0
    if pc_:
        part = pc.dictionary_encode(pc.cast(t.column(pc_), pa.string())).combine_chunks().indices.to_numpy(zero_copy_only=False)
        brk |= np.diff(part) != 0
    first = np.concatenate([[0], np.flatnonzero(brk) + 1]) if len(code) else np.zeros(0, np.int64)
    start = np.concatenate([first, [len(code)]]).astype(np.int64)
    ids = [names_[c] for c in code[first].tolist()]
    return ids, start, np.column_stack([x, y])


def _wkb_outer_rings(buf):
    """Outer rings of a WKB Polygon / MultiPolygon (2-D, or with Z / M which are dropped) and
    the number of holes left out."""
    mv = memoryview(buf)
    rings, holes = [], 0

    def header(off):
        fmt = "<" if mv[off] == 1 else ">"
        t = struct.unpack_from(fmt + "I", mv, off + 1)[0]
        off += 5
        hz, hm = bool(t & 0x80000000), bool(t & 0x40000000)
        if t & 0x20000000:
            off += 4
        base = t & 0x0FFFFFFF
        if base >= 1000:
            hz, hm, base = hz or base // 1000 in (1, 3), hm or base // 1000 in (2, 3), base % 1000
        return fmt, base, 2 + hz + hm, off

    def polygon(fmt, dim, off):
        nonlocal holes
        nr = struct.unpack_from(fmt + "I", mv, off)[0]
        off += 4
        for i in range(nr):
            k = struct.unpack_from(fmt + "I", mv, off)[0]
            off += 4
            if i == 0:
                rings.append(np.frombuffer(mv, dtype=fmt + "f8", count=k * dim, offset=off).reshape(k, dim)[:, :2])
            else:
                holes += 1
            off += 8 * k * dim
        return off

    if len(mv) < 9:
        return rings, holes
    fmt, base, dim, off = header(0)
    if base == 3:
        polygon(fmt, dim, off)
    elif base == 6:
        n = struct.unpack_from(fmt + "I", mv, off)[0]
        off += 4
        for _ in range(n):
            f2, b2, d2, off = header(off)
            if b2 != 3:
                break
            off = polygon(f2, d2, off)
    return rings, holes


def _from_ring_list(ids, rings):
    cnt = np.array([len(r) for r in rings], np.int64)
    start = np.concatenate([[0], np.cumsum(cnt)])
    xy = np.concatenate(rings) if rings else np.zeros((0, 2))
    return ids, start, np.asarray(xy, np.float64)


def read_wkb_table(path, id_col=None, geom_col=None, z_col=None, z=None):
    """A parquet table with one WKB geometry per row (MERSCOPE cell_boundaries.parquet and the
    Vizgen post-processing tool's *_micron_space.parquet; any GeoParquet). With a z column
    (`ZIndex`) each cell keeps the row whose plane is nearest to `z` (the middle plane present
    when `z` is None) among its rows that hold a polygon.
    Returns (ids per polygon, start, xy, info)."""
    import pyarrow as pa
    import pyarrow.compute as pc
    import pyarrow.parquet as pq
    pf = pq.ParquetFile(str(path))
    names = list(pf.schema_arrow.names)
    idc = id_col or _pick(names, "EntityID", *ID_COLUMNS)
    gc = geom_col or _pick(names, "Geometry", "geometry", "geom", "wkb")
    zc = z_col or _pick(names, "ZIndex", "z_index", "z")
    if not (idc and gc):
        raise ValueError(f"{Path(path).name}: needs an id column and a WKB geometry column (found {names[:8]})")
    t = pf.read(columns=[idc, gc] + ([zc] if zc else []))
    ids = t.column(idc).combine_chunks()
    key = ids.cast(pa.string()).to_pylist()
    geom = t.column(gc).combine_chunks()
    size = np.asarray(pc.binary_length(geom).fill_null(0).to_numpy(zero_copy_only=False), np.int64)
    info = {}
    order = np.arange(len(key))
    if zc:
        zz = t.column(zc).to_numpy().astype(np.float64)
        planes = np.unique(zz[np.isfinite(zz)])
        z0 = float(z) if z is not None and len(planes) and z in planes else float(planes[len(planes) // 2]) if len(planes) else 0.0
        info = {"zPlane": int(z0) if float(z0).is_integer() else z0, "zPlanes": [int(p) if float(p).is_integer() else float(p) for p in planes]}
        code = pc.dictionary_encode(ids).indices.to_numpy(zero_copy_only=False).astype(np.int64)
        has = size > 9
        o = np.lexsort((np.abs(zz - z0), ~has, code))
        firsts = o[np.r_[True, code[o][1:] != code[o][:-1]]] if len(o) else o
        order = np.sort(firsts[has[firsts]])
        info["otherPlane"] = int((zz[order] != z0).sum())
    out_ids, rings, holes = [], [], 0
    blobs = geom.to_pylist()
    for i in order.tolist():
        b = blobs[i]
        if b is None:
            continue
        rr, h = _wkb_outer_rings(b)
        holes += h
        for r in rr:
            out_ids.append(key[i])
            rings.append(r)
    info["holes"] = holes
    return (*_from_ring_list(out_ids, rings), info)


def read_geojson(path, id_prop=None):
    """Polygons / MultiPolygons of a GeoJSON file; the cell id is the property `id_prop` (else
    the first of cell_id, cell, EntityID, id, label, name) or the feature's own "id".
    Returns (ids per polygon, start, xy, info)."""
    gj = json.loads(Path(path).read_text(encoding="utf-8"))
    feats = gj.get("features", []) if isinstance(gj, dict) and gj.get("type") == "FeatureCollection" else [gj]
    ids, rings, holes, no_id = [], [], 0, 0
    for f in feats:
        g = (f or {}).get("geometry") or {}
        p = f.get("properties") or {}
        cid = p.get(id_prop) if id_prop else next((p[k] for k in ID_COLUMNS if k in p), f.get("id"))
        polys = [g.get("coordinates")] if g.get("type") == "Polygon" else g.get("coordinates", []) if g.get("type") == "MultiPolygon" else []
        if polys and cid is None:
            no_id += 1
            continue
        for poly in polys:
            if not poly:
                continue
            holes += len(poly) - 1
            ids.append(str(cid))
            rings.append(np.asarray(poly[0], np.float64)[:, :2])
    return (*_from_ring_list(ids, rings), {"holes": holes, "featuresWithoutId": no_id})


def xenium_cell_ids(cell_id) -> list[str]:
    """Xenium cell ids as text: (prefix, suffix) pairs become the eight letters a-p of the
    prefix's hex digits + "-" + the suffix; a plain number stays a number."""
    a = np.asarray(cell_id)
    if a.ndim == 2 and a.shape[1] == 2:
        tr = str.maketrans("0123456789abcdef", "abcdefghijklmnop")
        return [f"{int(p):08x}".translate(tr) + f"-{int(s)}" for p, s in a.tolist()]
    return [str(v.decode() if isinstance(v, bytes) else v) for v in a.ravel().tolist()]


def read_xenium_zarr(path) -> dict:
    """Outlines from Xenium's cells.zarr.zip: {"cell": (ids, start, xy), "nucleus": ...} in µm.
    Newer runs: polygon_sets/<i>/{cell_index, num_vertices, vertices} (set names in the
    `polygon_set_names` attribute, else 0 = nucleus, 1 = cell); older runs: polygon_vertices
    (2, cells, 2 x vertices) + polygon_num_vertices, 0 = nucleus, 1 = cell."""
    import zarr
    p = Path(path)
    store = None
    if p.is_dir():
        root = zarr.open(str(p), mode="r")
    else:
        store = zarr.ZipStore(str(p), mode="r")
        root = zarr.open(store, mode="r")
    try:
        ids = xenium_cell_ids(root["cell_id"][:]) if "cell_id" in root else None
        out = {}

        def rings(v, nv, ci):
            v = np.asarray(v, np.float64).reshape(len(v), -1, 2)
            nv = np.clip(np.asarray(nv, np.int64).ravel(), 0, v.shape[1])
            mask = np.arange(v.shape[1])[None, :] < nv[:, None]
            start = np.concatenate([[0], np.cumsum(nv)])
            who = [ids[i] if ids else str(i) for i in np.asarray(ci, np.int64).ravel().tolist()]
            return who, start, v[mask]
        if "polygon_sets" in root:
            ps = root["polygon_sets"]
            names = ps.attrs.get("polygon_set_names") or root.attrs.get("polygon_set_names") or ["nucleus", "cell"]
            for i, name in enumerate(names):
                if str(i) not in ps:
                    continue
                g = ps[str(i)]
                n = len(g["num_vertices"])
                ci = g["cell_index"][:] if "cell_index" in g else np.arange(n)
                out["nucleus" if "nucle" in str(name).lower() else "cell"] = rings(g["vertices"][:], g["num_vertices"][:], ci)
        elif "polygon_vertices" in root:
            pv, pn = root["polygon_vertices"][:], root["polygon_num_vertices"][:]
            for i, name in enumerate(("nucleus", "cell")):
                out[name] = rings(pv[i], pn[i], np.arange(pv.shape[1]))
        return out
    finally:
        if store is not None:
            store.close()


def read_merscope_hdf5(folder, z=None):
    """Older MERSCOPE runs: cell_boundaries/feature_data_<fov>.hdf5 with
    featuredata/<cell id>/zIndex_<k>/p_<i>/coordinates (1, n, 2) in µm. Each cell's plane
    nearest to `z` (the middle one when None) is read. Returns (ids, start, xy, info)."""
    import h5py
    ids, rings, planes, other = [], [], set(), 0
    files = sorted(Path(folder).glob("feature_data_*.hdf5"))
    for f in files:
        with h5py.File(str(f), "r") as h:
            fd = h.get("featuredata")
            if fd is None:
                continue
            for cid in fd:
                g = fd[cid]
                zs = sorted(int(k.split("_")[1]) for k in g if k.startswith("zIndex_") and len(g[k]))
                if not zs:
                    continue
                planes.update(zs)
                z0 = z if z is not None else zs[len(zs) // 2]
                zk = min(zs, key=lambda v: (abs(v - z0), v))
                other += zk != z0
                for pk in g[f"zIndex_{zk}"]:
                    c = np.asarray(g[f"zIndex_{zk}"][pk]["coordinates"], np.float64).reshape(-1, 2)
                    ids.append(str(cid))
                    rings.append(c)
    zp = z if z is not None else (sorted(planes)[len(planes) // 2] if planes else None)
    return (*_from_ring_list(ids, rings), {"files": len(files), "zPlane": zp, "zPlanes": sorted(planes), "otherPlane": int(other)})


def read_cell_border(path):
    """Stereo-seq cell bin GEF: cellBin/cellBorder, (cells, 32, 2) int16 offsets from the cell's
    x, y (DNB), padded with 32767. Returns (cell x, cell y, start, offsets) in the file's cell
    order, or None when the file has no borders."""
    import h5py
    with h5py.File(str(path), "r") as f:
        cb = f.get("cellBin")
        if cb is None or "cellBorder" not in cb:
            return None
        b = np.asarray(cb["cellBorder"][()])
        cell = cb["cell"][()]
    ok = (b[:, :, 0] != 32767) & (b[:, :, 1] != 32767)
    start = np.concatenate([[0], np.cumsum(ok.sum(1))])
    return cell["x"].astype(np.int64), cell["y"].astype(np.int64), start, b[ok].astype(np.float64)


# ------------------------------------------------------------------------------- writing
def _layer(ds_dir: Path, key: str):
    """(cache xy float32 with NaN = no position, a function giving the barcodes)."""
    import zarr
    g = zarr.open_group(str(Path(ds_dir) / "expr" / f"{key}.zarr"), mode="r")
    xy = np.asarray(g["xy"])

    def barcodes():
        if "barcodes" not in g:
            return None
        return [b.decode() if isinstance(b, bytes) else str(b) for b in g["barcodes"][()]]
    return xy, barcodes


def _pct(a, q):
    return [round(float(v), 3) for v in np.percentile(a, q)] if len(a) else []


def _prepare(s: ShapeSet, xy_cache, barcodes, um: float, tol_um: float, chunk_px: int, log) -> dict | None:
    """Match, clean, simplify and check one set. Returns what the tile writer needs, or None."""
    tag = f"[outlines] {s.key}"
    n_in = len(s.start) - 1
    if not n_in:
        log(f"{tag}: no polygons in {s.source}")
        return None
    on = ~np.isnan(xy_cache[:, 0])
    row_of = np.cumsum(on) - 1
    pos = xy_cache[on].astype(np.float64)
    n_rows = int(on.sum())
    rows = s.rows
    if rows is None:
        bc = barcodes()
        if bc is None:
            log(f"{tag}: the layer was prepared without cell ids (an older cache): outlines cannot be matched; skipped")
            return None
        where = {b: i for i, b in enumerate(bc)}
        rows = np.fromiter((where.get(i, -1) for i in s.ids), np.int64, n_in)
    rows = np.where((rows >= 0) & (rows < len(on)), rows, -1)
    not_in = int((rows < 0).sum())
    if n_in - not_in < MATCH_MIN * n_in:
        log(f"{tag}: only {n_in - not_in:,} of {n_in:,} polygons belong to a cell of the count matrix: these are not this "
            f"dataset's outlines (another run or segmentation?). Not written.")
        return None
    ok = rows >= 0
    ok[ok] = on[rows[ok]]
    no_pos = int(ok.size - ok.sum()) - not_in
    start, xy = _keep_rings(s.start, s.xy, ok)
    rows = row_of[rows[ok]]
    v_in = int(len(xy))
    start, xy = clean_rings(start, xy)
    tol_px = max(tol_um / um if um else 0.0, s.tol_px_min if tol_um > 0 else 0.0)
    start, xy, raised = simplify(start, xy, tol_px)
    area = signed_area(start, xy)
    xy = orient(start, xy, area)
    area = np.abs(area)
    cnt = np.diff(start)
    # degenerate, absurdly large, and crumbs of a cell in several parts
    bad = (cnt < 3) | ~(area > 0)
    ring = _ring_of(start)
    ext = np.zeros(len(cnt))
    if len(xy):
        nz = cnt > 0
        lo = np.full((len(cnt), 2), np.inf)
        hi = np.full((len(cnt), 2), -np.inf)
        lo[nz] = np.minimum.reduceat(xy, start[:-1][nz], axis=0)
        hi[nz] = np.maximum.reduceat(xy, start[:-1][nz], axis=0)
        ext[nz] = (hi[nz] - lo[nz]).max(1)
    large = ~bad & (ext > 2 * chunk_px)
    tot = np.bincount(rows[~bad & ~large], weights=area[~bad & ~large], minlength=n_rows)
    crumb = ~bad & ~large & (area < MIN_PART * tot[rows])
    keep = ~(bad | large | crumb)
    start, xy = _keep_rings(start, xy, keep)
    rows, area, cnt = rows[keep], area[keep], cnt[keep]
    if not len(rows):
        log(f"{tag}: no usable polygon in {s.source}; not written")
        return None
    # the frame check: a cell's own centre lies on its own outline
    cx, cy = pos[rows, 0], pos[rows, 1]
    inside = contains(start, xy, cx, cy)
    dist = np.where(inside, 0.0, distance(start, xy, cx, cy))
    best = np.full(n_rows, np.inf)
    np.minimum.at(best, rows, dist)
    has = np.isfinite(best)
    near_um = 1.0
    in_share = float((best[has] == 0).mean())
    near_share = float((best[has] * (um or 1.0) <= near_um).mean())
    med_d = float(np.sqrt(np.median(np.bincount(rows, weights=area, minlength=n_rows)[has]) / np.pi) * 2)   # px, an outline across
    by_cell = float((best[has] <= med_d).mean())
    nucleus = s.key == "nucleus"
    score = by_cell if nucleus else near_share
    what = (f"{in_share:.1%} of the cell centres lie inside their nucleus outline, {by_cell:.1%} within a nucleus width of it" if nucleus else
            f"{in_share:.2%} of the cell centres lie inside their own outline ({near_share:.2%} within {near_um:g} µm of it)")
    if score < MATCH_MIN:
        log(f"{tag}: {what}: the outlines do not sit on the cells (wrong units, a flipped axis, or another run). Not written.")
        return None
    if score < MATCH_WARN:
        log(f"{tag}: WARNING {what}: check the alignment before trusting these outlines")
    crossing = self_intersecting(start, xy)
    parts = np.bincount(rows, minlength=n_rows)
    stats = {"polygonsRead": n_in, "polygonsNotInMatrix": not_in, "polygonsOfCellsWithoutPosition": no_pos,
             "dropped": {"degenerate": int(bad.sum()), "tooLarge": int(large.sum()), "smallParts": int(crumb.sum()),
                         "holes": int(s.info.get("holes", 0))},
             "rowsInSeveralParts": int((parts > 1).sum()), "selfIntersecting": int(crossing.sum()),
             "centreInside": round(in_share, 5), "centreWithin1um": round(near_share, 5),
             "centreWithinOneDiameter": round(by_cell, 5),
             "verticesRead": v_in, "verticesPerPolygon": _pct(cnt, [0, 50, 100]), "meanVertices": round(float(cnt.mean()), 2),
             "areaUm2": _pct(area * um * um, [25, 50, 75]) if um else None,
             "simplify": {"toleranceUm": tol_um, "tolerancePx": round(tol_px, 4), "maxVertices": MAX_VERTICES,
                          "toleranceRaised": int(raised)}}
    stats.update({k: v for k, v in s.info.items() if k != "holes"})
    log(f"{tag}: {len(rows):,} polygons for {int(has.sum()):,} of {n_rows:,} cells from {s.source} "
        f"({not_in:,} polygons of cells not in the matrix, {no_pos:,} of cells without a position, "
        f"{int(bad.sum()) + int(large.sum()) + int(crumb.sum()):,} dropped); {what}; "
        f"{cnt.mean():.1f} vertices a polygon ({v_in / max(1, n_in - not_in - no_pos):.1f} read), {int(crossing.sum()):,} cross themselves")
    if int((~has).sum()):
        log(f"{tag}:   {int((~has).sum()):,} cells have no {'nucleus outline' if nucleus else 'outline and stay dots'}")
    return {"rows": rows, "start": start, "xy": xy, "pos": pos, "n_rows": n_rows, "stats": stats, "has": has}


def _write_tiles(d: Path, P: dict, chunk_px: int, complete: bool) -> dict:
    """Tiles + what index.json says about them."""
    rows, start, xy, pos, n_rows = P["rows"], P["start"], P["xy"], P["pos"], P["n_rows"]
    cnt = np.diff(start)
    e_row, e_cnt, e_start = rows, cnt, start[:-1]
    if complete:
        miss = np.flatnonzero(~P["has"])
        e_row = np.concatenate([rows, miss])
        e_cnt = np.concatenate([cnt, np.zeros(len(miss), np.int64)])
        e_start = np.concatenate([start[:-1], np.zeros(len(miss), np.int64)])
    p32 = pos[e_row].astype(np.float32).astype(np.float64)          # the float32 values positions.bin holds
    tx = np.floor(p32[:, 0] / chunk_px).astype(np.int64)
    ty = np.floor(p32[:, 1] / chunk_px).astype(np.int64)
    order = np.lexsort((np.arange(len(e_row)), e_row, ty, tx))
    e_row, e_cnt, e_start, tx, ty = e_row[order], e_cnt[order], e_start[order], tx[order], ty[order]
    v_end = np.cumsum(e_cnt)
    v_beg = v_end - e_cnt
    vxy = xy[np.repeat(e_start - v_beg, e_cnt) + np.arange(int(v_end[-1]) if len(v_end) else 0)]
    split = np.flatnonzero((np.diff(tx) != 0) | (np.diff(ty) != 0)) + 1
    t0 = np.concatenate([[0], split])
    t1 = np.concatenate([split, [len(e_row)]])
    # each tile's own origin; the step is the finest power of two under which every tile fits in 16 bits
    origin, extent, over = [], 0.0, 0.0
    for a, b in zip(t0, t1):
        va, vb = int(v_beg[a]), int(v_end[b - 1])
        cx, cy = int(tx[a]) * chunk_px, int(ty[a]) * chunk_px
        if vb > va:
            lo, hi = vxy[va:vb].min(0), vxy[va:vb].max(0)
            o = np.floor(lo)
            extent = max(extent, float((hi - o).max()))
            over = max(over, cx - lo[0], cy - lo[1], hi[0] - (cx + chunk_px), hi[1] - (cy + chunk_px), 0.0)
        else:
            o = np.array([cx, cy], np.float64)
        origin.append(o)
    step = STEP_MIN
    while extent / step > Q:
        step *= 2
    d.mkdir(parents=True, exist_ok=True)
    for old in d.glob("*.bin"):
        old.unlink()
    chunks, total = [], 0
    for (a, b), o in zip(zip(t0, t1), origin):
        va, vb = int(v_beg[a]), int(v_end[b - 1])
        q = np.clip(np.rint((vxy[va:vb] - o) / step), 0, Q).astype("<u2")
        buf = b"".join([np.array([b - a, vb - va], "<u4").tobytes(), o.astype("<f4").tobytes(),
                        e_row[a:b].astype("<u4").tobytes(), e_cnt[a:b].astype("<u2").tobytes(), q.tobytes()])
        name = f"{int(tx[a])}_{int(ty[a])}.bin"
        tmp = d / (name + ".part")
        tmp.write_bytes(buf)
        tmp.replace(d / name)
        chunks.append({"tx": int(tx[a]), "ty": int(ty[a]), "n": int(b - a), "nv": int(vb - va), "bytes": len(buf)})
        total += len(buf)
    return {"step": step, "overhangPx": round(float(over) + step, 3), "chunks": chunks, "bytes": total,
            "nPolygons": int(len(rows)), "nVertices": int(len(xy))}


def write_layer(ds_dir: Path, layer_key: str, sets: list, pixel_size_um: float | None, tol_um: float | None = None,
                log=print) -> dict | None:
    """Write shapes/<layer>/<set>/ in the cache for every set that passes its checks. Returns the
    manifest block for the layer's expr[] entry (`shapes`), or None when nothing was written."""
    t0 = time.time()
    ds_dir = Path(ds_dir)
    tol_um = TOLERANCE_UM if tol_um is None else max(0.0, float(tol_um))
    um = float(pixel_size_um) if pixel_size_um else None
    xy_cache, barcodes = _layer(ds_dir, layer_key)
    chunk_px = int(round(CHUNK_UM / um)) if um else 1024
    base = ds_dir / DIR / layer_key
    blocks, seen = [], set()
    for s in sets:
        if s.key in seen:
            continue
        P = _prepare(s, xy_cache, barcodes, um, tol_um, chunk_px, log)
        if P is None:
            continue
        seen.add(s.key)
        complete = not blocks
        T = _write_tiles(base / s.key, P, chunk_px, complete)
        n_with = int(P["has"].sum())
        index = {"format": FORMAT, "version": VERSION, "set": s.key, "label": s.label,
                 "units": "full-resolution image pixels", "pixelSizeUm": um, "chunkPx": chunk_px,
                 "chunkUm": CHUNK_UM if um else None, "step": T["step"], "overhangPx": T["overhangPx"],
                 "rings": "open (first vertex not repeated), clockwise on screen, outer rings only",
                 "nRows": P["n_rows"], "nRowsWithShape": n_with, "complete": complete,
                 "nPolygons": T["nPolygons"], "nVertices": T["nVertices"], "bytes": T["bytes"],
                 "source": s.source, "stats": P["stats"], "chunks": T["chunks"]}
        deid.assert_clean({k: v for k, v in index.items() if k != "chunks"}, f"shapes/{s.key}/index.json")
        tmp = base / s.key / "index.json.part"
        tmp.write_text(json.dumps(index, separators=(",", ":"), allow_nan=False), encoding="utf-8")
        tmp.replace(base / s.key / "index.json")
        blocks.append({"key": s.key, "label": s.label, "nPolygons": T["nPolygons"], "nRowsWithShape": n_with,
                       "nVertices": T["nVertices"], "chunks": len(T["chunks"]), "bytes": T["bytes"], "source": s.source})
        log(f"[outlines] {s.key}: {len(T['chunks'])} tiles of {chunk_px} px, {T['bytes'] / 1e6:.1f} MB "
            f"({T['bytes'] / max(1, n_with):.0f} bytes a cell), vertex step {T['step']:g} px, reach outside a tile "
            f"{T['overhangPx']:.0f} px  ({time.time() - t0:.0f}s)")
    if not blocks:
        return None                                           # nothing passed: what an earlier run wrote is left alone
    for stale in [p for p in base.iterdir() if p.is_dir() and p.name not in seen]:
        shutil.rmtree(stale)                                  # a set of an earlier run that is not part of this one
    return {"version": VERSION, "format": FORMAT, "dir": DIR, "chunkPx": chunk_px, "sets": blocks}


def attach(ds_dir: Path, meta: dict, sets, pixel_size_um=None, args=None, log=print) -> None:
    """Write `sets` for the layer of manifest entry `meta` and note them on it (meta["shapes"]).
    `sets` empty or None: a message, nothing written. A failure here never stops the import."""
    meta.pop("shapes", None)
    if not sets:
        return
    try:
        block = write_layer(ds_dir, meta["key"], sets, pixel_size_um or meta.get("micronsPerPixel"),
                            tol_um=getattr(args, "outline_tolerance_um", None), log=log)
    except (ValueError, KeyError, OSError, IndexError) as e:
        log(f"[outlines] WARNING could not write the outlines ({type(e).__name__}: {e}); the cells stay dots")
        return
    if block:
        meta["shapes"] = block


def read_set(d: Path):
    """A set back as (index, rows per polygon, start, xy float64 image pixels); entries with no
    polygon come back as empty rings. For tests and checks."""
    d = Path(d)
    index = json.loads((d / "index.json").read_text(encoding="utf-8"))
    rows, cnts, xys = [], [], []
    for c in index["chunks"]:
        buf = (d / f"{c['tx']}_{c['ty']}.bin").read_bytes()
        n, nv = np.frombuffer(buf, "<u4", 2)
        o = np.frombuffer(buf, "<f4", 2, 8).astype(np.float64)
        rows.append(np.frombuffer(buf, "<u4", n, 16))
        cnts.append(np.frombuffer(buf, "<u2", n, 16 + 4 * n))
        q = np.frombuffer(buf, "<u2", 2 * nv, 16 + 6 * n).reshape(-1, 2).astype(np.float64)
        xys.append(o + q * index["step"])
    cnt = np.concatenate(cnts).astype(np.int64) if cnts else np.zeros(0, np.int64)
    return (index, np.concatenate(rows).astype(np.int64) if rows else np.zeros(0, np.int64),
            np.concatenate([[0], np.cumsum(cnt)]), np.concatenate(xys) if xys else np.zeros((0, 2)))


def export_shapes(ds_dir: Path, out: Path, m: dict, log=print) -> dict:
    """Copy the cache's shapes/ into the export (the cache format is the static format).
    Returns {layer key: manifest block} for the layers whose outlines were copied; a set whose
    rows are not the rows of the exported positions.bin is left out with a message."""
    done = {}
    for layer in m.get("expr", []):
        block, key = layer.get("shapes"), layer.get("key")
        if not block:
            continue
        pos = Path(out) / "expr" / key / "positions.bin"
        n_rows = pos.stat().st_size // 8 if pos.exists() else -1
        sets = []
        for s in block.get("sets", []):
            src = Path(ds_dir) / DIR / key / s["key"]
            if not (src / "index.json").exists():
                log(f"[export] outlines {key}/{s['key']}: not in the cache, left out")
                continue
            index = json.loads((src / "index.json").read_text(encoding="utf-8"))
            if index.get("nRows") != n_rows:
                log(f"[export] outlines {key}/{s['key']}: made for {index.get('nRows'):,} rows, the layer has {n_rows:,}: left out "
                    f"(run prepare.py outlines again)")
                continue
            dst = Path(out) / "expr" / key / DIR / s["key"]
            dst.mkdir(parents=True, exist_ok=True)
            names = {f"{c['tx']}_{c['ty']}.bin" for c in index["chunks"]}
            for old in dst.glob("*.bin"):
                if old.name not in names:
                    old.unlink()
            for name in sorted(names) + ["index.json"]:
                shutil.copyfile(src / name, dst / name)
            sets.append(s)
        if sets:
            done[key] = {**block, "sets": sets}
            log(f"[export] outlines {key}: " + ", ".join(f"{s['key']} {s['nPolygons']:,} polygons in {s['chunks']} tiles "
                                                         f"({s['bytes'] / 1e6:.1f} MB)" for s in sets))
    return done


def manifest_layers(layers: list, done: dict) -> list:
    """The export manifest's expr[] with `shapes` only where export_shapes copied them."""
    return [({**L, "shapes": done[L["key"]]} if L.get("key") in done else {k: v for k, v in L.items() if k != "shapes"})
            for L in layers]


# ------------------------------------------------------------------------- platform glue
def _first(folder: Path, *names):
    for n in names:
        if (Path(folder) / n).exists():
            return Path(folder) / n
    return None


def xenium_sets(folder: Path, pixel_size_um: float, log=print) -> list:
    """cell_boundaries / nucleus_boundaries (.parquet, .csv.gz, .csv) of a Xenium output folder,
    else cells.zarr.zip. Xenium's µm frame starts at the image's top-left corner: pixel = µm /
    pixel size, as for the centroids."""
    folder, out = Path(folder), []
    um = float(pixel_size_um)
    zarr_sets = None
    for key in ("cell", "nucleus"):
        f = _first(folder, f"{key}_boundaries.parquet", f"{key}_boundaries.csv.gz", f"{key}_boundaries.csv")
        if f is not None:
            ids, start, xy = read_vertex_table(f, id_col="cell_id", x_col="vertex_x", y_col="vertex_y")
            out.append(ShapeSet(key, start, xy / um, ids=ids, source=f.name))
            continue
        if zarr_sets is None:
            z = _first(folder, "cells.zarr.zip", "cells.zarr")
            try:
                zarr_sets = read_xenium_zarr(z) if z is not None else {}
            except (KeyError, ValueError, OSError) as e:
                log(f"[outlines] WARNING {z.name}: outlines not read ({type(e).__name__}: {e})")
                zarr_sets = {}
        if key in zarr_sets:
            ids, start, xy = zarr_sets[key]
            out.append(ShapeSet(key, start, xy / um, ids=ids, source="cells.zarr.zip"))
    if not out:
        log("[outlines] no cell outlines: cell_boundaries.parquet / .csv.gz or cells.zarr.zip not found; cells are drawn as dots")
    return out


_MERSCOPE_FILES = re.compile(r"^(cell_boundaries|.*micron_space)\.parquet$", re.I)


def merscope_find(inp: dict):
    """The boundaries of a MERSCOPE folder: ("parquet", file) or ("hdf5", folder) or None. Looked
    for beside the cell tables, then beside the transcripts and one folder up."""
    seen, dirs = set(), []
    for k in ("cells", "matrix", "transcripts"):
        p = inp.get(k)
        for d in ([Path(p).parent, Path(p).parent.parent] if p else []):
            if d not in seen:
                seen.add(d)
                dirs.append(d)
    for d in dirs:
        try:
            files = sorted((f for f in d.iterdir() if f.is_file() and _MERSCOPE_FILES.match(f.name)),
                           key=lambda f: (not f.name.lower().startswith("cell_boundaries"), f.name.lower()))
        except OSError:
            continue
        if files:
            return "parquet", files[0]
    for d in dirs:
        h = d / "cell_boundaries"
        if h.is_dir() and any(h.glob("feature_data_*.hdf5")):
            return "hdf5", h
    return None


def merscope_sets(inp: dict, frame, z=None, log=print) -> list:
    """Cell boundaries of a MERSCOPE region in µm, through the same micron -> mosaic pixel
    transform as the centroids (`frame.to_px`)."""
    found = merscope_find(inp)
    if not found:
        log("[outlines] no cell outlines: cell_boundaries.parquet, *_micron_space.parquet or cell_boundaries/*.hdf5 not "
            "found; cells are drawn as dots")
        return []
    kind, path = found
    if kind == "parquet":
        ids, start, xy, info = read_wkb_table(path, z=z)
        source = "cell_boundaries.parquet" if path.name.lower().startswith("cell_boundaries") else "micron_space.parquet"
    else:
        ids, start, xy, info = read_merscope_hdf5(path, z=z)
        source = "cell_boundaries/*.hdf5"
    px, py = frame.to_px(xy[:, 0], xy[:, 1])
    if info.get("otherPlane"):
        log(f"[outlines] {info['otherPlane']:,} cells have no outline in z-plane {info.get('zPlane')}: their nearest plane is used")
    return [ShapeSet("cell", start, np.column_stack([px, py]), ids=ids, source=source, info=info)]


def stereoseq_sets(path, file_offset, offset, image_scale: float, log=print) -> list:
    """cellBin/cellBorder of a cell bin GEF. A border point is a DNB, placed like the cell itself:
    (cell + offset in the border + file offset + 0.5 - image offset) x image pixels per DNB.
    Rows are the file's cells in order, as stereoseq.build_cells writes them."""
    got = read_cell_border(path)
    if got is None:
        log("[outlines] no cell outlines: the cell bin GEF has no cellBin/cellBorder; cells are drawn as dots")
        return []
    cx, cy, start, off = got
    cnt = np.diff(start)
    x = (np.repeat(cx, cnt) + off[:, 0] + file_offset[0] + 0.5 - offset[0]) * image_scale
    y = (np.repeat(cy, cnt) + off[:, 1] + file_offset[1] + 0.5 - offset[1]) * image_scale
    return [ShapeSet("cell", start, np.column_stack([x, y]), rows=np.arange(len(cx)), source="cellBin/cellBorder",
                     info={"note": "border points are DNB centres: the outline runs through the outermost DNBs of the cell"})]


def cosmx_sets(P, cells: dict, log=print) -> list:
    """CosMx outlines: the *polygons file of newer exports when there is one, else traced from
    the CellLabels images. `P` is cosmx.plan()'s Plan (the FOV frame and the mosaic), `cells`
    cosmx.read_cells()'s result (fov and cell_ID of every matrix row)."""
    import flatfiles as ff
    key = np.asarray(cells["fov"], np.int64) * (1 << 32) + np.asarray(cells["cell_ID"], np.int64)
    order = np.argsort(key, kind="stable")
    skey = key[order]

    def row_of(fov, cid):
        k = np.asarray(fov, np.int64) * (1 << 32) + np.asarray(cid, np.int64)
        if not len(skey):
            return np.full(len(k), -1, np.int64)
        p = np.minimum(np.searchsorted(skey, k), len(skey) - 1)
        return np.where(skey[p] == k, order[p], -1)
    poly = P.inp.get("polygons")
    if poly:
        names = ff.header(poly)
        fc, cc = ff.pick(names, "fov", "FOV"), ff.pick(names, "cellID", "cell_ID", "cell_id")
        xc, yc = ff.pick(names, "x_global_px", "x_slide_px"), ff.pick(names, "y_global_px", "y_slide_px")
        if fc and cc and xc and yc:
            t = ff.read_columns(poly, ints=(fc, cc), floats=(xc, yc))
            k = t[fc] * (1 << 32) + t[cc]
            first = np.concatenate([[0], np.flatnonzero(np.diff(k) != 0) + 1]) if len(k) else np.zeros(0, np.int64)
            start = np.concatenate([first, [len(k)]])
            px, py = P.frame.to_px(t[xc], t[yc])
            return [ShapeSet("cell", start, np.column_stack([px, py]), rows=row_of(t[fc][first], t[cc][first]),
                             source="polygons file")]
        log(f"[outlines] {Path(poly).name}: no fov / cellID / x_global_px / y_global_px columns (found {names[:8]}); "
            f"trying the CellLabels images")
    labels = P.inp["images"].get("labels") or {}
    fov = getattr(P, "fov", None)
    if not labels or not fov:
        log("[outlines] no cell outlines: no polygons file and no CellLabels images; cells are drawn as dots")
        return []
    import tifffile
    t0 = time.time()
    (sx, sy), (fx, fy), origins = fov["signs"], fov["flip"], fov["origins"]
    mos = getattr(P, "mosaic", None) or {}
    tiles = mos.get("tiles") or {}
    rows, starts, xys, holes, nv = [], [], [], 0, 0
    tol = TRACE_TOL_PX * abs(P.frame.ax)                          # mask pixels -> image pixels
    for f in sorted(labels):
        if f not in origins:
            continue
        lab = tifffile.imread(str(labels[f]))
        if lab.ndim != 2:
            lab = lab.reshape(lab.shape[-2], lab.shape[-1])
        H, W = lab.shape
        lb, st, xy, h = trace_labels(lab)
        holes += h
        if f in tiles and list(mos.get("fovSize") or []) == [W, H]:
            X, Y = tiles[f][1] + xy[:, 0], tiles[f][2] + xy[:, 1]          # as the FOV's image is pasted
        else:
            lx = W - xy[:, 0] if fx else xy[:, 0]
            ly = H - xy[:, 1] if fy else xy[:, 1]
            X, Y = P.frame.to_px(sx * lx + origins[f][0], sy * ly + origins[f][1])
        st, pts, _ = simplify(st, np.column_stack([X, Y]), tol)
        rows.append(row_of(np.full(len(lb), f), lb))
        starts.append(st[1:] + nv)
        xys.append(pts)
        nv += len(pts)
    if not rows:
        log("[outlines] no cell outlines: no CellLabels image belongs to a FOV with a position")
        return []
    start = np.concatenate([[0]] + starts)
    log(f"[outlines] traced {len(start) - 1:,} outlines from the CellLabels images of {len(rows)} FOVs  ({time.time() - t0:.0f}s)")
    return [ShapeSet("cell", start, np.concatenate(xys), rows=np.concatenate(rows), source="CellLabels images (traced)",
                     tol_px_min=tol, info={"holes": holes, "traced": True, "traceTolerancePx": TRACE_TOL_PX})]


def read_boundaries(path, id_col=None):
    """A lab's own boundaries file by its type: GeoJSON, a parquet with a WKB geometry column, or
    a vertex table (.csv, .csv.gz, .parquet). Returns (ids, start, xy, info, what it was read as)."""
    path = Path(path)
    low = path.name.lower()
    if low.endswith((".geojson", ".json")):
        return (*read_geojson(path, id_col), "GeoJSON polygons")
    if low.endswith(".parquet"):
        import pyarrow as pa
        import pyarrow.parquet as pq
        sch = pq.ParquetFile(str(path)).schema_arrow
        if any(pa.types.is_binary(f.type) or pa.types.is_large_binary(f.type) for f in sch):
            return (*read_wkb_table(path, id_col=id_col), "parquet with WKB geometry")
    if low.endswith((".parquet", ".csv", ".csv.gz")):
        return (*read_vertex_table(path, id_col=id_col), {}, "vertex table")
    raise ValueError(f"{path.name}: expected .geojson, .parquet, .csv or .csv.gz")


def user_sets(cell_file, nucleus_file, units: str, manifest: dict, id_col=None, log=print) -> list:
    """Outlines from a lab's own file(s). `units`: "px" (full-resolution image pixels), "um"
    (µm from the image's top-left corner) or "frame" (the platform's own coordinates, through
    the manifest's `frame` affine)."""
    um = manifest.get("pixelSizeUm")
    A = (manifest.get("frame") or {}).get("affine")
    if units == "um" and not um:
        raise SystemExit("[outlines] --boundary-units um needs the dataset's pixel size")
    if units == "frame" and not A:
        raise SystemExit("[outlines] --boundary-units frame: this dataset's manifest has no `frame` (CosMx / MERSCOPE have one)")
    out = []
    for key, f in (("cell", cell_file), ("nucleus", nucleus_file)):
        if not f:
            continue
        if not Path(f).exists():
            log(f"[outlines] {f}: not found, skipped")
            continue
        ids, start, xy, info, what = read_boundaries(f, id_col)
        if units == "um":
            xy = xy / float(um)
        elif units == "frame":
            xy = np.column_stack([A[0][0] * xy[:, 0] + A[0][1] * xy[:, 1] + A[0][2], A[1][0] * xy[:, 0] + A[1][1] * xy[:, 1] + A[1][2]])
        out.append(ShapeSet(key, start, xy, ids=ids, source=f"{what} (own file)", info=info))
    return out


# -------------------------------------------------------------------------- prepare.py glue
def add_prepare_args(p) -> None:
    """Options of the platform subcommands (xenium, cosmx, merscope, stereoseq)."""
    p.add_argument("--no-outlines", action="store_true", help="skip the cell outlines (cells are drawn as dots)")
    p.add_argument("--outline-tolerance-um", type=float,
                   help=f"outlines are simplified so that no point moves by more than this many µm (default "
                        f"{TOLERANCE_UM}; 0 keeps every vertex of the file)")


def add_args(p) -> None:
    """`prepare.py outlines`: add outlines to a dataset that is already prepared (and exported)."""
    p.add_argument("--id", required=True)
    p.add_argument("--out", required=True, help="cache root that prepare wrote into")
    p.add_argument("--dest", help="export root holding <dest>/<id>/: the outlines are copied there too and its manifest updated")
    p.add_argument("--layer", default="cells", help="the expression layer the outlines belong to (default cells)")
    p.add_argument("--src", help="the platform's output folder (Xenium, MERSCOPE) or cell bin GEF (Stereo-seq) to read the "
                                 "outlines from; default: the one the dataset was prepared from")
    p.add_argument("--boundaries", help="your own cell outlines instead: GeoJSON polygons with the cell id as a property, a "
                                        "table with cell_id, vertex_x, vertex_y (.csv, .csv.gz, .parquet), or a parquet with "
                                        "a WKB Geometry column and an id column")
    p.add_argument("--nuclei", help="your own nucleus outlines, same file types (optional)")
    p.add_argument("--boundary-units", choices=["px", "um", "frame"], default="px",
                   help="coordinates of --boundaries / --nuclei: px = full-resolution image pixels (default), um = µm from "
                        "the image's top-left corner, frame = the platform's own coordinates (manifest `frame`)")
    p.add_argument("--boundary-id", help="the column / property holding the cell id (default: cell_id, cell, EntityID, id ...)")
    p.add_argument("--outline-tolerance-um", type=float,
                   help=f"simplification tolerance in µm (default {TOLERANCE_UM}; 0 keeps every vertex)")


def run(args, log=print) -> None:
    """prepare.py outlines."""
    ds_dir = Path(args.out) / args.id
    mp = ds_dir / "manifest.json"
    if not mp.exists():
        raise SystemExit(f"[outlines] {mp} not found: prepare the dataset first")
    m = json.loads(mp.read_text(encoding="utf-8"))
    layer = next((L for L in m.get("expr", []) if L.get("key") == args.layer), None)
    if layer is None or layer.get("kind") != "cells":
        raise SystemExit(f"[outlines] {args.id} has no cells layer '{args.layer}' (layers: {[L.get('key') for L in m.get('expr', [])]})")
    um, kind, src = m.get("pixelSizeUm"), m.get("kind"), m.get("source") or {}
    if args.boundaries or args.nuclei:
        sets = user_sets(args.boundaries, args.nuclei, args.boundary_units, m, id_col=args.boundary_id, log=log)
    elif kind == "xenium":
        folder = Path(args.src or src.get("folder") or "")
        if folder.is_file():                                     # a zipped Xenium output (xenium.py / zipinput.py)
            import xenium
            zi = xenium.open_zip(folder, Path(args.out) / "_unzipped" / args.id, log=log)
            try:
                sets = xenium.outline_sets(zi, um, log=log)
            finally:
                zi.cleanup()
                zi.close()
        else:
            if not folder.is_dir():
                raise SystemExit(f"[outlines] Xenium folder {folder} not found: give --src")
            if not (folder / "experiment.xenium").exists() and (folder / "data").is_dir():
                folder = folder / "data"
            sets = xenium_sets(folder, um, log=log)
    elif kind == "merscope":
        import merscope
        folder = Path(args.src or src.get("folder") or "")
        if not folder.is_dir():
            raise SystemExit(f"[outlines] MERSCOPE folder {folder} not found: give --src")
        inp = merscope.find_inputs(folder)
        if not any(inp.get(k) for k in ("cells", "matrix", "transcripts")):
            inp["cells"] = folder / "cell_metadata.csv"           # only its folder is used
        fr = m.get("frame") or {}
        if not fr.get("affine"):
            raise SystemExit("[outlines] this MERSCOPE cache has no `frame` in its manifest: prepare it again")
        sets = merscope_sets(inp, merscope.Frame(fr["affine"]), fr.get("zPlane"), log=log)
    elif kind == "stereoseq":
        import h5py
        path = args.src or src.get("cellbin")
        if not path or not Path(path).exists():
            raise SystemExit(f"[outlines] cell bin GEF {path} not found: give --src <SN>.cellbin.gef")
        st = m.get("stereoseq") or {}
        with h5py.File(str(path), "r") as f:
            fo = [int(np.ravel(f.attrs.get(k, 0))[0]) for k in ("offsetX", "offsetY")]
        sets = stereoseq_sets(path, fo, st.get("offsetDnb") or [0, 0], float(st.get("imageScale") or 1.0), log=log)
        if sets and len(sets[0].start) - 1 != layer.get("nBins"):
            raise SystemExit(f"[outlines] {Path(path).name} has {len(sets[0].start) - 1:,} cells, the layer {layer.get('nBins'):,}: "
                             f"not the file this dataset was prepared from")
    else:
        raise SystemExit(f"[outlines] a {kind} dataset: give your outlines with --boundaries (CosMx outlines are traced "
                         f"from the CellLabels images by prepare.py cosmx itself)")
    before = layer.get("shapes")
    attach(ds_dir, layer, sets, um, args=args, log=log)
    if not layer.get("shapes"):
        if before:
            layer["shapes"] = before
        log("[outlines] nothing written" + ("; the outlines this dataset already had are kept" if before else ""))
        return
    mp.write_text(json.dumps(m, indent=1), encoding="utf-8")
    if args.dest:
        out = Path(args.dest) / m["id"]
        ep = out / "manifest.json"
        if not ep.exists():
            raise SystemExit(f"[outlines] {ep} not found: export the dataset first (the cache now has the outlines, so a "
                             f"plain export will include them)")
        em = json.loads(ep.read_text(encoding="utf-8"))
        done = export_shapes(ds_dir, out, m, log=log)
        em["expr"] = [({**L, "shapes": done[L["key"]]} if L.get("key") in done else L) for L in em.get("expr", [])]
        ep.write_text(json.dumps(em, indent=1), encoding="utf-8")
    log(f"[outlines] done -> {ds_dir / DIR}" + (f" and {Path(args.dest) / m['id']}" if args.dest else ""))
