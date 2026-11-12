"""CosMx SMI flat-file export -> image mosaic + cells layer + molecules layer.

Reads, when present (each part is optional), from the folder given as --src or up to two
folders below it:

    *exprMat_file.csv[.gz]         fov, cell_ID, [cell,] one column per target (counts per cell)
    *metadata_file.csv[.gz]        fov, cell_ID, Area, CenterX/Y_local_px, CenterX/Y_global_px, ...
    *fov_positions_file.csv[.gz]   fov, x_px, y_px  (older)  |  FOV, x_global_px, y_global_px  |  ..._mm
    *tx_file.csv[.gz]              fov, cell_ID, x_global_px, y_global_px, x_local_px, y_local_px, z, target, CellComp
    Morphology2D/*_F<fov>.TIF      raw morphology stains, one multi-page 16-bit TIFF per FOV
    RawMorphologyImages/ | ResegmentationImages/ *_F<fov>_Z<z>.TIF   the same per z-plane
    CellComposite/CellComposite_F<fov>.jpg   the vendor's rendered RGB composite per FOV
    CellLabels/CellLabels_F<fov>.tif         cell label images, used only to check orientation

(.parquet copies of the four tables are read too.)

Coordinates. CosMx positions are instrument pixels. "local" is within one FOV image, "global"
is over the slide. Two things differ between exports and are therefore measured, not assumed:

  * global = sign * local + FOV origin, per axis. The sign and every FOV's origin are fitted
    from the cell table (which has both); the FOV positions file gives the FOVs with no cells.
  * where local (0, 0) sits in the FOV image. In the 2021 NSCLC release local y runs UP from
    the bottom edge (image row = height - y_local); this is checked against CellLabels when
    they are there (the label under each centroid must be the cell's own id) and can be forced
    with --fov-y up|down.

The mosaic places every FOV image unflipped at its origin, so mosaic x, y = +-global + constant;
cells and transcripts go through the same two numbers (manifest `frame`).

Cells become an ordinary expression layer (expr.build_cell_layer), transcripts the molecules
layer (molecules.py). Left out of both: negative probes (NegPrb*, Negative*), SystemControl*
and FalseCode* codes; from the cells also cell_ID 0, which is each FOV's background row
(transcripts in no cell), not a cell.
"""
from __future__ import annotations

import re
import time
from collections import OrderedDict
from pathlib import Path

import numpy as np

import flatfiles as ff
from pyramid import Source

PLATFORM = "cosmx"
# µm per instrument pixel by FOV image size (vendor documentation): the 2021 instruments
# (5472 x 3648) and the commercial SMI (4256 x 4256)
PIXEL_UM_BY_FOV = {(5472, 3648): 0.18, (4256, 4256): 0.12028}
DEFAULT_PIXEL_UM = 0.12028
CONTROL_PREFIXES = ("negprb", "negative", "systemcontrol", "falsecode")
IMAGE_DIRS = {"morphology2d": "morphology", "rawmorphologyimages": "zstack", "resegmentationimages": "zstack",
              "cellcomposite": "composite", "celllabels": "labels"}
_FOV_FILE = re.compile(r"(?:^|_)F(\d{1,5})(?:_Z(\d{1,4}))?\.(?:jpe?g|png|tiff?)$", re.I)


def is_control(name: str) -> bool:
    return str(name).strip().lower().startswith(CONTROL_PREFIXES)


# ------------------------------------------------------------------------------- finding
def find_inputs(folder: Path) -> dict:
    """The CosMx files under `folder` (None / {} for each part that is missing)."""
    folder = Path(folder)
    skip = tuple(IMAGE_DIRS) + ("compartmentlabels", "celloverlay")
    inp = {k: ff.find_table(folder, pat, skip_dirs=skip) for k, pat in
           (("matrix", r".*exprmat_file"), ("cells", r".*metadata_file"),
            ("fov_positions", r".*fov_positions_file"), ("transcripts", r".*tx_file"))}
    images: dict[str, dict] = {"morphology": {}, "zstack": {}, "composite": {}, "labels": {}}
    level = [folder]
    for _ in range(5):
        nxt = []
        for d in level:
            try:
                subs = sorted((p for p in d.iterdir() if p.is_dir()), key=lambda p: p.name.lower())
            except OSError:
                continue
            for s in subs:
                kind = IMAGE_DIRS.get(s.name.lower())
                if kind is None:
                    nxt.append(s)
                    continue
                for f in sorted(s.iterdir()):
                    m = _FOV_FILE.search(f.name)
                    if not m:
                        continue
                    fov = int(m.group(1))
                    if kind == "zstack" or (kind == "morphology" and m.group(2)):
                        images["zstack"].setdefault(fov, {})[int(m.group(2) or 0)] = f
                    else:
                        images[kind].setdefault(fov, f)
        level = nxt
    inp["images"] = images
    return inp


def _image_size(path: Path) -> tuple[int, int, int, str]:
    """(width, height, channels, dtype name) of one FOV image without decoding it."""
    if path.suffix.lower() in (".tif", ".tiff"):
        import tifffile
        with tifffile.TiffFile(str(path)) as tf:
            p = tf.pages[0]
            n = len(tf.pages) if p.samplesperpixel == 1 else int(p.samplesperpixel)
            return int(p.imagewidth), int(p.imagelength), n, str(p.dtype)
    from PIL import Image
    with Image.open(path) as im:
        return im.size[0], im.size[1], len(im.getbands()), "uint8"


# ------------------------------------------------------------------------------ geometry
def fit_fov_frame(fov, lx, ly, gx, gy):
    """global = sign * local + origin, per axis: the sign (one for all FOVs) and each FOV's origin.

    Returns ((sx, sy), {fov: (ox, oy)}, worst residual in pixels)."""
    fov = np.asarray(fov)
    order = np.argsort(fov, kind="stable")
    f_sorted = fov[order]
    starts = np.r_[0, np.flatnonzero(np.diff(f_sorted)) + 1]
    ends = np.r_[starts[1:], len(f_sorted)]
    signs, origins, worst = [], {}, 0.0
    for l, g in ((lx, gx), (ly, gy)):
        l, g = np.asarray(l, np.float64)[order], np.asarray(g, np.float64)[order]
        spread = {+1: 0.0, -1: 0.0}
        for s in (+1, -1):
            d = g - s * l
            spread[s] = float(sum(np.ptp(d[a:b]) for a, b in zip(starts, ends)))
        s = +1 if spread[+1] <= spread[-1] else -1
        signs.append(s)
        d = g - s * l
        for a, b in zip(starts, ends):
            o = float(np.median(d[a:b]))
            worst = max(worst, float(np.abs(d[a:b] - o).max()))
            origins.setdefault(int(f_sorted[a]), []).append(o)
    return (signs[0], signs[1]), {f: (o[0], o[1]) for f, o in origins.items()}, worst


def tx_frame_sample(path: Path, max_rows: int = 3_000_000):
    """(fov, local x, local y, global x, global y) of the first transcripts, or None when the
    file does not have both local and global coordinates. Used to fit the FOV frame when the
    cell table cannot give it (no table, or no global centroids)."""
    names = ff.header(path)
    cols = [ff.pick(names, "fov", "FOV"), ff.pick(names, "x_local_px", "x_FOV_px"), ff.pick(names, "y_local_px", "y_FOV_px"),
            ff.pick(names, "x_global_px", "x_slide_px"), ff.pick(names, "y_global_px", "y_slide_px")]
    if not all(cols):
        return None
    parts, n = [], 0
    for b in ff.batches(path, ints=cols[:1], floats=cols[1:]):
        parts.append([ff.col(b, c) for c in cols])
        n += b.num_rows
        if n >= max_rows:
            break
    if not n:
        return None
    return tuple(np.concatenate([p[i] for p in parts]) for i in range(5))


def read_fov_positions(path: Path, pixel_um: float | None) -> dict:
    """{fov: (x, y)} in global pixels from any of the FOV positions layouts."""
    names = ff.header(path)
    fcol = ff.pick(names, "fov", "FOV", "Slide_FOV")
    xcol = ff.pick(names, "x_global_px", "x_px", "X_px", "x")
    ycol = ff.pick(names, "y_global_px", "y_px", "Y_px", "y")
    scale = 1.0
    if not (xcol and ycol):
        xcol = ff.pick(names, "x_global_mm", "x_mm", "X_mm")
        ycol = ff.pick(names, "y_global_mm", "y_mm", "Y_mm")
        scale = 1000.0 / pixel_um if pixel_um else None
    if not (fcol and xcol and ycol and scale):
        return {}
    t = ff.read_columns(path, ints=(fcol,), floats=(xcol, ycol))
    return {int(f): (float(x) * scale, float(y) * scale) for f, x, y in zip(t[fcol], t[xcol], t[ycol])}


def mm_pixel_size(path: Path) -> float | None:
    """µm per pixel when the FOV positions file gives both pixels and millimetres."""
    names = ff.header(path)
    px = ff.pick(names, "x_global_px", "x_px", "X_px")
    mm = ff.pick(names, "x_global_mm", "x_mm", "X_mm")
    if not (px and mm):
        return None
    t = ff.read_columns(path, floats=(px, mm))
    dp, dm = np.ptp(t[px]), np.ptp(t[mm])
    return round(float(dm * 1000.0 / dp), 6) if dp > 0 and dm > 0 else None


def check_labels(label_path: Path, cell_id, lx, ly) -> tuple[bool, bool, float, dict]:
    """Where local (0, 0) sits in the FOV image, from a CellLabels image: the label under each
    centroid must be the cell's id. Returns (flip_x, flip_y, share matched, {name: share})."""
    import tifffile
    lab = tifffile.imread(str(label_path))
    if lab.ndim != 2:
        lab = lab.reshape(lab.shape[-2], lab.shape[-1])
    H, W = lab.shape
    x = np.rint(np.asarray(lx)).astype(np.int64)
    y = np.rint(np.asarray(ly)).astype(np.int64)
    cid = np.asarray(cell_id).astype(np.int64)
    res = {}
    for fx in (False, True):
        for fy in (False, True):
            cc = np.clip(W - x if fx else x, 0, W - 1)
            rr = np.clip(H - y if fy else y, 0, H - 1)
            res[(fx, fy)] = float((lab[rr, cc] == cid).mean()) if len(cid) else 0.0
    (fx, fy), best = max(res.items(), key=lambda kv: kv[1])
    names = {f"x {'from the right' if a else 'from the left'}, y {'up from the bottom' if b else 'down from the top'}": v
             for (a, b), v in res.items()}
    return fx, fy, best, names


class Frame:
    """mosaic pixel = a * global + b per axis (a = +-k, k = image pixels per instrument pixel)."""

    def __init__(self, ax, bx, ay, by):
        self.ax, self.bx, self.ay, self.by = float(ax), float(bx), float(ay), float(by)

    def to_px(self, gx, gy):
        return np.asarray(gx, np.float64) * self.ax + self.bx, np.asarray(gy, np.float64) * self.ay + self.by

    def affine(self):
        return [[self.ax, 0.0, self.bx], [0.0, self.ay, self.by]]


class MosaicSource(Source):
    """Per-FOV images pasted unflipped at their mosaic origins, read one band at a time.

    `tiles`: {fov: (path, x0, y0)} with integer mosaic origins; `planes`: per output channel the
    (page or sample index, "page" | "sample") to read from each file. A decoded plane is kept
    for as long as bands still touch it (least recently used, one row of FOVs and a spare)."""

    def __init__(self, tiles: dict, fov_size, names, planes, dtype, width, height, is_rgb=False):
        self.tiles = tiles
        self.fw, self.fh = fov_size
        self.names = list(names)
        self.planes = planes
        self.dtype = np.dtype(dtype)
        self.width, self.height = int(width), int(height)
        self.is_rgb = is_rgb
        self._cache: OrderedDict = OrderedDict()
        ys = sorted(y0 for _p, _x, y0 in tiles.values())
        per_band = max((sum(1 for y in ys if y < y0 + 1024 and y + self.fh > y0) for y0 in ys), default=1)
        self._cap = per_band + 2
        self._mm: dict = {}
        self._decoded = None                          # (fov, whole decoded RGB image)

    def _plane(self, fov: int, c: int) -> np.ndarray:
        key = (fov, c)
        if key in self._cache:
            self._cache.move_to_end(key)
            return self._cache[key]
        path = self.tiles[fov][0]
        idx, how = self.planes[c]
        if how == "page":
            import tifffile
            try:
                a = tifffile.memmap(str(path), page=idx, mode="r")      # uncompressed: rows on demand
            except Exception:
                a = tifffile.imread(str(path), key=idx)
            if a.ndim == 3:
                a = a[..., 0]
        else:
            if self._decoded is None or self._decoded[0] != fov:
                if path.suffix.lower() in (".tif", ".tiff"):
                    import tifffile
                    full = tifffile.imread(str(path))
                else:
                    from PIL import Image
                    Image.MAX_IMAGE_PIXELS = None
                    with Image.open(path) as im:
                        full = np.asarray(im.convert("RGB"))
                self._decoded = (fov, full)
            a = np.ascontiguousarray(self._decoded[1][..., idx])
        if a.shape != (self.fh, self.fw):
            raise ValueError(f"{path.name}: {a.shape[1]}x{a.shape[0]} px, the other FOV images are {self.fw}x{self.fh}")
        self._cache[key] = a
        while len(self._cache) > self._cap:
            self._cache.popitem(last=False)
        return a

    def read_band(self, c, y0, y1):
        band = np.zeros((y1 - y0, self.width), self.dtype)
        for fov in sorted(self.tiles):
            _p, x0, ty0 = self.tiles[fov]
            a0, a1 = max(y0, ty0), min(y1, ty0 + self.fh)
            if a1 <= a0:
                continue
            xa, xb = max(x0, 0), min(x0 + self.fw, self.width)
            if xb <= xa:
                continue
            src = self._plane(fov, c)
            band[a0 - y0:a1 - y0, xa:xb] = src[a0 - ty0:a1 - ty0, xa - x0:xb - x0]
        return band

    def close(self):
        self._cache.clear()
        self._decoded = None


def overlap_report(tiles: dict, fw: int, fh: int) -> dict:
    """How neighbouring FOV images meet in the mosaic: the largest overlap, and the gaps
    between each FOV and its nearest neighbour to the right / below (pixels)."""
    if len(tiles) < 2:
        return {"maxOverlapPx": 0, "gapX": None, "gapY": None}
    xy = np.array([[x, y] for _p, x, y in tiles.values()], np.float64)
    dx = xy[:, None, 0] - xy[None, :, 0]
    dy = xy[:, None, 1] - xy[None, :, 1]
    ox, oy = fw - np.abs(dx), fh - np.abs(dy)
    both = (ox > 0) & (oy > 0)
    np.fill_diagonal(both, False)
    max_ov = float(np.minimum(ox, oy)[both].max()) if both.any() else 0.0
    same_row = (np.abs(dy) < fh / 2) & (dx > 0)
    same_col = (np.abs(dx) < fw / 2) & (dy > 0)
    gx = np.where(same_row, dx, np.inf).min(axis=0)
    gy = np.where(same_col, dy, np.inf).min(axis=0)
    gx, gy = gx[np.isfinite(gx)] - fw, gy[np.isfinite(gy)] - fh
    return {"maxOverlapPx": max_ov, "gapX": [float(gx.min()), float(gx.max())] if len(gx) else None,
            "gapY": [float(gy.min()), float(gy.max())] if len(gy) else None}


# ---------------------------------------------------------------------------------- cells
def read_cell_table(path: Path) -> dict:
    """The metadata file's columns that are used, under fixed names."""
    names = ff.header(path)
    want = {"fov": ff.pick(names, "fov", "FOV"), "cell_ID": ff.pick(names, "cell_ID", "cell_id"),
            "cell": ff.pick(names, "cell"), "area": ff.pick(names, "Area", "area"),
            "lx": ff.pick(names, "CenterX_local_px", "x_FOV_px"), "ly": ff.pick(names, "CenterY_local_px", "y_FOV_px"),
            "gx": ff.pick(names, "CenterX_global_px", "x_slide_px"),
            "gy": ff.pick(names, "CenterY_global_px", "y_slide_px")}
    if not (want["fov"] and want["cell_ID"]):
        raise SystemExit(f"[cosmx] {path.name}: no fov / cell_ID columns (found {names[:8]} ...)")
    t = ff.read_columns(path, ints=[want[k] for k in ("fov", "cell_ID")],
                        floats=[want[k] for k in ("area", "lx", "ly", "gx", "gy") if want[k]],
                        strings=[want["cell"]] if want["cell"] else [])
    out = {k: t[v] for k, v in want.items() if v}
    out["mean_columns"] = [n[5:] for n in names if n.startswith("Mean.")]
    return out


def cell_key(fov, cell_id) -> np.ndarray:
    return np.asarray(fov, np.int64) * (1 << 32) + np.asarray(cell_id, np.int64)


def read_cells(inp: dict, table: dict | None, log=print) -> dict | None:
    """The count matrix (cells only, controls out) with each row's cell-table row, or None."""
    if not inp.get("matrix") or table is None:
        missing = ([] if inp.get("matrix") else ["*exprMat_file.csv"]) + ([] if table is not None else ["*metadata_file.csv"])
        log(f"[cosmx] no cells layer: {' and '.join(missing)} not found")
        return None
    t0 = time.time()
    X, ids, feats = ff.read_count_matrix(inp["matrix"], int_ids=("fov", "cell_ID", "FOV", "cell_id"), str_ids=("cell",))
    fcol = "fov" if "fov" in ids else "FOV"
    ccol = "cell_ID" if "cell_ID" in ids else "cell_id"
    if fcol not in ids or ccol not in ids:
        raise SystemExit(f"[cosmx] {Path(inp['matrix']).name}: no fov / cell_ID columns")
    fov, cid = ids[fcol], ids[ccol]
    is_cell = cid != 0
    n_bg = int((~is_cell).sum())
    keep = np.array([not is_control(g) for g in feats])
    X = X[np.flatnonzero(is_cell)][:, np.flatnonzero(keep)].tocsr()
    fov, cid = fov[is_cell], cid[is_cell]
    genes = [g for g, k in zip(feats, keep) if k]
    tkey = cell_key(table["fov"], table["cell_ID"])
    order = np.argsort(tkey, kind="stable")
    mkey = cell_key(fov, cid)
    pos = np.searchsorted(tkey[order], mkey)
    pos = np.minimum(pos, len(order) - 1) if len(order) else np.zeros(len(mkey), np.int64)
    hit = (tkey[order][pos] == mkey) if len(order) else np.zeros(len(mkey), bool)
    row = np.where(hit, order[pos] if len(order) else 0, -1)
    if "cell" in ids and any(ids["cell"]):
        bc = [b for b, k in zip(ids["cell"], is_cell) if k]
    else:
        bc = [f"{c}_{f}" for f, c in zip(fov.tolist(), cid.tolist())]      # cell_ID_fov, as Seurat and squidpy name them
    log(f"[cosmx] cells: {len(mkey):,} in the matrix ({n_bg} background rows with cell_ID 0 left out), "
        f"{len(tkey):,} in {Path(inp['cells']).name}, {int(hit.sum()):,} matched; {len(genes)} genes "
        f"({int((~keep).sum())} control targets left out), nnz {X.nnz:,}  ({time.time() - t0:.0f}s)")
    if int((~hit).sum()):
        log(f"[cosmx]   WARNING {int((~hit).sum())} matrix cells are not in the cell table: no position")
    return {"X": X, "genes": genes, "row": row, "barcodes": bc, "nBackgroundRows": n_bg,
            "controls": [g for g, k in zip(feats, keep) if not k], "fov": fov, "cell_ID": cid}


# ----------------------------------------------------------------------------------- plan
class Plan:
    """What prepare.py needs: the image source and pixel size now, the layers after the pyramid."""

    def __init__(self):
        self.source = None
        self.pixel_size_um = None
        self.manifest_source = {}
        self.frame = None
        self.background = None
        self._writer = None
        self._mol_meta = None

    def finish(self, ds_dir: Path, manifest: dict, args, log=print):
        from expr import build_cell_layer
        um = manifest["pixelSizeUm"]
        size = (manifest["image"]["width"], manifest["image"]["height"])
        manifest["frame"] = {"from": "x_global_px, y_global_px (instrument pixels)", "affine": self.frame.affine(),
                             "instrumentPixelUm": self.instrument_um}
        if self.background:
            manifest["image"]["background"] = self.background
            manifest["image"]["channels"][0]["color"] = "#b4b4b4"
        cells = self.cells
        if cells is not None and not args.no_cells:
            t = self.table
            gx, gy = self.cell_global
            px, py = self.frame.to_px(gx, gy)
            xy = np.full((len(cells["row"]), 2), np.nan, np.float64)
            ok = cells["row"] >= 0
            xy[ok, 0] = px[cells["row"][ok]] * um
            xy[ok, 1] = py[cells["row"][ok]] * um
            diam = None
            if "area" in t and len(t["area"]):
                a = float(np.nanmedian(t["area"])) * self.instrument_um ** 2      # px^2 -> µm^2
                diam = 2.0 * (a / np.pi) ** 0.5 if a > 0 else None
            meta = build_cell_layer(ds_dir, xy, um, cells["X"], cells["genes"], key="cells", label="Cells",
                                    cell_diameter_um=diam, log=log, extra={"platform": PLATFORM},
                                    barcodes=cells["barcodes"])
            manifest["expr"].append(meta)
        if not args.no_transcripts:
            if self._writer is None and self.inp.get("transcripts"):
                self._stream(size, um, log)
            if self._writer is not None:
                w = self._writer
                manifest["molecules"] = w.finish(ds_dir, meta={**self._mol_meta, "read": int(w.rows_read)})
            elif not self.inp.get("transcripts"):
                log("[cosmx] no molecules layer: *tx_file.csv not found")

    # transcripts -> MoleculeWriter (pass 1 for the target list, pass 2 for the molecules)
    def _tx_columns(self):
        names = ff.header(self.inp["transcripts"])
        g = ff.pick(names, "target", "gene", "target_name")
        x = ff.pick(names, "x_global_px", "x_slide_px")
        y = ff.pick(names, "y_global_px", "y_slide_px")
        if not (g and x and y):
            raise SystemExit(f"[cosmx] {Path(self.inp['transcripts']).name}: needs target, x_global_px and "
                             f"y_global_px columns (found {names})")
        return g, x, y

    def _scan(self, log):
        if self.scan is None:
            g, x, y = self._tx_columns()
            self.scan = ff.scan_transcripts(self.inp["transcripts"], g, x, y, log=log, tag="cosmx")
        return self.scan

    def _stream(self, size, um, log, density=None):
        g, x, y = self._tx_columns()
        sc = self._scan(log)
        names = sc["names"]
        keep = np.array([bool(n) and not is_control(n) for n in names])
        log(f"[cosmx] transcripts from {Path(self.inp['transcripts']).name}: {int(keep.sum())} genes, "
            f"{int((~keep).sum())} control targets left out")
        self._writer = ff.write_molecules(self.inp["transcripts"], g, x, y, names, keep, self.frame.to_px, size, um,
                                          log=log, tag="cosmx", density=density)
        self._mol_meta = {"platform": PLATFORM, "source": "tx_file"}


def plan(args, log=print) -> Plan:
    """Find the files, fix the coordinate frame and the pixel size, and build the image source."""
    folder = Path(args.src)
    if not folder.is_dir():
        raise SystemExit(f"[cosmx] --src {folder} is not a folder")
    inp = find_inputs(folder)
    P = Plan()
    P.inp, P.scan = inp, None
    for k, what in (("matrix", "*exprMat_file"), ("cells", "*metadata_file"), ("fov_positions", "*fov_positions_file"),
                    ("transcripts", "*tx_file")):
        log(f"[cosmx] {what}: {inp[k].name if inp[k] else 'not found'}")
    img = inp["images"]
    log("[cosmx] images: " + ", ".join(f"{k} {len(v)} FOVs" for k, v in img.items()))

    # ---- the cell table (positions; also fits the FOV frame)
    P.table = table = read_cell_table(inp["cells"]) if inp["cells"] else None
    has_local = table is not None and all(k in table for k in ("lx", "ly"))
    has_global = table is not None and all(k in table for k in ("gx", "gy"))

    # ---- which images
    want = getattr(args, "images", "auto") or "auto"
    z_pick = getattr(args, "z", None)
    stack = {}
    for fov, zs in img["zstack"].items():
        zz = sorted(zs)
        stack[fov] = zs[z_pick] if z_pick in zs else zs[zz[len(zz) // 2]]
    raw = img["morphology"] or stack
    raw_kind = "Morphology2D" if img["morphology"] else ("z-stack images" if stack else None)
    comp = img["composite"]
    if want == "none":
        use, files = None, {}
    elif want == "morphology":
        if not raw:
            raise SystemExit("[cosmx] --images morphology: no Morphology2D / RawMorphologyImages / ResegmentationImages TIFFs found")
        use, files = "morphology", raw
    elif want == "composite":
        if not comp:
            raise SystemExit("[cosmx] --images composite: no CellComposite images found")
        use, files = "composite", comp
    elif raw and len(raw) >= len(comp):
        use, files = "morphology", raw
    elif comp:
        use, files = "composite", comp
        if raw:
            log(f"[cosmx] raw morphology images cover {len(raw)} of {len(comp)} FOVs: using CellComposite "
                f"(--images morphology uses the raw ones and leaves the other FOVs black)")
    else:
        use, files = None, {}
    fov_size = None
    if files:
        w, h, nch, dt = _image_size(files[min(files)])
        fov_size = (w, h)

    # ---- pixel size (µm per instrument pixel)
    um, why = None, ""
    if args.pixel_size:
        um, why = float(args.pixel_size), "--pixel-size"
    if um is None and inp["fov_positions"]:
        um = mm_pixel_size(inp["fov_positions"])
        why = "mm and pixel columns of the FOV positions file"
    if um is None:
        size = fov_size
        if size is None and has_local and len(table["lx"]):
            size = (5472, 3648) if float(np.nanmax(table["lx"])) > 4300 else (4256, 4256)
        if size in PIXEL_UM_BY_FOV:
            um, why = PIXEL_UM_BY_FOV[size], f"{size[0]} x {size[1]} px FOVs"
    if um is None:
        um, why = DEFAULT_PIXEL_UM, "default"
        log(f"[cosmx] WARNING pixel size unknown, assuming {um} µm per pixel: give --pixel-size if this is wrong")
    else:
        log(f"[cosmx] pixel size {um} µm per instrument pixel ({why})")
    P.instrument_um = um

    # ---- FOV frame: global = sign * local + origin
    file_pos = read_fov_positions(inp["fov_positions"], um) if inp["fov_positions"] else {}
    signs, origins = (+1, +1), {}
    fit, fit_from = None, ""
    if has_local and has_global and len(table["fov"]):
        fit, fit_from = (table["fov"], table["lx"], table["ly"], table["gx"], table["gy"]), "the cell table"
    elif inp["transcripts"]:
        fit, fit_from = tx_frame_sample(inp["transcripts"]), "the first transcripts"
    if fit is not None:
        signs, origins, worst = fit_fov_frame(*fit)
        log(f"[cosmx] FOV frame from {fit_from}: global x = {'+' if signs[0] > 0 else '-'}local x + origin, "
            f"global y = {'+' if signs[1] > 0 else '-'}local y + origin; {len(origins)} FOVs, worst residual {worst:.2f} px")
        common = sorted(set(origins) & set(file_pos))
        if common:
            d = np.array([[origins[f][0] - file_pos[f][0], origins[f][1] - file_pos[f][1]] for f in common])
            delta = np.median(d, axis=0)
            log(f"[cosmx] FOV positions file vs fitted origins: offset ({delta[0]:.1f}, {delta[1]:.1f}) px, "
                f"spread {float(np.abs(d - delta).max()):.2f} px over {len(common)} FOVs")
            for f, (x, y) in file_pos.items():
                origins.setdefault(f, (x + float(delta[0]), y + float(delta[1])))
    elif file_pos:
        origins = dict(file_pos)
        log(f"[cosmx] FOV frame from the FOV positions file alone ({len(origins)} FOVs): ASSUMING global = "
            f"local + position on both axes (no table with both local and global coordinates to check)")
    if has_global:
        P.cell_global = (table["gx"], table["gy"])
    elif has_local and origins:
        ox = np.array([origins.get(int(f), (np.nan, np.nan))[0] for f in table["fov"]])
        oy = np.array([origins.get(int(f), (np.nan, np.nan))[1] for f in table["fov"]])
        P.cell_global = (signs[0] * table["lx"] + ox, signs[1] * table["ly"] + oy)
        log("[cosmx] cell table has no global centroids: computed from local centroids and FOV origins")
    else:
        P.cell_global = None
        if table is not None:
            log("[cosmx] WARNING the cell table has no usable centroid columns: no cells layer")

    # ---- where local (0, 0) sits in the FOV image
    fov_y = getattr(args, "fov_y", "auto") or "auto"
    flip_x, flip_y, how = False, True, "assumed: as in the 2021 NSCLC release"
    if fov_y in ("up", "down"):
        flip_y, how = fov_y == "up", "--fov-y"
    elif img["labels"] and has_local:
        fovs, counts = np.unique(table["fov"], return_counts=True)
        cand = [int(f) for f, c in zip(fovs, counts) if c >= 20 and int(f) in img["labels"]]
        if cand:
            f0 = cand[0]
            sel = table["fov"] == f0
            flip_x, flip_y, share, alls = check_labels(img["labels"][f0], table["cell_ID"][sel], table["lx"][sel], table["ly"][sel])
            how = f"CellLabels of FOV {f0}: the cell's own label is under {share:.1%} of its centroids"
            if share < 0.8:
                log(f"[cosmx] WARNING CellLabels do not settle the orientation ({alls}); check the result")
    log(f"[cosmx] local y runs {'up from the bottom' if flip_y else 'down from the top'} edge of the FOV image"
        f"{', local x from the right edge' if flip_x else ''} ({how})")
    tx = signs[0] * (-1 if flip_x else 1)
    ty = signs[1] * (-1 if flip_y else 1)

    P.cells = None
    if not args.no_cells and P.cell_global is not None:
        P.cells = read_cells(inp, table, log=log)
    elif not args.no_cells and table is None:
        read_cells(inp, None, log=log)

    P.manifest_source = {"folder": str(folder), **{k: str(inp[k]) for k in ("matrix", "cells", "fov_positions", "transcripts") if inp[k]}}
    if use:
        fw, fh = fov_size
        tiles_f = {}
        for f in sorted(set(files) | set(origins)):
            if f not in origins:
                continue
            ox, oy = origins[f]
            tiles_f[f] = (tx * (ox + signs[0] * (fw if flip_x else 0)), ty * (oy + signs[1] * (fh if flip_y else 0)))
        no_origin = sorted(set(files) - set(origins))
        if no_origin:
            log(f"[cosmx] WARNING {len(no_origin)} FOV images have no position (not in the cell table or the "
                f"FOV positions file) and are left out: {no_origin[:8]}")
        if not tiles_f:
            raise SystemExit("[cosmx] no FOV has both an image and a position: need *fov_positions_file.csv or "
                             "*metadata_file.csv next to the images (or --images none)")
        bx = -min(v[0] for v in tiles_f.values())
        by = -min(v[1] for v in tiles_f.values())
        tiles = {f: (files[f], int(round(v[0] + bx)), int(round(v[1] + by))) for f, v in tiles_f.items() if f in files}
        W = int(round(max(v[0] for v in tiles_f.values()) + bx)) + fw
        H = int(round(max(v[1] for v in tiles_f.values()) + by)) + fh
        P.frame = Frame(tx, bx, ty, by)
        P.pixel_size_um = um
        if use == "morphology":
            given = [c.strip() for c in (getattr(args, "channels", None) or "").split(",") if c.strip()]
            means = table["mean_columns"] if table is not None else []
            if len(given) == nch:
                names, src_names = given, "--channels"
            elif len(means) == nch:
                names, src_names = list(means), "the Mean.* columns of the cell table, in their order"
            else:
                names, src_names = [f"ch{i}" for i in range(nch)], "no names found; --channels sets them"
            planes = [(i, "page") for i in range(nch)]
            log(f"[cosmx] image: raw morphology ({raw_kind}), {nch} channels named from {src_names}: {', '.join(names)}")
            P.source = MosaicSource(tiles, fov_size, names, planes, dt, W, H)
        else:
            names = ["CellComposite red", "CellComposite green", "CellComposite blue"][:max(nch, 1)]
            planes = [(i, "sample") for i in range(len(names))]
            log(f"[cosmx] image: CellComposite (the vendor's rendered 8-bit colour composite of the stains), "
                f"{len(tiles)} FOVs")
            P.source = MosaicSource(tiles, fov_size, names, planes, "uint8", W, H, is_rgb=True)
        rep = overlap_report(tiles, fw, fh)
        log(f"[cosmx] mosaic {W} x {H} px from {len(tiles)} FOV images of {fw} x {fh}; largest overlap "
            f"{rep['maxOverlapPx']:.0f} px, gap to the next FOV x {rep['gapX']}, y {rep['gapY']} px "
            f"(a later FOV is drawn over an earlier one where they overlap)")
        P.manifest_source["images"] = use
        P.mosaic = {"tiles": {f: (str(p), x, y) for f, (p, x, y) in tiles.items()}, "fovSize": [fw, fh], **rep}
        return P

    # ---- no images: a neutral background from the transcripts (else the cells)
    lims = []
    if inp["transcripts"] and not args.no_transcripts:
        sc = P._scan(log)
        if sc["n"]:
            lims.append((sc["xlim"], sc["ylim"]))
    if P.cell_global is not None and len(P.cell_global[0]):
        gx, gy = P.cell_global
        lims.append(([float(np.nanmin(gx)), float(np.nanmax(gx))], [float(np.nanmin(gy)), float(np.nanmax(gy))]))
    if not lims:
        raise SystemExit(f"[cosmx] nothing to import from {folder}: no images, no transcripts and no cell positions")
    xl = sorted((tx * min(l[0][0] for l in lims), tx * max(l[0][1] for l in lims)))
    yl = sorted((ty * min(l[1][0] for l in lims), ty * max(l[1][1] for l in lims)))
    bg_um = getattr(args, "background_um", None) or ff.BACKGROUND_UM
    W, H, bg_um, x0, y0 = ff.background_grid(xl, yl, 1.0 / um, bg_um=bg_um)
    k = um / bg_um                                            # background pixels per instrument pixel
    P.frame = Frame(tx * k, -x0 * k, ty * k, -y0 * k)
    P.pixel_size_um = bg_um
    dens = ff.Density(W, H)
    if inp["transcripts"] and not args.no_transcripts:
        P._stream((W, H), bg_um, log, density=dens)
        label, sigma = "Transcript density (no image)", 1.5 / bg_um
        P.background = "transcript density"
    else:
        dens.add(*P.frame.to_px(*P.cell_global))
        label, sigma = "Cell density (no image)", 5.0 / bg_um
        P.background = "cell density"
    log(f"[cosmx] no images{' used (--images none)' if want == 'none' else ' found'}: background = "
        f"{P.background} at {bg_um:g} µm per pixel, {W} x {H} px ({dens.n:,} points)")
    P.source = ff.ArraySource([label], dens.image(sigma))
    P.manifest_source["images"] = None
    return P
