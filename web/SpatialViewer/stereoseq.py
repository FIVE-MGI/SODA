"""Stereo-seq (STOmics / BGI, SAW pipeline) output -> square-bin layers + a cells layer.

Reads, straight with h5py / pyarrow (no stereopy, no gefpy):

    <SN>.raw.gef | <SN>.tissue.gef | <SN>.gef     bin GEF (HDF5): geneExp/binN/{gene, expression}
    <SN>.cellbin.gef                              cell bin GEF:   cellBin/{cell, cellExp, gene, cellBorder}
    <SN>.gem | <SN>.tissue.gem.gz                 the same counts as text: geneID [geneName] x y MIDCount
    <SN>_<stain>_regist.tif                       the registered H&E / ssDNA image (built by prepare.py)
    <SN>_<stain>_tissue_cut.tif                   optional tissue mask (same size as the image)

Coordinates. A Stereo-seq chip is a grid of DNB spots `resolution` nm apart (500 nm on every
chip so far: 0.5 µm). The registered image has one pixel per DNB with pixel (0, 0) = chip DNB
(0, 0), x along the columns and y down the rows (no flip). x and y in a file are DNB indices
counted from the file's own offset: STOmics' writer (geftools) stores x - minX, y - minY and
keeps the offset beside them, as `minX` / `minY` on a bin GEF's `expression`, `offsetX` /
`offsetY` on a cell bin GEF's root and `#OffsetX` / `#OffsetY` in a GEM header. A whole-chip
file has offset 0; a file cut to a region (e.g. the tissue's bounding box) has the region's
corner. So

    chip DNB                  = x + file offset
    image pixel (column, row) = (chip DNB + 0.5 - image offset) * imageScale
    µm per image pixel        = resolution / 1000 / imageScale

with image offset 0 and imageScale 1 for a SAW registered image (`--offset` is for an image
that was itself cropped to start at another DNB).

Bin layers. "binN" is an N x N block of DNBs: bin20 = 10 µm, bin50 = 25 µm, bin100 = 50 µm,
bin200 = 100 µm. Bins are counted from the file's origin, bin (i, j) = file x in [i*N, (i+1)*N),
y in [j*N, (j+1)*N), as SAW and Stereopy do (x // N * N), so a bin here is the same bin there.
A GEF level that is already stored (geneExp/binN) is used as it is; a missing one is summed
from the largest stored level that divides it (bin1 in a raw GEF). The GEF is gene-major, so a
pass over it in whole-gene chunks gives the genes x bins CSR directly; the reduced chunks go to
a temporary file per level, so memory is bounded by the largest finished layer, never by the
270 M-row bin1 table. The GEM text is streamed the same way.

Tissue. `raw.gef` covers the whole chip; `tissue.gef` / `tissue.gem.gz` hold only the DNBs
under the tissue mask. For a whole-chip source the cut is made here, from (first found): a
tissue mask TIFF, the GEF's own contour/tissueContour polygons. At bin1 the cut is per DNB
(what SAW's tissue cut does); at a stored coarser level a bin is kept when any of its DNBs is
tissue (exact for a tissue.gef, slightly generous at the border for a raw one).

Layers are written with expr._write_layer / expr.build_cell_layer: positions at bin centres
(cell centroids) in image pixels, `binSizeUm`, and a barcode per row ("x_y" of the bin's first
DNB as the file writes it, e.g. "5300_9500"; the cell id for cells) so the lab's own clusters
can be attached.
"""
from __future__ import annotations

import gzip
import json
import shutil
import time
from pathlib import Path

import numpy as np

DEFAULT_BINS = (20, 50, 100)
DNB_UM_DEFAULT = 0.5             # 500 nm DNB pitch, used only when no file states its resolution
CHUNK = 8_000_000                # records per pass chunk (whole genes; one gene may exceed it)
_GEM_COUNT_COLS = ("midcount", "midcounts", "umicount", "umicounts", "umi_count", "count", "counts", "values")


def add_args(p) -> None:
    """The `prepare.py stereoseq` arguments (the common ones are added by prepare.py)."""
    p.add_argument("--image", required=True, help="the registered image, <SN>_<stain>_regist.tif (H&E RGB or ssDNA grey), "
                                                  "one pixel per DNB")
    p.add_argument("--gef", help="bin GEF: <SN>.raw.gef (whole chip, cut to tissue here), <SN>.tissue.gef or <SN>.gef")
    p.add_argument("--gem", help="GEM text (.gem / .gem.gz), used for the bin layers when no --gef is given")
    p.add_argument("--cellbin", help="cell bin GEF, <SN>.cellbin.gef: the cells layer")
    p.add_argument("--bins", default=",".join(str(b) for b in DEFAULT_BINS),
                   help="comma list of bin sizes in DNBs (default 20,50,100 = 10, 25 and 50 µm at 0.5 µm per DNB; "
                        "200 = 100 µm). Sizes under 8 µm are written but get no QC / UMAP and are not offered "
                        "in the viewer.")
    p.add_argument("--tissue", choices=["auto", "mask", "contour", "none"], default="auto",
                   help="cut the counts to the tissue: auto = --tissue-mask if given, else the GEF's own "
                        "contour/tissueContour, else nothing (a tissue.gef / tissue.gem is already cut); "
                        "none = keep the whole chip")
    p.add_argument("--tissue-mask", help="tissue mask TIFF, <SN>_<stain>_tissue_cut.tif (same size as the image, "
                                         "non-zero = tissue)")
    p.add_argument("--image-scale", type=float, default=None,
                   help="image pixels per DNB (default 1: a SAW registered image)")
    p.add_argument("--offset", help="X,Y: the chip DNB at the image's first pixel (default 0,0: a registered image "
                                    "covers the chip from DNB 0,0). The files' own offsets (minX / minY, "
                                    "offsetX / offsetY, #OffsetX) are read and applied without this.")
    p.add_argument("--dnb-um", type=float, default=None,
                   help="µm between DNBs when no file states it (GEM only; default 0.5)")
    p.add_argument("--barcodes", choices=["xy", "stereopy"], default="xy",
                   help="row names kept for attaching your own clusters: xy = \"x_y\" of the bin's first DNB and "
                        "the cell id for cells; stereopy = the names Stereopy gives bins and cells, x * 2^32 + y")
    p.add_argument("--no-cells", action="store_true", help="skip the cells layer")


def parse_bins(text) -> list[int]:
    try:
        bins = sorted({int(t) for t in str(text).replace(" ", "").split(",") if t})
    except ValueError:
        raise SystemExit(f"--bins wants whole numbers, e.g. 20,50,100 (got {text!r})")
    if not bins or min(bins) < 1:
        raise SystemExit(f"--bins wants whole numbers >= 1 (got {text!r})")
    return bins


# ------------------------------------------------------------------------------ GEF basics
def _s(v) -> str:
    return v.decode("utf-8", "replace") if isinstance(v, (bytes, np.bytes_)) else str(v)


def _attr(obj, name, default=None):
    """An HDF5 attribute as a plain Python value (SAW writes scalars as arrays of one)."""
    if name not in obj.attrs:
        return default
    v = obj.attrs[name]
    if isinstance(v, np.ndarray):
        v = v.reshape(-1)
        if v.size == 0:
            return default
        v = v[0] if v.size == 1 else v.tolist()
    if isinstance(v, (bytes, np.bytes_)):
        return _s(v)
    return v.item() if isinstance(v, np.generic) else v


def gef_kind(path) -> str:
    """"bin", "cellbin" or "unknown" (an HDF5 file with neither group) / "not hdf5"."""
    import h5py
    try:
        with h5py.File(str(path), "r") as f:
            if "cellBin" in f:
                return "cellbin"
            if "geneExp" in f:
                return "bin"
    except OSError:
        return "not hdf5"
    return "unknown"


def bin_levels(f) -> list[int]:
    """Bin sizes stored in a bin GEF (geneExp/binN with gene + expression), ascending."""
    out = []
    for k in f.get("geneExp", {}):
        if k.startswith("bin") and k[3:].isdigit() and "expression" in f["geneExp"][k] and "gene" in f["geneExp"][k]:
            out.append(int(k[3:]))
    return sorted(out)


def gene_labels(ids: list[str], names: list[str]) -> list[str]:
    """One distinct label per gene: its name (the id where there is none). A name several ids
    share - or that differs from another only by case, which a Windows / macOS folder cannot
    hold as two gene files - keeps the plain name for the first id and becomes "Name_<id>" for
    the others, so no gene is merged or lost."""
    out, seen = [], set()
    for gid, name in zip(ids, names):
        lab = name.strip() if name and name.strip() not in ("", "*", "-", "NA", "nan") else gid
        if lab.casefold() in seen:
            lab = f"{lab}_{gid}" if gid and gid != lab else lab
            k = 2
            while lab.casefold() in seen:
                lab = f"{lab}.{k}"
                k += 1
        seen.add(lab.casefold())
        out.append(lab)
    return out


def read_gene_table(ds) -> dict:
    """geneExp/binN/gene or cellBin/gene: ids, names, offsets and row counts.
    GEF version 4: geneID + geneName; earlier versions: one `gene` field holding the name."""
    t = ds[()]
    fields = t.dtype.names or ()
    name_f = "geneName" if "geneName" in fields else ("gene" if "gene" in fields else None)
    id_f = "geneID" if "geneID" in fields else name_f
    if id_f is None:
        raise SystemExit(f"[stereoseq] {ds.name}: no geneID / geneName / gene field (fields: {list(fields)})")
    ids = [_s(v) for v in t[id_f]]
    names = [_s(v) for v in t[name_f]] if name_f else list(ids)
    n_f = next((k for k in ("count", "cellCount") if k in fields), None)
    return {"ids": ids, "names": names, "offset": t["offset"].astype(np.int64),
            "count": t[n_f].astype(np.int64) if n_f else None}


# ----------------------------------------------------------------------------- tissue mask
class TissueMask:
    """Which DNBs are tissue: one bit per DNB (rows packed along x), so a whole 23,520 x 23,520
    chip is 69 MB."""

    def __init__(self, bits: np.ndarray, width: int, height: int, source: str):
        self.bits, self.width, self.height, self.source = bits, int(width), int(height), source

    @property
    def area(self) -> int:
        return int(sum(int(np.unpackbits(self.bits[y0:y0 + 2048]).sum()) for y0 in range(0, self.height, 2048)))

    def contains(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        x = np.asarray(x, np.int64)
        y = np.asarray(y, np.int64)
        ok = (x >= 0) & (x < self.width) & (y >= 0) & (y < self.height)
        out = np.zeros(len(x), bool)
        xs, ys = x[ok], y[ok]
        out[ok] = (self.bits[ys, xs >> 3] >> (7 - (xs & 7)).astype(np.uint8)) & 1
        return out

    def blocks(self, n: int, x0: int = 0, y0: int = 0) -> np.ndarray:
        """(rows, columns) grid of n x n DNB blocks counted from DNB (x0, y0): True where any DNB of
        the block is tissue."""
        nbx, nby = max(0, -(-(self.width - x0) // n)), max(0, -(-(self.height - y0) // n))
        out = np.zeros((nby, nbx), bool)
        starts = np.arange(x0, self.width, n)
        for j in range(nby):
            rows = np.unpackbits(self.bits[y0 + j * n:y0 + (j + 1) * n], axis=1)[:, :self.width]
            out[j] = np.logical_or.reduceat(rows.any(axis=0), starts)
        return out

    @classmethod
    def from_contour(cls, text, width: int, height: int, strip: int = 2048) -> "TissueMask":
        """Fill SAW's contour/tissueContour (GeoJSON-like: geometries -> rings of [x, y] DNBs).

        The rings are the border DNBs of the tissue mask (OpenCV contours), so a ring is filled
        including its outline. A ring that starts on tissue already filled by a larger ring is a
        hole: its inside comes out, its outline (tissue DNBs around the hole) stays."""
        from PIL import Image, ImageDraw
        j = json.loads(_s(text))
        geoms = j.get("geometries") or [f.get("geometry", f) for f in j.get("features", [])] or [j]
        rings = []
        for g in geoms:
            co = g.get("coordinates") or []
            polys = co if g.get("type") == "MultiPolygon" else [co]
            for poly in polys:
                for ring in poly:
                    a = np.asarray(ring, np.int64).reshape(-1, 2)
                    if len(a) >= 3:
                        rings.append(a)
        if not rings:
            raise ValueError("no polygon in the tissue contour")

        def shoelace(a):
            return abs(float(np.dot(a[:, 0], np.roll(a[:, 1], -1)) - np.dot(a[:, 1], np.roll(a[:, 0], -1)))) / 2

        rings.sort(key=shoelace, reverse=True)
        mask = np.zeros((height, width), bool) if height * width <= 64_000_000 else None
        bits = np.zeros((height, (width + 7) // 8), np.uint8)

        def get(y0, y1):
            if mask is not None:
                return mask[y0:y1]
            return np.unpackbits(bits[y0:y1], axis=1)[:, :width].astype(bool)

        def put(y0, y1, m):
            if mask is not None:
                mask[y0:y1] = m
            else:
                bits[y0:y1] = np.packbits(m, axis=1)

        for a in rings:
            x0, y0 = (int(v) for v in a.min(0))
            x1, y1 = (int(v) for v in a.max(0))
            fx, fy = int(a[0, 0]), int(a[0, 1])
            inside = (0 <= fx < width and 0 <= fy < height and bool(get(fy, fy + 1)[0, fx]))
            for s0 in range(max(0, y0), min(height, y1 + 1), strip):
                s1 = min(height, y1 + 1, s0 + strip)
                w = min(width, x1 + 1) - max(0, x0)
                if w <= 0:
                    continue
                ox = max(0, x0)
                pts = [(int(px) - ox, int(py) - s0) for px, py in a]
                im = Image.new("1", (w, s1 - s0), 0)
                ImageDraw.Draw(im).polygon(pts, fill=1, outline=1)
                full = np.asarray(im, dtype=bool)
                cur = get(s0, s1).copy()
                if inside:
                    ed = Image.new("1", (w, s1 - s0), 0)
                    ImageDraw.Draw(ed).polygon(pts, fill=0, outline=1)
                    cur[:, ox:ox + w] &= ~(full & ~np.asarray(ed, dtype=bool))
                else:
                    cur[:, ox:ox + w] |= full
                put(s0, s1, cur)
        if mask is not None:
            bits = np.packbits(mask, axis=1)
        return cls(bits, width, height, f"GEF contour ({len(rings)} ring{'s' if len(rings) != 1 else ''})")

    @classmethod
    def from_image(cls, path, strip: int = 1024) -> "TissueMask":
        """A mask TIFF (non-zero = tissue), read in bands of rows."""
        import tifffile
        import zarr
        with tifffile.TiffFile(str(path)) as tf:
            s = tf.series[0]
            ax = s.axes
            H, W = s.shape[ax.index("Y")], s.shape[ax.index("X")]
            bits = np.zeros((H, (W + 7) // 8), np.uint8)
            z = zarr.open(s.aszarr(), mode="r")
            for y0 in range(0, H, strip):
                idx = tuple(slice(y0, min(H, y0 + strip)) if a == "Y" else slice(None) for a in ax)
                band = np.asarray(z[idx])
                for i, a in reversed(list(enumerate(ax))):
                    if a not in "YX":
                        band = band.max(axis=i)
                bits[y0:y0 + band.shape[0]] = np.packbits(band > 0, axis=1)
        return cls(bits, W, H, "mask image")           # no file name: this text goes into the export manifest


# --------------------------------------------------------------------- bin accumulation
class _Level:
    def __init__(self, n: int, src: int, nbx: int, nby: int, tmp: Path):
        self.n, self.src, self.nbx, self.nby = n, src, nbx, nby
        self.p_idx, self.p_val = tmp / f"bin{n}.idx", tmp / f"bin{n}.val"
        self.f_idx, self.f_val = open(self.p_idx, "wb"), open(self.p_val, "wb")
        self.chunks: list[tuple[int, np.ndarray]] = []      # (first gene, nnz per gene) per reduced chunk
        self.g_hi = -1
        self.ordered = True
        self.nnz = 0
        self.total = 0

    def close(self):
        for f in (self.f_idx, self.f_val):
            if not f.closed:
                f.close()


class BinAccumulator:
    """Sum (gene, x, y, count) records into genes x bins matrices for several bin sizes at once.

    `add` takes a chunk of records (x, y = first DNB of the source bin), reduces it to one
    value per (gene, bin) with a sparse sum and appends the result to a file per level; chunks
    that come gene by gene (a GEF, a gene-sorted GEM) only need concatenating at the end."""

    def __init__(self, bins: list[int], src_of: dict[int, int], width: int, height: int, tmp: Path):
        self.width, self.height = int(width), int(height)
        self.tmp = Path(tmp)
        self.tmp.mkdir(parents=True, exist_ok=True)
        self.levels: dict[int, _Level] = {}
        for n in bins:
            nbx, nby = -(-self.width // n), -(-self.height // n)
            if nbx * nby >= 2 ** 31:
                raise SystemExit(f"[stereoseq] bin{n}: {nbx} x {nby} bins is more than a layer can index")
            self.levels[n] = _Level(n, src_of[n], nbx, nby, self.tmp)

    def add(self, src: int, gene: np.ndarray, x: np.ndarray, y: np.ndarray, count: np.ndarray) -> None:
        from scipy import sparse
        if not len(gene):
            return
        g0, g1 = int(gene.min()), int(gene.max())
        local = (gene - g0).astype(np.int32)
        for L in self.levels.values():
            if L.src != src:
                continue
            b = ((y // L.n) * L.nbx + (x // L.n)).astype(np.int32)
            m = sparse.csr_matrix((count.astype(np.uint32), (local, b)), shape=(g1 - g0 + 1, L.nbx * L.nby))
            m.sum_duplicates()
            m.sort_indices()
            if g0 <= L.g_hi:
                L.ordered = False
            L.g_hi = max(L.g_hi, g1)
            L.chunks.append((g0, np.diff(m.indptr).astype(np.int64)))
            m.indices.astype("<i4").tofile(L.f_idx)
            m.data.astype("<f4").tofile(L.f_val)
            L.nnz += int(m.nnz)
            L.total += int(m.data.sum(dtype=np.int64))

    def finish(self, n: int, n_genes: int):
        """(indptr int64 over genes, bin index on the level's grid int32, value float32)."""
        L = self.levels[n]
        L.close()
        idx = np.fromfile(L.p_idx, dtype="<i4")
        val = np.fromfile(L.p_val, dtype="<f4")
        per_gene = np.zeros(n_genes, np.int64)
        for g0, nn in L.chunks:
            per_gene[g0:g0 + len(nn)] += nn
        if not L.ordered:
            # records did not come gene by gene (an unsorted GEM): sort the reduced chunks by
            # (gene, bin) and sum what several chunks hold for the same pair. Costs memory.
            gene = np.concatenate([np.repeat(np.arange(g0, g0 + len(nn), dtype=np.int64), nn) for g0, nn in L.chunks])
            key = gene * (L.nbx * L.nby) + idx
            del gene
            order = np.argsort(key, kind="stable")
            key, val = key[order], val[order]
            del order
            first = np.flatnonzero(np.concatenate([[True], key[1:] != key[:-1]])) if len(key) else np.zeros(0, np.int64)
            val = np.add.reduceat(val.astype(np.float64), first).astype(np.float32) if len(key) else val
            key = key[first]
            idx = (key % (L.nbx * L.nby)).astype(np.int32)
            per_gene = np.bincount(key // (L.nbx * L.nby), minlength=n_genes).astype(np.int64)
            del key
        indptr = np.zeros(n_genes + 1, np.int64)
        np.cumsum(per_gene, out=indptr[1:])
        for p in (L.p_idx, L.p_val):
            p.unlink(missing_ok=True)
        return indptr, idx, val

    def cleanup(self):
        for L in self.levels.values():
            L.close()
        shutil.rmtree(self.tmp, ignore_errors=True)


def plan_levels(bins: list[int], stored: list[int]) -> dict[int, int]:
    """For each wanted bin size, the stored level to read: itself, else the largest that divides it."""
    out = {}
    for n in bins:
        src = [m for m in stored if n % m == 0]
        if not src:
            raise SystemExit(f"[stereoseq] bin{n} cannot be made: the file holds bin sizes {stored}, none of which "
                             f"divides {n}. Give the raw / tissue GEF (bin1) or change --bins.")
        out[n] = max(src)
    return out


# ---------------------------------------------------------------------------- GEF source
def gef_offset(path) -> tuple[int, int]:
    """(minX, minY) of a bin GEF's finest level: what its stored x, y are counted from."""
    import h5py
    with h5py.File(str(path), "r") as f:
        lv = bin_levels(f)
        if not lv:
            return (0, 0)
        e = f["geneExp"][f"bin{lv[0]}"]["expression"]
        return (int(_attr(e, "minX", 0) or 0), int(_attr(e, "minY", 0) or 0))


def read_contour(path):
    import h5py
    with h5py.File(str(path), "r") as f:
        if "contour" in f and "tissueContour" in f["contour"]:
            v = f["contour"]["tissueContour"][()]
            return v if not isinstance(v, np.ndarray) else v.reshape(-1)[0]
    return None


def gef_chunks(f, level: int, chunk: int = CHUNK):
    """(gene index, x, y, count) arrays for whole genes, about `chunk` records at a time."""
    grp = f["geneExp"][f"bin{level}"]
    gt = read_gene_table(grp["gene"])
    ds = grp["expression"]
    off, cnt = gt["offset"], gt["count"]
    order = np.argsort(off, kind="stable")
    a = 0
    G = len(order)
    while a < G:
        b, n = a, 0
        while b < G and (n == 0 or n + cnt[order[b]] <= chunk):
            n += int(cnt[order[b]])
            b += 1
        genes = order[a:b]
        r0 = int(off[genes[0]])
        r1 = int(off[genes[-1]] + cnt[genes[-1]])
        if r1 - r0 != n:                               # gaps between the genes' row ranges
            parts = [ds[int(off[g]):int(off[g] + cnt[g])] for g in genes if cnt[g]]
            rec = np.concatenate(parts) if parts else ds[0:0]
        else:
            rec = ds[r0:r1]
        a = b
        if not len(rec):
            continue
        yield (np.repeat(genes.astype(np.int32), cnt[genes]), rec["x"].astype(np.int64),
               rec["y"].astype(np.int64), rec["count"].astype(np.int64))


# ---------------------------------------------------------------------------- GEM source
def gem_header(path) -> dict:
    """The `#Key=Value` lines and the column row of a GEM file."""
    p = Path(path)
    op = gzip.open if p.suffix.lower() == ".gz" else open
    comments, n_skip, cols = {}, 0, None
    with op(p, "rt", encoding="utf-8", errors="replace") as f:
        for ln in f:
            if ln.startswith("#"):
                n_skip += 1
                if "=" in ln:
                    k, v = ln[1:].strip().split("=", 1)
                    comments[k.strip()] = v.strip()
                continue
            cols = ln.rstrip("\r\n").split("\t")
            break
    if not cols or len(cols) < 4:
        raise SystemExit(f"[stereoseq] {p.name}: no tab-separated GEM column row (geneID x y MIDCount) found")
    low = [c.strip().lower() for c in cols]
    try:
        ix, iy = low.index("x"), low.index("y")
    except ValueError:
        raise SystemExit(f"[stereoseq] {p.name}: the GEM column row has no x / y columns: {cols}")
    ic = next((low.index(c) for c in _GEM_COUNT_COLS if c in low), None)
    if ic is None:
        raise SystemExit(f"[stereoseq] {p.name}: the GEM column row has no MIDCount / UMICount column: {cols}")
    iid = low.index("geneid") if "geneid" in low else 0
    iname = low.index("genename") if "genename" in low else None

    def num(key, default=0):
        try:
            return int(float(comments.get(key, default)))
        except ValueError:
            return default
    return {"comments": comments, "skip": n_skip, "columns": cols, "id": cols[iid],
            "name": cols[iname] if iname is not None else None, "x": cols[ix], "y": cols[iy], "count": cols[ic],
            "binSize": max(1, num("BinSize", 1)), "offset": [num("OffsetX"), num("OffsetY")],
            "binType": comments.get("BinType")}


def gem_chunks(path, hdr: dict, genes: dict, block_bytes: int = 64 << 20):
    """(gene index, x, y, count) from a GEM, streamed `block_bytes` of text at a time. `genes`
    is filled as ids appear: {"ids": [], "names": [], "index": {}}. Records of the gene a batch
    ends on are held back for the next batch, so a gene-sorted GEM arrives gene by gene (no
    final sort needed)."""
    import pyarrow as pa
    import pyarrow.csv as pcsv
    want = [c for c in (hdr["id"], hdr["name"], hdr["x"], hdr["y"], hdr["count"]) if c]
    types = {hdr["x"]: pa.int64(), hdr["y"]: pa.int64(), hdr["count"]: pa.int64(), hdr["id"]: pa.string()}
    if hdr["name"]:
        types[hdr["name"]] = pa.string()
    rd = pcsv.open_csv(str(path), read_options=pcsv.ReadOptions(skip_rows=hdr["skip"], block_size=int(block_bytes)),
                       parse_options=pcsv.ParseOptions(delimiter="\t"),
                       convert_options=pcsv.ConvertOptions(include_columns=want, column_types=types))
    index = genes["index"]
    carry = None
    for batch in rd:
        if not batch.num_rows:
            continue
        ids = batch.column(hdr["id"]).dictionary_encode()
        d_ids = ids.dictionary.to_pylist()
        code = ids.indices.to_numpy(zero_copy_only=False)
        names = None
        if hdr["name"]:
            nm = batch.column(hdr["name"]).dictionary_encode()
            first = np.full(len(d_ids), -1, np.int64)
            first[code[::-1]] = np.arange(len(code) - 1, -1, -1)       # first row of each id in the batch
            d_nm = nm.dictionary.to_pylist()
            nm_code = nm.indices.to_numpy(zero_copy_only=False)
            names = [d_nm[nm_code[r]] if r >= 0 else None for r in first]
        lut = np.empty(len(d_ids), np.int32)
        for i, gid in enumerate(d_ids):
            k = index.get(gid)
            if k is None:
                k = index[gid] = len(genes["ids"])
                genes["ids"].append(gid)
                genes["names"].append((names[i] if names else None) or gid)
            lut[i] = k
        rec = (lut[code], batch.column(hdr["x"]).to_numpy(zero_copy_only=False),
               batch.column(hdr["y"]).to_numpy(zero_copy_only=False),
               batch.column(hdr["count"]).to_numpy(zero_copy_only=False))
        if carry is not None:
            rec = tuple(np.concatenate([c, r]) for c, r in zip(carry, rec))
            carry = None
        g = rec[0]
        if len(g) > 1 and bool((g[1:] >= g[:-1]).all()):
            cut = int(np.searchsorted(g, g[-1], side="left"))          # start of the last gene's rows
            carry = tuple(r[cut:] for r in rec)
            rec = tuple(r[:cut] for r in rec)
        if len(rec[0]):
            yield rec
    if carry is not None and len(carry[0]):
        yield carry


# -------------------------------------------------------------------------- the bin layers
def _mid(v: int) -> str:
    return f"{int(v):,}"


def build_bin_layers(ds_dir: Path, bins: list[int], gef=None, gem=None, mask: TissueMask | None = None,
                     dnb_um: float = DNB_UM_DEFAULT, image_scale: float = 1.0, offset=(0, 0), image_size=None,
                     barcodes: str = "xy", log=print) -> tuple[list[dict], dict]:
    """Write expr/binNNN.zarr for every wanted bin size from a bin GEF (preferred) or a GEM.
    Returns (manifest layer entries, a summary block of counts)."""
    from scipy import sparse
    from expr import _write_layer
    t0 = time.time()
    ds_dir = Path(ds_dir)
    tmp = ds_dir / "expr" / "_stereoseq_tmp"
    shutil.rmtree(tmp, ignore_errors=True)
    stats = {"records": 0, "mid": 0, "recordsKept": 0, "midKept": 0}
    info: dict = {}
    f = None
    try:
        if gef:
            import h5py
            f = h5py.File(str(gef), "r")
            stored = bin_levels(f)
            src_of = plan_levels(bins, stored)
            e0 = f["geneExp"][f"bin{stored[0]}"]["expression"]
            fo = (int(_attr(e0, "minX", 0) or 0), int(_attr(e0, "minY", 0) or 0))
            for m in sorted(set(src_of.values())):
                em = f["geneExp"][f"bin{m}"]["expression"]
                if (int(_attr(em, "minX", 0) or 0), int(_attr(em, "minY", 0) or 0)) != fo:
                    raise SystemExit(f"[stereoseq] {Path(gef).name}: bin{m} and bin{stored[0]} have different minX / minY")
            ext = (int(_attr(e0, "maxX", 0)) + stored[0], int(_attr(e0, "maxY", 0)) + stored[0])
            gt = read_gene_table(f["geneExp"][f"bin{min(src_of.values())}"]["gene"])
            genes = {"ids": gt["ids"], "names": gt["names"]}
            info = {"source": "gef", "gefVersion": _attr(f, "version"), "levelsStored": stored}
            log(f"[stereoseq] bin GEF version {_attr(f, 'version')}: stored bin sizes {stored}; "
                f"{len(gt['ids']):,} genes; extent x {_attr(e0, 'minX', 0)}-{_attr(e0, 'maxX')}, "
                f"y {_attr(e0, 'minY', 0)}-{_attr(e0, 'maxY')} (DNB), resolution {_attr(e0, 'resolution')} nm"
                + (f"; file offset minX, minY = {list(fo)}: added to every x, y" if fo != (0, 0) else ""))
        else:
            hdr = gem_header(gem)
            if hdr["binType"] and str(hdr["binType"]).lower() not in ("bin", "bins"):
                log(f"[stereoseq] WARNING {Path(gem).name}: #BinType={hdr['binType']} (a square-bin GEM is expected)")
            stored = [hdr["binSize"]]
            src_of = plan_levels(bins, stored)
            genes = {"ids": [], "names": [], "index": {}}
            ext = None
            fo = (int(hdr["offset"][0]), int(hdr["offset"][1]))
            info = {"source": "gem", "gemFormat": hdr["comments"].get("FileFormat"), "levelsStored": stored}
            log(f"[stereoseq] GEM {hdr['comments'].get('FileFormat', '')}: bin size {hdr['binSize']}, columns "
                f"{', '.join(hdr['columns'])}; header #OffsetX, #OffsetY = {list(fo)}"
                + (": added to every x, y" if fo != (0, 0) else ""))
        # Chip DNBs [lo, hi) are on the image (one pixel per DNB); records off it are left out. Bins
        # are counted in the file's own coordinates (chip DNB - file offset), on a grid that reaches
        # the image's far edge.
        lo = (0, 0)
        if image_size and image_scale:
            hi = (int(np.ceil(image_size[0] / image_scale + offset[0])), int(np.ceil(image_size[1] / image_scale + offset[1])))
            lo = (max(0, int(offset[0])), max(0, int(offset[1])))
        elif ext:
            hi = ext
        else:
            raise SystemExit("[stereoseq] a GEM without an image has no known extent")
        width, height = hi[0] - fo[0], hi[1] - fo[1]
        if width <= 0 or height <= 0:
            raise SystemExit(f"[stereoseq] the file's offset {list(fo)} lies beyond the image ({hi[0]} x {hi[1]} DNB): "
                             f"is --image the registered image of this chip?")
        if mask is not None and (mask.width, mask.height) != hi:
            log(f"[stereoseq] WARNING tissue mask is {mask.width} x {mask.height} DNB, the image covers {hi[0]} x {hi[1]}")
        acc = BinAccumulator(bins, src_of, width, height, tmp)
        block_masks = {m: mask.blocks(m, fo[0], fo[1]) for m in sorted(set(src_of.values())) if mask is not None and m > 1}
        for n in bins:
            log(f"[stereoseq] bin{n} ({n * dnb_um:g} µm): " + ("stored level, used as it is" if src_of[n] == n
                else f"summed from bin{src_of[n]}"))

        def feed(src, it):
            n_chunks = 0
            for g, x, y, c in it:
                stats["records"] += len(g)
                stats["mid"] += int(c.sum())
                ax, ay = x + fo[0], y + fo[1]                         # chip DNB
                oob = (x < 0) | (y < 0) | (ax < lo[0]) | (ay < lo[1]) | (ax >= hi[0]) | (ay >= hi[1])
                keep = ~oob
                stats["outside"] = stats.get("outside", 0) + int(oob.sum())
                if mask is not None:
                    if src == 1:
                        keep &= mask.contains(ax, ay)
                    else:
                        bm = block_masks[src]
                        bj, bi = y // src, x // src
                        ok = keep & (bj < bm.shape[0]) & (bi < bm.shape[1])
                        keep[:] = False
                        keep[ok] = bm[bj[ok], bi[ok]]
                if not keep.all():
                    g, x, y, c = g[keep], x[keep], y[keep], c[keep]
                stats["recordsKept"] += len(g)
                stats["midKept"] += int(c.sum())
                acc.add(src, g, x, y, c)
                n_chunks += 1
                if n_chunks % 8 == 0:
                    log(f"[stereoseq]   {stats['records']:,} records read  ({time.time() - t0:.0f}s)")

        per_src = {}
        if gef:
            for src in sorted(set(src_of.values())):
                before = dict(stats)
                feed(src, gef_chunks(f, src))
                per_src[src] = {k: stats[k] - before.get(k, 0) for k in ("records", "mid", "recordsKept", "midKept")}
        else:
            src = stored[0]
            feed(src, gem_chunks(gem, hdr, genes))
            per_src[src] = {k: stats[k] for k in ("records", "mid", "recordsKept", "midKept")}
        for src, s in sorted(per_src.items()):
            log(f"[stereoseq] bin{src} source: {s['records']:,} records, {_mid(s['mid'])} MID read; "
                + (f"{s['recordsKept']:,} records, {_mid(s['midKept'])} MID inside the tissue ({mask.source})"
                   if mask is not None else "no tissue cut applied (whole file kept)"))
        if stats.get("outside"):
            log(f"[stereoseq] WARNING {stats['outside']:,} records lie outside the image (chip DNB {lo[0]}-{hi[0] - 1} x "
                f"{lo[1]}-{hi[1] - 1}) and were left out: is --image the registered image of this chip (or is an "
                f"--offset needed)?")
        labels = gene_labels(genes["ids"], genes["names"])
        n_dup = sum(1 for lab, nm in zip(labels, genes["names"]) if lab != nm and nm)
        if n_dup:
            log(f"[stereoseq] {n_dup} genes share a name with another gene id: shown as Name_<geneID>")
        G = len(labels)
        metas = []
        for n in bins:
            L = acc.levels[n]
            total = L.total
            indptr, idx, val = acc.finish(n, G)
            if not len(idx):
                log(f"[stereoseq] bin{n}: no counts left (nothing in the file, or nothing under the tissue cut): layer skipped")
                continue
            occupied = np.zeros(L.nbx * L.nby, bool)
            for a in range(0, len(idx), 1 << 24):
                occupied[idx[a:a + (1 << 24)]] = True
            cells = np.flatnonzero(occupied)
            lut = np.cumsum(occupied, dtype=np.int64) - 1
            del occupied
            for a in range(0, len(idx), 1 << 24):
                idx[a:a + (1 << 24)] = lut[idx[a:a + (1 << 24)]]
            del lut
            B = len(cells)
            bx, by = (cells % L.nbx).astype(np.int64), (cells // L.nbx).astype(np.int64)
            xy = np.empty((B, 2), np.float32)
            xy[:, 0] = ((bx + 0.5) * n + fo[0] - offset[0]) * image_scale
            xy[:, 1] = ((by + 0.5) * n + fo[1] - offset[1]) * image_scale
            if barcodes == "stereopy":
                names = [str(v) for v in ((bx * n) << 32 | (by * n)).tolist()]
            else:
                names = [f"{a}_{b}" for a, b in zip((bx * n).tolist(), (by * n).tolist())]
            # scipy wants one index type for indices and indptr
            csr = sparse.csr_matrix((val, idx, indptr.astype(np.int32) if len(idx) < 2 ** 31 - 1 else indptr),
                                    shape=(G, B), copy=False)
            um = n * dnb_um
            key = f"bin{n:03d}"
            meta = {"key": key, "label": f"bin{n} ({um:g} µm)", "nBins": int(B), "nGenes": int(G),
                    "binSizeUm": float(um), "spotDiameterPx": float(n * image_scale),
                    "micronsPerPixel": float(dnb_um / image_scale), "imageScale": float(image_scale),
                    "platform": "stereoseq", "binDnb": int(n), "midCount": int(total)}
            _write_layer(ds_dir, xy, csr, labels, meta, names)
            log(f"[stereoseq] {key}: {B:,} bins x {G:,} genes, nnz {len(idx):,}, {_mid(total)} MID  "
                f"({time.time() - t0:.0f}s)")
            metas.append(meta)
            del csr, idx, val, indptr
        if stats.get("outside"):
            info["recordsOffImage"] = int(stats["outside"])
        base = per_src[min(per_src)]              # the finest level read (each stored level read is one pass)
        info.update({"midRead": int(base["mid"]), "midInTissue": int(base["midKept"]),
                     "recordsRead": int(base["records"]), "recordsInTissue": int(base["recordsKept"]),
                     "tissueCut": mask.source if mask is not None else "none (source kept whole)",
                     "nGenes": int(G), "fileOffsetDnb": [int(fo[0]), int(fo[1])]})
        if len(per_src) > 1:
            info["perSourceLevel"] = {f"bin{k}": v for k, v in sorted(per_src.items())}
        return metas, info
    finally:
        if f is not None:
            f.close()
        if "acc" in locals():
            acc.cleanup()
        shutil.rmtree(tmp, ignore_errors=True)


# ------------------------------------------------------------------------------- cells
def read_cellbin(path, log=print) -> dict:
    """cellBin/{cell, cellExp, gene} of a cell bin GEF -> cells x genes CSR, centroids (DNB),
    ids and areas (DNB^2), in the file's cell order."""
    import h5py
    from scipy import sparse
    t0 = time.time()
    with h5py.File(str(path), "r") as f:
        if "cellBin" not in f:
            raise SystemExit(f"[stereoseq] {Path(path).name}: no cellBin group (is this a bin GEF? give it as --gef)")
        cb = f["cellBin"]
        cell = cb["cell"][()]
        gt = read_gene_table(cb["gene"])
        ce = cb["cellExp"][()]
        res = _attr(f, "resolution")
        head = {"version": _attr(f, "version"), "offset": [_attr(f, "offsetX", 0), _attr(f, "offsetY", 0)],
                "hasBorders": "cellBorder" in cb}
    n, G = len(cell), len(gt["ids"])
    off = cell["offset"].astype(np.int64)
    ng = cell["geneCount"].astype(np.int64)
    gidx = ce["geneID"].astype(np.int32)
    data = ce["count"].astype(np.float32)
    del ce
    if n and bool((off[1:] == off[:-1] + ng[:-1]).all()) and off[0] == 0 and off[-1] + ng[-1] == len(gidx):
        indptr = np.concatenate([off, [len(gidx)]])
        X = sparse.csr_matrix((data, gidx, indptr), shape=(n, G))
    else:                                              # rows not laid out cell after cell: gather them
        rows = np.concatenate([np.arange(o, o + k) for o, k in zip(off, ng)]) if n else np.zeros(0, np.int64)
        X = sparse.csr_matrix((data[rows], gidx[rows], np.concatenate([[0], np.cumsum(ng)])), shape=(n, G))
    if len(gidx) and int(gidx.max()) >= G:
        raise SystemExit(f"[stereoseq] {Path(path).name}: cellExp names gene {int(gidx.max())} but there are {G} genes")
    area = cell["area"].astype(np.float64) if "area" in cell.dtype.names else np.full(n, np.nan)
    out = {"X": X, "ids": cell["id"].astype(np.int64), "x": cell["x"].astype(np.int64), "y": cell["y"].astype(np.int64),
           "area": area, "geneIds": gt["ids"], "geneNames": gt["names"], "resolutionNm": res,
           "mid": int(X.data.sum(dtype=np.float64)), **head}
    log(f"[stereoseq] cell bin GEF version {head['version']}: {n:,} cells x {G:,} genes, nnz {X.nnz:,}, "
        f"{_mid(out['mid'])} MID; resolution {res} nm; offset {head['offset']}  ({time.time() - t0:.0f}s)")
    return out


def build_cells(ds_dir: Path, path, dnb_um: float, image_scale: float = 1.0, offset=(0, 0), barcodes: str = "xy",
                image_size=None, log=print, outlines_args=None) -> dict:
    from expr import build_cell_layer
    c = read_cellbin(path, log)
    px_um = dnb_um / image_scale
    fo = (int(c["offset"][0] or 0), int(c["offset"][1] or 0))      # cell x, y are counted from offsetX, offsetY
    if fo != (0, 0):
        log(f"[stereoseq] cell bin GEF offsetX, offsetY = {list(fo)}: added to every cell x, y")
    xy_um = np.column_stack([(c["x"] + fo[0] + 0.5 - offset[0]) * dnb_um, (c["y"] + fo[1] + 0.5 - offset[1]) * dnb_um])
    if image_size:
        out = ((xy_um[:, 0] < 0) | (xy_um[:, 1] < 0) | (xy_um[:, 0] > image_size[0] * px_um) |
               (xy_um[:, 1] > image_size[1] * px_um))
        if out.any():
            log(f"[stereoseq] WARNING {int(out.sum()):,} cells lie outside the image and get no position: is --image "
                f"the registered image of this chip (or is an --offset needed)?")
            xy_um[out] = np.nan
    med = float(np.nanmedian(c["area"])) if np.isfinite(c["area"]).any() else None
    diam = 2.0 * (med * dnb_um * dnb_um / np.pi) ** 0.5 if med else None
    ids = ([str(v) for v in ((c["x"] << 32) | c["y"]).tolist()] if barcodes == "stereopy"
           else [str(v) for v in c["ids"].tolist()])
    meta = build_cell_layer(ds_dir, xy_um, px_um, c["X"], gene_labels(c["geneIds"], c["geneNames"]), key="cells",
                            label="Cells (cell bin)", cell_diameter_um=diam, log=log,
                            extra={"platform": "stereoseq", "midCount": c["mid"],
                                   "medianAreaUm2": round(med * dnb_um * dnb_um, 2) if med else None,
                                   "fileOffsetDnb": [fo[0], fo[1]]},
                            barcodes=ids)
    if not getattr(outlines_args, "no_outlines", False):
        import cellshapes
        cellshapes.attach(ds_dir, meta, cellshapes.stereoseq_sets(path, fo, offset, image_scale, log=log), px_um,
                          args=outlines_args, log=log)
    return meta


# ------------------------------------------------------------------------- prepare.py glue
def resolution_um(args, log=print) -> tuple[float, str]:
    """µm per DNB: the GEF's `resolution` (nm), else the cell bin GEF's, else --dnb-um / 0.5."""
    import h5py
    for path, where in ((args.gef, "geneExp"), (getattr(args, "cellbin", None), "root")):
        if not path or not Path(path).exists():
            continue
        try:
            with h5py.File(str(path), "r") as f:
                if where == "geneExp":
                    lv = bin_levels(f)
                    r = _attr(f["geneExp"][f"bin{lv[0]}"]["expression"], "resolution") if lv else None
                else:
                    r = _attr(f, "resolution")
            if r:
                return float(r) / 1000.0, f"resolution {r} nm in {Path(path).name}"
        except OSError:
            pass
    if args.dnb_um:
        return float(args.dnb_um), "--dnb-um"
    return DNB_UM_DEFAULT, "no file states it: the Stereo-seq DNB pitch of 500 nm is assumed (--dnb-um to change)"


def precheck(args, log=print) -> None:
    """Before the (slow) image pyramid: drop inputs that are missing or of the wrong kind, each with
    a message, and stop now if the wanted bin sizes cannot be made from the file given."""
    if getattr(args, "_stereoseq_checked", False):
        return
    for name in ("gef", "gem", "cellbin", "tissue_mask"):
        p = getattr(args, name, None)
        if p and not Path(p).exists():
            log(f"[stereoseq] --{name.replace('_', '-')} {p}: not found, skipped")
            setattr(args, name, None)
    kinds = {name: gef_kind(getattr(args, name)) for name in ("gef", "cellbin") if getattr(args, name, None)}
    if kinds.get("gef") == "cellbin" and kinds.get("cellbin", "bin") == "bin":
        log(f"[stereoseq] --gef {Path(args.gef).name} is a cell bin GEF"
            + (f" and --cellbin {Path(args.cellbin).name} a bin GEF: swapped" if args.cellbin else ": used as --cellbin"))
        args.gef, args.cellbin = args.cellbin, args.gef
        kinds = {"gef": kinds.get("cellbin"), "cellbin": "cellbin"}
    if args.gef and kinds.get("gef") != "bin":
        log(f"[stereoseq] --gef {Path(args.gef).name} is not a bin GEF (no geneExp group; {kinds.get('gef')}): skipped")
        args.gef = None
    if args.cellbin and kinds.get("cellbin") != "cellbin":
        log(f"[stereoseq] --cellbin {Path(args.cellbin).name} is not a cell bin GEF (no cellBin group; "
            f"{kinds.get('cellbin')}): skipped")
        args.cellbin = None
    bins = parse_bins(args.bins)
    if args.gef:
        import h5py
        with h5py.File(str(args.gef), "r") as f:
            plan_levels(bins, bin_levels(f))
    elif args.gem:
        plan_levels(bins, [gem_header(args.gem)["binSize"]])
    if args.offset:
        try:
            if len([int(v) for v in str(args.offset).split(",")]) != 2:
                raise ValueError
        except ValueError:
            raise SystemExit(f"--offset wants two whole numbers X,Y (got {args.offset!r})")
    args._stereoseq_checked = True


def build(args, ds_dir: Path, manifest: dict, log=print) -> None:
    """Everything after the image pyramid: pixel size, tissue cut, bin layers, cells layer."""
    t0 = time.time()
    ds_dir = Path(ds_dir)
    W, H = manifest["image"]["width"], manifest["image"]["height"]
    precheck(args, log)
    dnb_um, why = resolution_um(args, log)
    scale = args.image_scale
    if scale is None:
        scale = dnb_um / args.pixel_size if args.pixel_size else 1.0
    offset = tuple(int(v) for v in str(args.offset).split(",")) if args.offset else (0, 0)
    px_um = dnb_um / scale
    manifest["pixelSizeUm"] = px_um
    if len(manifest["image"]["channels"]) == 1:        # a grey ssDNA / DAPI image under the bins: shown white, not tinted
        manifest["image"]["channels"][0]["color"] = "#ffffff"
    log(f"[stereoseq] {dnb_um:g} µm per DNB ({why}); image {W} x {H} px at {scale:g} px per DNB "
        f"-> {px_um:g} µm per pixel; offset {list(offset)}")
    bins = parse_bins(args.bins)
    block = {"dnbUm": dnb_um, "imageScale": scale, "offsetDnb": list(offset)}

    if args.gef or args.gem:
        mask = None
        if args.tissue != "none":
            if args.tissue_mask and args.tissue in ("auto", "mask"):
                mask = TissueMask.from_image(args.tissue_mask)
            elif args.tissue == "mask":
                log("[stereoseq] --tissue mask needs --tissue-mask <tif>: no tissue cut applied")
            elif args.gef and args.tissue == "auto" and gef_offset(args.gef) != (0, 0):
                log(f"[stereoseq] {Path(args.gef).name} is cut to a region (minX, minY = {list(gef_offset(args.gef))}): taken as "
                    f"already cut to the tissue, its contour is not applied (--tissue contour applies it, in chip DNBs)")
            elif args.gef:
                contour = read_contour(args.gef)
                if contour is not None:
                    gw, gh = int(np.ceil(W / scale + offset[0])), int(np.ceil(H / scale + offset[1]))
                    try:
                        mask = TissueMask.from_contour(contour, gw, gh)
                    except (ValueError, KeyError, TypeError) as e:
                        log(f"[stereoseq] WARNING the GEF's tissue contour could not be read ({e}): no tissue cut")
                elif args.tissue == "contour":
                    log(f"[stereoseq] --tissue contour: {Path(args.gef).name} has no contour/tissueContour: no tissue cut")
        if mask is not None:
            which = f" {Path(args.tissue_mask).name}" if mask.source == "mask image" else ""
            log(f"[stereoseq] tissue from {mask.source}{which}: {mask.area:,} of {mask.width * mask.height:,} DNBs  "
                f"({time.time() - t0:.0f}s)")
        else:
            log("[stereoseq] no tissue mask or contour: every record of the source is kept (a tissue.gef / "
                "tissue.gem is already cut to the tissue; a raw.gef is the whole chip)")
        if args.gef and args.gem:
            log(f"[stereoseq] both --gef and --gem given: counts are read from the GEF, {Path(args.gem).name} is not used"
                + ("" if mask is not None or args.tissue == "none" else
                   " (if the GEF is a whole chip and the GEM the tissue, leave out --gef to read the GEM instead)"))
        metas, info = build_bin_layers(ds_dir, bins, gef=args.gef, gem=None if args.gef else args.gem, mask=mask,
                                       dnb_um=dnb_um, image_scale=scale, offset=offset, image_size=(W, H),
                                       barcodes=args.barcodes, log=log)
        manifest["expr"].extend(metas)
        block.update(info)
    else:
        log("[stereoseq] no --gef / --gem: no bin layers")

    if args.cellbin and not args.no_cells:
        try:
            manifest["expr"].append(build_cells(ds_dir, args.cellbin, dnb_um, scale, offset, args.barcodes,
                                                image_size=(W, H), log=log, outlines_args=args))
        except KeyError as e:
            log(f"[stereoseq] WARNING {Path(args.cellbin).name}: dataset {e} missing, no cells layer")
    elif not args.cellbin:
        log("[stereoseq] no --cellbin: no cells layer")
    manifest["stereoseq"] = block
    log(f"[stereoseq] layers done in {time.time() - t0:.0f}s")
