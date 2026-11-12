"""Shared pieces of the flat-file readers of imaging platforms (cosmx.py, merscope.py).

CosMx and MERSCOPE both export plain tables: a dense cells x genes count matrix, a cell table
with centroids, and one row per decoded transcript (tens to hundreds of millions of rows), as
.csv, .csv.gz or .parquet. This module reads them in chunks with pyarrow (bounded memory, no
pandas), turns the count matrix into a sparse one, streams transcripts into
molecules.MoleculeWriter, and builds the neutral background used when a dataset has no images.

    find_table(folder, pattern)          the first <pattern>.csv | .csv.gz | .parquet under a folder
    header(path)                         column names (an unnamed index column becomes "_c0")
    batches(path, ints, floats, strings) record batches of the named columns only
    read_columns(path, ...)              the same columns, whole (cell tables)
    read_count_matrix(path, id_columns)  (X cells x features CSR float32, id columns, feature names)
    scan_transcripts(...)                pass 1 over a transcript table: features, counts, extent
    write_molecules(...)                 pass 2: transcripts -> MoleculeWriter (+ a density image)
    Density / ArraySource                transcript- or cell-density image as a pyramid source
"""
from __future__ import annotations

import csv
import gzip
import re
import time
from pathlib import Path

import numpy as np

from pyramid import Source

TABLE_SUFFIXES = (".csv", ".csv.gz", ".parquet")
# CSV read block. pyarrow parses several blocks ahead in parallel, so peak memory follows this,
# not the file size: 8 MB kept a 26-million-row transcript table under 0.6 GB (32 MB: 1.6 GB) at
# the same speed.
BLOCK_BYTES = 8 << 20
PARQUET_ROWS = 1_000_000
BACKGROUND_UM = 1.0               # µm per pixel of the neutral background (no images)
BACKGROUND_MAX_PX = 250e6         # above this many pixels the background gets coarser pixels


def is_parquet(path) -> bool:
    return str(path).lower().endswith(".parquet")


def find_table(folder: Path, pattern: str, depth: int = 2, skip_dirs=()) -> Path | None:
    """First file whose name matches `pattern` (a regex for the part before .csv / .csv.gz /
    .parquet, case-insensitive), looking in `folder` first and then up to `depth` folders down."""
    rx = re.compile(rf"(?:{pattern})\.(?:parquet|csv\.gz|csv)$", re.I)
    level = [Path(folder)]
    for _ in range(depth + 1):
        nxt = []
        for d in level:
            try:
                entries = sorted(d.iterdir(), key=lambda p: p.name.lower())
            except OSError:
                continue
            hits = [p for p in entries if p.is_file() and rx.search(p.name)]
            if hits:
                # parquet before csv.gz before csv when a folder holds several copies
                rank = {".parquet": 0, ".gz": 1, ".csv": 2}
                return sorted(hits, key=lambda p: (rank.get(p.suffix.lower(), 3), p.name.lower()))[0]
            nxt += [p for p in entries if p.is_dir() and p.name.lower() not in skip_dirs]
        level = nxt
    return None


def header(path: Path) -> list[str]:
    """Column names. Empty names (the unnamed index column pandas writes) become "_c<i>", and a
    repeated name gets "_<i>" appended, so every column can be asked for by name."""
    path = Path(path)
    if is_parquet(path):
        import pyarrow.parquet as pq
        names = list(pq.ParquetFile(str(path)).schema_arrow.names)
    else:
        op = gzip.open if path.name.lower().endswith(".gz") else open
        with op(path, "rt", encoding="utf-8-sig", newline="") as f:
            names = next(csv.reader(f), [])
    out, seen = [], set()
    for i, n in enumerate(names):
        n = str(n).strip()
        if not n:
            n = f"_c{i}"
        if n in seen:
            n = f"{n}_{i}"
        seen.add(n)
        out.append(n)
    return out


def is_index_column(name: str) -> bool:
    """A row-number column left by a re-save (pandas `to_csv` without index=False)."""
    return bool(re.match(r"^(_c\d+|unnamed: ?\d+(_\d+)?|index|__index_level_\d+__)$", name.strip(), re.I))


def pick(names: list[str], *candidates: str) -> str | None:
    """The first of `candidates` present in `names`, compared case-insensitively."""
    low = {n.lower(): n for n in names}
    for c in candidates:
        if c.lower() in low:
            return low[c.lower()]
    return None


def batches(path: Path, ints=(), floats=(), strings=(), block_bytes: int = BLOCK_BYTES):
    """Record batches holding only the named columns, with fixed types, of a CSV (.gz ok) or
    Parquet table. Memory is one block of the file at a time."""
    import pyarrow as pa
    path = Path(path)
    types = {**{c: pa.int64() for c in ints}, **{c: pa.float64() for c in floats},
             **{c: pa.string() for c in strings}}
    cols = list(types)
    if is_parquet(path):
        import pyarrow.parquet as pq
        pf = pq.ParquetFile(str(path))
        schema = pa.schema([(c, types[c]) for c in cols])
        for b in pf.iter_batches(batch_size=PARQUET_ROWS, columns=cols):
            b = pa.RecordBatch.from_arrays([b.column(b.schema.get_field_index(c)) for c in cols], names=cols)
            yield b.cast(schema) if b.schema != schema else b
        return
    import pyarrow.csv as pcsv
    rd = pcsv.open_csv(str(path),
                       read_options=pcsv.ReadOptions(column_names=header(path), skip_rows=1, block_size=block_bytes),
                       convert_options=pcsv.ConvertOptions(include_columns=cols, column_types=types,
                                                           strings_can_be_null=False))
    try:
        for b in rd:
            if b.num_rows:
                yield b
    finally:
        rd.close()


def col(batch, name: str) -> np.ndarray:
    """One numeric column of a batch as numpy (nulls -> NaN for floats)."""
    return batch.column(batch.schema.get_field_index(name)).to_numpy(zero_copy_only=False)


def strs(batch, name: str) -> list[str]:
    return ["" if v is None else v for v in batch.column(batch.schema.get_field_index(name)).to_pylist()]


def read_columns(path: Path, ints=(), floats=(), strings=()) -> dict:
    """Whole columns of a (cell) table: numpy arrays for numbers, lists for strings."""
    parts: dict[str, list] = {c: [] for c in (*ints, *floats, *strings)}
    for b in batches(path, ints, floats, strings):
        for c in (*ints, *floats):
            parts[c].append(col(b, c))
        for c in strings:
            parts[c].append(strs(b, c))
    out = {}
    for c in ints:
        out[c] = np.concatenate(parts[c]).astype(np.int64) if parts[c] else np.zeros(0, np.int64)
    for c in floats:
        out[c] = np.concatenate(parts[c]).astype(np.float64) if parts[c] else np.zeros(0, np.float64)
    for c in strings:
        out[c] = [v for p in parts[c] for v in p]
    return out


def read_count_matrix(path: Path, int_ids=(), str_ids=(), block_bytes: int = BLOCK_BYTES):
    """A dense cells x features table -> (X CSR float32, {id column: values}, feature names).

    `int_ids` / `str_ids`: the columns that identify a row (read as integers / text); every other
    column is a feature, except row-number columns left by a re-save. Read block by block, so
    the dense form of only one block is ever in memory."""
    import pyarrow as pa
    from scipy import sparse
    path = Path(path)
    names = header(path)
    int_ids = [c for c in int_ids if c in names]
    str_ids = [c for c in str_ids if c in names]
    ids = set(int_ids) | set(str_ids)
    feats = [n for n in names if n not in ids and not is_index_column(n)]
    types = {**{c: pa.int64() for c in int_ids}, **{c: pa.string() for c in str_ids},
             **{c: pa.float32() for c in feats}}
    if is_parquet(path):
        import pyarrow.parquet as pq
        it = pq.ParquetFile(str(path)).iter_batches(batch_size=20_000, columns=list(types))
    else:
        import pyarrow.csv as pcsv
        it = pcsv.open_csv(str(path),
                           read_options=pcsv.ReadOptions(column_names=names, skip_rows=1, block_size=block_bytes),
                           convert_options=pcsv.ConvertOptions(include_columns=list(types), column_types=types))
    parts, id_parts = [], {c: [] for c in ids}
    for b in it:
        n = b.num_rows
        if not n:
            continue
        where = {name: i for i, name in enumerate(b.schema.names)}
        a = np.empty((n, len(feats)), np.float32)
        for j, g in enumerate(feats):
            a[:, j] = b.column(where[g]).to_numpy(zero_copy_only=False)
        np.nan_to_num(a, copy=False)
        parts.append(sparse.csr_matrix(a))
        for c in int_ids:
            id_parts[c].append(b.column(where[c]).to_numpy(zero_copy_only=False).astype(np.int64))
        for c in str_ids:
            id_parts[c].append(["" if v is None else str(v) for v in b.column(where[c]).to_pylist()])
    X = sparse.vstack(parts, format="csr") if parts else sparse.csr_matrix((0, len(feats)), dtype=np.float32)
    out = {c: (np.concatenate(id_parts[c]) if id_parts[c] else np.zeros(0, np.int64)) for c in int_ids}
    out.update({c: [v for p in id_parts[c] for v in p] for c in str_ids})
    return X.astype(np.float32), out, feats


# ---------------------------------------------------------------------------- transcripts
def _codes(arr, index: dict) -> np.ndarray:
    """Feature ids of a string column through `index` (name -> id); unknown names -> -1."""
    import pyarrow.compute as pc
    d = pc.dictionary_encode(arr)
    if hasattr(d, "combine_chunks"):
        d = d.combine_chunks()
    lut = np.array([index.get(n, -1) for n in d.dictionary.to_pylist()] + [-1], np.int64)
    idx = pc.fill_null(d.indices, len(lut) - 1).to_numpy(zero_copy_only=False)
    return lut[idx]


def scan_transcripts(path: Path, gene_col: str, x_col: str, y_col: str, log=print, tag="") -> dict:
    """Pass 1 over a transcript table: every feature name with its row count, and the extent
    of the coordinates. Returns {"names": sorted names, "counts": {name: rows}, "n", "xlim", "ylim"}."""
    import pyarrow.compute as pc
    t0 = time.time()
    counts: dict[str, int] = {}
    n = 0
    lim = [np.inf, -np.inf, np.inf, -np.inf]
    for b in batches(path, floats=(x_col, y_col), strings=(gene_col,)):
        vc = pc.value_counts(b.column(b.schema.get_field_index(gene_col)))
        for name, c in zip(vc.field(0).to_pylist(), vc.field(1).to_pylist()):
            counts[name or ""] = counts.get(name or "", 0) + int(c)
        x, y = col(b, x_col), col(b, y_col)
        if len(x):
            lim = [min(lim[0], float(np.nanmin(x))), max(lim[1], float(np.nanmax(x))),
                   min(lim[2], float(np.nanmin(y))), max(lim[3], float(np.nanmax(y)))]
        n += b.num_rows
    log(f"[{tag}] transcripts, pass 1: {n:,} rows, {len(counts)} features  ({time.time() - t0:.0f}s)")
    return {"names": sorted(counts), "counts": counts, "n": n, "xlim": lim[:2], "ylim": lim[2:]}


def write_molecules(path: Path, gene_col: str, x_col: str, y_col: str, names, keep, to_px, image_size,
                    pixel_size_um: float, log=print, tag="", density=None):
    """Pass 2: stream a transcript table into a molecules.MoleculeWriter (returned, not yet
    finished). `to_px(x, y)` maps the table's coordinates to full-resolution image pixels;
    `density`, when given, also counts the kept molecules per background pixel."""
    import molecules
    t0 = time.time()
    index = {n: i for i, n in enumerate(names)}
    keep = np.asarray(keep, bool)
    w = molecules.MoleculeWriter(names, image_size, pixel_size_um, keep=keep, log=log)
    n = 0
    for i, b in enumerate(batches(path, floats=(x_col, y_col), strings=(gene_col,))):
        g = _codes(b.column(b.schema.get_field_index(gene_col)), index)
        px, py = to_px(col(b, x_col), col(b, y_col))
        ok = np.isfinite(px) & np.isfinite(py)
        g, px, py = g[ok], px[ok], py[ok]
        w.add(g, px * pixel_size_um, py * pixel_size_um)        # the writer takes µm in the image frame
        if density is not None:
            k = (g >= 0) & keep[np.maximum(g, 0)]
            density.add(px[k], py[k])
        n += b.num_rows
        if (i + 1) % 100 == 0:
            log(f"[{tag}]   {n:,} transcripts read  ({time.time() - t0:.0f}s)")
    w.rows_read = n
    return w


# ------------------------------------------------------------------------------ background
class ArraySource(Source):
    """Channels held in memory (C, H, W) as a pyramid source."""

    def __init__(self, names: list[str], data: np.ndarray):
        data = np.asarray(data)
        if data.ndim == 2:
            data = data[None]
        self.names = list(names)
        self._a = data
        self.height, self.width = int(data.shape[1]), int(data.shape[2])
        self.dtype = data.dtype

    def read_band(self, c, y0, y1):
        return self._a[c, y0:y1]


class Density:
    """Counts of points per pixel of a (width x height) image, for the neutral background."""

    def __init__(self, width: int, height: int):
        self.W, self.H = int(width), int(height)
        self.a = np.zeros(self.W * self.H, np.uint32)
        self.n = 0

    def add(self, x_px, y_px):
        ix = np.floor(np.asarray(x_px, np.float64)).astype(np.int64)
        iy = np.floor(np.asarray(y_px, np.float64)).astype(np.int64)
        ok = (ix >= 0) & (ix < self.W) & (iy >= 0) & (iy < self.H)
        u, c = np.unique(iy[ok] * self.W + ix[ok], return_counts=True)
        self.a[u] += c.astype(np.uint32)
        self.n += int(ok.sum())

    def image(self, sigma_px: float) -> np.ndarray:
        """uint8 (H, W): the counts smoothed with a Gaussian, 255 = the 99.9th percentile of
        the non-empty pixels."""
        from scipy import ndimage as ndi
        d = ndi.gaussian_filter(self.a.reshape(self.H, self.W).astype(np.float32), float(sigma_px))
        nz = d[d > 0]
        top = float(np.percentile(nz, 99.9)) if nz.size else 1.0
        return np.clip(d * (255.0 / max(top, 1e-12)), 0, 255).astype(np.uint8)


def background_grid(xlim, ylim, units_per_um: float, margin_um: float = 10.0, bg_um: float = BACKGROUND_UM):
    """Size of the neutral background over an extent given in source units.

    Returns (width, height, bg_um, x0, y0): pixels of `bg_um` µm each, covering the extent plus a
    margin, with (x0, y0) the source coordinates of the low edges. `bg_um` doubles until the
    image is under BACKGROUND_MAX_PX pixels."""
    w_um = (xlim[1] - xlim[0]) / units_per_um + 2 * margin_um
    h_um = (ylim[1] - ylim[0]) / units_per_um + 2 * margin_um
    if not (np.isfinite(w_um) and np.isfinite(h_um)) or w_um <= 0 or h_um <= 0:
        raise ValueError("no finite coordinates to build a background from")
    while (w_um / bg_um) * (h_um / bg_um) > BACKGROUND_MAX_PX:
        bg_um *= 2
    width, height = int(np.ceil(w_um / bg_um)), int(np.ceil(h_um / bg_um))
    return max(width, 2), max(height, 2), bg_um, xlim[0] - margin_um * units_per_um, ylim[0] - margin_um * units_per_um
