r"""A small, self-contained reader for zarr v3 stores (zarr.json), the subset OME-NGFF 0.5 images and
AnnData / SpatialData stores written with zarr-python 3 use. Read-only, local folders only.

SODA's own zarr is held at 2.x (its image pyramids are written with it), and zarr-python 2 cannot read
the v3 format, so the v3 layout is read here directly, to the spec
(https://zarr-specs.readthedocs.io/en/latest/v3/core/index.html):

    open_node(path)            -> Group | Array (by zarr.json's node_type)
    Group                      .attrs (dict), .keys(), "a/b" in g, g["a/b"]  (the zarr 2 API subset
                                annstore.py and ngff.py use)
    Array                      .shape, .dtype, .ndim, .chunks, .fill_value, .attrs, a[ints / slices / ...]
    check_store(path)          the first array below `path` SODA cannot read, with the reason (or None)

What is read:
    metadata         zarr.json of groups and arrays; attributes (NGFF 0.5 puts its metadata under
                     attributes["ome"]); dimension_names
    chunk grid       regular
    chunk keys       "default" (c/0/1, separator "/" or ".") and "v2" (0.1, separator "." or "/")
    data types       bool, int8 .. int64, uint8 .. uint64, float16 / 32 / 64, complex64 / 128,
                     variable-length strings ("string", vlen-utf8) and bytes (vlen-bytes),
                     fixed_length_utf32
    codecs           transpose (array -> array); bytes (endianness), vlen-utf8, vlen-bytes and
                     sharding_indexed (index at the end or the start; inner chunks read one by one
                     with a seek, never the whole shard) (array -> bytes); blosc, gzip, zstd (through
                     numcodecs) and crc32c (checked on shard indexes) (bytes -> bytes)
    fill_value       numbers, "NaN" / "Infinity" / "-Infinity", "0x..." bit patterns, booleans, "";
                     a chunk that is not stored (no file, or an empty shard index entry) is the fill value
Refused (V3Error, a SystemExit, with the reason and the codec / data type named): any other codec,
data type, chunk grid, chunk key encoding or a storage transformer; a negative slice step.
"""
from __future__ import annotations

import itertools
import json
import math
import struct
import threading
from pathlib import Path

import numpy as np

CODECS = ("bytes", "transpose", "blosc", "gzip", "zstd", "crc32c", "sharding_indexed", "vlen-utf8", "vlen-bytes")
_NUM = {"bool": "?", "int8": "i1", "int16": "i2", "int32": "i4", "int64": "i8", "uint8": "u1", "uint16": "u2",
        "uint32": "u4", "uint64": "u8", "float16": "f2", "float32": "f4", "float64": "f8",
        "complex64": "c8", "complex128": "c16"}
_STR = ("string", "variable_length_utf8", "str")
_BYTES = ("bytes", "variable_length_bytes")
EMPTY = 2 ** 64 - 1                      # a shard index entry for a chunk that is not stored


class V3Error(SystemExit):
    """A zarr v3 store (or part of one) outside what this reader reads; the message says what."""


def _refuse(msg: str):
    raise V3Error(f"[zarr v3] {msg}. SODA reads zarr v3 stores with the codecs {', '.join(CODECS)}; write the store "
                  f"again without it (or as zarr v2: anndata.settings.zarr_write_format = 2 before write_zarr; "
                  f"OME-NGFF 0.4 for images), or give the data as .h5ad / OME-TIFF")


def read_meta(folder) -> dict | None:
    p = Path(folder) / "zarr.json"
    if not p.exists():
        return None
    try:
        m = json.loads(p.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        raise V3Error(f"[zarr v3] a zarr.json that is not JSON")
    if not isinstance(m, dict) or m.get("zarr_format") != 3:
        raise V3Error(f"[zarr v3] a zarr.json without zarr_format 3")
    return m


def attrs_of(folder) -> dict:
    """The attributes of a v3 node ({} when there is none)."""
    try:
        m = read_meta(folder)
    except V3Error:
        return {}
    return dict((m or {}).get("attributes") or {})


# ---------------------------------------------------------------------------------------
# data types and fill values
# ---------------------------------------------------------------------------------------
def _dtype(dt) -> tuple[np.dtype, str]:
    """(numpy dtype, kind): kind "num", "str" (variable-length text), "bytes" (variable-length bytes), "U" (fixed utf32)."""
    name, conf = (dt, {}) if isinstance(dt, str) else (str(dt.get("name")), dt.get("configuration") or {})
    if name in _NUM:
        return np.dtype(_NUM[name]), "num"
    if name in _STR:
        return np.dtype(object), "str"
    if name in _BYTES:
        return np.dtype(object), "bytes"
    if name == "fixed_length_utf32":
        n = int(conf.get("length_bytes", 4)) // 4
        return np.dtype(f"<U{max(n, 1)}"), "U"
    _refuse(f"the data type {name!r} is not read")


def _fill(fv, dtype: np.dtype, kind: str):
    if kind in ("str", "bytes"):
        return "" if kind == "str" and fv in (None, 0) else (fv if fv is not None else b"")
    if kind == "U":
        return fv or ""
    if fv is None:
        return np.zeros((), dtype)[()]
    if dtype.kind == "b":
        return bool(fv)
    if dtype.kind == "c" and isinstance(fv, (list, tuple)) and len(fv) == 2:
        return complex(_fill(fv[0], np.dtype("f8"), "num"), _fill(fv[1], np.dtype("f8"), "num"))
    if isinstance(fv, str):
        if fv == "NaN":
            return np.nan
        if fv == "Infinity":
            return np.inf
        if fv == "-Infinity":
            return -np.inf
        if fv.startswith("0x"):
            bits = int(fv, 16)
            size = dtype.itemsize if dtype.kind != "c" else dtype.itemsize // 2
            raw = bits.to_bytes(size, "big")
            ft = {2: ">f2", 4: ">f4", 8: ">f8"}.get(size)
            if dtype.kind in "fc" and ft:
                return float(np.frombuffer(raw, ft)[0])
            return int(np.frombuffer(raw, np.dtype(dtype).newbyteorder(">"))[0])
        _refuse(f"the fill value {fv!r} is not understood")
    return np.asarray(fv).astype(dtype)[()] if dtype.kind != "b" else bool(fv)


# ---------------------------------------------------------------------------------------
# crc32c (Castagnoli), for the shard indexes
# ---------------------------------------------------------------------------------------
def _crc_table():
    t = []
    for i in range(256):
        c = i
        for _ in range(8):
            c = (c >> 1) ^ 0x82F63B78 if c & 1 else c >> 1
        t.append(c)
    return t


_CRC = _crc_table()


def crc32c(buf: bytes) -> int:
    c = 0xFFFFFFFF
    t = _CRC
    for b in buf:
        c = t[(c ^ b) & 0xFF] ^ (c >> 8)
    return c ^ 0xFFFFFFFF


# ---------------------------------------------------------------------------------------
# codecs
# ---------------------------------------------------------------------------------------
def _codec(c) -> tuple[str, dict]:
    if isinstance(c, str):
        return c, {}
    return str(c.get("name")), dict(c.get("configuration") or {})


class Pipeline:
    """A codec chain: array -> array codecs, one array -> bytes codec, bytes -> bytes codecs."""

    def __init__(self, codecs: list, dtype: np.dtype, kind: str, fill, chunk_shape: tuple):
        self.dtype, self.kind, self.fill = dtype, kind, fill
        self.aa, self.ab, self.bb = [], None, []
        for c in codecs or []:
            name, conf = _codec(c)
            if name == "transpose":
                if self.ab is not None:
                    _refuse("a transpose codec after the array -> bytes codec")
                order = conf.get("order")
                nd = len(chunk_shape)
                if order == "C" or order is None:
                    order = list(range(nd))
                elif order == "F":
                    order = list(range(nd))[::-1]
                order = [int(v) for v in order]
                if sorted(order) != list(range(nd)):
                    _refuse(f"a transpose order {order} for {nd} dimensions")
                self.aa.append(order)
            elif name in ("bytes", "endian", "vlen-utf8", "vlen-bytes", "sharding_indexed"):
                if self.ab is not None:
                    _refuse("two array -> bytes codecs")
                self.ab = (name, conf)
            elif name in ("blosc", "gzip", "zstd", "crc32c"):
                if self.ab is None:
                    _refuse(f"the codec {name!r} before the array -> bytes codec")
                self.bb.append((name, conf))
            else:
                _refuse(f"the codec {name!r} is not read")
        if self.ab is None:
            _refuse("no array -> bytes codec")
        name, conf = self.ab
        if name == "vlen-utf8" and kind != "str" or name == "vlen-bytes" and kind != "bytes":
            _refuse(f"the codec {name!r} for a data type it does not encode")
        if name in ("bytes", "endian") and kind in ("str", "bytes"):
            _refuse("variable-length strings with the bytes codec")
        self.shape_in = tuple(int(v) for v in chunk_shape)
        s = self.shape_in
        for order in self.aa:
            s = tuple(s[i] for i in order)
        self.shape_ab = s                       # the shape the array -> bytes codec sees
        self.shard = None
        if name == "sharding_indexed":
            self.shard = Shard(conf, dtype, kind, fill, self.shape_ab)
        self._nc = {}

    def _numcodec(self, name, conf):
        c = self._nc.get(name)
        if c is None:
            import numcodecs
            if name == "blosc":
                c = numcodecs.Blosc()
            elif name == "gzip":
                c = numcodecs.GZip(level=int(conf.get("level", 5)))
            elif name == "zstd":
                c = numcodecs.Zstd()
            self._nc[name] = c
        return c

    def bytes_decode(self, buf: bytes) -> bytes:
        for name, conf in reversed(self.bb):
            if name == "crc32c":
                if len(buf) < 4:
                    raise V3Error("[zarr v3] a crc32c chunk shorter than its checksum")
                if len(buf) <= (1 << 16) and crc32c(buf[:-4]) != struct.unpack("<I", buf[-4:])[0]:
                    raise V3Error("[zarr v3] a chunk whose crc32c checksum does not match (a damaged file)")
                buf = buf[:-4]
            elif name == "gzip":
                import gzip
                import zlib
                try:
                    buf = gzip.decompress(buf)
                except OSError:
                    buf = zlib.decompress(buf)
            else:
                buf = bytes(self._numcodec(name, conf).decode(buf))
        return buf

    def decode(self, buf: bytes) -> np.ndarray:
        """Stored bytes of one chunk -> the chunk (chunk_shape, dtype)."""
        buf = self.bytes_decode(buf)
        name, conf = self.ab
        n = int(np.prod(self.shape_ab)) if self.shape_ab else 1
        if name == "sharding_indexed":
            a = self.shard.decode_all(buf)
        elif name == "vlen-utf8" or name == "vlen-bytes":
            a = _vlen(buf, name == "vlen-utf8", n).reshape(self.shape_ab)
        else:
            if self.dtype.itemsize > 1 or self.dtype.kind in "U":
                end = conf.get("endian", "little")
                dt = self.dtype.newbyteorder("<" if end == "little" else ">")
            else:
                dt = self.dtype
            if len(buf) != n * dt.itemsize:
                raise V3Error(f"[zarr v3] a chunk of {len(buf)} bytes for {n} values of {dt.itemsize} bytes")
            a = np.frombuffer(buf, dt).reshape(self.shape_ab).astype(self.dtype.newbyteorder("="), copy=False)
        for order in reversed(self.aa):
            a = a.transpose(np.argsort(order))
        return a


def _vlen(buf: bytes, text: bool, n: int) -> np.ndarray:
    """numcodecs' VLenUTF8 / VLenBytes layout: u32 count, then u32 length + the bytes of each item."""
    mv = memoryview(buf)
    if len(buf) < 4:
        raise V3Error("[zarr v3] a variable-length chunk without its header")
    cnt = struct.unpack_from("<I", mv, 0)[0]
    if cnt != n:
        raise V3Error(f"[zarr v3] a variable-length chunk of {cnt} items for {n}")
    out = np.empty(n, object)
    off = 4
    for i in range(n):
        ln = struct.unpack_from("<I", mv, off)[0]
        off += 4
        b = bytes(mv[off:off + ln])
        off += ln
        out[i] = b.decode("utf-8", "replace") if text else b
    return out


class Shard:
    """sharding_indexed: a shard is inner chunks + an index of (offset, nbytes) per inner chunk."""

    def __init__(self, conf: dict, dtype, kind, fill, shard_shape: tuple):
        self.inner_shape = tuple(int(v) for v in conf.get("chunk_shape") or [])
        if len(self.inner_shape) != len(shard_shape) or any(s <= 0 or ss % s for s, ss in zip(self.inner_shape, shard_shape)):
            _refuse(f"a sharding_indexed inner chunk shape {list(self.inner_shape)} that does not divide the shard "
                    f"{list(shard_shape)}")
        self.shard_shape = shard_shape
        self.per = tuple(ss // s for s, ss in zip(self.inner_shape, shard_shape))
        self.n = int(np.prod(self.per)) if self.per else 1
        self.inner = Pipeline(conf.get("codecs") or [{"name": "bytes", "configuration": {"endian": "little"}}],
                              dtype, kind, fill, self.inner_shape)
        self.dtype, self.kind, self.fill = dtype, kind, fill
        self.loc = str(conf.get("index_location", "end"))
        if self.loc not in ("end", "start"):
            _refuse(f"a shard index location {self.loc!r}")
        self.icodecs = [_codec(c) for c in conf.get("index_codecs") or [{"name": "bytes"}, {"name": "crc32c"}]]
        for name, c in self.icodecs:
            if name not in ("bytes", "crc32c"):
                _refuse(f"the shard index codec {name!r} is not read")
        self.icrc = any(name == "crc32c" for name, _ in self.icodecs)
        ab = next((c for name, c in self.icodecs if name == "bytes"), {})
        self.iend = "<" if ab.get("endian", "little") == "little" else ">"
        self.isize = 16 * self.n + (4 if self.icrc else 0)

    def index(self, raw: bytes) -> np.ndarray:
        if len(raw) != self.isize:
            raise V3Error("[zarr v3] a shard index of the wrong size (a damaged file)")
        if self.icrc:
            if crc32c(raw[:-4]) != struct.unpack("<I", raw[-4:])[0]:
                raise V3Error("[zarr v3] a shard index whose crc32c checksum does not match (a damaged file)")
            raw = raw[:-4]
        return np.frombuffer(raw, self.iend + "u8").reshape(self.per + (2,)).astype(np.uint64)

    def index_of_file(self, fh, size: int) -> np.ndarray:
        if size < self.isize:
            raise V3Error("[zarr v3] a shard file shorter than its index (a damaged file)")
        fh.seek(0 if self.loc == "start" else size - self.isize)
        return self.index(fh.read(self.isize))

    def inner_chunk(self, entry, read) -> np.ndarray | None:
        off, nb = int(entry[0]), int(entry[1])
        if off == EMPTY and nb == EMPTY:
            return None
        return self.inner.decode(read(off, nb))

    def decode_all(self, buf: bytes) -> np.ndarray:
        """A whole shard from its bytes (used when bytes -> bytes codecs wrap the shard)."""
        raw = buf[:self.isize] if self.loc == "start" else buf[len(buf) - self.isize:]
        idx = self.index(raw)
        out = _full(self.shard_shape, self.dtype, self.kind, self.fill)
        for pos in itertools.product(*[range(p) for p in self.per]):
            ch = self.inner_chunk(idx[pos], lambda o, n: buf[o:o + n])
            if ch is not None:
                sl = tuple(slice(p * s, (p + 1) * s) for p, s in zip(pos, self.inner_shape))
                out[sl] = ch
        return out


def _full(shape, dtype, kind, fill) -> np.ndarray:
    if kind in ("str", "bytes"):
        a = np.empty(shape, object)
        a.fill(fill)
        return a
    return np.full(shape, fill, dtype)


# ---------------------------------------------------------------------------------------
# nodes
# ---------------------------------------------------------------------------------------
class Group:
    def __init__(self, path: Path, meta: dict):
        self.path, self.meta = Path(path), meta
        self.attrs = dict(meta.get("attributes") or {})

    def keys(self) -> list[str]:
        return sorted(d.name for d in self.path.iterdir() if d.is_dir() and (d / "zarr.json").exists())

    def __iter__(self):
        return iter(self.keys())

    def _resolve(self, key: str) -> Path | None:
        parts = [p for p in str(key).replace("\\", "/").split("/") if p not in ("", ".")]
        if any(p == ".." for p in parts):
            return None
        q = self.path.joinpath(*parts) if parts else self.path
        return q if (q / "zarr.json").exists() else None

    def __contains__(self, key) -> bool:
        return self._resolve(key) is not None

    def __getitem__(self, key):
        q = self._resolve(key)
        if q is None:
            raise KeyError(key)
        return open_node(q)


class Array:
    def __init__(self, path: Path, meta: dict):
        self.path, self.meta = Path(path), meta
        self.attrs = dict(meta.get("attributes") or {})
        if meta.get("storage_transformers"):
            _refuse(f"a storage transformer ({_codec(meta['storage_transformers'][0])[0]!r})")
        self.shape = tuple(int(v) for v in meta.get("shape") or [])
        self.ndim = len(self.shape)
        self.dtype, self.kind = _dtype(meta.get("data_type"))
        grid = meta.get("chunk_grid") or {}
        gname, gconf = _codec(grid)
        if gname != "regular":
            _refuse(f"the chunk grid {gname!r} is not read")
        self.chunk_shape = tuple(int(v) for v in gconf.get("chunk_shape") or [])
        if len(self.chunk_shape) != self.ndim or any(c <= 0 for c in self.chunk_shape):
            raise V3Error(f"[zarr v3] a chunk shape {list(self.chunk_shape)} for a {self.ndim}-D array")
        kname, kconf = _codec(meta.get("chunk_key_encoding") or {"name": "default"})
        if kname not in ("default", "v2"):
            _refuse(f"the chunk key encoding {kname!r} is not read")
        self.key_v2 = kname == "v2"
        self.sep = str(kconf.get("separator", "." if self.key_v2 else "/"))
        if self.sep not in ("/", "."):
            _refuse(f"the chunk key separator {self.sep!r}")
        self.fill_value = _fill(meta.get("fill_value"), self.dtype, self.kind)
        self.pipe = Pipeline(meta.get("codecs") or [], self.dtype, self.kind, self.fill_value, self.chunk_shape)
        self.dimension_names = meta.get("dimension_names")
        # what one read returns: the inner chunk when the shards can be read piece by piece
        sh = self.pipe.shard
        self._partial = sh is not None and not self.pipe.aa and not self.pipe.bb
        self.chunks = sh.inner_shape if self._partial else self.chunk_shape
        self._lock = threading.Lock()
        self._idx_cache: dict = {}
        self._last = None

    def __len__(self):
        if not self.shape:
            raise TypeError("len() of a 0-d array")
        return self.shape[0]

    @property
    def size(self):
        return int(np.prod(self.shape)) if self.shape else 1

    def _key(self, coords) -> Path:
        if self.key_v2:
            rel = self.sep.join(str(c) for c in coords) if coords else "0"
        else:
            rel = self.sep.join(["c"] + [str(c) for c in coords])
        return self.path.joinpath(*rel.split("/"))

    def _stored_chunk(self, coords) -> np.ndarray | None:
        """A whole stored chunk (shard), or None (not stored)."""
        f = self._key(coords)
        if not f.exists():
            return None
        last = self._last
        if last is not None and last[0] == coords:
            return last[1]
        a = self.pipe.decode(f.read_bytes())
        self._last = (coords, a)
        return a

    def _inner(self, coords) -> np.ndarray | None:
        """One inner chunk of a sharded array, read with a seek (the shard's index is cached)."""
        sh = self.pipe.shard
        scoords = tuple(c // p for c, p in zip(coords, sh.per))
        pos = tuple(c % p for c, p in zip(coords, sh.per))
        f = self._key(scoords)
        if not f.exists():
            return None
        with open(f, "rb") as fh:
            with self._lock:
                idx = self._idx_cache.get(scoords)
            if idx is None:
                idx = sh.index_of_file(fh, f.stat().st_size)
                with self._lock:
                    if len(self._idx_cache) > 4096:
                        self._idx_cache.clear()
                    self._idx_cache[scoords] = idx

            def read(o, n):
                fh.seek(o)
                b = fh.read(n)
                if len(b) != n:
                    raise V3Error("[zarr v3] a shard shorter than its index says (a damaged file)")
                return b
            return sh.inner_chunk(idx[pos], read)

    def _chunk(self, coords):
        return self._inner(coords) if self._partial else self._stored_chunk(coords)

    def _norm(self, key):
        if not isinstance(key, tuple):
            key = (key,)
        if any(k is Ellipsis for k in key):
            i = next(j for j, k in enumerate(key) if k is Ellipsis)
            fill = self.ndim - (len(key) - 1)
            key = key[:i] + (slice(None),) * fill + key[i + 1:]
        if len(key) > self.ndim:
            raise IndexError(f"{len(key)} indices for a {self.ndim}-D array")
        key = key + (slice(None),) * (self.ndim - len(key))
        out = []
        for k, n in zip(key, self.shape):
            if isinstance(k, (int, np.integer)):
                i = int(k) + (n if k < 0 else 0)
                if not 0 <= i < n:
                    raise IndexError(f"index {int(k)} out of range for {n}")
                out.append((i, i + 1, 1, True))
            elif isinstance(k, slice):
                a, b, s = k.indices(n)
                if s <= 0:
                    raise V3Error("[zarr v3] a negative slice step is not read")
                out.append((a, max(a, b), s, False))
            else:
                raise IndexError(f"index {k!r}: only integers and slices")
        return out

    def __getitem__(self, key):
        if not self.ndim:
            ch = self._chunk(())
            v = _full((), self.dtype, self.kind, self.fill_value) if ch is None else np.asarray(ch).reshape(())
            return v[()]
        sel = self._norm(key)
        box = [b - a for a, b, _, _ in sel]
        out = _full(tuple(box), self.dtype, self.kind, self.fill_value)
        if all(n > 0 for n in box):
            cs = self.chunks
            ranges = [range(a // c, (b - 1) // c + 1) for (a, b, _, _), c in zip(sel, cs)]
            for coords in itertools.product(*ranges):
                ch = self._chunk(tuple(coords))
                if ch is None:
                    continue
                src, dst = [], []
                for (a, b, _, _), c, k in zip(sel, cs, coords):
                    c0 = k * c
                    lo, hi = max(a, c0), min(b, c0 + c)
                    src.append(slice(lo - c0, hi - c0))
                    dst.append(slice(lo - a, hi - a))
                out[tuple(dst)] = ch[tuple(src)]
        step = tuple(slice(None, None, s) if s != 1 else slice(None) for _, _, s, _ in sel)
        if any(s != 1 for _, _, s, _ in sel):
            out = out[step]
        drop = tuple(0 if d else slice(None) for *_, d in sel)
        if any(d for *_, d in sel):
            out = out[drop]
        return out


def open_node(path):
    """A v3 Group or Array (refused with the reason when it is outside the subset read here)."""
    p = Path(path)
    m = read_meta(p)
    if m is None:
        raise V3Error("[zarr v3] no zarr.json here")
    t = m.get("node_type")
    if t == "group":
        return Group(p, m)
    if t == "array":
        return Array(p, m)
    raise V3Error(f"[zarr v3] a node of type {t!r}")


def open_array(path) -> Array:
    n = open_node(path)
    if not isinstance(n, Array):
        raise V3Error("[zarr v3] a group where an array was expected")
    return n


def open_group(path) -> Group:
    n = open_node(path)
    if not isinstance(n, Group):
        raise V3Error("[zarr v3] an array where a group was expected")
    return n


def check_store(path, limit: int = 5000) -> str | None:
    """The reason the first array below `path` cannot be read (its codec, data type ...), or None. Reads the
    zarr.json files only (at most `limit` of them)."""
    seen = 0
    stack = [Path(path)]
    while stack and seen < limit:
        d = stack.pop()
        try:
            m = read_meta(d)
        except V3Error as e:
            return str(e)
        if m is None:
            continue
        seen += 1
        if m.get("node_type") == "array":
            try:
                Array(d, m)
            except V3Error as e:
                return str(e)
            continue
        try:
            stack += sorted((c for c in d.iterdir() if c.is_dir() and (c / "zarr.json").exists()), reverse=True)
        except OSError:
            continue
    return None
