r"""AnnData files read with h5py (.h5ad) or zarr (AnnData-Zarr, also the tables/ of a SpatialData store).

No anndata, scanpy or pandas: the on-disk encoding is read directly. Used by stdformats.py.

    AnnFile(path)              an .h5ad file or an AnnData .zarr folder (zarr v2 with zarr 2; zarr v3 with zarrv3.py)
      .obs_names / .var_names  the dataframe indexes
      .column(df, name)        one obs / var column: Column(values) or Column(codes, categories)
      .columns(df)             the column names in their stored order
      .obsm(key)               an obsm array (n_obs, k) float64
      .matrix(where)           X, raw/X or layers/<name> as a scipy CSR (n_obs x n_var), float32
      .dict(path)              a uns group as plain Python (scalars, strings, lists, nested dicts)
      .image(path)             a uns image array (H, W[, 3 / 4])

What is understood (anndata's encoding-type attributes, anndata 0.7 to 0.11):
    dataframe (_index, column-order), categorical (codes + categories), string-array, array,
    nullable-integer / -boolean / -string-array (values + mask: masked entries become NA), csr_matrix,
    csc_matrix (also the older h5sparse_format attributes), dense arrays, dict, numeric-scalar, string;
    the anndata 0.7 categorical (a codes dataset whose "categories" attribute points at obs/__categories/).
    A zarr v3 store (anndata 0.11+ with zarr-python 3: zarr.json, strings as vlen-utf8) is read the same way.
Refused, with the reason: an .h5ad of anndata older than 0.7 (obs / var as one compound table), a zarr v3
array outside the subset zarrv3.py reads (the codec is named), anything that is not an AnnData file.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

import ngff

ROW_BLOCK = 20_000               # dense X is read in blocks of this many rows, each made sparse


def _s(x) -> str:
    return x.decode("utf-8", "replace") if isinstance(x, (bytes, bytearray, np.bytes_)) else str(x)


class Column:
    """One dataframe column: plain values (`values`), or categories with integer codes (-1 = NA)."""

    def __init__(self, values=None, codes=None, categories=None, ordered=False, encoding=""):
        self.values, self.codes, self.categories, self.ordered = values, codes, categories, ordered
        self.encoding = encoding

    @property
    def is_categorical(self) -> bool:
        return self.codes is not None

    @property
    def is_text(self) -> bool:
        return self.values is not None and self.values.dtype.kind in "OUS"

    def strings(self) -> list[str]:
        """Every row as text ("NA" for a missing value)."""
        if self.is_categorical:
            cats = [_s(c) for c in self.categories]
            return [cats[c] if c >= 0 else "NA" for c in self.codes.tolist()]
        v = self.values
        if v.dtype.kind == "f":
            return ["NA" if not np.isfinite(x) else (str(int(x)) if float(x).is_integer() else repr(float(x))) for x in v.tolist()]
        if v.dtype.kind == "b":
            return ["True" if x else "False" for x in v.tolist()]
        return [_s(x) for x in v.tolist()]

    def __len__(self):
        return len(self.codes) if self.is_categorical else len(self.values)


class AnnFile:
    def __init__(self, path, what: str = "anndata"):
        self.path, self.what = Path(path), what
        self._h5 = None
        if self.path.is_dir():
            ver = ngff.zarr_version(self.path)
            if ver is None:
                raise SystemExit(f"[{what}] not a zarr folder (no .zgroup / zarr.json): the AnnData-Zarr store is a "
                                 f"folder holding .zgroup (or zarr.json), obs/, var/ and X")
            try:
                self.root = ngff.open_group(self.path)
            except SystemExit as e:
                raise SystemExit(f"[{what}] {str(e).replace('[zarr v3] ', 'zarr v3: ', 1)}")
            self.kind = "zarr"
            self.zarr_format = ver
        else:
            import h5py
            try:
                self._h5 = h5py.File(str(self.path), "r")
            except OSError:
                raise SystemExit(f"[{what}] the file cannot be opened as HDF5: not an .h5ad file")
            self.root = self._h5
            self.kind = "h5ad"
        if "obs" not in self.root or "var" not in self.root:
            raise SystemExit(f"[{what}] no obs / var: not an AnnData file")
        for k in ("obs", "var"):
            if not self.is_group(self.root[k]):
                raise SystemExit(f"[{what}] {k} is one table (a compound array): written by anndata older than 0.7. "
                                 f"Open it with a current anndata / scanpy and save it again")
        self.obs_names = self.index(self.root["obs"])
        self.var_names = self.index(self.root["var"])
        self.n_obs, self.n_var = len(self.obs_names), len(self.var_names)

    def close(self):
        if self._h5 is not None:
            self._h5.close()
            self._h5 = None

    def __enter__(self):
        return self

    def __exit__(self, *a):
        self.close()

    # ---- nodes ---------------------------------------------------------------------------
    def is_group(self, n) -> bool:
        if self.kind == "zarr":
            return ngff.is_group(n)
        import h5py
        return isinstance(n, h5py.Group)

    def attrs(self, n) -> dict:
        out = {}
        for k, v in n.attrs.items():
            if isinstance(v, (bytes, np.bytes_)):
                v = _s(v)
            elif isinstance(v, np.ndarray):
                v = [_s(x) if isinstance(x, (bytes, np.bytes_)) else (x.item() if hasattr(x, "item") else x) for x in v.tolist()] \
                    if v.dtype.kind in "OSU" else v.tolist()
            elif isinstance(v, np.generic):
                v = v.item()
            out[k] = v
        return out

    def enc(self, n) -> str:
        return str(self.attrs(n).get("encoding-type", "") or "")

    def get(self, path: str):
        try:
            return self.root[path] if path in self.root else None
        except (KeyError, ValueError):
            return None

    def keys(self, g) -> list[str]:
        return list(g.keys()) if g is not None and self.is_group(g) else []

    def array(self, n) -> np.ndarray:
        """A dataset as numpy; text as an object array of str."""
        if self.kind == "h5ad":
            import h5py
            if h5py.check_string_dtype(n.dtype) is not None:
                v = n.asstr()[()]
                return np.asarray(v, dtype=object) if np.ndim(v) else v
            v = n[()]
        else:
            v = n[...] if n.shape else n[()]
        if isinstance(v, np.ndarray) and v.dtype.kind in "OSU":
            return np.array([_s(x) for x in v.ravel()], dtype=object).reshape(v.shape)
        if isinstance(v, (bytes, np.bytes_)):
            return _s(v)
        return v

    # ---- dataframes ----------------------------------------------------------------------
    def index(self, df) -> list[str]:
        a = self.attrs(df)
        name = a.get("_index", "_index")
        if name not in df:
            for alt in ("_index", "index", "obs_names", "var_names"):
                if alt in df:
                    name = alt
                    break
            else:
                raise SystemExit(f"[{self.what}] a dataframe without its index (_index)")
        return [_s(x) for x in np.asarray(self.array(df[name])).ravel().tolist()]

    def columns(self, df) -> list[str]:
        a = self.attrs(df)
        idx = a.get("_index", "_index")
        order = a.get("column-order")
        names = [str(c) for c in order] if isinstance(order, list) and order else \
            [k for k in self.keys(df) if k not in (idx, "__categories")]
        return [c for c in names if c in df and c != idx]

    def column(self, df, name: str) -> Column:
        n = df[name]
        e = self.enc(n)
        if self.is_group(n):
            if e == "categorical" or ("codes" in n and "categories" in n):
                codes = np.asarray(self.array(n["codes"])).astype(np.int64)
                cats = np.asarray(self.array(n["categories"]))
                return Column(codes=codes, categories=list(cats.tolist()), ordered=bool(self.attrs(n).get("ordered")),
                              encoding="categorical")
            if e.startswith("nullable") or ("values" in n and "mask" in n):
                vals = np.asarray(self.array(n["values"]))
                mask = np.asarray(self.array(n["mask"])).astype(bool)
                if vals.dtype.kind in "OSU":
                    vals = vals.astype(object)
                    vals[mask] = "NA"
                else:
                    vals = vals.astype(np.float64)
                    vals[mask] = np.nan
                return Column(values=vals, encoding=e or "nullable")
            raise ValueError(f"column of encoding {e!r}")
        a = self.attrs(n)
        cref = n.attrs.get("categories") if self.kind == "h5ad" else None
        if cref is not None and not isinstance(cref, (str, bytes, np.ndarray, list)):   # anndata 0.7: an object reference
            cats = np.asarray(self.array(self.root[cref]))
            return Column(codes=np.asarray(n[()]).astype(np.int64), categories=list(cats.tolist()), encoding="categorical-0.7")
        v = np.asarray(self.array(n))
        if v.ndim != 1:
            raise ValueError(f"a {v.ndim}-D column")
        return Column(values=v, encoding=a.get("encoding-type", "") or "array")

    # ---- arrays and matrices -------------------------------------------------------------
    def obsm_keys(self) -> list[str]:
        return self.keys(self.get("obsm"))

    def obsm(self, key: str) -> np.ndarray | None:
        g = self.get("obsm")
        if g is None or key not in g:
            return None
        n = g[key]
        if self.is_group(n):                           # a dataframe in obsm (spatialdata-io sometimes)
            cols = self.columns(n)
            return np.column_stack([np.asarray(self.column(n, c).values, np.float64) for c in cols]) if cols else None
        return np.asarray(self.array(n), np.float64)

    def matrix(self, where: str = "X"):
        """X, raw/X or layers/<name> as CSR float32 (n_obs x its n_var)."""
        from scipy import sparse
        n = self.get(where)
        if n is None:
            raise SystemExit(f"[{self.what}] no {where} in the AnnData file")
        if self.is_group(n):
            a = self.attrs(n)
            e = str(a.get("encoding-type") or a.get("h5sparse_format") or "")
            shape = a.get("shape") or a.get("h5sparse_shape")
            if not shape or not all(k in n for k in ("data", "indices", "indptr")):
                raise SystemExit(f"[{self.what}] {where}: a group that is not a sparse matrix (encoding {e!r})")
            shape = tuple(int(v) for v in shape)
            data = np.asarray(self.array(n["data"])).astype(np.float32, copy=False)
            indices = np.asarray(self.array(n["indices"]))
            indptr = np.asarray(self.array(n["indptr"]))
            if "csc" in e:
                return sparse.csc_matrix((data, indices, indptr), shape=shape).tocsr()
            return sparse.csr_matrix((data, indices, indptr), shape=shape)
        if len(n.shape) != 2:
            raise SystemExit(f"[{self.what}] {where}: a {len(n.shape)}-D array, not a matrix")
        blocks = []
        for r0 in range(0, n.shape[0], ROW_BLOCK):
            blocks.append(sparse.csr_matrix(np.asarray(n[r0:r0 + ROW_BLOCK], dtype=np.float32)))
        return sparse.vstack(blocks).tocsr() if blocks else sparse.csr_matrix((0, n.shape[1]), dtype=np.float32)

    def dict(self, path: str, depth: int = 6):
        """A uns entry as plain Python (dicts of scalars, strings and small lists)."""
        n = self.get(path)
        return None if n is None else self._plain(n, depth)

    def _plain(self, n, depth):
        if self.is_group(n):
            if depth <= 0:
                return {}
            return {k: self._plain(n[k], depth - 1) for k in self.keys(n)}
        v = self.array(n)
        if isinstance(v, np.ndarray):
            if v.ndim == 0:
                return v.item()
            if v.size <= 64:
                return [x.item() if hasattr(x, "item") else x for x in v.ravel().tolist()] if v.ndim == 1 else v
            return v
        return v.item() if hasattr(v, "item") else v

    def image(self, path: str) -> np.ndarray | None:
        n = self.get(path)
        if n is None or self.is_group(n):
            return None
        return np.asarray(self.array(n))


def looks_like_counts(X, sample: int = 200_000) -> bool:
    """Non-negative whole numbers (the first `sample` stored values)."""
    d = np.asarray(X.data[:sample], np.float64)
    return bool(d.size == 0 or (np.all(d >= 0) and np.all(np.abs(d - np.round(d)) < 1e-6)))
