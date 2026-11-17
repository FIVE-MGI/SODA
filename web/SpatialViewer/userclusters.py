r"""Bring your own UMAP and clusters: a lab's Seurat / scanpy results on a spatial expression layer.

SODA computes its own whole-section UMAP and Leiden clusters for every spatial layer
(sectionumap.py). A lab that already clustered the same bins / spots / cells, and named the
clusters, can attach those results to one layer instead of (or next to) SODA's:

    --user-layer 008um                       the layer the results belong to (required)
    --user-umap  umap.csv                    barcode, then two numeric columns (UMAP_1, UMAP_2 ...)
    --user-obs   clusters.csv  (repeatable)  barcode, then one column per clustering / annotation
    --user-h5ad  analysis.h5ad               obs columns + an obsm embedding (+ uns/<col>_colors)
    --user-embedding X_umap                  the obsm key (h5ad only)
    --user-columns seurat_clusters,cell_type  which columns become clusterings (default: every usable one)
    --user-barcode-prefix S1_ / --user-barcode-suffix _1
                                             keep only the barcodes of one sample of a merged object
    --user-replace                           drop SODA's own clusterings (default: keep both, user's first)

The results are written into the layer's umap/ folder in the format the viewer already reads
(docs/DATA_FORMAT.md, "umap/"): each user clustering becomes a clusters_user<N>.bin (uint16 per
positions.bin row, 0xFFFF = not assigned) and an index.json entry with "id", "name",
"source": "user", the category names ("categories") and, when given, their colours ("colors").
A user embedding replaces coords.bin; without one SODA's coordinates stay. SODA's own index.json
and coords.bin are kept as index_soda.json / coords_soda.bin, so a later export without user
files puts SODA's back (restore) and the suggested cell types (celltypes.py) are always scored on
SODA's population of bins. Markers per user cluster come from sectionumap._markers on the same
log1p CP10k values, so the legends, suggested names, volcano group pickers and QC all work.

Rows are matched by barcode. The prepare cache keeps each layer's barcodes in the zarr
(expr/<key>.zarr/barcodes, in xy order, written by expr.build_expr); for caches made before
that, they are read again from the Space Ranger outs named in the cache manifest. Any layer
that keeps barcodes this way (Visium spots, Visium HD bins, cell layers) works the same.
Matching: exact first, then the usual Seurat / scanpy variants: a different or missing
"-1"-style suffix ("-1", "_1", "-1-0", ".1") and a sample prefix ("S1_", "slice1_", "S1:").
Below MIN_MATCH of the file's barcodes found, the import stops; a partial match is a warning.

De-identification: column names and category values go through deid.py exactly as
`prepare.py rnaseq` does (a sample-identity column has every value replaced by a number);
index.json is checked for leftover IDs and barcodes before it is written. No barcode and no
input file name is ever written to the export.
"""
from __future__ import annotations

import csv
import gzip
import json
import re
import time
from pathlib import Path

import numpy as np
from scipy import sparse

import deid
import sectionumap
from hvgchunks import _BARCODE

NO_CLUSTER = 0xFFFF
MIN_MATCH = 0.5               # below this share of the file's barcodes found: stop
MAX_CATEGORIES = 200          # as rnaseq.py: more distinct values is an id column, not a clustering
MAX_INT_CATEGORIES = 30       # integer columns with at most this many values are categories ...
# ... and up to MAX_CATEGORIES when the name says clustering (seurat_clusters, leiden, SCT_snn_res.0.8) or the
# column was asked for with --user-columns: a fine Seurat clustering easily has more than 30 clusters
_CLUSTER_NAME = re.compile(r"cluster|leiden|louvain|seurat|ident|snn_res|_res\.|^res\.|celltype|cell_type|annot", re.I)
SODA_INDEX, SODA_COORDS = "index_soda.json", "coords_soda.bin"
_NA = {"", "na", "nan", "none", "null", "<na>"}
_SUFFIX = re.compile(r"(?:[-_.]\d{1,3})+$")          # -1, _1, -1-0, -1_2, .1
_SEPS = re.compile(r"[_:#|]")
_HEX = re.compile(r"^#?([0-9a-fA-F]{6})(?:[0-9a-fA-F]{2})?$")
_COLOR_COL = re.compile(r"^(.*?)[._]colou?rs?$", re.I)
_BC_HEADER = {"barcode", "barcodes", "cell", "cells", "cell_id", "cellid", "spot", "spot_id", "bin", "bin_id"}


# ---------------------------------------------------------------------------------------
# command line (shared by `prepare.py export` and `prepare.py userclusters`)
# ---------------------------------------------------------------------------------------
def add_args(p):
    g = p.add_argument_group("your own UMAP / clusters (Seurat, scanpy) for one spatial layer (userclusters.py)")
    g.add_argument("--user-layer", help="layer key the files below belong to, e.g. 008um, 016um or spots")
    g.add_argument("--user-umap", help="CSV: barcode, then two numeric columns (the embedding)")
    g.add_argument("--user-obs", action="append", default=[],
                   help="CSV: barcode, then one column per clustering / annotation (repeatable)")
    g.add_argument("--user-h5ad", help="AnnData with the clusterings in obs and the embedding in obsm")
    g.add_argument("--user-embedding", default=None, help="obsm key of the embedding in --user-h5ad (default X_umap)")
    g.add_argument("--user-columns", help="comma list of the columns to use (default: every categorical column)")
    g.add_argument("--user-barcode-prefix", help="keep only barcodes starting with this (one sample of a merged object)")
    g.add_argument("--user-barcode-suffix", help="keep only barcodes ending with this (e.g. _2 from Seurat merge)")
    g.add_argument("--user-replace", action="store_true",
                   help="leave SODA's own Leiden clusterings out (default: keep both, yours shown first)")


def specs_from_args(args) -> list[dict]:
    """[spec] from the --user-* options, or [] when none was given."""
    files = args.user_umap or args.user_obs or args.user_h5ad
    if not files:
        if args.user_layer:
            raise SystemExit("--user-layer given without --user-umap, --user-obs or --user-h5ad")
        return []
    if not args.user_layer:
        raise SystemExit("--user-umap / --user-obs / --user-h5ad need --user-layer (e.g. 008um): "
                         "say which layer the barcodes belong to")
    return [make_spec(layer=args.user_layer, umap=args.user_umap, obs=args.user_obs, h5ad=args.user_h5ad,
                      embedding=args.user_embedding,
                      columns=args.user_columns.split(",") if args.user_columns else None,
                      prefix=args.user_barcode_prefix, suffix=args.user_barcode_suffix,
                      replace=args.user_replace)]


def make_spec(layer: str, umap=None, obs=(), h5ad=None, embedding=None, columns=None, prefix=None,
              suffix=None, replace=False) -> dict:
    if isinstance(obs, (str, Path)):
        obs = [obs]
    spec = {"layer": str(layer), "umap": Path(umap) if umap else None, "obs": [Path(p) for p in obs or []],
            "h5ad": Path(h5ad) if h5ad else None, "embedding": embedding,
            "columns": [c.strip() for c in columns if c.strip()] if columns else None,
            "prefix": prefix or None, "suffix": suffix or None, "replace": bool(replace)}
    if spec["h5ad"] and (spec["umap"] or spec["obs"]):
        raise SystemExit("give either --user-h5ad or the CSV files (--user-umap / --user-obs), not both")
    if embedding and not spec["h5ad"]:
        raise SystemExit("--user-embedding names an obsm key: it needs --user-h5ad")
    for p in [spec["umap"], spec["h5ad"], *spec["obs"]]:
        if p is not None and not p.exists():
            raise SystemExit(f"[userclusters] file not found: {p}")
    return spec


def check_layers(specs: list[dict], manifest: dict):
    """Fail before any work when a spec names a layer the dataset does not have."""
    keys = [L["key"] for L in manifest.get("expr", [])]
    seen = set()
    for s in specs or []:
        if s["layer"] not in keys:
            raise SystemExit(f"[userclusters] --user-layer {s['layer']}: this dataset has the layers {keys}")
        if s["layer"] in seen:
            raise SystemExit(f"[userclusters] layer {s['layer']} given twice")
        seen.add(s["layer"])


# ---------------------------------------------------------------------------------------
# reading the user's files -> tables {barcodes, embedding, columns{name: (values, order)}, colors}
# ---------------------------------------------------------------------------------------
def _open_text(p: Path):
    return gzip.open(p, "rt", encoding="utf-8-sig", newline="") if p.suffix == ".gz" \
        else open(p, encoding="utf-8-sig", newline="")


def read_table(path: Path) -> tuple[list[str], list[str], list[list[str]]]:
    """(barcodes, column names, columns as lists of strings) of a CSV / TSV (.gz ok).

    The barcode column is the one headed barcode / cell / cell_id ..., else the first column
    (R's write.csv and pandas' to_csv put the row names there, often with an empty header).
    """
    with _open_text(path) as fh:
        first = fh.readline()
        delim = "\t" if first.count("\t") > first.count(",") else ","
        rows = list(csv.reader([first], delimiter=delim)) + list(csv.reader(fh, delimiter=delim))
    rows = [r for r in rows if r and any(x.strip() for x in r)]
    if len(rows) < 2:
        raise SystemExit(f"[userclusters] {path.name}: no data rows")
    header = [h.strip() for h in rows[0]]
    body = rows[1:]
    # R's write.csv(row.names = TRUE) writes one header cell fewer than data cells
    if all(len(r) == len(header) + 1 for r in body[:100]):
        header = [""] + header
    bi = next((i for i, h in enumerate(header) if h.lower() in _BC_HEADER), 0)
    width = len(header)
    bad = sum(1 for r in body if len(r) != width)
    if bad:
        raise SystemExit(f"[userclusters] {path.name}: {bad} rows do not have {width} fields")
    bcs = [r[bi].strip() for r in body]
    names = [h for i, h in enumerate(header) if i != bi]
    cols = [[r[i].strip() for r in body] for i in range(width) if i != bi]
    _no_duplicates(bcs, path.name)
    return bcs, names, cols


def _no_duplicates(bcs, what: str):
    if len(set(bcs)) != len(bcs):
        n = len(bcs) - len(set(bcs))
        raise SystemExit(f"[userclusters] {what}: {n} barcodes appear more than once. One row per barcode: "
                         f"for a merged object pick one sample with --user-barcode-prefix / --user-barcode-suffix "
                         f"or export that sample alone.")


def _as_float(v: list[str]) -> np.ndarray | None:
    """float64 array (NaN for NA) if every non-NA value is a number, else None."""
    out = np.empty(len(v), np.float64)
    for i, x in enumerate(v):
        if x.lower() in _NA:
            out[i] = np.nan
            continue
        try:
            out[i] = float(x)
        except ValueError:
            return None
    return out


def read_embedding_csv(path: Path) -> dict:
    bcs, names, cols = read_table(path)
    num = [(n, a) for n, c in zip(names, cols) if (a := _as_float(c)) is not None]
    if len(num) < 2:
        raise SystemExit(f"[userclusters] {path.name}: needs two numeric columns after the barcode "
                         f"(found {[n for n, _ in num]} among {names})")
    (n1, x), (n2, y) = num[:2]
    if len(num) > 2:
        print(f"[userclusters] {path.name}: using the first two numeric columns '{n1}', '{n2}' as the embedding")
    return {"barcodes": bcs, "xy": np.column_stack([x, y]).astype(np.float32), "embName": f"{n1} / {n2}",
            "what": "embedding CSV"}


def read_obs_csv(path: Path) -> dict:
    bcs, names, cols = read_table(path)
    return {"barcodes": bcs, "columns": {n: (c, None) for n, c in zip(names, cols)}, "colors": {},
            "what": f"obs CSV {path.name}"}


def read_h5ad(path: Path, embedding: str | None) -> dict:
    import h5py
    from rnaseq import _h5_column, _s
    with h5py.File(path, "r") as f:
        og = f["obs"]
        idx = _s(og.attrs.get("_index", "_index"))
        bcs = [_s(b).strip() for b in og[idx][()]]
        _no_duplicates(bcs, path.name)
        columns, colors = {}, {}
        for name in og.keys():
            if name in (idx, "__categories"):
                continue
            obj = og[name]
            try:
                vals = _h5_column(og, name)
            except Exception:
                continue
            order = None
            if isinstance(obj, h5py.Group) and "categories" in obj:
                order = [_s(c) for c in obj["categories"][()]]
                pal = f.get(f"uns/{name}_colors")
                if pal is not None:
                    pal = [_s(c) for c in pal[()]]
                    if len(pal) == len(order):
                        colors[name] = {c: p for c, p in zip(order, pal)}
            if getattr(vals, "dtype", None) is not None and vals.dtype.kind in "fiub":
                vals = [("nan" if (isinstance(x, float) and np.isnan(x)) else str(x.item() if hasattr(x, "item") else x))
                        for x in vals]
            else:
                vals = [_s(x) for x in vals]
            columns[name] = (vals, order)
        xy, emb_name = None, None
        obsm = f.get("obsm")
        keys = [embedding] if embedding else ["X_umap", "X_UMAP", "umap", "UMAP"]
        if obsm is not None:
            for k in keys:
                if k in obsm and isinstance(obsm[k], h5py.Dataset):
                    xy = np.asarray(obsm[k][()], dtype=np.float32)[:, :2]
                    emb_name = k
                    break
        if xy is None and embedding:
            raise SystemExit(f"[userclusters] {path.name}: no '{embedding}' in obsm "
                             f"(found: {list(obsm.keys()) if obsm is not None else 'no obsm'})")
        if xy is None:
            print(f"[userclusters] {path.name}: no UMAP in obsm - only the clusterings are used, on SODA's UMAP")
    return {"barcodes": bcs, "columns": columns, "colors": colors, "xy": xy, "embName": emb_name,
            "what": f"h5ad {path.name}"}


def _palettes(columns: dict) -> dict:
    """Colour columns of a CSV: '<col>_color' / '<col>.colour' / '<col>_colors' holding #rrggbb per row,
    or a bare 'color' column next to exactly one other column (Xenium Explorer's cell groups CSV:
    cell_id,group,color)."""
    out = {}
    for name, (vals, _) in list(columns.items()):
        m = _COLOR_COL.match(name)
        target = m.group(1) if m else None
        if m is None and name.lower() in ("color", "colour", "colors", "colours") and len(columns) == 2:
            target = next(k for k in columns if k != name)
        if target is None or target not in columns:
            continue
        cats = columns[target][0]
        pal = {}
        for c, v in zip(cats, vals):
            h = _HEX.match(v.strip())
            if h and c.lower() not in _NA:
                pal.setdefault(c, "#" + h.group(1).lower())
        out[target] = pal
        del columns[name]
    return out


def categorical(values: list[str], order: list[str] | None,
                ints_up_to: int = MAX_INT_CATEGORIES) -> tuple[list[str] | None, np.ndarray | None, str]:
    """(categories, codes with -1 = NA, why-not) of one column; categories None when it is not one.

    Text columns, and integer columns with at most `ints_up_to` values (seurat_clusters), are
    categories; other numeric columns are not (yet) shown. 2 .. MAX_CATEGORIES values.
    """
    na = np.array([v.lower() in _NA for v in values])
    num = _as_float(values)
    vals = list(values)
    if num is not None:
        u = np.unique(num[~np.isnan(num)])
        if len(u) > ints_up_to or not np.all(np.mod(u, 1) == 0):
            return None, None, "numeric"
        vals = ["NA" if np.isnan(x) else str(int(x)) for x in num]
        if order is not None:
            order = [str(int(float(o))) for o in order if o.lower() not in _NA]
    present = sorted({v for v, n in zip(vals, na) if not n})
    if order is not None:
        cats = [c for c in dict.fromkeys(order) if c in set(present)]
        cats += [c for c in present if c not in set(cats)]
    else:
        cats = sorted(present, key=lambda c: (not c.lstrip("-").isdigit(),
                                              int(c) if c.lstrip("-").isdigit() else 0, c.lower(), c))
    if len(cats) < 2:
        return None, None, f"{len(cats)} value"
    if len(cats) > MAX_CATEGORIES:
        return None, None, f"{len(cats)} values (more than {MAX_CATEGORIES}: an id column?)"
    pos = {c: i for i, c in enumerate(cats)}
    codes = np.array([-1 if n else pos[v] for v, n in zip(vals, na)], np.int64)
    return cats, codes, ""


def deidentify(name: str, cats: list[str], codes: np.ndarray, colors: dict | None, ctx: str, log):
    """(safe name, categories, codes, colours) as rnaseq.write_obs does it: a sample-identity column
    has every value replaced by a number, any other value goes through the ID pattern. Values that
    become equal are merged."""
    safe_name = deid.anonymize(str(name), form="text", context=f"user obs column {ctx}")
    if deid.is_sample_column(name, cats):
        anon = [deid.sample_value(c, context=f"user sample column {ctx}") for c in cats]
        log(f"[deid] user column '{safe_name}': {len(cats)} sample values replaced by numbers")
    else:
        anon = [deid.anonymize(c, form="value", context=f"user obs {ctx}") for c in cats]
    new = list(dict.fromkeys(anon))
    pos = {c: i for i, c in enumerate(new)}
    remap = np.array([pos[a] for a in anon], np.int64)
    codes = np.where(codes >= 0, remap[np.maximum(codes, 0)], -1)
    cols = None
    if colors:
        first = {}
        for c, a in zip(cats, anon):
            h = _HEX.match(str(colors.get(c, "")).strip())
            if h:
                first.setdefault(a, "#" + h.group(1).lower())
        if all(a in first for a in new):
            cols = [first[a] for a in new]
    return safe_name, new, codes, cols


# ---------------------------------------------------------------------------------------
# barcodes of a layer, and matching
# ---------------------------------------------------------------------------------------
def layer_barcodes(er, ds_dir: Path | None = None, log=print) -> list[str]:
    """Barcodes of the layer in xy order: from the cache (expr.build_expr keeps them), else from the
    Space Ranger outs the cache manifest names (caches made before barcodes were kept)."""
    key = er.meta.get("key")
    n = len(er.xy)
    if "barcodes" in er.g:
        bcs = [b.decode() if isinstance(b, bytes) else str(b) for b in np.asarray(er.g["barcodes"])]
    else:
        bcs = None
        src = {}
        if ds_dir is not None and (Path(ds_dir) / "manifest.json").exists():
            src = json.loads((Path(ds_dir) / "manifest.json").read_text(encoding="utf-8")).get("source") or {}
        if src.get("outs"):
            import h5py
            from expr import find_visium_layers
            for L in find_visium_layers(Path(src["outs"])):
                if L["key"] == key:
                    with h5py.File(L["h5"], "r") as f:
                        bcs = [b.decode() for b in f["matrix"]["barcodes"][()]]
                    log(f"[userclusters] {key}: barcodes read from the Space Ranger outs (this cache predates "
                        f"barcodes in the cache)")
                    break
        if bcs is None:
            raise SystemExit(f"[userclusters] {key}: the cache has no barcodes for this layer and its source files "
                             f"are not reachable. Re-run `prepare.py visium` (it now keeps them in the cache).")
    if len(bcs) != n:
        raise SystemExit(f"[userclusters] {key}: {len(bcs)} barcodes but {n} bins in the cache")
    return bcs


def match_barcodes(user: list[str], layer: list[str], row_of: np.ndarray, prefix=None, suffix=None,
                   what: str = "", log=print) -> tuple[np.ndarray, dict]:
    """positions.bin row per user barcode (-1 = not found, -2 = found but off the image / cropped).

    Tries, per barcode: exact; the layer barcode without its "-1"-style suffix vs the user's without
    theirs; then the same after dropping a prefix up to one of the first three separators (_ : # |).
    `prefix` / `suffix`: only barcodes with it count (the rest are another sample and ignored), and
    it is removed before matching. Fails when several user barcodes land on one row from different
    sample tags (a merged object): pick a sample with prefix / suffix then.
    """
    t0 = time.time()
    exact = {b: i for i, b in enumerate(layer)}
    core: dict[str, int] = {}
    clash = set()
    for i, b in enumerate(layer):
        c = _SUFFIX.sub("", b)
        if c in core:
            clash.add(c)
        core[c] = i
    for c in clash:                                   # a stripped form two barcodes share is no evidence
        del core[c]
    nu = len(user)
    idx = np.full(nu, -1, np.int64)
    method = np.zeros(nu, np.int8)                    # 0 none, 1 exact, 2 suffix, 3 prefix, 4 prefix + suffix
    tag = [""] * nu
    ignored = 0
    for k, u0 in enumerate(user):
        u = u0
        if prefix:
            if not u.startswith(prefix):
                method[k] = -1; ignored += 1; continue
            u = u[len(prefix):]
        if suffix:
            if not u.endswith(suffix):
                method[k] = -1; ignored += 1; continue
            u = u[: len(u) - len(suffix)]
        j = exact.get(u)
        if j is not None:
            idx[k], method[k] = j, 1
            continue
        cu = _SUFFIX.sub("", u)
        j = core.get(cu)
        if j is not None:
            idx[k], method[k], tag[k] = j, 2, u[len(cu):]
            continue
        for m in list(_SEPS.finditer(u))[:3]:
            rest = u[m.end():]
            j = exact.get(rest)
            if j is not None:
                idx[k], method[k], tag[k] = j, 3, u[:m.end()]
                break
            cr = _SUFFIX.sub("", rest)
            j = core.get(cr)
            if j is not None:
                idx[k], method[k], tag[k] = j, 4, u[:m.end()] + "…" + rest[len(cr):]
                break
    considered = nu - ignored
    found = idx >= 0
    rows = np.where(found, row_of[np.maximum(idx, 0)], -1)
    off = found & (rows < 0)
    rows[off] = -2
    # one row claimed by several barcodes
    on = np.flatnonzero(rows >= 0)
    u_rows, first, counts = np.unique(rows[on], return_index=True, return_counts=True)
    dup_rows = set(u_rows[counts > 1].tolist())
    n_dup = int(sum(counts[counts > 1]) - len(dup_rows))
    if dup_rows:
        dk = [k for k in on if rows[k] in dup_rows]
        tags = {}
        for k in dk:
            tags[tag[k] or "(none)"] = tags.get(tag[k] or "(none)", 0) + 1
        top = sorted(tags.items(), key=lambda kv: -kv[1])[:6]
        if len(tags) > 1:
            raise SystemExit(
                f"[userclusters] {what}: {n_dup:,} barcodes land on a bin another barcode already took, with "
                f"different sample tags {top}. The file holds several samples (a merged object): keep one with "
                f"--user-barcode-prefix or --user-barcode-suffix (e.g. the prefix or suffix of this section).")
        keep = np.zeros(nu, bool)
        keep[on[first]] = True
        drop = on[~keep[on]]
        rows[drop] = -1
        log(f"[userclusters] WARNING {what}: {n_dup:,} barcodes matched a bin already matched; the first kept")
    n_found = int(found.sum())
    rate = n_found / considered if considered else 0.0
    rep = {"inFile": int(nu), "considered": int(considered), "matched": n_found, "rate": round(rate, 4),
           "offImage": int(off.sum()),
           "methods": {"exact": int((method == 1).sum()), "suffix": int((method == 2).sum()),
                       "prefix": int((method == 3).sum()), "prefixSuffix": int((method == 4).sum())}}
    if ignored:
        rep["otherSample"] = int(ignored)
    log(f"[userclusters] {what}: {n_found:,} of {considered:,} barcodes matched ({rate:.1%}): "
        + ", ".join(f"{k} {v:,}" for k, v in rep["methods"].items() if v)
        + (f"; {rep['offImage']:,} of them are off the image (cropped)" if rep["offImage"] else "")
        + (f"; {ignored:,} of another sample ignored" if ignored else "") + f"  {time.time() - t0:.1f}s")
    if rate < MIN_MATCH:
        miss = [user[k] for k in range(nu) if method[k] == 0][:3]
        hint = ""
        hd = {m.group(1) for b in user[:2000] for m in [re.search(r"s_(\d{3})um_", b)] if m}
        lk = {m.group(1) for b in layer[:5] for m in [re.search(r"s_(\d{3})um_", b)] if m}
        if hd and lk and hd != lk:
            hint = f" Your barcodes look like {'/'.join(sorted(hd))} µm bins, this layer has {'/'.join(sorted(lk))} µm: check --user-layer."
        raise SystemExit(f"[userclusters] {what}: only {n_found:,} of {considered:,} barcodes ({rate:.1%}) match "
                         f"this layer (at least {MIN_MATCH:.0%} needed). Examples not found: {miss}; the layer's "
                         f"barcodes look like {layer[:2]}.{hint}")
    if n_found < considered:
        log(f"[userclusters] WARNING {what}: {considered - n_found:,} barcodes ({1 - rate:.1%}) not found in this "
            f"layer, e.g. {[user[k] for k in range(nu) if method[k] == 0][:3]}")
    return rows, rep


# ---------------------------------------------------------------------------------------
# the umap/ folder
# ---------------------------------------------------------------------------------------
def _write(ud: Path, name: str, buf: bytes):
    tmp = ud / (name + ".part")
    tmp.write_bytes(buf)
    tmp.replace(ud / name)


def _soda_base(ud: Path, n_rows: int, log) -> tuple[dict | None, np.ndarray | None]:
    """SODA's own index.json and coords.bin: kept aside as index_soda.json / coords_soda.bin the first
    time, read from there once index.json carries user results. (None, None) without a SODA UMAP."""
    ip = ud / "index.json"
    if ip.exists():
        cur = json.loads(ip.read_text(encoding="utf-8"))
        if "userImport" not in cur:                      # sectionumap's own, fresh or reused
            _write(ud, SODA_COORDS, (ud / "coords.bin").read_bytes())
            _write(ud, SODA_INDEX, ip.read_bytes())
    if not (ud / SODA_INDEX).exists() or not (ud / SODA_COORDS).exists():
        return None, None
    base = json.loads((ud / SODA_INDEX).read_text(encoding="utf-8"))
    co = np.fromfile(ud / SODA_COORDS, dtype="<f4").reshape(-1, 2)
    if len(co) != n_rows or base.get("nBins") != n_rows:
        log(f"[userclusters] {ud}: SODA's saved UMAP does not fit this layer ({len(co)} rows): not used")
        return None, None
    return base, co


def manifest_entry(index: dict) -> dict:
    """The manifest's view of an index.json, naming each clustering."""
    cl = []
    for c in index["clusterings"]:
        e = {k: c[k] for k in ("id", "resolution", "name", "source") if c.get(k) is not None}
        e["n"] = c["n"]
        cl.append(e)
    out = {"nEmbedded": index["nEmbedded"], "clusterings": cl, "default": index["default"]}
    if (index.get("embedding") or {}).get("source") == "user":
        out["embedding"] = "user"
    return out


def restore(d: Path, log=print) -> dict | None:
    """Put SODA's own umap/ back after an earlier import (an export run without user files)."""
    ud = Path(d) / "umap"
    ip = ud / "index.json"
    if not ip.exists() or "userImport" not in json.loads(ip.read_text(encoding="utf-8")):
        # sectionumap wrote a fresh index.json: leftovers of an older import are not listed in it
        for p in [*ud.glob("clusters_user*.bin"), ud / SODA_INDEX, ud / SODA_COORDS]:
            if ip.exists() and p.exists():
                p.unlink()
        return None
    if not ((ud / SODA_INDEX).exists() and (ud / SODA_COORDS).exists()):
        log(f"[userclusters] {ud}: earlier user results but no SODA UMAP saved - left as it is")
        return manifest_entry(json.loads(ip.read_text(encoding="utf-8")))
    _write(ud, "coords.bin", (ud / SODA_COORDS).read_bytes())
    _write(ud, "index.json", (ud / SODA_INDEX).read_bytes())
    # our own files: the folder is again exactly what sectionumap wrote
    for p in [*ud.glob("clusters_user*.bin"), ud / SODA_INDEX, ud / SODA_COORDS]:
        p.unlink()
    log(f"[userclusters] {ud.parent.name}: no user files given - SODA's own UMAP and clusters restored")
    return manifest_entry(json.loads(ip.read_text(encoding="utf-8")))


def _log_cp10k(er, h5_idx: np.ndarray) -> tuple[sparse.csr_matrix, np.ndarray]:
    """log1p CP10k (rows x genes) of the given bins (h5 / xy order), and their library sizes."""
    G, B = len(er.genes), len(er.xy)
    gm = sparse.csr_matrix((np.asarray(er.g["data"], dtype=np.float32), np.asarray(er.g["indices"]),
                            er.indptr), shape=(G, B))
    X = gm.T.tocsr()[h5_idx]
    del gm
    lib = np.asarray(X.sum(1)).ravel().astype(np.float64)
    scale = np.divide(1e4, lib, out=np.zeros_like(lib), where=lib > 0)
    Xn = (sparse.diags(scale) @ X).tocsr().astype(np.float32)
    Xn.data = np.log1p(Xn.data)
    return Xn, lib


def apply(er, d: Path, spec: dict, ds_dir: Path | None = None, log=print) -> dict:
    """Write the user's embedding / clusterings into <d>/umap/ (d = the export's expr/<key>/).
    Returns the manifest's "umap" entry."""
    t0 = time.time()
    d = Path(d)
    key = er.meta.get("key")
    ud = d / "umap"
    ud.mkdir(parents=True, exist_ok=True)
    on_image = ~np.isnan(er.xy[:, 0])
    row_of = np.where(on_image, np.cumsum(on_image) - 1, -1)
    n_rows = int(on_image.sum())
    if (d / "positions.bin").exists() and (d / "positions.bin").stat().st_size != 8 * n_rows:
        raise SystemExit(f"[userclusters] {key}: positions.bin does not have the cache's {n_rows} rows")
    h5_of_row = np.flatnonzero(on_image)
    layer = layer_barcodes(er, ds_dir, log)

    # ---- read ----------------------------------------------------------------------------
    tables = []
    if spec.get("h5ad"):
        tables.append(read_h5ad(spec["h5ad"], spec.get("embedding")))
    if spec.get("umap"):
        tables.append(read_embedding_csv(spec["umap"]))
    for p in spec.get("obs") or []:
        tables.append(read_obs_csv(p))
    for t in tables:
        if t.get("columns"):
            t["colors"] = {**_palettes(t["columns"]), **t.get("colors", {})}

    # ---- match ---------------------------------------------------------------------------
    reports = []
    for t in tables:
        t["rows"], rep = match_barcodes(t["barcodes"], layer, row_of, spec.get("prefix"), spec.get("suffix"),
                                        f"{key} {t['what']}", log)
        reports.append(rep)

    base, base_co = _soda_base(ud, n_rows, log)
    emb_t = next((t for t in tables if t.get("xy") is not None), None)
    if emb_t is None and base is None:
        raise SystemExit(f"[userclusters] {key}: clusters without an embedding go onto SODA's UMAP, and this "
                         f"layer has none (too few bins, --no-umap, or the clustering packages are missing). "
                         f"Give your UMAP too (--user-umap, or obsm in --user-h5ad).")

    # ---- clusterings ---------------------------------------------------------------------
    want = spec.get("columns")
    found_cols = {}
    for t in tables:
        for name, (vals, order) in (t.get("columns") or {}).items():
            if name in found_cols:
                log(f"[userclusters] {key}: column '{name}' is in more than one file - the first is used")
                continue
            found_cols[name] = (t, vals, order)
    if want:
        missing = [c for c in want if c not in found_cols]
        if missing:
            raise SystemExit(f"[userclusters] {key}: --user-columns {missing} not found; columns are "
                             f"{sorted(found_cols)}")
    user_cl, labels = [], []
    for name in (want or list(found_cols)):
        t, vals, order = found_cols[name]
        named = bool(want) or bool(_CLUSTER_NAME.search(name))
        cats, codes, why = categorical(vals, order, MAX_CATEGORIES if named else MAX_INT_CATEGORIES)
        if cats is None:
            msg = f"[userclusters] {key}: column '{deid.anonymize(name, 'text', 'user obs column')}' skipped ({why})"
            if want:
                raise SystemExit(msg + " - it was asked for with --user-columns")
            log(msg)
            continue
        safe, cats, codes, colors = deidentify(name, cats, codes, (t.get("colors") or {}).get(name), key, log)
        lab = np.full(n_rows, -1, np.int64)
        r = t["rows"]
        ok = (r >= 0) & (codes >= 0)
        lab[r[ok]] = codes[ok]
        i = len(user_cl) + 1
        entry = {"id": f"user{i}", "name": safe, "source": "user", "file": f"clusters_user{i}.bin",
                 "n": len(cats), "categories": cats}
        if colors:
            entry["colors"] = colors
        user_cl.append(entry)
        labels.append(lab)
    if not user_cl and emb_t is None:
        raise SystemExit(f"[userclusters] {key}: nothing to import - no usable categorical column and no embedding")

    # markers per user cluster: log1p CP10k as sectionumap, over the assigned rows with counts
    if user_cl:
        assigned = np.flatnonzero(np.any(np.stack(labels) >= 0, axis=0))
        Xn, lib = _log_cp10k(er, h5_of_row[assigned])
        has = lib > 0
        Xn = Xn[np.flatnonzero(has)]
        rows_m = assigned[has]
        for c, lab in zip(user_cl, labels):
            la = lab[rows_m]
            m = la >= 0
            c["sizes"] = [int(v) for v in np.bincount(lab[lab >= 0], minlength=c["n"])]
            c["nAssigned"] = int((lab >= 0).sum())
            with np.errstate(divide="ignore", invalid="ignore"):
                c["markers"] = sectionumap._markers(Xn[np.flatnonzero(m)], la[m], c["n"], er.genes)
            log(f"[userclusters] {key}: '{c['name']}': {c['n']} groups, {c['nAssigned']:,} of {n_rows:,} rows "
                f"assigned  {time.time() - t0:.0f}s")
        del Xn

    # ---- coordinates ---------------------------------------------------------------------
    if emb_t is not None:
        coords = np.full((n_rows, 2), np.nan, np.float32)
        r = emb_t["rows"]
        ok = r >= 0
        coords[r[ok]] = emb_t["xy"][ok]
        bad = ~np.isfinite(coords).all(1)
        coords[bad] = np.nan
        embedding = {"source": "user", "name": deid.anonymize(str(emb_t["embName"]), "text", "user embedding")}
    else:
        coords = base_co
        embedding = {"source": "soda"}
    fin = np.isfinite(coords).all(1)
    if not fin.any():
        raise SystemExit(f"[userclusters] {key}: no matched row has a place in the embedding")
    e = coords[fin]

    soda_cl = []
    if base is not None and not spec.get("replace"):
        soda_cl = [{**c, "name": c.get("name") or f"Leiden {c['resolution']:g}", "source": "soda"}
                   for c in base["clusterings"]]
    index = dict(base) if base is not None else {"nBins": n_rows}
    index.update({
        "nBins": n_rows, "nEmbedded": int(fin.sum()),
        "bounds": [round(float(v), 4) for v in (e[:, 0].min(), e[:, 1].min(), e[:, 0].max(), e[:, 1].max())],
        "embedding": embedding,
        "clusterings": user_cl + soda_cl,
        "default": 0 if user_cl else (base or {}).get("default", 0),
        "userImport": {"layer": key, "kind": "h5ad" if spec.get("h5ad") else "csv", "barcodes": reports,
                       "sodaClusterings": "replaced" if spec.get("replace") else "kept" if base else "none",
                       "markers": "as sectionumap: log2 fold change of mean CP10k vs rest, pct >= "
                                  f"{sectionumap.MIN_PCT:g}, lfc > 0, over the assigned rows with counts"},
    })
    if base is None:
        index.setdefault("params", {})
        index.setdefault("clusterings", [])
    if not index["clusterings"]:
        raise SystemExit(f"[userclusters] {key}: no clustering left (SODA's replaced, none of yours usable)")
    deid.assert_clean(index, f"umap/index.json of {key}")
    txt = json.dumps(index, separators=(",", ":"), allow_nan=False)
    if _BARCODE.search(txt):
        raise SystemExit("[userclusters] REFUSING to write index.json: it contains a barcode")

    files = {"coords.bin": np.ascontiguousarray(coords, dtype="<f4").tobytes()}
    for c, lab in zip(user_cl, labels):
        files[c["file"]] = np.where(lab >= 0, lab, NO_CLUSTER).astype("<u2").tobytes()
    for p in ud.glob("clusters_user*.bin"):          # our own files from an earlier import
        if p.name not in files:
            p.unlink()
    files["index.json"] = txt.encode("utf-8")        # last: its presence means the folder is complete
    for name, buf in files.items():
        _write(ud, name, buf)
    log(f"[userclusters] {key}: {len(user_cl)} user clustering(s)"
        + (f" + {len(soda_cl)} of SODA's" if soda_cl else "")
        + f", embedding {embedding['source']} ({int(fin.sum()):,} rows placed)  {time.time() - t0:.0f}s")
    return manifest_entry(index)


def for_layer(specs: list[dict] | None, er, d: Path, ds_dir: Path | None, umap_entry: dict | None,
              log=print) -> dict | None:
    """The export's hook: apply the spec for this layer, or put SODA's results back after an earlier
    import when no spec names it. Returns the manifest's umap entry (or `umap_entry` unchanged)."""
    key = er.meta.get("key")
    spec = next((s for s in specs or [] if s["layer"] == key), None)
    if spec is not None and spec.get("auto"):     # kept from a standard-format file (stdformats.py): never fatal
        try:
            return apply(er, d, spec, ds_dir, log)
        except SystemExit as e:
            log(f"[userclusters] {key}: the file's own clusterings / UMAP were not attached ({e})")
            ud = Path(d) / "umap"
            if ud.is_dir() and not any(ud.iterdir()):    # apply() made it and stopped before writing anything
                ud.rmdir()
            return restore(d, log) or umap_entry if umap_entry is not None else umap_entry
    if spec is not None:
        return apply(er, d, spec, ds_dir, log)
    if umap_entry is not None:
        return restore(d, log) or umap_entry
    return umap_entry


# ---------------------------------------------------------------------------------------
# `prepare.py userclusters`: the same on an existing export, without re-running it
# ---------------------------------------------------------------------------------------
def apply_to_export(ds_dir: Path, dest: Path, specs: list[dict], markers: Path | None = None,
                    celltypes: bool = True, log=print) -> Path:
    import celltypes as celltypes_mod
    from expr import ExprReader
    ds_dir = Path(ds_dir)
    m = json.loads((ds_dir / "manifest.json").read_text(encoding="utf-8"))
    out = Path(dest) / m["id"]
    mp = out / "manifest.json"
    if not mp.exists():
        raise SystemExit(f"[userclusters] no export at {out}: run `prepare.py export` first")
    sm = json.loads(mp.read_text(encoding="utf-8"))
    check_layers(specs, sm)
    for spec in specs:
        key = spec["layer"]
        er = ExprReader(ds_dir / "expr" / f"{key}.zarr")
        d = out / "expr" / key
        u = apply(er, d, spec, ds_dir, log)
        if celltypes:
            mk = markers
            old = d / "umap" / "celltypes.json"
            if mk is None and old.exists():           # the list the export used before
                f = (json.loads(old.read_text(encoding="utf-8")).get("markerSet") or {}).get("file") or ""
                cand = Path(__file__).resolve().parent / f
                mk = cand if f.startswith("markers/") and cand.exists() else None
                if mk is None and f:                  # e.g. an --markers file from outside markers/
                    log(f"[userclusters] WARNING {key}: the export's marker list ({f}) is not in "
                        f"web/SpatialViewer/markers/, so suggested names fall back to the default list. "
                        f"Pass the same --markers file as the export to keep its names.")
            c = celltypes_mod.write_celltypes(er, d, markers=mk, log=log)
            if c:
                u = {**u, "celltypes": c}
        sm["expr"] = [({**L, "umap": u} if L["key"] == key else L) for L in sm.get("expr", [])]
    mp.write_text(json.dumps(sm, indent=1), encoding="utf-8")
    log(f"[userclusters] done -> {mp}")
    return out
