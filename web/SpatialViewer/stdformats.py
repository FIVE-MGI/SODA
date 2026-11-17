r"""Standard formats as pipeline input: OME-Zarr images, spatial AnnData (.h5ad / AnnData-Zarr) and SpatialData.

    python prepare.py omezarr     --src image.ome.zarr            --id x --out <cache> [--series 0] [--z N]
    python prepare.py anndata     --src section.h5ad | table.zarr --id x --out <cache> [--image fullres.tif]
    python prepare.py spatialdata --src sdata.zarr                --id x --out <cache> [--table / --region /
                                  --coordinate-system / --image-element / --points]

Each one is turned into an ordinary SODA cache (image.zarr, expr/<key>.zarr, shapes/, molecules/, manifest.json),
exactly what the instrument readers write, so `prepare.py export` and the viewer need nothing new. The files
are read with h5py / zarr 2 / pyarrow (annstore.py, ngff.py): anndata, scanpy, spatialdata and pandas are not
needed. Reading these formats in the browser itself is not done (a possible later step).

OME-Zarr (kind omezarr): the image (level 0 of multiscales; a bioformats2raw layout's series 0 or --series),
channel names and colours from omero, the pixel size from the level-0 scale when the x / y axes carry a unit,
the middle z-plane (or --z). As `prepare.py tiff`: the export adds the protein grid layer.

Spatial AnnData (kind anndata), the Scanpy / Squidpy layout:
    obsm['spatial']          positions (full-resolution image pixels for Visium; --spatial-key for another key)
    X / layers['counts'] / raw/X   the counts: layers['counts'] first, then raw/X, then X (--matrix picks);
                             values that are not whole numbers are shown as given, with a warning
    var                      gene names: the index, or gene_symbols / feature_name / gene_name ... when present;
                             feature_types other than Gene Expression and control probes are left out
    uns['spatial'][library]  images (hires, else lowres; --image-res) and scalefactors (tissue_*_scalef,
                             spot_diameter_fullres, microns_per_pixel, bin_size_um); one library only
    obs categorical columns, obsm['X_umap'], uns['<column>_colors']
                             -> "your own clusters / UMAP" on the layer (userclusters.py, at export)
    --image <TIFF / OME-Zarr>  the full-resolution image instead (positions as they are, x --image-scale)
    no image at all          a neutral cell-density background; --coords-um says how many µm one unit is
SpatialData (kind spatialdata), a SpatialData .zarr (zarr v2):
    tables/<name>            an AnnData each (X, obs, var, obsm as above); uns/spatialdata_attrs names the
                             region element, region_key and instance_key. Every table whose region is a shapes
                             element becomes one layer (--table picks one, --region one region of a table)
    shapes/<region>          positions: circle centres (spot / dot size from radius) or polygon centroids
                             (polygons become the cell outlines; equal axis-aligned squares are bins)
    shapes/<other>           polygons with the same ids -> more outlines ("nucle..." = nucleus outlines)
    images/<name>            the image under it (the largest one in the coordinate system, or --image-element)
    points/<name>            molecules (x, y, the feature_key column; qv >= 20 when there is a qv column;
                             control probes left out); needs the pixel size
    coordinate systems       --coordinate-system, else "global" when every element has it, else the one they
                             share. Each element -> that system -> the image's pixels: identity, scale,
                             translation, affine, sequence (ngff.to_affine); anything else is refused
Refused, each with its reason and never with a sample id, barcode or file name in the message: zarr v3 stores,
an AnnData with several image libraries (save one section per file), a SpatialData table whose region is a
labels element, regions / tables / coordinate systems that leave a choice open (the message names the option),
time series, HCS plates.

Privacy: obs names, instance ids and library / element names stay in the cache. Category values and column
names go through deid.py at export (userclusters.py); the manifest carries no file or element name. Messages
show a name only when deid.find_ids sees no sample ID in it ([sample ID] otherwise).
"""
from __future__ import annotations

import csv
import json
import re
import time
from pathlib import Path

import numpy as np

import annstore
import deid
import ngff

CONTROL_PREFIXES = ("negcontrolprobe", "negcontrolcodeword", "blank", "unassignedcodeword", "deprecatedcodeword",
                    "intergenic", "antisense", "genomiccontrol", "negprb", "negative", "systemcontrol", "falsecode")
GENE_COLUMNS = ("gene_symbols", "feature_name", "gene_name", "gene_names", "features", "name")
_VISIUM_BC = re.compile(r"[ACGT]{16}-\d{1,3}\Z")
VISIUM_SPOT_UM = 55.0
QV_MIN = 20.0
USER_DIR = "userobs"                         # in the cache: the obs columns / UMAP for userclusters.py at export
MAX_TABLES = 16


def log_name(name) -> str:
    """A file / element / library name for a message: as it is, or [sample ID] when it holds one."""
    s = str(name)
    return "[sample ID]" if deid.find_ids(s) else s


def is_control(name: str) -> bool:
    return str(name).strip().lower().startswith(CONTROL_PREFIXES)


def unique(names: list[str]) -> list[str]:
    """Gene names made unique the way scanpy's var_names_make_unique does (A, A-1, A-2)."""
    seen, out = {}, []
    taken = set(names)
    for n in names:
        if n not in seen:
            seen[n] = 0
            out.append(n)
            continue
        k = seen[n]
        while True:
            k += 1
            cand = f"{n}-{k}"
            if cand not in taken:
                break
        seen[n] = k
        taken.add(cand)
        out.append(cand)
    return out


def median_spacing(xy: np.ndarray) -> float:
    """Median distance to the nearest other point (positions in any unit; at most 50,000 sampled)."""
    from scipy.spatial import cKDTree
    p = xy[np.isfinite(xy).all(axis=1)]
    if len(p) < 2:
        return 1.0
    if len(p) > 50_000:
        p = p[np.random.default_rng(0).choice(len(p), 50_000, replace=False)]
    d, _ = cKDTree(p).query(p, k=2)
    d = d[:, 1]
    d = d[d > 0]
    return float(np.median(d)) if d.size else 1.0


# ---------------------------------------------------------------------------------------
# AnnData pieces shared by `anndata` and the tables of `spatialdata`
# ---------------------------------------------------------------------------------------
def read_counts(ad: annstore.AnnFile, which: str | None, tag: str, log):
    """(X n_obs x genes CSR, gene names, where, counts?) with feature types / control probes left out."""
    from scipy import sparse
    layers = ad.keys(ad.get("layers"))
    raw_ok = ad.get("raw/X") is not None and ad.get("raw/var") is not None
    if which in (None, "", "auto"):
        where = "layers/counts" if "counts" in layers else ("raw/X" if raw_ok else "X")
    elif which in ("X", "raw"):
        where = "raw/X" if which == "raw" else "X"
        if which == "raw" and not raw_ok:
            raise SystemExit(f"[{tag}] --matrix raw: the file has no raw/X")
    else:
        if which not in layers:
            raise SystemExit(f"[{tag}] --matrix {which}: no such layer (layers: {', '.join(layers) or 'none'})")
        where = f"layers/{which}"
    X = ad.matrix(where)
    var = ad.get("raw/var") if where == "raw/X" else ad.get("var")
    names = ad.index(var)
    for col in GENE_COLUMNS:
        if col in var:
            try:
                c = ad.column(var, col)
            except ValueError:
                continue
            vals = c.strings()
            if (c.is_categorical or c.is_text) and all(v and v != "NA" for v in vals):
                names = vals
                break
    if X.shape[1] != len(names):
        raise SystemExit(f"[{tag}] {where} has {X.shape[1]} columns for {len(names)} genes")
    keep = np.ones(len(names), bool)
    if "feature_types" in var:
        try:
            ft = [v.lower() for v in ad.column(var, "feature_types").strings()]
            if "gene expression" in ft:
                keep &= np.array([v == "gene expression" for v in ft])
        except ValueError:
            pass
    ctrl = np.array([is_control(n) for n in names])
    keep &= ~ctrl
    if not keep.all():
        log(f"[{tag}] {int((~keep).sum())} features left out ({int(ctrl.sum())} control probes, the rest not Gene Expression)")
        X = X[:, np.flatnonzero(keep)]
        names = [n for n, k in zip(names, keep) if k]
    X = sparse.csr_matrix(X, dtype=np.float32)
    X.eliminate_zeros()
    counts = annstore.looks_like_counts(X)
    log(f"[{tag}] counts from {where}: {X.shape[0]:,} rows x {X.shape[1]:,} genes, nnz {X.nnz:,}"
        + ("" if counts else " - WARNING these are not whole numbers (normalised values?): SODA's QC, UMAP and "
                             "marker scores assume counts; the values are shown as they are. --matrix picks another"))
    return X, unique([str(n) for n in names]), where, counts


def obs_clusterings(ad: annstore.AnnFile, rows: np.ndarray | None, skip=()) -> tuple[dict, dict]:
    """({column: values per row as text}, {column: {category: #colour}}) of the obs columns that can be
    clusterings: categorical and text columns. The rest (numbers, ids) is left to userclusters' own rules."""
    obs = ad.get("obs")
    idx_name = ad.attrs(obs).get("_index", "_index")
    out, colors = {}, {}
    for name in ad.columns(obs):
        if name in skip or name == idx_name:
            continue
        try:
            c = ad.column(obs, name)
        except (ValueError, KeyError, TypeError):
            continue
        if not (c.is_categorical or c.is_text):
            continue
        vals = c.strings()
        if rows is not None:
            vals = [vals[i] for i in rows.tolist()]
        out[name] = vals
        pal = ad.get(f"uns/{name}_colors")
        if c.is_categorical and pal is not None and not ad.is_group(pal):
            p = [annstore._s(x) for x in np.asarray(ad.array(pal)).ravel().tolist()]
            cats = [annstore._s(x) for x in c.categories]
            if len(p) == len(cats):
                colors[name] = dict(zip(cats, p))
    return out, colors


def write_user_tables(ds_dir: Path, key: str, barcodes: list[str], columns: dict, colors: dict, umap,
                      log, tag: str) -> dict | None:
    """The obs clusterings and the UMAP as the CSVs userclusters.py reads (into the cache, never the export).
    Returns the manifest's userObs entry, or None when there is nothing."""
    if not columns and umap is None:
        return None
    d = Path(ds_dir) / USER_DIR / key
    d.mkdir(parents=True, exist_ok=True)
    entry = {"layer": key, "obs": None, "umap": None}
    if columns:
        names = list(columns)
        used = set(names)
        col_names = {}
        for n in names:                               # a colour column for each column with a palette
            if n in colors:
                cn = f"{n}_color"
                while cn in used:
                    cn += "_"
                used.add(cn)
                col_names[n] = cn
        with open(d / "obs.csv", "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["barcode"] + names + [col_names[n] for n in names if n in col_names])
            for i, b in enumerate(barcodes):
                w.writerow([b] + [columns[n][i] for n in names]
                           + [colors[n].get(columns[n][i], "") for n in names if n in col_names])
        entry["obs"] = f"{USER_DIR}/{key}/obs.csv"
    if umap is not None:
        with open(d / "umap.csv", "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["barcode", "UMAP_1", "UMAP_2"])
            for b, (x, y) in zip(barcodes, np.asarray(umap, np.float64)[:, :2].tolist()):
                w.writerow([b, repr(x), repr(y)])
        entry["umap"] = f"{USER_DIR}/{key}/umap.csv"
    log(f"[{tag}] {key}: {len(columns)} obs column(s) kept as your own clusterings"
        + (" and obsm X_umap as your own UMAP" if umap is not None else "") + " (attached at export, userclusters.py)")
    return entry


def user_specs(ds_dir: Path, m: dict, given: list[dict] | None) -> list[dict]:
    """export: the specs of the obs columns a standard-format prepare kept (manifest userObs), for every layer
    the user did not give files for (their own --user-* files win)."""
    import userclusters
    out = list(given or [])
    have = {s["layer"] for s in out}
    for e in m.get("userObs") or []:
        if e.get("layer") in have:
            continue
        obs = [Path(ds_dir) / e["obs"]] if e.get("obs") else []
        umap = Path(ds_dir) / e["umap"] if e.get("umap") else None
        if not obs and umap is None:
            continue
        if not all(p.exists() for p in obs + ([umap] if umap else [])):
            continue
        spec = userclusters.make_spec(layer=e["layer"], umap=umap, obs=obs)
        spec["auto"] = True                      # a problem with it is a message, not a failed export
        out.append(spec)
    return out


def write_layer(ds_dir: Path, key: str, label: str, xy_px: np.ndarray, X, genes: list[str], barcodes: list[str],
                spot_px: float, mpp: float | None, kind: str | None = None, bin_um: float | None = None,
                extra: dict | None = None, log=print, tag: str = "") -> dict:
    """One expression layer in the cache format (expr._write_layer) from n x genes counts."""
    from expr import _write_layer
    t0 = time.time()
    gm = X.T.tocsr()
    gm.sort_indices()
    n = X.shape[0]
    meta = {"key": key, "label": label, "nBins": int(n), "nGenes": int(len(genes)), "binSizeUm": bin_um,
            "spotDiameterPx": float(spot_px), "micronsPerPixel": float(mpp) if mpp else None, "imageScale": 1.0}
    if kind == "cells":
        meta = {"key": key, "label": label, "kind": "cells", "nBins": int(n), "nGenes": int(len(genes)), "binSizeUm": None,
                "cellDiameterUm": round(float(spot_px) * float(mpp), 2) if mpp else None,
                "spotDiameterPx": float(spot_px), "micronsPerPixel": float(mpp) if mpp else None, "imageScale": 1.0}
    meta.update(extra or {})
    _write_layer(Path(ds_dir), np.asarray(xy_px, np.float32), gm, list(genes), meta, barcodes)
    on = int(np.isfinite(xy_px).all(axis=1).sum())
    log(f"[{tag}] {key}: {n:,} {'cells' if kind == 'cells' else 'spots / bins'} ({on:,} placed) x {len(genes):,} genes, "
        f"nnz {gm.nnz:,}, dot {spot_px:.1f} px ({time.time() - t0:.1f}s)")
    return meta


def density_background(xy: np.ndarray, units_um: float | None, spacing: float, bg_um: float | None, what: str, log, tag):
    """A neutral density image under positions given in some unit. Returns (source, to_px(x, y), µm per px or None,
    note). With `units_um` (µm per unit): pixels of bg_um µm (default 1, coarser above 250 Mpx); without: pixels of
    a quarter of the median spacing."""
    import flatfiles as ff
    ok = np.isfinite(xy).all(axis=1)
    if not ok.any():
        raise SystemExit(f"[{tag}] no finite positions to draw a background from")
    xl = (float(xy[ok, 0].min()), float(xy[ok, 0].max()))
    yl = (float(xy[ok, 1].min()), float(xy[ok, 1].max()))
    if units_um:
        W, H, bg, x0, y0 = ff.background_grid(xl, yl, 1.0 / units_um, bg_um=bg_um or ff.BACKGROUND_UM)
        step = bg / units_um                                   # units per pixel
        um = bg
    else:
        step = max(spacing / 4.0, 1e-9)
        margin = 4 * spacing
        while ((xl[1] - xl[0] + 2 * margin) / step) * ((yl[1] - yl[0] + 2 * margin) / step) > ff.BACKGROUND_MAX_PX:
            step *= 2
        W = max(2, int(np.ceil((xl[1] - xl[0] + 2 * margin) / step)))
        H = max(2, int(np.ceil((yl[1] - yl[0] + 2 * margin) / step)))
        x0, y0, um = xl[0] - margin, yl[0] - margin, None

    def to_px(x, y):
        return (np.asarray(x, np.float64) - x0) / step, (np.asarray(y, np.float64) - y0) / step
    dens = ff.Density(W, H)
    dens.add(*to_px(xy[ok, 0], xy[ok, 1]))
    sigma = max(1.0, 0.5 * spacing / step)
    src = ff.ArraySource([f"{what.capitalize()} density (no image)"], dens.image(sigma))
    log(f"[{tag}] no image: background = {what} density, {W} x {H} px"
        + (f" at {um:g} µm per pixel" if um else " (unit of the positions unknown: no scale bar; --coords-um sets it)"))
    return src, to_px, um, f"{what} density"


class Plan:
    """What prepare.py needs: the image source and pixel size before the pyramid, the layers after it."""

    def __init__(self, tag: str):
        self.tag = tag
        self.source = None
        self.pixel_size_um = None
        self.manifest_source = {}
        self.background = None
        self.colors = None
        self.imaging = False
        self.steps = []                         # callables (ds_dir, manifest, args, log) run by finish
        self.block = {}

    def finish(self, ds_dir: Path, manifest: dict, args, log=print):
        if self.background:
            manifest["image"]["background"] = self.background
            manifest["image"]["channels"][0]["color"] = "#b4b4b4"
        if self.imaging:
            manifest["imaging"] = True
        for step in self.steps:
            step(Path(ds_dir), manifest, args, log)
        if self.block:
            manifest[self.tag] = self.block


def _image_source(path: Path, args, log, tag):
    """--image: a TIFF / OME-TIFF or an OME-Zarr folder."""
    import pyramid
    if path.is_dir():
        img, note = ngff.open_image(path, getattr(args, "series", None), what=tag)
        if note:
            log(f"[{tag}] image: {note}")
        return img.source(z=getattr(args, "z", None), log=log), img.colors
    return pyramid.TiffStackSource(path, rgb_split=True), None


def _is_photo(src) -> bool:
    return bool(getattr(src, "is_rgb", False)) or [n.lower() for n in src.names] == ["red", "green", "blue"]


# ---------------------------------------------------------------------------------------
# omezarr
# ---------------------------------------------------------------------------------------
def add_args_omezarr(p):
    p.add_argument("--src", required=True, help="the OME-Zarr image (a .zarr / .ome.zarr folder, NGFF 0.4 or older, zarr v2)")
    p.add_argument("--series", help="bioformats2raw layout: the series to read (default 0)")
    p.add_argument("--z", type=int, help="z-stack: the plane to read (default the middle one)")
    p.add_argument("--channels", help="comma-separated channel names (override the omero labels)")


def plan_omezarr(args, log=print) -> Plan:
    tag = "omezarr"
    img, note = ngff.open_image(Path(args.src), args.series, what=tag)
    names = [s.strip() for s in args.channels.split(",")] if args.channels else None
    P = Plan(tag)
    P.source = img.source(z=args.z, names=names, log=log)
    P.colors = None if img.is_rgb else img.colors
    P.pixel_size_um = img.pixel_size_um
    P.manifest_source = {"folder": str(args.src)}
    P.block = {"ngffVersion": img.version or None, "channels": img.nC, "zPlanes": img.nZ, "zPlane": P.source.z if img.nZ > 1 else None,
               "levels": len(img.level_paths), "dtype": img.dtype.name}
    log(f"[{tag}] NGFF {img.version or '?'}{', ' + note if note else ''}: {img.width} x {img.height} px, {img.nC} channel(s) "
        f"({', '.join(log_name(n) for n in P.source.names[:12])}{' ...' if img.nC > 12 else ''}), {img.dtype.name}, "
        f"{len(img.level_paths)} level(s)" + (f", z-plane {P.source.z} of {img.nZ}" if img.nZ > 1 else ""))
    if img.pixel_size_um:
        log(f"[{tag}] pixel size {img.pixel_size_um:.4f} µm from the level-0 scale")
    else:
        log(f"[{tag}] no pixel size in the file ({img.pixel_note}): give --pixel-size for the scale bar and protein squares")
    return P


# ---------------------------------------------------------------------------------------
# anndata
# ---------------------------------------------------------------------------------------
def add_args_anndata(p):
    p.add_argument("--src", required=True, help="the AnnData file: .h5ad, or an AnnData .zarr folder (zarr v2)")
    p.add_argument("--image", help="the full-resolution image the positions are in (TIFF / OME-TIFF / OME-Zarr); "
                                   "default: the image in uns['spatial'], else a cell-density background")
    p.add_argument("--image-scale", type=float, default=1.0, help="with --image: multiply the positions by this")
    p.add_argument("--image-res", choices=["hires", "lowres"], default="hires",
                   help="which uns['spatial'] image to use (default hires; lowres when there is no hires)")
    p.add_argument("--spatial-key", default="spatial", help="the obsm key of the positions (default spatial)")
    p.add_argument("--library", help="uns['spatial'] library to use when there are several: its name or number (1 = first); "
                                     "the rows of the other libraries are left out by obs['library_id']")
    p.add_argument("--matrix", help="the counts: X, raw, or a layer name (default: layers['counts'], else raw/X, else X)")
    p.add_argument("--coords-um", type=float,
                   help="µm per unit of obsm['spatial'] when the file cannot say (no scale factors): sets the scale bar")
    p.add_argument("--layer-kind", choices=["auto", "spots", "cells"], default="auto",
                   help="what the rows are: spots (Visium-like bins) or cells; auto: spots when the scale factors give a "
                        "spot size, cells otherwise")
    p.add_argument("--no-user-clusters", action="store_true",
                   help="do not keep the obs columns / X_umap as your own clusterings / UMAP")
    p.add_argument("--background-um", type=float, help="no image, --coords-um given: µm per pixel of the background (default 1)")


def plan_anndata(args, log=print) -> Plan:
    tag = "anndata"
    src = Path(args.src)
    ad = annstore.AnnFile(src, what=tag)
    P = Plan(tag)
    P.manifest_source = {"file": str(src)}
    try:
        xy = ad.obsm(args.spatial_key)
        if xy is None:
            raise SystemExit(f"[{tag}] no obsm['{args.spatial_key}'] (positions): obsm holds "
                             f"{', '.join(ad.obsm_keys()) or 'nothing'}. Without positions this is single-cell data: "
                             f"prepare it with --kind rnaseq, or name the key with --spatial-key")
        if xy.ndim != 2 or xy.shape[1] < 2 or len(xy) != ad.n_obs:
            raise SystemExit(f"[{tag}] obsm['{args.spatial_key}'] has shape {list(xy.shape)}: not (cells, 2)")
        xy = xy[:, :2].astype(np.float64)
        rows = np.arange(ad.n_obs)
        # ---- libraries (uns['spatial'])
        sp = ad.get("uns/spatial")
        libs = [k for k in ad.keys(sp) if ad.is_group(sp[k])] if sp is not None else []
        lib = None
        if len(libs) > 1 and not args.library:
            raise SystemExit(f"[{tag}] uns['spatial'] holds {len(libs)} image libraries (sections merged into one object). "
                             f"SODA prepares one section at a time: give --library <name or number> (its rows are taken "
                             f"by obs['library_id']), or save each section to its own file")
        if libs:
            if args.library:
                lib = libs[int(args.library) - 1] if str(args.library).isdigit() and 0 < int(args.library) <= len(libs) \
                    else args.library
                if lib not in libs:
                    raise SystemExit(f"[{tag}] --library {log_name(args.library)}: not one of the {len(libs)} libraries "
                                     f"(give a number 1 .. {len(libs)})")
            else:
                lib = libs[0]
            if len(libs) > 1:
                if "library_id" not in ad.get("obs"):
                    raise SystemExit(f"[{tag}] several libraries but no obs['library_id'] to say which rows belong to "
                                     f"which: save each section to its own file")
                lid = ad.column(ad.get("obs"), "library_id").strings()
                rows = np.array([i for i, v in enumerate(lid) if v == lib], np.int64)
                if not len(rows):
                    raise SystemExit(f"[{tag}] no row has obs['library_id'] == the library asked for")
                log(f"[{tag}] library {libs.index(lib) + 1} of {len(libs)}: {len(rows):,} of {ad.n_obs:,} rows")
        sf = (ad.dict(f"uns/spatial/{lib}/scalefactors") or {}) if lib else {}
        sf = {k: (float(v[0]) if isinstance(v, list) and v else float(v)) for k, v in sf.items()
              if isinstance(v, (int, float)) or (isinstance(v, list) and v and isinstance(v[0], (int, float)))}
        barcodes = [ad.obs_names[i] for i in rows.tolist()]
        xy = xy[rows]
        visium_bc = bool(barcodes) and sum(1 for b in barcodes[:2000] if _VISIUM_BC.match(b)) >= 0.95 * min(len(barcodes), 2000)
        # ---- µm per full-resolution pixel
        mpp_full, how = None, ""
        if sf.get("microns_per_pixel"):
            mpp_full, how = sf["microns_per_pixel"], "scalefactors microns_per_pixel"
        elif visium_bc and sf.get("spot_diameter_fullres"):
            mpp_full, how = VISIUM_SPOT_UM / sf["spot_diameter_fullres"], \
                "classic Visium assumed (Visium barcodes and a spot diameter in the scale factors): 55 µm spots"
        elif args.coords_um:
            mpp_full, how = args.coords_um, "--coords-um"
        # ---- the image
        scale = 1.0
        if args.image:
            P.source, P.colors = _image_source(Path(args.image), args, log, tag)
            scale = float(args.image_scale or 1.0)
            P.manifest_source["image"] = str(args.image)
            to_px = (lambda x, y: (np.asarray(x) * scale, np.asarray(y) * scale))
            log(f"[{tag}] image: --image, {P.source.width} x {P.source.height} px; positions x {scale:g}")
        else:
            img, res = None, None
            if lib:
                for r in ([args.image_res] + [x for x in ("hires", "lowres") if x != args.image_res]):
                    img = ad.image(f"uns/spatial/{lib}/images/{r}")
                    if img is not None:
                        res = r
                        break
            if img is not None:
                sfk = f"tissue_{res}_scalef"
                if not sf.get(sfk):
                    raise SystemExit(f"[{tag}] the {res} image has no scalefactors {sfk}: its pixels cannot be matched "
                                     f"to obsm['spatial']; give the full-resolution image with --image")
                scale = sf[sfk]
                P.source = _array_source(img, log, tag)
                P.manifest_source["image"] = f"uns spatial {res}"
                to_px = (lambda x, y: (np.asarray(x) * scale, np.asarray(y) * scale))
                log(f"[{tag}] image: uns['spatial'] {res} image, {P.source.width} x {P.source.height} px, positions x {scale:g} "
                    f"({sfk})" + (" - the other one is not used" if res != args.image_res else ""))
            else:
                spacing = median_spacing(xy)
                P.source, to_px, um, P.background = density_background(
                    xy, mpp_full, spacing, args.background_um, "spot" if sf.get("spot_diameter_fullres") else "cell", log, tag)
                scale = (1.0 / (um / mpp_full)) if (um and mpp_full) else None
        if args.pixel_size:
            P.pixel_size_um, how = float(args.pixel_size), "--pixel-size"
            log(f"[{tag}] pixel size {P.pixel_size_um:.4f} µm (--pixel-size)")
        elif mpp_full and scale:
            P.pixel_size_um = mpp_full / scale
            log(f"[{tag}] pixel size {P.pixel_size_um:.4f} µm ({how})")
        elif getattr(P.source, "pixel_size_um", None) and args.image:
            P.pixel_size_um, how = float(P.source.pixel_size_um), "the --image file"
            log(f"[{tag}] pixel size {P.pixel_size_um:.4f} µm (from the --image file)")
        elif not P.background:
            log(f"[{tag}] no pixel size from the file (no microns_per_pixel, not classic Visium): give --pixel-size "
                f"or --coords-um for the scale bar")
        P.imaging = not _is_photo(P.source) and not P.background
        # ---- counts, genes, clusterings: read now (the file is open), written after the pyramid
        X, genes, where, counts = read_counts(ad, args.matrix, tag, log)
        X = X[rows]
        cols, colors = ({}, {}) if args.no_user_clusters else obs_clusterings(ad, rows if len(rows) != ad.n_obs else None)
        umap = None
        if not args.no_user_clusters:
            u = ad.obsm("X_umap")
            if u is not None and u.ndim == 2 and u.shape[1] >= 2 and len(u) == ad.n_obs:
                umap = u[rows, :2]
    finally:
        ad.close()
    xpx, ypx = to_px(xy[:, 0], xy[:, 1])
    xy_px = np.column_stack([xpx, ypx])
    # ---- what the rows are
    hd_bin = sf.get("bin_size_um")
    kind = args.layer_kind
    if kind == "auto":
        kind = "spots" if (sf.get("spot_diameter_fullres") or hd_bin) else "cells"
    if sf.get("spot_diameter_fullres") and scale:
        spot_px = sf["spot_diameter_fullres"] * scale
    elif hd_bin and P.pixel_size_um:
        spot_px = hd_bin / P.pixel_size_um
    else:
        spot_px = median_spacing(xy_px)
    if kind == "spots" and hd_bin:
        key, label, bin_um = f"{int(round(hd_bin)):03d}um", f"{hd_bin:g} µm bins", float(hd_bin)
    elif kind == "spots":
        key, bin_um = "spots", None
        label = "Visium spots (55 µm)" if how.startswith("classic Visium") else "Spots"
    else:
        key, label, bin_um = "cells", "Cells", None
    plat = "visium" if (visium_bc and key == "spots" and not hd_bin) else "anndata"
    P.block = {"source": "h5ad" if ad.kind == "h5ad" else "zarr", "counts": where, "countsAreIntegers": counts,
               "positions": f"obsm {args.spatial_key}", "libraries": len(libs), "pixelSize": how or None,
               "rows": int(len(rows)), "rowsTotal": int(ad.n_obs)}

    def step(ds_dir, manifest, a, lg):
        meta = write_layer(ds_dir, key, label, xy_px, X, genes, barcodes, spot_px, manifest.get("pixelSizeUm"),
                           kind="cells" if kind == "cells" else None, bin_um=bin_um, extra={"platform": plat},
                           log=lg, tag=tag)
        manifest["expr"].append(meta)
        e = write_user_tables(ds_dir, key, barcodes, cols, colors, umap, lg, tag)
        if e:
            manifest["userObs"] = [e]
    P.steps.append(step)
    return P


def _array_source(img: np.ndarray, log, tag):
    """A uns image (H, W), (H, W, 3 / 4), float 0..1 or integers, as a pyramid source."""
    import flatfiles as ff
    a = np.asarray(img)
    if a.dtype.kind == "f":
        top = 1.0 if np.nanmax(a) <= 1.0 + 1e-6 else float(np.nanmax(a))
        a = np.clip(np.nan_to_num(a) / top * 255.0 + 0.5, 0, 255).astype(np.uint8)
    elif a.dtype not in (np.uint8, np.uint16):
        a = np.clip(a, 0, 65535).astype(np.uint16)
    if a.ndim == 3 and a.shape[2] in (3, 4):
        s = ff.ArraySource(["Red", "Green", "Blue"], np.moveaxis(a[:, :, :3], 2, 0).copy())
        s.is_rgb = True
        return s
    if a.ndim == 2:
        return ff.ArraySource(["Image"], a)
    raise SystemExit(f"[{tag}] the uns['spatial'] image has shape {list(a.shape)}: not an image")


# ---------------------------------------------------------------------------------------
# spatialdata
# ---------------------------------------------------------------------------------------
ELEMENT_KINDS = ("images", "labels", "shapes", "points", "tables")


def add_args_spatialdata(p):
    p.add_argument("--src", required=True, help="the SpatialData store (a .zarr folder, zarr v2)")
    p.add_argument("--table", help="the table to read (default: every table that annotates shapes, one layer each)")
    p.add_argument("--region", help="a table that annotates several regions: the one to read")
    p.add_argument("--coordinate-system", help="the coordinate system to draw in (default global, else the one the "
                                                "elements share)")
    p.add_argument("--image-element", help="the image to draw under the data (default the largest in the coordinate system)")
    p.add_argument("--points", help="the points element for the molecules (default: the only one, or the one named "
                                    "transcripts)")
    p.add_argument("--matrix", help="the counts of the table(s): X, raw, or a layer name (default layers['counts'], else raw/X, else X)")
    p.add_argument("--z", type=int, help="an image with a z-stack: the plane to read (default the middle one)")
    p.add_argument("--no-transcripts", action="store_true", help="skip the molecules layer")
    p.add_argument("--no-user-clusters", action="store_true",
                   help="do not keep the tables' obs columns / X_umap as your own clusterings / UMAP")
    p.add_argument("--coords-um", type=float,
                   help="no image: µm per unit of the coordinate system (sets the scale bar and lets the molecules "
                        "be drawn); with an image its own pixel size (or --pixel-size) is used")
    p.add_argument("--background-um", type=float, help="no image, --coords-um given: µm per pixel of the background (default 1)")


def sd_elements(root: Path) -> dict[str, dict[str, Path]]:
    """{images|labels|shapes|points|tables: {name: folder}} of a SpatialData store (an older store's single
    table/ group counts as a table)."""
    out = {k: {} for k in ELEMENT_KINDS}
    for k in ELEMENT_KINDS:
        d = root / k
        if d.is_dir():
            for e in sorted(d.iterdir(), key=lambda q: q.name.lower()):
                if e.is_dir() and not e.name.startswith("."):
                    out[k][e.name] = e
    t = root / "table"
    if t.is_dir() and not out["tables"]:
        if ngff.read_attrs(t).get("encoding-type") == "anndata" or (t / "obs").is_dir():
            out["tables"]["table"] = t
        else:
            for e in sorted(t.iterdir()):
                if e.is_dir() and (e / "obs").is_dir():
                    out["tables"][e.name] = e
    return out


def is_spatialdata(path) -> bool:
    p = Path(path)
    if not p.is_dir() or ngff.zarr_version(p) is None:
        return False
    a = ngff.read_attrs(p)
    if "spatialdata_attrs" in a:
        return True
    if "multiscales" in a or a.get("encoding-type") == "anndata":
        return False
    return any((p / k).is_dir() and ngff.zarr_version(p / k) for k in ELEMENT_KINDS + ("table",))


def table_attrs(ad: annstore.AnnFile) -> dict:
    """uns/spatialdata_attrs of a table: {"region": [names], "region_key", "instance_key"}."""
    a = ad.dict("uns/spatialdata_attrs") or {}
    reg = a.get("region")
    if isinstance(reg, str):
        reg = [reg]
    elif isinstance(reg, np.ndarray):
        reg = [str(x) for x in reg.ravel().tolist()]
    return {"region": [str(r) for r in reg or []], "region_key": a.get("region_key"), "instance_key": a.get("instance_key")}


def _pandas_index(pf) -> list | None:
    md = pf.schema_arrow.metadata or {}
    raw = md.get(b"pandas")
    if not raw:
        return None
    try:
        return json.loads(raw.decode("utf-8")).get("index_columns")
    except ValueError:
        return None


def _wkb_points(blobs) -> np.ndarray:
    import struct
    out = np.full((len(blobs), 2), np.nan)
    for i, b in enumerate(blobs):
        if b is None or len(b) < 21:
            continue
        fmt = "<" if b[0] == 1 else ">"
        t = struct.unpack_from(fmt + "I", b, 1)[0]
        off = 5 + (4 if t & 0x20000000 else 0)
        if (t & 0x0FFFFFFF) % 1000 != 1:
            continue
        out[i] = struct.unpack_from(fmt + "dd", b, off)
    return out


def read_shapes(folder: Path, tag: str) -> dict:
    """A shapes element (GeoParquet, spatialdata shapes format 0.2): {"ids": [str], "kind": "circles" | "polygons",
    "xy": centres (n, 2), "radius" | ("start", "ring_xy", "ring_ids"), "attrs"}."""
    import pyarrow.parquet as pq
    import cellshapes
    attrs = ngff.read_attrs(folder)
    f = folder / "shapes.parquet"
    if not f.exists():
        if (folder / "coords").is_dir() or (folder / "coords" / ".zarray").exists() or "geos" in attrs:
            raise SystemExit(f"[{tag}] shapes {log_name(folder.name)}: the older zarr-array shapes layout (spatialdata "
                             f"shapes format 0.1). SODA reads the GeoParquet layout (format 0.2): save the store again with "
                             f"a current spatialdata")
        raise SystemExit(f"[{tag}] shapes {log_name(folder.name)}: no shapes.parquet")
    pf = pq.ParquetFile(str(f))
    names = pf.schema_arrow.names
    gcol = "geometry" if "geometry" in names else next((n for n in names if n.lower() in ("geom", "wkb")), None)
    if gcol is None:
        raise SystemExit(f"[{tag}] shapes {log_name(folder.name)}: no geometry column")
    idx = _pandas_index(pf)
    want = [gcol] + (["radius"] if "radius" in names else [])
    icol = None
    if idx and isinstance(idx[0], str) and idx[0] in names:
        icol = idx[0]
        want.append(icol)
    t = pf.read(columns=want)
    n = t.num_rows
    if icol:
        ids = [str(v) for v in t.column(icol).to_pylist()]
    elif idx and isinstance(idx[0], dict) and idx[0].get("kind") == "range":
        r = idx[0]
        ids = [str(v) for v in range(int(r.get("start", 0)), int(r.get("stop", n)), int(r.get("step", 1)))][:n]
    else:
        ids = [str(i) for i in range(n)]
    blobs = t.column(gcol).to_pylist()
    first = next((b for b in blobs if b), None)
    if first is None:
        raise SystemExit(f"[{tag}] shapes {log_name(folder.name)}: no geometry")
    fmt = "<" if first[0] == 1 else ">"
    import struct
    gtype = struct.unpack_from(fmt + "I", first, 1)[0] & 0x0FFFFFFF
    gtype %= 1000
    out = {"ids": ids, "attrs": attrs, "n": n}
    if gtype == 1:
        out["kind"] = "circles"
        out["xy"] = _wkb_points(blobs)
        if "radius" not in names:
            raise SystemExit(f"[{tag}] shapes {log_name(folder.name)}: points without a radius column")
        out["radius"] = np.asarray(t.column("radius").to_numpy(zero_copy_only=False), np.float64)
        return out
    if gtype not in (3, 6):
        raise SystemExit(f"[{tag}] shapes {log_name(folder.name)}: geometry type {gtype} (SODA reads points with a radius, "
                         f"polygons and multipolygons)")
    out["kind"] = "polygons"
    ring_ids, rings, cx, cy, area = [], [], np.full(n, np.nan), np.full(n, np.nan), np.zeros(n)
    nvert = np.zeros(n, np.int64)
    for i, b in enumerate(blobs):
        if b is None:
            continue
        rr, _holes = cellshapes._wkb_outer_rings(b)
        A = sx = sy = 0.0
        for r in rr:
            r = np.asarray(r, np.float64)
            if len(r) > 1 and np.allclose(r[0], r[-1]):
                r = r[:-1]
            if len(r) < 3:
                continue
            x, y = r[:, 0], r[:, 1]
            x1, y1 = np.roll(x, -1), np.roll(y, -1)
            cr = x * y1 - x1 * y
            a = cr.sum() / 2.0
            A += a
            sx += ((x + x1) * cr).sum() / 6.0
            sy += ((y + y1) * cr).sum() / 6.0
            ring_ids.append(ids[i])
            rings.append(r)
            nvert[i] += len(r)
        if abs(A) > 0:
            cx[i], cy[i], area[i] = sx / A, sy / A, abs(A)
    out["xy"] = np.column_stack([cx, cy])
    out["area"] = area
    out["nvert"] = nvert
    out["ring_ids"] = ring_ids
    cnt = np.array([len(r) for r in rings], np.int64)
    out["start"] = np.concatenate([[0], np.cumsum(cnt)])
    out["ring_xy"] = np.concatenate(rings) if rings else np.zeros((0, 2))
    return out


def _squares(sh: dict, M: np.ndarray) -> float | None:
    """The side (in image pixels) when every polygon is the same axis-aligned square (Visium HD bins), else None."""
    if sh["kind"] != "polygons" or not sh["n"]:
        return None
    if np.mean(sh["nvert"] == 4) < 0.99 or abs(M[0, 1]) > 1e-9 or abs(M[1, 0]) > 1e-9:
        return None
    a = sh["area"][sh["area"] > 0]
    if not a.size or np.std(a) / np.mean(a) > 0.01:
        return None
    return float(np.sqrt(np.median(a)) * np.sqrt(abs(M[0, 0] * M[1, 1])))


def _pick_cs(want: str | None, systems: list[dict], tag: str) -> str:
    """The coordinate system every element in `systems` (each {cs: matrix}) has."""
    common = set(systems[0]) if systems else set()
    for s in systems[1:]:
        common &= set(s)
    if want:
        if want not in common:
            raise SystemExit(f"[{tag}] --coordinate-system {log_name(want)}: not shared by the elements read "
                             f"(they share: {', '.join(log_name(c) for c in sorted(common)) or 'none'})")
        return want
    if "global" in common:
        return "global"
    if len(common) == 1:
        return next(iter(common))
    if not common:
        raise SystemExit(f"[{tag}] the elements read share no coordinate system: SODA cannot place them on one image")
    raise SystemExit(f"[{tag}] the elements share {len(common)} coordinate systems ("
                     f"{', '.join(log_name(c) for c in sorted(common))}): say which with --coordinate-system")


def plan_spatialdata(args, log=print) -> Plan:
    tag = "spatialdata"
    root = Path(args.src)
    ngff.refuse_v3(root, tag)
    if not is_spatialdata(root):
        raise SystemExit(f"[{tag}] not a SpatialData store (a zarr folder with images/, shapes/, points/ or tables/)")
    el = sd_elements(root)
    P = Plan(tag)
    P.manifest_source = {"folder": str(root)}
    log(f"[{tag}] elements: " + ", ".join(f"{len(v)} {k}" for k, v in el.items() if v))
    if not el["tables"]:
        raise SystemExit(f"[{tag}] no table (tables/): SODA needs a table of counts. An image on its own can be "
                         f"prepared with --kind omezarr --src <the store>/images/<name>")
    tables = list(el["tables"])
    if args.table:
        if args.table not in el["tables"]:
            raise SystemExit(f"[{tag}] --table {log_name(args.table)}: not one of the {len(tables)} tables")
        tables = [args.table]
    if len(tables) > MAX_TABLES:
        raise SystemExit(f"[{tag}] {len(tables)} tables: give the one to read with --table")
    # ---- each table: its region (a shapes element), its rows
    layers = []
    for tn in tables:
        ad = annstore.AnnFile(el["tables"][tn], what=tag)
        try:
            ta = table_attrs(ad)
            regs = ta["region"]
            if not regs:
                raise SystemExit(f"[{tag}] table {log_name(tn)}: no region in uns['spatialdata_attrs']: it annotates no "
                                 f"element, so its rows have no position")
            if args.region:
                if args.region not in regs:
                    if args.table:
                        raise SystemExit(f"[{tag}] --region {log_name(args.region)}: table {log_name(tn)} annotates "
                                         f"{', '.join(log_name(r) for r in regs)}")
                    continue
                reg = args.region
            elif len(regs) > 1:
                raise SystemExit(f"[{tag}] table {log_name(tn)} annotates {len(regs)} regions ("
                                 f"{', '.join(log_name(r) for r in regs)}): give the one to read with --region")
            else:
                reg = regs[0]
            if reg in el["labels"]:
                raise SystemExit(f"[{tag}] table {log_name(tn)} annotates the labels element {log_name(reg)} (a segmentation "
                                 f"mask). SODA reads positions from shapes (circles or polygons) only; a labels region is "
                                 f"not read yet")
            if reg not in el["shapes"]:
                if len(tables) > 1 and not args.table:
                    log(f"[{tag}] table {log_name(tn)}: its region is not a shapes element, skipped")
                    continue
                raise SystemExit(f"[{tag}] table {log_name(tn)}: its region {log_name(reg)} is not a shapes element of the store")
            obs = ad.get("obs")
            ik, rk = ta["instance_key"], ta["region_key"]
            if not ik or ik not in obs:
                raise SystemExit(f"[{tag}] table {log_name(tn)}: no instance_key column in obs")
            inst = ad.column(obs, ik).strings()
            rows = np.arange(ad.n_obs)
            if rk and rk in obs:
                rv = ad.column(obs, rk).strings()
                rows = np.array([i for i, v in enumerate(rv) if v == reg], np.int64)
            inst = [_norm_id(inst[i]) for i in rows.tolist()]
            X, genes, where, counts = read_counts(ad, args.matrix, tag, log)
            X = X[rows]
            skip = {c for c in (ik, rk) if c}
            cols, colors = ({}, {}) if args.no_user_clusters else obs_clusterings(ad, rows if len(rows) != ad.n_obs else None, skip)
            umap = None
            if not args.no_user_clusters:
                u = ad.obsm("X_umap")
                if u is not None and u.ndim == 2 and u.shape[1] >= 2 and len(u) == ad.n_obs:
                    umap = u[rows, :2]
            layers.append({"table": tn, "region": reg, "ids": inst, "X": X, "genes": genes, "where": where, "counts": counts,
                           "cols": cols, "colors": colors, "umap": umap})
        finally:
            ad.close()
    if not layers:
        raise SystemExit(f"[{tag}] no table annotates a shapes element" + (" with that --region" if args.region else ""))
    # ---- regions, transforms, the coordinate system, the image
    shapes = {}
    for L in layers:
        if L["region"] not in shapes:
            shapes[L["region"]] = read_shapes(el["shapes"][L["region"]], tag)
    tr = {n: ngff.element_transforms(sh["attrs"], f"shapes {log_name(n)}") for n, sh in shapes.items()}
    for n, t in tr.items():
        if not t:
            raise SystemExit(f"[{tag}] shapes {log_name(n)}: no coordinate transformation to any coordinate system")
    img_tr = {}
    for n, f in el["images"].items():
        a = ngff.read_attrs(f)
        try:
            img_tr[n] = ngff.element_transforms(a, f"image {log_name(n)}")
        except SystemExit as e:
            if args.image_element == n:
                raise
            log(f"[{tag}] image {log_name(n)} left out: {e}")
    cs = _pick_cs(args.coordinate_system, list(tr.values()), tag)
    image, img = None, None
    if args.image_element:
        if args.image_element not in el["images"]:
            raise SystemExit(f"[{tag}] --image-element {log_name(args.image_element)}: not one of the {len(el['images'])} images")
        if cs not in img_tr.get(args.image_element, {}):
            raise SystemExit(f"[{tag}] --image-element {log_name(args.image_element)}: not in the coordinate system {log_name(cs)}")
        image = args.image_element
    else:
        best = None
        for n in el["images"]:
            if cs not in img_tr.get(n, {}):
                continue
            try:
                im = ngff.NgffImage(el["images"][n], ngff.read_attrs(el["images"][n]), what=tag)
            except SystemExit as e:
                log(f"[{tag}] image {log_name(n)} left out: {e}")
                continue
            px = im.width * im.height
            if best is None or px > best[0]:
                best = (px, n, im)
        if best:
            image, img = best[1], best[2]
            others = [n for n in el["images"] if n != image]
            if others:
                log(f"[{tag}] image: {log_name(image)}, the largest in {log_name(cs)}; not used: "
                    f"{', '.join(log_name(o) for o in others)} (--image-element picks another)")
    if image is not None:
        if img is None:
            img = ngff.NgffImage(el["images"][image], ngff.read_attrs(el["images"][image]), what=tag)
        M_img = img_tr[image][cs]
        if abs(np.linalg.det(M_img[:2, :2])) < 1e-15:
            raise SystemExit(f"[{tag}] the image's transformation cannot be inverted")
        to_img = np.linalg.inv(M_img)
        P.source = img.source(z=args.z, log=log)
        P.colors = None if img.is_rgb else img.colors
        P.imaging = not img.is_rgb
        P.pixel_size_um = img.pixel_size_um
        P.manifest_source["image"] = str(el["images"][image])
        log(f"[{tag}] image {log_name(image)}: {img.width} x {img.height} px, {img.nC} channel(s)"
            + (f", z-plane {P.source.z} of {img.nZ}" if img.nZ > 1 else "")
            + (f"; pixel size {img.pixel_size_um:.4f} µm from its axes" if img.pixel_size_um else
               "; no pixel size in it" + (", --pixel-size used" if getattr(args, "pixel_size", None) else
                                          " (give --pixel-size for the scale bar and the molecules)")))
    else:
        to_img = None
        log(f"[{tag}] no image in {log_name(cs)}: a cell-density background is drawn")
    # every region's positions in the coordinate system, then in image pixels
    pos_cs = {}
    for n, sh in shapes.items():
        M = tr[n][cs] if cs in tr[n] else None
        if M is None:
            raise SystemExit(f"[{tag}] shapes {log_name(n)}: not in the coordinate system {log_name(cs)}")
        pos_cs[n] = ngff.apply(M, sh["xy"][:, 0], sh["xy"][:, 1])
    if to_img is None:
        allxy = np.concatenate([np.column_stack(v) for v in pos_cs.values()])
        spacing = median_spacing(allxy)
        P.source, bg_to_px, um, P.background = density_background(allxy, args.coords_um, spacing,
                                                                   args.background_um, "cell", log, tag)
        P.pixel_size_um = um
        # density_background's px from cs units: (x - x0) / step = an affine
        x0p, y0p = bg_to_px(0.0, 0.0)
        x1p, _ = bg_to_px(1.0, 0.0)
        s = float(x1p - x0p)
        to_img = np.array([[s, 0, float(x0p)], [0, s, float(y0p)], [0, 0, 1.0]])
        P.imaging = False
    um = P.pixel_size_um
    if getattr(args, "pixel_size", None) and image is not None:
        um = P.pixel_size_um = float(args.pixel_size)
    P.block = {"coordinateSystem": "global" if cs == "global" else "chosen", "tables": len(layers),
               "imageFrom": "images" if image is not None else "density background",
               "transforms": "element -> coordinate system -> image pixels (ngff.to_affine)"}
    # ---- layers (after the pyramid)
    used_keys = set()
    for L in layers:
        sh = shapes[L["region"]]
        M = to_img @ tr[L["region"]][cs]
        idx = {v: i for i, v in enumerate(sh["ids"])}
        hit = np.array([idx.get(i, -1) for i in L["ids"]], np.int64)
        miss = int((hit < 0).sum())
        if miss:
            log(f"[{tag}] WARNING {miss:,} of {len(hit):,} table rows have no shape in their region: no position")
        if miss == len(hit):
            raise SystemExit(f"[{tag}] table {log_name(L['table'])}: none of its instance ids is in its region: not the same data")
        px, py = ngff.apply(M, sh["xy"][np.maximum(hit, 0), 0], sh["xy"][np.maximum(hit, 0), 1])
        xy_px = np.column_stack([px, py])
        xy_px[hit < 0] = np.nan
        sc = ngff.scale_of(M)
        side = _squares(sh, M)
        if sh["kind"] == "circles":
            r = sh["radius"][np.maximum(hit, 0)][hit >= 0] * sc
            spot_px = float(np.median(2 * r)) if r.size else median_spacing(xy_px)
            spots = r.size and float(np.std(r) / max(np.mean(r), 1e-12)) < 0.01
            kind = None if spots else "cells"
        elif side is not None:
            spot_px, kind = side, None
        else:
            a = sh["area"][np.maximum(hit, 0)][hit >= 0] * sc * sc
            spot_px = float(np.median(2 * np.sqrt(a / np.pi))) if a.size else median_spacing(xy_px)
            kind = "cells"
        if kind == "cells":
            key, label, bin_um = "cells", "Cells", None
        elif side is not None:
            bin_um = side * um if um else None
            key, label = (f"{int(round(bin_um)):03d}um", f"{round(bin_um):g} µm bins") if bin_um else ("bins", "Bins")
        else:
            key, label, bin_um = "spots", "Spots", None
        k0, i = key, 2
        while key in used_keys:
            key = f"{k0}{i}"
            label = f"{label.split(' (')[0]} ({i})"
            i += 1
        used_keys.add(key)
        L.update(key=key, label=label, kind=kind, bin_um=bin_um, spot_px=spot_px, xy_px=xy_px)
        L["outlines"] = []
        if kind == "cells":
            L["outlines"] = _outline_sets(L, sh, shapes, el, tr, cs, to_img, tag, log)
    P.block["layers"] = [{"key": L["key"], "rows": len(L["ids"]), "counts": L["where"], "countsAreIntegers": L["counts"],
                          "region": "circles" if shapes[L["region"]]["kind"] == "circles" else "polygons"} for L in layers]
    if P.pixel_size_um is None and um:
        P.pixel_size_um = um

    def step(ds_dir, manifest, a, lg):
        mpp = manifest.get("pixelSizeUm")
        user = []
        for L in layers:
            meta = write_layer(ds_dir, L["key"], L["label"], L["xy_px"], L["X"], L["genes"], L["ids"], L["spot_px"], mpp,
                               kind=L["kind"], bin_um=L["bin_um"], extra={"platform": "spatialdata"}, log=lg, tag=tag)
            if L["outlines"] and not getattr(a, "no_outlines", False):
                import cellshapes
                cellshapes.attach(ds_dir, meta, L["outlines"], mpp, args=a, log=lg)
            manifest["expr"].append(meta)
            e = write_user_tables(ds_dir, L["key"], L["ids"], L["cols"], L["colors"], L["umap"], lg, tag)
            if e:
                user.append(e)
        if user:
            manifest["userObs"] = user
        if not a.no_transcripts:
            mol = _points(ds_dir, el, cs, to_img, (manifest["image"]["width"], manifest["image"]["height"]), mpp, a, lg, tag)
            if mol:
                manifest["molecules"] = mol
    P.steps.append(step)
    return P


def _norm_id(v: str) -> str:
    """An instance id as text, the same for 12, "12" and 12.0."""
    s = str(v).strip()
    try:
        f = float(s)
        if f.is_integer() and re.fullmatch(r"[-+]?\d+(\.0*)?", s):
            return str(int(f))
    except ValueError:
        pass
    return s


def _outline_sets(L, sh, shapes, el, tr, cs, to_img, tag, log) -> list:
    """Cell outlines of a cells layer: its region's polygons, and other polygon shapes with the same ids
    ("nucle..." in the name = nucleus outlines)."""
    import cellshapes
    sets = []
    want = set(L["ids"])

    def make(name, s, key):
        M = to_img @ ngff.element_transforms(s["attrs"], f"shapes {log_name(name)}")[cs]
        x, y = ngff.apply(M, s["ring_xy"][:, 0], s["ring_xy"][:, 1])
        return cellshapes.ShapeSet(key, s["start"], np.column_stack([x, y]), ids=[_norm_id(i) for i in s["ring_ids"]],
                                   source="SpatialData shapes")
    if sh["kind"] == "polygons":
        sets.append(make(L["region"], sh, "cell"))
    for name, f in el["shapes"].items():
        if name == L["region"]:
            continue
        try:
            s = shapes.get(name) or read_shapes(f, tag)
        except SystemExit as e:
            log(f"[{tag}] shapes {log_name(name)} not read for outlines: {e}")
            continue
        if s["kind"] != "polygons":
            continue
        t = ngff.element_transforms(s["attrs"], f"shapes {log_name(name)}")
        if cs not in t:
            continue
        ids = {_norm_id(i) for i in s["ids"]}
        if len(ids & want) < 0.5 * min(len(want), len(ids)):
            continue
        key = "nucleus" if "nucle" in name.lower() else "cell"
        if any(x.key == key for x in sets):
            log(f"[{tag}] shapes {log_name(name)}: a second set of {key} outlines, not used")
            continue
        sets.append(make(name, s, key))
        log(f"[{tag}] {key} outlines from shapes {log_name(name)} ({len(s['ids']):,} polygons)")
    sets.sort(key=lambda s: 0 if s.key == "cell" else 1)
    return sets


def _points(ds_dir, el, cs, to_img, size, mpp, args, log, tag):
    """The molecules layer from a points element (parquet parts), or None with the reason."""
    import pyarrow.parquet as pq
    import molecules
    pts = el["points"]
    if not pts:
        return None
    if args.points:
        if args.points not in pts:
            raise SystemExit(f"[{tag}] --points {log_name(args.points)}: not one of the {len(pts)} points elements")
        name = args.points
    elif len(pts) == 1:
        name = next(iter(pts))
    else:
        hits = [n for n in pts if "transcript" in n.lower()]
        if len(hits) != 1:
            log(f"[{tag}] {len(pts)} points elements and none (or several) named transcripts: no molecules layer "
                f"(--points picks one)")
            return None
        name = hits[0]
    if not mpp:
        log(f"[{tag}] no molecules layer: the molecules need the pixel size (give --pixel-size)")
        return None
    f = pts[name]
    attrs = ngff.read_attrs(f)
    t = ngff.element_transforms(attrs, f"points {log_name(name)}")
    if cs not in t:
        log(f"[{tag}] points {log_name(name)} are not in the coordinate system {log_name(cs)}: no molecules layer")
        return None
    M = to_img @ t[cs]
    sa = attrs.get("spatialdata_attrs") or {}
    fk = sa.get("feature_key")
    pq_dir = f / "points.parquet"
    parts = sorted(pq_dir.glob("*.parquet")) if pq_dir.is_dir() else ([pq_dir] if pq_dir.exists() else [])
    if not parts:
        log(f"[{tag}] points {log_name(name)}: no points.parquet, no molecules layer")
        return None
    cols = pq.ParquetFile(str(parts[0])).schema_arrow.names
    if not fk or fk not in cols or "x" not in cols or "y" not in cols:
        log(f"[{tag}] points {log_name(name)}: needs x, y and the feature_key column (spatialdata_attrs): no molecules layer")
        return None
    qv = "qv" if "qv" in cols else None
    t0 = time.time()
    names: dict[str, int] = {}
    # pass 1: the feature names (dictionary-encoded or plain text)
    for p in parts:
        for b in pq.ParquetFile(str(p)).iter_batches(batch_size=500_000, columns=[fk]):
            col = b.column(0)
            vals = col.dictionary.to_pylist() if hasattr(col, "dictionary") else col.unique().to_pylist()
            for v in vals:
                if v is not None and str(v) not in names:
                    names[str(v)] = len(names)
    gene_list = list(names)
    keep = np.array([not is_control(g) for g in gene_list])
    w = molecules.MoleculeWriter(gene_list, size, mpp, keep=keep, log=log)
    low_q = read = 0
    for p in parts:
        for b in pq.ParquetFile(str(p)).iter_batches(batch_size=500_000, columns=["x", "y", fk] + ([qv] if qv else [])):
            x = np.asarray(b.column(0).to_numpy(zero_copy_only=False), np.float64)
            y = np.asarray(b.column(1).to_numpy(zero_copy_only=False), np.float64)
            col = b.column(2)
            if hasattr(col, "dictionary"):
                d = np.array([names.get(str(v), -1) for v in col.dictionary.to_pylist()], np.int64)
                ix = np.asarray(col.indices.to_numpy(zero_copy_only=False))
                g = np.where(ix >= 0, d[np.maximum(ix, 0)], -1)
            else:
                g = np.array([names.get(str(v), -1) if v is not None else -1 for v in col.to_pylist()], np.int64)
            read += len(x)
            if qv:
                q = np.asarray(b.column(3).to_numpy(zero_copy_only=False), np.float64)
                ok = q >= QV_MIN
                low_q += int((~ok).sum())
                x, y, g = x[ok], y[ok], g[ok]
            px, py = ngff.apply(M, x, y)
            w.add(g, px * mpp, py * mpp)
    log(f"[{tag}] points {log_name(name)}: {read:,} read" + (f", {low_q:,} with qv < {QV_MIN:g} left out" if qv else "")
        + f", {int((~keep).sum())} control features left out ({time.time() - t0:.0f}s)")
    meta = {"platform": "spatialdata", "source": "points"}
    if qv:
        meta["qvMin"] = QV_MIN
    return w.finish(Path(ds_dir), meta=meta)


# ---------------------------------------------------------------------------------------
# prepare.py hooks
# ---------------------------------------------------------------------------------------
PLANS = {"omezarr": plan_omezarr, "anndata": plan_anndata, "spatialdata": plan_spatialdata}
ADD_ARGS = {"omezarr": add_args_omezarr, "anndata": add_args_anndata, "spatialdata": add_args_spatialdata}
HELP = {"omezarr": "an OME-Zarr image (OME-NGFF 0.4 or older, zarr v2): the image layer, as `tiff` does for an OME-TIFF",
        "anndata": "a spatial AnnData (.h5ad or AnnData-Zarr, Scanpy / Squidpy layout): obsm['spatial'] positions, counts, "
                   "the uns['spatial'] image, obs clusterings and X_umap as your own (stdformats.py)",
        "spatialdata": "a SpatialData store (.zarr, zarr v2): the image, every table on its shapes region (circles / "
                       "polygons -> outlines), points as molecules (stdformats.py)"}


def plan(kind: str, args, log=print) -> Plan:
    return PLANS[kind](args, log=log)
