"""Xenium Onboard Analysis output folder (or a zip of it) -> cells layer + molecules layer.

Reads, when present (each part is optional; a folder with stains only still works, and so does one with
cells and / or transcripts but no stain images):

    experiment.xenium                       pixel_size (µm per pixel of the morphology images)
    morphology_focus/*.ome.tif              the stain images; without them the dataset is drawn on a neutral
                                            density background (transcripts, else cells), see plan()
    cells.parquet | cells.csv.gz            cell_id, x_centroid, y_centroid (µm), cell_area
    cell_feature_matrix.h5 | cell_feature_matrix/ (barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz)
    transcripts.parquet | transcripts.zarr.zip | transcripts.zarr/

--src may be the folder, or a zip of it (Xenium's output bundle, or a zip of the whole output folder, with
the dataset folder at the top of the zip or a few folders down): zipinput.py reads the members in place
where that is cheap and copies out only the stain images (and a compressed transcripts.zarr.zip) into a
work folder beside the cache (<out>/_unzipped/<id>/, --unzip-to), never beside the zip.

Xenium coordinates are µm in the frame of the morphology images, origin at the top-left
corner of the full-resolution image, so pixel = µm / pixel_size (the same transform Xenium
Explorer and FIVE Atlas use; no affine). Without images the background keeps that frame: its pixels
are a whole number of morphology pixels (k x pixel_size, about --background-um), origin at Xenium's 0, 0.

Cells become an ordinary expression layer (expr.build_cell_layer), so the whole-section
UMAP, Leiden clusters, suggested names, bin QC, circle -> UMAP and the volcano panel run on
them unchanged. Transcripts become the molecules layer (molecules.py): quality value >= 20
(Xenium's own "high quality" cut, used for its cell-feature matrix), control probes and
codewords left out.
"""
from __future__ import annotations

import ast
import gzip
import io
import json
import math
import time
from pathlib import Path

import numpy as np

QV_MIN = 20.0
# names of codewords / probes that are not genes, in every Xenium version seen so far
CONTROL_PREFIXES = ("negcontrolprobe", "negcontrolcodeword", "blank", "unassignedcodeword",
                    "deprecatedcodeword", "intergenic", "antisense", "genomiccontrol")
GENE_FEATURE_TYPES = ("gene expression",)
DEFAULT_PIXEL_SIZE_UM = 0.2125      # Xenium's morphology images (every Xenium Onboard Analysis version so far)
BACKGROUND_UM = 1.0                 # no stain images: about this many µm per background pixel
MARGIN_UM = 10.0
# names that make a folder (or a folder inside a zip) a Xenium output folder
ROOT_NAMES = {"experiment.xenium", "morphology_focus", "transcripts.parquet", "transcripts.zarr.zip", "transcripts.zarr",
              "cells.parquet", "cells.csv.gz", "cell_feature_matrix.h5", "cell_feature_matrix", "cells.zarr.zip"}


def is_control(name: str) -> bool:
    return str(name).strip().lower().startswith(CONTROL_PREFIXES)


def is_root(names) -> bool:
    return bool(set(names) & ROOT_NAMES)


# A dataset folder is a Path or a zipinput.ZipInput; a file in it a Path or a zipinput.Member (both have
# .name and .suffix). These helpers read either.
def _child(folder, name: str):
    if isinstance(folder, Path):
        return folder / name
    from zipinput import Member
    if isinstance(folder, Member):                     # a folder inside the zip
        return folder.zi.member(f"{folder.rel}/{name}")
    return folder.member(name)


def _exists(p) -> bool:
    return p is not None and (p.exists() if isinstance(p, Path) else True)


def _is_dir(p) -> bool:
    return p is not None and p.is_dir()


def _first(folder, *names):
    for n in names:
        p = _child(folder, n)
        if _exists(p):
            return p
    return None


def _rand(p):
    """Something pyarrow / h5py can seek in."""
    return str(p) if isinstance(p, Path) else p.random_access()


def _stream(p):
    return open(p, "rb") if isinstance(p, Path) else p.open()


def _read_text(p) -> str:
    return p.read_text(encoding="utf-8") if isinstance(p, Path) else p.read_bytes().decode("utf-8")


def find_inputs(folder) -> dict:
    """The Xenium files present under `folder` (a Path or a zipinput.ZipInput; None for each one that is missing)."""
    folder = Path(folder) if isinstance(folder, str) else folder
    mtx_dir = _child(folder, "cell_feature_matrix")
    has_mtx = _is_dir(mtx_dir) and (_first(mtx_dir, "matrix.mtx.gz", "matrix.mtx") is not None)
    return {
        "experiment": _first(folder, "experiment.xenium"),
        "cells": _first(folder, "cells.parquet", "cells.csv.gz", "cells.csv"),
        "matrix_h5": _first(folder, "cell_feature_matrix.h5"),
        "matrix_mtx": mtx_dir if has_mtx else None,
        "transcripts": _first(folder, "transcripts.parquet", "transcripts.zarr.zip", "transcripts.zarr"),
        "cell_boundaries": _first(folder, "cell_boundaries.parquet", "cell_boundaries.csv.gz"),
    }


def pixel_size(folder) -> float | None:
    p = _first(Path(folder) if isinstance(folder, str) else folder, "experiment.xenium")
    if p is None:
        return None
    try:
        v = json.loads(_read_text(p)).get("pixel_size")
        return float(v) if v else None
    except (ValueError, TypeError, OSError):
        return None


# ---------------------------------------------------------------------------------- cells
def _table(path, columns: list[str]) -> dict:
    """Columns of a CSV (optionally gzipped) or Parquet table as numpy / lists."""
    if path.suffix == ".parquet":
        import pyarrow.parquet as pq
        pf = pq.ParquetFile(_rand(path))
        have = set(pf.schema_arrow.names)
        t = pf.read(columns=[c for c in columns if c in have])
    else:
        import pyarrow as pa
        import pyarrow.csv as pcsv
        if isinstance(path, Path):
            t = pcsv.read_csv(str(path))
        else:                                             # a member of a zip, read in order
            with path.open() as f:
                t = pcsv.read_csv(pa.input_stream(f, compression="gzip" if path.suffix == ".gz" else None))
    out = {}
    for c in columns:
        if c in t.column_names:
            col = t.column(c)
            out[c] = col.to_pylist() if c == "cell_id" else col.to_numpy()
    return out


def _cell_ids(v) -> list[str]:
    return [x.decode() if isinstance(x, bytes) else str(x) for x in v]


def read_matrix(inp: dict):
    """(X cells x features CSR float32, barcodes, feature names, feature types)."""
    from scipy import io as sio, sparse
    if inp.get("matrix_h5"):
        import h5py
        with h5py.File(_rand(inp["matrix_h5"]), "r") as f:
            m = f["matrix"]
            G, B = (int(v) for v in m["shape"][()])
            X = sparse.csc_matrix((m["data"][()].astype(np.float32), m["indices"][()], m["indptr"][()]),
                                  shape=(G, B)).T.tocsr()
            bc = _cell_ids(m["barcodes"][()])
            names = _cell_ids(m["features/name"][()])
            types = _cell_ids(m["features/feature_type"][()])
        return X, bc, names, types
    d = inp["matrix_mtx"]

    def opened(p):
        raw = _stream(p)
        return gzip.GzipFile(fileobj=raw, mode="rb") if p.suffix == ".gz" else raw

    def lines(stem):
        with io.TextIOWrapper(opened(_first(d, stem + ".gz", stem)), encoding="utf-8") as f:
            return [ln.rstrip("\n\r").split("\t") for ln in f if ln.strip()]
    bc = [r[0] for r in lines("barcodes.tsv")]
    feats = lines("features.tsv")
    names = [r[1] if len(r) > 1 else r[0] for r in feats]
    types = [r[2] if len(r) > 2 else "Gene Expression" for r in feats]
    mp = _first(d, "matrix.mtx.gz", "matrix.mtx")
    if isinstance(mp, Path):
        M = sio.mmread(str(mp))                                      # features x cells
    else:
        with opened(mp) as f:
            M = sio.mmread(io.BytesIO(f.read()))
    X = sparse.csr_matrix(M.T, dtype=np.float32)
    return X, bc, names, types


def read_cells(inp: dict, log=print) -> dict | None:
    """Cell centroids (µm) in count-matrix order + the gene matrix, or None when either is missing."""
    if not inp.get("cells") or not (inp.get("matrix_h5") or inp.get("matrix_mtx")):
        missing = (([] if inp.get("cells") else ["cells.parquet / cells.csv.gz"]) +
                   ([] if (inp.get("matrix_h5") or inp.get("matrix_mtx")) else ["cell_feature_matrix"]))
        log(f"[xenium] no cells layer: {' and '.join(missing)} not found")
        return None
    t0 = time.time()
    X, bc, names, types = read_matrix(inp)
    keep = np.array([t.strip().lower() in GENE_FEATURE_TYPES and not is_control(n) for n, t in zip(names, types)])
    n_ctrl = int((~keep).sum())
    X = X[:, np.flatnonzero(keep)].tocsr()
    genes = [n for n, k in zip(names, keep) if k]
    tab = _table(inp["cells"], ["cell_id", "x_centroid", "y_centroid", "cell_area"])
    ids = _cell_ids(tab["cell_id"])
    row = {c: i for i, c in enumerate(ids)}
    xy = np.full((len(bc), 2), np.nan, np.float64)
    area = np.full(len(bc), np.nan, np.float64)
    hit = np.array([row.get(b, -1) for b in bc], np.int64)
    ok = hit >= 0
    xy[ok, 0] = np.asarray(tab["x_centroid"], np.float64)[hit[ok]]
    xy[ok, 1] = np.asarray(tab["y_centroid"], np.float64)[hit[ok]]
    if "cell_area" in tab:
        area[ok] = np.asarray(tab["cell_area"], np.float64)[hit[ok]]
    log(f"[xenium] cells: {len(bc):,} in the matrix, {len(ids):,} in {inp['cells'].name}, "
        f"{int(ok.sum()):,} matched; {len(genes)} genes ({n_ctrl} control features left out), "
        f"nnz {X.nnz:,}  ({time.time() - t0:.0f}s)")
    if int((~ok).sum()):
        log(f"[xenium]   WARNING {int((~ok).sum())} matrix cells have no centroid: no position")
    med_area = float(np.nanmedian(area)) if np.isfinite(area).any() else None
    return {"X": X, "genes": genes, "xy_um": xy, "barcodes": list(bc), "nCellsTable": len(ids), "medianAreaUm2": med_area,
            "nControlFeatures": n_ctrl}


# ---------------------------------------------------------------------------- transcripts
def _attr_list(v):
    return v if isinstance(v, (list, dict)) else ast.literal_eval(v)


def _open_zarr(path, log=print):
    """(root group, a function that closes the store)."""
    import zarr
    if not isinstance(path, Path):                     # inside a zip (zipinput.py)
        import zipinput
        store, close = zipinput.zarr_store(path, log=log)
        return zarr.open(store, mode="r"), close
    if path.is_dir():
        return zarr.open(str(path), mode="r"), (lambda: None)
    store = zarr.ZipStore(str(path), mode="r")
    return zarr.open(store, mode="r"), store.close


def zarr_genes(root) -> tuple[list[str], np.ndarray]:
    """Gene names indexed by gene_identity, and which of them are genes (not controls)."""
    names = [str(n) for n in _attr_list(root.attrs["gene_names"])]
    keep = np.array([not is_control(n) for n in names])
    if "gene_category" in root:                                   # (genes x flags), flag 0 = is_gene
        try:
            cols = _attr_list(root["gene_category"].attrs.get("column_names", []))
            cat = np.asarray(root["gene_category"][:])
            if cat.ndim == 2 and cat.shape[0] == len(names) and "is_gene" in cols:
                keep &= cat[:, cols.index("is_gene")].astype(bool)
        except Exception:
            pass
    return names, keep


def zarr_batches(root, qv_min: float, stats: dict):
    """(gene_identity, x µm, y µm) per level-0 grid chunk, quality value >= qv_min."""
    g0 = root["grids/0"]
    keys = list(g0.group_keys())
    for k in keys:
        c = g0[k]
        if "gene_identity" not in c:
            continue
        gid = np.asarray(c["gene_identity"][:]).reshape(-1)
        if not len(gid):
            continue
        loc = np.asarray(c["location"][:])
        q = np.asarray(c["quality_score"][:]).reshape(-1)
        ok = q >= qv_min
        if "valid" in c:
            ok &= np.asarray(c["valid"][:]).reshape(-1).astype(bool)
        stats["read"] += len(gid)
        stats["lowQv"] += int((q < qv_min).sum())
        yield gid[ok], loc[ok, 0], loc[ok, 1]


def parquet_batches(path, names: list[str], qv_min: float, stats: dict):
    import pyarrow.parquet as pq
    pf = pq.ParquetFile(str(path) if isinstance(path, Path) else path)
    idx = {n: i for i, n in enumerate(names)}
    cols = ["feature_name", "x_location", "y_location", "qv"]
    for b in pf.iter_batches(batch_size=4_000_000, columns=cols):
        fn = b.column(0).to_pylist()
        x = b.column(1).to_numpy()
        y = b.column(2).to_numpy()
        q = b.column(3).to_numpy()
        fn = [f.decode() if isinstance(f, bytes) else str(f) for f in fn]
        gid = np.fromiter((idx.get(f, -1) for f in fn), np.int64, len(fn))
        ok = q >= qv_min
        stats["read"] += len(gid)
        stats["lowQv"] += int((~ok).sum())
        yield gid[ok], x[ok], y[ok]


def parquet_genes(path) -> list[str]:
    import pyarrow.compute as pc
    import pyarrow.parquet as pq
    t = pq.ParquetFile(str(path) if isinstance(path, Path) else path).read(columns=["feature_name"])
    u = pc.unique(t.column(0)).to_pylist()
    return sorted(f.decode() if isinstance(f, bytes) else str(f) for f in u)


def open_transcripts(inp: dict, qv_min: float, stats: dict, log=print):
    """(gene names, keep mask, batches of (gene id, x µm, y µm), kind, close)."""
    src = inp["transcripts"]
    if src.suffix == ".parquet":
        f = _rand(src)
        names = parquet_genes(f)
        if hasattr(f, "seek"):
            f.seek(0)
        keep = np.array([not is_control(n) for n in names])
        return names, keep, parquet_batches(f, names, qv_min, stats), "transcripts.parquet", (lambda: None)
    root, close = _open_zarr(src, log=log)
    names, keep = zarr_genes(root)
    return names, keep, zarr_batches(root, qv_min, stats), "transcripts.zarr", close


def transcript_extent(inp: dict, log=print) -> tuple[float, float] | None:
    """(x max, y max) in µm of the transcripts, without reading them all: the zarr's grid of chunks (an
    upper bound, at most one chunk too far), or the parquet's column statistics (x / y read when it has none)."""
    src = inp.get("transcripts")
    if not src:
        return None
    if src.suffix == ".parquet":
        import pyarrow.compute as pc
        import pyarrow.parquet as pq
        pf = pq.ParquetFile(_rand(src))
        names = pf.schema_arrow.names
        mx = []
        for col in ("x_location", "y_location"):
            j = names.index(col)
            st = [pf.metadata.row_group(r).column(j).statistics for r in range(pf.metadata.num_row_groups)]
            if st and all(x is not None and x.has_min_max for x in st):
                mx.append(max(float(x.max) for x in st))
            else:
                mx.append(float(pc.max(pf.read(columns=[col]).column(0)).as_py()))
        return mx[0], mx[1]
    root, close = _open_zarr(src, log=log)
    try:
        size = _attr_list(root["grids"].attrs.get("grid_size", [250.0])) if "grids" in root else [250.0]
        size = float(size[0] if isinstance(size, (list, tuple)) else size)
        keys = [k for k in root["grids/0"].group_keys() if "," in k]
        if not keys:
            return None
        ij = np.array([[int(v) for v in k.split(",")[:2]] for k in keys])
        return float((ij[:, 0].max() + 1) * size), float((ij[:, 1].max() + 1) * size)
    finally:
        close()


def write_transcripts(inp: dict, ds_dir: Path, image_size, pixel_size_um: float, qv_min: float = QV_MIN,
                      log=print, density=None, finish: bool = True):
    """molecules/ from transcripts.parquet or transcripts.zarr(.zip); returns the manifest block. With
    `density` (flatfiles.Density over the same pixels) every molecule kept is also counted there; with
    finish=False the writer is returned unfinished (finish_transcripts writes it after the pyramid)."""
    import molecules
    src = inp.get("transcripts")
    if not src:
        log("[xenium] no molecules layer: transcripts.parquet / transcripts.zarr.zip not found")
        return None
    t0 = time.time()
    stats = {"read": 0, "lowQv": 0}
    close = None
    try:
        names, keep, batches, kind, close = open_transcripts(inp, qv_min, stats, log=log)
        log(f"[xenium] transcripts from {kind}: {int(keep.sum())} genes, {int((~keep).sum())} control codewords "
            f"left out, quality >= {qv_min:g}")
        w = molecules.MoleculeWriter(names, image_size, pixel_size_um, keep=keep, log=log)
        for i, (g, x, y) in enumerate(batches):
            w.add(g, x, y)
            if density is not None:
                gi = np.asarray(g).astype(np.int64, copy=False).ravel()
                k = (gi >= 0) & (gi < len(keep))
                k[k] = keep[gi[k]]
                density.add(np.asarray(x, np.float64).ravel()[k] / pixel_size_um,
                            np.asarray(y, np.float64).ravel()[k] / pixel_size_um)
            if (i + 1) % 200 == 0:
                log(f"[xenium]   {stats['read']:,} transcripts read  ({time.time() - t0:.0f}s)")
    finally:
        if close is not None:
            close()
    log(f"[xenium] transcripts: {stats['read']:,} read, {stats['lowQv']:,} below quality {qv_min:g}  "
        f"({time.time() - t0:.0f}s)")
    w.xenium_meta = {"platform": "xenium", "source": kind, "qvMin": qv_min, "read": int(stats["read"]),
                     "lowQv": int(stats["lowQv"])}
    return finish_transcripts(w, ds_dir) if finish else w


def finish_transcripts(w, ds_dir: Path) -> dict:
    return w.finish(ds_dir, meta=w.xenium_meta)


# ---------------------------------------------------------------------------------- plan
def add_args(p) -> None:
    p.add_argument("--src", required=True, help="the Xenium output folder (the one holding morphology_focus/, or "
                                                "cells / transcripts without stain images), or a zip of it")
    p.add_argument("--no-cells", action="store_true", help="skip the cells layer")
    p.add_argument("--no-transcripts", action="store_true", help="skip the molecules layer")
    p.add_argument("--qv-min", type=float, default=20.0,
                   help="molecules: keep transcripts with quality value >= this (default 20, Xenium's own cut)")
    p.add_argument("--background-um", type=float,
                   help=f"no stain images: about this many µm per pixel of the density background (default "
                        f"{BACKGROUND_UM:g}; rounded to a whole number of morphology pixels)")
    p.add_argument("--unzip-to", help="a zipped --src: the work folder for the files copied out of the zip (default "
                                      "<out>/_unzipped/<id>; never next to the zip)")
    p.add_argument("--keep-unzipped", action="store_true",
                   help="a zipped --src: keep the files copied out of the zip after a successful prepare (default: "
                        "they are removed; the zip itself is never changed)")


def resolve_folder(src: Path) -> Path:
    """The Xenium folder for --src: itself, or its data/ folder (an archive's <name>/data/ layout)."""
    try:
        here = {p.name for p in src.iterdir()}
    except OSError:
        return src
    if "morphology_focus" not in here and (src / "data" / "morphology_focus").is_dir():
        return src / "data"
    if not is_root(here) and (src / "data").is_dir() and is_root(p.name for p in (src / "data").iterdir()):
        return src / "data"
    return src


class Plan:
    """What prepare.py needs: the image source and pixel size now, the layers after the pyramid."""

    def __init__(self):
        self.source = None
        self.pixel_size_um = None
        self.manifest_source = {}
        self.background = None
        self.zip = None
        self.folder = None
        self.inp = {}
        self.cells = None
        self._writer = None
        self.xenium_px = None

    def finish(self, ds_dir: Path, manifest: dict, args, log=print):
        from expr import build_cell_layer
        import cellshapes
        um = manifest["pixelSizeUm"]
        size = (manifest["image"]["width"], manifest["image"]["height"])
        if self.background:
            manifest["image"]["background"] = self.background
            manifest["image"]["channels"][0]["color"] = "#b4b4b4"
            manifest["frame"] = {"from": "Xenium x, y (µm)", "affine": [[1.0 / um, 0.0, 0.0], [0.0, 1.0 / um, 0.0]],
                                 "xeniumPixelSizeUm": self.xenium_px}
        if not um:
            log("[xenium] WARNING no pixel size (experiment.xenium / OME): cells and transcripts skipped; give --pixel-size")
            return
        cells = self.cells if self.background else (None if args.no_cells else read_cells(self.inp, log=log))
        if cells:
            area = cells["medianAreaUm2"]
            diam = 2.0 * (area / math.pi) ** 0.5 if area else None
            meta = build_cell_layer(ds_dir, cells["xy_um"], um, cells["X"], cells["genes"], key="cells",
                                    label="Cells", cell_diameter_um=diam, log=log,
                                    extra={"platform": "xenium"}, barcodes=cells["barcodes"])
            if not args.no_outlines:
                cellshapes.attach(ds_dir, meta, outline_sets(self.folder, um, log=log), um, args=args, log=log)
            manifest["expr"].append(meta)
        if args.no_transcripts:
            return
        if self._writer is not None:
            manifest["molecules"] = finish_transcripts(self._writer, ds_dir)
        elif not self.background:
            mol = write_transcripts(self.inp, ds_dir, size, um, qv_min=args.qv_min, log=log)
            if mol:
                manifest["molecules"] = mol

    def cleanup(self, args, log=print):
        if self.zip is not None:
            self.zip.cleanup(keep=getattr(args, "keep_unzipped", False))
            self.zip.close()


def outline_sets(folder, um: float, log=print) -> list:
    """cellshapes.xenium_sets for a folder, or for a zip (the boundary files it would read copied out first)."""
    import cellshapes
    if isinstance(folder, Path):
        return cellshapes.xenium_sets(folder, um, log=log)
    want = []
    for key in ("cell", "nucleus"):
        m = _first(folder, f"{key}_boundaries.parquet", f"{key}_boundaries.csv.gz", f"{key}_boundaries.csv")
        if m is not None:
            want.append(m.rel)
    if len(want) < 2:
        m = _first(folder, "cells.zarr.zip")
        if m is not None and not m.is_dir():
            want.append(m.rel)
    if want:
        folder.extract(want, why="the cell outlines")
    return cellshapes.xenium_sets(folder.work, um, log=log)


def open_zip(path: Path, work: Path, log=print):
    import zipinput
    zi = zipinput.ZipInput(path, work, log=log, root_test=is_root, tag="xenium")
    log(f"[xenium] reading the zip {path.name} as it is ({zipinput._fmt_gb(path.stat().st_size)}; {len(zi._files):,} "
        f"files in the dataset folder, which is {zi.where()})")
    return zi


def plan(args, log=print) -> Plan:
    """Find the files (folder or zip), the stain images or a density background, and the pixel size."""
    import pyramid
    P = Plan()
    src = Path(args.src)
    if src.is_file():
        work = Path(args.unzip_to) if getattr(args, "unzip_to", None) else Path(args.out) / "_unzipped" / args.id
        P.zip = folder = open_zip(src, work, log=log)
        P.manifest_source = {"zip": str(src), "folder": str(src)}
    elif src.is_dir():
        folder = resolve_folder(src)
        P.manifest_source = {"folder": str(folder)}
    else:
        raise SystemExit(f"[xenium] --src {src} not found")
    P.folder = folder
    P.inp = inp = find_inputs(folder)
    if isinstance(folder, Path):
        stains = sorted((folder / "morphology_focus").glob("*.ome.tif*"))
    else:
        stains = folder.glob("morphology_focus", r"\.ome\.tiff?$")
    xpx = pixel_size(folder)

    if stains:
        import re as _re
        import tifffile as _tf
        if not isinstance(folder, Path):               # random access into big compressed TIFFs: copy them out
            P.manifest_source["files"] = [m.rel for m in stains]       # in the zip's dataset folder (copies go)
            stains = folder.extract([m.rel for m in stains], why="the stain images, read tile by tile")
        with _tf.TiffFile(str(stains[0])) as tf0:           # the OME block lists every channel
            ome = tf0.ome_metadata or ""
        names = _re.findall(r'<Channel[^>]*?\sName="([^"]*)"', ome)
        px = _re.findall(r'PhysicalSizeX="([^"]*)"', ome)
        pairs = []
        for f in stains:
            m = _re.match(r"ch(\d+)_", f.name)
            idx = int(m.group(1)) if m else len(pairs)
            pairs.append((names[idx] if idx < len(names) else f.stem.split("_", 1)[-1], f))
        P.source = pyramid.TiffFileSource(pairs)
        um = args.pixel_size
        if um is None and px:
            um = float(px[0])
        if args.pixel_size is None and xpx:
            if um and abs(um - xpx) / xpx > 0.01:
                log(f"[xenium] pixel size: experiment.xenium says {xpx:.4f} µm/px, the OME metadata {um:.4f}; "
                    f"using experiment.xenium's (cells and transcripts are placed with it)")
            um = xpx
        P.pixel_size_um = um
        ann = _first(folder, "Custom_annotation.geojson")
        if ann is not None and not args.geojson:
            args.geojson = [f"Annotations={ann if isinstance(ann, Path) else ann.local()}"]
        P.manifest_source.setdefault("files", [str(p) for _, p in pairs])
        return P

    # ---- no stain images: a neutral density background in Xenium's own frame
    if args.pixel_size:
        xpx = float(args.pixel_size)
        log(f"[xenium] morphology pixel size {xpx:g} µm (--pixel-size)")
    elif xpx:
        log(f"[xenium] morphology pixel size {xpx:g} µm (experiment.xenium)")
    else:
        xpx = DEFAULT_PIXEL_SIZE_UM
        log(f"[xenium] no experiment.xenium: the morphology pixel size is taken as {xpx:g} µm, Xenium's default "
            f"(--pixel-size sets it); it sets the background's pixel grid only, as positions are in µm")
    P.xenium_px = xpx
    tx = inp.get("transcripts") if not args.no_transcripts else None
    P.cells = None if args.no_cells else read_cells(inp, log=log)
    if not tx and P.cells is None:
        raise SystemExit(f"[xenium] no stain images (morphology_focus/*.ome.tif), and no transcripts or cells "
                         f"(cells table + cell_feature_matrix) to draw instead"
                         + (" (--no-cells / --no-transcripts left them out)" if args.no_cells or args.no_transcripts else "")
                         + ": nothing to prepare")
    xmax = ymax = 0.0
    if tx:
        e = transcript_extent(inp, log=log)
        if e:
            xmax, ymax = max(xmax, e[0]), max(ymax, e[1])
    if P.cells is not None:
        xy = P.cells["xy_um"]
        if np.isfinite(xy).any():
            xmax, ymax = max(xmax, float(np.nanmax(xy[:, 0]))), max(ymax, float(np.nanmax(xy[:, 1])))
    if not (xmax > 0 and ymax > 0):
        raise SystemExit("[xenium] no positions to draw a background from")
    import flatfiles as ff
    want = getattr(args, "background_um", None) or BACKGROUND_UM
    k = max(1, int(round(want / xpx)))
    while ((xmax + MARGIN_UM) / (k * xpx)) * ((ymax + MARGIN_UM) / (k * xpx)) > ff.BACKGROUND_MAX_PX:
        k *= 2
    bg = k * xpx
    W, H = max(2, int(math.ceil((xmax + MARGIN_UM) / bg))), max(2, int(math.ceil((ymax + MARGIN_UM) / bg)))
    dens = ff.Density(W, H)
    P.pixel_size_um = bg
    if tx:
        P._writer = write_transcripts(inp, None, (W, H), bg, qv_min=args.qv_min, log=log, density=dens, finish=False)
        label, sigma, P.background = "Transcript density (no image)", 1.5 / bg, "transcript density"
    else:
        xy = P.cells["xy_um"]
        ok = np.isfinite(xy).all(axis=1)
        dens.add(xy[ok, 0] / bg, xy[ok, 1] / bg)
        label, sigma, P.background = "Cell density (no image)", 5.0 / bg, "cell density"
    log(f"[xenium] no stain images: background = {P.background}, {W} x {H} px of {bg:.4g} µm ({k} x {xpx:g} µm "
        f"morphology pixels), origin at Xenium's 0, 0 ({dens.n:,} points)")
    P.source = ff.ArraySource([label], dens.image(sigma))
    P.manifest_source["images"] = None
    return P
