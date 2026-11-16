r"""What kind of dataset is in this folder?

    python prepare.py detect <folder or file> [--json] [--kind <kind>]
    python prepare.py auto --src <folder or file> --id <id> --out <cache> [--kind <kind>] [options of that kind]

`detect()` looks at a folder (or one file) and says which `prepare.py` kind it holds, how sure
that is, what was found, what is missing and the exact arguments `prepare.py` would get. It
only reads names, sizes and small headers (a TIFF's first directory, a PNG's size, the first
lines of a CSV / sheet, the key names of an HDF5 file, a zip's file list). It never loads a
matrix or an image and never writes anything.

The signatures are what each reader needs (prepare.py, xenium.py, cosmx.py, merscope.py,
stereoseq.py, geomx.py, rnaseq.py, expr.find_visium_layers):

    visium     a Space Ranger folder: binned_outputs/square_NNNum/ (Visium HD) or spatial/ with
               tissue positions (classic), each layer with filtered_feature_bc_matrix.h5; plus the
               full-resolution image given to spaceranger, which is NOT in that folder: a TIFF
               next to it (its pixel size is checked against tissue_hires_image.png / scale factor)
    xenium     morphology_focus/*.ome.tif (required); experiment.xenium, cells + cell_feature_matrix,
               transcripts.parquet / .zarr.zip (each optional)
    cosmx      *_exprMat_file, *_metadata_file, *_fov_positions_file, *_tx_file (csv / csv.gz /
               parquet) and the FOV image folders (Morphology2D, CellComposite ...)
    merscope   detected_transcripts, cell_by_gene + cell_metadata, images/mosaic_<stain>_z<k>.tif
               with micron_to_mosaic_pixel_transform.csv
    stereoseq  <SN>_<stain>_regist.tif (required) + <SN>.tissue.gef / .raw.gef / .gem(.gz) /
               <SN>.cellbin.gef
    geomx      *.dcc files (or a zip of them), *.pkc, an annotation sheet with a Sample_ID column
    akoya      a folder of one TIFF per marker (NN_Marker_..._gray.tif)
    tiff       one image file: QPTIFF, OME-TIFF, a TIFF stack, an RGB (H&E) image, one channel
    rnaseq     an .h5ad with a UMAP in obsm, or a 10x matrix folder + a UMAP CSV

Verdicts (the rule: when it is not sure it does not guess, it shows what it found):

    certain     one kind fits, every file that kind requires is there, no choice is open
    ambiguous   more than one kind fits, or a choice is open (two candidate images, two
                segmentation runs, several datasets in one folder ...)
    incomplete  a kind is recognised but a required file is missing (it says which), or the
                data is still zipped
    unknown     nothing recognised, an empty folder, a path that does not exist or cannot be read

Files of a platform win over loose image files: a Xenium folder that also holds derived TIFFs
is a Xenium dataset, and the TIFFs are listed as "not used".
"""
from __future__ import annotations

import json
import os
import re
import time
import zipfile
from pathlib import Path

PLATFORM_KINDS = ("visium", "xenium", "cosmx", "merscope", "stereoseq", "geomx")
GENERIC_KINDS = ("akoya", "tiff", "rnaseq")
KINDS = PLATFORM_KINDS + GENERIC_KINDS
# other names people (and the notebook's older input folders) use for a kind
ALIASES = {"proteomics_file": "tiff", "proteomics_markers": "akoya", "scrnaseq": "rnaseq", "qptiff": "tiff",
           "ome-tiff": "tiff", "ometiff": "tiff", "he": "tiff", "h&e": "tiff", "visiumhd": "visium",
           "visium-hd": "visium", "visium_hd": "visium", "h5ad": "rnaseq", "10x": "rnaseq", "codex": "tiff",
           "phenocycler": "tiff", "stereo-seq": "stereoseq", "saw": "stereoseq"}
KIND_WORDS = {"visium": "Visium (Space Ranger output)", "xenium": "Xenium output folder", "cosmx": "CosMx flat files",
              "merscope": "MERSCOPE region folder", "stereoseq": "Stereo-seq (SAW output)", "geomx": "GeoMx DSP study",
              "akoya": "one TIFF per marker", "tiff": "one image file", "rnaseq": "single-cell RNA-seq"}
NEEDS = {
    "visium": "Space Ranger's outs folder (binned_outputs/square_NNNum/ for Visium HD, or filtered_feature_bc_matrix.h5 "
              "+ spatial/tissue_positions for classic Visium) and the full-resolution image given to spaceranger",
    "xenium": "a Xenium output folder with morphology_focus/*.ome.tif (cells and transcripts are optional)",
    "cosmx": "CosMx flat files: *_exprMat_file + *_metadata_file and / or *_tx_file (csv, csv.gz or parquet)",
    "merscope": "a MERSCOPE region folder: detected_transcripts, cell_by_gene + cell_metadata and / or images/mosaic_*.tif",
    "stereoseq": "the registered image <SN>_<stain>_regist.tif and a .gef / .gem(.gz) / .cellbin.gef file",
    "geomx": "*.dcc files (or a zip of them), a *.pkc probe kit file and an annotation sheet with a Sample_ID column",
    "akoya": "a folder of TIFF files, one per marker, all the same size",
    "tiff": "one image file (.qptiff, .ome.tif, .tif, .tiff, .btf)",
    "rnaseq": "an .h5ad with a UMAP in obsm, or a 10x matrix folder (matrix.mtx, features.tsv, barcodes.tsv) with a UMAP CSV",
}
# the options of each kind that detection fills in (in the order prepare.py gets them)
KIND_ARGS = {"visium": ("image", "outs"), "xenium": ("src",), "cosmx": ("src",), "merscope": ("src", "cells_dir"),
             "stereoseq": ("image", "gef", "gem", "cellbin", "tissue_mask"),
             "geomx": ("dcc", "pkc", "annotation", "image"), "akoya": ("src",), "tiff": ("src",),
             "rnaseq": ("h5ad", "mtx", "umap", "obs")}
_REPEAT = {"annotation", "image", "obs"}      # given once per value: --annotation a --annotation b
_NARGS = {"pkc"}                              # given once with every value: --pkc a b
# options that are not files but change what is required
OVERRIDE_FLAGS = ("images", "embedding")

MAX_DEPTH = 8                  # folders below the one given
MAX_ENTRIES = 200_000          # files + folders listed before giving up
MAX_SECONDS = 20.0             # listing time before giving up (slow network shares)
MAX_ZIPS = 12                  # zip file lists read per call
MAX_PEEKS = 60                 # TIFF headers opened per call
TIFF_EXT = (".tif", ".tiff", ".btf", ".qptiff")
_SKIP_DIRS = {"__pycache__", "node_modules", "$recycle.bin", "system volume information"}
_TABLE_EXT = r"\.(?:parquet|csv\.gz|csv)$"
_SQUARE = re.compile(r"^square_(\d+)um$")
_POSITIONS = ("tissue_positions.parquet", "tissue_positions.csv", "tissue_positions_list.csv")
_MOSAIC = re.compile(r"^mosaic_(.+)_z(\d+)\.tiff?$", re.I)                     # merscope.py
_CX_TABLE = re.compile(r"^(.*?)(exprmat_file|metadata_file|fov_positions_file|tx_file)" + _TABLE_EXT, re.I)
_CX_ROLE = {"exprmat_file": "matrix", "metadata_file": "cells", "fov_positions_file": "fov_positions", "tx_file": "transcripts"}
_CX_IMG = {"morphology2d": "morphology", "rawmorphologyimages": "zstack", "resegmentationimages": "zstack",
           "cellcomposite": "composite", "celllabels": "labels"}                # cosmx.IMAGE_DIRS
_CX_OTHER = ("compartmentlabels", "celloverlay")
_CX_FOV = re.compile(r"(?:^|_)F(\d{1,5})(?:_Z(\d{1,4}))?\.(?:jpe?g|png|tiff?)$", re.I)
_ST_IMG = re.compile(r"^(?P<sn>.+?)_(?P<stain>[^_]+)_regist\.tiff?$", re.I)
_ST_MASK = re.compile(r"^(?P<sn>.+?)_(?P<stain>[^_]+)_tissue_cut\.tiff?$", re.I)
_ST_STAINS = ("he", "h&e", "ssdna", "dapi")                                     # the order the docs list them


def norm_kind(kind):
    """A kind name or one of its other spellings -> the prepare.py kind, or None."""
    if not kind:
        return None
    k = str(kind).strip().lower()
    k = ALIASES.get(k, k)
    return k if k in KINDS else None


# ---------------------------------------------------------------------------------------
# listing: names and sizes only
# ---------------------------------------------------------------------------------------
def _join(rel: str, *names: str) -> str:
    out = rel
    for n in names:
        out = f"{out}/{n}" if out else n
    return out


def _parent(rel: str) -> str:
    return rel.rpartition("/")[0]


def _depth(rel: str) -> int:
    return 0 if not rel else rel.count("/") + 1


def _under(rel: str, root: str) -> bool:
    return root == "" or rel == root or rel.startswith(root + "/")


def _lca(rels) -> str:
    parts = [r.split("/") if r else [] for r in rels]
    out = []
    for xs in zip(*parts):
        if any(x != xs[0] for x in xs):
            break
        out.append(xs[0])
    return "/".join(out)


class Tree:
    """The folders under one root: files[rel] = {name: size}, dirs[rel] = [subfolder names].
    `rel` is the folder's path below the root with forward slashes ("" = the root itself)."""

    def __init__(self, root, inside: str | None = None):
        self.root = Path(root) if root is not None else None
        self.inside = inside                 # a zip's name when this is a zip's file list
        self.files: dict[str, dict[str, int]] = {}
        self.dirs: dict[str, list[str]] = {}
        self.unreadable: list[str] = []
        self.truncated: str | None = None
        self.deep = False                    # folders below the depth limit exist
        self.n_files = 0
        self.seconds = 0.0
        self.single = None                   # the file name when a single file was given

    @property
    def real(self) -> bool:
        return self.inside is None

    def name(self, rel: str) -> str:
        if rel:
            return rel.rpartition("/")[2]
        return self.root.name if self.root is not None else ""

    def path(self, rel: str, name: str | None = None) -> str:
        rel = _join(rel, name) if name else rel
        if not self.real:
            return f"{self.inside}!/{rel}"
        return str(self.root / rel) if rel else str(self.root)

    def f(self, rel: str) -> dict[str, int]:
        return self.files.get(rel, {})

    def d(self, rel: str) -> list[str]:
        return self.dirs.get(rel, [])

    def below(self, root: str, max_depth: int | None = None):
        base = _depth(root)
        for rel in self.files:
            if _under(rel, root) and (max_depth is None or _depth(rel) - base <= max_depth):
                yield rel

    def _ensure(self, rel: str):
        if rel in self.files:
            return
        self.files[rel] = {}
        self.dirs.setdefault(rel, [])
        if rel:
            par = _parent(rel)
            self._ensure(par)
            self.dirs[par].append(rel.rpartition("/")[2])


def _listdir(path) -> tuple[dict, list] | None:
    """({file name: size}, [folder names]) of one folder, or None when it cannot be read."""
    files, dirs = {}, []
    try:
        with os.scandir(path) as it:
            for e in it:
                try:
                    if e.is_dir(follow_symlinks=False):
                        dirs.append(e.name)
                    elif e.is_file():
                        files[e.name] = e.stat().st_size
                except OSError:
                    continue
    except OSError:
        return None
    return files, sorted(dirs, key=str.lower)


def scan(root: Path, max_depth: int = MAX_DEPTH, max_entries: int = MAX_ENTRIES,
         max_seconds: float = MAX_SECONDS) -> Tree:
    """List `root` breadth-first, to a depth, entry and time limit. Zarr stores and hidden
    folders are not looked into."""
    t = Tree(root)
    t0 = time.time()
    level, n = [""], 0
    for depth in range(max_depth + 1):
        nxt = []
        for rel in level:
            got = _listdir(t.root / rel if rel else t.root)
            if got is None:
                t.unreadable.append(rel or ".")
                if not rel:
                    t.seconds = time.time() - t0
                    return t
                continue
            files, dirs = got
            t.files[rel] = files
            t.dirs[rel] = dirs
            n += len(files) + len(dirs)
            t.n_files += len(files)
            if ".zgroup" in files or ".zarray" in files:           # a zarr store: thousands of chunk files
                continue
            for dn in dirs:
                low = dn.lower()
                if low.startswith(".") or low in _SKIP_DIRS or low.endswith(".zarr"):
                    continue
                nxt.append(_join(rel, dn))
            if n > max_entries or time.time() - t0 > max_seconds:
                t.truncated = (f"stopped after {n:,} entries" if n > max_entries
                               else f"stopped after {max_seconds:.0f} s of listing")
                t.seconds = time.time() - t0
                return t
        level = nxt
        if not level:
            break
    else:
        t.deep = bool(level)
    t.seconds = time.time() - t0
    return t


def zip_tree(path: Path, max_names: int = MAX_ENTRIES) -> Tree | None:
    """The file list of a zip as a Tree (nothing is unpacked), or None when it is not a readable zip."""
    try:
        with zipfile.ZipFile(str(path)) as z:
            infos = z.infolist()
    except (zipfile.BadZipFile, OSError, ValueError, NotImplementedError):
        return None
    t = Tree(None, inside=Path(path).name)
    t._ensure("")
    for i in infos[:max_names]:
        name = i.filename.replace("\\", "/").strip("/")
        if not name or name.startswith("__MACOSX/"):
            continue
        if i.is_dir():
            t._ensure(name)
            continue
        d, _, fn = name.rpartition("/")
        t._ensure(d)
        t.files[d][fn] = i.file_size
        t.n_files += 1
    return t


# ---------------------------------------------------------------------------------------
# small headers
# ---------------------------------------------------------------------------------------
def tiff_info(path) -> dict | None:
    """Size and channel layout of a TIFF from its directories (no pixel is read). None when
    tifffile is not installed; {"error": ...} when the file is not a readable TIFF."""
    try:
        import tifffile
    except ImportError:
        return None
    try:
        with tifffile.TiffFile(str(path)) as tf:
            p = tf.pages[0]
            samples = int(p.samplesperpixel)
            info = {"width": int(p.imagewidth), "height": int(p.imagelength), "rgb": samples >= 3,
                    "dtype": str(p.dtype), "ome": bool(tf.is_ome), "qptiff": bool(getattr(tf, "is_qpi", False)),
                    "channels": 1}
            if info["ome"]:                                    # the OME block of the first page lists the channels
                m = re.search(r'SizeC="(\d+)"', p.description or "")
                info["channels"] = int(m.group(1)) if m else 1
                if info["rgb"]:
                    info["channels"] = max(1, info["channels"] // samples)
            else:
                s = tf.series[0]
                n = 1
                for a, k in zip(s.axes, s.shape):
                    if a not in "YXS":
                        n *= int(k)
                info["channels"] = n
            return info
    except Exception as e:
        return {"error": type(e).__name__}


def tiff_variant(path, info: dict | None) -> tuple[str, str]:
    """(variant, plain words) of one image file."""
    low = str(path).lower()
    if info is None or "error" in info:
        return ("qptiff", "QPTIFF scan (not opened)") if low.endswith(".qptiff") else ("tiff", "TIFF image (not opened)")
    size = f"{info['width']:,} x {info['height']:,} px"
    if info["qptiff"] or low.endswith(".qptiff"):
        return "qptiff", f"QPTIFF scan, {info['channels']} channels, {size}"
    if info["rgb"] and info["channels"] == 1:
        return "he", f"RGB image (H&E or other brightfield), {size}"
    if info["channels"] > 1:
        return ("ome-tiff", f"OME-TIFF, {info['channels']} channels, {size}") if info["ome"] \
            else ("tiff-stack", f"TIFF stack, {info['channels']} channels, {size}")
    return "single-channel", f"one grey channel, {size}"


def png_size(path) -> tuple[int, int] | None:
    try:
        with open(path, "rb") as fh:
            b = fh.read(24)
    except OSError:
        return None
    if len(b) < 24 or b[:8] != b"\x89PNG\r\n\x1a\n":
        return None
    return int.from_bytes(b[16:20], "big"), int.from_bytes(b[20:24], "big")


def csv_first_cell(path) -> str | None:
    """The first header cell of a CSV, lower case (quotes removed), or None when it cannot be read."""
    try:
        with open(path, "r", encoding="utf-8-sig", errors="replace") as fh:
            line = fh.readline(4096)
    except OSError:
        return None
    return line.split(",")[0].strip().strip('"').strip().lower()


def h5ad_obsm(path) -> dict | None:
    """{"obsm": [keys]} of an AnnData file, {"error": ...} when it is not one, None without h5py."""
    try:
        import h5py
    except ImportError:
        return None
    try:
        with h5py.File(str(path), "r") as f:
            if "X" not in f and "obs" not in f:
                return {"error": "no X / obs group: not an AnnData file"}
            return {"obsm": sorted(f["obsm"].keys()) if "obsm" in f else []}
    except Exception as e:
        return {"error": f"not a readable HDF5 file ({type(e).__name__})"}


def geomx_sheet_roles(path) -> set | None:
    """The GeoMx column roles (sample_id, roi, segment ...) in the header of an annotation file,
    read with geomx.py's own header rules from the top of the file. None: could not be checked."""
    try:
        import geomx
    except Exception:
        return None
    path = Path(path)
    try:
        if path.suffix.lower() in (".xlsx", ".xlsm"):
            best: set = set()
            for _name, member in geomx.xlsx_sheets(path):
                top = geomx.xlsx_rows(path, member, max_rows=60)
                h = geomx._header_row(top)
                if h is None:
                    continue
                roles = set(geomx._roles([geomx._clean(c) for c in top[h]]))
                if "sample_id" in roles:
                    return roles
                best = best or roles
            return best
        with open(path, "r", encoding="utf-8-sig", errors="replace", newline="") as fh:
            lines = [ln.rstrip("\r\n") for _, ln in zip(range(400), fh)]
        delim = "\t" if sum(ln.count("\t") for ln in lines[:60]) > sum(ln.count(",") for ln in lines[:60]) else ","
        import csv
        rows = list(csv.reader(lines, delimiter=delim))
        h = geomx._header_row(rows)
        return set(geomx._roles([geomx._clean(c) for c in rows[h]])) if h is not None else set()
    except (SystemExit, Exception):
        return set()


# ---------------------------------------------------------------------------------------
# candidates
# ---------------------------------------------------------------------------------------
class Ctx:
    def __init__(self, tree: Tree, overrides: dict | None = None):
        self.t = tree
        self.ov = overrides or {}
        self.peeks = 0

    def tiff(self, path) -> dict | None:
        if not self.t.real or self.peeks >= MAX_PEEKS:
            return None
        self.peeks += 1
        return tiff_info(path)

    def tiff_ok(self, path) -> dict | None:
        """The header when the file could be opened, else None."""
        info = self.tiff(path)
        return info if info and "error" not in info else None


def _cand(kind: str, variant: str, label: str, root: str, is_file: bool = False) -> dict:
    return {"kind": kind, "variant": variant, "label": label, "root": root, "status": "complete",
            "evidence": [], "missing": [], "choices": [], "notes": [], "args": {}, "info": {},
            "_file": is_file, "_claim_dirs": [], "_claim_files": []}


def _ev(c: dict, what: str, where: str, detail: str = ""):
    c["evidence"].append({"what": what, "where": where or ".", "detail": detail})


def _miss(c: dict, what: str, arg: str | None = None, fix: str = ""):
    c["missing"].append({"what": what, "arg": arg, "fix": fix})


def _choice(c: dict, arg: str | None, why: str, options: list[dict]):
    c["choices"].append({"arg": arg, "why": why, "options": options})


def _mb(n: int) -> str:
    return f"{n / 1e9:.2f} GB" if n >= 1e9 else f"{n / 1e6:.1f} MB" if n >= 1e5 else f"{n:,} bytes"


def _first(files: dict, *names: str) -> str | None:
    return next((n for n in names if n in files), None)


def _table(files: dict, stem: str) -> str | None:
    """A table the way flatfiles.find_table picks it: parquet before csv.gz before csv."""
    rx = re.compile(rf"(?:{stem}){_TABLE_EXT}", re.I)
    hits = [n for n in files if rx.search(n)]
    rank = {".parquet": 0, ".gz": 1, ".csv": 2}
    return sorted(hits, key=lambda n: (rank.get(Path(n).suffix.lower(), 3), n.lower()))[0] if hits else None


def _is_tiff(name: str) -> bool:
    return name.lower().endswith(TIFF_EXT)


# ----------------------------------------------------------------------------- xenium
def find_xenium(cx: Ctx) -> list[dict]:
    t, out = cx.t, []
    for rel, f in t.files.items():
        d = t.d(rel)
        stains = sorted(n for n in t.f(_join(rel, "morphology_focus")) if ".ome.tif" in n.lower()) \
            if "morphology_focus" in d else []
        cells = _first(f, "cells.parquet", "cells.csv.gz", "cells.csv")
        mdir = t.f(_join(rel, "cell_feature_matrix")) if "cell_feature_matrix" in d else {}
        matrix = "cell_feature_matrix.h5" if "cell_feature_matrix.h5" in f else \
            ("cell_feature_matrix/" if _first(mdir, "matrix.mtx.gz", "matrix.mtx") else None)
        tx = _first(f, "transcripts.parquet", "transcripts.zarr.zip") or ("transcripts.zarr" if "transcripts.zarr" in d else None)
        old = _first(f, "morphology_focus.ome.tif", "morphology_mip.ome.tif")
        strong = ("experiment.xenium" in f or matrix or tx in ("transcripts.zarr.zip", "transcripts.zarr")
                  or "cells.zarr.zip" in f or old)
        weak = sum(1 for n in ("transcripts.parquet", "cells.parquet", "cells.csv.gz", "cell_boundaries.parquet",
                               "nucleus_boundaries.parquet", "morphology.ome.tif", "gene_panel.json") if n in f)
        if not (stains or strong or weak >= 2):
            continue
        has_cells = bool(cells and matrix)
        parts = (["cells"] if has_cells else []) + (["molecules"] if tx else [])
        variant = "+".join(parts) or "stains"
        c = _cand("xenium", variant, "Xenium (" + ", ".join(["stains"] * bool(stains) + parts) + ")", rel)
        if stains:
            _ev(c, "stain images", _join(rel, "morphology_focus") + "/",
                f"{len(stains)} file(s): " + ", ".join(stains[:4]) + (" ..." if len(stains) > 4 else ""))
            c["_claim_dirs"].append(_join(rel, "morphology_focus"))
        else:
            why = (f" ({old} is the single-file layout of Xenium Onboard Analysis 1.x: the xenium reader needs the "
                   f"morphology_focus/ folder of 2.0 and later; the file can be prepared on its own with --kind tiff)"
                   if old else "")
            _miss(c, "the stain images morphology_focus/*.ome.tif" + why, None,
                  "copy the morphology_focus folder of the Xenium output into this folder")
        if "experiment.xenium" in f:
            _ev(c, "pixel size", _join(rel, "experiment.xenium"))
        else:
            c["notes"].append("no experiment.xenium: the pixel size comes from the image's own OME metadata")
        if has_cells:
            _ev(c, "cells layer", rel, f"{cells} + {matrix}")
        else:
            lack = [n for n, v in (("cells.parquet / cells.csv.gz", cells), ("cell_feature_matrix", matrix)) if not v]
            c["notes"].append("no " + " and no ".join(lack) + ": no cells layer")
        if "cell_feature_matrix" in d:
            c["_claim_dirs"].append(_join(rel, "cell_feature_matrix"))
        if tx:
            _ev(c, "molecules layer", _join(rel, tx))
        else:
            c["notes"].append("no transcripts file: no molecules layer")
        if "Custom_annotation.geojson" in f:
            _ev(c, "overlay", _join(rel, "Custom_annotation.geojson"))
        c["args"]["src"] = t.path(rel)
        out.append(c)
    return out


# ------------------------------------------------------------------------------ visium
def _visium_expect(cx: Ctx, outs: str, squares: list[str]):
    """The full-resolution image size Space Ranger recorded: tissue_hires_image.png / tissue_hires_scalef."""
    t = cx.t
    if not t.real:
        return None
    for sp in [_join(outs, "spatial")] + [_join(outs, "binned_outputs", s, "spatial") for s in squares]:
        f = t.f(sp)
        if "scalefactors_json.json" not in f:
            continue
        png = sp if "tissue_hires_image.png" in f else (_join(outs, "spatial") if "tissue_hires_image.png" in t.f(_join(outs, "spatial")) else None)
        if png is None:
            continue
        try:
            sf = float(json.loads(Path(t.path(sp, "scalefactors_json.json")).read_text(encoding="utf-8"))["tissue_hires_scalef"])
        except (OSError, ValueError, KeyError, TypeError):
            continue
        size = png_size(t.path(png, "tissue_hires_image.png"))
        if size and sf > 0:
            return {"width": size[0] / sf, "height": size[1] / sf, "tol": 1.0 / sf + 2.0}
    return None


def _visium_images(cx: Ctx, outs: str, droot: str, expect, other_outs: set) -> list[list[dict]]:
    """TIFFs that could be the image given to spaceranger, nearest first, as levels of options.

    Free levels (any TIFF there is a candidate): the dataset's own folder, i.e. next to an `outs`
    folder and in the folders next to it, or the top of a Space Ranger folder that is not called
    outs. Strict levels (marked `_outside`: only a TIFF of the size Space Ranger recorded counts,
    so they are skipped when that size is unknown): everything further out, which may belong to
    other datasets - one level further up, and the parent of the folder given."""
    t = cx.t
    levels: list[list[dict]] = []
    named_outs = t.name(outs).lower() == "outs"

    def opts(files: dict, base_path, rel, strict: bool):
        out = []
        for n in sorted(files, key=str.lower):
            if _is_tiff(n):
                o = {"path": str(base_path / n) if base_path is not None else t.path(rel, n), "size": files[n],
                     "_rel": rel, "_name": n}
                if strict:
                    o["_outside"] = True
                out.append(o)
        return out

    def holds_outs(rel: str) -> bool:
        return any(_under(o, rel) for o in other_outs | {outs})

    def listed(base: str, strict: bool, own_only: bool = False):
        levels.append(opts(t.f(base), None, base, strict))
        if not own_only:
            levels.append([o for dn in t.d(base) if not holds_outs(_join(base, dn))
                           for o in opts(t.f(_join(base, dn)), None, _join(base, dn), strict)])

    def unlisted(folder: Path, skip: str, strict_files: bool):
        """The parent of the folder given: one listing, plus one per folder in it (at most 40)."""
        got = _listdir(folder) if t.real else None
        if got is None:
            return
        files, dirs = got
        levels.append(opts(files, folder, None, strict_files))
        sub = []
        for dn in dirs[:40]:
            if dn == skip or dn.lower().startswith("."):
                continue
            g2 = _listdir(folder / dn)
            if g2 and "binned_outputs" not in g2[1] and "outs" not in g2[1]:
                sub += opts(g2[0], folder / dn, None, True)
        levels.append(sub)

    can_check = expect is not None
    if outs and named_outs:                          # <dataset>/outs: the dataset folder, then further out
        listed(_parent(outs), strict=False)
        if can_check and _parent(outs):
            listed(_parent(_parent(outs)), strict=True)
        elif can_check and t.root is not None:
            unlisted(t.root.parent, t.root.name, strict_files=True)
    elif outs:                                       # a Space Ranger folder under another name
        listed(outs, strict=False, own_only=True)
        if can_check:
            listed(_parent(outs), strict=True)
    elif named_outs:                                 # the folder given IS outs: the TIFFs beside it
        unlisted(t.root.parent, t.root.name, strict_files=False)
        if not can_check:
            levels.pop()                             # folders beside it: only with a size to check
    else:                                            # the folder given is a Space Ranger folder under another name
        listed("", strict=False, own_only=True)
        if can_check and t.root is not None:
            unlisted(t.root.parent, t.root.name, strict_files=True)
    return levels


def find_visium(cx: Ctx) -> list[dict]:
    t, found = cx.t, []
    for rel in sorted(t.files, key=_depth):
        if any(_under(rel, o) for o, _, _ in found):
            continue                                              # square_NNNum/, segmented_outputs/ of an outs
        d = t.d(rel)
        squares = sorted(s for s in t.d(_join(rel, "binned_outputs")) if _SQUARE.match(s)) if "binned_outputs" in d else []
        sp = t.f(_join(rel, "spatial")) if "spatial" in d else {}
        classic = any(n in sp for n in _POSITIONS + ("scalefactors_json.json", "tissue_hires_image.png"))
        if squares or classic:
            found.append((rel, squares, sp))
    all_outs = {o for o, _, _ in found}
    out = []
    for rel, squares, sp in found:
        f, d, name = t.f(rel), t.d(rel), t.name(rel)
        par_name = t.name(_parent(rel)) if rel else (t.root.parent.name if t.real else "")
        if _SQUARE.match(name) and par_name == "binned_outputs" and not squares:
            c = _cand("visium", "visium-hd", "one bin-size folder of a Visium HD run", rel)
            _ev(c, "one Visium HD bin folder", rel + "/" if rel else name + "/")
            _miss(c, "this is binned_outputs/" + name + ", one bin size of a Visium HD run: point at the outs folder "
                     "(two levels up) so that every bin size is read", None, "give the outs folder instead")
            out.append(c)
            continue
        layers, bad = [], []
        for s in squares:
            base = _join(rel, "binned_outputs", s)
            lack = ([] if "filtered_feature_bc_matrix.h5" in t.f(base) else ["filtered_feature_bc_matrix.h5"]) + \
                   ([] if "tissue_positions.parquet" in t.f(_join(base, "spatial")) else ["spatial/tissue_positions.parquet"])
            if lack:
                bad.append(f"binned_outputs/{s}: no " + " and no ".join(lack))
            else:
                layers.append(f"{int(_SQUARE.match(s).group(1)):03d}um")
        pos = _first(sp, *_POSITIONS)
        if "filtered_feature_bc_matrix.h5" in f and pos:
            layers.append("spots")
        elif not squares:
            if "filtered_feature_bc_matrix.h5" not in f:
                bad.append("no filtered_feature_bc_matrix.h5"
                           + (" (only the matrix folder filtered_feature_bc_matrix/ is here: the Visium reader reads the .h5 file)"
                              if "filtered_feature_bc_matrix" in d else "")
                           + (" (raw_feature_bc_matrix.h5 is not read)" if "raw_feature_bc_matrix.h5" in f else ""))
            if not pos:
                bad.append("no spatial/tissue_positions.parquet (or tissue_positions.csv / tissue_positions_list.csv)")
        hd = any(k != "spots" for k in layers) or bool(squares)
        variant = "visium-hd" if hd else "visium"
        droot = _parent(rel) if (rel and name.lower() == "outs") else rel
        c = _cand("visium", variant, "Visium HD" if hd else "Visium (classic, 55 µm spots)", droot)
        c["info"] = {"layers": layers, "outs": rel}
        c["_claim_dirs"].append(rel)
        if layers:
            _ev(c, "Space Ranger output", (rel + "/") if rel else ".", "layers " + ", ".join(layers))
            c["notes"] += [b + ": that layer is skipped" for b in bad]
        else:
            _ev(c, "Space Ranger output", (rel + "/") if rel else ".", "no complete count layer")
            for b in bad:
                _miss(c, b)
        if "filtered_feature_bc_matrix" in d and "analysis" in d:
            c["notes"].append("the matrix folder and analysis/ here would also load as a UMAP-only single-cell dataset "
                              "(--kind rnaseq); with spatial/ present this is read as Visium")
        c["args"]["outs"] = t.path(rel)
        # ---- the image given to spaceranger: not in the Space Ranger folder
        expect = _visium_expect(cx, rel, squares)
        levels = _visium_images(cx, rel, droot, expect, all_outs - {rel})
        if "image" in cx.ov:                                       # given by the user: nothing to decide
            c["_claim_files"] += [(o["_rel"], o["_name"]) for lv in levels for o in lv if o["_rel"] is not None]
        else:
            def check(o):
                if expect is None or "match" in o:
                    return
                info = cx.tiff_ok(o["path"])
                o["pixels"] = f"{info['width']:,} x {info['height']:,}" if info else None
                o["match"] = None if info is None else (abs(info["width"] - expect["width"]) <= expect["tol"]
                                                        and abs(info["height"] - expect["height"]) <= expect["tol"])
            seen, chosen = [], None
            for lv in levels:
                if not lv:
                    continue
                for o in lv[:12]:
                    check(o)
                if lv[0].get("_outside"):
                    lv = [o for o in lv if o.get("match")]
                    if not lv:
                        continue
                seen.append(lv)
                if expect is None or any(o.get("match") is not False for o in lv):
                    chosen = lv
                    break
            cyt = _join(rel, "spatial")
            if expect is not None and not (chosen and any(o.get("match") for o in chosen)) and "cytassist_image.tiff" in t.f(cyt):
                o = {"path": t.path(cyt, "cytassist_image.tiff"), "size": t.f(cyt)["cytassist_image.tiff"], "_rel": cyt,
                     "_name": "cytassist_image.tiff"}
                check(o)
                if o.get("match"):                                 # a run without a microscope image
                    chosen = [o]
                    seen.append(chosen)
            options = chosen if chosen is not None else (seen[0] if seen else [])
            for lv in seen:
                c["_claim_files"] += [(o["_rel"], o["_name"]) for o in lv if o["_rel"] is not None]
            exp_txt = f"{expect['width']:,.0f} x {expect['height']:,.0f} px" if expect else ""

            def public(o):
                d_ = {"path": o["path"], "size": o["size"], "sizeText": _mb(o["size"])}
                if o.get("pixels"):
                    d_["pixels"] = o["pixels"]
                if expect is not None:
                    d_["note"] = ("the size Space Ranger recorded" if o.get("match") else
                                  "could not be opened" if o.get("match") is None else
                                  f"not the size Space Ranger recorded ({exp_txt})")
                return d_
            good = [o for o in options if o.get("match")]
            pick = None
            if not options:
                _miss(c, "the full-resolution image given to `spaceranger count --image` (it is not inside the Space "
                         "Ranger folder, and no TIFF was found next to it)", "image",
                      "add --image <that image>; not tissue_hires_image.png, which is downscaled")
            elif expect is not None and len(good) == 1:
                pick = good[0]
                why = f"its size {pick['pixels']} px is the one Space Ranger recorded" + \
                    ("; it is outside the folder given" if pick.get("_outside") else "")
            elif expect is not None and len(good) > 1:
                _choice(c, "image", f"{len(good)} images have the size Space Ranger recorded ({exp_txt}); "
                                    f"which one was given to spaceranger?", [public(o) for o in good])
            elif expect is not None and all(o.get("match") is False for o in options):
                _choice(c, "image", f"no TIFF next to the Space Ranger folder has the size it recorded ({exp_txt}). "
                                    f"A resized copy also needs --image-scale", [public(o) for o in options])
            elif len(options) == 1:
                pick, why = options[0], "the only TIFF next to the Space Ranger folder" + \
                    (" (its size could not be checked)" if expect is not None else
                     " (no tissue_hires_image.png here to check its size against)")
            else:
                _choice(c, "image", f"{len(options)} TIFF files are next to the Space Ranger folder; which one was "
                                    f"given to spaceranger?", [public(o) for o in options])
            if pick is not None:
                c["args"]["image"] = pick["path"]
                _ev(c, "full-resolution image", pick["path"], f"{_mb(pick['size'])}; {why}")
                rest = [o for o in options if o is not pick]
                if rest:
                    c["notes"].append(f"{len(rest)} other TIFF file(s) next to it are not the image: "
                                      + ", ".join(Path(o["path"]).name for o in rest[:4]) + (" ..." if len(rest) > 4 else ""))
        c["args"] = {k: c["args"][k] for k in ("image", "outs") if k in c["args"]}
        out.append(c)
    return out


# ---------------------------------------------------------------------------- merscope
def find_merscope(cx: Ctx) -> list[dict]:
    t = cx.t
    skip = ("images", "result_tiles", "cell_boundaries")
    tx_at, mos_at, cell_at = {}, {}, {}
    for rel, f in t.files.items():
        n = _table(f, "detected_transcripts")
        if n:
            tx_at[rel] = n
        mos = {}
        for fn in f:
            m = _MOSAIC.match(fn)
            if m:
                mos.setdefault(int(m.group(2)), {})[m.group(1)] = fn
        if mos:
            root = _parent(rel) if (rel and t.name(rel) == "images") else rel
            if root not in mos_at or t.name(rel) == "images":
                mos_at[root] = (rel, mos)
        cbg, cmeta = _table(f, "cell_by_gene"), _table(f, "cell_metadata")
        if cbg or cmeta:
            cell_at[rel] = (cbg, cmeta)
    roots = set(tx_at) | set(mos_at)
    tx_of = {r: (r, tx_at[r]) for r in tx_at}
    for r in sorted(tx_at):                           # the reader also looks one folder down for the transcripts
        p = _parent(r)
        if r and p in mos_at and p not in tx_at and t.name(r).lower() not in skip \
                and sum(1 for q in tx_at if q and _parent(q) == p) == 1:
            roots.discard(r)
            tx_of[p] = (r, tx_at[r])
            tx_of.pop(r, None)
    runs: dict[str, list[str]] = {}
    for cd in sorted(cell_at, key=_depth):            # cell tables: the region folder itself or up to two folders down
        anc, owner = cd, None
        for _ in range(3):
            if anc in roots:
                owner = anc
                break
            if not anc:
                break
            anc = _parent(anc)
        if owner is None:
            roots.add(cd)
            owner = cd
        if not any(part.lower() in skip for part in cd[len(owner):].strip("/").split("/") if part):
            runs.setdefault(owner, []).append(cd)
    out = []
    for r in sorted(roots):
        img_rel, mos = mos_at.get(r, (None, {}))
        tdir = img_rel if img_rel is not None else _join(r, "images")
        transform = "micron_to_mosaic_pixel_transform.csv" in t.f(tdir)
        tx = tx_of.get(r)
        mine = runs.get(r, [])
        full = [cd for cd in mine if all(cell_at[cd])]
        cells_dir, cell_note = None, None
        c = _cand("merscope", "", "", r)
        if "cells_dir" in cx.ov:
            cells_dir = "given"
        elif r in mine and all(cell_at[r]):
            cells_dir = r
            others = [cd for cd in full if cd != r]
            if others:
                cell_note = ("other segmentation runs here are not used: " + ", ".join(o[len(r):].strip("/") for o in others)
                             + " (--cells-dir picks one)")
        elif len(full) == 1:
            cells_dir = full[0]
        elif len(full) > 1:
            _choice(c, "cells_dir", f"{len(full)} folders hold cell_by_gene + cell_metadata (segmentation runs); which one?",
                    [{"path": t.path(cd), "size": sum(t.f(cd).values()), "sizeText": _mb(sum(t.f(cd).values()))} for cd in full])
        use_images = bool(mos) and transform and cx.ov.get("images") != "none"
        has_cells = cells_dir is not None or bool(c["choices"])
        parts = (["images"] if use_images else []) + (["cells"] if has_cells else []) + (["molecules"] if tx else [])
        c["variant"] = "+".join(parts) or "nothing"
        c["label"] = "MERSCOPE (" + (", ".join(parts) or "nothing usable") + ")"
        if mos:
            zs = sorted(mos)
            _ev(c, "mosaic images", (img_rel + "/") if img_rel else ".",
                f"{sum(len(v) for v in mos.values())} file(s), z-planes {', '.join(str(z) for z in zs)}, stains "
                + ", ".join(sorted({s for v in mos.values() for s in v})))
            if img_rel != r:
                c["_claim_dirs"].append(img_rel)
            else:
                c["_claim_files"] += [(img_rel, fn) for v in mos.values() for fn in v.values()]
            if transform:
                _ev(c, "pixel positions", _join(tdir, "micron_to_mosaic_pixel_transform.csv"))
            elif cx.ov.get("images") != "none":
                _miss(c, f"micron_to_mosaic_pixel_transform.csv next to the mosaic images ({tdir}/): without it the "
                         f"images cannot be placed", "images",
                      "copy it into the images folder, or add --images none to go without the images")
        else:
            c["notes"].append("no images/mosaic_<stain>_z<k>.tif: a grey transcript-density background is used")
        if tx:
            _ev(c, "molecules layer", _join(tx[0], tx[1]), _mb(t.f(tx[0])[tx[1]]))
        else:
            c["notes"].append("no detected_transcripts file: no molecules layer")
        if cells_dir is not None and cells_dir != "given":
            _ev(c, "cells layer", cells_dir or ".", " + ".join(cell_at[cells_dir]))
        elif not c["choices"] and "cells_dir" not in cx.ov:
            half = [cd for cd in mine if any(cell_at[cd])]
            c["notes"].append(("only one of cell_by_gene / cell_metadata is here" if half else "no cell_by_gene + cell_metadata")
                              + ": no cells layer")
        if cell_note:
            c["notes"].append(cell_note)
        if not (tx or has_cells or use_images) and not c["missing"]:
            _miss(c, "detected_transcripts, or cell_by_gene + cell_metadata, or mosaic images with their transform file")
        c["args"]["src"] = t.path(r)
        out.append(c)
    return out


# ------------------------------------------------------------------------------- cosmx
def find_cosmx(cx: Ctx) -> list[dict]:
    t = cx.t
    tables: dict[str, dict[str, list[str]]] = {}        # folder -> role -> file names
    prefix_of: dict[str, set] = {}
    img_dirs: dict[str, tuple[str, int]] = {}           # folder -> (kind, FOVs)
    for rel, f in t.files.items():
        for fn in f:
            m = _CX_TABLE.match(fn)
            if m:
                tables.setdefault(rel, {}).setdefault(_CX_ROLE[m.group(2).lower()], []).append(fn)
                prefix_of.setdefault(rel, set()).add(m.group(1).lower())
        kind = _CX_IMG.get(t.name(rel).lower()) if rel else None
        if kind:
            fovs = {int(m.group(1)) for m in (_CX_FOV.search(fn) for fn in f) if m}
            if fovs:
                img_dirs[rel] = (kind, len(fovs))
    groups: list[set] = []                              # folders that hold tables of the same slide: one folder,
    for rel in sorted(tables):                          # or folders whose tables share a name and complete each other
        hit = next((g for g in groups if any(prefix_of[rel] & prefix_of[o] - {""} and not set(tables[rel]) & set(tables[o])
                                             for o in g)), None)
        if hit is not None:
            hit.add(rel)
        else:
            groups.append({rel})
    out, used_img = [], set()
    bases = [_lca(sorted(g)) for g in groups]
    for g, base in zip(groups, bases):
        def images_under(root):
            return {r: v for r, v in img_dirs.items() if _under(r, root) and r != root}
        if not images_under(base):                      # image folders one or two levels up belong to the slide
            anc = base
            for _ in range(2):
                if not anc:
                    break
                anc = _parent(anc)
                others = [b for b in bases if b != base and _under(b, anc)]
                if images_under(anc) and not others and all(_depth(r) - _depth(anc) <= 2 for r in g):
                    base = anc
                    break
        roles: dict[str, list[tuple[str, str]]] = {}
        for rel in sorted(g, key=lambda r: (_depth(r), r.lower())):
            for role, names in tables[rel].items():
                stems = {}
                for n in names:
                    stems.setdefault(re.sub(_TABLE_EXT, "", n, flags=re.I).lower(), []).append(n)
                for stem in sorted(stems):
                    roles.setdefault(role, []).append((rel, _table({n: 0 for n in stems[stem]}, ".*")))
        imgs = images_under(base)
        used_img |= set(imgs)
        have = {k: roles.get(k, [None])[0] for k in ("matrix", "cells", "fov_positions", "transcripts")}
        kinds = {}
        for r, (k, nf) in imgs.items():
            kinds[k] = max(kinds.get(k, 0), nf)
        pictures = any(k in kinds for k in ("morphology", "zstack", "composite"))
        has_cells = bool(have["matrix"] and have["cells"])
        placed = pictures and bool(have["fov_positions"] or have["cells"])
        parts = (["images"] if placed else []) + (["cells"] if has_cells else []) + (["molecules"] if have["transcripts"] else [])
        c = _cand("cosmx", "+".join(parts) or "nothing", "CosMx (" + (", ".join(parts) or "nothing usable") + ")", base)
        words = {"matrix": "counts per cell", "cells": "cell table", "fov_positions": "FOV positions", "transcripts": "molecules layer"}
        for k in ("matrix", "cells", "fov_positions", "transcripts"):
            if have[k]:
                _ev(c, words[k], _join(*have[k]), _mb(t.f(have[k][0])[have[k][1]]))
            if len(roles.get(k, [])) > 1:
                _choice(c, None, f"{len(roles[k])} *{[n for n, v in _CX_ROLE.items() if v == k][0]} tables are in this "
                                 f"folder and the CosMx reader would take the first: keep one slide per folder",
                        [{"path": t.path(*x), "size": t.f(x[0])[x[1]], "sizeText": _mb(t.f(x[0])[x[1]])} for x in roles[k]])
        if kinds:
            _ev(c, "FOV images", base or ".", ", ".join(f"{k} ({n} FOVs)" for k, n in sorted(kinds.items())))
        for r in t.below(base, 5):
            if t.name(r).lower() in tuple(_CX_IMG) + _CX_OTHER and r != base:
                c["_claim_dirs"].append(r)
        if not pictures:
            c["notes"].append("no FOV images (Morphology2D / CellComposite ...): a grey transcript-density background is used")
        elif not placed:
            c["notes"].append("FOV images but no *_fov_positions_file or *_metadata_file to place them: they are left out")
        if not has_cells:
            c["notes"].append("no cells layer (needs *_exprMat_file AND *_metadata_file)")
        if not have["transcripts"]:
            c["notes"].append("no *_tx_file: no molecules layer")
        if not (have["transcripts"] or has_cells or placed):
            _miss(c, "*_tx_file, or *_exprMat_file with *_metadata_file, or FOV images with *_fov_positions_file")
        c["args"]["src"] = t.path(base)
        out.append(c)
    loose: dict[str, dict] = {}
    for r, (k, nf) in img_dirs.items():                 # image folders with no table anywhere near
        if r not in used_img:
            loose.setdefault(_parent(r), {})[k] = nf
    for base, kinds in sorted(loose.items()):
        c = _cand("cosmx", "images-only", "CosMx FOV images without tables", base)
        _ev(c, "FOV images", base or ".", ", ".join(f"{k} ({n} FOVs)" for k, n in sorted(kinds.items())))
        _miss(c, "the flat files (*_fov_positions_file or *_metadata_file place the FOV images; *_exprMat_file and "
                 "*_tx_file give cells and molecules)")
        for r in t.below(base, 1):
            if t.name(r).lower() in tuple(_CX_IMG) + _CX_OTHER and r != base:
                c["_claim_dirs"].append(r)
        c["args"]["src"] = t.path(base)
        out.append(c)
    return out


# --------------------------------------------------------------------------- stereoseq
def find_stereoseq(cx: Ctx) -> list[dict]:
    t = cx.t
    items: dict[str, list[tuple]] = {}                  # serial number -> (role, folder, file, stain)
    for rel, f in t.files.items():
        for fn in f:
            low = fn.lower()
            m = _ST_IMG.match(fn)
            if m:
                items.setdefault(m.group("sn"), []).append(("image", rel, fn, m.group("stain")))
                continue
            m = _ST_MASK.match(fn)
            if m:
                items.setdefault(m.group("sn"), []).append(("mask", rel, fn, m.group("stain")))
            elif low.endswith(".gef"):
                items.setdefault(fn.split(".")[0], []).append(("cellbin" if "cellbin" in low else "gef", rel, fn, ""))
            elif low.endswith((".gem", ".gem.gz")):
                items.setdefault(fn.split(".")[0], []).append(("gem", rel, fn, ""))
    items = {sn: v for sn, v in items.items() if any(r != "mask" for r, *_ in v)}
    by_root: dict[str, list[str]] = {}
    for sn, v in items.items():
        by_root.setdefault(_lca(sorted({rel for _, rel, _, _ in v})), []).append(sn)
    for root, sns in by_root.items():                   # renamed files: an image-only and a counts-only name in one folder
        only_img = [s for s in sns if all(r in ("image", "mask") for r, *_ in items[s])]
        only_cnt = [s for s in sns if not any(r == "image" for r, *_ in items[s])]
        if len(sns) == 2 and len(only_img) == 1 and len(only_cnt) == 1:
            items[only_cnt[0]] += items.pop(only_img[0])
            sns.remove(only_img[0])
    out = []
    for root, sns in sorted(by_root.items()):
        for sn in sorted(sns):
            v = items[sn]
            c = _cand("stereoseq", "", "", root)
            c["info"]["sn"] = sn

            def pick(role, rank, arg, what):
                hits = [(rank(fn, st), rel, fn) for r, rel, fn, st in v if r == role]
                if not hits:
                    return None
                best = min(h[0] for h in hits)
                top = sorted((h for h in hits if h[0] == best), key=lambda h: (_depth(h[1]), h[2].lower()))
                if len(top) > 1 and arg not in cx.ov:
                    _choice(c, arg, f"{len(top)} files could be the {what}; which one?",
                            [{"path": t.path(rel, fn), "size": t.f(rel)[fn], "sizeText": _mb(t.f(rel)[fn])} for _, rel, fn in top])
                    return "open"
                rest = [fn for _, _, fn in sorted(hits) if fn != top[0][2]]
                return top[0][1], top[0][2], rest

            def stain_rank(fn, st):
                return _ST_STAINS.index(st.lower()) if st.lower() in _ST_STAINS else len(_ST_STAINS)

            def name_rank(*suffixes):
                want = [f"{sn}{s}".lower() for s in suffixes]
                return lambda fn, st: next((i for i, w in enumerate(want) if fn.lower() == w or fn.lower() == w + ".gz"), len(want))
            image = pick("image", stain_rank, "image", "registered image")
            gef = pick("gef", name_rank(".tissue.gef", ".raw.gef", ".gef"), "gef", "bin GEF")
            gem = pick("gem", name_rank(".tissue.gem", ".gem"), "gem", "GEM file") if gef is None else None
            cellbin = pick("cellbin", name_rank(".cellbin.gef", ".adjusted.cellbin.gef"), "cellbin", "cell bin GEF")
            stain = (_ST_IMG.match(image[1]).group("stain").lower() if image and image != "open" else "")
            mask = pick("mask", lambda fn, st: 0 if st.lower() == stain else 1, "tissue_mask", "tissue mask")
            for r, rel, fn, _ in v:
                if r in ("image", "mask"):
                    c["_claim_files"].append((rel, fn))
            got = {"image": image, "gef": gef, "gem": gem, "cellbin": cellbin, "tissue_mask": mask}
            words = {"image": "registered image", "gef": "bin layers", "gem": "bin layers (GEM text: slower than a GEF)",
                     "cellbin": "cells layer", "tissue_mask": "tissue mask"}
            for k, val in got.items():
                if val and val != "open":
                    c["args"][k] = t.path(val[0], val[1])
                    extra = ""
                    if k == "gef" and ".raw." in val[1].lower():
                        extra = "whole chip: cut to the tissue"
                    _ev(c, words[k], _join(val[0], val[1]), extra)
                    if val[2]:
                        c["notes"].append(f"not used ({words[k].split(' (')[0]}): " + ", ".join(val[2][:4])
                                          + f"; --{k.replace('_', '-')} picks another")
            if gef and any(r == "gem" for r, *_ in v):
                c["notes"].append("the GEM file(s) are not needed: the GEF is read")
            bins = bool(gef or gem)
            cells = bool(cellbin)
            if not image and "image" not in cx.ov:
                _miss(c, f"the registered image {sn}_HE_regist.tif (or {sn}_ssDNA_regist.tif / {sn}_DAPI_regist.tif): "
                         f"one pixel per DNB, not the raw microscope image", "image", "add --image <the registered image>")
            if not (bins or cells) and not any(k in cx.ov for k in ("gef", "gem", "cellbin")):
                _miss(c, f"a count file: {sn}.tissue.gef / {sn}.raw.gef, {sn}.tissue.gem.gz or {sn}.cellbin.gef", "gef",
                      "add --gef, --gem or --cellbin")
            if not bins:
                c["notes"].append("no .gef / .gem: no bin layers")
            if not cells:
                c["notes"].append("no .cellbin.gef: no cells layer")
            parts = (["bins"] if bins else []) + (["cells"] if cells else [])
            c["variant"] = "+".join(parts) or "image-only"
            c["label"] = f"Stereo-seq ({', '.join(parts) or 'image only'})"
            out.append(c)
    return out


# ------------------------------------------------------------------------------- geomx
def find_geomx(cx: Ctx, zips: dict | None = None) -> list[dict]:
    t = cx.t
    dcc_dirs = {rel: sum(1 for n in f if n.lower().endswith(".dcc")) for rel, f in t.files.items()}
    dcc_dirs = {r: n for r, n in dcc_dirs.items() if n}
    pkc_files = [(rel, n) for rel, f in t.files.items() for n in sorted(f) if n.lower().endswith(".pkc")]
    dcc_zips = [(rel, n) for (rel, n), z in (zips or {}).items() if z.get("dcc")]
    pkc_zips = [(rel, n) for (rel, n), z in (zips or {}).items() if z.get("pkc") and not z.get("dcc")]
    if not (dcc_dirs or pkc_files or dcc_zips):
        return []

    def has_pkc(root):
        return any(_under(r, root) for r, _ in pkc_files + pkc_zips)
    roots: dict[str, dict] = {}
    for src in sorted(set(dcc_dirs) | {r for r, _ in dcc_zips}):
        anc = src
        while not has_pkc(anc) and anc:
            anc = _parent(anc)
        root = anc if has_pkc(anc) else src
        roots.setdefault(root, {"dirs": [], "zips": []})
        if src in dcc_dirs:
            roots[root]["dirs"].append(src)
    for rel, n in dcc_zips:
        root = next((r for r in sorted(roots, key=_depth, reverse=True) if _under(rel, r)), None)
        if root is not None:
            roots[root]["zips"].append((rel, n))
    if not roots:                                       # probe kit files only
        for rel, n in pkc_files:
            roots.setdefault(rel, {"dirs": [], "zips": []})
    out = []
    sheet_rx = re.compile(r"\.(xlsx|xlsm|csv|tsv)$|labworksheet.*\.txt$", re.I)
    for root in sorted(roots):
        info = roots[root]

        def scans(r):
            return [(rel, n) for rel in t.below(r, 4) for n in sorted(t.f(rel)) if re.search(r"\.ome\.tiff?$", n, re.I)]
        top = root
        if not scans(top):                               # the scans often sit in a folder next to the counts
            anc = root
            for _ in range(3):
                if not anc:
                    break
                anc = _parent(anc)
                if any(_under(r, anc) and not _under(r, root) for r in roots):
                    break
                if scans(anc):
                    top = anc
                    break
        c = _cand("geomx", "", "", top)
        n_dcc = sum(dcc_dirs[d] for d in info["dirs"])
        if info["dirs"]:
            dcc_rel = _lca(sorted(info["dirs"]))
            c["args"]["dcc"] = t.path(dcc_rel)
            _ev(c, "DCC count files", (dcc_rel + "/") if dcc_rel else ".", f"{n_dcc} file(s)")
            if info["zips"]:
                c["notes"].append("the DCC zip here is not needed: the unzipped .dcc files are read")
        elif len(info["zips"]) == 1:
            c["args"]["dcc"] = t.path(*info["zips"][0])
            _ev(c, "DCC count files", _join(*info["zips"][0]), "a zip of .dcc files (read as it is)")
        elif len(info["zips"]) > 1 and "dcc" not in cx.ov:
            _choice(c, "dcc", f"{len(info['zips'])} zips of .dcc files; which one?",
                    [{"path": t.path(r, n), "size": t.f(r)[n], "sizeText": _mb(t.f(r)[n])} for r, n in info["zips"]])
        elif "dcc" not in cx.ov:
            _miss(c, "the *.dcc count files (a folder or a zip of them)", "dcc", "add --dcc <folder or zip>")
        pk = [(r, n) for r, n in pkc_files if _under(r, top)]
        if pk:
            c["args"]["pkc"] = [t.path(r, n) for r, n in pk]
            _ev(c, "probe kit", _lca(sorted({r for r, _ in pk})) or ".", ", ".join(n for _, n in pk))
        else:
            pz = [(r, n) for r, n in pkc_zips if _under(r, top)]
            if pz:
                c["args"]["pkc"] = [t.path(r, n) for r, n in pz]
                _ev(c, "probe kit", _join(*pz[0]), "a zip of .pkc files (read as it is)")
            elif "pkc" not in cx.ov:
                _miss(c, "the probe kit file(s) *.pkc of the panel(s) on the slide (from the vendor's PKC files page)",
                      "pkc", "add --pkc <file.pkc>")
        # annotation sheets: Sample_ID is what names each AOI's DCC file
        cands = [(rel, n) for rel in t.below(top, 5) for n in sorted(t.f(rel))
                 if sheet_rx.search(n) and not n.startswith("~$") and t.f(rel)[n] < 50e6]
        cands.sort(key=lambda x: (0 if re.search(r"annot|labworksheet|segment", x[1], re.I) else 1, _depth(x[0]), x[1].lower()))
        bases, joins, unchecked = [], [], []
        for rel, n in cands[:16]:
            roles = geomx_sheet_roles(t.path(rel, n)) if t.real else None
            if roles is None:
                if _depth(rel) - _depth(top) == 0 and re.search(r"\.(xlsx|csv|tsv)$|labworksheet", n, re.I):
                    unchecked.append((rel, n))
            elif "sample_id" in roles and roles & {"roi", "segment", "scan", "slide", "aoi"}:
                bases.append((rel, n))
            elif {"roi", "segment"} <= roles:
                joins.append((rel, n))
        if "annotation" in cx.ov:
            pass
        elif bases or unchecked:
            use = (bases or unchecked) + (joins if len(joins) == 1 else [])
            c["args"]["annotation"] = [t.path(r, n) for r, n in use]
            _ev(c, "annotation sheet", _lca(sorted({r for r, _ in use})) or ".",
                ", ".join(n for _, n in use) + ("" if bases else " (header not checked)"))
            if len(joins) > 1:
                c["notes"].append(f"{len(joins)} sheets with ROI + segment columns but no Sample_ID were not added "
                                  f"(--annotation adds the one you want): " + ", ".join(n for _, n in joins[:4]))
        else:
            _miss(c, "an annotation sheet with a Sample_ID column (*.xlsx, *.csv or the readout package's "
                     "*_LabWorksheet.txt)" + (f"; {len(joins)} sheet(s) here have ROI + segment columns but no Sample_ID"
                                              if joins else ""), "annotation", "add --annotation <sheet>")
        sc = scans(top)
        if sc and "image" not in cx.ov:
            c["args"]["image"] = [t.path(r, n) for r, n in sc]
            _ev(c, "scan images", _lca(sorted({r for r, _ in sc})) or ".", ", ".join(n for _, n in sc[:6]))
            c["_claim_files"] += sc
        elif not sc:
            c["notes"].append("no scan image (*.ome.tiff): the ROIs are drawn on a labelled grid, not at their real positions")
        c["variant"] = "scans" if sc else "grid"
        c["label"] = "GeoMx DSP (" + (f"{len(sc)} scan image(s)" if sc else "no scan image: ROIs on a grid") + ")"
        c["args"] = {k: c["args"][k] for k in KIND_ARGS["geomx"] if k in c["args"]}
        out.append(c)
    return out


# ------------------------------------------------------------------- generic: rnaseq
def _is_mtx_dir(f: dict) -> bool:
    return bool(_first(f, "matrix.mtx", "matrix.mtx.gz") and _first(f, "barcodes.tsv", "barcodes.tsv.gz")
                and _first(f, "features.tsv", "features.tsv.gz", "genes.tsv", "genes.tsv.gz"))


def find_rnaseq(cx: Ctx, skip_dir, skip_file) -> list[dict]:
    t, out = cx.t, []
    emb = cx.ov.get("embedding") or "X_umap"
    for rel, f in t.files.items():
        if skip_dir(rel):
            continue
        for n in sorted(f):
            if not n.lower().endswith(".h5ad") or skip_file(rel, n):
                continue
            c = _cand("rnaseq", "h5ad", "single-cell RNA-seq (AnnData .h5ad)", _join(rel, n), is_file=True)
            info = h5ad_obsm(t.path(rel, n)) if t.real else None
            _ev(c, "AnnData file", _join(rel, n), _mb(f[n]))
            if info is None:
                c["notes"].append("the file was not opened: its UMAP (obsm) is not checked")
            elif "error" in info:
                _miss(c, f"{n}: {info['error']}")
            else:
                keys = info["obsm"]
                if any(k in keys for k in (emb, emb.upper(), "X_UMAP", "umap", "UMAP")):
                    _ev(c, "UMAP", f"obsm of {n}", ", ".join(keys))
                else:
                    _miss(c, f"an embedding in obsm ('{emb}' not found; obsm holds: {', '.join(keys) or 'nothing'})", "embedding",
                          "add --embedding <one of those keys>, or save the UMAP into the file from scanpy / Seurat first")
            c["args"]["h5ad"] = t.path(rel, n)
            out.append(c)
    mtx = [rel for rel, f in t.files.items() if _is_mtx_dir(f) and not skip_dir(rel)]
    by_parent: dict[str, list[str]] = {}
    for rel in mtx:
        by_parent.setdefault(_parent(rel) if rel else "", []).append(rel)
    for par, dirs in sorted(by_parent.items()):
        own = dirs == [""]                               # the folder given IS the matrix folder
        c = _cand("rnaseq", "10x-matrix", "single-cell RNA-seq (10x matrix folder + UMAP CSV)", par)
        names = {t.name(r).lower(): r for r in dirs}
        if len(dirs) == 1:
            use = dirs[0]
        elif "filtered_feature_bc_matrix" in names:
            use = names["filtered_feature_bc_matrix"]
            c["notes"].append("the filtered matrix is used; not used: " + ", ".join(t.name(r) for r in dirs if r != use))
        else:
            use = None
            if "mtx" not in cx.ov:
                _choice(c, "mtx", f"{len(dirs)} matrix folders are here; which one?",
                        [{"path": t.path(r), "size": sum(t.f(r).values()), "sizeText": _mb(sum(t.f(r).values()))} for r in dirs])
        if use is not None:
            c["args"]["mtx"] = t.path(use)
            _ev(c, "10x matrix", (use + "/") if use else ".", "matrix.mtx, features.tsv, barcodes.tsv")
        beside = t.f(par) if not own else (_listdir(t.root.parent) or ({}, []))[0] if t.real else {}
        base_path = (lambda *p: t.path(_join(par, *p))) if not own else (lambda *p: str(t.root.parent.joinpath(*p)))
        umap = None
        udir = _join(par, "analysis", "umap")
        if not own:
            proj = [d for d in t.d(udir) if "projection.csv" in t.f(_join(udir, d))]
            pref = [d for d in proj if d == "gene_expression_2_components"] or proj
            if len(pref) == 1:
                umap = base_path("analysis", "umap", pref[0], "projection.csv")
        if umap is None and "umap.csv" in beside:
            umap = base_path("umap.csv")
        if umap:
            c["args"]["umap"] = umap
            _ev(c, "UMAP", umap)
        elif "umap" not in cx.ov:
            _miss(c, "the UMAP: analysis/umap/gene_expression_2_components/projection.csv or a umap.csv "
                     "(barcode,UMAP1,UMAP2) next to the matrix folder", "umap", "add --umap <csv>")
        obs = []
        if not own and "clusters.csv" in t.f(_join(par, "analysis", "clustering", "gene_expression_graphclust")):
            obs.append(base_path("analysis", "clustering", "gene_expression_graphclust", "clusters.csv"))
        for n in sorted(beside):
            if n.lower().endswith(".csv") and n != "umap.csv":
                first = csv_first_cell(base_path(n)) if t.real else None
                if first in ("barcode", "barcodes", ""):
                    obs.append(base_path(n))
                else:
                    c["notes"].append(f"{n} is not used as colour-by columns (its first column is not 'barcode')")
        if obs and "obs" not in cx.ov:
            c["args"]["obs"] = obs
            _ev(c, "colour-by columns", par or ".", ", ".join(Path(o).name for o in obs))
        if _first(beside, "filtered_feature_bc_matrix.h5", "raw_feature_bc_matrix.h5") and not own:
            c["notes"].append("the .h5 matrix files are not read for single-cell data: the matrix folder is")
        c["args"] = {k: c["args"][k] for k in ("mtx", "umap", "obs") if k in c["args"]}
        out.append(c)
    return out


# -------------------------------------------------------------- generic: image files
def find_images(cx: Ctx, skip_dir, skip_file, want: str | None = None) -> list[dict]:
    """Loose image files. `want` ("akoya" / "tiff"): the user said which, so the naming is not asked for."""
    t, out = cx.t, []
    for rel in sorted(t.files):
        if skip_dir(rel):
            continue
        f = t.f(rel)
        tifs = sorted((n for n in f if _is_tiff(n) and not skip_file(rel, n)), key=str.lower)
        if not tifs:
            continue

        def one(n):
            info = cx.tiff(t.path(rel, n))
            variant, words = tiff_variant(n, info)
            c = _cand("tiff", variant, words, _join(rel, n), is_file=True)
            c["info"] = info or {}
            _ev(c, "image file", _join(rel, n), _mb(f[n]))
            if info and "error" in info:
                _miss(c, f"{n} is not a TIFF that can be opened ({info['error']})")
            if variant == "single-channel":
                c["notes"].append("one grey channel: if it is one marker of several, put the marker files together "
                                  "in a folder (one TIFF per marker) and use --kind akoya")
            c["args"]["src"] = t.path(rel, n)
            return c
        if want == "tiff":
            out += [one(n) for n in tifs[:40]]
            continue
        scans = [n for n in tifs if n.lower().endswith(".qptiff")] if want != "akoya" else []
        rest = [n for n in tifs if n not in scans]
        out += [one(n) for n in scans]
        if want == "akoya":
            ak = _cand("akoya", "per-marker", f"one TIFF per marker ({len(rest)} files)", rel)
            ak["args"]["src"] = t.path(rel)
            _ev(ak, "marker files", (rel + "/") if rel else ".", f"{len(rest)} TIFF file(s): " + ", ".join(rest[:4])
                + (" ..." if len(rest) > 4 else ""))
            sizes = {(i["width"], i["height"]) for i in (cx.tiff_ok(t.path(rel, n)) for n in rest[:6]) if i}
            if len(sizes) > 1:
                _miss(ak, "marker files of one size: " + " and ".join(f"{w} x {h}" for w, h in sorted(sizes))
                      + " px were found, and the per-marker reader needs every file to be the same width and height")
            out.append(ak)
            continue
        if len(rest) == 1:
            out.append(one(rest[0]))
        if len(rest) < 2:
            continue
        indexed = sum(1 for n in rest if re.match(r"^\d+[_ ]", n))
        grey = sum(1 for n in rest if re.search(r"_(gray|grey)\.tiff?$", n, re.I))
        ak = _cand("akoya", "per-marker", f"one TIFF per marker ({len(rest)} files)", rel)
        ak["args"]["src"] = t.path(rel)
        _ev(ak, "marker files", (rel + "/") if rel else ".", f"{len(rest)} TIFF file(s): " + ", ".join(rest[:4])
            + (" ..." if len(rest) > 4 else ""))
        probe = [rest[0], rest[-1]] if len(rest) > 2 else rest
        if max(indexed, grey) >= 0.8 * len(rest):        # NN_Marker_..._gray.tif: the Akoya export's own naming
            infos = [cx.tiff_ok(t.path(rel, n)) for n in probe]
            sizes = {(i["width"], i["height"]) for i in infos if i}
            if len(sizes) > 1:
                _miss(ak, "marker files of one size: " + " and ".join(f"{w} x {h}" for w, h in sorted(sizes))
                      + " px were found, and the per-marker reader needs every file to be the same width and height")
            ak["notes"].append("these files rarely store a pixel size: add --pixel-size (µm per pixel) for the scale bar "
                               "and the protein squares")
            out.append(ak)
            continue
        infos = {n: cx.tiff_ok(t.path(rel, n)) for n in rest[:4]}
        known = [i for i in infos.values() if i]
        same = len({(i["width"], i["height"]) for i in known}) == 1
        if known and len(known) == len(infos) and same and all(i["channels"] == 1 and not i["rgb"] for i in known):
            # grey files of one size: the markers of one image, or separate images? Do not guess.
            ak["notes"].append("the names do not follow NN_Marker_..._gray.tif, so this was not assumed")
            out.append(ak)
            alt = _cand("tiff", "single-channel", f"{len(rest)} separate one-channel images", rel)
            _ev(alt, "image files", (rel + "/") if rel else ".", f"{len(rest)} grey TIFF files of {known[0]['width']:,} x "
                                                                 f"{known[0]['height']:,} px")
            _choice(alt, "src", "each file would be its own dataset: give one file as --src",
                    [{"path": t.path(rel, n), "size": f[n], "sizeText": _mb(f[n])} for n in rest[:20]])
            out.append(alt)
        else:
            out += [one(n) for n in rest[:40]]
            if len(rest) > 40:
                out[-1]["notes"].append(f"{len(rest) - 40} more TIFF files in this folder are not listed")
    return out


# ---------------------------------------------------------------------------------------
# putting it together
# ---------------------------------------------------------------------------------------
PLATFORM_FINDERS = {"visium": find_visium, "xenium": find_xenium, "cosmx": find_cosmx, "merscope": find_merscope,
                    "stereoseq": find_stereoseq, "geomx": find_geomx}


def _zip_listing(path: Path) -> dict | None:
    """What a zip holds, from its file list: {"kind", "variant", "label"} of the platform inside,
    and whether it is a zip of .dcc / .pkc files (those the GeoMx reader reads as they are)."""
    zt = zip_tree(path)
    if zt is None:
        return None
    names = [n.lower() for f in zt.files.values() for n in f]
    out = {"dcc": bool(names) and sum(n.endswith(".dcc") for n in names) >= max(1, 0.5 * len(names)),
           "pkc": any(n.endswith(".pkc") for n in names), "files": zt.n_files}
    zc = Ctx(zt)
    best = None
    for kind in ("xenium", "visium", "cosmx", "merscope", "stereoseq", "geomx"):
        for c in PLATFORM_FINDERS[kind](zc):
            if c["kind"] == "geomx" and not any(n.endswith(".dcc") for n in names):
                continue
            score = (len(c["evidence"]), -_depth(c["root"]))
            if best is None or score > best[0]:
                best = (score, c)
    if best:
        out.update(kind=best[1]["kind"], variant=best[1]["variant"], label=best[1]["label"], inner=best[1]["root"])
    return out


def _argv(kind: str, args: dict) -> list[str]:
    out = []
    for k in KIND_ARGS[kind]:
        if k not in args or args[k] in (None, [], ""):
            continue
        opt, v = "--" + k.replace("_", "-"), args[k]
        if k in _REPEAT:
            for x in (v if isinstance(v, (list, tuple)) else [v]):
                out += [opt, str(x)]
        elif k in _NARGS:
            out += [opt, *[str(x) for x in (v if isinstance(v, (list, tuple)) else [v])]]
        else:
            out += [opt, str(v[-1] if isinstance(v, (list, tuple)) else v)]
    return out


def _finish(c: dict, ov: dict) -> dict:
    """Apply what the user gave explicitly, then settle the candidate's status and command."""
    for k, v in ov.items():
        if k in KIND_ARGS[c["kind"]]:
            c["args"][k] = v
            _ev(c, f"--{k.replace('_', '-')} (given by you)", ", ".join(str(x) for x in v) if isinstance(v, (list, tuple)) else str(v))
        if k in KIND_ARGS[c["kind"]] or k in OVERRIDE_FLAGS:
            c["missing"] = [m for m in c["missing"] if m["arg"] != k]
            c["choices"] = [ch for ch in c["choices"] if ch["arg"] != k]
            for x in (v if isinstance(v, (list, tuple)) else [v]) if k in KIND_ARGS[c["kind"]] else []:
                p = str(x).split("=", 1)[-1] if (c["kind"] == "geomx" and k == "image" and not Path(str(x)).exists()) else str(x)
                if not Path(p).exists():
                    _miss(c, f"--{k.replace('_', '-')} {x}: not found", k)
    if c["kind"] == "stereoseq" and "gef" in c["args"] and "gem" in c["args"] and "gem" in ov and "gef" not in ov:
        c["args"].pop("gef")                             # the user's GEM wins over a GEF found here
    c["args"] = {k: c["args"][k] for k in KIND_ARGS[c["kind"]] if k in c["args"]}
    c["argv"] = [c["kind"]] + _argv(c["kind"], c["args"])
    c["status"] = "incomplete" if c["missing"] else "complete"
    return c


def _public(c: dict) -> dict:
    return {k: v for k, v in c.items() if not k.startswith("_")}


def _quote(s: str) -> str:
    s = str(s)
    return f'"{s}"' if (" " in s or "&" in s or "(" in s or "#" in s or not s) else s


def command(c: dict, ds_id: str = "<id>", out: str = "<cache>", dest: str = "<export>", extra=()) -> str:
    """The prepare.py command of a candidate, as one line to paste."""
    return "python prepare.py " + " ".join(_quote(a) for a in prepare_argv(c, ds_id, out, dest, extra))


def prepare_argv(c: dict, ds_id: str, out=None, dest=None, extra=()) -> list[str]:
    """prepare.py's arguments for a candidate: the kind, its files, --id, --out (or --dest for
    single-cell data, which is written straight to the export folder), then `extra` as given."""
    where = ["--dest", str(dest)] if c["kind"] == "rnaseq" else ["--out", str(out)]
    return list(c["argv"]) + ["--id", str(ds_id)] + where + [str(x) for x in extra]


def _pick_options(c: dict) -> list[str]:
    """The exact options that would settle a candidate's open choices."""
    out = []
    for ch in c["choices"]:
        for o in ch["options"]:
            if ch["arg"] is None:
                continue
            if ch["arg"] == "src":
                out.append(f"--src {_quote(o['path'])}")
            else:
                out.append(f"--kind {c['kind']} --{ch['arg'].replace('_', '-')} {_quote(o['path'])}")
    return out


def _unit(cands: list[dict], extras: list[str]) -> dict:
    """The verdict for the candidates that share one dataset folder."""
    plat = [c for c in cands if c["kind"] in PLATFORM_KINDS]
    main = plat or cands
    kinds = sorted({c["kind"] for c in main})
    rep = {"root": main[0]["root"], "path": "", "kind": None, "variant": None, "label": None,
           "candidates": [_public(c) for c in main + [c for c in cands if c not in main]], "extras": extras, "suggestions": []}
    if len(kinds) > 1:
        rep.update(verdict="ambiguous", reason="two-kinds",
                   message="files of more than one kind are in this folder (" + ", ".join(KIND_WORDS[k] for k in kinds)
                           + "): say which one to prepare")
        for c in main:
            picks = _pick_options(c)
            rep["suggestions"] += picks if (picks and c["kind"] != "tiff") else ([f"--kind {c['kind']}"] if c["kind"] != "tiff" else picks)
        return rep
    if len(main) > 1:
        rep.update(verdict="ambiguous", reason="several-in-one-folder", kind=kinds[0],
                   message=f"{len(main)} {KIND_WORDS[kinds[0]]} datasets share this folder: give each its own folder, "
                           f"or name the files of the one you want")
        for c in main:
            rep["suggestions"].append(f"--kind {c['kind']} " + " ".join(_quote(a) for a in c["argv"][1:]))
        return rep
    c = main[0]
    rep.update(kind=c["kind"], variant=c["variant"], label=c["label"])
    if c["missing"]:
        zipped = c.get("zipped")
        rep.update(verdict="incomplete", reason="zipped" if zipped else "missing-files",
                   message=(f"this is a zipped {KIND_WORDS[c['kind']]}: unzip it into a folder first, then point at that folder"
                            if zipped else f"this looks like {c['label']}, but something it needs is missing: "
                            + "; ".join(m["what"] for m in c["missing"])))
        rep["suggestions"] = [m["fix"] for m in c["missing"] if m["fix"]]
        alts = [a for a in cands if a is not c and a["status"] == "complete" and a["kind"] in GENERIC_KINDS]
        if len(alts) > 6:
            rep["extras"] = list(extras) + [f"{len(alts)} loose image / AnnData files are also in this folder; the first 6 are listed"]
            rep["candidates"] = [_public(x) for x in [c] + alts[:6]]
        for alt in alts[:6]:                             # --kind as well: the file's name alone may say "platform"
            src = alt["args"].get("src") or alt["args"].get("h5ad")
            rep["suggestions"].append(f"--kind {alt['kind']}" + (f" --src {_quote(src)}" if alt["_file"] and src else "")
                                      + f"   (prepares {alt['label']} on its own)")
    elif c["choices"]:
        rep.update(verdict="ambiguous", reason="open-choice",
                   message=f"this is {c['label']}, but a choice is open: " + "; ".join(ch["why"] for ch in c["choices"]))
        rep["suggestions"] = _pick_options(c)
    else:
        rep.update(verdict="certain", reason="one-kind", message=f"this is {c['label']}")
    return rep


def detect(path, kind: str | None = None, hint: str | None = None, overrides: dict | None = None,
           max_depth: int = MAX_DEPTH, max_entries: int = MAX_ENTRIES, max_seconds: float = MAX_SECONDS) -> dict:
    """Look at `path` (a folder or one file) and report what dataset(s) it holds.

    kind       only that kind is considered (the "someone selects it" path); the report is
               `incomplete` / `unknown` with the reason when its files are not there
    hint       a kind suggested by where the folder sits (inputs/visium/...): used when the files
               agree with it; a clash is reported as `ambiguous`, never silently resolved
    overrides  options the user already gave ({"image": path, "cells_dir": path, "images": "none",
               "embedding": "X_tsne" ...}): they settle the matching missing file or open choice

    Returns a dict (JSON-serialisable): verdict, reason, message, kind, variant, candidates
    (each with evidence, missing, choices, notes, args, argv), datasets (one entry per dataset
    found: several mean the folder holds several datasets), extras, suggestions, scan."""
    t0 = time.time()
    p = Path(path)
    ov = dict(overrides or {})
    want = None
    if kind:
        want = norm_kind(kind)
        if want is None:
            raise ValueError(f"unknown kind {kind!r}: one of {', '.join(KINDS)}")
    hint = norm_kind(hint)
    rep = {"path": str(p), "verdict": "unknown", "reason": "", "message": "", "kind": None, "variant": None, "label": None,
           "candidates": [], "datasets": [], "extras": [], "suggestions": [], "forcedKind": want, "hint": hint,
           "scan": {"files": 0, "folders": 0, "seconds": 0.0, "truncated": None, "unreadable": []}}

    def done(**kw):
        rep.update(kw)
        rep["scan"]["seconds"] = round(time.time() - t0, 3)
        return rep
    try:
        exists = p.exists()
    except OSError:
        exists = False
    if not exists:
        return done(reason="not-found", message="this path does not exist (check the spelling, and that the drive or share is connected)")
    if p.is_file():
        t = Tree(p.parent)
        t.single = p.name
        try:
            t.files[""] = {p.name: p.stat().st_size}
        except OSError:
            return done(reason="unreadable", message="this file cannot be read (permissions?)")
        t.dirs[""] = []
        t.n_files = 1
    else:
        t = scan(p, max_depth, max_entries, max_seconds)
        if "" not in t.files:
            return done(reason="unreadable", message="this folder cannot be read (permissions, or a share that is not connected)")
        if not t.files[""] and not t.dirs[""]:
            return done(reason="empty", message="this folder is empty")
    rep["scan"].update(files=t.n_files, folders=len(t.files), truncated=t.truncated, unreadable=t.unreadable[:20],
                       depthLimit=max_depth if t.deep else None)
    cx = Ctx(t, ov)

    # ---- zips: read the file list only
    zips: dict[tuple, dict] = {}
    for rel in sorted(t.files, key=_depth):
        for n in sorted(t.f(rel)):
            if n.lower().endswith(".zip") and not n.lower().endswith(".zarr.zip") and len(zips) < MAX_ZIPS:
                z = _zip_listing(Path(t.path(rel, n)))
                if z is not None:
                    zips[(rel, n)] = z

    # ---- platform kinds
    cands: list[dict] = []
    for k, finder in PLATFORM_FINDERS.items():
        if want and want != k:
            continue
        cands += finder(cx, zips) if k == "geomx" else finder(cx)
    for c in cands:
        _finish(c, ov)
    complete_roots = [c["root"] for c in cands if c["status"] == "complete"]      # an open choice still owns its folder
    claim_dirs = [d for c in cands for d in c["_claim_dirs"]]
    claim_files = {x for c in cands if c["status"] == "complete" for x in c["_claim_files"]}
    forced_generic = want in GENERIC_KINDS

    def skip_dir(rel):
        if forced_generic:
            return False
        return any(_under(rel, d) for d in claim_dirs) or any(_under(rel, r) for r in complete_roots)

    def skip_file(rel, n):
        return (not forced_generic) and (rel, n) in claim_files

    # ---- loose image files, AnnData files, 10x matrix folders
    generic: list[dict] = []
    if not want or want in ("tiff", "akoya"):
        generic += [c for c in find_images(cx, skip_dir, skip_file, want) if not want or c["kind"] == want]
    if not want or want == "rnaseq":
        generic += find_rnaseq(cx, skip_dir, skip_file)
    for c in generic:
        _finish(c, ov)

    # ---- what a complete platform dataset leaves unused
    extras_of: dict[str, list[str]] = {}
    if not forced_generic:
        for root in complete_roots:
            n_tif, eg, n_h5 = 0, [], 0
            for rel in t.below(root):
                if any(_under(rel, d) for d in claim_dirs):
                    continue
                inner = [r for r in complete_roots if r != root and _under(rel, r) and _under(r, root)]
                if inner:
                    continue
                for n in t.f(rel):
                    if (rel, n) in claim_files:
                        continue
                    if _is_tiff(n):
                        n_tif += 1
                        eg.append(n) if len(eg) < 4 else None
                    elif n.lower().endswith(".h5ad"):
                        n_h5 += 1
            if n_tif:
                extras_of.setdefault(root, []).append(f"{n_tif} other TIFF file(s) in this folder are not used (" + ", ".join(eg)
                                                      + (" ..." if n_tif > len(eg) else "") + "); one can be prepared on its own with --src <file>")
            if n_h5:
                extras_of.setdefault(root, []).append(f"{n_h5} .h5ad file(s) in this folder are not used; one can be prepared "
                                                      f"on its own with --src <file>")

    # ---- zipped instrument output
    zipped = []
    for (rel, n), z in sorted(zips.items()):
        if z.get("kind") and not z["dcc"]:
            if want and z["kind"] != want:
                continue
            c = _cand(z["kind"], z["variant"], f"{z['label']}, zipped", _join(rel, n), is_file=True)
            c["zipped"] = t.path(rel, n)
            _ev(c, "zip file list", _join(rel, n), f"{_mb(t.f(rel)[n])}, {z['files']:,} files inside; holds {z['label']}")
            _miss(c, f"the files are inside {n}: prepare.py does not read zipped instrument output", None,
                  f"unzip {n} into a folder, then point at that folder" + (f" (the dataset is in {z['inner']}/ inside it)" if z.get("inner") else ""))
            zipped.append(_finish(c, {}))
    other_archives = sorted(n for f in t.files.values() for n in f
                            if n.lower().endswith((".tar.gz", ".tgz", ".tar", ".7z", ".rar")))

    # ---- group into datasets (one per folder / file)
    units: dict[str, list[dict]] = {}
    for c in cands:
        units.setdefault(c["root"], []).append(c)
    incomplete_roots = sorted((r for r, cs in units.items() if all(x["status"] != "complete" for x in cs)), key=_depth, reverse=True)
    for c in generic:
        home = next((r for r in incomplete_roots if _under(c["root"], r)), None) if not forced_generic else None
        units.setdefault(home if home is not None else c["root"], []).append(c)
    real_units = bool(units)
    notes = []
    if zipped and not real_units:
        for c in zipped:
            units.setdefault(c["root"], []).append(c)
    elif zipped:
        notes += [f"{Path(c['zipped']).name} holds a zipped {KIND_WORDS[c['kind']]} (not looked into further)" for c in zipped]
    if other_archives:
        notes.append(f"{len(other_archives)} archive(s) were not looked into (only .zip file lists are read): "
                     + ", ".join(other_archives[:4]) + (" ..." if len(other_archives) > 4 else ""))
    lone = sorted(n for f in t.files.values() for n in f if n.lower().endswith("transcripts.zarr.zip") and n != "transcripts.zarr.zip")
    if lone:
        notes.append(f"{len(lone)} file(s) named like {lone[0]} look like a Xenium transcripts store on its own: the Xenium "
                     f"reader needs the output folder with morphology_focus/ and the file named transcripts.zarr.zip")
    n_zip = sum(1 for f in t.files.values() for n in f if n.lower().endswith(".zip") and not n.lower().endswith(".zarr.zip"))
    if n_zip > len(zips):
        notes.append(f"{n_zip - len(zips)} zip file(s) were not looked into (the first {MAX_ZIPS} are)")
    if t.truncated:
        notes.append(f"the listing was cut short ({t.truncated}): point at one dataset folder for a full answer")
    if t.unreadable:
        notes.append(f"{len(t.unreadable)} folder(s) could not be read: " + ", ".join(t.unreadable[:4]))

    datasets = []
    for root, cs in sorted(units.items(), key=lambda kv: (_depth(kv[0]), kv[0].lower())):
        extras = list(extras_of.get(root, []))
        clash = False
        if hint and not want:                            # the folder's own name as a hint (inputs/visium/<id>/)
            keep = [c for c in cs if c["kind"] == hint]
            if keep and len(keep) != len(cs):
                extras.append(f"the folder name says {hint}: files of other kinds here were left out")
            clash = not keep
            cs = keep or cs
        d = _unit(cs, extras)
        if clash:
            d.update(verdict="ambiguous", reason="kind-mismatch",
                     message=f"the folder name says {KIND_WORDS[hint]}, but the files look like "
                             + (d["label"] or ", ".join(sorted({c["label"] for c in d["candidates"]}))) + ": say which it is",
                     suggestions=[f"--kind {k}" for k in sorted({c["kind"] for c in d["candidates"]})])
        d["path"] = t.path(root) if not cs[0]["_file"] or root != cs[0]["root"] else t.path(cs[0]["root"])
        datasets.append(d)
    rep["datasets"] = datasets
    rep["extras"] = notes
    if not datasets:
        if want:
            msg = f"no {KIND_WORDS[want]} was found here. It needs: {NEEDS[want]}"
            return done(verdict="incomplete", reason="kind-not-found", message=msg, kind=want)
        exts: dict[str, int] = {}
        for f in t.files.values():
            for n in f:
                e = "".join(Path(n).suffixes[-2:]).lower() or "(no extension)"
                exts[e] = exts.get(e, 0) + 1
        top = ", ".join(f"{n} {e}" for e, n in sorted(exts.items(), key=lambda kv: -kv[1])[:6])
        what = "this file is" if t.single else "nothing in this folder is"
        return done(reason="nothing-recognised",
                    message=f"{what} recognised as a dataset SODA can prepare" + (f" (files here: {top})" if top else "")
                            + ". If you know what it is, say so with --kind")
    if len(datasets) == 1:
        d = datasets[0]
        rep.update({k: d[k] for k in ("verdict", "reason", "message", "kind", "variant", "label", "candidates", "suggestions")})
        rep["extras"] = d["extras"] + notes
        if d["root"] and not t.single and d["verdict"] == "certain":
            is_file = d["candidates"][0]["kind"] in ("tiff", "rnaseq") and Path(d["path"]).is_file()
            rep["extras"].insert(0, f"the dataset is the file {d['root']} in the folder given" if is_file
                                 else f"the dataset is in {d['root']}/ below the folder given")
        return done()
    kinds = sorted({d["kind"] for d in datasets if d["kind"]})
    return done(verdict="ambiguous", reason="several-datasets",
                message=f"this folder holds {len(datasets)} datasets" + (f" ({', '.join(kinds)})" if kinds else "")
                        + ": prepare each one on its own",
                suggestions=[f"--src {_quote(d['path'])}" for d in datasets[:40]])


# ---------------------------------------------------------------------------------------
# plain words
# ---------------------------------------------------------------------------------------
def _rel_to(base: str, where: str) -> str:
    try:
        return Path(where).relative_to(base).as_posix() if Path(where).is_absolute() else where
    except ValueError:
        return where


def _format_candidate(c: dict, indent: str, base: str) -> list[str]:
    out = []
    for e in c["evidence"]:
        out.append(f"{indent}found    {e['what']}: {_rel_to(base, e['where'])}" + (f"  ({e['detail']})" if e["detail"] else ""))
    for m in c["missing"]:
        out.append(f"{indent}MISSING  {m['what']}")
    for ch in c["choices"]:
        out.append(f"{indent}CHOOSE   {ch['why']}")
        for o in ch["options"]:
            extra = ", ".join(x for x in (o.get("sizeText"), (o.get("pixels") or "") + (" px" if o.get("pixels") else ""), o.get("note")) if x)
            out.append(f"{indent}           {o['path']}  ({extra})")
    for n in c["notes"]:
        out.append(f"{indent}note     {n}")
    return out


def format_report(rep: dict, ds_id: str = "<id>", out: str = "<cache>", dest: str = "<export>") -> str:
    """The report in plain words (what `prepare.py detect` prints)."""
    L = [str(rep["path"])]

    def one(d: dict, ind: str, src: str, brief: bool = False):
        L.append(f"{ind}{d['verdict'].upper()}: {d['message']}")
        cs = d["candidates"]
        for i, c in enumerate(cs):
            if len(cs) > 1:
                L.append(f"{ind}  candidate {i + 1}: {c['label']}  [{c['kind']}]")
            L.extend(_format_candidate(c, ind + ("    " if len(cs) > 1 else "  "), src))
        for x in d.get("extras", []):
            L.append(f"{ind}  also     {x}")
        if d["verdict"] == "certain" and brief:
            L.append(f"{ind}  prepare it with: python prepare.py auto --src {_quote(src)} --id {ds_id} "
                     + (f"--dest {dest}" if d["kind"] == "rnaseq" else f"--out {out}"))
        elif d["verdict"] == "certain":
            L.append(f"{ind}  prepare it with:")
            L.append(f"{ind}    python prepare.py auto --src {_quote(src)} --id {ds_id} "
                     + (f"--dest {dest}" if d["kind"] == "rnaseq" else f"--out {out}"))
            L.append(f"{ind}  which runs:")
            L.append(f"{ind}    {command(cs[0], ds_id, out, dest)}")
        elif d.get("suggestions"):
            L.append(f"{ind}  to go on, add one of these to `prepare.py auto`:" if d["verdict"] == "ambiguous"
                     else f"{ind}  to go on:")
            for s in d["suggestions"]:
                L.append(f"{ind}    {s}")

    if rep["datasets"] and len(rep["datasets"]) > 1:
        L.append(f"  {rep['verdict'].upper()}: {rep['message']}")
        for i, d in enumerate(rep["datasets"]):
            L.append(f"  [{i + 1}] {d['path']}")
            one(d, "      ", d["path"], brief=True)
        for x in rep["extras"]:
            L.append(f"  also     {x}")
    elif rep["datasets"]:
        d = dict(rep["datasets"][0])
        d["extras"] = rep["extras"]
        one(d, "  ", str(rep["path"]) if rep["verdict"] != "certain" or not d["root"] else d["path"])
    else:
        L.append(f"  {rep['verdict'].upper()}: {rep['message']}")
        for x in rep["extras"]:
            L.append(f"  also     {x}")
    s = rep["scan"]
    L.append(f"  (looked at {s['files']:,} files in {s['folders']:,} folders in {s['seconds']:.2f} s"
             + (f"; {s['truncated']}" if s["truncated"] else "") + ")")
    return "\n".join(L)


def option_lines(rep: dict, ds_id: str, name: str = "OPTIONS") -> list[str]:
    """For the notebook: what to type into its options dictionary to settle a report that is
    not `certain` (the same choices `prepare.py auto` prints as --kind ... options)."""
    def entry(**kv):
        body = ", ".join(f'"{k}": ' + (f'r"{v}"' if k != "kind" else f'"{v}"') for k, v in kv.items())
        return f'{name}["{ds_id}"] = {{{body}}}'
    out = []
    if rep["reason"] == "several-datasets":
        return ["this folder holds several datasets: give each its own folder under inputs/"]
    if rep["reason"] in ("two-kinds", "kind-mismatch"):
        return ["pick one:"] + ["  " + entry(kind=k) for k in dict.fromkeys(c["kind"] for c in rep["candidates"])]
    if rep["reason"] == "several-in-one-folder":
        return ["several datasets of one kind share this folder: move each into its own folder under inputs/"]
    if not rep["candidates"]:
        return ["if you know what it is: " + entry(kind="<" + " | ".join(KINDS) + ">")] if rep["reason"] == "nothing-recognised" else []
    c = rep["candidates"][0]
    for ch in c["choices"]:
        out.append(ch["why"])
        for o in ch["options"]:
            extra = ", ".join(x for x in (o.get("sizeText"), o.get("note")) if x)
            out.append("  " + (entry(kind=c["kind"], **{ch["arg"]: o["path"]}) if ch["arg"] and ch["arg"] != "src"
                               else o["path"]) + (f"   # {extra}" if extra else ""))
    for m in c["missing"]:
        out.append("missing: " + m["what"])
        if c.get("zipped"):
            out.append("  " + m["fix"])
        elif m["arg"] and m["arg"] in KIND_ARGS[c["kind"]]:
            out.append("  or point at it: " + entry(kind=c["kind"], **{m["arg"]: "<path>"}))
        elif m["arg"]:
            out.append("  " + m["fix"].replace("add --" + m["arg"], f'set "{m["arg"]}" in {name}["{ds_id}"]:'))
    return out


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description="Say what kind of dataset a folder holds (read-only).")
    ap.add_argument("path")
    ap.add_argument("--kind", help="only consider this kind: " + ", ".join(KINDS))
    ap.add_argument("--json", action="store_true")
    a = ap.parse_args(argv)
    rep = detect(a.path, kind=a.kind)
    print(json.dumps(rep, indent=1) if a.json else format_report(rep))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
