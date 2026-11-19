r"""Prepare and export a whole folder of datasets in one go, with a log of any that fail.

    python prepare.py batch --src <folder of dataset folders> --out <cache folder> --dest <export folder>
           [--options options.json] [--jobs 1] [--rebuild] [--only id1,id2] [--dry-run]
           [--fmt jpg] [--workers 8] [--bins 008um,016um] [--markers-genes x.json] [--markers-proteins y.json]
           [--title "Our atlas"] [--key <key file>] [--check-ids] [--no-ids] [--logs <folder>] [--threads N]

    (python batch.py ... is the same)

This is the notebook (notebooks/SODA_pipeline.ipynb, Sections 2 to 5) without the notebook:

1. Each folder under --src is one dataset, and its name becomes the dataset id (letters, digits,
   _ and -; anything else becomes _). One image file or .h5ad, a zipped Xenium output or a Xenium
   transcripts store (<name>_transcripts.zarr.zip) on its own counts too; other loose files are ignored,
   and so is any name starting with _ or . (a place to keep what is not a dataset). x.zip beside a folder
   x/ is a second dataset, x_zip, when it can be read as it is (the log says they may be the same data).
   The notebook's older inputs/<kind>/<id>/ layout works as well. The notebook (Section 2) and the
   desktop app find datasets with this same code (find_datasets). Each one is
   looked at with detect.py (`prepare.py detect`). Only a "certain" verdict is run. Anything else
   is skipped and the log says why and what to add to the options file (or the --kind ... flag of
   `prepare.py auto`) to settle it.
2. `prepare.py auto` (the cache), then `prepare.py export` (the export folder), each in its own
   process, so a dataset that fails or runs out of memory does not stop the others.
3. Resumable: a dataset whose <export>/<id>/manifest.json exists is skipped, and one whose
   <cache>/<id>/manifest.json exists is not prepared again (export then fills in what is missing).
   --rebuild does both again (old tiles stay: delete <export>/<id>/ first for a clean export).
4. catalog.py writes <export>/datasets.json, then the sample ID check (deid.py) reads the export
   and the logs. With a key file it is a gate (exit code 3), as in the notebook; without one it
   only warns.
5. The log: <logs>/batch_<date>_<time>.md (a table to read), .json and .csv (for scripts), and per
   dataset <id>_prepare.log and <id>_export.log, all in <logs> (default: a folder "logs" next to
   the export folder, never inside it). Anything that looks like a sample ID is written as
   [sample ID].

The options file is the notebook's OPTIONS dictionary as JSON: one entry per dataset id, with the
same keys (kind, label, src, image, pixel_size, crop, flip, user_layer, no_transcripts ...):

    {"nerve_section1": {"label": "Normal nerve"},
     "two_runs":       {"kind": "merscope", "cells_dir": "D:/data/two_runs/vpt-run-a"},
     "tonsil_codex":   {"cells": "segmentation/tonsil_cp_masks.tif"}}

("cells" / "nuclei" / "cells_units": a lab's own cell segmentation of a protein image, prepare.py --cells.)

Relative paths in it are taken from the folder the options file is in (user_umap, user_obs and
user_h5ad from the dataset's folder, as in the notebook).

Barcodes / cell ids: each export holds every layer's barcodes / cell ids (expr/<key>/ids.txt.gz) by
default, as `prepare.py export` writes them: a layer whose ids look like sample IDs, or do not fit a
known platform's id format, is left out with a *** line in its export log (layerids.py). --no-ids leaves every
list out. You decide what you publish; the sample ID check below reads the lists too.

Exit code: 0 every dataset that ran is fine (skipped ones are listed, not errors); 1 a dataset or
catalog.py failed; 2 nothing to do or a wrong option; 3 the sample ID check is a gate and found
something. Only the standard library is used here; the work is done by prepare.py in child processes.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import subprocess
import sys
import threading
import time
import traceback
from collections import deque
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import deid        # noqa: E402  (only find_ids / check_paths: nothing is written)
import detect      # noqa: E402

PREPARE = HERE / "prepare.py"
CATALOG = HERE / "catalog.py"
TAIL_LINES = 15                 # last lines of a failed step kept in the log
WHOLE_SECTION_GB = 8.5          # largest measured peak (README, Requirements): a whole Stereo-seq chip's export

# The notebook's OPTIONS keys (Section 2), plus "markers" (a marker list for this dataset only).
# Keys that settle which files are read (detect.py overrides, passed to `prepare.py auto`):
FILE_OPTIONS = ("image", "outs", "gef", "gem", "cellbin", "tissue_mask", "dcc", "pkc", "annotation", "cells_dir",
                "h5ad", "mtx", "umap", "obs", "images",
                # standard formats: choices, not files (detect.OVERRIDE_FLAGS)
                "library", "series", "table", "region", "coordinate_system", "image_element",
                "embedding")                    # single-cell .h5ad: which obsm key is the map (default X_umap)
CHOICE_OPTIONS = {"images", "library", "series", "table", "region", "coordinate_system", "image_element", "embedding"}
# Options of the kind's own command, and the kinds that take them:
KIND_OPTIONS = {"pixel_size": None, "channels": ("tiff", "omezarr"), "image_scale": ("visium", "anndata"),
                "no_cells": ("xenium", "cosmx", "merscope"),
                "no_transcripts": ("xenium", "cosmx", "merscope", "spatialdata"),
                "qv_min": ("xenium",), "z": ("cosmx", "merscope", "omezarr", "spatialdata"), "bins": ("stereoseq",),
                "matrix": ("anndata", "spatialdata"), "spatial_key": ("anndata",), "image_res": ("anndata",),
                "layer_kind": ("anndata",), "coords_um": ("anndata", "spatialdata"), "points": ("spatialdata",),
                "no_user_clusters": ("anndata", "spatialdata"),
                "sheet": ("geomx",), "layout": ("geomx",), "keep_slide_names": ("geomx",),
                "crop": ("tiff", "akoya", "visium"), "flip": ("tiff", "akoya", "visium"),
                # a lab's own cell segmentation of a protein image (protcells.py): a label mask TIFF or a GeoJSON
                "cells": ("tiff", "akoya", "omezarr"), "nuclei": ("tiff", "akoya", "omezarr"),
                "cells_units": ("tiff", "akoya", "omezarr")}
SEGMENTATION_OPTIONS = ("cells", "nuclei")      # files, read relative to the options file like the others
USER_OPTIONS = ("user_layer", "user_umap", "user_obs", "user_h5ad", "user_embedding", "user_columns",
                "user_barcode_prefix", "user_barcode_suffix", "user_replace")
OTHER_OPTIONS = ("kind", "label", "group", "src", "markers")
KNOWN_OPTIONS = set(FILE_OPTIONS) | set(KIND_OPTIONS) | set(USER_OPTIONS) | set(OTHER_OPTIONS)
PATH_OPTIONS = set(FILE_OPTIONS) - CHOICE_OPTIONS | {"src", "markers"} | set(SEGMENTATION_OPTIONS)
# the notebook's names for three kinds (also the folder names of its older inputs/<kind>/<id>/ layout)
KIND_FOLDERS = {"visium", "xenium", "stereoseq", "cosmx", "merscope", "geomx", "proteomics_file",
                "proteomics_markers", "scrnaseq"}
PROTEIN_KINDS = ("tiff", "akoya", "omezarr")   # these get the protein marker list, every other kind the gene list

_print_lock = threading.Lock()


def say(*a):
    with _print_lock:
        print(*a, flush=True)


# ---------------------------------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------------------------------
def scrub(text):
    """`text` with anything deid.py would call a sample ID replaced by [sample ID] (logs never hold one)."""
    if not isinstance(text, str) or not text:
        return text
    out, pos = [], 0
    for m in deid._ID.finditer(text):
        if m.start() < pos or m.group(0) in deid._DTYPES or deid._is_gene(text, m):
            continue
        out.append(text[pos:m.start()])
        out.append("[sample ID]")
        pos = m.end()
        while True:                                   # shorthand that continues an ID, as deid.anonymize reads it
            c = deid._CONT.match(text[pos:])
            if not c or not ((len(c.group(2)) == 4 and c.group(1) != "-") or (c.group(1) == "-" and len(c.group(2)) <= 2)):
                break
            out.append(c.group(1) + "[sample ID]")
            pos += c.end()
    out.append(text[pos:])
    return "".join(out)


def scrub_obj(obj):
    if isinstance(obj, dict):
        return {scrub(str(k)): scrub_obj(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [scrub_obj(v) for v in obj]
    return scrub(obj) if isinstance(obj, str) else obj


def safe_id(name) -> str:
    """The notebook's rule: a folder name becomes the dataset id (a store's .zarr / .ome.zarr is dropped)."""
    name = re.sub(r"(?:\.ome)?\.zarr$", "", str(name), flags=re.I)
    return re.sub(r"[^A-Za-z0-9_-]+", "_", name).strip("_") or "dataset"


def secs(s) -> str:
    if s is None:
        return ""
    s = int(round(s))
    if s < 60:
        return f"{s} s"
    if s < 3600:
        return f"{s // 60} min {s % 60:02d} s"
    return f"{s // 3600} h {s % 3600 // 60:02d} min"


def gb(n) -> str:
    return "" if not n else f"{n / 2**30:.2f} GB"


def size_text(files, nbytes) -> str:
    if not files:
        return ""
    return f"{files:,} files, " + (f"{nbytes / 2**30:.2f} GB" if nbytes >= 2**30 else f"{nbytes / 2**20:.1f} MB")


def folder_size(p: Path) -> tuple[int, int]:
    n = b = 0
    for root, _dirs, files in os.walk(p):
        for f in files:
            try:
                b += os.path.getsize(os.path.join(root, f))
                n += 1
            except OSError:
                pass
    return n, b


def memory_gb():
    """(total, free) memory in GB, or None when this system does not say (the notebook's check)."""
    try:
        if sys.platform == "win32":
            import ctypes

            class MemStatus(ctypes.Structure):
                _fields_ = [("dwLength", ctypes.c_ulong), ("dwMemoryLoad", ctypes.c_ulong)] + [
                    (n, ctypes.c_ulonglong) for n in ("ullTotalPhys", "ullAvailPhys", "ullTotalPageFile",
                                                      "ullAvailPageFile", "ullTotalVirtual", "ullAvailVirtual",
                                                      "ullAvailExtendedVirtual")]
            m = MemStatus()
            m.dwLength = ctypes.sizeof(MemStatus)
            if ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(m)):
                return m.ullTotalPhys / 2**30, m.ullAvailPhys / 2**30
        elif Path("/proc/meminfo").exists():
            kb = {}
            for line in open("/proc/meminfo"):
                parts = line.split()
                if len(parts) >= 2 and parts[1].isdigit():
                    kb[parts[0].rstrip(":")] = int(parts[1])
            return kb["MemTotal"] / 2**20, kb.get("MemAvailable", kb.get("MemFree", 0)) / 2**20
        elif sys.platform == "darwin":
            total = int(subprocess.run(["sysctl", "-n", "hw.memsize"], capture_output=True, text=True).stdout)
            vm = subprocess.run(["vm_stat"], capture_output=True, text=True).stdout
            page = int(re.search(r"page size of (\d+)", vm).group(1))
            pages = sum(int(n) for _, n in re.findall(r"Pages (free|inactive|speculative):\s+(\d+)", vm))
            return total / 2**30, pages * page / 2**30
    except Exception:
        pass
    return None


# Peak memory of a step. The step reports its own: on Windows a venv's python.exe is a small launcher
# that starts the real interpreter as a second process, so the parent cannot see the right one.
CHILD = "--child"
PEAK_MARK = "[batch] peak memory of this step:"


def own_peak_bytes() -> int | None:
    """This process's peak memory: the peak working set on Windows, the maximum resident set elsewhere."""
    try:
        if sys.platform == "win32":
            import ctypes
            from ctypes import wintypes

            class PMC(ctypes.Structure):
                _fields_ = [("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD)] + [
                    (n, ctypes.c_size_t) for n in ("PeakWorkingSetSize", "WorkingSetSize", "QuotaPeakPagedPoolUsage",
                                                   "QuotaPagedPoolUsage", "QuotaPeakNonPagedPoolUsage",
                                                   "QuotaNonPagedPoolUsage", "PagefileUsage", "PeakPagefileUsage")]
            k32 = ctypes.WinDLL("kernel32")
            k32.GetCurrentProcess.restype = wintypes.HANDLE
            k32.K32GetProcessMemoryInfo.argtypes = (wintypes.HANDLE, ctypes.POINTER(PMC), wintypes.DWORD)
            c = PMC()
            c.cb = ctypes.sizeof(PMC)
            return int(c.PeakWorkingSetSize) if k32.K32GetProcessMemoryInfo(k32.GetCurrentProcess(), ctypes.byref(c), c.cb) else None
        import resource
        return int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) * (1 if sys.platform == "darwin" else 1024)
    except Exception:
        return None


def run_child(argv: list[str]):
    """python batch.py --child <script> <args>: run the script as the main program (prepare.py) and print
    its peak memory when it ends, whether it finished or failed."""
    import atexit
    import runpy
    atexit.register(lambda: print(f"{PEAK_MARK} {own_peak_bytes() or 0}", flush=True))
    sys.argv = list(argv)
    runpy.run_path(argv[0], run_name="__main__")


# ---------------------------------------------------------------------------------------------
# finding the datasets: the one copy, used here, by the notebook (Section 2) and by the desktop app
# ---------------------------------------------------------------------------------------------
def load_options(path: Path | None) -> tuple[dict, list[str]]:
    """The options file: {dataset id: {key: value}}. Returns (options, warnings)."""
    if not path:
        return {}, []
    try:
        raw = json.loads(Path(path).read_text(encoding="utf-8-sig"))
    except (OSError, ValueError) as e:
        raise SystemExit(f"[batch] the options file {path} could not be read: {e}")
    if not isinstance(raw, dict):
        raise SystemExit(f"[batch] the options file {path} must hold one {{...}} with an entry per dataset id")
    base = Path(path).resolve().parent
    out, warn = {}, []
    for ds_id, o in raw.items():
        if str(ds_id).startswith("_"):                 # "_comment": ... and the like
            continue
        if not isinstance(o, dict):
            raise SystemExit(f"[batch] options file: the entry for {ds_id!r} must be a {{...}} of settings")
        o = dict(o)
        unknown = sorted(k for k in o if k not in KNOWN_OPTIONS)
        if unknown:
            warn.append(f"options for {ds_id}: {', '.join(unknown)} not known, left out "
                        f"(known: {', '.join(sorted(KNOWN_OPTIONS))})")
            for k in unknown:
                o.pop(k)
        for k in PATH_OPTIONS & set(o):
            def absolute(v):
                s = str(v)
                if k == "image" and "=" in s and not Path(s).exists():     # GeoMx: "scan name=path"
                    name, p = s.split("=", 1)
                    return f"{name}={absolute(p)}"
                q = Path(s).expanduser()
                return str(q if q.is_absolute() else (base / q).resolve())
            o[k] = [absolute(v) for v in o[k]] if isinstance(o[k], list) else absolute(o[k])
        out[str(ds_id)] = o
    return out, warn


def _json_entry(ds_id: str, **kv) -> str:
    return json.dumps(ds_id) + ": " + json.dumps(kv)


def hint_lines(rep: dict, ds_id: str) -> list[str]:
    """What to add to settle a report that is not certain: the options-file entry (the notebook's
    option_lines, as JSON) and the flags of `prepare.py auto` (detect.py's suggestions)."""
    out = []
    c = rep["candidates"][0] if rep["candidates"] else None
    reason = rep["reason"]
    if reason in ("several-datasets", "several-in-one-folder"):
        out.append("several datasets share this folder: move each into its own folder under --src")
    elif reason in ("two-kinds", "kind-mismatch"):
        out.append("say which kind it is, in the options file (one of):")
        out += ["  " + _json_entry(ds_id, kind=k) for k in dict.fromkeys(x["kind"] for x in rep["candidates"])]
    elif not c:
        if reason == "nothing-recognised" or rep["verdict"] == "unknown":
            out.append("if you know what it is, say so in the options file: "
                       + _json_entry(ds_id, kind="<" + " | ".join(detect.KINDS) + ">"))
    else:
        for ch in c["choices"]:
            out.append(ch["why"] + " Add one of these to the options file:")
            for o in ch["options"]:
                if ch["arg"] and ch["arg"] != "src":
                    out.append("  " + _json_entry(ds_id, kind=c["kind"], **{ch["arg"]: o["path"]}))
                else:
                    out.append("  " + _json_entry(ds_id, src=o["path"]))
        for m in c["missing"]:
            out.append("missing: " + m["what"])
            if c.get("zipped"):
                out.append("  " + m["fix"])
            elif m["arg"] and m["arg"] in detect.KIND_ARGS[c["kind"]]:
                out.append("  or point at it in the options file: " + _json_entry(ds_id, kind=c["kind"], **{m["arg"]: "<path>"}))
            elif m["fix"]:
                out.append("  " + m["fix"])
    flags = [s for s in rep.get("suggestions") or [] if s.startswith("--")]
    if flags:
        out.append("(the same as adding to `prepare.py auto`: " + "  or  ".join(flags[:6]) + ")")
    return out


def _user_files(d: dict, o: dict) -> None:
    """Visium: the lab's own UMAP / clusters (the notebook's user_files): OPTIONS user_* keys, else
    my_umap.csv, my_clusters*.csv, my_analysis.h5ad in the dataset folder."""
    here = d["path"] if d["path"].is_dir() else d["path"].parent

    def path(p):
        p = Path(p)
        return p if p.is_absolute() else here / p
    umap = path(o["user_umap"]) if "user_umap" in o else here / "my_umap.csv"
    obs = o.get("user_obs")
    obs = ([path(p) for p in ([obs] if isinstance(obs, (str, Path)) else obs)] if obs is not None
           else [here / "my_clusters.csv"] + sorted(here.glob("my_clusters_*.csv")))
    h5ad = path(o["user_h5ad"]) if "user_h5ad" in o else here / "my_analysis.h5ad"
    u = {"umap": umap if umap.exists() else None, "obs": [p for p in obs if p.exists()],
         "h5ad": h5ad if h5ad.exists() else None}
    for k, p in (("user_umap", umap), ("user_h5ad", h5ad)):
        if k in o and not p.exists():
            d["problems"].append(f"{k}: not found: {p}")
    if "user_obs" in o:
        d["problems"] += [f"user_obs: not found: {p}" for p in obs if not p.exists()]
    if not (u["umap"] or u["obs"] or u["h5ad"]):
        return
    if u["h5ad"] and (u["umap"] or u["obs"]):
        d["problems"].append("your own results: give either the .h5ad or the CSVs, not both")
    if not o.get("user_layer"):
        d["problems"].append('your own UMAP / clusters found: set "user_layer" (e.g. "008um") in the options file')
        d["hint"].append(_json_entry(d["id"], user_layer="008um"))
    d["user"] = u
    d["notes"].append(f"your own results for layer {o.get('user_layer')}: "
                      + ", ".join(p.name for p in [u["umap"], u["h5ad"], *u["obs"]] if p))


def check_dataset(ds_id: str, src: Path, options: dict, bins: str | None, hint: str | None = None,
                  expand: bool = True) -> list[dict]:
    """One folder (or file): detect.py's report, then the checks of the notebook's Section 2."""
    o = options.get(ds_id, {})
    src = Path(o["src"]) if "src" in o else Path(src)
    ov = {k: o[k] for k in FILE_OPTIONS if k in o}
    try:
        rep = detect.detect(src, kind=o.get("kind"), hint=hint, overrides=ov)
    except ValueError as e:                              # a "kind" that is not one of the kinds
        rep = detect.detect(src, hint=hint)
        rep.update(verdict="ambiguous", reason="kind-mismatch", message=str(e))
    if rep["reason"] == "several-datasets" and expand:   # a folder of datasets: one row each, never merged
        subs = [Path(s["path"]) for s in rep["datasets"]]
        return [r for s in subs for r in check_dataset(safe_id(s.name if s.is_dir() else s.stem), s, options, bins,
                                                        hint, expand=False)]
    c = rep["candidates"][0] if rep["candidates"] else None
    kinds = list(dict.fromkeys(x["kind"] for x in rep["candidates"]))
    d = {"id": ds_id, "src": src, "path": Path(rep["datasets"][0]["path"]) if rep.get("datasets") else src,
         "verdict": rep["verdict"], "kind": rep["kind"] or ("/".join(kinds) + "?" if kinds else "?"),
         "variant": rep["variant"] or "", "label": o.get("label", ds_id), "opts": o,
         "forceKind": bool(o.get("kind") or hint), "problems": [], "hint": [], "notes": [],
         "user": None, "status": None, "reason": "", "steps": {},
         "report": rep}                                  # detect.py's report (the notebook shows it; not in the logs)
    if rep["verdict"] != "certain":
        d["problems"].append(f"not {'sure' if rep['verdict'] == 'ambiguous' else rep['verdict']}: {rep['message']}")
        d["hint"] = hint_lines(rep, ds_id)
        return [d]
    d["notes"] += list(c["notes"]) + list(rep["extras"])
    for k in SEGMENTATION_OPTIONS:                      # a lab's own segmentation: protein images only, and it must exist
        if o.get(k) and d["kind"] not in KIND_OPTIONS[k]:
            d["problems"].append(f"{k}: a lab's own cell segmentation is read for protein images only "
                                 f"({', '.join(KIND_OPTIONS[k])}), not for {d['kind']}")
        elif o.get(k) and not Path(o[k]).is_file():
            d["problems"].append(f"{k}: not found: {o[k]}")
    if o.get("nuclei") and not o.get("cells"):
        d["problems"].append('nuclei: needs "cells" too (the nucleus outlines belong to its cells)')
    if o.get("cells"):
        d["notes"].append(f"your own cell segmentation: {Path(o['cells']).name}"
                          + (f", nuclei {Path(o['nuclei']).name}" if o.get("nuclei") else ""))
    if d["kind"] == "akoya" and "pixel_size" not in o:
        d["notes"].append("no pixel_size option: no scale bar or protein squares unless the files carry it")
    if d["kind"] == "visium":
        layers = c["info"].get("layers", [])
        lack = [b for b in (bins or "").split(",") if b and b != "all" and b not in layers] \
            if d["variant"] == "visium-hd" else []
        if lack:
            d["problems"].append(f"--bins asks for {', '.join(lack)}, but outs/ has only {', '.join(layers)}")
        _user_files(d, o)
    return [d]


ZIP_SUFFIX = "_zip"      # the id of x.zip when the folder x/ beside it is a dataset too


def find_datasets(src: Path, options: dict, bins: str | None) -> tuple[list[dict], list[str]]:
    """Every dataset under `src`: batch.py's, the notebook's (Section 2 calls this) and the desktop app's one
    finder. Returns (datasets, warnings).

    Names starting with "." or "_" are left out (a place to keep what is not a dataset), and so is README.txt.
    x.zip beside a folder x/: both are datasets when the zip can be prepared as it is (a zipped Xenium output), the
    zip with the id x_zip, and a warning says they may be the same data; a zip that cannot (it must be unzipped
    first) is left out with a warning, as the folder beside it is taken to be its unzipped copy."""
    rows, warn = [], []

    def is_kind_folder(p: Path) -> bool:              # inputs/visium/<id>/ ... of the older layout
        return p.name in KIND_FOLDERS and not any(s["root"] == "" for s in detect.detect(p)["datasets"])

    entries = sorted(src.iterdir(), key=lambda q: q.name.lower())
    folders = {q.name.lower() for q in entries if q.is_dir()}
    parked = sum(1 for p in entries if p.name.startswith("_"))
    if parked:
        warn.append(f"{parked} name(s) under the inputs folder start with _ and are left out (not datasets)")
    for p in entries:
        if p.name.startswith((".", "_")) or p.name == "README.txt":
            continue
        if (p.is_file() and p.suffix.lower() == ".zip" and p.stem.lower() in folders
                and safe_id(p.stem) in {r["id"] for r in rows}):          # x/ sorts before x.zip
            zid = safe_id(p.stem) + ZIP_SUFFIX
            got = check_dataset(zid, p, options, bins)
            if got and all(r["verdict"] == "certain" for r in got):
                rows += got
                warn.append(f"{p.name} and the folder {p.stem} beside it are prepared as two datasets, "
                            f"{safe_id(p.stem)} and {zid}: they may be the same data (the folder an unzipped copy of "
                            f"the zip). If they are, keep one: move the other out of the inputs folder")
            else:
                warn.append(f"{p.name} is not prepared: it cannot be read as it is, and the folder {p.stem} beside it "
                            f"is prepared (an unzipped copy, it seems)")
            continue
        if p.is_dir() and is_kind_folder(p):
            for sub in sorted(q for q in p.iterdir() if q.is_dir()):
                rows += check_dataset(safe_id(sub.name), sub, options, bins, hint=p.name)
        elif p.is_dir():
            rows += check_dataset(safe_id(p.name), p, options, bins)
        else:                                    # a file on its own; stray files are not datasets
            rows += [r for r in check_dataset(safe_id(p.stem), p, options, bins) if r["verdict"] != "unknown"]
    have = {r["id"] for r in rows}
    for ds_id, o in options.items():             # datasets that only exist as a `src` somewhere else
        if "src" in o and ds_id not in have:
            rows += check_dataset(ds_id, Path(o["src"]), options, bins, expand=False)
        elif ds_id not in have:
            warn.append(f"the options file names {ds_id}, but there is no such dataset under --src (and no \"src\")")
    seen = set()
    for r in rows:
        if r["id"] != safe_id(r["id"]):
            r["problems"].append("dataset id may only use letters, digits, _ and -")
        if deid.find_ids(r["id"]):
            r["problems"].append("the folder name looks like a sample ID, and it would become the dataset id (in every "
                                 "link and file name): rename the folder")
        if "label" in r["opts"] and deid.find_ids(str(r["label"])):
            r["problems"].append("its label (options file) looks like a sample ID: change it")
        if r["id"] in seen:
            r["problems"].append(f"two datasets would get the id {r['id']}: rename one folder")
        seen.add(r["id"])
    return rows, warn


# ---------------------------------------------------------------------------------------------
# the commands: the only place they are built (the notebook and the desktop app call these too)
# ---------------------------------------------------------------------------------------------
def _flag(key: str) -> str:
    return "--" + key.replace("_", "-")


def gene_file_argv(a) -> list[str]:
    """--gene-files / --shard-mb as given to the batch (genepack.py): for every export, and for
    single-cell RNA-seq, which is written straight to the export by its prepare step."""
    out = ["--gene-files", a.gene_files] if getattr(a, "gene_files", None) else []
    return out + (["--shard-mb", str(a.shard_mb)] if getattr(a, "shard_mb", None) else [])


def prepare_argv(d: dict, a) -> list[str]:
    o, k = d["opts"], d["kind"]
    argv = ["auto", "--src", str(d["src"]), "--id", d["id"], "--out", str(a.out), "--dest", str(a.dest),
            "--label", str(o.get("label", d["id"]))]
    if d["forceKind"]:
        argv += ["--kind", k]
    for key in FILE_OPTIONS:                      # files the options settle (detect.py overrides)
        if key not in o or o[key] in (None, "", []):
            continue
        vals = o[key] if isinstance(o[key], list) else [o[key]]
        if key == "pkc":
            argv += [_flag(key), *[str(v) for v in vals]]
        elif key in ("annotation", "image", "obs"):
            for v in vals:
                argv += [_flag(key), str(v)]
        else:
            argv += [_flag(key), str(vals[-1])]
    for key, kinds in KIND_OPTIONS.items():
        v = o.get(key)
        if (kinds is None or k in kinds) and v is not None and v is not False and v != "":
            if isinstance(v, list):
                v = ",".join(str(x) for x in v)
            argv += [_flag(key)] + ([] if v is True else [str(v)])
    if k == "visium" and d["variant"] == "visium-hd" and a.bins:
        argv += ["--bins", a.bins]
    if k == "rnaseq":
        argv += gene_file_argv(a)
    return argv


def export_argv(d: dict, a) -> list[str]:
    """The export command. Per-row barcodes / cell ids are written by `prepare.py export` by default (a
    layer that fails layerids.py's checks is left out with a message, the export goes on); --no-ids passes
    on to leave them all out. It never passes --ids (which would stop a dataset at a failed check); the
    options file cannot add either (only KNOWN_OPTIONS pass)."""
    o = d["opts"]
    argv = ["export", "--id", d["id"], "--out", str(a.out), "--dest", str(a.dest), "--fmt", a.fmt,
            "--workers", str(a.workers)] + gene_file_argv(a) + (["--no-ids"] if getattr(a, "no_ids", False) else [])
    markers = o.get("markers") or (a.markers_proteins if d["kind"] in PROTEIN_KINDS else a.markers_genes)
    if markers:
        argv += ["--markers", str(markers)]
    u = d["user"]
    if u:                                         # your own UMAP / clusters
        argv += ["--user-layer", str(o["user_layer"])]
        argv += ["--user-umap", str(u["umap"])] if u["umap"] else []
        for p in u["obs"]:
            argv += ["--user-obs", str(p)]
        argv += ["--user-h5ad", str(u["h5ad"])] if u["h5ad"] else []
        for key in ("user_embedding", "user_columns", "user_barcode_prefix", "user_barcode_suffix"):
            if o.get(key):
                argv += [_flag(key), str(o[key])]
        argv += ["--user-replace"] if o.get("user_replace") else []
    return argv


def run_step(argv: list[str], log_path: Path, threads: int | None, running: set) -> dict:
    """One prepare.py command in its own process; its output (sample IDs scrubbed) goes to log_path."""
    args = [*argv] + (["--threads", str(threads)] if threads else [])
    cmd = [sys.executable, str(Path(__file__).resolve()), CHILD, str(PREPARE), *args]
    t0 = time.time()
    tail, warned, peak = deque(maxlen=TAIL_LINES), False, None
    with open(log_path, "a", encoding="utf-8") as lg:
        lg.write(scrub(f"\n=== {datetime.now().isoformat(timespec='seconds')} batch: python prepare.py "
                       + " ".join(detect._quote(x) for x in args) + "\n"))
        lg.flush()
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, encoding="utf-8",
                             errors="replace", env={**os.environ, "PYTHONUNBUFFERED": "1", "PYTHONIOENCODING": "utf-8"})
        running.add(p)
        try:
            for line in p.stdout:
                if line.startswith(PEAK_MARK):
                    try:
                        peak = int(line[len(PEAK_MARK):].strip()) or None
                    except ValueError:
                        pass
                    continue
                line = scrub(line)
                lg.write(line)
                tail.append(line.rstrip("\n"))
                warned = warned or "***" in line
            rc = p.wait()
        finally:
            running.discard(p)
        lg.write(f"=== exit {rc} after {time.time() - t0:.0f} s" + (f", peak memory {gb(peak)}" if peak else "") + "\n")
    return {"exit": rc, "seconds": round(time.time() - t0, 1), "peakBytes": peak, "log": str(log_path),
            "tail": list(tail) if rc != 0 else [], "warning": warned}


def process(d: dict, a, n: int, total: int, logs: Path, threads: int | None, running: set, stop: threading.Event):
    """Prepare then export one dataset; fills in d["status"], d["reason"], d["steps"]."""
    tag = f"[batch] [{n}/{total}] {d['id']}:"
    t0 = time.time()
    out_dir, cache_dir = Path(a.dest) / d["id"], Path(a.out) / d["id"]
    try:
        if (out_dir / "manifest.json").exists() and not a.rebuild:
            d["status"], d["reason"] = "skipped", "already exported (its manifest.json is there); --rebuild to do it again"
            say(f"{tag} skipped, already exported")
            return
        rnaseq = d["kind"] == "rnaseq"
        if (cache_dir / "manifest.json").exists() and not a.rebuild and not rnaseq:
            d["steps"]["prepare"] = {"skipped": "already in the cache"}
            say(f"{tag} prepare: already in the cache")
        else:
            say(f"{tag} prepare ({d['kind']}{', ' + d['variant'] if d['variant'] else ''}) ...")
            st = d["steps"]["prepare"] = run_step(prepare_argv(d, a), logs / f"{d['id']}_prepare.log", threads, running)
            if st["exit"] != 0:
                d["status"], d["reason"] = "failed", f"prepare failed (exit {st['exit']})"
                return
        if stop.is_set():
            d["status"], d["reason"] = "failed", "stopped (Ctrl+C) before the export"
            return
        if not rnaseq:                            # single-cell RNA-seq is written straight to the export
            say(f"{tag} export ...")
            st = d["steps"]["export"] = run_step(export_argv(d, a), logs / f"{d['id']}_export.log", threads, running)
            if st["exit"] != 0:
                d["status"], d["reason"] = "failed", f"export failed (exit {st['exit']})"
                return
            if st["warning"]:
                d["notes"].append(f"the export log has a *** line: read {st['log']} before uploading")
        if not (out_dir / "manifest.json").exists():
            d["status"], d["reason"] = "failed", f"finished, but {out_dir / 'manifest.json'} is not there"
            return
        d["status"] = "ok"
    except Exception as e:                        # never let one dataset stop the others
        d["status"], d["reason"] = "failed", f"batch error: {e}"
        d["steps"]["error"] = {"tail": traceback.format_exc().splitlines()[-TAIL_LINES:]}
    finally:
        d["seconds"] = round(time.time() - t0, 1)
        peaks = [s.get("peakBytes") for s in d["steps"].values() if s.get("peakBytes")]
        d["peakBytes"] = max(peaks) if peaks else None
        if d["status"] == "ok":
            d["exportFiles"], d["exportBytes"] = folder_size(out_dir)
        if d["status"] in ("ok", "failed"):
            say(f"{tag} {d['status'].upper() if d['status'] == 'failed' else d['status']} in {secs(d['seconds'])}"
                + (f" ({d['reason']})" if d["reason"] else "")
                + (f", peak memory {gb(d['peakBytes'])}" if d.get("peakBytes") else "")
                + (f", {size_text(d.get('exportFiles'), d.get('exportBytes', 0))}" if d.get("exportFiles") else ""))


# ---------------------------------------------------------------------------------------------
# the log
# ---------------------------------------------------------------------------------------------
def error_tail(d: dict) -> list[str]:
    for k in ("error", "export", "prepare"):
        s = d["steps"].get(k) or {}
        if s.get("tail"):
            return s["tail"]
    return []


def failed_step(d: dict) -> dict:
    for k in ("export", "prepare"):
        s = d["steps"].get(k) or {}
        if s.get("exit"):
            return s
    return {}


def kind_text(d: dict) -> str:
    return d["kind"] + (f" ({d['variant']})" if d.get("variant") else "")


def row_json(d: dict) -> dict:
    return {"id": d["id"], "kind": d["kind"], "variant": d.get("variant") or None, "status": d["status"],
            "reason": d["reason"], "seconds": d.get("seconds"),
            "peakMemoryBytes": d.get("peakBytes"), "exportFiles": d.get("exportFiles"),
            "exportBytes": d.get("exportBytes"), "source": str(d["src"]),
            "hint": d["hint"], "problems": d["problems"], "notes": d["notes"],
            "steps": {k: {x: v for x, v in s.items() if x != "tail"} for k, s in d["steps"].items()},
            "errorTail": error_tail(d)}


def write_reports(stem: Path, rows: list[dict], meta: dict) -> list[Path]:
    """<stem>.md (to read), .json and .csv (for scripts). Everything goes through scrub()."""
    md = [f"# SODA batch export, {meta['started']}", ""]
    for k in ("command", "source", "cache", "export", "options", "logs"):
        if meta.get(k):
            md.append(f"- {k}: `{meta[k]}`")
    md += [f"- finished: {meta['finished']} ({secs(meta['seconds'])}), exit code {meta['exitCode']}",
           f"- datasets: {meta['counts']}", ""]
    md += ["| id | kind | status | time | peak memory | export size | note |", "|---|---|---|---|---|---|---|"]
    for d in rows:
        md.append("| " + " | ".join(str(x).replace("|", "/") for x in (
            d["id"], kind_text(d), d["status"], secs(d.get("seconds")) if d["status"] != "skipped" else "",
            gb(d.get("peakBytes")), size_text(d.get("exportFiles"), d.get("exportBytes", 0)), d["reason"])) + " |")
    failed = [d for d in rows if d["status"] == "failed"]
    if failed:
        md += ["", "## Failed", ""]
        for d in failed:
            st = failed_step(d)
            md += [f"### {d['id']}: {d['reason']}", ""]
            if st.get("log"):
                md.append(f"full log: `{st['log']}`")
            md += ["", "last lines:", "", "```"] + (error_tail(d) or ["(no output)"]) + ["```", ""]
    skipped = [d for d in rows if d["status"] == "skipped" and d["problems"]]
    if skipped:
        md += ["", "## Skipped: needs a choice or a fix", ""]
        for d in skipped:
            md += [f"### {d['id']} ({kind_text(d)})", f"folder: `{d['src']}`", ""]
            md += [f"- {p}" for p in d["problems"]]
            if d["hint"]:
                md += ["", "```"] + d["hint"] + ["```"]
            md.append("")
    noted = [d for d in rows if d["notes"] and d["status"] == "ok"]
    if noted:
        md += ["", "## Notes", ""]
        md += [f"- {d['id']}: {n}" for d in noted for n in d["notes"]]
    md += ["", "## After the datasets", "", f"- catalog: {meta['catalog']['text']}",
           f"- sample ID check: {meta['deid']['text']}"]
    md += [f"- warning: {w}" for w in meta.get("warnings", [])]
    paths = [stem.with_suffix(".md"), stem.with_suffix(".json"), stem.with_suffix(".csv")]
    paths[0].write_text(scrub("\n".join(md) + "\n"), encoding="utf-8")
    paths[1].write_text(json.dumps(scrub_obj({**meta, "datasets": [row_json(d) for d in rows]}), indent=1,
                                   ensure_ascii=False) + "\n", encoding="utf-8")
    with open(paths[2], "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["id", "kind", "variant", "status", "reason", "seconds", "peak_memory_mb", "export_files",
                    "export_mb", "log", "hint", "error_tail"])
        for d in rows:
            st = failed_step(d)
            w.writerow([scrub(str(x)) for x in (
                d["id"], d["kind"], d.get("variant") or "", d["status"], d["reason"] or "; ".join(d["problems"]),
                d.get("seconds") if d.get("seconds") is not None else "",
                round(d["peakBytes"] / 2**20) if d.get("peakBytes") else "",
                d.get("exportFiles") or "", round(d["exportBytes"] / 2**20, 1) if d.get("exportBytes") else "",
                st.get("log", ""), " / ".join(d["hint"]), " | ".join(error_tail(d)))])
    return paths


# ---------------------------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="prepare.py batch", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--src", required=True, help="the folder holding one folder (or file) per dataset")
    ap.add_argument("--out", required=True, help="the cache folder (stays on this computer)")
    ap.add_argument("--dest", required=True, help="the export folder (what you upload)")
    ap.add_argument("--options", help="JSON file: settings per dataset id, the notebook's OPTIONS keys")
    ap.add_argument("--jobs", type=int, default=1,
                    help="datasets at once (default 1). Each whole-section dataset may need about "
                         f"{WHOLE_SECTION_GB} GB of memory, so N jobs may need N times that")
    ap.add_argument("--rebuild", action="store_true",
                    help="prepare and export again even when the cache / export already has the dataset "
                         "(old tiles stay: delete <export>/<id>/ first for a clean export)")
    ap.add_argument("--only", help="comma list of dataset ids: run just these")
    ap.add_argument("--dry-run", action="store_true", help="say what would run, run nothing, write nothing")
    ap.add_argument("--fmt", choices=["jpg", "png"], default="jpg",
                    help="image tiles: jpg (default, as the notebook: several times smaller) or png (lossless)")
    ap.add_argument("--workers", type=int, default=8, help="tile writers per export (default 8)")
    ap.add_argument("--gene-files", choices=["pack", "loose", "both"],
                    help="per-gene data: pack (prepare.py's default: a few ~64 MB packs per layer, read with HTTP "
                         "Range requests), loose (one file per gene, for a host that ignores Range) or both")
    ap.add_argument("--shard-mb", type=float, help="gene pack size in MB (prepare.py's default 64)")
    ap.add_argument("--bins", help="Visium HD layers to build, e.g. 008um,016um (the default) or all")
    ap.add_argument("--markers-genes", help="marker list for gene data (default markers/rna_generic.json)")
    ap.add_argument("--markers-proteins", help="marker list for protein images (default markers/protein_generic.json)")
    ap.add_argument("--title", help="a name for the whole collection, for datasets.json (catalog.py --title)")
    ap.add_argument("--key", help="the de-identification key file (default: SPATIALVIZ_KEY, else "
                                  "~/private/soda/sample_id_key.csv as in the notebook). When it exists the "
                                  "sample ID check is a gate")
    ap.add_argument("--check-ids", action="store_true",
                    help="make the sample ID check a gate even without a key file (the notebook's ALWAYS_CHECK_IDS)")
    ap.add_argument("--no-ids", action="store_true",
                    help="write no barcode / cell id lists (expr/<key>/ids.txt.gz). By default each export holds them; "
                         "a layer whose ids fail the id-format or sample-ID check is left out with a message")
    ap.add_argument("--logs", help="folder for the logs (default: 'logs' next to the export folder)")
    ap.add_argument("--threads", type=int,
                    help="BLAS / OpenMP threads per dataset (default: as prepare.py, split between --jobs)")
    return ap


def main(argv=None) -> int:
    for stream in (sys.stdout, sys.stderr):      # a redirected Windows console (cp1252) must not stop the run
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(errors="replace")
    a = build_parser().parse_args(argv)
    if a.jobs < 1 or a.workers < 1 or (a.threads is not None and a.threads < 1):
        print("[batch] --jobs, --workers and --threads must be 1 or more", file=sys.stderr)
        return 2
    src, a.out, a.dest = Path(a.src).resolve(), Path(a.out).resolve(), Path(a.dest).resolve()
    for p in (a.markers_genes, a.markers_proteins):
        if p and not Path(p).exists():
            print(f"[batch] marker list not found: {p}", file=sys.stderr)
            return 2
    a.markers_genes = str(Path(a.markers_genes).resolve()) if a.markers_genes else None
    a.markers_proteins = str(Path(a.markers_proteins).resolve()) if a.markers_proteins else None
    if not src.is_dir():
        print(f"[batch] --src is not a folder: {src}", file=sys.stderr)
        return 2
    logs = Path(a.logs).resolve() if a.logs else a.dest.parent / "logs"
    if logs == a.dest or a.dest in logs.parents:
        print(f"[batch] keep the logs out of the export folder (they hold local paths): {logs}", file=sys.stderr)
        return 2

    options, warnings = load_options(a.options)
    datasets, more = find_datasets(src, options, a.bins)
    warnings += more
    if a.only:
        want = [x.strip() for x in a.only.split(",") if x.strip()]
        missing = [x for x in want if x not in {d["id"] for d in datasets}]
        if missing:
            warnings.append(f"--only: no dataset {', '.join(missing)} under --src")
        datasets = [d for d in datasets if d["id"] in want]
    for w in warnings:
        say(scrub(f"[batch] WARNING {w}"))
    if not datasets:
        say(f"[batch] no datasets under {src}. Put one folder per dataset there.")
        return 2

    ready = [d for d in datasets if not d["problems"]]
    for d in datasets:
        if d["problems"]:
            d["status"], d["reason"] = "skipped", d["problems"][0]
    say(scrub(f"[batch] {len(datasets)} dataset(s) under {src}: {len(ready)} ready, "
              f"{len(datasets) - len(ready)} need a choice or a fix"))
    for d in datasets:
        say(scrub(f"  {'READY' if not d['problems'] else 'SKIP '} {kind_text(d):34s} {d['id']}"))
        for line in (d["problems"] + d["hint"]) if d["problems"] else []:
            say(scrub(f"        {line}"))

    if a.dry_run:
        for d in ready:
            done = (a.dest / d["id"] / "manifest.json").exists() and not a.rebuild
            say(scrub(f"\n{d['id']}: " + ("already exported, would be skipped" if done else "would run:")))
            if not done:
                say(scrub("  python prepare.py " + " ".join(detect._quote(x) for x in prepare_argv(d, a))))
                if d["kind"] != "rnaseq":
                    say(scrub("  python prepare.py " + " ".join(detect._quote(x) for x in export_argv(d, a))))
        say("\n[batch] dry run: nothing was run or written")
        return 0

    # the de-identification key, as the notebook: --key, else SPATIALVIZ_KEY, else the private default
    key = Path(a.key or os.environ.get("SPATIALVIZ_KEY") or Path.home() / "private" / "soda" / "sample_id_key.csv")
    if a.key or os.environ.get("SPATIALVIZ_KEY") or key.exists():
        os.environ["SPATIALVIZ_KEY"] = str(key)       # every prepare.py started from here uses this key
    gate = key.exists() or a.check_ids

    mem = memory_gb()
    if mem:
        need = WHOLE_SECTION_GB * a.jobs
        say(f"[batch] memory: {mem[0]:.1f} GB in total, {mem[1]:.1f} GB free. A whole-section dataset can need "
            f"about {WHOLE_SECTION_GB} GB" + (f", so {a.jobs} at once about {need:.0f} GB" if a.jobs > 1 else ""))
        if a.jobs > 1 and mem[1] < need:
            say(f"[batch] WARNING only {mem[1]:.1f} GB free for {a.jobs} jobs: big datasets may run out of memory "
                "(each then fails on its own and is logged). Use --jobs 1 for whole sections.")
    threads = a.threads
    if threads is None and a.jobs > 1:
        import threadcap
        threads = max(1, threadcap.cap() // a.jobs)

    for p in (a.out, a.dest, logs):
        p.mkdir(parents=True, exist_ok=True)
    started = datetime.now()
    stem = logs / f"batch_{started.strftime('%Y%m%d_%H%M%S')}"
    i = 1
    while stem.with_suffix(".json").exists():
        stem = logs / f"batch_{started.strftime('%Y%m%d_%H%M%S')}_{i}"
        i += 1
    say(f"[batch] {len(ready)} to run, {a.jobs} at a time; logs in {logs}")

    t0 = time.time()
    running, stop = set(), threading.Event()
    interrupted = False
    with ThreadPoolExecutor(max_workers=a.jobs) as ex:
        futs = [ex.submit(process, d, a, n, len(ready), logs, threads, running, stop) for n, d in enumerate(ready, 1)]
        try:
            for f in as_completed(futs):
                f.result()
        except KeyboardInterrupt:
            interrupted = True
            stop.set()
            say("\n[batch] stopping: Ctrl+C")
            for f in futs:
                f.cancel()
            for p in list(running):
                try:
                    p.kill()
                except OSError:
                    pass
    for d in ready:
        if d["status"] is None:
            d["status"], d["reason"] = "failed", "not run (stopped with Ctrl+C)"

    # datasets.json (catalog.py), as the notebook's Section 4
    cat = {"exit": None, "text": "not run: no dataset in the export folder"}
    if any((p / "manifest.json").exists() for p in a.dest.iterdir() if p.is_dir()):
        cmd = [sys.executable, str(CATALOG), str(a.dest)] + (["--title", a.title] if a.title else [])
        r = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8", errors="replace",
                           env={**os.environ, "PYTHONIOENCODING": "utf-8"})
        lines = [scrub(x) for x in (r.stdout + r.stderr).strip().splitlines()]
        for x in lines:
            say(x)
        cat = {"exit": r.returncode, "lines": lines,
               "text": (f"{a.dest / 'datasets.json'} written" if r.returncode == 0 else
                        f"FAILED (exit {r.returncode}): " + (lines[-1] if lines else "no output"))}

    # the sample ID check (the notebook's Section 5) on the export and this run's own logs
    run_logs = sorted({Path(s["log"]) for d in ready for s in d["steps"].values() if s.get("log")})
    probs = deid.check_paths([a.dest])
    log_probs = deid.check_paths(run_logs) if run_logs else []
    clean = not probs and not log_probs
    if clean:
        dtext = "clean" + (" (a gate: key file found)" if key.exists() else " (a gate: --check-ids)" if gate
                           else f" (no key file at {key}, so not a gate)")
    else:
        dtext = (f"NOT clean: {len(probs)} in the export, {len(log_probs)} in the logs"
                 + (" (a gate: fix them before uploading)" if gate else
                    f" (no key file at {key}, so a warning: if your lab de-identifies samples, fix them)")
                 + f". To see them: python {HERE / 'deid.py'} check {a.dest}")
    deid_info = {"clean": clean, "gate": gate, "keyFile": str(key), "keyFound": key.exists(),
                 "exportHits": len(probs), "logHits": len(log_probs), "text": dtext}
    for p in (probs + log_probs)[:40]:
        say(scrub("[deid] " + p))

    failed = [d for d in datasets if d["status"] == "failed"]
    code = 3 if (gate and not clean) else 1 if (failed or (cat["exit"] not in (None, 0))) else 0
    if interrupted:
        code = 130
    counts = {s: sum(d["status"] == s for d in datasets) for s in ("ok", "failed", "skipped")}
    done_before = sum(d["status"] == "skipped" and d["reason"].startswith("already exported") for d in datasets)
    count_text = (f"{counts['ok']} ok, {counts['failed']} failed, {counts['skipped']} skipped"
                  + (f" ({done_before} already exported, {counts['skipped'] - done_before} need a choice or a fix)"
                     if counts["skipped"] else ""))
    meta = {"started": started.isoformat(timespec="seconds"), "finished": datetime.now().isoformat(timespec="seconds"),
            "seconds": round(time.time() - t0, 1), "exitCode": code,
            "command": "python prepare.py batch " + " ".join(detect._quote(x) for x in (argv if argv is not None else sys.argv[1:])),
            "source": str(src), "cache": str(a.out), "export": str(a.dest), "logs": str(logs),
            "options": str(Path(a.options).resolve()) if a.options else "", "jobs": a.jobs, "rebuild": a.rebuild,
            "counts": count_text, "catalog": cat, "deid": deid_info, "warnings": warnings}
    paths = write_reports(stem, datasets, meta)
    leaked = deid.check_paths(paths)               # the scrub above makes this a formality; never skip it
    if leaked:
        say(f"[batch] *** the log itself holds {len(leaked)} sample ID(s): delete {stem}.* before sharing it")
        code = 3

    say("")
    say(f"[batch] {count_text}, in {secs(meta['seconds'])}" + (" (stopped with Ctrl+C)" if interrupted else ""))
    for d in failed:
        st = failed_step(d)
        say(scrub(f"[batch]   FAILED  {d['id']}: {d['reason']}" + (f"; log: {st['log']}" if st.get("log") else "")))
    for d in datasets:
        if d["status"] == "skipped" and d["problems"]:
            say(scrub(f"[batch]   skipped {d['id']}: {d['problems'][0]}"))
            for h in d["hint"][:4]:
                say(scrub(f"[batch]           {h}"))
    say(scrub(f"[batch] catalog: {cat['text']}"))
    say(scrub(f"[batch] sample ID check: {dtext}"))
    if a.rebuild:
        say("[batch] note: --rebuild leaves old tiles in place; delete <export>/<id>/ first for a clean export")
    say(f"[batch] export: {a.dest}")
    say(f"[batch] cache:  {a.out}")
    say(f"[batch] log:    {paths[0]}")
    say(f"[batch]         {paths[1]}")
    say(f"[batch]         {paths[2]}")
    say(f"[batch] exit code {code}")
    return code


if __name__ == "__main__":
    if sys.argv[1:2] == [CHILD]:
        run_child(sys.argv[2:])
    else:
        sys.exit(main())
