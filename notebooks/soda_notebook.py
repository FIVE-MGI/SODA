"""The code behind notebooks/SODA_pipeline.ipynb: one function per step, so the notebook itself is
only the settings you edit and one call per section.

    import soda_notebook as soda
    soda.setup()                             # 0. the repository, Python and the packages
    soda.set_work(WORK, CLOUD_URL)           # 0. the work folder and the cloud folder's link
    soda.advanced(...)                       # 0. the optional settings (all have defaults)
    soda.list_inputs()                       # 1. what is in WORK/inputs
    soda.check_inputs(OPTIONS)               # 2. what each dataset is (batch.find_datasets)
    soda.inspect_images()                    # 2. optional: look inside proteomics / H&E files
    soda.memory_check()                      # 3. memory before the slow step (only warns)
    soda.prepare_all()                       # 3. prepare.py, per dataset, into cache/
    soda.export_all()                        # 4. prepare.py export into export/, then catalog.py
    soda.deid_check()                        # 5. deid.py check on export/
    soda.sizes(); soda.preview(); soda.stop_preview()   # 6. look before you upload
    soda.share(IS_THIS_DATA_PUBLIC)          # 7. the upload command and the link to share
    soda.test_link()                         # 7. after the upload: can the SODA page read it?

It runs the same scripts with the same arguments as the command lines in docs/PREPARING_DATA.md:
web/SpatialViewer/prepare.py, catalog.py, deid.py and static_serve.py. Finding the datasets is
batch.find_datasets (the same code as `prepare.py batch` and the desktop app) and what each folder
is, is detect.py: there is no file finding of its own here.

Only the standard library is imported at the top, so this module loads before setup() has
checked the packages; detect.py and batch.py are imported when Section 2 runs. What one step
leaves for the next (the datasets found, the prepare results, the sample ID check's answer) is
kept in `state`.
"""
import atexit
import html
import importlib
import json
import os
import re
import socket
import subprocess
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path
from types import SimpleNamespace
from urllib.parse import quote, urlsplit

REPO = Path(__file__).resolve().parent.parent
SV = REPO / "web" / "SpatialViewer"          # the scripts: prepare.py, catalog.py, deid.py, static_serve.py
PY = sys.executable                          # the same Python this notebook runs on
if not (SV / "prepare.py").exists():
    raise FileNotFoundError("soda_notebook.py must stay in the notebooks/ folder of the SODA repository "
                            "(the folder that has web/SpatialViewer/prepare.py).")
if str(SV) not in sys.path:
    sys.path.insert(0, str(SV))

SODA_PAGE = "https://soda.mgifive.org/"   # the shared SODA page (GitHub Pages of FIVE-MGI/SODA, custom domain)
PAGE_PLACEHOLDER = "<the SODA page's address>"   # shown in the link when SODA_PAGE is set to ""
DEFAULT_KEY_FILE = Path(os.environ.get("SPATIALVIZ_KEY") or Path.home() / "private" / "soda" / "sample_id_key.csv")

# The settings (the notebook's Section 0 sets them) and what each step leaves for the next.
state = SimpleNamespace(
    work=None, inputs=None, cache=None, export=None, cloud_url="",
    key_file=DEFAULT_KEY_FILE, always_check_ids=False, soda_page=SODA_PAGE, title="", read_pass="",
    tile_format="jpg", visium_bins="008um,016um", workers=8, markers={"visium": "", "proteomics": ""},
    rebuild=False, include_ids=True,
    options={}, datasets=None, ready=None, results=None, deid_status=None, servers=[], link=None)

# this notebook's names for three of prepare.py's kinds (also the folder names of the older inputs/<kind>/ layout)
NB_KIND = {"tiff": "proteomics_file", "akoya": "proteomics_markers", "rnaseq": "scrnaseq"}
GROUP = {"proteomics_file": "proteomics", "proteomics_markers": "proteomics", "visium": "transcriptomics",
         "xenium": "transcriptomics", "stereoseq": "transcriptomics", "cosmx": "transcriptomics", "merscope": "transcriptomics", "scrnaseq": "scrnaseq", "geomx": "transcriptomics",
         "anndata": "transcriptomics", "spatialdata": "transcriptomics", "omezarr": "proteomics"}   # picker headings (Section 7)
HOW_SURE = {"certain": "certain", "ambiguous": "NOT SURE: pick one (see To do)",
            "incomplete": "recognised, but something is missing", "unknown": "not recognised"}
# batch.py's messages speak of its command line; here the same things are OPTIONS, inputs/ and VISIUM_BINS
WORDS = (("(options file)", "(OPTIONS)"), ("the options file", "OPTIONS"), ("under --src", "under inputs/"),
         ("--bins asks", "VISIUM_BINS asks"))


def _need(value, message):
    if value is None:
        raise RuntimeError(message)
    return value


def _need_work():
    _need(state.work, "Run the Settings cell (Section 0) first: it sets WORK.")


# ---------------------------------------------------------------------------------------------
# 0. Setup
# ---------------------------------------------------------------------------------------------
def check_packages():
    """Print ok / MISSING for every package SODA needs, and WRONG VERSION for zarr 3 or a new tifffile."""
    need = {"numpy": "core", "tifffile": "core", "zarr": "core", "imagecodecs": "core", "PIL": "core",
            "h5py": "core", "pyarrow": "core", "scipy": "core", "cryptography": "core",
            "sklearn": "UMAP (optional)", "umap": "UMAP (optional)", "igraph": "UMAP (optional)",
            "leidenalg": "UMAP (optional)"}
    missing = []
    for mod, what in need.items():
        try:
            m = importlib.import_module(mod)
            print(f"  ok       {mod:12s} {getattr(m, '__version__', '')}")
        except ImportError:
            print(f"  MISSING  {mod:12s} ({what})")
            missing.append(mod)

    def major(mod):
        try:
            return int(importlib.import_module(mod).__version__.split(".")[0])
        except Exception:
            return None
    wrong = []
    if (major("zarr") or 0) >= 3:
        wrong.append("zarr must be 2.x (found 3+): images cannot be written")
    if (major("tifffile") or 0) >= 2025:
        wrong.append("tifffile must be older than 2025 (it needs zarr 3)")
    for w in wrong:
        print("  WRONG VERSION:", w)
    if wrong or any(need[m] == "core" for m in missing):
        print("\nThis Python is not set up for SODA. Make the soda-env above, install the requirements in it,"
              "\nand run this notebook with that environment. Stop here until this cell shows no MISSING core or WRONG lines.")


def setup():
    """Section 0: where the repository is, which Python runs, and its packages."""
    print("repository:", REPO)
    print("python:    ", PY, sys.version.split()[0])
    check_packages()


def _make_folders():
    state.inputs, state.cache, state.export = state.work / "inputs", state.work / "cache", state.work / "export"
    for d in (state.inputs, state.cache, state.export):
        d.mkdir(parents=True, exist_ok=True)
    os.environ["SPATIALVIZ_KEY"] = str(state.key_file)  # every script started from here uses this key


def set_work(work, cloud_url=""):
    """Section 0, Settings: WORK (inputs/, cache/ and export/ go inside) and CLOUD_URL (may stay empty until Section 7)."""
    state.work = Path(work)
    state.cloud_url = cloud_url or ""
    _make_folders()
    print("work folder: ", state.work)
    print("cloud folder:", state.cloud_url.strip().split("?", 1)[0] or "(not set yet: paste it before Section 7)")


def advanced(SODA_PAGE=SODA_PAGE, KEY_FILE=DEFAULT_KEY_FILE, ALWAYS_CHECK_IDS=False, TITLE="", READ_PASS="",
             TILE_FORMAT="jpg", VISIUM_BINS="008um,016um", WORKERS=8, MARKERS=None, REBUILD=False, INCLUDE_IDS=True):
    """Section 0, Advanced settings (optional). The defaults here are the notebook's defaults; see
    docs/ADVANCED.md, "Every notebook setting"."""
    _need_work()
    state.soda_page, state.key_file, state.always_check_ids = SODA_PAGE, Path(KEY_FILE), ALWAYS_CHECK_IDS
    state.title, state.read_pass = TITLE, READ_PASS
    state.tile_format, state.visium_bins, state.workers = TILE_FORMAT, VISIUM_BINS, WORKERS
    state.markers = dict(MARKERS) if MARKERS is not None else {"visium": "", "proteomics": ""}
    state.rebuild, state.include_ids = REBUILD, INCLUDE_IDS
    _make_folders()
    key = state.key_file
    print("key file:   ", key, "(found: Section 5 checks for sample IDs)" if key.exists()
          else "(none: Section 5 skips the sample ID check)")
    if key.resolve().is_relative_to(state.work.resolve()) or key.resolve().is_relative_to(REPO.resolve()):
        print("WARNING: keep the key file outside the work folder and the repository.")
    print("SODA page:  ", SODA_PAGE if SODA_PAGE.strip() and "<" not in SODA_PAGE else "(not set: Section 7 prints the link with a placeholder)")


# ---------------------------------------------------------------------------------------------
# 1. Put your files here
# ---------------------------------------------------------------------------------------------
def list_inputs():
    """Section 1: write inputs/README.txt once and list what is in inputs/."""
    _need_work()
    inputs = state.inputs
    readme = inputs / "README.txt"
    if not readme.exists():
        readme.write_text("Put each dataset in its own folder here. Any name works; the name becomes the dataset id.\n"
                          "SODA works out what each folder holds: run Section 2 of notebooks/SODA_pipeline.ipynb.\n"
                          "What each instrument's folder should contain: docs/PREPARING_DATA.md.\n", encoding="utf-8")
    here = sorted(p.name + ("/" if p.is_dir() else "") for p in inputs.iterdir() if p.name != "README.txt")
    print(f"{inputs}:")
    print("  " + ("\n  ".join(here) if here else "empty. Put one folder per dataset in it, then run Section 2."))


# ---------------------------------------------------------------------------------------------
# 2. Check the inputs
# ---------------------------------------------------------------------------------------------
def nb_text(s):
    for a, b in WORDS:
        s = s.replace(a, b)
    return s


def nb_row(r):
    # One dataset as batch.py found it, in this notebook's words: its kind name, the files found (d["found"],
    # also your own UMAP / clusters as found["user"]), the prepare.py arguments (d["argv"]) and the To do list.
    import detect
    rep = r["report"]
    c = rep["candidates"][0] if rep["candidates"] else {"args": {}, "argv": [], "variant": "", "notes": []}
    found = {k: [Path(x) for x in v] if isinstance(v, list) else Path(v) for k, v in c["args"].items()}
    if r["user"]:
        found["user"] = r["user"]
    if rep["verdict"] != "certain":       # what detection says, then the OPTIONS line(s) that would settle it
        todo = [rep["message"]] + detect.option_lines(rep, r["id"]) + r["problems"][1:]
        notes = list(c["notes"]) + list(rep["extras"])
    else:
        todo, notes = r["problems"], r["notes"]
    return dict(id=r["id"], src=r["src"], report=rep, kind=NB_KIND.get(rep["kind"], rep["kind"]) or "?",
                variant=c["variant"], found=found, argv=list(c["argv"]),
                problems=[nb_text(p) for p in todo], notes=[nb_text(n) for n in notes])


def find_datasets():
    # Every dataset under inputs/: batch.find_datasets, the one finder that `prepare.py batch` and the desktop
    # app use too. Names starting with _ or . are left out; x.zip beside a folder x/ becomes a second dataset, x_zip.
    import batch
    rows, warnings = batch.find_datasets(state.inputs, state.options, state.visium_bins)
    for w in warnings:
        print("Note:", nb_text(w))
    return [nb_row(r) for r in rows]


def _in_ipython():
    try:
        from IPython import get_ipython
        return get_ipython() is not None
    except Exception:
        return False


def show(d):
    """One dataset's table: what it is, how sure, what was found, what to do."""
    rep = d["report"]
    c = rep["candidates"][0] if rep["candidates"] else None
    found = [f"{e['what']}: {e['where']}" + (f" ({e['detail']})" if e["detail"] else "") for e in (c["evidence"] if c else [])]
    rows = [("Dataset id", d["id"]), ("Folder", str(d["src"])),
            ("Detected", rep["label"] or ", ".join(x["label"] for x in rep["candidates"]) or "-"),
            ("How sure", HOW_SURE[rep["verdict"]]), ("Found", found or ["-"]),
            ("To do", d["problems"] or ["nothing: ready"]), ("Notes", d["notes"] or ["-"])]
    if _in_ipython():
        from IPython.display import display, HTML
        cell = lambda v: "<br>".join(html.escape(str(x)) for x in (v if isinstance(v, list) else [v]))  # noqa: E731
        display(HTML("<table>" + "".join(
            f"<tr><th style='text-align:left;vertical-align:top;white-space:nowrap'>{k}</th>"
            f"<td style='text-align:left'>{cell(v)}</td></tr>" for k, v in rows) + "</table>"))
    else:
        for k, v in rows:
            for i, x in enumerate(v if isinstance(v, list) else [v]):
                print(f"  {k if i == 0 else '':11s} {x}")
        print()


def check_inputs(options=None):
    """Section 2: look at every folder under inputs/ (nothing is converted), one table per dataset,
    then the READY / FIX list. `options` is the notebook's OPTIONS dictionary."""
    _need(state.inputs, "Run the Settings cell (Section 0) first: it sets WORK.")
    import batch
    import detect
    importlib.reload(detect)
    importlib.reload(batch)
    state.options = dict(options or {})
    state.datasets = datasets = find_datasets()
    if not datasets:
        print(f"No datasets yet. Put one folder per dataset under {state.inputs} (Section 1), then run this cell again.")
    for d in datasets:
        show(d)
    state.ready = ready = [d for d in datasets if not d["problems"]]
    for d in datasets:
        print(f"{'READY' if not d['problems'] else 'FIX  '} {d['kind']:19s} {d['id']}")
    print(f"\n{len(ready)} of {len(datasets)} dataset(s) ready." + ("" if len(ready) == len(datasets)
          else " For the others: see To do in their table, change the folder or OPTIONS, and run this cell again."))


def _ready():
    return _need(state.ready, "Run Section 2 (Check the inputs) first.")


def inspect_images():
    """Section 2, optional: axes, size, levels and channel names of each proteomics / H&E file. Writes nothing."""
    for d in _ready():
        if d["kind"] == "proteomics_file":
            subprocess.run([PY, str(SV / "prepare.py"), "inspect", "--src", str(d["found"]["src"])])


# ---------------------------------------------------------------------------------------------
# 3. Prepare
# ---------------------------------------------------------------------------------------------
def memory_gb():
    """(total, free) memory in GB (2**30 bytes), or None when this system does not say."""
    try:
        import psutil
        vm = psutil.virtual_memory()
        return vm.total / 2**30, vm.available / 2**30
    except ImportError:
        pass
    import batch                          # the same check without psutil: Windows, Linux, macOS
    return batch.memory_gb()


def files_size(p):
    p = Path(p)
    if p.is_file():
        return p.stat().st_size
    if p.is_dir():
        return sum(q.stat().st_size for q in p.rglob("*") if q.is_file())
    return 0


def input_gb(d):
    # About the bytes prepare.py reads for this dataset: for Visium only the image and the count
    # matrices it uses; for a folder input every file in it (an upper bound)
    f, bins = d["found"], state.visium_bins
    if d["kind"] == "visium" and "outs" in f:
        paths = [f.get("image")]
        binned = f["outs"] / "binned_outputs"
        if binned.is_dir():
            for layer in binned.glob("square_*um"):
                if not bins or layer.name[len("square_"):] in bins.split(","):
                    paths += [layer / "filtered_feature_bc_matrix.h5", layer / "spatial" / "tissue_positions.parquet"]
        else:
            paths += [f["outs"] / "filtered_feature_bc_matrix.h5"]
    else:
        paths = [p for k, v in f.items() if k != "user" for p in (v if isinstance(v, list) else [v])]
    return sum(files_size(p) for p in paths if isinstance(p, (str, Path))) / 2**30


def memory_check():
    """Section 3, before the slow step: total and free memory, the biggest dataset's input files, and the
    measured peaks from the README. It only warns; nothing is stopped."""
    ready, bins = _ready(), state.visium_bins
    try:
        mem = memory_gb()
        if mem:
            print(f"memory: {mem[0]:.1f} GB in total, {mem[1]:.1f} GB free now")
        else:
            print("memory: this system does not report it here; compare with the numbers below yourself")
        sizes = sorted(((input_gb(d), d) for d in ready), key=lambda t: t[0], reverse=True)
        if sizes:
            gb, big = sizes[0]
            print(f"largest input: {big['id']} ({big['kind']}), about {gb:.2f} GB of input files")
            # Measured peaks (README, Requirements): a CosMx slide's counts 0.76 GB;
            # a whole mouse brain Stereo-seq chip 2.7 GB to prepare (3.9 GB from tissue.gem.gz) and
            # 8.5 GB to export, the largest measured. Whole-section data: have more free than 8.5 GB.
            whole = gb >= 1.0                       # 1 GB of input or more: treat as a whole section
            need = 8.5 if whole else 1.0
            print(f"suggested free memory: more than {need} GB "
                  + ("(whole-section data: the largest measured peak, a whole Stereo-seq chip, was 8.5 GB)"
                     if whole else "(small data: a CosMx slide's counts peaked at 0.76 GB)"))
            hd2 = [d["id"] for d in ready if d["kind"] == "visium" and "outs" in d["found"]
                   and (d["found"]["outs"] / "binned_outputs" / "square_002um").is_dir()
                   and (not bins or "002um" in bins.split(","))]
            if hd2:
                print("note: the 2 um Visium HD layer is on for " + ", ".join(hd2) + "; it has not been measured "
                      "and is expected to need the most memory. Take 002um out of VISIUM_BINS if it runs out.")
            if mem and mem[1] < need:
                print(f"\nWARNING: only {mem[1]:.1f} GB free, less than the {need} GB suggested. Prepare or export may "
                      "run out of memory (a MemoryError, or a very slow computer). Close other programs, use a "
                      "computer with more memory, or prepare the biggest dataset on its own. Nothing is stopped: "
                      "run the next cell when ready.")
            elif mem:
                print("ok: more free memory than suggested")
        else:
            print("no ready datasets (Section 2)")
    except Exception as e:                                             # a check, never a blocker
        print("memory check skipped:", e)


def run(args, logfile=None):
    """Run one script, print its output as it comes, keep a log. Returns True on success."""
    print("$", " ".join(str(a) for a in args))
    t0 = time.time()
    log = open(logfile, "a", encoding="utf-8") if logfile else None
    p = subprocess.Popen([str(a) for a in args], cwd=str(REPO), stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace",
                         env={**os.environ, "PYTHONUNBUFFERED": "1", "PYTHONIOENCODING": "utf-8"})
    for line in p.stdout:
        print(line, end="")
        if log:
            log.write(line)
    p.wait()
    if log:
        log.close()
    print(f"-> exit {p.returncode} in {time.time() - t0:.0f} s\n")
    return p.returncode == 0


def prepare_args(d):
    """The prepare.py command for one dataset: Section 2 found the files (d["argv"] is the kind and its
    file arguments); OPTIONS add the rest."""
    o, k = state.options.get(d["id"], {}), d["kind"]
    a = [PY, SV / "prepare.py"] + d["argv"] + ["--id", d["id"], "--out", state.cache, "--label", o.get("label", d["id"])]

    def opt(key, *kinds):
        v = o.get(key)
        if (not kinds or k in kinds) and v is not None and v is not False and v != "":
            a.extend(["--" + key.replace("_", "-")] + ([] if v is True else [v]))
    opt("pixel_size")
    opt("channels", "proteomics_file", "omezarr")
    if k == "visium" and d["variant"] == "visium-hd":
        a += ["--bins", state.visium_bins or "all"]      # "" = every layer (prepare.py alone keeps 8 and 16 um)
    opt("image_scale", "visium")
    opt("no_cells", "xenium", "cosmx", "merscope")
    opt("no_transcripts", "xenium", "cosmx", "merscope", "spatialdata")
    opt("qv_min", "xenium")
    opt("images", "cosmx", "merscope")
    opt("z", "cosmx", "merscope", "omezarr", "spatialdata")
    # standard formats (stdformats.py)
    opt("image_scale", "anndata")
    for key in ("library", "spatial_key", "image_res", "layer_kind"):
        opt(key, "anndata")
    for key in ("matrix", "coords_um", "no_user_clusters"):
        opt(key, "anndata", "spatialdata")
    for key in ("table", "region", "coordinate_system", "image_element", "points"):
        opt(key, "spatialdata")
    opt("series", "omezarr")
    opt("bins", "stereoseq")
    for key in ("sheet", "layout", "keep_slide_names"):
        opt(key, "geomx")
    for key in ("crop", "flip"):
        opt(key, "proteomics_file", "proteomics_markers", "visium")
    for key in ("cells", "nuclei", "cells_units"):   # a lab's own cell segmentation of a protein image (--cells)
        opt(key, "proteomics_file", "proteomics_markers", "omezarr")
    return a


def rnaseq_args(d):
    """The prepare.py command for single-cell RNA-seq: straight into export/ (no cache)."""
    o = state.options.get(d["id"], {})
    a = [PY, SV / "prepare.py"] + d["argv"] + ["--id", d["id"], "--label", o.get("label", d["id"]), "--dest", state.export]
    if o.get("embedding"):                             # OPTIONS: which obsm key is the map (else X_umap, a t-SNE, then PCA)
        a += ["--embedding", o["embedding"]]
    return a


def prepare_all():
    """Section 3: prepare every ready dataset into cache/ (single-cell RNA-seq straight into export/).
    A dataset already there is skipped unless REBUILD = True."""
    ready = _ready()
    (state.work / "logs").mkdir(exist_ok=True)
    state.results = results = {}
    for d in ready:
        log = state.work / "logs" / f"{d['id']}_prepare.log"
        if d["kind"] == "scrnaseq":
            if (state.export / d["id"] / "manifest.json").exists() and not state.rebuild:
                print(f"skip {d['id']}: already in export/"); results[d["id"]] = True; continue  # noqa: E702
            results[d["id"]] = run(rnaseq_args(d), log)
        else:
            if (state.cache / d["id"] / "manifest.json").exists() and not state.rebuild:
                print(f"skip {d['id']}: already in cache/"); results[d["id"]] = True; continue  # noqa: E702
            results[d["id"]] = run(prepare_args(d), log)
    print("prepared:", [k for k, v in results.items() if v])
    print("FAILED:  ", [k for k, v in results.items() if not v] or "none")


# ---------------------------------------------------------------------------------------------
# 4. Export
# ---------------------------------------------------------------------------------------------
def export_args(d):
    """The prepare.py export command for one prepared dataset."""
    a = [PY, SV / "prepare.py", "export", "--id", d["id"], "--out", state.cache, "--dest", state.export,
         "--fmt", state.tile_format, "--workers", state.workers]
    # gene data (Visium, Xenium, CosMx, MERSCOPE, Stereo-seq, GeoMx) takes the gene list,
    # protein images the protein list
    group = "proteomics" if d["kind"].startswith("proteomics") or d["kind"] == "omezarr" else "visium"
    if state.markers.get(group):
        a += ["--markers", state.markers[group]]
    if not state.include_ids:                    # barcode / cell id lists: written unless turned off
        a += ["--no-ids"]
    u, opt = d["found"].get("user"), state.options.get(d["id"], {})
    if u:                                        # your own UMAP / clusters (Section 1)
        a += ["--user-layer", opt["user_layer"]]
        a += ["--user-umap", u["umap"]] if u["umap"] else []
        for p in u["obs"]:
            a += ["--user-obs", p]
        a += ["--user-h5ad", u["h5ad"]] if u["h5ad"] else []
        for k in ("user_embedding", "user_columns", "user_barcode_prefix", "user_barcode_suffix"):
            if opt.get(k):
                a += ["--" + k.replace("_", "-"), opt[k]]
        a += ["--user-replace"] if opt.get("user_replace") else []
    return a


def export_all():
    """Section 4: export every prepared dataset into export/<id>/, then write export/datasets.json (catalog.py)."""
    ready = _ready()
    results = _need(state.results, "Run Section 3 (Prepare) first: it skips the datasets that are already in cache/.")
    exported = {}
    for d in ready:
        if d["kind"] == "scrnaseq" or not results.get(d["id"]):
            continue
        log = state.work / "logs" / f"{d['id']}_export.log"
        before = log.stat().st_size if log.exists() else 0      # run() appends: read only this run's lines
        exported[d["id"]] = run(export_args(d), log)
        with open(log, "rb") as f:
            f.seek(before)
            txt = f.read().decode("utf-8", errors="replace")
        if "***" in txt:
            print(f"WARNING {d['id']}: the export log has a *** line; read it before uploading.")
    print("exported:", [k for k, v in exported.items() if v])
    print("FAILED:  ", [k for k, v in exported.items() if not v] or "none")

    # The list of every dataset in export/ (datasets.json): the link you share opens this list.
    export = state.export
    if any((p / "manifest.json").exists() for p in export.iterdir() if p.is_dir()):
        print()
        if not run([PY, SV / "catalog.py", export] + (["--title", state.title] if state.title.strip() else [])):
            print("datasets.json was NOT written: read the [catalog] line above, fix it and run this cell again.")
    else:
        print("No dataset in export/ yet, so no datasets.json either.")


# ---------------------------------------------------------------------------------------------
# 5. De-identification check
# ---------------------------------------------------------------------------------------------
def deid_check():
    """Section 5: deid.py check on export/. With a key file (or ALWAYS_CHECK_IDS) it is a gate for Section 7;
    without one it only warns."""
    _need_work()
    key = state.key_file
    r = subprocess.run([PY, str(SV / "deid.py"), "check", str(state.export)], cwd=str(REPO),
                       capture_output=True, text=True, encoding="utf-8", errors="replace")
    found_clean = r.returncode == 0 and "[deid] clean" in r.stdout
    if key.exists() or state.always_check_ids:
        # this lab de-identifies: the check is a gate
        print(r.stdout[-5000:] or r.stderr[-2000:])
        clean = found_clean
        state.deid_status = "clean" if clean else "not clean"
        print("\nOK to upload." if clean else "\nNOT clean: fix the listed names before uploading. Section 7 stops until this says clean.")
    else:
        # no key file: the check still runs, but only warns (its patterns are one lab's ID formats)
        state.deid_status = "skipped"
        if found_clean:
            print(f"Sample ID check: nothing found (no key file at {key}, so this is not a gate).")
        else:
            print(r.stdout[-5000:] or r.stderr[-2000:])
            print(f"\nWARNING: the names above look like sample IDs. No key file at {key}, so this does not stop the upload.\n"
                  "If your lab de-identifies samples, set KEY_FILE in the advanced settings, run that cell and this one again,\n"
                  "and fix the names before uploading. If these are not sample IDs, carry on.")


# ---------------------------------------------------------------------------------------------
# 6. Look before you upload
# ---------------------------------------------------------------------------------------------
def folder_size(p):
    n = b = 0
    for f in p.rglob("*"):
        if f.is_file():
            n += 1
            b += f.stat().st_size
    return n, b


def sizes():
    """Section 6: files and GB per dataset in export/ (many hosts limit both)."""
    _need_work()
    tot_n = tot_b = 0
    for ds in sorted(p for p in state.export.iterdir() if p.is_dir()):
        n, b = folder_size(ds)
        tot_n += n
        tot_b += b
        print(f"{ds.name:40s} {n:>9,d} files {b / 1e9:8.2f} GB")
    print(f"{'TOTAL':40s} {tot_n:>9,d} files {tot_b / 1e9:8.2f} GB")


def free_port():
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


def preview():
    """Section 6: two small local servers, one for export/ (static_serve.py) and one for the viewer page,
    and the links that open them. stop_preview() stops them."""
    _need_work()
    export = state.export
    stop_preview(quiet=True)                     # a second run replaces the servers of the first
    data_port, page_port = free_port(), free_port()
    state.servers = [
        subprocess.Popen([PY, str(SV / "static_serve.py"), str(export), str(data_port)], cwd=str(SV)),
        subprocess.Popen([PY, "-m", "http.server", str(page_port), "--bind", "127.0.0.1", "--directory", str(SV)]),
    ]
    atexit.register(stop_preview, True)         # the servers end with the notebook's kernel
    time.sleep(1.5)
    local_folder = f"http://127.0.0.1:{data_port}/"
    print("All datasets, the way your shared link will open them:")
    print(f"  http://127.0.0.1:{page_port}/index.html?data={quote(local_folder, safe='')}")
    if not (export / "datasets.json").exists():
        print("  (no datasets.json yet: run Section 4 to the end first)")
    print("\nOne dataset at a time (also works when your config.js sets a password):")
    for ds in sorted(p.name for p in export.iterdir() if p.is_dir()):
        print(f"  http://127.0.0.1:{page_port}/index.html?static=http://127.0.0.1:{data_port}&ds={ds}")
    print("\nOpen a link and click through every dataset. The servers run until the notebook stops "
          "(or soda.stop_preview()).")


def stop_preview(quiet=False):
    """Section 6: stop the servers preview() started."""
    for s in state.servers:
        s.terminate()
    state.servers = []
    if not quiet:
        print("servers stopped")


# ---------------------------------------------------------------------------------------------
# 7. Upload and share
# ---------------------------------------------------------------------------------------------
def page_address(page):
    """The SODA page's index.html, or None when SODA_PAGE is not set."""
    p = (page or "").strip()
    if not p or "<" in p:
        return None
    if not p.startswith(("https://", "http://")):
        p = "https://" + p
    return p if p.split("?", 1)[0].endswith(".html") else p.rstrip("/") + "/index.html"


def share(is_this_data_public):
    """Section 7: needs the answer to "is this data public?" ("yes" or "no") and a clean (or skipped)
    Section 5. Prints the upload command, THE LINK TO SHARE and the config.js lines. Runs and uploads nothing."""
    answer = {True: "yes", False: "no"}.get(is_this_data_public, str(is_this_data_public).strip().lower())
    if answer not in ("yes", "no"):
        raise ValueError('Answer the question first: set IS_THIS_DATA_PUBLIC = "yes" or "no" at the top of this cell.')
    mode = "public" if answer == "yes" else "private"

    if state.deid_status is None:
        raise RuntimeError("Run Section 5 (the de-identification check) first.")
    if state.deid_status == "not clean":
        raise RuntimeError("Section 5 found sample IDs. Fix them and run Section 5 again until it says clean.")
    if state.deid_status == "skipped":
        print("Note: the sample ID check was skipped (no key file). Read the labels and channel names yourself.\n")
    export = state.export
    if not (export / "datasets.json").exists():
        print("Note: export/ has no datasets.json. Run Section 4 to the end before uploading.\n")

    url = state.cloud_url.strip()
    if "?" in url:                              # a pasted pass (SAS) is a secret: drop it, never print it
        url = url.split("?", 1)[0]
        print("Note: removed the '?...' part of CLOUD_URL. Paste only the folder address there "
              "(a read-only pass goes in READ_PASS).\n")
    base_url = url.rstrip("/")
    if not base_url.startswith(("https://", "http://")) or "<" in base_url:
        raise ValueError("Set CLOUD_URL (Section 0) to your cloud folder's link, e.g. https://myaccount.blob.core.windows.net/spatial")
    on_azure = ".blob.core.windows.net/" in base_url + "/"

    read_pass = state.read_pass.strip().lstrip("?")
    folder_link = base_url + "/" + ("?" + read_pass if read_pass else "")
    page = page_address(state.soda_page)
    share_link = (page or PAGE_PLACEHOLDER + "/index.html") + "?data=" + quote(folder_link, safe="")
    state.link = SimpleNamespace(mode=mode, page=page, folder_link=folder_link, share_link=share_link, read_pass=read_pass)

    lines = []
    # a folder without manifest.json is an export that stopped part-way: catalog.py leaves it out of datasets.json too
    for ds in sorted(p for p in export.iterdir() if p.is_dir() and (p / "manifest.json").is_file()):
        man = json.loads((ds / "manifest.json").read_text(encoding="utf-8"))
        d = next((x for x in state.datasets or [] if x["id"] == ds.name), None)
        o = state.options.get(ds.name, {})
        group = o.get("group") or (GROUP.get(d["kind"], "transcriptomics") if d else "transcriptomics")
        label = o.get("label") or man.get("label") or ds.name
        entry = f'{ds.name}: {{ group: "{group}", label: {json.dumps(label)}'
        entry += f', base: "{base_url}/{ds.name}/" }},' if mode == "public" else " },"
        lines.append("    " + entry)

    print(f"This data is {mode.upper()}.\n")
    print("# 1. Upload export/ to", base_url)
    if on_azure:
        print("# Make a short-lived WRITE pass (SAS) for the container first")
        print("# (Azure portal or Storage Explorer: container > Generate SAS, permissions Write + Create + List).")
        print('$sas = "<paste the SAS token here, starting with sv=>"')
        print(f'azcopy copy "{export}\\*" "{base_url}?$sas" --recursive')
    else:
        print(f"# Copy everything inside {export} to that folder with your host's upload tool.")
    if mode == "public":
        print("# Public data: allow anonymous read access to the folder.")
    print("# 2. Then run the next cell: it tests that the SODA page may read the folder (CORS).")

    if mode == "public" or read_pass:
        print("\n# 3. THE LINK TO SHARE:")
        print(share_link)
        if not page:
            print("# (SODA_PAGE is not set: replace <the SODA page's address> with the page's address,")
            print("#  or set SODA_PAGE in the advanced settings and run that cell and this one again.)")
        if read_pass:
            print("# This link carries a read pass: anyone holding it can read the data until the pass")
            print("# expires. Send it only to people who may see the data, like a password.")
    else:
        print("\n# 3. No link to share: private data opens on your password-protected SODA site")
        print("#    (deploy/azure/README.md), with the config.js lines below. For a link instead,")
        print("#    set READ_PASS in the advanced settings to a READ-only pass and run this cell again.")

    print("\n// Only for your own SODA site: paste into web/SpatialViewer/config.js (copy from config.example.js first)")
    print("datasets: {\n" + "\n".join(lines) + "\n},")
    print('tokenEndpoint: "/api/token",' if mode == "private" else 'tokenEndpoint: "",')


def storage_kind(link):
    host = (urlsplit(link).hostname or "").lower()
    if host.endswith(".blob.core.windows.net"):
        return "azure"
    if host.endswith(".amazonaws.com") or host.endswith(".cloudfront.net"):
        return "s3"
    if host in ("storage.googleapis.com", "storage.cloud.google.com") or host.endswith(".storage.googleapis.com"):
        return "gcs"
    if host in ("localhost", "127.0.0.1", "::1"):
        return "local"
    return "web"


def cors_fix(link, origin):
    """The setting to change, per storage provider (as in docs/SHARE_YOUR_DATA.md, Step 4)."""
    o = origin or "<the SODA page's address, e.g. https://soda.mgifive.org>"
    kind = storage_kind(link)
    if kind == "azure":
        return (f"Azure: Storage account -> Resource sharing (CORS) -> Blob service: Allowed origins {o},\n"
                "  Allowed methods GET, HEAD, OPTIONS, Allowed headers *, Exposed headers *, Max age 3600. Save.\n"
                "  (If you list headers instead of *: allow Range, expose Content-Range and Content-Length.)")
    if kind == "s3":
        return ("Amazon S3: Bucket -> Permissions -> Cross-origin resource sharing (CORS) -> Edit, paste:\n"
                f'  [{{ "AllowedOrigins": ["{o}"], "AllowedMethods": ["GET", "HEAD"],\n'
                '     "AllowedHeaders": ["*"], "ExposeHeaders": ["Content-Range", "Content-Length"], "MaxAgeSeconds": 3600 }]'
                + ("\n  (CloudFront: also forward the Origin header and the CORS response headers.)"
                   if (urlsplit(link).hostname or "").endswith(".cloudfront.net") else ""))
    if kind == "gcs":
        u = urlsplit(link)
        host = u.hostname or ""
        bucket = host[:-len(".storage.googleapis.com")] if host.endswith(".storage.googleapis.com") \
            else (u.path.strip("/").split("/") or ["<bucket>"])[0] or "<bucket>"
        return ("Google Cloud Storage: save this as cors.json\n"
                f'  [{{"origin": ["{o}"], "method": ["GET", "HEAD"], "responseHeader": ["Content-Type", "Range", "Content-Range", "Content-Length"], "maxAgeSeconds": 3600}}]\n'
                f"  and run: gcloud storage buckets update gs://{bucket} --cors-file=cors.json")
    if kind == "local":
        return ("A local test server: serve the folder with web/SpatialViewer/static_serve.py "
                "(it allows every page); python -m http.server does not.")
    return (f"A web server: send the header 'Access-Control-Allow-Origin: {o}' with every file of the folder,\n"
            "  'Access-Control-Expose-Headers: Content-Range, Content-Length', allow the Range header in a preflight,\n"
            "  and answer Range requests (the nginx example is in docs/HOSTING.md, 'Host notes').")


def range_test(root, qs, cat, origin, page, timeout=30):
    # The genes of an export are packed (a few large files per layer) and the page reads one gene
    # with an HTTP Range request. Ask for the first 12 bytes of one pack the way the page would:
    # the answer must be 206 with exactly those bytes, readable from the page's origin.
    def get(url, headers, method="GET"):
        req = urllib.request.Request(url, method=method, headers={"Origin": origin, "User-Agent": "SODA-notebook-link-test", **headers})
        try:
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.status, r.headers, r.read(1 << 16)
        except urllib.error.HTTPError as e:
            return e.code, e.headers, b""
    q = "?" + qs if qs else ""
    for d in cat.get("datasets", []):
        st, _, body = get(f"{root}/{d['id']}/manifest.json{q}", {})
        if st != 200:
            continue
        try:
            m = json.loads(body.decode("utf-8-sig"))
        except ValueError:
            continue
        packs = ([("", m["genePack"])] if m.get("genePack") else []) + \
                [(f"expr/{L['key']}/", L["genePack"]) for L in m.get("expr", []) if L.get("genePack")] + \
                ([(m["molecules"].get("dir", "molecules") + "/", m["molecules"]["genePack"])]
                 if (m.get("molecules") or {}).get("genePack") else [])
        if not packs:
            continue
        sub, pk = packs[0]
        path = f"{d['id']}/{sub}{pk.get('dir', 'genepack')}/{pk['shards'][0]['file']}"
        shown = f"{root}/{path}" + ("?<read pass>" if qs else "")
        st, h, body = get(f"{root}/{path}{q}", {"Range": "bytes=0-11"})
        if st == 200:
            return [f"{shown} came back whole (HTTP 200): this host ignores Range requests, so the page cannot read "
                    "one gene out of a pack. Use blob storage / S3 / a web server that answers Range, or export again "
                    "with --gene-files loose for this host."]
        if st != 206 or len(body) != 12:
            return [f"A Range request for {shown} answered HTTP {st} with {len(body)} bytes (want 206 and 12 bytes)."]
        expose = (h.get("Access-Control-Expose-Headers") or "").lower()
        if expose != "*" and "content-range" not in expose:
            return ["Range requests work, but the page cannot see Content-Range (CORS). Fix:\n  " + cors_fix(root, origin if page else None)]
        st, h, _ = get(f"{root}/{path}{q}", {"Access-Control-Request-Method": "GET", "Access-Control-Request-Headers": "range"}, "OPTIONS")
        allow = (h.get("Access-Control-Allow-Headers") or "").lower()
        if st >= 400 or (allow != "*" and "range" not in allow):
            return [f"Browsers that ask before sending a Range header would be refused (CORS preflight: HTTP {st}, "
                    f"allowed headers {allow or 'none'}). Fix:\n  " + cors_fix(root, origin if page else None)]
        return []
    return []                                  # nothing packed (exports with --gene-files loose): nothing to test


def cors_test(link, origin, page, timeout=30):
    """Fetch <folder>/datasets.json as a browser on `origin` would. Returns (ok, problems)."""
    root, _, qs = link.partition("?")
    url = root.rstrip("/") + "/datasets.json" + ("?" + qs if qs else "")
    shown = root.rstrip("/") + "/datasets.json" + ("?<read pass>" if qs else "")   # never print the pass
    req = urllib.request.Request(url, headers={"Origin": origin, "User-Agent": "SODA-notebook-link-test"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            code, headers, body = r.status, r.headers, r.read()
    except urllib.error.HTTPError as e:
        code, headers, body = e.code, e.headers, b""
    except (urllib.error.URLError, OSError, ValueError) as e:
        return False, [f"Could not reach {shown}: {getattr(e, 'reason', e)}. Check CLOUD_URL."]
    problems = []
    if code == 404:
        problems.append(f"{shown} was not found (HTTP 404). Upload the whole export/ folder, with datasets.json, "
                        "to exactly this folder. (A private Azure container also answers 404 without a read pass.)")
    elif code in (401, 403, 409):
        problems.append(f"The storage refused {shown} (HTTP {code}): the folder is private. Public data: allow "
                        "anonymous read access. Private data: set READ_PASS to a read-only pass (or a new one if it expired).")
    elif code != 200:
        problems.append(f"{shown} answered HTTP {code}.")
    else:
        try:
            cat = json.loads(body.decode("utf-8-sig"))
            if not (isinstance(cat, dict) and isinstance(cat.get("datasets"), list) and cat["datasets"]):
                problems.append(f"{shown} lists no datasets. Run Section 4 to the end and upload datasets.json again.")
            elif (headers.get("Access-Control-Allow-Origin") or "").strip() in ("*", origin):
                problems += range_test(root.rstrip("/"), qs, cat, origin, page, timeout)
        except ValueError:
            problems.append(f"{shown} is not a SODA dataset list (is the link the export folder itself?).")
        except (urllib.error.URLError, OSError) as e:   # the Range test's requests: a timeout or a dropped connection
            problems.append(f"The Range test of one gene pack got no answer ({getattr(e, 'reason', e)}). "
                            "Check the connection and run this cell again.")
    allow = (headers.get("Access-Control-Allow-Origin") or "").strip()
    if allow not in ("*", origin):
        problems.append(f"The storage does not let the page on {origin} read it (CORS: "
                        + (f"it allows only {allow}" if allow else "no Access-Control-Allow-Origin header") + "). Fix:\n  "
                        + cors_fix(link, origin if page else None))
    return not problems, problems


def test_link():
    """Section 7, after the upload: read datasets.json and one gene pack the way the SODA page would, and
    print OK and THE LINK TO SHARE, or the exact setting to change."""
    s = _need(state.link, "Run the cell above (Section 7, is this data public?) first.")
    if s.mode == "private" and not s.read_pass:
        print("Private data without a read pass: only your password-protected site can read the folder,")
        print("so this test cannot. Set the CORS rule for your site's address as in deploy/azure/README.md.")
        return
    origin = "{0.scheme}://{0.netloc}".format(urlsplit(s.page)) if s.page else "https://example.org"
    if not s.page:
        print("SODA_PAGE is not set, so this tests as a page on https://example.org.")
        print("Set SODA_PAGE in the advanced settings to test with the real page's address.\n")
    ok, problems = cors_test(s.folder_link, origin, s.page)
    if ok:
        print(f"OK: a page on {origin} can read your folder.\n")
        print("THE LINK TO SHARE:")
        print(s.share_link)
        if not s.page:
            print("(replace <the SODA page's address> with the page's address first)")
        print("\nOpen it yourself first (a private browser window is best): you should see every dataset.")
    else:
        print("NOT READY yet:\n")
        for p in problems:
            print("- " + p + "\n")
        print("Change that, wait a minute, and run this cell again.")
