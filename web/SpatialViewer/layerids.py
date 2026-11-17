"""Per-layer id list: the instrument's own barcode / cell id of every bin, for the viewer's
"Save barcodes (CSV)" (selections.js), as Loupe Browser, Xenium Explorer and Vitessce keep them.

    <dest>/<id>/expr/<key>/ids.txt.gz    gzip of UTF-8 text, one id per line, one line per row of
                                         positions.bin, in the same order (no header)

and the layer's manifest entry gets `"ids": {"file": "ids.txt.gz", "n": rows, "bytes": size of the file,
"column": "barcode" | "cell_id"}`. Without that entry the export has no id list.

On by default: `prepare.py export` (and batch.py, the notebook and the desktop app, which run it) writes the
lists; `--no-ids` leaves them out. `prepare.py ids --id X --out <cache> --dest <export root>` adds them to an
export made without them (it reads the barcodes from the prepare cache, writes the files and updates
manifest.json; nothing else changes). `--ids` is still accepted: it asks for the lists explicitly, which makes
a failed check stop the export (below). The data owner decides what is published: the checks below catch the
usual mistakes, they do not make an id list safe to publish.

The ids go into the folder that gets published, so two gates come before anything is written, for every
layer of the dataset:

1. Format, for the platforms in FORMATS below (taken from what the readers in this folder put into the cache):
   every id must FULLY match its platform's format. Anything else - a merge prefix such as "S1_" or a sample
   name glued to a barcode, a run or slide name in a cell id - fails it. Ids are never rewritten or stripped
   here: fix the input, or use --no-ids. A platform with no format here (an AnnData or SpatialData table's
   obs names, ...) has no such check: its ids pass when they hold no line break or control character (one id
   per line must stay one row) and gate 2 finds nothing.
2. deid.find_ids over the whole list: the lab's own sample-ID patterns, for every platform.

What a failed gate does depends on how the lists were asked for:
  - by default (no flag): that layer's list is left out with a *** message, and the export goes on (other
    layers keep theirs);
  - explicitly (`--ids`, or `prepare.py ids`): the whole dataset is refused and nothing is written at all.
A message names the layer, the platform, the reason and how many ids failed; it never prints an id (the log
of a run would otherwise carry the very IDs it left out).

GeoMx ROI layers (their row names carry slide and sample names) and protein grid layers (squares have no
id) get none.

Size: gzip of the text. Measured on the smoke-test demo data, see docs/DATA_FORMAT.md.
"""
from __future__ import annotations

import gzip
import json
import re
from pathlib import Path

import numpy as np

import deid

FILE = "ids.txt.gz"
SKIP_KINDS = {"rois", "protein"}

# The id formats the readers write into the cache, per platform (each pattern must match the WHOLE id).
#   visium_hd  expr.build_expr from Space Ranger's filtered_feature_bc_matrix.h5: s_<bin µm, 3 digits>um_<row 5>_<col 5>-1
#   visium     the same for classic spots: 16 nucleotides and the GEM well suffix (-1; aggr gives -2, -3 ...)
#   xenium     xenium.read_matrix: cell_feature_matrix barcodes = cells.csv cell_id; since Xenium Ranger 1.3
#              8 letters a-p and a dataset suffix (aaabbbcc-1), integers before that
#   cosmx      cosmx.read_cells: the matrix's "cell" column when it has one (AtoMx: c_<slide>_<fov>_<cell>, all
#              integers; the slide number is the instrument's slide index, kept as it is because it is how
#              AtoMx and the lab's Seurat objects name the cell), else <cell_ID>_<fov> made from two integers
#   merscope   merscope.read_cells: EntityID, an integer (up to 39 digits in the demo and in Vizgen output)
#   stereoseq  stereoseq.py: bins "<x>_<y>" (DNB coordinates) or with --barcodes stereopy x * 2^32 + y; cells the
#              cell bin GEF's integer id (or x * 2^32 + y)
FORMATS = {
    "visium_hd": r"s_\d{3}um_\d{5}_\d{5}-1",
    "visium": r"[ACGT]{16}-\d{1,3}",
    "xenium": r"[a-p]{8}-\d{1,3}|\d{1,12}",
    "cosmx": r"c_\d{1,6}_\d{1,6}_\d{1,9}|\d{1,9}_\d{1,6}",
    "merscope": r"\d{1,40}",
    "stereoseq": r"\d{1,9}_\d{1,9}|\d{1,20}",
}
_FULL = {k: re.compile(rf"(?:{v})\Z") for k, v in FORMATS.items()}


def column_for(layer: dict) -> str:
    return "cell_id" if layer.get("kind") == "cells" else "barcode"


def platform_of(dataset_kind: str | None, layer: dict) -> str | None:
    """The FORMATS key for a layer of a dataset (cache manifest `kind`), or None when there is none."""
    kind = (dataset_kind or "").lower()
    plat = str(layer.get("platform") or "").lower()
    if plat in FORMATS:
        return plat
    if kind == "visium":
        return "visium_hd" if layer.get("binSizeUm") else "visium"
    if kind in FORMATS:
        return kind
    return None


def problem(ids: list[str], platform: str | None) -> tuple[str, str] | None:
    """Why this list may not be published, as (kind, reason), or None when the gates pass.
    kind: "mismatch" (ids outside the platform's format), "shape" (a line break or control character in an id of
    a platform without a format), "sample" (deid.py's sample-ID patterns). The reason names counts and the
    format; it never contains an id."""
    if platform in _FULL:
        rx = _FULL[platform]
        bad = sum(1 for i in ids if not rx.match(i))
        if bad:
            return ("mismatch", f"{bad:,} of {len(ids):,} ids do not match the {platform} id format ({FORMATS[platform]}); "
                                f"a sample or run name in the ids, or a merge prefix, would be published")
    else:
        bad = sum(1 for i in ids if _CTRL.search(i))
        if bad:
            return ("shape", f"{bad:,} of {len(ids):,} ids hold a line break or another control character, so the "
                             f"one-id-per-line list would not line up with the rows ({platform or 'no known id format'})")
    hits = len(deid.find_ids("\n".join(ids)))          # one id per line: a match never spans two ids
    if hits:
        return ("sample", f"{hits:,} ids look like sample IDs (deid.py patterns)")
    return None


_CTRL = re.compile(r"[\x00-\x1f\x7f]")


def check_ids(ids: list[str], platform: str | None, key: str) -> None:
    """Refuse (SystemExit) unless the list passes the gates: the lists were asked for explicitly (--ids,
    prepare.py ids). The message names the layer and counts; it never contains an id."""
    why = problem(ids, platform)
    if why is None:
        return
    raise SystemExit(f"[ids] REFUSING to write {key}/{FILE}: {why[1]}. Nothing was written. Without --ids the "
                     f"export leaves this list out and goes on; --no-ids writes no id list; or fix the ids at the source.")


def skip_message(key: str, platform: str | None, why: tuple[str, str]) -> str:
    """The line logged when a default export leaves a layer's list out (never an id): a *** line, which the
    notebook and batch.py flag for reading."""
    return (f"[ids] *** {key}: barcode / cell id list LEFT OUT of the export: {why[1]}. The export goes on without "
            f"it. Fix the ids at the source, or use --no-ids to leave every list out.")


def layer_ids(er, ds_dir: Path | None = None, log=print) -> list[str] | None:
    """The ids of the positions.bin rows of a layer (rows without a position left out), or None."""
    import userclusters
    try:
        bcs = userclusters.layer_barcodes(er, ds_dir, log=log)
    except SystemExit:
        log(f"[ids] {er.meta.get('key')}: no barcodes in the cache, no id list")
        return None
    on = ~np.isnan(np.asarray(er.xy)[:, 0])
    return [b for b, keep in zip(bcs, on.tolist()) if keep]


def plan(ds_dir: Path, m: dict, log=print, strict: bool = True) -> dict:
    """Read and check the ids of every layer that gets a list, before anything is written: {key: ids}.
    strict (--ids, prepare.py ids): one refused layer stops the whole dataset (SystemExit), so no partial set of
    lists is published. Not strict (the default export): a layer that fails is left out with a message
    (skip_message) and the others keep theirs."""
    from expr import ExprReader
    out = {}
    for layer in m.get("expr", []):
        key = layer.get("key")
        if layer.get("kind") in SKIP_KINDS:
            continue
        zp = Path(ds_dir) / "expr" / f"{key}.zarr"
        if not zp.exists():
            continue
        er = ExprReader(zp)
        ids = layer_ids(er, ds_dir, log=log)
        if ids is None:
            continue
        platform = platform_of(m.get("kind"), layer)
        if strict:
            check_ids(ids, platform, key)
        else:
            why = problem(ids, platform)
            if why is not None:
                log(skip_message(key, platform, why))
                continue
        out[key] = ids
    return out


def write_ids(ids: list[str], d: Path, layer: dict, log=print) -> dict:
    """Write expr/<key>/ids.txt.gz from ids that `plan` checked; returns the manifest entry."""
    key = layer.get("key")
    txt = "\n".join(ids) + "\n"
    data = gzip.compress(txt.encode("utf-8"), compresslevel=9, mtime=0)
    d = Path(d)
    d.mkdir(parents=True, exist_ok=True)
    tmp = d / (FILE + ".part")
    tmp.write_bytes(data)
    tmp.replace(d / FILE)
    n = len(ids)
    log(f"[ids] {key}: {n:,} {column_for(layer)}s, {len(data) / 1e3:.1f} kB ({len(data) / max(n, 1):.2f} bytes per row; "
        f"{len(txt.encode('utf-8')) / max(n, 1):.1f} as plain text)")
    return {"file": FILE, "n": n, "bytes": len(data), "column": column_for(layer)}


def add_to_export(ds_dir: Path, dest: Path, log=print) -> dict:
    """`prepare.py ids`: write the id lists of every layer into an existing export and name them in its manifest."""
    ds_dir = Path(ds_dir)
    m = json.load(open(ds_dir / "manifest.json", encoding="utf-8"))
    out = Path(dest) / m["id"]
    mp = out / "manifest.json"
    if not mp.exists():
        raise SystemExit(f"[ids] no export at {out}: run prepare.py export first")
    sm = json.load(open(mp, encoding="utf-8"))
    planned = plan(ds_dir, m, log=log)                   # every check first: a refusal writes nothing
    for key, ids in planned.items():
        p = out / "expr" / key / "positions.bin"
        if not p.exists():
            raise SystemExit(f"[ids] {key}: the export has no positions.bin for this layer")
        if p.stat().st_size // 8 != len(ids):
            raise SystemExit(f"[ids] {key}: the cache has {len(ids):,} rows, the export's positions.bin "
                             f"{p.stat().st_size // 8:,}: this export was made from another cache. Nothing was written.")
    done = {}
    for layer in sm.get("expr", []):
        key = layer.get("key")
        if key in planned:
            layer["ids"] = done[key] = write_ids(planned[key], out / "expr" / key, layer, log=log)
    tmp = mp.with_suffix(".json.part")
    tmp.write_text(json.dumps(sm, indent=1), encoding="utf-8")
    tmp.replace(mp)
    log(f"[ids] {m['id']}: id lists for {', '.join(done) or 'no layer'} -> {out}")
    return done
