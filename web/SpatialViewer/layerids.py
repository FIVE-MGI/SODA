"""Optional per-layer id list: the instrument's own barcode / cell id of every bin, for the viewer's
"Save barcodes (CSV)" (selections.js).

    <dest>/<id>/expr/<key>/ids.txt.gz    gzip of UTF-8 text, one id per line, one line per row of
                                         positions.bin, in the same order (no header)

and the layer's manifest entry gets `"ids": {"file": "ids.txt.gz", "n": rows, "bytes": size of the file,
"column": "barcode" | "cell_id"}`. Without that entry the export has no id list (the default).

Off unless asked: `prepare.py export --ids`, or `prepare.py ids --id X --out <cache> --dest <export root>` to
add the lists to an export made before (it reads the barcodes from the prepare cache, writes the files and
updates manifest.json; nothing else changes). The notebook, the README examples and catalog.py never ask.

The ids go into the folder that gets published, so two gates come before anything is written, for every
layer of the dataset (a refusal writes nothing at all):

1. Allow-list. Every id must FULLY match the id format of its platform (FORMATS below, taken from what the
   readers in this folder put into the cache). Anything else - a merge prefix such as "S1_" or a sample name
   glued to a barcode, a run or slide name in a cell id, a GeoMx Sample_ID - stops it. Ids are never
   rewritten or stripped here: fix the input, or leave --ids off.
2. deid.find_ids over the whole list: the lab's own sample-ID patterns.

A refusal names the layer, the platform and how many ids failed; it never prints an id (the log of a
refused run would otherwise carry the very IDs it refused). A platform with no format here gets no list.

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


def check_ids(ids: list[str], platform: str | None, key: str) -> None:
    """Refuse (SystemExit) unless every id fits the platform's format and none looks like a sample ID.
    The message names the layer and counts; it never contains an id."""
    where = f"{key}/{FILE}"
    if platform not in _FULL:
        raise SystemExit(f"[ids] REFUSING to write {where}: no known id format for this platform "
                         f"({platform or 'unknown'}); id lists are written only for {', '.join(sorted(FORMATS))}")
    rx = _FULL[platform]
    bad = sum(1 for i in ids if not rx.match(i))
    if bad:
        raise SystemExit(f"[ids] REFUSING to write {where}: {bad:,} of {len(ids):,} ids do not match the {platform} "
                         f"id format ({FORMATS[platform]}); a sample or run name in the ids, or a merge prefix, "
                         f"would be published. Nothing was written. Leave --ids off, or prepare from the "
                         f"instrument's own files.")
    hits = len(deid.find_ids("\n".join(ids)))          # one id per line: a match never spans two ids
    if hits:
        raise SystemExit(f"[ids] REFUSING to write {where}: {hits:,} ids look like sample IDs (deid.py patterns). "
                         f"Nothing was written.")


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


def plan(ds_dir: Path, m: dict, log=print) -> dict:
    """Read and check the ids of every layer that gets a list, before anything is written: {key: ids}.
    One refused layer stops the whole dataset (SystemExit), so no partial set of lists is published."""
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
        check_ids(ids, platform_of(m.get("kind"), layer), key)
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
