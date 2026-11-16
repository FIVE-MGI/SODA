"""Optional per-layer id list: the instrument's own barcode / cell id of every bin, for the viewer's
"Save barcodes (CSV)" (selections.js).

    <dest>/<id>/expr/<key>/ids.txt.gz    gzip of UTF-8 text, one id per line, one line per row of
                                         positions.bin, in the same order (no header)

and the layer's manifest entry gets `"ids": {"file": "ids.txt.gz", "n": rows, "bytes": size of the file,
"column": "barcode" | "cell_id"}`. Without that entry the export has no id list (the default).

Off unless asked: `prepare.py export --ids`, or `prepare.py ids --id X --out <cache> --dest <export root>` to
add the lists to an export made before (it reads the barcodes from the prepare cache, writes the files and
updates manifest.json; nothing else changes).

The ids are what the cache keeps per layer (expr.build_expr / build_cell_layer: Visium spot and HD bin
barcodes, Xenium cell_id, CosMx <cell_ID>_<fov>, MERSCOPE EntityID, Stereo-seq bin x_y / cell id). They are
not sample IDs; every list is still run through deid.find_ids first and the export stops if one is found.
GeoMx ROI layers (their row names carry slide and sample names) and protein grid layers (squares have no
id) get none.

Size: gzip of the text. Measured on the smoke-test demo data, see docs/DATA_FORMAT.md.
"""
from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np

import deid

FILE = "ids.txt.gz"
SKIP_KINDS = {"rois", "protein"}


def column_for(layer: dict) -> str:
    return "cell_id" if layer.get("kind") == "cells" else "barcode"


def layer_ids(er, ds_dir: Path | None = None, log=print) -> list[str] | None:
    """The ids of the positions.bin rows of a layer (rows without a position left out), or None."""
    import userclusters
    try:
        bcs = userclusters.layer_barcodes(er, ds_dir, log=log)
    except SystemExit as e:
        log(f"[ids] {er.meta.get('key')}: no barcodes in the cache, no id list ({e})")
        return None
    on = ~np.isnan(np.asarray(er.xy)[:, 0])
    return [b for b, keep in zip(bcs, on.tolist()) if keep]


def write_layer_ids(er, d: Path, layer: dict, ds_dir: Path | None = None, log=print) -> dict | None:
    """Write expr/<key>/ids.txt.gz for one layer; returns its manifest entry (or None: no ids for it)."""
    key = layer.get("key")
    if layer.get("kind") in SKIP_KINDS:
        return None
    ids = layer_ids(er, ds_dir, log=log)
    if ids is None:
        return None
    n = int(np.count_nonzero(~np.isnan(np.asarray(er.xy)[:, 0])))
    if len(ids) != n:
        raise SystemExit(f"[ids] {key}: {len(ids)} ids for {n} rows of positions.bin")
    bad = [i for i in ids if "\n" in i or "\r" in i]
    if bad:
        raise SystemExit(f"[ids] {key}: ids with a line break, e.g. {bad[0]!r}")
    txt = "\n".join(ids) + "\n"
    hits = deid.find_ids(txt)
    if hits:
        raise SystemExit(f"[ids] REFUSING to write {key}/{FILE}: it would contain sample IDs: {sorted(set(hits))[:5]}")
    data = gzip.compress(txt.encode("utf-8"), compresslevel=9, mtime=0)
    d = Path(d)
    d.mkdir(parents=True, exist_ok=True)
    tmp = d / (FILE + ".part")
    tmp.write_bytes(data)
    tmp.replace(d / FILE)
    log(f"[ids] {key}: {n:,} {column_for(layer)}s, {len(data) / 1e3:.1f} kB ({len(data) / max(n, 1):.2f} bytes per row; "
        f"{len(txt.encode('utf-8')) / max(n, 1):.1f} as plain text)")
    return {"file": FILE, "n": n, "bytes": len(data), "column": column_for(layer)}


def add_to_export(ds_dir: Path, dest: Path, log=print) -> dict:
    """`prepare.py ids`: write the id lists of every layer into an existing export and name them in its manifest."""
    from expr import ExprReader
    ds_dir = Path(ds_dir)
    m = json.load(open(ds_dir / "manifest.json", encoding="utf-8"))
    out = Path(dest) / m["id"]
    mp = out / "manifest.json"
    if not mp.exists():
        raise SystemExit(f"[ids] no export at {out}: run prepare.py export first")
    sm = json.load(open(mp, encoding="utf-8"))
    cache_keys = {L["key"] for L in m.get("expr", [])}
    done = {}
    for layer in sm.get("expr", []):
        key = layer.get("key")
        if key not in cache_keys or not (out / "expr" / key / "positions.bin").exists():
            continue
        er = ExprReader(ds_dir / "expr" / f"{key}.zarr")
        npos = (out / "expr" / key / "positions.bin").stat().st_size // 8
        entry = write_layer_ids(er, out / "expr" / key, layer, ds_dir, log=log)
        if entry is None:
            continue
        if entry["n"] != npos:
            (out / "expr" / key / FILE).unlink()
            raise SystemExit(f"[ids] {key}: the cache has {entry['n']} rows, the export's positions.bin {npos}: "
                             f"this export was made from another cache")
        layer["ids"] = entry
        done[key] = entry
    tmp = mp.with_suffix(".json.part")
    tmp.write_text(json.dumps(sm, indent=1), encoding="utf-8")
    tmp.replace(mp)
    log(f"[ids] {m['id']}: id lists for {', '.join(done) or 'no layer'} -> {out}")
    return done
