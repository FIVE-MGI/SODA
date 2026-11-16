"""Export a prepared dataset as a plain static file tree (no server needed).

    <dest>/<id>/manifest.json
    <dest>/<id>/tiles/<c>/<z>/<x>/<y>.png       one channel, 8-bit grey (or .jpg)
    <dest>/<id>/expr/<key>/genes.json           [{n, t}]
    <dest>/<id>/expr/<key>/positions.bin        float32 [x,y]*
    <dest>/<id>/expr/<key>/gene/<NAME>.bin      float32 [x,y,v]*   (NAME percent-encoded like encodeURIComponent)
    <dest>/<id>/expr/<key>/hvg/...              bin-major top-variable genes in spatial chunks (hvgchunks.py)
    <dest>/<id>/expr/<key>/binqc/...            per-bin total counts, genes detected, % MT (binqc.py)
    <dest>/<id>/expr/<key>/umap/...             whole-section UMAP + Leiden clusters (sectionumap.py)
    <dest>/<id>/expr/<key>/umap/celltype*       suggested cell type per bin / cluster from marker lists (celltypes.py)
    <dest>/<id>/expr/<key>/umap/clusters_user*  the lab's own clusterings / UMAP, matched by barcode (userclusters.py)
    <dest>/<id>/expr/grid<NNN>um/...            image datasets (Akoya / multiplex TIFF): marker means per grid
                                                square as a "protein" expr layer, with qc.json and umap/
                                                (protgrid.py)
    <dest>/<id>/molecules/index.json, gene/<NAME>.bin   imaging platforms (Xenium): every transcript,
                                                uint16 x, y per molecule, one file per gene (molecules.py)
    <dest>/<id>/expr/<key>/shapes/<set>/...     cell (and nucleus) outlines in spatial tiles (cellshapes.py)
    <dest>/<id>/expr/rois/...                   GeoMx ROI segments (geomx.py): gene files, binqc/ with the GeoMx
                                                QC values, hvg/ with every target, umap/ with the annotation
                                                columns as groupings and no embedding
    <dest>/<id>/overlays/<file>.geojson

Upload the <id> folder to any static host (Azure blob, S3, a web server) and
point config.js (or a `<id>.link.txt` in web/SpatialViewer) at its URL. Re-running skips files that already exist, so an
interrupted export can be resumed.
"""
from __future__ import annotations

import json
import math
import shutil
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import quote

import numpy as np

import binqc
import cellshapes
import celltypes as celltypes_mod
import hvgchunks
import molecules
import protgrid
import regions
import sectionumap
import userclusters
from expr import ExprReader
from pyramid import PyramidReader


def _encode(name: str) -> str:
    # identical to JavaScript encodeURIComponent
    return quote(name, safe="-_.!~*'()")


def verify_tiles(r: PyramidReader, out: Path, ext: str, log=print) -> list[int]:
    """Re-derive one tile per channel and compare it with what is on disk.

    One real dataset shipped with twelve of its thirty-two channels wrong - eight
    written as solid black, four saturated to white - while the pyramid they came
    from was perfectly intact. Nothing noticed until someone opened the dataset and
    saw markers missing, days later. The tile loop has no way to tell a blank tile
    from a blank marker, so the only honest check is to recompute a tile and look.

    One tile per channel at the coarsest level: that is a single small read each,
    it covers the whole slide, and it catches any channel that came out uniformly
    wrong - which is the shape every failure of this kind has taken.
    """
    from PIL import Image

    bad = []
    z = 0                                            # coarsest level is z 0 on disk
    for c in range(len(r.channels)):
        want = r.tile_u8(c, r.nlevels - 1, 0, 0)
        p = out / "tiles" / str(c) / str(z) / "0" / f"0.{ext}"
        if want is None:
            continue
        if not p.exists():
            bad.append(c); continue
        got = np.asarray(Image.open(p).convert("L"), dtype=np.float32)
        w = float(want.mean())
        # jpeg moves the mean of a near-black tile by a grey level or so; a channel
        # that is blank-when-it-should-not-be, or saturated, is off by far more
        if abs(w - float(got.mean())) > max(1.5, 0.15 * w):
            bad.append(c)
    if bad:
        names = ", ".join(f"{c} ({r.channels[c]['name']})" for c in bad)
        log(f"[export] *** {len(bad)} channel(s) do not match the pyramid: {names}")
        log(f"[export] *** delete {out / 'tiles'} for those channels and re-run - "
            f"an export skips tiles that already exist, so re-running alone will not fix it")
    else:
        log(f"[export] verified {len(r.channels)} channel(s) against the pyramid")
    return bad


def load_region_overlay(ds_dir: Path, m: dict, outline_names=regions.DEFAULT_OUTLINE_NAMES, log=print):
    """The one region set used for Tier A linking: the first overlay with polygons.

    Returns (file, geojson, feats, feat_id) or (None, ...). Multiple region sets per
    dataset is a later extension; one covers the Visium / Xenium case.
    """
    for ov in m.get("overlays", []):
        src = ds_dir / "overlays" / ov["file"]
        if not src.exists():
            continue
        gj = json.load(open(src, encoding="utf-8"))
        feats = regions.load_features(gj)
        feat_id = regions.assign_region_ids(feats, outline_names=outline_names)
        if any(rid for rid in feat_id.values()):
            n_out = sum(1 for rid in feat_id.values() if rid == 0)
            log(f"[export] region overlay '{ov['file']}': {max(feat_id.values())} regions"
                + (f", {n_out} outline feature(s) excluded" if n_out else ""))
            return ov["file"], gj, feats, feat_id
    return None, None, None, None


def export_static(ds_dir: Path, dest: Path, fmt: str = "png", workers: int = 8, log=print,
                  verify: bool = True, outline_names=regions.DEFAULT_OUTLINE_NAMES, hvg: bool = True,
                  umap: bool = True, celltypes: bool = True, markers: Path | None = None,
                  grid_um: float = protgrid.DEFAULT_GRID_UM,
                  tissue_smooth_um: float = protgrid.DEFAULT_TISSUE_SMOOTH_UM,
                  umap_sample: int | None = None, seam_correct: bool = False,
                  user_clusters: list[dict] | None = None) -> Path:
    """`markers`: the marker list for suggested cluster names. Default by modality:
    markers/rna_generic.json for RNA (Visium) layers, markers/protein_generic.json for the
    protein grid layer of an image dataset. `grid_um` (image datasets only): grid square size
    in µm for the protein layer; 0 = no grid layer. `tissue_smooth_um`: moving-average window
    (µm) over the square means before the grid's tissue threshold; 0 = raw means.
    `umap_sample`: the grid layer's UMAP is fitted on a sample of this many squares above it
    (None = sectionumap.UMAP_SAMPLE, 0 = all). `seam_correct`: correct the grid layer's square
    means for per-tile differences where tile seams are found (off: detected and reported only).
    `user_clusters`: the lab's own UMAP / clusterings per layer (userclusters.make_spec)."""
    ds_dir = Path(ds_dir)
    m = json.load(open(ds_dir / "manifest.json", encoding="utf-8"))
    userclusters.check_layers(user_clusters, m)
    out = Path(dest) / m["id"]
    out.mkdir(parents=True, exist_ok=True)
    r = PyramidReader(ds_dir, m)
    T = r.tile
    ext = "jpg" if fmt == "jpg" else "png"

    # ---- tiles -------------------------------------------------------------
    jobs = []
    for L, (H, W) in enumerate(r.shapes):
        z = r.nlevels - 1 - L
        nx, ny = math.ceil(W / T), math.ceil(H / T)
        for c in range(len(r.channels)):
            for ty in range(ny):
                for tx in range(nx):
                    p = out / "tiles" / str(c) / str(z) / str(tx) / f"{ty}.{ext}"
                    if not p.exists():
                        jobs.append((c, L, tx, ty, p))
    log(f"[export] {m['id']}: {len(jobs):,} tiles to write ({len(r.channels)} ch x {r.nlevels} levels) -> {out}")
    t0 = time.time()
    done = [0]

    def work(job):
        c, L, tx, ty, p = job
        data = r.tile_bytes(c, L, tx, ty, fmt=ext)
        if data is None:
            return
        p.parent.mkdir(parents=True, exist_ok=True)
        tmp = p.with_suffix(p.suffix + ".part")
        tmp.write_bytes(data)
        tmp.replace(p)
        done[0] += 1
        if done[0] % 2000 == 0:
            log(f"[export]   {done[0]:,}/{len(jobs):,} tiles  {time.time() - t0:.0f}s")

    with ThreadPoolExecutor(max_workers=workers) as ex:
        list(ex.map(work, jobs))
    log(f"[export] tiles done in {time.time() - t0:.0f}s")
    if verify:
        verify_tiles(r, out, ext, log=log)

    # ---- region overlay (Tier A): one region set, assigned once per dataset ----
    region_file, region_geojson, region_feats, region_feat_id = load_region_overlay(
        ds_dir, m, outline_names=outline_names, log=log)
    layer_regions: dict[str, list] = {}
    layer_hvg: dict[str, dict] = {}
    layer_umap: dict[str, dict] = {}
    layer_binqc: dict[str, dict] = {}

    # ---- expression --------------------------------------------------------
    for layer in m.get("expr", []):
        key = layer["key"]
        er = ExprReader(ds_dir / "expr" / f"{key}.zarr")
        d = out / "expr" / key
        (d / "gene").mkdir(parents=True, exist_ok=True)
        (d / "genes.json").write_text(json.dumps(er.gene_list(), separators=(",", ":")), encoding="utf-8")
        pos = er.positions()                         # same bins, same order, as region.bin below
        (d / "positions.bin").write_bytes(np.ascontiguousarray(pos, dtype="<f4").tobytes())
        if region_feats is not None:
            ids, table = regions.index_bins(pos, region_feats, region_feat_id)
            if table:
                (d / "region.bin").write_bytes(np.ascontiguousarray(ids, dtype="<u2").tobytes())
                layer_regions[key] = table
                log(f"[export] expr {key}: {len(table)} regions, "
                    f"{int(np.count_nonzero(ids)):,}/{len(ids):,} bins assigned")
        # GeoMx ROI segments (geomx.py): a few dozen rows. Their own QC values, every target in hvg/,
        # the annotation columns as groupings; no UMAP, Leiden, cell types or variable-gene selection.
        roi = layer.get("kind") == "rois"
        if roi:
            import geomx
            if any(s["layer"] == key for s in user_clusters or []):
                raise SystemExit(f"[geomx] --user-layer {key}: an ROI layer takes its groups from the annotation sheet. "
                                 f"Add your columns there (or give a second --annotation file) and run prepare.py geomx again.")
            rq, rh, ru = geomx.export_layer(er, d, ds_dir, layer, log=log)
            for store, val in ((layer_binqc, rq), (layer_hvg, rh), (layer_umap, ru)):
                if val:
                    store[key] = val
        # per-bin QC (total counts, genes detected, % MT) for the QC filter; same 8 um floor
        if not roi and (not layer.get("binSizeUm") or layer["binSizeUm"] >= 8):
            q = binqc.write_bin_qc(er, d, log=log)
            if q:
                layer_binqc[key] = q
        # circle -> UMAP: bin-major top-variable genes. Not for the 2 um grid: too sparse to embed.
        if hvg and not roi and (not layer.get("binSizeUm") or layer["binSizeUm"] >= 8):
            h = hvgchunks.write_hvg_chunks(er, d, layer.get("micronsPerPixel"), log=log)
            if h:
                layer_hvg[key] = h
        # whole-section UMAP + Leiden clusters of the same bins (same 8 um floor)
        u = None
        if umap and not roi and (not layer.get("binSizeUm") or layer["binSizeUm"] >= 8):
            u = sectionumap.write_section_umap(er, d, log=log)
        # the lab's own UMAP / clusters for this layer (or SODA's put back after an earlier import)
        u = userclusters.for_layer(user_clusters, er, d, ds_dir, u, log=log) if not roi else None
        if u:
            layer_umap[key] = u
            # suggested cluster names from the marker lists, on the same embedded bins
            if celltypes:
                c = celltypes_mod.write_celltypes(er, d, markers=markers, log=log)
                if c:
                    layer_umap[key] = {**u, "celltypes": c}
        # slice in memory: reading each gene through zarr decompresses whole chunks
        indices = np.asarray(er.g["indices"])
        data = np.asarray(er.g["data"], dtype=np.float32)
        indptr = er.indptr
        xy = er.xy
        n = 0
        for gi, name in enumerate(er.genes):
            p = d / "gene" / f"{_encode(name)}.bin"
            if p.exists():
                continue
            a, b = int(indptr[gi]), int(indptr[gi + 1])
            if b <= a:
                continue
            ids = indices[a:b]
            xyv = np.empty((b - a, 3), np.float32)
            xyv[:, :2] = xy[ids]
            xyv[:, 2] = data[a:b]
            xyv = xyv[~np.isnan(xyv[:, 0])]
            if len(xyv) == 0:
                continue
            p.write_bytes(np.ascontiguousarray(xyv, dtype="<f4").tobytes())
            n += 1
        log(f"[export] expr {key}: {n:,} gene files")

    # ---- molecules (Xenium ...): per-gene point files, already in their static format ----
    if m.get("molecules"):
        molecules.export_molecules(ds_dir, out, log=log)

    # ---- cell outlines (cellshapes.py): tiles of polygons, already in their static format ----
    layer_shapes = cellshapes.export_shapes(ds_dir, out, m, log=log)

    # ---- protein grid layer (image datasets): marker means per square ----------
    grid_layer = None
    if grid_um and grid_um > 0 and protgrid.applies_to(m):
        grid_layer = protgrid.write_grid_layer(r, m, out, grid_um=grid_um, markers=markers,
                                               umap=umap, celltypes=celltypes, log=log,
                                               tissue_smooth_um=tissue_smooth_um, umap_sample=umap_sample,
                                               seam_correct=seam_correct)
        if grid_layer and region_feats is not None:
            gk = grid_layer["key"]
            pos = np.fromfile(out / "expr" / gk / "positions.bin", dtype="<f4").reshape(-1, 2)
            ids, table = regions.index_bins(pos, region_feats, region_feat_id)
            if table:
                (out / "expr" / gk / "region.bin").write_bytes(np.ascontiguousarray(ids, dtype="<u2").tobytes())
                grid_layer = {**grid_layer, "regions": table, "regionOverlay": region_file}
                log(f"[export] expr {gk}: {len(table)} regions, "
                    f"{int(np.count_nonzero(ids)):,}/{len(ids):,} squares assigned")

    # ---- overlays + manifest -----------------------------------------------
    for ov in m.get("overlays", []):
        src = ds_dir / "overlays" / ov["file"]
        if not src.exists():
            continue
        (out / "overlays").mkdir(exist_ok=True)
        if ov["file"] == region_file:
            # write the region ids back onto the features so a click reads one straight off
            regions.annotate_geojson(region_geojson, region_feats, region_feat_id)
            (out / "overlays" / ov["file"]).write_text(
                json.dumps(region_geojson, separators=(",", ":")), encoding="utf-8")
        else:
            shutil.copyfile(src, out / "overlays" / ov["file"])
    sm = dict(m)
    sm["static"] = True
    sm["tileFormat"] = ext
    sm.pop("source", None)          # local paths are nobody else's business
    if layer_regions:               # Tier A: per-layer region tables for the viewer
        sm["expr"] = [({**layer, "regions": layer_regions[layer["key"]], "regionOverlay": region_file}
                       if layer["key"] in layer_regions else dict(layer))
                      for layer in m.get("expr", [])]
    if layer_hvg:
        sm["expr"] = [({**layer, "hvg": layer_hvg[layer["key"]]} if layer["key"] in layer_hvg else layer)
                      for layer in sm.get("expr", [])]
    if layer_binqc:
        sm["expr"] = [({**layer, "binQc": layer_binqc[layer["key"]]} if layer["key"] in layer_binqc else layer)
                      for layer in sm.get("expr", [])]
    if layer_umap:
        sm["expr"] = [({**layer, "umap": layer_umap[layer["key"]]} if layer["key"] in layer_umap else layer)
                      for layer in sm.get("expr", [])]
    sm["expr"] = cellshapes.manifest_layers(sm.get("expr", []), layer_shapes)
    if grid_layer:
        sm["expr"] = [L for L in sm.get("expr", []) if L.get("key") != grid_layer["key"]] + [grid_layer]
    (out / "manifest.json").write_text(json.dumps(sm, indent=1), encoding="utf-8")
    log(f"[export] done -> {out}")
    return out
