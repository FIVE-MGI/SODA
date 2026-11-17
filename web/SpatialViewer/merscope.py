"""MERSCOPE (Vizgen) region folder -> mosaic stains + cells layer + molecules layer.

Reads, when present (each part is optional), from the folder given as --src:

    images/mosaic_<stain>_z<k>.tif             one 16-bit mosaic per stain and z-plane (DAPI, PolyT,
                                               Cellbound1-3, ...); one z-plane is imported
    images/micron_to_mosaic_pixel_transform.csv   3 x 3 matrix: mosaic pixel = A . (x µm, y µm, 1)
    cell_by_gene.csv[.gz]                      cell, one column per gene and Blank-* barcode (counts)
    cell_metadata.csv[.gz]                     cell id (EntityID or the unnamed first column), center_x,
                                               center_y, min_x, max_x, min_y, max_y (µm), volume
    cell_boundaries.parquet | *_micron_space.parquet | cell_boundaries/feature_data_<fov>.hdf5
                                               the cells' outlines in µm (cellshapes.py)
    detected_transcripts.csv[.gz] | .parquet   barcode_id, global_x, global_y (µm), global_z, x, y, fov,
                                               gene, transcript_id [, cell_id]

The two cell tables may also sit one folder down (a Vizgen Post-processing Tool run such as
`vpt-run/`); --cells-dir picks the folder when there are several. .parquet copies of the
tables are read too, and a re-saved table with extra row-number columns in front is fine:
columns are found by name.

MERSCOPE coordinates are µm over the region ("global"). The transform file maps them to mosaic
pixels (a scale of about 9.26 px/µm = 0.108 µm per pixel and a shift, y pointing down as in the
image), so cells and transcripts need nothing else. Without images (some deposits leave them
out because of their size) the dataset gets a neutral background instead: the transcript
density at 1 µm per pixel, in the same orientation.

Cells become an ordinary expression layer (expr.build_cell_layer), transcripts the molecules
layer (molecules.py). The Blank-* barcodes (MERSCOPE's negative controls) are left out of both.
"""
from __future__ import annotations

import json
import re
import time
from pathlib import Path

import numpy as np

import flatfiles as ff

PLATFORM = "merscope"
CONTROL_PREFIXES = ("blank",)
_MOSAIC = re.compile(r"^mosaic_(.+)_z(\d+)\.tiff?$", re.I)
TRANSFORM = "micron_to_mosaic_pixel_transform.csv"
# stains in the order people look for them; anything else follows by name
_STAIN_ORDER = ("dapi", "polyt")


def is_control(name: str) -> bool:
    return str(name).strip().lower().startswith(CONTROL_PREFIXES)


def find_inputs(folder: Path, cells_dir: Path | None = None) -> dict:
    """The MERSCOPE files under `folder` (None / {} for each part that is missing)."""
    folder = Path(folder)
    skip = ("images", "result_tiles", "cell_boundaries")
    cd = Path(cells_dir) if cells_dir else folder
    depth = 0 if cells_dir else 2
    inp = {"transcripts": ff.find_table(folder, r"detected_transcripts", depth=1, skip_dirs=skip),
           "matrix": ff.find_table(cd, r"cell_by_gene", depth=depth, skip_dirs=skip),
           "cells": ff.find_table(cd, r"cell_metadata", depth=depth, skip_dirs=skip)}
    if inp["matrix"] and inp["cells"] and inp["matrix"].parent != inp["cells"].parent:
        same = ff.find_table(inp["matrix"].parent, r"cell_metadata", depth=0)
        if same:
            inp["cells"] = same                       # both tables of one segmentation run
    mosaics: dict[int, dict[str, Path]] = {}
    idir = next((d for d in (folder / "images", folder) if d.is_dir() and any(_MOSAIC.match(p.name) for p in d.iterdir())), None)
    if idir:
        for p in sorted(idir.iterdir()):
            m = _MOSAIC.match(p.name)
            if m:
                mosaics.setdefault(int(m.group(2)), {})[m.group(1)] = p
    inp["mosaics"] = mosaics
    tdir = idir or folder / "images"
    inp["transform"] = (tdir / TRANSFORM) if (tdir / TRANSFORM).exists() else None
    inp["image_manifest"] = (tdir / "manifest.json") if (tdir / "manifest.json").exists() else None
    return inp


def read_transform(path: Path) -> np.ndarray:
    """The 3 x 3 micron -> mosaic pixel matrix (rows of numbers separated by spaces or commas)."""
    rows = [[float(v) for v in re.split(r"[\s,]+", ln.strip()) if v]
            for ln in Path(path).read_text(encoding="utf-8").splitlines() if ln.strip()]
    A = np.array(rows, np.float64)
    if A.shape != (3, 3):
        raise SystemExit(f"[merscope] {path.name}: expected a 3 x 3 matrix, got {A.shape}")
    return A


class Frame:
    """image pixel = A . (x µm, y µm, 1): the rows of a 2 x 3 affine."""

    def __init__(self, A):
        self.A = np.asarray(A, np.float64)[:2, :3]

    def to_px(self, x, y):
        x, y = np.asarray(x, np.float64), np.asarray(y, np.float64)
        a = self.A
        return a[0, 0] * x + a[0, 1] * y + a[0, 2], a[1, 0] * x + a[1, 1] * y + a[1, 2]

    def affine(self):
        return [[float(v) for v in r] for r in self.A]


def read_cell_table(path: Path) -> dict:
    names = ff.header(path)
    idc = ff.pick(names, "EntityID", "cell", "cell_id") or names[0]
    cols = {k: ff.pick(names, k) for k in ("center_x", "center_y", "min_x", "max_x", "min_y", "max_y", "volume")}
    if not (cols["center_x"] and cols["center_y"]):
        raise SystemExit(f"[merscope] {path.name}: no center_x / center_y columns (found {names[:8]} ...)")
    t = ff.read_columns(path, floats=[v for v in cols.values() if v], strings=[idc])
    out = {k: t[v] for k, v in cols.items() if v}
    out["id"] = [s.strip() for s in t[idc]]
    return out


def cell_diameter_um(table: dict) -> tuple[float | None, str]:
    """A dot size for the cells: MERSCOPE reports a volume and a bounding box, not an area, so
    the diameter of the circle with the bounding box's inscribed-ellipse area is used."""
    if all(k in table for k in ("min_x", "max_x", "min_y", "max_y")):
        w = np.asarray(table["max_x"]) - np.asarray(table["min_x"])
        h = np.asarray(table["max_y"]) - np.asarray(table["min_y"])
        d = np.sqrt(np.clip(w, 0, None) * np.clip(h, 0, None))
        d = d[np.isfinite(d) & (d > 0)]
        if len(d):
            return float(np.median(d)), "median of sqrt(bounding-box width x height)"
    return None, "default"


def read_cells(inp: dict, log=print) -> dict | None:
    """Cell centroids (µm) in count-matrix order + the gene matrix, or None when either is missing."""
    if not inp.get("matrix") or not inp.get("cells"):
        missing = ([] if inp.get("matrix") else ["cell_by_gene.csv"]) + ([] if inp.get("cells") else ["cell_metadata.csv"])
        log(f"[merscope] no cells layer: {' and '.join(missing)} not found")
        return None
    t0 = time.time()
    names = ff.header(inp["matrix"])
    idc = ff.pick(names, "cell", "EntityID", "cell_id") or names[0]
    X, ids, feats = ff.read_count_matrix(inp["matrix"], str_ids=(idc,))
    bc = [s.strip() for s in ids[idc]]
    keep = np.array([not is_control(g) for g in feats])
    X = X[:, np.flatnonzero(keep)].tocsr()
    genes = [g for g, k in zip(feats, keep) if k]
    tab = read_cell_table(inp["cells"])
    row = {c: i for i, c in enumerate(tab["id"])}
    hit = np.array([row.get(b, -1) for b in bc], np.int64)
    ok = hit >= 0
    xy = np.full((len(bc), 2), np.nan, np.float64)
    xy[ok, 0] = np.asarray(tab["center_x"], np.float64)[hit[ok]]
    xy[ok, 1] = np.asarray(tab["center_y"], np.float64)[hit[ok]]
    diam, how = cell_diameter_um(tab)
    log(f"[merscope] cells: {len(bc):,} in {Path(inp['matrix']).name}, {len(tab['id']):,} in "
        f"{Path(inp['cells']).name}, {int(ok.sum()):,} matched; {len(genes)} genes ({int((~keep).sum())} "
        f"Blank barcodes left out), nnz {X.nnz:,}  ({time.time() - t0:.0f}s)")
    if int((~ok).sum()):
        log(f"[merscope]   WARNING {int((~ok).sum())} matrix cells are not in the cell table: no position")
    if diam:
        log(f"[merscope] cell dot size {diam:.1f} µm ({how})")
    return {"X": X, "genes": genes, "xy_um": xy, "barcodes": bc, "nCellsTable": len(tab["id"]),
            "diameterUm": diam, "controls": [g for g, k in zip(feats, keep) if not k]}


class Plan:
    """What prepare.py needs: the image source and pixel size now, the layers after the pyramid."""

    def __init__(self):
        self.source = None
        self.pixel_size_um = None
        self.manifest_source = {}
        self.frame = None
        self.background = None
        self.scan = None
        self._writer = None

    def _tx_columns(self):
        names = ff.header(self.inp["transcripts"])
        g, x, y = ff.pick(names, "gene"), ff.pick(names, "global_x"), ff.pick(names, "global_y")
        if not (g and x and y):
            raise SystemExit(f"[merscope] {Path(self.inp['transcripts']).name}: needs gene, global_x and global_y "
                             f"columns (found {names})")
        return g, x, y

    def _scan(self, log):
        if self.scan is None:
            g, x, y = self._tx_columns()
            self.scan = ff.scan_transcripts(self.inp["transcripts"], g, x, y, log=log, tag="merscope")
        return self.scan

    def _stream(self, size, um, log, density=None):
        g, x, y = self._tx_columns()
        names = self._scan(log)["names"]
        keep = np.array([bool(n) and not is_control(n) for n in names])
        log(f"[merscope] transcripts from {Path(self.inp['transcripts']).name}: {int(keep.sum())} genes, "
            f"{int((~keep).sum())} Blank barcodes left out")
        self._writer = ff.write_molecules(self.inp["transcripts"], g, x, y, names, keep, self.frame.to_px, size, um,
                                          log=log, tag="merscope", density=density)

    def finish(self, ds_dir: Path, manifest: dict, args, log=print):
        from expr import build_cell_layer
        um = manifest["pixelSizeUm"]
        size = (manifest["image"]["width"], manifest["image"]["height"])
        manifest["frame"] = {"from": "global_x, global_y (µm)", "affine": self.frame.affine()}
        if self.z is not None:
            manifest["frame"]["zPlane"] = self.z
        if self.background:
            manifest["image"]["background"] = self.background
            manifest["image"]["channels"][0]["color"] = "#b4b4b4"
        cells = self.cells
        if cells is not None:
            px, py = self.frame.to_px(cells["xy_um"][:, 0], cells["xy_um"][:, 1])
            meta = build_cell_layer(ds_dir, np.column_stack([px, py]) * um, um, cells["X"], cells["genes"],
                                    key="cells", label="Cells", cell_diameter_um=cells["diameterUm"], log=log,
                                    extra={"platform": PLATFORM}, barcodes=cells["barcodes"])
            if not getattr(args, "no_outlines", False):
                import cellshapes
                cellshapes.attach(ds_dir, meta, cellshapes.merscope_sets(self.inp, self.frame, self.z, log=log), um,
                                  args=args, log=log)
            manifest["expr"].append(meta)
        if not args.no_transcripts:
            if self._writer is None and self.inp.get("transcripts"):
                self._stream(size, um, log)
            if self._writer is not None:
                w = self._writer
                manifest["molecules"] = w.finish(ds_dir, meta={"platform": PLATFORM, "source": "detected_transcripts",
                                                               "read": int(w.rows_read)})
            else:
                log("[merscope] no molecules layer: detected_transcripts.csv not found")


def plan(args, log=print) -> Plan:
    """Find the files, fix the coordinate frame and the pixel size, and build the image source."""
    import pyramid
    folder = Path(args.src)
    if not folder.is_dir():
        raise SystemExit(f"[merscope] --src {folder} is not a folder")
    inp = find_inputs(folder, getattr(args, "cells_dir", None))
    P = Plan()
    P.inp, P.z = inp, None

    def rel(p):
        try:
            return Path(p).relative_to(folder).as_posix()
        except ValueError:
            return str(p)
    for k, what in (("transcripts", "detected_transcripts"), ("matrix", "cell_by_gene"), ("cells", "cell_metadata"),
                    ("transform", TRANSFORM)):
        log(f"[merscope] {what}: {rel(inp[k]) if inp[k] else 'not found'}")
    mos = inp["mosaics"]
    log("[merscope] mosaic images: " + (", ".join(f"z{z}: {' '.join(sorted(s))}" for z, s in sorted(mos.items())) or "none"))

    P.cells = None if args.no_cells else read_cells(inp, log=log)
    P.manifest_source = {"folder": str(folder), **{k: str(inp[k]) for k in ("matrix", "cells", "transcripts") if inp[k]}}

    want = getattr(args, "images", "auto") or "auto"
    if want != "none" and mos and not inp["transform"]:
        log(f"[merscope] WARNING mosaic images found but no {TRANSFORM}: the images cannot be placed and are "
            f"left out (copy the file into the images folder)")
    if want != "none" and mos and inp["transform"]:
        zs = sorted(mos)
        z = getattr(args, "z", None)
        if z is None:
            z = zs[len(zs) // 2]
        elif z not in mos:
            raise SystemExit(f"[merscope] --z {z}: no mosaic for that plane (found z = {zs})")
        stains = mos[z]
        order = sorted(stains, key=lambda s: (_STAIN_ORDER.index(s.lower()) if s.lower() in _STAIN_ORDER
                                              else len(_STAIN_ORDER), s.lower()))
        A = read_transform(inp["transform"])
        if abs(A[0, 1]) > 1e-9 or abs(A[1, 0]) > 1e-9 or A[0, 0] <= 0 or A[1, 1] <= 0:
            log(f"[merscope] WARNING the transform is not a plain scale + shift ({A[:2].tolist()}): positions use "
                f"it as given; check the alignment")
        um = 1.0 / float(np.sqrt(abs(A[0, 0] * A[1, 1] - A[0, 1] * A[1, 0])))
        if args.pixel_size and abs(args.pixel_size - um) / um > 0.01:
            log(f"[merscope] --pixel-size {args.pixel_size} ignored: the transform file says {um:.4f} µm per pixel")
        P.frame = Frame(A)
        P.pixel_size_um = um
        P.z = int(z)
        P.source = pyramid.TiffFileSource([(s, stains[s]) for s in order])
        log(f"[merscope] image: z-plane {z} of {zs} ({', '.join(order)}), {P.source.width} x {P.source.height} px, "
            f"{um:.4f} µm per pixel from {TRANSFORM}")
        if inp["image_manifest"]:
            try:
                mj = json.loads(inp["image_manifest"].read_text(encoding="utf-8"))
                mw, mh = mj.get("mosaic_width_pixels"), mj.get("mosaic_height_pixels")
                if mw and mh and (int(mw), int(mh)) != (P.source.width, P.source.height):
                    log(f"[merscope] WARNING images/manifest.json says {mw} x {mh} px, the mosaic is "
                        f"{P.source.width} x {P.source.height}")
            except (ValueError, OSError):
                pass
        P.manifest_source["images"] = [str(stains[s]) for s in order]
        return P

    # ---- no images: a neutral background from the transcripts (else the cells), y down as the mosaic
    lims = []
    if inp["transcripts"] and not args.no_transcripts:
        sc = P._scan(log)
        if sc["n"]:
            lims.append((sc["xlim"], sc["ylim"]))
    if P.cells is not None and np.isfinite(P.cells["xy_um"]).any():
        xy = P.cells["xy_um"]
        lims.append(([float(np.nanmin(xy[:, 0])), float(np.nanmax(xy[:, 0]))],
                     [float(np.nanmin(xy[:, 1])), float(np.nanmax(xy[:, 1]))]))
    if not lims:
        raise SystemExit(f"[merscope] nothing to import from {folder}: no mosaic images, no transcripts and no cell positions")
    xl = (min(l[0][0] for l in lims), max(l[0][1] for l in lims))
    yl = (min(l[1][0] for l in lims), max(l[1][1] for l in lims))
    bg_um = getattr(args, "background_um", None) or ff.BACKGROUND_UM
    W, H, bg_um, x0, y0 = ff.background_grid(xl, yl, 1.0, bg_um=bg_um)
    P.frame = Frame([[1.0 / bg_um, 0.0, -x0 / bg_um], [0.0, 1.0 / bg_um, -y0 / bg_um]])
    P.pixel_size_um = bg_um
    dens = ff.Density(W, H)
    if inp["transcripts"] and not args.no_transcripts:
        P._stream((W, H), bg_um, log, density=dens)
        label, sigma = "Transcript density (no image)", 1.5 / bg_um
        P.background = "transcript density"
    else:
        xy = P.cells["xy_um"]
        dens.add(*P.frame.to_px(xy[:, 0], xy[:, 1]))
        label, sigma = "Cell density (no image)", 5.0 / bg_um
        P.background = "cell density"
    log(f"[merscope] no images{' used (--images none)' if want == 'none' else ' found'}: background = "
        f"{P.background} at {bg_um:g} µm per pixel, {W} x {H} px ({dens.n:,} points)")
    P.source = ff.ArraySource([label], dens.image(sigma))
    P.manifest_source["images"] = None
    return P
