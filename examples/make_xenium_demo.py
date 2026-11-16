#!/usr/bin/env python
r"""Make a small SYNTHETIC Xenium-shaped output folder (three variants) for the smoke test.

Nothing here comes from a tissue, an instrument or a run: every pixel, cell and molecule
is drawn by a seeded random generator from a made-up layout, so the same --seed always
gives the same files.

    python examples/make_xenium_demo.py --out demo/raw

writes

    <out>/xenium/           everything prepare.py xenium reads:
        experiment.xenium                   {"pixel_size": 0.5, ...}
        morphology_focus/ch0000_dapi.ome.tif, ch0001_boundary.ome.tif
                                            1600 x 1200 uint16, one channel per file; the first
                                            file's OME block lists both channels (as Xenium does)
        cells.csv.gz                        cell_id, x_centroid, y_centroid (µm), ..., cell_area
                                            (rows in a different order from the matrix)
        cell_feature_matrix/{barcodes,features}.tsv.gz, matrix.mtx.gz
                                            14 genes + 3 control features, counts = molecules
                                            assigned to the cell with quality >= 20
        transcripts.zarr.zip                grids/0/<x>,<y>/{location, gene_identity,
                                            quality_score, valid}, 250 µm chunks, gene_names
                                            and gene_category on the root
        cell_boundaries.csv.gz, nucleus_boundaries.csv.gz
                                            cell_id, vertex_x, vertex_y (µm), label_id; closed rings in
                                            cells.csv order, every other one anticlockwise. Two cells
                                            have no boundary, one has 80 vertices, three polygons belong
                                            to cells that are not in the matrix; two cells have two
                                            nuclei, four have none
    <out>/xenium_pq/        stains + transcripts.parquet only (no cell files): the newer
                            transcript table, same molecules
    <out>/xenium_nt/        stains + cells.parquet + cell_feature_matrix.h5 only (no transcripts),
                            and the same outlines as cell_boundaries.parquet / nucleus_boundaries.parquet
    <out>/../outline_files/ outlines alone, beside <out> (not datasets, so not among the dataset folders):
        xenium_zarr/        cells.zarr.zip: the same outlines as polygon_sets/{0 nucleus, 1 cell}
                            /{cell_index, num_vertices, vertices} with cell_id (prefix, suffix)
        xenium_zarr_v1/     cells.zarr.zip in the older layout (polygon_vertices, 13 vertices at most)
        xenium_own/         "a lab's own" outlines for the same cells: cells.geojson (µm; one cell in two
                            parts, one with a hole, one with a crumb) and nuclei.csv (a plain vertex table)
    <out>/truth_xenium.json, truth_xenium.npz   what the test checks against

Layout: an elliptical "tissue" of 800 x 600 µm (0.5 µm/px) with ~940 cells on a jittered
20 µm grid, three cell types in vertical bands (left / middle / right), each with four
genes of its own (SynA1-4, SynB1-4, SynC1-4) plus two shared ones (SynHk1-2). Each counted
molecule sits within 5 µm of its cell's centre; 15 % have a quality below 20; a few hundred
lie between cells (in no cell); a few are control codewords; two lie just off the image.
Nuclei are bright disks on the DAPI stain at the cell centres, so the test can check that
cells and molecules land on them.
"""
from __future__ import annotations

import argparse
import gzip
import json
import shutil
from pathlib import Path

import numpy as np

PX_UM = 0.5
W, H = 1600, 1200
GENES = [f"Syn{t}{i}" for t in "ABC" for i in range(1, 5)] + ["SynHk1", "SynHk2"]
CONTROLS = [("NegControlProbe_00001", "Negative Control Probe"),
            ("NegControlCodeword_0500", "Negative Control Codeword"),
            ("UnassignedCodeword_0001", "Unassigned Codeword")]
CHUNK_UM = 250.0
QV_MIN = 20.0


def cell_id(i: int) -> str:
    """Xenium-style id: base-16 digits written as the letters a-p, then '-1'."""
    s = ""
    for _ in range(8):
        s = chr(ord("a") + i % 16) + s
        i //= 16
    return s + "-1"


def ome_xml(names, w, h):
    ch = "".join(f'<Channel ID="Channel:{i}" Name="{n}" SamplesPerPixel="1"/>' for i, n in enumerate(names))
    td = "".join(f'<TiffData FirstZ="0" FirstT="0" FirstC="{i}" PlaneCount="1"/>' for i in range(len(names)))
    return ('<?xml version="1.0" encoding="UTF-8"?>'
            '<OME xmlns="http://www.openmicroscopy.org/Schemas/OME/2016-06">'
            '<Image ID="Image:0" Name="SYNTHETIC demo xenium">'
            f'<Pixels ID="Pixels:0" DimensionOrder="XYZCT" Type="uint16" SizeX="{w}" SizeY="{h}" SizeZ="1" '
            f'SizeC="{len(names)}" SizeT="1" PhysicalSizeX="{PX_UM}" PhysicalSizeXUnit="&#181;m" '
            f'PhysicalSizeY="{PX_UM}" PhysicalSizeYUnit="&#181;m">{ch}{td}</Pixels></Image></OME>')


def write_stains(folder: Path, xy_um: np.ndarray, rng):
    import tifffile
    d = folder / "morphology_focus"
    d.mkdir(parents=True, exist_ok=True)
    yy, xx = np.mgrid[0:H, 0:W]
    dapi = rng.normal(300, 40, (H, W))
    bnd = rng.normal(250, 30, (H, W))
    r_nuc, r_cell = 3.0 / PX_UM, 8.0 / PX_UM
    for x, y in xy_um / PX_UM:
        x0, x1, y0, y1 = int(max(0, x - 20)), int(min(W, x + 21)), int(max(0, y - 20)), int(min(H, y + 21))
        dd = np.hypot(xx[y0:y1, x0:x1] - x, yy[y0:y1, x0:x1] - y)
        dapi[y0:y1, x0:x1] += 6000 * (dd <= r_nuc)
        bnd[y0:y1, x0:x1] += 2500 * (np.abs(dd - r_cell) <= 1.0)
    names = ["DAPI", "Boundary"]
    for i, (img, stem) in enumerate([(dapi, "ch0000_dapi"), (bnd, "ch0001_boundary")]):
        tifffile.imwrite(str(d / f"{stem}.ome.tif"), np.clip(img, 0, 65535).astype(np.uint16),
                         photometric="minisblack", tile=(256, 256), compression="zlib", metadata=None,
                         description=ome_xml(names, W, H) if i == 0 else None,
                         software="make_xenium_demo.py (synthetic)")


def ring(centre, radius, k, phase, wobble=0.0, wphase=0.0):
    """k vertices round `centre` (µm), angle increasing: clockwise on screen (y down)."""
    a = phase + 2 * np.pi * np.arange(k) / k
    r = radius * (1 + wobble * np.sin(3 * a + wphase))
    return np.asarray(centre) + np.column_stack([r * np.cos(a), r * np.sin(a)])


def make_outlines(xy, area, seed):
    """The truth: one 12-gon per cell (mean radius from its area), one 8-gon nucleus of 3 µm (the
    DAPI disk); the special cases are listed in the returned dict."""
    rb = np.random.default_rng([seed, 5])
    n = len(xy)
    rad = np.sqrt(area / np.pi)
    ph, wp = rb.uniform(0, 2 * np.pi, n), rb.uniform(0, 2 * np.pi, n)
    cell = np.stack([ring(xy[i], rad[i], 12, ph[i], 0.06, wp[i]) for i in range(n)])
    nuc = np.stack([ring(xy[i], 3.0, 8, ph[i]) for i in range(n)])
    sp = {"noBoundary": [7, 11], "dense": 5, "denseN": 80, "binucleate": [20, 21], "noNucleus": [30, 31, 32, 33],
          "notInMatrix": 3, "ownMulti": 40, "ownHole": 41, "ownCrumb": 42}
    dense = ring(xy[sp["dense"]], rad[sp["dense"]], sp["denseN"], 0.3)
    twin = {i: [ring(xy[i] + [-2.5, 0], 2.0, 8, 0.0), ring(xy[i] + [2.5, 0], 2.0, 8, 0.5)] for i in sp["binucleate"]}
    return {"cell": cell, "nucleus": nuc, "dense": dense, "twin": twin, "special": sp}


def outline_polygons(O, barcodes, order):
    """(cell polygons, nucleus polygons) as lists of (cell id, ring) in file order."""
    sp = O["special"]
    cells, nuclei = [], []
    for j, i in enumerate(order):
        i = int(i)
        if i not in sp["noBoundary"]:
            r = O["dense"] if i == sp["dense"] else O["cell"][i]
            cells.append((barcodes[i], r[::-1] if j % 2 else r))
        if i in sp["binucleate"]:
            nuclei += [(barcodes[i], r) for r in O["twin"][i]]
        elif i not in sp["noNucleus"]:
            nuclei.append((barcodes[i], O["nucleus"][i][::-1] if j % 2 else O["nucleus"][i]))
    for k in range(sp["notInMatrix"]):                           # cells the matrix does not hold
        cells.append((cell_id(900000 + k), ring([50.0 + 20 * k, 50.0], 6.0, 12, 0.0)))
    return cells, nuclei


def write_boundary_tables(folder: Path, polys, stem: str, parquet: bool):
    cid, vx, vy, lab = [], [], [], []
    for k, (c, r) in enumerate(polys):
        r = np.vstack([r, r[:1]])                                # closed, as Xenium writes them
        cid += [c] * len(r)
        vx += r[:, 0].tolist()
        vy += r[:, 1].tolist()
        lab += [k + 1] * len(r)
    if parquet:
        import pyarrow as pa
        import pyarrow.parquet as pq
        t = pa.table({"cell_id": pa.array(cid), "vertex_x": pa.array(np.array(vx, np.float32)),
                      "vertex_y": pa.array(np.array(vy, np.float32)), "label_id": pa.array(np.array(lab, np.int32))})
        pq.write_table(t.replace_schema_metadata({"note": "SYNTHETIC"}), str(folder / f"{stem}.parquet"))
        return
    with gzip.open(folder / f"{stem}.csv.gz", "wt", encoding="utf-8", newline="") as f:
        f.write('"cell_id","vertex_x","vertex_y","label_id"\n')
        f.write("".join(f'"{c}",{x:.4f},{y:.4f},{k}\n' for c, x, y, k in zip(cid, vx, vy, lab)))


def write_cells_zarr(folder: Path, O, barcodes, old_layout: bool):
    """cells.zarr.zip with the outlines only, in matrix order (cell_index = row)."""
    import zarr
    folder.mkdir(parents=True, exist_ok=True)
    sp = O["special"]
    n = len(barcodes)
    tr = str.maketrans("abcdefghijklmnop", "0123456789abcdef")
    ids = np.array([[int(b.split("-")[0].translate(tr), 16), int(b.split("-")[1])] for b in barcodes], np.uint32)
    store = zarr.ZipStore(str(folder / "cells.zarr.zip"), mode="w")
    root = zarr.group(store=store)
    root.attrs.update({"major_version": 1 if old_layout else 6, "minor_version": 0, "number_cells": n,
                       "spatial_units": "microns", "note": "SYNTHETIC"})
    root.create_dataset("cell_id", data=ids)
    cells = [(i, O["cell"][i] if old_layout or i != sp["dense"] else O["dense"]) for i in range(n) if i not in sp["noBoundary"]]
    nuclei = []
    for i in range(n):
        if i in sp["binucleate"] and not old_layout:
            nuclei += [(i, r) for r in O["twin"][i]]
        elif i not in sp["noNucleus"] and i not in sp["binucleate"]:
            nuclei.append((i, O["nucleus"][i]))

    def packed(polys, vmax):
        v = np.zeros((len(polys), 2 * vmax), np.float32)
        nv = np.zeros(len(polys), np.int32)
        for k, (_i, r) in enumerate(polys):
            r = np.vstack([r, r[:1]])
            v[k, :2 * len(r)] = r.ravel()
            nv[k] = len(r)
        return v, nv
    if old_layout:                                               # (2, cells, 26): 0 = nucleus, 1 = cell; no polygon = 0 vertices
        pv, pn = np.zeros((2, n, 26), np.float32), np.zeros((2, n), np.int32)
        for s, polys in ((0, nuclei), (1, cells)):
            v, nv = packed(polys, 13)
            rows = [i for i, _r in polys]
            pv[s, rows], pn[s, rows] = v, nv
        root.create_dataset("polygon_vertices", data=pv)
        root.create_dataset("polygon_num_vertices", data=pn)
    else:
        ps = root.create_group("polygon_sets")
        ps.attrs.update({"number_polygon_sets": 2, "polygon_set_names": ["nucleus", "cell"],
                         "polygon_set_display_names": ["Nucleus boundaries", "Cell boundaries"]})
        for s, polys in ((0, nuclei), (1, cells)):
            g = ps.create_group(str(s))
            v, nv = packed(polys, max(len(r) for _i, r in polys) + 1)
            g.create_dataset("cell_index", data=np.array([i for i, _r in polys], np.uint32))
            g.create_dataset("method", data=np.zeros(len(polys), np.uint32))
            g.create_dataset("num_vertices", data=nv)
            g.create_dataset("vertices", data=v)
    store.close()


def write_own_outlines(folder: Path, O, barcodes):
    """A lab's own files: GeoJSON cell polygons (µm) and a plain vertex table of nuclei."""
    folder.mkdir(parents=True, exist_ok=True)
    sp = O["special"]
    feats = []
    for i, b in enumerate(barcodes):
        r = O["cell"][i]
        closed = np.vstack([r, r[:1]]).round(4).tolist()
        c = r.mean(0)
        if i == sp["ownMulti"]:                                  # a second part, 14 % of the cell's area
            sat = ring(c + [14.0, 0], 3.0, 8, 0.0)
            geom = {"type": "MultiPolygon", "coordinates": [[closed], [np.vstack([sat, sat[:1]]).round(4).tolist()]]}
        elif i == sp["ownCrumb"]:                                # a part of 1 %: dropped
            sat = ring(c + [14.0, 0], 0.8, 8, 0.0)
            geom = {"type": "MultiPolygon", "coordinates": [[closed], [np.vstack([sat, sat[:1]]).round(4).tolist()]]}
        elif i == sp["ownHole"]:
            hole = ring(c, 2.0, 8, 0.0)[::-1]
            geom = {"type": "Polygon", "coordinates": [closed, np.vstack([hole, hole[:1]]).round(4).tolist()]}
        else:
            geom = {"type": "Polygon", "coordinates": [closed]}
        feats.append({"type": "Feature", "properties": {"cell_id": b, "other_run": f"x{i}"}, "geometry": geom})
    (folder / "cells.geojson").write_text(json.dumps({"type": "FeatureCollection", "note": "SYNTHETIC", "features": feats}),
                                          encoding="utf-8")
    with open(folder / "nuclei.csv", "w", encoding="utf-8", newline="") as f:
        f.write("cell_id,vertex_x,vertex_y\n")
        for i, b in enumerate(barcodes):
            r = np.vstack([O["nucleus"][i], O["nucleus"][i][:1]])
            f.write("".join(f"{b},{x:.4f},{y:.4f}\n" for x, y in r))


def write_experiment(folder: Path, n_cells: int):
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "experiment.xenium").write_text(json.dumps(
        {"major_version": 6, "minor_version": 0, "pixel_size": PX_UM, "num_cells": n_cells,
         "panel_organism": "SYNTHETIC", "run_name": "SYNTHETIC demo"}, indent=4), encoding="utf-8")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="demo/raw")
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args(argv)
    out = Path(args.out)
    rng = np.random.default_rng(args.seed)

    # ---- cells: jittered grid inside an ellipse, three types in vertical bands
    gx, gy = np.meshgrid(np.arange(10, 800, 20.0), np.arange(10, 600, 20.0))
    xy = np.column_stack([gx.ravel(), gy.ravel()]) + rng.uniform(-4, 4, (gx.size, 2))
    inside = ((xy[:, 0] - 400) / 380) ** 2 + ((xy[:, 1] - 300) / 280) ** 2 <= 1
    xy = xy[inside]
    n = len(xy)
    ctype = np.digitize(xy[:, 0], [800 / 3, 1600 / 3])            # 0 left, 1 middle, 2 right
    names = GENES + [c[0] for c in CONTROLS]
    G, Gg = len(names), len(GENES)
    lam = np.full((n, G), 0.4)
    for t in range(3):
        lam[ctype == t, 4 * t:4 * t + 4] = 9.0
    lam[:, 12:14] = 5.0
    lam[:, Gg:] = 0.08                                            # control codewords, rare
    counts_all = rng.poisson(lam)                                 # molecules per cell and gene (any quality)

    mol_g, mol_xy, mol_q, mol_cell = [], [], [], []
    for i in range(n):
        for g in np.flatnonzero(counts_all[i]):
            k = int(counts_all[i, g])
            ang = rng.uniform(0, 2 * np.pi, k)
            rad = np.minimum(np.abs(rng.normal(0, 2.5, k)), 5.0)
            mol_xy.append(xy[i] + np.column_stack([rad * np.cos(ang), rad * np.sin(ang)]))
            mol_g.append(np.full(k, g))
            mol_cell.append(np.full(k, i))
            mol_q.append(np.where(rng.random(k) < 0.85, rng.uniform(20, 40, k), rng.uniform(0, 19.9, k)))
    m_free = 400                                                  # between cells: no cell
    fx = np.column_stack([rng.uniform(40, 760, m_free), rng.uniform(40, 560, m_free)])
    mol_xy.append(fx); mol_g.append(rng.integers(0, Gg, m_free)); mol_cell.append(np.full(m_free, -1))
    mol_q.append(rng.uniform(20, 40, m_free))
    mol_xy.append(np.array([[-3.0, 100.0], [-1.5, 200.0]]))       # just off the image: dropped
    mol_g.append(np.array([0, 1])); mol_cell.append(np.array([-1, -1])); mol_q.append(np.array([30.0, 30.0]))
    mg = np.concatenate(mol_g).astype(np.int64)
    mxy = np.concatenate(mol_xy).astype(np.float32)
    mq = np.concatenate(mol_q).astype(np.float32)
    mc = np.concatenate(mol_cell).astype(np.int64)
    mz = rng.uniform(5, 15, len(mg)).astype(np.float32)

    # the count matrix, as Xenium makes it: molecules of quality >= 20 assigned to the cell
    good = mq >= QV_MIN
    cnt = np.zeros((n, G), np.int64)
    np.add.at(cnt, (mc[good & (mc >= 0)], mg[good & (mc >= 0)]), 1)
    barcodes = [cell_id(37 * i + 11) for i in range(n)]

    area = np.pi * rng.uniform(7.0, 9.0, n) ** 2
    full = out / "xenium"
    if full.exists():
        shutil.rmtree(full)
    write_experiment(full, n)
    write_stains(full, xy, np.random.default_rng([args.seed, 1]))
    # cells.csv.gz: real columns, rows in a different order from the matrix
    perm = np.random.default_rng([args.seed, 2]).permutation(n)
    gene_tot = cnt[:, :Gg].sum(1)
    ctrl = cnt[:, Gg:]
    with gzip.open(full / "cells.csv.gz", "wt", encoding="utf-8", newline="") as f:
        f.write('"cell_id","x_centroid","y_centroid","transcript_counts","control_probe_counts",'
                '"genomic_control_counts","control_codeword_counts","unassigned_codeword_counts",'
                '"deprecated_codeword_counts","total_counts","cell_area","nucleus_area","nucleus_count",'
                '"segmentation_method"\n')
        for i in perm:
            f.write(f'"{barcodes[i]}",{xy[i, 0]:.4f},{xy[i, 1]:.4f},{gene_tot[i]},{ctrl[i, 0]},0,{ctrl[i, 1]},'
                    f'{ctrl[i, 2]},0,{cnt[i].sum()},{area[i]:.4f},{np.pi * 9:.4f},1,"SYNTHETIC"\n')
    # outlines, in cells.csv order (not the matrix's)
    O = make_outlines(xy, area, args.seed)
    cell_polys, nuc_polys = outline_polygons(O, barcodes, perm)
    write_boundary_tables(full, cell_polys, "cell_boundaries", False)
    write_boundary_tables(full, nuc_polys, "nucleus_boundaries", False)
    md = full / "cell_feature_matrix"
    md.mkdir()
    feats = [(f"SYN{j:06d}", g, "Gene Expression") for j, g in enumerate(GENES)] + \
            [(f"CTRL{j:06d}", c, t) for j, (c, t) in enumerate(CONTROLS)]
    with gzip.open(md / "features.tsv.gz", "wt", encoding="utf-8") as f:
        f.write("".join(f"{a}\t{b}\t{c}\n" for a, b, c in feats))
    with gzip.open(md / "barcodes.tsv.gz", "wt", encoding="utf-8") as f:
        f.write("".join(b + "\n" for b in barcodes))
    nz = np.argwhere(cnt.T > 0)                                   # (feature, cell)
    with gzip.open(md / "matrix.mtx.gz", "wt", encoding="utf-8") as f:
        f.write("%%MatrixMarket matrix coordinate integer general\n%metadata_json: {\"SYNTHETIC\": true}\n")
        f.write(f"{G} {n} {len(nz)}\n")
        f.write("".join(f"{a + 1} {b + 1} {cnt[b, a]}\n" for a, b in nz))

    # transcripts.zarr.zip: Xenium's grid layout (level 0 only), molecules sorted by gene in a chunk
    import zarr
    zp = full / "transcripts.zarr.zip"
    store = zarr.ZipStore(str(zp), mode="w")
    root = zarr.group(store=store)
    root.attrs.update({"gene_names": names, "number_genes": G, "major_version": 5, "minor_version": 2,
                       "spatial_units": "micron", "name": "RnaDataset", "number_rnas": int(len(mg)),
                       "note": "SYNTHETIC"})
    cat = np.zeros((G, 9), bool)
    cat[:Gg, 0] = True
    cat[Gg + 0, 2] = cat[Gg + 1, 1] = cat[Gg + 2, 3] = True
    gc = root.create_dataset("gene_category", data=cat)
    gc.attrs["column_names"] = ["is_gene", "is_negative_control_codeword", "is_negative_control_probe",
                                "is_unassigned_codeword", "is_predesigned_gene", "is_deprecated_codeword",
                                "is_custom_gene", "is_genomic_control", "is_boosted_gene"]
    grids = root.create_group("grids")
    g0 = grids.create_group("0")
    # chunk = floor(µm / 250); the two off-image molecules (x < 0) are stored in column 0
    cx = np.maximum(np.floor(mxy[:, 0] / CHUNK_UM).astype(int), 0)
    cy = np.maximum(np.floor(mxy[:, 1] / CHUNK_UM).astype(int), 0)
    keys = []
    for a, b in sorted(set(zip(cx.tolist(), cy.tolist()))):
        sel = np.flatnonzero((cx == a) & (cy == b))
        sel = sel[np.argsort(mg[sel], kind="stable")]
        key = f"{a},{b}"
        c = g0.create_group(key)
        arrs = {"location": np.column_stack([mxy[sel], mz[sel]]).astype(np.float32),
                "gene_identity": mg[sel].astype(np.uint16).reshape(-1, 1),
                "quality_score": mq[sel].reshape(-1, 1), "valid": np.ones((len(sel), 1), np.uint8)}
        for k, v in arrs.items():
            c.create_dataset(k, data=v, chunks=v.shape)
        keys.append(key)
    grids.attrs.update({"grid_keys": [keys], "grid_size": [CHUNK_UM], "number_levels": 1})
    store.close()

    # variant 2: stains + transcripts.parquet only
    import pyarrow as pa
    import pyarrow.parquet as pq
    pqd = out / "xenium_pq"
    if pqd.exists():
        shutil.rmtree(pqd)
    write_experiment(pqd, 0)
    shutil.copytree(full / "morphology_focus", pqd / "morphology_focus")
    t = pa.table({"transcript_id": pa.array(np.arange(len(mg), dtype=np.uint64) + 281474976710656),
                  "cell_id": pa.array([barcodes[c] if c >= 0 else "UNASSIGNED" for c in mc]),
                  "overlaps_nucleus": pa.array(np.zeros(len(mg), np.uint8)),
                  "feature_name": pa.array([names[g] for g in mg]),
                  "x_location": pa.array(mxy[:, 0]), "y_location": pa.array(mxy[:, 1]),
                  "z_location": pa.array(mz), "qv": pa.array(mq),
                  "is_gene": pa.array(mg < Gg)})
    t = t.replace_schema_metadata({"note": "SYNTHETIC"})
    pq.write_table(t, str(pqd / "transcripts.parquet"), row_group_size=4000)

    # variant 3: stains + cells.parquet + cell_feature_matrix.h5, no transcripts
    import h5py
    ntd = out / "xenium_nt"
    if ntd.exists():
        shutil.rmtree(ntd)
    write_experiment(ntd, n)
    shutil.copytree(full / "morphology_focus", ntd / "morphology_focus")
    pq.write_table(pa.table({"cell_id": pa.array([barcodes[i] for i in perm]),
                             "x_centroid": pa.array(xy[perm, 0]), "y_centroid": pa.array(xy[perm, 1]),
                             "cell_area": pa.array(area[perm])}), str(ntd / "cells.parquet"))
    write_boundary_tables(ntd, cell_polys, "cell_boundaries", True)
    write_boundary_tables(ntd, nuc_polys, "nucleus_boundaries", True)
    extra = out.parent / "outline_files"
    for d in (extra / "xenium_zarr", extra / "xenium_zarr_v1", extra / "xenium_own"):
        if d.exists():
            shutil.rmtree(d)
    write_cells_zarr(extra / "xenium_zarr", O, barcodes, False)
    write_cells_zarr(extra / "xenium_zarr_v1", O, barcodes, True)
    write_own_outlines(extra / "xenium_own", O, barcodes)
    from scipy import sparse
    csc = sparse.csc_matrix(cnt.T.astype(np.int32))               # features x cells, CSC by cell
    with h5py.File(ntd / "cell_feature_matrix.h5", "w") as f:
        f.attrs["note"] = "SYNTHETIC"
        m = f.create_group("matrix")
        m.create_dataset("data", data=csc.data)
        m.create_dataset("indices", data=csc.indices.astype(np.int64))
        m.create_dataset("indptr", data=csc.indptr.astype(np.int64))
        m.create_dataset("shape", data=np.array([G, n], np.int32))
        m.create_dataset("barcodes", data=np.array([b.encode() for b in barcodes]))
        ft = m.create_group("features")
        ft.create_dataset("id", data=np.array([a.encode() for a, _, _ in feats]))
        ft.create_dataset("name", data=np.array([b.encode() for _, b, _ in feats]))
        ft.create_dataset("feature_type", data=np.array([c.encode() for _, _, c in feats]))

    # truth for the test (numpy only)
    on_img = (mxy[:, 0] >= 0) & (mxy[:, 0] <= W * PX_UM) & (mxy[:, 1] >= 0) & (mxy[:, 1] <= H * PX_UM)
    keep = good & (mg < Gg) & on_img
    np.savez_compressed(out / "truth_xenium.npz", xy_um=xy, ctype=ctype, counts=cnt[:, :Gg],
                        mol_gene=mg[keep], mol_xy_um=mxy[keep], mol_cell=mc[keep],
                        cell_poly_um=O["cell"], nucleus_poly_um=O["nucleus"], dense_poly_um=O["dense"],
                        twin_poly_um=np.stack([np.stack(O["twin"][i]) for i in O["special"]["binucleate"]]),
                        cell_area_um2=area)
    truth = {"pixelSizeUm": PX_UM, "width": W, "height": H, "nCells": int(n), "genes": GENES,
             "controls": [c for c, _ in CONTROLS],
             "moleculesPerGene": {g: int((keep & (mg == j)).sum()) for j, g in enumerate(GENES)},
             "controlMolecules": int((good & (mg >= Gg) & on_img).sum()),
             "lowQuality": int((~good).sum()), "outsideImage": int((good & (mg < Gg) & ~on_img).sum()),
             "read": int(len(mg)), "outlines": O["special"]}
    (out / "truth_xenium.json").write_text(json.dumps(truth, indent=1), encoding="utf-8")
    print(f"[xenium demo] {n} cells, {len(mg):,} molecules ({int(keep.sum()):,} kept), {G} features -> {out}")


if __name__ == "__main__":
    main()
