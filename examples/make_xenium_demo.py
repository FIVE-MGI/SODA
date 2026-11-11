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
    <out>/xenium_pq/        stains + transcripts.parquet only (no cell files): the newer
                            transcript table, same molecules
    <out>/xenium_nt/        stains + cells.parquet + cell_feature_matrix.h5 only (no transcripts)
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
                        mol_gene=mg[keep], mol_xy_um=mxy[keep], mol_cell=mc[keep])
    truth = {"pixelSizeUm": PX_UM, "width": W, "height": H, "nCells": int(n), "genes": GENES,
             "controls": [c for c, _ in CONTROLS],
             "moleculesPerGene": {g: int((keep & (mg == j)).sum()) for j, g in enumerate(GENES)},
             "controlMolecules": int((good & (mg >= Gg) & on_img).sum()),
             "lowQuality": int((~good).sum()), "outsideImage": int((good & (mg < Gg) & ~on_img).sum()),
             "read": int(len(mg))}
    (out / "truth_xenium.json").write_text(json.dumps(truth, indent=1), encoding="utf-8")
    print(f"[xenium demo] {n} cells, {len(mg):,} molecules ({int(keep.sum()):,} kept), {G} features -> {out}")


if __name__ == "__main__":
    main()
