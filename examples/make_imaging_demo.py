#!/usr/bin/env python
r"""Make small SYNTHETIC CosMx- and MERSCOPE-shaped folders for the smoke test.

Nothing here comes from a tissue, an instrument or a run: every pixel, cell and molecule
is drawn by a seeded random generator from a made-up layout, so the same --seed always
gives the same files.

    python examples/make_imaging_demo.py --out demo/raw

writes

    <out>/cosmx/            everything prepare.py cosmx reads, with the older column names:
        Demo_exprMat_file.csv           fov, cell_ID, 14 genes + NegPrb1, NegPrb2, SystemControl1;
                                        one cell_ID 0 row per FOV (transcripts in no cell)
        Demo_metadata_file.csv          fov, cell_ID, Area, CenterX/Y_local_px, CenterX/Y_global_px,
                                        Mean.* ... (rows in a different order from the matrix)
        Demo_fov_positions_file.csv     fov, x_px, y_px
        Demo_tx_file.csv                fov, cell_ID, x/y_global_px, x/y_local_px, z, target, CellComp
        CellComposite/CellComposite_F001..F004.jpg   800 x 600 RGB per FOV (nuclei in blue)
        CellLabels/CellLabels_F001..F004.tif         uint16 labels = cell_ID (a disk of 7 µm per cell, cut at
                                                     its own FOV's edge): the outlines are traced from these
                            four FOVs on a 2 x 2 grid, 2 px (x) and 3 px (y) apart, global y pointing
                            UP and y_local measured from the BOTTOM edge of each FOV image
    <out>/cosmx_new/        the newer names and containers, no cell labels:
        Demo_exprMat_file.csv.gz        fov, cell_ID, cell, genes + Negative1, Negative2, SystemControl1
        Demo_metadata_file.csv.gz       ..., cell (c_1_<fov>_<cell_ID>), no global centroids
        Demo_fov_positions_file.csv.gz  FOV, x_global_px, y_global_px, x_global_mm, y_global_mm
        Demo_tx_file.csv.gz             the same molecules, plus a cell column
        Morphology2D/20260101_000000_S1_C902_P99_N99_F001..F004.TIF   two uint16 pages (membrane, DAPI)
        Demo-polygons.csv.gz            fov, cellID, cell, x/y_local_px, x/y_global_px: a 12-gon of 7 µm per cell
    <out>/cosmx_tx/         Demo_tx_file.parquet + Demo_fov_positions_file.csv only (no images, no cells)

    <out>/merscope/         everything prepare.py merscope reads, older layout:
        images/mosaic_DAPI_z0..z2.tif, mosaic_PolyT_z0..z2.tif   1600 x 1200 uint16 (z1 is in focus)
        images/micron_to_mosaic_pixel_transform.csv, images/manifest.json
        cell_by_gene.csv                cell, 14 genes + Blank-1, Blank-2
        cell_metadata.csv               (unnamed id), fov, volume, center_x, center_y, min_x, max_x, min_y, max_y
        detected_transcripts.csv        (unnamed index), barcode_id, global_x, global_y, global_z, x, y, fov,
                                        gene, transcript_id
        cell_boundaries/feature_data_7.hdf5   featuredata/<id>/zIndex_<0..2>/p_0/coordinates (µm): a 16-gon of
                                        the cell's radius in z1, 0.8 x that in z0 and z2; cell 3 only in z0,
                                        cell 4 not at all
    <out>/merscope_noimg/   no images; the cell tables one folder down with the newer names
                            (vpt/cell_by_gene.csv, vpt/cell_metadata.csv with EntityID) and
                            detected_transcripts.csv.gz re-saved with a second index column and cell_id;
                            vpt/cellpose_micron_space.parquet (EntityID, ZIndex, WKB Geometry: the same
                            outlines; cell 0 in two parts, cell 1 with a hole, cell 2 with a crumb, cell 3
                            empty in z1, cell 4 missing) and a cellpose_mosaic_space.parquet that must not be used
    <out>/merscope_pq/      images + detected_transcripts.parquet only (no cell tables)
    <out>/truth_imaging.json, truth_imaging.npz   what the test checks against

Layout (shared by both platforms): an elliptical "tissue" of 800 x 600 µm at 0.5 µm per pixel
with ~840 cells on a jittered 20 µm grid, three cell types in vertical bands, each with four
genes of its own (SynA1-4, SynB1-4, SynC1-4) plus two shared ones (SynHk1-2). Each counted
molecule sits within 5 µm of its cell's centre; a few hundred lie between cells (in no cell); a
few are control codes. Nuclei are bright disks at the cell centres, so the test can check
that cells and molecules land on them.
"""
from __future__ import annotations

import argparse
import gzip
import json
import shutil
from pathlib import Path

import numpy as np

PX_UM = 0.5
W0, H0 = 1600, 1200                      # the layout canvas (pixels)
GENES = [f"Syn{t}{i}" for t in "ABC" for i in range(1, 5)] + ["SynHk1", "SynHk2"]
# CosMx: 2 x 2 FOVs of 800 x 600 px, 2 px apart in x and 3 px in y; global origin far from zero
FOV_W, FOV_H, GAP_X, GAP_Y = 800, 600, 2, 3
COSMX_W, COSMX_H = 2 * FOV_W + GAP_X, 2 * FOV_H + GAP_Y
COSMX_GX0, COSMX_GY0 = 1000.25, 50000.5
COSMX_CONTROLS = ["NegPrb1", "NegPrb2", "SystemControl1"]
COSMX_CONTROLS_NEW = ["Negative1", "Negative2", "SystemControl1"]
# MERSCOPE: µm = origin + pixel * 0.5; mosaic pixel = 2 * µm - 2 * origin
MER_X0, MER_Y0 = 5000.0, 2000.0
MER_CONTROLS = ["Blank-1", "Blank-2"]


def layout(seed: int) -> dict:
    """Cells (canvas pixels), their types, and every molecule (gene id, canvas pixels, cell or -1).
    Gene ids >= len(GENES) are control codes (three of them)."""
    rng = np.random.default_rng(seed)
    gx, gy = np.meshgrid(np.arange(10, 800, 20.0), np.arange(10, 600, 20.0))
    xy = np.column_stack([gx.ravel(), gy.ravel()]) + rng.uniform(-4, 4, (gx.size, 2))
    inside = ((xy[:, 0] - 400) / 380) ** 2 + ((xy[:, 1] - 300) / 280) ** 2 <= 1
    xy = xy[inside]
    n = len(xy)
    ctype = np.digitize(xy[:, 0], [800 / 3, 1600 / 3])
    G, Gg = len(GENES) + 3, len(GENES)
    lam = np.full((n, G), 0.4)
    for t in range(3):
        lam[ctype == t, 4 * t:4 * t + 4] = 9.0
    lam[:, 12:14] = 5.0
    lam[:, Gg:] = 0.08
    counts = rng.poisson(lam)
    mol_g, mol_xy, mol_cell = [], [], []
    for i in range(n):
        for g in np.flatnonzero(counts[i]):
            k = int(counts[i, g])
            ang = rng.uniform(0, 2 * np.pi, k)
            rad = np.minimum(np.abs(rng.normal(0, 2.5, k)), 5.0)
            mol_xy.append(xy[i] + np.column_stack([rad * np.cos(ang), rad * np.sin(ang)]))
            mol_g.append(np.full(k, g))
            mol_cell.append(np.full(k, i))
    m_free = 400
    mol_xy.append(np.column_stack([rng.uniform(40, 760, m_free), rng.uniform(40, 560, m_free)]))
    mol_g.append(rng.integers(0, Gg, m_free))
    mol_cell.append(np.full(m_free, -1))
    mg = np.concatenate(mol_g).astype(np.int64)
    mxy = np.concatenate(mol_xy) / PX_UM                         # canvas pixels
    mc = np.concatenate(mol_cell).astype(np.int64)
    order = rng.permutation(len(mg))                             # files are not sorted by gene
    return {"xy_px": xy / PX_UM, "ctype": ctype, "counts": counts, "mol_gene": mg[order], "mol_px": mxy[order],
            "mol_cell": mc[order], "radius_um": rng.uniform(7.0, 9.0, n), "n": n}


def draw(xy_px: np.ndarray, rng, w: int, h: int, bright: float = 1.0):
    """(nuclear stain, cell-boundary stain) as float images: disks of 3 µm and rings at 8 µm."""
    yy, xx = np.mgrid[0:h, 0:w]
    dapi = rng.normal(300, 40, (h, w))
    bnd = rng.normal(250, 30, (h, w))
    r_nuc, r_cell = 3.0 / PX_UM, 8.0 / PX_UM
    for x, y in xy_px:
        x0, x1, y0, y1 = int(max(0, x - 20)), int(min(w, x + 21)), int(max(0, y - 20)), int(min(h, y + 21))
        dd = np.hypot(xx[y0:y1, x0:x1] - x, yy[y0:y1, x0:x1] - y)
        dapi[y0:y1, x0:x1] += bright * 6000 * (dd <= r_nuc)
        bnd[y0:y1, x0:x1] += bright * 2500 * (np.abs(dd - r_cell) <= 1.0)
    return dapi, bnd


def write_csv(path: Path, header: list[str], rows, gz: bool = False):
    op = gzip.open if gz else open
    with op(path, "wt", encoding="utf-8", newline="") as f:
        f.write(",".join(header) + "\n")
        f.write("".join(",".join(str(v) for v in r) + "\n" for r in rows))


# ---------------------------------------------------------------------------------- CosMx
def cosmx_frame(px: np.ndarray):
    """Canvas pixels -> (mosaic pixels, fov 1..4, local x, local y, global x, global y).

    The canvas is cut into four FOVs; the right / lower ones move by the gap. FOV 1 is the
    one with the lowest global y, i.e. the BOTTOM-left of the mosaic (global y points up)."""
    col = (px[:, 0] >= FOV_W).astype(int)
    row = (px[:, 1] >= FOV_H).astype(int)                        # 0 = top row of the mosaic
    mx = px[:, 0] + col * GAP_X
    my = px[:, 1] + row * GAP_Y
    fov = (1 - row) * 2 + col + 1
    gx = COSMX_GX0 + mx
    gy = COSMX_GY0 + (COSMX_H - my)
    lx = mx - col * (FOV_W + GAP_X)
    ly = (row * (FOV_H + GAP_Y) + FOV_H) - my                    # from the bottom edge of the FOV image
    return np.column_stack([mx, my]), fov, lx, ly, gx, gy


def fov_origin(fov: int):
    """Global pixel of local (0, 0) of a FOV (its bottom-left corner)."""
    col, row = (fov - 1) % 2, 1 - (fov - 1) // 2
    return COSMX_GX0 + col * (FOV_W + GAP_X), COSMX_GY0 + COSMX_H - (row * (FOV_H + GAP_Y) + FOV_H)


def make_cosmx(out: Path, L: dict, seed: int) -> dict:
    import tifffile
    from PIL import Image
    n, Gg = L["n"], len(GENES)
    cell_m, cell_fov, clx, cly, cgx, cgy = cosmx_frame(L["xy_px"])
    cid = np.zeros(n, np.int64)                                  # cell_ID restarts at 1 in every FOV
    for f in range(1, 5):
        sel = np.flatnonzero(cell_fov == f)
        cid[sel] = np.arange(1, len(sel) + 1)
    mol_m, mol_fov, mlx, mly, mgx, mgy = cosmx_frame(L["mol_px"])
    mc, mg = L["mol_cell"], L["mol_gene"]
    # a molecule belongs to its cell only inside the cell's own FOV (as the instrument segments per FOV)
    in_cell = (mc >= 0) & (mol_fov == cell_fov[np.maximum(mc, 0)])
    mol_cid = np.where(in_cell, cid[np.maximum(mc, 0)], 0)
    G = Gg + 3
    cnt = np.zeros((n, G), np.int64)
    np.add.at(cnt, (mc[in_cell], mg[in_cell]), 1)
    bg = np.zeros((5, G), np.int64)                              # the cell_ID 0 row of each FOV
    np.add.at(bg, (mol_fov[~in_cell], mg[~in_cell]), 1)
    area_px = np.pi * (L["radius_um"] / PX_UM) ** 2

    dapi, bnd = draw(cell_m, np.random.default_rng([seed, 1]), COSMX_W, COSMX_H)
    lab = np.zeros((COSMX_H, COSMX_W), np.uint16)
    yy, xx = np.mgrid[0:COSMX_H, 0:COSMX_W]
    for i, (x, y) in enumerate(cell_m):
        # a cell is segmented inside its own FOV only (its id means nothing in the next FOV's image)
        bx, by = ((cell_fov[i] - 1) % 2) * (FOV_W + GAP_X), (1 - (cell_fov[i] - 1) // 2) * (FOV_H + GAP_Y)
        x0, x1 = int(max(bx, x - 16)), int(min(bx + FOV_W, x + 17))
        y0, y1 = int(max(by, y - 16)), int(min(by + FOV_H, y + 17))
        dd = np.hypot(xx[y0:y1, x0:x1] - x, yy[y0:y1, x0:x1] - y)
        lab[y0:y1, x0:x1][dd <= 7.0 / PX_UM] = cid[i]
    blue = np.clip(dapi / 6300.0 * 255, 0, 255).astype(np.uint8)
    green = np.clip(bnd / 2750.0 * 200, 0, 255).astype(np.uint8)
    red = np.clip(np.random.default_rng([seed, 2]).normal(12, 3, blue.shape), 0, 255).astype(np.uint8)
    rgb = np.dstack([red, green, blue])

    def fov_box(f):
        col, row = (f - 1) % 2, 1 - (f - 1) // 2
        x0, y0 = col * (FOV_W + GAP_X), row * (FOV_H + GAP_Y)
        return slice(y0, y0 + FOV_H), slice(x0, x0 + FOV_W)

    genes_old = GENES + COSMX_CONTROLS
    genes_new = GENES + COSMX_CONTROLS_NEW
    meta_cols = ["Mean.MembraneStain", "Max.MembraneStain", "Mean.DAPI", "Max.DAPI"]
    perm = np.random.default_rng([seed, 3]).permutation(n)
    comp = np.where(in_cell, "Nuclear", "0")

    def tables(folder: Path, new: bool):
        folder.mkdir(parents=True, exist_ok=True)
        names = genes_new if new else genes_old
        gz = new
        ext = ".csv.gz" if gz else ".csv"
        q = (lambda s: s) if new else (lambda s: f'"{s}"')        # the older exports quote their headers
        cell_name = [f"c_1_{f}_{c}" for f, c in zip(cell_fov, cid)]
        rows = []
        for f in range(1, 5):                                    # by FOV: the background row, then its cells
            rows.append([f, 0] + (["c_1_%d_0" % f] if new else []) + bg[f].tolist())
            for i in np.flatnonzero(cell_fov == f):
                rows.append([f, cid[i]] + ([cell_name[i]] if new else []) + cnt[i].tolist())
        write_csv(folder / f"Demo_exprMat_file{ext}", [q(c) for c in ["fov", "cell_ID"] + (["cell"] if new else []) + names],
                  rows, gz)
        head = ["fov", "cell_ID"] + (["cell"] if new else []) + ["Area", "AspectRatio", "CenterX_local_px", "CenterY_local_px"] + \
               ([] if new else ["CenterX_global_px", "CenterY_global_px"]) + ["Width", "Height"] + meta_cols
        rows = []
        for i in perm:
            r = [cell_fov[i], cid[i]] + ([cell_name[i]] if new else []) + \
                [f"{area_px[i]:.0f}", "1.0", f"{clx[i]:.3f}", f"{cly[i]:.3f}"] + \
                ([] if new else [f"{cgx[i]:.6f}", f"{cgy[i]:.6f}"]) + [32, 32, 2500, 5000, 3000, 6300]
            rows.append(r)
        write_csv(folder / f"Demo_metadata_file{ext}", [q(c) for c in head], rows, gz)
        if new:
            write_csv(folder / f"Demo_fov_positions_file{ext}",
                      ["FOV", "x_global_px", "y_global_px", "x_global_mm", "y_global_mm"],
                      [[f, f"{fov_origin(f)[0]:.6f}", f"{fov_origin(f)[1]:.6f}", f"{fov_origin(f)[0] * PX_UM / 1000:.9f}",
                        f"{fov_origin(f)[1] * PX_UM / 1000:.9f}"] for f in range(1, 5)], gz)
        else:
            write_csv(folder / f"Demo_fov_positions_file{ext}", [q(c) for c in ["fov", "x_px", "y_px"]],
                      [[f, f"{fov_origin(f)[0]:.6f}", f"{fov_origin(f)[1]:.6f}"] for f in range(1, 5)], gz)
        head = ["fov", "cell_ID"] + (["cell"] if new else []) + ["x_global_px", "y_global_px", "x_local_px", "y_local_px",
                                                                  "z", "target", "CellComp"]
        rows = ([mol_fov[j], mol_cid[j]] + ([f"c_1_{mol_fov[j]}_{mol_cid[j]}"] if new else []) +
                [f"{mgx[j]:.6f}", f"{mgy[j]:.6f}", f"{mlx[j]:.4f}", f"{mly[j]:.4f}", j % 8, (q(names[mg[j]])), q(comp[j])]
                for j in range(len(mg)))
        write_csv(folder / f"Demo_tx_file{ext}", [q(c) for c in head], rows, gz)

    old = out / "cosmx"
    new = out / "cosmx_new"
    txo = out / "cosmx_tx"
    for d in (old, new, txo):
        if d.exists():
            shutil.rmtree(d)
    tables(old, False)
    (old / "CellComposite").mkdir()
    (old / "CellLabels").mkdir()
    tables(new, True)
    # the newer exports' polygons file: a 12-gon of 7 µm round every cell, in local and global pixels
    ang = 2 * np.pi * np.arange(12) / 12
    poly_m = cell_m[:, None, :] + (7.0 / PX_UM) * np.stack([np.cos(ang), np.sin(ang)], 1)[None]     # mosaic pixels
    rows = []
    for i in range(n):
        f, c = cell_fov[i], cid[i]
        col, row = (f - 1) % 2, 1 - (f - 1) // 2
        for mx, my in poly_m[i]:
            rows.append([f, c, f"c_1_{f}_{c}", f"{mx - col * (FOV_W + GAP_X):.3f}", f"{row * (FOV_H + GAP_Y) + FOV_H - my:.3f}",
                         f"{COSMX_GX0 + mx:.4f}", f"{COSMX_GY0 + (COSMX_H - my):.4f}"])
    write_csv(new / "Demo-polygons.csv.gz", ["fov", "cellID", "cell", "x_local_px", "y_local_px", "x_global_px", "y_global_px"],
              rows, True)
    (new / "Morphology2D").mkdir()
    for f in range(1, 5):
        ys, xs = fov_box(f)
        Image.fromarray(np.ascontiguousarray(rgb[ys, xs])).save(old / "CellComposite" / f"CellComposite_F{f:03d}.jpg", quality=92)
        tifffile.imwrite(str(old / "CellLabels" / f"CellLabels_F{f:03d}.tif"), lab[ys, xs], compression="zlib",
                         software="make_imaging_demo.py (synthetic)")
        stack = np.stack([np.clip(bnd[ys, xs], 0, 65535), np.clip(dapi[ys, xs], 0, 65535)]).astype(np.uint16)
        tifffile.imwrite(str(new / "Morphology2D" / f"20260101_000000_S1_C902_P99_N99_F{f:03d}.TIF"), stack,
                         photometric="minisblack", metadata=None, software="make_imaging_demo.py (synthetic)")
    # transcripts as Parquet + FOV positions only
    import pyarrow as pa
    import pyarrow.parquet as pq
    txo.mkdir(parents=True)
    shutil.copyfile(old / "Demo_fov_positions_file.csv", txo / "Demo_fov_positions_file.csv")
    t = pa.table({"fov": pa.array(mol_fov), "cell_ID": pa.array(mol_cid), "x_global_px": pa.array(mgx),
                  "y_global_px": pa.array(mgy), "x_local_px": pa.array(mlx), "y_local_px": pa.array(mly),
                  "z": pa.array(np.arange(len(mg)) % 8), "target": pa.array([genes_old[g] for g in mg]),
                  "CellComp": pa.array(comp.tolist())}).replace_schema_metadata({"note": "SYNTHETIC"})
    pq.write_table(t, str(txo / "Demo_tx_file.parquet"), row_group_size=5000)

    return {"cells_px": cell_m, "fov": cell_fov, "cell_ID": cid, "counts": cnt[:, :Gg], "mol_px": mol_m, "poly_px": poly_m,
            "labels": lab,
            "mol_global": np.column_stack([mgx, mgy]), "control": mg >= Gg,
            "areaPx2": area_px, "nBackgroundRows": 4}


# ------------------------------------------------------------------------------- MERSCOPE
def make_merscope(out: Path, L: dict, seed: int) -> dict:
    import tifffile
    n, Gg = L["n"], len(GENES)
    mg_all, mc_all = L["mol_gene"], L["mol_cell"]
    use = mg_all < Gg + 2                                        # two Blank barcodes (the third control is CosMx's)
    mg, mc, mpx = mg_all[use], mc_all[use], L["mol_px"][use]
    names = GENES + MER_CONTROLS
    G = len(names)
    cnt = np.zeros((n, G), np.int64)
    np.add.at(cnt, (mc[mc >= 0], mg[mc >= 0]), 1)
    xy_um = np.column_stack([MER_X0, MER_Y0]) + L["xy_px"] * PX_UM
    mol_um = np.column_stack([MER_X0, MER_Y0]) + mpx * PX_UM
    r = L["radius_um"]
    ids_old = [str(110883424764611924400221639916314253469 + 7919 * i) for i in range(n)]     # 39-digit ids
    ids_new = [str(4057517800002100001 + 3 * i) for i in range(n)]
    perm = np.random.default_rng([seed, 4]).permutation(n)

    def cell_tables(folder: Path, ids, new: bool):
        folder.mkdir(parents=True, exist_ok=True)
        write_csv(folder / "cell_by_gene.csv", ["cell"] + names, ([ids[i]] + cnt[i].tolist() for i in range(n)))
        head = (["EntityID"] if new else [""]) + ["fov", "volume", "center_x", "center_y", "min_x", "max_x", "min_y", "max_y"]
        write_csv(folder / "cell_metadata.csv", head,
                  ([ids[i], "" if new else 7, f"{4 / 3 * np.pi * r[i] ** 3:.3f}", f"{xy_um[i, 0]:.6f}", f"{xy_um[i, 1]:.6f}",
                    f"{xy_um[i, 0] - r[i]:.6f}", f"{xy_um[i, 0] + r[i]:.6f}", f"{xy_um[i, 1] - r[i]:.6f}",
                    f"{xy_um[i, 1] + r[i]:.6f}"] for i in perm))

    def transcripts(path: Path, resaved: bool):
        head = ["", "barcode_id", "global_x", "global_y", "global_z", "x", "y", "fov", "gene", "transcript_id"]
        if resaved:                                              # a second row-number column in front, and cell_id
            head = ["", "Unnamed: 0"] + head[1:] + ["cell_id"]
        rows = (([j, j % 97] if resaved else [j]) + [mg[j], f"{mol_um[j, 0]:.4f}", f"{mol_um[j, 1]:.4f}", float(j % 7),
                                                     f"{mpx[j, 0] % 400:.3f}", f"{mpx[j, 1] % 400:.3f}", 7, names[mg[j]],
                                                     f"SYNT{mg[j]:011d}"] + ([ids_new[mc[j]] if mc[j] >= 0 else -1] if resaved else [])
                for j in range(len(mg)))
        write_csv(path, head, rows, gz=path.name.endswith(".gz"))

    def images(folder: Path):
        d = folder / "images"
        d.mkdir(parents=True)
        for z, bright in ((0, 0.25), (1, 1.0), (2, 0.25)):       # z1 is the plane in focus
            dapi, poly = draw(L["xy_px"], np.random.default_rng([seed, 10 + z]), W0, H0, bright)
            for stain, img in (("DAPI", dapi), ("PolyT", poly)):
                tifffile.imwrite(str(d / f"mosaic_{stain}_z{z}.tif"), np.clip(img, 0, 65535).astype(np.uint16),
                                 photometric="minisblack", metadata=None, software="make_imaging_demo.py (synthetic)")
        s = 1.0 / PX_UM
        np.savetxt(d / "micron_to_mosaic_pixel_transform.csv", np.array([[s, 0, -MER_X0 * s], [0, s, -MER_Y0 * s], [0, 0, 1]]))
        (d / "manifest.json").write_text(json.dumps(
            {"bbox_microns": [MER_X0, MER_Y0, MER_X0 + W0 * PX_UM, MER_Y0 + H0 * PX_UM], "mosaic_width_pixels": W0,
             "mosaic_height_pixels": H0, "note": "SYNTHETIC"}, indent=1), encoding="utf-8")

    # outlines: a 16-gon of the cell's radius in the middle z-plane, 0.8 x that in the two others
    ang = 2 * np.pi * np.arange(16) / 16
    poly_um = xy_um[:, None, :] + r[:, None, None] * np.stack([np.cos(ang), np.sin(ang)], 1)[None]

    def at_z(i, z):
        return xy_um[i] + (poly_um[i] - xy_um[i]) * (1.0 if z == 1 else 0.8)

    def boundaries_hdf5(folder: Path):
        import h5py
        d = folder / "cell_boundaries"
        d.mkdir(parents=True)
        with h5py.File(d / "feature_data_7.hdf5", "w") as f:
            f.attrs["note"] = "SYNTHETIC"
            fd = f.create_group("featuredata")
            for i in range(n):
                if i == 4:
                    continue
                g = fd.create_group(ids_old[i])
                for z in ((0,) if i == 3 else (0, 1, 2)):
                    p = at_z(i, z)
                    g.create_group(f"zIndex_{z}").create_group("p_0").create_dataset("coordinates", data=np.vstack([p, p[:1]])[None])

    def wkb(parts):
        """MultiPolygon; a part is a list of rings (the first the outline, the rest holes)."""
        import struct
        b = struct.pack("<BII", 1, 6, len(parts))
        for rings in parts:
            b += struct.pack("<BII", 1, 3, len(rings))
            for rg in rings:
                rg = np.vstack([rg, rg[:1]])
                b += struct.pack("<I", len(rg)) + rg.astype("<f8").tobytes()
        return b

    def boundaries_parquet(folder: Path):
        import pyarrow as pa
        import pyarrow.parquet as pq
        small = np.stack([np.cos(ang[::2]), np.sin(ang[::2])], 1)
        rows = []
        for i in range(n):
            if i == 4:
                continue
            for z in (0, 1, 2):
                p = at_z(i, z)
                parts = [[p]]
                if i == 0:
                    parts.append([xy_um[i] + [1.6 * r[i], 0] + 0.5 * r[i] * small])       # a second part, a quarter of the area
                elif i == 1:
                    parts = [[p, (xy_um[i] + 0.25 * r[i] * small)[::-1]]]                 # a hole
                elif i == 2:
                    parts.append([xy_um[i] + [1.6 * r[i], 0] + 0.08 * r[i] * small])      # a crumb: dropped
                elif i == 3 and z == 1:
                    parts = []                                                            # nothing in the middle plane
                rows.append((int(ids_new[i]), z, wkb(parts)))
        order = np.random.default_rng([seed, 6]).permutation(len(rows))
        for name, scale in (("cellpose_micron_space.parquet", 1.0), ("cellpose_mosaic_space.parquet", None)):
            geo = [rows[j][2] for j in order]
            if scale is None:                                    # the same table in mosaic pixels: not the one to read
                geo = [wkb([[np.array([[0, 0], [9, 0], [9, 9], [0, 9]], float) + 37.0 * (j % 11)]]) for j in order]
            t = pa.table({"ID": pa.array(np.arange(len(rows), dtype=np.int64)),
                          "EntityID": pa.array(np.array([rows[j][0] for j in order], np.int64)),
                          "ZIndex": pa.array(np.array([rows[j][1] for j in order], np.int64)),
                          "Geometry": pa.array(geo, pa.binary()), "Type": pa.array(["cell"] * len(rows)),
                          "ZLevel": pa.array(np.array([1.5 * (rows[j][1] + 1) for j in order]))})
            pq.write_table(t.replace_schema_metadata({"note": "SYNTHETIC"}), str(folder / name))

    full, noimg, pqd = out / "merscope", out / "merscope_noimg", out / "merscope_pq"
    for d in (full, noimg, pqd):
        if d.exists():
            shutil.rmtree(d)
    cell_tables(full, ids_old, False)
    boundaries_hdf5(full)
    transcripts(full / "detected_transcripts.csv", False)
    images(full)
    cell_tables(noimg / "vpt", ids_new, True)
    boundaries_parquet(noimg / "vpt")
    transcripts(noimg / "detected_transcripts.csv.gz", True)
    import pyarrow as pa
    import pyarrow.parquet as pq
    pqd.mkdir(parents=True)
    shutil.copytree(full / "images", pqd / "images")
    t = pa.table({"barcode_id": pa.array(mg), "global_x": pa.array(mol_um[:, 0].astype(np.float32)),
                  "global_y": pa.array(mol_um[:, 1].astype(np.float32)), "global_z": pa.array((np.arange(len(mg)) % 7).astype(np.float32)),
                  "fov": pa.array(np.full(len(mg), 7)), "gene": pa.array([names[g] for g in mg]),
                  "transcript_id": pa.array([f"SYNT{g:011d}" for g in mg])}).replace_schema_metadata({"note": "SYNTHETIC"})
    pq.write_table(t, str(pqd / "detected_transcripts.parquet"), row_group_size=5000)
    return {"cells_um": xy_um, "counts": cnt[:, :Gg], "mol_um": mol_um, "mol_gene": mg, "ids_old": ids_old, "poly_um": poly_um,
            "ids_new": ids_new, "radius_um": r}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="demo/raw")
    ap.add_argument("--seed", type=int, default=11)
    args = ap.parse_args(argv)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    L = layout(args.seed)
    Gg = len(GENES)
    c = make_cosmx(out, L, args.seed)
    m = make_merscope(out, L, args.seed)
    np.savez_compressed(out / "truth_imaging.npz", ctype=L["ctype"], mol_gene=L["mol_gene"],
                        cosmx_cells_px=c["cells_px"], cosmx_fov=c["fov"], cosmx_cell_id=c["cell_ID"],
                        cosmx_counts=c["counts"], cosmx_mol_px=c["mol_px"], cosmx_mol_global=c["mol_global"],
                        cosmx_area_px2=c["areaPx2"], cosmx_poly_px=c["poly_px"], cosmx_labels=c["labels"],
                        merscope_poly_um=m["poly_um"],
                        merscope_cells_um=m["cells_um"], merscope_counts=m["counts"], merscope_mol_um=m["mol_um"],
                        merscope_mol_gene=m["mol_gene"], merscope_radius_um=m["radius_um"])
    truth = {"pixelSizeUm": PX_UM, "genes": GENES, "nCells": int(L["n"]),
             "cosmx": {"width": COSMX_W, "height": COSMX_H, "fovSize": [FOV_W, FOV_H], "controls": COSMX_CONTROLS,
                       "controlsNew": COSMX_CONTROLS_NEW, "nBackgroundRows": 4,
                       "read": int(len(L["mol_gene"])), "controlMolecules": int(c["control"].sum()),
                       "moleculesPerGene": {g: int((L["mol_gene"] == j).sum()) for j, g in enumerate(GENES)}},
             "merscope": {"width": W0, "height": H0, "origin": [MER_X0, MER_Y0], "controls": MER_CONTROLS,
                          "read": int(len(m["mol_gene"])), "controlMolecules": int((m["mol_gene"] >= Gg).sum()),
                          "moleculesPerGene": {g: int((m["mol_gene"] == j).sum()) for j, g in enumerate(GENES)},
                          "idsOld": m["ids_old"], "idsNew": m["ids_new"]}}
    (out / "truth_imaging.json").write_text(json.dumps(truth, indent=1), encoding="utf-8")
    print(f"[imaging demo] {L['n']} cells, {len(L['mol_gene']):,} molecules -> {out} "
          f"(cosmx, cosmx_new, cosmx_tx, merscope, merscope_noimg, merscope_pq)")


if __name__ == "__main__":
    main()
