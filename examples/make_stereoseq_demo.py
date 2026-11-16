#!/usr/bin/env python
r"""Make a small SYNTHETIC Stereo-seq (STOmics SAW) shaped output folder for the smoke test.

Nothing here comes from a tissue, a chip or a run: every pixel, DNB count and cell is drawn
by a seeded random generator from a made-up layout, so the same --seed always gives the same
files. The HDF5 group / dataset / field names and attributes follow the GEF files SAW writes
(bin GEF version 4 and the older one-`gene`-field layout; cell bin GEF version 4).

    python examples/make_stereoseq_demo.py --out demo/raw

writes

    <out>/stereoseq/
        SYNCHIP_HE_regist.tif       1400 x 1200 RGB uint8, one row per strip, deflate (as SAW's
                                    registered image): one pixel per DNB, white glass, pink tissue,
                                    a dark nucleus at every cell centre
        SYNCHIP_ssDNA_regist.tif    the same chip as a grey nuclear stain, uint16 (bright nuclei)
        SYNCHIP_HE_tissue_cut.tif   the tissue mask, uint8, 255 = tissue
        SYNCHIP.raw.gef             bin GEF version 4, bin1 only, WHOLE CHIP (tissue + background):
                                      geneExp/bin1/expression   (x u4, y u4, count u1), gene-major
                                      geneExp/bin1/gene         (geneID S64, geneName S64, offset, count)
                                      geneExp/bin1/exon         u1
                                      contour/tissueContour     JSON: the mask's border rings
                                    attrs: version, bin_type, omics, sn; minX.., resolution = 500
        SYNCHIP.tissue.gef          the older layout, tissue only, stored levels bin10 and bin50
                                    (no bin1): geneExp/binN/gene is (gene S32, offset, count).
                                    Cut to a region as STOmics' writer does it: x and y are
                                    stored minus the region's corner, kept as minX = 100, minY = 0
        SYNCHIP.tissue.gem.gz       GEMv0.2 text of the tissue DNBs, gene after gene:
                                    geneID geneName x y MIDCount ExonCount; cut the same way,
                                    with #OffsetX=100 in the header
        SYNCHIP.raw.shuffled.gem    whole chip, rows in random order, the old four columns
                                    (geneID x y MIDCount, the gene NAME in the first column)
        SYNCHIP.cellbin.gef         cell bin GEF version 4: cellBin/{cell, cellExp, gene, geneExp,
                                    cellBorder, cellTypeList, blockIndex, blockSize}
        SYNCHIP.cut.cellbin.gef     the same cells as made from the cut bin GEF: cell x is stored
                                    minus 100, with offsetX = 100 on the root
    <out>/truth_stereoseq.json, truth_stereoseq.npz   what the test checks against

Layout: a 1400 x 1200 DNB chip at 0.5 µm per DNB (700 x 600 µm). The tissue is an octagon
with an octagonal hole and a small separate island; three vertical bands (x < 500, 500 to
900, >= 900 DNB, so no bin of 20, 50 or 100 DNBs straddles two) are "cell types" TypeA / TypeB / TypeC, each with four genes of its own (SynA1-4,
SynB1-4, SynC1-4), two shared genes (SynHk1-2), a mitochondrial one (mt-Syn1), one with the
name of a real mouse gene that looks like a sample ID to deid.py (H2-T23), and two gene ids that
share the name SynDup. Background DNBs off the tissue carry a few counts too, as on
a real chip. Cells sit on a jittered 20 µm grid; a cell's counts are the DNB counts inside
its disk.
"""
from __future__ import annotations

import argparse
import gzip
import json
import shutil
from pathlib import Path

import numpy as np

SN = "SYNCHIP"
W, H = 1400, 1200
RES_NM = 500
TYPES = ["typea", "typeb", "typec"]
NAMES = ([f"Syn{t}{i}" for t in "ABC" for i in range(1, 5)] + ["SynHk1", "SynHk2", "mt-Syn1", "H2-T23", "SynDup", "SynDup"])
IDS = [f"SYNG{i + 1:011d}" for i in range(len(NAMES))]
# tissue: |dx| <= a, |dy| <= b, |dx| + |dy| <= c around (cx, cy); the hole and the island likewise
OUTER = (700, 610, 560, 470, 840)
HOLE = (520, 470, 60, 60, 90)          # ring on the tissue DNBs around the hole: its inside is cut out
ISLAND = (1290, 110, 50, 50, 75)
CELL_R = 13                            # DNBs
GRID = 40                              # DNBs between cell centres (20 µm)
CUT = (100, 0)                         # the corner the "tissue" files are cut to (a multiple of every bin size)
BAND_EDGES = (500, 900)                # x where TypeA | TypeB | TypeC meet (multiples of every bin size used)


def octagon(cx, cy, a, b, c, shrink=0):
    yy, xx = np.mgrid[0:H, 0:W]
    dx, dy = np.abs(xx - cx), np.abs(yy - cy)
    return (dx <= a - shrink) & (dy <= b - shrink) & (dx + dy <= c - shrink)


def ring(cx, cy, a, b, c):
    """The octagon's corner DNBs in drawing order (what an OpenCV contour of it holds)."""
    p, q = c - b, c - a
    return [[cx - p, cy - b], [cx + p, cy - b], [cx + a, cy - q], [cx + a, cy + q],
            [cx + p, cy + b], [cx - p, cy + b], [cx - a, cy + q], [cx - a, cy - q]]


def tissue_mask():
    return (octagon(*OUTER) & ~octagon(*HOLE, shrink=1)) | octagon(*ISLAND)


def region_of(x):
    return np.searchsorted(np.asarray(BAND_EDGES), np.asarray(x), side="right").astype(np.int64)


def draw_molecules(mask, rng):
    """(gene, x, y, count) records of the whole chip, one per (gene, DNB), gene-major."""
    G = len(NAMES)
    prof = np.full((3, G), 0.004)
    for t in range(3):
        prof[t, 4 * t:4 * t + 4] = 0.17
    prof[:, 12:14] = 0.06            # SynHk1-2
    prof[:, 14] = 0.03               # mt-Syn1
    prof[:, 15] = 0.03               # H2-T23
    prof[:, 16] = 0.02               # SynDup (first id)
    prof[:, 17] = 0.008              # SynDup (second id)
    prof /= prof.sum(1, keepdims=True)
    ty, tx = np.nonzero(mask)
    n_t = int(0.55 * len(tx))
    pick = rng.integers(0, len(tx), n_t)
    mx, my = tx[pick], ty[pick]
    reg = region_of(mx)
    mg = np.empty(n_t, np.int64)
    for t in range(3):
        sel = reg == t
        mg[sel] = rng.choice(G, int(sel.sum()), p=prof[t])
    by, bx = np.nonzero(~mask)
    n_b = int(0.02 * len(bx))
    pick = rng.integers(0, len(bx), n_b)
    gx = np.concatenate([mx, bx[pick]])
    gy = np.concatenate([my, by[pick]])
    gg = np.concatenate([mg, rng.integers(0, G, n_b)])
    key = (gg * H + gy) * W + gx
    uk, cnt = np.unique(key, return_counts=True)
    return (uk // (W * H)).astype(np.int64), (uk % W).astype(np.int64), ((uk // W) % H).astype(np.int64), cnt.astype(np.int64)


def write_image(folder: Path, mask, cells_xy, rng):
    import tifffile
    img = np.empty((H, W, 3), np.float64)
    img[:] = (236, 234, 238)
    img[mask] = (226, 150, 190)
    img += rng.normal(0, 4, img.shape)
    yy, xx = np.mgrid[0:H, 0:W]
    for x, y in cells_xy:
        x0, x1, y0, y1 = max(0, x - 8), min(W, x + 9), max(0, y - 8), min(H, y + 9)
        d = np.hypot(xx[y0:y1, x0:x1] - x, yy[y0:y1, x0:x1] - y) <= 5
        img[y0:y1, x0:x1][d] = (70, 40, 120)
    tifffile.imwrite(str(folder / f"{SN}_HE_regist.tif"), np.clip(img, 0, 255).astype(np.uint8), photometric="rgb",
                     rowsperstrip=1, compression="zlib", software="make_stereoseq_demo.py (synthetic)")
    grey = rng.normal(300, 30, (H, W)) + 500 * mask + 9000 * (img[..., 0] < 120)
    tifffile.imwrite(str(folder / f"{SN}_ssDNA_regist.tif"), np.clip(grey, 0, 65535).astype(np.uint16),
                     photometric="minisblack", rowsperstrip=16, compression="zlib",
                     software="make_stereoseq_demo.py (synthetic)")
    tifffile.imwrite(str(folder / f"{SN}_HE_tissue_cut.tif"), (mask * 255).astype(np.uint8), photometric="minisblack",
                     rowsperstrip=64, compression="zlib", software="make_stereoseq_demo.py (synthetic)")


def _a(v, dt):
    return np.array([v], dtype=dt)


def gene_major(g, x, y, c, G):
    order = np.lexsort((x, y, g))
    g, x, y, c = g[order], x[order], y[order], c[order]
    cnt = np.bincount(g, minlength=G)
    off = np.concatenate([[0], np.cumsum(cnt)[:-1]])
    return g, x, y, c, off, cnt


def write_bin_group(f, level, g, x, y, c, G, old_layout, rng, origin=(0, 0)):
    """geneExp/bin<level> from bin1 records (summed to the level's bins). x, y are counted from
    `origin`, which goes into minX / minY (maxX / maxY stay chip coordinates), as geftools writes."""
    bx, by = (x // level) * level, (y // level) * level
    key = (g * (H + level) + by) * (W + level) + bx
    uk, inv = np.unique(key, return_inverse=True)
    cs = np.bincount(inv, weights=c).astype(np.int64)
    gg, yy, xx = uk // ((W + level) * (H + level)), (uk // (W + level)) % (H + level), uk % (W + level)
    gg, xx, yy, cs, off, cnt = gene_major(gg, xx, yy, cs, G)
    grp = f.require_group("geneExp").create_group(f"bin{level}")
    cdt = "u1" if cs.max() < 256 else ("<u2" if cs.max() < 65536 else "<u4")
    e = np.empty(len(gg), dtype=[("x", "<u4"), ("y", "<u4"), ("count", cdt)])
    e["x"], e["y"], e["count"] = xx, yy, cs
    d = grp.create_dataset("expression", data=e)
    for k, v in (("minX", origin[0]), ("minY", origin[1]), ("maxX", origin[0] + int(xx.max())),
                 ("maxY", origin[1] + int(yy.max()))):
        d.attrs[k] = _a(v, "<i4")
    d.attrs["maxExp"] = _a(int(cs.max()), "<u4")
    d.attrs["resolution"] = _a(RES_NM, "<u4")
    if old_layout:
        t = np.empty(G, dtype=[("gene", "S32"), ("offset", "<u4"), ("count", "<u4")])
        t["gene"] = [n.encode() for n in NAMES_UNIQUE]
    else:
        t = np.empty(G, dtype=[("geneID", "S64"), ("geneName", "S64"), ("offset", "<u4"), ("count", "<u4")])
        t["geneID"], t["geneName"] = [i.encode() for i in IDS], [n.encode() for n in NAMES]
        ex = grp.create_dataset("exon", data=np.minimum(cs, rng.integers(0, 3, len(cs))).astype("u1"))
        ex.attrs["maxExon"] = _a(int(cs.max()), "<u4")
    t["offset"], t["count"] = off, cnt
    grp.create_dataset("gene", data=t)


# the older layout has one name per gene: the second SynDup id gets a name of its own there
NAMES_UNIQUE = NAMES[:-1] + ["SynDup2"]


def write_gefs(folder: Path, rec, mask, rng):
    import h5py
    g, x, y, c = rec
    G = len(NAMES)
    with h5py.File(folder / f"{SN}.raw.gef", "w") as f:
        f.attrs["version"] = _a(4, "<u4")
        f.attrs["bin_type"] = np.array([b"Bin"], dtype="S32")
        f.attrs["omics"] = np.array([b"Transcriptomics"], dtype="S32")
        f.attrs["sn"] = SN
        gg, xx, yy, cc, off, cnt = gene_major(g, x, y, c, G)
        grp = f.create_group("geneExp").create_group("bin1")
        e = np.empty(len(gg), dtype=[("x", "<u4"), ("y", "<u4"), ("count", "u1")])
        e["x"], e["y"], e["count"] = xx, yy, cc
        d = grp.create_dataset("expression", data=e)
        for k, v in (("minX", 0), ("minY", 0), ("maxX", W - 1), ("maxY", H - 1)):       # whole chip
            d.attrs[k] = _a(v, "<i4")
        d.attrs["maxExp"] = _a(int(cc.max()), "<u4")
        d.attrs["resolution"] = _a(RES_NM, "<u4")
        t = np.empty(G, dtype=[("geneID", "S64"), ("geneName", "S64"), ("offset", "<u4"), ("count", "<u4")])
        t["geneID"], t["geneName"] = [i.encode() for i in IDS], [n.encode() for n in NAMES]
        t["offset"], t["count"] = off, cnt
        grp.create_dataset("gene", data=t)
        ex = grp.create_dataset("exon", data=np.minimum(cc, rng.integers(0, 3, len(cc))).astype("u1"))
        ex.attrs["maxExon"] = _a(int(cc.max()), "<u4")
        contour = {"type": "GeometryCollection", "geometries": [{
            "type": "Polygon", "properties": {"label": "tissue", "uID": 0},
            # the order SAW writes: small rings first, the tissue outline last
            "coordinates": [ring(*HOLE), ring(*ISLAND), ring(*OUTER)]}]}
        f.create_group("contour").create_dataset("tissueContour", data=json.dumps(contour, separators=(",", ":")),
                                                 dtype=h5py.string_dtype())
    inside = mask[y, x]
    with h5py.File(folder / f"{SN}.tissue.gef", "w") as f:
        f.attrs["version"] = _a(2, "<u4")
        f.attrs["bin_type"] = np.array([b"Bin"], dtype="S32")
        for level in (10, 50):
            write_bin_group(f, level, g[inside], x[inside] - CUT[0], y[inside] - CUT[1], c[inside], G, True, rng, origin=CUT)


def write_gems(folder: Path, rec, mask, rng):
    g, x, y, c = rec
    inside = mask[y, x]
    gg, xx, yy, cc, _, _ = gene_major(g[inside], x[inside], y[inside], c[inside], len(NAMES))
    head = (f"#FileFormat=GEMv0.2\n#SortedBy=None\n#BinType=Bin\n#BinSize=1\n#Omics=Transcriptomics\n"
            f"#Stereo-seqChip={SN}\n#OffsetX={CUT[0]}\n#OffsetY={CUT[1]}\n")
    xx, yy = xx - CUT[0], yy - CUT[1]
    with gzip.open(folder / f"{SN}.tissue.gem.gz", "wt", encoding="utf-8", newline="\n") as f:
        f.write(head + "geneID\tgeneName\tx\ty\tMIDCount\tExonCount\n")
        f.write("".join(f"{IDS[a]}\t{NAMES[a]}\t{b}\t{d}\t{e}\t{min(e, 1)}\n" for a, b, d, e in
                        zip(gg.tolist(), xx.tolist(), yy.tolist(), cc.tolist())))
    order = rng.permutation(len(g))
    with open(folder / f"{SN}.raw.shuffled.gem", "w", encoding="utf-8", newline="\n") as f:
        f.write("#FileFormat=GEMv0.1\n#BinSize=1\n#OffsetX=0\n#OffsetY=0\ngeneID\tx\ty\tMIDCount\n")
        f.write("".join(f"{NAMES_UNIQUE[a]}\t{b}\t{d}\t{e}\n" for a, b, d, e in
                        zip(g[order].tolist(), x[order].tolist(), y[order].tolist(), c[order].tolist())))


def make_cells(mask, rec, rng):
    """Cell centres on a jittered grid inside the tissue, and the DNB records each disk holds."""
    g, x, y, c = rec
    cells = []
    for gy in range(GRID // 2, H, GRID):
        for gx in range(GRID // 2, W, GRID):
            cx, cy = gx + int(rng.integers(-5, 6)), gy + int(rng.integers(-5, 6))
            if (CELL_R < cx < W - CELL_R and CELL_R < cy < H - CELL_R and mask[cy, cx]
                    and mask[cy - CELL_R, cx] and mask[cy + CELL_R, cx] and mask[cy, cx - CELL_R] and mask[cy, cx + CELL_R]):
                cells.append((cx, cy))
    cells = np.array(cells, np.int64)
    label = np.full((H, W), -1, np.int64)
    yy, xx = np.mgrid[-CELL_R:CELL_R + 1, -CELL_R:CELL_R + 1]
    disk = (xx * xx + yy * yy) <= CELL_R * CELL_R
    for i, (cx, cy) in enumerate(cells):
        sub = label[cy - CELL_R:cy + CELL_R + 1, cx - CELL_R:cx + CELL_R + 1]
        sub[disk & (sub < 0)] = i
    lab = label[y, x]
    ok = lab >= 0
    counts = np.zeros((len(cells), len(NAMES)), np.int64)
    np.add.at(counts, (lab[ok], g[ok]), c[ok])
    dnb = np.bincount(lab[ok], minlength=len(cells))
    area = np.bincount(label[label >= 0], minlength=len(cells))
    keep = counts.sum(1) > 0
    return cells[keep], counts[keep], dnb[keep], area[keep]


def write_cellbin(folder: Path, cells, counts, dnb, area, rng, name=f"{SN}.cellbin.gef", origin=(0, 0)):
    """`origin`: the corner of the bin GEF the cells were cut from; cell x, y are counted from it."""
    import h5py
    n, G = counts.shape
    cells = cells - np.asarray(origin)
    with h5py.File(folder / name, "w") as f:
        f.attrs["version"] = _a(4, "<u4")
        f.attrs["bin_type"] = np.array([b"CellBin"], dtype="S32")
        f.attrs["omics"] = np.array([b"Transcriptomics"], dtype="S32")
        f.attrs["sn"] = SN
        f.attrs["geftool_ver"] = np.array([1, 2, 2], dtype="<u4")
        for k, v in (("resolution", RES_NM), ("offsetX", origin[0]), ("offsetY", origin[1]), ("maxX", W - 1), ("maxY", H - 1)):
            f.attrs[k] = _a(v, "<u4" if k == "resolution" else "<i4")
        cb = f.create_group("cellBin")
        ci, gi = np.nonzero(counts)
        vals = counts[ci, gi]
        ng = np.bincount(ci, minlength=n)
        cell = np.empty(n, dtype=[("id", "<u4"), ("x", "<i4"), ("y", "<i4"), ("offset", "<u4"), ("geneCount", "<u2"),
                                  ("expCount", "<u2"), ("dnbCount", "<u2"), ("area", "<u2"), ("cellTypeID", "<u2"),
                                  ("clusterID", "<u2")])
        cell["id"], cell["x"], cell["y"] = np.arange(n), cells[:, 0], cells[:, 1]
        cell["offset"] = np.concatenate([[0], np.cumsum(ng)[:-1]])
        cell["geneCount"], cell["expCount"], cell["dnbCount"], cell["area"] = ng, counts.sum(1), dnb, area
        cell["cellTypeID"], cell["clusterID"] = 0, 0
        d = cb.create_dataset("cell", data=cell)
        for k, v in (("minX", cells[:, 0].min()), ("maxX", cells[:, 0].max()), ("minY", cells[:, 1].min()),
                     ("maxY", cells[:, 1].max())):
            d.attrs[k] = _a(int(v), "<i4")
        d.attrs["medianArea"] = _a(float(np.median(area)), "<f4")
        ce = np.empty(len(ci), dtype=[("geneID", "<u4"), ("count", "<u2")])
        ce["geneID"], ce["count"] = gi, vals
        cb.create_dataset("cellExp", data=ce).attrs["maxCount"] = _a(int(vals.max()), "<u2")
        order = np.lexsort((ci, gi))
        ge = np.empty(len(ci), dtype=[("cellID", "<u4"), ("count", "<u2")])
        ge["cellID"], ge["count"] = ci[order], vals[order]
        cb.create_dataset("geneExp", data=ge).attrs["maxCount"] = _a(int(vals.max()), "<u2")
        cc = np.bincount(gi, minlength=G)
        gene = np.empty(G, dtype=[("geneID", "S64"), ("geneName", "S64"), ("offset", "<u4"), ("cellCount", "<u4"),
                                  ("expCount", "<u4"), ("maxMIDcount", "<u2")])
        gene["geneID"], gene["geneName"] = [i.encode() for i in IDS], [nm.encode() for nm in NAMES]
        gene["offset"] = np.concatenate([[0], np.cumsum(cc)[:-1]])
        gene["cellCount"], gene["expCount"], gene["maxMIDcount"] = cc, counts.sum(0), counts.max(0)
        cb.create_dataset("gene", data=gene)
        border = np.full((n, 32, 2), 32767, np.int16)                 # offsets from the centre, padded
        ang = np.linspace(0, 2 * np.pi, 16, endpoint=False)
        border[:, :16, 0] = np.rint(CELL_R * np.cos(ang))
        border[:, :16, 1] = np.rint(CELL_R * np.sin(ang))
        cb.create_dataset("cellBorder", data=border)
        cb.create_dataset("cellTypeList", data=np.array([b"default"], dtype="S32"))
        cb.create_dataset("blockIndex", data=np.array([0, n], dtype="<u4"))
        cb.create_dataset("blockSize", data=np.array([W, H, 1, 1], dtype="<u4"))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="demo/raw")
    ap.add_argument("--seed", type=int, default=11)
    args = ap.parse_args(argv)
    out = Path(args.out)
    folder = out / "stereoseq"
    if folder.exists():
        shutil.rmtree(folder)
    folder.mkdir(parents=True)
    rng = np.random.default_rng(args.seed)
    mask = tissue_mask()
    rec = draw_molecules(mask, rng)
    cells, counts, dnb, area = make_cells(mask, rec, rng)
    write_image(folder, mask, cells, rng)
    write_gefs(folder, rec, mask, rng)
    write_gems(folder, rec, mask, rng)
    write_cellbin(folder, cells, counts, dnb, area, rng)
    write_cellbin(folder, cells, counts, dnb, area, rng, name=f"{SN}.cut.cellbin.gef", origin=CUT)
    g, x, y, c = rec
    inside = mask[y, x]
    np.savez_compressed(out / "truth_stereoseq.npz", gene=g, x=x, y=y, count=c, inside=inside, mask=np.packbits(mask),
                        cell_xy=cells, cell_counts=counts, cell_area=area,
                        cell_type=region_of(cells[:, 0]))
    truth = {"note": "SYNTHETIC. Made by examples/make_stereoseq_demo.py; no tissue, chip or run is behind it.",
             "seed": args.seed, "sn": SN, "width": W, "height": H, "dnbUm": RES_NM / 1000, "geneIds": IDS,
             "geneNames": NAMES, "labels": NAMES[:-1] + [f"SynDup_{IDS[-1]}"], "oldLayoutNames": NAMES_UNIQUE,
             "types": TYPES, "bandEdges": list(BAND_EDGES), "island": list(ISLAND), "cut": list(CUT), "records": int(len(g)), "mid": int(c.sum()), "recordsInTissue": int(inside.sum()),
             "midInTissue": int(c[inside].sum()), "tissueDnb": int(mask.sum()), "nCells": int(len(cells)), "cellRadiusDnb": CELL_R,
             "cellMid": int(counts.sum())}
    (out / "truth_stereoseq.json").write_text(json.dumps(truth, indent=1), encoding="utf-8")
    print(f"synthetic Stereo-seq folder -> {folder}")
    print(f"  chip {W} x {H} DNB at {RES_NM} nm; tissue {int(mask.sum()):,} DNB; {len(g):,} records, {int(c.sum()):,} MID "
          f"({int(c[inside].sum()):,} in the tissue); {len(cells)} cells, {int(counts.sum()):,} MID in cells")


if __name__ == "__main__":
    main()
