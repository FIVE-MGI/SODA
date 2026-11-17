#!/usr/bin/env python
r"""Make small SYNTHETIC files in the standard formats SODA reads (stdformats.py) for the smoke test.

Nothing here comes from a tissue, an instrument or a run: every pixel, spot, cell and molecule is
drawn by a seeded random generator from a made-up layout. The files are written with h5py, zarr 2
and pyarrow directly, in the on-disk layouts anndata (0.8+), OME-NGFF 0.4 and spatialdata (0.1 /
shapes 0.2) write; anndata, spatialdata and ome-zarr are not needed.

    python examples/make_stdformats_demo.py --out demo/raw [--cases demo/stdformats_cases]

writes into --out (every one is a dataset detect.py must call certain):

    anndata_visium.h5ad     the Squidpy / scanpy.read_visium layout of a classic Visium section: 240 spots
                            on a hexagonal 100 µm grid (barcodes ACGT x 16 + "-1"), X a CSR matrix of
                            counts (20 genes + 2 Antibody Capture features + one BLANK control), var
                            gene_ids index + gene_symbols, obs: cluster (categorical, uns cluster_colors),
                            sample (categorical; its value is a made-up sample ID), in_tissue, total; obsm
                            spatial (full-resolution pixels) and X_umap; uns spatial/<library>/images/hires
                            (RGB uint8) + lowres (float) and scalefactors (hires 0.25, spot 89 px). The
                            library name holds a made-up sample ID too (it must not reach a message or the export)
    anndata_cells.zarr      an AnnData-Zarr of 300 cells (dense X, string index, a leiden column, obsm
                            spatial in µm), no image (prepare with --coords-um 1)
    image.ome.zarr          an OME-Zarr image, NGFF 0.4: c, y, x with micrometer units, 0.5 µm pixels,
                            four uint16 channels (DAPI, CD45, PanCK, CD3 with omero colours), three levels
    sdata_xenium.zarr       a SpatialData store shaped like spatialdata-io's Xenium reader: images/morphology
                            (2 uint16 channels, identity to global, also in an "aligned" system), shapes
                            cell_circles (radius from the cell size; scale µm -> pixels), cell_boundaries
                            (16-gons, scale), nucleus_boundaries (8-gons, an affine), points/transcripts
                            (two parquet parts, one with the gene as a dictionary column, qv, a sequence
                            transformation; two control codes), tables/table (counts of the cells, region
                            cell_circles, instance_key cell_id, cell_type categorical)
    sdata_visiumhd.zarr     a SpatialData store shaped like a Visium HD 16 µm layer: images/hires (RGB uint8,
                            a scale of 4 to global), shapes/square_016um (equal squares, a range index),
                            tables/square_016um
    truth_stdformats.npz / .json   what the smoke test checks against

and with --cases, folders / files SODA must refuse or ask about (never in --out, which detect.py
must list as certain datasets only):

    v3.ome.zarr             a zarr v3 store (zarr.json, OME-NGFF 0.5)
    plate.ome.zarr          an HCS plate
    two_libraries.h5ad      uns spatial holds two libraries (two sections in one object)
    labels_region.zarr      a SpatialData table annotating a labels element
    mapaxis.zarr            a SpatialData shapes element with a mapAxis transformation
"""
from __future__ import annotations

import argparse
import json
import struct
from pathlib import Path

import numpy as np

FAKE_ID = "INT" + "9961"            # a made-up sample ID (never a real one), built so this file holds no ID  (deid:ok)
LIBRARY = FAKE_ID + "_A1"           # the Visium library name: a sample ID, as Space Ranger run names often are
GENES = [f"SynA{i}" for i in range(1, 6)] + [f"SynB{i}" for i in range(1, 6)] + \
        [f"SynC{i}" for i in range(1, 6)] + [f"SynHk{i}" for i in range(1, 6)]
SPOT_PX = 89.0                      # spot_diameter_fullres
HIRES = 0.25                        # tissue_hires_scalef
LOWRES = 0.05
XE_PX = 0.5                         # µm per morphology pixel of the Xenium-like store
HD_IMG_SCALE = 4.0                  # global (full-resolution) pixels per hires pixel of the Visium HD store
HD_FULL_UM = 0.5                    # µm per global pixel there (so a 16 µm bin is 32 global / 8 hires pixels)


# ---------------------------------------------------------------------------------------
# AnnData writing (anndata 0.8+ encoding), the same code for h5py and zarr
# ---------------------------------------------------------------------------------------
class Writer:
    def __init__(self, root, kind):
        self.root, self.kind = root, kind

    def group(self, parent, name, enc=None, ver="0.1.0", **attrs):
        g = parent.create_group(name)
        if enc:
            g.attrs["encoding-type"] = enc
            g.attrs["encoding-version"] = ver
        for k, v in attrs.items():
            g.attrs[k] = v
        return g

    def strings(self, parent, name, values, enc=True):
        vals = [str(v) for v in values]
        if self.kind == "h5":
            import h5py
            d = parent.create_dataset(name, data=np.array(vals, dtype=object), dtype=h5py.string_dtype("utf-8"))
        else:
            import numcodecs
            d = parent.create_dataset(name, data=np.array(vals, dtype=object), dtype=object,
                                      object_codec=numcodecs.VLenUTF8(), shape=(len(vals),))
        if enc:
            d.attrs["encoding-type"] = "string-array"
            d.attrs["encoding-version"] = "0.2.0"
        return d

    def array(self, parent, name, a):
        a = np.asarray(a)
        if self.kind == "h5":
            d = parent.create_dataset(name, data=a)
        else:
            d = parent.create_dataset(name, data=a, shape=a.shape, dtype=a.dtype)
        d.attrs["encoding-type"] = "array"
        d.attrs["encoding-version"] = "0.2.0"
        return d

    def scalar(self, parent, name, v):
        if isinstance(v, str):
            if self.kind == "h5":
                import h5py
                d = parent.create_dataset(name, data=v, dtype=h5py.string_dtype("utf-8"))
            else:
                import numcodecs
                d = parent.create_dataset(name, data=np.array(v, dtype=object), dtype=object,
                                          object_codec=numcodecs.VLenUTF8(), shape=())
            d.attrs["encoding-type"] = "string"
        else:
            if self.kind == "h5":
                d = parent.create_dataset(name, data=np.float64(v))
            else:
                d = parent.create_dataset(name, data=np.array(v, np.float64), shape=())
            d.attrs["encoding-type"] = "numeric-scalar"
        d.attrs["encoding-version"] = "0.2.0"
        return d

    def categorical(self, parent, name, values, colors_into=None, palette=None):
        cats = sorted(set(values))
        g = self.group(parent, name, "categorical", "0.2.0", ordered=False)
        self.array(g, "codes", np.array([cats.index(v) for v in values], np.int8))
        self.strings(g, "categories", cats)
        if colors_into is not None and palette:
            self.strings(colors_into, f"{name}_colors", [palette[c] for c in cats])
        return g

    def dataframe(self, parent, name, index, columns: dict, uns=None, palettes=None):
        g = self.group(parent, name, "dataframe", "0.2.0", _index="_index", **{"column-order": list(columns)})
        self.strings(g, "_index", index)
        for col, v in columns.items():
            if isinstance(v, tuple) and v[0] == "cat":
                self.categorical(g, col, v[1], uns, (palettes or {}).get(col))
            elif isinstance(v, tuple) and v[0] == "str":
                self.strings(g, col, v[1])
            else:
                self.array(g, col, v)
        return g

    def csr(self, parent, name, X):
        from scipy import sparse
        X = sparse.csr_matrix(X)
        g = self.group(parent, name, "csr_matrix", "0.1.0", shape=list(X.shape))
        self.array(g, "data", X.data.astype(np.float32))
        self.array(g, "indices", X.indices.astype(np.int32))
        self.array(g, "indptr", X.indptr.astype(np.int64))
        return g

    def dict_(self, parent, name, d: dict):
        g = self.group(parent, name, "dict", "0.1.0")
        for k, v in d.items():
            if isinstance(v, dict):
                self.dict_(g, k, v)
            elif isinstance(v, np.ndarray):
                self.array(g, k, v)
            else:
                self.scalar(g, k, v)
        return g


def anndata(path: Path, kind: str, X, obs_names, var_names, obs: dict, var: dict, obsm: dict, uns: dict,
            palettes=None, dense=False):
    """An AnnData file / store: X (CSR, or dense), obs / var dataframes, obsm arrays, uns (nested dicts)."""
    if kind == "h5":
        import h5py
        f = h5py.File(str(path), "w")
    else:
        import zarr
        f = zarr.open_group(str(path), mode="w")
    f.attrs["encoding-type"] = "anndata"
    f.attrs["encoding-version"] = "0.1.0"
    w = Writer(f, kind)
    ug = w.group(f, "uns", "dict")
    w.dataframe(f, "obs", obs_names, obs, ug, palettes)
    w.dataframe(f, "var", var_names, var)
    if dense:
        w.array(f, "X", np.asarray(X, np.float32))
    else:
        w.csr(f, "X", X)
    og = w.group(f, "obsm", "dict")
    for k, v in obsm.items():
        w.array(og, k, v)
    for k, v in uns.items():
        if isinstance(v, dict):
            w.dict_(ug, k, v)
        else:
            w.scalar(ug, k, v)
    w.group(f, "layers", "dict")
    w.group(f, "varm", "dict")
    w.group(f, "obsp", "dict")
    w.group(f, "varp", "dict")
    if kind == "h5":
        f.close()


# ---------------------------------------------------------------------------------------
# the layouts
# ---------------------------------------------------------------------------------------
def counts_for(types: np.ndarray, rng, n_genes=20, base=6.0):
    """Counts with a block per type (A / B / C genes) and shared housekeeping genes."""
    n = len(types)
    lam = np.full((n, n_genes), 0.3)
    for t in range(3):
        lam[types == t, t * 5:(t + 1) * 5] = base
    lam[:, 15:20] = base / 2
    return rng.poisson(lam).astype(np.float32)


def write_anndata_visium(out: Path, rng, truth: dict, two_libraries: bool = False):
    from scipy import sparse
    pts = []
    for r in range(12):
        for q in range(20):
            pts.append((400 + q * 160 + (80 if r % 2 else 0), 400 + r * 139))      # 100 µm grid, 160 px apart
    xy = np.array(pts, np.float64)
    n = len(xy)
    types = np.minimum((xy[:, 0] - 400) // 1100, 2).astype(int)
    X = counts_for(types, rng)
    extra = np.column_stack([rng.poisson(3.0, n), rng.poisson(3.0, n), rng.poisson(1.0, n)]).astype(np.float32)
    Xall = sparse.csr_matrix(np.hstack([X, extra]))
    letters = np.array(list("ACGT"))
    bcs = ["".join(letters[rng.integers(0, 4, 16)]) + "-1" for _ in range(n)]
    gene_ids = [f"SYNG{i:05d}" for i in range(len(GENES) + 3)]
    symbols = GENES + ["AbCd45", "AbCd3", "BLANK_0001"]
    ftypes = ["Gene Expression"] * len(GENES) + ["Antibody Capture"] * 2 + ["Gene Expression"]
    names = ["Region A", "Region B", "Region C"]
    cluster = [names[t] for t in types]
    palette = {"Region A": "#1f77b4", "Region B": "#ff7f0e", "Region C": "#2ca02c"}
    W, H = 4000, 2400
    hires = np.full((int(H * HIRES), int(W * HIRES), 3), 235, np.uint8)
    for (x, y), t in zip(xy, types):                       # a pink dot under each spot, darker by type
        cx, cy = int(x * HIRES), int(y * HIRES)
        hires[max(cy - 8, 0):cy + 9, max(cx - 8, 0):cx + 9] = (200 - 40 * t, 120, 180)
    lowres = (hires[::5, ::5].astype(np.float32) / 255.0)
    sf = {"tissue_hires_scalef": HIRES, "tissue_lowres_scalef": LOWRES, "spot_diameter_fullres": SPOT_PX,
          "fiducial_diameter_fullres": 134.0}
    lib = {"images": {"hires": hires, "lowres": lowres}, "scalefactors": sf,
           "metadata": {"chemistry_description": "Spatial 3' v1", "software_version": "spaceranger-2.1.0"}}
    spatial = {LIBRARY: lib}
    obs = {"in_tissue": np.ones(n, np.int64), "cluster": ("cat", cluster), "sample": ("cat", [f"{FAKE_ID}_A{1 + (i >= n // 2)}" for i in range(n)]),
           "total": np.asarray(Xall.sum(1)).ravel().astype(np.float32)}
    if two_libraries:
        other = FAKE_ID + "_B1"
        spatial[other] = lib
        obs["library_id"] = ("cat", [LIBRARY if i < n // 2 else other for i in range(n)])
    umap = np.column_stack([types * 4.0 + rng.normal(0, 0.5, n), rng.normal(0, 0.5, n)])
    path = out / ("two_libraries.h5ad" if two_libraries else "anndata_visium.h5ad")
    anndata(path, "h5", Xall, bcs, gene_ids, obs,
            {"gene_symbols": ("str", symbols), "feature_types": ("cat", ftypes)},
            {"spatial": xy, "X_umap": umap}, {"spatial": spatial}, palettes={"cluster": palette})
    if not two_libraries:
        truth["av_xy"] = xy
        truth["av_counts"] = X
        truth["av_types"] = types
        truth["av_bcs"] = np.array(bcs)


def write_anndata_cells(out: Path, rng, truth: dict):
    n = 300
    xy = np.column_stack([rng.uniform(20, 380, n), rng.uniform(20, 280, n)])        # µm
    types = np.minimum((xy[:, 0] - 20) // 120, 2).astype(int)
    X = counts_for(types, rng, base=4.0)
    ids = [f"cell{i:04d}" for i in range(n)]
    anndata(out / "anndata_cells.zarr", "zarr", X, ids, GENES, {"leiden": ("cat", [str(t) for t in types])},
            {"highly_variable": np.ones(len(GENES), bool)}, {"spatial": xy}, {}, dense=True)
    truth["ac_xy"] = xy
    truth["ac_counts"] = X
    truth["ac_types"] = types


# ---------------------------------------------------------------------------------------
# OME-Zarr (NGFF 0.4) and SpatialData
# ---------------------------------------------------------------------------------------
def _axes(names, units=None):
    t = {"c": "channel", "t": "time", "z": "space", "y": "space", "x": "space"}
    out = []
    for n in names:
        a = {"name": n, "type": t[n]}
        if units and t[n] == "space":
            a["unit"] = units
        out.append(a)
    return out


def ngff_image(path: Path, data: np.ndarray, names, colors, px_um=None, cs=None, sdata=False):
    """A multiscale c, y, x image (levels halved twice). px_um: the level-0 scale with micrometer units.
    cs: SpatialData {coordinate system: transformation dict} put on the multiscales."""
    import zarr
    g = zarr.open_group(str(path), mode="w")
    levels = [data]
    for _ in range(2):
        a = levels[-1]
        levels.append(a[:, ::2, ::2].copy())
    datasets = []
    for i, a in enumerate(levels):
        g.create_dataset(str(i), data=a, shape=a.shape, dtype=a.dtype, chunks=(1, 128, 128), dimension_separator="/")
        s = (px_um or 1.0) * 2 ** i
        datasets.append({"path": str(i), "coordinateTransformations": [{"type": "scale", "scale": [1.0, s, s]}]})
    ms = {"version": "0.4", "name": path.stem, "axes": _axes("cyx", "micrometer" if px_um else None), "datasets": datasets}
    if cs:
        ms["coordinateTransformations"] = [{**t, "input": {"name": "cyx", "axes": _axes("cyx")},
                                            "output": {"name": name, "axes": _axes("cyx")}} for name, t in cs.items()]
    attrs = {"multiscales": [ms],
             "omero": {"name": path.stem, "version": "0.4",
                       "channels": [{"label": n, "color": c, "active": True, "window": {"start": 0, "end": 4000, "min": 0, "max": 65535}}
                                    for n, c in zip(names, colors)]}}
    if sdata:
        attrs["spatialdata_attrs"] = {"version": "0.1"}
    g.attrs.update(attrs)


def _wkb_point(x, y) -> bytes:
    return struct.pack("<BIdd", 1, 1, float(x), float(y))


def _wkb_polygon(ring) -> bytes:
    ring = list(ring) + [ring[0]]
    return struct.pack("<BII", 1, 3, 1) + struct.pack("<I", len(ring)) + b"".join(struct.pack("<dd", *p) for p in ring)


def _xy_axes():
    return [{"name": "x", "type": "space", "unit": "unit"}, {"name": "y", "type": "space", "unit": "unit"}]


def sd_transforms(ts: dict) -> list:
    return [{**t, "input": {"name": "xy", "axes": _xy_axes()}, "output": {"name": name, "axes": _xy_axes()}}
            for name, t in ts.items()]


def sd_shapes(folder: Path, geoms: list, ids, transforms: dict, radius=None, range_index=False):
    """shapes/<name>: shapes.parquet (GeoParquet: WKB geometry, radius, the index) + .zattrs."""
    import pyarrow as pa
    import pyarrow.parquet as pq
    import zarr
    zarr.open_group(str(folder), mode="w")
    cols = {"geometry": pa.array(geoms, pa.binary())}
    if radius is not None:
        cols["radius"] = pa.array(np.asarray(radius, np.float64))
    if range_index:
        index = [{"kind": "range", "name": None, "start": 0, "stop": len(geoms), "step": 1}]
    else:
        cols["cell_id"] = pa.array([str(i) for i in ids])
        index = ["cell_id"]
    t = pa.table(cols)
    geo = {"version": "1.0.0", "primary_column": "geometry",
           "columns": {"geometry": {"encoding": "WKB", "geometry_types": ["Point" if radius is not None else "Polygon"]}}}
    pandas = {"index_columns": index, "column_indexes": [], "columns": [], "pandas_version": "2.2.0"}
    t = t.replace_schema_metadata({b"geo": json.dumps(geo).encode(), b"pandas": json.dumps(pandas).encode()})
    pq.write_table(t, str(folder / "shapes.parquet"))
    zarr.open_group(str(folder), mode="r+").attrs.update({
        "axes": ["x", "y"], "coordinateTransformations": sd_transforms(transforms), "encoding-type": "ngff:shapes",
        "spatialdata_attrs": {"version": "0.2"}})


def ring(cx, cy, r, k, phase=0.0):
    a = phase + np.arange(k) * 2 * np.pi / k
    return list(zip((cx + r * np.cos(a)).tolist(), (cy + r * np.sin(a)).tolist()))


def write_sdata_xenium(out: Path, rng, truth: dict):
    import pyarrow as pa
    import pyarrow.parquet as pq
    import zarr
    root = out / "sdata_xenium.zarr"
    g = zarr.open_group(str(root), mode="w")
    g.attrs["spatialdata_attrs"] = {"version": "0.1"}
    W, H = 400, 300                                       # morphology pixels (0.5 µm): 200 x 150 µm
    n = 80
    centers = np.column_stack([rng.uniform(15, 185, n), rng.uniform(15, 135, n)])     # µm
    rad = rng.uniform(4.0, 7.0, n)                        # µm
    types = np.minimum((centers[:, 0] - 15) // 57, 2).astype(int)
    ids = ["".join(rng.choice(list("abcdefghijklmnop"), 8)) + "-1" for _ in range(n)]
    yy, xx = np.mgrid[0:H, 0:W]
    dapi = np.zeros((H, W), np.float64)
    for (cx, cy), r in zip(centers / XE_PX, rad / XE_PX):
        dapi += 3000 * ((xx - cx) ** 2 + (yy - cy) ** 2 < (0.5 * r) ** 2)
    img = np.stack([dapi, 0.5 * dapi + 200]).astype(np.uint16)
    g.create_group("images")
    ngff_image(root / "images" / "morphology", img, ["DAPI", "Boundary"], ["0000FF", "00FF00"], sdata=True,
               cs={"global": {"type": "identity"}, "aligned": {"type": "translation", "translation": [0.0, 5.0, 5.0]}})
    g.create_group("shapes")
    s = 1.0 / XE_PX
    sd_shapes(root / "shapes" / "cell_circles", [_wkb_point(x, y) for x, y in centers], ids,
              {"global": {"type": "scale", "scale": [s, s]}, "aligned": {"type": "scale", "scale": [s, s]}}, radius=rad)
    sd_shapes(root / "shapes" / "cell_boundaries", [_wkb_polygon(ring(x, y, r, 16)) for (x, y), r in zip(centers, rad)], ids,
              {"global": {"type": "scale", "scale": [s, s]}})
    sd_shapes(root / "shapes" / "nucleus_boundaries", [_wkb_polygon(ring(x, y, 0.5 * r, 8, 0.2)) for (x, y), r in zip(centers, rad)],
              ids, {"global": {"type": "affine", "affine": [[s, 0.0, 0.0], [0.0, s, 0.0], [0.0, 0.0, 1.0]]}})
    # transcripts: per cell its counts, within the cell; a few controls and low-quality ones
    X = counts_for(types, rng, base=5.0)
    gx, gy, gn, qv = [], [], [], []
    for i in range(n):
        for j in range(len(GENES)):
            k = int(X[i, j])
            if not k:
                continue
            a = rng.uniform(0, 2 * np.pi, k)
            d = rng.uniform(0, rad[i] * 0.9, k)
            gx += (centers[i, 0] + d * np.cos(a)).tolist()
            gy += (centers[i, 1] + d * np.sin(a)).tolist()
            gn += [GENES[j]] * k
            qv += [30.0] * k
    n_good = len(gx)
    for name in ("NegControlProbe_00001", "BLANK_0002"):
        gx += rng.uniform(10, 190, 15).tolist()
        gy += rng.uniform(10, 140, 15).tolist()
        gn += [name] * 15
        qv += [30.0] * 15
    low = 40
    gx += rng.uniform(10, 190, low).tolist()
    gy += rng.uniform(10, 140, low).tolist()
    gn += [GENES[0]] * low
    qv += [10.0] * low
    order = rng.permutation(len(gx))
    gx, gy, gn, qv = (np.asarray(v)[order] for v in (gx, gy, gn, qv))
    pdir = root / "points" / "transcripts"
    zarr.open_group(str(root / "points"), mode="w")
    zarr.open_group(str(pdir), mode="w")
    (pdir / "points.parquet").mkdir()
    half = len(gx) // 2
    for part, sl, dictionary in ((0, slice(0, half), True), (1, slice(half, None), False)):
        feat = pa.array(gn[sl].tolist())
        if dictionary:
            feat = feat.dictionary_encode()
        t = pa.table({"x": pa.array(gx[sl]), "y": pa.array(gy[sl]), "z": pa.array(np.full(len(gx[sl]), 5.0)),
                      "feature_name": feat, "qv": pa.array(qv[sl])})
        pq.write_table(t, str(pdir / "points.parquet" / f"part.{part}.parquet"))
    zarr.open_group(str(pdir), mode="r+").attrs.update({
        "axes": ["x", "y"], "encoding-type": "ngff:points",
        "coordinateTransformations": sd_transforms({"global": {"type": "sequence", "transformations": [
            {"type": "scale", "scale": [s, s]}, {"type": "translation", "translation": [0.0, 0.0]}]}}),
        "spatialdata_attrs": {"feature_key": "feature_name", "instance_key": "cell_id", "version": "0.1"}})
    # the table: rows in another order than the shapes
    perm = rng.permutation(n)
    zarr.open_group(str(root / "tables"), mode="w")
    names = ["Type A", "Type B", "Type C"]
    anndata(root / "tables" / "table", "zarr", X[perm], [str(i) for i in range(n)], GENES,
            {"cell_id": ("str", [ids[i] for i in perm]), "region": ("cat", ["cell_circles"] * n),
             "cell_type": ("cat", [names[types[i]] for i in perm]), "transcript_counts": X[perm].sum(1)},
            {"feature_types": ("cat", ["Gene Expression"] * len(GENES))}, {"spatial": centers[perm]},
            {"spatialdata_attrs": {"region": "cell_circles", "region_key": "region", "instance_key": "cell_id"}})
    truth["sx_xy_px"] = centers[perm] / XE_PX
    truth["sx_counts"] = X[perm]
    truth["sx_ids"] = np.array([ids[i] for i in perm])
    truth["sx_types"] = types[perm]
    truth["sx_rad_px"] = rad[perm] / XE_PX
    truth["sx_molecules"] = np.array([n_good])
    truth["sx_mol_per_gene"] = np.array([int((gn[qv >= 20] == gname).sum()) for gname in GENES])


def write_sdata_visiumhd(out: Path, rng, truth: dict):
    import zarr
    from scipy import sparse
    root = out / "sdata_visiumhd.zarr"
    g = zarr.open_group(str(root), mode="w")
    g.attrs["spatialdata_attrs"] = {"version": "0.1"}
    side = 16.0 / HD_FULL_UM                              # global pixels per bin
    nx, ny = 14, 10
    x0, y0 = 40.0, 48.0
    cells = [(i, j) for j in range(ny) for i in range(nx)]
    n = len(cells)
    geoms, centres = [], []
    for i, j in cells:
        x, y = x0 + i * side, y0 + j * side
        geoms.append(_wkb_polygon([(x, y), (x + side, y), (x + side, y + side), (x, y + side)]))
        centres.append((x + side / 2, y + side / 2))
    centres = np.array(centres)
    types = np.array([min(i * 3 // nx, 2) for i, _ in cells])
    X = counts_for(types, rng, base=8.0)
    Wg, Hg = 560, 448                                     # global (full-resolution) extent
    hi = np.full((int(Hg / HD_IMG_SCALE), int(Wg / HD_IMG_SCALE), 3), 230, np.uint8)
    for (cx, cy), t in zip(centres / HD_IMG_SCALE, types):
        hi[int(cy) - 2:int(cy) + 3, int(cx) - 2:int(cx) + 3] = (180 - 50 * t, 100, 170)
    g.create_group("images")
    ngff_image(root / "images" / "hires_image", np.moveaxis(hi, 2, 0).copy(), ["r", "g", "b"], ["FF0000", "00FF00", "0000FF"],
               sdata=True, cs={"global": {"type": "scale", "scale": [1.0, HD_IMG_SCALE, HD_IMG_SCALE]}})
    g.create_group("shapes")
    sd_shapes(root / "shapes" / "square_016um", geoms, None, {"global": {"type": "identity"}}, range_index=True)
    zarr.open_group(str(root / "tables"), mode="w")
    bcs = [f"s_016um_{j:05d}_{i:05d}-1" for i, j in cells]
    perm = rng.permutation(n)
    anndata(root / "tables" / "square_016um", "zarr", sparse.csr_matrix(X[perm]), [bcs[k] for k in perm], GENES,
            {"location_id": np.asarray(perm, np.int64), "region": ("cat", ["square_016um"] * n)},
            {}, {}, {"spatialdata_attrs": {"region": "square_016um", "region_key": "region", "instance_key": "location_id"}})
    truth["sh_xy_px"] = centres[perm] / HD_IMG_SCALE
    truth["sh_counts"] = X[perm]
    truth["sh_bcs"] = np.array([bcs[k] for k in perm])


def write_omezarr(out: Path, rng, truth: dict):
    H, W = 320, 480
    yy, xx = np.mgrid[0:H, 0:W]
    dapi = (800 + 3000 * (((xx // 40) + (yy // 40)) % 2)).astype(np.uint16)
    cd45 = (100 + 20 * xx // 4).astype(np.uint16)
    panck = (100 + 20 * yy // 4).astype(np.uint16)
    cd3 = (100 + 2500 * ((xx - 240) ** 2 + (yy - 160) ** 2 < 90 ** 2)).astype(np.uint16)
    img = np.stack([dapi, cd45, panck, cd3])
    ngff_image(out / "image.ome.zarr", img, ["DAPI", "CD45", "PanCK", "CD3"], ["0000FF", "FF0000", "00FF00", "FFFF00"],
               px_um=0.5)
    truth["oz_image"] = img


# ---------------------------------------------------------------------------------------
# cases SODA must refuse or ask about
# ---------------------------------------------------------------------------------------
def write_cases(cases: Path, rng):
    import zarr
    cases.mkdir(parents=True, exist_ok=True)
    v3 = cases / "v3.ome.zarr"
    v3.mkdir(exist_ok=True)
    (v3 / "zarr.json").write_text(json.dumps({"zarr_format": 3, "node_type": "group", "attributes": {"ome": {
        "version": "0.5", "multiscales": [{"axes": _axes("cyx"), "datasets": [{"path": "0", "coordinateTransformations": [
            {"type": "scale", "scale": [1, 1, 1]}]}]}]}}}), encoding="utf-8")
    p = zarr.open_group(str(cases / "plate.ome.zarr"), mode="w")
    p.attrs["plate"] = {"version": "0.4", "rows": [{"name": "A"}], "columns": [{"name": "1"}], "wells": [{"path": "A/1", "rowIndex": 0, "columnIndex": 0}]}
    write_anndata_visium(cases, rng, {}, two_libraries=True)
    # a table on a labels element
    root = cases / "labels_region.zarr"
    g = zarr.open_group(str(root), mode="w")
    g.attrs["spatialdata_attrs"] = {"version": "0.1"}
    zarr.open_group(str(root / "labels"), mode="w")
    lab = np.zeros((1, 64, 64), np.uint16)
    lab[0, 10:20, 10:20] = 1
    ngff_image(root / "labels" / "cells", lab, ["labels"], ["FFFFFF"], sdata=True, cs={"global": {"type": "identity"}})
    zarr.open_group(str(root / "labels" / "cells"), mode="r+").attrs["image-label"] = {"version": "0.4"}
    zarr.open_group(str(root / "tables"), mode="w")
    anndata(root / "tables" / "table", "zarr", np.ones((1, 3), np.float32), ["0"], GENES[:3],
            {"instance_id": np.array([1], np.int64), "region": ("cat", ["cells"])}, {}, {},
            {"spatialdata_attrs": {"region": "cells", "region_key": "region", "instance_key": "instance_id"}})
    # a shapes element with a mapAxis transformation
    root = cases / "mapaxis.zarr"
    g = zarr.open_group(str(root), mode="w")
    g.attrs["spatialdata_attrs"] = {"version": "0.1"}
    zarr.open_group(str(root / "shapes"), mode="w")
    sd_shapes(root / "shapes" / "spots", [_wkb_point(10, 10), _wkb_point(20, 20)], ["a", "b"],
              {"global": {"type": "mapAxis", "mapAxis": {"x": "y", "y": "x"}}}, radius=[2.0, 2.0])
    zarr.open_group(str(root / "tables"), mode="w")
    anndata(root / "tables" / "table", "zarr", np.ones((2, 3), np.float32), ["0", "1"], GENES[:3],
            {"spot": ("str", ["a", "b"]), "region": ("cat", ["spots", "spots"])}, {}, {},
            {"spatialdata_attrs": {"region": "spots", "region_key": "region", "instance_key": "spot"}})


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--cases", help="also write the cases SODA must refuse or ask about into this folder")
    ap.add_argument("--seed", type=int, default=11)
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(a.seed)
    truth: dict = {}
    write_anndata_visium(out, rng, truth)
    write_anndata_cells(out, rng, truth)
    write_omezarr(out, rng, truth)
    write_sdata_xenium(out, rng, truth)
    write_sdata_visiumhd(out, rng, truth)
    np.savez_compressed(out / "truth_stdformats.npz", **truth)
    (out / "truth_stdformats.json").write_text(json.dumps({
        "genes": GENES, "spotPx": SPOT_PX, "hires": HIRES, "xeniumPx": XE_PX, "hdImageScale": HD_IMG_SCALE,
        "hdFullUm": HD_FULL_UM, "library": "[made-up sample ID]_A1"}, indent=1), encoding="utf-8")
    if a.cases:
        write_cases(Path(a.cases), np.random.default_rng(a.seed + 1))
    print(f"wrote the standard-format fixtures to {out}" + (f" and the cases to {a.cases}" if a.cases else ""))


if __name__ == "__main__":
    main()
