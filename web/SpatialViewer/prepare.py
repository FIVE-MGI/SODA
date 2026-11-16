r"""Build the viewer cache for one dataset.

    python prepare.py akoya  --src <folder of NN_Marker_Sample_gray.tif> --id sample1 --out <cache folder>
    python prepare.py tiff   --src <stack.tif|qptiff|ome.tif>            --id x        --out <cache folder>
    python prepare.py visium --image <hires.tif> --outs <spaceranger outs> --id sample1_visiumhd --out <cache folder>
    python prepare.py xenium --src <xenium output folder> --id x1 --out <cache folder>
           (stains; cells -> a "cells" expression layer; transcripts -> the molecules layer)
    python prepare.py geomx  --dcc <folder of .dcc> --pkc <panel.pkc ...> --annotation <xlsx|csv|LabWorksheet.txt>
           [--image <scan.ome.tiff> ...] --id g1 --out <cache folder>
           (GeoMx DSP: one row per ROI segment in a "rois" layer, ROI outlines, annotation columns as groups)
    python prepare.py cosmx    --src <CosMx flat-file folder>  --id c1 --out <cache folder>
    python prepare.py merscope --src <MERSCOPE region folder>  --id m1 --out <cache folder>
           (FOV mosaic / mosaic stains, cells and molecules as for xenium; cosmx.py, merscope.py)
    python prepare.py stereoseq --image <SN>_HE_regist.tif --gef <SN>.raw.gef --cellbin <SN>.cellbin.gef --id s1 --out <cache folder>
           (Stereo-seq / SAW: bin20 / bin50 / bin100 layers cut to the tissue + a "cells" layer; --gem instead of --gef works too)

Common options: --label, --pixel-size (µm/px), --tile, --geojson "Label=path" (repeatable),
--rgb-split (keep RGB files as three channels instead of max-projecting).
The result is <out>/<id>/ with image.zarr, expr/, overlays/ and manifest.json;
the server just scans <out> for manifest.json files.

    python prepare.py export --id sample1 --out <cache folder> --dest <static folder> [--fmt jpg]
           [--user-layer 008um --user-umap umap.csv --user-obs clusters.csv | --user-h5ad analysis.h5ad]
    python prepare.py userclusters --id sample1 --out <cache folder> --dest <static folder> --user-layer 008um ...

writes the static tree. For an image dataset (akoya / tiff / xenium) it also adds a 16 µm grid
of marker means (from full-resolution pixels) as a "protein" expression layer with a QC report,
a UMAP, Leiden clusters and suggested names (--grid-um N to change the square size, 0 to skip;
--markers picks the marker list, default markers/protein_generic.json for protein and
markers/rna_generic.json for Visium genes; tissue-specific lists such as markers/nerve.json,
markers/brain.json and markers/akoya_nerve.json are passed with --markers). The --user-* options attach the lab's own
UMAP and clusterings (Seurat / scanpy, matched by barcode) to one spatial layer (userclusters.py);
`userclusters` does the same on an existing export without re-running it.
"""
from __future__ import annotations

import argparse
import json
import re
import shutil
import sys
import time
from datetime import datetime
from pathlib import Path

import threadcap                                   # first: caps the BLAS / OpenMP threads before numpy loads
threadcap.apply(threadcap.from_argv(sys.argv))     # --threads N, else SODA_THREADS, else min(8, usable CPUs)

import pyramid
import stereoseq
import userclusters
from expr import build_cell_layer, build_expr, find_visium_layers


def log(*a):
    print(*a, flush=True)


VISIUM_DEFAULT_BINS = ("008um", "016um")   # Visium HD layers kept without --bins; --bins all keeps every one


def select_visium_layers(layers: list[dict], bins: str | None) -> list[dict]:
    """The Space Ranger layers `prepare.py visium` builds, logging what is kept and left out.
    bins None/"": Visium HD 8 and 16 µm (all HD layers if neither exists), classic spots always;
    "all": every layer; otherwise a comma list of layer keys, taken literally."""
    present = [L["key"] for L in layers]
    want = (bins or "").strip()
    if want.lower() == "all":
        keep, how = set(present), "--bins all"
    elif want:
        keep, how = {b.strip() for b in want.split(",") if b.strip()}, "--bins " + want
        unknown = sorted(keep - set(present))
        if unknown:
            log(f"[visium] WARNING --bins: no layer {', '.join(unknown)} under outs (present: "
                f"{', '.join(present) or 'none'})")
    else:
        hd = [L["key"] for L in layers if L["binSizeUm"] is not None]
        hit = set(VISIUM_DEFAULT_BINS) & set(hd)
        keep = (hit or set(hd)) | {L["key"] for L in layers if L["binSizeUm"] is None}   # + classic spots
        how = "default" if hit or not hd else "default, but no 8 or 16 µm layer, so all of them"
    chosen = [L for L in layers if L["key"] in keep]
    left = [k for k in present if k not in keep]
    log(f"[visium] layers kept ({how}): {', '.join(L['key'] for L in chosen) or 'none'}")
    if left:
        log(f"[visium] left out: {', '.join(left)}. To build them too: --bins all, or list the ones you want, "
            f"e.g. --bins {','.join(sorted(set(left) | {L['key'] for L in chosen}))}")
    return chosen


def akoya_files(folder: Path):
    """Marker files in an Akoya export folder, ordered by their leading index.
    `NN_Marker_<Sample>_gray.tif` is a plain grey channel; `_blue` etc. are
    RGB renderings and only used when no grey twin exists."""
    folder = Path(folder)
    tifs = [p for p in folder.iterdir() if p.suffix.lower() in (".tif", ".tiff")]
    sample = folder.name
    entries = {}
    for p in tifs:
        stem = p.stem
        m = re.match(r"^(.*)_(gray|grey|blue|red|green|rgb)$", stem, re.I)
        variant = m.group(2).lower() if m else "gray"
        base = m.group(1) if m else stem
        if base.endswith("_" + sample):
            base = base[: -len(sample) - 1]
        m2 = re.match(r"^(\d+)[_ ]+(.*)$", base)
        order = int(m2.group(1)) if m2 else 10 ** 6
        name = (m2.group(2) if m2 else base).strip() or stem
        cur = entries.get(name)
        pref = 0 if variant in ("gray", "grey") else 1
        if cur is None or pref < cur[0]:
            entries[name] = (pref, order, p)
    return [(n, e[2]) for n, e in sorted(entries.items(), key=lambda kv: (kv[1][1], kv[0]))]


def parse_geojson_args(items, ds_dir: Path):
    out = []
    for i, it in enumerate(items or []):
        if "=" in it:
            label, path = it.split("=", 1)
        else:
            label, path = Path(it).stem, it
        src = Path(path)
        if not src.exists():
            log(f"[overlay] missing: {src}")
            continue
        dst = ds_dir / "overlays" / f"{i:02d}_{re.sub(r'[^A-Za-z0-9_.-]+', '_', label)}.geojson"
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src, dst)
        out.append({"key": dst.stem, "label": label, "file": dst.name})
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="kind", required=True)

    def common(p):
        p.add_argument("--id", required=True, help="dataset id (folder name in the cache, URL-safe)")
        p.add_argument("--out", required=True, help="cache root folder")
        p.add_argument("--label")
        p.add_argument("--pixel-size", type=float, help="µm per full-res pixel (scale bar)")
        p.add_argument("--tile", type=int, default=pyramid.TILE)
        p.add_argument("--geojson", action="append", help='"Label=path.geojson" overlay in image pixel coords')
        p.add_argument("--rgb-split", action="store_true", help="RGB files -> three channels")

    a = sub.add_parser("akoya", help="folder of per-marker TIFFs")
    a.add_argument("--src", required=True)
    a.add_argument("--flip", help="h, v or hv - mirror the image")
    a.add_argument("--crop", help="x0,y0,x1,y1 in full-res pixels: build only that rectangle. Run once per region with a different --id to split a whole slide into samples.")
    common(a)
    t = sub.add_parser("tiff", help="one multi-channel TIFF / qptiff / OME-TIFF")
    t.add_argument("--src", required=True)
    t.add_argument("--channels", help="comma-separated channel names (override)")
    t.add_argument("--flip", help="h, v or hv - mirror the image")
    t.add_argument("--crop", help="x0,y0,x1,y1 in full-res pixels: build only that rectangle. Run once per region with a different --id to split a whole slide into samples.")
    common(t)
    e = sub.add_parser("export", help="write a prepared dataset as static files for blob/web hosting (no server)")
    e.add_argument("--id", required=True)
    e.add_argument("--out", required=True, help="cache root that prepare wrote into")
    e.add_argument("--dest", required=True, help="folder to write <dest>/<id>/ into (then upload it)")
    e.add_argument("--fmt", choices=["png", "jpg"], default="png")
    e.add_argument("--workers", type=int, default=8)
    e.add_argument("--markers", help="marker-list JSON for suggested cluster names. Default by what the dataset "
                                     "holds: markers/rna_generic.json (genes) for Visium / RNA layers, "
                                     "markers/protein_generic.json (\"modality\": \"protein\") for the grid layer "
                                     "of an image dataset. Tissue-specific lists: markers/nerve.json, "
                                     "markers/brain.json, markers/akoya_nerve.json. A list of the other modality is skipped with a message.")
    e.add_argument("--grid-um", type=float, default=16.0,
                   help="image datasets (Akoya / multiplex TIFF, not Visium or H&E): size in µm of the grid squares "
                        "whose marker means become a 'protein' expression layer with a UMAP and clusters "
                        "(protgrid.py). Default 16 (as the Visium HD 16 µm bins); 8 is about 4x the squares and "
                        "files. 0 = no grid layer. Needs the dataset's pixel size.")
    e.add_argument("--tissue-smooth-um", type=float, default=100.0,
                   help="grid layer: the tissue mask thresholds the DAPI square means after a moving average "
                        "over this many µm (default 100 = 7x7 squares of 16 µm); 0 = threshold the raw means "
                        "(keeps only nucleus-dense squares).")
    e.add_argument("--seam-correct", action="store_true",
                   help="grid layer: where camera-tile seams are found, correct the square means for "
                        "brightness / background differences between tiles (seams.py). Seams are always "
                        "detected and reported in the slide QC; shading shared by every tile is not corrected.")
    e.add_argument("--umap-sample", type=int, default=None,
                   help="grid layer: above this many tissue squares the section UMAP is fitted on a seeded random "
                        "sample of this size and the other squares are placed into it (default 200000; 0 = always "
                        "fit on every square, slower on big slides).")
    e.add_argument("--no-celltypes", action="store_true", help="skip the suggested cell types")
    e.add_argument("--no-umap", action="store_true", help="skip the whole-section UMAP + Leiden clusters (and so the cell types)")
    e.add_argument("--no-hvg", action="store_true", help="skip the gene chunks the circle -> UMAP tool reads")
    userclusters.add_args(e)
    uc = sub.add_parser("userclusters", help="add your own UMAP / clusters (Seurat, scanpy) to a layer of an existing export")
    uc.add_argument("--id", required=True)
    uc.add_argument("--out", required=True, help="cache root that prepare wrote into")
    uc.add_argument("--dest", required=True, help="export root holding <dest>/<id>/")
    uc.add_argument("--markers", help="marker list for the suggested names (default: the one the export used)")
    uc.add_argument("--no-celltypes", action="store_true", help="skip the suggested cell types")
    userclusters.add_args(uc)
    r = sub.add_parser("rnaseq", help="single-cell / nucleus RNA-seq -> static UMAP dataset (written straight to --dest)")
    r.add_argument("--id", required=True)
    r.add_argument("--label")
    r.add_argument("--dest", required=True, help="static export root")
    r.add_argument("--h5ad", help="AnnData file with an embedding in obsm")
    r.add_argument("--embedding", default="X_umap")
    r.add_argument("--mtx", help="10x matrix folder (matrix.mtx, features.tsv, barcodes.tsv)")
    r.add_argument("--umap", help="CSV: barcode, UMAP1, UMAP2 (with --mtx)")
    r.add_argument("--obs", action="append", default=[], help="CSV: barcode, column... (repeatable, with --mtx)")
    r.add_argument("--tmp", help="local folder for the temporary gene-major arrays of big files (needs ~8 bytes per non-zero)")
    i = sub.add_parser("inspect", help="print channels / levels / dtype of a TIFF or qptiff without building")
    i.add_argument("--src", required=True)
    x = sub.add_parser("xenium", help="Xenium output folder: morphology_focus stains, cells (cells.csv.gz + "
                                      "cell_feature_matrix) as a 'cells' expression layer, transcripts as the "
                                      "molecules layer (+ Custom_annotation.geojson if present)")
    x.add_argument("--src", required=True, help="the Xenium output folder (the one holding morphology_focus/)")
    x.add_argument("--no-cells", action="store_true", help="skip the cells layer")
    x.add_argument("--no-transcripts", action="store_true", help="skip the molecules layer")
    x.add_argument("--qv-min", type=float, default=20.0,
                   help="molecules: keep transcripts with quality value >= this (default 20, Xenium's own cut)")
    common(x)
    gx = sub.add_parser("geomx", help="GeoMx DSP: DCC counts + PKC probe kit(s) + annotation sheet (+ scan OME-TIFFs "
                                      "with the ROI polygons) -> a 'rois' layer, one row per ROI segment (geomx.py)")
    import geomx as _geomx
    _geomx.add_args(gx)
    common(gx)
    cx = sub.add_parser("cosmx", help="CosMx SMI flat files: per-FOV images stitched into one mosaic, cells "
                                      "(exprMat + metadata) as a 'cells' expression layer, tx_file as the "
                                      "molecules layer (cosmx.py)")
    cx.add_argument("--src", required=True, help="folder holding *_exprMat_file.csv, *_metadata_file.csv, "
                                                 "*_fov_positions_file.csv, *_tx_file.csv and the image folders")
    cx.add_argument("--images", choices=["auto", "morphology", "composite", "none"], default="auto",
                    help="auto: raw morphology TIFFs when every FOV has one, else CellComposite, else a "
                         "transcript-density background; none: the background even when images exist")
    cx.add_argument("--z", type=int, help="raw morphology z-stacks (*_F<fov>_Z<z>.TIF): the plane to use (default the middle one)")
    cx.add_argument("--channels", help="comma-separated names for the raw morphology channels (default: the "
                                       "Mean.* columns of the metadata file, in their order)")
    cx.add_argument("--fov-y", choices=["auto", "up", "down"], default="auto",
                    help="whether y_local runs up from the bottom edge of a FOV image or down from the top; "
                         "auto checks against CellLabels when present, else up")
    mc = sub.add_parser("merscope", help="MERSCOPE region folder: mosaic stains (one z-plane), cells (cell_by_gene + "
                                         "cell_metadata) as a 'cells' expression layer, detected_transcripts as the "
                                         "molecules layer (merscope.py)")
    mc.add_argument("--src", required=True, help="the region folder (holding detected_transcripts.csv and images/)")
    mc.add_argument("--images", choices=["auto", "none"], default="auto",
                    help="none: a transcript-density background even when mosaic images exist")
    mc.add_argument("--z", type=int, help="z-plane of the mosaic images (default the middle one present)")
    mc.add_argument("--cells-dir", help="folder holding cell_by_gene.csv and cell_metadata.csv when it is not "
                                        "--src itself (e.g. one of several segmentation runs)")
    for p_ in (cx, mc):
        p_.add_argument("--no-cells", action="store_true", help="skip the cells layer")
        p_.add_argument("--no-transcripts", action="store_true", help="skip the molecules layer")
        p_.add_argument("--background-um", type=float,
                        help="no images: µm per pixel of the density background (default 1)")
        common(p_)
    v = sub.add_parser("visium", help="Space Ranger outs + the full-res image it was run on")
    v.add_argument("--image", required=True, help="the --image TIFF given to spaceranger (RGB or grey)")
    v.add_argument("--outs", required=True)
    v.add_argument("--bins", help="Visium HD layers to build: a comma list of layer keys, e.g. 002um,008um,016um, "
                                  "or 'all'. Default 008um,016um (the 2 µm layer is by far the largest); "
                                  "classic Visium spots are built by default too")
    v.add_argument("--crop", help="x0,y0,x1,y1 in full-res pixels; bin positions move with it")
    v.add_argument("--flip", help="h, v or hv - mirror the image and the bins together")
    v.add_argument("--image-scale", type=float, default=1.0,
                   help="multiply spaceranger positions by this if --image differs in scale")
    common(v)
    ss = sub.add_parser("stereoseq", help="Stereo-seq (STOmics SAW) output: registered image + bin GEF or GEM "
                                          "(square-bin layers) + cell bin GEF (a 'cells' layer)")
    stereoseq.add_args(ss)
    common(ss)
    for p in sub.choices.values():
        p.add_argument("--threads", type=int,
                       help="BLAS / OpenMP threads for numpy and scipy (default: min(8, CPUs this process may use); "
                            "also the SODA_THREADS environment variable). See threadcap.py.")
    args = ap.parse_args(argv)
    if args.threads is not None and args.threads < 1:
        ap.error("--threads must be 1 or more")

    if args.kind == "rnaseq":
        import rnaseq
        if args.h5ad and Path(args.h5ad).stat().st_size > 2e9:          # large files: stream, never load
            rnaseq.export_h5ad_streaming(Path(args.h5ad), Path(args.dest), args.id, args.label or args.id,
                                         embedding=args.embedding, tmpdir=args.tmp, log=log)
            return
        if args.h5ad:
            X, genes, emb, obs = rnaseq.read_h5ad(Path(args.h5ad), args.embedding)
        elif args.mtx and args.umap:
            X, genes, emb, obs = rnaseq.read_10x(Path(args.mtx), Path(args.umap), [Path(p) for p in args.obs])
        else:
            sys.exit("rnaseq needs --h5ad, or --mtx together with --umap")
        rnaseq.export_rnaseq(X, genes, emb, obs, Path(args.dest), args.id, args.label or args.id, log=log)
        return

    if args.kind == "export":
        from export import export_static
        log(f"[export] {threadcap.describe()}")
        export_static(Path(args.out) / args.id, Path(args.dest), fmt=args.fmt, workers=args.workers, log=log,
                      hvg=not args.no_hvg, umap=not args.no_umap,
                      celltypes=not (args.no_celltypes or args.no_umap), markers=Path(args.markers) if args.markers else None,
                      grid_um=args.grid_um, tissue_smooth_um=args.tissue_smooth_um, umap_sample=args.umap_sample,
                      seam_correct=args.seam_correct, user_clusters=userclusters.specs_from_args(args))
        return

    if args.kind == "userclusters":
        specs = userclusters.specs_from_args(args)
        if not specs:
            sys.exit("userclusters needs --user-layer and --user-umap / --user-obs / --user-h5ad")
        userclusters.apply_to_export(Path(args.out) / args.id, Path(args.dest), specs,
                                     markers=Path(args.markers) if args.markers else None,
                                     celltypes=not args.no_celltypes, log=log)
        return

    if args.kind == "inspect":
        import tifffile
        with tifffile.TiffFile(args.src) as tf:
            s0 = tf.series[0]
            print(f"file: {args.src}")
            print(f"  bigtiff={tf.is_bigtiff} qpi={getattr(tf, 'is_qpi', None)} ome={tf.is_ome} series={len(tf.series)}")
            print(f"  series[0]: axes={s0.axes} shape={s0.shape} dtype={s0.dtype} levels={len(s0.levels)}")
            for L, lv in enumerate(s0.levels[:10]):
                print(f"    level {L}: {lv.shape}")
            pg = s0.pages[0]
            print(f"  page0: tiled={pg.is_tiled} tile={pg.tilewidth}x{pg.tilelength} compression={pg.compression.name}")
        src = pyramid.TiffStackSource(Path(args.src), rgb_split=True)
        print(f"  channels ({src.nchannels}): " + ", ".join(src.names))
        src.close()
        return

    if args.kind == "geomx":
        _geomx.prepare(args, log=log, geojson=parse_geojson_args)
        return

    ds_dir = Path(args.out) / args.id
    ds_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    manifest = {"id": args.id, "label": args.label or args.id, "kind": args.kind,
                "pixelSizeUm": args.pixel_size, "created": datetime.now().isoformat(timespec="seconds"),
                "source": {}, "expr": [], "overlays": []}

    if args.kind == "akoya":
        files = akoya_files(Path(args.src))
        if not files:
            sys.exit(f"no TIFFs in {args.src}")
        log(f"[akoya] {len(files)} markers: " + ", ".join(n for n, _ in files))
        src = pyramid.TiffFileSource(files, rgb_split=args.rgb_split)
        manifest["source"] = {"folder": str(args.src), "files": [str(p) for _, p in files]}
    elif args.kind == "xenium":
        import re as _re
        import tifffile as _tf
        folder = Path(args.src)
        if not (folder / "morphology_focus").is_dir() and (folder / "data" / "morphology_focus").is_dir():
            folder = folder / "data"                      # an archive's <name>/data/ layout
        files = sorted((folder / "morphology_focus").glob("*.ome.tif*"))
        if not files:
            sys.exit(f"no morphology_focus/*.ome.tif in {folder}")
        with _tf.TiffFile(str(files[0])) as tf0:           # the OME block lists every channel
            ome = tf0.ome_metadata or ""
        names = _re.findall(r'<Channel[^>]*?\sName="([^"]*)"', ome)
        px = _re.findall(r'PhysicalSizeX="([^"]*)"', ome)
        pairs = []
        for f in files:
            m = _re.match(r"ch(\d+)_", f.name)
            idx = int(m.group(1)) if m else len(pairs)
            pairs.append((names[idx] if idx < len(names) else f.stem.split("_", 1)[-1], f))
        src = pyramid.TiffFileSource(pairs)
        if args.pixel_size is None and px:
            manifest["pixelSizeUm"] = float(px[0])
        import xenium as _xen
        xpx = _xen.pixel_size(folder)
        if args.pixel_size is None and xpx:
            old = manifest["pixelSizeUm"]
            if old and abs(old - xpx) / xpx > 0.01:
                log(f"[xenium] pixel size: experiment.xenium says {xpx:.4f} µm/px, the OME metadata {old:.4f}; "
                    f"using experiment.xenium's (cells and transcripts are placed with it)")
            manifest["pixelSizeUm"] = xpx
        ann = folder / "Custom_annotation.geojson"
        if ann.exists() and not args.geojson:
            args.geojson = [f"Annotations={ann}"]
        manifest["source"] = {"folder": str(folder), "files": [str(p) for _, p in pairs]}
    elif args.kind in ("cosmx", "merscope"):
        import importlib
        plat = importlib.import_module(args.kind).plan(args, log=log)      # cosmx.py / merscope.py
        src = plat.source
        manifest["pixelSizeUm"] = plat.pixel_size_um
        manifest["source"] = plat.manifest_source
    elif args.kind == "tiff":
        names = args.channels.split(",") if args.channels else None
        src = pyramid.TiffStackSource(Path(args.src), names=names, rgb_split=True)
        manifest["source"] = {"file": str(args.src)}
    elif args.kind == "stereoseq":
        stereoseq.precheck(args, log=log)
        src = pyramid.TiffStackSource(Path(args.image), rgb_split=True)
        if src.nchannels == 1:
            src.names = [next((n for n in ("ssDNA", "DAPI") if n.lower() in Path(args.image).name.lower()), "Image")]
        manifest["source"] = {"image": str(args.image), "gef": args.gef, "gem": args.gem, "cellbin": args.cellbin}
    else:
        src = pyramid.TiffStackSource(Path(args.image), rgb_split=True)
        manifest["source"] = {"image": str(args.image), "outs": str(args.outs)}

    if getattr(args, "crop", None):
        try:
            box = [int(v) for v in str(args.crop).replace(" ", "").split(",")]
            if len(box) != 4:
                raise ValueError
        except ValueError:
            raise SystemExit(f"--crop wants four numbers x0,y0,x1,y1 (got {args.crop!r})")
        full_w, full_h = src.width, src.height
        src = pyramid.CropSource(src, box)
        manifest["crop"] = {"box": list(src.box), "ofFullSize": [full_w, full_h]}
        log(f"[image] crop {src.box} of {full_w}x{full_h} -> {src.width}x{src.height}")

    flip = (getattr(args, "flip", "") or "").lower()
    if flip:
        if set(flip) - set("hv"):
            raise SystemExit(f"--flip takes h, v or hv (got {args.flip!r})")
        src = pyramid.FlipSource(src, flip)
        manifest["flip"] = flip
        log(f"[image] flip {flip} ({'horizontal ' if 'h' in flip else ''}{'vertical' if 'v' in flip else ''})")

    if manifest["pixelSizeUm"] is None and getattr(src, "pixel_size_um", None):
        manifest["pixelSizeUm"] = float(src.pixel_size_um)
        log(f"[image] pixel size from file: {manifest['pixelSizeUm']:.4f} µm/px")
    try:
        manifest["image"] = pyramid.build_image(src, ds_dir, tile=args.tile, log=log)
    finally:
        src.close()

    if args.kind == "xenium":
        import xenium as _xen
        inp = _xen.find_inputs(folder)
        um = manifest["pixelSizeUm"]
        size = (manifest["image"]["width"], manifest["image"]["height"])
        if not um:
            log("[xenium] WARNING no pixel size (experiment.xenium / OME): cells and transcripts skipped; "
                "give --pixel-size")
        else:
            if not args.no_cells:
                cells = _xen.read_cells(inp, log=log)
                if cells:
                    area = cells["medianAreaUm2"]
                    diam = 2.0 * (area / 3.141592653589793) ** 0.5 if area else None
                    meta = build_cell_layer(ds_dir, cells["xy_um"], um, cells["X"], cells["genes"], key="cells",
                                            label="Cells", cell_diameter_um=diam, log=log,
                                            extra={"platform": "xenium"}, barcodes=cells["barcodes"])
                    manifest["expr"].append(meta)
            if not args.no_transcripts:
                mol = _xen.write_transcripts(inp, ds_dir, size, um, qv_min=args.qv_min, log=log)
                if mol:
                    manifest["molecules"] = mol

    if args.kind in ("cosmx", "merscope"):
        plat.finish(ds_dir, manifest, args, log=log)                       # cells layer, molecules
    if args.kind == "stereoseq":
        stereoseq.build(args, ds_dir, manifest, log=log)

    if args.kind == "visium":
        layers = find_visium_layers(Path(args.outs))
        for layer in select_visium_layers(layers, args.bins):
            meta = build_expr(layer, ds_dir, image_scale=args.image_scale, log=log,
                              crop_box=(manifest.get("crop") or {}).get("box"),
                              flip=manifest.get("flip", ""),
                              image_size=(manifest["image"]["width"], manifest["image"]["height"]))
            manifest["expr"].append(meta)
            # Space Ranger's microns_per_pixel is measured on this exact image grid, so it wins
            # over the TIFF's own resolution tag; only --pixel-size overrides it.
            sr_um = meta.get("micronsPerPixel")
            if sr_um and args.pixel_size is None:
                old = manifest["pixelSizeUm"]
                if old and abs(old - sr_um) / sr_um > 0.05:
                    log(f"[visium] pixel size: Space Ranger says {sr_um:.4f} µm/px, the image file "
                        f"{old:.4f}; using Space Ranger's")
                manifest["pixelSizeUm"] = sr_um
        if not layers:
            log("[visium] WARNING: no count layers found under", args.outs)

    manifest["overlays"] = parse_geojson_args(args.geojson, ds_dir)
    with open(ds_dir / "manifest.json", "w") as f:
        json.dump(manifest, f, indent=1)
    log(f"[done] {ds_dir}\\manifest.json  ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
