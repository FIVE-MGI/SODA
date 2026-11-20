"""Single molecules (decoded transcripts) -> one compact point file per gene.

Imaging platforms (Xenium, CosMx, MERSCOPE) locate every transcript. The viewer draws a
picked gene's molecules as small dots, several genes at once, loading each gene's file only
when it is picked. This module is platform-neutral: an adapter (xenium.py) streams batches of
(gene id, x µm, y µm) that already passed its own quality filter into MoleculeWriter, which
writes

    molecules/index.json          {"format": "xy-uint16", "n", "nGenes", "bounds", "scale", ...,
                                   "genes": [{"n": name, "c": count}, ...]}
    molecules/gene/<NAME>.bin     uint16 pairs qx, qy per molecule, little-endian, no header

Coordinates are full-resolution image pixels (µm / pixel size, the same frame as the tiles and
expr/ positions), quantised to 16 bits over the image: x = qx * scale[0], y = qy * scale[1],
with scale = image size / 65535. For a 50,000 px wide Xenium image that is 0.76 px (0.16 µm)
per step, under the platform's localisation precision, at 4 bytes a molecule (float32 x, y
would be 8).

Within a file the molecules are in a seeded random order, so any prefix is a uniform random
sample of the gene. The viewer draws a prefix when the whole section is in view (a few million
dots in all) and the whole file once zoomed in, without a second file or a second request.
Molecules outside the image are dropped and counted. Same input, same bytes.

That is the cache's layout (server.py serves it as it is). The static export packs the gene
files like expr/ (genepack.py): molecules/genepack/<hash>-<i>.pack, each gene's place as
"s", "o", "b" on its index.json entry and the shard list as "genePack" on the manifest's
"molecules" block.
"""
from __future__ import annotations

import json
import time
from pathlib import Path
from urllib.parse import quote

import numpy as np

import deid

DIR = "molecules"
FORMAT = "xy-uint16"
Q = 65535
SEED = 42
# consolidate each gene's pending batches into one array after this many batches, so a
# 200-million-molecule section does not hold a million small arrays
_CONSOLIDATE_EVERY = 64


def encode(name: str) -> str:
    """File name of a gene, identical to JavaScript encodeURIComponent."""
    return quote(name, safe="-_.!~*'()")


class MoleculeWriter:
    """Collect molecules batch by batch, then write one shuffled uint16 file per gene.

    `gene_names`: every gene id's name (ids index this list). `keep`: bool per gene id, False
    for ids never written (controls); their molecules are counted as `controls`.
    `image_size`: (width, height) in full-resolution pixels; `pixel_size_um`: µm per pixel.
    Memory: 4 bytes per kept molecule until `finish`.
    """

    def __init__(self, gene_names, image_size, pixel_size_um: float, keep=None, log=print):
        self.names = [str(g) for g in gene_names]
        self.G = len(self.names)
        self.keep = np.ones(self.G, bool) if keep is None else np.asarray(keep, bool)
        if len(self.keep) != self.G:
            raise ValueError("keep must have one entry per gene name")
        self.W, self.H = float(image_size[0]), float(image_size[1])
        if not (self.W > 0 and self.H > 0 and pixel_size_um and pixel_size_um > 0):
            raise ValueError(f"molecules need the image size and pixel size (got {image_size}, {pixel_size_um})")
        self.px = float(pixel_size_um)
        self.sx, self.sy = self.W / Q, self.H / Q
        self.parts: list[list[np.ndarray]] = [[] for _ in range(self.G)]
        self.counts = np.zeros(self.G, np.int64)
        self.stats = {"read": 0, "controls": 0, "outside": 0, "badGene": 0}
        self.batches = 0
        self.log = log
        self.t0 = time.time()

    def add(self, gene_ids, x_um, y_um):
        """One batch of molecules that passed the platform's quality filter."""
        g = np.asarray(gene_ids).astype(np.int64, copy=False).ravel()
        x = np.asarray(x_um, dtype=np.float64).ravel() / self.px
        y = np.asarray(y_um, dtype=np.float64).ravel() / self.px
        self.stats["read"] += len(g)
        ok = (g >= 0) & (g < self.G)
        self.stats["badGene"] += int((~ok).sum())
        g, x, y = g[ok], x[ok], y[ok]
        k = self.keep[g]
        self.stats["controls"] += int((~k).sum())
        g, x, y = g[k], x[k], y[k]
        inside = (x >= 0) & (x <= self.W) & (y >= 0) & (y <= self.H)
        self.stats["outside"] += int((~inside).sum())
        g, x, y = g[inside], x[inside], y[inside]
        qx = np.minimum(np.rint(x / self.sx), Q).astype(np.uint32)
        qy = np.minimum(np.rint(y / self.sy), Q).astype(np.uint32)
        packed = qx | (qy << 16)                   # little-endian: qx then qy as two uint16
        order = np.argsort(g, kind="stable")
        g, packed = g[order], packed[order]
        bounds = np.flatnonzero(np.diff(g)) + 1
        starts = np.r_[0, bounds]
        ends = np.r_[bounds, len(g)]
        for s, e in zip(starts, ends):
            if e > s:
                gi = int(g[s])
                self.parts[gi].append(packed[s:e].copy())
                self.counts[gi] += e - s
        self.batches += 1
        if self.batches % _CONSOLIDATE_EVERY == 0:
            for gi, p in enumerate(self.parts):
                if len(p) > 1:
                    self.parts[gi] = [np.concatenate(p)]

    def finish(self, ds_dir: Path, meta: dict | None = None) -> dict:
        """Write molecules/ under `ds_dir`; returns the manifest `molecules` block."""
        d = Path(ds_dir) / DIR
        (d / "gene").mkdir(parents=True, exist_ok=True)
        rng = np.random.default_rng(SEED)
        genes, total, nbytes = [], 0, 0
        for gi, name in enumerate(self.names):
            if not self.keep[gi]:
                continue
            p = self.parts[gi]
            a = np.concatenate(p) if len(p) > 1 else (p[0] if p else np.zeros(0, np.uint32))
            self.parts[gi] = []
            n = len(a)
            genes.append({"n": name, "c": int(n)})
            if not n:
                continue                           # listed with c = 0, no file (as expr/ genes)
            a = a[rng.permutation(n)]
            buf = a.astype("<u4").tobytes()        # = uint16 qx, qy per molecule
            f = d / "gene" / f"{encode(name)}.bin"
            tmp = f.with_suffix(".bin.part")
            tmp.write_bytes(buf)
            tmp.replace(f)
            total += n
            nbytes += len(buf)
        deid.assert_clean([g["n"] for g in genes], "molecules gene names")
        index = {"format": FORMAT, "n": int(total), "nGenes": len(genes),
                 "bounds": [0.0, 0.0, self.W, self.H], "scale": [self.sx, self.sy],
                 "units": "full-resolution image pixels", "pixelSizeUm": self.px,
                 "order": "seeded random: any prefix is a uniform sample of the gene",
                 "excluded": {"controls": self.stats["controls"], "outsideImage": self.stats["outside"],
                              "unknownGeneId": self.stats["badGene"]},
                 **(meta or {}), "genes": genes}
        deid.assert_clean({k: v for k, v in index.items() if k != "genes"}, "molecules/index.json")
        tmp = d / "index.json.part"
        tmp.write_text(json.dumps(index, separators=(",", ":"), allow_nan=False), encoding="utf-8")
        tmp.replace(d / "index.json")
        block = {"dir": DIR, "format": FORMAT, "n": int(total), "nGenes": len(genes), "bytes": int(nbytes),
                 **{k: v for k, v in (meta or {}).items() if k in ("qvMin", "platform")}}
        top = sorted(genes, key=lambda g: -g["c"])[:3]
        self.log(f"[molecules] {total:,} molecules of {len(genes)} genes, {nbytes / 1e6:.1f} MB "
                 f"({self.stats['controls']:,} control, {self.stats['outside']:,} outside the image left out); "
                 f"largest " + ", ".join(f"{g['n']} {g['c']:,}" for g in top)
                 + f"  ({time.time() - self.t0:.0f}s)")
        return block


def write_molecules(ds_dir: Path, gene_names, batches, image_size, pixel_size_um: float,
                    keep=None, meta: dict | None = None, log=print) -> dict:
    """Generic entry point for any platform adapter.

    `batches` yields (gene_ids, x_um, y_um) arrays, already quality-filtered by the adapter
    (or pass a single tuple in a list). Returns the manifest `molecules` block."""
    w = MoleculeWriter(gene_names, image_size, pixel_size_um, keep=keep, log=log)
    for g, x, y in batches:
        w.add(g, x, y)
    return w.finish(ds_dir, meta=meta)


def export_molecules(ds_dir: Path, out: Path, log=print) -> tuple[int, dict | None]:
    """The cache's molecules/ into the export: the gene files packed in index order (genepack.py),
    index.json with each gene's s / o / b. Returns (genes written, the genePack block)."""
    import genepack
    src = Path(ds_dir) / DIR
    if not (src / "index.json").exists():
        return 0, None
    dst = Path(out) / DIR
    dst.mkdir(parents=True, exist_ok=True)
    index = json.loads((src / "index.json").read_text(encoding="utf-8"))
    w = genepack.PackWriter(dst, unit=4, log=log)
    where = {}
    try:
        for g in index["genes"]:
            f = src / "gene" / f"{encode(g['n'])}.bin"
            if g.get("c") and f.exists() and g["n"] not in where:
                where[g["n"]] = w.add(g["n"], f.read_bytes())
    except BaseException:
        w.abort()
        raise
    pack = w.finish()
    index = {**index, "genes": genepack.with_index(index["genes"], where)}
    n = len(where)
    genepack.note_loose(dst, log=log)
    log(f"[export] molecules: {n:,} genes -> "
        + (f"{len(pack['shards'])} pack shard(s), {pack['bytes'] / 1e6:.1f} MB" if pack else "no pack (no molecules)"))
    deid.assert_clean([g["n"] for g in index["genes"]], "molecules gene names")
    genepack.atomic_json(dst / "index.json", index, allow_nan=False)
    return n, pack
