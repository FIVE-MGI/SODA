r"""Gene packs: the per-gene files of an export as a few large shards, read with HTTP Range requests.

An export used to write one small file per gene: expr/<key>/gene/<NAME>.bin for every expressed
gene of every layer (16,000-18,000 per Visium HD layer), molecules/gene/<NAME>.bin for imaging
platforms, gene/<NAME>.bin for an RNA-seq set. A hundred datasets were about 6 million files,
and writing them took about half of an export. A pack holds the same bytes end to end in shards
of about 64 MB; the viewer reads one gene with one Range request ("bytes=o-(o+b-1)").

Layout, version 1 (one "group" = the gene files of one expression layer, of the molecules, or of
an RNA-seq set):

    <group>/genepack/<hash>-<i>.pack    the payloads of the group's genes, end to end, byte for
                                        byte what <group>/gene/<NAME>.bin holds in the loose layout,
                                        in listing order; a gene is never split across shards
    listing entry   {"n": ..., "s": shard, "o": byte offset in the shard, "b": bytes}
                    (genes.json, or the "genes" of molecules/index.json); a gene with no counts
                    has no s / o / b, as it had no file
    manifest        "genePack": {"version": 1, "dir": "genepack", "unit": 12, "genes": n,
                                 "bytes": total, "shards": [{"file": "<hash>-0.pack", "bytes": N}, ...]}
                    on the expr layer, on the "molecules" block, or at the top of an RNA-seq manifest

<hash> is a hash of every name and payload of the group, so a listing can never point into
shards of another export run; the shard sizes in the manifest are checked by the viewer too.

    python genepack.py check  <export>/<id> ...                 index vs shards, every byte accounted for
    python genepack.py repack <export>/<id> ... [--delete-loose] [--shard-mb 64]
                       turn an existing export's gene files into packs (no cache needed), verify
                       every gene against its file, then (only with --delete-loose) delete the files
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
from pathlib import Path
from urllib.parse import quote

VERSION = 1
DIR = "genepack"
SHARD_MB = 64
MODES = ("pack", "loose", "both")
DEFAULT_MODE = "pack"
# fields a listing entry may carry, per kind of group (anything else in a published listing is
# refused by `check`: a listing holds gene names and numbers, nothing about a sample)
FIELDS = {"expr": {"n", "t", "s", "o", "b"}, "molecules": {"n", "c", "s", "o", "b"},
          "rnaseq": {"n", "p", "m", "s", "o", "b"}}


def encode(name: str) -> str:
    """A gene's loose file name: identical to JavaScript encodeURIComponent."""
    return quote(name, safe="-_.!~*'()")


class PackWriter:
    """Write a group's gene payloads into shards, in the order they are added.

    add() returns the index fields {s, o, b} of the gene; finish() renames the shards to their
    final names and returns the manifest's genePack block (None when nothing was added). Shards
    are written as <pid>-<i>.pack.part and renamed only at the end, so a crash never leaves a
    short shard under a name a listing points at."""

    def __init__(self, group_dir: Path, shard_mb: float = SHARD_MB, unit: int | None = None, log=print):
        self.dir = Path(group_dir) / DIR
        self.limit = max(1, int(shard_mb * (1 << 20)))
        self.unit = unit
        self.log = log
        self.h = hashlib.blake2b(digest_size=8)
        self.parts: list[Path] = []
        self.sizes: list[int] = []
        self.f = None
        self.n = 0

    def _open(self):
        if self.f is not None:
            self.f.close()
        self.dir.mkdir(parents=True, exist_ok=True)
        p = self.dir / f"{os.getpid()}-{len(self.parts)}.pack.part"
        self.f = open(p, "wb", buffering=8 << 20)
        self.parts.append(p)
        self.sizes.append(0)

    def add(self, name: str, data) -> dict:
        mv = memoryview(data).cast("B")
        b = mv.nbytes
        if b == 0:
            raise ValueError(f"{name}: empty payload (a gene with no counts gets no pack entry)")
        if self.unit and b % self.unit:
            raise ValueError(f"{name}: {b} bytes is not a multiple of {self.unit}")
        if self.f is None or (self.sizes[-1] > 0 and self.sizes[-1] + b > self.limit):
            self._open()
        o = self.sizes[-1]
        self.f.write(mv)
        self.sizes[-1] += b
        nb = name.encode("utf-8")
        self.h.update(len(nb).to_bytes(4, "little") + nb + b.to_bytes(8, "little"))
        self.h.update(mv)
        self.n += 1
        return {"s": len(self.parts) - 1, "o": o, "b": b}

    def abort(self):
        if self.f is not None:
            self.f.close()
            self.f = None
        for p in self.parts:
            try:
                p.unlink()
            except OSError:
                pass

    def finish(self) -> dict | None:
        if self.f is not None:
            self.f.close()
            self.f = None
        if not self.parts:
            return None
        tag = self.h.hexdigest()
        shards = []
        for i, (p, size) in enumerate(zip(self.parts, self.sizes)):
            final = self.dir / f"{tag}-{i}.pack"
            os.replace(p, final)
            shards.append({"file": final.name, "bytes": int(size)})
        keep = {s["file"] for s in shards}
        stale = sorted(q.name for q in self.dir.iterdir()
                       if q.is_file() and (q.suffix == ".pack" or q.name.endswith(".pack.part")) and q.name not in keep)
        if stale:
            self.log(f"[genepack] {self.dir}: {len(stale)} shard(s) of an earlier export are not used any more "
                     f"({', '.join(stale[:3])}{', ...' if len(stale) > 3 else ''}); delete them before uploading")
        pack = {"version": VERSION, "dir": DIR, "genes": self.n, "bytes": int(sum(self.sizes)), "shards": shards}
        if self.unit:
            pack["unit"] = self.unit
        return pack


def with_index(listing: list[dict], index: dict[str, dict]) -> list[dict]:
    """Listing entries with the {s, o, b} of the genes that were packed (on the first entry of a
    name listed twice: one payload, one entry pointing at it)."""
    out, seen = [], set()
    for e in listing:
        if e["n"] in index and e["n"] not in seen:
            seen.add(e["n"])
            out.append({**e, **index[e["n"]]})
        else:
            out.append(e)
    return out


def atomic_json(p: Path, obj, **kw):
    p = Path(p)
    tmp = p.with_name(p.name + ".part")
    tmp.write_text(json.dumps(obj, separators=(",", ":"), **kw), encoding="utf-8")
    os.replace(tmp, p)


def note_loose(group_dir: Path, log=print):
    """A gene/ folder left by an earlier export in the loose layout: not read any more."""
    g = Path(group_dir) / "gene"
    if g.is_dir():
        n = sum(1 for _ in g.iterdir())
        if n:
            log(f"[genepack] {g}: {n:,} loose gene files of an earlier export are not used any more "
                f"(the layer is packed now); delete that folder before uploading to save the files")


# ---- reading: groups of an export, one gene, checks --------------------------------------------
def groups(m: dict) -> list[dict]:
    """The gene-file groups of an exported dataset's manifest: where the listing is, the listing's
    gene entries (a list, or the "genes" of a dict), the pack block (None = loose files)."""
    out = []
    if m.get("kind") == "rnaseq":
        out.append({"kind": "rnaseq", "dir": "", "listing": "genes.json", "key": None,
                    "pack": m.get("genePack"), "what": "RNA-seq genes"})
    for L in m.get("expr", []) or []:
        out.append({"kind": "expr", "dir": f"expr/{L['key']}", "listing": "genes.json", "key": None,
                    "pack": L.get("genePack"), "what": f"expr {L['key']}"})
    mo = m.get("molecules")
    if mo:
        out.append({"kind": "molecules", "dir": mo.get("dir", "molecules"), "listing": "index.json", "key": "genes",
                    "pack": mo.get("genePack"), "what": "molecules"})
    return out


def _gdir(ds_dir: Path, g: dict) -> Path:
    return Path(ds_dir) / g["dir"] if g["dir"] else Path(ds_dir)


def read_listing(ds_dir: Path, g: dict) -> list[dict]:
    p = _gdir(ds_dir, g) / g["listing"]
    obj = json.loads(p.read_text(encoding="utf-8").replace(":NaN", ":null").replace(":Infinity", ":null"))
    return obj[g["key"]] if g["key"] else obj


def gene_bytes(ds_dir: Path, g: dict, entry: dict) -> bytes | None:
    """One gene's payload, from the pack when the group is packed, else from its loose file.
    None = the gene has none (no counts)."""
    d = _gdir(ds_dir, g)
    pk = g.get("pack")
    if pk:
        if "s" not in entry:
            return None
        with open(d / pk["dir"] / pk["shards"][entry["s"]]["file"], "rb") as f:
            f.seek(entry["o"])
            return f.read(entry["b"])
    p = d / "gene" / f"{encode(entry['n'])}.bin"
    return p.read_bytes() if p.exists() else None


def check_group(ds_dir: Path, g: dict, verify_hash: bool = False) -> list[str]:
    """A packed group, structurally: the manifest's shards exist with their sizes, every listing
    entry points inside one shard, the entries tile the shards with no gap and no overlap (so a
    pack holds nothing but the listed genes' payloads), and the listing carries only the fields of
    its kind with names as plain strings. `verify_hash`: also re-read every shard and recompute
    the content hash in the file names. Returns the problems (empty = fine)."""
    pk = g.get("pack")
    if not pk:
        return []
    what = g["what"]
    d = _gdir(ds_dir, g)
    probs = []
    if pk.get("version") != VERSION:
        return [f"{what}: genePack version {pk.get('version')!r} is not {VERSION}"]
    shards = pk.get("shards") or []
    tags = set()
    for i, s in enumerate(shards):
        f = d / pk["dir"] / s["file"]
        tag, _, rest = s["file"].rpartition("-")
        tags.add(tag)
        if rest != f"{i}.pack":
            probs.append(f"{what}: shard {i} is named {s['file']}")
        if not f.is_file():
            probs.append(f"{what}: shard {s['file']} is missing")
        elif f.stat().st_size != s["bytes"]:
            probs.append(f"{what}: shard {s['file']} is {f.stat().st_size} bytes, the manifest says {s['bytes']}")
    if len(tags) > 1:
        probs.append(f"{what}: shards of different export runs: {sorted(tags)}")
    try:
        entries = read_listing(ds_dir, g)
    except (OSError, ValueError, KeyError, TypeError) as e:
        return probs + [f"{what}: listing {g['listing']} unreadable: {e}"]
    allowed = FIELDS[g["kind"]]
    spans = [[] for _ in shards]
    n = 0
    for e in entries:
        extra = set(e) - allowed
        if extra:
            probs.append(f"{what}: listing entry {e.get('n')!r} has fields {sorted(extra)} (allowed {sorted(allowed)})")
        if not isinstance(e.get("n"), str):
            probs.append(f"{what}: listing entry without a gene name: {e!r}"[:200])
        if "s" not in e:
            continue
        s, o, b = e.get("s"), e.get("o"), e.get("b")
        if not all(isinstance(v, int) for v in (s, o, b)) or not 0 <= s < len(shards) or o < 0 or b <= 0:
            probs.append(f"{what}: {e.get('n')!r}: bad pack fields s={s!r} o={o!r} b={b!r}")
            continue
        if pk.get("unit") and b % pk["unit"]:
            probs.append(f"{what}: {e['n']!r}: {b} bytes is not a multiple of {pk['unit']}")
        spans[s].append((o, o + b, e["n"]))
        n += 1
    if n != pk.get("genes"):
        probs.append(f"{what}: {n} packed genes in the listing, the manifest says {pk.get('genes')}")
    total = 0
    for i, sp in enumerate(spans):
        sp.sort()
        at = 0
        for a, z, name in sp:
            if a != at:
                probs.append(f"{what}: shard {i}: {'gap' if a > at else 'overlap'} at byte {min(a, at)} (before {name!r})")
                break
            at = z
        if i < len(shards) and at != shards[i]["bytes"]:
            probs.append(f"{what}: shard {i}: the listing covers {at} of its {shards[i]['bytes']} bytes")
        total += at
    if total != pk.get("bytes"):
        probs.append(f"{what}: the listing covers {total} bytes, the manifest says {pk.get('bytes')}")
    if verify_hash and not probs:
        h = hashlib.blake2b(digest_size=8)
        order = sorted((e for e in entries if "s" in e), key=lambda e: (e["s"], e["o"]))
        handles = {}
        try:
            for e in order:
                f = handles.get(e["s"])
                if f is None:
                    f = handles[e["s"]] = open(d / pk["dir"] / shards[e["s"]]["file"], "rb")
                f.seek(e["o"])
                buf = f.read(e["b"])
                nb = e["n"].encode("utf-8")
                h.update(len(nb).to_bytes(4, "little") + nb + len(buf).to_bytes(8, "little"))
                h.update(buf)
        finally:
            for f in handles.values():
                f.close()
        if shards and h.hexdigest() not in tags:
            probs.append(f"{what}: content hash {h.hexdigest()} does not match the shard names {sorted(tags)}")
    return probs


def check_dataset(ds_dir: Path, verify_hash: bool = False) -> list[str]:
    m = json.loads((Path(ds_dir) / "manifest.json").read_text(encoding="utf-8"))
    probs = []
    for g in groups(m):
        probs += [f"{Path(ds_dir).name}: {p}" for p in check_group(ds_dir, g, verify_hash=verify_hash)]
    return probs


# ---- repack: an existing loose export -> packs ------------------------------------------------
UNITS = {"expr": 12, "molecules": 4, "rnaseq": None}


def repack(ds_dir: Path, shard_mb: float = SHARD_MB, delete_loose: bool = False, log=print) -> dict:
    """Pack the loose gene files of an existing export in place: packs and the new listings and
    manifest are written, every gene is read back from its pack and compared with its file, and
    only then (with delete_loose) are the files deleted. Already packed groups are left alone.
    Returns {group: (genes, shards)}."""
    ds_dir = Path(ds_dir)
    mp = ds_dir / "manifest.json"
    m = json.loads(mp.read_text(encoding="utf-8"))
    done = {}
    new_packs = {}
    for g in groups(m):
        if g["pack"]:
            log(f"[genepack] {ds_dir.name} {g['what']}: already packed")
            continue
        d = _gdir(ds_dir, g)
        if not (d / "gene").is_dir():
            continue
        lp = d / g["listing"]
        obj = json.loads(lp.read_text(encoding="utf-8").replace(":NaN", ":null").replace(":Infinity", ":null"))
        entries = obj[g["key"]] if g["key"] else obj
        w = PackWriter(d, shard_mb=shard_mb, unit=UNITS[g["kind"]], log=log)
        index = {}
        try:
            for e in entries:
                f = d / "gene" / f"{encode(e['n'])}.bin"
                if f.exists() and f.stat().st_size and e["n"] not in index:
                    index[e["n"]] = w.add(e["n"], f.read_bytes())
        except BaseException:
            w.abort()
            raise
        pack = w.finish()
        if not pack:
            continue
        g2 = {**g, "pack": pack}
        new_entries = with_index(entries, index)
        bad = [e["n"] for e in new_entries if "s" in e and
               gene_bytes(ds_dir, g2, e) != (d / "gene" / f"{encode(e['n'])}.bin").read_bytes()]
        if bad:
            raise SystemExit(f"[genepack] {ds_dir.name} {g['what']}: {len(bad)} gene(s) differ after packing "
                             f"({bad[:5]}); nothing changed in the listing or manifest")
        if g["key"]:
            obj[g["key"]] = new_entries
            atomic_json(lp, obj, allow_nan=False)
        else:
            atomic_json(lp, new_entries, allow_nan=False)
        new_packs[(g["kind"], g["dir"])] = pack
        done[g["what"]] = (pack["genes"], len(pack["shards"]))
        log(f"[genepack] {ds_dir.name} {g['what']}: {pack['genes']:,} genes -> {len(pack['shards'])} shard(s), "
            f"{pack['bytes'] / 1e6:.1f} MB, every gene verified against its file")
    if new_packs:
        if m.get("kind") == "rnaseq" and ("rnaseq", "") in new_packs:
            m["genePack"] = new_packs[("rnaseq", "")]
        for L in m.get("expr", []) or []:
            if ("expr", f"expr/{L['key']}") in new_packs:
                L["genePack"] = new_packs[("expr", f"expr/{L['key']}")]
        mo = m.get("molecules")
        if mo and ("molecules", mo.get("dir", "molecules")) in new_packs:
            mo["genePack"] = new_packs[("molecules", mo.get("dir", "molecules"))]
        tmp = mp.with_name("manifest.json.part")
        tmp.write_text(json.dumps(m, indent=1), encoding="utf-8")
        os.replace(tmp, mp)
        probs = check_dataset(ds_dir, verify_hash=True)
        if probs:
            raise SystemExit("[genepack] check failed after repack:\n  " + "\n  ".join(probs))
        if delete_loose:
            for kind, gd in new_packs:
                gdir = ds_dir / gd if gd else ds_dir
                n = 0
                for f in (gdir / "gene").glob("*.bin"):
                    f.unlink()
                    n += 1
                try:
                    (gdir / "gene").rmdir()
                except OSError:
                    pass
                log(f"[genepack] {ds_dir.name} {gd or '.'}: {n:,} loose gene files deleted")
    return done


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("check", help="check the packs of exported datasets against their listings")
    c.add_argument("datasets", nargs="+", help="<export>/<id> folders")
    c.add_argument("--hash", action="store_true", help="also re-read every shard and check its content hash")
    r = sub.add_parser("repack", help="pack the loose gene files of exported datasets in place")
    r.add_argument("datasets", nargs="+", help="<export>/<id> folders")
    r.add_argument("--shard-mb", type=float, default=SHARD_MB)
    r.add_argument("--delete-loose", action="store_true",
                   help="delete the loose gene files once every gene has been verified in its pack")
    a = ap.parse_args(argv)
    if a.cmd == "check":
        probs = []
        for ds in a.datasets:
            probs += check_dataset(Path(ds), verify_hash=a.hash)
        for p in probs:
            print(p)
        print(f"[genepack] {len(probs)} problem(s)" if probs else "[genepack] packs consistent")
        return 1 if probs else 0
    for ds in a.datasets:
        repack(Path(ds), shard_mb=a.shard_mb, delete_loose=a.delete_loose)
    return 0


if __name__ == "__main__":
    sys.exit(main())
