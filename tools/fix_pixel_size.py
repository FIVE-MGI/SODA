#!/usr/bin/env python
r"""Correct the scale-bar pixel size of Visium datasets prepared before the fix.

Until 2026-11-11, prepare.py took the image's TIFF resolution tag over Space Ranger's
microns_per_pixel. Many microscope exports carry a print setting there (300 dpi), which gave
pixelSizeUm = 84.6667 and a scale bar ~200x too long. Only manifest.json's top-level
`pixelSizeUm` was wrong: bins, genes, UMAPs and the circle tool use each layer's own
`micronsPerPixel`, which was always right. So no re-export is needed, only this field.

    python tools/fix_pixel_size.py <folder> [<folder> ...]           # report only
    python tools/fix_pixel_size.py <folder> [<folder> ...] --write   # change manifest.json

A <folder> is one dataset (it has manifest.json) or a folder of datasets (cache or export
root). Run it on the cache too, or the next export copies the old value back. Re-upload each
changed export manifest.json afterwards.
"""
import json
import sys
from pathlib import Path


def manifests(paths):
    for p in map(Path, paths):
        if (p / "manifest.json").exists():
            yield p / "manifest.json"
        else:
            yield from sorted(p.glob("*/manifest.json"))


def main(argv):
    write = "--write" in argv
    paths = [a for a in argv if a != "--write"]
    if not paths:
        print(__doc__)
        return 2
    changed = 0
    for mf in manifests(paths):
        m = json.loads(mf.read_text(encoding="utf-8"))
        sr = next((e.get("micronsPerPixel") for e in m.get("expr", []) if e.get("micronsPerPixel")), None)
        px = m.get("pixelSizeUm")
        if m.get("kind") != "visium" or not sr:
            continue
        if px is not None and abs(px - sr) / sr <= 0.05:
            print(f"ok      {mf.parent.name:32s} {px:.4f} um/px")
            continue
        print(f"{'FIXED' if write else 'WRONG'}   {mf.parent.name:32s} {px} -> {sr:.4f} um/px (Space Ranger)")
        if write:
            m["pixelSizeUm"] = sr
            text = json.dumps(m, indent=1 if mf.read_text(encoding='utf-8').startswith('{\n ') else None)
            mf.write_text(text, encoding="utf-8")
        changed += 1
    print(f"\n{changed} manifest(s) {'changed' if write else 'to change (run again with --write)'}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
