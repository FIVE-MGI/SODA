#!/usr/bin/env python3
r"""Does any binary inside SODA.app need a newer macOS than its Info.plist promises?

    python3 app/check_min_os.py <path to SODA.app>

Reads LSMinimumSystemVersion from Contents/Info.plist and, with `otool -l`, the minimum macOS of
every Mach-O file in the bundle (LC_BUILD_VERSION minos, or LC_VERSION_MIN_MACOSX). Exit 1 and a
list when any is newer than the plist: such an app launches on the older macOS and then dies
loading that library (FIVE Atlas shipped exactly this once, with numpy). macOS only.
"""
import plistlib
import re
import subprocess
import sys
from pathlib import Path


def key(v: str):
    return tuple(int(x) for x in (re.findall(r"\d+", v) + ["0", "0"])[:3])


def min_os(path: Path):
    try:
        out = subprocess.run(["otool", "-l", str(path)], capture_output=True, text=True, timeout=60).stdout
    except Exception:
        return None
    found = re.findall(r"cmd LC_BUILD_VERSION.*?minos (\S+)", out, re.S) \
        + re.findall(r"cmd LC_VERSION_MIN_MACOSX.*?version (\S+)", out, re.S)
    return max(found, key=key) if found else None


def main() -> int:
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    app = Path(sys.argv[1])
    plist = plistlib.loads((app / "Contents" / "Info.plist").read_bytes())
    promised = plist.get("LSMinimumSystemVersion", "0")
    worst, too_new, n = "0", [], 0
    for f in (app / "Contents").rglob("*"):
        if not f.is_file() or f.is_symlink():
            continue
        with open(f, "rb") as h:
            magic = h.read(4)
        if magic not in (b"\xcf\xfa\xed\xfe", b"\xce\xfa\xed\xfe", b"\xca\xfe\xba\xbe"):
            continue
        n += 1
        v = min_os(f)
        if not v:
            continue
        if key(v) > key(worst):
            worst = v
        if key(v) > key(promised):
            too_new.append(f"{v}  {f.relative_to(app)}")
    print(f"{n} Mach-O files; Info.plist promises macOS {promised}; the newest any of them needs: {worst}")
    if too_new:
        print(f"{len(too_new)} need a newer macOS than {promised} (set SODA_MIN_MACOS={worst} or pin older wheels):")
        for line in sorted(too_new, key=lambda s: key(s.split()[0]), reverse=True)[:40]:
            print("  " + line)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
