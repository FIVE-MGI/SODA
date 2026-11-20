r"""Put the page's two outside libraries (deck.gl, umap-js) in vendor/ next to index.html, so SODA
needs no CDN: on a network without internet, in the desktop app, or for a site that would rather
not depend on cdn.jsdelivr.net.

    python web/SpatialViewer/vendor.py                 (download the pinned files, check them, write vendor/)
    python web/SpatialViewer/prepare.py vendor         (the same)
    python web/SpatialViewer/vendor.py --check         (only check what is in vendor/; downloads nothing)
    python web/SpatialViewer/vendor.py --from D:/usb/vendor
                                                       (copy from a folder laid out the same way, e.g. another
                                                        computer's vendor/, instead of downloading)
    python web/SpatialViewer/vendor.py --dest D:/site  (the folder that holds index.html; default: this one)
    python web/SpatialViewer/vendor.py --source https://unpkg.com/
                                                       (another npm mirror; the hash check is the same)

The versions and SHA-384 hashes come from libs.js (the JSON between its LIBS markers), the one list
the page, lasso_worker.js and sw.js read too. A file is written only when its size and SHA-384 are
exactly the pinned ones, so a wrong or damaged download never lands in vendor/; the page checks the
same hash again before it runs the file (Subresource Integrity) and otherwise uses the CDN.

Layout (the CDN's own paths under https://cdn.jsdelivr.net/npm/):

    vendor/deck.gl@9.1.12/dist.min.js
    vendor/umap-js@1.4.0/lib/umap-js.min.js

vendor/ is git-ignored: the files are not part of the repository. Copy the folder along with the
page when you publish it (or run this command on the server) to serve them from your own site.
Exit code 0 when every file is in place and correct, 1 otherwise. Standard library only.
"""
from __future__ import annotations

import argparse
import base64
import hashlib
import json
import os
import re
import sys
import tempfile
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
CDN = "https://cdn.jsdelivr.net/npm/"


def read_libs(libs_js: Path = HERE / "libs.js") -> dict:
    """The pinned libraries from libs.js: {key: {name, version, path, bytes, sri}}."""
    text = libs_js.read_text(encoding="utf-8")
    m = re.search(r"/\* LIBS-JSON \*/(.*?)/\* END-LIBS-JSON \*/", text, re.S)
    if not m:
        raise SystemExit(f"{libs_js}: the LIBS-JSON block is missing")
    libs = json.loads(m.group(1))
    for k, v in libs.items():
        if not (isinstance(v.get("path"), str) and re.fullmatch(r"[\w.@-]+(/[\w.@-]+)+", v["path"]) and ".." not in v["path"]):
            raise SystemExit(f"libs.js: {k} has an unexpected path {v.get('path')!r}")
        if not re.fullmatch(r"sha384-[A-Za-z0-9+/]{64}", v.get("sri", "")):
            raise SystemExit(f"libs.js: {k} has no sha384 hash")
    return libs


def sri(data: bytes) -> str:
    return "sha384-" + base64.b64encode(hashlib.sha384(data).digest()).decode("ascii")


def verdict(lib: dict, data: bytes | None) -> str:
    """'' when the bytes are the pinned file, else what is wrong."""
    if data is None:
        return "missing"
    if lib.get("bytes") and len(data) != lib["bytes"]:
        return f"{len(data):,} bytes, expected {lib['bytes']:,}"
    got = sri(data)
    return "" if got == lib["sri"] else f"SHA-384 {got}, expected {lib['sri']}"


def download(url: str, timeout: float = 60) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "SODA vendor.py", "Accept-Encoding": "identity"})
    with urllib.request.urlopen(req, timeout=timeout) as r:      # noqa: S310 (a fixed https address)
        if r.status != 200:
            raise OSError(f"HTTP {r.status}")
        return r.read()


def write_atomic(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=".part-")
    try:
        with os.fdopen(fd, "wb") as f:
            f.write(data)
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dest", default=str(HERE), help="the folder that holds index.html (vendor/ goes inside it)")
    ap.add_argument("--check", action="store_true", help="only check the files already in vendor/; download nothing")
    ap.add_argument("--from", dest="src", metavar="FOLDER", help="copy from this folder (laid out as vendor/) instead of downloading")
    ap.add_argument("--source", default=CDN, help=f"npm mirror to download from (default {CDN})")
    args = ap.parse_args(argv)
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(errors="replace")
    libs = read_libs()
    vendor = Path(args.dest) / "vendor"
    bad = 0
    for key, lib in libs.items():
        target = vendor / lib["path"]
        have = target.read_bytes() if target.is_file() else None
        why = verdict(lib, have)
        if not why:
            print(f"ok       {lib['name']} {lib['version']}: {target}")
            continue
        if args.check:
            print(f"{'MISSING' if have is None else 'WRONG':8} {lib['name']} {lib['version']}: {target} ({why})")
            bad += 1
            continue
        if args.src:
            src = Path(args.src) / lib["path"]
            where = str(src)
            try:
                data = src.read_bytes()
            except OSError as e:
                print(f"FAILED   {lib['name']} {lib['version']}: cannot read {where} ({e})")
                bad += 1
                continue
        else:
            where = args.source.rstrip("/") + "/" + lib["path"]
            if not where.startswith("https://"):
                print(f"FAILED   {lib['name']}: --source must start with https://")
                bad += 1
                continue
            try:
                data = download(where)
            except Exception as e:                    # noqa: BLE001 (network errors of every kind)
                print(f"FAILED   {lib['name']} {lib['version']}: {where} ({e})")
                bad += 1
                continue
        why = verdict(lib, data)
        if why:
            print(f"REFUSED  {lib['name']} {lib['version']}: {where} is not the pinned file ({why}); nothing written")
            bad += 1
            continue
        write_atomic(target, data)
        print(f"written  {lib['name']} {lib['version']}: {target} ({len(data):,} bytes, {lib['sri']})")
    if bad:
        print(f"\n{bad} of {len(libs)} libraries are not in place: the page uses the CDN for "
              f"{'it' if bad == 1 else 'them'}.")
        return 1
    print(f"\nvendor/ is complete ({len(libs)} libraries): the page loads them from {vendor} instead of the CDN.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
