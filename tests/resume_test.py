#!/usr/bin/env python
r"""A prepare or export that fails part-way leaves no manifest.json, so the next batch run does it again.

    python tests/resume_test.py [--work <folder>] [--keep]        (from the repo root, with soda-env's python)
    python -m pytest tests/resume_test.py                         (same test, if pytest is installed)

No data needed: a small two-channel TIFF made here (about a minute in all). Steps:
  1. `prepare.py batch` on it: ok, with <cache>/tiny/manifest.json and <export>/tiny/manifest.json.
  2. The export again, in this process, made to fail after its tiles (export.verify_tiles raises): the old
     <export>/tiny/manifest.json must be gone.
  3. `prepare.py batch` (no --rebuild): the export is run again (the prepare is not: the cache is whole), ok.
  4. <export>/tiny/ deleted (a clean export), then the prepare again, in this process, made to fail after its image
     (pyramid.build_image raises once it has rewritten image.zarr): no <cache>/tiny/manifest.json may be left.
  5. `prepare.py batch` (no --rebuild): prepare and export both run again (before, the half-rebuilt cache was
     taken as whole and exported), ok, both manifests back.
  6. The failure reasons of batch.py (failure_reason): out of memory from the message or the exit code, else the
     step's last error line, never a traceback frame.
Exit code 0 = every check passed.
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
APP = REPO / "web" / "SpatialViewer"
DS = "tiny"


class Simulated(RuntimeError):
    pass


def make_input(root: Path) -> Path:
    import numpy as np
    import tifffile
    src = root / "inputs" / DS
    src.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(0)
    img = (rng.random((2, 700, 900)) * 4000).astype(np.uint16)
    img[0, 200:500, 300:700] += 20000                        # some "tissue" for the protein grid
    tifffile.imwrite(src / f"{DS}.ome.tif", img, metadata={"axes": "CYX", "Channel": {"Name": ["DAPI", "CD3"]},
                                                              "PhysicalSizeX": 0.5, "PhysicalSizeY": 0.5}, ome=True)
    return src / f"{DS}.ome.tif"


def run_batch(root: Path) -> dict:
    cmd = [sys.executable, str(APP / "prepare.py"), "batch", "--src", str(root / "inputs"), "--out",
           str(root / "cache"), "--dest", str(root / "export"), "--logs", str(root / "logs"), "--threads", "2"]
    r = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8", errors="replace")
    if r.returncode != 0:
        print(r.stdout[-4000:], r.stderr[-4000:])
        raise AssertionError(f"batch exit {r.returncode}")
    newest = max((root / "logs").glob("batch_*.json"), key=lambda p: p.stat().st_mtime)
    rows = {d["id"]: d for d in json.loads(newest.read_text(encoding="utf-8"))["datasets"]}
    return rows[DS]


def check(cond, what):
    print(("  ok    " if cond else "  FAIL  ") + what, flush=True)
    if not cond:
        raise AssertionError(what)


def export_failing(root: Path):
    import export
    real = export.verify_tiles

    def boom(*a, **k):
        real(*a, **k)
        raise Simulated("export stopped part-way (test)")
    export.verify_tiles = boom
    try:
        import prepare
        prepare.main(["export", "--id", DS, "--out", str(root / "cache"), "--dest", str(root / "export"),
                      "--no-umap", "--threads", "2"])
    except Simulated:
        return
    finally:
        export.verify_tiles = real
    raise AssertionError("the simulated export failure did not happen")


def prepare_failing(root: Path, tif: Path):
    import pyramid
    real = pyramid.build_image

    def boom(*a, **k):
        real(*a, **k)                                       # image.zarr rewritten (mode "w"), then the crash
        raise Simulated("prepare stopped part-way (test)")
    pyramid.build_image = boom
    try:
        import prepare
        prepare.main(["tiff", "--src", str(tif), "--id", DS, "--out", str(root / "cache"), "--threads", "2"])
    except Simulated:
        return
    finally:
        pyramid.build_image = real
    raise AssertionError("the simulated prepare failure did not happen")


def reasons():
    import batch
    tb = ["Traceback (most recent call last):", '  File "x.py", line 3, in <module>', "    g.create_dataset(1)",
          "    ^^^^^^^^^^^^^^^^", "AttributeError: 'Group' object has no attribute 'create_dataset'"]
    r = batch.failure_reason("prepare", {"exit": 1, "tail": tb})
    check(r == "prepare failed (exit 1): AttributeError: 'Group' object has no attribute 'create_dataset'", r)
    r = batch.failure_reason("prepare", {"exit": 1, "tail": ["[xenium] reading cells", "[xenium] --src x not found"]})
    check(r == "prepare failed (exit 1): [xenium] --src x not found", r)
    mem = ["[image] 2 channels", "numpy._core._exceptions._ArrayMemoryError: Unable to allocate 41.0 GiB for an array"]
    r = batch.failure_reason("export", {"exit": 1, "tail": mem, "peakBytes": 30 * 2**30})
    check(r.startswith("export ran out of memory - try fewer datasets at once") and "peak" in r, r)
    check(batch.failure_reason("prepare", {"exit": 1, "tail": ["MemoryError"]}).startswith("prepare ran out of memory"),
          "a bare MemoryError")
    check(batch.failure_reason("prepare", {"exit": -9, "tail": []}).startswith("prepare ran out of memory"),
          "killed with SIGKILL (exit -9)")
    r = batch.failure_reason("prepare", {"exit": 0xC0000005, "tail": []})
    check(r == "prepare failed (exit 3221225477 = 0xC0000005)", r)
    r = batch.failure_reason("prepare", {"exit": 1, "tail": ["[batch] WARNING may run out of memory", "ValueError: x"]})
    check(r == "prepare failed (exit 1): ValueError: x", "a log line about memory is not an out-of-memory: " + r)


def run(work: Path | None = None, keep: bool = False) -> int:
    root = Path(work) if work else Path(tempfile.mkdtemp(prefix="soda-resume-"))
    if str(APP) not in sys.path:
        sys.path.insert(0, str(APP))
    try:
        for sub in ("inputs", "cache", "export", "logs"):
            shutil.rmtree(root / sub, ignore_errors=True)
        tif = make_input(root)
        cm, em = root / "cache" / DS / "manifest.json", root / "export" / DS / "manifest.json"

        print("1. batch on a new dataset", flush=True)
        d = run_batch(root)
        check(d["status"] == "ok", f"status ok ({d['status']}: {d['reason']})")
        check(cm.is_file() and em.is_file(), "both manifest.json written")

        print("2. export fails part-way", flush=True)
        export_failing(root)
        check(cm.is_file(), "the cache manifest.json is still there (the cache is whole)")
        check(not em.exists(), "no export manifest.json left behind")

        print("3. batch again: the export is done again", flush=True)
        d = run_batch(root)
        check(d["status"] == "ok", f"status ok ({d['status']}: {d['reason']})")
        check(d["steps"]["prepare"].get("skipped") == "already in the cache", "prepare skipped (the cache is whole)")
        check(d["steps"].get("export", {}).get("exit") == 0, "export run again")
        check(em.is_file(), "export manifest.json back")

        print("4. a clean export (its folder deleted, as the docs say) whose prepare fails part-way", flush=True)
        shutil.rmtree(root / "export" / DS)
        prepare_failing(root, tif)
        check(not cm.exists(), "no cache manifest.json left behind")
        check((root / "cache" / DS / "image.zarr").is_dir(), "(image.zarr was rewritten before the failure)")

        print("5. batch again: prepare and export both done again", flush=True)
        d = run_batch(root)
        check(d["status"] == "ok", f"status ok ({d['status']}: {d['reason']})")
        check(d["steps"]["prepare"].get("exit") == 0, "prepare run again (not skipped)")
        check(d["steps"].get("export", {}).get("exit") == 0, "export run again")
        check(cm.is_file() and em.is_file(), "both manifest.json back")

        print("6. failure reasons", flush=True)
        reasons()
        print("PASS", flush=True)
        return 0
    finally:
        if not keep and not work:
            shutil.rmtree(root, ignore_errors=True)


def test_resume():
    assert run() == 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--work", help="folder to work in (default: a new temporary folder, removed at the end)")
    ap.add_argument("--keep", action="store_true", help="keep the temporary folder")
    a = ap.parse_args()
    try:
        sys.exit(run(Path(a.work) if a.work else None, a.keep))
    except AssertionError as e:
        print(f"FAILED: {e}")
        sys.exit(1)
