#!/usr/bin/env python
r"""Quick check that SODA's Python packages import and have the right versions.

    python tools/check_install.py

install.ps1 and install.sh run it at the end. It imports every package the tools need,
prints its version, and checks the version rules the code depends on (zarr 2.x, tifffile
older than 2025, Python 3.12 or 3.13). It takes about 5-10 seconds, up to a minute right
after a fresh install (Python compiles every package once, and umap-learn some code with
numba). Exit code 0 means everything is fine.
"""
import importlib
import sys
import time

# (module to import, package name people know it by, what it is for)
PACKAGES = [
    ("numpy", "numpy", "core"),
    ("scipy", "scipy", "core"),
    ("tifffile", "tifffile", "core"),
    ("imagecodecs", "imagecodecs", "core"),
    ("zarr", "zarr", "core"),
    ("numcodecs", "numcodecs", "core"),
    ("PIL", "pillow", "core"),
    ("h5py", "h5py", "core"),
    ("pyarrow", "pyarrow", "core"),
    ("cryptography", "cryptography", "core"),
    ("fastapi", "fastapi", "tile server"),
    ("uvicorn", "uvicorn", "tile server"),
    ("sklearn", "scikit-learn", "UMAP and clusters"),
    ("umap", "umap-learn", "UMAP and clusters"),
    ("igraph", "igraph", "UMAP and clusters"),
    ("leidenalg", "leidenalg", "UMAP and clusters"),
]


def major(version):
    try:
        return int(str(version).split(".")[0])
    except ValueError:
        return None


def main():
    t0 = time.time()
    problems = []
    py = sys.version_info
    print(f"Python {py.major}.{py.minor}.{py.micro} at {sys.executable}")
    if (py.major, py.minor) not in ((3, 12), (3, 13)):
        problems.append(f"Python {py.major}.{py.minor}: the pinned packages are for Python 3.12 or 3.13")

    versions = {}
    for module, name, what in PACKAGES:
        try:
            m = importlib.import_module(module)
        except Exception as e:  # ImportError, or a broken binary (OSError, ...)
            print(f"  MISSING  {name:14s} ({what}): {type(e).__name__}: {e}")
            problems.append(f"{name} does not import")
            continue
        versions[module] = getattr(m, "__version__", "")
        print(f"  ok       {name:14s} {versions[module]}")

    if "zarr" in versions and (major(versions["zarr"]) or 0) >= 3:
        problems.append(f"zarr {versions['zarr']}: must be 2.x, images cannot be written with zarr 3")
    if "tifffile" in versions and (major(versions["tifffile"]) or 0) >= 2025:
        problems.append(f"tifffile {versions['tifffile']}: must be older than 2025 (newer needs zarr 3)")
    if "numpy" in versions and (major(versions["numpy"]) or 0) < 2:
        problems.append(f"numpy {versions['numpy']}: must be 2.x")

    print(f"  ({time.time() - t0:.1f} s)")
    if problems:
        print("\nNOT READY:")
        for p in problems:
            print("  -", p)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
