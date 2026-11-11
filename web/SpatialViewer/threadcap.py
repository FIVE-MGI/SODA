"""Cap the BLAS / OpenMP thread pools used by numpy and scipy.

numpy and scipy each bring an OpenBLAS that starts one thread per logical CPU of the
machine, whatever share of them this process may actually use. With the process pinned
to 8 of 24 logical CPUs (an affinity mask, a container CPU limit, a job scheduler), or
with the export's own thread pools busy, those threads spend their time waiting on each
other: the 2000 x 2000 eigensolve in sectionumap._pca took 192 s instead of about 1.5 s.

    import threadcap; threadcap.apply()      # before numpy is first imported

sets OPENBLAS_NUM_THREADS, OMP_NUM_THREADS, MKL_NUM_THREADS (and the other names BLAS
builds read) to min(8, CPUs this process may use). The libraries read them once, when they
load, so this has to run first; child processes inherit them through the environment.

Override: SODA_THREADS=<n> in the environment, or `prepare.py ... --threads <n>` (which
sets SODA_THREADS before anything is imported). A variable that is already set (e.g.
OMP_NUM_THREADS=4 by a job scheduler) is left alone unless SODA_THREADS / --threads is
given, in which case that number wins for all of them.

When numpy was imported before apply() could run (a notebook, another program importing
export.py), the environment is too late: `limit()` caps the pools that are already
loaded for the duration of a `with` block, through threadpoolctl (installed with
scikit-learn, which umap-learn needs). sectionumap wraps its PCA in it.

The cap changes how many threads work on a product or a factorisation, not what is
computed. See tests/smoke_test.py (check_threadcap) for the output comparison.

Only the standard library is imported here, so importing this module never loads numpy.
"""
from __future__ import annotations

import contextlib
import os
import sys

ENV = "SODA_THREADS"
DEFAULT_MAX = 8
# every name a BLAS / OpenMP build is known to read at load time
VARS = ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
        "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS")


def usable_cpus() -> int:
    """Logical CPUs this process may run on: the affinity mask where there is one
    (os.cpu_count() reports the whole machine even when the process is pinned to a few)."""
    try:
        return max(1, len(os.sched_getaffinity(0)))            # Linux
    except (AttributeError, OSError):
        pass
    if sys.platform == "win32":
        try:
            import ctypes
            from ctypes import wintypes
            k32 = ctypes.WinDLL("kernel32", use_last_error=True)
            k32.GetCurrentProcess.restype = wintypes.HANDLE
            k32.GetProcessAffinityMask.argtypes = [wintypes.HANDLE, ctypes.POINTER(ctypes.c_size_t),
                                                   ctypes.POINTER(ctypes.c_size_t)]
            k32.GetProcessAffinityMask.restype = wintypes.BOOL
            proc, system = ctypes.c_size_t(), ctypes.c_size_t()
            if k32.GetProcessAffinityMask(k32.GetCurrentProcess(), ctypes.byref(proc), ctypes.byref(system)):
                n = bin(proc.value).count("1")
                if n:
                    return n
        except Exception:                                       # no kernel32 (Wine, odd builds): fall through
            pass
    return max(1, os.cpu_count() or 1)


def _parse(text) -> int | None:
    try:
        n = int(str(text).strip())
    except (TypeError, ValueError):
        return None
    return n if n >= 1 else None


def requested() -> int | None:
    """The user's own number (SODA_THREADS), or None."""
    return _parse(os.environ.get(ENV))


def cap() -> int:
    """Threads the BLAS / OpenMP pools get: SODA_THREADS if set, else min(8, usable CPUs)."""
    return requested() or min(DEFAULT_MAX, usable_cpus())


def from_argv(argv) -> int | None:
    """`--threads N` / `--threads=N` read straight from the command line: argparse runs
    after the imports, which is too late for the environment variables."""
    for i, a in enumerate(argv):
        if a == "--threads" and i + 1 < len(argv):
            return _parse(argv[i + 1])
        if a.startswith("--threads="):
            return _parse(a.split("=", 1)[1])
    return None


def apply(n: int | None = None) -> int:
    """Set the thread-count variables for this process and its children; returns the cap.

    `n` (from --threads) or SODA_THREADS overrides everything. Otherwise variables the
    user already set are kept. Call before numpy is imported: a BLAS that is already
    loaded keeps the thread count it started with (use limit() for that case).
    """
    if n is not None and _parse(n):
        os.environ[ENV] = str(_parse(n))
    force = requested() is not None
    k = cap()
    for v in VARS:
        if force or not os.environ.get(v):
            os.environ[v] = str(k)
    return k


@contextlib.contextmanager
def limit(n: int | None = None):
    """Cap the BLAS / OpenMP pools that are already loaded while the block runs.

    Does nothing when they are within the cap already (the usual case after apply())
    or when threadpoolctl is not installed.
    """
    k = n or cap()
    try:
        from threadpoolctl import threadpool_info, threadpool_limits
    except ImportError:
        yield k
        return
    if all(p.get("num_threads", 1) <= k for p in threadpool_info()):
        yield k
        return
    with threadpool_limits(limits=k):
        yield k


def describe() -> str:
    """One line for the export log: the cap and where it came from."""
    src = f"{ENV}" if requested() else f"min({DEFAULT_MAX}, {usable_cpus()} usable CPUs)"
    return f"{cap()} BLAS/OpenMP threads ({src})"
