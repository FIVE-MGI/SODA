r"""Time one SODA preparing step (prepare.py <kind> / auto / export / rnaseq ...) and record its peak memory.

    python tools/perf/bench_prepare.py --label visiumhd_pub --platform "Visium HD" --step export \
        --csv results.csv --logdir logs [--input PATH ...] [--measure DIR] [--cond "dest=local"] \
        -- <python> web/SpatialViewer/prepare.py export --id x --out cache --dest export

Everything after `--` is the command, run as is (its own Python, its own working directory: --cwd).
While it runs, the script samples every 0.5 s the resident memory (RSS) of the command and all of its
child processes, and at the end it records:
    wall_s         wall-clock time of the whole command
    peak_rss_mb    the highest sum of RSS over the process tree seen in a sample
    peak_wset_mb   Windows only: the largest peak working set any one process of the tree reported
                   (the OS keeps this peak itself, so it catches spikes shorter than the sampling interval)
    cpu_s          user + system CPU seconds of the tree (at the last sample)
    input_bytes / input_files     the --input paths (files or folders, summed)
    output_bytes / output_files   the --measure folder after the run (the export, or the cache entry)
The log (<logdir>/<label>_<step>.log, or _2, _3 ... when that name is taken) has each output line prefixed with seconds since the start;
<logdir>/<label>_<step>.json has the summary and <logdir>/<label>_<step>_mem.csv the samples.
Rows go to --csv in the long format of the performance report (one row per measured value):
    date,category,dataset,platform,step,condition,metric,value,unit,rep,cpu,notes

CPU affinity: this script does not set any. Start it under the affinity you want and the command
inherits it, e.g. on Windows (logical CPUs 16-23 = the E-cores of an i7-13700KF):
    cmd /c "start /affinity FF0000 /wait /b python tools/perf/bench_prepare.py ..."
The affinity the command actually had is recorded in the JSON and in the `cpu` column.

Needs psutil in the Python that runs this script (not in the one that runs the command).
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import os
import platform
import subprocess
import sys
import threading
import time

import psutil


def tree_size(paths):
    """(bytes, files) of files and folders, recursively."""
    b = n = 0
    for p in paths or []:
        if os.path.isfile(p):
            b += os.path.getsize(p); n += 1
            continue
        for root, _dirs, files in os.walk(p):
            for f in files:
                try:
                    b += os.path.getsize(os.path.join(root, f)); n += 1
                except OSError:
                    pass
    return b, n


def cpu_name():
    if sys.platform == "win32":
        try:
            out = subprocess.run(["powershell", "-NoProfile", "-Command", "(Get-CimInstance Win32_Processor).Name"],
                                 capture_output=True, text=True, timeout=30).stdout.strip()
            if out:
                return out
        except Exception:
            pass
    return platform.processor() or platform.machine()


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if "--" not in argv:
        sys.exit("give the command after --")
    k = argv.index("--")
    own, cmd = argv[:k], argv[k + 1:]
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--label", required=True, help="dataset label for the report (no sample IDs)")
    ap.add_argument("--platform", default="")
    ap.add_argument("--step", required=True, help="the step timed, for the report: prepare, export, ...")
    ap.add_argument("--cond", default="", help="condition text for the CSV, e.g. 'cache=D: dest=C:'")
    ap.add_argument("--csv", required=True)
    ap.add_argument("--logdir", required=True)
    ap.add_argument("--cwd", default=None)
    ap.add_argument("--input", action="append", default=[], help="input file/folder (repeatable), for its size")
    ap.add_argument("--measure", action="append", default=[], help="output folder to size after the run (repeatable)")
    ap.add_argument("--notes", default="")
    ap.add_argument("--interval", type=float, default=0.5)
    a = ap.parse_args(own)
    os.makedirs(a.logdir, exist_ok=True)
    prefix = base = os.path.join(a.logdir, f"{a.label}_{a.step}")
    k = 1
    while os.path.exists(prefix + ".json") or os.path.exists(prefix + ".log"):   # never overwrite an earlier run
        k += 1; prefix = f"{base}_{k}"

    in_b, in_n = tree_size(a.input)
    env = dict(os.environ, PYTHONUNBUFFERED="1")
    t0 = time.time()
    p = subprocess.Popen(cmd, cwd=a.cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                         encoding="utf-8", errors="replace", env=env)
    root = psutil.Process(p.pid)
    try:
        affinity = root.cpu_affinity()
    except Exception:
        affinity = None
    samples, stop = [], threading.Event()
    peak = {"rss": 0.0, "wset": 0.0, "cpu": 0.0}
    seen = {}

    def sample():
        while not stop.is_set():
            try:
                procs = [root] + root.children(recursive=True)
            except psutil.Error:
                break
            rss = 0.0; cpu = 0.0; wset = 0.0
            for q in procs:
                try:
                    mi = q.memory_info(); ct = q.cpu_times()
                except psutil.Error:
                    continue
                rss += mi.rss
                wset = max(wset, getattr(mi, "peak_wset", mi.rss))
                seen[q.pid] = ct.user + ct.system
            cpu = sum(seen.values())
            t = time.time() - t0
            samples.append((round(t, 1), round(rss / 2**20, 1), round(wset / 2**20, 1), round(cpu, 1), len(procs)))
            peak["rss"] = max(peak["rss"], rss / 2**20); peak["wset"] = max(peak["wset"], wset / 2**20); peak["cpu"] = cpu
            stop.wait(a.interval)

    th = threading.Thread(target=sample, daemon=True); th.start()
    lines = []
    with open(prefix + ".log", "w", encoding="utf-8") as log:
        log.write("# " + subprocess.list2cmdline(cmd) + "\n")
        log.write(f"# affinity {affinity}\n")
        for line in p.stdout:
            t = time.time() - t0
            log.write(f"{t:9.1f}  {line}"); log.flush()
            lines.append((round(t, 1), line.rstrip()[:300]))
    rc = p.wait()
    wall = time.time() - t0
    stop.set(); th.join(5)
    out_b, out_n = tree_size(a.measure)

    with open(prefix + "_mem.csv", "w", encoding="utf-8") as f:
        f.write("t_s,rss_mb,peak_wset_mb,cpu_s,n_procs\n")
        for s in samples:
            f.write(",".join(map(str, s)) + "\n")
    ncpu = len(affinity) if affinity else os.cpu_count()
    cpu = f"{cpu_name()}; {ncpu} of {os.cpu_count()} logical CPUs" + (f" ({min(affinity)}-{max(affinity)})" if affinity else "")
    summary = {"label": a.label, "platform": a.platform, "step": a.step, "condition": a.cond, "cmd": cmd, "cwd": a.cwd,
               "returncode": rc, "wall_s": round(wall, 1), "peak_rss_mb": round(peak["rss"], 1),
               "peak_wset_mb": round(peak["wset"], 1), "cpu_s": round(peak["cpu"], 1), "affinity": affinity, "cpu": cpu,
               "input": a.input, "input_bytes": in_b, "input_files": in_n,
               "measure": a.measure, "output_bytes": out_b, "output_files": out_n,
               "started": dt.datetime.fromtimestamp(t0).isoformat(timespec="seconds"),
               "host_python": sys.version.split()[0], "lines": lines[-60:]}
    with open(prefix + ".json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=1)

    rows = [("wall_s", round(wall, 1), "s"), ("peak_rss_mb", round(peak["rss"], 1), "MB"),
            ("peak_wset_mb", round(peak["wset"], 1), "MB"), ("cpu_s", round(peak["cpu"], 1), "s")]
    if a.input:
        rows += [("input_bytes", in_b, "B"), ("input_files", in_n, "files")]
    if a.measure:
        rows += [("output_bytes", out_b, "B"), ("output_files", out_n, "files")]
    if rc != 0:
        rows = [("returncode", rc, "")] + rows
    new = not os.path.exists(a.csv)
    with open(a.csv, "a", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        if new:
            w.writerow(["date", "category", "dataset", "platform", "step", "condition", "metric", "value", "unit", "rep", "cpu", "notes"])
        for m, v, u in rows:
            w.writerow([summary["started"], "prepare", a.label, a.platform, a.step, a.cond, m, v, u, 1, cpu, a.notes])
    print(json.dumps({k: v for k, v in summary.items() if k not in ("lines", "cmd")}))
    return rc


if __name__ == "__main__":
    sys.exit(main())
