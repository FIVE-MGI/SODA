r"""Turn the rows written by bench_prepare.py and bench_open.js into Markdown tables (medians over repetitions).

    python tools/perf/summarize.py results.csv [more.csv ...] > tables.md

Columns of the input (long format, one row per measured value):
    date,category,dataset,platform,step,condition,metric,value,unit,rep,cpu,notes
category "prepare" -> one table of wall time, peak memory, output size and file count per dataset and step;
category "open"    -> one table per network condition: first tile, whole first view, bytes, gene, memory,
                      cold and warm, median with min-max over the repetitions;
category "region"  -> Draw a region -> Run UMAP time per region size and Leiden setting.
Sizes on disk and over the network are decimal (1 MB = 10^6 bytes); memory is binary as the OS
reports it (1 MB = 2^20 bytes). Standard library only.
"""
from __future__ import annotations

import csv
import statistics as st
import sys
from collections import OrderedDict, defaultdict


def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def fmt_s(x):
    if x is None:
        return "–"
    return f"{x:.0f} s" if x >= 100 else f"{x:.1f} s"


def fmt_mb(b):
    if b is None:
        return "–"
    mb = b / 1e6
    return f"{mb / 1000:.2f} GB" if mb >= 1000 else (f"{mb:.0f} MB" if mb >= 10 else f"{mb:.2f} MB")


def fmt_mem(mib):
    """memory as the OS reports it, in MiB/GiB (written MB/GB)"""
    return f"{mib / 1024:.1f} GB" if mib >= 1024 else f"{mib:.0f} MB"


def med_rng(vals, f, rng=True):
    vals = [v for v in vals if v is not None]
    if not vals:
        return "–"
    m = st.median(vals)
    if not rng or len(vals) == 1 or min(vals) == max(vals):
        return f(m)
    return f"{f(m)} ({f(min(vals))}–{f(max(vals))})"


def ms_to_s(v):
    return f"{v / 1000:.2f} s" if v < 10000 else f"{v / 1000:.1f} s"


def load(paths):
    rows = []
    for p in paths:
        with open(p, newline="", encoding="utf-8") as f:
            rows += list(csv.DictReader(f))
    return rows


def prepare_table(rows):
    g = OrderedDict()
    for r in rows:
        if r["category"] != "prepare":
            continue
        k = (r["dataset"], r["platform"], r["step"], r["condition"])
        g.setdefault(k, defaultdict(list))[r["metric"]].append(num(r["value"]))
    out = ["| Dataset | Platform | Step | Condition | Wall time | Peak memory (RSS) | Input | Output | Files |",
           "|---|---|---|---|---|---|---|---|---|"]
    for (ds, plat, step, cond), m in g.items():
        def one(k):
            return (m.get(k) or [None])[-1]
        fail = " **(failed)**" if one("returncode") else ""
        rss = one("peak_rss_mb")
        files = one("output_files")
        out.append(f"| {ds} | {plat} | {step}{fail} | {cond} | {fmt_s(one('wall_s'))} "
                   f"| {fmt_mem(rss) if rss else '–'} | {fmt_mb(one('input_bytes')) if one('input_bytes') else '–'} "
                   f"| {fmt_mb(one('output_bytes')) if one('output_bytes') is not None else '–'} "
                   f"| {f'{int(files):,}' if files is not None else '–'} |")
    return "\n".join(out)


def open_tables(rows):
    g = OrderedDict()
    for r in rows:
        if r["category"] != "open":
            continue
        visit, _, net = r["condition"].partition("; ")
        g.setdefault(net, OrderedDict()).setdefault((r["dataset"], r["platform"]), defaultdict(lambda: defaultdict(list)))
        g[net][(r["dataset"], r["platform"])][visit][(r["step"], r["metric"])].append(num(r["value"]))
    parts = []
    for net, dsets in g.items():
        parts.append(f"**{net}** (median over repetitions, min–max in brackets)\n")
        parts.append("| Dataset | Visit | First tile | First view complete | Bytes to first view | Bytes, whole open | Gene (top) | Gene (median) | JS heap | Browser memory |")
        parts.append("|---|---|---|---|---|---|---|---|---|---|")
        for (ds, plat), visits in dsets.items():
            for visit, m in visits.items():
                gt = m.get(("gene (top)", "gene_ms"), [])
                gm = m.get(("gene (median)", "gene_ms"), [])
                bt = m.get(("gene (top)", "bytes_gene"), [])
                parts.append(
                    f"| {ds} | {visit} | {med_rng(m.get(('open', 'first_tile_ms'), []), ms_to_s)} | {med_rng(m.get(('open', 'view_ms'), []), ms_to_s)} "
                    f"| {med_rng(m.get(('open', 'bytes_view'), []), fmt_mb, False)} "
                    f"| {med_rng(m.get(('open', 'bytes_open'), []), fmt_mb, False)} "
                    f"| {med_rng(gt, ms_to_s)}{' / ' + med_rng(bt, fmt_mb, False) if bt else ''} | {med_rng(gm, ms_to_s)} "
                    f"| {med_rng(m.get(('open', 'js_heap_mb'), []), lambda v: f'{v:.0f} MB', False)} "
                    f"| {med_rng(m.get(('open', 'browser_private_mb'), []), lambda v: f'{v:.0f} MB', False)} |")
        parts.append("")
    return "\n".join(parts)


def region_table(rows):
    g = OrderedDict()
    for r in rows:
        if r["category"] != "region":
            continue
        g.setdefault((r["dataset"], r["step"]), defaultdict(list))[r["metric"]].append(num(r["value"]))
    out = ["| Dataset | Region | Run UMAP: started to shown | Browser memory peak |", "|---|---|---|---|"]
    for (ds, step), m in g.items():
        out.append(f"| {ds} | {step} | {med_rng(m.get('region_ms', []), ms_to_s)} | {med_rng(m.get('browser_private_peak_mb', []), lambda v: f'{v:.0f} MB', False)} |")
    return "\n".join(out)


def main(argv=None):
    paths = (argv if argv is not None else sys.argv[1:])
    if not paths:
        sys.exit(__doc__)
    sys.stdout.reconfigure(encoding="utf-8")
    rows = load(paths)
    print("## Preparing\n\n" + prepare_table(rows) + "\n")
    print("## Opening\n\n" + open_tables(rows))
    print("## Draw a region\n\n" + region_table(rows) + "\n")


if __name__ == "__main__":
    main()
