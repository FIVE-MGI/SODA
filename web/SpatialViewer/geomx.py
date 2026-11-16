r"""GeoMx DSP (NanoString / Bruker) -> a SODA dataset: scan image, ROI outlines, one "rois" layer.

    python prepare.py geomx --dcc <folder or .zip of *.dcc> --pkc <file.pkc ...> \
           --annotation <annotations.xlsx | .csv | *_LabWorksheet.txt> [--annotation <more>] \
           [--image <scan.ome.tiff> ...] --id study1 --out <cache folder>

GeoMx is not cells or bins. A person draws a few dozen regions of interest (ROIs) on a scan,
sometimes split into segments (AOIs: e.g. PanCK+ and PanCK- inside one ROI); each AOI goes
into one well and gives ONE count vector. So the layer written here (key "rois",
kind "rois") has one row per AOI:

    position      centroid of the AOI's mask in the scan's OME-XML (else of the ROI polygon),
                  moved inside the polygon if it falls outside; full-resolution pixels.
                  Rows are in the scans' left-to-right order, then annotation order
    counts        DCC probe counts collapsed to targets with the PKC, as GeomxTools
                  aggregateCounts: one probe -> its count; several -> their geometric mean,
                  0.5 added to all of them first if any is 0. Negative-probe targets are not
                  genes: they give the background (NegGeoMean) and the limit of
                  quantification LOQ = max(2, NegGeoMean * NegGeoSD^2) per AOI
    barcodes      the Sample_ID (= DCC file name); cache only, never exported
    overlay       rois.geojson: one polygon per ROI, the dataset's region overlay
    groupings     every categorical annotation column, de-identified as userclusters.py does,
                  written at export as named clusterings with "source": "user"
    q3            GeomxTools' Q3 normalisation (Q3_RULE): each AOI's Q3 in the cache, the normFactors
                  per exported row in the layer's manifest entry ("q3") at export; the viewer
                  multiplies a count by 1 / normFactor (its default for plots, volcano, colours)

Several --image scans are put side by side on one canvas (a study is usually 2-4 slides and
its comparisons run across them). Without an image (or with --layout grid) the AOIs are drawn
on a labelled grid so the dataset still opens; those positions are not real.

What is NOT written for such a layer, on purpose: the whole-section UMAP, Leiden clusters,
suggested cell types and the variable-gene selection. A few dozen hand-picked regions are not
a sample of the tissue, and the viewer hides those tools for kind "rois". The hvg/ chunks the
volcano reads hold every target instead (export_layer).

Annotation files: the GeoMx "annotation" workbook used with GeomxTools (one sheet, a header
row with Sample_ID, slide name, scan name, panel, roi, segment, aoi, area, then the study's
own columns), a CSV / TSV of the same, the *_LabWorksheet.txt of the readout package (the same
table after a preamble), and a DSP Data Analysis suite export with a SegmentProperties sheet
(no Sample_ID: joined by scan, ROI and segment). .xlsx is read with the standard library only
(zip + XML): no openpyxl needed.

Needs only numpy, scipy, zarr, tifffile and Pillow (all in requirements.txt).
"""
from __future__ import annotations

import base64
import csv
import io
import json
import math
import re
import time
import zipfile
import xml.etree.ElementTree as ET
from datetime import datetime
from pathlib import Path

import numpy as np

import deid

LAYER_KEY = "rois"
OVERLAY_FILE = "rois.geojson"
CACHE_DIR = "geomx"                 # <cache>/<id>/geomx/rows.json: annotation values and QC per row (cache only)
NTC_NAME = "notemplatecontrol"
MIN_MAPPED = 0.8                    # share of a study's DCC counts the PKC(s) must explain
MAX_MISSING_DCC = 0.2               # share of annotation rows allowed to have no DCC file
MAX_UNPLACED = 0.5                  # share of a scan's rows allowed to have no ROI in its OME-XML
LOQ_MIN, LOQ_SD = 2.0, 2.0          # GeomxTools: LOQ = max(2, NegGeoMean * NegGeoSD ^ 2)
ZERO_SHIFT = 0.5                    # GeomxTools thresholdValues: + 0.5 on a target's probes if any is 0
MOSAIC_GAP = 512                    # px of empty canvas between scans
GRID_CELL, GRID_HEAD, GRID_COLS = 256, 96, 12
DISC = 0.6                          # drawn disc: this share of the diameter of a circle with the median AOI area,
                                    # so that it sits inside most outlines (ROIs are rarely round)
COLLAPSE_RULE = ("probes -> targets as GeomxTools aggregateCounts: a target with one probe keeps its count; "
                 "with several, the geometric mean of its probes, 0.5 added to all of them first when any is 0")
LOQ_RULE = "LOQ = max(2, NegGeoMean * NegGeoSD^2) of the module's negative probes, per AOI"
# Q3 (upper-quartile) normalisation, as GeomxTools defines it. Source: GeomxTools (Bioconductor,
# R/NanoStringGeoMxSet-normalize.R), quantileNorm, called by normalize(object, norm_method = "quant",
# desiredQuantile = 0.75):
#     qs <- apply(exprs(object), 2, function(x) stats::quantile(x, desiredQuantile))
#     normFactors <- qs / ngeoMean(qs);   normalised = counts / normFactors   (per segment)
# ngeoMean (R/utils.R) = geometric mean after thresholdValues (+0.5 on every value when the minimum
# is under 0.5): _geo below. stats::quantile's default (type 7) = numpy's linear percentile.
# exprs() holds whatever targets are left at that point; the GeoMx workflow vignette (Bioconductor
# GeoMxWorkflows, "GeomxTools_RNA-NGS_Analysis") first keeps the targets above the LOQ in >= 10% of
# the segments plus the negative-probe targets ("DetectionRate >= 0.1 | TargetName %in% neg_probes"),
# and calls Q3 "typically the preferred normalization strategy for most DSP-NGS RNA studies".
# Without negative probes (no LOQ) every target is used.
Q3_QUANTILE = 75.0
Q3_MIN_DETECTION = 0.10
Q3_RULE = ("Q3 normalisation as GeomxTools normalize(norm_method = 'quant', desiredQuantile = 0.75): per AOI, Q3 = the "
           "75th percentile (R quantile type 7) of its target counts; normFactor = Q3 / geometric mean of the Q3s of the "
           "AOIs (GeomxTools ngeoMean); value = count / normFactor, so values stay count-like. Targets: those above the "
           "LOQ in >= 10% of the AOIs plus the negative-probe targets, as the GeoMx workflow vignette filters before "
           "normalising (every target when the panel has no negative probes)")
NO_EMBED_WHY = ("a few dozen hand-drawn regions: no UMAP or Leiden clusters are computed for an ROI layer; "
                "the groupings are the annotation columns")


# ---------------------------------------------------------------------------------------
# command line (prepare.py geomx)
# ---------------------------------------------------------------------------------------
def add_args(p):
    p.add_argument("--dcc", required=True, help="folder (searched recursively) or .zip of the *.dcc files; "
                                                 "files of other scans are ignored")
    p.add_argument("--pkc", required=True, nargs="+", help="probe kit file(s) *.pkc (or a folder / .zip of them)")
    p.add_argument("--annotation", required=True, action="append",
                   help="annotation workbook (.xlsx), CSV / TSV, or the *_LabWorksheet.txt of the readout package. "
                        "Repeatable: more files (e.g. a DSP Data Analysis export with a SegmentProperties sheet) "
                        "are joined by Sample_ID, else by scan + ROI + segment")
    p.add_argument("--sheet", help="sheet name in the .xlsx (default: the first sheet with a Sample_ID header)")
    p.add_argument("--image", action="append", default=[],
                   help="scan OME-TIFF exported from the GeoMx (its OME-XML holds the ROI polygons). Repeatable: "
                        "several scans go side by side. 'ScanName=path' if the file's own scan name is not the "
                        "one in the annotation")
    p.add_argument("--layout", choices=["auto", "grid"], default="auto",
                   help="grid = draw the AOIs on a labelled grid even when an image is given (auto: the scan "
                        "image when there is one, else the grid)")
    p.add_argument("--keep-slide-names", action="store_true",
                   help="keep the slide / scan names as group names (through the sample-ID pattern only). "
                        "Default: they are replaced by Slide 1, Slide 2 ..., as these free-text names often "
                        "carry a patient or block number")


# ---------------------------------------------------------------------------------------
# tables: .xlsx (standard library), CSV / TSV, LabWorksheet.txt
# ---------------------------------------------------------------------------------------
_X = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
_XR = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"
_EXCEL_TEXT = re.compile(r'^="(.*)"$')          # ="001": how the GeoMx writes a number that must stay text


def _clean(v) -> str:
    s = "" if v is None else str(v).strip()
    m = _EXCEL_TEXT.match(s)
    return m.group(1).strip() if m else s


def _col_index(ref: str) -> int:
    n = 0
    for ch in ref:
        if not ch.isalpha():
            break
        n = n * 26 + (ord(ch.upper()) - 64)
    return n - 1


def xlsx_sheets(path: Path) -> list[tuple[str, str]]:
    """[(sheet name, member path)] of a workbook, in workbook order."""
    with zipfile.ZipFile(path) as z:
        wb = ET.fromstring(z.read("xl/workbook.xml"))
        rels = {r.get("Id"): r.get("Target") for r in ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))}
        out = []
        for s in wb.iter(_X + "sheet"):
            t = rels.get(s.get(_XR + "id"), "")
            t = t.lstrip("/") if t.startswith("/") else "xl/" + t
            out.append((s.get("name") or "", t))
    return out


def xlsx_rows(path: Path, member: str, max_rows: int | None = None) -> list[list[str]]:
    """Rows of one sheet as lists of strings (shared and inline strings resolved, gaps filled).
    Streamed, so one sheet of a large workbook can be read without the others."""
    with zipfile.ZipFile(path) as z:
        strings = []
        if "xl/sharedStrings.xml" in z.namelist():
            with z.open("xl/sharedStrings.xml") as f:
                for _, el in ET.iterparse(f):
                    if el.tag == _X + "si":
                        # plain text or rich-text runs; the phonetic guides (rPh) some locales add are not text
                        strings.append("".join((r.text or "") if r.tag == _X + "t" else
                                               "".join(t.text or "" for t in r.iter(_X + "t"))
                                               for r in el if r.tag in (_X + "t", _X + "r")))
                        el.clear()
        rows = []
        with z.open(member) as f:
            for _, el in ET.iterparse(f):
                if el.tag != _X + "row":
                    continue
                cells = []
                for c in el.iter(_X + "c"):
                    i = _col_index(c.get("r") or "") if c.get("r") else len(cells)
                    t, v = c.get("t"), c.find(_X + "v")
                    if t == "s":
                        val = strings[int(v.text)] if v is not None and v.text else ""
                    elif t == "inlineStr":
                        val = "".join(x.text or "" for x in c.iter(_X + "t"))
                    else:
                        val = v.text if v is not None and v.text is not None else ""
                    while len(cells) < i:
                        cells.append("")
                    cells.append(val)
                rows.append(cells)
                el.clear()
                if max_rows and len(rows) >= max_rows:
                    break
    return rows


def delimited_rows(path: Path) -> list[list[str]]:
    text = Path(path).read_text(encoding="utf-8-sig", errors="replace")
    head = text[:20000]
    delim = "\t" if head.count("\t") > head.count(",") else ","
    return [r for r in csv.reader(io.StringIO(text), delimiter=delim)]


def _norm(h: str) -> str:
    return re.sub(r"[^a-z0-9]", "", str(h).lower())


# role -> header spellings (normalised): the GeomxTools annotation sheet, the lab worksheet, DSPDA's SegmentProperties
ROLES = {
    "sample_id": ("sampleid", "dccfile", "dccfilename"),
    "slide": ("slidename", "slide"),
    "scan": ("scanname", "scanlabel", "scan"),
    "roi": ("roi", "roilabel", "roiname"),
    "segment": ("segment", "segmentlabel", "segmentname"),
    "aoi": ("aoi", "aoiname"),
    "area": ("area", "aoisurfacearea", "aoiarea"),
    "nuclei": ("nuclei", "aoinucleicount", "nucleicount", "nucleus", "aoinuclei"),
    "panel": ("panel",),
}


def _roles(header: list[str]) -> dict[str, int]:
    n = [_norm(h) for h in header]
    out = {}
    for role, names in ROLES.items():
        for nm in names:
            if nm in n:
                out[role] = n.index(nm)
                break
    return out


def _header_row(rows: list[list[str]]) -> int | None:
    """Index of the header: the first row naming Sample_ID, else ROI + segment columns."""
    for want in (("sample_id",), ("roi", "segment")):
        for i, r in enumerate(rows[:400]):
            ro = _roles([_clean(c) for c in r])
            if all(w in ro for w in want) and len([c for c in r if _clean(c)]) >= 2:
                return i
    return None


def read_table(path: Path, sheet: str | None = None, log=print) -> dict:
    """{"name", "header": [...], "rows": [[str]], "roles": {role: column}} of one annotation file."""
    path = Path(path)
    if not path.exists():
        raise SystemExit(f"[geomx] annotation file not found: {path}")
    if path.suffix.lower() in (".xlsx", ".xlsm"):
        sheets = xlsx_sheets(path)
        if sheet:
            pick = [s for s in sheets if s[0].lower() == sheet.lower()]
            if not pick:
                raise SystemExit(f"[geomx] {path.name}: no sheet '{sheet}' (sheets: {[s[0] for s in sheets]})")
        else:
            # a sheet with Sample_ID first, then one with ROI + segment columns; only the top is read to decide
            scored = []
            for s in sheets:
                top = xlsx_rows(path, s[1], max_rows=60)
                h = _header_row(top)
                if h is not None:
                    scored.append((0 if "sample_id" in _roles([_clean(c) for c in top[h]]) else 1, s))
            if not scored:
                raise SystemExit(f"[geomx] {path.name}: no sheet has a header row with Sample_ID (or ROI and "
                                 f"segment columns). Sheets: {[s[0] for s in sheets]}")
            pick = [sorted(scored, key=lambda t: t[0])[0][1]]
        rows = xlsx_rows(path, pick[0][1])
        what = f"{path.name} [{pick[0][0]}]"
    elif path.suffix.lower() == ".xls":
        raise SystemExit(f"[geomx] {path.name}: the old .xls format is not read. Save it as .xlsx or .csv.")
    else:
        rows = delimited_rows(path)
        what = path.name
    h = _header_row(rows)
    if h is None:
        raise SystemExit(f"[geomx] {what}: no header row with Sample_ID (or ROI and segment columns) found")
    header = [_clean(c) for c in rows[h]]
    while header and not header[-1]:
        header.pop()
    body = []
    for r in rows[h + 1:]:
        r = [_clean(c) for c in r[:len(header)]] + [""] * max(0, len(header) - len(r))
        if any(r):
            body.append(r)
    roles = _roles(header)
    log(f"[geomx] annotation {what}: {len(body)} rows, {len(header)} columns"
        + ("" if "sample_id" in roles else " (no Sample_ID: joined by scan + ROI + segment)"))
    return {"name": what, "header": header, "rows": body, "roles": roles}


def _roi_key(v: str) -> str:
    v = str(v).strip()
    return str(int(v)) if v.isdigit() else v.lower()


def _join_key(t: dict, r: list[str]):
    ro = t["roles"]
    if not all(k in ro for k in ("scan", "roi", "segment")):
        return None
    return (r[ro["scan"]].strip().lower(), _roi_key(r[ro["roi"]]), r[ro["segment"]].strip().lower())


def read_annotation(paths, sheet: str | None = None, log=print) -> dict:
    """The annotation as one table: the file with Sample_ID is the base, the others add columns."""
    tables = [read_table(Path(p), sheet if i == 0 else None, log) for i, p in enumerate(paths)]
    base = next((t for t in tables if "sample_id" in t["roles"]), None)
    if base is None:
        raise SystemExit("[geomx] no annotation file has a Sample_ID column: it is what names each AOI's DCC file. "
                         "Give the annotation workbook or the *_LabWorksheet.txt of the readout package.")
    sid = base["roles"]["sample_id"]
    seen, dup = set(), []
    for r in base["rows"]:
        if r[sid] in seen:
            dup.append(r[sid])
        seen.add(r[sid])
    if dup:
        raise SystemExit(f"[geomx] {base['name']}: Sample_ID appears more than once: {dup[:5]}")
    for t in tables:
        if t is base:
            continue
        have = {_norm(h) for h in base["header"]}
        role_of = {v: k for k, v in t["roles"].items()}
        # not a second copy of a column the base already has under another spelling (ScanLabel = scan name)
        add = [j for j, h in enumerate(t["header"]) if h and _norm(h) not in have
               and role_of.get(j) not in base["roles"]]
        if "sample_id" in t["roles"]:
            idx = {r[t["roles"]["sample_id"]]: r for r in t["rows"]}
            find = lambda r: idx.get(r[sid])
            how = "Sample_ID"
        else:
            idx = {_join_key(t, r): r for r in t["rows"]}
            if None in idx or not all(k in base["roles"] for k in ("scan", "roi", "segment")):
                raise SystemExit(f"[geomx] {t['name']}: needs Sample_ID, or scan, ROI and segment columns in both files, "
                                 f"to be joined")
            find = lambda r: idx.get(_join_key(base, r))
            how = "scan + ROI + segment"
        hit = 0
        for r in base["rows"]:
            o = find(r)
            hit += o is not None
            r.extend([(o[j] if o is not None else "") for j in add])
        base["header"].extend(t["header"][j] for j in add)
        base["roles"] = {**_roles(base["header"]), **base["roles"]}
        log(f"[geomx] annotation {t['name']}: {len(add)} columns added to {hit} of {len(base['rows'])} rows (by {how})")
    return base


# ---------------------------------------------------------------------------------------
# DCC files
# ---------------------------------------------------------------------------------------
class DccStore:
    """The *.dcc files of a folder (recursive) or a .zip, by file stem (= Sample_ID)."""

    def __init__(self, src: Path):
        self.src = Path(src)
        self._zip = None
        self.files: dict[str, object] = {}
        if self.src.is_dir():
            for p in sorted(self.src.rglob("*.dcc")):
                self.files.setdefault(p.stem, p)
        elif self.src.suffix.lower() == ".zip" and self.src.exists():
            self._zip = zipfile.ZipFile(self.src)
            for n in sorted(self._zip.namelist()):
                if n.lower().endswith(".dcc"):
                    self.files.setdefault(Path(n).stem, n)
        else:
            raise SystemExit(f"[geomx] --dcc: {self.src} is not a folder or a .zip")
        if not self.files:
            raise SystemExit(f"[geomx] no .dcc files under {self.src}")

    def text(self, stem: str) -> str:
        f = self.files[stem]
        if self._zip is not None:
            return self._zip.read(f).decode("utf-8", errors="replace")
        return Path(f).read_text(encoding="utf-8", errors="replace")

    def close(self):
        if self._zip is not None:
            self._zip.close()


def parse_dcc(text: str) -> dict:
    """{"id", "plate", "well", "raw", "trimmed", "stitched", "aligned", "umiQ30", "rtsQ30", "codes": {RTS_ID: count}}."""
    out = {"codes": {}}
    section = None
    num = {"Raw": "raw", "Trimmed": "trimmed", "Stitched": "stitched", "Aligned": "aligned",
           "umiQ30": "umiQ30", "rtsQ30": "rtsQ30"}
    codes = out["codes"]
    for ln in text.splitlines():
        ln = ln.strip()
        if not ln:
            continue
        if ln.startswith("<"):
            section = None if ln.startswith("</") else ln.strip("<>")
            continue
        k, _, v = ln.partition(",")
        if section == "Code_Summary":
            try:
                codes[k] = codes.get(k, 0.0) + float(v)
            except ValueError:
                pass
        elif section == "Scan_Attributes":
            if k == "ID":
                out["id"] = v.strip('"')
            elif k == "Plate_ID":
                out["plate"] = v.strip('"')
            elif k == "Well":
                out["well"] = v.strip('"')
        elif section == "NGS_Processing_Attributes" and k in num:
            try:
                out[num[k]] = float(v)
            except ValueError:
                pass
    return out


# ---------------------------------------------------------------------------------------
# PKC files and the probe -> target collapse
# ---------------------------------------------------------------------------------------
def _pkc_texts(paths) -> list[tuple[str, str]]:
    out = []
    for p in map(Path, paths):
        if p.is_dir():
            out += [(q.name, q.read_text(encoding="utf-8-sig")) for q in sorted(p.rglob("*.pkc"))]
        elif p.suffix.lower() == ".zip" and p.exists():
            with zipfile.ZipFile(p) as z:
                out += [(Path(n).name, z.read(n).decode("utf-8-sig")) for n in sorted(z.namelist())
                        if n.lower().endswith(".pkc")]
        elif p.exists():
            out.append((p.name, p.read_text(encoding="utf-8-sig")))
        else:
            raise SystemExit(f"[geomx] --pkc: {p} not found")
    if not out:
        raise SystemExit("[geomx] --pkc: no .pkc file found")
    return out


def read_pkcs(paths, log=print) -> dict:
    """The probe kits as flat arrays.

    {"modules": [{"name", "version", "targets", "probes", "negativeProbes"}],
     "rts": {RTS_ID: probe index}, "probe_target": int array (probe -> target),
     "targets": [{"name", "codeClass", "module", "negative"}]}
    A target name used by two modules is kept once per module (the later one as "NAME (module)").
    """
    modules, rts, probe_target, targets = [], {}, [], []
    names = set()
    for fname, text in _pkc_texts(paths):
        try:
            j = json.loads(text)
        except ValueError as e:
            raise SystemExit(f"[geomx] {fname} is not a PKC file (JSON expected): {e}")
        if "Targets" not in j:
            raise SystemExit(f"[geomx] {fname}: no 'Targets' in the PKC")
        mod = str(j.get("Name") or Path(fname).stem)
        if any(m["name"] == mod for m in modules):
            raise SystemExit(f"[geomx] two PKC files are the panel '{mod}': give one version of each panel")
        mi = len(modules)
        n_neg = n_pr = 0
        for t in j["Targets"]:
            cc = str(t.get("CodeClass") or "")
            neg = cc.lower().startswith("negative")
            name = str(t.get("DisplayName") or "").strip() or f"target{len(targets)}"
            if name in names:
                name = f"{name} ({mod})"
            names.add(name)
            ti = len(targets)
            targets.append({"name": name, "codeClass": re.sub(r"\d+$", "", cc), "module": mi, "negative": neg})
            for p in t.get("Probes") or []:
                r = str(p.get("RTS_ID") or "")
                if not r:
                    continue
                if r in rts:
                    raise SystemExit(f"[geomx] probe {r} is in two PKC files: give one version of each panel")
                rts[r] = len(probe_target)
                probe_target.append(ti)
                n_pr += 1
                n_neg += neg
        modules.append({"name": mod, "version": str(j.get("Version", "")), "targets": sum(1 for t in targets if t["module"] == mi),
                        "probes": n_pr, "negativeProbes": n_neg})
        log(f"[geomx] PKC {fname}: '{mod}' v{j.get('Version', '?')}, {modules[-1]['targets']:,} targets, "
            f"{n_pr:,} probes ({n_neg} negative)")
    return {"modules": modules, "rts": rts, "probe_target": np.asarray(probe_target, np.int64), "targets": targets}


def _geo(v: np.ndarray, axis: int = 0) -> tuple[np.ndarray, np.ndarray]:
    """Geometric mean and geometric SD (sample SD of the logs) along `axis`, after GeomxTools'
    thresholdValues: 0.5 is added to every value of a vector whose minimum is under 0.5."""
    v = np.asarray(v, np.float64)
    low = v.min(axis=axis, keepdims=True) < ZERO_SHIFT
    lg = np.log(v + ZERO_SHIFT * low)
    gm = np.exp(lg.mean(axis=axis))
    gsd = np.exp(lg.std(axis=axis, ddof=1)) if v.shape[axis] > 1 else np.ones_like(gm)
    return gm, gsd


def collapse(P: np.ndarray, probe_target: np.ndarray, n_targets: int) -> np.ndarray:
    """Probe counts (probes x AOIs) -> target values (targets x AOIs), COLLAPSE_RULE."""
    P = np.asarray(P, np.float64)
    order = np.argsort(probe_target, kind="stable")
    pt = probe_target[order]
    counts = np.bincount(pt, minlength=n_targets)
    T = np.zeros((n_targets, P.shape[1]), np.float64)
    start = np.concatenate([[0], np.cumsum(counts)])
    single = np.flatnonzero(counts == 1)
    T[single] = P[order[start[single]]]
    for t in np.flatnonzero(counts > 1):
        T[t] = _geo(P[order[start[t]:start[t + 1]]], axis=0)[0]
    return T


def q3_values(G: np.ndarray, neg_T: np.ndarray, detected: np.ndarray | None, n: int, log=print) -> dict:
    """Each AOI's Q3 for Q3_RULE (cache only; the factors are made at export, over the AOIs exported).

    G: gene targets x AOIs; neg_T: negative-probe targets x AOIs (collapsed as the genes);
    detected: per gene target, the number of AOIs where it is above its module's LOQ (None: no LOQ).
    """
    if detected is not None and n:
        keep = detected / n >= Q3_MIN_DETECTION
        Q = np.vstack([G[keep], neg_T])
        n_neg = int(len(neg_T))
    else:
        keep = np.ones(len(G), bool)
        Q, n_neg = G, 0
    q = np.percentile(Q, Q3_QUANTILE, axis=0) if len(Q) else np.zeros(n)
    log(f"[geomx] Q3 normalisation: Q3 over {int(keep.sum()):,} of {len(G):,} targets"
        + (f" (above LOQ in >= {Q3_MIN_DETECTION:.0%} of the AOIs) + {n_neg} negative-probe target(s)" if detected is not None
           else " (no negative probes, so no LOQ: every target)")
        + (f"; median Q3 {np.median(q):.2f}" if len(q) else "")
        + (f"; WARNING {int((q <= 0).sum())} AOI(s) with Q3 = 0 cannot be Q3-normalised" if (q <= 0).any() else ""))
    return {"q3": [float(x) for x in q], "genesKept": int(keep.sum()), "genes": int(len(G)), "negativeTargets": n_neg,
            "loqFilter": detected is not None}


def q3_entry(rows: dict, on: np.ndarray, log=print) -> dict | None:
    """The layer's "q3" entry: GeomxTools' normFactors (Q3 / ngeoMean of the Q3s) per exported row,
    in positions.bin order. The browser multiplies a count by 1 / normFactor. None for a cache
    prepared before Q3 was written (the viewer then keeps counts per 10,000)."""
    qn = rows.get("q3norm")
    if not qn or len(qn.get("q3", [])) != rows["n"]:
        log("[geomx] no Q3 values in this cache (prepared by an older version): the plots stay on counts per 10,000. "
            "Run prepare.py geomx again for Q3 normalisation.")
        return None
    q = np.array([np.nan if x is None else x for x in qn["q3"]], np.float64)[on]
    if not len(q) or not np.isfinite(q).all() or not (q > 0).any():
        log("[geomx] Q3 normalisation left out: no AOI has a Q3 above 0")
        return None
    gm = float(_geo(q)[0])                            # GeomxTools ngeoMean(qs): thresholdValues, then the geometric mean
    f = q / gm
    zero = int((q <= 0).sum())
    log(f"[geomx] Q3 normFactors for {len(q)} AOIs: geometric mean of the Q3s {gm:.3f}, factors {f.min():.3f} .. {f.max():.3f}"
        + (f"; {zero} AOI(s) with Q3 = 0 get no factor (left out of the Q3 plots)" if zero else ""))
    return {"method": "Q3", "quantile": Q3_QUANTILE / 100, "normFactors": [float(f"{x:.7g}") for x in f],
            "geoMeanQ3": round(gm, 6), "genesKept": qn.get("genesKept"), "genes": qn.get("genes"),
            "negativeTargets": qn.get("negativeTargets"), "loqFilter": qn.get("loqFilter"),
            "minDetection": Q3_MIN_DETECTION if qn.get("loqFilter") else None, "zeroQ3": zero,
            "default": True, "rule": Q3_RULE}


# ---------------------------------------------------------------------------------------
# geometry: OME-XML ROIs, polygons, masks
# ---------------------------------------------------------------------------------------
def _tag(el) -> str:
    return el.tag.rsplit("}", 1)[-1]


def poly_area_centroid(p: np.ndarray) -> tuple[float, float, float]:
    """(area, cx, cy) of a simple polygon (n, 2); the vertex mean when the area is 0."""
    x, y = p[:, 0], p[:, 1]
    x1, y1 = np.roll(x, -1), np.roll(y, -1)
    cr = x * y1 - x1 * y
    a = cr.sum() / 2.0
    if abs(a) < 1e-9:
        return 0.0, float(x.mean()), float(y.mean())
    return abs(a), float(((x + x1) * cr).sum() / (6 * a)), float(((y + y1) * cr).sum() / (6 * a))


def point_in_poly(x: float, y: float, p: np.ndarray) -> bool:
    xi, yi = p[:, 0], p[:, 1]
    xj, yj = np.roll(xi, -1), np.roll(yi, -1)
    cond = (yi > y) != (yj > y)
    with np.errstate(divide="ignore", invalid="ignore"):
        xint = xi + (xj - xi) * (y - yi) / (yj - yi)
    return bool(np.count_nonzero(cond & (x < xint)) % 2)


def chord(p: np.ndarray, x: float, y: float) -> tuple[float, float]:
    """The horizontal run of the polygon's inside at height y that holds x (else its widest run)."""
    xi, yi = p[:, 0], p[:, 1]
    xj, yj = np.roll(xi, -1), np.roll(yi, -1)
    cond = (yi > y) != (yj > y)
    if not cond.any():
        return float(x), float(x)
    xs = np.sort(xi[cond] + (xj[cond] - xi[cond]) * (y - yi[cond]) / (yj[cond] - yi[cond]))
    runs = [(float(xs[k]), float(xs[k + 1])) for k in range(0, len(xs) - 1, 2)]
    for a, b in runs:
        if a <= x <= b:
            return a, b
    return max(runs, key=lambda r: r[1] - r[0])


def interior_point(p: np.ndarray, x: float, y: float) -> tuple[float, float]:
    """(x, y) if it is inside the polygon, else the middle of the widest inside run at a height
    where the polygon is wide (its centroid's height, then its bounding box's middle)."""
    if point_in_poly(x, y, p):
        return float(x), float(y)
    y0, y1 = float(p[:, 1].min()), float(p[:, 1].max())
    for yy in (y if y0 < y < y1 else None, (y0 + y1) / 2, y0 + (y1 - y0) / 3, y0 + 2 * (y1 - y0) / 3):
        if yy is None:
            continue
        a, b = chord(p, float(p[:, 0].mean()), yy)
        if b > a and point_in_poly((a + b) / 2, yy, p):
            return (a + b) / 2, float(yy)
    return float(p[:, 0].mean()), float(p[:, 1].mean())


def _ellipse(cx, cy, rx, ry, n=64) -> np.ndarray:
    a = np.arange(n) * (2 * math.pi / n)
    return np.column_stack([cx + rx * np.cos(a), cy + ry * np.sin(a)])


def decode_mask(el, bin_el) -> dict | None:
    """Area (pixels) and centroid (image pixels) of an OME Mask: BinData is one bit per pixel,
    row by row, most significant bit first, no padding at the row ends."""
    try:
        w, h = int(round(float(el.get("Width")))), int(round(float(el.get("Height"))))
        x0, y0 = float(el.get("X")), float(el.get("Y"))
        raw = base64.b64decode(bin_el.text or "")
        comp = (bin_el.get("Compression") or "none").lower()
        if comp == "zlib":
            import zlib
            raw = zlib.decompress(raw)
        elif comp == "bzip2":
            import bz2
            raw = bz2.decompress(raw)
        elif comp != "none":
            return None
        bits = np.unpackbits(np.frombuffer(raw, np.uint8), bitorder="big")
        if w <= 0 or h <= 0 or len(bits) < w * h:
            return None
        m = bits[: w * h].reshape(h, w)
        n = int(m.sum())
        if not n:
            return None
        ys, xs = np.nonzero(m)
        return {"areaPx": n, "cx": x0 + float(xs.mean()) + 0.5, "cy": y0 + float(ys.mean()) + 0.5,
                "bbox": (x0, y0, x0 + w, y0 + h)}
    except Exception:
        return None


def ome_xml_of(tiff: Path) -> str:
    """The OME-XML of a scan: the block inside the TIFF, else a .ome.xml file next to it."""
    import tifffile
    xml = None
    try:
        with tifffile.TiffFile(str(tiff)) as tf:
            xml = tf.ome_metadata
    except Exception:
        xml = None
    if not xml or "<ROI" not in xml:
        for side in (Path(str(tiff).rsplit(".ome.", 1)[0] + ".ome.xml"), Path(tiff).with_suffix(".xml")):
            if side.exists() and side != Path(tiff):
                t = side.read_text(encoding="utf-8", errors="replace")
                if "<ROI" in t or not xml:
                    xml = t
                break
    return xml or ""


def read_scan(tiff: Path, log=print) -> dict:
    """Scan name, pixel size, channel names and the ROIs of a GeoMx scan OME-TIFF.

    {"scan", "px": (x µm, y µm) | None, "channels": [...], "size": (W, H) | None,
     "rois": {roi key: {"label", "poly": (n, 2) px, "masks": {segment lower: {"name", "areaPx", "cx", "cy"}}}}}
    """
    xml = ome_xml_of(tiff)
    out = {"scan": None, "px": None, "channels": [], "size": None, "rois": {}}
    if not xml.strip():
        log(f"[geomx] {Path(tiff).name}: no OME-XML in or next to the file")
        return out
    root = ET.fromstring(xml.encode("utf-8"))
    img = next((e for e in root.iter() if _tag(e) == "Image"), None)
    if img is not None:
        out["scan"] = img.get("Name")
        px = next((e for e in img.iter() if _tag(e) == "Pixels"), None)
        if px is not None:
            try:
                sx, sy = float(px.get("PhysicalSizeX")), float(px.get("PhysicalSizeY") or px.get("PhysicalSizeX"))
                unit = (px.get("PhysicalSizeXUnit") or "µm").lower()
                k = {"nm": 1e-3, "mm": 1e3}.get(unit, 1.0)
                out["px"] = (sx * k, sy * k)
            except (TypeError, ValueError):
                pass
            try:
                out["size"] = (int(px.get("SizeX")), int(px.get("SizeY")))
            except (TypeError, ValueError):
                pass
            for ch in px:
                if _tag(ch) == "Channel":
                    fl, nm = (ch.get("Fluor") or "").strip(), (ch.get("Name") or "").strip()
                    out["channels"].append(f"{fl} ({nm})" if fl and nm and fl != nm else (fl or nm or f"ch{len(out['channels'])}"))
    for e in root.iter():          # the "ScanName" comment annotation wins over the image name when both exist
        if _tag(e) == "CommentAnnotation":
            d = next((c.text for c in e if _tag(c) == "Description"), None)
            v = next((c.text for c in e if _tag(c) == "Value"), None)
            if d and v and d.strip().lower() == "scanname":
                out["scan"] = v.strip()
    n_shapes = {"Polygon": 0, "Ellipse": 0, "Rectangle": 0, "mask only": 0}
    for roi in root.iter():
        if _tag(roi) != "ROI":
            continue
        label, poly, masks = None, None, {}
        for sh in roi.iter():
            t = _tag(sh)
            if t == "Label" and label is None:
                label = (sh.get("Text") or "").strip()
            elif t == "Polygon" and sh.get("Points"):
                try:
                    poly = np.array([[float(v) for v in q.split(",")[:2]] for q in sh.get("Points").split()], np.float64)
                    n_shapes["Polygon"] += 1
                except ValueError:
                    poly = None
            elif t == "Ellipse":
                poly = _ellipse(float(sh.get("X")), float(sh.get("Y")), float(sh.get("RadiusX")), float(sh.get("RadiusY")))
                n_shapes["Ellipse"] += 1
            elif t == "Rectangle":
                x, y, w, h = (float(sh.get(k)) for k in ("X", "Y", "Width", "Height"))
                poly = np.array([[x, y], [x + w, y], [x + w, y + h], [x, y + h]], np.float64)
                n_shapes["Rectangle"] += 1
            elif t == "Mask":
                b = next((c for c in sh if _tag(c) == "BinData"), None)
                m = decode_mask(sh, b) if b is not None else None
                if m:
                    name = (sh.get("Text") or "").strip()
                    masks[name.lower()] = {"name": name, **m}
        if not label:
            continue
        if poly is None or len(poly) < 3:
            if not masks:
                continue
            x0 = min(m["bbox"][0] for m in masks.values()); y0 = min(m["bbox"][1] for m in masks.values())
            x1 = max(m["bbox"][2] for m in masks.values()); y1 = max(m["bbox"][3] for m in masks.values())
            poly = np.array([[x0, y0], [x1, y0], [x1, y1], [x0, y1]], np.float64)
            n_shapes["mask only"] += 1
        out["rois"][_roi_key(label)] = {"label": label, "poly": poly, "masks": masks}
    log(f"[geomx] scan {Path(tiff).name}: {len(out['rois'])} ROIs in the OME-XML ("
        + ", ".join(f"{v} {k}" for k, v in n_shapes.items() if v) + ")"
        + (f", {out['px'][0]:.4f} µm/px" if out["px"] else ", no pixel size"))
    return out


# ---------------------------------------------------------------------------------------
# images: scans side by side, or a drawn grid
# ---------------------------------------------------------------------------------------
class MosaicSource:
    """Several sources left to right on one canvas (pyramid.Source interface), `gap` px apart,
    tops aligned, the rest zero. `offsets[i]` is the x of source i's left edge."""

    def __init__(self, sources: list, gap: int = MOSAIC_GAP):
        self._src = sources
        first = sources[0]
        self.dtype = np.dtype(first.dtype)
        self.names = list(first.names)
        for s in sources[1:]:
            if np.dtype(s.dtype) != self.dtype or len(s.names) != len(self.names):
                raise SystemExit("[geomx] the scans differ in channel count or bit depth: they cannot share one "
                                 "canvas. Prepare one dataset per scan.")
        self.offsets, x = [], 0
        for s in sources:
            self.offsets.append(x)
            x += s.width + gap
        self.width = x - gap
        self.height = max(s.height for s in sources)
        self.is_rgb = False
        self.pixel_size_um = getattr(first, "pixel_size_um", None)

    @property
    def nchannels(self) -> int:
        return len(self.names)

    def read_band(self, c: int, y0: int, y1: int) -> np.ndarray:
        out = np.zeros((y1 - y0, self.width), self.dtype)
        for s, x in zip(self._src, self.offsets):
            if y0 >= s.height:
                continue
            b = s.read_band(c, y0, min(y1, s.height))
            out[: b.shape[0], x:x + s.width] = b
        return out

    def close(self):
        for s in self._src:
            s.close()


class ArraySource:
    """An in-memory (H, W) uint8 picture as a one-channel source."""

    def __init__(self, a: np.ndarray, name: str):
        self._a = a
        self.names = [name]
        self.height, self.width = a.shape
        self.dtype = np.dtype(np.uint8)
        self.is_rgb = False
        self.pixel_size_um = None

    @property
    def nchannels(self) -> int:
        return 1

    def read_band(self, c, y0, y1):
        return self._a[y0:y1]

    def close(self):
        pass


def grid_layout(slide_of: list[int], labels: list[tuple[str, str]], n_slides: int):
    """A drawn grid for AOIs without positions: one block per slide, one cell per AOI.

    Returns (picture uint8 (H, W), centres (n, 2), circles [(cx, cy, r)]).
    `labels[i]` = (ROI label, segment), already de-identified; slide headers are "Slide N".
    """
    from PIL import Image, ImageDraw, ImageFont
    n = len(slide_of)
    cols = min(GRID_COLS, max(4, int(math.ceil(math.sqrt(max(1, n))))))
    blocks, y = [], 0
    for s in range(n_slides):
        rows = [i for i in range(n) if slide_of[i] == s]
        nr = max(1, int(math.ceil(len(rows) / cols)))
        blocks.append((s, rows, y))
        y += GRID_HEAD + nr * GRID_CELL
    W, H = cols * GRID_CELL, max(y, GRID_CELL)
    im = Image.new("L", (W, H), 0)
    dr = ImageDraw.Draw(im)

    def font(size):
        try:
            return ImageFont.load_default(size=size)
        except TypeError:                      # Pillow < 10.1: one small bitmap font
            return ImageFont.load_default()
    f_head, f_cell = font(44), font(26)
    centres = np.zeros((n, 2), np.float64)
    circles = [None] * n
    r = 0.36 * GRID_CELL
    for s, rows, y0 in blocks:
        dr.text((16, y0 + 22), f"Slide {s + 1}" + (f"  ({len(rows)} AOIs)" if rows else ""), fill=235, font=f_head)
        dr.line([(0, y0 + GRID_HEAD - 2), (W, y0 + GRID_HEAD - 2)], fill=90, width=2)
        for k, i in enumerate(rows):
            cx = (k % cols + 0.5) * GRID_CELL
            cy = y0 + GRID_HEAD + (k // cols + 0.5) * GRID_CELL
            centres[i] = (cx, cy)
            circles[i] = (cx, cy, r)
            dr.rectangle([cx - GRID_CELL / 2 + 3, cy - GRID_CELL / 2 + 3, cx + GRID_CELL / 2 - 3, cy + GRID_CELL / 2 - 3],
                         outline=45, width=1)
            roi, seg = labels[i]
            dr.text((cx - GRID_CELL / 2 + 10, cy - GRID_CELL / 2 + 8), str(roi)[:14], fill=220, font=f_cell)
            if seg:
                dr.text((cx - GRID_CELL / 2 + 10, cy + GRID_CELL / 2 - 38), str(seg)[:16], fill=150, font=f_cell)
    return np.asarray(im, np.uint8), centres, circles


# ---------------------------------------------------------------------------------------
# prepare: everything into <cache>/<id>/
# ---------------------------------------------------------------------------------------
def _float(v: str) -> float:
    try:
        return float(str(v).replace(",", ""))
    except ValueError:
        return float("nan")


def _is_ntc(t: dict, r: list[str]) -> bool:
    ro = t["roles"]
    if "slide" in ro and _norm(r[ro["slide"]]) == NTC_NAME:
        return True
    return all(not r[ro[k]] for k in ("roi", "segment", "area") if k in ro) and any(k in ro for k in ("roi", "segment"))


def _scan_args(items) -> list[tuple[str | None, Path]]:
    out = []
    for it in items or []:
        name, path = None, it
        if "=" in it and not Path(it).exists():
            name, path = it.split("=", 1)
        p = Path(path)
        if not p.exists():
            raise SystemExit(f"[geomx] --image: {p} not found")
        out.append((name.strip() if name else None, p))
    return out


def prepare(args, log=print, geojson=None) -> Path:
    """`prepare.py geomx`: build <out>/<id>/ (image.zarr, expr/rois.zarr, overlays/rois.geojson,
    geomx/rows.json, manifest.json). `geojson`: prepare.parse_geojson_args, for extra --geojson overlays."""
    import pyramid
    from expr import _write_layer
    from scipy import sparse

    t0 = time.time()
    ds_dir = Path(args.out) / args.id
    ds_dir.mkdir(parents=True, exist_ok=True)

    # ---- annotation, DCC, PKC -------------------------------------------------------------
    ann = read_annotation(args.annotation, args.sheet, log)
    ro, H = ann["roles"], ann["header"]
    sid = ro["sample_id"]
    ntc_rows = [r for r in ann["rows"] if _is_ntc(ann, r)]
    rows = [r for r in ann["rows"] if not _is_ntc(ann, r)]
    col = lambda role, r: (r[ro[role]] if role in ro else "")
    scans = _scan_args(args.image)
    use_grid = args.layout == "grid" or not scans
    if scans and "roi" not in ro:
        raise SystemExit(f"[geomx] the annotation has no 'roi' column to find each AOI's outline in the scan by "
                         f"(columns: {H})")
    scan_info = []
    if scans:
        for given, path in scans:
            info = read_scan(path, log)
            info["path"] = path
            info["name"] = given or info["scan"] or re.sub(r"\.ome$", "", path.stem, flags=re.I)
            scan_info.append(info)
        names = {s["name"].strip().lower(): i for i, s in enumerate(scan_info)}
        if len(names) != len(scan_info):
            raise SystemExit(f"[geomx] two --image files are the same scan ({[s['name'] for s in scan_info]})")
        if "scan" not in ro:
            if len(scan_info) > 1:
                raise SystemExit("[geomx] several --image scans need a 'scan name' column in the annotation")
            keep = rows
        else:
            keep = [r for r in rows if r[ro["scan"]].strip().lower() in names]
            if not keep:
                have = sorted({r[ro["scan"]] for r in rows})
                raise SystemExit(f"[geomx] no annotation row has the scan name of the image(s) "
                                 f"{[s['name'] for s in scan_info]}. Scan names in the annotation: {have[:12]}. "
                                 f"Use --image 'ScanName=path' to say which scan a file is.")
        if len(keep) < len(rows):
            log(f"[geomx] {len(rows) - len(keep)} annotation rows belong to scans without an --image: left out")
        rows = keep
        if "scan" in ro:                 # rows in the scans' left-to-right order, so "Slide 1" is the leftmost scan
            rows = sorted(rows, key=lambda r: names[r[ro["scan"]].strip().lower()])
    if not rows:
        raise SystemExit("[geomx] the annotation has no AOI rows")

    store = DccStore(Path(args.dcc))
    missing = [r[sid] for r in rows if r[sid] not in store.files]
    if missing:
        share = len(missing) / len(rows)
        if share > MAX_MISSING_DCC or len(missing) == len(rows):
            raise SystemExit(f"[geomx] {len(missing)} of {len(rows)} annotation rows have no <Sample_ID>.dcc under "
                             f"{args.dcc} (e.g. {missing[:3]}; files there look like {sorted(store.files)[:2]})")
        log(f"[geomx] WARNING {len(missing)} annotation rows have no DCC file and are left out: {missing[:5]}")
        rows = [r for r in rows if r[sid] in store.files]
    n = len(rows)
    log(f"[geomx] {n} AOIs to read; {len(store.files)} DCC files under --dcc, {len(store.files) - n} not used "
        f"(other scans, controls)")

    pk = read_pkcs(args.pkc, log)
    n_probes, n_targets = len(pk["probe_target"]), len(pk["targets"])
    P = np.zeros((n_probes, n), np.float64)
    seq = {k: np.full(n, np.nan) for k in ("raw", "trimmed", "stitched", "aligned", "umiQ30", "rtsQ30")}
    dedup = np.zeros(n)
    unmapped = np.zeros(n)
    plates = []
    for j, r in enumerate(rows):
        d = parse_dcc(store.text(r[sid]))
        for k in seq:
            seq[k][j] = d.get(k, np.nan)
        plates.append(d.get("plate") or "")
        for code, cnt in d["codes"].items():
            i = pk["rts"].get(code)
            if i is None:
                unmapped[j] += cnt
            else:
                P[i, j] += cnt
        dedup[j] = sum(d["codes"].values())
    tot_all = float(dedup.sum())
    mapped = 1.0 - float(unmapped.sum()) / tot_all if tot_all else 0.0
    log(f"[geomx] DCC: {n} files read, {tot_all:,.0f} deduplicated counts, {mapped:.2%} on probes of the PKC(s)  "
        f"({time.time() - t0:.0f}s)")
    if mapped < MIN_MAPPED:
        raise SystemExit(f"[geomx] only {mapped:.1%} of the DCC counts are on probes of the PKC file(s) given "
                         f"(at least {MIN_MAPPED:.0%} needed): wrong or missing --pkc?")
    ntc = []
    for r in ntc_rows:                                  # the No Template Control of each plate the AOIs came from
        if r[sid] in store.files:
            d = parse_dcc(store.text(r[sid]))
            if (d.get("plate") or "") in set(plates):
                ntc.append({"plate": d.get("plate") or "", "total": float(sum(d["codes"].values()))})
    store.close()

    T = collapse(P, pk["probe_target"], n_targets)
    neg_t = [i for i, t in enumerate(pk["targets"]) if t["negative"]]
    gene_t = [i for i, t in enumerate(pk["targets"]) if not t["negative"]]
    genes = [pk["targets"][i]["name"] for i in gene_t]
    G = T[gene_t]                                           # genes x AOIs
    # background and LOQ per module that has negative probes
    neg = []
    for mi, mod in enumerate(pk["modules"]):
        pr = np.flatnonzero(np.array([pk["targets"][t]["negative"] and pk["targets"][t]["module"] == mi
                                      for t in pk["probe_target"]]))
        if not len(pr):
            continue
        gm, gsd = _geo(P[pr], axis=0)
        neg.append({"module": mi, "geoMean": gm, "geoSD": gsd, "loq": np.maximum(LOQ_MIN, gm * gsd ** LOQ_SD)})
    neg.sort(key=lambda b: -pk["modules"][b["module"]]["targets"])     # the main panel first: its values are the QC metrics
    mod_of_gene = np.array([pk["targets"][i]["module"] for i in gene_t])
    above = np.zeros(n)
    detected = np.zeros(len(gene_t))                     # per target: in how many AOIs it is above its module's LOQ
    if neg:
        loq_by_mod = {b["module"]: b["loq"] for b in neg}
        first = neg[0]["loq"]
        for mi in np.unique(mod_of_gene):
            sel = mod_of_gene == mi
            hit = G[sel] > loq_by_mod.get(int(mi), first)[None, :]
            above += hit.sum(0)
            detected[sel] = hit.sum(1)
    q3 = np.percentile(G, 75, axis=0) if len(G) else np.zeros(n)
    q3norm = q3_values(G, T[neg_t], detected if neg else None, n, log)
    area = np.array([_float(col("area", r)) for r in rows])
    nuclei = np.array([_float(col("nuclei", r)) for r in rows])
    log(f"[geomx] {len(genes):,} targets ({len(neg_t)} negative-probe target(s) kept out of the genes), "
        f"median {np.median(G.sum(0)):,.0f} counts per AOI"
        + (f", background {np.median(neg[0]['geoMean']):.2f}, LOQ {np.median(neg[0]['loq']):.2f} (medians)" if neg else ""))

    # ---- slides, scans, labels --------------------------------------------------------------
    def ordinal(values):
        first = {}
        for v in values:
            first.setdefault(v, len(first))
        return [first[v] for v in values], len(first)
    slide_of, n_slides = ordinal([col("slide", r) or col("scan", r) for r in rows])
    roi_lab = [deid.anonymize(col("roi", r), "text", "geomx roi label") for r in rows]
    seg_lab = [deid.anonymize(col("segment", r), "text", "geomx segment") for r in rows]

    # ---- positions, ROI polygons, image -----------------------------------------------------
    xy = np.full((n, 2), np.nan, np.float64)
    feats, px_um, scan_meta = [], None, []
    area_px = np.full(n, np.nan)
    if not use_grid:
        sources = []
        for s in scan_info:
            src = pyramid.TiffStackSource(s["path"], names=s["channels"] or None, rgb_split=True)
            if s["channels"] and len(s["channels"]) != len(src.names):
                src.close()
                src = pyramid.TiffStackSource(s["path"], rgb_split=True)
            sources.append(src)
            s["um"] = (s["px"][0] if s["px"] else None) or getattr(src, "pixel_size_um", None)
        ums = [s["um"] for s in scan_info if s["um"]]
        if ums and max(ums) / min(ums) > 1.01:
            raise SystemExit(f"[geomx] the scans have different pixel sizes ({ums}): one dataset per scan")
        px_um = args.pixel_size or (ums[0] if ums else None)
        src = MosaicSource(sources) if len(sources) > 1 else sources[0]
        offsets = src.offsets if len(sources) > 1 else [0]
        idx_of = {s["name"].strip().lower(): i for i, s in enumerate(scan_info)}
        by_roi: dict[tuple, list[int]] = {}
        unplaced = []
        for j, r in enumerate(rows):
            si = idx_of[r[ro["scan"]].strip().lower()] if "scan" in ro else 0
            g = scan_info[si]["rois"].get(_roi_key(col("roi", r)))
            if g is None:
                unplaced.append(j)
                continue
            by_roi.setdefault((si, _roi_key(col("roi", r))), []).append(j)
        for (si, rk), members in by_roi.items():
            g = scan_info[si]["rois"][rk]
            poly = g["poly"] + np.array([offsets[si], 0.0])
            a_px, cx, cy = poly_area_centroid(poly)
            pts = []
            for j in members:
                m = g["masks"].get(col("segment", rows[j]).strip().lower())
                if m is None and len(g["masks"]) == 1 and len(members) == 1:
                    m = next(iter(g["masks"].values()))
                if m:
                    area_px[j] = m["areaPx"]
                    pts.append(interior_point(poly, m["cx"] + offsets[si], m["cy"]))
                else:
                    area_px[j] = a_px / len(members)
                    pts.append(interior_point(poly, cx, cy))
            # AOIs of one ROI must not share a place (the viewer finds a row by its coordinates)
            if len(members) > 1 and min(math.dist(pts[a], pts[b]) for a in range(len(pts)) for b in range(a)) < 2.0:
                ix, iy = interior_point(poly, cx, cy)
                a, b = chord(poly, ix, iy)
                pts = [(a + (k + 1) * (b - a) / (len(members) + 1), iy) for k in range(len(members))]
            for j, p in zip(members, pts):
                xy[j] = p
            lab = deid.anonymize(g["label"], "text", "geomx roi label")
            name = f"ROI {lab}" if len(scan_info) == 1 else f"Slide {slide_of[members[0]] + 1} · ROI {lab}"
            ring = [[round(float(x), 2), round(float(y), 2)] for x, y in poly]
            feats.append({"type": "Feature", "properties": {"name": name, "roi": lab, "slide": slide_of[members[0]] + 1,
                                                            "aois": len(members)},
                          "geometry": {"type": "Polygon", "coordinates": [ring + [ring[0]]]}})
        for si, s in enumerate(scan_info):
            mine = [j for j, r in enumerate(rows) if ("scan" not in ro or idx_of[r[ro["scan"]].strip().lower()] == si)]
            miss = [j for j in mine if j in set(unplaced)]
            if mine and len(miss) / len(mine) > MAX_UNPLACED:
                raise SystemExit(f"[geomx] {len(miss)} of {len(mine)} AOIs of scan '{s['name']}' have no ROI with their "
                                 f"label in {s['path'].name} (labels there: {sorted(x['label'] for x in s['rois'].values())[:8]}): "
                                 f"wrong scan file, or an export without ROIs? Use --layout grid to open it without positions.")
            # the masks' areas against the annotation's: says at once if this is another scan's file
            # (ROI numbers repeat from scan to scan) or the pixel size is off
            um2 = (s["px"][0] * s["px"][1]) if s["px"] else (px_um * px_um if px_um else None)
            rt = np.array([area_px[j] * um2 / area[j] for j in mine
                           if um2 and np.isfinite(area_px[j]) and np.isfinite(area[j]) and area[j] > 0])
            ratio = float(np.median(rt)) if len(rt) else None
            if ratio is not None and not 0.9 <= ratio <= 1.1:
                log(f"[geomx] WARNING scan '{s['name']}': the ROI areas in {s['path'].name} are {ratio:.2f} x the areas in "
                    f"the annotation (median over {len(rt)} AOIs). Another scan's file, or a wrong pixel size? "
                    f"The ROIs may be drawn in the wrong place.")
            scan_meta.append({"scan": si + 1, "offsetPx": [int(offsets[si]), 0],
                              "size": [int(sources[si].width), int(sources[si].height)],
                              "rois": len(s["rois"]), "aois": len(mine), "unplaced": len(miss),
                              "maskAreaOverAnnotation": round(ratio, 4) if ratio is not None else None})
        if unplaced:
            log(f"[geomx] WARNING {len(unplaced)} AOIs have no ROI in the scan's OME-XML: no position, not drawn")
        image = pyramid.build_image(src, ds_dir, tile=args.tile, log=log)
        src.close()
        if px_um and not np.isfinite(area).all():            # no area column: from the mask / polygon
            um2 = px_um * (scan_info[0]["px"][1] if scan_info[0]["px"] else px_um)
            area = np.where(np.isfinite(area), area, area_px * um2)
        med = float(np.nanmedian(area)) if np.isfinite(area).any() else None
        spot_px = DISC * ((2.0 * math.sqrt(med / math.pi) / px_um) if (med and px_um) else
                          (2.0 * math.sqrt(float(np.nanmedian(area_px)) / math.pi) if np.isfinite(area_px).any() else 200.0))
    else:
        pic, centres, circles = grid_layout(slide_of, list(zip(roi_lab, seg_lab)), n_slides)
        xy = centres
        src = ArraySource(pic, "ROI grid (no scan image)")
        image = pyramid.build_image(src, ds_dir, tile=args.tile, colors=["#ffffff"], log=log)
        spot_px = 0.5 * GRID_CELL
        for j, (cx, cy, rad) in enumerate(circles):
            ring = [[round(float(x), 2), round(float(y), 2)] for x, y in _ellipse(cx, cy, rad, rad, 48)]
            nm = f"ROI {roi_lab[j]}" + (f" · {seg_lab[j]}" if seg_lab[j] else "") + (f" · Slide {slide_of[j] + 1}" if n_slides > 1 else "")
            feats.append({"type": "Feature", "properties": {"name": nm, "roi": roi_lab[j], "slide": slide_of[j] + 1, "aois": 1},
                          "geometry": {"type": "Polygon", "coordinates": [ring + [ring[0]]]}})
        log(f"[geomx] no scan image: {n} AOIs on a drawn grid ({pic.shape[1]} x {pic.shape[0]} px, {n_slides} slide block(s)); "
            f"positions are a layout, not the tissue")

    placed = np.isfinite(xy[:, 0])
    k = np.round(xy[placed] * 100).astype(np.int64)
    if len(np.unique(k, axis=0)) != int(placed.sum()):
        raise SystemExit("[geomx] two AOIs ended up at the same position (internal error: please report the files)")

    # ---- the layer --------------------------------------------------------------------------
    gm = sparse.csr_matrix(G.astype(np.float32))
    gm.eliminate_zeros()
    gm.sort_indices()
    meta = {"key": LAYER_KEY, "label": "GeoMx ROIs", "kind": "rois", "platform": "geomx", "nBins": int(n),
            "nGenes": int(len(genes)), "binSizeUm": None, "spotDiameterPx": round(float(spot_px), 2),
            "micronsPerPixel": float(px_um) if px_um else None, "imageScale": 1.0,
            "unit": "counts", "nRois": len(feats) if not use_grid else len({(s, r) for s, r in zip(slide_of, roi_lab)}),
            "layout": "grid" if use_grid else "scan"}
    _write_layer(ds_dir, xy.astype(np.float32), gm, genes, meta, [r[sid] for r in rows])
    log(f"[expr] {LAYER_KEY}: {n} AOIs x {len(genes):,} targets, nnz {gm.nnz:,}, dot {spot_px:.0f} px")

    # ---- overlay ----------------------------------------------------------------------------
    od = ds_dir / "overlays"
    od.mkdir(exist_ok=True)
    gj = {"type": "FeatureCollection", "features": feats}
    deid.assert_clean(gj, OVERLAY_FILE)
    (od / OVERLAY_FILE).write_text(json.dumps(gj, separators=(",", ":")), encoding="utf-8")
    overlays = [{"key": "rois", "label": "ROI outlines", "file": OVERLAY_FILE}]
    if geojson and getattr(args, "geojson", None):
        overlays += geojson(args.geojson, ds_dir)

    # ---- per-row values the export needs (cache only: raw annotation text lives here) ---------
    def role_of_name(h):                # by spelling, so a second scan / slide column is treated as one too
        return next((k for k, names in ROLES.items() if _norm(h) in names), None)
    columns = [{"name": H[c], "role": role_of_name(H[c]), "values": [r[c] for r in rows]}
               for c in range(len(H)) if H[c] and role_of_name(H[c]) not in ("sample_id", "aoi", "area", "nuclei", "roi")]
    qc = {"rawReads": seq["raw"], "pctAligned": 100.0 * seq["aligned"] / np.maximum(seq["raw"], 1.0),
          "saturation": 100.0 * (1.0 - dedup / np.maximum(seq["aligned"], 1.0)),
          "area": area, "nuclei": nuclei, "q3": q3}
    if neg:
        qc.update({"negGeoMean": neg[0]["geoMean"], "loq": neg[0]["loq"], "genesAboveLoq": above,
                   "q3OverNeg": q3 / np.maximum(neg[0]["geoMean"], 1e-9)})
    cd = ds_dir / CACHE_DIR
    cd.mkdir(exist_ok=True)
    (cd / "rows.json").write_text(json.dumps({
        "n": n, "keepSlideNames": bool(args.keep_slide_names), "slide": slide_of, "columns": columns,
        "qc": {k: [None if not np.isfinite(x) else float(x) for x in v] for k, v in qc.items()},
        "q3norm": q3norm}), encoding="utf-8")

    # ---- manifest ---------------------------------------------------------------------------
    plate_ord, _ = ordinal(plates)
    ntc_plate = {p: i for p, i in zip(plates, plate_ord)}
    info = {
        "nAoi": int(n), "nRoi": meta["nRois"], "nSlides": int(n_slides), "layout": meta["layout"],
        "panels": [{"name": deid.anonymize(m["name"], "text", "geomx panel"), "version": m["version"],
                    "targets": m["targets"], "probes": m["probes"], "negativeProbes": m["negativeProbes"]}
                   for m in pk["modules"]],
        "collapse": COLLAPSE_RULE, "loq": LOQ_RULE if neg else None,
        "normalisation": "stored: target counts as collapsed, not normalised; plots and volcano: log1p of the Q3-normalised "
                         "counts (GeomxTools Q3, factors in the layer's 'q3' entry at export; the default), or log1p counts per "
                         "10,000; the colours on the slide: Q3-normalised counts, or the raw counts",
        "dcc": {"used": int(n), "inFolder": len(store.files), "sharePkcMapped": round(mapped, 5)},
        "scans": scan_meta,
        "ntc": [{"plate": int(ntc_plate[x["plate"]]) + 1, "total": x["total"]} for x in ntc],
        "background": [{"panel": b["module"] + 1, "medianNegGeoMean": round(float(np.median(b["geoMean"])), 4),
                        "medianLoq": round(float(np.median(b["loq"])), 4)} for b in neg],
        "unplaced": int((~placed).sum()),
    }
    manifest = {"id": args.id, "label": args.label or args.id, "kind": "geomx", "pixelSizeUm": float(px_um) if px_um else None,
                "created": datetime.now().isoformat(timespec="seconds"),
                "source": {"dcc": str(args.dcc), "pkc": [str(p) for p in args.pkc],
                           "annotation": [str(p) for p in args.annotation], "images": [str(p) for _, p in scans]},
                "expr": [meta], "overlays": overlays, "image": image, "geomx": info}
    deid.assert_clean({k: v for k, v in manifest.items() if k != "source"}, "manifest.json")
    (ds_dir / "manifest.json").write_text(json.dumps(manifest, indent=1), encoding="utf-8")
    log(f"[done] {ds_dir}\\manifest.json  ({time.time() - t0:.0f}s)")
    return ds_dir


# ---------------------------------------------------------------------------------------
# export: QC metrics, all-target chunks for the volcano, the annotation columns as groupings
# ---------------------------------------------------------------------------------------
QC_METRICS = [   # key, label, file, unit, scale
    ("rawReads", "Raw reads", "raw_reads.bin", "reads", "log"),
    ("pctAligned", "% of raw reads aligned", "pct_aligned.bin", "%", "linear"),
    ("saturation", "Sequencing saturation %", "saturation.bin", "%", "linear"),
    ("area", "AOI area (µm²)", "area.bin", "µm²", "log"),
    ("nuclei", "Nuclei", "nuclei.bin", "nuclei", "log"),
    ("negGeoMean", "Background (negative probes, geometric mean)", "neg_geomean.bin", "counts", "linear"),
    ("loq", "Limit of quantification (LOQ)", "loq.bin", "counts", "linear"),
    ("genesAboveLoq", "Targets above LOQ", "genes_above_loq.bin", "targets", "log"),
    ("q3", "Q3 of the target counts", "q3.bin", "counts", "log"),
    ("q3OverNeg", "Signal to background (Q3 / negative probes)", "q3_over_neg.bin", "ratio", "linear"),
]
_SLIDE_ROLES = {"slide": "Slide", "scan": "Scan"}
# columns of a DSP Data Analysis SegmentProperties sheet that describe the scan or the sequencing run, not the tissue
_TECHNICAL = {"scanwidth", "scanheight", "scanoffsetx", "scanoffsety", "roicoordinatex", "roicoordinatey", "roiid",
              "segmentid", "segmentdisplayname", "umiq30", "rtsq30", "rawreads", "alignedreads", "deduplicatedreads",
              "trimmedreads", "stitchedreads", "sequencingsaturation", "sequencingsetid", "geomxngspipelineversion"}


def _cache_rows(ds_dir: Path, er) -> tuple[dict, np.ndarray]:
    p = Path(ds_dir) / CACHE_DIR / "rows.json"
    if not p.exists():
        raise SystemExit(f"[geomx] {p} is missing: run `prepare.py geomx` again for this dataset")
    rows = json.loads(p.read_text(encoding="utf-8"))
    if rows["n"] != len(er.xy):
        raise SystemExit(f"[geomx] {p} has {rows['n']} rows, the layer {len(er.xy)}: run `prepare.py geomx` again")
    return rows, ~np.isnan(er.xy[:, 0])


def write_qc(er, d: Path, rows: dict, on: np.ndarray, log=print) -> dict | None:
    """binqc/ as binqc.py writes it (total counts, targets detected, % MT) plus the GeoMx values per AOI."""
    import binqc
    block = binqc.write_bin_qc(er, d, log=log)
    if not block:
        return None
    hd = Path(d) / block["dir"]
    for key, label, fname, unit, scale in QC_METRICS:
        v = rows["qc"].get(key)
        if v is None:
            continue
        a = np.array([np.nan if x is None else x for x in v], np.float64)[on]
        if not len(a) or not np.isfinite(a).all():
            continue                                   # a value missing for some AOI: the metric is left out
        (hd / fname).write_bytes(np.ascontiguousarray(a, dtype="<f4").tobytes())
        block["metrics"].append({"key": key, "label": label, "file": fname, "dtype": "float32", "unit": unit,
                                 "scale": scale, **binqc._stats(a.astype(np.float32))})
    for m in block["metrics"]:
        if m["key"] == "genes":
            m["label"] = "Targets detected"
    deid.assert_clean(block, "binqc/ manifest block")
    log(f"[geomx] binqc: {len(block['metrics'])} metrics per AOI ({', '.join(m['key'] for m in block['metrics'])})")
    return block


def write_target_chunks(er, d: Path, mpp: float | None, chunk_um: float = 500.0, log=print) -> dict | None:
    """hvg/ in hvgchunks.py's format, but holding EVERY target with a count (no variable-gene
    selection, no minimum number of rows): the volcano reads these chunks, and for a few dozen
    ROIs it should test every target. Counts are rounded and clipped at 65,535 (uint16)."""
    from scipy import sparse
    xy = er.xy
    on = ~np.isnan(xy[:, 0])
    row_of = np.cumsum(on) - 1
    Gn, B = len(er.genes), len(xy)
    gm = sparse.csr_matrix((np.asarray(er.g["data"], dtype=np.float32), np.asarray(er.g["indices"]), er.indptr),
                           shape=(Gn, B))
    X = gm.T.tocsr()
    X.sort_indices()
    lib = np.asarray(X.sum(1)).ravel().astype(np.float32)
    cols = np.flatnonzero(np.asarray(er.totals) > 0)
    if len(cols) > 65535:
        cols = np.sort(cols[np.argsort(-np.asarray(er.totals)[cols], kind="stable")[:65535]])
    if not len(cols) or not on.any():
        return None
    genes = [er.genes[i] for i in cols]
    deid.assert_clean(genes, "hvg/genes.json")
    Xs = X[:, cols].tocsr()
    n_clip = int((Xs.data > 65535).sum())
    n_frac = int((np.abs(Xs.data - np.rint(Xs.data)) > 1e-6).sum())
    Xs.data = np.minimum(np.rint(Xs.data), 65535).astype(np.float32)
    Xs.eliminate_zeros()
    Xs.sort_indices()
    chunk_px = int(round(chunk_um / mpp)) if mpp else 1024
    keep = np.flatnonzero(on)
    xy64 = xy[keep].astype(np.float64)
    tx = np.floor(xy64[:, 0] / chunk_px).astype(np.int64)
    ty = np.floor(xy64[:, 1] / chunk_px).astype(np.int64)
    hd = Path(d) / "hvg"
    hd.mkdir(parents=True, exist_ok=True)
    for old in hd.glob("*.bin"):
        old.unlink()
    order = np.lexsort((keep, ty, tx))
    keys = np.stack([tx[order], ty[order]], 1)
    split = np.flatnonzero(np.any(np.diff(keys, axis=0) != 0, axis=1)) + 1
    chunks, total = [], 0
    for grp in np.split(order, split):
        bins = keep[grp]
        sub = Xs[bins]
        buf = b"".join([np.uint32(len(bins)).astype("<u4").tobytes(), row_of[bins].astype("<u4").tobytes(),
                        lib[bins].astype("<f4").tobytes(), sub.indptr.astype("<u4").tobytes(),
                        sub.indices.astype("<u2").tobytes(), sub.data.astype("<u2").tobytes()])
        cx, cy = int(tx[grp[0]]), int(ty[grp[0]])
        (hd / f"{cx}_{cy}.bin").write_bytes(buf)
        chunks.append({"tx": cx, "ty": cy, "n": int(len(bins)), "bytes": len(buf)})
        total += len(buf)
    index = {"chunkPx": chunk_px, "chunkUm": chunk_um if mpp else None, "nGenes": len(genes), "minUmi": 10,
             "selection": "all targets with a count (ROI layer: no variable-gene selection)",
             "rounded": n_frac, "clipped": n_clip, "chunks": chunks}
    (hd / "genes.json").write_text(json.dumps(genes, separators=(",", ":")), encoding="utf-8")
    (hd / "index.json").write_text(json.dumps(index, separators=(",", ":")), encoding="utf-8")
    log(f"[geomx] hvg/: all {len(genes):,} targets, {len(chunks)} chunks, {total / 1e6:.1f} MB"
        + (f"; {n_frac:,} values rounded to whole counts (multi-probe geometric means)" if n_frac else "")
        + (f"; WARNING {n_clip:,} values above 65,535 clipped (volcano only; gene files are exact)" if n_clip else ""))
    return {"nGenes": len(genes), "chunkPx": chunk_px, "chunks": len(chunks), "bytes": total, "allTargets": True}


def write_groupings(er, d: Path, rows: dict, on: np.ndarray, log=print) -> dict | None:
    """umap/ with the annotation columns as named clusterings ("source": "user") and NO embedding.

    Same files and index.json fields as userclusters.py writes, same de-identification
    (userclusters.deidentify + deid.assert_clean). Slide / scan names become "Slide 1" ...
    unless prepare was run with --keep-slide-names. Returns the manifest's "umap" entry.
    """
    import sectionumap
    import userclusters as uc
    key = er.meta.get("key")
    n_rows = int(on.sum())
    h5_of_row = np.flatnonzero(on)
    ud = Path(d) / "umap"
    ud.mkdir(parents=True, exist_ok=True)
    # shown in this order: the study's own columns, then segment / panel, slide and scan last. Looked at in the
    # other order, so that of two columns with the same groups the better-known one (slide, segment) is kept.
    rank = lambda c: 2 if c.get("role") in _SLIDE_ROLES else 1 if c.get("role") else 0
    cols = sorted([c for c in rows["columns"] if _norm(c["name"]) not in _TECHNICAL], key=lambda c: -rank(c))
    found, seen = [], []
    for c in cols:
        role, name = c.get("role"), c["name"]
        vals = [str(v) for v in c["values"]]
        if role in _SLIDE_ROLES and not rows.get("keepSlideNames"):
            first = {}
            for v in vals:
                if v.lower() not in uc._NA:
                    first.setdefault(v, len(first) + 1)
            vals = [f"{_SLIDE_ROLES[role]} {first[v]}" if v in first else "" for v in vals]
            name = _SLIDE_ROLES[role]
        cats, codes, why = uc.categorical(vals, None)
        safe_col = deid.anonymize(name, "text", "geomx annotation column")
        if cats is None:
            log(f"[geomx] annotation column '{safe_col}': not a grouping ({why})")
            continue
        if len(cats) > max(2, n_rows // 2):
            log(f"[geomx] annotation column '{safe_col}': not a grouping ({len(cats)} values for {n_rows} AOIs: an id)")
            continue
        safe, cats, codes, _ = uc.deidentify(name, cats, codes, None, key, log)
        lab = np.where(codes[h5_of_row] >= 0, codes[h5_of_row], -1).astype(np.int64)
        if len(set(lab[lab >= 0].tolist())) < 2:
            continue
        part = tuple(lab.tolist())
        canon = tuple({v: i for i, v in enumerate(dict.fromkeys(part))}[v] for v in part)
        if canon in seen:
            log(f"[geomx] annotation column '{safe}': the same groups as an earlier column, left out")
            continue
        seen.append(canon)
        found.append((rank(c), len(found), safe, cats, lab))
    out, labels = [], []
    for _, _, safe, cats, lab in sorted(found, key=lambda t: (t[0], t[1])):
        i = len(out) + 1
        out.append({"id": f"user{i}", "name": safe, "source": "user", "file": f"clusters_user{i}.bin",
                    "n": len(cats), "categories": cats})
        labels.append(lab)
    if not out:
        log(f"[geomx] {key}: no annotation column is a usable grouping: no legend, no group pickers")
        return None
    assigned = np.flatnonzero(np.any(np.stack(labels) >= 0, axis=0))
    Xn, lib = uc._log_cp10k(er, h5_of_row[assigned])
    has = lib > 0
    Xn = Xn[np.flatnonzero(has)]
    rows_m = assigned[has]
    for c, lab in zip(out, labels):
        la = lab[rows_m]
        m = la >= 0
        c["sizes"] = [int(v) for v in np.bincount(lab[lab >= 0], minlength=c["n"])]
        c["nAssigned"] = int((lab >= 0).sum())
        with np.errstate(divide="ignore", invalid="ignore"):
            c["markers"] = sectionumap._markers(Xn[np.flatnonzero(m)], la[m], c["n"], er.genes)
        log(f"[geomx] grouping '{c['name']}': {c['n']} groups, sizes {c['sizes'][:12]}{' ...' if c['n'] > 12 else ''}")
    index = {"nBins": n_rows, "nEmbedded": 0, "bounds": [0, 0, 1, 1],
             "embedding": {"source": "none", "why": NO_EMBED_WHY},
             "params": {"markers": f"log2 fold change of mean CP10k vs the other AOIs, pct >= {sectionumap.MIN_PCT:g}, lfc > 0"},
             "clusterings": out, "default": 0}
    deid.assert_clean(index, f"umap/index.json of {key}")
    files = {"coords.bin": np.full((n_rows, 2), np.nan, "<f4").tobytes()}
    for c, lab in zip(out, labels):
        files[c["file"]] = np.where(lab >= 0, lab, uc.NO_CLUSTER).astype("<u2").tobytes()
    for p in ud.glob("clusters_*.bin"):
        if p.name not in files:
            p.unlink()
    files["index.json"] = json.dumps(index, separators=(",", ":"), allow_nan=False).encode("utf-8")
    for name, buf in files.items():
        uc._write(ud, name, buf)
    entry = uc.manifest_entry(index)
    entry["embedding"] = "none"
    return entry


def export_layer(er, d: Path, ds_dir: Path, layer: dict, log=print) -> tuple[dict | None, dict | None, dict | None]:
    """export.py's hook for a layer of kind "rois": (binQc block, hvg entry, umap entry).
    No section UMAP, no Leiden, no suggested cell types, no variable-gene selection.
    Also sets layer["q3"] (q3_entry) on the manifest's layer dict itself: export.py writes the
    exported manifest's "expr" entries from these same dicts."""
    rows, on = _cache_rows(ds_dir, er)
    log(f"[geomx] {layer.get('key')}: ROI layer ({int(on.sum())} AOIs): no UMAP, Leiden clusters or cell-type "
        f"suggestions are computed; groupings come from the annotation")
    q3 = q3_entry(rows, on, log)
    layer.pop("q3", None)
    if q3:
        deid.assert_clean(q3, f"q3 entry of {layer.get('key')}")
        layer["q3"] = q3
    q = write_qc(er, Path(d), rows, on, log)
    h = write_target_chunks(er, Path(d), layer.get("micronsPerPixel"), log=log)
    u = write_groupings(er, Path(d), rows, on, log)
    return q, h, u
