"""Suggested cluster names from marker lists: score every embedded bin for known cell types.

A marker set (markers/nerve.json: {name, note, version, types: [{name, genes}]}) lists
genes per cell type. Each bin the whole-section UMAP embedded gets a score per type,
its best type, and each Leiden cluster a suggested name with a confidence. The browser
scores the circle tool's bins with the same rule (lasso_worker.js), so keep them identical.

    expr/<key>/umap/celltypes.json   {markerSet, types: [{name, genes, missing: [{gene, reason}], scored, reason?}], rule,
                                      clusterings: {"<res>": [{cluster, suggest, top, pct, second, pct2,
                                                               confidence, meanScore}]}}
    expr/<key>/umap/celltype.bin     uint8 per positions.bin row: index into types,
                                     254 = embedded but no type, 255 = not embedded

Rule, on the bins embedded in the UMAP (>= MIN_UMI UMI):
  a marker gene is missing with reason "absent" when the layer has no gene of that name
  (case-insensitive), "rare" when it is counted (> 0) in < MIN_DETECT of those bins - a single
  UMI of a gene almost never seen z-scores so high that it alone would decide the bin;
  x = log1p(count / libSize * 1e4) per remaining gene;
  z = (x - mean) / sd per gene over those bins (population sd);
  score[bin, type] = mean z over the type's genes; a type with < MIN_GENES genes left is not scored.
  Bin's type = argmax score (first in list order on ties) if >= BIN_MIN, else none.
  Cluster: pct[t] = share of the cluster's bins (none included) whose type is t; top/second by pct
  (stable, list order on ties). high: pct >= HIGH and >= 2x second; medium: >= MEDIUM; low: >= LOW;
  otherwise none -> "Unclear".

Reads umap/coords.bin and clusters_r*.bin as written by sectionumap.py: it never
recomputes the UMAP, and takes seconds.
"""
from __future__ import annotations

import json
import re
import time
from pathlib import Path

import numpy as np

import deid
from hvgchunks import _BARCODE, MIN_UMI

DEFAULT_MARKERS = Path(__file__).resolve().parent / "markers" / "rna_generic.json"
MIN_GENES = 2
MIN_DETECT = 0.005          # share of the scored bins a marker gene must be counted in
# 0.5 left 44% of the 8 um bins of s2 without a type and called clear MBP/MPZ clusters "maybe"
BIN_MIN = 0.25
HIGH, MEDIUM, LOW = 0.5, 0.3, 0.15
NO_TYPE = 254
NOT_EMBEDDED = 255
NO_CLUSTER = 0xFFFF


def rule_text() -> str:
    return (f"Each embedded bin (>= {MIN_UMI} UMI): log1p CP10k of each marker gene, z-scored over the "
            f"embedded bins (population sd), averaged per type; the bin's type is the best-scoring one if "
            f"its score >= {BIN_MIN:g}, else none. A marker gene counted in fewer than {MIN_DETECT:.1%} of "
            f"those bins is left out as rare. A type needs >= {MIN_GENES} genes left to be scored. "
            f"Per cluster, pct = share of its bins of each type (none counts in the denominator). "
            f"high: top pct >= {HIGH:g} and >= 2x the second; medium: >= {MEDIUM:g}; low: >= {LOW:g} "
            f"(shown as 'maybe'); otherwise Unclear.")


def load_markers(path: Path | None = None) -> tuple[dict, str]:
    """The parsed marker set and the file name the manifest/viewer refer to it by."""
    p = Path(path) if path else DEFAULT_MARKERS
    ms = json.loads(p.read_text(encoding="utf-8"))
    rel = f"markers/{p.name}"
    return ms, rel


def match_genes(types: list[dict], genes: list[str]) -> list[dict]:
    """Per type: column indices of its genes found (case-insensitive, first match), names found / missing."""
    first: dict[str, int] = {}
    for i, g in enumerate(genes):
        first.setdefault(g.upper(), i)
    out = []
    for t in types:
        cols, found, missing = [], [], []
        for g in t["genes"]:
            j = first.get(g.upper())
            if j is None:
                missing.append({"gene": g, "reason": "absent"})
            else:
                cols.append(j)
                found.append(genes[j])
        out.append({"name": t["name"], "cols": cols, "genes": found, "missing": missing})
    return out


def score_bins(cnt: dict[int, np.ndarray], lib: np.ndarray, matched: list[dict]) -> tuple[np.ndarray, np.ndarray]:
    """score (bins x types, NaN for types not scored) and scored flag per type.

    `cnt[col]` is the dense count vector of gene column `col` over the scored bins.
    A gene counted in < MIN_DETECT of the bins (sd 0 included) is dropped as rare;
    a type left with < MIN_GENES is not scored.
    """
    n = len(lib)
    z = {}
    for col, c in cnt.items():
        if np.count_nonzero(c) < MIN_DETECT * n:
            continue
        x = np.log1p(c / lib * 1e4)
        sd = x.std()
        if sd > 0:
            z[col] = (x - x.mean()) / sd
    score = np.full((n, len(matched)), np.nan)
    scored = np.zeros(len(matched), bool)
    for ti, t in enumerate(matched):
        cols = [c for c in t["cols"] if c in z]
        t["used"] = [g for g, c in zip(t["genes"], t["cols"]) if c in z]
        t["rare"] = [g for g, c in zip(t["genes"], t["cols"]) if c not in z]
        if len(cols) >= MIN_GENES:
            score[:, ti] = np.mean([z[c] for c in cols], axis=0)
            scored[ti] = True
    return score, scored


def best_type(score: np.ndarray, scored: np.ndarray, bin_min: float = BIN_MIN) -> np.ndarray:
    """Index of each bin's best type, NO_TYPE when its best score is < bin_min (or nothing is scored)."""
    best = np.full(len(score), NO_TYPE, np.int64)
    if not scored.any():
        return best
    s = np.where(scored[None, :], score, -np.inf)
    arg = np.argmax(s, axis=1)                         # first max: list order breaks ties
    ok = s[np.arange(len(s)), arg] >= bin_min
    best[ok] = arg[ok]
    return best


def summarise_clusters(lab: np.ndarray, best: np.ndarray, score: np.ndarray, scored: np.ndarray,
                       names: list[str]) -> list[dict]:
    """One suggestion per cluster id 0..max(lab)."""
    T = len(names)
    k = int(lab.max()) + 1
    out = []
    for ci in range(k):
        m = lab == ci
        nb = int(m.sum())
        bt = best[m]
        pct = np.bincount(bt[bt < T], minlength=T).astype(np.float64) / max(nb, 1)
        order = np.argsort(-pct, kind="stable")
        top, sec = int(order[0]), int(order[1]) if T > 1 else None
        p1 = float(pct[top])
        p2 = float(pct[sec]) if sec is not None else 0.0
        if p1 >= HIGH and p1 >= 2 * p2:
            conf = "high"
        elif p1 >= MEDIUM:
            conf = "medium"
        elif p1 >= LOW:
            conf = "low"
        else:
            conf = "none"
        suggest = names[top] if conf in ("high", "medium") else f"maybe {names[top]}" if conf == "low" else "Unclear"
        ms = [round(float(score[m, t].mean()), 3) if scored[t] and nb else None for t in range(T)]
        out.append({"cluster": ci, "suggest": suggest, "top": names[top], "pct": round(p1, 3),
                    "second": names[sec] if sec is not None else None, "pct2": round(p2, 3),
                    "confidence": conf, "meanScore": ms})
    return out


def write_celltypes(er, d, markers: Path | None = None, log=print) -> dict | None:
    """Write umap/celltypes.json + umap/celltype.bin for one layer (ExprReader `er`, output dir `d`).

    Needs umap/index.json from sectionumap.write_section_umap. Returns the manifest's
    "celltypes" entry, or None when skipped.
    """
    key = er.meta.get("key")
    ud = d / "umap"
    ip = ud / "index.json"
    if not ip.exists():
        log(f"[celltypes] {key}: no umap/index.json, skipped")
        return None
    t0 = time.time()
    index = json.loads(ip.read_text(encoding="utf-8"))
    ms, rel = load_markers(markers)
    if ms.get("modality", "rna") != "rna":
        log(f"[celltypes] {key}: marker list {rel} is for {ms.get('modality')!r}, this layer is RNA - skipped "
            f"(pass an RNA list with --markers, e.g. markers/nerve.json)")
        return None
    matched = match_genes(ms["types"], er.genes)
    names = [t["name"] for t in matched]
    if len(names) >= NO_TYPE:
        raise SystemExit(f"[celltypes] {len(names)} types: celltype.bin holds at most {NO_TYPE}")

    xy = er.xy
    on_image = ~np.isnan(xy[:, 0])
    # positions.bin keeps only bins with a position, in order: that is the row space
    row_of = np.cumsum(on_image) - 1
    n_rows = int(on_image.sum())
    coords = np.fromfile(ud / "coords.bin", dtype="<f4").reshape(-1, 2)
    if len(coords) != n_rows:
        raise SystemExit(f"[celltypes] {key}: coords.bin has {len(coords)} rows, positions.bin {n_rows}")
    emb_rows = np.flatnonzero(~np.isnan(coords[:, 0]))

    # library size over ALL genes, as sectionumap/hvgchunks: one pass over the count matrix
    indices = np.asarray(er.g["indices"])
    data = np.asarray(er.g["data"], dtype=np.float32)
    lib_all = np.bincount(indices, weights=data, minlength=len(xy)).astype(np.float32)
    del indices, data
    bins = np.flatnonzero(on_image & (lib_all >= MIN_UMI))
    if not np.array_equal(row_of[bins], emb_rows):
        raise SystemExit(f"[celltypes] {key}: embedded rows in coords.bin do not match bins >= {MIN_UMI} UMI")
    pos_of = np.full(len(xy), -1, np.int64)
    pos_of[bins] = np.arange(len(bins))
    lib = lib_all[bins].astype(np.float64)

    cnt: dict[int, np.ndarray] = {}
    for col in sorted({c for t in matched for c in t["cols"]}):
        a, b = int(er.indptr[col]), int(er.indptr[col + 1])
        v = np.zeros(len(bins))
        if b > a:
            ids = np.asarray(er.g["indices"][a:b])
            vals = np.asarray(er.g["data"][a:b], dtype=np.float64)
            p = pos_of[ids]
            keep = p >= 0
            v[p[keep]] = vals[keep]
        cnt[col] = v
    score, scored = score_bins(cnt, lib, matched)
    best = best_type(score, scored)

    types = []
    for t, s in zip(matched, scored):
        e = {"name": t["name"], "genes": t["used"],
             "missing": t["missing"] + [{"gene": g, "reason": "rare"} for g in t["rare"]], "scored": bool(s)}
        if not s:
            e["reason"] = (f"{len(t['used'])} of {len(t['used']) + len(e['missing'])} genes measurable "
                           f"here, needs {MIN_GENES}")
        types.append(e)
    clusterings = {}
    for c in index["clusterings"]:
        cl = np.fromfile(ud / c["file"], dtype="<u2")
        lab = cl[emb_rows].astype(np.int64)
        if len(cl) != n_rows or (lab == NO_CLUSTER).any() or int(lab.max()) + 1 != c["n"]:
            raise SystemExit(f"[celltypes] {key}: {c['file']} does not match coords.bin / index.json")
        clusterings[f"{c['resolution']:g}"] = summarise_clusters(lab, best, score, scored, names)

    out = {"markerSet": {"name": ms.get("name"), "version": ms.get("version"), "file": rel},
           "types": types, "rule": rule_text(), "clusterings": clusterings}
    deid.assert_clean(out, "umap/celltypes.json")
    txt = json.dumps(out, separators=(",", ":"))
    if _BARCODE.search(txt):
        raise SystemExit("[celltypes] REFUSING to write celltypes.json: it contains a barcode")
    ct = np.full(n_rows, NOT_EMBEDDED, np.uint8)
    ct[emb_rows] = best
    for name, buf in (("celltype.bin", ct.tobytes()), ("celltypes.json", txt.encode("utf-8"))):
        tmp = ud / (name + ".part")
        tmp.write_bytes(buf)
        tmp.replace(ud / name)
    n_typed = int((best < NO_TYPE).sum())
    log(f"[celltypes] {key}: {int(scored.sum())}/{len(names)} types scored, {n_typed:,}/{len(bins):,} bins "
        f"typed, " + ", ".join(f"r{r}: {sum(x['confidence'] == 'high' for x in v)} high"
                               for r, v in clusterings.items()) + f"  {time.time() - t0:.1f}s")
    return {"nTypes": len(names), "markerSet": ms.get("name")}


# ---------------------------------------------------------------------------------------
# Protein rule: grid squares of a multiplexed immunofluorescence image (protgrid.py)
# ---------------------------------------------------------------------------------------
DEFAULT_PROTEIN_MARKERS = Path(__file__).resolve().parent / "markers" / "protein_generic.json"
PROTEIN_MIN_CHANNELS = 1     # a list without "minChannels" (akoya_nerve.json): one channel per type is enough
PROTEIN_BIN_MIN = 0.5
ARTEFACT_SHARE = 0.5         # a cluster with more than this share of "everything bright" squares ...
ARTEFACT_NAME = "Artefact? (fold / edge)"   # ... is suggested as this, whatever its markers say
FLAG_BRIGHT = 1              # = protgrid.FLAG_BRIGHT

# Channel-name cleanup for matching (= lasso_worker.js normName): a leading index ("01_",
# "12 - "), then trailing fluorophore / block tags after a separator ("-488", "_Opal520",
# " (Atto550)", " AF647", " FFPE"), then case, spaces, "-", "_" and "." are ignored.
_LEAD = re.compile(r"^\s*\d+\s*[_\-.:)]\s*")
_TAG = re.compile(r"(?:[\s_\-]+|\s*\()(?:opal\s*\d+|atto\s*\d+|af\s*\d+|alexa(?:\s*fluor)?\s*\d+|cy\s*\d(?:\.\d)?|"
                  r"fitc|tritc|pe|apc|ffpe|405|488|514|520|532|546|555|568|570|594|620|647|650|660|680|690|700|750|780|790)\)?\s*$", re.I)
_SEP = re.compile(r"[\s_\-.]")


def norm_name(s: str) -> str:
    s = _LEAD.sub("", str(s))
    while True:
        t = _TAG.sub("", s)
        if t == s or not t.strip():
            break
        s = t
    return _SEP.sub("", s).lower()


def match_channels(types: list[dict], channels: list[str], aliases: dict | None = None) -> list[dict]:
    """match_genes for protein panels: names compared after norm_name, each list name also
    under its aliases (the marker set's "aliases": {name: [other names]}); first channel wins."""
    first: dict[str, int] = {}
    for i, c in enumerate(channels):
        first.setdefault(norm_name(c), i)
    aliases = aliases or {}
    out = []
    for t in types:
        cols, found, missing = [], [], []
        for g in t["genes"]:
            j = next((first[k] for k in [norm_name(g)] + [norm_name(a) for a in aliases.get(g, [])] if k in first), None)
            if j is None or j in cols:
                missing.append({"gene": g, "reason": "absent"})
            else:
                cols.append(j)
                found.append(channels[j])
        out.append({"name": t["name"], "cols": cols, "genes": found, "missing": missing})
    return out


def type_min(ms: dict, t: dict) -> int:
    """Channels a type needs on the panel to be scored: its own "min", else the list's
    "minChannels", else PROTEIN_MIN_CHANNELS."""
    return int(t.get("min", ms.get("minChannels", PROTEIN_MIN_CHANNELS)))


def rule_text_protein(ms: dict | None = None) -> str:
    mins = sorted({type_min(ms or {}, t) for t in (ms or {}).get("types", [])}) or [PROTEIN_MIN_CHANNELS]
    need = (f">= {mins[0]} channel{'s' if mins[0] > 1 else ''} present" if len(mins) == 1 else
            f">= {ms.get('minChannels', PROTEIN_MIN_CHANNELS)} channels present (some types, marked with their own "
            f"minimum, need fewer)")
    return (f"Protein rule. Each square: x = its mean intensity per marker channel with the slide's background "
            f"taken off, divided by the channel's 99.5th percentile over the tissue and put through asinh "
            f"(0 = background, 1 = the 99.5th percentile); z-scored over the squares (population sd) and "
            f"averaged over the type's channels (no log / CP10k, no rare-channel drop). "
            f"The square's type is the best-scoring one if its score >= {PROTEIN_BIN_MIN:g}, else none. "
            f"Channel names are matched ignoring case, number prefixes, fluorophore tags and listed aliases. "
            f"A channel the panel lacks is 'absent'; one that is flat over the tissue is left out as 'flat'. "
            f"A type needs {need} to be scored. "
            f"Per cluster, pct = share of its squares of each type (none counts in the denominator). "
            f"high: top pct >= {HIGH:g} and >= 2x the second; medium: >= {MEDIUM:g}; low: >= {LOW:g} "
            f"(shown as 'maybe'); otherwise Unclear. A cluster with more than {ARTEFACT_SHARE:.0%} of its squares "
            f"flagged 'everything bright' (folds, bubbles, edges) is suggested as '{ARTEFACT_NAME}'.")


def flag_artefacts(summary: list[dict], lab: np.ndarray, flags) -> list[dict]:
    """Per cluster: the share of its squares flagged everything-bright; above ARTEFACT_SHARE the
    suggestion becomes ARTEFACT_NAME (confidence "low": shown as a warning, never auto-accepted)."""
    if flags is None:
        return summary
    br = (np.asarray(flags) & FLAG_BRIGHT).astype(bool)
    for s in summary:
        m = lab == s["cluster"]
        share = float(br[m].mean()) if m.any() else 0.0
        s["artefact"] = round(share, 3)
        if share > ARTEFACT_SHARE:
            s.update(suggest=ARTEFACT_NAME, confidence="low", artefactRule=True)
    return summary


def write_celltypes_protein(d, xp: np.ndarray, channels: list[str], usable, markers: Path | None = None,
                            key: str = "", log=print, flags=None) -> dict | None:
    """umap/celltypes.json + umap/celltype.bin for a protein layer (protgrid.py).

    `xp` = squares x channels of x (every positions.bin row is embedded), `channels`
    their names, `usable[j]` False for a channel that is flat over the tissue, `flags`
    protgrid's qc_flags (uint8 per square) for the artefact rule, or None.
    Needs umap/index.json from sectionumap.write_feature_umap.
    """
    d = Path(d)
    ud = d / "umap"
    ip = ud / "index.json"
    if not ip.exists():
        log(f"[celltypes] {key}: no umap/index.json, skipped")
        return None
    t0 = time.time()
    index = json.loads(ip.read_text(encoding="utf-8"))
    ms, rel = load_markers(markers or DEFAULT_PROTEIN_MARKERS)
    if ms.get("modality") != "protein":
        log(f"[celltypes] {key}: marker list {rel} has no \"modality\": \"protein\" - skipped "
            f"(pass a protein list with --markers, e.g. markers/akoya_nerve.json)")
        return None
    matched = match_channels(ms["types"], channels, ms.get("aliases"))
    mins = [type_min(ms, t) for t in ms["types"]]
    names = [t["name"] for t in matched]
    if len(names) >= NO_TYPE:
        raise SystemExit(f"[celltypes] {len(names)} types: celltype.bin holds at most {NO_TYPE}")
    n = len(xp)
    if index.get("nBins") != n or index.get("nEmbedded") != n:
        raise SystemExit(f"[celltypes] {key}: umap/index.json does not cover the {n} squares")

    z = {}
    for j in sorted({c for t in matched for c in t["cols"]}):
        x = xp[:, j].astype(np.float64)
        sd = x.std()
        if usable[j] and sd > 0:
            z[j] = (x - x.mean()) / sd
    score = np.full((n, len(matched)), np.nan)
    scored = np.zeros(len(matched), bool)
    types = []
    for ti, t in enumerate(matched):
        cols = [c for c in t["cols"] if c in z]
        used = [g for g, c in zip(t["genes"], t["cols"]) if c in z]
        flat = [g for g, c in zip(t["genes"], t["cols"]) if c not in z]
        if len(cols) >= mins[ti]:
            score[:, ti] = np.mean([z[c] for c in cols], axis=0)
            scored[ti] = True
        e = {"name": t["name"], "genes": used,
             "missing": t["missing"] + [{"gene": g, "reason": "flat"} for g in flat], "scored": bool(scored[ti]),
             "min": mins[ti]}
        if not scored[ti]:
            e["reason"] = (f"{len(used)} of {len(used) + len(e['missing'])} channels measurable here, "
                           f"needs {mins[ti]}")
        types.append(e)
    best = best_type(score, scored, PROTEIN_BIN_MIN)

    clusterings = {}
    for c in index["clusterings"]:
        lab = np.fromfile(ud / c["file"], dtype="<u2").astype(np.int64)
        if len(lab) != n or (lab == NO_CLUSTER).any() or int(lab.max()) + 1 != c["n"]:
            raise SystemExit(f"[celltypes] {key}: {c['file']} does not match index.json")
        clusterings[f"{c['resolution']:g}"] = flag_artefacts(summarise_clusters(lab, best, score, scored, names),
                                                            lab, flags)

    out = {"markerSet": {"name": ms.get("name"), "version": ms.get("version"), "file": rel},
           "modality": "protein", "types": types, "rule": rule_text_protein(ms), "clusterings": clusterings}
    deid.assert_clean(out, "umap/celltypes.json")
    txt = json.dumps(out, separators=(",", ":"))
    ct = best.astype(np.uint8)                       # every square is embedded: no 255 rows
    for name, buf in (("celltype.bin", ct.tobytes()), ("celltypes.json", txt.encode("utf-8"))):
        tmp = ud / (name + ".part")
        tmp.write_bytes(buf)
        tmp.replace(ud / name)
    n_typed = int((best < NO_TYPE).sum())
    log(f"[celltypes] {key}: protein, {int(scored.sum())}/{len(names)} types scored, {n_typed:,}/{n:,} squares "
        f"typed, " + ", ".join(f"r{r}: {sum(x['confidence'] == 'high' for x in v)} high"
                               for r, v in clusterings.items()) + f"  {time.time() - t0:.1f}s")
    return {"nTypes": len(names), "markerSet": ms.get("name"), "modality": "protein"}
