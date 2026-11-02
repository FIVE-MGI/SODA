"""Tier A: assign every expression bin to the overlay region that contains it.

This runs ONCE, at export, so the viewer never does geometry. For each bin centre
in `positions.bin` order we find the smallest overlay polygon that contains it and
write a parallel `region.bin` (uint16 region id, 0 = none). The same ids are written
back into the exported GeoJSON as `properties.__region_id`, so a click in the viewer
reads the id straight off the feature with no matching.

The link is frame-free: bin centres (expr.py, `pxl_*_in_fullres`) and the overlay
polygons are already in the same full-res image pixel grid, so no registration step.

Containment rule: assign to the SMALLEST-area containing polygon, and treat any
whole-slide outline as "no region" - otherwise a `hemi`-style outline that wraps
every other region swallows every bin. Outlines are detected by name and by area.
"""
from __future__ import annotations

import numpy as np

DEFAULT_OUTLINE_NAMES = ("hemi",)
OUTLINE_AREA_FRAC = 0.85          # a feature >= this fraction of the overlay extent is an outline


def _rings_of(geom: dict) -> list:
    """GeoJSON geometry -> list of polygons, each [exterior, hole, ...] of coord lists."""
    t = (geom or {}).get("type")
    c = (geom or {}).get("coordinates")
    if t == "Polygon":
        return [c]
    if t == "MultiPolygon":
        return c
    return []


def _shoelace(ring: np.ndarray) -> float:
    x, y = ring[:, 0], ring[:, 1]
    return 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def _feat_name(props: dict) -> str:
    if not props:
        return ""
    cl = props.get("classification")
    if isinstance(cl, dict) and cl.get("name"):
        return str(cl["name"])
    for k in ("name", "label", "region", "id"):
        if props.get(k):
            return str(props[k])
    return ""


def load_features(geojson: dict) -> list:
    """Parse a FeatureCollection into [{feature, name, polys, area, bbox} | None].

    `None` marks a non-polygon feature (points, lines) so indices line up with
    geojson['features'] and the id written back lands on the right feature.
    """
    out = []
    for f in (geojson.get("features") or []):
        polys = _rings_of(f.get("geometry") or {})
        parsed, area = [], 0.0
        minx = miny = np.inf
        maxx = maxy = -np.inf
        for poly in polys:
            if not poly:
                continue
            ext = np.asarray(poly[0], dtype=np.float64)
            if ext.ndim != 2 or len(ext) < 3:
                continue
            holes = [np.asarray(h, dtype=np.float64) for h in poly[1:] if len(h) >= 3]
            parsed.append((ext, holes))
            area += _shoelace(ext) - sum(_shoelace(h) for h in holes)
            pts = np.vstack([ext, *holes]) if holes else ext
            minx, miny = min(minx, pts[:, 0].min()), min(miny, pts[:, 1].min())
            maxx, maxy = max(maxx, pts[:, 0].max()), max(maxy, pts[:, 1].max())
        if not parsed:
            out.append(None)
            continue
        out.append({"feature": f, "name": _feat_name(f.get("properties")),
                    "polys": parsed, "area": area, "bbox": (minx, miny, maxx, maxy)})
    return out


def assign_region_ids(feats: list, outline_names=DEFAULT_OUTLINE_NAMES,
                      outline_area_frac: float = OUTLINE_AREA_FRAC) -> dict:
    """feature index -> region id (1..K), with outlines mapped to 0.

    Ids are stable across every expression layer of the dataset: they depend only on
    feature order and the overlay's own geometry, never on a layer's bins.
    """
    polys = [(i, fe) for i, fe in enumerate(feats) if fe]
    if not polys:
        return {}
    xs = [fe["bbox"][0] for _, fe in polys] + [fe["bbox"][2] for _, fe in polys]
    ys = [fe["bbox"][1] for _, fe in polys] + [fe["bbox"][3] for _, fe in polys]
    extent = max((max(xs) - min(xs)) * (max(ys) - min(ys)), 1.0)
    onames = {n.lower() for n in outline_names}
    feat_id, nxt = {}, 1
    for i, fe in polys:
        outline = fe["name"].lower() in onames or fe["area"] >= outline_area_frac * extent
        if outline:
            feat_id[i] = 0
        else:
            feat_id[i] = nxt
            nxt += 1
    return feat_id


def _inside(pts: np.ndarray, polys: list) -> np.ndarray:
    """Even-odd point-in-polygon for (N,2) pts against one feature's polygons+holes."""
    res = np.zeros(len(pts), dtype=bool)
    for ext, holes in polys:
        m = _ray(pts, ext)
        for h in holes:
            m &= ~_ray(pts, h)
        res |= m
    return res


def _ray(pts: np.ndarray, ring: np.ndarray) -> np.ndarray:
    x, y = pts[:, 0], pts[:, 1]
    xi, yi = ring[:, 0], ring[:, 1]
    xj, yj = np.roll(xi, -1), np.roll(yi, -1)
    inside = np.zeros(len(pts), dtype=bool)
    for k in range(len(ring)):
        a, b, c, d = xi[k], yi[k], xj[k], yj[k]
        if b == d:                                   # horizontal edge: no crossing
            continue
        cond = (b > y) != (d > y)
        xint = a + (c - a) * (y - b) / (d - b)
        inside ^= cond & (x < xint)
    return inside


def index_bins(xy: np.ndarray, feats: list, feat_id: dict):
    """Region id per bin (uint16, 0 = none) + the layer's region table.

    Assigns smallest-area-first so an inner region wins over a larger one that also
    contains the bin. A bbox prefilter keeps the work near the bins each region covers.
    """
    n = len(xy)
    ids = np.zeros(n, dtype=np.uint16)
    if not n or not feat_id:
        return ids, []
    order = sorted((i for i, rid in feat_id.items() if rid), key=lambda i: feats[i]["area"])
    free = np.ones(n, dtype=bool)
    gx, gy = xy[:, 0], xy[:, 1]
    for i in order:
        fe = feats[i]
        minx, miny, maxx, maxy = fe["bbox"]
        cand = free & (gx >= minx) & (gx <= maxx) & (gy >= miny) & (gy <= maxy)
        idx = np.where(cand)[0]
        if not len(idx):
            continue
        hit = idx[_inside(xy[idx], fe["polys"])]
        ids[hit] = feat_id[i]
        free[hit] = False
    table = []
    for i in order:
        rid = feat_id[i]
        table.append({"id": int(rid), "name": feats[i]["name"] or f"region {rid}",
                      "nBins": int(np.count_nonzero(ids == rid))})
    table.sort(key=lambda r: r["id"])
    return ids, [r for r in table if r["nBins"] > 0]


def annotate_geojson(geojson: dict, feats: list, feat_id: dict) -> dict:
    """Write the assigned region id into each feature's properties, in place."""
    for i, fe in enumerate(feats):
        if fe is None:
            continue
        f = fe["feature"]
        props = f.setdefault("properties", {}) or {}
        f["properties"] = props
        rid = feat_id.get(i, 0)
        if rid:
            props["__region_id"] = int(rid)
        else:
            props["__region_outline"] = True
    return geojson
