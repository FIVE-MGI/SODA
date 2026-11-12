#!/usr/bin/env python
r"""Make a small SYNTHETIC GeoMx DSP study for the smoke test (tests/smoke_test.py).

Nothing here comes from a tissue, an instrument, a probe kit or a run: every name, count,
pixel and ROI is made up and drawn by a seeded generator, so the same --seed gives the same
files.

    python examples/make_geomx_demo.py --out demo/raw

writes <out>/geomx/:

    dcc/DSP-1001000000001-A-A01.dcc ... one DCC per well: Header, Scan_Attributes, NGS_Processing_Attributes,
                             Code_Summary (RTS_ID,count; zero counts absent; one code no PKC knows).
                             A01 = No Template Control; 13 AOIs of two scans; 2 AOIs of a third scan
                             that has no image.
    pkc/Demo_Panel_v1.0.pkc  18 targets: single-probe genes, 4 genes with 3 probes, one control target,
                             one negative-probe target (6 probes). Names include "H2-T23" (a gene symbol
                             that looks like a sample ID) and "TRB1/2" (a name that needs encoding).
    pkc/Demo_Addon_v1.0.pkc  a second panel: 3 targets (one named like a target of the first panel)
                             and its own negative-probe target (3 probes).
    annotations.xlsx         the annotation workbook: Sample_ID, slide name, scan name, panel, roi,
                             segment, aoi, area, tissue, genotype, nuclei. The slide names carry a
                             made-up sample ID, to check that none reaches the export.
    annotations.csv          the same table as CSV
    Demo_LabWorksheet.txt    the same as the readout package writes it: a preamble, then the table
                             with Sample_ID .. area and ="001"-style ROI numbers
    segment_properties.xlsx  a DSP Data Analysis style sheet "SegmentProperties" without Sample_ID
                             (SlideName, ScanLabel, ROILabel, SegmentLabel, AOINucleiCount, Region)
    demoA.ome.tiff, demoB.ome.tiff   two scans (2 channels, uint16, 0.5 µm/px) whose OME-XML holds the
                             ROIs: Label, Polygon / Ellipse / Rectangle, and one bit Mask per segment.
                             ROI 003 of demoA has two segments; ROI 005 is a "C" whose centroid lies
                             outside it.
    <out>/truth_geomx.json   what the test checks against: target values (collapsed here by plain
                             loops), background, LOQ, positions, areas, groups.
"""
from __future__ import annotations

import argparse
import base64
import json
import math
import zipfile
from pathlib import Path
from xml.sax.saxutils import escape

import numpy as np

PX_UM = 0.5
PLATE = "1001000000001"
FAKE_ID = "INT" + "9951"            # a made-up sample ID (never a real one), built so this file holds no ID  (deid:ok)
SCANS = {"demoA": (1200, 900), "demoB": (1000, 800)}
SLIDE = {"demoA": f"{FAKE_ID} block A", "demoB": f"{FAKE_ID} block B", "demoC": "block C"}


def targets():
    """[(name, code class, n probes)] of the two panels."""
    main = [(f"Gene{i}", "Endogenous01", 1) for i in range(1, 9)]
    main += [("H2-T23", "Endogenous01", 1), ("TRB1/2", "Endogenous01", 1), ("Shared", "Endogenous01", 1)]
    main += [(f"Multi{i}", "Endogenous01", 3) for i in range(1, 5)]
    main += [("mt-Demo1", "Endogenous01", 1), ("HK1", "Control01", 2), ("NegProbe-DEMO", "Negative01", 6)]
    addon = [("Addon1", "Endogenous02", 1), ("Addon2", "Endogenous02", 4), ("Shared", "Endogenous02", 1),
             ("Addon Negative", "Negative02", 3)]
    return {"Demo Panel": main, "Demo Addon": addon}


def aois():
    """One dict per AOI: well, scan, roi, segment, shape, groups."""
    rows = []
    def add(scan, roi, seg, shape, tissue, geno):
        rows.append({"scan": scan, "roi": roi, "segment": seg, "shape": shape, "tissue": tissue, "genotype": geno})
    add("demoA", "001", "Full ROI", ("ellipse", 200, 200, 90, 90), "cortex", "WT")
    add("demoA", "002", "Full ROI", ("rect", 380, 120, 180, 150), "cortex", "KO")
    add("demoA", "003", "PanCK+", ("rect", 650, 100, 240, 160), "tumour", "WT")       # two segments of one ROI
    add("demoA", "003", "PanCK-", ("rect", 650, 100, 240, 160), "tumour", "WT")
    add("demoA", "004", "Full ROI", ("poly", [(150, 500), (330, 480), (360, 640), (240, 720), (130, 650)]), "tumour", "KO")
    add("demoA", "005", "Full ROI", ("poly", [(500, 450), (760, 450), (760, 520), (580, 520), (580, 680), (760, 680),
                                              (760, 750), (500, 750)]), "stroma", "WT")   # a "C": centroid outside
    add("demoA", "006", "Full ROI", ("ellipse", 980, 600, 120, 80), "stroma", "KO")
    add("demoB", "001", "Full ROI", ("ellipse", 180, 180, 100, 100), "cortex", "WT")
    add("demoB", "002", "Full ROI", ("rect", 400, 80, 200, 200), "cortex", "KO")
    add("demoB", "003", "Full ROI", ("poly", [(700, 100), (900, 150), (880, 330), (720, 300)]), "tumour", "WT")
    add("demoB", "004", "Full ROI", ("ellipse", 250, 550, 140, 110), "tumour", "KO")
    add("demoB", "005", "Full ROI", ("rect", 520, 450, 160, 220), "stroma", "WT")
    add("demoB", "006", "Full ROI", ("ellipse", 820, 600, 90, 120), "stroma", "KO")
    add("demoC", "001", "Full ROI", None, "cortex", "WT")                               # a scan with no image
    add("demoC", "002", "Full ROI", None, "tumour", "KO")
    for i, r in enumerate(rows):
        w = i + 2                                                                         # A01 is the NTC
        r["well"] = f"{'ABCDEFGH'[(w - 1) // 12]}{(w - 1) % 12 + 1:02d}"
        r["sample_id"] = f"DSP-{PLATE}-A-{r['well']}"
    return rows


def raster(shape, W, H) -> np.ndarray:
    """Boolean mask (H, W) of a shape, pixel centres tested."""
    yy, xx = np.mgrid[0:H, 0:W]
    x, y = xx + 0.5, yy + 0.5
    kind = shape[0]
    if kind == "ellipse":
        _, cx, cy, rx, ry = shape
        return ((x - cx) / rx) ** 2 + ((y - cy) / ry) ** 2 <= 1.0
    if kind == "rect":
        _, x0, y0, w, h = shape
        return (x >= x0) & (x < x0 + w) & (y >= y0) & (y < y0 + h)
    pts = shape[1]
    m = np.zeros((H, W), bool)
    n = len(pts)
    for i in range(n):
        (x1, y1), (x2, y2) = pts[i], pts[(i + 1) % n]
        if y1 == y2:
            continue
        cond = (y1 > y) != (y2 > y)
        m ^= cond & (x < x1 + (x2 - x1) * (y - y1) / (y2 - y1))
    return m


def mask_xml(m: np.ndarray, text: str, sid: int) -> tuple[str, dict]:
    ys, xs = np.nonzero(m)
    x0, x1, y0, y1 = xs.min(), xs.max() + 1, ys.min(), ys.max() + 1
    sub = m[y0:y1, x0:x1]
    b64 = base64.b64encode(np.packbits(sub.ravel(), bitorder="big").tobytes()).decode()
    xml = (f'<Mask FillColor="-1" Height="{y1 - y0}.0" ID="Shape:{sid}" Text="{escape(text)}" Width="{x1 - x0}.0" '
           f'X="{x0}.0" Y="{y0}.0"><BinData BigEndian="false" Compression="none" Length="{len(b64)}">{b64}</BinData></Mask>')
    return xml, {"areaPx": int(m.sum()), "cx": float(xs.mean()) + 0.5, "cy": float(ys.mean()) + 0.5}


def write_scan(path: Path, name: str, rows: list[dict], rng) -> dict:
    """A 2-channel OME-TIFF with the ROIs in its OME-XML. Returns {(roi, segment): mask facts}."""
    import tifffile
    W, H = SCANS[name]
    img = rng.integers(20, 60, size=(2, H, W)).astype(np.uint16)
    facts, rois_xml, refs, sid = {}, [], [], 0
    by_roi: dict[str, list[dict]] = {}
    for r in rows:
        by_roi.setdefault(r["roi"], []).append(r)
    for k, (roi, members) in enumerate(by_roi.items()):
        shape = members[0]["shape"]
        full = raster(shape, W, H)
        img[0][full] += 900                                    # "nuclei" where the ROI is
        img[1][full] += 300
        parts = [f'<Label FontSize="32" ID="Shape:{sid}" Text="{roi}" X="0.0" Y="0.0"/>']
        sid += 1
        for i, r in enumerate(members):
            m = full
            if len(members) > 1:                               # segments: left / right part of the ROI
                xs = np.nonzero(full.any(0))[0]
                cut = xs.min() + (xs.max() + 1 - xs.min()) * 0.4
                half = (np.arange(W) + 0.5 < cut)[None, :]
                m = full & (half if i == 0 else ~half)
            xml, f = mask_xml(m, r["segment"], sid)
            sid += 1
            parts.append(xml)
            facts[(roi, r["segment"])] = f
        if shape[0] == "ellipse":
            parts.append(f'<Ellipse ID="Shape:{sid}" RadiusX="{shape[3]}" RadiusY="{shape[4]}" X="{shape[1]}" Y="{shape[2]}"/>')
        elif shape[0] == "rect":
            parts.append(f'<Rectangle ID="Shape:{sid}" Height="{shape[4]}" Width="{shape[3]}" X="{shape[1]}" Y="{shape[2]}"/>')
        else:
            parts.append(f'<Polygon ID="Shape:{sid}" Points="{" ".join(f"{x},{y}" for x, y in shape[1])}"/>')
        sid += 1
        rois_xml.append(f'<ROI ID="ROI:{k}"><Union>{"".join(parts)}</Union></ROI>')
        refs.append(f'<ROIRef ID="ROI:{k}"/>')
    ch = ('<Channel ID="Channel:0:0" Name="FITC/525nm" Fluor="SYTO 13" SamplesPerPixel="1"/>'
          '<Channel ID="Channel:0:1" Name="Cy3/568nm" Fluor="Alexa 532" SamplesPerPixel="1"/>')
    td = "".join(f'<TiffData FirstC="{c}" FirstT="0" FirstZ="0" IFD="{c}"/>' for c in range(2))
    xml = ('<?xml version="1.0" encoding="UTF-8"?>'
           '<OME xmlns="http://www.openmicroscopy.org/Schemas/OME/2016-06">'
           f'<Image ID="Image:0" Name="{name}"><Description>SYNTHETIC demo scan - not from an instrument</Description>'
           f'<Pixels ID="Pixels:0" DimensionOrder="XYCZT" Type="uint16" SizeX="{W}" SizeY="{H}" SizeC="2" SizeZ="1" SizeT="1" '
           f'PhysicalSizeX="{PX_UM}" PhysicalSizeXUnit="&#181;m" PhysicalSizeY="{PX_UM}" PhysicalSizeYUnit="&#181;m">{ch}{td}</Pixels>'
           f'{"".join(refs)}</Image>{"".join(rois_xml)}</OME>')
    tifffile.imwrite(str(path), img, photometric="minisblack", tile=(256, 256), compression="zlib",
                     metadata=None, description=xml, software="make_geomx_demo.py (synthetic)")
    return facts


def write_xlsx(path: Path, sheet: str, rows: list[list]):
    """A minimal workbook: text as shared strings (one cell as an inline string), numbers as numbers."""
    strings, idx = [], {}

    def cell(r, c, v):
        ref = f"{'ABCDEFGHIJKLMNOPQRSTUVWXYZ'[c]}{r + 1}"
        if isinstance(v, (int, float)):
            return f'<c r="{ref}"><v>{v!r}</v></c>'
        v = str(v)
        if v == "":
            return ""                                           # an empty cell is simply absent, as Excel writes it
        if r == 2 and c == 1:
            return f'<c r="{ref}" t="inlineStr"><is><t>{escape(v)}</t></is></c>'
        if v not in idx:
            idx[v] = len(strings)
            strings.append(v)
        return f'<c r="{ref}" t="s"><v>{idx[v]}</v></c>'
    body = "".join(f'<row r="{r + 1}">{"".join(cell(r, c, v) for c, v in enumerate(row))}</row>' for r, row in enumerate(rows))
    ns = 'xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"'
    rel = 'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"'
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("[Content_Types].xml",
                   '<?xml version="1.0" encoding="UTF-8"?><Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
                   '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
                   '<Default Extension="xml" ContentType="application/xml"/>'
                   '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
                   '<Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
                   '<Override PartName="/xl/sharedStrings.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sharedStrings+xml"/>'
                   '</Types>')
        z.writestr("_rels/.rels",
                   '<?xml version="1.0" encoding="UTF-8"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
                   '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>'
                   '</Relationships>')
        z.writestr("xl/workbook.xml", f'<?xml version="1.0" encoding="UTF-8"?><workbook {ns} {rel}><sheets>'
                                      f'<sheet name="{escape(sheet)}" sheetId="1" r:id="rId1"/></sheets></workbook>')
        z.writestr("xl/_rels/workbook.xml.rels",
                   '<?xml version="1.0" encoding="UTF-8"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
                   '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet1.xml"/>'
                   '<Relationship Id="rId2" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/sharedStrings" Target="sharedStrings.xml"/>'
                   '</Relationships>')
        z.writestr("xl/worksheets/sheet1.xml", f'<?xml version="1.0" encoding="UTF-8"?><worksheet {ns}><sheetData>{body}</sheetData></worksheet>')
        z.writestr("xl/sharedStrings.xml", f'<?xml version="1.0" encoding="UTF-8"?><sst {ns} count="{len(strings)}" uniqueCount="{len(strings)}">'
                                           + "".join(f"<si><t>{escape(s)}</t></si>" for s in strings) + "</sst>")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True, help="folder to write geomx/ and truth_geomx.json into")
    ap.add_argument("--seed", type=int, default=11)
    args = ap.parse_args(argv)
    out = Path(args.out)
    d = out / "geomx"
    (d / "dcc").mkdir(parents=True, exist_ok=True)
    (d / "pkc").mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    rows = aois()

    # ---- probe kits ------------------------------------------------------------------------
    panels = targets()
    probes, rts_n = {}, 9000000                       # (panel, target) -> [RTS_ID]
    for pname, tl in panels.items():
        T = []
        for name, cc, k in tl:
            ids = []
            for i in range(k):
                rts_n += 1
                ids.append(f"RTS{rts_n:07d}")
            probes[(pname, name)] = ids
            T.append({"AnalyteType": "RNA", "DisplayName": name, "CodeClass": cc,
                      "Probes": [{"RTS_ID": r, "DisplayName": f"{name}_{i + 1:02d}", "ProbeID": f"demo{r[-4:]}",
                                  "SystematicName": [name], "Accession": [], "GeneID": []} for i, r in enumerate(ids)]})
        j = {"Name": pname, "Version": 1.0, "AnalyteType": "RNA", "MinArea": 16000, "MinNuclei": 200,
             "Note": "SYNTHETIC demo panel - not a NanoString / Bruker probe kit", "ProbeGroups": [], "Targets": T}
        (d / "pkc" / f"{pname.replace(' ', '_')}_v1.0.pkc").write_text(json.dumps(j, indent=1), encoding="utf-8")

    # ---- counts: tissue and genotype shift some genes; multi-probe targets get zeros ---------
    tissue_gene = {"cortex": "Gene1", "tumour": "Gene2", "stroma": "Gene3"}
    counts = {}                                       # sample -> {RTS: count}
    for r in rows:
        depth = float(rng.uniform(0.6, 1.8))
        c = {}
        for (pname, name), ids in probes.items():
            neg = name.startswith("NegProbe") or name.endswith("Negative")
            base = 3.0 if neg else 40.0
            if name == tissue_gene[r["tissue"]]:
                base = 600.0
            if name == "Gene4" and r["genotype"] == "KO":
                base = 4.0
            if name == "Gene5":
                base = 0.0                            # a target never counted: no gene file
            if name.startswith("Multi"):
                base = 6.0
            for k, rid in enumerate(ids):
                v = int(rng.poisson(base * depth))
                if name == "Multi1" and k == 0:
                    v = 0                             # one probe always 0: the + 0.5 rule on every AOI
                if name == "Multi2":
                    v = 0                             # all probes 0: geometric mean of 0.5s = 0.5
                c[rid] = v
        c["RTS9999999"] = 7                           # a code in no PKC
        counts[r["sample_id"]] = c
    ntc_id = f"DSP-{PLATE}-A-A01"
    counts[ntc_id] = {probes[("Demo Panel", "Gene1")][0]: 2, probes[("Demo Panel", "Gene2")][0]: 1}
    seq = {}
    for sid, c in counts.items():
        ded = sum(c.values())
        aligned = int(ded * float(rng.uniform(1.8, 3.2))) + 5
        raw = int(aligned / float(rng.uniform(0.80, 0.95))) + 5
        seq[sid] = {"raw": raw, "trimmed": int(raw * 0.99), "stitched": int(raw * 0.97), "aligned": aligned, "dedup": ded}
        code = "\n".join(f"{k},{v}" for k, v in sorted(c.items()) if v > 0)
        (d / "dcc" / f"{sid}.dcc").write_text(
            f'<Header>\nFileVersion,0.02\nSoftwareVersion,"synthetic demo"\nDate,2026-01-01\n</Header>\n\n'
            f"<Scan_Attributes>\nID,{sid}\nPlate_ID,{PLATE}\nWell,{sid[-3:]}\n</Scan_Attributes>\n\n"
            f"<NGS_Processing_Attributes>\nSeqSetId,DEMO:1:SYNTHETIC\nRaw,{raw}\nTrimmed,{seq[sid]['trimmed']}\n"
            f"Stitched,{seq[sid]['stitched']}\nAligned,{aligned}\numiQ30,0.9900\nrtsQ30,0.9900\n</NGS_Processing_Attributes>\n\n"
            f"<Code_Summary>\n{code}\n</Code_Summary>\n", encoding="utf-8")

    # ---- scans -----------------------------------------------------------------------------
    facts = {}
    for scan in SCANS:
        f = write_scan(d / f"{scan}.ome.tiff", scan, [r for r in rows if r["scan"] == scan], rng)
        facts.update({(scan, roi, seg): v for (roi, seg), v in f.items()})
    for r in rows:
        f = facts.get((r["scan"], r["roi"], r["segment"]))
        r["area"] = round(f["areaPx"] * PX_UM * PX_UM, 2) if f else 20000.0
        r["nuclei"] = int(r["area"] / 110)

    # ---- annotation files --------------------------------------------------------------------
    head = ["Sample_ID", "slide name", "scan name", "panel", "roi", "segment", "aoi", "area", "tissue", "genotype", "nuclei"]
    panel_txt = "(v1.0) Demo Panel, (v1.0) Demo Addon"
    table = [head, [ntc_id, "No Template Control"] + [""] * 9]
    for r in rows:
        table.append([r["sample_id"], SLIDE[r["scan"]], r["scan"], panel_txt, r["roi"], r["segment"],
                      f"{r['segment']}-aoi-001", r["area"], r["tissue"], r["genotype"], r["nuclei"]])
    write_xlsx(d / "annotations.xlsx", "Annotation", table)
    q = lambda v: '"' + str(v).replace('"', '""') + '"' if isinstance(v, str) and ("," in v or '"' in v) else str(v)
    (d / "annotations.csv").write_text("\n".join(",".join(q(v) for v in row) for row in table) + "\n", encoding="utf-8")
    pre = ["Experiment Summary", "Readout group name\tDEMO", "Readout mode\tNGS", f"Number Of AOIs\t{len(rows)}", "",
           "Library Prep Summary", f'Library Prep Plate\tCollection Plate\nDEMO-A\t"=""{PLATE}"""', "", "Annotations"]
    ws = ["\t".join(head[:8])] + [f"{ntc_id}\tNo Template Control"]
    ws += ["\t".join([r["sample_id"], SLIDE[r["scan"]], r["scan"], panel_txt, f'"=""{r["roi"]}"""', r["segment"],
                      f"{r['segment']}-aoi-001", f"{r['area']:.2f}"]) for r in rows]
    (d / "Demo_LabWorksheet.txt").write_text("\n".join(pre + ws) + "\n", encoding="utf-8")
    region = {"cortex": "outer", "tumour": "core", "stroma": "outer"}
    sp = [["SlideName", "ScanLabel", "ROILabel", "SegmentLabel", "SegmentDisplayName", "AOISurfaceArea", "AOINucleiCount", "Region"]]
    sp += [[SLIDE[r["scan"]], r["scan"], r["roi"], r["segment"], f"{r['scan']} | {r['roi']} | {r['segment']}", r["area"],
            r["nuclei"], region[r["tissue"]]] for r in rows]
    write_xlsx(d / "segment_properties.xlsx", "SegmentProperties", sp)

    # ---- truth: the collapse by plain loops, background, LOQ ---------------------------------
    names, used = [], set()
    for pname in sorted(panels, key=lambda p: p.replace(" ", "_")):     # the order prepare reads the .pkc files in
        for name, cc, k in panels[pname]:
            if cc.startswith("Negative"):
                continue
            nm = name if name not in used else f"{name} ({pname})"
            used.add(nm)
            names.append((nm, pname, name))

    def geo(v):
        if min(v) < 0.5:
            v = [x + 0.5 for x in v]
        lg = [math.log(x) for x in v]
        m = sum(lg) / len(lg)
        sd = math.sqrt(sum((x - m) ** 2 for x in lg) / (len(lg) - 1)) if len(lg) > 1 else 0.0
        return math.exp(m), math.exp(sd)
    values, bg = {}, {}
    for r in rows:
        c = counts[r["sample_id"]]
        col = []
        for nm, pname, name in names:
            v = [float(c[x]) for x in probes[(pname, name)]]
            col.append(v[0] if len(v) == 1 else geo(v)[0])
        values[r["sample_id"]] = col
        gm, gsd = geo([float(c[x]) for x in probes[("Demo Panel", "NegProbe-DEMO")]])
        bg[r["sample_id"]] = {"negGeoMean": gm, "loq": max(2.0, gm * gsd ** 2)}
    truth = {
        "note": "SYNTHETIC", "pixelSizeUm": PX_UM, "scans": {k: list(v) for k, v in SCANS.items()}, "gap": 512,
        "genes": [n for n, _, _ in names], "negativeTargets": ["NegProbe-DEMO", "Addon Negative"],
        "fakeId": FAKE_ID, "slideNames": sorted(set(SLIDE.values())), "ntc": ntc_id, "ntcTotal": 3,
        "aois": [{"sample_id": r["sample_id"], "scan": r["scan"], "roi": r["roi"], "segment": r["segment"],
                  "tissue": r["tissue"], "genotype": r["genotype"], "area": r["area"], "nuclei": r["nuclei"],
                  "values": values[r["sample_id"]], **bg[r["sample_id"]], **seq[r["sample_id"]],
                  "mask": facts.get((r["scan"], r["roi"], r["segment"]))} for r in rows],
    }
    (out / "truth_geomx.json").write_text(json.dumps(truth), encoding="utf-8")
    print(f"synthetic GeoMx study -> {d}  ({len(rows)} AOIs + 1 NTC, {len(names)} targets, 2 scans)")


if __name__ == "__main__":
    main()
