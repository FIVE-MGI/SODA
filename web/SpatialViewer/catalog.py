"""Write datasets.json: the list of every dataset in an export folder.

The SODA page can open a whole export folder from one link,

    https://<soda page>/index.html?data=<link to the export folder>

and list every dataset in it, with no config.js to edit. To do that it reads one
small file at the top of the folder, datasets.json, which this script writes:

    python catalog.py <export folder> [--title "Our nerve atlas"]

Run it after the last `prepare.py export` and before uploading. Running it again is
safe: labels, groups and order that you changed by hand in datasets.json are kept,
datasets that were added are appended, and ones whose folder is gone are dropped
(each change is printed).

datasets.json (everything marked optional is yours to add or change by hand, in any
text editor; this script keeps it):

    {"soda": 1,
     "title": "Our nerve atlas",                    # optional: the page title
     "byline": "Smith lab, 2026",                   # optional: the line under the title
     "paper": "https://doi.org/...",                # optional: the "Link to paper" button
     "default": "demo_visiumhd",                    # optional: the dataset opened first
     "datasets": [
       {"id": "demo_visiumhd",                      # the folder name; do not change
        "label": "Demo section",                    # what the list shows
        "group": "Spatial transcriptomics",         # the heading it sits under
        "kind": "visium",                           # filled in by this script
        "hideChannels": ["CD11c test"],             # optional: stains to leave out of the list
        "defaultChannels": ["DAPI", "CD3"],         # optional: stains switched on at first
        "labelPositions": {"T cells": [1.2, -3.4]}  # optional: hand-placed UMAP cluster labels
       }, ...]}

Group order on the page is the order groups first appear in this list.

Only the standard library is used, so this runs on any Python 3.8+.
"""
import argparse
import json
import os
import sys
import time

CATALOG = "datasets.json"
FORMAT_VERSION = 1

# What a reader would call each kind of dataset; the picker groups by this.
GROUP_BY_KIND = {
    "visium": "Spatial transcriptomics",
    "xenium": "Spatial transcriptomics",
    "cosmx": "Spatial transcriptomics",
    "merscope": "Spatial transcriptomics",
    "stereoseq": "Spatial transcriptomics",
    "geomx": "Spatial transcriptomics",
    "rnaseq": "Single-cell RNA-seq",
    "tiff": "Images",
    "akoya": "Images",
}

# Files that must never be uploaded. Finding one inside the export folder stops the
# script: it means something private was copied next to the data by mistake.
NEVER_UPLOAD = {"sample_id_key.csv"}


def group_for(manifest):
    kind = manifest.get("kind", "")
    if kind in ("tiff", "akoya") and any(e.get("kind") == "protein" for e in manifest.get("expr") or []):
        return "Protein images"
    return GROUP_BY_KIND.get(kind, "Other")


def read_manifest(folder):
    path = os.path.join(folder, "manifest.json")
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def count_files(folder):
    n = size = 0
    for root, _dirs, files in os.walk(folder):
        for name in files:
            n += 1
            try:
                size += os.path.getsize(os.path.join(root, name))
            except OSError:
                pass
    return n, size


def find_private_files(root):
    hits = []
    for dirpath, _dirs, files in os.walk(root):
        for name in files:
            if name.lower() in NEVER_UPLOAD:
                hits.append(os.path.join(dirpath, name))
    return hits


def looks_like_cache(folder):
    """A prepare cache (raw full-resolution zarr) instead of an export."""
    try:
        return any(n.endswith(".zarr") for n in os.listdir(folder))
    except OSError:
        return False


def build(root, title=None, count=False, fresh=False, log=print):
    root = os.path.abspath(root)
    if not os.path.isdir(root):
        raise SystemExit(f"[catalog] not a folder: {root}")

    private = find_private_files(root)
    if private:
        raise SystemExit("[catalog] *** private file inside the export folder, remove it before uploading:\n  "
                         + "\n  ".join(private))

    found = {}
    for name in sorted(os.listdir(root)):
        folder = os.path.join(root, name)
        if not os.path.isdir(folder) or not os.path.isfile(os.path.join(folder, "manifest.json")):
            continue
        if looks_like_cache(folder):
            raise SystemExit(f"[catalog] *** {name} looks like a prepare cache (it holds .zarr folders), "
                             "not an export. Point this at the folder `prepare.py export --dest` wrote.")
        try:
            m = read_manifest(folder)
        except (OSError, ValueError) as e:
            log(f"[catalog] skipped {name}: manifest.json unreadable ({e})")
            continue
        if m.get("id") and m["id"] != name:
            log(f"[catalog] note: folder {name} holds a manifest with id {m['id']!r}; "
                "the folder name is what the page uses")
        found[name] = m

    if not found:
        raise SystemExit(f"[catalog] no dataset folders (with a manifest.json) in {root}")

    path = os.path.join(root, CATALOG)
    old = {"datasets": []}
    if os.path.isfile(path) and not fresh:
        try:
            # utf-8-sig: Notepad and PowerShell on Windows save JSON with a byte-order mark
            with open(path, encoding="utf-8-sig") as f:
                old = json.load(f)
        except (OSError, ValueError) as e:
            # never overwrite it: it may hold labels someone typed by hand
            raise SystemExit(f"[catalog] *** {path} could not be read ({e}).\n"
                             "  Fix the typo it points to, or run again with --fresh to replace it "
                             "(labels and groups typed into it are then lost).")

    entries = []
    kept = set()
    # first the datasets already listed, in the order (and with the labels) the host left them
    for e in old.get("datasets", []):
        did = e.get("id")
        if did in found and did not in kept:
            e = dict(e)
            e["kind"] = found[did].get("kind", e.get("kind", ""))
            entries.append(e)
            kept.add(did)
        elif did not in found:
            log(f"[catalog] dropped {did}: its folder is gone")
    # then any new ones
    for did, m in found.items():
        if did in kept:
            continue
        entries.append({"id": did, "label": m.get("label") or did, "group": group_for(m), "kind": m.get("kind", "")})
        log(f"[catalog] added {did} ({entries[-1]['group']})")

    if count:
        total_n = total_b = 0
        for e in entries:
            n, b = count_files(os.path.join(root, e["id"]))
            e["files"], e["bytes"] = n, b
            total_n += n
            total_b += b
        log(f"[catalog] {total_n:,} files, {total_b / 1e9:.2f} GB to upload")
    else:
        for e in entries:              # counts from an earlier --count run would be stale now
            e.pop("files", None)
            e.pop("bytes", None)

    out = {"soda": FORMAT_VERSION}
    t = title if title is not None else old.get("title")
    if t:
        out["title"] = t
    # anything else typed in at the top (byline, paper, default, ...) stays as it was
    for k, v in old.items():
        if k not in ("soda", "title", "datasets", "updated"):
            out[k] = v
    if out.get("default") and out["default"] not in {e["id"] for e in entries}:
        log(f"[catalog] note: \"default\" names {out['default']!r}, which is not in this folder; "
            "the page will open the first dataset instead")
    out["datasets"] = entries
    out["updated"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
        f.write("\n")
    log(f"[catalog] wrote {path} ({len(entries)} dataset(s))")
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("root", help="the export folder (the one holding one folder per dataset)")
    ap.add_argument("--title", help="a name for the whole collection, shown as the page title")
    ap.add_argument("--count", action="store_true", help="also count files and bytes per dataset (slow on network drives)")
    ap.add_argument("--fresh", action="store_true", help="ignore an existing datasets.json instead of keeping its edits")
    a = ap.parse_args(argv)
    build(a.root, a.title, a.count, a.fresh)


if __name__ == "__main__":
    sys.exit(main())
