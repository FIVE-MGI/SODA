# De-identification

Anything you upload can end up being seen by people you did not plan for. Sample or
patient identifiers must not be in it: not in dataset ids, labels, folder names, column
names, category values or `config.js`.

`web/SpatialViewer/deid.py` helps with two things:

1. **Replacing** sample IDs with plain numbers (`Sample 1`, `Sample 2`, ...), the same
   number for the same ID every time, using a private key file.
2. **Checking** files and folders for anything that still looks like a sample ID, as a
   gate before upload.

It is a pattern matcher, not a guarantee. Read [Limits](#limits) before relying on it.

## What it recognises

An ID, for `deid.py`, is one of the prefixes `INT`, `TFI`, `UTD`, `IO`, `NM`, `T` or `I`
(any case), optionally followed by `-`, `_` or a space, then 2 to 4 digits, optionally
one trailing letter. It must not be glued to a letter or digit before it, and no digit may
follow. Shorthand lists after an ID are expanded too: an ID followed by `_` and four
more digits counts as two IDs, and an ID followed by `-` and one or two digits replaces
its last digits (a range written in short form).

These prefixes are the identifier formats of the lab that wrote the tool. **Your
identifiers probably look different.** If so, the pattern does nothing for them until you
change it: edit `_ID` near the top of `deid.py` to match your formats, and test it with
`deid.py text` (below) before trusting it.

Each ID is replaced by its number in one of three forms:

| Form | Used for | Result |
|---|---|---|
| `text` | labels, free text | `Sample 3` (several in a row: `Samples 3 / 4`) |
| `value` | table cells | `3` |
| `slug` | ids and folder names | `s3` |

A letter suffix is kept: an ID ending in `b` becomes `3b`.

## Where it runs automatically

- **`prepare.py rnaseq`** passes the dataset id, label, every `obs` column name and every
  category value through it. A column that *is* the sample identity (named `sample`,
  `sample_id`, `donor`, `patient`, `subject`, `individual`, `orig.ident`, `library`,
  `batch`, `specimen` and similar, or whose values mostly contain IDs) has **every** value
  replaced by a number, whatever its format. The manifest is checked again before it is
  written; the export stops rather than write a leftover ID. Cell barcodes are never
  written.
- **`prepare.py export`** refuses to write the highly variable gene list or the UMAP
  index if they contain an ID or a Space Ranger barcode, and drops the local file paths
  from the manifest.

It does **not** rewrite the labels, channel names or overlay names of image datasets.
For those, `deid.py check` is the gate.

## The key file

The mapping from original IDs to numbers is kept in one CSV file with the columns
`number, original_id, first_seen, context`. New IDs are appended the first time they are
seen, so numbering stays consistent across datasets and runs.

**This file re-identifies your data.** Treat it like the source data:

- keep it **outside** the repository and outside every export folder;
- never upload it, commit it or paste it into a ticket;
- back it up privately: if it is lost, the next run numbers everything afresh and the
  numbers will no longer match what you published before.

Its location comes from the `SPATIALVIZ_KEY` environment variable. **If that is not set,
the file is `sample_id_key.csv` in whatever folder you run the command from**, which is
easy to commit by accident. Set it before running `prepare.py` or `deid.py`.

macOS / Linux (add the line to `~/.bashrc` or `~/.zshrc` to keep it):

```bash
export SPATIALVIZ_KEY="$HOME/private/spatial-viewer/sample_id_key.csv"
```

Windows PowerShell (current window only):

```powershell
$env:SPATIALVIZ_KEY = "$HOME\private\spatial-viewer\sample_id_key.csv"
```

Windows, permanently for new terminal windows:

```powershell
setx SPATIALVIZ_KEY "$HOME\private\spatial-viewer\sample_id_key.csv"
```

Use a folder only you (or your project's approved staff) can read.

## Command line

```bash
python web/SpatialViewer/deid.py text  "<some label>"     # replace IDs, print the result (also: value, slug)
python web/SpatialViewer/deid.py check <file or folder> [...]
python web/SpatialViewer/deid.py key                      # print the key: do not do this on a shared screen
```

`text`, `value` and `slug` add any new IDs to the key file.

## `deid.py check` as a gate

```bash
python web/SpatialViewer/deid.py check demo/export/demo_visiumhd web/SpatialViewer/config.js
```

It looks at:

- the **names** of every file and folder given, recursively;
- the **contents** of text files (`.js .json .html .txt .md .css .geojson .csv .py .ps1
  .yml .yaml`) under 50 MB, line by line.

It skips `.venv`, `__pycache__`, `.git`, `tiles` and `gene` folders (image tiles and
per-gene files: binary, named by numbers and gene symbols), and `*.link.txt` files
(encrypted, random-looking text). A line containing the marker `deid:ok` is skipped: it is
meant for documented, made-up examples in source code. Do not use it to silence a real hit.

It prints every hit as `name <path> -> <id>` or `content <path>:<line> -> <id>`, then
`[deid] N sample ID(s) found` and exits with status 1, or prints `[deid] clean` and
exits with 0. So it can guard an upload:

macOS / Linux:

```bash
python web/SpatialViewer/deid.py check demo/export/demo_visiumhd && azcopy copy "demo/export/demo_visiumhd" "https://<account>.blob.core.windows.net/<container>?<UPLOAD-SAS>" --recursive
```

Windows PowerShell:

```powershell
python web/SpatialViewer/deid.py check demo/export/demo_visiumhd
if ($LASTEXITCODE -eq 0) { azcopy copy "demo/export/demo_visiumhd" "https://<account>.blob.core.windows.net/<container>?<UPLOAD-SAS>" --recursive }
```

## Limits

- **It only knows its pattern.** Names, initials, dates of birth, medical record or
  accession numbers, study codes in other formats, and free-text notes are not detected.
- **It does not look at images.** A slide label, barcode sticker or handwriting that is
  visible in the scanned image is published with the tiles. Crop it out with `--crop`
  when you prepare the dataset, and look at the whole slide at low zoom before uploading.
- **It does not read binary files** (`.bin`, tiles). Those hold numbers derived from the
  data, not text, but a de-identification review should know that.
- **Image datasets are only checked, not rewritten.** Labels, channel names and overlay /
  region names come from your files and `--label`; choose them carefully.
- **False positives happen.** A gene, marker or file name can match the pattern. Read what
  it reports instead of renaming blindly.
- It is a helper, not a certification. Whether a dataset may be shared, and how, is decided
  by your IRB, data-use agreements and institutional policy.

**Always also review by eye**: the dataset labels and groups in `config.js`, the `label`
and channel names in each `manifest.json`, overlay and region names, RNA-seq category
values, `about.html`, and the URL of the site itself.
