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

**It only knows our lab's ID formats.** Before you use it, edit the `_ID` pattern near the
top of `web/SpatialViewer/deid.py` to match your sample IDs, and test it with `deid.py text`
(see [What it recognises](#what-it-recognises)). Until then it does not find your IDs.

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
  from the manifest. **Barcodes and cell ids are written by default** (`expr/<key>/ids.txt.gz`,
  as Loupe Browser and Xenium Explorer keep them); `--no-ids` (the notebook's
  `INCLUDE_IDS = False`, the desktop app's "Include barcode / cell id lists" box) leaves them
  out. Before a list is written, no id may match the sample-ID patterns and, for the platforms
  with a known format (Visium / Visium HD barcodes, Xenium, CosMx, MERSCOPE and Stereo-seq cell
  and bin ids, see `ids.txt.gz` in [DATA_FORMAT.md](DATA_FORMAT.md)), every id must fully match
  the instrument's id format. AnnData and SpatialData ids have no format check, only the
  sample-ID one. By default a failing layer gets no list (a `***` line in the export log) and
  the export goes on; with `--ids` (or `prepare.py ids`) one failing id stops the dataset and no
  list is written for any layer. Either way the message gives the layer and the count but never
  the ids, so a log holds none. Ids are never stripped or rewritten to make them pass. The
  checks know only the patterns in `deid.py`: you are responsible for what you publish.

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
python web/SpatialViewer/deid.py check export/my_dataset web/SpatialViewer/config.js
```

It looks at:

- the **names** of every file and folder given, recursively;
- the **contents** of text files (`.js .json .html .txt .md .css .geojson .csv .py .ps1
  .yml .yaml`) under 50 MB, line by line;
- the **contents** of gzip files under 50 MB that unpack to text (the `ids.txt.gz` barcode
  lists), line by line as a stream. A hit there is printed as
  `content <path>:<line> -> looks like a sample ID (not shown)`, never the ID (10 lines per
  file, then a count).

It skips `soda-env`, `.venv`, `__pycache__`, `.git`, `tiles` and `gene` folders (image tiles and
per-gene files: binary, named by numbers and gene symbols). A line containing the marker `deid:ok` is skipped: it is
meant for documented, made-up examples in source code. Do not use it to silence a real hit.

Gene packs (`genepack/*.pack`, the per-gene data of an export, see
[DATA_FORMAT.md](DATA_FORMAT.md#gene-packs-genepackhash-ipack)) are binary and are not read as
text. Instead, for every `manifest.json` it finds, the check reads each packed group's listing
(`genes.json`, `molecules/index.json`): every listed name is checked like any other text
(`gene <path> -> <id>`), a listing entry may carry only the gene fields (name, totals, shard,
offset, length), and every byte of every shard the manifest names must belong to a listed gene,
with no gap and no overlap (`pack <path> -> <problem>`). So a pack cannot carry anything a gene
listing does not account for. (Shards an earlier export left in the same folder are not named by
the manifest and never read; the export log lists them for deleting before the upload.)

It prints every hit as `name <path> -> <id>` or `content <path>:<line> -> <id>`, then
`[deid] N sample ID(s) found` and exits with status 1, or prints `[deid] clean` and
exits with 0. So it can guard an upload:

macOS / Linux:

```bash
python web/SpatialViewer/deid.py check export/my_dataset && azcopy copy "export/my_dataset" "https://<account>.blob.core.windows.net/<container>?<UPLOAD-SAS>" --recursive
```

Windows PowerShell:

```powershell
python web/SpatialViewer/deid.py check export/my_dataset
if ($LASTEXITCODE -eq 0) { azcopy copy "export/my_dataset" "https://<account>.blob.core.windows.net/<container>?<UPLOAD-SAS>" --recursive }
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
  One family is built in as an exception: the mouse MHC genes `H2-T3` ... `H2-T24`. A `T` plus
  digits directly after `H2-` is read as the gene symbol, neither replaced nor reported
  (`_GENE_PREFIX` in `deid.py`); without it the gene list of any mouse whole-transcriptome
  dataset fails the check.
- It is a helper, not a certification. Whether a dataset may be shared, and how, is decided
  by your IRB, data-use agreements and institutional policy.

**Always also review by eye**: the dataset labels and groups in `config.js`, the `label`
and channel names in each `manifest.json`, overlay and region names, RNA-seq category
values, `about.html`, and the URL of the site itself.
