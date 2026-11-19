# Getting started: publish your spatial omics data with SODA

For people who have never used the command line, Python or cloud storage. Follow the steps in
order. Each one says what to click or type and what you should see; if you see something else,
stop there and look at [If something goes wrong](#if-something-goes-wrong).

---

## What SODA does

SODA shows spatial omics data (a tissue image with genes or proteins measured across it) in a
web browser. Readers zoom into the tissue like a map, colour it by any gene or protein, and see
which cells or spots are alike (the UMAP and clusters). **Readers install nothing**: they open a
link.

See it before you start: **[the example datasets](https://soda.mgifive.org/index.html?data=https://fivestore.blob.core.windows.net/soda/examples/)**.
Pick a dataset from the **Dataset** menu, scroll to zoom, type a gene into **Search gene**, and
try **Tissue | Side by side | UMAP** at the top.

To publish your own data you do three things, once per collection of datasets:

```
1. PREPARE  (your computer)        2. UPLOAD  (your cloud storage)      3. SHARE  (one link)
   the files your instrument  ───►    a folder online that you   ───►     soda.mgifive.org/?data=<your folder>
   wrote become web-ready files       own and pay for                     opens your data for anyone
```

Your data never goes onto SODA's website: the page only reads it from **your** storage while
someone is looking.

**What you need**
- A Windows or Mac computer with free disk space of about 3 times the size of your raw data.
- Your instrument's output folder(s): Visium / Visium HD (Space Ranger `outs` + the microscope
  image), Xenium, CosMx, MERSCOPE, Stereo-seq, GeoMx DSP, PhenoCycler / CODEX / other OME-TIFF
  and QPTIFF images, or single-cell RNA-seq (`.h5ad`). Exactly what each needs:
  [PREPARING_DATA.md](PREPARING_DATA.md).
- Cloud storage that can serve files to a web page. This guide uses **Microsoft Azure** (Step 8);
  Amazon S3 and Google Cloud work too ([SHARE_YOUR_DATA.md, Step 3](SHARE_YOUR_DATA.md#step-3-upload-the-whole-export-folder)).
  Box, Google Drive, Dropbox and OneDrive do **not** work.
- About 1 to 2 hours the first time, most of it waiting.

**No data of your own yet?** Practise with a small public dataset: the MERSCOPE human fetal
colon sample (275 MB, Zenodo record 19450420, CC BY 4.0). Download link and command:
[REPRODUCE.md, section 9](../REPRODUCE.md#9-the-public-datasets) (`merscope_2`). Unzipping it
gives a folder `AAN4`: in Step 5, put that folder in `inputs/` as your dataset.

## Words you will meet

| Word | Meaning |
|---|---|
| **Terminal** | a window where you type commands. Windows: Start menu → type `PowerShell` → open it. Mac: Spotlight (⌘ Space) → type `Terminal` → open it. Paste with right-click (Windows) or ⌘ V (Mac); press Enter to run. |
| **Repository** | SODA's folder of code, kept on GitHub. |
| **Notebook** | a page in your browser made of boxes ("cells"). Each box runs one step; its answer appears under it. |
| **Run a cell** | click inside the box and press **Shift + Enter**. `[*]` beside it means "working"; a number means "done". |
| **Export** | the web-ready copy of your data that SODA writes. Only this is uploaded. |
| **Container** | a folder in Azure storage. |
| **CORS** | a setting on your storage that lets the SODA page read your files. Set once. |

## Before you publish: is this data allowed to be public?

Anyone with your link can view **and download** the exported data. Publish only data that may be
public: published or example data, with nothing that identifies a patient or donor. If not, see
[If your data must stay private](#if-your-data-must-stay-private) instead of Part 3.

---

## Part 1. Install (once, about 15 minutes)

### Step 1. Install Python

SODA needs **Python 3.13** (3.12 also works).

- **Windows:** go to https://www.python.org/downloads/windows/, download the latest
  **Python 3.13** "Windows installer (64-bit)", run it, tick **Add python.exe to PATH** on the
  first screen, then **Install Now**.
- **Mac:** go to https://www.python.org/downloads/macos/, download the latest **Python 3.13**
  "macOS 64-bit universal2 installer", open it and click through.

### Step 2. Get SODA

On https://github.com/FIVE-MGI/SODA click the green **Code** button → **Download ZIP**. Unzip it
to a short path, and rename the folder to `soda`:
- Windows: `C:\soda`
- Mac: your home folder, `~/soda`

(If you already use Git: `git clone https://github.com/FIVE-MGI/SODA.git soda` does the same.)

### Step 3. Install SODA's packages

Open a terminal and paste these lines one at a time (Enter after each).

**Windows (PowerShell):**
```powershell
cd C:\soda
```
```powershell
powershell -ExecutionPolicy Bypass -File install.ps1
```
```powershell
.venv\Scripts\python -m pip install jupyter
```

**Mac (Terminal):**
```bash
cd ~/soda
```
```bash
bash install.sh
```
```bash
.venv/bin/python -m pip install jupyter
```

✅ The second line ends with **`SODA is installed`**; the third ends with `Successfully installed`.
Everything goes into the `.venv` folder inside `soda`; nothing else on your computer changes.

---

## Part 2. Prepare your data (on your computer)

### Step 4. Open the SODA notebook

In the same terminal:

- **Windows:** `.venv\Scripts\jupyter lab`
- **Mac:** `.venv/bin/jupyter lab`

A browser tab called **JupyterLab** opens. (Keep the terminal open; closing it stops the
notebook.) On the left, double-click **notebooks** → **SODA_pipeline.ipynb**.

Run the first code cell (under **0. Setup**). ✅ No line says `MISSING` or `WRONG VERSION`.

### Step 5. Make a work folder and put your data in it

Make a folder on a disk with plenty of space, e.g. `D:\SODA_work` (Windows) or
`~/SODA_work` (Mac). Inside it make a folder **`inputs`**, and inside that **one folder per
dataset**:

```
SODA_work/
  inputs/
    mouse_brain_hd/        ← a Space Ranger outs/ folder + the full-resolution microscope image
    tumour_xenium/         ← a Xenium output folder, exactly as the instrument wrote it
```

Copy (don't move) your instrument's files into them. **The folder name is shown on the website**:
use a plain name such as `mouse_brain_hd`, never a patient ID, a sample number or a date.

### Step 6. Settings

1. In the cell that starts with `WORK =` (under **Settings: edit this cell**), put your work folder:
   ```python
   WORK = r"D:\SODA_work"      # Mac: WORK = "/Users/<you>/SODA_work"
   CLOUD_URL = ""              # leave empty for now; Step 9 fills it in
   ```
   Run it.
2. Run the **Advanced settings** cell as it is. (Optional: `TITLE = "My lab's atlas"` names
   the collection on the website.)

### Step 7. Check, prepare, export, look

Run these sections in order; each says when it is done.

1. **1. Put your files here**: ✅ lists your dataset folders.
2. **2. Check the inputs**: one table per dataset. ✅ **certain**.
   - **NOT SURE: pick one**: the table prints a line to copy into the `OPTIONS` cell; paste it, run again.
   - **missing**: a file the instrument writes is not in the folder; the table names it.
   - To give a dataset a nicer name, add to `OPTIONS`: `"mouse_brain_hd": {"label": "Mouse brain, Visium HD"},`
3. **3. Prepare**: the slow step (minutes to an hour per dataset). ✅ no `FAILED`.
4. **4. Export**: ✅ `SODA_work/export/` now holds one folder per dataset and a `datasets.json`.
5. **5. De-identification check**: scans the export for anything that looks like a sample ID.
   ✅ ends with `[deid] clean` (or says it was skipped: then read your names yourself). If it lists
   names, rename the folder or label, delete that dataset's folder in `export/` and `cache/`,
   and run 3 to 5 again.
6. **6. Look before you upload**: run `sizes()`, then `preview()`, and click the **first link**.
   Your data opens in your browser, **only on your computer**. ✅ every dataset looks right.
   Then run `stop_preview()` (delete the `#` at its start first).

---

## Part 3. Upload and share (public data)

### Step 8. Make your storage on Azure (once)

You need an Azure account (many universities provide one; otherwise https://azure.microsoft.com/free).
Azure charges a little for storage and for what readers download; small collections cost very little.

1. Go to https://portal.azure.com and sign in.
2. **Create the storage account.** Search bar → **Storage accounts** → **+ Create**:
   - Resource group: **Create new** → `soda`
   - Storage account name: lowercase letters and digits, unique worldwide, e.g. `smithlabsoda`
   - Region: near you. Performance: **Standard**. Redundancy: **LRS**.
   - **Advanced** tab: tick **Allow enabling anonymous access on individual containers**.
   - **Review + create** → **Create**. Wait for "Your deployment is complete" → **Go to resource**.
3. **Create the container.** Left menu **Data storage → Containers** → **+ Container**:
   - Name: `soda`
   - Anonymous access level: **Blob (anonymous read access for blobs only)**
   - **Create**.
4. **Let the SODA page read it (CORS).** Left menu **Settings → Resource sharing (CORS)** →
   **Blob service** tab → fill one row → **Save**:

   | Allowed origins | Allowed methods | Allowed headers | Exposed headers | Max age |
   |---|---|---|---|---|
   | `https://soda.mgifive.org` | `GET`, `HEAD`, `OPTIONS` | `*` | `*` | `3600` |

### Step 9. Upload the export folder

The easiest way is **Azure Storage Explorer**, a free drag-and-drop app:

1. Download it from https://azure.microsoft.com/products/storage/storage-explorer/, install,
   open, and **Sign in with Azure** (the person icon on the left).
2. In the left tree: your subscription → **Storage Accounts** → your account → **Blob Containers** → `soda`.
3. Click **Upload → Upload Folder…** → choose your `SODA_work/export` folder → **Upload**.
   It appears in the container as a folder named `export`.
   ✅ The activity bar at the bottom says the transfer completed, with no failures.

Upload **only** `export`. Never `cache`, `inputs` or any key file.

Your folder's address is now:
```
https://<your storage account name>.blob.core.windows.net/soda/export
```

### Step 10. Get your link and test it

1. In the notebook's **Settings** cell, put that address in `CLOUD_URL` and run the cell:
   ```python
   CLOUD_URL = "https://smithlabsoda.blob.core.windows.net/soda/export"
   ```
2. Under **7. Upload and share**, set `IS_THIS_DATA_PUBLIC = "yes"` and run the cell.
   It prints **THE LINK TO SHARE**. (It also prints an upload command: you already uploaded in
   Step 9, so skip it.)
3. Run **`test_link()`** below it. ✅ **OK**. If it reports a CORS problem, check Step 8.4 (the
   origin must be exactly `https://soda.mgifive.org`), wait a minute, and run it again.
4. Open the link in a **private browser window**. ✅ your datasets are in the **Dataset** menu
   and genes load.

**Done.** That link is what you put in your paper, your lab page or an email. It looks like:
```
https://soda.mgifive.org/?data=https%3A%2F%2Fsmithlabsoda.blob.core.windows.net%2Fsoda%2Fexport%2F
```

### Changing names later

Edit `datasets.json` in `SODA_work/export` with any text editor (names, groups, order, title),
then upload just that one file again (Storage Explorer: open the `export` folder → **Upload →
Upload Files**, overwrite). Reload the page. Every setting:
[SHARE_YOUR_DATA.md, What you can change in datasets.json](SHARE_YOUR_DATA.md#what-you-can-change-in-datasetsjson).

### Adding more datasets later

Put the new folder in `SODA_work/inputs`, run the notebook's Sections 1 to 6 again (datasets
already done are skipped), and upload the `export` folder again (Storage Explorer asks about
files that already exist: **Replace**). The same link now lists the new dataset too.

---

## Optional: your own SODA website

Your own copy of the SODA page at your own address, with your title and logos, opening your data
without the `?data=` part: [SHARE_YOUR_DATA.md, 2. Make it your own](SHARE_YOUR_DATA.md#2-make-it-your-own)
(about 30 minutes in the browser, with a free GitHub account). Add your site's address to the CORS
row from Step 8.4.

## If your data must stay private

Do Parts 1 and 2 the same way. Then, instead of Part 3:
- **A few people, for a while:** keep the container **Private**, and share a link carrying a
  read-only pass (Azure SAS, set as `READ_PASS` in the notebook's advanced settings). Anyone with
  that link can read the data until the pass expires: send it like a password, never put it in a paper.
- **Your lab, long term:** a password-protected SODA site on Azure: [deploy/azure/README.md](../deploy/azure/README.md).
- **Only you:** upload nothing; Section 6's preview is enough.

Details: [SHARE_YOUR_DATA.md, Start here](SHARE_YOUR_DATA.md#start-here-who-may-see-the-data).
Sample IDs and the key file: [DEIDENTIFICATION.md](DEIDENTIFICATION.md).

## If something goes wrong

| What you see | What to do |
|---|---|
| Step 3: "SODA needs 64-bit Python 3.12 or 3.13" | install Python 3.13 (Step 1), close the terminal, open a new one, run Step 3 again |
| Step 3: red text, no `SODA is installed` | run the same line again (it is safe); if it fails again, copy the last 20 lines into a [GitHub issue](https://github.com/FIVE-MGI/SODA/issues) |
| Step 4: `jupyter` not found | Step 3's third line did not finish: run it again |
| Step 7.2: NOT SURE / missing | follow what the table says; what each instrument needs: [PREPARING_DATA.md](PREPARING_DATA.md) |
| Step 7.3: the computer runs out of memory | prepare one dataset at a time; [ADVANCED.md](ADVANCED.md) |
| Step 7.5: names listed | rename and redo, as Step 7.5 says |
| The link says "could not reach" your folder | CORS (Step 8.4), or the address in `CLOUD_URL` is not the folder that holds `datasets.json` |
| The link says the folder is private / 403 | the container's access level must be **Blob** (Step 8.3), and the account must allow it (Step 8.2, Advanced tab) |
| A dataset is missing from the list | its folder was not uploaded, or `datasets.json` was uploaded before Section 4 finished |

More: [HOSTING.md](HOSTING.md), [SHARE_YOUR_DATA.md](SHARE_YOUR_DATA.md), and the full README.
