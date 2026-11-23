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

See it before you start: **[the example datasets](https://sidolans01.mgifive.org/SODA-examples/)**.
Pick a dataset from the **Dataset** menu, scroll to zoom, type a gene into **Search gene**, and
try **Tissue | Side by side | UMAP** at the top.

To publish your own data you do three things, once per collection of datasets:

```
1. PREPARE  (your computer)        2. UPLOAD  (your cloud storage)      3. YOUR SITE  (GitHub, free)
   the files your instrument  ───►    a folder online that you   ───►     your copy of SODA, which
   wrote become web-ready files       own and pay for                     opens that folder for anyone
```

This guide covers step 1. Steps 2 and 3 are in [SHARE_YOUR_DATA.md](SHARE_YOUR_DATA.md). Your
data never goes onto GitHub: your site only reads it from **your** storage while someone is
looking.

**What you need**
- A Windows or Mac computer with free disk space of about 3 times the size of your raw data.
- Your instrument's output folder(s): Visium / Visium HD (Space Ranger `outs` + the microscope
  image), Xenium, CosMx, MERSCOPE, Stereo-seq, GeoMx DSP, PhenoCycler / CODEX / other OME-TIFF
  and QPTIFF images, or single-cell RNA-seq (`.h5ad`). Exactly what each needs:
  [PREPARING_DATA.md](PREPARING_DATA.md).
- To publish: cloud storage that can serve files to a web page (Microsoft Azure, Amazon S3 or
  Google Cloud; Box, Google Drive, Dropbox and OneDrive do **not** work) and a free GitHub
  account. Both are explained in [SHARE_YOUR_DATA.md](SHARE_YOUR_DATA.md).
- About 1 to 2 hours the first time, most of it waiting.

**Two ways to prepare.** Both run the same code and make the same export:
- **The SODA app** (easiest, no code): [Part A](#part-a-the-soda-app-no-code), then Part 3.
- **The notebook** (more control: every setting, in a page you can read and change): Parts 1
  and 2, then Part 3.

**How much memory?** Most datasets need a few GB; the largest (Visium HD with 2 µm bins, a 5K
Xenium panel, a whole Stereo-seq chip) need 8 to 16 GB. The table:
[DESKTOP_APP.md, How much memory do I need?](DESKTOP_APP.md#how-much-memory-do-i-need).

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

## Before you publish: is this data allowed to be public?

Anyone with your link can view **and download** the exported data. Publish only data that may be
public: published or example data, with nothing that identifies a patient or donor. If not, see
[SHARE_YOUR_DATA.md, Before you start](SHARE_YOUR_DATA.md#before-you-start-who-may-see-the-data).

---

## Part A. The SODA app (no code)

### Step A1. Download and install

Open https://github.com/FIVE-MGI/SODA/releases/latest and download, **on the computer that will
run it**:
- **Windows 10 or 11:** `SODA-<version>-setup.exe`. Double-click it. The first time, Windows says
  "Windows protected your PC": click **More info**, then **Run anyway**. No administrator rights
  are needed. ✅ SODA is in the Start menu.
- **Mac with Apple Silicon (M1 or newer), macOS 14 or newer:** `SODA-<version>-macos-arm64.dmg`.
  Double-click it and drag SODA onto **Applications**. The first time, macOS refuses it once
  ("damaged" or "could not verify"): follow **README.txt** in the .dmg (one Terminal line, or
  System Settings → Privacy & Security → **Open Anyway**). There is no Intel Mac version.

### Step A2. Put your data in one folder

Make one folder for your data (for example `D:\SODA_data`) and put **each dataset in its own
folder** inside it, as your instrument wrote it (any names). Make a second, empty folder for the
output on a disk with room (about 3 times the data).

### Step A3. Prepare

Start SODA (Start menu, or Applications). After a few seconds its page opens in your browser at
`http://127.0.0.1:8770/app/` (a black window opens too on Windows: leave it open, it is SODA).
Under **Prepare a folder of datasets**, pick the data folder and the output folder (**Browse...**
or paste the path) and click **Prepare**. Progress shows on the page; large datasets take a while.

✅ A table lists each dataset as ok, and a **green** banner says "clean: the export folder can
be shared". A dataset marked failed or skipped says why and what to add. A **red** "NOT clean"
banner means sample IDs were found: do not upload; see
[DEIDENTIFICATION.md](DEIDENTIFICATION.md).

### Step A4. Look at it

Click **Open viewer**: your datasets open in a new tab, **only on your computer**. ✅ every
dataset looks right. To stop SODA: close its black window (Windows) or press **Quit SODA** on its
page (Mac).

Then go to [Part 3](#part-3-put-it-online) with the `export` folder inside your
output folder. More about the app: [DESKTOP_APP.md](DESKTOP_APP.md).

---

## Part 1. Or the notebook: install (once, about 15 minutes)

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
soda-env\Scripts\python -m pip install jupyter
```

**Mac (Terminal):**
```bash
cd ~/soda
```
```bash
bash install.sh
```
```bash
soda-env/bin/python -m pip install jupyter
```

✅ The second line ends with **`SODA is installed`**; the third ends with `Successfully installed`.
Everything goes into the `soda-env` folder inside `soda`; nothing else on your computer changes.

---

## Part 2. Prepare your data (on your computer)

### Step 4. Open the SODA notebook

Pick **one** of these two ways.

**A. In your browser (JupyterLab).** In the same terminal:

- **Windows:** `soda-env\Scripts\jupyter lab`
- **Mac:** `soda-env/bin/jupyter lab`

A browser tab called **JupyterLab** opens. (Keep the terminal open; closing it stops the
notebook.) On the left, double-click **notebooks** → **SODA_pipeline.ipynb**.

**B. In VS Code** (if you use it):

1. **File → Open Folder…** → choose the `soda` folder.
2. On the left, open **notebooks** → **SODA_pipeline.ipynb**. If VS Code offers to install
   the **Python** and **Jupyter** extensions, say yes.
3. Top right of the notebook: **Select Kernel** → **Python Environments** → the one marked
   **`soda-env`** (or **Jupyter Kernel** → **SODA**). (Not listed? Close the notebook and open it again.)

Either way: the text cells explain each step; the grey boxes with code are the ones you run.
Commands that start with `powershell`, `bash` or `soda-env` (like Step 3's) go in a **terminal**,
never in the notebook. In VS Code, **Terminal → New Terminal** opens one already in the `soda`
folder.

Run the first code cell (under **0. Setup**): click in it, press **Shift + Enter**.
✅ No line says `MISSING` or `WRONG VERSION`.

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
   CLOUD_URL = ""              # leave empty for now (Part 3)
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

## Part 3. Put it online

Your `export` folder (inside the app's output folder, or `SODA_work/export` from the notebook)
is ready. Follow [SHARE_YOUR_DATA.md](SHARE_YOUR_DATA.md) from
[Step 2](SHARE_YOUR_DATA.md#step-2-make-your-cloud-storage-once):

- **Step 2-4:** make your cloud storage (public read), upload the `export` folder with an
  upload tool, and set the CORS rule that lets your site read it.
- **Step 5-7:** click **Use this template** on https://github.com/FIVE-MGI/SODA, turn on
  GitHub Pages, and paste your storage folder's address into `web/SpatialViewer/config.js`.

✅ Your site, `https://<your GitHub name>.github.io/<your copy>/SODA/`, lists your datasets and
genes load.

**Notebook users:** Section 7 can print the upload command and test your storage the way a
browser would. Put your storage folder's address in `CLOUD_URL` (Settings cell) and your site's
address in `SODA_PAGE` (Advanced settings), run both cells, then Section 7 and `test_link()`.

**Adding datasets later:** put the new folder with the others, prepare again (datasets already
done are skipped), and upload the `export` folder again (replace existing files). Your site
lists the new dataset by itself.

## If something goes wrong

| What you see | What to do |
|---|---|
| Step 3: "SODA needs 64-bit Python 3.12 or 3.13" | install Python 3.13 (Step 1), close the terminal, open a new one, run Step 3 again |
| Step 3: red text, no `SODA is installed` | run the same line again (it is safe); if it fails again, copy the last 20 lines into a [GitHub issue](https://github.com/FIVE-MGI/SODA/issues) |
| Step 4: `jupyter` not found | Step 3's third line did not finish: run it again |
| Step 7.2: NOT SURE / missing | follow what the table says; what each instrument needs: [PREPARING_DATA.md](PREPARING_DATA.md) |
| Step 7.3 or the app: a dataset failed with "ran out of memory" | prepare one dataset at a time, or use a computer with more memory ([how much](DESKTOP_APP.md#how-much-memory-do-i-need)) |
| Step 7.5: names listed | rename and redo, as Step 7.5 says |
| Your site shows a message, or no data | [SHARE_YOUR_DATA.md, If something goes wrong](SHARE_YOUR_DATA.md#if-something-goes-wrong) |

More: [SHARE_YOUR_DATA.md](SHARE_YOUR_DATA.md) and the [README](../README.md#troubleshooting).
