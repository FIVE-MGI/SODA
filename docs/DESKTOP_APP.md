# The SODA desktop app (Windows and Mac)

The goal: someone with no Python, no terminal and no notebook downloads one file, double-clicks
it, picks a folder of datasets, and gets an export they can look at in the viewer and, once the
sample ID check is clean, upload. It is built and released the way FIVE Atlas is, so the lab has
one way of shipping desktop apps.

Status: the Windows build is made and tested on a lab PC (numbers below). The macOS build script
and both GitHub Actions workflows are written but have not run yet: they need a Mac or a GitHub
runner (see [What is still open](#what-is-still-open)).

## How FIVE Atlas is packaged (the pattern to follow)

FIVE Atlas is **not** an Avalonia/.NET app, though it is sometimes described as one. A chat on
2026-08-05 ("Five Atlas on Avalonia for macOS") looked at Avalonia and decided against it: it
would mean rewriting a Python + web app in C#, and the scientific Python stack has no .NET
equivalent. FIVE Atlas is:

| | FIVE Atlas |
|---|---|
| What it is | a Python (FastAPI) server plus a built React page, frozen with **PyInstaller 6.21.0** as a **one-folder** build (not one-file: one-file re-extracts hundreds of MB on every launch) |
| UI | no window of its own: the app opens the **system browser** at `http://127.0.0.1:8050` |
| Windows | `FiveAtlas.exe` with a console window (it shows the address and any error; closing it stops the app). Shipped as `FiveAtlas-<ver>-win64.zip`; `installer.iss` makes a per-user Inno Setup installer when Inno Setup is installed (it is not on the build machines, so in practice: the zip). **Not code-signed**: SmartScreen says "Windows protected your PC" -> More info -> Run anyway |
| macOS | `FiveAtlas.app` (no console, `LSUIElement`: no Dock icon, a log file, an idle watchdog) in `FiveAtlas-<ver>-macos-arm64.dmg` and `.zip`. **Apple Silicon only** (the Intel runner pool never started). **Ad-hoc signed only, not notarised**: the README in the .dmg walks people through `xattr -dr com.apple.quarantine` or System Settings -> Open Anyway |
| Build | `build_release.ps1` / `build_release_mac.sh` -> one `FiveAtlas.spec` for both. The version comes from the git tag and is baked into the app at build time. A launch check runs the built binary before packaging |
| CI | `.github/workflows/build-windows.yml` (windows-latest) and `build-macos.yml` (macos-15, arm64): run by hand or by pushing a `v*` tag; actions pinned to commit SHAs; the built app is smoke-tested headless (health, page, `/api/selftest`, version), uploaded as a workflow artifact, and on a tag attached to **one GitHub Release** with both downloads and `SHA256SUMS` |
| Release notes | the tag's section of `CHANGELOG.md`, plus install instructions, on the Release page |
| Updates | the app checks `releases/latest` on GitHub (a few times a day), shows "a newer version is out"; from 1.3.3 it also downloads the zip, checks it against `SHA256SUMS`, and swaps itself in on restart (`selfupdate.py`) |
| Where people get it | `https://github.com/jessicahelms/FiveAtlas/releases/latest` (a public repository) |

## Options for SODA

SODA is a Python pipeline (`prepare.py`, `batch.py`, `detect.py`, `deid.py` and the modules they
import, pinned in `requirements-lock.txt`) plus a static page (`web/SpatialViewer/index.html`).
Three ways to put that in an app:

**(a) A small Avalonia/.NET shell** with an embedded Python, the locked packages and the page;
"Prepare a folder" runs `prepare.py batch`, "Open viewer" starts a local server and opens the
browser or an embedded web view.
- For: a real native window with buttons; an embedded web view would keep the viewer inside the app.
- Against: two runtimes to ship and keep in step (.NET 10 *and* Python 3.13, about 70 MB more);
  a second language in a Python + JavaScript repository; an embedded Python on macOS must still
  be relocatable and signed binary by binary (the same problem PyInstaller already solves);
  Avalonia's web view is young, and the viewer needs WebGL2 and workers, which the system
  browser already does well. It is not "the same way as FIVE Atlas": FIVE Atlas tried this
  idea and rejected it. It also needs Avalonia and its NuGet packages, which are not on this
  machine (the local NuGet cache has no Avalonia package), so it could not even be built here
  without downloads.

**(b) PyInstaller one-folder build of a Python launcher** (what FIVE Atlas does).
- For: the identical toolchain, version (6.21.0), workflow layout, release page, install
  instructions and Gatekeeper story as FIVE Atlas; one language; the launcher is a few hundred
  lines of standard library; SODA's scripts run unchanged.
- Against: about 1 GB on disk (the scientific stack); a console window on Windows (as FIVE
  Atlas); the viewer opens in the browser, not inside an app window.

**(c) Other**: a conda-constructor installer or "embedded Python + .bat" (lighter to build but
not relocatable on macOS and no single double-click app), Briefcase (packages a native app, but
its macOS support for numba/llvmlite is less proven than PyInstaller's), or Electron/Tauri
around the page (adds Node or Rust plus still has to ship Python).

**Recommendation: (b), PyInstaller, built exactly like FIVE Atlas.** It is the lab's proven
recipe, it is the one that needs nothing new on the build machines, and the user experience is
the same as FIVE Atlas's (download, double-click, the browser opens). This is what is built.

## What the user sees (no coding)

1. **Download** from the Release page: `SODA-<ver>-win64.zip` (Windows 10/11) or
   `SODA-<ver>-macos-arm64.dmg` (Apple Silicon Mac).
2. **Install**: Windows: unzip, copy the `SODA` folder to the local disk, double-click `SODA.exe`
   (first time: More info -> Run anyway). Mac: drag SODA to Applications; first time, the
   one-line Terminal fix or Open Anyway (until the app is notarised, see below).
3. **The SODA home page opens in the browser** (`http://127.0.0.1:8770/app/`) with three steps:
   - **Prepare a folder of datasets**: a data folder (one sub-folder per dataset, or loose
     OME-TIFF / QPTIFF / .h5ad files) and an output folder, each with a Browse... button (a
     native folder picker) or a pasted path. Prepare runs `prepare.py batch` on it: every
     dataset is recognised (`detect.py`), prepared and exported, one process per dataset, then
     `datasets.json` and the **sample ID check**. Progress streams onto the page. The result is
     a table of datasets (ok / failed / skipped, with the reason and what to add), and a
     green "clean: the export folder can be shared" or a red "NOT clean: do not upload or share
     the export folder yet" banner. Nothing in the data folder is changed. Output:
     `<output>/export` (what gets uploaded), `<output>/cache`, `<output>/logs`.
   - **Open viewer**: the viewer opens in a new tab on the export folder, served by the app on
     127.0.0.1 (`?data=` on the folder's `datasets.json`, every dataset in the list).
   - **Check a folder for sample IDs**: `deid.py check` on any folder, with each hit shown by
     place and the ID itself hidden.
4. **Quit**: close the console window (Windows) or press Quit SODA on the page (Mac; it also
   stops by itself after 30 minutes with no page open, as FIVE Atlas does).

The same program also works from a terminal, for scripts and for support:

```
SODA prepare --src <data folder> --work <output folder> [any prepare.py batch option]
SODA view <export folder> [--port N] [--no-browser]
SODA check <folder> ...
SODA selftest [--quick]
SODA prepare.py detect <folder>     (any of SODA's scripts, as `python <script>` would run it)
```

## How it is built

| Piece | File |
|---|---|
| launcher: home page server, `prepare` / `view` / `check` / `selftest`, python stand-in | `app/soda_app.py` |
| the home page | `app/home.html` |
| PyInstaller spec (Windows and macOS) | `app/SODA.spec` |
| build scripts | `app/build_release.ps1`, `app/build_release_mac.sh` |
| Windows installer (optional, needs Inno Setup 6) | `app/installer.iss` |
| README inside the zip / dmg | `app/readme_windows.txt`, `app/readme_mac.txt` |
| macOS: hardened-runtime entitlements, min-OS check | `app/entitlements.plist`, `app/check_min_os.py` |
| CI | `.github/workflows/build-windows.yml`, `.github/workflows/build-macos.yml` |
| tests (from source or a built app; no window opened) | `tests/app_test.py`, `tests/app_browser_test.js` |

Design points, each a deliberate difference from a plain PyInstaller app:

- **SODA's own scripts ship as source.** `web/SpatialViewer` is copied whole to
  `_internal/soda/` (minus `config.js`, `logos/`, `.password` and key files), and SODA's modules
  are taken back out of PyInstaller's archive. So `prepare.py`, `batch.py`, `deid.py` run
  byte-for-byte as in the repository, `markers/*.json` are found next to them, and the page is
  served from the same folder. The third-party libraries are still found by analysing those
  modules, so a new import in `web/SpatialViewer` needs no spec change.
- **The app is also a python.** `batch.py` starts every dataset in its own process with
  `sys.executable`; frozen, that is `SODA.exe`, so `SODA <script>.py args` runs the script
  like python would. That keeps batch's isolation (one dataset running out of memory does not
  stop the others) without changing batch.py.
- **umap-learn and pynndescent ship as source** and numba's cache goes to the user's app-data
  folder: numba's `cache=True` needs real source files and a writable cache.
- **One local server** for the home page, the viewer page and the export folders (FIVE Atlas
  also serves everything from one port). It listens on 127.0.0.1 only, refuses other Host
  headers (DNS rebinding), needs a per-launch token for every action (start, stop, pick, check,
  quit), never serves the `.py` tools, never serves outside a chosen export folder, and keeps
  no request log (paths can name datasets). The viewer gets an empty `config.js`, so no site
  configuration is needed for local viewing.
- **Version** from the git tag, baked in at build time (`app/_buildinfo.py`, git-ignored), the
  same number in the file name, `SODA version`, the home page and the Mac `Info.plist`.

### Build it yourself

Windows (from the repository folder, after `install.ps1`):

```
pip install pyinstaller==6.21.0          (into .venv)
powershell -ExecutionPolicy Bypass -File app\build_release.ps1 -Version 0.1.0 -OutRoot D:\builds\soda
```

macOS (on an Apple Silicon Mac, after `install.sh`):

```
.venv/bin/pip install pyinstaller==6.21.0
./app/build_release_mac.sh --version 0.1.0
```

Then `python tests/app_test.py --exe <...>/SODA/SODA.exe` (or
`.../SODA.app/Contents/MacOS/SODA`) checks the built app without opening a window.

## Size and start-up (measured on Windows)

MEASURED_TABLE

## The sample ID check stays in the flow

The notebook ends with the sample ID check (Section 5) and makes it a gate when a
de-identification key exists. The app is stricter, because its users are the ones least likely
to read a warning in a log:

- Every Prepare runs `prepare.py batch` with `--check-ids`, so the check is **always a gate**,
  with or without a key file (batch exit code 3 when it finds something). There is no switch to
  turn it off in the app.
- The result is the most visible thing on the page: green "clean, can be shared" or red "NOT
  clean, do not upload or share the export folder yet", with where to look.
- Dataset folders whose names look like sample IDs are skipped by batch before anything is
  written (as in batch.py); the key file is used from the same places as the notebook
  (`SPATIALVIZ_KEY`, else `~/private/soda/sample_id_key.csv`).
- Nothing the app writes holds a sample ID: progress lines, the app's `prepare.log` and the
  Check results go through `batch.scrub` (IDs become `[sample ID]`); the server keeps no request
  log; folder paths are not remembered between launches.
- The app never uploads. Publishing stays a separate, deliberate step (docs/SHARE_YOUR_DATA.md,
  docs/HOSTING.md) for a folder the check has passed.

## macOS: Gatekeeper, notarisation, Apple Silicon and Intel

- **As FIVE Atlas ships today: ad-hoc signed, not notarised.** It runs, but the first launch is
  refused once ("damaged" or "cannot verify"); the README in the .dmg has the fix (the
  `xattr` line, or Open Anyway in System Settings on macOS 15). Do not re-sign the bundle with
  `codesign --deep`: FIVE Atlas found that this makes the kernel kill the app at launch.
- **Notarised (recommended for people outside the lab)** needs an Apple Developer Program
  membership (USD 99 a year, in the lab's or the university's name) and a "Developer ID
  Application" certificate. `build_release_mac.sh` then signs every binary with the hardened
  runtime and `app/entitlements.plist` (numba/llvmlite need JIT and unsigned executable memory),
  submits the .dmg to Apple's notary service and staples the ticket. The workflow does this
  when the `MACOS_*` secrets exist (listed at the top of `build-macos.yml`); without them it
  falls back to the ad-hoc build. The app then opens with a normal double-click.
- **Apple Silicon (arm64) only**, as FIVE Atlas: PyInstaller cannot make a universal app unless
  every wheel is universal (imagecodecs, llvmlite are not). An Intel build would need the
  macos-13 runner (which FIVE Atlas found unusable) or an Intel Mac running the same script.
- **Minimum macOS**: `LSMinimumSystemVersion` is 12.0 by default (`SODA_MIN_MACOS`);
  `app/check_min_os.py` fails the CI build if any library inside needs newer, so the plist never
  promises an older macOS than the binaries support.

## Windows: SmartScreen and signing

Unsigned, as FIVE Atlas: "Windows protected your PC" -> More info -> Run anyway, once. A
code-signing certificate (an OV or EV certificate from a CA, or Azure Trusted Signing, about USD
10 a month) would remove that; `signtool sign` would go after PyInstaller in
`build_release.ps1`. The installer (`installer.iss`) is per-user, so no administrator rights are
needed on managed lab PCs; it is built only where Inno Setup 6 is installed.

## Updates

- **Now**: the build bakes in the repository it was built from (`-UpdateRepo` / the workflow's
  `github.repository`). At start-up the app asks that repository's GitHub `releases/latest` once
  (4-second timeout, only its version is sent) and the home page shows "SODA x is out: Download
  it". Offline, or a private repository, simply means no banner.
- **Next, as FIVE Atlas 1.3.3**: one-click self-update (download the zip/dmg in the background,
  check it against the release's SHA256SUMS, swap on restart, roll back on failure). FIVE Atlas's
  `backend/selfupdate.py` can be adapted almost as is; it needs the downloads to be on a public
  Release (a private repository's assets need a token).
- Prepared data is never touched by an update: it lives in the folders people chose.

## What is still open

- **macOS build**: not run. It needs either a GitHub Actions run (Actions tab -> "Build SODA app
  (macOS)" -> Run workflow; in FIVE-MGI repositories the workflows only run by hand) or an Apple
  Silicon Mac with `install.sh` + `pip install pyinstaller==6.21.0` + `app/build_release_mac.sh`.
- **Where releases live**: FIVE Atlas releases on a public repository. If SODA's repository
  stays private, Release downloads need a GitHub login with access, and the update check sees
  nothing. Options: make releases from a public repository, or attach the zip/dmg to a lab share.
- **Notarisation and Windows signing**: need the certificates above.
- **Offline viewing**: the viewer loads deck.gl (and umap-js for the circle tool) from
  cdn.jsdelivr.net, so "Open viewer" needs internet, as the hosted page does. Bundling those
  two files into the app would make it fully offline.
- **Self-update**, and an app icon (both `icon=None` today, as FIVE Atlas).
