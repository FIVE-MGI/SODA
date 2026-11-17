# The SODA desktop app (Windows and Mac)

The goal: someone with no Python, no terminal and no notebook downloads one file, double-clicks
it, picks a folder of datasets, and gets an export they can look at in the viewer and, once the
sample ID check is clean, upload. It is built and released the way FIVE Atlas is, so the lab has
one way of shipping desktop apps.

Status: the Windows build and its installer are made and tested on a lab PC (numbers below; the
installer installed, tested and uninstalled silently). The macOS workflow
ran once on GitHub (it built; the minimum-macOS check then failed, fixed since: macOS 14 or
newer) and the Windows workflow has not run yet (see [What is still open](#what-is-still-open)).

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

1. **Download** from the Release page: `SODA-<ver>-setup.exe` (Windows 10/11, the installer) or
   `SODA-<ver>-win64.zip` (the same app, no installer), or `SODA-<ver>-macos-arm64.dmg` (Apple
   Silicon Mac, macOS 14 Sonoma or newer).
2. **Install**: Windows: double-click `SODA-<ver>-setup.exe` (first time: More info -> Run
   anyway; no administrator rights needed), then SODA in the Start menu. Or unzip, copy the
   `SODA` folder to the local disk and double-click `SODA.exe`. See
   [Windows: installer or zip](#windows-installer-or-zip-smartscreen-signing). Mac: drag SODA
   to Applications; first time, the one-line Terminal fix or Open Anyway (until the app is
   notarised, see below).
3. **The SODA home page opens in the browser** (`http://127.0.0.1:8770/app/`) with three steps:
   - **Prepare a folder of datasets**: a data folder (one sub-folder per dataset, or loose
     OME-TIFF / QPTIFF / .h5ad files, or OME-Zarr / SpatialData / AnnData `.zarr` folders) and an output folder, each with a Browse... button (a
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
| Windows installer (per-user; built when Inno Setup 6 is found, always in CI) | `app/installer.iss` |
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
  also serves everything from one port). Files are served by `static_serve.py`'s own handler,
  so HTTP Range (206 with Content-Range, which the viewer uses to read one gene out of a gene
  pack) behaves exactly as with `python static_serve.py` and on blob storage. The app's
  subclass only routes `/viewer/` and `/data/<n>/` and tightens it: it listens on 127.0.0.1
  only, refuses other Host headers (DNS rebinding), sends no CORS header (the page is
  same-origin here, and an open one would let any website read the export while SODA runs),
  lists no folders, never serves the `.py` tools or anything outside a chosen export folder,
  needs a per-launch token for every action (start, stop, pick, check, quit), and keeps no
  request log (paths can name datasets). The viewer gets the app's own two-line `config.js`
  (no site configuration is needed for local viewing) with `offlineKeep = false`: the Save tab's
  Offline use section offers no "Keep this dataset for offline use", because the export is
  already a folder on this computer. The page still keeps itself for offline use as on any site.
- **deck.gl and umap-js inside the app** (release builds): the workflows run
  `python web/SpatialViewer/prepare.py vendor` before PyInstaller, which fetches the versions
  pinned in `libs.js` and writes a file only when its size and SHA-384 match; the build scripts
  (`-RequireVendor` / `--require-vendor`) then refuse to build without a complete `vendor/`, and
  the app test (`--expect-vendor`) checks the built app serves both with the pinned hashes and,
  with `--browser`, that the viewer loads deck.gl with the CDN blocked. So Open viewer needs no
  internet. A local build without `vendor/` still works; its viewer loads them from the CDN.
- **Version** from the git tag, baked in at build time (`app/_buildinfo.py`, git-ignored), the
  same number in the file name, `SODA version`, the home page and the Mac `Info.plist`. SODA
  stays at 0.x (`0.x.y`, and `0.0.0-dev` / `0.x.y-dev` for test builds, the default in the build
  scripts and workflows) until the paper; the first public release is 1.0.0.

### Build it yourself

Windows (from the repository folder, after `install.ps1`):

```
pip install pyinstaller==6.21.0          (into .venv)
python web\SpatialViewer\prepare.py vendor   (deck.gl and umap-js into the app; about 1.5 MB from cdn.jsdelivr.net)
powershell -ExecutionPolicy Bypass -File app\build_release.ps1 -Version 0.1.0 -OutRoot D:\builds\soda -RequireVendor
```

The installer is made too when Inno Setup 6 is on the computer. `build_release.ps1` finds
`ISCC.exe` in an all-users install (Program Files or Program Files (x86)), a per-user install
(`%LOCALAPPDATA%\Programs\Inno Setup 6`, what Inno Setup's own installer offers without admin
rights), on PATH, through Inno Setup's uninstall entry in the registry, or where `-Iscc <path>`
or `%ISCC%` says. Without it the build makes the zip only and says so; `-RequireInstaller` makes
that an error instead (the workflow uses it).

macOS (on an Apple Silicon Mac, after `install.sh`):

```
.venv/bin/pip install pyinstaller==6.21.0
.venv/bin/python web/SpatialViewer/prepare.py vendor
./app/build_release_mac.sh --version 0.1.0 --require-vendor
```

**CI.** Besides the two build workflows (run by hand, or on a `v*` tag in a copy of SODA), the
`smoke` workflow (every push and pull request, Ubuntu) runs `tests/app_test.py` from source after
`tests/smoke_test.py`: the launcher's selftest, Prepare on the demo data with the sample ID
check, Check, the viewer server and the home page's API, in about 6 minutes, with no browser.

Then `python tests/app_test.py --exe <...>/SODA/SODA.exe` (or
`.../SODA.app/Contents/MacOS/SODA`) checks the built app without opening a window (it sets
`SODA_NO_BROWSER=1`, so no browser opens, and `SODA_WORKDIR` to its scratch folder, so the real
`%LOCALAPPDATA%\SODA` is not touched). To check the installer the same way, install it silently
into a scratch folder, run the test against the installed `SODA.exe`, and uninstall silently:

```
SODA-0.1.0-setup.exe /VERYSILENT /SUPPRESSMSGBOXES /NORESTART /CURRENTUSER /NOICONS /DIR="D:\scratch\soda-installed"
python tests/app_test.py --exe D:\scratch\soda-installed\SODA.exe
D:\scratch\soda-installed\unins000.exe /VERYSILENT /SUPPRESSMSGBOXES /NORESTART
```

## Size and start-up (measured on Windows)

Built and tested on the lab PC (i7-13700KF, run on the E-cores; Python 3.13.9; PyInstaller
6.21.0), with the build output on D:, which is a network share. FIVE Atlas's figures for
comparison are from its handoff notes.

| | SODA 0.1.0-dev, Windows x64 | FIVE Atlas |
|---|---|---|
| app folder | 448 MB, 2,470 files | about 1,100 files |
| download: zip | 176 MB | |
| download: installer (`setup.exe`, LZMA2) | 117 MB; installed 453 MB, 2,472 files (with the uninstaller) | |
| biggest parts | llvmlite 115 MB (numba's compiler, for UMAP), pyarrow 83 MB, scipy 70 MB, imagecodecs 44 MB, numpy 28 MB; SODA itself 2 MB | |
| build time (`build_release.ps1`) | 7.8 to 9.1 min without the installer; 12.8 to 13.2 min with it (Inno Setup's LZMA2 compression) | "a few minutes" |
| `SODA.exe version` | 1.1 s first time, 0.3 s after | |
| double-click to home page served | 0.4 to 0.6 s (then the browser opens) | about 2 s |
| `selftest` (every library + a UMAP) | 91 s the first time (numba compiles UMAP's code into its cache), seconds after | |
| Prepare the README demo (Visium HD + 6-channel OME-TIFF, two datasets), through the app | 132 to 180 s; from source with python, 144 to 146 s | |

Notes: the first UMAP on a new computer is slow once while numba compiles (the cache then
lives in `%LOCALAPPDATA%\SODA\numba-cache`). Each dataset Prepare runs in its own copy of
`SODA.exe`, which starts in under a second, so the frozen app costs nothing measurable over
python. A macOS build should be about the same size (the same libraries, arm64 builds).

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
  log; folder paths are not remembered between launches. (One known false alarm: running
  Check on the app's own data folder flags two file names in `numba-cache\`, where numba names
  compiled code after a function and a line number, and `tau_rand_int-<line>` reads like an ID.
  That folder is a local cache and is never shared.)
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
- **Apple Silicon (arm64) only, no Intel build** (decided), as FIVE Atlas: PyInstaller cannot
  make a universal app unless every wheel is universal (imagecodecs, llvmlite are not), and an
  Intel build would need the macos-13 runner (which FIVE Atlas found unusable) or an Intel Mac.
- **Minimum macOS: 14 (Sonoma) or newer.** `LSMinimumSystemVersion` is 14.0 by default
  (`SODA_MIN_MACOS`, in `SODA.spec` and `build_release_mac.sh`): the first GitHub build (run
  44114170391, main b7eb968) found that SciPy's arm64 wheels in the lock file need macOS 14.0,
  so the earlier 12.0 promised too much. `app/check_min_os.py` fails the CI build if any library
  inside needs newer, so the plist never promises an older macOS than the binaries support.
  Older wheels are not pinned to reach further back: the lock file stays as it is.

## Windows: installer or zip, SmartScreen, signing

Each Windows release has the same app in two forms:

| | `SODA-<ver>-setup.exe` (installer) | `SODA-<ver>-win64.zip` |
|---|---|---|
| size (0.1.0-dev) | 117 MB download | 176 MB download |
| installed | 453 MB in `%LOCALAPPDATA%\Programs\SODA` | 448 MB wherever the folder is copied |
| admin rights | none (per-user install, no UAC prompt) | none |
| start it | Start menu -> SODA (a desktop shortcut if ticked) | double-click `SODA.exe` in the folder |
| update | run the new setup.exe: it installs over the old version (same folder, one entry in Settings -> Apps) | replace the folder |
| remove | Settings -> Apps -> SODA -> Uninstall, or Start menu -> Uninstall SODA | delete the folder |
| for | most people; IT deployment (silent flags below) | no-install use, a USB stick, people who cannot run installers |

What the installer does and does not do (`app/installer.iss`):

- **Per-user** (`PrivilegesRequired=lowest`): installs into `%LOCALAPPDATA%\Programs\SODA`, writes
  only to the user's own registry hive (the Settings -> Apps entry), needs no administrator
  rights, so it works on managed lab PCs without IT. `/ALLUSERS` on the command line installs for
  everyone into `Program Files\SODA` instead (that one asks for admin rights).
- **Start menu**: SODA and Uninstall SODA, in a SODA folder. A desktop shortcut only if ticked
  (off by default, so a silent install adds none).
- **The uninstaller removes the program folder only.** Prepared data lives in the folders people
  chose; the app's own logs and numba's compiled-code cache in `%LOCALAPPDATA%\SODA`; neither is
  inside the program folder, and neither is touched (checked: see below).
- **Version from the build**: `build_release.ps1` passes the same version that is baked into
  `SODA.exe` (`/DMyAppVersion`); it is the file name, the Settings -> Apps version and the
  setup.exe's file properties.
- **The AppId is fixed** (`{3B0E6C55-...}`, not FIVE Atlas's): every version installs over the
  previous one. Never change it. An update first clears the old `_internal` library folder (only
  if that folder holds a `SODA.exe`), so no file of the old version lingers.
- **Nothing machine-specific inside**: the staging and output folders are build-time defines,
  Inno Setup stores no source paths, and the script names no dataset or sample.
- **Not code-signed**: SmartScreen shows "Windows protected your PC" for the setup.exe (and for
  `SODA.exe` from the zip): click **More info**, then **Run anyway**, once. Browsers may also say
  the download "is not commonly downloaded": Keep. A code-signing certificate (an OV or EV
  certificate from a CA, or Azure Trusted Signing, about USD 10 a month) would remove both;
  `signtool sign` would go after PyInstaller (on `SODA.exe`) and after Inno Setup (on the
  setup.exe, or through Inno's `SignTool=` directive) in `build_release.ps1`.

**Silent install for IT** (no window at all; exit code 0 = installed):

```
SODA-<ver>-setup.exe /VERYSILENT /SUPPRESSMSGBOXES /NORESTART
    /DIR="<folder>"          install somewhere else (default %LOCALAPPDATA%\Programs\SODA)
    /ALLUSERS                for every user, into Program Files (run elevated)
    /CURRENTUSER             for this user only (the default)
    /NOICONS                 no Start-menu entry
    /MERGETASKS=desktopicon  also a desktop shortcut
    /LOG="<file>"            write a log
"<install folder>\unins000.exe" /VERYSILENT /SUPPRESSMSGBOXES /NORESTART      silent uninstall
```

A silent install never starts SODA (the "Start SODA" box at the end is skipped). The uninstall
command is also in the registry, under
`HKCU\Software\Microsoft\Windows\CurrentVersion\Uninstall\{3B0E6C55-41A7-4F7E-9C1D-5D2B8E0A6F31}_is1`
(`HKLM\...` for an `/ALLUSERS` install), as `QuietUninstallString` (that one has `/SILENT`, which
shows a progress bar; use the line above for nothing on screen). An update remembers the first
install's choices: its folder (unless `/DIR=` is given) and whether there is a Start-menu entry
(a later `setup.exe` without `/NOICONS` does not add one to an install made with it).

**Tested** on the lab PC (0.1.0-dev, Inno Setup 6.7.3, all silent, nothing shown on screen):

1. `/VERYSILENT /SUPPRESSMSGBOXES /NORESTART /CURRENTUSER /NOICONS /DIR=<scratch>`: installed in
   12 to 14 s, exit 0; no Start-menu folder, no desktop shortcut; the Settings -> Apps entry
   (HKCU only, none in HKLM) says SODA 0.1.0-dev; the setup.exe's file properties say 0.1.0 /
   0.1.0-dev.
2. `tests/app_test.py --exe <scratch>\SODA.exe --full-selftest --browser` against the installed
   app: 58 passed, 0 failed (including the viewer and home page in headless Chromium, and no
   "Keep this dataset" in the viewer). The app wrote no file into its install folder.
3. The same setup.exe again over it (the update path): a stray file put in `_internal\` was
   removed, the rest reinstalled.
4. `unins000.exe /VERYSILENT /SUPPRESSMSGBOXES /NORESTART`: exit 0, returned after the work was
   done (2 s), the install folder and the Settings -> Apps entry gone; a file put in
   `%LOCALAPPDATA%\SODA` and the app test's 805 files of prepared data were still there.
5. A fresh install without `/NOICONS`: Start menu -> SODA folder with "SODA" and "Uninstall
   SODA" pointing into the install folder; the silent uninstall removed the folder and both
   shortcuts.

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

- **macOS build**: run once on GitHub (FIVE-MGI/SODA run 44114170391, main b7eb968): the Apple
  Silicon build succeeded, then "Check bundle minimum macOS and version" failed (Info.plist said
  12.0; SciPy's wheels need 14.0) and the app test was skipped. The minimum is now 14.0; the
  workflow needs a rerun (Actions tab -> "Build SODA app (macOS)" -> Run workflow; in FIVE-MGI
  repositories the workflows only run by hand) to get past that check and test the app.
- **Windows workflow**: not run yet either. GitHub's `windows-latest` (Windows Server 2025)
  image lists Inno Setup 6.7.1, so the workflow uses it and only falls back to
  `choco install innosetup` if a future image drops it; it then installs, starts and uninstalls
  the setup.exe silently and attaches it to the Release next to the zip.
- **Where releases live**: FIVE Atlas releases on a public repository. If SODA's repository
  stays private, Release downloads need a GitHub login with access, and the update check sees
  nothing. Options: make releases from a public repository, or attach the zip/dmg to a lab share.
- **Notarisation and Windows signing**: need the certificates above.
- **Offline viewing in a local build**: the release workflows now put deck.gl and umap-js in
  the app (above), but the 0.1.0-dev build and installer tested on the lab PC were made without
  `vendor/` (the pinned files were not fetched on that machine), so their viewer loads deck.gl
  from cdn.jsdelivr.net, and keeps it after one "Open viewer" with internet (its service worker,
  see [ADVANCED.md, Offline use](ADVANCED.md#9-offline-use)). The `--expect-vendor` and
  CDN-blocked checks have therefore not run yet against a real vendored build.
- **Self-update**, and an app icon (both `icon=None` today, as FIVE Atlas).
