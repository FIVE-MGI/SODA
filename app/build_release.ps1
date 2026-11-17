<#
    Build the shippable SODA desktop app for Windows (the counterpart of FIVE Atlas's
    build_release.ps1; see docs/DESKTOP_APP.md).

    Produces, in <OutRoot>\release\ :
      SODA\                      the app folder (double-click SODA.exe)
      SODA-<ver>-win64.zip       that folder, zipped: this is what you send
      SODA-<ver>-setup.exe       a per-user installer, only if Inno Setup is installed
      SHA256SUMS-windows.txt     checksums of the above

    Everything the build writes goes to OutRoot (about 2 GB of scratch), outside the repository,
    so build trees never show up in git status.

    Usage (from the repository folder):
        powershell -ExecutionPolicy Bypass -File app\build_release.ps1
        powershell -ExecutionPolicy Bypass -File app\build_release.ps1 -Version 0.1.0 -OutRoot D:\builds\soda
        ... -Python C:\path\to\python.exe        build with this Python (CI: the runner's own)
        ... -ExtraSitePath D:\tools\site           a folder put first on PYTHONPATH while building,
                                                   e.g. a PyInstaller copied from another venv
        ... -UpdateRepo owner/name                 where the app looks for newer releases (default: none)
#>
param(
    [string]$Version = "",
    [string]$OutRoot = "$env:LOCALAPPDATA\SODA_build",
    [string]$Python = "",
    [string]$ExtraSitePath = "",
    [string]$UpdateRepo = ""
)

$ErrorActionPreference = "Stop"
$ProgressPreference = "SilentlyContinue"      # Compress-Archive is many times slower with a progress bar
$App = $PSScriptRoot
$Repo = Split-Path $App -Parent

# No -Version: the git tag this tree is at (as FIVE Atlas), else 0.0.0-dev.
if (-not $Version) {
    try { $Version = (& git -C $Repo describe --tags --dirty --always 2>$null) } catch { $Version = $null }
    if ($Version) { $Version = ($Version -replace '^v', '').Trim() }
    if (-not $Version) { $Version = "0.0.0-dev" }
}
if ($Version -notmatch '^[0-9A-Za-z][0-9A-Za-z.+-]*$') { throw "version '$Version' is not a plain version string" }

# The Python to build with: -Python, else the repository's .venv (what install.ps1 makes).
if ($Python) {
    $Py = (Get-Command $Python -ErrorAction SilentlyContinue).Source
    if (-not $Py) { $Py = $Python }
} else {
    $Py = Join-Path $Repo ".venv\Scripts\python.exe"
}
if (-not (Test-Path $Py) -and -not (Get-Command $Py -ErrorAction SilentlyContinue)) {
    throw "Python not found: $Py. Run install.ps1 first, or pass -Python."
}

$oldPP = $env:PYTHONPATH
if ($ExtraSitePath) { $env:PYTHONPATH = "$ExtraSitePath;$env:PYTHONPATH" }

$Work = Join-Path $OutRoot "work"
$DistDir = Join-Path $OutRoot "dist"
$Scratch = Join-Path $OutRoot "temp"
$Release = Join-Path $OutRoot "release"
New-Item -ItemType Directory -Force -Path $OutRoot, $Work, $DistDir, $Scratch, $Release | Out-Null

Write-Host "== SODA $Version (Windows) ==" -ForegroundColor Cyan
Write-Host "   build output -> $OutRoot" -ForegroundColor DarkGray
Write-Host "   python       -> $Py" -ForegroundColor DarkGray
$t0 = Get-Date

try {
    # 1. the environment ---------------------------------------------------------------
    Write-Host "[1/5] checking the build environment..." -ForegroundColor Yellow
    & $Py (Join-Path $Repo "tools\check_install.py")
    if ($LASTEXITCODE -ne 0) { throw "tools\check_install.py failed: the packages do not match requirements-lock.txt" }
    $piv = & $Py -m PyInstaller --version
    if ($LASTEXITCODE -ne 0) { throw "PyInstaller is not installed for $Py (pip install pyinstaller==6.21.0)" }
    Write-Host "   PyInstaller $piv" -ForegroundColor DarkGray
    # deck.gl and umap-js, optional: with web\SpatialViewer\vendor\ (python web\SpatialViewer\prepare.py vendor
    # before the build) they go into the app (SODA.spec takes the whole page folder) and Open viewer needs no
    # internet; without it the viewer loads them from cdn.jsdelivr.net. A vendor\ that is there must be complete.
    if (Test-Path (Join-Path $Repo "web\SpatialViewer\vendor")) {
        & $Py (Join-Path $Repo "web\SpatialViewer\vendor.py") --check
        if ($LASTEXITCODE -ne 0) { throw "web\SpatialViewer\vendor\ is incomplete or not the pinned files: run python web\SpatialViewer\prepare.py vendor again" }
        Write-Host "   vendor\: deck.gl and umap-js included (the viewer works with no internet)" -ForegroundColor DarkGray
    } else {
        Write-Host "   no web\SpatialViewer\vendor\: the viewer will load deck.gl from the CDN (prepare.py vendor includes it)" -ForegroundColor DarkGray
    }

    # 2. PyInstaller ---------------------------------------------------------------------
    Write-Host "[2/5] running PyInstaller (several minutes)..." -ForegroundColor Yellow
    Get-Process SODA -ErrorAction SilentlyContinue | Where-Object { $_.Path -like "$OutRoot*" } |
        Stop-Process -Force -ErrorAction SilentlyContinue
    $oldTemp, $oldTmp = $env:TEMP, $env:TMP
    try {
        $env:TEMP = $Scratch; $env:TMP = $Scratch
        $env:SODA_VERSION = $Version
        $env:SODA_UPDATE_REPO = $UpdateRepo
        & $Py -m PyInstaller --noconfirm --clean --log-level WARN `
            --workpath $Work --distpath $DistDir (Join-Path $App "SODA.spec")
        if ($LASTEXITCODE -ne 0) { throw "PyInstaller failed" }
    } finally {
        $env:TEMP = $oldTemp; $env:TMP = $oldTmp
    }
    $Built = Join-Path $DistDir "SODA"
    $Exe = Join-Path $Built "SODA.exe"
    if (-not (Test-Path $Exe)) { throw "SODA.exe was not produced" }

    # 3. does it start? --------------------------------------------------------------------
    # The version inside must be the one asked for, and every library must load (a missing
    # hidden import or DLL fails here, not on a colleague's laptop).
    Write-Host "[3/5] launch check..." -ForegroundColor Yellow
    $got = (& $Exe version) -join ""
    if ($got -ne "SODA $Version") { throw "SODA.exe says '$got', expected 'SODA $Version'" }
    $env:SODA_WORKDIR = Join-Path $Scratch "launch-check-workdir"
    & $Exe selftest --quick | Out-File -Encoding utf8 (Join-Path $Scratch "selftest.json")
    $rc = $LASTEXITCODE
    Remove-Item Env:\SODA_WORKDIR
    if ($rc -ne 0) {
        Get-Content (Join-Path $Scratch "selftest.json") | Write-Host
        throw "SODA.exe selftest --quick failed"
    }
    Write-Host "   starts, version and selftest OK" -ForegroundColor DarkGray

    # 4. stage + zip ------------------------------------------------------------------------
    Write-Host "[4/5] staging and zipping..." -ForegroundColor Yellow
    $Staged = Join-Path $Release "SODA"
    if (Test-Path $Staged) { Remove-Item -Recurse -Force $Staged }
    Copy-Item -Recurse $Built $Staged
    $readme = [IO.File]::ReadAllText((Join-Path $App "readme_windows.txt")).Replace("{VERSION}", $Version)
    [IO.File]::WriteAllText((Join-Path $Staged "README.txt"), $readme, (New-Object Text.UTF8Encoding $false))
    $Zip = Join-Path $Release "SODA-$Version-win64.zip"
    if (Test-Path $Zip) { Remove-Item -Force $Zip }
    Compress-Archive -Path $Staged -DestinationPath $Zip -CompressionLevel Optimal

    # 5. installer (optional, as FIVE Atlas) -----------------------------------------------
    Write-Host "[5/5] installer..." -ForegroundColor Yellow
    $Iscc = $null
    foreach ($c in @("ISCC.exe", "C:\Program Files (x86)\Inno Setup 6\ISCC.exe", "C:\Program Files\Inno Setup 6\ISCC.exe")) {
        $found = if ($c -eq "ISCC.exe") { (Get-Command $c -ErrorAction SilentlyContinue).Source } elseif (Test-Path $c) { $c } else { $null }
        if ($found) { $Iscc = $found; break }
    }
    $sums = @($Zip)
    if ($Iscc) {
        & $Iscc "/DMyAppVersion=$Version" "/DStageDir=$Staged" "/DOutDir=$Release" (Join-Path $App "installer.iss")
        if ($LASTEXITCODE -ne 0) { throw "Inno Setup failed" }
        $sums += (Join-Path $Release "SODA-$Version-setup.exe")
    } else {
        Write-Host "      Inno Setup not found: zip only." -ForegroundColor DarkGray
    }
    $lines = foreach ($f in $sums) { "{0}  {1}" -f (Get-FileHash -Algorithm SHA256 $f).Hash.ToLower(), (Split-Path $f -Leaf) }
    [IO.File]::WriteAllLines((Join-Path $Release "SHA256SUMS-windows.txt"), [string[]]$lines)
} finally {
    $env:PYTHONPATH = $oldPP
}

$mins = ((Get-Date) - $t0).TotalMinutes
Write-Host ""
Write-Host ("Done in {0:N1} min. In {1} :" -f $mins, $Release) -ForegroundColor Green
Get-ChildItem $Release | ForEach-Object {
    $size = if ($_.PSIsContainer) {
        $files = Get-ChildItem $_.FullName -Recurse -File
        "{0:N0} MB, {1:N0} files (folder)" -f (($files | Measure-Object Length -Sum).Sum / 1MB), $files.Count
    } else { "{0:N0} MB" -f ($_.Length / 1MB) }
    "  {0,-34} {1}" -f $_.Name, $size
}
