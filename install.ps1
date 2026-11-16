<#
SODA one-command install for Windows.

From the repository folder, in PowerShell or Command Prompt:

    powershell -ExecutionPolicy Bypass -File install.ps1

("-ExecutionPolicy Bypass" lets this one script run; it changes no setting on the computer.)

What it does:
  1. finds Python 3.12 or 3.13 (the py launcher first, then python on the PATH)
  2. makes a virtual environment in .venv inside the repository (or reuses it)
  3. upgrades pip in it
  4. installs the exact package versions in web\SpatialViewer\requirements-lock.txt
  5. checks that every package imports (tools\check_install.py, about 10 seconds; up to a
     minute right after a fresh install)

No administrator rights are needed; everything goes into .venv. Safe to run again: it reuses
.venv and only changes packages whose version differs from the lock file.

Options:
  -Python <path to python.exe>   use this Python instead of searching for one
  -Venv <folder>                 make the environment here instead of .venv
#>
param(
    [string]$Python = "",
    [string]$Venv = ""
)

$ErrorActionPreference = "Continue"
$Repo = $PSScriptRoot
$Lock = Join-Path $Repo "web\SpatialViewer\requirements-lock.txt"
$Check = Join-Path $Repo "tools\check_install.py"
if (-not $Venv) { $Venv = Join-Path $Repo ".venv" }
$Supported = @("3.12", "3.13")

function Say($text) { Write-Host $text }
function Step($text) { Write-Host ""; Write-Host "== $text" -ForegroundColor Cyan }
function Fail($text) {
    Write-Host ""
    Write-Host "STOPPED: $text" -ForegroundColor Red
    exit 1
}

# Ask a Python for its version, bitness and path. Returns $null if it does not run.
function Probe-Python($exe, $pre) {
    $code = "import sys, struct; print('%d.%d|%d|%s' % (sys.version_info[0], sys.version_info[1], struct.calcsize('P') * 8, sys.executable))"
    try {
        $out = & $exe @pre -c $code 2>$null
    } catch {
        return $null
    }
    if ($LASTEXITCODE -ne 0 -or -not $out) { return $null }
    $parts = ([string]($out | Select-Object -Last 1)).Trim().Split("|")
    if ($parts.Count -lt 3) { return $null }
    return [pscustomobject]@{ Version = $parts[0]; Bits = $parts[1]; Path = $parts[2]; Exe = $exe; Pre = $pre }
}

if (-not (Test-Path $Lock)) {
    Fail "Cannot find $Lock. Run this script from the SODA repository folder (it must stay next to the web folder)."
}

# ---------------------------------------------------------------- 1. Python
Step "1/5 Looking for Python 3.12 or 3.13"
$candidates = @()
if ($Python) {
    $candidates += ,@($Python, @())
} else {
    if (Get-Command py -ErrorAction SilentlyContinue) {
        $candidates += ,@("py", @("-3.13"))
        $candidates += ,@("py", @("-3.12"))
    }
    foreach ($name in @("python", "python3")) {
        $cmd = Get-Command $name -ErrorAction SilentlyContinue | Select-Object -First 1
        # The WindowsApps "python" is only a shortcut to the Microsoft Store, not a Python.
        if ($cmd -and $cmd.Source -notlike "*\WindowsApps\*") { $candidates += ,@($cmd.Source, @()) }
    }
}

$found = @()
$chosen = $null
foreach ($c in $candidates) {
    $p = Probe-Python $c[0] $c[1]
    if (-not $p) { continue }
    $found += $p
    if (($Supported -contains $p.Version) -and $p.Bits -eq "64") { $chosen = $p; break }
}

if (-not $chosen) {
    if ($found.Count -gt 0) {
        Say "Found these Pythons, but none is a 64-bit Python 3.12 or 3.13:"
        foreach ($p in $found) { Say ("  Python {0} ({1}-bit) at {2}" -f $p.Version, $p.Bits, $p.Path) }
    } elseif ($Python) {
        Say "Could not run $Python."
    } else {
        Say "No Python found."
    }
    Fail ("SODA needs 64-bit Python 3.12 or 3.13 (newer and older versions do not have all the packages yet).`n" +
          "         Install Python 3.13 from https://www.python.org/downloads/windows/ (the default 64-bit installer;`n" +
          "         no administrator rights needed), then run this script again. If it is installed somewhere`n" +
          "         unusual: powershell -ExecutionPolicy Bypass -File install.ps1 -Python C:\path\to\python.exe")
}
Say ("Using Python {0} at {1}" -f $chosen.Version, $chosen.Path)

# ---------------------------------------------------------------- 2. venv
Step "2/5 Virtual environment: $Venv"
$VenvPy = Join-Path $Venv "Scripts\python.exe"
if (Test-Path $Venv) {
    $existing = $null
    if (Test-Path $VenvPy) { $existing = Probe-Python $VenvPy @() }
    if (-not $existing) {
        Fail "The folder $Venv exists but has no working Python in it. Rename or delete that folder, then run this script again."
    }
    if (-not ($Supported -contains $existing.Version)) {
        Fail ("$Venv was made with Python {0}; SODA needs 3.12 or 3.13. Rename or delete that folder, then run this script again." -f $existing.Version)
    }
    Say ("Already there (Python {0}); reusing it." -f $existing.Version)
} else {
    $pre = @($chosen.Pre | Where-Object { $_ })
    & $chosen.Exe @pre -m venv $Venv
    if ($LASTEXITCODE -ne 0 -or -not (Test-Path $VenvPy)) {
        Fail "Could not make the virtual environment in $Venv (see the message above)."
    }
    Say "Made it."
}

# ---------------------------------------------------------------- 3. pip
Step "3/5 Upgrading pip"
& $VenvPy -m pip install --disable-pip-version-check --upgrade pip
if ($LASTEXITCODE -ne 0) {
    Fail "pip could not be upgraded. Check the internet connection (and proxy settings), then run this script again."
}

# ---------------------------------------------------------------- 4. packages
Step "4/5 Installing the pinned packages (a few minutes the first time)"
& $VenvPy -m pip install --disable-pip-version-check -r $Lock
if ($LASTEXITCODE -ne 0) {
    Fail ("Installing the packages failed; the lines above name the package. Check the internet connection`n" +
          "         and run this script again (it continues where it stopped).")
}
& $VenvPy -m pip check --disable-pip-version-check
if ($LASTEXITCODE -ne 0) {
    Fail "pip check found packages that do not fit together (listed above). Please report this with the output."
}

# ---------------------------------------------------------------- 5. check
Step "5/5 Checking that everything imports"
& $VenvPy $Check
if ($LASTEXITCODE -ne 0) {
    Fail "The import check failed (see above). Please report this with the output."
}

$rel = $Venv
if ($Venv.StartsWith($Repo, [System.StringComparison]::OrdinalIgnoreCase)) {
    $rel = $Venv.Substring($Repo.Length).TrimStart("\")
}
Write-Host ""
Write-Host "SODA is installed." -ForegroundColor Green
Say ""
Say "Next, from the repository folder, turn the environment on in your window:"
Say "  PowerShell:      $rel\Scripts\Activate.ps1"
Say "  Command Prompt:  $rel\Scripts\activate.bat"
Say "(If PowerShell refuses, run first: Set-ExecutionPolicy -Scope Process Bypass  - it lasts only for that window.)"
Say ""
Say "Then make and prepare the demo data (README, 'Try it in 5 minutes', step 2):"
Say "  python examples/make_demo.py --out demo/raw"
exit 0
