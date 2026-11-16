#!/usr/bin/env bash
# SODA one-command install for macOS and Linux.
#
# From the repository folder, in Terminal:
#
#     bash install.sh
#
# What it does:
#   1. finds Python 3.12 or 3.13 (python3.13, python3.12, python3, python on the PATH)
#   2. makes a virtual environment in .venv inside the repository (or reuses it)
#   3. upgrades pip in it
#   4. installs the exact package versions in web/SpatialViewer/requirements-lock.txt
#   5. checks that every package imports (tools/check_install.py, about 10 seconds; up to a
#      minute right after a fresh install)
#
# No administrator rights (sudo) are needed; everything goes into .venv. Safe to run again:
# it reuses .venv and only changes packages whose version differs from the lock file.
#
# Options:
#   --python <path>   use this Python instead of searching for one (or set PYTHON=<path>)
#   --venv <folder>   make the environment here instead of .venv (or set VENV=<folder>)
#
# The lock file covers macOS 12+ on Apple Silicon (M1 and newer) and Linux on x86_64 (glibc
# 2.28+). Not yet run on a Mac by us. Intel Macs
# and Linux on ARM are not covered: there, use the manual install in the README.

set -u

REPO="$(cd "$(dirname "$0")" && pwd)"
LOCK="$REPO/web/SpatialViewer/requirements-lock.txt"
CHECK="$REPO/tools/check_install.py"
PYTHON="${PYTHON:-}"
VENV="${VENV:-}"

while [ $# -gt 0 ]; do
  case "$1" in
    --python) PYTHON="${2:-}"; shift 2 ;;
    --venv)   VENV="${2:-}";   shift 2 ;;
    -h|--help) sed -n '2,/^$/p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "Unknown option: $1 (see: bash install.sh --help)"; exit 2 ;;
  esac
done
[ -n "$VENV" ] || VENV="$REPO/.venv"

step() { printf '\n== %s\n' "$1"; }
fail() { printf '\nSTOPPED: %s\n' "$1"; exit 1; }

# Prints "major.minor|bits|machine|platform|path" for a Python, or nothing if it does not run.
probe() {
  "$@" -c 'import sys, struct, platform; print("%d.%d|%d|%s|%s|%s" % (sys.version_info[0], sys.version_info[1], struct.calcsize("P") * 8, platform.machine(), sys.platform, sys.executable))' 2>/dev/null | tail -n 1
}
supported() { [ "$1" = "3.12" ] || [ "$1" = "3.13" ]; }

[ -f "$LOCK" ] || fail "Cannot find $LOCK. Run this script from the SODA repository folder."

# ---------------------------------------------------------------- 1. Python
step "1/5 Looking for Python 3.12 or 3.13"
chosen=""; chosen_info=""; found=""
try_python() {
  [ -z "$chosen" ] || return 0
  command -v "$1" >/dev/null 2>&1 || [ -x "$1" ] || return 0
  info="$(probe "$1")"
  [ -n "$info" ] || return 0
  ver="${info%%|*}"
  found="$found
  Python $ver at ${info##*|}"
  if supported "$ver"; then chosen="$1"; chosen_info="$info"; fi
}
if [ -n "$PYTHON" ]; then
  try_python "$PYTHON"
else
  for c in python3.13 python3.12 python3 python; do try_python "$c"; done
fi

if [ -z "$chosen" ]; then
  if [ -n "$found" ]; then
    echo "Found these Pythons, but none is Python 3.12 or 3.13:$found"
  else
    echo "No Python found."
  fi
  echo
  echo "SODA needs Python 3.12 or 3.13 (newer and older versions do not have all the packages yet)."
  echo "  macOS: install Python 3.13 from https://www.python.org/downloads/macos/"
  echo "         (or with Homebrew: brew install python@3.13)"
  echo "  Linux: install python3.12 or python3.13 and its venv module with your package manager"
  echo "Then run this script again. For a Python somewhere else: bash install.sh --python /path/to/python3.13"
  fail "no suitable Python."
fi

IFS='|' read -r ver bits machine plat path <<EOF
$chosen_info
EOF
echo "Using Python $ver at $path ($plat, $machine)"

[ "$bits" = "64" ] || fail "This Python is $bits-bit; SODA needs a 64-bit Python."
case "$plat|$machine" in
  darwin\|arm64|linux\|x86_64|linux\|amd64) ;;
  win32\|*|cygwin\|*|msys\|*)
    fail "This is Windows: run install.ps1 instead, in PowerShell or Command Prompt:
         powershell -ExecutionPolicy Bypass -File install.ps1" ;;
  darwin\|x86_64)
    fail "This is an Intel Mac, or an Intel (x86_64) Python running under Rosetta on an Apple Silicon Mac.
         The pinned packages have no ready-made versions for Intel Macs. On an Apple Silicon Mac,
         install the arm64 Python from python.org and run this again; on an Intel Mac use the manual
         install in README.md (Requirements)." ;;
  *)
    fail "The pinned packages are tested for macOS on Apple Silicon and Linux on x86_64, not $plat on $machine.
         Use the manual install in README.md (Requirements)." ;;
esac

# ---------------------------------------------------------------- 2. venv
step "2/5 Virtual environment: $VENV"
VPY="$VENV/bin/python"
if [ -e "$VENV" ]; then
  vinfo=""
  [ -x "$VPY" ] && vinfo="$(probe "$VPY")"
  [ -n "$vinfo" ] || fail "The folder $VENV exists but has no working Python in it. Rename or delete that folder, then run this script again."
  vver="${vinfo%%|*}"
  supported "$vver" || fail "$VENV was made with Python $vver; SODA needs 3.12 or 3.13. Rename or delete that folder, then run this script again."
  echo "Already there (Python $vver); reusing it."
else
  if ! "$chosen" -m venv "$VENV"; then
    fail "Could not make the virtual environment in $VENV (see the message above).
         On Debian or Ubuntu the venv module is a separate package: sudo apt install python$ver-venv
         (this one step needs administrator rights), then run this script again."
  fi
  [ -x "$VPY" ] || fail "Could not make the virtual environment in $VENV."
  echo "Made it."
fi

# ---------------------------------------------------------------- 3. pip
step "3/5 Upgrading pip"
"$VPY" -m pip install --disable-pip-version-check --upgrade pip \
  || fail "pip could not be upgraded. Check the internet connection (and proxy settings), then run this script again."

# ---------------------------------------------------------------- 4. packages
step "4/5 Installing the pinned packages (a few minutes the first time)"
"$VPY" -m pip install --disable-pip-version-check -r "$LOCK" \
  || fail "Installing the packages failed; the lines above name the package. Check the internet connection
         and run this script again (it continues where it stopped)."
"$VPY" -m pip check --disable-pip-version-check \
  || fail "pip check found packages that do not fit together (listed above). Please report this with the output."

# ---------------------------------------------------------------- 5. check
step "5/5 Checking that everything imports"
"$VPY" "$CHECK" || fail "The import check failed (see above). Please report this with the output."

case "$VENV" in
  "$REPO"/*) rel="${VENV#"$REPO"/}" ;;
  *) rel="$VENV" ;;
esac
echo
echo "SODA is installed."
echo
echo "Next, from the repository folder, turn the environment on in your window:"
echo "  source $rel/bin/activate"
echo
echo "Then make and prepare the demo data (README, 'Try it in 5 minutes', step 2):"
echo "  python examples/make_demo.py --out demo/raw"
exit 0
