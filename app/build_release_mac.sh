#!/usr/bin/env bash
#
# Build the shippable SODA desktop app for macOS. The counterpart of app/build_release.ps1, and
# of FIVE Atlas's build_release_mac.sh (see docs/DESKTOP_APP.md).
#
# Produces, in <out-root>/release/ :
#   SODA.app                          the app (drag to /Applications)
#   SODA-<ver>-macos-<arch>.dmg       what you send people
#   SODA-<ver>-macos-<arch>.zip       the same, for people who prefer a zip
#   SHA256SUMS                        checksums of both
#
# PyInstaller cannot cross-compile: this must run ON a Mac, and the app is for that Mac's CPU.
# Apple Silicon Mac -> arm64 app. (Intel is not built by default; see the workflow.)
#
# Signing (all optional; without them the app is ad-hoc signed, as FIVE Atlas ships):
#   CODESIGN_ID="Developer ID Application: ..."   sign with hardened runtime + app/entitlements.plist
#   NOTARY_PROFILE=<keychain profile>              also notarise and staple the .dmg
#                                                  (xcrun notarytool store-credentials <profile> ...)
#
# Usage:
#     ./app/build_release_mac.sh
#     ./app/build_release_mac.sh --version 0.2.0 --out-root /Volumes/Scratch/soda --python .venv/bin/python
#
set -euo pipefail

APP="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(dirname "$APP")"
VERSION="$(git -C "$REPO" describe --tags --dirty --always 2>/dev/null | sed 's/^v//')"
VERSION="${VERSION:-0.0.0-dev}"
OUT_ROOT="${HOME}/Library/Caches/SODA_build"
PY=""
UPDATE_REPO="${SODA_UPDATE_REPO:-}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --version)     VERSION="$2"; shift 2 ;;
        --out-root)    OUT_ROOT="$2"; shift 2 ;;
        --python)      PY="$2"; shift 2 ;;
        --update-repo) UPDATE_REPO="$2"; shift 2 ;;
        -h|--help)     sed -n '2,25p' "$0"; exit 0 ;;
        *) echo "unknown option: $1" >&2; exit 2 ;;
    esac
done
printf '%s' "$VERSION" | grep -Eq '^[0-9A-Za-z][0-9A-Za-z.+-]*$' || { echo "version '$VERSION' is not a plain version string" >&2; exit 2; }

if [[ -z "$PY" ]]; then
    if [[ -x "$REPO/.venv/bin/python" ]]; then PY="$REPO/.venv/bin/python"; else PY="$(command -v python3 || true)"; fi
fi
[[ -n "$PY" ]] || { echo "no python found: run install.sh first, or pass --python" >&2; exit 1; }
ARCH="$(uname -m)"
MIN_MACOS="${SODA_MIN_MACOS:-12.0}"

WORK="$OUT_ROOT/work"; DIST_DIR="$OUT_ROOT/dist"; RELEASE="$OUT_ROOT/release"; SCRATCH="$OUT_ROOT/temp"
mkdir -p "$WORK" "$DIST_DIR" "$RELEASE" "$SCRATCH"

echo "== SODA $VERSION (macOS $ARCH) =="
echo "   build output -> $OUT_ROOT"
echo "   python       -> $PY"

# 1. the environment ----------------------------------------------------------------------
echo "[1/6] checking the build environment..."
"$PY" "$REPO/tools/check_install.py"
"$PY" -m PyInstaller --version >/dev/null || { echo "PyInstaller is not installed for $PY (pip install pyinstaller==6.21.0)" >&2; exit 1; }

# 2. PyInstaller ----------------------------------------------------------------------------
echo "[2/6] running PyInstaller (several minutes)..."
pkill -f "SODA.app/Contents/MacOS/SODA" 2>/dev/null || true
export SODA_VERSION="$VERSION" SODA_UPDATE_REPO="$UPDATE_REPO" SODA_MIN_MACOS="$MIN_MACOS"
"$PY" -m PyInstaller --noconfirm --clean --log-level WARN --workpath "$WORK" --distpath "$DIST_DIR" "$APP/SODA.spec"
BUILT="$DIST_DIR/SODA.app"
[[ -d "$BUILT" ]] || { echo "SODA.app was not produced" >&2; exit 1; }

# 3. sign -----------------------------------------------------------------------------------
# PyInstaller has already ad-hoc signed every Mach-O. Do NOT "re-sign to be safe" with
# `codesign --force --deep --sign -`: FIVE Atlas learned that this clobbers the per-binary
# signatures and the kernel then SIGKILLs the app at launch while `codesign --verify` still passes.
echo "[3/6] signing..."
if [[ -n "${CODESIGN_ID:-}" ]]; then
    # Developer ID: hardened runtime, with the entitlements numba / llvmlite's JIT needs.
    # Inner binaries first (no --deep, it is deprecated), then the bundle.
    while IFS= read -r -d '' f; do
        if file -b "$f" | grep -q "Mach-O"; then
            codesign --force --timestamp --options runtime --entitlements "$APP/entitlements.plist" \
                --sign "$CODESIGN_ID" "$f"
        fi
    done < <(find "$BUILT/Contents" -type f \( -name "*.so" -o -name "*.dylib" -o -perm -u+x \) -print0)
    codesign --force --timestamp --options runtime --entitlements "$APP/entitlements.plist" \
        --sign "$CODESIGN_ID" "$BUILT"
    echo "      signed with: $CODESIGN_ID"
else
    echo "      leaving PyInstaller's ad-hoc signature alone (no CODESIGN_ID)"
fi
if codesign --verify --strict "$BUILT"; then echo "      signature verifies"
else echo "the bundle signature does not verify: refusing to package it" >&2; exit 1; fi

# 4. does it start? ---------------------------------------------------------------------------
# Verifying is not running: start the real binary. version must match, every library must load.
echo "[4/6] launch check..."
EXE="$BUILT/Contents/MacOS/SODA"
got="$("$EXE" version 2>&1 || true)"
[[ "$got" == "SODA $VERSION" ]] || { echo "the app says '$got', expected 'SODA $VERSION'" >&2; exit 1; }
rc=0
SODA_WORKDIR="$SCRATCH/launch-check-workdir" "$EXE" selftest --quick >"$SCRATCH/selftest.json" 2>"$SCRATCH/selftest.err" || rc=$?
if [[ $rc -ne 0 ]]; then
    echo "SODA selftest failed (exit $rc)" >&2
    [[ $rc -eq 137 || $rc -eq 9 ]] && echo "  137/9 = SIGKILL: the kernel is rejecting the bundle's signature" >&2
    cat "$SCRATCH/selftest.json" "$SCRATCH/selftest.err" >&2 || true
    exit 1
fi
echo "      starts, version and selftest OK"

STAGED="$RELEASE/SODA.app"
rm -rf "$STAGED"
ditto "$BUILT" "$STAGED"            # ditto, not cp: keeps symlinks and the signature
sed -e "s/{VERSION}/$VERSION/g" -e "s/{ARCH}/$ARCH/g" -e "s/{MIN_MACOS}/$MIN_MACOS/g" \
    "$APP/readme_mac.txt" > "$RELEASE/README.txt"

# 5. dmg (+ notarise) -------------------------------------------------------------------------
echo "[5/6] building dmg..."
DMG="$RELEASE/SODA-$VERSION-macos-$ARCH.dmg"
rm -f "$DMG"
STAGE="$(mktemp -d)"
trap 'rm -rf "$STAGE"' EXIT
ditto "$STAGED" "$STAGE/SODA.app"
cp "$RELEASE/README.txt" "$STAGE/README.txt"
ln -s /Applications "$STAGE/Applications"
hdiutil create -volname "SODA $VERSION" -srcfolder "$STAGE" -ov -format UDZO "$DMG" >/dev/null
hdiutil verify "$DMG" >/dev/null || { echo "dmg does not verify" >&2; exit 1; }
if [[ -n "${CODESIGN_ID:-}" ]]; then
    codesign --force --timestamp --sign "$CODESIGN_ID" "$DMG"
fi
if [[ -n "${CODESIGN_ID:-}" && -n "${NOTARY_PROFILE:-}" ]]; then
    echo "      notarising (a few minutes)..."
    xcrun notarytool submit "$DMG" --keychain-profile "$NOTARY_PROFILE" --wait
    xcrun stapler staple "$DMG"
    xcrun stapler staple "$STAGED" || true
    spctl --assess --type open --context context:primary-signature -v "$DMG" || true
fi

# 6. zip --------------------------------------------------------------------------------------
echo "[6/6] zipping..."
ZIP="$RELEASE/SODA-$VERSION-macos-$ARCH.zip"
rm -f "$ZIP"
ZSTAGE="$(mktemp -d)/SODA-$VERSION"
mkdir -p "$ZSTAGE"
ditto "$STAGED" "$ZSTAGE/SODA.app"
cp "$RELEASE/README.txt" "$ZSTAGE/README.txt"
ditto -c -k --sequesterRsrc --keepParent "$ZSTAGE" "$ZIP"
rm -rf "$(dirname "$ZSTAGE")"

( cd "$RELEASE" && shasum -a 256 "$(basename "$DMG")" "$(basename "$ZIP")" > SHA256SUMS )

echo
echo "Done. In $RELEASE :"
du -sh "$RELEASE"/* 2>/dev/null | sed 's/^/  /'
