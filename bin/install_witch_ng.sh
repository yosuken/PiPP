#!/usr/bin/env bash
# Download the witch-ng binary into the directory this script lives in.
#
# witch-ng is the default --aligner in PiPP but is not packaged on Bioconda
# or conda-forge, so we fetch the prebuilt binary from upstream releases.
#   https://github.com/RuneBlaze/WITCH-NG/releases
#
# Usage:
#   ./bin/install_witch_ng.sh                # latest (currently v0.0.4)
#   WITCH_NG_VERSION=v0.0.4 ./bin/install_witch_ng.sh
#
# License: witch-ng is GPLv3. We do not redistribute the binary; we just
# fetch it on the user's behalf.

set -euo pipefail

VERSION="${WITCH_NG_VERSION:-v0.0.4}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEST="$HERE/witch-ng"

uname_s=$(uname -s)
uname_m=$(uname -m)

case "$uname_s/$uname_m" in
  Linux/x86_64)   asset="witch-ng-x86_64-unknown-linux-gnu.tar.gz" ;;
  Linux/aarch64)  asset="witch-ng-aarch64-unknown-linux-gnu.tar.gz" ;;
  Darwin/x86_64)  asset="witch-ng-x86_64-apple-darwin.tar.gz" ;;
  Darwin/arm64)   asset="witch-ng-aarch64-apple-darwin.tar.gz" ;;
  *)
    echo "Unsupported platform: $uname_s/$uname_m" >&2
    echo "Supported: Linux/x86_64, Linux/aarch64, Darwin/x86_64, Darwin/arm64." >&2
    echo "Build from source instead: https://github.com/RuneBlaze/WITCH-NG" >&2
    exit 2
    ;;
esac

url="https://github.com/RuneBlaze/WITCH-NG/releases/download/${VERSION}/${asset}"
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

echo "==> downloading $asset (${VERSION})"
curl -fSL -o "$tmp/$asset" "$url"

echo "==> extracting"
tar -xzf "$tmp/$asset" -C "$tmp"

# Locate the binary inside the extracted tree (release tarballs put it at
# the top level for some versions, in a subdir for others).
src=$(find "$tmp" -maxdepth 3 -type f -name 'witch-ng' -perm -u+x 2>/dev/null \
       | head -n 1 || true)
if [[ -z "$src" ]]; then
  src=$(find "$tmp" -maxdepth 3 -type f -name 'witch-ng' | head -n 1 || true)
fi
if [[ -z "$src" ]]; then
  echo "Could not find witch-ng binary inside $asset" >&2
  ls -la "$tmp" >&2
  exit 1
fi

install -m 0755 "$src" "$DEST"
echo "==> installed $DEST"
"$DEST" --version 2>&1 | head -1 || true
