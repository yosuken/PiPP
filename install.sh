#!/usr/bin/env bash
# Install PiPP with pixi: provision the runtime environment (locked in
# pixi.lock) and build the bundled Rust binary (pipp_util).
#
#   $ ./install.sh                     # pixi install + build pipp_util
#   $ ./install.sh --skip-rust         # don't build pipp_util
#   $ ./install.sh --skip-env          # don't provision the runtime env
#
# Notes:
#  - Needs `pixi` on PATH (https://pixi.sh). The runtime env (bio tools + apples
#    etc.) is pinned by pixi.lock. The Rust toolchain lives in a separate `build`
#    environment so it never co-solves with the runtime deps.
#  - duckdb is NOT a dependency: pipp_util bundles its own DuckDB. A duckdb CLI
#    (>=1.0) is optional, only for querying the output DBs.
#  - Prefer conda? environment.yaml is kept for that; see README.md
#    ("conda / micromamba (compatible)").

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

skip_env=0
skip_rust=0
for arg in "$@"; do
  case "$arg" in
    --skip-env)  skip_env=1 ;;
    --skip-rust) skip_rust=1 ;;
    -h|--help)
      sed -n '2,16p' "$0" | sed 's/^# \?//'
      exit 0
      ;;
    *)
      echo "unknown arg: $arg" >&2
      exit 2
      ;;
  esac
done

if ! command -v pixi >/dev/null 2>&1; then
  echo "error: pixi not found on PATH. Install it from https://pixi.sh" >&2
  echo "       (or use the conda/micromamba path in README.md)" >&2
  exit 1
fi

if [[ $skip_env -eq 0 ]]; then
  echo "==> provisioning runtime env from pixi.lock (pixi install)"
  pixi install
fi

if [[ $skip_rust -eq 0 ]]; then
  echo "==> building pipp_util in the isolated build env (pixi run -e build build)"
  pixi run -e build build
  echo
  echo "Binary built at: rust/target/release/pipp_util"
fi

echo
echo "Run pipp inside the env, e.g. from the repo dir:"
echo "  pixi run ./pipp -q <query.fa> -r <refpkg> -o <out>"
echo "Or put a launcher on PATH so pipp works from anywhere (see README.md)."
