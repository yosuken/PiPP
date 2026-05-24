#!/usr/bin/env bash
# Install PiPP dependencies into a conda environment and build the bundled
# Rust binary. This is a thin wrapper around `environment.yaml` and
# `cargo build`; you can run those two commands directly if you prefer.
#
#   $ ./install.sh                     # use micromamba (default)
#   $ CONDA=mamba ./install.sh         # or mamba / conda
#   $ ./install.sh --skip-rust         # don't build pipp_util (e.g. CI without Rust)
#   $ ./install.sh --skip-env          # don't create the conda env
#
# Notes:
#  - The runtime env is created with `--channel-priority flexible`. PiPP pulls
#    packages (e.g. pplacer, fasttree) whose deps cannot be reconciled with the
#    latest conda-forge toolchain under strict priority, so a strict solve fails
#    regardless of your global condarc. flexible is required.
#  - rust is NOT in environment.yaml (build-time only). This script builds
#    pipp_util with a system `cargo` if one is on PATH (rustup recommended);
#    otherwise it creates a throwaway conda env just for the build.
#  - duckdb is NOT in environment.yaml either (pipp_util bundles its own DuckDB).
#    A duckdb CLI (>=1.0) is optional, only for querying the output DBs.

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

CONDA="${CONDA:-micromamba}"
skip_env=0
skip_rust=0
for arg in "$@"; do
  case "$arg" in
    --skip-env)  skip_env=1 ;;
    --skip-rust) skip_rust=1 ;;
    -h|--help)
      sed -n '2,20p' "$0" | sed 's/^# \?//'
      exit 0
      ;;
    *)
      echo "unknown arg: $arg" >&2
      exit 2
      ;;
  esac
done

if [[ $skip_env -eq 0 ]]; then
  echo "==> creating runtime conda env from environment.yaml (via $CONDA, flexible priority)"
  "$CONDA" env create -f environment.yaml --channel-priority flexible
  env_name="$(awk '/^name:/{print $2; exit}' environment.yaml)"
  echo
  echo "Activate the env with:"
  echo "  $CONDA activate $env_name"
fi

if [[ $skip_rust -eq 0 ]]; then
  echo "==> building pipp_util (cargo build --release)"
  if command -v cargo >/dev/null 2>&1; then
    echo "    using system cargo: $(command -v cargo)"
    cargo build --release --manifest-path rust/Cargo.toml
  else
    echo "    no system cargo on PATH; creating a throwaway conda Rust env just to build"
    build_root="$(mktemp -d)"
    build_env="$build_root/rust"
    trap '"$CONDA" env remove -y -p "$build_env" >/dev/null 2>&1 || true; rm -rf "$build_root"' EXIT
    "$CONDA" create -y -p "$build_env" -c conda-forge 'rust>=1.70'
    "$CONDA" run -p "$build_env" cargo build --release --manifest-path rust/Cargo.toml
  fi
  echo
  echo "Binary built at: rust/target/release/pipp_util"
fi
