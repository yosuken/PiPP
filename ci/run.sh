#!/usr/bin/env bash
# CI entry point. Safe to run locally:
#   $ ci/run.sh                # run all phases
#   $ ci/run.sh ruby           # only Ruby syntax
#   $ ci/run.sh rust           # only Rust fmt+clippy+build
#   $ ci/run.sh smoke          # only pipp_util smoke test
#
# Each phase exits non-zero on the first failure. Set CI=1 to add
# extra noise reduction (e.g. cargo --color always).

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

bold() { printf "\033[1m%s\033[0m\n" "$*"; }
green() { printf "\033[1;32m%s\033[0m\n" "$*"; }
red() { printf "\033[1;31m%s\033[0m\n" "$*" >&2; }

phase_ruby() {
  # -w adds warnings (uninitialized vars, ambiguous syntax, etc.) on top of
  # plain `-c` syntax check. We deliberately do not run a full linter
  # (rubocop/standardrb) to avoid style wars; `ruby -wc` catches real bugs.
  bold "=== ruby syntax (ruby -wc) ==="
  ruby -wc PiPP
  ruby -wc PiPP.rake
  for f in script/*.rb; do
    ruby -wc "$f"
  done
  green "ruby OK"
}

phase_rust() {
  bold "=== rust fmt ==="
  cargo fmt --manifest-path rust/Cargo.toml --check

  bold "=== rust clippy ==="
  cargo clippy --manifest-path rust/Cargo.toml --release --locked -- -D warnings

  bold "=== rust build ==="
  cargo build --manifest-path rust/Cargo.toml --release --locked

  bold "=== rust test ==="
  cargo test --manifest-path rust/Cargo.toml --release --locked

  green "rust OK"
}

phase_smoke() {
  bold "=== pipp_util smoke test ==="
  "$ROOT/ci/smoke_pipp_util.sh"
  green "smoke OK"
}

phases=("${@:-ruby rust smoke}")
for p in ${phases[@]}; do
  case "$p" in
    ruby)  phase_ruby ;;
    rust)  phase_rust ;;
    smoke) phase_smoke ;;
    *)     red "unknown phase: $p"; exit 2 ;;
  esac
done

green "ALL OK"
