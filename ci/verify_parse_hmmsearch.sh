#!/usr/bin/env bash
# Compare the Ruby (script/parse_hmmsearch.rb) and Rust
# (pipp_util parse-hmmsearch) implementations on one PPP run.
#
# Usage:
#   ci/verify_parse_hmmsearch.sh <run_dir>
#
# run_dir is something like
#   /aptmp/yosuke/www_db/uniparc/data/2024-02-06/PPP/001/rhodopsin_2020-11-17
# i.e. the directory that contains `prefilter/hmmsearch/<sub>/sub*.out`
# and `prefilter/query/*.fa`.
#
# Outputs both implementations into:
#   $TMPDIR/<run_basename>/{ruby,rust}/
# and writes a diff summary to stdout.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RUBY_SCRIPT="$ROOT/script/parse_hmmsearch.rb"
RUST_BIN="${PIPP_UTIL_BIN:-$ROOT/rust/target/release/pipp_util}"

run_dir="${1:?usage: verify_parse_hmmsearch.sh <run_dir>}"
[[ -d "$run_dir" ]] || { echo "no such dir: $run_dir" >&2; exit 1; }

# Locate inputs
hmm_dir=$(find "$run_dir/prefilter/hmmsearch" -mindepth 1 -maxdepth 1 -type d | head -n1)
[[ -n "$hmm_dir" ]] || { echo "no prefilter/hmmsearch/<sub> under $run_dir" >&2; exit 1; }
query_fa=$(find "$run_dir/prefilter/query" -maxdepth 1 -name '*.fa' | head -n1)
[[ -n "$query_fa" ]] || { echo "no prefilter/query/*.fa under $run_dir" >&2; exit 1; }

# PiPP defaults
EVALUE=1e-5            # gene-level (passed to -ge / --gene-evalue)
EVALUEDOM=1e-2         # domain-level (passed to -e / --evalue)

# Working area
tmp_base="${VERIFY_TMP:-/tmp/pipp_verify}"
run_name=$(basename "$(dirname "$hmm_dir")")_$(basename "$run_dir")    # e.g. uniparc_active_p1.UO2X_rhodopsin...
tmp="$tmp_base/$(basename "$(dirname "$run_dir")")_${run_name}"
rm -rf "$tmp"
mkdir -p "$tmp/ruby" "$tmp/rust"

# Concatenate sub*.out (same as pipp.rake 01-2b)
concat="$tmp/hmmsearch_concat.out"
cat "$hmm_dir"/sub*.out > "$concat"
n_subout=$(ls "$hmm_dir"/sub*.out | wc -l)
n_lines=$(wc -l < "$concat")
echo "[$(date +%T)] inputs: hmmsearch_concat from $n_subout sub*.out ($n_lines lines); query $(du -h "$query_fa" | cut -f1)"

# Common args
common="--create-evalue-table \
  --min-hmm-len-dom 0 --min-hmm-cov-dom 0.0 --min-ali-len-dom 0 --min-ali-cov-dom 0.0 \
  --min-hmm-len 0     --min-hmm-cov 0.0     --min-ali-len 0     --min-ali-cov 0.0"

# Ruby run
echo "[$(date +%T)] ruby ..."
t0=$(date +%s)
ruby "$RUBY_SCRIPT" \
  -ge "$EVALUE" -e "$EVALUEDOM" $common \
  -i "$concat" -f "$query_fa" -o "$tmp/ruby" >"$tmp/ruby.log" 2>&1
t_ruby=$(( $(date +%s) - t0 ))
echo "[$(date +%T)] ruby done in ${t_ruby}s"

# Rust run
echo "[$(date +%T)] rust ..."
t0=$(date +%s)
"$RUST_BIN" parse-hmmsearch \
  --gene-evalue "$EVALUE" --evalue "$EVALUEDOM" $common \
  -i "$concat" -f "$query_fa" -o "$tmp/rust" >"$tmp/rust.log" 2>&1
t_rust=$(( $(date +%s) - t0 ))
echo "[$(date +%T)] rust done in ${t_rust}s"

# ---- compare ----
status=0
echo "----- diff summary -----"
for f in all-hit.tsv best-hit.tsv best-hit.whole.fa best-hit.fa evalues.tsv; do
  ruby_p="$tmp/ruby/$f"
  rust_p="$tmp/rust/$f"
  if [[ ! -f "$ruby_p" && ! -f "$rust_p" ]]; then
    printf "  %-20s -- both missing\n" "$f"
    continue
  fi
  if [[ ! -f "$ruby_p" || ! -f "$rust_p" ]]; then
    printf "  %-20s ONLY %s exists  FAIL\n" "$f" "$([[ -f "$ruby_p" ]] && echo ruby || echo rust)"
    status=1
    continue
  fi

  r_n=$(wc -l < "$ruby_p")
  s_n=$(wc -l < "$rust_p")

  # Try byte-exact first, then sorted line-set equality (robust to ordering).
  if cmp -s "$ruby_p" "$rust_p"; then
    printf "  %-20s rows=%s  BYTE-IDENTICAL\n" "$f" "$r_n"
  elif diff -q <(sort "$ruby_p") <(sort "$rust_p") >/dev/null; then
    printf "  %-20s rows=%s  IDENTICAL after sort\n" "$f" "$r_n"
  else
    printf "  %-20s rows ruby=%s rust=%s  DIFF\n" "$f" "$r_n" "$s_n"
    # Save unified diff (first 50 lines)
    diff -u <(sort "$ruby_p") <(sort "$rust_p") | head -50 > "$tmp/${f}.diff" || true
    echo "    (first 50 lines of diff saved to $tmp/${f}.diff)"
    status=1
  fi
done

echo "outputs kept in $tmp"
exit "$status"
