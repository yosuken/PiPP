#!/usr/bin/env bash
# Driver that runs ci/verify_parse_hmmsearch.sh across every PPP run
# under /aptmp/yosuke/www_db/uniparc/data/2024-02-06/PPP/*/rhodopsin_2020-11-17,
# in parallel, and summarises PASS / DIFF / FAIL counts.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
GLOB="${VERIFY_GLOB:-/aptmp/yosuke/www_db/uniparc/data/2024-02-06/PPP/*/rhodopsin_2020-11-17}"
JOBS="${JOBS:-4}"
LOGDIR="${VERIFY_LOG:-/tmp/pipp_verify_all}"

mkdir -p "$LOGDIR"
joblog="$LOGDIR/joblog"
results="$LOGDIR/results.tsv"
: > "$results"

runs=( $GLOB )
echo "verifying ${#runs[@]} runs with -j$JOBS"

# Per-run wrapper: run verification, capture status, write a one-line summary
# into $results. On PASS, remove the temp dir to keep disk small.
run_one() {
  local d="$1"
  local stamp
  stamp=$(date +%s)
  local log="$LOGDIR/$(basename "$(dirname "$d")")_${stamp}.log"
  if "$ROOT/ci/verify_parse_hmmsearch.sh" "$d" >"$log" 2>&1; then
    # tail finds the summary block to record byte-identical / sort-equal status
    local n_bi n_sort
    n_bi=$(grep -c BYTE-IDENTICAL "$log" || true)
    n_sort=$(grep -c 'IDENTICAL after sort' "$log" || true)
    echo -e "PASS\t$d\t bi=$n_bi sortequal=$n_sort" >> "$results"
    # Clean the temp tree (PASS case only — DIFF saves them in /tmp/pipp_verify)
    base=$(basename "$d")
    parent=$(basename "$(dirname "$d")")
    rm -rf "/tmp/pipp_verify/${parent}_"* 2>/dev/null || true
    rm -f "$log"
  else
    echo -e "FAIL\t$d" >> "$results"
    echo "  log: $log" >> "$results"
  fi
}
export -f run_one
export ROOT LOGDIR results

printf '%s\n' "${runs[@]}" \
  | parallel -j "$JOBS" --joblog "$joblog" --bar run_one {}

# Summarise
n_pass=$(grep -c '^PASS' "$results" || true)
n_fail=$(grep -c '^FAIL' "$results" || true)
echo
echo "============================ summary ============================"
echo "  total runs:     ${#runs[@]}"
echo "  PASS  (5/5 byte-identical): $n_pass"
echo "  FAIL  (some file differed): $n_fail"
echo
if [[ "$n_fail" -gt 0 ]]; then
  echo "FAILED runs:"
  grep '^FAIL' "$results" | sed 's/^/    /'
fi
echo "Results: $results"
echo "joblog:  $joblog"

exit "$n_fail"
