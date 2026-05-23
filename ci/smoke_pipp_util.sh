#!/usr/bin/env bash
# Hand-crafted smoke test for the pipp_util Rust binary.
#
# Writes synthetic versions of the three TSVs that PiPP produces
# (assign/per_query.tsv, feature/aa/feature.tsv,
# alignment/aligned_position.tsv) into a temp dir, runs
# `pipp_util import`, and verifies the resulting DuckDB file contains
# the expected row counts via the duckdb CLI.
#
# Requires:
#   - rust/target/release/pipp_util (run `cargo build --release` first)
#   - duckdb (CLI) on $PATH for verification
#
# Exits non-zero on any mismatch.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${PIPP_UTIL_BIN:-$ROOT/rust/target/release/pipp_util}"

if [[ ! -x "$BIN" ]]; then
  echo "pipp_util not found at $BIN. Build it first:" >&2
  echo "  cargo build --release --manifest-path rust/Cargo.toml" >&2
  exit 1
fi

if ! command -v duckdb >/dev/null 2>&1; then
  echo "duckdb CLI not found on \$PATH (needed to verify the imported DB)." >&2
  echo "  install: see https://duckdb.org/docs/installation/ or use the bundled conda env" >&2
  exit 1
fi

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

REFPKG="$TMP/result/synth_refpkg"
mkdir -p "$REFPKG/assign" "$REFPKG/feature/aa" "$REFPKG/alignment"

# --- assign/per_query.tsv : (gappa examine assign output) ---
cat >"$REFPKG/assign/per_query.tsv" <<'TSV'
name	LWR	fract	aLWR	afract	taxopath
seq1	0	0	1	1	Eukaryota
seq1	0.4	0.4	1	1	Eukaryota;Animalia
seq1	0.6	0.6	0.6	0.6	Eukaryota;Animalia;Chordata
seq2	0.5	0.5	1	1	Bacteria
TSV
EXP_ASSIGN=4

# --- feature/aa/feature.tsv : (01-6a.aa_feature output, 28 cols) ---
cat >"$REFPKG/feature/aa/feature.tsv" <<'TSV'
gene	len	len_of_std_aa	avg_MW	N-ARSC	C-ARSC	S-ARSC	K	R	H	D	E	N	Q	S	T	Y	A	V	L	I	P	F	M	W	G	C	others
seq1	100	100	120.5	0.5	0.3	0.1	5	5	2	8	8	5	3	7	6	3	7	7	8	6	5	4	2	1	5	3	0
seq2	150	148	118.2	0.4	0.35	0.12	8	7	3	10	12	7	4	10	9	5	10	10	11	8	7	6	3	2	8	5	2
TSV
EXP_AA=2

# --- alignment/aligned_position.tsv : (01-4h.aligned_position output, dynamic columns) ---
cat >"$REFPKG/alignment/aligned_position.tsv" <<'TSV'
query	TM2_D51	TM3_C78	TM3_E102	fract	taxpath
seq1	N	R	E	1.0	Eukaryota;Animalia
seq2	K	C	D	0.8	Bacteria
TSV
# 2 queries x 3 positions = 6 rows
EXP_POS=6

# --- run import ---
"$BIN" import "$REFPKG" --overwrite --refpkg synth_refpkg

DB="$REFPKG/pipp.duckdb"
[[ -f "$DB" ]] || { echo "DB not produced at $DB" >&2; exit 1; }

# --- verify row counts ---
got_assign=$(duckdb "$DB" -noheader -list -c "SELECT COUNT(*) FROM assignments")
got_aa=$(duckdb "$DB"     -noheader -list -c "SELECT COUNT(*) FROM aa_features")
got_pos=$(duckdb "$DB"    -noheader -list -c "SELECT COUNT(*) FROM aligned_positions")

fail=0
check() {
  local label="$1" got="$2" exp="$3"
  if [[ "$got" == "$exp" ]]; then
    printf "  %-18s %s rows  OK\n" "$label" "$got"
  else
    printf "  %-18s got=%s expected=%s  FAIL\n" "$label" "$got" "$exp" >&2
    fail=1
  fi
}
check assignments       "$got_assign" "$EXP_ASSIGN"
check aa_features       "$got_aa"     "$EXP_AA"
check aligned_positions "$got_pos"    "$EXP_POS"

# refpkg column sanity
got_refpkg=$(duckdb "$DB" -noheader -list -c "SELECT DISTINCT refpkg FROM assignments")
if [[ "$got_refpkg" != "synth_refpkg" ]]; then
  echo "  refpkg column   got=$got_refpkg expected=synth_refpkg  FAIL" >&2
  fail=1
fi

# aligned_positions schema sanity (long format)
got_pos_label=$(duckdb "$DB" -noheader -list -c "SELECT residues FROM aligned_positions WHERE query_name='seq1' AND pos_label='TM3_C78'")
if [[ "$got_pos_label" != "R" ]]; then
  echo "  long format     got=$got_pos_label expected=R (seq1 / TM3_C78)  FAIL" >&2
  fail=1
fi

exit "$fail"
