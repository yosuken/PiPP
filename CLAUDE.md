# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

PiPP (Pipeline for Phylogenetic Placement) is a bioinformatics tool for phylogenetic placement of a single query protein FASTA onto a clade- or taxonomy-defined reference phylogenetic tree. It is implemented as a Ruby + Rake orchestration layer plus a bundled Rust binary (`pipp_util`).

## Build and Test Commands

### One-shot CI (same as GitHub Actions)
```bash
./ci/run.sh             # ruby -wc + rust fmt/clippy/build/test + pipp_util smoke
./ci/run.sh ruby        # only ruby
./ci/run.sh rust        # only rust (fmt, clippy, build, test)
./ci/run.sh smoke       # only the pipp_util import smoke test
```

There is no `rake test` target. The whole `test/` directory is gitignored (it holds large local fixtures plus unimplemented Minitest scaffolds); automated tests live in `cargo test` and `ci/`.

### Building the bundled Rust binary
```bash
cargo build --release --manifest-path rust/Cargo.toml
```

The Rake pipeline locates `pipp_util` in this order: `$PIPP_UTIL_BIN` → `which pipp_util` → `rust/target/release/pipp_util`.

### Running the Main Pipeline
```bash
# Basic usage (single query FASTA only since v0.4.0)
./pipp -q <query_fasta> -r <refpkg_dir> -o <output_dir>

# Example with options
./pipp -q queries/sample.fa -r refpkgs/refpkg1 -o results --ncpus 4 --evalue 1e-10 --keep-intermediate
```

## Architecture Overview

### Core Components

1. **`pipp` (Ruby entry point)**: parses CLI options, sets ENV, hands off to Rake. (`PiPP` is a backward-compat alias.)
2. **`pipp.rake`**: orchestrates the sequential task graph.
3. **`script/` (Ruby helpers)**: small per-task scripts (validation, alignment post-processing, feature extraction).
4. **`rust/` (`pipp_util` crate)**: bundled Rust binary providing the `import` (TSV → DuckDB) and `parse-hmmsearch` subcommands. Hot-path scripts are being ported here over time.

### Pipeline Workflow

Numbered task sequence (01-1a through 01-7a):
- **01-1x**: Query and refpkg validation
- **01-2x**: HMM-based prefiltering (`hmmsearch` + `pipp_util parse-hmmsearch`)
- **01-3x**: Chunked alignment (witch-ng / mafft-add) and placement (pplacer / apples-2 / epa-ng)
- **01-4x**: Placement analysis (info, lwr-list, assign, graft, aligned_position)
- **01-6x**: AA feature extraction
- **01-7x**: DuckDB import (`pipp_util import`)

### Key Dependencies

External tools required:
- hmmer (≥3.0) for sequence similarity detection
- mafft (7.453+) for sequence alignment
- gappa (0.6.0+) for placement processing and chunking
- pplacer (1.1.alpha19+) for phylogenetic placement
- witch-ng (default `--aligner`; not on conda, install manually — see README)
- GNU parallel for parallelization
- Rust toolchain (for building the bundled `pipp_util` binary)

### Directory Structure

- `script/`: Ruby processing utilities (per Rake task)
- `rust/`: `pipp_util` crate (subcommands: `import`, `parse-hmmsearch`)
- `ci/`: local-runnable CI (`run.sh`, `smoke_pipp_util.sh`)
- `.github/workflows/`: GitHub Actions
- `bin/`: helper scripts and the manually-placed `witch-ng` binary
- Output structure: `result/<refpkg>/{seq,alignment,placement,assign,graft,feature,...}/`
- DuckDB output: `result/<refpkg>/pipp.duckdb` (tables: `assignments`, `aa_features`, `aligned_positions`)
- Intermediate dirs (`prefilter/`, `chunks/`, `batch/`, `log/tasks/`) are removed at end of run unless `--keep-intermediate` is passed

### Configuration

The pipeline accepts extensive command-line configuration including:
- E-value thresholds for prefiltering (`-e`, `--evaluedom`)
- Aligner choice (`witch-ng` (default) or `mafft-add`)
- witch-ng eHMM decomposition size (`--witch-ng-hmm-size-lb`, default: backbone leaf count / 20)
- Placer choice (`pplacer` (default), `apples-2`, `epa-ng`)
- Chunk sizes for parallelization (`-c`)
- Refpkg derived-file handling (`--copy-refpkg`): by default a valid `<refpkg>/derived/` cache (backbone.mfa, FastTree min-evo/gamma trees, taxit package, witch-ng eHMMs) is referenced in place so the run's `refpkg/` dir stays tiny; `--copy-refpkg` materializes a full isolated copy under `refpkg/<name>/` (tens of MB per refpkg, mostly eHMMs)

## Development Notes

- Ruby (≥2.7) + Rake for orchestration; Rust (≥1.70, edition 2021) for the bundled binary.
- `ruby -wc` is used in CI for syntax + warnings (no rubocop/standardrb — explicit choice to avoid style wars).
- `pipp_util` ships its own DuckDB (via the duckdb-rs `bundled` feature); CI installs the duckdb CLI separately for verification.
- Hot-path Ruby scripts are being ported to Rust one subcommand at a time (parse-hmmsearch already; others candidates: 01-3e.unchunkify_alignment, 01-6a.aa_feature).
