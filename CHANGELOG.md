# Changelog

All notable changes to PiPP are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project aims to adhere to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).
Pre-1.0 minor bumps may still include breaking changes.

## [Unreleased]

### Added
- `pipp_util` Rust binary that ingests the per-refpkg TSV outputs into a
  DuckDB file (`result/<refpkg>/pipp.duckdb`) at the end of a run. Three
  tables: `assignments`, `aa_features`, `aligned_positions` (long format).
- New Rake task `01-7a.duckdb_import` (runs as the final pipeline step).
- `--keep-intermediate` flag on `PiPP` to retain `prefilter/`, `chunks/`,
  `batch/`, and `log/tasks/` after a run.
- `environment.yaml` for conda environment setup (`PiPP_v0.4.0`).
- `meta.yaml` skeleton for a future Bioconda recipe.
- `ci/run.sh` (local-runnable) and `.github/workflows/ci.yml` covering
  Ruby syntax, Rust fmt/clippy/build, and a self-contained `pipp_util`
  smoke test using synthetic TSVs.
- `install.sh` rewritten as a thin wrapper around `environment.yaml` +
  `cargo build`.

### Changed
- `-q` now accepts a single FASTA only. Glob patterns and
  comma-separated lists are no longer supported.
- Output layout flattened to `result/<refpkg>/<task>/`. The `all/`,
  `each/`, and per-query subdirectories are gone; file names within each
  task directory are fixed.
- Intermediate directories (`prefilter/`, `chunks/`, `batch/`,
  `log/tasks/`) are removed by default after a successful run. Pass
  `--keep-intermediate` to retain them.
- `.gitignore` narrowed: only `test/example1/` is excluded under `test/`
  (large fixtures); test code under `test/1/` is no longer ignored. Added
  the usual editor/OS/Python/scratch patterns.
- README bash status badges replaced with a real CI badge.

### Removed
- Obsolete `PiPP.sh` bash wrapper (the Ruby `PiPP` has been the entry
  point for a while).
- The `merge_jplace.rb` invocation in `01-3d.unchunkify`. With a single
  query there is nothing to merge; `gappa unchunkify` now writes the
  final `.jplace` directly into `result/<refpkg>/placement/`.

### Migration notes (v0.3.x → 0.4.0)
- Downstream tooling that read `result/<refpkg>/all/...` or
  `result/<refpkg>/each/<query>/...` paths must be updated. The plain
  `result/<refpkg>/<task>/` layout is the only supported output.
- Mid-pipeline failures cannot be resumed in place when the default
  cleanup runs. Pass `--keep-intermediate` if you might need to resume.

## [0.3.0] - 2025-09-27

First tagged release. Single-FASTA / multi-FASTA support, witch-ng,
pplacer/apples-2/epa-ng placers, gappa-based chunkify pipeline.

[Unreleased]: https://github.com/yosuken/PiPP/compare/v0.3.0...HEAD
[0.3.0]: https://github.com/yosuken/PiPP/releases/tag/v0.3.0
