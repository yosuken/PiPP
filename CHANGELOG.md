# Changelog

All notable changes to PiPP are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project aims to adhere to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).
Pre-1.0 minor bumps may still include breaking changes.

## Release process

When cutting a release `vX.Y.Z`:

1. Rename `## [Unreleased]` below to `## [X.Y.Z] - YYYY-MM-DD` and add a
   fresh empty `## [Unreleased]` block above it.
2. Bump the version in `pipp` (Ruby `VERSION`/`MODIDATE`), `rust/Cargo.toml`,
   `environment.yaml` (`name: PiPP_vX.Y.Z`), `meta.yaml`
   (`{% set version = "X.Y.Z" %}`), and README's `### PiPP ver X.Y.Z` block.
3. Update the link references at the bottom of this file
   (`[Unreleased]`, `[X.Y.Z]`).
4. Commit, then tag and push:

   ```
   git tag -a vX.Y.Z -m "Release X.Y.Z"
   git push origin main vX.Y.Z
   ```

   The release workflow (`.github/workflows/release.yml`) builds prebuilt
   `pipp_util` binaries for Linux x86_64 and macOS x86_64/aarch64 and
   attaches them to the GitHub Release.

## [Unreleased]

### Added
- `pipp_util` Rust binary that ingests the per-refpkg TSV outputs into a
  DuckDB file (`result/<refpkg>/pipp.duckdb`) at the end of a run. Three
  tables: `assignments`, `aa_features`, `aligned_positions` (long format).
- New Rake task `01-7a.duckdb_import` (runs as the final pipeline step).
- `--keep-intermediate` flag on `PiPP` to retain `prefilter/`, `chunks/`,
  `batch/`, and `log/tasks/` after a run.
- pixi as the primary environment manager: `pixi.toml` + `pixi.lock` pin the
  conda runtime (bio tools) and PyPI deps (apples, taxtastic) in one lockfile;
  a separate `build` environment isolates the Rust toolchain that compiles
  `pipp_util` (`pixi run -e build build`). `install.sh` now wraps
  `pixi install` + the build task.
- `environment.yaml` for conda environment setup (`PiPP_v0.4.0`), kept as a
  compat path for conda/micromamba users (pixi is preferred).
- `meta.yaml` skeleton for a future Bioconda recipe.
- `ci/run.sh` (local-runnable) and `.github/workflows/ci.yml` covering
  Ruby syntax, Rust fmt/clippy/build, and a self-contained `pipp_util`
  smoke test using synthetic TSVs.
- The entry point and Rake file are now lowercase (`pipp`, `pipp.rake`); `pipp`
  is the canonical command. `PiPP` is kept as a backward-compat alias created at
  install time (PATH launcher, or a guarded `$PREFIX/bin` symlink for Bioconda),
  not committed to the repo — `PiPP` and `pipp` case-fold to the same name and
  would collide on case-insensitive filesystems (macOS/Windows), so only the
  lowercase files live in git. The entry point resolves symlinks
  (`File.realpath`) so invocation through the `PiPP` alias still finds
  `pipp.rake`. Project branding remains "PiPP".

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
- Ruby constraint relaxed to `>=3` (was `>=3,<4`); validated on 3.4 and 4.0.
  The pixi.lock stays on the battle-tested 3.4.x until `pixi update ruby`.

### Removed
- `environment.lock.txt` and `requirements.lock.txt` (manual conda/pip locks),
  superseded by `pixi.lock` which pins conda + PyPI in one file.
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

### Release-time TODO (drop when merging this branch to main)
- README's CI badge currently pins `?branch=feature/duckdb`. Once this
  branch lands on `main`, drop the query string so the badge tracks the
  default branch.

## [0.3.0] - 2025-09-27

First tagged release. Single-FASTA / multi-FASTA support, witch-ng,
pplacer/apples-2/epa-ng placers, gappa-based chunkify pipeline.

[Unreleased]: https://github.com/yosuken/PiPP/compare/v0.3.0...HEAD
[0.3.0]: https://github.com/yosuken/PiPP/releases/tag/v0.3.0
