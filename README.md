
# PiPP - a Pipeline for Phylogenetic Placement

[![CI](https://github.com/yosuken/PiPP/actions/workflows/ci.yml/badge.svg?branch=feature/duckdb)](https://github.com/yosuken/PiPP/actions/workflows/ci.yml?query=branch%3Afeature%2Fduckdb)

## currently PiPP is beta version. Any specification might be changed in a future version.
PiPP is developed as a tool for phylogenetic placement onto a clade or taxonomy defined phylogenetic tree through procedures below. A large query file is acceptable.

## install

PiPP uses [pixi](https://pixi.sh) as its primary environment manager: one
`pixi.toml` / `pixi.lock` pins both the conda bio tools and the PyPI deps
(apples, taxtastic), and a separate `build` environment holds the Rust
toolchain that compiles `pipp_util`.

### 1. one-liner (pixi)

```
$ ./install.sh                       # pixi install + build pipp_util
```

`install.sh` wraps two pixi steps. You can run them by hand if you prefer:

```
### [1] runtime env (all bio tools + apples/taxtastic), pinned by pixi.lock
$ pixi install

### [2] build the bundled Rust binary `pipp_util` in the isolated build env
$ pixi run -e build build            # = cargo build --release (rust-only env)
```

Then run `pipp` through pixi (it activates the env and keeps your cwd, so
relative `-q`/`-r`/`-o` paths resolve as expected):

```
$ pixi run ./pipp -q <query.fa> -r <refpkg> -o <out>
```

The canonical command is `pipp` (lowercase); `PiPP` is kept as a backward-compatible
alias. To call them from anywhere, drop a launcher on your `PATH`:

```
$ cat > ~/.local/bin/pipp <<'EOF'
#!/usr/bin/env bash
exec pixi run --manifest-path /ABS/PATH/TO/PiPP/pixi.toml /ABS/PATH/TO/PiPP/pipp "$@"
EOF
$ chmod +x ~/.local/bin/pipp
$ ln -sf pipp ~/.local/bin/PiPP        # compatible alias (skip on case-insensitive FS)
```

> `pipp` and `PiPP` case-fold to the same name, so only the lowercase `pipp`
> lives in the repo; the uppercase alias is created here at install time (and
> guarded in the Bioconda recipe) to avoid breaking clones on case-insensitive
> filesystems (macOS/Windows).

`pipp_util` ingests the per-refpkg TSV outputs into a DuckDB file at the end of the pipeline.

> **Note on `rust` and `duckdb`:** Rust is build-time only — it lives in pixi's
> `build` environment (never co-solved with the runtime deps) and just compiles
> `pipp_util`. `duckdb` is not a runtime dependency at all — `pipp_util` bundles
> its own DuckDB. A `duckdb` CLI (>=1.0) is optional, only for querying the
> output DBs (use any system install / mise / conda).

#### witch-ng binary (default `--aligner`)

`witch-ng` (GPLv3) is the default aligner but is not packaged on Bioconda or conda-forge. A helper fetches the right prebuilt binary into `bin/`:

```
$ bin/install_witch_ng.sh                       # detects platform, fetches v0.0.4
$ WITCH_NG_VERSION=v0.0.4 bin/install_witch_ng.sh
```

If you'd rather skip witch-ng, pass `--aligner mafft-add` on every run.


### 2. conda / micromamba

If you'd rather not use pixi, `environment.yaml` is kept for conda-family tools.
A strict solve fails (pplacer/fasttree deps vs latest conda-forge), so flexible
priority is required, and you build `pipp_util` with a system `cargo`:

```
$ micromamba env create -f environment.yaml --channel-priority flexible
$ micromamba activate PiPP_v0.4.0
$ cargo build --release --manifest-path rust/Cargo.toml
```

The pipeline locates the binary in this order:

1. `$PIPP_UTIL_BIN` if set and executable
2. `pipp_util` found on `$PATH` (e.g. via Bioconda)
3. `rust/target/release/pipp_util` (after the `cargo build` above)

#### witch-ng binary (default `--aligner`)

See above.

#### apples (only used with `--placer apples-2`)

The default placer is `pplacer`; `apples` (GPL-3.0, PyPI-only) is used only when
you pass `--placer apples-2`. No separate `pip install apples` is needed — it is
already bundled in **both** install paths (`pixi install` via `pixi.toml`, and
`environment.yaml` for the conda path). If `apples` is missing the run still
works with the default placer; PiPP's startup tool check just notes it as
`(not found)`.

### (future) Bioconda

A `meta.yaml` recipe is included in the repo as a skeleton for a future Bioconda package. The plan is for `pipp` to install the Ruby orchestration and the Rust `pipp_util` binary in one shot. witch-ng / apples will likely remain user-installed (see notes above).

## usage 
```
### PiPP ver 0.4.0 (2026-05-23) ###

PiPP - Pipeline for phylogenetic placement.
PiPP is developed as a tool for phylogenetic placement onto a clade or taxonomy defined phylogenetic tree through procedures below.

1. prefilter query sequences by similarity detection (hmmsearch). Queries and references should be protein sequences at this moment.
2. align query sequences to a given reference alignment ('witch-ng', 'mafft --add', or 'mafft --addfragments')
3. perform phylogenetic placement ('pplacer', 'apples-2', or 'epa-ng') with efficient parallelization using 'gappa prepare chunkify' and 'gappa prepare unchunkify'
4. analysis of placed sequences
  a. assign clade/taxonomy and generate statistics ('gappa examine assign')
  b. extract placed sequences and placement file for each clade/taxonomy ('gappa prepare extract')
  c. extract placed sequences and placement file for each clade/taxonomy ('gappa prepare extract')

[usage]
$ pipp [options] -q <query fasta> -r <refpkg dir(s)> -o <output dir>

[dependencies]
- ruby (ver >= 2.0)
- hmmer (ver >= 3.0)
- mafft (tested by 7.453 and 7.520)
- gappa (tested by 0.6.0)
- pplacer (tested by 1.1.alpha19)
- witch-ng (ver >= 0.0.4)

[output files]
  result/<refpkg name>/{seq,alignment,placement,assign,graft,feature,...} -- result of placement and further analysis
  result/<refpkg name>/pipp.duckdb -- DuckDB file with assignments / aa_features / aligned_positions tables (produced by pipp_util)

[options]
[File/directory]
    -q, --query FILE                 Query sequence file (protein fasta, can be gzipped) [required]
    -r, --refpkg DIR(S)              Reference package(s) made by taxtastic [required]
    -o, --outdir PATH                Output directory [required]
        --[no-]overwrite             Overwrite output directory (default: overwrite)
        --keep-intermediate          Keep intermediate files (prefilter/, chunks/, batch/, log/tasks/) after run (default: removed)
        --copy-refpkg                Copy refpkg derived files (backbone/trees/eHMM) into the run's refpkg/ dir (default: reference the <refpkg>/derived cache in place)

[Task]
        --only-detect                Only detect homologous regions of input sequences using hmmsearch

[Prefilter (result cutoffs)]
    -e, --evalue NUM                 E-value threshold of hmmsearch (default: 1e-5)
        --minseqlen INT              set a cutoff of minimum amino acid length of input sequences (default: 0)
        --minhmmlen INT              Minimum hmm hit length in linked result of hmmsearch (default: 0)
        --minhmmcov FLOAT            Minimum fraction of hmm length in linked result of hmmsearch (default: 0)
        --minalilen INT              Minimum hmm hit length in linked result of hmmsearch (default: 0)
        --minalicov FLOAT            Minimum fraction of hmm length in linked result of hmmsearch (default: 0)

[Prefilter (domain-level cutoffs)]
        --evaluedom NUM              Domain E-value threshold of hmmsearch (default: 1e-2)
        --minhmmlendom INT           Minimum hmm hit length in domain-level result of hmmsearch (default: 0)
        --minhmmcovdom FLOAT         Minimum fraction of hmm length in domain-level result of hmmsearch (default: 0)
        --minalilendom INT           Minimum hmm hit length in domain-level result of hmmsearch (default: 0)
        --minalicovdom FLOAT         Minimum fraction of hmm length in domain-level result of hmmsearch (default: 0)

[Alignment]
        --aligner OPTION             query sequence aligner (default: witch-ng)
        --mafft-method METHOD        MAFFT add method (default: E-INS-i)
        --witch-ng-hmm-size-lb INT   witch-ng eHMM decomposition size lower bound (default: backbone leaf count / 20)

[Placement]
        --placer OPTION              phylogenetic placement tool (default: pplacer) [pplacer|apples-2|epa-ng]
        --epa-ng-model MODEL         model for epa-ng, either model name (e.g., LG, PROTGTR, ...) or tree log file (compatible with RAxML 8.x and IQ-TREE)
                                     [required when '--placer epa-ng' is selected]
                                     Please refer to epa-ng document. [https://github.com/pierrebarbera/epa-ng?tab=readme-ov-file#setting-the-model-parameters]

[Computation]
    -c, --chunk-size INT             Chunk size of 'gappa prepare chunkify' (default: 10000)
    -n, --ncpus INT                  Number of CPUs to use (default: 1)

[General]
    -h, --help                       Show this help message
    -v, --version                    Show version
```

## DuckDB output

At the end of a run, `pipp_util import` loads three TSVs from each `result/<refpkg>/` into `result/<refpkg>/pipp.duckdb`:

| table              | source                                   | grain                                |
|--------------------|------------------------------------------|--------------------------------------|
| `assignments`      | `assign/per_query.tsv`                   | one row per (query, taxopath)        |
| `aa_features`      | `feature/aa/feature.tsv`                 | one row per query                    |
| `aligned_positions`| `alignment/aligned_position.tsv`         | one row per (query, position_label); `pos_index` preserves the source TSV column order so the wide layout can be reconstructed |
| `jplace_clamps`    | `placement/*.clamp.tsv`                  | one row per sanitized jplace (only when clamp-jplace fixed non-finite/negative-branch values) |
| `refpkgs`          | `refpkg/<name>/backbone.json`            | one row per refpkg (identity + source provenance: refpkg_dir, hmmname, hmmlen, aln/tree/hmm source paths) |
| `run_params`       | `run_params.json`                        | one row per (refpkg, option) — all run-time options as key/value, plus pipp_version, run_datetime, command_line, query |
| `software`         | `run_params.json`                        | one row per external tool (name, resolved path, version) |
| `prefilter_hits`   | `prefilter/.../best-hit.tsv`             | one row per detected region for this refpkg's HMM (hmmsearch i-Evalue, score, hmm/ali/env coords, …) |
| `prefilter_evalues`| `prefilter/.../evalues.tsv`              | one row per (seq) with a non-empty best i-Evalue against this refpkg's HMM |
| `query_whole`      | `seq/whole.fa`                           | one row per detected protein — whole (ungapped) sequence |
| `query_aligned`    | `alignment/aligned_wo_ref.fa`            | one row per region — aligned (gapped) query in backbone columns |

The other sequence files are **derived**, not stored, to avoid redundancy:

- `query_region` (a view) reconstructs `seq/region.fa` = `query_whole` sliced at each hit's `prefilter_hits.ali_fm..ali_to`.
- the full `alignment/aligned.fa` = the refpkg backbone alignment (`<refpkg>/derived/backbone.mfa`) + `query_aligned`.

`aligned_positions` is stored in long format: the dynamic per-position columns from the TSV (driven by the refpkg's `position.tsv`) become rows keyed by `pos_label`.

Each table has a `refpkg` column so multiple DBs can be `ATTACH`-ed and unioned. Quick examples:

```
$ duckdb result/opsin/pipp.duckdb -c "SELECT taxopath, COUNT(*) FROM assignments GROUP BY 1 ORDER BY 2 DESC LIMIT 10"
$ duckdb result/opsin/pipp.duckdb -c "SELECT * FROM aligned_positions WHERE pos_label='Lys296' AND residues LIKE 'K%'"
```

Standalone usage:

```
$ pipp_util import result/opsin/         # writes result/opsin/pipp.duckdb
$ pipp_util import result/opsin/ --db custom.duckdb --refpkg opsin --overwrite
```

## development / CI

The same checks run on GitHub Actions can be run locally:

```
$ ci/run.sh             # ruby syntax + rust fmt/clippy/build + pipp_util smoke
$ ci/run.sh ruby        # only ruby
$ ci/run.sh rust        # only rust
$ ci/run.sh smoke       # only the pipp_util import smoke test
```

The smoke test (`ci/smoke_pipp_util.sh`) writes synthetic TSVs to a temp dir, runs `pipp_util import`, and verifies row counts via the `duckdb` CLI. It does **not** require any of the heavy bio tools (hmmsearch, mafft, gappa, pplacer, witch-ng).

## migration note (v0.3.x → v0.4.0)

- `-q` now accepts a single query FASTA only. Glob patterns and comma-separated lists are no longer supported.
- Output directory layout has been flattened: `result/<refpkg>/<task>/` (no more `all/`, `each/`, or per-query subdirectories).
- Intermediate directories (`prefilter/`, `chunks/`, `batch/`, `log/tasks/`) are removed by default after a successful run. Pass `--keep-intermediate` to retain them (useful for debugging or resuming a failed run).
- New `01-7a.duckdb_import` task produces `result/<refpkg>/pipp.duckdb` from the per-task TSVs. Requires the bundled Rust binary `pipp_util` (see install section).

## citation
```
```
