
# PiPP - a Pipeline for Phylogenetic Placement

[![CI](https://github.com/yosuken/PiPP/actions/workflows/ci.yml/badge.svg)](https://github.com/yosuken/PiPP/actions/workflows/ci.yml)

## currently PiPP is beta version. Any specification might be changed in a future version.
PiPP is developed as a tool for phylogenetic placement onto a clade or taxonomy defined phylogenetic tree through procedures below. A large query file is acceptable.

## install

### one-liner

```
$ ./install.sh                       # micromamba; use CONDA=mamba or CONDA=conda otherwise
$ micromamba activate PiPP_v0.4.0
```

`install.sh` wraps two steps. You can run them by hand if you prefer:

```
### [1] conda env (all bio tools + rust toolchain)
$ micromamba env create -f environment.yaml
$ micromamba activate PiPP_v0.4.0

### [2] build the bundled Rust binary `pipp_util`
$ cargo build --release --manifest-path rust/Cargo.toml
```

`pipp_util` ingests the per-refpkg TSV outputs into a DuckDB file at the end of the pipeline.

The pipeline locates the binary in this order:

1. `$PIPP_UTIL_BIN` if set and executable
2. `pipp_util` found on `$PATH` (e.g. via Bioconda)
3. `rust/target/release/pipp_util` (after the `cargo build` above)

### (future) Bioconda

A `meta.yaml` recipe is included in the repo as a skeleton for a future Bioconda package. The plan is for `pipp` to install both the Ruby orchestration and the Rust `pipp_util` binary in one shot.

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
$ PiPP [options] -q <query fasta> -r <refpkg dir(s)> -o <output dir>

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

[Placement]
        --placer OPTION              query sequence aligner (default: pplacer) [pplacer|apples-2|epa-ng]
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

| table              | source TSV                              | grain                                |
|--------------------|------------------------------------------|--------------------------------------|
| `assignments`      | `assign/per_query.tsv`                   | one row per (query, taxopath)        |
| `aa_features`      | `feature/aa/feature.tsv`                 | one row per query                    |
| `aligned_positions`| `alignment/aligned_position.tsv`         | one row per (query, position_label)  |

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
