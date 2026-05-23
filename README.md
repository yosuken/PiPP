
# PiPP - a Pipeline for Phylogenetic Placement

[![license](https://img.shields.io/github/license/mashape/apistatus.svg)](/LICENSE)
[![size](https://img.shields.io/github/size/webcaetano/craft/build/phaser-craft.min.js.svg)]()

## currently PiPP is beta version. Any specification might be changed in a future version.
PiPP is developed as a tool for phylogenetic placement onto a clade or taxonomy defined phylogenetic tree through procedures below. Very large queries are acceptable.

## install (use conda environment)

### [1] make conda environment and install packages
```
$ conda create -n PiPP -y && conda activate PiPP
```

### [2] install packages
```
### create environment (name: PiPP_v0.3.0) and install packages
### change 'micromamba' to 'mamba' or 'conda' if you use mamba or conda instead of micromamba
$ v=PiPP_v0.3.0 && micromamba create -n $v -c conda-forge -c bioconda ruby=3.4.5 hmmer=3.4 parallel=20250822 gappa=0.9.0 pplacer=1.1.alpha19 mafft=7.520 fasttree=2.2.0 epa-ng=0.3.8 python=3.12 -y && micromamba activate $v && pip install apples taxtastic
```

### [3] place witch-ng binary (default aligner)

`witch-ng` is the default `--aligner` in PiPP but is not available on Bioconda / conda-forge, so you need to place the binary yourself. Pre-built binaries (Linux x86_64 / macOS) are published on the WITCH-NG releases page:

  https://github.com/RuneBlaze/WITCH-NG/releases

Download the binary for your platform, then put it under PiPP's `bin/` and mark it executable:

```
### example for Linux x86_64 (witch-ng v0.0.4)
$ curl -L -o bin/witch-ng https://github.com/RuneBlaze/WITCH-NG/releases/download/v0.0.4/witch-ng-x86_64-linux
$ chmod +x bin/witch-ng
```

If you'd rather not install witch-ng, pass `--aligner mafft-add` on every run.

## usage 
```
### PiPP ver 0.3.0 (2025-09-27) ###

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
$ PiPP [options] -q <query fasta(s)> -r <refpkg dir(s)> -o <output dir>

[dependencies]
- ruby (ver >= 2.0)
- hmmer (ver >= 3.0)
- mafft (tested by 7.453 and 7.520)
- gappa (tested by 0.6.0)
- pplacer (tested by 1.1.alpha19)
- witch-ng (ver >= 0.0.4)

[output files]
  result/<refpkg name>/{all,each}/{seq,alignment,placement,assign,extract,...} -- result of placement and further analysis

[options]
[File/directory]
    -q, --query FILE(S)              Query sequence file(s) (protein fasta, can be gzipped) [required]
    -r, --refpkg DIR(S)              Reference package(s) made by taxtastic [required]
    -o, --outdir PATH                Output directory [required]
        --[no-]overwrite             Overwrite output directory (default: overwrite)

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

## citation
```
```
