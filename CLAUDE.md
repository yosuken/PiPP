# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

PiPP (Pipeline for Phylogenetic Placement) is a Ruby-based bioinformatics tool for phylogenetic placement of query sequences onto reference phylogenetic trees. The pipeline accepts a single (large) query FASTA file and runs it through a multi-step workflow involving sequence prefiltering, alignment, and phylogenetic placement analysis.

## Build and Test Commands

### Running Tests
```bash
rake test
# Or using the default task
rake
```

### Running the Main Pipeline
```bash
# Basic usage
./PiPP -q <query_fasta> -r <refpkg_dir> -o <output_dir>

# Example with options
./PiPP -q queries/sample.fa -r "refpkgs/*" -o results --ncpus 4 --evalue 1e-10
```

### Development Commands
The project uses Rake as its build system. Key Rake tasks are defined in `PiPP.rake`.

## Architecture Overview

### Core Components

1. **Main Pipeline Script (`PiPP`)**: Bash wrapper that validates dependencies, parses command-line arguments, and invokes the Rake-based workflow

2. **Rake Workflow (`PiPP.rake`)**: Orchestrates the entire pipeline through sequential tasks:
   - Query validation and preprocessing
   - Reference package validation  
   - HMM-based sequence prefiltering (hmmsearch)
   - Sequence alignment (MAFFT)
   - Phylogenetic placement (pplacer with gappa chunking)
   - Analysis and visualization

3. **Ruby Processing Scripts (`script/`)**: Specialized utilities for:
   - Sequence validation and trimming
   - HMM search result parsing
   - Alignment processing and format conversion
   - Feature extraction

### Pipeline Workflow

The pipeline follows a numbered task sequence (01-1a through 01-6a):
- **01-1x**: Query validation and preprocessing
- **01-2x**: HMM-based prefiltering using hmmsearch
- **01-3x**: Chunked alignment and phylogenetic placement
- **01-4x**: Placement analysis (assignment, grafting, visualization)
- **01-5x**: Comparative analysis (KRD, EdgePCA)
- **01-6x**: Feature extraction

### Key Dependencies

External tools required:
- hmmer (≥3.0) for sequence similarity detection
- mafft (7.453+) for sequence alignment  
- gappa (0.6.0+) for placement processing and chunking
- pplacer (1.1.alpha19+) for phylogenetic placement
- GNU parallel for parallelization

### Directory Structure

- `script/`: Ruby processing utilities
- `test/`: Minitest-based test suite
- Output structure: `result/<refpkg>/{seq,alignment,placement,assign,graft,feature,...}/`
- Intermediate dirs (`prefilter/`, `chunks/`, `batch/`, `log/tasks/`) are removed at end of run unless `--keep-intermediate` is passed

### Configuration

The pipeline accepts extensive command-line configuration including:
- E-value thresholds for prefiltering
- Alignment methods (FFT-NS-2, FFT-NS-i, E-INS-i)  
- Chunk sizes for parallelization
- Extraction levels for taxonomic analysis

## Development Notes

- Uses Ruby 2.7+ with Rake for task orchestration
- Implements custom Range extensions for overlap detection
- Batch job generation for parallel execution
- Comprehensive logging and error handling throughout pipeline