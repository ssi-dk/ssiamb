# ssiamb — Ambiguous Sites Counter

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.12+](https://img.shields.io/badge/python-3.12+-blue.svg)](https://www.python.org/downloads/)

**Author:** Povilas Matusevicius <pmat@ssi.dk>  
**Repository:** [https://github.com/ssi-dk/ssiamb](https://github.com/ssi-dk/ssiamb)  
**License:** MIT  
**Minimum Python:** 3.12+

## Overview

`ssiamb` computes an "ambiguous sites" metric for bacterial whole genome sequencing (WGS) as a measure of within-sample heterogeneity. This tool modernizes and standardizes the lab's prior definition while providing robust packaging, CLI interface, and Galaxy integration capabilities.

## What are "Ambiguous Sites"?

An ambiguous site is a genomic position with:
 
- **Sufficient coverage**: Depth ≥ `dp_min` (default: 10)
- **Minor-allele signal**: Minor-allele fraction (MAF) ≥ `maf_min` (default: 0.10)

These metrics are determined from variant calls after normalization and atomization, counting **once per locus** (multi-allelic sites count once if any ALT passes the thresholds).

## Supported Modes

### Self-mapping Mode (`ssiamb self`)

- **Input**: Reads → Sample's own assembly
- **Use case**: Analyze heterogeneity against the sample's assembled genome
- **Mapping space**: Uses the assembly as reference

### Reference-mapped Mode (`ssiamb ref`)

- **Input**: Reads → Species canonical reference
- **Use case**: Compare against standardized reference genomes
- **Reference selection**: Via admin directory or user override

## Key Features

- **Flexible mapping**: Support for minimap2 (default) and bwa-mem2
- **Multiple variant callers**: BBTools (default) and bcftools
- **Comprehensive outputs**: Summary TSV (always) + optional VCF, BED, matrices, per-contig analysis
- **Depth analysis**: Using mosdepth (default) or samtools
- **Reusable workflows**: Accept pre-computed BAM/VCF files
- **Galaxy integration**: Designed for workflow environments
- **Quality control**: Configurable thresholds with sensible defaults

## Installation

This project is under active development. You can install the package for development (editable) from the repository so changes to the source are immediately available:

```bash
# Install in editable/development mode (recommended for contributors)
pip install -e .

# After editable install you can run the CLI via the console script or module:
ssiamb --help
# or
python -m ssiamb --help
```

 
When a stable release is published the package will also be available via PyPI and Bioconda (example future commands):

```bash
# Future installation via pip (PyPI)
pip install ssiamb

# Future installation via conda (Bioconda)
conda install -c bioconda ssiamb
```

## Quick Start

```bash
# Self-mapping mode: analyze reads against sample's own assembly
ssiamb self --r1 sample_R1.fastq.gz --r2 sample_R2.fastq.gz --assembly sample.fna

# Reference-mapped mode: analyze against species reference
ssiamb ref --r1 sample_R1.fastq.gz --r2 sample_R2.fastq.gz --species "Escherichia coli"

# With custom thresholds
ssiamb self --r1 reads_R1.fastq.gz --r2 reads_R2.fastq.gz --assembly assembly.fna \
  --dp-min 15 --maf-min 0.05

# Output to stdout (no files written)
ssiamb self --r1 reads_R1.fastq.gz --r2 reads_R2.fastq.gz --assembly assembly.fna --stdout
```

## Output

### Primary Output

- **`ambiguous_summary.tsv`**: Single-row summary with ambiguous site counts and quality metrics

### Optional Outputs (via flags)

- **`--emit-vcf`**: Variant calls with ambiguity annotations
- **`--emit-bed`**: BED file of ambiguous sites
- **`--emit-matrix`**: Depth×MAF cumulative count matrix
- **`--emit-per-contig`**: Per-contig breakdown
- **`--emit-provenance`**: Analysis provenance and parameters
- **`--emit-multiqc`**: MultiQC-compatible reports

## Dependencies

### Required External Tools

- **Mapping**: `minimap2` or `bwa-mem2`
- **Variant calling**: `BBTools` (callvariants.sh) or `bcftools`
- **Depth analysis**: `mosdepth` or `samtools`
- **VCF processing**: `bcftools` (for normalization)

### Python Dependencies

- Python 3.12+
- `typer[all]`, `rich`, `pyyaml`, `pandas`, `numpy`, `pysam`, `biopython`

 
Install Python deps in a virtual environment (example):

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e .
```

### Running tests

The project includes a pytest test-suite. Run all tests with:

```bash
python -m pytest tests/ -v
```

Contributors should install in editable mode and run the tests before opening PRs.

Test dependencies:

- numpy
- pysam
- biopython


## Development Status

This project is currently in active development. The implementation follows a structured approach:

1. ✅ **Planning & Specification** - Comprehensive requirements defined
2. 🚧 **Repository Bootstrap** - Setting up package structure
3. ⏳ **Core Implementation** - CLI, models, and processing pipelines
4. ⏳ **External Tool Integration** - Mapping and variant calling
5. ⏳ **Testing & Validation** - Unit tests and integration testing
6. ⏳ **Packaging & Distribution** - Bioconda, containers, Galaxy tools

 
## Contributing

This project is developed by the SSI team. For questions or contributions, please contact:

- Povilas Matusevicius <pmat@ssi.dk>

## Citation

> **Note**: Citation information will be provided upon publication.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
