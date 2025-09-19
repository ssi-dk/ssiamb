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
# Check what would be done (dry run)
ssiamb --dry-run self --r1 sample_R1.fastq.gz --r2 sample_R2.fastq.gz --assembly sample.fna

# Self-mapping mode: analyze reads against sample's own assembly
ssiamb self --r1 sample_R1.fastq.gz --r2 sample_R2.fastq.gz --assembly sample.fna

# Reference-mapped mode: analyze against species reference
ssiamb ref --r1 sample_R1.fastq.gz --r2 sample_R2.fastq.gz --species "Escherichia coli"

# Summarize existing VCF and BAM files
ssiamb summarize --vcf sample.vcf.gz --bam sample.bam

# With custom thresholds and optional outputs
ssiamb self --r1 reads_R1.fastq.gz --r2 reads_R2.fastq.gz --assembly assembly.fna \
  --dp-min 15 --maf-min 0.05 --emit-vcf --emit-bed

# Output to stdout (no files written)
ssiamb self --r1 reads_R1.fastq.gz --r2 reads_R2.fastq.gz --assembly assembly.fna --stdout
```

## Error Codes

`ssiamb` follows a structured exit code system for programmatic handling:

- **0**: Success
- **1**: CLI/input errors (missing files, invalid sample names, bad arguments)
- **2**: Reference mode selection errors (species not found, Bracken failures)
- **3**: Reuse compatibility errors (VCF/BAM mismatch with reference)
- **4**: External tool failures (missing tools, tool execution errors)
- **5**: QC failures (only when `--qc-action fail` is enabled)

Errors include helpful suggestions and available options when applicable.

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

## Release Process

This project uses automated publishing to PyPI, Bioconda, and Galaxy ToolShed. The release process is as follows:

### 1. Version Update

1. Update version in `pyproject.toml`:

   ```toml
   [project]
   version = "1.0.0"  # Update this
   ```

2. Update version in `recipes/ssiamb/meta.yaml`:

   ```yaml
   {% set version = "1.0.0" %}  # Update this
   ```

3. Update version in `galaxy/ssiamb.xml`:

   ```xml
   <tool id="ssiamb" name="Ambiguous Sites Counter" version="1.0.0+galaxy0">
   ```

### 2. Create Release

1. Commit version changes:

   ```bash
   git add pyproject.toml recipes/ssiamb/meta.yaml galaxy/ssiamb.xml
   git commit -m "Bump version to v1.0.0"
   git push origin main
   ```

2. Create and push tag:

   ```bash
   git tag v1.0.0
   git push origin v1.0.0
   ```

### 3. Automated Publishing

#### PyPI Publishing (Automatic)

- GitHub Actions automatically publishes to PyPI on tag push
- Uses PyPI Trusted Publishing (OIDC) - no tokens needed
- Creates signed GitHub release with artifacts

#### Bioconda Publishing (Manual)

1. Wait for PyPI release to complete
2. Update `recipes/ssiamb/meta.yaml` with correct SHA256:

   ```bash
   # Get SHA256 from PyPI release
   pip download ssiamb==1.0.0 --no-deps
   shasum -a 256 ssiamb-1.0.0.tar.gz
   ```

3. Fork [bioconda/bioconda-recipes](https://github.com/bioconda/bioconda-recipes)
4. Copy `recipes/ssiamb/` to `recipes/ssiamb/` in the fork
5. Create pull request to bioconda-recipes
6. Address review feedback and wait for merge

#### Galaxy ToolShed Publishing (Manual)

1. Install planemo: `pip install planemo`
2. Test wrapper: `planemo test galaxy/ssiamb.xml` (may fail until bioconda is available)
3. Create account on [Galaxy ToolShed](https://toolshed.g2.bx.psu.edu/)
4. Upload wrapper:

   ```bash
   cd galaxy/
   planemo shed_upload --shed_target toolshed
   ```

### 4. Post-Release

1. Verify all distributions:
   - PyPI: <https://pypi.org/project/ssiamb/>
   - Bioconda: <https://anaconda.org/bioconda/ssiamb>
   - Galaxy ToolShed: <https://toolshed.g2.bx.psu.edu/>
   - BioContainers: <https://quay.io/repository/biocontainers/ssiamb>

2. Update documentation if needed
3. Announce release

### Version Numbering

- Use semantic versioning: `MAJOR.MINOR.PATCH`
- Galaxy wrapper versions: `SOFTWARE_VERSION+galaxy0` (increment galaxy# for wrapper-only changes)
- Pre-releases: `1.0.0rc1`, `1.0.0a1`, etc.

### Troubleshooting

See `PYPI_SETUP.md` for PyPI Trusted Publishing configuration details.

## Citation

> **Note**: Citation information will be provided upon publication.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
