# SSIAMB Documentation

Welcome to the SSIAMB documentation!

SSIAMB (SSI Ambiguous Site Detection Tool) is a comprehensive bioinformatics tool for detecting and analyzing ambiguous sites in bacterial genomic sequences.

## Overview

SSIAMB provides comprehensive analysis of ambiguous sites in bacterial WGS data by:

- **Self-mapping mode**: Analyzing reads against sample's own assembly
- **Reference-mapping mode**: Comparing against species canonical references  
- **Automatic reference management**: Download and index references from NCBI RefSeq
- **Flexible workflows**: Supporting multiple mappers and variant callers
- **Rich outputs**: Summary statistics, VCF/BED files, matrices, and provenance

## Quick Start

### Installation

```bash
# Clone and setup development environment
git clone https://github.com/ssi-dk/ssiamb.git
cd ssiamb
conda env create -f environment.yml
conda activate ssiamb
pip install -e .
```

### Basic Usage

```bash
# Self-mapping mode
ssiamb self --r1 sample_R1.fastq.gz --r2 sample_R2.fastq.gz --assembly sample.fna

# Reference-mapping mode
export SSIAMB_REF_DIR=/path/to/references
ssiamb ref --r1 sample_R1.fastq.gz --r2 sample_R2.fastq.gz --species "Escherichia coli"

# Download reference genomes
python -m ssiamb.refseq download --species "Escherichia coli" --output-dir $SSIAMB_REF_DIR
```

## Key Features

- **🧬 Automatic Reference Management**: Download and index references from NCBI RefSeq
- **🔧 Multiple Mappers**: minimap2 (default) and bwa-mem2 support
- **📊 Multiple Callers**: BBTools (default) and bcftools variant calling
- **📋 Comprehensive Outputs**: TSV summaries, VCF/BED files, depth matrices
- **♻️ Workflow Integration**: Reuse existing BAM/VCF files
- **🧪 Robust Testing**: Comprehensive test suite with smoke testing
- **🌟 Galaxy Integration**: Available as a Galaxy tool wrapper

## What are Ambiguous Sites?

An ambiguous site is a genomic position with:

- **Sufficient coverage**: Depth ≥ `dp_min` (default: 10)
- **Minor-allele signal**: Minor-allele fraction (MAF) ≥ `maf_min` (default: 0.10)

These metrics provide a measure of within-sample heterogeneity, useful for quality control and population genetics analysis.

## Documentation Structure

- **[Main README](../README.md)** - Comprehensive installation and usage guide
- **[API Reference](reference/index.md)** - Complete API documentation  
- **[Galaxy Wrapper](https://github.com/ssi-dk/ssiamb/tree/main/galaxy)** - Galaxy tool wrapper
- **[Specification](../spec.md)** - Detailed technical specification

## Support and Development

- **[GitHub Repository](https://github.com/ssi-dk/ssiamb)** - Source code and development
- **[GitHub Issues](https://github.com/ssi-dk/ssiamb/issues)** - Bug reports and feature requests
- **Contact**: Povilas Matusevicius <pmat@ssi.dk>

The project is feature-complete and actively maintained by the SSI team.
