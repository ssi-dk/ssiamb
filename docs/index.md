# SSIAMB

Welcome to the SSIAMB documentation!

SSIAMB (SSI Ambiguous Site Detection Tool) is a bioinformatics tool for detecting and analyzing ambiguous sites in bacterial genomic sequences.

## Overview

SSIAMB provides comprehensive analysis of ambiguous sites in genomic data by:

- Processing BAM alignment files to identify ambiguous sites
- Calculating depth and quality metrics
- Generating detailed reports and summaries
- Supporting multiple alignment and variant calling workflows

## Quick Start

### Installation

```bash
pip install ssiamb
```

### Basic Usage

```bash
# Analyze a BAM file for ambiguous sites
ssiamb analyze sample.bam --reference ref.fasta --output results/

# Process with custom thresholds
ssiamb analyze sample.bam --reference ref.fasta --min-depth 10 --min-quality 20
```

## Features

- **Multi-format support**: Works with BAM alignment files
- **Flexible workflows**: Compatible with BWA and Minimap2 aligners
- **Quality control**: Comprehensive depth and quality filtering
- **Rich output**: Detailed reports and summary statistics
- **Galaxy integration**: Available as a Galaxy tool

## Documentation Structure

- [API Reference](reference/index.md) - Complete API documentation
- [Galaxy Wrapper](https://github.com/ssi-dk/ssiamb/tree/main/galaxy) - Galaxy tool wrapper

## Support

- [GitHub Issues](https://github.com/ssi-dk/ssiamb/issues) - Bug reports and feature requests
- [Repository](https://github.com/ssi-dk/ssiamb) - Source code and development