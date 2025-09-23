# ssiamb — Ambiguous Sites Counter

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.12+](https://img.shields.io/badge/python-3.12+-blue.svg)](https://www.python.org/downloads/)

**Author:** Povilas Matusevicius <pmat@ssi.dk>  
**Repository:** [https://github.com/ssi-dk/ssiamb](https://github.com/ssi-dk/ssiamb)  
**License:** MIT  
**Minimum Python:** 3.12+

## Overview

`ssiamb` computes an "ambiguous sites" metric for bacterial whole genome sequencing (WGS) as a measure of within-sample heterogeneity. This tool modernizes and standardizes the lab's prior definition while providing robust packaging, CLI interface, and Galaxy integration capabilities.

### What are "Ambiguous Sites"?

An ambiguous site is a genomic position with:
 
- **Sufficient coverage**: Depth ≥ `dp_min` (default: 10)
- **Minor-allele signal**: Minor-allele fraction (MAF) ≥ `maf_min` (default: 0.10)

These metrics are determined from variant calls after normalization and atomization, counting **once per locus** (multi-allelic sites count once if any ALT passes the thresholds).

### Supported Modes

#### Self-mapping Mode (`ssiamb self`)
- **Input**: Reads → Sample's own assembly
- **Use case**: Analyze heterogeneity against the sample's assembled genome
- **Mapping space**: Uses the assembly as reference

#### Reference-mapped Mode (`ssiamb ref`)
- **Input**: Reads → Species canonical reference
- **Use case**: Compare against standardized reference genomes  
- **Reference selection**: Via admin directory, user override, or Bracken classification

#### Summarize Mode (`ssiamb summarize`)
- **Input**: Pre-computed VCF + BAM files
- **Use case**: Reanalyze existing variant calls with different thresholds
- **Speed**: Fast analysis without remapping or variant calling

### Key Features

- **🧬 Automatic Reference Management**: Download and index references from NCBI RefSeq
- **🔧 Flexible Mapping**: Support for minimap2 (default) and bwa-mem2
- **📊 Multiple Variant Callers**: BBTools (default) and bcftools
- **📋 Comprehensive Outputs**: Summary TSV (always) + optional VCF, BED, matrices, per-contig analysis
- **📏 Depth Analysis**: Using mosdepth (default) or samtools
- **♻️ Reusable Workflows**: Accept pre-computed BAM/VCF files
- **🧪 Galaxy Integration**: Designed for workflow environments
- **✅ Quality Control**: Configurable thresholds with sensible defaults
- **🧪 Robust Testing**: Comprehensive test suite with smoke testing

## Installation

### Development Installation (Recommended)

For development or to get the latest features:

```bash
# Clone the repository
git clone https://github.com/ssi-dk/ssiamb.git
cd ssiamb

# Create conda environment with dependencies
conda env create -f environment.yml
conda activate ssiamb

# Install in editable mode
pip install -e .

# Verify installation
ssiamb --help
```

### Conda Environment Setup

The project includes an `environment.yml` file with all required dependencies:

```bash
# Create environment
conda env create -f environment.yml

# Activate environment  
conda activate ssiamb

# Install ssiamb
pip install -e .
```

### Future Distribution Methods

When stable releases are published, the package will also be available via:

```bash
# Future installation via pip (PyPI)
pip install ssiamb

# Future installation via conda (Bioconda)
conda install -c bioconda ssiamb
```

### External Tool Dependencies

`ssiamb` requires several external bioinformatics tools. These are included in the conda environment:

- **Mapping**: `minimap2` and/or `bwa-mem2`
- **Variant calling**: BBTools and `bcftools`
- **Depth analysis**: `mosdepth` and `samtools`
- **VCF processing**: `bcftools` (for normalization)

## Admin Reference Directory Setup

For reference-mapped mode, you need an admin reference directory containing indexed reference genomes. `ssiamb` includes a built-in downloader to automatically fetch and index references from NCBI RefSeq.

### Setting Up References

1. **Set environment variable** (recommended):
   ```bash
   export SSIAMB_REF_DIR=/path/to/references
   ```

2. **Download common bacterial references**:
   ```bash
   # Download single species
   python -m ssiamb.refseq download --species "Escherichia coli" --output-dir $SSIAMB_REF_DIR
   
   # Download multiple species  
   python -m ssiamb.refseq download --species "Salmonella enterica" --output-dir $SSIAMB_REF_DIR
   python -m ssiamb.refseq download --species "Staphylococcus aureus" --output-dir $SSIAMB_REF_DIR
   ```

3. **Verify setup**:
   ```bash
   ls $SSIAMB_REF_DIR
   # Should show: Escherichia_coli.fna, Escherichia_coli.fna.mmi, Escherichia_coli.fna.bwa.*, etc.
   ```

### Reference Downloader Features

- **Automatic Selection**: Chooses best RefSeq reference genome (complete > chromosome > scaffold)
- **Index Generation**: Creates both minimap2 (`.mmi`) and bwa-mem2 (`.bwa.*`) indexes
- **Species Normalization**: Handles common name variations and aliases
- **Progress Reporting**: Shows download progress with rich progress bars
- **Fallback Logic**: Tries multiple genomes if primary choice fails

### Common Species Examples

```bash
# Popular bacterial pathogens
python -m ssiamb.refseq download --species "Escherichia coli" --output-dir $SSIAMB_REF_DIR
python -m ssiamb.refseq download --species "Salmonella enterica" --output-dir $SSIAMB_REF_DIR  
python -m ssiamb.refseq download --species "Staphylococcus aureus" --output-dir $SSIAMB_REF_DIR
python -m ssiamb.refseq download --species "Streptococcus pneumoniae" --output-dir $SSIAMB_REF_DIR
python -m ssiamb.refseq download --species "Klebsiella pneumoniae" --output-dir $SSIAMB_REF_DIR
```

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

### Basic Usage

```bash
# Check what would be done (dry run)
ssiamb self --r1 sample_R1.fastq.gz --r2 sample_R2.fastq.gz --assembly sample.fna --dry-run

# Self-mapping mode: analyze reads against sample's own assembly
ssiamb self --r1 sample_R1.fastq.gz --r2 sample_R2.fastq.gz --assembly sample.fna

# Reference-mapped mode: analyze against species reference
ssiamb ref --r1 sample_R1.fastq.gz --r2 sample_R2.fastq.gz --species "Escherichia coli"

# Summarize existing VCF and BAM files
ssiamb summarize --vcf sample.vcf.gz --bam sample.bam
```

### Comprehensive Examples

#### Self-Mapping Mode

```bash
# Basic self-mapping
ssiamb self \
  --r1 reads_R1.fastq.gz \
  --r2 reads_R2.fastq.gz \
  --assembly assembly.fna \
  --sample MySample \
  --outdir results/

# With custom thresholds and optional outputs
ssiamb self \
  --r1 reads_R1.fastq.gz \
  --r2 reads_R2.fastq.gz \
  --assembly assembly.fna \
  --dp-min 15 \
  --maf-min 0.05 \
  --emit-vcf \
  --emit-bed \
  --emit-matrix \
  --threads 8

# Using bwa-mem2 instead of minimap2
ssiamb self \
  --r1 reads_R1.fastq.gz \
  --r2 reads_R2.fastq.gz \
  --assembly assembly.fna \
  --mapper bwa-mem2 \
  --caller bcftools
```

#### Reference-Mapping Mode

```bash
# Using species name (requires admin reference directory)
export SSIAMB_REF_DIR=/path/to/references
ssiamb ref \
  --r1 reads_R1.fastq.gz \
  --r2 reads_R2.fastq.gz \
  --species "Escherichia coli" \
  --sample MySample

# Using direct reference file
ssiamb ref \
  --r1 reads_R1.fastq.gz \
  --r2 reads_R2.fastq.gz \
  --reference reference_genome.fna \
  --sample MySample

# Using Bracken classification results
ssiamb ref \
  --r1 reads_R1.fastq.gz \
  --r2 reads_R2.fastq.gz \
  --bracken sample.bracken \
  --ref-dir /path/to/references \
  --min-bracken-frac 0.8
```

#### Advanced Usage

```bash
# Output to stdout (no files written)
ssiamb self \
  --r1 reads_R1.fastq.gz \
  --r2 reads_R2.fastq.gz \
  --assembly assembly.fna \
  --stdout

# Reuse existing BAM file
ssiamb self \
  --r1 reads_R1.fastq.gz \
  --r2 reads_R2.fastq.gz \
  --assembly assembly.fna \
  --bam existing_alignment.bam

# Append to existing TSV instead of overwriting
ssiamb self \
  --r1 reads_R1.fastq.gz \
  --r2 reads_R2.fastq.gz \
  --assembly assembly.fna \
  --tsv-mode append

# Comprehensive output with provenance
ssiamb self \
  --r1 reads_R1.fastq.gz \
  --r2 reads_R2.fastq.gz \
  --assembly assembly.fna \
  --emit-vcf \
  --emit-bed \
  --emit-matrix \
  --emit-per-contig \
  --emit-provenance \
  --emit-multiqc
```

### Testing Your Installation

Run the built-in smoke test to verify everything works:

```bash
# Quick test (skips downloads)
python smoke_test.py --skip-downloads

# Full test including reference downloads (slower)
python smoke_test.py
```

## Output Files and Formats

### Primary Output

#### `ambiguous_summary.tsv`
Always generated. Single-row summary with comprehensive metrics:

| Column | Description | Example |
|--------|-------------|---------|
| `sample` | Sample identifier | `MySample` |
| `mode` | Analysis mode | `self`, `ref` |
| `mapper` | Alignment tool used | `minimap2`, `bwa-mem2` |
| `caller` | Variant caller used | `bbtools`, `bcftools` |
| `dp_min` | Minimum depth threshold | `10` |
| `maf_min` | Minimum MAF threshold | `0.1` |
| `ambiguous_snv_count` | Number of ambiguous SNVs | `42` |
| `ambiguous_snv_per_mb` | SNVs per megabase | `15.23` |
| `callable_bases` | Bases with adequate coverage | `4651234` |
| `genome_length` | Total genome length | `4800000` |
| `breadth_10x` | Fraction covered at 10x+ | `0.9691` |
| `ref_label` | Reference identifier | `Escherichia_coli|GCF_000005825.2` |
| `runtime_sec` | Analysis runtime | `245.67` |

### Optional Outputs

Enable with command-line flags:

#### `--emit-vcf`: Variant Call Format
- **File**: `{SAMPLE}.ambiguous_sites.vcf.gz` + `.tbi` index
- **Content**: Normalized, atomized variants passing thresholds
- **Annotations**: Custom INFO fields (MAF, AMBIG flag, etc.)
- **Format**: BGzip compressed, tabix indexed

#### `--emit-bed`: Browser Extensible Data
- **File**: `{SAMPLE}.ambiguous_sites.bed.gz` + `.tbi` index  
- **Content**: Genomic intervals of ambiguous sites
- **Columns**: `chrom`, `start`, `end`, `name`, `score`, `strand`, `sample`, `variant_class`, `ref`, `alt`, `maf`, `dp`
- **Format**: BGzip compressed, tabix indexed

#### `--emit-matrix`: Depth×MAF Matrix
- **File**: `{SAMPLE}.ambiguity_matrix.tsv.gz`
- **Content**: 100×51 cumulative count matrix
- **Rows**: Depth thresholds (1-100)
- **Columns**: MAF bins (0-50, representing 0.00-0.50)
- **Values**: Cumulative variant counts

#### `--emit-per-contig`: Per-Contig Summary
- **File**: `{SAMPLE}.per_contig_summary.tsv`
- **Content**: Breakdown by chromosome/contig
- **Columns**: `sample`, `contig`, `length`, `callable_bases_10x`, `breadth_10x`, `ambiguous_snv_count`, etc.

#### `--emit-provenance`: Analysis Provenance
- **File**: `run_provenance.json`
- **Content**: Complete analysis parameters, tool versions, runtime info
- **Format**: JSON array (one entry per sample)

#### `--emit-multiqc`: MultiQC Integration
- **File**: `{SAMPLE}.multiqc.tsv`
- **Content**: Curated metrics for MultiQC reporting
- **Use case**: Integration with MultiQC pipelines

### Output Directory Structure

```
results/
├── ambiguous_summary.tsv                    # Always generated
├── MySample.ambiguous_sites.vcf.gz          # --emit-vcf
├── MySample.ambiguous_sites.vcf.gz.tbi      # VCF index
├── MySample.ambiguous_sites.bed.gz          # --emit-bed  
├── MySample.ambiguous_sites.bed.gz.tbi      # BED index
├── MySample.ambiguity_matrix.tsv.gz         # --emit-matrix
├── MySample.per_contig_summary.tsv          # --emit-per-contig
├── MySample.multiqc.tsv                     # --emit-multiqc
└── run_provenance.json                      # --emit-provenance
```

## CLI Reference

### Global Options

| Option | Default | Description |
|--------|---------|-------------|
| `--threads` | `4` | Number of CPU threads |
| `--outdir` | `.` | Output directory |
| `--sample` | *inferred* | Sample name (required if auto-inference fails) |
| `--dp-min` | `10` | Minimum depth for ambiguous sites |
| `--maf-min` | `0.1` | Minimum minor allele frequency |
| `--dp-cap` | `100` | Maximum depth cap (clipped to 100) |
| `--mapper` | `minimap2` | Alignment tool (`minimap2`, `bwa-mem2`) |
| `--caller` | `bbtools` | Variant caller (`bbtools`, `bcftools`) |
| `--depth-tool` | `mosdepth` | Depth analysis tool (`mosdepth`, `samtools`) |
| `--require-pass` | `False` | Only use PASS variants |
| `--tsv-mode` | `overwrite` | TSV handling (`overwrite`, `append`, `fail`) |
| `--stdout` | `False` | Write summary to stdout only |

### Command-Specific Options

#### `ssiamb self`
| Option | Required | Description |
|--------|----------|-------------|
| `--r1` | ✅ | Forward reads (FASTQ, gzipped OK) |
| `--r2` | ✅ | Reverse reads (FASTQ, gzipped OK) |
| `--assembly` | ✅ | Assembly FASTA file |
| `--vcf` | ❌ | Reuse existing VCF file |
| `--bam` | ❌ | Reuse existing BAM file |

#### `ssiamb ref`
| Option | Required | Description |
|--------|----------|-------------|
| `--r1` | ✅ | Forward reads (FASTQ, gzipped OK) |
| `--r2` | ✅ | Reverse reads (FASTQ, gzipped OK) |
| `--reference` | ❌* | Direct reference FASTA |
| `--species` | ❌* | Species name for lookup |
| `--bracken` | ❌* | Bracken classification file |
| `--ref-dir` | ❌ | Admin reference directory |
| `--min-bracken-frac` | `0.70` | Minimum Bracken fraction |
| `--min-bracken-reads` | `100000` | Minimum Bracken reads |
| `--on-fail` | `error` | Bracken failure action (`error`, `self`) |

*One of `--reference`, `--species`, or `--bracken` is required.

#### `ssiamb summarize`
| Option | Required | Description |
|--------|----------|-------------|
| `--vcf` | ✅ | VCF file to summarize |
| `--bam` | ✅ | BAM file for denominator |
| `--output` | ❌ | Output file path |

## Error Codes and Troubleshooting

### Exit Codes

`ssiamb` follows a structured exit code system for programmatic handling:

- **0**: Success
- **1**: CLI/input errors (missing files, invalid sample names, bad arguments)
- **2**: Reference mode selection errors (species not found, Bracken failures)
- **3**: Reuse compatibility errors (VCF/BAM mismatch with reference)
- **4**: External tool failures (missing tools, tool execution errors)
- **5**: QC failures (only when `--qc-action fail` is enabled)

### Common Issues and Solutions

#### Missing External Tools

**Error**: `Command not found: minimap2`

**Solution**: Install dependencies via conda:
```bash
conda install -c bioconda minimap2 bwa-mem2 bcftools samtools mosdepth bbmap
```

#### Reference Directory Issues

**Error**: `Species 'Escherichia coli' not found in admin directory`

**Solutions**:
```bash
# Option 1: Download the species reference
python -m ssiamb.refseq download --species "Escherichia coli" --output-dir $SSIAMB_REF_DIR

# Option 2: Use direct reference file
ssiamb ref --reference /path/to/ecoli.fna --r1 reads_R1.fastq.gz --r2 reads_R2.fastq.gz

# Option 3: Set environment variable
export SSIAMB_REF_DIR=/path/to/your/references
```

#### Index File Issues

**Error**: `Reference found but indexes missing for 'Escherichia_coli': minimap2 index`

**Solution**: Regenerate indexes:
```bash
cd $SSIAMB_REF_DIR
minimap2 -d Escherichia_coli.fna.mmi Escherichia_coli.fna
bwa-mem2 index Escherichia_coli.fna
```

#### Sample Name Issues

**Error**: `Sample name could not be safely inferred`

**Solution**: Explicitly provide sample name:
```bash
ssiamb self --r1 reads_R1.fastq.gz --r2 reads_R2.fastq.gz --assembly assembly.fna --sample MySample
```

#### Memory Issues

**Error**: `BBTools out of memory`

**Solutions**:
```bash
# Reduce threads
ssiamb self --threads 2 --r1 reads_R1.fastq.gz --r2 reads_R2.fastq.gz --assembly assembly.fna

# Set BBTools memory limit  
ssiamb self --bbtools-mem 4g --r1 reads_R1.fastq.gz --r2 reads_R2.fastq.gz --assembly assembly.fna
```

#### Permission Issues

**Error**: `Cannot write to output directory`

**Solutions**:
```bash
# Create directory with proper permissions
mkdir -p /path/to/output
chmod 755 /path/to/output

# Or use different output directory
ssiamb self --outdir ~/results --r1 reads_R1.fastq.gz --r2 reads_R2.fastq.gz --assembly assembly.fna
```

### Debugging Tips

1. **Use dry-run mode** to validate inputs:
   ```bash
   ssiamb self --dry-run --r1 reads_R1.fastq.gz --r2 reads_R2.fastq.gz --assembly assembly.fna
   ```

2. **Check tool versions**:
   ```bash
   minimap2 --version
   bcftools --version
   mosdepth --version
   ```

3. **Validate input files**:
   ```bash
   # Check FASTQ files
   zcat reads_R1.fastq.gz | head -4
   
   # Check FASTA files
   head assembly.fna
   ```

4. **Run smoke test**:
   ```bash
   python smoke_test.py --skip-downloads
   ```

### Getting Help

- **CLI help**: `ssiamb --help` or `ssiamb COMMAND --help`
- **GitHub Issues**: [Report bugs or request features](https://github.com/ssi-dk/ssiamb/issues)
- **Contact**: Povilas Matusevicius <pmat@ssi.dk>

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

## Testing

The project includes comprehensive testing:

### Running Tests

```bash
# Run all unit tests
python -m pytest tests/ -v

# Run specific test modules
python -m pytest tests/test_refseq.py -v
python -m pytest tests/test_refdir.py -v

# Run smoke test (integration test)
python smoke_test.py --skip-downloads  # Fast version
python smoke_test.py                   # Full version with downloads
```

### Test Coverage

The test suite covers:

- **Unit tests**: Core algorithms, species normalization, reference downloading
- **Integration tests**: Full pipeline validation via smoke test
- **Edge cases**: Error handling, malformed inputs, missing dependencies
- **Mock testing**: External API calls and tool dependencies

### Test Dependencies

Install test dependencies:

```bash
# Via conda environment (recommended)
conda env create -f environment.yml

# Or manually
pip install pytest numpy pysam biopython requests
```

## Development Status

This project has completed its major development milestones:

- ✅ **Planning & Specification** - Comprehensive requirements defined
- ✅ **Repository Bootstrap** - Package structure, CI/CD, documentation
- ✅ **Core Implementation** - CLI, models, and processing pipelines
- ✅ **External Tool Integration** - Mapping and variant calling workflows
- ✅ **Reference Management** - Automatic RefSeq download and indexing
- ✅ **Testing & Validation** - Unit tests, integration testing, smoke tests
- 🚧 **Packaging & Distribution** - Bioconda, containers, Galaxy tools

### Recent Features

- **RefSeq Integration**: Automatic reference genome downloading from NCBI
- **Robust Indexing**: Automatic minimap2 and bwa-mem2 index generation
- **Enhanced Testing**: Comprehensive unit tests and smoke testing
- **Improved CLI**: Better help text, error messages, and validation
- **Output Flexibility**: Centralized directory handling with fallbacks

 
## Multi-Sample Processing

`ssiamb` supports batch processing via manifest files for analyzing multiple samples efficiently.

### Manifest Files

Create TSV files listing samples and their inputs:

#### Self-Mode Manifest

```tsv
sample	r1	r2	assembly	bam	vcf
Sample1	data/Sample1_R1.fastq.gz	data/Sample1_R2.fastq.gz	assemblies/Sample1.fna		
Sample2	data/Sample2_R1.fastq.gz	data/Sample2_R2.fastq.gz	assemblies/Sample2.fna	existing/Sample2.bam	
Sample3	data/Sample3_R1.fastq.gz	data/Sample3_R2.fastq.gz	assemblies/Sample3.fna		existing/Sample3.vcf.gz
```

#### Reference-Mode Manifest

```tsv
sample	r1	r2	bracken	reference	species	bam	vcf
Sample1	data/Sample1_R1.fastq.gz	data/Sample1_R2.fastq.gz	Sample1.bracken			
Sample2	data/Sample2_R1.fastq.gz	data/Sample2_R2.fastq.gz		ref/custom.fna		
Sample3	data/Sample3_R1.fastq.gz	data/Sample3_R2.fastq.gz			Escherichia coli	
```

### Running with Manifests

```bash
# Process self-mode manifest
ssiamb self --manifest samples.tsv --outdir results/

# Process reference-mode manifest  
ssiamb ref --manifest samples.tsv --ref-dir $SSIAMB_REF_DIR --outdir results/

# With custom settings
ssiamb self --manifest samples.tsv --dp-min 15 --emit-vcf --threads 8
```

### Manifest Features

- **Relative paths**: Resolved relative to manifest file location
- **Optional columns**: Empty cells skip optional inputs (bam, vcf)
- **Comments**: Lines starting with `#` are ignored
- **Sequential processing**: Samples processed one at a time
- **Consolidated output**: Single `ambiguous_summary.tsv` with all samples

## Performance Considerations

### Resource Usage

- **CPU**: Scales with `--threads` parameter (default: 4)
- **Memory**: 
  - BBTools: 4-8GB (adjustable with `--bbtools-mem`)
  - bwa-mem2: 2-4GB 
  - minimap2: 1-2GB
- **Disk**: ~2-5x input size for intermediate files
- **Network**: Required for RefSeq downloads only

### Optimization Tips

1. **Use appropriate thread count**:
   ```bash
   # For high-memory systems
   ssiamb self --threads 16 --r1 reads_R1.fastq.gz --r2 reads_R2.fastq.gz --assembly assembly.fna
   ```

2. **Reuse intermediate files**:
   ```bash
   # First run - saves BAM
   ssiamb self --r1 reads_R1.fastq.gz --r2 reads_R2.fastq.gz --assembly assembly.fna
   
   # Subsequent runs - reuses BAM
   ssiamb self --r1 reads_R1.fastq.gz --r2 reads_R2.fastq.gz --assembly assembly.fna --bam sample.sorted.bam --dp-min 15
   ```

3. **Pre-download references**:
   ```bash
   # Download all needed references first
   python -m ssiamb.refseq download --species "Escherichia coli" --output-dir $SSIAMB_REF_DIR
   python -m ssiamb.refseq download --species "Salmonella enterica" --output-dir $SSIAMB_REF_DIR
   ```

4. **Use faster mappers for large datasets**:
   ```bash
   # minimap2 is generally faster
   ssiamb self --mapper minimap2 --r1 reads_R1.fastq.gz --r2 reads_R2.fastq.gz --assembly assembly.fna
   ```

### Typical Runtimes

| Dataset Size | Mode | Mapper | Runtime (4 cores) |
|-------------|------|--------|--------------------|
| 1M PE reads | self | minimap2 | 2-5 minutes |
| 1M PE reads | ref | minimap2 | 3-7 minutes |
| 5M PE reads | self | minimap2 | 8-15 minutes |
| 5M PE reads | ref | bwa-mem2 | 15-25 minutes |

*Times vary based on genome size, coverage, and hardware*

## Contributing

This project is developed by the SSI team. The codebase is now feature-complete with comprehensive testing.

### Development Setup

```bash
# Clone repository
git clone https://github.com/ssi-dk/ssiamb.git
cd ssiamb

# Create development environment
conda env create -f environment.yml
conda activate ssiamb

# Install in editable mode
pip install -e .

# Run tests
python -m pytest tests/ -v
python smoke_test.py --skip-downloads
```

### Code Quality

Before submitting changes:

1. **Run the test suite**:
   ```bash
   python -m pytest tests/ -v
   ```

2. **Run smoke tests**:
   ```bash
   python smoke_test.py
   ```

3. **Check code style** (if configured):
   ```bash
   black src/ tests/
   isort src/ tests/
   ```

### Contact

For questions, contributions, or issues:

- **Primary Contact**: Povilas Matusevicius <pmat@ssi.dk>
- **GitHub Issues**: [Report bugs or request features](https://github.com/ssi-dk/ssiamb/issues)
- **Repository**: [https://github.com/ssi-dk/ssiamb](https://github.com/ssi-dk/ssiamb)

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
