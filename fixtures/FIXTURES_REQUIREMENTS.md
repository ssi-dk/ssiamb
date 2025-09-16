# ssiamb Test Fixtures Requirements

## Overview
This document outlines the test fixtures needed for comprehensive ssiamb testing. Fixtures should be **small but realistic** to enable fast testing while covering all functionality.

## 🧬 Core Test Data Requirements

### 1. Primary Test Dataset (Small Bacterial Genome/Plasmid)
**Target**: Small bacterial sequence (~10-50kb) with known variants

#### Required Files:
- **`small_genome.fna`** - Reference/assembly sequence (10-50kb)
- **`sample_R1.fastq.gz`** - Forward reads (5000-10000 reads)  
- **`sample_R2.fastq.gz`** - Reverse reads (5000-10000 reads)
- **`sample_expected_summary.tsv`** - Expected ssiamb output for validation

#### Characteristics Needed:
- **Known ambiguous sites**: Include positions with 10-30% minor allele frequency
- **Good coverage**: ~20-50x average depth
- **Mixed variants**: SNVs, small indels
- **Realistic error profile**: Typical Illumina error patterns

---

## 📁 Test Files by Category

### 2. Self-Mapping Mode Files
```
fixtures/self_mode/
├── assembly.fna              # Small assembly (~10-50kb)
├── reads_R1.fastq.gz         # PE reads that map to assembly  
├── reads_R2.fastq.gz
├── expected_summary.tsv      # Expected ambiguous_summary.tsv output
└── expected_variants.vcf     # Expected variant calls (for --emit-vcf testing)
```

### 3. Reference-Mapped Mode Files
```
fixtures/ref_mode/
├── reference.fna             # Reference genome (same sequence as assembly above)
├── reads_R1.fastq.gz         # Same reads as self-mode (for comparison)
├── reads_R2.fastq.gz  
├── expected_summary.tsv      # Expected output (should match self-mode)
└── species_test/
    ├── bracken_output.tsv    # Sample Bracken classification output
    └── species_mapping.txt   # Species name → reference file mapping
```

### 4. Admin Reference Directory Structure
```
fixtures/admin_refs/
├── Escherichia_coli.fna      # Example reference
├── Escherichia_coli.fna.fai  # samtools index
├── Escherichia_coli.mmi      # minimap2 index  
├── Escherichia_coli.bwt      # bwa-mem2 index files
├── Escherichia_coli.pac
├── Escherichia_coli.0123
└── Salmonella_enterica.fna   # Second species for testing
```

### 5. Reuse/Integration Testing Files
```
fixtures/reuse/
├── pre_mapped.bam           # Pre-computed BAM (single-sample)
├── pre_mapped.bai           # BAM index
├── pre_called.vcf.gz        # Pre-computed VCF (single-sample)  
├── pre_called.vcf.gz.tbi    # VCF index
├── mismatched_sample.bam    # BAM with different sample name (for error testing)
└── multi_sample.vcf.gz      # Multi-sample VCF (should cause error)
```

### 6. Edge Cases & Error Testing
```
fixtures/edge_cases/
├── empty_reads_R1.fastq.gz     # Empty FASTQ files
├── empty_reads_R2.fastq.gz
├── malformed.vcf               # Malformed VCF file
├── low_coverage_R1.fastq.gz    # Very low coverage reads
├── low_coverage_R2.fastq.gz
├── no_variants_R1.fastq.gz     # Reads with no ambiguous sites
├── no_variants_R2.fastq.gz
└── large_indels_R1.fastq.gz    # Reads with complex variants
```

### 7. Expected Outputs for Validation
```
fixtures/expected_outputs/
├── minimal_summary.tsv         # Minimal expected output
├── full_output_with_vcf.vcf    # Expected VCF with --emit-vcf
├── depth_matrix.tsv            # Expected --emit-matrix output
├── per_contig_breakdown.tsv    # Expected --emit-per-contig output
└── run_provenance.json         # Expected --emit-provenance output
```

---

## 🎯 Specific File Requirements

### FASTQ Files Specifications:
- **Format**: Illumina paired-end, gzipped
- **Read length**: 150bp (typical modern Illumina)
- **Coverage**: 20-50x for main tests, 5-10x for low coverage tests
- **Quality scores**: Realistic Phred33 quality scores
- **Insert size**: ~300-500bp typical insert size distribution

### FASTA Files Specifications:
- **Size**: 10-50kb (small enough for fast testing)
- **Content**: Realistic bacterial sequence (could be a plasmid or small genome fragment)
- **Headers**: Simple, clean headers (e.g., `>test_sequence`)

### VCF Files Specifications:
- **Format**: VCF 4.2+ format
- **Single-sample**: Must contain exactly one sample
- **Variants**: Mix of SNVs and small indels
- **Annotations**: Basic INFO fields (DP, AF, etc.)

---

## 🔬 Suggested Test Organism/Sequence

### Option 1: Small Bacterial Plasmid
- **Advantages**: Small size, well-characterized, realistic
- **Example**: pUC19 (2.7kb) or similar cloning vector
- **Simulate**: Add some SNVs/indels for ambiguous sites

### Option 2: Genome Fragment  
- **Advantages**: More realistic than plasmid
- **Example**: 20-50kb fragment from E. coli or Salmonella
- **Source**: Extract from reference genome

### Option 3: Synthetic Sequence
- **Advantages**: Full control over variants and characteristics
- **Generate**: Use tools like ART or PBSIM to simulate reads

---

## 📋 Creation Checklist

When providing fixtures, please ensure:

- [ ] **Small file sizes** (< 1MB each for fast CI/testing)
- [ ] **Realistic data** (actual biological sequences)
- [ ] **Known truth set** (documented expected results)
- [ ] **Gzipped FASTQ** files (space efficiency)
- [ ] **Consistent naming** (follow the structure above)
- [ ] **Sample metadata** (document what organism, coverage, etc.)

---

## 🚀 Priority Order

1. **HIGH PRIORITY** - Basic self-mode test set (assembly + reads + expected output)
2. **HIGH PRIORITY** - Reference-mode test set (reference + same reads)  
3. **MEDIUM PRIORITY** - Edge cases (empty files, low coverage)
4. **MEDIUM PRIORITY** - Reuse testing (pre-computed BAM/VCF)
5. **LOW PRIORITY** - Admin directory structure examples

Would you like me to help you create any of these files, or do you have access to suitable test datasets we can adapt?