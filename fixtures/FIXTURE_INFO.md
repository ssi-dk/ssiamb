# ssiamb Test Fixtures - Listeria monocytogenes listeria_test_sample

## Overview

These fixtures are based on **Listeria monocytogenes** sample **listeria_test_sample** from SSI's production pipeline.

## Sample Information

- **Organism**: *Listeria monocytogenes* 
- **Sample ID**: listeria_test_sample
- **Source**: SSI production pipeline
- **Species Confidence**: 99.52% *L. monocytogenes* (Bracken classification)
- **Assembly Size**: 2.8 MB (realistic bacterial genome size)
- **Sequencing**: Illumina paired-end, 150bp reads

## Fixture Files

### Self-Mapping Mode (`fixtures/self_mode/`)
- **`assembly.fna`** (2.8M) - Complete assembly of listeria_test_sample
- **`reads_R1.fastq.gz`** (702K) - 10,000 subsampled forward reads  
- **`reads_R2.fastq.gz`** (742K) - 10,000 subsampled reverse reads

### Reference-Mapping Mode (`fixtures/ref_mode/`)  
- **`reference.fna`** (2.8M) - L. monocytogenes EGD-e reference genome (NC_003210.1)
- **`reads_R1.fastq.gz`** (702K) - Identical reads for comparison
- **`reads_R2.fastq.gz`** (742K) - Identical reads for comparison

### Supporting Files (`fixtures/supporting_files/`)
- **`bracken_listeria.tsv`** (979B) - Species classification output

## Usage

### Self-Mapping Test
```bash
ssiamb self \
  --r1 fixtures/self_mode/reads_R1.fastq.gz \
  --r2 fixtures/self_mode/reads_R2.fastq.gz \
  --assembly fixtures/self_mode/assembly.fna
```

### Reference-Mapping Test  
```bash
ssiamb ref \
  --r1 fixtures/ref_mode/reads_R1.fastq.gz \
  --r2 fixtures/ref_mode/reads_R2.fastq.gz \
  --reference fixtures/ref_mode/reference.fna
```

## Expected Characteristics

### Coverage
- **Target Coverage**: ~10-20x (estimated from 10k reads × 150bp / 2.8Mb)
- **Read Length**: 150bp Illumina PE
- **Insert Size**: ~300-500bp (typical)

### Expected Ambiguous Sites
- **Realistic heterozygosity**: Sample bacterial isolate
- **Known variants**: Should have some ambiguous sites from sequencing/assembly process  
- **Mode comparison**: Self-mode vs ref-mode will show different results due to reference genome differences

## File Generation

Files were created from original data using:

```bash
# Assembly (already assembled)
cp original/listeria_test_sample.fasta fixtures/*/assembly.fna

# Subsampled reads (seed 123 for reproducibility)
seqtk sample -s 123 original_R1.fastq.gz 10000 | gzip > reads_R1.fastq.gz
seqtk sample -s 123 original_R2.fastq.gz 10000 | gzip > reads_R2.fastq.gz

# Species classification
cp original/bracken.txt supporting_files/bracken_listeria.tsv
```

## Quality Control

- ✅ **File sizes**: Manageable for testing (~700KB reads, 2.8MB assembly)
- ✅ **Realistic data**: Actual production sample from SSI pipeline  
- ✅ **Species purity**: 99.52% single species
- ✅ **Paired reads**: Same seed ensures proper pairing
- ✅ **Complete**: Both self and reference mode files available

## Next Steps

1. **Validate assembly quality**: Check N50, contiguity
2. **Generate expected results**: Run existing tools to create truth sets
3. **Add more edge cases**: Low coverage, complex variants, etc.
4. **Create small reference directory**: For admin ref testing

---
*Generated from SSI production data on September 16, 2025*