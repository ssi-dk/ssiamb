# Placeholder expected outputs - to be generated when ssiamb tool is implemented

This directory contains expected output files that ssiamb should produce:

## Generated Files (TODO)
- `self_mode_summary.tsv` - Expected ambiguous_summary.tsv for self mode
- `ref_mode_summary.tsv` - Expected ambiguous_summary.tsv for ref mode  
- `self_mode_sites.vcf.gz` - Expected VCF with --emit-vcf (self mode)
- `ref_mode_sites.vcf.gz` - Expected VCF with --emit-vcf (ref mode)
- `variant_matrix.tsv.gz` - Expected --emit-matrix output
- `per_contig_summary.tsv` - Expected --emit-per-contig output
- `run_provenance.json` - Expected --emit-provenance output

## Usage
These files will be used for:
1. **Regression testing** - Ensure outputs don't change unexpectedly
2. **Validation** - Verify tool produces correct results
3. **CI/CD** - Automated testing in GitHub Actions

## Generation
Run ssiamb on fixture data to generate these files:
```bash
ssiamb self --emit-vcf --emit-matrix --emit-per-contig --emit-provenance \
  --r1 fixtures/self_mode/reads_R1.fastq.gz \
  --r2 fixtures/self_mode/reads_R2.fastq.gz \
  --assembly fixtures/self_mode/assembly.fna
```