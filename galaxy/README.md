# Galaxy Wrapper for ssiamb

This directory contains the Galaxy ToolShed wrapper for the Ambiguous Sites Counter (ssiamb).

## Files

- `ssiamb.xml` - Main Galaxy tool wrapper
- `test-data/` - Test datasets for functional testing
- `README.md` - This file

## Tool Information

- **Owner**: pmat
- **Category**: Sequence Analysis  
- **Auto-links**: Bioconda package and BioContainer
- **Resource requirements**: 4 threads, 16GB memory, 2h walltime

## Testing

To test the wrapper locally:

```bash
# Install planemo
pip install planemo

# Lint the wrapper
planemo lint ssiamb.xml

# Test the wrapper (requires external tools)
planemo test ssiamb.xml
```

## Submission to ToolShed

The wrapper will be submitted to the Galaxy ToolShed under owner "pmat" in the "Sequence Analysis" category.

## Dependencies

The wrapper requires the following Conda packages:

- ssiamb (main tool)
- minimap2, bwa-mem2 (mappers)
- samtools, bcftools, htslib (BAM/VCF processing)
- bbmap (alternative variant caller)
- mosdepth (depth analysis)
- bedtools (coordinate processing)

These are automatically installed via Bioconda when the tool is used in Galaxy.
