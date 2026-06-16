# ssiamb

`ssiamb` counts ambiguous SNV sites from trimmed paired-end reads mapped back to
their own assembly.

The first supported workflow matches production use at SSI:

```bash
ssiamb self \
    --r1 reads_R1.fastq.gz \
    --r2 reads_R2.fastq.gz \
    --assembly shovill_assembly.fasta \
    --sample SAMPLE \
    --outdir outputs
```

The tool assumes reads are already trimmed and the assembly has already been
created, for example by Shovill.

## Method

The default flow maps reads to the supplied assembly with `minimap2`, sorts the
alignment with `samtools`, calls variants with BBTools `callvariants.sh`, then
groups VCF records by `(CHROM, POS)` before counting ambiguous loci.

Default ambiguity criteria:

- `DP >= 10`
- `MAF >= 0.10`
- same-position ALT allele frequencies are summed before converting to MAF
- high ALT frequencies are converted with `MAF = min(AF, 1 - AF)`
- filtered caller records remain visible in diagnostics and are not hidden from
  the primary candidate set

## Outputs

For sample `SAMPLE`, the CLI writes:

- `SAMPLE.summary.tsv`: one-row summary
- `SAMPLE.loci.tsv`: grouped locus diagnostics
- `SAMPLE.matrix.tsv`: cumulative depth/MAF matrix
- `SAMPLE.variants.vcf`: raw called variants
- `SAMPLE.provenance.json`: command, thresholds, inputs, outputs, and tool info

## External Tools

The end-to-end `self` workflow requires these commands on `PATH`:

- `minimap2`
- `samtools`
- `callvariants.sh` from BBTools/BBMap

Core VCF grouping and counting logic is pure Python and covered by unit tests.
