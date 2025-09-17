#!/usr/bin/env python3
"""
ssiamb - SSI Ambiguous Site Detection Tool

Command-line interface for detecting ambiguous sites in bacterial genomes
through mapping and variant calling approaches.
"""

from pathlib import Path
from typing import Optional, Annotated
import typer
from rich.console import Console
from rich.table import Table

from .version import __version__
from .models import Mode, TSVMode, DepthTool
from .runner import create_run_plan, execute_plan, run_summarize

app = typer.Typer(
    name="ssiamb",
    help="SSI Ambiguous Site Detection Tool",
    add_completion=False,
    no_args_is_help=True,
)

console = Console()


def version_callback(value: bool) -> None:
    """Print version and exit."""
    if value:
        console.print(f"ssiamb version {__version__}")
        raise typer.Exit()


@app.callback()
def main(
    version: Annotated[
        Optional[bool],
        typer.Option(
            "--version",
            "-V",
            help="Show version and exit",
            callback=version_callback,
            is_eager=True,
        ),
    ] = None,
    config: Annotated[
        Optional[Path],
        typer.Option(
            "--config",
            "-c",
            help="Path to custom configuration file",
            exists=True,
            file_okay=True,
            dir_okay=False,
            readable=True,
        ),
    ] = None,
    verbose: Annotated[
        bool,
        typer.Option(
            "--verbose", "-v", help="Enable verbose output"
        ),
    ] = False,
    quiet: Annotated[
        bool,
        typer.Option(
            "--quiet", "-q", help="Suppress non-error output"
        ),
    ] = False,
    dry_run: Annotated[
        bool,
        typer.Option(
            "--dry-run", help="Show what would be done without executing"
        ),
    ] = False,
) -> None:
    """
    SSI Ambiguous Site Detection Tool.
    
    Detect ambiguous sites in bacterial genomes through mapping and variant calling.
    """
    if verbose and quiet:
        console.print("[red]Error: --verbose and --quiet cannot be used together[/red]")
        raise typer.Exit(1)
    
    # Load configuration if specified
    if config:
        from .config import load_config
        try:
            load_config(config)
            if verbose:
                console.print(f"[green]Loaded configuration from: {config}[/green]")
        except Exception as e:
            console.print(f"[red]Error loading config file {config}: {e}[/red]")
            raise typer.Exit(1)


@app.command()
def self(
    ctx: typer.Context,
    r1: Annotated[
        Path,
        typer.Option(
            "--r1", help="Forward reads (FASTQ, can be gzipped)"
        ),
    ],
    r2: Annotated[
        Path,
        typer.Option(
            "--r2", help="Reverse reads (FASTQ, can be gzipped)"
        ),
    ],
    assembly: Annotated[
        Path,
        typer.Option(
            "--assembly", help="Assembly FASTA file"
        ),
    ],
    sample: Annotated[
        Optional[str],
        typer.Option(
            "--sample", help="Sample name (inferred from filenames if not provided)"
        ),
    ] = None,
    output_dir: Annotated[
        Optional[Path],
        typer.Option(
            "--outdir", "-o", help="Output directory (default: current)"
        ),
    ] = None,
    threads: Annotated[
        int,
        typer.Option(
            "--threads", "-t", help="Number of threads"
        ),
    ] = 4,
    mapper: Annotated[
        str,
        typer.Option(
            "--mapper", help="Mapper to use"
        ),
    ] = "minimap2",
    caller: Annotated[
        str,
        typer.Option(
            "--caller", help="Variant caller to use"
        ),
    ] = "bbtools",
    bbtools_mem: Annotated[
        Optional[str],
        typer.Option(
            "--bbtools-mem", help="BBTools heap memory (e.g., '4g', '8g')"
        ),
    ] = None,
    dp_min: Annotated[
        int,
        typer.Option(
            "--dp-min", help="Minimum depth for ambiguous sites"
        ),
    ] = 10,
    maf_min: Annotated[
        float,
        typer.Option(
            "--maf-min", help="Minimum minor allele frequency"
        ),
    ] = 0.1,
    mapq: Annotated[
        int,
        typer.Option(
            "--mapq", help="Minimum mapping quality for depth analysis"
        ),
    ] = 30,
    dp_cap: Annotated[
        int,
        typer.Option(
            "--dp-cap", help="Maximum depth cap for variant analysis"
        ),
    ] = 100,
    depth_tool: Annotated[
        DepthTool,
        typer.Option(
            "--depth-tool", help="Tool for depth analysis"
        ),
    ] = DepthTool.MOSDEPTH,
    require_pass: Annotated[
        bool,
        typer.Option(
            "--require-pass", help="Only consider variants that pass caller filters"
        ),
    ] = False,
    emit_vcf: Annotated[
        bool,
        typer.Option(
            "--emit-vcf", help="Emit VCF file with ambiguous sites"
        ),
    ] = False,
    emit_bed: Annotated[
        bool,
        typer.Option(
            "--emit-bed", help="Emit BED file with ambiguous sites"
        ),
    ] = False,
    emit_matrix: Annotated[
        bool,
        typer.Option(
            "--emit-matrix", help="Emit variant matrix in TSV format"
        ),
    ] = False,
    emit_per_contig: Annotated[
        bool,
        typer.Option(
            "--emit-per-contig", help="Emit per-contig summary statistics"
        ),
    ] = False,
    emit_multiqc: Annotated[
        bool,
        typer.Option(
            "--emit-multiqc", help="Emit MultiQC-compatible metrics"
        ),
    ] = False,
    emit_provenance: Annotated[
        bool,
        typer.Option(
            "--emit-provenance", help="Emit JSON provenance file"
        ),
    ] = False,
    tsv_mode: Annotated[
        TSVMode,
        typer.Option(
            "--tsv-mode", help="TSV output mode"
        ),
    ] = TSVMode.OVERWRITE,
    stdout: Annotated[
        bool,
        typer.Option(
            "--stdout", help="Write summary to stdout instead of file"
        ),
    ] = False,
) -> None:
    """
    Self-mapping mode: map reads to their own assembly.
    
    Maps paired-end reads to the provided assembly to identify ambiguous sites
    where the assembly may not represent the true sequence.
    """
    try:
        # Validate emit flags with stdout
        if stdout and any([emit_vcf, emit_bed, emit_matrix, emit_per_contig, emit_multiqc, emit_provenance]):
            console.print("[red]Error: --stdout cannot be used with --emit-* flags[/red]")
            raise typer.Exit(1)
        
        # Get global options from context
        dry_run = ctx.parent.params.get("dry_run", False) if ctx.parent else False
        
        # Create execution plan
        plan = create_run_plan(
            mode=Mode.SELF,
            r1=r1,
            r2=r2,
            assembly=assembly,
            sample=sample,
            output_dir=output_dir,
            threads=threads,
            mapper=mapper,
            caller=caller,
            bbtools_mem=bbtools_mem,
            dp_min=dp_min,
            maf_min=maf_min,
            dp_cap=dp_cap,
            mapq=mapq,
            depth_tool=depth_tool.value,
            require_pass=require_pass,
            dry_run=dry_run,
            to_stdout=stdout,
            emit_vcf=emit_vcf,
            emit_bed=emit_bed,
            emit_matrix=emit_matrix,
            emit_per_contig=emit_per_contig,
            emit_multiqc=emit_multiqc,
            emit_provenance=emit_provenance,
            tsv_mode=tsv_mode,
        )
        
        # Execute plan
        result, provenance_record = execute_plan(plan)
        
        # Handle provenance output
        if emit_provenance and provenance_record:
            from .provenance import write_provenance_json
            out_dir = output_dir or Path.cwd()
            provenance_path = out_dir / "run_provenance.json"
            write_provenance_json([provenance_record], provenance_path)
            console.print(f"Provenance written to {provenance_path}")
        
        if not stdout:
            console.print(f"[green]Self-mapping completed for sample {result.sample}[/green]")
            total_sites = result.ambiguous_snv_count + result.ambiguous_indel_count + result.ambiguous_del_count
            console.print(f"Found {total_sites} ambiguous sites")
        
    except Exception as e:
        console.print(f"[red]Error: {e}[/red]")
        raise typer.Exit(1)


@app.command()
def ref(
    ctx: typer.Context,
    r1: Annotated[
        Path,
        typer.Option(
            "--r1", help="Forward reads (FASTQ, can be gzipped)"
        ),
    ],
    r2: Annotated[
        Path,
        typer.Option(
            "--r2", help="Reverse reads (FASTQ, can be gzipped)"
        ),
    ],
    reference: Annotated[
        Optional[Path],
        typer.Option(
            "--reference", help="Reference genome FASTA file"
        ),
    ] = None,
    species: Annotated[
        Optional[str],
        typer.Option(
            "--species", help="Species name for reference lookup"
        ),
    ] = None,
    bracken: Annotated[
        Optional[Path],
        typer.Option(
            "--bracken", help="Bracken classification file"
        ),
    ] = None,
    sample: Annotated[
        Optional[str],
        typer.Option(
            "--sample", help="Sample name (inferred from filenames if not provided)"
        ),
    ] = None,
    output_dir: Annotated[
        Optional[Path],
        typer.Option(
            "--outdir", "-o", help="Output directory (default: current)"
        ),
    ] = None,
    threads: Annotated[
        int,
        typer.Option(
            "--threads", "-t", help="Number of threads"
        ),
    ] = 4,
    mapper: Annotated[
        str,
        typer.Option(
            "--mapper", help="Mapper to use"
        ),
    ] = "minimap2",
    caller: Annotated[
        str,
        typer.Option(
            "--caller", help="Variant caller to use"
        ),
    ] = "bbtools",
    bbtools_mem: Annotated[
        Optional[str],
        typer.Option(
            "--bbtools-mem", help="BBTools heap memory (e.g., '4g', '8g')"
        ),
    ] = None,
    dp_min: Annotated[
        int,
        typer.Option(
            "--dp-min", help="Minimum depth for ambiguous sites"
        ),
    ] = 10,
    maf_min: Annotated[
        float,
        typer.Option(
            "--maf-min", help="Minimum minor allele frequency"
        ),
    ] = 0.1,
    mapq: Annotated[
        int,
        typer.Option(
            "--mapq", help="Minimum mapping quality for depth analysis"
        ),
    ] = 30,
    dp_cap: Annotated[
        int,
        typer.Option(
            "--dp-cap", help="Maximum depth cap for variant analysis"
        ),
    ] = 100,
    depth_tool: Annotated[
        DepthTool,
        typer.Option(
            "--depth-tool", help="Tool for depth analysis"
        ),
    ] = DepthTool.MOSDEPTH,
    require_pass: Annotated[
        bool,
        typer.Option(
            "--require-pass", help="Only consider variants that pass caller filters"
        ),
    ] = False,
    min_bracken_frac: Annotated[
        float,
        typer.Option(
            "--min-bracken-frac", help="Minimum Bracken abundance fraction"
        ),
    ] = 0.70,
    min_bracken_reads: Annotated[
        int,
        typer.Option(
            "--min-bracken-reads", help="Minimum Bracken read count"
        ),
    ] = 100000,
    ref_dir: Annotated[
        Optional[Path],
        typer.Option(
            "--ref-dir", help="Admin reference directory"
        ),
    ] = None,
    on_fail: Annotated[
        str,
        typer.Option(
            "--on-fail", help="Action when reference resolution fails"
        ),
    ] = "error",
    emit_vcf: Annotated[
        bool,
        typer.Option(
            "--emit-vcf", help="Emit VCF file with ambiguous sites"
        ),
    ] = False,
    emit_bed: Annotated[
        bool,
        typer.Option(
            "--emit-bed", help="Emit BED file with ambiguous sites"
        ),
    ] = False,
    emit_matrix: Annotated[
        bool,
        typer.Option(
            "--emit-matrix", help="Emit variant matrix in TSV format"
        ),
    ] = False,
    emit_per_contig: Annotated[
        bool,
        typer.Option(
            "--emit-per-contig", help="Emit per-contig summary statistics"
        ),
    ] = False,
    emit_multiqc: Annotated[
        bool,
        typer.Option(
            "--emit-multiqc", help="Emit MultiQC-compatible metrics"
        ),
    ] = False,
    emit_provenance: Annotated[
        bool,
        typer.Option(
            "--emit-provenance", help="Emit JSON provenance file"
        ),
    ] = False,
    tsv_mode: Annotated[
        TSVMode,
        typer.Option(
            "--tsv-mode", help="TSV output mode"
        ),
    ] = TSVMode.OVERWRITE,
    stdout: Annotated[
        bool,
        typer.Option(
            "--stdout", help="Write summary to stdout instead of file"
        ),
    ] = False,
) -> None:
    """
    Reference-mapping mode: map reads to a reference genome.
    
    Maps paired-end reads to a reference genome to identify ambiguous sites.
    Reference can be provided directly, looked up by species, or selected from Bracken.
    """
    try:
        # Validate emit flags with stdout
        if stdout and any([emit_vcf, emit_bed, emit_matrix, emit_per_contig, emit_multiqc, emit_provenance]):
            console.print("[red]Error: --stdout cannot be used with --emit-* flags[/red]")
            raise typer.Exit(1)
        
        # Get global options from context
        dry_run = ctx.parent.params.get("dry_run", False) if ctx.parent else False
        
        # Create execution plan
        plan = create_run_plan(
            mode=Mode.REF,
            r1=r1,
            r2=r2,
            reference=reference,
            sample=sample,
            output_dir=output_dir,
            threads=threads,
            mapper=mapper,
            caller=caller,
            bbtools_mem=bbtools_mem,
            dp_min=dp_min,
            maf_min=maf_min,
            dp_cap=dp_cap,
            mapq=mapq,
            depth_tool=depth_tool.value,
            require_pass=require_pass,
            dry_run=dry_run,
            to_stdout=stdout,
            emit_vcf=emit_vcf,
            emit_bed=emit_bed,
            emit_matrix=emit_matrix,
            emit_per_contig=emit_per_contig,
            emit_multiqc=emit_multiqc,
            emit_provenance=emit_provenance,
            species=species,
            bracken=bracken,
            ref_dir=ref_dir,
            on_fail=on_fail,
            tsv_mode=tsv_mode,
            min_bracken_frac=min_bracken_frac,
            min_bracken_reads=min_bracken_reads,
        )
        
        # Execute plan
        result, provenance_record = execute_plan(plan, 
                            species=species, 
                            bracken=bracken, 
                            ref_dir=ref_dir, 
                            on_fail=on_fail,
                            min_bracken_frac=min_bracken_frac,
                            min_bracken_reads=min_bracken_reads)
        
        # Handle provenance output
        if emit_provenance and provenance_record:
            from .provenance import write_provenance_json
            out_dir = output_dir or Path.cwd()
            provenance_path = out_dir / "run_provenance.json"
            write_provenance_json([provenance_record], provenance_path)
            console.print(f"Provenance written to {provenance_path}")
        
        if not stdout:
            console.print(f"[green]Reference-mapping completed for sample {result.sample}[/green]")
            total_sites = result.ambiguous_snv_count + result.ambiguous_indel_count + result.ambiguous_del_count
            console.print(f"Found {total_sites} ambiguous sites")
        
    except Exception as e:
        console.print(f"[red]Error: {e}[/red]")
        raise typer.Exit(1)


@app.command()
def summarize(
    vcf: Annotated[
        Path,
        typer.Option(
            "--vcf", help="VCF file to summarize"
        ),
    ],
    bam: Annotated[
        Path,
        typer.Option(
            "--bam", help="BAM file for denominator calculation"
        ),
    ],
    output: Annotated[
        Optional[Path],
        typer.Option(
            "--output", "-o", help="Output summary file"
        ),
    ] = None,
    dp_min: Annotated[
        int,
        typer.Option(
            "--dp-min", help="Minimum depth for ambiguous sites"
        ),
    ] = 10,
    maf_min: Annotated[
        float,
        typer.Option(
            "--maf-min", help="Minimum minor allele frequency"
        ),
    ] = 0.1,
    dp_cap: Annotated[
        int,
        typer.Option(
            "--dp-cap", help="Maximum depth cap for variant analysis"
        ),
    ] = 100,
    require_pass: Annotated[
        bool,
        typer.Option(
            "--require-pass", help="Only consider variants that pass caller filters"
        ),
    ] = False,
    emit_vcf: Annotated[
        bool,
        typer.Option(
            "--emit-vcf", help="Emit VCF file with ambiguous sites"
        ),
    ] = False,
    emit_bed: Annotated[
        bool,
        typer.Option(
            "--emit-bed", help="Emit BED file with ambiguous sites"
        ),
    ] = False,
    emit_matrix: Annotated[
        bool,
        typer.Option(
            "--emit-matrix", help="Emit variant matrix in TSV format"
        ),
    ] = False,
    emit_per_contig: Annotated[
        bool,
        typer.Option(
            "--emit-per-contig", help="Emit per-contig summary statistics"
        ),
    ] = False,
    emit_multiqc: Annotated[
        bool,
        typer.Option(
            "--emit-multiqc", help="Emit MultiQC-compatible metrics"
        ),
    ] = False,
    emit_provenance: Annotated[
        bool,
        typer.Option(
            "--emit-provenance", help="Emit JSON provenance file"
        ),
    ] = False,
    mode: Annotated[
        str,
        typer.Option(
            "--mode", help="Summary mode"
        ),
    ] = "combined",
    stdout: Annotated[
        bool,
        typer.Option(
            "--stdout", help="Write summary to stdout instead of file"
        ),
    ] = False,
) -> None:
    """
    Summarize VCF and BAM files to generate ambiguous site summary.
    
    Analyzes a VCF file with BAM for denominator calculation to produce
    ambiguous site statistics.
    """
    try:
        console.print("[green]Running summarize mode...[/green]")
        console.print(f"VCF file: {vcf}")
        console.print(f"BAM file: {bam}")
        console.print(f"Output: {output or 'stdout' if stdout else 'auto'}")
        console.print(f"Mode: {mode}")
        console.print(f"Thresholds: dp_min={dp_min}, maf_min={maf_min}")
        
        # Call the summarize logic
        results = run_summarize(
            vcf=vcf,
            bam=bam,
            output=output,
            dp_min=dp_min,
            maf_min=maf_min,
            dp_cap=dp_cap,
            require_pass=require_pass,
            emit_vcf=emit_vcf,
            emit_bed=emit_bed,
            emit_matrix=emit_matrix,
            emit_per_contig=emit_per_contig,
            emit_multiqc=emit_multiqc,
            emit_provenance=emit_provenance,
            to_stdout=stdout,
        )
        
        if not stdout:
            console.print(f"[green]Summarize completed with {len(results)} results[/green]")
        
    except Exception as e:
        console.print(f"[red]Error: {e}[/red]")
        raise typer.Exit(1)


if __name__ == "__main__":
    app()