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
from .models import Mode
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
            "--output-dir", "-o", help="Output directory (default: current)"
        ),
    ] = None,
    threads: Annotated[
        int,
        typer.Option(
            "--threads", "-t", help="Number of threads"
        ),
    ] = 1,
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
        if stdout and any([emit_vcf, emit_bed, emit_matrix, emit_per_contig, emit_multiqc]):
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
            dp_min=dp_min,
            maf_min=maf_min,
            mapq=mapq,
            dry_run=dry_run,
            to_stdout=stdout,
            emit_vcf=emit_vcf,
            emit_bed=emit_bed,
            emit_matrix=emit_matrix,
            emit_per_contig=emit_per_contig,
            emit_multiqc=emit_multiqc,
        )
        
        # Execute plan
        result = execute_plan(plan)
        
        if not stdout:
            console.print(f"[green]Self-mapping completed for sample {result.sample}[/green]")
            console.print(f"Found {result.ambiguous_sites} ambiguous sites")
        
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
            "--output-dir", "-o", help="Output directory (default: current)"
        ),
    ] = None,
    threads: Annotated[
        int,
        typer.Option(
            "--threads", "-t", help="Number of threads"
        ),
    ] = 1,
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
        if stdout and any([emit_vcf, emit_bed, emit_matrix, emit_per_contig, emit_multiqc]):
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
            dp_min=dp_min,
            maf_min=maf_min,
            mapq=mapq,
            dry_run=dry_run,
            to_stdout=stdout,
            emit_vcf=emit_vcf,
            emit_bed=emit_bed,
            emit_matrix=emit_matrix,
            emit_per_contig=emit_per_contig,
            emit_multiqc=emit_multiqc,
            species=species,
            bracken=bracken,
            ref_dir=ref_dir,
            on_fail=on_fail,
        )
        
        # Execute plan
        result = execute_plan(plan, 
                            species=species, 
                            bracken=bracken, 
                            ref_dir=ref_dir, 
                            on_fail=on_fail)
        
        if not stdout:
            console.print(f"[green]Reference-mapping completed for sample {result.sample}[/green]")
            console.print(f"Found {result.ambiguous_sites} ambiguous sites")
        
    except Exception as e:
        console.print(f"[red]Error: {e}[/red]")
        raise typer.Exit(1)


@app.command()
def summarize(
    input_files: Annotated[
        list[Path],
        typer.Argument(help="Input TSV files to summarize")
    ],
    output: Annotated[
        Optional[Path],
        typer.Option(
            "--output", "-o", help="Output summary file"
        ),
    ] = None,
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
    Summarize multiple ssiamb results.
    
    Combines results from multiple ssiamb runs into a consolidated summary.
    """
    console.print("[green]Running summarize mode...[/green]")
    console.print(f"Input files: {len(input_files)}")
    for f in input_files:
        console.print(f"  - {f}")
    console.print(f"Output: {output or 'stdout' if stdout else 'auto'}")
    console.print(f"Mode: {mode}")
    
    # TODO: Implement summarize logic
    console.print("[yellow]Implementation pending...[/yellow]")


if __name__ == "__main__":
    app()