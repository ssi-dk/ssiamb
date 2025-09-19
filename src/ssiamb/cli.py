#!/usr/bin/env python3
"""
ssiamb - SSI Ambiguous Site Detection Tool

Command-line interface for detecting ambiguous sites in bacterial genomes
through mapping and variant calling approaches.
"""

from pathlib import Path
from typing import Optional, Annotated
import sys
import os
import typer
from rich.console import Console

from .version import __version__
from .models import Mode, TSVMode, DepthTool
from .runner import create_run_plan, execute_plan, run_summarize
from .errors import handle_exception_with_exit

app = typer.Typer(
    name="ssiamb",
    help="SSI Ambiguous Site Detection Tool",
    add_completion=False,
    no_args_is_help=True,
)


# Configure console for better test compatibility
def _is_test_environment() -> bool:
    """Detect if we're running in a test environment."""
    return (
        "pytest" in sys.modules
        or os.getenv("PYTEST_CURRENT_TEST") is not None
        or os.getenv("CI") is not None  # CI environments often set CI=true or CI=1
        or os.getenv("GITHUB_ACTIONS") is not None
        or os.getenv("RUNNER_OS") is not None  # GitHub Actions specific
        or "unittest" in sys.modules
    )


# Set NO_COLOR environment variable for test environments
if _is_test_environment():
    os.environ["NO_COLOR"] = "1"
    os.environ["FORCE_COLOR"] = "0"

console = Console(
    force_terminal=not _is_test_environment(),
    no_color=_is_test_environment(),  # Completely disable colors in test environments
    width=80 if _is_test_environment() else None,  # Fixed width for consistent output
    legacy_windows=False,  # Disable legacy Windows console behavior
    _environ={"NO_COLOR": "1"} if _is_test_environment() else None,  # Force NO_COLOR
)


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
        typer.Option("--verbose", "-v", help="Enable verbose output"),
    ] = False,
    quiet: Annotated[
        bool,
        typer.Option("--quiet", "-q", help="Suppress non-error output"),
    ] = False,
    dry_run: Annotated[
        bool,
        typer.Option("--dry-run", help="Show what would be done without executing"),
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
        typer.Option("--r1", help="Forward reads (FASTQ, can be gzipped)"),
    ],
    r2: Annotated[
        Path,
        typer.Option("--r2", help="Reverse reads (FASTQ, can be gzipped)"),
    ],
    assembly: Annotated[
        Path,
        typer.Option("--assembly", help="Assembly FASTA file"),
    ],
    sample: Annotated[
        Optional[str],
        typer.Option(
            "--sample", help="Sample name (inferred from filenames if not provided)"
        ),
    ] = None,
    output_dir: Annotated[
        Optional[Path],
        typer.Option("--outdir", "-o", help="Output directory (default: current)"),
    ] = None,
    threads: Annotated[
        int,
        typer.Option("--threads", "-t", help="Number of threads"),
    ] = 4,
    mapper: Annotated[
        str,
        typer.Option("--mapper", help="Mapper to use"),
    ] = "minimap2",
    caller: Annotated[
        str,
        typer.Option("--caller", help="Variant caller to use"),
    ] = "bbtools",
    bbtools_mem: Annotated[
        Optional[str],
        typer.Option("--bbtools-mem", help="BBTools heap memory (e.g., '4g', '8g')"),
    ] = None,
    dp_min: Annotated[
        int,
        typer.Option("--dp-min", help="Minimum depth for ambiguous sites"),
    ] = 10,
    maf_min: Annotated[
        float,
        typer.Option("--maf-min", help="Minimum minor allele frequency"),
    ] = 0.1,
    mapq: Annotated[
        int,
        typer.Option("--mapq", help="Minimum mapping quality for depth analysis"),
    ] = 30,
    dp_cap: Annotated[
        int,
        typer.Option("--dp-cap", help="Maximum depth cap for variant analysis"),
    ] = 100,
    depth_tool: Annotated[
        DepthTool,
        typer.Option("--depth-tool", help="Tool for depth analysis"),
    ] = DepthTool.MOSDEPTH,
    require_pass: Annotated[
        bool,
        typer.Option(
            "--require-pass", help="Only consider variants that pass caller filters"
        ),
    ] = False,
    emit_vcf: Annotated[
        bool,
        typer.Option("--emit-vcf", help="Emit VCF file with ambiguous sites"),
    ] = False,
    emit_bed: Annotated[
        bool,
        typer.Option("--emit-bed", help="Emit BED file with ambiguous sites"),
    ] = False,
    emit_matrix: Annotated[
        bool,
        typer.Option("--emit-matrix", help="Emit variant matrix in TSV format"),
    ] = False,
    emit_per_contig: Annotated[
        bool,
        typer.Option("--emit-per-contig", help="Emit per-contig summary statistics"),
    ] = False,
    emit_multiqc: Annotated[
        bool,
        typer.Option("--emit-multiqc", help="Emit MultiQC-compatible metrics"),
    ] = False,
    emit_provenance: Annotated[
        bool,
        typer.Option("--emit-provenance", help="Emit JSON provenance file"),
    ] = False,
    tsv_mode: Annotated[
        TSVMode,
        typer.Option("--tsv-mode", help="TSV output mode"),
    ] = TSVMode.OVERWRITE,
    stdout: Annotated[
        bool,
        typer.Option("--stdout", help="Write summary to stdout instead of file"),
    ] = False,
) -> None:
    """
    Self-mapping mode: map reads to their own assembly.

    Maps paired-end reads to the provided assembly to identify ambiguous sites
    where the assembly may not represent the true sequence.
    """
    try:
        # Validate required parameters
        if not r1:
            console.print("[red]Error: Missing option '--r1'[/red]")
            raise typer.Exit(2)
        if not r2:
            console.print("[red]Error: Missing option '--r2'[/red]")
            raise typer.Exit(2)
        if not assembly:
            console.print("[red]Error: Missing option '--assembly'[/red]")
            raise typer.Exit(2)

        # Validate emit flags with stdout
        if stdout and any(
            [
                emit_vcf,
                emit_bed,
                emit_matrix,
                emit_per_contig,
                emit_multiqc,
                emit_provenance,
            ]
        ):
            console.print(
                "[red]Error: --stdout cannot be used with --emit-* flags[/red]"
            )
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

        # If this was a dry run, we're done
        if dry_run:
            return

        # Handle provenance output
        if emit_provenance and provenance_record:
            from .provenance import write_provenance_json

            out_dir = output_dir or Path.cwd()
            provenance_path = out_dir / "run_provenance.json"
            write_provenance_json([provenance_record], provenance_path)
            console.print(f"Provenance written to {provenance_path}")

        if not stdout:
            console.print(
                f"[green]Self-mapping completed for sample {result.sample}[/green]"
            )
            total_sites = (
                result.ambiguous_snv_count
                + result.ambiguous_indel_count
                + result.ambiguous_del_count
            )
            console.print(f"Found {total_sites} ambiguous sites")

    except Exception as e:
        handle_exception_with_exit(e, "Self-mapping mode failed")


@app.command()
def ref(
    ctx: typer.Context,
    r1: Annotated[
        Path,
        typer.Option("--r1", help="Forward reads (FASTQ, can be gzipped)"),
    ],
    r2: Annotated[
        Path,
        typer.Option("--r2", help="Reverse reads (FASTQ, can be gzipped)"),
    ],
    reference: Annotated[
        Optional[Path],
        typer.Option("--reference", help="Reference genome FASTA file"),
    ] = None,
    species: Annotated[
        Optional[str],
        typer.Option("--species", help="Species name for reference lookup"),
    ] = None,
    bracken: Annotated[
        Optional[Path],
        typer.Option("--bracken", help="Bracken classification file"),
    ] = None,
    sample: Annotated[
        Optional[str],
        typer.Option(
            "--sample", help="Sample name (inferred from filenames if not provided)"
        ),
    ] = None,
    output_dir: Annotated[
        Optional[Path],
        typer.Option("--outdir", "-o", help="Output directory (default: current)"),
    ] = None,
    threads: Annotated[
        int,
        typer.Option("--threads", "-t", help="Number of threads"),
    ] = 4,
    mapper: Annotated[
        str,
        typer.Option("--mapper", help="Mapper to use"),
    ] = "minimap2",
    caller: Annotated[
        str,
        typer.Option("--caller", help="Variant caller to use"),
    ] = "bbtools",
    bbtools_mem: Annotated[
        Optional[str],
        typer.Option("--bbtools-mem", help="BBTools heap memory (e.g., '4g', '8g')"),
    ] = None,
    dp_min: Annotated[
        int,
        typer.Option("--dp-min", help="Minimum depth for ambiguous sites"),
    ] = 10,
    maf_min: Annotated[
        float,
        typer.Option("--maf-min", help="Minimum minor allele frequency"),
    ] = 0.1,
    mapq: Annotated[
        int,
        typer.Option("--mapq", help="Minimum mapping quality for depth analysis"),
    ] = 30,
    dp_cap: Annotated[
        int,
        typer.Option("--dp-cap", help="Maximum depth cap for variant analysis"),
    ] = 100,
    depth_tool: Annotated[
        DepthTool,
        typer.Option("--depth-tool", help="Tool for depth analysis"),
    ] = DepthTool.MOSDEPTH,
    require_pass: Annotated[
        bool,
        typer.Option(
            "--require-pass", help="Only consider variants that pass caller filters"
        ),
    ] = False,
    min_bracken_frac: Annotated[
        float,
        typer.Option("--min-bracken-frac", help="Minimum Bracken abundance fraction"),
    ] = 0.70,
    min_bracken_reads: Annotated[
        int,
        typer.Option("--min-bracken-reads", help="Minimum Bracken read count"),
    ] = 100000,
    ref_dir: Annotated[
        Optional[Path],
        typer.Option("--ref-dir", help="Admin reference directory"),
    ] = None,
    on_fail: Annotated[
        str,
        typer.Option("--on-fail", help="Action when reference resolution fails"),
    ] = "error",
    emit_vcf: Annotated[
        bool,
        typer.Option("--emit-vcf", help="Emit VCF file with ambiguous sites"),
    ] = False,
    emit_bed: Annotated[
        bool,
        typer.Option("--emit-bed", help="Emit BED file with ambiguous sites"),
    ] = False,
    emit_matrix: Annotated[
        bool,
        typer.Option("--emit-matrix", help="Emit variant matrix in TSV format"),
    ] = False,
    emit_per_contig: Annotated[
        bool,
        typer.Option("--emit-per-contig", help="Emit per-contig summary statistics"),
    ] = False,
    emit_multiqc: Annotated[
        bool,
        typer.Option("--emit-multiqc", help="Emit MultiQC-compatible metrics"),
    ] = False,
    emit_provenance: Annotated[
        bool,
        typer.Option("--emit-provenance", help="Emit JSON provenance file"),
    ] = False,
    tsv_mode: Annotated[
        TSVMode,
        typer.Option("--tsv-mode", help="TSV output mode"),
    ] = TSVMode.OVERWRITE,
    stdout: Annotated[
        bool,
        typer.Option("--stdout", help="Write summary to stdout instead of file"),
    ] = False,
) -> None:
    """
    Reference-mapping mode: map reads to a reference genome.

    Maps paired-end reads to a reference genome to identify ambiguous sites.
    Reference can be provided directly, looked up by species, or selected from Bracken.
    """
    try:
        # Validate required parameters
        if not r1:
            console.print("[red]Error: Missing option '--r1'[/red]")
            raise typer.Exit(2)
        if not r2:
            console.print("[red]Error: Missing option '--r2'[/red]")
            raise typer.Exit(2)

        # At least one of reference, species, or bracken must be provided
        if not any([reference, species, bracken]):
            console.print(
                "[red]Error: Must provide one of --reference, --species, or --bracken[/red]"
            )
            raise typer.Exit(2)

        # Validate emit flags with stdout
        if stdout and any(
            [
                emit_vcf,
                emit_bed,
                emit_matrix,
                emit_per_contig,
                emit_multiqc,
                emit_provenance,
            ]
        ):
            console.print(
                "[red]Error: --stdout cannot be used with --emit-* flags[/red]"
            )
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
        result, provenance_record = execute_plan(
            plan,
            species=species,
            bracken=bracken,
            ref_dir=ref_dir,
            on_fail=on_fail,
            min_bracken_frac=min_bracken_frac,
            min_bracken_reads=min_bracken_reads,
        )

        # If this was a dry run, we're done
        if dry_run:
            return

        # Handle provenance output
        if emit_provenance and provenance_record:
            from .provenance import write_provenance_json

            out_dir = output_dir or Path.cwd()
            provenance_path = out_dir / "run_provenance.json"
            write_provenance_json([provenance_record], provenance_path)
            console.print(f"Provenance written to {provenance_path}")

        if not stdout:
            console.print(
                f"[green]Reference-mapping completed for sample {result.sample}[/green]"
            )
            total_sites = (
                result.ambiguous_snv_count
                + result.ambiguous_indel_count
                + result.ambiguous_del_count
            )
            console.print(f"Found {total_sites} ambiguous sites")

    except Exception as e:
        handle_exception_with_exit(e, "Reference-mapping mode failed")


@app.command()
def summarize(
    ctx: typer.Context,
    vcf: Annotated[
        Path,
        typer.Option("--vcf", help="VCF file to summarize"),
    ],
    bam: Annotated[
        Path,
        typer.Option("--bam", help="BAM file for denominator calculation"),
    ],
    output: Annotated[
        Optional[Path],
        typer.Option("--output", "-o", help="Output summary file"),
    ] = None,
    dp_min: Annotated[
        int,
        typer.Option("--dp-min", help="Minimum depth for ambiguous sites"),
    ] = 10,
    maf_min: Annotated[
        float,
        typer.Option("--maf-min", help="Minimum minor allele frequency"),
    ] = 0.1,
    dp_cap: Annotated[
        int,
        typer.Option("--dp-cap", help="Maximum depth cap for variant analysis"),
    ] = 100,
    require_pass: Annotated[
        bool,
        typer.Option(
            "--require-pass", help="Only consider variants that pass caller filters"
        ),
    ] = False,
    emit_vcf: Annotated[
        bool,
        typer.Option("--emit-vcf", help="Emit VCF file with ambiguous sites"),
    ] = False,
    emit_bed: Annotated[
        bool,
        typer.Option("--emit-bed", help="Emit BED file with ambiguous sites"),
    ] = False,
    emit_matrix: Annotated[
        bool,
        typer.Option("--emit-matrix", help="Emit variant matrix in TSV format"),
    ] = False,
    emit_per_contig: Annotated[
        bool,
        typer.Option("--emit-per-contig", help="Emit per-contig summary statistics"),
    ] = False,
    emit_multiqc: Annotated[
        bool,
        typer.Option("--emit-multiqc", help="Emit MultiQC-compatible metrics"),
    ] = False,
    emit_provenance: Annotated[
        bool,
        typer.Option("--emit-provenance", help="Emit JSON provenance file"),
    ] = False,
    mode: Annotated[
        str,
        typer.Option("--mode", help="Summary mode"),
    ] = "combined",
    stdout: Annotated[
        bool,
        typer.Option("--stdout", help="Write summary to stdout instead of file"),
    ] = False,
) -> None:
    """
    Summarize VCF and BAM files to generate ambiguous site summary.

    Analyzes a VCF file with BAM for denominator calculation to produce
    ambiguous site statistics.
    """
    try:
        # Validate required parameters
        if not vcf:
            console.print("[red]Error: Missing option '--vcf'[/red]")
            raise typer.Exit(2)
        if not bam:
            console.print("[red]Error: Missing option '--bam'[/red]")
            raise typer.Exit(2)

        # Get global options from context
        dry_run = ctx.parent.params.get("dry_run", False) if ctx.parent else False

        if dry_run:
            console.print("[yellow]DRY RUN - Summarize mode plan:[/yellow]")
            console.print("  [bold]Input validation:[/bold]")
            console.print(
                f"    VCF file: {vcf} {'✓' if vcf.exists() else '✗ (missing)'}"
            )
            console.print(
                f"    BAM file: {bam} {'✓' if bam.exists() else '✗ (missing)'}"
            )
            console.print(
                f"  [bold]Sample name:[/bold] {vcf.stem.split('.')[0] if vcf.exists() else 'unknown'}"
            )
            console.print("  [bold]Analysis plan:[/bold]")
            console.print(
                "    1. Run mosdepth on BAM (MAPQ≥30, depth≥10, exclude duplicates)"
            )
            console.print("    2. Parse VCF for ambiguous sites")
            console.print(
                f"    3. Count SNVs: dp_min={dp_min}, maf_min={maf_min}, dp_cap={dp_cap}"
            )
            console.print("    4. Count indels and deletions for secondary metrics")
            console.print("    5. Calculate mapping rate from BAM")
            console.print("  [bold]Outputs planned:[/bold]")
            if stdout:
                console.print("    Summary: stdout")
            else:
                console.print(f"    Summary: {output or 'ambiguous_summary.tsv'}")
            if emit_vcf:
                console.print(f"    VCF: {vcf.stem}.ambiguous_sites.vcf.gz")
            if emit_bed:
                console.print(f"    BED: {vcf.stem}.ambiguous_sites.bed.gz")
            if emit_matrix:
                console.print(f"    Matrix: {vcf.stem}.variant_matrix.tsv.gz")
            if emit_per_contig:
                console.print(f"    Per-contig: {vcf.stem}.per_contig_summary.tsv")
            if emit_multiqc:
                console.print(f"    MultiQC: {vcf.stem}.multiqc.tsv")
            if emit_provenance:
                console.print("    Provenance: run_provenance.json")
            console.print(
                f"  [bold]Filters:[/bold] {'PASS-only' if require_pass else 'All variants'}"
            )
            console.print("[green]Dry run completed - no files written[/green]")
            return

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
            console.print(
                f"[green]Summarize completed with {len(results)} results[/green]"
            )

    except Exception as e:
        handle_exception_with_exit(e, "Summarize mode failed")


if __name__ == "__main__":
    app()
