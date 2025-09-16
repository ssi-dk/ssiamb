"""
Main orchestration logic for ssiamb commands.

This module contains the high-level functions that coordinate the execution
of self, ref, and summarize modes by calling appropriate modules.
"""

from pathlib import Path
from typing import Optional, List
import logging

from .models import (
    Mode, Mapper, Caller, DepthTool, SummaryRow, RunPlan, Paths, Thresholds,
    TSVMode, OnFail
)
from .io_utils import (
    infer_sample_name, validate_sample_name, write_tsv_summary, 
    write_tsv_to_stdout, SampleNameError
)
from .version import __version__

logger = logging.getLogger(__name__)


def create_run_plan(
    mode: Mode,
    r1: Path,
    r2: Path,
    assembly: Optional[Path] = None,
    reference: Optional[Path] = None,
    sample: Optional[str] = None,
    output_dir: Optional[Path] = None,
    threads: int = 1,
    mapper: str = "minimap2",
    caller: str = "bbtools",
    dp_min: int = 10,
    maf_min: float = 0.1,
    dry_run: bool = False,
    to_stdout: bool = False,
    **kwargs,
) -> RunPlan:
    """
    Create a validated execution plan from CLI arguments.
    
    Args:
        mode: Analysis mode
        r1: Forward reads file
        r2: Reverse reads file
        assembly: Assembly file (for self mode)
        reference: Reference file (for ref mode)
        sample: Sample name (inferred if not provided)
        output_dir: Output directory
        threads: Number of threads
        mapper: Mapper tool name
        caller: Caller tool name
        dp_min: Minimum depth threshold
        maf_min: Minimum MAF threshold
        dry_run: Dry run mode
        to_stdout: Write to stdout
        **kwargs: Additional arguments
        
    Returns:
        Validated RunPlan object
        
    Raises:
        ValueError: If configuration is invalid
        SampleNameError: If sample name is invalid
        FileNotFoundError: If required files don't exist
    """
    # Validate and resolve sample name
    if sample is None:
        sample = infer_sample_name(r1, r2)
    else:
        sample = validate_sample_name(sample)
    
    # Create paths object
    paths = Paths(
        r1=r1,
        r2=r2,
        assembly=assembly,
        reference=reference,
        output_dir=output_dir or Path("."),
        sample=sample,
    )
    
    # Create thresholds
    thresholds = Thresholds(
        dp_min=dp_min,
        maf_min=maf_min,
    )
    
    # Validate mapper and caller
    try:
        mapper_enum = Mapper(mapper)
    except ValueError:
        valid_mappers = [m.value for m in Mapper]
        raise ValueError(f"Invalid mapper '{mapper}'. Valid options: {valid_mappers}")
    
    try:
        caller_enum = Caller(caller)
    except ValueError:
        valid_callers = [c.value for c in Caller]
        raise ValueError(f"Invalid caller '{caller}'. Valid options: {valid_callers}")
    
    # Mode-specific validation
    if mode == Mode.SELF and assembly is None:
        raise ValueError("Self mode requires --assembly")
    
    if mode == Mode.REF and all(x is None for x in [reference, kwargs.get("species"), kwargs.get("bracken")]):
        raise ValueError("Ref mode requires one of: --reference, --species, or --bracken")
    
    # Create run plan
    plan = RunPlan(
        mode=mode,
        sample=sample,
        paths=paths,
        thresholds=thresholds,
        mapper=mapper_enum,
        caller=caller_enum,
        threads=threads,
        dry_run=dry_run,
        to_stdout=to_stdout,
        # Extract output flags from kwargs
        emit_vcf=kwargs.get("emit_vcf", False),
        emit_bed=kwargs.get("emit_bed", False),
        emit_matrix=kwargs.get("emit_matrix", False),
        emit_per_contig=kwargs.get("emit_per_contig", False),
        emit_provenance=kwargs.get("emit_provenance", False),
        emit_multiqc=kwargs.get("emit_multiqc", False),
    )
    
    logger.info(f"Created run plan for {mode.value} mode, sample {sample}")
    return plan


def run_self(plan: RunPlan) -> SummaryRow:
    """
    Execute self-mapping mode.
    
    Args:
        plan: Validated execution plan
        
    Returns:
        Summary row with results
    """
    logger.info(f"Running self mode for sample {plan.sample}")
    
    if plan.dry_run:
        logger.info("DRY RUN - would execute self-mapping pipeline")
        logger.info(f"  Map {plan.paths.r1} + {plan.paths.r2} to {plan.paths.assembly}")
        logger.info(f"  Using {plan.mapper.value} mapper and {plan.caller.value} caller")
        logger.info(f"  Thresholds: dp_min={plan.thresholds.dp_min}, maf_min={plan.thresholds.maf_min}")
    
    # TODO: Implement actual self-mapping pipeline
    # For now, return placeholder data
    return SummaryRow(
        sample=plan.sample,
        mode=plan.mode.value,
        ref_label=plan.paths.assembly.name if plan.paths.assembly else "unknown",
        mapper=plan.mapper.value,
        caller=plan.caller.value,
        dp_min=plan.thresholds.dp_min,
        maf_min=plan.thresholds.maf_min,
        callable_bases=2800000,  # ~2.8Mb typical bacterial genome
        genome_length=2900000,
        breadth_10x=0.95,
        mean_depth=25.0,
        mapping_rate=0.98,
        ambiguous_sites=42,
        ref_calls=2799950,
        alt_calls=8,
        no_call=50,
        qc_warnings="",
    )


def run_ref(plan: RunPlan, **kwargs) -> SummaryRow:
    """
    Execute reference-mapping mode.
    
    Args:
        plan: Validated execution plan
        **kwargs: Additional arguments (species, bracken, ref_dir, on_fail)
        
    Returns:
        Summary row with results
    """
    logger.info(f"Running ref mode for sample {plan.sample}")
    
    # TODO: Resolve reference from precedence: file > species > bracken
    ref_source = "unknown"
    ref_label = "unknown"
    
    if plan.paths.reference:
        ref_source = "file"
        ref_label = plan.paths.reference.name
    elif kwargs.get("species"):
        ref_source = "species"
        ref_label = kwargs["species"]
    elif kwargs.get("bracken"):
        ref_source = "bracken"
        ref_label = "bracken-selected"
    
    if plan.dry_run:
        logger.info("DRY RUN - would execute reference-mapping pipeline")
        logger.info(f"  Resolve reference from {ref_source}")
        logger.info(f"  Map {plan.paths.r1} + {plan.paths.r2} to reference")
        logger.info(f"  Using {plan.mapper.value} mapper and {plan.caller.value} caller")
        logger.info(f"  Thresholds: dp_min={plan.thresholds.dp_min}, maf_min={plan.thresholds.maf_min}")
    
    # TODO: Implement actual ref-mapping pipeline
    # For now, return placeholder data
    return SummaryRow(
        sample=plan.sample,
        mode=plan.mode.value,
        ref_label=ref_label,
        mapper=plan.mapper.value,
        caller=plan.caller.value,
        dp_min=plan.thresholds.dp_min,
        maf_min=plan.thresholds.maf_min,
        callable_bases=2750000,  # Reference genome size
        genome_length=2800000,
        breadth_10x=0.92,
        mean_depth=23.5,
        mapping_rate=0.95,
        ambiguous_sites=67,
        ref_calls=2749920,
        alt_calls=13,
        no_call=80,
        qc_warnings="low_mapping_rate",
    )


def run_summarize(
    input_files: List[Path],
    output: Optional[Path] = None,
    mode: str = "combined",
    to_stdout: bool = False,
) -> List[SummaryRow]:
    """
    Execute summarize mode.
    
    Args:
        input_files: List of TSV files to summarize
        output: Output file (auto-generated if not provided)
        mode: Summary mode
        to_stdout: Write to stdout
        
    Returns:
        List of summary rows
    """
    logger.info(f"Running summarize mode on {len(input_files)} files")
    
    # TODO: Implement actual summarization logic
    # For now, return empty list
    return []


def execute_plan(plan: RunPlan, **kwargs) -> SummaryRow:
    """
    Execute a run plan and handle output.
    
    Args:
        plan: Validated execution plan
        **kwargs: Mode-specific arguments
        
    Returns:
        Summary row with results
    """
    # Execute the appropriate mode
    if plan.mode == Mode.SELF:
        result = run_self(plan)
    elif plan.mode == Mode.REF:
        result = run_ref(plan, **kwargs)
    else:
        raise ValueError(f"Unsupported mode: {plan.mode}")
    
    # Handle output
    if plan.to_stdout:
        write_tsv_to_stdout([result])
    else:
        output_file = plan.paths.output_dir / f"{plan.get_sample_prefix()}_ambiguous_summary.tsv"
        write_tsv_summary(output_file, [result], TSVMode.OVERWRITE)
        logger.info(f"Results written to {output_file}")
    
    return result