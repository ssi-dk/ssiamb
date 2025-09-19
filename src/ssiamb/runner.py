"""
Main orchestration logic for ssiamb commands.

This module contains the high-level functions that coordinate the execution
of self, ref, and summarize modes by calling appropriate modules.
"""

from pathlib import Path
from typing import Optional, List, Tuple, Any
import logging
import time
import subprocess
import re

from .models import (
    Mode,
    Mapper,
    Caller,
    DepthTool,
    SummaryRow,
    RunPlan,
    Paths,
    Thresholds,
    TSVMode,
)
from .io_utils import (
    infer_sample_name,
    validate_sample_name,
    write_tsv_summary,
    write_tsv_to_stdout,
)
from .refdir import (
    resolve_reference_directory,
    resolve_species_to_fasta,
    parse_accession_from_fasta_header,
    ReferenceResolutionError,
)
from .mapping import (
    ensure_indexes_self,
    map_fastqs,
    check_external_tools,
    MappingError,
    ExternalToolError,
)
from .depth import analyze_depth, DepthAnalysisError, check_mosdepth_available
from .calling import call_variants, VariantCallingError, caller_tools_available
from .vcf_ops import (
    normalize_and_split,
    count_ambiguous_sites,
    VCFOperationError,
    VariantClass,
    AmbigGrid,
    emit_vcf,
    emit_bed,
    emit_matrix,
    emit_per_contig,
    emit_multiqc,
)
from .mapping import calculate_mapping_rate
from .reuse import has_duplicate_flags, run_markdup_for_depth, CompatibilityError
from .qc import check_qc_metrics, log_qc_warnings
from .provenance import (
    create_provenance_record,
    ProvenanceRecord,
    ProvenanceMappingStats,
    ProvenanceCounts,
)
from .version import __version__
import pysam

logger = logging.getLogger(__name__)


def get_tool_versions(plan: RunPlan) -> str:
    """
    Detect versions of external tools used in the workflow.

    Args:
        plan: Run plan containing tool choices

    Returns:
        Formatted string with tool versions
    """
    versions = {}

    # Map plan tools to version commands
    version_commands = {
        plan.mapper.value: _get_mapper_version(plan.mapper.value),
        plan.caller.value: _get_caller_version(plan.caller.value),
        "samtools": ["samtools", "--version"],
        "mosdepth": ["mosdepth", "--version"],
        "bcftools": ["bcftools", "--version"],
    }

    for tool, cmd in version_commands.items():
        try:
            if cmd:
                result = subprocess.run(cmd, capture_output=True, text=True, timeout=5)
                if result.returncode == 0:
                    version = _extract_version(result.stdout)
                    if version:
                        versions[tool] = version
        except (
            subprocess.TimeoutExpired,
            subprocess.CalledProcessError,
            FileNotFoundError,
        ):
            pass

    # Format as tool:version pairs
    if versions:
        return ";".join(f"{tool}:{ver}" for tool, ver in sorted(versions.items()))
    else:
        return "NA"


def _get_mapper_version(mapper: str) -> Optional[List[str]]:
    """Get version command for mapper."""
    if mapper == "minimap2":
        return ["minimap2", "--version"]
    elif mapper == "bwa-mem2":
        return ["bwa-mem2", "version"]
    return None


def _get_caller_version(caller: str) -> Optional[List[str]]:
    """Get version command for variant caller."""
    if caller == "bbtools":
        return ["bbmap.sh", "version"]
    elif caller == "bcftools":
        return ["bcftools", "--version"]
    return None


def _extract_version(output: str) -> Optional[str]:
    """Extract version number from tool output."""
    # Common version patterns
    patterns = [
        r"version\s+([0-9]+\.[0-9]+(?:\.[0-9]+)?)",
        r"v?([0-9]+\.[0-9]+(?:\.[0-9]+)?)",
        r"([0-9]+\.[0-9]+(?:\.[0-9]+)?)",
    ]

    for pattern in patterns:
        match = re.search(pattern, output, re.IGNORECASE)
        if match:
            return match.group(1)
    return None


def extract_ref_accession(fasta_path: Path) -> str:
    """
    Extract the reference accession from a FASTA file header.

    Args:
        fasta_path: Path to the FASTA file

    Returns:
        Reference accession or "NA" if not found
    """
    try:
        with pysam.FastaFile(str(fasta_path)) as fasta:
            # Get the first reference sequence header
            if len(fasta.references) > 0:
                first_ref = fasta.references[0]
                # Try to extract accession from header format like >NC_123456.1 description
                if first_ref.startswith(("NC_", "NZ_", "AC_", "CP_")):
                    # Extract up to first space or dot
                    accession = first_ref.split(".")[0]
                    return accession
                else:
                    return first_ref[:50]  # First 50 chars of header
    except Exception as e:
        logger.debug(f"Could not extract ref accession from {fasta_path}: {e}")

    return "NA"


def detect_reuse_from_plan(plan: RunPlan) -> Tuple[bool, bool]:
    """
    Detect if BAM or VCF files are being reused based on plan paths.

    Args:
        plan: The execution plan

    Returns:
        Tuple of (reused_bam, reused_vcf) booleans
    """
    reused_bam = False
    reused_vcf = False

    # Check if BAM path is provided (indicating reuse)
    if plan.paths.bam is not None and plan.paths.bam.exists():
        reused_bam = True

    # Check if VCF path is provided (indicating reuse)
    if plan.paths.vcf is not None and plan.paths.vcf.exists():
        reused_vcf = True

    return reused_bam, reused_vcf


def create_summary_row(
    sample: str,
    mode: str,
    mapper: str,
    caller: str,
    dp_min: int,
    maf_min: float,
    dp_cap: int,
    callable_bases: int,
    genome_length: int,
    ambiguous_snv_count: int,
    ambiguous_indel_count: int = 0,
    ambiguous_del_count: int = 0,
    ref_label: str = "NA",
    ref_accession: str = "NA",
    bracken_species: str = "NA",
    bracken_frac: float = 0.0,
    bracken_reads: int = 0,
    alias_used: str = "NA",
    reused_bam: bool = False,
    reused_vcf: bool = False,
    runtime_sec: float = 0.0,
    denom_policy: str = "exclude_dups",
    mapping_rate: Optional[float] = None,
) -> SummaryRow:
    """
    Create a SummaryRow with proper field calculations and formatting.

    Args:
        Basic parameters for SummaryRow creation
        mapping_rate: Optional mapping rate for QC checks

    Returns:
        Fully populated SummaryRow object
    """
    # Calculate derived fields
    breadth_10x = callable_bases / genome_length if genome_length > 0 else 0.0
    ambiguous_snv_per_mb = (
        (ambiguous_snv_count * 1_000_000) / callable_bases
        if callable_bases > 0
        else 0.0
    )

    # Perform QC checks and log warnings
    from .models import Mode

    mode_enum = Mode(mode) if mode in [m.value for m in Mode] else None
    qc_warnings = check_qc_metrics(
        breadth_10x=breadth_10x,
        callable_bases=callable_bases,
        mapping_rate=mapping_rate,
        mode=mode_enum,
    )

    # Log QC warnings
    log_qc_warnings(qc_warnings)

    return SummaryRow(
        sample=sample,
        mode=mode,
        mapper=mapper,
        caller=caller,
        dp_min=dp_min,
        maf_min=maf_min,
        dp_cap=dp_cap,
        denom_policy=denom_policy,
        callable_bases=callable_bases,
        genome_length=genome_length,
        breadth_10x=breadth_10x,
        ambiguous_snv_count=ambiguous_snv_count,
        ambiguous_snv_per_mb=ambiguous_snv_per_mb,
        ambiguous_indel_count=ambiguous_indel_count,
        ambiguous_del_count=ambiguous_del_count,
        ref_label=ref_label,
        ref_accession=ref_accession,
        bracken_species=bracken_species,
        bracken_frac=bracken_frac,
        bracken_reads=bracken_reads,
        alias_used=alias_used,
        reused_bam=reused_bam,
        reused_vcf=reused_vcf,
        runtime_sec=runtime_sec,
        tool_version=__version__,
    )


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
    bbtools_mem: Optional[str] = None,
    dp_min: int = 10,
    maf_min: float = 0.1,
    dp_cap: int = 100,
    mapq: int = 30,
    depth_tool: str = "mosdepth",
    require_pass: bool = False,
    dry_run: bool = False,
    to_stdout: bool = False,
    emit_vcf: bool = False,
    emit_bed: bool = False,
    emit_matrix: bool = False,
    emit_per_contig: bool = False,
    emit_multiqc: bool = False,
    emit_provenance: bool = False,
    **kwargs: Any,
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
        mapq: Minimum mapping quality for depth analysis
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
        dp_cap=dp_cap,
        mapq_min=mapq,
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

    if mode == Mode.REF and all(
        x is None for x in [reference, kwargs.get("species"), kwargs.get("bracken")]
    ):
        raise ValueError(
            "Ref mode requires one of: --reference, --species, or --bracken"
        )

    # Validate depth_tool
    try:
        depth_tool_enum = DepthTool(depth_tool)
    except ValueError:
        valid_depth_tools = [d.value for d in DepthTool]
        raise ValueError(
            f"Invalid depth_tool '{depth_tool}'. Valid options: {valid_depth_tools}"
        )

    # Create run plan
    plan = RunPlan(
        mode=mode,
        sample=sample,
        paths=paths,
        thresholds=thresholds,
        mapper=mapper_enum,
        caller=caller_enum,
        depth_tool=depth_tool_enum,
        bbtools_mem=bbtools_mem,
        require_pass=require_pass,
        threads=threads,
        dry_run=dry_run,
        to_stdout=to_stdout,
        # Extract output flags
        emit_vcf=emit_vcf,
        emit_bed=emit_bed,
        emit_matrix=emit_matrix,
        emit_per_contig=emit_per_contig,
        emit_provenance=emit_provenance,
        emit_multiqc=emit_multiqc,
        tsv_mode=kwargs.get("tsv_mode", TSVMode.OVERWRITE),
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

    Raises:
        ExternalToolError: If required tools are not available
        MappingError: If mapping fails
    """
    logger.info(f"Running self mode for sample {plan.sample}")
    start_time = time.time()

    # Initialize variables for error handling
    ambiguous_snv_count = 0
    all_variants_count = 0
    qc_warnings = ""

    # Ensure thresholds are properly initialized
    assert plan.thresholds.dp_min is not None, "dp_min threshold not set"
    assert plan.thresholds.maf_min is not None, "maf_min threshold not set"
    assert plan.paths.assembly is not None, "Assembly path required for self mode"

    if plan.dry_run:
        print("🔍 DRY RUN - Self-mapping pipeline plan:")
        print(f"  📋 Sample: {plan.sample}")
        print(
            f"  📁 Inputs: {plan.paths.r1.name} + {plan.paths.r2.name} → {plan.paths.assembly.name}"
        )
        print(
            f"  🔧 Tools: {plan.mapper.value} mapper, {plan.caller.value} caller, mosdepth"
        )
        print(
            f"  📊 Thresholds: dp_min={plan.thresholds.dp_min}, maf_min={plan.thresholds.maf_min}, mapq≥{plan.thresholds.mapq_min}"
        )
        print("  ⚙️  Steps planned:")
        print(
            f"    1. Check tool availability ({plan.mapper.value}, samtools, mosdepth, {plan.caller.value})"
        )
        print(f"    2. Ensure/build {plan.mapper.value} indexes for assembly")
        print(f"    3. Map reads to assembly → {plan.sample}.sorted.bam")
        print(f"    4. Run depth analysis → {plan.sample}.depth.mosdepth.summary.txt")
        print(f"    5. Call variants → {plan.sample}.vcf")
        print(f"    6. Normalize/split VCF → {plan.sample}.normalized.vcf.gz")
        print("    7. Count ambiguous sites and generate summary")
        if plan.emit_vcf:
            print(f"    8. Emit filtered VCF → {plan.sample}.ambiguous_sites.vcf.gz")
        if plan.emit_bed:
            print(f"    8. Emit BED → {plan.sample}.ambiguous_sites.bed.gz")
        if plan.emit_matrix:
            print(f"    8. Emit matrix → {plan.sample}.variant_matrix.tsv.gz")
        print(f"  📤 Output: {plan.paths.output_dir}/ambiguous_summary.tsv")
        print("✅ Dry run completed - no files would be written")

        # Detect reuse and extract ref accession
        reused_bam, reused_vcf = detect_reuse_from_plan(plan)
        ref_accession = extract_ref_accession(plan.paths.assembly)

        # Return placeholder data for dry run
        return create_summary_row(
            sample=plan.sample,
            mode=plan.mode.value,
            mapper=plan.mapper.value,
            caller=plan.caller.value,
            dp_min=plan.thresholds.dp_min,
            maf_min=plan.thresholds.maf_min,
            dp_cap=(
                plan.thresholds.dp_cap if plan.thresholds.dp_cap is not None else 100
            ),
            callable_bases=2800000,
            genome_length=2900000,
            ambiguous_snv_count=42,
            ambiguous_indel_count=0,
            ambiguous_del_count=0,
            ref_label=plan.paths.assembly.name if plan.paths.assembly else "NA",
            ref_accession=ref_accession,
            reused_bam=reused_bam,
            reused_vcf=reused_vcf,
            runtime_sec=0.0,
        )

    # Check tool availability
    tools = check_external_tools()
    tool_name = plan.mapper.value  # minimap2 or bwa-mem2
    mosdepth_available = check_mosdepth_available()
    caller_available = caller_tools_available(plan.caller)

    if not tools.get(tool_name, {}).get("available", False) or not tools.get(
        "samtools", {}
    ).get("available", False):
        missing = []
        for t in [tool_name, "samtools"]:
            if not tools.get(t, {}).get("available", False):
                version_info = tools.get(t, {}).get("version", "unknown")
                missing.append(f"{t} ({version_info})")
        raise ExternalToolError(f"Required tools not available: {missing}")

    if not mosdepth_available:
        raise ExternalToolError(
            "mosdepth not found in PATH - required for depth analysis"
        )

    if not caller_available:
        raise ExternalToolError(
            f"Variant caller {plan.caller.value} tools not found in PATH"
        )

    try:
        # Step 1: Ensure indexes exist for self-mode
        logger.info(f"Ensuring indexes for {plan.paths.assembly}")
        ensure_indexes_self(plan.paths.assembly, plan.mapper)

        # Step 2: Map reads to assembly
        logger.info(f"Mapping reads to assembly using {plan.mapper.value}")
        bam_path = map_fastqs(
            mapper=plan.mapper,
            fasta_path=plan.paths.assembly,
            r1_path=plan.paths.r1,
            r2_path=plan.paths.r2,
            sample_name=plan.sample,
            threads=plan.threads,
            output_path=plan.paths.output_dir / f"{plan.sample}.sorted.bam",
        )

        logger.info(f"Mapping completed: {bam_path}")

        # Calculate mapping rate from BAM statistics
        try:
            mapping_rate = calculate_mapping_rate(bam_path)
            logger.info(f"Mapping rate: {mapping_rate:.4f}")
        except (MappingError, ExternalToolError) as e:
            logger.warning(f"Could not calculate mapping rate: {e}")
            mapping_rate = 0.0  # Use fallback value

        # Step 3: Check for duplicates and run markdup if needed
        logger.info("Checking for duplicate flags in BAM")
        final_bam_path = bam_path
        try:
            if not has_duplicate_flags(bam_path):
                logger.info(
                    "No duplicate flags detected, running samtools markdup for depth analysis"
                )
                final_bam_path = run_markdup_for_depth(bam_path, plan.paths.output_dir)
            else:
                logger.info("Duplicate flags already present")
        except CompatibilityError as e:
            logger.warning(f"Could not check/mark duplicates: {e}")
            # Continue with original BAM file
            final_bam_path = bam_path

        # Step 4: Run depth analysis using mosdepth
        logger.info("Running depth analysis with mosdepth")
        try:
            # Ensure thresholds are properly initialized
            assert plan.thresholds.mapq_min is not None, "mapq_min threshold not set"

            depth_summary = analyze_depth(
                bam_path=final_bam_path,  # Use markdup BAM if created
                output_dir=plan.paths.output_dir,
                sample_name=plan.sample,
                mapq_threshold=plan.thresholds.mapq_min,  # Use CLI/config MAPQ threshold
                depth_threshold=plan.thresholds.dp_min,  # Use CLI dp_min
                threads=plan.threads,
            )
            logger.info(
                f"Depth analysis completed: {depth_summary.genome_length:,} bp genome, "
                f"{depth_summary.breadth_10x:.2%} breadth ≥{plan.thresholds.dp_min}x"
            )
        except DepthAnalysisError as e:
            logger.error(f"Depth analysis failed: {e}")
            # Use placeholder values if depth analysis fails
            depth_summary = None
            qc_warnings = f"depth_analysis_failed:{e}"

        # Step 5: Run variant calling
        logger.info(f"Running variant calling with {plan.caller.value}")
        try:
            # Ensure thresholds are properly initialized
            assert plan.thresholds.mapq_min is not None, "mapq_min threshold not set"
            assert plan.thresholds.baseq_min is not None, "baseq_min threshold not set"

            vcf_path = plan.paths.output_dir / f"{plan.sample}.vcf"
            variant_result = call_variants(
                bam_path=bam_path,  # Use original BAM for variant calling (not markdup)
                reference_path=plan.paths.assembly,
                output_vcf=vcf_path,
                caller=plan.caller,
                sample_name=plan.sample,
                threads=plan.threads,
                mapq_min=plan.thresholds.mapq_min,
                baseq_min=plan.thresholds.baseq_min,
                minallelefraction=0.0,  # Get all variants for analysis
                bbtools_mem=plan.bbtools_mem,
            )
            logger.info(f"Variant calling completed: {variant_result.vcf_path}")
        except VariantCallingError as e:
            logger.error(f"Variant calling failed: {e}")
            raise ExternalToolError(f"Variant calling failed: {e}")

        # Step 5: Normalize VCF and count ambiguous sites
        logger.info("Normalizing VCF and counting ambiguous sites")
        normalized_vcf_path = None  # Initialize to avoid NameError in exception handler
        try:
            # Normalize and split VCF
            normalized_vcf_path = normalize_and_split(
                vcf_in=variant_result.vcf_path,
                reference=plan.paths.assembly,
                output_dir=plan.paths.output_dir,
            )

            # Get included contigs from depth analysis
            included_contigs = None
            if depth_summary:
                # Use the same contig filtering as depth analysis (≥500 bp)
                from .depth import list_included_contigs

                mosdepth_summary = (
                    plan.paths.output_dir / f"{plan.sample}.mosdepth.summary.txt"
                )
                if mosdepth_summary.exists():
                    included_contigs = list_included_contigs(
                        mosdepth_summary, min_len=500
                    )

            # Count ambiguous sites (SNVs only for primary count)
            ambiguous_snv_count, grid = count_ambiguous_sites(
                vcf_path=normalized_vcf_path,
                dp_min=plan.thresholds.dp_min,
                maf_min=plan.thresholds.maf_min,
                dp_cap=100,  # As per spec
                included_contigs=included_contigs,
                variant_classes=[VariantClass.SNV],
            )

            # Count total variant sites by class
            all_variants_count, _ = count_ambiguous_sites(
                vcf_path=normalized_vcf_path,
                dp_min=0,  # Count all variants regardless of thresholds
                maf_min=0.0,
                dp_cap=100,
                included_contigs=included_contigs,
                variant_classes=[VariantClass.SNV, VariantClass.INS, VariantClass.DEL],
            )

            logger.info(
                f"Found {ambiguous_snv_count} ambiguous SNVs, {all_variants_count} total variant sites"
            )

            # Optional emitters (only create files if flags are set)
            if plan.emit_vcf:
                vcf_output_path = (
                    plan.paths.output_dir / f"{plan.sample}.ambiguous_sites.vcf.gz"
                )
                logger.info(f"Emitting VCF to {vcf_output_path}")
                emit_vcf(
                    normalized_vcf_path=normalized_vcf_path,
                    output_path=vcf_output_path,
                    dp_min=plan.thresholds.dp_min,
                    maf_min=plan.thresholds.maf_min,
                    sample_name=plan.sample,
                    require_pass=False,  # TODO: Add this to plan if needed
                    included_contigs=included_contigs,
                )

            if plan.emit_bed:
                bed_output_path = (
                    plan.paths.output_dir / f"{plan.sample}.ambiguous_sites.bed.gz"
                )
                logger.info(f"Emitting BED to {bed_output_path}")
                emit_bed(
                    normalized_vcf_path=normalized_vcf_path,
                    output_path=bed_output_path,
                    dp_min=plan.thresholds.dp_min,
                    maf_min=plan.thresholds.maf_min,
                    sample_name=plan.sample,
                    included_contigs=included_contigs,
                )

            if plan.emit_matrix:
                matrix_output_path = (
                    plan.paths.output_dir / f"{plan.sample}.ambiguity_matrix.tsv.gz"
                )
                logger.info(f"Emitting ambiguity matrix to {matrix_output_path}")
                emit_matrix(
                    grid=grid, output_path=matrix_output_path, sample_name=plan.sample
                )

            if plan.emit_per_contig:
                per_contig_output_path = (
                    plan.paths.output_dir / f"{plan.sample}.per_contig.tsv"
                )
                logger.info(f"Emitting per-contig summary to {per_contig_output_path}")
                mosdepth_summary = (
                    plan.paths.output_dir / f"{plan.sample}.mosdepth.summary.txt"
                )
                emit_per_contig(
                    depth_summary_path=mosdepth_summary,
                    normalized_vcf_path=normalized_vcf_path,
                    output_path=per_contig_output_path,
                    dp_min=plan.thresholds.dp_min,
                    maf_min=plan.thresholds.maf_min,
                    sample_name=plan.sample,
                    included_contigs=included_contigs,
                )

        except VCFOperationError as e:
            logger.error(f"VCF analysis failed: {e}")
            # Use placeholder values if VCF analysis fails
            ambiguous_snv_count = 0
            all_variants_count = 0
            if qc_warnings:
                qc_warnings += ";vcf_analysis_failed"
            else:
                qc_warnings = "vcf_analysis_failed"

        # Step 6: Extract metrics for SummaryRow
        if depth_summary:
            callable_bases = depth_summary.callable_bases
            genome_length = depth_summary.genome_length
            breadth_10x = depth_summary.breadth_10x
            qc_warnings = ""

            # Add QC warning if no contigs meet length threshold
            if depth_summary.included_contigs == 0:
                qc_warnings = "no_contigs_ge_500bp"
            elif depth_summary.included_contigs < depth_summary.total_contigs:
                excluded = depth_summary.total_contigs - depth_summary.included_contigs
                qc_warnings = f"excluded_{excluded}_short_contigs"
        else:
            # Fallback values if depth analysis failed
            callable_bases = 0
            genome_length = 0
            breadth_10x = 0.0
            qc_warnings = qc_warnings or "depth_analysis_unavailable"

        # TODO: Add mapping rate calculation from BAM stats
        # Calculate ref_calls (sites that are reference-like, i.e., not variant)
        (
            callable_bases - all_variants_count
            if callable_bases > all_variants_count
            else 0
        )

        # Emit MultiQC metrics if requested
        if plan.emit_multiqc:
            multiqc_output_path = plan.paths.output_dir / f"{plan.sample}.multiqc.tsv"
            logger.info(f"Emitting MultiQC metrics to {multiqc_output_path}")
            emit_multiqc(
                sample_name=plan.sample,
                ambiguous_snv_count=ambiguous_snv_count,
                breadth_10x=breadth_10x,
                callable_bases=callable_bases,
                genome_length=genome_length,
                dp_min=plan.thresholds.dp_min,
                maf_min=plan.thresholds.maf_min,
                mapper=plan.mapper.value,
                caller=plan.caller.value,
                mode=plan.mode.value,
                output_path=multiqc_output_path,
            )

        # Calculate variant counts by type (if VCF normalization succeeded)
        ambiguous_indel_count = 0
        ambiguous_del_count = 0
        ref_grid: Optional[AmbigGrid] = None

        if normalized_vcf_path is not None:
            _, ref_grid = count_ambiguous_sites(
                vcf_path=normalized_vcf_path,
                dp_min=plan.thresholds.dp_min,
                maf_min=plan.thresholds.maf_min,
                dp_cap=(
                    plan.thresholds.dp_cap
                    if plan.thresholds.dp_cap is not None
                    else 100
                ),
                included_contigs=included_contigs,
                variant_classes=[VariantClass.SNV],
            )

            # Count indels and deletions separately
            ambiguous_indel_count, _ = count_ambiguous_sites(
                vcf_path=normalized_vcf_path,
                dp_min=plan.thresholds.dp_min,
                maf_min=plan.thresholds.maf_min,
                dp_cap=(
                    plan.thresholds.dp_cap
                    if plan.thresholds.dp_cap is not None
                    else 100
                ),
                included_contigs=included_contigs,
                variant_classes=[VariantClass.INS],
            )

            ambiguous_del_count, _ = count_ambiguous_sites(
                vcf_path=normalized_vcf_path,
                dp_min=plan.thresholds.dp_min,
                maf_min=plan.thresholds.maf_min,
                dp_cap=(
                    plan.thresholds.dp_cap
                    if plan.thresholds.dp_cap is not None
                    else 100
                ),
                included_contigs=included_contigs,
                variant_classes=[VariantClass.DEL],
            )
        else:
            logger.warning(
                "Skipping variant type counting due to VCF normalization failure"
            )

        # Detect reuse and extract ref accession
        reused_bam, reused_vcf = detect_reuse_from_plan(plan)
        ref_accession = extract_ref_accession(plan.paths.assembly)

        return create_summary_row(
            sample=plan.sample,
            mode=plan.mode.value,
            mapper=plan.mapper.value,
            caller=plan.caller.value,
            dp_min=plan.thresholds.dp_min,
            maf_min=plan.thresholds.maf_min,
            dp_cap=(
                plan.thresholds.dp_cap if plan.thresholds.dp_cap is not None else 100
            ),
            callable_bases=callable_bases,
            genome_length=genome_length,
            ambiguous_snv_count=ambiguous_snv_count,
            ambiguous_indel_count=ambiguous_indel_count,
            ambiguous_del_count=ambiguous_del_count,
            ref_label=plan.paths.assembly.name if plan.paths.assembly else "NA",
            ref_accession=ref_accession,
            reused_bam=reused_bam,
            reused_vcf=reused_vcf,
            runtime_sec=time.time() - start_time,
            mapping_rate=mapping_rate,
        )

    except (ExternalToolError, MappingError) as e:
        logger.error(f"Self-mapping failed: {e}")
        raise
    except DepthAnalysisError as e:
        logger.error(f"Depth analysis failed: {e}")
        # Detect reuse and extract ref accession
        reused_bam, reused_vcf = detect_reuse_from_plan(plan)
        ref_accession = extract_ref_accession(plan.paths.assembly)

        # Return summary with depth analysis failure
        return create_summary_row(
            sample=plan.sample,
            mode=plan.mode.value,
            mapper=plan.mapper.value,
            caller=plan.caller.value,
            dp_min=plan.thresholds.dp_min,
            maf_min=plan.thresholds.maf_min,
            dp_cap=(
                plan.thresholds.dp_cap if plan.thresholds.dp_cap is not None else 100
            ),
            callable_bases=0,
            genome_length=0,
            ambiguous_snv_count=0,
            ambiguous_indel_count=0,
            ambiguous_del_count=0,
            ref_label=plan.paths.assembly.name if plan.paths.assembly else "unknown",
            ref_accession=ref_accession,
            reused_bam=reused_bam,
            reused_vcf=reused_vcf,
            runtime_sec=0.0,
        )


def run_ref(
    plan: RunPlan,
    species: Optional[str] = None,
    bracken: Optional[Path] = None,
    ref_dir: Optional[Path] = None,
    on_fail: str = "error",
    min_bracken_frac: float = 0.70,
    min_bracken_reads: int = 100000,
) -> SummaryRow:
    """
    Execute reference-mapping mode.

    Args:
        plan: Validated execution plan
        **kwargs: Additional arguments (species, bracken, ref_dir, on_fail)

    Returns:
        Summary row with results
    """
    logger.info(f"Running ref mode for sample {plan.sample}")
    start_time = time.time()

    # Initialize variables for error handling
    ambiguous_snv_count = 0
    all_variants_count = 0
    qc_warnings = ""

    # Ensure thresholds are properly initialized
    assert plan.thresholds.dp_min is not None, "dp_min threshold not set"
    assert plan.thresholds.maf_min is not None, "maf_min threshold not set"

    # Resolve reference using precedence: file > species > bracken
    ref_source = "unknown"
    ref_label = "unknown"
    reference_path = None
    bracken_species = "NA"
    bracken_frac = 0.0
    bracken_reads = 0

    try:
        if plan.paths.reference:
            # Direct reference file provided
            ref_source = "file"
            ref_label = plan.paths.reference.name
            reference_path = plan.paths.reference

        elif species:
            # Species lookup in admin directory
            ref_source = "species"
            ref_dir = resolve_reference_directory(ref_dir)

            # Resolve species to FASTA file
            reference_path = resolve_species_to_fasta(ref_dir, species)
            ref_label = reference_path.name

            # Try to get accession for better labeling
            accession = parse_accession_from_fasta_header(reference_path)
            if accession:
                ref_label = f"{species} ({accession})"
            else:
                ref_label = species

        elif bracken:
            # Bracken-based species selection
            from .bracken import (
                select_species_from_file,
                BrackenThresholds,
                BrackenSelectionError,
            )

            ref_source = "bracken"
            bracken_path = Path(bracken)

            # Get thresholds from configuration or defaults
            min_frac = min_bracken_frac
            min_reads = min_bracken_reads
            thresholds = BrackenThresholds(min_frac=min_frac, min_reads=min_reads)

            try:
                # Select species from Bracken
                selection = select_species_from_file(bracken_path, thresholds)

                if selection is None:
                    raise ReferenceResolutionError(
                        f"No species in Bracken file {bracken_path} meet selection criteria. "
                        f"Try using --species or --reference instead."
                    )

                # Store bracken information
                bracken_species = selection.species.name
                bracken_frac = selection.species.fraction_total_reads
                bracken_reads = selection.species.new_est_reads

                # Extract species name from Bracken (e.g., "Listeria monocytogenes" from potentially complex name)
                bracken_name = selection.species.name
                parts = bracken_name.strip().split()
                if len(parts) >= 2:
                    species_for_resolution = (
                        f"{parts[0]} {parts[1]}"  # "Listeria monocytogenes"
                    )
                else:
                    species_for_resolution = bracken_name

                # Resolve species to reference file
                ref_dir = resolve_reference_directory(ref_dir)
                reference_path = resolve_species_to_fasta(
                    ref_dir, species_for_resolution
                )

                # Create detailed label with Bracken info
                from .refdir import normalize_species_name

                accession = parse_accession_from_fasta_header(reference_path)
                normalized_name = normalize_species_name(
                    species_for_resolution
                )  # For label only
                if accession:
                    ref_label = f"{normalized_name} ({accession}) [Bracken: {selection.species.name}]"
                else:
                    ref_label = f"{normalized_name} [Bracken: {selection.species.name}]"

            except BrackenSelectionError as e:
                # Convert BrackenSelectionError to ReferenceResolutionError to use existing --on-fail logic
                raise ReferenceResolutionError(str(e))

        else:
            raise ReferenceResolutionError(
                "No reference specified. Provide one of: --reference, --species, or --bracken"
            )

    except ReferenceResolutionError as e:
        if on_fail == "error":
            raise
        elif on_fail == "warn":
            logger.warning(f"Reference resolution failed: {e}")
            # Continue with placeholder data
            ref_source = "failed"
            ref_label = "resolution_failed"
        elif on_fail == "skip":
            logger.info(f"Skipping due to reference resolution failure: {e}")
            # Detect reuse (though unlikely in skip case)
            reused_bam, reused_vcf = detect_reuse_from_plan(plan)

            # Return a minimal result indicating skip
            return create_summary_row(
                sample=plan.sample,
                mode=plan.mode.value,
                mapper=plan.mapper.value,
                caller=plan.caller.value,
                dp_min=plan.thresholds.dp_min,
                maf_min=plan.thresholds.maf_min,
                dp_cap=(
                    plan.thresholds.dp_cap
                    if plan.thresholds.dp_cap is not None
                    else 100
                ),
                callable_bases=0,
                genome_length=0,
                ambiguous_snv_count=0,
                ambiguous_indel_count=0,
                ambiguous_del_count=0,
                ref_label="SKIPPED",
                ref_accession="NA",
                bracken_species="NA",  # No bracken data for skipped samples
                bracken_frac=0.0,
                bracken_reads=0,
                reused_bam=reused_bam,
                reused_vcf=reused_vcf,
                runtime_sec=0.0,
            )

    if plan.dry_run:
        print("🔍 DRY RUN - Reference-mapping pipeline plan:")
        print(f"  📋 Sample: {plan.sample}")
        print(f"  📁 Inputs: {plan.paths.r1.name} + {plan.paths.r2.name}")
        print(f"  🧬 Reference source: {ref_source}")
        print(f"  📖 Reference: {reference_path.name if reference_path else 'TBD'}")
        print(
            f"  🔧 Tools: {plan.mapper.value} mapper, {plan.caller.value} caller, mosdepth"
        )
        print(
            f"  📊 Thresholds: dp_min={plan.thresholds.dp_min}, maf_min={plan.thresholds.maf_min}, mapq≥{plan.thresholds.mapq_min}"
        )
        print("  ⚙️  Steps planned:")
        print(f"    1. Resolve reference from {ref_source}")
        print(
            f"    2. Check tool availability ({plan.mapper.value}, samtools, mosdepth, {plan.caller.value})"
        )
        print(f"    3. Map reads to reference → {plan.sample}.sorted.bam")
        print(f"    4. Run depth analysis → {plan.sample}.depth.mosdepth.summary.txt")
        print(f"    5. Call variants → {plan.sample}.vcf")
        print(f"    6. Normalize/split VCF → {plan.sample}.normalized.vcf.gz")
        print("    7. Count ambiguous sites and generate summary")
        if plan.emit_vcf:
            print(f"    8. Emit filtered VCF → {plan.sample}.ambiguous_sites.vcf.gz")
        if plan.emit_bed:
            print(f"    8. Emit BED → {plan.sample}.ambiguous_sites.bed.gz")
        if plan.emit_matrix:
            print(f"    8. Emit matrix → {plan.sample}.variant_matrix.tsv.gz")
        print(f"  📤 Output: {plan.paths.output_dir}/ambiguous_summary.tsv")
        print("✅ Dry run completed - no files would be written")

        # Detect reuse and extract ref accession
        reused_bam, reused_vcf = detect_reuse_from_plan(plan)
        ref_accession = "NA"
        if reference_path:
            ref_accession = extract_ref_accession(reference_path)

        # Return placeholder data for dry run
        return create_summary_row(
            sample=plan.sample,
            mode=plan.mode.value,
            mapper=plan.mapper.value,
            caller=plan.caller.value,
            dp_min=plan.thresholds.dp_min,
            maf_min=plan.thresholds.maf_min,
            dp_cap=(
                plan.thresholds.dp_cap if plan.thresholds.dp_cap is not None else 100
            ),
            callable_bases=2750000,  # Reference genome size
            genome_length=2800000,
            ambiguous_snv_count=67,
            ambiguous_indel_count=0,
            ambiguous_del_count=0,
            ref_label=ref_label,
            ref_accession=ref_accession,
            bracken_species=bracken_species,
            bracken_frac=bracken_frac,
            bracken_reads=bracken_reads,
            reused_bam=reused_bam,
            reused_vcf=reused_vcf,
            runtime_sec=0.0,
        )

    # Update the plan's reference info
    plan.ref_source = ref_source
    plan.ref_label = ref_label
    if reference_path:
        plan.paths.reference = reference_path

    # Check tool availability
    tools = check_external_tools()
    tool_name = plan.mapper.value  # minimap2 or bwa-mem2
    mosdepth_available = check_mosdepth_available()
    caller_available = caller_tools_available(plan.caller)

    if not tools.get(tool_name, {}).get("available", False) or not tools.get(
        "samtools", {}
    ).get("available", False):
        missing = []
        for t in [tool_name, "samtools"]:
            if not tools.get(t, {}).get("available", False):
                version_info = tools.get(t, {}).get("version", "unknown")
                missing.append(f"{t} ({version_info})")
        raise ExternalToolError(f"Required tools not available: {missing}")

    if not mosdepth_available:
        raise ExternalToolError(
            "mosdepth not found in PATH - required for depth analysis"
        )

    if not caller_available:
        raise ExternalToolError(
            f"Variant caller {plan.caller.value} tools not found in PATH"
        )

    # Only proceed if we have a valid reference
    if ref_source == "failed" or not reference_path:
        # Detect reuse
        reused_bam, reused_vcf = detect_reuse_from_plan(plan)

        return create_summary_row(
            sample=plan.sample,
            mode=plan.mode.value,
            mapper=plan.mapper.value,
            caller=plan.caller.value,
            dp_min=plan.thresholds.dp_min,
            maf_min=plan.thresholds.maf_min,
            dp_cap=(
                plan.thresholds.dp_cap if plan.thresholds.dp_cap is not None else 100
            ),
            callable_bases=0,
            genome_length=0,
            ambiguous_snv_count=0,
            ambiguous_indel_count=0,
            ambiguous_del_count=0,
            ref_label=ref_label,
            ref_accession="NA",
            reused_bam=reused_bam,
            reused_vcf=reused_vcf,
            runtime_sec=0.0,
        )

    try:
        # Step 1: Map reads to reference
        logger.info(f"Mapping reads to reference using {plan.mapper.value}")
        bam_path = map_fastqs(
            mapper=plan.mapper,
            fasta_path=reference_path,
            r1_path=plan.paths.r1,
            r2_path=plan.paths.r2,
            sample_name=plan.sample,
            threads=plan.threads,
            output_path=plan.paths.output_dir / f"{plan.sample}.sorted.bam",
        )

        logger.info(f"Mapping completed: {bam_path}")

        # Calculate mapping rate from BAM statistics
        try:
            mapping_rate = calculate_mapping_rate(bam_path)
            logger.info(f"Mapping rate: {mapping_rate:.4f}")
        except (MappingError, ExternalToolError) as e:
            logger.warning(f"Could not calculate mapping rate: {e}")
            mapping_rate = 0.0  # Use fallback value

        # Step 2: Check for duplicates and run markdup if needed
        logger.info("Checking for duplicate flags in BAM")
        final_bam_path = bam_path
        try:
            if not has_duplicate_flags(bam_path):
                logger.info(
                    "No duplicate flags detected, running samtools markdup for depth analysis"
                )
                final_bam_path = run_markdup_for_depth(bam_path, plan.paths.output_dir)
            else:
                logger.info("Duplicate flags already present")
        except CompatibilityError as e:
            logger.warning(f"Could not check/mark duplicates: {e}")
            # Continue with original BAM file
            final_bam_path = bam_path

        # Step 3: Run depth analysis using mosdepth
        logger.info("Running depth analysis with mosdepth")
        try:
            # Ensure thresholds are properly initialized
            assert plan.thresholds.mapq_min is not None, "mapq_min threshold not set"

            depth_summary = analyze_depth(
                bam_path=final_bam_path,  # Use markdup BAM if created
                output_dir=plan.paths.output_dir,
                sample_name=plan.sample,
                mapq_threshold=plan.thresholds.mapq_min,  # Use CLI/config MAPQ threshold
                depth_threshold=plan.thresholds.dp_min,  # Use CLI dp_min
                threads=plan.threads,
            )
            logger.info(
                f"Depth analysis completed: {depth_summary.genome_length:,} bp genome, "
                f"{depth_summary.breadth_10x:.2%} breadth ≥{plan.thresholds.dp_min}x"
            )
        except DepthAnalysisError as e:
            logger.error(f"Depth analysis failed: {e}")
            # Use placeholder values if depth analysis fails
            depth_summary = None
            qc_warnings = f"depth_analysis_failed:{e}"

        # Step 4: Run variant calling
        logger.info(f"Running variant calling with {plan.caller.value}")
        try:
            # Ensure thresholds are properly initialized
            assert plan.thresholds.mapq_min is not None, "mapq_min threshold not set"
            assert plan.thresholds.baseq_min is not None, "baseq_min threshold not set"

            vcf_path = plan.paths.output_dir / f"{plan.sample}.vcf"
            variant_result = call_variants(
                bam_path=bam_path,
                reference_path=reference_path,
                output_vcf=vcf_path,
                caller=plan.caller,
                sample_name=plan.sample,
                threads=plan.threads,
                mapq_min=plan.thresholds.mapq_min,
                baseq_min=plan.thresholds.baseq_min,
                minallelefraction=0.0,  # Get all variants for analysis
                bbtools_mem=plan.bbtools_mem,
            )
            logger.info(f"Variant calling completed: {variant_result.vcf_path}")
        except VariantCallingError as e:
            logger.error(f"Variant calling failed: {e}")
            raise ExternalToolError(f"Variant calling failed: {e}")

        # Step 4: Normalize VCF and count ambiguous sites
        logger.info("Normalizing VCF and counting ambiguous sites")
        try:
            # Normalize and split VCF
            normalized_vcf_path = normalize_and_split(
                vcf_in=variant_result.vcf_path,
                reference=reference_path,
                output_dir=plan.paths.output_dir,
            )

            # Get included contigs from depth analysis
            included_contigs = None
            if depth_summary:
                # Use the same contig filtering as depth analysis (≥500 bp)
                from .depth import list_included_contigs

                mosdepth_summary = (
                    plan.paths.output_dir / f"{plan.sample}.mosdepth.summary.txt"
                )
                if mosdepth_summary.exists():
                    included_contigs = list_included_contigs(
                        mosdepth_summary, min_len=500
                    )

            # Count ambiguous sites (SNVs only for primary count)
            ambiguous_snv_count, grid = count_ambiguous_sites(
                vcf_path=normalized_vcf_path,
                dp_min=plan.thresholds.dp_min,
                maf_min=plan.thresholds.maf_min,
                dp_cap=100,  # As per spec
                included_contigs=included_contigs,
                variant_classes=[VariantClass.SNV],
            )

            # Count total variant sites by class
            all_variants_count, _ = count_ambiguous_sites(
                vcf_path=normalized_vcf_path,
                dp_min=0,  # Count all variants regardless of thresholds
                maf_min=0.0,
                dp_cap=100,
                included_contigs=included_contigs,
                variant_classes=[VariantClass.SNV, VariantClass.INS, VariantClass.DEL],
            )

            logger.info(
                f"Found {ambiguous_snv_count} ambiguous SNVs, {all_variants_count} total variant sites"
            )

            # Optional emitters (only create files if flags are set)
            if plan.emit_vcf:
                vcf_output_path = (
                    plan.paths.output_dir / f"{plan.sample}.ambiguous_sites.vcf.gz"
                )
                logger.info(f"Emitting VCF to {vcf_output_path}")
                emit_vcf(
                    normalized_vcf_path=normalized_vcf_path,
                    output_path=vcf_output_path,
                    dp_min=plan.thresholds.dp_min,
                    maf_min=plan.thresholds.maf_min,
                    sample_name=plan.sample,
                    require_pass=False,  # TODO: Add this to plan if needed
                    included_contigs=included_contigs,
                )

            if plan.emit_bed:
                bed_output_path = (
                    plan.paths.output_dir / f"{plan.sample}.ambiguous_sites.bed.gz"
                )
                logger.info(f"Emitting BED to {bed_output_path}")
                emit_bed(
                    normalized_vcf_path=normalized_vcf_path,
                    output_path=bed_output_path,
                    dp_min=plan.thresholds.dp_min,
                    maf_min=plan.thresholds.maf_min,
                    sample_name=plan.sample,
                    included_contigs=included_contigs,
                )

            if plan.emit_matrix:
                matrix_output_path = (
                    plan.paths.output_dir / f"{plan.sample}.ambiguity_matrix.tsv.gz"
                )
                logger.info(f"Emitting ambiguity matrix to {matrix_output_path}")
                emit_matrix(
                    grid=grid, output_path=matrix_output_path, sample_name=plan.sample
                )

            if plan.emit_per_contig:
                per_contig_output_path = (
                    plan.paths.output_dir / f"{plan.sample}.per_contig.tsv"
                )
                logger.info(f"Emitting per-contig summary to {per_contig_output_path}")
                mosdepth_summary = (
                    plan.paths.output_dir / f"{plan.sample}.mosdepth.summary.txt"
                )
                emit_per_contig(
                    depth_summary_path=mosdepth_summary,
                    normalized_vcf_path=normalized_vcf_path,
                    output_path=per_contig_output_path,
                    dp_min=plan.thresholds.dp_min,
                    maf_min=plan.thresholds.maf_min,
                    sample_name=plan.sample,
                    included_contigs=included_contigs,
                )

        except VCFOperationError as e:
            logger.error(f"VCF analysis failed: {e}")
            # Use placeholder values if VCF analysis fails
            ambiguous_snv_count = 0
            all_variants_count = 0
            if qc_warnings:
                qc_warnings += ";vcf_analysis_failed"
            else:
                qc_warnings = "vcf_analysis_failed"

        # Step 5: Extract metrics for SummaryRow
        if depth_summary:
            callable_bases = depth_summary.callable_bases
            genome_length = depth_summary.genome_length
            breadth_10x = depth_summary.breadth_10x
            qc_warnings = ""

            # Add QC warning if no contigs meet length threshold
            if depth_summary.included_contigs == 0:
                qc_warnings = "no_contigs_ge_500bp"
            elif depth_summary.included_contigs < depth_summary.total_contigs:
                excluded = depth_summary.total_contigs - depth_summary.included_contigs
                qc_warnings = f"excluded_{excluded}_short_contigs"
        else:
            # Fallback values if depth analysis failed
            callable_bases = 0
            genome_length = 0
            breadth_10x = 0.0
            qc_warnings = qc_warnings or "depth_analysis_unavailable"

        # TODO: Add mapping rate calculation from BAM stats
        # Calculate ref_calls (sites that are reference-like, i.e., not variant)
        (
            callable_bases - all_variants_count
            if callable_bases > all_variants_count
            else 0
        )

        # Emit MultiQC metrics if requested
        if plan.emit_multiqc:
            multiqc_output_path = plan.paths.output_dir / f"{plan.sample}.multiqc.tsv"
            logger.info(f"Emitting MultiQC metrics to {multiqc_output_path}")
            emit_multiqc(
                sample_name=plan.sample,
                ambiguous_snv_count=ambiguous_snv_count,
                breadth_10x=breadth_10x,
                callable_bases=callable_bases,
                genome_length=genome_length,
                dp_min=plan.thresholds.dp_min,
                maf_min=plan.thresholds.maf_min,
                mapper=plan.mapper.value,
                caller=plan.caller.value,
                mode=plan.mode.value,
                output_path=multiqc_output_path,
            )

        # Calculate variant counts by type
        _, grid = count_ambiguous_sites(
            vcf_path=normalized_vcf_path,
            dp_min=plan.thresholds.dp_min,
            maf_min=plan.thresholds.maf_min,
            dp_cap=(
                plan.thresholds.dp_cap if plan.thresholds.dp_cap is not None else 100
            ),
            included_contigs=included_contigs,
            variant_classes=[VariantClass.SNV],
        )

        # Count indels and deletions separately
        ambiguous_indel_count, _ = count_ambiguous_sites(
            vcf_path=normalized_vcf_path,
            dp_min=plan.thresholds.dp_min,
            maf_min=plan.thresholds.maf_min,
            dp_cap=(
                plan.thresholds.dp_cap if plan.thresholds.dp_cap is not None else 100
            ),
            included_contigs=included_contigs,
            variant_classes=[VariantClass.INS],
        )

        ambiguous_del_count, _ = count_ambiguous_sites(
            vcf_path=normalized_vcf_path,
            dp_min=plan.thresholds.dp_min,
            maf_min=plan.thresholds.maf_min,
            dp_cap=(
                plan.thresholds.dp_cap if plan.thresholds.dp_cap is not None else 100
            ),
            included_contigs=included_contigs,
            variant_classes=[VariantClass.DEL],
        )

        # Detect reuse and extract ref accession
        reused_bam, reused_vcf = detect_reuse_from_plan(plan)
        ref_accession = (
            extract_ref_accession(reference_path) if reference_path else "NA"
        )

        return create_summary_row(
            sample=plan.sample,
            mode=plan.mode.value,
            mapper=plan.mapper.value,
            caller=plan.caller.value,
            dp_min=plan.thresholds.dp_min,
            maf_min=plan.thresholds.maf_min,
            dp_cap=(
                plan.thresholds.dp_cap if plan.thresholds.dp_cap is not None else 100
            ),
            callable_bases=callable_bases,
            genome_length=genome_length,
            ambiguous_snv_count=ambiguous_snv_count,
            ambiguous_indel_count=ambiguous_indel_count,
            ambiguous_del_count=ambiguous_del_count,
            ref_label=ref_label,
            ref_accession=ref_accession,
            bracken_species=bracken_species,
            bracken_frac=bracken_frac,
            bracken_reads=bracken_reads,
            reused_bam=reused_bam,
            reused_vcf=reused_vcf,
            runtime_sec=time.time() - start_time,
            mapping_rate=mapping_rate,
        )

    except (ExternalToolError, MappingError) as e:
        logger.error(f"Reference-mapping failed: {e}")
        raise
    except DepthAnalysisError as e:
        logger.error(f"Depth analysis failed: {e}")
        # Detect reuse and extract ref accession
        reused_bam, reused_vcf = detect_reuse_from_plan(plan)
        ref_accession = (
            extract_ref_accession(reference_path) if reference_path else "NA"
        )

        # Return summary with depth analysis failure
        return create_summary_row(
            sample=plan.sample,
            mode=plan.mode.value,
            mapper=plan.mapper.value,
            caller=plan.caller.value,
            dp_min=plan.thresholds.dp_min,
            maf_min=plan.thresholds.maf_min,
            dp_cap=(
                plan.thresholds.dp_cap if plan.thresholds.dp_cap is not None else 100
            ),
            callable_bases=0,
            genome_length=0,
            ambiguous_snv_count=0,
            ambiguous_indel_count=0,
            ambiguous_del_count=0,
            ref_label=ref_label,
            ref_accession=ref_accession,
            bracken_species=bracken_species,
            bracken_frac=bracken_frac,
            bracken_reads=bracken_reads,
            reused_bam=reused_bam,
            reused_vcf=reused_vcf,
            runtime_sec=0.0,
            mapping_rate=0.0,  # Use fallback for failed depth analysis
        )


def run_summarize(
    vcf: Path,
    bam: Path,
    output: Optional[Path] = None,
    dp_min: int = 10,
    maf_min: float = 0.1,
    dp_cap: int = 100,
    require_pass: bool = False,
    emit_vcf: bool = False,
    emit_bed: bool = False,
    emit_matrix: bool = False,
    emit_per_contig: bool = False,
    emit_multiqc: bool = False,
    emit_provenance: bool = False,
    to_stdout: bool = False,
) -> List[SummaryRow]:
    """
    Execute summarize mode.

    Requires both VCF and BAM files. Reuses denominator pipeline from BAM
    and counting logic from VCF to produce summary statistics.

    Args:
        vcf: VCF file to summarize (required)
        bam: BAM file for denominator calculation (required)
        output: Output file (auto-generated if not provided)
        dp_min: Minimum depth threshold
        maf_min: Minimum MAF threshold
        dp_cap: Maximum depth cap for variant analysis
        require_pass: Only consider variants that pass caller filters
        emit_vcf: Emit filtered VCF
        emit_bed: Emit BED file
        emit_matrix: Emit variant matrix
        emit_per_contig: Emit per-contig summary
        emit_multiqc: Emit MultiQC metrics
        emit_provenance: Emit JSON provenance file
        to_stdout: Write to stdout instead of file

    Returns:
        List of summary rows

    Raises:
        FileNotFoundError: If VCF or BAM files are missing
        ValueError: If inputs are invalid
    """
    from .depth import analyze_depth, list_included_contigs
    from .vcf_ops import (
        count_ambiguous_sites,
        emit_vcf as emit_vcf_file,
        emit_bed as emit_bed_file,
        emit_matrix as emit_matrix_file,
        VariantClass,
    )
    from .io_utils import (
        infer_sample_name,
        validate_sample_name,
        write_tsv_summary,
        write_tsv_to_stdout,
    )
    from .models import TSVMode

    logger.info(f"Running summarize mode on VCF: {vcf}, BAM: {bam}")
    start_time = time.time()

    try:
        # Step 1: Validate required inputs
        if not vcf.exists():
            raise FileNotFoundError(f"VCF file not found: {vcf}")
        if not bam.exists():
            raise FileNotFoundError(f"BAM file not found: {bam}")

        # Validate that stdout mode doesn't conflict with emit flags
        if to_stdout and any(
            [
                emit_vcf,
                emit_bed,
                emit_matrix,
                emit_per_contig,
                emit_multiqc,
                emit_provenance,
            ]
        ):
            raise ValueError("Cannot use --stdout with --emit-* flags")

        # Step 2: Infer and validate sample name from inputs
        sample_name = infer_sample_name(
            r1=Path(
                "dummy.fastq"
            ),  # Required by function but not used for VCF inference
            vcf=vcf,
            bam=bam,
        )
        validate_sample_name(sample_name)
        logger.info(f"Using sample name: {sample_name}")

        # Step 3: Denominator calculation from BAM using depth analysis
        logger.info("Calculating denominator from BAM using depth analysis")
        output_dir = output.parent if output else Path.cwd()
        temp_depth_dir = output_dir / "temp_depth"
        depth_summary = analyze_depth(
            bam_path=bam,
            output_dir=temp_depth_dir,
            sample_name=sample_name,
            mapq_threshold=30,
            depth_threshold=10,
            threads=4,
        )

        # Get included contigs (>=500 bp) for consistent filtering
        summary_file = temp_depth_dir / f"{sample_name}.depth.mosdepth.summary.txt"
        included_contigs = list_included_contigs(summary_file, min_len=500)

        # Step 4: Count ambiguous SNVs from VCF
        logger.info("Counting ambiguous SNVs from VCF")
        ambiguous_snv_count, grid = count_ambiguous_sites(
            vcf_path=vcf,
            dp_min=dp_min,
            maf_min=maf_min,
            dp_cap=dp_cap,
            included_contigs=included_contigs,
            variant_classes=[VariantClass.SNV],
        )

        # Step 5: Count indels and deletions for secondary metrics
        logger.info("Counting indels and deletions")
        indel_count, _ = count_ambiguous_sites(
            vcf_path=vcf,
            dp_min=dp_min,
            maf_min=maf_min,
            dp_cap=dp_cap,
            included_contigs=included_contigs,
            variant_classes=[VariantClass.INS],
        )

        del_count, _ = count_ambiguous_sites(
            vcf_path=vcf,
            dp_min=dp_min,
            maf_min=maf_min,
            dp_cap=dp_cap,
            included_contigs=included_contigs,
            variant_classes=[VariantClass.DEL],
        )

        # Step 6: Calculate mapping rate from BAM
        calculate_mapping_rate(bam)

        # Step 7: Create summary row
        summary = SummaryRow(
            sample=sample_name,
            mode="summarize",
            mapper="unknown",  # Can't reliably determine from existing files
            caller="unknown",  # Can't reliably determine from existing files
            dp_min=dp_min,
            maf_min=maf_min,
            dp_cap=dp_cap,
            denom_policy="reused_bam",
            callable_bases=depth_summary.callable_bases,
            genome_length=depth_summary.genome_length,
            breadth_10x=depth_summary.breadth_10x,
            ambiguous_snv_count=ambiguous_snv_count,
            ambiguous_snv_per_mb=(
                ambiguous_snv_count / (depth_summary.callable_bases / 1_000_000)
                if depth_summary.callable_bases > 0
                else 0.0
            ),
            ambiguous_indel_count=indel_count,
            ambiguous_del_count=del_count,
            ref_label="unknown",
            ref_accession="unknown",
            bracken_species="NA",
            bracken_frac=0.0,
            bracken_reads=0,
            alias_used="NA",
            reused_bam=True,
            reused_vcf=True,
            runtime_sec=time.time() - start_time,
            tool_version=__version__,
        )

        # Step 8: Handle outputs
        if to_stdout:
            write_tsv_to_stdout([summary])
        else:
            # Write main summary TSV
            if not output:
                output = output_dir / "ambiguous_summary.tsv"
            write_tsv_summary(output, [summary], TSVMode.OVERWRITE)
            logger.info(f"Summary written to {output}")

        # Step 9: Emit optional outputs if requested
        if emit_vcf and not to_stdout:
            vcf_output = output_dir / f"{sample_name}.ambiguous_sites.vcf.gz"
            emit_vcf_file(
                normalized_vcf_path=vcf,
                output_path=vcf_output,
                dp_min=dp_min,
                maf_min=maf_min,
                sample_name=sample_name,
                require_pass=require_pass,
                included_contigs=included_contigs,
            )
            logger.info(f"VCF written to {vcf_output}")

        if emit_bed and not to_stdout:
            bed_output = output_dir / f"{sample_name}.ambiguous_sites.bed.gz"
            emit_bed_file(
                normalized_vcf_path=vcf,
                output_path=bed_output,
                dp_min=dp_min,
                maf_min=maf_min,
                sample_name=sample_name,
                included_contigs=included_contigs,
            )
            logger.info(f"BED written to {bed_output}")

        if emit_matrix and not to_stdout:
            matrix_output = output_dir / f"{sample_name}.variant_matrix.tsv.gz"
            emit_matrix_file(grid, matrix_output, sample_name)
            logger.info(f"Matrix written to {matrix_output}")

        if emit_per_contig and not to_stdout:
            per_contig_output = output_dir / f"{sample_name}.per_contig_summary.tsv"
            # Note: emit_per_contig implementation would need to be verified
            # For now, we'll skip this emission and just log that it was requested
            logger.info(f"Per-contig summary requested for {per_contig_output}")

        if emit_multiqc and not to_stdout:
            multiqc_output = output_dir / f"{sample_name}.multiqc.tsv"
            # Note: emit_multiqc implementation would need proper parameters
            # For now, we'll skip this emission and just log that it was requested
            logger.info(f"MultiQC metrics requested for {multiqc_output}")

        # Step 10: Handle provenance if requested
        if emit_provenance and not to_stdout:
            provenance_output = output_dir / "run_provenance.json"
            # Note: Full provenance implementation would go here
            logger.info(f"Provenance requested for {provenance_output}")

        # Clean up temporary depth directory
        import shutil

        if temp_depth_dir.exists():
            shutil.rmtree(temp_depth_dir, ignore_errors=True)

        logger.info(
            f"Summarize mode completed: {ambiguous_snv_count} ambiguous SNVs found"
        )
        return [summary]

    except Exception as e:
        logger.error(f"Summarize mode failed: {e}")
        raise


def execute_plan(
    plan: RunPlan,
    species: Optional[str] = None,
    bracken: Optional[Path] = None,
    ref_dir: Optional[Path] = None,
    on_fail: str = "error",
    min_bracken_frac: float = 0.1,
    min_bracken_reads: int = 10,
) -> Tuple[SummaryRow, Optional[ProvenanceRecord]]:
    """
    Execute a run plan and handle output.

    Args:
        plan: Validated execution plan
        **kwargs: Mode-specific arguments

    Returns:
        Tuple of (summary row with results, provenance record if emit_provenance is True)
    """
    from datetime import datetime

    # Handle dry run early - just execute the mode functions which will print plans and return placeholders
    if plan.dry_run:
        # Execute the appropriate mode (which will print dry-run info and return placeholder data)
        if plan.mode == Mode.SELF:
            result = run_self(plan)
        elif plan.mode == Mode.REF:
            result = run_ref(
                plan,
                species=species,
                bracken=bracken,
                ref_dir=ref_dir,
                on_fail=on_fail,
                min_bracken_frac=min_bracken_frac,
                min_bracken_reads=min_bracken_reads,
            )
        else:
            raise ValueError(f"Unsupported mode: {plan.mode}")

        # Return placeholder data without writing any files
        logger.info("Dry run completed - no files written")
        return result, None

    # Record start time
    start_time = datetime.now()

    # Get tool versions before execution
    tool_version = get_tool_versions(plan)

    # Execute the appropriate mode
    if plan.mode == Mode.SELF:
        result = run_self(plan)
    elif plan.mode == Mode.REF:
        result = run_ref(
            plan,
            species=species,
            bracken=bracken,
            ref_dir=ref_dir,
            on_fail=on_fail,
            min_bracken_frac=min_bracken_frac,
            min_bracken_reads=min_bracken_reads,
        )
    else:
        raise ValueError(f"Unsupported mode: {plan.mode}")

    # Record end time
    end_time = datetime.now()

    # Calculate runtime
    runtime_sec = (end_time - start_time).total_seconds()

    # Update the result with runtime and tool version
    result.runtime_sec = runtime_sec
    result.tool_version = tool_version

    # Handle output
    if plan.to_stdout:
        write_tsv_to_stdout([result])
    else:
        output_file = (
            plan.paths.output_dir / f"{plan.get_sample_prefix()}_ambiguous_summary.tsv"
        )
        write_tsv_summary(output_file, [result], plan.tsv_mode)
        logger.info(f"Results written to {output_file}")

    # Create provenance record if requested
    provenance_record = None
    if plan.emit_provenance:
        # Create mapping stats
        mapping_stats = ProvenanceMappingStats(
            map_rate=None,  # We don't store mapping_rate in SummaryRow currently
            breadth_10x=result.breadth_10x,
        )

        # Create counts
        counts = ProvenanceCounts(
            ambiguous_snv_count=result.ambiguous_snv_count,
            ambiguous_indel_count=result.ambiguous_indel_count,
            ambiguous_del_count=result.ambiguous_del_count,
            callable_bases=result.callable_bases,
            genome_length=result.genome_length,
        )

        # Extract reference path and species info
        reference_path = None
        reference_species = None
        if plan.mode == Mode.REF:
            # Try to get reference info from args or plan
            reference_path = plan.paths.assembly
            reference_species = species

        provenance_record = create_provenance_record(
            sample=result.sample,
            mode=plan.mode,
            started_at=start_time,
            finished_at=end_time,
            threads=plan.threads,
            mapper=plan.mapper.value,
            caller=plan.caller.value,
            dp_min=plan.thresholds.dp_min or 10,  # Use default if None
            maf_min=plan.thresholds.maf_min or 0.1,  # Use default if None
            dp_cap=plan.thresholds.dp_cap or 100,  # Use default if None
            denom_policy=result.denom_policy,  # Get from result, not plan
            depth_tool=plan.depth_tool.value,
            emit_vcf=plan.emit_vcf,
            emit_bed=plan.emit_bed,
            emit_matrix=plan.emit_matrix,
            emit_per_contig=plan.emit_per_contig,
            r1=plan.paths.r1,
            r2=plan.paths.r2,
            assembly=plan.paths.assembly,
            reference=plan.paths.reference,
            bam=plan.paths.bam,
            vcf=plan.paths.vcf,
            species=species,
            mapping_stats=mapping_stats,
            counts=counts,
            reference_path=reference_path,
            reference_species=reference_species,
        )

    return result, provenance_record
