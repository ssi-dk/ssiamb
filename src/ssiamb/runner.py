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
from .refdir import (
    resolve_reference_directory, resolve_species_to_fasta, 
    parse_accession_from_fasta_header, ReferenceResolutionError
)
from .mapping import (
    ensure_indexes_self, map_fastqs, check_external_tools,
    MappingError, ExternalToolError
)
from .depth import (
    analyze_depth, DepthAnalysisError, check_mosdepth_available
)
from .calling import (
    call_variants, VariantCallingError, check_caller_tools, VariantCallResult
)
from .vcf_ops import (
    normalize_and_split, count_ambiguous_sites, VCFOperationError,
    VariantClass, check_vcf_tools, emit_vcf, emit_bed, emit_matrix, emit_per_contig, emit_multiqc
)
from .mapping import calculate_mapping_rate
from .reuse import has_duplicate_flags, run_markdup_for_depth, CompatibilityError
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
    mapq: int = 30,
    dry_run: bool = False,
    to_stdout: bool = False,
    emit_vcf: bool = False,
    emit_bed: bool = False,
    emit_matrix: bool = False,
    emit_per_contig: bool = False,
    emit_multiqc: bool = False,
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
        # Extract output flags
        emit_vcf=emit_vcf,
        emit_bed=emit_bed,
        emit_matrix=emit_matrix,
        emit_per_contig=emit_per_contig,
        emit_provenance=kwargs.get("emit_provenance", False),
        emit_multiqc=emit_multiqc,
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
    
    # Initialize variables for error handling
    ambiguous_snv_count = 0
    all_variants_count = 0
    qc_warnings = ""
    
    # Ensure thresholds are properly initialized
    assert plan.thresholds.dp_min is not None, "dp_min threshold not set"
    assert plan.thresholds.maf_min is not None, "maf_min threshold not set"
    assert plan.paths.assembly is not None, "Assembly path required for self mode"
    
    if plan.dry_run:
        logger.info("DRY RUN - would execute self-mapping pipeline")
        logger.info(f"  Check tool availability: {plan.mapper.value}, samtools, mosdepth, {plan.caller.value}")
        logger.info(f"  Ensure indexes for {plan.paths.assembly} ({plan.mapper.value})")
        logger.info(f"  Map {plan.paths.r1} + {plan.paths.r2} to {plan.paths.assembly}")
        logger.info(f"  Run depth analysis with mosdepth (MAPQ≥{plan.thresholds.mapq_min}, depth≥{plan.thresholds.dp_min})")
        logger.info(f"  Call variants with {plan.caller.value} (MAPQ≥{plan.thresholds.mapq_min}, BASEQ≥{plan.thresholds.baseq_min})")
        logger.info(f"  Normalize VCF with bcftools norm (atomize and split multi-allelic sites)")
        logger.info(f"  Count ambiguous sites: dp_min={plan.thresholds.dp_min}, maf_min={plan.thresholds.maf_min}")
        logger.info(f"  Output: {plan.sample}.sorted.bam, {plan.sample}.mosdepth.summary.txt, {plan.sample}.vcf, {plan.sample}.normalized.vcf.gz")
        
        # Return placeholder data for dry run
        return SummaryRow(
            sample=plan.sample,
            mode=plan.mode.value,
            ref_label=plan.paths.assembly.name,
            mapper=plan.mapper.value,
            caller=plan.caller.value,
            dp_min=plan.thresholds.dp_min,
            maf_min=plan.thresholds.maf_min,
            callable_bases=2800000,
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
    
    # Check tool availability
    tools = check_external_tools()
    tool_name = plan.mapper.value.replace('-', '')  # minimap2 or bwamem2
    mosdepth_available = check_mosdepth_available()
    caller_available = check_caller_tools(plan.caller)
    
    if not tools.get(tool_name) or not tools.get('samtools'):
        missing = [t for t, avail in tools.items() if not avail and t in [tool_name, 'samtools']]
        raise ExternalToolError(f"Required tools not available: {missing}")
    
    if not mosdepth_available:
        raise ExternalToolError("mosdepth not found in PATH - required for depth analysis")
    
    if not caller_available:
        raise ExternalToolError(f"Variant caller {plan.caller.value} tools not found in PATH")
    
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
            output_path=plan.paths.output_dir / f"{plan.sample}.sorted.bam"
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
                logger.info("No duplicate flags detected, running samtools markdup for depth analysis")
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
                threads=plan.threads
            )
            logger.info(f"Depth analysis completed: {depth_summary.genome_length:,} bp genome, "
                       f"{depth_summary.breadth_10x:.2%} breadth ≥{plan.thresholds.dp_min}x")
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
            )
            logger.info(f"Variant calling completed: {variant_result.vcf_path}")
        except VariantCallingError as e:
            logger.error(f"Variant calling failed: {e}")
            raise ExternalToolError(f"Variant calling failed: {e}")
        
        # Step 5: Normalize VCF and count ambiguous sites
        logger.info("Normalizing VCF and counting ambiguous sites")
        try:
            # Normalize and split VCF
            normalized_vcf_path = normalize_and_split(
                vcf_in=variant_result.vcf_path,
                reference=plan.paths.assembly,
                output_dir=plan.paths.output_dir
            )
            
            # Get included contigs from depth analysis
            included_contigs = None
            if depth_summary:
                # Use the same contig filtering as depth analysis (≥500 bp)
                from .depth import list_included_contigs
                mosdepth_summary = plan.paths.output_dir / f"{plan.sample}.mosdepth.summary.txt"
                if mosdepth_summary.exists():
                    included_contigs = list_included_contigs(mosdepth_summary, min_len=500)
            
            # Count ambiguous sites (SNVs only for primary count)
            ambiguous_snv_count, grid = count_ambiguous_sites(
                vcf_path=normalized_vcf_path,
                dp_min=plan.thresholds.dp_min,
                maf_min=plan.thresholds.maf_min,
                dp_cap=100,  # As per spec
                included_contigs=included_contigs,
                variant_classes=[VariantClass.SNV]
            )
            
            # Count total variant sites by class
            all_variants_count, _ = count_ambiguous_sites(
                vcf_path=normalized_vcf_path,
                dp_min=0,  # Count all variants regardless of thresholds
                maf_min=0.0,
                dp_cap=100,
                included_contigs=included_contigs,
                variant_classes=[VariantClass.SNV, VariantClass.INS, VariantClass.DEL]
            )
            
            logger.info(f"Found {ambiguous_snv_count} ambiguous SNVs, {all_variants_count} total variant sites")
            
            # Optional emitters (only create files if flags are set)
            if plan.emit_vcf:
                vcf_output_path = plan.paths.output_dir / f"{plan.sample}.ambiguous_sites.vcf.gz"
                logger.info(f"Emitting VCF to {vcf_output_path}")
                emit_vcf(
                    normalized_vcf_path=normalized_vcf_path,
                    output_path=vcf_output_path,
                    dp_min=plan.thresholds.dp_min,
                    maf_min=plan.thresholds.maf_min,
                    sample_name=plan.sample,
                    require_pass=False,  # TODO: Add this to plan if needed
                    included_contigs=included_contigs
                )
            
            if plan.emit_bed:
                bed_output_path = plan.paths.output_dir / f"{plan.sample}.ambiguous_sites.bed.gz"
                logger.info(f"Emitting BED to {bed_output_path}")
                emit_bed(
                    normalized_vcf_path=normalized_vcf_path,
                    output_path=bed_output_path,
                    dp_min=plan.thresholds.dp_min,
                    maf_min=plan.thresholds.maf_min,
                    sample_name=plan.sample,
                    included_contigs=included_contigs
                )
            
            if plan.emit_matrix:
                matrix_output_path = plan.paths.output_dir / f"{plan.sample}.ambiguity_matrix.tsv.gz"
                logger.info(f"Emitting ambiguity matrix to {matrix_output_path}")
                emit_matrix(
                    grid=grid,
                    output_path=matrix_output_path,
                    sample_name=plan.sample
                )
            
            if plan.emit_per_contig:
                per_contig_output_path = plan.paths.output_dir / f"{plan.sample}.per_contig.tsv"
                logger.info(f"Emitting per-contig summary to {per_contig_output_path}")
                mosdepth_summary = plan.paths.output_dir / f"{plan.sample}.mosdepth.summary.txt"
                emit_per_contig(
                    depth_summary_path=mosdepth_summary,
                    normalized_vcf_path=normalized_vcf_path,
                    output_path=per_contig_output_path,
                    dp_min=plan.thresholds.dp_min,
                    maf_min=plan.thresholds.maf_min,
                    sample_name=plan.sample,
                    included_contigs=included_contigs
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
            mean_depth = depth_summary.mean_depth
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
            mean_depth = 0.0
            qc_warnings = qc_warnings or "depth_analysis_unavailable"
        
        # TODO: Add mapping rate calculation from BAM stats
        # Calculate ref_calls (sites that are reference-like, i.e., not variant)
        ref_calls = callable_bases - all_variants_count if callable_bases > all_variants_count else 0
        
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
                output_path=multiqc_output_path
            )
        
        return SummaryRow(
            sample=plan.sample,
            mode=plan.mode.value,
            ref_label=plan.paths.assembly.name,
            mapper=plan.mapper.value,
            caller=plan.caller.value,
            dp_min=plan.thresholds.dp_min,
            maf_min=plan.thresholds.maf_min,
            callable_bases=callable_bases,
            genome_length=genome_length,
            breadth_10x=breadth_10x,
            mean_depth=mean_depth,
            mapping_rate=mapping_rate,       # Calculated from BAM stats
            ambiguous_sites=ambiguous_snv_count,
            ref_calls=ref_calls,
            alt_calls=all_variants_count,
            no_call=genome_length - callable_bases if genome_length > callable_bases else 0,
            qc_warnings=qc_warnings,
        )
        
    except (ExternalToolError, MappingError) as e:
        logger.error(f"Self-mapping failed: {e}")
        raise
    except DepthAnalysisError as e:
        logger.error(f"Depth analysis failed: {e}")
        # Return summary with depth analysis failure
        return SummaryRow(
            sample=plan.sample,
            mode=plan.mode.value,
            ref_label=plan.paths.assembly.name if plan.paths.assembly else "unknown",
            mapper=plan.mapper.value,
            caller=plan.caller.value,
            dp_min=plan.thresholds.dp_min,
            maf_min=plan.thresholds.maf_min,
            callable_bases=0,
            genome_length=0,
            breadth_10x=0.0,
            mean_depth=0.0,
            mapping_rate=0.0,
            ambiguous_sites=0,
            ref_calls=0,
            alt_calls=0,
            no_call=0,
            qc_warnings=f"depth_analysis_failed:{e}",
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
    
    try:
        if plan.paths.reference:
            # Direct reference file provided
            ref_source = "file"
            ref_label = plan.paths.reference.name
            reference_path = plan.paths.reference
            
        elif kwargs.get("species"):
            # Species lookup in admin directory
            ref_source = "species"
            species = kwargs["species"]
            ref_dir = resolve_reference_directory(kwargs.get("ref_dir"))
            
            # Resolve species to FASTA file
            reference_path = resolve_species_to_fasta(ref_dir, species)
            ref_label = reference_path.name
            
            # Try to get accession for better labeling
            accession = parse_accession_from_fasta_header(reference_path)
            if accession:
                ref_label = f"{species} ({accession})"
            else:
                ref_label = species
                
        elif kwargs.get("bracken"):
            # Bracken-based species selection
            from .bracken import select_species_from_file, BrackenThresholds, BrackenSelectionError
            
            ref_source = "bracken"
            bracken_path = Path(kwargs["bracken"])
            
            # Get thresholds from configuration or defaults
            min_frac = kwargs.get("min_bracken_frac", 0.70)
            min_reads = kwargs.get("min_bracken_reads", 100000)
            thresholds = BrackenThresholds(min_frac=min_frac, min_reads=min_reads)
            
            try:
                # Select species from Bracken
                selection = select_species_from_file(bracken_path, thresholds)
                
                if selection is None:
                    raise ReferenceResolutionError(
                        f"No species in Bracken file {bracken_path} meet selection criteria. "
                        f"Try using --species or --reference instead."
                    )
                
                # Extract species name from Bracken (e.g., "Listeria monocytogenes" from potentially complex name)
                bracken_name = selection.species.name
                parts = bracken_name.strip().split()
                if len(parts) >= 2:
                    species_for_resolution = f"{parts[0]} {parts[1]}"  # "Listeria monocytogenes"
                else:
                    species_for_resolution = bracken_name
                
                # Resolve species to reference file
                ref_dir = resolve_reference_directory(kwargs.get("ref_dir"))
                reference_path = resolve_species_to_fasta(ref_dir, species_for_resolution)
                
                # Create detailed label with Bracken info
                from .refdir import normalize_species_name
                accession = parse_accession_from_fasta_header(reference_path)
                normalized_name = normalize_species_name(species_for_resolution)  # For label only
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
        on_fail_action = kwargs.get("on_fail", "error")
        
        if on_fail_action == "error":
            raise
        elif on_fail_action == "warn":
            logger.warning(f"Reference resolution failed: {e}")
            # Continue with placeholder data
            ref_source = "failed"
            ref_label = "resolution_failed"
        elif on_fail_action == "skip":
            logger.info(f"Skipping due to reference resolution failure: {e}")
            # Return a minimal result indicating skip
            return SummaryRow(
                sample=plan.sample,
                mode=plan.mode.value,
                ref_label="SKIPPED",
                mapper=plan.mapper.value,
                caller=plan.caller.value,
                dp_min=plan.thresholds.dp_min,
                maf_min=plan.thresholds.maf_min,
                callable_bases=0,
                genome_length=0,
                breadth_10x=0.0,
                mean_depth=0.0,
                mapping_rate=0.0,
                ambiguous_sites=0,
                ref_calls=0,
                alt_calls=0,
                no_call=0,
                qc_warnings="reference_resolution_failed",
            )
    
    if plan.dry_run:
        logger.info("DRY RUN - would execute reference-mapping pipeline")
        logger.info(f"  Reference source: {ref_source}")
        logger.info(f"  Reference: {reference_path}")
        logger.info(f"  Check tool availability: {plan.mapper.value}, samtools, mosdepth, {plan.caller.value}")
        logger.info(f"  Map {plan.paths.r1} + {plan.paths.r2} to reference")
        logger.info(f"  Run depth analysis with mosdepth (MAPQ≥{plan.thresholds.mapq_min}, depth≥{plan.thresholds.dp_min})")
        logger.info(f"  Call variants with {plan.caller.value} (MAPQ≥{plan.thresholds.mapq_min}, BASEQ≥{plan.thresholds.baseq_min})")
        logger.info(f"  Normalize VCF with bcftools norm (atomize and split multi-allelic sites)")
        logger.info(f"  Count ambiguous sites: dp_min={plan.thresholds.dp_min}, maf_min={plan.thresholds.maf_min}")
        
        # Return placeholder data for dry run
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
            qc_warnings="",
        )
    
    # Update the plan's reference info
    plan.ref_source = ref_source
    plan.ref_label = ref_label
    if reference_path:
        plan.paths.reference = reference_path
    
    # Check tool availability
    tools = check_external_tools()
    tool_name = plan.mapper.value.replace('-', '')  # minimap2 or bwamem2
    mosdepth_available = check_mosdepth_available()
    caller_available = check_caller_tools(plan.caller)
    
    if not tools.get(tool_name) or not tools.get('samtools'):
        missing = [t for t, avail in tools.items() if not avail and t in [tool_name, 'samtools']]
        raise ExternalToolError(f"Required tools not available: {missing}")
    
    if not mosdepth_available:
        raise ExternalToolError("mosdepth not found in PATH - required for depth analysis")
    
    if not caller_available:
        raise ExternalToolError(f"Variant caller {plan.caller.value} tools not found in PATH")
    
    # Only proceed if we have a valid reference
    if ref_source == "failed" or not reference_path:
        return SummaryRow(
            sample=plan.sample,
            mode=plan.mode.value,
            ref_label=ref_label,
            mapper=plan.mapper.value,
            caller=plan.caller.value,
            dp_min=plan.thresholds.dp_min,
            maf_min=plan.thresholds.maf_min,
            callable_bases=0,
            genome_length=0,
            breadth_10x=0.0,
            mean_depth=0.0,
            mapping_rate=0.0,
            ambiguous_sites=0,
            ref_calls=0,
            alt_calls=0,
            no_call=0,
            qc_warnings="reference_resolution_failed",
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
            output_path=plan.paths.output_dir / f"{plan.sample}.sorted.bam"
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
                logger.info("No duplicate flags detected, running samtools markdup for depth analysis")
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
                threads=plan.threads
            )
            logger.info(f"Depth analysis completed: {depth_summary.genome_length:,} bp genome, "
                       f"{depth_summary.breadth_10x:.2%} breadth ≥{plan.thresholds.dp_min}x")
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
                output_dir=plan.paths.output_dir
            )
            
            # Get included contigs from depth analysis
            included_contigs = None
            if depth_summary:
                # Use the same contig filtering as depth analysis (≥500 bp)
                from .depth import list_included_contigs
                mosdepth_summary = plan.paths.output_dir / f"{plan.sample}.mosdepth.summary.txt"
                if mosdepth_summary.exists():
                    included_contigs = list_included_contigs(mosdepth_summary, min_len=500)
            
            # Count ambiguous sites (SNVs only for primary count)
            ambiguous_snv_count, grid = count_ambiguous_sites(
                vcf_path=normalized_vcf_path,
                dp_min=plan.thresholds.dp_min,
                maf_min=plan.thresholds.maf_min,
                dp_cap=100,  # As per spec
                included_contigs=included_contigs,
                variant_classes=[VariantClass.SNV]
            )
            
            # Count total variant sites by class
            all_variants_count, _ = count_ambiguous_sites(
                vcf_path=normalized_vcf_path,
                dp_min=0,  # Count all variants regardless of thresholds
                maf_min=0.0,
                dp_cap=100,
                included_contigs=included_contigs,
                variant_classes=[VariantClass.SNV, VariantClass.INS, VariantClass.DEL]
            )
            
            logger.info(f"Found {ambiguous_snv_count} ambiguous SNVs, {all_variants_count} total variant sites")
            
            # Optional emitters (only create files if flags are set)
            if plan.emit_vcf:
                vcf_output_path = plan.paths.output_dir / f"{plan.sample}.ambiguous_sites.vcf.gz"
                logger.info(f"Emitting VCF to {vcf_output_path}")
                emit_vcf(
                    normalized_vcf_path=normalized_vcf_path,
                    output_path=vcf_output_path,
                    dp_min=plan.thresholds.dp_min,
                    maf_min=plan.thresholds.maf_min,
                    sample_name=plan.sample,
                    require_pass=False,  # TODO: Add this to plan if needed
                    included_contigs=included_contigs
                )
            
            if plan.emit_bed:
                bed_output_path = plan.paths.output_dir / f"{plan.sample}.ambiguous_sites.bed.gz"
                logger.info(f"Emitting BED to {bed_output_path}")
                emit_bed(
                    normalized_vcf_path=normalized_vcf_path,
                    output_path=bed_output_path,
                    dp_min=plan.thresholds.dp_min,
                    maf_min=plan.thresholds.maf_min,
                    sample_name=plan.sample,
                    included_contigs=included_contigs
                )
            
            if plan.emit_matrix:
                matrix_output_path = plan.paths.output_dir / f"{plan.sample}.ambiguity_matrix.tsv.gz"
                logger.info(f"Emitting ambiguity matrix to {matrix_output_path}")
                emit_matrix(
                    grid=grid,
                    output_path=matrix_output_path,
                    sample_name=plan.sample
                )
            
            if plan.emit_per_contig:
                per_contig_output_path = plan.paths.output_dir / f"{plan.sample}.per_contig.tsv"
                logger.info(f"Emitting per-contig summary to {per_contig_output_path}")
                mosdepth_summary = plan.paths.output_dir / f"{plan.sample}.mosdepth.summary.txt"
                emit_per_contig(
                    depth_summary_path=mosdepth_summary,
                    normalized_vcf_path=normalized_vcf_path,
                    output_path=per_contig_output_path,
                    dp_min=plan.thresholds.dp_min,
                    maf_min=plan.thresholds.maf_min,
                    sample_name=plan.sample,
                    included_contigs=included_contigs
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
            mean_depth = depth_summary.mean_depth
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
            mean_depth = 0.0
            qc_warnings = qc_warnings or "depth_analysis_unavailable"
        
        # TODO: Add mapping rate calculation from BAM stats
        # Calculate ref_calls (sites that are reference-like, i.e., not variant)
        ref_calls = callable_bases - all_variants_count if callable_bases > all_variants_count else 0
        
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
                output_path=multiqc_output_path
            )
        
        return SummaryRow(
            sample=plan.sample,
            mode=plan.mode.value,
            ref_label=ref_label,
            mapper=plan.mapper.value,
            caller=plan.caller.value,
            dp_min=plan.thresholds.dp_min,
            maf_min=plan.thresholds.maf_min,
            callable_bases=callable_bases,
            genome_length=genome_length,
            breadth_10x=breadth_10x,
            mean_depth=mean_depth,
            mapping_rate=mapping_rate,       # Calculated from BAM stats
            ambiguous_sites=ambiguous_snv_count,
            ref_calls=ref_calls,
            alt_calls=all_variants_count,
            no_call=genome_length - callable_bases if genome_length > callable_bases else 0,
            qc_warnings=qc_warnings,
        )
        
    except (ExternalToolError, MappingError) as e:
        logger.error(f"Reference-mapping failed: {e}")
        raise
    except DepthAnalysisError as e:
        logger.error(f"Depth analysis failed: {e}")
        # Return summary with depth analysis failure
        return SummaryRow(
            sample=plan.sample,
            mode=plan.mode.value,
            ref_label=ref_label,
            mapper=plan.mapper.value,
            caller=plan.caller.value,
            dp_min=plan.thresholds.dp_min,
            maf_min=plan.thresholds.maf_min,
            callable_bases=0,
            genome_length=0,
            breadth_10x=0.0,
            mean_depth=0.0,
            mapping_rate=0.0,
            ambiguous_sites=0,
            ref_calls=0,
            alt_calls=0,
            no_call=0,
            qc_warnings=f"depth_analysis_failed:{e}",
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