"""
Variant calling module for ssiamb.

This module implements variant calling using BBTools and bcftools pipelines.
Supports both bacterial genome analysis with appropriate ploidy and quality settings.
"""

import logging
import subprocess
import shutil
from pathlib import Path
from typing import Optional, List, Dict, Any
from dataclasses import dataclass

from .models import Caller
from .tool_config import get_tool_path, get_tool_path_optional
from .config import get_config

logger = logging.getLogger(__name__)


class VariantCallingError(Exception):
    """Raised when variant calling fails."""

    pass


@dataclass
class VariantCallResult:
    """Result of variant calling operation."""

    vcf_path: Path
    caller: Caller
    success: bool
    error_message: Optional[str] = None
    runtime_seconds: Optional[float] = None


def check_caller_tools_detailed(caller: Caller) -> Dict[str, Dict[str, Any]]:
    """
    Check availability and versions of required tools for the specified caller.

    Args:
        caller: Variant caller to check

    Returns:
        Dictionary mapping tool names to availability and version info.
        Each tool entry contains:
        - 'available': boolean status
        - 'version': version string if available, or error message
    """
    tools = {}

    if caller == Caller.BBTOOLS:
        # Check for BBTools executables
        for tool in ["pileup.sh", "callvariants.sh"]:
            tools[tool] = {"available": False, "version": "unknown"}

            tool_path = get_tool_path_optional(tool)
            if tool_path:
                tools[tool]["available"] = True
                tools[tool]["version"] = "BBTools (version check not supported)"
                logger.debug(f"Tool {tool}: found")
            else:
                tools[tool]["version"] = "not found in PATH"
                logger.debug(f"Tool {tool}: not found")
            logger.debug(f"Tool {tool}: available")

    elif caller == Caller.BCFTOOLS:
        # Check for bcftools
        tools["bcftools"] = {"available": False, "version": "unknown"}

        if shutil.which("bcftools") is None:
            tools["bcftools"]["version"] = "not found in PATH"
            logger.debug("Tool bcftools: not found")
        else:
            # Try to get bcftools version
            try:
                result = subprocess.run(
                    ["bcftools", "--version"], capture_output=True, text=True, timeout=5
                )
                if result.returncode == 0:
                    tools["bcftools"]["available"] = True
                    # Extract version from output (first line typically)
                    version_line = (
                        result.stdout.strip().split("\n")[0]
                        if result.stdout
                        else "version check succeeded"
                    )
                    tools["bcftools"]["version"] = version_line
                    logger.debug(f"Tool bcftools: available, version: {version_line}")
                else:
                    tools["bcftools"][
                        "version"
                    ] = f"version check failed (exit {result.returncode})"
                    logger.debug("Tool bcftools: found but version check failed")
            except subprocess.TimeoutExpired:
                tools["bcftools"]["version"] = "version check timed out"
                logger.debug("Tool bcftools: found but version check timed out")
            except Exception as e:
                tools["bcftools"]["version"] = f"version check error: {e}"
                logger.debug(f"Tool bcftools: found but version check error: {e}")
    else:
        raise ValueError(f"Unknown caller: {caller}")

    return tools


def check_caller_tools(caller: Caller) -> bool:
    """
    Check if required tools for the specified caller are available.

    Args:
        caller: Variant caller to check

    Returns:
        True if all required tools are available, False otherwise
    """
    if caller == Caller.BBTOOLS:
        # Check for BBTools executables
        return (
            get_tool_path_optional("pileup.sh") is not None
            and get_tool_path_optional("callvariants.sh") is not None
        )
    elif caller == Caller.BCFTOOLS:
        # Check for bcftools
        return get_tool_path_optional("bcftools") is not None
    else:
        raise ValueError(f"Unknown caller: {caller}")


def caller_tools_available(caller: Caller) -> bool:
    """
    Check if all required tools for the specified caller are available.

    Args:
        caller: Variant caller to check

    Returns:
        True if all required tools are available, False otherwise
    """
    tools = check_caller_tools_detailed(caller)
    return all(tool_info.get("available", False) for tool_info in tools.values())


def run_bbtools_calling(
    bam_path: Path,
    reference_path: Path,
    output_vcf: Path,
    sample_name: str,
    threads: Optional[int] = None,
    mapq_min: Optional[int] = None,
    baseq_min: Optional[int] = None,
    minallelefraction: Optional[float] = None,
    bbtools_mem: Optional[str] = None,
) -> VariantCallResult:
    """
    Run BBTools variant calling pipeline.

    Executes:
    1. callvariants.sh directly with BAM input (no pileup step needed)

    Args:
        bam_path: Input BAM file
        reference_path: Reference genome FASTA
        output_vcf: Output VCF file path
        sample_name: Sample name for VCF header
        threads: Number of threads to use
        mapq_min: Minimum mapping quality
        baseq_min: Minimum base quality
        minallelefraction: Minimum allele fraction for variant calling
        bbtools_mem: BBTools heap memory (e.g., '4g', '8g')

    Returns:
        VariantCallResult with execution details
    """
    import time

    # Get config defaults for None parameters, but only if clearfilters is not active
    config = get_config()
    bbtools_config = config.tools.get("bbtools", {})
    clearfilters_enabled = bbtools_config.get("clearfilters", False)

    if threads is None:
        threads = config.get_calling_setting("default_threads", 1)

    # If clearfilters is on, we want mapq/baseq to be None so they can be set to 0 later
    # Otherwise, if they are None, set them to their default values.
    if not clearfilters_enabled:
        if mapq_min is None:
            mapq_min = config.get_calling_setting("default_mapq_min", 20)
        if baseq_min is None:
            baseq_min = config.get_calling_setting("default_baseq_min", 20)

    if minallelefraction is None:
        minallelefraction = config.get_calling_setting("default_minallelefraction", 0.0)

    start_time = time.time()

    try:
        # Ensure output directory exists
        output_vcf.parent.mkdir(parents=True, exist_ok=True)

        # Run callvariants.sh directly with BAM input (simpler approach)
        callvariants_path = get_tool_path("callvariants.sh")
        bbtools_config = config.tools.get("bbtools", {})

        # Use configured parameters with fallbacks to passed values
        ploidy = bbtools_config.get("ploidy", 1)
        clearfilters_enabled = bbtools_config.get("clearfilters", False)
        clearfilters = "t" if clearfilters_enabled else "f"
        configured_maf = bbtools_config.get("minallelefraction", 0.0)
        configured_mapq = bbtools_config.get("minavgmapq", 20)
        configured_baseq = bbtools_config.get("minquality", 20)

        # Use passed parameters if provided, otherwise use configured values
        final_maf = (
            minallelefraction if minallelefraction is not None else configured_maf
        )

        # When clearfilters is enabled and parameters are None (from runner),
        # use 0 to truly clear filters instead of config defaults
        if clearfilters_enabled:
            final_mapq = mapq_min if mapq_min is not None else 0
            final_baseq = baseq_min if baseq_min is not None else 0
        else:
            final_mapq = mapq_min if mapq_min is not None else configured_mapq
            final_baseq = baseq_min if baseq_min is not None else configured_baseq

        callvariants_cmd = [
            callvariants_path,
            f"in={bam_path}",
            f"ref={reference_path}",
            f"vcf={output_vcf}",
            f"ploidy={ploidy}",
            f"clearfilters={clearfilters}",
            f"threads={threads}",
        ]

        # Only add quality filter parameters if clearfilters is not enabled
        # When clearfilters=t, BBTools ignores filters placed after it,
        # but we want to honor the clearfilters intent to truly clear all filters
        if not clearfilters_enabled:
            callvariants_cmd.extend(
                [
                    f"minallelefraction={final_maf}",
                    f"minavgmapq={final_mapq}",
                    f"minquality={final_baseq}",
                ]
            )
        else:
            # When clearfilters is enabled, only add explicit non-zero parameters
            # to allow users to selectively re-enable specific filters after clearfilters
            if final_maf > 0:
                callvariants_cmd.append(f"minallelefraction={final_maf}")
            if final_mapq > 0:
                callvariants_cmd.append(f"minavgmapq={final_mapq}")
            if final_baseq > 0:
                callvariants_cmd.append(f"minquality={final_baseq}")

        # Add memory setting if provided
        if bbtools_mem:
            callvariants_cmd.append(f"-Xmx{bbtools_mem}")

        # Add extra BBTools args from config
        bbtools_extra_args = config.tools.get("bbtools", {}).get("extra_args", [])
        callvariants_cmd.extend(bbtools_extra_args)

        logger.info(
            f"Running BBTools callvariants: {' '.join(map(str, callvariants_cmd))}"
        )

        result = subprocess.run(
            callvariants_cmd,
            capture_output=True,
            text=True,
            check=False,
            timeout=3600,  # 1 hour timeout
        )

        if result.returncode != 0:
            error_msg = f"BBTools callvariants failed (exit {result.returncode}): {result.stderr}"
            logger.error(error_msg)
            return VariantCallResult(
                vcf_path=output_vcf,
                caller=Caller.BBTOOLS,
                success=False,
                error_message=error_msg,
                runtime_seconds=time.time() - start_time,
            )

        # Log stdout for debugging
        if result.stdout:
            logger.debug(f"BBTools callvariants stdout: {result.stdout.strip()}")

        # Verify output VCF was created
        if not output_vcf.exists() or output_vcf.stat().st_size == 0:
            error_msg = "BBTools callvariants completed but no VCF output found"
            logger.error(error_msg)
            return VariantCallResult(
                vcf_path=output_vcf,
                caller=Caller.BBTOOLS,
                success=False,
                error_message=error_msg,
                runtime_seconds=time.time() - start_time,
            )

        logger.info(f"BBTools variant calling completed: {output_vcf}")
        return VariantCallResult(
            vcf_path=output_vcf,
            caller=Caller.BBTOOLS,
            success=True,
            runtime_seconds=time.time() - start_time,
        )

    except subprocess.TimeoutExpired:
        error_msg = "BBTools variant calling timed out"
        logger.error(error_msg)
        return VariantCallResult(
            vcf_path=output_vcf,
            caller=Caller.BBTOOLS,
            success=False,
            error_message=error_msg,
            runtime_seconds=time.time() - start_time,
        )
    except Exception as e:
        error_msg = f"BBTools variant calling failed with exception: {e}"
        logger.error(error_msg)
        return VariantCallResult(
            vcf_path=output_vcf,
            caller=Caller.BBTOOLS,
            success=False,
            error_message=error_msg,
            runtime_seconds=time.time() - start_time,
        )


def run_bcftools_pileup(
    bam_path: Path,
    reference_path: Path,
    output_vcf: Path,
    sample_name: str,
    threads: Optional[int] = None,
    mapq_min: Optional[int] = None,
    baseq_min: Optional[int] = None,
) -> VariantCallResult:
    """
    Run BCFtools pileup generation pipeline.

    Executes:
    1. bcftools mpileup -q20 -Q20 -B -a AD,ADF,ADR,DP
    2. bcftools view -Oz (keeps all sites, compress to VCF.gz)

    This creates candidate sites (pileup) rather than called variants,
    suitable for downstream VCF normalization and ambiguous site detection.
    Normalization is now handled separately in step #4.

    Args:
        bam_path: Input BAM file
        reference_path: Reference genome FASTA
        output_vcf: Output VCF file path (will be .pileup.vcf.gz)
        sample_name: Sample name for VCF header
        threads: Number of threads to use
        mapq_min: Minimum mapping quality
        baseq_min: Minimum base quality

    Returns:
        VariantCallResult with execution details
    """
    import time

    start_time = time.time()

    try:
        # Ensure output directory exists
        output_vcf.parent.mkdir(parents=True, exist_ok=True)

        # Convert output to .pileup.vcf.gz naming convention
        if not str(output_vcf).endswith(".pileup.vcf.gz"):
            pileup_vcf = output_vcf.with_suffix(".pileup.vcf.gz")
        else:
            pileup_vcf = output_vcf

        # Get BCFtools configuration
        config = get_config()
        bcftools_config = config.tools.get("bcftools", {})

        # Use configured parameters with fallbacks to passed values
        configured_mapq = bcftools_config.get("mapq_min", 20)
        configured_baseq = bcftools_config.get("baseq_min", 20)

        # Use passed parameters if provided, otherwise use configured values
        final_mapq = mapq_min if mapq_min is not None else configured_mapq
        final_baseq = baseq_min if baseq_min is not None else configured_baseq

        # Build mpileup command using configuration
        bcftools_path = get_tool_path("bcftools")
        mpileup_cmd = [bcftools_path, "mpileup"]

        # Add output type from config
        output_type = bcftools_config.get("mpileup_output_type", "-Ou")
        mpileup_cmd.append(output_type)

        # Add standard arguments
        mpileup_cmd.extend(
            [
                "--threads",
                str(threads),
                "-q",
                str(final_mapq),  # Minimum mapping quality
                "-Q",
                str(final_baseq),  # Minimum base quality
            ]
        )

        # Add BAQ disable flag if configured
        if bcftools_config.get("disable_baq", True):
            mpileup_cmd.append("-B")

        # Add max depth from config
        max_depth = bcftools_config.get("max_depth", 100000)
        mpileup_cmd.extend(["--max-depth", str(max_depth)])

        # Add format annotations from config
        annotations = bcftools_config.get("annotations", "FORMAT/AD,ADF,ADR,DP")
        mpileup_cmd.extend(["-a", annotations])

        # Add existing mpileup args from config
        mpileup_args = bcftools_config.get("mpileup_args", [])
        mpileup_cmd.extend(mpileup_args)

        # Add reference and BAM files
        mpileup_cmd.extend(
            [
                "-f",
                str(reference_path),  # Reference FASTA
                str(bam_path),  # Input BAM
            ]
        )

        # Build view command using configuration (outputs compressed VCF directly)
        view_cmd = [bcftools_path, "view"]
        view_cmd.extend(["--threads", str(threads)])

        # Add view output type from config (should be -Oz for compressed VCF)
        view_output_type = bcftools_config.get("view_output_type", "-Oz")
        view_cmd.append(view_output_type)

        # Add output file directly to view command
        view_cmd.extend(["-o", str(pileup_vcf)])

        # Add view args from config
        view_args = bcftools_config.get("view_args", [])
        view_cmd.extend(view_args)

        logger.info(
            f"Running bcftools pileup pipeline: {' '.join(mpileup_cmd)} | {' '.join(view_cmd)}"
        )

        # Run mpileup | view pipeline
        mpileup_proc = subprocess.Popen(
            mpileup_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )

        view_proc = subprocess.Popen(
            view_cmd,
            stdin=mpileup_proc.stdout,
            stderr=subprocess.PIPE,
            text=True,
        )

        # Close upstream stdout to allow SIGPIPE
        if mpileup_proc.stdout:
            mpileup_proc.stdout.close()

        # Wait for both processes
        view_stderr = view_proc.communicate()[1]
        mpileup_stderr = mpileup_proc.communicate()[1]

        # Check return codes
        if mpileup_proc.returncode != 0:
            error_msg = f"bcftools mpileup failed (exit {mpileup_proc.returncode}): {mpileup_stderr}"
            logger.error(error_msg)
            return VariantCallResult(
                vcf_path=pileup_vcf,
                caller=Caller.BCFTOOLS,
                success=False,
                error_message=error_msg,
                runtime_seconds=time.time() - start_time,
            )

        if view_proc.returncode != 0:
            error_msg = (
                f"bcftools view failed (exit {view_proc.returncode}): {view_stderr}"
            )
            logger.error(error_msg)
            return VariantCallResult(
                vcf_path=pileup_vcf,
                caller=Caller.BCFTOOLS,
                success=False,
                error_message=error_msg,
                runtime_seconds=time.time() - start_time,
            )

        # Verify output VCF was created
        if not pileup_vcf.exists():
            error_msg = "bcftools view completed but no VCF output found"
            logger.error(error_msg)
            return VariantCallResult(
                vcf_path=pileup_vcf,
                caller=Caller.BCFTOOLS,
                success=False,
                error_message=error_msg,
                runtime_seconds=time.time() - start_time,
            )

        logger.info(f"bcftools pileup generation completed: {pileup_vcf}")
        return VariantCallResult(
            vcf_path=pileup_vcf,
            caller=Caller.BCFTOOLS,
            success=True,
            runtime_seconds=time.time() - start_time,
        )

    except Exception as e:
        error_msg = f"bcftools pileup generation failed with exception: {e}"
        logger.error(error_msg)
        return VariantCallResult(
            vcf_path=output_vcf,
            caller=Caller.BCFTOOLS,
            success=False,
            error_message=error_msg,
            runtime_seconds=time.time() - start_time,
        )


def call_variants(
    bam_path: Path,
    reference_path: Path,
    output_vcf: Path,
    caller: Caller,
    sample_name: str,
    threads: int = 1,
    mapq_min: Optional[int] = None,
    baseq_min: Optional[int] = None,
    minallelefraction: float = 0.0,
    bbtools_mem: Optional[str] = None,
) -> VariantCallResult:
    """
    Run variant calling with the specified caller.

    Args:
        bam_path: Input BAM file
        reference_path: Reference genome FASTA
        output_vcf: Output VCF file path
        caller: Variant caller to use
        sample_name: Sample name for VCF header
        threads: Number of threads to use
        mapq_min: Minimum mapping quality
        baseq_min: Minimum base quality
        minallelefraction: Minimum allele fraction (BBTools only)
        bbtools_mem: BBTools heap memory (e.g., '4g', '8g')

    Returns:
        VariantCallResult with execution details

    Raises:
        VariantCallingError: If caller tools are not available or calling fails
    """
    # Check tool availability
    if not caller_tools_available(caller):
        raise VariantCallingError(
            f"Required tools for {caller.value} are not available"
        )

    # Validate inputs
    if not bam_path.exists():
        raise VariantCallingError(f"BAM file not found: {bam_path}")
    if not reference_path.exists():
        raise VariantCallingError(f"Reference file not found: {reference_path}")

    # Run appropriate caller
    if caller == Caller.BBTOOLS:
        result = run_bbtools_calling(
            bam_path=bam_path,
            reference_path=reference_path,
            output_vcf=output_vcf,
            sample_name=sample_name,
            threads=threads,
            mapq_min=mapq_min,
            baseq_min=baseq_min,
            minallelefraction=minallelefraction,
            bbtools_mem=bbtools_mem,
        )
    elif caller == Caller.BCFTOOLS:
        result = run_bcftools_pileup(
            bam_path=bam_path,
            reference_path=reference_path,
            output_vcf=output_vcf,
            sample_name=sample_name,
            threads=threads,
            mapq_min=mapq_min,
            baseq_min=baseq_min,
        )
    else:
        raise VariantCallingError(f"Unsupported caller: {caller}")

    # Raise exception if calling failed
    if not result.success:
        raise VariantCallingError(
            result.error_message or f"{caller.value} variant calling failed"
        )

    return result


def get_available_callers() -> List[Caller]:
    """
    Get list of available variant callers based on tool availability.

    Returns:
        List of available Caller enum values
    """
    available = []

    for caller in Caller:
        if caller_tools_available(caller):
            available.append(caller)

    return available
