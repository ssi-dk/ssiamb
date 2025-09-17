"""
Variant calling module for ssiamb.

This module implements variant calling using BBTools and bcftools pipelines.
Supports both bacterial genome analysis with appropriate ploidy and quality settings.
"""

import logging
import subprocess
import shutil
from pathlib import Path
from typing import Optional, List, Tuple
from dataclasses import dataclass

from .models import Caller

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
        return (shutil.which("pileup.sh") is not None and 
                shutil.which("callvariants.sh") is not None)
    elif caller == Caller.BCFTOOLS:
        # Check for bcftools
        return shutil.which("bcftools") is not None
    else:
        raise ValueError(f"Unknown caller: {caller}")


def run_bbtools_calling(
    bam_path: Path,
    reference_path: Path,
    output_vcf: Path,
    sample_name: str,
    threads: int = 1,
    mapq_min: int = 20,
    baseq_min: int = 20,
    minallelefraction: float = 0.0,
) -> VariantCallResult:
    """
    Run BBTools variant calling pipeline.
    
    Executes:
    1. pileup.sh minmapq=X
    2. callvariants.sh ploidy=1 clearfilters=t minallelefraction=Z minavgmapq=X minquality=Y
    
    Args:
        bam_path: Input BAM file
        reference_path: Reference genome FASTA
        output_vcf: Output VCF file path
        sample_name: Sample name for VCF header
        threads: Number of threads to use
        mapq_min: Minimum mapping quality
        baseq_min: Minimum base quality
        minallelefraction: Minimum allele fraction for variant calling
        
    Returns:
        VariantCallResult with execution details
    """
    import time
    start_time = time.time()
    
    try:
        # Ensure output directory exists
        output_vcf.parent.mkdir(parents=True, exist_ok=True)
        
        # Step 1: Run pileup.sh
        pileup_out = output_vcf.with_suffix('.var')
        
        pileup_cmd = [
            "pileup.sh",
            f"in={bam_path}",
            f"ref={reference_path}",
            f"out={pileup_out}",
            f"minmapq={mapq_min}",
            f"threads={threads}",
        ]
        
        logger.info(f"Running BBTools pileup: {' '.join(map(str, pileup_cmd))}")
        
        result = subprocess.run(
            pileup_cmd,
            capture_output=True,
            text=True,
            check=False
        )
        
        if result.returncode != 0:
            error_msg = f"BBTools pileup failed (exit {result.returncode}): {result.stderr}"
            logger.error(error_msg)
            return VariantCallResult(
                vcf_path=output_vcf,
                caller=Caller.BBTOOLS,
                success=False,
                error_message=error_msg,
                runtime_seconds=time.time() - start_time
            )
        
        # Step 2: Run callvariants.sh
        callvariants_cmd = [
            "callvariants.sh",
            f"in={pileup_out}",
            f"ref={reference_path}",
            f"out={output_vcf}",
            "ploidy=1",
            "clearfilters=t",  # Clear all filters to get raw variants
            f"minallelefraction={minallelefraction}",
            f"minavgmapq={mapq_min}",
            f"minquality={baseq_min}",
            f"threads={threads}",
        ]
        
        logger.info(f"Running BBTools callvariants: {' '.join(map(str, callvariants_cmd))}")
        
        result = subprocess.run(
            callvariants_cmd,
            capture_output=True,
            text=True,
            check=False
        )
        
        if result.returncode != 0:
            error_msg = f"BBTools callvariants failed (exit {result.returncode}): {result.stderr}"
            logger.error(error_msg)
            return VariantCallResult(
                vcf_path=output_vcf,
                caller=Caller.BBTOOLS,
                success=False,
                error_message=error_msg,
                runtime_seconds=time.time() - start_time
            )
        
        # Clean up intermediate file
        try:
            if pileup_out.exists():
                pileup_out.unlink()
        except OSError:
            # Ignore cleanup errors
            pass
        
        # Verify output VCF was created
        if not output_vcf.exists():
            error_msg = "BBTools callvariants completed but no VCF output found"
            logger.error(error_msg)
            return VariantCallResult(
                vcf_path=output_vcf,
                caller=Caller.BBTOOLS,
                success=False,
                error_message=error_msg,
                runtime_seconds=time.time() - start_time
            )
        
        logger.info(f"BBTools variant calling completed: {output_vcf}")
        return VariantCallResult(
            vcf_path=output_vcf,
            caller=Caller.BBTOOLS,
            success=True,
            runtime_seconds=time.time() - start_time
        )
        
    except Exception as e:
        error_msg = f"BBTools variant calling failed with exception: {e}"
        logger.error(error_msg)
        return VariantCallResult(
            vcf_path=output_vcf,
            caller=Caller.BBTOOLS,
            success=False,
            error_message=error_msg,
            runtime_seconds=time.time() - start_time
        )


def run_bcftools_calling(
    bam_path: Path,
    reference_path: Path,
    output_vcf: Path,
    sample_name: str,
    threads: int = 1,
    mapq_min: int = 20,
    baseq_min: int = 20,
) -> VariantCallResult:
    """
    Run bcftools variant calling pipeline.
    
    Executes:
    1. bcftools mpileup -q20 -Q20 -B -a AD,ADF,ADR,DP
    2. bcftools call -m --ploidy 1
    
    Args:
        bam_path: Input BAM file
        reference_path: Reference genome FASTA
        output_vcf: Output VCF file path
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
        
        # Combine mpileup and call in a pipeline
        mpileup_cmd = [
            "bcftools", "mpileup",
            f"--threads", str(threads),
            f"-q", str(mapq_min),      # Minimum mapping quality
            f"-Q", str(baseq_min),     # Minimum base quality
            "-B",                      # Disable BAQ computation
            "-a", "AD,ADF,ADR,DP",     # Annotations to include
            "-f", str(reference_path),  # Reference FASTA
            str(bam_path)              # Input BAM
        ]
        
        call_cmd = [
            "bcftools", "call",
            f"--threads", str(threads),
            "-m",                      # Multiallelic caller
            "--ploidy", "1",           # Haploid organism
            "-v",                      # Output only variants
            "-o", str(output_vcf)      # Output file
        ]
        
        logger.info(f"Running bcftools pipeline: {' '.join(mpileup_cmd)} | {' '.join(call_cmd)}")
        
        # Run mpileup | call pipeline
        mpileup_proc = subprocess.Popen(
            mpileup_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
        call_proc = subprocess.Popen(
            call_cmd,
            stdin=mpileup_proc.stdout,
            stderr=subprocess.PIPE,
            text=True
        )
        
        # Close mpileup stdout to allow it to receive SIGPIPE
        if mpileup_proc.stdout:
            mpileup_proc.stdout.close()
        
        # Wait for both processes
        call_stderr = call_proc.communicate()[1]
        mpileup_stderr = mpileup_proc.communicate()[1]
        
        # Check return codes
        if mpileup_proc.returncode != 0:
            error_msg = f"bcftools mpileup failed (exit {mpileup_proc.returncode}): {mpileup_stderr}"
            logger.error(error_msg)
            return VariantCallResult(
                vcf_path=output_vcf,
                caller=Caller.BCFTOOLS,
                success=False,
                error_message=error_msg,
                runtime_seconds=time.time() - start_time
            )
        
        if call_proc.returncode != 0:
            error_msg = f"bcftools call failed (exit {call_proc.returncode}): {call_stderr}"
            logger.error(error_msg)
            return VariantCallResult(
                vcf_path=output_vcf,
                caller=Caller.BCFTOOLS,
                success=False,
                error_message=error_msg,
                runtime_seconds=time.time() - start_time
            )
        
        # Verify output VCF was created
        if not output_vcf.exists():
            error_msg = "bcftools call completed but no VCF output found"
            logger.error(error_msg)
            return VariantCallResult(
                vcf_path=output_vcf,
                caller=Caller.BCFTOOLS,
                success=False,
                error_message=error_msg,
                runtime_seconds=time.time() - start_time
            )
        
        logger.info(f"bcftools variant calling completed: {output_vcf}")
        return VariantCallResult(
            vcf_path=output_vcf,
            caller=Caller.BCFTOOLS,
            success=True,
            runtime_seconds=time.time() - start_time
        )
        
    except Exception as e:
        error_msg = f"bcftools variant calling failed with exception: {e}"
        logger.error(error_msg)
        return VariantCallResult(
            vcf_path=output_vcf,
            caller=Caller.BCFTOOLS,
            success=False,
            error_message=error_msg,
            runtime_seconds=time.time() - start_time
        )


def call_variants(
    bam_path: Path,
    reference_path: Path,
    output_vcf: Path,
    caller: Caller,
    sample_name: str,
    threads: int = 1,
    mapq_min: int = 20,
    baseq_min: int = 20,
    minallelefraction: float = 0.0,
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
        
    Returns:
        VariantCallResult with execution details
        
    Raises:
        VariantCallingError: If caller tools are not available or calling fails
    """
    # Check tool availability
    if not check_caller_tools(caller):
        raise VariantCallingError(f"Required tools for {caller.value} are not available")
    
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
        )
    elif caller == Caller.BCFTOOLS:
        result = run_bcftools_calling(
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
        raise VariantCallingError(result.error_message or f"{caller.value} variant calling failed")
    
    return result


def get_available_callers() -> List[Caller]:
    """
    Get list of available variant callers based on tool availability.
    
    Returns:
        List of available Caller enum values
    """
    available = []
    
    for caller in Caller:
        if check_caller_tools(caller):
            available.append(caller)
    
    return available