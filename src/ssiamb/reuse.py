"""
Reuse and compatibility checking for BAM/VCF inputs.

This module provides functionality for:
- Checking compatibility between BAM/VCF files and reference FASTA
- Validating contig overlap and length thresholds per spec §7
- Auto-detecting duplicate flags in BAM files
- Running transient samtools markdup for depth analysis if needed
"""

import logging
import subprocess
import shutil
from pathlib import Path
from typing import Dict, Set, Optional
from dataclasses import dataclass
import pysam

logger = logging.getLogger(__name__)


class CompatibilityError(Exception):
    """Raised when BAM/VCF compatibility checks fail."""

    pass


@dataclass
class ContigInfo:
    """Information about a contig."""

    name: str
    length: int


@dataclass
class CompatibilityResult:
    """Result of compatibility checking."""

    is_compatible: bool
    shared_length: int
    total_fasta_length: int
    coverage_fraction: float
    per_contig_diffs: Dict[str, float]  # contig -> length diff percentage
    total_length_diff: float  # percentage difference in total length
    missing_contigs: Set[str]  # contigs in FASTA but not in BAM/VCF
    extra_contigs: Set[str]  # contigs in BAM/VCF but not in FASTA
    error_message: Optional[str] = None


def get_fasta_contigs(fasta_path: Path, min_length: int = 500) -> Dict[str, int]:
    """
    Extract contig names and lengths from FASTA file.

    Args:
        fasta_path: Path to FASTA file
        min_length: Minimum contig length to include (default 500bp as per spec)

    Returns:
        Dictionary mapping contig names to lengths

    Raises:
        CompatibilityError: If FASTA file cannot be read
    """
    try:
        contigs = {}
        with pysam.FastaFile(str(fasta_path)) as fasta:
            for contig_name in fasta.references:
                length = fasta.get_reference_length(contig_name)
                if length >= min_length:
                    contigs[contig_name] = length

        logger.debug(f"Found {len(contigs)} contigs ≥{min_length}bp in {fasta_path}")
        return contigs

    except Exception as e:
        raise CompatibilityError(f"Failed to read FASTA file {fasta_path}: {e}") from e


def get_bam_contigs(bam_path: Path) -> Dict[str, int]:
    """
    Extract contig names and lengths from BAM header.

    Args:
        bam_path: Path to BAM file

    Returns:
        Dictionary mapping contig names to lengths

    Raises:
        CompatibilityError: If BAM file cannot be read
    """
    try:
        contigs = {}
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            for contig_name, length in zip(bam.references, bam.lengths):
                contigs[contig_name] = length

        logger.debug(f"Found {len(contigs)} contigs in BAM header: {bam_path}")
        return contigs

    except Exception as e:
        raise CompatibilityError(f"Failed to read BAM file {bam_path}: {e}") from e


def get_vcf_contigs(vcf_path: Path) -> Dict[str, int]:
    """
    Extract contig names and lengths from VCF header.

    Args:
        vcf_path: Path to VCF file

    Returns:
        Dictionary mapping contig names to lengths

    Raises:
        CompatibilityError: If VCF file cannot be read
    """
    try:
        contigs = {}
        with pysam.VariantFile(str(vcf_path)) as vcf:
            for contig_name, record in vcf.header.contigs.items():
                # Access length attribute directly
                if hasattr(record, "length") and record.length is not None:
                    contigs[str(contig_name)] = int(record.length)

        logger.debug(f"Found {len(contigs)} contigs in VCF header: {vcf_path}")
        return contigs

    except Exception as e:
        raise CompatibilityError(f"Failed to read VCF file {vcf_path}: {e}") from e


def check_compatibility(
    fasta_contigs: Dict[str, int],
    other_contigs: Dict[str, int],
    coverage_threshold: float = 0.95,
    per_contig_threshold: float = 0.01,
    total_threshold: float = 0.02,
) -> CompatibilityResult:
    """
    Check compatibility between FASTA and BAM/VCF contigs.

    As per spec §7:
    - ≥95% of analyzed FASTA length covered by shared contigs (name match)
    - Per-contig length diff ≤1%
    - Total diff ≤2%

    Args:
        fasta_contigs: Contigs from FASTA file (name -> length)
        other_contigs: Contigs from BAM/VCF file (name -> length)
        coverage_threshold: Minimum fraction of FASTA length that must be covered
        per_contig_threshold: Maximum per-contig length difference (fraction)
        total_threshold: Maximum total length difference (fraction)

    Returns:
        CompatibilityResult with detailed information
    """
    # Find shared contigs
    shared_names = set(fasta_contigs.keys()) & set(other_contigs.keys())
    missing_contigs = set(fasta_contigs.keys()) - set(other_contigs.keys())
    extra_contigs = set(other_contigs.keys()) - set(fasta_contigs.keys())

    # Calculate shared length and coverage
    shared_length = sum(fasta_contigs[name] for name in shared_names)
    total_fasta_length = sum(fasta_contigs.values())
    coverage_fraction = (
        shared_length / total_fasta_length if total_fasta_length > 0 else 0.0
    )

    # Check per-contig length differences
    per_contig_diffs = {}
    per_contig_violations = []

    for contig_name in shared_names:
        fasta_len = fasta_contigs[contig_name]
        other_len = other_contigs[contig_name]

        if fasta_len > 0:
            diff_fraction = abs(fasta_len - other_len) / fasta_len
            per_contig_diffs[contig_name] = diff_fraction

            if diff_fraction > per_contig_threshold:
                per_contig_violations.append(
                    f"{contig_name}: {diff_fraction:.3f} ({fasta_len} vs {other_len})"
                )

    # Calculate total length difference
    total_fasta_shared = sum(fasta_contigs[name] for name in shared_names)
    total_other_shared = sum(other_contigs[name] for name in shared_names)
    total_length_diff = (
        abs(total_fasta_shared - total_other_shared) / total_fasta_shared
        if total_fasta_shared > 0
        else 0.0
    )

    # Determine compatibility
    violations = []

    if coverage_fraction < coverage_threshold:
        violations.append(
            f"Coverage {coverage_fraction:.3f} < {coverage_threshold:.3f} "
            f"({shared_length}/{total_fasta_length} bp)"
        )

    if per_contig_violations:
        violations.append(
            f"Per-contig length violations: {'; '.join(per_contig_violations)}"
        )

    if total_length_diff > total_threshold:
        violations.append(
            f"Total length diff {total_length_diff:.3f} > {total_threshold:.3f}"
        )

    is_compatible = len(violations) == 0
    error_message = "; ".join(violations) if violations else None

    logger.info(
        f"Compatibility check: coverage={coverage_fraction:.3f}, "
        f"per_contig_max={max(per_contig_diffs.values()) if per_contig_diffs else 0:.3f}, "
        f"total_diff={total_length_diff:.3f}, compatible={is_compatible}"
    )

    return CompatibilityResult(
        is_compatible=is_compatible,
        shared_length=shared_length,
        total_fasta_length=total_fasta_length,
        coverage_fraction=coverage_fraction,
        per_contig_diffs=per_contig_diffs,
        total_length_diff=total_length_diff,
        missing_contigs=missing_contigs,
        extra_contigs=extra_contigs,
        error_message=error_message,
    )


def check_bam_compatibility(bam_path: Path, fasta_path: Path) -> CompatibilityResult:
    """
    Check compatibility between BAM and FASTA files.

    Args:
        bam_path: Path to BAM file
        fasta_path: Path to FASTA file

    Returns:
        CompatibilityResult

    Raises:
        CompatibilityError: If files cannot be read
    """
    fasta_contigs = get_fasta_contigs(fasta_path)
    bam_contigs = get_bam_contigs(bam_path)

    return check_compatibility(fasta_contigs, bam_contigs)


def check_vcf_compatibility(vcf_path: Path, fasta_path: Path) -> CompatibilityResult:
    """
    Check compatibility between VCF and FASTA files.

    Args:
        vcf_path: Path to VCF file
        fasta_path: Path to FASTA file

    Returns:
        CompatibilityResult

    Raises:
        CompatibilityError: If files cannot be read
    """
    fasta_contigs = get_fasta_contigs(fasta_path)
    vcf_contigs = get_vcf_contigs(vcf_path)

    return check_compatibility(fasta_contigs, vcf_contigs)


def has_duplicate_flags(bam_path: Path) -> bool:
    """
    Check if BAM file has duplicate flags set.

    Samples a subset of reads to determine if duplicate marking has been performed.

    Args:
        bam_path: Path to BAM file

    Returns:
        True if duplicate flags are detected, False otherwise

    Raises:
        CompatibilityError: If BAM file cannot be read
    """
    try:
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            # Check first 10000 reads for duplicate flags
            duplicate_count = 0
            total_count = 0

            for read in bam.fetch():
                total_count += 1
                if read.is_duplicate:
                    duplicate_count += 1

                # Stop after checking enough reads
                if total_count >= 10000:
                    break

            if total_count == 0:
                logger.warning(f"No reads found in BAM file: {bam_path}")
                return False

            duplicate_fraction = duplicate_count / total_count
            has_dups = duplicate_count > 0

            logger.debug(
                f"Duplicate check: {duplicate_count}/{total_count} "
                f"({duplicate_fraction:.3f}) reads marked as duplicates"
            )

            return has_dups

    except Exception as e:
        raise CompatibilityError(
            f"Failed to check duplicates in BAM file {bam_path}: {e}"
        ) from e


def run_markdup_for_depth(bam_path: Path, output_dir: Path) -> Path:
    """
    Run samtools markdup transiently to mark duplicates for depth analysis.

    Args:
        bam_path: Path to input BAM file
        output_dir: Directory for temporary output

    Returns:
        Path to temporary BAM file with duplicates marked

    Raises:
        CompatibilityError: If samtools markdup fails
    """
    if not shutil.which("samtools"):
        raise CompatibilityError("samtools not found in PATH")

    # Create temporary output path
    temp_bam = output_dir / f"{bam_path.stem}.markdup.bam"

    try:
        # First sort by queryname for markdup
        queryname_bam = output_dir / f"{bam_path.stem}.queryname.bam"

        sort_cmd = [
            "samtools",
            "sort",
            "-n",  # Sort by queryname
            "-o",
            str(queryname_bam),
            str(bam_path),
        ]

        logger.debug(f"Sorting by queryname: {bam_path} -> {queryname_bam}")

        subprocess.run(sort_cmd, capture_output=True, text=True, check=True)

        # Then run fixmate to add MC tags
        fixmate_bam = output_dir / f"{bam_path.stem}.fixmate.bam"

        fixmate_cmd = [
            "samtools",
            "fixmate",
            "-m",  # Add mate score tags
            str(queryname_bam),
            str(fixmate_bam),
        ]

        logger.debug(f"Running samtools fixmate: {queryname_bam} -> {fixmate_bam}")

        subprocess.run(fixmate_cmd, capture_output=True, text=True, check=True)

        # Sort again by coordinate for markdup
        coord_sorted_bam = output_dir / f"{bam_path.stem}.coord_sorted.bam"

        coord_sort_cmd = [
            "samtools",
            "sort",
            "-o",
            str(coord_sorted_bam),
            str(fixmate_bam),
        ]

        logger.debug(f"Sorting by coordinate: {fixmate_bam} -> {coord_sorted_bam}")

        subprocess.run(coord_sort_cmd, capture_output=True, text=True, check=True)

        # Finally run samtools markdup on the coordinate-sorted fixmate output
        cmd = ["samtools", "markdup", str(coord_sorted_bam), str(temp_bam)]

        logger.info(f"Running samtools markdup: {coord_sorted_bam} -> {temp_bam}")

        subprocess.run(cmd, capture_output=True, text=True, check=True)

        logger.debug("samtools markdup completed successfully")

        # Clean up intermediate files
        for temp_file in [queryname_bam, fixmate_bam, coord_sorted_bam]:
            if temp_file.exists():
                temp_file.unlink()

        # Index the marked BAM file
        index_cmd = ["samtools", "index", str(temp_bam)]
        subprocess.run(index_cmd, capture_output=True, text=True, check=True)

        return temp_bam

    except subprocess.CalledProcessError as e:
        raise CompatibilityError(f"samtools markdup failed: {e.stderr}") from e
    except Exception as e:
        raise CompatibilityError(f"Failed to run markdup: {e}") from e
