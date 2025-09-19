"""
Depth analysis using mosdepth for computing denominator metrics.

This module provides functionality for:
- Running mosdepth on BAM files to compute depth statistics
- Parsing mosdepth summary output files
- Computing callable bases, genome length, breadth, and mean depth
- Filtering contigs by size (≥500 bp) for consistent numerator/denominator
"""

import subprocess
import logging
from pathlib import Path
from typing import List, Set
from dataclasses import dataclass

logger = logging.getLogger(__name__)


class DepthAnalysisError(Exception):
    """Raised when depth analysis operations fail."""

    pass


@dataclass
class ContigDepthStats:
    """Depth statistics for a single contig."""

    name: str
    length: int
    bases_covered: int
    mean_depth: float
    breadth_10x: float

    @property
    def is_long_enough(self) -> bool:
        """Check if contig meets minimum length threshold (≥500 bp)."""
        return self.length >= 500


@dataclass
class DepthSummary:
    """Summary of depth analysis results."""

    callable_bases: int  # Total bases with depth ≥10 in contigs ≥500bp
    genome_length: int  # Total length of contigs ≥500bp
    breadth_10x: float  # Fraction of genome with depth ≥10x
    mean_depth: float  # Mean depth across all callable bases
    total_contigs: int  # Total number of contigs
    included_contigs: int  # Number of contigs ≥500bp
    contig_stats: List[ContigDepthStats]  # Per-contig statistics

    @property
    def included_contig_names(self) -> Set[str]:
        """Set of contig names that meet length threshold."""
        return {contig.name for contig in self.contig_stats if contig.is_long_enough}


def check_mosdepth_available() -> bool:
    """
    Check if mosdepth is available in PATH.

    Returns:
        True if mosdepth is available, False otherwise
    """
    try:
        subprocess.run(["mosdepth", "--version"], capture_output=True, check=True)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False


def run_mosdepth(
    bam_path: Path, output_prefix: Path, mapq_threshold: int = 30, threads: int = 4
) -> Path:
    """
    Run mosdepth on a BAM file to compute depth statistics.

    Args:
        bam_path: Path to input BAM file
        output_prefix: Output prefix for mosdepth files
        mapq_threshold: Minimum mapping quality (default: 30)
        threads: Number of threads to use

    Returns:
        Path to the generated summary file (.mosdepth.summary.txt)

    Raises:
        DepthAnalysisError: If mosdepth is not available or execution fails
        FileNotFoundError: If BAM file doesn't exist
    """
    if not bam_path.exists():
        raise FileNotFoundError(f"BAM file not found: {bam_path}")

    if not check_mosdepth_available():
        raise DepthAnalysisError("mosdepth not found in PATH")

    # Construct mosdepth command
    cmd = [
        "mosdepth",
        "--mapq",
        str(mapq_threshold),
        "--no-per-base",  # Skip per-base output for efficiency
        "--threads",
        str(threads),
        str(output_prefix),
        str(bam_path),
    ]

    logger.info(f"Running mosdepth: {' '.join(cmd)}")

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)

        # Log any stderr output for debugging
        if result.stderr.strip():
            logger.debug(f"mosdepth stderr: {result.stderr.strip()}")

    except subprocess.CalledProcessError as e:
        error_msg = f"mosdepth failed with return code {e.returncode}"
        if e.stderr:
            error_msg += f": {e.stderr.strip()}"
        raise DepthAnalysisError(error_msg)

    # Check that expected output file was created
    summary_file = Path(str(output_prefix) + ".mosdepth.summary.txt")
    if not summary_file.exists():
        raise DepthAnalysisError(
            f"Expected mosdepth summary file not created: {summary_file}"
        )

    logger.info(f"mosdepth completed successfully: {summary_file}")
    return summary_file


def parse_mosdepth_summary(summary_file: Path) -> DepthSummary:
    """
    Parse mosdepth summary file to extract depth statistics.

    The mosdepth summary file has format:
    chrom	length	bases	mean	min	max
    contig1	5000	4850	25.3	0	100
    ...
    total	50000	48500	24.8	0	100

    Args:
        summary_file: Path to .mosdepth.summary.txt file

    Returns:
        DepthSummary object with computed statistics

    Raises:
        DepthAnalysisError: If file parsing fails
        FileNotFoundError: If summary file doesn't exist
    """
    if not summary_file.exists():
        raise FileNotFoundError(f"mosdepth summary file not found: {summary_file}")

    logger.info(f"Parsing mosdepth summary: {summary_file}")

    try:
        contig_stats = []

        with open(summary_file, "r") as f:
            # Skip header line
            header = f.readline().strip()
            if not header.startswith("chrom"):
                raise DepthAnalysisError(f"Unexpected header format: {header}")

            for line_num, line in enumerate(f, start=2):
                line = line.strip()
                if not line:
                    continue

                parts = line.split("\t")
                if len(parts) < 4:
                    raise DepthAnalysisError(
                        f"Invalid line format at line {line_num}: {line}"
                    )

                chrom, length_str, bases_str, mean_str = parts[:4]

                try:
                    length = int(length_str)
                    bases_covered = int(bases_str)
                    mean_depth = float(mean_str)
                except ValueError as e:
                    raise DepthAnalysisError(
                        f"Invalid numeric value at line {line_num}: {e}"
                    )

                if chrom == "total":
                    pass
                else:
                    # Calculate breadth for this contig
                    breadth_10x = bases_covered / length if length > 0 else 0.0

                    contig_stats.append(
                        ContigDepthStats(
                            name=chrom,
                            length=length,
                            bases_covered=bases_covered,
                            mean_depth=mean_depth,
                            breadth_10x=breadth_10x,
                        )
                    )

        if not contig_stats:
            raise DepthAnalysisError("No contig data found in summary file")

        # Compute summary statistics for contigs ≥500bp
        long_contigs = [c for c in contig_stats if c.is_long_enough]

        if not long_contigs:
            logger.warning("No contigs ≥500bp found - all contigs are too short")
            genome_length = 0
            callable_bases = 0
            breadth_10x = 0.0
            mean_depth = 0.0
        else:
            genome_length = sum(c.length for c in long_contigs)
            callable_bases = sum(c.bases_covered for c in long_contigs)
            breadth_10x = callable_bases / genome_length if genome_length > 0 else 0.0

            # Weighted mean depth across long contigs
            total_depth = sum(c.mean_depth * c.length for c in long_contigs)
            mean_depth = total_depth / genome_length if genome_length > 0 else 0.0

        summary = DepthSummary(
            callable_bases=callable_bases,
            genome_length=genome_length,
            breadth_10x=breadth_10x,
            mean_depth=mean_depth,
            total_contigs=len(contig_stats),
            included_contigs=len(long_contigs),
            contig_stats=contig_stats,
        )

        logger.info(
            f"Parsed depth summary: {summary.included_contigs}/{summary.total_contigs} contigs ≥500bp, "
            f"{summary.callable_bases:,} callable bases, {summary.mean_depth:.1f}x mean depth"
        )

        return summary

    except Exception as e:
        if isinstance(e, DepthAnalysisError):
            raise
        raise DepthAnalysisError(f"Failed to parse mosdepth summary: {e}")


def analyze_depth(
    bam_path: Path,
    output_dir: Path,
    sample_name: str,
    mapq_threshold: int = 30,
    depth_threshold: int = 10,
    threads: int = 4,
) -> DepthSummary:
    """
    Complete depth analysis workflow: run mosdepth and parse results.

    Args:
        bam_path: Path to input BAM file
        output_dir: Directory for output files
        sample_name: Sample name for output prefix
        mapq_threshold: Minimum mapping quality (default: 30)
        depth_threshold: Depth threshold for breadth calculation (default: 10)
        threads: Number of threads to use

    Returns:
        DepthSummary object with computed statistics

    Raises:
        DepthAnalysisError: If analysis fails at any step
    """
    logger.info(f"Starting depth analysis for {sample_name}")

    # Create output directory if needed
    output_dir.mkdir(parents=True, exist_ok=True)

    # Run mosdepth
    output_prefix = output_dir / f"{sample_name}.depth"
    summary_file = run_mosdepth(
        bam_path=bam_path,
        output_prefix=output_prefix,
        mapq_threshold=mapq_threshold,
        threads=threads,
    )

    # Parse results
    summary = parse_mosdepth_summary(summary_file)

    logger.info(
        f"Depth analysis complete for {sample_name}: "
        f"{summary.genome_length:,} bp genome, "
        f"{summary.breadth_10x:.2%} breadth ≥{depth_threshold}x"
    )

    return summary


def get_depth_from_existing_summary(summary_file: Path) -> DepthSummary:
    """
    Parse existing mosdepth summary file without running mosdepth.

    Useful for reusing existing depth analysis results.

    Args:
        summary_file: Path to existing .mosdepth.summary.txt file

    Returns:
        DepthSummary object with computed statistics
    """
    return parse_mosdepth_summary(summary_file)


def list_included_contigs(summary_file: Path, min_len: int = 500) -> Set[str]:
    """
    Get set of contig names that meet minimum length threshold.

    Args:
        summary_file: Path to .mosdepth.summary.txt file
        min_len: Minimum contig length threshold

    Returns:
        Set of contig names that meet the length threshold
    """
    if not summary_file.exists():
        logger.warning(f"Summary file not found: {summary_file}")
        return set()

    try:
        included_contigs = set()

        with open(summary_file, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue

                fields = line.split("\t")
                if len(fields) < 4:
                    continue

                contig_name = fields[0]
                try:
                    contig_length = int(fields[1])
                    if contig_length >= min_len:
                        included_contigs.add(contig_name)
                except (ValueError, IndexError):
                    continue

        logger.debug(f"Found {len(included_contigs)} contigs >= {min_len} bp")
        return included_contigs

    except Exception as e:
        logger.error(f"Error reading mosdepth summary {summary_file}: {e}")
        return set()
