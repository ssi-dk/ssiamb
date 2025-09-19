"""
Mapping and index handling for minimap2 and bwa-mem2.

This module provides functionality for:
- Building and managing indexes for reference sequences
- Mapping paired-end FASTQ files to references
- Generating sorted BAM files with proper read groups
"""

import subprocess
import logging
from pathlib import Path
from typing import Dict, List, Optional, Any
import shutil

from .models import Mapper


logger = logging.getLogger(__name__)


class MappingError(Exception):
    """Raised when mapping operations fail."""

    pass


class ExternalToolError(Exception):
    """Raised when external tools are missing or fail."""

    pass


def check_external_tools() -> Dict[str, Dict[str, Any]]:
    """
    Check availability and versions of external mapping tools.

    Returns:
        Dictionary mapping tool names to availability and version info.
        Each tool entry contains:
        - 'available': boolean status
        - 'version': version string if available, or error message

    Examples:
        >>> tools = check_external_tools()
        >>> if not tools['minimap2']['available']:
        ...     raise ExternalToolError("minimap2 not found")
        >>> print(f"minimap2 version: {tools['minimap2']['version']}")
    """
    tools = {}

    # Define tools and their version commands
    tool_commands = {
        "minimap2": ["minimap2", "--version"],
        "bwa-mem2": ["bwa-mem2", "version"],
        "samtools": ["samtools", "--version"],
    }

    for tool, version_cmd in tool_commands.items():
        tools[tool] = {"available": False, "version": "unknown"}

        # Check if tool is in PATH
        if shutil.which(tool) is None:
            tools[tool]["version"] = "not found in PATH"
            logger.debug(f"Tool {tool}: not found")
            continue

        # Try to get version
        try:
            result = subprocess.run(
                version_cmd,
                capture_output=True,
                text=True,
                timeout=5,  # Prevent hanging
            )
            if result.returncode == 0:
                tools[tool]["available"] = True
                # Extract version from output (first line typically)
                version_line = (
                    result.stdout.strip().split("\n")[0]
                    if result.stdout
                    else "version check succeeded"
                )
                tools[tool]["version"] = version_line
                logger.debug(f"Tool {tool}: available, version: {version_line}")
            else:
                tools[tool][
                    "version"
                ] = f"version check failed (exit {result.returncode})"
                logger.debug(f"Tool {tool}: found but version check failed")
        except subprocess.TimeoutExpired:
            tools[tool]["version"] = "version check timed out"
            logger.debug(f"Tool {tool}: found but version check timed out")
        except Exception as e:
            tools[tool]["version"] = f"version check error: {e}"
            logger.debug(f"Tool {tool}: found but version check error: {e}")

    return tools


def get_index_files(fasta_path: Path, mapper: Mapper) -> List[Path]:
    """
    Get expected index file paths for a given FASTA and mapper.

    Args:
        fasta_path: Path to reference FASTA file
        mapper: Mapper type (minimap2 or bwa-mem2)

    Returns:
        List of expected index file paths

    Raises:
        ValueError: If mapper is not supported
    """
    if mapper == Mapper.MINIMAP2:
        return [fasta_path.with_suffix(".mmi")]

    elif mapper == Mapper.BWA_MEM2:
        # BWA-MEM2 creates index files with the original filename plus extensions
        extensions = [".0123", ".amb", ".ann", ".pac", ".bwt.2bit.64"]
        return [Path(str(fasta_path) + ext) for ext in extensions]

    else:
        raise ValueError(f"Unsupported mapper: {mapper}")


def indexes_exist(fasta_path: Path, mapper: Mapper) -> bool:
    """
    Check if all required index files exist for a given FASTA and mapper.

    Args:
        fasta_path: Path to reference FASTA file
        mapper: Mapper type

    Returns:
        True if all index files exist, False otherwise
    """
    index_files = get_index_files(fasta_path, mapper)
    return all(idx_file.exists() for idx_file in index_files)


def index_bam(bam_path: Path) -> None:
    """
    Index a BAM file using samtools.

    Args:
        bam_path: Path to BAM file to index

    Raises:
        MappingError: If indexing fails
    """
    if not bam_path.exists():
        raise MappingError(f"BAM file not found: {bam_path}")

    index_path = bam_path.with_suffix(".bam.bai")

    logger.debug(f"Indexing BAM file: {bam_path}")

    try:
        cmd = ["samtools", "index", str(bam_path)]
        subprocess.run(cmd, capture_output=True, text=True, check=True)

        if not index_path.exists():
            raise MappingError(f"BAM index was not created: {index_path}")

        logger.debug(f"Created BAM index: {index_path}")

    except subprocess.CalledProcessError as e:
        raise MappingError(f"Failed to index BAM {bam_path}: {e.stderr}")


def ensure_indexes_self(fasta_path: Path, mapper: Mapper) -> None:
    """
    Ensure index files exist for self-mode mapping, building them if missing.

    This function builds index files next to the FASTA file if they don't exist.
    For minimap2, creates .mmi file. For bwa-mem2, creates .0123, .amb, .ann,
    .pac, .bwt.2bit.64 files.

    Args:
        fasta_path: Path to reference FASTA file
        mapper: Mapper type (minimap2 or bwa-mem2)

    Raises:
        FileNotFoundError: If FASTA file doesn't exist
        ExternalToolError: If required mapper tool is not available
        MappingError: If index building fails

    Examples:
        >>> ensure_indexes_self(Path("ref.fasta"), Mapper.MINIMAP2)
        # Creates ref.mmi if it doesn't exist
    """
    if not fasta_path.exists():
        raise FileNotFoundError(f"Reference FASTA not found: {fasta_path}")

    # Check if indexes already exist
    if indexes_exist(fasta_path, mapper):
        logger.info(f"Index files already exist for {fasta_path} ({mapper.value})")
        return

    # Check tool availability
    tools = check_external_tools()
    tool_name = mapper.value  # Use the actual tool name without modification
    if not tools.get(tool_name, {}).get("available", False):
        version_info = tools.get(tool_name, {}).get("version", "unknown")
        raise ExternalToolError(f"{tool_name} not available: {version_info}")

    logger.info(f"Building {mapper.value} index for {fasta_path}")

    try:
        if mapper == Mapper.MINIMAP2:
            _build_minimap2_index(fasta_path)
        elif mapper == Mapper.BWA_MEM2:
            _build_bwa_mem2_index(fasta_path)
        else:
            raise ValueError(f"Unsupported mapper: {mapper}")

    except subprocess.CalledProcessError as e:
        raise MappingError(f"Failed to build {mapper.value} index: {e}")


def _build_minimap2_index(fasta_path: Path) -> None:
    """Build minimap2 index (.mmi file)."""
    index_path = fasta_path.with_suffix(".mmi")

    cmd = [
        "minimap2",
        "-d",
        str(index_path),  # Output index file
        str(fasta_path),  # Input FASTA
    ]

    logger.debug(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True, capture_output=True, text=True)

    if not index_path.exists():
        raise MappingError(f"minimap2 index was not created: {index_path}")

    logger.info(f"Created minimap2 index: {index_path}")


def _build_bwa_mem2_index(fasta_path: Path) -> None:
    """Build bwa-mem2 index files."""
    cmd = ["bwa-mem2", "index", str(fasta_path)]

    logger.debug(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True, capture_output=True, text=True)

    # Verify all expected index files were created
    index_files = get_index_files(fasta_path, Mapper.BWA_MEM2)
    missing_files = [f for f in index_files if not f.exists()]

    if missing_files:
        raise MappingError(f"bwa-mem2 index files were not created: {missing_files}")

    logger.info(f"Created bwa-mem2 index files for: {fasta_path}")


def map_fastqs(
    mapper: Mapper,
    fasta_path: Path,
    r1_path: Path,
    r2_path: Path,
    sample_name: str,
    threads: int = 4,
    output_path: Optional[Path] = None,
) -> Path:
    """
    Map paired-end FASTQ files to reference and return sorted BAM.

    Args:
        mapper: Mapper type (minimap2 or bwa-mem2)
        fasta_path: Path to reference FASTA file
        r1_path: Path to R1 FASTQ file
        r2_path: Path to R2 FASTQ file
        sample_name: Sample name for read group
        threads: Number of threads to use
        output_path: Output BAM path (auto-generated if None)

    Returns:
        Path to sorted BAM file

    Raises:
        FileNotFoundError: If input files don't exist
        ExternalToolError: If required tools are not available
        MappingError: If mapping fails

    Examples:
        >>> bam_path = map_fastqs(
        ...     Mapper.MINIMAP2,
        ...     Path("ref.fasta"),
        ...     Path("sample_R1.fastq.gz"),
        ...     Path("sample_R2.fastq.gz"),
        ...     "sample123"
        ... )
    """
    # Validate input files
    for file_path in [fasta_path, r1_path, r2_path]:
        if not file_path.exists():
            raise FileNotFoundError(f"Input file not found: {file_path}")

    # Check tool availability
    tools = check_external_tools()
    tool_name = mapper.value  # Use the actual tool name without modification
    if not tools.get(tool_name, {}).get("available", False) or not tools.get(
        "samtools", {}
    ).get("available", False):
        missing = [
            t
            for t, info in tools.items()
            if not info.get("available", False) and t in [tool_name, "samtools"]
        ]
        raise ExternalToolError(f"Required tools not available: {missing}")

    # Ensure indexes exist
    ensure_indexes_self(fasta_path, mapper)

    # Generate output path if not provided
    if output_path is None:
        output_path = Path(f"{sample_name}.sorted.bam")

    logger.info(
        f"Mapping {r1_path.name} and {r2_path.name} to {fasta_path.name} using {mapper.value}"
    )

    try:
        if mapper == Mapper.MINIMAP2:
            _map_with_minimap2(
                fasta_path, r1_path, r2_path, sample_name, threads, output_path
            )
        elif mapper == Mapper.BWA_MEM2:
            _map_with_bwa_mem2(
                fasta_path, r1_path, r2_path, sample_name, threads, output_path
            )
        else:
            raise ValueError(f"Unsupported mapper: {mapper}")

    except subprocess.CalledProcessError as e:
        raise MappingError(f"Mapping with {mapper.value} failed: {e}")

    if not output_path.exists():
        raise MappingError(f"Output BAM was not created: {output_path}")

    # Index the BAM file
    index_bam(output_path)

    logger.info(f"Created sorted BAM: {output_path}")
    return output_path


def _map_with_minimap2(
    fasta_path: Path,
    r1_path: Path,
    r2_path: Path,
    sample_name: str,
    threads: int,
    output_path: Path,
) -> None:
    """Map with minimap2 and create sorted BAM."""
    index_path = fasta_path.with_suffix(".mmi")

    # Create read group string
    read_group = f"@RG\\tID:{sample_name}\\tSM:{sample_name}\\tPL:ILLUMINA"

    # minimap2 command: map to SAM, pipe to samtools for sorting
    minimap2_cmd = [
        "minimap2",
        "-x",
        "sr",  # Short read preset
        "-a",  # Output alignments (required by spec)
        "-t",
        str(threads),  # Threads
        "-R",
        read_group,  # Read group
        str(index_path),  # Index file
        str(r1_path),  # R1 FASTQ
        str(r2_path),  # R2 FASTQ
    ]

    samtools_cmd = [
        "samtools",
        "sort",
        "-@",
        str(threads),  # Threads
        "-o",
        str(output_path),  # Output file
    ]

    logger.debug(f"Running: {' '.join(minimap2_cmd)} | {' '.join(samtools_cmd)}")

    # Run minimap2 and pipe to samtools sort
    minimap2_proc = subprocess.Popen(
        minimap2_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )

    samtools_proc = subprocess.Popen(
        samtools_cmd, stdin=minimap2_proc.stdout, stderr=subprocess.PIPE, text=True
    )

    # Close minimap2 stdout so samtools gets EOF
    if minimap2_proc.stdout:
        minimap2_proc.stdout.close()

    # Wait for both processes
    samtools_stderr = samtools_proc.communicate()[1]
    minimap2_stderr = minimap2_proc.communicate()[1]

    # Check return codes
    if minimap2_proc.returncode != 0:
        raise subprocess.CalledProcessError(
            minimap2_proc.returncode, minimap2_cmd, stderr=minimap2_stderr
        )

    if samtools_proc.returncode != 0:
        raise subprocess.CalledProcessError(
            samtools_proc.returncode, samtools_cmd, stderr=samtools_stderr
        )


def _map_with_bwa_mem2(
    fasta_path: Path,
    r1_path: Path,
    r2_path: Path,
    sample_name: str,
    threads: int,
    output_path: Path,
) -> None:
    """Map with bwa-mem2 and create sorted BAM."""
    # Create read group string
    read_group = f"@RG\\tID:{sample_name}\\tSM:{sample_name}\\tPL:ILLUMINA"

    # bwa-mem2 command: map to SAM, pipe to samtools for sorting
    bwa_cmd = [
        "bwa-mem2",
        "mem",
        "-t",
        str(threads),  # Threads
        "-R",
        read_group,  # Read group
        str(fasta_path),  # Reference FASTA
        str(r1_path),  # R1 FASTQ
        str(r2_path),  # R2 FASTQ
    ]

    samtools_cmd = [
        "samtools",
        "sort",
        "-@",
        str(threads),  # Threads
        "-o",
        str(output_path),  # Output file
    ]

    logger.debug(f"Running: {' '.join(bwa_cmd)} | {' '.join(samtools_cmd)}")

    # Run bwa-mem2 and pipe to samtools sort
    bwa_proc = subprocess.Popen(
        bwa_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )

    samtools_proc = subprocess.Popen(
        samtools_cmd, stdin=bwa_proc.stdout, stderr=subprocess.PIPE, text=True
    )

    # Close bwa stdout so samtools gets EOF
    if bwa_proc.stdout:
        bwa_proc.stdout.close()

    # Wait for both processes
    samtools_stderr = samtools_proc.communicate()[1]
    bwa_stderr = bwa_proc.communicate()[1]

    # Check return codes
    if bwa_proc.returncode != 0:
        raise subprocess.CalledProcessError(
            bwa_proc.returncode, bwa_cmd, stderr=bwa_stderr
        )

    if samtools_proc.returncode != 0:
        raise subprocess.CalledProcessError(
            samtools_proc.returncode, samtools_cmd, stderr=samtools_stderr
        )


def calculate_mapping_rate(bam_path: Path) -> float:
    """
    Calculate mapping rate from BAM file using samtools stats.

    Args:
        bam_path: Path to sorted BAM file

    Returns:
        Mapping rate as fraction (0.0-1.0)

    Raises:
        MappingError: If BAM file doesn't exist or samtools fails
        ExternalToolError: If samtools is not available
    """
    if not shutil.which("samtools"):
        raise ExternalToolError("samtools not found in PATH")

    if not bam_path.exists():
        raise MappingError(f"BAM file not found: {bam_path}")

    logger.debug(f"Calculating mapping rate for {bam_path}")

    try:
        # Run samtools stats to get read counts
        cmd = ["samtools", "stats", str(bam_path)]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)

        total_reads = 0
        mapped_reads = 0

        # Parse samtools stats output
        for line in result.stdout.splitlines():
            if line.startswith("SN\t"):
                parts = line.split("\t")
                if len(parts) >= 3:
                    metric = parts[1]
                    value_str = parts[2].strip()

                    # Extract numbers (some lines have additional text after the number)
                    try:
                        value = int(value_str.split()[0])
                    except (ValueError, IndexError):
                        continue

                    if "raw total sequences:" in metric:
                        total_reads = value
                    elif "reads mapped:" in metric:
                        mapped_reads = value

        if total_reads == 0:
            logger.warning(f"No reads found in BAM file {bam_path}")
            return 0.0

        mapping_rate = mapped_reads / total_reads
        logger.debug(f"Mapping rate: {mapped_reads}/{total_reads} = {mapping_rate:.4f}")

        return mapping_rate

    except subprocess.CalledProcessError as e:
        raise MappingError(f"samtools stats failed: {e.stderr}")
    except Exception as e:
        raise MappingError(f"Failed to calculate mapping rate: {e}")
