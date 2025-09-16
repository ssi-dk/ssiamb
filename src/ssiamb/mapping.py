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
from typing import Dict, List, Optional
import shutil

from .models import Mapper


logger = logging.getLogger(__name__)


class MappingError(Exception):
    """Raised when mapping operations fail."""
    pass


class ExternalToolError(Exception):
    """Raised when external tools are missing or fail."""
    pass


def check_external_tools() -> Dict[str, bool]:
    """
    Check availability of external mapping tools.
    
    Returns:
        Dictionary mapping tool names to availability status
        
    Examples:
        >>> tools = check_external_tools()
        >>> if not tools['minimap2']:
        ...     raise ExternalToolError("minimap2 not found")
    """
    tools = {}
    
    for tool in ['minimap2', 'bwa-mem2', 'samtools']:
        tools[tool] = shutil.which(tool) is not None
        logger.debug(f"Tool {tool}: {'available' if tools[tool] else 'missing'}")
    
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
    base_path = fasta_path.with_suffix('')  # Remove .fasta/.fa/.fna extension
    
    if mapper == Mapper.MINIMAP2:
        return [fasta_path.with_suffix('.mmi')]
    
    elif mapper == Mapper.BWA_MEM2:
        extensions = ['.0123', '.amb', '.ann', '.pac', '.bwt.2bit.64']
        return [base_path.with_suffix(ext) for ext in extensions]
    
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
    tool_name = mapper.value.replace('-', '')  # minimap2 or bwamem2
    if not tools.get(tool_name):
        raise ExternalToolError(f"{tool_name} not found in PATH")
    
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
    index_path = fasta_path.with_suffix('.mmi')
    
    cmd = [
        'minimap2',
        '-d', str(index_path),  # Output index file
        str(fasta_path)         # Input FASTA
    ]
    
    logger.debug(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    
    if not index_path.exists():
        raise MappingError(f"minimap2 index was not created: {index_path}")
    
    logger.info(f"Created minimap2 index: {index_path}")


def _build_bwa_mem2_index(fasta_path: Path) -> None:
    """Build bwa-mem2 index files."""
    cmd = [
        'bwa-mem2',
        'index',
        str(fasta_path)
    ]
    
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
    output_path: Optional[Path] = None
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
    tool_name = mapper.value.replace('-', '')
    if not tools.get(tool_name) or not tools.get('samtools'):
        missing = [t for t, avail in tools.items() if not avail and t in [tool_name, 'samtools']]
        raise ExternalToolError(f"Required tools not found: {missing}")
    
    # Ensure indexes exist
    ensure_indexes_self(fasta_path, mapper)
    
    # Generate output path if not provided
    if output_path is None:
        output_path = Path(f"{sample_name}.sorted.bam")
    
    logger.info(f"Mapping {r1_path.name} and {r2_path.name} to {fasta_path.name} using {mapper.value}")
    
    try:
        if mapper == Mapper.MINIMAP2:
            _map_with_minimap2(fasta_path, r1_path, r2_path, sample_name, threads, output_path)
        elif mapper == Mapper.BWA_MEM2:
            _map_with_bwa_mem2(fasta_path, r1_path, r2_path, sample_name, threads, output_path)
        else:
            raise ValueError(f"Unsupported mapper: {mapper}")
            
    except subprocess.CalledProcessError as e:
        raise MappingError(f"Mapping with {mapper.value} failed: {e}")
    
    if not output_path.exists():
        raise MappingError(f"Output BAM was not created: {output_path}")
    
    logger.info(f"Created sorted BAM: {output_path}")
    return output_path


def _map_with_minimap2(
    fasta_path: Path,
    r1_path: Path,
    r2_path: Path,
    sample_name: str,
    threads: int,
    output_path: Path
) -> None:
    """Map with minimap2 and create sorted BAM."""
    index_path = fasta_path.with_suffix('.mmi')
    
    # Create read group string
    read_group = f"@RG\\tID:{sample_name}\\tSM:{sample_name}\\tPL:ILLUMINA"
    
    # minimap2 command: map to SAM, pipe to samtools for sorting
    minimap2_cmd = [
        'minimap2',
        '-ax', 'sr',           # Short read preset
        '-t', str(threads),    # Threads
        '-R', read_group,      # Read group
        str(index_path),       # Index file
        str(r1_path),          # R1 FASTQ
        str(r2_path)           # R2 FASTQ
    ]
    
    samtools_cmd = [
        'samtools', 'sort',
        '-@', str(threads),    # Threads
        '-o', str(output_path) # Output file
    ]
    
    logger.debug(f"Running: {' '.join(minimap2_cmd)} | {' '.join(samtools_cmd)}")
    
    # Run minimap2 and pipe to samtools sort
    minimap2_proc = subprocess.Popen(
        minimap2_cmd, 
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    
    samtools_proc = subprocess.Popen(
        samtools_cmd,
        stdin=minimap2_proc.stdout,
        stderr=subprocess.PIPE,
        text=True
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
            minimap2_proc.returncode, 
            minimap2_cmd,
            stderr=minimap2_stderr
        )
    
    if samtools_proc.returncode != 0:
        raise subprocess.CalledProcessError(
            samtools_proc.returncode,
            samtools_cmd, 
            stderr=samtools_stderr
        )


def _map_with_bwa_mem2(
    fasta_path: Path,
    r1_path: Path,
    r2_path: Path,
    sample_name: str,
    threads: int,
    output_path: Path
) -> None:
    """Map with bwa-mem2 and create sorted BAM."""
    # Create read group string
    read_group = f"@RG\\tID:{sample_name}\\tSM:{sample_name}\\tPL:ILLUMINA"
    
    # bwa-mem2 command: map to SAM, pipe to samtools for sorting
    bwa_cmd = [
        'bwa-mem2', 'mem',
        '-t', str(threads),    # Threads
        '-R', read_group,      # Read group
        str(fasta_path),       # Reference FASTA
        str(r1_path),          # R1 FASTQ
        str(r2_path)           # R2 FASTQ
    ]
    
    samtools_cmd = [
        'samtools', 'sort',
        '-@', str(threads),    # Threads
        '-o', str(output_path) # Output file
    ]
    
    logger.debug(f"Running: {' '.join(bwa_cmd)} | {' '.join(samtools_cmd)}")
    
    # Run bwa-mem2 and pipe to samtools sort
    bwa_proc = subprocess.Popen(
        bwa_cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    
    samtools_proc = subprocess.Popen(
        samtools_cmd,
        stdin=bwa_proc.stdout,
        stderr=subprocess.PIPE,
        text=True
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
            bwa_proc.returncode,
            bwa_cmd,
            stderr=bwa_stderr
        )
    
    if samtools_proc.returncode != 0:
        raise subprocess.CalledProcessError(
            samtools_proc.returncode,
            samtools_cmd,
            stderr=samtools_stderr
        )