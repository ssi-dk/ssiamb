"""
Core data models for ssiamb.

This module defines the primary data structures used throughout the application,
including configuration objects, output records, and enums for various options.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Optional, Dict, Any


class Mode(Enum):
    """Analysis mode."""
    SELF = "self"
    REF = "ref"
    SUMMARIZE = "summarize"


class Mapper(Enum):
    """Read mapping tool."""
    MINIMAP2 = "minimap2"
    BWA_MEM2 = "bwa-mem2"


class Caller(Enum):
    """Variant calling tool."""
    BBTOOLS = "bbtools"
    BCFTOOLS = "bcftools"


class DepthTool(Enum):
    """Depth analysis tool."""
    MOSDEPTH = "mosdepth"
    SAMTOOLS = "samtools"


class TSVMode(Enum):
    """TSV output mode."""
    OVERWRITE = "overwrite"
    APPEND = "append"
    FAIL = "fail"


class OnFail(Enum):
    """Action when reference resolution fails."""
    ERROR = "error"
    WARN = "warn"
    SKIP = "skip"


@dataclass
class Thresholds:
    """Thresholds for ambiguous site detection."""
    dp_min: int = 10
    maf_min: float = 0.1
    dp_cap: int = 100
    mapq_min: int = 20
    baseq_min: int = 20
    
    def __post_init__(self) -> None:
        """Validate threshold values."""
        if self.dp_min < 1:
            raise ValueError("dp_min must be >= 1")
        if not 0.0 <= self.maf_min <= 1.0:
            raise ValueError("maf_min must be between 0.0 and 1.0")
        if self.dp_cap < self.dp_min:
            raise ValueError("dp_cap must be >= dp_min")
        if self.mapq_min < 0:
            raise ValueError("mapq_min must be >= 0")
        if self.baseq_min < 0:
            raise ValueError("baseq_min must be >= 0")


@dataclass
class SummaryRow:
    """
    Single row of the ambiguous_summary.tsv output.
    
    Based on spec.md §4.1 - Primary Output Format.
    """
    sample: str
    mode: str
    ref_label: str
    mapper: str
    caller: str
    dp_min: int
    maf_min: float
    
    # Denominator metrics
    callable_bases: int
    genome_length: int
    breadth_10x: float
    mean_depth: float
    mapping_rate: float
    
    # Numerator counts
    ambiguous_sites: int
    ref_calls: int
    alt_calls: int
    no_call: int
    
    # Optional QC flags (comma-separated warnings)
    qc_warnings: str = ""
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for TSV output."""
        return {
            "sample": self.sample,
            "mode": self.mode,
            "ref_label": self.ref_label,
            "mapper": self.mapper,
            "caller": self.caller,
            "dp_min": self.dp_min,
            "maf_min": self.maf_min,
            "callable_bases": self.callable_bases,
            "genome_length": self.genome_length,
            "breadth_10x": self.breadth_10x,
            "mean_depth": self.mean_depth,
            "mapping_rate": self.mapping_rate,
            "ambiguous_sites": self.ambiguous_sites,
            "ref_calls": self.ref_calls,
            "alt_calls": self.alt_calls,
            "no_call": self.no_call,
            "qc_warnings": self.qc_warnings,
        }


@dataclass
class Paths:
    """File paths for input and output."""
    r1: Path
    r2: Path
    assembly: Optional[Path] = None
    reference: Optional[Path] = None
    output_dir: Path = Path(".")
    sample: Optional[str] = None
    
    # Intermediate file paths (computed)
    bam: Optional[Path] = None
    vcf: Optional[Path] = None
    depth_summary: Optional[Path] = None
    
    def __post_init__(self) -> None:
        """Validate paths and resolve output directory."""
        if not self.r1.exists():
            raise FileNotFoundError(f"R1 file not found: {self.r1}")
        if not self.r2.exists():
            raise FileNotFoundError(f"R2 file not found: {self.r2}")
        
        # Ensure output directory exists
        self.output_dir.mkdir(parents=True, exist_ok=True)


@dataclass
class RunPlan:
    """
    Complete execution plan for a ssiamb run.
    
    This captures all resolved inputs, tools, and output specifications
    before execution begins.
    """
    mode: Mode
    sample: str
    paths: Paths
    thresholds: Thresholds
    
    # Tool selection
    mapper: Mapper = Mapper.MINIMAP2
    caller: Caller = Caller.BBTOOLS
    depth_tool: DepthTool = DepthTool.MOSDEPTH
    
    # Reference information (for ref mode)
    ref_source: str = "unknown"  # "file", "species", "bracken"
    ref_label: str = "unknown"
    
    # Output options
    emit_vcf: bool = False
    emit_bed: bool = False
    emit_matrix: bool = False
    emit_per_contig: bool = False
    emit_provenance: bool = False
    emit_multiqc: bool = False
    to_stdout: bool = False
    
    # Runtime options
    threads: int = 1
    dry_run: bool = False
    
    def get_sample_prefix(self) -> str:
        """Get prefix for output files."""
        return f"{self.sample}_{self.mode.value}"


@dataclass
class Provenance:
    """
    Provenance information for reproducibility.
    
    Tracks tool versions, command lines, input file hashes, etc.
    """
    ssiamb_version: str
    command_line: str
    timestamp: str
    
    # Input file information
    input_files: Dict[str, str] = field(default_factory=dict)  # filename -> md5
    tool_versions: Dict[str, str] = field(default_factory=dict)  # tool -> version
    
    # Runtime environment
    hostname: str = ""
    working_dir: str = ""
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON output."""
        return {
            "ssiamb_version": self.ssiamb_version,
            "command_line": self.command_line,
            "timestamp": self.timestamp,
            "input_files": self.input_files,
            "tool_versions": self.tool_versions,
            "hostname": self.hostname,
            "working_dir": self.working_dir,
        }