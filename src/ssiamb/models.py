"""
Core data models for ssiamb.

This module defines the primary data structures used throughout the application,
including configuration objects, output records, and enums for various options.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Optional, Dict, Any, TYPE_CHECKING

if TYPE_CHECKING:
    from .config import SsiambConfig


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

    dp_min: Optional[int] = None
    maf_min: Optional[float] = None
    dp_cap: Optional[int] = None
    mapq_min: Optional[int] = None
    baseq_min: Optional[int] = None

    def __post_init__(self) -> None:
        """Load defaults from configuration and validate threshold values."""
        # Import here to avoid circular imports
        from .config import get_config

        # Load defaults from config if not explicitly set
        config = get_config()
        if self.dp_min is None:
            self.dp_min = config.get_threshold("dp_min", 10)
        if self.maf_min is None:
            self.maf_min = config.get_threshold("maf_min", 0.1)
        if self.dp_cap is None:
            self.dp_cap = config.get_threshold("dp_cap", 100)
        if self.mapq_min is None:
            self.mapq_min = config.get_threshold("mapq_min", 20)
        if self.baseq_min is None:
            self.baseq_min = config.get_threshold("baseq_min", 20)

        # Validate threshold values (now they should all be set)
        assert self.dp_min is not None
        assert self.maf_min is not None
        assert self.dp_cap is not None
        assert self.mapq_min is not None
        assert self.baseq_min is not None

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

    @classmethod
    def from_config(
        cls, config: Optional["SsiambConfig"] = None, **overrides: Any
    ) -> "Thresholds":
        """
        Create Thresholds from configuration with optional overrides.

        Args:
            config: Configuration object (uses global if None)
            **overrides: Explicit threshold values to override

        Returns:
            Thresholds instance with config defaults and overrides applied
        """
        if config is None:
            from .config import get_config

            config = get_config()

        # Start with config defaults
        kwargs = {
            "dp_min": config.get_threshold("dp_min", 10),
            "maf_min": config.get_threshold("maf_min", 0.1),
            "dp_cap": config.get_threshold("dp_cap", 100),
            "mapq_min": config.get_threshold("mapq_min", 20),
            "baseq_min": config.get_threshold("baseq_min", 20),
        }

        # Apply overrides
        kwargs.update(overrides)

        return cls(**kwargs)


@dataclass
class SummaryRow:
    """
    Single row of the ambiguous_summary.tsv output.

    Based on spec.md §4.1 - Primary Output Format.
    Schema: sample, mode, mapper, caller, dp_min, maf_min, dp_cap, denom_policy, callable_bases,
    genome_length, breadth_10x, ambiguous_snv_count, ambiguous_snv_per_mb,
    ambiguous_indel_count, ambiguous_del_count,
    ref_label, ref_accession, bracken_species, bracken_frac, bracken_reads,
    alias_used, reused_bam, reused_vcf, runtime_sec, tool_version
    """

    # Core run parameters
    sample: str
    mode: str  # self, ref, summarize
    mapper: str  # minimap2, bwa-mem2
    caller: str  # bbtools, bcftools
    dp_min: int
    maf_min: float
    dp_cap: int
    denom_policy: str  # e.g. "exclude_dups"

    # Denominator metrics
    callable_bases: int
    genome_length: int
    breadth_10x: float  # 4 decimal places per spec

    # Numerator counts
    ambiguous_snv_count: int
    ambiguous_snv_per_mb: float  # 2 decimal places per spec
    ambiguous_indel_count: int
    ambiguous_del_count: int

    # Reference information
    ref_label: str  # SpeciesKey|Accession or just SpeciesKey
    ref_accession: str  # GCF_/GCA_ accession or "NA"

    # Bracken fields (NA if not used)
    bracken_species: str
    bracken_frac: float
    bracken_reads: int

    # Reuse information
    alias_used: str  # "NA" if no alias
    reused_bam: bool
    reused_vcf: bool

    # Runtime information
    runtime_sec: float
    tool_version: str

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for TSV output."""
        return {
            "sample": self.sample,
            "mode": self.mode,
            "mapper": self.mapper,
            "caller": self.caller,
            "dp_min": self.dp_min,
            "maf_min": self.maf_min,
            "dp_cap": self.dp_cap,
            "denom_policy": self.denom_policy,
            "callable_bases": self.callable_bases,
            "genome_length": self.genome_length,
            "breadth_10x": f"{self.breadth_10x:.4f}",
            "ambiguous_snv_count": self.ambiguous_snv_count,
            "ambiguous_snv_per_mb": f"{self.ambiguous_snv_per_mb:.2f}",
            "ambiguous_indel_count": self.ambiguous_indel_count,
            "ambiguous_del_count": self.ambiguous_del_count,
            "ref_label": self.ref_label,
            "ref_accession": self.ref_accession,
            "bracken_species": self.bracken_species,
            "bracken_frac": self.bracken_frac,
            "bracken_reads": self.bracken_reads,
            "alias_used": self.alias_used,
            "reused_bam": self.reused_bam,
            "reused_vcf": self.reused_vcf,
            "runtime_sec": self.runtime_sec,
            "tool_version": self.tool_version,
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

    # Tool options
    bbtools_mem: Optional[str] = None

    # Calling options
    require_pass: bool = False

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
    tsv_mode: TSVMode = TSVMode.OVERWRITE

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
