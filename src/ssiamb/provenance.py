"""
Provenance tracking for ssiamb runs.

This module handles collection and output of detailed per-sample provenance
information when --emit-provenance is used.
"""

import json
import hashlib
import os
import sys
from dataclasses import dataclass, asdict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Any

from .models import Mode
from .qc import QCWarning


@dataclass
class ProvenanceInput:
    """Input file information with MD5 hash."""

    path: str
    md5: Optional[str] = None


@dataclass
class ProvenanceReferenceInfo:
    """Reference genome information."""

    species_requested: Optional[str] = None
    alias_applied: Optional[str] = None
    species_final: Optional[str] = None
    fasta_path: Optional[str] = None
    fasta_md5: Optional[str] = None
    source: Optional[Dict[str, Any]] = None


@dataclass
class ProvenanceSpeciesSelection:
    """Species selection from Bracken."""

    bracken_species: Optional[str] = None
    bracken_frac: Optional[float] = None
    bracken_reads: Optional[int] = None
    thresholds: Optional[Dict[str, Any]] = None
    on_fail: Optional[str] = None


@dataclass
class ProvenanceMappingStats:
    """Mapping statistics."""

    total_reads: Optional[int] = None
    mapped_reads: Optional[int] = None
    map_rate: Optional[float] = None
    mean_depth: Optional[float] = None
    breadth_1x: Optional[float] = None
    breadth_10x: Optional[float] = None


@dataclass
class ProvenanceDuplicatePolicy:
    """Duplicate handling policy."""

    denominator_excludes_dups: bool = True
    bam_had_dups_flag: Optional[bool] = None


@dataclass
class ProvenanceContigFilters:
    """Contig filtering settings."""

    min_contig_len: int = 500


@dataclass
class ProvenanceCallerParams:
    """Variant caller parameters."""

    exact_cmdlines: List[str]


@dataclass
class ProvenanceCounts:
    """Analysis counts."""

    ambiguous_snv_count: Optional[int] = None
    ambiguous_indel_count: Optional[int] = None
    ambiguous_del_count: Optional[int] = None
    callable_bases: Optional[int] = None
    genome_length: Optional[int] = None


@dataclass
class ProvenanceGridCell:
    """Grid cell used for analysis."""

    depth: int
    maf_bin: int


@dataclass
class ProvenanceExtrasEmitted:
    """Optional outputs emitted."""

    vcf: bool = False
    bed: bool = False
    matrix: bool = False
    per_contig: bool = False


@dataclass
class ProvenanceRecord:
    """Complete provenance record for a sample."""

    # Tool and environment info
    tool_version: str
    python_version: str
    conda_env: Optional[str]

    # Timing
    started_at: str
    finished_at: str
    runtime_sec: float
    threads: int

    # Sample and analysis
    sample: str
    mode: str
    mapper: str
    caller: str

    # Analysis parameters
    thresholds: Dict[str, Any]
    denom_policy: str
    depth_tool: str

    # Inputs
    inputs: Dict[str, Any]

    # Reference information
    reference_info: Dict[str, Any]

    # Species selection (for ref mode)
    species_selection: Optional[Dict[str, Any]] = None

    # Results
    mapping_stats: Optional[Dict[str, Any]] = None
    duplicate_policy: Optional[Dict[str, Any]] = None
    contig_filters: Optional[Dict[str, Any]] = None
    caller_params: Optional[Dict[str, Any]] = None
    counts: Optional[Dict[str, Any]] = None
    grid_cell_used: Optional[Dict[str, Any]] = None
    extras_emitted: Optional[Dict[str, Any]] = None
    warnings: Optional[List[str]] = None


def calculate_md5(file_path: Path) -> Optional[str]:
    """Calculate MD5 hash of a file."""
    hash_md5 = hashlib.md5()
    try:
        with open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)
        return hash_md5.hexdigest()
    except (FileNotFoundError, PermissionError):
        return None


def get_conda_env() -> Optional[str]:
    """Get the name of the current conda environment."""
    conda_env = os.environ.get("CONDA_DEFAULT_ENV")
    if conda_env and conda_env != "base":
        return conda_env
    return None


def get_tool_version() -> str:
    """Get the current tool version."""
    try:
        from importlib.metadata import version

        return version("ssiamb")
    except Exception:
        return "unknown"


def create_provenance_record(
    sample: str,
    mode: Mode,
    started_at: datetime,
    finished_at: datetime,
    threads: int,
    mapper: str,
    caller: str,
    dp_min: int,
    maf_min: float,
    dp_cap: int,
    denom_policy: str,
    depth_tool: str,
    emit_vcf: bool,
    emit_bed: bool,
    emit_matrix: bool,
    emit_per_contig: bool,
    r1: Optional[Path] = None,
    r2: Optional[Path] = None,
    assembly: Optional[Path] = None,
    reference: Optional[Path] = None,
    bam: Optional[Path] = None,
    vcf: Optional[Path] = None,
    species: Optional[str] = None,
    mapping_stats: Optional[ProvenanceMappingStats] = None,
    species_selection: Optional[ProvenanceSpeciesSelection] = None,
    counts: Optional[ProvenanceCounts] = None,
    warnings: Optional[List[QCWarning]] = None,
    reference_path: Optional[Path] = None,
    reference_species: Optional[str] = None,
) -> ProvenanceRecord:
    """Create a complete provenance record for a sample."""

    # Calculate runtime
    runtime_sec = (finished_at - started_at).total_seconds()

    # Prepare inputs
    inputs = {}
    if r1:
        inputs["r1"] = {"path": str(r1), "md5": calculate_md5(r1)}
    if r2:
        inputs["r2"] = {"path": str(r2), "md5": calculate_md5(r2)}
    if assembly:
        inputs["assembly"] = {"path": str(assembly), "md5": calculate_md5(assembly)}
    if reference:
        inputs["reference"] = {"path": str(reference), "md5": calculate_md5(reference)}
    if bam:
        inputs["bam"] = {"path": str(bam), "md5": calculate_md5(bam)}
    if vcf:
        inputs["vcf"] = {"path": str(vcf), "md5": calculate_md5(vcf)}

    # Reference info
    reference_info = {}
    if mode == Mode.REF:
        reference_info["species_requested"] = species
        reference_info["species_final"] = reference_species
        if reference_path:
            reference_info["fasta_path"] = str(reference_path)
            reference_info["fasta_md5"] = calculate_md5(reference_path)

    # Grid cell used
    grid_cell = {
        "depth": dp_min,
        "maf_bin": int(maf_min * 100),  # Floor of 100 * maf_min
    }

    # Extras emitted
    extras = {
        "vcf": emit_vcf,
        "bed": emit_bed,
        "matrix": emit_matrix,
        "per_contig": emit_per_contig,
    }

    # Convert warnings to string list
    warning_messages = []
    if warnings:
        warning_messages = [f"{w.metric}: {w.message}" for w in warnings]

    return ProvenanceRecord(
        tool_version=get_tool_version(),
        python_version=f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}",
        conda_env=get_conda_env(),
        started_at=started_at.isoformat(),
        finished_at=finished_at.isoformat(),
        runtime_sec=runtime_sec,
        threads=threads,
        sample=sample,
        mode=mode.value,
        mapper=mapper,
        caller=caller,
        thresholds={"dp_min": dp_min, "maf_min": maf_min, "dp_cap": dp_cap},
        denom_policy=denom_policy,
        depth_tool=depth_tool,
        inputs=inputs,
        reference_info=reference_info,
        species_selection=asdict(species_selection) if species_selection else None,
        mapping_stats=asdict(mapping_stats) if mapping_stats else None,
        duplicate_policy=(
            asdict(ProvenanceDuplicatePolicy()) if mode == Mode.SELF else None
        ),
        contig_filters=asdict(ProvenanceContigFilters(min_contig_len=500)),
        caller_params={"exact_cmdlines": []},  # TODO: Collect actual command lines
        counts=asdict(counts) if counts else {},
        grid_cell_used=grid_cell,
        extras_emitted=extras,
        warnings=warning_messages,
    )


def write_provenance_json(records: List[ProvenanceRecord], output_path: Path) -> None:
    """Write provenance records to JSON file."""
    with open(output_path, "w") as f:
        json.dump([asdict(record) for record in records], f, indent=2)
