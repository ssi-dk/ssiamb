from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class VcfRecord:
    chrom: str
    pos: int
    record_id: str
    ref: str
    alts: tuple[str, ...]
    qual: str
    filt: str
    info: dict[str, str | bool]
    format_keys: tuple[str, ...] = ()
    sample_values: dict[str, str] = field(default_factory=dict)
    line_number: int = 0

    @property
    def is_snv(self) -> bool:
        return len(self.ref) == 1 and any(len(alt) == 1 for alt in self.alts)


@dataclass(frozen=True)
class RecordEvidence:
    record: VcfRecord
    depth: int | None
    alt_af: float | None
    af_source: str
    variant_type: str


@dataclass(frozen=True)
class LocusSummary:
    sample: str
    chrom: str
    pos: int
    depth: int
    summed_alt_af: float
    maf: float
    counted: bool
    filters: tuple[str, ...]
    evidence_sources: tuple[str, ...]
    record_count: int
    variant_types: tuple[str, ...]


@dataclass(frozen=True)
class AmbiguityResult:
    sample: str
    dp_min: int
    maf_min: float
    dp_cap: int
    loci: tuple[LocusSummary, ...]

    @property
    def candidate_loci(self) -> int:
        return sum("snv" in locus.variant_types for locus in self.loci)

    @property
    def counted_loci(self) -> int:
        return sum(locus.counted for locus in self.loci)
