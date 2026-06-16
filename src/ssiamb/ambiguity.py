from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Iterable

from .models import AmbiguityResult, LocusSummary, RecordEvidence, VcfRecord


def _first_float(value: str | bool | None) -> float | None:
    if value is None or isinstance(value, bool):
        return None
    first = value.split(",", 1)[0]
    if first in {"", "."}:
        return None
    try:
        return float(first)
    except ValueError:
        return None


def _float_list(value: str | bool | None) -> list[float]:
    if value is None or isinstance(value, bool):
        return []
    floats: list[float] = []
    for item in value.split(","):
        if item in {"", "."}:
            continue
        try:
            floats.append(float(item))
        except ValueError:
            continue
    return floats


def _int_list(value: str | bool | None) -> list[int]:
    if value is None or isinstance(value, bool):
        return []
    ints: list[int] = []
    for item in value.replace(",", ":").split(":"):
        if item in {"", "."}:
            continue
        try:
            ints.append(int(float(item)))
        except ValueError:
            continue
    return ints


def extract_depth(record: VcfRecord) -> int | None:
    info_dp = _first_float(record.info.get("DP"))
    if info_dp is not None:
        return max(0, int(info_dp))

    sample_dp = _first_float(record.sample_values.get("DP"))
    if sample_dp is not None:
        return max(0, int(sample_dp))

    ad_values = _int_list(record.sample_values.get("AD"))
    if ad_values:
        return max(0, sum(ad_values))

    dp4_values = _int_list(record.info.get("DP4"))
    if dp4_values:
        return max(0, sum(dp4_values))

    return None


def extract_alt_af(record: VcfRecord) -> tuple[float | None, str]:
    info_af = _float_list(record.info.get("AF"))
    if info_af:
        return min(1.0, max(0.0, sum(info_af))), "INFO/AF"

    sample_af = _float_list(record.sample_values.get("AF"))
    if sample_af:
        return min(1.0, max(0.0, sum(sample_af))), "FORMAT/AF"

    ad_values = _int_list(record.sample_values.get("AD"))
    if len(ad_values) >= 2:
        total = sum(ad_values)
        if total > 0:
            return min(1.0, max(0.0, sum(ad_values[1:]) / total)), "FORMAT/AD"

    dp4_values = _int_list(record.info.get("DP4"))
    if len(dp4_values) >= 4:
        total = sum(dp4_values[:4])
        if total > 0:
            return min(1.0, max(0.0, sum(dp4_values[2:4]) / total)), "INFO/DP4"

    return None, "missing"


def evidence_from_record(record: VcfRecord) -> RecordEvidence:
    if record.is_snv:
        variant_type = "snv"
    elif any(len(alt) != len(record.ref) for alt in record.alts):
        variant_type = "indel"
    else:
        variant_type = "complex"

    alt_af, af_source = extract_alt_af(record)
    return RecordEvidence(
        record=record,
        depth=extract_depth(record),
        alt_af=alt_af,
        af_source=af_source,
        variant_type=variant_type,
    )


def summarize_records(
    records: Iterable[VcfRecord],
    *,
    sample: str,
    dp_min: int = 10,
    maf_min: float = 0.10,
    dp_cap: int = 100,
) -> AmbiguityResult:
    grouped: dict[tuple[str, int], list[RecordEvidence]] = defaultdict(list)
    for record in records:
        grouped[(record.chrom, record.pos)].append(evidence_from_record(record))

    loci: list[LocusSummary] = []
    for (chrom, pos), evidences in sorted(grouped.items()):
        snv_evidences = [
            evidence for evidence in evidences if evidence.variant_type == "snv"
        ]
        depth_values = [
            evidence.depth for evidence in evidences if evidence.depth is not None
        ]
        depth = max(depth_values) if depth_values else 0

        summed_alt_af = min(
            1.0,
            sum(
                evidence.alt_af
                for evidence in snv_evidences
                if evidence.alt_af is not None
            ),
        )
        maf = min(summed_alt_af, 1.0 - summed_alt_af)
        has_snv_evidence = bool(snv_evidences)
        counted = has_snv_evidence and depth >= dp_min and (maf + 1e-12) >= maf_min

        loci.append(
            LocusSummary(
                sample=sample,
                chrom=chrom,
                pos=pos,
                depth=depth,
                summed_alt_af=summed_alt_af,
                maf=maf,
                counted=counted,
                filters=tuple(sorted({evidence.record.filt for evidence in evidences})),
                evidence_sources=tuple(
                    sorted({evidence.af_source for evidence in evidences})
                ),
                record_count=len(evidences),
                variant_types=tuple(
                    sorted({evidence.variant_type for evidence in evidences})
                ),
            )
        )

    return AmbiguityResult(
        sample=sample,
        dp_min=dp_min,
        maf_min=maf_min,
        dp_cap=dp_cap,
        loci=tuple(loci),
    )


def matrix_counts(result: AmbiguityResult) -> list[dict[str, int]]:
    snv_loci = [locus for locus in result.loci if "snv" in locus.variant_types]
    rows: list[dict[str, int]] = []

    for dp_threshold in range(1, result.dp_cap + 1):
        row: dict[str, int] = {"dp_threshold": dp_threshold}
        for maf_threshold in range(0, 51):
            row[f"maf_{maf_threshold:02d}"] = 0
        rows.append(row)

    for locus in snv_loci:
        depth_bin = max(1, min(result.dp_cap, locus.depth))
        maf_bin = max(0, min(50, math.floor((locus.maf * 100.0) + 1e-12)))
        for row in rows[:depth_bin]:
            for maf_threshold in range(0, maf_bin + 1):
                row[f"maf_{maf_threshold:02d}"] += 1

    return rows
