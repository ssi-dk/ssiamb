from __future__ import annotations

import csv
import json
import sys
from pathlib import Path
from typing import Any

from .ambiguity import matrix_counts
from .models import AmbiguityResult


def write_summary(
    path: Path,
    result: AmbiguityResult,
    *,
    r1: Path | None,
    r2: Path | None,
    assembly: Path | None,
) -> None:
    fields = [
        "sample",
        "ambiguous_sites",
        "dp_min",
        "maf_min",
        "candidate_loci",
        "counted_loci",
        "assembly",
        "r1",
        "r2",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerow(
            {
                "sample": result.sample,
                "ambiguous_sites": result.counted_loci,
                "dp_min": result.dp_min,
                "maf_min": f"{result.maf_min:.4g}",
                "candidate_loci": result.candidate_loci,
                "counted_loci": result.counted_loci,
                "assembly": str(assembly) if assembly else "",
                "r1": str(r1) if r1 else "",
                "r2": str(r2) if r2 else "",
            }
        )


def write_loci(path: Path, result: AmbiguityResult) -> None:
    fields = [
        "sample",
        "chrom",
        "pos",
        "depth",
        "summed_alt_af",
        "maf",
        "counted",
        "filters",
        "evidence_sources",
        "record_count",
        "variant_types",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for locus in result.loci:
            writer.writerow(
                {
                    "sample": locus.sample,
                    "chrom": locus.chrom,
                    "pos": locus.pos,
                    "depth": locus.depth,
                    "summed_alt_af": f"{locus.summed_alt_af:.6g}",
                    "maf": f"{locus.maf:.6g}",
                    "counted": "true" if locus.counted else "false",
                    "filters": ",".join(locus.filters),
                    "evidence_sources": ",".join(locus.evidence_sources),
                    "record_count": locus.record_count,
                    "variant_types": ",".join(locus.variant_types),
                }
            )


def write_matrix(path: Path, result: AmbiguityResult) -> None:
    rows = matrix_counts(result)
    fields = ["dp_threshold", *[f"maf_{i:02d}" for i in range(0, 51)]]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_provenance(path: Path, payload: dict[str, Any]) -> None:
    payload = {
        "python": sys.version.split()[0],
        **payload,
    }
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
