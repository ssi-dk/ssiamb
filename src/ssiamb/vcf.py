from __future__ import annotations

import gzip
from collections.abc import Iterator
from pathlib import Path

from .models import VcfRecord


class VcfParseError(ValueError):
    """Raised when a VCF line cannot be parsed."""


def open_text(path: Path) -> Iterator[str]:
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            yield from handle
    else:
        with path.open("r", encoding="utf-8") as handle:
            yield from handle


def parse_info(info_text: str) -> dict[str, str | bool]:
    if info_text in {"", "."}:
        return {}

    info: dict[str, str | bool] = {}
    for item in info_text.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
        else:
            info[item] = True
    return info


def parse_record(line: str, line_number: int = 0) -> VcfRecord:
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 8:
        raise VcfParseError(f"VCF line {line_number} has fewer than 8 columns")

    chrom, pos_text, record_id, ref, alt_text, qual, filt, info_text = fields[:8]
    try:
        pos = int(pos_text)
    except ValueError as exc:
        raise VcfParseError(
            f"VCF line {line_number} has invalid POS: {pos_text}"
        ) from exc

    format_keys: tuple[str, ...] = ()
    sample_values: dict[str, str] = {}
    if len(fields) >= 10:
        format_keys = tuple(fields[8].split(":"))
        values = fields[9].split(":")
        sample_values = dict(zip(format_keys, values, strict=False))

    alts = tuple(alt for alt in alt_text.split(",") if alt and alt != ".")
    return VcfRecord(
        chrom=chrom,
        pos=pos,
        record_id=record_id,
        ref=ref,
        alts=alts,
        qual=qual,
        filt=filt,
        info=parse_info(info_text),
        format_keys=format_keys,
        sample_values=sample_values,
        line_number=line_number,
    )


def read_vcf(path: Path) -> Iterator[VcfRecord]:
    for line_number, line in enumerate(open_text(path), start=1):
        if not line.strip() or line.startswith("#"):
            continue
        yield parse_record(line, line_number=line_number)
