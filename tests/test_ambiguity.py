from __future__ import annotations

from pathlib import Path

from ssiamb.ambiguity import matrix_counts, summarize_records
from ssiamb.vcf import read_vcf


def write_vcf(tmp_path: Path, body: str) -> Path:
    path = tmp_path / "input.vcf"
    path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                body.strip(),
                "",
            ]
        ),
        encoding="utf-8",
    )
    return path


def summarize(tmp_path: Path, body: str):
    return summarize_records(read_vcf(write_vcf(tmp_path, body)), sample="S")


def test_af_010_dp_10_counts(tmp_path: Path) -> None:
    result = summarize(tmp_path, "ctg1\t10\t.\tA\tG\t.\tPASS\tDP=10;AF=0.10")
    assert result.counted_loci == 1
    assert result.loci[0].maf == 0.10


def test_af_0095_does_not_round_into_count(tmp_path: Path) -> None:
    result = summarize(tmp_path, "ctg1\t10\t.\tA\tG\t.\tPASS\tDP=10;AF=0.095")
    assert result.counted_loci == 0
    assert result.loci[0].maf == 0.095


def test_high_alt_af_converts_to_maf(tmp_path: Path) -> None:
    result = summarize(tmp_path, "ctg1\t10\t.\tA\tG\t.\tPASS\tDP=10;AF=0.90")
    assert result.counted_loci == 1
    assert round(result.loci[0].maf, 6) == 0.10


def test_same_position_af_is_summed_before_counting(tmp_path: Path) -> None:
    result = summarize(
        tmp_path,
        "\n".join(
            [
                "ctg1\t10\t.\tA\tG\t.\tPASS\tDP=10;AF=0.06",
                "ctg1\t10\t.\tA\tT\t.\tPASS\tDP=11;AF=0.05",
            ]
        ),
    )
    assert result.counted_loci == 1
    assert result.loci[0].depth == 11
    assert round(result.loci[0].summed_alt_af, 6) == 0.11


def test_same_position_af_is_capped(tmp_path: Path) -> None:
    result = summarize(
        tmp_path,
        "\n".join(
            [
                "ctg1\t10\t.\tA\tG\t.\tPASS\tDP=10;AF=0.80",
                "ctg1\t10\t.\tA\tT\t.\tPASS\tDP=10;AF=0.70",
            ]
        ),
    )
    assert result.counted_loci == 0
    assert result.loci[0].summed_alt_af == 1.0
    assert result.loci[0].maf == 0.0


def test_filtered_records_stay_in_diagnostics(tmp_path: Path) -> None:
    result = summarize(tmp_path, "ctg1\t10\t.\tA\tG\t.\tq10\tDP=10;AF=0.10")
    assert result.counted_loci == 1
    assert result.loci[0].filters == ("q10",)


def test_low_depth_is_reported_but_not_counted(tmp_path: Path) -> None:
    result = summarize(tmp_path, "ctg1\t10\t.\tA\tG\t.\tPASS\tDP=9;AF=0.20")
    assert result.counted_loci == 0
    assert result.candidate_loci == 1
    assert result.loci[0].depth == 9


def test_missing_info_af_falls_back_to_ad(tmp_path: Path) -> None:
    path = tmp_path / "input.vcf"
    path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS",
                "ctg1\t10\t.\tA\tG\t.\tPASS\tDP=10\tAD\t9,1",
                "",
            ]
        ),
        encoding="utf-8",
    )
    result = summarize_records(read_vcf(path), sample="S")
    assert result.counted_loci == 1
    assert result.loci[0].evidence_sources == ("FORMAT/AD",)


def test_matrix_uses_floor_bins(tmp_path: Path) -> None:
    result = summarize(tmp_path, "ctg1\t10\t.\tA\tG\t.\tPASS\tDP=10;AF=0.095")
    matrix = matrix_counts(result)
    row_10 = matrix[9]
    assert row_10["maf_09"] == 1
    assert row_10["maf_10"] == 0
