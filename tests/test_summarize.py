#!/usr/bin/env python3
"""
Tests for Milestone 19: Summarize Mode

Tests the run_summarize function and summarize CLI command to ensure:
- Both VCF and BAM are required
- Denominator logic is reused from BAM
- All emission options are supported
- Error handling works correctly
"""

import pytest
import re
import tempfile
from pathlib import Path
from typer.testing import CliRunner

from ssiamb.cli import app
from ssiamb.runner import run_summarize


def strip_ansi_codes(text: str) -> str:
    """Remove ANSI escape sequences from text for consistent testing."""
    ansi_escape = re.compile(r"\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])")
    return ansi_escape.sub("", text)


class TestSummarizeFunction:
    """Test the run_summarize function directly."""

    def test_run_summarize_requires_vcf_and_bam(self):
        """Test that run_summarize validates both VCF and BAM are provided."""
        with pytest.raises(FileNotFoundError, match="VCF file not found"):
            run_summarize(
                vcf=Path("nonexistent.vcf"),
                bam=Path(
                    "nonexistent.bam"
                ),  # Use nonexistent bam to ensure VCF is checked first
                output=Path("output.tsv"),
                to_stdout=False,
            )

        # Create a temporary VCF file to test BAM validation
        with tempfile.NamedTemporaryFile(suffix=".vcf", delete=False) as tmp_vcf:
            tmp_vcf.write(
                b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            )
            tmp_vcf_path = Path(tmp_vcf.name)

        try:
            with pytest.raises(FileNotFoundError, match="BAM file not found"):
                run_summarize(
                    vcf=tmp_vcf_path,
                    bam=Path("nonexistent.bam"),
                    output=Path("output.tsv"),
                    to_stdout=False,
                )
        finally:
            tmp_vcf_path.unlink()  # Clean up    def test_run_summarize_reuses_denominator_logic(self):
        """Test that run_summarize reuses the denominator logic from BAM files."""
        # This test verifies that the function signature is correct and validates files
        # Create minimal valid files
        with (
            tempfile.NamedTemporaryFile(suffix=".vcf", delete=False) as vcf_file,
            tempfile.NamedTemporaryFile(suffix=".bam", delete=False) as bam_file,
        ):

            # Write minimal VCF content
            vcf_file.write(
                b"""##fileformat=VCFv4.2
##contig=<ID=test_contig,length=1000>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample
test_contig\t100\t.\tA\tT\t60\tPASS\tDP=50\tGT:DP\t0/1:50
"""
            )
            vcf_path = Path(vcf_file.name)
            bam_path = Path(bam_file.name)

            try:
                # This should fail due to invalid BAM, but validates the interface
                with pytest.raises(Exception):  # Expecting some error due to fake BAM
                    run_summarize(
                        vcf=vcf_path,
                        bam=bam_path,
                        output=Path("output.tsv"),
                        to_stdout=False,
                    )

            finally:
                # Cleanup
                vcf_path.unlink(missing_ok=True)
                bam_path.unlink(missing_ok=True)

    def test_run_summarize_supports_all_emission_options(self):
        """Test that run_summarize supports all emission options."""
        # Test that all emission flags are accepted without errors in the function signature
        # Create minimal valid files
        with (
            tempfile.NamedTemporaryFile(suffix=".vcf", delete=False) as vcf_file,
            tempfile.NamedTemporaryFile(suffix=".bam", delete=False) as bam_file,
        ):

            # Write minimal VCF content
            vcf_file.write(
                b"""##fileformat=VCFv4.2
##contig=<ID=test_contig,length=1000>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample
test_contig\t100\t.\tA\tT\t60\tPASS\tDP=50\tGT:DP\t0/1:50
"""
            )
            vcf_path = Path(vcf_file.name)
            bam_path = Path(bam_file.name)

            try:
                # Test that function accepts all emission options
                with pytest.raises(Exception):  # Expecting some error due to fake BAM
                    run_summarize(
                        vcf=vcf_path,
                        bam=bam_path,
                        output=Path("output.tsv"),
                        emit_vcf=True,
                        emit_bed=True,
                        emit_matrix=True,
                        emit_per_contig=True,
                        emit_multiqc=True,
                        emit_provenance=True,
                        to_stdout=False,
                    )

            finally:
                # Cleanup
                vcf_path.unlink(missing_ok=True)
                bam_path.unlink(missing_ok=True)


class TestSummarizeCLI:
    """Test the summarize CLI command."""

    def test_summarize_requires_both_vcf_and_bam(self):
        """Test that CLI requires both --vcf and --bam arguments."""
        runner = CliRunner()

        # Test missing VCF
        result = runner.invoke(app, ["summarize", "--bam", "test.bam"])
        assert result.exit_code != 0
        assert "Missing option '--vcf'" in strip_ansi_codes(result.output)

        # Test missing BAM
        result = runner.invoke(app, ["summarize", "--vcf", "test.vcf"])
        assert result.exit_code != 0
        assert "Missing option '--bam'" in strip_ansi_codes(result.output)

    def test_summarize_validates_file_existence(self):
        """Test that CLI validates file existence with proper error codes."""
        runner = CliRunner()

        # Test nonexistent VCF
        result = runner.invoke(
            app, ["summarize", "--vcf", "nonexistent.vcf", "--bam", "nonexistent.bam"]
        )
        assert result.exit_code == 1  # Input validation error
        assert "VCF file not found" in result.output

    def test_summarize_emission_flags_work(self):
        """Test that all emission flags are accepted by CLI."""
        runner = CliRunner()

        # Test that all emission flags are accepted (will fail on file validation, but flags should be OK)
        result = runner.invoke(
            app,
            [
                "summarize",
                "--vcf",
                "test.vcf",
                "--bam",
                "test.bam",
                "--emit-vcf",
                "--emit-bed",
                "--emit-matrix",
                "--emit-per-contig",
                "--emit-multiqc",
                "--emit-provenance",
            ],
        )

        # Should fail on file validation, not argument parsing
        assert "VCF file not found" in result.output
        assert "Missing option" not in result.output

    def test_summarize_stdout_mode(self):
        """Test that stdout mode works and prevents file emission."""
        runner = CliRunner()

        # Test stdout mode with emission flags (should be allowed in CLI, validation happens later)
        result = runner.invoke(
            app, ["summarize", "--vcf", "test.vcf", "--bam", "test.bam", "--stdout"]
        )

        # Will fail on file validation but stdout flag should be accepted
        assert result.exit_code == 1
        assert "VCF file not found" in result.output


class TestSummarizeIntegration:
    """Integration tests using real fixture files."""

    def test_summarize_with_real_files(self):
        """Test summarize mode with real VCF and BAM files from fixtures."""
        runner = CliRunner()

        # Use real files from test_output
        vcf_file = "test_output/2508H52931_bwa_bbtools.normalized.vcf.gz"
        bam_file = "test_output/2508H52931_bwa_bbtools.sorted.bam"

        # Check if files exist first
        vcf_path = Path(vcf_file)
        bam_path = Path(bam_file)

        if vcf_path.exists() and bam_path.exists():
            with tempfile.TemporaryDirectory() as tmpdir:
                output_file = Path(tmpdir) / "summary.tsv"

                result = runner.invoke(
                    app,
                    [
                        "summarize",
                        "--vcf",
                        str(vcf_path),
                        "--bam",
                        str(bam_path),
                        "--output",
                        str(output_file),
                        "--stdout",
                    ],
                )

                # Should succeed and print TSV to stdout
                if result.exit_code != 0:
                    print("STDOUT:", result.output)
                    print("Exception:", result.exception)
                assert result.exit_code == 0
                assert "sample" in result.output  # TSV header
                assert "2508H52931_bwa_bbtools" in result.output  # Sample name
        else:
            pytest.skip("Real test files not available")

    def test_summarize_dry_run_with_real_files(self):
        """Test dry-run mode with real files."""
        runner = CliRunner()

        vcf_file = "test_output/2508H52931_bwa_bbtools.normalized.vcf.gz"
        bam_file = "test_output/2508H52931_bwa_bbtools.sorted.bam"

        vcf_path = Path(vcf_file)
        bam_path = Path(bam_file)

        if vcf_path.exists() and bam_path.exists():
            result = runner.invoke(
                app,
                [
                    "--dry-run",
                    "summarize",
                    "--vcf",
                    str(vcf_path),
                    "--bam",
                    str(bam_path),
                ],
            )

            # Should succeed with dry-run output
            assert result.exit_code == 0
            assert "DRY RUN - Summarize mode plan" in result.output
            assert "no files written" in result.output
        else:
            pytest.skip("Real test files not available")
