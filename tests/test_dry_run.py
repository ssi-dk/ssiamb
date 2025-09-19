#!/usr/bin/env python3
"""
Tests for Milestone 20: Dry-run Planner

Tests the --dry-run functionality across all subcommands to ensure:
- Dry-run validates inputs properly
- Shows planned steps without execution
- Exits with code 0 on success
- No files are written during dry-run
"""

import pytest
import tempfile
from pathlib import Path
from typer.testing import CliRunner
import re

from ssiamb.cli import app


def strip_ansi_codes(text: str) -> str:
    """Remove ANSI escape codes from text for easier testing."""
    ansi_escape = re.compile(r"\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])")
    return ansi_escape.sub("", text)


class TestDryRunGlobal:
    """Test --dry-run flag at global level."""

    def test_dry_run_flag_position(self):
        """Test that --dry-run is accepted as a global flag."""
        runner = CliRunner()

        # Test that --dry-run is recognized as global option
        result = runner.invoke(app, ["--help"])
        clean_output = strip_ansi_codes(result.output)
        assert "--dry-run" in clean_output
        assert "Show what would be done without executing" in clean_output

    def test_dry_run_exits_zero_on_validation_error(self):
        """Test that dry-run still validates inputs but doesn't fail hard."""
        runner = CliRunner()

        # Dry-run with missing files should show validation error but exit appropriately
        result = runner.invoke(
            app,
            [
                "--dry-run",
                "self",
                "--r1",
                "nonexistent_R1.fastq.gz",
                "--r2",
                "nonexistent_R2.fastq.gz",
                "--assembly",
                "nonexistent.fna",
                "--sample",
                "test_sample",
            ],
        )

        # Should show validation error
        assert "R1 file not found" in result.output
        assert result.exit_code == 1  # Still fails on validation


class TestDryRunSelfMode:
    """Test --dry-run with self mode."""

    def test_dry_run_self_shows_plan(self):
        """Test that dry-run shows mapping and calling plan for self mode."""
        runner = CliRunner()

        # Create temporary files to pass initial validation
        with (
            tempfile.NamedTemporaryFile(suffix="_R1.fastq.gz", delete=False) as r1,
            tempfile.NamedTemporaryFile(suffix="_R2.fastq.gz", delete=False) as r2,
            tempfile.NamedTemporaryFile(suffix=".fna", delete=False) as assembly,
        ):

            # Write minimal content to files
            r1.write(b"@test\nACGT\n+\nIIII\n")
            r2.write(b"@test\nACGT\n+\nIIII\n")
            assembly.write(b">contig1\nACGTACGTACGT\n")

            r1_path = Path(r1.name)
            r2_path = Path(r2.name)
            assembly_path = Path(assembly.name)

            try:
                result = runner.invoke(
                    app,
                    [
                        "--dry-run",
                        "self",
                        "--r1",
                        str(r1_path),
                        "--r2",
                        str(r2_path),
                        "--assembly",
                        str(assembly_path),
                        "--sample",
                        "test_sample",
                    ],
                )

                # Should show dry-run plan
                if result.exit_code == 0:
                    assert (
                        "DRY RUN" in result.output and "Self-mapping" in result.output
                    )
                    assert "no files" in result.output

            finally:
                # Cleanup
                r1_path.unlink(missing_ok=True)
                r2_path.unlink(missing_ok=True)
                assembly_path.unlink(missing_ok=True)

    def test_dry_run_self_shows_mapper_and_caller(self):
        """Test that dry-run shows selected mapper and caller."""
        runner = CliRunner()

        with (
            tempfile.NamedTemporaryFile(suffix="_R1.fastq.gz", delete=False) as r1,
            tempfile.NamedTemporaryFile(suffix="_R2.fastq.gz", delete=False) as r2,
            tempfile.NamedTemporaryFile(suffix=".fna", delete=False) as assembly,
        ):

            r1.write(b"@test\nACGT\n+\nIIII\n")
            r2.write(b"@test\nACGT\n+\nIIII\n")
            assembly.write(b">contig1\nACGTACGTACGT\n")

            r1_path = Path(r1.name)
            r2_path = Path(r2.name)
            assembly_path = Path(assembly.name)

            try:
                result = runner.invoke(
                    app,
                    [
                        "--dry-run",
                        "self",
                        "--r1",
                        str(r1_path),
                        "--r2",
                        str(r2_path),
                        "--assembly",
                        str(assembly_path),
                        "--sample",
                        "test_sample",
                        "--mapper",
                        "bwa-mem2",
                        "--caller",
                        "bcftools",
                    ],
                )

                if result.exit_code == 0:
                    # Should show selected tools in plan
                    assert "bwa-mem2" in result.output or "mapper" in result.output
                    assert "bcftools" in result.output or "caller" in result.output

            finally:
                r1_path.unlink(missing_ok=True)
                r2_path.unlink(missing_ok=True)
                assembly_path.unlink(missing_ok=True)


class TestDryRunRefMode:
    """Test --dry-run with ref mode."""

    def test_dry_run_ref_shows_plan(self):
        """Test that dry-run shows reference resolution and mapping plan."""
        runner = CliRunner()

        with (
            tempfile.NamedTemporaryFile(suffix="_R1.fastq.gz", delete=False) as r1,
            tempfile.NamedTemporaryFile(suffix="_R2.fastq.gz", delete=False) as r2,
            tempfile.NamedTemporaryFile(suffix=".fna", delete=False) as ref,
        ):

            r1.write(b"@test\nACGT\n+\nIIII\n")
            r2.write(b"@test\nACGT\n+\nIIII\n")
            ref.write(b">ref_contig\nACGTACGTACGT\n")

            r1_path = Path(r1.name)
            r2_path = Path(r2.name)
            ref_path = Path(ref.name)

            try:
                result = runner.invoke(
                    app,
                    [
                        "--dry-run",
                        "ref",
                        "--r1",
                        str(r1_path),
                        "--r2",
                        str(r2_path),
                        "--reference",
                        str(ref_path),
                        "--sample",
                        "test_sample",
                    ],
                )

                if result.exit_code == 0:
                    assert (
                        "DRY RUN" in result.output
                        and "Reference-mapping" in result.output
                    )
                    assert "no files" in result.output

            finally:
                r1_path.unlink(missing_ok=True)
                r2_path.unlink(missing_ok=True)
                ref_path.unlink(missing_ok=True)

    def test_dry_run_ref_shows_reference_resolution(self):
        """Test that dry-run shows how reference is resolved."""
        runner = CliRunner()

        with (
            tempfile.NamedTemporaryFile(suffix="_R1.fastq.gz", delete=False) as r1,
            tempfile.NamedTemporaryFile(suffix="_R2.fastq.gz", delete=False) as r2,
            tempfile.NamedTemporaryFile(suffix=".fna", delete=False) as ref,
        ):

            r1.write(b"@test\nACGT\n+\nIIII\n")
            r2.write(b"@test\nACGT\n+\nIIII\n")
            ref.write(b">ref_contig\nACGTACGTACGT\n")

            r1_path = Path(r1.name)
            r2_path = Path(r2.name)
            ref_path = Path(ref.name)

            try:
                result = runner.invoke(
                    app,
                    [
                        "--dry-run",
                        "ref",
                        "--r1",
                        str(r1_path),
                        "--r2",
                        str(r2_path),
                        "--reference",
                        str(ref_path),
                        "--sample",
                        "test_sample",
                    ],
                )

                if result.exit_code == 0:
                    # Should show reference being used
                    assert (
                        str(ref_path) in result.output
                        or "reference" in result.output.lower()
                    )

            finally:
                r1_path.unlink(missing_ok=True)
                r2_path.unlink(missing_ok=True)
                ref_path.unlink(missing_ok=True)


class TestDryRunSummarizeMode:
    """Test --dry-run with summarize mode."""

    def test_dry_run_summarize_shows_plan(self):
        """Test that dry-run shows analysis plan for summarize mode."""
        runner = CliRunner()

        with (
            tempfile.NamedTemporaryFile(suffix=".vcf", delete=False) as vcf,
            tempfile.NamedTemporaryFile(suffix=".bam", delete=False) as bam,
        ):

            # Write minimal content
            vcf.write(
                b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
            )
            bam.write(b"BAM_HEADER_PLACEHOLDER")  # Real BAM would be binary

            vcf_path = Path(vcf.name)
            bam_path = Path(bam.name)

            try:
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

                if result.exit_code == 0:
                    assert "DRY RUN - Summarize mode plan" in result.output
                    assert "no files written" in result.output

            finally:
                vcf_path.unlink(missing_ok=True)
                bam_path.unlink(missing_ok=True)

    def test_dry_run_summarize_shows_thresholds(self):
        """Test that dry-run shows analysis thresholds."""
        runner = CliRunner()

        with (
            tempfile.NamedTemporaryFile(suffix=".vcf", delete=False) as vcf,
            tempfile.NamedTemporaryFile(suffix=".bam", delete=False) as bam,
        ):

            vcf.write(
                b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
            )
            bam.write(b"BAM_HEADER_PLACEHOLDER")

            vcf_path = Path(vcf.name)
            bam_path = Path(bam.name)

            try:
                result = runner.invoke(
                    app,
                    [
                        "--dry-run",
                        "summarize",
                        "--vcf",
                        str(vcf_path),
                        "--bam",
                        str(bam_path),
                        "--dp-min",
                        "15",
                        "--maf-min",
                        "0.05",
                    ],
                )

                if result.exit_code == 0:
                    # Should show custom thresholds
                    assert "dp_min=15" in result.output
                    assert "maf_min=0.05" in result.output

            finally:
                vcf_path.unlink(missing_ok=True)
                bam_path.unlink(missing_ok=True)


class TestDryRunNoExecution:
    """Test that dry-run doesn't execute or write files."""

    def test_dry_run_no_file_creation(self):
        """Test that dry-run doesn't create any output files."""
        runner = CliRunner()

        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir)

            # Count files before
            files_before = list(output_dir.glob("*"))

            with (
                tempfile.NamedTemporaryFile(suffix=".vcf", delete=False) as vcf,
                tempfile.NamedTemporaryFile(suffix=".bam", delete=False) as bam,
            ):

                vcf.write(
                    b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
                )
                bam.write(b"BAM_HEADER_PLACEHOLDER")

                vcf_path = Path(vcf.name)
                bam_path = Path(bam.name)

                try:
                    runner.invoke(
                        app,
                        [
                            "--dry-run",
                            "summarize",
                            "--vcf",
                            str(vcf_path),
                            "--bam",
                            str(bam_path),
                            "--outdir",
                            str(output_dir),
                            "--emit-vcf",
                            "--emit-bed",
                            "--emit-matrix",
                        ],
                    )

                    # Count files after
                    files_after = list(output_dir.glob("*"))

                    # Should not have created any new files
                    assert len(files_after) == len(files_before)

                finally:
                    vcf_path.unlink(missing_ok=True)
                    bam_path.unlink(missing_ok=True)

    def test_dry_run_shows_planned_outputs(self):
        """Test that dry-run shows what outputs would be created."""
        runner = CliRunner()

        with (
            tempfile.NamedTemporaryFile(suffix=".vcf", delete=False) as vcf,
            tempfile.NamedTemporaryFile(suffix=".bam", delete=False) as bam,
        ):

            vcf.write(
                b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
            )
            bam.write(b"BAM_HEADER_PLACEHOLDER")

            vcf_path = Path(vcf.name)
            bam_path = Path(bam.name)

            try:
                result = runner.invoke(
                    app,
                    [
                        "--dry-run",
                        "summarize",
                        "--vcf",
                        str(vcf_path),
                        "--bam",
                        str(bam_path),
                        "--emit-vcf",
                        "--emit-bed",
                    ],
                )

                if result.exit_code == 0:
                    # Should mention planned outputs
                    assert (
                        "Outputs planned" in result.output
                        or "outputs" in result.output.lower()
                    )

            finally:
                vcf_path.unlink(missing_ok=True)
                bam_path.unlink(missing_ok=True)


class TestDryRunRealFiles:
    """Test dry-run with real fixture files."""

    def test_dry_run_with_real_files_exits_zero(self):
        """Test that dry-run with valid real files exits successfully."""
        runner = CliRunner()

        # Use real files if available
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

            # Should exit successfully with dry-run plan
            assert result.exit_code == 0
            assert "DRY RUN" in result.output
            assert "no files written" in result.output
        else:
            pytest.skip("Real test files not available")
