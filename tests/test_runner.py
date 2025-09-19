"""
Tests for runner.py workflow orchestration.

Tests basic functions in runner.py that can be easily tested.
"""

import pytest
import tempfile
from pathlib import Path
from unittest.mock import patch, MagicMock

from src.ssiamb.runner import (
    _extract_version,
    extract_ref_accession,
    detect_reuse_from_plan,
    create_run_plan,
)
from src.ssiamb.models import Mode, Mapper, Caller


class TestVersionExtraction:
    """Test version extraction functionality."""

    def test_extract_version_simple_version(self):
        """Test version extraction from simple version strings."""
        test_cases = [
            ("2.30", "2.30"),
            ("1.17", "1.17"),
            ("version 1.2.3", "1.2.3"),
            ("v2.5.0", "2.5.0"),
            ("no version here", None),
        ]

        for output, expected in test_cases:
            result = _extract_version(output)
            assert result == expected, f"Failed for input: {output}"


class TestReferenceAccession:
    """Test reference accession extraction."""

    def test_extract_ref_accession_file_not_found(self):
        """Test reference accession extraction with non-existent file."""
        non_existent_path = Path("/tmp/non_existent.fasta")

        result = extract_ref_accession(non_existent_path)
        assert result == "NA"  # Based on the actual return value from the error


class TestReuseDetection:
    """Test reuse detection functionality."""

    def test_detect_reuse_from_plan_no_reuse(self):
        """Test reuse detection when no files are being reused."""
        mock_plan = MagicMock()
        mock_plan.paths.bam = None
        mock_plan.paths.vcf = None

        reuse_mapping, reuse_calling = detect_reuse_from_plan(mock_plan)
        assert reuse_mapping is False
        assert reuse_calling is False

    def test_detect_reuse_from_plan_bam_reuse(self):
        """Test reuse detection when BAM file exists."""
        with tempfile.TemporaryDirectory() as tmpdir:
            bam_path = Path(tmpdir) / "test.bam"
            bam_path.write_text("")  # Create empty file

            mock_plan = MagicMock()
            mock_plan.paths.bam = bam_path
            mock_plan.paths.vcf = None

            reuse_mapping, reuse_calling = detect_reuse_from_plan(mock_plan)
            assert reuse_mapping is True
            assert reuse_calling is False

    def test_detect_reuse_from_plan_vcf_reuse(self):
        """Test reuse detection when VCF file exists."""
        with tempfile.TemporaryDirectory() as tmpdir:
            vcf_path = Path(tmpdir) / "test.vcf"
            vcf_path.write_text("")  # Create empty file

            mock_plan = MagicMock()
            mock_plan.paths.bam = None
            mock_plan.paths.vcf = vcf_path

            reuse_mapping, reuse_calling = detect_reuse_from_plan(mock_plan)
            assert reuse_mapping is False
            assert reuse_calling is True


class TestRunPlanCreation:
    """Test run plan creation functionality."""

    def test_create_run_plan_self_mode_minimal(self):
        """Test creating run plan for self mode with minimal arguments."""
        with tempfile.TemporaryDirectory() as tmpdir:
            r1_path = Path(tmpdir) / "test_R1.fastq"
            r2_path = Path(tmpdir) / "test_R2.fastq"
            assembly_path = Path(tmpdir) / "assembly.fasta"

            # Create dummy files
            r1_path.write_text("@read1\nACGT\n+\nIIII\n")
            r2_path.write_text("@read2\nACGT\n+\nIIII\n")
            assembly_path.write_text(">contig1\nACGTACGT\n")

            with (
                patch(
                    "src.ssiamb.runner.infer_sample_name", return_value="test_sample"
                ),
                patch(
                    "src.ssiamb.runner.validate_sample_name", return_value="test_sample"
                ),
            ):

                plan = create_run_plan(
                    mode=Mode.SELF,
                    r1=r1_path,
                    r2=r2_path,
                    assembly=assembly_path,
                    sample=None,  # Will be inferred
                    output_dir=None,  # Will default to current
                    threads=4,
                    mapper="minimap2",
                    caller="bbtools",
                )

                assert plan.mode == Mode.SELF
                assert plan.sample == "test_sample"
                assert plan.mapper == Mapper.MINIMAP2
                assert plan.caller == Caller.BBTOOLS
                assert plan.threads == 4
                assert plan.paths.r1 == r1_path
                assert plan.paths.r2 == r2_path
                assert plan.paths.assembly == assembly_path

    def test_create_run_plan_ref_mode_defaults(self):
        """Test creating run plan for ref mode using defaults."""
        with tempfile.TemporaryDirectory() as tmpdir:
            r1_path = Path(tmpdir) / "sample_R1.fastq"
            r2_path = Path(tmpdir) / "sample_R2.fastq"
            reference_path = Path(tmpdir) / "reference.fasta"
            output_dir = Path(tmpdir) / "output"

            # Create dummy files and directories
            r1_path.write_text("@read1\nACGT\n+\nIIII\n")
            r2_path.write_text("@read2\nACGT\n+\nIIII\n")
            reference_path.write_text(">ref\nACGTACGT\n")
            output_dir.mkdir()

            plan = create_run_plan(
                mode=Mode.REF,
                r1=r1_path,
                r2=r2_path,
                reference=reference_path,
                sample="custom_sample",
                output_dir=output_dir,
                threads=8,
                mapper="bwa-mem2",
                caller="bcftools",
            )

            assert plan.mode == Mode.REF
            assert plan.sample == "custom_sample"
            assert plan.mapper == Mapper.BWA_MEM2
            assert plan.caller == Caller.BCFTOOLS
            assert plan.threads == 8
            # Use default thresholds since we didn't override them
            assert plan.thresholds.dp_min == 10  # Default value
            assert plan.thresholds.maf_min == 0.1  # Default value
            assert plan.thresholds.dp_cap == 100  # Default value
            assert plan.paths.output_dir == output_dir

    def test_create_run_plan_sample_validation_error(self):
        """Test run plan creation with invalid sample name."""
        with tempfile.TemporaryDirectory() as tmpdir:
            r1_path = Path(tmpdir) / "test_R1.fastq"
            r2_path = Path(tmpdir) / "test_R2.fastq"
            assembly_path = Path(tmpdir) / "assembly.fasta"

            r1_path.write_text("@read1\nACGT\n+\nIIII\n")
            r2_path.write_text("@read2\nACGT\n+\nIIII\n")
            assembly_path.write_text(">contig1\nACGTACGT\n")

            from src.ssiamb.io_utils import SampleNameError

            with patch(
                "src.ssiamb.runner.validate_sample_name",
                side_effect=SampleNameError("Invalid sample name"),
            ):
                with pytest.raises(SampleNameError):
                    create_run_plan(
                        mode=Mode.SELF,
                        r1=r1_path,
                        r2=r2_path,
                        assembly=assembly_path,
                        sample="invalid sample name",
                        threads=4,
                        mapper="minimap2",
                        caller="bbtools",
                    )
