"""
Unit tests for the calling module.

Tests BBTools and bcftools variant calling functionality with mocked subprocess calls.
"""

import pytest
from unittest.mock import Mock, patch

from src.ssiamb.calling import (
    check_caller_tools,
    run_bbtools_calling,
    call_variants,
    get_available_callers,
    VariantCallingError,
    VariantCallResult,
)
from src.ssiamb.models import Caller


class TestCallerToolChecks:
    """Test tool availability checking."""

    @patch("src.ssiamb.calling.shutil.which")
    def test_bbtools_tools_available(self, mock_which):
        """Test BBTools tool availability check when tools are present."""
        mock_which.side_effect = lambda tool: (
            f"/usr/bin/{tool}" if tool in ["pileup.sh", "callvariants.sh"] else None
        )

        assert check_caller_tools(Caller.BBTOOLS) is True
        assert mock_which.call_count == 2

    @patch("src.ssiamb.calling.shutil.which")
    def test_bbtools_tools_missing(self, mock_which):
        """Test BBTools tool availability check when tools are missing."""
        mock_which.return_value = None

        assert check_caller_tools(Caller.BBTOOLS) is False

    @patch("src.ssiamb.calling.shutil.which")
    def test_bcftools_available(self, mock_which):
        """Test bcftools tool availability check when tool is present."""
        mock_which.side_effect = lambda tool: (
            "/usr/bin/bcftools" if tool == "bcftools" else None
        )

        assert check_caller_tools(Caller.BCFTOOLS) is True
        mock_which.assert_called_once_with("bcftools")

    @patch("src.ssiamb.calling.shutil.which")
    def test_bcftools_missing(self, mock_which):
        """Test bcftools tool availability check when tool is missing."""
        mock_which.return_value = None

        assert check_caller_tools(Caller.BCFTOOLS) is False


class TestBBToolsCalling:
    """Test BBTools variant calling functionality."""

    @patch("src.ssiamb.calling.subprocess.run")
    def test_successful_bbtools_calling(self, mock_run, tmp_path):
        """Test successful BBTools variant calling pipeline."""
        # Setup paths
        bam_path = tmp_path / "test.bam"
        ref_path = tmp_path / "ref.fasta"
        output_vcf = tmp_path / "output.vcf"

        # Create real files so Path.exists() and Path.stat() work properly
        bam_path.write_bytes(b"dummy bam content")
        ref_path.write_text("dummy fasta content")

        # Mock successful subprocess call - also create the output file
        def mock_subprocess(*args, **kwargs):
            # Create the output VCF file as if the tool did it
            output_vcf.write_text(
                "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            )
            return Mock(returncode=0, stderr="", stdout="")

        mock_run.side_effect = mock_subprocess

        # Run BBTools calling
        result = run_bbtools_calling(
            bam_path=bam_path,
            reference_path=ref_path,
            output_vcf=output_vcf,
            sample_name="test_sample",
            threads=2,
            mapq_min=25,
            baseq_min=25,
        )

        # Check result
        assert result.success is True
        assert result.caller == Caller.BBTOOLS
        assert result.vcf_path == output_vcf
        assert result.error_message is None
        assert result.runtime_seconds is not None

        # Check that callvariants was called once (new implementation)
        assert mock_run.call_count == 1

    @patch("src.ssiamb.calling.subprocess.run")
    def test_bbtools_pileup_failure(self, mock_run, tmp_path):
        """Test BBTools calling when pileup.sh fails."""
        # Setup paths
        bam_path = tmp_path / "test.bam"
        ref_path = tmp_path / "ref.fasta"
        output_vcf = tmp_path / "output.vcf"

        # Mock failed pileup
        mock_run.return_value = Mock(returncode=1, stderr="Pileup failed")

        # Run BBTools calling
        result = run_bbtools_calling(
            bam_path=bam_path,
            reference_path=ref_path,
            output_vcf=output_vcf,
            sample_name="test_sample",
        )

        # Check result indicates failure
        assert result.success is False
        assert result.error_message is not None
        assert "BBTools callvariants failed" in result.error_message


class TestCallVariants:
    """Test the main call_variants function."""

    @patch("src.ssiamb.calling.check_caller_tools")
    @patch("src.ssiamb.calling.run_bbtools_calling")
    @patch("src.ssiamb.calling.Path.exists")
    def test_call_variants_bbtools(
        self, mock_exists, mock_bbtools, mock_check_tools, tmp_path
    ):
        """Test call_variants with BBTools caller."""
        # Setup
        bam_path = tmp_path / "test.bam"
        ref_path = tmp_path / "ref.fasta"
        output_vcf = tmp_path / "output.vcf"

        mock_exists.return_value = True
        mock_check_tools.return_value = True
        mock_bbtools.return_value = VariantCallResult(
            vcf_path=output_vcf, caller=Caller.BBTOOLS, success=True
        )

        # Call function
        result = call_variants(
            bam_path=bam_path,
            reference_path=ref_path,
            output_vcf=output_vcf,
            caller=Caller.BBTOOLS,
            sample_name="test_sample",
        )

        # Check result
        assert result.success is True
        assert result.caller == Caller.BBTOOLS
        mock_bbtools.assert_called_once()

    @patch("src.ssiamb.calling.caller_tools_available")
    def test_call_variants_tools_unavailable(self, mock_tools_available, tmp_path):
        """Test call_variants when caller tools are not available."""
        # Setup - create dummy files
        bam_path = tmp_path / "test.bam"
        ref_path = tmp_path / "ref.fasta"
        output_vcf = tmp_path / "output.vcf"

        # Create dummy files
        bam_path.write_bytes(b"dummy bam")
        ref_path.write_text("dummy fasta")

        mock_tools_available.return_value = False

        # Should raise VariantCallingError
        with pytest.raises(VariantCallingError, match="Required tools.*not available"):
            call_variants(
                bam_path=bam_path,
                reference_path=ref_path,
                output_vcf=output_vcf,
                caller=Caller.BBTOOLS,
                sample_name="test_sample",
            )


class TestGetAvailableCallers:
    """Test the get_available_callers function."""

    @patch("src.ssiamb.calling.caller_tools_available")
    def test_get_available_callers_all_available(self, mock_tools_available):
        """Test when all callers are available."""
        mock_tools_available.return_value = True

        available = get_available_callers()

        assert len(available) == 2
        assert Caller.BBTOOLS in available
        assert Caller.BCFTOOLS in available

    @patch("src.ssiamb.calling.caller_tools_available")
    def test_get_available_callers_none(self, mock_tools_available):
        """Test when no callers are available."""
        mock_tools_available.return_value = False

        available = get_available_callers()

        assert len(available) == 0
