"""
Unit tests for reuse and compatibility checking module.

Tests BAM/VCF compatibility checks, duplicate detection,
and compatibility thresholds per spec requirements.
"""

import pytest
import tempfile
import subprocess
from pathlib import Path
from unittest.mock import patch, MagicMock

from src.ssiamb.reuse import (
    CompatibilityError,
    get_fasta_contigs,
    get_bam_contigs,
    get_vcf_contigs,
    check_compatibility,
    has_duplicate_flags,
    run_markdup_for_depth,
)


class TestCompatibilityChecking:
    """Test compatibility checking functionality."""

    def test_check_compatibility_basic_pass(self):
        """Test basic compatibility check that passes all thresholds."""
        fasta_contigs = {"chr1": 1000, "chr2": 2000, "chr3": 500}

        bam_contigs = {
            "chr1": 1005,  # 0.5% diff, within 1% threshold
            "chr2": 1990,  # 0.5% diff, within 1% threshold
            "chr3": 495,  # 1.0% diff, exactly at threshold
        }

        result = check_compatibility(fasta_contigs, bam_contigs)

        assert result.is_compatible is True
        assert result.coverage_fraction >= 0.95  # Should be 100%
        assert result.total_length_diff <= 0.02
        assert result.error_message is None

    def test_check_compatibility_coverage_fail(self):
        """Test compatibility check that fails coverage threshold."""
        fasta_contigs = {
            "chr1": 1000,
            "chr2": 2000,
            "chr3": 3000,  # Large contig missing from BAM
        }

        bam_contigs = {
            "chr1": 1000,
            "chr2": 2000,
            # chr3 missing - only 50% coverage
        }

        result = check_compatibility(fasta_contigs, bam_contigs)

        assert result.is_compatible is False
        assert result.coverage_fraction < 0.95
        assert result.error_message is not None
        assert "Coverage" in result.error_message
        assert "chr3" in result.missing_contigs

    def test_check_compatibility_per_contig_fail(self):
        """Test compatibility check that fails per-contig threshold."""
        fasta_contigs = {"chr1": 1000, "chr2": 2000}

        bam_contigs = {
            "chr1": 1000,
            "chr2": 1800,  # 10% difference, exceeds 1% threshold
        }

        result = check_compatibility(fasta_contigs, bam_contigs)

        assert result.is_compatible is False
        assert result.error_message is not None
        assert "Per-contig length violations" in result.error_message
        assert "chr2" in result.error_message

    def test_check_compatibility_total_length_fail(self):
        """Test compatibility check that fails total length threshold."""
        fasta_contigs = {"chr1": 1000, "chr2": 2000, "chr3": 3000}

        bam_contigs = {
            "chr1": 1000,
            "chr2": 2000,
            "chr3": 2800,  # 6.7% difference in chr3, total diff > 2%
        }

        result = check_compatibility(fasta_contigs, bam_contigs)

        assert result.is_compatible is False
        assert result.error_message is not None
        assert "Total length diff" in result.error_message

    def test_check_compatibility_edge_cases(self):
        """Test compatibility edge cases."""
        # Empty FASTA
        result = check_compatibility({}, {"chr1": 1000})
        assert result.coverage_fraction == 0.0
        assert result.is_compatible is False

        # Empty BAM
        result = check_compatibility({"chr1": 1000}, {})
        assert result.coverage_fraction == 0.0
        assert result.is_compatible is False

        # Perfect match
        contigs = {"chr1": 1000, "chr2": 2000}
        result = check_compatibility(contigs, contigs)
        assert result.is_compatible is True
        assert result.coverage_fraction == 1.0
        assert result.total_length_diff == 0.0


class TestContigExtraction:
    """Test contig extraction from various file formats."""

    def test_get_fasta_contigs_basic(self):
        """Test FASTA contig extraction."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            # Create test FASTA
            fasta_content = (
                ">chr1\n"
                + "A" * 1000
                + "\n>chr2\n"
                + "C" * 2000
                + "\n>chr3\n"
                + "G" * 300
                + "\n"
            )
            fasta_file = tmpdir / "test.fa"
            fasta_file.write_text(fasta_content)

            # Index the FASTA for pysam
            with patch("pysam.FastaFile") as mock_fasta:
                mock_fasta.return_value.__enter__.return_value.references = [
                    "chr1",
                    "chr2",
                    "chr3",
                ]
                mock_fasta.return_value.__enter__.return_value.get_reference_length.side_effect = [
                    1000,
                    2000,
                    300,
                ]

                contigs = get_fasta_contigs(fasta_file, min_length=500)

                # Only chr1 and chr2 should be included (≥500bp)
                assert len(contigs) == 2
                assert contigs["chr1"] == 1000
                assert contigs["chr2"] == 2000
                assert "chr3" not in contigs  # Too short

    def test_get_bam_contigs_basic(self):
        """Test BAM contig extraction."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            bam_file = tmpdir / "test.bam"

            with patch("pysam.AlignmentFile") as mock_bam:
                mock_bam.return_value.__enter__.return_value.references = [
                    "chr1",
                    "chr2",
                ]
                mock_bam.return_value.__enter__.return_value.lengths = [1000, 2000]

                contigs = get_bam_contigs(bam_file)

                assert len(contigs) == 2
                assert contigs["chr1"] == 1000
                assert contigs["chr2"] == 2000

    def test_get_vcf_contigs_basic(self):
        """Test VCF contig extraction."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            vcf_file = tmpdir / "test.vcf"

            # Mock VCF header contigs
            mock_contig1 = MagicMock()
            mock_contig1.length = 1000
            mock_contig2 = MagicMock()
            mock_contig2.length = 2000
            mock_contig3 = MagicMock()
            mock_contig3.length = None  # Missing length

            with patch("pysam.VariantFile") as mock_vcf:
                mock_vcf.return_value.__enter__.return_value.header.contigs.items.return_value = [
                    ("chr1", mock_contig1),
                    ("chr2", mock_contig2),
                    ("chr3", mock_contig3),
                ]

                contigs = get_vcf_contigs(vcf_file)

                # Only chr1 and chr2 should be included (have length)
                assert len(contigs) == 2
                assert contigs["chr1"] == 1000
                assert contigs["chr2"] == 2000
                assert "chr3" not in contigs  # No length


class TestDuplicateDetection:
    """Test duplicate flag detection functionality."""

    def test_has_duplicate_flags_with_duplicates(self):
        """Test duplicate detection when duplicates are present."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            bam_file = tmpdir / "test.bam"

            # Mock reads with some duplicates
            mock_read1 = MagicMock()
            mock_read1.is_duplicate = False
            mock_read2 = MagicMock()
            mock_read2.is_duplicate = True
            mock_read3 = MagicMock()
            mock_read3.is_duplicate = False

            with patch("pysam.AlignmentFile") as mock_bam:
                mock_bam.return_value.__enter__.return_value.fetch.return_value = [
                    mock_read1,
                    mock_read2,
                    mock_read3,
                ]

                result = has_duplicate_flags(bam_file)

                assert result is True

    def test_has_duplicate_flags_without_duplicates(self):
        """Test duplicate detection when no duplicates are present."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            bam_file = tmpdir / "test.bam"

            # Mock reads without duplicates
            mock_reads = [MagicMock() for _ in range(5)]
            for read in mock_reads:
                read.is_duplicate = False

            with patch("pysam.AlignmentFile") as mock_bam:
                mock_bam.return_value.__enter__.return_value.fetch.return_value = (
                    mock_reads
                )

                result = has_duplicate_flags(bam_file)

                assert result is False

    def test_has_duplicate_flags_empty_bam(self):
        """Test duplicate detection with empty BAM file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            bam_file = tmpdir / "test.bam"

            with patch("pysam.AlignmentFile") as mock_bam:
                mock_bam.return_value.__enter__.return_value.fetch.return_value = []

                result = has_duplicate_flags(bam_file)

                assert result is False


class TestMarkdupFunctionality:
    """Test samtools markdup functionality."""

    @patch("shutil.which")
    @patch("subprocess.run")
    def test_run_markdup_for_depth_success(self, mock_run, mock_which):
        """Test successful markdup operation."""
        mock_which.return_value = "/usr/bin/samtools"
        mock_run.return_value = MagicMock(returncode=0)

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            input_bam = tmpdir / "input.bam"
            input_bam.touch()  # Create dummy file

            result = run_markdup_for_depth(input_bam, tmpdir)

            assert str(result).endswith(".markdup.bam")
            assert (
                mock_run.call_count == 5
            )  # queryname sort + fixmate + coord sort + markdup + index

    @patch("shutil.which")
    def test_run_markdup_for_depth_missing_samtools(self, mock_which):
        """Test markdup failure when samtools is missing."""
        mock_which.return_value = None

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            input_bam = tmpdir / "input.bam"

            with pytest.raises(CompatibilityError, match="samtools not found"):
                run_markdup_for_depth(input_bam, tmpdir)

    @patch("shutil.which")
    @patch("subprocess.run")
    def test_run_markdup_for_depth_command_failure(self, mock_run, mock_which):
        """Test markdup failure when samtools command fails."""
        mock_which.return_value = "/usr/bin/samtools"
        mock_run.side_effect = subprocess.CalledProcessError(
            1, "samtools", stderr="Error"
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            input_bam = tmpdir / "input.bam"

            with pytest.raises(CompatibilityError, match="samtools markdup failed"):
                run_markdup_for_depth(input_bam, tmpdir)


class TestErrorHandling:
    """Test error handling in reuse module."""

    def test_compatibility_error_inheritance(self):
        """Test that CompatibilityError is properly defined."""
        error = CompatibilityError("test message")
        assert isinstance(error, Exception)
        assert str(error) == "test message"

    @patch("pysam.FastaFile")
    def test_get_fasta_contigs_error_handling(self, mock_fasta):
        """Test error handling in FASTA contig extraction."""
        mock_fasta.side_effect = Exception("File not found")

        with pytest.raises(CompatibilityError, match="Failed to read FASTA file"):
            get_fasta_contigs(Path("nonexistent.fa"))

    @patch("pysam.AlignmentFile")
    def test_get_bam_contigs_error_handling(self, mock_bam):
        """Test error handling in BAM contig extraction."""
        mock_bam.side_effect = Exception("File not found")

        with pytest.raises(CompatibilityError, match="Failed to read BAM file"):
            get_bam_contigs(Path("nonexistent.bam"))

    @patch("pysam.VariantFile")
    def test_get_vcf_contigs_error_handling(self, mock_vcf):
        """Test error handling in VCF contig extraction."""
        mock_vcf.side_effect = Exception("File not found")

        with pytest.raises(CompatibilityError, match="Failed to read VCF file"):
            get_vcf_contigs(Path("nonexistent.vcf"))


if __name__ == "__main__":
    pytest.main([__file__])
