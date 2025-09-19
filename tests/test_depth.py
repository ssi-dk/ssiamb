"""
Unit tests for depth analysis functionality.

Tests mosdepth execution, parsing, error handling, and contig filtering logic.
"""

import pytest
import tempfile
from pathlib import Path
from unittest.mock import Mock, patch

from src.ssiamb.depth import (
    DepthSummary,
    ContigDepthStats,
    DepthAnalysisError,
    check_mosdepth_available,
    run_mosdepth,
    parse_mosdepth_summary,
    analyze_depth,
)


class TestContigDepthStats:
    """Test ContigDepthStats dataclass functionality."""

    def test_contig_meets_length_threshold(self):
        """Test that contigs ≥500bp meet the threshold."""
        contig = ContigDepthStats(
            name="contig1",
            length=1000,
            bases_covered=950,
            mean_depth=25.0,
            breadth_10x=0.95,
        )
        assert contig.is_long_enough is True

    def test_contig_below_length_threshold(self):
        """Test that contigs <500bp don't meet the threshold."""
        contig = ContigDepthStats(
            name="contig2",
            length=499,
            bases_covered=450,
            mean_depth=20.0,
            breadth_10x=0.90,
        )
        assert contig.is_long_enough is False


class TestDepthSummary:
    """Test DepthSummary dataclass functionality."""

    def test_included_contig_names_property(self):
        """Test that included_contig_names returns contigs ≥500bp."""
        contig_stats = [
            ContigDepthStats("short1", 300, 280, 15.0, 0.93),  # Too short
            ContigDepthStats("long1", 1000, 950, 25.0, 0.95),  # Long enough
            ContigDepthStats("long2", 800, 720, 20.0, 0.90),  # Long enough
        ]

        summary = DepthSummary(
            callable_bases=1670,
            genome_length=1800,
            breadth_10x=0.925,
            mean_depth=22.5,
            total_contigs=3,
            included_contigs=2,
            contig_stats=contig_stats,
        )

        expected_names = {"long1", "long2"}
        assert summary.included_contig_names == expected_names


class TestMosdepthAvailability:
    """Test mosdepth availability checking."""

    @patch("subprocess.run")
    def test_mosdepth_available(self, mock_run):
        """Test successful mosdepth detection."""
        mock_run.return_value = Mock(returncode=0)
        assert check_mosdepth_available() is True

    @patch("subprocess.run")
    def test_mosdepth_not_available(self, mock_run):
        """Test mosdepth not available."""
        mock_run.side_effect = FileNotFoundError()
        assert check_mosdepth_available() is False


class TestRunMosdepth:
    """Test mosdepth execution functionality."""

    def test_bam_file_not_found(self):
        """Test error when BAM file doesn't exist."""
        with tempfile.TemporaryDirectory() as temp_dir:
            bam_path = Path(temp_dir) / "missing.bam"
            output_prefix = Path(temp_dir) / "test_output"

            with pytest.raises(FileNotFoundError, match="BAM file not found"):
                run_mosdepth(bam_path, output_prefix)

    @patch("src.ssiamb.depth.check_mosdepth_available")
    def test_mosdepth_not_available(self, mock_check):
        """Test error when mosdepth is not available."""
        with tempfile.TemporaryDirectory() as temp_dir:
            bam_path = Path(temp_dir) / "test.bam"
            bam_path.touch()
            output_prefix = Path(temp_dir) / "test_output"

            mock_check.return_value = False

            with pytest.raises(DepthAnalysisError, match="mosdepth not found in PATH"):
                run_mosdepth(bam_path, output_prefix)

    @patch("src.ssiamb.depth.check_mosdepth_available")
    @patch("subprocess.run")
    def test_successful_mosdepth_execution(self, mock_run, mock_check):
        """Test successful mosdepth execution."""
        with tempfile.TemporaryDirectory() as temp_dir:
            bam_path = Path(temp_dir) / "test.bam"
            bam_path.touch()
            output_prefix = Path(temp_dir) / "test_output"

            mock_check.return_value = True
            mock_run.return_value = Mock(returncode=0, stderr="")

            # Create expected output file
            expected_summary = output_prefix.with_suffix(".mosdepth.summary.txt")
            expected_summary.touch()

            result = run_mosdepth(bam_path, output_prefix)

            assert result == expected_summary
            mock_run.assert_called_once()

            # Check command construction
            call_args = mock_run.call_args[0][0]
            assert call_args[0] == "mosdepth"
            assert "--mapq" in call_args
            assert "--no-per-base" in call_args


class TestParseMosdepthSummary:
    """Test mosdepth summary parsing functionality."""

    def test_file_not_found(self):
        """Test error when summary file doesn't exist."""
        with tempfile.TemporaryDirectory() as temp_dir:
            non_existent = Path(temp_dir) / "missing.txt"

            with pytest.raises(
                FileNotFoundError, match="mosdepth summary file not found"
            ):
                parse_mosdepth_summary(non_existent)

    def test_successful_parsing_with_filtering(self):
        """Test successful parsing with contig length filtering."""
        with tempfile.TemporaryDirectory() as temp_dir:
            summary_file = Path(temp_dir) / "test.mosdepth.summary.txt"

            content = """chrom\tlength\tbases\tmean\tmin\tmax
short_contig\t300\t280\t15.5\t0\t50
long_contig1\t1000\t950\t25.0\t0\t100
long_contig2\t800\t720\t20.0\t0\t80
total\t2100\t1950\t20.5\t0\t100
"""
            summary_file.write_text(content)

            result = parse_mosdepth_summary(summary_file)

            # Check summary statistics
            assert result.total_contigs == 3
            assert result.included_contigs == 2  # Only contigs ≥500bp
            assert result.genome_length == 1800  # 1000 + 800
            assert result.callable_bases == 1670  # 950 + 720
            assert result.breadth_10x == pytest.approx(1670 / 1800, rel=1e-6)

            # Check included contig names
            expected_included = {"long_contig1", "long_contig2"}
            assert result.included_contig_names == expected_included

    def test_no_contigs_above_threshold(self):
        """Test parsing when no contigs meet length threshold."""
        with tempfile.TemporaryDirectory() as temp_dir:
            summary_file = Path(temp_dir) / "test.mosdepth.summary.txt"

            content = """chrom\tlength\tbases\tmean\tmin\tmax
small1\t300\t280\t15.5\t0\t50
small2\t400\t360\t12.0\t0\t45
total\t700\t640\t13.5\t0\t50
"""
            summary_file.write_text(content)

            result = parse_mosdepth_summary(summary_file)

            # Check that all metrics are zero when no contigs are long enough
            assert result.total_contigs == 2
            assert result.included_contigs == 0
            assert result.genome_length == 0
            assert result.callable_bases == 0
            assert result.breadth_10x == 0.0
            assert result.mean_depth == 0.0


class TestAnalyzeDepth:
    """Test complete depth analysis workflow."""

    @patch("src.ssiamb.depth.run_mosdepth")
    @patch("src.ssiamb.depth.parse_mosdepth_summary")
    def test_successful_workflow(self, mock_parse, mock_run):
        """Test successful complete depth analysis workflow."""
        with tempfile.TemporaryDirectory() as temp_dir:
            bam_path = Path(temp_dir) / "test.bam"
            bam_path.touch()

            # Mock mosdepth execution
            summary_file = Path(temp_dir) / "test.depth.mosdepth.summary.txt"
            mock_run.return_value = summary_file

            # Mock parsing results
            expected_summary = DepthSummary(
                callable_bases=1670,
                genome_length=1800,
                breadth_10x=0.925,
                mean_depth=22.5,
                total_contigs=4,
                included_contigs=2,
                contig_stats=[],
            )
            mock_parse.return_value = expected_summary

            result = analyze_depth(
                bam_path=bam_path, output_dir=Path(temp_dir), sample_name="test_sample"
            )

            assert result == expected_summary

            # Check that functions were called with correct parameters
            mock_run.assert_called_once_with(
                bam_path=bam_path,
                output_prefix=Path(temp_dir) / "test_sample.depth",
                mapq_threshold=30,
                threads=4,
            )
            mock_parse.assert_called_once_with(summary_file)


if __name__ == "__main__":
    pytest.main([__file__])
