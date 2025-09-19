"""
Integration tests for runner.py utility functions.

Tests utility functions like extract_ref_accession, version detection,
and other testable helper functions from the runner module.
"""

import tempfile
from pathlib import Path

from src.ssiamb.runner import (
    extract_ref_accession,
    _extract_version,
    _get_mapper_version,
    _get_caller_version,
)


class TestRefAccessionExtraction:
    """Test reference accession extraction from FASTA files."""

    def test_extract_ref_accession_success(self):
        """Test successful reference accession extraction."""
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        ) as tmp_file:
            tmp_file.write(
                ">NC_000913.3 Escherichia coli str. K-12 substr. MG1655, complete genome\n"
            )
            tmp_file.write("ACGTACGTACGT\n")
            tmp_file.flush()

            accession = extract_ref_accession(Path(tmp_file.name))
            assert accession == "NC_000913"  # Note: function extracts base accession

        # Clean up
        Path(tmp_file.name).unlink()

    def test_extract_ref_accession_genbank_format(self):
        """Test accession extraction from GenBank format header."""
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        ) as tmp_file:
            tmp_file.write(
                ">gi|49175990|ref|NC_000913.2| Escherichia coli str. K-12 substr. MG1655 chromosome, complete genome\n"
            )
            tmp_file.write("ACGTACGTACGT\n")
            tmp_file.flush()

            accession = extract_ref_accession(Path(tmp_file.name))
            assert accession == "gi|49175990|ref|NC_000913.2|"  # Returns full header

        # Clean up
        Path(tmp_file.name).unlink()

    def test_extract_ref_accession_simple_id(self):
        """Test accession extraction from simple identifier."""
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        ) as tmp_file:
            tmp_file.write(">CP000819.1\n")
            tmp_file.write("ACGTACGTACGT\n")
            tmp_file.flush()

            accession = extract_ref_accession(Path(tmp_file.name))
            assert accession == "CP000819.1"  # CP accessions keep version

        # Clean up
        Path(tmp_file.name).unlink()

    def test_extract_ref_accession_no_header(self):
        """Test reference accession extraction failure with no header."""
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        ) as tmp_file:
            tmp_file.write("ACGTACGTACGT\n")  # No header
            tmp_file.flush()

            accession = extract_ref_accession(Path(tmp_file.name))
            assert accession == "NA"  # Returns NA for invalid files

        # Clean up
        Path(tmp_file.name).unlink()

    def test_extract_ref_accession_malformed_header(self):
        """Test reference accession extraction with malformed header."""
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        ) as tmp_file:
            tmp_file.write(">malformed_header_no_accession\n")
            tmp_file.write("ACGTACGTACGT\n")
            tmp_file.flush()

            accession = extract_ref_accession(Path(tmp_file.name))
            assert accession == "malformed_header_no_accession"  # Returns header as-is

        # Clean up
        Path(tmp_file.name).unlink()

    def test_extract_ref_accession_empty_file(self):
        """Test reference accession extraction from empty file."""
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        ) as tmp_file:
            # Empty file
            tmp_file.flush()

            accession = extract_ref_accession(Path(tmp_file.name))
            assert accession == "NA"  # Returns NA for empty files

        # Clean up
        Path(tmp_file.name).unlink()

    def test_extract_ref_accession_missing_file(self):
        """Test reference accession extraction from missing file."""
        non_existent_file = Path("/non/existent/file.fasta")

        accession = extract_ref_accession(non_existent_file)
        assert accession == "NA"  # Returns NA for missing files


class TestVersionExtraction:
    """Test version extraction from tool output."""

    def test_extract_version_minimap2_format(self):
        """Test version extraction from minimap2-style output."""
        output = "2.24-r1122"
        version = _extract_version(output)
        assert version == "2.24"  # Extracts base version only

    def test_extract_version_bbtools_format(self):
        """Test version extraction from BBTools-style output."""
        output = "BBMap version 38.90"
        version = _extract_version(output)
        assert version == "38.90"

    def test_extract_version_bcftools_format(self):
        """Test version extraction from bcftools-style output."""
        output = "bcftools 1.15.1\nUsing htslib 1.15.1"
        version = _extract_version(output)
        assert version == "1.15.1"

    def test_extract_version_no_match(self):
        """Test version extraction when no version pattern is found."""
        output = "Some tool output without version information"
        version = _extract_version(output)
        assert version is None

    def test_extract_version_empty_output(self):
        """Test version extraction from empty output."""
        output = ""
        version = _extract_version(output)
        assert version is None

    def test_extract_version_complex_output(self):
        """Test version extraction from complex multi-line output."""
        output = """
        Tool Name: SomeMapper
        Version: 2.1.4-beta
        Build: 2023-01-15
        Options: --help
        """
        version = _extract_version(output)
        assert version == "2.1.4"  # Extracts base version only


class TestToolVersions:
    """Test tool version detection functions."""

    def test_get_mapper_version_minimap2(self):
        """Test getting minimap2 version command."""
        result = _get_mapper_version("minimap2")
        assert result == ["minimap2", "--version"]

    def test_get_mapper_version_bwa_mem2(self):
        """Test getting bwa-mem2 version command."""
        result = _get_mapper_version("bwa-mem2")
        assert result == ["bwa-mem2", "version"]

    def test_get_mapper_version_failure(self):
        """Test mapper version command for unknown mapper."""
        result = _get_mapper_version("nonexistent_mapper")
        assert result is None

    def test_get_caller_version_bbtools(self):
        """Test getting BBTools caller version command."""
        result = _get_caller_version("bbtools")
        assert result == ["bbmap.sh", "version"]

    def test_get_caller_version_bcftools(self):
        """Test getting bcftools version command."""
        result = _get_caller_version("bcftools")
        assert result == ["bcftools", "--version"]

    def test_get_caller_version_failure(self):
        """Test caller version command for unknown caller."""
        result = _get_caller_version("nonexistent_caller")
        assert result is None


class TestRunnerErrorHandling:
    """Test error handling in runner utility functions."""

    def test_extract_ref_accession_file_permission_error(self):
        """Test reference accession extraction with permission error."""
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        ) as tmp_file:
            tmp_file.write(">NC_000913.3 E. coli\n")
            tmp_file.write("ACGTACGT\n")
            tmp_file.flush()

            # Change permissions to make file unreadable
            import os

            os.chmod(tmp_file.name, 0o000)

            try:
                accession = extract_ref_accession(Path(tmp_file.name))
                # Should return "NA" on permission error
                assert accession == "NA"
            finally:
                # Restore permissions for cleanup
                os.chmod(tmp_file.name, 0o644)
                Path(tmp_file.name).unlink()

    def test_extract_ref_accession_with_spaces(self):
        """Test reference accession extraction with spaces in header."""
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        ) as tmp_file:
            tmp_file.write(">NC_000913.3 Escherichia coli K-12 MG1655\n")
            tmp_file.write("ACGTACGT\n")
            tmp_file.flush()

            accession = extract_ref_accession(Path(tmp_file.name))
            assert accession == "NC_000913"

        # Clean up
        Path(tmp_file.name).unlink()

    def test_version_extraction_edge_cases(self):
        """Test version extraction with edge cases."""
        # Test with just numbers
        assert _extract_version("1.0") == "1.0"

        # Test with alpha version
        assert _extract_version("tool version 2.1a") == "2.1"

        # Test with multiple version-like strings (should get first)
        assert _extract_version("old version 1.0, new version 2.0") == "1.0"

        # Test with no digits
        assert _extract_version("no version here") is None
