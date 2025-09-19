"""
Tests for mapping.py functions.

Tests the core mapping functionality including tool detection, index management,
and mapping rate calculation that are most critical for coverage.
"""

import pytest
import tempfile
import subprocess
from pathlib import Path
from unittest.mock import patch

from src.ssiamb.mapping import (
    get_index_files,
    indexes_exist,
    calculate_mapping_rate,
    MappingError,
)
from src.ssiamb.models import Mapper


class TestIndexFiles:
    """Test index file path generation and detection."""

    def test_get_index_files_minimap2(self):
        """Test index file paths for minimap2."""
        fasta_path = Path("/tmp/reference.fasta")

        index_files = get_index_files(fasta_path, Mapper.MINIMAP2)

        assert len(index_files) == 1
        assert index_files[0] == Path("/tmp/reference.mmi")

    def test_get_index_files_bwa_mem2(self):
        """Test index file paths for bwa-mem2."""
        fasta_path = Path("/tmp/reference.fasta")

        index_files = get_index_files(fasta_path, Mapper.BWA_MEM2)

        expected_files = [
            Path("/tmp/reference.fasta.0123"),
            Path("/tmp/reference.fasta.amb"),
            Path("/tmp/reference.fasta.ann"),
            Path("/tmp/reference.fasta.pac"),
            Path("/tmp/reference.fasta.bwt.2bit.64"),
        ]

        assert len(index_files) == 5
        assert set(index_files) == set(expected_files)

    def test_get_index_files_fasta_extension_variants(self):
        """Test index file generation with different FASTA extensions."""
        for ext in [".fasta", ".fa", ".fna"]:
            fasta_path = Path(f"/tmp/reference{ext}")

            index_files = get_index_files(fasta_path, Mapper.MINIMAP2)

            assert index_files[0] == Path("/tmp/reference.mmi")

    def test_indexes_exist_all_present(self):
        """Test index existence check when all files are present."""
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = Path(tmpdir) / "reference.fasta"
            fasta_path.touch()

            # Create minimap2 index
            index_file = Path(tmpdir) / "reference.mmi"
            index_file.touch()

            assert indexes_exist(fasta_path, Mapper.MINIMAP2) is True

    def test_indexes_exist_some_missing(self):
        """Test index existence check when some files are missing."""
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = Path(tmpdir) / "reference.fasta"
            fasta_path.touch()

            # Create only some bwa-mem2 index files
            (Path(tmpdir) / "reference.0123").touch()
            (Path(tmpdir) / "reference.amb").touch()
            # Missing .ann, .pac, .bwt.2bit.64

            assert indexes_exist(fasta_path, Mapper.BWA_MEM2) is False

    def test_indexes_exist_none_present(self):
        """Test index existence check when no files are present."""
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = Path(tmpdir) / "reference.fasta"
            fasta_path.touch()

            assert indexes_exist(fasta_path, Mapper.MINIMAP2) is False
            assert indexes_exist(fasta_path, Mapper.BWA_MEM2) is False


class TestMappingRate:
    """Test mapping rate calculation functionality."""

    @patch("src.ssiamb.mapping.subprocess.run")
    @patch("src.ssiamb.mapping.shutil.which")
    def test_calculate_mapping_rate_success(self, mock_which, mock_run):
        """Test successful mapping rate calculation."""
        mock_which.return_value = "/usr/bin/samtools"

        # Mock samtools stats output with correct format
        mock_run.return_value.returncode = 0
        mock_run.return_value.stdout = """
SN	raw total sequences:	1000
SN	reads mapped:	850
SN	reads unmapped:	150
"""

        with tempfile.NamedTemporaryFile(suffix=".bam") as tmp_bam:
            bam_path = Path(tmp_bam.name)

            mapping_rate = calculate_mapping_rate(bam_path)

            assert mapping_rate == 0.85  # 850/1000
            mock_run.assert_called_once()

    @patch("src.ssiamb.mapping.subprocess.run")
    @patch("src.ssiamb.mapping.shutil.which")
    def test_calculate_mapping_rate_zero_reads(self, mock_which, mock_run):
        """Test mapping rate calculation with zero reads."""
        mock_which.return_value = "/usr/bin/samtools"

        # Mock samtools stats output with zero reads
        mock_run.return_value.returncode = 0
        mock_run.return_value.stdout = """
SN	raw total sequences:	0
SN	reads mapped:	0
SN	reads unmapped:	0
"""

        with tempfile.NamedTemporaryFile(suffix=".bam") as tmp_bam:
            bam_path = Path(tmp_bam.name)

            mapping_rate = calculate_mapping_rate(bam_path)

            assert mapping_rate == 0.0

    @patch("src.ssiamb.mapping.subprocess.run")
    @patch("src.ssiamb.mapping.shutil.which")
    def test_calculate_mapping_rate_samtools_error(self, mock_which, mock_run):
        """Test mapping rate calculation when samtools fails."""
        mock_which.return_value = "/usr/bin/samtools"
        # Mock subprocess.CalledProcessError for non-zero return code
        mock_run.side_effect = subprocess.CalledProcessError(
            1, ["samtools", "stats"], stderr="samtools: error reading file"
        )

        with tempfile.NamedTemporaryFile(suffix=".bam") as tmp_bam:
            bam_path = Path(tmp_bam.name)

            with pytest.raises(MappingError, match="samtools stats failed"):
                calculate_mapping_rate(bam_path)


class TestInputValidation:
    """Test input validation for mapping functions."""

    def test_get_index_files_different_extensions(self):
        """Test index file generation handles various FASTA extensions correctly."""
        extensions = [".fasta", ".fa", ".fna", ".fas"]

        for ext in extensions:
            fasta_path = Path(f"/data/reference{ext}")

            # Test minimap2
            index_files = get_index_files(fasta_path, Mapper.MINIMAP2)
            assert len(index_files) == 1
            assert index_files[0].suffix == ".mmi"
            assert index_files[0].stem == "reference"

            # Test bwa-mem2
            index_files = get_index_files(fasta_path, Mapper.BWA_MEM2)
            assert len(index_files) == 5
            # All should have the same stem (reference)
            assert all(f.stem.startswith("reference") for f in index_files)

    def test_indexes_exist_edge_cases(self):
        """Test index existence checking with edge cases."""
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = Path(tmpdir) / "test.fasta"
            fasta_path.touch()

            # Test with empty directory (no index files)
            assert indexes_exist(fasta_path, Mapper.MINIMAP2) is False
            assert indexes_exist(fasta_path, Mapper.BWA_MEM2) is False

            # Test with partial index files for bwa-mem2
            partial_files = [".fasta.0123", ".fasta.amb"]  # Missing 3 files
            for ext in partial_files:
                (Path(tmpdir) / f"test{ext}").touch()

            assert (
                indexes_exist(fasta_path, Mapper.BWA_MEM2) is False
            )  # Not all present

            # Complete the index files
            remaining_files = [".fasta.ann", ".fasta.pac", ".fasta.bwt.2bit.64"]
            for ext in remaining_files:
                (Path(tmpdir) / f"test{ext}").touch()

            assert indexes_exist(fasta_path, Mapper.BWA_MEM2) is True  # All present now


class TestMappingUtilities:
    """Test utility functions for mapping operations."""

    @patch("src.ssiamb.mapping.subprocess.run")
    @patch("src.ssiamb.mapping.shutil.which")
    def test_calculate_mapping_rate_complex_output(self, mock_which, mock_run):
        """Test mapping rate calculation with complex samtools output."""
        mock_which.return_value = "/usr/bin/samtools"

        # Mock realistic samtools stats output
        mock_run.return_value.returncode = 0
        mock_run.return_value.stdout = """
CHK	PSC	0d41d8cd98f00b204e9800998ecf8427e	.
SN	raw total sequences:	10000
SN	filtered sequences:	0
SN	sequences:	10000
SN	is sorted:	1
SN	1st fragments:	5000
SN	last fragments:	5000
SN	reads mapped:	9234
SN	reads mapped and paired:	9000
SN	reads unmapped:	766
SN	reads properly paired:	8800
SN	reads paired:	10000
SN	reads duplicated:	45
SN	reads MQ0:	234
SN	reads QC failed:	0
SN	non-primary alignments:	0
SN	total length:	1500000
SN	total first fragment length:	750000
SN	total last fragment length:	750000
SN	bases mapped:	1385100
SN	bases mapped (cigar):	1385100
SN	bases trimmed:	0
SN	bases duplicated:	6750
SN	mismatches:	2356
SN	error rate:	1.700450e-03
SN	average length:	150
SN	average first fragment length:	150
SN	average last fragment length:	150
SN	maximum length:	150
SN	maximum first fragment length:	150
SN	maximum last fragment length:	150
SN	average quality:	37.0
SN	insert size average:	450.5
SN	insert size standard deviation:	123.4
SN	inward oriented pairs:	4400
SN	outward oriented pairs:	0
SN	pairs with other orientation:	0
SN	pairs on different chromosomes:	0
SN	percentage of properly paired reads (%):	88.0
"""

        with tempfile.NamedTemporaryFile(suffix=".bam") as tmp_bam:
            bam_path = Path(tmp_bam.name)

            mapping_rate = calculate_mapping_rate(bam_path)

            # 9234 mapped out of 10000 total
            assert mapping_rate == 0.9234
            mock_run.assert_called_once()

    @patch("src.ssiamb.mapping.subprocess.run")
    @patch("src.ssiamb.mapping.shutil.which")
    def test_calculate_mapping_rate_malformed_output(self, mock_which, mock_run):
        """Test mapping rate calculation with malformed samtools output."""
        mock_which.return_value = "/usr/bin/samtools"

        # Mock malformed output (missing expected fields)
        mock_run.return_value.returncode = 0
        mock_run.return_value.stdout = """
SN	some other field:	1000
SN	different field:	850
CHK	checksum:	abc123
"""

        with tempfile.NamedTemporaryFile(suffix=".bam") as tmp_bam:
            bam_path = Path(tmp_bam.name)

            mapping_rate = calculate_mapping_rate(bam_path)

            # Should return 0.0 when no expected fields are found
            assert mapping_rate == 0.0
