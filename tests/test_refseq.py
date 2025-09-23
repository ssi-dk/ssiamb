"""
Unit tests for RefSeq downloader functionality.

Tests NCBI API integration, genome selection policy,
download logic, and index generation.
"""

import tempfile
from pathlib import Path
from unittest.mock import patch, MagicMock

from src.ssiamb.refseq import (
    RefSeqDownloader,
    RefSeqGenome,
)


class TestRefSeqGenome:
    """Test RefSeqGenome data model."""

    def test_genome_creation(self):
        """Test basic genome object creation."""
        genome = RefSeqGenome(
            accession="GCF_000005825.2",
            organism_name="Escherichia coli str. K-12 substr. MG1655",
            category="reference",
            assembly_level="Complete Genome",
            assembly_type="RefSeq",
            release_date="2013/09/26",
            contig_count=1,
            total_length=4641652,
            ftp_path="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/825/GCF_000005825.2_ASM582v2",
        )

        assert genome.accession == "GCF_000005825.2"
        assert genome.organism_name == "Escherichia coli str. K-12 substr. MG1655"
        assert genome.assembly_level == "Complete Genome"
        assert genome.category == "reference"

    def test_genome_selection_priority(self):
        """Test genome selection priority scoring."""
        # Reference genome should have higher priority (lower score) than representative
        ref_genome = RefSeqGenome(
            accession="GCF_000005825.2",
            organism_name="Escherichia coli",
            category="reference",
            assembly_level="Complete Genome",
            assembly_type="RefSeq",
            release_date="2013/09/26",
            contig_count=1,
            total_length=4641652,
            ftp_path="https://ftp.example.com/ref",
        )

        rep_genome = RefSeqGenome(
            accession="GCF_000006745.1",
            organism_name="Escherichia coli",
            category="representative",
            assembly_level="Complete Genome",
            assembly_type="RefSeq",
            release_date="2013/09/26",
            contig_count=1,
            total_length=4641652,
            ftp_path="https://ftp.example.com/rep",
        )

        # Reference should have lower priority value (higher priority)
        assert ref_genome.get_selection_priority() < rep_genome.get_selection_priority()

    def test_genome_ordering(self):
        """Test genome sorting by selection priority."""
        genomes = [
            RefSeqGenome(
                "GCF_1",
                "E. coli",
                "na",
                "Scaffold",
                "RefSeq",
                "2020/01/01",
                100,
                4000000,
                "ftp1",
            ),
            RefSeqGenome(
                "GCF_2",
                "E. coli",
                "reference",
                "Complete Genome",
                "RefSeq",
                "2020/01/01",
                1,
                4641652,
                "ftp2",
            ),
            RefSeqGenome(
                "GCF_3",
                "E. coli",
                "representative",
                "Complete Genome",
                "RefSeq",
                "2020/01/01",
                1,
                4641652,
                "ftp3",
            ),
        ]

        sorted_genomes = sorted(genomes, key=lambda g: g.get_selection_priority())

        # Should be ordered: reference > representative > other
        assert sorted_genomes[0].category == "reference"
        assert sorted_genomes[1].category == "representative"
        assert sorted_genomes[2].category == "na"


class TestRefSeqDownloader:
    """Test RefSeqDownloader functionality."""

    def test_downloader_initialization(self):
        """Test downloader initialization."""
        with tempfile.TemporaryDirectory() as temp_dir:
            downloader = RefSeqDownloader(
                Path(temp_dir), verbose=True, create_indexes=False
            )
            assert downloader.output_dir == Path(temp_dir)
            assert downloader.verbose is True
            assert downloader.create_indexes is False
            assert downloader.session is not None

    def test_species_normalization(self):
        """Test species name normalization."""
        with tempfile.TemporaryDirectory() as temp_dir:
            downloader = RefSeqDownloader(Path(temp_dir))

            result = downloader.normalize_species_name("Escherichia coli")
            assert result == "Escherichia_coli"

            result = downloader.normalize_species_name(
                "Salmonella enterica subsp. enterica"
            )
            assert result == "Salmonella_enterica"

    @patch("src.ssiamb.refseq.requests.Session.get")
    def test_search_assemblies_success(self, mock_get):
        """Test successful assembly search."""
        # Mock NCBI API response
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "esearchresult": {"idlist": ["1234", "5678"], "count": "2"}
        }
        mock_get.return_value = mock_response

        with tempfile.TemporaryDirectory() as temp_dir:
            downloader = RefSeqDownloader(Path(temp_dir))

            result = downloader.search_assemblies("Escherichia coli")
            assert result == ["1234", "5678"]

            # Verify API call was made
            mock_get.assert_called_once()
            call_args = mock_get.call_args
            assert "esearch.fcgi" in call_args[0][0]

    @patch("src.ssiamb.refseq.requests.Session.get")
    def test_search_assemblies_no_results(self, mock_get):
        """Test assembly search with no results."""
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "esearchresult": {"idlist": [], "count": "0"}
        }
        mock_get.return_value = mock_response

        with tempfile.TemporaryDirectory() as temp_dir:
            downloader = RefSeqDownloader(Path(temp_dir))

            result = downloader.search_assemblies("Nonexistent species")
            assert result == []

    def test_index_generation_missing_file(self):
        """Test index generation with missing FASTA file."""
        with tempfile.TemporaryDirectory() as temp_dir:
            downloader = RefSeqDownloader(Path(temp_dir))
            nonexistent_file = Path(temp_dir) / "nonexistent.fna"

            result = downloader.generate_indexes(
                nonexistent_file, "test_species", dry_run=False
            )
            assert result is False

    def test_index_generation_dry_run(self):
        """Test index generation in dry run mode."""
        with tempfile.TemporaryDirectory() as temp_dir:
            downloader = RefSeqDownloader(Path(temp_dir))
            fasta_file = Path(temp_dir) / "test.fna"

            result = downloader.generate_indexes(
                fasta_file, "test_species", dry_run=True
            )
            assert result is True  # Should always succeed in dry run

    @patch("subprocess.run")
    def test_generate_minimap2_index(self, mock_run):
        """Test minimap2 index generation."""
        mock_run.return_value = MagicMock(returncode=0)

        with tempfile.TemporaryDirectory() as temp_dir:
            # Create a mock FASTA file
            fasta_file = Path(temp_dir) / "test.fna"
            fasta_file.write_text(">test\\nACGT\\n")

            downloader = RefSeqDownloader(Path(temp_dir))
            result = downloader._generate_minimap2_index(fasta_file, "test_species")

            assert result is True
            mock_run.assert_called_once()

            # Check command
            call_args = mock_run.call_args[0][0]
            assert call_args[0] == "minimap2"
            assert "-d" in call_args

    @patch("subprocess.run")
    def test_generate_bwa_mem2_index(self, mock_run):
        """Test bwa-mem2 index generation."""
        mock_run.return_value = MagicMock(returncode=0)

        with tempfile.TemporaryDirectory() as temp_dir:
            # Create a mock FASTA file
            fasta_file = Path(temp_dir) / "test.fna"
            fasta_file.write_text(">test\\nACGT\\n")

            # Create expected index files
            expected_suffixes = [".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"]
            for suffix in expected_suffixes:
                index_file = fasta_file.with_suffix(f"{fasta_file.suffix}{suffix}")
                index_file.write_text("mock index content")

            downloader = RefSeqDownloader(Path(temp_dir))
            result = downloader._generate_bwa_mem2_index(fasta_file, "test_species")

            assert result is True
            mock_run.assert_called_once()

            # Check command
            call_args = mock_run.call_args[0][0]
            assert call_args[0] == "bwa-mem2"
            assert "index" in call_args

    @patch("subprocess.run")
    def test_index_generation_failure(self, mock_run):
        """Test handling of index generation failure."""
        mock_run.return_value = MagicMock(returncode=1, stderr="Error")

        with tempfile.TemporaryDirectory() as temp_dir:
            fasta_file = Path(temp_dir) / "test.fna"
            fasta_file.write_text(">test\\nACGT\\n")

            downloader = RefSeqDownloader(Path(temp_dir))
            result = downloader._generate_minimap2_index(fasta_file, "test_species")

            assert result is False
