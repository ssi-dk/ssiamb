"""
Essential tests for Bracken parsing and species selection.

Focused on core functionality and critical edge cases.
"""

import pytest
import tempfile
from pathlib import Path

from src.ssiamb.bracken import (
    BrackenSpecies,
    BrackenThresholds,
    BrackenSelectionError,
    parse_bracken_file,
    select_species,
    select_species_from_file,
    normalize_species_name_for_resolution,
)


class TestBrackenParsing:
    """Test core Bracken file parsing."""

    def create_test_bracken_file(self, content: str) -> Path:
        """Helper to create temporary Bracken file."""
        temp_file = tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False)
        temp_file.write(content)
        temp_file.close()
        return Path(temp_file.name)

    def test_parse_valid_file(self):
        """Test parsing a valid Bracken file."""
        content = """name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction_total_reads
Listeria monocytogenes\t1639\tS\t1718713\t145289\t1864002\t0.99520
Listeria innocua\t1642\tS\t5984\t764\t6748\t0.00360
Escherichia coli\t562\tG\t1000\t500\t1500\t0.00080"""

        file_path = self.create_test_bracken_file(content)

        try:
            species_list = parse_bracken_file(file_path)

            # Should only get species-level entries (taxonomy_lvl == "S")
            assert len(species_list) == 2
            assert species_list[0].name == "Listeria monocytogenes"
            assert species_list[0].new_est_reads == 1864002
            assert species_list[0].fraction_total_reads == 0.99520

        finally:
            file_path.unlink()

    def test_parse_file_errors(self):
        """Test error handling for invalid files."""
        # Missing file
        with pytest.raises(FileNotFoundError):
            parse_bracken_file(Path("nonexistent.tsv"))

        # Empty file
        empty_file = self.create_test_bracken_file("")
        try:
            with pytest.raises(ValueError, match="empty"):
                parse_bracken_file(empty_file)
        finally:
            empty_file.unlink()


class TestSpeciesSelection:
    """Test species selection algorithm."""

    def create_test_species(self) -> list[BrackenSpecies]:
        """Create test species data."""
        return [
            BrackenSpecies(
                name="Listeria monocytogenes",
                taxonomy_id=1639,
                taxonomy_lvl="S",
                kraken_assigned_reads=1718713,
                added_reads=145289,
                new_est_reads=1864002,
                fraction_total_reads=0.99520,
            ),
            BrackenSpecies(
                name="Escherichia coli",
                taxonomy_id=562,
                taxonomy_lvl="S",
                kraken_assigned_reads=45000,
                added_reads=55000,
                new_est_reads=100000,
                fraction_total_reads=0.70000,
            ),
            BrackenSpecies(
                name="Salmonella enterica",
                taxonomy_id=28901,
                taxonomy_lvl="S",
                kraken_assigned_reads=30000,
                added_reads=20000,
                new_est_reads=50000,
                fraction_total_reads=0.35000,
            ),
        ]

    def test_basic_selection(self):
        """Test basic species selection with default thresholds."""
        species_list = self.create_test_species()
        selection = select_species(species_list, BrackenThresholds())

        assert selection is not None
        assert selection.species.name == "Listeria monocytogenes"
        assert selection.rank == 1

    def test_no_species_meet_thresholds(self):
        """Test when no species meet thresholds."""
        species_list = self.create_test_species()
        strict_thresholds = BrackenThresholds(min_frac=0.999, min_reads=2000000)

        selection = select_species(species_list, strict_thresholds)
        assert selection is None

    def test_tie_breaking(self):
        """Test tie-breaking logic: fraction > reads > name."""
        # Equal read counts, different fractions
        species_list = [
            BrackenSpecies(
                name="Species A",
                taxonomy_id=1,
                taxonomy_lvl="S",
                kraken_assigned_reads=50000,
                added_reads=50000,
                new_est_reads=100000,
                fraction_total_reads=0.75000,
            ),
            BrackenSpecies(
                name="Species B",
                taxonomy_id=2,
                taxonomy_lvl="S",
                kraken_assigned_reads=50000,
                added_reads=50000,
                new_est_reads=100000,
                fraction_total_reads=0.70000,
            ),
        ]

        selection = select_species(species_list, BrackenThresholds())
        assert selection is not None
        assert selection.species.name == "Species A"  # Higher fraction wins

    def test_custom_thresholds(self):
        """Test selection with custom thresholds."""
        species_list = self.create_test_species()
        loose_thresholds = BrackenThresholds(min_frac=0.30, min_reads=40000)

        selection = select_species(species_list, loose_thresholds)
        assert selection is not None
        assert selection.species.name == "Listeria monocytogenes"


class TestFileToSelectionWorkflow:
    """Test complete file-to-selection workflow."""

    def create_test_bracken_file(self, content: str) -> Path:
        """Helper to create temporary Bracken file."""
        temp_file = tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False)
        temp_file.write(content)
        temp_file.close()
        return Path(temp_file.name)

    def test_successful_selection_from_file(self):
        """Test successful species selection from file."""
        content = """name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction_total_reads
Listeria monocytogenes\t1639\tS\t1718713\t145289\t1864002\t0.99520"""

        file_path = self.create_test_bracken_file(content)

        try:
            selection = select_species_from_file(file_path)
            assert selection is not None
            assert selection.species.name == "Listeria monocytogenes"
        finally:
            file_path.unlink()

    def test_selection_failure_with_helpful_error(self):
        """Test error message when no species meet thresholds."""
        content = """name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction_total_reads
Escherichia coli\t562\tS\t45000\t5000\t50000\t0.50000"""

        file_path = self.create_test_bracken_file(content)

        try:
            with pytest.raises(BrackenSelectionError) as exc_info:
                select_species_from_file(file_path)

            error_msg = str(exc_info.value)
            assert "No species" in error_msg
            assert "0.70 fraction" in error_msg
            assert "100000 reads" in error_msg
            assert "Consider using --species or --reference instead" in error_msg
        finally:
            file_path.unlink()


class TestThresholdValidation:
    """Test threshold validation."""

    def test_invalid_thresholds(self):
        """Test validation of invalid threshold values."""
        with pytest.raises(ValueError, match="min_frac must be between 0.0 and 1.0"):
            BrackenThresholds(min_frac=1.5)

        with pytest.raises(ValueError, match="min_reads must be >= 0"):
            BrackenThresholds(min_reads=-1)


class TestSpeciesNameNormalization:
    """Test species name normalization."""

    def test_basic_normalization(self):
        """Test basic species name normalization."""
        result = normalize_species_name_for_resolution("Listeria monocytogenes")
        assert result == "Listeria_monocytogenes"

    def test_strain_removal(self):
        """Test removal of strain/serovar information."""
        result = normalize_species_name_for_resolution(
            "Salmonella enterica serovar Typhimurium"
        )
        assert result == "Salmonella_enterica"

    def test_invalid_name_format(self):
        """Test error for invalid species format."""
        with pytest.raises(
            ValueError, match="Species name must contain at least genus and species"
        ):
            normalize_species_name_for_resolution("Bacteria")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
