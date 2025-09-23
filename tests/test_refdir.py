"""
Unit tests for reference directory operations.

Tests species name normalization, reference resolution,
and directory handling functionality.
"""

import pytest
import tempfile
from pathlib import Path
from unittest.mock import patch

from src.ssiamb.refdir import (
    normalize_species_name,
    resolve_reference_directory,
    ReferenceResolutionError,
)


class TestSpeciesNormalization:
    """Test species name normalization functionality."""

    def test_basic_normalization(self):
        """Test basic genus species normalization."""
        result = normalize_species_name("Listeria monocytogenes")
        assert result == "Listeria_monocytogenes"

    def test_capitalization_normalization(self):
        """Test proper capitalization handling."""
        result = normalize_species_name("listeria MONOCYTOGENES")
        assert result == "Listeria_monocytogenes"

    def test_subspecies_removal(self):
        """Test removal of subspecies designations."""
        # Test subsp. pattern
        result = normalize_species_name("Salmonella enterica subsp. enterica")
        assert result == "Salmonella_enterica"

        # Test subspecies pattern
        result = normalize_species_name("Bacillus subtilis subspecies spizizenii")
        assert result == "Bacillus_subtilis"

        # Test serovar pattern
        result = normalize_species_name("Salmonella enterica serovar Typhimurium")
        assert result == "Salmonella_enterica"

        # Test strain pattern
        result = normalize_species_name("Escherichia coli strain K-12")
        assert result == "Escherichia_coli"

        # Test var pattern
        result = normalize_species_name("Pseudomonas aeruginosa var. aeruginosa")
        assert result == "Pseudomonas_aeruginosa"

    def test_strain_info_removal(self):
        """Test removal of strain and serotype information."""
        result = normalize_species_name("Escherichia coli O157:H7")
        assert result == "Escherichia_coli"

        result = normalize_species_name("Staphylococcus aureus MRSA")
        assert result == "Staphylococcus_aureus"

    def test_punctuation_cleaning(self):
        """Test removal of unwanted punctuation."""
        result = normalize_species_name("Clostridium.difficile test")
        assert result == "Clostridiumdifficile_test"

        result = normalize_species_name("Acinetobacter.baumannii other")
        assert result == "Acinetobacterbaumannii_other"

    def test_hyphen_in_single_word_error(self):
        """Test that hyphens without spaces cause insufficient words error."""
        with pytest.raises(ValueError, match="must contain at least genus and species"):
            normalize_species_name("Clostridium-difficile")

    def test_hyphen_preservation(self):
        """Test that hyphens in species names are preserved."""
        result = normalize_species_name("Mycobacterium avium-intracellulare")
        assert result == "Mycobacterium_avium-intracellulare"

    def test_extra_words_removal(self):
        """Test removal of extra words beyond genus and species."""
        result = normalize_species_name("Klebsiella pneumoniae complex group")
        assert result == "Klebsiella_pneumoniae"

    def test_whitespace_handling(self):
        """Test handling of various whitespace patterns."""
        result = normalize_species_name("  Enterococcus   faecalis  ")
        assert result == "Enterococcus_faecalis"

        result = normalize_species_name("Streptococcus\tpyogenes\n")
        assert result == "Streptococcus_pyogenes"

    def test_empty_input_error(self):
        """Test error handling for empty input."""
        with pytest.raises(ValueError, match="Species name cannot be empty"):
            normalize_species_name("")

        with pytest.raises(ValueError, match="Species name cannot be empty"):
            normalize_species_name("   ")

    def test_insufficient_words_error(self):
        """Test error handling for insufficient words."""
        with pytest.raises(ValueError, match="must contain at least genus and species"):
            normalize_species_name("Monogenus")

    def test_invalid_characters_error(self):
        """Test error handling for completely invalid characters."""
        with pytest.raises(ValueError, match="Invalid genus or species after cleaning"):
            normalize_species_name("!@# $%^")

    def test_numeric_species_names(self):
        """Test that numeric species names are allowed."""
        result = normalize_species_name("123 456")
        assert result == "123_456"

    def test_edge_cases(self):
        """Test various edge cases."""
        # Numbers in species names
        result = normalize_species_name("Candidatus Liberibacter solanacearum")
        assert result == "Candidatus_liberibacter"

        # Single letter species
        result = normalize_species_name("Hepatitis B")
        assert result == "Hepatitis_b"

        # Hyphenated genus
        result = normalize_species_name("Alpha-proteobacterium species")
        assert result == "Alpha-proteobacterium_species"


class TestReferenceDirectoryResolution:
    """Test reference directory resolution functionality."""

    def test_resolve_existing_directory(self):
        """Test resolution of existing directory."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            result = resolve_reference_directory(temp_path)
            assert result == temp_path

    def test_resolve_nonexistent_directory_error(self):
        """Test error for non-existent directory."""
        nonexistent = Path("/this/does/not/exist")
        with pytest.raises(ReferenceResolutionError, match="does not exist"):
            resolve_reference_directory(nonexistent)

    def test_resolve_file_instead_of_directory_error(self):
        """Test error when path points to file instead of directory."""
        with tempfile.NamedTemporaryFile() as temp_file:
            file_path = Path(temp_file.name)
            with pytest.raises(ReferenceResolutionError, match="not a directory"):
                resolve_reference_directory(file_path)

    @patch.dict("os.environ", {}, clear=True)  # Clear SSIAMB_REF_DIR
    def test_resolve_none_uses_default(self):
        """Test that None input uses default directory."""
        with tempfile.TemporaryDirectory() as temp_dir:
            # Create refs subdirectory
            refs_dir = Path(temp_dir) / "refs"
            refs_dir.mkdir()

            # Change to temp directory and test default
            import os

            old_cwd = os.getcwd()
            try:
                os.chdir(temp_dir)
                result = resolve_reference_directory(None)
                assert result == Path("refs")
            finally:
                os.chdir(old_cwd)

    def test_resolve_string_path(self):
        """Test resolution with Path object input."""
        with tempfile.TemporaryDirectory() as temp_dir:
            result = resolve_reference_directory(Path(temp_dir))
            assert result == Path(temp_dir)
            assert isinstance(result, Path)

    def test_resolve_environment_variable(self):
        """Test resolution using environment variable."""
        with tempfile.TemporaryDirectory() as temp_dir:
            with patch.dict("os.environ", {"SSIAMB_REF_DIR": temp_dir}):
                result = resolve_reference_directory()
                assert result == Path(temp_dir)

    def test_resolve_default_path(self):
        """Test default resolution to 'refs' directory."""
        with tempfile.TemporaryDirectory() as temp_dir:
            # Create a refs subdirectory
            refs_dir = Path(temp_dir) / "refs"
            refs_dir.mkdir()

            # Change to temp directory and test default
            import os

            old_cwd = os.getcwd()
            try:
                os.chdir(temp_dir)
                result = resolve_reference_directory()
                assert result == Path("refs")
            finally:
                os.chdir(old_cwd)


class TestReferenceDirectoryPermissions:
    """Test reference directory permission handling."""

    @patch("src.ssiamb.refdir.Path.iterdir")
    def test_unreadable_directory_error(self, mock_iterdir):
        """Test error for unreadable directory."""
        mock_iterdir.side_effect = PermissionError("Permission denied")

        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            with pytest.raises(
                ReferenceResolutionError, match="Cannot read reference directory"
            ):
                resolve_reference_directory(temp_path)
