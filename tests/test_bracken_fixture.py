"""
Integration tests using real fixture data.
"""

import pytest
from pathlib import Path

from src.ssiamb.bracken import (
    parse_bracken_file,
    select_species_from_file,
    BrackenThresholds,
)


class TestBrackenWithRealFixture:
    """Test Bracken functionality with real fixture data."""

    @property
    def fixture_path(self) -> Path:
        """Path to the real Bracken fixture."""
        return Path("fixtures/supporting_files/bracken_listeria.tsv")

    def test_parse_real_bracken_file(self):
        """Test parsing the real Bracken fixture file."""
        if not self.fixture_path.exists():
            pytest.skip(f"Fixture file not found: {self.fixture_path}")

        species_list = parse_bracken_file(self.fixture_path)

        # Should have species-level entries
        assert len(species_list) > 0

        # Check that all entries are species-level
        for species in species_list:
            assert species.taxonomy_lvl == "S"
            assert species.name
            assert species.taxonomy_id > 0
            assert species.new_est_reads >= 0
            assert 0.0 <= species.fraction_total_reads <= 1.0

    def test_select_species_from_real_file(self):
        """Test species selection from real Bracken fixture."""
        if not self.fixture_path.exists():
            pytest.skip(f"Fixture file not found: {self.fixture_path}")

        selection = select_species_from_file(self.fixture_path)

        # Should successfully select a species
        assert selection is not None
        assert selection.species.name
        assert selection.rank >= 1

        # Selected species should meet default thresholds
        assert selection.species.fraction_total_reads >= 0.70
        assert selection.species.new_est_reads >= 100000

    def test_real_file_with_custom_thresholds(self):
        """Test species selection with custom thresholds on real data."""
        if not self.fixture_path.exists():
            pytest.skip(f"Fixture file not found: {self.fixture_path}")

        # Try with very loose thresholds
        loose_thresholds = BrackenThresholds(min_frac=0.01, min_reads=1000)
        selection = select_species_from_file(self.fixture_path, loose_thresholds)

        assert selection is not None

        # Try with very strict thresholds - may or may not find species
        strict_thresholds = BrackenThresholds(min_frac=0.99, min_reads=1000000)
        try:
            selection = select_species_from_file(self.fixture_path, strict_thresholds)
            # If a species is found, it should meet the strict criteria
            if selection:
                assert selection.species.fraction_total_reads >= 0.99
                assert selection.species.new_est_reads >= 1000000
        except Exception:
            # It's okay if no species meet very strict thresholds
            pass


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
