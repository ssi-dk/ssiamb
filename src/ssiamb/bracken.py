"""
Bracken output parsing and species selection for ssiamb.

This module handles:
- Parsing Bracken abundance estimation files
- Species-level selection based on fraction and read count thresholds
- Tie-breaking logic for species selection
"""

from __future__ import annotations
import pandas as pd
from pathlib import Path
from typing import Optional, NamedTuple, List
from dataclasses import dataclass
import logging

logger = logging.getLogger(__name__)


class BrackenSelectionError(Exception):
    """Raised when species selection from Bracken fails."""

    pass


@dataclass
class BrackenSpecies:
    """
    A single species entry from Bracken output.

    Represents a species-level taxonomy assignment with abundance estimates.
    """

    name: str  # Species name from Bracken
    taxonomy_id: int  # NCBI taxonomy ID
    taxonomy_lvl: str  # Taxonomy level (should be "S" for species)
    kraken_assigned_reads: int  # Reads assigned to this taxon by Kraken
    added_reads: int  # Additional reads assigned by Bracken
    new_est_reads: int  # Total estimated reads for this taxon
    fraction_total_reads: float  # Fraction of total reads in sample

    @classmethod
    def from_row(cls, row: pd.Series) -> "BrackenSpecies":
        """Create BrackenSpecies from a pandas DataFrame row."""
        return cls(
            name=str(row["name"]),
            taxonomy_id=int(row["taxonomy_id"]),
            taxonomy_lvl=str(row["taxonomy_lvl"]),
            kraken_assigned_reads=int(row["kraken_assigned_reads"]),
            added_reads=int(row["added_reads"]),
            new_est_reads=int(row["new_est_reads"]),
            fraction_total_reads=float(row["fraction_total_reads"]),
        )


class BrackenSelection(NamedTuple):
    """Result of species selection from Bracken."""

    species: BrackenSpecies
    rank: int  # 1-based rank (1 = top species)


@dataclass
class BrackenThresholds:
    """Thresholds for species selection from Bracken."""

    min_frac: float = 0.70  # Minimum fraction of total reads
    min_reads: int = 100000  # Minimum number of reads

    def __post_init__(self) -> None:
        """Validate thresholds."""
        if not 0.0 <= self.min_frac <= 1.0:
            raise ValueError("min_frac must be between 0.0 and 1.0")
        if self.min_reads < 0:
            raise ValueError("min_reads must be >= 0")


def parse_bracken_file(path: Path) -> List[BrackenSpecies]:
    """
    Parse a Bracken abundance estimation file.

    Args:
        path: Path to Bracken output file

    Returns:
        List of species-level entries from the file

    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If file format is invalid
    """
    if not path.exists():
        raise FileNotFoundError(f"Bracken file not found: {path}")

    try:
        # Read the Bracken file (tab-separated)
        df = pd.read_csv(path, sep="\t")

        # Validate required columns
        required_columns = [
            "name",
            "taxonomy_id",
            "taxonomy_lvl",
            "kraken_assigned_reads",
            "added_reads",
            "new_est_reads",
            "fraction_total_reads",
        ]

        missing_columns = [col for col in required_columns if col not in df.columns]
        if missing_columns:
            raise ValueError(
                f"Missing required columns in Bracken file: {missing_columns}"
            )

        # Filter to species-level entries only
        species_df = df[df["taxonomy_lvl"] == "S"].copy()

        if species_df.empty:
            logger.warning(f"No species-level entries found in Bracken file: {path}")
            return []

        # Convert to BrackenSpecies objects
        species_list = [
            BrackenSpecies.from_row(row) for _, row in species_df.iterrows()
        ]

        logger.info(f"Parsed {len(species_list)} species-level entries from {path}")
        return species_list

    except pd.errors.EmptyDataError:
        raise ValueError(f"Bracken file is empty: {path}")
    except pd.errors.ParserError as e:
        raise ValueError(f"Failed to parse Bracken file {path}: {e}")
    except Exception as e:
        raise ValueError(f"Error reading Bracken file {path}: {e}")


def select_species(
    species_list: List[BrackenSpecies], thresholds: BrackenThresholds
) -> Optional[BrackenSelection]:
    """
    Select the best species from Bracken results based on thresholds.

    Selection algorithm:
    1. Filter species that meet both fraction and read count thresholds
    2. Tie-breaking: highest fraction, then highest reads, then alphabetical

    Args:
        species_list: List of species from Bracken parsing
        thresholds: Selection thresholds

    Returns:
        BrackenSelection if a species meets criteria, None otherwise
    """
    if not species_list:
        logger.warning("No species provided for selection")
        return None

    # Filter species that meet thresholds
    candidates = [
        sp
        for sp in species_list
        if (
            sp.fraction_total_reads >= thresholds.min_frac
            and sp.new_est_reads >= thresholds.min_reads
        )
    ]

    if not candidates:
        logger.info(
            f"No species meet thresholds (min_frac={thresholds.min_frac}, "
            f"min_reads={thresholds.min_reads})"
        )
        return None

    # Sort by tie-breaking criteria:
    # 1. Highest fraction (descending)
    # 2. Highest read count (descending)
    # 3. Alphabetical by name (ascending, for stability)
    candidates.sort(
        key=lambda sp: (-sp.fraction_total_reads, -sp.new_est_reads, sp.name)
    )

    selected = candidates[0]
    rank = species_list.index(selected) + 1  # 1-based rank in original list

    logger.info(
        f"Selected species: {selected.name} "
        f"(frac={selected.fraction_total_reads:.3f}, "
        f"reads={selected.new_est_reads}, rank={rank})"
    )

    return BrackenSelection(species=selected, rank=rank)


def select_species_from_file(
    path: Path, thresholds: Optional[BrackenThresholds] = None
) -> Optional[BrackenSelection]:
    """
    Parse Bracken file and select the best species.

    Args:
        path: Path to Bracken output file
        thresholds: Selection thresholds (uses defaults if None)

    Returns:
        BrackenSelection if successful, None if no species meets criteria

    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If file format is invalid
        BrackenSelectionError: If selection fails due to thresholds
    """
    if thresholds is None:
        thresholds = BrackenThresholds()

    # Parse the file
    species_list = parse_bracken_file(path)

    if not species_list:
        raise BrackenSelectionError(
            f"No species-level entries found in Bracken file: {path}"
        )

    # Select species
    selection = select_species(species_list, thresholds)

    if selection is None:
        # Provide helpful feedback about what was available
        top_species = sorted(
            species_list,
            key=lambda sp: (-sp.fraction_total_reads, -sp.new_est_reads, sp.name),
        )[
            :3
        ]  # Show top 3

        suggestions = []
        for sp in top_species:
            suggestions.append(
                f"  {sp.name}: {sp.fraction_total_reads:.3f} fraction, "
                f"{sp.new_est_reads} reads"
            )

        suggestion_text = "\n".join(suggestions)
        raise BrackenSelectionError(
            f"No species in {path} meet selection criteria:\n"
            f"  Required: {thresholds.min_frac:.2f} fraction, "
            f"{thresholds.min_reads} reads\n"
            f"Top candidates:\n{suggestion_text}\n"
            f"Consider using --species or --reference instead."
        )

    return selection


def normalize_species_name_for_resolution(bracken_name: str) -> str:
    """
    Normalize a Bracken species name for resolution in the admin reference directory.

    Bracken names often include strain information or other details that need
    to be stripped for reference resolution.

    Args:
        bracken_name: Species name from Bracken (e.g., "Escherichia coli O157:H7")

    Returns:
        Normalized name suitable for reference resolution (e.g., "Escherichia coli")
    """
    # Import here to avoid circular imports
    from .refdir import normalize_species_name

    # Remove strain/serovar information that commonly appears in Bracken
    # Examples: "Escherichia coli O157:H7" -> "Escherichia coli"
    #           "Salmonella enterica serovar Typhimurium" -> "Salmonella enterica"

    # Split and take first two words (genus + species)
    parts = bracken_name.strip().split()
    if len(parts) >= 2:
        genus_species = f"{parts[0]} {parts[1]}"
    else:
        genus_species = bracken_name

    # Use the standard normalization function
    return normalize_species_name(genus_species)
