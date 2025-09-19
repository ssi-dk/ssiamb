"""
Reference directory management and species resolution for ssiamb.

This module handles:
- Resolution of admin reference directories
- Species name normalization and aliasing
- Mapping species names to reference FASTA files
- Parsing accession numbers from FASTA headers
"""

import os
import re
from pathlib import Path
from typing import Optional, List
import logging

logger = logging.getLogger(__name__)


class ReferenceResolutionError(Exception):
    """Raised when reference resolution fails."""

    pass


def resolve_reference_directory(ref_dir: Optional[Path] = None) -> Path:
    """
    Resolve the admin reference directory path.

    Precedence:
    1. Provided ref_dir parameter
    2. SSIAMB_REF_DIR environment variable
    3. Default: ./refs/ (current working directory)

    Args:
        ref_dir: Explicit reference directory path

    Returns:
        Resolved reference directory path

    Raises:
        ReferenceResolutionError: If directory doesn't exist or isn't readable
    """
    if ref_dir is not None:
        resolved_dir = ref_dir
    elif "SSIAMB_REF_DIR" in os.environ:
        resolved_dir = Path(os.environ["SSIAMB_REF_DIR"])
    else:
        resolved_dir = Path("refs")

    # Check if directory exists and is readable
    if not resolved_dir.exists():
        raise ReferenceResolutionError(
            f"Reference directory does not exist: {resolved_dir}\\n"
            f"Please create the directory or set SSIAMB_REF_DIR environment variable."
        )

    if not resolved_dir.is_dir():
        raise ReferenceResolutionError(
            f"Reference path is not a directory: {resolved_dir}"
        )

    # Check if we can read the directory
    try:
        list(resolved_dir.iterdir())
    except PermissionError:
        raise ReferenceResolutionError(
            f"Cannot read reference directory: {resolved_dir}\\n"
            f"Please check permissions."
        )

    logger.info(f"Using reference directory: {resolved_dir}")
    return resolved_dir


def normalize_species_name(species: str) -> str:
    """
    Normalize species name to standardized format.

    Converts "Genus species ..." to "Genus_species" by:
    - Replacing spaces with underscores
    - Removing subspecies designations (subsp., subspecies)
    - Removing strain/serovar info
    - Removing punctuation except hyphens in species names
    - Taking only first two words (genus + species)

    Args:
        species: Raw species name

    Returns:
        Normalized species name in format "Genus_species"

    Examples:
        "Listeria monocytogenes" -> "Listeria_monocytogenes"
        "Salmonella enterica subsp. enterica" -> "Salmonella_enterica"
        "Escherichia coli O157:H7" -> "Escherichia_coli"
    """
    if not species or not species.strip():
        raise ValueError("Species name cannot be empty")

    # Start with the input, strip whitespace
    name = species.strip()

    # Remove common subspecies patterns
    subspecies_patterns = [
        r"\\s+subsp\\.?\\s+\\S+",
        r"\\s+subspecies\\s+\\S+",
        r"\\s+serovar\\s+\\S+",
        r"\\s+strain\\s+\\S+",
        r"\\s+var\\.?\\s+\\S+",
    ]

    for pattern in subspecies_patterns:
        name = re.sub(pattern, "", name, flags=re.IGNORECASE)

    # Split into words and take first two (genus + species)
    words = name.split()
    if len(words) < 2:
        raise ValueError(
            f"Species name must contain at least genus and species: '{species}'"
        )

    genus = words[0]
    species_epithet = words[1]

    # Clean up individual words - remove non-alphanumeric except hyphens
    genus = re.sub(r"[^a-zA-Z0-9-]", "", genus)
    species_epithet = re.sub(r"[^a-zA-Z0-9-]", "", species_epithet)

    if not genus or not species_epithet:
        raise ValueError(f"Invalid genus or species after cleaning: '{species}'")

    # Ensure proper capitalization
    genus = genus.capitalize()
    species_epithet = species_epithet.lower()

    normalized = f"{genus}_{species_epithet}"
    logger.debug(f"Normalized '{species}' -> '{normalized}'")

    return normalized


def resolve_species_alias(species: str) -> str:
    """
    Resolve species alias to canonical form.

    Args:
        species: Species name (normalized)

    Returns:
        Canonical species name (may be same as input)
    """
    # Import here to avoid circular imports
    from .config import get_config

    config = get_config()
    canonical = config.get_species_alias(species)
    if canonical != species:
        logger.info(f"Resolved species alias: '{species}' -> '{canonical}'")
    return canonical


def list_available_species(ref_dir: Path) -> List[str]:
    """
    List all available species in the reference directory.

    Looks for files matching pattern "Genus_species.fna"

    Args:
        ref_dir: Reference directory path

    Returns:
        List of available species names
    """
    species_list = []

    try:
        for fasta_file in ref_dir.glob("*.fna"):
            # Extract species name from filename
            name = fasta_file.stem
            # Validate it looks like a species name (Genus_species)
            if re.match(r"^[A-Z][a-z]+_[a-z]+$", name):
                species_list.append(name)
    except Exception as e:
        logger.warning(f"Error listing species in {ref_dir}: {e}")

    return sorted(species_list)


def resolve_species_to_fasta(
    ref_dir: Path, species_or_alias: str, require_indexes: bool = True
) -> Path:
    """
    Resolve species name to reference FASTA file.

    Args:
        ref_dir: Reference directory path
        species_or_alias: Species name or alias
        require_indexes: Whether to require mapping indexes exist

    Returns:
        Path to reference FASTA file

    Raises:
        ReferenceResolutionError: If species not found or indexes missing
    """
    # Normalize and resolve aliases
    try:
        normalized = normalize_species_name(species_or_alias)
    except ValueError as e:
        raise ReferenceResolutionError(f"Invalid species name: {e}")

    canonical = resolve_species_alias(normalized)

    # Look for FASTA file
    fasta_path = ref_dir / f"{canonical}.fna"

    if not fasta_path.exists():
        # List available species for helpful error message
        available = list_available_species(ref_dir)
        raise ReferenceResolutionError(
            f"Reference not found for species '{species_or_alias}' (normalized: '{canonical}')\\n"
            f"Expected: {fasta_path}\\n"
            f"Available species in {ref_dir}: {', '.join(available) if available else 'none'}"
        )

    # Check for required indexes if requested
    if require_indexes:
        missing_indexes = []

        # Check for minimap2 index
        mmi_path = fasta_path.with_suffix(".fna.mmi")
        if not mmi_path.exists():
            missing_indexes.append(f"minimap2 index ({mmi_path.name})")

        # Check for bwa-mem2 indexes
        bwa_extensions = [".0123", ".amb", ".ann", ".pac", ".bwt.2bit.64"]
        for ext in bwa_extensions:
            bwa_path = fasta_path.with_suffix(f".fna{ext}")
            if not bwa_path.exists():
                missing_indexes.append(f"bwa-mem2 index ({bwa_path.name})")
                break  # Only report first missing bwa index

        if missing_indexes:
            raise ReferenceResolutionError(
                f"Reference found but indexes missing for '{canonical}': {', '.join(missing_indexes)}\\n"
                f"Reference: {fasta_path}\\n"
                f"Please ensure all mapping indexes are present."
            )

    logger.info(f"Resolved species '{species_or_alias}' to {fasta_path}")
    return fasta_path


def parse_accession_from_fasta_header(fasta_path: Path) -> Optional[str]:
    """
    Parse accession number from FASTA header.

    Looks for common accession patterns in the first line:
    - RefSeq: NC_*, NZ_*, etc.
    - GenBank: CP*, AP*, etc.
    - Simple identifiers

    Args:
        fasta_path: Path to FASTA file

    Returns:
        Accession number if found, None otherwise
    """
    try:
        with open(fasta_path, "r") as f:
            first_line = f.readline().strip()

        if not first_line.startswith(">"):
            logger.warning(f"Invalid FASTA format in {fasta_path}")
            return None

        # Remove '>' and split by whitespace to get first field
        header = first_line[1:].split()[0] if first_line[1:].split() else ""

        # Look for common accession patterns
        accession_patterns = [
            r"(N[CTZW]_\d+(?:\.\d+)?)",  # RefSeq
            r"([A-Z]{2}\d+(?:\.\d+)?)",  # GenBank
            r"(\w+)",  # Fallback: any word
        ]

        for pattern in accession_patterns:
            match = re.search(pattern, header)
            if match:
                accession = match.group(1)
                logger.debug(f"Parsed accession '{accession}' from {fasta_path}")
                return accession

        logger.warning(f"Could not parse accession from header: {first_line}")
        return None

    except Exception as e:
        logger.warning(f"Error reading FASTA header from {fasta_path}: {e}")
        return None
