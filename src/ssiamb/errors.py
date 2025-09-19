#!/usr/bin/env python3
"""
Centralized error handling and exit code mapping for ssiamb.

Exit code specification:
- 0: success
- 1: CLI/input errors (missing args/files, unsafe sample name, invalid manifest)
- 2: ref-mode selection errors
- 3: reuse mismatch too severe
- 4: external tool failure
- 5: QC failure (only if user opts into --qc-action fail)
"""

import typer
from pathlib import Path
from typing import List, Optional
from rich.console import Console

# Import all custom exceptions for mapping
from .io_utils import SampleNameError, TSVWriteError
from .mapping import MappingError, ExternalToolError
from .calling import VariantCallingError
from .vcf_ops import VCFOperationError
from .depth import DepthAnalysisError
from .reuse import CompatibilityError
from .bracken import BrackenSelectionError
from .refdir import ReferenceResolutionError

console = Console()


class ExitCode:
    """Exit code constants."""

    SUCCESS = 0
    INPUT_ERROR = 1
    REF_MODE_ERROR = 2
    REUSE_MISMATCH = 3
    EXTERNAL_TOOL_ERROR = 4
    QC_FAILURE = 5


def map_exception_to_exit_code(exception: Exception) -> int:
    """
    Map exception types to appropriate exit codes.

    Args:
        exception: The exception to map

    Returns:
        Appropriate exit code (1-5)
    """
    # CLI/Input errors (exit code 1)
    if isinstance(
        exception,
        (
            ValueError,  # Invalid CLI arguments
            FileNotFoundError,  # Missing input files
            SampleNameError,  # Invalid sample names
            TSVWriteError,  # File writing issues
            typer.BadParameter,  # Typer CLI parameter issues
        ),
    ):
        return ExitCode.INPUT_ERROR

    # Reference mode selection errors (exit code 2)
    if isinstance(
        exception,
        (
            ReferenceResolutionError,  # Can't resolve reference
            BrackenSelectionError,  # Bracken selection failed
        ),
    ):
        return ExitCode.REF_MODE_ERROR

    # Reuse compatibility errors (exit code 3)
    if isinstance(exception, CompatibilityError):
        return ExitCode.REUSE_MISMATCH

    # External tool failures (exit code 4)
    if isinstance(
        exception,
        (
            ExternalToolError,  # External tool not found or failed
            MappingError,  # Mapping tool failure
            VariantCallingError,  # Calling tool failure
            DepthAnalysisError,  # Depth analysis tool failure
            VCFOperationError,  # VCF processing tool failure
        ),
    ):
        return ExitCode.EXTERNAL_TOOL_ERROR

    # QC failures would be exit code 5, but are warn-only by default
    # This would only be used if --qc-action fail is implemented later

    # Default to input error for unknown exceptions
    return ExitCode.INPUT_ERROR


def format_helpful_error_message(exception: Exception) -> str:
    """
    Format exception with helpful suggestions.

    Args:
        exception: The exception to format

    Returns:
        User-friendly error message with suggestions
    """
    error_msg = str(exception)

    # Add helpful suggestions based on exception type
    if isinstance(exception, FileNotFoundError):
        missing_file = error_msg
        suggestions = [
            "Check that the file path is correct",
            "Ensure the file exists and is readable",
            "Use absolute paths to avoid confusion",
        ]
        return f"{missing_file}\n\nSuggestions:\n" + "\n".join(
            f"  • {s}" for s in suggestions
        )

    elif isinstance(exception, SampleNameError):
        suggestions = [
            "Use --sample to provide a valid sample name",
            "Sample names must match pattern: ^[A-Za-z0-9._-]{1,64}$",
            "Remove special characters or spaces from filenames",
        ]
        return f"{error_msg}\n\nSuggestions:\n" + "\n".join(
            f"  • {s}" for s in suggestions
        )

    elif isinstance(exception, ReferenceResolutionError):
        suggestions = [
            "Use --reference to provide a direct reference file",
            "Check that $SSIAMB_REF_DIR is set and contains reference files",
            "Use --species with a valid species name from the admin directory",
            "List available species with the admin directory contents",
        ]
        return f"{error_msg}\n\nSuggestions:\n" + "\n".join(
            f"  • {s}" for s in suggestions
        )

    elif isinstance(exception, BrackenSelectionError):
        suggestions = [
            "Use --species to bypass Bracken and select reference directly",
            "Use --reference to provide a specific reference file",
            "Check Bracken file format and species-level classifications",
            "Adjust --min-bracken-frac or --min-bracken-reads thresholds",
        ]
        return f"{error_msg}\n\nSuggestions:\n" + "\n".join(
            f"  • {s}" for s in suggestions
        )

    elif isinstance(exception, CompatibilityError):
        suggestions = [
            "Check that VCF and BAM are from the same reference",
            "Ensure contigs in VCF match those in the reference",
            "Use --force-lenient if minor mismatches are acceptable (when implemented)",
            "Regenerate VCF/BAM with consistent reference",
        ]
        return f"{error_msg}\n\nSuggestions:\n" + "\n".join(
            f"  • {s}" for s in suggestions
        )

    elif isinstance(exception, ExternalToolError):
        suggestions = [
            "Ensure required tools are installed and in PATH",
            "Check tool versions are compatible",
            "Use --dry-run to see which tools would be called",
            "Install tools via conda: conda install minimap2 bwa-mem2 samtools bcftools bbmap mosdepth",
        ]
        return f"{error_msg}\n\nSuggestions:\n" + "\n".join(
            f"  • {s}" for s in suggestions
        )

    elif isinstance(exception, (MappingError, VariantCallingError, DepthAnalysisError)):
        suggestions = [
            "Check input file formats and integrity",
            "Ensure sufficient disk space and memory",
            "Try with different --threads setting",
            "Check tool-specific logs for detailed error messages",
        ]
        return f"{error_msg}\n\nSuggestions:\n" + "\n".join(
            f"  • {s}" for s in suggestions
        )

    # Default message for unknown exceptions
    return error_msg


def handle_exception_with_exit(exception: Exception, context: str = "") -> None:
    """
    Handle exception with appropriate exit code and helpful message.

    Args:
        exception: The exception to handle
        context: Additional context for the error
    """
    exit_code = map_exception_to_exit_code(exception)
    error_message = format_helpful_error_message(exception)

    # Add context if provided
    if context:
        full_message = f"{context}: {error_message}"
    else:
        full_message = error_message

    # Print error with appropriate styling
    console.print(f"[red]Error: {full_message}[/red]")

    # Exit with appropriate code
    raise typer.Exit(exit_code)


def suggest_species_from_directory(ref_dir: Path) -> List[str]:
    """
    List available species in reference directory for helpful error messages.

    Args:
        ref_dir: Path to reference directory

    Returns:
        List of available species names
    """
    species_list = []
    if ref_dir and ref_dir.exists():
        try:
            for fna_file in ref_dir.glob("*.fna"):
                # Extract species name from filename (remove .fna extension)
                species_name = fna_file.stem
                species_list.append(species_name)
        except Exception:
            pass  # Ignore errors when listing directory

    return sorted(species_list)


def create_species_suggestion_message(ref_dir: Optional[Path]) -> str:
    """
    Create helpful message with available species list.

    Args:
        ref_dir: Path to reference directory

    Returns:
        Formatted message with available species
    """
    if not ref_dir:
        return "Set $SSIAMB_REF_DIR or use --ref-dir to specify reference directory."

    species_list = suggest_species_from_directory(ref_dir)

    if not species_list:
        return f"No species references found in {ref_dir}. Check directory contents."

    species_msg = "Available species:\n" + "\n".join(
        f"  • {species}" for species in species_list[:10]
    )
    if len(species_list) > 10:
        species_msg += f"\n  ... and {len(species_list) - 10} more"

    return species_msg
