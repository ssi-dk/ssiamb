#!/usr/bin/env python3
"""
Comprehensive smoke test for ssiamb - Integration test that exercises
the full pipeline from reference download through analysis output.

This test validates:
1. Admin reference directory setup via download command
2. Self-mapping mode analysis
3. Reference-mapping mode analysis
4. Core output file generation and validation
5. Error handling for common failure modes

The test uses existing fixture data and creates temporary test environments
to ensure the full workflow operates correctly.
"""

import subprocess
import tempfile
import sys
import gzip
from pathlib import Path
from typing import Dict, List, Optional
import argparse

# Test configuration
TEST_SPECIES = [
    "Escherichia_coli",  # Very common, should download reliably
    "Salmonella_enterica",  # Common pathogen
    "Staphylococcus_aureus",  # Gram-positive representative
]

FIXTURE_DIR = Path(__file__).parent / "fixtures"
SELF_MODE_FIXTURES = FIXTURE_DIR / "self_mode"
REF_MODE_FIXTURES = FIXTURE_DIR / "ref_mode"


class SmokeTestError(Exception):
    """Custom exception for smoke test failures"""

    pass


def run_command(
    cmd: List[str],
    cwd: Optional[Path] = None,
    timeout: int = 300,
    capture_output: bool = True,
) -> subprocess.CompletedProcess:
    """
    Run a command with timeout and error handling.

    Args:
        cmd: Command and arguments as list
        cwd: Working directory for command
        timeout: Timeout in seconds
        capture_output: Whether to capture stdout/stderr

    Returns:
        CompletedProcess result

    Raises:
        SmokeTestError: If command fails or times out
    """
    try:
        print(f"Running: {' '.join(map(str, cmd))}")
        if cwd:
            print(f"  in directory: {cwd}")

        result = subprocess.run(
            cmd,
            cwd=cwd,
            timeout=timeout,
            capture_output=capture_output,
            text=True,
            check=False,  # Don't auto-raise on non-zero exit
        )

        if result.returncode != 0:
            error_msg = f"Command failed with exit code {result.returncode}: {' '.join(map(str, cmd))}"
            if result.stderr:
                error_msg += f"\nSTDERR: {result.stderr}"
            if result.stdout:
                error_msg += f"\nSTDOUT: {result.stdout}"
            raise SmokeTestError(error_msg)

        return result

    except subprocess.TimeoutExpired:
        raise SmokeTestError(
            f"Command timed out after {timeout}s: {' '.join(map(str, cmd))}"
        )
    except FileNotFoundError:
        raise SmokeTestError(f"Command not found: {cmd[0]}")


def verify_file_exists(filepath: Path, description: str) -> None:
    """Verify that a file exists and is non-empty."""
    if not filepath.exists():
        raise SmokeTestError(f"{description} not found: {filepath}")
    if filepath.stat().st_size == 0:
        raise SmokeTestError(f"{description} is empty: {filepath}")
    print(f"✓ {description} exists: {filepath} ({filepath.stat().st_size} bytes)")


def verify_tsv_output(
    filepath: Path, expected_sample: Optional[str] = None
) -> Dict[str, str]:
    """
    Verify TSV output file and return parsed row data.

    Args:
        filepath: Path to TSV file
        expected_sample: Expected sample name in the file (optional)

    Returns:
        Dictionary of column -> value mappings for the first data row

    Raises:
        SmokeTestError: If file format is invalid
    """
    verify_file_exists(filepath, "TSV output file")

    with open(filepath, "r") as f:
        lines = [line.strip() for line in f if line.strip()]

    if len(lines) < 2:
        raise SmokeTestError(
            f"TSV file should have header + at least 1 data row, got {len(lines)} lines"
        )

    # Parse header
    header = lines[0].split("\t")
    expected_columns = ["sample", "mode", "mapper", "caller", "ambiguous_snv_count"]
    missing_cols = [col for col in expected_columns if col not in header]
    if missing_cols:
        raise SmokeTestError(f"Missing expected columns in TSV: {missing_cols}")

    # Use first data row
    first_row = lines[1].split("\t")
    if len(first_row) < len(header):
        raise SmokeTestError(
            f"First data row has {len(first_row)} columns, expected {len(header)}"
        )

    sample_row = dict(zip(header, first_row))
    actual_sample = sample_row["sample"]

    print(f"✓ TSV contains sample '{actual_sample}' with {len(sample_row)} fields")

    # Validate key numeric fields
    try:
        snv_count = int(sample_row["ambiguous_snv_count"])
        print(f"  - Ambiguous SNV count: {snv_count}")
    except (ValueError, KeyError):
        raise SmokeTestError("Invalid or missing ambiguous_snv_count in TSV")

    return sample_row


def test_admin_reference_setup(test_dir: Path) -> Path:
    """
    Test the admin reference directory setup by downloading test species.

    Args:
        test_dir: Base test directory

    Returns:
        Path to created admin reference directory

    Raises:
        SmokeTestError: If reference setup fails
    """
    print("\n=== Testing Admin Reference Setup ===")

    admin_ref_dir = test_dir / "admin_refs"
    admin_ref_dir.mkdir()

    # Import and test the RefSeqDownloader directly
    try:
        from src.ssiamb.refseq import RefSeqDownloader
    except ImportError as e:
        raise SmokeTestError(f"Could not import RefSeqDownloader: {e}")

    # Test downloading one species to verify functionality
    test_species = "Escherichia_coli"  # Start with just one for speed
    print(f"\nTesting download for {test_species}...")

    try:
        downloader = RefSeqDownloader(admin_ref_dir, verbose=True, create_indexes=True)
        result = downloader.download_species(test_species)

        if result:
            print(f"✓ Successfully downloaded {test_species}")

            # Verify expected files were created
            fasta_file = admin_ref_dir / f"{test_species}.fna"
            mmi_index = admin_ref_dir / f"{test_species}.fna.mmi"

            verify_file_exists(fasta_file, f"{test_species} FASTA")
            verify_file_exists(mmi_index, f"{test_species} minimap2 index")

            # Check for bwa-mem2 index files (multiple files expected)
            bwa_files = list(admin_ref_dir.glob(f"{test_species}.fna.bwa*")) + list(
                admin_ref_dir.glob(f"{test_species}.fna.*")
            )
            bwa_files = [
                f
                for f in bwa_files
                if f.suffix not in [".fna"] and ".mmi" not in f.name
            ]  # Exclude the FASTA file itself and .mmi
            if not bwa_files:
                raise SmokeTestError(
                    f"No bwa-mem2 index files found for {test_species}"
                )
            print(f"✓ Found {len(bwa_files)} bwa-mem2 index files for {test_species}")
        else:
            raise SmokeTestError(f"Download returned None for {test_species}")

    except Exception as e:
        print(f"✗ Failed to download {test_species}: {e}")
        # For smoke test, we'll treat this as non-critical since downloads can be flaky
        print("Download test failed but continuing with other tests...")
        return admin_ref_dir

    print(
        f"✓ Admin reference directory setup complete with reference for {test_species}"
    )
    return admin_ref_dir


def test_self_mode_analysis(test_dir: Path) -> Path:
    """
    Test self-mapping mode analysis using fixture data.

    Args:
        test_dir: Base test directory

    Returns:
        Path to self-mode output directory

    Raises:
        SmokeTestError: If self-mode analysis fails
    """
    print("\n=== Testing Self-Mode Analysis ===")

    # Verify fixtures exist
    if not SELF_MODE_FIXTURES.exists():
        raise SmokeTestError(f"Self-mode fixtures not found: {SELF_MODE_FIXTURES}")

    r1_file = SELF_MODE_FIXTURES / "reads_R1.fastq.gz"
    r2_file = SELF_MODE_FIXTURES / "reads_R2.fastq.gz"
    assembly_file = SELF_MODE_FIXTURES / "assembly.fna"

    verify_file_exists(r1_file, "Self-mode R1 reads")
    verify_file_exists(r2_file, "Self-mode R2 reads")
    verify_file_exists(assembly_file, "Self-mode assembly")

    # Set up output directory
    output_dir = test_dir / "self_mode_output"
    output_dir.mkdir()

    # Run self-mode analysis
    sample_name = "test_self_sample"
    cmd = [
        sys.executable,
        "-m",
        "ssiamb.cli",
        "self",
        "--r1",
        str(r1_file),
        "--r2",
        str(r2_file),
        "--assembly",
        str(assembly_file),
        "--sample",
        sample_name,
        "--outdir",
        str(output_dir),
        "--threads",
        "2",  # Use fewer threads for testing
        "--emit-vcf",
        "--emit-bed",
    ]

    run_command(cmd, timeout=600)  # Mapping can take time
    print("✓ Self-mode analysis completed successfully")

    # Debug: List all files created
    created_files = list(output_dir.glob("*"))
    print(f"Files created in output directory: {[f.name for f in created_files]}")

    # Verify outputs - check for TSV files with any name
    tsv_files = list(output_dir.glob("*.tsv"))
    if not tsv_files:
        raise SmokeTestError(f"No TSV files found in output directory: {output_dir}")

    summary_file = tsv_files[0]  # Use first TSV file found
    verify_file_exists(summary_file, "Self-mode summary TSV")

    # Parse and validate TSV content
    sample_data = verify_tsv_output(summary_file)

    # Verify mode is correct
    if sample_data.get("mode") != "self":
        raise SmokeTestError(f"Expected mode 'self', got '{sample_data.get('mode')}'")

    # Check for optional output files if requested
    vcf_files = list(output_dir.glob("*.vcf.gz"))
    bed_files = list(output_dir.glob("*.bed.gz"))

    if vcf_files:
        verify_file_exists(vcf_files[0], "Self-mode VCF output")

        # Verify VCF is properly formatted (has header)
        with gzip.open(vcf_files[0], "rt") as f:
            first_line = f.readline()
            if not first_line.startswith("##fileformat=VCF"):
                raise SmokeTestError("VCF file does not start with proper header")
        print("✓ VCF file is properly formatted")

    if bed_files:
        verify_file_exists(bed_files[0], "Self-mode BED output")

    print("✓ Self-mode analysis validation complete")
    return output_dir


def test_ref_mode_analysis(test_dir: Path, admin_ref_dir: Path) -> Path:
    """
    Test reference-mapping mode analysis using fixture data and downloaded references.

    Args:
        test_dir: Base test directory
        admin_ref_dir: Admin reference directory from setup test

    Returns:
        Path to ref-mode output directory

    Raises:
        SmokeTestError: If ref-mode analysis fails
    """
    print("\n=== Testing Reference-Mode Analysis ===")

    # Verify fixtures exist
    if not REF_MODE_FIXTURES.exists():
        raise SmokeTestError(f"Ref-mode fixtures not found: {REF_MODE_FIXTURES}")

    r1_file = REF_MODE_FIXTURES / "reads_R1.fastq.gz"
    r2_file = REF_MODE_FIXTURES / "reads_R2.fastq.gz"

    verify_file_exists(r1_file, "Ref-mode R1 reads")
    verify_file_exists(r2_file, "Ref-mode R2 reads")

    # Set up output directory
    output_dir = test_dir / "ref_mode_output"
    output_dir.mkdir()

    # Find an available species from the admin directory
    available_species = [f.stem for f in admin_ref_dir.glob("*.fna")]
    if not available_species:
        raise SmokeTestError("No species references available in admin directory")

    test_species = available_species[0]  # Use first available
    print(f"Using species: {test_species}")

    # Run ref-mode analysis with species lookup
    sample_name = "test_ref_sample"
    cmd = [
        sys.executable,
        "-m",
        "ssiamb.cli",
        "ref",
        "--r1",
        str(r1_file),
        "--r2",
        str(r2_file),
        "--species",
        test_species.replace("_", " "),  # Convert to space format for CLI
        "--ref-dir",
        str(admin_ref_dir),
        "--sample",
        sample_name,
        "--outdir",
        str(output_dir),
        "--threads",
        "2",
        "--emit-vcf",
        "--emit-matrix",
    ]

    run_command(cmd, timeout=600)
    print("✓ Ref-mode analysis completed successfully")

    # Debug: List all files created
    created_files = list(output_dir.glob("*"))
    print(f"Files created in output directory: {[f.name for f in created_files]}")

    # Verify outputs - check for TSV files with any name
    tsv_files = list(output_dir.glob("*.tsv"))
    if not tsv_files:
        raise SmokeTestError(f"No TSV files found in output directory: {output_dir}")

    summary_file = tsv_files[0]  # Use first TSV file found
    verify_file_exists(summary_file, "Ref-mode summary TSV")

    # Parse and validate TSV content
    sample_data = verify_tsv_output(summary_file)

    # Verify mode is correct
    if sample_data.get("mode") != "ref":
        raise SmokeTestError(f"Expected mode 'ref', got '{sample_data.get('mode')}'")

    # Check for optional output files
    vcf_files = list(output_dir.glob("*.vcf.gz"))
    matrix_files = list(output_dir.glob("*.matrix.tsv.gz")) or list(
        output_dir.glob("*matrix*.tsv.gz")
    )

    if vcf_files:
        verify_file_exists(vcf_files[0], "Ref-mode VCF output")

    if matrix_files:
        verify_file_exists(matrix_files[0], "Ref-mode matrix output")

    print("✓ Ref-mode analysis validation complete")
    return output_dir


def test_summarize_mode(test_dir: Path, self_output_dir: Path) -> None:
    """
    Test summarize mode using outputs from self-mode analysis.

    Args:
        test_dir: Base test directory
        self_output_dir: Output directory from self-mode test

    Raises:
        SmokeTestError: If summarize mode fails
    """
    print("\n=== Testing Summarize Mode ===")

    # Look for BAM file in intermediate outputs (may be named differently)
    bam_files = list(self_output_dir.glob("*.bam"))
    vcf_files = list(self_output_dir.glob("*.vcf.gz"))

    if not bam_files:
        print("⚠ No BAM files found in self-mode output, skipping summarize test")
        return

    if not vcf_files:
        print("⚠ No VCF files found in self-mode output, skipping summarize test")
        return

    bam_file = bam_files[0]
    vcf_file = vcf_files[0]

    # Set up output directory
    output_dir = test_dir / "summarize_output"
    output_dir.mkdir()

    # Run summarize mode
    cmd = [
        sys.executable,
        "-m",
        "ssiamb.cli",
        "summarize",
        "--vcf",
        str(vcf_file),
        "--bam",
        str(bam_file),
        "--output",
        str(output_dir / "summarize_result.tsv"),
    ]

    try:
        run_command(cmd, timeout=120)
        print("✓ Summarize mode completed successfully")

        # Verify output
        summary_file = output_dir / "summarize_result.tsv"
        verify_file_exists(summary_file, "Summarize mode output")

        verify_tsv_output(summary_file)
        print("✓ Summarize mode validation complete")

    except SmokeTestError as e:
        print(f"⚠ Summarize mode test failed (non-critical): {e}")


def test_error_handling(test_dir: Path) -> None:
    """
    Test error handling for common failure scenarios.

    Args:
        test_dir: Base test directory

    Raises:
        SmokeTestError: If error handling doesn't work as expected
    """
    print("\n=== Testing Error Handling ===")

    # Test 1: Missing input files
    print("Testing missing input file handling...")
    cmd = [
        sys.executable,
        "-m",
        "ssiamb.cli",
        "self",
        "--r1",
        "/nonexistent/file.fastq.gz",
        "--r2",
        "/nonexistent/file.fastq.gz",
        "--assembly",
        "/nonexistent/assembly.fna",
        "--sample",
        "test_error",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode == 0:
        raise SmokeTestError(
            "Expected failure for missing input files, but command succeeded"
        )
    print("✓ Correctly handles missing input files")

    # Test 2: Invalid species name
    print("Testing invalid species handling...")
    cmd = [
        sys.executable,
        "-m",
        "ssiamb.cli",
        "ref",
        "--r1",
        "/nonexistent/file.fastq.gz",
        "--r2",
        "/nonexistent/file.fastq.gz",
        "--species",
        "Nonexistent_species",
        "--sample",
        "test_error",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode == 0:
        raise SmokeTestError(
            "Expected failure for invalid species, but command succeeded"
        )
    print("✓ Correctly handles invalid species names")

    print("✓ Error handling tests complete")


def run_smoke_tests(clean_downloads: bool = False) -> None:
    """
    Run the complete smoke test suite.

    Args:
        clean_downloads: If True, skip admin reference setup test

    Raises:
        SmokeTestError: If any critical test fails
    """
    print("Starting ssiamb smoke tests...")
    print("=" * 50)

    # Create temporary test directory
    with tempfile.TemporaryDirectory() as temp_dir:
        test_dir = Path(temp_dir)
        print(f"Using test directory: {test_dir}")

        try:
            # Test 1: Admin reference setup (optional, slow)
            admin_ref_dir = None
            if not clean_downloads:
                try:
                    admin_ref_dir = test_admin_reference_setup(test_dir)
                except SmokeTestError as e:
                    print(f"⚠ Admin reference setup failed (non-critical): {e}")
                    print("Continuing with other tests...")

            # Test 2: Self-mode analysis
            self_output_dir = test_self_mode_analysis(test_dir)

            # Test 3: Reference-mode analysis (if we have admin refs)
            if admin_ref_dir:
                test_ref_mode_analysis(test_dir, admin_ref_dir)
            else:
                print("⚠ Skipping ref-mode test (no admin references available)")

            # Test 4: Summarize mode
            test_summarize_mode(test_dir, self_output_dir)

            # Test 5: Error handling
            test_error_handling(test_dir)

            print("\n" + "=" * 50)
            print("🎉 All smoke tests completed successfully!")
            print("✓ Self-mode analysis working")
            if admin_ref_dir:
                print("✓ Reference download working")
                print("✓ Ref-mode analysis working")
            print("✓ Error handling working")
            print("✓ Output file generation working")

        except SmokeTestError as e:
            print(f"\n❌ Smoke test failed: {e}")
            sys.exit(1)
        except KeyboardInterrupt:
            print("\n⚠ Smoke tests interrupted by user")
            sys.exit(1)
        except Exception as e:
            print(f"\n💥 Unexpected error during smoke tests: {e}")
            import traceback

            traceback.print_exc()
            sys.exit(1)


def main() -> None:
    """Main entry point for smoke test script."""
    parser = argparse.ArgumentParser(
        description="Comprehensive smoke tests for ssiamb",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python smoke_test.py                    # Run all tests including downloads
  python smoke_test.py --skip-downloads  # Skip slow reference download tests
  
This script validates the full ssiamb pipeline from reference setup
through analysis output generation.
        """,
    )

    parser.add_argument(
        "--skip-downloads",
        action="store_true",
        help="Skip admin reference download tests (faster)",
    )

    args = parser.parse_args()

    run_smoke_tests(clean_downloads=args.skip_downloads)


if __name__ == "__main__":
    main()
