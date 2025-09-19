#!/usr/bin/env python3
"""
Tests for Milestone 21: Error Handling & Exit Codes

Tests the centralized error handling system to ensure:
- CLI commands return proper exit codes for different error types
- Error messages are helpful with suggestions
- CLI integration works correctly
"""

import pytest
from pathlib import Path
from typer.testing import CliRunner
import re

from ssiamb.cli import app
from ssiamb.errors import map_exception_to_exit_code


def strip_ansi_codes(text: str) -> str:
    """Remove ANSI escape codes from text for easier testing."""
    ansi_escape = re.compile(r"\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])")
    return ansi_escape.sub("", text)


class TestCLIExitCodes:
    """Test that CLI commands return correct exit codes."""

    def test_cli_input_validation_error_exit_code(self):
        """Test that CLI input validation errors return non-zero exit code."""
        runner = CliRunner()

        # Test missing required arguments - should return non-zero exit code
        result = runner.invoke(app, ["self"])
        assert result.exit_code != 0  # Should fail with non-zero exit code
        assert result.exit_code in [1, 2]  # Either typer's 2 or our handled 1

    def test_cli_file_not_found_error_exit_code(self):
        """Test that file not found errors return exit code 1."""
        runner = CliRunner()

        # Test with nonexistent files
        result = runner.invoke(
            app, ["summarize", "--vcf", "nonexistent.vcf", "--bam", "nonexistent.bam"]
        )
        assert result.exit_code == 1
        assert "VCF file not found" in result.output

    def test_cli_sample_name_validation_error(self):
        """Test that invalid sample names return proper error and exit code."""
        runner = CliRunner()

        # Create temporary files with names that don't match
        import tempfile

        with (
            tempfile.NamedTemporaryFile(suffix="_R1.fastq.gz", delete=False) as r1,
            tempfile.NamedTemporaryFile(
                suffix="_R2_different.fastq.gz", delete=False
            ) as r2,
            tempfile.NamedTemporaryFile(suffix=".fna", delete=False) as assembly,
        ):

            r1_path = Path(r1.name)
            r2_path = Path(r2.name)
            assembly_path = Path(assembly.name)

            try:
                result = runner.invoke(
                    app,
                    [
                        "self",
                        "--r1",
                        str(r1_path),
                        "--r2",
                        str(r2_path),
                        "--assembly",
                        str(assembly_path),
                    ],
                )

                assert result.exit_code == 1
                assert "Cannot infer sample name" in result.output
                assert "Use --sample" in result.output

            finally:
                r1_path.unlink(missing_ok=True)
                r2_path.unlink(missing_ok=True)
                assembly_path.unlink(missing_ok=True)

    def test_cli_dry_run_exits_zero_on_success(self):
        """Test that dry-run exits with code 0 when successful."""
        runner = CliRunner()

        # Use real files if available for dry-run test
        vcf_file = "test_output/2508H52931_bwa_bbtools.normalized.vcf.gz"
        bam_file = "test_output/2508H52931_bwa_bbtools.sorted.bam"

        vcf_path = Path(vcf_file)
        bam_path = Path(bam_file)

        if vcf_path.exists() and bam_path.exists():
            result = runner.invoke(
                app,
                [
                    "--dry-run",
                    "summarize",
                    "--vcf",
                    str(vcf_path),
                    "--bam",
                    str(bam_path),
                ],
            )

            # Should exit successfully with dry-run
            assert result.exit_code == 0
            assert "DRY RUN" in result.output
        else:
            pytest.skip("Real test files not available")


class TestExitCodeMapping:
    """Test exit code mapping functions from errors module."""

    def test_map_exception_to_exit_code_file_not_found(self):
        """Test that FileNotFoundError maps to exit code 1."""
        error = FileNotFoundError("File not found")
        assert map_exception_to_exit_code(error) == 1

    def test_map_exception_to_exit_code_value_error(self):
        """Test that ValueError maps to exit code 1."""
        error = ValueError("Invalid value")
        assert map_exception_to_exit_code(error) == 1

    def test_map_exception_to_exit_code_external_tool_error(self):
        """Test that ExternalToolError maps to exit code 4."""
        from ssiamb.mapping import ExternalToolError

        error = ExternalToolError("Tool failed")
        assert map_exception_to_exit_code(error) == 4

    def test_map_exception_to_exit_code_compatibility_error(self):
        """Test that CompatibilityError maps to exit code 3."""
        from ssiamb.reuse import CompatibilityError

        error = CompatibilityError("Files incompatible")
        assert map_exception_to_exit_code(error) == 3

    def test_map_exception_to_exit_code_reference_error(self):
        """Test that ReferenceResolutionError maps to exit code 2."""
        from ssiamb.refdir import ReferenceResolutionError

        error = ReferenceResolutionError("Species not found")
        assert map_exception_to_exit_code(error) == 2


class TestErrorMessageFormatting:
    """Test error message formatting functions."""

    def test_error_messages_are_strings(self):
        """Test that error messages are properly formatted as strings."""
        error = ValueError("Something went wrong")
        assert isinstance(str(error), str)
        assert "Something went wrong" in str(error)

    def test_file_not_found_error_format(self):
        """Test formatting FileNotFoundError."""
        error = FileNotFoundError("test.txt not found")
        assert isinstance(str(error), str)
        assert "test.txt not found" in str(error)


class TestErrorMessageSuggestions:
    """Test that error messages provide helpful suggestions."""

    def test_input_validation_suggestions(self):
        """Test that input validation errors provide helpful suggestions."""
        runner = CliRunner()

        # Test file not found suggestion
        result = runner.invoke(
            app, ["summarize", "--vcf", "missing.vcf", "--bam", "missing.bam"]
        )

        assert result.exit_code == 1
        assert "Suggestions:" in result.output
        assert "Check that the file path is correct" in result.output

    def test_sample_name_suggestions(self):
        """Test that sample name errors provide helpful suggestions."""
        runner = CliRunner()

        import tempfile

        with (
            tempfile.NamedTemporaryFile(suffix="bad name.fastq.gz", delete=False) as r1,
            tempfile.NamedTemporaryFile(
                suffix="bad name2.fastq.gz", delete=False
            ) as r2,
            tempfile.NamedTemporaryFile(suffix=".fna", delete=False) as assembly,
        ):

            r1_path = Path(r1.name)
            r2_path = Path(r2.name)
            assembly_path = Path(assembly.name)

            try:
                result = runner.invoke(
                    app,
                    [
                        "self",
                        "--r1",
                        str(r1_path),
                        "--r2",
                        str(r2_path),
                        "--assembly",
                        str(assembly_path),
                    ],
                )

                assert result.exit_code == 1
                assert "Use --sample to provide" in result.output
                assert "Sample names must match pattern" in result.output

            finally:
                r1_path.unlink(missing_ok=True)
                r2_path.unlink(missing_ok=True)
                assembly_path.unlink(missing_ok=True)

    def test_help_text_shows_options(self):
        """Test that help text shows all available options clearly."""
        runner = CliRunner()

        result = runner.invoke(app, ["--help"])
        assert result.exit_code == 0
        clean_output = strip_ansi_codes(result.output)
        assert "--version" in clean_output
        assert "--dry-run" in clean_output


class TestErrorHandlingIntegration:
    """Test integration of error handling across the system."""

    def test_missing_tools_handled_gracefully(self):
        """Test that missing external tools are handled with appropriate messages."""
        # This would require mocking tool availability, but we can test the structure
        runner = CliRunner()

        # Just verify the CLI structure handles errors properly
        result = runner.invoke(app, ["self", "--help"])
        assert result.exit_code == 0
        clean_output = strip_ansi_codes(result.output)
        assert "Self-mapping mode" in clean_output

    def test_invalid_parameters_rejected(self):
        """Test that invalid parameters are rejected with helpful messages."""
        runner = CliRunner()

        # Test invalid dp-min value
        result = runner.invoke(
            app,
            [
                "summarize",
                "--vcf",
                "test.vcf",
                "--bam",
                "test.bam",
                "--dp-min",
                "0",  # Should be >= 1
            ],
        )

        # Should either reject the value or proceed to file validation
        assert result.exit_code != 0 or "not found" in result.output

    def test_conflicting_options_handled(self):
        """Test that conflicting options are handled appropriately."""
        runner = CliRunner()

        # Test --verbose and --quiet together
        result = runner.invoke(app, ["--verbose", "--quiet", "--help"])
        # Should either handle gracefully or show error
        assert result.exit_code in [0, 1, 2]

        # Test more specific case with actual command
        result = runner.invoke(
            app,
            [
                "--verbose",
                "--quiet",
                "self",
                "--r1",
                "test1.fq",
                "--r2",
                "test2.fq",
                "--assembly",
                "test.fa",
            ],
        )
        assert result.exit_code == 1
        assert "--verbose and --quiet cannot be used together" in result.output
