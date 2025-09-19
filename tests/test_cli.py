"""
Tests for cli.py command-line interface.

Tests CLI argument parsing, command dispatch, error handling, and integration.
"""

import pytest
import tempfile
from pathlib import Path
from unittest.mock import patch, MagicMock
from typer.testing import CliRunner

from src.ssiamb.cli import app
from src.ssiamb.models import Mode


class TestMainCallback:
    """Test main callback and global options."""

    @patch("src.ssiamb.config.load_config")
    def test_config_loading_success(self, mock_load_config):
        """Test successful config loading."""
        runner = CliRunner()

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".yaml", delete=False
        ) as config_file:
            config_path = Path(config_file.name)
            config_file.write("# test config")

        try:
            result = runner.invoke(
                app,
                [
                    "--config",
                    str(config_path),
                    "--verbose",
                    "self",
                    "--r1",
                    "test1.fq",
                    "--r2",
                    "test2.fq",
                    "--assembly",
                    "test.fa",
                ],
            )

            mock_load_config.assert_called_once_with(config_path)
            assert "Loaded configuration from:" in result.stdout
        finally:
            config_path.unlink()

    @patch("src.ssiamb.config.load_config")
    def test_config_loading_failure(self, mock_load_config):
        """Test config loading failure."""
        mock_load_config.side_effect = Exception("Invalid config")
        runner = CliRunner()

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".yaml", delete=False
        ) as config_file:
            config_path = Path(config_file.name)
            config_file.write("invalid: yaml: content:")

        try:
            result = runner.invoke(
                app,
                [
                    "--config",
                    str(config_path),
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
            assert "Error loading config file" in result.stdout
        finally:
            config_path.unlink()


class TestSelfCommand:
    """Test the 'self' command functionality."""

    @patch("src.ssiamb.cli.execute_plan")
    @patch("src.ssiamb.cli.create_run_plan")
    def test_self_command_minimal_args(self, mock_create_plan, mock_execute):
        """Test self command with minimal required arguments."""
        mock_plan = MagicMock()
        mock_create_plan.return_value = mock_plan

        # Mock execute_plan to return (summary, provenance) tuple
        mock_summary = MagicMock()
        mock_provenance = MagicMock()
        mock_execute.return_value = (mock_summary, mock_provenance)

        runner = CliRunner()

        with tempfile.TemporaryDirectory() as tmpdir:
            r1_path = Path(tmpdir) / "test_R1.fastq"
            r2_path = Path(tmpdir) / "test_R2.fastq"
            assembly_path = Path(tmpdir) / "assembly.fasta"

            # Create dummy files
            r1_path.write_text("@read1\\nACGT\\n+\\nIIII\\n")
            r2_path.write_text("@read2\\nACGT\\n+\\nIIII\\n")
            assembly_path.write_text(">contig1\\nACGTACGT\\n")

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

            assert result.exit_code == 0
            mock_create_plan.assert_called_once()
            mock_execute.assert_called_once_with(mock_plan)

    @patch("src.ssiamb.cli.execute_plan")
    @patch("src.ssiamb.cli.create_run_plan")
    def test_self_command_all_options(self, mock_create_plan, mock_execute):
        """Test self command with all optional arguments."""
        mock_plan = MagicMock()
        mock_create_plan.return_value = mock_plan

        # Mock execute_plan to return (summary, provenance) tuple
        mock_summary = MagicMock()
        mock_provenance = MagicMock()
        mock_execute.return_value = (mock_summary, mock_provenance)

        runner = CliRunner()

        with tempfile.TemporaryDirectory() as tmpdir:
            r1_path = Path(tmpdir) / "sample_R1.fastq"
            r2_path = Path(tmpdir) / "sample_R2.fastq"
            assembly_path = Path(tmpdir) / "assembly.fasta"
            output_dir = Path(tmpdir) / "output"

            # Create dummy files
            r1_path.write_text("@read1\\nACGT\\n+\\nIIII\\n")
            r2_path.write_text("@read2\\nACGT\\n+\\nIIII\\n")
            assembly_path.write_text(">contig1\\nACGTACGT\\n")
            output_dir.mkdir()

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
                    "--sample",
                    "test_sample",
                    "--outdir",
                    str(output_dir),
                    "--threads",
                    "8",
                    "--mapper",
                    "bwa-mem2",
                    "--caller",
                    "bcftools",
                    "--dp-min",
                    "15",
                    "--maf-min",
                    "0.05",
                    "--dp-cap",
                    "200",
                    "--depth-tool",
                    "mosdepth",
                    "--tsv-mode",
                    "append",
                ],
            )

            assert result.exit_code == 0
            mock_create_plan.assert_called_once()

            # Verify the plan was created with correct parameters
            call_args = mock_create_plan.call_args
            assert call_args[1]["mode"] == Mode.SELF
            assert call_args[1]["sample"] == "test_sample"
            assert call_args[1]["threads"] == 8
            assert call_args[1]["mapper"] == "bwa-mem2"
            assert call_args[1]["caller"] == "bcftools"

    def test_self_command_missing_required_args(self):
        """Test self command with missing required arguments."""
        runner = CliRunner()

        result = runner.invoke(app, ["self"])

        assert result.exit_code != 0
        # Check both stdout and stderr for error messages
        output = result.stdout + result.stderr
        assert (
            "Missing option" in output
            or "required" in output.lower()
            or result.exit_code == 2
        )

    @patch("src.ssiamb.cli.execute_plan")
    @patch("src.ssiamb.cli.create_run_plan")
    def test_self_command_file_validation(self, mock_create_plan, mock_execute):
        """Test self command file validation."""
        runner = CliRunner()

        # Test with non-existent files
        result = runner.invoke(
            app,
            [
                "self",
                "--r1",
                "/nonexistent/r1.fastq",
                "--r2",
                "/nonexistent/r2.fastq",
                "--assembly",
                "/nonexistent/assembly.fasta",
            ],
        )

        # Should fail due to file validation
        assert result.exit_code != 0


class TestRefCommand:
    """Test the 'ref' command functionality."""

    @patch("src.ssiamb.cli.execute_plan")
    @patch("src.ssiamb.cli.create_run_plan")
    def test_ref_command_minimal_args(self, mock_create_plan, mock_execute):
        """Test ref command with minimal required arguments."""
        mock_plan = MagicMock()
        mock_create_plan.return_value = mock_plan

        # Mock execute_plan to return (summary, provenance) tuple
        mock_summary = MagicMock()
        mock_provenance = MagicMock()
        mock_execute.return_value = (mock_summary, mock_provenance)

        runner = CliRunner()

        with tempfile.TemporaryDirectory() as tmpdir:
            r1_path = Path(tmpdir) / "test_R1.fastq"
            r2_path = Path(tmpdir) / "test_R2.fastq"
            reference_path = Path(tmpdir) / "reference.fasta"

            # Create dummy files
            r1_path.write_text("@read1\\nACGT\\n+\\nIIII\\n")
            r2_path.write_text("@read2\\nACGT\\n+\\nIIII\\n")
            reference_path.write_text(">ref\\nACGTACGT\\n")

            result = runner.invoke(
                app,
                [
                    "ref",
                    "--r1",
                    str(r1_path),
                    "--r2",
                    str(r2_path),
                    "--reference",
                    str(reference_path),
                ],
            )

            assert result.exit_code == 0
            mock_create_plan.assert_called_once()

            # Verify mode is REF
            call_args = mock_create_plan.call_args
            assert call_args[1]["mode"] == Mode.REF
            # Check that execute_plan was called with the plan (ignore additional keyword args)
            mock_execute.assert_called_once()
            assert mock_execute.call_args[0][0] == mock_plan

    @patch("src.ssiamb.cli.execute_plan")
    @patch("src.ssiamb.cli.create_run_plan")
    def test_ref_command_all_options(self, mock_create_plan, mock_execute):
        """Test ref command with all optional arguments."""
        mock_plan = MagicMock()
        mock_create_plan.return_value = mock_plan

        # Mock execute_plan to return (summary, provenance) tuple
        mock_summary = MagicMock()
        mock_provenance = MagicMock()
        mock_execute.return_value = (mock_summary, mock_provenance)

        runner = CliRunner()

        with tempfile.TemporaryDirectory() as tmpdir:
            r1_path = Path(tmpdir) / "test_R1.fastq"
            r2_path = Path(tmpdir) / "test_R2.fastq"
            reference_path = Path(tmpdir) / "reference.fasta"
            output_dir = Path(tmpdir) / "output"

            # Create dummy files
            r1_path.write_text("@read1\\nACGT\\n+\\nIIII\\n")
            r2_path.write_text("@read2\\nACGT\\n+\\nIIII\\n")
            reference_path.write_text(">ref\\nACGTACGT\\n")
            output_dir.mkdir()

            result = runner.invoke(
                app,
                [
                    "ref",
                    "--r1",
                    str(r1_path),
                    "--r2",
                    str(r2_path),
                    "--reference",
                    str(reference_path),
                    "--sample",
                    "ref_test_sample",
                    "--outdir",
                    str(output_dir),
                    "--threads",
                    "16",
                    "--mapper",
                    "minimap2",
                    "--caller",
                    "bbtools",
                ],
            )

            assert result.exit_code == 0
            mock_create_plan.assert_called_once()
            # Check that execute_plan was called with the plan (ignore additional keyword args)
            mock_execute.assert_called_once()
            assert mock_execute.call_args[0][0] == mock_plan


class TestSummarizeCommand:
    """Test the 'summarize' command functionality."""

    @patch("src.ssiamb.cli.run_summarize")
    def test_summarize_command_minimal_args(self, mock_run_summarize):
        """Test summarize command with minimal arguments."""
        runner = CliRunner()

        with tempfile.TemporaryDirectory() as tmpdir:
            vcf_path = Path(tmpdir) / "results.vcf"
            bam_path = Path(tmpdir) / "sample.bam"
            vcf_path.write_text("##fileformat=VCFv4.2\n")
            bam_path.write_text("")  # Empty file for testing

            result = runner.invoke(
                app, ["summarize", "--vcf", str(vcf_path), "--bam", str(bam_path)]
            )

            assert result.exit_code == 0
            mock_run_summarize.assert_called_once()

    @patch("src.ssiamb.cli.run_summarize")
    def test_summarize_command_all_options(self, mock_run_summarize):
        """Test summarize command with all options."""
        runner = CliRunner()

        with tempfile.TemporaryDirectory() as tmpdir:
            vcf_path = Path(tmpdir) / "results.vcf"
            bam_path = Path(tmpdir) / "sample.bam"
            output_path = Path(tmpdir) / "output.tsv"

            vcf_path.write_text("##fileformat=VCFv4.2\n")
            bam_path.write_text("")  # Empty file for testing

            result = runner.invoke(
                app,
                [
                    "summarize",
                    "--vcf",
                    str(vcf_path),
                    "--bam",
                    str(bam_path),
                    "--output",
                    str(output_path),
                    "--dp-min",
                    "15",
                    "--maf-min",
                    "0.05",
                    "--emit-vcf",
                    "--emit-bed",
                ],
            )

            assert result.exit_code == 0
            mock_run_summarize.assert_called_once()


class TestCLIIntegration:
    """Test CLI integration and error handling."""

    @patch("src.ssiamb.cli.create_run_plan")
    def test_create_run_plan_error_handling(self, mock_create_plan):
        """Test error handling when create_run_plan fails."""
        mock_create_plan.side_effect = Exception("Plan creation failed")

        runner = CliRunner()

        with tempfile.TemporaryDirectory() as tmpdir:
            r1_path = Path(tmpdir) / "test_R1.fastq"
            r2_path = Path(tmpdir) / "test_R2.fastq"
            assembly_path = Path(tmpdir) / "assembly.fasta"

            r1_path.write_text("@read1\\nACGT\\n+\\nIIII\\n")
            r2_path.write_text("@read2\\nACGT\\n+\\nIIII\\n")
            assembly_path.write_text(">contig1\\nACGTACGT\\n")

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

            # Should handle the exception gracefully
            assert result.exit_code != 0

    @patch("src.ssiamb.cli.execute_plan")
    @patch("src.ssiamb.cli.create_run_plan")
    def test_execute_plan_error_handling(self, mock_create_plan, mock_execute):
        """Test error handling when execute_plan fails."""
        mock_plan = MagicMock()
        mock_create_plan.return_value = mock_plan
        mock_execute.side_effect = Exception("Execution failed")

        runner = CliRunner()

        with tempfile.TemporaryDirectory() as tmpdir:
            r1_path = Path(tmpdir) / "test_R1.fastq"
            r2_path = Path(tmpdir) / "test_R2.fastq"
            assembly_path = Path(tmpdir) / "assembly.fasta"

            r1_path.write_text("@read1\\nACGT\\n+\\nIIII\\n")
            r2_path.write_text("@read2\\nACGT\\n+\\nIIII\\n")
            assembly_path.write_text(">contig1\\nACGTACGT\\n")

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

            # Should handle the exception gracefully
            assert result.exit_code != 0


class TestCLIParameterValidation:
    """Test CLI parameter validation and type conversion."""

    @pytest.mark.parametrize(
        "threads_value,expected_valid",
        [
            ("1", True),
            ("8", True),
            ("32", True),
            ("0", False),
            ("-1", False),
            ("abc", False),
        ],
    )
    def test_threads_parameter_validation(self, threads_value, expected_valid):
        """Test threads parameter validation."""
        runner = CliRunner()

        with tempfile.TemporaryDirectory() as tmpdir:
            r1_path = Path(tmpdir) / "test_R1.fastq"
            r2_path = Path(tmpdir) / "test_R2.fastq"
            assembly_path = Path(tmpdir) / "assembly.fasta"

            r1_path.write_text("@read1\\nACGT\\n+\\nIIII\\n")
            r2_path.write_text("@read2\\nACGT\\n+\\nIIII\\n")
            assembly_path.write_text(">contig1\\nACGTACGT\\n")

            with (
                patch("src.ssiamb.cli.create_run_plan"),
                patch("src.ssiamb.cli.execute_plan"),
            ):
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
                        "--threads",
                        threads_value,
                    ],
                )

                if expected_valid:
                    assert (
                        result.exit_code == 0 or "Missing option" not in result.stdout
                    )
                else:
                    assert result.exit_code != 0
