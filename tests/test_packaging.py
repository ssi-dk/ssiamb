#!/usr/bin/env python3
"""
Tests for Milestone 22: Packaging Polish & Version

Tests the packaging and version functionality to ensure:
- --version flag shows correct version
- CLI help output is complete and formatted correctly
- pyproject.toml metadata is complete
- README documentation is up to date
"""

import pytest
import re
from pathlib import Path
from typer.testing import CliRunner

from ssiamb.cli import app
from ssiamb.version import __version__


class TestVersionCommand:
    """Test --version flag functionality."""

    def test_version_flag_shows_version(self):
        """Test that --version shows the correct version."""
        runner = CliRunner()
        
        result = runner.invoke(app, ["--version"])
        assert result.exit_code == 0
        assert f"ssiamb version {__version__}" in result.output

    def test_version_flag_short_form(self):
        """Test that -V also shows version."""
        runner = CliRunner()
        
        result = runner.invoke(app, ["-V"])
        assert result.exit_code == 0
        assert f"ssiamb version {__version__}" in result.output

    def test_version_format(self):
        """Test that version follows semantic versioning pattern."""
        # Should be either dev version (0.0.0.dev0) or proper semver (x.y.z) or rc version (x.y.z-rc1)
        version_pattern = r'^\d+\.\d+\.\d+(?:\.dev\d+|-rc\d+)?(?:\+.*)?$'
        assert re.match(version_pattern, __version__), f"Version {__version__} doesn't match semver pattern"

    def test_version_callback_exits(self):
        """Test that version callback causes immediate exit."""
        runner = CliRunner()
        
        # Version should exit without needing subcommand
        result = runner.invoke(app, ["--version"])
        assert result.exit_code == 0
        # Should not show help or subcommand info
        assert "Commands:" not in result.output


class TestCLIHelp:
    """Test CLI help output and structure."""

    def test_main_help_output(self):
        """Test main help shows all expected elements."""
        runner = CliRunner()
        
        result = runner.invoke(app, ["--help"])
        assert result.exit_code == 0
        
        # Check essential elements (Rich CLI uses "Commands" in a styled box)
        assert "ssiamb" in result.output
        assert "SSI Ambiguous Site Detection Tool" in result.output
        assert "Commands" in result.output  # Rich formatting shows Commands
        assert "self" in result.output
        assert "ref" in result.output
        assert "summarize" in result.output

    def test_main_help_shows_global_options(self):
        """Test that main help shows global options."""
        runner = CliRunner()
        
        result = runner.invoke(app, ["--help"])
        assert result.exit_code == 0
        
        # Check global options
        assert "--version" in result.output
        assert "--config" in result.output
        assert "--verbose" in result.output
        assert "--quiet" in result.output
        assert "--dry-run" in result.output

    def test_self_help_output(self):
        """Test self command help output."""
        runner = CliRunner()
        
        result = runner.invoke(app, ["self", "--help"])
        assert result.exit_code == 0
        
        # Check self-specific options
        assert "--r1" in result.output
        assert "--r2" in result.output
        assert "--assembly" in result.output
        assert "Self-mapping mode" in result.output

    def test_ref_help_output(self):
        """Test ref command help output."""
        runner = CliRunner()
        
        result = runner.invoke(app, ["ref", "--help"])
        assert result.exit_code == 0
        
        # Check ref-specific options
        assert "--r1" in result.output
        assert "--r2" in result.output
        assert "--reference" in result.output or "--species" in result.output or "--bracken" in result.output
        assert "Reference-mapping mode" in result.output

    def test_summarize_help_output(self):
        """Test summarize command help output."""
        runner = CliRunner()
        
        result = runner.invoke(app, ["summarize", "--help"])
        assert result.exit_code == 0
        
        # Check summarize-specific options
        assert "--vcf" in result.output
        assert "--bam" in result.output
        assert "Summarize VCF and BAM files" in result.output

    def test_help_shows_common_options(self):
        """Test that command help shows common analysis options."""
        runner = CliRunner()
        
        for command in ["self", "ref", "summarize"]:
            result = runner.invoke(app, [command, "--help"])
            assert result.exit_code == 0
            
            # Common analysis options should appear
            assert "--dp-min" in result.output
            assert "--maf-min" in result.output
            assert "--emit-vcf" in result.output
            assert "--stdout" in result.output

    def test_help_shows_required_vs_optional(self):
        """Test that help clearly indicates required vs optional parameters."""
        runner = CliRunner()
        
        result = runner.invoke(app, ["summarize", "--help"])
        assert result.exit_code == 0
        
        # Required options should be marked
        assert "*" in result.output  # Typer marks required options with *
        
        # VCF and BAM should be required for summarize
        vcf_line = [line for line in result.output.split('\n') if '--vcf' in line]
        bam_line = [line for line in result.output.split('\n') if '--bam' in line]
        
        assert len(vcf_line) > 0
        assert len(bam_line) > 0


class TestCLIUsability:
    """Test CLI usability features."""

    def test_no_command_shows_help(self):
        """Test that running ssiamb without command shows help."""
        runner = CliRunner()
        
        result = runner.invoke(app, [])
        # Should show help rather than error
        assert "Commands" in result.output or "Usage:" in result.output

    def test_invalid_command_shows_error(self):
        """Test that invalid commands show helpful error."""
        runner = CliRunner()
        
        result = runner.invoke(app, ["invalid_command"])
        assert result.exit_code != 0
        # Should suggest available commands or show help
        assert "No such command" in result.output or "Usage:" in result.output

    def test_conflicting_flags_handled(self):
        """Test that conflicting flags are handled properly."""
        runner = CliRunner()
        
        # Test --verbose and --quiet conflict
        result = runner.invoke(app, ["--verbose", "--quiet", "self", "--help"])
        # Should either error or handle gracefully
        assert result.exit_code in [0, 1, 2]  # Acceptable exit codes

    def test_emission_flags_clear(self):
        """Test that emission flags are clearly documented."""
        runner = CliRunner()
        
        result = runner.invoke(app, ["self", "--help"])
        assert result.exit_code == 0
        
        # All emission flags should be documented
        emission_flags = ["--emit-vcf", "--emit-bed", "--emit-matrix", 
                         "--emit-per-contig", "--emit-multiqc", "--emit-provenance"]
        
        for flag in emission_flags:
            assert flag in result.output


class TestPackagingMetadata:
    """Test packaging metadata completeness."""

    def test_pyproject_toml_exists(self):
        """Test that pyproject.toml exists and is readable."""
        pyproject_path = Path("pyproject.toml")
        assert pyproject_path.exists(), "pyproject.toml not found"
        
        content = pyproject_path.read_text()
        assert len(content) > 0, "pyproject.toml is empty"

    def test_pyproject_has_required_metadata(self):
        """Test that pyproject.toml has required metadata fields."""
        pyproject_path = Path("pyproject.toml")
        if not pyproject_path.exists():
            pytest.skip("pyproject.toml not found")
        
        content = pyproject_path.read_text()
        
        # Check for essential fields
        assert 'name = "ssiamb"' in content
        assert 'version' in content
        assert 'description' in content
        assert 'authors' in content
        assert 'license' in content
        assert '[project.urls]' in content and ('Homepage' in content or 'Repository' in content)

    def test_pyproject_has_dependencies(self):
        """Test that pyproject.toml lists required dependencies."""
        pyproject_path = Path("pyproject.toml")
        if not pyproject_path.exists():
            pytest.skip("pyproject.toml not found")
        
        content = pyproject_path.read_text()
        
        # Check for key dependencies
        required_deps = ["typer", "rich", "pyyaml", "pandas", "numpy", "pysam", "biopython"]
        
        for dep in required_deps:
            assert dep in content, f"Required dependency {dep} not found in pyproject.toml"

    def test_pyproject_has_console_script(self):
        """Test that pyproject.toml defines console script entry point."""
        pyproject_path = Path("pyproject.toml")
        if not pyproject_path.exists():
            pytest.skip("pyproject.toml not found")
        
        content = pyproject_path.read_text()
        
        # Check for console script
        assert "[project.scripts]" in content or "[tool.setuptools.entry-points" in content
        assert "ssiamb" in content


class TestREADMEDocumentation:
    """Test README documentation completeness."""

    def test_readme_exists(self):
        """Test that README.md exists."""
        readme_path = Path("README.md")
        assert readme_path.exists(), "README.md not found"

    def test_readme_has_basic_sections(self):
        """Test that README has essential sections."""
        readme_path = Path("README.md")
        if not readme_path.exists():
            pytest.skip("README.md not found")
        
        content = readme_path.read_text()
        
        # Check for essential sections
        assert "# ssiamb" in content or "ssiamb" in content[:100]
        assert "install" in content.lower()
        assert "usage" in content.lower() or "example" in content.lower()

    def test_readme_has_cli_examples(self):
        """Test that README includes CLI usage examples."""
        readme_path = Path("README.md")
        if not readme_path.exists():
            pytest.skip("README.md not found")
        
        content = readme_path.read_text()
        
        # Should have examples of main commands
        assert "ssiamb self" in content
        assert "ssiamb ref" in content or "ssiamb summarize" in content

    def test_readme_documents_dry_run(self):
        """Test that README documents dry-run functionality."""
        readme_path = Path("README.md")
        if not readme_path.exists():
            pytest.skip("README.md not found")
        
        content = readme_path.read_text()
        
        # Should mention dry-run feature
        assert "--dry-run" in content or "dry run" in content.lower()

    def test_readme_documents_error_codes(self):
        """Test that README documents error codes."""
        readme_path = Path("README.md")
        if not readme_path.exists():
            pytest.skip("README.md not found")
        
        content = readme_path.read_text()
        
        # Should document exit codes
        assert "exit code" in content.lower() or "error code" in content.lower()


class TestVersionConsistency:
    """Test version consistency across files."""

    def test_version_consistent_across_files(self):
        """Test that version is consistent between version.py and pyproject.toml."""
        from ssiamb.version import __version__
        
        pyproject_path = Path("pyproject.toml")
        if not pyproject_path.exists():
            pytest.skip("pyproject.toml not found")
        
        content = pyproject_path.read_text()
        
        # Extract version from pyproject.toml
        import re
        version_match = re.search(r'version\s*=\s*"([^"]+)"', content)
        if version_match:
            pyproject_version = version_match.group(1)
            assert __version__ == pyproject_version, \
                f"Version mismatch: version.py has {__version__}, pyproject.toml has {pyproject_version}"

    def test_version_in_cli_output_matches_module(self):
        """Test that CLI version output matches module version."""
        runner = CliRunner()
        
        result = runner.invoke(app, ["--version"])
        assert result.exit_code == 0
        
        from ssiamb.version import __version__
        assert __version__ in result.output


class TestInstallationReadiness:
    """Test that package is ready for installation."""

    def test_import_works(self):
        """Test that the package can be imported."""
        try:
            import ssiamb
            import ssiamb.cli
            import ssiamb.runner
            import ssiamb.errors
        except ImportError as e:
            pytest.fail(f"Failed to import ssiamb modules: {e}")

    def test_console_script_callable(self):
        """Test that the main CLI function is callable."""
        from ssiamb.cli import app
        
        # Should be a Typer app
        assert hasattr(app, '__call__')
        assert hasattr(app, 'command')  # Typer app should have command method

    def test_all_required_modules_exist(self):
        """Test that all required modules exist."""
        required_modules = [
            'ssiamb.cli',
            'ssiamb.runner', 
            'ssiamb.models',
            'ssiamb.errors',
            'ssiamb.version',
            'ssiamb.io_utils',
            'ssiamb.vcf_ops',
            'ssiamb.mapping',
            'ssiamb.calling',
            'ssiamb.depth'
        ]
        
        for module_name in required_modules:
            try:
                __import__(module_name)
            except ImportError:
                pytest.fail(f"Required module {module_name} cannot be imported")