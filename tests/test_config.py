"""
Tests for config.py configuration management.

Tests configuration loading, merging, validation, environment variables,
and global config management.
"""

import pytest
import tempfile
import os
import yaml
from pathlib import Path
from unittest.mock import patch

from src.ssiamb.config import SsiambConfig, get_config, set_config, load_config


class TestSsiambConfig:
    """Test SsiambConfig class functionality."""

    def test_load_defaults_when_no_file(self):
        """Test loading defaults when defaults.yaml doesn't exist."""
        with patch("src.ssiamb.config.Path.exists", return_value=False):
            config = SsiambConfig._load_defaults()

            assert "thresholds" in config
            assert config["thresholds"]["dp_min"] == 10
            assert config["thresholds"]["maf_min"] == 0.1
            assert config["thresholds"]["dp_cap"] == 100
            assert config["thresholds"]["mapq_min"] == 20
            assert config["thresholds"]["baseq_min"] == 20
            assert config["species_aliases"] == {}
            assert config["tools"] == {}
            assert config["output"] == {}

    def test_load_defaults_from_file(self):
        """Test loading defaults from existing defaults.yaml file."""
        defaults_content = {
            "thresholds": {"dp_min": 15, "maf_min": 0.05, "dp_cap": 150},
            "species_aliases": {"ecoli": "Escherichia coli"},
            "tools": {"minimap2": {"threads": 8}},
        }

        with (
            patch("src.ssiamb.config.Path.exists", return_value=True),
            patch(
                "src.ssiamb.config.SsiambConfig._load_yaml",
                return_value=defaults_content,
            ),
        ):
            config = SsiambConfig._load_defaults()

            assert config["thresholds"]["dp_min"] == 15
            assert config["thresholds"]["maf_min"] == 0.05
            assert config["species_aliases"]["ecoli"] == "Escherichia coli"
            assert config["tools"]["minimap2"]["threads"] == 8

    def test_load_yaml_success(self):
        """Test successful YAML loading."""
        test_config = {"thresholds": {"dp_min": 20}, "tools": {"bwa": {"mem": True}}}

        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
            yaml.dump(test_config, f)
            config_path = Path(f.name)

        try:
            result = SsiambConfig._load_yaml(config_path)
            assert result == test_config
        finally:
            config_path.unlink()

    def test_load_yaml_invalid_yaml(self):
        """Test YAML loading with invalid YAML content."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
            f.write("invalid: yaml: content: [")  # Invalid YAML
            config_path = Path(f.name)

        try:
            with pytest.raises(ValueError, match="Invalid YAML"):
                SsiambConfig._load_yaml(config_path)
        finally:
            config_path.unlink()

    def test_load_yaml_file_not_found(self):
        """Test YAML loading with non-existent file."""
        non_existent_path = Path("/tmp/non_existent_config.yaml")

        with pytest.raises(FileNotFoundError, match="Could not read config file"):
            SsiambConfig._load_yaml(non_existent_path)

    def test_load_yaml_empty_file(self):
        """Test YAML loading with empty file."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
            f.write("")  # Empty file
            config_path = Path(f.name)

        try:
            result = SsiambConfig._load_yaml(config_path)
            assert result == {}
        finally:
            config_path.unlink()

    def test_merge_configs_simple(self):
        """Test simple configuration merging."""
        base = {
            "thresholds": {"dp_min": 10, "maf_min": 0.1},
            "tools": {"minimap2": {"version": "2.30"}},
        }

        overlay = {
            "thresholds": {"dp_min": 20},  # Override
            "output": {"format": "tsv"},  # New section
        }

        result = SsiambConfig._merge_configs(base, overlay)

        assert result["thresholds"]["dp_min"] == 20  # Overridden
        assert result["thresholds"]["maf_min"] == 0.1  # Preserved
        assert result["tools"]["minimap2"]["version"] == "2.30"  # Preserved
        assert result["output"]["format"] == "tsv"  # New

    def test_merge_configs_nested(self):
        """Test nested configuration merging."""
        base = {
            "tools": {
                "minimap2": {"threads": 4, "preset": "sr"},
                "bwa": {"algorithm": "mem"},
            }
        }

        overlay = {
            "tools": {
                "minimap2": {"threads": 8},  # Override threads, keep preset
                "samtools": {"sort_mem": "2G"},  # New tool
            }
        }

        result = SsiambConfig._merge_configs(base, overlay)

        assert result["tools"]["minimap2"]["threads"] == 8
        assert result["tools"]["minimap2"]["preset"] == "sr"
        assert result["tools"]["bwa"]["algorithm"] == "mem"
        assert result["tools"]["samtools"]["sort_mem"] == "2G"

    def test_merge_configs_replace_non_dict(self):
        """Test merging when overlay value is not a dict."""
        base = {"thresholds": {"dp_min": 10, "maf_min": 0.1}, "enabled": True}

        overlay = {
            "thresholds": "disabled",  # Replace dict with string
            "enabled": False,  # Replace bool
        }

        result = SsiambConfig._merge_configs(base, overlay)

        assert result["thresholds"] == "disabled"
        assert result["enabled"] is False


class TestEnvironmentOverrides:
    """Test environment variable override functionality."""

    def test_apply_env_overrides_no_env_vars(self):
        """Test applying environment overrides when no env vars are set."""
        config = {"thresholds": {"dp_min": 10, "maf_min": 0.1}}

        result = SsiambConfig._apply_env_overrides(config)

        assert result == config  # Unchanged

    @patch.dict(
        os.environ,
        {"SSIAMB_DP_MIN": "25", "SSIAMB_MAF_MIN": "0.05", "SSIAMB_DP_CAP": "200"},
    )
    def test_apply_env_overrides_with_values(self):
        """Test applying environment overrides with set values."""
        config = {"thresholds": {"dp_min": 10, "maf_min": 0.1, "baseq_min": 20}}

        result = SsiambConfig._apply_env_overrides(config)

        assert result["thresholds"]["dp_min"] == 25  # Overridden
        assert result["thresholds"]["maf_min"] == 0.05  # Overridden
        assert result["thresholds"]["dp_cap"] == 200  # Added
        assert result["thresholds"]["baseq_min"] == 20  # Preserved

    @patch.dict(os.environ, {"SSIAMB_DP_MIN": "invalid_number"})
    def test_apply_env_overrides_invalid_value(self):
        """Test environment override with invalid value."""
        config = {"thresholds": {"dp_min": 10}}

        with pytest.raises(ValueError, match="Invalid value for SSIAMB_DP_MIN"):
            SsiambConfig._apply_env_overrides(config)

    @patch.dict(os.environ, {"SSIAMB_DP_MIN": "15"})
    def test_apply_env_overrides_create_section(self):
        """Test environment override creating missing section."""
        config = {}  # No thresholds section

        result = SsiambConfig._apply_env_overrides(config)

        assert result["thresholds"]["dp_min"] == 15


class TestSpeciesAliasNormalization:
    """Test species alias normalization functionality."""

    @patch("src.ssiamb.refdir.normalize_species_name")
    def test_normalize_species_aliases(self, mock_normalize):
        """Test species alias normalization."""
        mock_normalize.side_effect = lambda x: x.lower().replace(" ", "_")

        aliases = {
            "E. coli": "Escherichia coli",
            "S aureus": "Staphylococcus aureus",
            "listeria": "Listeria monocytogenes",
        }

        result = SsiambConfig._normalize_species_aliases(aliases)

        assert result["e._coli"] == "Escherichia coli"
        assert result["s_aureus"] == "Staphylococcus aureus"
        assert result["listeria"] == "Listeria monocytogenes"


class TestSsiambConfigLoad:
    """Test SsiambConfig.load() method integration."""

    @patch("src.ssiamb.config.SsiambConfig._load_defaults")
    @patch("src.ssiamb.config.SsiambConfig._apply_env_overrides")
    @patch("src.ssiamb.config.SsiambConfig._normalize_species_aliases")
    def test_load_defaults_only(self, mock_normalize, mock_env, mock_defaults):
        """Test loading with defaults only (no user config)."""
        mock_defaults.return_value = {
            "thresholds": {"dp_min": 10},
            "species_aliases": {"ecoli": "E. coli"},
            "tools": {},
            "output": {},
        }
        mock_env.return_value = mock_defaults.return_value
        mock_normalize.return_value = {"ecoli": "E. coli"}

        config = SsiambConfig.load()

        assert config.thresholds == {"dp_min": 10}
        assert config.species_aliases == {"ecoli": "E. coli"}
        assert config.tools == {}
        assert config.output == {}

        mock_defaults.assert_called_once()
        mock_env.assert_called_once()
        mock_normalize.assert_called_once()

    def test_load_with_user_config(self):
        """Test loading with user config file."""
        defaults = {
            "thresholds": {"dp_min": 10, "maf_min": 0.1},
            "species_aliases": {},
            "tools": {},
            "output": {},
        }

        user_config = {
            "thresholds": {"dp_min": 20},
            "tools": {"minimap2": {"threads": 8}},
        }

        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
            yaml.dump(user_config, f)
            config_path = Path(f.name)

        try:
            with (
                patch(
                    "src.ssiamb.config.SsiambConfig._load_defaults",
                    return_value=defaults,
                ),
                patch(
                    "src.ssiamb.config.SsiambConfig._apply_env_overrides",
                    side_effect=lambda x: x,
                ),
                patch(
                    "src.ssiamb.config.SsiambConfig._normalize_species_aliases",
                    side_effect=lambda x: x,
                ),
            ):

                config = SsiambConfig.load(config_path)

                assert config.thresholds["dp_min"] == 20  # Overridden
                assert config.thresholds["maf_min"] == 0.1  # From defaults
                assert config.tools["minimap2"]["threads"] == 8  # From user config
        finally:
            config_path.unlink()

    def test_load_with_nonexistent_user_config(self):
        """Test loading with non-existent user config file."""
        non_existent_path = Path("/tmp/non_existent.yaml")

        with (
            patch("src.ssiamb.config.SsiambConfig._load_defaults") as mock_defaults,
            patch(
                "src.ssiamb.config.SsiambConfig._apply_env_overrides",
                side_effect=lambda x: x,
            ),
            patch(
                "src.ssiamb.config.SsiambConfig._normalize_species_aliases",
                side_effect=lambda x: x,
            ),
        ):

            mock_defaults.return_value = {
                "thresholds": {},
                "species_aliases": {},
                "tools": {},
                "output": {},
            }

            # Should not raise an error, just use defaults
            SsiambConfig.load(non_existent_path)

            mock_defaults.assert_called_once()


class TestConfigGetterMethods:
    """Test configuration getter methods."""

    def test_get_threshold(self):
        """Test get_threshold method."""
        config = SsiambConfig(
            thresholds={"dp_min": 15, "maf_min": 0.05},
            species_aliases={},
            tools={},
            output={},
        )

        assert config.get_threshold("dp_min") == 15
        assert config.get_threshold("maf_min") == 0.05
        assert config.get_threshold("nonexistent") is None
        assert config.get_threshold("nonexistent", "default") == "default"

    def test_get_species_alias(self):
        """Test get_species_alias method."""
        config = SsiambConfig(
            thresholds={},
            species_aliases={
                "ecoli": "Escherichia coli",
                "staph": "Staphylococcus aureus",
            },
            tools={},
            output={},
        )

        assert config.get_species_alias("ecoli") == "Escherichia coli"
        assert config.get_species_alias("staph") == "Staphylococcus aureus"
        assert config.get_species_alias("unknown") == "unknown"  # Returns original

    def test_get_tool_setting(self):
        """Test get_tool_setting method."""
        config = SsiambConfig(
            thresholds={},
            species_aliases={},
            tools={
                "minimap2": {"threads": 8, "preset": "sr"},
                "bwa": {"algorithm": "mem"},
            },
            output={},
        )

        assert config.get_tool_setting("minimap2", "threads") == 8
        assert config.get_tool_setting("minimap2", "preset") == "sr"
        assert config.get_tool_setting("bwa", "algorithm") == "mem"
        assert config.get_tool_setting("minimap2", "nonexistent") is None
        assert config.get_tool_setting("unknown_tool", "setting") is None
        assert (
            config.get_tool_setting("minimap2", "nonexistent", "default") == "default"
        )

    def test_get_output_setting(self):
        """Test get_output_setting method."""
        config = SsiambConfig(
            thresholds={},
            species_aliases={},
            tools={},
            output={"format": "tsv", "precision": 3},
        )

        assert config.get_output_setting("format") == "tsv"
        assert config.get_output_setting("precision") == 3
        assert config.get_output_setting("nonexistent") is None
        assert config.get_output_setting("nonexistent", "default") == "default"


class TestGlobalConfigManagement:
    """Test global configuration management functions."""

    def teardown_method(self):
        """Reset global config after each test."""
        global _config
        import src.ssiamb.config

        src.ssiamb.config._config = None

    @patch("src.ssiamb.config.SsiambConfig.load")
    def test_get_config_lazy_loading(self, mock_load):
        """Test get_config lazy loading."""
        mock_config = SsiambConfig(
            thresholds={}, species_aliases={}, tools={}, output={}
        )
        mock_load.return_value = mock_config

        # First call should trigger loading
        config1 = get_config()
        assert config1 is mock_config
        mock_load.assert_called_once()

        # Second call should return cached instance
        config2 = get_config()
        assert config2 is mock_config
        assert config2 is config1
        mock_load.assert_called_once()  # Still only called once

    def test_set_config(self):
        """Test set_config function."""
        test_config = SsiambConfig(
            thresholds={"dp_min": 25}, species_aliases={}, tools={}, output={}
        )

        set_config(test_config)

        retrieved_config = get_config()
        assert retrieved_config is test_config
        assert retrieved_config.thresholds["dp_min"] == 25

    @patch("src.ssiamb.config.SsiambConfig.load")
    def test_load_config(self, mock_load):
        """Test load_config function."""
        test_config = SsiambConfig(
            thresholds={"dp_min": 30}, species_aliases={}, tools={}, output={}
        )
        mock_load.return_value = test_config

        config_path = Path("/tmp/test.yaml")
        result = load_config(config_path)

        assert result is test_config
        mock_load.assert_called_once_with(config_path)

        # Should be set as global config
        assert get_config() is test_config


class TestConfigIntegration:
    """Test configuration integration scenarios."""

    def test_full_config_loading_workflow(self):
        """Test complete configuration loading workflow."""
        # Create a realistic user config
        user_config_content = {
            "thresholds": {"dp_min": 20, "maf_min": 0.02},
            "species_aliases": {
                "e_coli": "Escherichia coli",
                "s_aureus": "Staphylococcus aureus",
            },
            "tools": {
                "minimap2": {"threads": 16, "preset": "sr"},
                "bbtools": {"min_depth": 5},
            },
            "output": {"format": "tsv", "include_headers": True},
        }

        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
            yaml.dump(user_config_content, f)
            config_path = Path(f.name)

        try:
            with (
                patch("src.ssiamb.config.SsiambConfig._load_defaults") as mock_defaults,
                patch.dict(os.environ, {"SSIAMB_DP_MIN": "25"}),
                patch(
                    "src.ssiamb.refdir.normalize_species_name",
                    side_effect=lambda x: x.lower(),
                ),
            ):

                # Mock reasonable defaults
                mock_defaults.return_value = {
                    "thresholds": {"dp_min": 10, "maf_min": 0.1, "baseq_min": 20},
                    "species_aliases": {},
                    "tools": {"samtools": {"sort_mem": "1G"}},
                    "output": {"precision": 2},
                }

                config = SsiambConfig.load(config_path)

                # Check that user config overrode defaults
                assert config.get_threshold("maf_min") == 0.02  # User override
                assert config.get_threshold("baseq_min") == 20  # From defaults

                # Check that environment variable overrode user config
                assert config.get_threshold("dp_min") == 25  # Environment override

                # Check species aliases
                assert config.get_species_alias("e_coli") == "Escherichia coli"

                # Check tool settings
                assert config.get_tool_setting("minimap2", "threads") == 16
                assert config.get_tool_setting("minimap2", "preset") == "sr"
                assert (
                    config.get_tool_setting("samtools", "sort_mem") == "1G"
                )  # From defaults

                # Check output settings
                assert config.get_output_setting("format") == "tsv"
                assert config.get_output_setting("precision") == 2  # From defaults
        finally:
            config_path.unlink()
