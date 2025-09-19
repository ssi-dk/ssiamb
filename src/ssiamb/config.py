"""
Configuration management for ssiamb.

This module handles loading and merging configuration from:
1. Built-in defaults (config/defaults.yaml)
2. User-specified config files (--config)
3. Environment variables
4. Command-line overrides
"""

from __future__ import annotations
import os
import yaml
from pathlib import Path
from typing import Dict, Any, Optional
from dataclasses import dataclass


@dataclass
class SsiambConfig:
    """
    Complete ssiamb configuration.

    This holds all configurable values that were previously hardcoded,
    allowing users to customize behavior via config files.
    """

    # Analysis thresholds
    thresholds: Dict[str, Any]

    # Species aliases for reference resolution
    species_aliases: Dict[str, str]

    # Tool-specific settings
    tools: Dict[str, Any]

    # Output formatting
    output: Dict[str, Any]

    @classmethod
    def load(cls, config_path: Optional[Path] = None) -> SsiambConfig:
        """
        Load configuration from files and environment.

        Args:
            config_path: Optional path to user config file

        Returns:
            Merged configuration object
        """
        # Start with built-in defaults
        config = cls._load_defaults()

        # Overlay user config if provided
        if config_path and config_path.exists():
            user_config = cls._load_yaml(config_path)
            config = cls._merge_configs(config, user_config)

        # Apply environment variable overrides
        config = cls._apply_env_overrides(config)

        # Normalize species alias keys for consistent lookup
        if "species_aliases" in config:
            config["species_aliases"] = cls._normalize_species_aliases(
                config["species_aliases"]
            )

        return cls(
            thresholds=config.get("thresholds", {}),
            species_aliases=config.get("species_aliases", {}),
            tools=config.get("tools", {}),
            output=config.get("output", {}),
        )

    @staticmethod
    def _load_defaults() -> Dict[str, Any]:
        """Load built-in default configuration."""
        defaults_path = Path(__file__).parent / "config" / "defaults.yaml"
        if not defaults_path.exists():
            # Fallback to minimal defaults if file doesn't exist
            return {
                "thresholds": {
                    "dp_min": 10,
                    "maf_min": 0.1,
                    "dp_cap": 100,
                    "mapq_min": 20,
                    "baseq_min": 20,
                },
                "species_aliases": {},
                "tools": {},
                "output": {},
            }
        return SsiambConfig._load_yaml(defaults_path)

    @staticmethod
    def _load_yaml(path: Path) -> Dict[str, Any]:
        """Load YAML configuration file."""
        try:
            with open(path, "r") as f:
                return yaml.safe_load(f) or {}
        except yaml.YAMLError as e:
            raise ValueError(f"Invalid YAML in config file {path}: {e}")
        except Exception as e:
            raise FileNotFoundError(f"Could not read config file {path}: {e}")

    @staticmethod
    def _merge_configs(base: Dict[str, Any], overlay: Dict[str, Any]) -> Dict[str, Any]:
        """Recursively merge two configuration dictionaries."""
        result = base.copy()

        for key, value in overlay.items():
            if (
                key in result
                and isinstance(result[key], dict)
                and isinstance(value, dict)
            ):
                result[key] = SsiambConfig._merge_configs(result[key], value)
            else:
                result[key] = value

        return result

    @staticmethod
    def _apply_env_overrides(config: Dict[str, Any]) -> Dict[str, Any]:
        """Apply environment variable overrides to configuration."""
        # Support SSIAMB_* environment variables
        env_mappings = {
            "SSIAMB_DP_MIN": ("thresholds", "dp_min", int),
            "SSIAMB_MAF_MIN": ("thresholds", "maf_min", float),
            "SSIAMB_DP_CAP": ("thresholds", "dp_cap", int),
            "SSIAMB_MAPQ_MIN": ("thresholds", "mapq_min", int),
            "SSIAMB_BASEQ_MIN": ("thresholds", "baseq_min", int),
        }

        for env_var, (section, key, type_func) in env_mappings.items():
            value = os.environ.get(env_var)
            if value is not None:
                try:
                    if section not in config:
                        config[section] = {}
                    config[section][key] = type_func(value)
                except ValueError:
                    raise ValueError(f"Invalid value for {env_var}: {value}")

        return config

    @staticmethod
    def _normalize_species_aliases(aliases: Dict[str, str]) -> Dict[str, str]:
        """
        Normalize species alias keys for consistent lookup.

        This ensures that alias keys match the normalized form used during
        species name resolution.
        """
        # Import here to avoid circular imports
        from .refdir import normalize_species_name

        normalized_aliases = {}
        for key, value in aliases.items():
            normalized_key = normalize_species_name(key)
            normalized_aliases[normalized_key] = value

        return normalized_aliases

    def get_threshold(self, key: str, default: Any = None) -> Any:
        """Get a threshold value with fallback."""
        return self.thresholds.get(key, default)

    def get_species_alias(self, species: str) -> str:
        """Get species alias, returning original name if no alias exists."""
        return self.species_aliases.get(species, species)

    def get_tool_setting(self, tool: str, key: str, default: Any = None) -> Any:
        """Get a tool-specific setting."""
        return self.tools.get(tool, {}).get(key, default)

    def get_output_setting(self, key: str, default: Any = None) -> Any:
        """Get an output formatting setting."""
        return self.output.get(key, default)


# Global configuration instance
_config: Optional[SsiambConfig] = None


def get_config() -> SsiambConfig:
    """Get the global configuration instance."""
    global _config
    if _config is None:
        _config = SsiambConfig.load()
    return _config


def set_config(config: SsiambConfig) -> None:
    """Set the global configuration instance."""
    global _config
    _config = config


def load_config(config_path: Optional[Path] = None) -> SsiambConfig:
    """Load and set configuration from file."""
    config = SsiambConfig.load(config_path)
    set_config(config)
    return config
