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

    # Quality control thresholds
    qc: Dict[str, Any]

    # Depth analysis parameters
    depth: Dict[str, Any]

    # VCF processing parameters
    vcf: Dict[str, Any]

    # RefSeq API configuration
    refseq: Dict[str, Any]

    # Variant calling defaults
    calling: Dict[str, Any]

    # Resource management
    resources: Dict[str, Any]

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
            qc=config.get("qc", {}),
            depth=config.get("depth", {}),
            vcf=config.get("vcf", {}),
            refseq=config.get("refseq", {}),
            calling=config.get("calling", {}),
            resources=config.get("resources", {}),
            species_aliases=config.get("species_aliases", {}),
            tools=config.get("tools", {}),
            output=config.get("output", {}),
        )

    @staticmethod
    def _load_defaults() -> Dict[str, Any]:
        """Load built-in default configuration."""
        defaults_path = Path(__file__).parent / "config" / "defaults.yaml"
        tools_path = Path(__file__).parent / "config" / "tools.yaml"

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
                "qc": {
                    "min_breadth_10x": 0.80,
                    "min_callable_bases": 1000000,
                    "min_mapping_rate_ref": 0.70,
                },
                "depth": {
                    "min_contig_length": 500,
                    "default_mapq_threshold": 30,
                    "default_threads": 4,
                },
                "vcf": {
                    "dp_cap": 100,
                    "maf_bins": 51,
                    "maf_max": 0.50,
                },
                "refseq": {
                    "rate_limit_delay": 0.4,
                    "api_timeout": 30,
                    "batch_timeout": 60,
                    "user_agent": "ssiamb/0.8.0 (https://github.com/ssi-dk/ssiamb)",
                    "max_retries": 3,
                },
                "calling": {
                    "default_threads": 1,
                    "default_mapq_min": 20,
                    "default_baseq_min": 20,
                    "default_minallelefraction": 0.0,
                    "ploidy": 1,
                },
                "resources": {
                    "default_threads": 4,
                    "default_memory": "16g",
                    "timeout_short": 30,
                    "timeout_long": 300,
                    "walltime_hint": "2h",
                },
                "species_aliases": {},
                "tools": {},
                "output": {},
            }

        # Load main defaults
        config = SsiambConfig._load_yaml(defaults_path)

        # Load and merge tools configuration if it exists
        if tools_path.exists():
            tools_config = SsiambConfig._load_yaml(tools_path)
            config["tools"] = tools_config
        else:
            config["tools"] = {}

        return config

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

    def get_qc_setting(self, key: str, default: Any = None) -> Any:
        """Get a QC threshold setting."""
        return self.qc.get(key, default)

    def get_depth_setting(self, key: str, default: Any = None) -> Any:
        """Get a depth analysis setting."""
        return self.depth.get(key, default)

    def get_vcf_setting(self, key: str, default: Any = None) -> Any:
        """Get a VCF processing setting."""
        return self.vcf.get(key, default)

    def get_refseq_setting(self, key: str, default: Any = None) -> Any:
        """Get a RefSeq API setting."""
        return self.refseq.get(key, default)

    def get_calling_setting(self, key: str, default: Any = None) -> Any:
        """Get a variant calling setting."""
        return self.calling.get(key, default)

    def get_resource_setting(self, key: str, default: Any = None) -> Any:
        """Get a resource management setting."""
        return self.resources.get(key, default)


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
