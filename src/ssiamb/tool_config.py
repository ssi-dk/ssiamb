"""
Tool configuration management for ssiamb.

This module handles the configuration of external bioinformatics tools,
allowing users to specify custom paths and providing helpful error messages.
"""

import os
import shutil
import subprocess
from pathlib import Path
from typing import Optional, Dict, Tuple, Set
import yaml
from rich.console import Console
from rich.table import Table

console = Console()


class ToolConfigError(Exception):
    """Exception raised when tool configuration fails."""

    pass


class ToolConfig:
    """Manages configuration for external bioinformatics tools."""

    DEFAULT_TOOLS = {
        "minimap2": {"required_for": ["self", "ref"], "conda_package": "minimap2"},
        "bwa-mem2": {"required_for": ["self", "ref"], "conda_package": "bwa-mem2"},
        "samtools": {
            "required_for": ["self", "ref", "summarize"],
            "conda_package": "samtools",
        },
        "mosdepth": {
            "required_for": ["self", "ref", "summarize"],
            "conda_package": "mosdepth",
        },
        "bcftools": {
            "required_for": ["self", "ref", "summarize"],
            "conda_package": "bcftools",
        },
        "bgzip": {
            "required_for": ["self", "ref", "summarize"],
            "conda_package": "htslib",
        },
        "pileup.sh": {"required_for": ["self", "ref"], "conda_package": "bbmap"},
        "callvariants.sh": {"required_for": ["self", "ref"], "conda_package": "bbmap"},
        "tabix": {
            "required_for": ["self", "ref", "summarize"],
            "conda_package": "htslib",
        },
    }

    def __init__(self, config_path: Optional[Path] = None):
        """Initialize tool configuration."""
        if config_path is None:
            config_path = self._get_default_config_path()

        self.config_path = config_path
        self.config = self._load_config()

    def _get_default_config_path(self) -> Path:
        """Get the default configuration file path within the package."""
        # Get the path to this module's directory
        package_dir = Path(__file__).parent
        config_dir = package_dir / "config"
        config_dir.mkdir(parents=True, exist_ok=True)
        return config_dir / "tools.yaml"

    def _load_config(self) -> Dict:
        """Load configuration from file or create default."""
        if self.config_path.exists():
            try:
                with open(self.config_path, "r") as f:
                    config = yaml.safe_load(f) or {}
            except Exception as e:
                console.print(
                    f"[yellow]Warning: Could not load config file {self.config_path}: {e}[/yellow]"
                )
                config = {}
        else:
            # Create default config structure
            config = {"tools": {tool: "" for tool in self.DEFAULT_TOOLS}}
            # Don't save immediately, let it be created when first modified

        # Ensure all tools have entries
        if "tools" not in config:
            config["tools"] = {}

        for tool in self.DEFAULT_TOOLS:
            if tool not in config["tools"]:
                config["tools"][tool] = ""

        return config

    def save_config(self) -> None:
        """Save configuration to file."""
        try:
            self.config_path.parent.mkdir(parents=True, exist_ok=True)
            with open(self.config_path, "w") as f:
                yaml.dump(self.config, f, default_flow_style=False, sort_keys=True)
        except Exception as e:
            console.print(f"[red]Error saving config to {self.config_path}: {e}[/red]")
            raise

    def find_tool_in_path(self, tool_name: str) -> Optional[str]:
        """Find tool in system PATH."""
        return shutil.which(tool_name)

    def get_tool_path(self, tool_name: str) -> Optional[str]:
        """Get the configured or auto-detected path for a tool."""
        if tool_name not in self.DEFAULT_TOOLS:
            return None

        # Get tool configuration (could be string or dict)
        tool_config = self.config["tools"].get(tool_name, "")

        # Handle legacy string format and new dict format
        if isinstance(tool_config, dict):
            # New format: look for 'path' key, or special keys for bbtools
            if tool_name == "callvariants.sh":
                user_path = tool_config.get("callvariants_path", "")
            elif tool_name == "pileup.sh":
                user_path = tool_config.get("pileup_path", "")
            else:
                user_path = tool_config.get("path", "")
        else:
            # Legacy string format
            user_path = str(tool_config) if tool_config else ""

        # If user has configured a specific path, use it
        if user_path and user_path.strip():
            return str(user_path)

        # Otherwise, try to find in PATH
        return self.find_tool_in_path(tool_name)

    def require_tool_path(self, tool_name: str) -> str:
        """Get the path for a tool or raise an error if not found."""
        path = self.get_tool_path(tool_name)
        if path is None:
            tool_info = self.DEFAULT_TOOLS.get(tool_name, {})
            conda_pkg = tool_info.get("conda_package", tool_name)
            raise ToolConfigError(
                f"Tool '{tool_name}' not found in PATH or configuration.\n"
                f"To install: conda install {conda_pkg}\n"
                f"To configure: ssiamb config set {tool_name} /path/to/{tool_name}"
            )
        return path

    def set_tool_path(self, tool_name: str, path: str) -> None:
        """Set the path for a specific tool."""
        if tool_name not in self.DEFAULT_TOOLS:
            raise ValueError(
                f"Unknown tool: {tool_name}. Available tools: {list(self.DEFAULT_TOOLS.keys())}"
            )

        # Validate that the path exists and is executable
        tool_path = Path(path)
        if not tool_path.exists():
            raise FileNotFoundError(f"Tool path does not exist: {path}")

        if not tool_path.is_file():
            raise ValueError(f"Tool path is not a file: {path}")

        if not os.access(tool_path, os.X_OK):
            raise PermissionError(f"Tool is not executable: {path}")

        # Test that the tool actually works
        try:
            subprocess.run([str(tool_path), "--help"], capture_output=True, timeout=10)
            # Most tools return 0 or 1 for --help, but shouldn't crash
        except (subprocess.TimeoutExpired, subprocess.CalledProcessError):
            pass  # Some tools may not support --help nicely
        except FileNotFoundError:
            raise ValueError(f"Tool cannot be executed: {path}")

        # Update configuration - handle both legacy and new format
        if isinstance(self.config["tools"].get(tool_name), dict):
            # New format: update the appropriate path key
            if tool_name == "callvariants.sh":
                self.config["tools"][tool_name]["callvariants_path"] = str(tool_path)
            elif tool_name == "pileup.sh":
                self.config["tools"][tool_name]["pileup_path"] = str(tool_path)
            else:
                self.config["tools"][tool_name]["path"] = str(tool_path)
        else:
            # Legacy format or create new dict format
            if tool_name in ["callvariants.sh", "pileup.sh"]:
                # Create dict format for bbtools
                if tool_name not in self.config["tools"] or not isinstance(
                    self.config["tools"][tool_name], dict
                ):
                    self.config["tools"][tool_name] = {}
                if tool_name == "callvariants.sh":
                    self.config["tools"][tool_name]["callvariants_path"] = str(
                        tool_path
                    )
                else:
                    self.config["tools"][tool_name]["pileup_path"] = str(tool_path)
            else:
                # Create dict format for other tools
                if tool_name not in self.config["tools"] or not isinstance(
                    self.config["tools"][tool_name], dict
                ):
                    self.config["tools"][tool_name] = {}
                self.config["tools"][tool_name]["path"] = str(tool_path)

        # Update related indexing and normalization configurations
        self._update_related_tool_paths(tool_name, str(tool_path))

        self.save_config()
        console.print(f"[green]✓[/green] Set {tool_name} path to: {path}")

    def _update_related_tool_paths(self, tool_name: str, tool_path: str) -> None:
        """Update paths for related indexing and normalization configurations."""
        # Map main tools to their related configurations
        related_configs = {
            "minimap2": ["minimap2_index"],
            "bwa-mem2": ["bwa-mem2_index"],
            "bcftools": ["vcf_normalization"],
        }

        if tool_name in related_configs:
            for related_config in related_configs[tool_name]:
                # Ensure the related config exists as a dict
                if related_config not in self.config["tools"]:
                    self.config["tools"][related_config] = {}
                elif not isinstance(self.config["tools"][related_config], dict):
                    self.config["tools"][related_config] = {}

                # Set the path to match the main tool
                self.config["tools"][related_config]["path"] = tool_path
                console.print(
                    f"[green]✓[/green] Also set {related_config} path to: {tool_path}"
                )

    def _reset_related_tool_paths(self, tool_name: str) -> None:
        """Reset paths for related indexing and normalization configurations."""
        # Map main tools to their related configurations
        related_configs = {
            "minimap2": ["minimap2_index"],
            "bwa-mem2": ["bwa-mem2_index"],
            "bcftools": ["vcf_normalization"],
        }

        if tool_name in related_configs:
            for related_config in related_configs[tool_name]:
                if related_config in self.config["tools"]:
                    if isinstance(self.config["tools"][related_config], dict):
                        self.config["tools"][related_config]["path"] = ""
                    else:
                        self.config["tools"][related_config] = ""
                    console.print(
                        f"[green]✓[/green] Also reset {related_config} to auto-detection"
                    )

    def reset_tool_config(self, tool_name: Optional[str] = None) -> None:
        """Reset tool configuration to auto-detection."""
        if tool_name:
            if tool_name not in self.DEFAULT_TOOLS:
                raise ValueError(f"Unknown tool: {tool_name}")

            # Reset path in tool config while preserving other parameters
            if isinstance(self.config["tools"].get(tool_name), dict):
                if tool_name == "callvariants.sh":
                    self.config["tools"][tool_name]["callvariants_path"] = ""
                elif tool_name == "pileup.sh":
                    self.config["tools"][tool_name]["pileup_path"] = ""
                else:
                    self.config["tools"][tool_name]["path"] = ""
            else:
                self.config["tools"][tool_name] = ""

            # Reset related configurations
            self._reset_related_tool_paths(tool_name)
            console.print(f"[green]✓[/green] Reset {tool_name} to auto-detection")
        else:
            # Reset all tools
            for tool in self.DEFAULT_TOOLS:
                if isinstance(self.config["tools"].get(tool), dict):
                    if tool == "callvariants.sh":
                        self.config["tools"][tool]["callvariants_path"] = ""
                    elif tool == "pileup.sh":
                        self.config["tools"][tool]["pileup_path"] = ""
                    else:
                        self.config["tools"][tool]["path"] = ""
                else:
                    self.config["tools"][tool] = ""
                # Reset related configurations for each tool
                self._reset_related_tool_paths(tool)
            console.print("[green]✓[/green] Reset all tools to auto-detection")

        self.save_config()

    def check_tool_availability(
        self, tool_name: str
    ) -> Tuple[bool, Optional[str], str]:
        """Check if a tool is available and return status info."""
        path = self.get_tool_path(tool_name)

        if path and Path(path).exists():
            try:
                # Try to get version info
                result = subprocess.run(
                    [path, "--version"], capture_output=True, text=True, timeout=5
                )
                if result.returncode == 0 and result.stdout:
                    version_info = result.stdout.strip().split("\n")[0]
                else:
                    version_info = "Available"
                return True, path, version_info
            except Exception:
                return True, path, "Available (version unknown)"

        return False, None, "Not found"

    def check_all_tools(
        self, required_for: Optional[str] = None
    ) -> Dict[str, Tuple[bool, Optional[str], str]]:
        """Check availability of all tools, optionally filtered by requirement."""
        results = {}

        for tool_name, tool_info in self.DEFAULT_TOOLS.items():
            if required_for is None or required_for in tool_info["required_for"]:
                results[tool_name] = self.check_tool_availability(tool_name)

        return results

    def get_missing_tools_message(self, mode: str) -> str:
        """Get a helpful message about missing tools for a specific mode."""
        missing_tools = []
        tool_status = self.check_all_tools(required_for=mode)

        for tool_name, (available, path, status) in tool_status.items():
            if not available:
                missing_tools.append(tool_name)

        if not missing_tools:
            return ""

        message_parts = [
            f"Missing required tools for '{mode}' mode: {', '.join(missing_tools)}",
            "",
            "To fix this issue, you can:",
            "",
            "1. Install via conda (recommended):",
        ]

        # Group by conda package
        conda_packages: Set[str] = set()
        for tool in missing_tools:
            conda_packages.add(str(self.DEFAULT_TOOLS[tool]["conda_package"]))

        message_parts.append(
            f"   conda install -c bioconda {' '.join(sorted(conda_packages))}"
        )
        message_parts.extend(
            [
                "",
                "2. Or set custom tool paths:",
            ]
        )

        for tool in missing_tools:
            message_parts.append(f"   ssiamb config set {tool} /path/to/{tool}")

        message_parts.extend(
            [
                "",
                "3. Check current tool configuration:",
                "   ssiamb config check",
                "",
                "For more help: ssiamb config --help",
            ]
        )

        return "\n".join(message_parts)

    def print_config_table(self) -> None:
        """Print a formatted table of current tool configuration."""
        table = Table(title="SSI Ambiguous Site Detection Tool Configuration")
        table.add_column("Tool", style="cyan", no_wrap=True)
        table.add_column("Status", style="green")
        table.add_column("Path", style="yellow")
        table.add_column("Required For", style="blue")
        table.add_column("Conda Package", style="magenta")

        for tool_name, tool_info in self.DEFAULT_TOOLS.items():
            available, path, status = self.check_tool_availability(tool_name)

            status_display = "✓ Available" if available else "✗ Missing"
            path_display = path or "Not found"
            required_for = ", ".join(tool_info["required_for"])
            conda_package = str(tool_info["conda_package"])

            table.add_row(
                tool_name, status_display, path_display, required_for, conda_package
            )

        console.print(table)
        console.print(f"\nConfiguration file: {self.config_path}")


# Global instance
config = ToolConfig()


def get_tool_path(tool_name: str) -> str:
    """Get the configured or auto-detected path for a tool, raising error if not found."""
    return config.require_tool_path(tool_name)


def get_tool_path_optional(tool_name: str) -> Optional[str]:
    """Get the configured or auto-detected path for a tool, returning None if not found."""
    return config.get_tool_path(tool_name)


def check_required_tools_for_mode(mode: str) -> None:
    """Check that all tools required for a given mode are available."""
    missing_tools = []
    for tool_name, tool_info in config.DEFAULT_TOOLS.items():
        if mode in tool_info.get("required_for", []):
            if not config.get_tool_path(tool_name):
                missing_tools.append(tool_name)

    if missing_tools:
        error_message = config.get_missing_tools_message(mode)
        from .errors import ExternalToolError

        raise ExternalToolError(error_message)
