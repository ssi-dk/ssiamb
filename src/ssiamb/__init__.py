"""Ambiguous-site counting for reads mapped to their own assembly."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("ssiamb")
except PackageNotFoundError:  # pragma: no cover - source tree without install metadata
    __version__ = "0+unknown"
