"""
Entry point for python -m ssiamb.

This allows the package to be executed as a module:
    python -m ssiamb --help
"""

from .cli import app

if __name__ == "__main__":
    app()
