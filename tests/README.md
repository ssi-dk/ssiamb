# Tests for ssiamb

This directory contains unit tests for the ssiamb package.

## Running Tests

```bash
# Install test dependencies
pip install pytest

# Run all tests
pytest tests/

# Run with verbose output
pytest -v tests/

# Run specific test file
pytest tests/test_bracken.py
```

## Test Structure

- `test_bracken.py` - Tests for Bracken parsing and species selection
- `test_refdir.py` - Tests for reference directory and species resolution
- `test_models.py` - Tests for data models and configuration
- `test_io_utils.py` - Tests for I/O utilities and sample validation