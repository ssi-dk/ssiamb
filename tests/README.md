# Tests for ssiamb

This directory contains comprehensive unit and integration tests for the ssiamb package.

## Running Tests

```bash
# Install test dependencies (automatically handled by pyproject.toml)
pip install -e .[test]

# Run all tests
pytest tests/

# Run with verbose output
pytest -v tests/

# Run specific test file
pytest tests/test_bracken.py

# Run tests with coverage
pytest --cov=ssiamb tests/

# Run only unit tests (exclude integration tests)
pytest tests/ -k "not integration"

# Run only integration tests
pytest tests/ -k "integration"
```

## Test Structure

### Core Functionality Tests

- `test_bracken.py` - Tests for Bracken parsing and species selection
- `test_refdir.py` - Tests for reference directory and species resolution
- `test_io_utils.py` - Tests for I/O utilities and sample validation
- `test_config.py` - Tests for configuration management and tool paths
- `test_provenance.py` - Tests for provenance tracking and metadata
- `test_qc.py` - Tests for quality control metrics and warnings

### Analysis Module Tests

- `test_mapping.py` - Tests for read mapping functionality
- `test_calling.py` - Tests for variant calling operations
- `test_vcf_ops.py` - Tests for VCF processing and analysis
- `test_depth.py` - Tests for depth analysis with mosdepth
- `test_reuse.py` - Tests for file reuse and compatibility checking
- `test_runner.py` - Tests for pipeline execution and workflow management
- `test_summarize.py` - Tests for result summarization

### Integration Tests

- `test_mapping_integration.py` - End-to-end mapping pipeline tests
- `test_calling_integration.py` - End-to-end variant calling tests
- `test_vcf_ops_integration.py` - Full VCF processing workflow tests
- `test_depth_integration.py` - Complete depth analysis pipeline tests
- `test_reuse_integration.py` - File reuse workflow tests
- `test_runner_integration.py` - Full pipeline integration tests

### CLI and Interface Tests

- `test_cli.py` - Tests for command-line interface
- `test_dry_run.py` - Tests for dry-run functionality
- `test_error_handling.py` - Tests for error handling and user feedback

### Specialized Tests

- `test_mapping_core.py` - Core mapping algorithm tests
- `test_refseq.py` - Tests for RefSeq integration and downloads
- `test_packaging.py` - Tests for package distribution and installation
- `test_bracken_fixture.py` - Bracken test data and fixtures

## Test Categories

- **Unit Tests**: Fast, isolated tests for individual functions and classes
- **Integration Tests**: Slower tests that verify component interactions and workflows
- **CLI Tests**: Tests for command-line interface behavior and user interactions
- **Fixture Tests**: Tests using real-world data and scenarios

## Test Data

Test fixtures and sample data are stored in the `fixtures/` directory at the project root. These include reference genomes, sample reads, and expected outputs for validation.
