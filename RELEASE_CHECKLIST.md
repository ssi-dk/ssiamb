# Release Checklist

This document outlines the steps that must be completed before creating a tagged release that will be automatically uploaded to PyPI.

## Pre-Release Checklist

### 1. Version Update

- [ ] Update version in `src/ssiamb/version.py`
- [ ] Update version in `pyproject.toml`
- [ ] Ensure both version numbers match exactly

### 2. Code Quality Checks

Run all code quality checks to ensure the codebase is ready for release:

```bash
# Install/update development dependencies
pip install -e ".[dev]"

# Format code with black
black src/ tests/

# Lint code with ruff
ruff check src/ tests/

# Type checking with mypy
mypy src/

# Run full test suite
pytest tests/ -v --cov=src/ssiamb --cov-report=term-missing
```

All checks must pass before proceeding:

- [ ] Black formatting: ✅ No issues
- [ ] Ruff linting: ✅ All checks passed  
- [ ] MyPy type checking: ✅ No issues found
- [ ] Pytest test suite: ✅ All tests pass
- [ ] Test coverage: ✅ Adequate coverage (aim for >90%)

### 3. Documentation Updates

- [ ] Update `README.md` if needed (new features, installation, usage)
- [ ] Update `CHANGELOG.md` or release notes with changes in this version
- [ ] Verify all example commands in documentation still work

### 4. Dependencies and Environment

- [ ] Ensure `requirements.txt` is up to date (if used)
- [ ] Verify `environment.yml` includes all necessary packages for development
- [ ] Check that `pyproject.toml` dependencies are current and minimal

### 5. Final Validation

Run a comprehensive validation:

```bash
# Clean build environment
rm -rf dist/ build/ src/*.egg-info/

# Create fresh environment and test installation
python -m venv test_env
source test_env/bin/activate  # or test_env\Scripts\activate on Windows
pip install .

# Test CLI functionality
ssiamb --version
ssiamb --help

# Clean up
deactivate
rm -rf test_env
```

- [ ] Package builds successfully
- [ ] CLI installs and runs correctly
- [ ] Version output matches expected version

## Creating the Release

Once all checklist items are completed:

### 1. Commit and Tag

```bash
# Stage and commit all changes
git add .
git commit -m "Release version X.Y.Z"

# Create and push tag
git tag vX.Y.Z
git push origin main
git push origin vX.Y.Z
```

### 2. GitHub Release

1. Go to the GitHub repository releases page
2. Click "Create a new release"
3. Select the tag you just created (vX.Y.Z)
4. Fill in release title and description
5. Publish the release

### 3. Automated PyPI Upload

The GitHub Actions workflow will automatically:

- Build the package
- Run quality checks
- Upload to PyPI (via trusted publishing)

Monitor the GitHub Actions workflow to ensure successful deployment.

## Post-Release

- [ ] Verify package is available on PyPI
- [ ] Test installation from PyPI: `pip install ssiamb==X.Y.Z`
- [ ] Update development version (e.g., bump to X.Y.Z+1-dev) if using development versions

## Emergency Release Process

If critical fixes are needed:

1. Create hotfix branch from main
2. Make minimal necessary changes
3. Follow abbreviated checklist (focus on critical checks)
4. Create patch release (increment patch version)
5. Merge hotfix back to main

## Notes

- All releases use semantic versioning (MAJOR.MINOR.PATCH)
- GitHub Actions will handle the actual PyPI upload via trusted publishing
- Never upload manually to PyPI once trusted publishing is set up
- Always test in a clean environment before releasing
