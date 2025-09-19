# PyPI Trusted Publishing Setup

This document outlines the steps to set up PyPI Trusted Publishing for the ssiamb project.

## Prerequisites

1. **PyPI Account**: You need a PyPI account with maintainer privileges for the ssiamb project
2. **GitHub Repository**: The repository must be public or you need GitHub Pro/Enterprise

## Setting Up PyPI Trusted Publishing

### 1. Create PyPI Project (First Time Only)

If the project doesn't exist on PyPI yet:

1. Build and upload manually first time:

   ```bash
   python -m pip install build twine
   python -m build
   python -m twine upload dist/*
   ```

### 2. Configure Trusted Publishing on PyPI

1. Go to <https://pypi.org/manage/project/ssiamb/settings/publishing/>
2. Click "Add a new pending publisher"
3. Fill in the form:
   - **PyPI Project Name**: `ssiamb`
   - **Owner**: `ssi-dk`
   - **Repository name**: `ssiamb`
   - **Workflow filename**: `publish.yml`
   - **Environment name**: `pypi` (optional but recommended)

### 3. Configure GitHub Environment (Recommended)

1. Go to your GitHub repository settings
2. Navigate to Environments
3. Create a new environment named `pypi`
4. Add protection rules:
   - Required reviewers (optional)
   - Deployment branches: Only protected branches
5. **No secrets needed!** Trusted publishing uses OIDC authentication automatically

### 4. Create Release

Once everything is set up, create releases by:

1. Update version in `pyproject.toml`
2. Commit changes
3. Create and push a tag:

   ```bash
   git tag v1.0.0
   git push origin v1.0.0
   ```

The GitHub Action will automatically:

- Build the package
- Publish to PyPI using OIDC
- Create a GitHub release
- Sign artifacts with Sigstore

## Troubleshooting

### Common Issues

1. **"Trusted publishing is not configured"**:
   - Ensure the pending publisher is correctly configured on PyPI
   - Check that repository owner, name, and workflow filename match exactly

2. **"Environment protection rules"**:
   - If using environment protection, ensure the tag push is from a protected branch
   - Or remove environment protection for automated releases

3. **"Package already exists"**:
   - Ensure version number is incremented in `pyproject.toml`
   - PyPI doesn't allow overwriting existing versions

### Verifying Setup

Test the setup with a pre-release version:

```bash
# Update version to something like 0.1.0rc1
git tag v0.1.0rc1
git push origin v0.1.0rc1
```

Check the GitHub Actions tab for any errors in the workflow.

## Security Benefits

Trusted Publishing provides:

- No need to store PyPI tokens in GitHub secrets
- Short-lived, per-release tokens
- Cryptographic verification of the publishing source
- Sigstore signing for supply chain security
