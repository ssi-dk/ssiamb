# Bioconda Submission Guide for ssiamb

## Steps to Submit ssiamb to Bioconda

### 1. Fork the bioconda-recipes repository
- Go to https://github.com/bioconda/bioconda-recipes
- Click "Fork" to create your own fork

### 2. Clone your fork locally
```bash
git clone https://github.com/YOUR_USERNAME/bioconda-recipes.git
cd bioconda-recipes
```

### 3. Create a new branch for ssiamb
```bash
git checkout -b add-ssiamb
```

### 4. Copy the recipe
Copy the contents of `recipes/ssiamb/meta.yaml` to:
```
bioconda-recipes/recipes/ssiamb/meta.yaml
```

### 5. Test the recipe locally (optional but recommended)
```bash
cd bioconda-recipes
conda-build recipes/ssiamb
```

### 6. Commit and push the changes
```bash
git add recipes/ssiamb/
git commit -m "Add ssiamb recipe"
git push origin add-ssiamb
```

### 7. Create a Pull Request
- Go to your fork on GitHub
- Click "Compare & pull request"
- Title: "Add ssiamb recipe"
- Description should include:
  - Brief description of the tool
  - Link to the original repository
  - Any relevant notes about dependencies

### 8. Wait for Review
- Bioconda maintainers will review the PR
- They may request changes or ask questions
- The automated tests will run to ensure the recipe works

## Current Recipe Status
- ✅ Recipe validates locally with `conda-build --check`
- ✅ Version 0.1.0 with correct SHA256 hash
- ✅ All dependencies specified
- ✅ Maintainer: PovilasMat

## Notes
- The recipe is ready for submission
- Once merged, it will be available via `conda install -c bioconda ssiamb`
- BioContainer will be automatically built from the Bioconda package