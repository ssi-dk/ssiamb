# Fixtures

This directory contains test fixtures for ssiamb development and testing.

## 📋 Fixture Requirements

See **[FIXTURES_REQUIREMENTS.md](FIXTURES_REQUIREMENTS.md)** for a comprehensive list of all test files needed.

## 🎯 Current Status

- [x] **Core test dataset** (Listeria monocytogenes listeria_test_sample)
- [x] **Self-mapping mode fixtures** (assembly + subsampled reads)
- [x] **Reference-mapped mode fixtures** (EGD-e reference + same reads)
- [x] **Admin reference directory examples** (with proper indexes)
- [x] **Reuse/integration test files** (BAM/VCF + error cases)
- [x] **Edge cases and error testing files** (empty, low coverage, malformed)
- [x] **Expected output validation files** (placeholder templates)
- [x] **Manifest files** (multi-sample testing)

## ✅ **Ready for Testing**

**Primary test files are now available:**
- **Species**: *Listeria monocytogenes* (99.52% purity)
- **Assembly**: 2.9 MB, N50=509kb, 15 contigs
- **Reads**: 10,000 paired-end reads (150bp)
- **Coverage**: ~10-20x estimated
- **File sizes**: Manageable for CI/testing (<1MB each)

## 📁 Directory Structure

Once populated, this directory will contain:

```
fixtures/
├── FIXTURES_REQUIREMENTS.md     # This requirements document
├── README.md                    # This file
├── self_mode/                   # Self-mapping test files
├── ref_mode/                    # Reference-mapping test files  
├── admin_refs/                  # Admin reference directory examples
├── reuse/                       # Pre-computed BAM/VCF files
├── edge_cases/                  # Error testing and edge cases
└── expected_outputs/            # Expected results for validation
```

## 🚀 Next Steps

1. Review the requirements in `FIXTURES_REQUIREMENTS.md`
2. Identify or create suitable small test datasets  
3. Populate fixture directories with test files
4. Document expected results for validation testing

For questions about fixture requirements, see the detailed specifications in `FIXTURES_REQUIREMENTS.md`.
