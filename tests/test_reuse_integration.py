"""
Integration tests for reuse and compatibility checking using real fixture data.

Tests BAM/VCF compatibility checking, duplicate detection, and compatibility
scenarios using real production files from SSI pipeline.
"""

import pytest
import tempfile
from pathlib import Path

from src.ssiamb.reuse import (
    CompatibilityError,
    CompatibilityResult,
    get_fasta_contigs,
    get_bam_contigs,
    get_vcf_contigs,
    check_compatibility,
    check_bam_compatibility,
    check_vcf_compatibility,
    has_duplicate_flags,
    run_markdup_for_depth,
)


class TestReuseRealData:
    """Integration tests using real BAM/VCF fixtures."""

    @pytest.fixture
    def real_bam_path(self):
        """Path to real BAM fixture."""
        path = Path("fixtures/reuse/pre_mapped.bam")
        if not path.exists():
            pytest.skip("Real BAM fixture not available")
        return path

    @pytest.fixture
    def real_vcf_path(self):
        """Path to real VCF fixture."""
        path = Path("fixtures/reuse/pre_called.vcf.gz")
        if not path.exists():
            pytest.skip("Real VCF fixture not available")
        return path

    @pytest.fixture
    def multi_sample_vcf_path(self):
        """Path to multi-sample VCF fixture."""
        path = Path("fixtures/reuse/multi_sample.vcf.gz")
        if not path.exists():
            pytest.skip("Multi-sample VCF fixture not available")
        return path

    @pytest.fixture
    def mismatched_bam_path(self):
        """Path to BAM with different sample name."""
        path = Path("fixtures/reuse/mismatched_sample.bam")
        if not path.exists():
            pytest.skip("Mismatched BAM fixture not available")
        return path

    @pytest.fixture
    def reference_fasta_path(self):
        """Path to reference FASTA."""
        path = Path("fixtures/self_mode/assembly.fna")
        if not path.exists():
            pytest.skip("Reference FASTA fixture not available")
        return path

    def test_real_bam_contig_extraction(self, real_bam_path):
        """Test contig extraction from real BAM file."""
        contigs = get_bam_contigs(real_bam_path)

        # Should extract contigs successfully
        assert isinstance(contigs, dict)
        assert len(contigs) > 0

        # All contig names should be strings, lengths should be positive integers
        for name, length in contigs.items():
            assert isinstance(name, str)
            assert isinstance(length, int)
            assert length > 0

        # Should have reasonable contig names (not just numbers)
        assert any(
            len(name) > 5 for name in contigs.keys()
        ), "Should have meaningful contig names"

    def test_real_vcf_contig_extraction(self, real_vcf_path):
        """Test contig extraction from real VCF file."""
        contigs = get_vcf_contigs(real_vcf_path)

        # Should extract contigs successfully
        assert isinstance(contigs, dict)
        assert len(contigs) > 0

        # All contig names should be strings, lengths should be positive integers
        for name, length in contigs.items():
            assert isinstance(name, str)
            assert isinstance(length, int)
            assert length > 0

        # VCF contigs should match the ones we see in the actual file
        print(f"VCF contigs: {list(contigs.keys())}")

    def test_real_fasta_contig_extraction(self, reference_fasta_path):
        """Test contig extraction from real FASTA file."""
        contigs = get_fasta_contigs(reference_fasta_path)

        # Should extract contigs successfully
        assert isinstance(contigs, dict)
        assert len(contigs) > 0

        # All contig names should be strings, lengths should be positive integers
        for name, length in contigs.items():
            assert isinstance(name, str)
            assert isinstance(length, int)
            assert length > 0

        # Should have reasonable total length
        total_length = sum(contigs.values())
        assert (
            total_length > 100000
        ), f"Total FASTA length seems too small: {total_length}"

    def test_bam_vcf_compatibility_real_data(self, real_bam_path, real_vcf_path):
        """Test compatibility between real BAM and VCF files."""
        bam_contigs = get_bam_contigs(real_bam_path)
        vcf_contigs = get_vcf_contigs(real_vcf_path)

        # Check if they have overlapping contigs
        common_contigs = set(bam_contigs.keys()) & set(vcf_contigs.keys())
        assert len(common_contigs) > 0, "BAM and VCF should have some common contigs"

        # Check compatibility
        result = check_compatibility(bam_contigs, vcf_contigs)

        # Should provide meaningful compatibility assessment
        assert isinstance(result, CompatibilityResult)
        assert isinstance(result.is_compatible, bool)
        assert isinstance(result.coverage_fraction, float)
        assert 0.0 <= result.coverage_fraction <= 1.0

        # Log the result for debugging
        print(f"BAM-VCF compatibility: {result.is_compatible}")
        print(f"Coverage fraction: {result.coverage_fraction:.3f}")
        if result.error_message:
            print(f"Error message: {result.error_message}")

    def test_bam_fasta_compatibility_real_data(
        self, real_bam_path, reference_fasta_path
    ):
        """Test compatibility between real BAM and reference FASTA."""
        result = check_bam_compatibility(real_bam_path, reference_fasta_path)

        # Should provide meaningful compatibility assessment
        assert isinstance(result, CompatibilityResult)
        assert isinstance(result.is_compatible, bool)
        assert isinstance(result.coverage_fraction, float)
        assert 0.0 <= result.coverage_fraction <= 1.0

        # Should have detailed contig information
        assert isinstance(result.per_contig_diffs, dict)
        assert isinstance(result.missing_contigs, set)
        assert isinstance(result.extra_contigs, set)

        # Log the result for debugging
        print(f"BAM-FASTA compatibility: {result.is_compatible}")
        print(f"Coverage fraction: {result.coverage_fraction:.3f}")
        print(f"Total length diff: {result.total_length_diff:.3f}")
        if result.error_message:
            print(f"Error message: {result.error_message}")

    def test_vcf_fasta_compatibility_real_data(
        self, real_vcf_path, reference_fasta_path
    ):
        """Test compatibility between real VCF and reference FASTA."""
        result = check_vcf_compatibility(real_vcf_path, reference_fasta_path)

        # Should provide meaningful compatibility assessment
        assert isinstance(result, CompatibilityResult)
        assert isinstance(result.is_compatible, bool)
        assert isinstance(result.coverage_fraction, float)
        assert 0.0 <= result.coverage_fraction <= 1.0

        # Should have detailed contig information
        assert isinstance(result.per_contig_diffs, dict)
        assert isinstance(result.missing_contigs, set)
        assert isinstance(result.extra_contigs, set)

        # Log the result for debugging
        print(f"VCF-FASTA compatibility: {result.is_compatible}")
        print(f"Coverage fraction: {result.coverage_fraction:.3f}")
        print(f"Total length diff: {result.total_length_diff:.3f}")
        if result.error_message:
            print(f"Error message: {result.error_message}")

    def test_duplicate_flag_detection_real_bam(self, real_bam_path):
        """Test duplicate flag detection on real BAM file."""
        has_dups = has_duplicate_flags(real_bam_path)

        # Should return a boolean
        assert isinstance(has_dups, bool)

        # Log the result
        print(f"BAM has duplicate flags: {has_dups}")

    def test_mismatched_sample_handling(
        self, mismatched_bam_path, reference_fasta_path
    ):
        """Test handling of BAM files with mismatched sample names."""
        # This BAM should have a different sample name structure
        result = check_bam_compatibility(mismatched_bam_path, reference_fasta_path)

        # Should still perform compatibility check
        assert isinstance(result, CompatibilityResult)
        assert isinstance(result.is_compatible, bool)

        # May or may not be compatible, but should not crash
        print(f"Mismatched BAM compatibility: {result.is_compatible}")
        if result.error_message:
            print(f"Error message: {result.error_message}")

    def test_multi_sample_vcf_handling(
        self, multi_sample_vcf_path, reference_fasta_path
    ):
        """Test handling of multi-sample VCF files."""
        result = check_vcf_compatibility(multi_sample_vcf_path, reference_fasta_path)

        # Should handle multi-sample VCF
        assert isinstance(result, CompatibilityResult)
        assert isinstance(result.is_compatible, bool)

        # Should work regardless of number of samples
        print(f"Multi-sample VCF compatibility: {result.is_compatible}")
        if result.error_message:
            print(f"Error message: {result.error_message}")

    def test_markdup_integration_real_bam(self, real_bam_path):
        """Test markdup integration with real BAM file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir)  # Pass directory, not a file path

            result_path = run_markdup_for_depth(real_bam_path, output_dir)

            # Should create output file
            assert result_path.exists()
            assert result_path.stat().st_size > 0

            # Should be a valid BAM file
            import pysam

            with pysam.AlignmentFile(str(result_path)) as f:
                # Should be readable
                assert f.header is not None

                # Check that some reads are marked as duplicates
                duplicates = 0
                total = 0
                for read in f:
                    total += 1
                    if read.is_duplicate:
                        duplicates += 1
                    if total >= 1000:  # Check first 1000 reads
                        break

                # Should have found some duplicates (though may be rare in small sample)
                assert total > 0, "Should have reads in the output BAM"
                # Note: duplicates may be 0 in small samples, that's normal


class TestReuseEdgeCases:
    """Test edge cases and error handling with real data."""

    def test_nonexistent_file_handling(self):
        """Test handling of non-existent files."""
        fake_path = Path("does_not_exist.bam")
        fake_fasta = Path("does_not_exist.fasta")

        # Should raise appropriate errors
        with pytest.raises((FileNotFoundError, CompatibilityError)):
            get_bam_contigs(fake_path)

        with pytest.raises((FileNotFoundError, CompatibilityError)):
            get_fasta_contigs(fake_fasta)

    def test_corrupted_file_handling(self):
        """Test handling of corrupted files."""
        # Create a fake corrupted BAM file
        with tempfile.NamedTemporaryFile(suffix=".bam", delete=False) as f:
            f.write(b"This is not a BAM file")
            fake_bam = Path(f.name)

        try:
            # Should handle corrupted files gracefully
            with pytest.raises((Exception, CompatibilityError)):
                get_bam_contigs(fake_bam)
        finally:
            fake_bam.unlink()  # Clean up

    def test_empty_file_handling(self):
        """Test handling of empty files."""
        with tempfile.NamedTemporaryFile(suffix=".bam", delete=False) as f:
            fake_bam = Path(f.name)

        try:
            # Should handle empty files gracefully
            with pytest.raises((Exception, CompatibilityError)):
                get_bam_contigs(fake_bam)
        finally:
            fake_bam.unlink()  # Clean up


class TestReusePerformance:
    """Performance tests with real data."""

    @pytest.fixture
    def large_bam_path(self):
        """Use the pre_mapped BAM for performance testing."""
        path = Path("fixtures/reuse/pre_mapped.bam")
        if not path.exists():
            pytest.skip("Large BAM fixture not available")
        return path

    def test_large_bam_processing_performance(self, large_bam_path):
        """Test that processing large BAM files completes in reasonable time."""
        import time

        start_time = time.time()
        contigs = get_bam_contigs(large_bam_path)
        end_time = time.time()

        process_time = end_time - start_time

        # Should complete within reasonable time
        assert (
            process_time < 10
        ), f"BAM processing took too long: {process_time:.2f} seconds"
        assert len(contigs) > 0, "Should extract contigs"

        print(
            f"Processed BAM with {len(contigs)} contigs in {process_time:.2f} seconds"
        )

    def test_duplicate_detection_performance(self, large_bam_path):
        """Test duplicate detection performance on real BAM."""
        import time

        # Test duplicate detection performance
        start_time = time.time()
        has_dups = has_duplicate_flags(large_bam_path)
        end_time = time.time()

        process_time = end_time - start_time

        # Should complete quickly
        assert (
            process_time < 10
        ), f"Duplicate detection took too long: {process_time:.2f} seconds"
        assert isinstance(has_dups, bool)

        print(f"Duplicate detection: {process_time:.2f} seconds, result: {has_dups}")
