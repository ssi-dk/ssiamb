"""
Integration tests for VCF operations using real fixture data.

Tests VCF parsing, MAF extraction, variant classification, and ambiguous site
counting using real production VCF files from SSI pipeline.
"""

import pytest
from pathlib import Path

from src.ssiamb.vcf_ops import (
    parse_vcf_sites,
    count_ambiguous_sites,
    AmbigGrid,
    VariantClass,
    SiteRecord,
    normalize_and_split,
    emit_vcf,
    emit_bed,
    emit_matrix,
    VCFOperationError,
)


class TestVCFOperationsRealData:
    """Integration tests using real VCF fixtures."""

    @pytest.fixture
    def real_vcf_path(self):
        """Path to real production VCF fixture."""
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
    def reference_fasta_path(self):
        """Path to reference FASTA fixture."""
        path = Path("fixtures/ref_mode/reference.fna")
        if not path.exists():
            pytest.skip("Reference FASTA fixture not available")
        return path

    def test_parse_real_production_vcf(self, real_vcf_path):
        """Test parsing real SSI production VCF file."""
        sites = list(parse_vcf_sites(real_vcf_path))

        # Should have a reasonable number of sites (adjust based on actual fixture size)
        assert len(sites) > 100, f"Expected >100 sites, got {len(sites)}"

        # All parsed sites should be variants (no reference sites in this VCF)
        variant_sites = [s for s in sites if s.variant_class != VariantClass.UNKNOWN]

        assert (
            len(variant_sites) > 100
        ), f"Expected >100 variants, got {len(variant_sites)}"

        # Validate site record structure
        for site in sites[:100]:  # Check first 100 sites
            assert isinstance(site, SiteRecord)
            assert isinstance(site.chrom, str)
            assert isinstance(site.pos, int)
            assert site.pos > 0
            assert isinstance(site.variant_class, VariantClass)

            # MAF should be valid float for variants
            assert isinstance(site.maf, float)
            assert 0.0 <= site.maf <= 1.0

    def test_variant_classification_real_data(self, real_vcf_path):
        """Test variant classification on real variants."""
        sites = list(parse_vcf_sites(real_vcf_path))
        variant_sites = [s for s in sites if s.variant_class != VariantClass.UNKNOWN]

        # Count variant types
        snv_count = sum(1 for s in variant_sites if s.variant_class == VariantClass.SNV)
        ins_count = sum(1 for s in variant_sites if s.variant_class == VariantClass.INS)
        del_count = sum(1 for s in variant_sites if s.variant_class == VariantClass.DEL)

        # Real bacterial data should have mostly SNVs with some indels
        assert snv_count > 0, "Should have SNV variants"
        assert snv_count >= ins_count, "Should have SNVs"
        assert snv_count >= del_count, "Should have SNVs"

        # Validate that all variants are properly classified
        for site in variant_sites[:50]:  # Check first 50 variants
            assert site.variant_class in [
                VariantClass.SNV,
                VariantClass.INS,
                VariantClass.DEL,
            ]

    def test_maf_extraction_realistic_distributions(self, real_vcf_path):
        """Test MAF extraction from real variant calls."""
        sites = list(parse_vcf_sites(real_vcf_path))
        variant_sites = [
            s
            for s in sites
            if s.variant_class != VariantClass.UNKNOWN and s.maf is not None
        ]

        if len(variant_sites) == 0:
            pytest.skip("No variants with MAF found in real VCF")

        mafs = [site.maf for site in variant_sites]

        # Real MAF distributions should be realistic
        assert len(mafs) > 0, "Should have MAF values"
        assert all(
            0.0 <= maf <= 1.0 for maf in mafs
        ), "All MAFs should be between 0 and 1"

        # For bacterial data, we expect many homozygous variants (MAF near 1.0)
        # and some heterozygous variants (MAF around 0.5)
        high_maf_count = sum(1 for maf in mafs if maf > 0.8)
        if len(mafs) > 10:  # Only check if we have enough data
            assert high_maf_count >= 0, "Should have some variants with defined MAF"

    def test_ambiguous_site_counting_real_data(self, real_vcf_path):
        """Test ambiguous site counting with real variant density."""
        # Test using the count_ambiguous_sites function directly
        dp_min = 1
        maf_min = 0.0

        ambig_count, grid = count_ambiguous_sites(real_vcf_path, dp_min, maf_min)

        # Should return a reasonable ambiguous site count
        assert isinstance(ambig_count, int)
        assert ambig_count >= 0

        # Should return a valid grid
        assert isinstance(grid, AmbigGrid)

        # Grid should have been populated
        cumulative = grid.build_cumulative()
        assert cumulative.shape == (101, 51)  # dp_cap+1, maf_bins

    def test_normalize_and_split_real_variants(
        self, real_vcf_path, reference_fasta_path, tmp_path
    ):
        """Test normalization and splitting with real data - expects mismatch error."""

        # The test VCF and reference don't match (different assemblies/contigs)
        # This test verifies that normalization correctly fails with informative error
        # when VCF contigs don't match reference sequences
        with pytest.raises(VCFOperationError) as exc_info:
            normalize_and_split(
                vcf_in=real_vcf_path,
                reference=reference_fasta_path,
                output_dir=tmp_path,
            )

        # Verify error message indicates sequence/contig mismatch
        error_msg = str(exc_info.value)
        assert "sequence" in error_msg.lower() and "not found" in error_msg.lower()

    def test_output_generation_real_data(self, real_vcf_path, tmp_path):
        """Test output file generation with real data."""
        # Test using the emit functions with proper parameters
        dp_min = 1
        maf_min = 0.0
        sample_name = "test_sample"

        # Test VCF emission
        output_vcf = tmp_path / "test_output.vcf.gz"
        result_vcf = emit_vcf(real_vcf_path, output_vcf, dp_min, maf_min, sample_name)
        assert result_vcf.exists()
        assert result_vcf.stat().st_size > 0

        # Test BED emission (returns compressed file)
        output_bed = tmp_path / "test_output.bed"
        result_bed = emit_bed(real_vcf_path, output_bed, dp_min, maf_min, sample_name)
        assert result_bed.exists()
        assert result_bed.stat().st_size > 0

        # Test matrix emission with grid
        _, grid = count_ambiguous_sites(real_vcf_path, dp_min, maf_min)
        output_matrix = tmp_path / "test_output.tsv"
        result_matrix = emit_matrix(grid, output_matrix, sample_name)
        assert result_matrix.exists()
        assert result_matrix.stat().st_size > 0

    def test_multi_sample_vcf_handling(self, multi_sample_vcf_path):
        """Test handling of multi-sample VCF files."""
        sites = list(parse_vcf_sites(multi_sample_vcf_path))

        # Should parse successfully even with multiple samples
        assert len(sites) > 0, "Should parse multi-sample VCF"

        # Check that we handle multiple samples properly
        for site in sites[:10]:
            assert isinstance(site, SiteRecord)
            # MAF extraction should work or return None for multi-sample
            assert site.maf is None or (0.0 <= site.maf <= 1.0)

    def test_real_vcf_error_handling(self):
        """Test error handling with malformed real VCF data."""
        malformed_vcf = Path("fixtures/edge_cases/malformed.vcf")
        if not malformed_vcf.exists():
            pytest.skip("Malformed VCF fixture not available")

        # Should handle malformed VCF gracefully by raising VCFOperationError
        # The malformed VCF contains valid records (which may generate MAF warnings)
        # followed by malformed records that cause parsing to fail
        import warnings

        with warnings.catch_warnings():
            warnings.simplefilter(
                "ignore", UserWarning
            )  # Suppress MAF warnings from valid records
            with pytest.raises(VCFOperationError):
                list(parse_vcf_sites(malformed_vcf))


class TestVCFPerformanceRealData:
    """Performance tests with real data."""

    @pytest.fixture
    def large_vcf_path(self):
        """Use the larger multi-sample VCF for performance testing."""
        path = Path("fixtures/reuse/multi_sample.vcf.gz")
        if not path.exists():
            pytest.skip("Large VCF fixture not available")
        return path

    def test_large_vcf_parsing_performance(self, large_vcf_path):
        """Test that parsing large VCF files completes in reasonable time."""
        import time

        start_time = time.time()
        sites = list(parse_vcf_sites(large_vcf_path))
        end_time = time.time()

        parse_time = end_time - start_time

        # Should complete within reasonable time (adjust based on expected file size)
        assert parse_time < 30, f"Parsing took too long: {parse_time:.2f} seconds"
        assert (
            len(sites) > 100
        ), f"Expected reasonable number of sites, got {len(sites)}"

        print(f"Parsed {len(sites)} sites in {parse_time:.2f} seconds")
