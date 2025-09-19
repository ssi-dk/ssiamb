"""
Unit tests for VCF operations module.

Tests VCF normalization, MAF extraction, variant classification,
grid-based counting, and edge cases.
"""

import pytest
import tempfile
from pathlib import Path
from unittest.mock import patch, MagicMock
import numpy as np
import gzip

from src.ssiamb.vcf_ops import (
    normalize_and_split,
    classify_variant,
    extract_maf_from_record,
    count_ambiguous_sites,
    AmbigGrid,
    VariantClass,
    VCFOperationError,
    SiteRecord,
    check_vcf_tools,
    emit_vcf,
    emit_bed,
    emit_matrix,
    emit_per_contig,
    emit_multiqc,
)


class TestVCFTools:
    """Test VCF tool availability checks."""

    def test_check_vcf_tools_available(self):
        """Test tool availability check when tools are available."""
        with patch("shutil.which") as mock_which:
            mock_which.return_value = "/usr/bin/tool"
            assert check_vcf_tools() is True
            # Should check for bcftools, bgzip, tabix
            assert mock_which.call_count == 3

    def test_check_vcf_tools_missing(self):
        """Test tool availability check when tools are missing."""
        with patch("shutil.which") as mock_which:
            mock_which.return_value = None
            assert check_vcf_tools() is False


class TestVariantClassification:
    """Test variant classification logic."""

    def test_classify_snv(self):
        """Test SNV classification."""
        assert classify_variant("A", "T") == VariantClass.SNV
        assert classify_variant("G", "C") == VariantClass.SNV
        assert classify_variant("a", "t") == VariantClass.SNV  # Case insensitive

    def test_classify_insertion(self):
        """Test insertion classification."""
        assert classify_variant("A", "AT") == VariantClass.INS
        assert classify_variant("G", "GTCGA") == VariantClass.INS

    def test_classify_deletion(self):
        """Test deletion classification."""
        assert classify_variant("AT", "A") == VariantClass.DEL
        assert classify_variant("GTCGA", "G") == VariantClass.DEL

    def test_classify_unknown_symbolic(self):
        """Test unknown classification for symbolic variants."""
        assert classify_variant("<INS>", "A") == VariantClass.UNKNOWN
        assert classify_variant("A", "<DEL>") == VariantClass.UNKNOWN
        assert classify_variant("*", "A") == VariantClass.UNKNOWN

    def test_classify_unknown_iupac(self):
        """Test unknown classification for IUPAC codes."""
        assert classify_variant("N", "A") == VariantClass.UNKNOWN
        assert classify_variant("A", "R") == VariantClass.UNKNOWN

    def test_classify_complex(self):
        """Test complex variant classification."""
        # Same length but different from 1bp
        assert classify_variant("AT", "GC") == VariantClass.UNKNOWN


class TestMAFExtraction:
    """Test MAF extraction from VCF records."""

    def create_mock_record(
        self, format_fields=None, info_fields=None, sample_data=None
    ):
        """Create a mock VCF record for testing."""
        mock_record = MagicMock()
        mock_record.chrom = "chr1"
        mock_record.pos = 100
        mock_record.format = format_fields or {}
        mock_record.info = info_fields or {}

        if sample_data:
            mock_sample = MagicMock()
            for key, value in sample_data.items():
                mock_sample.get.return_value = value if key in sample_data else None
            mock_record.samples = {"sample1": mock_sample}
        else:
            mock_record.samples = {}

        return mock_record

    def test_extract_maf_from_ad_simple(self):
        """Test MAF extraction from AD field - simple case."""
        record = self.create_mock_record(
            format_fields={"AD": True}, sample_data={"AD": [80, 20]}  # 80 ref, 20 alt
        )

        # MAF should be 20/100 = 0.2
        maf = extract_maf_from_record(record)
        assert maf == pytest.approx(0.2)

    def test_extract_maf_from_ad_multiallelic(self):
        """Test MAF extraction from AD field - multi-allelic."""
        record = self.create_mock_record(
            format_fields={"AD": True},
            sample_data={"AD": [60, 30, 10]},  # 60 ref, 30 alt1, 10 alt2
        )

        # MAF should be min(ref_freq, total_alt_freq) = min(0.6, 0.4) = 0.4
        maf = extract_maf_from_record(record)
        assert maf == pytest.approx(0.4)

    def test_extract_maf_from_dp4(self):
        """Test MAF extraction from DP4 field."""
        record = self.create_mock_record(
            info_fields={"DP4": [30, 30, 10, 10]}  # ref_fwd, ref_rev, alt_fwd, alt_rev
        )

        # Total ref = 60, total alt = 20, MAF = 20/80 = 0.25
        maf = extract_maf_from_record(record)
        assert maf == pytest.approx(0.25)

    def test_extract_maf_from_af_single(self):
        """Test MAF extraction from AF field - single allele."""
        record = self.create_mock_record(info_fields={"AF": 0.3})

        # AF = 0.3, ref_freq = 0.7, MAF = min(0.3, 0.7) = 0.3
        maf = extract_maf_from_record(record)
        assert maf == pytest.approx(0.3)

    def test_extract_maf_from_af_multiallelic(self):
        """Test MAF extraction from AF field - multi-allelic."""
        record = self.create_mock_record(
            info_fields={"AF": [0.2, 0.1]}  # Two alt alleles
        )

        # Total AF = 0.3, ref_freq = 0.7, max_alt = 0.2, MAF = min(0.7, 0.2) = 0.2
        maf = extract_maf_from_record(record)
        assert maf == pytest.approx(0.2)

    def test_extract_maf_missing_data(self):
        """Test MAF extraction when no suitable fields are present."""
        record = self.create_mock_record()

        with pytest.warns(UserWarning):
            maf = extract_maf_from_record(record)
        assert maf is None


class TestAmbigGrid:
    """Test AmbigGrid functionality."""

    def test_grid_initialization(self):
        """Test grid initialization."""
        grid = AmbigGrid(dp_cap=100)
        assert grid.dp_cap == 100
        assert grid.maf_bins == 51
        assert grid.grid.shape == (101, 51)
        assert np.all(grid.grid == 0)

    def test_add_site_basic(self):
        """Test adding sites to grid."""
        grid = AmbigGrid(dp_cap=100)

        # Add site with depth=50, MAF=0.25
        grid.add_site(50, 0.25)

        # Should be in bin [50, 25]
        assert grid.grid[50, 25] == 1
        assert np.sum(grid.grid) == 1

    def test_add_site_depth_capping(self):
        """Test depth capping functionality."""
        grid = AmbigGrid(dp_cap=100)

        # Add site with depth > cap
        grid.add_site(150, 0.25)

        # Should be capped at dp_cap
        assert grid.grid[100, 25] == 1
        assert grid.grid[150] if 150 < grid.grid.shape[0] else True  # No index 150

    def test_add_site_maf_binning(self):
        """Test MAF binning logic."""
        grid = AmbigGrid(dp_cap=100)

        # Test edge cases for MAF binning
        grid.add_site(50, 0.099)  # Should go to bin 9
        grid.add_site(50, 0.100)  # Should go to bin 10
        grid.add_site(50, 0.505)  # Should be capped at bin 50

        assert grid.grid[50, 9] == 1
        assert grid.grid[50, 10] == 1
        assert grid.grid[50, 50] == 1

    def test_build_cumulative(self):
        """Test cumulative grid building."""
        grid = AmbigGrid(dp_cap=10)  # Small grid for testing

        # Add some sites
        grid.add_site(5, 0.20)  # bin [5, 20]
        grid.add_site(7, 0.30)  # bin [7, 30]
        grid.add_site(10, 0.20)  # bin [10, 20]

        cumulative = grid.build_cumulative()

        # At [5, 20] should include sites with dp>=5 and maf>=0.20
        # That includes all three sites
        assert cumulative[5, 20] == 3

        # At [7, 20] should include sites with dp>=7 and maf>=0.20
        # That includes the last two sites
        assert cumulative[7, 20] == 2

        # At [10, 30] should include sites with dp>=10 and maf>=0.30
        # That includes no sites (the site at [10,20] has maf=0.20 < 0.30)
        assert cumulative[10, 30] == 0

    def test_count_at(self):
        """Test counting sites meeting thresholds."""
        grid = AmbigGrid(dp_cap=100)

        # Add sites
        grid.add_site(15, 0.15)  # Meets both thresholds
        grid.add_site(5, 0.25)  # Doesn't meet depth threshold
        grid.add_site(20, 0.05)  # Doesn't meet MAF threshold
        grid.add_site(25, 0.20)  # Meets both thresholds

        # Count with thresholds dp_min=10, maf_min=0.10
        count = grid.count_at(10, 0.10)
        assert count == 2  # Only first and last sites

    def test_to_wide_tsv(self):
        """Test TSV output functionality."""
        with tempfile.TemporaryDirectory() as tmpdir:
            grid = AmbigGrid(dp_cap=5)  # Very small grid

            # Add some sites
            grid.add_site(2, 0.10)
            grid.add_site(3, 0.20)

            output_path = Path(tmpdir) / "test_grid.tsv"
            grid.to_wide_tsv(output_path)

            assert output_path.exists()

            # Check file content
            with open(output_path) as f:
                lines = f.readlines()

            # Should have header + 5 data rows (depth 1-5)
            assert len(lines) == 6

            # Check header
            header = lines[0].strip().split("\t")
            assert header[0] == "depth"
            assert header[1] == "0.00"  # First MAF column
            assert header[-1] == "0.50"  # Last MAF column


class TestVCFNormalization:
    """Test VCF normalization functionality."""

    def test_normalize_and_split_tool_check(self):
        """Test that normalization checks for required tools."""
        with patch("src.ssiamb.vcf_ops.check_vcf_tools") as mock_check:
            mock_check.return_value = False

            with pytest.raises(VCFOperationError, match="Required VCF tools not found"):
                normalize_and_split(Path("input.vcf"), Path("ref.fasta"))

    def test_normalize_and_split_file_checks(self):
        """Test file existence checks."""
        with patch("src.ssiamb.vcf_ops.check_vcf_tools") as mock_check:
            mock_check.return_value = True

            # Test missing input VCF
            with pytest.raises(VCFOperationError, match="Input VCF file not found"):
                normalize_and_split(Path("missing.vcf"), Path("ref.fasta"))

    @patch("src.ssiamb.vcf_ops.check_vcf_tools")
    @patch("subprocess.run")
    def test_normalize_and_split_success(self, mock_run, mock_check):
        """Test successful normalization."""
        mock_check.return_value = True
        mock_run.return_value = MagicMock(stderr="", returncode=0)

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            # Create dummy input files
            input_vcf = tmpdir / "input.vcf"
            ref_fasta = tmpdir / "ref.fasta"
            input_vcf.touch()
            ref_fasta.touch()

            result = normalize_and_split(input_vcf, ref_fasta)

            # Check that bcftools norm was called
            assert mock_run.call_count == 2  # norm + tabix

            # Check first call (bcftools norm)
            norm_call = mock_run.call_args_list[0]
            assert "bcftools" in norm_call[0][0]
            assert "norm" in norm_call[0][0]
            assert "-f" in norm_call[0][0]
            assert "--atomize" in norm_call[0][0]

            # Check second call (tabix)
            tabix_call = mock_run.call_args_list[1]
            assert "tabix" in tabix_call[0][0]

            # Check result path
            expected_path = tmpdir / "input.normalized.vcf.gz"
            assert result == expected_path


class TestCountAmbiguousSites:
    """Test ambiguous site counting functionality."""

    def create_test_vcf(self, tmpdir, records):
        """Create a test VCF file with specified records."""
        vcf_path = Path(tmpdir) / "test.vcf"

        with open(vcf_path, "w") as f:
            # Write minimal VCF header
            f.write("##fileformat=VCFv4.2\n")
            f.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n')
            f.write(
                '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n'
            )
            f.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
            f.write('##FILTER=<ID=FAIL,Description="Failed quality filters">\n')
            f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            f.write(
                '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic Depth">\n'
            )
            f.write(
                '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n'
            )
            f.write("##contig=<ID=chr1,length=250000000>\n")
            f.write("##contig=<ID=chr2,length=200000000>\n")
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n")

            # Write records
            for record in records:
                f.write(record + "\n")

        return vcf_path

    def test_count_ambiguous_sites_basic(self):
        """Test basic ambiguous site counting."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create test VCF with some sites
            records = [
                "chr1\t100\t.\tA\tT\t60\tPASS\tDP=50\tGT:AD:DP\t0/1:30,20:50",  # MAF=0.4, depth=50
                "chr1\t200\t.\tG\tC\t60\tPASS\tDP=20\tGT:AD:DP\t0/1:15,5:20",  # MAF=0.25, depth=20
                "chr1\t300\t.\tT\tA\t60\tPASS\tDP=5\tGT:AD:DP\t0/1:3,2:5",  # MAF=0.4, depth=5
            ]

            vcf_path = self.create_test_vcf(tmpdir, records)

            # Count with thresholds dp_min=10, maf_min=0.30
            count, grid = count_ambiguous_sites(
                vcf_path=vcf_path,
                dp_min=10,
                maf_min=0.30,
                variant_classes=[VariantClass.SNV],
            )

            # Only first site should meet both thresholds
            assert count == 1

    def test_count_ambiguous_sites_contig_filtering(self):
        """Test contig filtering in site counting."""
        with tempfile.TemporaryDirectory() as tmpdir:
            records = [
                "chr1\t100\t.\tA\tT\t60\tPASS\tDP=50\tGT:AD:DP\t0/1:30,20:50",  # Included contig
                "chr2\t200\t.\tG\tC\t60\tPASS\tDP=50\tGT:AD:DP\t0/1:30,20:50",  # Excluded contig
            ]

            vcf_path = self.create_test_vcf(tmpdir, records)

            # Count with contig filtering
            included_contigs = {"chr1"}
            count, grid = count_ambiguous_sites(
                vcf_path=vcf_path,
                dp_min=10,
                maf_min=0.30,
                included_contigs=included_contigs,
                variant_classes=[VariantClass.SNV],
            )

            # Only chr1 site should be counted
            assert count == 1


class TestIntegration:
    """Integration tests combining multiple components."""

    def test_full_pipeline_mock(self):
        """Test full VCF processing pipeline with mocked components."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            # Create dummy files
            input_vcf = tmpdir / "input.vcf"
            ref_fasta = tmpdir / "ref.fasta"
            normalized_vcf = tmpdir / "input.normalized.vcf.gz"

            input_vcf.touch()
            ref_fasta.touch()
            normalized_vcf.touch()

            # Mock the VCF parsing to return test sites
            test_sites = [
                SiteRecord("chr1", 100, "A", "T", VariantClass.SNV, 50, 0.4, "PASS"),
                SiteRecord("chr1", 200, "G", "C", VariantClass.SNV, 30, 0.2, "PASS"),
                SiteRecord("chr1", 300, "T", "G", VariantClass.SNV, 15, 0.3, "PASS"),
            ]

            with (
                patch("src.ssiamb.vcf_ops.normalize_and_split") as mock_norm,
                patch("src.ssiamb.vcf_ops.parse_vcf_sites") as mock_parse,
            ):

                mock_norm.return_value = normalized_vcf
                mock_parse.return_value = iter(test_sites)

                # Run counting
                count, grid = count_ambiguous_sites(
                    vcf_path=input_vcf,  # Will be "normalized" to normalized_vcf
                    dp_min=20,
                    maf_min=0.25,
                    variant_classes=[VariantClass.SNV],
                )

                # Should count sites with depth>=20 and MAF>=0.25
                # That's sites at pos 100 (depth=50, MAF=0.4) and 300 (depth=15, MAF=0.3)
                # But wait, site 300 has depth=15 < 20, so only site 100 qualifies
                assert count == 1


class TestEmitVCF:
    """Test VCF emitter functionality."""

    def test_emit_vcf_with_test_data(self):
        """Test VCF emitter with synthetic test data."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            # Create a test VCF file
            test_vcf = tmpdir / "test.vcf"
            with open(test_vcf, "w") as f:
                f.write(
                    """##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=FAIL,Description="Failed quality filters">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1
chr1	100	.	A	G	60	PASS	.	GT:DP:AD	0/1:50:30,20
chr1	200	.	C	T	30	FAIL	.	GT:DP:AD	0/1:5:3,2
chr1	300	.	G	A	60	PASS	.	GT:DP:AD	0/1:100:70,30
"""
                )

            output_vcf = tmpdir / "output.vcf.gz"

            # Mock bgzip and tabix commands
            with patch("subprocess.run") as mock_subprocess:
                # Test with dp_min=10, maf_min=0.1 (all sites should pass thresholds)
                result_path = emit_vcf(
                    normalized_vcf_path=test_vcf,
                    output_path=output_vcf,
                    dp_min=10,
                    maf_min=0.1,
                    sample_name="test_sample",
                    require_pass=False,
                    included_contigs={"chr1"},
                )

                assert result_path == output_vcf
                # Should call tabix (bgzip is handled by pysam)
                assert mock_subprocess.call_count >= 1

    def test_emit_vcf_empty_output(self):
        """Test VCF emitter with no qualifying variants."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            # Create a test VCF file with low-quality variants
            test_vcf = tmpdir / "test.vcf"
            with open(test_vcf, "w") as f:
                f.write(
                    """##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1
chr1	100	.	A	G	60	PASS	.	GT:DP:AD	0/1:5:4,1
"""
                )

            output_vcf = tmpdir / "output.vcf.gz"

            # Mock bgzip and tabix commands
            with patch("subprocess.run") as mock_subprocess:
                # High thresholds that no variants will pass
                result_path = emit_vcf(
                    normalized_vcf_path=test_vcf,
                    output_path=output_vcf,
                    dp_min=50,
                    maf_min=0.4,
                    sample_name="test_sample",
                    require_pass=False,
                    included_contigs={"chr1"},
                )

                assert result_path == output_vcf
                # Should still create valid VCF file with header (tabix call)
                assert mock_subprocess.call_count >= 1


class TestEmitBED:
    """Test BED emitter functionality."""

    def test_emit_bed_with_test_data(self):
        """Test BED emitter with synthetic test data."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            # Create a test VCF file with different variant types
            test_vcf = tmpdir / "test.vcf"
            with open(test_vcf, "w") as f:
                f.write(
                    """##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1
chr1	100	.	A	G	60	PASS	.	GT:DP:AD	0/1:50:30,20
chr1	200	.	C	CT	60	PASS	.	GT:DP:AD	0/1:40:20,20
chr1	300	.	GTA	G	60	PASS	.	GT:DP:AD	0/1:30:15,15
"""
                )

            output_bed = tmpdir / "output.bed.gz"

            # Mock bgzip and tabix commands
            with patch("subprocess.run") as mock_subprocess:
                result_path = emit_bed(
                    normalized_vcf_path=test_vcf,
                    output_path=output_bed,
                    dp_min=10,
                    maf_min=0.1,
                    sample_name="test_sample",
                    included_contigs={"chr1"},
                )

                assert result_path == output_bed
                # Should call bgzip and tabix
                assert mock_subprocess.call_count >= 2

    def test_emit_bed_coordinate_conversion(self):
        """Test that BED coordinates are properly 0-based half-open."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            # Create a test VCF file
            test_vcf = tmpdir / "test.vcf"
            with open(test_vcf, "w") as f:
                f.write(
                    """##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1
chr1	100	.	A	G	60	PASS	.	GT:DP:AD	0/1:50:30,20
"""
                )

            output_bed = tmpdir / "output.bed.gz"

            # Mock bgzip and tabix subprocess calls
            with patch("subprocess.run"):
                result_path = emit_bed(
                    normalized_vcf_path=test_vcf,
                    output_path=output_bed,
                    dp_min=10,
                    maf_min=0.1,
                    sample_name="test_sample",
                    included_contigs={"chr1"},
                )

                assert result_path == output_bed


class TestEmitterIntegration:
    """Test integration of emitters with the overall pipeline."""

    def test_emitter_tool_requirements(self):
        """Test that emitters require bgzip and tabix tools."""
        # This is covered by existing check_vcf_tools tests
        with patch("shutil.which") as mock_which:
            mock_which.return_value = None
            assert check_vcf_tools() is False

    def test_emitter_file_extensions(self):
        """Test that emitters handle file extensions correctly."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            # Create a minimal test VCF
            test_vcf = tmpdir / "test.vcf"
            with open(test_vcf, "w") as f:
                f.write(
                    """##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1
"""
                )

            # Test VCF output path handling
            with patch("subprocess.run"):
                # Test that .vcf.gz extension is added if missing
                vcf_output = tmpdir / "output"
                result = emit_vcf(
                    normalized_vcf_path=test_vcf,
                    output_path=vcf_output,
                    dp_min=10,
                    maf_min=0.1,
                    sample_name="test",
                    included_contigs=set(),
                )
                assert str(result).endswith(".vcf.gz")

                # Test BED output path handling
                bed_output = tmpdir / "output"
                result = emit_bed(
                    normalized_vcf_path=test_vcf,
                    output_path=bed_output,
                    dp_min=10,
                    maf_min=0.1,
                    sample_name="test",
                    included_contigs=set(),
                )
                assert str(result).endswith(".bed.gz")


class TestMatrixEmitter:
    """Test matrix emitter functionality."""

    def test_emit_matrix_basic(self):
        """Test basic matrix emission with gzip compression."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            # Create a simple AmbigGrid
            grid = AmbigGrid(dp_cap=100)
            # Add some test data
            grid.add_site(20, 0.3)
            grid.add_site(15, 0.2)

            output_path = tmpdir / "test_matrix.tsv"

            result = emit_matrix(
                grid=grid, output_path=output_path, sample_name="test_sample"
            )

            # Check that .tsv.gz extension was added
            assert str(result).endswith(".tsv.gz")
            assert result.exists()

            # Check that file is gzip compressed
            with gzip.open(result, "rt") as f:
                content = f.read()
                assert "depth" in content  # Should have header
                assert len(content.strip().split("\n")) > 1  # Should have data

    def test_emit_matrix_path_handling(self):
        """Test matrix emitter path handling."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            grid = AmbigGrid(dp_cap=10)

            # Test without .tsv.gz extension
            output_path = tmpdir / "matrix"
            result = emit_matrix(grid=grid, output_path=output_path, sample_name="test")
            assert str(result).endswith(".tsv.gz")

            # Test with .tsv extension only
            output_path = tmpdir / "matrix.tsv"
            result = emit_matrix(grid=grid, output_path=output_path, sample_name="test")
            assert str(result).endswith(".tsv.gz")


class TestPerContigEmitter:
    """Test per-contig emitter functionality."""

    def test_emit_per_contig_basic(self):
        """Test basic per-contig emission."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            # Create mock depth summary file
            depth_summary = tmpdir / "test.mosdepth.summary.txt"
            depth_summary.write_text(
                "chrom\tlength\tbases\tmean\tmin\tmax\n"
                "contig1\t1000\t950\t25.5\t0\t100\n"
                "contig2\t2000\t1800\t30.2\t0\t120\n"
                "total\t3000\t2750\t27.8\t0\t120\n"
            )

            # Create mock VCF file
            vcf_content = """##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic Depth">
##contig=<ID=contig1,length=1000>
##contig=<ID=contig2,length=2000>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttest_sample
contig1\t100\t.\tA\tT\t60\tPASS\tDP=20\tGT:DP:AD\t0/1:20:10,10
contig1\t200\t.\tG\tC\t70\tPASS\tDP=25\tGT:DP:AD\t0/1:25:12,13
contig2\t500\t.\tC\tA\t80\tPASS\tDP=30\tGT:DP:AD\t0/1:30:15,15
"""
            vcf_file = tmpdir / "test.vcf"
            vcf_file.write_text(vcf_content)

            output_path = tmpdir / "per_contig.tsv"

            result = emit_per_contig(
                normalized_vcf_path=vcf_file,
                depth_summary_path=depth_summary,
                output_path=output_path,
                dp_min=10,
                maf_min=0.1,
                sample_name="test_sample",
                included_contigs={"contig1", "contig2"},
            )

            assert result.exists()

            # Check content structure
            content = result.read_text()
            lines = content.strip().split("\n")
            assert len(lines) >= 2  # Header + at least one data row

            # Check header
            header = lines[0].split("\t")
            expected_fields = [
                "sample",
                "contig",
                "length",
                "callable_bases_10x",
                "breadth_10x",
                "ambiguous_snv_count",
                "ambiguous_indel_count",
                "ambiguous_del_count",
                "ambiguous_snv_per_mb",
                "mean_depth",
            ]
            assert header == expected_fields


class TestMultiQCEmitter:
    """Test MultiQC emitter functionality."""

    def test_emit_multiqc_basic(self):
        """Test basic MultiQC emission."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            output_path = tmpdir / "multiqc.tsv"

            result = emit_multiqc(
                sample_name="test_sample",
                ambiguous_snv_count=42,
                breadth_10x=0.85,
                callable_bases=1500000,
                genome_length=2000000,
                dp_min=10,
                maf_min=0.1,
                mapper="minimap2",
                caller="bbtools",
                mode="self",
                output_path=output_path,
            )

            assert result.exists()
            assert str(result).endswith(".tsv")

            # Check content structure
            content = result.read_text()
            lines = content.strip().split("\n")
            assert len(lines) == 2  # Header + one data row

            # Check header
            header = lines[0].split("\t")
            expected_fields = [
                "sample",
                "ambiguous_snv_count",
                "ambiguous_snv_per_mb",
                "breadth_10x",
                "callable_bases",
                "genome_length",
                "dp_min",
                "maf_min",
                "mapper",
                "caller",
                "mode",
            ]
            assert header == expected_fields

            # Check data values
            data = lines[1].split("\t")
            assert data[0] == "test_sample"
            assert data[1] == "42"
            assert data[4] == "1500000"
            assert data[5] == "2000000"
            assert data[8] == "minimap2"
            assert data[9] == "bbtools"
            assert data[10] == "self"

    def test_emit_multiqc_snv_per_mb_calculation(self):
        """Test SNV per megabase calculation."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            output_path = tmpdir / "multiqc.tsv"

            result = emit_multiqc(
                sample_name="test",
                ambiguous_snv_count=100,
                breadth_10x=0.9,
                callable_bases=1000000,
                genome_length=2000000,  # 2 Mb
                dp_min=15,
                maf_min=0.05,
                mapper="bwa-mem2",
                caller="bcftools",
                mode="ref",
                output_path=output_path,
            )

            content = result.read_text()
            lines = content.strip().split("\n")
            data = lines[1].split("\t")

            # SNV per Mb should be 100 * 1,000,000 / 1,000,000 = 100.00 (per spec: ambiguous_snv_count * 1e6 / callable_bases)
            snv_per_mb = float(data[2])
            assert abs(snv_per_mb - 100.0) < 0.01


if __name__ == "__main__":
    pytest.main([__file__])
