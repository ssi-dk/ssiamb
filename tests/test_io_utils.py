"""
Tests for io_utils.py functions.

Tests I/O functionality including sample name validation/inference, MD5 calculation, and TSV writing.
"""

import pytest
import tempfile
import hashlib
from pathlib import Path
from unittest.mock import patch

from src.ssiamb.io_utils import (
    validate_sample_name,
    infer_sample_name,
    compute_md5,
    write_tsv_summary,
    write_tsv_to_stdout,
    file_lock,
    SampleNameError,
    TSVWriteError,
    TSVMode,
)
from src.ssiamb.models import SummaryRow


class TestSampleNameValidation:
    """Test sample name validation functionality."""

    @pytest.mark.parametrize(
        "valid_name",
        [
            "sample1",
            "test_sample",
            "sample-123",
            "S123",
            "abc_def_123",
            "a",  # Minimum valid length
            "a" * 64,  # Maximum valid length
            "test._-123",  # Special characters that should be allowed
        ],
    )
    def test_validate_sample_name_valid(self, valid_name):
        """Test validation of valid sample names."""
        result = validate_sample_name(valid_name)
        assert result == valid_name

    @pytest.mark.parametrize(
        "invalid_name,reason",
        [
            ("", "Empty"),
            ("sample with spaces", "Spaces"),
            ("sample/slash", "Slash"),
            ("sample\\backslash", "Backslash"),
            ("sample:colon", "Colon"),
            ("sample*asterisk", "Asterisk"),
            ("sample?question", "Question mark"),
            ('sample"quote', "Quote"),
            ("sample<less", "Less than"),
            ("sample>greater", "Greater than"),
            ("sample|pipe", "Pipe"),
            ("a" * 65, "Too long (65 chars)"),
        ],
    )
    def test_validate_sample_name_invalid(self, invalid_name, reason):
        """Test validation of invalid sample names."""
        with pytest.raises(SampleNameError):
            validate_sample_name(invalid_name)


class TestSampleNameInference:
    """Test sample name inference from file paths."""

    def test_infer_sample_name_paired_reads(self):
        """Test sample name inference from paired reads."""
        r1 = Path("sample123_R1.fastq.gz")
        r2 = Path("sample123_R2.fastq.gz")

        inferred = infer_sample_name(r1, r2)
        assert inferred == "sample123"

    def test_infer_sample_name_single_read(self):
        """Test sample name inference from single read file."""
        r1 = Path("sample456_R1.fastq")

        inferred = infer_sample_name(r1)
        assert inferred == "sample456"

    def test_infer_sample_name_complex_patterns(self):
        """Test sample name inference with complex filename patterns."""
        test_cases = [
            (Path("SRR123456_1.fastq.gz"), Path("SRR123456_2.fastq.gz"), "SRR123456"),
            (
                Path("test_sample.R1.fastq.gz"),
                Path("test_sample.R2.fastq.gz"),
                "test_sample",
            ),
        ]

        for r1, r2, expected in test_cases:
            inferred = infer_sample_name(r1, r2)
            assert inferred == expected, f"Expected {expected}, got {inferred}"

        # Test cases that should raise errors due to mismatched stems
        error_cases = [
            (Path("sample.forward.fq"), Path("sample.reverse.fq")),
            (
                Path("data_L001_R1_001.fastq"),
                Path("data_L001_R2_001.fastq"),
            ),  # Doesn't handle multiple suffixes
        ]

        for r1, r2 in error_cases:
            with pytest.raises(SampleNameError):
                infer_sample_name(r1, r2)

    def test_infer_sample_name_no_common_prefix(self):
        """Test sample name inference when files have no common prefix."""
        r1 = Path("fileA_R1.fastq")
        r2 = Path("fileB_R2.fastq")

        # Should raise error since stems don't match
        with pytest.raises(SampleNameError):
            infer_sample_name(r1, r2)

    def test_infer_sample_name_from_vcf(self):
        """Test sample name inference from VCF filename."""
        vcf = Path("test_sample.vcf")

        inferred = infer_sample_name(Path("dummy_R1.fastq"), vcf=vcf)
        # Should use dummy_R1 since vcf isn't used for inference
        assert inferred == "dummy"

    def test_infer_sample_name_from_bam(self):
        """Test sample name inference from BAM filename."""
        bam = Path("sample_mapped.bam")

        inferred = infer_sample_name(Path("test_R1.fastq"), bam=bam)
        # Should use test_R1 since bam isn't used for inference
        assert inferred == "test"


class TestFileOperations:
    """Test file operation utilities."""

    def test_compute_md5_simple_file(self):
        """Test MD5 computation for a simple file."""
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as tmp_file:
            tmp_file.write("Hello, world!")
            tmp_file.flush()
            file_path = Path(tmp_file.name)

        try:
            # Compute MD5
            md5_hash = compute_md5(file_path)

            # Verify it's a valid MD5 hash (32 hex characters)
            assert len(md5_hash) == 32
            assert all(c in "0123456789abcdef" for c in md5_hash)

            # Verify it matches expected value
            expected_md5 = hashlib.md5(b"Hello, world!").hexdigest()
            assert md5_hash == expected_md5
        finally:
            file_path.unlink()

    def test_compute_md5_empty_file(self):
        """Test MD5 computation for an empty file."""
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as tmp_file:
            file_path = Path(tmp_file.name)

        try:
            md5_hash = compute_md5(file_path)
            expected_md5 = hashlib.md5(b"").hexdigest()
            assert md5_hash == expected_md5
        finally:
            file_path.unlink()

    def test_compute_md5_binary_file(self):
        """Test MD5 computation for a binary file."""
        with tempfile.NamedTemporaryFile(mode="wb", delete=False) as tmp_file:
            binary_data = b"\x00\x01\x02\x03\xff\xfe\xfd"
            tmp_file.write(binary_data)
            tmp_file.flush()
            file_path = Path(tmp_file.name)

        try:
            md5_hash = compute_md5(file_path)
            expected_md5 = hashlib.md5(binary_data).hexdigest()
            assert md5_hash == expected_md5
        finally:
            file_path.unlink()

    def test_compute_md5_large_file(self):
        """Test MD5 computation for a larger file."""
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as tmp_file:
            # Write 1MB of data
            data = "A" * (1024 * 1024)
            tmp_file.write(data)
            tmp_file.flush()
            file_path = Path(tmp_file.name)

        try:
            md5_hash = compute_md5(file_path)
            expected_md5 = hashlib.md5(data.encode()).hexdigest()
            assert md5_hash == expected_md5
        finally:
            file_path.unlink()

    def test_compute_md5_nonexistent_file(self):
        """Test MD5 computation for nonexistent file."""
        nonexistent_path = Path("/tmp/nonexistent_file_12345.txt")

        with pytest.raises(FileNotFoundError):
            compute_md5(nonexistent_path)


class TestFileLocking:
    """Test file locking functionality."""

    def test_file_lock_context_manager(self):
        """Test file locking as context manager."""
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as tmp_file:
            file_path = Path(tmp_file.name)

        try:
            # Use file lock
            with file_lock(file_path):
                # Lock should be active
                lock_path = file_path.with_suffix(file_path.suffix + ".lock")
                assert lock_path.exists()

            # Lock should be released
            assert not lock_path.exists()
        finally:
            file_path.unlink()
            # Clean up any remaining lock file
            lock_path = file_path.with_suffix(file_path.suffix + ".lock")
            if lock_path.exists():
                lock_path.unlink()

    def test_file_lock_creates_lock_file(self):
        """Test that file lock creates expected lock file."""
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as tmp_file:
            file_path = Path(tmp_file.name)

        try:
            lock_path = file_path.with_suffix(file_path.suffix + ".lock")

            with file_lock(file_path):
                assert lock_path.exists()
                # Lock file might be empty or contain PID - just verify it exists
                assert lock_path.is_file()
        finally:
            file_path.unlink()
            if lock_path.exists():
                lock_path.unlink()


class TestTSVWriting:
    """Test TSV writing functionality."""

    def test_write_tsv_summary_basic(self):
        """Test basic TSV writing functionality."""
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False
        ) as tmp_file:
            output_path = Path(tmp_file.name)

        # Create test data using correct SummaryRow fields
        rows = [
            SummaryRow(
                sample="test1",
                mode="self",
                mapper="minimap2",
                caller="bbtools",
                dp_min=10,
                maf_min=0.1,
                dp_cap=100,
                denom_policy="exclude_dups",
                callable_bases=1000,
                genome_length=1200,
                breadth_10x=0.8333,
                ambiguous_snv_count=45,
                ambiguous_snv_per_mb=37.5,
                ambiguous_indel_count=5,
                ambiguous_del_count=2,
                ref_label="E.coli|NC_000913",
                ref_accession="NC_000913",
                bracken_species="NA",
                bracken_frac=0.0,
                bracken_reads=0,
                alias_used="NA",
                reused_bam=False,
                reused_vcf=False,
                runtime_sec=120.5,
                tool_version="1.0.0",
            )
        ]

        # Write TSV
        write_tsv_summary(output_path, rows, TSVMode.OVERWRITE)

        # Verify file was created and has content
        assert output_path.exists()
        content = output_path.read_text()
        assert "sample" in content  # Header
        assert "test1" in content  # Data

        # Clean up
        output_path.unlink()

    def test_write_tsv_summary_overwrite_mode(self):
        """Test TSV writing in overwrite mode."""
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False
        ) as tmp_file:
            output_path = Path(tmp_file.name)
            # Write initial content
            tmp_file.write("initial content")
            tmp_file.flush()

        rows = [
            SummaryRow(
                sample="overwrite_test",
                mode="self",
                mapper="minimap2",
                caller="bbtools",
                dp_min=10,
                maf_min=0.1,
                dp_cap=100,
                denom_policy="exclude_dups",
                callable_bases=100,
                genome_length=100,
                breadth_10x=1.0,
                ambiguous_snv_count=10,
                ambiguous_snv_per_mb=100.0,
                ambiguous_indel_count=0,
                ambiguous_del_count=0,
                ref_label="E.coli|NC_000913",
                ref_accession="NC_000913",
                bracken_species="NA",
                bracken_frac=0.0,
                bracken_reads=0,
                alias_used="NA",
                reused_bam=False,
                reused_vcf=False,
                runtime_sec=60.0,
                tool_version="1.0.0",
            )
        ]

        # Overwrite the file
        write_tsv_summary(output_path, rows, TSVMode.OVERWRITE)

        # Verify original content is gone
        content = output_path.read_text()
        assert "initial content" not in content
        assert "overwrite_test" in content

        # Clean up
        output_path.unlink()

    def test_write_tsv_summary_append_mode(self):
        """Test TSV writing in append mode."""
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False
        ) as tmp_file:
            output_path = Path(tmp_file.name)
            # Write TSV header and one row
            tmp_file.write(
                "sample\tmode\tmapper\tcaller\tdp_min\tmaf_min\tdp_cap\tdenom_policy\tcallable_bases\tgenome_length\tbreadth_10x\tambiguous_snv_count\tambiguous_snv_per_mb\tambiguous_indel_count\tambiguous_del_count\tref_label\tref_accession\tbracken_species\tbracken_frac\tbracken_reads\talias_used\treused_bam\treused_vcf\truntime_sec\ttool_version\n"
            )
            tmp_file.write(
                "existing\tself\tminimap2\tbbtools\t10\t0.1\t100\texclude_dups\t500\t500\t1.0\t5\t10.0\t0\t0\tE.coli|NC_000913\tNC_000913\tNA\t0.0\t0\tNA\tFalse\tFalse\t50.0\t1.0.0\n"
            )
            tmp_file.flush()

        rows = [
            SummaryRow(
                sample="appended_test",
                mode="self",
                mapper="minimap2",
                caller="bbtools",
                dp_min=10,
                maf_min=0.1,
                dp_cap=100,
                denom_policy="exclude_dups",
                callable_bases=200,
                genome_length=200,
                breadth_10x=1.0,
                ambiguous_snv_count=20,
                ambiguous_snv_per_mb=100.0,
                ambiguous_indel_count=0,
                ambiguous_del_count=0,
                ref_label="E.coli|NC_000913",
                ref_accession="NC_000913",
                bracken_species="NA",
                bracken_frac=0.0,
                bracken_reads=0,
                alias_used="NA",
                reused_bam=False,
                reused_vcf=False,
                runtime_sec=75.0,
                tool_version="1.0.0",
            )
        ]

        # Append to the file
        write_tsv_summary(output_path, rows, TSVMode.APPEND)

        # Verify both old and new content exist
        content = output_path.read_text()
        assert "existing" in content
        assert "appended_test" in content

        # Clean up
        output_path.unlink()

    def test_write_tsv_summary_fail_mode(self):
        """Test TSV writing in fail mode when file exists."""
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False
        ) as tmp_file:
            output_path = Path(tmp_file.name)
            tmp_file.write("existing content")
            tmp_file.flush()

        rows = [
            SummaryRow(
                sample="fail_test",
                mode="self",
                mapper="minimap2",
                caller="bbtools",
                dp_min=10,
                maf_min=0.1,
                dp_cap=100,
                denom_policy="exclude_dups",
                callable_bases=100,
                genome_length=100,
                breadth_10x=1.0,
                ambiguous_snv_count=5,
                ambiguous_snv_per_mb=50.0,
                ambiguous_indel_count=0,
                ambiguous_del_count=0,
                ref_label="E.coli|NC_000913",
                ref_accession="NC_000913",
                bracken_species="NA",
                bracken_frac=0.0,
                bracken_reads=0,
                alias_used="NA",
                reused_bam=False,
                reused_vcf=False,
                runtime_sec=30.0,
                tool_version="1.0.0",
            )
        ]

        # Should raise error when file exists and mode is FAIL
        with pytest.raises(FileExistsError):
            write_tsv_summary(output_path, rows, TSVMode.FAIL)

        # Clean up
        output_path.unlink()

    def test_write_tsv_summary_empty_rows(self):
        """Test TSV writing with empty rows list."""
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False
        ) as tmp_file:
            output_path = Path(tmp_file.name)

        # Should raise error for empty rows
        with pytest.raises(TSVWriteError):
            write_tsv_summary(output_path, [], TSVMode.OVERWRITE)

        # Clean up
        output_path.unlink()


class TestTSVStdout:
    """Test TSV stdout writing functionality."""

    @patch("sys.stdout.write")
    def test_write_tsv_to_stdout_basic(self, mock_write):
        """Test writing TSV to stdout."""
        rows = [
            SummaryRow(
                sample="stdout_test",
                mode="self",
                mapper="minimap2",
                caller="bbtools",
                dp_min=10,
                maf_min=0.1,
                dp_cap=100,
                denom_policy="exclude_dups",
                callable_bases=300,
                genome_length=300,
                breadth_10x=1.0,
                ambiguous_snv_count=25,
                ambiguous_snv_per_mb=83.33,
                ambiguous_indel_count=0,
                ambiguous_del_count=0,
                ref_label="E.coli|NC_000913",
                ref_accession="NC_000913",
                bracken_species="NA",
                bracken_frac=0.0,
                bracken_reads=0,
                alias_used="NA",
                reused_bam=False,
                reused_vcf=False,
                runtime_sec=90.0,
                tool_version="1.0.0",
            )
        ]

        write_tsv_to_stdout(rows)

        # Verify stdout.write was called
        assert mock_write.called

        # Reconstruct the written content
        written_content = "".join(call.args[0] for call in mock_write.call_args_list)
        assert "sample" in written_content  # Header
        assert "stdout_test" in written_content  # Data

    @patch("sys.stdout.write")
    def test_write_tsv_to_stdout_empty(self, mock_write):
        """Test writing empty TSV to stdout."""
        write_tsv_to_stdout([])

        # With empty list, nothing should be written
        assert not mock_write.called


class TestErrorHandling:
    """Test error handling in I/O operations."""

    def test_sample_name_error_creation(self):
        """Test SampleNameError exception creation."""
        error = SampleNameError("Invalid name")
        assert str(error) == "Invalid name"
        assert isinstance(error, ValueError)

    def test_tsv_write_error_creation(self):
        """Test TSVWriteError exception creation."""
        error = TSVWriteError("Write failed")
        assert str(error) == "Write failed"
        assert isinstance(error, Exception)

    @patch("builtins.open", side_effect=PermissionError("Access denied"))
    def test_compute_md5_permission_error(self, mock_open):
        """Test MD5 computation with permission error."""
        file_path = Path("/protected/file.txt")

        with pytest.raises(PermissionError):
            compute_md5(file_path)

    def test_write_tsv_summary_invalid_path(self):
        """Test TSV writing with invalid output path."""
        invalid_path = Path("/root/protected/file.tsv")  # Likely no permission
        rows = [
            SummaryRow(
                sample="invalid_path_test",
                mode="self",
                mapper="minimap2",
                caller="bbtools",
                dp_min=10,
                maf_min=0.1,
                dp_cap=100,
                denom_policy="exclude_dups",
                callable_bases=100,
                genome_length=100,
                breadth_10x=1.0,
                ambiguous_snv_count=5,
                ambiguous_snv_per_mb=50.0,
                ambiguous_indel_count=0,
                ambiguous_del_count=0,
                ref_label="E.coli|NC_000913",
                ref_accession="NC_000913",
                bracken_species="NA",
                bracken_frac=0.0,
                bracken_reads=0,
                alias_used="NA",
                reused_bam=False,
                reused_vcf=False,
                runtime_sec=30.0,
                tool_version="1.0.0",
            )
        ]

        # Should handle permission/path errors gracefully
        with pytest.raises((PermissionError, OSError, TSVWriteError)):
            write_tsv_summary(invalid_path, rows, TSVMode.OVERWRITE)
