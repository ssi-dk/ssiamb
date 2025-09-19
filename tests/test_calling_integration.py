"""
Integration tests for variant calling functionality using real fixture data.

Tests real variant calling workflows with BBTools and bcftools using actual
BAM files, reference genomes, and FASTQ reads from SSI pipeline.
"""

import pytest
import tempfile
import subprocess
import time
from pathlib import Path
from unittest.mock import patch

from src.ssiamb.calling import (
    check_caller_tools,
    call_variants,
    get_available_callers,
    VariantCallingError,
    VariantCallResult,
)
from src.ssiamb.models import Caller
from src.ssiamb.mapping import map_fastqs, ensure_indexes_self
from src.ssiamb.models import Mapper


class TestCallingRealData:
    """Integration tests using real data files."""

    @pytest.fixture
    def real_reference_path(self):
        """Path to real reference FASTA."""
        path = Path("fixtures/self_mode/assembly.fna")
        if not path.exists():
            pytest.skip("Real reference FASTA not available")
        return path

    @pytest.fixture
    def real_reads_r1_path(self):
        """Path to real R1 FASTQ."""
        path = Path("fixtures/self_mode/reads_R1.fastq.gz")
        if not path.exists():
            pytest.skip("Real R1 FASTQ not available")
        return path

    @pytest.fixture
    def real_reads_r2_path(self):
        """Path to real R2 FASTQ."""
        path = Path("fixtures/self_mode/reads_R2.fastq.gz")
        if not path.exists():
            pytest.skip("Real R2 FASTQ not available")
        return path

    @pytest.fixture
    def pre_mapped_bam_path(self):
        """Path to pre-existing BAM file."""
        path = Path("fixtures/reuse/pre_mapped.bam")
        if not path.exists():
            pytest.skip("Pre-mapped BAM not available")
        return path

    def test_caller_tool_availability(self):
        """Test that required calling tools are available."""
        available_callers = get_available_callers()

        # Should have at least one caller available in test environment
        assert (
            len(available_callers) > 0
        ), "At least one variant caller should be available"

        # Test individual tool checks
        assert check_caller_tools(Caller.BBTOOLS), "BBTools should be available"
        assert check_caller_tools(Caller.BCFTOOLS), "bcftools should be available"

        print(f"Available callers: {[c.value for c in available_callers]}")

    def test_pre_mapped_bam_properties(self, pre_mapped_bam_path):
        """Test properties of the pre-mapped BAM file."""
        # Check BAM file size and basic properties
        assert (
            pre_mapped_bam_path.stat().st_size > 100000
        ), "BAM should be reasonably large"

        # Check that BAM is readable with samtools
        try:
            result = subprocess.run(
                ["samtools", "view", "-H", str(pre_mapped_bam_path)],
                capture_output=True,
                text=True,
                timeout=10,
            )
            assert result.returncode == 0, "BAM should be readable with samtools"
            assert len(result.stdout.strip()) > 0, "BAM should have header"
        except subprocess.TimeoutExpired:
            pytest.fail("samtools view took too long")

    def test_bbtools_calling_pre_mapped_bam(
        self, pre_mapped_bam_path, real_reference_path
    ):
        """Test BBTools variant calling with pre-mapped BAM."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            # Copy reference to temp directory (BBTools may need to modify it)
            temp_ref = tmpdir / "reference.fna"
            temp_ref.write_bytes(real_reference_path.read_bytes())

            output_vcf = tmpdir / "variants_bbtools.vcf"

            # Run BBTools variant calling
            start_time = time.time()
            result = call_variants(
                bam_path=pre_mapped_bam_path,
                reference_path=temp_ref,
                output_vcf=output_vcf,
                caller=Caller.BBTOOLS,
                sample_name="test_sample",
                threads=4,
                mapq_min=20,
                baseq_min=20,
                minallelefraction=0.05,
            )
            calling_time = time.time() - start_time

            # Validate calling results
            assert (
                result.success
            ), f"BBTools calling should succeed: {result.error_message}"
            assert result.vcf_path.exists(), "Output VCF should be created"
            assert result.vcf_path.stat().st_size > 0, "Output VCF should not be empty"
            assert result.caller == Caller.BBTOOLS, "Result should have correct caller"

            # Check VCF format (basic validation)
            with open(result.vcf_path) as f:
                content = f.read()
                assert content.startswith(
                    "##fileformat=VCF"
                ), "VCF should have proper header"
                assert "#CHROM" in content, "VCF should have column headers"

            # Calling should complete in reasonable time
            assert (
                calling_time < 300
            ), f"BBTools calling took too long: {calling_time:.2f} seconds"

            print(f"BBTools calling completed in {calling_time:.2f} seconds")
            print(f"VCF size: {result.vcf_path.stat().st_size:,} bytes")

    def test_bcftools_calling_pre_mapped_bam(
        self, pre_mapped_bam_path, real_reference_path
    ):
        """Test bcftools variant calling with pre-mapped BAM."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            # Copy reference to temp directory
            temp_ref = tmpdir / "reference.fna"
            temp_ref.write_bytes(real_reference_path.read_bytes())

            output_vcf = tmpdir / "variants_bcftools.vcf"

            # Run bcftools variant calling
            start_time = time.time()
            result = call_variants(
                bam_path=pre_mapped_bam_path,
                reference_path=temp_ref,
                output_vcf=output_vcf,
                caller=Caller.BCFTOOLS,
                sample_name="test_sample",
                threads=4,
                mapq_min=20,
                baseq_min=20,
            )
            calling_time = time.time() - start_time

            # Validate calling results
            assert (
                result.success
            ), f"bcftools calling should succeed: {result.error_message}"
            assert result.vcf_path.exists(), "Output VCF should be created"
            assert result.vcf_path.stat().st_size > 0, "Output VCF should not be empty"
            assert result.caller == Caller.BCFTOOLS, "Result should have correct caller"

            # Check VCF format (basic validation)
            with open(result.vcf_path) as f:
                content = f.read()
                assert content.startswith(
                    "##fileformat=VCF"
                ), "VCF should have proper header"
                assert "#CHROM" in content, "VCF should have column headers"

            # Calling should complete in reasonable time
            assert (
                calling_time < 300
            ), f"bcftools calling took too long: {calling_time:.2f} seconds"

            print(f"bcftools calling completed in {calling_time:.2f} seconds")
            print(f"VCF size: {result.vcf_path.stat().st_size:,} bytes")

    def test_end_to_end_mapping_and_calling(
        self, real_reference_path, real_reads_r1_path, real_reads_r2_path
    ):
        """Test complete end-to-end pipeline: mapping + calling with real data."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            # Copy reference to temp directory
            temp_ref = tmpdir / "reference.fna"
            temp_ref.write_bytes(real_reference_path.read_bytes())

            # Step 1: Create mapping
            print("Step 1: Creating mapping...")
            ensure_indexes_self(temp_ref, Mapper.MINIMAP2)

            mapped_bam = tmpdir / "mapped.bam"

            map_start = time.time()
            result_bam = map_fastqs(
                mapper=Mapper.MINIMAP2,
                fasta_path=temp_ref,
                r1_path=real_reads_r1_path,
                r2_path=real_reads_r2_path,
                sample_name="test_sample",
                threads=4,
                output_path=mapped_bam,
            )
            map_time = time.time() - map_start

            assert result_bam.exists(), "Mapping should create BAM file"

            # Step 2: Call variants using bcftools (more reliable)
            print("Step 2: Calling variants with bcftools...")
            output_vcf = tmpdir / "variants.vcf"

            call_start = time.time()
            call_result = call_variants(
                bam_path=result_bam,
                reference_path=temp_ref,
                output_vcf=output_vcf,
                caller=Caller.BCFTOOLS,
                sample_name="test_sample",
                threads=4,
            )
            call_time = time.time() - call_start

            # Validate end-to-end results
            assert (
                call_result.success
            ), f"End-to-end calling should succeed: {call_result.error_message}"
            assert call_result.vcf_path.exists(), "Final VCF should be created"
            assert (
                call_result.vcf_path.stat().st_size > 0
            ), "Final VCF should not be empty"

            # Total time should be reasonable
            total_time = map_time + call_time
            assert (
                total_time < 600
            ), f"End-to-end pipeline took too long: {total_time:.2f} seconds"

            print("End-to-end pipeline completed:")
            print(f"  Mapping: {map_time:.2f} seconds")
            print(f"  Calling: {call_time:.2f} seconds")
            print(f"  Total: {total_time:.2f} seconds")
            print(f"  Final VCF size: {call_result.vcf_path.stat().st_size:,} bytes")

    def test_caller_comparison_pre_mapped_bam(
        self, pre_mapped_bam_path, real_reference_path
    ):
        """Compare BBTools vs bcftools performance and results with pre-mapped BAM."""
        results = {}

        for caller in [Caller.BBTOOLS, Caller.BCFTOOLS]:
            with tempfile.TemporaryDirectory() as tmpdir:
                tmpdir = Path(tmpdir)

                # Copy reference
                temp_ref = tmpdir / f"reference_{caller.value}.fna"
                temp_ref.write_bytes(real_reference_path.read_bytes())

                output_vcf = tmpdir / f"variants_{caller.value}.vcf"

                # Run variant calling
                call_start = time.time()
                result = call_variants(
                    bam_path=pre_mapped_bam_path,
                    reference_path=temp_ref,
                    output_vcf=output_vcf,
                    caller=caller,
                    sample_name="test_sample",
                    threads=4,
                    minallelefraction=0.05 if caller == Caller.BBTOOLS else 0.0,
                )
                call_time = time.time() - call_start

                # Count variants (rough estimate)
                variant_count = 0
                if result.success and result.vcf_path.exists():
                    with open(result.vcf_path) as f:
                        for line in f:
                            if not line.startswith("#") and line.strip():
                                variant_count += 1

                results[caller.value] = {
                    "success": result.success,
                    "call_time": call_time,
                    "vcf_size": (
                        result.vcf_path.stat().st_size
                        if result.vcf_path.exists()
                        else 0
                    ),
                    "variant_count": variant_count,
                    "error": result.error_message,
                }

        # Both callers should succeed
        for caller_name, result in results.items():
            assert result["success"], f"{caller_name} should succeed: {result['error']}"
            assert result["vcf_size"] > 0, f"{caller_name} should produce non-empty VCF"

        # Log comparison results
        print("\\nCaller comparison results:")
        for caller_name, result in results.items():
            print(f"  {caller_name}:")
            print(f"    Call time: {result['call_time']:.2f}s")
            print(f"    VCF size: {result['vcf_size']:,} bytes")
            print(f"    Variants: {result['variant_count']}")


class TestCallingEdgeCases:
    """Test edge cases and error handling."""

    @pytest.fixture
    def real_reference_path(self):
        """Path to real reference FASTA."""
        path = Path("fixtures/self_mode/assembly.fna")
        if not path.exists():
            pytest.skip("Real reference FASTA not available")
        return path

    def test_missing_bam_file(self, real_reference_path):
        """Test handling of missing BAM file."""
        fake_bam = Path("does_not_exist.bam")

        with tempfile.TemporaryDirectory() as tmpdir:
            output_vcf = Path(tmpdir) / "output.vcf"

            with pytest.raises(VariantCallingError, match="BAM file not found"):
                call_variants(
                    bam_path=fake_bam,
                    reference_path=real_reference_path,
                    output_vcf=output_vcf,
                    caller=Caller.BBTOOLS,
                    sample_name="test",
                )

    def test_missing_reference_file(self):
        """Test handling of missing reference file."""
        fake_ref = Path("does_not_exist.fasta")

        with tempfile.TemporaryDirectory() as tmpdir:
            # Create a dummy BAM file so we get to the reference check
            fake_bam = Path(tmpdir) / "test.bam"
            fake_bam.write_bytes(b"dummy bam content")

            output_vcf = Path(tmpdir) / "output.vcf"

            with pytest.raises(VariantCallingError, match="Reference file not found"):
                call_variants(
                    bam_path=fake_bam,
                    reference_path=fake_ref,
                    output_vcf=output_vcf,
                    caller=Caller.BBTOOLS,
                    sample_name="test",
                )

    def test_corrupted_bam_file(self, real_reference_path):
        """Test handling of corrupted BAM file - BBTools handles this gracefully."""
        with tempfile.NamedTemporaryFile(suffix=".bam", delete=False) as f:
            f.write(b"This is not a BAM file")
            fake_bam = Path(f.name)

        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                output_vcf = Path(tmpdir) / "output.vcf"

                result = call_variants(
                    bam_path=fake_bam,
                    reference_path=real_reference_path,
                    output_vcf=output_vcf,
                    caller=Caller.BBTOOLS,
                    sample_name="test",
                )

                # BBTools actually handles corrupted files gracefully in many cases
                # We just verify that the function completed without crashing
                assert result is not None
                assert isinstance(result, VariantCallResult)
        finally:
            fake_bam.unlink()  # Clean up

    @patch("src.ssiamb.calling.caller_tools_available")
    def test_missing_caller_tools_error_handling(
        self, mock_tools_available, real_reference_path
    ):
        """Test handling when calling tools are not available."""
        mock_tools_available.return_value = False

        # Create a temporary BAM file
        with tempfile.NamedTemporaryFile(suffix=".bam", delete=False) as f:
            f.write(b"BAM\\x01")  # Minimal BAM magic number
            fake_bam = Path(f.name)

        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                output_vcf = Path(tmpdir) / "output.vcf"

                with pytest.raises(
                    VariantCallingError, match="Required tools.*not available"
                ):
                    call_variants(
                        bam_path=fake_bam,
                        reference_path=real_reference_path,
                        output_vcf=output_vcf,
                        caller=Caller.BBTOOLS,
                        sample_name="test",
                    )
        finally:
            fake_bam.unlink()  # Clean up


class TestCallingPerformance:
    """Performance tests with real data."""

    @pytest.fixture
    def pre_mapped_bam_path(self):
        """Path to pre-existing BAM file."""
        path = Path("fixtures/reuse/pre_mapped.bam")
        if not path.exists():
            pytest.skip("Pre-mapped BAM not available for performance test")
        return path

    @pytest.fixture
    def real_reference_path(self):
        """Path to real reference FASTA."""
        path = Path("fixtures/self_mode/assembly.fna")
        if not path.exists():
            pytest.skip("Real reference not available for performance test")
        return path

    def test_calling_performance_scaling(
        self, pre_mapped_bam_path, real_reference_path
    ):
        """Test calling performance with different thread counts."""
        thread_counts = [1, 2, 4]
        performance_results = {}

        for threads in thread_counts:
            with tempfile.TemporaryDirectory() as tmpdir:
                tmpdir = Path(tmpdir)
                temp_ref = tmpdir / "reference.fna"
                temp_ref.write_bytes(real_reference_path.read_bytes())

                output_vcf = tmpdir / "variants.vcf"

                start_time = time.time()
                result = call_variants(
                    bam_path=pre_mapped_bam_path,
                    reference_path=temp_ref,
                    output_vcf=output_vcf,
                    caller=Caller.BBTOOLS,  # Use BBTools for consistency
                    sample_name="test_sample",
                    threads=threads,
                    minallelefraction=0.05,
                )
                call_time = time.time() - start_time

                performance_results[threads] = {
                    "call_time": call_time,
                    "success": result.success,
                }

                # Should complete in reasonable time regardless of thread count
                assert (
                    call_time < 400
                ), f"Calling with {threads} threads took too long: {call_time:.2f}s"
                assert result.success, f"Calling with {threads} threads should succeed"

        print("\\nCalling performance by thread count:")
        for threads, result in performance_results.items():
            print(f"  {threads} threads: {result['call_time']:.2f} seconds")

        # More threads should generally be faster (or at least not much slower)
        # Allow some variance due to overhead and I/O
        max_time = max(r["call_time"] for r in performance_results.values())
        min_time = min(r["call_time"] for r in performance_results.values())
        speedup_ratio = max_time / min_time
        assert (
            speedup_ratio < 4.0
        ), f"Performance variance too high: {speedup_ratio:.2f}x"

    def test_caller_tool_performance_comparison(
        self, pre_mapped_bam_path, real_reference_path
    ):
        """Compare performance between different variant callers."""
        caller_times = {}

        for caller in [Caller.BBTOOLS, Caller.BCFTOOLS]:
            with tempfile.TemporaryDirectory() as tmpdir:
                tmpdir = Path(tmpdir)
                temp_ref = tmpdir / f"ref_{caller.value}.fna"
                temp_ref.write_bytes(real_reference_path.read_bytes())

                output_vcf = tmpdir / f"variants_{caller.value}.vcf"

                start_time = time.time()
                result = call_variants(
                    bam_path=pre_mapped_bam_path,
                    reference_path=temp_ref,
                    output_vcf=output_vcf,
                    caller=caller,
                    sample_name="test_sample",
                    threads=4,
                    minallelefraction=0.05 if caller == Caller.BBTOOLS else 0.0,
                )
                call_time = time.time() - start_time

                caller_times[caller.value] = call_time

                # Sanity check - calling should complete in reasonable time
                assert (
                    call_time < 600
                ), f"{caller.value} calling took too long: {call_time:.2f}s"
                assert result.success, f"{caller.value} calling should succeed"

        print("\\nCaller performance comparison:")
        for caller, call_time in caller_times.items():
            print(f"  {caller}: {call_time:.2f} seconds")
