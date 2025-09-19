"""
Integration tests for mapping functionality using real fixture data.

Tests real mapping workflows with minimap2 and bwa-mem2 using actual
FASTQ reads and reference genome from SSI pipeline.
"""

import pytest
import tempfile
import time
from pathlib import Path
from unittest.mock import patch

from src.ssiamb.mapping import (
    check_external_tools,
    get_index_files,
    indexes_exist,
    ensure_indexes_self,
    map_fastqs,
    calculate_mapping_rate,
    MappingError,
    ExternalToolError,
)
from src.ssiamb.models import Mapper


class TestMappingRealData:
    """Integration tests using real FASTQ and reference data."""

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

    def test_tool_availability(self):
        """Test that required mapping tools are available."""
        tools = check_external_tools()

        # Should have all tools available in test environment
        assert tools["minimap2"], "minimap2 should be available"
        assert tools["bwa-mem2"], "bwa-mem2 should be available"
        assert tools["samtools"], "samtools should be available"

    def test_real_reference_properties(self, real_reference_path):
        """Test properties of the real reference genome."""
        # Check reference file size and structure
        assert (
            real_reference_path.stat().st_size > 1000000
        ), "Reference should be reasonably large"

        # Check that the reference has proper FASTA format
        with open(real_reference_path) as f:
            first_line = f.readline().strip()
            assert first_line.startswith(
                ">"
            ), "Reference should start with FASTA header"

        # Check if index exists
        fai_path = real_reference_path.with_suffix(".fna.fai")
        assert fai_path.exists(), "Reference should have samtools faidx index"

    def test_real_reads_properties(self, real_reads_r1_path, real_reads_r2_path):
        """Test properties of the real FASTQ reads."""
        # Check read file sizes
        assert (
            real_reads_r1_path.stat().st_size > 100000
        ), "R1 reads should be reasonably large"
        assert (
            real_reads_r2_path.stat().st_size > 100000
        ), "R2 reads should be reasonably large"

        # Basic FASTQ format check (first few lines)
        import gzip

        with gzip.open(real_reads_r1_path, "rt") as f:
            header = f.readline().strip()
            assert header.startswith("@"), "FASTQ should start with @ header"

        with gzip.open(real_reads_r2_path, "rt") as f:
            header = f.readline().strip()
            assert header.startswith("@"), "FASTQ should start with @ header"

    def test_minimap2_index_creation_real_data(self, real_reference_path):
        """Test minimap2 index creation with real reference."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Copy reference to temp directory to avoid modifying original
            temp_ref = Path(tmpdir) / "reference.fna"
            temp_ref.write_bytes(real_reference_path.read_bytes())

            # Test index creation
            start_time = time.time()
            ensure_indexes_self(temp_ref, Mapper.MINIMAP2)
            index_time = time.time() - start_time

            # Should create index files
            assert indexes_exist(
                temp_ref, Mapper.MINIMAP2
            ), "Minimap2 index should be created"

            # Index creation should complete in reasonable time
            assert (
                index_time < 30
            ), f"Index creation took too long: {index_time:.2f} seconds"

            # Check that index file exists and has reasonable size
            index_files = get_index_files(temp_ref, Mapper.MINIMAP2)
            for index_file in index_files:
                assert index_file.exists(), f"Index file {index_file} should exist"
                assert (
                    index_file.stat().st_size > 1000
                ), f"Index file {index_file} should not be empty"

            print(f"Minimap2 indexing completed in {index_time:.2f} seconds")

    def test_bwa_mem2_index_creation_real_data(self, real_reference_path):
        """Test bwa-mem2 index creation with real reference."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Copy reference to temp directory
            temp_ref = Path(tmpdir) / "reference.fna"
            temp_ref.write_bytes(real_reference_path.read_bytes())

            # Test index creation
            start_time = time.time()
            ensure_indexes_self(temp_ref, Mapper.BWA_MEM2)
            index_time = time.time() - start_time

            # Should create index files
            assert indexes_exist(
                temp_ref, Mapper.BWA_MEM2
            ), "BWA-MEM2 index should be created"

            # Index creation should complete in reasonable time
            assert (
                index_time < 60
            ), f"Index creation took too long: {index_time:.2f} seconds"

            # Check that all expected index files exist
            index_files = get_index_files(temp_ref, Mapper.BWA_MEM2)
            for index_file in index_files:
                assert index_file.exists(), f"Index file {index_file} should exist"
                assert (
                    index_file.stat().st_size > 0
                ), f"Index file {index_file} should not be empty"

            print(f"BWA-MEM2 indexing completed in {index_time:.2f} seconds")

    def test_minimap2_mapping_real_data(
        self, real_reference_path, real_reads_r1_path, real_reads_r2_path
    ):
        """Test minimap2 mapping with real data."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            # Copy reference to temp directory
            temp_ref = tmpdir / "reference.fna"
            temp_ref.write_bytes(real_reference_path.read_bytes())

            # Ensure indexes exist
            ensure_indexes_self(temp_ref, Mapper.MINIMAP2)

            # Set up output path
            output_bam = tmpdir / "mapped.bam"

            # Run mapping
            start_time = time.time()
            result_bam = map_fastqs(
                mapper=Mapper.MINIMAP2,
                fasta_path=temp_ref,
                r1_path=real_reads_r1_path,
                r2_path=real_reads_r2_path,
                sample_name="test_sample",
                threads=4,
                output_path=output_bam,
            )
            mapping_time = time.time() - start_time

            # Validate mapping results
            assert result_bam.exists(), "Output BAM should be created"
            assert result_bam.stat().st_size > 1000, "Output BAM should not be empty"

            # Check mapping rate
            mapping_rate = calculate_mapping_rate(result_bam)
            assert 0.0 <= mapping_rate <= 1.0, "Mapping rate should be between 0 and 1"
            assert (
                mapping_rate > 0.5
            ), f"Mapping rate should be reasonable: {mapping_rate:.3f}"

            # Mapping should complete in reasonable time
            assert (
                mapping_time < 120
            ), f"Mapping took too long: {mapping_time:.2f} seconds"

            print(f"Minimap2 mapping completed in {mapping_time:.2f} seconds")
            print(f"Mapping rate: {mapping_rate:.3f}")

    def test_bwa_mem2_mapping_real_data(
        self, real_reference_path, real_reads_r1_path, real_reads_r2_path
    ):
        """Test bwa-mem2 mapping with real data."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            # Copy reference to temp directory
            temp_ref = tmpdir / "reference.fna"
            temp_ref.write_bytes(real_reference_path.read_bytes())

            # Ensure indexes exist
            ensure_indexes_self(temp_ref, Mapper.BWA_MEM2)

            # Set up output path
            output_bam = tmpdir / "mapped.bam"

            # Run mapping
            start_time = time.time()
            result_bam = map_fastqs(
                mapper=Mapper.BWA_MEM2,
                fasta_path=temp_ref,
                r1_path=real_reads_r1_path,
                r2_path=real_reads_r2_path,
                sample_name="test_sample",
                threads=4,
                output_path=output_bam,
            )
            mapping_time = time.time() - start_time

            # Validate mapping results
            assert result_bam.exists(), "Output BAM should be created"
            assert result_bam.stat().st_size > 1000, "Output BAM should not be empty"

            # Check mapping rate
            mapping_rate = calculate_mapping_rate(result_bam)
            assert 0.0 <= mapping_rate <= 1.0, "Mapping rate should be between 0 and 1"
            assert (
                mapping_rate > 0.5
            ), f"Mapping rate should be reasonable: {mapping_rate:.3f}"

            # Mapping should complete in reasonable time
            assert (
                mapping_time < 120
            ), f"Mapping took too long: {mapping_time:.2f} seconds"

            print(f"BWA-MEM2 mapping completed in {mapping_time:.2f} seconds")
            print(f"Mapping rate: {mapping_rate:.3f}")

    def test_mapping_rate_calculation_realistic(
        self, real_reference_path, real_reads_r1_path, real_reads_r2_path
    ):
        """Test mapping rate calculation with realistic data."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            # Copy reference and create minimap2 mapping (faster for this test)
            temp_ref = tmpdir / "reference.fna"
            temp_ref.write_bytes(real_reference_path.read_bytes())
            ensure_indexes_self(temp_ref, Mapper.MINIMAP2)

            output_bam = tmpdir / "mapped.bam"

            # Create BAM file
            map_fastqs(
                mapper=Mapper.MINIMAP2,
                fasta_path=temp_ref,
                r1_path=real_reads_r1_path,
                r2_path=real_reads_r2_path,
                sample_name="test_sample",
                threads=4,
                output_path=output_bam,
            )

            # Test mapping rate calculation
            start_time = time.time()
            mapping_rate = calculate_mapping_rate(output_bam)
            calc_time = time.time() - start_time

            # Validate mapping rate
            assert isinstance(mapping_rate, float), "Mapping rate should be a float"
            assert 0.0 <= mapping_rate <= 1.0, "Mapping rate should be between 0 and 1"

            # For real data, expect decent mapping rate
            assert (
                mapping_rate > 0.3
            ), f"Real data should have reasonable mapping rate: {mapping_rate:.3f}"

            # Calculation should be fast
            assert (
                calc_time < 10
            ), f"Mapping rate calculation took too long: {calc_time:.2f} seconds"

            print(
                f"Mapping rate calculation: {mapping_rate:.3f} in {calc_time:.3f} seconds"
            )

    def test_mapper_comparison_real_data(
        self, real_reference_path, real_reads_r1_path, real_reads_r2_path
    ):
        """Compare minimap2 vs bwa-mem2 performance and mapping rates with real data."""
        results = {}

        for mapper in [Mapper.MINIMAP2, Mapper.BWA_MEM2]:
            with tempfile.TemporaryDirectory() as tmpdir:
                tmpdir = Path(tmpdir)

                # Copy reference
                temp_ref = tmpdir / f"reference_{mapper.value}.fna"
                temp_ref.write_bytes(real_reference_path.read_bytes())

                # Create indexes
                index_start = time.time()
                ensure_indexes_self(temp_ref, mapper)
                index_time = time.time() - index_start

                # Run mapping
                output_bam = tmpdir / "mapped.bam"

                map_start = time.time()
                result_bam = map_fastqs(
                    mapper=mapper,
                    fasta_path=temp_ref,
                    r1_path=real_reads_r1_path,
                    r2_path=real_reads_r2_path,
                    sample_name="test_sample",
                    threads=4,
                    output_path=output_bam,
                )
                map_time = time.time() - map_start

                # Calculate mapping rate
                mapping_rate = calculate_mapping_rate(result_bam)

                results[mapper.value] = {
                    "index_time": index_time,
                    "map_time": map_time,
                    "mapping_rate": mapping_rate,
                    "bam_size": result_bam.stat().st_size,
                }

        # Both mappers should produce reasonable results
        for mapper_name, result in results.items():
            assert (
                result["mapping_rate"] > 0.3
            ), f"{mapper_name} should have decent mapping rate"
            assert (
                result["bam_size"] > 1000
            ), f"{mapper_name} should produce non-empty BAM"

        # Log comparison results
        print("\\nMapper comparison results:")
        for mapper_name, result in results.items():
            print(f"  {mapper_name}:")
            print(f"    Index time: {result['index_time']:.2f}s")
            print(f"    Map time: {result['map_time']:.2f}s")
            print(f"    Mapping rate: {result['mapping_rate']:.3f}")
            print(f"    BAM size: {result['bam_size']:,} bytes")


class TestMappingEdgeCases:
    """Test edge cases and error handling."""

    @pytest.fixture
    def real_reference_path(self):
        """Path to real reference FASTA."""
        path = Path("fixtures/self_mode/assembly.fna")
        if not path.exists():
            pytest.skip("Real reference FASTA not available")
        return path

    def test_missing_reference_file(self):
        """Test handling of missing reference file."""
        fake_ref = Path("does_not_exist.fasta")

        with pytest.raises((FileNotFoundError, MappingError)):
            ensure_indexes_self(fake_ref, Mapper.MINIMAP2)

    def test_missing_fastq_files(self, real_reference_path):
        """Test handling of missing FASTQ files."""
        fake_r1 = Path("does_not_exist_R1.fastq.gz")
        fake_r2 = Path("does_not_exist_R2.fastq.gz")

        with tempfile.TemporaryDirectory() as tmpdir:
            output_bam = Path(tmpdir) / "output.bam"

            with pytest.raises((FileNotFoundError, MappingError)):
                map_fastqs(
                    mapper=Mapper.MINIMAP2,
                    fasta_path=real_reference_path,
                    r1_path=fake_r1,
                    r2_path=fake_r2,
                    sample_name="test",
                    output_path=output_bam,
                )

    def test_corrupted_bam_mapping_rate(self):
        """Test mapping rate calculation with corrupted BAM."""
        with tempfile.NamedTemporaryFile(suffix=".bam", delete=False) as f:
            f.write(b"This is not a BAM file")
            fake_bam = Path(f.name)

        try:
            with pytest.raises((Exception, MappingError)):
                calculate_mapping_rate(fake_bam)
        finally:
            fake_bam.unlink()  # Clean up

    @patch("src.ssiamb.mapping.check_external_tools")
    def test_missing_tools_error_handling(self, mock_check):
        """Test handling when mapping tools are not available."""
        mock_check.return_value = {
            "minimap2": {"available": False, "version": "not found in PATH"},
            "bwa-mem2": {"available": False, "version": "not found in PATH"},
            "samtools": {"available": False, "version": "not found in PATH"},
        }

        # Create a real temporary file to avoid the FileNotFoundError
        with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False) as f:
            f.write(b">test\nACGT\n")
            fake_ref = Path(f.name)

        try:
            # Should raise appropriate errors when tools are missing
            with pytest.raises((ExternalToolError, MappingError)):
                ensure_indexes_self(fake_ref, Mapper.MINIMAP2)
        finally:
            fake_ref.unlink()  # Clean up


class TestMappingPerformance:
    """Performance tests with real data."""

    @pytest.fixture
    def small_reference_path(self):
        """Use the real reference for performance testing."""
        path = Path("fixtures/self_mode/assembly.fna")
        if not path.exists():
            pytest.skip("Real reference not available for performance test")
        return path

    def test_indexing_performance_comparison(self, small_reference_path):
        """Compare indexing performance between mappers."""
        index_times = {}

        for mapper in [Mapper.MINIMAP2, Mapper.BWA_MEM2]:
            with tempfile.TemporaryDirectory() as tmpdir:
                temp_ref = Path(tmpdir) / f"ref_{mapper.value}.fna"
                temp_ref.write_bytes(small_reference_path.read_bytes())

                start_time = time.time()
                ensure_indexes_self(temp_ref, mapper)
                index_time = time.time() - start_time

                index_times[mapper.value] = index_time

                # Sanity check - indexing should complete in reasonable time
                assert (
                    index_time < 120
                ), f"{mapper.value} indexing took too long: {index_time:.2f}s"

        print("\\nIndexing performance:")
        for mapper, index_time in index_times.items():
            print(f"  {mapper}: {index_time:.2f} seconds")

    def test_mapping_performance_scaling(self, small_reference_path):
        """Test mapping performance with different thread counts."""
        reads_r1 = Path("fixtures/self_mode/reads_R1.fastq.gz")
        reads_r2 = Path("fixtures/self_mode/reads_R2.fastq.gz")

        if not reads_r1.exists() or not reads_r2.exists():
            pytest.skip("FASTQ files not available for performance test")

        thread_counts = [1, 2, 4]
        performance_results = {}

        for threads in thread_counts:
            with tempfile.TemporaryDirectory() as tmpdir:
                tmpdir = Path(tmpdir)
                temp_ref = tmpdir / "reference.fna"
                temp_ref.write_bytes(small_reference_path.read_bytes())

                ensure_indexes_self(temp_ref, Mapper.MINIMAP2)

                output_bam = tmpdir / "mapped.bam"

                start_time = time.time()
                map_fastqs(
                    mapper=Mapper.MINIMAP2,
                    fasta_path=temp_ref,
                    r1_path=reads_r1,
                    r2_path=reads_r2,
                    sample_name="test_sample",
                    threads=threads,
                    output_path=output_bam,
                )
                map_time = time.time() - start_time

                performance_results[threads] = map_time

                # Should complete in reasonable time regardless of thread count
                assert (
                    map_time < 180
                ), f"Mapping with {threads} threads took too long: {map_time:.2f}s"

        print("\\nMapping performance by thread count:")
        for threads, map_time in performance_results.items():
            print(f"  {threads} threads: {map_time:.2f} seconds")

        # More threads should generally be faster (or at least not much slower)
        # Allow some variance due to overhead and small dataset size
        max_time = max(performance_results.values())
        min_time = min(performance_results.values())
        speedup_ratio = max_time / min_time
        assert (
            speedup_ratio < 3.0
        ), f"Performance variance too high: {speedup_ratio:.2f}x"
