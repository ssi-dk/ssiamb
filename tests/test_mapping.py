"""
Tests for mapping.py functions.

Tests mapping workflow including tool detection, index management,
mapping execution, and mapping rate calculation.
"""

import pytest
import tempfile
import subprocess
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


class TestExternalToolDetection:
    """Test external tool availability detection."""

    @pytest.mark.parametrize(
        "availability_scenario,expected_tools",
        [
            # All tools available
            (
                {
                    "minimap2": "/usr/bin/minimap2",
                    "bwa-mem2": "/usr/bin/bwa-mem2",
                    "samtools": "/usr/bin/samtools",
                },
                {"minimap2": True, "bwa-mem2": True, "samtools": True},
            ),
            # Some tools missing (only minimap2 and samtools available)
            (
                {
                    "minimap2": "/usr/bin/minimap2",
                    "bwa-mem2": None,
                    "samtools": "/usr/bin/samtools",
                },
                {"minimap2": True, "bwa-mem2": False, "samtools": True},
            ),
            # All tools missing
            (
                {"minimap2": None, "bwa-mem2": None, "samtools": None},
                {"minimap2": False, "bwa-mem2": False, "samtools": False},
            ),
            # Only samtools available
            (
                {"minimap2": None, "bwa-mem2": None, "samtools": "/usr/bin/samtools"},
                {"minimap2": False, "bwa-mem2": False, "samtools": True},
            ),
        ],
    )
    @patch("src.ssiamb.mapping.shutil.which")
    def test_check_external_tools(
        self, mock_which, availability_scenario, expected_tools
    ):
        """Test tool detection across different availability scenarios."""

        def mock_which_side_effect(tool):
            return availability_scenario.get(tool)

        mock_which.side_effect = mock_which_side_effect

        tools = check_external_tools()

        for tool, expected_available in expected_tools.items():
            assert tools[tool]["available"] is expected_available
        assert len(tools) == 3


class TestIndexFiles:
    """Test index file path generation and detection."""

    def test_get_index_files_minimap2(self):
        """Test index file paths for minimap2."""
        fasta_path = Path("/tmp/reference.fasta")

        index_files = get_index_files(fasta_path, Mapper.MINIMAP2)

        assert len(index_files) == 1
        assert index_files[0] == Path("/tmp/reference.mmi")

    def test_get_index_files_bwa_mem2(self):
        """Test index file paths for bwa-mem2."""
        fasta_path = Path("/tmp/reference.fasta")

        index_files = get_index_files(fasta_path, Mapper.BWA_MEM2)

        expected_files = [
            Path("/tmp/reference.fasta.0123"),
            Path("/tmp/reference.fasta.amb"),
            Path("/tmp/reference.fasta.ann"),
            Path("/tmp/reference.fasta.pac"),
            Path("/tmp/reference.fasta.bwt.2bit.64"),
        ]

        assert len(index_files) == 5
        assert set(index_files) == set(expected_files)

    def test_get_index_files_fasta_extension_variants(self):
        """Test index file generation with different FASTA extensions."""
        for ext in [".fasta", ".fa", ".fna"]:
            fasta_path = Path(f"/tmp/reference{ext}")

            index_files = get_index_files(fasta_path, Mapper.MINIMAP2)

            assert index_files[0] == Path("/tmp/reference.mmi")

    # def test_get_index_files_unsupported_mapper(self):
    #     """Test error handling for unsupported mapper."""
    #     # This test is disabled due to type checking constraints
    #     # The function handles this case, but we can't easily test it
    #     # without bypassing type checking
    #     pass

    def test_indexes_exist_all_present(self):
        """Test index existence check when all files are present."""
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = Path(tmpdir) / "reference.fasta"
            fasta_path.touch()

            # Create minimap2 index
            index_file = Path(tmpdir) / "reference.mmi"
            index_file.touch()

            assert indexes_exist(fasta_path, Mapper.MINIMAP2) is True

    def test_indexes_exist_some_missing(self):
        """Test index existence check when some files are missing."""
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = Path(tmpdir) / "reference.fasta"
            fasta_path.touch()

            # Create only some bwa-mem2 index files
            (Path(tmpdir) / "reference.0123").touch()
            (Path(tmpdir) / "reference.amb").touch()
            # Missing .ann, .pac, .bwt.2bit.64

            assert indexes_exist(fasta_path, Mapper.BWA_MEM2) is False

    def test_indexes_exist_none_present(self):
        """Test index existence check when no files are present."""
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = Path(tmpdir) / "reference.fasta"
            fasta_path.touch()

            assert indexes_exist(fasta_path, Mapper.MINIMAP2) is False
            assert indexes_exist(fasta_path, Mapper.BWA_MEM2) is False


class TestIndexBuilding:
    """Test index building functionality."""

    @patch("pathlib.Path.exists")
    @patch("src.ssiamb.mapping.indexes_exist")
    def test_ensure_indexes_self_already_exist(
        self, mock_indexes_exist, mock_path_exists
    ):
        """Test index building when indexes already exist."""
        mock_indexes_exist.return_value = True
        mock_path_exists.return_value = True

        fasta_path = Path("/tmp/reference.fasta")

        # Should not raise any errors
        ensure_indexes_self(fasta_path, Mapper.MINIMAP2)

        mock_indexes_exist.assert_called_once_with(fasta_path, Mapper.MINIMAP2)

    @patch("pathlib.Path.exists")
    @patch("src.ssiamb.mapping.check_external_tools")
    @patch("src.ssiamb.mapping.indexes_exist")
    def test_ensure_indexes_self_tool_missing(
        self, mock_indexes_exist, mock_check_tools, mock_path_exists
    ):
        """Test index building when required tool is missing."""
        mock_indexes_exist.return_value = False
        mock_check_tools.return_value = {
            "minimap2": {"available": False, "version": "not found in PATH"},
            "bwa-mem2": {"available": False, "version": "not found in PATH"},
            "samtools": {"available": True, "version": "1.0"},
        }
        mock_path_exists.return_value = True

        fasta_path = Path("/tmp/reference.fasta")

        with pytest.raises(ExternalToolError, match="minimap2 not available"):
            ensure_indexes_self(fasta_path, Mapper.MINIMAP2)

    @patch("pathlib.Path.exists")
    @patch("src.ssiamb.mapping._build_minimap2_index")
    @patch("src.ssiamb.mapping.check_external_tools")
    @patch("src.ssiamb.mapping.indexes_exist")
    def test_ensure_indexes_self_minimap2_build(
        self, mock_indexes_exist, mock_check_tools, mock_build, mock_path_exists
    ):
        """Test minimap2 index building."""
        mock_indexes_exist.return_value = False
        mock_check_tools.return_value = {
            "minimap2": {"available": True, "version": "2.30"},
            "bwa-mem2": {"available": True, "version": "2.0"},
            "samtools": {"available": True, "version": "1.0"},
        }
        mock_path_exists.return_value = True

        fasta_path = Path("/tmp/reference.fasta")

        ensure_indexes_self(fasta_path, Mapper.MINIMAP2)

        mock_build.assert_called_once_with(fasta_path)

    @patch("pathlib.Path.exists")
    @patch("src.ssiamb.mapping._build_bwa_mem2_index")
    @patch("src.ssiamb.mapping.check_external_tools")
    @patch("src.ssiamb.mapping.indexes_exist")
    def test_ensure_indexes_self_bwa_mem2_build(
        self, mock_indexes_exist, mock_check_tools, mock_build, mock_path_exists
    ):
        """Test bwa-mem2 index building."""
        mock_indexes_exist.return_value = False
        mock_check_tools.return_value = {
            "minimap2": {"available": True, "version": "2.30"},
            "bwa-mem2": {"available": True, "version": "2.0"},
            "samtools": {"available": True, "version": "1.0"},
        }
        mock_path_exists.return_value = True

        fasta_path = Path("/tmp/reference.fasta")

        ensure_indexes_self(fasta_path, Mapper.BWA_MEM2)

        mock_build.assert_called_once_with(fasta_path)


class TestMappingExecution:
    """Test mapping execution functionality."""

    def test_map_fastqs_input_validation(self):
        """Test input file validation for mapping."""
        fasta_path = Path("/non/existent/reference.fasta")
        r1_path = Path("/non/existent/R1.fastq")
        r2_path = Path("/non/existent/R2.fastq")

        with pytest.raises(FileNotFoundError):
            map_fastqs(Mapper.MINIMAP2, fasta_path, r1_path, r2_path, "test_sample")

    @patch("src.ssiamb.mapping.check_external_tools")
    def test_map_fastqs_tool_missing(self, mock_check_tools):
        """Test mapping when required tools are missing."""
        mock_check_tools.return_value = {
            "minimap2": {"available": False, "version": "not found in PATH"},
            "bwa-mem2": {"available": False, "version": "not found in PATH"},
            "samtools": {"available": False, "version": "not found in PATH"},
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            # Create input files
            fasta_path = Path(tmpdir) / "reference.fasta"
            r1_path = Path(tmpdir) / "R1.fastq"
            r2_path = Path(tmpdir) / "R2.fastq"

            fasta_path.write_text(">ref\nACGT\n")
            r1_path.write_text("@read1\nACGT\n+\nIIII\n")
            r2_path.write_text("@read2\nACGT\n+\nIIII\n")

            with pytest.raises(ExternalToolError, match="Required tools not available"):
                map_fastqs(Mapper.MINIMAP2, fasta_path, r1_path, r2_path, "test_sample")

    @patch("src.ssiamb.mapping.index_bam")
    @patch("src.ssiamb.mapping._map_with_minimap2")
    @patch("src.ssiamb.mapping.ensure_indexes_self")
    @patch("src.ssiamb.mapping.check_external_tools")
    def test_map_fastqs_minimap2_success(
        self, mock_check_tools, mock_ensure_indexes, mock_map, mock_index
    ):
        """Test successful mapping with minimap2."""
        mock_check_tools.return_value = {
            "minimap2": {"available": True, "version": "2.30"},
            "bwa-mem2": {"available": True, "version": "2.0"},
            "samtools": {"available": True, "version": "1.0"},
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            # Create input files
            fasta_path = Path(tmpdir) / "reference.fasta"
            r1_path = Path(tmpdir) / "R1.fastq"
            r2_path = Path(tmpdir) / "R2.fastq"
            expected_bam = Path("test_sample.sorted.bam")

            fasta_path.write_text(">ref\nACGT\n")
            r1_path.write_text("@read1\nACGT\n+\nIIII\n")
            r2_path.write_text("@read2\nACGT\n+\nIIII\n")

            # Mock the mapping function to create the expected BAM file
            def create_bam_file(*args, **kwargs):
                expected_bam.touch()

            mock_map.side_effect = create_bam_file

            result = map_fastqs(
                Mapper.MINIMAP2, fasta_path, r1_path, r2_path, "test_sample"
            )

            assert result == expected_bam
            mock_ensure_indexes.assert_called_once_with(fasta_path, Mapper.MINIMAP2)
            mock_map.assert_called_once()
            mock_index.assert_called_once_with(expected_bam)

    @patch("src.ssiamb.mapping.index_bam")
    @patch("src.ssiamb.mapping._map_with_bwa_mem2")
    @patch("src.ssiamb.mapping.ensure_indexes_self")
    @patch("src.ssiamb.mapping.check_external_tools")
    def test_map_fastqs_bwa_mem2_success(
        self, mock_check_tools, mock_ensure_indexes, mock_map, mock_index
    ):
        """Test successful mapping with bwa-mem2."""
        mock_check_tools.return_value = {
            "minimap2": {"available": True, "version": "2.30"},
            "bwa-mem2": {"available": True, "version": "2.0"},
            "samtools": {"available": True, "version": "1.0"},
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            # Create input files
            fasta_path = Path(tmpdir) / "reference.fasta"
            r1_path = Path(tmpdir) / "R1.fastq"
            r2_path = Path(tmpdir) / "R2.fastq"
            expected_bam = Path("test_sample.sorted.bam")

            fasta_path.write_text(">ref\nACGT\n")
            r1_path.write_text("@read1\nACGT\n+\nIIII\n")
            r2_path.write_text("@read2\nACGT\n+\nIIII\n")

            # Mock the mapping function to create the expected BAM file
            def create_bam_file(*args, **kwargs):
                expected_bam.touch()

            mock_map.side_effect = create_bam_file

            result = map_fastqs(
                Mapper.BWA_MEM2, fasta_path, r1_path, r2_path, "test_sample"
            )

            assert result == expected_bam
            mock_ensure_indexes.assert_called_once_with(fasta_path, Mapper.BWA_MEM2)
            mock_map.assert_called_once()
            mock_index.assert_called_once_with(expected_bam)

    def test_map_fastqs_custom_output_path(self):
        """Test mapping with custom output path."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create input files
            fasta_path = Path(tmpdir) / "reference.fasta"
            r1_path = Path(tmpdir) / "R1.fastq"
            r2_path = Path(tmpdir) / "R2.fastq"
            custom_output = Path(tmpdir) / "custom_output.bam"

            fasta_path.write_text(">ref\nACGT\n")
            r1_path.write_text("@read1\nACGT\n+\nIIII\n")
            r2_path.write_text("@read2\nACGT\n+\nIIII\n")

            with (
                patch("src.ssiamb.mapping.check_external_tools") as mock_check_tools,
                patch("src.ssiamb.mapping.ensure_indexes_self"),
                patch("src.ssiamb.mapping._map_with_minimap2") as mock_map,
                patch("src.ssiamb.mapping.index_bam") as mock_index,
            ):

                mock_check_tools.return_value = {
                    "minimap2": {"available": True, "version": "2.30"},
                    "bwa-mem2": {"available": True, "version": "2.0"},
                    "samtools": {"available": True, "version": "1.0"},
                }

                # Mock the mapping function to create the expected BAM file
                def create_bam_file(*args, **kwargs):
                    custom_output.touch()

                mock_map.side_effect = create_bam_file

                result = map_fastqs(
                    Mapper.MINIMAP2,
                    fasta_path,
                    r1_path,
                    r2_path,
                    "test_sample",
                    output_path=custom_output,
                )

                assert result == custom_output
                mock_index.assert_called_once_with(custom_output)


class TestMappingRate:
    """Test mapping rate calculation functionality."""

    @patch("src.ssiamb.mapping.subprocess.run")
    @patch("src.ssiamb.mapping.shutil.which")
    def test_calculate_mapping_rate_success(self, mock_which, mock_run):
        """Test successful mapping rate calculation."""
        mock_which.return_value = "/usr/bin/samtools"

        # Mock samtools stats output with correct format
        mock_run.return_value.returncode = 0
        mock_run.return_value.stdout = """
SN	raw total sequences:	1000
SN	reads mapped:	850
SN	reads unmapped:	150
"""

        with tempfile.NamedTemporaryFile(suffix=".bam") as tmp_bam:
            bam_path = Path(tmp_bam.name)

            mapping_rate = calculate_mapping_rate(bam_path)

            assert mapping_rate == 0.85  # 850/1000
            mock_run.assert_called_once()

    @patch("src.ssiamb.mapping.subprocess.run")
    @patch("src.ssiamb.mapping.shutil.which")
    def test_calculate_mapping_rate_zero_reads(self, mock_which, mock_run):
        """Test mapping rate calculation with zero reads."""
        mock_which.return_value = "/usr/bin/samtools"

        # Mock samtools stats output with zero reads
        mock_run.return_value.returncode = 0
        mock_run.return_value.stdout = """
SN	raw total sequences:	0
SN	reads mapped:	0
SN	reads unmapped:	0
"""

        with tempfile.NamedTemporaryFile(suffix=".bam") as tmp_bam:
            bam_path = Path(tmp_bam.name)

            mapping_rate = calculate_mapping_rate(bam_path)

            assert mapping_rate == 0.0

    @patch("src.ssiamb.mapping.subprocess.run")
    @patch("src.ssiamb.mapping.shutil.which")
    def test_calculate_mapping_rate_samtools_error(self, mock_which, mock_run):
        """Test mapping rate calculation when samtools fails."""
        mock_which.return_value = "/usr/bin/samtools"
        # Mock subprocess.CalledProcessError for non-zero return code
        mock_run.side_effect = subprocess.CalledProcessError(
            1, ["samtools", "stats"], stderr="samtools: error reading file"
        )

        with tempfile.NamedTemporaryFile(suffix=".bam") as tmp_bam:
            bam_path = Path(tmp_bam.name)

            with pytest.raises(MappingError, match="samtools stats failed"):
                calculate_mapping_rate(bam_path)


class TestMappingErrorHandling:
    """Test error handling in mapping functions."""

    @pytest.mark.parametrize(
        "error_class,message",
        [
            (MappingError, "Test mapping error"),
            (ExternalToolError, "Tool not found"),
            (MappingError, "Failed to build index"),
            (ExternalToolError, "samtools not available"),
        ],
    )
    def test_error_creation(self, error_class, message):
        """Test error exception creation for various error types."""
        error = error_class(message)
        assert str(error) == message
        assert isinstance(error, Exception)

    @pytest.mark.parametrize(
        "scenario,mock_setup,expected_error,expected_message",
        [
            # Tool not available scenarios
            (
                "samtools_missing",
                {
                    "which_return": None,
                    "setup_patches": ["src.ssiamb.mapping.shutil.which"],
                },
                ExternalToolError,
                "samtools not found",
            ),
            # File not found scenarios
            (
                "bam_file_missing",
                {
                    "which_return": "/usr/bin/samtools",
                    "bam_path": "/non/existent/file.bam",
                    "setup_patches": ["src.ssiamb.mapping.shutil.which"],
                },
                MappingError,
                "BAM file not found",
            ),
        ],
    )
    @patch("src.ssiamb.mapping.subprocess.run")
    @patch("src.ssiamb.mapping.shutil.which")
    def test_function_error_scenarios(
        self,
        mock_which,
        mock_run,
        scenario,
        mock_setup,
        expected_error,
        expected_message,
    ):
        """Test error scenarios across different mapping functions."""
        # Setup mocks based on scenario
        mock_which.return_value = mock_setup.get("which_return")

        # Test the appropriate function based on scenario
        bam_path = Path(mock_setup.get("bam_path", "/tmp/test.bam"))

        with pytest.raises(expected_error, match=expected_message):
            calculate_mapping_rate(bam_path)

    @patch("src.ssiamb.mapping.subprocess.run")
    @patch("src.ssiamb.mapping.shutil.which")
    def test_subprocess_failure_error(self, mock_which, mock_run):
        """Test subprocess failure scenario separately with proper file setup."""
        mock_which.return_value = "/usr/bin/samtools"
        mock_run.side_effect = subprocess.CalledProcessError(1, "samtools")

        # Create a temporary BAM file for this test
        with tempfile.NamedTemporaryFile(suffix=".bam", delete=False) as tmp_bam:
            bam_path = Path(tmp_bam.name)
            tmp_bam.write(b"dummy bam content")

        try:
            with pytest.raises(MappingError, match="samtools stats failed"):
                calculate_mapping_rate(bam_path)
        finally:
            bam_path.unlink()  # Clean up the temporary file

    @patch("src.ssiamb.mapping.subprocess.run")
    def test_index_building_subprocess_error(self, mock_run):
        """Test index building when subprocess fails."""
        mock_run.side_effect = subprocess.CalledProcessError(1, "minimap2")

        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = Path(tmpdir) / "reference.fasta"
            fasta_path.write_text(">ref\nACGT\n")

            with (
                patch("src.ssiamb.mapping.check_external_tools") as mock_check_tools,
                patch("src.ssiamb.mapping.indexes_exist") as mock_indexes_exist,
            ):

                mock_check_tools.return_value = {
                    "minimap2": {"available": True, "version": "2.30"},
                    "bwa-mem2": {"available": True, "version": "2.0"},
                    "samtools": {"available": True, "version": "1.0"},
                }
                mock_indexes_exist.return_value = False

                with pytest.raises(MappingError, match="Failed to build"):
                    ensure_indexes_self(fasta_path, Mapper.MINIMAP2)
