"""
Tests for provenance tracking system (milestone 17).

Tests provenance record creation, JSON output formatting, MD5 calculation,
and integration with different analysis modes.
"""

import json
import tempfile
import hashlib
from pathlib import Path
from datetime import datetime
from unittest.mock import patch

from src.ssiamb.provenance import (
    ProvenanceRecord,
    ProvenanceInput,
    ProvenanceMappingStats,
    ProvenanceCounts,
    create_provenance_record,
    write_provenance_json,
    calculate_md5,
    get_conda_env,
    get_tool_version,
)
from src.ssiamb.models import Mode
from src.ssiamb.qc import QCWarning


class TestProvenanceDataClasses:
    """Test provenance data structure classes."""

    def test_provenance_input_creation(self):
        """Test ProvenanceInput creation and attributes."""
        input_file = ProvenanceInput(path="/path/to/file.fastq", md5="abc123")
        assert input_file.path == "/path/to/file.fastq"
        assert input_file.md5 == "abc123"

    def test_provenance_mapping_stats(self):
        """Test ProvenanceMappingStats creation."""
        stats = ProvenanceMappingStats(
            total_reads=1000000,
            mapped_reads=950000,
            map_rate=0.95,
            mean_depth=45.2,
            breadth_1x=0.98,
            breadth_10x=0.85,
        )
        assert stats.map_rate == 0.95
        assert stats.breadth_10x == 0.85

    def test_provenance_counts(self):
        """Test ProvenanceCounts creation."""
        counts = ProvenanceCounts(
            ambiguous_snv_count=42,
            ambiguous_indel_count=5,
            ambiguous_del_count=3,
            callable_bases=2500000,
            genome_length=2750000,
        )
        assert counts.ambiguous_snv_count == 42
        assert counts.callable_bases == 2500000


class TestMD5Calculation:
    """Test MD5 hash calculation for files."""

    def test_calculate_md5_success(self):
        """Test successful MD5 calculation."""
        test_content = b"Hello, World!"
        expected_md5 = hashlib.md5(test_content).hexdigest()

        with tempfile.NamedTemporaryFile(delete=False) as tmp_file:
            tmp_file.write(test_content)
            tmp_file.flush()

            result_md5 = calculate_md5(Path(tmp_file.name))
            assert result_md5 == expected_md5

        # Clean up
        Path(tmp_file.name).unlink()

    def test_calculate_md5_missing_file(self):
        """Test MD5 calculation for missing file."""
        non_existent_file = Path("/non/existent/file.txt")
        result = calculate_md5(non_existent_file)
        assert result is None

    def test_calculate_md5_permission_error(self):
        """Test MD5 calculation with permission error."""
        with patch("builtins.open", side_effect=PermissionError):
            result = calculate_md5(Path("/some/file.txt"))
            assert result is None


class TestEnvironmentInfo:
    """Test environment information gathering."""

    @patch.dict("os.environ", {"CONDA_DEFAULT_ENV": "ssiamb-env"})
    def test_get_conda_env_with_env(self):
        """Test getting conda environment name."""
        env_name = get_conda_env()
        assert env_name == "ssiamb-env"

    @patch.dict("os.environ", {"CONDA_DEFAULT_ENV": "base"})
    def test_get_conda_env_base(self):
        """Test that base environment returns None."""
        env_name = get_conda_env()
        assert env_name is None

    @patch.dict("os.environ", {}, clear=True)
    def test_get_conda_env_no_env(self):
        """Test getting conda environment when not set."""
        env_name = get_conda_env()
        assert env_name is None

    def test_get_tool_version(self):
        """Test getting tool version."""
        version = get_tool_version()
        # Should return either a version string or "unknown"
        assert isinstance(version, str)
        assert len(version) > 0


class TestProvenanceRecordCreation:
    """Test creation of complete provenance records."""

    def test_create_basic_provenance_record(self):
        """Test creation of basic provenance record."""
        start_time = datetime(2025, 9, 17, 10, 0, 0)
        end_time = datetime(2025, 9, 17, 10, 30, 0)

        counts = ProvenanceCounts(
            ambiguous_snv_count=42, callable_bases=2500000, genome_length=2750000
        )

        record = create_provenance_record(
            sample="test_sample",
            mode=Mode.SELF,
            started_at=start_time,
            finished_at=end_time,
            threads=4,
            mapper="minimap2",
            caller="bbtools",
            dp_min=10,
            maf_min=0.1,
            dp_cap=100,
            denom_policy="exclude_dups",
            depth_tool="mosdepth",
            emit_vcf=True,
            emit_bed=False,
            emit_matrix=True,
            emit_per_contig=False,
            counts=counts,
        )

        # Check basic fields
        assert record.sample == "test_sample"
        assert record.mode == "self"
        assert record.mapper == "minimap2"
        assert record.caller == "bbtools"
        assert record.threads == 4
        assert record.runtime_sec == 1800.0  # 30 minutes

        # Check thresholds
        assert record.thresholds["dp_min"] == 10
        assert record.thresholds["maf_min"] == 0.1
        assert record.thresholds["dp_cap"] == 100

        # Check grid cell
        assert record.grid_cell_used["depth"] == 10
        assert record.grid_cell_used["maf_bin"] == 10  # floor(100 * 0.1)

        # Check extras emitted
        assert record.extras_emitted["vcf"] is True
        assert record.extras_emitted["bed"] is False
        assert record.extras_emitted["matrix"] is True
        assert record.extras_emitted["per_contig"] is False

    def test_create_provenance_with_inputs(self):
        """Test provenance record creation with input files."""
        start_time = datetime.now()
        end_time = datetime.now()

        # Create temporary files for testing
        with (
            tempfile.NamedTemporaryFile(suffix=".fastq") as r1_file,
            tempfile.NamedTemporaryFile(suffix=".fastq") as r2_file,
            tempfile.NamedTemporaryFile(suffix=".fasta") as assembly_file,
        ):

            # Write some content to calculate MD5
            r1_file.write(b"@read1\nACGT\n+\nIIII\n")
            r2_file.write(b"@read2\nTGCA\n+\nIIII\n")
            assembly_file.write(b">contig1\nACGTACGT\n")
            r1_file.flush()
            r2_file.flush()
            assembly_file.flush()

            record = create_provenance_record(
                sample="test_sample",
                mode=Mode.SELF,
                started_at=start_time,
                finished_at=end_time,
                threads=4,
                mapper="minimap2",
                caller="bbtools",
                dp_min=10,
                maf_min=0.1,
                dp_cap=100,
                denom_policy="exclude_dups",
                depth_tool="mosdepth",
                emit_vcf=False,
                emit_bed=False,
                emit_matrix=False,
                emit_per_contig=False,
                r1=Path(r1_file.name),
                r2=Path(r2_file.name),
                assembly=Path(assembly_file.name),
            )

        # Check that inputs are recorded with MD5s
        assert "r1" in record.inputs
        assert "r2" in record.inputs
        assert "assembly" in record.inputs

        assert record.inputs["r1"]["path"] == r1_file.name
        assert record.inputs["r2"]["path"] == r2_file.name
        assert record.inputs["assembly"]["path"] == assembly_file.name

        # MD5s should be calculated
        assert record.inputs["r1"]["md5"] is not None
        assert record.inputs["r2"]["md5"] is not None
        assert record.inputs["assembly"]["md5"] is not None
        assert len(record.inputs["r1"]["md5"]) == 32  # MD5 hex length

    def test_create_provenance_ref_mode(self):
        """Test provenance record for ref mode."""
        start_time = datetime.now()
        end_time = datetime.now()

        mapping_stats = ProvenanceMappingStats(map_rate=0.92, breadth_10x=0.85)

        record = create_provenance_record(
            sample="test_sample",
            mode=Mode.REF,
            started_at=start_time,
            finished_at=end_time,
            threads=4,
            mapper="minimap2",
            caller="bbtools",
            dp_min=10,
            maf_min=0.1,
            dp_cap=100,
            denom_policy="exclude_dups",
            depth_tool="mosdepth",
            emit_vcf=False,
            emit_bed=False,
            emit_matrix=False,
            emit_per_contig=False,
            species="Listeria monocytogenes",
            mapping_stats=mapping_stats,
            reference_species="Listeria monocytogenes",
        )

        # Check ref mode specific fields
        assert record.mode == "ref"
        assert record.reference_info["species_requested"] == "Listeria monocytogenes"
        assert record.reference_info["species_final"] == "Listeria monocytogenes"
        assert record.mapping_stats is not None
        assert record.mapping_stats["map_rate"] == 0.92
        assert record.mapping_stats["breadth_10x"] == 0.85

    def test_create_provenance_with_warnings(self):
        """Test provenance record creation with QC warnings."""
        start_time = datetime.now()
        end_time = datetime.now()

        warnings = [
            QCWarning("breadth_10x", 0.75, 0.80, "Low breadth: 0.750 < 0.800"),
            QCWarning(
                "callable_bases",
                500_000,
                1_000_000,
                "Low callable bases: 500,000 < 1,000,000",
            ),
        ]

        record = create_provenance_record(
            sample="test_sample",
            mode=Mode.SELF,
            started_at=start_time,
            finished_at=end_time,
            threads=4,
            mapper="minimap2",
            caller="bbtools",
            dp_min=10,
            maf_min=0.1,
            dp_cap=100,
            denom_policy="exclude_dups",
            depth_tool="mosdepth",
            emit_vcf=False,
            emit_bed=False,
            emit_matrix=False,
            emit_per_contig=False,
            warnings=warnings,
        )

        # Check warnings are formatted correctly
        assert record.warnings is not None
        assert len(record.warnings) == 2
        assert "breadth_10x: Low breadth: 0.750 < 0.800" in record.warnings
        assert (
            "callable_bases: Low callable bases: 500,000 < 1,000,000" in record.warnings
        )


class TestProvenanceJSONOutput:
    """Test JSON output generation."""

    def test_write_provenance_json_single_record(self):
        """Test writing single provenance record to JSON."""
        record = ProvenanceRecord(
            tool_version="1.0.0",
            python_version="3.12.0",
            conda_env="test-env",
            started_at="2025-09-17T10:00:00",
            finished_at="2025-09-17T10:30:00",
            runtime_sec=1800.0,
            threads=4,
            sample="test_sample",
            mode="self",
            mapper="minimap2",
            caller="bbtools",
            thresholds={"dp_min": 10, "maf_min": 0.1, "dp_cap": 100},
            denom_policy="exclude_dups",
            depth_tool="mosdepth",
            inputs={"r1": {"path": "/path/to/r1.fastq", "md5": "abc123"}},
            reference_info={},
            contig_filters={"min_contig_len": 500},
            caller_params={"exact_cmdlines": []},
            counts={"ambiguous_snv_count": 42},
            grid_cell_used={"depth": 10, "maf_bin": 10},
            extras_emitted={
                "vcf": False,
                "bed": False,
                "matrix": False,
                "per_contig": False,
            },
            warnings=[],
        )

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".json", delete=False
        ) as tmp_file:
            output_path = Path(tmp_file.name)

        try:
            write_provenance_json([record], output_path)

            # Read back and verify
            with open(output_path) as f:
                data = json.load(f)

            assert isinstance(data, list)
            assert len(data) == 1

            record_data = data[0]
            assert record_data["sample"] == "test_sample"
            assert record_data["mode"] == "self"
            assert record_data["tool_version"] == "1.0.0"
            assert record_data["runtime_sec"] == 1800.0
            assert record_data["thresholds"]["dp_min"] == 10

        finally:
            output_path.unlink()

    def test_write_provenance_json_multiple_records(self):
        """Test writing multiple provenance records to JSON."""
        records = []
        for i in range(3):
            record = ProvenanceRecord(
                tool_version="1.0.0",
                python_version="3.12.0",
                conda_env=None,
                started_at=f"2025-09-17T10:{i:02d}:00",
                finished_at=f"2025-09-17T10:{i+1:02d}:00",
                runtime_sec=60.0,
                threads=4,
                sample=f"sample_{i}",
                mode="self",
                mapper="minimap2",
                caller="bbtools",
                thresholds={"dp_min": 10, "maf_min": 0.1, "dp_cap": 100},
                denom_policy="exclude_dups",
                depth_tool="mosdepth",
                inputs={},
                reference_info={},
                contig_filters={"min_contig_len": 500},
                caller_params={"exact_cmdlines": []},
                counts={"ambiguous_snv_count": i * 10},
                grid_cell_used={"depth": 10, "maf_bin": 10},
                extras_emitted={
                    "vcf": False,
                    "bed": False,
                    "matrix": False,
                    "per_contig": False,
                },
                warnings=[],
            )
            records.append(record)

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".json", delete=False
        ) as tmp_file:
            output_path = Path(tmp_file.name)

        try:
            write_provenance_json(records, output_path)

            # Read back and verify
            with open(output_path) as f:
                data = json.load(f)

            assert isinstance(data, list)
            assert len(data) == 3

            for i, record_data in enumerate(data):
                assert record_data["sample"] == f"sample_{i}"
                assert record_data["counts"]["ambiguous_snv_count"] == i * 10

        finally:
            output_path.unlink()

    def test_json_output_format_validation(self):
        """Test that JSON output follows expected schema."""
        record = ProvenanceRecord(
            tool_version="1.0.0",
            python_version="3.12.0",
            conda_env="test-env",
            started_at="2025-09-17T10:00:00",
            finished_at="2025-09-17T10:30:00",
            runtime_sec=1800.0,
            threads=4,
            sample="test_sample",
            mode="ref",
            mapper="minimap2",
            caller="bbtools",
            thresholds={"dp_min": 10, "maf_min": 0.1, "dp_cap": 100},
            denom_policy="exclude_dups",
            depth_tool="mosdepth",
            inputs={"r1": {"path": "/path/to/r1.fastq", "md5": "abc123"}},
            reference_info={"species_requested": "Listeria monocytogenes"},
            species_selection={
                "bracken_species": "Listeria monocytogenes",
                "bracken_frac": 0.85,
            },
            mapping_stats={"map_rate": 0.92, "breadth_10x": 0.85},
            contig_filters={"min_contig_len": 500},
            caller_params={"exact_cmdlines": ["bcftools call -mv"]},
            counts={"ambiguous_snv_count": 42, "callable_bases": 2500000},
            grid_cell_used={"depth": 10, "maf_bin": 10},
            extras_emitted={
                "vcf": True,
                "bed": False,
                "matrix": True,
                "per_contig": False,
            },
            warnings=["breadth_10x: Low coverage"],
        )

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".json", delete=False
        ) as tmp_file:
            output_path = Path(tmp_file.name)

        try:
            write_provenance_json([record], output_path)

            # Read back and validate structure
            with open(output_path) as f:
                data = json.load(f)

            record_data = data[0]

            # Check all required top-level fields are present
            required_fields = [
                "tool_version",
                "python_version",
                "started_at",
                "finished_at",
                "runtime_sec",
                "threads",
                "sample",
                "mode",
                "mapper",
                "caller",
                "thresholds",
                "denom_policy",
                "depth_tool",
                "inputs",
                "reference_info",
                "grid_cell_used",
                "extras_emitted",
            ]

            for field in required_fields:
                assert field in record_data, f"Missing required field: {field}"

            # Check nested structure
            assert isinstance(record_data["thresholds"], dict)
            assert isinstance(record_data["inputs"], dict)
            assert isinstance(record_data["extras_emitted"], dict)
            assert isinstance(record_data["warnings"], list)

        finally:
            output_path.unlink()


class TestProvenanceIntegration:
    """Test provenance system integration scenarios."""

    def test_complete_provenance_workflow(self):
        """Test complete provenance workflow from creation to JSON output."""
        # Create a realistic provenance record
        start_time = datetime(2025, 9, 17, 10, 0, 0)
        end_time = datetime(2025, 9, 17, 10, 45, 30)  # 45 minutes 30 seconds

        mapping_stats = ProvenanceMappingStats(map_rate=0.89, breadth_10x=0.82)
        counts = ProvenanceCounts(
            ambiguous_snv_count=67,
            ambiguous_indel_count=5,
            ambiguous_del_count=2,
            callable_bases=2750000,
            genome_length=2800000,
        )
        warnings = [QCWarning("breadth_10x", 0.82, 0.85, "Low breadth: 0.820 < 0.850")]

        record = create_provenance_record(
            sample="integration_test",
            mode=Mode.REF,
            started_at=start_time,
            finished_at=end_time,
            threads=8,
            mapper="bwa-mem2",
            caller="bcftools",
            dp_min=15,
            maf_min=0.05,
            dp_cap=150,
            denom_policy="include_dups",
            depth_tool="mosdepth",
            emit_vcf=True,
            emit_bed=True,
            emit_matrix=False,
            emit_per_contig=True,
            species="Escherichia coli",
            mapping_stats=mapping_stats,
            counts=counts,
            warnings=warnings,
            reference_species="Escherichia coli",
        )

        # Write to JSON
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".json", delete=False
        ) as tmp_file:
            output_path = Path(tmp_file.name)

        try:
            write_provenance_json([record], output_path)

            # Read back and validate complete workflow
            with open(output_path) as f:
                data = json.load(f)

            record_data = data[0]

            # Validate timing calculation
            assert record_data["runtime_sec"] == 2730.0  # 45 min 30 sec

            # Validate parameter capture
            assert record_data["sample"] == "integration_test"
            assert record_data["mode"] == "ref"
            assert record_data["threads"] == 8
            assert record_data["thresholds"]["dp_min"] == 15
            assert record_data["thresholds"]["maf_min"] == 0.05

            # Validate results capture
            assert record_data["counts"]["ambiguous_snv_count"] == 67
            assert record_data["mapping_stats"]["map_rate"] == 0.89

            # Validate flags capture
            assert record_data["extras_emitted"]["vcf"] is True
            assert record_data["extras_emitted"]["bed"] is True
            assert record_data["extras_emitted"]["matrix"] is False

            # Validate warnings capture
            assert len(record_data["warnings"]) == 1
            assert "breadth_10x: Low breadth" in record_data["warnings"][0]

        finally:
            output_path.unlink()
