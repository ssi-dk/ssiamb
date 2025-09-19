"""
Tests for QC warning system (milestone 16).

Tests the QC threshold validation, warning generation, and integration
with the summary reporting system.
"""

import logging

from src.ssiamb.qc import (
    QCThresholds,
    QCWarning,
    check_qc_metrics,
    log_qc_warnings,
    format_qc_warnings_for_summary,
)
from src.ssiamb.models import Mode


class TestQCThresholds:
    """Test QC threshold configuration."""

    def test_default_thresholds(self):
        """Test default QC threshold values."""
        thresholds = QCThresholds()
        assert thresholds.min_breadth_10x == 0.80
        assert thresholds.min_callable_bases == 1_000_000
        assert thresholds.min_mapping_rate_ref == 0.70

    def test_custom_thresholds(self):
        """Test custom QC threshold configuration."""
        thresholds = QCThresholds(
            min_breadth_10x=0.90,
            min_callable_bases=2_000_000,
            min_mapping_rate_ref=0.80,
        )
        assert thresholds.min_breadth_10x == 0.90
        assert thresholds.min_callable_bases == 2_000_000
        assert thresholds.min_mapping_rate_ref == 0.80


class TestQCWarning:
    """Test QC warning data structure."""

    def test_qc_warning_creation(self):
        """Test QCWarning creation and attributes."""
        warning = QCWarning(
            metric="breadth_10x",
            value=0.75,
            threshold=0.80,
            message="Low breadth of coverage: 0.750 < 0.800",
        )

        assert warning.metric == "breadth_10x"
        assert "Low breadth of coverage" in warning.message
        assert warning.threshold == 0.80
        assert warning.value == 0.75


class TestQCMetricChecking:
    """Test QC metric checking logic."""

    def test_no_warnings_good_metrics(self):
        """Test that good metrics generate no warnings."""
        warnings = check_qc_metrics(
            breadth_10x=0.90,  # Above threshold
            callable_bases=2_000_000,  # Above threshold
            mapping_rate=0.80,  # Above threshold
            mode=Mode.REF,
        )

        assert len(warnings) == 0

    def test_breadth_warning_generation(self):
        """Test breadth_10x warning generation."""
        warnings = check_qc_metrics(
            breadth_10x=0.75,  # Below threshold (0.80)
            callable_bases=2_000_000,  # Good
            mapping_rate=0.80,  # Good
            mode=Mode.REF,
        )

        assert len(warnings) == 1
        assert warnings[0].metric == "breadth_10x"
        assert "0.750" in warnings[0].message
        assert "0.800" in warnings[0].message
        assert warnings[0].threshold == 0.80
        assert warnings[0].value == 0.75

    def test_callable_bases_warning_generation(self):
        """Test callable_bases warning generation."""
        warnings = check_qc_metrics(
            breadth_10x=0.90,  # Good
            callable_bases=500_000,  # Below threshold (1e6)
            mapping_rate=0.80,  # Good
            mode=Mode.REF,
        )

        assert len(warnings) == 1
        assert warnings[0].metric == "callable_bases"
        assert "500,000" in warnings[0].message
        assert "1,000,000" in warnings[0].message
        assert warnings[0].threshold == 1_000_000
        assert warnings[0].value == 500_000

    def test_mapping_rate_warning_ref_mode(self):
        """Test mapping rate warning in ref mode."""
        warnings = check_qc_metrics(
            breadth_10x=0.90,  # Good
            callable_bases=2_000_000,  # Good
            mapping_rate=0.65,  # Below threshold (0.70) for ref mode
            mode=Mode.REF,
        )

        assert len(warnings) == 1
        assert warnings[0].metric == "mapping_rate"
        assert "0.650" in warnings[0].message
        assert "0.700" in warnings[0].message
        assert warnings[0].threshold == 0.70
        assert warnings[0].value == 0.65

    def test_mapping_rate_no_warning_self_mode(self):
        """Test that mapping rate warnings don't apply to self mode."""
        warnings = check_qc_metrics(
            breadth_10x=0.90,  # Good
            callable_bases=2_000_000,  # Good
            mapping_rate=0.65,  # Would be low for ref mode
            mode=Mode.SELF,  # But self mode doesn't check mapping rate
        )

        assert len(warnings) == 0

    def test_mapping_rate_none_value(self):
        """Test handling of None mapping rate."""
        warnings = check_qc_metrics(
            breadth_10x=0.90,
            callable_bases=2_000_000,
            mapping_rate=None,  # Missing mapping rate
            mode=Mode.REF,
        )

        # Should not crash, but also no mapping rate warning
        breadth_warnings = [w for w in warnings if w.metric == "mapping_rate"]
        assert len(breadth_warnings) == 0

    def test_multiple_warnings(self):
        """Test multiple warnings generated simultaneously."""
        warnings = check_qc_metrics(
            breadth_10x=0.70,  # Below threshold
            callable_bases=800_000,  # Below threshold
            mapping_rate=0.60,  # Below threshold for ref mode
            mode=Mode.REF,
        )

        assert len(warnings) == 3

        # Check that all three warnings are present
        metrics = {w.metric for w in warnings}
        assert metrics == {"breadth_10x", "callable_bases", "mapping_rate"}

    def test_custom_thresholds(self):
        """Test QC checking with custom thresholds."""
        custom_thresholds = QCThresholds(
            min_breadth_10x=0.90,
            min_callable_bases=2_000_000,
            min_mapping_rate_ref=0.80,
        )

        warnings = check_qc_metrics(
            breadth_10x=0.85,  # Would be OK with default, but low with custom
            callable_bases=1_500_000,  # Would be OK with default, but low with custom
            mapping_rate=0.75,  # Would be OK with default, but low with custom
            mode=Mode.REF,
            thresholds=custom_thresholds,
        )

        assert len(warnings) == 3

        # Check that custom thresholds are used
        breadth_warning = next(w for w in warnings if w.metric == "breadth_10x")
        assert breadth_warning.threshold == 0.90


class TestQCWarningFormatting:
    """Test QC warning formatting and logging."""

    def test_format_qc_warnings_for_summary(self):
        """Test formatting warnings for summary output."""
        warnings = [
            QCWarning(
                metric="breadth_10x",
                value=0.75,
                threshold=0.80,
                message="Low breadth: 0.750 < 0.800",
            ),
            QCWarning(
                metric="callable_bases",
                value=500_000,
                threshold=1_000_000,
                message="Low callable bases: 500,000 < 1,000,000",
            ),
        ]

        formatted = format_qc_warnings_for_summary(warnings)
        assert "low_breadth_0.75" in formatted
        assert "low_callable_5e+05" in formatted
        assert ";" in formatted  # Multiple warnings separated by semicolon

    def test_format_empty_warnings(self):
        """Test formatting empty warning list."""
        formatted = format_qc_warnings_for_summary([])
        assert formatted == ""

    def test_format_single_warning(self):
        """Test formatting single warning."""
        warnings = [
            QCWarning(
                metric="breadth_10x",
                value=0.75,
                threshold=0.80,
                message="Low breadth: 0.750 < 0.800",
            )
        ]

        formatted = format_qc_warnings_for_summary(warnings)
        assert formatted == "low_breadth_0.75"
        assert ";" not in formatted  # No semicolon for single warning


class TestQCLogging:
    """Test QC warning logging integration."""

    def test_log_qc_warnings(self, caplog):
        """Test that QC warnings are properly logged."""
        warnings = [
            QCWarning(
                metric="breadth_10x",
                value=0.75,
                threshold=0.80,
                message="Low breadth: 0.750 < 0.800",
            ),
            QCWarning(
                metric="callable_bases",
                value=500_000,
                threshold=1_000_000,
                message="Low callable bases: 500,000 < 1,000,000",
            ),
        ]

        with caplog.at_level(logging.WARNING):
            log_qc_warnings(warnings)

        # Check that warnings were logged
        assert len(caplog.records) >= 1
        assert any(
            "QC warnings detected" in record.message for record in caplog.records
        )
        assert any(
            "breadth_10x" in record.message or "0.750" in record.message
            for record in caplog.records
        )

    def test_log_no_warnings(self, caplog):
        """Test logging when there are no warnings."""
        with caplog.at_level(logging.INFO):
            log_qc_warnings([])

        # Should log info message about passing thresholds
        assert len(caplog.records) >= 1
        assert any(
            "passed warning thresholds" in record.message for record in caplog.records
        )


class TestQCIntegration:
    """Test QC system integration scenarios."""

    def test_full_qc_workflow(self):
        """Test complete QC workflow from metrics to formatted output."""
        # Simulate realistic QC scenario
        warnings = check_qc_metrics(
            breadth_10x=0.75,  # Low
            callable_bases=800_000,  # Low
            mapping_rate=0.85,  # OK
            mode=Mode.REF,
        )

        # Should have 2 warnings
        assert len(warnings) == 2

        # Format for summary
        summary_text = format_qc_warnings_for_summary(warnings)
        assert "low_breadth_0.75" in summary_text
        assert "low_callable_8e+05" in summary_text
        assert "mapping_rate" not in summary_text  # This one was OK

        # Should be properly formatted for downstream use
        assert isinstance(summary_text, str)
        assert len(summary_text) > 0

    def test_edge_case_values(self):
        """Test QC checking with edge case values."""
        # Test exactly at thresholds
        warnings = check_qc_metrics(
            breadth_10x=0.80,  # Exactly at threshold
            callable_bases=1_000_000,  # Exactly at threshold
            mapping_rate=0.70,  # Exactly at threshold
            mode=Mode.REF,
        )

        # Should not generate warnings for values at threshold
        assert len(warnings) == 0

        # Test just below thresholds
        warnings = check_qc_metrics(
            breadth_10x=0.799,  # Just below
            callable_bases=999_999,  # Just below
            mapping_rate=0.699,  # Just below
            mode=Mode.REF,
        )

        # Should generate all warnings
        assert len(warnings) == 3
