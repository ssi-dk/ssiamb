"""
Quality control policies and warning thresholds for ssiamb.

This module implements QC policies that warn but do not fail runs,
as specified in spec.md §8.
"""

from dataclasses import dataclass
from typing import List, Optional
import logging

from .models import Mode

logger = logging.getLogger(__name__)


@dataclass
class QCWarning:
    """A QC warning with message and metric values."""

    metric: str
    value: float
    threshold: float
    message: str


@dataclass
class QCThresholds:
    """QC warning thresholds as specified in spec.md §8."""

    min_breadth_10x: float = 0.80
    min_callable_bases: int = 1_000_000  # 1e6
    min_mapping_rate_ref: float = 0.70  # Only for ref mode


def check_qc_metrics(
    breadth_10x: float,
    callable_bases: int,
    mapping_rate: Optional[float] = None,
    mode: Optional[Mode] = None,
    thresholds: Optional[QCThresholds] = None,
) -> List[QCWarning]:
    """
    Check QC metrics against warning thresholds.

    Args:
        breadth_10x: Breadth of coverage at 10x depth (0.0-1.0)
        callable_bases: Number of callable bases
        mapping_rate: Mapping rate (0.0-1.0), required for ref mode
        mode: Analysis mode (used to determine if mapping rate applies)
        thresholds: QC thresholds (uses defaults if None)

    Returns:
        List of QC warnings
    """
    if thresholds is None:
        thresholds = QCThresholds()

    warnings = []

    # Check breadth_10x threshold
    if breadth_10x < thresholds.min_breadth_10x:
        warnings.append(
            QCWarning(
                metric="breadth_10x",
                value=breadth_10x,
                threshold=thresholds.min_breadth_10x,
                message=f"Low breadth of coverage at 10x: {breadth_10x:.3f} < {thresholds.min_breadth_10x:.3f}",
            )
        )

    # Check callable bases threshold
    if callable_bases < thresholds.min_callable_bases:
        warnings.append(
            QCWarning(
                metric="callable_bases",
                value=callable_bases,
                threshold=thresholds.min_callable_bases,
                message=f"Low callable bases: {callable_bases:,} < {thresholds.min_callable_bases:,}",
            )
        )

    # Check mapping rate threshold (ref mode only)
    if mode == Mode.REF and mapping_rate is not None:
        if mapping_rate < thresholds.min_mapping_rate_ref:
            warnings.append(
                QCWarning(
                    metric="mapping_rate",
                    value=mapping_rate,
                    threshold=thresholds.min_mapping_rate_ref,
                    message=f"Low mapping rate in ref mode: {mapping_rate:.3f} < {thresholds.min_mapping_rate_ref:.3f}",
                )
            )

    return warnings


def log_qc_warnings(warnings: List[QCWarning]) -> None:
    """
    Log QC warnings to the logger.

    Args:
        warnings: List of QC warnings to log
    """
    if not warnings:
        logger.info("All QC metrics passed warning thresholds")
        return

    logger.warning(f"QC warnings detected ({len(warnings)} issues):")
    for warning in warnings:
        logger.warning(f"  - {warning.message}")


def format_qc_warnings_for_summary(warnings: List[QCWarning]) -> str:
    """
    Format QC warnings for inclusion in summary TSV qc_warnings field.

    Args:
        warnings: List of QC warnings

    Returns:
        Semicolon-separated string of warning codes, or empty string if no warnings
    """
    if not warnings:
        return ""

    # Create concise warning codes
    codes = []
    for warning in warnings:
        if warning.metric == "breadth_10x":
            codes.append(f"low_breadth_{warning.value:.2f}")
        elif warning.metric == "callable_bases":
            codes.append(f"low_callable_{warning.value:.0e}")
        elif warning.metric == "mapping_rate":
            codes.append(f"low_maprate_{warning.value:.2f}")

    return ";".join(codes)
