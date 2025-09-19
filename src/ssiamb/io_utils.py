"""
I/O utilities for ssiamb.

This module provides file handling utilities including sample name validation,
TSV writing with atomic operations, and file hashing.
"""

import csv
import hashlib
import os
import re
import tempfile
from pathlib import Path
from typing import List, Optional, Iterator
import fcntl
from contextlib import contextmanager

from .models import TSVMode, SummaryRow


# Sample name validation regex
SAMPLE_NAME_PATTERN = re.compile(r"^[A-Za-z0-9._-]{1,64}$")


class SampleNameError(ValueError):
    """Raised when sample name is invalid."""

    pass


class TSVWriteError(Exception):
    """Raised when TSV writing fails."""

    pass


def validate_sample_name(sample: str) -> str:
    """
    Validate sample name according to ssiamb rules.

    Args:
        sample: Sample name to validate

    Returns:
        Validated sample name

    Raises:
        SampleNameError: If sample name is invalid
    """
    if not sample:
        raise SampleNameError("Sample name cannot be empty")

    if not SAMPLE_NAME_PATTERN.match(sample):
        raise SampleNameError(
            f"Sample name '{sample}' is invalid. "
            "Must be 1-64 characters containing only letters, numbers, "
            "periods, underscores, and hyphens."
        )

    return sample


def infer_sample_name(
    r1: Path,
    r2: Optional[Path] = None,
    vcf: Optional[Path] = None,
    bam: Optional[Path] = None,
) -> str:
    """
    Infer sample name from input filenames.

    Tries to extract a common prefix from R1/R2 filenames, removing
    common suffixes like _R1, _R2, .fastq, .gz, etc.

    Args:
        r1: Forward reads file
        r2: Reverse reads file (optional)
        vcf: VCF file (optional)
        bam: BAM file (optional)

    Returns:
        Inferred sample name

    Raises:
        SampleNameError: If sample name cannot be inferred
    """
    # Start with R1 filename (remove .gz first if present)
    name = r1.name
    if name.endswith(".gz"):
        name = name[:-3]

    # Remove .fastq extension
    if name.endswith(".fastq"):
        name = name[:-6]
    elif name.endswith(".fq"):
        name = name[:-3]

    # Remove R1 suffix
    for suffix in ["_R1", "_1", ".R1", ".1"]:
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break

    # If R2 provided, ensure it has a compatible name
    if r2:
        r2_name = r2.name
        if r2_name.endswith(".gz"):
            r2_name = r2_name[:-3]

        # Remove .fastq extension
        if r2_name.endswith(".fastq"):
            r2_name = r2_name[:-6]
        elif r2_name.endswith(".fq"):
            r2_name = r2_name[:-3]

        # Remove R2 suffix
        for suffix in ["_R2", "_2", ".R2", ".2"]:
            if r2_name.endswith(suffix):
                r2_name = r2_name[: -len(suffix)]
                break

        # Check if R1 and R2 have the same base name
        if name != r2_name:
            raise SampleNameError(
                f"Cannot infer sample name: R1 stem '{name}' and "
                f"R2 stem '{r2_name}' do not match"
            )

    # Validate the inferred name
    if not name:
        raise SampleNameError("Cannot infer sample name: filename too short")

    try:
        validate_sample_name(name)
        return name
    except SampleNameError as e:
        raise SampleNameError(f"Cannot infer valid sample name from '{r1}': {e}")


@contextmanager
def file_lock(file_path: Path) -> Iterator[None]:
    """
    Context manager for file locking (best-effort cross-platform).

    Args:
        file_path: Path to lock
    """
    lock_file = file_path.with_suffix(file_path.suffix + ".lock")

    try:
        with open(lock_file, "w") as f:
            try:
                fcntl.flock(f.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
                yield
            except (OSError, IOError):
                # Lock failed - either file is locked or system doesn't support locks
                # For now, just proceed (best effort)
                yield
    except (OSError, IOError):
        # Can't create lock file - just proceed
        yield
    finally:
        try:
            lock_file.unlink(missing_ok=True)
        except (OSError, IOError):
            pass


def compute_md5(file_path: Path) -> str:
    """
    Compute MD5 hash of a file.

    Args:
        file_path: Path to file

    Returns:
        MD5 hash as hex string
    """
    hash_md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def write_tsv_summary(
    output_path: Path,
    rows: List[SummaryRow],
    mode: TSVMode = TSVMode.OVERWRITE,
) -> None:
    """
    Write summary rows to TSV file with atomic operations.

    Args:
        output_path: Output TSV file path
        rows: List of SummaryRow objects to write
        mode: Write mode (overwrite, append, or fail if exists)

    Raises:
        TSVWriteError: If writing fails
        FileExistsError: If file exists and mode is FAIL
    """
    if not rows:
        raise TSVWriteError("No rows to write")

    # Check file existence and mode
    file_exists = output_path.exists()

    if mode == TSVMode.FAIL and file_exists:
        raise FileExistsError(f"Output file already exists: {output_path}")

    # Determine if we need to write header
    write_header = (mode == TSVMode.OVERWRITE) or (
        mode == TSVMode.APPEND and not file_exists
    )

    # For atomic writes, use temporary file
    with file_lock(output_path):
        if mode == TSVMode.APPEND and file_exists:
            # Direct append
            with open(output_path, "a", newline="") as f:
                writer = csv.DictWriter(
                    f, fieldnames=rows[0].to_dict().keys(), delimiter="\t"
                )
                for row in rows:
                    writer.writerow(row.to_dict())
        else:
            # Atomic write using temporary file
            temp_file = None
            try:
                # Create temporary file in same directory
                temp_fd, temp_path = tempfile.mkstemp(
                    suffix=".tmp",
                    prefix=output_path.name + ".",
                    dir=output_path.parent,
                )
                temp_file = Path(temp_path)

                with os.fdopen(temp_fd, "w", newline="") as f:
                    writer = csv.DictWriter(
                        f, fieldnames=rows[0].to_dict().keys(), delimiter="\t"
                    )

                    if write_header:
                        writer.writeheader()

                    for row in rows:
                        writer.writerow(row.to_dict())

                # Atomic rename
                temp_file.replace(output_path)
                temp_file = None  # Successfully moved

            except Exception as e:
                if temp_file and temp_file.exists():
                    temp_file.unlink()
                raise TSVWriteError(f"Failed to write TSV: {e}") from e


def write_tsv_to_stdout(rows: List[SummaryRow]) -> None:
    """
    Write summary rows to stdout in TSV format.

    Args:
        rows: List of SummaryRow objects to write
    """
    import sys

    if not rows:
        return

    writer = csv.DictWriter(
        sys.stdout, fieldnames=rows[0].to_dict().keys(), delimiter="\t"
    )
    writer.writeheader()

    for row in rows:
        writer.writerow(row.to_dict())
