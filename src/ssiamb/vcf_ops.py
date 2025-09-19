"""
VCF operations module for ssiamb.

This module implements VCF normalization, atomization, MAF extraction,
variant classification, and grid-based counting for ambiguous sites.
"""

import logging
import subprocess
import shutil
import warnings
from pathlib import Path
from typing import Optional, List, Tuple, Iterator, Set, Any
from dataclasses import dataclass
from enum import Enum

import pysam
import numpy as np

logger = logging.getLogger(__name__)


class VCFOperationError(Exception):
    """Raised when VCF operations fail."""

    pass


class VariantClass(Enum):
    """Variant classification types."""

    SNV = "SNV"
    INS = "INS"
    DEL = "DEL"
    UNKNOWN = "UNKNOWN"


@dataclass
class SiteRecord:
    """Record for a genomic site with variant information."""

    chrom: str
    pos: int
    ref: str
    alt: str
    variant_class: VariantClass
    depth: int
    maf: float
    original_filter: str


@dataclass
class VCFNormalizationResult:
    """Result of VCF normalization operation."""

    normalized_vcf_path: Path
    success: bool
    error_message: Optional[str] = None
    records_processed: int = 0


def check_vcf_tools() -> bool:
    """
    Check if required VCF processing tools are available.

    Returns:
        True if all required tools are available, False otherwise
    """
    required_tools = ["bcftools", "bgzip", "tabix"]
    return all(shutil.which(tool) is not None for tool in required_tools)


def normalize_and_split(
    vcf_in: Path, reference: Path, output_dir: Optional[Path] = None
) -> Path:
    """
    Normalize VCF to reference and decompose into primitive records.

    Uses bcftools norm to:
    1. Normalize variants to reference (-f REF)
    2. Split multiallelic sites (-m -both)
    3. Atomize variants (--atomize)
    4. Compress and index output

    Args:
        vcf_in: Input VCF file path
        reference: Reference FASTA file path
        output_dir: Output directory (defaults to input directory)

    Returns:
        Path to normalized, compressed VCF file

    Raises:
        VCFOperationError: If normalization fails
    """
    if not check_vcf_tools():
        raise VCFOperationError(
            "Required VCF tools not found. Need: bcftools, bgzip, tabix"
        )

    if not vcf_in.exists():
        raise VCFOperationError(f"Input VCF file not found: {vcf_in}")

    if not reference.exists():
        raise VCFOperationError(f"Reference file not found: {reference}")

    # Determine output path
    if output_dir is None:
        output_dir = vcf_in.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    normalized_vcf = output_dir / f"{vcf_in.stem}.normalized.vcf.gz"

    try:
        # Build bcftools norm command
        cmd = [
            "bcftools",
            "norm",
            "-f",
            str(reference),
            "-m",
            "-both",
            "--atomize",
            "-cw",  # Check REF alleles and warn about issues (instead of error)
            "-d",
            "exact",  # Remove exact duplicates (replaces deprecated -D)
            "-O",
            "z",  # Compress output
            "-o",
            str(normalized_vcf),
            str(vcf_in),
        ]

        logger.debug(f"Running bcftools norm: {' '.join(cmd)}")

        result = subprocess.run(cmd, capture_output=True, text=True, check=True)

        if result.stderr:
            logger.debug(f"bcftools norm stderr: {result.stderr}")

        # Index the compressed VCF
        tabix_cmd = ["tabix", "-p", "vcf", str(normalized_vcf)]
        logger.debug(f"Indexing VCF: {' '.join(tabix_cmd)}")

        subprocess.run(tabix_cmd, capture_output=True, text=True, check=True)

        logger.info(f"VCF normalized and indexed: {normalized_vcf}")
        return normalized_vcf

    except subprocess.CalledProcessError as e:
        error_msg = f"VCF normalization failed: {e.stderr if e.stderr else str(e)}"
        logger.error(error_msg)
        raise VCFOperationError(error_msg) from e
    except Exception as e:
        error_msg = f"Unexpected error during VCF normalization: {str(e)}"
        logger.error(error_msg)
        raise VCFOperationError(error_msg) from e


def classify_variant(ref: str, alt: str) -> VariantClass:
    """
    Classify variant type based on REF and ALT alleles.

    Args:
        ref: Reference allele
        alt: Alternative allele

    Returns:
        VariantClass enum value
    """
    # Skip symbolic alleles and IUPAC codes
    if any(c in ref + alt for c in ["<", ">", "*", "N"]):
        return VariantClass.UNKNOWN

    # Skip non-ATCG characters
    valid_bases = set("ATCG")
    if not (set(ref.upper()) <= valid_bases and set(alt.upper()) <= valid_bases):
        return VariantClass.UNKNOWN

    ref_len = len(ref)
    alt_len = len(alt)

    if ref_len == alt_len == 1:
        return VariantClass.SNV
    elif alt_len > ref_len:
        return VariantClass.INS
    elif alt_len < ref_len:
        return VariantClass.DEL
    else:
        # Complex variants
        return VariantClass.UNKNOWN


def extract_maf_from_record(record: Any) -> Optional[float]:
    """
    Extract minor allele frequency from VCF record using precedence: AD → DP4 → AF.

    Args:
        record: pysam VariantRecord

    Returns:
        MAF as float between 0 and 1, or None if cannot be determined
    """
    try:
        # Precedence 1: FORMAT/AD (Allelic Depth)
        if "AD" in record.format:
            samples = list(record.samples)
            if samples:
                sample = record.samples[samples[0]]  # Single sample
                ad = sample.get("AD")
                dp = sample.get("DP")  # Total depth

                if ad is not None and dp is not None:
                    # Handle different AD formats
                    if isinstance(ad, (list, tuple)) and len(ad) >= 2:
                        # Standard format: [ref_depth, alt_depth, ...]
                        ref_depth = int(ad[0])
                        alt_depth = sum(int(x) for x in ad[1:])  # Sum all alt alleles
                        total_depth = ref_depth + alt_depth
                    elif isinstance(ad, (int, str)):
                        # BBTools format: AD is just alt_depth, use DP for total
                        alt_depth = int(ad)
                        total_depth = int(dp)
                        ref_depth = total_depth - alt_depth
                    else:
                        # Skip if AD format is unexpected
                        alt_depth = ref_depth = total_depth = 0

                    if total_depth > 0 and ref_depth >= 0 and alt_depth >= 0:
                        ref_freq = ref_depth / total_depth
                        alt_freq = alt_depth / total_depth
                        return min(ref_freq, alt_freq)  # MAF is the minor frequency

        # Precedence 2: INFO/DP4 (legacy format from bcftools/samtools)
        if "DP4" in record.info:
            dp4 = record.info["DP4"]
            if len(dp4) == 4:
                ref_depth = int(dp4[0]) + int(dp4[1])  # Forward + reverse ref
                alt_depth = int(dp4[2]) + int(dp4[3])  # Forward + reverse alt
                total_depth = ref_depth + alt_depth
                if total_depth > 0:
                    ref_freq = ref_depth / total_depth
                    alt_freq = alt_depth / total_depth
                    return min(ref_freq, alt_freq)

        # Precedence 3: INFO/AF (Allele Frequency)
        if "AF" in record.info:
            af = record.info["AF"]
            if isinstance(af, (list, tuple)):
                # Multi-allelic: compute site-level MAF
                if len(af) > 0:
                    max_af = float(max(af))
                    ref_freq = 1.0 - sum(af)
                    return min(ref_freq, max_af)
            else:
                # Single alt allele
                af_val = float(af)
                ref_freq = 1.0 - af_val
                return min(ref_freq, af_val)

        # No MAF data available
        warnings.warn(f"No MAF data found for variant at {record.chrom}:{record.pos}")
        return None

    except (ValueError, TypeError, KeyError) as e:
        warnings.warn(f"Error extracting MAF from {record.chrom}:{record.pos}: {e}")
        return None


def parse_vcf_sites(
    vcf_path: Path, included_contigs: Optional[Set[str]] = None
) -> Iterator[SiteRecord]:
    """
    Parse normalized VCF and yield site records with MAF and variant classification.

    Args:
        vcf_path: Path to VCF file (can be compressed)
        included_contigs: Set of contig names to include (None = include all)

    Yields:
        SiteRecord objects for each variant site
    """
    try:
        with pysam.VariantFile(str(vcf_path)) as vcf:
            for record in vcf:
                # Filter by included contigs
                if (
                    included_contigs is not None
                    and record.chrom not in included_contigs
                ):
                    continue

                # Skip if no ALT alleles or no REF
                if not record.alts or not record.ref:
                    continue

                # For multi-allelic sites (after normalization should be rare),
                # process each ALT allele separately
                for alt_allele in record.alts:
                    # At this point record.ref is guaranteed to be non-None
                    assert record.ref is not None
                    variant_class = classify_variant(record.ref, alt_allele)

                    # Skip unknown/symbolic variants
                    if variant_class == VariantClass.UNKNOWN:
                        continue

                    # Extract MAF
                    maf = extract_maf_from_record(record)
                    if maf is None:
                        continue

                    # Get depth from DP field or calculate from AD
                    depth = 0
                    if "DP" in record.info:
                        depth = record.info["DP"]
                    elif "AD" in record.format:
                        samples = list(record.samples)
                        if samples:
                            sample = record.samples[samples[0]]
                            ad = sample.get("AD")
                            if ad is not None:
                                depth = sum(ad)

                    # Get original filter
                    orig_filter = (
                        ";".join(str(f) for f in record.filter.keys())
                        if record.filter.keys()
                        else "PASS"
                    )

                    yield SiteRecord(
                        chrom=record.chrom,
                        pos=record.pos,
                        ref=record.ref,
                        alt=alt_allele,
                        variant_class=variant_class,
                        depth=depth,
                        maf=maf,
                        original_filter=orig_filter,
                    )

    except Exception as e:
        error_msg = f"Error parsing VCF file {vcf_path}: {str(e)}"
        logger.error(error_msg)
        raise VCFOperationError(error_msg) from e


class AmbigGrid:
    """
    100×51 cumulative grid for ambiguous site counting.

    Tracks sites by depth (0-100, with capping) and MAF bins (0-50, representing 0.00-0.50).
    MAF binning: bin = floor(100 * MAF), capped at 50.
    """

    def __init__(self, dp_cap: int = 100):
        """
        Initialize grid.

        Args:
            dp_cap: Maximum depth value (higher values are capped)
        """
        self.dp_cap = dp_cap
        self.maf_bins = 51  # 0.00 to 0.50 in 0.01 increments

        # Grid: [depth][maf_bin] = count
        self.grid = np.zeros((dp_cap + 1, self.maf_bins), dtype=int)

    def add_site(self, depth: int, maf: float) -> None:
        """
        Add a site to the grid.

        Args:
            depth: Site depth (will be capped at dp_cap)
            maf: Minor allele frequency (0.0 to 1.0)
        """
        # Cap depth
        if depth > self.dp_cap:
            depth = self.dp_cap
            logger.debug(f"Depth capped at {self.dp_cap}")

        # Calculate MAF bin (floor of 100 * MAF)
        maf_bin = int(np.floor(100 * maf))

        # Cap MAF bin at 50 (representing 0.50)
        if maf_bin >= self.maf_bins:
            maf_bin = self.maf_bins - 1

        # Ensure non-negative
        depth = max(0, depth)
        maf_bin = max(0, maf_bin)

        self.grid[depth, maf_bin] += 1

    def build_cumulative(self) -> np.ndarray:
        """
        Build cumulative counts matrix.

        Returns:
            Cumulative grid where each cell contains count of sites
            with depth >= row_index and MAF >= col_index/100
        """
        # Start from bottom-right and work backwards
        cumulative = np.zeros_like(self.grid)

        for dp in range(self.dp_cap, -1, -1):
            for maf_bin in range(self.maf_bins - 1, -1, -1):
                cumulative[dp, maf_bin] = self.grid[dp, maf_bin]

                # Add counts from higher depth (same MAF bin)
                if dp < self.dp_cap:
                    cumulative[dp, maf_bin] += cumulative[dp + 1, maf_bin]

                # Add counts from higher MAF bin (same depth)
                if maf_bin < self.maf_bins - 1:
                    cumulative[dp, maf_bin] += cumulative[dp, maf_bin + 1]

                # Subtract double-counted intersection
                if dp < self.dp_cap and maf_bin < self.maf_bins - 1:
                    cumulative[dp, maf_bin] -= cumulative[dp + 1, maf_bin + 1]

        return cumulative

    def count_at(self, dp_min: int, maf_min: float) -> int:
        """
        Count sites meeting minimum thresholds.

        Args:
            dp_min: Minimum depth threshold
            maf_min: Minimum MAF threshold

        Returns:
            Count of sites with depth >= dp_min and MAF >= maf_min
        """
        cumulative = self.build_cumulative()

        # Convert MAF to bin (floor of 100 * MAF)
        maf_bin = int(np.floor(100 * maf_min))
        maf_bin = min(maf_bin, self.maf_bins - 1)

        # Cap depth
        dp_min = min(dp_min, self.dp_cap)

        return int(cumulative[dp_min, maf_bin])

    def to_wide_tsv(self, output_path: Path) -> None:
        """
        Write cumulative grid to wide TSV format.

        Per spec.md §4.4: "Wide table: rows depth=1..100; columns maf_0..maf_50"

        Args:
            output_path: Output file path
        """
        cumulative = self.build_cumulative()

        # Create header: depth, then MAF columns (0.00, 0.01, ..., 0.50)
        maf_headers = [f"{i/100:.2f}" for i in range(self.maf_bins)]
        header = ["depth"] + maf_headers

        with open(output_path, "w") as f:
            # Write header
            f.write("\t".join(header) + "\n")

            # Write data rows - start from depth=1, not depth=0 (per spec §4.4)
            for dp in range(1, self.dp_cap + 1):
                row = [str(dp)] + [
                    str(cumulative[dp, maf_bin]) for maf_bin in range(self.maf_bins)
                ]
                f.write("\t".join(row) + "\n")

        logger.info(f"Cumulative grid written to {output_path}")


def count_ambiguous_sites(
    vcf_path: Path,
    dp_min: int,
    maf_min: float,
    dp_cap: int = 100,
    included_contigs: Optional[Set[str]] = None,
    variant_classes: Optional[List[VariantClass]] = None,
) -> Tuple[int, AmbigGrid]:
    """
    Count ambiguous sites from normalized VCF.

    Args:
        vcf_path: Path to normalized VCF file
        dp_min: Minimum depth threshold
        maf_min: Minimum MAF threshold
        dp_cap: Maximum depth (higher values capped)
        included_contigs: Set of contigs to include
        variant_classes: List of variant classes to count (default: [SNV])

    Returns:
        Tuple of (ambiguous_count, grid_object)
    """
    if variant_classes is None:
        variant_classes = [VariantClass.SNV]

    grid = AmbigGrid(dp_cap=dp_cap)

    # Process all sites and build grid
    for site in parse_vcf_sites(vcf_path, included_contigs):
        if site.variant_class in variant_classes:
            grid.add_site(site.depth, site.maf)

    # Count ambiguous sites
    ambiguous_count = grid.count_at(dp_min, maf_min)

    logger.info(
        f"Found {ambiguous_count} ambiguous sites (dp>={dp_min}, maf>={maf_min})"
    )

    return ambiguous_count, grid


def emit_vcf(
    normalized_vcf_path: Path,
    output_path: Path,
    dp_min: int,
    maf_min: float,
    sample_name: str,
    require_pass: bool = False,
    included_contigs: Optional[Set[str]] = None,
) -> Path:
    """
    Emit VCF output with only records passing ambiguous site thresholds.

    Args:
        normalized_vcf_path: Path to normalized input VCF
        output_path: Output VCF path (will be bgzipped)
        dp_min: Minimum depth threshold
        maf_min: Minimum MAF threshold
        sample_name: Sample name for output
        require_pass: If True, include ORIG_FILTER info field
        included_contigs: Set of contigs to include

    Returns:
        Path to compressed, indexed VCF file

    Raises:
        VCFOperationError: If emission fails
    """
    try:
        # Ensure output path has .vcf.gz extension
        if not str(output_path).endswith(".vcf.gz"):
            output_path = output_path.with_suffix(".vcf.gz")

        output_path.parent.mkdir(parents=True, exist_ok=True)

        with pysam.VariantFile(str(normalized_vcf_path)) as input_vcf:
            # Create output VCF with enhanced header
            header = input_vcf.header.copy()

            # Add custom INFO fields
            header.add_line(
                '##INFO=<ID=AMBIG,Number=0,Type=Flag,Description="Site passes ambiguous thresholds">'
            )
            header.add_line(
                '##INFO=<ID=MAF,Number=1,Type=Float,Description="Minor allele frequency">'
            )
            header.add_line(
                '##INFO=<ID=MAF_BIN,Number=1,Type=Integer,Description="MAF bin (floor(100*MAF))">'
            )
            header.add_line(
                '##INFO=<ID=DP_CAP,Number=1,Type=Integer,Description="Depth cap used (100)">'
            )
            header.add_line(
                '##INFO=<ID=VARIANT_CLASS,Number=1,Type=String,Description="Variant class (SNV/INS/DEL)">'
            )

            if not require_pass:
                header.add_line(
                    '##INFO=<ID=ORIG_FILTER,Number=1,Type=String,Description="Original FILTER value">'
                )

            # Add FILTER definitions
            header.add_line('##FILTER=<ID=PASS,Description="All filters passed">')

            with pysam.VariantFile(str(output_path), "w", header=header) as output_vcf:
                records_written = 0

                for record in input_vcf:
                    # Filter by included contigs
                    if (
                        included_contigs is not None
                        and record.chrom not in included_contigs
                    ):
                        continue

                    # Skip if no ALT alleles
                    if not record.alts:
                        continue

                    # Process each ALT allele (though normalization should make these single)
                    for alt_allele in record.alts:
                        # Skip if ref or alt is None
                        if record.ref is None or alt_allele is None:
                            continue

                        variant_class = classify_variant(record.ref, alt_allele)

                        # Skip unknown/symbolic variants
                        if variant_class == VariantClass.UNKNOWN:
                            continue

                        # Extract MAF and depth
                        maf = extract_maf_from_record(record)
                        if maf is None:
                            continue

                        # Get depth
                        depth = 0
                        if "DP" in record.info:
                            depth = record.info["DP"]
                        elif "AD" in record.format:
                            samples = list(record.samples)
                            if samples:
                                sample = record.samples[samples[0]]
                                ad = sample.get("AD")
                                if ad is not None:
                                    depth = sum(ad)

                        # Check if passes thresholds
                        passes_thresholds = depth >= dp_min and maf >= maf_min

                        if passes_thresholds:
                            # Create output record
                            new_record = output_vcf.new_record(
                                contig=record.chrom,
                                start=record.start,
                                stop=record.stop,
                                alleles=(record.ref, alt_allele),
                                id=record.id,
                                qual=record.qual,
                            )

                            # Copy existing INFO fields
                            for key, value in record.info.items():
                                new_record.info[key] = value

                            # Add our custom INFO fields
                            new_record.info["AMBIG"] = True
                            new_record.info["MAF"] = round(maf, 4)
                            new_record.info["MAF_BIN"] = int(np.floor(100 * maf))
                            new_record.info["DP_CAP"] = 100
                            new_record.info["VARIANT_CLASS"] = variant_class.value

                            if not require_pass:
                                filter_keys = list(record.filter.keys())
                                orig_filter = (
                                    ";".join(str(k) for k in filter_keys)
                                    if filter_keys
                                    else "PASS"
                                )
                                new_record.info["ORIG_FILTER"] = orig_filter

                            # Set FILTER to PASS
                            new_record.filter.clear()
                            new_record.filter.add("PASS")

                            # Copy FORMAT fields
                            for sample_id in record.samples:
                                sample_name = str(sample_id)
                                for key in record.format:
                                    new_record.samples[sample_name][key] = (
                                        record.samples[sample_id][key]
                                    )

                            output_vcf.write(new_record)
                            records_written += 1

        # Index the compressed VCF
        tabix_cmd = ["tabix", "-p", "vcf", str(output_path)]
        subprocess.run(tabix_cmd, capture_output=True, text=True, check=True)

        logger.info(
            f"Emitted VCF with {records_written} ambiguous sites: {output_path}"
        )
        return output_path

    except Exception as e:
        error_msg = f"VCF emission failed: {str(e)}"
        logger.error(error_msg)
        raise VCFOperationError(error_msg) from e


def emit_bed(
    normalized_vcf_path: Path,
    output_path: Path,
    dp_min: int,
    maf_min: float,
    sample_name: str,
    included_contigs: Optional[Set[str]] = None,
) -> Path:
    """
    Emit BED output with ambiguous sites in 0-based half-open coordinates.

    BED format columns:
    chrom, start, end, name, score, strand, sample, variant_class, ref, alt, maf, dp, maf_bin, dp_cap

    Args:
        normalized_vcf_path: Path to normalized input VCF
        output_path: Output BED path (will be bgzipped)
        dp_min: Minimum depth threshold
        maf_min: Minimum MAF threshold
        sample_name: Sample name
        included_contigs: Set of contigs to include

    Returns:
        Path to compressed, indexed BED file

    Raises:
        VCFOperationError: If emission fails
    """
    try:
        # Ensure output path has .bed.gz extension
        if not str(output_path).endswith(".bed.gz"):
            output_path = output_path.with_suffix(".bed.gz")

        output_path.parent.mkdir(parents=True, exist_ok=True)

        bed_records = []

        with pysam.VariantFile(str(normalized_vcf_path)) as input_vcf:
            for record in input_vcf:
                # Filter by included contigs
                if (
                    included_contigs is not None
                    and record.chrom not in included_contigs
                ):
                    continue

                # Skip if no ALT alleles
                if not record.alts:
                    continue

                # Process each ALT allele
                for alt_allele in record.alts:
                    # Skip if ref or alt is None
                    if record.ref is None or alt_allele is None:
                        continue

                    variant_class = classify_variant(record.ref, alt_allele)

                    # Skip unknown/symbolic variants
                    if variant_class == VariantClass.UNKNOWN:
                        continue

                    # Extract MAF and depth
                    maf = extract_maf_from_record(record)
                    if maf is None:
                        continue

                    # Get depth
                    depth = 0
                    if "DP" in record.info:
                        depth = record.info["DP"]
                    elif "AD" in record.format:
                        samples = list(record.samples)
                        if samples:
                            sample = record.samples[samples[0]]
                            ad = sample.get("AD")
                            if ad is not None:
                                depth = sum(ad)

                    # Check if passes thresholds
                    passes_thresholds = depth >= dp_min and maf >= maf_min

                    if passes_thresholds:
                        # Convert to 0-based coordinates
                        # VCF is 1-based, BED is 0-based half-open
                        start = record.start  # pysam already converts to 0-based

                        if variant_class == VariantClass.SNV:
                            # SNV: 1bp interval
                            end = start + 1
                        elif variant_class == VariantClass.INS:
                            # Insertion: 1bp anchor position
                            end = start + 1
                        elif variant_class == VariantClass.DEL:
                            # Deletion: span the deleted bases
                            if record.ref is not None:
                                end = start + len(record.ref)
                            else:
                                end = start + 1
                        else:
                            end = start + 1

                        # Create BED record
                        name = f"{record.chrom}:{record.pos}_{record.ref}>{alt_allele}"
                        score = min(1000, int(1000 * maf))  # Scale MAF to 0-1000
                        maf_bin = int(np.floor(100 * maf))

                        bed_record = [
                            record.chrom,  # chrom
                            str(start),  # start (0-based)
                            str(end),  # end (0-based half-open)
                            name,  # name
                            str(score),  # score (0-1000)
                            ".",  # strand (not applicable)
                            sample_name,  # sample
                            variant_class.value,  # variant_class
                            record.ref,  # ref
                            alt_allele,  # alt
                            f"{maf:.4f}",  # maf
                            str(depth),  # dp
                            str(maf_bin),  # maf_bin
                            "100",  # dp_cap
                        ]

                        bed_records.append(bed_record)

        # Sort records by chromosome and position
        bed_records.sort(key=lambda x: (x[0], int(x[1])))

        # Write BED file
        temp_bed = output_path.with_suffix(".bed")
        with open(temp_bed, "w") as f:
            # Write header comment
            f.write(
                "# chrom\tstart\tend\tname\tscore\tstrand\tsample\tvariant_class\tref\talt\tmaf\tdp\tmaf_bin\tdp_cap\n"
            )

            # Write records
            for bed_record in bed_records:
                f.write("\t".join(bed_record) + "\n")

        # Compress with bgzip
        bgzip_cmd = ["bgzip", "-c", str(temp_bed)]
        with open(output_path, "w") as f:
            subprocess.run(bgzip_cmd, stdout=f, check=True)

        # Remove temporary file
        temp_bed.unlink()

        # Index with tabix
        tabix_cmd = ["tabix", "-p", "bed", str(output_path)]
        subprocess.run(tabix_cmd, capture_output=True, text=True, check=True)

        logger.info(
            f"Emitted BED with {len(bed_records)} ambiguous sites: {output_path}"
        )
        return output_path

    except Exception as e:
        error_msg = f"BED emission failed: {str(e)}"
        logger.error(error_msg)
        raise VCFOperationError(error_msg) from e


def emit_matrix(grid: AmbigGrid, output_path: Path, sample_name: str) -> Path:
    """
    Emit variant matrix as compressed TSV.

    Args:
        grid: Populated AmbigGrid with variant counts
        output_path: Output path (will be ensured to end with .tsv.gz)
        sample_name: Sample name for logging

    Returns:
        Path to compressed matrix file

    Raises:
        VCFOperationError: If matrix emission fails
    """
    try:
        # Ensure output path has .tsv.gz extension
        if not str(output_path).endswith(".tsv.gz"):
            output_path = output_path.with_suffix(".tsv.gz")

        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Write to temporary uncompressed file first
        temp_tsv = output_path.with_suffix(".tsv")
        grid.to_wide_tsv(temp_tsv)

        # Compress with gzip
        import gzip

        with open(temp_tsv, "rb") as f_in:
            with gzip.open(output_path, "wb") as f_out:
                f_out.write(f_in.read())

        # Remove temporary file
        temp_tsv.unlink()

        logger.info(f"Emitted variant matrix for {sample_name}: {output_path}")
        return output_path

    except Exception as e:
        error_msg = f"Matrix emission failed: {str(e)}"
        logger.error(error_msg)
        raise VCFOperationError(error_msg) from e


def emit_per_contig(
    normalized_vcf_path: Path,
    depth_summary_path: Path,
    output_path: Path,
    dp_min: int,
    maf_min: float,
    sample_name: str,
    included_contigs: Optional[Set[str]] = None,
) -> Path:
    """
    Emit per-contig summary statistics.

    Args:
        normalized_vcf_path: Path to normalized VCF file
        depth_summary_path: Path to mosdepth summary file
        output_path: Output TSV path
        dp_min: Minimum depth threshold
        maf_min: Minimum MAF threshold
        sample_name: Sample name
        included_contigs: Set of contigs to include

    Returns:
        Path to per-contig summary file

    Raises:
        VCFOperationError: If emission fails
    """
    try:
        from .depth import parse_mosdepth_summary

        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Parse depth summary for contig-level stats
        depth_summary = parse_mosdepth_summary(depth_summary_path)

        # Count variants per contig
        contig_variant_counts = {}

        with pysam.VariantFile(str(normalized_vcf_path)) as vcf:
            for record in vcf:
                # Filter by included contigs
                if (
                    included_contigs is not None
                    and record.chrom not in included_contigs
                ):
                    continue

                # Skip if no ALT alleles
                if not record.alts:
                    continue

                # Initialize contig counts if needed
                if record.chrom not in contig_variant_counts:
                    contig_variant_counts[record.chrom] = {
                        "snv": 0,
                        "indel": 0,
                        "del": 0,
                    }

                # Process each ALT allele
                for alt_allele in record.alts:
                    if record.ref is None or alt_allele is None:
                        continue

                    variant_class = classify_variant(record.ref, alt_allele)

                    if variant_class == VariantClass.UNKNOWN:
                        continue

                    # Extract MAF and depth
                    maf = extract_maf_from_record(record)
                    if maf is None:
                        continue

                    # Get depth
                    depth = 0
                    if "DP" in record.info:
                        depth = int(record.info["DP"])
                    elif "AD" in record.format:
                        samples = list(record.samples)
                        if samples:
                            sample = record.samples[samples[0]]
                            ad = sample.get("AD")
                            if ad is not None:
                                depth = sum(int(x) for x in ad)

                    # Check if passes thresholds
                    if depth >= dp_min and maf >= maf_min:
                        if variant_class == VariantClass.SNV:
                            contig_variant_counts[record.chrom]["snv"] += 1
                        elif variant_class == VariantClass.INS:
                            contig_variant_counts[record.chrom]["indel"] += 1
                        elif variant_class == VariantClass.DEL:
                            contig_variant_counts[record.chrom]["del"] += 1

        # Write per-contig summary
        with open(output_path, "w") as f:
            # Write header
            header = [
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
            f.write("\t".join(header) + "\n")

            # Write data for each contig in depth summary
            for contig_stat in depth_summary.contig_stats:
                contig_name = contig_stat.name
                length = contig_stat.length
                callable_bases = contig_stat.bases_covered
                mean_depth = contig_stat.mean_depth
                breadth_10x = contig_stat.breadth_10x

                # Get variant counts for this contig
                counts = contig_variant_counts.get(
                    contig_name, {"snv": 0, "indel": 0, "del": 0}
                )
                snv_count = counts["snv"]
                indel_count = counts["indel"]
                del_count = counts["del"]

                # Calculate SNV per Mb - use callable_bases, not total length
                # Per spec.md §4.1: ambiguous_snv_per_mb = ambiguous_snv_count * 1e6 / callable_bases
                snv_per_mb = (
                    (snv_count * 1_000_000) / callable_bases
                    if callable_bases > 0
                    else 0.0
                )

                row = [
                    sample_name,
                    contig_name,
                    str(length),
                    str(callable_bases),
                    f"{breadth_10x:.4f}",
                    str(snv_count),
                    str(indel_count),
                    str(del_count),
                    f"{snv_per_mb:.2f}",
                    f"{mean_depth:.2f}",
                ]
                f.write("\t".join(row) + "\n")

        logger.info(f"Emitted per-contig summary for {sample_name}: {output_path}")
        return output_path

    except Exception as e:
        error_msg = f"Per-contig summary emission failed: {str(e)}"
        logger.error(error_msg)
        raise VCFOperationError(error_msg) from e


def emit_multiqc(
    sample_name: str,
    ambiguous_snv_count: int,
    breadth_10x: float,
    callable_bases: int,
    genome_length: int,
    dp_min: int,
    maf_min: float,
    mapper: str,
    caller: str,
    mode: str,
    output_path: Path,
) -> Path:
    """
    Emit MultiQC-compatible metrics TSV.

    As per spec.md §4.6: sample, ambiguous_snv_count, ambiguous_snv_per_mb, breadth_10x,
    callable_bases, genome_length, dp_min, maf_min, mapper, caller, mode

    Args:
        sample_name: Sample identifier
        ambiguous_snv_count: Count of ambiguous SNVs
        breadth_10x: Fraction of genome with ≥10x coverage
        callable_bases: Number of callable bases
        genome_length: Total genome length
        dp_min: Minimum depth threshold used
        maf_min: Minimum MAF threshold used
        mapper: Mapper tool name
        caller: Caller tool name
        mode: Analysis mode (self or ref)
        output_path: Output path (will be ensured to end with .tsv)

    Returns:
        Path to emitted MultiQC TSV file

    Raises:
        VCFOperationError: If MultiQC emission fails
    """
    try:
        # Ensure output path has .tsv extension
        if not str(output_path).endswith(".tsv"):
            output_path = output_path.with_suffix(".tsv")

        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Calculate ambiguous SNVs per megabase
        # Per spec.md §4.1: ambiguous_snv_per_mb = ambiguous_snv_count * 1e6 / callable_bases
        ambiguous_snv_per_mb = (
            (ambiguous_snv_count * 1_000_000) / callable_bases
            if callable_bases > 0
            else 0.0
        )

        with open(output_path, "w") as f:
            # Write header
            header = [
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
            f.write("\t".join(header) + "\n")

            # Write data row
            row = [
                sample_name,
                str(ambiguous_snv_count),
                f"{ambiguous_snv_per_mb:.2f}",
                f"{breadth_10x:.4f}",
                str(callable_bases),
                str(genome_length),
                str(dp_min),
                f"{maf_min:.3f}",
                mapper,
                caller,
                mode,
            ]
            f.write("\t".join(row) + "\n")

        logger.info(f"Emitted MultiQC metrics for {sample_name}: {output_path}")
        return output_path

    except Exception as e:
        error_msg = f"MultiQC emission failed: {str(e)}"
        logger.error(error_msg)
        raise VCFOperationError(error_msg) from e
