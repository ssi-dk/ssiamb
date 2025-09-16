"""
VCF operations module for ssiamb.

This module implements VCF normalization, atomization, MAF extraction,
variant classification, and grid-based counting for ambiguous sites.
"""

import logging
import subprocess
import shutil
import tempfile
import warnings
from pathlib import Path
from typing import Optional, List, Dict, Tuple, Iterator, Set
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


def normalize_and_split(vcf_in: Path, reference: Path, 
                       output_dir: Optional[Path] = None) -> Path:
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
            "bcftools", "norm",
            "-f", str(reference),
            "-m", "-both",
            "--atomize",
            "-O", "z",  # Compress output
            "-o", str(normalized_vcf),
            str(vcf_in)
        ]
        
        logger.debug(f"Running bcftools norm: {' '.join(cmd)}")
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )
        
        if result.stderr:
            logger.debug(f"bcftools norm stderr: {result.stderr}")
        
        # Index the compressed VCF
        tabix_cmd = ["tabix", "-p", "vcf", str(normalized_vcf)]
        logger.debug(f"Indexing VCF: {' '.join(tabix_cmd)}")
        
        subprocess.run(
            tabix_cmd,
            capture_output=True,
            text=True,
            check=True
        )
        
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


def extract_maf_from_record(record) -> Optional[float]:
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
                if ad is not None and len(ad) >= 2:
                    ref_depth = ad[0]
                    alt_depth = sum(ad[1:])  # Sum all alt alleles
                    total_depth = ref_depth + alt_depth
                    if total_depth > 0:
                        ref_freq = ref_depth / total_depth
                        alt_freq = alt_depth / total_depth
                        return min(ref_freq, alt_freq)  # MAF is the minor frequency
        
        # Precedence 2: INFO/DP4 (legacy format from bcftools/samtools)
        if "DP4" in record.info:
            dp4 = record.info["DP4"]
            if len(dp4) == 4:
                ref_depth = dp4[0] + dp4[1]  # Forward + reverse ref
                alt_depth = dp4[2] + dp4[3]  # Forward + reverse alt
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
                    max_af = max(af)
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


def parse_vcf_sites(vcf_path: Path, included_contigs: Optional[Set[str]] = None) -> Iterator[SiteRecord]:
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
                if included_contigs is not None and record.chrom not in included_contigs:
                    continue
                
                # Skip if no ALT alleles
                if not record.alts:
                    continue
                
                # For multi-allelic sites (after normalization should be rare),
                # process each ALT allele separately
                for alt_allele in record.alts:
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
                    orig_filter = ";".join(record.filter.keys()) if record.filter.keys() else "PASS"
                    
                    yield SiteRecord(
                        chrom=record.chrom,
                        pos=record.pos,
                        ref=record.ref,
                        alt=alt_allele,
                        variant_class=variant_class,
                        depth=depth,
                        maf=maf,
                        original_filter=orig_filter
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
        
        Args:
            output_path: Output file path
        """
        cumulative = self.build_cumulative()
        
        # Create header: depth, then MAF columns (0.00, 0.01, ..., 0.50)
        maf_headers = [f"{i/100:.2f}" for i in range(self.maf_bins)]
        header = ["depth"] + maf_headers
        
        with open(output_path, 'w') as f:
            # Write header
            f.write('\t'.join(header) + '\n')
            
            # Write data rows
            for dp in range(self.dp_cap + 1):
                row = [str(dp)] + [str(cumulative[dp, maf_bin]) for maf_bin in range(self.maf_bins)]
                f.write('\t'.join(row) + '\n')
        
        logger.info(f"Cumulative grid written to {output_path}")


def count_ambiguous_sites(vcf_path: Path, dp_min: int, maf_min: float, 
                         dp_cap: int = 100, 
                         included_contigs: Optional[Set[str]] = None,
                         variant_classes: Optional[List[VariantClass]] = None) -> Tuple[int, AmbigGrid]:
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
    
    logger.info(f"Found {ambiguous_count} ambiguous sites (dp>={dp_min}, maf>={maf_min})")
    
    return ambiguous_count, grid