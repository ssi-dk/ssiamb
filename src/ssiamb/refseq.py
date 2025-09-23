"""
RefSeq genome download functionality.

This module handles downloading reference genomes from NCBI RefSeq
according to the specified selection policy.
"""

import re
import time
import subprocess
from pathlib import Path
from typing import List, Dict, Optional, Tuple, Any
from dataclasses import dataclass
from datetime import datetime
import requests
import gzip
import shutil
from rich.console import Console
from rich.progress import Progress, BarColumn, TextColumn, TimeRemainingColumn

console = Console()


@dataclass
class RefSeqGenome:
    """Represents a RefSeq genome with metadata for selection."""

    accession: str
    organism_name: str
    category: str  # reference, representative, na
    assembly_level: str  # Complete Genome, Chromosome, Scaffold, Contig
    assembly_type: str  # RefSeq or GenBank
    release_date: str
    contig_count: int
    total_length: int
    ftp_path: str
    has_plasmids: bool = False

    def get_selection_priority(self) -> Tuple[int, int, int, str, int]:
        """
        Return a tuple for sorting by selection policy.
        Lower values = higher priority.

        Returns: (category_priority, level_priority, type_priority, date_desc, contig_count)
        """
        # Category priority: reference=0, representative=1, other=2
        category_map = {"reference": 0, "representative": 1, "na": 2}
        category_priority = category_map.get(self.category.lower(), 2)

        # Assembly level priority: Complete=0, Chromosome=1, Scaffold=2, Contig=3
        level_map = {"complete genome": 0, "chromosome": 1, "scaffold": 2, "contig": 3}
        level_priority = level_map.get(self.assembly_level.lower(), 3)

        # Type priority: RefSeq=0, GenBank=1
        type_priority = 0 if self.assembly_type == "RefSeq" else 1

        # Date (descending - newer first), convert to negative for sorting
        try:
            date_obj = datetime.strptime(self.release_date, "%Y/%m/%d")
            date_desc = date_obj.strftime("%Y%m%d")
            date_desc = "99999999" + date_desc[::-1]  # Reverse for desc sort
        except ValueError:
            date_desc = "00000000"

        return (
            category_priority,
            level_priority,
            type_priority,
            date_desc,
            self.contig_count,
        )


class RefSeqDownloader:
    """Handles RefSeq genome downloading with best selection policy."""

    ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    ASSEMBLY_URL = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=assembly"
    )

    def __init__(
        self, output_dir: Path, verbose: bool = False, create_indexes: bool = True
    ):
        """Initialize downloader with output directory."""
        self.output_dir = Path(output_dir)
        self.verbose = verbose
        self.create_indexes = create_indexes
        self.session = requests.Session()
        # Set user agent for NCBI API
        self.session.headers.update(
            {"User-Agent": "ssiamb/0.1.0 (https://github.com/ssi-dk/ssiamb)"}
        )

    def normalize_species_name(self, species: str) -> str:
        """Normalize species name to Genus_species format."""
        # Remove any non-alphanumeric characters except underscore and space
        cleaned = re.sub(r"[^\w\s]", "", species.strip())
        # Replace spaces with underscores
        normalized = re.sub(r"\s+", "_", cleaned)
        # Ensure format is Genus_species (capitalize first letter of each part)
        parts = normalized.split("_")
        if len(parts) >= 2:
            return f"{parts[0].capitalize()}_{parts[1].lower()}"
        return normalized.capitalize()

    def search_assemblies(self, species: str) -> List[str]:
        """Search for assembly IDs for a given species."""
        # Convert Genus_species to "Genus species" for NCBI search
        search_term = species.replace("_", " ")

        # First try to find reference genomes specifically
        params_ref: Dict[str, Any] = {
            "db": "assembly",
            "term": f'"{search_term}"[Organism] AND "reference genome"[RefSeq Category]',
            "retmode": "json",
            "retmax": 50,
        }

        # Then try representative genomes
        params_rep: Dict[str, Any] = {
            "db": "assembly",
            "term": f'"{search_term}"[Organism] AND "representative genome"[RefSeq Category]',
            "retmode": "json",
            "retmax": 50,
        }

        # Finally try all genomes with some quality filters
        params_all: Dict[str, Any] = {
            "db": "assembly",
            "term": f'"{search_term}"[Organism] AND "complete genome"[Assembly Level]',
            "retmode": "json",
            "retmax": 100,
        }

        all_ids = []

        try:
            # Try reference genomes first
            for params in [params_ref, params_rep, params_all]:
                time.sleep(0.4)  # Rate limiting
                response = self.session.get(self.ESEARCH_URL, params=params, timeout=30)
                response.raise_for_status()
                data = response.json()

                if "esearchresult" in data and "idlist" in data["esearchresult"]:
                    assembly_ids = data["esearchresult"]["idlist"]
                    all_ids.extend(assembly_ids)
                    if self.verbose:
                        console.print(
                            f"Found {len(assembly_ids)} assemblies with search: {params['term']}"
                        )

                    # If we found some reference or representative genomes, stop here
                    if assembly_ids and (
                        "reference genome" in str(params["term"])
                        or "representative genome" in str(params["term"])
                    ):
                        break

            # Remove duplicates while preserving order
            seen = set()
            unique_ids = []
            for id in all_ids:
                if id not in seen:
                    seen.add(id)
                    unique_ids.append(id)

            if self.verbose:
                console.print(
                    f"Found {len(unique_ids)} total unique assemblies for {species}"
                )

            return unique_ids

        except Exception as e:
            console.print(f"[red]Error searching for {species}: {e}[/red]")
            return []

    def get_assembly_metadata(self, assembly_ids: List[str]) -> List[RefSeqGenome]:
        """Get detailed metadata for assembly IDs."""
        if not assembly_ids:
            return []

        # NCBI API URL length limit is around 8000 chars, so use smaller batches
        # Each ID is about 8 chars plus comma, so ~50 IDs per batch is safe
        genomes = []
        for i in range(0, len(assembly_ids), 50):
            batch_ids = assembly_ids[i : i + 50]
            batch_genomes = self._get_batch_metadata(batch_ids)
            genomes.extend(batch_genomes)

        return genomes

    def _get_batch_metadata(self, assembly_ids: List[str]) -> List[RefSeqGenome]:
        """Get metadata for a batch of assembly IDs."""
        params = {"db": "assembly", "id": ",".join(assembly_ids), "retmode": "json"}

        try:
            # Add delay to respect NCBI rate limits (max 3 requests per second)
            time.sleep(0.4)

            response = self.session.get(self.ESUMMARY_URL, params=params, timeout=60)
            response.raise_for_status()
            data = response.json()

            genomes = []
            if "result" in data:
                for assembly_id in assembly_ids:
                    if assembly_id in data["result"]:
                        genome_data = data["result"][assembly_id]
                        genome = self._parse_genome_data(genome_data)
                        if genome:
                            genomes.append(genome)

            return genomes

        except Exception as e:
            console.print(f"[red]Error getting assembly metadata: {e}[/red]")
            return []

    def _parse_genome_data(self, data: Dict) -> Optional[RefSeqGenome]:
        """Parse genome data from NCBI assembly summary."""
        try:
            # Check if this is a RefSeq genome (has RefSeq accession)
            accession = data.get("assemblyaccession", "")
            if not accession.startswith("GCF_"):
                return None  # Skip GenBank-only assemblies

            # Extract metadata using correct field names from NCBI API
            organism_name = data.get("organism", data.get("speciesname", ""))

            # Assembly level is in assemblystatus
            assembly_status = data.get("assemblystatus", "").lower()
            assembly_level = (
                assembly_status  # Will be like "complete genome", "chromosome", etc.
            )

            release_date = data.get(
                "asmreleasedate_refseq", data.get("asmreleasedate_genbank", "")
            )

            # Use contign50 as a proxy for assembly quality (higher N50 = fewer contigs for same size)
            contig_n50 = int(data.get("contign50", 0))
            scaffold_n50 = int(data.get("scaffoldn50", 0))
            # Use the better of the two N50 values, prefer scaffold N50
            best_n50 = max(contig_n50, scaffold_n50)
            # Convert N50 to approximate contig count (inverse relationship)
            # For sorting, we'll use negative N50 so higher N50 (fewer contigs) sorts first
            contig_count = -best_n50 if best_n50 > 0 else 999999

            total_length = 0  # Not directly available in this API response
            ftp_path = data.get("ftppath_refseq", "")

            # Determine category (reference/representative/other)
            category = "na"
            refseq_category = data.get("refseq_category", "").lower()
            if "reference" in refseq_category:
                category = "reference"
            elif "representative" in refseq_category:
                category = "representative"

            # Check for plasmids in the assembly name
            assembly_name = data.get("assemblyname", "")
            has_plasmids = "plasmid" in assembly_name.lower()

            return RefSeqGenome(
                accession=accession,
                organism_name=organism_name,
                category=category,
                assembly_level=assembly_level,
                assembly_type="RefSeq",
                release_date=release_date,
                contig_count=contig_count,
                total_length=total_length,
                ftp_path=ftp_path,
                has_plasmids=has_plasmids,
            )

        except (ValueError, KeyError) as e:
            if self.verbose:
                console.print(
                    f"[yellow]Warning: Could not parse genome data: {e}[/yellow]"
                )
            return None

    def select_best_genome(self, genomes: List[RefSeqGenome]) -> Optional[RefSeqGenome]:
        """Select the best genome according to the selection policy."""
        if not genomes:
            return None

        # Filter out genomes without FTP paths
        valid_genomes = [g for g in genomes if g.ftp_path]
        if not valid_genomes:
            return None

        if self.verbose:
            console.print(
                f"Found {len(valid_genomes)} valid RefSeq genomes, selecting best..."
            )
            for g in valid_genomes[:5]:  # Show top 5
                n50_val = -g.contig_count if g.contig_count < 0 else g.contig_count
                console.print(
                    f"  {g.accession}: {g.category}, {g.assembly_level}, N50={n50_val}"
                )

        # Sort by selection criteria
        # Priority: category, assembly_level, type, date (desc), contig_count (asc)
        sorted_genomes = sorted(valid_genomes, key=lambda g: g.get_selection_priority())

        best_genome = sorted_genomes[0]

        if self.verbose:
            n50_display = (
                -best_genome.contig_count if best_genome.contig_count < 0 else "unknown"
            )
            console.print(
                f"Selected: {best_genome.accession} ({best_genome.category}, "
                f"{best_genome.assembly_level}, N50={n50_display})"
            )

        return best_genome

    def download_genome(
        self, genome: RefSeqGenome, species: str, dry_run: bool = False
    ) -> Optional[Path]:
        """Download a genome FASTA file."""
        if not genome.ftp_path:
            console.print(f"[red]No FTP path available for {genome.accession}[/red]")
            return None

        # Convert FTP URL to HTTPS URL for downloading
        ftp_path = genome.ftp_path
        if ftp_path.startswith("ftp://ftp.ncbi.nlm.nih.gov"):
            https_path = ftp_path.replace(
                "ftp://ftp.ncbi.nlm.nih.gov", "https://ftp.ncbi.nlm.nih.gov"
            )
        else:
            https_path = ftp_path

        output_file = self.output_dir / f"{species}.fna"

        # Try different filename patterns that NCBI uses
        # Extract assembly name from the FTP path
        assembly_name = https_path.split("/")[-1]  # e.g., "GCF_000005845.2_ASM584v2"

        # Possible filename patterns
        filename_patterns = [
            f"{assembly_name}_genomic.fna.gz",  # Most common pattern
            f"{genome.accession}_genomic.fna.gz",  # Alternative pattern
        ]

        if dry_run:
            fasta_url = f"{https_path}/{filename_patterns[0]}"
            console.print(
                f"[yellow]DRY RUN: Would download {fasta_url} to {output_file}[/yellow]"
            )
            return output_file

        # Try each filename pattern
        for filename in filename_patterns:
            fasta_url = f"{https_path}/{filename}"

            try:
                # First check if the file exists (HEAD request)
                head_response = self.session.head(fasta_url, timeout=30)
                if head_response.status_code == 200:
                    # Found the file, proceed with download
                    console.print(f"Downloading {genome.accession} for {species}...")

                    # Download with progress bar
                    response = self.session.get(fasta_url, stream=True, timeout=300)
                    response.raise_for_status()

                    total_size = int(response.headers.get("content-length", 0))

                    with Progress(
                        TextColumn("[progress.description]{task.description}"),
                        BarColumn(),
                        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
                        TimeRemainingColumn(),
                        console=console,
                    ) as progress:
                        task = progress.add_task(
                            f"Downloading {species}", total=total_size
                        )

                        # Download to temporary file first (atomic write)
                        temp_file = output_file.with_suffix(".tmp.gz")
                        with open(temp_file, "wb") as f:
                            for chunk in response.iter_content(chunk_size=8192):
                                if chunk:
                                    f.write(chunk)
                                    progress.update(task, advance=len(chunk))

                        # Decompress the file
                        with gzip.open(temp_file, "rb") as f_in:
                            with open(output_file, "wb") as f_out:
                                shutil.copyfileobj(f_in, f_out)

                        # Clean up temporary file
                        temp_file.unlink()

                        console.print(
                            f"[green]Successfully downloaded {species} to {output_file}[/green]"
                        )
                        return output_file

                elif self.verbose:
                    console.print(
                        f"[dim]Tried {filename}: status {head_response.status_code}[/dim]"
                    )

            except Exception as e:
                if self.verbose:
                    console.print(f"[dim]Error trying {filename}: {e}[/dim]")
                continue

        # If we get here, none of the filename patterns worked
        if self.verbose:
            console.print(
                f"[yellow]No downloadable files found for {genome.accession}[/yellow]"
            )

        # Clean up any partial files
        if output_file.exists():
            output_file.unlink()
        temp_file = output_file.with_suffix(".tmp.gz")
        if temp_file.exists():
            temp_file.unlink()

        return None

    def download_species(self, species: str, dry_run: bool = False) -> Optional[Path]:
        """Download the best RefSeq genome for a species."""
        normalized_species = self.normalize_species_name(species)

        console.print(f"Processing {normalized_species}...")

        # Search for assemblies
        assembly_ids = self.search_assemblies(normalized_species)
        if not assembly_ids:
            console.print(
                f"[yellow]No assemblies found for {normalized_species}[/yellow]"
            )
            return None

        # Get metadata
        genomes = self.get_assembly_metadata(assembly_ids)
        if not genomes:
            console.print(
                f"[yellow]No RefSeq genomes found for {normalized_species}[/yellow]"
            )
            return None

        # Sort all genomes by selection criteria (best first)
        valid_genomes = [g for g in genomes if g.ftp_path]
        if not valid_genomes:
            console.print(
                f"[yellow]No genomes with valid FTP paths found for {normalized_species}[/yellow]"
            )
            return None

        sorted_genomes = sorted(valid_genomes, key=lambda g: g.get_selection_priority())

        if self.verbose:
            console.print(
                f"Found {len(sorted_genomes)} valid RefSeq genomes, trying in order..."
            )

        # Try downloading genomes in order of preference
        for i, genome in enumerate(sorted_genomes[:5]):  # Try top 5 candidates
            if self.verbose:
                n50_val = (
                    -genome.contig_count
                    if genome.contig_count < 0
                    else genome.contig_count
                )
                console.print(
                    f"  Trying #{i+1}: {genome.accession} ({genome.category}, {genome.assembly_level}, N50={n50_val})"
                )

            result = self.download_genome(genome, normalized_species, dry_run)
            if result:
                # Generate indexes after successful download (if enabled)
                if self.create_indexes:
                    if self.generate_indexes(result, normalized_species, dry_run):
                        console.print(
                            f"[green]✓ Successfully downloaded and indexed {normalized_species}[/green]"
                        )
                    else:
                        console.print(
                            f"[yellow]⚠ Downloaded {normalized_species} but failed to generate some indexes[/yellow]"
                        )
                else:
                    console.print(
                        f"[green]✓ Successfully downloaded {normalized_species}[/green]"
                    )
                return result
            elif not dry_run:
                console.print(
                    f"[yellow]Download failed for {genome.accession}, trying next option...[/yellow]"
                )

        console.print(
            f"[red]All download attempts failed for {normalized_species}[/red]"
        )
        return None

    def generate_indexes(
        self, fasta_file: Path, species: str, dry_run: bool = False
    ) -> bool:
        """Generate minimap2 and bwa-mem2 indexes for a FASTA file."""
        if dry_run:
            console.print(
                f"[yellow]DRY RUN: Would generate minimap2 and bwa-mem2 indexes for {fasta_file}[/yellow]"
            )
            return True

        if not fasta_file.exists():
            console.print(f"[red]FASTA file not found: {fasta_file}[/red]")
            return False

        success = True

        # Generate minimap2 index
        success &= self._generate_minimap2_index(fasta_file, species)

        # Generate bwa-mem2 index
        success &= self._generate_bwa_mem2_index(fasta_file, species)

        return success

    def _generate_minimap2_index(self, fasta_file: Path, species: str) -> bool:
        """Generate minimap2 index."""
        index_file = (
            fasta_file.parent / f"{fasta_file.name}.mmi"
        )  # Use .fna.mmi pattern

        try:
            console.print(f"Generating minimap2 index for {species}...")

            # minimap2 -d output.mmi input.fasta
            cmd = ["minimap2", "-d", str(index_file), str(fasta_file)]

            result = subprocess.run(
                cmd, capture_output=True, text=True, timeout=300  # 5 minute timeout
            )

            if result.returncode == 0:
                console.print(
                    f"[green]✓ Generated minimap2 index: {index_file}[/green]"
                )
                return True
            else:
                console.print(
                    f"[red]✗ Failed to generate minimap2 index: {result.stderr}[/red]"
                )
                return False

        except subprocess.TimeoutExpired:
            console.print(f"[red]✗ Minimap2 indexing timed out for {species}[/red]")
            return False
        except Exception as e:
            console.print(f"[red]✗ Error generating minimap2 index: {e}[/red]")
            return False

    def _generate_bwa_mem2_index(self, fasta_file: Path, species: str) -> bool:
        """Generate bwa-mem2 index."""
        try:
            console.print(f"Generating bwa-mem2 index for {species}...")

            # bwa-mem2 index input.fasta
            cmd = ["bwa-mem2", "index", str(fasta_file)]

            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=600,  # 10 minute timeout (bwa-mem2 can be slower)
            )

            if result.returncode == 0:
                # Check if index files were created
                expected_suffixes = [".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"]
                index_files = [
                    fasta_file.with_suffix(f"{fasta_file.suffix}{suffix}")
                    for suffix in expected_suffixes
                ]

                if all(f.exists() for f in index_files):
                    console.print(
                        f"[green]✓ Generated bwa-mem2 index files for {species}[/green]"
                    )
                    return True
                else:
                    console.print(
                        "[yellow]⚠ bwa-mem2 completed but some index files missing[/yellow]"
                    )
                    return False
            else:
                console.print(
                    f"[red]✗ Failed to generate bwa-mem2 index: {result.stderr}[/red]"
                )
                return False

        except subprocess.TimeoutExpired:
            console.print(f"[red]✗ bwa-mem2 indexing timed out for {species}[/red]")
            return False
        except Exception as e:
            console.print(f"[red]✗ Error generating bwa-mem2 index: {e}[/red]")
            return False
