from __future__ import annotations

import os
import re
import shutil
import subprocess
import time
from dataclasses import dataclass
from os import PathLike
from pathlib import Path

from .ambiguity import summarize_records
from .reports import write_loci, write_matrix, write_provenance, write_summary
from .vcf import VcfParseError, read_vcf


class PipelineError(RuntimeError):
    """Raised when the external self workflow fails."""


@dataclass(frozen=True)
class OutputPaths:
    summary: Path
    loci: Path
    matrix: Path
    vcf: Path
    provenance: Path
    bam: Path


def safe_sample_name(sample: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]", "_", sample.strip())
    return cleaned or "sample"


def require_tool(name: str) -> str:
    path = shutil.which(name)
    if path is None:
        raise PipelineError(f"Required tool not found on PATH: {name}")
    return path


def tool_version(command: list[str]) -> str:
    try:
        result = subprocess.run(
            command,
            check=False,
            capture_output=True,
            text=True,
            timeout=10,
        )
    except Exception as exc:  # pragma: no cover - defensive provenance path
        return f"version unavailable: {exc}"
    output = (result.stdout or result.stderr).strip().splitlines()
    return output[0] if output else "version unavailable"


def env_with_conda_java(tool_path: str | PathLike[str]) -> dict[str, str]:
    env = os.environ.copy()
    prefix = Path(tool_path).parent.parent
    conda_java_home = prefix / "lib" / "jvm"
    conda_java = conda_java_home / "bin" / "java"
    if conda_java.exists():
        env["JAVA_HOME"] = str(conda_java_home)
        env["PATH"] = f"{conda_java_home / 'bin'}:{env.get('PATH', '')}"
    return env


def output_paths(outdir: Path, sample: str) -> OutputPaths:
    prefix = outdir / safe_sample_name(sample)
    return OutputPaths(
        summary=Path(f"{prefix}.summary.tsv"),
        loci=Path(f"{prefix}.loci.tsv"),
        matrix=Path(f"{prefix}.matrix.tsv"),
        vcf=Path(f"{prefix}.variants.vcf"),
        provenance=Path(f"{prefix}.provenance.json"),
        bam=Path(f"{prefix}.sorted.bam"),
    )


def run_mapping(
    *,
    minimap2: str,
    samtools: str,
    r1: Path,
    r2: Path,
    assembly: Path,
    bam: Path,
    threads: int,
) -> None:
    minimap_cmd = [
        minimap2,
        "-t",
        str(threads),
        "--MD",
        "-ax",
        "sr",
        str(assembly),
        str(r1),
        str(r2),
    ]
    sort_cmd = [samtools, "sort", "-@", str(threads), "-o", str(bam), "-"]
    minimap_stderr_path = bam.with_suffix(".minimap2.stderr.log")

    with (
        minimap_stderr_path.open("wb") as minimap_stderr_handle,
        subprocess.Popen(
            minimap_cmd,
            stdout=subprocess.PIPE,
            stderr=minimap_stderr_handle,
            text=False,
        ) as minimap_proc,
    ):
        sort_proc = subprocess.run(
            sort_cmd,
            stdin=minimap_proc.stdout,
            capture_output=True,
            text=True,
            check=False,
        )
        if minimap_proc.stdout is not None:
            minimap_proc.stdout.close()
        minimap_return = minimap_proc.wait()

    minimap_stderr = minimap_stderr_path.read_text(
        encoding="utf-8",
        errors="replace",
    )
    minimap_stderr_path.unlink(missing_ok=True)

    if minimap_return != 0:
        raise PipelineError(
            f"minimap2 failed with exit {minimap_return}: {minimap_stderr}"
        )
    if sort_proc.returncode != 0:
        raise PipelineError(
            f"samtools sort failed with exit {sort_proc.returncode}: {sort_proc.stderr}"
        )


def run_callvariants(
    *,
    callvariants: str,
    bam: Path,
    assembly: Path,
    vcf: Path,
    threads: int,
    bbtools_mem: str | None,
) -> None:
    cmd = [callvariants]
    if bbtools_mem:
        cmd.append(f"-Xmx{bbtools_mem}")
    cmd.extend(
        [
            f"in={bam}",
            f"ref={assembly}",
            f"vcf={vcf}",
            "ploidy=1",
            "clearfilters=t",
            f"threads={threads}",
        ]
    )
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        check=False,
        env=env_with_conda_java(callvariants),
    )
    if result.returncode != 0:
        raise PipelineError(
            f"callvariants.sh failed with exit {result.returncode}: {result.stderr}"
        )


def validate_input(path: Path, label: str) -> None:
    if not path.exists():
        raise PipelineError(f"{label} does not exist: {path}")
    if not path.is_file():
        raise PipelineError(f"{label} is not a file: {path}")


def run_self(
    *,
    r1: Path,
    r2: Path,
    assembly: Path,
    sample: str,
    outdir: Path,
    threads: int,
    dp_min: int,
    maf_min: float,
    dp_cap: int,
    keep_intermediates: bool,
    bbtools_mem: str | None = None,
) -> OutputPaths:
    validate_input(r1, "R1 reads")
    validate_input(r2, "R2 reads")
    validate_input(assembly, "Assembly")
    if threads < 1:
        raise PipelineError("--threads must be >= 1")
    if dp_min < 1:
        raise PipelineError("--dp-min must be >= 1")
    if not 0 <= maf_min <= 0.5:
        raise PipelineError("--maf-min must be between 0 and 0.5")
    if dp_cap < 1:
        raise PipelineError("--dp-cap must be >= 1")

    outdir.mkdir(parents=True, exist_ok=True)
    paths = output_paths(outdir, sample)

    minimap2 = require_tool("minimap2")
    samtools = require_tool("samtools")
    callvariants = require_tool("callvariants.sh")

    started = time.time()
    run_mapping(
        minimap2=minimap2,
        samtools=samtools,
        r1=r1,
        r2=r2,
        assembly=assembly,
        bam=paths.bam,
        threads=threads,
    )
    run_callvariants(
        callvariants=callvariants,
        bam=paths.bam,
        assembly=assembly,
        vcf=paths.vcf,
        threads=threads,
        bbtools_mem=bbtools_mem,
    )

    try:
        result = summarize_records(
            read_vcf(paths.vcf),
            sample=sample,
            dp_min=dp_min,
            maf_min=maf_min,
            dp_cap=dp_cap,
        )
    except VcfParseError as exc:
        raise PipelineError(f"Failed to parse callvariants VCF: {exc}") from exc
    write_summary(paths.summary, result, r1=r1, r2=r2, assembly=assembly)
    write_loci(paths.loci, result)
    write_matrix(paths.matrix, result)
    write_provenance(
        paths.provenance,
        {
            "sample": sample,
            "safe_sample": safe_sample_name(sample),
            "inputs": {"r1": str(r1), "r2": str(r2), "assembly": str(assembly)},
            "thresholds": {"dp_min": dp_min, "maf_min": maf_min, "dp_cap": dp_cap},
            "outputs": {
                "summary": str(paths.summary),
                "loci": str(paths.loci),
                "matrix": str(paths.matrix),
                "vcf": str(paths.vcf),
                "bam": str(paths.bam) if keep_intermediates else "",
            },
            "tools": {
                "minimap2": {
                    "path": minimap2,
                    "version": tool_version([minimap2, "--version"]),
                },
                "samtools": {
                    "path": samtools,
                    "version": tool_version([samtools, "--version"]),
                },
                "callvariants.sh": {"path": callvariants, "version": "BBTools script"},
            },
            "runtime_seconds": round(time.time() - started, 3),
        },
    )

    if not keep_intermediates and paths.bam.exists():
        paths.bam.unlink()

    return paths
