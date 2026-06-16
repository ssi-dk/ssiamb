from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__
from .pipeline import PipelineError, run_self


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="ssiamb",
        description="Count ambiguous sites from reads mapped back to their assembly.",
    )
    parser.add_argument("--version", action="version", version=f"ssiamb {__version__}")

    subparsers = parser.add_subparsers(dest="command", required=True)
    self_parser = subparsers.add_parser(
        "self",
        help=(
            "Map trimmed paired reads to a supplied assembly and count ambiguous sites."
        ),
    )
    self_parser.add_argument(
        "--r1", required=True, type=Path, help="Forward reads FASTQ"
    )
    self_parser.add_argument(
        "--r2", required=True, type=Path, help="Reverse reads FASTQ"
    )
    self_parser.add_argument(
        "--assembly", required=True, type=Path, help="Assembly FASTA"
    )
    self_parser.add_argument("--sample", required=True, help="Sample name")
    self_parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("."),
        help="Output directory (default: current directory)",
    )
    self_parser.add_argument("--threads", type=int, default=1, help="Worker threads")
    self_parser.add_argument("--dp-min", type=int, default=10, help="Minimum depth")
    self_parser.add_argument(
        "--maf-min",
        type=float,
        default=0.10,
        help="Minimum minor allele frequency",
    )
    self_parser.add_argument(
        "--dp-cap",
        type=int,
        default=100,
        help="Depth cap for matrix/reporting",
    )
    self_parser.add_argument(
        "--bbtools-mem",
        help="Optional BBTools heap size without -Xmx prefix, for example 4g",
    )
    self_parser.add_argument(
        "--keep-intermediates",
        action="store_true",
        help="Keep sorted BAM intermediate",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "self":
        try:
            paths = run_self(
                r1=args.r1,
                r2=args.r2,
                assembly=args.assembly,
                sample=args.sample,
                outdir=args.outdir,
                threads=args.threads,
                dp_min=args.dp_min,
                maf_min=args.maf_min,
                dp_cap=args.dp_cap,
                keep_intermediates=args.keep_intermediates,
                bbtools_mem=args.bbtools_mem,
            )
        except PipelineError as exc:
            print(f"ssiamb: error: {exc}", file=sys.stderr)
            return 2
        print(paths.summary)
        return 0

    parser.error(f"unknown command: {args.command}")
    return 2
