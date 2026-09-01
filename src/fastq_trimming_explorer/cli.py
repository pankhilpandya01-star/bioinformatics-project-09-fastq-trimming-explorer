"""Command-line interface for the FASTQ Trimming Explorer."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from fastq_trimming_explorer import __version__
from fastq_trimming_explorer.workflow import (
    DEFAULT_ADAPTER,
    WorkflowConfig,
    WorkflowError,
    run_workflow,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="fastq-trimming-explorer",
        description=(
            "Run traceable single-end quality, adapter, and minimum-length trimming."
        ),
    )
    parser.add_argument("--fastq", required=True, type=Path, help="input FASTQ or FASTQ.GZ")
    parser.add_argument("--adapter", default=DEFAULT_ADAPTER, help="3-prime adapter sequence")
    parser.add_argument(
        "--quality-cutoff",
        type=float,
        default=20.0,
        help="3-prime Cutadapt quality cutoff (default: 20)",
    )
    parser.add_argument(
        "--minimum-length",
        type=int,
        default=75,
        help="minimum retained read length (default: 75)",
    )
    parser.add_argument(
        "--minimum-overlap",
        type=int,
        default=8,
        help="minimum adapter overlap (default: 8)",
    )
    parser.add_argument(
        "--error-rate",
        type=float,
        default=0.1,
        help="maximum adapter error rate (default: 0.1)",
    )
    parser.add_argument("--output-dir", required=True, type=Path, help="new output directory")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    config = WorkflowConfig(
        fastq=args.fastq,
        adapter=args.adapter,
        quality_cutoff=args.quality_cutoff,
        minimum_length=args.minimum_length,
        minimum_overlap=args.minimum_overlap,
        error_rate=args.error_rate,
        output_dir=args.output_dir,
    )
    try:
        result = run_workflow(config)
    except WorkflowError as error:
        print(f"error: {error}", file=sys.stderr)
        return 2

    print(f"Output directory: {result.output_dir}")
    print(f"Input reads: {result.input_reads}")
    print(f"Retained reads: {result.retained_reads}")
    print(f"Discarded reads: {result.discarded_reads}")
    print(f"Quality-trimmed bases: {result.quality_bases_removed}")
    print(f"Adapter-trimmed bases: {result.adapter_bases_removed}")
    print(f"Total trimmed bases: {result.total_bases_removed}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
