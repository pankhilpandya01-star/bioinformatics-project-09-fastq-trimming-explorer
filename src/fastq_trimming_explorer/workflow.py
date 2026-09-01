"""Cutadapt orchestration and exact per-read accounting."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from dataclasses import asdict, dataclass
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

from fastq_trimming_explorer.analysis import AnalysisError, create_analysis_outputs
from fastq_trimming_explorer.fastq import (
    FastqRecord,
    FastqValidationError,
    read_fastq,
    records_by_identifier,
)


EXPECTED_CUTADAPT_VERSION = "5.2"
DEFAULT_ADAPTER = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPTER_ALPHABET = frozenset("ACGTURYSWKMBDHVN")


class WorkflowError(RuntimeError):
    """Raised when validation, Cutadapt, or accounting fails."""


@dataclass(frozen=True, slots=True)
class WorkflowConfig:
    fastq: Path
    adapter: str = DEFAULT_ADAPTER
    quality_cutoff: float = 20.0
    minimum_length: int = 75
    minimum_overlap: int = 8
    error_rate: float = 0.1
    output_dir: Path = Path("results/run")


@dataclass(frozen=True, slots=True)
class WorkflowResult:
    output_dir: Path
    input_reads: int
    retained_reads: int
    discarded_reads: int
    input_bases: int
    quality_bases_removed: int
    adapter_bases_removed: int

    @property
    def total_bases_removed(self) -> int:
        return self.quality_bases_removed + self.adapter_bases_removed


def _resolved(path: Path) -> Path:
    return path.expanduser().resolve()


def validate_config(config: WorkflowConfig) -> WorkflowConfig:
    fastq = _resolved(config.fastq)
    output_dir = _resolved(config.output_dir)
    adapter = config.adapter.strip().upper()

    if not fastq.exists():
        raise WorkflowError(f"input FASTQ does not exist: {fastq}")
    if not fastq.is_file():
        raise WorkflowError(f"input FASTQ is not a file: {fastq}")
    if output_dir.exists():
        raise WorkflowError(
            f"output directory already exists; choose a new path: {output_dir}"
        )
    if not adapter:
        raise WorkflowError("adapter must not be empty")
    invalid_adapter_symbols = sorted(set(adapter) - ADAPTER_ALPHABET)
    if invalid_adapter_symbols:
        symbols = "".join(invalid_adapter_symbols)
        raise WorkflowError(f"adapter contains unsupported symbols: {symbols}")
    if not 0 <= config.quality_cutoff <= 93:
        raise WorkflowError("quality cutoff must be between 0 and 93")
    if config.minimum_length < 1:
        raise WorkflowError("minimum length must be at least 1")
    if config.minimum_overlap < 1:
        raise WorkflowError("minimum overlap must be at least 1")
    if config.minimum_overlap > len(adapter):
        raise WorkflowError("minimum overlap cannot exceed adapter length")
    if not 0 <= config.error_rate < 1:
        raise WorkflowError("error rate must be at least 0 and less than 1")

    return WorkflowConfig(
        fastq=fastq,
        adapter=adapter,
        quality_cutoff=config.quality_cutoff,
        minimum_length=config.minimum_length,
        minimum_overlap=config.minimum_overlap,
        error_rate=config.error_rate,
        output_dir=output_dir,
    )


def _cutadapt_version() -> str:
    try:
        installed_version = version("cutadapt")
    except PackageNotFoundError as error:
        raise WorkflowError(
            "Cutadapt is not installed; install this project before running the workflow"
        ) from error
    if installed_version != EXPECTED_CUTADAPT_VERSION:
        raise WorkflowError(
            f"Cutadapt {EXPECTED_CUTADAPT_VERSION} is required; found {installed_version}"
        )
    return installed_version


def _format_number(value: float | int) -> str:
    return format(value, "g")


def _run_cutadapt(arguments: list[str], stage_name: str, staging_dir: Path) -> None:
    command = [sys.executable, "-m", "cutadapt", *arguments]
    completed = subprocess.run(command, capture_output=True, text=True, check=False)

    log_dir = staging_dir / "logs"
    log_dir.mkdir(exist_ok=True)
    (log_dir / f"{stage_name}.stdout.txt").write_text(
        completed.stdout, encoding="utf-8"
    )
    (log_dir / f"{stage_name}.stderr.txt").write_text(
        completed.stderr, encoding="utf-8"
    )
    if completed.returncode != 0:
        detail = completed.stderr.strip() or completed.stdout.strip() or "no diagnostic output"
        raise WorkflowError(
            f"Cutadapt stage {stage_name!r} failed with exit code "
            f"{completed.returncode}: {detail}"
        )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _load_json(path: Path) -> dict[str, Any]:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise WorkflowError(f"could not read Cutadapt JSON report {path.name}: {error}") from error


def _replace_path_fragments(value: Any, replacements: list[tuple[str, str]]) -> Any:
    """Replace run-specific absolute paths inside nested report values."""

    if isinstance(value, str):
        for original, portable in replacements:
            value = value.replace(original, portable)
        return value
    if isinstance(value, list):
        return [_replace_path_fragments(item, replacements) for item in value]
    if isinstance(value, dict):
        return {
            key: _replace_path_fragments(item, replacements)
            for key, item in value.items()
        }
    return value


def _make_reports_portable(
    staging_dir: Path,
    input_fastq: Path,
    report_paths: tuple[Path, ...],
) -> None:
    """Remove machine-specific paths from reports and captured tool logs."""

    replacements = sorted(
        [
            (str(staging_dir), "."),
            (str(input_fastq), input_fastq.name),
        ],
        key=lambda item: len(item[0]),
        reverse=True,
    )

    for report_path in report_paths:
        report = _replace_path_fragments(_load_json(report_path), replacements)
        report_path.write_text(
            json.dumps(report, indent=2) + "\n",
            encoding="utf-8",
        )

    for log_path in (staging_dir / "logs").glob("*.txt"):
        log_text = log_path.read_text(encoding="utf-8")
        for original, portable in replacements:
            log_text = log_text.replace(original, portable)
        log_path.write_text(log_text, encoding="utf-8")


def _assert_same_identifiers(
    expected: dict[str, FastqRecord],
    observed: dict[str, FastqRecord],
    stage_name: str,
) -> None:
    if expected.keys() != observed.keys():
        missing = sorted(expected.keys() - observed.keys())
        extra = sorted(observed.keys() - expected.keys())
        raise WorkflowError(
            f"identifier mismatch after {stage_name}: "
            f"missing={missing[:5]}, extra={extra[:5]}"
        )


def _adapter_trimmed_bases(adapter_report: dict[str, Any]) -> int:
    total = 0
    for adapter in adapter_report.get("adapters_read1") or []:
        end = adapter.get("three_prime_end") or adapter.get("five_prime_end") or {}
        for entry in end.get("trimmed_lengths") or []:
            total += int(entry["len"]) * sum(int(count) for count in entry["counts"])
    return total


def _validate_reports(
    original: list[FastqRecord],
    quality: list[FastqRecord],
    adapter: list[FastqRecord],
    retained: list[FastqRecord],
    discarded: list[FastqRecord],
    quality_report: dict[str, Any],
    adapter_report: dict[str, Any],
    length_report: dict[str, Any],
) -> None:
    original_bases = sum(record.length for record in original)
    quality_bases = sum(record.length for record in quality)
    adapter_bases = sum(record.length for record in adapter)
    retained_bases = sum(record.length for record in retained)
    discarded_bases = sum(record.length for record in discarded)

    checks = {
        "quality input reads": (quality_report["read_counts"]["input"], len(original)),
        "quality output reads": (quality_report["read_counts"]["output"], len(quality)),
        "quality input bases": (quality_report["basepair_counts"]["input"], original_bases),
        "quality output bases": (quality_report["basepair_counts"]["output"], quality_bases),
        "adapter input reads": (adapter_report["read_counts"]["input"], len(quality)),
        "adapter output reads": (adapter_report["read_counts"]["output"], len(adapter)),
        "adapter input bases": (adapter_report["basepair_counts"]["input"], quality_bases),
        "adapter output bases": (adapter_report["basepair_counts"]["output"], adapter_bases),
        "length input reads": (length_report["read_counts"]["input"], len(adapter)),
        "length retained reads": (length_report["read_counts"]["output"], len(retained)),
        "length discarded reads": (
            length_report["read_counts"]["filtered"]["too_short"] or 0,
            len(discarded),
        ),
        "length input bases": (length_report["basepair_counts"]["input"], adapter_bases),
        "length retained bases": (length_report["basepair_counts"]["output"], retained_bases),
    }
    for label, (reported, observed) in checks.items():
        if reported != observed:
            raise WorkflowError(
                f"accounting mismatch for {label}: report={reported}, observed={observed}"
            )

    if adapter_bases != retained_bases + discarded_bases:
        raise WorkflowError(
            "base accounting mismatch: adapter output does not equal retained plus discarded"
        )
    report_quality_removed = quality_report["basepair_counts"]["quality_trimmed"]
    if report_quality_removed != original_bases - quality_bases:
        raise WorkflowError("quality-trimmed base count disagrees with FASTQ lengths")
    if _adapter_trimmed_bases(adapter_report) != quality_bases - adapter_bases:
        raise WorkflowError("adapter-trimmed base count disagrees with FASTQ lengths")


def _write_per_read_csv(
    path: Path,
    original: list[FastqRecord],
    quality: dict[str, FastqRecord],
    adapter: dict[str, FastqRecord],
    retained: dict[str, FastqRecord],
    discarded: dict[str, FastqRecord],
) -> None:
    fieldnames = [
        "read_id",
        "original_length",
        "post_quality_length",
        "post_adapter_length",
        "final_length",
        "quality_bases_removed",
        "adapter_bases_removed",
        "total_bases_removed",
        "status",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for original_record in original:
            read_id = original_record.identifier
            quality_length = quality[read_id].length
            adapter_length = adapter[read_id].length
            if read_id in retained:
                final_length = retained[read_id].length
                status = "retained"
            else:
                final_length = discarded[read_id].length
                status = "discarded_too_short"
            quality_removed = original_record.length - quality_length
            adapter_removed = quality_length - adapter_length
            writer.writerow(
                {
                    "read_id": read_id,
                    "original_length": original_record.length,
                    "post_quality_length": quality_length,
                    "post_adapter_length": adapter_length,
                    "final_length": final_length,
                    "quality_bases_removed": quality_removed,
                    "adapter_bases_removed": adapter_removed,
                    "total_bases_removed": quality_removed + adapter_removed,
                    "status": status,
                }
            )


def run_workflow(config: WorkflowConfig) -> WorkflowResult:
    """Run the complete workflow and atomically publish a verified output directory."""

    config = validate_config(config)
    cutadapt_version = _cutadapt_version()
    try:
        original_records = read_fastq(config.fastq)
    except FastqValidationError as error:
        raise WorkflowError(str(error)) from error

    input_checksum = _sha256(config.fastq)
    config.output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging_dir = Path(
        tempfile.mkdtemp(
            prefix=f".project09-{config.output_dir.name}-",
            dir=config.output_dir.parent,
        )
    )

    quality_fastq = staging_dir / "01_quality_trimmed.fastq"
    adapter_fastq = staging_dir / "02_adapter_trimmed.fastq"
    retained_fastq = staging_dir / "retained_reads.fastq"
    discarded_fastq = staging_dir / "discarded_too_short.fastq"
    reports_dir = staging_dir / "reports"
    reports_dir.mkdir()
    quality_json = reports_dir / "01_quality.cutadapt.json"
    adapter_json = reports_dir / "02_adapter.cutadapt.json"
    length_json = reports_dir / "03_length_filter.cutadapt.json"

    try:
        _run_cutadapt(
            [
                "--cores",
                "1",
                "--quality-cutoff",
                _format_number(config.quality_cutoff),
                "--json",
                str(quality_json),
                "--output",
                str(quality_fastq),
                str(config.fastq),
            ],
            "01_quality",
            staging_dir,
        )
        _run_cutadapt(
            [
                "--cores",
                "1",
                "--adapter",
                config.adapter,
                "--overlap",
                str(config.minimum_overlap),
                "--error-rate",
                _format_number(config.error_rate),
                "--json",
                str(adapter_json),
                "--output",
                str(adapter_fastq),
                str(quality_fastq),
            ],
            "02_adapter",
            staging_dir,
        )
        _run_cutadapt(
            [
                "--cores",
                "1",
                "--minimum-length",
                str(config.minimum_length),
                "--too-short-output",
                str(discarded_fastq),
                "--json",
                str(length_json),
                "--output",
                str(retained_fastq),
                str(adapter_fastq),
            ],
            "03_length_filter",
            staging_dir,
        )

        quality_records = read_fastq(quality_fastq, allow_empty_sequences=True)
        adapter_records = read_fastq(adapter_fastq, allow_empty_sequences=True)
        retained_records = read_fastq(retained_fastq, allow_empty=True)
        discarded_records = read_fastq(
            discarded_fastq,
            allow_empty=True,
            allow_empty_sequences=True,
        )

        original_by_id = records_by_identifier(original_records)
        quality_by_id = records_by_identifier(quality_records)
        adapter_by_id = records_by_identifier(adapter_records)
        retained_by_id = records_by_identifier(retained_records)
        discarded_by_id = records_by_identifier(discarded_records)
        _assert_same_identifiers(original_by_id, quality_by_id, "quality trimming")
        _assert_same_identifiers(original_by_id, adapter_by_id, "adapter trimming")
        if retained_by_id.keys() & discarded_by_id.keys():
            raise WorkflowError("retained and discarded outputs contain overlapping identifiers")
        partition = {**retained_by_id, **discarded_by_id}
        _assert_same_identifiers(original_by_id, partition, "length filtering")

        quality_report = _load_json(quality_json)
        adapter_report = _load_json(adapter_json)
        length_report = _load_json(length_json)
        _validate_reports(
            original_records,
            quality_records,
            adapter_records,
            retained_records,
            discarded_records,
            quality_report,
            adapter_report,
            length_report,
        )
        _make_reports_portable(
            staging_dir,
            config.fastq,
            (quality_json, adapter_json, length_json),
        )

        per_read_csv = staging_dir / "per_read_trimming.csv"
        _write_per_read_csv(
            per_read_csv,
            original_records,
            quality_by_id,
            adapter_by_id,
            retained_by_id,
            discarded_by_id,
        )

        analysis_result = create_analysis_outputs(
            original_records=original_records,
            retained_records=retained_records,
            discarded_records=discarded_records,
            per_read_csv=per_read_csv,
            output_dir=staging_dir,
            input_label=config.fastq.name,
        )

        quality_bases_removed = sum(
            record.length - quality_by_id[record.identifier].length
            for record in original_records
        )
        adapter_bases_removed = sum(
            quality_by_id[record.identifier].length
            - adapter_by_id[record.identifier].length
            for record in original_records
        )
        result = WorkflowResult(
            output_dir=config.output_dir,
            input_reads=len(original_records),
            retained_reads=len(retained_records),
            discarded_reads=len(discarded_records),
            input_bases=sum(record.length for record in original_records),
            quality_bases_removed=quality_bases_removed,
            adapter_bases_removed=adapter_bases_removed,
        )
        manifest = {
            "workflow_version": "0.1.0",
            "cutadapt_version": cutadapt_version,
            "input": {
                "path": config.fastq.name,
                "sha256": input_checksum,
            },
            "parameters": {
                "adapter": config.adapter,
                "quality_cutoff": config.quality_cutoff,
                "minimum_length": config.minimum_length,
                "minimum_overlap": config.minimum_overlap,
                "error_rate": config.error_rate,
            },
            "counts": {
                **asdict(result),
                "output_dir": ".",
                "total_bases_removed": result.total_bases_removed,
            },
            "analysis": asdict(analysis_result),
        }
        (staging_dir / "run_manifest.json").write_text(
            json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
        )

        if _sha256(config.fastq) != input_checksum:
            raise WorkflowError("input FASTQ changed while the workflow was running")
        os.replace(staging_dir, config.output_dir)
        return result
    except (AnalysisError, FastqValidationError, KeyError, OSError) as error:
        raise WorkflowError(f"workflow output validation failed: {error}") from error
    finally:
        if staging_dir.exists():
            shutil.rmtree(staging_dir)
