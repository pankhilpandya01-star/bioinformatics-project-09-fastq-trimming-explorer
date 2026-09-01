"""Before/after tables and dashboard generation."""

from __future__ import annotations

import csv
import os
import statistics
import tempfile
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(tempfile.gettempdir()) / "fastq-trimming-explorer-matplotlib"),
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from fastq_trimming_explorer.fastq import FastqRecord


class AnalysisError(RuntimeError):
    """Raised when analysis artifacts cannot be generated consistently."""


@dataclass(frozen=True, slots=True)
class AnalysisResult:
    before_mean_quality: float
    after_mean_quality: float
    before_q20_percent: float
    after_q20_percent: float
    before_q30_percent: float
    after_q30_percent: float
    discarded_bases: int
    quality_changed_reads: int
    adapter_changed_reads: int
    any_changed_reads: int


@dataclass(slots=True)
class _DatasetStatistics:
    read_count: int
    total_bases: int
    lengths: list[int]
    mean_length: float
    median_length: float
    minimum_length: int
    maximum_length: int
    mean_quality: float
    q20_percent: float
    q30_percent: float
    positional_rows: list[dict[str, int | float]]


def _phred_scores(record: FastqRecord) -> list[int]:
    return [ord(character) - 33 for character in record.quality]


def _dataset_statistics(records: list[FastqRecord]) -> _DatasetStatistics:
    lengths = [record.length for record in records]
    total_bases = sum(lengths)
    all_quality_sum = 0
    q20_bases = 0
    q30_bases = 0
    maximum_length = max(lengths, default=0)
    position_scores: list[list[int]] = [[] for _ in range(maximum_length)]

    for record in records:
        scores = _phred_scores(record)
        all_quality_sum += sum(scores)
        q20_bases += sum(score >= 20 for score in scores)
        q30_bases += sum(score >= 30 for score in scores)
        for index, score in enumerate(scores):
            position_scores[index].append(score)

    positional_rows: list[dict[str, int | float]] = []
    for index, scores in enumerate(position_scores, start=1):
        positional_rows.append(
            {
                "position": index,
                "reads_at_position": len(scores),
                "mean_quality": statistics.fmean(scores),
                "median_quality": statistics.median(scores),
                "q20_percent": 100 * sum(score >= 20 for score in scores) / len(scores),
                "q30_percent": 100 * sum(score >= 30 for score in scores) / len(scores),
            }
        )

    return _DatasetStatistics(
        read_count=len(records),
        total_bases=total_bases,
        lengths=lengths,
        mean_length=statistics.fmean(lengths) if lengths else 0.0,
        median_length=statistics.median(lengths) if lengths else 0.0,
        minimum_length=min(lengths, default=0),
        maximum_length=maximum_length,
        mean_quality=all_quality_sum / total_bases if total_bases else 0.0,
        q20_percent=100 * q20_bases / total_bases if total_bases else 0.0,
        q30_percent=100 * q30_bases / total_bases if total_bases else 0.0,
        positional_rows=positional_rows,
    )


def _write_positional_csv(path: Path, rows: list[dict[str, int | float]]) -> None:
    fieldnames = [
        "position",
        "reads_at_position",
        "mean_quality",
        "median_quality",
        "q20_percent",
        "q30_percent",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    **row,
                    "mean_quality": f"{float(row['mean_quality']):.4f}",
                    "median_quality": f"{float(row['median_quality']):.4f}",
                    "q20_percent": f"{float(row['q20_percent']):.4f}",
                    "q30_percent": f"{float(row['q30_percent']):.4f}",
                }
            )


def _relative_change(before: float | int, after: float | int) -> str:
    if before == 0:
        return ""
    return f"{100 * (after - before) / before:.6f}"


def _write_before_after_summary(
    path: Path,
    before: _DatasetStatistics,
    after: _DatasetStatistics,
) -> None:
    metrics: list[tuple[str, float | int, float | int, str]] = [
        ("read_count", before.read_count, after.read_count, "reads"),
        ("total_bases", before.total_bases, after.total_bases, "bases"),
        ("mean_read_length", before.mean_length, after.mean_length, "nt"),
        ("median_read_length", before.median_length, after.median_length, "nt"),
        ("minimum_read_length", before.minimum_length, after.minimum_length, "nt"),
        ("maximum_read_length", before.maximum_length, after.maximum_length, "nt"),
        ("mean_base_quality", before.mean_quality, after.mean_quality, "Phred"),
        ("q20_bases", before.q20_percent, after.q20_percent, "percent"),
        ("q30_bases", before.q30_percent, after.q30_percent, "percent"),
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            ["metric", "before", "after", "absolute_change", "relative_change_percent", "unit"]
        )
        for metric, before_value, after_value, unit in metrics:
            writer.writerow(
                [
                    metric,
                    f"{before_value:.6f}" if isinstance(before_value, float) else before_value,
                    f"{after_value:.6f}" if isinstance(after_value, float) else after_value,
                    f"{after_value - before_value:.6f}",
                    _relative_change(before_value, after_value),
                    unit,
                ]
            )


def _read_per_read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _write_trimming_summary(
    path: Path,
    *,
    before: _DatasetStatistics,
    after: _DatasetStatistics,
    discarded: _DatasetStatistics,
    per_read_rows: list[dict[str, str]],
) -> dict[str, int]:
    quality_changed = sum(int(row["quality_bases_removed"]) > 0 for row in per_read_rows)
    adapter_changed = sum(int(row["adapter_bases_removed"]) > 0 for row in per_read_rows)
    any_changed = sum(int(row["total_bases_removed"]) > 0 for row in per_read_rows)
    quality_removed = sum(int(row["quality_bases_removed"]) for row in per_read_rows)
    adapter_removed = sum(int(row["adapter_bases_removed"]) for row in per_read_rows)
    metrics: list[tuple[str, int | float, str]] = [
        ("input_reads", before.read_count, "reads"),
        ("retained_reads", after.read_count, "reads"),
        ("discarded_too_short_reads", discarded.read_count, "reads"),
        ("quality_changed_reads", quality_changed, "reads"),
        ("adapter_changed_reads", adapter_changed, "reads"),
        ("any_trimmed_reads", any_changed, "reads"),
        ("input_bases", before.total_bases, "bases"),
        ("quality_trimmed_bases", quality_removed, "bases"),
        ("adapter_trimmed_bases", adapter_removed, "bases"),
        ("total_trimmed_bases", quality_removed + adapter_removed, "bases"),
        ("retained_bases", after.total_bases, "bases"),
        ("bases_in_discarded_reads", discarded.total_bases, "bases"),
        ("trimmed_bases_percent", 100 * (quality_removed + adapter_removed) / before.total_bases, "percent"),
        ("retained_bases_percent", 100 * after.total_bases / before.total_bases, "percent"),
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["metric", "value", "unit"])
        for metric, value, unit in metrics:
            rendered = f"{value:.6f}" if isinstance(value, float) else value
            writer.writerow([metric, rendered, unit])

    return {
        "quality_changed_reads": quality_changed,
        "adapter_changed_reads": adapter_changed,
        "any_changed_reads": any_changed,
        "quality_removed_bases": quality_removed,
        "adapter_removed_bases": adapter_removed,
    }


def _annotate_bars(axis: Any) -> None:
    for patch in axis.patches:
        height = patch.get_height()
        axis.annotate(
            f"{int(height):,}",
            (patch.get_x() + patch.get_width() / 2, height),
            ha="center",
            va="bottom",
            xytext=(0, 5),
            textcoords="offset points",
            fontsize=9,
            fontweight="bold",
        )


def _create_dashboard(
    path: Path,
    *,
    before: _DatasetStatistics,
    after: _DatasetStatistics,
    discarded: _DatasetStatistics,
    trimming: dict[str, int],
    input_label: str,
) -> None:
    plt.style.use("seaborn-v0_8-whitegrid")
    navy = "#17324D"
    teal = "#00A6A6"
    orange = "#F28E2B"
    red = "#D1495B"
    gray = "#6B7280"
    figure, axes = plt.subplots(2, 2, figsize=(14, 10), constrained_layout=True)
    figure.patch.set_facecolor("#F7F9FC")
    figure.suptitle(
        "Project 09 · FASTQ Trimming Explorer",
        fontsize=22,
        fontweight="bold",
        color=navy,
    )
    figure.text(
        0.5,
        0.955,
        f"{input_label} · Cutadapt 5.2 · Q20 · minimum length 75 nt",
        ha="center",
        fontsize=11,
        color=gray,
    )

    disposition = axes[0, 0]
    disposition.bar(
        ["Retained", "Too short"],
        [after.read_count, discarded.read_count],
        color=[teal, red],
        width=0.62,
    )
    disposition.set_title("A · Read disposition", loc="left", fontweight="bold", color=navy)
    disposition.set_ylabel("Reads")
    _annotate_bars(disposition)

    lengths_axis = axes[0, 1]
    before_lengths = Counter(before.lengths)
    after_lengths = Counter(after.lengths)
    lengths_axis.plot(
        sorted(before_lengths),
        [before_lengths[length] for length in sorted(before_lengths)],
        marker="o",
        linewidth=2.2,
        drawstyle="steps-mid",
        color=navy,
        label="Before: all raw reads",
    )
    if after_lengths:
        lengths_axis.plot(
            sorted(after_lengths),
            [after_lengths[length] for length in sorted(after_lengths)],
            marker="o",
            linewidth=2.2,
            drawstyle="steps-mid",
            color=teal,
            label="After: retained reads",
        )
    lengths_axis.set_title("B · Read-length distribution", loc="left", fontweight="bold", color=navy)
    lengths_axis.set_xlabel("Read length (nt)")
    lengths_axis.set_ylabel("Reads (log scale)")
    lengths_axis.set_yscale("log")
    lengths_axis.legend(frameon=False)

    quality_axis = axes[1, 0]
    quality_axis.plot(
        [row["position"] for row in before.positional_rows],
        [row["mean_quality"] for row in before.positional_rows],
        color=navy,
        linewidth=2.2,
        label="Before: all raw reads",
    )
    if after.positional_rows:
        quality_axis.plot(
            [row["position"] for row in after.positional_rows],
            [row["mean_quality"] for row in after.positional_rows],
            color=teal,
            linewidth=2.2,
            label="After: retained reads",
        )
    quality_axis.axhline(20, color=orange, linestyle="--", linewidth=1.5, label="Q20")
    quality_axis.set_title("C · Mean quality by position", loc="left", fontweight="bold", color=navy)
    quality_axis.set_xlabel("Position (nt)")
    quality_axis.set_ylabel("Mean Phred score")
    quality_axis.legend(frameon=False)

    consequences = axes[1, 1]
    consequence_labels = ["Quality\ntrimmed", "Adapter\ntrimmed", "Routed with\ntoo-short reads"]
    consequence_values = [
        trimming["quality_removed_bases"],
        trimming["adapter_removed_bases"],
        discarded.total_bases,
    ]
    consequences.bar(
        consequence_labels,
        consequence_values,
        color=[orange, navy, red],
        width=0.62,
    )
    consequences.set_title("D · Base-level consequences", loc="left", fontweight="bold", color=navy)
    consequences.set_ylabel("Bases")
    _annotate_bars(consequences)
    consequences.text(
        0.02,
        0.98,
        "Too-short bases are preserved in a separate FASTQ.",
        transform=consequences.transAxes,
        va="top",
        fontsize=9,
        color=gray,
    )

    for axis in axes.flat:
        axis.set_facecolor("white")
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)

    figure.savefig(path, dpi=180, facecolor=figure.get_facecolor())
    plt.close(figure)


def create_analysis_outputs(
    *,
    original_records: list[FastqRecord],
    retained_records: list[FastqRecord],
    discarded_records: list[FastqRecord],
    per_read_csv: Path,
    output_dir: Path,
    input_label: str,
) -> AnalysisResult:
    """Write verified before/after CSVs and a four-panel PNG dashboard."""

    try:
        before = _dataset_statistics(original_records)
        after = _dataset_statistics(retained_records)
        discarded = _dataset_statistics(discarded_records)
        per_read_rows = _read_per_read_rows(per_read_csv)
        if len(per_read_rows) != before.read_count:
            raise AnalysisError(
                "per-read table row count does not match the input FASTQ"
            )
        if after.read_count + discarded.read_count != before.read_count:
            raise AnalysisError("retained and discarded read counts do not partition the input")

        _write_before_after_summary(
            output_dir / "before_after_summary.csv",
            before,
            after,
        )
        _write_positional_csv(
            output_dir / "before_positional_quality.csv",
            before.positional_rows,
        )
        _write_positional_csv(
            output_dir / "after_positional_quality.csv",
            after.positional_rows,
        )
        trimming = _write_trimming_summary(
            output_dir / "trimming_summary.csv",
            before=before,
            after=after,
            discarded=discarded,
            per_read_rows=per_read_rows,
        )
        expected_retained_bases = (
            before.total_bases
            - trimming["quality_removed_bases"]
            - trimming["adapter_removed_bases"]
            - discarded.total_bases
        )
        if expected_retained_bases != after.total_bases:
            raise AnalysisError(
                "analysis base totals do not reconcile with retained and discarded FASTQs"
            )
        _create_dashboard(
            output_dir / "comparison_dashboard.png",
            before=before,
            after=after,
            discarded=discarded,
            trimming=trimming,
            input_label=input_label,
        )
        return AnalysisResult(
            before_mean_quality=before.mean_quality,
            after_mean_quality=after.mean_quality,
            before_q20_percent=before.q20_percent,
            after_q20_percent=after.q20_percent,
            before_q30_percent=before.q30_percent,
            after_q30_percent=after.q30_percent,
            discarded_bases=discarded.total_bases,
            quality_changed_reads=trimming["quality_changed_reads"],
            adapter_changed_reads=trimming["adapter_changed_reads"],
            any_changed_reads=trimming["any_changed_reads"],
        )
    except AnalysisError:
        raise
    except Exception as error:
        raise AnalysisError(f"could not generate analysis outputs: {error}") from error
