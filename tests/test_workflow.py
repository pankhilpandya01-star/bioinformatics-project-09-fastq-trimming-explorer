from __future__ import annotations

import csv
import hashlib
import io
import json
import subprocess
import tempfile
import unittest
from contextlib import redirect_stderr
from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

from fastq_trimming_explorer.cli import main as cli_main
from fastq_trimming_explorer.fastq import read_fastq
from fastq_trimming_explorer.workflow import (
    DEFAULT_ADAPTER,
    WorkflowConfig,
    WorkflowError,
    run_workflow,
    validate_config,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
CONTROLLED_FASTQ = PROJECT_ROOT / "tests" / "fixtures" / "controlled.fastq"
PUBLIC_FASTQ = PROJECT_ROOT / "data" / "SRR13921545_read1_first5000.fastq"
PUBLIC_FASTQ_SHA256 = (
    "c3d708ba4b1c7eb4bb95dbae6f7189f291250d129514f4d26b50aad272fecc15"
)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def fastq_ids(path: Path, *, allow_empty: bool = False) -> list[str]:
    return [record.identifier for record in read_fastq(path, allow_empty=allow_empty)]


class WorkflowValidationTests(unittest.TestCase):
    def test_rejects_invalid_configuration_values(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            fastq = root / "valid.fastq"
            fastq.write_text("@read1\nACGT\n+\nIIII\n", encoding="ascii")
            base = WorkflowConfig(fastq=fastq, output_dir=root / "output")
            cases = {
                "empty_adapter": replace(base, adapter=""),
                "invalid_adapter": replace(base, adapter="ACGT*"),
                "negative_quality": replace(base, quality_cutoff=-1),
                "quality_above_phred_range": replace(base, quality_cutoff=94),
                "zero_minimum_length": replace(base, minimum_length=0),
                "zero_overlap": replace(base, minimum_overlap=0),
                "overlap_above_adapter_length": replace(base, minimum_overlap=34),
                "negative_error_rate": replace(base, error_rate=-0.01),
                "error_rate_one_is_an_absolute_count_in_cutadapt": replace(base, error_rate=1),
                "error_rate_above_one": replace(base, error_rate=1.01),
                "missing_input": replace(base, fastq=root / "missing.fastq"),
            }
            existing_output = root / "existing"
            existing_output.mkdir()
            cases["existing_output"] = replace(base, output_dir=existing_output)

            for name, config in cases.items():
                with self.subTest(case=name):
                    with self.assertRaises(WorkflowError):
                        validate_config(config)

    def test_cli_returns_two_and_creates_nothing_for_invalid_value(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            output = Path(temporary_directory) / "must_not_exist"
            error_stream = io.StringIO()
            with redirect_stderr(error_stream):
                exit_code = cli_main(
                    [
                        "--fastq",
                        str(CONTROLLED_FASTQ),
                        "--quality-cutoff",
                        "-1",
                        "--output-dir",
                        str(output),
                    ]
                )
            self.assertEqual(2, exit_code)
            self.assertIn("quality cutoff", error_stream.getvalue())
            self.assertFalse(output.exists())

    def test_subprocess_failure_leaves_no_partial_output(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            fastq = root / "input.fastq"
            fastq.write_text("@read1\n" + "A" * 100 + "\n+\n" + "I" * 100 + "\n")
            output = root / "failed_output"
            failed_process = subprocess.CompletedProcess(
                args=["cutadapt"],
                returncode=17,
                stdout="",
                stderr="controlled subprocess failure",
            )
            with patch(
                "fastq_trimming_explorer.workflow.subprocess.run",
                return_value=failed_process,
            ):
                with self.assertRaisesRegex(WorkflowError, "controlled subprocess failure"):
                    run_workflow(WorkflowConfig(fastq=fastq, output_dir=output))
            self.assertFalse(output.exists())
            self.assertEqual([], list(root.glob(".project09-*")))


class WorkflowIntegrationTests(unittest.TestCase):
    def test_known_quality_adapter_variable_length_and_exact_accounting(self) -> None:
        expected = {
            "unchanged": (100, 100, 0, 0, "retained"),
            "quality_tail": (100, 90, 10, 0, "retained"),
            "adapter_retained": (113, 80, 0, 33, "retained"),
            "adapter_discarded": (103, 70, 0, 33, "discarded_too_short"),
            "quality_discarded": (80, 70, 10, 0, "discarded_too_short"),
        }
        checksum_before = hashlib.sha256(CONTROLLED_FASTQ.read_bytes()).hexdigest()
        with tempfile.TemporaryDirectory() as temporary_directory:
            output = Path(temporary_directory) / "controlled"
            result = run_workflow(
                WorkflowConfig(fastq=CONTROLLED_FASTQ, output_dir=output)
            )
            rows = read_csv(output / "per_read_trimming.csv")
            rows_by_id = {row["read_id"]: row for row in rows}

            self.assertEqual(5, result.input_reads)
            self.assertEqual(3, result.retained_reads)
            self.assertEqual(2, result.discarded_reads)
            self.assertEqual(20, result.quality_bases_removed)
            self.assertEqual(66, result.adapter_bases_removed)
            self.assertEqual(expected.keys(), rows_by_id.keys())
            for read_id, (
                original_length,
                final_length,
                quality_removed,
                adapter_removed,
                status,
            ) in expected.items():
                row = rows_by_id[read_id]
                self.assertEqual(original_length, int(row["original_length"]))
                self.assertEqual(final_length, int(row["final_length"]))
                self.assertEqual(quality_removed, int(row["quality_bases_removed"]))
                self.assertEqual(adapter_removed, int(row["adapter_bases_removed"]))
                self.assertEqual(status, row["status"])

            input_ids = fastq_ids(CONTROLLED_FASTQ)
            quality_ids = fastq_ids(output / "01_quality_trimmed.fastq")
            adapter_ids = fastq_ids(output / "02_adapter_trimmed.fastq")
            retained_ids = fastq_ids(output / "retained_reads.fastq")
            discarded_ids = fastq_ids(output / "discarded_too_short.fastq")
            self.assertEqual(input_ids, quality_ids)
            self.assertEqual(input_ids, adapter_ids)
            self.assertFalse(set(retained_ids) & set(discarded_ids))
            self.assertEqual(set(input_ids), set(retained_ids) | set(discarded_ids))

            reports = [
                json.loads((output / "reports" / name).read_text(encoding="utf-8"))
                for name in [
                    "01_quality.cutadapt.json",
                    "02_adapter.cutadapt.json",
                    "03_length_filter.cutadapt.json",
                ]
            ]
            self.assertEqual(5, reports[0]["read_counts"]["output"])
            self.assertEqual(2, reports[1]["read_counts"]["read1_with_adapter"])
            self.assertEqual(3, reports[2]["read_counts"]["output"])
            self.assertEqual(2, reports[2]["read_counts"]["filtered"]["too_short"])
            self.assertTrue((output / "comparison_dashboard.png").is_file())
            self.assertEqual(9, len(read_csv(output / "before_after_summary.csv")))

        self.assertEqual(checksum_before, hashlib.sha256(CONTROLLED_FASTQ.read_bytes()).hexdigest())

    def test_identifier_descriptions_are_preserved_and_no_discard_is_supported(self) -> None:
        records = [
            ("alpha", "first description", "A" * 100),
            ("beta", "second description", "C" * 90),
        ]
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            fastq = root / "descriptions.fastq"
            content = "".join(
                f"@{read_id} {description}\n{sequence}\n"
                f"+{read_id} {description}\n{'I' * len(sequence)}\n"
                for read_id, description, sequence in records
            )
            fastq.write_text(content, encoding="ascii")
            output = root / "output"
            run_workflow(WorkflowConfig(fastq=fastq, output_dir=output))
            self.assertEqual(["alpha", "beta"], fastq_ids(output / "retained_reads.fastq"))
            self.assertEqual([], fastq_ids(output / "discarded_too_short.fastq", allow_empty=True))
            self.assertEqual(
                ["alpha", "beta"],
                [row["read_id"] for row in read_csv(output / "per_read_trimming.csv")],
            )

    def test_all_reads_can_be_routed_to_too_short_output(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            fastq = root / "short.fastq"
            fastq.write_text("@short\n" + "A" * 20 + "\n+\n" + "I" * 20 + "\n")
            output = root / "output"
            result = run_workflow(WorkflowConfig(fastq=fastq, output_dir=output))
            self.assertEqual(0, result.retained_reads)
            self.assertEqual(1, result.discarded_reads)
            self.assertEqual([], fastq_ids(output / "retained_reads.fastq", allow_empty=True))
            self.assertEqual(["short"], fastq_ids(output / "discarded_too_short.fastq"))
            self.assertTrue((output / "comparison_dashboard.png").is_file())

    def test_full_public_dataset_execution_and_cross_artifact_totals(self) -> None:
        self.assertTrue(PUBLIC_FASTQ.is_file(), "the public test dataset must be committed")
        self.assertEqual(
            PUBLIC_FASTQ_SHA256,
            hashlib.sha256(PUBLIC_FASTQ.read_bytes()).hexdigest(),
        )
        with tempfile.TemporaryDirectory() as temporary_directory:
            output = Path(temporary_directory) / "full_dataset"
            result = run_workflow(
                WorkflowConfig(
                    fastq=PUBLIC_FASTQ,
                    adapter=DEFAULT_ADAPTER,
                    quality_cutoff=20,
                    minimum_length=75,
                    minimum_overlap=8,
                    error_rate=0.1,
                    output_dir=output,
                )
            )
            rows = read_csv(output / "per_read_trimming.csv")
            retained = read_fastq(output / "retained_reads.fastq")
            discarded = read_fastq(output / "discarded_too_short.fastq")
            original = read_fastq(PUBLIC_FASTQ)
            manifest = json.loads((output / "run_manifest.json").read_text(encoding="utf-8"))
            trimming_summary = {
                row["metric"]: row["value"]
                for row in read_csv(output / "trimming_summary.csv")
            }

            self.assertEqual((5000, 4986, 14), (result.input_reads, result.retained_reads, result.discarded_reads))
            self.assertEqual((71, 924), (result.quality_bases_removed, result.adapter_bases_removed))
            self.assertEqual(5000, len(rows))
            self.assertEqual(5000, len({row["read_id"] for row in rows}))
            original_ids = {record.identifier for record in original}
            retained_ids = {record.identifier for record in retained}
            discarded_ids = {record.identifier for record in discarded}
            self.assertFalse(retained_ids & discarded_ids)
            self.assertEqual(original_ids, retained_ids | discarded_ids)
            self.assertEqual(71, sum(int(row["quality_bases_removed"]) for row in rows))
            self.assertEqual(924, sum(int(row["adapter_bases_removed"]) for row in rows))
            retained_bases = sum(record.length for record in retained)
            discarded_bases = sum(record.length for record in discarded)
            self.assertEqual(498107, retained_bases)
            self.assertEqual(898, discarded_bases)
            self.assertEqual(500000, result.total_bases_removed + retained_bases + discarded_bases)
            self.assertEqual("498107", trimming_summary["retained_bases"])
            self.assertEqual("898", trimming_summary["bases_in_discarded_reads"])
            self.assertEqual(PUBLIC_FASTQ_SHA256, manifest["input"]["sha256"])
            self.assertEqual(PUBLIC_FASTQ.name, manifest["input"]["path"])
            self.assertEqual(".", manifest["counts"]["output_dir"])
            self.assertEqual(102, manifest["analysis"]["any_changed_reads"])
            quality_report = json.loads(
                (output / "reports" / "01_quality.cutadapt.json").read_text(encoding="utf-8")
            )
            adapter_report = json.loads(
                (output / "reports" / "02_adapter.cutadapt.json").read_text(encoding="utf-8")
            )
            length_report = json.loads(
                (output / "reports" / "03_length_filter.cutadapt.json").read_text(encoding="utf-8")
            )
            self.assertEqual(71, quality_report["basepair_counts"]["quality_trimmed"])
            self.assertEqual(42, adapter_report["read_counts"]["read1_with_adapter"])
            self.assertEqual(14, length_report["read_counts"]["filtered"]["too_short"])
            portable_artifacts = [
                output / "run_manifest.json",
                *(output / "logs").glob("*.txt"),
                *(output / "reports").glob("*.json"),
            ]
            portable_text = "\n".join(
                path.read_text(encoding="utf-8") for path in portable_artifacts
            )
            self.assertNotIn(str(PUBLIC_FASTQ.resolve()), portable_text)
            self.assertNotIn(str(output.resolve()), portable_text)
            self.assertNotIn(".project09-", portable_text)
            before_positions = read_csv(output / "before_positional_quality.csv")
            after_positions = read_csv(output / "after_positional_quality.csv")
            self.assertEqual(100, len(before_positions))
            self.assertEqual(100, len(after_positions))
            self.assertEqual("5000", before_positions[-1]["reads_at_position"])
            self.assertEqual("4898", after_positions[-1]["reads_at_position"])
            self.assertGreater((output / "comparison_dashboard.png").stat().st_size, 100_000)

        self.assertEqual(
            PUBLIC_FASTQ_SHA256,
            hashlib.sha256(PUBLIC_FASTQ.read_bytes()).hexdigest(),
        )


if __name__ == "__main__":
    unittest.main()
