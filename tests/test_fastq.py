from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from fastq_trimming_explorer.fastq import FastqValidationError, read_fastq


class FastqValidationTests(unittest.TestCase):
    def test_rejects_malformed_fastq_records(self) -> None:
        cases = {
            "empty": "",
            "truncated": "@read1\nACGT\n+\n",
            "bad_header": "read1\nACGT\n+\nIIII\n",
            "bad_separator": "@read1\nACGT\n-\nIIII\n",
            "mismatched_lengths": "@read1\nACGT\n+\nIII\n",
            "empty_sequence": "@read1\n\n+\n\n",
            "invalid_sequence_symbol": "@read1\nACZT\n+\nIIII\n",
            "invalid_quality_character": "@read1\nA\n+\n \n",
            "mismatched_repeated_id": "@read1\nACGT\n+read2\nIIII\n",
            "duplicate_identifier": (
                "@read1 first\nACGT\n+\nIIII\n"
                "@read1 second\nTGCA\n+\nIIII\n"
            ),
        }
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            for name, content in cases.items():
                with self.subTest(case=name):
                    path = root / f"{name}.fastq"
                    path.write_text(content, encoding="ascii")
                    with self.assertRaises(FastqValidationError):
                        read_fastq(path)

    def test_accepts_iupac_bases_descriptions_and_repeated_identifier(self) -> None:
        content = (
            "@read1 description\nACGTRYSWKMBDHVN\n"
            "+read1 description\nIIIIIIIIIIIIIII\n"
        )
        with tempfile.TemporaryDirectory() as temporary_directory:
            path = Path(temporary_directory) / "valid.fastq"
            path.write_text(content, encoding="ascii")
            records = read_fastq(path)
        self.assertEqual(["read1"], [record.identifier for record in records])
        self.assertEqual("ACGTRYSWKMBDHVN", records[0].sequence)
        self.assertEqual("IIIIIIIIIIIIIII", records[0].quality)

    def test_allows_empty_cutadapt_output_when_explicitly_requested(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            path = Path(temporary_directory) / "empty.fastq"
            path.write_text("", encoding="ascii")
            self.assertEqual([], read_fastq(path, allow_empty=True))


if __name__ == "__main__":
    unittest.main()
