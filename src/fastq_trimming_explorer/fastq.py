"""Strict FASTQ parsing used for validation and read accounting."""

from __future__ import annotations

import gzip
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO


FASTQ_SEQUENCE_ALPHABET = frozenset("ACGTURYSWKMBDHVN")


class FastqValidationError(ValueError):
    """Raised when an input is not a valid four-line FASTQ file."""


@dataclass(frozen=True, slots=True)
class FastqRecord:
    identifier: str
    sequence: str
    quality: str

    @property
    def length(self) -> int:
        return len(self.sequence)


@contextmanager
def open_fastq(path: Path) -> Iterator[TextIO]:
    """Open an uncompressed or gzip-compressed FASTQ as ASCII text."""

    if path.suffix.lower() == ".gz":
        with gzip.open(path, mode="rt", encoding="ascii", newline=None) as handle:
            yield handle
    else:
        with path.open(mode="rt", encoding="ascii", newline=None) as handle:
            yield handle


def _strip_newline(value: str) -> str:
    return value.rstrip("\r\n")


def read_fastq(
    path: Path,
    *,
    allow_empty: bool = False,
    allow_empty_sequences: bool = False,
) -> list[FastqRecord]:
    """Read and validate a standard four-line FASTQ file.

    Identifiers are the first whitespace-delimited token after ``@``. They must be
    unique because they are the stable keys used for cross-stage accounting.
    """

    records: list[FastqRecord] = []
    identifiers: set[str] = set()

    try:
        with open_fastq(path) as handle:
            record_number = 0
            while True:
                header_line = handle.readline()
                if header_line == "":
                    break

                record_number += 1
                sequence_line = handle.readline()
                plus_line = handle.readline()
                quality_line = handle.readline()
                if "" in (sequence_line, plus_line, quality_line):
                    raise FastqValidationError(
                        f"record {record_number} is incomplete; FASTQ records require four lines"
                    )

                header = _strip_newline(header_line)
                sequence = _strip_newline(sequence_line)
                plus = _strip_newline(plus_line)
                quality = _strip_newline(quality_line)

                if not header.startswith("@"):
                    raise FastqValidationError(
                        f"record {record_number} header must start with '@'"
                    )
                header_body = header[1:].strip()
                identifier = header_body.split(maxsplit=1)[0] if header_body else ""
                if not identifier:
                    raise FastqValidationError(
                        f"record {record_number} has an empty identifier"
                    )
                if identifier in identifiers:
                    raise FastqValidationError(
                        f"duplicate FASTQ identifier: {identifier!r}"
                    )
                if not plus.startswith("+"):
                    raise FastqValidationError(
                        f"record {record_number} separator must start with '+'"
                    )
                plus_body = plus[1:].strip()
                if plus_body and plus_body != header_body:
                    raise FastqValidationError(
                        f"record {record_number} repeats a different description on its '+' line"
                    )
                if not sequence and not allow_empty_sequences:
                    raise FastqValidationError(
                        f"record {record_number} has an empty sequence"
                    )
                if any(character.isspace() for character in sequence):
                    raise FastqValidationError(
                        f"record {record_number} contains whitespace in its sequence"
                    )
                invalid_sequence_symbols = sorted(
                    set(sequence.upper()) - FASTQ_SEQUENCE_ALPHABET
                )
                if invalid_sequence_symbols:
                    symbols = "".join(invalid_sequence_symbols)
                    raise FastqValidationError(
                        f"record {record_number} contains non-IUPAC sequence symbols: {symbols}"
                    )
                if len(sequence) != len(quality):
                    raise FastqValidationError(
                        f"record {record_number} has sequence length {len(sequence)} "
                        f"but quality length {len(quality)}"
                    )
                if any(not 33 <= ord(character) <= 126 for character in quality):
                    raise FastqValidationError(
                        f"record {record_number} contains a quality character outside Phred+33"
                    )

                identifiers.add(identifier)
                records.append(
                    FastqRecord(
                        identifier=identifier,
                        sequence=sequence,
                        quality=quality,
                    )
                )
    except (OSError, UnicodeError) as error:
        raise FastqValidationError(f"could not read FASTQ file: {error}") from error

    if not records and not allow_empty:
        raise FastqValidationError("FASTQ file contains no records")
    return records


def records_by_identifier(records: list[FastqRecord]) -> dict[str, FastqRecord]:
    return {record.identifier: record for record in records}
