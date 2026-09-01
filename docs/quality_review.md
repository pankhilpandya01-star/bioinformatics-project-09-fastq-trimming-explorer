# Milestone 4: Quality Review

## Status

The automated quality suite, clean package-build check, dependency check, and local
Python 3.12 execution all pass. A GitHub Actions workflow is ready to run the same
suite on `ubuntu-latest` with Python 3.12 once the repository is published.

## Automated coverage

| Acceptance area | Automated evidence |
|---|---|
| Malformed FASTQ | Rejects empty, truncated, invalid header/separator, mismatched sequence and quality lengths, empty sequences, non-IUPAC symbols, invalid Phred+33 characters, mismatched repeated descriptions, and duplicate identifiers. |
| Invalid CLI/configuration | Rejects empty or invalid adapters, Q-scores outside 0–93, invalid minimum lengths and overlaps, error rates outside the half-open interval 0–1, missing inputs, and existing output directories. The CLI returns exit code 2 and creates no output for an invalid Q-score. |
| Known quality trimming | A controlled Q0 tail removes exactly 10 bases from each of two reads. |
| Known adapter trimming | Two exact 33-base adapter tails each remove exactly 33 bases. |
| Variable lengths | Controlled original lengths of 80, 100, 103, and 113 nt produce verified final lengths of 70, 80, 90, and 100 nt. |
| Subprocess failures | A mocked Cutadapt exit code 17 raises a diagnostic workflow error, removes the staging directory, and leaves no public output directory. |
| Identifier preservation | Input order is preserved through quality and adapter intermediates; retained and discarded identifiers form a disjoint exact input partition. Headers with descriptions are also exercised. |
| Empty routing branches | Tests cover both no discarded reads and all reads discarded as too short. |
| Exact accounting | Controlled CSV, FASTQ, JSON, manifest, positional-quality, and dashboard artifacts agree exactly. |
| Full dataset | The authentic 5,000-read SRR13921545 subset is processed end to end and reproduces 4,986 retained reads, 14 too-short reads, 71 quality-trimmed bases, 924 adapter-trimmed bases, 498,107 retained bases, and 898 preserved too-short bases. |
| Input immutability | The full-dataset SHA-256 is verified before and after execution. |

## Defect found during review

The first suite run exposed a FASTQ-description validation gap. The parser initially
allowed a `+` line to repeat only the identifier token while Cutadapt requires a
non-empty repeated description to match the complete header description. The
parser now enforces Cutadapt's rule before starting any subprocess. An unrelated
test assertion also expected five before/after metrics even though the file
correctly contains nine; that assertion was corrected.

The next full run passed all tests.

The final documentation audit identified one additional Cutadapt boundary: a
value of exactly `1` for `-e` means one allowed error rather than a 100% error
rate. Because this project's option is explicitly named `--error-rate`, validation
now requires a value from 0 inclusive to 1 exclusive, and the suite covers that
boundary.

## Local verification

The project was built as a non-editable wheel and reinstalled before the final test
run, ensuring the suite exercised the packaged code rather than only the source
tree.

```text
Ran 10 tests in 6.009s

OK
No broken requirements found.
```

Verification command:

```bash
python -m unittest discover -s tests -v
```

## GitHub Actions

`.github/workflows/tests.yml` runs on pushes, pull requests, and manual dispatches.
It grants read-only repository permissions, uses `ubuntu-latest`, installs Python
3.12, installs the packaged project from `pyproject.toml`, and runs the same
standard-library test suite with a 15-minute timeout. Matplotlib is forced to its
non-interactive backend and a runner-temporary cache directory.

The workflow uses the current official major versions `actions/checkout@v7` and
`actions/setup-python@v7`. Remote CI has not been triggered because the agreed
publication boundary still prohibits creating or pushing the GitHub repository.

## Review boundary

Milestone 5 will complete the README, methods, limitations, provenance, citations,
portfolio navigation, and final repository audit. It will pause again before any
GitHub repository is created or published.
