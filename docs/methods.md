# Methods

## Study design

Project 09 is a single-end preprocessing study. It applies three sequential
Cutadapt 5.2 operations to an authentic public FASTQ subset, retains every
intermediate FASTQ, and constructs an independent per-read accounting table from
the resulting record lengths. The workflow does not align reads or use a reference
genome.

The published configuration is:

| Setting | Value |
|---|---:|
| Quality cutoff | Q20 at the 3-prime end |
| Adapter | `AGATCGGAAGAGCACACGTCTGAACTCCAGTCA` |
| Minimum adapter overlap | 8 nt |
| Maximum adapter error rate | 0.1 |
| Minimum retained length | 75 nt |
| Cutadapt cores per stage | 1 |

The adapter is an Illumina Read 1 trimming sequence documented for TruSeq-based
kits. Cutadapt uses regular 3-prime adapter matching for `-a`, including partial
matches at the read end and error-tolerant matching.

## Dataset selection gate

The first candidate was the same 5,000-read SRR13921545 Read 1 subset used in
Project 08. It would be retained if at least 1% of reads changed or at least 0.5%
of input bases were removed. Otherwise, the planned fallback was complete public
run DRR063281.

The candidate produced 102 changed reads (2.04%) and 995 trimmed bases (0.199%),
so it passed on the changed-read criterion. The fallback was not downloaded or
analyzed. Gate evidence is recorded in `docs/dataset_review.md` and
`results/gate/`.

## Input validation

Before creating an output directory, the program validates:

1. the input exists and is a regular file;
2. the requested output directory does not already exist;
3. the adapter is non-empty and contains supported IUPAC symbols;
4. the quality cutoff is between 0 and 93 inclusive;
5. minimum length and overlap are positive integers;
6. overlap does not exceed adapter length;
7. error rate is in the half-open interval `[0, 1)`;
8. each FASTQ record contains exactly four lines;
9. headers and separators begin with `@` and `+`;
10. a non-empty repeated `+` description equals the complete header description;
11. identifier tokens are non-empty and unique;
12. sequences are non-empty IUPAC strings; and
13. sequence and Phred+33 quality strings have equal length.

Gzip-compressed input is detected from the `.gz` suffix. Empty sequences are
rejected in raw input but supported in intermediate and discarded outputs because
trimming can legitimately reduce a read to length zero.

## Stage 1: quality trimming

Cutadapt runs with `--quality-cutoff 20`. With one cutoff, Cutadapt applies quality
trimming at the 3-prime end and assumes ASCII-encoded Phred+33 qualities. The
quality-trimming algorithm is the same algorithm used by BWA.

The stage produces:

- `01_quality_trimmed.fastq`;
- `reports/01_quality.cutadapt.json`; and
- captured standard-output and standard-error logs.

All input reads remain present after this modification stage.

## Stage 2: adapter trimming

The post-quality FASTQ is processed with:

```text
--adapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
--overlap 8
--error-rate 0.1
```

Indels remain enabled, which is Cutadapt's default. At an eight-base partial
overlap, a 10% maximum rate permits no errors because one error would have rate
1/8 = 0.125. Longer matches can tolerate errors when the calculated rate remains
at or below 0.1.

The stage produces `02_adapter_trimmed.fastq`, its JSON report, and logs. No reads
are filtered at this stage.

## Stage 3: minimum-length routing

The post-adapter FASTQ is processed with `--minimum-length 75`. Reads below 75 nt
are redirected with `--too-short-output` instead of being deleted.

The stage produces:

- `retained_reads.fastq` for reads of at least 75 nt;
- `discarded_too_short.fastq` for reads below 75 nt;
- `reports/03_length_filter.cutadapt.json`; and
- logs.

Cutadapt specifies that a read is written to at most one filtering output. The
workflow additionally verifies that retained and too-short identifier sets are
disjoint and their union equals the input identifier set.

## Per-read accounting

The program parses the raw, post-quality, post-adapter, retained, and too-short
FASTQ files and joins them by the first whitespace-delimited identifier token.
For every input read it records:

```text
quality bases removed = original length - post-quality length
adapter bases removed = post-quality length - post-adapter length
total bases removed   = quality removed + adapter removed
final length          = length in retained or too-short FASTQ
```

The final status is `retained` or `discarded_too_short`. Records remain in input
order in `per_read_trimming.csv`.

## Cross-artifact verification

Direct FASTQ record and base counts are compared with every Cutadapt JSON report.
The workflow also checks:

```text
adapter-stage bases = retained bases + too-short bases
input bases = quality-trimmed bases
            + adapter-trimmed bases
            + retained bases
            + too-short bases
```

The input SHA-256 is measured before processing and again before publication.

## Atomic publication

All outputs are generated in a temporary sibling directory. Any validation,
Cutadapt, JSON, accounting, analysis, or filesystem failure removes this staging
directory. On success, the completed directory is renamed to the requested output
path in one operation. Existing output directories are rejected.

## Before/after analysis

"Before" includes all raw reads. "After" includes only the analysis-ready retained
reads; too-short reads are summarized separately.

For each set, the workflow calculates read count, total bases, mean/median/minimum/
maximum length, mean Phred quality, percent of bases at Q20 or higher, and percent
at Q30 or higher.

For every sequenced position, the positional tables report:

- number of reads reaching the position;
- mean and median Phred score;
- percent of bases at Q20 or higher; and
- percent of bases at Q30 or higher.

The dashboard contains read disposition, before/after length distribution, mean
quality by position, and base-level consequences. The final panel distinguishes
trimmed bases from bases preserved in too-short reads.

## Software and reproducibility

- Python 3.12
- Cutadapt 5.2
- Matplotlib 3.11.1 using the non-interactive Agg backend
- one Cutadapt core per stage
- standard-library CSV, JSON, FASTQ validation, hashing, subprocess, and testing

Exact runtime parameters, versions, checksum, counts, and analysis metrics are
stored in `run_manifest.json`.

Before outputs are published, machine-specific absolute paths in Cutadapt JSON
reports and captured logs are replaced with portable filenames. The manifest
also records the input filename and `.` as the output location so a run does not
disclose or depend on the directory layout of the computer that produced it.
