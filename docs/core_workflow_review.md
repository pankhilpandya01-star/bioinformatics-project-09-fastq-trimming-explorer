# Milestone 2: Core Workflow Review

## Status

The validated command-line workflow is implemented and demonstrated on controlled
single-end reads. Real-data analysis has not started.

## Interface

The installed command is `fastq-trimming-explorer` and exposes every planned
parameter:

```text
--fastq
--adapter
--quality-cutoff
--minimum-length
--minimum-overlap
--error-rate
--output-dir
```

The default adapter is
`AGATCGGAAGAGCACACGTCTGAACTCCAGTCA`; the remaining defaults are Q20, minimum
length 75 nt, minimum overlap 8 nt, and maximum error rate 0.1.

## Workflow architecture

Cutadapt 5.2 runs in three independent subprocesses:

1. `01_quality`: 3-prime quality trimming.
2. `02_adapter`: 3-prime adapter trimming.
3. `03_length_filter`: routing into retained and too-short FASTQ files.

Each stage writes its own JSON report and captured standard-output/error logs.
Keeping the length filter separate means adapter-trimmed bases are not confused
with the bases that remain in reads routed to the too-short file.

Before Cutadapt runs, the workflow checks the input path, output path, adapter,
numeric ranges, FASTQ structure, unique identifiers, IUPAC sequence symbols,
sequence/quality lengths, and Phred+33 characters. It requires the installed
Cutadapt version to be exactly 5.2.

The workflow builds all files in a temporary sibling directory. It publishes the
requested output directory only after Cutadapt succeeds, all identifiers are
accounted for, direct FASTQ totals agree with the JSON reports, and the input
checksum is unchanged. A failed run removes its staging directory and does not
leave a partial output directory.

## Output contract

Every successful run produces:

- `01_quality_trimmed.fastq`
- `02_adapter_trimmed.fastq`
- `retained_reads.fastq`
- `discarded_too_short.fastq`
- `per_read_trimming.csv`
- `run_manifest.json`
- `reports/01_quality.cutadapt.json`
- `reports/02_adapter.cutadapt.json`
- `reports/03_length_filter.cutadapt.json`
- one standard-output and standard-error log per stage under `logs/`

The per-read table preserves input order and records the original, post-quality,
post-adapter, and final length; bases removed by each trimming stage; total bases
removed; and retained or discarded status.

## Controlled demonstration

The five-read fixture deliberately contains an unchanged read, a read with a
low-quality tail, an adapter-contaminated retained read, an adapter-contaminated
too-short read, and a quality-trimmed too-short read.

| Read | Original | Post-quality | Final | Quality removed | Adapter removed | Status |
|---|---:|---:|---:|---:|---:|---|
| `unchanged` | 100 | 100 | 100 | 0 | 0 | retained |
| `quality_tail` | 100 | 90 | 90 | 10 | 0 | retained |
| `adapter_retained` | 113 | 113 | 80 | 0 | 33 | retained |
| `adapter_discarded` | 103 | 103 | 70 | 0 | 33 | discarded too short |
| `quality_discarded` | 80 | 70 | 70 | 10 | 0 | discarded too short |

Verified totals:

- 5 input reads and 496 input bases
- 3 retained reads and 2 too-short reads
- 20 quality-trimmed bases
- 66 adapter-trimmed bases
- 86 total trimmed bases
- 5 unique per-read rows
- all three Cutadapt JSON reports agree with direct FASTQ accounting
- the controlled input checksum remained unchanged

An invalid Q-score demonstration returned exit code 2 and created no output
directory. Package installation, the console entry point, byte-code compilation,
and dependency consistency checks also passed.

## Review boundary

Milestone 3 will run this workflow on the authentic 5,000-read dataset and add the
before/after summary tables, positional-quality tables, four-panel dashboard, and
biological interpretation. Those analysis artifacts are intentionally deferred
until this milestone is approved.
