# Milestone 3: Authentic Dataset Analysis

## Status

The complete workflow was run on the authentic 5,000-read SRR13921545 Read 1
subset selected in Milestone 1. The FASTQ outputs, Cutadapt reports, per-read
records, summary tables, positional-quality tables, manifest, and four-panel
dashboard were generated together in one verified atomic run.

## Comparison definition

"Before" means all 5,000 raw input reads. "After" means the 4,986 analysis-ready
reads retained at a minimum length of 75 nt. The 14 too-short reads are not included
in after-quality or after-length statistics, but remain available in
`discarded_too_short.fastq` and are reported separately throughout the analysis.

## Primary results

| Measure | Result |
|---|---:|
| Input reads | 5,000 |
| Retained reads | 4,986 (99.72%) |
| Too-short reads | 14 (0.28%) |
| Reads changed by Q20 trimming | 60 (1.20%) |
| Reads changed by adapter trimming | 42 (0.84%) |
| Reads changed by either stage | 102 (2.04%) |
| Quality-trimmed bases | 71 |
| Adapter-trimmed bases | 924 |
| Total trimmed bases | 995 (0.199% of input) |
| Bases retained for analysis | 498,107 (99.6214% of input) |
| Bases preserved in too-short reads | 898 |

The quality-changed and adapter-changed read sets did not overlap. All 14 too-short
reads were adapter-trimmed; the other 28 adapter-trimmed reads remained at least
75 nt and were retained.

## Before/after quality and length

| Measure | Before | After retained reads |
|---|---:|---:|
| Mean read length | 100.000 nt | 99.901 nt |
| Median read length | 100 nt | 100 nt |
| Minimum read length | 100 nt | 77 nt |
| Maximum read length | 100 nt | 100 nt |
| Mean base quality | Q36.4305 | Q36.4335 |
| Bases at Q20 or higher | 99.0076% | 99.0205% |
| Bases at Q30 or higher | 96.4124% | 96.4225% |
| Mean quality at position 100 | Q36.2224 | Q36.5272 |

The small quality increase is consistent with removing only a limited set of
low-quality or adapter-associated tails. It should not be interpreted as a broad
quality problem in the raw subset: nearly all bases were already high quality, the
median read length stayed at 100 nt, and 4,898 retained reads remained exactly
100 nt long.

The 42 adapter matches provide a measurable preprocessing signal. They are
consistent with adapter sequence occurring at the 3-prime ends of a small subset
of reads, but this workflow does not establish why each match occurred or claim
that every removed base was biologically invalid. The 14 reads that fell below
75 nt were separated rather than destroyed so the filtering decision remains
auditable.

## Generated artifacts

- `01_quality_trimmed.fastq`: post-Q20 intermediate
- `02_adapter_trimmed.fastq`: post-adapter intermediate before length routing
- `retained_reads.fastq`: 4,986 analysis-ready reads
- `discarded_too_short.fastq`: 14 preserved too-short reads
- `per_read_trimming.csv`: 5,000 input-ordered accounting records
- `before_after_summary.csv`: before/retained comparison metrics
- `trimming_summary.csv`: trimming and routing totals
- `before_positional_quality.csv`: raw positional quality through position 100
- `after_positional_quality.csv`: retained-read positional quality through position 100
- `comparison_dashboard.png`: four-panel comparison dashboard
- `run_manifest.json`: input checksum, versions, parameters, counts, and analysis metrics
- `reports/`: one Cutadapt JSON report per processing stage
- `logs/`: captured Cutadapt standard output and error for each stage

## Consistency verification

- The original SHA-256 remained
  `c3d708ba4b1c7eb4bb95dbae6f7189f291250d129514f4d26b50aad272fecc15`.
- All 5,000 input identifiers appear once in the per-read table.
- Retained and discarded identifiers are disjoint and their union is the input.
- Per-read removal totals are 71 quality bases and 924 adapter bases.
- The three Cutadapt JSON reports agree with the FASTQ record and base counts.
- 500,000 input bases reconcile exactly as 995 trimmed bases, 498,107 retained
  bases, and 898 bases preserved in too-short reads.
- Both positional-quality CSVs contain all 100 sequenced positions.
- The dashboard was visually inspected after generation.

## Review boundary

Milestone 4 will expand automated coverage for malformed FASTQ, invalid parameter
values, known trimming outcomes, variable lengths, subprocess failures, identifier
preservation, exact accounting, and the full public dataset. It will also add the
Python 3.12 GitHub Actions workflow. That quality-review work is intentionally
deferred until this milestone is approved.
