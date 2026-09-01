# Milestone 1: Dataset Review

## Decision

Keep the Project 08 test dataset for Project 09. The suitability gate passed because
102 of 5,000 reads (2.04%) changed during the two trimming stages, exceeding the
1% changed-read threshold. The fallback run, DRR063281, is therefore not needed.

## Source and integrity

- Accession: SRR13921545
- Organism: *Escherichia coli* K-12 MG1655
- Library strategy: whole-genome sequencing (WGS)
- Source layout: paired-end; Project 09 uses Read 1 only
- Instrument: Illumina NovaSeq 6000
- Subset: first 5,000 spots, Read 1
- Source command: `fastq-dump --split-files --skip-technical -X 5000 SRR13921545`
- FASTQ records: 5,000
- Original bases: 500,000
- SHA-256: `c3d708ba4b1c7eb4bb95dbae6f7189f291250d129514f4d26b50aad272fecc15`

The checksum matched the Project 08 provenance metadata before trimming and was
unchanged after the gate run. The original FASTQ was used only as an input; all
gate outputs were written under `results/gate/`.

## Gate configuration

Cutadapt 5.2 was run as two separate stages:

1. 3-prime quality trimming at Q20.
2. 3-prime adapter trimming with adapter
   `AGATCGGAAGAGCACACGTCTGAACTCCAGTCA`, minimum overlap 8, maximum error
   rate 0.1, and minimum retained length 75 nt. Too-short reads were preserved
   separately.

The dataset would pass if either at least 1% of reads changed or at least 0.5% of
input bases were removed.

## Observed trimming signal

| Measure | Result |
|---|---:|
| Input reads | 5,000 |
| Reads changed by quality trimming | 60 (1.20%) |
| Reads changed by adapter trimming | 42 (0.84%) |
| Reads changed by either stage | 102 (2.04%) |
| Quality-trimmed bases | 71 |
| Adapter-trimmed bases | 924 |
| Total trimmed bases | 995 (0.199% of input bases) |
| Retained reads (at least 75 nt) | 4,986 |
| Too-short reads preserved separately | 14 |
| Shortest retained read | 77 nt |
| Longest discarded read | 74 nt |

The quality and adapter stages changed disjoint sets of reads in this gate run,
so the changed-read total is 60 + 42 = 102.

## Accounting checks

- All 5,000 input identifiers were preserved through the quality stage.
- The retained and too-short outputs were disjoint.
- Their identifier union exactly matched the 5,000 input identifiers.
- The quality report recorded 500,000 input bases, 71 quality-trimmed bases, and
  499,929 post-quality bases.
- The adapter report recorded 42 adapter matches, 14 too-short reads, and 4,986
  retained reads.
- Adapter-trimmed lengths summed to 924 bases.
- The adapter report's retained-output base count excludes the remaining bases of
  the 14 reads routed to the too-short output. Those reads are preserved, not lost.

## Interpretation

This is a modest but measurable preprocessing signal. Most reads remain unchanged,
which is expected for a high-quality short-read subset, while the authentic data
still exercises quality trimming, adapter removal, and minimum-length routing.
Passing the gate does not imply that every removed base is biologically invalid;
it establishes that the dataset can demonstrate and quantify the workflow's
effects without synthetic reads in the published analysis.
