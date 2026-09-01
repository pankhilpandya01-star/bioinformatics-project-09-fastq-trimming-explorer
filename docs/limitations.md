# Limitations

## Biological scope

- The analysis uses the first 5,000 Read 1 records from one *E. coli* run. It is a
  reproducible subset, not a random sample and not evidence about the complete run,
  other organisms, library preparations, or sequencing instruments.
- The source experiment is paired-end, but this project intentionally processes
  only Read 1. Filtering one mate independently can desynchronize paired data, so
  this workflow must not be substituted for a paired-end pipeline.
- Adapter matches are computational matches under the configured overlap and error
  rules. They are consistent with adapter sequence but do not prove the cause of
  every match or that every removed base was biologically invalid.
- The 75 nt threshold is an educational project decision. An appropriate threshold
  depends on the downstream aligner, assay, read length, and research question.
- The quality and adapter settings were not optimized against downstream mapping,
  assembly, or variant-calling performance.

## Quality-trimming scope

- The source instrument is an Illumina NovaSeq 6000, which uses two-color chemistry.
  Cutadapt documents a specialized `--nextseq-trim` mode for high-quality erroneous
  terminal `G` calls from two-color instruments. This project follows the planned
  regular Q20 3-prime trimming design and does not evaluate that alternative.
- The parser assumes Phred+33 encoding. Legacy Phred+64 FASTQ is not supported.
- A Phred threshold describes sequencing confidence; it does not directly determine
  whether a base is biologically correct or useful downstream.

## Adapter-model scope

- The default sequence is the documented Illumina Read 1 trimming sequence used by
  multiple TruSeq-based kits, but the SRA metadata bundled here does not include the
  original sample sheet. The workflow therefore tests an explicit adapter hypothesis
  rather than reconstructing the exact laboratory configuration.
- Cutadapt's default indel-aware matching is retained. Results could differ with
  `--no-indels`, anchoring, non-internal matching, a different overlap, or another
  adapter sequence.
- The workflow searches one 3-prime adapter. It does not trim 5-prime adapters,
  poly-A tails, primers, UMIs, or multiple adapter families.

## Comparison scope

- "Before" uses all raw reads, whereas "after" uses only retained reads. Changes in
  after-quality metrics therefore reflect both trimming and the exclusion of 14
  too-short reads.
- Positional statistics at later positions use fewer retained reads because trimmed
  reads no longer reach those positions. `reads_at_position` must be considered when
  interpreting the quality curves.
- The observed quality change is very small. The dashboard is descriptive and does
  not include confidence intervals or inferential statistics.

## Software scope

- Input must use unwrapped four-line FASTQ records with unique first-token
  identifiers. Multi-line FASTQ variants and duplicate identifier tokens are
  rejected so exact per-read accounting remains possible.
- Cutadapt is intentionally pinned to version 5.2. A different version is rejected
  rather than silently changing report schemas or behavior.
- Output directories must be new. The CLI does not provide an overwrite flag.
- One Cutadapt core is used for deterministic, easy-to-audit project runs rather
  than maximum throughput.
- The project does not perform FastQC/MultiQC integration, contamination screening,
  reference alignment, assembly, taxonomic classification, duplicate removal,
  variant analysis, or downstream biological validation.

