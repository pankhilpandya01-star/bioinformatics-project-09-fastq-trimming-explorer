# Data Provenance

## Public source

The authentic analysis dataset derives from public NCBI Sequence Read Archive run
[`SRR13921545`](https://www.ncbi.nlm.nih.gov/sra/SRR13921545).

| Field | Recorded value |
|---|---|
| Organism | *Escherichia coli* K-12 MG1655 |
| Library strategy | WGS |
| Library layout | Paired |
| Platform | Illumina |
| Instrument | Illumina NovaSeq 6000 |
| Selection used here | Read 1 from spots 1–5,000 |
| Bundled records | 5,000 |
| Bundled bases | 500,000 |
| Retrieval date | 2026-08-05 |
| SRA Toolkit version | 3.4.1 |

Only Read 1 is analyzed. The project does not claim to process the paired experiment
as a synchronized pair.

## Acquisition

The subset was generated with the NCBI SRA Toolkit command recorded by Project 08:

```powershell
fastq-dump --split-files --skip-technical -X 5000 SRR13921545
```

The bundled file is:

```text
data/SRR13921545_read1_first5000.fastq
```

Its SHA-256 is:

```text
C3D708BA4B1C7EB4BB95DBAE6F7189F291250D129514F4D26B50AAD272FECC15
```

Machine-readable metadata is stored in `data/source_metadata.csv`. The FASTQ is a
byte-for-byte copy of the verified Project 08 input so Projects 08 and 09 examine
the same records from diagnosis through preprocessing.

## Dataset selection

The selection rule was fixed before the workflow was implemented:

- keep this subset if at least 1% of reads changed; or
- keep it if at least 0.5% of bases were removed; otherwise
- use complete public single-end run DRR063281.

Under Q20 quality trimming, the specified Read 1 adapter, eight-base overlap, 10%
maximum error rate, and 75 nt minimum length, 102 reads changed (2.04%) and 995
bases were trimmed (0.199%). The read criterion passed, so DRR063281 was not used.

The gate is documented in `docs/dataset_review.md`; its Cutadapt reports and FASTQ
outputs are in `results/gate/`.

## Input immutability

The workflow opens the source FASTQ only for reading. It hashes the input before
processing and again immediately before atomically publishing results. Tests also
verify the public checksum before and after a complete run.

The verified checksum remained unchanged throughout all milestones.

## Authentic and constructed data

- All published biological results in `results/full_dataset/` use authentic public
  SRR13921545 reads.
- `tests/fixtures/controlled.fastq` contains five deliberately constructed reads
  used only to verify known quality, adapter, and length-routing outcomes.
- Controlled fixtures are not mixed into the biological summaries or dashboard.

## Derived data

`results/full_dataset/` contains derived FASTQ, CSV, JSON, log, manifest, and PNG
artifacts. Every input identifier appears once in `per_read_trimming.csv` and once
in either the retained or too-short FASTQ. The source FASTQ remains separate and
unchanged under `data/`.

