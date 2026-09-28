# SRR6468562 example

This example contains paired non-human metagenomic reads and the compact result tables produced by CARDCell.

## Input files

- `input/SRR6468562.nonhuman_R1.fastq.gz`
- `input/SRR6468562.nonhuman_R2.fastq.gz`

Each mate contains 469,775 reads. The combined mean read length is 96.212222872652 bp.

SHA-256 checksums:

```text
99a72d9f47083ba3764bd5f6f7ba567d785d05157260a00985f1e85e24770ddb  SRR6468562.nonhuman_R1.fastq.gz
05ad9f34994335d765ca0521c0e39f9a37cdc3cfcd2159770aac05be5f440766  SRR6468562.nonhuman_R2.fastq.gz
```

## Reproduce the example

Run from the repository root:

```bash
bin/cardoap_kma16s \
  --r1 examples/SRR6468562/input/SRR6468562.nonhuman_R1.fastq.gz \
  --r2 examples/SRR6468562/input/SRR6468562.nonhuman_R2.fastq.gz \
  --sample SRR6468562 \
  --outdir results \
  --threads 8 \
  --card-local-db /path/to/localDB \
  --metaxa-dir /path/to/Metaxa2_2.2.3 \
  --metaxa-env-bin /path/to/metaxa2_dependency_bin
```

The published compact outputs are stored under `expected_results/`. Large SAM/BAM files and temporary files are intentionally excluded.

Drug classes are assigned directly from the repository's bundled `CARD_OAP_full_6059_annotation.csv` by exact ARO accession. The validation sample has 15/15 exact annotation matches.

Expected headline results:

```text
Bacterial 16S reads: 406
Cell-equivalent coverage: 5.455958682164
CARD references detected: 15
Total KMA ARG Depth: 5.600000000000
Total ARG copies/cell: 1.026400734722
```
