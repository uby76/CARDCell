# CARDCell: CARD/KMA ARG Quantification with 16S-Based Cell Normalization

CARDCell is a standalone, resumable pipeline for quantifying ARG abundance in paired-end metagenomic data.

## Clone the repository

```bash
git clone git@github.com:uby76/CARDCell.git
cd CARDCell
chmod +x bin/cardoap_kma16s
```

Requirements include Python 3.9 or later, RGI with its KMA, samtools, and bamtools runtime dependencies, a loaded canonical CARD local database, and Metaxa2 2.2.3 with its SSU database. The rrnDB 5.10 pan-taxa NCBI statistics table is included. CARD and Metaxa2 databases are not downloaded or replaced automatically.

## Method overview

```text
ARG numerator: FASTQ → RGI-bwt/KMA → CARD → KMA-reported CARD reference Depth
Cell denominator: FASTQ → Metaxa2 bacterial 16S → 16S coverage → taxonomy → rrnDB → cell-equivalent coverage
Final result: ARG copies/cell = KMA Depth / cell-equivalent coverage
```

The pipeline combines **CARD/KMA ARG profiling** with **[Zhu et al. 2025](https://doi.org/10.1038/s41467-025-59019-3)-style bacterial cell normalization based on 16S and rrnDB**. The complete pipeline should not be described as the original Zhu et al. 2025 workflow; only the cell denominator follows its 16S normalization concept.

## Core formulas

Read lengths are calculated directly from every FASTQ record. No fixed read length such as 100 or 150 bp is assumed.

```text
16S_coverage = N_16S_total × combined_mean_read_length / 1432

weighted_mean_16S_copy_number
  = Σ(relative_abundance_i × rrn_copy_number_i)

cell_equivalent_coverage
  = 16S_coverage / weighted_mean_16S_copy_number

ARG_copies_per_cell
  = KMA-reported_CARD_reference_Depth / cell_equivalent_coverage
```

The 1,432 bp reference length is the average bacterial 16S rRNA gene length used by [Zhu et al. 2025](https://doi.org/10.1038/s41467-025-59019-3). The ARG numerator is taken exclusively from the `Depth` field in the RGI-bwt/KMA `allele_mapping_data.txt` output. It is explicitly referred to as **KMA-reported CARD reference Depth**.

The following are excluded from the primary calculation: BAM-derived mean depth, ARGs-OAP KO30 nCell, the experimental KMA-KO30 database, assembly, ORF prediction, DIAMOND, BLASTX, DeepARG, and Kraken2/Bracken. RGI internally generates BAM files and invokes samtools/bamtools, but the CARDCell calculation code does not read BAM files.

## rrnDB matching

Metaxa2 identifies bacterial SSU/16S reads independently in R1 and R2. Community relative abundance is calculated using all bacterial 16S reads. Reads assigned to genus are summarized at genus level; otherwise, the deepest reliable taxonomic assignment is retained.

The first exact rank/name match in rrnDB is used in this order:

```text
genus → family → order → class → phylum → domain
```

Every fallback is recorded in `04_rrnDB_copy_number.tsv`. No arbitrary fixed value such as 4 or 4.1 is assigned to missing taxa. If unmatched abundance is present, the matched abundance is renormalized before calculating the weighted mean. The matched fraction, unmatched fraction, and renormalization status are reported.

## Usage

```bash
bin/cardoap_kma16s \
  --r1 sample_R1.fastq.gz \
  --r2 sample_R2.fastq.gz \
  --sample SAMPLE_ID \
  --outdir results \
  --threads 8 \
  --card-local-db /path/to/localDB \
  --metaxa-dir /path/to/Metaxa2_2.2.3 \
  --metaxa-env-bin /path/to/metaxa2_dependency_bin
```

By default, the pipeline uses `resources/rrnDB-5.10_pantaxa_stats_NCBI.tsv`. A different current official pan-taxa NCBI statistics table can be supplied with `--rrndb FILE`. CARD is never downloaded or replaced automatically.

If RGI, Metaxa2, or their databases are not available in the default locations, specify them with `--rgi-bin`, `--card-local-db`, `--metaxa-dir`, and `--metaxa-env-bin`. To list all options:

```bash
bin/cardoap_kma16s --help
```

Common options:

- `--force`: rebuild only the selected sample inside the specified output directory without modifying input data or other workflows.
- `--keep-temp`: retain decompressed temporary FASTQ files.
- Without `--force`, valid non-empty RGI and Metaxa2 outputs are detected and reused automatically.
- `--comparison-ko30`, `--comparison-previous-16s`, and `--comparison-previous-16s-count`: write validation comparisons only; these values never enter the primary calculation.

## Quality control

Default thresholds produce transparent flags and never delete ARG detections:

- mapped reads `< 3`: `LOW_MAPPED_READS_lt_3`
- percent coverage `< 10%`: `LOW_PERCENT_COVERAGE_lt_10`
- bacterial 16S reads `< 100`: sample-level `LOW_16S_READS_lt_100`

These thresholds can be adjusted with `--qc-low-mapped-reads`, `--qc-low-percent-coverage`, and `--qc-low-16s-reads`. Both raw and QC-flagged tables are retained.

## Per-sample outputs

```text
results/SAMPLE_ID/
├── 01_read_statistics.tsv
├── 02_CARD_KMA_depth.tsv
├── 02_CARD_KMA_depth_QC.tsv
├── 03_Metaxa2_16S.tsv
├── 04_rrnDB_copy_number.tsv
├── 05_cell_coverage.tsv
├── 06_ARG_copies_per_cell.tsv
├── 06_ARG_copies_per_cell_QC.tsv
├── 07_sample_summary.tsv
├── 07_abundance_by_annotation.tsv
├── validation_comparison.tsv
├── rgi_bwt/
├── metaxa2/
└── logs/
```

`02_CARD_KMA_depth.tsv` retains the CARD model, database, allele source, completely/flanking/all mapped-read counts, coverage metrics, reference length, KMA Depth, and functional annotations. The `06` tables contain the primary per-reference quantification results.

In `07_abundance_by_annotation.tsv`, values are summed directly by ARO term. For gene family, drug class, and resistance mechanism fields containing multiple labels, each reference value is allocated equally among its labels to prevent duplicated abundance.

## SRR6468562 validation example

The repository includes the paired SRR6468562 FASTQ files and compact expected results under `examples/SRR6468562/`. The validation results were generated by rereading the FASTQ files and independently rerunning RGI-bwt/KMA and Metaxa2.

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

The validation command, output, and numerical audit are documented in `workflow.log`. Software and database versions are listed in `environment_versions.txt`.

## Reference

Zhu, C., Wu, L., Ning, D. *et al.* Global diversity and distribution of antibiotic resistance genes in human wastewater treatment systems. *Nature Communications* **16**, 4006 (2025). [https://doi.org/10.1038/s41467-025-59019-3](https://doi.org/10.1038/s41467-025-59019-3)
