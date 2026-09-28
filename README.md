# CARDCell: CARD/KMA ARG Quantification with 16S-Based Cell Normalization

CARDCell is a standalone, resumable pipeline for quantifying ARG abundance in paired-end metagenomic data.

## Clone the repository

```bash
git clone git@github.com:uby76/CARDCell.git
cd CARDCell
chmod +x bin/cardoap_kma16s
bash scripts/install.sh
```

The only prerequisite for automatic installation is Conda or Mamba. `scripts/install.sh` creates two isolated, repository-local environments and downloads the official databases into `.cardcell/`:

```text
.cardcell/
├── envs/rgi/                 RGI 6.0.8, KMA, samtools, bamtools
├── envs/metaxa2/             Perl, HMMER, BLAST, MAFFT, VSEARCH
├── software/Metaxa2_2.2.3/  Metaxa2 and its bundled SSU database
└── databases/card/localDB/  current official canonical CARD database
```

Two environments are intentional. The Bioconda `metaxa=2.2.3` recipe pins obsolete BLAST/HMMER releases that are not resolvable on current Apple Silicon systems; installing the official Metaxa2 2.2.3 program and database with current compatible dependencies is reproducible across macOS and Linux. The pipeline finds these paths automatically. `.cardcell/` is ignored by Git because generated environments and databases are platform-specific.

The normal total network transfer is below 1 GB: the pinned official Metaxa2 archive is 47.8 MB, while Conda packages and the current CARD canonical archive vary by platform/version. The installer prints this information before downloading. It reuses an existing valid installation and does not overwrite it. To verify all components:

```bash
bash scripts/check_installation.sh
```

To install software first and defer CARD download, use `bash scripts/install.sh --skip-card-db`. To refresh CARD later while preserving the previous database as a timestamped backup, use `bash scripts/install.sh --force-card-db`.

The rrnDB 5.10 pan-taxa NCBI statistics table and the 6,059-entry CARD-OAP class annotation database are already included in the repository. Users do not need to download or reorganize either table.

## Method overview

```text
ARG numerator: FASTQ → RGI-bwt/KMA → CARD → KMA-reported CARD reference Depth
Cell denominator: FASTQ → Metaxa2 bacterial 16S → 16S coverage → taxonomy → rrnDB → cell-equivalent coverage
Final result: ARG copies/cell = KMA Depth / cell-equivalent coverage
ARG class: ARO accession → bundled CARD-OAP database → one Drug_Class per reference
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
  --threads 8
```

By default, the pipeline uses `resources/rrnDB-5.10_pantaxa_stats_NCBI.tsv`. A different current official pan-taxa NCBI statistics table can be supplied with `--rrndb FILE`. CARD is downloaded only by the explicit installation command; analysis runs never download or replace databases.

The pipeline also uses `resources/CARD_OAP_full_6059_annotation.csv` automatically. Its `Class` field directly replaces the original multi-label CARD/RGI `Drug Class` field in all primary result tables. Matching is performed by exact ARO accession. The bundled database contains one unique class for each of 6,059 ARO accessions, so no fractional allocation across multiple drug classes is used. An alternative compatible table can be supplied with `--card-oap-annotation FILE`, but no separate download is required for normal use.

For a custom external installation, override the automatic paths with `--rgi-bin`, `--card-local-db`, `--metaxa-dir`, and `--metaxa-env-bin`. To list all options:

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
├── 08_Drug_Class_abundance.tsv
├── 09_CARD_OAP_annotation_audit.tsv
├── validation_comparison.tsv
├── rgi_bwt/
├── metaxa2/
└── logs/
```

`02_CARD_KMA_depth.tsv` retains the CARD model, database, allele source, completely/flanking/all mapped-read counts, coverage metrics, reference length, KMA Depth, and functional annotations. The `Drug_Class` column is the unique `Class` obtained from the bundled CARD-OAP annotation database, not the original multi-label RGI field. The `06` tables contain the primary per-reference quantification results.

In `07_abundance_by_annotation.tsv`, values are summed directly by ARO term and by the unique bundled `Drug_Class`. Gene-family and resistance-mechanism fields containing multiple labels are still allocated equally among their labels to prevent duplicated abundance. `08_Drug_Class_abundance.tsv` is the concise class-level result. `09_CARD_OAP_annotation_audit.tsv` records the exact ARO match for every detected reference.

## SRR6468562 validation example

The repository includes the paired SRR6468562 FASTQ files and compact expected results under `examples/SRR6468562/`. The validation results were generated by rereading the FASTQ files and independently rerunning RGI-bwt/KMA and Metaxa2.

```bash
bin/cardoap_kma16s \
  --r1 examples/SRR6468562/input/SRR6468562.nonhuman_R1.fastq.gz \
  --r2 examples/SRR6468562/input/SRR6468562.nonhuman_R2.fastq.gz \
  --sample SRR6468562 \
  --outdir results \
  --threads 8
```

The validation command, output, and numerical audit are documented in `workflow.log`. Software and database versions are listed in `environment_versions.txt`.

## Reference

Zhu, C., Wu, L., Ning, D. *et al.* Global diversity and distribution of antibiotic resistance genes in human wastewater treatment systems. *Nature Communications* **16**, 4006 (2025). [https://doi.org/10.1038/s41467-025-59019-3](https://doi.org/10.1038/s41467-025-59019-3)
