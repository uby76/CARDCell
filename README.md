# CARDCell：基于 CARD/KMA 与 16S 的 ARG 细胞归一化流程

CARDCell 是一套独立、可续跑的双端宏基因组定量流程。

## 获取仓库

```bash
git clone git@github.com:uby76/CARDCell.git
cd CARDCell
chmod +x bin/cardoap_kma16s
```

运行前需要：Python 3.9 或更高版本、RGI（含 KMA、samtools、bamtools 运行依赖）、已加载的 CARD 标准本地数据库、Metaxa2 2.2.3 及其 SSU 数据库。仓库已包含 rrnDB 5.10 泛分类群 NCBI 统计表；CARD 和 Metaxa2 数据库不会自动下载或替换。

## 方法一句话说明

```text
ARG 分子：FASTQ → RGI-bwt/KMA → CARD → KMA 报告的 CARD 参考序列深度
细胞分母：FASTQ → Metaxa2 细菌 16S → 16S 覆盖度 → 分类组成 → rrnDB → 细胞等效覆盖度
最终结果：每细胞 ARG 拷贝数 = KMA 深度 / 细胞等效覆盖度
```

本流程把 **CARD/KMA ARG 分析**与 **Zhu 等（2025）风格的 16S + rrnDB 细菌细胞归一化**结合。不要把整套流程称为 Zhu 等（2025）的原始流程；只有细胞分母采用该文描述的 16S 归一化思想。

## 核心公式

实际读长由输入 FASTQ 逐条统计，不预设 100 或 150 bp。

```text
16S_coverage = N_16S_total × combined_mean_read_length / 1432

weighted_mean_16S_copy_number
  = Σ(relative_abundance_i × rrn_copy_number_i)

cell_equivalent_coverage
  = 16S_coverage / weighted_mean_16S_copy_number

ARG_copies_per_cell
  = KMA-reported_CARD_reference_Depth / cell_equivalent_coverage
```

`1432 bp` 是 Zhu 等（2025）使用的平均细菌 16S rRNA 基因长度。ARG 分子只取 RGI-bwt/KMA 输出 `allele_mapping_data.txt` 的 `Depth` 字段。本流程将其明确称为 **KMA 报告的 CARD 参考序列深度**。

以下项目不参与主计算：BAM 计算的平均深度、ARGs-OAP KO30 nCell、实验性 KMA-KO30 数据库、组装、ORF 预测、DIAMOND、BLASTX、DeepARG、Kraken2/Bracken。RGI 内部会生成 BAM 并调用 samtools/bamtools，但本流程的计算脚本不读取 BAM。

## rrnDB 匹配规则

Metaxa2 在 R1、R2 上分别识别细菌 SSU/16S 序列。群落相对丰度以全部细菌 16S 序列为分母：可鉴定到属时按属汇总；不能鉴定到属时保留最深的可靠分类层级。

rrnDB 使用第一个精确的 rank/name 匹配，顺序为：

```text
属 → 科 → 目 → 纲 → 门 → 域
```

每个回退层级写入 `04_rrnDB_copy_number.tsv`。不使用固定 4 或 4.1 替代缺失值。若存在未匹配丰度，只用匹配部分重新归一化来计算加权均值，并在结果中报告匹配/未匹配比例及是否发生重新归一化。

## 运行方法

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

默认使用本仓库的 `resources/rrnDB-5.10_pantaxa_stats_NCBI.tsv`。也可用 `--rrndb FILE` 指定当前官方 pan-taxa NCBI 统计表。CARD 不会自动下载或替换。

RGI、Metaxa2 或数据库不在 PATH/默认位置时，分别用 `--rgi-bin`、`--card-local-db`、`--metaxa-dir` 和 `--metaxa-env-bin` 指定。查看全部参数：

```bash
bin/cardoap_kma16s --help
```

常用控制参数：

- `--force`：只重建指定样本在本输出目录中的结果，不触碰输入或旧流程。
- `--keep-temp`：保留解压后的临时 FASTQ。
- 无 `--force` 时，如果 RGI 和 Metaxa2 的必需输出完整且非空，则自动续跑。
- `--comparison-ko30`、`--comparison-previous-16s`、`--comparison-previous-16s-count`：只写验证比较表，绝不进入主公式。

## QC 规则

默认阈值均为透明标记，不删除任何 ARG：

- 比对序列数 `< 3`：`LOW_MAPPED_READS_lt_3`
- 覆盖百分比 `< 10%`：`LOW_PERCENT_COVERAGE_lt_10`
- 细菌 16S 序列数 `< 100`：样本级标记 `LOW_16S_READS_lt_100`

阈值分别可用 `--qc-low-mapped-reads`、`--qc-low-percent-coverage`、`--qc-low-16s-reads` 调整。原始表与 QC 标记表同时保留。

## 每个样本的输出

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

`02_CARD_KMA_depth.tsv` 保留 CARD 模型、数据库、等位基因来源、完整/带侧翼/全部比对序列数、覆盖度、参考序列长度、KMA 深度和功能注释。`06` 是逐 CARD 参考序列的主定量结果。

`07_abundance_by_annotation.tsv` 对 ARO 条目直接求和；含多个标签的基因家族、药物类别或耐药机制采用**等比例分配**，避免同一参考序列因多个标签而被重复累计。

## SRR6468562 验证

仓库提供可公开复现的 SRR6468562 双端 FASTQ 和结果实例，见 `examples/SRR6468562/`。测试结果来自重新读取 FASTQ、重新运行 RGI-bwt/KMA 和重新运行 Metaxa2。

示例命令：

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

验证命令、运行输出和数值审计见 `workflow.log`；软件及数据库版本见 `environment_versions.txt`。
