# SRR6468562 示例

本示例包含一组去除人源序列后的双端宏基因组数据，以及 CARDCell 生成的精简结果表。

## 输入文件

- `input/SRR6468562.nonhuman_R1.fastq.gz`
- `input/SRR6468562.nonhuman_R2.fastq.gz`

R1 和 R2 各包含 469,775 条序列，合并平均读长为 96.212222872652 bp。

SHA-256 校验值：

```text
99a72d9f47083ba3764bd5f6f7ba567d785d05157260a00985f1e85e24770ddb  SRR6468562.nonhuman_R1.fastq.gz
05ad9f34994335d765ca0521c0e39f9a37cdc3cfcd2159770aac05be5f440766  SRR6468562.nonhuman_R2.fastq.gz
```

## 重新运行

在仓库根目录运行：

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

已发布的精简结果位于 `expected_results/`。大型 SAM/BAM 和临时文件未纳入仓库。

预期核心结果：

```text
细菌 16S 序列数：406
细胞等效覆盖度：5.455958682164
检出的 CARD 参考序列数：15
KMA ARG 总深度：5.600000000000
每细胞 ARG 总拷贝数：1.026400734722
```
