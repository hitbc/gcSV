# gcSV: 全面的结构变异检测统一框架


## 概述

gcSV 是一个全面的结构变异检测统一框架，支持多种测序数据类型和分析模式，能够准确检测和基因分型基因组中的结构变异。

主要特点：
- 支持纯长读长测序(LRS)数据的结构变异检测和基因分型
- 支持纯短读长测序(SRS)数据的结构变异检测和基因分型
- 支持结合低深度LRS数据和SRS数据进行混合结构变异检测
- 通过局部序列组装精确重建插入序列
- 支持极低深度(深度小于5X)的三代数据结构变异检测
- 支持多种测序平台数据，包括ASM(T2T)、ONT、PacBio、Illumina、BGI-T7等
- 提供面向群体样本复杂结构变异（CSV）分析策略

目录：

* [简介](#简介)
* [安装](#安装)
  * [1. 直接获取静态链接二进制文件并部署](#1-直接获取静态链接二进制文件并部署)
  * [2. 从源代码编译gcSV](#2-从源代码编译gcSV)
* [使用方法](#使用方法)
  * [1. 使用纯 LRS 数据检测和基因分型结构变异](#1-使用纯-lrs-数据检测和基因分型结构变异)
    * [全基因组 SV 检测](#全基因组-sv-检测)
    * [执行1号染色体 SV 检测](#执行1号染色体-sv-检测)
    * [执行特定区域检测 SV（chr2: 11,500,000–18,500,000）](#执行特定区域检测-svchr2-1150000018500000)
  * [2. 使用纯 SRS 数据检测和基因分型结构变异](#2-使用纯-srs-数据检测和基因分型结构变异)
    * [执行全基因组检测](#执行全基因组检测)
    * [特定区域检测 SV（chrX: 11,500,000–18,500,000）](#特定区域检测-svchrX-1150000018500000)
  * [3. 结合低深度 LRS 数据和 SRS 数据进行混合 SV 检测](#3-结合低深度-lrs-数据和-srs-数据进行混合-sv-检测)
  * [4. 多进程处理](#4-多进程处理)
    * [将整个基因组划分为不同区域，同时进行变异检测](#将整个基因组划分为不同区域同时进行变异检测)
    * [将变异检测结果文件合并为一个完整文件](#将变异检测结果文件合并为一个完整文件)
  * [5. LRS 数据分析的高级参数](#5-lrs-数据分析的高级参数)
  * [6. 针对指定区间的局部组装](#6-针对指定区间的局部组装)
    * [生成VCF格式结果](#生成vcf格式结果)
    * [生成BED格式结果](#生成bed格式结果)
  * [7. 面向群体样本的结构变异分析](#7-面向群体样本的结构变异分析)
    * [7.1 群体样本复杂结构变异（CSV）聚类分析](#71-群体样本复杂结构变异csv聚类分析)
    * [7.2 群体样本复杂结构变异区域-MEI-TR阵列检测](#72-群体样本复杂结构变异区域-mei-tr阵列检测)
    * [7.3 群体样本复杂结构变异区域-嵌套变异检测](#73-群体样本复杂结构变异区域-嵌套变异检测)
* [演示](#演示)
* [许可证](#许可证)


## 简介

- gcSV 是一款基于序列比对的结构变异检测工具(以BAM文件为输入)。
- 使用纯 LRS（长读长测序）数据检测和基因分型结构变异
- 使用纯 SRS（短读长测序）数据检测和基因分型结构变异
- 通过结合低深度 LRS 数据和 SRS 数据进行混合结构变异检测和基因分型
- 通过局部序列组装精确重建插入序列
- gcSV支持极低深度（深度小于5X）的三代数据结构变异检测

gcSV支持常见的二代测序平台的数据（Illumina、BGI-T7等）；也支持常见的三代测序平台的测序数据，包括ONT、PacBio的测序数据，或者全局组装结果数据（例如ASM、T2T）。
gcSV在使用高质量三代测序数据（例如ONT-Q26、PacBio HiFi等）表现出较高的检测率和基因分型准确性；同时，gcSV也能分析错误率较高的三代测序数据（错误率约3%-5%的三代测序数据）。

## 安装

gcSV 提供预编译的静态链接二进制文件。使用以下命令进行快速部署和执行：

### 1. 直接获取静态链接二进制文件并部署

```bash
wget https://github.com/hitbc/gcSV/Release/gcSV
./gcSV --help
```

### 2. 从源代码编译gcSV

```bash
git clone https://github.com/hitbc/gcSV/
cd ./Release
make clean
make all -j 12
./gcSV --help
```

## 使用方法

gcSV 支持多种数据的输出，包括基于纯 LRS 数据、基于纯 SRS 数据 或者 结合（低深度） LRS 数据和 SRS 数据进行混合 SV 检测。

另外，可以通过设定待检测的局部区域，来快速获取特定区域的 SV 检测结果。gcSV不提供内置的多线程并行，用户通过将全染色体划分为多个不同的区域，并行运行多个gcSV实例来实现多进程并行。

针对不同类型的三代测序数据，gcSV 提供了不同的预设参数设定，以优化检测和基因分型准确性。"ONT_Q20" 以及 "HIFI" 可以用于错误率低于1%的数据，"ERR_PRONE" 可以用于错误率在1%到6%之间的数据。

### 1. 使用纯 LRS 数据检测和基因分型结构变异

通过设定-l 参数，引入LRS数据或组装结果数据（ASM）进行检测。

#### 全基因组 SV 检测

```bash
gcSV call -l sample.LRS.bam -r ref.fa -o output.vcf 2> /dev/null
```

#### 执行1号染色体 SV 检测

参数-S 和 -E 分别表示起始和结束的染色体 ID，注意染色体 ID 是 0 基的。

```bash
gcSV call -S 0 -E 0 -l sample.LRS.bam -r ref.fa -o output.vcf 2> /dev/null
```

#### 执行特定区域检测 SV（chr2: 11,500,000–18,500,000）

```bash
gcSV call -S 1 -E 1 -s 11500000 -F 18500000 -l sample.LRS.bam -r ref.fa -o output.vcf 2> /dev/null
```

注意：染色体 ID 和位置是 0 基的。默认情况下，执行全基因组 SV 检测。

### 2. 使用纯 SRS 数据检测和基因分型结构变异

在分析前，我们建议使用 `gcSV ngs_fa_stat` 为参考基因组预先计算局部重复复杂性信息。常见人类参考基因组 grch38 和 hg37d5 的预构建结果可在 [GRCh38.stat.txt.gz](https://github.com/hitbc/gcSV/blob/main/demo_gcSV/GRCh38.stat.txt.gz) 和 [hs37d5.stat.txt.gz](https://github.com/hitbc/gcSV/blob/main/demo_gcSV/hs37d5.stat.txt.gz) 获取，使用前需要解压缩。

当进行大规模序列研究或分析完整基因组的时候，使用 `gcSV ngs_trans_reads` 对于SRS数据进行预处理，以减少 I/O 需求并加速分析；该指令不会影响变异检测结果，但是会显著提高分析速度。

#### 执行全基因组检测

```bash
gcSV ngs_fa_stat ref.fa > ref.stat.txt
gcSV ngs_trans_reads ref.fa sample.SRS.bam TL.bam 
samtools sort --output-fmt=BAM -o TL.sort.bam TL.bam
samtools index TL.sort.bam
gcSV call -n sample.SRS.bam -L TL.sort.bam -r ref.fa -I ref.stat.txt -o output.vcf 2> /dev/null
```

#### 特定区域检测 SV（chrX: 11,500,000–18,500,000）

```bash
gcSV call -S 22 -E 22 -s 11500000 -F 18500000 -n sample.SRS.bam -I ref.stat.txt -r ref.fa -o output.vcf 2> /dev/null
```

### 3. 结合低深度 LRS 数据和 SRS 数据进行混合 SV 检测

在一次 gcSV-call 中同时提供 LRS 和 SRS 数据集，以利用它们的互补优势：

```bash
gcSV call -l sample.LRS.bam -n sample.SRS.bam -L TL.sort.bam -r ref.fa -I ref.stat.txt -o output.vcf 2> /dev/null
```

### 4. 多进程处理

gcSV 不直接支持多线程变异检测。可以通过划分染色体的不同区域，并行运行多个gcSV实例来实现多进程并行。

脚本gcSV_multy_process.sh给定了一个示例，展示了如何将整个基因组划分为24个不同的区域，并行运行24个gcSV实例来实现多进程并行变异检测。

```bash
bash gcSV_multy_process.sh sample.LRS.bam ./WORK_DIR ONT_Q26 ref.fa
```

使用如下多进程处理方法：

#### 将整个基因组划分为不同区域，同时进行变异检测

```bash
gcSV call -S 0 -E 0 -l sample.LRS.bam -r ref.fa -o chr1.vcf 2> /dev/null
.....
gcSV call -S 22 -E 22 -H -l sample.LRS.bam -r ref.fa -o chrX.vcf 2> /dev/null
gcSV call -S 23 -E 23 -H -l sample.LRS.bam -r ref.fa -o chrY.vcf 2> /dev/null
```

#### 将变异检测结果文件合并为一个完整文件

```bash
cat chr1.vcf > output.vcf
.....
cat chrX.vcf >> output.vcf
cat chrY.vcf >> output.vcf
```

这种方法适用于所有分析类型（纯 LRS、纯 SRS、混合）。

### 5. LRS 数据分析的高级参数

- `--LRS_preset ONT_Q20` 或 `-p ONT_Q20`：用于预设输入的三代测序数据的类型。默认值是高质量的ONT数据集（ONT_Q20），其他可设定的预设参数包括 "ASM"、"HIFI" 以及 "ERR_PRONE"。"ASM"用于经过组装的结果数据，"ONT_Q20" 以及 "HIFI" 可以用于错误率低于1%的数据，"ERR_PRONE" 可以用于输入错误率在1%到6%之间的数据。

- `--random_phasing 1` 或 `-f 1`：用于对变异进行随机定相（0/1-->0|1 或 1|0 随机；1/1-->1|1），位于同一局部单倍型上的杂合变异将共享相同的随机定相结果。默认开启，设定为 `-f 0` 则关闭随机定相。

- `--Small_var 1` 或 `-v 1`：可以输出与 SV 位于同一局部单倍型上的小变异（SNP 或 INDEL）。默认开启，设定为 `-v 0` 则关闭输出小变异。当分析高错误率数据的时候建议关闭输出小变异。

### 6. 针对指定区间的局部组装

当基于如下参数给出指定区间的时候（-b --FC_BED），gcSV生成特定区间的局部组装序列。该功能仅仅适用于ASM以及LRS数据集，暂时不能应用与SRS数据集，并且每个指定的区间长度不能超过50,000 bp。

#### 生成VCF格式结果

生成的结果存储为VCF格式，contig存储在INFO字段中：

```bash
gcSV call -b region.bed -l sample.LRS.bam -r ref.fa -o output.vcf 2> /dev/null
```

#### 生成BED格式结果

生成的局部contig的结果存储为BED格式，针对完成组装的数据进行分析：

```bash
gcSV call -B BED -b region.bed -l sample.ASM.bam -r ref.fa -o output.bed 2> /dev/null
```

### 7. 面向群体样本的结构变异分析

该功能用于分析特定染色体区间内，群体样本的结构变异结构。特别是为了解决复杂结构变异（CSV）的检测问题而设立。gcSV通过直接分析群体样本的局部单倍型序列信息（局部contig），而非基于vcf，避免了生成VCF的时候，一些细节丢失的问题，从而更好检测和分析群体样本的复杂的结构变异构成。该功能仅仅适用于ASM以及LRS数据集，暂时不能应用与SRS数据集，并且每个指定的区间长度不能超过50,000 bp。

#### 7.1 群体样本复杂结构变异（CSV）聚类分析

大致步骤：

1. 执行全体样本的单样本变异检测；
2. 基于单样本变异检测结果，寻找感兴趣的基因组区间；
3. 使用（-b --FC_BED）参数，在每个感兴趣的基因组区间，对每个样本进行局部组装，生成局部contig；
4. 基于TRF算法与repeat Masker算法对所有局部contig进行contig结构注释；
5. 基于相近congtig的聚类算法，将相近的contig聚在一起，并分析群体样本的复杂的结构变异构成。

完整的工作流程参阅“cohort_csv_analysis.md”。

#### 7.2 群体样本复杂结构变异区域-MEI-TR阵列检测

1. 基于BED文件，对每个样本重新召回局部contig序列；
2. contig 序列 RM、TRF注释；
3. SV-TR阵列检测；

完整的工作流程 参阅“py/科研-MEI阵列-单样本检测.py”。

基于流程7.1, 获取 repeat Masker的变异注释结果： [region].fa.out，然后执行：

```bash
python3 ./py/科研-MEI阵列-单样本检测.py region.fa.out
```

#### 7.3 群体样本复杂结构变异区域-嵌套变异检测

1. 基于BED文件，对每个样本重新召回局部contig序列；
2. contig 序列 RM、TRF注释；
3. 嵌套变异检测；

基于流程7.1的结果文件“nestedStructures.log”, 获取嵌套变异检测结果。

## 演示

演示数据提供了来自 GIAB 的 HG002 真实数据集的一个片段（1:869000-870000），包括演示 Pacbio HiFi 数据集（平均 10x 深度）和演示 Illumina 数据（平均 60x 深度）。

由于演示数据规模较小，无法有效计算 Illumina 数据集的插入片段大小分布。因此，提供了一个额外的文件 `SRS_HG002_stat.json` 来指示 Illumina 数据集的深度和插入片段大小分布统计信息（-T SRS_HG002_stat.json）。在一般使用场景中（例如 WGS），不需要提供此文件或参数。

参数 `-S 0 -E 0 -s 0 -F 1000000` 将变异检测的实际区域限制为 1:1-1000000。

```bash
cd ./demo_gcSV/
# 获取参考分析文件
zcat hs37d5.stat.txt.gz > hs37d5.stat.txt
# LRS SV 检测
gcSV call -S 0 -E 0 -s 0 -F 1000000 -r hs37d5_1_0_1000000.fa -o lrs_demo.vcf -p HIFI -l LRS_HG002_1_869000_870000_10X_demo.bam
# SRS SV 检测
gcSV call -S 0 -E 0 -s 0 -F 1000000 -T SRS_HG002_stat.json -r hs37d5_1_0_1000000.fa -o srs_demo.vcf -I hs37d5.stat.txt -n  SRS_HG002_1_869000_870000_60X_demo.bam
# 混合 SV 检测
gcSV call -S 0 -E 0 -s 0 -F 1000000 -T SRS_HG002_stat.json -r  hs37d5_1_0_1000000.fa -o hybrid_demo.vcf -I hs37d5.stat.txt -n SRS_HG002_1_869000_870000_60X_demo.bam -l LRS_HG002_1_869000_870000_10X_demo.bam -p HIFI
```

## 许可证

该项目受 [GNU GPL v3](LICENSE) 许可证保护。