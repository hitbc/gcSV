# gcSV: A Comprehensive Unified Framework for Structural Variant Detection

## Overview

gcSV is a comprehensive unified framework for structural variant detection that supports multiple sequencing data types and analysis modes, enabling accurate detection and genotyping of structural variants in genomes.

Key Features:
- Supports structural variant detection and genotyping using pure long-read sequencing (LRS) data
- Supports structural variant detection and genotyping using pure short-read sequencing (SRS) data
- Supports hybrid structural variant detection by combining low-depth LRS data with SRS data
- Precisely reconstructs inserted sequences through local sequence assembly
- Supports structural variant detection from ultra-low depth (less than 5X) third-generation sequencing data
- Supports multiple sequencing platforms, including ASM (T2T), ONT, PacBio, Illumina, BGI-T7, etc.
- Provides strategies for complex structural variant (CSV) analysis of population samples

Table of Contents:

* [Introduction](#introduction)
* [Installation](#installation)
  * [1. Directly Obtain and Deploy Precompiled Static Binary](#1-directly-obtain-and-deploy-precompiled-static-binary)
  * [2. Compile gcSV from Source Code](#2-compile-gcsv-from-source-code)
* [Usage](#usage)
  * [1. Detect and Genotype Structural Variants Using Pure LRS Data](#1-detect-and-genotype-structural-variants-using-pure-lrs-data)
    * [Whole Genome SV Detection](#whole-genome-sv-detection)
    * [Chromosome 1 SV Detection](#chromosome-1-sv-detection)
    * [Specific Region SV Detection (chr2: 11,500,000–18,500,000)](#specific-region-sv-detection-chr2-1150000018500000)
  * [2. Detect and Genotype Structural Variants Using Pure SRS Data](#2-detect-and-genotype-structural-variants-using-pure-srs-data)
    * [Whole Genome Detection](#whole-genome-detection)
    * [Specific Region SV Detection (chrX: 11,500,000–18,500,000)](#specific-region-sv-detection-chrx-1150000018500000)
  * [3. Hybrid SV Detection Combining Low-depth LRS Data and SRS Data](#3-hybrid-sv-detection-combining-low-depth-lrs-data-and-srs-data)
  * [4. Multi-process Processing](#4-multi-process-processing)
    * [Divide the Entire Genome into Different Regions for Simultaneous Variant Detection](#divide-the-entire-genome-into-different-regions-for-simultaneous-variant-detection)
    * [Merge Variant Detection Result Files into a Complete File](#merge-variant-detection-result-files-into-a-complete-file)
  * [5. Advanced Parameters for LRS Data Analysis](#5-advanced-parameters-for-lrs-data-analysis)
  * [6. Local Assembly for Specified Intervals](#6-local-assembly-for-specified-intervals)
    * [Generate VCF Format Results](#generate-vcf-format-results)
    * [Generate BED Format Results](#generate-bed-format-results)
  * [7. Structural Variant Analysis for Population Samples](#7-structural-variant-analysis-for-population-samples)
    * [7.1 Population Sample Complex Structural Variant (CSV) Clustering Analysis](#71-population-sample-complex-structural-variant-csv-clustering-analysis)
    * [7.2 Population Sample Complex Structural Variant Region - MEI-TR Array Detection](#72-population-sample-complex-structural-variant-region-mei-tr-array-detection)
    * [7.3 Population Sample Complex Structural Variant Region - Nested Variant Detection](#73-population-sample-complex-structural-variant-region-nested-variant-detection)
* [Demo](#demo)
* [License](#license)


## Introduction

- gcSV is a sequence alignment-based structural variant detection tool (taking BAM files as input).
- Detects and genotypes structural variants using pure LRS (long-read sequencing) data or assembly result data (ASM)
- Detects and genotypes structural variants using pure SRS (short-read sequencing) data
- Performs hybrid structural variant detection and genotyping by combining (low-depth) LRS data with SRS data
- Precisely reconstructs inserted sequences through local sequence assembly
- Supports structural variant detection from ultra-low depth (less than 5X) third-generation sequencing data

gcSV supports data from common second-generation sequencing platforms (Illumina, BGI-T7, etc.) as well as data from common third-generation sequencing platforms, including ONT, PacBio sequencing data, or global assembly result data (such as ASM, T2T).
gcSV demonstrates high detection rates and genotyping accuracy when using high-quality third-generation sequencing data (such as ONT-Q26, PacBio HiFi, etc.); at the same time, gcSV can also analyze third-generation sequencing data with higher error rates (error rates of approximately 3%-5%).

## Installation
gcSV provides precompiled static linked binary files. Use the following commands for quick deployment and execution:

### 1. Directly Obtain and Deploy Precompiled Static Binary

```bash
wget https://github.com/hitbc/gcSV/Release/gcSV
./gcSV --help
```

### 2. Compile gcSV from Source Code

```bash
git clone https://github.com/hitbc/gcSV/
cd ./Release
make clean
make all -j 12
./gcSV --help
```

## Usage

gcSV supports multiple data inputs for SV detection, including pure LRS data, pure SRS data, or a combination of (low-depth) LRS data and SRS data.

In addition, you can quickly obtain SV detection results for specific regions by setting the local region to be detected. gcSV does not provide built-in multi-threaded parallelization; users can divide entire chromosomes into multiple different regions and run multiple gcSV instances in parallel to achieve multi-process parallelization.

For different types of third-generation sequencing data, gcSV provides different preset parameter settings to optimize detection and genotyping accuracy. "ONT_Q20" and "HIFI" can be used for data with error rates below 1%, and "ERR_PRONE" can be used for data with error rates between 1% and 6%.

### 1. Detect and Genotype Structural Variants Using Pure LRS Data

Introduce LRS data or assembly result data (ASM) for detection by setting the -l parameter.

#### Whole Genome SV Detection

```bash
gcSV call -l sample.LRS.bam -r ref.fa -o output.vcf 2> /dev/null
```

#### Chromosome 1 SV Detection

Parameters -S and -E represent the start and end chromosome IDs, respectively. Note that chromosome IDs are 0-based.

```bash
gcSV call -S 0 -E 0 -l sample.LRS.bam -r ref.fa -o output.vcf 2> /dev/null
```

#### Specific Region SV Detection (chr2: 11,500,000–18,500,000)

```bash
gcSV call -S 1 -E 1 -s 11500000 -F 18500000 -l sample.LRS.bam -r ref.fa -o output.vcf 2> /dev/null
```

Note: Chromosome IDs and positions are 0-based. By default, whole-genome SV detection is performed.

### 2. Detect and Genotype Structural Variants Using Pure SRS Data

Before analysis, we recommend using `gcSV ngs_fa_stat` to pre-calculate local repeat complexity information for the reference genome. Pre-built results for common human reference genomes grch38 and hg37d5 are available at [GRCh38.stat.txt.gz](https://github.com/hitbc/gcSV/blob/main/demo_gcSV/GRCh38.stat.txt.gz) and [hs37d5.stat.txt.gz](https://github.com/hitbc/gcSV/blob/main/demo_gcSV/hs37d5.stat.txt.gz), which need to be decompressed before use.

When conducting large-scale sequence research or analyzing complete genomes, use `gcSV ngs_trans_reads` to preprocess SRS data to reduce I/O requirements and accelerate analysis; this command does not affect variant detection results but significantly improves analysis speed.

#### Whole Genome Detection

```bash
gcSV ngs_fa_stat ref.fa > ref.stat.txt
gcSV ngs_trans_reads ref.fa sample.SRS.bam TL.bam 
samtools sort --output-fmt=BAM -o TL.sort.bam TL.bam
samtools index TL.sort.bam
gcSV call -n sample.SRS.bam -L TL.sort.bam -r ref.fa -I ref.stat.txt -o output.vcf 2> /dev/null
```

#### Specific Region SV Detection (chrX: 11,500,000–18,500,000)

```bash
gcSV call -S 22 -E 22 -s 11500000 -F 18500000 -n sample.SRS.bam -I ref.stat.txt -r ref.fa -o output.vcf 2> /dev/null
```

### 3. Hybrid SV Detection Combining Low-depth LRS Data and SRS Data

Provide both LRS and SRS datasets in a single gcSV-call to leverage their complementary advantages:

```bash
gcSV call -l sample.LRS.bam -n sample.SRS.bam -L TL.sort.bam -r ref.fa -I ref.stat.txt -o output.vcf 2> /dev/null
```

### 4. Multi-process Processing

gcSV does not directly support multi-threaded variant detection. Multi-process parallelization can be achieved by dividing chromosomes into different regions and running multiple gcSV instances in parallel.

The script gcSV_multy_process.sh provides an example showing how to divide the entire genome into 24 different regions and run 24 gcSV instances in parallel to achieve multi-process parallel variant detection.

```bash
bash gcSV_multy_process.sh sample.LRS.bam ./WORK_DIR ONT_Q26 ref.fa
```

Using the following multi-process processing method:

#### Divide the Entire Genome into Different Regions for Simultaneous Variant Detection

```bash
gcSV call -S 0 -E 0 -l sample.LRS.bam -r ref.fa -o chr1.vcf 2> /dev/null
.....
gcSV call -S 22 -E 22 -H -l sample.LRS.bam -r ref.fa -o chrX.vcf 2> /dev/null
gcSV call -S 23 -E 23 -H -l sample.LRS.bam -r ref.fa -o chrY.vcf 2> /dev/null
```

#### Merge Variant Detection Result Files into a Complete File

```bash
cat chr1.vcf > output.vcf
.....
cat chrX.vcf >> output.vcf
cat chrY.vcf >> output.vcf
```

This method is applicable to all analysis types (pure LRS, pure SRS, hybrid).

### 5. Advanced Parameters for LRS Data Analysis

- `--LRS_preset ONT_Q20` or `-p ONT_Q20`: Used to preset the type of input third-generation sequencing data. The default value is high-quality ONT dataset (ONT_Q20), and other available preset parameters include "ASM", "HIFI", and "ERR_PRONE". "ASM" is used for assembled result data, "ONT_Q20" and "HIFI" can be used for data with error rates below 1%, and "ERR_PRONE" can be used for data with error rates between 1% and 6%.

- `--random_phasing 1` or `-f 1`: Used for random phasing of variants (0/1 --> 0|1 or 1|0 randomly; 1/1 --> 1|1). Heterozygous variants located on the same local haplotype will share the same random phasing result. It is enabled by default; setting it to `-f 0` will disable random phasing.

- `--Small_var 1` or `-v 1`: Can output small variants (SNPs or INDELs) located on the same local haplotype as SVs. It is enabled by default; setting it to `-v 0` will disable the output of small variants. It is recommended to disable the output of small variants when analyzing high-error-rate data.

### 6. Local Assembly for Specified Intervals

When specifying intervals with the following parameters (-b --FC_BED), gcSV generates local assembly sequences for specific intervals. This function only applies to ASM and LRS datasets, not to SRS datasets, and the length of each specified interval cannot exceed 50,000 bp.

#### Generate VCF Format Results

The generated results are stored in VCF format, with contigs stored in the INFO field:

```bash
gcSV call -b region.bed -l sample.LRS.bam -r ref.fa -o output.vcf 2> /dev/null
```

#### Generate BED Format Results

The generated local contig results are stored in BED format for analysis of completed assembly data:

```bash
gcSV call -B BED -b region.bed -l sample.ASM.bam -r ref.fa -o output.bed 2> /dev/null
```

### 7. Structural Variant Analysis for Population Samples

This function is used to analyze the structural variant structure of population samples within specific chromosomal intervals. It is particularly designed to address the detection of complex structural variants (CSV). gcSV directly analyzes the local haplotype sequence information (local contigs) of population samples instead of being based on VCF, avoiding the loss of some details when generating VCF, thereby better detecting and analyzing the complex structural variant composition of population samples. This function only applies to ASM and LRS datasets, not to SRS datasets, and the length of each specified interval cannot exceed 50,000 bp.

#### 7.1 Population Sample Complex Structural Variant (CSV) Clustering Analysis

General steps:

1. Perform single-sample variant detection for all samples;
2. Based on the single-sample variant detection results, identify genomic intervals of interest;
3. Use the (-b --FC_BED) parameter to perform local assembly on each sample in each genomic interval of interest to generate local contigs;
4. Annotate the structure of all local contigs using the TRF algorithm and repeat Masker algorithm;
5. Based on the clustering algorithm for similar contigs, cluster similar contigs together and analyze the complex structural variant composition of population samples.

Refer to "cohort_csv_analysis.md" for the complete workflow.

#### 7.2 Population Sample Complex Structural Variant Region - MEI-TR Array Detection

1. Re-call local contig sequences for each sample based on BED files;
2. Contig sequence RM and TRF annotation;
3. SV-TR array detection;

Refer to "py/科研-MEI阵列-单样本检测.py" (科研-MEI阵列-单样本检测.py translates to: Research-MEI Array-Single Sample Detection.py) for the complete workflow.

Based on workflow 7.1, obtain the variant annotation results from repeat Masker: [region].fa.out, then execute:

```bash
python3 ./py/科研-MEI阵列-单样本检测.py region.fa.out
```

#### 7.3 Population Sample Complex Structural Variant Region - Nested Variant Detection

1. Re-call local contig sequences for each sample based on BED files;
2. Contig sequence RM and TRF annotation;
3. Nested variant detection;

Obtain nested variant detection results based on the result file "nestedStructures.log" from workflow 7.1.

## Demo

Demo data provides a fragment (1:869000-870000) from GIAB's HG002 real dataset, including demo Pacbio HiFi dataset (average 10x depth) and demo Illumina data (average 60x depth).

Due to the small size of the demo data, the insertion fragment size distribution of the Illumina dataset cannot be effectively calculated. Therefore, an additional file `SRS_HG002_stat.json` is provided to indicate the depth and insertion fragment size distribution statistics of the Illumina dataset (-T SRS_HG002_stat.json). In general usage scenarios (such as WGS), there is no need to provide this file or parameter.

The parameters `-S 0 -E 0 -s 0 -F 1000000` limit the actual region for variant detection to 1:1-1000000.

```bash
cd ./demo_gcSV/
# Obtain reference analysis files
zcat hs37d5.stat.txt.gz > hs37d5.stat.txt
# LRS SV detection
gcSV call -S 0 -E 0 -s 0 -F 1000000 -r hs37d5_1_0_1000000.fa -o lrs_demo.vcf -p HIFI -l LRS_HG002_1_869000_870000_10X_demo.bam
# SRS SV detection
gcSV call -S 0 -E 0 -s 0 -F 1000000 -T SRS_HG002_stat.json -r hs37d5_1_0_1000000.fa -o srs_demo.vcf -I hs37d5.stat.txt -n  SRS_HG002_1_869000_870000_60X_demo.bam
# Hybrid SV detection
gcSV call -S 0 -E 0 -s 0 -F 1000000 -T SRS_HG002_stat.json -r  hs37d5_1_0_1000000.fa -o hybrid_demo.vcf -I hs37d5.stat.txt -n SRS_HG002_1_869000_870000_60X_demo.bam -l LRS_HG002_1_869000_870000_10X_demo.bam -p HIFI
```

## License

This project is protected by the [GNU GPL v3](LICENSE) license.
