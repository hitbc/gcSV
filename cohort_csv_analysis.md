# Complete Workflow for Cohort CSV Analysis

## Introduction

This document introduces a complete workflow for cohort CSV (Complex Structural Variation) analysis using the gcSV tool. This workflow is suitable for Long Read Sequencing (LRS) data and implements a complete process from single-sample variant detection to population CSV region analysis, filtering, precise detection, filtering, and clustering analysis through a series of steps.

The workflow mainly includes the following core functions:
- Preliminary variant detection at the single-sample level
- CSV region identification based on population data
- Screening of medium and high-frequency CSV regions
- High-precision Force Calling variant detection
- Filtering and storage of Contig sequences
- Contig clustering analysis based on RM-TRF annotation

This workflow is suitable for structural variation analysis of large-scale population samples and can help researchers efficiently identify structurally variant regions of important biological significance in the population.

## Table of Contents

1. [Workflow Overview](#workflow-overview)
2. [Data Preparation](#data-preparation-before-starting)
3. [Single-Sample Variant Detection](#1-single-sample-variant-detection)
   - [Create Single-Sample Variant Detection Script](#11-create-single-sample-variant-detection-script-calldefaultsh)
   - [Batch Generate Single-Sample Variant Detection Tasks](#12-batch-generate-and-submit-single-sample-variant-detection-tasks)
4. [Single-Sample CSV Detection and Population CSV Region Detection](#2-single-sample-csv-detection-and-population-csv-region-detection)
   - [Single-Sample CSV Region Analysis](#21-single-sample-csv-region-analysis)
   - [Batch Submit Single-Sample CSV Region Analysis Jobs](#22-batch-submit-single-sample-csv-region-analysis-jobs)
   - [Merge Population CSV Region Analysis Results](#23-merge-population-csv-region-analysis-results)
5. [Screening Medium and High-Frequency Population CSV](#3-screening-medium-and-high-frequency-population-csv-generate-bed-file)
   - [Screen CSV Regions with Medium and High Allele Counts](#31-screen-csv-regions-with-medium-and-high-allele-counts-ac≥56)
   - [Sort and Merge BED Files](#32-sort-and-merge-bed-files)
6. [Recall Local Contig Sequences for Each Sample Based on BED File](#4-recall-local-contig-sequences-for-each-sample-based-on-bed-file)
   - [Create Single-Sample Force Calling Script](#41-create-single-sample-fc-script-singlesamplesh)
   - [Batch Submit Force Calling Jobs](#42-batch-submit-force-calling-jobs)
7. [Contig Filtering](#5-contig-filtering)
   - [Create Filtering Script](#51-create-filtering-script-singlesamplefiltersh)
   - [Batch Submit Filtering Jobs](#52-batch-submit-filtering-jobs)
8. [Store Variants in Corresponding Regions](#6-store-variants-in-corresponding-regions)
   - [Prepare Input File List](#61-prepare-input-file-list)
   - [Execute Variant Storage](#62-execute-variant-storage)
9. [Contig Clustering Analysis](#7-contig-clustering-analysis)
   - [Create Annotation and Analysis Script](#71-create-annotation-and-analysis-script-runsingleannosh)
   - [Batch Submit Annotation Jobs](#72-batch-submit-annotation-jobs)

## Workflow Overview

1. **Single-Sample Variant Detection**: Perform preliminary variant detection for each sample, generate VCF files
2. **Population CSV Region Detection**: Analyze VCF files of all samples to identify potential CSV regions
3. **Screening Medium and High-Frequency Population CSV**: Filter out CSV regions with sufficient sample support based on allele frequency
4. **Force Calling**: Perform high-precision variant detection again for each sample based on the filtered CSV regions, generate contig sequences
5. **Contig Filtering**: Filter out contigs with insufficient read support
6. **Variant Region Storage**: Classify and store filtered contigs according to regions
7. **Contig Clustering Analysis**: Perform RM-TRF annotation and further clustering analysis on contigs

## Data Preparation Before Starting:

BAM.list: A file describing the storage paths of BAM results for all samples, and storing sample names.
PRESET: Describes the type of LRS sequencing (ONT_Q26, ASM, ERR_PRONE, etc.)

## 1. Single-Sample Variant Detection

### 1.1 Create Single-Sample Variant Detection Script [call_default.sh]

```bash
#!/bin/bash
INPUT=$1
WORK_DIR=$2
PRESET=$3

#Fixed parameters
GCSV_TOOL="gcSVL call "
REF=GRCh38.fa

mkdir ${WORK_DIR}
cd ${WORK_DIR}

BAM_ORI=${INPUT}

VCF=${WORK_DIR}/D.vcf
LOG=${WORK_DIR}/D.log

${GCSV_TOOL} -p ${PRESET} -l ${BAM_ORI} -r ${REF} -o ${VCF} 2> ${LOG}
```

### 1.2 Batch Generate and Submit Single-Sample Variant Detection Tasks
```bash
cat BAM.list | awk '{printf "sbatch %s ./RST/%s ERR_PRONE \n",$1,$2}' > run_all_SV_calling.sh
bash run_all_SV_calling.sh
```

## 2. Single-Sample CSV Detection and Population CSV Region Detection

### 2.1 Single-Sample CSV Region Analysis
```bash
#!/bin/bash

#Parameters for different samples
SAMPLE_NAME=$1

IN_VCF=./${SAMPLE_NAME}/D.vcf
WORK_DIR=./${SAMPLE_NAME}/

#Fixed parameters
GCSV_TOOL="gcSV "
REF=GRCh38.fa

mkdir ${WORK_DIR}
cd ${WORK_DIR}

LOG=${WORK_DIR}/D.log

${GCSV_TOOL} tools csv_region_500w ${IN_VCF}  500 50 > ${SAMPLE_NAME}.main.bed 2> ${LOG}

cat ${SAMPLE_NAME}.main.bed | grep "HAP1" | grep "500W" | sort | uniq > ${SAMPLE_NAME}.hap1.bed
cat ${SAMPLE_NAME}.main.bed | grep "HAP2" | grep "500W" | sort | uniq > ${SAMPLE_NAME}.hap2.bed
```

### 2.2 Batch Submit Single-Sample CSV Region Analysis Jobs

```bash
mkdir ./RST_ANA_CSV/
cat BAM.list | awk '{printf "sbatch csv_region_500w_ana.sh %s \n",$2}' > run_all_csv_region_500w_ana.sh
bash run_all_csv_region_500w_ana.sh > run_all_csv_region_500w_ana.bash.log &
```

### 2.3 Merge Population CSV Region Analysis Results

```bash
cd ./RST_ANA_CSV/
cat */*hap*bed > ALL_single_hap.bed
cat ALL_single_hap.bed | awk '{printf "%s\n",$1;}'
cat ALL_single_hap.bed | awk '{printf "%s:%s-%s\n",$1,$2,$3;}' | sort | uniq -c > ALL_single_hap_summary.txt
```

## 3. Screening Medium and High-Frequency Population CSV, Generate BED File

### 3.1 Screen CSV Regions with Medium and High Allele Counts (AC≥56)

```bash
cat ALL_single_hap_summary.txt  | awk '{split($2, A, ":"); split(A[2], B, "-"); if($1 >=56) printf "%s\t%s\t%s\t%s\n", A[1],B[1],B[2],$1}' > AF_0.05_CSV_AF.bed
```

### 3.2 Sort and Merge BED Files
```bash
cd ./RST_ANA_CSV/
bedtools sort -i AF_0.05_CSV_AF.bed > AF_0.05_CSV_AF.sort.bed
bedtools merge -i AF_0.05_CSV_AF.sort.bed > AF_0.05_CSV_AF.sort.merge.bed
```

## 4. Recall Local Contig Sequences for Each Sample Based on BED File

### 4.1 Create Single-Sample FC Script (single_sample.sh)
```bash
#!/bin/bash

#Parameters
BAM_FILE=$1
OUT_DIR=$2
BED_FILE=$3
PRESET=$4

#Fixed parameters
GCSV_TOOL="gcSVL call "
REF=GRCh38.fa

mkdir -p ${OUT_DIR}
cd ${OUT_DIR}

${GCSV_TOOL} -p ${PRESET} -b ${BED_FILE} -l ${BAM_FILE} -r ${REF} -o ./gcSV_FC_CSV_RST.txt 2> ./gcSV_FC_CSV_RST.log
```

### 4.2 Batch Submit Force Calling Jobs
```bash
cd ./
cat BAM.list | awk '{printf "sbatch ./single_sample.sh %s ./RST/%s/ ./AF_0.05_CSV_AF.sort.merge.bed ERR_PRONE \n",$1,$2;}' >  ./single_FC_gcSV_run_ALL.sh
bash ./single_FC_gcSV_run_ALL.sh > ./single_FC_gcSV_run_ALL.log &

```

## 5. Contig Filtering

Filter as needed, here we filtered out contigs with less than 2 supporting reads

### 5.1 Create Filtering Script (single_sample_FILTER.sh)
```bash
#!/bin/bash

#Parameters
IN_DIR=$1

cd ${IN_DIR}

#Delete results with read support equal to 1
cat gcSV_FC_CSV_RST.txt | awk '{split($4,A,"="); split(A[2],B,":"); if(B[1] > 1) print $0}' > gcSV_FC_CSV_RST.FILTER.txt
```

### 5.2 Batch Submit Filtering Jobs
```bash
cd ./
cat BAM.list | awk '{printf "sbatch ./single_sample_FILTER.sh ./RST/%s/\n",$2;}' > ./single_sample_FILTER_run_ALL.sh
bash ./single_sample_FILTER_run_ALL.sh > ./single_sample_FILTER_run_ALL.log
```

## 6. Store Variants in Corresponding Regions

### 6.1 Prepare Input File List
```bash
cd ./
cat BAM.list | awk '{printf "./RST/%s/gcSV_FC_CSV_RST.FILTER.txt\n",$2}' > ont_contig_string_data_ALL_FL.txt
```

### 6.2 Execute Variant Storage
```bash
gcSVL tools contig_file_split ont_contig_string_data_ALL_FL.txt ./CONTIG_IN_REGION/ GRCh38.fa
```

## 7. Contig Clustering Analysis

### 7.1 Create Annotation and Analysis Script (run_single_ANNO.sh)

```bash
#!/bin/bash

PREFIX=$1
WORD_DIR=$2
FC_DIR=$3

GCSV_TOOL="gcSVL call "
REF_38=GRCh38.fa
RM_TOOL=RepeatMasker
TRF_TOOL=trf

cd ${WORD_DIR}
mkdir ${PREFIX}
cd ${PREFIX}
cp ${FC_DIR}/${PREFIX} ./${PREFIX}.fa

#Execute RM annotation:
${RM_TOOL} ${PREFIX}.fa > rm_log.txt 2>  rm_log2.txt
#Execute TRF annotation:
${TRF_TOOL} ${PREFIX}.fa 2 7 7 80 10 10 500 -h -d -ngs 1> ${PREFIX}.fa.trf.out 2>  trf_log2.txt
#Execute population CSV clustering analysis:
${GCSV_TOOL} tools fc_joint_ana ${REGION_ID} .
```

#### 7.2 Batch Submit Annotation Jobs
```bash
cd ./RST/

cat AF_0.05_CSV_AF.sort.merge.bed | awk '{printf "sbatch run_single_ANNO.sh  %s_%s_%s ONT_567/RST/ ./CONTIG_IN_REGION/ \n",$6,$7,$8;}' > ONT_567/run_all.sh
bash ONT_567/run_all.sh > ONT_567/run_all.log &
```