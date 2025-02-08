# gcSV: a unified framework for comprehensive structural variant detection

## Overview

* [Introduction](#introduction)
* [Installation](#installation)
* [Usage](#usage)
* [License](#license)


## Introduction
- gcSV is a sequence alignment-based structural variation detection tool.
- Detects and genotypes structural variations using pure LRS data
- Detects and genotypes structural variations using pure SRS data
- Performs hybrid structural variation detection and genotyping by combining low-depth LRS data with SRS data
- Precisely reconstructs inserted sequences through local sequence assembly

This project is under active development. Support for specific datasets (e.g., ONT datasets) is currently being implemented. Feedback and suggestions of all kinds are welcome.

##  Installation
gcSV provides pre-compiled statically linked binaries. Use the following command for quick deployment and execution:

### 1: Get statically linked binaries
```bash
> wget https://github.com/hitbc/gcSV/gcSV_release_v1.0.0
> ./gcSV_release_v1.0.0 --help
```

### 2: Compile from source code
```bash
> git clone https://github.com/hitbc/gcSV/
> cd ./Release
> make clean
> make all -j 12
> ./gcSV --help
```

## Usage

### 1. Detects and genotypes structural variations using pure LRS data
Support for ONT datasets is under development; currently only PacBio HiFi datasets are supported.

Whole-genome SV detection

```bash
> gcSV call -c ccs_input.bam -r ref.fa -o output.vcf 2> /dev/null
```

Chromosome 1 SV detection
```bash
> gcSV call -S 0 -E 0 -c ccs_input.bam -r ref.fa -o output.vcf 2> /dev/null
```

Region-specific detection SV (chr2: 11,500,000–18,500,000)
```bash
> gcSV call -S 1 -E 1 -s 11500000 -F 18500000 -c ccs_input.bam -r ref.fa -o output.vcf 2> /dev/null
```
Note: Chromosome IDs and positions are 0-based. By default, whole-genome SV detection is performed.

### 2. Detects and genotypes structural variations using pure SRS data


Before analysis, we recommend precomputing local repeat complexity information using `ngs_fa_stat` for the reference genome. Pre-built results for common human reference genomes grch38 and hg37d5 are available at https://github.com/hitbc/gcSV/ref_stat/.

Preprocess alignment data using `ngs_trans_reads` to reduce I/O demands and accelerate analysis (recommended for WGS or large-scale studies):

```bash
> gcSV ngs_fa_stat ref.fa > ref.stat.txt
> gcSV ngs_trans_reads ref.fa Ill_input.bam TL.bam 
> samtools sort --output-fmt=BAM -o TL.sort.bam TL.bam
> gcSV call -n Ill_input.bam -L TL.sort.bam -r ref.fa -I ref.stat.txt -o output.vcf 2> /dev/null
```

Region-specific detection SV (chrX: 11,500,000–18,500,000)
```bash
> gcSV call -S 22 -E 22 -s 11500000 -F 18500000 -n Ill_input.bam -r ref.fa -o output.vcf 2> /dev/null
```

### 3. Hybrid SV detection by combining low-depth LRS data with SRS data
Provide both LRS and SRS datasets in one gcSV-call to leverage their complementary advantages:

```bash
gcSV call -c ccs_input.bam -n Ill_input.bam -L TL.sort.bam -r ref.fa -I ref.stat.txt -o output.vcf 2> /dev/null
```

4. Multiprocessing
gcSV does not natively support multithreading. For parallel execution, use multiprocessing as follow:

# Divide the whole genome into different regions and perform variant detection simultaneously
```bash
> gcSV call -S 0 -E 0 -c ccs_input.bam -r ref.fa -o chr1.vcf 2> /dev/null
> gcSV call -S 1 -E 1 -H -c ccs_input.bam -r ref.fa -o chr2.vcf 2> /dev/null
> gcSV call -S 2 -E 2 -H -c ccs_input.bam -r ref.fa -o chr3.vcf 2> /dev/null
.....
> gcSV call -S 22 -E 22 -H -c ccs_input.bam -r ref.fa -o chrX.vcf 2> /dev/null
> gcSV call -S 23 -E 23 -H -c ccs_input.bam -r ref.fa -o chrY.vcf 2> /dev/null
```

# Merge the variant detection result files into a complete file.
```bash
> cat chr1.vcf > output.vcf
> cat chr2.vcf >> output.vcf
> cat chr3.vcf >> output.vcf
.....
> cat chrX.vcf >> output.vcf
> cat chrY.vcf >> output.vcf
```

This approach works for all analysis types (pure LRS, pure SRS, hybrid).

5. Advanced parameters for LRS data analysis

Optional flags to optimize output for specific use cases (e.g., contig reconstruction):


 `--random_phasing` or `-f`, Randomly phasing all SVs randomly (0/1-->0|1 or 1|0; 1/1-->1|1)
 `--Small_var` or  `-v` 

The `--random_phasing` or `-f` parameter is used to perform random phasing (0/1-->0|1 or 1|0 randomly; 1/1-->1|1) on variants, and heterozygous variants located on the same local haplotype will share the same random phasing result.
The `--Small_var` or  `-v`  parameter can output small variants (SNPs or INDELs) that are on the same local haplotype as the SV.

## License
This project is covered under the <a href="LICENSE">GNU GPL v3</a> license.

