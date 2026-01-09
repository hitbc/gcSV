#!/bin/bash

INPUT=$1
WORK_DIR=$2
PRESET=$3
REF=$4

GCSV_TOOL="./gcSV call "

mkdir ${WORK_DIR}
cd ${WORK_DIR}

BAM_ORI=${INPUT}

VCF=${WORK_DIR}/D.vcf
LOG=${WORK_DIR}/D.log

#call_single_chr
task_call_single_chr() {
    echo "Starting task task_call_single_chr-$1"
    ${GCSV_TOOL} -S $1 -E $1 -p ${PRESET} -l ${BAM_ORI} -r ${REF} -o ${WORK_DIR}/PART_$1_D.vcf 2> /dev/null
    echo "Finished task task_call_single_chr-$1"
}

#combine
task_combine() {
    echo "Starting task combine"
    cat ${WORK_DIR}/PART_0_D.vcf > ${VCF}
    for i in {1..23}; do
        echo ${WORK_DIR}/PART_${i}_D.vcf
        cat ${WORK_DIR}/PART_${i}_D.vcf | grep -v "#" >> ${VCF}
    done
    echo "Finished task combine"
