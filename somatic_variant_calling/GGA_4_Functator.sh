#!/bin/sh
set -e
SAMPLE1=$1
SAMPLE2=$2
config_path=$3
source ${config_path}/config.sh
echo ${vcf_path}/${SAMPLE1}_${SAMPLE2}_GGA_Filter.vcf.gz
zcat ${vcf_path}/${SAMPLE1}_${SAMPLE2}_GGA_Filter.vcf.gz > ${vcf_path}/${SAMPLE1}_${SAMPLE2}_GGA_Filter.vcf
${gatk4} --java-options "-Xmx10g" Funcotator \
--variant ${vcf_path}/${SAMPLE1}_${SAMPLE2}_GGA_Filter.vcf \
--reference ${seq_ref}/GRCh37.fa \
--ref-version hg19 \
--data-sources-path ${funcotator_ref} \
--output ${maf_path}/${SAMPLE1}_${SAMPLE2}_GGA_Filter_funcotated.maf \
--output-file-format MAF \
--remove-filtered-variants true \
--transcript-selection-mode CANONICAL \
--force-b37-to-hg19-reference-contig-conversion \
2>  ${log_path}/${SAMPLE1}_${SAMPLE2}_GGA_Filter_funcotated.log
