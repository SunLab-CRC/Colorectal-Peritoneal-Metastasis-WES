#!/bin/sh
source <PROJECT_ROOT>
set -e
SAMPLE1=$1
SAMPLE2=$2
echo ${SAMPLE1}_${SAMPLE2}
${gatk4} --java-options "-Xmx10g" Funcotator \
--variant ${vcf_path}/${SAMPLE1}_${SAMPLE2}_filtered.vcf \
--reference ${seq_ref}/GRCh37.fa \
--ref-version hg19 \
--data-sources-path ${funcotator_ref} \
--output ${vcf_path}/${SAMPLE1}_${SAMPLE2}_funcotated.maf \
--output-file-format MAF \
--remove-filtered-variants true \
--transcript-selection-mode CANONICAL \
--force-b37-to-hg19-reference-contig-conversion \
2>  ${log_path}/${SAMPLE1}_${SAMPLE2}_Funcotator.log
