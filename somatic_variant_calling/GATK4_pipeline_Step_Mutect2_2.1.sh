#!/bin/sh
sample=$1
config_path=$2
source ${config_path}/config.sh
set -e
${gatk4} --java-options "-Xmx20g" GetPileupSummaries \
-I ${bam_path}/${sample}-recal.bam \
-V ${variants_ref}/af-only-gnomad.raw.sites.b37.vcf.gz \
-L ${target_bed} \
-O ${vcf_path}/${sample}_getpileupsummaries.table \
2>  ${log_path}/${sample}_getpileupsummaries.log
