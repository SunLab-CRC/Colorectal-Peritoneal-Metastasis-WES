# Colorectal Cancer Peritoneal Metastasis WES

This repository contains the core analysis scripts for whole-exome sequencing of colorectal cancer with peritoneal metastasis.

## Repository structure

* `somatic_variant_calling/`: read trimming, BWA alignment, BQSR, Mutect2 calling, filtering, GGA cross-sample recall, and Funcotator annotation.
* `analysis/sample_information/`: sample-quality summary.
* `analysis/global_genomic_features/`: tumour purity, ploidy, copy-number burden, and mutation burden analyses.
* `analysis/mutation_landscape/`: MAF preparation, driver-gene analysis, mutation landscape, and shared-mutation analyses.
* `analysis/metastatic_driver_genes/`: paired and external-cohort analyses of somatic mutation and copy-number driver genes. `CNV_paired_test/` contains helper functions for consensus copy-number events.
* `analysis/drug_targets/`: somatic mutation and copy-number analyses of clinically actionable targets.

## Somatic variant calling

The calling workflow uses GRCh37/hg19 and is based on Trimmomatic, BWA-MEM, GATK 4, Picard, samtools, bcftools, bedtools, bgzip/tabix, and GATK Funcotator. The scripts require a local `config.sh` that defines tool locations, reference files, input paths, output paths, and compute settings.

Run the scripts in the following order:

1. `GATK4_pipeline_trim.sh`
2. `GATK4_pipeline_bwamap.sh.sh`
3. `GATK4_pipeline_bqsr.sh`
4. `GATK4_pipeline_Step_Mutect2_1.1.sh` and `GATK4_pipeline_Step_Mutect2_1.2.sh`
5. `GATK4_pipeline_Step_Mutect2_1.3.sh`, `GATK4_pipeline_Step_Mutect2_2.1.sh`, `GATK4_pipeline_Step_Mutect2_2.2.sh`, `GATK4_pipeline_Step_Mutect2_2.3.sh`, and `GATK4_pipeline_Step_Mutect2_2.4.sh`
6. `GGA_1_CombineVcf.sh`, `GGA_2_Recall.sh`, `GGA_3_Filter.sh`, and `GGA_4_Functator.sh`

GGA recall calls candidate mutations across all tumour samples from the same patient. Consequently, an output locus can have `t_alt_count = 0` in a particular sample.

## Data availability

The de-identified somatic mutation table is distributed separately through OMIX. This repository does not contain sequencing data, clinical information, laboratory sample identifiers, ID-mapping tables, intermediate results, or large output files.

## Reuse

Before running the workflow on a new cohort, configure all paths and reference resources in `config.sh` and provide the required input files for each downstream script. Please cite the associated study when using this code or data.
