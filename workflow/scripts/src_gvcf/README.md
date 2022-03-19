See also [issue 107](https://github.com/sinarueeger/G2G-HBV/issues/107).

This directory contains all the scripts needed to create a PLINK file from all the single gvcf files in data/raw/gilead_20181126/wes_gvcf/*vcf.

Pipeline to create the plink file
.
|-- 1_tabix_gvcf.sh         ## this file never really ran, bc tabix, combine and genotype were together once
|-- 2_combine_gvcf.sh       ## array job (1-30)
|-- 3_genotype_gvcf.sh      ## genotyping 
|-- 4_filter_gvcf.sh        ## hard filtering
|-- 5_annotation_gvcf.sh    ## adding gene + function
|-- 6_misingness_filter.sh  ## missingness bcftools
|-- 7_plink_gvcf.sh         ## transform into a plink format
`-- _archive/

_archive contains (1) and (2).

(1)
TRIAL pipeline: There are three scripts that only use gvcf files
.
|-- TRIAL_1_tabix_combine_genotype_gvcf.sh  ## old file that used to tabix, combine and genotype
|-- TRIAL_2_combine_gvcf.sh     ## combine using an array job (new, replaces the TRIAL_1* one)
|-- TRIAL_2_filter_gvcf.sh      ## filter won't work on 23 ppl only (only works if hard filter is used from file combine_gvcf_dylan.sh)
|-- TRIAL_3_annotation_gvcf.sh  ## run annotation (perl does not work right now)
`-- TRIAL_4_plink_gvcf.sh       ## convert into plink

(2)
File that dylan sent me with instructions
.
`-- combine_gvcf_dylan.sh       ## combine, genotype, hard filter, annotation
`-- vep_gcvf_dylan.sh           ## soft filter, annotation (perl won't run though)
