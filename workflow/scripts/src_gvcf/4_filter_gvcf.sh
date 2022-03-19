#!/bin/sh
#SBATCH --workdir /work/gr-fe/rueger/G2G-HBV/data/src_gvcf/
#SBATCH --job-name 4_filter_gvcf
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 1
#SBATCH --time 12:00:00
#SBATCH --mem 60gb ## 160gb on helvetios, #SBATCH --mem=64gb on deneb2 when reservation
#SBATCH --out /work/gr-fe/rueger/G2G-HBV/data/log/4_filter_gvcf.%J.out
#SBATCH --error /work/gr-fe/rueger/G2G-HBV/data/log/4_filter_gvcf.%J.err
#SBATCH --mail-type ALL
#SBATCH --mail-user sina.rueeger@epfl.ch
## #SBATCH --reservation gr-fe
 

# Load modules required
module load gcc/6.4.0
module load gcc/7.3.0
module load intel/18.0.2

# Directories
DIROUT=/scratch/rueger/gvcf_tmp
export TMPDIR=/scratch/rueger/$SLURM_JOB_ID

## if grfepc4
#HOME=/home/rueger/G2G-HBV
HOME=/work/gr-fe/rueger/G2G-HBV

# Tools
VEP=/work/gr-fe/lawless/ref/variant_effect_predictor
vcfhacks=/work/gr-fe/lawless/tool/vcfhacks

# References
hg1kv37=/work/gr-fe/lawless/ref/b37/human_g1k_v37.fasta
GATK=/work/gr-fe/lawless/tool/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar
K1GINDEL=/work/gr-fe/lawless/ref/b37/1000G_phase1.indels.b37.vcf
MILLS1000G=/work/gr-fe/lawless/ref/b37/Mills_and_1000G_gold_standard.indels.b37.sites.vcf
DBSNP=/work/gr-fe/lawless/ref/b37/dbSnp146.b37.vcf.gz
SSV5=/work/gr-fe/lawless/ref/SureSelectAllExonV5/S04380110_Regions_b37.bed
HAPMAPSITES=/work/gr-fe/lawless/ref/b37/hapmap_3.3.b37.sites.vcf
K1GOMNISITES=/work/gr-fe/lawless/ref/b37/1000G_omni2.5.b37.sites.vcf
K1GSNP=/work/gr-fe/lawless/ref/b37/1000G_phase1.snps.high_confidence.b37.vcf


## START
echo started at `date`


# 4 ----------- SOFT filter for SNPs and Indels
# # Each done seperately, then join the two datasets (select + filter SNPs, select + filter Indels)


# Options to filter for quality scores. 
# There are two options:
# Option 1 - Variant quality score recalibration. In general, said to be better and based on assessing the whole cohort and recalibrating low quality calls. https://gatkforums.broadinstitute.org/gatk/discussion/39/variant-quality-score-recalibration-vqsr
# Option 2 - Hard filter, simpler algorith, search for SNPs and Indels seperately then combine both together. See this page and the first link "why you would use hard filters" https://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set  
# I cannot remember which is better for single cases in a cohort, versus assossiation type senario with higher % carying a variant. 

#filter=soft ##option 1
filter=hard ## option 2

# # Hard filter for SNPs and Indels
# # Each done seperately, then join the two datasets (select + filter SNPs, select + filter Indels)

if [ $filter = "hard" ]
then
echo applying hard filter


  java -jar $GATK \
  -T SelectVariants \
  -R $hg1kv37 \
  -selectType SNP \
  --variant $DIROUT/hbv_gilead_joint_genotype.vcf \
  -o $DIROUT/hbv_gilead_joint_genotype_raw-snps.vcf

 
  java -jar $GATK \
  -T SelectVariants \
  -R $hg1kv37 \
  --variant $DIROUT/hbv_gilead_joint_genotype.vcf \
  -selectType INDEL -selectType MNP \
  -o $DIROUT/hbv_gilead_joint_genotype_raw-indels.vcf
 
  java -jar $GATK \
  -T VariantFiltration \
  -R $hg1kv37 \
  -V $DIROUT/hbv_gilead_joint_genotype_raw-snps.vcf \
  --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MappingQualityRankSum < -12.5 || ReadPosRankSum < -8.0" \
  --filterName "snp_hard_filter" \
  -o $DIROUT/hbv_gilead_joint_genotype_raw-snps.filtered.snvs.vcf 

  java -jar $GATK \
  -T VariantFiltration \
  -R $hg1kv37 \
  -V $DIROUT/hbv_gilead_joint_genotype_raw-indels.vcf \
  --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
  --filterName "indel_hard_filter" \
  -o $DIROUT/hbv_gilead_joint_genotype_raw-indels.filtered.indels.vcf
 
 
  # Now join the two datasets back together
  java -jar $GATK \
  -T CombineVariants \
  -R $hg1kv37 \
  --variant $DIROUT/hbv_gilead_joint_genotype_raw-snps.filtered.snvs.vcf \
  --variant $DIROUT/hbv_gilead_joint_genotype_raw-indels.filtered.indels.vcf \
  -o $DIROUT/hbv_gilead_joint_genotype.fltd-combinedvars.vcf \
  --genotypemergeoption UNSORTED

fi

## OPTION 1 - variant quality recalibration
## ----------------------------------------

if [ $filter = "soft" ]
then
echo applying soft filter


# Pick SNPs and make tranche file
java -Xmx24g -Djava.io.tmpdir=/scratch/rueger/${SLURM_JOBID} \
  -jar $GATK \
  -T VariantRecalibrator \
  -R $hg1kv37 \
  -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $HAPMAPSITES \
  -resource:omni,known=false,training=true,truth=false,prior=12.0 $K1GOMNISITES \  ## truth-false
  -resource:1000G,known=false,training=true,truth=false,prior=10.0 $K1GSNP \
  -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $DBSNP \
  -an QD -an MQRankSum -an ReadPosRankSum -an FS \
  -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
  -mode SNP \
  -input $DIROUT/hbv_gilead_joint_genotype.vcf \
  -recalFile $DIROUT/hbv_gilead_joint_genotype.snvrecal \
  -tranchesFile $DIROUT/hbv_gilead_joint_genotype.snvrecal.tranches \
  -rscriptFile $DIROUT/hbv_gilead_joint_genotype.snvrecal.plots.R -nt 1 && \


# create a tranche file of Indels, again on the ORIGIN input vcf, not the output of the last step
java -Xmx24g -Djava.io.tmpdir=/scratch/rueger/${SLURM_JOBID} \
  -jar $GATK \
  -T VariantRecalibrator \
  -R $hg1kv37 \
  -resource:mills,known=false,training=true,truth=true,prior=12.0 $MILLS1000G \
  -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $DBSNP \
  -an FS -an ReadPosRankSum -an MQRankSum -mode INDEL \
  -input $DIROUT/hbv_gilead_joint_genotype.vcf \
  -recalFile $DIROUT/hbv_gilead_joint_genotype.indelrecal \
  -tranchesFile $DIROUT/hbv_gilead_joint_genotype.indelrecal.tranches \
  -rscriptFile $DIROUT/hbv_gilead_joint_genotype.indelrecal.plots.R -nt 1 && \

# Apply that recalibration of SNPs
java -Xmx24g -Djava.io.tmpdir=/scratch/rueger/${SLURM_JOBID} \
  -jar $GATK \
  -T ApplyRecalibration \
  -R $hg1kv37 \
  -input $DIROUT/hbv_gilead_joint_genotype.vcf \
  -recalFile $DIROUT/hbv_gilead_joint_genotype.snvrecal \
  -tranchesFile $DIROUT/hbv_gilead_joint_genotype.snvrecal.tranches \
  -o $DIROUT/hbv_gilead_joint_genotype.snvts99pt9.vcf \
  -mode SNP --ts_filter_level 99.9 -nt 1 && \
 
# Apply that recalibration of Indels to the previous output
java -Xmx24g -Djava.io.tmpdir=/scratch/rueger/${SLURM_JOBID} \
  -jar $GATK \
  -T ApplyRecalibration \
  -R $hg1kv37 \
  -input $DIROUT/hbv_gilead_joint_genotype.snvts99pt9.vcf \
  -recalFile $DIROUT/hbv_gilead_joint_genotype.indelrecal \
  -tranchesFile $DIROUT/hbv_gilead_joint_genotype.indelrecal.tranches  \
  -o $DIROUT/hbv_gilead_joint_genotype.ts99pt9.vcf \
  -mode INDEL --ts_filter_level 99.9 -nt 1 && \

# Skip to keep common variants
# filter variants in dbSNP >/= 1% and not listed as pothogenic by ClinVar
perl $vcfhacks/annotateSnps.pl \
  -d $DBSNP /work/gr-fe/lawless/ref/b37/clinvar_20160531.vcf.gz -f 1 -pathogenic \
  -i $DIROUT/hbv_gilead_joint_genotype.ts99pt9.vcf \
  -o $DIROUT/hbv_gilead_joint_genotype.1pcdbsnp.vcf -t 16 && \
 
# Skip to keep common variants
# filter variants in EVS (exome variants server) greater >/= 1%
perl $vcfhacks/filterOnEvsMaf.pl \
  -d /work/gr-fe/lawless/ref/evs/ -f 1 --progress \
  -i $DIROUT/hbv_gilead_joint_genotype.ts99pt9.1pcdbsnp.vcf \
  -o $DIROUT/hbv_gilead_joint_genotype.ts99pt9.1pcdbsnp.1pcEVS.vcf -t 16 && \

fi

echo ended at `date`
exit
