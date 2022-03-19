#!/bin/sh
#SBATCH --workdir /work/gr-fe/rueger/G2G-HBV/data/src_gvcf/
#SBATCH --job-name 6_bcftools_gvcf
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --time 10:00:00
#SBATCH --mem 8G
#SBATCH --out /work/gr-fe/rueger/G2G-HBV/data/log/6_bcftools_gvcf.%J.out
#SBATCH --error /work/gr-fe/rueger/G2G-HBV/data/log/6_bcftools_gvcf.%J.err
#SBATCH --mail-type ALL
#SBATCH --mail-user sina.rueeger@epfl.ch
#SBATCH --account gr-fe

BCFTOOLS=/work/gr-fe/rueger/G2G-HBV/bin/bcftools-1.9/bcftools
## number of samples
## bcftools query -l file.bcf | wc -l
#BCFTOOLS=/home/rueger/G2G-HBV/bin/bcftools-1.9/bin/bcftools
DIROUT=/scratch/rueger/gvcf_tmp
mkdir -p $DIROUT

## bcftools query -l file.vcf


INPUT=$DIROUT/hbv_gilead_joint_genotype.fltd-combinedvars.snpEff
OUT=$DIROUT/hbv_gilead_joint_genotype.fltd-combinedvars.snpEff.bcftools

## run BCFTOOLS (script sent by Christian)
$BCFTOOLS filter -i 'QUAL>=30 & INFO/DP>=20 & FMT/DP>=10 & FMT/GQ>=20' -S . -Oz -o $OUT.vcf.gz --threads 30 $INPUT.vcf


echo ended at `date`
exit
