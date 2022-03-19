#!/bin/sh
#SBATCH --workdir /work/gr-fe/rueger/G2G-HBV/data/src_gvcf/
#SBATCH --job-name 7_plink_gvcf
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 1
#SBATCH --time 07:00:00
#SBATCH --mem 60gb ## 160gb on helvetios, #SBATCH --mem=64gb on deneb2 when reservation
#SBATCH --out /work/gr-fe/rueger/G2G-HBV/data/log/7_plink_gvcf.%J.out
#SBATCH --error /work/gr-fe/rueger/G2G-HBV/data/log/7_plink_gvcf.%J.err
#SBATCH --mail-type ALL
#SBATCH --mail-user sina.rueeger@epfl.ch
#SBATCH --account gr-fe
 


# Load all modules required
module load gcc/6.4.0
module load gcc/7.3.0
module load intel/18.0.2
module load samtools/1.8
module load picard/2.18.14

## if grfepc4
#HOME=/home/rueger/G2G-HBV
HOME=/work/gr-fe/rueger/G2G-HBV

# create folders
DIROUT=/scratch/rueger/gvcf_tmp
DIRFINAL=$HOME/data/raw/gilead_20181126/wes_plink

# Directories
GVCF=$HOME/data/raw/gilead_20181126/wes_gvcf
BIN=$HOME/bin
export TMPDIR=/scratch/rueger/$SLURM_JOB_ID

## START
echo started at `date`

# 5 ----------- turn g.vcf into plink file
## turn vcf into plink (vcf-filter keeps only the FILTER = PASS ones)
##  --vcf-filter
## vcf-filter, keeps only high quality snps (PASS): grep "LowQual" hbv_gilead.vcf

$BIN/plink1.9_linux --vcf $DIROUT/hbv_gilead_joint_genotype.fltd-combinedvars.snpEff.bcftools.vcf.gz --set-missing-var-ids @:#:\$1:\$2 --keep-allele-order --recode --make-bed --threads 20 --out $DIRFINAL/hbv_gilead

#$BIN/plink1.9_linux --vcf $DIROUT/hbv_gilead_joint_genotype.fltd-combinedvars.snpEff.vcf --set-missing-var-ids @:#:\$1:\$2 --keep-allele-order --recode --make-bed --threads 20 --out $DIRFINAL/hbv_gilead

# 6 ------------ cleaning up

#rm -rf $DIROUT/*.list

# 7 ------------- checks

## should be 643 individuals


echo ended at `date`
exit
