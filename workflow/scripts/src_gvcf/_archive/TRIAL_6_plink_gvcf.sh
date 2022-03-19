#!/bin/sh
#SBATCH --workdir /work/gr-fe/rueger/G2G-HBV/
#SBATCH --job-name 6_plink_gvcf_TRIAL
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 1
#SBATCH --time 01:00:00
#SBATCH --mem 60gb ## 160gb on helvetios, #SBATCH --mem=64gb on deneb2 when reservation
#SBATCH --out /work/gr-fe/rueger/G2G-HBV/data/log/6_plink_gvcf_TRIAL.%J.out
#SBATCH --error /work/gr-fe/rueger/G2G-HBV/data/log/6_plink_gvcf_TRIAL.%J.err
#SBATCH --mail-type ALL
#SBATCH --mail-user sina.rueeger@epfl.ch
## #SBATCH --reservation=gr-fe
 


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
DIROUT=$HOME/data/raw/gilead_20181126/wes_plink_tmp
mkdir -p $DIROUT

# Directories
GVCF=$HOME/data/raw/gilead_20181126/wes_gvcf_tmp
BIN=$HOME/bin
export TMPDIR=/scratch/rueger/$SLURM_JOB_ID

# Tool
GATK=/work/gr-fe/lawless/tool/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar
DBSNP=/work/gr-fe/lawless/ref/b37/dbSnp146.b37.vcf.gz
hg1kv37=/work/gr-fe/lawless/ref/b37/human_g1k_v37.fasta

## START
echo started at `date`


# 7 ----------- turn g.vcf into plink file
## turn vcf into plink (vcf-filter keeps only the FILTER = PASS ones)
## vcf-filter, keeps only high quality snps (PASS): grep "LowQual" hbv_gilead.vcf
$BIN/plink1.9_linux --vcf $DIROUT/hbv_gilead_joint_genotype.ts99pt9.1pcdbsnp.1pcEVS.vep.vcf --vcf-filter --recode --make-bed --threads 20 --out $DIROUT/hbv_gilead


# 8 ------------ cleaning up

#rm -rf $DIROUT/*.list

# 9 ------------- checks

## should be 643 individuals

echo ended at `date`
exit
