#!/bin/sh
#SBATCH --workdir /work/gr-fe/rueger/G2G-HBV/data/src_gvcf/
#SBATCH --job-name 1_tabix_gvcf
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 1
#SBATCH --time 02:00:00
#SBATCH --mem 60gb ## 160gb on helvetios, #SBATCH --mem=64gb on deneb2 when reservation
#SBATCH --out /work/gr-fe/rueger/G2G-HBV/data/log/1_tabix_gvcf.%J.out
#SBATCH --error /work/gr-fe/rueger/G2G-HBV/data/log/1_tabix_gvcf.%J.err
#SBATCH --mail-type ALL
#SBATCH --mail-user sina.rueeger@epfl.ch
## #SBATCH --reservation gr-fe
 


# Load all modules required
module load gcc/6.4.0
module load gcc/7.3.0
module load intel/18.0.2
module load samtools/1.8
module load picard/2.18.14

## if grfepc4
#HOME=/home/rueger/G2G-HBV
HOME=/work/gr-fe/rueger/G2G-HBV

# Directories
GVCF=$HOME/data/raw/gilead_20181126/wes_gvcf
BIN=$HOME/bin
export TMPDIR=/scratch/rueger/$SLURM_JOB_ID

# Tool
GATK=/work/gr-fe/lawless/tool/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar
DBSNP=/work/gr-fe/lawless/ref/b37/dbSnp146.b37.vcf.gz
hg1kv37=/work/gr-fe/lawless/ref/b37/human_g1k_v37.fasta

## parameters
## how many datasets
NDAT=10

## START
echo started at `date`

## 1 -------- export stuff for perl and tabix
echo generate tabix files

export PERL5LIB=$BIN/vcftools_0.1.13/perl/
export PATH=${PATH}:$BIN/tabix-0.2.6/

# create tbi files
for file in $GVCF/*g.vcf.gz
do
  $BIN/tabix-0.2.6/tabix -p vcf -f $file
done


echo ended at `date`
exit
