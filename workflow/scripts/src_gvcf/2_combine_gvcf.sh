#!/bin/sh
#SBATCH --workdir /work/gr-fe/rueger/G2G-HBV/data/src_gvcf/
#SBATCH --job-name 2_combine_gvcf
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 1
#SBATCH --time 3-0
#SBATCH --mem 60gb ## 160gb on helvetios, #SBATCH --mem=64gb on deneb2 when reservation
#SBATCH --out /work/gr-fe/rueger/G2G-HBV/data/log/2_combine_gvcf.%A.%a.out
#SBATCH --error /work/gr-fe/rueger/G2G-HBV/data/log/2_combine_gvcf.%A.%a.err
#SBATCH --mail-type ALL
#SBATCH --mail-user sina.rueeger@epfl.ch
## #SBATCH --reservation gr-fe
#SBATCH --array 1-30



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
mkdir -p $DIROUT

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
NDAT=30

## START
echo started at `date`


# 2a ----------- create lists for NDAT chunks
echo split files in batches

ls $GVCF/*g.vcf.gz > $DIROUT/vcfs_combine.list
NSAMP=$(wc -l $DIROUT/vcfs_combine.list | cut -f1 -d' ')
SIZE=$(($NSAMP / $NDAT))

## create files with SIZE lines
split -l $SIZE --numeric-suffixes --additional-suffix=.list $DIROUT/vcfs_combine.list $DIROUT/vcfs_combine_split


## clean up all the ones that did not find a place in one of the lists 
## (no clue why this is not working with the above command)

## find all IDs that are not in any other file
## and for simplicity - use R for this
module load gcc/6.4.0 openblas/0.2.20-openmp
module load r/3.5.0
Rscript 2_combine_gvcf_extra.R


# 2b ----------- create these NDAT chunks
echo loop through batches and combine in each batch all files

INPUTALL=($(ls $DIROUT/vcfs_combine_split*.list))

INPUT=${INPUTALL[${SLURM_ARRAY_TASK_ID}]}

INPUTROOT="${INPUT##*/}"
echo starting with files $INPUTROOT

java -jar $GATK \
    -T CombineGVCFs \
    -R $hg1kv37 \
    --variant $INPUT \
    -o $DIROUT/hbv_gilead_$INPUTROOT.g.vcf


echo ended at `date`

