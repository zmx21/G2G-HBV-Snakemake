#!/bin/sh
#SBATCH --workdir /work/gr-fe/rueger/G2G-HBV/data/src_gvcf/
#SBATCH --job-name 3_genotype_gvcf
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 16
#SBATCH --time 48:00:00
#SBATCH --mem 60gb ## 160gb on helvetios, #SBATCH --mem=64gb on deneb2 when reservation
#SBATCH --out /work/gr-fe/rueger/G2G-HBV/data/log/3_genotype_gvcf.%J.out
#SBATCH --error /work/gr-fe/rueger/G2G-HBV/data/log/3_genotype_gvcf.%J.err
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

# create folders
DIROUT=/scratch/rueger/gvcf_tmp

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


## START
echo started at `date`


# 3 ----------- merge NDAT chunks to 1 file

echo joint genotype a cohort

# Load the batches of combined gVCF into genotypeGVCFs and joint genotype a cohort. 
# Better accuracy per variant when the quality per nucleotide for the  whole cohort is compared, rather than individually. 

ls $DIROUT/hbv_gilead_*split*g.vcf > $DIROUT/vcfs_genotype.list

java -Xmx50g -Djava.io.tmpdir=/scratch/rueger/${SLURM_JOBID} -XX:ParallelGCThreads=16 \
    -jar $GATK \
   -T GenotypeGVCFs \
   -R $hg1kv37 \
   -D $DBSNP \
    -stand_call_conf 30 \
    --variant $DIROUT/vcfs_genotype.list \
    -o $DIROUT/hbv_gilead_joint_genotype.vcf -nda --showFullBamList -nt 16 && \
    
    ## KEEP THIS FILE


echo ended at `date`
exit
