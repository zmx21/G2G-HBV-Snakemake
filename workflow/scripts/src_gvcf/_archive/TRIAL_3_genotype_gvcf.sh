#!/bin/sh
#SBATCH --workdir /work/gr-fe/rueger/G2G-HBV/
#SBATCH --job-name=3_genotype_gvcf_TRIAL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-6
#SBATCH --mem=60gb ## 160gb on helvetios, #SBATCH --mem=64gb on deneb2 when reservation
#SBATCH --out=/work/gr-fe/rueger/G2G-HBV/data/log/3_genotype_gvcf_TRIAL.%J.out
#SBATCH --error=/work/gr-fe/rueger/G2G-HBV/data/log/3_genotype_gvcf_TRIAL.%J.err
#SBATCH --mail-type ALL
#SBATCH --mail-user sina.rueeger@epfl.ch
# #SBATCH --reservation=gr-fe
 


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
#DIROUT=/work/gr-fe/rueger/DATA/HBV/gilead_20181126/wes_plink_tmp
DIROUT=$HOME/data/raw/gilead_20181126/wes_plink
mkdir -p $DIROUT

# Directories
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



# 4 ----------- merge NDAT chunks to 1 file

echo joint genotype a cohort

# Load the batches of combined gVCF into genotypeGVCFs and joint genotype a cohort. 
# Better accuracy per variant when the quality per nucleotide for the  whole cohort is compared, rather than individually. 

java -Xmx50g -Djava.io.tmpdir=/scratch/rueger/${SLURM_JOBID} -XX:ParallelGCThreads=1 \
    -jar $GATK \
    -T GenotypeGVCFs \
    -R $hg1kv37 \
    -D $DBSNP \
    -stand_call_conf 30 \
    -V $DIROUT/hbv_gilead_vcfs_combine_split00.list.g.vcf \
    -o $DIROUT/hbv_gilead_vcfs_combine_genotyped.vcf -nda --showFullBamList -nt 1 && \




echo ended at `date`
