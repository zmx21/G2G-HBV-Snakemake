#!/bin/sh
#SBATCH --workdir /work/gr-fe/rueger/G2G-HBV/
#SBATCH --job-name=1_tabix_combine_genotype_gvcf_TRIAL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-6
#SBATCH --mem=60gb ## 160gb on helvetios, #SBATCH --mem=64gb on deneb2 when reservation
#SBATCH --out=/work/gr-fe/rueger/G2G-HBV/data/log/1_tabix_combine_genotype_gvcf_TRIAL.%J.out
#SBATCH --error=/work/gr-fe/rueger/G2G-HBV/data/log/1_tabix_combine_genotype_gvcf_TRIAL.%J.err
#SBATCH --mail-type ALL
#SBATCH --mail-user sina.rueeger@epfl.ch
#SBATCH --reservation=gr-fe
 


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

## parameters
## how many datasets
NDAT=2

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


# 2 ----------- create lists for NDAT chunks
echo split files in batches

ls $GVCF/*g.vcf.gz > $DIROUT/vcfs_combine.list
NSAMP=$(wc -l $DIROUT/vcfs_combine.list | cut -f1 -d' ')
SIZE=$(($NSAMP / $NDAT))

## create files with SIZE lines
split -l $SIZE --numeric-suffixes --additional-suffix=.list $DIROUT/vcfs_combine.list $DIROUT/vcfs_combine_split


# 3 ----------- create these NDAT chunks
echo loop through batches and combine in each batch all files

for INPUT in $(ls $DIROUT/vcfs_combine_split*.list)
do

INPUTROOT="${INPUT##*/}"
echo starting with files $INPUTROOT

## might not need that first line# java -jar $GATK \

java -jar $GATK \
    -T CombineGVCFs \
    -R $hg1kv37 \
    --variant $INPUT \
    -o $DIROUT/hbv_gilead_$INPUTROOT.g.vcf

done

# 4 ----------- merge NDAT chunks to 1 file

echo joint genotype a cohort

# Load the batches of combined gVCF into genotypeGVCFs and joint genotype a cohort. 
# Better accuracy per variant when the quality per nucleotide for the  whole cohort is compared, rather than individually. 

ls $DIROUT/hbv_gilead_*split*g.vcf > $DIROUT/vcfs_genotype.list

java -Xmx50g -Djava.io.tmpdir=/scratch/rueger/${SLURM_JOBID} -XX:ParallelGCThreads=1 \
    -jar $GATK \
   -T GenotypeGVCFs \
   -R $hg1kv37 \
   -D $DBSNP \
    -stand_call_conf 30 \
    --variant $DIROUT/vcfs_genotype.list \
    -o $DIROUT/hbv_gilead.g.vcf -nda --showFullBamList -nt 1 && \

# from dylan
# java -Xmx150g -Djava.io.tmpdir=/scratch/lawless/${SLURM_JOBID} -XX:ParallelGCThreads=24 \
# -jar $GATK \
# -T GenotypeGVCFs \
# -R $hg1kv37 \
# -D $DBSNP \
# -stand_call_conf 30 \
# -V $GENO/cohort.batch.1.combined.g.vcf \
# -V $GENO/cohort.batch.2.combined.g.vcf \
# -V $GENO/cohort.batch.3.combined.g.vcf \
# -o $GENO/cohort.joint.genotype.vcf -nda --showFullBamList -nt  && \

    


echo ended at `date`
