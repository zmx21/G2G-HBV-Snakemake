#!/bin/sh
#SBATCH --workdir /work/gr-fe/rueger/G2G-HBV/data/
#SBATCH --job-name 5_annotation_gvcf_TRIAL
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 16
#SBATCH --time 03:00:00
#SBATCH --mem 60gb ## 160gb on helvetios, #SBATCH --mem=64gb on deneb2 when reservation
#SBATCH --out /work/gr-fe/rueger/G2G-HBV/data/log/5_annotation_gvcf_TRIAL.%J.out
#SBATCH --error /work/gr-fe/rueger/G2G-HBV/data/log/5_annotation_gvcf_TRIAL.%J.err
#SBATCH --mail-type ALL
#SBATCH --mail-user sina.rueeger@epfl.ch
## #SBATCH --reservation=gr-fe


# Load all modules required
module load intel/18.0.2
module load htslib/1.8
module load perl/5.24.1

export LANGUAGE=en_US.UTF-8
export LANG=en_US.UTF-8
export LC_ALL=en_US.UTF-8

cpanm DBI
cpanm CGI

## if grfepc4
#HOME=/home/rueger/G2G-HBV
HOME=/work/gr-fe/rueger/G2G-HBV

# create folders
DIROUT=$HOME/data/raw/gilead_20181126/wes_plink_tmp

# Directories
GVCF=$HOME/data/raw/gilead_20181126/wes_gvcf_tmp
BIN=$HOME/bin
export TMPDIR=/scratch/rueger/$SLURM_JOB_ID


# Tool
GATK=/work/gr-fe/lawless/tool/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar
DBSNP=/work/gr-fe/lawless/ref/b37/dbSnp146.b37.vcf.gz
hg1kv37=/work/gr-fe/lawless/ref/b37/human_g1k_v37.fasta
VEP=/work/gr-fe/lawless/ref/variant_effect_predictor

## START
echo started at `date`

# 6 ----------- Annotation

# annotate with Variant Effect Predictor
perl $VEP/variant_effect_predictor.pl \
--offline --vcf --everything \
--dir_cache $VEP/vep_cache \
--dir_plugins $VEP/vep_cache/Plugins \
--plugin Condel,$VEP/vep_cache/Plugins/config/Condel/config/ \
--plugin ExAC,/home/gr-fe/lawless/ref/ExAC/ExAC.r0.3.sites.vep.vcf.gz \
--plugin SpliceConsensus \
--fasta $VEP/cache_b37/homo_sapiens/92_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz
-i $DIROUT/hbv_gilead_joint_genotype.vcf \
-o $DIROUT/hbv_gilead_joint_genotype.ts99pt9.1pcdbsnp.1pcEVS.vep.vcf \
--fork 1 && \


echo ended at `date`
exit
