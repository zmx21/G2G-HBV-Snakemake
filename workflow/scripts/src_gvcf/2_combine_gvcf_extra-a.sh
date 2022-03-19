#!/bin/sh
#SBATCH --workdir /work/gr-fe/rueger/G2G-HBV/data/src_gvcf/
#SBATCH --job-name 2b_combine_gvcf_extra
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 1
#SBATCH --time 48:00:00
#SBATCH --mem 60gb ## 160gb on helvetios, #SBATCH --mem=64gb on deneb2 when reservation
#SBATCH --out /work/gr-fe/rueger/G2G-HBV/data/log/2b_combine_gvcf_extra.%A.%a.out
#SBATCH --error /work/gr-fe/rueger/G2G-HBV/data/log/2b_combine_gvcf_extra.%A.%a.err
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
HOME=/work/gr-fe/rueger/G2G-HBV


# Directories
GVCF=$HOME/data/raw/gilead_20181126/wes_gvcf
BIN=$HOME/bin
export TMPDIR=/scratch/rueger/$SLURM_JOB_ID

# Tool
GATK=/work/gr-fe/lawless/tool/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar
hg1kv37=/work/gr-fe/lawless/ref/b37/human_g1k_v37.fasta

DIROUT=/scratch/rueger/gvcf_tmp

java -jar $GATK \
    -T CombineGVCFs \
    -R $hg1kv37 \
    --variant $DIROUT/vcfs_combine_vcfs_split31.list \
    -o $DIROUT/hbv_gilead_combine_vcfs_split31.list.g.vcf
