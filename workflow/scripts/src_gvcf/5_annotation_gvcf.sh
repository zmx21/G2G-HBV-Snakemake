#!/bin/sh
#SBATCH --workdir /work/gr-fe/rueger/G2G-HBV/src_gvcf/
#SBATCH --job-name 5_annotation_gvcf
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --time 24:00:00
#SBATCH --mem 8G
#SBATCH --out /work/gr-fe/rueger/G2G-HBV/data/log/5_annotation_gvcf.%J.out
#SBATCH --error /work/gr-fe/rueger/G2G-HBV/data/log/5_annotation_gvcf.%J.err
#SBATCH --mail-type ALL
#SBATCH --mail-user sina.rueeger@epfl.ch
#SBATCH --account gr-fe

## took 7 hours last time

# 5 ----------- Annotation

DIROUT=/scratch/rueger/gvcf_tmp
mkdir -p $DIROUT

FILE=$DIROUT/hbv_gilead_joint_genotype.fltd-combinedvars  #$1

SNPEFF=/work/gr-fe/rueger/G2G-HBV/bin/snpEff/snpEff.jar

# Annotate
java -Xmx8g -jar $SNPEFF \
        -c /work/gr-fe/rueger/G2G-HBV/bin/snpEff/snpEff.config GRCh37.75 -v -o gatk $FILE.vcf > $FILE.snpEff.vcf
        
## Store the annotation in a separate file
sed '/^##/ d' $FILE.snpEff.vcf | awk 'BEGIN {OFS ="," ; FS = "\t"};{print $1, $2, $3, $8}' > $DIROUT/snpEff_annotation.csv

echo ended at `date`
exit
