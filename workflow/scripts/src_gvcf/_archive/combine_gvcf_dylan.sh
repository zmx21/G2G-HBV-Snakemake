#!/bin/sh
#SBATCH --mail-type END
#SBATCH --mail-user email@epfl.ch
#SBATCH -J genotypeGVCF
#SBATCH -o ./log/genotypeGVCF.%J.out
#SBATCH -e ./log/genotypeGVCF.%J.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task 1
#SBATCH --time=14:00:00
#SBATCH --mem=160gb
#SBATCH --reservation=gr-fe

# helvetios has 36 cores per node, 1 thread per core and ~190GB ram.
# some GATK commands can use multithreading (documentation options list either -nt or -nct). 
# However, there is 1 thread per core on scitas. 
# So I think the options for genotypeGVCFs on 36 threads (-nt 36) would be #SBATCH --cpus-per-task 36

# set -e

echo started at `date`

# swap scratch username, not sure if you can write to my scratch dir
export TMPDIR=/scratch/lawless/$SLURM_JOB_ID

# Load modules required
module load gcc/6.4.0
module load gcc/7.3.0
module load intel/18.0.2

# Directories, list you own
GVCF=/work/gr-fe/.....
GENO=/work/gr-fe/.....

# References, copy to your own or read from mine
hg1kv37=/work/gr-fe/lawless/ref/b37/human_g1k_v37.fasta
GATK=/work/gr-fe/lawless/tool/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar
DBSNP=/work/gr-fe/lawless/ref/b37/dbSnp146.b37.vcf.gz

# # Combine gGVCFs into batches, so that genotypeGVCFs can load a smaller number of files.

# java -jar $GATK \
# -T CombineGVCFs \
# -R $hg1kv37 \
# -V $GVCF/001.HC.g.vcf \
# -V $GVCF/002.HC.g.vcf \
# -V $GVCF/004.HC.g.vcf \
# -V $GVCF/005.HC.g.vcf \
# -V $GVCF/007.HC.g.vcf \
# -V $GVCF/008.HC.g.vcf \
# -V $GVCF/099.HC.g.vcf \
# -o $GENO/cohort.batch.1.combined.g.vcf && \




# Load the batches of combined gVCF into genotypeGVCFs and joint genotype a cohort. 
# Better accuracy per variant when the quality per nucleotide for the  whole cohort is compared, rather than individually. 

# java -Xmx150g -Djava.io.tmpdir=/scratch/lawless/${SLURM_JOBID} -XX:ParallelGCThreads=24 \
# -jar $GATK \
# -T GenotypeGVCFs \
# -R $hg1kv37 \
# -D $DBSNP \
# -stand_call_conf 30 \
# -V $GENO/cohort.batch.1.combined.g.vcf \
# -V $GENO/cohort.batch.2.combined.g.vcf \
# -V $GENO/cohort.batch.3.combined.g.vcf \
# -o $GENO/cohort.joint.genotype.vcf -nda --showFullBamList -nt 24 && \

# # Hard filter for SNPs and Indels
# # Each done seperately, then join the two datasets (select + filter SNPs, select + filter Indels)

# java -jar $GATK \
# -T SelectVariants \
# -R $hg1kv37 \
# -selectType SNP \
# --variant $GENO/cohort.joint.genotype.vcf \
# -o $GENO/cohort.joint.genotype.raw-snps.vcf && \
 
# java -jar $GATK \
# -T SelectVariants \
# -R $hg1kv37 \
# --variant $GENO/cohort.joint.genotype.vcf \
# -selectType INDEL -selectType MNP \
# -o $GENO/cohort.joint.genotype.raw-indels.vcf && \
 
# java -jar $GATK \
# -T VariantFiltration \
# -R $hg1kv37 \
# -V $GENO/cohort.joint.genotype.raw-snps.vcf \
# --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MappingQualityRankSum < -12.5 || ReadPosRankSum < -8.0" \
# --filterName "snp_hard_filter" \
# -o $GENO/cohort.joint.genotype.raw-snps.filtered.snvs.vcf && \
 
# java -jar $GATK \
# -T VariantFiltration \
# -R $hg1kv37 \
# -V $GENO/cohort.joint.genotype.raw-indels.vcf \
# --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
# --filterName "indel_hard_filter" \
# -o $GENO/cohort.joint.genotype.raw-indels.filtered.indels.vcf && \
 
# # Now join the two datasets back together
# java -jar $GATK \
# -T CombineVariants \
# -R $hg1kv37 \
# --variant $GENO/cohort.joint.genotype.raw-snps.filtered.snvs.vcf \
# --variant $GENO/cohort.joint.genotype.raw-indels.filtered.indels.vcf \
# -o $GENO/cohort.joint.genotype.fltd-combinedvars.vcf \
# --genotypemergeoption UNSORTED && \

# # Annotate with Variant Effect Predictor

# perl /work/gr-fe/lawless/ref/variant_effect_predictor/variant_effect_predictor.pl \
# --offline --vcf --everything \
# --dir_cache /home/variant_effect_predictor/vep_cache \
# --dir_plugins /home/variant_effect_predictor/vep_cache/Plugins \
# --plugin Condel,/home/variant_effect_predictor/vep_cache/Plugins/config/Condel/config/ \
# --plugin ExAC,/home/ref/ExAC/ExAC.r0.3.sites.vep.vcf.gz \
# --plugin SpliceConsensus \
# --fasta /home/variant_effect_predictor/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
# -i $GENO/cohort.joint.genotype.fltd-combinedvars.vcf \
# -o $GENO/cohort.joint.genotype.fltd-combinedvars.vep.vcf \
# --fork 16 && \

exit
