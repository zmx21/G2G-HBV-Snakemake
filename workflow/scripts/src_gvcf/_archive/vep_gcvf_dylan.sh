# Directories
GENO=/work/gr-fe/....

# Tools
VEP=/work/gr-fe/lawless/ref/variant_effect_predictor
vcfhacks=/work/gr-fe/lawless/tool/vcfhacks

# References
hg1kv37=/work/gr-fe/lawless/ref/b37/human_g1k_v37.fasta
GATK=/work/gr-fe/lawless/tool/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar
K1GINDEL=/work/gr-fe/lawless/ref/b37/1000G_phase1.indels.b37.vcf
MILLS1000G=/work/gr-fe/lawless/ref/b37/Mills_and_1000G_gold_standard.indels.b37.sites.vcf
DBSNP=/work/gr-fe/lawless/ref/b37/dbSnp146.b37.vcf.gz
SSV5=/work/gr-fe/lawless/ref/SureSelectAllExonV5/S04380110_Regions_b37.bed
HAPMAPSITES=/work/gr-fe/lawless/ref/b37/hapmap_3.3.b37.sites.vcf
K1GOMNISITES=/work/gr-fe/lawless/ref/b37/1000G_omni2.5.b37.sites.vcf
K1GSNP=/work/gr-fe/lawless/ref/b37/1000G_phase1.snps.high_confidence.b37.vcf

# Options to filter for quality scores. 
# There are two options:
# Option 1 - Variant quality score recalibration. In general, said to be better and based on assessing the whole cohort and recalibrating low quality calls. https://gatkforums.broadinstitute.org/gatk/discussion/39/variant-quality-score-recalibration-vqsr
# Option 2 - Hard filter, simpler algorith, search for SNPs and Indels seperately then combine both together. See this page and the first link "why you would use hard filters" https://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set  
# I cannot remember which is better for single cases in a cohort, versus assossiation type senario with higher % carying a variant. 

# # The next 4 commands make up the VQSR method
# # Pick SNPs and make tranche file
# java -Xmx24g -Djava.io.tmpdir=/scratch/lawless/${SLURM_JOBID} \
# -jar $GATK \
# -T VariantRecalibrator \
# -R $hg1kv37 \
# -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $HAPMAPSITES \
# -resource:omni,known=false,training=true,truth=true,prior=12.0 $K1GOMNISITES \
# -resource:1000G,known=false,training=true,truth=false,prior=10.0 $K1GSNP \
# -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $DBSNP \
# -an QD -an MQRankSum -an ReadPosRankSum -an FS \
# -mode SNP \
# -input $GENO/pri.leeds.hiv.combined.genotype.vcf \
# -recalFile $GENO/pri.leeds.hiv.combined.genotype.snvrecal \
# -tranchesFile $GENO/pri.leeds.hiv.combined.genotype.snvrecal.tranches \
# -rscriptFile $GENO/pri.leeds.hiv.combined.genotype.snvrecal.plots.R -nt 16 && \
#
# # create a tranche file of Indels, again on the ORIGIN input vcf, not the output of the last step
# java -Xmx24g -Djava.io.tmpdir=/scratch/lawless/${SLURM_JOBID} \
# -jar $GATK \
# -T VariantRecalibrator \
# -R $hg1kv37 \
# -resource:mills,known=false,training=true,truth=true,prior=12.0 $MILLS1000G \
# -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $DBSNP \
# -an FS -an ReadPosRankSum -an MQRankSum -mode INDEL \
# -input $GENO/pri.leeds.hiv.combined.genotype.vcf \
# -recalFile $GENO/pri.leeds.hiv.combined.genotype.indelrecal \
# -tranchesFile $GENO/pri.leeds.hiv.combined.genotype.indelrecal.tranches \
# -rscriptFile $GENO/pri.leeds.hiv.combined.genotype.indelrecal.plots.R -nt 16 && \
#
# # Apply that recalibration of SNPs
# java -Xmx24g -Djava.io.tmpdir=/scratch/lawless/${SLURM_JOBID} \
# -jar $GATK \
# -T ApplyRecalibration \
# -R $hg1kv37 \
# -input $GENO/pri.leeds.hiv.combined.genotype.vcf \
# -recalFile $GENO/pri.leeds.hiv.combined.genotype.snvrecal \
# -tranchesFile $GENO/pri.leeds.hiv.combined.genotype.snvrecal.tranches \
# -o $GENO/pri.leeds.hiv.combined.genotype.snvts99pt9.vcf \
# -mode SNP --ts_filter_level 99.9 -nt 16 && \
# 
# # Apply that recalibration of Indels to the previous output
# java -Xmx24g -Djava.io.tmpdir=/scratch/lawless/${SLURM_JOBID} \
# -jar $GATK \
# -T ApplyRecalibration \
# -R $hg1kv37 \
# -input $GENO/pri.leeds.hiv.combined.genotype.snvts99pt9.vcf \
# -recalFile $GENO/pri.leeds.hiv.combined.genotype.indelrecal \
# -tranchesFile $GENO/pri.leeds.hiv.combined.genotype.indelrecal.tranches  \
# -o $GENO/pri.leeds.hiv.combined.genotype.ts99pt9.vcf \
# -mode INDEL --ts_filter_level 99.9 -nt 16 && \

# # Skip to keep common variants
# # filter variants in dbSNP >/= 1% and not listed as pothogenic by ClinVar
# perl $vcfhacks/annotateSnps.pl \
# -d $DBSNP /work/gr-fe/lawless/ref/b37/clinvar_20160531.vcf.gz -f 1 -pathogenic \
# -i $GENO/cohort.combined.genotype.ts99pt9.vcf \
# -o $GENO/cohort.combined.genotype.ts99pt9.1pcdbsnp.vcf -t 16 && \
# 
# # Skip to keep common variants
# # filter variants in EVS (exome variants server) greater >/= 1%
# perl $vcfhacks/filterOnEvsMaf.pl \
# -d /work/gr-fe/lawless/ref/evs/ -f 1 --progress \
# -i $GENO/cohort.combined.genotype.ts99pt9.1pcdbsnp.vcf \
# -o $GENO/cohort.combined.genotype.ts99pt9.1pcdbsnp.1pcEVS.vcf -t 16 && \

# annotate with Variant Effect Predictor
perl $VEP/variant_effect_predictor.pl \
--offline --vcf --everything \
--dir_cache $VEP/vep_cache \
--dir_plugins $VEP/vep_cache/Plugins \
--plugin Condel,$VEP/vep_cache/Plugins/config/Condel/config/ \
--plugin ExAC,/home/gr-fe/lawless/ref/ExAC/ExAC.r0.3.sites.vep.vcf.gz \
--plugin SpliceConsensus \
--fasta $VEP/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
-i $GENO/cohort.combined.genotype.ts99pt9.1pcdbsnp.1pcEVS.vcf \
-o $GENO/cohort.combined.genotype.ts99pt9.1pcdbsnp.1pcEVS.vep.vcf \
--fork 16 && \
