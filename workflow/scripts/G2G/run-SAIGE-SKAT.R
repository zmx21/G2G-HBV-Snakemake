library(SAIGE)
CreateGRMGenoFileBURDEN <- function(results_path,n_cores_input,max_maf = 0.05){
  cur_dir <- here::here()
  load(glue::glue("{results_path}/prep-data.rda"))
  #Write out PLINK files without Chr 6
  system(glue::glue("{PLINK1} --bfile {DIR_PROCESSED}/hbv_gilead_QC --mac 1 -make-bed --out {DIR_PROCESSED}/hbv_gilead_QC_mac_1"))
  system(glue::glue("{PLINK1} --not-chr 6 --bfile {DIR_PROCESSED}/hbv_gilead_QC_mac_1 -make-bed --out {DIR_PROCESSED}/hbv_gilead_QC_SAIGE_no_chr6"))
  
  #Write out VCF File for genotype.
  system(glue::glue("{PLINK} --bfile {DIR_PROCESSED}/hbv_gilead_QC_mac_1 --mac 1 --max-maf {max_maf} --export vcf --out {DIR_PROCESSED}/hbv_gilead_QC_SAIGE"))
  system(glue::glue("./bin/bgzip -f {DIR_PROCESSED}/hbv_gilead_QC_SAIGE.vcf"))
  system(glue::glue("{BCFTOOLS} index -f {DIR_PROCESSED}/hbv_gilead_QC_SAIGE.vcf.gz"))
  
  #Make sure ordering same as fam file
  fam_file <- data.table::fread(glue::glue("{DIR_PROCESSED}/hbv_gilead_QC_mac_1.fam"),header = F)
  write(fam_file$V1,file = glue::glue("{DIR_PROCESSED}/hbv_gilead_QC_SAIGE.vcf.gz.sample"))
  
  plinkFile <- glue::glue("{DIR_PROCESSED}/hbv_gilead_QC_SAIGE_no_chr6")
  
  #Create File for LOF and MISSENSE variants
  pos_file <- data.table::fread(text = system(glue::glue("{cur_dir}/bin/bcftools-1.9/bin/bcftools query -f '%ID %CHROM %POS %REF %ALT\n' {DIR_PROCESSED}/hbv_gilead_QC_SAIGE.vcf.gz"),intern = T)) %>%
    dplyr::select(ID=V1,CHROM=V2,POS=V3,REF=V4,ALT=V5)
  
  #RUN VEP beforehand, and save as hbv_gilead_QC_SAIGE.vep
  loftee_HC <- data.table::fread(cmd = glue::glue("grep 'LoF=HC' {DIR_PROCESSED}/hbv_gilead_QC_SAIGE.vep"),header = F) %>% dplyr::select(ID=V1,GENE=V4,Type = V7) %>% 
    dplyr::distinct(ID,GENE) %>% dplyr::left_join(pos_file)
  #Write out gene set files
  unique_genes <- unique(loftee_HC$GENE)
  system(glue::glue('mkdir -p {DIR_PROCESSED}/LOF_HC_variants/'))
  out_path <- glue::glue('{DIR_PROCESSED}/LOF_HC_variants/LOF_HC_variants_')
  log <- pbmclapply(1:length(unique_genes),function(i){
    cur_df <- dplyr::filter(loftee_HC,GENE == unique_genes[i])
    if(nrow(cur_df) > 1){
      ids <- paste0(cur_df$CHR,':',cur_df$POS,'_',cur_df$REF,'/',cur_df$ALT)
      string_out <- paste0(unique_genes[i],'\t',paste(ids,collapse = '\t'))
      write(string_out,file = glue::glue("{out_path}{unique_genes[i]}.txt"),append = F)
    }
  },mc.cores = n_cores_input,ignore.interactive = T)
  
  loftee_LC <- data.table::fread(cmd = glue::glue("grep 'LoF=LC' {DIR_PROCESSED}/hbv_gilead_QC_SAIGE.vep"),header = F) %>% dplyr::select(ID=V1,GENE=V4,Type = V7) %>% 
    dplyr::distinct(ID,GENE) %>% dplyr::left_join(pos_file)
  loftee_LoF <- rbind(loftee_HC,loftee_LC) %>% dplyr::distinct(ID,GENE,.keep_all = T) 
  #Write out gene set files
  unique_genes <- unique(loftee_LoF$GENE)
  system(glue::glue('mkdir -p {DIR_PROCESSED}/LOF_variants/'))
  out_path <- glue::glue('{DIR_PROCESSED}/LOF_variants/LOF_variants_')
  log <- pbmclapply(1:length(unique_genes),function(i){
    cur_df <- dplyr::filter(loftee_LoF,GENE == unique_genes[i]) 
    if(nrow(cur_df) > 1){
      ids <- paste0(cur_df$CHR,':',cur_df$POS,'_',cur_df$REF,'/',cur_df$ALT)
      string_out <- paste0(unique_genes[i],'\t',paste(ids,collapse = '\t'))
      write(string_out,file = glue::glue("{out_path}{unique_genes[i]}.txt"),append = F)
    }
  },mc.cores = n_cores_input,ignore.interactive = T)
  
  VEP_missense <- data.table::fread(cmd = glue::glue("grep -v '#' {DIR_PROCESSED}/hbv_gilead_QC_SAIGE.vep | grep -v 'IMPACT=LOW' | grep -v 'IMPACT=MODIFIER' "),header = F)
  VEP_missense <- VEP_missense %>% dplyr::select(ID=V1,GENE=V4,Type = V7) %>% 
    dplyr::distinct(ID,GENE,.keep_all = T) %>% dplyr::left_join(pos_file)
  #Write out gene set files
  unique_genes <- unique(VEP_missense$GENE)
  system(glue::glue('mkdir -p {DIR_PROCESSED}/missense_variants/'))
  out_path <- glue::glue('{DIR_PROCESSED}/missense_variants/missense_variants_')
  log <- pbmclapply(1:length(unique_genes),function(i){
    cur_df <- dplyr::filter(VEP_missense,GENE == unique_genes[i]) 
    if(nrow(cur_df) > 1){
      ids <- paste0(cur_df$CHR,':',cur_df$POS,'_',cur_df$REF,'/',cur_df$ALT)
      string_out <- paste0(unique_genes[i],'\t',paste(ids,collapse = '\t'))
      write(string_out,file = glue::glue("{out_path}{unique_genes[i]}.txt"),append = F)
    }
  },mc.cores = n_cores_input,ignore.interactive = T)
  
  saveRDS(list(loftee_HC=loftee_HC,loftee_LC=loftee_LC,VEP_missense=VEP_missense),file = glue::glue("{DIR_PROCESSED}/variants.rds"))
  
  #Create sparse grm
  SAIGE::createSparseGRM(plinkFile = plinkFile,
                         outputPrefix = glue::glue("{DIR_PROCESSED}/sparseGRM"),
                         nThreads = n_cores)
}
RunSAIGEBurden <- function(results_path,AA_residue=NULL,LOCO=T,n_cores_input,run_missense = F,run_LoF = T){
  cur_dir <- here::here()
  load(glue::glue("{results_path}/prep-data.rda"))
  
  if(LOCO){
    plinkFile <- glue::glue("{DIR_PROCESSED}/hbv_gilead_QC")
    
  }else{
    plinkFile <- glue::glue("{DIR_PROCESSED}/hbv_gilead_QC_SAIGE_no_chr6")
  }
  
  #Generate df with same order as IDs in PLINK file
  covar <- dplyr::left_join(ids_order %>% dplyr::select(ID=X1),pc %>% dplyr::select(-FID),by=c('ID'='ID')) %>%
    dplyr::left_join(pPC,by=c('ID'='ID')) %>%
    dplyr::left_join(dat_covars_raw %>% dplyr::select(ID=host_id,AGE,SEX),by=c('ID'='ID')) %>% as.data.frame(.)
  
  #Remap Sex as binary
  covar$SEX <- ifelse(covar$SEX == 'M',0,1)
  #Make sure ordering same as fam file
  fam_file <- data.table::fread(glue::glue("{DIR_PROCESSED}/hbv_gilead_QC.fam"),header = F)
  covar <- covar[match(fam_file$V1,covar$ID),]
  dat_pathogen <- dat_pathogen[match(fam_file$V1,dat_pathogen$ID),]
  dat_pathogen_AA_only <- dat_pathogen[,-1]
  #Write Pheno fIle
  pheno_file_path <- glue::glue('{DIR_PROCESSED}/hbv_gilead_QC_SAIGE_pheno.txt')
  data.table::fwrite(cbind(covar,dat_pathogen_AA_only),file = pheno_file_path,sep = '\t',col.names = T,na = 'NA',quote = F)
  
  cur_AA <- AA_residue
  #Prepare NULL Matrix
  if(!(glue::glue('{DIR_SCRATCH}/SAIGE_LOCO_{LOCO}_{cur_AA}.varianceRatio.txt') %in% dir(DIR_SCRATCH))){
    nullGLMM_stat <- tryCatch(expr = {
      SAIGE::fitNULLGLMM(plinkFile = plinkFile,
                         phenoFile = pheno_file_path,
                         phenoCol = cur_AA,
                         traitType = 'binary',
                         invNormalize = F,
                         covarColList = setdiff(colnames(covar),c("ID")),
                         sampleIDColinphenoFile = "ID",
                         LOCO = LOCO,
                         nThreads = n_cores_input,
                         IsSparseKin = T,
                         sparseGRMFile = glue::glue("{DIR_PROCESSED}/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx"),
                         sparseGRMSampleIDFile = glue::glue("{DIR_PROCESSED}/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"),
                         isCateVarianceRatio = T,
                         outputPrefix = glue::glue('{DIR_SCRATCH}/SAIGE_LOCO_{LOCO}_{cur_AA}'))
      T
    },error = function(e){
      write("NULL MODEL ERROR",glue::glue('{DIR_SCRATCH}/SAIGE_LOCO_{LOCO}_{cur_AA}.error'),append = T)
      write(as.character(e),glue::glue('{DIR_SCRATCH}/SAIGE_LOCO_{LOCO}_{cur_AA}.error'),append = T)
      return(F)},
    warning = function(w) {
      write("NULL MODEL WARNING",glue::glue('{DIR_SCRATCH}/SAIGE_LOCO_{LOCO}_{cur_AA}.warning'),append = T)
      write(as.character(w),glue::glue('{DIR_SCRATCH}/SAIGE_LOCO_{LOCO}_{cur_AA}.error'),append = T)
      return(F)})
    if(!nullGLMM_stat){
      return()
    }
  }
  if(run_LoF & !(glue::glue('{DIR_SCRATCH}/SAIGE_LOCO_{LOCO}_{cur_AA}.error') %in% dir(DIR_SCRATCH))){
    #Run LOF Variants
    system(glue::glue("mkdir -p {DIR_SCRATCH}/SAIGE_BURDEN_LOF_LOCO_{LOCO}_{cur_AA}/"))
    lof_files <- setdiff(dir(glue::glue("{DIR_PROCESSED}/LOF_variants/")),"LOF_variants.txt")
    
    log <- mclapply(1:length(lof_files),function(i){
      cur_prefix <- gsub(pattern = '.txt',replacement = '',x = glue::glue('{DIR_SCRATCH}/SAIGE_BURDEN_LOF_LOCO_{LOCO}_{cur_AA}/{lof_files[i]}'))
      SPA_test_stat <- tryCatch(expr = {
        log_file <- file(glue::glue('{cur_prefix}.log'))
        sink(file = log_file)
        sink(file = log_file,type = 'message')
        SAIGE::SPAGMMATtest(vcfFile = glue::glue("{DIR_PROCESSED}/hbv_gilead_QC_SAIGE.vcf.gz"),
                            vcfFileIndex = glue::glue("{DIR_PROCESSED}/hbv_gilead_QC_SAIGE.vcf.gz.tbi"),
                            vcfField = "GT",
                            minMAC = 0.5,
                            chrom = '1',
                            IsDropMissingDosages = F,
                            IsOutputAFinCaseCtrl = T,
                            IsOutputNinCaseCtrl = T,
                            minMAF = 0,
                            LOCO = LOCO,
                            groupFile = glue::glue(glue::glue('{DIR_PROCESSED}/LOF_variants/{lof_files[i]}')),
                            numLinesOutput = .Machine$integer.max,
                            sampleFile= glue::glue("{DIR_PROCESSED}/hbv_gilead_QC_SAIGE.vcf.gz.sample"),
                            GMMATmodelFile= glue::glue('{DIR_SCRATCH}/SAIGE_LOCO_{LOCO}_{cur_AA}.rda'),
                            varianceRatioFile=glue::glue('{DIR_SCRATCH}/SAIGE_LOCO_{LOCO}_{cur_AA}.varianceRatio.txt'),
                            sparseSigmaFile = glue::glue("{DIR_SCRATCH}/SAIGE_LOCO_{LOCO}_{cur_AA}.varianceRatio.txt_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseSigma.mtx"),
                            SAIGEOutputFile= glue::glue("{cur_prefix}.out"))
        sink()
        sink(type = 'message')
        close(log_file)
        rm(log_file)
        T
      },error = function(e){
        write(glue::glue("ERROR"),glue::glue('{cur_prefix}.error'),append = T)
        write(as.character(e),glue::glue('{cur_prefix}.error'),append = T)
        sink()
        sink(type = 'message')
        close(log_file)
        rm(log_file)
      },
      warning = function(w){
        write(glue::glue("WARNING"),glue::glue('{cur_prefix}.warning'),append = T)
        write(as.character(w),glue::glue('{cur_prefix}.warning'),append = T)
        sink()
        sink(type = 'message')
        close(log_file)
        rm(log_file)
      })
    },mc.cores = n_cores_input,mc.silent = T,mc.preschedule = F)
    # system(glue::glue("find {DIR_SCRATCH}/SAIGE_BURDEN_LOF_LOCO_{LOCO}_{cur_AA}/ -name *.out | xargs cat | head -n 1 > {DIR_SCRATCH}/SAIGE_BURDEN_LOF_LOCO_{LOCO}_{cur_AA}.txt"))
    # system(glue::glue("find {DIR_SCRATCH}/SAIGE_BURDEN_LOF_LOCO_{LOCO}_{cur_AA}/ -name *.out | xargs cat | grep -v 'SKAT' >> {DIR_SCRATCH}/SAIGE_BURDEN_LOF_LOCO_{LOCO}_{cur_AA}.txt"))
    # system(glue::glue("find {DIR_SCRATCH}/SAIGE_BURDEN_LOF_LOCO_{LOCO}_{cur_AA}/ -name *.out | xargs rm"))
    
    #Run LOF_HC Variants
    system(glue::glue("mkdir -p {DIR_SCRATCH}/SAIGE_BURDEN_LOF_HC_LOCO_{LOCO}_{cur_AA}/"))
    LOF_HC_files <- setdiff(dir(glue::glue("{DIR_PROCESSED}/LOF_HC_variants/")),"LOF_HC_variants.txt")
    
    log <- mclapply(1:length(LOF_HC_files),function(i){
      cur_prefix <- gsub(pattern = '.txt',replacement = '',x = glue::glue('{DIR_SCRATCH}/SAIGE_BURDEN_LOF_HC_LOCO_{LOCO}_{cur_AA}/{LOF_HC_files[i]}'))
      SPA_test_stat <- tryCatch(expr = {
        log_file <- file(glue::glue('{cur_prefix}.log'))
        sink(file = log_file)
        sink(file = log_file,type = 'message')
        SAIGE::SPAGMMATtest(vcfFile = glue::glue("{DIR_PROCESSED}/hbv_gilead_QC_SAIGE.vcf.gz"),
                            vcfFileIndex = glue::glue("{DIR_PROCESSED}/hbv_gilead_QC_SAIGE.vcf.gz.tbi"),
                            vcfField = "GT",
                            chrom = '1',
                            minMAC = 0.5,
                            IsDropMissingDosages = F,
                            IsOutputAFinCaseCtrl = T,
                            IsOutputNinCaseCtrl = T,
                            minMAF = 0,
                            LOCO = LOCO,
                            groupFile = glue::glue(glue::glue('{DIR_PROCESSED}/LOF_HC_variants/{LOF_HC_files[i]}')),
                            numLinesOutput = .Machine$integer.max,
                            sampleFile= glue::glue("{DIR_PROCESSED}/hbv_gilead_QC_SAIGE.vcf.gz.sample"),
                            GMMATmodelFile= glue::glue('{DIR_SCRATCH}/SAIGE_LOCO_{LOCO}_{cur_AA}.rda'),
                            varianceRatioFile=glue::glue('{DIR_SCRATCH}/SAIGE_LOCO_{LOCO}_{cur_AA}.varianceRatio.txt'),
                            sparseSigmaFile = glue::glue("{DIR_SCRATCH}/SAIGE_LOCO_{LOCO}_{cur_AA}.varianceRatio.txt_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseSigma.mtx"),
                            SAIGEOutputFile= glue::glue("{cur_prefix}.out"))
        sink()
        sink(type = 'message')
        close(log_file)
        rm(log_file)
        T
      },error = function(e){
        write(glue::glue("ERROR"),glue::glue('{cur_prefix}.error'),append = T)
        write(as.character(e),glue::glue('{cur_prefix}.error'),append = T)
        sink()
        sink(type = 'message')
        close(log_file)
        rm(log_file)
      },
      warning = function(w){
        write(glue::glue("WARNING"),glue::glue('{cur_prefix}.warning'),append = T)
        write(as.character(w),glue::glue('{cur_prefix}.warning'),append = T)
        sink()
        sink(type = 'message')
        close(log_file)
        rm(log_file)
      })
    },mc.cores = n_cores_input,mc.silent = T,mc.preschedule = F)
    # system(glue::glue("find {DIR_SCRATCH}/SAIGE_BURDEN_LOF_HC_LOCO_{LOCO}_{cur_AA}/ -name *.out | xargs cat | head -n 1 > {DIR_SCRATCH}/SAIGE_BURDEN_LOF_HC_LOCO_{LOCO}_{cur_AA}.txt"))
    # system(glue::glue("find {DIR_SCRATCH}/SAIGE_BURDEN_LOF_HC_LOCO_{LOCO}_{cur_AA}/ -name *.out | xargs cat | grep -v 'SKAT' >> {DIR_SCRATCH}/SAIGE_BURDEN_LOF_HC_LOCO_{LOCO}_{cur_AA}.txt"))
    # system(glue::glue("find {DIR_SCRATCH}/SAIGE_BURDEN_LOF_HC_LOCO_{LOCO}_{cur_AA}/ -name *.out | xargs rm"))
    
  }
  #Run Missense Variants
  if(run_missense & !(glue::glue('{DIR_SCRATCH}/SAIGE_LOCO_{LOCO}_{cur_AA}.error') %in% dir(DIR_SCRATCH))){
    system(glue::glue("mkdir -p {DIR_SCRATCH}/SAIGE_BURDEN_MISSENSE_LOCO_{LOCO}_{cur_AA}/"))
    missense_files <- setdiff(dir(glue::glue("{DIR_PROCESSED}/missense_variants/")),"missense_variants.txt")

    log <- mclapply(1:length(missense_files),function(i) {
      cur_prefix <- gsub(pattern = '.txt',replacement = '',x = glue::glue('{DIR_SCRATCH}/SAIGE_BURDEN_MISSENSE_LOCO_{LOCO}_{cur_AA}/{missense_files[i]}'))
      SPA_test_stat <- tryCatch(expr = {
        log_file <- file(glue::glue('{cur_prefix}.log'))
        sink(file = log_file)
        sink(file = log_file,type = 'message')
        SAIGE::SPAGMMATtest(vcfFile = glue::glue("{DIR_PROCESSED}/hbv_gilead_QC_SAIGE.vcf.gz"),
                            vcfFileIndex = glue::glue("{DIR_PROCESSED}/hbv_gilead_QC_SAIGE.vcf.gz.tbi"),
                            vcfField = "GT",
                            minMAC = 0.5,
                            chrom = '1',
                            IsDropMissingDosages = F,
                            IsOutputAFinCaseCtrl = T,
                            IsOutputNinCaseCtrl = T,
                            minMAF = 0,
                            LOCO = LOCO,
                            groupFile = glue::glue(glue::glue('{DIR_PROCESSED}/missense_variants/{missense_files[i]}')),
                            numLinesOutput = .Machine$integer.max,
                            sampleFile= glue::glue("{DIR_PROCESSED}/hbv_gilead_QC_SAIGE.vcf.gz.sample"),
                            GMMATmodelFile= glue::glue('{DIR_SCRATCH}/SAIGE_LOCO_{LOCO}_{cur_AA}.rda'),
                            varianceRatioFile=glue::glue('{DIR_SCRATCH}/SAIGE_LOCO_{LOCO}_{cur_AA}.varianceRatio.txt'),
                            sparseSigmaFile = glue::glue("{DIR_SCRATCH}/SAIGE_LOCO_{LOCO}_{cur_AA}.varianceRatio.txt_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseSigma.mtx"),
                            SAIGEOutputFile= glue::glue("{cur_prefix}.out"))
        sink()
        sink(type = 'message')
        close(log_file)
        rm(log_file)
      },error = function(e){
        write(glue::glue("ERROR"),glue::glue('{cur_prefix}.error'),append = T)
        write(as.character(e),glue::glue('{cur_prefix}.error'),append = T)
        sink()
        sink(type = 'message')
        close(log_file)
        rm(log_file)
      },
      warning = function(w){
        write(glue::glue("WARNING"),glue::glue('{cur_prefix}.warning'),append = T)
        write(as.character(w),glue::glue('{cur_prefix}.warning'),append = T)
        sink()
        sink(type = 'message')
        close(log_file)
        rm(log_file)
      })
    },mc.cores = n_cores_input,mc.silent = T,mc.preschedule = F)
    
    # system(glue::glue("find {DIR_SCRATCH}/SAIGE_BURDEN_MISSENSE_LOCO_{LOCO}_{cur_AA}/ -name *.out | xargs cat | head -n 1 > {DIR_SCRATCH}/SAIGE_BURDEN_MISSENSE_LOCO_{LOCO}_{cur_AA}.txt"))
    # system(glue::glue("find {DIR_SCRATCH}/SAIGE_BURDEN_MISSENSE_LOCO_{LOCO}_{cur_AA}/ -name *.out | xargs cat | grep -v 'SKAT' >> {DIR_SCRATCH}/SAIGE_BURDEN_MISSENSE_LOCO_{LOCO}_{cur_AA}.txt"))
    # system(glue::glue("find {DIR_SCRATCH}/SAIGE_BURDEN_MISSENSE_LOCO_{LOCO}_{cur_AA}/ -name *.out | xargs rm"))
  }
}