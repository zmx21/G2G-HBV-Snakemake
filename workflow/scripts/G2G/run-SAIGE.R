library(SAIGE)
CreateGRMGenoFile <- function(results_path){
  cur_dir <- here::here()
  load(glue::glue("{results_path}/prep-data.rda"))
  
  #Write out PLINK files without Chr 6
  system(glue::glue("{PLINK1} --not-chr 6 --bfile {DIR_PROCESSED}/hbv_gilead_QC -make-bed --out {DIR_PROCESSED}/hbv_gilead_QC_SAIGE_no_chr6"))
  
  #Write out VCF File for genotype.
  system(glue::glue("{PLINK} --bfile {DIR_PROCESSED}/hbv_gilead_QC --export vcf --out {DIR_PROCESSED}/hbv_gilead_QC_SAIGE"))
  system(glue::glue("./bin/bgzip -f {DIR_PROCESSED}/hbv_gilead_QC_SAIGE.vcf"))
  system(glue::glue("{BCFTOOLS} index -f -t {DIR_PROCESSED}/hbv_gilead_QC_SAIGE.vcf.gz"))
  
  #Make sure ordering same as fam file
  fam_file <- data.table::fread(glue::glue("{DIR_PROCESSED}/hbv_gilead_QC.fam"),header = F)
  write(fam_file$V1,file = glue::glue("{DIR_PROCESSED}/hbv_gilead_QC_SAIGE.vcf.gz.sample"))
  
}

RunSAIGE <- function(results_path,AA_residue=NULL,LOCO=T,n_cores_input = 1){
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
  
  if(is.null(AA_residue)){
    for(i in 1:ncol(dat_pathogen_AA_only)){
      cur_AA <- colnames(dat_pathogen_AA_only)[i]
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
                           outputPrefix = glue::glue('{DIR_SCRATCH}/SAIGE_LOCO_{LOCO}_{cur_AA}'))
        return(T)
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
      err_counter <- mclapply(1:22,function(j){
        SPA_test_stat <- tryCatch(expr = {
          SAIGE::SPAGMMATtest(vcfFile = glue::glue("{DIR_PROCESSED}/hbv_gilead_QC_SAIGE.vcf.gz"),
                              vcfFileIndex = glue::glue("{DIR_PROCESSED}/hbv_gilead_QC_SAIGE.vcf.gz.tbi"),
                              vcfField = "GT",
                              chrom = as.character(j),
                              minMAC = 4,
                              IsDropMissingDosages = T,
                              IsOutputAFinCaseCtrl = T,
                              IsOutputNinCaseCtrl = T,
                              minMAF = 0.05,
                              LOCO = LOCO,
                              numLinesOutput = .Machine$integer.max,
                              sampleFile= glue::glue("{DIR_PROCESSED}/hbv_gilead_QC_SAIGE.vcf.gz.sample"),
                              GMMATmodelFile= glue::glue('{DIR_SCRATCH}/SAIGE_LOCO_{LOCO}_{cur_AA}.rda'),
                              varianceRatioFile=glue::glue('{DIR_SCRATCH}/SAIGE_LOCO_{LOCO}_{cur_AA}.varianceRatio.txt'),
                              SAIGEOutputFile= glue::glue('{DIR_SCRATCH}/SAIGE_SPA_TEST_LOCO_{LOCO}_chr{j}_{cur_AA}'))
          return(T)}
          ,error = function(e){
            write(glue::glue("CHR {j} ERROR"),glue::glue('{DIR_SCRATCH}/SAIGE_LOCO_{LOCO}_{cur_AA}.error'),append = T)
            write(as.character(e),glue::glue('{DIR_SCRATCH}/SAIGE_LOCO_{LOCO}_{cur_AA}.error'),append = T)
            return(F)},
          warning = function(w){
            write(glue::glue("CHR {j} WARNING"),glue::glue('{DIR_SCRATCH}/SAIGE_LOCO_{LOCO}_{cur_AA}.warning'),append = T)
            write(as.character(w),glue::glue('{DIR_SCRATCH}/SAIGE_LOCO_{LOCO}_{cur_AA}.warning'),append = T)
            return(F)})
        return(SPA_test_stat)
      },mc.preschedule = F,mc.cores = n_cores_input)
      err_counter <- unlist(err_counter)
      if(sum(err_counter) == 22){
        system(glue::glue("head -n 1 {DIR_SCRATCH}/SAIGE_SPA_TEST_LOCO_{LOCO}_chr{which(!err_counter)[1]}_{cur_AA} > {DIR_SCRATCH}/SAIGE_SPA_TEST_LOCO_{LOCO}_{cur_AA}"))
        system(glue::glue("cat {DIR_SCRATCH}/SAIGE_SPA_TEST_LOCO_{LOCO}_chr*_{cur_AA} | grep -v 'CHR' >> {DIR_SCRATCH}/SAIGE_SPA_TEST_LOCO_{LOCO}_{cur_AA}"))
        system(glue::glue("rm {DIR_SCRATCH}/SAIGE_SPA_TEST_LOCO_{LOCO}_chr*_{cur_AA}"))
      }else{
        return()
      }
    }
  }else{
    cur_AA <- AA_residue
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
    err_counter <- mclapply(1:22,function(j){
      SPA_test_stat <- tryCatch(expr = {
        SAIGE::SPAGMMATtest(vcfFile = glue::glue("{DIR_PROCESSED}/hbv_gilead_QC_SAIGE.vcf.gz"),
                            vcfFileIndex = glue::glue("{DIR_PROCESSED}/hbv_gilead_QC_SAIGE.vcf.gz.tbi"),
                            vcfField = "GT",
                            chrom = as.character(j),
                            minMAC = 4,
                            IsDropMissingDosages = T,
                            IsOutputAFinCaseCtrl = T,
                            IsOutputNinCaseCtrl = T,
                            minMAF = 0.05,
                            LOCO = LOCO,
                            numLinesOutput = .Machine$integer.max,
                            sampleFile= glue::glue("{DIR_PROCESSED}/hbv_gilead_QC_SAIGE.vcf.gz.sample"),
                            GMMATmodelFile= glue::glue('{DIR_SCRATCH}/SAIGE_LOCO_{LOCO}_{cur_AA}.rda'),
                            varianceRatioFile=glue::glue('{DIR_SCRATCH}/SAIGE_LOCO_{LOCO}_{cur_AA}.varianceRatio.txt'),
                            SAIGEOutputFile= glue::glue('{DIR_SCRATCH}/SAIGE_SPA_TEST_LOCO_{LOCO}_chr{j}_{cur_AA}'))
        T
      }
      ,error = function(e){
        write(glue::glue("CHR {j} ERROR"),glue::glue('{DIR_SCRATCH}/SAIGE_LOCO_{LOCO}_{cur_AA}.error'),append = T)
        write(as.character(e),glue::glue('{DIR_SCRATCH}/SAIGE_LOCO_{LOCO}_{cur_AA}.error'),append = T)
        return(F)},
      warning = function(w){
        write(glue::glue("CHR {j} WARNING"),glue::glue('{DIR_SCRATCH}/SAIGE_LOCO_{LOCO}_{cur_AA}.warning'),append = T)
        write(as.character(w),glue::glue('{DIR_SCRATCH}/SAIGE_LOCO_{LOCO}_{cur_AA}.warning'),append = T)
        return(F)})
      return(SPA_test_stat)
    },mc.preschedule = F,mc.cores = n_cores_input,mc.silent = T)
    err_counter <- unlist(err_counter)
    if(sum(err_counter) == 22){
      system(glue::glue("head -n 1 {DIR_SCRATCH}/SAIGE_SPA_TEST_LOCO_{LOCO}_chr{which(err_counter)[1]}_{cur_AA} > {DIR_SCRATCH}/SAIGE_SPA_TEST_LOCO_{LOCO}_{cur_AA}"))
      system(glue::glue("cat {DIR_SCRATCH}/SAIGE_SPA_TEST_LOCO_{LOCO}_chr*_{cur_AA} | grep -v 'CHR' >> {DIR_SCRATCH}/SAIGE_SPA_TEST_LOCO_{LOCO}_{cur_AA}"))
      system(glue::glue("rm {DIR_SCRATCH}/SAIGE_SPA_TEST_LOCO_{LOCO}_chr*_{cur_AA}"))
    }else{
      return()
    }
  }
}

RunSAIGEPheno <- function(results_path,cur_outcome,LOCO=T,n_cores_input = 1){
  cur_dir <- here::here()
  load(glue::glue("{results_path}/prep-data.rda"))
  system(glue::glue("mkdir -p {DIR_SCRATCH}/GWAS/"))
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
  dat_alternative_outcome <- dat_alternative_outcome[match(fam_file$V1,dat_alternative_outcome$host_id),]
  #Write Pheno fIle
  pheno_file_path <- glue::glue('{DIR_PROCESSED}/hbv_gilead_QC_SAIGE_pheno_alt.txt')
  data.table::fwrite(dplyr::inner_join(covar,dat_alternative_outcome,by=c('ID' = 'host_id')),file = pheno_file_path,sep = '\t',col.names = T,na = 'NA',quote = F)
  if(length(unique(as.vector(t(dat_alternative_outcome[,cur_outcome])))) == 2){
    mdl_type = 'binary'
  }else{
    mdl_type = 'quantitative'
  }
  nullGLMM_stat <- tryCatch(expr = {
    SAIGE::fitNULLGLMM(plinkFile = plinkFile,
                       phenoFile = pheno_file_path,
                       phenoCol = cur_outcome,
                       traitType = mdl_type,
                       invNormalize = F,
                       covarColList = setdiff(colnames(covar),c("ID")),
                       sampleIDColinphenoFile = "ID",
                       LOCO = LOCO,
                       nThreads = n_cores_input,
                       outputPrefix = glue::glue('{DIR_SCRATCH}/GWAS/SAIGE_LOCO_{LOCO}_{cur_outcome}'))
    T
  },error = function(e){
    write("NULL MODEL ERROR",glue::glue('{DIR_SCRATCH}/GWAS/SAIGE_LOCO_{LOCO}_{cur_outcome}.error'),append = T)
    write(as.character(e),glue::glue('{DIR_SCRATCH}/GWAS/SAIGE_LOCO_{LOCO}_{cur_outcome}.error'),append = T)
    return(F)},
  warning = function(w) {
    write("NULL MODEL WARNING",glue::glue('{DIR_SCRATCH}/GWAS/SAIGE_LOCO_{LOCO}_{cur_outcome}.warning'),append = T)
    write(as.character(w),glue::glue('{DIR_SCRATCH}/GWAS/SAIGE_LOCO_{LOCO}_{cur_outcome}.error'),append = T)
    return(F)})
  if(!nullGLMM_stat){
    return()
  }
  err_counter <- mclapply(1:22,function(j){
    SPA_test_stat <- tryCatch(expr = {
      SAIGE::SPAGMMATtest(vcfFile = glue::glue("{DIR_PROCESSED}/hbv_gilead_QC_SAIGE.vcf.gz"),
                          vcfFileIndex = glue::glue("{DIR_PROCESSED}/hbv_gilead_QC_SAIGE.vcf.gz.tbi"),
                          vcfField = "GT",
                          chrom = as.character(j),
                          minMAC = 4,
                          IsDropMissingDosages = T,
                          IsOutputAFinCaseCtrl = T,
                          IsOutputNinCaseCtrl = T,
                          minMAF = 0.05,
                          LOCO = LOCO,
                          numLinesOutput = .Machine$integer.max,
                          sampleFile= glue::glue("{DIR_PROCESSED}/hbv_gilead_QC_SAIGE.vcf.gz.sample"),
                          GMMATmodelFile= glue::glue('{DIR_SCRATCH}/GWAS/SAIGE_LOCO_{LOCO}_{cur_outcome}.rda'),
                          varianceRatioFile=glue::glue('{DIR_SCRATCH}/GWAS/SAIGE_LOCO_{LOCO}_{cur_outcome}.varianceRatio.txt'),
                          SAIGEOutputFile= glue::glue('{DIR_SCRATCH}/GWAS/SAIGE_SPA_TEST_LOCO_{LOCO}_chr{j}_{cur_outcome}'))
      T
    }
    ,error = function(e){
      write(glue::glue("CHR {j} ERROR"),glue::glue('{DIR_SCRATCH}/GWAS/SAIGE_LOCO_{LOCO}_{cur_outcome}.error'),append = T)
      write(as.character(e),glue::glue('{DIR_SCRATCH}/GWAS/SAIGE_LOCO_{LOCO}_{cur_outcome}.error'),append = T)
      return(F)},
    warning = function(w){
      write(glue::glue("CHR {j} WARNING"),glue::glue('{DIR_SCRATCH}/GWAS/SAIGE_LOCO_{LOCO}_{cur_outcome}.warning'),append = T)
      write(as.character(w),glue::glue('{DIR_SCRATCH}/GWAS/SAIGE_LOCO_{LOCO}_{cur_outcome}.warning'),append = T)
      return(F)})
    return(SPA_test_stat)
  },mc.preschedule = F,mc.cores = n_cores_input,mc.silent = T)
  err_counter <- unlist(err_counter)
  if(sum(err_counter) == 22){
    system(glue::glue("head -n 1 {DIR_SCRATCH}/GWAS/SAIGE_SPA_TEST_LOCO_{LOCO}_chr{which(err_counter)[1]}_{cur_outcome} > {DIR_SCRATCH}/GWAS/SAIGE_SPA_TEST_LOCO_{LOCO}_{cur_outcome}"))
    system(glue::glue("cat {DIR_SCRATCH}/GWAS/SAIGE_SPA_TEST_LOCO_{LOCO}_chr*_{cur_outcome} | grep -v 'CHR' >> {DIR_SCRATCH}/GWAS/SAIGE_SPA_TEST_LOCO_{LOCO}_{cur_outcome}"))
    system(glue::glue("rm {DIR_SCRATCH}/GWAS/SAIGE_SPA_TEST_LOCO_{LOCO}_chr*_{cur_outcome}"))
  }else{
    return()
  }
  
}
RunSAIGE_HLA <- function(results_path,AA_residue=NULL,LOCO=F,n_cores_input = 1,OUT_DIR,PLINK_HLA,run_AA){
  load(glue::glue("{results_path}/prep-data.rda"))
  #Write HLA Plink file to VCF
  #Write out VCF File for genotype.
  cur_dir <- here::here()
  system(paste0("awk '{ print $1\"\t\"$1 }' < ",DIR_PROCESSED,"/hbv_gilead_QC_SAIGE.vcf.gz.sample > ",OUT_DIR,'hla_samples.txt'))
  if(!run_AA){
    system(glue::glue("{cur_dir}/bin/plink2_linux --bfile {PLINK_HLA} --chr 6 --from-bp 28477797 --to-bp 33448354 --geno 0.2 --keep {OUT_DIR}hla_samples.txt --export vcf --out {PLINK_HLA}"))
    system(glue::glue("{cur_dir}/bin/bcftools-1.9/bin/bcftools reheader -s {DIR_PROCESSED}/hbv_gilead_QC_SAIGE.vcf.gz.sample {PLINK_HLA}.vcf > {PLINK_HLA}.rehead.vcf"))
    
  }else{
    system(glue::glue("{cur_dir}/bin/plink2_linux --bfile {PLINK_HLA} --chr 6 --from-bp 28477797 --to-bp 33448354  --keep {OUT_DIR}hla_samples.txt --export vcf --out {PLINK_HLA}"))
    system(glue::glue("{cur_dir}/bin/bcftools-1.9/bin/bcftools reheader -s {DIR_PROCESSED}/hbv_gilead_QC_SAIGE.vcf.gz.sample {PLINK_HLA}.vcf > {PLINK_HLA}.rehead.vcf"))
  }
  system(glue::glue("/home/zmxu/Software/bgzip -f {PLINK_HLA}.rehead.vcf"))
  system(glue::glue("{cur_dir}/bin/bcftools-1.9/bin/bcftools index -f -t {PLINK_HLA}.rehead.vcf.gz"))
  
  
  #Generate df with same order as IDs in PLINK file
  covar <- dplyr::left_join(ids_order %>% dplyr::select(ID=X1),pc %>% dplyr::select(-FID),by=c('ID'='ID')) %>%
    dplyr::left_join(pPC,by=c('ID'='ID')) %>%
    dplyr::left_join(dat_covars_raw %>% dplyr::select(ID=host_id,AGE,SEX),by=c('ID'='ID')) %>% as.data.frame(.)
  
  #Remap Sex as binary
  covar$SEX <- ifelse(covar$SEX == 'M',0,1)
  
  #Use existing NULL model file to run new single variant analysis
  tryCatch(expr = {SAIGE::fitNULLGLMM(plinkFile = glue::glue("{DIR_PROCESSED}/hbv_gilead_QC_SAIGE_no_chr6"),
                                      phenoFile = glue::glue('{DIR_PROCESSED}/hbv_gilead_QC_SAIGE_pheno.txt'),
                                      phenoCol = AA_residue,
                                      traitType = 'binary',
                                      invNormalize = F,
                                      covarColList = setdiff(colnames(covar),c("ID")),
                                      sampleIDColinphenoFile = "ID",
                                      LOCO = LOCO,
                                      nThreads = n_cores_input,
                                      outputPrefix = glue::glue('{DIR_SCRATCH}/SAIGE_LOCO_{LOCO}_{AA_residue}'))
  },error = function(e) print(e))
  
  tryCatch(expr = {SAIGE::SPAGMMATtest(vcfFile = glue::glue("{PLINK_HLA}.rehead.vcf.gz"),
                                       vcfFileIndex = glue::glue("{PLINK_HLA}.rehead.vcf.gz.tbi"),
                                       vcfField = "GT",
                                       chrom = '6',
                                       minMAC = 4,
                                       IsDropMissingDosages = T,
                                       IsOutputAFinCaseCtrl = T,
                                       IsOutputNinCaseCtrl = T,
                                       minMAF = 0.01,
                                       LOCO = LOCO,
                                       numLinesOutput = .Machine$integer.max,
                                       sampleFile= glue::glue("{DIR_PROCESSED}/hbv_gilead_QC_SAIGE.vcf.gz.sample"),
                                       GMMATmodelFile= glue::glue('{DIR_SCRATCH}/SAIGE_LOCO_{LOCO}_{AA_residue}.rda'),
                                       varianceRatioFile=glue::glue('{DIR_SCRATCH}/SAIGE_LOCO_{LOCO}_{AA_residue}.varianceRatio.txt'),
                                       SAIGEOutputFile= glue::glue('{PLINK_HLA}_HBV_{LOCO}_{AA_residue}'))},error = function(e) return(F))
  saige_result_no_cond <- data.table::fread(glue::glue("{PLINK_HLA}_HBV_{LOCO}_{AA_residue}"))
  prev_snp_result <- saige_result_no_cond
  prev_aa_result <- saige_result_no_cond
  
  max_iter <- 3
  cur_iter <- 1
  min_snp_nuid <- c()
  min_snp_id<- c()
  min_add_nuid <- c()
  min_aa_id <- c()
  while (cur_iter <= max_iter) {
    if(!is.null(prev_snp_result)){
      #Condition on top SNP
      snp_results <- dplyr::filter(prev_snp_result,grepl('rs',SNPID))
      if(cur_iter == 1){
        min_snp <- snp_results[which.min(snp_results$p.value),]
      }else{
        min_snp <- snp_results[which.min(snp_results$p.value_cond),]
      }
      min_snp_nuid <- c(min_snp_nuid,paste0(min_snp$CHR,':',min_snp$POS,"_",min_snp$Allele1,'/',min_snp$Allele2))
      min_snp_id <- c(min_snp_id,min_snp$SNPID)
      min_snp_id_string <- paste0(min_snp_id,collapse = '_')
      tryCatch(expr = {SAIGE::SPAGMMATtest(vcfFile = glue::glue("{PLINK_HLA}.rehead.vcf.gz"),
                                           vcfFileIndex = glue::glue("{PLINK_HLA}.rehead.vcf.gz.tbi"),
                                           vcfField = "GT",
                                           chrom = '6',
                                           minMAC = 4,
                                           IsDropMissingDosages = T,
                                           IsOutputAFinCaseCtrl = T,
                                           IsOutputNinCaseCtrl = T,
                                           minMAF = 0.01,
                                           LOCO = LOCO,
                                           numLinesOutput = .Machine$integer.max,
                                           condition = min_snp_nuid,
                                           sampleFile= glue::glue("{DIR_PROCESSED}/hbv_gilead_QC_SAIGE.vcf.gz.sample"),
                                           GMMATmodelFile= glue::glue('{DIR_SCRATCH}/SAIGE_LOCO_{LOCO}_{AA_residue}.rda'),
                                           varianceRatioFile=glue::glue('{DIR_SCRATCH}/SAIGE_LOCO_{LOCO}_{AA_residue}.varianceRatio.txt'),
                                           SAIGEOutputFile= glue::glue('{PLINK_HLA}_HBV_{LOCO}_{AA_residue}_cond_snp_{min_snp_id_string}'))},error = function(e) print(e))
      if(!file.exists(glue::glue('{PLINK_HLA}_HBV_{LOCO}_{AA_residue}_cond_snp_{min_snp_id_string}'))){
        prev_snp_result <- NULL
      }else{
        prev_snp_result <- data.table::fread(glue::glue('{PLINK_HLA}_HBV_{LOCO}_{AA_residue}_cond_snp_{min_snp_id_string}'))
      }
    }
    if(run_AA & !is.null(prev_aa_result)){
      aa_results <- dplyr::filter(prev_aa_result,grepl('AA',SNPID)) 
      if(cur_iter == 1){
        min_aa <- aa_results[which.min(aa_results$p.value),]
      }else{
        min_aa <- aa_results[which.min(aa_results$p.value_cond),]
      }
      min_add_nuid <- c(min_add_nuid,paste0(min_aa$CHR,':',min_aa$POS,"_",min_aa$Allele1,'/',min_aa$Allele2))
      min_aa_id <- c(min_aa_id,min_aa$SNPID)
      min_add_id_string <- paste0(min_aa_id,collapse = '_')
      tryCatch(expr = {SAIGE::SPAGMMATtest(vcfFile = glue::glue("{PLINK_HLA}.rehead.vcf.gz"),
                                           vcfFileIndex = glue::glue("{PLINK_HLA}.rehead.vcf.gz.tbi"),
                                           vcfField = "GT",
                                           chrom = '6',
                                           minMAC = 4,
                                           IsDropMissingDosages = T,
                                           IsOutputAFinCaseCtrl = T,
                                           IsOutputNinCaseCtrl = T,
                                           minMAF = 0.01,
                                           LOCO = LOCO,
                                           numLinesOutput = .Machine$integer.max,
                                           condition = min_add_nuid,
                                           sampleFile= glue::glue("{DIR_PROCESSED}/hbv_gilead_QC_SAIGE.vcf.gz.sample"),
                                           GMMATmodelFile= glue::glue('{DIR_SCRATCH}/SAIGE_LOCO_{LOCO}_{AA_residue}.rda'),
                                           varianceRatioFile=glue::glue('{DIR_SCRATCH}/SAIGE_LOCO_{LOCO}_{AA_residue}.varianceRatio.txt'),
                                           SAIGEOutputFile= glue::glue('{PLINK_HLA}_HBV_{LOCO}_{AA_residue}_cond_aa_{min_add_id_string}'))},error = function(e) print(e))
      if(!file.exists(glue::glue('{PLINK_HLA}_HBV_{LOCO}_{AA_residue}_cond_aa_{min_add_id_string}'))){
        prev_snp_result <- NULL
      }else{
        prev_snp_result <- data.table::fread(glue::glue('{PLINK_HLA}_HBV_{LOCO}_{AA_residue}_cond_aa_{min_add_id_string}'))
      }
    }
    cur_iter <- cur_iter + 1
  }
}



RunSAIGE_HLA_PyHLA <- function(results_path,AA_residue=NULL,LOCO=F,n_cores_input = 1,out_dir,vcf_path,create_vcf = F){
  load(glue::glue("{results_path}/prep-data.rda"))
  vcf_path_alllele <- gsub(pattern = 'assocAA',replacement = 'assocAllele',x = vcf_path)
  #Write HLA Plink file to VCF
  #Write out VCF File for genotype.
  cur_dir <- here::here()
  if(create_vcf){
    system(paste0("awk '{ print $1\"\t\"$1 }' < ",DIR_PROCESSED,"/hbv_gilead_QC_SAIGE.vcf.gz.sample > ",out_dir,'hla_samples.txt'))
    system(glue::glue("{cur_dir}/bin/bcftools-1.9/bin/bcftools reheader -s {DIR_PROCESSED}/hbv_gilead_QC_SAIGE.vcf.gz.sample {vcf_path} > {gsub('.vcf','.rehead.vcf',vcf_path)}"))
    
    system(glue::glue("/home/zmxu/Software/bgzip -f {gsub('.vcf','.rehead.vcf',vcf_path)}"))
    system(glue::glue("{cur_dir}/bin/bcftools-1.9/bin/bcftools index -f -t {gsub('.vcf','.rehead.vcf.gz',vcf_path)}"))
    
    system(glue::glue("{cur_dir}/bin/bcftools-1.9/bin/bcftools reheader -s {DIR_PROCESSED}/hbv_gilead_QC_SAIGE.vcf.gz.sample {vcf_path_alllele} > {gsub('.vcf','.rehead.vcf',vcf_path_alllele)}"))
    system(glue::glue("/home/zmxu/Software/bgzip -f {gsub('.vcf','.rehead.vcf',vcf_path_alllele)}"))
    system(glue::glue("{cur_dir}/bin/bcftools-1.9/bin/bcftools index -f -t {gsub('.vcf','.rehead.vcf.gz',vcf_path_alllele)}"))
  }
  
  #Generate df with same order as IDs in PLINK file
  covar <- dplyr::left_join(ids_order %>% dplyr::select(ID=X1),pc %>% dplyr::select(-FID),by=c('ID'='ID')) %>%
    dplyr::left_join(pPC,by=c('ID'='ID')) %>%
    dplyr::left_join(dat_covars_raw %>% dplyr::select(ID=host_id,AGE,SEX),by=c('ID'='ID')) %>% as.data.frame(.)
  
  #Remap Sex as binary
  covar$SEX <- ifelse(covar$SEX == 'M',0,1)
  
  #Use existing NULL model file to run new single variant analysis
  tryCatch(expr = {SAIGE::fitNULLGLMM(plinkFile = glue::glue("{DIR_PROCESSED}/hbv_gilead_QC_SAIGE_no_chr6"),
                                      phenoFile = glue::glue('{DIR_PROCESSED}/hbv_gilead_QC_SAIGE_pheno.txt'),
                                      phenoCol = AA_residue,
                                      traitType = 'binary',
                                      invNormalize = F,
                                      covarColList = setdiff(colnames(covar),c("ID")),
                                      sampleIDColinphenoFile = "ID",
                                      LOCO = LOCO,
                                      nThreads = n_cores_input,
                                      outputPrefix = glue::glue('{DIR_SCRATCH}/SAIGE_LOCO_{LOCO}_{AA_residue}'))
  },error = function(e) print(e))
  
  #get all HBV genes (Chr)
  all_chr <- unique(system(glue::glue("{cur_dir}/bin/bcftools-1.9/bin/bcftools query -f '%CHROM\n' {gsub('.vcf','.rehead.vcf.gz',vcf_path)}"),intern = T))
  
  #Run SAIGE for HLA AA
  for(i in 1:length(all_chr)){
    tryCatch(expr = {SAIGE::SPAGMMATtest(vcfFile = glue::glue("{gsub('.vcf','.rehead.vcf.gz',vcf_path)}"),
                                         vcfFileIndex = glue::glue("{gsub('.vcf','.rehead.vcf.gz',vcf_path)}.tbi"),
                                         vcfField = "GT",
                                         chrom = all_chr[i],
                                         minMAC = 4,
                                         IsDropMissingDosages = T,
                                         IsOutputAFinCaseCtrl = T,
                                         IsOutputNinCaseCtrl = T,
                                         minMAF = 0.01,
                                         LOCO = LOCO,
                                         numLinesOutput = .Machine$integer.max,
                                         sampleFile= glue::glue("{DIR_PROCESSED}/hbv_gilead_QC_SAIGE.vcf.gz.sample"),
                                         GMMATmodelFile= glue::glue('{DIR_SCRATCH}/SAIGE_LOCO_{LOCO}_{AA_residue}.rda'),
                                         varianceRatioFile=glue::glue('{DIR_SCRATCH}/SAIGE_LOCO_{LOCO}_{AA_residue}.varianceRatio.txt'),
                                         SAIGEOutputFile= glue::glue("{gsub('.vcf','',vcf_path)}_HBV_{AA_residue}_chr{all_chr[i]}"))},error = function(e) print(e))
  }
  system(glue::glue("head -n 1 {gsub('.vcf','',vcf_path)}_HBV_{AA_residue}_chr{all_chr[1]} > {gsub('.vcf','',vcf_path)}_HBV_{AA_residue}"))
  system(glue::glue("cat {gsub('.vcf','',vcf_path)}_HBV_{AA_residue}_chr* | grep -v 'CHR' >> {gsub('.vcf','',vcf_path)}_HBV_{AA_residue}"))
  system(glue::glue("rm {gsub('.vcf','',vcf_path)}_HBV_{AA_residue}_chr*"))
  saige_result_no_cond <- data.table::fread(glue::glue("{gsub('.vcf','',vcf_path)}_HBV_{AA_residue}"))
  if(nrow(saige_result_no_cond) == 0){
    system(glue::glue("rm {gsub('.vcf','',vcf_path)}_HBV_{AA_residue}"))
  }
  prev_aa_result <- saige_result_no_cond
  
  #Run Saige for HLA Alleles
  for(i in 1:length(all_chr)){
    tryCatch(expr = {SAIGE::SPAGMMATtest(vcfFile = glue::glue("{gsub('.vcf','.rehead.vcf.gz',vcf_path_alllele)}"),
                                         vcfFileIndex = glue::glue("{gsub('.vcf','.rehead.vcf.gz',vcf_path_alllele)}.tbi"),
                                         vcfField = "GT",
                                         chrom = all_chr[i],
                                         minMAC = 4,
                                         IsDropMissingDosages = T,
                                         IsOutputAFinCaseCtrl = T,
                                         IsOutputNinCaseCtrl = T,
                                         minMAF = 0.01,
                                         LOCO = LOCO,
                                         numLinesOutput = .Machine$integer.max,
                                         sampleFile= glue::glue("{DIR_PROCESSED}/hbv_gilead_QC_SAIGE.vcf.gz.sample"),
                                         GMMATmodelFile= glue::glue('{DIR_SCRATCH}/SAIGE_LOCO_{LOCO}_{AA_residue}.rda'),
                                         varianceRatioFile=glue::glue('{DIR_SCRATCH}/SAIGE_LOCO_{LOCO}_{AA_residue}.varianceRatio.txt'),
                                         SAIGEOutputFile= glue::glue("{gsub('.vcf','',vcf_path_alllele)}_HBV_{AA_residue}_chr{all_chr[i]}"))},error = function(e) print(e))
  }
  system(glue::glue("head -n 1 {gsub('.vcf','',vcf_path_alllele)}_HBV_{AA_residue}_chr{all_chr[1]} > {gsub('.vcf','',vcf_path_alllele)}_HBV_{AA_residue}"))
  system(glue::glue("cat {gsub('.vcf','',vcf_path_alllele)}_HBV_{AA_residue}_chr* | grep -v 'CHR' >> {gsub('.vcf','',vcf_path_alllele)}_HBV_{AA_residue}"))
  system(glue::glue("rm {gsub('.vcf','',vcf_path_alllele)}_HBV_{AA_residue}_chr*"))
  saige_result_no_cond <- data.table::fread(glue::glue("{gsub('.vcf','',vcf_path_alllele)}_HBV_{AA_residue}"))
  if(nrow(saige_result_no_cond) == 0){
    system(glue::glue("rm {gsub('.vcf','',vcf_path_alllele)}_HBV_{AA_residue}"))
  }
  prev_allele_result <- saige_result_no_cond
  
  #Run Conditional Analysis
  max_iter <- 5
  cur_iter <- 1
  min_aa_nuid <- c()
  min_aa_id <- c()
  
  min_allele_nuid <- c()
  min_allele_id <- c()
  
  while (cur_iter <= max_iter) {
    #Run AA cond analysis
    continue <- T
    if((nrow(prev_aa_result) != 0) & !is.null(nrow(prev_aa_result))){
      aa_results <- prev_aa_result
      if(cur_iter == 1){
        min_aa <- aa_results[which.min(aa_results$p.value),]
      }else{
        min_aa <- aa_results[which.min(aa_results$p.value_cond),]
        if(length(min_aa_id) == 0){
          continue <- F
        }else{
          if(min_aa$SNPID %in% min_aa_id){
            continue <- F
          }
        }
      }
      if(continue){
        min_aa_nuid <- c(min_aa_nuid,paste0(min_aa$CHR,':',min_aa$POS,"_",min_aa$Allele1,'/',min_aa$Allele2))
        min_aa_id <- c(min_aa_id,min_aa$SNPID)
        min_aa_id_string <- paste0(min_aa_id,collapse = '.')
        prev_min_aa <- min_aa
        false_counter <- 0
        for(i in 1:length(all_chr)){
          tryCatch(expr = {SAIGE::SPAGMMATtest(vcfFile = glue::glue("{gsub('.vcf','.rehead.vcf.gz',vcf_path)}"),
                                               vcfFileIndex = glue::glue("{gsub('.vcf','.rehead.vcf.gz',vcf_path)}.tbi"),
                                               vcfField = "GT",
                                               chrom = all_chr[i],
                                               minMAC = 4,
                                               IsDropMissingDosages = T,
                                               IsOutputAFinCaseCtrl = T,
                                               IsOutputNinCaseCtrl = T,
                                               minMAF = 0.01,
                                               LOCO = LOCO,
                                               numLinesOutput = .Machine$integer.max,
                                               condition = min_aa_nuid,
                                               sampleFile= glue::glue("{DIR_PROCESSED}/hbv_gilead_QC_SAIGE.vcf.gz.sample"),
                                               GMMATmodelFile= glue::glue('{DIR_SCRATCH}/SAIGE_LOCO_{LOCO}_{AA_residue}.rda'),
                                               varianceRatioFile=glue::glue('{DIR_SCRATCH}/SAIGE_LOCO_{LOCO}_{AA_residue}.varianceRatio.txt'),
                                               SAIGEOutputFile= glue::glue("{gsub('.vcf','',vcf_path)}_HBV_{AA_residue}_cond_aa_{min_aa_id_string}_chr{all_chr[i]}"))
          },error = function(e) print(e))
        }
        system(glue::glue("head -n 1 {gsub('.vcf','',vcf_path)}_HBV_{AA_residue}_cond_aa_{min_aa_id_string}_chr{all_chr[1]} > {gsub('.vcf','',vcf_path)}_HBV_{AA_residue}_cond_aa_{min_aa_id_string}"))
        system(glue::glue("cat {gsub('.vcf','',vcf_path)}_HBV_{AA_residue}_cond_aa_{min_aa_id_string}_chr* | grep -v 'CHR' >> {gsub('.vcf','',vcf_path)}_HBV_{AA_residue}_cond_aa_{min_aa_id_string}"))
        system(glue::glue("rm {gsub('.vcf','',vcf_path)}_HBV_{AA_residue}_cond_aa_{min_aa_id_string}_chr*"))
        
        if(!file.exists(glue::glue("{gsub('.vcf','',vcf_path)}_HBV_{AA_residue}_cond_aa_{min_aa_id_string}"))){
          prev_aa_result <- NULL
        }else{
          prev_aa_result <- data.table::fread(glue::glue("{gsub('.vcf','',vcf_path)}_HBV_{AA_residue}_cond_aa_{min_aa_id_string}"))
        }
      }
    }
    
    #Run alelle conditional analysis
    continue <- T
    if((nrow(prev_allele_result) != 0) & !is.null(nrow(prev_allele_result))){
      allele_results <- prev_allele_result
      if(cur_iter == 1){
        min_allele <- allele_results[which.min(allele_results$p.value),]
      }else{
        min_allele <- allele_results[which.min(allele_results$p.value_cond),]
        if(length(min_allele_id) == 0){
          continue <- F
        }else{
          if(min_allele$SNPID %in% min_allele_id){
            continue <- F
          }
        }
      }
      if(continue){
        min_allele_nuid <- c(min_allele_nuid,paste0(min_allele$CHR,':',min_allele$POS,"_",min_allele$Allele1,'/',min_allele$Allele2))
        min_allele_id <- c(min_allele_id,min_allele$SNPID)
        min_allele_id_string <- paste0(min_allele_id,collapse = '.')
        prev_min_allele <- min_allele
        false_counter <- 0
        for(i in 1:length(all_chr)){
          tryCatch(expr = {SAIGE::SPAGMMATtest(vcfFile = glue::glue("{gsub('.vcf','.rehead.vcf.gz',vcf_path_alllele)}"),
                                               vcfFileIndex = glue::glue("{gsub('.vcf','.rehead.vcf.gz',vcf_path_alllele)}.tbi"),
                                               vcfField = "GT",
                                               chrom = all_chr[i],
                                               minMAC = 4,
                                               IsDropMissingDosages = T,
                                               IsOutputAFinCaseCtrl = T,
                                               IsOutputNinCaseCtrl = T,
                                               minMAF = 0.01,
                                               LOCO = LOCO,
                                               numLinesOutput = .Machine$integer.max,
                                               condition = min_allele_nuid,
                                               sampleFile= glue::glue("{DIR_PROCESSED}/hbv_gilead_QC_SAIGE.vcf.gz.sample"),
                                               GMMATmodelFile= glue::glue('{DIR_SCRATCH}/SAIGE_LOCO_{LOCO}_{AA_residue}.rda'),
                                               varianceRatioFile=glue::glue('{DIR_SCRATCH}/SAIGE_LOCO_{LOCO}_{AA_residue}.varianceRatio.txt'),
                                               SAIGEOutputFile= glue::glue("{gsub('.vcf','',vcf_path_alllele)}_HBV_{AA_residue}_cond_allele_{min_allele_id_string}_chr{all_chr[i]}"))
          },error = function(e) print(e))
        }
        system(glue::glue("head -n 1 {gsub('.vcf','',vcf_path_alllele)}_HBV_{AA_residue}_cond_allele_{min_allele_id_string}_chr{all_chr[1]} > {gsub('.vcf','',vcf_path_alllele)}_HBV_{AA_residue}_cond_allele_{min_allele_id_string}"))
        system(glue::glue("cat {gsub('.vcf','',vcf_path_alllele)}_HBV_{AA_residue}_cond_allele_{min_allele_id_string}_chr* | grep -v 'CHR' >> {gsub('.vcf','',vcf_path_alllele)}_HBV_{AA_residue}_cond_allele_{min_allele_id_string}"))
        system(glue::glue("rm {gsub('.vcf','',vcf_path_alllele)}_HBV_{AA_residue}_cond_allele_{min_allele_id_string}_chr*"))
        
        if(!file.exists(glue::glue("{gsub('.vcf','',vcf_path_alllele)}_HBV_{AA_residue}_cond_allele_{min_allele_id_string}"))){
          prev_allele_result <- NULL
        }else{
          prev_allele_result <- data.table::fread(glue::glue("{gsub('.vcf','',vcf_path_alllele)}_HBV_{AA_residue}_cond_allele_{min_allele_id_string}"))
        }
      }
    }
    
    cur_iter <- cur_iter + 1
  }
  
}
RunSAIGEPheno <- function(results_path,cur_outcome,LOCO=T,n_cores_input = 1){
  cur_dir <- here::here()
  load(glue::glue("{results_path}/prep-data.rda"))
  system(glue::glue("mkdir -p {DIR_SCRATCH}/GWAS/"))
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
  dat_alternative_outcome <- dat_alternative_outcome[match(fam_file$V1,dat_alternative_outcome$host_id),]
  #Write Pheno fIle
  pheno_file_path <- glue::glue('{DIR_PROCESSED}/hbv_gilead_QC_SAIGE_pheno_alt.txt')
  data.table::fwrite(dplyr::inner_join(covar,dat_alternative_outcome,by=c('ID' = 'host_id')),file = pheno_file_path,sep = '\t',col.names = T,na = 'NA',quote = F)
  if(length(unique(as.vector(t(dat_alternative_outcome[,cur_outcome])))) == 2){
    mdl_type = 'binary'
  }else{
    mdl_type = 'quantitative'
  }
  nullGLMM_stat <- tryCatch(expr = {
    SAIGE::fitNULLGLMM(plinkFile = plinkFile,
                       phenoFile = pheno_file_path,
                       phenoCol = cur_outcome,
                       traitType = mdl_type,
                       invNormalize = F,
                       covarColList = setdiff(colnames(covar),c("ID")),
                       sampleIDColinphenoFile = "ID",
                       LOCO = LOCO,
                       nThreads = n_cores_input,
                       outputPrefix = glue::glue('{DIR_SCRATCH}/GWAS/SAIGE_LOCO_{LOCO}_{cur_outcome}'))
    T
  },error = function(e){
    write("NULL MODEL ERROR",glue::glue('{DIR_SCRATCH}/GWAS/SAIGE_LOCO_{LOCO}_{cur_outcome}.error'),append = T)
    write(as.character(e),glue::glue('{DIR_SCRATCH}/GWAS/SAIGE_LOCO_{LOCO}_{cur_outcome}.error'),append = T)
    return(F)},
  warning = function(w) {
    write("NULL MODEL WARNING",glue::glue('{DIR_SCRATCH}/GWAS/SAIGE_LOCO_{LOCO}_{cur_outcome}.warning'),append = T)
    write(as.character(w),glue::glue('{DIR_SCRATCH}/GWAS/SAIGE_LOCO_{LOCO}_{cur_outcome}.error'),append = T)
    return(F)})
  if(!nullGLMM_stat){
    return()
  }
  err_counter <- mclapply(1:22,function(j){
    SPA_test_stat <- tryCatch(expr = {
      SAIGE::SPAGMMATtest(vcfFile = glue::glue("{DIR_PROCESSED}/hbv_gilead_QC_SAIGE.vcf.gz"),
                          vcfFileIndex = glue::glue("{DIR_PROCESSED}/hbv_gilead_QC_SAIGE.vcf.gz.tbi"),
                          vcfField = "GT",
                          chrom = as.character(j),
                          minMAC = 4,
                          IsDropMissingDosages = T,
                          IsOutputAFinCaseCtrl = T,
                          IsOutputNinCaseCtrl = T,
                          minMAF = 0.05,
                          LOCO = LOCO,
                          numLinesOutput = .Machine$integer.max,
                          sampleFile= glue::glue("{DIR_PROCESSED}/hbv_gilead_QC_SAIGE.vcf.gz.sample"),
                          GMMATmodelFile= glue::glue('{DIR_SCRATCH}/GWAS/SAIGE_LOCO_{LOCO}_{cur_outcome}.rda'),
                          varianceRatioFile=glue::glue('{DIR_SCRATCH}/GWAS/SAIGE_LOCO_{LOCO}_{cur_outcome}.varianceRatio.txt'),
                          SAIGEOutputFile= glue::glue('{DIR_SCRATCH}/GWAS/SAIGE_SPA_TEST_LOCO_{LOCO}_chr{j}_{cur_outcome}'))
      T
    }
    ,error = function(e){
      write(glue::glue("CHR {j} ERROR"),glue::glue('{DIR_SCRATCH}/GWAS/SAIGE_LOCO_{LOCO}_{cur_outcome}.error'),append = T)
      write(as.character(e),glue::glue('{DIR_SCRATCH}/GWAS/SAIGE_LOCO_{LOCO}_{cur_outcome}.error'),append = T)
      return(F)},
    warning = function(w){
      write(glue::glue("CHR {j} WARNING"),glue::glue('{DIR_SCRATCH}/GWAS/SAIGE_LOCO_{LOCO}_{cur_outcome}.warning'),append = T)
      write(as.character(w),glue::glue('{DIR_SCRATCH}/GWAS/SAIGE_LOCO_{LOCO}_{cur_outcome}.warning'),append = T)
      return(F)})
    return(SPA_test_stat)
  },mc.preschedule = F,mc.cores = n_cores_input,mc.silent = T)
  err_counter <- unlist(err_counter)
  if(sum(err_counter) == 22){
    system(glue::glue("head -n 1 {DIR_SCRATCH}/GWAS/SAIGE_SPA_TEST_LOCO_{LOCO}_chr{which(err_counter)[1]}_{cur_outcome} > {DIR_SCRATCH}/GWAS/SAIGE_SPA_TEST_LOCO_{LOCO}_{cur_outcome}"))
    system(glue::glue("cat {DIR_SCRATCH}/GWAS/SAIGE_SPA_TEST_LOCO_{LOCO}_chr*_{cur_outcome} | grep -v 'CHR' >> {DIR_SCRATCH}/GWAS/SAIGE_SPA_TEST_LOCO_{LOCO}_{cur_outcome}"))
    system(glue::glue("rm {DIR_SCRATCH}/GWAS/SAIGE_SPA_TEST_LOCO_{LOCO}_chr*_{cur_outcome}"))
  }else{
    return()
  }
}