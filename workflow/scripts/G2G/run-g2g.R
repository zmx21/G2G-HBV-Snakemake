# 2019, Aug 8th
# (c) Sina Rüeger

## Aim: This script runs the G2G analysis for HBV
## It contains all the steps needed from data preparation onwards
## NOTE: Some pre-datapreparation (vcf to plink conversion), was done in data/combine_vcfs.sh

RunG2G <- function(POP,GT,PCA_type,LOCO,DIR_SCRATCH,n_cores,debugging = F,sigHitsOnly=F,burden=F,run_LoF = T,run_missense = F){
  ##////////////////////////////////////////////////////////////////
  ##                        Dependencies                          //
  ##////////////////////////////////////////////////////////////////
  library(dplyr)
  library(pbmcapply)
  library(data.table)
  library(tidyr)
  library(here)
  
  source("./scripts/G2G/functions-data-pathogen.R")
  source("./scripts/G2G/run-SAIGE.R")
  source("./scripts/G2G/run-SAIGE-SKAT.R")

  #Load Setup stored Variables
  load(glue::glue("{DIR_SCRATCH}/prep-data.rda"))
  
  ## add string to debugging
  if (debugging) {
    NAM <- glue::glue("{NAM}_debugging")
  }
  
  ##////////////////////////////////////////////////////////////////
  ##                    Load data                                 //
  ##////////////////////////////////////////////////////////////////
  ## Assemble plan --------------------------------------------------------------
  sm <- summarise_pathogen(data = dat_pathogen)
  plan <-
    data.frame(
      global.outcome = NAM,
      outcome = sm$aa,
      outcome.freq = sm$maf,
      outcome.count = sm$mac,
      covars = NA,
      str.out = NA,
      include = NA
    )
  
  if (debugging)
  {
    ## if debugging, just select first three
    plan <- plan %>% filter(between(outcome.freq, 0.2, 0.8)) %>% slice(1:3)
  }
  
  plan <- plan %>% dplyr::rename(outcome.freq.all = outcome.freq, outcome.count.all = outcome.count) %>%
    left_join(sm %>% dplyr::rename(outcome.freq = maf, outcome.count = mac), by = c("outcome" = "aa"))
  
  ## write out plan
  data.table::fwrite(plan,file = glue::glue("{DIR_SCRATCH}/plan_{NAM}.txt"),sep = ' ')
  # readr::write_delim(plan,
  #             path = glue::glue("{DIR_SCRATCH}/plan_{NAM}.txt"),
  #             delim = " ")
  
  ##////////////////////////////////////////////////////////////////
  ##                      Decide which Outcomes                   //
  ##////////////////////////////////////////////////////////////////
  
  
  
  if(debugging)
  {
    OUTCOME <- plan %>% pull(outcome)
    plan <- plan %>% mutate(include.gcta = outcome %in% OUTCOME)
    
  }else{
    OUTCOME <- plan %>% filter(outcome.count >= mac_path_variant_thresh) %>% pull(outcome)
    plan <- plan %>% mutate(include.gcta = outcome %in% OUTCOME)
    
    ## TODO: why? 
    OUTCOME <- intersect(names(dat_pathogen[,-1]), OUTCOME)
    
  }
  
  n_outcome <- length(OUTCOME)
  
  dat_pathogen <- dat_pathogen[,c("ID", OUTCOME)]
  
  
  
  ##///////////////////////////////////////////////////////////////
  ##                    Harmonise Identifiers                    //
  ##///////////////////////////////////////////////////////////////
  ## join with fam file
  
  # Write out pathogen data -----------------------------------
  
  dat_pathogen <-
    left_join(ids_order %>% select(X1, X2),
              dat_pathogen,
              by = c("X1" = "ID")) 
  ## :::
  ## there should be an overlap of 725 individuals!, see ../munge/investigate-ids.R
  ## :::
  
  stopifnot(identical(as.character(ids_order$X1), as.character(dat_pathogen$X1)))
  data.table::fwrite(dat_pathogen,file = FILE_PATHOGEN_OUT,sep = ' ',col.names = F)
  
  # readr::write_delim(dat_pathogen,
  #             path = FILE_PATHOGEN_OUT,
  #             col_names = FALSE,
  #             delim = " ")
  
  # Write out covars data -----------------------------------
  
  ## turn into same order than genotyped data
  dat_covars_numeric <-
    left_join(ids_order %>% select(X1, X2),
              dat_covars_numeric,
              by = c("X1" = "host_id")) 
  
  
  
  dat_covars_discrete <-
    left_join(ids_order %>% select(X1, X2),
              dat_covars_discrete,
              by = c("X1" = "host_id")) 
  
  
  ## write out
  stopifnot(identical(as.character(ids_order$X1), as.character(dat_covars_discrete$X1)))
  stopifnot(identical(as.character(ids_order$X1), as.character(dat_covars_numeric$X1)))
  data.table::fwrite(dat_covars_discrete,file = FILE_COVARS_DISCRETE_OUT,sep = ' ',col.names = F)
  
  # readr::write_delim(dat_covars_discrete,
  #             FILE_COVARS_DISCRETE_OUT,
  #             delim = " ",
  #             col_names = FALSE)
  
  dat_covars_numeric_sub1 <- dat_covars_numeric %>% select(X1, X2, AGE, pPC1, pPC2, pPC3, pPC4, pPC5, pPC6)
  data.table::fwrite(dat_covars_numeric_sub1,file = FILE_COVARS_NUMERIC_OUT,sep = ' ',col.names = F)
  
  # readr::write_delim(dat_covars_numeric_sub1,
  #             FILE_COVARS_NUMERIC_OUT,
  #             delim = " ",
  #             col_names = FALSE)
  
  
  
  ## plink covars
  dat_covars_plink <- dat_covars_discrete %>% left_join(dat_covars_numeric) %>% mutate(SEX_num = case_when(SEX == "F" ~ 1, 
                                                                                                           SEX == "M" ~ 0, 
                                                                                                           TRUE ~ NA_real_)) %>%
    select(-SEX)
  FILE_COVARS_PLINK <- glue::glue("{DIR_PROCESSED}/covars-plink.covar")
  stopifnot(identical(as.character(ids_order$X1), as.character(dat_covars_plink$X1)))
  data.table::fwrite(dat_covars_plink,file = FILE_COVARS_PLINK,sep = ' ',col.names = F)
  
  # readr::write_delim(dat_covars_plink,
  #             FILE_COVARS_PLINK,
  #             delim = " ",
  #             col_names = FALSE)
  
  
  ##////////////////////////////////////////////////////////////////
  ##                           Run PLINK                          //
  ##////////////////////////////////////////////////////////////////
  OUTCOME = colnames(dat_pathogen[,c(-1,-2)])
  
  if(sigHitsOnly){
    OUTCOME_sig <- c()
    if('Asian' %in% POP){
      OUTCOME_sig <-c(OUTCOME_sig,'gene_S_pos_0017_A') #c(OUTCOME_sig,"gene_PC_C_pos_0160_A",'gene_PC_C_pos_0160_P','gene_Pol_pos_0049_N','gene_Pol_pos_0197_C','gene_Pol_pos_0215_Q','gene_S_pos_0017_A','gene_S_pos_0032_L','gene_S_pos_0035_R','gene_S_pos_0051_P')
    }
    if('European' %in% POP){
      OUTCOME_sig <- c(OUTCOME_sig,'gene_PC_C_pos_0067_Y')
    }
  }else{
    OUTCOME_sig <- OUTCOME
  }

  OUTCOME_sig_which <- which(OUTCOME %in% OUTCOME_sig)
  if(!burden){
    for (counter_within_outcome in OUTCOME_sig_which) {
      
      str_out_plink <-
        glue::glue("{NAM_PLINK}_{OUTCOME[counter_within_outcome]}")
      
      if (debugging) {
        ## DEBUGGING ------------------------------------------
        ## no covars
        
        system(
          glue::glue(
            "{PLINK} --bfile {FILE_HOST_OUT} --threads {n_cores} --no-sex --logistic --pheno {FILE_PATHOGEN_OUT} --pheno-col-nums {counter_within_outcome + 2}  --out {DIR_SCRATCH}/{str_out_plink} --1 "
          )
        )
        
      } else {
        
        
        system(
          glue::glue(
            "{PLINK} --bfile {FILE_HOST_OUT} --no-sex --threads {n_cores} --glm firth --pheno {FILE_PATHOGEN_OUT} --pheno-col-nums {counter_within_outcome + 2}  --out {DIR_SCRATCH}/{str_out_plink} --covar {FILE_COVARS_PLINK} --1 --ci 0.95"
          )
        )
        
        
        system(
          glue::glue(
            "{PLINK1} --bfile {FILE_HOST_OUT} --no-sex --threads {n_cores} --freq case-control --pheno {FILE_PATHOGEN_OUT} --mpheno {counter_within_outcome}  --out {DIR_SCRATCH}/{str_out_plink} --1 "
          )
        )
        
        
      }
    }
    
  }
  
  ##////////////////////////////////////////////////////////////////
  ##                           Run SAIGE                           //
  ##////////////////////////////////////////////////////////////////
  
  existing_sage_files <- dir(DIR_SCRATCH)
  if(!burden){
    existing_sage_files <- existing_sage_files[grepl(glue::glue("SAIGE_SPA_TEST_LOCO_{as.character(LOCO)}"),existing_sage_files)]
    existing_sage_pos <- sapply(existing_sage_files,function(x) strsplit(x,split = glue::glue('SAIGE_SPA_TEST_LOCO_{as.character(LOCO)}_'))[[1]][2])
  }else{
    if(run_LoF){
      existing_sage_files <- existing_sage_files[grepl(glue::glue("SAIGE_BURDEN_LOF_LOCO_{as.character(LOCO)}"),existing_sage_files)]
      existing_sage_files <- existing_sage_files[!grepl(glue::glue(".txt"),existing_sage_files)]
      existing_sage_pos <- sapply(existing_sage_files,function(x) strsplit(x,split = glue::glue('SAIGE_BURDEN_LOF_LOCO_{as.character(LOCO)}_'))[[1]][2])
    }else if(run_missense){
      existing_sage_files <- existing_sage_files[grepl(glue::glue("SAIGE_BURDEN_MISSENSE_LOCO_{as.character(LOCO)}"),existing_sage_files)]
      existing_sage_files <- existing_sage_files[!grepl(glue::glue(".txt"),existing_sage_files)]
      existing_sage_pos <- sapply(existing_sage_files,function(x) strsplit(x,split = glue::glue('SAIGE_BURDEN_MISSENSE_LOCO_{as.character(LOCO)}_'))[[1]][2])
    }
  }

  sum(!sapply(setdiff(OUTCOME,existing_sage_pos),function(x) any(grepl(pattern = x,x = dir(DIR_SCRATCH)))))
  if(!burden){
    CreateGRMGenoFile(DIR_SCRATCH)
  }else{
    if(!('missense_variants' %in% dir(DIR_PROCESSED) & 'LOF_variants' %in% dir(DIR_PROCESSED))){
      CreateGRMGenoFileBURDEN(DIR_SCRATCH,n_cores)
    }
  }
  print(c(POP,GT,PCA_type,LOCO,n_cores,debugging))
  if(!sigHitsOnly){
    if(!burden){
      log <- lapply(which(OUTCOME %in% setdiff(OUTCOME,c(existing_sage_pos,OUTCOME_sig))), function(counter_within_outcome){
        cur_outcome <- OUTCOME[counter_within_outcome]
        RunSAIGE(DIR_SCRATCH,AA_residue=cur_outcome,LOCO=LOCO,n_cores_input = n_cores)
      })
    }else{
      log <- lapply(which(OUTCOME %in% setdiff(OUTCOME,c(existing_sage_pos,OUTCOME_sig))), function(counter_within_outcome){
        cur_outcome <- OUTCOME[counter_within_outcome]
        RunSAIGEBurden(DIR_SCRATCH,AA_residue=cur_outcome,LOCO=LOCO,n_cores_input = n_cores,run_LoF = run_LoF,run_missense = run_missense)
      })
    }

    for (counter_within_outcome in which(OUTCOME %in% setdiff(OUTCOME,c(existing_sage_pos,OUTCOME_sig)))){
      cur_outcome <- OUTCOME[counter_within_outcome]
      RunSAIGE(DIR_SCRATCH,AA_residue=cur_outcome,LOCO=LOCO,n_cores = n_cores)
    }
  }else{
    if(!burden){
      log <- lapply(which(OUTCOME %in% setdiff(OUTCOME_sig,existing_sage_pos)), function(counter_within_outcome){
        cur_outcome <- OUTCOME[counter_within_outcome]
        RunSAIGE(DIR_SCRATCH,AA_residue=cur_outcome,LOCO=LOCO,n_cores_input = n_cores)
      })
    }else{
      log <- lapply(which(OUTCOME %in% setdiff(OUTCOME_sig,existing_sage_pos)), function(counter_within_outcome){
        cur_outcome <- OUTCOME[counter_within_outcome]
        RunSAIGEBurden(DIR_SCRATCH,AA_residue=cur_outcome,LOCO=LOCO,n_cores_input = n_cores,run_LoF = run_LoF,run_missense = run_missense)
      })
    }

    for (counter_within_outcome in which(OUTCOME %in% setdiff(OUTCOME_sig,existing_sage_pos))) {
      cur_outcome <- OUTCOME[counter_within_outcome]
      RunSAIGE(DIR_SCRATCH,AA_residue=cur_outcome,LOCO=LOCO,n_cores = n_cores)
    }
  }
  
  # ///////////////
  # write out plan
  # ///////////////
  data.table::fwrite(plan,file = glue::glue("{DIR_SCRATCH}/plan_{NAM}.txt"),sep = ' ',col.names = F)
  
  # readr::write_delim(plan,
  #             path = glue::glue("{DIR_SCRATCH}/plan_{NAM}.txt"),
  #             delim = " ")
  saveRDS(log,file = glue::glue("{DIR_SCRATCH}/g2g-log.rds"))
}

args <- commandArgs(trailingOnly = TRUE)
POP <- strsplit(args[[1]],split = '_')[[1]]
PCA_type <- args[[2]]
GT <- strsplit(args[[3]],split = '_')[[1]]
LOCO <- as.logical(args[[4]])
burden <- as.logical(args[[5]])
sigHitsOnly <- as.logical(args[[6]])
DIR_RUN <- args[[7]]
n_cores <- as.numeric(args[[8]])

RunG2G(POP=POP,GT=GT,PCA_type=PCA_type,LOCO=LOCO,burden=burden,DIR_SCRATCH = DIR_RUN,n_cores=n_cores,sigHitsOnly=sigHitsOnly)



# if(length(args) > 0){
#   RunG2G(POP = unlist(strsplit(args[1],split = ',')),GT =unlist(strsplit(args[2],split = ',')),PCA_type = args[3],LOCO = as.logical(args[4]),n_cores = as.numeric(args[5]),debugging = ifelse(args[6]=='forreal',F,T),sigHitsOnly = as.logical(args[7]),burden = as.logical(args[8]))
# }
# RunG2G(POP = c('asian'),GT = c('B'),PCA_type = 'pPCA',LOCO = F,n_cores = 22,debugging = F,sigHitsOnly = T,burden = F)
# RunG2G(POP = c('asian'),GT = c('C'),PCA_type = 'pPCA',LOCO = F,n_cores = 22,debugging = F,sigHitsOnly = T,burden = F)
# RunG2G(POP = c('asian'),GT = c('A','C','D','B'),PCA_type = 'pPCA',LOCO = F,n_cores = 22,debugging = F,sigHitsOnly = T,burden = F)

# RunG2G(POP = c('european'),GT = c('A','C','D','F','H'),PCA_type = 'pPCA',LOCO = F,n_cores = 22,debugging = F,sigHitsOnly = F,burden = F)
# RunG2G(POP = c('european'),GT = c('D'),PCA_type = 'pPCA',LOCO = F,n_cores = 22,debugging = F,sigHitsOnly = T,burden = F)

# RunG2G(POP = c('asian'),GT = c('A','C','D','B'),PCA_type = 'pPCA',LOCO = F,n_cores = 22,debugging = F,sigHitsOnly = F,burden = F)
# RunG2G(POP = c('european'),GT = c('A','C','D','F','H'),PCA_type = 'pPCA',LOCO = F,n_cores = 22,debugging = F,sigHitsOnly = F,burden = F)

# RunG2G(POP = c('asian'),GT = c('A','C','D','B'),PCA_type = 'pPCA',LOCO = F,n_cores = 22,debugging = F,sigHitsOnly = F,burden = T)
# RunG2G(POP = c('european'),GT = c('A','C','D','F','H'),PCA_type = 'pPCA',LOCO = F,n_cores = 22,debugging = F,sigHitsOnly = F,burden = T)

# RunG2G(POP = c('asian','european'),GT = c('A','C','D','B','F','H'),PCA_type = 'pPCA',LOCO = F,n_cores = 22,debugging = F,sigHitsOnly = F,burden = F)

# RunG2G(POP = c('asian'),GT = c('A','C','D','B'),PCA_type = 'pPCA',LOCO = F,n_cores = 100,debugging = F,sigHitsOnly = F,burden = T,run_LoF = F,run_missense = T)
# RunG2G(POP = c('european'),GT = c('A','C','D','F','H'),PCA_type = 'pPCA',LOCO = F,n_cores = 15,debugging = F,sigHitsOnly = F,burden = T,run_LoF = T,run_missense = F)
