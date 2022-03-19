## This script loads the covariate data, the pathogen data and the plink data
## and does some QC on it. 
# 2019, Aug 8th
# (c) Sina Rüeger

RunPrepData <- function(POP,GT,PCA_type,LOCO,burden,DIR_RUN,n_cores,debugging=F){
  library(glue)
  library(dplyr)
  library(magrittr)
  library(data.table)
  library(stringr)
  ## also need to be installed
  ## phylobase
  ## ape
  ## adephylo
  library(ggplot2)
  library(tidyr)
  library(entropy) ## for shannon entropy
  library(cluster) ## for hierarchical clustering with agnes
  library(furrr) ## future purrr
  library(methods)  ## this needs to be here otherwise phylobase won't work (no clue why)
  source('./scripts/GWAS.utils/R/inv_normal.R')
  
  ##////////////////////////////////////////////////////////////////
  ##                        Dependencies                          //
  ##////////////////////////////////////////////////////////////////
  source("./scripts/G2G/functions-data-pathogen.R")

  load(glue::glue("{DIR_RUN}/setup.rda"))

  ##////////////////////////////////////////////////////////////////
  ##                      Prepare covar data                      //
  ##////////////////////////////////////////////////////////////////
  
  # Read data ---------------------------------------------------------------
  ## USUBJID,STUDYID,SUBJID,Screening ID: these IDs are used internally in clinical trial; can be safely ignored for your analysis
  ## IGM_ID: this ID can be used to identify VCF file for each subject.
  ## gilead_id: this ID can be used to identify viral sequence files for each subject.
  
  dat_covars_raw <-
    #readr::read_delim(FILE_COVARS_IN, delim = ",") %>%
    as_tibble(data.table::fread(FILE_COVARS_IN,sep = ',',header = T)) %>%
    dplyr::rename(host_id = IGM_ID, pathogen_id = gilead_id) %>%
    dplyr::select(-USUBJID,-STUDYID,-SUBJID,-`Screening ID`)
  
  # IDs ---------------------------------------------------------------------
  
  dat_ids <- dat_covars_raw  %>% dplyr::select(host_id, pathogen_id) 
  

  # Covars selection --------------------------------------------------------
  
  dat_covars_discrete <- dat_covars_raw %>% dplyr::select(host_id, SEX)
  dat_covars_numeric <- dat_covars_raw %>% dplyr::select(host_id, AGE)
  ## formatting see details: cnsgenomics.com/software/gcta/GCTA_UserManual_v1.24.pdf

  ## old stuff
  ## OAV_EXPERIENCErefers to prior oral antiviral treatment experience, ie the patient has had a treatment in the past (antiviral, interferon treatment would not be included here).
  ## use same covars as nimisha: "To prevent spurious associations due to host and viral stratification, we included human principal components and viral genotype as covariates"
  #dat_covars_discrete$RACE <- forcats::fct_recode(dat_covars_discrete$RACE, "OTHER"= "BLACK OR AFRICAN AMERICAN", "OTHER"= "NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER", "OTHER" = "WHITE")
  
  # Alternative outcomes --------------------------------------------------------
  # only for src/gGWAS.R

  ## BASELINE_ALT_U/L:           serum alanine aminotransferase levels       
  ## BASELINE_HBSAG_log10_IU/mL: surface antigen
  ## BASELINE_HBVDNA_IU/mL:      hep b dna levels
  dat_alternative_outcome <- dat_covars_raw %>% 
    dplyr::select(host_id, `BASELINE_ALT_U/L`, `BASELINE_HBSAG_log10_IU/mL`, `BASELINE_HBVDNA_IU/mL`, `BASELINE_HBVDNA_Dil_IU/mL`, BASELINE_HBEAG_STATUS) %>% 
    dplyr::mutate(
      viralload_dil = as.numeric(`BASELINE_HBVDNA_Dil_IU/mL`),
      viralload = as.numeric(`BASELINE_HBVDNA_IU/mL`),
      viralload_merged = case_when(
        !is.na(viralload_dil) ~ viralload_dil,
        is.na(viralload_dil) ~ viralload,
        TRUE ~ NA_real_
      ),
      viralload_dil_log10 = log10(viralload_dil),
      viralload_log10 = log10(viralload),
      viralload_merged_log10 = log10(viralload_merged),
      viralload_merged_invnormal = trans_inv_normal(viralload_merged),
      surface_antigen_invnormal = trans_inv_normal(`BASELINE_HBSAG_log10_IU/mL`),
      ALT_level_invnormal = trans_inv_normal(`BASELINE_ALT_U/L`)
    )  %>% 
    dplyr::mutate(antigen_status = case_when(
      BASELINE_HBEAG_STATUS == "Negative" ~ 0, 
      BASELINE_HBEAG_STATUS == "Positive" ~ 1
    )) %>% 
    dplyr::select(
      -BASELINE_HBEAG_STATUS,
      -`BASELINE_HBVDNA_IU/mL`,
      -`BASELINE_HBVDNA_Dil_IU/mL`,
      -`BASELINE_ALT_U/L`,
      -`BASELINE_HBSAG_log10_IU/mL`,
      -viralload,
      -viralload_dil,
      -viralload_merged,
      -viralload_log10,
      -viralload_dil_log10,
      -viralload_merged_log10
    )
  ##////////////////////////////////////////////////////////////////
  ##                       Prepare HBV data                       //
  ##////////////////////////////////////////////////////////////////
  
  
  # Read file --------------------------------------------------------
  
  dat_pathogen_raw <-
    read_data_pathogen(FILE_PATHOGEN_AA_IN, DAT_IDS = dat_ids)
  
  # Remove Individuals with High Missingness in AA Sequences--------------------------------------------------------
  
  dat_pathogen <- apply_qc_pathogen(data = dat_pathogen_raw, 
                                    NA_ids = callrate_ids_thresh, 
                                    NA_variants = 1, 
                                    MAC = 0)
  
  
  # Calculate shannon entropy as an outcome -----------------------------
  dat_pathogen_long <-
    dat_pathogen %>% tidyr::gather("aminoacid", "dosage",-ID)
  entropy_id <-
    dat_pathogen_long %>% group_by(ID) %>% summarise(
      entropy_mm = entropy::entropy(na.omit(dosage), method = "MM"),
      entropy_shrink = entropy::entropy(na.omit(dosage), method = "shrink"),
      entropy_mm_log = log(entropy_mm)
    )
  
  #<!----------Using entropy is not a good idea in a mixed setting, because C and D have hihger frequency of cor promoter and pre-S mutations than gt A and B. -------->
  ## add to alternative outcome
  ## not adding entropy
  ##  dat_alternative_outcome %<>% full_join(entropy_id %>% select(-entropy_mm), by = c("host_id"= "ID"))
  
  
  # Write out IDs -----------------------------------------------------------
  
  idlist <- data.frame(id = dat_pathogen$ID, id2 = dat_pathogen$ID)
  data.table::fwrite(idlist,file = glue::glue("{DIR_PROCESSED}/tmp_idlist.txt"),col.names = F,sep = ' ')
  # write_delim(
  #   idlist,
  #   path = glue::glue("{DIR_PROCESSED}/tmp_idlist.txt"),
  #   col_names = FALSE,
  #   delim = " "
  # )
  
  
  ##///////////////////////////////////////////////////////////////
  ##                      Prepare host data                      //
  ##///////////////////////////////////////////////////////////////
  
  ## use either a vector of POP or HBV_GT to select individuals. 
  ## if you want to use ALL OF THEM, use all levels, e.g. POP <- c("asian", "european", "remaining"). 
  
  ## VCF files have been prepared using a script in data/
  ## this file turns all VCF files into a plink file
  
  ## select clusters (POP) ------------------------------------------
  #Make plink file keeping IDs which are avaliable
  system(
    glue::glue(
      "{PLINK1} --bfile {FILE_HOST_IN} --geno {callrate_snps_thresh} --maf 0.01 --hwe {hwe_thresh} --make-bed --out {FILE_HOST_OUT}_cluster --keep {DIR_PROCESSED}/tmp_idlist.txt"
    )
  )
  
  #Make GRM
  system(
    glue::glue(
      "{GCTA} --bfile {FILE_HOST_OUT}_cluster --autosome --make-grm --out {FILE_GRM_OUT}_cluster --thread-num {n_cores}"
    )
  )
  system(glue::glue("{GCTA} --grm {FILE_GRM_OUT}_cluster  --pca 8 --out {FILE_GRM_OUT}_cluster"))
  pc <- as_tibble(data.table::fread(glue::glue("{FILE_GRM_OUT}_cluster.eigenvec"), header = F)) %>% dplyr::rename(ID = V1, FID = V2)
  #pc <- readr::read_table2(glue::glue("{FILE_GRM_OUT}_cluster.eigenvec"), col_names = FALSE) %>% dplyr::rename(ID = X1, FID = X2)
  cl <- cluster::agnes((pc %>% dplyr::select(-ID, -FID)), method = "ward", diss = FALSE, metric = "euclidean")
  gr_cluster <- cutree(cl, k = 3)
  
  ## choose clusters -------------------------
  
  tmp_cluster <-
    dat_covars_raw %>% dplyr::right_join(pc %>% dplyr::select(-FID) %>% cbind(gr_cluster), by = c("host_id" = "ID"))
  dat_covars_raw$GT[dat_covars_raw$GT == 'Mixed genotype detected'] <- 'Mixed' 
  ## cluster 1: ASIAN
  ## cluster 2: EUROPEAN
  ## cluster 3: AFRICAN
  
  tmp_cluster$ancestry <- dplyr::case_when(
    tmp_cluster$gr_cluster == 1 ~ "asian",
    tmp_cluster$gr_cluster == 2 ~ "european",
    TRUE ~ "remaining"
  )
  
  #Keep intersect between genetically obtained ancestry and self reported ancestry
  tmp_cluster$SELF_REPORTED_ANCESTRY <- ifelse(tmp_cluster$RACE == 'ASIAN','asian',ifelse(tmp_cluster$RACE == 'WHITE','european','other'))
  dat_ids_pop <- tmp_cluster %>% dplyr::filter(tolower(ancestry) %in% tolower(POP) & tolower(SELF_REPORTED_ANCESTRY) %in% tolower(POP))
  
  #Filter based on HBV Genotype
  dat_ids_pop <- dat_ids_pop %>% dplyr::select(GT, host_id, pathogen_id) %>% filter(GT %in% HBV_GT)
  
  idlist_pop <-
    data.frame(id = dat_ids_pop$host_id, id2 = dat_ids_pop$host_id)
  data.table::fwrite(idlist_pop,file = glue::glue("{DIR_PROCESSED}/tmp_idlist_pop.txt"),col.names = F,sep = ' ')
  # readr::write_delim(
  #   idlist_pop,
  #   path = glue::glue("{DIR_PROCESSED}/tmp_idlist_pop.txt"),
  #   col_names = FALSE,
  #   delim = " "
  # )
  
  fs::file_copy(glue::glue("{DIR_PROCESSED}/tmp_idlist_pop.txt"), glue::glue("{DIR_SCRATCH}/idlist_pop.txt"), overwrite = TRUE)
  
  idlist_pop_path <- as.character(dat_ids_pop$pathogen_id)
  
  # revisit pathogen data, remove individuals ------------------------------------
  dat_pathogen <- dat_pathogen %>% filter(ID %in% idlist_pop$id)
  
  dat_pathogen <- apply_qc_pathogen(data = dat_pathogen, 
                                    NA_ids = callrate_path_ids_thresh, 
                                    NA_variants = callrate_path_variant_thresh, 
                                    MAC = mac_path_variant_thresh)
  
  # PLINK ------------------------------------------------------------------------
  # debugging just selects only chr 22
  
  if(debugging)
  {
    
    FILE_HOST_tmp0 <- glue::glue("{FILE_HOST_tmp0}_debugging")
    FILE_HOST_tmp1 <- glue::glue("{FILE_HOST_tmp1}_debugging")
    FILE_HOST_OUT <- glue::glue("{FILE_HOST_OUT}_debugging")
    
    ## apply callrate snps
    ## maf
    ## hwe
    ## rename missing snps
    system(
      glue::glue(
        "{PLINK1} --bfile {FILE_HOST_IN} --geno {callrate_snps_thresh} --maf {maf_thresh} --hwe {hwe_thresh} --make-bed --out {FILE_HOST_tmp0} --keep {DIR_PROCESSED}/tmp_idlist_pop.txt --chr 20-22"
      )
    )
    
    
  }else{
    
    
    ## do clustering
    
    ## apply callrate snps
    ## maf
    ## hwe
    ## rename missing snps: set SNP == . to chr:pos
    ##  --set-missing-var-ids @:#:\\$1:\\$2
    if(burden){
      system(
        glue::glue(
          "{PLINK1} --bfile {FILE_HOST_IN} --geno {callrate_snps_thresh} --maf {maf_thresh} --max-maf {maf_max} --hwe {hwe_thresh} --make-bed --out {FILE_HOST_tmp0} --autosome --keep {DIR_PROCESSED}/tmp_idlist_pop.txt" 
        )
      )
      
    }else{
      system(
        glue::glue(
          "{PLINK1} --bfile {FILE_HOST_IN} --geno {callrate_snps_thresh} --maf {maf_thresh} --hwe {hwe_thresh} --make-bed --out {FILE_HOST_tmp0} --autosome --keep {DIR_PROCESSED}/tmp_idlist_pop.txt" 
        )
      )
    }
  }
  
  
  ## apply callrate individuals
  ## maf, hwe
  system(
    glue::glue(
      "{PLINK1} --bfile {FILE_HOST_tmp0} --mind {callrate_ids_thresh} --geno {callrate_snps_thresh} --maf {maf_thresh} --hwe {hwe_thresh} --make-bed --out {FILE_HOST_tmp1}"
    )
  )
  
  
  ## KING -------------------------------------------
  system(
   glue::glue(
     "{KING} -b {FILE_HOST_IN}.bed --related --duplicate" # --reml-bendV
   )
  )
  ## no duplicate
  
  # ## KING -------------------------------------------
  system(
    glue::glue(
      "{KING} -b {FILE_HOST_tmp1}.bed --unrelated" # --reml-bendV
    )
  )
  #KING1     igm160324,igm160612
  #KING2     igm160325,igm160613
  #KING3     igm160030,igm160031
  #KING4     igm160870,igm160883
  
  # system(glue::glue("mv kingunrelated.txt {DIR_PROCESSED}/kingunrelated.txt"))
  fs::file_move("kingunrelated.txt", glue::glue("{DIR_PROCESSED}/kingunrelated.txt"))
  fs::file_delete("kingallsegs.txt")
  fs::file_delete("kingunrelated_toberemoved.txt")

  # apply king
  system(
    glue::glue(
      "{PLINK1} --bfile {FILE_HOST_tmp1} --make-bed --out {FILE_HOST_OUT} --keep {DIR_PROCESSED}/kingunrelated.txt"
    )
  )
  
  ids_unrelated <- as_tibble(data.table::fread(glue::glue("{DIR_PROCESSED}/kingunrelated.txt"),header = F)) %>% dplyr::rename(X1=V1,X2=V2) #read_tsv(glue::glue("{DIR_PROCESSED}/kingunrelated.txt"), col_names = FALSE)
  # data.table::fwrite(ids_unrelated,glue::glue("{DIR_PROCESSED}/kingunrelated.txt"),sep = ' ',col.names = F,row.names = F,quote = F)
  
  dat_ids %<>% filter(host_id %in% ids_unrelated$X1)
  
  ## calc frequencies
  system(
    glue::glue(
      "{PLINK1} --bfile {FILE_HOST_OUT} --out {DIR_SCRATCH}/hbv_gilead_QC --freq" 
    )
  )
  
  # Read fam ----------------------------------------------------------------
  
  ids_order <-
    as_tibble(data.table::fread(glue::glue("{FILE_HOST_OUT}.fam"),header = F) %>% dplyr::select(X1=V1,X2=V2,X3=V3,X4=V4,X5=V5,X6=V6))
    # readr::read_delim(glue::glue("{FILE_HOST_OUT}.fam"),
    #            delim = " ",
    #            col_names = FALSE)
    # 
  
  ##///////////////////////////////////////////////////////////////
  ##                          Estimate GRM                       //
  ##///////////////////////////////////////////////////////////////
  ## note: this is not needed when mlma-loco
  ## only to generate PCs
  
  system(
    glue::glue(
      "{PLINK1} --bfile {FILE_HOST_OUT} --not-chr 6 --make-bed --out {FILE_HOST_OUT}_grm"
    )
  )
  
  
  
  system(
    glue::glue(
      "{GCTA} --bfile {FILE_HOST_OUT}_grm --autosome --maf 0.05 --make-grm --out {FILE_GRM_OUT} --thread-num {n_cores}"
    )
  )
  ## suggested here: http://gcta.freeforums.net/thread/243/variance-covaraince-matrix-positive-definite
  
  
  
  ##///////////////////////////////////////////////////////////////
  ##                Add host PCs                                 //
  ##///////////////////////////////////////////////////////////////
  
  system(glue::glue("{GCTA} --grm {FILE_GRM_OUT} --pca {n_pc} --out {FILE_GRM_OUT}"))#  --keep {DIR_PROCESSED}/tmp_idlist_pop.txt"))
  
  ## extract pcs
  pc <-
    as_tibble(data.table::fread(glue::glue("{FILE_GRM_OUT}.eigenvec"),header = F))
    # readr::read_delim(glue::glue("{FILE_GRM_OUT}.eigenvec"),
    #            delim = " ",
    #            col_names = FALSE)
  names(pc) <- c("ID", "FID", paste0("host_PC", 1:n_pc))
  dat_covars_numeric <-
    dat_covars_numeric %>% left_join(pc %>% dplyr::select(-FID), by = c("host_id" = "ID"))
  
  
  ## loadings ------------------
  if(FALSE)
  {
    system(
      glue::glue(
        "{GCTA} --bfile {FILE_HOST_OUT} --pc-loading {FILE_GRM_OUT} --out {FILE_GRM_OUT}"
      )
    )
    
    ## extract pcs
    pcl <-
      readr::read_delim(glue::glue("{FILE_GRM_OUT}.pcl"),
                 delim = "\t",
                 col_names = TRUE)
    
    #which SNPs contribute a lot?
    thresh <- 0.001
    pcl$driver1 <-
      (pcl$pc1_loading < quantile(pcl$pc1_loading, c(thresh))) |
      (pcl$pc1_loading > quantile(pcl$pc1_loading, c(1 - thresh)))
    pcl$driver2 <-
      (pcl$pc2_loading < quantile(pcl$pc2_loading, c(thresh))) |
      (pcl$pc2_loading > quantile(pcl$pc2_loading, c(1 - thresh)))
    pcl$driver3 <-
      (pcl$pc3_loading < quantile(pcl$pc3_loading, c(thresh))) |
      (pcl$pc3_loading > quantile(pcl$pc3_loading, c(1 - thresh)))
    pcl$driver4 <-
      (pcl$pc4_loading < quantile(pcl$pc4_loading, c(thresh))) |
      (pcl$pc4_loading > quantile(pcl$pc4_loading, c(1 - thresh)))
    
    map <-
      readr::read_delim(glue::glue("{FILE_HOST_OUT}_grm.bim"),
                 delim = "\t",
                 col_names = FALSE)
    
    
    pcl_freq <-
      pcl %>% filter(driver4) %>% left_join(map, by = c("SNP" = "X2"))
    
    table(pcl_freq$X1)
    
  }
  
  
  ##///////////////////////////////////////////////////////////////
  ##                Pathogen principal components                //
  ##///////////////////////////////////////////////////////////////
  if (debugging) {
    ## use ordinary pca when debugging
    pPC <- pathogen_pca(
      path_tree = FILE_PATHOGEN_TREE,
      path_pathogen = FILE_PATHOGEN_AA_IN,
      id.list = idlist_pop_path,
      pca = TRUE,
      ppca = FALSE,
      n.pc = n_ppc,
      loadings = FALSE
    )
    
    
  } else{
    if(PCA_type == 'pPCA'){
      pPC <- pathogen_pca(
        path_tree = FILE_PATHOGEN_TREE,
        path_pathogen = FILE_PATHOGEN_AA_IN,
        id.list = idlist_pop_path,
        pca = FALSE,
        ppca = TRUE,
        n.pc = n_ppc,
        loadings = FALSE
      )
    }else if(PCA_type == 'PCA'){
      pPC <- pathogen_pca(
        path_tree = FILE_PATHOGEN_TREE,
        path_pathogen = FILE_PATHOGEN_AA_IN,
        id.list = idlist_pop_path,
        pca = TRUE,
        ppca = FALSE,
        n.pc = n_ppc,
        loadings = FALSE)
    }
    
    
  }
  
  ## add the correct IDs from host
  pPC <- dat_ids %>%
    left_join(pPC %>% mutate_if(is.factor, as.character),
              c("pathogen_id" = "ID")) %>%
    dplyr::select(-pathogen_id) %>%
    dplyr::rename(ID = host_id)
  
  dat_covars_numeric <-
    dat_covars_numeric %>% left_join(pPC, by = c("host_id" = "ID"))
  save(list = ls(all.names = TRUE), file =  glue("{DIR_SCRATCH}/prep-data.rda", envir = 
                                                   environment()))
  
}
args <-  commandArgs(trailingOnly = T)
POP <- strsplit(args[[1]],split = '_')[[1]]
PCA_type <- args[[2]]
GT <- strsplit(args[[3]],split = '_')[[1]]
LOCO <- as.logical(args[[4]])
burden <- as.logical(args[[5]])
DIR_RUN <- args[[6]]
n_cores <- as.numeric(args[[7]])

RunPrepData(POP,GT,PCA_type,LOCO,burden,DIR_RUN,n_cores)