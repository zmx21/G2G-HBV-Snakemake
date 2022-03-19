# 2019, Aug 8th
# (c) Sina Rüeger

## Aim: define all "global" variables and loads all necessary packages

##////////////////////////////////////////////////////////////////
##                      load libraries                          //
##////////////////////////////////////////////////////////////////
RunSetUp <- function(POP,GT,PCA_type,LOCO,burden,n_pc,n_ppc,DIR_RAW,n_cores){
  # library(remotes)
  # remotes::install_github("sinarueeger/ggGWAS")
  # remotes::install_github("sinarueeger/GWAS.utils")
  # 
  # library(dplyr)
  # library(magrittr)
  # library(here)
  # library(readr)
  # library(stringr)
  # ## also need to be installed
  # ## phylobase
  # ## ape
  # ## adephylo
  # library(ggplot2)
  # library(ggGWAS) ## remotes::install_github("sinarueeger/ggGWAS")
  # library(GWAS.utils) ## remotes::install_github("sinarueeger/GWAS.utils")
  # library(tidyr)
  # library(entropy) ## for shannon entropy
  # library(cluster) ## for hierarchical clustering with agnes
  # library(furrr) ## future purrr
  # library(methods)  ## this needs to be here otherwise phylobase won't work (no clue why)
  # 
  # ## from https://github.com/clauswilke/colorblindr/blob/master/R/palettes.R
  # palette_OkabeIto_black <- c("#000000", "#009E73", "#0072B2", "#D55E00", "#56B4E9", "#CC79A7", "#000000", "#F0E442")
  # library(qqman)
  library(glue)
  library(fs)
  ##////////////////////////////////////////////////////////////////
  ##                      define paths                            //
  ##////////////////////////////////////////////////////////////////
  
  DIR_HOST <- paste0(DIR_RAW,'wes_plink') #here::here("data", "raw", "gilead_20181126", "wes_plink")
  DIR_PATHOGEN <- paste0(DIR_RAW,'viral_seq') #here::here("data", "raw", "gilead_20181126", "viral_seq")
  DIR_COVAR <- DIR_RAW
  if(burden){
    SUFFIX <- glue("BURDEN_POP_{paste(POP,collapse = '_')}_GT_{paste(GT,collapse = '_')}_TYPE_{PCA_type}_LOCO_{LOCO}")
  }else{
    SUFFIX <- glue("POP_{paste(POP,collapse = '_')}_GT_{paste(GT,collapse = '_')}_TYPE_{PCA_type}_LOCO_{LOCO}")
  }

  ## PROCESSED: this is for all processed raw data, can be deleted
  DIR_PROCESSED <- glue("../results/processed_{SUFFIX}") #here::here("data", glue("processed_{SUFFIX}"))
  system(glue('mkdir -p {DIR_PROCESSED}'))
  
  ## OUTPUT: these are all shiny figures and summarised results
  DIR_OUTPUT <- glue("../results/output_{SUFFIX}") #here::here("data",glue("output_{SUFFIX}")) 
  system(glue('mkdir -p {DIR_OUTPUT}'))

  ## this is where all the written results go
  DIR_SCRATCH <-glue("../results/results_{SUFFIX}")
  system(glue('mkdir -p {DIR_SCRATCH}'))
  
  # Filenames ---------------------------------------------------------------
  
  ## HOST
  FILE_HOST_IN <- glue::glue("{DIR_HOST}/hbv_gilead")
  HOST <- fs::path_file(FILE_HOST_IN)
  FILE_HOST_tmp0 <- glue::glue("{DIR_PROCESSED}/{HOST}_tmp0")
  FILE_HOST_tmp1 <- glue::glue("{DIR_PROCESSED}/{HOST}_tmp1")
  FILE_HOST_OUT <- glue::glue("{DIR_PROCESSED}/{HOST}_QC")
  
  ## PATHOGEN
  FILE_PATHOGEN_AA_IN <-
    glue::glue("{DIR_PATHOGEN}/OUT_aa_binary_table.txt")
  FILE_PATHOGEN_NT_IN <-
    glue::glue("{DIR_PATHOGEN}/OUT_nt_binary_table.txt") ## nucleotide
  FILE_PATHOGEN_OUT <- glue::glue("{DIR_PROCESSED}/hbv-data.pheno")
  FILE_PATHOGEN_TREE <- glue::glue("{DIR_PATHOGEN}/HBV_WG+og_rr.nw")
  
  FILE_GRM_OUT <- glue::glue("{DIR_PROCESSED}/host-grm")
  
  ## COVARS
  FILE_COVARS_IN <- glue::glue("{DIR_RAW}/to_fellay_manifest.csv")
  FILE_COVARS_NUMERIC_OUT <-
    glue::glue("{DIR_PROCESSED}/covars-numeric.covar")
  FILE_COVARS_DISCRETE_OUT <-
    glue::glue("{DIR_PROCESSED}/covars-discrete.covar")
  
  NAM_PLINK <- paste0("plink_aa")#, fs::path_file(FILE_PATHOGEN_AA_IN))
  NAM <- "gcta_aa"
  
  ##////////////////////////////////////////////////////////////////
  ##                      define binaries                         //
  ##////////////////////////////////////////////////////////////////
  
  
  VCFTOOLS <- 'vcftools' #here::here("bin", "vcftools_0.1.13", "bin", "vcftools")
  VCFMERGE <- 'vcf-merge' #here::here("bin", "vcftools_0.1.13", "bin", "vcf-merge")
  BCFTOOLS <- 'bcftools' #here::here("bin", "bcftools-1.9", "bin", "bcftools")
  
  PLINK <- 'plink2' #here::here("bin", "plink2_linux")
  PLINK1 <- 'plink' #here::here("bin", "plink1.9_linux")
  KING <- './bin/king'
  
  ## gcta:
  ## help: https://cnsgenomics.com/software/gcta/#ManipulatingtheGRM
  GCTA <- 'gcta64' # here::here("bin", "gcta64_1.91.7_linux.beta")
  
  
  ##////////////////////////////////////////////////////////////////
  ##                      define parameters                       //
  ##////////////////////////////////////////////////////////////////
  
  ## how many pcs
  n_pc <- 4 ## keep at 2
  n_ppc <- 6
  
  ## thresholds discussed with Christian Thorball
  if(burden){
    maf_thresh <- 0.00001 
    maf_max <- 0.499
  }else{
    maf_thresh <- 0.05 ## if 0.01 many high odds ratios.
  }
  hwe_thresh <- 1e-6
  
  ## setting SNP threshold low after using an extra filtering step with BCFTOOLS
  callrate_snps_thresh <- 0.1
  ## but leavign the individual one kinda lenient
  callrate_ids_thresh <- 0.1
  
  #--geno [maximum per-variant]
  #--mind [maximum per-sample]
  #--geno filters out all variants with missing call rates exceeding the provided value (default 0.1) to be removed, while --mind does the same for samples.
  
  
  callrate_path_variant_thresh <- 0.1
  callrate_path_ids_thresh <- 0.1
  
  mac_path_variant_thresh_tmp <- 0 ## used in run-multiskat.R
  #If only testign one genotype, we lower the threshold to ensure that the specified variants are tested. 
  if(length(GT) == 1){
    mac_path_variant_thresh <- 1
  }else{
    mac_path_variant_thresh <- 20
  }

  p.THRESH.lenient <- 5e-8 ## checked with jacques
  
  HBV_GT <- GT ## sets GT in prep-data.R
  POP <- POP
  
  ##////////////////////////////////////////////////////////////////
  ##                     custom functions                         //
  ##////////////////////////////////////////////////////////////////
  
  
  extract_annotation_vcf <- function(path_to_vcfs, extract = c("fix","info","gt"))
  {
    dat_annotation <- path_to_vcfs %>% purrr::map_dfr(extract_vcf_custom, extract = extract) %>% dplyr::distinct()# %>% ## assembles all annotation, makes it distinct 
    #   assertr::assert(assertr::is_uniq, POS) ## checks if there are duplicated positions
    
    return(dat_annotation)
  }
  
  
  
  
  
  #' Turn VCF into DATA FRAME
  #' Extract per VCF file a subset of columns
  #'
  #' @param filename filename to vcf or vcf.gz file
  #'
  #' @return dataframe with variants by row and and annotation as column (but no genotype data)
  #'
  #' @examples
  #' ## get G04860.var.homo.SNPs.vcf from Sinergia switch drive
  #' dat.vcf <- extract_vcf_custom("data/raw/Mtb_VCF_files/G04860.var.homo.SNPs.vcf")
  #' 
  #' 
  extract_vcf_custom <- function(filename, extract = c("header","annotation","gt"))
  {
    ## turn data into a R readtable vcf format
    vcf <- vcfR::read.vcfR(filename)
    
    
    if (extract %in% "gt")
    {
      
      ## GT data ---------------------------------------------
      gt <- vcfR::extract_gt_tidy(vcf) %>%
        dplyr::select(Indiv, gt_GT, gt_GT_alleles, gt_FREQ) %>% ## only select ID, gt_GT, gt_GT in alleles and GT freq 
        ## strwrap(vcf@meta) shows that FREQ is the variant allele frequency (the one of the alternative allele, not the reference allele)
        dplyr::mutate(gt_GT_binary = dplyr::case_when(
          gt_GT == "1/1" ~ 1,
          gt_GT == "0/0" ~ 0,
          TRUE ~ NA_real_
        )) ## turn 1/1 into 1 and 0/0 into 0
      
      ## combine ---------------------------------------------
      if(exists("out")){
        out <- cbind(out, gt)
      }else{
        out <- gt
      }
      
    }
    
    if(extract %in% "info")
    { 
      ## SNP info data ---------------------------------------
      ann_cols <- vcfR::queryMETA(vcf, element = 'INFO=<ID=ANN')[[1]][4] %>% 
        stringr::str_remove_all("\\t") %>% ## remove all tabs
        stringr::str_remove_all("Description=Functionalannotations:'") %>% ## remove pretext
        stringr::str_remove_all("'") %>% ## remove all apostrophes
        stringr::str_split(pattern = "\\|") %>% # split by pipe
        unlist()
      info <- vcfR::extract_info_tidy(vcf) %>% dplyr::select(ANN, LOF, NMD) 
      
      info_ann <- info$ANN %>% stringr::str_split(",") %>%
        purrr::map(1) %>% ## first, only the relevant ones, the first ones
        stringr::str_split("\\|") %>% ## split by pipe
        do.call(rbind.data.frame, .) %>% ## turn into dataframe
        magrittr::set_colnames(ann_cols) %>% ## add colnames ann_cols
        dplyr::mutate_if(is.factor, as.character)
      ## make sure no factor!!!
      
      info %<>% dplyr::select(-ANN) %>% dplyr::bind_cols(info_ann)
      ## queryMETA(vcf, element = 'INFO=<ID=HOM')
      ## removing ADP, WT, HET, HOM and NC (e.g. all samples will be homozygous)
      ## queryMETA(vcf, element = 'INFO=<ID=LOF')
      ## "##INFO=<ID=LOF,Number=.,Type=String,Description=\"Predicted loss of function effects for"    
      ## "this variant. Format: 'Gene_Name | Gene_ID | Number_of_transcripts_in_gene |"                
      ## "Percent_of_transcripts_affected' \">"                                                        
      ## "##INFO=<ID=NMD,Number=.,Type=String,Description=\"Predicted nonsense mediated decay effects" 
      ## "for this variant. Format: 'Gene_Name | Gene_ID | Number_of_transcripts_in_gene |"            
      ## "Percent_of_transcripts_affected' \">
      
      ## combine ---------------------------------------------
      if(exists("out")){
        out <- cbind(out, info)
      }else{
        out <- info
      }
      
    }
    
    if(extract %in% "fix")
    { 
      ## FIX body ---------------------------------------------
      fix <- vcfR::getFIX(vcf)[, c("CHROM", "POS", "ID", "REF", "ALT")] %>% 
        as.data.frame() %>% 
        dplyr::mutate(POS = as.numeric(as.character(POS))) %>%
        dplyr::mutate_if(is.factor, as.character)
      
      ## combine ---------------------------------------------
      if(exists("out")){
        out <- cbind(out, fix)
      }else{
        out <- fix
      }
    }
    
    return(out)
  }
  
  
  str_ignore <- function(string, pattern) {
    string[!str_detect(string, pattern)]
  }
  
  
  
  
  sample_size_extract <-
    function(files = NULL, 
             phenotype = "plink_aa_OUT_aa_binary_table.txt")
    {
      DIR <- fs::path_dir(files[1])
      
      ## find those that errored
      not_run <- purrr::map(files, function(x)
        system(glue::glue("grep 'error occurs' {x}"), intern = TRUE)) 
      not_run <- lapply(not_run, function(x) case_when(
        length(x) == 0 ~ FALSE,
        length(x) > 0 ~ TRUE))
      not_run_df <- not_run %>% unlist() %>% data.frame() 
      names(not_run_df) <- "error_run"
      not_run_df$files <- files
      not_run_df$outcome <- files %>% 
        stringr::str_replace(DIR, "") %>% 
        stringr::str_replace(glue::glue("{phenotype}_"), "") %>% 
        stringr::str_replace(".log", "")
      
      
      ## sample size  
      files_sub <- not_run_df %>% filter(!error_run) %>% pull(files)
      sample_size <-
        purrr::map(files_sub, function(x)
          system(glue::glue("grep -F 'observations' {x}"), intern = TRUE) %>% first() %>% stringr::str_split(" ") %>% purrr::map(1)) %>% unlist() %>% as.numeric()
      
      sample_size_df <-
        tibble::tibble(files = files_sub, 
                       n = sample_size) %>% mutate(outcome = files %>% stringr::str_replace(DIR, "") %>% stringr::str_replace(glue::glue("/{phenotype}_"), "") %>% stringr::str_replace(".log", ""))
      
      ## assemble it together
      
      sample_size_not_run_df <- full_join(sample_size_df %>% select(-files), not_run_df %>% select(-files), by = "outcome") %>% select(outcome, n, error_run)
      
      return(sample_size_not_run_df)
    }
  
  sample_size_extract_gcta <-
    function(files = NULL, 
             phenotype = "gcta_aa_OUT_aa_binary_table.txt")
    {
      DIR <- fs::path_dir(files[1])
      
      ## find those that errored
      not_run <- purrr::map(files, function(x)
        system(glue::glue("grep 'error occurs' {x}"), intern = TRUE)) 
      not_run <- lapply(not_run, function(x) case_when(
        length(x) == 0 ~ FALSE,
        length(x) > 0 ~ TRUE))
      not_run_df <- not_run %>% unlist() %>% data.frame() 
      names(not_run_df) <- "error_run"
      not_run_df$files <- files
      not_run_df$outcome <- files %>% 
        stringr::str_replace(DIR, "") %>% 
        stringr::str_replace(glue::glue("{phenotype}_"), "") %>% 
        stringr::str_replace(".log", "")
      
      
      ## sample size  
      files_sub <- not_run_df %>% filter(!error_run) %>% pull(files)
      sample_size <-
        purrr::map(files, function(x)
          system(glue::glue("grep -F 'observations' {x}"), intern = TRUE) %>% first() %>% stringr::str_split(" ") %>% purrr::map(1)) %>% unlist() %>% as.numeric()
      
      sample_size_df <-
        tibble::tibble(files = files_sub, 
                       n = sample_size) %>% mutate(outcome = files %>% stringr::str_replace(DIR, "") %>% stringr::str_replace(glue::glue("{phenotype}_"), "") %>% stringr::str_replace(".log", ""))
      
      ## assemble it together
      
      sample_size_not_run_df <- full_join(sample_size_df %>% select(-files), not_run_df %>% select(-files), by = "outcome") %>% select(outcome, n, error_run)
      
      return(sample_size_not_run_df)
    }
  
  ## //////////////////
  ## PLOTs
  ## //////////////////
  qqplot_wrapper <- function(x, title = "", subtitle = "")
  {
    N <- length(x)
    
    ## calculate the expected axis
    expected <-
      sort(-log10((1:N) / N - 1 / (2 * N)))
    observed <-
      sort(-log10(x))
    
    qp <- ggplot() +
      geom_abline(intercept = 0,
                  slope = 1,
                  alpha = I(0.5)) +
      geom_point(aes(expected, observed), shape = 1) +
      #   geom_point(aes(expected, observed), shape = 16, alpha = I(0.3)) +
      xlab(expression(Expected ~ ~ -log[10](italic(p)))) +
      ylab(expression(Observed ~ ~ -log[10](italic(p)))) +
      labs(title = title, subtitle = subtitle)
    
    return(qp)
  }
  
  mhtplot_wrapper <-
    function(data , title, subtitle, THRESH1, THRESH2, color_highlight = "red") {
      ## chr
      ## pos
      ## p
      
      
      ## add small space
      ## equidistance
      data_max <-
        data %>% group_by(chr) %>% slice(which.max(pos)) %>% mutate(pos_max = TRUE)
      data <-
        left_join(data, data_max) %>% mutate(pos_max = case_when(is.na(pos_max) ~ FALSE, TRUE ~ pos_max))
      
      CONST <- 10000
      data <- data %>% dplyr::arrange(chr, pos) %>%
        dplyr::mutate(tmp = 1) %>%
        mutate(tmp = case_when(pos_max ~ tmp + CONST,
                               TRUE ~ tmp)) %>%
        mutate(cumsum.tmp = cumsum(tmp))
      
      ## real distance
      med.dat <-
        data %>% group_by(chr) %>% summarise(median.x = median(cumsum.tmp))
      
      ## signif snps
      CONST <- 4e5 ## same as in src/run_locuszoomplot
      data_snps_signif <-
        data %>% filter(p <= THRESH1) %>% select(chr, pos)
      data_snps_signif_range <-
        purrr::map2_dfr(data_snps_signif$chr, data_snps_signif$pos, function(x, y)
          data %>% filter(chr == x & pos < (y + CONST) & pos > (y - CONST)))
      
      
      if (nrow(data_snps_signif) > 0)
      {
        qp <- ggplot() +
          geom_point(data = data,
                     aes(cumsum.tmp,-log10(p), color = as.factor(chr)),
                     shape = 1) +
          #geom_point(data = data, aes(cumsum.tmp, -log10(p), color = as.factor(chr)),  shape = 16, alpha = I(0.5)) +
          labs(title = title) + xlab("Chromosomal Position") +
          ylab(expression(-log[10](italic(p)))) +
          scale_x_continuous(breaks = med.dat$median.x, labels = med.dat$chr,limits = c(0,max(med.dat$median.x + 10))) +
          geom_hline(yintercept = -log10(THRESH1),
                     linetype = 2) +
          geom_hline(
            yintercept = -log10(THRESH2),
            alpha = I(0.5),
            linetype = 2
          ) +
          scale_color_manual(values = rep(c("black", gray(0.4)), 11)) + theme(legend.position = "none") +
          geom_point(aes(cumsum.tmp,-log10(p)),
                     color = I(color_highlight),
                     data = data_snps_signif_range)
        
        
      } else{
        qp <- ggplot() +
          geom_point(data = data,
                     aes(cumsum.tmp,-log10(p), color = as.factor(chr)),
                     shape = 1) +
          #geom_point(data = data, aes(cumsum.tmp, -log10(p), color = as.factor(chr)),  shape = 16, alpha = I(0.5)) +
          labs(title = title) + xlab("Chromosomal Position") +
          ylab(expression(-log[10](italic(p)))) +
          scale_x_continuous(breaks = med.dat$median.x, labels = med.dat$chr) +
          geom_hline(yintercept = -log10(THRESH1),
                     linetype = 2) +
          geom_hline(
            yintercept = -log10(THRESH2),
            alpha = I(0.5),
            linetype = 2
          ) +
          scale_color_manual(values = rep(c("black", gray(0.4)), 11)) + theme(legend.position = "none")
        
      }
      return(qp)
      
      
    }
  
  
  
  #' EQTLs from EUGENE
  #'
  #' @param rsid
  #'
  #' @details from Christian Thorball and https://genepi.qimr.edu.au/staff/manuelF/eugene/download.html
  #'
  #' @return dataframe SNP, eugene.gene, eugene.pvalue, eugene.tissue
  #'
  #' @note only rsid possible
  #'
  #' @examples
  #' extract_eugene(rsid = "rs12597715")
  #'
  extract_eugene <- function(rsid)
  {
    if (str_detect(rsid, "rs"))
    {
      ## Tissue, Pvalue, Gene
      tmp <-
        system(glue::glue("grep -w {rsid} /home/rueger/EQTL_DB/EUGENE/eqtl*.cis"),
               intern = TRUE) %>%
        str_split(" ") %>%
        unlist() %>%
        data.frame() %>%
        t() %>%
        data.frame()
      
      if (nrow(tmp) > 0)
      {
        ## melt
        n_long <- 7 ## file, gene, pval, tissue, type, allele, something
        tmp_vec <- sapply(tmp, as.character)
        tmp_long <- matrix(tmp_vec, byrow = TRUE, ncol= n_long) %>% data.frame() %>%
          mutate_if(is.factor, as.character)
        tmp_long$SNP <- rsid
        names(tmp_long) <- c(
          "file",
          "eugene.gene",
          "eugene.pvalue",
          "eugene.tissue",
          "type",
          "allele",
          "something",
          "SNP")
        
        ## to numeric
        tmp_long$eugene.pvalue <- as.numeric(tmp_long$eugene.pvalue)
        tmp_long$something <- as.numeric(tmp_long$something)
        
        tmp_long %<>% dplyr::select(SNP, eugene.gene, eugene.pvalue, eugene.tissue, type) %>% arrange(eugene.pvalue)
        
        return(tmp_long)
        
      }
    }
    
  }
  save(list = ls(all.names = TRUE), file =  glue("{DIR_SCRATCH}/setup.rda", envir = 
         environment()))
}
args <-  commandArgs(trailingOnly = T)
POP <- strsplit(args[[1]],split = '_')[[1]]
PCA_type <- args[[2]]
GT <- strsplit(args[[3]],split = '_')[[1]]
LOCO <- as.logical(args[[4]])
burden <- as.logical(args[[5]])
n_pc <- as.numeric(args[[6]])
n_ppc <- as.numeric(args[[7]])
DIR_RAW <- args[[8]]
n_cores <- as.numeric(args[[9]])

RunSetUp(POP,GT,PCA_type,LOCO,burden,n_pc,n_ppc,DIR_RAW,n_cores)