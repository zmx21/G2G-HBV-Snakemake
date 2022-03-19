RunResultsG2G <- function(POP,GT,PCA_type,LOCO,n_cores,debugging = F,sigHitsOnly=F,burden = F,p.thresh=1e-4){
  ##////////////////////////////////////////////////////////////////
  ##                        Dependencies                          //
  ##////////////////////////////////////////////////////////////////
  library(pbmcapply)
  source(here::here("src", "setup.R"))
  #Load Setup stored Variables
  DIR_SCRATCH_New <- RunSetUp(POP,GT,PCA_type,LOCO,n_cores,burden)
  load(glue("{DIR_SCRATCH_New}/setup.rda"))
  load(glue("{DIR_SCRATCH_New}/prep-data.rda"))
  DIR_SCRATCH <- DIR_SCRATCH_New
  ## add string to debugging
  if (debugging) {
    NAM <- glue::glue("{NAM}_debugging")
  }
  
  ## files ----------------------------------------------------------------
  

    files <- fs::dir_ls(DIR_SCRATCH) %>%
      str_subset("SAIGE_SPA_TEST_")
    n_tests <- length(files) ## run munge/eff-number-phenotypes.R
    THRESH <- p.THRESH.lenient / (n_tests)
  
  
  ##////////////////////////////////////////////////////////////////
  ##                  Get Summary Stats                           //
  ##////////////////////////////////////////////////////////////////
  
  
  ## read in all data -------------------------------------------------------
  all_results  <- pbmclapply(files,function(x) data.table::fread(x,header = T) %>% dplyr::filter(p.value <= p.thresh),mc.cores = 30)
  names(all_results) <- sapply(strsplit(files,split = glue::glue("{LOCO}_")),function(x) x[2])
  saveRDS(all_results,file = glue::glue("{DIR_SCRATCH}/all_results.rds"))
  below_thresh <- sapply(all_results,function(x) nrow(x %>% dplyr::filter(p.value <= 5e-8)) != 0)
  if(all(!below_thresh)){
    return(NULL)
  }
  all_results_filt <- all_results[which(below_thresh)]
  all_results_filt <- lapply(all_results_filt,function(x) dplyr::filter(x,p.value < 5e-8))
  # dat_anno <- data.table::fread(glue::glue("{DIR_RAW}/wes_plink/snpEff_annotation_processed.csv"))
  
  all_results_filt_merged <- lapply(all_results_filt,function(dat_select){
    ## add annotation
    dat_select <- dat_select #%>% left_join(dat_anno, by = c("CHR" = "#CHROM", "POS" = "POS"))

    ## add eqtl
    # annotation_raw_eugene <- parallel::mclapply(dat_select$SNPID, function(x) if(str_detect(x, "rs")) extract_eugene(rsid = x), mc.cores = 12)
    # annotation_eugene <- tryCatch(annotation_raw_eugene %>% 
    #                                 bind_rows() %>% 
    #                                 group_by(SNP) %>%
    #                                 slice(which.min(eugene.pvalue)) %>% 
    #                                 summarize(eugene = paste(glue::glue("{eugene.gene}-{type} in {eugene.tissue} (P={eugene.pvalue})"), collapse = ", ") ) %>%
    #                                 ungroup(),error = function(e) data.frame(SNP = dat_select$SNPID,eugene = NA))
    # dat_select <- dat_select %>% left_join(annotation_eugene,by = c('SNPID' = 'SNP'))
    return(dat_select)
  })
  
  ##////////////////////////////////////////////////////////////////
  ##                         3a. Q-Q Plots                        //
  ##                         3b. MHT Plots                        //
  ##////////////////////////////////////////////////////////////////
  
  ## loop through all results to produce QQ and MHT plots
  
  for (i in 1:length(all_results_filt_merged)) {
    VARIANT <- names(all_results_filt_merged)[i]

    out <-
      data.table::fread(files[which(grepl(VARIANT,files))])  %>% 
      rename(SNP = SNPID,chr = CHR, pos = POS,p=p.value) %>% 
      filter(!is.na(p))
    
    
    
    ##----------------------------
    ## Lambda QC adaptation
    ##----------------------------
    
    lambda_gc <- GWAS.utils::genomic_inflation(P = na.omit(out$p))
    
    ##----------------------------
    ##  Q-Q plot                --
    ##----------------------------
    
    qp.qq <- qqplot_wrapper(x = out %>% pull(p), subtitle = glue::glue("GC lambda = {round(lambda_gc, 2)}"), title = glue::glue("{VARIANT}"))
    
    
    png(glue::glue(
      "{DIR_SCRATCH}/results_top_qqplot_SAIGE_{VARIANT}.png"
    ))
    print(qp.qq)
    dev.off()
    
    ##----------------------------
    ##  MHT plot                --
    ##----------------------------
    
    png(glue::glue("{DIR_SCRATCH}/results_top_manhattan_SAIGE_{VARIANT}.png"), width = 900/9*8, height = 480/9*8)
    if (debugging) {
      qqman::manhattan(out,
                       chr = "chr",
                       bp = "pos",
                       p = "p")
      
    } else {
      qp <- mhtplot_wrapper(data = out %>% select(SNP, chr, pos, p), title = glue::glue("{VARIANT}"), THRESH1 = THRESH, THRESH2 = 5e-08)
      print(qp)
    }
    dev.off()
  }
}
# args <- commandArgs(trailingOnly = TRUE)
# if(length(args) > 0){
#   RunResultsG2G(POP = unlist(strsplit(args[1],split = ',')),GT =unlist(strsplit(args[2],split = ',')),PCA_type = args[3],LOCO = as.logical(args[4]),n_cores = as.numeric(args[5]),debugging = ifelse(args[6]=='forreal',F,T),sigHitsOnly = as.logical(args[7]))
# }
# RunResultsG2G(POP = c('asian','european'),GT = c('A','B','C','D','F','H'),PCA_type = 'pPCA',LOCO = F,n_cores = 1,debugging = F,sigHitsOnly = F,burden = F)
RunResultsG2G(POP = 'asian',GT = 'C',PCA_type = 'pPCA',LOCO = F,n_cores = 1,debugging = F,sigHitsOnly = F,burden = F)
RunResultsG2G(POP = 'asian',GT = 'B',PCA_type = 'pPCA',LOCO = F,n_cores = 1,debugging = F,sigHitsOnly = F,burden = F,p.thresh=1)
RunResultsG2G(POP = 'asian',GT = c('A','C','D','B'),PCA_type = 'pPCA',LOCO = F,n_cores = 1,debugging = F,sigHitsOnly = F,burden = F)
# RunResultsG2G(POP = 'european',GT = 'D',PCA_type = 'pPCA',LOCO = F,n_cores = 1,debugging = F,sigHitsOnly = F,burden = F)
RunResultsG2G(POP = 'european',GT = c('A','C','D','F','H'),PCA_type = 'pPCA',LOCO = F,n_cores = 1,debugging = F,sigHitsOnly = F,burden = F)
