

## Data preparation HBV data
## - load data
## - apply QC
## - store in processed/


#' QC for Pathogen data
#'
#' @param data dataset with id and genetic variants as columns
#' @param NA_ids callrate threshold for individuals, e.g. 0.03
#' @param NA_variants callrate threshold for SNPs, e.g. 0.01
#' @param MAC minor allele count threshold, e.g. 10, removes all variants with counts equal or less than 10
#'
#' @return cleaned dataset
#'
#' @examples
#' dat_pathogen_raw <- read_data_pathogen("OUT_nt_binary_table.txt", dat_ids = x)
#' dat_pathogen <- apply_qc_pathogen(dat_pathogen_raw, NA_id = callrate_ids_thresh, NA_variants = callrate_snps_thresh)

apply_qc_pathogen <-
  function(data = NULL,
           NA_ids = 0,
           NA_variants = 0,
           MAC = 0)
  {
    ## snps
    NA_freq_by_snps <-
      apply(is.na(data), 2, sum, na.rm = TRUE) / nrow(data)
    rm_snps <- which(NA_freq_by_snps >= NA_variants)
    if (length(rm_snps) > 0) 
    {
      data <- data[,-rm_snps]
    }
    cat("rm ", length(rm_snps), "HBV variants \n")
    
    
    ## ids
    NA_freq_by_ids <-
      apply(is.na(data), 1, sum, na.rm = TRUE) / ncol(data)
    rm_ids <- which(NA_freq_by_ids >= NA_ids)
    if (length(rm_ids) > 0) 
    {
      data <- data[-rm_ids,]
    }
    cat("rm ", length(rm_ids), "individuals \n")
    
    
    ## snps
    counts_by_snps_0 <- apply(data[, -1] == 0, 2, sum, na.rm = TRUE)
    counts_by_snps_1 <- apply(data[, -1] == 1, 2, sum, na.rm = TRUE)
    
    rm_c_snps <- which(counts_by_snps_0 < MAC | counts_by_snps_1 < MAC)
    if (length(rm_c_snps) > 0)
    {
      data <-
        data[,-(rm_c_snps + 1)] ## cause shifted by one bc no ID included
    }
    cat("rm ", length(rm_c_snps), "HBV variants bc of low counts \n")
    
    
    return(data)
  }



#' Read Pathogen data
#' using the pathogen file, indicated with PATH
#' and matching it with an ID file, which has two columns: host_id, pathogen_id
#'
#' @param PATH path to pathogen data
#' @param DAT_IDS dataframe with host_id and pathogen_id as column
#'
#' @return dataset with ids that will match covars and host data
#'
#' @examples
#' dat_ids <- data.frame(host_id = c("igm160234", "igm160214", "igm160231", "igm160235", "igm160219", "igm160228"), pathogen_id = c("GS-US-320-0108-1001", "GS-US-320-0108-1002", "GS-US-320-0108-1003", "GS-US-320-0108-1005", "GS-US-320-0108-1006", "GS-US-320-0108-1007"))
#' out <- read_data_pathogen(glue::glue("{DIR_PATHOGEN}/OUT_nt_binary_table.txt"), DAT_IDS = dat_ids)

read_data_pathogen <- function(PATH, DAT_IDS) {
  data <- read.table(PATH, row.names = 1) %>%
    # as_tibble() %>%
    tibble::rownames_to_column("ID")
  
  ## add the correct IDs from host
  data <- DAT_IDS %>%
    left_join(data, c("pathogen_id" = "ID")) %>%
    dplyr::select(-pathogen_id) %>%
    dplyr::rename(ID = host_id)
  
  
  ## Sanitychecks
  if (!all(unique(unlist(data[, -1])) %in% c(NA, 0, 1)))
    stop("Other values than 0 and 1 in dataset")
  
  ## return
  return(data)
}



#In terms of viral data format, the numbering for all sequences is relative to HBV GT C (3215 bp), and the AA and nucleotide frequencies for variants are reported as binary 0 or 1, indicating the absence or  presence of a certain variant, respectively. When there is an NA in the frequency tables ("-" represents a gap) in the consensus sequence, that indicates that a variant could not be called at that position due to low coverage (<100 read coverage). Also, in the consensus sequences, we included mixed bases (IUPAC code) at a position when variants were detected at 15% or greater. In the binary tables, each Nuc column is formatted according to “pos_position123_variant” and each AA column is formatted according to “gene_genenameXYZ_pos_position123_variant.” Please note that the HBV genome has overlapping genes. For the consensus sequences reported, they are provided in an alignment and named according to USUBJID.


#' Extract PCs
#' This is a function to compute the phylogenetic PCs for the HBV project
#' Guide on constructing pPCs: https://gist.github.com/sinarueeger/718c7ec3cee473ad0beed4feea4d0e24
#'
#' @param path_tree path to tree
#' @param path_pathogen path to pathogen data
#' @param ppca logical, returning phylogenetic PCs
#' @param pca logical, returning PCs
#' @param n.pc number of PCs returned
#' @param loadings logical, returining loadings in case of pca = TRUE or ppca = TRUE
#' @param eigenvalues logical, returining eigenvalues in case of pca = TRUE
#' @param filter_threshold maf filter, everything < filter_threshold will be removed
#' @param id.list character vector of IDs to be incorporated
#'
#' @return dataframe with ID and PCs as columns
#' @examples
#' pathogen_pca(
#' path_tree = glue::glue("{DIR_PATHOGEN}/HBV_WG+og_rr.nw"),
#' path_pathogen = glue::glue("{DIR_PATHOGEN}/OUT_aa_binary_table.txt"),
#' pca = TRUE
#' )



pathogen_pca <- function(path_tree,
                         path_pathogen,
                         ppca = FALSE,
                         pca = TRUE,
                         n.pc = 4,
                         loadings = FALSE,
                         eigenvalues = FALSE,
                         filter_threshold = 0.05,
                         id.list = NULL)
{
  
  ## 1) TREE --------------------------------------
  ## ----------------------------------------------
  tree <- ape::read.tree(path_tree)
  #plot(tree)
  
  ## check if matrix is singular
  #solve(ape::vcv.phylo(tree))
  
  ## 2) Y matrix ----------------------------------
  ## ----------------------------------------------
  
  data_pathogen <- read.table(path_pathogen, row.names = 1) %>%
    tibble::rownames_to_column("ID")
  
  if (!is.null(id.list))
  {
    data_pathogen <- data_pathogen %>%
      filter(ID %in% id.list)
    
  }
  
  tree_subset <-
    ape::keep.tip(tree, as.character(data_pathogen$ID))
  #    plot(tree_subset)
  
  ## reorder rows, same tiplabels
  id_tree <- tibble(ID = tree_subset$tip.label)
  
  ## turn data into matrix and reorder rows according to tips
  Y <- left_join(id_tree,data_pathogen)
  
  stopifnot(identical((tree_subset$tip.label), as.character(Y$ID)))
  
  
  id_pathogen <- Y$ID
  
  ## remove ID
  Y_no_id <- Y %>% dplyr::select(-ID) %>% as.matrix()
  
  ## remove columns with only wild type
  Y_no_id <-
    Y_no_id[, which(apply(Y_no_id, 2, function(x)
      mean(x, na.rm = TRUE)) != 0)]
  
  ## Remove variant between 0.05 and 0.95 allele frequency
  Y_no_id.freq <-
    parallel::mclapply(1:ncol(Y_no_id), function(x) {
      mean(Y_no_id[, x], na.rm = T)
    }, mc.cores = 20) %>% unlist()
  Y_no_id <- Y_no_id[, which(Y_no_id.freq > filter_threshold & Y_no_id.freq < (1-filter_threshold)) ]
  
  ## only non varying ones
  Y_no_id.var <-
    parallel::mclapply(1:ncol(Y_no_id), function(x) {
      var(Y_no_id[, x], na.rm = T)
    }, mc.cores = 20) %>% unlist()
  Y_no_id <- Y_no_id[, which(Y_no_id.var != 0)]
  
  ## impute
  message("Imputing PCA matrix... (will take some time)")
  # Y_no_id4logitPCA <- round(missMDA::imputePCA(Y_no_id)$completeObs)
  

  
  ## 3) pPCA --------------------------------------
  ## ----------------------------------------------
  
  ## initiate out
  out <- data.frame(ID = id_pathogen)
  
  #print(sessionInfo())
  
  if (ppca)
  {
    ## from Nimisha
    vir_4d <- phylobase::phylo4d(tree_subset, Y_no_id)
    
    vir_pca <- adephylo::ppca(
      vir_4d,
      scale = TRUE,
      ## if scaled, very similar to pca
      scannf = FALSE,
      nfposi = 16,
      method = "Abouheif"
    )
    
    
    if (loadings)
    {
      ## return loadings
      ppca_out <- vir_pca$c1[, 1:n.pc]
      rownames(ppca_out) <- colnames(Y_no_id)
      return(ppca_out)
    } else{
      if (eigenvalues)
      {
        ppca_out <- vir_pca$eig
        return(ppca_out)
        
      } else {
        ## return ppca
        ppca_out <- vir_pca$li %>% dplyr::select(1:n.pc)
        names(ppca_out) <- paste0("pPC", 1:n.pc)
        
      }
    }
    
    out <- cbind(out, ppca_out)
  }
  
  
  ## 4) as comparison: PCA ------------------------
  ## ----------------------------------------------
  if (pca)
  {

    ## FactoMineR::PCA ---------------------------
    message("Estimating logistic PCA... (will take some time)")
    pca_logistic_raw <-
      logisticPCA::convexLogisticPCA(Y_no_id, k = n.pc)
    ## PCs (PC scores)
    ## U = loadings
    
    if (loadings)
    {
      ## return loadings
      out <-
        pca_logistic_raw$U %>% as.data.frame()  %>% dplyr::select(1:n.pc)
      
    } else{
      if (eigenvalues) {
        ## only for eigenvalues, a bit dodgy tho to compute this when logistic PCA was used
        pca_raw <-
          FactoMineR::PCA(Y_no_id4logitPCA,
                          scale.unit = FALSE,
                          graph = FALSE)
        
        ## return eigenvalues
        pca_out_eigenvalues <- (pca_raw$eig)#[1:n.pc]
        out <- pca_out_eigenvalues
        
      } else{
        pca_out <-
          pca_logistic_raw$PCs %>% as.data.frame() %>% dplyr::select(1:n.pc)
        names(pca_out) <- paste0("pPC", 1:n.pc)
        out <- cbind(out, pca_out)
      }
    }
    
  }
  
  
  return(out)
  
  
}



#' Summarise pathogen data with maf for amino acids
#'
#' @param data dataframe with pathogen data (ID as identifier)
#'
#' @return tibble
#'
#' @examples
#' summarise_pathogen(data = dat_pathogen)
#'
summarise_pathogen <- function(data = NULL) {
  data_long <- tidyr::gather(data, aa, value,-ID)
  
  counts_by_snps_0 <- apply(data[,-1] == 0, 2, sum, na.rm = TRUE)
  counts_by_snps_1 <- apply(data[,-1] == 1, 2, sum, na.rm = TRUE)
  
  
  ## calculate allele frequency, and allele counts
  out <-
    data_long %>% group_by(aa) %>% dplyr::summarise(
      af = mean(value, na.rm = TRUE),
      ac0 = sum(value == 0, na.rm = TRUE),
      ac1 = sum(value == 1, na.rm = TRUE)
    )
  
  ## turn ac into mac
  ## turn af into maf
  out <- out %>% mutate(
    maf = case_when(af > 0.5 ~ 1 - af,
                    af <= 0.5 ~ af,
                    TRUE ~ NA_real_),
    mac = case_when(ac1 > ac0 ~ ac0,
                    ac1 <= ac0 ~ ac1,
                    TRUE ~ NA_integer_)
  ) %>% dplyr::select(-af, -ac0, -ac1)
  
  return(out)
}
