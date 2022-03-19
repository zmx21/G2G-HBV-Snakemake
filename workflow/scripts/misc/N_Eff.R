#Calculates number of tests done based on output results

PruneAA <- function(dat_pathogen,prune_thresh = 1){
  cor_mat <- cor(dat_pathogen[,-1],use = 'complete.obs')
  
  locus_to_test <- colnames(dat_pathogen[,-1])
  locus_to_keep <- c()
  locus_to_remove <- c()
  while(length(locus_to_test) > 0){
    cur_locus <- locus_to_test[1]
    cor_with_locus <- cor_mat[cur_locus,]^2
    locus_to_remove <- c(locus_to_remove,setdiff(names(cor_with_locus)[cor_with_locus >= prune_thresh],cur_locus))
    locus_to_keep <- c(locus_to_keep,cur_locus)
    locus_to_test <- setdiff(locus_to_test[-1],locus_to_remove)
  }
  perf_pairs <- which(abs(cor_mat) == 1,arr.ind = T)
  perf_pairs <- perf_pairs[perf_pairs[,'row']!= perf_pairs[,'col'],]
  
  perf_pairs_df <- data.frame(SNP_1 = rownames(cor_mat)[perf_pairs[,'row']] ,SNP_2 = colnames(cor_mat)[perf_pairs[,'col']],stringsAsFactors = F)
  SNP_1_2 <- lapply(1:nrow(perf_pairs_df), function(i) sort(c(perf_pairs_df$SNP_1[i],perf_pairs_df$SNP_2[i])))
  perf_pairs_df <- perf_pairs_df[!duplicated(SNP_1_2),]
  return(list(dat_pathogen = dat_pathogen[,c('ID',locus_to_keep)],perf_pairs_df=perf_pairs_df))
}
asian_all_subtype_results_path <- '~/G2G-HBV/data/results_POP_asian_GT_A_C_D_B_TYPE_pPCA_LOCO_FALSE/'
eur_all_subtype_results_path <- '~/G2G-HBV/data/results_POP_european_GT_A_C_D_F_H_TYPE_pPCA_LOCO_FALSE/'

LOCO = 'FALSE'
sage_files_asn <- dir(asian_all_subtype_results_path)
sage_files_asn <- sage_files_asn[grepl(glue::glue("SAIGE_SPA_TEST_LOCO_{as.character(LOCO)}"),sage_files_asn)]
sage_genes_asn <- sapply(sage_files_asn,function(x) strsplit(x,split = glue::glue('SAIGE_SPA_TEST_LOCO_{as.character(LOCO)}_'))[[1]][2])

sage_files_eur <- dir(eur_all_subtype_results_path)
sage_files_eur <- sage_files_eur[grepl(glue::glue("SAIGE_SPA_TEST_LOCO_{as.character(LOCO)}"),sage_files_eur)]
sage_genes_eur <- sapply(sage_files_eur,function(x) strsplit(x,split = glue::glue('SAIGE_SPA_TEST_LOCO_{as.character(LOCO)}_'))[[1]][2])

load('~/G2G-HBV/data/results_POP_asian_GT_A_C_D_B_TYPE_pPCA_LOCO_FALSE/prep-data.rda')
dat_pathogen_pruned_asn <- PruneAA(dat_pathogen[,c('ID',intersect(sage_genes_asn,colnames(dat_pathogen)))],1)


load('~/G2G-HBV/data/results_POP_european_GT_A_C_D_F_H_TYPE_pPCA_LOCO_FALSE/prep-data.rda')
dat_pathogen_pruned_eur <- PruneAA(dat_pathogen[c('ID',intersect(sage_genes_eur,colnames(dat_pathogen)))],1)

N_Eff <- length(intersect(sage_genes_asn,colnames(dat_pathogen_pruned_asn$dat_pathogen))) + length(intersect(sage_genes_eur,colnames(dat_pathogen_pruned_eur$dat_pathogen)))
