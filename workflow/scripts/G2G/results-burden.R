library(parallel)
library(pbmcapply)

viral_path <- list.dirs('~/G2G-HBV/data/results_BURDEN_POP_asian_GT_A_C_D_B_TYPE_pPCA_LOCO_FALSE/')[-1]
gene_files <- pbmclapply(viral_path,function(x) dir(x)[grepl(dir(x),pattern = '.out') & !grepl(dir(x),pattern = '_single')],mc.cores = 20)
names(gene_files) <- viral_path

gene_files <- gene_files[!sapply(gene_files,function(x) length(x) == 0)]
res <- data.frame(Pathogen_Variant = character(),Host_Gene = character(),P = numeric())

for(i in 1:length(gene_files)){
  print(i)
  cur_dir <- names(gene_files)[i]
  cur_files <- gene_files[[i]]
  cur_variant <- strsplit(cur_dir,split = 'gene_')[[1]][2]
  
  gene_results <- mclapply(1:length(cur_files),function(k){
    cur_gene <- gsub(strsplit(cur_files[k],split = '_')[[1]][length(strsplit(cur_files[k],split = '_')[[1]])],pattern = '.out',replacement = '')
    cur_file <- readLines(paste0(cur_dir,'/',cur_files[k]))
    if(length(cur_file) > 0){
      p_value_numeric <- as.numeric(strsplit(cur_file[2],split = ' ')[[1]][2])
    }else{
      p_value_numeric <- NA
    }
    return(list(cur_gene=cur_gene,p_value_numeric=p_value_numeric))
  },mc.cores = 20)
  res <- rbind(res,data.frame(Pathogen_Variant=rep(cur_variant,length(cur_files)),Host_Gene = sapply(gene_results,function(x) x$cur_gene),P = sapply(gene_results,function(x) x$p_value_numeric)))
  
}
saveRDS(res,'~/G2G-HBV/data/results_BURDEN_POP_asian_GT_A_C_D_B_TYPE_pPCA_LOCO_FALSE/results.rds')