library(seqinr)
library(pbmcapply)
library(latex2exp)
library(dplyr)
library(ggplot2)
library(ggrepel)

#Extracts all possible k-mers which overlaps a specific mutation
ExtractKMers <- function(k,mut_pos,prot_seq,residue){
  prot_seq <- toupper(prot_seq)
  prot_seq[mut_pos] <- residue
  offset <- 0:(k-1)
  spanning_seq <- lapply(offset,function(i) paste0(prot_seq[(mut_pos-i):(mut_pos + (k - 1 - i))],collapse = ''))
  return(unlist(spanning_seq))
  
}

RunPredictBinding <- function(out_dir,suffix,k_mer_list,hla_allele){
  system(glue::glue("mkdir -p {out_dir}"))
  ref_file <- glue::glue("{out_dir}ref_{suffix}")
  mut_file <- glue::glue("{out_dir}mut_{suffix}")
  system(glue::glue('rm {ref_file}*'))
  system(glue::glue('rm {mut_file}*'))
  pred_methods <- c('netmhcpan_ba','netmhcpan_el')

  #Write mut file
  for(i in 1:length(k_mer_list)) {
    cur_mut <- k_mer_list[[i]]$mut
    header <- sapply(1:nchar(cur_mut[1]),function(x) paste0('pos_',x))
    residue <- strsplit(x = cur_mut[1],split = '')[[1]][1]
    write(paste0(mapply(function(x,y) paste0('>',y,'\n',x),cur_mut,header),collapse = '\n'),glue::glue("{mut_file}_{residue}_k_{nchar(cur_mut[1])}.fasta"),append = F)
    #Run netmhcpan
    lapply(pred_methods,function(x) system(glue::glue("python ~/G2G-HBV/mhc_i/src/predict_binding.py {x} {hla_allele} {nchar(cur_mut[1])} {mut_file}_{residue}_k_{nchar(cur_mut[1])}.fasta > {mut_file}_{residue}_k_{nchar(cur_mut[1])}.{x}.out")))
    #Run MixMHCPred
    #Use nearest proxy for HLA-A*33:03 (according to Prof. David Gfeller)
    if(hla_allele == "HLA-A*33:03"){
      hla_allele_mix_mhc_pred <- "A3301"
    }else{
      hla_prot <- strsplit(x=strsplit(hla_allele,split = '\\*')[[1]][1],split = '-')[[1]][2]
      hla_allele_mix_mhc_pred <- paste0(hla_prot,gsub(pattern = ':',replacement = '',x = strsplit(hla_allele,split = '\\*')[[1]][2]))
    }
    system(glue::glue("/home/zmxu/G2G-HBV/MixMHCpred/MixMHCpred -i {mut_file}_{residue}_k_{nchar(cur_mut[1])}.fasta -o {mut_file}_{residue}_k_{nchar(cur_mut[1])}.MixMHCpred.out -a {hla_allele_mix_mhc_pred}"))
  }
  #Write non-mut file
  for(i in 1:length(k_mer_list)){
    cur_non_mut <- k_mer_list[[i]]$non_mut
    for(j in 1:length(cur_non_mut)){
      cur_non_mut_spec_res <- cur_non_mut[[j]]
      header <- sapply(1:nchar(cur_non_mut_spec_res[1]),function(x) paste0('pos_',x))
      residue <- strsplit(x = cur_non_mut_spec_res[1],split = '')[[1]][1]
      write(paste0(mapply(function(x,y) paste0('>',y,'\n',x),cur_non_mut_spec_res,header),collapse = '\n'),glue::glue("{ref_file}_{residue}_k_{nchar(cur_non_mut_spec_res[1])}.fasta"),append = F)
      lapply(pred_methods,function(x) system(glue::glue("python ~/G2G-HBV/mhc_i/src/predict_binding.py {x} {hla_allele} {nchar(cur_non_mut_spec_res[1])} {ref_file}_{residue}_k_{nchar(cur_non_mut_spec_res[1])}.fasta > {ref_file}_{residue}_k_{nchar(cur_non_mut_spec_res[1])}.{x}.out")))
      system(glue::glue("~/G2G-HBV/MixMHCpred/MixMHCpred -i {ref_file}_{residue}_k_{nchar(cur_non_mut_spec_res[1])}.fasta -o {ref_file}_{residue}_k_{nchar(cur_non_mut_spec_res[1])}.MixMHCpred.out -a {hla_allele_mix_mhc_pred}"))
    }
  }
}
GetBestBinders <- function(prefix,aa,epitope_res_dir ='~/G2G-HBV/HLA_epitopes/',mhc_mix = T){
  epitope_res_files <- dir(epitope_res_dir)
  
  best_binder_per_k <- lapply(8:12,function(k) data.table::fread(glue::glue("{epitope_res_dir}{prefix}{k}.netmhcpan_el.out")) %>% dplyr::select(peptide,rank) %>% dplyr::arrange(rank))
  argmin_k <- which.min(sapply(best_binder_per_k,function(x) min(x$rank)))
  best_binder_el <- best_binder_per_k[[argmin_k]][1,] %>% dplyr::select(peptide = peptide,rank = rank)
  best_binder_peptide <- best_binder_el$peptide

  best_binder_ba <- data.table::fread(glue::glue("{epitope_res_dir}{prefix}{argmin_k + 7}.netmhcpan_ba.out"))
  best_binder_ba <- best_binder_ba[which(best_binder_ba$peptide == best_binder_peptide),] %>% dplyr::select(peptide,ic50,seq_num)
  best_binder_seq_num <- best_binder_ba$seq_num[which(best_binder_ba$peptide == best_binder_peptide)]
  best_binder_ba <- best_binder_ba %>% dplyr::select(-seq_num) 
  
  files_ba <- epitope_res_files[sapply(epitope_res_files, function(x) grepl(glue::glue("k_{argmin_k + 7}"),x) & grepl(glue::glue("_ba.out"),x) & grepl(glue::glue("ref"),x) & grepl(glue::glue("{aa}"),x))]
  ba_res <- lapply(files_ba,function(x) data.table::fread(paste0(epitope_res_dir,x)) %>% dplyr::filter(seq_num == as.character(best_binder_seq_num)) %>% dplyr::select(peptide,ic50))
  files_el <- epitope_res_files[sapply(epitope_res_files, function(x) grepl(glue::glue("k_{argmin_k + 7}"),x) & grepl(glue::glue("_el.out"),x) & grepl(glue::glue("ref"),x) & grepl(glue::glue("{aa}"),x))]
  el_res <- lapply(files_el,function(x) data.table::fread(paste0(epitope_res_dir,x)) %>% dplyr::filter(seq_num == as.character(best_binder_seq_num)) %>% dplyr::select(peptide,rank))
  
  if(mhc_mix){
    best_binder_mhcmix <- data.table::fread(glue::glue("{epitope_res_dir}{prefix}{argmin_k + 7}.MixMHCpred.out"),skip = 'Peptide') 
    best_binder_mhcmix <- best_binder_mhcmix[which(best_binder_mhcmix$Peptide == best_binder_peptide),] %>% dplyr::select(peptide = Peptide,perc_rank_MixMHCPred = `%Rank_bestAllele`)
    files_mhcmix <- epitope_res_files[sapply(epitope_res_files, function(x) grepl(glue::glue("k_{argmin_k + 7}"),x) & grepl(glue::glue("MixMHCpred.out"),x) & grepl(glue::glue("ref"),x) & grepl(glue::glue("{aa}"),x))]
    mhcmix_res <- lapply(files_mhcmix,function(x) data.table::fread(paste0(epitope_res_dir,x),skip = 'Peptide')[best_binder_seq_num,]  %>% dplyr::select(peptide = Peptide,perc_rank_MixMHCPred = `%Rank_bestAllele`))
    
    res <- rbind(dplyr::inner_join(best_binder_mhcmix,dplyr::inner_join(best_binder_ba,best_binder_el)),dplyr::inner_join(do.call(rbind,mhcmix_res),dplyr::inner_join(do.call(rbind,ba_res),do.call(rbind,el_res))))
    res$AA_res <- sapply(res$peptide,function(x) strsplit('',x=x)[[1]][best_binder_seq_num])
    
  }else{
    res <- rbind(dplyr::inner_join(best_binder_ba,best_binder_el),dplyr::inner_join(do.call(rbind,ba_res),do.call(rbind,el_res)))
    res$AA_res <- sapply(res$peptide,function(x) strsplit('',x=x)[[1]][best_binder_seq_num])
  }
  return(res)
}
GetEpitopeVsAlleleScores <- function(out_dir,epitopes){
  hla_alleles <- dir(out_dir)
  ba_epitope_allele_mat <- matrix(NA,nrow = length(hla_alleles),ncol = length(epitopes))
  rownames(ba_epitope_allele_mat) <- hla_alleles
  colnames(ba_epitope_allele_mat) <- epitopes
  el_epitope_allele_mat <- matrix(NA,nrow = length(hla_alleles),ncol = length(epitopes))
  rownames(el_epitope_allele_mat) <- hla_alleles
  colnames(el_epitope_allele_mat) <- epitopes
  
  for(i in 1:length(hla_alleles)){
    cur_hla_allele <- hla_alleles[i]
    cur_files <- dir(glue::glue("{out_dir}/{cur_hla_allele}"))

    cur_ba <- lapply(cur_files[grepl('ba',cur_files)],function(x) data.table::fread(glue::glue("{out_dir}/{cur_hla_allele}/{x}")))
    if(any(sapply(cur_ba,nrow) == 0)){
      next
    }
    cur_ba <- do.call(rbind,cur_ba) %>% dplyr::select(ic50,peptide) %>% dplyr::distinct()
    ba_epitope_allele_mat[i,] <- cur_ba$ic50[match(epitopes,cur_ba$peptide)]


    cur_el <- lapply(cur_files[grepl('el',cur_files)],function(x) data.table::fread(glue::glue("{out_dir}/{cur_hla_allele}/{x}")))
    cur_el <- do.call(rbind,cur_el) %>% dplyr::select(rank,peptide) %>% dplyr::distinct()
    el_epitope_allele_mat[i,] <- cur_el$rank[match(epitopes,cur_el$peptide)]
  }
  el_epitope_allele_mat <- el_epitope_allele_mat[apply(el_epitope_allele_mat,1,function(x) !any(is.na(x))),]
  ba_epitope_allele_mat <- ba_epitope_allele_mat[apply(ba_epitope_allele_mat,1,function(x) !any(is.na(x))),]

  return(list(ba_epitope_allele_mat=ba_epitope_allele_mat,el_epitope_allele_mat=el_epitope_allele_mat))
}

GetBestBindersAcrossAlleles <- function(AA_variant,best_kmers,all_kmers,hla_allele,hla_dosage){
  #Get All HLA Alleles in same gene
    hla_gene <- strsplit(hla_allele,split = '_')[[1]][1]
    hla_dosage <- hla_dosage[,grepl(pattern = paste0(hla_gene,'_'),colnames(hla_dosage))]
    all_alleles_raw <- colnames(hla_dosage)
    all_alleles <- sapply(all_alleles_raw,function(x) {
      allele <- strsplit(x=x,split = '_')[[1]]
      return(paste0('HLA-',allele[1],'*',gsub("[^0-9.]","",allele[2]),':',gsub("[^0-9.]","",allele[3])))
  })
  #Run Peptide prediction for peptide bound by G2G associated allele
  # lapply(1:length(ctl_alleles_raw),function(i){
  #   RunPredictBinding(out_dir = glue::glue('~/G2G-HBV/HLA_epitopes/best_binder_vs_alleles/{AA_variant}/{all_alleles_raw[i]}/'),all_alleles_raw[i],best_kmers,all_alleles[i])
  # })

  #Run peptide scan across all peptides (best in sliding window)
  # lapply(1:length(ctl_alleles_raw),function(i){
  #   RunPredictBinding(out_dir = glue::glue('~/G2G-HBV/HLA_epitopes/all_alleles/{AA_variant}/{all_alleles_raw[i]}/'),AA_variant,all_kmers,all_alleles[i])
  # })

  #Parse Peptide prediction
  allele_scores <- GetEpitopeVsAlleleScores(glue::glue('~/G2G-HBV/HLA_epitopes/best_binder_vs_alleles/{AA_variant}'),unlist(best_kmers))

  variant_allele <- strsplit(AA_variant,split = '_')[[1]][length(strsplit(AA_variant,split = '_')[[1]])]
  alleles_scores_sliding_window <- lapply(dir(glue::glue('~/G2G-HBV/HLA_epitopes/all_alleles/{AA_variant}')),function(x) {
    tryCatch(GetBestBinders(prefix = glue::glue("mut_{AA_variant}_{variant_allele}_k_"),aa = AA_variant,epitope_res_dir = glue::glue("~/G2G-HBV/HLA_epitopes/all_alleles/{AA_variant}/{x}/"),mhc_mix = F),error = function(e) return(NA))}
  )
  names(alleles_scores_sliding_window) <- dir(glue::glue('~/G2G-HBV/HLA_epitopes/all_alleles/{AA_variant}'))
  alleles_scores_sliding_window <- alleles_scores_sliding_window[sapply(alleles_scores_sliding_window,function(x) !all(is.na(x)))]

  el_epitope_allele_mat_sliding <- matrix(data = NA,nrow = length(alleles_scores_sliding_window),ncol = nrow(alleles_scores_sliding_window[[1]]))
  colnames(el_epitope_allele_mat_sliding) <- alleles_scores_sliding_window[[1]]$AA_res
  rownames(el_epitope_allele_mat_sliding) <- names(alleles_scores_sliding_window)
  ba_epitope_allele_mat_sliding <- el_epitope_allele_mat_sliding

  for(i in 1:nrow(el_epitope_allele_mat_sliding)){
    el_epitope_allele_mat_sliding[i,alleles_scores_sliding_window[[i]]$AA_res] <- alleles_scores_sliding_window[[i]]$rank
    ba_epitope_allele_mat_sliding[i,alleles_scores_sliding_window[[i]]$AA_res] <- alleles_scores_sliding_window[[i]]$ic50
  }
  
  peptides_sliding <- data.frame(HLA_Allele = character(),AA = character(),Peptide = character())
  for(i in 1:length(alleles_scores_sliding_window)){
    peptides_sliding <- rbind(peptides_sliding,data.frame(HLA_Allele = rep(names(alleles_scores_sliding_window)[i],nrow(alleles_scores_sliding_window[[i]])),
                                                          AA = alleles_scores_sliding_window[[i]]$AA_res,
                                                          Peptide = alleles_scores_sliding_window[[i]]$peptide))
  }
  
  return(list(allele_scores_el=allele_scores$el_epitope_allele_mat,allele_scores_ba=allele_scores$ba_epitope_allele_mat,allele_scores_el_sliding = el_epitope_allele_mat_sliding,allele_scores_ba_sliding = ba_epitope_allele_mat_sliding,peptides_sliding=peptides_sliding))
}


GetFreqByAllele <- function(hla_dosage,haplotype_df,hla_gene,target_allele,dat_env,binding_affinity,peptides_sliding,freq_cutoff=0.2,second_digit = F,show_other =F){
  alleles_to_keep <- which(grepl(pattern = paste0(hla_gene,'_'),colnames(hla_dosage)))
  hla_dosage_mat <- as.matrix(hla_dosage[,..alleles_to_keep])
  #Cleanup names
  colnames(hla_dosage_mat) <- sapply(colnames(hla_dosage_mat),function(x) gsub(gsub(x,pattern = 'Q',replacement = ''),pattern = 'N',replacement = ''))
  hla_dosage_id <- hla_dosage$IID
  hla_dosage_gilead_id <- dat_env$dat_ids$pathogen_id[match(hla_dosage_id,dat_env$dat_ids$host_id)]
  
  #Get HLA allele frequencies
  frq <- apply(hla_dosage_mat,2,function(x) sum(x!=0,na.rm = T)/sum(!is.na(x)))
  
  target_allele_mat <- list(hla_dosage_mat[,target_allele,drop=F])
  target_allele_name <- strsplit(target_allele,split = '_')[[1]] 
  names(target_allele_mat) <- paste0('G2G Associated:\n HLA-',target_allele_name[1],'*',target_allele_name[2],':',target_allele_name[3],' (',round(frq[which(names(frq) == target_allele)]*100,1),'%)')
  
  
  if(second_digit){
    hla_second_digit <- sapply(colnames(hla_dosage_mat),function(x) paste0('Other:\n HLA-',paste0(strsplit(x=x,split = '_')[[1]][1:2],collapse = '*')))
    uniq_digits <- unique(hla_second_digit)
    uniq_digits_freq <- sapply(uniq_digits,function(x) {
      matching_mat <- hla_dosage_mat[,hla_second_digit == x & colnames(hla_dosage_mat) != target_allele,drop=F]
      frq <- apply(matching_mat,2,function(x) sum(x!=0,na.rm = T)/sum(!is.na(x)))
      return(sum(frq))
    })
    
    digits_to_query <- as.list(uniq_digits[uniq_digits_freq > freq_cutoff])
    alleles_to_query_mat <- lapply(digits_to_query,function(x) hla_dosage_mat[,hla_second_digit %in% x & colnames(hla_dosage_mat) != target_allele,drop=F])
    names(alleles_to_query_mat) <- sapply(digits_to_query,function(x) paste0(x,collapse = '_'))
    
    alleles_to_query_mat <- c(target_allele_mat,alleles_to_query_mat)
    
  }else{
    hla_alleles <- sapply(1:ncol(hla_dosage_mat),function(x) paste0('HLA-',paste0(paste0(strsplit(x=colnames(hla_dosage_mat)[x],split = '_')[[1]][1:2],collapse = '*'),':',
                                                                                             strsplit(x=colnames(hla_dosage_mat)[x],split = '_')[[1]][3]),' (',round(frq[x]*100,1),'%)'))
    alleles_to_query <- as.list(setdiff(names(frq)[frq > freq_cutoff],target_allele))
    alleles_to_query_mat <- lapply(alleles_to_query,function(x) hla_dosage_mat[,x,drop=F])
    names(alleles_to_query_mat) <- setdiff(hla_alleles[frq > freq_cutoff],gsub(names(target_allele_mat),pattern = 'G2G Associated:\n ',replacement = ''))
    alleles_to_query_mat <- c(target_allele_mat,alleles_to_query_mat)
    
    if(show_other){
      other_alleles <- setdiff(names(frq),c(target_allele,alleles_to_query))
      other_alleles_mat <- list(hla_dosage_mat[,other_alleles,drop=F])
      names(other_alleles_mat) <- 'Other'
      alleles_to_query_mat <- c(alleles_to_query_mat,other_alleles_mat)
    }
  }

  freq_df <- data.frame(hla_group = character(),AA = character(),Frq = numeric(),Mixed = logical())
  
  uniq_AA <- unique(haplotype_df$AA)
  
  for(i in 1:length(alleles_to_query_mat)){
    cur_hla_allele <- names(alleles_to_query_mat)[i]
    cur_id <- hla_dosage_gilead_id[apply(alleles_to_query_mat[[i]],1,function(x) any(x!=0,na.rm = T))]
    cur_haplo_info <- dplyr::filter(haplotype_df,pathogen_id %in% cur_id)
    
    for(cur_AA in uniq_AA){
      dup_id <- cur_haplo_info$pathogen_id[duplicated(cur_haplo_info$pathogen_id)]
      #Get samples without mixed epitopes
      
      cur_AA_not_mixed <- nrow(dplyr::filter(cur_haplo_info,!pathogen_id %in% dup_id & AA == cur_AA))
      freq_df <- rbind(freq_df,data.frame(hla_group = cur_hla_allele,AA = cur_AA,Frq = cur_AA_not_mixed / length(unique(cur_haplo_info$pathogen_id)),Mixed=F))
      #Get samples with mixed epitopes
      cur_AA_mixed <- nrow(dplyr::filter(cur_haplo_info,pathogen_id %in% dup_id & AA == cur_AA))
      freq_df <- rbind(freq_df,data.frame(hla_group = cur_hla_allele,AA = cur_AA,Frq = cur_AA_mixed / length(unique(cur_haplo_info$pathogen_id)),Mixed=T))
    }
  }
  #Get weighted mean binding score
  allele_members <- lapply(alleles_to_query_mat,function(x) colnames(x))
  allele_freq <- lapply(alleles_to_query_mat,function(x) apply(x,2,function(x) sum(x!=0,na.rm = T)/sum(!is.na(x))))
  
  binding_affinity_by_hla_group <- data.frame(hla_group = character(),AA = character(),binding_score = numeric())
  for(i in 1:length(allele_freq)){
    cur_freq <- allele_freq[[i]]
    cur_common_alleles <- intersect(rownames(binding_affinity),names(cur_freq))
    
    cur_binding_affinity <- binding_affinity[cur_common_alleles,,drop=F]
    cur_freq <- cur_freq[cur_common_alleles]
    
    for(k in 1:ncol(cur_binding_affinity)){
      binding_affinity_by_hla_group <- rbind(binding_affinity_by_hla_group,data.frame(hla_group = names(allele_freq)[i],
                                                                                      AA = colnames(cur_binding_affinity)[k],
                                                                                      binding_score = as.vector(cur_freq %*% cur_binding_affinity[,k,drop=F] / sum(cur_freq))))
      
      
    }
    
  }
  
  freq_df <- dplyr::left_join(freq_df,binding_affinity_by_hla_group,by=c('hla_group'='hla_group','AA'='AA'))
  freq_df <- dplyr::mutate(freq_df,strong_binder = ifelse(round(binding_score,1) >= 1,F,T))
  
  freq_df$hla_group_raw <- sapply(freq_df$hla_group,function(x) strsplit(x=x,split = 'HLA-')[[1]][2])
  freq_df$hla_group_raw <- sapply(freq_df$hla_group_raw,function(x) gsub(x=gsub(x=x,pattern = '\\*',replacement = '_'),pattern = ':',replacement = '_'))
  freq_df <- dplyr::left_join(freq_df,peptides_sliding,by=c('hla_group_raw'='HLA_Allele','AA'='AA')) %>% dplyr::select(-hla_group_raw)
  
  peptides_sliding$second_digit <- sapply(peptides_sliding$HLA_Allele,function(x) paste0(strsplit(x=x,split = '_')[[1]][1:2],collapse = '_'))
  #Add peptides for non-target alleles
  for(i in 1:nrow(freq_df)){
    if(is.na(freq_df$Peptide[i])){
      cur_hla_group <- gsub(gsub(x=strsplit(strsplit(freq_df$hla_group[i],split = 'HLA-')[[1]][2],split = ' \\(')[[1]][1],pattern = '\\*',replacement = '_'),pattern = ':',replacement = '_')
      if(second_digit){
        cur_peptides <- dplyr::filter(peptides_sliding,second_digit == cur_hla_group & AA == freq_df$AA[i])
      }else{
        cur_peptides <- dplyr::filter(peptides_sliding,HLA_Allele == cur_hla_group & AA == freq_df$AA[i])
      }
      freq_df$Peptide[i] <- paste0(unique(cur_peptides$Peptide),collapse = '\n')
    }
  }
  
  return(freq_df)
}


#Load Asian Data
dat_env_asn <- new.env()
load('~/G2G-HBV/data/results_POP_asian_GT_A_C_D_B_TYPE_pPCA_LOCO_FALSE/prep-data.rda',envir = dat_env_asn)
asn_cluster <- dat_env_asn$ids_unrelated %>% dplyr::select(host_id = X1)

dat_pathogen_raw_asn <- as.data.frame(dat_env_asn$dat_pathogen_raw[match(asn_cluster$host_id,dat_env_asn$dat_pathogen_raw$ID),-1])
rownames(dat_pathogen_raw_asn) <- asn_cluster$host_id
all_AA <- colnames(dat_pathogen_raw_asn[,apply(dat_pathogen_raw_asn,2,function(x) min(table(x)[!is.na(table(x))]) >= 1)])
all_AA_prot <- sapply(all_AA,function(x) strsplit(x=strsplit(x=x,split = 'gene_')[[1]][2],split = '_pos')[[1]][1])
all_AA_pos <- as.numeric(sapply(all_AA,function(x) strsplit(strsplit(x,split = 'pos_')[[1]][2],split = '_')[[1]][1]))
all_AA_residual <- sapply(all_AA,function(x) strsplit(strsplit(x,split = 'pos_')[[1]][2],split = '_')[[1]][2])
all_AA_freq <- apply(dat_pathogen_raw_asn,2,function(x){
  tbl <- table(x)
  if("1" %in% names(tbl)){
    return(tbl["1"] / sum(tbl))
  }else{
    return(0)
  }
})
all_AA_df_asn <- data.frame(Prot = all_AA_prot,Pos = all_AA_pos,Residual = all_AA_residual,Freq = all_AA_freq,stringsAsFactors = F)
#Load HLA_Allele Dosage
hla_dosage_asn <- data.table::fread('~/G2G-HBV/PyHLA_Out/Asian/hla_allele_assocAllele.dosage')


#Reference genome for Genotype C
hbv_geno_c_prot_seq <- seqinr::read.fasta('~/G2G-HBV/data/raw/gilead_20181126/viral_seq/Ref_Genome/GQ924620_AA.fasta')
names(hbv_geno_c_prot_seq) <- c('Pol','S','X','PC_C','C')

#PC_C_160, run through all k-mer lengths and possible residues.
k_mer_length <- 8:14
mut_pos = 160
k_mers_PC_C_pos_160_A <- lapply(k_mer_length,function(i) list(mut = ExtractKMers(i,mut_pos,hbv_geno_c_prot_seq$PC_C,'A'),non_mut = lapply(setdiff(dplyr::filter(all_AA_df_asn,Prot == 'PC_C',Pos == mut_pos)$Residual,'A'),function(x) ExtractKMers(i,mut_pos,hbv_geno_c_prot_seq$PC_C,x))))
RunPredictBinding('~/G2G-HBV/HLA_epitopes/','PC_C_pos_0160_A',k_mers_PC_C_pos_160_A,"HLA-A*33:03")

#Pol_49
mut_pos = 49
k_mers_Pol_pos_49_N <- lapply(k_mer_length,function(i) list(mut = ExtractKMers(i,mut_pos,hbv_geno_c_prot_seq$Pol,'N'),non_mut = lapply(setdiff(dplyr::filter(all_AA_df_asn,Prot == 'Pol',Pos == mut_pos)$Residual,'N'),function(x) ExtractKMers(i,mut_pos,hbv_geno_c_prot_seq$Pol,x))))
RunPredictBinding('~/G2G-HBV/HLA_epitopes/','Pol_pos_0049_N',k_mers_Pol_pos_49_N,"HLA-A*02:06")


#European Cohort
dat_env_eur <- new.env()
load('~/G2G-HBV/data/results_POP_european_GT_A_C_D_F_H_TYPE_pPCA_LOCO_FALSE/prep-data.rda',envir = dat_env_eur)
eur_cluster <- dat_env_eur$ids_unrelated %>% dplyr::select(host_id = X1)

dat_pathogen_raw_eur <- as.data.frame(dat_env_eur$dat_pathogen_raw[match(eur_cluster$host_id,dat_env_eur$dat_pathogen_raw$ID),-1])
rownames(dat_pathogen_raw_eur) <- eur_cluster$host_id

all_AA <- colnames(dat_pathogen_raw_eur[,apply(dat_pathogen_raw_eur,2,function(x) min(table(x)[!is.na(table(x))]) >= 1)])
all_AA_prot <- sapply(all_AA,function(x) strsplit(x=strsplit(x=x,split = 'gene_')[[1]][2],split = '_pos')[[1]][1])
all_AA_pos <- as.numeric(sapply(all_AA,function(x) strsplit(strsplit(x,split = 'pos_')[[1]][2],split = '_')[[1]][1]))
all_AA_residual <- sapply(all_AA,function(x) strsplit(strsplit(x,split = 'pos_')[[1]][2],split = '_')[[1]][2])
all_AA_freq <- apply(dat_pathogen_raw_eur,2,function(x){
  tbl <- table(x)
  if("1" %in% names(tbl)){
    return(tbl["1"] / sum(tbl))
  }else{
    return(0)
  }
})
all_AA_df_eur <- data.frame(Prot = all_AA_prot,Pos = all_AA_pos,Residual = all_AA_residual,Freq = all_AA_freq,stringsAsFactors = F)
#Load HLA_Allele Dosage
hla_dosage_eur <- data.table::fread('~/G2G-HBV/PyHLA_Out/European/hla_allele_assocAllele.dosage')

#Reference genome for Genotype D
hbv_geno_d_prot_seq <- seqinr::read.fasta('~/G2G-HBV/data/raw/gilead_20181126/viral_seq/Ref_Genome/FJ356716_AA.fasta')[c(1,2,5,6)]
names(hbv_geno_d_prot_seq) <- c('Pol','S','X','PC_C')


mut_pos = 67
k_mers_PC_C_pos_67_Y <- lapply(k_mer_length,function(i) list(mut = ExtractKMers(i,mut_pos,hbv_geno_d_prot_seq$PC_C,'Y'),non_mut = lapply(setdiff(dplyr::filter(all_AA_df_eur,Prot == 'PC_C',Pos == mut_pos)$Residual,'Y'),function(x) ExtractKMers(i,mut_pos,hbv_geno_c_prot_seq$PC_C,x))))
RunPredictBinding('~/G2G-HBV/HLA_epitopes/','PC_C_pos_0067_Y',k_mers_PC_C_pos_67_Y,"HLA-A*01:01")