GetPopConsensusSeq <- function(Haplotypes,start_offset=33){
  sample_consensus <- Haplotypes %>% dplyr::group_by(ID) %>%
    dplyr::arrange(desc(prop_reads)) %>% dplyr::filter(row_number() == 1)
  #Get all haplotypes
  all_haplo <- sapply(sample_consensus$haplotype,function(x) paste0(strsplit(x,split = '')[[1]][(start_offset+1):(nchar(x)-1)],collapse = ''))
  #Align all haplotypes
  pop_aln <- msa(DNAStringSet(all_haplo))
  pop_consensus <- tolower(msaConsensusSequence(pop_aln,type = 'upperlower'))
  return(pop_consensus)
}

#From https://biostatmatt.com/archives/2902
round_preserve_sum <- function(x, digits = 0) {
  up <- 10 ^ digits
  x <- x * up
  y <- floor(x)
  indices <- tail(order(x-y), round(sum(x)) - sum(y))
  y[indices] <- y[indices] + 1
  y / up
}

CalculateDnDs <- function(Haplotypes,Ref_Seq = NULL,start_offset = 33,out_dir = '/home/zmxu/G2G-HBV/DnDs/',within_host = T){
  cur_wd <- getwd()
  uniq_IDs <- unique(Haplotypes$ID)
  res_df <- data.frame(ID = character(),dNN = numeric(),dSN = numeric(),dSS = numeric(),dNS = numeric())
  pure_samples <- 0
  raw_results <- vector(mode = 'list',length = length(uniq_IDs))
  names(raw_results) <- uniq_IDs
  for(i in 1:length(uniq_IDs)){
    cur_ID <- uniq_IDs[i]
    system(glue::glue('mkdir -p {out_dir}/{cur_ID}'))
    setwd(glue::glue('{out_dir}/{cur_ID}'))
    cur_Haplotypes <- dplyr::filter(Haplotypes,ID == cur_ID)
    #Extract MyrB from second start codon
    cur_Haplotypes$haplotype <- sapply(cur_Haplotypes$haplotype,function(x) paste0(strsplit(x,split = '')[[1]][(start_offset+1):(nchar(x)-1)],collapse = ''))
    #Collapse identical haplotypes
    cur_Haplotypes <- cur_Haplotypes %>% dplyr::group_by(haplotype) %>% dplyr::summarise(prop_reads=sum(prop_reads))
    if(nrow(cur_Haplotypes) <= 1){
      res_df <- rbind(res_df,data.frame(ID = cur_ID,dNN=NA,dSN=NA,dSS=NA,dNS=NA,dN=NA,dS=NA,stringsAsFactors = F))
      pure_samples <- pure_samples+1
    }else{
      if(!within_host){
        #Get Consensus within sample
        intra_host_consensus <- tolower(msaConsensusSequence(msa(DNAStringSet(rep(cur_Haplotypes$haplotype,round(10*cur_Haplotypes$prop_reads)))),type = 'upperlower'))
        
        for(k in 1:nrow(cur_Haplotypes)){
          if(is.null(Ref_Seq)){
            cur_seq <- list(intra_host_consensus,cur_Haplotypes$haplotype[k])
          }else{
            cur_seq <- list(Ref_Seq,cur_Haplotypes$haplotype[k])
          }
          write.fasta(sequences = cur_seq,names = c('Ref',paste0('Haplo',k)),file.out = glue::glue("{out_dir}/{cur_ID}/haplo{k}.fasta"))
          system(glue::glue('~/Software/OLGenie/OLGenie.pl --fasta_file=haplo{k}.fasta --frame=ss12 --output_file=haplo{k}.txt'))
        }
        dnds <- lapply(1:nrow(cur_Haplotypes),function(k) data.table::fread(glue::glue("{out_dir}/{cur_ID}/haplo{k}.txt")))
        dNN <-  sum(as.vector(cur_Haplotypes$prop_reads %*% do.call(rbind,lapply(dnds,function(x)x$NN_diffs))) )/ sum(dnds[[1]]$NN_sites)
        dSN <-  sum(as.vector(cur_Haplotypes$prop_reads %*% do.call(rbind,lapply(dnds,function(x)x$SN_diffs)))) / sum(dnds[[1]]$SN_sites)
        dSS <-  sum(as.vector(cur_Haplotypes$prop_reads %*% do.call(rbind,lapply(dnds,function(x)x$SS_diffs)))) / sum(dnds[[1]]$SS_sites)
        dNS <- sum(as.vector(cur_Haplotypes$prop_reads %*% do.call(rbind,lapply(dnds,function(x)x$NS_diffs)))) / sum(dnds[[1]]$NS_sites)
        res_df <- rbind(res_df,data.frame(ID = cur_ID,dNN=dNN,dSN=dSN,dSS=dSS,dNS=dNS,dN=dN,dS=dS,stringsAsFactors = F))
        raw_results[[i]] <- dnds
      }else{
        #Construct pseudo alignment, by up-sampling haplotypes to 1000 sequences
        up_sample_haplo <- rep(cur_Haplotypes$haplotype,round_preserve_sum(1000*cur_Haplotypes$prop_reads,0))
        haplo_names <- unlist(sapply(1:nrow(cur_Haplotypes),function(x) paste0('Haplo',x,'_',1:(round_preserve_sum(1000*cur_Haplotypes$prop_reads,0)[x]))))
        write.fasta(sequences = lapply(up_sample_haplo,function(x) strsplit(x=toupper(x),split = '')[[1]]),names = haplo_names,file.out = glue::glue("{out_dir}/{cur_ID}/pseudo_aln.fasta"))
        system(glue::glue('~/Software/OLGenie/OLGenie.pl --fasta_file=pseudo_aln.fasta --frame=ss12 --output_file=dnds.txt'))
        
        dnds <- data.table::fread(glue::glue("{out_dir}/{cur_ID}/dnds.txt"))
        dNN <-  sum(dnds$NN_diffs) / sum(dnds$NN_sites)
        dSN <-  sum(dnds$SN_diffs) / sum(dnds$SN_sites)
        dSS <-  sum(dnds$SS_diffs) / sum(dnds$SS_sites) 
        dNS <- sum(dnds$NS_diffs) / sum(dnds$NS_sites)  
        dN <- sum(dnds$NN_diffs + dnds$NS_diffs) / sum(dnds$NN_sites + dnds$NS_sites)
        dS <- sum(dnds$SS_diffs + dnds$SN_diffs) / sum(dnds$SS_sites + dnds$SN_sites)
        
        res_df <- rbind(res_df,data.frame(ID = cur_ID,dNN=dNN,dSN=dSN,dSS=dSS,dNS=dNS,dN=dN,dS=dS,stringsAsFactors = F))
        raw_results[[i]] <- dnds
      }
    }
  }
  raw_results_filt <- raw_results[sapply(raw_results,function(x) !is.null(x))]
  per_site_dNN <- do.call(rbind,lapply(raw_results_filt,function(x) x$NN_diffs / x$NN_sites))
  per_site_dSN <- do.call(rbind,lapply(raw_results_filt,function(x) x$SN_diffs / x$SN_sites))
  per_site_dSS <- do.call(rbind,lapply(raw_results_filt,function(x) x$SS_diffs / x$SS_sites))
  per_site_dNS <- do.call(rbind,lapply(raw_results_filt,function(x) x$NS_diffs / x$NS_sites))
  setwd(cur_wd)
  return(list(res_df=res_df,raw_results=raw_results,per_site = list(per_site_dNN=per_site_dNN,per_site_dSN=per_site_dSN,per_site_dSS=per_site_dSS,per_site_dNS=per_site_dNS,pure_samples=pure_samples)))
}
Geno_C_Ref <- read.fasta('~/G2G-HBV/data/raw/gilead_20181126/viral_seq/Ref_Genome/GQ924620_NUC.fasta')
# Geno_C_Ref_Myr_B <- paste0(Geno_C_Ref$GQ924620.1[2881:3024],collapse = '')

# geno_C_consensus <- GetPopConsensusSeq(rbind(prop_geno_C$results_shotgun_bind,prop_WT_geno_C$results_shotgun_bind))
geno_C_dnds_raw <- CalculateDnDs(Haplotypes = prop_geno_C$results_shotgun_bind,Ref_Seq = NULL,out_dir = '/home/zmxu/G2G-HBV/DnDs/')
geno_C_WT_dnds_raw <- CalculateDnDs(Haplotypes = prop_WT_geno_C$results_shotgun_bind,Ref_Seq = NULL,out_dir ='/home/zmxu/G2G-HBV/DnDs/') 

geno_C_dnds <-  geno_C_dnds_raw$res_df %>% dplyr::left_join(dat_ids_pop,by=c('ID'='pathogen_id')) %>% dplyr::left_join(dat_alternative_outcome) %>% dplyr::left_join(dat_covars_raw %>% dplyr::select(host_id,AVG_COVERAGE))
geno_C_WT_dnds <-  geno_C_WT_dnds_raw$res_df %>% dplyr::left_join(dat_ids_pop,by=c('ID'='pathogen_id')) %>% dplyr::left_join(dat_alternative_outcome) %>% dplyr::left_join(dat_covars_raw %>% dplyr::select(host_id,AVG_COVERAGE))

# Geno_B_Ref <- read.fasta('~/G2G-HBV/data/raw/gilead_20181126/viral_seq/Ref_Genome/AB219428_NUC.fasta')
# Geno_B_Ref_Myr_B <- paste0(Geno_B_Ref$AB219428.1[2881:3024],collapse = '')

# geno_B_consensus <- GetPopConsensusSeq(rbind(prop_geno_B$results_shotgun_bind,prop_WT_geno_B$results_shotgun_bind))
geno_B_dnds_raw <- CalculateDnDs(Haplotypes = prop_geno_B$results_shotgun_bind,Ref_Seq = NULL,out_dir = '/home/zmxu/G2G-HBV/DnDs/')
geno_B_WT_dnds_raw <- CalculateDnDs(Haplotypes = prop_WT_geno_B$results_shotgun_bind,Ref_Seq = NULL,out_dir = '/home/zmxu/G2G-HBV/DnDs/') 

geno_B_dnds <-  geno_B_dnds_raw$res_df %>% dplyr::left_join(dat_ids_pop,by=c('ID'='pathogen_id')) %>% dplyr::left_join(dat_alternative_outcome) %>% dplyr::left_join(dat_covars_raw %>% dplyr::select(host_id,AVG_COVERAGE))
geno_B_WT_dnds <-  geno_B_WT_dnds_raw$res_df %>% dplyr::left_join(dat_ids_pop,by=c('ID'='pathogen_id')) %>% dplyr::left_join(dat_alternative_outcome) %>% dplyr::left_join(dat_covars_raw %>% dplyr::select(host_id,AVG_COVERAGE))

WT_dnds <- rbind(geno_C_WT_dnds,geno_B_WT_dnds)
S267F_dnds <- rbind(geno_C_dnds,geno_B_dnds)
merged_dnds <- rbind(data.frame(WT_dnds,rs2296651 = 'GG'),data.frame(S267F_dnds,rs2296651 = 'GA'))
