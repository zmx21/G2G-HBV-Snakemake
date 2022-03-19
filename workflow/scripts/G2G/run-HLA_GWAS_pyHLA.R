library(pbmcapply)
#Parse omnibus results given omnibus directory
ParseOmnibus <- function(cond_files_AA,parse_cond = T){
  all_files <- dir(cond_files_AA)
  chap_files <- sapply(all_files,function(x) grepl('.chap',x))
  chap_files <- all_files[chap_files]
  cond_files <- chap_files[sapply(chap_files,function(x) grepl('cond',x))]
  cond_files <- chap_files[sapply(chap_files,function(x) grepl('cond',x))]
  uncond_files <- setdiff(chap_files,cond_files)
  
  uncond_chap_results <- lapply(uncond_files,function(x) {
    raw_str <- system(glue::glue("grep 'p =' {cond_files_AA}{x}"),intern = T)
    if(is.null(unlist(strsplit(raw_str,split = 'p = ')))){
      return(NA)
    }else{
      return(unlist(strsplit(raw_str,split = 'p = ')))
    }})
  uncond_chap_p <- unlist(sapply(uncond_chap_results,function(x) as.numeric(x[2])))
  uncond_ids <- sapply(uncond_files, function(x) strsplit(x=x,split = '\\.')[[1]][1])
  uncond_chrom <- sapply(uncond_ids,function(x) strsplit(x=x,split = '_')[[1]][1])
  
  if(parse_cond){
    cond_chap_results <- lapply(cond_files,function(x) {
      raw_str <- system(glue::glue("grep 'p =' {cond_files_AA}{x}"),intern = T)
      if(is.null(unlist(strsplit(raw_str,split = 'p = ')))){
        return(NA)
      }else{
        return(unlist(strsplit(raw_str,split = 'p = ')))
      }})
    cond_chap_p <- unlist(sapply(cond_chap_results,function(x) as.numeric(x[2])))
    
    #Get IDs and CHROM
    cond_ids <- sapply(cond_files, function(x) strsplit(x=x,split = '\\.')[[1]][1])
    cond_chrom <- sapply(cond_ids,function(x) strsplit(x=x,split = '_')[[1]][1])
    
    aa_results <- list(cond = data.frame(`#CHROM` = cond_chrom,ID = cond_ids,P = cond_chap_p,stringsAsFactors=F),uncond = data.frame(`#CHROM` = uncond_chrom,ID = uncond_ids, P = uncond_chap_p,stringsAsFactors=F))
    colnames(aa_results$cond) <- c('#CHROM','ID','P')
    colnames(aa_results$uncond) <- c('#CHROM','ID','P')
    return(aa_results)
  }
  aa_results <- list(uncond = data.frame(`#CHROM` = uncond_chrom,ID = uncond_ids, P = uncond_chap_p,stringsAsFactors=F))
  colnames(aa_results$uncond) <- c('#CHROM','ID','P')
  return(aa_results)
}

ParseAllele <- function(input_pattern,allele_list,digit=NULL){
  if(is.na(input_pattern)){
    return(NA)
  }
  #If allele exists in input pattern
  if(input_pattern %in% allele_list$Allele){
    return(input_pattern)
  }#If doesn't exist, strip last letter
  else{
    #Find locus (before *)
    hla_locus <- strsplit(input_pattern,split = '\\*')[[1]][1]
    #Get remaining part of string
    hla_non_locus <- strsplit(input_pattern,split = '\\*')[[1]][2]
    hla_digits <- strsplit(hla_non_locus,split = ':')[[1]]
    #Remove all letters from digits
    hla_digits <- sapply(hla_digits,function(x) gsub("[^0-9.-]", "", x))
    #Return parsed hla allele, if digit specified up to specified digit
    if(is.null(digit)){
      parsed_alelle <- paste0(hla_locus,'*',paste0(hla_digits,collapse = ':'))
      return(parsed_alelle)
    }else{
      if(!digit %in% c(2,4,6)){
        stop('Digit not in multiples of 2 or exceeded 6')
      }
      parsed_alelle <- paste0(hla_locus,'*',paste0(hla_digits[1:(digit/2)],collapse = ':'))
      return(parsed_alelle)
      
    }
    
  }
}
Parse_HLA_LA_Output <- function(path,samples_to_include = NULL,pheno,allele_list,QC = T){
  result_file <- data.table::fread(path)
  #Do QC based on Q1 metric and Average coverage
  if(QC){
    #Set Low probability calls as NA (Q1 < 0.7) 
    result_file$Allele[which(result_file$Q1 < 0.7)] <- NA
  }
  
  #Filter out samples to include if specified 
  if(!is.null(samples_to_include)){
    result_file <- dplyr::filter(result_file,sample_id %in% samples_to_include)
  }
  #Convert to format which pyhla accepts
  n_col <- 2*length(unique(result_file$Locus)) + 2 #2 for id and pheno, rest for hla alleles (2 alleles at each locus)
  col_names <- c('ID','PHENO',sort(c(paste0(unique(result_file$Locus),'_1'),paste0(unique(result_file$Locus),'_2'))))
  
  #Create data frame in pyhla format
  output <- data.frame(matrix(ncol = n_col, nrow = length(samples_to_include)))
  colnames(output) <- col_names
  output$ID <- samples_to_include
  output$PHENO <- as.vector(t(pheno))
  #2 for cases, 1 for controls. Set REF AA (1) as Control, and Non-Ref AA (0) as Cases
  output$PHENO[output$PHENO==0] <- 2
  
  for(i in 3:ncol(output)){
    cur_locus <- colnames(output)[i]
    cur_hla_locus <- strsplit(cur_locus,split = '_')[[1]][1]
    cur_chr <- strsplit(cur_locus,split = '_')[[1]][2]
    cur_result <- dplyr::filter(result_file,Chromosome == cur_chr,Locus == cur_hla_locus)
    #If some samples with missing allele imputation
    if(nrow(cur_result) != nrow(output)){
      match_ind <- match(output$ID,cur_result$sample_id)
      match_ind_pseudo <- match_ind
      match_ind_pseudo[is.na(match_ind)] <- 1 #Construct pseudo row which will be fixed later
      cur_result_pseudo <- cur_result[match_ind_pseudo,]
      cur_result_pseudo[is.na(match_ind),] <- NA #Set pseudo row to NA
      cur_result_pseudo[is.na(match_ind),]$sample_id <- output$ID[is.na(match_ind)]
      cur_alleles <- sapply(cur_result_pseudo$Allele,function(x) ParseAllele(x,allele_list))
      output[,i] <- cur_alleles
    }else{
      cur_result <- cur_result[match(output$ID,cur_result$sample_id),] #Make sure ordering is same
      cur_alleles <- sapply(cur_result$Allele,function(x) ParseAllele(x,allele_list))
      output[,i] <- cur_alleles
    }
  }
  return(output)
}

ConvertPyHLAToDosage <- function(input_file_path,output_file_path){
  #File with sample and hla alelles. Mapped already to match IMGTHLA contained in PyHLA
  input_file <- data.table::fread(input_file_path)
  #Output file from PyHLA
  output_file <- data.table::fread(output_file_path,sep = ' ')
  #Process last column, which is seperated by tab
  n_col <- ncol(output_file)
  last_col <- output_file[,..n_col]
  alleles_for_each_residue <- sapply(as.vector(t(last_col)),function(x) strsplit(x=x,split = '\t')[[1]][2])
  names(alleles_for_each_residue) <- NULL
  alleles_for_each_residue <- lapply(alleles_for_each_residue,function(x) strsplit(x=x,split = ',')[[1]])
  aa_residues <- output_file$ID
  #Setup Dictionary which maps an single Allele to AA profile
  all_alleles <- unique(as.vector(t(input_file[,-c(1,2)])))
  #Get non-mapped alleles
  non_mapped_alleles <- system(glue::glue("grep 'no sequence is avaiable for: ' {output_file_path}.log.final.out"),intern = T)
  non_mapped_alleles <- sapply(non_mapped_alleles,function(x) gsub('no sequence is avaiable for: ','',x))
  dict <- vector(mode = 'list',length = length(all_alleles))
  names(dict) <- all_alleles
  for(i in 1:length(all_alleles)){
    cur_allele <- all_alleles[i]
    if(!cur_allele %in% non_mapped_alleles){
      #for each AA residue, check if cur allele has the residue
      dict[[i]] <- sapply(alleles_for_each_residue,function(x) ifelse(cur_allele %in% x,1,0))
    }else{
      cur_locus_allele <- strsplit(cur_allele,split = '*')[[1]][1]
      loci_aa <- sapply(aa_residues,function(x) strsplit(x,split = '_')[[1]][1])
      #Set all aa residues in same loci as NA, every other locus as 0
      dict[[i]] <- rep(NA,length(alleles_for_each_residue))
      dict[[i]][which(loci_aa!=cur_locus_allele)] <- 0
    }
  }
  #Get AA profile for each sample
  aa_dosage_matrix <- matrix(NA,nrow = nrow(input_file),ncol = length(aa_residues))
  rownames(aa_dosage_matrix) <- input_file$V1
  colnames(aa_dosage_matrix) <- aa_residues
  
  aa_geno_matrix <- matrix(NA,nrow = nrow(input_file),ncol = length(aa_residues))
  rownames(aa_geno_matrix) <- input_file$V1
  colnames(aa_geno_matrix) <- aa_residues
  
  locus_for_each_allele <- apply(input_file[,-c(1,2)],2,function(x){
    unique_alleles <- unique(sapply(x,function(y) strsplit(x=y,split = '\\*')[[1]][1]))
    return(unique_alleles[!is.na(unique_alleles)])
  })
  for(i in 1:nrow(aa_dosage_matrix)){
    #All alleles which sample has
    cur_allele_profile <- as.vector(t(input_file[i,-c(1,2)]))
    if(all(is.na(cur_allele_profile))){
      aa_dosage_matrix[i,] <- rep(NA,ncol(aa_dosage_matrix))
      aa_geno_matrix[i,] <- rep(".|.",ncol(aa_geno_matrix))
    }else{
      #Extract corresponding AA residue profile for each allele which sample has
      cur_dict_mapping <- dict[cur_allele_profile]
      #For uncalled residues, set AA in locus as NA
      cur_dict_mapping_na_added <- cur_dict_mapping
      if(any(is.na(cur_allele_profile))){
        for(cntr in which(is.na(cur_allele_profile))){
          cur_locus <- locus_for_each_allele[cntr]
          loci_aa <- sapply(aa_residues,function(x) strsplit(x,split = '_')[[1]][1])
          mapping <- rep(NA,length(loci_aa))
          mapping[loci_aa != cur_locus] <- 0
          cur_dict_mapping_na_added[[cntr]] <- mapping
        }
      }
      cur_aa_profile_mat <- do.call(rbind,cur_dict_mapping_na_added)
      #Get AA dosage profile (copies of Allele which has the specific AA residue)
      cur_aa_profile_dosage <- apply(cur_aa_profile_mat,2,function(x) {
        #If any copy is missing, cannot assign dosage at specific loci 
        # (previously #Set all aa residues in same loci as NA, every other locus as 0)
        if(any(is.na(x))){
          return(NA)
        }#Calculate number of copies present
        else{
          return(sum(x==1,na.rm = T))
        }
      })
      #Get AA dosage profile with specific genotypes (chr 1 is always odd rows, chr 2 even rows)
      cur_aa_profile_geno <- apply(cur_aa_profile_mat,2,function(x) {
        #If any copy is missing, cannot assign dosage at specific loci 
        # (previously #Set all aa residues in same loci as NA, every other locus as 0)
        if(any(is.na(x))){
          return('.|.')
        }#Calculate number of copies present
        else{
          rows_which_mapped <- which(x==1)
          #Even rows are chr 2, odd rows are chr 1
          if(length(rows_which_mapped) == 0){
            return('0|0')
          }else if((length(rows_which_mapped) == 1)){
            if((rows_which_mapped %% 2 == 0)){
              return('0|1')
            }else{
              return('1|0')
            }
          }else{
            return('1|1')
          }
        }
      })
      aa_dosage_matrix[i,] <- cur_aa_profile_dosage
      aa_geno_matrix[i,] <- cur_aa_profile_geno
    }
  }
  return(list(aa_dosage_matrix=aa_dosage_matrix,aa_geno_matrix=aa_geno_matrix))
}

WriteVCFFromAADosageMatrix <- function(aa_geno_matrix,output){
  all_aa <- colnames(aa_geno_matrix)
  all_genes <- sapply(all_aa,function(x) strsplit(x=x,split = '_')[[1]][1])
  #Construct pseudo postion, as the element index
  all_pos <- unlist(lapply(table(all_genes),function(x) seq(1,x,by=1)))
  vcf_df <- data.frame("CHROM" = all_genes,POS = all_pos,ID=all_aa,REF = 'A',ALT = 'P',QUAL = '.',FILTER = '.',INFO = 'PR',FORMAT = 'GT',stringsAsFactors = F)
  geno_df <- as.data.frame(t(aa_geno_matrix))
  vcf_out <- cbind(vcf_df,geno_df)
  colnames(vcf_out)[1] <- "#CHROM"
  
  chrom_length <- sapply(unique(all_genes),function(x) max(dplyr::filter(vcf_df,CHROM == x)$POS))
  contigs_string <- ''
  for(i in 1:length(chrom_length)){
    cur_contig <- glue::glue("##contig=<ID={names(chrom_length)[i]},length={chrom_length[i]}>")
    if(i > 1){
      contigs_string <- paste0(contigs_string,'\n',cur_contig)
    }else{
      contigs_string <- paste0(contigs_string,cur_contig)
    }
  }
  header = '##fileformat=VCFv4.3\n##FILTER=<ID=PASS,Description="All filters passed">\n##INFO=<ID=PR,Number=0,Type=Flag,Description="Provisional reference allele, may not be based on real reference genome">\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
  con <- file(output,'w')
  write(paste0(header,contigs_string),file = con,append = F)
  close(con)
  data.table::fwrite(vcf_out,file = output,append = T,col.names = T,row.names = F,quote = F,sep = '\t')
}

WriteVCFFromAlleleDosageMatrix <- function(allele_geno_matrix,output){
  all_allele <- colnames(allele_geno_matrix)
  all_genes <- sapply(all_allele,function(x) strsplit(x=x,split = '_')[[1]][1])
  #Construct pseudo postion, as the element index
  all_pos <- unlist(lapply(table(all_genes),function(x) seq(1,x,by=1)))
  
  vcf_df <- data.frame("CHROM" = all_genes,POS = all_pos,ID=all_allele,REF = 'A',ALT = 'P',QUAL = '.',FILTER = '.',INFO = 'PR',FORMAT = 'GT',stringsAsFactors = F)
  geno_df <- as.data.frame(t(allele_geno_matrix))
  vcf_out <- cbind(vcf_df,geno_df)
  colnames(vcf_out)[1] <- "#CHROM"
  
  chrom_length <- sapply(unique(all_genes),function(x) max(dplyr::filter(vcf_df,CHROM == x)$POS))
  contigs_string <- ''
  for(i in 1:length(chrom_length)){
    cur_contig <- glue::glue("##contig=<ID={names(chrom_length)[i]},length={chrom_length[i]}>")
    if(i > 1){
      contigs_string <- paste0(contigs_string,'\n',cur_contig)
    }else{
      contigs_string <- paste0(contigs_string,cur_contig)
    }
  }
  header = '##fileformat=VCFv4.3\n##FILTER=<ID=PASS,Description="All filters passed">\n##INFO=<ID=PR,Number=0,Type=Flag,Description="Provisional reference allele, may not be based on real reference genome">\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
  con <- file(output,'w')
  write(paste0(header,contigs_string),file = con,append = F)
  close(con)
  data.table::fwrite(vcf_out,file = output,append = T,col.names = T,row.names = F,quote = F,sep = '\t')
}




RunHLAGWASPyHLA <- function(POP,GT,PCA_type,LOCO,n_cores,debugging = F,sigHitsOnly=F,burden=F,out_dir,PyHLAPath,OUTCOME_sig,alt_outcome = F,hla_la_output_path){
  
  ##////////////////////////////////////////////////////////////////
  ##                        Dependencies                          //
  ##////////////////////////////////////////////////////////////////
  source(here::here("src", "prep-data.R"))
  source(here::here("src","run-GMMAT.R"))
  # source(here::here("src","run-SAIGE.R"))
  # source(here::here("src","run-SAIGE-SKAT.R"))
  library(pbmcapply)
  #Load Setup stored Variables
  DIR_SCRATCH <- RunPrepData(POP,GT,PCA_type,LOCO,n_cores,debugging,burden)
  load(glue("{DIR_SCRATCH}/prep-data.rda"))
  
  ## add string to debugging
  if (debugging) {
    NAM <- glue::glue("{NAM}_debugging")
  }
  
  
  ##///////////////////////////////////////////////////////////////
  ##                    Harmonise Identifiers                    //
  ##///////////////////////////////////////////////////////////////
  
  
  ## join with fam file,
  
  # Write out alternative outcome data ------------------------
  dat_alternative_outcome <-
    left_join(ids_order %>% dplyr::select(X1, X2),
              dat_alternative_outcome,
              by = c("X1" = "host_id")) 
  stopifnot(identical(as.character(ids_order$X1), as.character(dat_alternative_outcome$X1)))
  
  FILE_ALTERNATIVE_OUTCOME <- glue::glue("{DIR_PROCESSED}/alternative-outcome.pheno")
  
  write_delim(dat_alternative_outcome,
              path = FILE_ALTERNATIVE_OUTCOME,
              col_names = FALSE,
              delim = " ")
  
  # Write out covars data -----------------------------------
  
  ## turn into same order than genotyped data
  dat_covars_numeric <-
    left_join(ids_order %>% dplyr::select(X1, X2),
              dat_covars_numeric,
              by = c("X1" = "host_id")) 
  
  
  
  dat_covars_discrete <-
    left_join(ids_order %>% dplyr::select(X1, X2),
              dat_covars_discrete,
              by = c("X1" = "host_id")) 
  
  
  ## write out
  stopifnot(identical(as.character(ids_order$X1), as.character(dat_covars_discrete$X1)))
  stopifnot(identical(as.character(ids_order$X1), as.character(dat_covars_numeric$X1)))
  write_delim(dat_covars_discrete,
              FILE_COVARS_DISCRETE_OUT,
              delim = " ",
              col_names = FALSE)
  
  write_delim(dat_covars_numeric,
              FILE_COVARS_NUMERIC_OUT,
              delim = " ",
              col_names = FALSE)
  
  #Against HBV AA
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
  write_delim(plan,
              path = glue::glue("{DIR_SCRATCH}/plan_{NAM}.txt"),
              delim = " ")
  
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
    OUTCOME <- intersect(names(dat_pathogen), OUTCOME)
    
  }
  
  n_outcome <- length(OUTCOME)
  dat_pathogen <- dat_pathogen[,c("ID", OUTCOME)]
  
  
  
  ##///////////////////////////////////////////////////////////////
  ##                    Harmonise Identifiers                    //
  ##///////////////////////////////////////////////////////////////
  ## join with fam file
  
  # Write out pathogen data -----------------------------------
  
  dat_pathogen <-
    left_join(ids_order %>% dplyr::select(X1, X2),
              dat_pathogen,
              by = c("X1" = "ID")) 
  stopifnot(identical(as.character(ids_order$X1), as.character(dat_pathogen$X1)))
  write_delim(dat_pathogen,
              path = FILE_PATHOGEN_OUT,
              col_names = FALSE,
              delim = " ")
  
  ##///////////////////////////////////////////////////////////////
  ##                    Run PyHLA for Sig AA                    //
  ##///////////////////////////////////////////////////////////////
  OUTCOME_alt <- colnames(dat_alternative_outcome)[-c(1,2)]
  if(alt_outcome){
    OUTCOME_sig_which <- which(OUTCOME_alt %in% OUTCOME_sig)
  }else{
    OUTCOME_sig_which <- which(OUTCOME %in% OUTCOME_sig)
  }
  allele_list <- data.table::fread('~/G2G-HBV/IMGTHLA/Allelelist.txt',skip = 6)
  
  #Obtain phenotypes
  cur_pheno_tbl <- data.table::fread(FILE_PATHOGEN_OUT,stringsAsFactors = F)
  counter_within_outcome <- OUTCOME_sig_which[1]
  ind <- as.integer(counter_within_outcome + 2)
  cur_pheno_vect <- cur_pheno_tbl[,..ind]
  #Parse HLA alelles
  parsed_for_pyhla <- Parse_HLA_LA_Output(path = hla_la_output_path,
                                          samples_to_include = cur_pheno_tbl$V1,
                                          pheno = cur_pheno_vect,
                                          allele_list=allele_list,
                                          QC=F)
  #Write input file for PyHLA
  system(paste0('mkdir -p ',out_dir))
  out_path <- paste0(out_dir,'hla_allele.txt')
  data.table::fwrite(parsed_for_pyhla,file = out_path,sep = ' ',col.names = F,row.names = F,na = 'NA',quote = F)
  #Run PyHLA for first pass, obtain alleles not in reference
  cur_dir <- getwd()
  setwd(PyHLAPath)
  system(glue::glue("python2 PyHLA.py --assoc-AA --input {out_path} --out {gsub('.txt','_assocAA',out_path)} --print --consensus > {gsub('.txt','_assocAA',out_path)}.log.out"),intern = T)
  non_mapped_alleles <- system(glue::glue("grep 'no sequence is avaiable for: ' {gsub('.txt','_assocAA',out_path)}.log.out"),intern = T)
  non_mapped_alleles <- sapply(non_mapped_alleles,function(x) gsub('no sequence is avaiable for: ','',x))
  #Replace non-mapped alleles to two level alleles
  for(j in 3:ncol(parsed_for_pyhla)){
    cur_parsed_col <- parsed_for_pyhla[,j]
    mismapped_ind <- which(cur_parsed_col %in% non_mapped_alleles)
    if(length(mismapped_ind) > 0){
      new_parsed_alleles <- sapply(cur_parsed_col[which(cur_parsed_col %in% non_mapped_alleles)],function(x)ParseAllele(x,allele_list,4))
      cur_parsed_col[which(cur_parsed_col %in% non_mapped_alleles)] <- new_parsed_alleles
      parsed_for_pyhla[,j] <- cur_parsed_col
    }else{
      next
    }
  }
  #Re-write input file and Re-run PyHLA
  data.table::fwrite(parsed_for_pyhla,file = out_path,sep = ' ',col.names = F,row.names = F,na = 'NA',quote = F)
  system(glue::glue("python2 PyHLA.py --assoc-AA --input {out_path} --out {gsub('.txt','_assocAA',out_path)} --print --consensus > {gsub('.txt','_assocAA',out_path)}.log.final.out"),intern = T)
  system(glue::glue("python2 PyHLA.py --assoc --input {out_path} --out {gsub('.txt','_assocAllele',out_path)} --test logistic --model additive --print --consensus > {gsub('.txt','_assocAllele',out_path)}.log.final.out"),intern = T)
  dosage_matrix <- ConvertPyHLAToDosage(input_file_path = out_path,
                                        output_file_path  = gsub('.txt','_assocAA',out_path))
  #Load allele dosage matrix
  allele_dosage_matrix_raw <- data.table::fread(glue::glue("{gsub('.txt','_assocAllele',out_path)}.dosage"))
  #Match sample ordering 
  allele_dosage_matrix_df <- allele_dosage_matrix_raw[match(rownames(dosage_matrix$aa_geno_matrix),allele_dosage_matrix_raw$IID),]
  allele_dosage_matrix <- as.matrix(allele_dosage_matrix_df[,-c(1,2)])
  rownames(allele_dosage_matrix) <- allele_dosage_matrix_df$IID
  allele_geno_matrix <- apply(allele_dosage_matrix,2,function(x){
    x[x==1] <- '0/1'
    x[is.na(x)] <- './.'
    x[x==2] <- '1/1'
    x[x==0] <- '0/0'
    return(x)
  })
  
  WriteVCFFromAADosageMatrix(aa_geno_matrix = dosage_matrix$aa_geno_matrix,gsub('.txt','_assocAA.vcf',out_path))
  WriteVCFFromAlleleDosageMatrix(allele_geno_matrix = allele_geno_matrix,gsub('.txt','_assocAllele.vcf',out_path))
  
  setwd(cur_dir)
  # 
  # RunSAIGE_HLA_PyHLA(results_path = DIR_SCRATCH,AA_residue=OUTCOME_sig[1],LOCO=F,n_cores_input = 1,out_dir = out_dir ,vcf_path = gsub('.txt','_assocAA.vcf',out_path),create_vcf=T)
  # 
  # pbmclapply(2:length(OUTCOME_sig),function(i) {
  #   RunSAIGE_HLA_PyHLA(results_path = DIR_SCRATCH,AA_residue=OUTCOME_sig[i],LOCO=F,n_cores_input = 1,out_dir = out_dir ,vcf_path = gsub('.txt','_assocAA.vcf',out_path),create_vcf=F)
  # },mc.cores = 30)
  
  mod_type = 'logistic'
  #If alt outcome specified, replace path to pheno file 
  if(alt_outcome){
    alt_outcome_path <- glue::glue("/home/zmxu/G2G-HBV/data/processed_POP_{POP}_GT_{paste0(GT,collapse = '_')}_TYPE_pPCA_LOCO_FALSE/alternative-outcome.pheno")
    FILE_PATHOGEN_OUT <- alt_outcome_path
    mod_type = 'linear'
  }
  
  
  #Run SNPs
  id_list_path <- glue::glue("/home/zmxu/G2G-HBV/data/processed_POP_{POP}_GT_{paste0(GT,collapse = '_')}_TYPE_pPCA_LOCO_FALSE/kingunrelated.txt")
  for (i in 1:length(OUTCOME_sig)) {
    #RunSAIGE_HLA(results_path = DIR_SCRATCH,AA_residue=OUTCOME_sig[i],LOCO=F,n_cores_input = 1,OUT_DIR = out_dir ,PLINK_HLA = PLINK_HLA,run_AA=run_AA)
    counter_within_outcome <- OUTCOME_sig_which[i]
    #extract HLA Region
    system(glue::glue("{here()}/bin/plink2_linux --bfile /home/zmxu/G2G-HBV/data/raw/gilead_20181126/wes_plink/hbv_gilead --maf 0.05 --geno 0.2 --chr 6 --from-bp 28477797 --to-bp 33448354 --keep {id_list_path} --make-bed --out {out_dir}wgs_hla"))
    
    FILE_COVARS_PLINK <- glue("{DIR_PROCESSED}/covars-plink.covar")
    system(
      glue::glue(
        "{PLINK} --bfile {out_dir}wgs_hla --no-sex --{mod_type} --pheno {FILE_PATHOGEN_OUT} --pheno-col-nums {counter_within_outcome + 2}  --out {out_dir}hla_snps_{OUTCOME_sig[i]} --covar {FILE_COVARS_PLINK} --1 "
      )
    )
    if(OUTCOME_sig[i] %in% c("viralload_merged_invnormal","surface_antigen_invnormal","ALT_level_invnormal")){
      cur_result <- data.table::fread(glue::glue("{DIR_SCRATCH}/GWAS/SAIGE_SPA_TEST_LOCO_FALSE_{OUTCOME_sig[i]}"))
    }else{
      cur_result <- data.table::fread(glue::glue("{DIR_SCRATCH}/SAIGE_SPA_TEST_LOCO_FALSE_{OUTCOME_sig[i]}"))
    }
    snp_to_cond <- cur_result$SNPID[which.min(cur_result$p.value)]
    snp_list_file <- tempfile()
    write(snp_to_cond,file = snp_list_file,append = F)
    
    max_iter <- 3
    for(j in 1:max_iter){
      system(
        glue::glue(
          "{PLINK} --bfile {out_dir}wgs_hla --no-sex --{mod_type} --pheno {FILE_PATHOGEN_OUT} --pheno-col-nums {counter_within_outcome + 2}  --out {out_dir}hla_snps_{OUTCOME_sig[i]}_cond_snp_{paste0(snp_to_cond,collapse = '_')} --covar {FILE_COVARS_PLINK} --1 --condition-list {snp_list_file}"
        )
      )
      cur_result <- data.table::fread(glue::glue("{out_dir}hla_snps_{OUTCOME_sig[i]}_cond_snp_{paste0(snp_to_cond,collapse = '_')}.PHENO1.glm.{mod_type}")) %>% dplyr::filter(TEST=='ADD')
      snp_to_cond <- c(snp_to_cond,cur_result$ID[which.min(cur_result$P)])
      write(snp_to_cond,file = snp_list_file,append = F)
      
    }
    system(glue::glue("rm {snp_list_file}"))
  }
  
  #Run allele 
  for (i in 1:length(OUTCOME_sig)) {
    counter_within_outcome <- OUTCOME_sig_which[i]
    FILE_COVARS_PLINK <- glue("{DIR_PROCESSED}/covars-plink.covar")
    system(
      glue::glue(
        "{PLINK} --vcf {gsub('.txt','_assocAllele.vcf',out_path)} --no-sex --{mod_type} dominant --pheno {FILE_PATHOGEN_OUT} --pheno-col-nums {counter_within_outcome + 2}  --out {gsub('.txt',paste0('_',OUTCOME_sig[i],'_assocAllele'),out_path)} --covar {FILE_COVARS_PLINK} --1 --allow-extra-chr --maf 0.05 --geno 0.2"
      )
    )
    system(
      glue::glue(
        "{PLINK} --vcf {gsub('.txt','_assocAllele.vcf',out_path)} --no-sex --{mod_type} --pheno {FILE_PATHOGEN_OUT} --pheno-col-nums {counter_within_outcome + 2}  --out {gsub('.txt',paste0('_',OUTCOME_sig[i],'_assocAllele_add'),out_path)} --covar {FILE_COVARS_PLINK} --1 --allow-extra-chr --maf 0.05 --geno 0.2"
      )
    )
    
    cur_result <- data.table::fread(glue::glue("{gsub('.txt',paste0('_',OUTCOME_sig[i],'_assocAllele'),out_path)}.PHENO1.glm.{mod_type}")) %>% dplyr::filter(TEST=='DOM')
    snp_to_cond <- cur_result$ID[which.min(cur_result$P)]
    snp_list_file <- tempfile()
    write(snp_to_cond,file = snp_list_file,append = F)
    #For HLA alleles, run dominant model
    max_iter <- 3
    for(j in 1:max_iter){
      system(
        glue::glue(
          "{PLINK} --vcf {gsub('.txt','_assocAllele.vcf',out_path)} --no-sex --{mod_type} dominant --pheno {FILE_PATHOGEN_OUT} --pheno-col-nums {counter_within_outcome + 2}  --out {gsub('.txt',paste0('_',OUTCOME_sig[i],'_assocAllele'),out_path)}_cond_snp_{paste0(snp_to_cond,collapse = '_')} --covar {FILE_COVARS_PLINK} --1 --condition-list {snp_list_file} dominant --allow-extra-chr --maf 0.05 --geno 0.2"
        )
      )
      system(
        glue::glue(
          "{PLINK} --vcf {gsub('.txt','_assocAllele.vcf',out_path)} --no-sex --{mod_type} --pheno {FILE_PATHOGEN_OUT} --pheno-col-nums {counter_within_outcome + 2}  --out {gsub('.txt',paste0('_',OUTCOME_sig[i],'_assocAllele_add'),out_path)}_cond_snp_{paste0(snp_to_cond,collapse = '_')} --covar {FILE_COVARS_PLINK} --1 --condition-list {snp_list_file} --allow-extra-chr --maf 0.05 --geno 0.2"
        )
      )
      
      cur_result <- data.table::fread(glue::glue("{gsub('.txt',paste0('_',OUTCOME_sig[i],'_assocAllele'),out_path)}_cond_snp_{paste0(snp_to_cond,collapse = '_')}.PHENO1.glm.{mod_type}")) %>% dplyr::filter(TEST=='DOM')
      snp_to_cond <- c(snp_to_cond,cur_result$ID[which.min(cur_result$P)])
      write(snp_to_cond,file = snp_list_file,append = F)
      
    }
    system(glue::glue("rm {snp_list_file}"))
  }
  
  #Run AA logistic
  for (i in 1:length(OUTCOME_sig)) {
    counter_within_outcome <- OUTCOME_sig_which[i]
    FILE_COVARS_PLINK <- glue("{DIR_PROCESSED}/covars-plink.covar")
    system(
      glue::glue(
        "{PLINK} --vcf {gsub('.txt','_assocAA.vcf',out_path)} --no-sex --{mod_type} --pheno {FILE_PATHOGEN_OUT} --pheno-col-nums {counter_within_outcome + 2}  --out {gsub('.txt',paste0('_',OUTCOME_sig[i],'_assocAA'),out_path)} --covar {FILE_COVARS_PLINK} --1 --allow-extra-chr --maf 0.05 --geno 0.2 --ci 0.95"
      )
    )
    cur_result <- data.table::fread(glue::glue("{gsub('.txt',paste0('_',OUTCOME_sig[i],'_assocAA'),out_path)}.PHENO1.glm.{mod_type}")) %>% dplyr::filter(TEST=='ADD')
    snp_to_cond <- cur_result$ID[which.min(cur_result$P)]
    snp_list_file <- tempfile()
    write(snp_to_cond,file = snp_list_file,append = F)
    
    max_iter <- 3
    for(j in 1:max_iter){
      system(
        glue::glue(
          "{PLINK} --vcf {gsub('.txt','_assocAA.vcf',out_path)} --no-sex --{mod_type} --pheno {FILE_PATHOGEN_OUT} --pheno-col-nums {counter_within_outcome + 2}  --out {gsub('.txt',paste0('_',OUTCOME_sig[i],'_assocAA'),out_path)}_cond_snp_{paste0(snp_to_cond,collapse = '_')} --covar {FILE_COVARS_PLINK} --1 --condition-list {snp_list_file} --allow-extra-chr --maf 0.05 --geno 0.2 --ci 0.95"
        )
      )
      cur_result <- data.table::fread(glue::glue("{gsub('.txt',paste0('_',OUTCOME_sig[i],'_assocAA'),out_path)}_cond_snp_{paste0(snp_to_cond,collapse = '_')}.PHENO1.glm.{mod_type}")) %>% dplyr::filter(TEST=='ADD')
      snp_to_cond <- c(snp_to_cond,cur_result$ID[which.min(cur_result$P)])
      write(snp_to_cond,file = snp_list_file,append = F)
      
    }
    system(glue::glue("rm {snp_list_file}"))
  }
  
  #Run AA omnibus
  system(glue::glue(
    "{PLINK} --vcf {gsub('.txt','_assocAA.vcf',out_path)} --make-bed --out {gsub('.txt','_assocAA',out_path)} --allow-extra-chr --maf 0.05 --geno 0.2"
  ))
  allele_freq <- data.table::fread(glue::glue("{gsub('.txt','_assocAA',out_path)}.afreq"))
  for (i in 1:length(OUTCOME_sig)) {
    aa_id <- data.table::fread(glue::glue("{gsub('.txt','_assocAA',out_path)}.bim"),header = F)$V2
    genes <- sapply(aa_id,function(x) strsplit(x=x,split = '_')[[1]][1])
    pos <- sapply(aa_id,function(x) strsplit(x=x,split = '_')[[1]][2])
    gene_pos_id <- paste0(genes,'_',pos)
    aa_id_grouped <- vector(mode = 'list',length = length(unique(gene_pos_id)))
    for(j in 1:length(unique(gene_pos_id))){
      cur_gene_pos_id <- unique(gene_pos_id)[j]
      cur_aa_id <- aa_id[which(gene_pos_id == cur_gene_pos_id)]
      #Sort by allele freq (Major allele as reference)
      cur_aa_id_sorted <- dplyr::filter(allele_freq, ID %in% cur_aa_id) %>% dplyr::arrange(ALT_FREQS) %>% dplyr::select(ID)
      aa_id_grouped[[j]] <- cur_aa_id_sorted$ID
    }
    counter_within_outcome <- OUTCOME_sig_which[i]
    FILE_COVARS_PLINK <- glue("{DIR_PROCESSED}/covars-plink.covar")
    system(glue::glue("mkdir -p {gsub('hla_allele.txt','omnibus',out_path)}"))
    pheno_file <- data.table::fread(FILE_PATHOGEN_OUT)
    ind_to_extract <- c(1,2,counter_within_outcome + 2)
    data.table::fwrite(pheno_file[,..ind_to_extract],file = gsub('hla_allele.txt',paste0('omnibus/',OUTCOME_sig[i],'.pheno'),out_path),sep = ' ',col.names = F,row.names = F,na = "-9",quote = F)
    bim_file <- data.table::fread(glue::glue("{gsub('.txt','_assocAA',out_path)}.bim"))
    chr_new <- c()
    for(k in 1:length(unique(bim_file$V1))){
      print(k)
      chr_new <- c(chr_new,rep(k,table(bim_file$V1)[k]))
    }
    bim_file$V1 <- chr_new
    bim_file$V4 <- pos
    data.table::fwrite(bim_file,file = glue::glue("{gsub('.txt','_assocAA',out_path)}.bim"),sep = '\t',col.names = F,row.names = F)
    
    system(glue::glue("mkdir -p {gsub('hla_allele.txt',paste0('omnibus/',OUTCOME_sig[i],'/'),out_path)}"))
    #Run Omnibus test
    for(j in 1:length(aa_id_grouped)){
      cur_hap <- aa_id_grouped[[j]]
      if(mod_type == 'logistic'){
        system(
          glue::glue(
            "{gsub('plink2_linux','plink-1.07-i686/plink',PLINK)} --noweb --bfile {gsub('.txt','_assocAA',out_path)} --no-sex --pheno {gsub('hla_allele.txt',paste0('omnibus/',OUTCOME_sig[i],'.pheno'),out_path)} --out {gsub('hla_allele.txt',paste0('omnibus/',OUTCOME_sig[i],'/',paste0(cur_hap,collapse='.')),out_path)} --covar {FILE_COVARS_PLINK} --hap-snps {paste0(cur_hap,collapse=',')} --chap --1"
          )
        )
      }else{
        system(
          glue::glue(
            "{gsub('plink2_linux','plink-1.07-i686/plink',PLINK)} --noweb --bfile {gsub('.txt','_assocAA',out_path)} --no-sex --pheno {gsub('hla_allele.txt',paste0('omnibus/',OUTCOME_sig[i],'.pheno'),out_path)} --out {gsub('hla_allele.txt',paste0('omnibus/',OUTCOME_sig[i],'/',paste0(cur_hap,collapse='.')),out_path)} --covar {FILE_COVARS_PLINK} --hap-snps {paste0(cur_hap,collapse=',')} --chap"
          )
        )
      }
    }
    #Condition on reference allele (highest AF)
    omnibus_results <- ParseOmnibus(gsub('hla_allele.txt',paste0('omnibus/',OUTCOME_sig[i],'/'),out_path),parse_cond=F)
    omnibus_uncond_min_p <- dplyr::arrange(omnibus_results$uncond,P)[1,]
    cond_aa <- strsplit(rownames(omnibus_uncond_min_p),split = '\\.')[[1]][1]
    # cur_result <- data.table::fread(glue::glue("{gsub('.txt',paste0('_',OUTCOME_sig[i],'_assocAA'),out_path)}.PHENO1.glm.{mod_type}")) %>% dplyr::filter(TEST=='ADD')
    # cur_result$AA_POS <- sapply(cur_result$ID,function(x) as.numeric(strsplit(x,split = '_')[[1]][2]))
    # cur_result_omnibus_min_p <- dplyr::filter(cur_result,`#CHROM` == omnibus_uncond_min_p$`#CHROM` & AA_POS == as.numeric(strsplit(omnibus_uncond_min_p$ID,split = '_')[[1]][2]))
    # cond_aa <- cur_result_omnibus_min_p$ID[which.min(cur_result_omnibus_min_p$P)]
    # cond_aa <- cur_result$ID[which.min(cur_result$P)]
    snp_list_file <- tempfile()
    write(cond_aa,file = snp_list_file,append = F)
    
    system(
      glue::glue(
        "{PLINK} --vcf {gsub('.txt','_assocAA.vcf',out_path)} --no-sex --{mod_type} --pheno {FILE_PATHOGEN_OUT} --pheno-col-nums {counter_within_outcome + 2}  --out {gsub('.txt',paste0('_',OUTCOME_sig[i],'_assocAA'),out_path)}_cond_snp_{paste0(cond_aa,collapse = '_')} --covar {FILE_COVARS_PLINK} --1 --condition-list {snp_list_file} --allow-extra-chr --maf 0.05 --geno 0.2 --ci 0.95"
      )
    )
    
    for(j in 1:length(aa_id_grouped)){
      #Run Conditional analysis
      cur_hap <- aa_id_grouped[[j]]
      if(mod_type == 'logistic'){
        system(
          glue::glue(
            "{gsub('plink2_linux','plink-1.07-i686/plink',PLINK)} --noweb --bfile {gsub('.txt','_assocAA',out_path)} --no-sex --pheno {gsub('hla_allele.txt',paste0('omnibus/',OUTCOME_sig[i],'.pheno'),out_path)} --out {gsub('hla_allele.txt',paste0('omnibus/',OUTCOME_sig[i],'/',paste0(cur_hap,collapse='.')),out_path)}_cond_{cond_aa} --covar {FILE_COVARS_PLINK} --hap-snps {paste0(cur_hap,collapse=',')} --chap --1 --condition {cond_aa}"
          )
        )
      }else{
        system(
          glue::glue(
            "{gsub('plink2_linux','plink-1.07-i686/plink',PLINK)} --noweb --bfile {gsub('.txt','_assocAA',out_path)} --no-sex --pheno {gsub('hla_allele.txt',paste0('omnibus/',OUTCOME_sig[i],'.pheno'),out_path)} --out {gsub('hla_allele.txt',paste0('omnibus/',OUTCOME_sig[i],'/',paste0(cur_hap,collapse='.')),out_path)}_cond_{cond_aa} --covar {FILE_COVARS_PLINK} --hap-snps {paste0(cur_hap,collapse=',')} --chap --condition {cond_aa}"
          )
        )
      }
    }
  }
  
}

# RunHLAGWASPyHLA(POP = 'asian',GT = c('A','C','D','B'),PCA_type = 'pPCA',LOCO = F,n_cores = 1,debugging = F,sigHitsOnly=F,burden=F,out_dir = '/home/zmxu/G2G-HBV/PyHLA_Out/Asian/',
#                 PyHLAPath = '/home/zmxu/G2G-HBV/bin/PyHLA/',  OUTCOME_sig = "gene_Pol_pos_0049_S",hla_la_output_path = '/home/zmxu/G2G-HBV/HLA-LA_Output.csv')

# RunHLAGWASPyHLA(POP = 'asian',GT = c('A','C','D','B'),PCA_type = 'pPCA',LOCO = F,n_cores = 1,debugging = F,sigHitsOnly=F,burden=F,out_dir = '/home/zmxu/G2G-HBV/PyHLA_Out/Asian/',
#                 PyHLAPath = '/home/zmxu/G2G-HBV/bin/PyHLA/',  OUTCOME_sig = c("gene_PC_C_pos_0160_A","gene_Pol_pos_0049_N"),hla_la_output_path = '/home/zmxu/G2G-HBV/HLA-LA_Output.csv')

# RunHLAGWASPyHLA(POP = 'asian',GT = c('A','C','D','B'),PCA_type = 'pPCA',LOCO = F,n_cores = 1,debugging = F,sigHitsOnly=F,burden=F,out_dir = '/home/mxu/G2G-HBV/PyHLA_Out/Asian/',
#                 PyHLAPath = '/home/mxu/G2G-HBV/bin/PyHLA/',  OUTCOME_sig = c("viralload_merged_invnormal","surface_antigen_invnormal","ALT_level_invnormal"),alt_outcome = T,hla_la_output_path = '~/G2G-HBV/HLA-LA_Output.csv')

RunHLAGWASPyHLA(POP = 'european',GT = c('A','C','D','F','H'),PCA_type = 'pPCA',LOCO = F,n_cores = 1,debugging = F,sigHitsOnly=F,burden=F,out_dir = '/home/zmxu/G2G-HBV/PyHLA_Out/European/',
                PyHLAPath = '/home/zmxu/G2G-HBV/bin/PyHLA/',OUTCOME_sig = c("gene_PC_C_pos_0067_Y"),hla_la_output_path = '/home/zmxu/G2G-HBV/HLA-LA_Output.csv')

# RunHLAGWASPyHLA(POP = 'european',GT = c('A','C','D','F','H'),PCA_type = 'pPCA',LOCO = F,n_cores = 1,debugging = F,sigHitsOnly=F,burden=F,out_dir = '/home/mxu/G2G-HBV/PyHLA_Out/European/',
#                 PyHLAPath = '/home/mxu/G2G-HBV/bin/PyHLA/',OUTCOME_sig = c("viralload_merged_invnormal","surface_antigen_invnormal","ALT_level_invnormal"),alt_outcome = T,hla_la_output_path = '~/G2G-HBV/HLA-LA_Output.csv')


ExtractAllHLAAlleles <- function(out_dir = '/home/zmxu/G2G-HBV/PyHLA_Out/All/',PyHLAPath = '/home/zmxu/G2G-HBV/bin/PyHLA/'){
  #Parse HLA alelles
  hla_la_output <- data.table::fread('/home/zmxu/G2G-HBV/HLA-LA_Output.csv')
  allele_list <- data.table::fread('~/G2G-HBV/IMGTHLA/Allelelist.txt',skip = 6)
  parsed_for_pyhla <- Parse_HLA_LA_Output(path = '/home/zmxu/G2G-HBV/HLA-LA_Output.csv',
                                          samples_to_include = unique(hla_la_output$sample_id),
                                          pheno = rep(0,length(unique(hla_la_output$sample_id))),
                                          allele_list=allele_list,
                                          QC=T)
  #Write input file for PyHLA
  system(paste0('mkdir -p ',out_dir))
  out_path <- paste0(out_dir,'hla_allele.txt')
  data.table::fwrite(parsed_for_pyhla,file = out_path,sep = ' ',col.names = F,row.names = F,na = 'NA',quote = F)
  #Run PyHLA for first pass, obtain alleles not in reference
  cur_dir <- getwd()
  setwd(PyHLAPath)
  system(glue::glue("python2 PyHLA.py --assoc-AA --input {out_path} --out {gsub('.txt','_assocAA',out_path)} --print --consensus > {gsub('.txt','_assocAA',out_path)}.log.out"),intern = T)
  non_mapped_alleles <- system(glue::glue("grep 'no sequence is avaiable for: ' {gsub('.txt','_assocAA',out_path)}.log.out"),intern = T)
  non_mapped_alleles <- sapply(non_mapped_alleles,function(x) gsub('no sequence is avaiable for: ','',x))
  #Replace non-mapped alleles to two level alleles
  for(j in 3:ncol(parsed_for_pyhla)){
    cur_parsed_col <- parsed_for_pyhla[,j]
    mismapped_ind <- which(cur_parsed_col %in% non_mapped_alleles)
    if(length(mismapped_ind) > 0){
      new_parsed_alleles <- sapply(cur_parsed_col[which(cur_parsed_col %in% non_mapped_alleles)],function(x)ParseAllele(x,allele_list,4))
      cur_parsed_col[which(cur_parsed_col %in% non_mapped_alleles)] <- new_parsed_alleles
      parsed_for_pyhla[,j] <- cur_parsed_col
    }else{
      next
    }
  }
  #Re-write input file and Re-run PyHLA
  data.table::fwrite(parsed_for_pyhla,file = out_path,sep = ' ',col.names = F,row.names = F,na = 'NA',quote = F)
  system(glue::glue("python2 PyHLA.py --assoc-AA --input {out_path} --out {gsub('.txt','_assocAA',out_path)} --print --consensus > {gsub('.txt','_assocAA',out_path)}.log.final.out"),intern = T)
  system(glue::glue("python2 PyHLA.py --assoc --input {out_path} --out {gsub('.txt','_assocAllele',out_path)} --test logistic --model additive --print --consensus > {gsub('.txt','_assocAllele',out_path)}.log.final.out"),intern = T)
  dosage_matrix <- ConvertPyHLAToDosage(input_file_path = out_path,
                                        output_file_path  = gsub('.txt','_assocAA',out_path))
}
ExtractAllHLAAlleles()
#Load allele dosage matrix
allele_dosage_matrix_raw <- data.table::fread("~/G2G-HBV/PyHLA_Out/All/hla_allele_assocAllele.dosage")
