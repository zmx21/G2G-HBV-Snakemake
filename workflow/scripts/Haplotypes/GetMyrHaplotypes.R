#Contructs and Calculates frequency of preS1 Haplotypes, calls Shorah.

library(dplyr)
library(glue)
library(seqinr)
library(pbmcapply)
library(stringr)
library(msa)
library(pheatmap)
library(ggpubr)
library(ggplot2)
library(EnvStats)
library(gridExtra)
library(grid)

Sys.setenv(PATH=paste(Sys.getenv("PATH"), "/usr/local/shorah/bin/", sep=":"))
load('~/G2G-HBV/data/results_POP_asian_GT_A_C_D_B_TYPE_pPCA_LOCO_FALSE/prep-data.rda')
DIR_PROCESSED <- gsub(DIR_PROCESSED,pattern = '/mnt/data1/grfepc4-data3/mxu/',replacement = '~/')
LoadNTCPDosage <- function(DIR_PROCESSED){
  #Load dosage of rs2296651
  saige_vcf <- glue::glue("{DIR_PROCESSED}/hbv_gilead_QC_SAIGE.vcf.gz")
  system(glue::glue("~/Software/bcftools query -f '%CHROM %POS %ID %REF %ALT\n' {saige_vcf} > {saige_vcf}.pos"))
  pos_info <- data.table::fread(glue::glue("{saige_vcf}.pos"))
  rs2296651_pos <- dplyr::filter(pos_info,V3 == "rs2296651")
  snp_dosage_data <- data.table::fread(text = system(glue::glue('~/Software/bcftools view -r {rs2296651_pos$V1}:{rs2296651_pos$V2} {saige_vcf} | ~/Software/bcftools +dosage -- -t "GT"'),intern = T))
  colnames(snp_dosage_data) <- sapply(colnames(snp_dosage_data),function(x) gsub(x=x,pattern = "\\[.*\\]", replacement = "",))
  #Remove descript columns
  ind_to_keep <- which(!colnames(snp_dosage_data) %in% c('#CHROM','POS','REF','ALT'))
  snp_dosage_data <- snp_dosage_data[,..ind_to_keep]
  #Parse IDs
  snp_dosage_samples <- sapply(colnames(snp_dosage_data),function(x) strsplit(x,split = '_')[[1]][1])
  snp_dosage_data <- as.vector(as.numeric(snp_dosage_data))
  names(snp_dosage_data) <- snp_dosage_samples
  return(snp_dosage_data)
}
GetPreS1Coordinates <- function(GT,full_length = T){
  #Get CDS of correct reference genome
  ref_genomes <- c(B='AB219428.gff3',
                   C='GQ924620.gff3')
  #Get large S coordinates
  coords <- data.table::fread(glue::glue("~/G2G-HBV/data/raw/gilead_20181126/viral_seq/Ref_Genome/",ref_genomes[GT]),header = F) %>% dplyr::filter(V3=='CDS')
  large_S_coord <- coords[grepl(pattern = 'large surface',x = coords$V9),]
  large_S_coord <- c(large_S_coord$V4,large_S_coord$V5)
  #Read in FASTA file of reference genome
  ref_fasta <- read.fasta(glue::glue("~/G2G-HBV/data/raw/gilead_20181126/viral_seq/Ref_Genome/{gsub(pattern = '.gff3',replacement = '_NUC.fasta',x = ref_genomes[GT])}"))[[1]]
  if(large_S_coord[2] > length(ref_fasta)){
    large_S_pos <- c(large_S_coord[1]:length(ref_fasta),1:(large_S_coord[2] - length(ref_fasta)))
  }else{
    large_S_pos <- large_S_coord[1]:large_S_coord[2]
  }
  #Translate S protein, and get start codons
  S_AA <- translate(ref_fasta[large_S_pos])
  #Third START codon should encode for the end of PreS1
  PreS1_start <- large_S_pos[1]
  PreS1_end <- PreS1_start + (3*(which(S_AA=='M')[3]-1))
  if(full_length){
    Myr_start <- PreS1_start 
    Myr_end <- PreS1_start + (3*(which(S_AA=='M')[2]-1)) + (48*3)
  }else{
    Myr_start <- PreS1_start + (3*(which(S_AA=='M')[2]-1))
    Myr_end <- Myr_start + (48*3)
    
  }
  return(list(PreS1 = c(large_S_pos[1],PreS1_end),Myr = c(Myr_start,Myr_end)))
}
RunShorah <- function(sample_df,data_dir,out_path,n_cores = 1){
  #Get CDS of correct reference genome
  ref_genome_path <- gsub(data_dir,pattern = 'SAM_Files',replacement = 'Ref_Genome')
  ref_genomes <- c(B='AB219428_NUC.fasta',
                   C='GQ924620_NUC.fasta')
  cur_dir <- getwd()
  # sample_df <- dplyr::filter(sample_df,pathogen_id == "GS-US-283-1062-3134")
  #lapply(1:nrow(sample_df),function(i){
  pbmclapply(1:nrow(sample_df),function(i){
    #Set up path
    cur_id <- sample_df$pathogen_id[i]
    cur_GT <- sample_df$GT[i]
    system(glue::glue('mkdir -p {out_path}{cur_id}'))
    system(glue::glue('mkdir -p {out_path}{cur_id}/amplicon'))
    system(glue::glue('mkdir -p {out_path}{cur_id}/shotgun'))
    system(glue::glue('mkdir -p {out_path}{cur_id}/amplicon_nodup'))
    system(glue::glue('mkdir -p {out_path}{cur_id}/shotgun_nodup'))
    
    bam_path <- glue::glue("{data_dir}{cur_id}.bam")

    #Get nucleotide positions which mycludex B covers
    myr_pos <- GetPreS1Coordinates(cur_GT)$Myr
    #Get chr name from sam file
    chr_name <- system(glue::glue("~/Software/samtools view {bam_path} | head -n 1"),intern = T,ignore.stderr = T)
    chr_name <- strsplit(chr_name,split = '\t')[[1]][3]
    
    #Shift Myr pos by 1 codon to allow for indels (Manually checked that's the max length)
    #myr_pos[1] <- myr_pos[1] - 4
    #myr_pos[2] <- myr_pos[2] + 3
    
    #Run shotgun analysis 
    setwd(glue::glue('{out_path}{cur_id}/shotgun'))
    system(glue::glue("shorah shotgun -r {chr_name}:{myr_pos[1]}-{myr_pos[2]} -b {bam_path} -f {ref_genome_path}{ref_genomes[cur_GT]} -R 999 -w {myr_pos[2] - myr_pos[1] + 1} -s {(myr_pos[2] - myr_pos[1] + 1) / 2} -c 100 > log.out 2>&1"))
    out_file <- glue::glue("w-{gsub(pattern = '_NUC.fasta',replacement = '',ref_genomes[cur_GT])}.1-{myr_pos[1]}-{myr_pos[2]}.reads-support.fas")
    if(glue::glue("{out_file}.gz") %in% dir(glue::glue("{out_path}{cur_id}/shotgun/support/"))){
      system(glue::glue('gunzip {out_path}{cur_id}/shotgun/support/{out_file}.gz'))
      system(glue::glue('cp {out_path}{cur_id}/shotgun/support/{out_file} {out_path}{cur_id}/shotgun_results.fas'))
    }
    
    #Run amplicon analysis
    setwd(glue::glue('{out_path}{cur_id}/amplicon'))
    system(glue::glue("shorah amplicon -r {chr_name}:{myr_pos[1]}-{myr_pos[2]} -b {bam_path} -f {ref_genome_path}{ref_genomes[cur_GT]} -R 999 > log.out 2>&1"))
    out_file <- glue::glue("w-{gsub(pattern = '_NUC.fasta',replacement = '',ref_genomes[cur_GT])}.1-{myr_pos[1]}-{myr_pos[2]}.reads-support.fas")
    if(out_file %in% dir(glue::glue('{out_path}{cur_id}/amplicon'))){
      system(glue::glue('cp {out_path}{cur_id}/amplicon/{out_file} {out_path}{cur_id}/amplicon_results.fas'))
    }
    
    #Run shotgun analysis without PCR duplicates
    system(glue::glue('~/Software/samtools markdup -r {bam_path} {out_path}{cur_id}/{cur_id}.nodup.bam'))
    setwd(glue::glue('{out_path}{cur_id}/shotgun_nodup'))
    system(glue::glue("shorah shotgun -r {chr_name}:{myr_pos[1]}-{myr_pos[2]} -b {out_path}{cur_id}/{cur_id}.nodup.bam -f {ref_genome_path}{ref_genomes[cur_GT]} -R 999 -w {myr_pos[2] - myr_pos[1] + 1} -s {(myr_pos[2] - myr_pos[1] + 1) / 2} -c 100 > log.out 2>&1"))
    out_file <- glue::glue("w-{gsub(pattern = '_NUC.fasta',replacement = '',ref_genomes[cur_GT])}.1-{myr_pos[1]}-{myr_pos[2]}.reads-support.fas")
    if(glue::glue("{out_file}.gz") %in% dir(glue::glue("{out_path}{cur_id}/shotgun_nodup/support/"))){
      system(glue::glue('gunzip {out_path}{cur_id}/shotgun_nodup/support/{out_file}.gz'))
      system(glue::glue('cp {out_path}{cur_id}/shotgun_nodup/support/{out_file} {out_path}{cur_id}/shotgun_nodup_results.fas'))
    }    
    #Run amplicon analysis without PCR duplicates
    setwd(glue::glue('{out_path}{cur_id}/amplicon_nodup'))
    system(glue::glue("shorah amplicon -r {chr_name}:{myr_pos[1]}-{myr_pos[2]} -b {out_path}{cur_id}/{cur_id}.nodup.bam -f {ref_genome_path}{ref_genomes[cur_GT]} -R 999 > log.out 2>&1"))
    out_file <- glue::glue("w-{gsub(pattern = '_NUC.fasta',replacement = '',ref_genomes[cur_GT])}.1-{myr_pos[1]}-{myr_pos[2]}.reads-support.fas")
    if(out_file %in% dir(glue::glue('{out_path}{cur_id}/amplicon_nodup'))){
      system(glue::glue('cp {out_path}{cur_id}/amplicon_nodup/{out_file} {out_path}{cur_id}/amplicon_nodup_results.fas'))
    }
    setwd(cur_dir)
  #})
  },mc.cores = n_cores)
}
PoolHaplotypes <- function(sample_df,results_dir,posterior_thresh = 0.95,ave_reads_thresh = 10,pos_of_interest = c(17,35,51),no_dup = T){
  ProcessFASTA <- function(path,posterior_thresh,freq_thresh,ave_reads_thresh,pos_of_interest){
    #Read in FASTA file
    fasta_file <- read.fasta(path)
    #Concatenate each character vector into string
    haplotypes <- sapply(fasta_file,function(x) paste(x,collapse = ''))
    names(haplotypes) <- NULL
    #Find correct ORF (Hard coding, since we shifted added 4 bases upstream and 3 bases downstream beforehand in RunShorah)
    # start_codons <- sapply(haplotypes, function(x) stri_locate_first(x,fixed = 'atg')[1])
    # haplotypes <- sapply(1:length(haplotypes), function(i) substr(haplotypes[i],start_codons[i],nchar(haplotypes[i]) - 3 + start_codons[i] - 5))
    
    #Get haplotype annotations
    annotation <- getAnnot(fasta_file)
    ID <- gsub(x=sapply(annotation,function(x) strsplit(x,split = '\\|posterior')[[1]][1]),pattern = '>',replacement = '')
    posterior <- as.numeric(sapply(annotation,function(x) strsplit(x=strsplit(x,split = 'posterior=')[[1]][2],split = ' ave_reads')[[1]][1]))
    ave_reads <- as.numeric(sapply(annotation,function(x) strsplit(x,split = 'ave_reads=')[[1]][2]))
    prop_reads <- ave_reads / sum(ave_reads)
    #Translate to Amino acids, remove gaps.
    AA_haplotypes <- lapply(haplotypes,function(x) seqinr::translate(strsplit(x,split = '')[[1]]))
    AA_haplotypes_String <- sapply(AA_haplotypes,function(x) paste(x,collapse = ''))
    names(AA_haplotypes_String) <- NULL
    #Get amino acid residue at positions of interest
    AA_pos_of_interest <- sapply(AA_haplotypes,function(x) paste(x[pos_of_interest],collapse = ''))
    results <- data.frame(ID=ID,posterior=posterior,ave_read=ave_reads,haplotype=haplotypes,AA_haplotype=AA_haplotypes_String,AA_pos_of_interest=AA_pos_of_interest,stringsAsFactors = F)
    results_filt <- dplyr::filter(results,posterior > posterior_thresh,ave_reads > ave_reads_thresh)
    results_filt$prop_reads <- results_filt$ave_read / sum(results_filt$ave_read)
    return(results_filt)
  }
  
  #Initialize results list
  results_shotgun <- list()
  for(i in 1:nrow(sample_df)){
    #print(i)
    #Set up path
    cur_id <- sample_df$pathogen_id[i]
    cur_path <- glue::glue("{results_dir}/{cur_id}")
    if(!no_dup){
      results_shotgun[[i]] <- tryCatch(ProcessFASTA(glue::glue("{cur_path}/shotgun_results.fas"),posterior_thresh,freq_thresh,ave_reads_thresh,pos_of_interest),error = function(e){NA})
    }else{
      results_shotgun[[i]] <- tryCatch(ProcessFASTA(glue::glue("{cur_path}/shotgun_nodup_results.fas"),posterior_thresh,freq_thresh,ave_reads_thresh,pos_of_interest),error = function(e){NA})
    }
  }
  names(results_shotgun) <- sample_df$pathogen_id
  results_shotgun <- results_shotgun[sapply(results_shotgun,function(x) !all(is.na(x)))]
  results_shotgun_bind <- do.call(rbind,results_shotgun)
  
  #Reads with deletions: Manually check and correct
  NTCP_WT_ID <- c('GS-US-320-0108-1365.3','GS-US-320-0108-1370','GS-US-320-0110-4714.2','GS-US-320-0110-4819.4',
                  'GS-US-320-0110-4869.3','GS-US-320-0110-5059','GS-US-320-0110-5186.4','GS-US-283-1062-3031.1','GS-US-283-1062-3068')
  NTCP_MT_ID <- c('GS-US-320-0110-4586.8','GS-US-320-0110-4630','GS-US-320-0110-5313','GS-US-283-1062-3046.6')
  NTCP_WT_ID_B <- 'GS-US-320-0108-1221.2'
  
  if(all(sapply(NTCP_WT_ID,function(x) any(grepl(x=rownames(results_shotgun_bind),pattern=x))))){
    results_shotgun_bind['GS-US-320-0108-1365.3',]$haplotype <- gsub(x=results_shotgun_bind['GS-US-320-0108-1365.3',]$haplotype,pattern = 'gaggcct--ca',replacement = 'gaggcaactca')
    results_shotgun_bind['GS-US-320-0108-1365.3',]$AA_haplotype <- paste(seqinr::translate(strsplit(x=results_shotgun_bind['GS-US-320-0108-1365.3',]$haplotype,split = '')[[1]]),collapse = '')
    results_shotgun_bind['GS-US-320-0108-1365.3',]$AA_pos_of_interest <- paste(strsplit(x = results_shotgun_bind['GS-US-320-0108-1365.3',]$AA_haplotype,split = '')[[1]][pos_of_interest],collapse = '')
    
    results_shotgun_bind[grepl(x=rownames(results_shotgun_bind),pattern='GS-US-320-0108-1370'),]$haplotype <- sapply(results_shotgun_bind[grepl(x=rownames(results_shotgun_bind),pattern='GS-US-320-0108-1370'),]$haplotype,function(x) gsub(x=x,pattern = 'gaggcct--ca',replacement = 'gaggcaactca'),USE.NAMES = F)
    results_shotgun_bind[grepl(x=rownames(results_shotgun_bind),pattern='GS-US-320-0108-1370'),]$AA_haplotype <- sapply(results_shotgun_bind[grepl(x=rownames(results_shotgun_bind),pattern='GS-US-320-0108-1370'),]$haplotype,function(x) paste(seqinr::translate(strsplit(x=x,split = '')[[1]]),collapse = ''))
    results_shotgun_bind[grepl(x=rownames(results_shotgun_bind),pattern='GS-US-320-0108-1370'),]$AA_pos_of_interest <- sapply(results_shotgun_bind[grepl(x=rownames(results_shotgun_bind),pattern='GS-US-320-0108-1370'),]$AA_haplotype,function(x) paste(strsplit(x = x,split = '')[[1]][pos_of_interest],collapse = ''))
    
    results_shotgun_bind['GS-US-320-0110-4714.2',]$haplotype <- gsub(x=results_shotgun_bind['GS-US-320-0110-4714.2',]$haplotype,pattern = 'gaggcct--ca',replacement = 'gaggcaactca')
    results_shotgun_bind['GS-US-320-0110-4714.2',]$AA_haplotype <- paste(seqinr::translate(strsplit(x=results_shotgun_bind['GS-US-320-0110-4714.2',]$haplotype,split = '')[[1]]),collapse = '')
    results_shotgun_bind['GS-US-320-0110-4714.2',]$AA_pos_of_interest <- paste(strsplit(x = results_shotgun_bind['GS-US-320-0110-4714.2',]$AA_haplotype,split = '')[[1]][pos_of_interest],collapse = '')
    
    results_shotgun_bind['GS-US-320-0110-4819.4',]$haplotype <- gsub(x=results_shotgun_bind['GS-US-320-0110-4819.4',]$haplotype,pattern = 'gaggcct--ca',replacement = 'gaggcaactca')
    results_shotgun_bind['GS-US-320-0110-4819.4',]$AA_haplotype <- paste(seqinr::translate(strsplit(x=results_shotgun_bind['GS-US-320-0110-4819.4',]$haplotype,split = '')[[1]]),collapse = '')
    results_shotgun_bind['GS-US-320-0110-4819.4',]$AA_pos_of_interest <- paste(strsplit(x = results_shotgun_bind['GS-US-320-0110-4819.4',]$AA_haplotype,split = '')[[1]][pos_of_interest],collapse = '')
    
    results_shotgun_bind['GS-US-320-0110-4869.3',]$haplotype <- gsub(x=results_shotgun_bind['GS-US-320-0110-4869.3',]$haplotype,pattern = 'gaggcct--ca',replacement = 'gaggcaactca')
    results_shotgun_bind['GS-US-320-0110-4869.3',]$AA_haplotype <- paste(seqinr::translate(strsplit(x=results_shotgun_bind['GS-US-320-0110-4869.3',]$haplotype,split = '')[[1]]),collapse = '')
    results_shotgun_bind['GS-US-320-0110-4869.3',]$AA_pos_of_interest <- paste(strsplit(x = results_shotgun_bind['GS-US-320-0110-4869.3',]$AA_haplotype,split = '')[[1]][pos_of_interest],collapse = '')
    
    results_shotgun_bind[grepl(x=rownames(results_shotgun_bind),pattern='GS-US-320-0110-5059'),]$haplotype <- sapply(results_shotgun_bind[grepl(x=rownames(results_shotgun_bind),pattern='GS-US-320-0110-5059'),]$haplotype,function(x) gsub(x=x,pattern = 'gaggcct--ca',replacement = 'gaggcaactca'),USE.NAMES = F)
    results_shotgun_bind[grepl(x=rownames(results_shotgun_bind),pattern='GS-US-320-0110-5059'),]$AA_haplotype <- sapply(results_shotgun_bind[grepl(x=rownames(results_shotgun_bind),pattern='GS-US-320-0110-5059'),]$haplotype,function(x) paste(seqinr::translate(strsplit(x=x,split = '')[[1]]),collapse = ''))
    results_shotgun_bind[grepl(x=rownames(results_shotgun_bind),pattern='GS-US-320-0110-5059'),]$AA_pos_of_interest <- sapply(results_shotgun_bind[grepl(x=rownames(results_shotgun_bind),pattern='GS-US-320-0110-5059'),]$AA_haplotype,function(x) paste(strsplit(x = x,split = '')[[1]][pos_of_interest],collapse = ''))
    
    results_shotgun_bind['GS-US-320-0110-5186.4',]$haplotype <- gsub(x=results_shotgun_bind['GS-US-320-0110-5186.4',]$haplotype,pattern = 'gaggcct--ca',replacement = 'gaggcaactca')
    results_shotgun_bind['GS-US-320-0110-5186.4',]$AA_haplotype <- paste(seqinr::translate(strsplit(x=results_shotgun_bind['GS-US-320-0110-5186.4',]$haplotype,split = '')[[1]]),collapse = '')
    results_shotgun_bind['GS-US-320-0110-5186.4',]$AA_pos_of_interest <- paste(strsplit(x = results_shotgun_bind['GS-US-320-0110-5186.4',]$AA_haplotype,split = '')[[1]][pos_of_interest],collapse = '')
    
    results_shotgun_bind['GS-US-283-1062-3031.1',]$haplotype <- gsub(x=results_shotgun_bind['GS-US-283-1062-3031.1',]$haplotype,pattern = 'gaggcct--cagg',replacement = 'gaggcaaatctgg')
    results_shotgun_bind['GS-US-283-1062-3031.1',]$AA_haplotype <- paste(seqinr::translate(strsplit(x=results_shotgun_bind['GS-US-283-1062-3031.1',]$haplotype,split = '')[[1]]),collapse = '')
    results_shotgun_bind['GS-US-283-1062-3031.1',]$AA_pos_of_interest <- paste(strsplit(x = results_shotgun_bind['GS-US-283-1062-3031.1',]$AA_haplotype,split = '')[[1]][pos_of_interest],collapse = '')
    
    
    results_shotgun_bind[grepl(x=rownames(results_shotgun_bind),pattern='GS-US-283-1062-3068'),]$haplotype <- sapply(results_shotgun_bind[grepl(x=rownames(results_shotgun_bind),pattern='GS-US-283-1062-3068'),]$haplotype,function(x) gsub(x=x,pattern = 'gaggcct--ca',replacement = 'gaggcaactca'),USE.NAMES = F)
    results_shotgun_bind[grepl(x=rownames(results_shotgun_bind),pattern='GS-US-283-1062-3068'),]$AA_haplotype <- sapply(results_shotgun_bind[grepl(x=rownames(results_shotgun_bind),pattern='GS-US-283-1062-3068'),]$haplotype,function(x) paste(seqinr::translate(strsplit(x=x,split = '')[[1]]),collapse = ''))
    results_shotgun_bind[grepl(x=rownames(results_shotgun_bind),pattern='GS-US-283-1062-3068'),]$AA_pos_of_interest <- sapply(results_shotgun_bind[grepl(x=rownames(results_shotgun_bind),pattern='GS-US-283-1062-3068'),]$AA_haplotype,function(x) paste(strsplit(x = x,split = '')[[1]][pos_of_interest],collapse = ''))
    
  }else if(all(sapply(NTCP_MT_ID,function(x) any(grepl(x=rownames(results_shotgun_bind),pattern=x))))){
    results_shotgun_bind['GS-US-320-0110-4586.8',]$haplotype <- gsub(x=results_shotgun_bind['GS-US-320-0110-4586.8',]$haplotype,pattern = 'gaggact--cagg',replacement = 'gaggtaactcagg')
    results_shotgun_bind['GS-US-320-0110-4586.8',]$AA_haplotype <- paste(seqinr::translate(strsplit(x=results_shotgun_bind['GS-US-320-0110-4586.8',]$haplotype,split = '')[[1]]),collapse = '')
    results_shotgun_bind['GS-US-320-0110-4586.8',]$AA_pos_of_interest <- paste(strsplit(x = results_shotgun_bind['GS-US-320-0110-4586.8',]$AA_haplotype,split = '')[[1]][pos_of_interest],collapse = '')
    
    results_shotgun_bind[grepl(x=rownames(results_shotgun_bind),pattern='GS-US-320-0110-4630'),]$haplotype <- sapply(results_shotgun_bind[grepl(x=rownames(results_shotgun_bind),pattern='GS-US-320-0110-4630'),]$haplotype,function(x) gsub(x=x,pattern = 'accc-gctt',replacement = 'acccgctat'),USE.NAMES = F)
    results_shotgun_bind[grepl(x=rownames(results_shotgun_bind),pattern='GS-US-320-0110-4630'),]$haplotype <- sapply(results_shotgun_bind[grepl(x=rownames(results_shotgun_bind),pattern='GS-US-320-0110-4630'),]$haplotype,function(x) gsub(x=x,pattern = 'ggc----ca',replacement = 'ggaaaacca'),USE.NAMES = F)
    results_shotgun_bind[grepl(x=rownames(results_shotgun_bind),pattern='GS-US-320-0110-4630'),]$AA_haplotype <- sapply(results_shotgun_bind[grepl(x=rownames(results_shotgun_bind),pattern='GS-US-320-0110-4630'),]$haplotype,function(x) paste(seqinr::translate(strsplit(x=x,split = '')[[1]]),collapse = ''))
    results_shotgun_bind[grepl(x=rownames(results_shotgun_bind),pattern='GS-US-320-0110-4630'),]$AA_pos_of_interest <- sapply(results_shotgun_bind[grepl(x=rownames(results_shotgun_bind),pattern='GS-US-320-0110-4630'),]$AA_haplotype,function(x) paste(strsplit(x = x,split = '')[[1]][pos_of_interest],collapse = ''))

    results_shotgun_bind[grepl(x=rownames(results_shotgun_bind),pattern='GS-US-320-0110-5313'),]$haplotype <- sapply(results_shotgun_bind[grepl(x=rownames(results_shotgun_bind),pattern='GS-US-320-0110-5313'),]$haplotype,function(x) gsub(x=x,pattern = 'accctg--ttc',replacement = 'acccattgttc'),USE.NAMES = F)
    results_shotgun_bind[grepl(x=rownames(results_shotgun_bind),pattern='GS-US-320-0110-5313'),]$AA_haplotype <- sapply(results_shotgun_bind[grepl(x=rownames(results_shotgun_bind),pattern='GS-US-320-0110-5313'),]$haplotype,function(x) paste(seqinr::translate(strsplit(x=x,split = '')[[1]]),collapse = ''))
    results_shotgun_bind[grepl(x=rownames(results_shotgun_bind),pattern='GS-US-320-0110-5313'),]$AA_pos_of_interest <- sapply(results_shotgun_bind[grepl(x=rownames(results_shotgun_bind),pattern='GS-US-320-0110-5313'),]$AA_haplotype,function(x) paste(strsplit(x = x,split = '')[[1]][pos_of_interest],collapse = ''))
    
    results_shotgun_bind['GS-US-283-1062-3046.6',]$haplotype <- gsub(x=results_shotgun_bind['GS-US-283-1062-3046.6',]$haplotype,pattern = 'ggcat--ca',replacement = 'gggacatca')
    results_shotgun_bind['GS-US-283-1062-3046.6',]$AA_haplotype <- paste(seqinr::translate(strsplit(x=results_shotgun_bind['GS-US-283-1062-3046.6',]$haplotype,split = '')[[1]]),collapse = '')
    results_shotgun_bind['GS-US-283-1062-3046.6',]$AA_pos_of_interest <- paste(strsplit(x = results_shotgun_bind['GS-US-283-1062-3046.6',]$AA_haplotype,split = '')[[1]][pos_of_interest],collapse = '')
    
  }else if(all(sapply(NTCP_WT_ID_B,function(x) any(grepl(x=rownames(results_shotgun_bind),pattern=x))))){
    results_shotgun_bind['GS-US-320-0108-1221.2',]$haplotype <- gsub(x=results_shotgun_bind['GS-US-320-0108-1221.2',]$haplotype,pattern = 'ctc--acaatc',replacement = 'ctcacaaactc')
    results_shotgun_bind['GS-US-320-0108-1221.2',]$AA_haplotype <- paste(seqinr::translate(strsplit(x=results_shotgun_bind['GS-US-320-0108-1221.2',]$haplotype,split = '')[[1]]),collapse = '')
    results_shotgun_bind['GS-US-320-0108-1221.2',]$AA_pos_of_interest <- paste(strsplit(x = results_shotgun_bind['GS-US-320-0108-1221.2',]$AA_haplotype,split = '')[[1]][pos_of_interest],collapse = '')
  }
  
  
  unique_pos_of_interest <- unique(results_shotgun_bind$AA_pos_of_interest)
  prop_pos_of_interest <- rep(NA,length(unique_pos_of_interest))
  names(prop_pos_of_interest) <- unique_pos_of_interest
  for(j in 1:length(unique_pos_of_interest)){
    cur_haplo <- dplyr::filter(results_shotgun_bind,AA_pos_of_interest == unique_pos_of_interest[j])
    prop_pos_of_interest[j] <- sum(cur_haplo$prop_read) / length(results_shotgun)
  }
  results_shotgun_bind$ID <- sapply(rownames(results_shotgun_bind),function(x) strsplit(x=x,split = '\\.')[[1]][1])
  return(list(prop_pos_of_interest = prop_pos_of_interest,results_shotgun_bind=results_shotgun_bind))
}
GetReadWidth <- function(sample_df,data_dir){
  width <- rep(NA,nrow(sample_df))
  for(i in 1:nrow(sample_df)){
    cur_id <- sample_df$pathogen_id[i]
    bam_path <- glue::glue("{data_dir}{cur_id}.bam")
    if(!glue::glue("{bam_path}.readwidth") %in% dir(data_dir)){
      system(glue::glue("/home/zmxu/Software/samtools stats {bam_path} | grep ^RL | cut -f 2-3 > {bam_path}.readwidth"))
    }
    cur_width <- data.table::fread(glue::glue("{bam_path}.readwidth"))
    width[i] <- weighted.mean(cur_width$V1,cur_width$V2)
  }
  return(width)
}


#Get NTCP Dosages, Limit to genotype B and C
NTCP_Dosage <- LoadNTCPDosage(DIR_PROCESSED)
#Get NTCP Mutants and Wildtypes
NTCP_WT <- dplyr::filter(dat_covars_raw,host_id %in% names(NTCP_Dosage[which(NTCP_Dosage == 0)])) %>% dplyr::filter(GT %in% c('B','C'))
NTCP_Mutant <- dplyr::filter(dat_covars_raw,host_id %in% names(NTCP_Dosage[which(NTCP_Dosage == 1)])) %>% dplyr::filter(GT %in% c('B','C'))

#Get Read Width
NTCP_WT_width <- GetReadWidth(NTCP_WT,'~/G2G-HBV/data/raw/gilead_20181126/viral_seq/SAM_Files/')
NTCP_Mutant_width <- GetReadWidth(NTCP_Mutant,'~/G2G-HBV/data/raw/gilead_20181126/viral_seq/SAM_Files/')


#Run Shorah for NTCP Mutants and WT
RunShorah(NTCP_Mutant,'~/G2G-HBV/data/raw/gilead_20181126/viral_seq/SAM_Files/','~/G2G-HBV/MyrcludexB_Full_Length_Haplotypes/',n_cores=20)
RunShorah(NTCP_WT,'~/G2G-HBV/data/raw/gilead_20181126/viral_seq/SAM_Files/','~/G2G-HBV/MyrcludexB_Full_Length_Haplotypes/',n_cores=20)

#Pool intra-host populations into Multiple Sequence Alignment, separately for NTCP WT and Mutants
prop_geno_C <- PoolHaplotypes(NTCP_Mutant %>% dplyr::filter(GT == 'C'),'~/G2G-HBV/MyrcludexB_Full_Length_Haplotypes/',no_dup = F)
prop_geno_B <- PoolHaplotypes(NTCP_Mutant %>% dplyr::filter(GT == 'B'),'~/G2G-HBV/MyrcludexB_Full_Length_Haplotypes/',no_dup = F)

prop_WT_geno_C <- PoolHaplotypes(NTCP_WT %>% dplyr::filter(GT == 'C'),'~/G2G-HBV/MyrcludexB_Full_Length_Haplotypes/',no_dup = F)
prop_WT_geno_B <- PoolHaplotypes(NTCP_WT %>% dplyr::filter(GT == 'B'),'~/G2G-HBV/MyrcludexB_Full_Length_Haplotypes/',no_dup = F)

ContructTable <- function(prop_vect,thresh = 0){
  prop_vect_df <- data.frame(Haplotype = names(prop_vect),Prop = signif(prop_vect,2),stringsAsFactors = F) %>% dplyr::filter(Prop > thresh)
  pos_17 <- sapply(prop_vect_df$Haplotype,function(x) strsplit(x=x,split = '')[[1]][1])
  if(any(pos_17 == 'A')){
    pos_17_A = prop_vect_df[pos_17 == 'A',] %>% dplyr::arrange(desc(Prop))
    pos_17_A$Rel_Prop <- signif(pos_17_A$Prop / sum(pos_17_A$Prop),2)
    rownames(pos_17_A) <- NULL
    pos_17_S = prop_vect_df[pos_17 == 'S',]  %>% dplyr::arrange(desc(Prop))
    pos_17_S$Rel_Prop <- signif(pos_17_S$Prop / sum(pos_17_S$Prop),2)
    rownames(pos_17_S) <- NULL
    return(list(pos_17_A=pos_17_A,pos_17_S = pos_17_S))
  }else{
    pos_17_S = prop_vect_df %>% dplyr::arrange(desc(Prop))
    rownames(pos_17_S) <- NULL
    return(list(pos_17_A=NA,pos_17_S = pos_17_S))
  }
}

prop_tbl_geno_C <- ContructTable(prop_geno_C$prop_pos_of_interest)
prop_tbl_geno_C_WT <- ContructTable(prop_WT_geno_C$prop_pos_of_interest)

prop_tbl_geno_B <- ContructTable(prop_geno_B$prop_pos_of_interest)
prop_tbl_geno_B_WT <- ContructTable(prop_WT_geno_B$prop_pos_of_interest)



arp_geno_C <- prop_geno_C$results_shotgun_bind %>% dplyr::filter(AA_pos_of_interest == 'ARP')
sgh_geno_C <- prop_WT_geno_C$results_shotgun_bind %>% dplyr::filter(AA_pos_of_interest == 'SGH')
sgh_geno_C_MT <- prop_geno_C$results_shotgun_bind %>% dplyr::filter(AA_pos_of_interest == 'SGH')

akn_geno_B <- prop_geno_B$results_shotgun_bind %>% dplyr::filter(AA_pos_of_interest == 'AKN')
skn_geno_B <- prop_WT_geno_B$results_shotgun_bind %>% dplyr::filter(AA_pos_of_interest == 'SKN')

GetPropBySample <- function(results){
  by_ID <- data.frame(ID=character(),prop_reads = numeric(),AA_pos_of_interest = character(),stringsAsFactors = F)
  for(i in 1:length(unique(results$ID))){
    cur_df <- dplyr::filter(results,ID == unique(results$ID)[i])
    cur_haplo <- unique(cur_df$AA_pos_of_interest)
    for(j in 1:length(cur_haplo)){
      cur_haplo_df <- dplyr::filter(cur_df,AA_pos_of_interest == cur_haplo[j])
      by_ID <- rbind(by_ID,data.frame(ID = unique(results$ID)[i],prop_reads = sum(cur_haplo_df$prop_reads),AA_pos_of_interest = cur_haplo[j],stringsAsFactors = F))
    }
  }
  return(by_ID)
}
indels <- unique(data.table::fread('~/G2G-HBV/with_indel.txt',header = F)$V1)
indels <- sapply(indels,function(x) strsplit(x = x,split = '/')[[1]][3])
indels <- sapply(indels,function(x) gsub(x=x,pattern = '.fasta',replacement = ''))

prop_geno_C_by_sample <- GetPropBySample(prop_geno_C$results_shotgun_bind)
prop_geno_C_WT_by_sample <- GetPropBySample(prop_WT_geno_C$results_shotgun_bind)

prop_geno_B_by_sample <- GetPropBySample(prop_geno_B$results_shotgun_bind)
prop_geno_B_WT_by_sample <- GetPropBySample(prop_WT_geno_B$results_shotgun_bind)

prop_geno_C_by_sample_ARP <- prop_geno_C_by_sample %>% dplyr::filter(AA_pos_of_interest == 'ARP' & prop_reads > 0.75)
prop_geno_C_by_sample_SGH <- prop_geno_C_by_sample %>% dplyr::filter(AA_pos_of_interest == 'SGH' & prop_reads > 0.75)
prop_geno_C_by_sample_SRP <- prop_geno_C_by_sample %>% dplyr::filter(AA_pos_of_interest == 'SRP' & prop_reads > 0.75)
prop_geno_C_by_sample_SGQ <- prop_geno_C_by_sample %>% dplyr::filter(AA_pos_of_interest == 'SGQ' & prop_reads > 0.75)

prop_geno_C_WT_by_sample_SGH <- prop_geno_C_WT_by_sample %>% dplyr::filter(AA_pos_of_interest == 'SGH' & prop_reads > 0.75 & !ID %in% indels)
prop_geno_C_WT_by_sample_SGQ <- prop_geno_C_WT_by_sample %>% dplyr::filter(AA_pos_of_interest == 'SGQ' & prop_reads > 0.75 & !ID %in% indels)

prop_geno_B_by_sample_AKN <- prop_geno_B_by_sample %>% dplyr::filter(AA_pos_of_interest == 'AKN' & prop_reads > 0.75)
prop_geno_B_by_sample_ARN <- prop_geno_B_by_sample %>% dplyr::filter(AA_pos_of_interest == 'ARN' & prop_reads > 0.75)
prop_geno_B_WT_by_sample_SKN <- prop_geno_B_WT_by_sample %>% dplyr::filter(AA_pos_of_interest == 'SKN' & prop_reads > 0.75)


Myr_ARP_aln <- msa(AAStringSet(rep(arp_geno_C$AA_haplotype,round(10 * arp_geno_C$prop_reads))),'ClustalW')
msaConsensusSequence(Myr_ARP_aln,type = 'upperlower')

write.fasta(as.list(rep(arp_geno_C$AA_haplotype,round(10 * arp_geno_C$prop_reads))),
            names = seq(1,sum(round(10 * arp_geno_C$prop_reads))),file.out = '~/G2G-HBV/MyrcludexB_Full_Length_Haplotypes/ARP.fasta')


Myr_SGH_aln <- msa(AAStringSet(rep(sgh_geno_C_MT$AA_haplotype,round(10 * sgh_geno_C_MT$prop_reads))),'ClustalW')
toupper(msaConsensusSequence(Myr_SGH_aln,type = 'upperlower'))

Myr_SGH_WT_aln <- msa(AAStringSet(rep(sgh_geno_C$AA_haplotype,round(10 * sgh_geno_C$prop_reads))),'ClustalW')
msaConsensusSequence(Myr_SGH_WT_aln,type = 'upperlower')

write.fasta(as.list(rep(sgh_geno_C$AA_haplotype,round(10 * sgh_geno_C$prop_reads))),
            names = seq(1,sum(round(10 * sgh_geno_C$prop_reads))),file.out = '~/G2G-HBV/MyrcludexB_Full_Length_Haplotypes/SGH.fasta')



srp_geno_C <- prop_geno_C$results_shotgun_bind %>% dplyr::filter(AA_pos_of_interest == 'SRP')
Myr_SRP_aln <- msa(AAStringSet(rep(srp_geno_C$AA_haplotype,round(10 * srp_geno_C$prop_reads))),'ClustalW')
msaConsensusSequence(Myr_SRP_aln,type = 'upperlower')

sgq_geno_C <- prop_geno_C$results_shotgun_bind %>% dplyr::filter(AA_pos_of_interest == 'SGQ')
Myr_sgq_aln <- msa(AAStringSet(rep(sgq_geno_C$AA_haplotype,round(10 * sgq_geno_C$prop_reads))),'ClustalW')
msaConsensusSequence(Myr_sgq_aln,type = 'upperlower')

sgq_WT_geno_C <- prop_WT_geno_C$results_shotgun_bind %>% dplyr::filter(AA_pos_of_interest == 'SGQ')
Myr_sgq_WT_aln <- msa(AAStringSet(rep(sgq_WT_geno_C$AA_haplotype,round(10 * sgq_WT_geno_C$prop_reads))),'ClustalW')
msaConsensusSequence(Myr_sgq_WT_aln,type = 'upperlower')

akn_geno_B <- prop_geno_B$results_shotgun_bind %>% dplyr::filter(AA_pos_of_interest == 'AKN')
Myr_akn_aln <- msa(AAStringSet(rep(akn_geno_B$AA_haplotype,round(10 * akn_geno_B$prop_reads))),'ClustalW')
toupper(msaConsensusSequence(Myr_akn_aln,type = 'upperlower'))


write.fasta(as.list(rep(akn_geno_B$AA_haplotype,round(10 * akn_geno_B$prop_reads))),
            names = seq(1,sum(round(10 * akn_geno_B$prop_reads))),file.out = '~/G2G-HBV/MyrcludexB_Full_Length_Haplotypes/AKN.fasta')



skn_WT_geno_B <- prop_WT_geno_B$results_shotgun_bind %>% dplyr::filter(AA_pos_of_interest == 'SKN')
Myr_skn_WT_aln <- msa(AAStringSet(rep(skn_WT_geno_B$AA_haplotype,round(10 * skn_WT_geno_B$prop_reads))),'ClustalW')
msaConsensusSequence(Myr_skn_WT_aln,type = 'upperlower')

write.fasta(as.list(rep(skn_WT_geno_B$AA_haplotype,round(10 * skn_WT_geno_B$prop_reads))),
            names = seq(1,sum(round(10 * skn_WT_geno_B$prop_reads))),file.out = '~/G2G-HBV/MyrcludexB_Full_Length_Haplotypes/SKN.fasta')


skn_geno_B <- prop_geno_B$results_shotgun_bind %>% dplyr::filter(AA_pos_of_interest == 'SKN')
Myr_skn_aln <- msa(AAStringSet(rep(skn_geno_B$AA_haplotype,round(10 * skn_geno_B$prop_reads))),'ClustalW')
msaConsensusSequence(Myr_skn_aln,type = 'upperlower')

arn_geno_B <- prop_geno_B$results_shotgun_bind %>% dplyr::filter(AA_pos_of_interest == 'ARN')
Myr_arn_aln <- msa(AAStringSet(rep(arn_geno_B$AA_haplotype,round(10 * arn_geno_B$prop_reads))),'ClustalW')
toupper(msaConsensusSequence(Myr_arn_aln,type = 'upperlower'))


msaPrettyPrint(msa(AAStringSet(c(toupper(msaConsensusSequence(Myr_akn_aln,type = 'upperlower'))
                                 ,akn_nuc_consensus$AA_Seq[1]))))

prop_geno_C_SGH_SGQ <- lapply(unique(prop_geno_C$results_shotgun_bind$ID),function(x) {
  # df <-  dplyr::filter(prop_geno_C$results_shotgun_bind,ID == x) %>% dplyr::group_by(AA_pos_of_interest) %>% dplyr::summarize(Sum = sum(prop_reads)) %>% dplyr::arrange(desc(Sum))
  # return(df$AA_pos_of_interest[1])
  WT_Haplo <-  prop_WT_geno_C$results_shotgun_bind 
  df <- dplyr::filter(prop_geno_C$results_shotgun_bind,ID == x) %>% dplyr::filter(!AA_pos_of_interest %in% unique(WT_Haplo)$AA_pos_of_interest)
  if(nrow(df) > 0){
    return(sum(df$prop_reads))
  }else{
    return(0)
  }
  })
prop_geno_C_WT_SGH_SGQ <- lapply(unique(prop_WT_geno_C$results_shotgun_bind$ID),function(x) {
  WT_Haplo <-  prop_WT_geno_C$results_shotgun_bind 
  df <- dplyr::filter(prop_WT_geno_C$results_shotgun_bind,ID == x) %>% dplyr::filter(!AA_pos_of_interest %in% unique(WT_Haplo)$AA_pos_of_interest)
  if(nrow(df) > 0){
    return(sum(df$prop_reads))
  }else{
    return(0)
  }
})

prop_geno_C_SGH_SGQ <- rbind(data.frame(ID = unique(prop_geno_C$results_shotgun_bind$ID),Freq = unlist(prop_geno_C_SGH_SGQ),rs2296651 = 'GA'),
                             data.frame(ID = unique(prop_WT_geno_C$results_shotgun_bind$ID),Freq = unlist(prop_geno_C_WT_SGH_SGQ),rs2296651 = 'GG'))

prop_geno_B_SKN <- lapply(unique(prop_geno_B$results_shotgun_bind$ID),function(x) {
  WT_Haplo <-  prop_WT_geno_B$results_shotgun_bind 
  
  df <- dplyr::filter(prop_geno_B$results_shotgun_bind,ID == x) %>% dplyr::filter(!AA_pos_of_interest %in% unique(WT_Haplo)$AA_pos_of_interest)
  if(nrow(df) > 0){
    return(sum(df$prop_reads))
  }else{
    return(0)
  }
})
prop_geno_B_WT_SKN <- lapply(unique(prop_WT_geno_B$results_shotgun_bind$ID),function(x) {
  WT_Haplo <-  prop_WT_geno_B$results_shotgun_bind 
  
  df <- dplyr::filter(prop_WT_geno_B$results_shotgun_bind,ID == x) %>% dplyr::filter(!AA_pos_of_interest %in% unique(WT_Haplo)$AA_pos_of_interest)
  if(nrow(df) > 0){
    return(sum(df$prop_reads))
  }else{
    return(0)
  }
})

prop_geno_B_SKN <- rbind(data.frame(ID = unique(prop_geno_B$results_shotgun_bind$ID),Freq = unlist(prop_geno_B_SKN),rs2296651 = 'GA'),
                         data.frame(ID = unique(prop_WT_geno_B$results_shotgun_bind$ID),Freq = unlist(prop_geno_B_WT_SKN),rs2296651 = 'GG'))