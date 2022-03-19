fig_label <- function(text, region="figure", pos="topleft", cex=NULL, ...) {
  
  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright", 
                          "left", "center", "right", 
                          "bottomleft", "bottom", "bottomright"))
  
  if(region %in% c("figure", "device")) {
    ds <- dev.size("in")
    # xy coordinates of device corners in user coordinates
    x <- grconvertX(c(0, ds[1]), from="in", to="user")
    y <- grconvertY(c(0, ds[2]), from="in", to="user")
    
    # fragment of the device we use to plot
    if(region == "figure") {
      # account for the fragment of the device that 
      # the figure is using
      fig <- par("fig")
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- x[1] + dx * fig[1:2]
      y <- y[1] + dy * fig[3:4]
    } 
  }
  
  # much simpler if in plotting region
  if(region == "plot") {
    u <- par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }
  
  sw <- strwidth(text, cex=cex) * 60/100
  sh <- strheight(text, cex=cex) * 60/100
  
  x1 <- switch(pos,
               topleft     =x[1] + sw, 
               left        =x[1] + sw,
               bottomleft  =x[1] + sw,
               top         =(x[1] + x[2])/2,
               center      =(x[1] + x[2])/2,
               bottom      =(x[1] + x[2])/2,
               topright    =x[2] - sw,
               right       =x[2] - sw,
               bottomright =x[2] - sw)
  
  y1 <- switch(pos,
               topleft     =y[2] - sh,
               top         =y[2] - sh,
               topright    =y[2] - sh,
               left        =(y[1] + y[2])/2,
               center      =(y[1] + y[2])/2,
               right       =(y[1] + y[2])/2,
               bottomleft  =y[1] + sh,
               bottom      =y[1] + sh,
               bottomright =y[1] + sh)
  
  old.par <- par(xpd=NA)
  on.exit(par(old.par))
  
  text(x1, y1, text, cex=cex, ...)
  return(invisible(c(x,y)))
}
library(dplyr)
library(scales)
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

PlotHLAAssoc <- function(cond_files_snps,cond_files_AA,cond_files_allele,omnibus= F,cond_aa_pos = NA,bonf_thresh){
  snp_results <- list(uncond = data.table::fread(cond_files_snps$uncond) %>% dplyr::select(`#CHROM`=CHR,POS,ID=SNPID,REF=Allele1,ALT=Allele2,A1=Allele2,N.Cases,N.Controls,BETA,SE=SE,T_STAT=Tstat,P=p.value) %>% dplyr::mutate(OR=exp(BETA),OBS_CT=N.Cases+N.Controls,TEST='ADD') %>% dplyr::select(-N.Cases,-N.Controls,-BETA) %>%
                        dplyr::filter(`#CHROM`==6 &  POS >= 28477797 & POS <= 33448354),
                      cond=data.table::fread(cond_files_snps$cond) %>% dplyr::filter(TEST == 'ADD'))
  allele_results <- list(uncond = data.table::fread(cond_files_allele$uncond) %>% dplyr::select(`#CHROM`=CHR,POS,ID=SNPID,REF=Allele1,ALT=Allele2,A1=Allele2,N.Cases,N.Controls,BETA,SE=SE,T_STAT=Tstat,P=p.value) %>% dplyr::mutate(OR=exp(BETA),OBS_CT=N.Cases+N.Controls,TEST='DOM') %>% dplyr::select(-N.Cases,-N.Controls,-BETA),
                         cond=data.table::fread(cond_files_allele$cond) %>% dplyr::filter(TEST == 'DOM'))
  
  if(!omnibus){
    aa_results <- lapply(cond_files_AA,function(x) data.table::fread(x) %>% dplyr::filter(TEST == 'ADD'))
  }else{
    aa_results <- ParseOmnibus(cond_files_AA)
    #Remove positions with no variation
    aa_results <- lapply(aa_results,function(x) x[sapply(rownames(aa_results$uncond),function(x) length(strsplit(x=x,split = '\\.')[[1]])) != 2,])
  }
  
  #Dowload GTF from:ftp://ftp.ensembl.org/pub/grch37/release-87/gtf/homo_sapiens/
  gtf <- rtracklayer::import('~/G2G-HBV/Homo_sapiens.GRCh37.87.gtf')
  gtf_df=as.data.frame(gtf)
  hla_gtf <- dplyr::filter(gtf_df,grepl('HLA',gene_name)) %>% dplyr::filter(seqnames == '6')
  
  genes <- data.frame(Gene = unique(hla_gtf$gene_name),Start = NA,End = NA,Strand=NA)
  for(i in 1:nrow(genes)){
    cur_gene <- dplyr::filter(hla_gtf,gene_name == genes$Gene[i] & type == 'gene')
    genes$Start[i] <- min(cur_gene$start)
    genes$End[i] <- max(cur_gene$end)
    genes$Strand[i] <- unique(as.character(cur_gene$strand))
    
  }
  genes <- dplyr::filter(genes, Gene %in% paste0('HLA-',unique(allele_results[[1]]$`#CHROM`)))
  minP <- max(-log10(min(do.call(rbind,aa_results)$P,na.rm = T)) + 1,max(c(-log10(min(do.call(rbind,snp_results)$P,na.rm = T)) + 1,-log10(min(do.call(rbind,allele_results)$P,na.rm = T)) + 1)))
  minY <- 2
  uncond_snp <- snp_results$uncond
  cond_snp <- snp_results$cond
  #Plot AAs
  uncond_aa <- aa_results$uncond %>% dplyr::mutate(Gene = paste0('HLA-',`#CHROM`)) %>% dplyr::left_join(genes,by=c('Gene'='Gene'))
  uncond_aa$AA_POS <- as.numeric(sapply(uncond_aa$ID,function(x) strsplit(x=x,split = '_')[[1]][2]))
  cond_aa <- aa_results$cond %>% dplyr::mutate(Gene = paste0('HLA-',`#CHROM`)) %>% dplyr::left_join(genes,by=c('Gene'='Gene')) 
  cond_aa$AA_POS <- as.numeric(sapply(cond_aa$ID,function(x) strsplit(x=x,split = '_')[[1]][2]))
  cond_aa$Nuc_Pos <- rep(NA,nrow(cond_aa))
  uncond_aa$Nuc_Pos <- rep(NA,nrow(uncond_aa))
  for(i in 1:nrow(cond_aa)){
    cur_hla_locus <- dplyr::filter(cond_aa,Gene == cond_aa$Gene[i])
    cur_aa_pos <- cond_aa$AA_POS[i]
    cur_start <- cond_aa$Start[i]
    cur_end <- cond_aa$End[i]
    
    cur_ratio_factor <- cur_start + ((cur_end -cur_start)/max(cur_hla_locus$AA_POS)*cur_aa_pos)
    cond_aa$Nuc_Pos[i] <- cur_ratio_factor
  }
  for(i in 1:nrow(uncond_aa)){
    cur_hla_locus <- dplyr::filter(uncond_aa,Gene == uncond_aa$Gene[i])
    cur_aa_pos <- uncond_aa$AA_POS[i]
    cur_start <- uncond_aa$Start[i]
    cur_end <- uncond_aa$End[i]
    
    cur_ratio_factor <- cur_start + ((cur_end -cur_start)/max(cur_hla_locus$AA_POS)*cur_aa_pos)
    uncond_aa$Nuc_Pos[i] <- cur_ratio_factor
  }
  
  #Plot alleles
  uncond_alleles <- allele_results$uncond %>% dplyr::mutate(Gene = paste0('HLA-',`#CHROM`)) %>% dplyr::left_join(genes,by=c('Gene'='Gene')) %>% dplyr::mutate(POS = round((Start + End) / 2))
  cond_alleles <- allele_results$cond %>% dplyr::mutate(Gene = paste0('HLA-',`#CHROM`)) %>% dplyr::left_join(genes,by=c('Gene'='Gene')) %>% dplyr::mutate(POS = round((Start + End) / 2))
  
  
  ##Plots
  plot(uncond_snp$POS/1000,-log10(uncond_snp$P),xlab=paste("Chr ",unique(uncond_snp$`#CHROM`)," Position (Mb)",sep=""),ylab="-log10(P)",ylim=c(-minY,minP),pch=21,col="black",bg="green",cex = 0.7,yaxt="n", xaxt="n",yaxt="n", xaxt="n")
  points(uncond_aa$Nuc_Pos/1000,-log10(uncond_aa$P),pch=22,col="black",bg="green",cex = 1.4,yaxt="n", xaxt="n")
  points(uncond_alleles$POS/1000,-log10(uncond_alleles$P),pch=23,col="black",bg="green",cex = 1.3,yaxt="n", xaxt="n")
  points(cond_aa$Nuc_Pos/1000,-log10(cond_aa$P),pch=22,col="black",bg=alpha("red", 0.6),cex = 1.3,yaxt="n", xaxt="n")
  points(cond_snp$POS/1000,-log10(cond_snp$P),pch=21,col="black",bg=alpha("red", 0.6),cex = 0.7)
  points(cond_alleles$POS/1000,-log10(cond_alleles$P),pch=23,col="black",bg=alpha("red", 0.6),cex = 1.3,yaxt="n", xaxt="n")
  
  tix <- round(seq(from = 0,to = minP + 1,by = 1))
  axis(2,at=tix,las=1)
  
  # X-axis ticks
  tax <- seq(min(uncond_snp$POS/1000),max(uncond_snp$POS/1000),by = 100)
  tux <- round(tax / 1000,1)
  axis(1, at=tax, labels=tux)
  
  ## Add genes
  if(dim(genes)[1]>0){
    genes$Start <- genes$Start/1000
    genes$End <- genes$End/1000
    
    seq <- rep(x = c(minY - 0.3 ,(minY - 0.3 )/ 3 * 2,(minY - 0.3 )/ 3),length.out = dim(genes)[1])
    for (i in 1:dim(genes)[1]){
      text(x=((genes$Start[i]+genes$End[i])/2),y=-(seq[i]),labels=genes$Gene[i],cex=0.8,pos=1, offset = 0.35)
    }
  }
  
  # add legend
  top_allele <- uncond_alleles$ID[which.min(uncond_alleles$P)]
  top_allele_split <- strsplit(top_allele,split = '_')[[1]]
  top_allele <- paste0('HLA-',top_allele_split[1],'*',top_allele_split[2],':',top_allele_split[3])
  points(c(0,.Machine$double.xmax),c(bonf_thresh,bonf_thresh),type = 'l',lwd = 2,lty = 'dashed')
  uncond_aa$ID[uncond_aa$ID == 'A_99_I'] <- 'HLA-A T99I'
  uncond_aa$ID[uncond_aa$ID == 'A_170_Q'] <- 'HLA-A K170Q'
  uncond_aa$ID[uncond_aa$ID == 'A_35_T'] <- 'HLA-A Y35T'
  
  text(x = uncond_aa$Nuc_Pos[which.min(uncond_aa$P)]/1000 - 800, y = -log10(uncond_aa$P[which.min(uncond_aa$P)]),labels = uncond_aa$ID[which.min(uncond_aa$P)])
  points(x= c(uncond_aa$Nuc_Pos[which.min(uncond_aa$P)]/1000 - 450,uncond_aa$Nuc_Pos[which.min(uncond_aa$P)]/1000 - 50),y=c(-log10(uncond_aa$P[which.min(uncond_aa$P)]),-log10(uncond_aa$P[which.min(uncond_aa$P)])),type = 'l',lty = 'dashed')
  text(x = uncond_snp$POS[which.min(uncond_snp$P)]/1000 + 400, y = -log10(uncond_snp$P[which.min(uncond_snp$P)]) + 0.3,labels = uncond_snp$ID[which.min(uncond_snp$P)])
  points(x= c(uncond_snp$POS[which.min(uncond_snp$P)]/1000  + 100,uncond_snp$POS[which.min(uncond_snp$P)]/1000 + 30),y=c(-log10(uncond_snp$P[which.min(uncond_snp$P)]) + 0.3,-log10(uncond_snp$P[which.min(uncond_snp$P)]) + 0.1),type = 'l',lty = 'dashed')
  text(x = uncond_alleles$POS[which.min(uncond_alleles$P)]/1000 - 800, y = -log10(uncond_alleles$P[which.min(uncond_alleles$P)]) - 0.05,labels = top_allele)
  points(x= c(uncond_alleles$POS[which.min(uncond_alleles$P)]/1000 - 450,uncond_alleles$POS[which.min(uncond_alleles$P)]/1000 - 50),y=c(-log10(uncond_alleles$P[which.min(uncond_alleles$P)]) - 0.05,-log10(uncond_alleles$P[which.min(uncond_alleles$P)]) - 0.05),type = 'l',lty = 'dashed')
  if(is.na(cond_aa_pos)){
    legend("topright", legend=c("SNP",paste0("SNP: Conditional on ",uncond_snp$ID[which.min(uncond_snp$P)]),
                                "Amino Acid",paste0("Amino Acid: Conditional on ",uncond_aa$ID[which.min(uncond_aa$P)]),
                                "Allele",paste0("Allele: Conditional on ",top_allele)), col="black", pt.bg=c("green","red",'green','red','green','red'), pch=c(21,21,22,22,23,23), cex=c(1,1,1,1,1,1))
  }else{
    legend("topright", legend=c("SNP",paste0("SNP: Conditional on ",uncond_snp$ID[which.min(uncond_snp$P)]),
                                "Amino Acid",paste0("Amino Acid: Conditional on ",uncond_aa$ID[which.min(uncond_aa$P)]),
                                "Allele",paste0("Allele: Conditional on ",top_allele)), col="black", pt.bg=c("green","red",'green','red','green','red'), pch=c(21,21,22,22,23,23), cex=c(1,1,1,1,1,1))
  }
  return(list(snp_results=snp_results,allele_results=allele_results,aa_results=aa_results))
}


PlotFracAssoc <- function(snp_file,bonf_thresh,snp_to_label){
  snp_results <- data.table::fread(snp_file) %>% dplyr::select(CHR,POS,ID,REF,ALT,BETA=Estimate,SE=StdError,T_STAT=t,P=p.value) %>% dplyr::mutate(OR=exp(BETA),TEST='ADD') %>% dplyr::select(-BETA) %>%
                        dplyr::filter(`CHR`==6 &  POS >= 28477797 & POS <= 33448354)
  
  #Dowload GTF from:ftp://ftp.ensembl.org/pub/grch37/release-87/gtf/homo_sapiens/
  gtf <- rtracklayer::import('~/G2G-HBV/Homo_sapiens.GRCh37.87.gtf')
  gtf_df=as.data.frame(gtf)
  hla_gtf <- dplyr::filter(gtf_df,grepl('HLA',gene_name)) %>% dplyr::filter(seqnames == '6')
  
  genes <- data.frame(Gene = unique(hla_gtf$gene_name),Start = NA,End = NA,Strand=NA)
  for(i in 1:nrow(genes)){
    cur_gene <- dplyr::filter(hla_gtf,gene_name == genes$Gene[i] & type == 'gene')
    genes$Start[i] <- min(cur_gene$start)
    genes$End[i] <- max(cur_gene$end)
    genes$Strand[i] <- unique(as.character(cur_gene$strand))
    
  }
  genes <- dplyr::filter(genes, Gene %in% c('HLA-A','HLA-G','HLA-E','HLA-B','HLA-C','HLA-DRB1','HLA-DQB1','HLA-DQA1','HLA-DPA1','HLA-DPB1'))
  minP <- -log10(min(snp_results$P,na.rm = T)) + 1
  minY <- 2
  uncond_snp <- snp_results

  ##Plots
  plot(uncond_snp$POS/1000,-log10(uncond_snp$P),xlab=paste("Chr ",unique(uncond_snp$`#CHROM`)," Position (Mb)",sep=""),ylab="-log10(P)",ylim=c(-minY,minP),pch=21,cex = 0.7,yaxt="n", xaxt="n",yaxt="n", xaxt="n")
  points(c(0,.Machine$double.xmax),c(bonf_thresh,bonf_thresh),type = 'l',lwd = 2,lty = 'dashed')
  text(x = uncond_snp$POS[which.min(uncond_snp$P)]/1000 + 400, y = -log10(uncond_snp$P[which.min(uncond_snp$P)]) + 0.3,labels = uncond_snp$ID[which.min(uncond_snp$P)])
  points(x= c(uncond_snp$POS[which.min(uncond_snp$P)]/1000  + 100,uncond_snp$POS[which.min(uncond_snp$P)]/1000 + 30),y=c(-log10(uncond_snp$P[which.min(uncond_snp$P)]) + 0.3,-log10(uncond_snp$P[which.min(uncond_snp$P)]) + 0.1),type = 'l',lty = 'dashed')
  
  text(x = uncond_snp$POS[which(uncond_snp$ID == snp_to_label)]/1000 + 400, y = -log10(uncond_snp$P[which(uncond_snp$ID == snp_to_label)]) + 0.3,labels = uncond_snp$ID[which(uncond_snp$ID == snp_to_label)])
  points(x= c(uncond_snp$POS[which(uncond_snp$ID == snp_to_label)]/1000  + 100,uncond_snp$POS[which(uncond_snp$ID == snp_to_label)]/1000 + 30),y=c(-log10(uncond_snp$P[which(uncond_snp$ID == snp_to_label)]) + 0.3,-log10(uncond_snp$P[which(uncond_snp$ID == snp_to_label)]) + 0.1),type = 'l',lty = 'dashed')
  
  tix <- round(seq(from = 0,to = minP + 1,by = 1))
  axis(2,at=tix,las=1)
  
  # X-axis ticks
  tax <- seq(min(uncond_snp$POS/1000),max(uncond_snp$POS/1000),by = 100)
  tux <- round(tax / 1000,1)
  axis(1, at=tax, labels=tux)
  
  ## Add genes
  if(dim(genes)[1]>0){
    genes$Start <- genes$Start/1000
    genes$End <- genes$End/1000
    
    seq <- rep(x = c(minY - 0.3 ,(minY - 0.3 )/ 3 * 2,(minY - 0.3 )/ 3),length.out = dim(genes)[1])
    for (i in 1:dim(genes)[1]){
      text(x=((genes$Start[i]+genes$End[i])/2),y=-(seq[i]),labels=genes$Gene[i],cex=0.8,pos=1, offset = 0.35)
    }
  }
}

source('./scripts/misc/N_Eff.R')
pdf(file='~/G2G-HBV/Final_Figures/HLA_G2G_Merged.pdf',width = 8.5,height = 11)
par(mfrow=c(3,1))
PC_C_pos_0160_A_hla_results <- PlotHLAAssoc(list(uncond = '~/G2G-HBV/data/results_POP_asian_GT_A_C_D_B_TYPE_pPCA_LOCO_FALSE/SAIGE_SPA_TEST_LOCO_FALSE_gene_PC_C_pos_0160_A',cond = '~/G2G-HBV/PyHLA_Out/Asian/hla_snps_gene_PC_C_pos_0160_A_cond_snp_rs41541212.PHENO1.glm.logistic'),
             '~/G2G-HBV/PyHLA_Out/Asian/omnibus/gene_PC_C_pos_0160_A/',
             list(uncond = '~/G2G-HBV/PyHLA_Out/Asian/hla_allele_assocAllele_HBV_gene_PC_C_pos_0160_A',cond = '~/G2G-HBV/PyHLA_Out/Asian/hla_allele_gene_PC_C_pos_0160_A_assocAllele_cond_snp_A_33_03.PHENO1.glm.logistic'),omnibus = T,bonf_thresh = -log10(5e-8 / N_Eff))
fig_label('A)',cex = 2)

Pol_pos_0049_N_hla_results <- PlotHLAAssoc(list(uncond = '~/G2G-HBV/data/results_POP_asian_GT_A_C_D_B_TYPE_pPCA_LOCO_FALSE/SAIGE_SPA_TEST_LOCO_FALSE_gene_Pol_pos_0049_N',cond = '~/G2G-HBV/PyHLA_Out/Asian/hla_snps_gene_Pol_pos_0049_N_cond_snp_rs3129697.PHENO1.glm.logistic'),
              '~/G2G-HBV/PyHLA_Out/Asian/omnibus/gene_Pol_pos_0049_N/',
              list(uncond = '~/G2G-HBV/PyHLA_Out/Asian/hla_allele_assocAllele_HBV_gene_Pol_pos_0049_N',cond = '~/G2G-HBV/PyHLA_Out/Asian/hla_allele_gene_Pol_pos_0049_N_assocAllele_cond_snp_A_02_06.PHENO1.glm.logistic'),omnibus = T,bonf_thresh = -log10(5e-8 / N_Eff))
fig_label('B)',cex = 2)

PC_C_pos_0067_Y_eur_hla_results <- PlotHLAAssoc(list(uncond = '~/G2G-HBV/data/results_POP_european_GT_A_C_D_F_H_TYPE_pPCA_LOCO_FALSE/SAIGE_SPA_TEST_LOCO_FALSE_gene_PC_C_pos_0067_Y',cond = '~/G2G-HBV/PyHLA_Out/European/hla_snps_gene_PC_C_pos_0067_Y.PHENO1.glm.logistic'),
                                                                    '~/G2G-HBV/PyHLA_Out/European/omnibus/gene_PC_C_pos_0067_Y/',
                                                                    list(uncond = '~/G2G-HBV/PyHLA_Out/European/hla_allele_assocAllele_HBV_gene_PC_C_pos_0067_Y',cond = '~/G2G-HBV/PyHLA_Out/European/hla_allele_gene_PC_C_pos_0067_Y_assocAllele_cond_snp_A_01_01.PHENO1.glm.logistic'),omnibus = T,bonf_thresh = -log10(5e-8 /N_Eff))
fig_label('C)',cex = 2)
dev.off()

saveRDS(list(PC_C_pos_0160_A_hla_results=PC_C_pos_0160_A_hla_results,
             Pol_pos_0049_N_hla_results=Pol_pos_0049_N_hla_results,
             PC_C_pos_0067_Y_eur_hla_results=PC_C_pos_0067_Y_eur_hla_results),file = '~/G2G-HBV/PyHLA_Out/HLA_result.rds')