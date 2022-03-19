library(dplyr,lib.loc = '~/R/Old/')
library(dbplyr,lib.loc = '~/R/Old/')
library(GenomicFeatures)
library(pbmcapply)
library(seqinr)
library(circlize)

asian_all_subtype_results_path <- '~/G2G-HBV/data/results_POP_asian_GT_A_C_D_B_TYPE_pPCA_LOCO_FALSE/'
eur_all_subtype_results_path <- '~/G2G-HBV/data/results_POP_european_GT_A_C_D_F_H_TYPE_pPCA_LOCO_FALSE/'

LOCO = 'FALSE'
sage_files_asn <- dir(asian_all_subtype_results_path)
sage_files_asn <- sage_files_asn[grepl(glue::glue("SAIGE_SPA_TEST_LOCO_{as.character(LOCO)}"),sage_files_asn)]
sage_genes_asn <- sapply(sage_files_asn,function(x) strsplit(x,split = glue::glue('SAIGE_SPA_TEST_LOCO_{as.character(LOCO)}_'))[[1]][2])

sage_files_eur <- dir(eur_all_subtype_results_path)
sage_files_eur <- sage_files_eur[grepl(glue::glue("SAIGE_SPA_TEST_LOCO_{as.character(LOCO)}"),sage_files_eur)]
sage_genes_eur <- sapply(sage_files_eur,function(x) strsplit(x,split = glue::glue('SAIGE_SPA_TEST_LOCO_{as.character(LOCO)}_'))[[1]][2])

asian_g2g_all_subtype_res <- pbmclapply(sage_files_asn,function(x) data.table::fread(paste0(asian_all_subtype_results_path,x)) %>% dplyr::filter(p.value < 5e-8) %>% dplyr::select(CHR,POS,p.value),mc.cores = 4)
names(asian_g2g_all_subtype_res) <- sage_genes_asn
eur_g2g_all_subtype_res <- pbmclapply(sage_files_eur,function(x) data.table::fread(paste0(eur_all_subtype_results_path,x)) %>% dplyr::filter(p.value < 5e-8) %>% dplyr::select(CHR,POS,p.value),mc.cores = 4)
names(eur_g2g_all_subtype_res) <- sage_genes_eur
sage_genes <- c(sage_genes_asn,sage_genes_eur)
g2g_all_subtype_res <- c(asian_g2g_all_subtype_res,eur_g2g_all_subtype_res)
  
gene_pos <- as.numeric(sapply(sage_genes,function(x) strsplit(strsplit(x,split = 'pos_')[[1]][2],split = '_')[[1]][1]))
gene_name <- sapply(sage_genes,function(x) strsplit(strsplit(x,split = 'gene_')[[1]][2],split = '_')[[1]][1])

#Get human chromosome information
all_human_chr <- sort(unique(do.call(rbind,asian_g2g_all_subtype_res)$CHR))
human_chr_length <- getChromInfoFromBiomart(biomart="ENSEMBL_MART_ENSEMBL",
                        dataset="hsapiens_gene_ensembl",
                        id_prefix="ensembl_",
                        host="grch37.ensembl.org",
                        port=80) %>% dplyr::filter(chrom %in% as.character(all_human_chr))
human_chr_length$chrom <- paste0('Chr ',human_chr_length$chrom)
# Create dummy data for human genome
human_chr_to_incl <- c('Chr 6','Chr 14')
human_data = data.frame(
  factor = c(rbind(as.character(human_chr_length$chrom),as.character(human_chr_length$chrom))),
  x = c(rbind(rep(0,nrow(human_chr_length)),human_chr_length$length)),
  y=rep(c(0,0.1),nrow(human_chr_length)),stringsAsFactors = F
)
human_data <- dplyr::filter(human_data,as.character(factor) %in% human_chr_to_incl)

# Create dummy data for pathogen (based on nucelotide sequences)
scale_factor_patho <- 2e5
pathogen_data = data.frame(
  factor = 'HBV',
  x = c(0,3215) * scale_factor_patho,
  y=c(0,0.1),stringsAsFactors = F
)


# Create dummy data for pathogen genes
scale_factor_patho <- 2e5
pathogen_gene_data = data.frame(
  factor = c(rep('Pol',2),rep("PreS1",2),rep("PreS2",2),rep("S",2),rep("X",4),rep("PreC/C",2)),
  x = c(2307 - 1814,
        1623 + 1401,
        2848 - 1814 + 3,
        2848 - 1814 + 3 + 119*3,
        2848 - 1814 + 3 + 119*3 + 1,
        2848 - 1814 + 3 + 119*3 + 1 + (174 - 120 + 1)*3,
        2848 - 1814 + 3 + 119*3 + 1 + (174 - 120 + 1)*3 + 1,
        2848 - 1814 + 3 + 119*3 + 1 + (174 - 120)*3 + 1 + (400 - (118) - ((174 - 120 + 1)))*3 + 1,
        1401 + 1374,
        3215,
        0,
        1838 - 1814,
        0,
        638) * scale_factor_patho,
  col = c(rep('green',2),rep('skyblue1',2),rep('dodgerblue',2),rep('blue',2),rep('firebrick',4),rep('orange',2)),
  y=c(0,0.1),stringsAsFactors = F
)

merged_data <- rbind(human_data,pathogen_data)
merged_factors <- unique(c(human_data$factor,pathogen_data$factor))

#Get Sig Interactions
g2g_all_subtype_res_recode <- lapply(g2g_all_subtype_res,function(x) x %>% dplyr::mutate(nuid = paste0('Chr ',CHR,':',POS)) %>% dplyr::select(nuid,p.value))
human_variants_nuid <- unique(do.call(rbind,g2g_all_subtype_res_recode)$nuid)
pathogen_variants <- names(g2g_all_subtype_res_recode)
pathogen_variants_gene <- sapply(pathogen_variants,function(x) strsplit(split = '_pos',strsplit(x,split = 'gene_')[[1]][2])[[1]][1])
pathogen_variants_pos <- as.numeric(sapply(pathogen_variants,function(x) strsplit(split = '_',strsplit(x,split = 'pos_')[[1]][2])[[1]][1]))
pathogen_variants_nuid <- mapply(function(x,y){
 if(x=='PC_C'){
    return(paste0('PreC/C',':',y))
 }else if(x == 'S'){
   return(paste0('PreS1',':',y))
   }else{
   return(paste0(x,':',y))
  }
},pathogen_variants_gene,pathogen_variants_pos)

interact_mat <- matrix(Inf,ncol = length(human_variants_nuid),nrow = length(pathogen_variants_nuid))
rownames(interact_mat) <- pathogen_variants
colnames(interact_mat) <- human_variants_nuid
for(i in 1:nrow(interact_mat)){
  print(i)
  cur_pathogen_variant <- rownames(interact_mat)[i]
  cur_int <- g2g_all_subtype_res_recode[[cur_pathogen_variant]]
  interact_mat[cur_pathogen_variant,cur_int$nuid] <- cur_int$p.value
}


bonf_threshold <- 5e-8 / 523
gwas_thresh <- 5e-8
sig_interactions <- as.data.frame(which(interact_mat < gwas_thresh,arr.ind = T))
sig_interactions_df <- data.frame(chr_name = character(),chr_pos = numeric(),hbv_pos = numeric(),patho_gene = character())
for(i in 1:nrow(sig_interactions)){
  cur_row <- sig_interactions$row[i]
  cur_col <- sig_interactions$col[i]
  patho_nuid <- pathogen_variants_nuid[cur_row]
  human_nuid <- colnames(interact_mat)[cur_col]
  
  patho_gene <- strsplit(patho_nuid,':')[[1]][1]
  patho_pos <- as.numeric(strsplit(patho_nuid,':')[[1]][2])
  human_chr <- strsplit(human_nuid,':')[[1]][1]
  human_pos <- as.numeric(strsplit(human_nuid,':')[[1]][2])
  sig_interactions_df <- rbind(sig_interactions_df,data.frame(chr_name = human_chr,chr_pos = human_pos,hbv_pos = patho_pos,patho_gene = patho_gene, p = interact_mat[cur_row,cur_col]))
  # circos.link(sector.index1 = human_chr, point1 = human_pos, sector.index2 = "HBV", point2 = patho_pos,col = 'red')
}
sig_interactions_df <- as.data.frame(sig_interactions_df %>% dplyr::group_by(patho_gene,hbv_pos)  %>% dplyr::arrange(p) %>% dplyr::filter(row_number() == 1))
sig_interactions_df <- sig_interactions_df %>% dplyr::filter(!(patho_gene == 'PreS1' & hbv_pos == 17 | hbv_pos == 32))

# Initialize the plot.
groups <- list(c('X'),c('PreC/C'),c('Pol'),c('PreS1','PreS2','S'))
circos.par(gap.after = c(2,20,20))
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
merged_factors <- factor(c(rbind(as.character(merged_factors),as.character(merged_factors))),levels = merged_factors,ordered = T)
circos.initialize(factors = merged_factors, x = merged_data$x )

for(i in 1:length(groups)){
  cur_group <- groups[[i]]
  if(i == 1){
    cur_regions <- dplyr::filter(pathogen_gene_data,factor %in% cur_group)
    cur_col <- c(rep('grey',length(unique(human_data$factor))),rgb(255,255,255,alpha = 100,maxColorValue = 255))
    circos.trackPlotRegion(factors = merged_factors, ylim = c(0,1), bg.col = cur_col,track.margin = c(0,0.05),
                           bg.border = NA, track.height = 0.05,
                           panel.fun = function(x, y) {
                             sector.index = get.cell.meta.data("sector.index")
                             name = get.cell.meta.data("sector.index")
                             i = get.cell.meta.data("sector.numeric.index")
                             xlim = get.cell.meta.data("xlim")
                             ylim = get.cell.meta.data("ylim")
                             
                             # #text direction (dd) and adjusmtents (aa)
                             # theta = circlize(mean(xlim), 1.3)[1, 1] %% 360
                             # dd <- ifelse(theta < 90 || theta > 270, "clockwise", "reverse.clockwise")
                             # aa = c(1, 0.5)
                             # if(theta < 90 || theta > 270)  aa = c(0, 0.5)
                             # 
                             # #plot labels
                             # circos.text(x=cur_regions$x[1], y=3.5, labels=name, facing = 'bending.outside', cex=1,niceFacing = T)
                             
                           })
    
    j=1
    while(j <= nrow(cur_regions)){
      circos.rect(cur_regions[j,]$x, 0, cur_regions[j+1,]$x,1, sector.index = unique(merged_factors)[3], col = cur_regions[j,]$col)
      major_tick_bp <- 30
      range <- c(cur_regions[j,]$x,cur_regions[j+1,]$x)
      n_ticks <- floor(max(c(2,(range[2] - range[1]) / scale_factor_patho / 3 / major_tick_bp)))
      ticks_to_plot <- seq(from = range[1],to = range[2],length.out = n_ticks)
      labels_to_plot <- round(seq(from = 0,to = round((range[2] - range[1])/scale_factor_patho/3),length.out = n_ticks),0)
      circos.axis(h = 'top',sector.index = unique(merged_factors)[3],labels.cex = 0.6,lwd = 0.5,labels.away.percentage = 1,labels = labels_to_plot,major.at = ticks_to_plot,major.tick.percentage = 0.1)
      
      cur_region_assoc <- dplyr::filter(sig_interactions_df,patho_gene == cur_regions[j,]$factor)
      if(nrow(cur_region_assoc) > 0){
        for(q in 1:nrow(cur_region_assoc)){
          circos.link(sector.index1 = cur_region_assoc$chr_name[q], point1 = cur_region_assoc$chr_pos[q], sector.index2 = "HBV", point2 = range[1] + cur_region_assoc$hbv_pos[q] * 3 * scale_factor_patho,col = 'red',rou1 = 0.9)
        }
      }
      j <- j + 2
    }
    
    
    major_tick_bp <- 2.5e7
    for(k in 1:length(unique(human_data$factor))){
      cur_factor <- unique(human_data$factor)[k]
      range <- max(dplyr::filter(human_data,factor == cur_factor)$x)
      n_ticks <- floor(range / major_tick_bp)
      ticks_to_plot <- seq(from = 0,to = range,length.out = n_ticks)
      circos.axis(h = 'top',sector.index = cur_factor,labels.cex = 0.6,lwd = 0.5,labels.away.percentage = 1,labels = c(paste0(round(ticks_to_plot[1] / 1e6),'Mb'),round(ticks_to_plot[2:(length(ticks_to_plot))] / 1e6)),major.at = ticks_to_plot,major.tick.percentage = 0.1)
    }
    
    
  }else{
    cur_col <- c(rep(rgb(255,255,255,alpha = 100,maxColorValue = 255),length(unique(human_data$factor))),rgb(255,255,255,alpha = 100,maxColorValue = 255))
    cur_regions <- dplyr::filter(pathogen_gene_data,factor %in% cur_group)
    circos.trackPlotRegion(factors = unique(merged_factors), ylim = c(0,1), bg.col = cur_col,track.margin = c(0,0.1),
                           bg.border = NA, track.height = 0.05,
                           panel.fun = function(x, y) {
                             sector.index = get.cell.meta.data("sector.index")
                             name = get.cell.meta.data("sector.index")
                             i = get.cell.meta.data("sector.numeric.index")
                             xlim = get.cell.meta.data("xlim")
                             ylim = get.cell.meta.data("ylim")
                             # #text direction (dd) and adjusmtents (aa)
                             # theta = circlize(mean(xlim), 1.3)[1, 1] %% 360
                             # dd <- ifelse(theta < 90 || theta > 270, "clockwise", "reverse.clockwise")
                             # aa = c(1, 0.5)
                             # if(theta < 90 || theta > 270)  aa = c(0, 0.5)
                             # 
                             # #plot labels
                             # circos.text(x=cur_regions$x[1], y=3.5, labels=name, facing = 'outside', cex=1,niceFacing = T)
                             
                           })
    j=1
    while(j <= nrow(cur_regions)){
      circos.rect(cur_regions[j,]$x, 0, cur_regions[j+1,]$x,1, sector.index = unique(merged_factors)[3], col = cur_regions[j,]$col)
      major_tick_bp <- 30
      range <- c(cur_regions[j,]$x,cur_regions[j+1,]$x)
      n_ticks <- floor(max(c(2,(range[2] - range[1]) / scale_factor_patho / 3 / major_tick_bp)))
      ticks_to_plot <- seq(from = range[1],to = range[2],length.out = n_ticks)
      labels_to_plot <- round(seq(from = 0,to = round((range[2] - range[1])/scale_factor_patho/3),length.out = n_ticks),0)
      circos.axis(h = 'top',sector.index = unique(merged_factors)[3],labels.cex = 0.6,lwd = 0.5,labels.away.percentage = 1,labels = labels_to_plot,major.at = ticks_to_plot,major.tick.percentage = 0.1)
      cur_region_assoc <- dplyr::filter(sig_interactions_df,patho_gene == cur_regions[j,]$factor)
      if(nrow(cur_region_assoc) > 0){
        for(q in 1:nrow(cur_region_assoc)){
          if(cur_region_assoc$p[q] < bonf_threshold){
            circos.link(sector.index1 = cur_region_assoc$chr_name[q], point1 = cur_region_assoc$chr_pos[q], sector.index2 = "HBV", point2 = range[1] + cur_region_assoc$hbv_pos[q] * 3 * scale_factor_patho,col = 'red',rou1 = 0.9)
            
          }else{
            if(cur_region_assoc$hbv_pos[q] == 67){
              circos.link(sector.index1 = cur_region_assoc$chr_name[q], point1 = cur_region_assoc$chr_pos[q], sector.index2 = "HBV", point2 = range[1] + cur_region_assoc$hbv_pos[q] * 3 * scale_factor_patho,col = 'blue',rou1 = 0.9,lty = 'dashed')
              
            }else{
              circos.link(sector.index1 = cur_region_assoc$chr_name[q], point1 = cur_region_assoc$chr_pos[q], sector.index2 = "HBV", point2 = range[1] + cur_region_assoc$hbv_pos[q] * 3 * scale_factor_patho,col = 'red',rou1 = 0.9,lty = 'dashed')
              
            }

                        
          }

        }
      }
      
      j <- j + 2
      
    }
    
  }
  
}
library(latex2exp)
legend("topleft", legend = c(TeX('East Asian (p < $9.6$ x $10^{-11}$)'),TeX('East Asian (p < $5$ x $10^{-8}$)'),TeX('European (p < $5$ x $10^{-8}$)')), col = c('red','red','blue'),
       ncol = 1, cex = 1, lwd = 1.3,lty = c('solid','dashed','dashed'))
