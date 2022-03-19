library(phytools)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggtree)
library(ape)
library(viridis)
library(scales)
source('~/G2G-HBV/src/functions-data-pathogen.R')
#load('~/G2G-HBV/data/results_POP_asian_GT_A_C_D_B_TYPE_pPCA_LOCO_FALSE/prep-data.rda')

PlotPCA <- function(data_path,gg_title,tree_file,aa_table){
  load(data_path)
  
  # tree <- phytools::read.newick(FILE_PATHOGEN_TREE)
  # host_id <- paste0(dat_covars_raw$host_id,"_",dat_covars_raw$GT)
  # names(host_id) <- dat_covars_raw$pathogen_id
  # tree_labl <- match(tree$tip.label,names(host_id))
  # new_tree_labl <- host_id[tree_labl]
  # new_tree_labl[which(tree$tip.label=='NC_028129')] <- 'Wooly Monkey HBV'
  # names(new_tree_labl) <- NULL
  # tree$tip.label <- new_tree_labl
  # pdf('~/G2G-HBV/hbv-phylo/tree.id.gt.pdf',height = 150,width = 10)
  # phytools::plotTree(tree)
  # dev.off()

  # GT <- dat_covars_raw$GT
  # names(GT) <- dat_covars_raw$pathogen_id
  # tree_labl <- match(tree$tip.label,names(GT))
  # new_tree_labl <- GT[tree_labl]
  # new_tree_labl[which(tree$tip.label=='NC_028129')] <- 'Wooly Monkey HBV'
  # names(new_tree_labl) <- NULL
  # tree$tip.label <- new_tree_labl
  # tree$tip.label[is.na(tree$tip.label)] <- 'None'
  # for(i in 1:length(tree$tip.label)){
  #   tree<-paintBranches(tree,edge=i,state=tree$tip.label[i])
  # }
  # 
  # cols<-setNames(c("blue","red",'green1','purple','yellow3','cyan','green3','magenta','red4','black','grey','black'),c("A","B","C","D","E","F",'G','H','Mixed','Wooly Monkey HBV','None',1))
  # tree$tip.label <- rep('',length(tree$tip.label))
  # plotSimmap(tree,type="fan",colors=cols)
  # add.simmap.legend(colors = cols,prompt = T)
  
  dat_pathogen <- dat_pathogen %>% dplyr::filter(ID %in% ids_order$X1)
  pPC <- pathogen_pca(
    path_tree = tree_file,
    path_pathogen = aa_table,
    id.list = idlist_pop_path,
    pca = FALSE,
    ppca = TRUE,
    n.pc = n_ppc,
    loadings = FALSE
  )
  pPC_eigen <- pathogen_pca(
    path_tree = tree_file,
    path_pathogen = aa_table,
    id.list = idlist_pop_path,
    eigenvalues = T,
    pca = FALSE,
    ppca = TRUE,
    n.pc = n_ppc,
    loadings = FALSE
  )
  
  dat_covar_pPC <- dplyr::inner_join(dat_covars_raw,pPC,by=c("pathogen_id"="ID"))
  pPC_data <- dat_covar_pPC %>% dplyr::filter(host_id %in% ids_order$X1)
  pathogen_pPCA_1_2 <- ggplot2::ggplot(data = pPC_data) + aes(x=pPC1,y=pPC2,colour = GT,label=host_id) + geom_text() + ggtitle(paste0(gg_title,':pPCA'))
  pathogen_pPCA_3_4 <- ggplot2::ggplot(data = pPC_data) + aes(x=pPC3,y=pPC4,colour = GT,label=host_id) + geom_text() + ggtitle(paste0(gg_title,':pPCA'))

  pca <- pathogen_pca(
    path_tree = tree_file,
    path_pathogen = aa_table,
    id.list = idlist_pop_path,
    pca = TRUE,
    ppca = FALSE,
    n.pc = n_ppc,
    loadings = FALSE)
  PC_data <- dplyr::inner_join(dat_covars_raw,data.frame(ID = pca$ID,PC1=pPC$pPC1,PC2 = pPC$pPC2,PC3 = pPC$pPC3,PC4 = pPC$pPC4),by=c("pathogen_id"="ID")) %>% dplyr::filter(host_id %in% ids_order$X1)
  pathogen_PCA_1_2 <- ggplot2::ggplot(data = PC_data) + aes(x=PC1,y=PC2,colour = GT,label=host_id) + geom_text()+ ggtitle(paste0(gg_title,':PCA'))
  pathogen_PCA_3_4 <- ggplot2::ggplot(data = PC_data) + aes(x=PC3,y=PC4,colour = GT,label=host_id) + geom_text()+ ggtitle(paste0(gg_title,':PCA'))

  return(list(plot=ggpubr::ggarrange(pathogen_pPCA_1_2,pathogen_pPCA_3_4,pathogen_PCA_1_2,pathogen_PCA_3_4,ncol = 2,nrow = 2),pPC=PC_data,pPC_eigen=pPC_eigen))
}
PlotTree <- function(tree_file,outgroup="NC_028129.1",outgrop_label="Woolly Monkey HBV",raw=F){
  #Get Metadata of all samples
  load('~/G2G-HBV/data/results_POP_asian_GT_A_C_D_B_TYPE_pPCA_LOCO_FALSE/prep-data.rda')
  asn_samples <- dplyr::filter(dat_covars_raw,host_id %in% dat_ids$host_id)
  load('~/G2G-HBV/data/results_POP_european_GT_A_C_D_F_H_TYPE_pPCA_LOCO_FALSE/prep-data.rda')
  eur_samples <- dplyr::filter(dat_covars_raw,host_id %in% dat_ids$host_id)
  all_samples <- rbind(asn_samples,eur_samples) %>% dplyr::select(ID=pathogen_id,GT)
  all_samples <- rbind(all_samples,data.frame(ID=outgroup,GT=outgrop_label))
  #Load Tree
  tree <- ape::read.tree(tree_file)
  #Re-root tree
  #tree <- ape::root(tree,outgroup = "NC_028129.1",edgelabel=T)
  if(raw){
    #Remove excluded samples
    tree_filt <- ape::drop.tip(tree,setdiff(tree$tip.label,c(all_samples$ID,outgroup)))
    
    # Plot with ggtree
    p1 <- ggtree(tree_filt) 
    p1 <- p1 %<+% all_samples + 
      geom_tippoint(aes(color = factor( GT )))
    return(p1)
    
  }else{
    #Remove excluded samples
    tree_filt <- ape::drop.tip(tree,setdiff(tree$tip.label,c(all_samples$ID,outgroup)))
    #Only labels nodes with depth less than 10
    node_depth <- node.depth(tree_filt,method = 2)
    internal_node_depth <- node_depth[which(node_depth!= 1)]
    tree_filt$node.label[internal_node_depth < 23 | as.numeric(tree_filt$node.label) < 70] <- ""
    # Plot with ggtree
    p1 <- ggtree(tree_filt) 
    p1 <- p1 %<+% all_samples + 
      geom_tippoint(aes(color = factor( GT )))
    
    p1$data$label[grepl(x=p1$data$label,'GS') | grepl(x=p1$data$label,'NC')] <- ''
    col <- hue_pal()(6)
    p1 <- p1 + ggrepel::geom_text_repel(aes(label=label),direction = 'both') + labs(color = "HBV Genotype") + scale_color_manual(values=c(col,'gray0'))
    return(p1)
    
  }
}
VisualizeG2G <- function(tree_file,data_path,AA_Variant,Host_SNP,ancestry = 'Asian',is_HLA=F,Host_Label,Pathogen_Label,outgroup="NC_028129.1"){
  load(data_path)
  if(is_HLA){
    #Load HLA_Allele Dosage
    host_dosage_raw <- data.table::fread(glue::glue('~/G2G-HBV/PyHLA_Out/{ancestry}/hla_allele_assocAllele.dosage'))
    host_dosage <- as.vector(t(host_dosage_raw[,..Host_SNP,drop = T]))
    host_dosage[host_dosage==0] <- 'Non-Carrier'
    host_dosage[host_dosage==1 | host_dosage==2] <- 'Carrier'
    names(host_dosage) <- host_dosage_raw$IID
  }else{
    bed_path <- gsub(data_path,pattern = 'results_',replacement = 'processed_')
    bed_path <- gsub(bed_path,pattern = 'prep-data.rda',replacement = 'hbv_gilead_QC')
    host_genotype <- snpStats::read.plink(bed = bed_path,select.snps = Host_SNP)
    host_dosage_raw <- as(host_genotype$genotypes, Class = 'numeric')
    host_dosage <- as.vector(host_dosage_raw)  
    names(host_dosage) <- rownames(host_dosage_raw)
    host_dosage[host_dosage == 2] <- paste0(host_genotype$map$allele.2,host_genotype$map$allele.2)
    host_dosage[host_dosage == 1] <- paste0(host_genotype$map$allele.2,host_genotype$map$allele.1)
    host_dosage[host_dosage == 0] <- paste0(host_genotype$map$allele.1,host_genotype$map$allele.1)
  }

  cluster <- dat_ids %>% dplyr::select(host_id = host_id)
  dat_pathogen_raw <- as.data.frame(dat_pathogen_raw[match(cluster$host_id,dat_pathogen_raw$ID),-1])
  rownames(dat_pathogen_raw) <- cluster$host_id
  pathogen_dosage <- dat_pathogen_raw[,AA_Variant,drop = T]
  pathogen_dosage <- pathogen_dosage[match(dat_ids$host_id,rownames(dat_pathogen_raw))]
  names(pathogen_dosage) <- dat_ids$pathogen_id
  pathogen_dosage[pathogen_dosage==0] <- 'Absent'
  pathogen_dosage[pathogen_dosage==1] <- 'Present'
  
  host_dosage <- host_dosage[match(dat_ids$host_id,names(host_dosage))]
  names(host_dosage) <- dat_ids$pathogen_id
  df <- data.frame(ID = names(host_dosage),host_dosage = host_dosage,pathogen_dosage = pathogen_dosage,GT = dat_ids_pop$GT[match(names(host_dosage),dat_ids_pop$pathogen_id)])
  
  GT <- unique(df$GT)
  p1 <- vector(mode = 'list',length = length(GT))
  for(i in 1:length(GT)){
    cur_GT <- GT[i]
    cur_df <- df %>% dplyr::filter(GT == cur_GT)
    #Load Tree
    tree <- ape::read.tree(tree_file)
    #Remove excluded samples
    tree_filt <- ape::drop.tip(tree,c(cur_df$ID[is.na(df$host_dosage)],cur_df$ID[is.na(cur_df$pathogen_dosage)],setdiff(tree$tip.label,cur_df$ID)))
    p1[[i]] <- ggtree(tree_filt) %<+% cur_df + 
      geom_tippoint(aes(color = factor(pathogen_dosage),shape = factor(host_dosage))) + scale_size_manual(values = 3) + labs(shape = Host_Label,color = Pathogen_Label) + ggtitle(paste0('Genotype:',cur_GT))
  }
  return(p1)
}
# p1 <- PlotTree('~/G2G-HBV/data/phylo_tree/HBV_WG+og_rr.nw',raw=T)
PC_C_pos_160_A <- VisualizeG2G('~/G2G-HBV/data/phylo_tree/HBV_WG+og_rr.nw','~/G2G-HBV/data/results_POP_asian_GT_A_C_D_B_TYPE_pPCA_LOCO_FALSE/prep-data.rda',AA_Variant = 'gene_PC_C_pos_0160_A',Host_SNP = 'A_33_03',Host_Label = 'HLA-A*33:03',Pathogen_Label = 'PreCore/Core: Position 160\nResidue A',is_HLA = T,ancestry = 'Asian')
ggpubr::ggarrange(plotlist = PC_C_pos_160_A[1:2],ncol = 2,common.legend = T)

Pol_pos_49_N <- VisualizeG2G('~/G2G-HBV/data/phylo_tree/HBV_WG+og_rr.nw','~/G2G-HBV/data/results_POP_asian_GT_A_C_D_B_TYPE_pPCA_LOCO_FALSE/prep-data.rda',AA_Variant = 'gene_Pol_pos_0049_N',Host_SNP = 'A_02_06',Host_Label = 'HLA-A*02:06',Pathogen_Label = 'Pol: Position 49\nResidue N',is_HLA = T,ancestry = 'Asian')
ggpubr::ggarrange(plotlist = Pol_pos_49_N[1:2],ncol=2,common.legend = T)

PC_C_pos_67_Y <- VisualizeG2G('~/G2G-HBV/data/phylo_tree/HBV_WG+og_rr.nw','~/G2G-HBV/data/results_POP_european_GT_A_C_D_F_H_TYPE_pPCA_LOCO_FALSE/prep-data.rda',AA_Variant = 'gene_PC_C_pos_0067_Y',Host_SNP = 'A_01_01',Host_Label = 'HLA-A*01:01',Pathogen_Label = 'PreCore/Core: Position 67\nResidue Y',is_HLA = T,ancestry = 'European')
ggpubr::ggarrange(plotlist = PC_C_pos_67_Y[1:2],ncol=2,common.legend = T)

S_pos_0017_A <- VisualizeG2G('~/G2G-HBV/data/phylo_tree/HBV_WG+og_rr.nw','/home/zmxu/G2G-HBV/data/results_POP_asian_GT_A_C_D_B_TYPE_pPCA_LOCO_FALSE/prep-data.rda',AA_Variant = 'gene_S_pos_0017_A',Host_SNP = 'rs2296651',Host_Label = 'rs2296651',Pathogen_Label = 'PreS1: Position 17\nResidue A',is_HLA = F,ancestry = 'Asian')
ggpubr::ggarrange(plotlist = S_pos_0017_A[1:2],ncol=2,common.legend = T)
