library(dplyr)
library(data.table)

## Load data
load('~/G2G-HBV/data/results_POP_asian_GT_A_C_D_B_TYPE_pPCA_LOCO_FALSE/setup.rda')
load('~/G2G-HBV/data/results_POP_asian_GT_A_C_D_B_TYPE_pPCA_LOCO_FALSE/prep-data.rda')

## Table 1 - Baseline Characteristics
tmp_cluster$GT[tmp_cluster$GT == 'Mixed genotype detected'] <- 'Mixed'
tmp_cluster$SEX[tmp_cluster$SEX == 'F'] <- 'Female'
tmp_cluster$SEX[tmp_cluster$SEX == 'M'] <- 'Male'

eur_cluster <- dplyr::filter(tmp_cluster,gr_cluster==2 & RACE == 'WHITE')
eur_unrelated <- data.table::fread('~/G2G-HBV/data/processed_POP_european_GT_A_C_D_F_H_TYPE_pPCA_LOCO_FALSE/kingunrelated.txt',header = F)
eur_cluster <- dplyr::filter(eur_cluster,host_id %in% eur_unrelated$V1)
eur_cluster_age <- glue("{median(eur_cluster$AGE)} ({paste(range(eur_cluster$AGE),collapse = '-')})")
names(eur_cluster_age) <- 'Median (range)'
eur_cluster_GT <- table(eur_cluster$GT)
eur_cluster_sex <- table(eur_cluster$SEX)
eur_cluster_sex <- sapply(eur_cluster_sex,function(x) glue('{x} ({signif(x/sum(eur_cluster_sex)*100,3)}%)'))
eur_cluster_self_reported_ancestry <- table(eur_cluster$RACE)
eur_cluster_self_reported_ancestry <- sapply(eur_cluster_self_reported_ancestry,function(x) glue('{x} ({signif(x/sum(eur_cluster_self_reported_ancestry)*100,3)}%)'))
eur_cluster_HBeAg <- table(eur_cluster$BASELINE_HBEAG_STATUS)
eur_cluster_HBeAg <- sapply(eur_cluster_HBeAg,function(x) glue('{x} ({signif(x/sum(eur_cluster_HBeAg)*100,3)}%)'))

asn_cluster <- dplyr::filter(tmp_cluster,gr_cluster==1 & RACE == 'ASIAN')
king_unrelated <- data.table::fread('~/G2G-HBV/data/processed_POP_asian_GT_A_C_D_B_TYPE_pPCA_LOCO_FALSE/kingunrelated.txt',header = F)
asn_cluster <- dplyr::filter(asn_cluster,host_id %in% king_unrelated$V1)
asn_cluster_age <- glue("{median(asn_cluster$AGE)} ({paste(range(asn_cluster$AGE),collapse = '-')})")
names(asn_cluster_age) <- "Median (range)"
asn_cluster_GT <- table(asn_cluster$GT)
merged_names <- union(names(eur_cluster_GT),names(asn_cluster_GT))
merged_names <- merged_names[!is.na(merged_names)]
asn_cluster_GT <- asn_cluster_GT[merged_names]
names(asn_cluster_GT) <- merged_names
eur_cluster_GT <- eur_cluster_GT[merged_names]
names(eur_cluster_GT) <- merged_names
asn_cluster_GT[is.na(asn_cluster_GT)] <- 0
eur_cluster_GT[is.na(eur_cluster_GT)] <- 0
eur_cluster_GT <- sapply(eur_cluster_GT,function(x) glue('{x} ({signif(x/sum(eur_cluster_GT)*100,3)}%)'))
asn_cluster_GT <- sapply(asn_cluster_GT,function(x) glue('{x} ({signif(x/sum(asn_cluster_GT)*100,3)}%)'))
asn_cluster_sex <- table(asn_cluster$SEX)[names(eur_cluster_sex)]
asn_cluster_sex <- sapply(asn_cluster_sex,function(x) glue('{x} ({signif(x/sum(asn_cluster_sex)*100,3)}%)'))
asn_cluster_self_reported_ancestry <- table(asn_cluster$RACE)[names(eur_cluster_self_reported_ancestry)]
asn_cluster_self_reported_ancestry[is.na(asn_cluster_self_reported_ancestry)] <- 0
asn_cluster_self_reported_ancestry <- sapply(asn_cluster_self_reported_ancestry,function(x) glue('{x} ({signif(x/sum(asn_cluster_self_reported_ancestry)*100,3)}%)'))
asn_cluster_HBeAg <- table(asn_cluster$BASELINE_HBEAG_STATUS)
asn_cluster_HBeAg <- sapply(asn_cluster_HBeAg,function(x) glue('{x} ({signif(x/sum(asn_cluster_HBeAg)*100,3)}%)'))

df_eur <- rbind(as.matrix(eur_cluster_age),as.matrix(eur_cluster_sex),as.matrix(eur_cluster_GT),as.matrix(eur_cluster_HBeAg))

df_asn <- rbind(as.matrix(asn_cluster_age),as.matrix(asn_cluster_sex),as.matrix(asn_cluster_GT),as.matrix(asn_cluster_HBeAg))
df_merged <- cbind(df_asn,df_eur)
colnames(df_merged) <- c(paste0('Asian(N=',nrow(asn_cluster),')'),paste0('European (N=',nrow(eur_cluster),')'))

kable(df_merged, format = "html", caption = "Population Characteristics of European and Asian Clusters", booktabs = T) %>%
  kable_styling(latex_options = "HOLD_position") %>%
  kableExtra::group_rows("Age", 1, 1) %>%
  kableExtra::group_rows("Sex", 2, 3) %>%
  kableExtra::group_rows("HBV Genotype", 4, 4 + length(eur_cluster_GT) - 1) %>%
  kableExtra::group_rows("HBeAg",4 + length(eur_cluster_GT) , 4 + length(eur_cluster_GT)  + 1)


dat_covars_raw_asn <- dplyr::filter(dat_covars_raw,host_id %in% asn_cluster$host_id)
dat_covars_raw_eur <- dplyr::filter(dat_covars_raw,host_id %in% eur_cluster$host_id)

table(dat_covars_raw_asn$COUNTRY)
table(dat_covars_raw_eur$COUNTRY)

## Table 2 - G2G SLC10A1 and preS1 association

#Get Stats from PLINK FIRTH output (p and OR)
BindPlinkPVal <- function(df,result_directory,AA_name){
  plink_result <- data.table::fread(glue("{result_directory}plink_aa_{AA_name}.PHENO1.glm.firth")) %>% dplyr::filter(TEST == 'ADD')
  plink_freq <- data.table::fread(glue("{result_directory}plink_aa_{AA_name}.frq.cc"))
  plink_df <- dplyr::left_join(plink_result,plink_freq,by=c('ID'='SNP')) %>% dplyr::select(ID,OR_Firth = OR,SE_Firth = SE,L95,U95,AF_Cases = MAF_A,
                                                                                           AF_Ctl = MAF_U,N_Cases = NCHROBS_A ,N_Ctl = NCHROBS_U ) %>%
    dplyr::mutate(N_Cases = N_Cases / 2,N_Ctl = N_Ctl / 2)
  return(dplyr::left_join(df,plink_df,by=c('SNPID' = 'ID')))
}
#Get Stats from SAIGE output (p and OR)
GetResults <- function(result_directory,hits,SNP){
  SAIGE_Results <- lapply(hits,function(x) tryCatch(data.table::fread(glue::glue("{result_directory}/SAIGE_SPA_TEST_LOCO_FALSE_{x}")) %>% dplyr::filter(SNPID == SNP) %>% dplyr::select(SNPID,p.value), error = function(e) return(data.frame(SNPID = SNP,p.value = NA))))
  names(SAIGE_Results) <- hits
  return(SAIGE_Results)
}
#Parse output into readable format
ParseSigHits <- function(hits_df){
  sig_hits <- do.call(rbind,hits_df) %>% dplyr::select(rsid=SNPID,OR_Firth,SE_Firth,p=p.value,N_Cases,N_Ctl,AF_Cases,AF_Ctl,L95,U95)
  sig_hits$AF_Cases <- signif(sig_hits$AF_Cases,3)
  sig_hits$AF_Ctl <- signif(sig_hits$AF_Ctl,3)
  
  a <- sig_hits$AF_Cases * sig_hits$N_Cases
  b <- sig_hits$AF_Ctl * sig_hits$N_Cases
  c <- sig_hits$N_Cases * (1 - sig_hits$AF_Cases)
  d <- sig_hits$N_Ctl * (1 - sig_hits$AF_Ctl)
  OR_Calc <- (a/c) / (b/d) #calculate OR based on Case vs CTL
  sig_hits$HBV_AA <- names(slc10a1_hits)
  sig_hits$OR_calc <- signif(OR_Calc,3)
  sig_hits$p <- as.character(signif(sig_hits$p,3))
  sig_hits <- sig_hits %>% dplyr::select(HBV_AA,`N (Alt AA)`=N_Cases,`N (Ref AA)`=N_Ctl,`SNP AF (Alt AA)`=AF_Cases,`SNP AF (Ref AA)`=AF_Ctl,p,OR_Firth,SE_Firth,L95,U95)%>% dplyr::arrange(as.numeric(p))
  return(sig_hits)
}
#Get All Hits
slc10a1_hits <- GetResults('~/G2G-HBV/data/results_POP_asian_GT_A_C_D_B_TYPE_pPCA_LOCO_FALSE/',hits = c('gene_S_pos_0017_A','gene_S_pos_0035_R','gene_S_pos_0051_P'),SNP = 'rs2296651')
slc10a1_hits_bind <- lapply(1:length(slc10a1_hits),function(i) BindPlinkPVal(slc10a1_hits[[i]],'~/G2G-HBV/data/results_POP_asian_GT_A_C_D_B_TYPE_pPCA_LOCO_FALSE/',names(slc10a1_hits)[i]))
#Get Hit in Genotype C subanalysis
slc10a1_hits_geno_c <- GetResults('~/G2G-HBV/data/results_POP_asian_GT_C_TYPE_pPCA_LOCO_FALSE/',hits = c('gene_S_pos_0017_A','gene_S_pos_0035_R','gene_S_pos_0051_P'),SNP = 'rs2296651')
slc10a1_hits_bind_geno_c <- lapply(1:length(slc10a1_hits_geno_c),function(i) BindPlinkPVal(slc10a1_hits_geno_c[[i]],'~/G2G-HBV/data/results_POP_asian_GT_C_TYPE_pPCA_LOCO_FALSE/',names(slc10a1_hits_geno_c)[i]))
#Get Hit in Genotype B subanalysis
slc10a1_hits_geno_b <- GetResults('~/G2G-HBV/data/results_POP_asian_GT_B_TYPE_pPCA_LOCO_FALSE/',hits = c('gene_S_pos_0017_A','gene_S_pos_0035_R','gene_S_pos_0051_P'),SNP = 'rs2296651')
slc10a1_hits_bind_geno_b <- lapply(1:length(slc10a1_hits_geno_b),function(i) BindPlinkPVal(slc10a1_hits_geno_b[[i]],'~/G2G-HBV/data/results_POP_asian_GT_B_TYPE_pPCA_LOCO_FALSE/',names(slc10a1_hits_geno_b)[i]))

#Parse data into readable format
sig_hits <- ParseSigHits(slc10a1_hits_bind)
sig_hits_c <- ParseSigHits(slc10a1_hits_bind_geno_c)
sig_hits_b <- ParseSigHits(slc10a1_hits_bind_geno_b)


merged_df <- rbind(sig_hits %>% dplyr::arrange(HBV_AA),sig_hits_c %>% dplyr::arrange(HBV_AA) ,sig_hits_b %>% dplyr::arrange(HBV_AA))
merged_df$HBV_Protein <- sapply(merged_df$HBV_AA,function(x) strsplit(x=strsplit(x=x,split = 'gene_')[[1]][2],split = '_pos')[[1]][1])
merged_df$HBV_Pos <- as.numeric(sapply(merged_df$HBV_AA,function(x) strsplit(strsplit(x,split = 'pos_')[[1]][2],split = '_')[[1]][1]))

dat_pathogen_raw_asn <- dplyr::filter(dat_pathogen_raw,ID %in% asn_cluster$host_id)
dat_pathogen_raw_asn <- dat_pathogen_raw_asn[,!apply(dat_pathogen_raw_asn,2,function(x) (sum(as.numeric(x)==0,na.rm = T) + sum(is.na(x))) == length(x))]
all_AA <- colnames(dat_pathogen_raw_asn[,apply(dat_pathogen_raw_asn,2,function(x) min(table(x)[!is.na(table(x))]) >= 1)])
all_AA_prot <- sapply(all_AA,function(x) strsplit(x=strsplit(x=x,split = 'gene_')[[1]][2],split = '_pos')[[1]][1])
all_AA_pos <- as.numeric(sapply(all_AA,function(x) strsplit(strsplit(x,split = 'pos_')[[1]][2],split = '_')[[1]][1]))
all_AA_residual <- sapply(all_AA,function(x) strsplit(strsplit(x,split = 'pos_')[[1]][2],split = '_')[[1]][2])
all_AA_df <- data.frame(Prot = all_AA_prot,Pos = all_AA_pos,Residual = all_AA_residual,stringsAsFactors = F)

dat_pathogen_raw_asn_GT_B <- dplyr::filter(dat_pathogen_raw,ID %in% asn_cluster$host_id[which(asn_cluster$GT == 'B')])
dat_pathogen_raw_asn_GT_B <- dat_pathogen_raw_asn_GT_B[,!apply(dat_pathogen_raw_asn_GT_B,2,function(x) (sum(as.numeric(x)==0,na.rm = T) + sum(is.na(x))) == length(x))]
all_AA_GT_B <- colnames(dat_pathogen_raw_asn_GT_B[,apply(dat_pathogen_raw_asn_GT_B,2,function(x) min(table(x)[!is.na(table(x))]) >= 1)])
all_AA_prot_GT_B <- sapply(all_AA_GT_B,function(x) strsplit(x=strsplit(x=x,split = 'gene_')[[1]][2],split = '_pos')[[1]][1])
all_AA_pos_GT_B <- as.numeric(sapply(all_AA_GT_B,function(x) strsplit(strsplit(x,split = 'pos_')[[1]][2],split = '_')[[1]][1]))
all_AA_residual_GT_B <- sapply(all_AA_GT_B,function(x) strsplit(strsplit(x,split = 'pos_')[[1]][2],split = '_')[[1]][2])

dat_pathogen_raw_asn_GT_C <- dplyr::filter(dat_pathogen_raw,ID %in% asn_cluster$host_id[which(asn_cluster$GT == 'C')])
dat_pathogen_raw_asn_GT_C <- dat_pathogen_raw_asn_GT_C[,!apply(dat_pathogen_raw_asn_GT_C,2,function(x) (sum(as.numeric(x)==0,na.rm = T) + sum(is.na(x))) == length(x))]
all_AA_GT_C <- colnames(dat_pathogen_raw_asn_GT_C[,apply(dat_pathogen_raw_asn_GT_C,2,function(x) min(table(x)[!is.na(table(x))]) >= 1)])
all_AA_prot_GT_C <- sapply(all_AA_GT_C,function(x) strsplit(x=strsplit(x=x,split = 'gene_')[[1]][2],split = '_pos')[[1]][1])
all_AA_pos_GT_C <- as.numeric(sapply(all_AA_GT_C,function(x) strsplit(strsplit(x,split = 'pos_')[[1]][2],split = '_')[[1]][1]))
all_AA_residual_GT_C <- sapply(all_AA_GT_C,function(x) strsplit(strsplit(x,split = 'pos_')[[1]][2],split = '_')[[1]][2])

all_AA_df <- data.frame(Prot = all_AA_prot,Pos = all_AA_pos,Residual = all_AA_residual,stringsAsFactors = F)
all_AA_df_GT_B <- data.frame(Prot = all_AA_prot_GT_B,Pos = all_AA_pos_GT_B,Residual = all_AA_residual_GT_B,stringsAsFactors = F)
all_AA_df_GT_C <- data.frame(Prot = all_AA_prot_GT_C,Pos = all_AA_pos_GT_C,Residual = all_AA_residual_GT_C,stringsAsFactors = F)

merged_df$HBV_VarAA <- rep(NA,nrow(merged_df))
merged_df$HBV_VarAA[1:nrow(sig_hits)] <- mapply(FUN = function(x,y) unique(dplyr::filter(all_AA_df,Prot==x & Pos == y)$Residual),merged_df$HBV_Protein[1:nrow(sig_hits)],merged_df$HBV_Pos[1:nrow(sig_hits)])

merged_df$HBV_VarAA[(nrow(sig_hits) + 1):(nrow(sig_hits) + nrow(sig_hits_c))] <- mapply(FUN = function(x,y) unique(dplyr::filter(all_AA_df_GT_C,Prot==x & Pos == y)$Residual),merged_df$HBV_Protein[(nrow(sig_hits) + 1):(nrow(sig_hits) + nrow(sig_hits_c))],merged_df$HBV_Pos[(nrow(sig_hits) + 1):(nrow(sig_hits) + nrow(sig_hits_c))])

merged_df$HBV_VarAA[(nrow(sig_hits) + nrow(sig_hits_c) + 1):(nrow(sig_hits) + nrow(sig_hits_c) + nrow(sig_hits_b))] <- mapply(FUN = function(x,y) unique(dplyr::filter(all_AA_df_GT_B,Prot==x & Pos == y)$Residual),merged_df$HBV_Protein[(nrow(sig_hits) + nrow(sig_hits_c) + 1):(nrow(sig_hits) + nrow(sig_hits_c) + nrow(sig_hits_b))],merged_df$HBV_Pos[(nrow(sig_hits) + nrow(sig_hits_c) + 1):(nrow(sig_hits) + nrow(sig_hits_c) + nrow(sig_hits_b))])

merged_df$AssocAA <- sapply(merged_df$HBV_AA,function(x) strsplit(strsplit(x,split = 'pos_')[[1]][2],split = '_')[[1]][2])

merged_df$HBV_VarAA <- sapply(1:nrow(merged_df),function(i) paste0(setdiff(merged_df$HBV_VarAA[[i]],merged_df$AssocAA[i]),collapse = ', '))


merged_df$AssocAA <- sapply(1:nrow(merged_df),function(i) paste0(merged_df$AssocAA[i],' (',round(merged_df$`N (Alt AA)`[i] / (merged_df$`N (Alt AA)`[i] + merged_df$`N (Ref AA)`[i]) * 100,1),'%)'))
merged_df$EAF_Ratio <- paste0(signif(merged_df$`SNP AF (Ref AA)`,2),' : ',signif(merged_df$`SNP AF (Alt AA)`,2))

merged_df$CI <- paste0(round(merged_df$L95,1),'-',round(merged_df$U95,1))

kable(merged_df %>% dplyr::select('PreS1 Position' = HBV_Pos,'Associated AA (%)' = AssocAA,'Other AA' = HBV_VarAA, 'rs2296651 MAF (Other AA : Associated AA)' = EAF_Ratio, p, 'OR' = OR_Firth,`95% CI` = CI) %>% dplyr::mutate(OR = round(OR,1)), format = "html", caption = "", booktabs = T) %>%
  kable_styling(latex_options = "HOLD_position") %>%
  kableExtra::group_rows("All HBV Genotypes", 1, nrow(sig_hits)) %>%
  kableExtra::group_rows("HBV Genotype C", nrow(sig_hits) + 1, nrow(sig_hits) + nrow(sig_hits_c)) %>%
  kableExtra::group_rows("HBV Genotype B", nrow(sig_hits) + nrow(sig_hits_c) + 1, nrow(sig_hits) + nrow(sig_hits_c) + nrow(sig_hits_b))


## Figure 1 - Circos Plot
source('./scripts/Figures/circos_plot_overlap.R')

## Figure 2 - NTCP Haplotypes
source('./scripts/Haplotypes/GetMyrHaplotypes.R')

GetAvgIntraHostProp <- function(prop_by_sample,prop_by_sample_WT){
  uniq_haplotypes <- unique(c(prop_by_sample$AA_pos_of_interest,prop_by_sample_WT$AA_pos_of_interest))
  avg_intra_host_prop <- lapply(uniq_haplotypes,function(x) {
    df <- dplyr::filter(prop_by_sample,AA_pos_of_interest==x & prop_reads > 0.15)
    if(nrow(df) == 0){
      return(data.frame(Avg_Prop_Pure = 0,Avg_Prop_Mixed = 0))
    }else{
      samples_present = df$ID
      mixed_samples <- samples_present[sapply(samples_present,function(k) nrow(dplyr::filter(prop_by_sample,ID == k & AA_pos_of_interest!=x & prop_reads > 0.15))>0)]
      pure_samples <- setdiff(samples_present,mixed_samples)
      return(data.frame(Avg_Prop_Pure = length(pure_samples)/length(unique(prop_by_sample$ID)),Avg_Prop_Mixed = length(mixed_samples)/length(unique(prop_by_sample$ID))))
    }
  })
  avg_intra_host_prop <- do.call(rbind,avg_intra_host_prop)
  
  avg_intra_host_prop_WT <- lapply(uniq_haplotypes,function(x) {
    df <- dplyr::filter(prop_by_sample_WT,AA_pos_of_interest==x & prop_reads > 0.15)
    if(nrow(df) == 0){
      return(data.frame(Avg_Prop_Pure = 0,Avg_Prop_Mixed = 0))
    }else{
      samples_present = df$ID
      mixed_samples <- samples_present[sapply(samples_present,function(k) nrow(dplyr::filter(prop_by_sample_WT,ID == k & AA_pos_of_interest!=x & prop_reads > 0.15))>0)]
      pure_samples <- setdiff(samples_present,mixed_samples)
      return(data.frame(Avg_Prop_Pure = length(pure_samples)/length(unique(prop_by_sample_WT$ID)),Avg_Prop_Mixed = length(mixed_samples)/length(unique(prop_by_sample_WT$ID))))
    }
  })
  avg_intra_host_prop_WT <- do.call(rbind,avg_intra_host_prop_WT)
  
  return(rbind(data.frame(Haplotype = uniq_haplotypes,avg_intra_host_prop,rs2296651 = 'GA',stringsAsFactors = F),
               data.frame(Haplotype = uniq_haplotypes,avg_intra_host_prop_WT,rs2296651 = 'GG',stringsAsFactors = F)))
}

avg_freq_geno_C_merged <- GetAvgIntraHostProp(prop_geno_C_by_sample,prop_geno_C_WT_by_sample) %>% dplyr::mutate(Avg_Prop = Avg_Prop_Pure + Avg_Prop_Mixed)

avg_freq_geno_C_merged$Haplotype <- factor(avg_freq_geno_C_merged$Haplotype,levels = avg_freq_geno_C_merged$Haplotype[order(avg_freq_geno_C_merged$Avg_Prop,decreasing = T)],
                                           labels = avg_freq_geno_C_merged$Haplotype[order(avg_freq_geno_C_merged$Avg_Prop,decreasing = T)])
rare_haplotypes_C <- avg_freq_geno_C_merged %>% dplyr::group_by(Haplotype) %>% dplyr::arrange(desc(Avg_Prop)) %>% dplyr::filter(row_number() == 1) %>% dplyr::filter(Avg_Prop < 0.01)
avg_freq_geno_C_merged_single <- data.frame(avg_freq_geno_C_merged %>% dplyr::select(Haplotype,Avg_Prop=Avg_Prop_Pure,rs2296651),type = 'Single') %>% dplyr::filter(!Haplotype %in% rare_haplotypes_C$Haplotype)
avg_freq_geno_C_merged_multiple <- data.frame(avg_freq_geno_C_merged %>% dplyr::select(Haplotype,Avg_Prop=Avg_Prop_Mixed,rs2296651),type = 'Mixed') %>% dplyr::filter(!Haplotype %in% rare_haplotypes_C$Haplotype)
avg_freq_geno_C_merged_multiple$Avg_Prop <- avg_freq_geno_C_merged_single$Avg_Prop + avg_freq_geno_C_merged_multiple$Avg_Prop

p1_geno_C <- ggplot2::ggplot() + 
  geom_bar(data= avg_freq_geno_C_merged_single,stat = 'identity',position='dodge',colour = "black",width = 0.7,aes(x=Haplotype,y=Avg_Prop,fill = rs2296651,alpha = type)) +
  geom_bar(data= avg_freq_geno_C_merged_multiple,stat = 'identity',position='dodge',width = 0.7,colour = "black",aes(x=Haplotype,y=Avg_Prop,fill = rs2296651,alpha = type)) +
  xlab("PreS1 Haplotype (Positions:17,35,51)") +
  ylab("Fraction of Samples\n(Intra-Host Proportion > 15%)") + ylim(0,1) +
  theme(plot.margin=unit(c(1,0.5,0.5,0.5),"cm")) + scale_alpha_discrete('Intra-Host Composition',range = c(0.35,1)) + ggtitle('Genotype C')

geno_C_haplotype_freq <- rep(prop_geno_C$results_shotgun_bind$AA_haplotype,round(10 * prop_geno_C$results_shotgun_bind$prop_reads))
p2_geno_C <- ggplot() + geom_logo(geno_C_haplotype_freq,method = 'prob',seq_type='aa') + 
  theme_logo() + annotate('text', x=17, y=1, label='*',size = 10)  + 
  annotate('text', x=35, y=1, label='*',size = 10) + theme(legend.position = "none")+
  annotate('text', x=51, y=1, label='*',size = 10) + scale_x_continuous(name ="PreS1 Position", 
                                                                        seq(1,59,2)) + ggtitle('Genotype C: rs2296651-GA (NTCP S267F)') + scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1,1.15),labels = c(0,0.25,0.5,0.75,1,''),limits = c(0,1.15))

geno_C_WT_haplotype_freq <- rep(prop_WT_geno_C$results_shotgun_bind$AA_haplotype,round(10 * prop_WT_geno_C$results_shotgun_bind$prop_reads))
p3_geno_C <- ggplot() + geom_logo(geno_C_WT_haplotype_freq,method = 'prob',seq_type='aa') + 
  theme_logo() + annotate('text', x=17, y=1, label='*',size = 10)  + 
  annotate('text', x=35, y=1, label='*',size = 10) + theme(legend.position = "none") +
  annotate('text', x=51, y=1, label='*',size = 10) + scale_x_continuous(name ="PreS1 Position", 
                                                                        seq(1,59,2)) + ggtitle('Genotype C: rs2296651-GG') + scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1,1.15),labels = c(0,0.25,0.5,0.75,1,''),limits = c(0,1.15))
geno_C_AA_Plot <- ggpubr::ggarrange(ggpubr::ggarrange(p2_geno_C,p3_geno_C,nrow = 2,common.legend = F,labels = c('A)','B)')),p1_geno_C,ncol = 2,widths = c(1.3,0.9),labels = c('','E)'))

avg_freq_geno_B_merged <- GetAvgIntraHostProp(prop_geno_B_by_sample,prop_geno_B_WT_by_sample) %>% dplyr::mutate(Avg_Prop = Avg_Prop_Pure + Avg_Prop_Mixed)

avg_freq_geno_B_merged$Haplotype <- factor(avg_freq_geno_B_merged$Haplotype,levels = avg_freq_geno_B_merged$Haplotype[order(avg_freq_geno_B_merged$Avg_Prop,decreasing = T)],
                                           labels = avg_freq_geno_B_merged$Haplotype[order(avg_freq_geno_B_merged$Avg_Prop,decreasing = T)])
rare_haplotypes_C <- avg_freq_geno_B_merged %>% dplyr::group_by(Haplotype) %>% dplyr::arrange(desc(Avg_Prop)) %>% dplyr::filter(row_number() == 1) %>% dplyr::filter(Avg_Prop < 0.01)
avg_freq_geno_B_merged_single <- data.frame(avg_freq_geno_B_merged %>% dplyr::select(Haplotype,Avg_Prop=Avg_Prop_Pure,rs2296651),type = 'Single') %>% dplyr::filter(!Haplotype %in% rare_haplotypes_C$Haplotype)
avg_freq_geno_B_merged_multiple <- data.frame(avg_freq_geno_B_merged %>% dplyr::select(Haplotype,Avg_Prop=Avg_Prop_Mixed,rs2296651),type = 'Mixed') %>% dplyr::filter(!Haplotype %in% rare_haplotypes_C$Haplotype)
avg_freq_geno_B_merged_multiple$Avg_Prop <- avg_freq_geno_B_merged_single$Avg_Prop + avg_freq_geno_B_merged_multiple$Avg_Prop

p1_geno_B <- ggplot2::ggplot() + 
  geom_bar(data= avg_freq_geno_B_merged_single,stat = 'identity',position='dodge',colour = "black",width = 0.7,aes(x=Haplotype,y=Avg_Prop,fill = rs2296651,alpha = type)) +
  geom_bar(data= avg_freq_geno_B_merged_multiple,stat = 'identity',position='dodge',colour = "black",width = 0.7,aes(x=Haplotype,y=Avg_Prop,fill = rs2296651,alpha = type)) +
  xlab("PreS1 Haplotype (Positions:17,35,51)") +
  ylab("Fraction of Samples\n(Intra-Host Proportion > 15%)") + ylim(0,1) +
  theme(plot.margin=unit(c(1,0.5,0.5,0.5),"cm")) + scale_alpha_discrete('Intra-Host Composition',range = c(0.35,1))+ ggtitle('Genotype B')

geno_B_haplotype_freq <- rep(prop_geno_B$results_shotgun_bind$AA_haplotype,round(10 * prop_geno_B$results_shotgun_bind$prop_reads))
p2_geno_B <- ggplot() + geom_logo(geno_B_haplotype_freq,method = 'prob',seq_type='aa') + 
  theme_logo() + annotate('text', x=17, y=1, label='*',size = 10)  + 
  annotate('text', x=35, y=1, label='*',size = 10) + theme(legend.position = "none")+
  annotate('text', x=51, y=1, label='*',size = 10) + scale_x_continuous(name ="PreS1 Position", 
                                                                        seq(1,59,2)) + ggtitle('Genotype B: rs2296651-GA (NTCP S267F)') +  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1,1.15),labels = c(0,0.25,0.5,0.75,1,''),limits = c(0,1.15))
geno_B_WT_haplotype_freq <- rep(prop_WT_geno_B$results_shotgun_bind$AA_haplotype,round(10 * prop_WT_geno_B$results_shotgun_bind$prop_reads))
p3_geno_B <- ggplot() + geom_logo(geno_B_WT_haplotype_freq,method = 'prob',seq_type='aa') + 
  theme_logo() + annotate('text', x=17, y=1, label='*',size = 10) + 
  annotate('text', x=35, y=1, label='*',size = 10) + theme(legend.position = "none") +
  annotate('text', x=51, y=1, label='*',size = 10) + scale_x_continuous(name ="PreS1 Position", 
                                                                        seq(1,59,2)) + ggtitle('Genotype B: rs2296651-GG') + scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1,1.15),labels = c(0,0.25,0.5,0.75,1,''),limits = c(0,1.15))

geno_B_AA_Plot <- ggpubr::ggarrange(ggpubr::ggarrange(p2_geno_B,p3_geno_B,nrow = 2,common.legend = F,labels = c('C)','D)')),p1_geno_B,ncol = 2,widths = c(1.3,0.9),labels = c('','F)'))

ggpubr::ggarrange(geno_C_AA_Plot,geno_B_AA_Plot,nrow = 2)

##Suppl Figure 7 - NTCP Heatmap
PlotHeatMap <- function(prop_geno_by_sample,title,row.names=F){
  by_sample_mat <- matrix(0,nrow = length(unique(prop_geno_by_sample$ID)),ncol = length(unique(prop_geno_by_sample$AA_pos_of_interest)))
  rownames(by_sample_mat) <- unique(prop_geno_by_sample$ID)
  colnames(by_sample_mat) <- unique(prop_geno_by_sample$AA_pos_of_interest)
  
  for(i in 1:nrow(by_sample_mat)){
    cur_df <- dplyr::filter(prop_geno_by_sample,ID == unique(prop_geno_by_sample$ID)[i])
    by_sample_mat[i,cur_df$AA_pos_of_interest] <- cur_df$prop_reads
  }
  # return(pheatmap(by_sample_mat,color = colorRampPalette(brewer.pal(n = 4, name ="BuGn"))(10),show_rownames = F,main = title)[[4]])
  
  return(pheatmap(by_sample_mat,color = c('grey',colorRampPalette(brewer.pal(n = 4, name ="BuGn"))(100)),breaks = c(0,0.01,seq(0.011,1,0.01)),show_rownames = row.names,main = title)[[4]])
}
geno_C_heat <- PlotHeatMap(prop_geno_C_by_sample,title = 'Genotype C (rs2296651-GA)')
geno_B_heat <- PlotHeatMap(prop_geno_B_by_sample,title = 'Genotype B (rs2296651-GA)')
geno_C_heat_WT <- PlotHeatMap(prop_geno_C_WT_by_sample %>% dplyr::filter(!grepl('X',AA_pos_of_interest)) ,title = 'Genotype C (rs2296651-GG)')
geno_B_heat_WT <- PlotHeatMap(prop_geno_B_WT_by_sample %>% dplyr::filter(!grepl('X',AA_pos_of_interest)),title = 'Genotype B (rs2296651-GG)')

p1 = grid.arrange(arrangeGrob(grobs= list(geno_C_heat,geno_C_heat_WT),ncol=2))
p2 = grid.arrange(arrangeGrob(grobs= list(geno_B_heat,geno_B_heat_WT),ncol=2))
ggarrange(p1,p2,nrow = 2,labels = c('A)','B)'))


## Table 3 - NTCP S267F Escape vs NTCP S267F NonEscape vs NTCP WT
#Get Nucleotide variants associated with each biomarker
nt_tbl <- data.table::fread('~/G2G-HBV/data/raw/gilead_20181126/viral_seq/OUT_nt_binary_table.txt') %>% dplyr::select(ID=V1,pos_1817_T,pos_1838_G,pos_1896_A,pos_1762_A,pos_1764_G)

outcomes <- dat_covars_raw %>% dplyr::select(host_id,HBeAg = BASELINE_HBEAG_STATUS,ALT = 'BASELINE_ALT_U/L',VL_Raw = 'BASELINE_HBVDNA_IU/mL',VL_Dil = 'BASELINE_HBVDNA_Dil_IU/mL',HBsAg = 'BASELINE_HBSAG_log10_IU/mL',GT) %>% dplyr::mutate(VL = ifelse(grepl(pattern = '>',x = VL_Raw) | is.na(VL_Raw),VL_Dil,VL_Raw)) %>% dplyr::select(-VL_Raw,-VL_Dil) %>% dplyr::mutate(HBeAg = factor(HBeAg),ALT = log10(ALT),VL = log10(as.numeric(VL)))

#Get NTCP Haplotype for each sample, join with covariates
WT_Prop <- rbind(prop_geno_C_SGH_SGQ,prop_geno_B_SKN) %>% dplyr::left_join(dat_ids,by=c('ID' = 'pathogen_id')) %>% 
  dplyr::left_join(dat_pathogen_raw %>% dplyr::select(ID,gene_PC_C_pos_0009_V),by=c('host_id'='ID')) %>%
  dplyr::left_join(nt_tbl) %>%
  dplyr::left_join(dat_covars_discrete) %>% dplyr::left_join(dat_covars_numeric) %>%
  dplyr::left_join(outcomes) %>% dplyr::filter(!is.na(VL) & !is.na(HBsAg) & !is.na(HBeAg) & !is.na(pos_1817_T) & !is.na(pos_1838_G) & !is.na(pos_1896_A)) %>%
  dplyr::left_join(dat_covars_raw %>% dplyr::select(host_id,BMI=BASELINE_BMI),by=c('host_id'='host_id'))
#Group according to rs2296651 and presence/absences of escape haplotypes (15%)
treatment_df <- WT_Prop %>% dplyr::mutate(treat = ifelse(Freq < 0.15 & rs2296651 == 'GG',0,ifelse(Freq < 0.15 & rs2296651 == 'GA',1,2)))

ctl_group <- dplyr::filter(treatment_df,treat == 0) #NTCP WT 
non_escape_group <- dplyr::filter(treatment_df,treat != 2) #NTCP WT vs NTCP S267F Non-Escape
escape_group <- dplyr::filter(treatment_df,treat != 0) #NTCP S267F Non-Escape vs NTCP S267F Escape
escape_group$treat <- escape_group$treat - 1

escape_vs_ctl_group <- dplyr::filter(treatment_df,treat == 0 | treat == 2)  #NTCP S267F Non-Escape vs NTCP WT
escape_vs_ctl_group$treat <- ifelse(escape_vs_ctl_group$treat==2,1,0)
#Fit LM for each group
fit1_reg.alt <- lm(ALT  ~ treat + AGE + SEX + HBeAg + pos_1838_G + host_PC1 + host_PC2 + host_PC3 + host_PC4 + pPC1 + pPC2 + pPC3 + pPC4 + pPC5 + pPC6, data = non_escape_group)
fit1_reg.vl <- lm(VL  ~ treat + AGE + SEX + HBeAg  + pos_1838_G + pos_1764_G + host_PC1 + host_PC2 + host_PC3 + host_PC4 + pPC1 + pPC2 + pPC3 + pPC4 + pPC5 + pPC6, data = non_escape_group)
fit1_reg.sa <- lm(HBsAg  ~ treat + AGE + SEX + HBeAg + pos_1764_G  + host_PC1 + host_PC2 + host_PC3 + host_PC4 + pPC1 + pPC2 + pPC3 + pPC4 + pPC5 + pPC6, data = non_escape_group)


fit2_reg.alt <- lm(ALT ~ treat + AGE + SEX + HBeAg + pos_1838_G  + host_PC1 + host_PC2 + host_PC3 + host_PC4 + pPC1 + pPC2 + pPC3 + pPC4 + pPC5 + pPC6, data = escape_group)
fit2_reg.vl <- lm(VL  ~ treat + AGE + SEX + HBeAg  + pos_1838_G + pos_1764_G + host_PC1 + host_PC2 + host_PC3 + host_PC4 + pPC1 + pPC2 + pPC3 + pPC4 + pPC5 + pPC6, data = escape_group)
fit2_reg.sa <- lm(HBsAg  ~ treat + AGE + SEX + HBeAg + pos_1764_G  + host_PC1 + host_PC2 + host_PC3 + host_PC4 + pPC1 + pPC2 + pPC3 + pPC4 + pPC5 + pPC6, data = escape_group)


fit3_reg.alt <- lm(ALT ~ treat + AGE + SEX + HBeAg  + pos_1838_G  + host_PC1 + host_PC2 + host_PC3 + host_PC4 + pPC1 + pPC2 + pPC3 + pPC4 + pPC5 + pPC6, data = escape_vs_ctl_group)
fit3_reg.vl <- lm(VL ~ treat + AGE + SEX + HBeAg  + pos_1838_G + pos_1764_G  + host_PC1 + host_PC2 + host_PC3 + host_PC4 + pPC1 + pPC2 + pPC3 + pPC4 + pPC5 + pPC6, data = escape_vs_ctl_group)
fit3_reg.sa <- lm(HBsAg ~ treat + AGE + SEX + HBeAg  + pos_1764_G  + host_PC1 + host_PC2 + host_PC3 + host_PC4 + pPC1 + pPC2 + pPC3 + pPC4 + pPC5 + pPC6, data = escape_vs_ctl_group)

#Get Mean/ Freq

tbl_groupA <- data.frame(VL = mean(ctl_group$VL),VL_SD = sd(ctl_group$VL),
                         ALT = mean(ctl_group$ALT),ALT_SD = sd(ctl_group$ALT),
                         HBsAg = mean(ctl_group$HBsAg),HBsAg_SD = sd(ctl_group$HBsAg))

groupB <- dplyr::filter(non_escape_group,treat == 1)
tbl_groupB <- data.frame(VL = mean(groupB$VL),VL_SD = sd(groupB$VL),
                         ALT = mean(groupB$ALT),ALT_SD = sd(groupB$ALT),
                         HBsAg = mean(groupB$HBsAg),HBsAg_SD = sd(groupB$HBsAg))
#Get p-values and Beta
p_groupB <- data.frame(VL.Beta = fit1_reg.vl$coefficients['treat'],
                       VL.SE = summary(fit1_reg.vl)$coefficients['treat','Std. Error'],
                       VL.p = summary(fit1_reg.vl)$coefficients['treat','Pr(>|t|)'],
                       ALT.Beta = fit1_reg.alt$coefficients['treat'],
                       ALT.SE = summary(fit1_reg.alt)$coefficients['treat','Std. Error'],
                       ALT.p = summary(fit1_reg.alt)$coefficients['treat','Pr(>|t|)'],
                       SA.Beta = fit1_reg.sa$coefficients['treat'],
                       SA.SE = summary(fit1_reg.sa)$coefficients['treat','Std. Error'],
                       SA.p = summary(fit1_reg.sa)$coefficients['treat','Pr(>|t|)'])

groupC <- dplyr::filter(escape_group,treat == 1)
tbl_groupC <- data.frame(VL = mean(groupC$VL),VL_SD = sd(groupC$VL),
                         ALT = mean(groupC$ALT),ALT_SD = sd(groupC$ALT),
                         HBsAg = mean(groupC$HBsAg),HBsAg_SD = sd(groupC$HBsAg))

p_groupC_B <- data.frame(VL.Beta = fit2_reg.vl$coefficients['treat'],
                         VL.SE = summary(fit2_reg.vl)$coefficients['treat','Std. Error'],
                         VL.p = summary(fit2_reg.vl)$coefficients['treat','Pr(>|t|)'],
                         ALT.Beta = fit2_reg.alt$coefficients['treat'],
                         ALT.SE = summary(fit2_reg.alt)$coefficients['treat','Std. Error'],
                         ALT.p = summary(fit2_reg.alt)$coefficients['treat','Pr(>|t|)'],
                         SA.Beta = fit2_reg.sa$coefficients['treat'],
                         SA.SE = summary(fit2_reg.sa)$coefficients['treat','Std. Error'],
                         SA.p = summary(fit2_reg.sa)$coefficients['treat','Pr(>|t|)'])

p_groupC_A <- data.frame(VL.Beta = fit3_reg.vl$coefficients['treat'],
                         VL.SE = summary(fit3_reg.vl)$coefficients['treat','Std. Error'],
                         VL.p = summary(fit3_reg.vl)$coefficients['treat','Pr(>|t|)'],
                         ALT.Beta = fit3_reg.alt$coefficients['treat'],
                         ALT.SE = summary(fit3_reg.alt)$coefficients['treat','Std. Error'],
                         ALT.p = summary(fit3_reg.alt)$coefficients['treat','Pr(>|t|)'],
                         SA.Beta = fit3_reg.sa$coefficients['treat'],
                         SA.SE = summary(fit3_reg.sa)$coefficients['treat','Std. Error'],
                         SA.p = summary(fit3_reg.sa)$coefficients['treat','Pr(>|t|)'])

## Figure 2 and Suppl Tbl 1-4 HLA G2G
source('~/G2G-HBV/src/plot-HLA_AssocPlink.R')
hla_result <- readRDS('~/G2G-HBV/PyHLA_Out/HLA_result.rds')
PC_C_pos_0160_hla_results <- hla_result$PC_C_pos_0160_hla_results
Pol_pos_0049_hla_results <- hla_result$Pol_pos_0049_hla_results
gene_PC_C_pos_0067_Y_eur_hla_results <- hla_result$gene_PC_C_pos_0067_Y_eur_hla_results
PC_C_pos_0160_single_association = data.table::fread('~/G2G-HBV/PyHLA_Out/Asian/hla_allele_gene_PC_C_pos_0160_A_assocAA.PHENO1.glm.logistic') %>% dplyr::filter(TEST == 'ADD')
PC_C_pos_0160_single_association$AA_POS <- sapply(PC_C_pos_0160_single_association$ID,function(x) strsplit(x=x,split = '_')[[1]][2])
PC_C_pos_0160_single_association$AA_Gene <- sapply(PC_C_pos_0160_single_association$ID,function(x) strsplit(x=x,split = '_')[[1]][1])
PC_C_pos_0160_single_association$AA_Residue <- sapply(PC_C_pos_0160_single_association$ID,function(x) strsplit(x=x,split = '_')[[1]][3])
PC_C_pos_0160_top_hit_results <- dplyr::filter(PC_C_pos_0160_single_association,AA_POS == 99 & AA_Gene == "A")
Pol_pos_0049_single_association = data.table::fread('~/G2G-HBV/PyHLA_Out/Asian/hla_allele_gene_Pol_pos_0049_N_assocAA_cond_snp_A_35_Y.PHENO1.glm.logistic') %>% dplyr::filter(TEST == 'ADD')
Pol_pos_0049_single_association$AA_POS <- sapply(Pol_pos_0049_single_association$ID,function(x) strsplit(x=x,split = '_')[[1]][2])
Pol_pos_0049_single_association$AA_Gene <- sapply(Pol_pos_0049_single_association$ID,function(x) strsplit(x=x,split = '_')[[1]][1])
Pol_pos_0049_single_association$AA_Residue <- sapply(Pol_pos_0049_single_association$ID,function(x) strsplit(x=x,split = '_')[[1]][3])
Pol_pos_0049_top_hit_results <- dplyr::filter(Pol_pos_0049_single_association,AA_POS == 131 & AA_Gene == "A")
Pol_pos_0067_single_association = data.table::fread('~/G2G-HBV/PyHLA_Out/European/hla_allele_gene_PC_C_pos_0067_Y_assocAA.PHENO1.glm.logistic') %>% dplyr::filter(TEST == 'ADD')
Pol_pos_0067_single_association$AA_POS <- sapply(Pol_pos_0067_single_association$ID,function(x) strsplit(x=x,split = '_')[[1]][2])
Pol_pos_0067_single_association$AA_Gene <- sapply(Pol_pos_0067_single_association$ID,function(x) strsplit(x=x,split = '_')[[1]][1])
Pol_pos_0067_single_association$AA_Residue <- sapply(Pol_pos_0067_single_association$ID,function(x) strsplit(x=x,split = '_')[[1]][3])
Pol_pos_0067_top_hit_results <- dplyr::filter(Pol_pos_0067_single_association,AA_POS == 153 & AA_Gene == "A")
system(glue::glue("../bin/plink2_linux --freq --bfile /home/mxu/G2G-HBV/PyHLA_Out/Asian/hla_allele_assocAA --out /home/mxu/G2G-HBV/PyHLA_Out/Asian/hla_allele_assocAA"))
AA_HLA_Freq <- data.table::fread('~/G2G-HBV/PyHLA_Out/Asian/hla_allele_assocAA.afreq')
system(glue::glue("../bin/plink2_linux --freq --bfile /home/mxu/G2G-HBV/PyHLA_Out/European/hla_allele_assocAA --out /home/mxu/G2G-HBV/PyHLA_Out/European/hla_allele_assocAA"))
AA_HLA_Freq_EUR <- data.table::fread('~/G2G-HBV/PyHLA_Out/European/hla_allele_assocAA.afreq')
merged_AA_HLA_results <- rbind(PC_C_pos_0160_top_hit_results %>% dplyr::select(-L95,-U95),Pol_pos_0049_top_hit_results %>% dplyr::select(-L95,-U95)) %>% dplyr::left_join(AA_HLA_Freq %>% dplyr::select(-OBS_CT),by = c('ID' = 'ID')) %>% dplyr::select(`Host Protein` = AA_Gene,`Protein Position` = AA_POS, `Test AA` = AA_Residue,OR,SE,Freq = ALT_FREQS,OBS_CT)
merged_AA_HLA_results <- rbind(merged_AA_HLA_results,Pol_pos_0067_top_hit_results %>% dplyr::select(-L95,-U95)%>% dplyr::left_join(AA_HLA_Freq %>% dplyr::select(-OBS_CT),by = c('ID' = 'ID')) %>% dplyr::select(`Host Protein` = AA_Gene,`Protein Position` = AA_POS, `Test AA` = AA_Residue,OR,SE,Freq = ALT_FREQS,OBS_CT))
merged_AA_HLA_results$OR <- paste0(round(merged_AA_HLA_results$OR,2),' (',round(exp(log(merged_AA_HLA_results$OR) - qt(p = 0.975,df=merged_AA_HLA_results$Freq*merged_AA_HLA_results$OBS_CT)*merged_AA_HLA_results$SE),2),',',round(exp(log(merged_AA_HLA_results$OR) + qt(p = 0.975,df=merged_AA_HLA_results$Freq*merged_AA_HLA_results$OBS_CT)*merged_AA_HLA_results$SE),2),')')
merged_AA_HLA_results <- merged_AA_HLA_results %>% dplyr::select(-SE,-OBS_CT)
merged_AA_HLA_results$`Host Protein` <- paste0('HLA-',merged_AA_HLA_results$`Host Protein`)
merged_AA_HLA_results$Freq <- round(merged_AA_HLA_results$Freq,2)
kable(merged_AA_HLA_results, format = "html", caption = "", booktabs = T) %>%
  kable_styling(latex_options = "HOLD_position") %>%
  kableExtra::group_rows(paste0("Asian - HBV Core Protein (Pos: 160, Residue: A), p = ",min(PC_C_pos_0160_hla_results$aa_results$uncond$P,na.rm = T)), 1, nrow(PC_C_pos_0160_top_hit_results)) %>%
  kableExtra::group_rows(paste0("Asian - HBV Polymerase (Pos: 49, Residue: N ), p = ",min(Pol_pos_0049_hla_results$aa_results$cond$P,na.rm = T)), nrow(PC_C_pos_0160_top_hit_results) + 1, nrow(PC_C_pos_0160_top_hit_results) + nrow(Pol_pos_0049_top_hit_results)) %>%
  kableExtra::group_rows(paste0("European - HBV Polymerase (Pos: 67, Residue: Y ), p = ",min(gene_PC_C_pos_0067_Y_eur_hla_results$aa_results$uncond$P,na.rm = T)), nrow(PC_C_pos_0160_top_hit_results) + nrow(Pol_pos_0049_top_hit_results) + 1, nrow(PC_C_pos_0160_top_hit_results) + nrow(Pol_pos_0049_top_hit_results) + nrow(Pol_pos_0067_top_hit_results))

hla_result <- readRDS('~/G2G-HBV/PyHLA_Out/HLA_result.rds')

hla_allele_results <- lapply(hla_result,function(x) x$allele_results$uncond %>% dplyr::filter(`#CHROM` == 'A') %>% 
                               dplyr::arrange(P) %>% dplyr::select(ID,P,OR,SE) %>% dplyr::mutate(P = scientific(P,digits = 3),OR = signif(OR,digits = 3),SE = signif(SE,digits = 3))) 

hla_allele_results$PC_C_pos_0160_A_hla_results <- hla_allele_results$PC_C_pos_0160_A_hla_results %>% dplyr::left_join(data.frame(ID = names(PC_C_160_allele_scores$freq$global_HLA_freq),Freq = signif(PC_C_160_allele_scores$freq$global_HLA_freq,3))) %>% dplyr::relocate(Freq)

hla_allele_results$Pol_pos_0049_N_hla_results <- hla_allele_results$Pol_pos_0049_N_hla_results

hla_allele_results_eur <- hla_allele_results$PC_C_pos_0067_Y_eur_hla_results %>% dplyr::left_join(data.frame(ID = names(PC_C_67_allele_scores$freq$global_HLA_freq),Freq = signif(PC_C_67_allele_scores$freq$global_HLA_freq,3))) %>% dplyr::relocate(Freq)


hla_allele_results_asn <- dplyr::full_join(hla_allele_results$PC_C_pos_0160_A_hla_results,hla_allele_results$Pol_pos_0049_N_hla_results,by=c('ID'='ID')) %>% dplyr::relocate(ID)

kable(hla_allele_results_asn %>% dplyr::arrange(desc(Freq)) %>% dplyr::mutate(ID = sapply(ID,function(x) {
  split_id <- strsplit(x=x,split = '_')[[1]]
  paste0('HLA-',split_id[1],'*',split_id[2],':',gsub("[^0-9.-]", "", split_id[3]))
})) %>% dplyr::rename(`HLA Allele` = ID,Frequency = Freq),format = "html", caption = "", booktabs = T) %>%
  kable_styling(latex_options = "HOLD_position")  %>%
  add_header_above(c(" " = 2, "PreCore/Core position 160 \n (Associated AA: A)" = 3, "Pol position 49 \n (Associated AA: N)" = 3))


kable(hla_allele_results_eur %>% dplyr::arrange(desc(Freq)) %>% dplyr::mutate(ID = sapply(ID,function(x) {
  split_id <- strsplit(x=x,split = '_')[[1]]
  paste0('HLA-',split_id[1],'*',split_id[2],':',gsub("[^0-9.-]", "", split_id[3]))
})) %>% dplyr::rename(`HLA Allele` = ID,Frequency = Freq) %>% dplyr::relocate(`HLA Allele`),format = "html", caption = "", booktabs = T) %>%
  kable_styling(latex_options = "HOLD_position")  %>%
  add_header_above(c(" " = 2, "PreCore/Core position 67 \n (Associated AA: Y)" = 3))


hla_aa_results <- lapply(hla_result,function(x) {
  df <- dplyr::left_join(x$aa_results$uncond,x$aa_results$cond,by = c('ID'='ID','#CHROM'='#CHROM')) %>% 
    dplyr::filter(`#CHROM` == 'A') %>% 
    dplyr::select(ID,P = P.x,P_Cond = P.y) %>% 
    dplyr::mutate(P = scientific(P,digits = 3)) %>% 
    dplyr::mutate(P_Cond = scientific(P_Cond,digits = 3))
  df$Pos <- as.numeric(sapply(df$ID,function(x) strsplit(x=x,split = '_')[[1]][2]))
  AA_String <- lapply(rownames(x$aa_results$uncond %>% dplyr::filter(`#CHROM` == 'A')),function(x) strsplit(x=x,split = '\\.')[[1]])
  AA_String_Parsed <- lapply(AA_String,function(x) sapply(x[x!='chap'],function(y) strsplit(x=y,split = '_')[[1]][3]))
  df$AA <- sapply(AA_String_Parsed,function(x) paste0(x[!is.na(x)],collapse = ','))
  return(df %>% dplyr::select(-ID) %>% dplyr::relocate(Pos,AA) %>% dplyr::arrange(Pos) %>% dplyr::filter(!is.na(P)) %>% dplyr::filter(grepl(',',AA)))
}) 


kable(dplyr::inner_join(hla_aa_results$PC_C_pos_0160_A_hla_results,hla_aa_results$Pol_pos_0049_N_hla_results,by=c('Pos'='Pos','AA'='AA')),format = "html", caption = "", booktabs = T) %>%
  kable_styling(latex_options = "HOLD_position")  %>%
  add_header_above(c(" " = 2, "Precore/core position 160 \n (Associated AA: A)" = 2, "Pol position 49 \n (Associated AA: N)" = 2))

kable(hla_aa_results$PC_C_pos_0067_Y_eur_hla_results,format = "html", caption = "", booktabs = T) %>%
  kable_styling(latex_options = "HOLD_position")  %>%
  add_header_above(c(" " = 2, "Precore/core position 67 \n (Associated AA: Y)" = 2))

## Figure 4 HLA Epitope

source('~/G2G-HBV/src/epitope_analysis.R')
PC_C_160_res <- GetBestBinders("mut_PC_C_pos_0160_A_A_k_","PC_C_pos_0160_A")
best_k_mers_PC_C_160 <- list(list(mut = PC_C_160_res$peptide[which.min(PC_C_160_res$ic50)],non_mut = list(PC_C_160_res$peptide[PC_C_160_res$ic50 != min(PC_C_160_res$ic50)])))
PC_C_160_allele_scores <- GetBestBindersAcrossAlleles("PC_C_pos_0160_A",best_k_mers_PC_C_160,k_mers_PC_C_pos_160_A,'A_33_03',hla_dosage_asn)
PC_C_160_Haplo <- readRDS('~/G2G-HBV/Epitope_Haplotypes/PC_C_160/shotgun_results.rds')
PC_C_160_freq <- GetFreqByAllele(hla_dosage_asn,PC_C_160_Haplo,'A','A_33_03',dat_env_asn,PC_C_160_allele_scores$allele_scores_el_sliding,PC_C_160_allele_scores$peptides_sliding)

PC_C_160_freq_raw <- GetFreqByAllele(hla_dosage_asn,PC_C_160_Haplo,'A','A_33_03',dat_env_asn,PC_C_160_allele_scores$allele_scores_el_sliding,PC_C_160_allele_scores$peptides_sliding,freq_cutoff=-Inf)


PC_C_160_freq_single <- PC_C_160_freq %>% dplyr::filter(Mixed == F,Frq > 0)
PC_C_160_freq_mixed <- PC_C_160_freq %>%  dplyr::filter(Mixed == T,Frq > 0)
PC_C_160_freq_mixed <- PC_C_160_freq_mixed %>% dplyr::left_join(PC_C_160_freq_single,by=c('hla_group'='hla_group','AA'='AA')) %>% 
  dplyr::mutate(Frq = Frq.x + Frq.y) %>% dplyr::select(hla_group,AA,Frq,binding_score = binding_score.x,Mixed = Mixed.x,strong_binder = strong_binder.x,Peptide=Peptide.x)

ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

p1 <- ggplot2::ggplot() + 
  geom_bar(data = PC_C_160_freq_single,aes(x=log10(binding_score),fill = strong_binder,y=Frq),alpha = 1,
           stat = 'identity',position='dodge',width = 0.05) + 
  geom_bar(data = PC_C_160_freq_mixed,aes(x=log10(binding_score),fill = strong_binder,y=Frq),alpha = 0.4,
           stat = 'identity',position='dodge',width = 0.05) +
  facet_grid(~hla_group) + geom_vline(xintercept = log10(1),linetype="dashed",size = 0.6,alpha =0.4) +  annotate("text", label = "Strong Binder",x = -1, y = 1.2, size = 2.5, colour = "black") + 
  annotate("text", label = "Weak Binder",x = 0.5, y = 1.2, size = 2.5, colour = "black") + ylim(0,1.2) + ylab('Fraction of Samples') + xlab('Predicted log10(Elution Percentile Rank)') + scale_y_continuous(breaks = seq(0, 1, len = 5)) + ggrepel::geom_text_repel(aes(x=log10(binding_score),y=Frq,label = Peptide),alpha = 1,direction = "both",size = 3,nudge_y = 0.03,data = rbind(PC_C_160_freq_mixed,dplyr::anti_join(PC_C_160_freq_single,PC_C_160_freq_mixed,by=c('hla_group'='hla_group','AA'='AA')))) + theme(legend.position = "none") + xlim(-2.1,1.1) + scale_fill_manual(values = ggplotColours(2)[2])

Pol_49_res <- GetBestBinders("mut_Pol_pos_0049_N_N_k_","Pol_pos_0049_N")
best_k_mers_Pol_49<- list(list(mut = Pol_49_res$peptide[which.min(Pol_49_res$ic50)],non_mut = list(Pol_49_res$peptide[Pol_49_res$ic50 != min(Pol_49_res$ic50)])))
Pol_49_allele_scores <- GetBestBindersAcrossAlleles("Pol_pos_0049_N",best_k_mers_Pol_49,k_mers_Pol_pos_49_N,'A_02_06',hla_dosage_asn)
Pol_49_Haplo <- readRDS('~/G2G-HBV/Epitope_Haplotypes/Pol_49/shotgun_results.rds')
Pol_49_freq <- GetFreqByAllele(hla_dosage_asn,Pol_49_Haplo,'A','A_02_06',dat_env_asn,Pol_49_allele_scores$allele_scores_el_sliding,Pol_49_allele_scores$peptides_sliding) 

Pol_49_freq_single <- Pol_49_freq %>% dplyr::filter(Mixed == F,Frq > 0)
Pol_49_freq_mixed <- Pol_49_freq %>%  dplyr::filter(Mixed == T,Frq > 0)
Pol_49_freq_mixed <- Pol_49_freq_mixed %>% dplyr::left_join(Pol_49_freq_single,by=c('hla_group'='hla_group','AA'='AA')) %>% 
  dplyr::mutate(Frq = Frq.x + Frq.y) %>% dplyr::select(hla_group,AA,Frq,binding_score = binding_score.x,Mixed = Mixed.x,strong_binder = strong_binder.x,Peptide = Peptide.x)

p2 <- ggplot2::ggplot() + 
  geom_bar(data = Pol_49_freq_single,aes(x=log10(binding_score),fill = strong_binder,y=Frq),alpha = 1,
           stat = 'identity',position='dodge',width = 0.05) + 
  geom_bar(data = Pol_49_freq_mixed,aes(x=log10(binding_score),fill = strong_binder,y=Frq),alpha = 0.4,
           stat = 'identity',position='dodge',width = 0.05) +
  facet_grid(~hla_group) + geom_vline(xintercept = log10(1),linetype="dashed",size = 0.6,alpha =0.4) +  annotate("text", label = "Strong Binder",x = -1, y = 1.2, size = 2.5, colour = "black") + 
  annotate("text", label = "Weak Binder",x = 0.5, y = 1.2, size = 2.5, colour = "black") + ylim(0,1.2) + ylab('Fraction of Samples') + xlab('Predicted log10(Elution Percentile Rank)') + scale_y_continuous(breaks = seq(0, 1, len = 5)) + ggrepel::geom_text_repel(aes(x=log10(binding_score),y=Frq,label = Peptide),alpha = 1,direction = "both",size = 3,nudge_y = 0.03,data = rbind(Pol_49_freq_mixed,dplyr::anti_join(Pol_49_freq_single,Pol_49_freq_mixed,by=c('hla_group'='hla_group','AA'='AA')))) + theme(legend.position = "none") + xlim(-2.1,1.1)



PC_C_67_res <- GetBestBinders("mut_PC_C_pos_0067_Y_Y_k_","PC_C_pos_0067_Y")
best_k_mers_PC_C_67 <- list(list(mut = PC_C_67_res$peptide[which.min(PC_C_67_res$ic50)],non_mut = list(PC_C_67_res$peptide[PC_C_67_res$ic50 != min(PC_C_67_res$ic50)])))
PC_C_67_allele_scores <- GetBestBindersAcrossAlleles("PC_C_pos_0067_Y",best_k_mers_PC_C_67,k_mers_PC_C_pos_67_Y,'A_01_01',hla_dosage_eur)
PC_C_67_Haplo <- readRDS('~/G2G-HBV/Epitope_Haplotypes/PC_C_67/shotgun_results.rds')
PC_C_67_freq <- GetFreqByAllele(hla_dosage_eur,PC_C_67_Haplo,'A','A_01_01',dat_env_eur,PC_C_67_allele_scores$allele_scores_el_sliding,PC_C_67_allele_scores$peptides_sliding)

PC_C_67_freq_single <- PC_C_67_freq %>% dplyr::filter(Mixed == F,Frq > 0)
PC_C_67_freq_mixed <- PC_C_67_freq %>%  dplyr::filter(Mixed == T,Frq > 0)
PC_C_67_freq_mixed <- PC_C_67_freq_mixed %>% dplyr::left_join(PC_C_67_freq_single,by=c('hla_group'='hla_group','AA'='AA')) %>% 
  dplyr::mutate(Frq = Frq.x + Frq.y) %>% dplyr::select(hla_group,AA,Frq,binding_score = binding_score.x,Mixed = Mixed.x,strong_binder = strong_binder.x,Peptide = Peptide.x)

p3 <- ggplot2::ggplot() + 
  geom_bar(data = PC_C_67_freq_single,aes(x=log10(binding_score),fill = strong_binder,y=Frq),alpha = 1,
           stat = 'identity',position='dodge',width = 0.05) + 
  geom_bar(data = PC_C_67_freq_mixed,aes(x=log10(binding_score),fill = strong_binder,y=Frq),alpha = 0.4,
           stat = 'identity',position='dodge',width = 0.05) +
  facet_grid(~hla_group) + geom_vline(xintercept = log10(1),linetype="dashed",size = 0.6,alpha =0.4) +  annotate("text", label = "Strong Binder",x = -1, y = 1.2, size = 2.5, colour = "black") + 
  annotate("text", label = "Weak Binder",x = 0.5, y = 1.2, size = 2.5, colour = "black") + ylim(0,1.2) + ylab('Fraction of Samples') + xlab('Predicted log10(Elution Percentile Rank)') + scale_y_continuous(breaks = seq(0, 1, len = 5)) + ggrepel::geom_text_repel(aes(x=log10(binding_score),y=Frq,label = Peptide),alpha = 1,direction = "both",size = 3,nudge_y = 0.03,data = rbind(PC_C_67_freq_mixed,dplyr::anti_join(PC_C_67_freq_single,PC_C_67_freq_mixed,by=c('hla_group'='hla_group','AA'='AA')))) + theme(legend.position = "none") + xlim(-2.1,1.1)


kable(rbind(PC_C_160_res,Pol_49_res,PC_C_67_res) %>% dplyr::filter(peptide %in% c('WIRTPPAYR','VSIPWTHKV','LLDTASALY')) %>% dplyr::select(`Peptide`=peptide,`MixMHCPred Elution Percentile Rank`=perc_rank_MixMHCPred,`NetMHCPan IC50`=ic50,`NetMHCPan Elution Percentile Rank`=rank),format = "html") %>%
  kableExtra::group_rows('Position 160 precore/core',1,1) %>%
  kableExtra::group_rows('Position 49 polymerase',2,2) %>% 
  kableExtra::group_rows('Position 67 precore/core',3,3)

## Suppl Fig 1 - 3 pPCA and PCA

## -------------------------------
## After running src/mapping-1kg.R
## -------------------------------
DIR_SCRATCH <- '~/G2G-HBV/data/results_POP_asian_GT_A_C_D_B_TYPE_pPCA_LOCO_FALSE/'
FILE_REF <- glue::glue("{DIR_SCRATCH}/mapping-1kg/1KG_reference") 
FILE_TARGET <- glue::glue("{DIR_SCRATCH}/mapping-1kg/hbv_gilead_QC_PCA_TRANS") 


idlist_pop_REF <- data.table::fread("/mnt/data2/naret/work/data/1KG/other_data/integrated_call_samples_v3.20130502.ALL.panel") %>% dplyr::select(sample, pop, super_pop)
idlist_pop_TARGET <- readr::read_delim(glue::glue("{DIR_SCRATCH}/idlist_pop.txt"), delim = " ", col_names = FALSE) 

pc_REF <- data.table::fread(glue::glue("{FILE_REF}_pca20.eigenvec")) %>% dplyr::select(-V1)
names(pc_REF) <- c("ID", paste0("PC", 1:20))
pc_TARGET <- data.table::fread(glue::glue("{FILE_TARGET}_pca20.proj.eigenvec"))%>% dplyr::select(-V1)
names(pc_TARGET) <- c("ID", paste0("PC", 1:21))

## merge labels to pc_ref ----------------------
pc_REF <- pc_REF %>% dplyr::left_join(idlist_pop_REF, by = c("ID" = "sample" )) %>% 
  dplyr::rename(`1KG`= super_pop)

## merge labels of asian cluster ---------------
pc_TARGET <- pc_TARGET %>% dplyr::mutate(hbv = ID %in% idlist_pop_TARGET$X1) %>% dplyr::left_join(dat_covars_raw %>% dplyr::select(host_id, GT,RACE), by = c("ID"= "host_id")) %>% dplyr::filter(!is.na(RACE)) %>% dplyr::filter(RACE %in% c('ASIAN','WHITE')) #%>% dplyr::filter(ID %in% ids_unrelated$X1)
pc_TARGET$`Self Reported Ancestry` <- forcats::fct_recode(as.factor(pc_TARGET$hbv), Asian = 'TRUE', `Non Asian` = 'FALSE' )

n_asian <- length(idlist_pop_TARGET$X1)

## plot ----------------------------------------
pc_TARGET <- pc_TARGET %>% dplyr::left_join(dplyr::select(tmp_cluster,ID=host_id,gr_cluster))
pc_TARGET$RACE[pc_TARGET$RACE == 'WHITE'] <- 'Self-Reported:European'
pc_TARGET$RACE[pc_TARGET$RACE == 'ASIAN'] <- 'Self-Reported:Asian'

pc_TARGET$RACE <- factor(pc_TARGET$RACE)
p1 <- ggplot(data = pc_TARGET) + 
  geom_point(data = pc_REF, aes(PC1, PC2, color = `1KG`)) + 
  geom_point(data = pc_TARGET, aes(PC1, PC2),  alpha = I(0.5)) +
  facet_wrap(~RACE,nrow = 2)
#theme(legend.position = "none")
p2 <- ggplot(data = pc_TARGET) + 
  geom_point(data = pc_REF, aes(PC3, PC4, color = `1KG`)) + 
  geom_point(data = pc_TARGET, aes(PC3, PC4),  alpha = I(0.5)) +
  facet_wrap(~RACE,nrow = 2)
p1_p2 <- ggpubr::ggarrange(p1,p2,common.legend = T)

#Get Table of percent variance explained 
FILE_GRM_OUT <- '~/G2G-HBV/data/processed_POP_asian_GT_A_C_D_B_TYPE_pPCA_LOCO_FALSE/host-grm'
eigen_val <- readr::read_table2(glue::glue("{FILE_GRM_OUT}_cluster.eigenval"), col_names = FALSE)
perc_var_explained <- round(t(eigen_val / sum(eigen_val) * 100),2)
load('~/G2G-HBV/data/results_POP_asian_european_GT_A_C_D_B_F_H_TYPE_pPCA_LOCO_FALSE/prep-data.rda')
asn_PCs <- dplyr::filter(dat_covars_numeric,host_id %in% asn_cluster$host_id)
eur_PCs <- dplyr::filter(dat_covars_numeric,host_id %in% eur_cluster$host_id)
asn_PCs$`Assigned Ancestry` <- 'East Asian'
eur_PCs$`Assigned Ancestry` <- 'European'
merged_PCs <- rbind(asn_PCs,eur_PCs)

p3 <- ggplot(data = merged_PCs) + aes(x=host_PC1,y=host_PC2,color = `Assigned Ancestry`) + geom_point() + xlab(paste0('PC1 (',perc_var_explained[1],'%)')) + ylab(paste0('PC2 (',perc_var_explained[2],'%)'))
p4 <- ggplot(data = merged_PCs) + aes(x=host_PC3,y=host_PC4,color = `Assigned Ancestry`) + geom_point() + xlab(paste0('PC3 (',perc_var_explained[3],'%)')) + ylab(paste0('PC4 (',perc_var_explained[4],'%)'))
p3_p4 <- ggpubr::ggarrange(p3,p4,common.legend = T)

host_PC_plot <- ggpubr::ggarrange(p1_p2,p3_p4,nrow=2,labels = c('A)','B)'),heights = c(1,0.6))

load('~/G2G-HBV/data/results_POP_asian_european_GT_A_C_D_B_F_H_TYPE_pPCA_LOCO_FALSE/prep-data.rda')
source('./scripts/misc/plot-phyla-tree.R')
pPC_eigen <- PlotPCA('~/G2G-HBV/data/results_POP_asian_european_GT_A_C_D_B_F_H_TYPE_pPCA_LOCO_FALSE/prep-data.rda','Asian_European','~/G2G-HBV/data/phylo_tree/HBV_WG+og_rr.nw','~/G2G-HBV/data/raw/gilead_20181126/viral_seq/OUT_aa_binary_table.txt')$pPC_eigen
var_explained <- round(pPC_eigen/sum(pPC_eigen)*100,2)
pPC_raw <- dat_covars_numeric %>% dplyr::filter(!is.na(pPC1)) %>% dplyr::left_join(dat_covars_raw)
pPC <- pPC_raw %>% dplyr::select(contains('pPC')) 
pPC$`HBV Genotype` <- pPC_raw$GT
p1 <- ggplot(data = pPC) + 
  geom_point(data = pPC, aes(pPC1, pPC2, color = `HBV Genotype`)) + xlab(paste0('pPC1 (',var_explained[1],'%)')) + ylab(paste0('pPC2 (',var_explained[2],'%)'))
p2 <- ggplot(data = pPC) + 
  geom_point(data = pPC, aes(pPC3, pPC4, color = `HBV Genotype`)) + xlab(paste0('pPC3 (',var_explained[3],'%)')) + ylab(paste0('pPC4 (',var_explained[4],'%)'))
p3 <- ggplot(data = pPC) + 
  geom_point(data = pPC, aes(pPC5, pPC6, color = `HBV Genotype`)) + xlab(paste0('pPC5 (',var_explained[5],'%)')) + ylab(paste0('pPC6 (',var_explained[6],'%)'))
ggpubr::ggarrange(plotlist = list(p1,p2,p3),nrow = 2,ncol=2,common.legend = T)

## Suppl Fig 4-5 QQ
source('./scripts/G2G/results-g2g.R')
library(grid)
library(png)
load('~/G2G-HBV/data/results_POP_asian_GT_A_C_D_B_TYPE_pPCA_LOCO_FALSE/prep-data.rda')
data_dir_asn <- '~/G2G-HBV/data/results_POP_asian_GT_A_C_D_B_TYPE_pPCA_LOCO_FALSE/'
data_dir_eur <- '~/G2G-HBV/data/results_POP_european_GT_A_C_D_F_H_TYPE_pPCA_LOCO_FALSE/'

ntcp_plots <- lapply(paste0(data_dir_asn,c('results_top_qqplot_SAIGE_gene_S_pos_0051_P.png',
                                           'results_top_qqplot_SAIGE_gene_S_pos_0035_R.png',
                                           'results_top_qqplot_SAIGE_gene_S_pos_0017_A.png',
                                           'results_top_qqplot_SAIGE_gene_Pol_pos_0197_C.png')),function(x){
                                             img <- as.raster(readPNG(x))
                                             rasterGrob(img, interpolate = FALSE)
                                           })

hla_plots <- lapply(c(paste0(data_dir_eur,'results_top_qqplot_SAIGE_gene_PC_C_pos_0067_Y.png'),
                      paste0(data_dir_asn,c('results_top_qqplot_SAIGE_gene_Pol_pos_0049_N.png',
                                            'results_top_qqplot_SAIGE_gene_PC_C_pos_0160_A.png'))),function(x){
                                              img <- as.raster(readPNG(x))
                                              rasterGrob(img, interpolate = FALSE)
                                            })

ntcp_manhattan <- lapply(paste0(data_dir_asn,c('results_top_manhattan_SAIGE_gene_S_pos_0051_P.png',
                                               'results_top_manhattan_SAIGE_gene_S_pos_0035_R.png',
                                               'results_top_manhattan_SAIGE_gene_S_pos_0017_A.png',
                                               'results_top_manhattan_SAIGE_gene_Pol_pos_0197_C.png')),function(x){
                                                 img <- as.raster(readPNG(x))
                                                 rasterGrob(img, interpolate = FALSE)
                                               })

hla_manhattan <- lapply(c(paste0(data_dir_eur,'results_top_manhattan_SAIGE_gene_PC_C_pos_0067_Y.png'),
                          paste0(data_dir_asn,c('results_top_manhattan_SAIGE_gene_Pol_pos_0049_N.png',
                                                'results_top_manhattan_SAIGE_gene_PC_C_pos_0160_A.png'))),function(x){
                                                  img <- as.raster(readPNG(x))
                                                  rasterGrob(img, interpolate = FALSE)
                                                })


library(gridExtra)

ggsave("~/G2G-HBV/Final_Figures/ntcp_qq.pdf",width=8, height=8, 
       marrangeGrob(grobs = ntcp_plots, nrow=2, ncol=2,top=NULL))
ggsave("~/G2G-HBV/Final_Figures/hla_qq.pdf",width=8, height=8, 
       marrangeGrob(grobs = hla_plots, nrow=2, ncol=2,top=NULL))
ggsave("~/G2G-HBV/Final_Figures/ntcp_manhattan.pdf",width=4, height=8, 
       marrangeGrob(grobs = ntcp_manhattan, nrow=4, ncol=1,top=NULL))
ggsave("~/G2G-HBV/Final_Figures/hla_manhattan.pdf",width=4, height=6, 
       marrangeGrob(grobs = hla_manhattan, nrow=3, ncol=1,top=NULL))

## Suppl Fig 6 DnDs
source('./scripts/DnDs/CalculateDnDs.R')
library(ggplot2)
library(latex2exp)
library(ggpmisc)
p1 <- ggplot2::ggplot(data = merged_dnds %>% dplyr::filter(!is.na(dSN))) + aes(x=rs2296651,y=dSN) + geom_boxplot(outlier.shape=NA) +
  geom_jitter(position=position_jitter(width=.1, height=0),alpha = 0.5) + ylab(TeX("$\\pi_{SN}$")) + stat_n_text() + stat_compare_means(method = 'wilcox.test',size = 3) +ylim(NA,0.115)
p2 <- ggplot2::ggplot(data = merged_dnds %>% dplyr::filter(!is.na(dSN))) + aes(x=rs2296651,y=dNN) + geom_boxplot(outlier.shape=NA) +
  geom_jitter(position=position_jitter(width=.1, height=0),alpha = 0.5) + ylab(TeX("$\\pi_{NN}$")) + stat_n_text() + stat_compare_means(method = 'wilcox.test',size = 3) +ylim(NA,0.07)


p3 <- ggplot2::ggplot(data = merged_dnds %>% dplyr::filter(!is.na(dSN) & dSN != 0)) + aes(x=rs2296651,y=dNN/dSN) + geom_boxplot(outlier.shape=NA) +
  geom_jitter(position=position_jitter(width=.1, height=0),alpha = 0.5) + ylab(TeX("$\\frac{\\pi_{NN}}{\\pi_{SN}}$")) + stat_n_text() + stat_compare_means(method = 'wilcox.test',size = 3) +ylim(NA,6)
p_piNN <- ggarrange(p1,p2,p3,ncol = 3)


p11 <- ggplot2::ggplot(data = merged_dnds %>% dplyr::filter(!is.na(dS))) + aes(x=rs2296651,y=dS) + geom_boxplot(outlier.shape=NA) +
  geom_jitter(position=position_jitter(width=.1, height=0),alpha = 0.5) + ylab(TeX("$\\pi_{S}$")) + stat_n_text() + stat_compare_means(method = 'wilcox.test',size = 3) +ylim(NA,0.115)
p22 <- ggplot2::ggplot(data = merged_dnds %>% dplyr::filter(!is.na(dS))) + aes(x=rs2296651,y=dN) + geom_boxplot(outlier.shape=NA) +
  geom_jitter(position=position_jitter(width=.1, height=0),alpha = 0.5) + ylab(TeX("$\\pi_{N}$")) + stat_n_text() + stat_compare_means(method = 'wilcox.test',size = 3) +ylim(NA,0.07)
p33 <- ggplot2::ggplot(data = merged_dnds %>% dplyr::filter(!is.na(dS) & dS != 0)) + aes(x=rs2296651,y=dN/dS) + geom_boxplot(outlier.shape=NA) +
  geom_jitter(position=position_jitter(width=.1, height=0),alpha = 0.5) + ylab(TeX("$\\frac{\\pi_{N}}{\\pi_{S}}$")) + stat_n_text() + stat_compare_means(method = 'wilcox.test',size = 3) +ylim(NA,6)

p_piN <- ggarrange(p11,p22,p33,ncol = 3)


ggarrange(p_piN+theme(plot.margin=unit(c(0.5,1,0,0.5),"cm")),
          p_piNN +theme(plot.margin=unit(c(0.5,1,0,0.5),"cm")),nrow = 2,labels = c('A)','B)'))

## Suppl Fig 8 Evolution Trajectory
library(phylotools)
library(stats)
library(dplyr)
library(ggplot2)
library(ggtree)
library(ggrepel)

PlotEvolTrajectory <- function(fasta_file,name){
  fasta <- seqinr::read.fasta(file=fasta_file)
  #Filter for Posterior > 0.95
  posterior <- as.numeric(sapply(names(fasta),function(x) strsplit(x=x,split = 'posterior=')[[1]][2]))
  reads <- round(as.numeric(sapply(phylotools::get.fasta.name(fasta_file),function(x) strsplit(x=x,split='reads=')[[1]][2]),1))
  fasta <- fasta[posterior > 0.95 & reads > 10]
  seqinr::write.fasta(sequences = fasta,names = phylotools::get.fasta.name(fasta_file)[posterior > 0.95 & reads > 10],file.out = gsub(fasta_file,pattern = '.fas',replacement = '_high_conf.fas'))
  #Get PreS1 haplotypes based on G2G associated residues
  haplotypes <- sapply(fasta,function(x) paste0(seqinr::translate(x)[c(17,35,51)],collapse = ''))
  dna <- adegenet::fasta2DNAbin(file=gsub(fasta_file,pattern = '.fas',replacement = '_high_conf.fas'))
  #Get distance matrix
  D <- ape::dist.dna(dna, model = "TN93")
  #Construct neighbour-joining tree
  tre <- nj(D)
  # dna2 <- as.phyDat(dna)
  # tree.pars <- optim.parsimony(tre, dna2)
  
  #Label tips
  tre$tip.label <- paste0(haplotypes,'(',round(reads[posterior > 0.95 & reads > 10] / sum(reads[posterior > 0.95 & reads > 10]) * 100,0),'%)')
  
  p <- ggtree(tre,layout = 'unrooted')  %<+% data.frame(ID = tre$tip.label,Haplotype = haplotypes) +
    geom_label_repel(aes(label=label),cex = 2.5) + geom_tippoint(aes(color=Haplotype)) + theme(legend.position = "none")
  return(p)
  # plot(tre,type = 'unrooted',cex = 0.7,no.margin = T)
  # myPal <- colorRampPalette(c("red",'orange',"yellow","green",'lightblue','violet'))
  # tiplabels('', bg = fac2col(factor(haplotypes), col.pal = myPal), cex = 0.2)
  
  
  # msaPrettyPrint(msa(DNAStringSet(sapply(fasta,function(x) paste0(x,collapse = '')))),file = paste0(name,'.pdf'))
}
p1 <- PlotEvolTrajectory('~/G2G-HBV/MyrcludexB_Full_Length_Haplotypes/GS-US-320-0110-4571/shotgun_results.fas',name = 'GS-US-320-0110-4571')
# PlotEvolTrajectory('~/G2G-HBV/MyrcludexB_Full_Length_Haplotypes/GS-US-320-0110-4901/shotgun_results.fas',name = 'GS-US-320-0110-4901')
p2 <- PlotEvolTrajectory('~/G2G-HBV/MyrcludexB_Full_Length_Haplotypes/GS-US-283-1062-3071/shotgun_results.fas',name = 'GS-US-283-1062-3071')
# PlotEvolTrajectory('~/G2G-HBV/MyrcludexB_Full_Length_Haplotypes/GS-US-320-0110-5266/shotgun_results.fas',name = 'GS-US-320-0110-5266')
p3 <- PlotEvolTrajectory('~/G2G-HBV/MyrcludexB_Full_Length_Haplotypes/GS-US-320-0110-4586/shotgun_results.fas',name = 'GS-US-320-0110-4586')
# PlotEvolTrajectory('~/G2G-HBV/MyrcludexB_Full_Length_Haplotypes/GS-US-320-0110-4663/shotgun_results.fas',name = 'GS-US-320-0110-4663')

p4 <- PlotEvolTrajectory('~/G2G-HBV/MyrcludexB_Full_Length_Haplotypes/GS-US-320-0108-1095/shotgun_results.fas',name = 'GS-US-320-0108-1095')
# PlotEvolTrajectory('~/G2G-HBV/MyrcludexB_Full_Length_Haplotypes/GS-US-320-0110-4529/shotgun_results.fas',name = 'GS-US-320-0110-4529')
ggpubr::ggarrange(p1,p2,p3,p4,nrow = 2,ncol = 2)

## Suppl Fig 9 - pGWAS
gwas_res <- readRDS('~/G2G-HBV/data/pGWAS/pGWAS_nt_pPCA_PCA_result.rds')
merged_pGWAS_result <- gwas_res$merged_pGWAS_result

threshold <- 0.05 / sum(!is.na(merged_pGWAS_result$pGWAS_surface_antigen$p))

pGWAS_sa <- merged_pGWAS_result$pGWAS_surface_antigen %>% dplyr::group_by(Pos) %>% dplyr::arrange(p) %>% dplyr::filter(row_number()==1)
pGWAS_vl <- merged_pGWAS_result$pGWAS_viral_load %>% dplyr::group_by(Pos) %>% dplyr::arrange(p) %>% dplyr::filter(row_number()==1)
pGWAS_ALT <- merged_pGWAS_result$pGWAS_ALT %>% dplyr::group_by(Pos) %>% dplyr::arrange(p) %>% dplyr::filter(row_number()==1)
pGWAS_HBeAg <- merged_pGWAS_result$pGWAS_HBEAG_df %>% dplyr::group_by(Pos) %>% dplyr::arrange(p) %>% dplyr::filter(row_number()==1)

#Pathogen GWAS
pGWAS <- rbind(data.frame(pGWAS_sa,Pheno = 'HBsAg Levels'),data.frame(pGWAS_vl,Pheno = 'Viral Load'),data.frame(pGWAS_ALT,Pheno = 'Serum ALT'),data.frame(pGWAS_HBeAg,Pheno = 'HBeAg Status'))
ggplot(pGWAS) + aes(x=Pos,y=-log10(p),color = r2) + geom_point() + geom_hline(colour = 'red',yintercept = -log10(threshold)) + ggrepel::geom_label_repel(data = pGWAS %>% dplyr::filter((Pheno != 'HBeAg Status' & p < threshold) | (Pheno == 'HBeAg Status' & Pos %in% c(1762,1764,1838,1817,1896)) ),
                                                                                                                                                         aes(label=as.factor(paste0(gsub('pos_','',SNP),' (beta=',round(beta,2),')'))),size=2.5,color = 'black') + facet_grid(rows = vars(Pheno),scales = 'free_y') + xlab('HBV Nucleotide Position') + scale_color_viridis()

