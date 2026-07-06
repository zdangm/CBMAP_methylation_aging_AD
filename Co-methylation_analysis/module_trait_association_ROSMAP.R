library(data.table)
library(dplyr)
library(rhdf5)
library(readxl)
library(WGCNA)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(stringr)
setwd("/methylation/methy_AD_revision2/results/WGCNA")

load('ROSMAP_pheno.RData')
pheno_df$projid <- as.character(str_pad(pheno_df$projid, width = 8, pad = "0", side = "left"))
hml_info <- read_excel('dataset_1495_cross-sectional_02-16-2025.xlsx')
metadata_long <- read_excel('dataset_1495_long_12-01-2024.xlsx') %>% as.data.frame()
PD_case <- metadata_long %>% filter(
  r_pd %in% c(1)) %>% pull(projid) %>% unique()
PD_0 <- unique(hml_info$projid[which(! hml_info$projid %chin% PD_case)])
LBD_0 <- hml_info %>% filter(
  dlbdx %in% c(0)) %>% pull(projid) %>% unique()
cogdx_0 <- hml_info %>% filter(
  cogdx %in% c(1,2,3)) %>% pull(projid) %>% unique()
case = hml_info %>% filter(
  ADNC_4level %in% c(1,2,3)) %>% pull(projid) %>% unique()
control = hml_info %>% filter(
  ADNC_4level %in% c(0) 
  & (projid %in% PD_0) 
  & (projid %in% LBD_0)
  & (projid %in% cogdx_0)) %>% pull(projid) %>% unique()
sample = intersect(c(case, control), pheno_df$projid) 
pheno_df = pheno_df[which(pheno_df$projid %chin% sample),]
pheno_df$age_death <- as.numeric(unlist(hml_info[match(pheno_df$projid, hml_info$projid), 'age_death']))
pheno_df$AD <- as.numeric(unlist(hml_info[match(pheno_df$projid, hml_info$projid), 'ADNC_4level']))
pheno_df$Braak.NFT.stage <- as.numeric(unlist(hml_info[match(pheno_df$projid, hml_info$projid), "braaksc"]))
pheno_df$A_beta_0_3 <- as.numeric(unlist(hml_info[match(pheno_df$projid, hml_info$projid), "a_score_st4"]))
pheno_df$C_cerad_0_3 <- as.numeric(unlist(hml_info[match(pheno_df$projid, hml_info$projid), "c_score"]))
pheno_df$msex <- as.factor(pheno_df$msex)
pheno_df$Sample_Plate <- as.factor(pheno_df$Sample_Plate)
pheno_df$batch <- as.factor(pheno_df$batch)
pheno_df$Study <- as.factor(pheno_df$Study)
pheno_df <- pheno_df[!is.na(pheno_df$AD),] 

idkey = read.csv('ROSMAP_IDkey.csv', header=T)
rin_info = read.csv('ROSMAP_assay_rnaSeq_metadata.csv',header=T)
rin_info = rin_info[,c("specimenID","RIN")]
rin_info$projid = unlist(idkey[match(rin_info$specimenID,idkey$rnaseq_id),"projid"])
rin_info$projid <- as.character(str_pad(rin_info$projid, width = 8, pad = "0", side = "left"))
pheno_df$rin = unlist(rin_info[match(pheno_df$projid,rin_info$projid),"RIN"]) 

cell.pro = read.table('ROSMAP_all_sample_methylation_ctp_with_clr_deconvolution.txt',header=T)
cell_aligned <- cell.pro[pheno_df$barcode_id, , drop = FALSE]
pheno_df <- cbind(pheno_df, cell_aligned)
rownames(pheno_df) = pheno_df$mwas_id

# prepare methylation matrix
cg <- h5read("/methylation/data/ROSMAP_methylaition/final_DNAm_invMdat.h5",'Mdat')
colnames(cg) = h5read("/methylation/data/ROSMAP_methylaition/final_DNAm_invMdat.h5", 'colnames')
rownames(cg) = h5read("/methylation/data/ROSMAP_methylaition/final_DNAm_invMdat.h5", 'rownames')
sample = intersect(rownames(cg), pheno_df$mwas_id) # 680 samples
cg = cg[sample,]
pheno_df = pheno_df[sample,]
sample_infor = pheno_df[,c("projid","mwas_id","msex","age_death","Sample_Plate","batch","Study","pmi","AD","Braak.NFT.stage","A_beta_0_3","C_cerad_0_3",     
                           "rin","Exc","Inh","NonN_Astro_FGF3R","NonN_Endo","NonN_Micro","NonN_Oligo_MBP","NonN_OPC")]
for (i in c(8,13:20)) {
  sample_infor[is.na(sample_infor[, i]), i] <- mean(sample_infor[, i], na.rm = TRUE)
}
identical(rownames(sample_infor), rownames(cg))

load('meta_WGCNA_res.RData')
adnc_lim_all <- list(
  net = net,
  power =  power_selected,
  moduleColors = moduleColors,
  color_table = color_table,
  MEs = MEs
)
cg = cg[,names(net$colors)]
adnc_lim_all$MEs <- moduleEigengenes(cg, net$colors)$eigengenes

lm_res = data.frame('me'=NA, 'trait'=NA, 'B'=NA, 'P.Value'=NA)
module_trait_cor <- function(data){
  me <- data$MEs[,2:ncol(data$MEs)]
  
  trait = sample_infor
  colnames(trait) = colnames(sample_infor)
  trait$Sample_Plate = factor(trait$Sample_Plate)
  trait$Study = factor(trait$Study)
  
  identical(rownames(me),rownames(trait))
  for (pheno in c('age_death','AD','Braak.NFT.stage','A_beta_0_3','C_cerad_0_3')) {
    for (i in 1:ncol(me)) {
      if (pheno %chin% c('AD','Braak.NFT.stage','A_beta_0_3','C_cerad_0_3')){
        fit <- lm(me[,i] ~ trait[,pheno] + msex + age_death + Study + Sample_Plate + pmi + rin + Exc + Inh + NonN_Astro_FGF3R + NonN_Endo + NonN_Micro + NonN_Oligo_MBP, data=trait)
      } else {
        fit <- lm(me[,i] ~ trait[,pheno] + msex + Study + Sample_Plate + pmi + rin + Exc + Inh + NonN_Astro_FGF3R + NonN_Endo + NonN_Micro + NonN_Oligo_MBP, data=trait)
      }
      res <- data.frame('me'=colnames(me)[i], 'trait'=pheno,
                        'B'=summary(fit)$coefficients[2,1],
                        'P.Value'=summary(fit)$coefficients[2,4])
      lm_res = rbind(lm_res, res)
    }
  }
  lm_res = lm_res[-1,]
  lm_res <- lm_res %>%
    group_by(trait) %>%
    mutate(adj.P.Val = p.adjust(P.Value, method = "BH")) %>%
    ungroup()
  return(lm_res)
}
module_trait_all <- module_trait_cor(adnc_lim_all)
module_trait_all$trait <-ifelse(module_trait_all$trait == 'AD','ADNC',
                                ifelse(module_trait_all$trait == 'age_death','Age',
                                       ifelse(module_trait_all$trait == 'A_beta_0_3','A-beta',
                                              ifelse(module_trait_all$trait == 'C_cerad_0_3','C-score',module_trait_all$trait))))
module_trait_all$pmarker <- ifelse(module_trait_all$P.Value < 0.05, "*", '')
module_trait_all$me <- factor(module_trait_all$me, levels = paste0('ME', 1:99))
module_trait_all$trait <- factor(module_trait_all$trait, levels = c('Age','ADNC','Braak.NFT.stage','A-beta','C-score'))
df <- module_trait_all %>%
  mutate(
    pmarker = case_when(
      P.Value < 0.001 ~ "***",
      P.Value < 0.01 ~ "**",
      P.Value < 0.05 ~ "*",
      TRUE ~ ""
    ),
    direction = ifelse(B > 0, "pos", ifelse(B < 0, "neg", "zero")),
    abs_B = abs(B)
  )
mat <- reshape2::acast(df, me ~ trait, value.var = "B")
label_mat <- reshape2::acast(df, me ~ trait, value.var = "pmarker")
mat_t <- t(mat)
label_mat_t <- t(label_mat)
col_fun <- colorRamp2(c(min(mat), 0, max(mat)), c("#4DBAD6", "white", "#E44A33"))
Heatmap(
  mat_t,
  name = "Estimate",
  col = col_fun,
  row_names_gp = gpar(fontsize = 13),
  column_names_gp = gpar(fontsize = 13),
  column_names_rot = 45,
  row_names_rot = 45,
  row_names_side = "left",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(label_mat_t[i, j], x, y, gp = gpar(fontsize = 13))
  },
  heatmap_legend_param = list(
    title = "Estimate",
    title_gp = gpar(fontsize = 15, fontface = "plain"), 
    labels_gp = gpar(fontsize = 12),                    
    legend_height = unit(5, "cm"),                     
    legend_width = unit(3, "cm")
  )
)
openxlsx::write.xlsx(module_trait_all,'/methylation/methy_AD_revision2/results/WGCNA/module_trait_association_ROSMAP.xlsx')
