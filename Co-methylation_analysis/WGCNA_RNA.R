library(rhdf5)
library(readxl)
library(data.table)
library(impute)
library(WGCNA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(writexl)
library(dplyr)
library(forcats)
library(stringr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(gridExtra)
library(scales)
library(AnnotationDbi)
library(openxlsx)
library(RColorBrewer)
library(Rcpp)
library(RcppArmadillo)
library(tidyverse)
library(foreach)
library(doParallel)
allowWGCNAThreads(nThreads = 30)
enableWGCNAThreads(nThreads = 30)

load('meta_hml4level.RData')
all.cpg = meta_final_ordered_df$cpg[which(! is.na(meta_final_ordered_df$pVal.final))]

setwd("/methylation/methy_AD_revision2/results/WGCNA")
load('/methylation/results/disease_meth_age/AD_meth/CBMAP_pheno.RData')
all_sample_info <- as.data.frame(read.csv("CBMAP_sample_info.csv",header=T)) 
selected_sample = all_sample_info[!is.na(all_sample_info[,'ADNC']),]
selected_sample <- selected_sample %>% 
  filter(
    diag_ALS == 0 & 
      diag_SCZ == 0 & 
      diag_epilepsy == 0 &
      diag_others == 0 &
      OTHER == 0) 
complete_rows <- complete.cases(selected_sample[, c("ADNC", "ADNC.L", "ADNC.M", "ADNC.H")])
selected_sample <- selected_sample[complete_rows, ]
selected_sample$marker <- ifelse(selected_sample$ADNC == 0, 0,
                                 ifelse(selected_sample$ADNC.L == 1, 1,
                                        ifelse(selected_sample$ADNC.M == 1, 2,
                                               ifelse(selected_sample$ADNC.H == 1, 3, NA))))
exclued_id <- selected_sample[which(selected_sample$marker == 0 & (selected_sample$LBD == 1 | selected_sample$ARTAG_Gray.matter == 1 |
                                                                     selected_sample$ARTAG_Perivascular == 1 | selected_sample$ARTAG_Subpial == 1 | selected_sample$diag_PD == 1 | selected_sample$diag_dementia == 1)),'id']
pheno_df$AD <- selected_sample[match(pheno_df$Sample_Name, selected_sample$id), 'marker']
pheno_df <- pheno_df[!is.na(pheno_df$AD),] 
pheno_df <- pheno_df[which(! pheno_df$Sample_Name %chin% exclued_id),]

# methylation matrix
Mdat <- h5read("/methylation/data/CBMAP/DNAm_processed/final_DNAm_Mdat.h5",'Mdat')
colnames(Mdat) = h5read("/methylation/data/CBMAP/DNAm_processed/final_DNAm_Mdat.h5",'colnames')
rownames(Mdat) = h5read("/methylation/data/CBMAP/DNAm_processed/final_DNAm_Mdat.h5",'rownames')
Mdat = Mdat[grep("cg",rownames(Mdat)),] 
cg <- 2^Mdat / (2^Mdat + 1)
cg <- t(cg)

sample = intersect(rownames(cg), rownames(pheno_df)) 
cg <- cg[sample,]
pheno_df <- pheno_df[sample,]
pheno_df$sex_male <- as.factor(pheno_df$sex_male)
pheno_df$Sample_Plate <- as.factor(pheno_df$Sample_Plate)
pheno_df$bank <- as.factor(pheno_df$bank)

cell.pro = read.table('CBMAP_all_sample_methylation_ctp_with_clr_deconvolution.txt',header=T)
cell_aligned <- cell.pro[pheno_df$barcode_id, , drop = FALSE]
pheno_df <- cbind(pheno_df, cell_aligned)
pheno_df[] <- lapply(pheno_df, function(x)
  if (is.numeric(x)) replace(x, is.na(x), median(x, na.rm = TRUE)) else x
)
identical(rownames(pheno_df), rownames(cg))

# extract high variance cpg
var_cpg <- apply(cg, 2, var)
threshold_lowest_20 <- quantile(var_cpg, 0.2) 
cpg_ids <- names(var_cpg)[var_cpg >= threshold_lowest_20]
cpg_ids <- unique(intersect(all.cpg,cpg_ids))
eqtm_res = fread('/eQTM/cis_eQTM_all_res_clr.txt',header=T)
eqtm_res0 = eqtm_res[fdr_global<0.05]
gene_ids = unique(eqtm_res0$gene_id[which(eqtm_res0$cpg %chin% cpg_ids)])
exp = read.table('/CBMAP_RNAseq_protein_coding/RNAseq_process/RNAseq_data_qn_comb.txt',check.names = F)
exp = t(exp)
exp = exp[,gene_ids]
sample = intersect(sample,rownames(exp))
exp = exp[sample,]

########### implement WGCNA ##############
gsg <- goodSamplesGenes(exp)
if (!gsg$allOK) {
  exp <- exp[gsg$goodSamples, gsg$goodGenes]
}

powers = c(1:20)
sft <- pickSoftThreshold(exp,
                         networkType = "signed",
                         corFnc = "bicor",
                         powerVector = powers,
                         verbose = 5)
pdf('/methylation/methy_AD_revision2/results/WGCNA/power_selected_RNA.pdf',width=9,height=5)
par(mfrow = c(1, 2))
cex1 = 0.9

plot(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  xlab = "Soft Threshold (power)",
  ylab = "Scale Free Topology Model Fit, signed R²",
  type = "n",
  main = "Scale Independence"
)
text(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  labels = powers,
  cex = cex1,
  col = "red"
)
abline(h = 0.8, col = "blue", lty = 2)  

plot(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  xlab = "Soft Threshold (power)",
  ylab = "Mean Connectivity",
  type = "n",
  main = "Mean Connectivity"
)
text(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  labels = powers,
  cex = cex1,
  col = "red"
)
abline(h = 500, col = "blue", lty = 2)  
dev.off()
save(sft,file='/methylation/methy_AD_revision2/results/WGCNA/sft_RNA.RData')

sft_data <- data.frame(Power = sft$fitIndices[, 1], R2 = -sign(sft$fitIndices[,3])*sft$fitIndices[,2], mean_connectivity = sft$fitIndices[, 5])
power_selected <- sft_data$Power[which(sft_data$R2 > 0.8 & sft_data$mean_connectivity < 500)[1]]
print(power_selected)

cor <- WGCNA::cor
net <- blockwiseModules(
  exp,
  power = power_selected,
  corType = "bicor",         
  networkType = "signed",     
  TOMType = "signed",        
  deepSplit = 2,             
  pamRespectsDendro = FALSE,   
  pamStage = FALSE,
  minModuleSize = 30,       
  mergeCutHeight = 0.25,     
  numericLabels = TRUE,      
  verbose = 3,
  nThreads = 30
)

table(net$colors)

MEs <- moduleEigengenes(exp, net$colors)$eigengenes
kME <- signedKME(exp, MEs, corFnc = "bicor")

net$new_colors <- net$colors
moduleColors <- labels2colors(net$colors)
table(moduleColors)
pdf('/methylation/methy_AD_revision2/results/WGCNA/dendrogram_RNA.pdf',width=9,height=5)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Modulecolors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

pro_num_module <- data.frame(table(moduleColors))

color_table <- data.frame(
  ModuleNumber = sort(unique(net$colors)),
  ModuleColor = labels2colors(sort(unique(net$colors)))
)
color_table$pro_num <- pro_num_module$Freq[match(color_table$ModuleColor, pro_num_module$moduleColors)]
color_table$module <- paste0('M',color_table$ModuleNumber,color_table$ModuleColor)
color_table <- color_table[-1,]
color_table$module <- paste0('M',color_table$ModuleNumber)
color_table$module <- factor(color_table$module,levels=paste0('M',1:17))

pdf('/methylation/methy_AD_revision2/results/WGCNA/count_per_module_RNA.pdf',width=9,height=18)
p <- ggplot(color_table, aes(y = pro_num, x = module, fill = ModuleColor)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = pro_num, y = pro_num + 20), 
            # angle = 45,
            # hjust = 0.5,
            size = 3.5) +
  scale_fill_identity() + 
  scale_y_continuous(expand=c(0,0),limits = c(0,580)) +
  labs(
    x = NULL,
    y = 'Count',
    title = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 13),
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = 'none'
  )
print(p)
dev.off()
save(sft,net,power_selected, moduleColors, color_table, MEs, kME, file='/methylation/methy_AD_revision2/results/WGCNA/res_RNA.RData')

load('/methylation/methy_AD_revision2/results/WGCNA/sft_RNA.RData')
sft_data <- data.frame(Power = sft$fitIndices[, 1], R2 = -sign(sft$fitIndices[,3])*sft$fitIndices[,2], mean_connectivity = sft$fitIndices[, 5])
power_selected <- sft_data$Power[which(sft_data$R2 > 0.8 & sft_data$mean_connectivity < 500)[1]]
load('/methylation/methy_AD_revision2/results/WGCNA/res_RNA.RData')
adnc_lim_all <- list(
  net = net,
  power =  power_selected,
  moduleColors = moduleColors,
  color_table = color_table,
  MEs = MEs,
  kME = kME
)

sample_infor = merge(pheno_df[,c("Sample_Name","Sample_Plate","bank","sex_male","age","PMD","RIN","AD","Exc","Inh","NonN_Astro_FGF3R","NonN_Endo","NonN_Micro","NonN_Oligo_MBP")], selected_sample[,c("id","Braak.NFT.stage","A_beta_0_3","C_cerad_0_3")], by.x='Sample_Name', by.y='id')
lm_res = data.frame('me'=NA, 'trait'=NA, 'B'=NA, 'P.Value'=NA)
module_trait_cor <- function(data){
  me <- data$MEs[,2:ncol(data$MEs)]

  trait = data.frame(matrix(ncol=ncol(sample_infor),nrow = nrow(me)))
  for (i in 1:ncol(sample_infor)) {
    trait[,i] = sample_infor[,i][match(rownames(me),sample_infor$Sample_Name)]
  }
  colnames(trait) = colnames(sample_infor)
  trait$Sample_Plate = factor(trait$Sample_Plate)
  trait$bank = factor(trait$bank)

  rownames(trait) = trait$Sample_Name
  identical(rownames(me),rownames(trait))
  for (pheno in c('age','AD','Braak.NFT.stage','A_beta_0_3','C_cerad_0_3')) {
    for (i in 1:ncol(me)) {
      if (pheno %chin% c('AD','Braak.NFT.stage','A_beta_0_3','C_cerad_0_3')){
        fit <- lm(me[,i] ~ trait[,pheno] + sex_male + age + bank + PMD + RIN + Exc + Inh + NonN_Astro_FGF3R + NonN_Endo + NonN_Micro + NonN_Oligo_MBP + Sample_Plate, data=trait)
      } else {
        fit <- lm(me[,i] ~ trait[,pheno] + sex_male + bank + PMD + RIN + Exc + Inh + NonN_Astro_FGF3R + NonN_Endo + NonN_Micro + NonN_Oligo_MBP + Sample_Plate, data=trait)
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
                                ifelse(module_trait_all$trait == 'age','Age',
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
openxlsx::write.xlsx(module_trait_all,'/methylation/methy_AD_revision2/results/WGCNA/module_trait_association_RNA.xlsx')

#### enrichment #####
setwd('/methylation/methy_AD_revision2/results/WGCNA')
datKME=signedKME(exp, adnc_lim_all$MEs, outputColumnName="kME_MM.")
hubcpg = list()
modules_of_interest <- c(1,5,6,15,18)
for (mod in modules_of_interest) {
  mod_col <- paste0("kME_MM.", mod)
  top_cpg <- rownames(datKME)[order(-abs(datKME[[mod_col]]))][1:50]
  hubcpg[[paste0("ME", mod)]] <- top_cpg
}
save(hubcpg,file='hub_cpg_RNA.RData')

hubcpg = list()
for (i in modules_of_interest) {
  hubcpg[[paste0("ME", i)]] = names(net$colors)[which(net$colors == i)]
}

exp = read.table('/CBMAP_RNAseq_protein_coding/RNAseq_process/RNAseq_data_qn_comb.txt',check.names = F)
bg_gene = rownames(exp);rm(exp)
bg_gene.df = bitr(bg_gene,fromType="ENSEMBL",toType="ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID

ego_ALL = list()
for (i in 1:length(hubcpg)) {
  gene = hubcpg[[i]]
  if (length(gene) == 0) {next} else {
    a <- enrichGO(
      gene = gene,
      OrgDb = org.Hs.eg.db,
      keyType = "ENSEMBL",
      ont = "ALL",
      universe = bg_gene,
      pAdjustMethod = "BH",
      minGSSize = 5,
      maxGSSize = 5000,
      pvalueCutoff = 1,
      qvalueCutoff = 0.05,
      readable = TRUE
    )
    if (nrow(a@result)==0) {next} else {
      ego_ALL[[i]] <- a
      ego_ALL[[i]] <- as.data.frame(ego_ALL[[i]])
      ego_ALL[[i]] <- ego_ALL[[i]][order(ego_ALL[[i]]$qvalue), ]
      ego_ALL[[i]]$module <- names(hubcpg)[i]
    }
  }
}
ego_res_ALL <- rbindlist(ego_ALL)
write.xlsx(ego_res_ALL, file='GOenrich_adnc_module_RNA.xlsx',colNames=T)
