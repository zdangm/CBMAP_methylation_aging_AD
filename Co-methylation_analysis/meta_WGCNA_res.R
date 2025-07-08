library(rhdf5)
library(readxl)
library(data.table)
library(impute)
library(WGCNA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(writexl)
library(dplyr)
library(tidyr)
library(forcats)
library(stringr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(gridExtra)
library(scales)
library(AnnotationDbi)
library(openxlsx)
library(patchwork)
library(RColorBrewer)

load('meta_hml4level.RData')
cpg = meta_final_ordered_df$cpg[which(meta_final_ordered_df$pVal.final<0.05)] 
load('CBMAP_pheno.RData')
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
cg <- h5read("final_DNAm_invMdat.h5",'Mdat')
colnames(cg) = h5read("final_DNAm_invMdat.h5",'colnames')
rownames(cg) = h5read("final_DNAm_invMdat.h5",'rownames')
cg = cg[,grep("cg",colnames(cg))] 
sample = intersect(rownames(cg), rownames(pheno_df)) 
cg <- cg[sample,cpg]

wgcna <- function(data){
  print(dim(data))
  
  gsg <- goodSamplesGenes(data)
  if (!gsg$allOK) {
    data <- data[gsg$goodSamples, gsg$goodGenes]
  }
  
  powers = c(1:20)
  sft <- pickSoftThreshold(data, 
                           networkType = "signed",
                           corFnc = "bicor",
                           powerVector = powers,
                           verbose = 5)
  pdf('meta_power_selected.pdf',width=9,height=5)
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
  
  sft_data <- data.frame(Power = sft$fitIndices[, 1], R2 = -sign(sft$fitIndices[,3])*sft$fitIndices[,2], mean_connectivity = sft$fitIndices[, 5])
  power_selected <- sft_data$Power[which(sft_data$R2 > 0.8 & sft_data$mean_connectivity < 500)[1]]
  print(power_selected)
  
  cor <- WGCNA::cor
  net <- blockwiseModules(
    data,
    power = power_selected,
    corType = "bicor",          
    networkType = "signed",     
    TOMType = "signed",         
    deepSplit = 4,              
    pamRespectsDendro = FALSE,   
    pamStage = FALSE,
    minModuleSize = 30,         
    mergeCutHeight = 0.15,      
    numericLabels = TRUE,       
    verbose = 3
  )
  
  table(net$colors)
  
  MEs <- moduleEigengenes(data, net$colors)$eigengenes
  kME <- signedKME(data, MEs, corFnc = "bicor")
  
  net$new_colors <- net$colors
  moduleColors <- labels2colors(net$colors)
  table(moduleColors)
  pdf('meta_dendrogram.pdf',width=9,height=5)
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
  
  pdf('meta_count_per_module.pdf',width=9,height=9)
  p <- ggplot(color_table, aes(x = pro_num, y = reorder(module, pro_num), fill = ModuleColor)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = pro_num), 
              hjust = -0.1,  
              size = 4) +    
    scale_fill_identity() +  
    labs(
      x = NULL,
      y = NULL,
      title = NULL
    ) +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 13),
      plot.title = element_text(size = 14, face = "bold")
    ) +
    xlim(0, max(color_table$pro_num) * 1.1)  
  print(p)
  dev.off()
  
  kWithin <- intramodularConnectivity(adjacency(data, power = power_selected), net$colors)
  
  hub_pro <- chooseTopHubInEachModule(
    data,
    colorh = moduleColors,
    power = power_selected,
    type = "signed"
  )
    
  return(list(
    net = net,
    power =  power_selected,
    moduleColors = moduleColors,
    color_table = color_table,
    MEs = MEs,
    kME = kME
    # ,hub_gene = hub_gene_module_info
  ))
}


adnc_lim_all <- wgcna(cg)
save(adnc_lim_all, file = "meta_WGCNA_res.RData")

sample_infor = merge(pheno_df[,c(2,5,7,8,12,14)], selected_sample[,c(1,33,49,51)], by.x='Sample_Name', by.y='id')
colnames(sample_infor) = c('id','Sample_Plate','sex_male','age','NeuN_pos','adnc_num','Braak.NFT.stage','A_beta_0_3','C_cerad_0_3')
lm_res = data.frame('me'=NA, 'trait'=NA, 'B'=NA, 'P.Value'=NA)
module_trait_cor <- function(data){
  me <- data$MEs[,2:ncol(data$MEs)]
  
  trait <- data.frame(
    id = sample_infor$id[match(rownames(me),sample_infor$id)],
    adnc = sample_infor$adnc_num[match(rownames(me),sample_infor$id)],
    age = sample_infor$age[match(rownames(me),sample_infor$id)],
    sex = sample_infor$sex_male[match(rownames(me),sample_infor$id)],
    braak = sample_infor$Braak.NFT.stage[match(rownames(me),sample_infor$id)],
    abeta = sample_infor$A_beta_0_3[match(rownames(me),sample_infor$id)],
    cscore = sample_infor$C_cerad_0_3[match(rownames(me),sample_infor$id)],
    Sample_Plate = sample_infor$Sample_Plate[match(rownames(me),sample_infor$id)],
    NeuN_pos = sample_infor$NeuN_pos[match(rownames(me),sample_infor$id)]
  )
  trait$Sample_Plate = factor(trait$Sample_Plate)
  # linear
  sample = intersect(rownames(me),trait$id)
  rownames(trait) = trait$id
  me = me[sample,]; trait = trait[sample,]
  identical(rownames(me),rownames(trait))
  for (pheno in c('age','adnc','braak','abeta','cscore')) {
    for (i in 1:ncol(me)) {
      if (pheno %chin% c('adnc','braak','abeta','cscore')){
        fit <- lm(me[,i] ~ trait[,pheno] + age + sex + Sample_Plate + NeuN_pos, data=trait)
      } else {
        fit <- lm(me[,i] ~ trait[,pheno], data=trait)
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
    mutate(adj.P.Val = p.adjust(P.Value, method = "fdr")) %>%
    ungroup()
  return(lm_res)
}
module_trait_all <- module_trait_cor(adnc_lim_all)
module_trait_all$trait <- ifelse(module_trait_all$trait=='age','Age',
                                 ifelse(module_trait_all$trait == 'adnc','ADNC',
                                               ifelse(module_trait_all$trait == 'braak','Braak.NFT.stage',
                                                      ifelse(module_trait_all$trait == 'abeta','A-beta',
                                                             ifelse(module_trait_all$trait == 'cscore','C-score',NA)))))
module_trait_all$pmarker <- ifelse(module_trait_all$adj.P.Val < 0.05, "*", '')
module_trait_all$me <- factor(module_trait_all$me, levels = paste0('ME', 1:39))
module_trait_all$trait <- factor(module_trait_all$trait, levels = c('Age','ADNC','Braak.NFT.stage','A-beta','C-score'))
