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
all.cpg = meta_final_ordered_df$cpg[which(! is.na(meta_final_ordered_df$pVal.final))] 
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
Mdat <- h5read("final_DNAm_Mdat.h5",'Mdat')
colnames(Mdat) = h5read("final_DNAm_Mdat.h5",'colnames')
rownames(Mdat) = h5read("final_DNAm_Mdat.h5",'rownames')
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
cg = cg[,cpg_ids]


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
  ))
}


adnc_lim_all <- wgcna(cg)
save(adnc_lim_all, file = "meta_WGCNA_res.RData")

sample_infor = merge(pheno_df[,c(2,5:8,10:11,14:21)], selected_sample[,c(1,33,49,51)], by.x='Sample_Name', by.y='id')
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
module_trait_all$trait <- ifelse(module_trait_all$trait=='age','Age',
                                 ifelse(module_trait_all$trait == 'adnc','ADNC',
                                               ifelse(module_trait_all$trait == 'braak','Braak.NFT.stage',
                                                      ifelse(module_trait_all$trait == 'abeta','A-beta',
                                                             ifelse(module_trait_all$trait == 'cscore','C-score',NA)))))
module_trait_all$pmarker <- ifelse(module_trait_all$adj.P.Val < 0.05, "*", '')
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
openxlsx::write.xlsx(module_trait_all,'/methylation/methy_AD_revision2/results/WGCNA/module_trait_association.xlsx')

#### enrichment #####
setwd('/methylation/methy_AD_revision2/results/WGCNA')
datKME=signedKME(cg, adnc_lim_all$MEs, outputColumnName="kME_MM.")
hubcpg = list()
modules_of_interest <- c(7,8,31,37,39,40)
for (mod in modules_of_interest) {
  mod_col <- paste0("kME_MM.", mod)
  top_cpg <- rownames(datKME)[order(-abs(datKME[[mod_col]]))][1:50]
  hubcpg[[paste0("ME", mod)]] <- top_cpg
}

hubcpg = list()
for (i in c(7,8,31,37,39,40)) {
  hubcpg[[paste0("ME", i)]] = names(net$colors)[which(net$colors == i)]
}
cpg_gene_res = read.table('/eQTM/cis_eQTM_all_res_clr.txt',header=T)
cpg_gene_res <- cpg_gene_res %>%
  group_by(cpg) %>%
  mutate(fdr = p.adjust(pvalue, method = "BH")) %>%
  ungroup()
cpg_gene_res = cpg_gene_res[which(cpg_gene_res$fdr < 0.05),]

exp = read.table('/CBMAP_RNAseq_protein_coding/RNAseq_process/RNAseq_data_qn_comb.txt',check.names = F)
bg_gene = rownames(exp);rm(exp)
bg_gene.df = bitr(bg_gene,fromType="ENSEMBL",toType="ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID

ego_ALL = list()
p = list()
for (i in 1:length(hubcpg)) {
  cpg = hubcpg[[i]]
  gene = unlist(cpg_gene_res[which(cpg_gene_res$cpg %chin% cpg),'gene_id'])
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
      
      # if (dim(ego_ALL[[i]])[1] > 15) {
      #   ego_ALL[[i]] <- ego_ALL[[i]][1:15, ]
      # }
      # ego_ALL[[i]]$GeneRatio <- round(as.integer(str_split_fixed(ego_ALL[[i]]$GeneRatio,'/',n=2)[,1]) / as.integer(str_split_fixed(ego_ALL[[i]]$GeneRatio,'/',n=2)[,2]), 4)
      # ego_ALL[[i]]$type_order <- factor(rev(1:nrow(ego_ALL[[i]])),labels=rev(ego_ALL[[i]]$Description))
      # ego_ALL[[i]]$type_order <- factor(ego_ALL[[i]]$type_order,
      #                                   levels = ego_ALL[[i]]$type_order[order(ego_ALL[[i]]$qvalue, decreasing = T)])
    }
    # p[[i]] = ggplot(ego_ALL[[i]],aes(x=-log10(qvalue),y=type_order,size=GeneRatio,color=qvalue))+
    #   geom_point()+
    #   scale_size(range=c(2, 8))+
    #   scale_color_gradient('P.adjust',low="#698b69",high = "#c1ffc1")+
    #   theme_bw(base_size=14)+
    #   labs(color=expression(qvalue,size="GeneRatio"),
    #        x="-Log10(P.adjust)",y="",title=names(hubcpg)[i])+
    #   #theme(legend.title=element_text(size=14), legend.text = element_text(size=14))+
    #   theme(axis.title.x = element_text(margin = margin(t = 10)),
    #         axis.text = element_text(color="black",size=13))
  }
}
# p[[1]] + p[[2]] + p[[3]] + p[[4]] + p[[5]] + p[[6]] + plot_layout(ncol = 2)
ego_res_ALL <- rbindlist(ego_ALL)
write.xlsx(ego_res_ALL, file='GOenrich_adnc_module.xlsx',colNames=T)

