library(rhdf5)
library(data.table)
library(dplyr)
library (tidyr)
library(ggplot2)
library(dendextend)
library(RColorBrewer)
library(cowplot)
library(patchwork)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)
library(stringr)
library(openxlsx)

load('meta_hml4level.RData')
pheno_df = read.csv('phenotype.csv',header=T)
cpg = meta_final_ordered_df[which(meta_final_ordered_df$fdr<0.05),'cpg']
Mdat <- h5read("final_DNAm_Mdat.h5",'Mdat')
colnames(Mdat) = h5read("final_DNAm_Mdat.h5",'colnames')
rownames(Mdat) = h5read("final_DNAm_Mdat.h5",'rownames')
cg <- 2^Mdat / (2^Mdat + 1)
cg <- t(cg)
cg = cg[,cpg]
cg <- data.frame(scale(cg))
cg$ID <- rownames(cg)
cg$age <- pheno_df[match(cg$ID, pheno_df$Sample_Name),'age']
cg$sex <- pheno_df[match(cg$ID, pheno_df$Sample_Name),'sex_male']
cg$sex <- factor(cg$sex)
cg <- cg %>% arrange(age)

# Transform data into a long format for plotting
cg.tf <- gather(cg, 1:(ncol(cg)-3), key = "cpg", value = "BetaValue")
inv_logit <- function(x) {
  exp(x) / (1 + exp(x))
}
cg.tf$BetaValue <- inv_logit(cg.tf$BetaValue)

res_summary <- NULL
# Fit loess model for each protein and extract fitted values
for (i in 1:(ncol(cg)-3)) {
  loess_model <- loess(cg[, i, drop = TRUE] ~ age, data = cg)
  fitted_values <- data.frame(loess_model$fitted)
  colnames(fitted_values) <- colnames(cg)[i]
  if (is.null(res_summary)) {
    res_summary <- fitted_values
  } else {
    res_summary <- cbind(res_summary, fitted_values)
  }
}
loess.mat <- t(res_summary)
colnames(loess.mat) <- cg$ID

# Perform kmeans clustering
loess.mat.scaled <- t(scale(t(loess.mat)))
set.seed(42) 
k.max <- 10
wss <- numeric(k.max)
for (k in 1:k.max) {
  km <- kmeans(
    loess.mat.scaled,
    centers = k,
    nstart = 12,
    iter.max = 100
  )
  wss[k] <- km$tot.withinss
}
pdf('Elbow.pdf',height = 5.70, width = 5.31)
plot(
  1:k.max, wss,
  type = "b",
  pch = 19,
  frame = FALSE,
  xlab = "Number of clusters (k)",
  ylab = "Total within-cluster sum of squares"
)
dev.off()
kmeans_result <- kmeans(loess.mat.scaled, centers = 2)
col6 <- brewer.pal(12, "Paired")[c(1,3,5,7,9,11)]
col62 <- brewer.pal(12, "Paired")[c(2,4,6,8,10,12)]
clustering.res.df <- data.frame(cpg = names(kmeans_result$cluster), Cluster = kmeans_result$cluster)
cg_df.cluster <- left_join(cg.tf, clustering.res.df, by = c("cpg" = "cpg"))
# Count cpgs in each cluster
cluster.count <- cg_df.cluster %>% dplyr::count(Cluster)
cluster.count$n <- cluster.count$n/1031
cg_df.cluster$Cluster <- as.factor(cg_df.cluster$Cluster)

# Prepare plots for each cluster
cluster_gg <- list()
for (cluster in 1:2) {
  cluster_data <- filter(cg_df.cluster, Cluster == cluster)
  cluster_gg[[cluster]] <- ggplot(cluster_data, aes(x = age, y = BetaValue)) +
    theme_classic(base_size = 9) +
    theme(aspect.ratio = 1,  legend.position = "none",
          plot.title = element_text(size = 8)) +
    geom_smooth(method = "loess", method.args = list(span = .75, degree = 1),aes(group = cpg), se = FALSE, color = col6[cluster], linewidth = 0.01) + 
    geom_smooth(method = "loess", method.args = list(span = .75, degree = 1),se = FALSE, color = col62[cluster], linewidth = 0.5) +
    ylab("Methylation level") + xlab("Chronological age") + 
    #coord_cartesian(ylim = c(-2.51, 2.51)) +
    scale_y_continuous(limits=c(0,1)) +
    ggtitle(paste0("Cluster ", cluster, " (", cluster.count$n[cluster], " CpGs)"))
}

# Combine and plot all cluster plots in a grid
pdf('meta_hml_loess.pdf',width=5.4,height=2.7)
p4 = cluster_gg[[1]]+cluster_gg[[2]]
print(p4)
dev.off()

# Enrichment analysis 
cpg_gene_res = fread('cis_eQTM_all_res.txt',header=T)
cpg_gene_res <- cpg_gene_res %>%
  group_by(cpg) %>%
  mutate(fdr = p.adjust(pvalue, method = "fdr")) %>%
  ungroup()
cpg_gene_res = cpg_gene_res[which(cpg_gene_res$fdr < 0.05),]
exp = read.table('RNAseq_data_qn_comb.txt',check.names = F)
bg_gene = rownames(exp);rm(exp)
bg_gene.df = bitr(bg_gene,fromType="ENSEMBL",toType="ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
cluster_cpg = list()
for(i in 1:2) {
  cluster_cpg[[i]] = list()
  cluster_cpg[[i]][['cpg']] = unique(clustering.res.df[which(clustering.res.df$Cluster==i),'cpg'])
  cluster_cpg[[i]][['ENSEMBL']] = unique(cpg_gene_res[which(cpg_gene_res$cpg %chin% cluster_cpg[[i]][['cpg']]),'gene_id'])
  cluster_cpg[[i]][['ENTREZID']] = bitr(cluster_cpg[[i]][['ENSEMBL']]$gene_id,fromType="ENSEMBL",toType="ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
}

# GO enrichment analysis
ego_ALL <- list()
for (j in 1:2) {
  ego_ALL[[j]] <- enrichGO(
    gene = cluster_cpg[[j]][['ENSEMBL']],
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
  ego_ALL[[j]] <- as.data.frame(ego_ALL[[j]])
  if (nrow(ego_ALL[[j]])==0) {next} 
  else {
    ego_ALL[[j]] <- ego_ALL[[j]][order(ego_ALL[[j]]$qvalue),]
    ego_ALL[[j]]$cluster <- paste0('Cluster ',j)
  }  
}
for (j in 1:2) {
  if (j == 1) {
    ego_res = ego_ALL[[j]]
  } else {
    ego_res = rbind(ego_res, ego_ALL[[j]])
  }
}
write.xlsx(ego_res, file='/enrichment/hml_loess_GOenrich.xlsx',colNames=T)

# KEGG enrichment analysis
kegg_ALL = list()
for (j in 1:2) {
  kegg_ALL[[j]] <- enrichKEGG(gene = cluster_cpg[[j]][['ENTREZID']],keyType = "kegg",organism= "human", minGSSize = 5, maxGSSize = 5000,universe=bg_gene.df, qvalueCutoff = 0.05, pvalueCutoff=1)
  kegg_ALL[[j]] <- as.data.frame(kegg_ALL[[j]])
  if (nrow(kegg_ALL[[j]])==0) {next} else {
    rownames(kegg_ALL[[j]]) <- 1:nrow(kegg_ALL[[j]])
    kegg_ALL[[j]] <- kegg_ALL[[j]][order(kegg_ALL[[j]]$qvalue), ]
    kegg_ALL[[j]]$set <- paste0('Cluster ', j)
  }  
}
for (j in 1:2) {
  if (j == 1) {
    kegg_res = kegg_ALL[[j]]
  } else {
    kegg_res = rbind(kegg_res, kegg_ALL[[j]])
  }
}
write.xlsx(kegg_res, file='/enrichment/hml_loess_KEGGenrich.xlsx',colNames=T)
