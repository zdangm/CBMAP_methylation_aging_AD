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

# Plot GA trajectories for CpGs
pdf('meta_hml_loess.pdf')
p1 = ggplot(cg.tf, aes(x = age, y = BetaValue, group = cpg)) +
  geom_smooth(method = "loess", method.args = list(span = 0.75, degree = 1), se = FALSE, linewidth = .2, color = "grey40") +
  theme_classic(base_size = 16) +
  labs(y="Methylation level",x="Chronological age (years)")+
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20))
print(p1)
dev.off()

#  cluster  
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

# heatmap plot
heatmap.df <- res_summary
heatmap.df$ID <- cg$ID
heatmap.df$age <- pheno_df[match(heatmap.df$ID, pheno_df$Sample_Name),'age']
heatmap.df <- heatmap.df %>% arrange(age)
heatmap.df.tf <- gather(heatmap.df, 1:(ncol(heatmap.df)-2), key = "cpg", value = "BetaValue")
heatmap.df.tf <- heatmap.df.tf[which(heatmap.df.tf$age>=14),]

all_ages <- seq(14, max(heatmap.df.tf$age), by = 1)
data_filled <- heatmap.df.tf %>%
  complete(age = all_ages, cpg, fill = list(Zscore = NA)) %>% 
  group_by(cpg) %>%
  mutate(BetaValue = zoo::na.approx(BetaValue, age, na.rm = FALSE))
cpg.order = data_filled[which(data_filled$age==min(data_filled$age)),]
cpg.order = unlist(cpg.order[order(cpg.order$BetaValue,decreasing=T),'cpg'])
data_filled$cpg = factor(data_filled$cpg,levels=cpg.order)

ggplot(data_filled, aes(x = age, y = cpg, fill = BetaValue)) +
  geom_tile() + 
  scale_fill_gradient2(low = "navy", mid = "white", high = "red4", midpoint = 0) +
  #stat_density_2d(aes(alpha=..level..), geom = "polygon", contour = TRUE) +
  theme_minimal(base_size = 13) +
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        panel.grid = element_blank()) +
  labs(x = "Chronological age", y = paste0(length(unique(cg.tf$cpg))," CpGs"), fill = "Z-score")

# Perform kmeans clustering
set.seed(42) 
kmeans_result <- kmeans(loess.mat, centers = 4)
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
for (cluster in 1:4) {
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
pdf('meta_hml_loess.pdf',width=8,height=2)
p4 = cluster_gg[[1]]+cluster_gg[[2]]+cluster_gg[[3]]+cluster_gg[[4]]+plot_layout(ncol=4)
print(p4)
dev.off()

# Enrichment analysis 
cpg_gene_res = read.table('gene_cpg_res2MB.txt',header=T)
cpg_gene_res = cpg_gene_res[which(cpg_gene_res$fdr<0.05),]

exp = read.table('RNAseq_data_qn_comb.txt',check.names = F)
load(gencode_v46_b38_annotation.RData')
gene_pos = gencode_v46_annotation[which(gencode_v46_annotation$type == 'gene'),c('gene_id','chr','start','end','gene_name')]
gene_pos$gene_id = unlist(str_split(gene_pos$gene_id,"[.]",simplify=T))[,1]
gene_pos = gene_pos[which(gene_pos$gene_id %chin% rownames(exp)),]
# delete_gene = gene_pos$gene_id[which(gene_pos$chr=='chr6' & gene_pos$start>=30500001 & gene_pos$end<=46200000)]
# cpg_gene_res = cpg_gene_res[which(! cpg_gene_res$gene %chin% delete_gene),] # 2619 rows  # 3727 rows

cluster_cpg = list()
for(i in 1:4) {
  cluster_cpg[[i]] = list()
  cluster_cpg[[i]][['cpg']] = unique(cg_df.cluster[which(cg_df.cluster$Cluster==i),'cpg'])
  cluster_cpg[[i]][['ENSEMBL']] = unique(cpg_gene_res[which(cpg_gene_res$cpg %chin% cluster_cpg[[i]][['cpg']]),'gene'])
  cluster_cpg[[i]][['ENTREZID']] = bitr(cluster_cpg[[i]][['ENSEMBL']],fromType="ENSEMBL",toType="ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
}
bg_gene = rownames(exp) 
bg_gene.df = bitr(bg_gene,fromType="ENSEMBL",toType="ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID

# GO enrichment analysis
ego_ALL <- list()
for (j in 1:4) {
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
  else if (nrow(ego_ALL[[j]])>0 & nrow(ego_ALL[[j]])<15) {
    ego_ALL[[j]] <- ego_ALL[[j]][order(ego_ALL[[j]]$qvalue),]
    ego_ALL[[j]]$cluster <- paste0('Cluster ',j)
  } else {
    ego_ALL[[j]] <- ego_ALL[[j]][order(ego_ALL[[j]]$qvalue),]
    ego_ALL[[j]] <- ego_ALL[[j]][1:15,]
    ego_ALL[[j]]$cluster <- paste0('Cluster ',j)
  }  
}
for (j in 1:4) {
  if (j == 1) {
    ego_res = ego_ALL[[j]]
  } else {
    ego_res = rbind(ego_res, ego_ALL[[j]])
  }
}
ego_res$type_order = factor(rev(1:nrow(ego_res)),labels=rev(ego_res$Description))
ego_res$GeneRatio = round(as.integer(str_split_fixed(ego_res$GeneRatio,'/',n=2)[,1]) / as.integer(str_split_fixed(ego_res$GeneRatio,'/',n=2)[,2]), 2)

ggplot(ego_res,aes(y=type_order,x=cluster,group=cluster))+
  geom_point(aes(size=GeneRatio,color=qvalue)) +
  scale_color_gradient('P.adjust',limits=c(0,0.05), breaks=c(0.01,0.02,0.03,0.04), low="#CDE8C3",high = "#375637")+
  labs(color=expression(qvalue,size="GeneRatio"), 
       x=NULL,y=NULL)+
  theme_bw()+
  theme(axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=11, vjust=0.6, hjust=.7, angle=30))

# KEGG enrichment analysis
kegg_ALL = list()
for (j in 1:4) {
  kegg_ALL[[j]] <- enrichKEGG(gene = cluster_cpg[[j]][['ENTREZID']],keyType = "kegg",organism= "human", minGSSize = 5, maxGSSize = 5000,universe=bg_gene.df, qvalueCutoff = 0.05, pvalueCutoff=1)
  kegg_ALL[[j]] <- as.data.frame(kegg_ALL[[j]])
  if (nrow(kegg_ALL[[j]])==0) {next} else {
    rownames(kegg_ALL[[j]]) <- 1:nrow(kegg_ALL[[j]])
    kegg_ALL[[j]] <- kegg_ALL[[j]][order(kegg_ALL[[j]]$qvalue), ]
    if (nrow(kegg_ALL[[j]])>20) {
      kegg_ALL[[j]] <- kegg_ALL[[j]][1:20, ]
    }
    # kegg_ALL[[j]]$type_order <- factor(rev(1:nrow(kegg_ALL[[j]])),labels=rev(kegg_ALL[[j]]$Description))
    # kegg_ALL[[j]]$GeneRatio <- round(as.integer(str_split_fixed(kegg_ALL[[j]]$GeneRatio,'/',n=2)[,1]) / as.integer(str_split_fixed(kegg_ALL[[j]]$GeneRatio,'/',n=2)[,2]), 2)
    for (i in 1:dim(kegg_ALL[[j]])[1]) {
      kegg_ALL[[j]]$geneID[i] <- str_replace_all(kegg_ALL[[j]]$geneID[i], "/", ",")
      gene <- as.integer(str_split(kegg_ALL[[j]]$geneID[i],',')[[1]])
      gene_name = bitr(gene ,fromType="ENTREZID",toType="SYMBOL", OrgDb = org.Hs.eg.db)$SYMBOL
      kegg_ALL[[j]]$geneID[i] <- paste(gene_name[1:10], collapse = ', ')
    }
    kegg_ALL[[j]]$set <- paste0('Cluster ', j)
  }  
}
for (j in 1:4) {
  if (j == 1) {
    kegg_res = kegg_ALL[[j]]
  } else {
    kegg_res = rbind(kegg_res, kegg_ALL[[j]])
  }
}
kegg_res$type_order = factor(rev(1:nrow(kegg_res)),labels=rev(kegg_res$Description))
kegg_res$GeneRatio = round(as.integer(str_split_fixed(kegg_res$GeneRatio,'/',n=2)[,1]) / as.integer(str_split_fixed(kegg_res$GeneRatio,'/',n=2)[,2]), 2)

ggplot(kegg_res,aes(y=type_order,x=set,group=set))+
  geom_point(aes(size=GeneRatio,color=qvalue)) +
  scale_color_gradient('P.adjust',limits=c(0,0.05), breaks=c(0.01,0.02,0.03,0.04), low="#CDE8C3",high = "#375637")+
  labs(color=expression(qvalue,size="GeneRatio"), 
       x="",y="")+
  theme_bw()

