library(rhdf5)
library(data.table)
library(dplyr)
library(Rfast)
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
library(matrixStats)
library(proxy)
library(writexl) 


load('CBMAP_pheno.RData')
Mdat <- h5read("final_DNAm_Mdat.h5",'Mdat')
colnames(Mdat) = h5read("final_DNAm_Mdat.h5",'colnames')
rownames(Mdat) = h5read("final_DNAm_Mdat.h5",'rownames')
cg <- 2^Mdat / (2^Mdat + 1)
cg <- t(cg)
sample = intersect(rownames(cg), rownames(pheno_df)) 
cg <- cg[sample,]
pheno_df <- pheno_df[sample,]

# Calculate the mean methylation level for group age
unique_ages <- sort(unique(pheno_df$age))
result_mat <- matrix(NA, nrow = length(unique_ages), ncol = ncol(cg))
rownames(result_mat) <- unique_ages
colnames(result_mat) <- colnames(cg)
for (i in seq_along(unique_ages)) {
  age_i <- unique_ages[i]
  idx <- which(pheno_df$age == age_i)
  result_mat[i, ] <- colMeans(cg[idx, , drop = FALSE], na.rm = TRUE)
}
var <- colVars(result_mat, na.rm = TRUE)

# Linear regression by cpgs Methylation
age = as.numeric(rownames(result_mat))
run_linear_regression <- function(methylation_data) {
  data <- data.frame(methylation_level = methylation_data)
  linear_model <- lm(methylation_level ~ age, data=data)
  beta_lm <- summary(linear_model)$coefficients[2,1]
  p_lm <- summary(linear_model)$coefficients[2,4]
  loess_model <- loess(methylation_level ~ age, data = data)
  y <- data$methylation_level
  y_pred_lm <- predict(linear_model)
  y_pred_loess <- predict(loess_model)

  rss_lm <- sum((y - y_pred_lm)^2)
  rss_loess <- sum((y - y_pred_loess)^2)
  tss <- sum((y - mean(y))^2)

  r2_lm <- 1 - rss_lm / tss
  r2_loess <- 1 - rss_loess / tss
  res <- data.frame('beta_lm'=beta_lm, 'p_lm'=p_lm, 'r2_lm'=r2_lm, 'r2_loess'=r2_loess)
  return(res)
}
analyze_methylation <- function(cg) {
  results <- apply(cg, 2, function(x) {
    methylation_data <- x
    result <- run_linear_regression(methylation_data)
    return(result)
  })
  return(results)
}
results = analyze_methylation(result_mat)
results = rbindlist(results)
results$cpg = colnames(cg)
results$fdr = p.adjust(results$p_lm,'BH')
results$var = var[match(results$cpg, names(var))]
save(results, file = 'r2.RData')

# Nonlinear
load('r2.RData')
results = results[which(results$fdr < 0.05),]
cpg_neg = results[results$r2_lm > 0.7 & results$beta_lm < 0,] 
cpg_pos = results[results$r2_lm > 0.7 & results$beta_lm > 0,] 
cpg_neg = cpg_neg[order(cpg_neg$p_lm)[1:500],'cpg']
cpg_pos = cpg_pos[order(cpg_pos$p_lm)[1:500],'cpg']
delete_cpg = c(cpg_neg$cpg, cpg_pos$cpg)
results = results[which(! results$cpg %chin% delete_cpg),]
cpg = results[which(results$r2_loess-results$r2_lm>0.4),]  

load('CBMAP_pheno.RData')
# methylation matrix
Mdat <- h5read("final_DNAm_Mdat.h5",'Mdat')
colnames(Mdat) = h5read("final_DNAm_Mdat.h5",'colnames')
rownames(Mdat) = h5read("final_DNAm_Mdat.h5",'rownames')
cg <- 2^Mdat / (2^Mdat + 1)
cg <- t(cg)
sample = intersect(rownames(cg), rownames(pheno_df)) 
cg <- cg[sample,]
pheno_df <- pheno_df[sample,]

# Calculate the mean methylation level for group age
unique_ages <- sort(unique(pheno_df$age))
result_mat <- matrix(NA, nrow = length(unique_ages), ncol = ncol(cg))
rownames(result_mat) <- unique_ages
colnames(result_mat) <- colnames(cg)
for (i in seq_along(unique_ages)) {
  age_i <- unique_ages[i]
  idx <- which(pheno_df$age == age_i)
  result_mat[i, ] <- colMeans(cg[idx, , drop = FALSE], na.rm = TRUE)
}
age = as.numeric(rownames(result_mat))

age_breaks <- c(seq(0, 100, by = 20), 106)
age_labels <- c("0-20", paste(seq(21, 81, by = 20), seq(40, 100, by = 20), sep = "-"), "101-106")
age_group <- cut(age, breaks = age_breaks, labels = age_labels, include.lowest = TRUE, right = TRUE)
dt <- data.table(age_group = age_group, idx = 1:length(age))
group_var <- lapply(split(dt$idx, dt$age_group), function(rows) {
  colVars(result_mat[rows, , drop = FALSE], na.rm = TRUE)
})
group_var_mat <- do.call(rbind, group_var)
colnames(group_var_mat) = colnames(result_mat)
group_var_mat <- group_var_mat[,which(! colnames(group_var_mat) %chin% delete_cpg)]
#group_var_mat <- group_var_mat[-c(1,11),]
selected_cols <- apply(group_var_mat, 2, function(x) any(x > 0.05))
sum(selected_cols==T) 
filtered_matrix <- group_var_mat[, selected_cols]

cg <- cg[,colnames(filtered_matrix)]
cg <- cg[,cpg$cpg]
cg <- data.frame(cg)
cg <- data.frame(scale(cg))

cg$ID <- rownames(cg)
cg$age <- pheno_df[match(cg$ID, pheno_df$Sample_Name),'age']
cg$sex <- pheno_df[match(cg$ID, pheno_df$Sample_Name),'sex_male']
cg$sex <- factor(cg$sex)
cg <- cg %>% arrange(age)
cg.tf <- gather(cg, 1:(ncol(cg)-3), key = "cpg", value = "BetaValue")

#cluster  
res_summary <- NULL
# Fit loess model for each CpG and extract fitted values
for (i in 1:(ncol(cg)-3)) {
  loess_model <- loess(cg[, i, drop = TRUE] ~ age, data = cg)
  fitted_values <- data.frame(loess_model$fitted)
  colnames(fitted_values) <- colnames(cg)[i]
  if (is.null(res_summary)) {
    res_summary <- fitted_values
  } else {
    res_summary <- cbind(res_summary, fitted_values)
  }
  if (i %% 10000 == 0) {
    cat(i, " Done\n")
  }
}
loess.mat <- t(res_summary)
colnames(loess.mat) <- cg$ID

set.seed(42) 
kmeans_result <- kmeans(loess.mat, centers = 4)
paired_colors <- list(
  c("#A6CEE3", "#1F78B4"),
  c("#B2DF8A", "#33A02C"),
  c("#FB9A99", "#E31A1C"),
  c("#FDBF6F", "#FF7F00"),
  c("#CAB2D6", "#6A3D9A"),
  c("#FFFF99", "#B15928"),
  c("#D9D9D9", "#525252"),
  c("#A6DBD0", "#01665E")
)
paired_vector <- unlist(paired_colors)
col6 <- paired_vector[c(1,3,5,7,9,11,13,15)]
col62 <- paired_vector[c(2,4,6,8,10,12,14,16)]

# Join clustering results with cpg data
clustering.res.df <- data.frame(cpg = names(kmeans_result$cluster), Cluster = kmeans_result$cluster)
cg_df.cluster <- left_join(cg.tf, clustering.res.df, by = c("cpg" = "cpg"))
save(clustering.res.df, cg_df.cluster, file='loess_kmeans4_cluster_res.RData')

# Count cpgs in each cluster
load('loess_kmeans4_cluster_res.RData')
cluster.count <- cg_df.cluster %>% dplyr::count(Cluster)
cluster.count$n <- cluster.count$n/1031

# Convert Cluster to factor for plotting
cg_df.cluster$Cluster <- as.factor(cg_df.cluster$Cluster)

sampled_cpg = list()
for (i in 1:4) {
  sampled_cpg[[i]] = unique(cg_df.cluster$cpg[which(cg_df.cluster$Cluster == i)])
  if (length(sampled_cpg[[i]]) > 200) {
    set.seed(1)
    sampled_cpg[[i]] = sample(sampled_cpg[[i]], 200, replace=F)
  }
}
sampled_cpg = unlist(sampled_cpg)
cg_df.cluster.sample = cg_df.cluster[which(cg_df.cluster$cpg %chin% sampled_cpg),]

inv_logit <- function(x) {
  exp(x) / (1 + exp(x))
}
cg_df.cluster.sample$BetaValue <- inv_logit(cg_df.cluster.sample$BetaValue)

# Prepare plots for each cluster
cluster_gg <- list()
for (cluster in 1:4) {
  cluster_data <- filter(cg_df.cluster.sample, Cluster == cluster)
  cluster_gg[[cluster]] <- ggplot(cluster_data, aes(x = age, y = BetaValue)) +
    theme_classic(base_size = 20) +
    theme(aspect.ratio = 1,  legend.position = "none",
          plot.title = element_text(size = 18)) +
    geom_smooth(method = "loess", method.args = list(span = .75, degree = 1),aes(group = cpg), se = FALSE, color = col6[cluster], linewidth = 0.01) +
    geom_smooth(method = "loess", method.args = list(span = .75, degree = 1),se = FALSE, color = col62[cluster], linewidth = 0.5) +
    ylab("Methylation level") + xlab("Chronological age") +
    ggtitle(paste0("Cluster ", cluster, " (", cluster.count$n[cluster], " CpGs)"))
}
# Combine and plot all cluster plots in a grid
pdf('loess_kmeans4_cluster.pdf',width=10,height=10)
p = (cluster_gg[[1]]|cluster_gg[[2]])/(cluster_gg[[3]]|cluster_gg[[4]])
print(p)
dev.off()



####### linear ###########
load('r2.RData')
results = results[which(results$fdr < 0.05),]
cpg_neg = results[results$r2_lm > 0.7 & results$beta_lm < 0,] 
cpg_pos = results[results$r2_lm > 0.7 & results$beta_lm > 0,] 
cpg_neg = cpg_neg[order(cpg_neg$p_lm)[1:500],'cpg']
cpg_pos = cpg_pos[order(cpg_pos$p_lm)[1:500],'cpg']

pheno_df = read.csv('phenotype.csv',header=T)
Mdat <- h5read("final_DNAm_Mdat.h5",'Mdat')
colnames(Mdat) = h5read("final_DNAm_Mdat.h5",'colnames')
rownames(Mdat) = h5read("final_DNAm_Mdat.h5",'rownames')
cg <- 2^Mdat / (2^Mdat + 1)
cg <- t(cg)
cg_pos <- cg[,cpg_pos$cpg]
cg_pos <- data.frame(cg_pos)
cg_pos <- data.frame(scale(cg_pos))
cg_neg <- cg[,cpg_neg$cpg]
cg_neg <- data.frame(cg_neg)

cg_pos$ID <- rownames(cg_pos)
cg_pos$age <- pheno_df[match(cg_pos$ID, pheno_df$Sample_Name),'age']
cg_pos$sex <- pheno_df[match(cg_pos$ID, pheno_df$Sample_Name),'sex_male']
cg_pos$sex <- factor(cg_pos$sex)
cg_pos <- cg_pos %>% arrange(age)
cg_neg$ID <- rownames(cg_neg)
cg_neg$age <- pheno_df[match(cg_neg$ID, pheno_df$Sample_Name),'age']
cg_neg$sex <- pheno_df[match(cg_neg$ID, pheno_df$Sample_Name),'sex_male']
cg_neg$sex <- factor(cg_neg$sex)
cg_neg <- cg_neg %>% arrange(age)

# Transform data into a long format for plotting
cg.tf_pos <- gather(cg_pos, 1:(ncol(cg_pos)-3), key = "cpg", value = "BetaValue")
cg.tf_neg <- gather(cg_neg, 1:(ncol(cg_neg)-3), key = "cpg", value = "BetaValue")

inv_logit <- function(x) {
  exp(x) / (1 + exp(x))
}
cg.tf_pos$BetaValue <- inv_logit(cg.tf_pos$BetaValue)
cg.tf_neg$BetaValue <- inv_logit(cg.tf_neg$BetaValue)


pos_plot <- ggplot(cg.tf_pos, aes(x = age, y = BetaValue)) +
  theme_classic(base_size = 20) +
  theme(aspect.ratio = 1,  legend.position = "none",
        plot.title = element_text(size = 18)) +
  geom_smooth(method = "lm", aes(group = cpg), se = FALSE, color = "#FB9A99", linewidth = 0.05) +
  geom_smooth(method = "lm", se = FALSE, color = "#E31A1C", linewidth = 0.5) +
  ylab("Methylation level") + xlab("Chronological age") +
  ggtitle("Age-up CpGs")
neg_plot <- ggplot(cg.tf_neg, aes(x = age, y = BetaValue)) +
  theme_classic(base_size = 20) +
  theme(aspect.ratio = 1,  legend.position = "none",
        plot.title = element_text(size = 18)) +
  geom_smooth(method = "lm", aes(group = cpg), se = FALSE, color = "#A6CEE3", linewidth = 0.05) +
  geom_smooth(method = "lm", se = FALSE, color = "#1F78B4", linewidth = 0.5) +
  ylab("Methylation level") + xlab("Chronological age") +
  ggtitle("Age-down CpGs")

# Combine and plot all cluster plots in a grid
pdf('linear_cpg_plot.pdf',width=10,height=5)
p4 = pos_plot|neg_plot
print(p4)
dev.off()