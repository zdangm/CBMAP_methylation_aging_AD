library(rhdf5)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(stringr)
library(Rfast)
library(writexl)

####calculate methylation residual####
load('CBMAP_pheno.RData')
cell.pro = read.table('CBMAP_all_sample_methylation_ctp_with_clr_deconvolution.txt',header=T)
cell_aligned <- cell.pro[pheno_df$barcode_id, , drop = FALSE]
pheno_df <- cbind(pheno_df, cell_aligned)

# methylation matrix
Mdat <- h5read("final_DNAm_Mdat.h5",'Mdat')
colnames(Mdat) = h5read("final_DNAm_Mdat.h5",'colnames')
rownames(Mdat) = h5read("final_DNAm_Mdat.h5",'rownames')
cg <- 2^Mdat / (2^Mdat + 1)
cg <- t(cg)

sample = intersect(rownames(cg), rownames(pheno_df))
cg <- cg[sample,]
pheno_df <- pheno_df[sample,]
pheno_df$sex_male <- as.factor(pheno_df$sex_male)
pheno_df$Sample_Plate <- as.factor(pheno_df$Sample_Plate)
pheno_df$bank <- as.factor(pheno_df$bank)
pheno_df[] <- lapply(pheno_df, function(x)
  if (is.numeric(x)) replace(x, is.na(x), median(x, na.rm = TRUE)) else x
)
identical(rownames(pheno_df), rownames(cg))

residuals_matrix <- matrix(NA, nrow = nrow(cg), ncol = ncol(cg))
colnames(residuals_matrix) <- colnames(cg)
rownames(residuals_matrix) <- rownames(cg)
for (cpg in colnames(cg)) {
  lm_model <- lm(cg[, cpg] ~ sex_male + bank + PMD + RIN + Exc + Inh + NonN_Astro_FGF3R + NonN_Endo + NonN_Micro + NonN_Oligo_MBP + NonN_OPC + Sample_Plate, data = pheno_df)
  residuals_matrix[, cpg] <- residuals(lm_model)
}
h5write(residuals_matrix, "residual_methy", "residual_methy")
h5write(rownames(residuals_matrix), "residual_methy", "rownames")
h5write(colnames(residuals_matrix), "residual_methy", "colnames")

# Calculate the mean methylation level for group age
unique_ages <- sort(unique(pheno_df$age))
result_mat <- matrix(NA, nrow = length(unique_ages), ncol = ncol(residual_methy))
rownames(result_mat) <- unique_ages
colnames(result_mat) <- colnames(residual_methy)
for (i in seq_along(unique_ages)) {
  age_i <- unique_ages[i]
  idx <- which(pheno_df$age == age_i)
  result_mat[i, ] <- colMeans(residual_methy[idx, , drop = FALSE], na.rm = TRUE)
}

####identify nonlinear cpgs####
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

cpg_neg = results[results$r2_lm > 0.5 & results$beta_lm < 0,]
cpg_pos = results[results$r2_lm > 0.5 & results$beta_lm > 0,]
delete_cpg = c(cpg_neg$cpg, cpg_pos$cpg)
results = results[which(! results$cpg %chin% delete_cpg),]
cpg = results[which(results$r2_loess-results$r2_lm>0.5),]
colnames(cpg) <- c('estimate_linear','p_linear','r2_linear','r2_loess','cpg','fdr_linear')
write.xlsx(cpg,'nonlinear_cpgs.xlsx')

####k-means_cluster####
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
cg <- cg[,cpg$cpg]
cg <- data.frame(cg)
cg <- data.frame(scale(cg))
cg$ID <- rownames(cg)
cg$age <- pheno_df[match(cg$ID, pheno_df$Sample_Name),'age']
cg$sex <- pheno_df[match(cg$ID, pheno_df$Sample_Name),'sex_male']
cg$sex <- factor(cg$sex)
cg <- cg %>% arrange(age)

cg.tf <- gather(cg, 1:(ncol(cg)-3), key = "cpg", value = "BetaValue")

# ###########  cluster  #################
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
}
loess.mat <- t(res_summary)
colnames(loess.mat) <- cg$ID

##elbow method
k_values <- 1:10
wss_values <- numeric(length(k_values))

set.seed(42) 
for (i in seq_along(k_values)) {
  k <- k_values[i]
  km <- kmeans(loess.mat, centers = k, nstart = 10) 
  wss_values[i] <- km$tot.withinss
  cat("k =", k, ", WSS =", wss_values[i], "\n")
}

elbow_data <- data.frame(
  k = k_values,
  WSS = wss_values
)
elbow_plot <- ggplot(elbow_data, aes(x = k, y = WSS)) +
  geom_point(size = 3, color = "#1F78B4") +
  geom_line(color = "#1F78B4", linewidth = 1) +
  scale_x_continuous(breaks = k_values) +
  theme_classic(base_size = 15) +
  labs(
    title = "",
    x = "Number of Clusters (k)",
    y = "Total Within-Cluster Sum of Squares (WSS)"
  ) +
  theme(axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18))

find_elbow_k <- function(k_values, wss_values) {
  k_norm <- (k_values - min(k_values)) / (max(k_values) - min(k_values))
  wss_norm <- (wss_values - min(wss_values)) / (max(wss_values) - min(wss_values))
  p1 <- c(k_norm[1], wss_norm[1])
  p2 <- c(k_norm[length(k_norm)], wss_norm[length(wss_norm)])
  A <- p1[2] - p2[2]
  B <- p2[1] - p1[1]
  C <- p1[1]*p2[2] - p2[1]*p1[2]
  distances <- numeric(length(k_values))
  for (i in seq_along(k_values)) {
    p0 <- c(k_norm[i], wss_norm[i])
    distances[i] <- abs(A*p0[1] + B*p0[2] + C) / sqrt(A^2 + B^2)
  }
  optimal_k <- k_values[which.max(distances)]
  return(optimal_k)
}

optimal_k <- find_elbow_k(elbow_data$k, elbow_data$WSS)
kmeans_result <- kmeans(loess.mat, centers = optimal_k)

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
cg <- cg[,cpg$cpg]
cg <- data.frame(cg)
# cg <- data.frame(scale(cg))
cg$ID <- rownames(cg)
cg$age <- pheno_df[match(cg$ID, pheno_df$Sample_Name),'age']
cg$sex <- pheno_df[match(cg$ID, pheno_df$Sample_Name),'sex_male']
cg$sex <- factor(cg$sex)
cg <- cg %>% arrange(age)

cg.tf <- gather(cg, 1:(ncol(cg)-3), key = "cpg", value = "BetaValue")
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

clustering.res.df <- data.frame(cpg = names(kmeans_result$cluster), Cluster = kmeans_result$cluster)
cg_df.cluster <- left_join(cg.tf, clustering.res.df, by = c("cpg" = "cpg"))

cluster.count <- cg_df.cluster %>% dplyr::count(Cluster)
cluster.count$n <- cluster.count$n/1031

# Convert Cluster to factor for plotting
cg_df.cluster$Cluster <- as.factor(cg_df.cluster$Cluster)
sampled_cpg = list()
for (i in 1:2) {
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

cluster_gg <- list()
for (cluster in 1:2) {
  cluster_data <- filter(cg_df.cluster.sample, Cluster == cluster)
  cluster_gg[[cluster]] <- ggplot(cluster_data, aes(x = age, y = BetaValue)) +
    theme_classic(base_size = 20) +
    theme(aspect.ratio = 1,  legend.position = "none",
          plot.title = element_text(size = 18)) +
    geom_smooth(method = "loess", method.args = list(span = .75, degree = 1),aes(group = cpg), se = FALSE, color = col6[cluster], linewidth = 0.1) +
    geom_smooth(method = "loess", method.args = list(span = .75, degree = 1),se = FALSE, color = col62[cluster], linewidth = 0.5) +
    ylab("Methylation level") + xlab("Chronological age") +
    ggtitle(paste0("Cluster ", cluster, " (", cluster.count$n[cluster], " CpGs)"))
}

p3 = cluster_gg[[1]]|cluster_gg[[2]]
ggsave(plot = p3,filename = 'r2diff2_kmeans_cluster.pdf')
