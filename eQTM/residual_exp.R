library(data.table)
library(stringr)
library(dplyr)
library(rhdf5)

exp = read.table('/RNAseq_process/RNAseq_data_qn_comb.txt',check.names = F)
exp = t(exp)

pca = fread('exp_pca_10.txt')
pca = as.data.frame(pca)
rownames(pca) = pca$V1
pca = pca[,c(2:7)] 
cell.pro = read.table('/CIBERSORT/CIBERSORT_prop_MIND.txt')
pheno_df = fread('/RNAseq_process/RNAseq_metadata.txt')
pheno_df = as.data.frame(pheno_df)
rownames(pheno_df) = pheno_df$sample_id
cell.pro = cell.pro[pheno_df$sample_id, , drop = FALSE]
pheno_df = cbind(pheno_df,cell.pro)
pca = pca[pheno_df$sample_id, , drop = FALSE]
pheno_df = cbind(pheno_df,pca)
pheno_df = pheno_df[pheno_df$SexMismatch == 0,]

sample = intersect(rownames(exp),pheno_df$sample_id)
exp = exp[sample,]
pheno_df = pheno_df[sample,]
pheno_df = pheno_df[,c(2:4,6:7,64:75)]
pheno_df$Bank = factor(pheno_df$Bank)
pheno_df$sex_male = factor(pheno_df$sex_male)
pheno_df[] <- lapply(pheno_df, function(x)
  if (is.numeric(x)) replace(x, is.na(x), median(x, na.rm = TRUE)) else x
)
identical(rownames(pheno_df),rownames(exp))

residuals_matrix <- matrix(NA, nrow = nrow(exp), ncol = ncol(exp))
colnames(residuals_matrix) <- colnames(exp)
rownames(residuals_matrix) <- rownames(exp)
for (gene in colnames(exp)) {
  lm_model <- lm(exp[, gene] ~ ., data = pheno_df)
  residuals_matrix[, gene] <- residuals(lm_model)
}

residuals_df <- as.data.frame(residuals_matrix)
fwrite(residuals_df,file='residual_exp.txt',col.names = T,row.names = T,quote=F,sep='\t')
