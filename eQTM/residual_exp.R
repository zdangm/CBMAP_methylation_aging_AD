library(data.table)
library(stringr)
library(dplyr)
library(rhdf5)
library(compositions)

exp = read.table('/CBMAP_RNAseq_protein_coding/RNAseq_process/RNAseq_data_qn_comb.txt',check.names = F)
exp = t(exp)
pca = fread('/eQTM/exp_pca_10.txt')
pca = as.data.frame(pca)
rownames(pca) = pca$V1
pca = pca[,c(2:7)] # top 6 PCs
cell.pro = read.table('/CBMAP_RNAseq_protein_coding/CIBERSORT/CIBERSORT_prop_MIND.txt')
cell.pro<-as.matrix(cell.pro)
cell.pro[which(cell.pro <= 0)] <- 1e-3
cell.pro<-as.data.frame(clr(cell.pro))
pheno_df = fread('/CBMAP_RNAseq_protein_coding/RNAseq_process/RNAseq_metadata.txt')
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
pheno_df = pheno_df[,c('Bank','sex_male','age','PMD','RIN','Astrocyte','Excitatory','Inhibitory','Oligo','Endothelial',paste0('V',2:7))] # V2-V7 columns represent the top 6 PCs of gene expression
pheno_df$Bank = factor(pheno_df$Bank)
pheno_df$sex_male = factor(pheno_df$sex_male)
miss_doc = apply(pheno_df, 2, function(x){sum(is.na(x))})
# Bank    sex_male         age         PMD         RIN   Astrocyte  Excitatory 
#    0           0           0           0           0          33          33 
# Inhibitory       Oligo Endothelial          V2          V3          V4          V5 
#         33          33          33           0           0           0           0 
# V6          V7 
#  0           0

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
fwrite(residuals_df,file='/eQTM/residual_exp.txt',col.names = T,row.names = T,quote=F,sep='\t')
