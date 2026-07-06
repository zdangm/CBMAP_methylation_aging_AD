library(data.table)
library(stringr)
library(dplyr)
library(rhdf5)

load('CBMAP_pheno.RData')

cell.pro = read.table('CBMAP_all_sample_methylation_ctp_with_clr_deconvolution.txt',header=T)
cell_aligned <- cell.pro[pheno_df$barcode_id, , drop = FALSE]
pheno_df <- cbind(pheno_df, cell_aligned)

pca = fread('/eQTM/DNAm_pca_10.txt')
pca = as.data.frame(pca)
rownames(pca) = pca$V1
pca = pca[,c(2:7)]
pca_aligned <- pca[pheno_df$Sample_Name, , drop = FALSE]
pheno_df <- cbind(pheno_df, pca_aligned)

# methylation matrix
cg <- h5read("/data/CBMAP_methylaition/final_DNAm_invMdat.h5",'Mdat')
colnames(cg) = h5read("/data/CBMAP_methylaition/final_DNAm_invMdat.h5",'colnames')
rownames(cg) = h5read("/data/CBMAP_methylaition/final_DNAm_invMdat.h5",'rownames')
cg = cg[,grep("cg",colnames(cg))] 
sample = intersect(rownames(cg), rownames(pheno_df)) 
cg <- cg[sample,]
pheno_df <- pheno_df[sample,]
pheno_df$sex_male <- as.factor(pheno_df$sex_male)
pheno_df$Sample_Plate <- as.factor(pheno_df$Sample_Plate)
pheno_df$bank <- as.factor(pheno_df$bank)
miss_doc = apply(pheno_df, 2, function(x){sum(is.na(x))})
# barcode_id      Sample_Name       Sentrix_ID Sentrix_Position     Sample_Plate 
#          0                0                0                0                0 
# bank         sex_male              age             year              PMD 
#    0                0                0                0                5 
# RIN         NeuN_pos         NeuN_neg              Exc              Inh 
#   0                0                0               22               22 
# NonN_Astro_FGF3R        NonN_Endo       NonN_Micro   NonN_Oligo_MBP         NonN_OPC 
#               22               22               22               22               22 
# NeuN_pos         NeuN_neg               V2               V3               V4 
#       22               22                0                0                0 
# V5               V6               V7 
#  0                0                0
# Note V2-V7 are the top 6 PCs of methylation
pheno_df[] <- lapply(pheno_df, function(x)
  if (is.numeric(x)) replace(x, is.na(x), median(x, na.rm = TRUE)) else x
)
pheno_df = pheno_df[,c('Sample_Plate','bank','sex_male','age','PMD','RIN','Exc','Inh','NonN_Astro_FGF3R','NonN_Endo','NonN_Micro','NonN_Oligo_MBP', paste0('V',2:7))] # V2-V7 columns represent the top 6 PCs of methylation
identical(rownames(pheno_df), rownames(cg))

residuals_matrix <- matrix(NA, nrow = nrow(cg), ncol = ncol(cg))
colnames(residuals_matrix) <- colnames(cg)
rownames(residuals_matrix) <- rownames(cg)
for (cpg in colnames(cg)) {
  lm_model <- lm(cg[, cpg] ~ ., data = pheno_df)
  residuals_matrix[, cpg] <- residuals(lm_model)
}

h5write(residuals_matrix, "residual_methy", "residual_methy")
h5write(rownames(residuals_matrix), "residual_methy", "rownames")
h5write(colnames(residuals_matrix), "residual_methy", "colnames")

