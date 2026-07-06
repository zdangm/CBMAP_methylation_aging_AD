library(data.table)
library(dplyr)
library(rhdf5)

setwd('/methylation/methy_AD_revision2/results/EWAS')
load('/methylation/data/GSE40279/processed/pheno_info.RData')
pheno_df$AGE <- as.numeric(pheno_df$AGE)
pheno_df = pheno_df[which(pheno_df$AGE > 60),] # Comment out this line if using all samples
rownames(pheno_df) <- pheno_df$sampleID
cg <- h5read("/methylation/data/GSE40279/processed/final_DNAm_invMdat.h5",'Mdat')
colnames(cg) = h5read("/methylation/data/GSE40279/processed/final_DNAm_invMdat.h5",'colnames')
rownames(cg) = h5read("/methylation/data/GSE40279/processed/final_DNAm_invMdat.h5",'rownames')
cg = cg[,grep("cg",colnames(cg))] 

sample = intersect(rownames(cg), rownames(pheno_df)) 
cg <- cg[sample,]
pheno_df <- pheno_df[sample,]
pheno_df$sex_male <- ifelse(pheno_df$Gender == 'M', 1, 0)
pheno_df$sex_male <- as.factor(pheno_df$sex_male)
pheno_df$Sample_Plate <- as.factor(pheno_df$Plate)
pheno_df$Ethnicity <- as.factor(pheno_df$Ethnicity)
identical(rownames(pheno_df), rownames(cg))
pheno_df = pheno_df %>% select(AGE, sex_male, Sample_Plate, Ethnicity, B, NK, CD4T, CD8T, Mono, Neutro, Eosino)  
apply(pheno_df, 2, function(x){sum(is.na(x))})
# AGE     sex_male Sample_Plate    Ethnicity            B           NK         CD4T 
#   0            0            0            0            0            0            0 
# CD8T         Mono       Neutro       Eosino 
#    0            0            0            0 

run_linear_regression <- function(methylation_data, covariates) {
  data <- data.frame(methylation_level = methylation_data, covariates)
  model <- lm(methylation_level ~ AGE + sex_male + Sample_Plate + Ethnicity + B + NK + CD4T + CD8T + Mono + Neutro, data = data)
  summary_model <- summary(model)
  return(summary_model$coefficients[2,])
}
analyze_methylation <- function(cg, pheno_df) {
  results <- apply(cg, 2, function(x) {
    methylation_data <- x  
    covariates <- pheno_df %>% select(AGE, sex_male, Sample_Plate, Ethnicity, B, NK, CD4T, CD8T, Mono, Neutro, Eosino)  
    result <- run_linear_regression(methylation_data, covariates)  
    return(result)
  })
  return(results)
}

results = analyze_methylation(cg, pheno_df)
res <- as.data.frame(t(results[c(1,2,4),]),stringsAsFactors = FALSE)
colnames(res) <- c("Estimate", "StdErr", "pValue")
res$cpg <- colnames(results)
res <- res[,c("cpg", "Estimate", "StdErr", "pValue")]
res$fdr <- p.adjust(res$pValue,'BH')
save(res, file='Hannum_meth_age_60years.RData')

