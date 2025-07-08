library(rhdf5)
library(data.table)
library(dplyr)

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
pheno_df$sex_male <- as.factor(pheno_df$sex_male)
pheno_df$Sample_Plate <- as.factor(pheno_df$Sample_Plate)
identical(rownames(pheno_df), rownames(cg))

# Linear regression by cpgs Methylation
run_linear_regression <- function(methylation_data, covariates) {
  data <- data.frame(methylation_level = methylation_data, covariates)
  model <- lm(methylation_level ~ age + sex_male + NeuN_pos + Sample_Plate, data = data)
  summary_model <- summary(model)
  return(summary_model$coefficients[2,])
}
analyze_methylation <- function(cg, pheno_df) {
  results <- apply(cg, 2, function(x) {
    methylation_data <- x  
    covariates <- pheno_df %>% select(age, sex_male, NeuN_pos, Sample_Plate)  
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
save(res, file='CBMAP_meth_age.RData')
