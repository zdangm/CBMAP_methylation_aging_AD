library(data.table)
library(dplyr)
library(rhdf5)
setwd('/methylation/methy_AD_revision2/results/EWAS')

pheno_df = fread('/methylation/data/GSE111629/sample_info.txt',header=T)
pheno_df$sample_name = paste0(pheno_df$geo_accession,'_',pheno_df$source_name)
pheno_df = as.data.frame(pheno_df)
pheno_df = pheno_df[which(pheno_df$age>60),] # Comment out this line if using all samples
rownames(pheno_df) <- pheno_df$sample_name
celltype = read.table('/methylation/data/GSE111629/celltype_with_clr.txt',header=T)
pheno_df = merge(pheno_df, celltype, by='row.names')
rownames(pheno_df) = pheno_df$Row.names
pheno_df = pheno_df[which(pheno_df$disease == "PD-free control"),] # Comment out this line if using all samples

cg <- h5read("/methylation/data/GSE111629/final_DNAm_invMdat.h5",'Mdat')
colnames(cg) = h5read("/methylation/data/GSE111629/final_DNAm_invMdat.h5",'colnames')
rownames(cg) = h5read("/methylation/data/GSE111629/final_DNAm_invMdat.h5",'rownames')
cg = cg[,grep("cg",colnames(cg))] 

sample = intersect(rownames(cg), rownames(pheno_df)) 
cg <- cg[sample,]
pheno_df <- pheno_df[sample,]
pheno_df$age <- as.numeric(pheno_df$age)
pheno_df$gender <- as.factor(pheno_df$gender)
pheno_df$ethnicity <- as.factor(pheno_df$ethnicity)
identical(rownames(pheno_df), rownames(cg))
pheno_df <- pheno_df %>% select(age, gender, ethnicity, B, NK, CD4T, CD8T, Mono, Neutro, Eosino)  
apply(pheno_df, 2, function(x){sum(is.na(x))})
# age    gender ethnicity         B        NK      CD4T      CD8T      Mono    Neutro 
#   0         0         0         0         0         0         0         0         0 
# Eosino
#      0

run_linear_regression <- function(methylation_data, covariates) {
  data <- data.frame(methylation_level = methylation_data, covariates)
  model <- lm(methylation_level ~ age + gender + ethnicity + B + NK + CD4T + CD8T + Mono + Neutro, data = data)
  summary_model <- summary(model)
  return(summary_model$coefficients[2,])
}
analyze_methylation <- function(cg, pheno_df) {
  results <- apply(cg, 2, function(x) {
    methylation_data <- x  
    covariates <- pheno_df %>% select(age, gender, ethnicity, B, NK, CD4T, CD8T, Mono, Neutro, Eosino)  
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
save(res, file='PEG1_meth_age_60years.RData')
