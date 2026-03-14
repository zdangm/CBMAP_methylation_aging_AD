library(data.table)
library(dplyr)
library(rhdf5)

pheno_df = fread('/GSE111629/sample_info.txt',header=T)
pheno_df$sample_name = paste0(pheno_df$geo_accession,'_',pheno_df$source_name)
pheno_df = as.data.frame(pheno_df)
pheno_df = pheno_df[which(pheno_df$age>60),]
rownames(pheno_df) <- pheno_df$sample_name
celltype = read.table('/GSE111629/celltype_with_clr.txt',header=T)
pheno_df = merge(pheno_df, celltype, by='row.names')
rownames(pheno_df) = pheno_df$Row.names

cg <- h5read("/GSE111629/final_DNAm_invMdat.h5",'Mdat')
colnames(cg) = h5read("/GSE111629/final_DNAm_invMdat.h5",'colnames')
rownames(cg) = h5read("/GSE111629/final_DNAm_invMdat.h5",'rownames')

sample = intersect(rownames(cg), rownames(pheno_df)) 
cg <- cg[sample,]
pheno_df <- pheno_df[sample,]
pheno_df$age <- as.numeric(pheno_df$age)
pheno_df$gender <- as.factor(pheno_df$gender)
pheno_df$ethnicity <- as.factor(pheno_df$ethnicity)
pheno_df$disease <- as.factor(pheno_df$disease)
identical(rownames(pheno_df), rownames(cg))

# Linear regression by cpgs Methylation
run_linear_regression <- function(methylation_data, covariates) {
  data <- data.frame(methylation_level = methylation_data, covariates)
  model <- lm(methylation_level ~ age + gender + ethnicity + B + NK + CD4T + CD8T + Mono + Neutro + Eosino, data = data)
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
save(res, file='PEG1_meth_age.RData')