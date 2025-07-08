library(rhdf5)
library(data.table)
library(dplyr)
library(readxl)
library(stringr)

load('ROSMAP_pheno.RData')
pheno_df$projid <- as.character(str_pad(pheno_df$projid, width = 8, pad = "0", side = "left"))
hml_info <- read_excel('dataset_1495_cross-sectional_02-16-2025.xlsx')
pheno_df$age_death <- as.numeric(unlist(hml_info[match(pheno_df$projid, hml_info$projid), 'age_death']))
pheno_df$msex <- as.factor(pheno_df$msex)
pheno_df$Sample_Plate <- as.factor(pheno_df$Sample_Plate)
pheno_df$batch <- as.factor(pheno_df$batch)
pheno_df$NeuN_pos <- as.numeric(pheno_df$NeuN_pos)

# prepare methylation matrix
Mdat <- h5read('BMIQ_Mdat_combat.h5', 'Mdat')
colnames(Mdat) = h5read('BMIQ_Mdat_combat.h5','colnames')
rownames(Mdat) = h5read('BMIQ_Mdat_combat.h5','rownames')
cg <- 2^Mdat / (2^Mdat + 1)
cg <- t(cg)
rownames(cg) <- pheno_df[match(rownames(cg), pheno_df$barcode_id), 'mwas_id']
sample = intersect(rownames(cg), pheno_df$mwas_id) 
cg = cg[sample,]
row.names(pheno_df) = pheno_df$mwas_id
pheno_df = pheno_df[sample,]
identical(row.names(pheno_df), rownames(cg))

# Linear regression by cpgs Methylation
run_linear_regression <- function(methylation_data, covariates) {
  data <- data.frame(methylation_level = methylation_data, covariates)
  model <- lm(methylation_level ~ age_death + msex + NeuN_pos + Sample_Plate + batch, data = data)
  summary_model <- summary(model)
  return(summary_model$coefficients[2,])
}
analyze_methylation <- function(cg, pheno_df) {
  results <- apply(cg, 2, function(x) {
    methylation_data <- x  
    covariates <- pheno_df %>% select(age_death, msex, NeuN_pos, Sample_Plate, batch)  
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
save(res, file='ROSMAP_meth_age.RData')