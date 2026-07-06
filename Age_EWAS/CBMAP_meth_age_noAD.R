library(rhdf5)
library(data.table)
library(dplyr)

setwd('/methylation/methy_AD_revision2/results/EWAS')
load('/methylation/results/disease_meth_age/AD_meth/CBMAP_pheno.RData')
cg <- h5read("/methylation/data/CBMAP/DNAm_processed/final_DNAm_invMdat.h5",'Mdat')
colnames(cg) = h5read("/methylation/data/CBMAP/DNAm_processed/final_DNAm_invMdat.h5",'colnames')
rownames(cg) = h5read("/methylation/data/CBMAP/DNAm_processed/final_DNAm_invMdat.h5",'rownames')
cg = cg[,grep("cg",colnames(cg))] 

cell.pro = read.table('/methylation/results/ctp_analysis/CBMAP_all_sample_methylation_ctp_with_clr_deconvolution.txt',header=T)
cell_aligned <- cell.pro[pheno_df$barcode_id, , drop = FALSE]
pheno_df <- cbind(pheno_df, cell_aligned)

sample = intersect(rownames(cg), rownames(pheno_df)) 
cg <- cg[sample,]
pheno_df <- pheno_df[sample,]
pheno_df$sex_male <- as.factor(pheno_df$sex_male)
pheno_df$Sample_Plate <- as.factor(pheno_df$Sample_Plate)
pheno_df$bank <- as.factor(pheno_df$bank)
pheno_df[] <- lapply(pheno_df, function(x)
  if (is.numeric(x)) replace(x, is.na(x), median(x, na.rm = TRUE)) else x
)

# filter samples 
all_sample_info <- as.data.frame(read.csv("CBMAP_sample_info.csv",header=T)) 
selected_sample = all_sample_info$id[which(all_sample_info$diag_dementia == 1)]
pheno_df = pheno_df[which(! pheno_df$Sample_Name %chin% selected_sample),] 
cg = cg[pheno_df$Sample_Name,]
identical(rownames(pheno_df), rownames(cg))

run_linear_regression <- function(methylation_data, covariates) {
  data <- data.frame(methylation_level = methylation_data, covariates)
  model <- lm(methylation_level ~ age + sex_male + bank + PMD + RIN + Exc + Inh + NonN_Astro_FGF3R + NonN_Endo + NonN_Micro + NonN_Oligo_MBP + Sample_Plate, data = data)
  summary_model <- summary(model)
  return(summary_model$coefficients[2,])
}
analyze_methylation <- function(cg, pheno_df) {
  results <- apply(cg, 2, function(x) {
    methylation_data <- x  
    covariates <- pheno_df %>% select(age, sex_male, bank, PMD, RIN, Exc, Inh, NonN_Astro_FGF3R, NonN_Endo, NonN_Micro, NonN_Oligo_MBP, NonN_OPC, Sample_Plate)  
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
save(res, file='CBMAP_meth_age_no_dementia_947.RData')
