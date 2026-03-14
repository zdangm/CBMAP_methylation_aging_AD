library(rhdf5)
library(data.table)
library(dplyr)
library(readxl)
library(stringr)

load('ROSMAP_pheno.RData')
pheno_df$projid <- as.character(str_pad(pheno_df$projid, width = 8, pad = "0", side = "left"))
hml_info <- read_excel('dataset_1495_cross-sectional_02-16-2025.xlsx')
cogdx_1 <- hml_info %>% filter(
  cogdx %in% c(4,5,6)) %>% pull(projid) %>% unique()
pheno_df$age_death <- as.numeric(unlist(hml_info[match(pheno_df$projid, hml_info$projid), 'age_death']))
pheno_df$msex <- as.factor(pheno_df$msex)
pheno_df$Sample_Plate <- as.factor(pheno_df$Sample_Plate)
pheno_df$batch <- as.factor(pheno_df$batch)
pheno_df$Study <- as.factor(pheno_df$Study)
pheno_df <- pheno_df[! pheno_df$projid %chin% cogdx_1,]

idkey = read.csv('ROSMAP_IDkey.csv', header=T)
rin_info = read.csv('ROSMAP_assay_rnaSeq_metadata.csv',header=T)
rin_info = rin_info[,c("specimenID","RIN")]
rin_info$projid = unlist(idkey[match(rin_info$specimenID,idkey$rnaseq_id),"projid"])
rin_info$projid <- as.character(str_pad(rin_info$projid, width = 8, pad = "0", side = "left"))
pheno_df$rin = unlist(rin_info[match(pheno_df$projid,rin_info$projid),"RIN"]) 

cell.pro = read.table('ROSMAP_all_sample_methylation_ctp_with_clr_deconvolution.txt',header=T)
cell_aligned <- cell.pro[pheno_df$barcode_id, , drop = FALSE]
pheno_df <- cbind(pheno_df, cell_aligned)
rownames(pheno_df) = pheno_df$mwas_id

# prepare methylatioSample_Plate# prepare methylation matrix
cg <- h5read("ROSMAP_methylaition/final_DNAm_invMdat.h5",'Mdat')
colnames(cg) = h5read("/ROSMAP_methylaition/final_DNAm_invMdat.h5", 'colnames')
rownames(cg) = h5read("/ROSMAP_methylaition/final_DNAm_invMdat.h5", 'rownames')

sample = intersect(rownames(cg), pheno_df$mwas_id) 
cg = cg[sample,]
pheno_df = pheno_df[sample,]
pheno_df[] <- lapply(pheno_df, function(x)
  if (is.numeric(x)) replace(x, is.na(x), median(x, na.rm = TRUE)) else x
)
identical(row.names(pheno_df), rownames(cg))

# Linear regression by cpgs Methylation
run_linear_regression <- function(methylation_data, covariates) {
  data <- data.frame(methylation_level = methylation_data, covariates)
  model <- lm(methylation_level ~ age_death + msex + Study + pmi + rin + Exc + Inh + NonN_Astro_FGF3R + NonN_Endo + NonN_Micro + NonN_Oligo_MBP + NonN_OPC + Sample_Plate + batch, data = data)
  summary_model <- summary(model)
  return(summary_model$coefficients[2,])
}
analyze_methylation <- function(cg, pheno_df) {
  results <- apply(cg, 2, function(x) {
    methylation_data <- x  
    covariates <- pheno_df %>% select(age_death, msex, Study, pmi, rin, Exc, Inh, NonN_Astro_FGF3R, NonN_Endo, NonN_Micro, NonN_Oligo_MBP, NonN_OPC, Sample_Plate, batch)  
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
save(res, file='ROSMAP_meth_age_noDementia.RData')
