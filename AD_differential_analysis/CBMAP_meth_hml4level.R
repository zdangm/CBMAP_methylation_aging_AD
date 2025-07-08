library(rhdf5)
library(data.table)
library(dplyr)
library(readxl)

load('CBMAP_pheno.RData')
# select case and control
all_sample_info <- as.data.frame(read.csv("CBMAP_sample_info.csv",header=T)) 
selected_sample = all_sample_info[!is.na(all_sample_info[,'ADNC']),]
selected_sample <- selected_sample %>% 
  filter(
    diag_ALS == 0 & 
      diag_SCZ == 0 & 
      diag_epilepsy == 0 &
      diag_others == 0 &
      OTHER == 0) 
complete_rows <- complete.cases(selected_sample[, c("ADNC", "ADNC.L", "ADNC.M", "ADNC.H")])
selected_sample <- selected_sample[complete_rows, ]
selected_sample$marker <- ifelse(selected_sample$ADNC == 0, 0,
                                 ifelse(selected_sample$ADNC.L == 1, 1,
                                        ifelse(selected_sample$ADNC.M == 1, 2,
                                               ifelse(selected_sample$ADNC.H == 1, 3, NA))))
exclued_id <- selected_sample[which(selected_sample$marker == 0 & (selected_sample$LBD == 1 | selected_sample$ARTAG_Gray.matter == 1 |
                                           selected_sample$ARTAG_Perivascular == 1 | selected_sample$ARTAG_Subpial == 1 | selected_sample$diag_PD == 1 | selected_sample$diag_dementia == 1)),'id']
pheno_df$AD <- selected_sample[match(pheno_df$Sample_Name, selected_sample$id), 'marker']
pheno_df <- pheno_df[!is.na(pheno_df$AD),] 
pheno_df <- pheno_df[which(! pheno_df$Sample_Name %chin% exclued_id),]
dim(pheno_df)
# methylation matrix
cg <- h5read("final_DNAm_invMdat.h5",'Mdat')
colnames(cg) = h5read("final_DNAm_invMdat.h5",'colnames')
rownames(cg) = h5read("final_DNAm_invMdat.h5",'rownames')
cg = cg[,grep("cg",colnames(cg))] # keep only probes that start with "cg"

sample = intersect(rownames(cg), rownames(pheno_df)) 
cg <- cg[sample,]
pheno_df <- pheno_df[sample,]
pheno_df$sex_male <- as.factor(pheno_df$sex_male)
pheno_df$Sample_Plate <- as.factor(pheno_df$Sample_Plate)
identical(row.names(pheno_df), rownames(cg))


# Linear regression by cpgs Methylation
run_linear_regression <- function(methylation_data, covariates) {
  data <- data.frame(methylation_level = methylation_data, covariates)
  model <- lm(methylation_level ~ AD + sex_male + age + NeuN_pos + Sample_Plate, data = data)
  summary_model <- summary(model)
  return(summary_model$coefficients[2,])
}
analyze_methylation <- function(cg, pheno_df) {
  results <- apply(cg, 2, function(x) {
    methylation_data <- x  
    covariates <- pheno_df %>% select(AD, age, sex_male, NeuN_pos, Sample_Plate)  
    result <- run_linear_regression(methylation_data, covariates)  
    return(result)
  })
  return(results)
}

results = analyze_methylation(cg, pheno_df)
save(results, file='CBMAP_meth_hml4level.RData')
