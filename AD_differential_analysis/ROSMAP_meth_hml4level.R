library(data.table)
library(dplyr)
library(rhdf5)
library(readxl)
library(stringr)

load('ROSMAP_pheno.RData')
pheno_df$projid <- as.character(str_pad(pheno_df$projid, width = 8, pad = "0", side = "left"))
hml_info <- read_excel('dataset_1495_cross-sectional_02-16-2025.xlsx')
metadata_long <- read_excel('dataset_1495_long_12-01-2024.xlsx') %>% as.data.frame()
PD_case <- metadata_long %>% filter(
  r_pd %in% c(1)) %>% pull(projid) %>% unique()
PD_0 <- unique(hml_info$projid[which(! hml_info$projid %chin% PD_case)])
LBD_0 <- hml_info %>% filter(
  dlbdx %in% c(0)) %>% pull(projid) %>% unique()
cogdx_0 <- hml_info %>% filter(
  cogdx %in% c(1,2,3)) %>% pull(projid) %>% unique()
case = hml_info %>% filter(
  ADNC_4level %in% c(1,2,3)) %>% pull(projid) %>% unique()
control = hml_info %>% filter(
  ADNC_4level %in% c(0) 
    & (projid %in% PD_0) 
    & (projid %in% LBD_0)
    & (projid %in% cogdx_0)) %>% pull(projid) %>% unique()
sample = intersect(c(case, control), pheno_df$projid) 
pheno_df = pheno_df[which(pheno_df$projid %chin% sample),]
pheno_df$age_death <- as.numeric(unlist(hml_info[match(pheno_df$projid, hml_info$projid), 'age_death']))
pheno_df$hml4 <- as.numeric(unlist(hml_info[match(pheno_df$projid, hml_info$projid), 'ADNC_4level']))
pheno_df$msex <- as.factor(pheno_df$msex)
pheno_df$Sample_Plate <- as.factor(pheno_df$Sample_Plate)
pheno_df$batch <- as.factor(pheno_df$batch)
pheno_df$NeuN_pos <- as.numeric(pheno_df$NeuN_pos)
pheno_df <- pheno_df[!is.na(pheno_df$hml4),] 

# prepare methylation matrix
cg <- h5read("final_DNAm_invMdat.h5",'Mdat')
colnames(cg) = h5read("final_DNAm_invMdat.h5", 'colnames')
rownames(cg) = h5read("final_DNAm_invMdat.h5", 'rownames')
sample = intersect(rownames(cg), pheno_df$mwas_id) 
cg = cg[sample,]
cg = cg[,grep("cg",colnames(cg))] # keep only probes that start with "cg"
row.names(pheno_df) = pheno_df$mwas_id
pheno_df = pheno_df[sample,]
identical(row.names(pheno_df), rownames(cg))

# Linear regression by cpgs Methylation
run_linear_regression <- function(methylation_data, covariates) {
  data <- data.frame(methylation_level = methylation_data, covariates)
  model <- lm(methylation_level ~ hml4 + msex + age_death + NeuN_pos + Sample_Plate + batch, data = data)
  summary_model <- summary(model)
  return(summary_model$coefficients[2,])
}
analyze_methylation <- function(cg, pheno_df) {
  results <- apply(cg, 2, function(x) {
    methylation_data <- x  
    covariates <- pheno_df %>% select(hml4, age_death, msex, NeuN_pos, Sample_Plate, batch)  
    result <- run_linear_regression(methylation_data, covariates)  
    return(result)
  })
  return(results)
}

results = analyze_methylation(cg, pheno_df)
save(results, file='ROSMAP_meth_hml4level.RData')