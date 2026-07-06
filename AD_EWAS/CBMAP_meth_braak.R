library(rhdf5)
library(data.table)
library(dplyr)
library(readxl)

load('CBMAP_pheno.RData')
all_sample_info <- as.data.frame(read.csv("CBMAP_sample_info.csv",header=T)) 
selected_sample = all_sample_info[!is.na(all_sample_info[,'Braak.NFT.stage']),]
selected_sample <- selected_sample %>% 
  filter(
    diag_ALS == 0 & 
      diag_SCZ == 0 & 
      diag_epilepsy == 0 &
      diag_others == 0 &
      OTHER == 0) 
braak_sample <- unlist(selected_sample %>% filter(!is.na(Braak.NFT.stage))%>% dplyr::select(id)) # 1039 samples

# methylation matrix
cg <- h5read("/data/CBMAP_methylaition/final_DNAm_invMdat.h5",'Mdat')
colnames(cg) = h5read("/data/CBMAP_methylaition/final_DNAm_invMdat.h5",'colnames')
rownames(cg) = h5read("/data/CBMAP_methylaition/final_DNAm_invMdat.h5",'rownames')
cg = cg[,grep("cg",colnames(cg))] 

sample = Reduce(intersect,list(rownames(cg),braak_sample,rownames(pheno_df))) 
cg <- cg[sample,]
pheno_df <- merge(pheno_df, selected_sample[,c('id','Braak.NFT.stage')], by.x='Sample_Name', by.y='id')
rownames(pheno_df) <- pheno_df$Sample_Name
pheno_df <- pheno_df[sample,]
pheno_df$Braak.NFT.stage <- as.numeric(pheno_df$Braak.NFT.stage) 
pheno_df$sex_male <- as.factor(pheno_df$sex_male)
pheno_df$Sample_Plate <- as.factor(pheno_df$Sample_Plate)
pheno_df$bank <- as.factor(pheno_df$bank)

cell.pro = read.table('CBMAP_all_sample_methylation_ctp_with_clr_deconvolution.txt',header=T)
cell_aligned <- cell.pro[pheno_df$barcode_id, , drop = FALSE]
pheno_df <- cbind(pheno_df, cell_aligned)
pheno_df <- pheno_df %>% select(Braak.NFT.stage, age, sex_male, bank, PMD, RIN, Exc, Inh, NonN_Astro_FGF3R, NonN_Endo, NonN_Micro, NonN_Oligo_MBP, NonN_OPC, Sample_Plate)  
apply(pheno_df, 2, function(x){sum(is.na(x))})
# Braak.NFT.stage              age         sex_male             bank              PMD 
#               0                0                0                0                3 
# RIN              Exc              Inh NonN_Astro_FGF3R        NonN_Endo 
#   0               10               10               10               10 
# NonN_Micro   NonN_Oligo_MBP         NonN_OPC     Sample_Plate 
#         10               10               10                0

pheno_df[] <- lapply(pheno_df, function(x)
  if (is.numeric(x)) replace(x, is.na(x), median(x, na.rm = TRUE)) else x
)
identical(rownames(pheno_df), rownames(cg))

run_linear_regression <- function(methylation_data, covariates) {
  data <- data.frame(methylation_level = methylation_data, covariates)
  model <- lm(methylation_level ~ Braak.NFT.stage + sex_male + age + bank + PMD + RIN + Exc + Inh + NonN_Astro_FGF3R + NonN_Endo + NonN_Micro + NonN_Oligo_MBP + Sample_Plate, data = data)
  summary_model <- summary(model)
  return(summary_model$coefficients[2,])
}
analyze_methylation <- function(cg, pheno_df) {
  results <- apply(cg, 2, function(x) {
    methylation_data <- x  
    covariates <- pheno_df %>% select(Braak.NFT.stage, age, sex_male, bank, PMD, RIN, Exc, Inh, NonN_Astro_FGF3R, NonN_Endo, NonN_Micro, NonN_Oligo_MBP, Sample_Plate)  
    result <- run_linear_regression(methylation_data, covariates)  
    return(result)
  })
  return(results)
}

results = analyze_methylation(cg, pheno_df)
save(results, file='CBMAP_meth_braak.RData')
