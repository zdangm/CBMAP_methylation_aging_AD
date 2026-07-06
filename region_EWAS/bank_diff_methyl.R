library(rhdf5)
library(data.table)
library(dplyr)
library(stringr)

# methylation matrix
cg <- h5read("final_DNAm_invMdat.h5",'Mdat')
colnames(cg) = h5read("final_DNAm_invMdat.h5",'colnames')
rownames(cg) = h5read("final_DNAm_invMdat.h5",'rownames')

# phenotype
load('CBMAP_pheno.RData')
cell.pro = read.table('CBMAP_all_sample_methylation_ctp_with_clr_deconvolution.txt',header=T)
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
identical(rownames(pheno_df), rownames(cg))

###pairwise_compare
pairwise_results <- apply(cg, 2, function(x){
  model <- aov(x ~ bank + age + sex_male + bank + PMD + RIN + Exc + Inh + NonN_Astro_FGF3R + NonN_Endo + NonN_Micro + NonN_Oligo_MBP + Sample_Plate, data=pheno_df)
  tukey <- TukeyHSD(model, which = "bank")
  results <- tukey$bank[, c('diff',"p adj")]
  return(results)
})
pairwise_results <- t(pairwise_results)
colnames(pairwise_results) <- c('diff_pumc-csu','diff_zju-csu','diff_zju-pumc','p_pumc-csu','p_zju-csu','p_zju-pumc')
save(pairwise_results,file='bank_pairwise_results.RData')
