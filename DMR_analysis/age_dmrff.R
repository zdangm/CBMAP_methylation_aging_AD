library(stringr)
library(readxl)
library(dplyr)
library(data.table)
library(rhdf5)
library(dmrff)

load("CBMAP_meth_age_no_dementia.RData")
# methylation matrix
cg <- h5read("final_DNAm_invMdat.h5",'Mdat')
colnames(cg) = h5read("final_DNAm_invMdat.h5",'colnames')
rownames(cg) = h5read("final_DNAm_invMdat.h5",'rownames')

#phenotype
load('CBMAP_pheno.RData')
cell.pro = read.table('CBMAP_all_sample_methylation_ctp_with_clr_deconvolution.txt',header=T)
cell_aligned <- cell.pro[pheno_df$barcode_id, , drop = FALSE]
pheno_df <- cbind(pheno_df, cell_aligned)

sample = intersect(rownames(cg), rownames(pheno_df)) # 1031 samples
cg <- cg[sample,]
pheno_df <- pheno_df[sample,]
pheno_df$sex_male <- as.factor(pheno_df$sex_male)
pheno_df$Sample_Plate <- as.factor(pheno_df$Sample_Plate)
pheno_df$bank <- as.factor(pheno_df$bank)
pheno_df[] <- lapply(pheno_df, function(x)
  if (is.numeric(x)) replace(x, is.na(x), median(x, na.rm = TRUE)) else x
)

# filter hc samples
all_sample_info <- as.data.frame(read.csv("CBMAP_sample_info.csv",header=T)) 
selected_sample = all_sample_info$id[which(all_sample_info$diag_dementia == 1)]
pheno_df = pheno_df[which(! pheno_df$Sample_Name %chin% selected_sample),]
cg = cg[pheno_df$Sample_Name,]
identical(rownames(pheno_df), rownames(cg))

anno_935k = read.csv('EPIC_annotation.csv',header=T)
anno_935k <- anno_935k[,c('Name','CHR','MAPINFO')]
result <- merge(res,anno_935k,by.x='cpg',by.y='Name')
cg <- t(cg)

###dmrff
age_dmr <- dmrff(estimate = result$Estimate,se=result$StdErr,p.value = as.vector(result$pvalue),methylation = cg,chr=result$CHR,pos=as.numeric(result$MAPINFO))
age_dmr_sig <- age_dmr[age_dmr$p.adjust < 0.05 & age_dmr$n >= 2, ]
save(age_dmr_sig,file = 'age_dmr_sig.RData')
