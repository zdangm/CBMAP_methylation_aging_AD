#devtools::install_github("MorganLevineLab/calcPCBrainAge")
library(calcPCBrainAge)
library(rhdf5)
library(data.table)

setwd('/methylation/results/disease_meth_age/DNAm_clock/PCBrainAge/pumctest')
load('/methylation/results/disease_meth_age/AD_meth/CBMAP_pheno.RData')
source('/methylation/script/disease_meth_age/DNAm_clock/PCBrainAge/calcPCBrainAge.R')
# methylation matrix
Mdat <- h5read("/methylation/data/CBMAP/DNAm_processed/final_DNAm_Mdat.h5",'Mdat')
colnames(Mdat) = h5read("/methylation/data/CBMAP/DNAm_processed/final_DNAm_Mdat.h5",'colnames')
rownames(Mdat) = h5read("/methylation/data/CBMAP/DNAm_processed/final_DNAm_Mdat.h5",'rownames')
Bdat <- 2^Mdat / (2^Mdat + 1)

pheno_df = pheno_df[which(pheno_df$bank == 'pumc'),]
sample = intersect(colnames(Bdat), rownames(pheno_df)) 
Bdat <- Bdat[,sample]
Bdat <- t(Bdat)
Bdat <- as.data.frame(Bdat)
pheno_df <- pheno_df[sample,]
identical(rownames(pheno_df), rownames(Bdat))
pheno_df = pheno_df[,c('Sample_Name','age')]

res = calcPCBrainAge(DNAm = Bdat, pheno = pheno_df, CpGImputation = imputeMissingBrainCpGs)
write.table(res,file='PCBrainAgePred.txt',row.names=F,col.names=T,quote=F)
