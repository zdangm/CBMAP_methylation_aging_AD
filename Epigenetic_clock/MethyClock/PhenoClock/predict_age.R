library(rhdf5)
library(data.table)

setwd('/methylation/results/disease_meth_age/DNAm_clock/PhenoClock/pumctest')
source('/methylation/script/disease_meth_age/DNAm_clock/PhenoClock/PhenoClock.R')
load('CBMAP_pheno.RData')
# methylation matrix
Mdat <- h5read("/methylation/data/CBMAP/DNAm_processed/final_DNAm_Mdat.h5",'Mdat')
colnames(Mdat) = h5read("/methylation/data/CBMAP/DNAm_processed/final_DNAm_Mdat.h5",'colnames')
rownames(Mdat) = h5read("/methylation/data/CBMAP/DNAm_processed/final_DNAm_Mdat.h5",'rownames')
Bdat <- 2^Mdat / (2^Mdat + 1)

pheno_df = pheno_df[which(pheno_df$bank == 'pumc'),]
sample = intersect(colnames(Bdat), rownames(pheno_df))
Bdat <- Bdat[,sample]
pheno_df <- pheno_df[sample,]
identical(rownames(pheno_df), colnames(Bdat))

dir <- '/methylation/script/disease_meth_age/DNAm_clock/PhenoClock/'
IDcol='Sample_Name'
Agecol='age'
pre_age <- PhenoClock(betas=Bdat, pheno=pheno_df, dir=dir, IDcol=IDcol, Agecol=Agecol)
