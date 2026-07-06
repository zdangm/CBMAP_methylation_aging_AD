library(data.table)
library(EpiDISH)
library(readxl)
library(bigmelon)
library(rhdf5)
data(centDHSbloodDMC.m)

## 1. estimate celltype proportion
cg = fread('/methylation/data/GSE40279/GSE40279_average_beta.txt.gz',header=T)
beta = cg[,2:ncol(cg)]
beta = as.matrix(beta)
rownames(beta) = cg$ID_REF
BloodFrac.m <- epidish(beta.m = beta, ref.m = centDHSbloodDMC.m, method = "RPC")$estF
write.table(BloodFrac.m,file='/methylation/data/GSE40279/processed/celltype.txt',row.names=T,col.names=T,quote=F,sep='\t')
boxplot(BloodFrac.m)

library(compositions)
ctp = read.table('/methylation/data/GSE40279/processed/celltype.txt',sep='\t',header=T)
ctp<-as.matrix(ctp)
ctp[which(ctp <= 0)] <- 1e-3
ctp<-clr(ctp)
write.table(ctp,"/methylation/data/GSE40279/processed/celltype_with_clr.txt",row.names = T,col.names = T,sep = "\t")



## 2. phenotype infomation
pheno_info <- read_excel('/methylation/data/GSE40279/raw_data/sample_info.xlsx')
celltype <- read.table('/methylation/data/GSE40279/processed/celltype_with_clr.txt')
pheno_df <- merge(pheno_info, celltype, by.x='sampleID', by.y='row.names')
save(pheno_df, file='/methylation/data/GSE40279/processed/pheno_info.RData')


## 3. methylation matrix
# inverse normalization transform
cg = fread('/methylation/data/GSE40279/raw_data/GSE40279_average_beta.txt.gz',header=T)
beta = cg[,2:ncol(cg)]
beta = as.matrix(beta)
rownames(beta) = cg$ID_REF
beta[beta < 0.01] <- 0.01
beta[beta > 0.99] <- 0.99
Mdat <- Beta2M(beta)
InvNorm <- function(x) {
  return(qnorm((rank(x, na.last="keep")-3/8)/(sum(!is.na(x))+1/4)))
}
inverse_DNAm = apply(Mdat, 1, InvNorm)
h5write(inverse_DNAm, '/methylation/data/GSE40279/processed/final_DNAm_invMdat.h5', "Mdat")
h5write(rownames(inverse_DNAm), '/methylation/data/GSE40279/processed/final_DNAm_invMdat.h5', "rownames")
h5write(colnames(inverse_DNAm), '/methylation/data/GSE40279/processed/final_DNAm_invMdat.h5', "colnames")
