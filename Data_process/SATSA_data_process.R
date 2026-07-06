library(data.table)
library(stringr)
library('ChAMP')
library(readxl)
library(dplyr)
suppressMessages(library(bigmelon))
#suppressMessages(library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19))
suppressMessages(library(IlluminaHumanMethylation450kanno.ilmn12.hg19))
suppressMessages(library(rhdf5))
suppressMessages(library(tidyverse))
suppressMessages(library(dplyr))
suppressMessages(library(wateRmelon))
suppressMessages(library(sva))
suppressMessages(library(corpcor))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(ewastools))

setwd('/methylation/data/SATSA')

## Step 1: read idat files
indir <- "/EUR_PBMC_methylation/"
gfile <- iadd2(indir, gds = 'DNAm.gds', chunksize = 100)
write.table(paste("dim_rawData:", paste(dim(betas(gfile)[,]), collapse="*")), file="log.DNAm_QC", row.names=F, ,quote=F, col.names=F, append=F)

## Step 2: read phenotype and batch information
pdat <- read_excel("/EUR_PBMC_methylation/sample_info.xlsx")
pdat <- as.data.frame(pdat)
pdat$Sentrix_ID <- str_split_fixed(pdat$`Assay Name`,'_',n=3)[,1]
pdat$Sentrix_Position <- str_split_fixed(pdat$`Assay Name`,'_',n=3)[,2]
pdat$msex = ifelse(pdat$sex == 'male',1,0)
pdat <- pdat[,c(1,3,4,10,13:16)]
colnames(pdat) = c('Source_Name','individualID','age','twin_pair','Assay_Name','Sentrix_ID', 'Sentrix_Position', 'msex')
pdat = distinct(pdat, Sentrix_ID, Sentrix_Position, .keep_all = T)
pdat = pdat[,c('Source_Name','individualID','age','twin_pair','Sentrix_ID', 'Sentrix_Position', 'msex')]
rownames(pdat) = paste0(pdat$Sentrix_ID,'_',pdat$Sentrix_Position)
write.table(pdat,'pheno_df.txt',col.names = T,row.names = T,quote=F,sep='\t')

## Step 3: filter
pfilter(
  gfile,
  perCount = 5, # rm probes with beadcount<3 in >5% (perCount) samples
  pnthresh = 0.01, # p value threshold
  perc = 1, # rm samples with detection p-value>0.01 (pnthresh) in >1% (perc) probes
  pthresh = 5 # rm probes with detection p-value>0.01 (pnthresh) in >5% (pthresh) samples
)
write.table(paste("dim_rmLowQual:", paste(dim(betas(gfile)[,]), collapse="*")), file="log.DNAm_QC", row.names=F, quote=F, col.names=F, append=T)


Probe_ID <- fData(gfile)$Probe_ID
nCpG_rmSNPprobe <- sum(Probe_ID %in% snp_probes)
nCpG_rmxreact <- sum(Probe_ID %in% cross_reactive)
nCpG_rm5bpSNP <- sum(Probe_ID %in% snp_1KG_SGCHS_SBE_sub)
gfile_filter(probe_rm, c())
write.table(paste("nCpG_rmSNPprobe:", nCpG_rmSNPprobe), file="log.DNAm_QC", row.names=F, quote=F, col.names=F, append=T)
write.table(paste("nCpG_rmxreact:", nCpG_rmxreact), file="log.DNAm_QC", row.names=F, quote=F, col.names=F, append=T)
write.table(paste("nCpG_rm5bpSNP:", nCpG_rm5bpSNP), file="log.DNAm_QC", row.names=F, quote=F, col.names=F, append=T)
write.table(paste("dim_rmSNP_xreact:", paste(dim(betas(gfile)[,]), collapse="*")), file="log.DNAm_QC", row.names=F, quote=F, col.names=F, append=T)


## Step 4: predict sex
idatlist <- sub("_Red.idat", "", paste0(indir, list.files(path=indir, pattern="*_Red.idat", recursive=T)))
XY <- t(sapply(
  idatlist, 
  function(infile){
    meth <- read_idats(infile, quiet=T)
    meth <- correct_dye_bias(meth)
    unlist(check_sex(meth))
  }
))
rownames(XY) <- basename(rownames(XY))
pdat2 <- cbind(pdat[rownames(XY),], XY)
pdat2$mSex <- toupper(predict_sex(pdat2$X, pdat2$Y, which(pdat2$msex==1), which(pdat2$msex==0)))
pdat2$mSex <- ifelse(pdat2$mSex=='M', 1, pdat2$mSex)
pdat2$mSex <- ifelse(pdat2$mSex=='F', 0, pdat2$mSex)
male <- which(pdat2$msex==1)
female <- which(pdat2$msex==0)
cutX = median(outer(pdat2$X[male], pdat2$X[female], "+"))/2 # Hodges-Lehmann estimators to separate male and female clusters
cutY = median(outer(pdat2$Y[male], pdat2$Y[female], "+"))/2
pdf("sexMatching.correctDye.pdf", width=5, height=5)
tmp1 = pdat2 %>% dplyr::filter(mSex==msex)
tmp2 = pdat2 %>% dplyr::filter(is.na(mSex))
tmp3 = pdat2 %>% dplyr::filter(mSex!=msex)
plot(Y ~ X, data=tmp1, pch=ifelse(tmp1$msex==0,1,4), asp=1, xlab="Normalized X chromosome intensities", ylab="Normalized Y chromosome intensities")
abline(h=cutY, lty=3, col="grey")
abline(v=cutX, lty=3, col="grey")
points(Y ~ X, data=tmp2, pch=ifelse(tmp2$msex==0,1,4), col="grey")
points(Y ~ X, data=tmp3, pch=ifelse(tmp3$msex==0,1,4), col=2)
legend("topright", pch=c(1,4), legend=c("female","male"))
dev.off()
write.table(pdat2, file="sexMatching.correctDye.txt", sep="\t", row.names=F, quote=F, col.names=T)

sp_rm <- rownames(pdat2 %>% dplyr::filter(mSex!=msex))
gfile_filter(c(), sp_rm)
write.table(paste("dim_rmSexMismatch:", paste(dim(betas(gfile)[,]), collapse="*")), file="log.DNAm_QC", row.names=F, quote=F, col.names=F, append=T)


## Step 5: normalization
normbeta <- read.gdsn(betas(gfile))
rownames(normbeta) <- fData(gfile)$Probe_ID
colnames(normbeta) <- basename(pData(gfile)$barcode)
bcfail <- read.gdsn(index.gdsn(gfile, 'NBeads')) < 3
pfail <- read.gdsn(index.gdsn(gfile, 'pvals')) > 0.01
normbeta[bcfail | pfail] <- NA # set NA for low quality signals
pdat <- read.table("/methylation/data/SATSA/pheno_df.txt", header=T,sep='\t')
pd = pdat[colnames(normbeta),]
normbeta = champ.impute(beta=normbeta, pd=pd)
normbeta = champ.norm(beta=normbeta$beta,
                      method="BMIQ",
                      plotBMIQ=FALSE,
                      arraytype="450k",
                      cores=1)
h5write(normbeta, "BMIQ_normbeta.h5", "normbeta")
h5write(rownames(normbeta), "BMIQ_normbeta.h5", "rownames")
h5write(colnames(normbeta), "BMIQ_normbeta.h5", "colnames") 

## Step 6: correct batch effect
normbeta[normbeta < 0.01] <- 0.01
normbeta[normbeta > 0.99] <- 0.99
Mdat <- Beta2M(normbeta)
#rownames(Mdat) <- fData(gfile)$Probe_ID
rm(normbeta); gc()
h5write(Mdat, "BMIQ_Mdat.h5", "Mdat")
h5write(rownames(Mdat), "BMIQ_Mdat.h5", "rownames")
h5write(colnames(Mdat), "BMIQ_Mdat.h5", "colnames") 

closefn.gds(gfile)


##### inverse normalization transform  #####
cg = h5read('/methylation/data/SATSA/BMIQ_Mdat.h5', 'Mdat')
colnames(cg) = h5read('/methylation/data/SATSA/BMIQ_Mdat.h5', 'colnames')
rownames(cg) = h5read('/methylation/data/SATSA/BMIQ_Mdat.h5', 'rownames')
InvNorm <- function(x) {
  return(qnorm((rank(x, na.last="keep")-3/8)/(sum(!is.na(x))+1/4)))
}
inverse_DNAm = apply(cg, 1, InvNorm)

pheno_info = read.table(file='/methylation/data/SATSA/pheno_df.txt', header=T, sep='\t')

h5write(inverse_DNAm, 'final_DNAm_invMdat.h5', "Mdat")
h5write(rownames(inverse_DNAm), 'final_DNAm_invMdat.h5', "rownames")
h5write(colnames(inverse_DNAm), 'final_DNAm_invMdat.h5', "colnames")


############# estimate cell proportion ##########
library(data.table)
library(EpiDISH)
library(readxl)
library(bigmelon)
library(rhdf5)
data(centDHSbloodDMC.m)

## 1. estimate celltype proportion
cg <- h5read("/methylation/data/SATSA/BMIQ_normbeta.h5",'normbeta')
colnames(cg) = h5read("/methylation/data/SATSA/BMIQ_normbeta.h5",'colnames')
rownames(cg) = h5read("/methylation/data/SATSA/BMIQ_normbeta.h5",'rownames')
BloodFrac.m <- epidish(beta.m = cg, ref.m = centDHSbloodDMC.m, method = "RPC")$estF
write.table(BloodFrac.m,file='/methylation/data/SATSA/celltype.txt',row.names=T,col.names=T,quote=F,sep='\t')
boxplot(BloodFrac.m)
#################### Perform clr-transform
library(compositions)
ctp = read.table('/methylation/data/SATSA/celltype.txt',sep='\t',header=T)
ctp<-as.matrix(ctp)
ctp[which(ctp <= 0)] <- 1e-3
ctp<-clr(ctp)
write.table(ctp,"/methylation/data/SATSA/celltype_with_clr.txt",row.names = T,col.names = T,sep = "\t")
