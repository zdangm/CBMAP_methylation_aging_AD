library(data.table)
library(stringr)
library('ChAMP')
library(readxl)
suppressMessages(library(bigmelon))
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

library(GEOquery)
library(stringr)

GSE_data = getGEO(filename = '/methylation/data/GSE111629/GSE111629_series_matrix.txt.gz',getGPL = F)
sample_info = pData(GSE_data)
sample_info = sample_info[,c("geo_accession","source_name_ch1","age:ch1","disease state:ch1","ethnicity:ch1","gender:ch1","tissue:ch1" )]
colnames(sample_info) = c("geo_accession","source_name","age","disease","ethnicity","gender","tissue")
sample_info$source_name = str_split_fixed(sample_info$source_name,'X',n=2)[,2]
sample_info$sample_name = paste0('raw_data/',sample_info$geo_accession,'_',sample_info$source_name)
sample_info$gender = ifelse(sample_info$gender == 'Male',1,
                            ifelse(sample_info$gender == 'Female',0,NA))
sample_info$age = as.integer(sample_info$age)
write.table(sample_info,file='/methylation/data/GSE111629/sample_info.txt',col.names = T,quote = F,row.names = F,sep='\t')


setwd('/methylation/data/GSE111629')
## Step 1: read idat files
indir <- "/methylation/data/GSE111629/raw_data/"
gfile <- iadd2('/methylation/data/GSE111629/raw_data', gds = 'GSE111629_DNAm.gds', chunksize = 100)
write.table(paste("dim_rawData:", paste(dim(betas(gfile)[,]), collapse="*")), file="GSE111629.DNAm_QC", row.names=F, ,quote=F, col.names=F, append=F)

## Step 2: read phenotype and batch information
pdat <- as.data.frame(fread("/methylation/data/GSE111629/sample_info.txt", header=T))
rownames(pdat) = paste0(pdat$geo_accession,'_',pdat$source_name)

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
pdat2$mSex <- toupper(predict_sex(pdat2$X, pdat2$Y, which(pdat2$gender==1), which(pdat2$gender==0)))
pdat2$mSex <- ifelse(pdat2$mSex=='M', 1, pdat2$mSex)
pdat2$mSex <- ifelse(pdat2$mSex=='F', 0, pdat2$mSex)
male <- which(pdat2$gender==1)
female <- which(pdat2$gender==0)
cutX = median(outer(pdat2$X[male], pdat2$X[female], "+"))/2 # Hodges-Lehmann estimators to separate male and female clusters
cutY = median(outer(pdat2$Y[male], pdat2$Y[female], "+"))/2
pdf("sexMatching.correctDye.pdf", width=5, height=5)
tmp1 = pdat2 %>% dplyr::filter(mSex==gender)
tmp2 = pdat2 %>% dplyr::filter(is.na(mSex))
tmp3 = pdat2 %>% dplyr::filter(mSex!=gender)
plot(Y ~ X, data=tmp1, pch=ifelse(tmp1$gender==0,1,4), asp=1, xlab="Normalized X chromosome intensities", ylab="Normalized Y chromosome intensities")
abline(h=cutY, lty=3, col="grey")
abline(v=cutX, lty=3, col="grey")
points(Y ~ X, data=tmp2, pch=ifelse(tmp2$gender==0,1,4), col="grey")
points(Y ~ X, data=tmp3, pch=ifelse(tmp3$gender==0,1,4), col=2)
legend("topright", pch=c(1,4), legend=c("female","male"))
dev.off()
write.table(pdat2, file="sexMatching.correctDye.txt", sep="\t", row.names=F, quote=F, col.names=T)

sp_rm <- rownames(pdat2 %>% dplyr::filter(mSex!=gender))
gfile_filter(c(), sp_rm)
write.table(paste("dim_rmSexMismatch:", paste(dim(betas(gfile)[,]), collapse="*")), file="log.DNAm_QC", row.names=F, quote=F, col.names=F, append=T)


## Step 5: normalization
normbeta <- read.gdsn(betas(gfile))
rownames(normbeta) <- fData(gfile)$Probe_ID
colnames(normbeta) <- basename(pData(gfile)$barcode)
bcfail <- read.gdsn(index.gdsn(gfile, 'NBeads')) < 3
pfail <- read.gdsn(index.gdsn(gfile, 'pvals')) > 0.01
normbeta[bcfail | pfail] <- NA # set NA for low quality signals
pd = pdat[,2:6]
normbeta = champ.impute(beta=normbeta, pd=pd)
normbeta = champ.norm(beta=normbeta$beta,
                      method="BMIQ",
                      plotBMIQ=FALSE,
                      arraytype="450k",
                      cores=1)
#h5write(normbeta, "BMIQ_normbeta.h5", "normbeta")
#h5write(rownames(normbeta), "BMIQ_normbeta.h5", "rownames")
#h5write(colnames(normbeta), "BMIQ_normbeta.h5", "colnames") 

## Step 6: correct batch effect
normbeta[normbeta < 0.01] <- 0.01
normbeta[normbeta > 0.99] <- 0.99
Mdat <- Beta2M(normbeta)
#rownames(Mdat) <- fData(gfile)$Probe_ID
rm(normbeta); gc()

# run combat, protecting age, sex
pdat = na.omit(pdat)
sample = intersect(colnames(Mdat),rownames(pdat))
pdat3 = pdat[sample,2:6]
pdat3 = factor(pdat3$disease)
Mdat = Mdat[,sample]
h5write(Mdat, "BMIQ_Mdat.h5", "Mdat")
h5write(rownames(Mdat), "BMIQ_Mdat.h5", "rownames")
h5write(colnames(Mdat), "BMIQ_Mdat.h5", "colnames") 
# no batch info
# modcombat = model.matrix(~ age + gender + disease, data=pdat3)
# combat = ComBat(dat=Mdat, batch=pdat3$Sample_Plate, mod=modcombat)
closefn.gds(gfile)


##### inverse normalization transform  #####
InvNorm <- function(x) {
  return(qnorm((rank(x, na.last="keep")-3/8)/(sum(!is.na(x))+1/4)))
}

inverse_DNAm = apply(Mdat, 1, InvNorm)

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
cg <- h5read("/methylation/data/GSE111629/BMIQ_Mdat.h5",'Mdat')
colnames(cg) = h5read("/methylation/data/GSE111629/BMIQ_Mdat.h5",'colnames')
rownames(cg) = h5read("/methylation/data/GSE111629/BMIQ_Mdat.h5",'rownames')
cg = M2Beta(cg)
BloodFrac.m <- epidish(beta.m = cg, ref.m = centDHSbloodDMC.m, method = "RPC")$estF
write.table(BloodFrac.m,file='/methylation/data/GSE111629/celltype.txt',row.names=T,col.names=T,quote=F,sep='\t')
boxplot(BloodFrac.m)
#################### Perform clr-transform
library(compositions)
ctp = read.table('/methylation/data/GSE111629/celltype.txt',sep='\t',header=T)
ctp<-as.matrix(ctp)
ctp[which(ctp <= 0)] <- 1e-3
ctp<-clr(ctp)
write.table(ctp,"/methylation/data/GSE111629/celltype_with_clr.txt",row.names = T,col.names = T,sep = "\t")
