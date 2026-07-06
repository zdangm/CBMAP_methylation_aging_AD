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

setwd('/methylation/data/ROSMAP_methylaition')
synapse_meta = fread('/ROSMAP/Data/Epigenetics/Epigenetics (DNA methylation array)/IDAT Files/SYNAPSE_METADATA_MANIFEST.tsv')[,c('path','specimenID')]
cov = fread('/ROSMAP/Data/Epigenetics/ROSMAP_arrayMethylation_covariates.tsv')
clinical_info = fread('/ROSMAP/Data/Metadata/ROSMAP_clinical.csv')
clinical_info$id = paste0(clinical_info$Study, '', clinical_info$projid)
idkey = read.csv('/ROSMAP/ROSMAP_IDkey.csv', header=T)

##phenotype information
mer1=merge(idkey[,c('projid','mwas_id')], clinical_info[,c('projid','msex','age_death')], by='projid')
pheno_info = merge(mer1, cov, by.x='mwas_id', by.y='Sample')
pheno_info$projid <- as.character(str_pad(pheno_info$projid, width = 8, pad = "0", side = "left"))
hml_info <- read_excel('/ROSMAP/Processed/metadata/dataset_1495_cross-sectional_02-16-2025.xlsx')
pheno_info$age_death <- unlist(hml_info[match(pheno_info$projid, hml_info$projid),'age_death'])
rm(mer1)
pheno_info$barcode_id = paste0(pheno_info$Sentrix_ID,'_',pheno_info$Sentrix_Position)
pheno_info = distinct(pheno_info)
rownames(pheno_info) = pheno_info$barcode_id
#5822038012_R02C01  TBI-AUTO73203-PT-35OA  has missing sex and age_death
pheno_info = pheno_info[-which(rownames(pheno_info)=='5822038012_R02C01'),]
write.table(pheno_info, file='/methylation/data/ROSMAP_methylaition/phenotype.txt', row.names=T, col.names=T, quote=F)


## Step 1: read idat files
indir <- "/data/ROSMAP_methylaition/idat_data/"
gfile <- iadd2(indir, gds = 'DNAm.gds', chunksize = 100)
write.table(paste("dim_rawData:", paste(dim(betas(gfile)[,]), collapse="*")), file="log.DNAm_QC", row.names=F, ,quote=F, col.names=F, append=F)

## Step 2: read phenotype and batch information
pdat <- read.table("/methylation/data/ROSMAP_methylaition/phenotype.txt", header=T)
pdat <- pdat[,c('barcode_id','msex','age_death','Sample_Plate')]

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
pdat <- read.table("/methylation/data/ROSMAP_methylaition/phenotype.txt", header=T)
pd = pdat[,c('barcode_id', 'projid', 'Sentrix_ID', 'Sentrix_Position', 'Sample_Plate', 'Sample_Well', 'Sample_Group')]
colnames(pd) = c('Sample_Name', 'Project', 'Pool_ID', 'Array', 'Sample_Plate', 'Sample_Well', 'Sample_Group')
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
#h5write(Mdat, "BMIQ_Mdat.h5", "Mdat")
#h5write(rownames(Mdat), "BMIQ_Mdat.h5", "rownames")
#h5write(colnames(Mdat), "BMIQ_Mdat.h5", "colnames") 


# run combat, protecting age, sex
pdat = na.omit(pdat)
sample = intersect(colnames(Mdat),rownames(pdat))
pdat3 = pdat[sample,c("projid","msex","age_death","Sentrix_ID", "Sentrix_Position" ,"Sample_Plate", "Sample_Well" ,"Sample_Group","batch")]
Mdat = Mdat[,sample]
modcombat = model.matrix(~ age_death + msex, data=pdat3)
combat = ComBat(dat=Mdat, batch=pdat3$Sample_Plate, mod=modcombat)

h5write(combat, "BMIQ_Mdat_combat.h5", "Mdat")
h5write(rownames(combat), "BMIQ_Mdat_combat.h5", "rownames")
h5write(colnames(combat), "BMIQ_Mdat_combat.h5", "colnames")

closefn.gds(gfile)


##### inverse normalization transform  #####
rosmap_cg = h5read('/methylation/data/ROSMAP_methylaition/BMIQ_Mdat_combat.h5', 'Mdat')
colnames(rosmap_cg) = h5read('/methylation/data/ROSMAP_methylaition/BMIQ_Mdat_combat.h5', 'colnames')
rownames(rosmap_cg) = h5read('/methylation/data/ROSMAP_methylaition/BMIQ_Mdat_combat.h5', 'rownames')
InvNorm <- function(x) {
  return(qnorm((rank(x, na.last="keep")-3/8)/(sum(!is.na(x))+1/4)))
}
inverse_DNAm = apply(rosmap_cg, 1, InvNorm)

pheno_info = read.table(file='/methylation/data/ROSMAP_methylaition/phenotype.txt', header=T)
rownames(inverse_DNAm) = pheno_info[match(rownames(inverse_DNAm), pheno_info$barcode_id), 'mwas_id']

h5write(inverse_DNAm, '/methylation/data/ROSMAP_methylaition/final_DNAm_invMdat.h5', "Mdat")
h5write(rownames(inverse_DNAm), '/methylation/data/ROSMAP_methylaition/final_DNAm_invMdat.h5', "rownames")
h5write(colnames(inverse_DNAm), '/methylation/data/ROSMAP_methylaition/final_DNAm_invMdat.h5', "colnames")

