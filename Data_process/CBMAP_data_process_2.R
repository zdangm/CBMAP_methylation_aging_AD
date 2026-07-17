library(data.table)
library(stringr)
library(ChAMP)
library(readxl)
suppressMessages(library(bigmelon))
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
suppressMessages(library(rhdf5))
suppressMessages(library(tidyverse))
suppressMessages(library(dplyr))
suppressMessages(library(wateRmelon))
suppressMessages(library(sva))
suppressMessages(library(corpcor))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(ewastools))

setwd('/methylation/data/CBMAP/dedup_mean_QC')
## Step 1: read idat files
gfile <- iadd2('/methylation/data/CBMAP/DNAm_processed/mirror_raw_data_DNAm/zk', gds = 'zk_DNAm_mean.gds', chunksize = 150)
write.table(paste("dim_rawData:", paste(dim(betas(gfile)[,]), collapse="*")), file="zk_mean_log.DNAm_QC", row.names=F, ,quote=F, col.names=F, append=F)

## Step 2: read phenotype and batch information
pdat <- read.csv('/methylation/data/CBMAP/DNAm_processed/phenotype.csv', header=T)
rownames(pdat) = pdat$barcode_id
pdat <- pdat[,c('barcode_id','sex_male','age','Sample_Plate','bank')]

## Step 3: filter
pfilter(
  gfile,
  perCount = 5, # rm probes with beadcount<3 in >5% (perCount) samples
  pnthresh = 0.01, # p value threshold
  perc = 1, # rm samples with detection p-value>0.01 (pnthresh) in >1% (perc) probes
  pthresh = 5 # rm probes with detection p-value>0.01 (pnthresh) in >5% (pthresh) samples
)
closefn.gds(gfile)
write.table(paste("dim_rmLowQual:", paste(dim(betas(gfile)[,]), collapse="*")), file="zk_mean_log.DNAm_QC", row.names=F, quote=F, col.names=F, append=T)

Probe_ID_prefix <- sub("_.*", "", fData(gfile)$Probe_ID)
nCpG_rmSNPprobe <- sum(Probe_ID_prefix %in% snp_probes) 
nCpG_rmxreact <- sum(Probe_ID_prefix %in% cross_reactive) 
nCpG_rm5bpSNP <- sum(Probe_ID_prefix %in% snp_1KG_SGCHS_SBE_sub) 
probe_rm <- fData(gfile)$Probe_ID[Probe_ID_prefix %in% probe_rm]  
gfile_filter(probe_rm, c())
write.table(paste("nCpG_rmSNPprobe:", nCpG_rmSNPprobe), file="zk_mean_log.DNAm_QC", row.names=F, quote=F, col.names=F, append=T)
write.table(paste("nCpG_rmxreact:", nCpG_rmxreact), file="zk_mean_log.DNAm_QC", row.names=F, quote=F, col.names=F, append=T)
write.table(paste("nCpG_rm5bpSNP:", nCpG_rm5bpSNP), file="zk_mean_log.DNAm_QC", row.names=F, quote=F, col.names=F, append=T)
write.table(paste("dim_rmSNP_xreact:", paste(dim(betas(gfile)[,]), collapse="*")), file="zk_mean_log.DNAm_QC", row.names=F, quote=F, col.names=F, append=T)


sample_ID = read.gdsn(index.gdsn(gfile,"pData/barcode"))
sample_ID = str_split_fixed(sample_ID,'/',n=2)[,2]
write.gdsn(index.gdsn(gfile, "pData/barcode"), sample_ID)
sample_rm = sample_ID[which(! sample_ID %chin% pdat$barcode_id)]
gfile_filter(c(), sample_rm)
write.table(paste("dim_rmFinalmissing:", paste(dim(betas(gfile)[,]), collapse="*")), file="zk_mean_log.DNAm_QC", row.names=F, quote=F, col.names=F, append=T)


## Step 4: predict sex
idatlist <- sub("_Red.idat", "", paste0('/methylation/data/CBMAP/DNAm_processed/mirror_raw_data_DNAm/zk/', list.files(path='/methylation/data/CBMAP/DNAm_processed/mirror_raw_data_DNAm/zk', pattern="*_Red.idat", recursive=T)))
XY <- t(sapply(
  idatlist, 
  function(infile){
    meth <- read_idats(infile, quiet=T)
    meth$manifest$chr <- str_split_fixed(meth$manifest$chr, 'chr', n=2)[,2]
    meth <- correct_dye_bias(meth)
    unlist(check_sex(meth))
  }
))
rownames(XY) <- basename(rownames(XY))
pdat2 <- cbind(pdat[rownames(XY),], XY)
pdat2$mSex <- toupper(predict_sex(pdat2$X, pdat2$Y, which(pdat2$sex_male==1), which(pdat2$sex_male==0)))
pdat2$mSex <- ifelse(pdat2$mSex=='M', 1, pdat2$mSex)
pdat2$mSex <- ifelse(pdat2$mSex=='F', 0, pdat2$mSex)
male <- which(pdat2$sex_male==1)
female <- which(pdat2$sex_male==0)
cutX = median(outer(pdat2$X[male], pdat2$X[female], "+"))/2 # Hodges-Lehmann estimators to separate male and female clusters
cutY = median(outer(pdat2$Y[male], pdat2$Y[female], "+"))/2
pdf("zk_sexMatching.correctDye.pdf", width=5, height=5)
tmp1 = pdat2 %>% dplyr::filter(mSex==sex_male)
tmp2 = pdat2 %>% dplyr::filter(is.na(mSex))
tmp3 = pdat2 %>% dplyr::filter(mSex!=sex_male)
plot(Y ~ X, data=tmp1, pch=ifelse(tmp1$sex_male==0,1,4), asp=1, xlab="Normalized X chromosome intensities", ylab="Normalized Y chromosome intensities")
abline(h=cutY, lty=3, col="grey")
abline(v=cutX, lty=3, col="grey")
points(Y ~ X, data=tmp2, pch=ifelse(tmp2$sex_male==0,1,4), col="grey")
points(Y ~ X, data=tmp3, pch=ifelse(tmp3$sex_male==0,1,4), col=2)
legend("topright", pch=c(1,4), legend=c("female","male"))
dev.off()
write.table(pdat2, file="zk_sexMatching.correctDye.txt", sep="\t", row.names=F, quote=F, col.names=T)

sp_rm <- rownames(pdat2 %>% dplyr::filter(mSex!=sex_male))
sp_rm = unique(c(sp_rm,sample_rm))
gfile_filter(probe_rm, sp_rm)
write.table(paste("dim_rmSexMismatch:", paste(dim(betas(gfile)[,]), collapse="*")), file="zk_mean_log.DNAm_QC", row.names=F, quote=F, col.names=F, append=T)

## Step 5: normalization
normbeta <- read.gdsn(betas(gfile))
rownames(normbeta) <- fData(gfile)$Probe_ID
colnames(normbeta) <- basename(pData(gfile)$barcode)
bcfail <- read.gdsn(index.gdsn(gfile, 'NBeads')) < 3
pfail <- read.gdsn(index.gdsn(gfile, 'pvals')) > 0.01
normbeta[bcfail | pfail] <- NA # set NA for low quality signals
normbeta = champ.impute(beta=normbeta, pd=pdat)

# Deduplication
normbeta = normbeta$beta
probe_prefix <- sub("_.*", "", rownames(normbeta))
dt <- as.data.table(normbeta)
dt[, Probe := probe_prefix]
dt_sizes <- dt[, .N, by = Probe]
probes_1 <- dt_sizes[N == 1, Probe]
probes_2 <- dt_sizes[N == 2, Probe]
probes_3_plus <- dt_sizes[N >= 3, Probe]

calc_mean_no_outliers_fast <- function(x) {
  x <- na.omit(x)
  n <- length(x)
  if (n <= 2) return(mean(x)) 
  
  Q1 <- quantile(x, 0.25) 
  Q3 <- quantile(x, 0.75)
  IQR_val <- Q3 - Q1
  
  lower_bound <- Q1 - 0.5 * IQR_val
  upper_bound <- Q3 + 0.5 * IQR_val
  
  x_filtered <- x[x >= lower_bound & x <= upper_bound]
  if (length(x_filtered) == 0) return(mean(x))
  
  return(mean(x_filtered))
}
res_1 <- dt[Probe %in% probes_1]
if (length(probes_2) > 0) {
  res_2 <- dt[Probe %in% probes_2, lapply(.SD, mean, na.rm = TRUE), by = Probe]
} else {
  res_2 <- data.table()
}
if (length(probes_3_plus) > 0) {
  res_3 <- dt[Probe %in% probes_3_plus, lapply(.SD, calc_mean_no_outliers_fast), by = Probe]
} else {
  res_3 <- data.table()
}
dt_dedup <- rbindlist(list(res_1, res_2, res_3), use.names = TRUE, fill = TRUE)
probe_names <- dt_dedup$Probe
dt_numeric <- dt_dedup[, !"Probe"]
normbeta_dedup <- as.matrix(dt_numeric)
rownames(normbeta_dedup) <- probe_names
rownames(normbeta_dedup) <- dt_dedup$Probe
normbeta = champ.norm(beta=normbeta_dedup,
                      method="BMIQ",
                      plotBMIQ=FALSE,
                      arraytype="EPIC")

## Step 6: correct batch effect
normbeta[normbeta < 0.01] <- 0.01
normbeta[normbeta > 0.99] <- 0.99
Mdat <- Beta2M(normbeta)
#rownames(Mdat) <- fData(gfile)$Probe_ID
rm(normbeta); gc()

# run combat, protecting age, sex
sample = intersect(colnames(Mdat),rownames(pdat))
modcombat = model.matrix(~ age + sex_male, data=pdat[sample,])
combat = ComBat(dat=Mdat[,sample], batch=pdat[sample, "Sample_Plate"], mod=modcombat)

h5write(combat, "zk_BMIQ_Mdat_combat.h5", "Mdat")
h5write(rownames(combat), "zk_BMIQ_Mdat_combat.h5", "rownames")
h5write(colnames(combat), "zk_BMIQ_Mdat_combat.h5", "colnames")

closefn.gds(gfile)
