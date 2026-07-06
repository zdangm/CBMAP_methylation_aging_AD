
#==============================================================================
#
# Estimate CTP using Houseman 2012 method + sequencing data (Luo et al. 2022)
#
#==============================================================================

#------------------------------------------------------------------------------
# FUNCTIONS: READ IN FIRST

# Houseman 2012 deconvolution implemented in minfi (Aryee et al. 2014 Bioinformatics)
# projectCellType (an internal function within minfi estimateCellCounts.R)
# (https://rdrr.io/bioc/minfi/src/R/estimateCellCounts.R)
# this function preferred here as requiring just methylation beta matrix 
# is more flexible

# Inputs:
# Y: bulk methylation beta matrix
# coefCellType: reference, matrix of m rows (CpG probes) x n columns (cell-types)

projectCellType <- function(Y, coefCellType, contrastCellType = NULL,
                            nonnegative = TRUE, lessThanOne = TRUE) {
  if (is.null(contrastCellType)) {
    Xmat <- coefCellType
  } else {
    Xmat <- tcrossprod(coefCellType, contrastCellType)
  }
  
  nCol <- dim(Xmat)[2]
  if (nCol == 2) {
    Dmat <- crossprod(Xmat)
    mixCoef <- t(
      apply(Y, 2, function(x) solve(Dmat, crossprod(Xmat, x))))
    colnames(mixCoef) <- colnames(Xmat)
    return(mixCoef)
  } else {
    nSubj <- dim(Y)[2]
    
    mixCoef <- matrix(0, nSubj, nCol)
    rownames(mixCoef) <- colnames(Y)
    colnames(mixCoef) <- colnames(Xmat)
    
    if (nonnegative) {
      if (lessThanOne) {
        Amat <- cbind(rep(-1, nCol), diag(nCol))
        b0vec <- c(-1, rep(0, nCol))
      } else {
        Amat <- diag(nCol)
        b0vec <- rep(0, nCol)
      }
      for (i in seq_len(nSubj)) {
        obs <- which(!is.na(Y[,i]))
        Dmat <- crossprod(Xmat[obs,])
        mixCoef[i,] <- solve.QP(
          Dmat = Dmat,
          dvec = crossprod(Xmat[obs,], Y[obs,i]),
          Amat = Amat,
          bvec = b0vec)$sol
      }
    } else {
      for (i in seq_len(nSubj)) {
        obs <- which(!is.na(Y[,i]))
        Dmat <- crossprod(Xmat[obs,])
        mixCoef[i,] <- solve(Dmat, t(Xmat[obs,]) %*% Y[obs,i])
      }
    }
    mixCoef
  }
}

# Seiler-Vellame et al. 2022 biorXiv https://www.biorxiv.org/content/10.1101/2022.06.15.496235v1

# Inputs:
# applyIndex: getErrorPerSample() implemented within sapply loop to iterate through each matrix row
# predictedIN: output cell-type proportions
# coefDataIN: reference cell-type profile matrix
# betasBulkIN: bulk methylation matrix

getErrorPerSample = function(applyIndex,
                             predictedIN = counts,
                             coefDataIN = coefs,
                             betasBulkIN = coefdat){
  
  trueBulk = matrix(ncol = 1, nrow = nrow(coefDataIN), data = 0)
  
  RMSE = function(m, o){
    sqrt(mean((m - o)^2))
  }
  
  for (i in 1:ncol(coefDataIN)){
    
    trueBulk[,1] = trueBulk[,1] + coefDataIN[,i]*predictedIN[applyIndex,i]
  }
  
  betasBulkIN = t(apply(betasBulkIN, 1, function(x){x[is.na(x)] = 0; return(x)}))
  
  error = RMSE(trueBulk, betasBulkIN[,applyIndex])
  return(error)
}
#------------------------------------------------------------------------------

setwd("/data/projects/China_Brain_MultiOmics/methylation/results/")

library(minfi)
library(quadprog)
library(csSAM)
library(data.table)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(ggsci)
library(gdsfmt)
library('ChAMP',lib="/home/zhaiwl/R/x86_64-pc-linux-gnu-library/4.3")
library(compositions)
suppressMessages(library(rhdf5))


# Coefficients from sequencing data - download from: brain_CTP_deconv/deconvolution/data
ref_dir <- "/data/projects/China_Brain_MultiOmics/methylation/data/ctp_reference/"
coefs_seq_dir <- paste(ref_dir, "Luo2020_extremes_dmr_ilmn450kepic_aggto100bp_c10_cell7_split6040all_beta.rds", sep = "")
coefs_neun_dir <- paste(ref_dir, "dlpfc_450k_guintivano_rowftest_cell2.rds", sep = "")

pdat <- read.csv('/data/projects/China_Brain_MultiOmics/methylation/data/CBMAP/DNAm_processed/phenotype.csv', header=T)
rownames(pdat) = pdat$barcode_id

clinical<-read.csv("/data/shared_data/China_Brain_MultiOmics/sample_information/final/CBMAP_sample_info_1187_final_20241101.csv",fileEncoding="GBK",check.names = T)

pheno_info<-merge(pdat,clinical,by.x="Sample_Name",by.y="id",all.x=T)

#pheno_info<-as.data.table(pheno_info)


############# read meth data
meth_dir<-"/data/projects/China_Brain_MultiOmics/methylation/results/ctp_analysis/OSCA/CBMAP/CBMAP_norm_beta_OSCA_input_aut_mask_t.txt"
meth.tmp <- fread(meth_dir)
meth <- as.matrix(meth.tmp[,2:ncol(meth.tmp)])
rownames(meth) <- data.frame(meth.tmp)[,1]
meth <- meth[order(rownames(meth)),]

############ order pheno sample
pheno_info<-pheno_info[match(colnames(meth),pheno_info$barcode_id),]
unique(colnames(meth)==pheno_info$barcode_id)


# Reference probes, order, and overlap with methylation matrix
# Formatted as a list here, as it gives flexibility to deconvolve 
# one methylation beta matrix with multiple different reference
# panels for comparison
coefs_seq <- readRDS(coefs_seq_dir)
coefs_neun <- readRDS(coefs_neun_dir)
coefs.ls <- list(coefs_seq, coefs_neun)
coefs.ls <- lapply(coefs.ls, function(x) x[order(rownames(x)),])
coefs.ls <- lapply(coefs.ls, function(x) x[which(rownames(x) %in% rownames(meth)),])

#==============================================================================
# DECONVOLUTION
#==============================================================================

# Deconvolve CTP
ctp.ls <- lapply(coefs.ls, function(x) projectCellType(meth[which(rownames(meth) %in% (rownames(x))), ], x, lessThanOne = TRUE))

# Add errors
ctp_error.ls <- lapply(1:length(ctp.ls), function(y) data.frame(IID = rownames(ctp.ls[[y]]), ctp.ls[[y]], error =
                                                                  sapply(1:nrow(ctp.ls[[y]]), function(x) getErrorPerSample(x, 
                                                                                                                            ctp.ls[[y]], 
                                                                                                                            coefs.ls[[y]][which(rownames(coefs.ls[[y]]) %in% rownames(meth)),], 
                                                                                                                            meth[which(rownames(meth) %in% (rownames(coefs.ls[[y]]))), ]))))
names(ctp_error.ls) <- names(ctp.ls)

ctp<-as.data.frame(ctp.ls)

# Write output
saveRDS(ctp_error.ls,  "/data/projects/China_Brain_MultiOmics/methylation/results/ctp_analysis/CBMAP_all_sample_methylation_ctp_deconvolution.rds")

write.table(ctp,"/data/projects/China_Brain_MultiOmics/methylation/results/ctp_analysis/CBMAP_all_sample_methylation_ctp_deconvolution.txt",row.names = T,col.names = T,sep = "\t")
#==============================================================================
# PLOT
#==============================================================================

ctp.long.ls <- lapply(ctp_error.ls, function(x) melt(x, id.var = "IID", variable.name = "celltype", value.name = "CTP"))
# - add comparison identifier
compar <- names(ctp.long.ls)
compar <- c("Luo2022_seq", "Guintivano2013_array")
tmp.ls <- mapply(cbind, ctp.long.ls, "comparison"=compar, SIMPLIFY=F)
# - append list elements
ctp.long.df <- do.call(rbind, tmp.ls)
ctp.gg <- ctp.long.df %>% 
  ggplot(aes(x = celltype, y = CTP, colour = celltype)) + 
  geom_boxplot() + 
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_colour_aaas() +
  facet_grid(~ comparison, scales = "free_x", space = "free_x") +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "gray20")
ggsave("/data/projects/China_Brain_MultiOmics/methylation/results/ctp_analysis/CBMAP_all_sample_methylation_ctp_deconvolution_boxwhiskers.png", ctp.gg)

#sum(pheno_info$AD_diagnosis,na.rm = T)

ctp<-read.delim("/data/projects/China_Brain_MultiOmics/methylation/results/ctp_analysis/CBMAP_all_sample_methylation_ctp_deconvolution.txt")
co_sample<-intersect(rownames(ctp),pheno_info$barcode_id)
ctp<-ctp[co_sample,]

############ order pheno sample
pheno_info<-pheno_info[match(rownames(ctp),pheno_info$barcode_id),]
unique(rownames(ctp)==pheno_info$barcode_id)

######################### Make minimum CTP = 1e-3
ctp<-as.matrix(ctp)
ctp[which(ctp <= 0)] <- 1e-3

#################### Perform clr-transform
ctp<-clr(ctp)

write.table(ctp,"/data/projects/China_Brain_MultiOmics/methylation/results/ctp_analysis/CBMAP_all_sample_methylation_ctp_with_clr_deconvolution.txt",row.names = T,col.names = T,sep = "\t")
