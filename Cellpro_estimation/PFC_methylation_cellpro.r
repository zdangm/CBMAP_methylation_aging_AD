rm(list = ls())
gc()
options(mc.cores = 10)

library(meffil)
library(data.table)
library(dplyr)
library(compositions)

############# create.samplesheet for CBMAP
indir_zk <- "/mirror_raw_data_DNAm/zk"
indir_sh <- "/mirror_raw_data_DNAm/sh"

samplesheet_zk <- meffil.create.samplesheet(indir_zk)

samplesheet_sh <- meffil.create.samplesheet(indir_sh)

############### 1057 samples
samplesheet<-rbind(samplesheet_zk,samplesheet_sh)

## Step 2: read phenotype and batch information
pdat <- read.csv('phenotype.csv', header=T)
rownames(pdat) = pdat$barcode_id


samplesheet<-merge(samplesheet,pdat,by.x="Sample_Name",by.y="barcode_id")

################### set sex info
samplesheet$Sex[samplesheet$sex_male==1]<-"M"
samplesheet$Sex[samplesheet$sex_male==0]<-"F"


# set qc.parameters
qc.parameters <- meffil.qc.parameters(
  beadnum.samples.threshold             = 0.1,
  detectionp.samples.threshold          = 0.1,
  detectionp.cpgs.threshold             = 0.1, 
  beadnum.cpgs.threshold                = 0.1,
  sex.outlier.sd                        = 5,
  snp.concordance.threshold             = 0.95,
  sample.genotype.concordance.threshold = 0.8
)

# Generate the qc.objects
qc.objects <- meffil.qc(samplesheet,
                        cell.type.reference = "guintivano dlpfc", verbose = TRUE)

# Generate the QC summary
qc.summary <- meffil.qc.summary(qc.objects, parameters = qc.parameters)

# Generate outlier
outlier <- qc.summary$bad.samples

# remove outlier
qc.objects <- meffil.remove.samples(qc.objects, outlier$sample.name)

# set qc.parameters
qc.parameters <- meffil.qc.parameters(
  beadnum.samples.threshold             = 0.1,
  detectionp.samples.threshold          = 0.1,
  detectionp.cpgs.threshold             = 0.1, 
  beadnum.cpgs.threshold                = 0.1,
  sex.outlier.sd                        = 5,
  snp.concordance.threshold             = 0.95,
  sample.genotype.concordance.threshold = 0.8
)

# Generate the QC summary
qc.summary <- meffil.qc.summary(qc.objects, parameters = qc.parameters)

# Perform quantile normalization
norm.objects <- meffil.normalize.quantiles(qc.objects, fixed.effects = "Sample_Plate", number.pcs = 10) 

# Generate normalized probe values
norm.beta <- meffil.normalize.samples(norm.objects, 
                                      cpglist.remove = qc.summary$bad.cpgs$name) 

save(norm.beta, file = paste("/ctp_analysis/CBMAP_meffil_", "norm.beta", ".Robj", sep = ""))

#==============================================================================
#
# Format files for OSCA
#
#==============================================================================

#------------------------------------------------------------------------------
# Directories
#------------------------------------------------------------------------------

norm_beta_dir <- paste("/ctp_analysis/CBMAP_meffil_", "norm.beta", ".Robj", sep = "")
out_dir <- "/ctp_analysis/OSCA/format_for_osca/CBMAP_norm_beta_OSCA_input.txt"

#------------------------------------------------------------------------------
# Read-in data and libraries
#------------------------------------------------------------------------------

print("Loading normalised beta table (from meffil)")
load(norm_beta_dir)

print("Tidying table")
norm.beta.df <- data.frame(IID = rownames(norm.beta), norm.beta, row.names = NULL)
colnames(norm.beta.df) <- gsub("X", "", colnames(norm.beta.df))

print("Writing table")
write.table(norm.beta.df, out_dir, 
            col.names = T, row.names = F, quote = F, sep = "\t")

################# process 935k EPIC for CBMAP
anno <- fread('/935k_annotation/EPIC-8v2-0_A1.csv')
colnames(anno) <- as.character(anno[8,])
anno <- anno[-c(1:8),]
anno <- anno[!duplicated(anno$Name, fromLast = TRUE), ]

# extract autosomal_probes
autosomal_probes <- anno$Name[anno$CHR %in% paste0("chr", 1:22)]
write.table(autosomal_probes, file="/ctp_analysis/OSCA/935k_EPIC_aut.probe", row.names=FALSE, col.names=FALSE, quote=FALSE)


# extract masked_probes
masked_probes <- anno$Name[!(anno$CHR %in% paste0("chr", 1:22))]
write.table(masked_probes, file="/ctp_analysis/OSCA/935k_EPIC_mask.probe", row.names=FALSE, col.names=FALSE, quote=FALSE)
