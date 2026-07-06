options(mc.cores = 10)

library(meffil)
library(data.table)
library(dplyr)
library(compositions)

############# create.samplesheet for CBMAP
indir_zk <- "/methylation/data/CBMAP/DNAm_processed/mirror_raw_data_DNAm/zk"
indir_sh <- "/methylation/data/CBMAP/DNAm_processed/mirror_raw_data_DNAm/sh"

samplesheet_zk <- meffil.create.samplesheet(indir_zk)

samplesheet_sh <- meffil.create.samplesheet(indir_sh)

############### 
samplesheet<-rbind(samplesheet_zk,samplesheet_sh)

## Step 2: read phenotype and batch information
pdat <- read.csv('/methylation/data/CBMAP/DNAm_processed/phenotype.csv', header=T)
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

save(norm.beta, file = paste("/methylation/results/ctp_analysis/CBMAP_meffil_", "norm.beta", ".Robj", sep = ""))

