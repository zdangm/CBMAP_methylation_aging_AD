options(mc.cores = 10)

library(meffil)
library(data.table)
library(dplyr)
library(compositions)

############# create.samplesheet
indir <- "/data/ROSMAP_methylaition/idat_data/"
samplesheet <- meffil.create.samplesheet(indir)

################## read clinical info
synapse_meta = fread('/ROSMAP/Data/Epigenetics/Epigenetics (DNA methylation array)/IDAT Files/SYNAPSE_METADATA_MANIFEST.tsv')[,c('path','specimenID')]
cov = fread('/ROSMAP/Data/Epigenetics/ROSMAP_arrayMethylation_covariates.tsv')
clinical_info = fread('/ROSMAP/Data/Metadata/ROSMAP_clinical.csv')
clinical_info$id = paste0(clinical_info$Study, '', clinical_info$projid)
idkey = read.csv('/ROSMAP/ROSMAP_IDkey.csv', header=T)

##phenotype information
mer1=merge(idkey[,c('projid','mwas_id')], clinical_info, by='projid')
pheno_info = merge(mer1, cov, by.x='mwas_id', by.y='Sample')
rm(mer1)
pheno_info$barcode_id = paste0(pheno_info$Sentrix_ID,'_',pheno_info$Sentrix_Position)
pheno_info = distinct(pheno_info)
rownames(pheno_info) = pheno_info$barcode_id
#5822038012_R02C01  TBI-AUTO73203-PT-35OA  has missing sex and age_death
pheno_info = pheno_info[-which(rownames(pheno_info)=='5822038012_R02C01'),]


################ add phenotype
samplesheet<-merge(samplesheet,pheno_info,by.x="Sample_Name",by.y="barcode_id")

################### set sex info
samplesheet$Sex[samplesheet$msex==1]<-"M"
samplesheet$Sex[samplesheet$msex==0]<-"F"

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

save(norm.beta, file = paste("/methylation/results/ctp_analysis/ROSMAP_meffil_", "norm.beta", ".Robj", sep = ""))

