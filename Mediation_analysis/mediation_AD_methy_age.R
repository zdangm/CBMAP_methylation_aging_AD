library(rhdf5)
library(dplyr)
library(data.table)
library(ggplot2)
library(bruceR)
library(patchwork)
library(ggrepel)
library(readxl)

# load meta-AD-meth results
load('meta_hml4level.RData')
hml_cpg = meta_final_ordered_df[which(meta_final_ordered_df$fdr<0.05),'cpg']  
load('CBMAP_pheno.RData')
# select case and control
all_sample_info <- as.data.frame(read.csv("CBMAP_sample_info.csv",header=T)) 
selected_sample = all_sample_info[!is.na(all_sample_info[,'ADNC']),]
selected_sample <- selected_sample %>% 
  filter(
    diag_ALS == 0 & 
      diag_SCZ == 0 & 
      diag_epilepsy == 0 &
      diag_others == 0 &
      OTHER == 0) 
dim(selected_sample)
complete_rows <- complete.cases(selected_sample[, c("ADNC", "ADNC.L", "ADNC.M", "ADNC.H")])
selected_sample <- selected_sample[complete_rows, ]
selected_sample$marker <- ifelse(selected_sample$ADNC == 0, 0,
                                 ifelse(selected_sample$ADNC.L == 1, 1,
                                        ifelse(selected_sample$ADNC.M == 1, 2,
                                               ifelse(selected_sample$ADNC.H == 1, 3, NA))))
exclued_id <- selected_sample[which(selected_sample$marker == 0 & (selected_sample$LBD == 1 | selected_sample$ARTAG_Gray.matter == 1 |
                                                                     selected_sample$ARTAG_Perivascular == 1 | selected_sample$ARTAG_Subpial == 1 | selected_sample$diag_PD == 1 | selected_sample$diag_dementia == 1)),'id']
pheno_df$AD <- selected_sample[match(pheno_df$Sample_Name, selected_sample$id), 'marker']
pheno_df <- pheno_df[!is.na(pheno_df$AD),] 
pheno_df <- pheno_df[which(! pheno_df$Sample_Name %chin% exclued_id),]
hml_sample <- pheno_df$Sample_Name 

# methylation matrix
cg <- h5read("final_DNAm_invMdat.h5",'Mdat')
colnames(cg) = h5read("final_DNAm_invMdat.h5",'colnames')
rownames(cg) = h5read("final_DNAm_invMdat.h5",'rownames')
sample = Reduce(intersect,list(rownames(cg),hml_sample,pheno_df$Sample_Name)) 
length(sample)
cg <- cg[sample,hml_cpg]
pheno_df <- pheno_df[sample,]
pheno_df$sex_male <- factor(pheno_df$sex_male)
pheno_df$NeuN_pos <- as.numeric(pheno_df$NeuN_pos)
pheno_df$AD <- as.numeric(pheno_df$AD) 
pheno_df$age <- as.numeric(pheno_df$age)
identical(rownames(pheno_df), rownames(cg))

##### mediation analysis #####
mediation_res = data.frame('cpg'=NA, 'Mediated.Prop'=NA, 'estimate'=NA, 'P'=NA)
o = 0
for (i in 1:length(hml_cpg)) {
  dat = data.frame('methy'=cg[,hml_cpg[i]], 'age'=pheno_df$age, 'AD'=pheno_df$AD, 'sex'=pheno_df$sex_male, 'NeuN_pos'=pheno_df$NeuN_pos)
  contcont <- PROCESS(data=dat, y='AD', x='age', meds='methy', covs=c('sex','NeuN_pos'), nsim=1000, seed=2, digits=7)
  o <- o + 1
  mediation_res[o,1] = hml_cpg[i]
  mediation_res[o,2]= contcont$results[[1]]$mediation[1,1] / contcont$results[[1]]$mediation[3,1]
  mediation_res[o,3]= contcont$results[[1]]$mediation[1,1]
  mediation_res[o,4]= contcont$results[[1]]$mediation[1,'pval']
}
mediation_res$fdr = p.adjust(mediation_res$P, 'BH')
write.table(mediation_res, file='mediation_hml_meth_age.txt', row.names=F, col.names=T, quote=F)