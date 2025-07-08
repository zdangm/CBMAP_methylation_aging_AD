library(data.table)
library(patchwork)
library(dplyr)
library(cowplot)
library(ggplot2)
library(readxl)

res = read.table('case_agePred.txt',header=T)
res = res[,c('ID','age','CBMAP_EN','CBMAP_PC','CorticalAge','PCBrainAge','MultiAge','BloodAge','HannumAge','PhenoAge')]
hc = read.table('hc_test_agePred.txt',header=T)
hc = hc[,c('ID','age','CBMAP_EN','CBMAP_PC','CorticalAge','PCBrainAge','MultiAge','BloodAge','HannumAge','PhenoAge')]
res = rbind(res, hc)
load('CBMAP_pheno.RData')
res = merge(res,pheno_df, by.x='ID',by.y='Sample_Name')
res$sex_male = factor(res$sex_male,levels=c(0,1))
source('DNAmAgeAccel_fun.R')

######## Associations between DNAmAgeAccel and diseases #########
disease = c('HML1','HML2','braak')
p=list()
for (i in 1:length(disease)) {
  p[[i]] = PreAge_plot(disease[i])
}
layout <- "
AABBBB
CCCCCC
"
p[[1]]+p[[2]]+p[[3]] + plot_layout(design = layout)