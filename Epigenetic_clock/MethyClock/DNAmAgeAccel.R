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
residual_type = c('Age-scaled','Unscaled')
p=list()
for (i in 1:length(disease)) {
  for (j in 1:length(residual_type)) {
    p[[length(p) + 1]] = PreAge_plot(disease[i], residual_type[j])
  }
}
layout <- "
AABBBB
CCCCCC
"
p[[1]]+p[[2]]+p[[3]] + plot_layout(design = layout)
p[[4]]+p[[5]]+p[[6]] + plot_layout(design = layout)