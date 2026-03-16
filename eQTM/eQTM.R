args = as.numeric(commandArgs(TRUE))
print(args)
n = args


library(data.table)
library(stringr)
library(dplyr)
library(rhdf5)

residual_methy <- h5read("residual_methy",'residual_methy')
colnames(residual_methy) = h5read("residual_methy",'colnames')
rownames(residual_methy) = h5read("residual_methy",'rownames')

residual_exp <- fread('residual_exp.txt',header=T)

load('cis-cpg.RData')
gene_cpg_list = gene_cpg_list[((n-1)*200+1):min(n*200,length(gene_cpg_list))]

common_samples <- intersect(rownames(residual_methy), residual_exp$V1) 
residual_methy <- residual_methy[common_samples, ]
residual_exp = as.data.frame(residual_exp)
rownames(residual_exp) = residual_exp$V1
residual_exp = residual_exp[,2:ncol(residual_exp)]
residual_exp = residual_exp[common_samples,]
identical(rownames(residual_exp),rownames(residual_methy))

res = data.frame('gene_id'=NA,'cpg'=NA,'Beta'=NA,'se'=NA,'pvalue'=NA)
for (i in 1:length(gene_cpg_list)) {
  gene = gene_cpg_list[[i]]$gene_id
  cpgs = gene_cpg_list[[i]]$cis_cpgs
  exp = residual_exp[,gene]
  cg = residual_methy[,cpgs]
  
  for (cpg in cpgs) {
    fit = lm(exp ~ cg[,cpg])
    fit = summary(fit)
    res0 = fit$coefficients[2,c(1,2,4)]
    res[nrow(res)+1,] = c(gene,cpg,res0)
  }
}
res = res[-1,]
fwrite(res, file=paste0('/eQTM_batch/batch_',n,'.txt'),
       col.names = T,row.names = F,quote=F,sep='\t')