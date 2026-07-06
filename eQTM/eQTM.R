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

load('/eQTM/cis-cpg.RData')
# batch_size <- 200
# total_genes <- length(gene_cpg_list)
# gene_batches <- split(gene_cpg_list, ceiling(seq_along(gene_cpg_list) / batch_size)) # 79 batches
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
fwrite(res, file=paste0('/eQTM/eQTM_batch/batch_clr_',n,'.txt'),
       col.names = T,row.names = F,quote=F,sep='\t')


### summary
library(data.table)
library(ggvenn)
setwd('/eQTM/eQTM_batch')
files = list.files('/eQTM/eQTM_batch')
files = files[grep("batch_clr_",files)]
res_list = list()
for (i in 1:length(files)) {
  res0 = fread(files[i])
  res_list[[length(res_list)+1]] = res0
}
res = rbindlist(res_list);rm(res_list)
res$fdr_global = p.adjust(res$pvalue,'BH')
res <- res %>%
  group_by(cpg) %>%
  mutate(fdr_per_cpg = p.adjust(pvalue, method = "fdr")) %>%
  ungroup()

load('gencode_v46_b38_annotation.RData')
gene_pos = gencode_v46_annotation[which(gencode_v46_annotation$type == 'gene'),c('gene_id','chr','start','end','gene_name')]
gene_pos$gene_id = unlist(str_split(gene_pos$gene_id,"[.]",simplify=T))[,1]
res = merge(res, gene_pos,by='gene_id')

anno_935k = read.csv('EPIC-8v2-0_A1.csv',header=T)
anno_935k = anno_935k[-c(1:6),]
colnames(anno_935k) = anno_935k[1,]
anno_935k = anno_935k[-1,]
anno_935k = anno_935k[!duplicated(anno_935k$Name), ]
anno_935k = anno_935k[,c('Name','CHR','MAPINFO')]
anno_935k$MAPINFO = as.integer(anno_935k$MAPINFO)
res=merge(res, anno_935k, by.x='cpg', by.y='Name')
colnames(res) = c('cpg','gene_id','Beta','se','pvalue','fdr_global','fdr_per_cpg','gene_chr','gene_start','gene_end','gene_name','cpg_chr','cpg_pos')
fwrite(res, file='/eQTM/cis_eQTM_all_res_clr.txt',quote = F,col.names = T,row.names = F, sep='\t')
