library(data.table)
library(dplyr)
library(org.Hs.eg.db)
library(stringr)
library(openxlsx)
library(rhdf5)
library(bruceR)
library(mediation)
load("r2diff2_kmeans_cluster_scale.RData")

# nonlinear_hypomethylated pattern as example
cpg_c2 <- rownames(clustering.res.df[clustering.res.df$Cluster == 2, ])

cpg_gene_res = as.data.frame(fread("cis_eQTM_all_res_clr.txt",header=T))
cpg_gene_res <- cpg_gene_res[cpg_gene_res$cpg %in% cpg_c2, ]
cpg_gene_res <- cpg_gene_res %>%
  group_by(cpg) %>%
  mutate(fdr = p.adjust(pvalue, method = "fdr")) %>%
  ungroup()
cpg_gene_res = cpg_gene_res[which(cpg_gene_res$fdr < 0.05),]
cpg_gene_res$gene_name <- mapIds(
  org.Hs.eg.db,
  keys = cpg_gene_res$gene_id,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# negative cpg-gene association as example
cpg_gene_neg <- cpg_gene_res[cpg_gene_res$Beta<0,]
gene <- cpg_gene_neg$gene_id

# age-gene
exp = read.table('RNAseq_data_qn_comb.txt',check.names = F)
exp = as.data.frame(t(exp))
exp <- exp[,gene]
exp$id = rownames(exp)

load('CBMAP_pheno.RData')
cell.pro = read.table('CBMAP_all_sample_methylation_ctp_with_clr_deconvolution.txt',header=T)
cell_aligned <- cell.pro[pheno_df$barcode_id, , drop = FALSE]
pheno_df <- cbind(pheno_df, cell_aligned)
pheno_df <- pheno_df[pheno_df$age<70,]

sample = intersect(rownames(exp), rownames(pheno_df)) 
exp <- exp[sample,]
pheno_df <- pheno_df[sample,]
pheno_df$sex_male <- as.factor(pheno_df$sex_male)
pheno_df$Sample_Plate <- as.factor(pheno_df$Sample_Plate)
pheno_df$bank <- as.factor(pheno_df$bank)
pheno_df[] <- lapply(pheno_df, function(x)
  if (is.numeric(x)) replace(x, is.na(x), median(x, na.rm = TRUE)) else x
)
identical(rownames(pheno_df), rownames(exp))

dat = merge(exp, pheno_df, by.x='id',by.y='Sample_Name')

ress = data.frame('Estimate'=NA, `Std. Error`=NA, `t value`=NA, `Pr(>|t|)`=NA)
for (i in 2:dim(exp)[2]) {
  model = lm(dat[,i] ~ age + sex_male + bank + PMD + RIN + Exc + Inh + NonN_Astro_FGF3R + NonN_Endo + NonN_Micro + NonN_Oligo_MBP + Sample_Plate, data=dat)
  res = summary(model)
  res = res$coefficients[2,]
  ress[i,] = res
}
ress = ress[-1,]
ress$gene = colnames(dat)[2:dim(exp)[2]]
ress$fdr <- p.adjust(ress$Pr...t..,'BH')
ress <- ress[ress$fdr<0.05 &ress$Estimate>0,]
ress <- ress[,c('gene','Estimate','Pr...t..','fdr')]
colnames(ress) <- c('gene','beta_age_gene','p_age_gene','fdr_age_gene')

cpg_gene_neg <- cpg_gene_neg[cpg_gene_neg$gene_id%in%ress$gene,]
cpg_gene_neg <- cpg_gene_neg[,-c(4,6,7)]
colnames(cpg_gene_neg)[c(3,4,11)] <- c('cpg_gene_estimate','cpg_gene_p','cpg_gene_fdr')
output <- merge(cpg_gene_neg,ress,by.x='gene_id',by.y='gene')
write.xlsx(output,'c2_neg_neg.xlsx')

####mediation analysis
data <- read.xlsx('c2_neg_neg.xlsx')
#cg
cg <- h5read("final_DNAm_Mdat.h5",'Mdat')
colnames(cg) = h5read("final_DNAm_Mdat.h5",'colnames')
rownames(cg) = h5read("final_DNAm_Mdat.h5",'rownames')
cg <- t(cg)

#exp
exp = read.table('RNAseq_data_qn_comb.txt',check.names = F)

#pheno
load('CBMAP_pheno.RData')
cell.pro = read.table('CBMAP_all_sample_methylation_ctp_with_clr_deconvolution.txt',header=T)
cell_aligned <- cell.pro[pheno_df$barcode_id, , drop = FALSE]
pheno_df <- cbind(pheno_df, cell_aligned)
pheno_df$sex_male <- as.factor(pheno_df$sex_male)
pheno_df$Sample_Plate <- as.factor(pheno_df$Sample_Plate)
pheno_df$bank <- as.factor(pheno_df$bank)
pheno_df[] <- lapply(pheno_df, function(x)
  if (is.numeric(x)) replace(x, is.na(x), median(x, na.rm = TRUE)) else x
)
pheno_df <- pheno_df[pheno_df$age<70,]

mediation_res <- data.frame('cpg' = rep(NA, nrow(data)),
                            'Mediated.Prop' = rep(NA, nrow(data)),
                            'estimate' = rep(NA, nrow(data)),
                            'P' = rep(NA, nrow(data)))
for (i in 1:nrow(data)) {
  cpg <- as.character(data[i, "cpg"])
  gene <- as.character(data[i, "gene_id"])
  
  cg1 <- as.data.frame(cg[,cpg,drop=F])
  cg1$id <- rownames(cg1)
  
  exp_subset <- exp[gene,]
  exp_subset <- as.data.frame(t(exp_subset))
  exp_subset$id <- rownames(exp_subset)
  
  common_ids <- intersect(cg1$id, exp_subset$id)
  cg1 <- cg1[cg1$id %in% common_ids, ]
  exp_subset <- exp_subset[exp_subset$id %in% common_ids, ]
  dat <- merge(x=cg1,y=exp_subset,by='id',all=F)
  common_ids <- intersect(dat$id, pheno_df$Sample_Name)
  dat <- dat[dat$id %in% common_ids, ]
  pheno_df_subset <- pheno_df[pheno_df$Sample_Name %in% common_ids, ]
  dat = merge(dat, pheno_df_subset, by.x='id',by.y='Sample_Name')
  
  dat = data.frame('methy' = dat[,2], 'age' = dat$age, 'exp' = dat[,3], 'sex' = dat$sex_male, 'rin' = dat$RIN, 'pmd' = dat$PMD,'bank'=dat$bank, 'Exc'=dat$Exc,
                   'Inh'=dat$Inh,  'NonN_Astro_FGF3R'=dat$NonN_Astro_FGF3R,'NonN_Endo'=dat$NonN_Endo,'NonN_Micro'=dat$NonN_Micro, 'NonN_Oligo_MBP'=dat$NonN_Oligo_MBP,
                   'Sample_Plate'=dat$Sample_Plate)
  
  contcont <- PROCESS(data = dat, y = 'exp', x = 'age', meds = 'methy', covs = c('sex', 'rin', 'pmd','bank','Exc', 'Inh',
                                                                                 'NonN_Astro_FGF3R','NonN_Endo','NonN_Micro', 
                                                                                 'NonN_Oligo_MBP', 'Sample_Plate'), nsim = 1000, seed = 2, digits = 7)
  mediation_res[i, 1] = colnames(cg1)[1]
  mediation_res[i, 2] = contcont$results[[1]]$mediation[1, 1] / contcont$results[[1]]$mediation[3, 1]
  mediation_res[i, 3] = contcont$results[[1]]$mediation[1, 1]
  mediation_res[i, 4] = contcont$results[[1]]$mediation[1, 'pval']
}

mediation_res$fdr <- p.adjust(mediation_res$P,'BH')
sum(mediation_res$fdr<0.05)
colnames(mediation_res)[3:5] <- c('med_beta','med_p','med_fdr')

identical(data$cpg,mediation_res$cpg)
res_table <- cbind(data,mediation_res)
res_table <- res_table[res_table$med_fdr<0.05,]
write.xlsx(res_table,'c2_neg_neg_med.xlsx')
