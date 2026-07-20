library(readxl)
library(org.Hs.eg.db)
library(rhdf5)
library(stringr)
library(openxlsx)
library(bruceR)
library(AnnotationDbi)

###histone genes
top2_enrich_gene <- read.xlsx("linear_GO.xlsx")
top2_genes <- top2_enrich_gene %>%
  filter(Description %in% c("structural constituent of chromatin", "nucleosome")) %>%
  pull(geneID) %>%
  strsplit(",") %>%
  unlist() %>%
  unique() %>%
  trimws()
histon_gene <- mapIds(
  org.Hs.eg.db,
  keys = top2_genes,
  column = "ENSEMBL",
  keytype = "SYMBOL",
  multiVals = "first"
)

###cpg-gene
cpg_gene_res = as.data.frame(fread("cis_eQTM_all_res_clr.txt",header=T))
cpg_gene_res <- cpg_gene_res[cpg_gene_res$gene_id %in% histon_gene, ]
cpg_gene_res <- cpg_gene_res %>%
  group_by(cpg) %>%
  mutate(fdr = p.adjust(pvalue, method = "fdr")) %>%
  ungroup()
cpg_gene_res = cpg_gene_res[which(cpg_gene_res$fdr < 0.05),]#1580

cpg_gene_res <- cpg_gene_res[, c("cpg", "gene_id", "Beta", "pvalue", "fdr")]
colnames(cpg_gene_res) <- c('cpg','gene_id','beta_cpg_gene','p_cpg_gene','fdr_cpg_gene')
cpg <- unique(cpg_gene_res$cpg)#533

####age-cpg
# methylation matrix
cg <- h5read("final_DNAm_Mdat.h5",'Mdat')
colnames(cg) = h5read("final_DNAm_Mdat.h5",'colnames')
rownames(cg) = h5read("final_DNAm_Mdat.h5",'rownames')
cg <- t(cg)

load('CBMAP_pheno.RData')
cell.pro = read.table('CBMAP_all_sample_methylation_ctp_with_clr_deconvolution.txt',header=T)
cell_aligned <- cell.pro[pheno_df$barcode_id, , drop = FALSE]
pheno_df <- cbind(pheno_df, cell_aligned)

sample = intersect(rownames(cg), rownames(pheno_df)) 
cg <- cg[sample,]
pheno_df <- pheno_df[sample,]
pheno_df$sex_male <- as.factor(pheno_df$sex_male)
pheno_df$Sample_Plate <- as.factor(pheno_df$Sample_Plate)
pheno_df$bank <- as.factor(pheno_df$bank)
pheno_df[] <- lapply(pheno_df, function(x)
  if (is.numeric(x)) replace(x, is.na(x), median(x, na.rm = TRUE)) else x
)
cg <- as.data.frame(cg)
cg <- cg[,cpg]
identical(rownames(pheno_df), rownames(cg))

cg$id <- rownames(cg)
dat = merge(cg, pheno_df, by.x='id',by.y='Sample_Name')
ress = data.frame('Estimate'=NA, `Std. Error`=NA, `t value`=NA, `Pr(>|t|)`=NA)
for (i in 2:dim(cg)[2]) {
  model = lm(dat[,i] ~ age + sex_male + bank + PMD + RIN + Exc + Inh + NonN_Astro_FGF3R + NonN_Endo + NonN_Micro + NonN_Oligo_MBP + Sample_Plate, data=dat)
  res = summary(model)
  res = res$coefficients[2,]
  ress[i,] = res
}
ress = ress[-1,]
ress$cpg = colnames(dat)[2:dim(cg)[2]]
ress$fdr <- p.adjust(ress$Pr...t..,'BH')
ress <- ress[ress$fdr<0.05,]
ress <- ress[,c("cpg","Estimate", "Pr...t..", "fdr")]
colnames(ress) <- c('cpg','beta_age_cpg','p_age_cpg','fdr_age_cpg')
combined_data <- merge(cpg_gene_res,ress,by='cpg')

###age-gene
gene <- unique(combined_data$gene_id)
exp = read.table('RNAseq_data_qn_comb.txt',check.names = F)
exp = as.data.frame(t(exp))
exp = exp[,colnames(exp) %in% gene]

load('CBMAP_pheno.RData')
cell.pro = read.table('CBMAP_all_sample_methylation_ctp_with_clr_deconvolution.txt',header=T)
cell_aligned <- cell.pro[pheno_df$barcode_id, , drop = FALSE]
pheno_df <- cbind(pheno_df, cell_aligned)

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

exp$id <- rownames(exp)
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
ress <- ress[ress$fdr<0.05,]

ress <- ress[,c('gene','Estimate','Pr...t..','fdr')]
colnames(ress) <- c('gene','beta_age_gene','p_age_gene','fdr_age_gene')
combined_all <- merge(combined_data,ress,by.x='gene_id',by.y='gene')

#direction filter
combined_all_filtered <- combined_all %>%
  mutate(
    same_sign = sign(beta_age_cpg) == sign(beta_cpg_gene),
    condition_met = case_when(
      same_sign & beta_age_gene > 0 ~ TRUE, 
      !same_sign & beta_age_gene < 0 ~ TRUE, 
      TRUE ~ FALSE  
    )
  ) %>%
  filter(condition_met) %>%
  dplyr::select(-same_sign, -condition_met)
write.csv(combined_all_filtered,'histone_association_direction.csv')

####mediation analysis
cg <- h5read("final_DNAm_Mdat.h5",'Mdat')
colnames(cg) = h5read("final_DNAm_Mdat.h5",'colnames')
rownames(cg) = h5read("final_DNAm_Mdat.h5",'rownames')
cg <- t(cg)

top2_enrich_gene <- read.csv('histone_association_direction.csv')

exp = read.table('RNAseq_data_qn_comb.txt',check.names = F)

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

mediation_res <- data.frame('cpg' = rep(NA, nrow(top2_enrich_gene)), 
                            'Mediated.Prop' = rep(NA, nrow(top2_enrich_gene)), 
                            'estimate' = rep(NA, nrow(top2_enrich_gene)), 
                            'P' = rep(NA, nrow(top2_enrich_gene)))
for (i in 1:nrow(top2_enrich_gene)) {
  cpg <- as.character(top2_enrich_gene[i, "cpg"])
  gene <- as.character(top2_enrich_gene[i, "gene_id"])
  
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
colnames(mediation_res)[3:5] <- c('med_beta','med_p','med_fdr')

original_data <- read.csv("histone_association_direction.csv")
original_data$gene_name <- mapIds(
  org.Hs.eg.db,
  keys = original_data$gene_id,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

identical(original_data$cpg,mediation_res$cpg)
res_table <- cbind(original_data,mediation_res)

med_sig <- res_table[res_table$med_fdr<0.05,]
write.xlsx(med_sig,'histone_med_sig.xlsx')

