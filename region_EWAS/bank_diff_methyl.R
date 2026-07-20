library(rhdf5)
library(data.table)
library(dplyr)
library(stringr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)
library(openxlsx)
library(limma)



# methylation matrix
cg <- h5read("final_DNAm_invMdat.h5",'Mdat')
colnames(cg) = h5read("final_DNAm_invMdat.h5",'colnames')
rownames(cg) = h5read("final_DNAm_invMdat.h5",'rownames')
cg = cg[,grep("cg",colnames(cg))] 

load('CBMAP_pheno.RData')
cell.pro = read.table('CBMAP_all_sample_methylation_ctp_with_clr_deconvolution.txt',header=T)
cell_aligned <- cell.pro[pheno_df$barcode_id, , drop = FALSE]
pheno_df <- cbind(pheno_df, cell_aligned)

sample = intersect(rownames(cg), rownames(pheno_df)) 
cg <- cg[sample,]
pheno_df <- pheno_df[sample,]
pheno_df$sex_male <- as.factor(pheno_df$sex_male)
pheno_df$Sample_Plate <- as.factor(pheno_df$Sample_Plate)
pheno_df[] <- lapply(pheno_df, function(x)
  if (is.numeric(x)) replace(x, is.na(x), median(x, na.rm = TRUE)) else x
)
identical(rownames(pheno_df), rownames(cg))

####pumc-zju
pheno_pumc_zju <- pheno_df[pheno_df$bank %in% c('pumc','zju'),]
pheno_pumc_zju$bank <- factor(pheno_pumc_zju$bank, levels = c('pumc','zju'))
cg_pumc_zju <- cg[rownames(pheno_pumc_zju),]
identical(rownames(cg_pumc_zju),rownames(pheno_pumc_zju))

cg_t <- t(cg_pumc_zju)

design <- model.matrix(~ bank + sex_male + age + PMD + RIN + 
                         Exc + Inh + NonN_Astro_FGF3R + NonN_Endo + 
                         NonN_Micro + NonN_Oligo_MBP + Sample_Plate, 
                       data = pheno_pumc_zju)
fit <- lmFit(cg_t, design)   
fit2 <- eBayes(fit)
res <- topTable(fit2, coef = "bankzju", number = Inf, sort.by = "none")
res_pumc_zju <- res[res$adj.P.Val<0.05,]#141,354

####pumc-csu
pheno_pumc_csu <- pheno_df[pheno_df$bank %in% c('pumc','csu'),]
pheno_pumc_csu$bank <- factor(pheno_pumc_csu$bank, levels = c('pumc','csu'))
cg_pumc_csu <- cg[rownames(pheno_pumc_csu),]
identical(rownames(cg_pumc_csu),rownames(pheno_pumc_csu))

cg_t <- t(cg_pumc_csu)

design <- model.matrix(~ bank + sex_male + age + PMD + RIN + 
                         Exc + Inh + NonN_Astro_FGF3R + NonN_Endo + 
                         NonN_Micro + NonN_Oligo_MBP + Sample_Plate, 
                       data = pheno_pumc_csu)
fit <- lmFit(cg_t, design)   
fit2 <- eBayes(fit)
res <- topTable(fit2, coef = "bankcsu", number = Inf, sort.by = "none")
res_pumc_csu <- res[res$adj.P.Val<0.05,]#27,790

####zju-csu
pheno_zju_csu <- pheno_df[pheno_df$bank %in% c('zju','csu'),]
pheno_zju_csu$bank <- factor(pheno_zju_csu$bank, levels = c('zju','csu'))
cg_zju_csu <- cg[rownames(pheno_zju_csu),]
identical(rownames(cg_zju_csu),rownames(pheno_zju_csu))

cg_t <- t(cg_zju_csu)

design <- model.matrix(~ bank + sex_male + age + PMD + RIN + 
                         Exc + Inh + NonN_Astro_FGF3R + NonN_Endo + 
                         NonN_Micro + NonN_Oligo_MBP + Sample_Plate, 
                       data = pheno_zju_csu)
fit <- lmFit(cg_t, design)   
fit2 <- eBayes(fit)
res <- topTable(fit2, coef = "bankcsu", number = Inf, sort.by = "none")
res_zju_csu <- res[res$P.Value>0.05,]#673,847

####significant region-cpgs
common_probes <- Reduce(intersect, list(
  rownames(res_pumc_zju),
  rownames(res_pumc_csu),
  rownames(res_zju_csu)
))#3,391
colnames(res_pumc_csu_sig) <- paste0('pumc_csu_',colnames(res_pumc_csu_sig))
colnames(res_pumc_zju_sig) <- paste0('pumc_zju_',colnames(res_pumc_zju_sig))
colnames(res_zju_csu_nonsig) <- paste0('zju_csu_',colnames(res_zju_csu_nonsig))
final_res <- cbind(res_pumc_csu_sig,res_pumc_zju_sig,res_zju_csu_nonsig)
final_res$cpg <- rownames(final_res)
write.xlsx(final_res,'significant_region_cpgs.xlsx')
