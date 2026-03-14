library(stringr)
library(data.table)
library(MASS)
library(rhdf5)
library(dplyr)
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
pheno_df$sex_male <- as.factor(pheno_df$sex_male)
pheno_df$Sample_Plate <- as.factor(pheno_df$Sample_Plate)
pheno_df$bank <- as.factor(pheno_df$bank)
pheno_df$AD <- factor(pheno_df$AD, levels=c(0,1,2,3), ordered = T)

cell.pro = read.table('CBMAP_all_sample_methylation_ctp_with_clr_deconvolution.txt',header=T)
cell_aligned <- cell.pro[pheno_df$barcode_id, , drop = FALSE]
pheno_df <- cbind(pheno_df, cell_aligned)
pheno_df[] <- lapply(pheno_df, function(x)
  if (is.numeric(x)) replace(x, is.na(x), median(x, na.rm = TRUE)) else x
)
identical(rownames(pheno_df), rownames(cg))
dim(pheno_df)[1]

###### meth~age #####
res_M = data.frame('cpg'=NA, 'beta'=NA,'se'=NA,'t'=NA,'p'=NA)
for (i in 1:ncol(cg)) {
  fit1 = lm(cg[,i] ~ age + sex_male + bank + PMD + RIN + Exc + Inh + NonN_Astro_FGF3R + NonN_Endo + NonN_Micro + NonN_Oligo_MBP + NonN_OPC + Sample_Plate, data = pheno_df)
  fit1 = summary(fit1)
  res_M[nrow(res_M)+1,1] = colnames(cg)[i]
  res_M[nrow(res_M),2:5] = fit1$coefficients[2,]
}
res_M = res_M[-1,]
##### AD ~ age #####
fit3 = polr(AD ~ age + sex_male + bank + PMD + RIN + Exc + Inh + NonN_Astro_FGF3R + NonN_Endo + NonN_Micro + NonN_Oligo_MBP + NonN_OPC + Sample_Plate, data = pheno_df,
            Hess = TRUE, method = 'logistic')
fit3 = summary(fit3)
c = fit3$coefficients['age',1]
se_c = fit3$coefficients['age',2]
##### AD ~ age + methy #####
res_b_age = data.frame('cpg'=NA, 'beta'=NA,'se'=NA)
res_b_methy = data.frame('cpg'=NA, 'beta'=NA,'se'=NA)
for (i in 1:ncol(cg)) {
  fit2 = polr(AD ~ age + cg[,i] + sex_male + bank + PMD + RIN + Exc + Inh + NonN_Astro_FGF3R + NonN_Endo + NonN_Micro + NonN_Oligo_MBP + NonN_OPC + Sample_Plate, data = pheno_df,
              Hess = TRUE, method = 'logistic')
  fit2 = summary(fit2)
  res_b_age[nrow(res_b_age)+1,1] = colnames(cg)[i]
  res_b_age[nrow(res_b_age),2:3] = fit2$coefficients[1,1:2]
  res_b_methy[nrow(res_b_methy)+1,1] = colnames(cg)[i]
  res_b_methy[nrow(res_b_methy),2:3] = fit2$coefficients[2,1:2]
}
res_b_age = res_b_age[-1,]
res_b_methy = res_b_methy[-1,]
###### standardized regression coefficients #####
sd_M = apply(cg,2,sd)
sd_X = sd(pheno_df$age)
var_X_age = var(pheno_df$age)
var_M_methy = apply(cg, 2, var)
cov_Xage_Methy = apply(cg, 2, function(x){cov(pheno_df$age, x)})
identical(names(var_M_methy),res_b_age$cpg)
identical(names(cov_Xage_Methy), res_b_methy$cpg)
var_Y = ((res_b_age$beta)^2) * var_X_age + ((res_b_methy$beta)^2) * var_M_methy + 2 * res_b_age$beta * (res_b_methy$beta) * cov_Xage_Methy + (pi^2)/3

b_std = res_b_methy$beta * sd_M / sqrt(var_Y)
se_b_std = res_b_methy$se * sd_M / sqrt(var_Y)

a_std = res_M$beta * sd_X / sd_M
se_a_std = res_M$se * sd_X / sd_M

cc_std = res_b_age$beta * sd_X / sqrt(var_Y)
c_std = c * sd_X / sqrt(c^2 * var_X_age + (pi^2)/3)

##### results #####
ab_std <- a_std * b_std
se_ab_std <- sqrt(
  (a_std^2) * (se_b_std^2) +
    (b_std^2) * (se_a_std^2)
)
z_value <- ab_std / se_ab_std
p_value <- 2 * (1 - pnorm(abs(z_value)))
sum(p_value<0.05)
CI_lower <- ab_std - qnorm(0.975) * se_ab_std
CI_upper <- ab_std + qnorm(0.975) * se_ab_std
medi_pro <- ab_std / (ab_std + cc_std)

res = data.frame('cpg'=colnames(cg), 'indirect_effect' = ab_std, 'indirect_se' = se_ab_std, 
                 'z_value' = z_value, 'p_value' = p_value,
                 'CI_lower' = CI_lower, 'CI_upper' = CI_upper, 'medi_pro' = medi_pro)
res$fdr = p.adjust(res$p_value,'BH')
sum(res$fdr < 0.05)
write.table(res,file='mediation_ADNC_res.txt',col.names = T,row.names = F,quote = F,sep='\t')
