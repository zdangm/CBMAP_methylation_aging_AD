library(stringr)
library(readxl)
library(dplyr)
library(data.table)
library(rhdf5)
library(dmrff)

####rosmap####
###rosmap_pheno
load('ROSMAP_pheno.RData')
pheno_df$projid <- as.character(str_pad(pheno_df$projid, width = 8, pad = "0", side = "left"))
hml_info <- read_excel('dataset_1495_cross-sectional_02-16-2025.xlsx')
pheno_df$age_death <- as.numeric(unlist(hml_info[match(pheno_df$projid, hml_info$projid), 'age_death']))
pheno_df$hml4 <- as.numeric(unlist(hml_info[match(pheno_df$projid, hml_info$projid), "a_score_st4"]))
pheno_df <- pheno_df[!is.na(pheno_df$hml4),] 

idkey = read.csv('ROSMAP_IDkey.csv', header=T)
rin_info = read.csv('ROSMAP_assay_rnaSeq_metadata.csv',header=T)
rin_info = rin_info[,c("specimenID","RIN")]
rin_info$projid = unlist(idkey[match(rin_info$specimenID,idkey$rnaseq_id),"projid"])
rin_info$projid <- as.character(str_pad(rin_info$projid, width = 8, pad = "0", side = "left"))
pheno_df$rin = unlist(rin_info[match(pheno_df$projid,rin_info$projid),"RIN"]) 

cell.pro = read.table('ROSMAP_all_sample_methylation_ctp_with_clr_deconvolution.txt',header=T)
cell_aligned <- cell.pro[pheno_df$barcode_id, , drop = FALSE]
pheno_df <- cbind(pheno_df, cell_aligned)
rownames(pheno_df) = pheno_df$mwas_id

###rosmap_cg
ros_cg <- h5read("ROSMAP_methylation/final_DNAm_invMdat.h5",'Mdat')
colnames(ros_cg) = h5read("ROSMAP_methylation/final_DNAm_invMdat.h5", 'colnames')
rownames(ros_cg) = h5read("ROSMAP_methylation/final_DNAm_invMdat.h5", 'rownames')
sample = intersect(rownames(ros_cg), pheno_df$mwas_id)
ros_cg = ros_cg[sample,]
pheno_df = pheno_df[sample,]
pheno_df[] <- lapply(pheno_df, function(x)
  if (is.numeric(x)) replace(x, is.na(x), median(x, na.rm = TRUE)) else x
)
identical(row.names(pheno_df), rownames(ros_cg))


load("ROSMAP_meth_Ascore.RData")
ros_result <- as.data.frame(t(results))
ros_result$fdr <- p.adjust(ros_result$`Pr(>|t|)`,'BH')
ros_result$cpg <- rownames(ros_result)


####cbmap####
###cbmap_pheno
load('CBMAP_pheno.RData')
all_sample_info <- as.data.frame(read.csv("CBMAP_sample_info.csv",header=T)) 
selected_sample = all_sample_info[!is.na(all_sample_info[,"A_beta_0_3"]),]
selected_sample <- selected_sample %>% 
  filter(
    diag_ALS == 0 & 
      diag_SCZ == 0 & 
      diag_epilepsy == 0 &
      diag_others == 0 &
      OTHER == 0) 
Ascore_sample <- unlist(selected_sample %>% filter(!is.na(A_beta_0_3))%>% dplyr::select(id)) 

###cbmap_cg
cb_cg <- h5read("CBMAP_methylation/final_DNAm_invMdat.h5",'Mdat')
colnames(cb_cg) = h5read("CBMAP_methylation/final_DNAm_invMdat.h5",'colnames')
rownames(cb_cg) = h5read("CBMAP_methylation/final_DNAm_invMdat.h5",'rownames')

sample = Reduce(intersect,list(rownames(cb_cg),Ascore_sample,rownames(pheno_df)))
cb_cg <- cb_cg[sample,]
pheno_df <- merge(pheno_df, selected_sample[,c('id','A_beta_0_3')], by.x='Sample_Name', by.y='id')
rownames(pheno_df) <- pheno_df$Sample_Name
pheno_df <- pheno_df[sample,]

cell.pro = read.table('CBMAP_all_sample_methylation_ctp_with_clr_deconvolution.txt',header=T)
cell_aligned <- cell.pro[pheno_df$barcode_id, , drop = FALSE]
pheno_df <- cbind(pheno_df, cell_aligned)
pheno_df[] <- lapply(pheno_df, function(x)
  if (is.numeric(x)) replace(x, is.na(x), median(x, na.rm = TRUE)) else x
)
identical(rownames(pheno_df), rownames(cb_cg))


load("CBMAP_meth_Ascore.RData")
cb_result <- as.data.frame(t(results))
cb_result$fdr <- p.adjust(cb_result$`Pr(>|t|)`,'BH')
cb_result$cpg <- rownames(cb_result)


common_cg <- intersect(colnames(cb_cg),colnames(ros_cg))
ros_cg <- ros_cg[,common_cg]
cb_cg <- cb_cg[,common_cg]
ros_cg <- t(ros_cg)
cb_cg <- t(cb_cg)
rownames(ros_result) <- ros_result$cpg
rownames(cb_result) <- cb_result$cpg
ros_result <- ros_result[common_cg,]
cb_result <- cb_result[common_cg,]

identical(rownames(ros_cg),ros_result$cpg)
identical(rownames(cb_cg),rownames(cb_result))


anno <- read.csv('EPIC_annotation.csv')
anno <- anno[,c('Name','CHR','MAPINFO')]
ros_result <- merge(ros_result,anno,by.x='cpg',by.y='Name')
cb_result <- merge(cb_result,anno,by.x='cpg',by.y='Name')


###dmrff.pre
ros_dmr_pre <- dmrff.pre(estimate = ros_result$Estimate,se=ros_result$`Std. Error`,methylation = ros_cg,chr = ros_result$CHR,pos = ros_result$MAPINFO)
cb_dmr_pre <- dmrff.pre(estimate = cb_result$Estimate,se=cb_result$`Std. Error`,methylation = cb_cg,chr = cb_result$CHR,pos = cb_result$MAPINFO)

####dmrff.meta
meta_dmr <- dmrff.meta(list(rosmap=ros_dmr_pre,cbmap=cb_dmr_pre))
save(ros_dmr_pre,cb_dmr_pre,meta_dmr,file = 'ascore_dmrff.RData')
dmrs <- meta_dmr$dmrs[which(meta_dmr$dmrs$p.adjust < 0.05 & meta_dmr$dmrs$n >= 3), ]
write.csv(dmrs,'ascore_sig_dmr.csv')

