library(writexl)
library(rhdf5)
library(stringr)
library(data.table)
###blood
load("meta_meth_age_blood_60years.RData")
blood_res = res
sum(blood_res$p_bacon_fdr<0.05,na.rm=T)

###CBMAP
load("CBMAP_meth_age_noDementia_60years.RData")
cbmap_res <- res
sum(cbmap_res$p_bacon_fdr<0.05,na.rm=T)

###ROSMAP
load("ROSMAP_meth_age_noDementia.RData")
rosmap_res <- res
sum(rosmap_res$p_bacon_fdr<0.05,na.rm=T)

cpg = Reduce(intersect,list(blood_res$cpg,cbmap_res$cpg,rosmap_res$cpg))

blood_res = blood_res[which(blood_res$cpg %chin% cpg & blood_res$p_bacon_fdr < 0.05),]
cbmap_res = cbmap_res[which(cbmap_res$cpg %chin% cpg & cbmap_res$p_bacon_fdr < 0.05),] 
rosmap_res = rosmap_res[which(rosmap_res$cpg %chin% cpg & rosmap_res$p_bacon_fdr < 0.05),] 
rownames(blood_res) <- blood_res$cpg
rownames(cbmap_res) <- cbmap_res$cpg
rownames(rosmap_res) <- rosmap_res$cpg

###brain_shared_cpgs
brain_shared_cpg <- intersect(cbmap_res$cpg,rosmap_res$cpg)#4692
dat_cb <- cbmap_res[brain_shared_cpg,]
dat_ros <- rosmap_res[brain_shared_cpg,]
identical(dat_cb$cpg,dat_ros$cpg)

colnames(dat_cb)[2:7] <- paste0("cbmap_", colnames(dat_cb)[2:7])
colnames(dat_ros)[2:7] <- paste0("rosmap_", colnames(dat_ros)[2:7])
brain_share_table <- cbind(dat_cb,dat_ros)
brain_share_table <- brain_share_table[,-8]
write_xlsx(brain_share_table,'brain_shared_cpg.xlsx')

###brain_specific_cpgs
load("meta_meth_age_blood_60years.RData")
blood_res = res
rownames(blood_res) <- blood_res$cpg
blood_p_sig <- blood_res[blood_res$p_bacon<0.05,'cpg']
brain_speci_cpg <- setdiff(brain_shared_cpg,blood_p_sig)#3123

#filter monash
monash_data <- read.csv("Blood_DMPs.csv")
brain_speci_cpg <- setdiff(brain_speci_cpg,monash_data$MarkerName)

dat_cb <- cbmap_res[brain_speci_cpg,]
dat_ros <- rosmap_res[brain_speci_cpg,]
dat_eur <- blood_res[brain_speci_cpg,]
identical(dat_eur$cpg,dat_cb$cpg);identical(dat_cb$cpg,dat_ros$cpg)

colnames(dat_cb)[2:7] <- paste0("cbmap_", colnames(dat_cb)[2:7])
colnames(dat_ros)[2:7] <- paste0("rosmap_", colnames(dat_ros)[2:7])
colnames(dat_eur)[2:23] <- paste0("eur_", colnames(dat_eur)[2:23])
brain_speci_table <- cbind(dat_cb,dat_ros,dat_eur)
brain_speci_table <- brain_speci_table[,-c(8,15)]
write.xlsx(brain_speci_table,'brain_specific_cpg.xlsx')

###CBMAP_specific
load("ROSMAP_meth_age_noDementia.RData")
rosmap_res <- res
rownames(rosmap_res) <- rosmap_res$cpg

ros_p_sig <- rosmap_res[rosmap_res$p_bacon<0.05,'cpg']
eas_speci_cpg <- setdiff(cbmap_res$cpg,ros_p_sig)#2792

dat_cb <- cbmap_res[eas_speci_cpg,]
dat_ros <- rosmap_res[eas_speci_cpg,]
identical(dat_cb$cpg,dat_ros$cpg)

colnames(dat_cb)[2:7] <- paste0("cbmap_", colnames(dat_cb)[2:7])
colnames(dat_ros)[2:7] <- paste0("rosmap_", colnames(dat_ros)[2:7])
eas_speci_table <- cbind(dat_cb,dat_ros)
eas_speci_table <- eas_speci_table[,-8]
write_xlsx(eas_speci_table,'eas_specific_cpg.xlsx')



