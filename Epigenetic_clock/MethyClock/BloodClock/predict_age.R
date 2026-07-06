library(rhdf5)
library(data.table)

setwd('/methylation/results/disease_meth_age/DNAm_clock/BloodClock/input/pumctest')
load('CBMAP_pheno.RData')
# methylation matrix
Mdat <- h5read("/methylation/data/CBMAP/DNAm_processed/final_DNAm_Mdat.h5",'Mdat')
colnames(Mdat) = h5read("/methylation/data/CBMAP/DNAm_processed/final_DNAm_Mdat.h5",'colnames')
rownames(Mdat) = h5read("/methylation/data/CBMAP/DNAm_processed/final_DNAm_Mdat.h5",'rownames')
Bdat <- 2^Mdat / (2^Mdat + 1)
Bdat <- t(Bdat)

pheno_df = pheno_df[which(pheno_df$bank == 'pumc'),]
sample = intersect(rownames(Bdat), rownames(pheno_df)) 
Bdat <- Bdat[sample,]
saveRDS(Bdat,file='betas.rds')
pheno_df <- pheno_df[sample,]
identical(rownames(pheno_df), rownames(Bdat))
real_age <- data.frame(ID=pheno_df$Sample_Name, age=pheno_df$age)
write.table(real_age,file='age.data.txt',quote=F,row.names=F,col.names=T)

##### Run these commands on the Slurm cluster
cd /methylation/script/disease_meth_age/DNAm_clock/BloodClock/DNAm-based-age-predictor-master
Rscript pred.R -i /methylation/results/disease_meth_age/DNAm_clock/BloodClock/input/pumctest/betas.rds -o /methylation/results/disease_meth_age/DNAm_clock/BloodClock/result/pumctest/result.txt -a /methylation/results/disease_meth_age/DNAm_clock/BloodClock/input/pumctest/age.data.txt


##### Script for plotting figure and calculating prediction accuracy #####
read.table(file='/methylation/results/disease_meth_age/DNAm_clock/BloodClock/result/pumctest/result.txt',stringsAsFactor=F,header=T) -> age
library(ggplot2)
library(reshape2)
colnames(age)<-c("ID","age.raw","Elastic Net","BLUP")
corv<-c(round(cor(age[,2],age[,3]),2),round(cor(age[,2],age[,4]),2))
rmse<-c(round(sqrt(mean((age[,2]-age[,3])^2)),2),round(sqrt(mean((age[,2]-age[,4])^2)),2))
age<-melt(age,measure.vars=c("Elastic Net","BLUP"))
df1 <- data.frame(variable=c('Elastic Net'),label=paste0("Corr = ",corv[1], " RMSE = ",rmse[1]))
df2 <- data.frame(variable=c('BLUP'),label=paste0("Corr = ",corv[2], " RMSE = ",rmse[2]))
ggplot(data=age,aes(x=age.raw,y=value))+
  facet_wrap(~variable)+
  geom_abline(intercept=0,slope=1)+
  geom_point()+
  xlab("Chronological Age")+
  ylab("Predicted Age")+
  theme(axis.title=element_text(size=12,face="bold"),
        legend.text=element_text(size=10,face="bold"),
        legend.title=element_text(face="bold"),
        axis.text=element_text(size=10,face="bold"),
        strip.text.x=element_text(face="bold")
        )+
  geom_text(data=df1, 
           aes(label = label), 
           size = 4, 
           x = min(age$age.raw)+1, 
           y = max(age$age.raw)-38,
           hjust=0
           )+
  geom_text(data=df2, 
            aes(label = label), 
            size = 4, 
            x = min(age$age.raw)+1, 
            y = max(age$age.raw)-38,
            hjust=0
  )
