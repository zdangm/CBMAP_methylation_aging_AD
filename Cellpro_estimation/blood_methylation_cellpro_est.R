library(data.table)
library(EpiDISH)
library(readxl)
library(bigmelon)
library(rhdf5)
# data(cent12CT450k.m)
data(centDHSbloodDMC.m)

# estimate celltype proportion
cg <- h5read("/GSE111629/BMIQ_Mdat.h5",'Mdat')
colnames(cg) = h5read("/GSE111629/BMIQ_Mdat.h5",'colnames')
rownames(cg) = h5read("/GSE111629/BMIQ_Mdat.h5",'rownames')
cg = M2Beta(cg)
BloodFrac.m <- epidish(beta.m = cg, ref.m = centDHSbloodDMC.m, method = "RPC")$estF
write.table(BloodFrac.m,file='/GSE111629/celltype.txt',row.names=T,col.names=T,quote=F,sep='\t')
boxplot(BloodFrac.m)

# Perform clr-transform
library(compositions)
ctp = read.table('/GSE111629/celltype.txt',sep='\t',header=T)
ctp<-as.matrix(ctp)
ctp[which(ctp <= 0)] <- 1e-3
ctp<-clr(ctp)
write.table(ctp,"/GSE111629/celltype_with_clr.txt",row.names = T,col.names = T,sep = "\t")

