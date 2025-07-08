library(glmnet)
library(rhdf5)
library(data.table)
library(ggplot2)
library(readxl)
library(dplyr)

load('CBMAP_pheno.RData')
# methylation matrix
Mdat <- h5read("final_DNAm_Mdat.h5",'Mdat')
colnames(Mdat) = h5read("final_DNAm_Mdat.h5",'colnames')
rownames(Mdat) = h5read("final_DNAm_Mdat.h5",'rownames')
Bdat <- 2^Mdat / (2^Mdat + 1)

all_sample_info <- as.data.frame(read.csv("CBMAP_sample_info.csv",header=T)) 
selected_sample = all_sample_info[!is.na(all_sample_info[,'ADNC']),]        
selected_sample <- selected_sample %>% 
  filter(
    diag_ALS == 0 & 
      diag_SCZ == 0 & 
      diag_epilepsy == 0 &
      diag_others == 0 &
      OTHER == 0) 
case <- selected_sample$id[which(selected_sample$ADNC == 1)] 
hc <- selected_sample$id[which(selected_sample$ADNC == 0 &
                                 selected_sample$diag_dementia == 0 &
                                 selected_sample$diag_PD == 0 & 
                                 selected_sample$LBD == 0 &
                                 selected_sample$ARTAG_Subpial == 0 &
                                 selected_sample$ARTAG_Gray.matter == 0 &
                                 selected_sample$ARTAG_Perivascular == 0)] 


# Separate training samples and test samples
train = pheno_df[which(pheno_df$Sample_Name %chin% hc & pheno_df$bank %chin% c('zju','csu')),]
test = pheno_df[which((pheno_df$Sample_Name %chin% hc & pheno_df$bank == 'pumc') | (pheno_df$Sample_Name %chin% case)),]
train = pheno_df[which(pheno_df$Sample_Name %chin% hc),] 
test = pheno_df[which(pheno_df$Sample_Name %chin% case),] 
train = pheno_df[which(pheno_df$bank %chin% c('zju','csu')),]
test = pheno_df[which(pheno_df$bank == 'pumc'),]


# training
betasTrain<-Bdat[,match(train$Sample_Name, colnames(Bdat))]
identical(colnames(betasTrain),train$Sample_Name)
betasTrain<-t(betasTrain)
# testing
betasTest<-Bdat[,match(test$Sample_Name, colnames(Bdat))]
identical(colnames(betasTest),test$Sample_Name)
betasTest<-t(betasTest)

#Perform PCA and projections. Remove last PC.
PCA = prcomp(betasTrain,scale.=F)
TrainPCData = PCA$x[,1:(dim(PCA$x)[2]-1)]
TestPCData = predict(PCA,betasTest)[,1:(dim(PCA$x)[2]-1)]

#Select phenotype to be predicted. For example, here we are predicting age.
TrainAge = train$age
TestAge = test$age

#Train PC clock. Can test different models using different alpha and lambda parameters (see glmnet documentation)
cv = cv.glmnet(TrainPCData, TrainAge, nfolds=10,alpha=0.5, family="gaussian")
fit = glmnet(TrainPCData, TrainAge, family="gaussian", alpha=0.5, nlambda=100)

#Examine full model
plot(TrainAge,predict(fit,TrainPCData,s = cv$lambda.min),xlab = "Age",ylab = "Predicted Age", main = "Training")
cor(TrainAge,predict(fit,TrainPCData,s = cv$lambda.min))
plot(TestAge,predict(fit,TestPCData,s = cv$lambda.min),xlab = "Age",ylab = "Predicted Age", main = "Testing")
cor(TestAge,predict(fit,TestPCData,s = cv$lambda.min))

#Examine sparse model
plot(TrainAge,predict(fit,TrainPCData,s = cv$lambda.1se),xlab = "Age",ylab = "Predicted Age", main = "Training")
cor(TrainAge,predict(fit,TrainPCData,s = cv$lambda.1se))
plot(TestAge,predict(fit,TestPCData,s = cv$lambda.1se),xlab = "Age",ylab = "Predicted Age", main = "Testing")
cor(TestAge,predict(fit,TestPCData,s = cv$lambda.1se))

#Most likely your final model will only use a small subset of PCs. Thus you can compress your model:
CalcPCAge <- vector(mode = "list",length = 0)
temp = as.matrix(coef(cv,s = cv$lambda.min))
CalcPCAge$model = temp[temp!=0,][-1]
CalcPCAge$intercept = temp[1,1]
CalcPCAge$center = PCA$center
CalcPCAge$rotation = PCA$rotation[,names(CalcPCAge$model)]
CpGs = rownames(Bdat)

#And save your model:
save(CalcPCAge,CpGs,file = "CalcPCAge.RData")

#Then you can calculate your new clock in test data using compressed code, which will now take less time to load and calculate PCs:
load(file = "CalcPCAge.RData")
PCAgetest <- sweep(betasTest,2,CalcPCAge$center) %*% CalcPCAge$rotation %*% CalcPCAge$model + CalcPCAge$intercept
testres <- data.frame(ID=rownames(betasTest),PCAge=PCAgetest)
testres <- merge(testres,test[,c('Sample_Name','age')],by.x='ID',by.y='Sample_Name')
write.table(testres,file='PCAge_test.txt',row.names=F,col.names=T,quote=F)
PCAgetrain <- sweep(betasTrain,2,CalcPCAge$center) %*% CalcPCAge$rotation %*% CalcPCAge$model + CalcPCAge$intercept
trainres <- data.frame(ID=rownames(betasTrain),PCAge=PCAgetrain)
trainres <- merge(trainres,train[,c('Sample_Name','age')],by.x='ID',by.y='Sample_Name')
write.table(trainres,file='PCAge_train.txt',row.names=F,col.names=T,quote=F)
plotAge<-function(data,filename){
  ggplot(data=data,aes(x=age,y=PCAge)) +
    geom_abline(intercept=0,slope=1) +
    geom_point()+
    xlab("Chronological Age")+
    ylab("Predicted Age")+
    ggtitle("CBMAP PCAge") +
    theme_minimal() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          axis.title=element_text(size=16,face="bold"),
          legend.text=element_text(size=14,face="bold"),
          plot.title=element_text(size=16, face="bold.italic"),
          axis.text=element_text(size=14,face="bold"),
          strip.text.x=element_text(face="bold",size=16))+
    annotate("text", label = paste0("Corr = ",round(cor(data[,2],data[,3]),2)," RMSE = ",round(sqrt(mean((data[,2]-data[,3])^2)),2)), size = 2, x = min(data$age)+2, y = max(data$PCAge)-5,hjust=0)
  ggsave(paste0(filename,'.pdf'), height= 5, width=6)
}

plotAge(testres, 'PCAge_test')
plotAge(trainres, 'PCAge_train')