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

trafo= function(x,adult.age=20) { x=(x+1)/(1+adult.age); y=ifelse(x<=1, log( x),x-1);y }
test$AgeT=trafo(test$age)
train$AgeT=trafo(train$age)
anti.trafo<-function(x,adult.age=20) { ifelse(x<0, (1+adult.age)*exp(x)-1, (1+adult.age)*x+adult.age) }


# training
betasTrain<-Bdat[,match(train$Sample_Name, colnames(Bdat))]
identical(colnames(betasTrain),train$Sample_Name)
# testing
betasTest<-Bdat[,match(test$Sample_Name, colnames(Bdat))]
identical(colnames(betasTest),test$Sample_Name)


# use 10 fold cross validation to estimate the lambda parameter 
# in the training data
alpha<-0.5
glmnet.Training.CV = cv.glmnet(t(betasTrain), train$AgeT,  nfolds=10,alpha=alpha,family="gaussian") 
# The definition of the lambda parameter:
lambda.glmnet.Training = glmnet.Training.CV$lambda.min
# Fit the elastic net predictor to the training data
glmnet.Training = glmnet(t(betasTrain),train$AgeT, family="gaussian", alpha=0.5, nlambda=100)
# Arrive at an estimate of of DNAmAge
DNAmAgeBrainTraining=predict(glmnet.Training,t(betasTrain),type="response",s=lambda.glmnet.Training)
DNAmAgeBrainTesting=predict(glmnet.Training,t(betasTest),type="response",s=lambda.glmnet.Training)


# transform ages back 
DNAmAgeBrainTraining[,1]<-anti.trafo(DNAmAgeBrainTraining[,1])
colnames(DNAmAgeBrainTraining) <- c('PreAge')
DNAmAgeBrainTesting[,1]<-anti.trafo(DNAmAgeBrainTesting[,1])
colnames(DNAmAgeBrainTesting) <- c('PreAge')
PreAgeTrain <- merge(DNAmAgeBrainTraining, pheno_df[,c('Sample_Name','age')], by.x='row.names', by.y='Sample_Name')
colnames(PreAgeTrain) <- c('Sample_Name','brainpred','Age')
PreAgeTest <- merge(DNAmAgeBrainTesting, pheno_df[,c('Sample_Name','age')], by.x='row.names', by.y='Sample_Name')
colnames(PreAgeTest) <- c('Sample_Name','brainpred','Age')

write.table(PreAgeTrain,"PreAgeTrain.txt", row.names=F,col.names=T,quote=F)
write.table(PreAgeTest,"PreAgeTest.txt", row.names=F,col.names=T,quote=F)
save(glmnet.Training,lambda.glmnet.Training,  file = "CBMAPclock.rda")

tmp_coeffs <- coef(glmnet.Training.CV, s = "lambda.min")
myCoeff<-data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)
write.table(myCoeff,"CBMAPclockCoefs.txt", row.names=F, col.names=T, quote=F)

plotAge<-function(data,filename){
  ggplot(data=data,aes(x=Age,y=brainpred)) +
    geom_abline(intercept=0,slope=1) +
    geom_point()+
    xlab("Chronological Age")+
    ylab("Predicted Age")+
    ggtitle("CBMAP Clock") +
    theme_minimal() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          axis.title=element_text(size=16,face="bold"),
          legend.text=element_text(size=14,face="bold"),
          plot.title=element_text(size=16, face="bold.italic"),
          axis.text=element_text(size=14,face="bold"),
          strip.text.x=element_text(face="bold",size=16))+
    annotate("text", label = paste0("Corr = ",round(cor(data[,2],data[,3]),2)," RMSE = ",round(sqrt(mean((data[,2]-data[,3])^2)),2)), size = 2, x = min(data$Age)+2, y = max(data$brainpred)-5,hjust=0)
  ggsave(paste0(filename,'.pdf'), height= 5, width=6)
}

plotAge(PreAgeTrain, 'PreAgeTrain')
plotAge(PreAgeTest, 'PreAgeTest')