PhenoClock<-function(betas, ## betas = betas matrix (rownames=cpgs, colnames=IDs)
                        pheno, ##  pheno file = file which contains IDs that match betas IDs and  contains actual Age col 
                        dir, ## directory where the coeffecients are saved 
                        IDcol, ## ID column which matches Betas IDs for your samples
                        Agecol){  ## Age column 
  
  ## subset pheno                    
  pheno<-pheno[c(IDcol,Agecol)]
  colnames(pheno)<-c("ID","Age")                 
  
  ## check ggplot2 is loaded - if not, load, if not installd then install
  pkgTest <- function(x)
  {
    if (!require(x,character.only = TRUE))
    {
      install.packages(x,dep=TRUE)
      if(!require(x,character.only = TRUE)) stop("Package not found")
    }
  }
  
  pkgTest("ggplot2")
  
  ## ensure colnames are ID and Age 
  colnames(pheno)<-c("ID", "Age")
  
  
  #########################################################
  ############# read in pheno clock coeffs #############
  #########################################################
  
  braincoef<-read.csv(paste0(dir,"aging-10-101414-s002.csv"),stringsAsFactor=F,header=T)[,c('CpG','Weight')]
  colnames(braincoef) <- c('probe','coef')
  braincoef$probe<-as.character(braincoef$probe)
  
  #########################################################
  ### find the overlap between the probes and your data ###
  #########################################################
  
  overlap<-braincoef[which(braincoef$probe %in% rownames(betas)),]
  betas <- betas[overlap$probe,]
  rownames(braincoef) <- braincoef$probe
  braincoef <- braincoef[overlap$probe,]
  
  
  
    
  ##############################################################################################
  ############# Age prediciton! - weighted sum of coefficients plus the intercept ##############
  ##############################################################################################
    
  braincoef<-braincoef[match(rownames(betas), braincoef$probe),]
  brainpred<-braincoef$coef%*%betas+60.664000
    
  ##############################################################################################
  ### anti transform the results (accounting for logarithmic relationship in ages 0-20)      ###
  ### Same as Horvath's: Horvath, S. (2013).                                                 ###
  ### DNA methylation age of human tissues and cell types. Genome biology, 14(10), 3156.     ###                                                              #######
  ##############################################################################################
    
  #anti.trafo<-function(x,adult.age=20) { ifelse(x<0, (1+adult.age)*exp(x)-1, (1+adult.age)*x+adult.age) }
  #brainpred<-anti.trafo(brainpred)
    
  #################################################
  #############  Save brain predictions ###########
  #################################################
  pheno<-pheno[match(colnames(betas), pheno$ID),]
  pheno$brainpred<-as.numeric(brainpred)
  pheno$Age<-as.numeric(pheno$Age)
  write.csv(pheno, "PhenoPred.csv")
    
  
  
  ############################################     
  ### get some stats for the predictions #####
  ############################################
  compareStats<-function(data){
    colnames(data)<-c("ID","Age","Pheno Clock")
    data<-data[complete.cases(data$Age),]
    corr<-round(cor(data[,2],data[,3]),2) ## correlation - pearsons r
    RMSE<-function(actualAge,estimatedAge){ ## root mean squared error (years)
      sqrt(mean((actualAge-estimatedAge)^2))
    }   
    rmse<-round(RMSE(data[,2],data[,3]),2)
    mad<-round(median(abs(data[,2]-data[,3])),2) ## mean absoluate deviation (years)
    stats<-matrix(ncol=1,nrow=3)
    colnames(stats)<-c("Pheno Clock")
    rownames(stats)<-c("Correlation (r)","RMSE (years)","MAD (years)")
    stats[1,]<-corr
    stats[2,]<-rmse
    stats[3,]<-mad
    print(stats)
    write.csv(stats, "accuracy_statistics.csv")
  }
  corr<-round(cor(pheno[,2],pheno[,3]),2) ## correlation - pearsons r
  RMSE<-function(actualAge,estimatedAge){ ## root mean squared error (years)
    sqrt(mean((actualAge-estimatedAge)^2))
  }   
  rmse<-round(RMSE(pheno[,2],pheno[,3]),2)
  
  compareStats(pheno)
  
  
  #####################
  ### make a plot #####
  #####################
  
  ## uses ggplot
  plotAge<-function(data){
    ggplot(data=data,aes(x=Age,y=brainpred)) +
      geom_abline(intercept=0,slope=1) +
      geom_point()+
      xlab("Chronological Age")+
      ylab("Predicted Age")+
      ggtitle("Pheno Clock") +
      theme_minimal() + 
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
            axis.title=element_text(size=16,face="bold"),
            legend.text=element_text(size=14,face="bold"),
            plot.title=element_text(size=16, face="bold.italic"),
            axis.text=element_text(size=14,face="bold"),
            strip.text.x=element_text(face="bold",size=16))+
      annotate("text", label = paste0("Corr = ",corr," RMSE = ",rmse), size = 2, x = min(data$Age)+2, y = max(data$brainpred)-5,hjust=0)
    
    ggsave("PhenoClockplot.pdf", height= 5, width=6)
  }
  
  plotAge(pheno)
  
  print("Your Pheno DNAm ages are saved in the file 'PhenoPred.csv'")
  print("Accuracy statistics comparing DNAm age and chronological age are saved in the file 'accuracy_statistics.csv'")
  print("Plot of Pheno DNAm age against age is saved in the file 'PhenoClockplot.pdf'")
  
}
