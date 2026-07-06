PreAge_plot = function(disease, residual_type) {
  residual = apply(res[,c("CBMAP_EN","CBMAP_PC","CorticalAge","PCBrainAge","MultiAge","BloodAge","HannumAge","PhenoAge")],2,function(x){
    model = lm(x~res$age.x+res$NeuN_pos+res$sex_male)
    res0 = residuals(model)
    res0 = scale(res0)
    if (residual_type == 'Age-scaled') {
      res = res0/res$age.x 
    } else if (residual_type == 'Unscaled') {
      res = res0
    }
    return(res)
  })
  residual = as.data.frame(residual)
  residual$ID = res$ID
  
  all_sample_info <- as.data.frame(read.csv("CBMAP_sample_info.csv",header=T))   
  selected_sample = all_sample_info[!is.na(all_sample_info[,'ADNC']),]
  selected_sample <- selected_sample %>% 
    filter(
      diag_ALS == 0 & 
        diag_SCZ == 0 & 
        diag_epilepsy == 0 &
        diag_others == 0 &
        OTHER == 0) 
  selected_sample$sex_male = factor(selected_sample$sex_male,levels=c(0,1))
  if (disease == 'HML1') {
    case1 <- unlist(selected_sample%>%filter(ADNC == 1)%>%select(id)) 
    contr <- unlist(selected_sample %>% filter(ADNC == 0 & LBD == 0 & 
                                                 ARTAG_Gray.matter ==0 & 
                                                 ARTAG_Perivascular == 0 & 
                                                 ARTAG_Subpial == 0 & 
                                                 diag_dementia == 0 & 
                                                 diag_PD == 0)
                    %>%select(id))
    selected_sample <- selected_sample %>%
      mutate(case = case_when(
        id %chin% case1 ~ 'AD',
        id %chin% contr ~ 'Not AD',
        TRUE ~ NA
      ))
    selected_sample <- selected_sample %>%
      mutate(y = case_when(
        id %chin% case1 ~ 1,
        id %chin% contr ~ 0,
        TRUE ~ NA_integer_
      ))
    selected_sample$case = factor(selected_sample$case, levels=c('Not AD', 'AD'))
  }
  if (disease == 'HML2') {
    case1 <- unlist(selected_sample%>%filter(ADNC == 1)%>%select(id)) 
    contr <- unlist(selected_sample %>% filter(ADNC == 0 & LBD == 0 & 
                                                 ARTAG_Gray.matter ==0 & 
                                                 ARTAG_Perivascular == 0 & 
                                                 ARTAG_Subpial == 0 & 
                                                 diag_dementia == 0 & 
                                                 diag_PD == 0)
                    %>%select(id))
    selected_sample <- selected_sample %>%
      mutate(case2 = case_when(
        id %chin% case1 ~ 'AD',
        id %chin% contr ~ 'Not AD',
        TRUE ~ NA
      ))
    selected_sample <- selected_sample %>%
      mutate(case = case_when(
        (case2 == 'AD' & ADNC.H) == 1 ~ 'ADNC.H',
        (case2 == 'AD' & ADNC.M) == 1 ~ 'ADNC.M',
        (case2 == 'AD' & ADNC.L) == 1 ~ 'ADNC.L',
        case2 == 'Not AD' ~ 'Not AD',
        TRUE ~ NA
      ))
    selected_sample$case = factor(selected_sample$case,levels=c('Not AD','ADNC.L','ADNC.M','ADNC.H'))
    selected_sample <- selected_sample %>%
      mutate(y = case_when(
        case == 'ADNC.H' ~ 3,
        case == 'ADNC.M' ~ 2,
        case == 'ADNC.L' ~ 1,
        case == 'Not AD' ~ 0,
        TRUE ~ NA_integer_
      ))
  }
  selected_sample <- selected_sample[!is.na(selected_sample$case),]
  residual_disease = merge(residual, selected_sample, by.x='ID',by.y='id')
  
  aov.mean<-apply(residual_disease[,c("CBMAP_EN","CBMAP_PC","CorticalAge","PCBrainAge","MultiAge","BloodAge","HannumAge","PhenoAge")],2,function(x){
    aggregate(x,by=list(type=residual_disease$case),FUN=mean)
  })
  aov.sd<-apply(residual_disease[,c("CBMAP_EN","CBMAP_PC","CorticalAge","PCBrainAge","MultiAge","BloodAge","HannumAge","PhenoAge")],2,function(x){
    aggregate(x,by=list(type=residual_disease$case),FUN=sd)
  })
  if (disease %chin% c('HML2','braak')) {
    case.lm<-apply(residual_disease[,c("CBMAP_EN","CBMAP_PC","CorticalAge","PCBrainAge","MultiAge","BloodAge","HannumAge","PhenoAge")],2,function(x){
      model<-lm(x~residual_disease$y)
      model<-summary(model)
      res = model$coefficients[2,c(1,2,4)]
      return(res)
    })
  } else {
    case.lm<-apply(residual_disease[,c("CBMAP_EN","CBMAP_PC","CorticalAge","PCBrainAge","MultiAge","BloodAge","HannumAge","PhenoAge")],2,function(x){
      model<-lm(x~residual_disease$y)
      model<-summary(model)
      res = model$coefficients[2,c(1,2,4)]
      return(res)
    })
  } 
  aov=vector("list", length = 8)
  for (i in 1:length(aov.mean)){
    aov[[i]]$dat <- merge(aov.mean[[i]], aov.sd[[i]], by='type')
    colnames(aov[[i]]$dat)<-c('type','mean','sd')
    aov[[i]]$title<-names(aov.mean)[i]
    aov[[i]]$disease<-disease
    aov[[i]]$p<-case.lm[3,i]
  }
  names(aov) = names(aov.mean)
  plot_fun(aov[[2]])
}

plot_fun = function(data) {
  ggplot(data=data$dat,aes(x=data$dat$type, y=data$dat$mean))+
    geom_bar(stat="identity",position="dodge",fill='#B7DAFE',color='black')+
    geom_errorbar(aes(ymax=data$dat$mean+data$dat$sd/2,ymin=data$dat$mean-data$dat$sd/2),position=position_dodge(0.9),width=0.15)+
    ggtitle(paste0('P=',format(data$p, scientific = TRUE, digits = 3))) +
    labs(x='ADNC Level',y="DNAm Age Acceleration")+
    theme_bw(base_size = 12) + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          plot.title = element_text(size=7)
    ) 
}

