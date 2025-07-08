library(bnlearn)
library(Rgraphviz)
library(rhdf5)
library(data.table)
library(dplyr)
library(tidyr)
library(gRain)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(readxl)


### clinical info
load('CBMAP_pheno.RData')
all_sample_info <- as.data.frame(read.csv("CBMAP_sample_info.csv",header=T)) 
selected_sample = all_sample_info[!is.na(all_sample_info[,'ADNC']),]
selected_sample <- selected_sample %>% 
  filter(
    diag_ALS == 0 & 
      diag_SCZ == 0 & 
      diag_epilepsy == 0 &
      diag_others == 0 &
      OTHER == 0) 
complete_rows <- complete.cases(selected_sample[, c("ADNC", "ADNC.L", "ADNC.M", "ADNC.H")])
selected_sample <- selected_sample[complete_rows, ]
selected_sample$marker <- ifelse(selected_sample$ADNC == 0, 0,
                                 ifelse(selected_sample$ADNC.L == 1, 1,
                                        ifelse(selected_sample$ADNC.M == 1, 2,
                                               ifelse(selected_sample$ADNC.H == 1, 3, NA))))
exclued_id <- selected_sample[which(selected_sample$marker == 0 & (selected_sample$LBD == 1 | selected_sample$ARTAG_Gray.matter == 1 |
                                                                     selected_sample$ARTAG_Perivascular == 1 | selected_sample$ARTAG_Subpial == 1 | selected_sample$diag_PD == 1 | selected_sample$diag_dementia == 1)),'id']
selected_sample <- selected_sample[which(! selected_sample$id %chin% exclued_id),] #955
selected_sample$AD <- selected_sample$marker

# selected_sample$AD <- selected_sample$ADNC
bioAge1 <- read.table('hc_train_test/PCAge_train.txt',header=T)
bioAge2 <- read.table('hc_train_test/PCAge_test.txt',header=T)
bioAge <- rbind(bioAge1, bioAge2)
rm(bioAge1);rm(bioAge2)
colnames(bioAge)[1:2] = c('Sample_Name','brainpred')
dat = merge(selected_sample, bioAge[,c('Sample_Name','brainpred')], by.x='id', by.y='Sample_Name')
dat = dat[,c('brainpred','CVD_CAA', 'LATE','ARTAG','AD')]
colnames(dat)[1] = 'BioAge'
dat = dat[-which(dat$BioAge>90),]
dat$BioAge <- cut(dat$BioAge, breaks = c(-Inf, 50, 60, 70, 80, 90), labels = c("<51", "[51,60]", "[61,70]", "[71,80]", "[81,90]"))
dat <- data.frame(lapply(dat, as.factor))
dat <- na.omit(dat)
arc.set = matrix(c("BioAge", "CVD_CAA",
                   "BioAge", "LATE",
                   "BioAge", "ARTAG",
                   "BioAge", "AD",
                   "CVD_CAA","LATE",
                   "CVD_CAA","ARTAG",
                   "LATE","AD",
                   "ARTAG","AD",
                   "CVD_CAA","AD"),
                 byrow=T, ncol=2,
                 dimnames=list(NULL,c("from","to")))
dag = empty.graph(c('BioAge','CVD_CAA','LATE','ARTAG', 'AD'))
arcs(dag) = arc.set
graphviz.plot(dag)
dag_hc <- hc(dat,start=dag)
graphviz.plot(dag_hc)

mydat <- dat[,nodes(dag_hc)]
bn.mle <- bn.fit(dag_hc, data=mydat, method='mle')
graphviz.chart(bn.mle, grid = TRUE, 
               main = "Original BN", 
               scale = c(0.75, 1.3), 
               col = "black", 
               bg = "white",
               bar.col = 'lightskyblue')
junction <- compile(as.grain(bn.mle))
jlate <- setEvidence(junction, nodes = "LATE", states = "1")
graphviz.chart(as.bn.fit(jlate, including.evidence = TRUE), grid = TRUE,
               bar.col = c(BioAge = "grey", CVD_CAA = "grey",AD = "lightskyblue",ARTAG = "grey",LATE = "lightskyblue"),
               strip.bg = c(BioAge = "white", CVD_CAA = "white", AD = "white",ARTAG = "white",LATE = "grey"),
               main = "BN with LATE Evidence")
jcvd <- setEvidence(junction, nodes = "CVD_CAA", states = "1")
graphviz.chart(as.bn.fit(jcvd, including.evidence = TRUE), grid = TRUE,
               bar.col = c(BioAge = "grey", CVD_CAA = "lightskyblue",AD = "lightskyblue",ARTAG = "grey",LATE = "grey"),
               strip.bg = c(BioAge = "white", CVD_CAA = "grey", AD = "white",ARTAG = "white",LATE = "white"),
               main = "BN with CVD_CAA Evidence")

bn.fit.barchart(bn.mle$AD, main = "Probability of ADNC conditional on CVD_CAA and LATE",
                xlab = "Pr(ADNC | CVD_CAA, LATE)", ylab = "")