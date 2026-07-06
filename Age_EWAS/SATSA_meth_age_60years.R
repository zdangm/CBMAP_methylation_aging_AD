library(data.table)
library(dplyr)
library(rhdf5)
library(lmerTest)
library(parallel)
ncore <- 40
setwd('/methylation/methy_AD_revision2/results/EWAS')

pheno_df = read.table('/methylation/data/SATSA/pheno_df.txt',header=T, sep='\t')
celltype = read.table('/methylation/data/SATSA/celltype_with_clr.txt',header=T)
pheno_df = merge(pheno_df, celltype, by='row.names')
rownames(pheno_df) = pheno_df$Row.names
pheno_df = pheno_df[which(pheno_df$age > 60),] # Comment out this line if using all samples

cg <- h5read("/methylation/data/SATSA/final_DNAm_invMdat.h5",'Mdat')
colnames(cg) = h5read("/methylation/data/SATSA/final_DNAm_invMdat.h5",'colnames')
rownames(cg) = h5read("/methylation/data/SATSA/final_DNAm_invMdat.h5",'rownames')
cg = cg[,grep("cg",colnames(cg))] 

sample = intersect(rownames(cg), rownames(pheno_df)) 
cg <- cg[sample,]
pheno_df <- pheno_df[sample,]
pheno_df$msex <- as.factor(pheno_df$msex)
pheno_df$twin_pair   <- factor(pheno_df$twin_pair)
pheno_df$individualID <- factor(pheno_df$individualID)
identical(rownames(pheno_df), rownames(cg))
pheno_df <- pheno_df %>% select(age, msex, B, NK, CD4T, CD8T, Mono, Neutro, twin_pair, individualID)
apply(pheno_df, 2, function(x){sum(is.na(x))})
# age         msex            B           NK         CD4T         CD8T         Mono 
#   0            0            0            0            0            0            0 
# Neutro    twin_pair individualID 
#      0            0            0

# Linear regression by cpgs Methylation
form <- DNAm ~ age + msex + B + NK + CD4T + CD8T + Mono + Neutro + (1|twin_pair / individualID)

res_list <- vector("list", ncol(cg))
names(res_list) <- colnames(cg)

res_list <- mclapply(
  X = seq_len(ncol(cg)),
  FUN = function(i) {
    dat <- pheno_df
    dat$DNAm <- cg[, i]
    
    fit <- lmer(form, data = dat, REML = FALSE)
    sm  <- summary(fit)$coefficients
    
    data.frame(
      CpG  = colnames(cg)[i],
      beta = sm["age", "Estimate"],
      se   = sm["age", "Std. Error"],
      p    = sm["age", "Pr(>|t|)"]
    )
  },
  mc.cores = ncore
)

res <- do.call(rbind, res_list)
res$FDR <- p.adjust(res$p, method = "BH")
save(res, file='SATSA_meth_age_60years.RData')
