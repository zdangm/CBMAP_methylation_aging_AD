library(data.table)
library(dplyr)
library(rhdf5)
library(lmerTest)
library(parallel)
ncore <- 40

pheno_df = read.table('/SATSA/pheno_df.txt',header=T, sep='\t')
celltype = read.table('/SATSA/celltype_with_clr.txt',header=T)
pheno_df = merge(pheno_df, celltype, by='row.names')
rownames(pheno_df) = pheno_df$Row.names
pheno_df = pheno_df[which(pheno_df$age > 60),] 

cg <- h5read("/SATSA/final_DNAm_invMdat.h5",'Mdat')
colnames(cg) = h5read("/SATSA/final_DNAm_invMdat.h5",'colnames')
rownames(cg) = h5read("/SATSA/final_DNAm_invMdat.h5",'rownames')

sample = intersect(rownames(cg), rownames(pheno_df)) 
cg <- cg[sample,]
pheno_df <- pheno_df[sample,]
pheno_df$msex <- as.factor(pheno_df$msex)
pheno_df$twin_pair   <- factor(pheno_df$twin_pair)
pheno_df$individualID <- factor(pheno_df$individualID)
identical(rownames(pheno_df), rownames(cg))

# Linear regression by cpgs Methylation
form <- DNAm ~ age + msex + B + NK + CD4T + CD8T + Mono + Neutro + Eosino + (1|twin_pair)

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

res_ewas <- do.call(rbind, res_list)
res_ewas$FDR <- p.adjust(res_ewas$p, method = "BH")
save(res_ewas, file='SATSA_meth_age.RData')