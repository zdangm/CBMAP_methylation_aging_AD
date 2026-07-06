library(methylGSA)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(openxlsx)

load("CBMAP_meth_age_no_dementia.RData")
cpg.pval <- as.vector(res$pValue)
names(cpg.pval) <- res$cpg

GS_types <- c("GO", "KEGG")
res <- vector("list", length(GS_types))
names(res) <- GS_types

# methylGSA::methylRRA
for ( gs in GS_types) {
  cat("Processing", gs, "\n")
  
  res1 <- methylRRA(
    cpg.pval = cpg.pval,
    array.type = "EPIC",  
    method = "GSEA",
    minsize = 5,
    maxsize = 500,
    GS.type = gs
  )
  res1.sig <- res1[res1$padj < 0.05, ]
  res[[gs]] <- res1.sig
}
write.xlsx(res,'age_methylRRA_GSEA.sig.xlsx')