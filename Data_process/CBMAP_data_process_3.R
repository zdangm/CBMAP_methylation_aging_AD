##### inverse normalization transform and merge sh and zk DNAm #####
library(rhdf5)
zk = h5read('/methylation/data/CBMAP/DNAm_processed/zk_DNAm_processed/BMIQ_Mdat_combat.h5', 'Mdat')
colnames(zk) = h5read('/methylation/data/CBMAP/DNAm_processed/zk_DNAm_processed/BMIQ_Mdat_combat.h5', 'colnames')
rownames(zk) = h5read('/methylation/data/CBMAP/DNAm_processed/zk_DNAm_processed/BMIQ_Mdat_combat.h5', 'rownames')

sh = h5read('/methylation/data/CBMAP/DNAm_processed/sh_DNAm_processed/BMIQ_Mdat_combat.h5', 'Mdat')
colnames(sh) = h5read('/methylation/data/CBMAP/DNAm_processed/sh_DNAm_processed/BMIQ_Mdat_combat.h5', 'colnames')
rownames(sh) = h5read('/methylation/data/CBMAP/DNAm_processed/sh_DNAm_processed/BMIQ_Mdat_combat.h5', 'rownames')

cg = intersect(rownames(zk), rownames(sh))
zk = zk[cg,]
sh = sh[cg,]
final_DNAm = cbind(zk, sh)
rm(zk);rm(sh)

pheno_info = read.csv(file='/methylation/data/CBMAP/DNAm_processed/phenotype.csv', header=T)
colnames(final_DNAm) = pheno_info[match(colnames(final_DNAm), pheno_info$barcode_id), 'Sample_Name']

h5write(final_DNAm, '/methylation/data/CBMAP/DNAm_processed/final_DNAm_Mdat.h5', "Mdat")
h5write(rownames(final_DNAm), '/methylation/data/CBMAP/DNAm_processed/final_DNAm_Mdat.h5', "rownames")
h5write(colnames(final_DNAm), '/methylation/data/CBMAP/DNAm_processed/final_DNAm_Mdat.h5', "colnames")
InvNorm <- function(x) {
  return(qnorm((rank(x, na.last="keep")-3/8)/(sum(!is.na(x))+1/4)))
}
inverse_DNAm = apply(final_DNAm, 1, InvNorm)
h5write(inverse_DNAm, '/methylation/data/CBMAP/final_DNAm_invMdat.h5', "Mdat")
h5write(rownames(inverse_DNAm), '/methylation/data/CBMAP/final_DNAm_invMdat.h5', "rownames")
h5write(colnames(inverse_DNAm), '/methylation/data/CBMAP/final_DNAm_invMdat.h5', "colnames")
