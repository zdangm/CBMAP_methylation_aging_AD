library(minfi)
library(data.table)
suppressMessages(library(IlluminaHumanMethylation450kanno.ilmn12.hg19))

data(Locations) 

# extract autosomal_probes
autosomal_probes <- rownames(Locations[Locations$chr %in% paste0("chr", 1:22), ])
write.table(autosomal_probes, file="/methylation/results/ctp_analysis/OSCA/450K_aut.probe", row.names=FALSE, col.names=FALSE, quote=FALSE)


# extract masked_probes
masked_probes <- rownames(Locations[!(Locations$chr %in% paste0("chr", 1:22)), ])
write.table(masked_probes, file="/methylation/results/ctp_analysis/OSCA/450K_mask.probe", row.names=FALSE, col.names=FALSE, quote=FALSE)


################# process 935k EPIC for CBMAP
anno <- fread('EPIC-8v2-0_A1.csv')
colnames(anno) <- as.character(anno[8,])
anno <- anno[-c(1:8),]
anno <- anno[!duplicated(anno$Name, fromLast = TRUE), ]

# extract autosomal_probes
autosomal_probes <- anno$Name[anno$CHR %in% paste0("chr", 1:22)]
write.table(autosomal_probes, file="/methylation/results/ctp_analysis/OSCA/935k_EPIC_aut.probe", row.names=FALSE, col.names=FALSE, quote=FALSE)


# extract masked_probes
masked_probes <- anno$Name[!(anno$CHR %in% paste0("chr", 1:22))]
write.table(masked_probes, file="/methylation/results/ctp_analysis/OSCA/935k_EPIC_mask.probe", row.names=FALSE, col.names=FALSE, quote=FALSE)
