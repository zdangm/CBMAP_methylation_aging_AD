library(methylGSA)
library(openxlsx)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

load('meta_hml4level.RData')
hml_cpg = meta_final_ordered_df$pVal.final[which(! is.na(meta_final_ordered_df$fdr))]
names(hml_cpg) = meta_final_ordered_df$cpg[which(! is.na(meta_final_ordered_df$fdr))]
load('meta_meth_braak.RData')
braak_cpg = meta_final_ordered_df$pVal.final[which(! is.na(meta_final_ordered_df$fdr))]
names(braak_cpg) = meta_final_ordered_df$cpg[which(! is.na(meta_final_ordered_df$fdr))]
load('meta_meth_Ascore.RData')
ascore_cpg = meta_final_ordered_df$pVal.final[which(! is.na(meta_final_ordered_df$fdr))]
names(ascore_cpg) = meta_final_ordered_df$cpg[which(! is.na(meta_final_ordered_df$fdr))]
cpg_list = list('hml'=hml_cpg, 'braak'=braak_cpg, 'ascore' = ascore_cpg)

# methylGSA::methylRRA
res = vector('list',6)
o = 0
for ( GS.type in c("GO","KEGG")) {
  for (j in 1:length(cpg_list)) {
    o = o + 1
    group = names(cpg_list)[j]
    sig.cpg = cpg_list[[j]]
    print(group)
    print(GS.type)
    
    res1 = methylRRA(cpg.pval = sig.cpg, array.type = "450K", method = "ORA", minsize = 5, maxsize = 500, GS.type = GS.type)
    res1.sig<-res1[res1$padj<0.05,]
    res[[o]] = res1.sig
    names(res)[o] = paste0(group,'_',GS.type)
  }
}
write.xlsx(res,'methylGSA_methylRRA_ORA.sig.xlsx')
