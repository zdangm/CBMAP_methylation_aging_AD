#adnc
load("/meta_hml4level.RData")
meta_final_ordered_df <- meta_final_ordered_df[!is.na(meta_final_ordered_df$fdr),]
adnc_cpg <- meta_final_ordered_df[meta_final_ordered_df$fdr<0.05,]

adnc_up_cpg <- adnc_cpg[adnc_cpg$CBMAP_estimate>0,]
adnc_down_cpg <- adnc_cpg[adnc_cpg$CBMAP_estimate<0,]  

#braak
load("meta_meth_braak.RData")
meta_final_ordered_df <- meta_final_ordered_df[!is.na(meta_final_ordered_df$fdr),]
braak_cpg <- meta_final_ordered_df[meta_final_ordered_df$fdr<0.05,]

braak_up_cpg <- braak_cpg[braak_cpg$CBMAP_estimate>0,]
braak_down_cpg <- braak_cpg[braak_cpg$CBMAP_estimate<0,]

#ascore
load("meta_meth_Ascore.RData")
meta_final_ordered_df <- meta_final_ordered_df[!is.na(meta_final_ordered_df$fdr),]
ascore_cpg <- meta_final_ordered_df[meta_final_ordered_df$fdr<0.05,]

ascore_up_cpg <- ascore_cpg[ascore_cpg$CBMAP_estimate>0,]
ascore_down_cpg <- ascore_cpg[ascore_cpg$CBMAP_estimate<0,]

cpg_sets <- list(
  adnc_up = adnc_up_cpg$cpg,
  adnc_down = adnc_down_cpg$cpg,
  braak_up = braak_up_cpg$cpg,
  braak_down = braak_down_cpg$cpg,
  ascore_up = ascore_up_cpg$cpg,
  ascore_down = ascore_down_cpg$cpg
)

bg_gene <- read.table("/CBMAP_RNAseq_protein_coding/RNAseq_process/RNAseq_data_qn_comb.txt")
bg_gene <- rownames(bg_gene)
bg_gene_entrez = bitr(bg_gene,fromType="ENSEMBL",toType="ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID


run_enrichment <- function(cpg_vec, prefix){
  cpg_gene_res <- as.data.frame(fread("/eQTM/cis_eQTM_all_res_clr.txt"))
  cpg_gene_res <- cpg_gene_res[cpg_gene_res$cpg %in% cpg_vec,]
  cpg_gene_res <- cpg_gene_res %>%
    group_by(cpg) %>%
    mutate(fdr = p.adjust(pvalue, method = "fdr")) %>%
    ungroup()
  cpg_gene_res <- cpg_gene_res[cpg_gene_res$fdr < 0.05,]
  gene_ens <- unique(cpg_gene_res$gene_id)
  
  gene_entrez <- bitr(
    gene_ens,
    fromType = "ENSEMBL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )$ENTREZID
  gene_entrez <- unique(na.omit(gene_entrez))
  
  ego <- enrichGO(
    gene = gene_ens,
    OrgDb = org.Hs.eg.db,
    keyType = "ENSEMBL",
    ont = "ALL",
    universe = bg_gene,
    pAdjustMethod = "BH",
    minGSSize = 5,
    maxGSSize = 5000,
    pvalueCutoff = 1,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  
  ego_df <- as.data.frame(ego)
  if(nrow(ego_df) > 0){
    ego_df <- ego_df[order(ego_df$qvalue), ]
    ego_df$geneID <- gsub("/", ", ", ego_df$geneID)
  }
  
  ekegg <- enrichKEGG(
    gene = gene_entrez,
    organism = "hsa",
    universe = bg_gene_entrez,
    qvalueCutoff = 0.05,
    pvalueCutoff = 1
  )
  
  ekegg_df <- as.data.frame(ekegg)
  if(nrow(ekegg_df) > 0){
    ekegg_df <- ekegg_df[order(ekegg_df$qvalue), ]
  }
  
  out_file <- paste0("/methylation/methy_AD_revision2/results/meta_enrichment_i/direction_stratified/",prefix, "_GO_KEGG.xlsx")
  write.xlsx(list(GO = ego_df, KEGG = ekegg_df),file = out_file)
  
  message("done: ", prefix)
  return(list(GO = ego_df, KEGG = ekegg_df))
}

results <- lapply(names(cpg_sets), function(x){
  run_enrichment(cpg_sets[[x]], x)
})
names(results) <- names(cpg_sets)


