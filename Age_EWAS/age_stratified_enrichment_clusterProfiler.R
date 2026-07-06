library(data.table)
library(dplyr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)
library(openxlsx)

cpg_gene_res = as.data.frame(fread("cis_eQTM_all_res_clr.txt",header=T))
bg_gene <- rownames(read.table("RNAseq_data_qn_comb.txt")) 

bg_gene_ensembl <- bg_gene
bg_gene_entrez <- bitr(bg_gene, 
                       fromType = "ENSEMBL", 
                       toType = "ENTREZID", 
                       OrgDb = org.Hs.eg.db)$ENTREZID

run_enrichment <- function(cpg_vec, prefix) {
  res <- cpg_gene_res %>%
    filter(cpg %in% cpg_vec) %>%
    group_by(cpg) %>%
    mutate(fdr = p.adjust(pvalue, method = "fdr")) %>%
    ungroup() %>%
    filter(fdr < 0.05)
  
  gene_ensembl <- unique(res$gene_id)
  
  if (length(gene_ensembl) == 0) {
    return(list(GO = NULL, KEGG = NULL))
  }
  
  # GO enrich
  ego <- enrichGO(
    gene          = gene_ensembl,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENSEMBL",
    ont           = "ALL",
    universe      = bg_gene_ensembl,
    pAdjustMethod = "BH",
    minGSSize     = 5,
    maxGSSize     = 5000,
    pvalueCutoff  = 1,
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
  ego_df <- as.data.frame(ego)
  ego_df <- ego_df[order(ego_df$qvalue), ]
  ego_df$geneID <- gsub("/", ", ", ego_df$geneID)
  write.xlsx(ego_df, file = paste0(prefix, "_GO.xlsx"))
  
  # KEGG enrich
  gene_entrez <- bitr(gene_ensembl,
                      fromType = "ENSEMBL",
                      toType = "ENTREZID",
                      OrgDb = org.Hs.eg.db)$ENTREZID
  
  if (length(gene_entrez) == 0) {
    kegg_df <- data.frame() 
  } else {
    kegg <- enrichKEGG(
      gene          = gene_entrez,
      keyType       = "kegg",
      organism      = "human",
      universe      = bg_gene_entrez,
      qvalueCutoff  = 0.05,
      pvalueCutoff  = 1
    )
    kegg_df <- as.data.frame(kegg)
    if (nrow(kegg_df) > 0) {
      kegg_df <- kegg_df[order(kegg_df$qvalue), ]
      rownames(kegg_df) <- NULL
    }
  }
  write.xlsx(kegg_df, file = paste0(prefix, "_KEGG.xlsx"))
  return(list(GO = ego_df, KEGG = kegg_df))
}

load("CBMAP_meth_age_no_dementia.RData")
up_cpgs <- res[res$p_bacon_fdr<0.05&res$Estimate>0,'cpg']
down_cpgs <- res[res$p_bacon_fdr<0.05&res$Estimate<0,'cpg']

res_up   <- run_enrichment(up_cpgs,   "age_upreg_cpgs")
res_down <- run_enrichment(down_cpgs, "age_downreg_cpgs")
