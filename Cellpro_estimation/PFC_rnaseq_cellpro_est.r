################## RNAseq deconvlution ###########


library(data.table)
library(stringr)
library(ggplot2)
library(sva)
library(edgeR)
library(MIND)

source("/tools/CIBERSORT/CIBERSORT.R")
load(paste0('/data/shared_data/annotate_file/gencode/gencode_v46_b38_annotation.RData'))   

########### read signature for cibersort
signature<-read.csv("/Signature_matrix_Darmanis.csv")
rownames(signature)<-signature$X
signature<-signature[,-1]

############## read metadata to get cbmap batch
metadata <- fread(paste0('RNAseq_metadata.txt')) 
########### only keep the sample with qc
final_sampelid<-read.table("final_sampleid.txt",header = F)
metadata <- metadata[metadata$sample_id %in% final_sampelid$V1,]

############ read cbmap raw counts
counts<-fread("raw_count")
rownames(counts)<-counts$gene_id

################ rename the cbmap first batch sampleid
colnames(counts)<-ifelse(str_detect(colnames(counts),"-"),paste0("20",colnames(counts)),colnames(counts))
colnames(counts)<-ifelse(str_detect(colnames(counts),"-"),str_replace(colnames(counts),"-","CBB"),colnames(counts))

########### only keep the sample with qc
counts<-as.data.frame(counts)
counts_keep<-counts[,colnames(counts) %in% metadata$sample_id]
rownames(counts_keep)<-counts$gene_id

###### keep the same sampleid order in metadata and counts_keep
metadata<-metadata[match(colnames(counts_keep),metadata$sample_id),]
unique(metadata$sample_id==colnames(counts_keep))
metadata<-as.data.frame(metadata)

########################### check the batch for cbmap
runPCA <- function(var) {
  
  p <- ggplot(as.data.frame(prcomp$x[,1:2])) +
    geom_point(aes(x = PC1, y = PC2, color = metadata[,var]), alpha = .6) + 
    theme_bw() +
    labs(color = var,
         title = var,
         x = paste("PC1 (", eigsP1, ")", sep = ""),
         y = paste("PC2 (", eigsP2, ")", sep = ""))
  
  return(p)
  
}


prcomp <- prcomp(t(counts_keep))
eigs <- prcomp$sdev^2
eigsP1 <- scales::percent(eigs[1] / sum(eigs))
eigsP2 <- scales::percent(eigs[2] / sum(eigs))


p1 <- runPCA('RIN')
p1
p2 <- runPCA('Bank')
p2
p3 <- runPCA("TestMethod")
p3
p4 <- runPCA('sequencingBatch')
p4

########################### corrected for batch. No covariates
counts_keep_combatseq <- ComBat_seq(counts = as.matrix(counts_keep), batch = metadata$TestMethod)
saveRDS(counts_keep_combatseq, file = "cbmap_qc_sample_counts_keep_combatseq.rds")

########### batch correction check
prcomp <- prcomp(t(counts_keep_combatseq))
eigs <- prcomp$sdev^2
eigsP1 <- scales::percent(eigs[1] / sum(eigs))
eigsP2 <- scales::percent(eigs[2] / sum(eigs))


p5 <- runPCA('RIN')
p5
p6 <- runPCA('Bank')
p6
p7 <- runPCA("TestMethod")
p7
p8 <- runPCA('sequencingBatch')
p8

################### analysis from correct batch to CIBERSORT input follow by devBrain
dge <- DGEList(counts = counts_keep_combatseq)
lib.size <- colSums(dge$counts)

########### normlize
cpm_matrix <- cpm(dge)
bulk.cpm <- as.data.frame(cpm_matrix)

############## get  signature geneid
sig_gene<-rownames(signature)

miss_sig_gene<-data.frame(gene=setdiff(sig_gene,gencode_v46_annotation$gene_name),gene_id="")
#[1] "MKL2"      "LRMP"      "ATPIF1"    "FAM65B"    "GPR98"     "LRRC16A"   "PPAP2B"    "GPR56"     "ZCCHC6"    "PTRF"      "H3F3C"    
#[12] "LHFP"      "SDPR"      "GPR116"    "KIAA0226L" "HMHA1"     "FAM105A"   "FYB"       "ZNF812"    "KIAA0247"  "ACPP"      "PQLC3"    
#[23] "HIST2H2BF" "EMR2"      "CLECL1"    "MB21D1"    "FAM46A"    "FAM175A"   "AGPAT9"  

miss_sig_gene$gene_id[which(miss_sig_gene$gene=="MKL2")]<-unique(gencode_v46_annotation$gene_id[which(gencode_v46_annotation$gene_name=="MRTFB")])
miss_sig_gene$gene_id[which(miss_sig_gene$gene=="LRMP")]<-unique(gencode_v46_annotation$gene_id[which(gencode_v46_annotation$gene_name=="IRAG2")])
miss_sig_gene$gene_id[which(miss_sig_gene$gene=="ATPIF1")]<-unique(gencode_v46_annotation$gene_id[which(gencode_v46_annotation$gene_name=="ATP5IF1")])
miss_sig_gene$gene_id[which(miss_sig_gene$gene=="FAM65B")]<-unique(gencode_v46_annotation$gene_id[which(gencode_v46_annotation$gene_name=="RIPOR2")])
miss_sig_gene$gene_id[which(miss_sig_gene$gene=="GPR98")]<-unique(gencode_v46_annotation$gene_id[which(gencode_v46_annotation$gene_name=="ADGRV1")])
miss_sig_gene$gene_id[which(miss_sig_gene$gene=="LRRC16A")]<-unique(gencode_v46_annotation$gene_id[which(gencode_v46_annotation$gene_name=="CARMIL1")])
miss_sig_gene$gene_id[which(miss_sig_gene$gene=="PPAP2B")]<-unique(gencode_v46_annotation$gene_id[which(gencode_v46_annotation$gene_name=="PLPP3")])
miss_sig_gene$gene_id[which(miss_sig_gene$gene=="GPR56")]<-unique(gencode_v46_annotation$gene_id[which(gencode_v46_annotation$gene_name=="ADGRG1")])
miss_sig_gene$gene_id[which(miss_sig_gene$gene=="ZCCHC6")]<-unique(gencode_v46_annotation$gene_id[which(gencode_v46_annotation$gene_name=="TUT7")])
miss_sig_gene$gene_id[which(miss_sig_gene$gene=="PTRF")]<-unique(gencode_v46_annotation$gene_id[which(gencode_v46_annotation$gene_name=="CAVIN1")])
miss_sig_gene$gene_id[which(miss_sig_gene$gene=="H3F3C")]<-unique(gencode_v46_annotation$gene_id[which(gencode_v46_annotation$gene_name=="H3-5")])
miss_sig_gene$gene_id[which(miss_sig_gene$gene=="LHFP")]<-unique(gencode_v46_annotation$gene_id[which(gencode_v46_annotation$gene_name=="LHFPL6")])
miss_sig_gene$gene_id[which(miss_sig_gene$gene=="SDPR")]<-unique(gencode_v46_annotation$gene_id[which(gencode_v46_annotation$gene_name=="CAVIN2")])
miss_sig_gene$gene_id[which(miss_sig_gene$gene=="GPR116")]<-unique(gencode_v46_annotation$gene_id[which(gencode_v46_annotation$gene_name=="ADGRF5")])
miss_sig_gene$gene_id[which(miss_sig_gene$gene=="KIAA0226L")]<-unique(gencode_v46_annotation$gene_id[which(gencode_v46_annotation$gene_name=="RUBCNL")])
miss_sig_gene$gene_id[which(miss_sig_gene$gene=="HMHA1")]<-unique(gencode_v46_annotation$gene_id[which(gencode_v46_annotation$gene_name=="ARHGAP45")])
miss_sig_gene$gene_id[which(miss_sig_gene$gene=="FAM105A")]<-unique(gencode_v46_annotation$gene_id[which(gencode_v46_annotation$gene_name=="OTULINL")])
miss_sig_gene$gene_id[which(miss_sig_gene$gene=="FYB")]<-unique(gencode_v46_annotation$gene_id[which(gencode_v46_annotation$gene_name=="FYB1")])
miss_sig_gene$gene_id[which(miss_sig_gene$gene=="ZNF812")]<-unique(gencode_v46_annotation$gene_id[which(gencode_v46_annotation$gene_name=="ZNF812P")])
miss_sig_gene$gene_id[which(miss_sig_gene$gene=="KIAA0247")]<-unique(gencode_v46_annotation$gene_id[which(gencode_v46_annotation$gene_name=="SUSD6")])
miss_sig_gene$gene_id[which(miss_sig_gene$gene=="ACPP")]<-unique(gencode_v46_annotation$gene_id[which(gencode_v46_annotation$gene_name=="ACP3")])
miss_sig_gene$gene_id[which(miss_sig_gene$gene=="PQLC3")]<-unique(gencode_v46_annotation$gene_id[which(gencode_v46_annotation$gene_name=="SLC66A3")])
miss_sig_gene$gene_id[which(miss_sig_gene$gene=="HIST2H2BF")]<-unique(gencode_v46_annotation$gene_id[which(gencode_v46_annotation$gene_name=="H2BC18")])
miss_sig_gene$gene_id[which(miss_sig_gene$gene=="EMR2")]<-unique(gencode_v46_annotation$gene_id[which(gencode_v46_annotation$gene_name=="ADGRE2")])
unique(gencode_v46_annotation$gene_id[which(gencode_v46_annotation$gene_name=="CLECL1P")]) 
#> unique(gencode_v46_annotation$gene_id[which(gencode_v46_annotation$gene_name=="CLECL1P")])
#[1] "ENSG00000293488.1"  "ENSG00000184293.10"
miss_sig_gene$gene_id[which(miss_sig_gene$gene=="CLECL1")]<-"ENSG00000184293.10"
miss_sig_gene$gene_id[which(miss_sig_gene$gene=="MB21D1")]<-unique(gencode_v46_annotation$gene_id[which(gencode_v46_annotation$gene_name=="CGAS")]) 
miss_sig_gene$gene_id[which(miss_sig_gene$gene=="FAM46A")]<-unique(gencode_v46_annotation$gene_id[which(gencode_v46_annotation$gene_name=="TENT5A")]) 
miss_sig_gene$gene_id[which(miss_sig_gene$gene=="FAM175A")]<-unique(gencode_v46_annotation$gene_id[which(gencode_v46_annotation$gene_name=="ABRAXAS1")]) 

############ combine sig_gene
sig_gene<-data.frame(gene=sig_gene)
sig_gene<-merge(sig_gene,gencode_v46_annotation,by.x="gene",by.y="gene_name")
sig_gene<-unique(sig_gene[,c("gene","gene_id")])
sig_gene<-rbind(sig_gene,miss_sig_gene)

########### split gene_id
split_count_geneid<-str_split(rownames(bulk.cpm),"\\.",simplify = T)[,1]
rownames(bulk.cpm)<-split_count_geneid

split_sig_geneid<-str_split(sig_gene$gene_id,"\\.",simplify = T)[,1]
sig_gene$gene_id<-split_sig_geneid

#################### only keep the signature gene in bulk.cpm
bulk.cpm_sig<-bulk.cpm[rownames(bulk.cpm) %in% sig_gene$gene_id,]

############ get genename for bulk.cpm_sig
rownames(bulk.cpm_sig)<-sig_gene$gene[match(rownames(bulk.cpm_sig),sig_gene$gene_id)]
sig<-signature

################### analysis from correct batch to CIBERSORT input follow by MIND
######## normlize
counts_keep_combatseq_cpm = apply(counts_keep_combatseq, 2, function(x) x/sum(x)*1e6)
split_count_geneid<-str_split(rownames(counts_keep_combatseq_cpm),"\\.",simplify = T)[,1]
rownames(counts_keep_combatseq_cpm)<-split_count_geneid

#################### only keep the signature gene in counts_keep_combatseq_cpm
counts_keep_combatseq_cpm_sig<-counts_keep_combatseq_cpm[rownames(counts_keep_combatseq_cpm) %in% sig_gene$gene_id,]

################# get genename for counts_keep_combatseq_cpm_sig
rownames(counts_keep_combatseq_cpm_sig)<-sig_gene$gene[match(rownames(counts_keep_combatseq_cpm_sig),sig_gene$gene_id)]

###################### run CIBERSORT
cbprop_test_MIND <- CIBERSORT(sig, counts_keep_combatseq_cpm_sig, perm=0, QN=TRUE, absolute=FALSE, abs_method='sig.score')

write.table(cbprop_test_MIND,"/CIBERSORT/CIBERSORT_prop_MIND.txt",row.names = T,col.names = T)
