library(data.table)
library(stringr)
library(dplyr)
library(rhdf5)

cg <- h5read("/data/CBMAP_methylaition/final_DNAm_invMdat.h5",'Mdat')
colnames(cg) = h5read("/data/CBMAP_methylaition/final_DNAm_invMdat.h5",'colnames')
rownames(cg) = h5read("/data/CBMAP_methylaition/final_DNAm_invMdat.h5",'rownames')
cg = cg[,grep("cg",colnames(cg))]
anno_935k = read.csv('EPIC-8v2-0_A1.csv',header=T)
anno_935k = anno_935k[-c(1:6),]
colnames(anno_935k) = anno_935k[1,]
anno_935k = anno_935k[-1,]
anno_935k = anno_935k[!duplicated(anno_935k$Name), ]
cpg = intersect(colnames(cg),anno_935k$Name)
anno_935k = anno_935k[which(anno_935k$Name %chin% cpg),c('Name','CHR','MAPINFO')]
anno_935k$MAPINFO = as.integer(anno_935k$MAPINFO)

exp = read.table('/CBMAP_RNAseq_protein_coding/RNAseq_process/RNAseq_data_qn_comb.txt',check.names = F)
dim(exp)
#[1] 15621   933
load('gencode_v46_b38_annotation.RData')
gene_pos = gencode_v46_annotation[which(gencode_v46_annotation$type == 'gene'),c('gene_id','chr','start','end','gene_name')]
gene_pos$gene_id = unlist(str_split(gene_pos$gene_id,"[.]",simplify=T))[,1]
gene_pos = gene_pos[which(gene_pos$gene_id %chin% rownames(exp)),]
cg = cg[,cpg]
exp = t(exp)

# exp PCA
exp_pca = prcomp(exp, center=T, scale=T)
exp_pca_final = as.data.frame(exp_pca$x)
fwrite(exp_pca_final[,1:10], row.names=T, col.names=F, quote=F, sep="\t",
       file='/eQTM/exp_pca_10.txt')
summary(exp_pca)
pdf('/eQTM/exp_pca.pdf')
p=screeplot(exp_pca, type='line', main=NULL,lwd=2)
plot(p)
dev.off()
# DNAm PCA
cg_pca = prcomp(cg, center=T, scale=T)
cg_pca_final = as.data.frame(cg_pca$x)
fwrite(cg_pca_final[,1:10], row.names=T, col.names=F, quote=F, sep="\t",
       file='/eQTM/DNAm_pca_10.txt')
summary(cg_pca)
pdf('/eQTM/dnam_pca.pdf')
p=screeplot(cg_pca, type='line', main=NULL,lwd=2)
plot(p)
dev.off()

# extract cis-CpGs
anno_935k = as.data.table()
setkey(gene_pos, chr)
setkey(anno_935k, CHR)
gene_cpg_list <- lapply(unique(gene_pos$chr), function(chr0) {
  genes_chr <- gene_pos[which(gene_pos$chr == chr0),]
  cpgs_chr <- anno_935k[CHR == chr0]

  lapply(1:nrow(genes_chr), function(i) {
    gene <- genes_chr[i,]
    nearby_cpgs <- cpgs_chr[
      abs(MAPINFO - gene$start) < 1e6 | abs(MAPINFO - gene$end) < 1e6,
      .(cpg_id = Name)
    ]
    if (nrow(nearby_cpgs) > 0) {
      return(list(gene_id = gene$gene_id, cis_cpgs = nearby_cpgs$cpg_id))
    } else {
      return(NULL)
    }
  })
})
gene_cpg_list <- unlist(gene_cpg_list, recursive = FALSE)
gene_cpg_list <- gene_cpg_list[!sapply(gene_cpg_list, is.null)]
save(gene_cpg_list,file='/eQTM/cis-cpg.RData')

