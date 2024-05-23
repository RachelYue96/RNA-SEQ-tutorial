library(DESeq2)
library(ggplot2)
library(pheatmap)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(dplyr)
library(gplots)
library(clusterProfiler)
library(pathview)
library('biomaRt')
library(ggfortify)
library(stringr)
library(enrichplot)
library(GOplot)
library(DOSE)
library(ggnewscale)
library(topGO)
library(circlize)
library(ComplexHeatmap)

# Set the working directory
wkdir <- "/share/result/sequencer/salus/video_example/RNA-seq/"
setwd(wkdir)

# Set the prefix for each output file name
outputPrefix <- "results/pipe_plots/"
dir.create(paste0(wkdir, outputPrefix))

countdata <- read.table("results/RNA_counts/allgene_count_matrix.txt", header=TRUE, row.names=1, check.names=FALSE, sep=",")
countdata <- as.matrix(countdata)

# Create a DataFrame that describes the condition for each sample
samples <- c('alt-1', 'alt-2', 'alt-3', 'con-1', 'con-2', 'con-3')
conditions <- c('Treated', 'Treated', 'Treated', 'Control', 'Control', 'Control')
coldata <- DataFrame(condition = factor(conditions), row.names = samples)

#guts
ddsHTSeq <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~ condition)
dds <- DESeq(ddsHTSeq)

# Normalize the data
# dds1 <- DESeq(dds,betaPrior = T)

res <- results(dds,contrast = c("condition","Control","Treated"),alpha = 0.05)
res <- res[order(res$padj),]

res_df <- as.data.frame(res)
gene_names <- mapIds(org.Hs.eg.db, keys=row.names(res_df), column="SYMBOL", keytype="ENSEMBL", multiVals="first")
# Add the gene_names to res_df
res_df <- res_df %>%
  dplyr::mutate(genes = gene_names[rownames(res_df)])
res_df <- res_df[order(res_df$padj),]

GO_database <- org.Hs.eg.db
KEGG_database <- 'hsa'
res_df_filter <- subset(res_df, 
                 !is.na(pvalue) & 
                 !is.na(padj))
gene <- bitr(res_df_filter$genes,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
gene <- dplyr::distinct(gene,ENTREZID,.keep_all=TRUE)

data2 <- res_df_filter %>%
  inner_join(gene,by=c("genes"="SYMBOL"))
data_sort <- data2 %>%
    arrange(desc(log2FoldChange))
data_sort_df <- as.data.frame(data_sort)
write.csv(data_sort_df, file = paste0(wkdir, "/results/pipe_plots/differential_gene_2exclude.csv"))
kegg_gene_list <- data_sort$log2FoldChange
names(kegg_gene_list) <- data_sort$ENTREZID

################################ Use dotplot to display #################################
library(DOSE)
generate_dotplot <- function(ont) {
  gse <- gseGO(geneList=kegg_gene_list, 
               ont = ont, 
               keyType = "ENTREZID", 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.05, 
               verbose = TRUE, 
               OrgDb = org.Hs.eg.db, 
               pAdjustMethod = "none")
  gse_df <- as.data.frame(gse)
  dir.create(paste0(wkdir, "/results/pipe_plots/GO/"))
  write.csv(gse_df, file = paste0(wkdir, "/results/pipe_plots/GO/", "GSE_GO_", ont, "_2exclude.csv"))
  # p <- dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign) + ggtitle(ont)
  p <- dotplot(gse, showCategory=10, split=".sign", color="pvalue") + facet_grid(.~.sign) + ggtitle(ont)
  ggsave(paste0(wkdir, "/results/pipe_plots/GO/", "GSE_GO_", ont, "_2exclude.png"), plot = p, width = 9, height = 12, dpi = 300)
}

# Generate the plots for 'BP', 'CC', and 'MF'
generate_dotplot('BP')
generate_dotplot('CC')
generate_dotplot('MF')

# KEGG
res <- gseKEGG(
  kegg_gene_list,    # 根据logFC排序的基因集
  organism = "hsa",    # 人的拉丁名缩写
  minGSSize = 3,
  pvalueCutoff = 0.05,
  pAdjustMethod = "none"
)
res_df <- as.data.frame(res)
dir.create(paste0(wkdir, "/results/pipe_plots/KEGG/"))
write.csv(res_df, file = paste0(wkdir, "/results/pipe_plots/KEGG/", "GSE_KEGG_2exclude.csv"))

k <- dotplot(
  res,
  showCategory=10,
  color="pvalue",
  split=".sign") + facet_grid(.~.sign) + ggtitle("KEGG")
ggsave(paste0(wkdir, "/results/pipe_plots/KEGG/", "GSE_KEGG_2exclude.png"), plot = k, width = 9, height = 12, dpi = 300)

