################################################################################################################################
# Got the tutorial from https://bioinformatics.uconn.edu/resources-and-events/tutorials-2/rna-seq-tutorial-with-reference-genome/
################################################################################################################################

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
library("RColorBrewer")
library("gplots")
library("genefilter")
library("ggplot2")
library("grDevices")
library("RColorBrewer")

# Set the working directory
directory <- "/share/result/sequencer/salus/video_example/RNA-seq/"
setwd(directory)

par(bg = "white")  # Ensure background is white

# Set the prefix for each output file name
outputPrefix <- "results/pipe_plots/"
dir.create(paste0(directory, outputPrefix))

countdata <- read.table("results/RNA_counts/allgene_count_matrix.txt", header=TRUE, row.names=1, check.names=FALSE, sep=",")
countdata <- as.matrix(countdata)

# Create a DataFrame that describes the condition for each sample
samples <- c('alt-1', 'alt-2', 'alt-3', 'con-1', 'con-2', 'con-3')
conditions <- c('Treated', 'Treated', 'Treated', 'Control', 'Control', 'Control')

coldata <- DataFrame(condition = factor(conditions), row.names = samples)

#guts
ddsHTSeq <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~ condition)
dds <- DESeq(ddsHTSeq)
res <- results(dds,contrast = c("condition","Control","Treated"),alpha = 0.05)

# copied from: https://benchtobioinformatics.wordpress.com/category/dexseq/
# order results by padj value (most significant to least)
res = subset(res, padj < 0.05 & abs(log2FoldChange) > 0.5)
res <- res[order(res$padj),]
# should see DataFrame of baseMean, log2Foldchange, stat, pval, padj

# save data results and normalized reads to csv
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata)[1] <- 'gene'
write.csv(resdata, file = paste0(outputPrefix, "results-with-normalized.csv"))

# send normalized counts to tab delimited file for GSEA, etc.
write.table(as.data.frame(counts(dds),normalized=T), file = paste0(outputPrefix, "normalized_counts.txt"), sep = '\t')

# produce DataFrame of results of statistical tests
mcols(res, use.names = T)
write.csv(as.data.frame(mcols(res, use.name = T)),file = paste0(outputPrefix, "test-conditions.csv"))

# replacing outlier value with estimated value as predicted by distrubution using
# "trimmed mean" approach. recommended if you have several replicates per treatment
# DESeq2 will automatically do this if you have 7 or more replicates

ddsClean <- replaceOutliersWithTrimmedMean(dds)
ddsClean <- DESeq(ddsClean)
tab <- table(initial = results(dds)$padj < 0.05,
             cleaned = results(ddsClean)$padj < 0.05)
addmargins(tab)
write.csv(as.data.frame(tab),file = paste0(outputPrefix, "replaceoutliers.csv"))
resClean <- results(ddsClean, contrast = c("condition","Control","Treated"),alpha = 0.05)
resClean = subset(res, padj<0.05)
resClean <- resClean[order(resClean$padj),]
write.csv(as.data.frame(resClean),file = paste0(outputPrefix, "replaceoutliers-results.csv"))

####################################################################################
# Exploratory data analysis of RNAseq data with DESeq2
#
# these next R scripts are for a variety of visualization, QC and other plots to
# get a sense of what the RNAseq data looks like based on DESEq2 analysis
#
# 1) MA plot
# 2) rlog stabilization and variance stabiliazation
# 3) variance stabilization plot
# 4) heatmap of clustering analysis
# 5) PCA plot
#
#
####################################################################################

# MA plot of RNAseq data for entire dataset
# http://en.wikipedia.org/wiki/MA_plot
# genes with padj < 0.1 are colored Red
Sys.setenv("DISPLAY"=":0.0")
capabilities()
Sys.getenv(c("DISPLAY"))
options(bitmapType='cairo')
png(paste0(outputPrefix, "MAplot_initial_analysis.png"), width = 800, height = 600, type = "cairo")
plotMA(dds, ylim=c(-8,8),main = "RNAseq experiment")
# dev.copy(png, paste0(outputPrefix, "MAplot_initial_analysis.png"))
dev.off()

# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.
rld <- rlogTransformation(dds, blind=T)
vsd <- varianceStabilizingTransformation(dds, blind=T)

# save normalized values
write.table(as.data.frame(assay(rld),file = paste0(outputPrefix, "rlog-transformed-counts.txt"), sep = '\t'))
write.table(as.data.frame(assay(vsd),file = paste0(outputPrefix, "vst-transformed-counts.txt"), sep = '\t'))

# plot to show effect of transformation
# axis is square root of variance over the mean for all samples
par(mai = ifelse(1:4 <= 2, par('mai'),0))
px <- counts(dds)[,1] / sizeFactors(dds)[1]
ord <- order(px)
ord <- ord[px[ord] < 150]
ord <- ord[seq(1,length(ord),length=50)]
last <- ord[length(ord)]
vstcol <- c('blue','black')
png(paste0(outputPrefix, "variance_stabilizing.png"), width = 800, height = 600, type = "cairo")
matplot(px[ord], cbind(assay(vsd)[,1], log2(px))[ord, ],type='l', lty = 1, col=vstcol, xlab = 'n', ylab = 'f(n)')
legend('bottomright',legend=c(expression('variance stabilizing transformation'), expression(log[2](n/s[1]))), fill=vstcol)
# dev.copy(png,paste0(outputPrefix, "variance_stabilizing.png"))
dev.off()

# clustering analysis
# excerpts from http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/

distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds), paste(condition, rownames(colData(dds)), sep=" : "))
#Or if you want conditions use:
#rownames(mat) <- colnames(mat) <- with(colData(dds),condition)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
# dev.copy(png, paste0(outputPrefix, "clustering.png"))
png(paste0(outputPrefix, "clustering.png"), width = 800, height = 600, type = "cairo")
heatmap.2(mat, trace = "none", col = rev(hmcol), margin = c(13,13))
dev.off()

#Principal components plot shows additional but rough clustering of samples

rv <- rowVars(assay(rld))
select <- order(rv, decreasing=T)[seq_len(min(500,length(rv)))]
pc <- prcomp(t(assay(vsd)[select,]))

# set condition
condition <- coldata
scores <- data.frame(pc$x, condition)

(pcaplot <- ggplot(scores, aes(x = PC1, y = PC2, col = (factor(condition))))
+ geom_point(size = 5)
+ ggtitle("Principal Components")
+ scale_colour_brewer(name = " ", palette = "Set1")
+ theme(
  plot.title = element_text(face = 'bold'),
  legend.position = c(.9,.2),
  legend.key = element_rect(fill = 'NA'),
  legend.text = element_text(size = 10, face = "bold"),
  axis.text.y = element_text(colour = "Black"),
  axis.text.x = element_text(colour = "Black"),
  axis.title.x = element_text(face = "bold"),
  axis.title.y = element_text(face = 'bold'),
  panel.grid.major.x = element_blank(),
  panel.grid.major.y = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.grid.minor.y = element_blank(),
  panel.background = element_rect(color = 'black',fill = NA)
))

ggsave(pcaplot,file=paste0(outputPrefix, "PCA.pdf"))

# scatter plot of rlog transformations between Sample conditions
# nice way to compare control and experimental samples
head(assay(rld))
# plot(log2(1+counts(dds,normalized=T)[,1:2]),col='black',pch=20,cex=0.3, main='Log2 transformed')
plot(assay(rld)[,1:3],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")
plot(assay(rld)[,2:4],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")
plot(assay(rld)[,6:5],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")

# heatmap of data
# 1000 top expressed genes with heatmap.2
select <- order(rowMeans(counts(ddsClean,normalized=T)),decreasing=T)[1:1000]
my_palette <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100)
png(paste0(outputPrefix, "HEATMAP.png"), width = 800, height = 600, type = "cairo")
heatmap.2(assay(vsd)[select,], col=my_palette,
          scale="row", key=T, keysize=1.5,
          symkey=T, density.info="none", trace="none",
          cexCol=1.2, labRow=F,
          main="1000 Top Expressed Genes Heatmap")
# dev.copy(png, paste0(outputPrefix, "HEATMAP.png"))
dev.off()