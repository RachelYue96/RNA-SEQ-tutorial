#####################################################
## Date: 2023/12/8
## Copyright (c) Salus BioMed
## Author: Rahcel
## Description: This file is for summarizing the output files from
##       featureCounts and form a matrix of feature for 
##       all samples.
#####################################################

library('biomaRt')

wkdir <- "/share/result/sequencer/salus/video_example/RNA-seq"

#################################### 1. Compute counts for genes and transcripts #################################
files_allg <- list.files(path=paste0(wkdir, "/results/hisat2_stringtie"), pattern="*_allgene.counts.txt$", full.names=TRUE, recursive=TRUE)
files_trans <- list.files(path=paste0(wkdir, "/results/hisat2_stringtie"), pattern="*_transcript.counts.txt", full.names=TRUE, recursive=TRUE)

names(files_allg) <- gsub("_allgene.counts.txt$", "", basename(files_allg))
exprs <- list()
exprs_trans <- list()

for (i in 1:length(files_allg)) {
  exprs[[i]] <- read.table(paste0(files_allg[i]), col.names = c("gene_id", names(files_allg)[i]), check.names = FALSE)
}

for (i in 1:length(files_trans)) {
  exprs_trans[[i]] <- read.table(paste0(files_trans[i]), col.names = c("gene_id", names(files_allg)[i]), check.names = FALSE)
}

merged_exprs <- Reduce(function(x, y) merge(x, y, by = "gene_id"), exprs)
merged_exprs <- merged_exprs[-nrow(merged_exprs), ]
write.csv(merged_exprs, file = paste0(wkdir, "/results/RNA_counts/allgene_count_matrix.txt"),row.names = FALSE)
merged_exprs_trans <- Reduce(function(x, y) merge(x, y, by = "gene_id"), exprs_trans)
merged_exprs_trans <- merged_exprs_trans[-nrow(merged_exprs_trans), ]
write.csv(merged_exprs_trans, file = paste0(wkdir, "/results/RNA_counts/transcript_count_matrix.txt"),row.names = FALSE)

###################################### 2. Calculate FPKM and TPM from raw counts ##################################
files_L_allg <- list.files(path=paste0(wkdir, "/results/hisat2_stringtie"), pattern="*_allgene.lengths.txt$", full.names=TRUE, recursive=TRUE)
files_L_trans = list.files(path=paste0(wkdir, "/results/hisat2_stringtie"), pattern="*_transcript.lengths.txt$", full.names=TRUE, recursive=TRUE)

## For gene-level length extraction
legnths <- read.table(paste0(files_L_allg[1]), header = TRUE)
legnths <- aggregate(legnths[, 2], list(legnths[, 1]), mean)
colnames(legnths) <- c("gene_id", "length")
new_merged_exprs <- merge(merged_exprs, legnths, by = "gene_id")
# write.csv(new_merged_exprs, file = paste0(wkdir, "/results/RNA_counts/temp.txt"), row.names = FALSE)

## For transcript-level length extraction
length_trans <- read.table(paste0(files_L_trans[1]), header = TRUE)
new_merged_exprs_trans <- merge(merged_exprs_trans, length_trans, by = "gene_id")

## mapped reads in million --> is the same for gene and transcript level
files_summary <- list.files(path=paste0(wkdir, "/results/hisat2_stringtie"), pattern="*.tran_summary.txt$", full.names=TRUE, recursive=TRUE)
total_mapped_reads <- list()

for (i in 1:length(files_summary)) {
  total_mapped_reads[[i]] <- as.numeric(strsplit(readLines(files_summary[i], n = 1)[1], " ")[[1]][1]) * 2 / 1e6
}

## calculation for fpkm
fpkm <- function(counts, lengths, total_mapped_reads) {
  counts <- as.numeric(counts)
  lengths <- as.numeric(lengths)
  rpm <- counts / total_mapped_reads
  (rpm) / (lengths / 1e3)
}

## in gene-level
fpkm_matrix <- data.frame(gene_id = new_merged_exprs$gene_id)
last_index <- 0
for (i in 1:(ncol(new_merged_exprs)-2)) {
  last_index <- i
  fpkm_values <- fpkm(new_merged_exprs[,i+1], new_merged_exprs$length, total_mapped_reads[[i]])
  fpkm_matrix <- cbind(fpkm_matrix, fpkm_values)
  colnames(fpkm_matrix)[i] <- colnames(new_merged_exprs)[i]
}
colnames(fpkm_matrix)[last_index+1] <- colnames(new_merged_exprs)[last_index+1]
write.csv(fpkm_matrix, file = paste0(wkdir, "/results/RNA_counts/allgene_fpkm_matrix.txt"), row.names = FALSE)

## in transcript level
fpkm_matrix_trans <- data.frame(gene_id = new_merged_exprs_trans$gene_id)
last_index <- 0
for (i in 1:(ncol(new_merged_exprs_trans)-2)) {
  last_index <- i
  fpkm_values_trans <- fpkm(new_merged_exprs_trans[,i+1], new_merged_exprs_trans$length, total_mapped_reads[[i]])
  fpkm_matrix_trans <- cbind(fpkm_matrix_trans, fpkm_values_trans)
  colnames(fpkm_matrix_trans)[i] <- colnames(new_merged_exprs_trans)[i]
}
colnames(fpkm_matrix_trans)[last_index+1] <- colnames(new_merged_exprs_trans)[last_index+1]
write.csv(fpkm_matrix_trans, file = paste0(wkdir, "/results/RNA_counts/transcript_fpkm_matrix.txt"), row.names = FALSE)

## calculation for TPM
tpm <- function(counts, lengths) {
  counts <- as.numeric(counts)
  lengths <- as.numeric(lengths)
  rpk <- counts / (lengths / 1e3)
  per_million_scaling_factor <- sum(rpk) / 1e6
  tpm_values <- rpk / per_million_scaling_factor
  return(tpm_values)
}

## in gene-level
tpm_matrix <- data.frame(gene_id = new_merged_exprs$gene_id)
last_index <- 0
for (i in 1:(ncol(new_merged_exprs)-2)) {
  last_index <- i
  tpm_values <- tpm(new_merged_exprs[,i+1], new_merged_exprs$length)
  tpm_matrix <- cbind(tpm_matrix, tpm_values)
  colnames(tpm_matrix)[i] <- colnames(new_merged_exprs)[i]
}
colnames(tpm_matrix)[last_index+1] <- colnames(new_merged_exprs)[last_index+1]
write.csv(tpm_matrix, file = paste0(wkdir, "/results/RNA_counts/allgene_tpm_matrix.txt"), row.names = FALSE)

## in transcript level
tpm_matrix_trans <- data.frame(gene_id = new_merged_exprs_trans$gene_id)
last_index <- 0
for (i in 1:(ncol(new_merged_exprs_trans)-2)) {
  last_index <- i
  tpm_values_trans <- tpm(new_merged_exprs_trans[,i+1], new_merged_exprs_trans$length)
  tpm_matrix_trans <- cbind(tpm_matrix_trans, tpm_values_trans)
  colnames(tpm_matrix_trans)[i] <- colnames(new_merged_exprs_trans)[i]
}
colnames(tpm_matrix_trans)[last_index+1] <- colnames(new_merged_exprs_trans)[last_index+1]
write.csv(tpm_matrix_trans, file = paste0(wkdir, "/results/RNA_counts/transcript_tpm_matrix.txt"), row.names = FALSE)

##################################################################################################################
