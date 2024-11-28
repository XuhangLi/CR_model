install.packages("Seurat")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")

library(stringr)
library(ggplot2)
library(ggpubr)
library(matrixStats)
library(Seurat)
library(SeuratDisk)
library(pheatmap)


# ORIGINAL EDA

# # Load the .h5ad file into a Seurat object
# Convert("./input/source_data/processed_data/K562_gwps_normalized_bulk_01.h5ad", dest = "h5seurat", overwrite = TRUE)
# seurat_object <- LoadH5Seurat("./input/source_data/processed_data/K562_gwps_normalized_bulk_01.h5seurat")
# 
# # Assuming your Seurat object is named seurat_object
# # get the normalzied expression matrix
# exp_mat <- seurat_object@assays$RNA@counts # same as scale.data so it is scaled data
# 
# # get the metadata for each perturbations 
# cell_annotations <- seurat_object@meta.data
# head(cell_annotations)  # View the first few rows of the cell annotations
# 
# # get metadata for each genes for future filtering 
# gene_annotations <- seurat_object@assays$RNA@meta.features
# head(gene_annotations)  # View the first few rows of the cell annotations
# 
# # subset metabolic genes and metabolic perturbations and visualize 
# metGenes = read.csv('input/human_model/genes.tsv', sep = '\t', header = T)
# # get the met perturbations
# cell_annotations$perturbed_gene = rownames(cell_annotations)
# cell_annotations$perturbed_gene = unlist(lapply(strsplit(cell_annotations$perturbed_gene, '_'),function(x){x[[4]]}))
# # subset matrix 
# exp_mat_met = as.data.frame(exp_mat[rownames(exp_mat) %in% metGenes$genes, rownames(cell_annotations)[cell_annotations$perturbed_gene %in% metGenes$genes]])
# exp_mat_met[exp_mat_met == Inf] = max(exp_mat_met[exp_mat_met!=Inf])
# 
# # focus on responsive ones
# exp_mat_met_resp = exp_mat_met[, colnames(exp_mat_met) %in% rownames(cell_annotations_met)[cell_annotations_met$anderson_darling_counts > 10]]
# 
# # visualization
# color_palette <- rev(colorRampPalette(c("red", "grey100", "blue"))(100))
# upperl = 2
# seq1 = seq(-upperl,upperl,length.out = 100)
# pheatmap(exp_mat_met_resp, color = color_palette,breaks = seq1
#          ,show_rownames = F,show_colnames = F
# )
# 
# # PDH
# hist(as.numeric(exp_mat_met['ENSG00000131828',]))
# hist(as.numeric(exp_mat_met[,'6175_PDHA1_P1P2_ENSG00000131828']))
# 
# exp_mat_met['ENSG00000131828','6175_PDHA1_P1P2_ENSG00000131828']
# 
# 
# # CS 
# hist(as.numeric(exp_mat_met['ENSG00000062485',]))
# hist(as.numeric(exp_mat_met[,'1919_CS_P1P2_ENSG00000062485']))
# 
# exp_mat_met['ENSG00000062485','1919_CS_P1P2_ENSG00000062485']
# 
# # HK1 
# hist(as.numeric(exp_mat_met['ENSG00000156515',]))
# hist(as.numeric(exp_mat_met[,'3803_HK1_P2_ENSG00000156515']))
# 
# exp_mat_met['ENSG00000156515','3803_HK1_P2_ENSG00000156515']
# 
# # cell ann met
# cell_annotations_met = cell_annotations[cell_annotations$perturbed_gene %in% metGenes$genes,]
# 
# 
# # SDHc
# hist(as.numeric(exp_mat_met['ENSG00000143252',]))
# hist(as.numeric(exp_mat_met[,'7727_SDHC_P1P2_ENSG00000143252']))
# 
# exp_mat_met['ENSG00000143252','7727_SDHC_P1P2_ENSG00000143252']
# 
# 
# # MED19
# hist(as.numeric(exp_mat['ENSG00000156603',]))
# hist(as.numeric(exp_mat[,'4941_MED19_P1P2_ENSG00000156603']))
# 
# exp_mat['ENSG00000156603','4941_MED19_P1P2_ENSG00000156603']
# 
# # next give it a go with FBA and visualization 
# # also will consider merge the essential gene matrix to increase coverage and signals


# generate some figures 
# number of up and down DEGs for metabolic and all DEGs and across thresholds 

# load DEG result
adDE_res = read.csv('input/source_data/processed_data/anderson-darling_p-values_BH-corrected.csv', row.names = 1)
# although not noted in the figshare page where the data is downloaded, the column number (# perturbations) is equal to 
# that of the K562_gwps (K562_gwps_normalized_bulk_01.h5seurat), indicating this is the result of the genome-wide perturb-seq dataset in K562 cells 

# we will use the DE result and skip the thresholding on pseudobulk data (because it is z-score anyways - not fold change)
# (we can use raw UMI count (like Count-per-Ten-Thousound)) but that is more noisy and has batch problems)
# Load the .h5ad file into a Seurat object for meta information and sign of changes 
seurat_object <- LoadH5Seurat("./input/source_data/processed_data/K562_gwps_normalized_bulk_01.h5seurat")

# Assuming your Seurat object is named seurat_object
# get the normalzied expression matrix
exp_mat <- seurat_object@assays$RNA@counts # same as scale.data so it is scaled data
colnames(exp_mat) = paste('X',colnames(exp_mat),sep = '')

# get the metadata for each perturbations 
cell_annotations <- seurat_object@meta.data
head(cell_annotations)  # View the first few rows of the cell annotations

# get metadata for each genes for future filtering 
gene_annotations <- seurat_object@assays$RNA@meta.features
head(gene_annotations)  # View the first few rows of the cell annotations

# subset metabolic genes and metabolic perturbations and visualize 
metGenes = read.csv('input/human_model/genes.tsv', sep = '\t', header = T)


# get the met perturbations
perturbed_gene = colnames(adDE_res)
perturbed_gene = unlist(lapply(strsplit(perturbed_gene, '_'),function(x){x[[4]]}))

# subset matrix for a metabolic GRN
padj_mat_met = adDE_res[, perturbed_gene %in% metGenes$genes]

# get the overlap (with normalized expression data)
cell_annotations <- seurat_object@meta.data
cell_annotations$perturbed_gene = perturbed_gene[match(rownames(cell_annotations), str_remove(colnames(adDE_res),'^X'))]
cell_annotations$perturbation_id = paste('X',rownames(cell_annotations),sep = '')
padj_mat_met = padj_mat_met[, colnames(padj_mat_met) %in% colnames(exp_mat)]

# the z-score matrix 
exp_mat_met = as.matrix(exp_mat[rownames(padj_mat_met), colnames(padj_mat_met)])

# mask change that is the RNAi targeted gene 
for (i in 1:ncol(padj_mat_met)){
  if (any(perturbed_gene[colnames(padj_mat_met)[i] == colnames(adDE_res)] == rownames(padj_mat_met))){
    padj_mat_met[perturbed_gene[colnames(padj_mat_met)[i] == colnames(adDE_res)],i] = 1 # mask to insignificant (1)
  }
}

# plot the number of DEGs up and down given threshold 
cutoffs = c(10^-(10:2), 0.05, seq(0.1,0.9,0.1))
N_up = c()
N_down = c()
N_up_met = c()
N_down_met = c()
for (i in 1:length(cutoffs)){
  N_up = c(N_up, sum(padj_mat_met < cutoffs[i] & exp_mat_met > 0))
  N_down = c(N_down, sum(padj_mat_met < cutoffs[i] & exp_mat_met < 0))
  N_up_met = c(N_up_met, sum(padj_mat_met[rownames(padj_mat_met) %in% metGenes$genes,] < cutoffs[i] & exp_mat_met[rownames(padj_mat_met) %in% metGenes$genes,] > 0))
  N_down_met = c(N_down_met, sum(padj_mat_met[rownames(padj_mat_met) %in% metGenes$genes,] < cutoffs[i] & exp_mat_met[rownames(padj_mat_met) %in% metGenes$genes,] < 0))
  
}
plot(log10(cutoffs), N_up/(N_up+N_down), type = 'l', col = 'red', 
     xlab = 'FDR cutoff', ylab = 'Number of DEGs', main = 'Number of DEGs up and down across FDR cutoffs',
     ylim = c(0.3, 0.7))
lines(log10(cutoffs), N_down/(N_up+N_down), col = 'blue')
lines(log10(cutoffs), N_up_met/(N_up_met+N_down_met), col = 'red', lty = 2)
lines(log10(cutoffs), N_down_met/(N_up_met+N_down_met), col = 'blue', lty = 2)

pdf('figures/DEG_direction.pdf')
hist(exp_mat_met[padj_mat_met<0.01],breaks = 100, xlab = 'Normalized gene expression')
# plot met genes
exp_mat_met2 = exp_mat_met[rownames(padj_mat_met) %in% metGenes$genes,]
padj_mat_met2 = padj_mat_met[rownames(padj_mat_met) %in% metGenes$genes,]
hist(exp_mat_met2[padj_mat_met2<0.01],breaks = 100, add = T, col = 'red', alpha = 0.5)
dev.off()

# also plot the percentage of increase in each RNA - more difficult to look at (because the pattern is not as striking as worm data)
# allRNAi = unique(colnames(padj_mat_met)[colSums(padj_mat_met<0.01)>9])
# UPperc = c()
# for (i in 1:length(allRNAi)){
#   UPperc = c(UPperc, sum(padj_mat_met[,allRNAi[i]]<0.01 & exp_mat_met[,allRNAi[i]] > 0) / sum(padj_mat_met[,allRNAi[i]]<0.01))
# }
# dev.off()
# pdf('figures/0_DE_summary/FC_distribution_UP_DE_percentage.pdf',width = 5,height = 5)
# hist(UPperc,xlab = 'Fraction of up-regulated DEGs',breaks = 30)
# dev.off()

# plot indegree distribution 
# first define a very rough GRN without too much filtering (remove zeros)
DEG_calls = padj_mat_met<0.01
library(tidyverse)
library(reshape2)

DEG_calls =  melt(DEG_calls)
DEG_calls = DEG_calls[DEG_calls$value,]

# calculate in-degree
inDegree = table(DEG_calls$Var1)
hist(log2(inDegree+1),breaks = 10,right = F,xlab = 'log2(in degree+1)')
# show the scale free distribution of in-degree
ct = table(inDegree)
# show the accumulated DE curve 
freq = sort(inDegree,decreasing = T)
freq_acc = c()
for (i in 1:length(freq)){
  freq_acc = c(freq_acc, sum(freq[1:i]))
}
freq_acc = freq_acc / sum(freq)

# add metabolic genes
metGenes = read.csv('input/human_model/genes.tsv', sep = '\t', header = T)
freq_icel = sort(freq[names(freq) %in% metGenes$genes],decreasing = T)
freq_acc_icel = c()
for (i in 1:length(freq_icel)){
  freq_acc_icel = c(freq_acc_icel, sum(freq_icel[1:i]))
}
freq_acc_icel = freq_acc_icel / sum(freq_icel)

pdf('figures/human_GRN_in_degree_accumulative.pdf',width = 5,height = 5)
plot( (1:length(freq_acc))/length(freq_acc),freq_acc, lty = 'solid', type = 'l',
      xlab = 'Fraction of genes in mGRN', ylab = 'Fraction of total interactions')
abline(v = 0.2, lty = 'dotted', col = 'red')
abline(h = 0.6, b = 1, lty = 'dotted', col = 'red')
dev.off()

