##ISTP Biostatistics Course Final 
##Sophie Son, Mahmoud Yousef, Joseph Cozzi 
library(dplyr)
library(ggplot2)
library(NMF)
library(SummarizedExperiment)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
library("DESeq2")
library("RColorBrewer") # Load a package giving more colors
library("pheatmap") # load a package for making heatmaps

#set working directory and import dataset
setwd("/Users/seuleson/Desktop/MSTPsummerstatistics-master/Biostats2020-master/")
raw.data <- as.matrix(read.csv("/Users/seuleson/Desktop/MSTPsummerstatistics-master/Biostats2020-master/GSE130683_LPS_IL4_supplemental_data.csv"))
rownames(raw.data) <- raw.data[,1]
raw.data <- raw.data[,2:18]
class(raw.data) <- "numeric"
log2.counts <- log2(raw.data +1)

##plot counts for outliers
boxplot(raw.data, main="Raw Counts")  #plot raw counts
hist(raw.data[,2])

boxplot(log2.counts, main="log2 transformed counts")
hist(log2.counts[,2])

##coordinating data into usable form
#compiling column data 
data.coldata <- colnames(read.csv("/Users/seuleson/Desktop/MSTPsummerstatistics-master/Biostats2020-master/GSE130683_LPS_IL4_supplemental_data.csv", row.names = 1))
data.coldata <- rbind(data.coldata, c("IL4", "IL4", "IL4", "IL4", "IL4", "IL4", "IL4", "IL4", "LPS", "LPS", "LPS", "LPS", "LPS", "PBS", "PBS", "PBS", "PBS"), c("anti-inflammatory", "anti-inflammatory", "anti-inflammatory", "anti-inflammatory", "anti-inflammatory", "anti-inflammatory", "anti-inflammatory", "anti-inflammatory", "pro-inflammatory", "pro-inflammatory", "pro-inflammatory","pro-inflammatory", "pro-inflammatory", "negative control", "negative control", "negative control", "negative control"))
rownames(data.coldata) <- c("sample", "treatment", "treatment type")
data.coldata<-t(data.coldata)

##construct a master DEseqdataset (dds)
data.dds <- DESeqDataSetFromMatrix(countData = raw.data, colData = data.coldata, design = ~ treatment)
#filter out extremely low read samples (only keep if rows has at least 10 reads)
keep.cond <- rowSums(counts(data.dds)) >= 10
data.dds <- data.dds[keep.cond,]
#set factor levels of comparison conditions
data.dds$treatment <- factor(data.dds$treatment, levels=c("IL4", "LPS", "PBS"))

#regularized log transformation of DESeq object 
rlog.dds <- rlogTransformation(data.dds)

#plotting effects of rlog compared to log2 transformed counts 
#observe smaller dispersion effects
par(mfrow =c(1,2))
plot(log2( 1 + counts(data.dds)[,1:2]),
     pch=16, cex=0.3, main = "log2")
plot(assay(rlog.dds)[ ,1:2], # The assay function returns the count values of rld in a matrix
     pch=16, cex=0.3, main = "rlog")

#plotting heatmap of counts 
distance.dds <- dist(t(assay(rlog.dds))) # Calculate distances using transformed (and normalized) counts
mat.dist <- as.matrix(distance.dds) # convert to matrix
rownames(mat.dist) <- colnames(mat.dist) <- with(colData(data.dds), paste(treatment)) # set rownames in the matrix
colnames(mat.dist) = NULL # remove column names
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) # set colors
pheatmap(mat.dist,
         clustering_distance_rows=distance.dds,
         clustering_distance_cols=distance.dds,
         col=colors)

#normalization of data: observe if estimate size factors are approximately 1
data.dds <- estimateSizeFactors(data.dds)
sizeFactors(data.dds)

##differential expression analysis (DEA)
deseq.dds <- DESeq(data.dds)
dea.results.il4 <- results(deseq.dds, contrast = c("treatment", "IL4", "PBS"))
dea.results.lps <- results(deseq.dds, contrast = c("treatment", "LPS", "PBS"))
dea.results.il4vlps <- results(deseq.dds, contrast = c("treatment", "IL4", "LPS"))
#results represents log2 fold changes, p values, and p adjusted values 

#reorder results based on p values (small -> large)
dea.results.il4.pval.order <- dea.results[order(dea.results.il4$pvalue),]
dea.results.lps.pval.order <- dea.results[order(dea.results.lps$pvalue),]

#filter significant genes based on p adjusted <0.01
signif.result.il4 <- subset(dea.results.il4, padj < 0.1)
signif.result.il4 
signif.result.lps <- subset(dea.results.lps, padj < 0.1)
signif.result.lps

##plotting results
#plot log2 fold changes of a given variable over mean of normalized counts for all samples in dataset
plotMA(dea.results, ylim=c(-2,2))

#heatmap of most significant differentially expressed genes
pheatmap()
