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
library("biomaRt")

#set working directory and import dataset
setwd("/Users/seuleson/Desktop/MSTPsummerstatistics-master/Biostats2020-master/")
raw.data <- as.matrix(read.csv("/Users/seuleson/Desktop/MSTPsummerstatistics-master/Biostats2020-master/GSE130683_LPS_IL4_supplemental_data.csv"))
metadata <- as.matrix(read.csv("/Users/seuleson/Desktop/MSTPsummerstatistics-master/Biostats2020-master/Supplement-metadata.csv"))[,1:10]

rownames(raw.data) <- raw.data[,1]
raw.data <- raw.data[,2:18]
class(raw.data) <- "numeric"
log2.counts <- log2(raw.data +1)

##plot counts for outliers
par(mfrow=c(2,1))
sample.labels <- gsub("_CRC38m", "", colnames(raw.data))

boxplot(raw.data, main="Raw Counts", xaxt="n") 
axis(1,tick= F, labels=FALSE)
text(x = 1:ncol(raw.data), y = par("usr")[3]-10, labels = sample.labels, xpd = TRUE, srt = 35, adj = 0.965, cex = 0.95)
#plot raw counts
hist(raw.data[,2], ylab = colnames(raw.data[,2]))

boxplot(log2.counts, main="log2 transformed counts", las=2)
hist(log2.counts[,2])
dev.off()

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

#regularized log transformation of DESeq object: moderates the variance across the mean to improve clustering 
  #log transform the everage across samples of each genes normalized count --> prevents effects of size on variance (RNA-seq data is not homoscedastic)
    #after shrinkage/transformation --> results are homoscedastic --> can use an ordinary t-test 
      # t-test assumes genes are independent (probably not true)
rlog.dds <- rlogTransformation(data.dds, blind= TRUE)

#plotting effects of rlog compared to log2 transformed counts 
#observe smaller dispersion effects
par(mfrow =c(1,2))
plot(log2( 1 + counts(data.dds)[,1:2]),
     pch=16, cex=0.3, main = "log2")
plot(assay(rlog.dds)[ ,1:2], # The assay function returns the count values of rld in a matrix
     pch=16, cex=0.3, main = "rlog")
dev.off()

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

#principal component analysis of significant differentially expressed genes
plotPCA(rlog.dds, intgroup = "treatment") 
pc.dds <- plotPCA(rlog.dds, intgroup = "treatment", returnData=T) 

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
dea.results.il4.pval.order <- dea.results.il4[order(dea.results.il4$pvalue),]
dea.results.lps.pval.order <- dea.results.lps[order(dea.results.lps$pvalue),]
dea.results.il4vlps.pval.order <- dea.results.il4vlps[order(dea.results.il4vlps$pvalue),]

#filter significant genes based on p adjusted <0.01
signif.result.il4 <- subset(dea.results.il4, padj < 0.05)
signif.result.il4 <- subset(signif.result.il4, log2FoldChange > 1)
signif.result.il4 
signif.result.lps <- subset(dea.results.lps, padj < 0.05)
signif.result.lps <- subset(signif.result.lps, log2FoldChange > 1)
signif.result.lps
signif.result.il4vlps <- subset(dea.results.il4vlps, padj <0.05)
signif.result.il4vlps <- subset(signif.results.il4vlps, log2FoldChange > 1)
signif.result.il4vlps

##plotting results
#plot log2 fold changes of a given variable over mean of normalized counts for all samples in dataset
plotMA(dea.results.il4, ylim=c(-10,10), main= "IL4 vs PBS")
plotMA(dea.results.lps, ylim=c(-10,10), main="LPS vs PBS")
plotMA(dea.results.il4vlps, ylim=c(-10,10), main ="IL4 vs LPS")

#Find the 10 most upregulated
head(signif.result.il4[order(signif.result.il4$log2FoldChange),], 10)
head(signif.result.lps[order(signif.result.lps$log2FoldChange),], 10)
head(signif.result.il4vlps[order(signif.result.il4vlps$log2FoldChange),], 10)
#Find the 10 most downregulated 
head(signif.result.il4[order(signif.result.il4$log2FoldChange),decreasing = TRUE], 10)
head(signif.result.lps[order(signif.result.lps$log2FoldChange),decreasing = TRUE], 10)
head(signif.result.il4vlps[order(signif.result.il4vlps$log2FoldChange),decreasing = TRUE], 10)

#Filter genes based on highest level of variance and plot with heatmap
library(genefilter)
library(gplots)
Var.result.dds <- head(order(rowVars(assay(rlog.dds)), decreasing =TRUE), 35)
pdf("heatmap.most.variance.pdf", width = 30, height = 40)
heatmap.2( assay(rlog.dds)[Var.result.dds,], scale="row",
           trace="none", dendrogram="column", margins = c(12,9),
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))
dev.off()


