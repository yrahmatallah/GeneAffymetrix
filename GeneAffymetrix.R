rm(list=ls())
# Install package BiocManager which enables installing Bioconductor packages if it has not been already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Install package GEOquery to download data from Gene Expression Omnibus (GEO) repository
if (!requireNamespace("GEOquery", quietly = TRUE)) BiocManager::install("GEOquery")

library(GEOquery)
# GEO series GSE15641 has 32 clear cell renal cell carcinoma kidney tumor samples 
# and 23 normal kidney samples. Additional samples are available but will not be used. 

# Download preprocessed data which contains needed sample labels
gse <- getGEO(GEO="GSE15641", destdir=getwd())

# Build phenotypic data table from extracted sample labels
tab <- pData(gse[[1]])[, c("title","geo_accession")]
head(tab)
labels <- tab[,"title"]
names(labels) <- tab[,"geo_accession"]

# Extract the part of labels before the keyword " kidney "
ss <- strsplit(x=labels, split=" kidney ", fixed=TRUE)
class.labels <- unlist(lapply(ss, FUN=function(x){return(x[1])}))

# Remove spaces in labels
class.labels <- gsub(class.labels, pattern=" ", replacement="")

# Create a dataframe of phenotypic data, including group labels and file names
df <- data.frame("Name"=names(class.labels), 
"FileName"=paste(names(class.labels), ".CEL", sep=""), "Target"=class.labels)
head(df)

# Download raw files as a compressed folder under a created folder carrying the name of the GEO series
filePaths <- getGEOSuppFiles("GSE15641")

# Set the working directory to the newly created folder where downloaded data is saved
setwd(paste(getwd(), "/GSE15641", sep=""))

# Decompress the download tar folder GSE15641_RAW.tar
untar("GSE15641_RAW.tar")

# Processing steps for raw files
# Delete the compressed files since no longer needed
file.remove(c("GSE15641_RAW.tar", "GSE15641_mas5_data.txt.gz"))
# find all raw (.CEL) files in the current working directory
fns <- list.files(path=getwd(), pattern=".CEL", full.names=TRUE)
# function gunzip from package R.utils will be needed to unzip files
if (!requireNamespace("R.utils", quietly = TRUE)) BiocManager::install("R.utils")
library(R.utils)
# unzip CEL files without keeping the original zipped files
sapply(fns, gunzip, remove=TRUE)

# Read raw files into one ExpressionSet object using package affy
# install package affy if it has not been already installed
if (!requireNamespace("affy", quietly = TRUE)) BiocManager::install("affy")
library(affy)

# Read the phenotypic data table into an annotated data table object
adf <- AnnotatedDataFrame(df)

# Get group labels
pd <- pData(adf)$Target
# Get the indices of only two groups: normal and clear cell renal cell carcinoma
ind <- which((pd == "normal") | (pd == "clearcellRCC"))
adf <- adf[ind]

# Get CEL file names that match file names in the annotated dataframe
cels <- list.files(pattern=".CEL")
cels <- intersect(cels, pData(adf)$FileName)

# Read CEL files and phenotypic data into one affybatch object
# Package will be needed. You can install the package or let it be installed automatically by command ReadAffy
library(hgu133acdf)
mybatch <- ReadAffy(filenames=cels, phenoData=adf)
mybatch

# Perform quality assessment using Normalized Unscaled Standard Error (NUSE)
# Install package affyPLM if it has not been already installed
if (!requireNamespace("affyPLM", quietly = TRUE)) BiocManager::install("affyPLM")
library(affyPLM)
dataPLM <- fitPLM(mybatch)

# Generate boxplots of calculated NUSE for samples
par(mfcol=c(1,1), mar=c(8,4,2,2)) # create 2-by-1 plotting panels
aa <- boxplot(dataPLM, main="NUSE", ylim=c(0.95, 1.15), outline=FALSE,
col="lightblue", las=3, whisklty=0, staplelty=0)
abline(h=1, lty=3)

# Identify poor quality samples as those deviating from the majority of the dataset
# We deem samples with median NUSE >1.02 as poor quality samples
badSamples.index <- which(aa$stats[3,] > 1.02)
# bad samples consist of 5 normal and 5 clearcellRCC samples
mybatch$Target[badSamples.index]

# Discard samples deemed of poor quality
mybatch.filtered <- mybatch[, -badSamples.index]

# Perform preprocessing steps: Background correction, normalization, and summarization
# using the Robust Multi-Array (RMA) method. The output is an ExpressionSet object
eSet <- rma(mybatch.filtered, background=TRUE, normalize=TRUE)
eSet

# Use package genefilter to summarize intensities from probes mapping to the same gene
# Install package genefilter if it has not been already installed
if (!requireNamespace("genefilter", quietly = TRUE)) BiocManager::install("genefilter")
# Install the annotation package for the Affynetrix chip if it has not been already installed
if (!requireNamespace("hgu133a.db", quietly = TRUE)) BiocManager::install("hgu133a.db")
# Install package org.Hs.eg.db if it has not been already installed (needed to map between identifiers)
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db")
library(genefilter)
library(hgu133a.db)
library(org.Hs.eg.db)

# Summarize probes pointing to the same gene by selecting the probe with largest Inter-Quartile Range (IQR)
filteredSet <- nsFilter(eSet, require.entrez=TRUE, remove.dupEntrez=TRUE, var.func=IQR, var.filter=TRUE, var.cutof=0.1)
fset <- filteredSet$eset
mat <- exprs(fset)
mat[1:4,1:4]

# Get gene SYMBOL annotations for the probe identifiers
# Install package AnnotationDbi if it has not been already installed
if (!requireNamespace("AnnotationDbi", quietly = TRUE)) BiocManager::install("AnnotationDbi")
library(AnnotationDbi)
annot <- select(hgu133a.db, keys=rownames(mat), columns=c("SYMBOL","GENENAME"), keytype="PROBEID")
rownames(mat) <- annot[,"SYMBOL"]
class.labels <- eSet$Target
mat[1:4,1:4]

# Generate PCA scatter plot
pca <- prcomp(t(mat), scale=T)
a <- pca$sdev
pc1 <- a[1]^2/sum(a^2)
pc2 <- a[2]^2/sum(a^2)
# Set point symbols and colors for the normal and clearcellRCC groups
colors <- ifelse(class.labels=="normal", "green", "red")
symbols <- ifelse(class.labels=="normal", 16, 17)
# Generate PCA plot
plot(pca$x[,1:2], pch=symbols, cex=1.5, col=colors, main="2D PCA plot",
xlab=paste("PC1 (", round(100*pc1, 2), "%)", sep=""), 
ylab=paste("PC2 (", round(100*pc2, 2), "%)", sep=""))
legend("topright", legend=c("Normal","clearcellRCC"), col=c("green","red"), pch=c(16,17), pt.cex=1.5)

# Perform differential expression (DE) analysis using package limma
# Install package limma if it has not been installed yet
if (!requireNamespace("limma", quietly = TRUE)) BiocManager::install("limma")
library(limma)
group.labels <- ifelse(class.labels=="normal", 1, 2)
design.matrix <- model.matrix(~group.labels)
fit <- lmFit(mat, design=design.matrix)
fit <- eBayes(fit)
limma.results <- topTable(fit, n=nrow(fit))
head(limma.results)

# Get indices of significant genes at the specified threshold
# In the limma results Table, column logFC is the log2 Fold-Change. Positive value
# indicates up-regulation in clearcellRCC compared to normal
# In the limma results Table, column adj.P.Val is the adjusted p-value using 
# Benjamini-Hochberg method which ensures control over False Discovery Rate (FDR)
logFC.thr <- 3
adj.pval.thr <- 1e-10
DE.ind <- which((limma.results$adj.P.Val < adj.pval.thr) & (abs(limma.results$logFC) > logFC.thr))

# Generate a volcano plot to visualize significant genes
colors <- rep("gray", nrow(limma.results))
# Highlight significant genes with red color
colors[DE.ind] <- "red"
plot(limma.results$logFC, -log10(limma.results$adj.P.Val), type="p", pch=16, col=colors,
cex=0.8, xlab="log2FC", ylab="-log10(FDR)", main="Normal vs clearcellRCC")
abline(v=c(-logFC.thr, logFC.thr), lty=3)
abline(h=-log10(adj.pval.thr), lty=3)

# Generate heatmap for the significant DE genes using package pheatmap
if(!requireNamespace("pheatmap", quietly=TRUE)) BiocManager::install("pheatmap")
library(pheatmap)
# Extract gene expressions for only significant genes
mydata <- mat[rownames(limma.results)[DE.ind],]
# Remove the exyension '.CEL' at the end of column names
colnames(mydata) <- gsub(x=colnames(mydata), pattern=".CEL", replacement="")

 Generate column annotations for the heatmap
annotation_col <- data.frame(Group=factor(class.labels))
rownames(annotation_col) = colnames(mydata)
# Generate row annotations for the heatmap
up_down <- limma.results[rownames(mydata),"logFC"]
up_down_col <- rep("Up_in_tumor", nrow(mydata))
up_down_col[up_down < 0] <- "Down_in_tumor"
annotation_row <- data.frame(Up_Down=factor(up_down_col))
rownames(annotation_row) <- rownames(mydata)
# Set annotation colors
ann_colors = list(Group=c(normal="gold", clearcellRCC="blue"), 
Up_Down=c(Down_in_tumor="darkgreen", Up_in_tumor="darkred"))

pheatmap(mydata, scale="row", annotation_col=annotation_col, annotation_row=annotation_row,
annotation_colors=ann_colors, show_rownames=T, show_colnames=T, fontsize_row=7, fontsize_col=7)

sessionInfo()
