# Gene Affymetrix Processing and Analysis

## Objective
This educational module demonstrates the steps needed to download raw microarray files from the Gene Expression Omnibus (GEO) repository, ensure quality control by filtering out biased samples, process data (background correction, between samples normalization, and probe summarization), perform differential expression (DE) analysis, and generate visualizations.

## Introduction
Both R code and the same R code executed within a Jupyter Notebook supported by a conda environment are available. The modeule uses the GSE15641* series downloaded from GEO. The dataset has clear cell renal cell carcinoma (ccRCC) samples and normal controls, among other groups of samples. The module will find genes that are significantly different between ccRCC and normal control samples.

*Jones J, Otu H, Spentzos D, Kolia S et al. Gene signatures of progression and metastasis in renal cell cancer. Clin Cancer Res 2005, 15;11(16):5730-9. PMID:16115910.

## Setup and Installation
The module was conducted using R version 4.4.2 in a Jupyter Notebook with R kernel. The associated Anaconda environment for the Jupyter Notebook is provided in file 'environment_R_4.4.2.yaml'. If you run the R code, ensure the proper version is installed on your machine along with the required libraries. To install Bioconductor packages, packge **BiocManager** is required. You can install it using the following line:

```R
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
```

Then use package BiocManager to install the Bioconductor packages needed in this module:

```R
if (!requireNamespace("GEOquery", quietly = TRUE)) BiocManager::install("GEOquery")
if (!requireNamespace("affy", quietly = TRUE)) BiocManager::install("affy")
if (!requireNamespace("R.utils", quietly = TRUE)) BiocManager::install("R.utils")
if (!requireNamespace("affyPLM", quietly = TRUE)) BiocManager::install("affyPLM")
if (!requireNamespace("genefilter", quietly = TRUE)) BiocManager::install("genefilter")
if (!requireNamespace("hgu133a.db", quietly = TRUE)) BiocManager::install("hgu133a.db")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db")
if (!requireNamespace("AnnotationDbi", quietly = TRUE)) BiocManager::install("AnnotationDbi")
if (!requireNamespace("limma", quietly = TRUE)) BiocManager::install("limma")
if (!requireNamespace("pheatmap", quietly = TRUE)) BiocManager::install("pheatmap")
```

If you prefer to work with Jupyter Notebook, import the yaml environment file into your Anaconda, install R and the needed packages as shown above. Alternatively, you can create your own conda environment. For example:
```
conda create -n R_4.4.2
```
Then activate the environment and install the desired version of R:
```
conda activate R_4.4.2
conda install -c conda-forge r-base=4.4.2 r-essentials -y 
```
To install the R kernel, you can start R in a terminal and use install.packages:
```
conda activate R_4.4.2
R
```
```R
> install.packages("IRkernel")
```

If installing any package directly from R using install.packages("package_name") proves to be problematic due to dependencies that require a specific package version, you may install such packages using conda after activating your environment. For example:
```
conda activate r_env
conda install -c r r-matrixStats=1.5.0
```

 All package versions used in this module are shown at the end of the Jupyter notebook.
 
