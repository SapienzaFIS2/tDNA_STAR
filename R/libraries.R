#### LIBRARIES #######
library(readr)
library(readxl)

library(dplyr, quietly = T)
library(tibble)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("GEOquery")
library(DESeq2, quietly = T)
library(GEOquery, quietly=T)


BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks, quietly = T)

library(SummarizedExperiment,quietly = T)
library(ggplot2)
install.packages("ggrepel")
library(ggrepel)
library(ggpubr)  
library(patchwork)


library(igraph)
library(network)










