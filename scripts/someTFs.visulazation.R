#!/usr/bin/env Rscript
library(tidyverse)
library(Seurat)
library(sctransform)

combined_data <- readRDS(file="int.SCT.after_clustering.rds")

meta_new <- read.table("metadata_for_some_TFs.tsv", sep="\t", header=T, row.names=1)
meta_new <- meta_new[rownames(combined_data@meta.data), ]
combined_data <- AddMetaData(combined_data, meta_new, col.name=colnames(meta_new))
rm(meta_new)

pdf(file="someTFs.1.pdf", width=7, height=7)
color <- c("Day1_Control"="Red", "Day3_Control"="Orange", "Day5_Control"="Brown", "others"="Gray")
DimPlot(combined_data, reduction="umap", group.by="Example_Control_only", order=c("Day1_Control", "Day3_Control", "Day5_Control"), cols=color)
dev.off()

pdf(file="someTFs.2.pdf", width=21, height=7)
DimPlot(combined_data, reduction="umap", group.by="iTF_of_interest", split.by="Day", cols=c("Control"="Purple", "Gata3"="Blue", "Fosl1"="Green", "Tbx5"="Cyan", "Cited1"="Red", "Id2"="Brown", "Snai1"="orange", "others"="gray"), order=c("Control", "Gata3", "Fosl1", "Tbx5", "Cited1", "Id2", "Snai1"))
dev.off()
