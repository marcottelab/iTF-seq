#!/usr/bin/env Rscript
#Description: Integrate scRNA-seq data from day1/3/5 using SCTransform
#Usage: [THIS SCRIPT]

library(tidyverse)
library(Seurat)
library(sctransform)

data_d1 <- Read10X(data.dir="80TFs1d/outs/filtered_feature_bc_matrix")
data_d3 <- Read10X(data.dir="80TFs3d/outs/filtered_feature_bc_matrix")
data_d5 <- Read10X(data.dir="80TFs5d/outs/filtered_feature_bc_matrix")

obj_d1 <- CreateSeuratObject(counts=data_d1, min.cells=0, min.features=0, project="Day1")
rm(data_d1)
obj_d3 <- CreateSeuratObject(counts=data_d3, min.cells=0, min.features=0, project="Day3")
rm(data_d3)
obj_d5 <- CreateSeuratObject(counts=data_d5, min.cells=0, min.features=0, project="Day5")
rm(data_d5)

meta_d1 <- read.table("d1.metadata.g1k_mt10_umi10k.umi3.txt", sep="\t", header=TRUE, row.names=1)
meta_d3 <- read.table("d3.metadata.g1k_mt10_umi10k.umi3.txt", sep="\t", header=TRUE, row.names=1)
meta_d5 <- read.table("d5.metadata.g1k_mt10_umi10k.umi3.txt", sep="\t", header=TRUE, row.names=1)

meta_d1 <- meta_d1[rownames(obj_d1@meta.data), ]
meta_d3 <- meta_d3[rownames(obj_d3@meta.data), ]
meta_d5 <- meta_d5[rownames(obj_d5@meta.data), ]

obj_d1 <- AddMetaData(obj_d1, meta_d1, col.name=colnames(meta_d1))
obj_d3 <- AddMetaData(obj_d3, meta_d3, col.name=colnames(meta_d3))
obj_d5 <- AddMetaData(obj_d5, meta_d5, col.name=colnames(meta_d5))
rm(meta_d1, meta_d3, meta_d5)

obj_d1[["percent.mt"]] <- PercentageFeatureSet(obj_d1, pattern="^mt-")
obj_d3[["percent.mt"]] <- PercentageFeatureSet(obj_d3, pattern="^mt-")
obj_d5[["percent.mt"]] <- PercentageFeatureSet(obj_d5, pattern="^mt-")

obj_d1 <- subset(obj_d1, subset=nFeature_RNA>1000 & percent.mt<10 & nCount_RNA>10000)
obj_d3 <- subset(obj_d3, subset=nFeature_RNA>1000 & percent.mt<10 & nCount_RNA>10000)
obj_d5 <- subset(obj_d5, subset=nFeature_RNA>1000 & percent.mt<10 & nCount_RNA>10000)

obj_d1 <- subset(obj_d1, subset=detected_iTF != "NOT_USED")
obj_d1 <- subset(obj_d1, subset=detected_iTF != "ZERO_BY_THR")

obj_d3 <- subset(obj_d3, subset=detected_iTF != "NOT_USED")
obj_d3 <- subset(obj_d3, subset=detected_iTF != "ZERO_BY_THR")

obj_d5 <- subset(obj_d5, subset=detected_iTF != "NOT_USED")
obj_d5 <- subset(obj_d5, subset=detected_iTF != "ZERO_BY_THR")

obj_d1 <- AddMetaData(obj_d1, metadata="Day 1", col.name="Day")
obj_d3 <- AddMetaData(obj_d3, metadata="Day 3", col.name="Day")
obj_d5 <- AddMetaData(obj_d5, metadata="Day 5", col.name="Day")

obj_d1 <- RenameCells(object=obj_d1, add.cell.id="d1")
obj_d3 <- RenameCells(object=obj_d3, add.cell.id="d3")
obj_d5 <- RenameCells(object=obj_d5, add.cell.id="d5")

my_list <- list(obj_d1, obj_d3, obj_d5)
my_list <- lapply(X=my_list, FUN=SCTransform, method="glmGamPoi", vars.to.regress="percent.mt", return.only.var.genes=FALSE)

features <- SelectIntegrationFeatures(object.list=my_list)
my_list <- PrepSCTIntegration(object.list=my_list, anchor.features=features)
my_anchors <- FindIntegrationAnchors(object.list=my_list, normalization.method="SCT", anchor.features=features)
combined_data <- IntegrateData(anchorset=my_anchors, normalization.method="SCT")

saveRDS(combined_data, file="int.SCT.after_integration.rds")

combined_data <- RunPCA(combined_data, npcs=50, verbose=F)
combined_data <- RunUMAP(combined_data, reduction="pca", dims=1:50)
combined_data <- FindNeighbors(combined_data, reduction="pca", dims=1:50)
combined_data <- FindClusters(combined_data, resolution=1.5)

saveRDS(combined_data, file="int.SCT.after_clustering.rds")

#combined_data <- readRDS(file="int.SCT.after_clustering.rds")
write.table(combined_data@meta.data, file="int.SCT.metadata.tsv", quote=F, sep="\t", row.names=T, col.names=T)

pdf(file="integrated_plot.pdf", width=21, height=7)
p1 <- DimPlot(combined_data, reduction="umap", group.by="Day", cols=c("Brown", "#FFC107", "#1A85FF"))
p2 <- DimPlot(combined_data, reduction="umap", label=TRUE)
p1 + p2
#DimPlot(combined_data, reduction="umap", split.by="Day")
DimPlot(combined_data, reduction="umap", group.by="Day", split.by="Day", cols=c("Brown", "#FFC107", "#1A85FF"))
dev.off()

