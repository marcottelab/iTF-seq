#!/usr/bin/env Rscript

library(tidyverse)
library(Seurat)
library(sctransform)

###
data_matrix <- Read10X(data.dir="outs/filtered_feature_bc_matrix")
seurat_obj <- CreateSeuratObject(counts=data_matrix, min.cells=0, min.features=0, project="iTF_seq")
rm(data_matrix)

###
meta_ <- read.table("d1.metadata.g1k_mt10_umi10k.umi3.txt", sep="\t", header=TRUE, row.names=1) 
meta_ <- meta_[rownames(seurat_obj@meta.data), ]
seurat_obj <- AddMetaData(seurat_obj, meta_, col.name = colnames(meta_))
rm(meta_)

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern="^mt-")
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 1000 & percent.mt < 10 & nCount_RNA > 10000)
seurat_obj <- subset(seurat_obj, subset = detected_iTF != "NOT_USED")
seurat_obj <- subset(seurat_obj, subset = detected_iTF != "ZERO_BY_THR")

seurat_obj <- SCTransform(seurat_obj, method="glmGamPoi", vars.to.regress="percent.mt", return.only.var.genes=FALSE)

target_values <- c("NO_iTF_TAGS", "Arid3a", "Arnt", "Cdx2", "Cebpb", "Cited1", "Cited2", "Creb3", "Creb3l2", "Dlx3", "Dlx4", "Dlx5", "Dnajb6", "Dpf1", "E2f5", "Elf5", "Erf", "Esrrg", "Ets2", "Etv4", "Fbxl19", "Fos", "Fosl1", "Fosl2", "Foxn2", "Gata3", "Gata4", "Gcm1", "Gli1", "Gsc", "Hand1", "Hdac1", "Hdac2", "Hmgn2", "Hopx", "Id1", "Id2", "Jun", "Maff", "Mafk", "Max", "Med26", "Mef2d", "Meis1", "Mllt1", "Nkx2-1", "Nkx2-5", "Npm1", "Otx2", "Pax9", "Pdx1", "Pole3", "Polr2e", "Pou2f1", "Psmc3", "Rbfox2", "Rbpj", "Rbpjl", "Rybp", "Smarcb1", "Smarce1", "Snai1", "Taf6", "Tbx20", "Tbx5", "Tcea1", "Tcf7", "Tfap2c", "Tfdp1", "Tgif1", "Tle3", "Trim25", "Tsc22d3", "Vgll3", "Vgll4", "Yy1", "Zfp13", "Zfp36l1")
print(length(target_values))

for(i in target_values){
    temp_set <- subset(seurat_obj, subset = detected_iTF == i)
	my_df <- as.data.frame((GetAssayData(object=temp_set[["SCT"]], slot="scale.data")))
	###
    write.table(my_df, file=paste("day1/", i, ".tsv", sep=""), sep="\t", row.names=TRUE)
}
