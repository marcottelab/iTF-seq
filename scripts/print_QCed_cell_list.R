library(Seurat)

# day 1
data_matrix <- Read10X(data.dir="outs/filtered_feature_bc_matrix")
seurat_obj <- CreateSeuratObject(counts=data_matrix, min.cells=0, min.features=0, project="iTF_seq")
rm(data_matrix)
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern="^mt-")
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 1000 & percent.mt < 10 & nCount_RNA > 10000)
write.table(seurat_obj@meta.data, "d1.g1k_mt10_umi10k.tsv", sep="\t", quote=F)

# day 3
data_matrix <- Read10X(data.dir="outs/filtered_feature_bc_matrix")
seurat_obj <- CreateSeuratObject(counts=data_matrix, min.cells=0, min.features=0, project="iTF_seq")
rm(data_matrix)
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern="^mt-")
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 1000 & percent.mt < 10 & nCount_RNA > 10000)
write.table(seurat_obj@meta.data, "d3.g1k_mt10_umi10k.tsv", sep="\t", quote=F)

# day 5
data_matrix <- Read10X(data.dir="outs/filtered_feature_bc_matrix")
seurat_obj <- CreateSeuratObject(counts=data_matrix, min.cells=0, min.features=0, project="iTF_seq")
rm(data_matrix)
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern="^mt-")
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 1000 & percent.mt < 10 & nCount_RNA > 10000)
write.table(seurat_obj@meta.data, "d5.g1k_mt10_umi10k.tsv", sep="\t", quote=F)
