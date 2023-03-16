# script to perform standard workflow steps to analyze single cell RNA-Seq data
# script: https://github.com/kpatel427/YouTubeTutorials/blob/main/singleCell_standard_workflow.R
# Youtube: https://youtu.be/5HBzgsz8qyk

# load libraries
library(Seurat)
library(tidyverse)
library(here)
setwd(here())

# Load the dataset
mtx_obj <- ReadMtx(mtx = "data/raw/GSE222510_raw_matrix.mtx.gz",
                   features = "data/raw/GSE222510_raw_genes.tsv.gz",
                   cells = "data/raw/GSE222510_raw_barcodes.tsv.gz", feature.column = 1)
seurat_mtx <- CreateSeuratObject(counts = mtx_obj, min.cells = 3, min.features = 200)
str(seurat_mtx)

# 1. QC -------
View(seurat_mtx@meta.data)

# % MT reads
seurat_mtx[["percent.mt"]] <- PercentageFeatureSet(seurat_mtx, pattern = "^MT-")
View(seurat_mtx@meta.data)

VlnPlot(seurat_mtx, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(seurat_mtx, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')

# 2. Filtering -----------------
seurat_mtx <- subset(seurat_mtx, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & 
                       percent.mt < 5)

# 3. Normalize data ----------
#nsclc.seurat.obj <- NormalizeData(nsclc.seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
# OR
seurat_mtx <- NormalizeData(seurat_mtx)

# 4. Identify highly variable features --------------
seurat_mtx <- FindVariableFeatures(seurat_mtx, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat_mtx), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat_mtx)
plot1 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1

# Save as png
png("figure/plot.png")
plot1
dev.off()
str(seurat_mtx)