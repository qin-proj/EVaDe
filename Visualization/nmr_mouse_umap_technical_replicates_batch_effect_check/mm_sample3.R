library(Seurat)
library(tidyverse)
setwd('bm_data/')
library(Matrix)
library(dplyr)
barcode.path <- "GSM6605115_mouse.BM.M_2_1_barcodes.tsv.gz"
features.path <- "GSM6605115_mouse.BM.M_2_1_genes.tsv.gz"
matrix.path <- "GSM6605115_mouse.BM.M_2_1_matrix.mtx.gz"
mat <- readMM(file = matrix.path)
feature.names <- read.table(features.path, header = FALSE, stringsAsFactors = FALSE)
barcode.names <- read.table(barcode.path, header = FALSE, stringsAsFactors = FALSE)
colnames(mat) <- barcode.names$V1
rownames(mat) <- feature.names$V1
mouse <- CreateSeuratObject(counts = mat)
current_cell_names <- Cells(mouse)
new_cell_names <- paste0(current_cell_names, "_mouse.BM.M_2.1")
mouse <- RenameCells(mouse, new.names = new_cell_names)
ct=read.csv('../mouse.csv',header=T)
cells_name=ct$cell_barcode
mouse_filtered <- subset(mouse, cells = cells_name)
rownames(ct)=ct$cell_barcode
mouse_filtered=AddMetaData(mouse_filtered,metadata=ct)
barcode.path <- "GSM6605116_mouse.BM.M_2_2_barcodes.tsv.gz"
features.path <- "GSM6605116_mouse.BM.M_2_2_genes.tsv.gz"
matrix.path <- "GSM6605116_mouse.BM.M_2_2_matrix.mtx.gz"
mat <- readMM(file = matrix.path)
feature.names <- read.table(features.path, header = FALSE, stringsAsFactors = FALSE)
barcode.names <- read.table(barcode.path, header = FALSE, stringsAsFactors = FALSE)
colnames(mat) <- barcode.names$V1
rownames(mat) <- feature.names$V1
mouse2 <- CreateSeuratObject(counts = mat)
current_cell_names <- Cells(mouse2)
new_cell_names <- paste0(current_cell_names, "_mouse.BM.M_2.2")
mouse2 <- RenameCells(mouse2, new.names = new_cell_names)
mouse2_filtered <- subset(mouse2, cells = cells_name)
mouse2_filtered=AddMetaData(mouse2_filtered,metadata=ct)
combined <- merge(mouse_filtered, y = mouse2_filtered,project = "combined")
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined)
combined <- ScaleData(combined)
combined <- RunPCA(combined)
combined <- RunUMAP(combined, dims = 1:30)
p1 <- DimPlot(combined, reduction = "umap", group.by = "cell_subset", label = TRUE, repel = TRUE) +
  ggtitle("UMAP by Cell Subset")
p2 <- DimPlot(combined, reduction = "umap", group.by = "sample", label = TRUE, repel = TRUE) +
  ggtitle("UMAP by Sample")
combined_plot <- p1 + p2
ggsave("mm4_combined_umap.png", combined_plot, width = 16, height = 8, dpi = 300)
saveRDS(combined,"mouse4.rds")
