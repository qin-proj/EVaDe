library(Seurat)
library(tidyverse)
setwd('bm_data/')
library(Matrix)
library(dplyr)
barcode.path <- "GSM6605103_nmr.BM.F_2_1_barcodes.tsv.gz"
features.path <- "GSM6605103_nmr.BM.F_2_1_genes.tsv.gz"
matrix.path <- "GSM6605103_nmr.BM.F_2_1_matrix.mtx.gz"
mat <- readMM(file = matrix.path)
feature.names <- read.table(features.path, header = FALSE, stringsAsFactors = FALSE)
barcode.names <- read.table(barcode.path, header = FALSE, stringsAsFactors = FALSE)
colnames(mat) <- barcode.names$V1
rownames(mat) <- feature.names$V1
nmr <- CreateSeuratObject(counts = mat)
current_cell_names <- Cells(nmr)
new_cell_names <- paste0(current_cell_names, "_nmr.BM.F_2.1")
nmr <- RenameCells(nmr, new.names = new_cell_names)
ct=read.csv('../nmr.csv',header=T)
cells_name=ct$cell_barcode
nmr_filtered <- subset(nmr, cells = cells_name)
rownames(ct)=ct$cell_barcode
nmr_filtered=AddMetaData(nmr_filtered,metadata=ct)
barcode.path <- "GSM6605104_nmr.BM.F_2_2_barcodes.tsv.gz"
features.path <- "GSM6605104_nmr.BM.F_2_2_genes.tsv.gz"
matrix.path <- "GSM6605104_nmr.BM.F_2_2_matrix.mtx.gz"
mat <- readMM(file = matrix.path)
feature.names <- read.table(features.path, header = FALSE, stringsAsFactors = FALSE)
barcode.names <- read.table(barcode.path, header = FALSE, stringsAsFactors = FALSE)
colnames(mat) <- barcode.names$V1
rownames(mat) <- feature.names$V1
nmr2 <- CreateSeuratObject(counts = mat)
current_cell_names <- Cells(nmr2)
new_cell_names <- paste0(current_cell_names, "_nmr.BM.F_2.2")
nmr2 <- RenameCells(nmr2, new.names = new_cell_names)
nmr2_filtered <- subset(nmr2, cells = cells_name)
nmr2_filtered=AddMetaData(nmr2_filtered,metadata=ct)
combined <- merge(nmr_filtered, y = nmr2_filtered,project = "combined")
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
ggsave("2combined_umap.png", combined_plot, width = 16, height = 8, dpi = 300)
saveRDS(combined,"nmr2.rds")
