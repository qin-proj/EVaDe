# Function: expression_variance_decomposition_overall
# Description: Perform expression variance decomposition for each cell type, analyzing both species and sample effects
# Parameters:
#   input_file: Path to input Seurat RDS file
#   output_dir: Directory to save output CSV files

expression_variance_decomposition_overall <- function(input_file, output_dir) {

  library(Seurat)
  library(tidyverse)
  
  # Read the input RDS file
  pbmc <- readRDS(input_file)
  DefaultAssay(pbmc) <- "RNA"
  
  # Get unique cell types from the data
  cell_types <- unique(pbmc$subclass)
  
  # Ensure output directory exists
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Process each cell type
  for(cl in cell_types) {
    message(paste("Processing cell type:", cl))
    
    # Extract data for current cell type
    vcell <- pbmc@assays$RNA@data[, pbmc$subclass == cl]
    
    # Prepare expression matrix
    exp <- vcell[, 1:ncol(vcell)]
    dimnames <- list(rownames(exp), colnames(exp))
    data <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
    data2 <- t(data)
    
    # Get sample and species information
    sam <- pbmc@meta.data[rownames(data2), ]$samplename
    sp <- pbmc@meta.data[rownames(data2), ]$species
    
    # Create data frame for analysis
    data3 <- data.frame(sam, data2)
    gene_list <- setdiff(colnames(data3), "sam")
    data4 <- data.frame(sp, data3)
    
    # Initialize results frame
    frame <- data.frame(aa = c(0, 0, 0))
    rownames(frame) <- c("sp", "sam", "Residuals")
    
    # Perform ANOVA for each gene
    for(ge in gene_list) {
        a1 <- aov(get(ge) ~ sp + sam, data4)
        a2 <- summary(a1)[[1]]["Mean Sq"]
        names(a2)[1] <- ge
        frame[ge] <- a2
    }
    
    # Save results
    output_file <- file.path(output_dir, paste0(cl, ".csv"))
    write_csv(frame, file = output_file)
    message(paste("Results saved to:", output_file))
  }
  
  message("Analysis completed successfully!")
}

# Usage example:
# expression_variance_decomposition_overall(
#   input_file = "path/to/pfc_rh.rds",
#   output_dir = "path/to/sum_frame_rh"
# )

		
