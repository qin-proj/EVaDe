# Function: expression_variance_decomposition_by_taxon
# Description: Perform expression variance decomposition for each cell type for a specified species
# Parameters:
#   input_file: Path to input Seurat RDS file
#   output_dir: Directory to save output CSV files
#   species: Species to analyze ("Human" or "Rhesus")

expression_variance_decomposition_by_taxon <- function(input_file, output_dir, species = c("Human", "Rhesus")) {

  library(Seurat)
  library(tidyverse)
  
  # Read the input Seurat RDS file
  pbmc <- readRDS(input_file)
  
  # Get unique cell types from the data
  cell_types <- unique(pbmc$subclass)
  
  # Ensure output directories exist for each species
  for (sp in species) {
    dir.create(file.path(output_dir, sp), recursive = TRUE, showWarnings = FALSE)
  }
  
  # Function to process data for a given cell type and species
  process_cell_type_species <- function(cell_type, sp) {
    # Subset data for the given cell type and species
    vcell <- pbmc@assays$RNA@data[, pbmc$subclass == cell_type & pbmc$species == sp]
    
    # Check if there are cells for this combination
    if (ncol(vcell) == 0) {
      message(paste("No cells found for", cell_type, "in species", sp))
      return()
    }
    
    # Prepare data for analysis
    exp <- vcell[, 1:ncol(vcell)]
    dimnames <- list(rownames(exp), colnames(exp))
    data <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
    data2 <- t(data)
    
    # Get sample names
    sam <- pbmc@meta.data[rownames(data2), ]$samplename
    
    # Create data frame for analysis
    data3 <- data.frame(sam, data2)
    
    # Get gene list
    gene_list <- setdiff(colnames(data3), "sam")
    
    # Initialize results frame
    frame <- data.frame(aa = c(0, 0))
    rownames(frame) <- c("sam", "Residuals")
    
    # Perform ANOVA for each gene
    for (ge in gene_list) {
      a1 <- aov(get(ge) ~ sam, data3)
      a2 <- summary(a1)[[1]]["Mean Sq"]
      names(a2)[1] <- ge
      frame[ge] <- a2
    }
    
    # Write results to CSV
    write_csv(frame, file = file.path(output_dir, sp, paste0(cell_type, ".csv")))
  }
  
  # Process each cell type for each specified species
  for (sp in species) {
    for (cl in cell_types) {
      process_cell_type_species(cl, sp)
    }
  }
}

# Usage example:
# expression_variance_decomposition_by_taxon(
#   input_file = "path/to/input/pfc_rh.rds",
#   output_dir = "path/to/output/",
#   species = c("Human", "Rhesus")  # Can specify one or both species
# )
