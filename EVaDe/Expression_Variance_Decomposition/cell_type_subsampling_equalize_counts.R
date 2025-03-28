
# Function: cell_type_subsampling_equalize_counts
# Description: Subsample cells to equalize counts across cell types
# Parameters:
#   input_file: Path to input RDS file
#   output_file: Path to output RDS file
#   cells_per_type: Number of cells to keep per cell type(The selection of the cell number threshold should consider the trade-off between reserving more cell types for downstream analyses and retaining more cells for each type after down-sampling.)
#   species: Vector of species to include in the analysis

cell_type_subsampling_equalize_counts <- function(input_file, output_file, cells_per_type = 1800, species = c("Rhesus", "Human")) {

  library(Seurat)
  library(tidyverse)
  
  # Read the input RDS file
  pbmc <- readRDS(input_file)
  
  pbmc <- pbmc[, pbmc$species %in% species]
  
  # Get unique celltype
  classes <- unique(pbmc$subclass)
  
  # Initialize a list to store cells to keep
  cells_to_keep <- list()
 
  for(class in classes) {
    # Get cells of current subclass
    cells <- names(pbmc$cell_name[pbmc$subclass == class])
    
    # If cell count meets or exceeds cells_per_type, subsample，If cell count is less than cells_per_type, this class will be excluded
    if(length(cells) >= cells_per_type) {
      cells_to_keep[[class]] <- sample(cells, cells_per_type, replace = FALSE)
    }
  }
  
  # Combine all cells to keep
  all_cells_to_keep <- unlist(cells_to_keep)
  
  # Subset the Seurat object to keep only the selected cells
  pbmc <- pbmc[, colnames(pbmc) %in% all_cells_to_keep]
  
  # Save the processed data to output file
  saveRDS(pbmc, file = output_file)
  
  # Return the processed object (optional)
  return(pbmc)
}

# Usage example:
# result <- cell_type_subsampling_equalize_counts(
#   input_file = "path/to/input/pfc.rds",
#   output_file = "path/to/output/pfc_rh.rds",
#   cells_per_type = 1800,
#   species = c("Rhesus", "Human")
# )