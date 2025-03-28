# Function: PFC_NC_FDR_significance_negative
# Description: Calculate FDR, and filter significant genes
# Parameters:
#   input_file: Path to input correlation p-value CSV file
#   output_all_genes: Path to output CSV file for all genes with FDR
#   output_key_genes: Path to output CSV file for key genes

PFC_NC_FDR_significance_negative <- function(input_file, output_all_genes, output_key_genes) {

  library(dplyr)
  
  # Read input correlation CSV file
  g01 <- read.csv(input_file, header = TRUE)
  
  # Calculate FDR for p1, p2, and p3
  calculate_fdr <- function(p_values) {
    p.adjust(p_values, method = "BH")
  }
  
  g04 <- g01 %>%
    mutate(
      FDR1 = calculate_fdr(p1),
      FDR2 = calculate_fdr(p2),
      FDR3 = calculate_fdr(p3)
    )
  
  # Write all genes with FDR to CSV
  write.csv(g04, file = output_all_genes, row.names = FALSE)
  
  # Filter significant genes
  x01 <- g04 %>% filter(FDR2 < 0.05 & Rps < 0)
  x02 <- g04 %>% filter(FDR1 < 0.05 & Rpm > 0)
  x03 <- g04 %>% filter(FDR3 < 0.05 & Rms < 0)
  
  # Find key genes using anti_join
  filtered_rows <- anti_join(x01, rbind(x02, x03))
  
  # Write key genes to CSV
  write.csv(filtered_rows, file = output_key_genes, row.names = FALSE)
  
  # Return the filtered rows (optional)
  return(filtered_rows)
}

# Usage example:
# result <- PFC_NC_FDR_significance_negative(
#   input_file = "path/to/input/gogene_H_allpvalue.csv",
#   output_all_genes = "path/to/output/all-genes-fdr.csv",
#   output_key_genes = "path/to/output/key-genes-fdr.csv",
# )
