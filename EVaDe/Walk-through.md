# Project Workflow Walk-through
#### This document provides a step-by-step guide on how to run the analysis pipeline of EVaDe, using primate Prefrontal Cortex (PFC) single-cell RNA-seq data as an example. 
---
## Step 1: Downsample the cells so that each cell type has an equal number of cells

**Script:** `Expression_Variance_Decomposition/cell_type_subsampling_equalize_counts.R`

**Usage Example:**

```R
# Assuming the R script is loaded or sourced, containing the function:
# cell_type_subsampling_equalize_counts(input_file, output_file, cells_per_type, species)

# Example command:
result <- cell_type_subsampling_equalize_counts(
  input_file = "/path/to/your/data/pfc.rds",         # Input seurat data
  output_file = "/path/to/your/output/pfc_rh.rds",   # Output seurat data
  cells_per_type = 1800,                             # Target cell number per type
  species = c("Human", "Rhesus")                     # Species to include
)
# Note:The selection of the cell number threshold should consider the trade-off between reserving more cell types for downstream analyses and retaining more cells for each type after down-sampling.
``` 

## Step 2: Expression Variance Decomposition in two species

**Script:** `Expression_Variance_Decomposition/expression_variance_decomposition_overall.R`

**Usage Example:**

```R
# Assuming the R script is loaded or sourced, containing the function:
# expression_variance_decomposition_overall(input_file, output_dir)

# Example command:
expression_variance_decomposition_overall(
  input_file = "/path/to/your/output/pfc_rh.rds",   # Input from Step 1
  output_dir = "/path/to/your/output/sum_frame_rh"  # Directory for results CSVs
)

# This will generate multiple CSV files (one per cell type) in the 'sum_frame_rh' directory,
# containing variance components (species, sample, residual) for each gene.
``` 

## Step 3: Expression Variance Decomposition for a specified species

**Script:** `Expression_Variance_Decomposition/expression_variance_decomposition_by_taxon.R`

**Usage Example:**

```R
# Assuming the R script is loaded or sourced, containing the function:
# expression_variance_decomposition_by_taxon(input_file, output_dir, species)

# Example command:
expression_variance_decomposition_by_taxon(
  input_file = "/path/to/your/output/pfc_rh.rds",   # Input from Step 1
  output_dir = "/path/to/your/output/sum_frame_rh", # Base directory for results
  species = c("Human", "Rhesus")                    # Species to analyze(1 or 2)
)

# This will generate CSV files in subdirectories 'sum_frame_rh/Human/' and 'sum_frame_rh/Rhesus/', containing a specified species variance components
# for each gene per cell type. The residual variance here represents the noise term.
``` 

## Step 4: Calculate cross-cell-type correlations for each gene (Identifing NC key genes step1)

**Script:** `Key_Gene_Identification/NC_calculate_correlation.py`

**Usage Example:**

```bash

python Key_Gene_Identification/NC_calculate_correlation.py \
    --cell_types_file /path/to/your/metadata/cell_types.csv \   #cell types CSV file
    --expression_file1 /path/to/expression/expression_Human.csv \   #species1 mean expression CSV file
    --expression_file2 /path/to/expression/expression_Rhesus.csv \  #species2 mean expression CSV file
    --anova_dir /path/to/your/output/sum_frame_rh/ \    #Directory containing 2species ANOVA results
    --anova_dir_species1 /path/to/your/output/sum_frame_rh/Human/ \ #Directory containing ANOVA results for species1
    --anova_dir_species2 /path/to/your/output/sum_frame_rh/Rhesus/ \    #Directory containing ANOVA results for species2
    --output_file1 /path/to/your/output/sum_frame_rh/Human/gogene_H_allpvalue.csv \ #Path to save the species1 output CSV file
    --output_file2 /path/to/your/output/sum_frame_rh/Rhesus/gogene_R_allpvalue.csv  #Path to save the species1 output CSV file

# The output files contain correlation coefficients and p-values for each gene.
``` 

## Step 5: Screen for significantly negatively correlated genes(Identifing NC key genes step2)

**Script:** `Key_Gene_Identification/PFC_NC_FDR_significance_negative.R`

**Usage Example:**

```R
# Assuming the R script is loaded or sourced, containing the function:
# NC_FDR_significance_negative(input_file, output_all_genes, output_key_genes)

# Example command (using Human results from Step 4):
result <- NC_FDR_significance_negative(
  input_file = "/path/to/your/output/sum_frame_rh/Human/gogene_H_allpvalue.csv", # Input from Step 4
  output_all_genes = "/path/to/your/output/NC_analysis/all-genes-fdr.csv",       # Output with all genes + FDR
  output_key_genes = "/path/to/your/output/NC_analysis/key-genes-fdr.csv"        # Output with significant NC genes
)

# This R script reads the output from Step 4, calculates the FDR for the correlation p-values.
# You might run this separately for the Rhesus results file as well.
``` 

## Step 6: Calculate the Divergence-to-Variation Ratio for each gene in each cell types and indentifies DVR key genes

**Script:** `Key_Gene_Identification/Calculate_DVR.py`

**Usage Example:**

```bash

python Key_Gene_Identification/Calculate_DVR.py \
  --cell_types_file /path/to/your/metadata/cell_types.csv \
  --expression_file1 /path/to/expression/human_expression.csv \
  --expression_file2 /path/to/expression/rhesus_expression.csv \
  --anova_dir /path/to/your/output/sum_frame_rh/ \
  --anova_dir_species1 /path/to/your/output/sum_frame_rh/Human/ \
  --anova_dir_species2 /path/to/your/output/sum_frame_rh/Rhesus/ \
  --output_file1 /path/to/your/output/DVR_analysis/human_dvr_results \
  --output_file2 /path/to/your/output/DVR_analysis/rhesus_dvr_results \
  --exp_threshold 0.01 \
  --dsp_v_threshold 0.005 \
  --dsam_v_threshold 0.05

# THRESHOLD PARAMERS:
# exp_threshold: Set appropriate cutoff to exclude lowly expressed genes, which may exhibit larger stochastic expression divergence or noise.
# dsp_v_threshold & dsam_v_threshold: Adjust to retain sufficient top-ranked genes for functional enrichment analysis.


# The output files contain DVR values of all genes and lists of identified DVR key genes.
``` 