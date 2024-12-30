# Unveiling cell-type-specific mode of evolution in comparative single-cell expression data
### General Overview
![图片描述](https://github.com/qin-proj/pic1/blob/main/fig1_1.jpg)
We propose an analytical framework, grounded in phenotypic evolution theory, for the comparative analysis of single-cell expression data. Given a single-cell expression dataset comprising two species or populations, our framework decomposes the total expression variance into the divergence and noise for each gene in each cell type. Comparing across cell types and genes, we introduce two distinct approaches to identify genes exhibiting large evolutionary divergence and small expression noise in certain cell types. The first approach, termed the Negative Correlation (NC) strategy, leverages the negative correlation between inter-taxon divergence($D_{sp}^{n}$) and noise in one specific taxon ($V_{taxon}^{n}$). The second approach, termed DVR strategy, is based on the ratio of D<sub>sp</sub> to V<sub>taxon</sub>.
We apply these methodologies to two single-cell expression datasets: the primate prefrontal cortex (PFC) and naked mole-rat bone marrow(BM).


### Folder Descriptions

**1. Expression_Variance_Decomposition**
   
This folder contains scripts for decomposing the variance of gene expression matrices from PFC and BM single-cell data. This includes:

*   Subsampling of cell types to ensure that the number of each cell type is consistent: `PFC_cell_type_subsampling_equalize_counts.R` (PFC) and `BM_cell_type_subsampling_equalize_counts.R` (BM).
*   Overall variance decomposition (ANOVA) for each cell type between-taxon: `PFC_expression_variance_decomposition_overall.R` (PFC) and `BM_expression_variance_decomposition_overall.R` (BM).
*   Variance decomposition (ANOVA) for each celltype in one specific focal taxon : `PFC_expression_variance_decomposition_by_taxon.R` (PFC) and `BM_expression_variance_decomposition_by_taxon.R` (BM).

**2. Key_Gene_Identification**

This folder contains scripts for identifying key genes in both PFC and BM datasets using NC and DVR strategies. This includes:

*   NC analysis:
    *   Correlation calculation: `PFC_NC_calculate_correlation.py` (PFC) and `BM_NC_calculate_correlation.py` (BM).
    *   FDR significance testing: `PFC_NC_FDR_significance_negative.py` (PFC) and `BM_NC_FDR_significance_negative.R` (BM).
*   DVR implementation: `PFC_calculate_DVR.py` (PFC) and `BM_calculate_DVR.py` (BM).

**3. RECNE_Nearest_Gene_Analysis**

This folder contains scripts to identify protein-coding genes located within 100kb flanking regions of rapidly evolving conserved noncoding elements (RECNEs). It includes:

*   `extract_overall_cds_regions_from_gff.pl`: Extracts coding sequence (CDS) regions from a GFF file.
*   `gene_id_matcher.pl`: Matches gene IDs between different annotation versions.
*   `output_gene_pairs_with_min_distances.pl`: Calculates and outputs gene pairs with minimum distances, specifically identifying genes near RECNEs.

**4. Visualization**

This folder contains scripts used for generating figures and visualizations for the project.
