"""
This script calculates the Divergence-to-Noise Ratio for each gene in each cell type and indentifies DVR key genes
"""


import argparse
import pandas as pd
import scipy
from scipy import stats
import numpy as np

def main(args):
    # Read cell types from file
    all_cell_types = pd.read_csv(args.cell_types_file)
    cp_list = all_cell_types['cell_type'].tolist()
    cp_list2 = cp_list[1:]  # All cell types except the first one

    # Read expression data
    exp = pd.read_csv(args.expression_file1, index_col=0)
    exp2 = pd.read_csv(args.expression_file2, index_col=0)

    # Read and process ANOVA results
    df0 = read_and_process_anova(args.anova_dir, cp_list[0], cp_list2)
    df1 = read_and_process_species_anova(args.anova_dir_species1, cp_list[0], cp_list2)
    df2 = read_and_process_species_anova(args.anova_dir_species2, cp_list[0], cp_list2)

    # Get gene list from dataframe columns
    gene_list = [x for x in df0.columns if x != "aa"]

    # Calculate correlations and save results for both species
    calculate_and_save_dvr(gene_list, cp_list, df0, df1, exp, args.output_file1, args.exp_threshold, args.dsp_v_threshold, args.dsam_v_threshold)
    calculate_and_save_dvr(gene_list, cp_list, df0, df2, exp2, args.output_file2, args.exp_threshold, args.dsp_v_threshold, args.dsam_v_threshold)

def read_and_process_anova(anova_dir, first_cp, other_cp_list):
    """
    Read and process ANOVA results for all cell types
    
    Returns:
    DataFrame: Processed ANOVA results with three rows per cell type
    """
    df = pd.read_csv(f"{anova_dir}/{first_cp}.csv", header=0)
    df = df.rename(index={0: f"{first_cp}--0", 1: f"{first_cp}--1", 2: f"{first_cp}--2"})
    
    for cp in other_cp_list:
        temp_df = pd.read_csv(f"{anova_dir}/{cp}.csv", header=0)
        temp_df = temp_df.rename(index={0: f"{cp}--0", 1: f"{cp}--1", 2: f"{cp}--2"})
        df = pd.concat([df, temp_df], axis=0, join='inner')
    
    return df

def read_and_process_species_anova(anova_dir, first_cp, other_cp_list):
    """
    Read and process species-specific ANOVA results for all cell types
    
    Returns:
    DataFrame: Processed ANOVA results with only the second row extracted and renamed
    """
    df = pd.read_csv(f"{anova_dir}/{first_cp}.csv", header=0)
    # Extract only the second row (index 1) and rename it
    df = df.loc[[1]]
    df = df.rename(index={1: f"{first_cp}--4"})
    
    for cp in other_cp_list:
        temp_df = pd.read_csv(f"{anova_dir}/{cp}.csv", header=0)
        # Extract only the second row (index 1) and rename it
        temp_df = temp_df.loc[[1]]
        temp_df = temp_df.rename(index={1: f"{cp}--4"})
        df = pd.concat([df, temp_df], axis=0, join='inner')
    
    return df

def calculate_and_save_dvr(gene_list, cp_list, df0, df_species, exp, output_file, exp_threshold=0.01, dsp_v_threshold=0.005,dsam_v_threshold=0.05):
    """
    Calculate Divergence-to-Variation Ratio and save results.

    """
    # Initialize lists to store results
    genes, cell_types, dsp_v_ratios, dsam_v_ratios = [], [], [], []
    
    # Calculate ratios for each gene and cell type
    for gene in gene_list:
        for cell_type in cp_list:
            # Check if gene meets criteria for analysis
            if (exp[gene][cell_type] > exp_threshold and 
                df0[gene][f'{cell_type}--0'] > 0 and 
                df_species[gene][f'{cell_type}--4'] > 0):
             
                genes.append(gene)
                cell_types.append(cell_type)
                
                # Calculate ratios
                # Dsp/V
                dsp_v_ratios.append(df0[gene][f'{cell_type}--0'] / df_species[gene][f'{cell_type}--4'])
                # Dsam/V
                dsam_v_ratios.append(df0[gene][f'{cell_type}--1'] / df_species[gene][f'{cell_type}--4'])
    
    # Create DataFrame with results
    result_df = pd.DataFrame({
        'gene': genes,
        'celltype': cell_types,
        'Dsp/V': dsp_v_ratios,
        'Dsam/V': dsam_v_ratios
    })
    
    # Add expression data
    expression_values = []
    for i, row in result_df.iterrows():
        gene = row['gene']
        cell_type = row['celltype']
        value = exp[gene][cell_type]
        expression_values.append(value)
    
    result_df['exp'] = expression_values
    
    # Save complete results
    result_df.to_csv(f"{output_file}_all.csv", index=False)
    
    # Sort and filter results
    sorted_by_dsp = result_df.sort_values(by='Dsp/V', ascending=False)
    sorted_by_dsam = result_df.sort_values(by='Dsam/V', ascending=False)
    
    # Get top 0.5% genes by Dsp/V
    top_percent_rows = int(len(sorted_by_dsp) * dsp_v_threshold)
    top_genes = sorted_by_dsp.head(top_percent_rows)
    
    # Get bottom 5% genes by Dsam/V
    bottom_percent_rows = int(len(sorted_by_dsam) * dsam_v_threshold)
    
    # Filter out genes that are in both top Dsp/V and bottom Dsam/V
    filtered_genes = top_genes[~top_genes.isin(sorted_by_dsam.head(bottom_percent_rows))].dropna()
    
    # Save filtered results
    filtered_genes.to_csv(f"{output_file}_filtered.csv", index=False)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Identify DVR key genes')
    
    # Input files
    parser.add_argument('--cell_types_file', required=True, help='Path to the cell types CSV file')
    parser.add_argument('--expression_file1', required=True, help='Path to the species1 mean expression CSV file')
    parser.add_argument('--expression_file2', required=True, help='Path to the species2 mean expression CSV file')
    
    # ANOVA directories
    parser.add_argument('--anova_dir', required=True, help='Directory containing 2species ANOVA results')
    parser.add_argument('--anova_dir_species1', required=True, help='Directory containing ANOVA results for species1')
    parser.add_argument('--anova_dir_species2', required=True, help='Directory containing ANOVA results for species2')
    
    # Output files
    parser.add_argument('--output_file1', required=True, help='Path to save the species1 output CSV file')
    parser.add_argument('--output_file2', required=True, help='Path to save the species2 output CSV file')
    
    # Parameters
    # exp_threshold: Set an appropriate cutoff to filter out lowly expressed genes potentially with larger stochastic expression divergence or noise.
    # dsp_v_threshold & dsam_v_threshold: Adjust to retain sufficient top-ranked genes for functional enrichment analysis

    parser.add_argument('--exp_threshold', type=float, default=0.01, help='Expression threshold (default: 0.01)')
    parser.add_argument('--dsp_v_threshold', type=float, default=0.005, help='Threshold for top Dsp/V ratio percentile (default: 0.005)')
    parser.add_argument('--dsam_v_threshold', type=float, default=0.05, help='Threshold for bottom Dsam/V ratio percentile (default: 0.05)')
    
    args = parser.parse_args()
    main(args)


