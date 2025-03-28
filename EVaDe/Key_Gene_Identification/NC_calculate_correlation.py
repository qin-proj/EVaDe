"""
Calculate correlation coefficients between inter-species variance, inter-sample variance, and within cell-type noise across different cell types for every gene.

This script calculates correlations between:
- Inter-species variance vs Inter-sample variance (Rpm)
- Inter-species variance vs Within cell-type noise (Rps)
- Inter-sample variance vs Within cell-type noise (Rms)

"""


import argparse
import pyreadr
import pandas as pd
import scipy
from scipy import stats
import numpy as np
from matplotlib import pyplot

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
    calculate_and_save_correlations(gene_list, cp_list, df0, df1, exp, args.output_file1)
    calculate_and_save_correlations(gene_list, cp_list, df0, df2, exp2, args.output_file2)

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

def calculate_and_save_correlations(gene_list, cp_list, df0, df1, exp, output_file):
    """
    Calculate correlations between different components and save results

    """
    ge_li, Rpm, Rps, Rms, pv1, pv2, pv3 = [], [], [], [], [], [], []

    for ge in gene_list:
        x, y, z = [], [], []
        for cp in cp_list:
            if (df0[ge][f'{cp}--0'] > 0) and (df1[ge][f'{cp}--4'] > 0) and (exp[ge][cp] > 0):
                x.append(df0[ge][f'{cp}--0'] / exp[ge][cp])
                y.append(df0[ge][f'{cp}--1'] / exp[ge][cp])
                z.append(df1[ge][f'{cp}--4'] / exp[ge][cp])
        
        # Calculate correlations if enough data points
        if len(x) > 5:
            # Calculate correlation between species and sample components
            slope, intercept, r1, p1, std_err = stats.linregress(x, y)
            # Calculate correlation between species and residual components
            slope, intercept, r2, p2, std_err = stats.linregress(x, z)
            # Calculate correlation between sample and residual components
            slope, intercept, r3, p3, std_err = stats.linregress(y, z)
            
            ge_li.append(ge)
            Rpm.append(r1)
            pv1.append(p1)
            Rps.append(r2)
            pv2.append(p2)
            Rms.append(r3)
            pv3.append(p3)

    # Create and save results dataframe
    corr_df = pd.DataFrame({
        'gene': ge_li,
        'Rpm': Rpm,  # Correlation between species and sample components
        'p1': pv1,   # p-value for Rpm
        'Rps': Rps,  # Correlation between species and residual components
        'p2': pv2,   # p-value for Rps
        'Rms': Rms,  # Correlation between sample and residual components
        'p3': pv3,   # p-value for Rms
    })

    corr_df.to_csv(output_file, index=False)

if __name__ == "__main__":
    # Set up command line argument parser
    parser = argparse.ArgumentParser(description="Calculate correlations")
    parser.add_argument("--cell_types_file", required=True, help="Path to the cell types CSV file")
    parser.add_argument("--expression_file1", required=True, help="Path to the species1 mean expression CSV file")
    parser.add_argument("--expression_file2", required=True, help="Path to the species2 mean expression CSV file")
    parser.add_argument("--anova_dir", required=True, help="Directory containing 2species ANOVA results")
    parser.add_argument("--anova_dir_species1", required=True, help="Directory containing ANOVA results for species1")
    parser.add_argument("--anova_dir_species2", required=True, help="Directory containing ANOVA results for species2")
    parser.add_argument("--output_file1", required=True, help="Path to save the species1 output CSV file")
    parser.add_argument("--output_file2", required=True, help="Path to save the species2 output CSV file")
    
    args = parser.parse_args()
    main(args)


