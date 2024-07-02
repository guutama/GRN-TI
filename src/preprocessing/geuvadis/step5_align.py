import sys
import os
from pathlib import Path
import pandas as pd
import numpy as np
import random
from utils import sort_dataframe_by_category, rank_inverse_normal_transformation
from path_config import COMPRESSED_PATH, ALIGNED_PATH


# Append directories to sys.path for module imports
directories_to_append = [
    os.path.abspath(os.path.join(os.path.dirname(__file__), '.')),
    os.path.abspath(os.path.join(os.path.dirname(__file__), '..')),
    # Add more directories as needed
]
for directory in directories_to_append:
    sys.path.append(directory)

def filter_genes(exp, map, have_eQTL):
    """
    Filter genes based on their presence in the eqtl mapping data.

    Parameters:
    exp (DataFrame): DataFrame containing expression data.
    map (DataFrame): DataFrame containing mapping data.
    have_eQTL (bool): Whether to filter genes that have eQTL or not.

    Returns:
    DataFrame: Filtered genes DataFrame.
    """
    if have_eQTL:
        genes = exp[exp["GENE_ID"].isin(map["GENE_ID"])]
        genes = genes.reset_index(drop=True)
    else:
        genes = exp[~exp["GENE_ID"].isin(map["GENE_ID"])]
        genes = genes.reset_index(drop=True)
    return genes

def align_data(in_path, out_path):
    """
    Align data from expression, genotype, and eqtl mapping files, and save the processed data.

    Parameters:
    in_path (str): Directory path where the input files are located.
    out_path (str): Directory path where the processed files will be saved.
    """
    exp_file = os.path.join(in_path, "expression.csv.gz")
    geno_file = os.path.join(in_path, 'genotypes_all.csv.gz')
    mapping_file = os.path.join(in_path, 'mappings_all.csv.gz')
    
    # Check if the mandatory files exist
    if os.path.exists(exp_file) and os.path.exists(geno_file):
        # Load the mandatory files
        expressions_df = pd.read_csv(exp_file, compression='gzip')
        genotypes_df = pd.read_csv(geno_file, compression='gzip')
        print("Mandatory files loaded successfully.")
    else:
        # Raise an exception if one or both mandatory files are missing
        if not os.path.exists(exp_file):
            raise FileNotFoundError("Error: expressions.csv.gz is missing.")
        if not os.path.exists(geno_file):
            raise FileNotFoundError("Error: genotypes.csv.gz is missing.")

    # Check and load the optional file
    if os.path.exists(mapping_file):
        # Load the optional file
        mapping_df = pd.read_csv(mapping_file, compression='gzip')
        print("Optional file loaded successfully.")
    else:
        mapping_df = None
        print("Optional file is missing, skipping loading of the optional file.")
    
    if mapping_df is not None:
        top_mapping = mapping_df.loc[mapping_df.groupby('GENE_ID')['rvalue'].idxmax()]
        top_mapping = top_mapping.reset_index(drop=True)
        top_mapping = top_mapping.drop_duplicates("SNP_ID")

        top_snp = genotypes_df[genotypes_df["SNP_ID"].isin(top_mapping["SNP_ID"])]
        top_mapping = top_mapping[top_mapping['SNP_ID'].isin(top_snp["SNP_ID"])]
        top_mapping = top_mapping.reset_index(drop=True)
        top_snp = top_snp.reset_index(drop=True)
        print(top_mapping)
        print(top_snp)

        top_genes = filter_genes(exp=expressions_df, map=top_mapping, have_eQTL=True)
        b_genes = filter_genes(exp=expressions_df, map=top_mapping, have_eQTL=False)

        top_mapping = sort_dataframe_by_category(dataframe=top_mapping, value_col="GENE_ID", categories_list=top_genes["GENE_ID"].to_list())
        top_snp = sort_dataframe_by_category(dataframe=top_snp, value_col="SNP_ID", categories_list=top_mapping["SNP_ID"].to_list())
        all_genes = pd.concat([top_genes, b_genes], ignore_index=True).drop_duplicates()  
        all_snp = genotypes_df.copy()
        all_mapping = mapping_df.copy()
    
    dataframes = {
        "top_snp": os.path.join(out_path, "top_snp.csv.gz"),
        "all_snp": os.path.join(out_path, "all_snp.csv.gz"),
        "top_mapping": os.path.join(out_path, "top_mapping.csv.gz"),
        "all_mapping": os.path.join(out_path, "all_mapping.csv.gz"),
        "top_genes": os.path.join(out_path, "top_genes.csv.gz"),
        "all_genes": os.path.join(out_path, "all_genes.csv.gz"),
    }
    for df_name, file_path in dataframes.items():
        vars()[df_name].to_csv(file_path, index=False, compression='gzip')

def main():
    """
    Main function to create the output directory and call data alignment functions.
    """
    # config = yaml.safe_load(open('src/param_config.yaml'))
    in_path = COMPRESSED_PATH
    out_path = ALIGNED_PATH
    if not out_path.exists():
        out_path.mkdir(parents=True)
    align_data(in_path=in_path, out_path=out_path)

if __name__ == '__main__':
    main()
