



import pandas as pd 
import numpy as np 
from scipy.stats import kruskal
import sys
import os
# Append directories to sys.path for module imports
directories_to_append = [
    os.path.abspath(os.path.join(os.path.dirname(__file__), '..')),
    # Add more directories as needed
]
for directory in directories_to_append:
    sys.path.append(directory)


from path_config import (
    DREAM_PATH,
    ALIGNED_PATH
) 


def get_pvalues(gene_expression_df, eqtl_data_df):
    """
    Computes the p-values for each gene using the Kruskal-Wallis H-test between gene expression data and eQTL data.

    Parameters:
    gene_expression_df (pd.DataFrame): A DataFrame containing gene expression data. Each row represents a gene, and each column represents a sample.
    eqtl_data_df (pd.DataFrame): A DataFrame containing eQTL data corresponding to the gene expression data. Each row represents a gene, and each column represents a sample.

    Returns:
    list: A list of p-values resulting from the Kruskal-Wallis H-test for each gene.
    """
    p_values = [] 
    for i, (gene_idx, gene_row) in enumerate(gene_expression_df.iterrows()):
        eqtl_row = eqtl_data_df.iloc[i]
        _, p_value = kruskal(gene_row, eqtl_row)
        p_values.append(p_value)
    return p_values 


def find_cis_eqtls(gene_expression_file, genotype_file):
    """
    Identify cis-eQTLs based on the strength of association between gene expression and genetic variants.

    Parameters:
    gene_expression_file (str): The path to the file containing gene expression data (in .tsv format).
    genotype_file (str): The path to the file containing genotype data (in .tsv format).

    Returns:
    tuple: A tuple containing the following elements:
        - cis_eqTL (pd.DataFrame): A DataFrame containing the cis-eQTL data.
        - A_genes (pd.DataFrame): A DataFrame containing the subset of gene expression data for cis-eQTL analysis.
        - all_genes (pd.DataFrame): A DataFrame containing all the gene expression data.
    """
    gene_expression_data = pd.read_csv(gene_expression_file, sep="\t", index_col=[-1])
    gene_expression_data = gene_expression_data.reset_index(drop=True)
    genotype = pd.read_csv(genotype_file, sep="\t", index_col=[-1])
    genotype = genotype.reset_index(drop=True)

    gene_expression_data = gene_expression_data.T
    genotype = genotype.T

    p = get_pvalues(gene_expression_data, genotype)
    sorted_ind = np.argsort(p) + 1
    number_value_to_get = int(len(sorted_ind) * 0.25)
    ind25 = sorted_ind[:number_value_to_get] 
    gene_id = np.array(['G' + str(x) for x in ind25])

    gene_expression_data["GENE_ID"] = gene_expression_data.index
    genotype["SNP_ID"] = genotype.index

    A_genes = gene_expression_data[gene_expression_data["GENE_ID"].isin(gene_id)]
    B_genes = gene_expression_data[~gene_expression_data["GENE_ID"].isin(gene_id)]
    all_genes = pd.concat([A_genes, B_genes], ignore_index=True).drop_duplicates()

    genotype = genotype.reindex(index=all_genes.GENE_ID)
    cis_eqTL = genotype.copy().iloc[:len(A_genes)]
    
    return cis_eqTL, A_genes, all_genes 



def process_dream_data(in_path: str, out_path: str) -> None:
    """
    Process DREAM data to identify cis-eQTLs and save the results to compressed CSV files.

    Parameters:
    in_path (str): The input directory path containing the DREAM data files.
    out_path (str): The output directory path where the processed data files will be saved.
    """
    dream_network = 5
  
    file_paths_exp = os.path.join(in_path, f'DREAM5_SysGenA999_Network{dream_network}_Expression.tsv')
    file_paths_geno = os.path.join(in_path, f'DREAM5_SysGenA999_Network{dream_network}_Genotype.tsv')

    geno, a_genes, all_genes = find_cis_eqtls(gene_expression_file=file_paths_exp,
                                              genotype_file=file_paths_geno)
   
    all_names = pd.Series(all_genes.index)
    a_names = pd.Series(a_genes.index)
    sample = pd.Series(all_genes.columns)
    
    dataframes = {
        "geno": geno,
        "a_genes": a_genes,
        "all_genes": all_genes,
        "all_names": all_names,
        "a_names": a_names,
        "sample": sample,
    }

    file_paths = {
        "geno": os.path.join(out_path, "geno_best_dream.csv.gz"),
        "a_genes": os.path.join(out_path, "a_genes_dream.csv.gz"),
        "all_genes": os.path.join(out_path, "all_genes_dream.csv.gz"),
        "all_names": os.path.join(out_path, "all_names_dream.csv.gz"),
        "a_names": os.path.join(out_path, "a_names_dream.csv.gz"),
        "sample": os.path.join(out_path, "sample_dream.csv.gz"),
    }

    # Save each dataframe to its corresponding compressed CSV file
    for df_name, file_path in file_paths.items():
        dataframes[df_name].to_csv(file_path, index=False, compression='gzip')




in_path = DREAM_PATH
out_path = ALIGNED_PATH
if not out_path.exists():
        out_path.mkdir(parents=True)

process_dream_data(in_path, out_path)