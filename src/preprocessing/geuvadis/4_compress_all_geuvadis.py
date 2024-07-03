import sys
import os
import pandas as pd
import numpy as np
from path_config import EXPRESSION_PATH, GENOTYPE_PATH, MAPPING_PATH, COMPRESSED_PATH

# Append directories to sys.path for module imports
directories_to_append = [
    os.path.abspath(os.path.join(os.path.dirname(__file__), '..')),
    # Add more directories as needed
]
for directory in directories_to_append:
    sys.path.append(directory)

def load_expression(in_path, out_path):
    """
    Load and process expression data, then save it as a compressed CSV file.

    Parameters:
    in_path (str): Directory path where the expression file is located.
    out_path (str): Directory path where the processed file will be saved.
    """
    exp_file = os.path.join(in_path, "GD462.GeneQuantRPKM.50FN.samplename.resk10.txt")
    try:
        df = pd.read_csv(exp_file, sep='\t')
        df.rename(columns={'TargetID': 'GENE_ID'}, inplace=True)
        file_path = os.path.join(out_path, "expression.csv.gz")
        df.to_csv(file_path, index=False, compression='gzip')
    except Exception as e:
        print(f"Error processing expression data: {e}")

def load_genotype(in_path, out_path, filename1='merged_all.csv'):
    """
    Load and process genotype data, then save it as a compressed CSV file.

    Parameters:
    in_path (str): Directory path where the genotype file is located.
    out_path (str): Directory path where the processed file will be saved.
    filename1 (str): Name of the genotype file to be processed. Default is 'merged_all.csv'.
    """
    file_path1 = os.path.join(in_path, filename1)
    try:
        df_all = pd.read_csv(file_path1)
        df_all.rename(columns={'SNP': 'SNP_ID'}, inplace=True)
        file_path = os.path.join(out_path, "genotypes_all.csv.gz")
        df_all.to_csv(file_path, index=False, compression='gzip')
    except Exception as e:
        print(f"Error processing genotype data: {e}")

def load_mapping(in_path, out_path, filename1='processed_all.rs137.csv'):
    """
    Load and process mapping data, then save it as a compressed CSV file.

    Parameters:
    in_path (str): Directory path where the mapping file is located.
    out_path (str): Directory path where the processed file will be saved.
    filename1 (str): Name of the mapping file to be processed. Default is 'processed_all.rs137.csv'.
    """
    file_path1 = os.path.join(in_path, filename1)
    try:
        df_all = pd.read_csv(file_path1)
        file_path = os.path.join(out_path, "mappings_all.csv.gz")
        df_all.to_csv(file_path, index=False, compression='gzip')
    except Exception as e:
        print(f"Error processing mapping data: {e}")

def main():
    """
    Main function to create the output directory and call data processing functions.
    """
    out_dir = COMPRESSED_PATH
    if not out_dir.exists():
        out_dir.mkdir(parents=True)

    load_expression(in_path=EXPRESSION_PATH, out_path=out_dir)
    load_genotype(in_path=GENOTYPE_PATH, out_path=out_dir)
    load_mapping(in_path=MAPPING_PATH, out_path=out_dir)

if __name__ == '__main__':
    main()
