import sys
import os
import pandas as pd
import yaml

directories_to_append = [
    os.path.abspath(os.path.join(os.path.dirname(__file__), '..')),
    os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')),
    # Add more directories as needed
]
for directory in directories_to_append:
    sys.path.append(directory)

from path_config import (
    MAPPING_PATH,
    MAPPING_PROCESSED 
)

def process_mapping_file(in_path, file_name, out_path):
    """
    Processes geuvadis eqtl mapping file by  adding headers if necessary, and saving it as csv.

    Parameters:
    in_path (str): The input directory path where the mapping file is located.
    file_name (str): The name of the mapping file to be processed.
    out_path (str): The output directory path where the processed file will be saved.

    Returns:
    None: The function saves the processed file to the specified output path.
    """
    mapping_file = os.path.join(in_path, file_name)
    header = ['SNP_ID', 'ID', 'GENE_ID', 'PROBE_ID', 'CHR_SNP', 'CHR_GENE', 'SNPpos', 'TSSpos', 'distance', 'rvalue', 'pvalue', 'log10pvalue']

    if mapping_file:
        if file_name == "EUR373.gene.cis.FDR5.best.rs137.txt":
            file_ext = 'processed_best.rs137.csv'
            df = pd.read_csv(mapping_file, sep='\t', header=None)
            df.columns = header
        else:
            df = pd.read_csv(mapping_file, sep='\t')
            file_ext = 'processed_all.rs137.csv'
        out_file = os.path.join(out_path,file_ext)
        df.to_csv(out_file, index=False)

# Example usage:
if __name__ == '__main__':
    config = yaml.safe_load(open('src/param_config.yaml'))
    preprocessing_config = config['preprocessing']
    file_name = preprocessing_config['geuvadis']['mapping_file_name']
    in_path = MAPPING_PATH

    out_path = MAPPING_PROCESSED
    if not out_path.exists():
        out_path.mkdir(parents=True)
    process_mapping_file(in_path, file_name, out_path)
  
    