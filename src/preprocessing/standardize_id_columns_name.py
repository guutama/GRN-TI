import pandas as pd
import sys
import os
import yaml
directories_to_append = [
    os.path.abspath(os.path.join(os.path.dirname(__file__), '.')),
    os.path.abspath(os.path.join(os.path.dirname(__file__), '..')),
    # Add more directories as needed
]
for directory in directories_to_append:
    sys.path.append(directory)

from path_config import (
    RAW_CSV_TXT_PATH,
    RAW_CSV_GZ_PATH
)

def read_csv_file(file_path, sep=None):
    """Read a CSV file and return the DataFrame."""
    return pd.read_csv(file_path, sep=sep)

def save_dataframe(df, output_path, filename, compression='gzip'):
    """Save DataFrame to a compressed CSV file."""
    output_file = os.path.join(output_path, filename)
    df.to_csv(output_file, index=False, compression=compression)
    print(f"Data saved to {output_file}.")

def check_and_rename_columns(df, old_id_name, new_id_name):
    """Check if old_id_name exists and rename it to new_id_name if they are not the same."""
    if old_id_name in df.columns:
        if old_id_name != new_id_name:
            df.rename(columns={old_id_name: new_id_name}, inplace=True)
        return True
    else:
        return False

def process_file(dir_path, out_path, filename, check_column, rename_column=None, new_name=None, file_type='csv'):
    """Process the file by checking, renaming, and saving the DataFrame."""
    file_path = os.path.join(dir_path, filename)
    sep = '\t' if file_type == 'tsv' else ','
    df = read_csv_file(file_path, sep=sep)
    print(df)
    # For mapping, check if both new_gene_id and new_snp_id exist
    if isinstance(check_column, list):
        if all(col in df.columns for col in check_column):
            save_dataframe(df, out_path, filename, compression='gzip')
        else:
            raise ValueError(f"Not all required columns ({check_column}) are present in the mapping file.")
    else:
        # For expression and genotype, check and rename old_id_name to new_id_name
        if check_and_rename_columns(df, check_column, rename_column if rename_column else check_column):
            save_dataframe(df, out_path, filename.replace('.csv', '.csv.gz'), compression='gzip')
        else:
            raise ValueError(f"Required column ({check_column}) is not present in the file.")






if __name__ == "__main__":
    # Example directory paths
    dir_path = RAW_CSV_TXT_PATH
    out_path = RAW_CSV_GZ_PATH
    config = yaml.safe_load(open('src/config.yaml'))
    preprocessing_params = config['preprocessing']

    # Updated variable names for clarity
    expression_gene_column = preprocessing_params['expression_gene_id_col_name']
    mapping_gene_column = preprocessing_params['mapping_gene_id_col_name']
    genotype_snp_column = preprocessing_params['genotype_snp_id_col_name']
    mapping_snp_column = preprocessing_params['mapping_snp_id_col_name']

    # Process each dataset with specific checks
    process_file(dir_path, out_path, 'mapping.csv', [mapping_gene_column, mapping_snp_column])
    process_file(dir_path, out_path, 'expression.tsv', expression_gene_column, mapping_gene_column,file_type='tsv')
    process_file(dir_path, out_path, 'genotype.csv', genotype_snp_column, mapping_snp_column)
    