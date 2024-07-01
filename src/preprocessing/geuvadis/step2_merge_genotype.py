
import os
import sys
import pandas as pd

directories_to_append = [
    os.path.abspath(os.path.join(os.path.dirname(__file__), '..')),
    os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')),
    # Add more directories as needed
]
for directory in directories_to_append:
    sys.path.append(directory)

from path_config import (
    RECODED_PATH,
    MERGED_PATH
)
import yaml
def merge_genotypes(dir_path, out_path, file_type):
    """
    Merges genotype files of a specified type from a directory into a single CSV file.

    Parameters:
    dir_path (str): Directory path where the genotype files are located.
    out_path (str): Directory path where the merged output file will be saved.
    file_type (str): There are two types of files ("best, all"). best file contains only the most significant eQTLs.

    Returns:
    None: Outputs a merged CSV file in the specified output directory.
    """
    # Create a list to hold dataframes
    df_list = []
    
    # Construct file glob pattern
    file_pattern = f"{file_type}.traw"
    
    # Iterate over files in dir_path that match the file_pattern
    for file in os.listdir(dir_path):
        if file.endswith(file_pattern):
            file_path = os.path.join(dir_path, file)
            df = pd.read_csv(file_path)
            df_list.append(df)
    
    # Merge all dataframes
    merged_df = pd.concat(df_list, ignore_index=True)
    
    # Output file
    output_file = os.path.join(out_path, f"merged_{file_type}.csv")
    merged_df.to_csv(output_file, index=False)
    print(f"Merged file {output_file} created.")

if __name__ == '__main__':
    dir_path = RECODED_PATH  # Directory path where the files are located
    out_path =  MERGED_PATH # Output directory path
    if not out_path.exists():
        out_path.mkdir(parents=True)
   
    config = yaml.safe_load(open('src/param_config.yaml'))
    preprocessing_config = config['preprocessing']
    file_type = preprocessing_config['geuvadis']['genotype_file']  # The file type we want to merge
    merge_genotypes(dir_path, out_path, file_type)