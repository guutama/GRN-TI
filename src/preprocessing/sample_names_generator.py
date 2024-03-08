import os,sys
import pandas as pd
directories_to_append = [
    os.path.abspath(os.path.join(os.path.dirname(__file__), '.')),
    os.path.abspath(os.path.join(os.path.dirname(__file__), '..')),
    # Add more directories as needed
]
for directory in directories_to_append:
    sys.path.append(directory) 
from path_config import (
    RAW_PATH,
    RAW_CSV_GZ_PATH
)

def get_samples(dir_path,out_path):
    file = os.path.join(dir_path, "E-GEUV-2.sdrf.txt")
    sample_info_data = pd.read_csv(file, sep="\t")
    populations = ["GBR", "FIN", "CEU", "TSI"]  # "YRI",
    unwanted_col = ['HG00107', 'HG00237', 'NA07000']
    # Filter sample_info_data based on populations
    filtered_data = sample_info_data[sample_info_data['Characteristics[population]'].isin(populations)]
    # Remove unwanted columns from Source Name
    filtered_data = filtered_data[~filtered_data['Source Name'].isin(unwanted_col)]

    # Extract unique values for cols_intersect
    cols_intersect = filtered_data['Source Name'].unique().tolist()
    
    # Convert cols_intersect to DataFrame
    cols_intersect_df = pd.DataFrame(cols_intersect, columns=['Source Name'])
    
    # Save cols_intersect_df to a compressed CSV file
    output_file_path = os.path.join(out_path, "sample_names.csv.gz")
    cols_intersect_df.to_csv(output_file_path, sep="\t", index=False, compression="gzip")
    return cols_intersect 


if __name__ =='__main__':
    in_path = RAW_PATH
    out_path = RAW_CSV_GZ_PATH
    get_samples(in_path,out_path)