from pathlib import Path

# --------------------------------------------------------------------------------
# CONFIGURATION PATHS - MODIFY THESE PATHS AS NEEDED
# --------------------------------------------------------------------------------

# Path to 'libfindr.so' - Update this path according to your environment
FINDR_PATH = "/cluster/projects/nn1015k/findr/libfindr.so"  # TODO: Needs to be updated

# Path to GRN-TI - Update this path according to your environment
MY_PATH = Path("/cluster/projects/nn1015k/GRN-TI")  # TODO: Needs to be updated








DATA_PATH = MY_PATH / "data"
RAW_CSV_TXT_PATH = DATA_PATH / "raw_csv_txt" 
RAW_CSV_GZ_PATH = DATA_PATH / "raw_csv_gz" 
ALIGNED_PATH = DATA_PATH / "aligned"
















