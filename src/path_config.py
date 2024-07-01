from pathlib import Path

# --------------------------------------------------------------------------------
# CONFIGURATION PATHS - MODIFY THESE PATHS AS NEEDED
# --------------------------------------------------------------------------------

# Path to 'libfindr.so' - Update this path according to your environment
FINDR_PATH = "/cluster/projects/nn1015k/findr/libfindr.so"  # TODO: Needs to be updated

# Path to GRN-TI - Update this path according to your environment
MY_PATH = Path("/cluster/projects/nn1015k/GRN-TI")  # TODO: Needs to be modified
DATA_PATH = MY_PATH / "data"
RAW_PATH = DATA_PATH / "raw"
PREPROCESSED_PATH = DATA_PATH / "preprocessed"
RECODED_PATH = PREPROCESSED_PATH / "recoded"
MERGED_PATH = PREPROCESSED_PATH / "merged"
ALIGNED_PATH = PREPROCESSED_PATH / "aligned"
SPLIT_PATH = PREPROCESSED_PATH / "split"
PAIRWISE_PROBABILITY_PATH =  DATA_PATH / "pairwise_probability"
NETWORKS_INFO_PATH = DATA_PATH / "network_info"
NETWORKS_PATH = DATA_PATH / "networks" 
FEATURIZED_PATH = DATA_PATH / "featurized"

MODELS_PATH = MY_PATH / 'models'
METRICES_PATH = MY_PATH / 'metrices'
















