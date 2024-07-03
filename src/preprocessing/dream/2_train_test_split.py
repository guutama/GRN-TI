import os,sys
import numpy as np
import pandas as pd
import yaml
from utils import rank_inverse_normal_transformation
from path_config import ALIGNED_PATH, SPLIT_PATH 


# Append directories to sys.path for module imports
directories_to_append = [
    os.path.abspath(os.path.join(os.path.dirname(__file__), '..')),
    os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')),
    # Add more directories as needed
]
for directory in directories_to_append:
    sys.path.append(directory)
for directory in directories_to_append:
    sys.path.append(directory)

# Load preprocessing parameters from YAML configuration file
params = yaml.safe_load(open("src/param_config.yaml"))["preprocessing"]
params_split = params['split']

# Define input and output paths
in_path = ALIGNED_PATH
out_path = SPLIT_PATH

# Define the test size for data splitting
test_size = float(params_split['test_size'])

# Create output directory if it doesn't exist
if not out_path.exists():
    out_path.mkdir(parents=True)

# Define file paths for input data
geno_best_path = os.path.join(in_path, "geno_best_dream.csv.gz")
a_genes_path = os.path.join(in_path, "a_genes_dream.csv.gz")
all_genes_path = os.path.join(in_path, "all_genes_dream.csv.gz")

# Load gene names from files
a_names = pd.read_csv(os.path.join(in_path, "a_names_dream.csv.gz")).iloc[:,0].values
all_names = pd.read_csv(os.path.join(in_path, "all_names_dream.csv.gz")).iloc[:,0].values

# Load data into DataFrames
geno_best = pd.read_csv(geno_best_path)
a_genes = pd.read_csv(a_genes_path)
all_genes = pd.read_csv(all_genes_path)

# Extract sample columns and gene expression data
samples = a_genes.drop("GENE_ID", axis=1).columns
exp_a = a_genes[samples]
exp_all = all_genes[samples]
eqtl_best = geno_best[samples]

# Apply rank inverse normal transformation to gene expression data
exp_a = rank_inverse_normal_transformation(exp_a)
exp_all = rank_inverse_normal_transformation(exp_all)

# Convert transformed data back to DataFrames with appropriate column names
exp_a = pd.DataFrame(exp_a, columns=samples)
exp_all = pd.DataFrame(exp_all, columns=samples)

# Add gene IDs back to the transformed data
exp_a["GENE_ID"] = a_names
exp_all["GENE_ID"] = all_names

# Shuffle samples
exp_a = exp_a.sample(axis=1, frac=1).reset_index(drop=True)
sample_shuffled = exp_a.drop("GENE_ID", axis=1).columns.to_list()

# Reorder other dataframes according to the shuffled samples
exp_all = exp_all[sample_shuffled]
eqtl_best = eqtl_best[sample_shuffled]

# Split data into training and testing sets
split_index = int(exp_a.shape[1] * (1 - test_size))
train_sample = sample_shuffled[:split_index]
test_sample = sample_shuffled[split_index:]

exp_a_train = exp_a.iloc[:, :split_index].copy()
exp_a_test = exp_a.iloc[:, split_index:].copy()
exp_all_train = exp_all.iloc[:, :split_index].copy()
exp_all_test = exp_all.iloc[:, split_index:].copy()
eqtl_best_train = eqtl_best.iloc[:, :split_index].copy()
eqtl_best_test = eqtl_best.iloc[:, split_index:].copy()

# Add gene IDs back to the split data
exp_a_train.loc[:, "GENE_ID"] = a_names
exp_a_test.loc[:, "GENE_ID"] = a_names
exp_all_train.loc[:, "GENE_ID"] = all_names
exp_all_test.loc[:, "GENE_ID"] = all_names
eqtl_best_train.loc[:, "SNP_ID"] = a_names
eqtl_best_test.loc[:, "SNP_ID"] = a_names

# Convert training and testing sample lists to Series
train_sample = pd.Series(train_sample)
test_sample = pd.Series(test_sample)

# Define paths for saving the split data
dataframes = {
    "exp_a_train": os.path.join(out_path, "top_exp_train_dream.csv.gz"),
    "exp_a_test": os.path.join(out_path, "top_exp_test_dream.csv.gz"),
    "exp_all_train": os.path.join(out_path, "all_exp_train_dream.csv.gz"),
    "exp_all_test": os.path.join(out_path, "all_exp_test_dream.csv.gz"),
    "eqtl_best_train": os.path.join(out_path, "top_eqtl_train_dream.csv.gz"),
    "eqtl_best_test": os.path.join(out_path, "top_eqtl_test_dream.csv.gz"),
    "train_sample": os.path.join(out_path, "train_sample_dream.csv.gz"),
    "test_sample": os.path.join(out_path, "test_sample_dream.csv.gz"),
}

# Save each DataFrame to a compressed CSV file
for df_name, file_path in dataframes.items():
    vars()[df_name].to_csv(file_path, index=False, compression='gzip')
