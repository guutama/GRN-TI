import sys 
import os 
import yaml 
import pandas as pan
from sklearn.model_selection import train_test_split
params = yaml.safe_load(open("src/param_config.yaml"))["preprocessing"]
params_split = params['split']


# Append necessary directories to the system path
directories_to_append = [
    os.path.abspath(os.path.join(os.path.dirname(__file__), '..')),
    # Add more directories as needed
]
for directory in directories_to_append:
    sys.path.append(directory)

from utils import get_numpy_array
from path_config import (
    COMPRESSED_PATH,
    SPLIT_PATH
) 

# Initialize variables for data dimensions (can be set to specific values if needed)
n_eqtls = None
n_genes = None
n_segregants = None 
test_size = float(params_split['test_size'])
# Define input and output paths
in_path = COMPRESSED_PATH
out_path = SPLIT_PATH

# Load genes-eqtls table
eqtl_path = os.path.join(in_path, 'strongest_eqtls_r_3columns.csv')
gene_eqtl_table = pan.read_csv(eqtl_path, nrows=n_genes, delimiter=',')

# Sort and remove duplicates based on the 'pmarker' column
gene_eqtl_table = gene_eqtl_table.sort_values(by='pmarker')
gene_eqtl_table = gene_eqtl_table.reset_index(drop=True)
gene_eqtl_table = gene_eqtl_table.drop_duplicates("pmarker") 

# Load genotype data for strongest eQTLs
geno_path = os.path.join(in_path, "genotypes_binary_strongest_eqtl.csv")
genotypes = pan.read_csv(geno_path, nrows=n_segregants, delimiter=',')  

# Load genotype data for all eQTLs
geno_all_path = os.path.join(in_path, "genotypes_binary_all.csv")
genotypes_all = pan.read_csv(geno_all_path, nrows=n_segregants, delimiter=',')

# Extract relevant columns from genotypes
eqtl_cols = gene_eqtl_table['pmarker']
genotypes = genotypes[eqtl_cols] 

# Load gene expression data
exp_path = os.path.join(in_path, "expression_statsmodels_linreg_residuals_01.csv")
yeast_exp = pan.read_table(exp_path, nrows=n_segregants, delimiter=',')

# Extract sample column and combine with genotypes
sample = yeast_exp['sample']
genotypes = pan.concat([sample, genotypes], axis=1)

# Identify gene columns that are present in the gene expression data
gene_cols = gene_eqtl_table['gene']
count_occuring_genes = sum([1 for col in gene_cols if col in yeast_exp.columns])

# Identify missing and extra genes
missing_genes = [col for col in gene_cols if col not in yeast_exp.columns]
extra_genes = [col for col in yeast_exp.columns if col not in list(gene_cols)]

# Sort columns of expression values
yeast_exp_with_eqtls = yeast_exp[gene_cols]
yeast_exp_with_eqtls = pan.concat([sample, yeast_exp_with_eqtls], axis=1)

yeast_exp_without_eqtls = yeast_exp[extra_genes]

# Combine the expression data with and without eQTLs
yeast_exp = pan.concat([yeast_exp_with_eqtls, yeast_exp_without_eqtls.drop('sample', axis=1)], sort=False, axis=1)
yeast_exp = yeast_exp.T.drop_duplicates(keep='first').T

# Split data into training and testing sets based on the 'sample' column
train_samples, test_samples = train_test_split(sample, test_size=test_size, random_state=42)

# Subset DataFrames based on training and testing samples
train_exp = yeast_exp_with_eqtls[yeast_exp_with_eqtls['sample'].isin(train_samples)]
test_exp = yeast_exp_with_eqtls[yeast_exp_with_eqtls['sample'].isin(test_samples)] 

train_exp_all = yeast_exp[yeast_exp['sample'].isin(train_samples)]
test_exp_all = yeast_exp[yeast_exp['sample'].isin(test_samples)]

train_geno = genotypes[genotypes['sample'].isin(train_samples)]
test_geno = genotypes[genotypes['sample'].isin(test_samples)]

train_geno_all = genotypes_all[genotypes_all['sample'].isin(train_samples)]
test_geno_all = genotypes_all[genotypes_all['sample'].isin(test_samples)]

# Define paths for saving the split data
dataframes = {
    "train_exp": os.path.join(out_path, "top_exp_train_yeast.csv.gz"),
    "test_exp": os.path.join(out_path, "top_exp_test_yeast.csv.gz"),
    "train_exp_all": os.path.join(out_path, "all_exp_train_yeast.csv.gz"),
    "test_exp_all": os.path.join(out_path, "all_exp_test_yeast.csv.gz"),
    "train_geno": os.path.join(out_path, "top_eqtl_train_yeast.csv.gz"),
    "test_geno": os.path.join(out_path, "top_eqtl_test_yeast.csv.gz"),
    "train_geno_all": os.path.join(out_path, "all_eqtl_train_yeast.csv.gz"),
    "test_geno_all": os.path.join(out_path, "all_eqtl_test_yeast.csv.gz"),
    "train_samples": os.path.join(out_path, "train_sample_yeast.csv.gz"),
    "test_samples": os.path.join(out_path, "test_sample_yeast.csv.gz"),
}

# Save each DataFrame to a compressed CSV file
for df_name, file_path in dataframes.items():
    vars()[df_name].to_csv(file_path, index=False, compression='gzip')
