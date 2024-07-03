import sys
import os
from pathlib import Path
import pandas as pd
import yaml
from utils import get_samples, rank_inverse_normal_transformation, transform_to_log
from path_config import ALIGNED_PATH, SPLIT_PATH

# Append directories to sys.path for module imports
directories_to_append = [
    os.path.abspath(os.path.join(os.path.dirname(__file__), '..')),
    # Add more directories as needed
]
for directory in directories_to_append:
    sys.path.append(directory)

def split(in_path, out_path, test_size):
    """
    Split the dataset into training and testing sets.

    Parameters:
    in_path (str): Directory path where the input files are located.
    out_path (str): Directory path where the split files will be saved.
    test_size (float): Proportion of the dataset to include in the test split.

    Returns:
    None
    """
    # Define file paths
    geno_best_path = os.path.join(in_path, "top_snp.csv.gz")
    geno_all_path = os.path.join(in_path, "all_snp.csv.gz")
    map_best_path = os.path.join(in_path, "top_mapping.csv.gz")
    a_genes_path = os.path.join(in_path, "top_genes.csv.gz")
    all_genes_path = os.path.join(in_path, "all_genes.csv.gz")
    samples_path = os.path.join(in_path, "sample_names.csv.gz")

    # Load each DataFrame
    geno_best = pd.read_csv(geno_best_path)
    geno_all = pd.read_csv(geno_all_path)
    map_best = pd.read_csv(map_best_path)
    a_genes = pd.read_csv(a_genes_path)
    all_genes = pd.read_csv(all_genes_path)
    samples = pd.read_csv(samples_path).iloc[:, 0].values

    # Assert that columns 'TargetID' in a_genes and 'GENE_ID' in map_best are the same
    assert a_genes['GENE_ID'].equals(map_best['GENE_ID']), "a_genes['GENE_ID'] and map_best['GENE_ID'] are not the same"

    # Assert that columns 'SNP_ID' in map_best and 'SNP' in geno_best are the same
    assert map_best['SNP_ID'].equals(geno_best["SNP_ID"]), "map_best['SNP_ID'] and geno_best['SNP'] are not the same"

    print("All specified columns are exactly the same.")
    print(a_genes)

    exp_a = a_genes[samples]
    exp_all = all_genes[samples]

    print(exp_a)
    exp_a = transform_to_log(exp_a)
    exp_all = transform_to_log(exp_all)
    print(exp_a)

    exp_sub_a = exp_a.to_numpy(dtype='float64')
    exp_sub_all = exp_all.to_numpy(dtype='float64')

    exp_a = pd.DataFrame(exp_a, columns=samples)
    exp_all = pd.DataFrame(exp_all, columns=samples)

    eqtl_best = geno_best[samples]
    eqtl_all = geno_all[samples]

    exp_a = exp_a.sample(axis=1, frac=1).reset_index(drop=True)  # shuffling
    sample_shuffled = exp_a.columns.to_list()

    exp_all = exp_all[sample_shuffled]
    eqtl_all = eqtl_all[sample_shuffled]
    eqtl_best = eqtl_best[sample_shuffled]

    split_index = int(exp_a.shape[1] * (1 - test_size))
    train_sample = sample_shuffled[:split_index]
    test_sample = sample_shuffled[split_index:]

    top_exp_train = exp_a.iloc[:, :split_index].copy()
    top_exp_test = exp_a.iloc[:, split_index:].copy()

    all_exp_train = exp_all.iloc[:, :split_index].copy()
    all_exp_test = exp_all.iloc[:, split_index:].copy()

    top_eqtl_train = eqtl_best.iloc[:, :split_index].copy()
    top_eqtl_test = eqtl_best.iloc[:, split_index:].copy()

    all_eqtl_train = eqtl_all.iloc[:, :split_index].copy()
    all_eqtl_test = eqtl_all.iloc[:, split_index:].copy()

    top_exp_train.loc[:, "GENE_ID"] = a_genes.loc[:, "GENE_ID"].copy()
    top_exp_test.loc[:, "GENE_ID"] = a_genes.loc[:, "GENE_ID"].copy()

    all_exp_train.loc[:, "GENE_ID"] = all_genes.loc[:, "GENE_ID"].copy()
    all_exp_test.loc[:, "GENE_ID"] = all_genes.loc[:, "GENE_ID"].copy()

    top_eqtl_train.loc[:, "SNP_ID"] = geno_best.loc[:, "SNP_ID"].copy()
    top_eqtl_test.loc[:, "SNP_ID"] = geno_best.loc[:, "SNP_ID"].copy()

    all_eqtl_train.loc[:, "SNP_ID"] = geno_all.loc[:, "SNP_ID"].copy()
    all_eqtl_test.loc[:, "SNP_ID"] = geno_all.loc[:, "SNP_ID"].copy()

    train_sample = pd.Series(train_sample)
    test_sample = pd.Series(test_sample)

    dataframes = {
        "top_exp_train": os.path.join(out_path, "top_exp_train.csv.gz"),
        "top_exp_test": os.path.join(out_path, "top_exp_test.csv.gz"),
        "all_exp_train": os.path.join(out_path, "all_exp_train.csv.gz"),
        "all_exp_test": os.path.join(out_path, "all_exp_test.csv.gz"),
        "top_eqtl_train": os.path.join(out_path, "top_eqtl_train.csv.gz"),
        "top_eqtl_test": os.path.join(out_path, "top_eqtl_test.csv.gz"),
        "all_eqtl_train": os.path.join(out_path, "all_eqtl_train.csv.gz"),
        "all_eqtl_test": os.path.join(out_path, "all_eqtl_test.csv.gz"),
        "train_sample": os.path.join(out_path, "train_sample.csv.gz"),
        "test_sample": os.path.join(out_path, "test_sample.csv.gz"),
    }

    for df_name, file_path in dataframes.items():
        vars()[df_name].to_csv(file_path, index=False, compression='gzip')

if __name__ == '__main__':
    """
    Main function to load parameters, create output directory, and call the split function.
    """
    params = yaml.safe_load(open("src/config_params.yaml"))["preprocessing"]
    params_split = params['split']

    in_path = ALIGNED_PATH
    out_path = SPLIT_PATH

    test_size = float(params_split['test_size'])
    if not out_path.exists():
        out_path.mkdir(parents=True)
    split(in_path=in_path, out_path=out_path, test_size=test_size)
