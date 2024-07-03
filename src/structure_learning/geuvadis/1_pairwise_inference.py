import sys 
import os 
import yaml
import pandas as pd
import numpy as np

import findr
from utils import get_numpy_array
from path_config import (
    FINDR_PATH,
    ALIGNED_PATH,
    PAIRWISE_INFERENCE_PATH,
)

# Append necessary directories to the system path
directories_to_append = [
    os.path.abspath(os.path.join(os.path.dirname(__file__), '..')),
    os.path.abspath(os.path.join(os.path.dirname(__file__), '../preprocessing')),
    # Add more directories as needed
]
for directory in directories_to_append:
    sys.path.append(directory)

def calculate_p_values(expression_A, expression_ALL, genotype, a_genes_only, method, n=None):
    """
    Calculate p-values using the Findr library.

    Parameters:
    expression_A (ndarray): Gene expression data for A genes.
    expression_ALL (ndarray): Gene expression data for all genes.
    genotype (ndarray): Genotype data.
    a_genes_only (bool): Flag to indicate if only A genes should be considered.
    method (object): Findr method object for calculating p-values.
    n (int): Number of genes to consider. Defaults to None.

    Returns:
    dict: Dictionary containing calculated posterior-values for different tests.
    """
    # Calculate p0
    p0_results = method.pij_rank(dt=expression_A, dt2=expression_ALL, nodiag=True)
    p0 = p0_results["p"][:, :n] if a_genes_only else p0_results["p"]

    # Calculate p_other_results
    p_other_results = method.pijs_gassist(dg=genotype, dt=expression_A, dt2=expression_ALL, nodiag=True)
    p2, p3, p4, p5 = [p_other_results[key] if not a_genes_only else p_other_results[key][:, :n] for key in ["p2", "p3", "p4", "p5"]]

    # Calculate combined probabilities
    p2p3 = p2 * p3
    p2p5 = p2 * p5
    p = 0.5 * (p2p5 + p4)

    # Return test results
    return {
        "p0": p0,
        "p2": p2,
        "p3": p3,
        "p4": p4,
        "p5": p5,
        "p2p3": p2p3,
        "p2p5": p2p5,
        "p": p,
        "random": p
    }

def pairwise_inference(in_path, out_path, a_genes_only=True):
    """
    Perform pairwise inference to calculate edge posteriors.

    Parameters:
    in_path (str): Input directory path containing gene and genotype data.
    out_path (str): Output directory path to save the results.
    a_genes_only (bool): Flag to indicate if only A genes should be considered. Defaults to True.
    """
    # Create output directory if it doesn't exist
    if not out_path.exists():
        out_path.mkdir(parents=True)

    # Initialize Findr library
    findr_lib = findr.lib(path=FINDR_PATH, loglv=6, rs=0, nth=0)

    # Load configuration parameters for structure learning
    params = yaml.safe_load(open('src/param_config.yaml'))['structure_learning']
    params_pairwise = params['pairwise_inference']
    network_types = ["p0", "p2", "p2p3", "p2p5", "p"]

    # Load each DataFrame individually
    exp_a_train = pd.read_csv(os.path.join(in_path, "top_exp_train.csv.gz"))
    exp_all_train = pd.read_csv(os.path.join(in_path, "all_exp_train.csv.gz"))
    eqtl_best_train = pd.read_csv(os.path.join(in_path, "top_eqtl_train.csv.gz"))
    sample_ids = pd.read_csv(os.path.join(in_path, "train_sample.csv.gz")).iloc[:, 0].values

    # Convert DataFrames to numpy arrays
    dt = get_numpy_array(dataframe=exp_a_train, colums_to_keep=sample_ids)
    dt2 = get_numpy_array(dataframe=exp_all_train, colums_to_keep=sample_ids)
    dg = get_numpy_array(dataframe=eqtl_best_train, colums_to_keep=sample_ids)
    n = dt.shape[0]

    # Calculate posterior probabilities
    posteriors = calculate_p_values(expression_A=dt,
                                    expression_ALL=dt2,
                                    genotype=dg,
                                    n=n,
                                    a_genes_only=a_genes_only,
                                    method=findr_lib)

    # Save edge posteriors to files
    for network_type in network_types:
        print(f'Using network {network_type}', flush=True)
        if network_type in posteriors.keys():
            edge_posteriors = posteriors[network_type]
        edge_posteriors = pd.DataFrame(edge_posteriors, columns=exp_a_train["GENE_ID"].values)
        out_file = os.path.join(out_path, f"edge_posteriors_{network_type}.csv.gz")
        edge_posteriors.to_csv(out_file, compression='gzip', index=False)
        print(edge_posteriors)

if __name__ == '__main__':
    in_path = ALIGNED_PATH
    out_path = PAIRWISE_INFERENCE_PATH
    if not out_path.exists():
        out_path.mkdir(parents=True)

    pairwise_inference(in_path=in_path,
                       out_path=out_path)