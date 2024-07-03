import os
import numpy as np
import pandas as pan
# import findr
# import roman 


import sys 
import os 
import yaml 
import findr
from utils import get_numpy_array


# Append necessary directories to the system path
directories_to_append = [
    os.path.abspath(os.path.join(os.path.dirname(__file__), '..')),
    # Add more directories as needed
]
for directory in directories_to_append:
    sys.path.append(directory) 


from path_config import (
    FINDR_PATH,
    SPLIT_PATH,
    PAIRWISE_INFERENCE_PATH

    
)


def calculate_p_values(expression_A, expression_ALL, genotype, a_genes_only, method, n=None):
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
        "random":p
    }



def main(in_path,out_path):
    """ Function doc
        Run findr P0: pij_rank, and causl tests: pij_gassist on yeast data and save results
    """
   

    # print('# 0.  Initialise findr.')
    l=findr.lib(path=FINDR_PATH,loglv = 6, rs=0,nth=0)
  

    yeast_exp_with_eqtls = pan.read_csv(os.path.join(in_path, "top_exp_train_yeast.csv.gz"))
    yeast_exp  = pan.read_csv(os.path.join(in_path, "all_exp_train_yeast.csv.gz"))
    genotypes= pan.read_csv(os.path.join(in_path, "top_eqtl_train_yeast.csv.gz"))
    # sample_ids = get_samples(dir_path=GEUV_INFO_PATH)
    sample_ids = pan.read_csv(os.path.join(in_path, "train_sample.csv.gz")).iloc[:,0].values

  

    yeast_exp_with_eqtls = pan.read_csv(os.path.join(in_path, "top_exp_train_yeast.csv.gz"))
    yeast_exp  = pan.read_csv(os.path.join(in_path, "all_exp_train_yeast.csv.gz"))
    genotypes= pan.read_csv(os.path.join(in_path, "top_eqtl_train_yeast.csv.gz"))
    # sample_ids = get_samples(dir_path=GEUV_INFO_PATH)
    sample_ids = pan.read_csv(os.path.join(in_path, "train_sample.csv.gz")).iloc[:,0].values



    gene_cols=  yeast_exp_with_eqtls.drop('sample', axis=1).columns
    print(gene_cols)
    array_yeast_exp =  yeast_exp.drop('sample', axis=1)
    array_genotypes = genotypes.drop('sample',axis=1)
    array_genes_with_eqtls =  yeast_exp_with_eqtls.drop('sample', axis=1)

   

    array_genotypes = np.transpose(array_genotypes).to_numpy(dtype=np.float64)
    array_yeast_exp = np.transpose(array_yeast_exp).to_numpy(dtype=np.float64)
    array_genes_with_eqtls = np.transpose(array_genes_with_eqtls).to_numpy(dtype=np.float64)
    n = array_genes_with_eqtls.shape[0]
     
    print(n)
    print(array_genes_with_eqtls.shape)
    print(array_genotypes.shape)
    print(array_yeast_exp.shape)

   
    posteriors = calculate_p_values(expression_A=array_genes_with_eqtls,
                                    expression_ALL=array_yeast_exp,
                                    genotype=array_genotypes,
                                    n=n,
                                    a_genes_only=True,
                                    method=l)
    network_types = ["p0","p2","p2p3","p2p5","p"] #"p0",

    for network_type in network_types:
        print(f'Using network {network_type}',flush=True)
        if network_type in posteriors.keys():
            edge_posteriors = posteriors[network_type]
        edge_posteriors = pan.DataFrame(edge_posteriors,columns=gene_cols)
        out_file = os.path.join(out_path, f"edge_posteriors_yeast_{network_type}.csv.gz")    
        edge_posteriors.to_csv(out_file,compression='gzip',index=False)
        print(edge_posteriors)


    print('# Done with Findr.')


# Run main:

if __name__ == '__main__':
    in_path = SPLIT_PATH
    out_path = PAIRWISE_INFERENCE_PATH
    main(in_path=in_path,out_path=out_path)