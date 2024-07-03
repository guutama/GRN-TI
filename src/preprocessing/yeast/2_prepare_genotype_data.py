#!/usr/bin/env python
# coding: utf-8


# This script is adapted from the original found at:
# https://github.com/michoel-lab/FindrCausalNetworkInferenceOnYeast
# It has been modified with permission to suit our specific needs.


#  Preprocessing the genotype data matrices with genes
#   in the same order as in the eqtl list.
#
#   NOTE: need to prepare eQTL list first !
#
# 0. imports first
import pandas as pan
import numpy as np

import sys 
import os 
import yaml
directories_to_append = [
    os.path.abspath(os.path.join(os.path.dirname(__file__), '..')),
    os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')),
    # Add more directories as needed
]
for directory in directories_to_append:
    sys.path.append(directory)
import pandas as pd





from path_config import (
    YEAST_PATH,
    COMPRESSED_PATH
)

def process_genotype(in_path,out_path):
    
    #  1. Load data:
    path_geno = os.path.join(in_path, "SI_Data_03_genotypes.txt.zip")

    genotypes = pan.read_csv( path_geno,sep='\t')

    genotypes  = genotypes.reset_index().rename(columns={'index': 'sample'})

   
    sample_col = genotypes[['sample']] 

    gt_binary = (genotypes.drop('sample', axis=1) > 0).astype(int)
    gt_binary['sample'] = genotypes['sample']
    gt_binary= pd.concat([sample_col, gt_binary], axis=1)
   

    out_file_geno_all = os.path.join(out_path, 'genotypes_binary_all.csv')
    gt_binary.to_csv(out_file_geno_all,index=False)


        #   load ordering from eqtls ...
    path_eqtls = COMPRESSED_PATH

    ser = pan.read_csv( os.path.join(path_eqtls, 'strongest_eqtls_r.csv'))
    sorted_ser = ser.sort_values('pmarker')

 
    eqtl_labels = sorted_ser.pmarker
    

    sel_gt = genotypes[eqtl_labels.values].copy()
    
 
    print( '#   Write output:' )
    #   transform format of genotypes to binary:
    sel_gt = (sel_gt > 0).astype(int)
    sel_gt = pd.concat([sample_col, sel_gt], axis=1)
    print(sel_gt)

    #   3. Write output
    #       save genotypes to file
    out_file = os.path.join(out_path, 'genotypes_binary_strongest_eqtl.csv')
    sel_gt.to_csv(out_file,index=False)
    #print(sel_gt)

    print( '# Done.' )
    # EOF



def main():
    out_path = COMPRESSED_PATH
    in_path = YEAST_PATH
    process_genotype(in_path=in_path,
                     out_path=out_path) 



    

if __name__ == '__main__':
    main()


# # EOF

