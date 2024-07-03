

# This script is adapted from the original found at:
# https://github.com/michoel-lab/FindrCausalNetworkInferenceOnYeast
# It has been modified with permission to suit our specific needs.


import os
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





from path_config import(
    YEAST_PATH,
    COMPRESSED_PATH,
    
)


def align(in_path,out_path):
    #  1. Load expression data:
    exp_path = os.path.join(YEAST_PATH, "SI_Data_01_expressionValues.txt.zip")
    expr_values = pan.read_csv(exp_path ,sep='\t')
    
    #  2. Load ordering from eqtls ...
    ser = pan.read_csv(os.path.join(in_path, 'strongest_eqtls_r.csv'))
    #  Sort them by 'alphabetical' order of eqtls
    sorted_ser = ser.sort_values('pmarker')
    #  Check statistics
    #    a) for eqtls
    eqtl_labels = sorted_ser.pmarker
 
    gene_labels = sorted_ser.gene
  
    list_genes = sorted_ser['gene'].tolist()
    list_eqtls = sorted_ser['pmarker'].tolist() 


    #   3. Build data frame from expression data
    sel_exp = expr_values[list_genes].copy()
    

    # find and add the other expression data 
    # and add a table that keeps the ordering information.
    other_gene_labels=set(expr_values.columns)-set(list_genes)

    numbers_of_genes = [len(other_gene_labels), len(expr_values.columns), len(list_genes)]
    print('       numbers_of_genes: ', numbers_of_genes)
    print('       genes with eQTLs: ', len(set(list_genes)))
    print('       all genes: ', numbers_of_genes[0]+numbers_of_genes[2])

    list_other_gene_labels=list(other_gene_labels)
    n_eqtls = len(list_genes) 
    sel_exp = pan.concat([sel_exp, expr_values[list_other_gene_labels]], axis=1)
    # for i in range(len(list_other_gene_labels)):
    #     sel_exp.insert( n_eqtls+i, list_other_gene_labels[i],
    #                 expr_values[list_other_gene_labels[i]] )
    





    check_columns = (sel_exp.columns == list_genes + list_other_gene_labels).all()
    print( '    Check if all expected columns are present:', check_columns)
    print( '    number of columns selected:', len(sel_exp.columns) )
    

    print('#    Reordering:')
    print('#    Writing output:')
    #   4. Write output:
    sel_exp.to_csv( os.path.join(out_path , 'expression_reordered_01.csv'))
    np.savetxt(os.path.join( out_path ,'list_genes_in_order_01.csv'), sel_exp.columns.to_numpy(), fmt='%s')

    #
    #   save columns in order
    #
    sel_exp.columns.to_frame().reset_index(drop=True).to_csv( os.path.join(out_path , 'columns_in_order_for_findr_b.csv'), header=1)

    print( '# Done.' )
    # EOF





def main():
    in_path = COMPRESSED_PATH
    out_path = COMPRESSED_PATH
    if not out_path.exists():
        out_path.mkdir(parents=True)
    align(in_path=in_path,
          out_path=out_path)


if __name__== '__main__':
    main()