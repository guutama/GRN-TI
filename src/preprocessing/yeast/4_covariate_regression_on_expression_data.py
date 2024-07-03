#!/usr/bin/env python
# coding: utf-8
#
#  Python script for linear regression on gene expression values with statsmodels.
#
# 0. imports first
#
import pandas as pan
import numpy as np
import seaborn as sns
import os

import matplotlib.pylab as plt

print('# Linear regression on gene expression:')

#out_path = 'plots_pca_expression_covariates/'
import sys 
import os 
import yaml
directories_to_append = [
    os.path.abspath(os.path.join(os.path.dirname(__file__), '..')),

    # Add more directories as needed
]
for directory in directories_to_append:
    sys.path.append(directory)



from path_config import (
    COMPRESSED_PATH,
    YEAST_PATH
)

def covariate_regression(in_path,out_path):
        #  1. Load sorted expression data:
    exp_path = os.path.join(in_path,"expression_reordered_01.csv")
    expr_data = pan.read_csv(exp_path,sep=',')

    # print(expr_data)
    # print(expr_data.columns)
    all_gene_labels = expr_data.columns[1:]

    #     Load genes-eqtls table:
    n_genes=None
    eqtl_path = os.path.join(COMPRESSED_PATH,'strongest_eqtls_r_3columns.csv')
    gene_eqtl_table = pan.read_csv(eqtl_path,
                                #header=0,
                                nrows=n_genes,
                                    delimiter=',')
    # print(gene_eqtl_table)
    # print(gene_eqtl_table.shape)
    n_eqtls = gene_eqtl_table.shape[0]

    print('#    Load covariates:')
    #     Load covariates:
    covar_path = os.path.join(YEAST_PATH,"SI_Data_02_covariates.xlsx")
    covariates = pan.read_excel(covar_path )
    
    covariates.describe()

    #   check ordering of covariates and samples labels
    sample_labels = expr_data[ expr_data.columns[0] ]
    #
    head_of_sample_labels = [ x.split('-')[0] for x in  sample_labels ]
    head_of_covariate_labels = [ x.split('-')[0] for x in  covariates['segregant'] ]
    print(' list comparison of labels :', 
        head_of_covariate_labels == head_of_sample_labels )
    print(' array comparison of labels :',
        (np.array(head_of_covariate_labels) == np.array(head_of_sample_labels)).all() )

    g = sns.jointplot( x='batch', y='OD_covariate', data=covariates,
                    marginal_kws=dict(bins=15),
                    #marginal_kws=dict(bins=15, rug=True),
                    # scatter_kws= {'s':1.2 ,'lw':1.0/8 , 'zorder':2 },
                    #'x_jitter':True
                    #alpha=0.6,
                    #s=1,
                    kind="reg") 
    plot_path = os.path.join(out_path, 'plot_covariates_batch_vs_od_scatter_plot_.png')
    plt.savefig(  plot_path) 



        #
    #   2. prepare data for regression
    #
    # print(covariates.columns)

    covariates['sample'] = covariates.apply (
        lambda row: row['segregant'].split('-')[0], axis=1)
    #   Make dict mapping sample label to batch
    covs=covariates.set_index('sample')['batch'].to_dict()

    #   Show contents of expression data frame:
    unl=expr_data.columns[0]
    #   Add column for sample to expr_data:
    expr_data['sample'] = expr_data.apply (
        lambda row: row[unl].split('-')[0], axis=1)
    #   Add column for batch to expr_data:
    expr_data['batch'] = expr_data.apply (
        lambda row: covs[row['sample'] ], axis=1)
    #   Add column for OD_covariate to expr_data:
    ods=covariates.set_index('sample')['OD_covariate'].to_dict()
    expr_data['OD_covariate'] = expr_data.apply (
        lambda row: ods[row['sample'] ], axis=1)

    print('#    Run regression with Statsmodels:')
    #
    #   3. Run regression with Statsmodels
    #       a) regression on multiple genes
    def my_map_number_to_letters(x):
        number_letter_dict = { str(i) : chr(ord('a') + i) for i in range(10) }
        sl=[ number_letter_dict[j] for j in str(x) ]
        return ''.join(sl)

    all_gene_labels[:200]
    rev_agl_dict = { i: all_gene_labels[i]  for i in range(0, len(all_gene_labels) ) }
    agl_dict = { all_gene_labels[i] : my_map_number_to_letters(i) for i in range(0, len(all_gene_labels) ) }

    import statsmodels.api as sm
    import statsmodels.formula.api as smf

    models = dict()
    residuals_dict = dict()
    df_for_fit = expr_data[ 'OD_covariate batch sample'.split() ]
    df_for_fit = df_for_fit.set_index('sample', drop=False)
    expr_data  = expr_data.set_index('sample', drop=False)

    for gl in list(all_gene_labels)[:]:
        df_for_fit = df_for_fit.assign( gex=expr_data[gl] )
        mod = smf.ols(formula =  ' gex ~ OD_covariate + C(batch)', data=df_for_fit)
        res = mod.fit()
        models[ gl ] = res
        residuals_dict[ gl ] = res.resid

    print(len(models))

    print('#    Writing output:')
    dfres = pan.DataFrame(residuals_dict)
    dfres.to_csv(os.path.join(out_path,'expression_statsmodels_linreg_residuals_01.csv'))
    # print(dfres)
    #   b) example on a single gene
    #mod = smf.ols(formula='YAL027W ~ OD_covariate + C(batch)', data=expr_data)
    #res = mod.fit()
    #print(res.summary())

    print('# Done with regression.')
    # EOF

def main():
    in_path = COMPRESSED_PATH
    covariate_regression(in_path=in_path,out_path=in_path)

if __name__ == '__main__':
    main()
    

