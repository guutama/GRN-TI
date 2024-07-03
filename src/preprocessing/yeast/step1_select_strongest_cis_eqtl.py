#!/usr/bin/env python
# coding: utf-8



# This script is adapted from the original found at:
# https://github.com/michoel-lab/FindrCausalNetworkInferenceOnYeast
# It has been modified with permission to suit our specific needs.


#
# Python script to load and process yeast eQTLs from Albert et al. (2018):
#
#   Albert, F. W., Bloom, J. S., Siegel, J., Day, L., & Kruglyak, L. (2018).
#       Genetics of trans-regulatory variation in gene expression. Elife, 7, e35471. doi:10.7554/eLife.35471
#
#   paper:  https://doi.org/10.7554/elife.35471
#   data:   https://figshare.com/s/83bddc1ddf3f97108ad4
#


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


import numpy as np 

from path_config import (
    YEAST_PATH,
    COMPRESSED_PATH
)


in_path = YEAST_PATH
out_path = COMPRESSED_PATH


def eqtl_selector(in_path):
    path = os.path.join(in_path,"SI_Data_04_eQTL.xlsx")
    df_eqtl_info = pd.read_excel(path)
    
    #  2. Select cis eQTLs
    cis_indices = df_eqtl_info['cis']==True
    df_cis=df_eqtl_info[cis_indices]
    df_cis =df_cis.drop_duplicates('pmarker')
    df_cis = df_cis.sort_values(by='pmarker')
    df_cis = df_cis.reset_index(drop=True)
 


    
    df_cis.to_csv(os.path.join(out_path,'cis_eqtls_sorted_by_r.csv'),index=False)
    #   write output to csv file:
    df_cis.to_csv(os.path.join(out_path,'cis_eqtls_sorted_by_r_threecolumns.csv'), columns=['gene','pmarker','r'],index=False)

 

    sel_max = df_cis.groupby('gene').apply(
    lambda x: x.loc[x['r'].abs().idxmax()])

    # Reset index if needed
    sel_max.reset_index(drop=True, inplace=True)

    sel_max  =  sel_max .sort_values(by='pmarker')
    sel_max = sel_max.reset_index(drop=True)
 
    #   write output csv files:
    sel_max.to_csv(os.path.join(out_path,'strongest_eqtls_r.csv'),index=False)
    sel_max.to_csv(os.path.join(out_path,'strongest_eqtls_r_3columns.csv'),columns=['gene','pmarker','r'],index=False)
    sel_max.to_csv(os.path.join(out_path,'strongest_eqtls_r_list.csv'),columns=['gene','pmarker'],index=False)

   

   

    # EOF






def main():
    in_path = YEAST_STAT_DIR
    eqtl_selector(in_path=in_path)


if __name__ == '__main__':
    main()
