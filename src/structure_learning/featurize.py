import os
import sys

directories_to_append = [
    os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    # Add more directories as needed
]
for directory in directories_to_append:
    sys.path.append(directory)

import pickle
import pandas as pd
import networkx as nx
import numpy as np
from path_config import (
   NETWORKS_PATH,
   SPLIT_PATH,
   ALIGNED_PATH,
   FEATURIZED_PATH


)
import yaml

def get_y_data(y_data_name, expr, sample):
    y = expr[expr["TargetID"] == y_data_name][sample].to_numpy(dtype='float64')
    y = y.reshape(-1,1)
    return y

def get_x_data_exp(x_data_names, expr,sample):
    x_gene_expression = expr[expr["TargetID"].isin(x_data_names)][sample].to_numpy(dtype='float64').T
    return x_gene_expression
  
      
    

def get_x_data_geno(gene_names, genotype, mapping, sample):
    # Ensure gene_name is a list
    if not isinstance(gene_names, list):
        gene_names = [gene_names]

    rna_seq_idx = mapping[mapping["GENE_ID"].isin(gene_names)]
    cis_e_QTL = genotype[genotype["SNP"].isin(rna_seq_idx["SNP_ID"])][sample].reset_index(drop=True).to_numpy(dtype='float64').T
    return cis_e_QTL


   
def get_parent_map(G):
    """
    Function to get the parent of each node in a directed graph based on topological sort.

    Parameters:
    G (nx.DiGraph): A directed graph.

    Returns:
    dict: A dictionary mapping each node to its parent. If a node has no parent, it maps to None.
    """
    # Perform topological sort
    topological_order = list(nx.topological_sort(G))

    # Map to store the parent of each node
    parent_map = {}

    # Iterate through each node in topological order
    for node in topological_order:
        # Get the parents (predecessors) of the node
        parents = list(G.predecessors(node))
        # Map each node to its first parent
        parent_map[node] = parents #if parents else []

    return parent_map

# Function to aggregate parent features
def parent_features(node, graph):
    parents = list(graph.predecessors(node))
    if not parents:
        return []
    parent_features = [graph.nodes[parent]['y'] for parent in parents]
    return parent_features
def featurize(graph_path,data_path,mapping_path, out_path):
    params = yaml.safe_load(open('src/params.yaml'))['structure_learning']
    params_pairwise = params['pairwise_inference']
    network_types = ['p']#params_pairwise['network_type']
    for network_type in network_types:
        print("Using: ",network_type)
        if not out_path.exists():
            out_path.mkdir(parents=True)
        graph_file = os.path.join(graph_path, f"dag_graph_fdr20_{network_type}.pickle")
        with open(graph_file, 'rb') as f:
            G = pickle.load(f)

        sorted_nodes = nx.topological_sort(G)#nx.lexicographical_t
        # create a new DiGraph with the sorted nodes
        G = G.subgraph(sorted_nodes).copy()
        train_best_eQTL_G = G.copy()
        test_best_eQTL_G = G.copy()

        train_all_eQTL_G = G.copy()
        test_all_eQTL_G = G.copy()

        # Load each DataFrame individually
        exp_a_train = pd.read_csv(os.path.join(data_path, "exp_a_train.csv.gz"))
        exp_a_test = pd.read_csv(os.path.join(data_path, "exp_a_test.csv.gz"))

        eqtl_best_train = pd.read_csv(os.path.join(data_path, "eqtl_best_train.csv.gz"))
        eqtl_all_train = pd.read_csv(os.path.join(data_path, "eqtl_all_train.csv.gz"))

        eqtl_best_test = pd.read_csv(os.path.join(data_path, "eqtl_best_test.csv.gz"))
        eqtl_all_test = pd.read_csv(os.path.join(data_path, "eqtl_all_test.csv.gz"))
        mapping = pd.read_csv(os.path.join(mapping_path, "map_best.csv.gz"))
        mapping_all = pd.read_csv(os.path.join(mapping_path, "map_all.csv.gz"))
        train_sample = pd.read_csv(os.path.join(data_path, "train_sample.csv.gz")).iloc[:,0].values
        test_sample = pd.read_csv(os.path.join(data_path, "test_sample.csv.gz")).iloc[:,0].values

        #genotype_data =exp_a_train[exp_a_train['TargetID'] == 'ENSG00000136237.12'][train_sample]
        # y = get_y_data(expr=exp_a_train,y_data_name='ENSG00000136237.12',sample=train_sample)
        

        for node in G.nodes():
            # Find the row in genotype data frame for this node
            obs_train = get_y_data(expr=exp_a_train,y_data_name=node,sample=train_sample)
            obs_test = get_y_data(expr=exp_a_test,y_data_name=node,sample=test_sample)
        
            # Find the row in gene expression data frame for this node
            eqtl_best_train_feature= get_x_data_geno(gene_names=node,genotype=eqtl_best_train,mapping=mapping,sample=train_sample)
            eqtl_best_test_feature= get_x_data_geno(gene_names=node,genotype=eqtl_best_test,mapping=mapping,sample=test_sample)

            eqtl_all_train_feature= get_x_data_geno(gene_names=node,genotype=eqtl_all_train,mapping=mapping_all,sample=train_sample)
            eqtl_all_test_feature= get_x_data_geno(gene_names=node,genotype=eqtl_all_test,mapping=mapping_all,sample=test_sample)
        
         
            

            parents = list(G.predecessors(node))
            num_parent_gene = len(parents)

            if parents:
                parent_train = get_x_data_exp(x_data_names=parents, expr=exp_a_train,sample=train_sample)
                parent_test = get_x_data_exp(x_data_names=parents, expr=exp_a_test,sample=test_sample)

                parent_eqtl_train = get_x_data_geno(gene_names=parents,genotype=eqtl_best_train,mapping=mapping,sample=train_sample)
                parent_eqtl_test = get_x_data_geno(gene_names=parents,genotype=eqtl_best_test,mapping=mapping,sample=test_sample)


                parent_eqtl__all_train = get_x_data_geno(gene_names=parents,genotype=eqtl_all_train,mapping=mapping_all,sample=train_sample)
                parent_eqtl_all_test = get_x_data_geno(gene_names=parents,genotype=eqtl_all_test,mapping=mapping_all,sample=test_sample)

                    # Assign this data to the node
                train_best_eQTL_G.nodes[node]['g_train'] = eqtl_best_train_feature
                train_best_eQTL_G.nodes[node]['g_train_parent']= parent_eqtl_train
                train_best_eQTL_G.nodes[node]['y_train'] = obs_train
                train_best_eQTL_G.nodes[node]['p_train'] = parent_train

                test_best_eQTL_G.nodes[node]['g_test'] = eqtl_best_test_feature
                test_best_eQTL_G.nodes[node]['g_test_parent'] = parent_eqtl_test
                test_best_eQTL_G.nodes[node]['y_test'] = obs_test
                test_best_eQTL_G.nodes[node]['p_test'] = parent_test

                # Assign this data to the node
                train_all_eQTL_G.nodes[node]['g_train'] = eqtl_all_train_feature
                train_all_eQTL_G.nodes[node]['g_train_parent']= parent_eqtl__all_train
                train_all_eQTL_G.nodes[node]['y_train'] = obs_train
                train_all_eQTL_G.nodes[node]['p_train'] = parent_train

                test_all_eQTL_G.nodes[node]['g_test'] = eqtl_all_test_feature
                test_all_eQTL_G.nodes[node]['g_test_parent'] = parent_eqtl_all_test
                test_all_eQTL_G.nodes[node]['y_test'] = obs_test
                test_all_eQTL_G.nodes[node]['p_test'] = parent_test
                # print(f'{eqtl_all_train_feature.shape=}',flush=True)
                print(f'{parent_eqtl__all_train.shape=}',flush=True)
                print()
                
            else:
                # Assign this data to the node
                train_best_eQTL_G.nodes[node]['g_train'] = eqtl_best_train_feature
                train_best_eQTL_G.nodes[node]['g_train_parent']= None
                train_best_eQTL_G.nodes[node]['y_train'] = obs_train
                train_best_eQTL_G.nodes[node]['p_train'] = None

                test_best_eQTL_G.nodes[node]['g_test'] = eqtl_best_test_feature
                test_best_eQTL_G.nodes[node]['g_test_parent'] = None
                test_best_eQTL_G.nodes[node]['y_test'] = obs_test
                test_best_eQTL_G.nodes[node]['p_test'] = None

                # Assign this data to the node
                train_all_eQTL_G.nodes[node]['g_train'] = eqtl_all_train_feature
                train_all_eQTL_G.nodes[node]['g_train_parent']= None
                train_all_eQTL_G.nodes[node]['y_train'] = obs_train
                train_all_eQTL_G.nodes[node]['p_train'] = None

                test_all_eQTL_G.nodes[node]['g_test'] = eqtl_all_test_feature
                test_all_eQTL_G.nodes[node]['g_test_parent'] = None
                test_all_eQTL_G.nodes[node]['y_test'] = obs_test
                test_all_eQTL_G.nodes[node]['p_test'] = None
                # print(f'{eqtl_all_train_feature.shape=}',flush=True)
                # print()
                
            


        



        graphs = {
        "train_best_eQTL_G": train_best_eQTL_G,
        "test_best_eQTL_G": test_best_eQTL_G,
        "train_all_eQTL_G": train_all_eQTL_G,
        "test_all_eQTL_G": test_all_eQTL_G
        }

        # save graph object to file in pickle format
        with open(os.path.join(out_path, f"graphs_fdr20_{network_type}.pkl"), "wb") as f:
            pickle.dump(graphs, f)

    
if __name__ == '__main__':
    featurize(graph_path=NETWORKS_PATH,
         data_path=SPLIT_PATH,
         mapping_path=ALIGNED_PATH,
         out_path=FEATURIZED_PATH)





