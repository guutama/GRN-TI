import sys 
import os 

directories_to_append = [
    os.path.abspath(os.path.join(os.path.dirname(__file__), '..')),

    # Add more directories as needed
]
for directory in directories_to_append:
    sys.path.append(directory)

import pandas as pd 
import numpy as np
import findr
from utils import melt_and_transform
import networkx as nx
import json
import gzip
import pickle

import networkx as nx
import yaml

import random
import json
from path_config import (
        FINDR_PATH,
        PAIRWISE_INFERENCE_DIR,
        NETWORK_INFERENCE_DIR,
        NETWORK_INFO_DIR) 



# Function to create a randomized DAG
def randomize_dag(G):
    # Get a topological sort of the original graph
    topo_sort = list(nx.topological_sort(G))
    
    # Shuffle the topological sort to create a new node order
    random.shuffle(topo_sort)
    
    # Create a mapping from old labels to new labels based on the shuffled order
    mapping = {node: topo_sort[i] for i, node in enumerate(G.nodes())}
    
    # Generate a new DAG using the mapping to relabel nodes
    randomized_G = nx.relabel_nodes(G, mapping)
    
    return randomized_G 


def randomize_dag_with_same_degrees(G):
    # Get the in-degree and out-degree of each node
    in_degrees = dict(G.in_degree())
    out_degrees = dict(G.out_degree())

    # Create a new empty directed graph
    new_G = nx.DiGraph()

    # Add the nodes to the new graph
    new_G.add_nodes_from(G.nodes())

    # Create a list of edges to be added
    potential_edges = []

    for node in G.nodes():
        potential_edges.extend([(node, target) for target in G.nodes() if node != target])

    # Randomly shuffle the potential edges
    random.shuffle(potential_edges)

    # Track the current in-degrees and out-degrees in the new graph
    current_in_degrees = {node: 0 for node in G.nodes()}
    current_out_degrees = {node: 0 for node in G.nodes()}

    # Add edges while maintaining the in-degree and out-degree constraints
    for edge in potential_edges:
        u, v = edge
        if current_out_degrees[u] < out_degrees[u] and current_in_degrees[v] < in_degrees[v]:
            new_G.add_edge(u, v)
            current_out_degrees[u] += 1
            current_in_degrees[v] += 1

            # Ensure the graph remains acyclic
            if not nx.is_directed_acyclic_graph(new_G):
                new_G.remove_edge(u, v)
                current_out_degrees[u] -= 1
                current_in_degrees[v] -= 1

    return new_G
def infer(in_path, out_path):
    if not out_path.exists():
        out_path.mkdir(parents=True)
    l = findr.lib(path=FINDR_PATH,loglv = 6, rs=0,nth=0)
    params = yaml.safe_load(open('src/param_config.yaml'))['structure_learning']
    
    params_network = params['network_inference']
    network_type = params_network['network_type']
    fdr_prior =  float(params_network['fdr_prior'])
    
    # Load the compressed data
    file = os.path.join(in_path, f'edge_posteriors_yeast_{network_type}.csv.gz') 
    edge_posteriors = pd.read_csv(file, compression='gzip')
   
    indices = edge_posteriors.columns.tolist()
    columns = edge_posteriors.columns.tolist()
    #
        
    adjacency = melt_and_transform(input_data=edge_posteriors,
                                    indices= indices,
                                    columns=columns,
                                    id_vars="A",
                                    var_name="B",
                                    value_name="p")

    adjacency = adjacency[adjacency["p"]>fdr_prior]
    adjacency = adjacency.reset_index(drop=True)
    print(f'Using network {network_type}',flush=True)
    print(f'Using threshold: {fdr_prior}',flush=True)
    means = adjacency["p"].mean()
    FDR = 1- means
    


    pairwise_graph = nx.from_pandas_edgelist(
        df=adjacency,
        source="A",
        target="B",
        edge_attr="p",
        create_using= nx.DiGraph
    )

    print(f'Number of nodes after threshold :{pairwise_graph.number_of_nodes()}',flush=True)
    adjacency_matrix = nx.adjacency_matrix(pairwise_graph,weight='p')
    weighted_adjacency = np.array(adjacency_matrix.todense())

    ans = l.netr_one_greedy(weighted_adjacency)
    dag_graph = nx.from_numpy_array(weighted_adjacency * ans['net'].astype(int), create_using=nx.DiGraph)


    nodes_list = list(pairwise_graph.nodes())
    label_mapping = {i: nodes_list[i] for i in range(len(pairwise_graph))}



    dag_graph = nx.relabel_nodes(dag_graph, label_mapping) 
    print("Nodes:", dag_graph,flush=True)
   
    # components = list(nx.weakly_connected_components(dag_graph))

    # # Find the largest component
    # largest_component = max(components, key=len)
    # H = dag_graph.subgraph(largest_component).copy()
    # dag_graph = H 

    out_file = os.path.join(out_path, f"graph_yeast_{network_type}.pickle")
    pickle.dump(dag_graph , open(out_file, 'wb'))





    visualize_path = NETWORK_INFO_DIR
    visualize_save_file = os.path.join(visualize_path,f"graph_yeast_{network_type}.csv")
    nx.write_edgelist(dag_graph,visualize_save_file,data=['attr1','attr2','attr3'])

  

    # # Add the mean in-degree weight to your DataFrame
    # degree_info['in_degree_fdr'] = degree_info['node'].map(mean_in_weights)
    total_nodes = dag_graph.number_of_nodes()
    total_edges = dag_graph.number_of_edges()
    mean_in_degree = sum(d for n, d in dag_graph.in_degree()) / total_nodes
    mean_out_degree = sum(d for n, d in dag_graph.out_degree()) / total_nodes

    # Counting root nodes, leaf nodes, and middle nodes
    root_nodes = sum(1 for n in dag_graph.nodes() if dag_graph.in_degree(n) == 0)
    leaf_nodes = sum(1 for n in dag_graph.nodes() if dag_graph.out_degree(n) == 0)
    middle_nodes = total_nodes - root_nodes - leaf_nodes

    graph_info = {'global_fdr': FDR, 
                  'total_nodes':total_nodes,
                  'total_edges':total_edges,
                  'mean_in_edges': mean_in_degree,
                  'mean_out_edges':mean_out_degree,
                  'num_root_nodes':root_nodes,
                  'num_leaf_nodes':leaf_nodes,
                  'num_intermidate_nodes':middle_nodes}

    # degree_info.to_csv(degree_info_path, index=False)
   
    print(f'{FDR=}',flush=True)
    print(dag_graph,flush=True)
    graph_info_path = os.path.join(visualize_path,f"graph_info_yeast_{network_type}.csv")
    with open(graph_info_path, 'w') as jsonfile:
        json.dump(graph_info, jsonfile)

  


    


#  $GEUV_PRO_PAIR_PATH $GEUV_PRO_DAG_PATH $prior_threshold
if __name__ == '__main__':
    in_path=PAIRWISE_INFERENCE_DIR
    out_path = NETWORK_INFERENCE_DIR


    if not out_path.exists():
        out_path.mkdir(parents=True) 
    if not NETWORK_INFO_DIR.exists():
        NETWORK_INFO_DIR.mkdir(parents=True)
    infer(in_path=in_path, 
          out_path =out_path)