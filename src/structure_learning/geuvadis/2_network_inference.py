import sys 
import os 
import pandas as pd 
import numpy as np
import findr
from utils import melt_and_transform
import networkx as nx
import json
import gzip
import pickle
import yaml
import random
from path_config import (
    FINDR_PATH,
    PAIRWISE_INFERENCE_PATH,
    NETWORKS_PATH,
    NETWORKS_INFO_PATH
) 

# Append necessary directories to the system path
directories_to_append = [
    os.path.abspath(os.path.join(os.path.dirname(__file__), '..')),
    # Add more directories as needed
]
for directory in directories_to_append:
    sys.path.append(directory)

def randomize_dag(G):
    """
    Create a randomized Directed Acyclic Graph (DAG) by shuffling the topological order of nodes.
    
    Parameters:
    G (networkx.DiGraph): Input DAG.
    
    Returns:
    networkx.DiGraph: Randomized DAG.
    """
    topo_sort = list(nx.topological_sort(G))
    random.shuffle(topo_sort)
    mapping = {node: topo_sort[i] for i, node in enumerate(G.nodes())}
    randomized_G = nx.relabel_nodes(G, mapping)
    return randomized_G 



def infer(in_path, out_path):
    """
    Perform network inference to calculate directed acyclic graph (DAG) from edge posteriors.
    
    Parameters:
    in_path (str): Input directory path containing edge posterior data.
    out_path (str): Output directory path to save the results.
    """
    # Create output directory if it doesn't exist
    if not out_path.exists():
        out_path.mkdir(parents=True)

    # Initialize Findr library
    l = findr.lib(path=FINDR_PATH, loglv=6, rs=0, nth=0)

    # Load configuration parameters for structure learning
    params = yaml.safe_load(open('src/config_params.yaml'))['structure_learning']
    params_network = params['network_inference']
    network_type = params_network['network_type']
    fdr_prior = float(params_network['fdr_prior'])
    
    # Load the compressed data
    file = os.path.join(in_path, f'edge_posteriors_{network_type}.csv.gz')
    edge_posteriors = pd.read_csv(file, compression='gzip')
   
    indices = edge_posteriors.columns.tolist()
    columns = edge_posteriors.columns.tolist()
    
    # Transform the data for adjacency matrix
    adjacency = melt_and_transform(input_data=edge_posteriors,
                                   indices=indices,
                                   columns=columns,
                                   id_vars="A",
                                   var_name="B",
                                   value_name="p")
    
    # Apply threshold based on FDR prior
    adjacency = adjacency[adjacency["p"] > fdr_prior]
    adjacency = adjacency.reset_index(drop=True)
    print(f'Using network {network_type}', flush=True)
    print(f'Using threshold: {fdr_prior}')
    means = adjacency["p"].mean()
    FDR = 1 - means
    
    # Create a directed graph from the adjacency matrix
    pairwise_graph = nx.from_pandas_edgelist(
        df=adjacency,
        source="A",
        target="B",
        edge_attr="p",
        create_using=nx.DiGraph
    )

    print(f'Number of nodes after threshold: {pairwise_graph.number_of_nodes()}')
    adjacency_matrix = nx.adjacency_matrix(pairwise_graph, weight='p')
    weighted_adjacency = np.array(adjacency_matrix.todense())

    # Perform greedy network inference
    ans = l.netr_one_greedy(weighted_adjacency)
    dag_graph = nx.from_numpy_array(weighted_adjacency * ans['net'].astype(int), create_using=nx.DiGraph)

    nodes_list = list(pairwise_graph.nodes())
    label_mapping = {i: nodes_list[i] for i in range(len(pairwise_graph))}

    # Relabel the nodes in the DAG
    dag_graph = nx.relabel_nodes(dag_graph, label_mapping)
    print("Nodes:", dag_graph)
   
    # Extract the largest weakly connected component
    components = list(nx.weakly_connected_components(dag_graph))
    largest_component = max(components, key=len)
    H = dag_graph.subgraph(largest_component).copy()
    dag_graph = H 

    # Save the graph object to file in pickle format
    out_file = os.path.join(out_path, f"graph_{network_type}.pickle")
    pickle.dump(dag_graph, open(out_file, 'wb'))
    print("Nodes:", dag_graph)

    visualize_path = NETWORKS_INFO_PATH
    visualize_save_file = os.path.join(visualize_path, f"graph_{network_type}.csv")
    nx.write_edgelist(dag_graph, visualize_save_file, data=['attr1', 'attr2', 'attr3'])

    # Calculate graph statistics
    total_nodes = dag_graph.number_of_nodes()
    total_edges = dag_graph.number_of_edges()
    mean_in_degree = sum(d for n, d in dag_graph.in_degree()) / total_nodes
    mean_out_degree = sum(d for n, d in dag_graph.out_degree()) / total_nodes

    # Counting root nodes, leaf nodes, and middle nodes
    root_nodes = sum(1 for n in dag_graph.nodes() if dag_graph.in_degree(n) == 0)
    leaf_nodes = sum(1 for n in dag_graph.nodes() if dag_graph.out_degree(n) == 0)
    middle_nodes = total_nodes - root_nodes - leaf_nodes

    # Store graph information
    graph_info = {
        'global_fdr': FDR, 
        'total_nodes': total_nodes,
        'total_edges': total_edges,
        'mean_in_edges': mean_in_degree,
        'mean_out_edges': mean_out_degree,
        'num_root_nodes': root_nodes,
        'num_leaf_nodes': leaf_nodes,
        'num_intermidate_nodes': middle_nodes
    }

    print(f'{FDR=}')
    print(dag_graph)
    graph_info_path = os.path.join(visualize_path, f"graph_info_{network_type}.json")
    with open(graph_info_path, 'w') as jsonfile:
        json.dump(graph_info, jsonfile)

if __name__ == '__main__':
    in_path = PAIRWISE_INFERENCE_PATH
    out_path = NETWORKS_PATH

    if not out_path.exists():
        out_path.mkdir(parents=True) 
    if not NETWORKS_INFO_PATH.exists():
        NETWORKS_INFO_PATH.mkdir(parents=True)
        
    infer(in_path=in_path, out_path=out_path)
