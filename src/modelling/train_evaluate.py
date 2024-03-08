
import math
import torch
import pyro
import gpytorch
from matplotlib import pyplot as plt

from sklearn.feature_selection import SelectFromModel


import os
import sys
import joblib
import json
import multiprocessing
from tqdm import tqdm
directories_to_append = [
    os.path.abspath(os.path.join(os.path.dirname(__file__), '..')),
    os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
    # Add more directories as needed
]
for directory in directories_to_append:
    sys.path.append(directory)

from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import r_regression
from sklearn.linear_model import BayesianRidge
import numpy as np
import networkx as nx
import pyro
import pyro.distributions as dist
import torch
import yaml
import pickle
import networkx as nx 
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.ensemble import RandomForestRegressor

from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import KFold, cross_val_score

from sklearn.metrics import (
    mean_absolute_percentage_error,
    mean_squared_error,
    r2_score,
    explained_variance_score,
    make_scorer
)
from path_config import (
 FEATURIZED_PATH
)
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import RFE
from sklearn.linear_model import LinearRegression
import random
from sklearn.utils import resample

from sklearn.model_selection import cross_val_score
import numpy as np

from numpy.random import default_rng

def plot_and_save_dual_distributions(data1, data2, save_path, title='Distribution Comparison', xlabel='Values', label1='y_true', label2='y_pred', figsize=(8, 7)):
    fig, ax = plt.subplots(figsize=figsize)
    sns.histplot(x=data1, kde=True, label=label1, ax=ax)
    sns.histplot(x=data2, kde=True, label=label2, ax=ax)
    ax.legend()
    ax.set(title=title, xlabel=xlabel)
    plt.savefig(save_path)
    plt.close()


    

def concat_arrays(arr1, arr2):
    return np.concatenate((arr1, arr2), axis=1)



        

def load_graph(graph_path,network_type='p'):
    graph_file = os.path.join(graph_path, f"graphs_{network_type}.pkl")
    with open(graph_file, "rb") as f:
        G = pickle.load(f)
    return G
   

   
def train_and_evaluate_model(X_train, y_train, X_test, y_test, node, input_type, n_splits, seed=0):
    # Initialize the BayesianRidge model
    model = BayesianRidge(tol=1e-6, compute_score=True, max_iter=10000)
    # Fit the model on the entire training dataset
    model.fit(X_train, y_train)

    # Evaluate the model on the test dataset
    pred = model.predict(X_test)
    test_score = explained_variance_score(y_test, pred)
    print(f"{node}: Explained Variance Score {input_type} (Test Set): {test_score}", flush=True)

    return model, pred, test_score

def scale_data(X_train, X_test):
    scaler = StandardScaler()
    scaler.fit(X_train)
    return scaler.transform(X_train), scaler.transform(X_test),scaler

from sklearn.feature_selection import f_regression
def kbest_method(X,Xt,y):
    k=X.shape[0]
    if X.shape[1] > k:
        selector = SelectKBest(score_func=f_regression,k=k)
        selector.fit(X,y)
        X_train = selector.transform(X)
        X_test = selector.transform(Xt)
        return X_train,X_test,selector
    else:
        return X,Xt,None
def process_node(node, G, model_type,seed):
    """
    Process a node in the graph and train the specified model.

    Parameters:
    - node: The node to process.
    - G: The graph.
    - model_type: Type of model to train ('eqtl', 'parent', 'both').
    - DIST_PLOT_PATH: The path to save distribution plots.

    Returns:
    - Trained model and R2 values for training and test sets, or None if Xboth or Xparent is None.
    """
    G_tr = G['train_all_eQTL_G']
    G_te = G['test_all_eQTL_G']
    
    y = np.array(G_tr.nodes[node]['y_train']).reshape(-1,)
    y_test = np.array(G_te.nodes[node]['y_test']).reshape(-1,)
    

    if model_type == 'eqtl_expression':
        Xparent = G_tr.nodes[node].get('p_train')
        if Xparent is None:
            X, X_test = G_tr.nodes[node]['g_train'], G_te.nodes[node]['g_test']
            # X, X_test,selector_g = kbest_method(X=X_g,Xt=X_test_g,y=y)
            scaler_p = None
            scaler_g = None
         
            
        else:
            X_g, X_test_g = G_tr.nodes[node]['g_train'], G_te.nodes[node]['g_test']
            X_p, X_test_p,scaler_p = scale_data(G_tr.nodes[node].get('p_train'), G_te.nodes[node]['p_test'])
            scaler_g = None
            # X_g, X_test_g,selector_g = kbest_method(X=X_g,Xt=X_test_g,y=y)
            # X_p, X_test_p,selector_p = kbest_method(X=X_p,Xt=X_test_p,y=y)

            X = concat_arrays(X_g, X_p)
            X_test = concat_arrays(X_test_g, X_test_p)
    else:
        raise ValueError("Invalid model_type. Choose from 'eqtl', 'parent', 'both'.")

    model,pred, score= train_and_evaluate_model(X_train=X, 
                                                y_train=y, 
                                                X_test=X_test, 
                                                y_test=y_test, 
                                                node=node, 
                                                input_type=model_type, 
                                                n_splits=3, 
                                                seed=seed)
    return model,pred,score, scaler_g,scaler_p



    
def process_node_parallel(node, G, counter, total_nodes, lock):
    result = process_node(node, G)
    with lock:
        counter.value += 1
        progress = (counter.value / total_nodes) * 100
        print(f"Processed {counter.value}/{total_nodes} nodes ({progress:.2f}% complete)",flush=True)

    return node, result
def save_scores_to_json(scores, out_scores_path):
    with open(out_scores_path, 'w') as file:
        json.dump(scores, file)

def save_partial_results(models, r2_test, out_models_path, out_scores_test):
    joblib.dump(models, out_models_path)

    with open(out_scores_test, 'w') as file:
        json.dump(r2_test, file)
if __name__ == '__main__':

    RANDOM_SEED = 8927
    rng = default_rng(RANDOM_SEED)


 
    graph_path=GEUV_SPLIT_FEA_PATH
    models_path=MODELS_PATH
    metrices_path = BRIDGE_R2
    # seed = 42
    # random.seed(seed)
    # np.random.seed(seed)
    if not models_path.exists():
        models_path.mkdir(parents=True)
    if not metrices_path.exists():
        metrices_path.mkdir(parents=True)
    if not DIST_PLOT_PATH.exists():
        DIST_PLOT_PATH.mkdir(parents=True)
 

        
    params =yaml.safe_load(open("src/params.yaml"))["modeling"]
    params_feature = params["feature_selection"]
    params_train = params['train']

    params_struc = yaml.safe_load(open('src/params.yaml'))['structure_learning']
    params_pairwise = params_struc['pairwise_inference']
    network_type = 'fdr20_p' #params_pairwise['network_type']


    number_of_features = params_feature["number_of_features"]
    learning_rate = params_train['learning_rate']
    num_iterations = params_train['num_iterations']
    weight_prior_std = params_train['weight_prior_std']
    num_parallel_processes = 4  
    n_splits = 10        
    print(network_type,flush=True)

    G = load_graph(graph_path,network_type=network_type)
    model_to_train = 'eqtl_expression'
    name_extention = f'bridge_{network_type}_{model_to_train}'

    models = {}
    r2_test = {}
    preds = {}
    scalers_g ={}
    scalers_p ={}
  
    pred_pred_score = {}
    out_scores_test = os.path.join(metrices_path, f"{name_extention}.json")
    out_scores_pred = os.path.join(metrices_path, f"{name_extention}_predicted.json")
    out_models_path = os.path.join(models_path, f"{name_extention}_models.pkl")
    # sorted_nodes = list(nx.topological_sort(G))

    un_sorted_nodes = list(G['train_all_eQTL_G'].nodes())


    


    total_nodes = len(un_sorted_nodes)
    
    for i, node in enumerate(un_sorted_nodes):
        model,pred,score, scaler_g,scaler_p= process_node(node, G,model_type=model_to_train,seed=RANDOM_SEED)
        
        # Update models, r2_train, r2_test as needed
        models[node] = model

        r2_test[node] = score
        preds[node] = pred
        scalers_g[node] =scaler_g
        scalers_p[node] =scaler_p
      



        # Save results periodically (adjust the frequency as needed)
        if (i + 1) % 10 == 0:
            save_partial_results(models, r2_test, out_models_path, out_scores_test)

        # Print progress
        progress = ((i + 1) / total_nodes) * 100
        print(f"Processed {i + 1}/{total_nodes} nodes ({progress:.2f}% complete)", flush=True)

    # Save final results
    save_partial_results(models, r2_test, out_models_path, out_scores_test)


    for i, node in enumerate(un_sorted_nodes):

           
        G_te = G['test_all_eQTL_G']
        Xparent = G_te.nodes[node].get('p_test')
        if Xparent is not None:
            model = models[node]
            parents = list(G_te.predecessors(node))
            scaler_g = scalers_g[node]
            scaler_p =scalers_p[node]
           
            y_test = np.array(G_te.nodes[node]['y_test']).reshape(-1,)
            
            parents_preds = [preds[parent] for parent in parents]
            parents_preds = [arr.reshape(-1,1) for arr in parents_preds]
            parents_preds = np.concatenate(parents_preds,axis=1)
            num_parents = len(parents)
            assert num_parents == Xparent.shape[1]
            parents_preds = parents_preds.reshape(-1, num_parents)

            X_g = G_te.nodes[node]['g_test']
            X_p = scaler_p.transform(parents_preds)

            # X_g = selector_g.transform(X_g) if selector_g is not None else X_g
            # X_p = selector_p.transform(X_p) if selector_p is not None else X_p
            X_test = concat_arrays(X_g, X_p)

            pred = model.predict(X_test)
            test_score = explained_variance_score(y_test, pred)
            print(f"{node}: Average R2  : {test_score}", flush=True)
            pred_pred_score[node]=test_score
    
        # Print R2 scores
    
            

            if (i + 1) % 10 == 0:
                with open(out_scores_pred , 'w') as file:
                    json.dump(pred_pred_score, file)

            progress = ((i + 1) / total_nodes) * 100
            print(f"Processed {i + 1}/{total_nodes} nodes ({progress:.2f}% complete)", flush=True)


        else:
            pred_pred_score[node]=None
    with open(out_scores_pred , 'w') as file:
        json.dump(pred_pred_score, file)




                



        
