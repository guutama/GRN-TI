
import os
import numpy as np
import pandas as pan
# import findr
# import roman 
import joblib
import json


import sys 
import os 
import yaml
# Append necessary directories to the system path
directories_to_append = [
    os.path.abspath(os.path.join(os.path.dirname(__file__), '..')),
    # Add more directories as needed
]
for directory in directories_to_append:
    sys.path.append(directory) 

import numpy as np
import networkx as nx



import yaml
import pickle
import networkx as nx 
import matplotlib.pyplot as plt
import seaborn as sns
from numpy.random import default_rng
from matplotlib import pyplot as plt

from sklearn.metrics import (
    mean_absolute_percentage_error,
    mean_squared_error,
    r2_score,
    explained_variance_score,
  
    make_scorer
)

from sklearn.model_selection import KFold, cross_val_score
from sklearn.linear_model import BayesianRidge
from path_config import (
    MODELS_PATH,
    METRICES_PATH,
    NETWORKS_PATH,
    SPLIT_PATH,
    ALIGNED_PATH
 
) 
from sklearn.preprocessing import StandardScaler

from sklearn.ensemble import BaggingRegressor 
from sklearn.svm import LinearSVR

from sklearn.feature_selection import SelectFromModel


from utils import (
    rank_inverse_normal_transformation
)
import pandas as pd
def load_graph(graph_path,network_type='p'):
    graph_file = os.path.join(graph_path, f"graph_geuv_{network_type}.pickle")
    with open(graph_file, "rb") as f:
        G = pickle.load(f)
    return G
   
def get_y_data(y_data_name, expr, sample):
    y = expr[expr["GENE_ID"] == y_data_name][sample].to_numpy(dtype='float64')
    y = y.reshape(-1,1)
    return y

def get_x_data_exp(x_data_names, expr,sample):
    x_gene_expression = expr[expr["GENE_ID"].isin(x_data_names)][sample].to_numpy(dtype='float64').T
    return x_gene_expression
  
      
    

def get_x_data_geno(gene_names, genotype, mapping, sample):
    # Ensure gene_name is a list
    if not isinstance(gene_names, list):
        gene_names = [gene_names]

    rna_seq_idx = mapping[mapping["GENE_ID"].isin(gene_names)]
    cis_e_QTL = genotype[genotype["SNP_ID"].isin(rna_seq_idx["SNP_ID"])][sample].reset_index(drop=True).to_numpy(dtype='float64').T
    return cis_e_QTL


def concat_arrays(arr1, arr2):
    return np.concatenate((arr1, arr2), axis=1) 


def save_partial_results(models, models_path, scores, scores_path):
    joblib.dump(models, models_path) 
    with open(scores_path, 'w') as file:
        json.dump(scores, file)

def find_parents_and_grandparents(graph, node):
    # Find parents of the node
    parents = list(graph.predecessors(node))
    
    # Initialize a set for parents and grandparents to avoid duplicates
    ancestors = set(parents)
    
    # Find parents of the parents (grandparents) only if they exist
    for parent in parents:
        grandparents = list(graph.predecessors(parent))
        if grandparents:  # Check if the parent has any predecessors
            ancestors.update(grandparents)
    
    return list(ancestors) 




from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline

from sklearn.model_selection import KFold, cross_val_score
from sklearn.linear_model import BayesianRidge
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
import numpy as np

model_name ='brr'
model_list = {'brr': BayesianRidge(tol=1e-6, compute_score=True, max_iter=100000)}
def cross_validation_model(X, y, n_splits=5, scaling=False):
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)
    
    
    base_model = model_list[model_name]
    if scaling:
        pipeline = Pipeline([
            ('scaler', StandardScaler()),
            ('bayesianridge', base_model)
        ])
    else:
        pipeline = Pipeline([
            ('bayesianridge', base_model)
        ])
    
    # Use cross_val_score to perform cross-validation with R² as the scoring metric
    scores = cross_val_score(pipeline, X, y, cv=kf, scoring='r2')
    
    # Fit the pipeline on the whole dataset (optional, based on your requirements)
    pipeline.fit(X, y)
    print(f"Mean cross-validation R² score: {np.mean(scores)}.3")
    # Return the fitted pipeline and cross-validation scores
    return pipeline




params=yaml.safe_load(open("param_config.yaml"))
param_train = params["modeling"]['train']
params_struc = params["structure_learning"]
params_pairwise = params_struc['pairwise_inference']
params_network =params_struc['network_inference'] 
network_type = params_network['network_type'] 



in_path = SPLIT_PATH
graph_path = NETWORKS_PATH
parameters_dir = MODELS_PATH
metrices_dir = METRICES_PATH

align_path = ALIGNED_PATH



if not parameters_dir.exists():
    parameters_dir.mkdir(parents=True)
if not metrices_dir.exists():
    metrices_dir.mkdir(parents=True) 




name_extention = f'{model_name}_all_eqtl_{network_type}'
out_scores_path = os.path.join(metrices_dir, f"{name_extention}_geuv.json")
out_models_path = os.path.join(parameters_dir, f"{name_extention}_geuv.pkl")

train_sample = pd.read_csv(os.path.join(in_path, "train_sample.csv.gz")).iloc[:,0].values 
test_sample = pd.read_csv(os.path.join(in_path, "test_sample.csv.gz")).iloc[:,0].values 

exp_train = pd.read_csv(os.path.join(in_path, "top_exp_train.csv.gz"))
exp_test = pd.read_csv(os.path.join(in_path, "top_exp_test.csv.gz"))





geno_train = pd.read_csv(os.path.join(in_path, "all_eqtl_train.csv.gz"))
geno_test = pd.read_csv(os.path.join(in_path, "all_eqtl_test.csv.gz"))






# geno = pd.read_csv(os.path.join(in_path, "all_snp.csv.gz"))


mapping = pd.read_csv(os.path.join(ALIGNED_PATH, "all_mapping.csv.gz"))
G = load_graph(graph_path,network_type=network_type)
    
sorted_nodes = list(nx.topological_sort(G))


exp_pred = exp_test.copy()
exp_pred[test_sample] = 0.0
exp_pred[test_sample] = exp_pred[test_sample].astype(float)


scores = {}

models = {}
print(f'{network_type=}')
for i, node in enumerate(sorted_nodes):
    train_data_target = get_y_data(expr=exp_train,y_data_name=node,sample=train_sample).ravel()
    train_data_geno = get_x_data_geno(gene_names=node,genotype=geno_train,mapping=mapping,sample=train_sample)
    parents = list(G.predecessors(node)) 



    model = cross_validation_model(X=train_data_geno,y=train_data_target,use_bagging=False)


        

    test_data_target = get_y_data(expr=exp_test,y_data_name=node,sample=test_sample).ravel()


    test_data_geno = get_x_data_geno(gene_names=node,genotype=geno_test,mapping=mapping,sample=test_sample)
    
    pred = model.predict(test_data_geno)
    pred = pred.astype(np.float64)

    score = r2_score(test_data_target, pred)
    condition = exp_pred["GENE_ID"]==node
    exp_pred.loc[condition,test_sample]=pred
    print(f'{i} out of {len(sorted_nodes)}',flush=True)
    print(f'{score=}',flush=True)
    print()
    print()

    models[node] = {'E': model}
    scores[node] = {'E': score}

    

    if parents:
     
        grandparents =find_parents_and_grandparents(graph=G,node=node)
        
        pred_train = model.predict(train_data_geno)
        residual_target = train_data_target - pred_train
        
        
        train_data_parent_exp = get_x_data_exp(x_data_names=parents, expr=exp_train,sample=train_sample)
        train_data_geno_parent = get_x_data_geno(gene_names=parents,genotype=geno_train,mapping=mapping,sample=train_sample)
        train_data_geno_grandparent = get_x_data_geno(gene_names=grandparents,genotype=geno_train,mapping=mapping,sample=train_sample)

      
        train_data_parent_exp  = train_data_parent_exp

        test_data_parent_exp = get_x_data_exp(x_data_names=parents, expr=exp_test,sample=test_sample)
        test_data_geno_parent = get_x_data_geno(gene_names=parents,genotype=geno_test,mapping=mapping,sample=test_sample)
        test_data_geno_grandparent = get_x_data_geno(gene_names=grandparents,genotype=geno_test,mapping=mapping,sample=test_sample) 

        test_data_parent_exp  = test_data_parent_exp 


        test_data_parent_pred = get_x_data_exp(x_data_names=parents, expr=exp_pred,sample=test_sample) 

        test_data_parent_pred = test_data_parent_pred
        

        model_exp= cross_validation_model(X=train_data_parent_exp,y=residual_target.ravel(),use_bagging=False,scaling=True)
        model_trans = cross_validation_model(X=train_data_geno_parent,y=residual_target.ravel(),use_bagging=False)
        model_grand = cross_validation_model(X=train_data_geno_grandparent,y=residual_target.ravel(),use_bagging=False)


        pred_res = model_exp.predict(test_data_parent_exp)
        pred_pred = model_exp.predict(test_data_parent_pred)
        pred_res_geno = model_trans.predict(test_data_geno_parent)
        pred_res_grand = model_grand.predict(test_data_geno_grandparent)

        pred_tot_exp= pred + pred_res 
        pred_tot_trans = pred + pred_res_geno
        pred_tot_pred = pred + pred_pred 
        pred_tot_grand = pred + pred_res_grand

        pred_tot_exp = pred_tot_exp.ravel()
        pred_tot_trans = pred_tot_trans.ravel()
        pred_tot_pred = pred_tot_pred.ravel()
        pred_tot_grand = pred_tot_grand.ravel()

        score_exp = r2_score(test_data_target,pred_tot_exp)
        score_trans = r2_score(test_data_target,pred_tot_trans)
        score_pred = r2_score(test_data_target,pred_tot_pred)
        score_grand = r2_score(test_data_target,pred_tot_grand) 

    

        print(f'{score_exp=}',flush=True)
        print(f'{score_trans=}',flush=True)
        print(f'{score_pred=}',flush=True)
        print(f'{score_grand=}',flush=True)
        print(flush=True)
        print(flush=True)
        

        models[node].update({'XP': model_exp, 'EP': model_trans})

        scores[node].update({
        'EplusXP': score_exp,
        'EplusXPP': score_pred,
        'EplusEP': score_trans,
        'EplusEGP': score_grand,
         }) # Assuming this is the same as one of the other residuals or needs separate calculation
    
    if i % 100 == 0:
        save_partial_results(models=models,
                                models_path=out_models_path,
                                scores=scores,
                                scores_path=out_scores_path)
save_partial_results(models=models,
                                models_path=out_models_path,
                                scores=scores,
                                scores_path=out_scores_path)
    
