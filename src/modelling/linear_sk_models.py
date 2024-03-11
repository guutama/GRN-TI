
from sklearn.linear_model import (
    BayesianRidge,
    LassoCV,
    RidgeCV,
    ElasticNetCV
)
from sklearn.ensemble import (
    RandomForestRegressor
)
from sklearn.model_selection import (
    KFold
)
import numpy as np 
import yaml 


params=yaml.safe_load(open("src/param_config.yaml"))["modelling"]['train']

n_splits = params['cross_val_split']
seed = params['random_state']
alphas = np.logspace(-3, 1, 50) 

kf = KFold(n_splits=n_splits, shuffle=True, random_state=seed)
  


elr  = ElasticNetCV(cv=kf, random_state=seed,n_jobs=-1,max_iter=100000,l1_ratio=0.5,alphas=alphas)
brr = BayesianRidge(tol=1e-6, compute_score=True, max_iter=10000)
rir = RidgeCV(alphas=alphas, cv=n_splits) 
lar = LassoCV(alphas=alphas, cv=kf, random_state=seed,max_iter=100000)
rfr = RandomForestRegressor(n_estimators=50,max_depth=10,min_samples_split=4,min_samples_leaf=2,random_state=seed)



regularized_ensemble_models = {
    'lasso':lar,
    'ridge':rir,
    'elastic':elr,
    'bayesian_ridge':brr,
    'forest':rfr
}