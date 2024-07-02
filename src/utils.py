import os 
import sys

import pandas as pd
import numpy as np
from scipy.stats import spearmanr
import time
from scipy.stats import norm,rankdata

def transform_to_log(arr):
    arr = np.maximum(arr, 0)       # Set negative values to zero
    arr = np.log1p(arr)      # Apply element-wise natural logarithm and transpose
    return arr
def get_filepaths(dir_path, exts):
    """
    Get the filepaths for all files with the given file extentions (exts) in the given directory (dir_path)

    return (str if lists)
    """

    if not isinstance(exts, list):
        exts = [exts]
        filepaths = [os.path.join(dir_path, f) for f in os.listdir(dir_path) if f.endswith(tuple(exts))]
        return sorted(filepaths)

        
def filter_file(filepaths,file_end):
    return (lambda x: x)(*[file for file in filepaths if file.endswith(file_end)])

# def get_samples(dir_path):
#     file      = os.path.join(dir_path,"E-GEUV-2.sdrf.txt")
#     sample_info_data = pd.read_csv(file, sep="\t")
#     populations = ["GBR", "FIN", "CEU", "TSI"] # "YRI",
#     unwanted_col = ['HG00107', 'HG00237', 'NA07000']
#     # Filter sample_info_data based on populations
#     filtered_data = sample_info_data[sample_info_data['Characteristics[population]'].isin(populations)]
#     # Remove unwanted columns from Source Name
#     filtered_data = filtered_data[~filtered_data['Source Name'].isin(unwanted_col)]

#     # Extract unique values for cols_intersect
#     cols_intersect = filtered_data['Source Name'].unique().tolist()
#     return cols_intersect
def get_samples(dir_path,out_path):
    file = os.path.join(dir_path, "E-GEUV-2.sdrf.txt")
    sample_info_data = pd.read_csv(file, sep="\t")
    populations = ["GBR", "FIN", "CEU", "TSI"]  # "YRI",
    unwanted_col = ['HG00107', 'HG00237', 'NA07000']
    # Filter sample_info_data based on populations
    filtered_data = sample_info_data[sample_info_data['Characteristics[population]'].isin(populations)]
    # Remove unwanted columns from Source Name
    filtered_data = filtered_data[~filtered_data['Source Name'].isin(unwanted_col)]

    # Extract unique values for cols_intersect
    cols_intersect = filtered_data['Source Name'].unique().tolist()
    
    # Convert cols_intersect to DataFrame
    cols_intersect_df = pd.DataFrame(cols_intersect, columns=['Source Name'])
    
    # Save cols_intersect_df to a compressed CSV file
    output_file_path = os.path.join(out_path, "sample_names.csv.gz")
    cols_intersect_df.to_csv(output_file_path, sep="\t", index=False, compression="gzip")
    return cols_intersect

def sort_dataframe_by_category(dataframe, value_col, categories_list):
    data = dataframe.copy()
    data["Sorter"] = pd.Categorical(values=data[value_col],categories=categories_list,ordered=True)
    data = data.sort_values("Sorter")
    data = data.reset_index(drop=True)
    del data["Sorter"]
    return data

def get_numpy_array(dataframe, colums_to_keep=None, datatype=np.float64):
    if colums_to_keep is not None:
        data_array = dataframe[colums_to_keep]
        data_array = data_array.to_numpy(dtype=datatype)
    else:
        data_array = dataframe.to_numpy(dtype=datatype)
    return data_array



def multiple_correlation(x, y,method="pearson"):
    print("Started the multiple correlation")
    start = time.time()
    """
    Calculates the multiple correlation coefficient between n independent variables (X) and one dependent variable (Y)
    Input: 
        - x: a 2D array of shape (n_samples, n_features), where n_samples is the number of observations and n_features is the number of independent variables
        - y: a 1D array of shape (n_samples,) representing the dependent variable
    Output: 
        - r: a float representing the multiple correlation coefficient between the independent variables (X) and the dependent variable (Y)
        - r_squared: a float representing the R-squared value of the multiple regression
        - adj_r_squared: a float representing the adjusted R-squared value of the multiple regression
    """
    n_samples, n_features = x.shape
    
    if method == "pearson":
        # Calculate the correlation matrix between all variables
        R = np.corrcoef(np.concatenate((x, y.reshape(-1, 1)), axis=1), rowvar=False)
    elif method == "spearman":
        if n_features == 1:
            R, _ = spearmanr(x.reshape(-1), y)
            return R**2
        else:
            R, _ = spearmanr(np.concatenate((x, y.reshape(-1, 1)), axis=1), axis=0)

   
    # Separate the correlation matrix into the submatrices
    R_xx = R[:-1, :-1]
    R_xy = R[:-1, -1]
    R_yx = R[-1, :-1]
    R_yy = R[-1, -1]

    
    # Add a small positive constant to the diagonal elements of R_xx to ensure that it is positive definite
    reg = 1e-10
    R_xx += np.eye(n_features) * reg

    # Calculate the eigenvalues and eigenvectors of R_xx
    eig_vals, eig_vecs = np.linalg.eig(R_xx)

    # Check if any eigenvalues are negative
    if np.any(eig_vals < 0):
    # Add a larger positive constant to the diagonal elements of R_xx to ensure that it is positive definite
        R_xx += np.eye(n_features) * (reg * 10)
    
    # Recalculate the eigenvalues and eigenvectors of R_xx
    eig_vals, eig_vecs = np.linalg.eig(R_xx)
    
    # Check again if any eigenvalues are negative
    if np.any(eig_vals < 0):
        raise ValueError('The correlation matrix R_xx is not positive definite.')

    if method == "pearson":
        # Try to calculate the inverse of R_xx
        try:
            R_xx_inv = np.linalg.inv(R_xx)
        except np.linalg.LinAlgError:
            # If the inverse cannot be computed, remove a collinear variable and try again
            corr_coeffs = np.abs(np.corrcoef(x, rowvar=False)[-1, :-1])
            idx_to_remove = np.argmax(corr_coeffs)
            x = np.delete(x, idx_to_remove, axis=1)
            return multiple_correlation(x, y)
    elif method == "spearman":
        # Try to calculate the inverse of R_xx
        try:
            R_xx_inv = np.linalg.inv(R_xx)
        except np.linalg.LinAlgError:
            # If the inverse cannot be computed, remove a collinear variable and try again
            corr_coeffs = np.abs(spearmanr(x, y, axis=0)[0][-1, :-1])
            idx_to_remove = np.argmax(corr_coeffs)
            x = np.delete(x, idx_to_remove, axis=1)
            return multiple_correlation(x, y)

    # Calculate the multiple correlation coefficient
    # print(f'{R_yy=}')
    print(f'{R_yy=}')
    r = (np.sqrt(R_yx.dot(R_xx_inv).dot(R_xy)))/R_yy

    # Calculate R-squared
    r_squared = r ** 2

    # Calculate adjusted R-squared
    adj_r_squared = 1 - ( (1 - r_squared) * (n_samples - 1) / (n_samples - n_features - 1) )

    print("total time taken: ", time.time() - start)
    return r_squared 



def rank_inverse_normal_transformation(X, c=3/8):
    """
    Computes the rank based inverse normal transform (INT) for a given dataset X.
    
    Parameters
    ----------
    X : array_like
        2-D array of data with samples in columns.
    c : float, optional
        Parameter that controls the degree of ties in the dataset. Default is 3/8.
        
    Returns
    -------
    ipt : array_like
        2-D array of IPT values with the same shape as X.
    """
    # Get the ranks of the data
    ranks = rankdata(X, axis=1, method='average')
#np.apply_along_axis(np.argsort, 0, np.argsort(X, axis=0))
    
    # Compute the IPT values
    n = X.shape[1]
    ipt = norm.ppf((ranks - c) / (n + 1 - 2*c))
   # $$\text{INT}(W_i) = \Phi^{-1} \left\{ \frac{\text{rank}(W_i) - c}{n+1-2c} \right\}, \quad c \in [0, \frac{1}{2}]$$
    
    return ipt 
def melt_and_transform(input_data, id_vars, indices, columns, var_name, value_name):
    # Check if input_data is a NumPy array and convert to DataFrame
    if isinstance(input_data, np.ndarray):
        input_data = pd.DataFrame(input_data, columns=columns)
    
    # If input_data is a DataFrame, make a copy to avoid modifying the original data
    elif isinstance(input_data, pd.DataFrame):
        input_data = input_data.copy()
    else:
        raise ValueError("Input must be a NumPy array or a pandas DataFrame.")

    # Insert id_vars as first column
    input_data.insert(loc=0, column=id_vars, value=indices)
    
    # Melt DataFrame
    melted_df = pd.melt(input_data, id_vars=id_vars, value_vars=columns, var_name=var_name, value_name=value_name, ignore_index=False)
    
    # Reset index and drop old index
    melted_df.reset_index(drop=True, inplace=True)
    
    return melted_df
