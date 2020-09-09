import numpy as np
import pandas as pd


def load_file(file_path, file_type = "csv", sample_axis = "rows"):
    if file_type == "tab":
        data_in = pd.read_table(file_path, sep='\t', index_col=0)
        
    elif file_type == "csv":
        data_in = pd.read_csv(file_path, index_col=0)
    
    else:
        raise Exception("Invalid file type!")
    
    if sample_axis == "rows":  # default: rows
        data_in = data_in.T
        
    print('Operating {} genes and {} samples.'.format(data_in.shape[1], data_in.shape[0]))
    return data_in


def preprocess(data, log1p = True, scale_data = True):
    if log1p:
        data = np.log1p(data)
        
    if scale_data:
        data = (data - data.mean(axis=0))/data.std(axis=0)

    return data
