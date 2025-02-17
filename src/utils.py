#
#
# Utils functions 
# 03.02.2025
# Made by Maxime Lecomte
#
#

import os
import pandas as pd
import csv
import arviz as az
from statistics import mean
import matplotlib.pyplot as plt
import cobra


def write_stats_file(trmt_ctrl,not_in_ctrl,not_in_trmt, dar,reac,files,dose,tot_reactions):

    files.write(f"\nProportion of molecules find in both conditions: {round(len(trmt_ctrl)/len(dar)*100,3)}\n")
    files.write(f"Proportion of molecules not find in the trmt condition but present in the control condition: {round(len(not_in_trmt)/len(dar)*100,3)}\n")
    files.write(f"Proportion of molecules not find in the control condition but present in the trmt condition: {round((len(not_in_ctrl)/ len(dar))*100,3)}\n")

    ### compute minimal number of activated reactions, maximal number of activated reaction, average number of activated reaction
    files.write(f"\nThe minimal number of activated reactions with a {dose} dose is {min(list(reac.values()))}\n")
    files.write(f"The maximal number of activated reactions with a {dose} dose is {max(list(reac.values()))}\n")
    files.write(f"The average number of activated reactions with a {dose} dose is {round(mean(list(reac.values())),0)}\n")
    dar_number = len(dar)
    files.write(f"\nThe number of DAR is {dar_number}\n")
    files.write(f"Total number of reactions: {tot_reactions}\n")
    files.write(f"Proportion of DAR: {round(dar_number/tot_reactions*100,3)}\n")

def df_to_dict(df):
    data={}
    tmp = df.to_dict('list')
    for k,v in tmp.items():
        if k != 'Unnamed: 0':
            data[k] = v
    return data


def read_sample_file(dir_path):
    files = [f for f in os.listdir(dir_path) if 'samples' in f]
    for f in files:
        df=pd.DataFrame()
        df = pd.read_csv(dir_path+'/'+f,sep='\t',index_col='Unnamed: 0')
        df_t = df.transpose()
        # data = df_to_dict(df)

    return df_t,df

def read_sample_file_imat(dir_path):
    file = [f for f in os.listdir(dir_path)][0]
    df=pd.DataFrame()
    df = pd.read_csv(dir_path+'/'+file,sep=',',index_col='Unnamed: 0').head(n=1000)
    df_t = df.transpose()
    # data = df_to_dict(df)

    return df_t,df


def read_dar_file(mol,reps,dose,dars,dar_file):
    with open(dar_file) as f:
        dar_file = csv.reader(f, delimiter="\t")
        for line in dar_file:
            if len(line) > 1:
                dars[mol][reps][dose].add(line[0])
    return dars
            
def plot_samples(data,path):
    print(len(data.keys()))
    az.plot_trace(data,compact=False)
    plt.show()

def remove_nan_from_list(list_with_nan):
    return [x for x in list(list_with_nan) if str(x) != 'nan']

def load_model(model_path):
    """ load and create cobra SBML model

    Args:
        model_path (str): Path of the SBML model

    Returns:
        cobra.core.Model
    """    
    return cobra.io.read_sbml_model(model_path)

def create_directory_if_not_exists(path):
    """ Create directory path from input path and check if the directory already exists

    Args:
        path (str): Path of the desired directory

    Returns:
        Bool
    """    

    if not os.path.exists(os.path.dirname(path)):
        os.makedirs(os.path.dirname(path))
        print(f"Le répertoire '{path}' a été créé.")
    else:
        print(f"Le répertoire '{path}' existe déjà.")
    return True


def get_fraction(path):
    """Get fraction of optimum objective value from previous maxfit optimization

    Args:
        path (str): Path of the desired directory

    Returns:
        int: Fraction of optimum objective value
    """    
    if os.path.exists(path):
        f = open(path)
        line = f.readline()
        return line.split(' ')[-1]
    else:
        print(f'No optimized fraction of optimum was find because the path {path} is invalid.\nFraction will be set to the defaut value of 0.8.')
        return 0.8
        