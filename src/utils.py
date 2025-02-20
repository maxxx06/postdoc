#
#
# Utils functions 
# 03.02.2025
# Made by Maxime Lecomte
#
#

import os
import pandas as pd
import cobra

def read_sample_file(dir_path,tool = str()):
    """
    Reads a sample file from the specified directory and returns two DataFrames.
    
    Based on the specified tool ('mana' or 'riptide'), the function identifies the sample file 
    in the directory, reads it into a DataFrame, transposes the DataFrame, and returns both 
    the transposed and original DataFrames.

    Args:
        dir_path (str): The directory path where the sample files are located.
        tool (str, optional): A string specifying the tool type ('mana' or 'riptide'). 
                              Defaults to an empty string.

    Returns:
        tuple: A tuple containing two pandas DataFrames:
            - df_t (DataFrame): The transposed version of the original DataFrame.
            - df (DataFrame): The original DataFrame read from the sample file.
    """
    if tool == 'mana':
        sample_file = [f for f in os.listdir(dir_path)][0]
    elif tool == 'riptide':
        sample_file = [f for f in os.listdir(dir_path) if 'samples' in f][0]

    df = pd.DataFrame()
    if tool == 'riptide':
        sep = '\t'
    elif tool == 'mana':
        sep=','
    df = pd.read_csv(dir_path+'/'+sample_file,sep=sep,index_col='Unnamed: 0')
    df_t = df.transpose()

    return df_t,df

def remove_nan_from_list(list_with_nan):
    """
    Removes any 'nan' values from the provided list.

    This function filters out elements from the input list that are equal to 'nan' (as a string).

    Args:
        list_with_nan (list): A list that may contain 'nan' values (represented as strings).

    Returns:
        list: A new list containing only the elements from the original list that are not 'nan'.
    
    Example:
        >>> remove_nan_from_list([1, 'nan', 3, 'nan', 5])
        [1, 3, 5]
    """
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
        