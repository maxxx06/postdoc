#
#
# Script to select desired microarray data 
# 15.01.2025
# Made by Maxime Lecomte
#
#

import os 
import pandas as pd

def data_filter(data_path,attribute_data_path,output_transcriptomic_path,sacrifice_period=str(),dose_level=str(),compound_name=str()):
    """_summary_

    Args:
        data_path (_type_): _description_
        attribute_data_path (_type_): _description_
        output_transcriptomic_path (_type_): _description_
        sacrifice_period (_type_, optional): _description_. Defaults to str().
        dose_level (_type_, optional): _description_. Defaults to str().
        compound_name (_type_, optional): _description_. Defaults to str().

    Returns:
        _type_: _description_
    """    
    if os.path.exists(output_transcriptomic_path):
        return pd.read_csv(output_transcriptomic_path)

    df = pd.read_csv(data_path,sep='\t')
    attribute_data = pd.read_csv(attribute_data_path, encoding_errors='ignore')
    attribute_data = attribute_data.loc[(attribute_data['SPECIES'] == 'Human') 
                                        & (attribute_data['ORGAN'] == 'Liver') 
                                        & (attribute_data['DOSE_LEVEL'] == dose_level) 
                                        & (attribute_data['SACRIFICE_PERIOD'] == sacrifice_period)
                                        & (attribute_data['COMPOUND_NAME'] == compound_name)
                                        ]
    
    attribute_data = attribute_data.drop(["SEX_TYPE","MATERIAL_ID","EXP_ID","GROUP_ID",'STRAIN_TYPE',"ADMINISTRATION_ROUTE_TYPE","ANIMAL_AGE(week)","ARR_DESIGN","EXP_TEST_TYPE","COMPOUND_NO","SINGLE_REPEAT_TYPE"],axis=1)

    cel_file = list(attribute_data["BARCODE"].values+'.CEL')
    cel_file+=["SYMBOL"]
    transcriptomic_data = df.loc[:,df.columns.isin(cel_file)]

    transcriptomic_data.to_csv(output_transcriptomic_path)

    return transcriptomic_data

def map_genes(transcriptomic_data_unmapped,mapped_genes=str(),output=str()):
    """_summary_

    Args:
        transcriptomic_data_unmapped (_type_): _description_
        mapped_genes (_type_, optional): _description_. Defaults to str().
        output (_type_, optional): _description_. Defaults to str().

    Returns:
        _type_: _description_
    """    
    if os.path.exists(output[1]):
        data = df_to_dict(output[0],output[1])
        return data
    
    mapped_genes_dict = {}
    with open(mapped_genes) as f:
        f.readline()
        for line in f:
            splitted_lines = line.split()
            if splitted_lines[0] not in mapped_genes_dict:
                mapped_genes_dict[splitted_lines[1]] = splitted_lines[0]

    for k,v in mapped_genes_dict.items():
        transcriptomic_data_unmapped.replace(k,v,inplace=True)
    transcriptomic_data_mapped = transcriptomic_data_unmapped
    transcriptomic_data_mapped.to_csv(output[0])

    data = df_to_dict(output[0],output[1])

    return data

def df_to_dict(path,output):
    """ Tranform pandas dataFrame to appropriate dictionnary

    Args:
        path (str): _description_
        output (_type_): _description_

    Returns:
        _type_: _description_
    """    
    df = pd.read_csv(path)
    df.set_index("SYMBOL",inplace=True)
    df = df.loc[:,df.columns.str.endswith('.CEL')]
    if output != '':
        df.to_csv(output)

    transcriptomic_data_dict = df.to_dict("index")
    
    data = {}
    for k,v in transcriptomic_data_dict.items():
        data[k] = []
        for kk,vv in v.items():
            data[k].append(vv)
            
    return data