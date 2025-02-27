#
#
# Script to select desired microarray data 
# 15.01.2025
# Made by Maxime Lecomte
#
#

import os 
import pandas as pd

def data_filter(data_path,attribute_data_path,output_transcriptomic_path,sacrifice_period=str(),dose_level=str(),compound_name=str(),tag=str()):
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
        return pd.read_csv(output_transcriptomic_path,index_col=['Unnamed: 0'])
    if tag:
        attribute_data = pd.read_excel(attribute_data_path,"Nomenclature échantillons")
        attribute_data = attribute_data.loc[ (attribute_data['Sample Code (nvlle nomenclature sept21)'].str.split('-',expand=True)[1] == sacrifice_period)
                                            & (attribute_data['Compound'] == compound_name)
                                            ]
        cel_file = [el.replace('/','').replace('-','_') for el in list(attribute_data['Sample Code (nvlle nomenclature sept21)'].values) if 'D'+dose_level+"_"+sacrifice_period in el.replace('/','').replace('-','_')]
        cel_file+=["Gene_Symbol",'ensembl_gene']
        df = pd.read_excel(data_path)

    else:
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
        
        df = pd.read_csv(data_path,sep='\t')

    transcriptomic_data = df.loc[:,df.columns.isin(cel_file)]
    transcriptomic_data.to_csv(output_transcriptomic_path)

    return transcriptomic_data

def map_genes(transcriptomic_data_unmapped,mapped_genes,output,model_version,tag=str()):
    """_summary_

    Args:
        transcriptomic_data_unmapped (_type_): _description_
        mapped_genes (_type_, optional): _description_. Defaults to str().
        output (_type_, optional): _description_. Defaults to str().

    Returns:
        _type_: _description_
    """    
    # if os.path.exists(output):
    #     data = df_to_dict(pd.read_csv(output))
    #     return data
    
    df = pd.read_csv(mapped_genes,sep='\t')
    if tag:
        transcriptomic_data_unmapped.rename(columns={"Gene_Symbol": 'Approved symbol'}, inplace = True)
        df = df[['Approved symbol','HGNC ID']]
    else:
        transcriptomic_data_unmapped.rename(columns={"SYMBOL": 'Approved symbol'}, inplace = True)
        df = df[['Ensembl ID(supplied by Ensembl)','Approved symbol','HGNC ID']]

    df.set_index("Approved symbol",inplace=True)
    transcriptomic_data_unmapped.set_index("Approved symbol",inplace=True)

    transcriptomic_data_mapped = transcriptomic_data_unmapped.join(df) #left join 
    if 'ensembl_gene' in transcriptomic_data_mapped.columns: 
        transcriptomic_data_mapped.rename(columns={"ensembl_gene": 'Ensembl ID(supplied by Ensembl)'}, inplace = True)
        transcriptomic_data_mapped = transcriptomic_data_mapped.iloc[:,[1,2,3,4,0]]
    transcriptomic_data_mapped.to_csv(output)
    data = df_to_dict(transcriptomic_data_mapped,model_version)
    return data

def df_to_dict(df,model_version):
    """ Tranform pandas dataFrame to appropriate dictionnary

    Args:
        path (str): _description_
        output (_type_): _description_

    Returns:
        _type_: _description_
    """   
    select_col_to_drop = ("Ensembl ID(supplied by Ensembl)","HGNC ID") if model_version == "recon2.2" else ("HGNC ID","Ensembl ID(supplied by Ensembl)")
    df.drop([select_col_to_drop[0]],axis=1,inplace=True)
    df.set_index(select_col_to_drop[1],inplace=True)
    df = df[~df.index.duplicated(keep='first')]

    transcriptomic_data_dict = df.to_dict("index")
    data = {}
    for idx,v in transcriptomic_data_dict.items():
        data[idx] = []
        for _,val in v.items():
            if '/' in str(val):
                val=val.split('/')[1]
            data[idx].append(val)
            
    return data