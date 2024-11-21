#
#
# Script to use riptide on microarray data 
# 19.11.2024
# Made by Maxime Lecomte
#
#

import cobra
import pandas as pd
import os 
import riptide

DOSE_LEVEL = ['Control','Low','Middle','High']
COMPOUND_NAME = ['valproic acid','amiodarone']

def data_filter(data_path,attribute_data_path,output_transcriptomic_path):

    if os.path.exists(output_transcriptomic_path):
        return pd.read_csv(output_transcriptomic_path)

    df = pd.read_csv(data_path,sep='\t')
    attribute_data = pd.read_csv(attribute_data_path, encoding_errors='ignore')
    attribute_data = attribute_data.loc[(attribute_data['SPECIES'] == 'Human') 
                                        & (attribute_data['ORGAN'] == 'Liver') 
                                        & (attribute_data['DOSE_LEVEL'].isin(DOSE_LEVEL)) 
                                        & (attribute_data['COMPOUND_NAME'].isin(COMPOUND_NAME))
                                        ]
    
    attribute_data = attribute_data.drop(["SEX_TYPE","MATERIAL_ID","EXP_ID","GROUP_ID",'STRAIN_TYPE',"ADMINISTRATION_ROUTE_TYPE","ANIMAL_AGE(week)","ARR_DESIGN","EXP_TEST_TYPE","COMPOUND_NO","SINGLE_REPEAT_TYPE"],axis=1)

    cel_file = list(attribute_data["BARCODE"].values+'.CEL')
    cel_file+=["SYMBOL"]
    transcriptomic_data = df.loc[:,df.columns.isin(cel_file)]

    transcriptomic_data.to_csv(output_transcriptomic_path)

    return transcriptomic_data

def load_model(model_path):
    return cobra.io.read_sbml_model(model_path)

def map_genes(transcriptomic_data_unmapped,mapped_genes=str(),output=str()):
    if os.path.exists(output):
        data = df_to_dict(output)
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

    transcriptomic_data_mapped.to_csv(output)
    data = df_to_dict(output)

    return data

def df_to_dict(path):
    df = pd.read_csv(path)
    df.set_index("SYMBOL",inplace=True)
    transcriptomic_data_dict = df.to_dict("index")
    
    data = {}
    for k,v in transcriptomic_data_dict.items():
        data[k] = []
        for kk,vv in v.items():
            data[k].append(vv)

    return data

     

def get_flux_samples(data_path,attribute_data_path,output_transcriptomic_path):
    transcriptomic_data = data_filter(data_path,attribute_data_path,output_transcriptomic_path)
    transcriptomic_data_mapped = map_genes(transcriptomic_data,mapped_genes='data/microarray/gene_with_protein_product.txt',output='data/microarray/transcriptomic_data_mapped.csv')
    model = load_model('data/metabolic_networks/recon2.2.xml')
    riptide_object = riptide.contextualize(model,transcriptome = transcriptomic_data_mapped,gpr=True,prune=True,objective=False,samples=1,fraction=0.001)
    riptide.save_output(riptide_obj=riptide_object, path="results/riptide_recon2.2",file_type='SBML')


get_flux_samples('data/microarray/annotated_data_uniq_high_sd_flagged.tsv','data/microarray/open_tggates_cel_file_attribute.csv','data/microarray/transcriptomic_data.csv')
