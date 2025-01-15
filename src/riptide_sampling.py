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
from pathlib import Path


# DOSE_LEVEL = ['Control','Low','Middle','High']
# COMPOUND_NAME = ['valproic acid','amiodarone']

def data_filter(data_path,attribute_data_path,output_transcriptomic_path,sacrifice_period=str(),dose_level=str(),compound_name=str()):

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

def load_model(model_path):
    return cobra.io.read_sbml_model(model_path)

def map_genes(transcriptomic_data_unmapped,mapped_genes=str(),output=str()):
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

def create_directory_if_not_exists(path):
    if is_valid_path(path):
        if not os.path.exists(path):
            os.makedirs(path)
            print(f"Le répertoire '{path}' a été créé.")
        else:
            print(f"Le répertoire '{path}' existe déjà.")
        return True
    else:
        return False

def is_valid_path(path):
    try:
        p = Path(path)
        p.resolve(strict=False)
        return True
    except (OSError, ValueError):
         return False

def get_flux_samples(data_path,attribute_data_path,sacrific_period=str(),dose_level=str(),compound_name=str(),sampling_coverage=False,replicates=int(),samples=int()):
    import dexom_python
    output_transcriptomic_path = 'data/microarray/'+compound_name+'/'+sacrific_period.split(' ')[0]+'_'+dose_level+'/'
    model = load_model('data/metabolic_networks/recon2v2_biomass_corrected_final.sbml')   
    model.name = 'data/metabolic_networks/recon2v2_biomass_corrected_final.sbml'
    # output_riptide_path = "results/riptide/recon2.2/"+compound_name+'/'+sacrific_period.split(' ')[0]+'_'+dose_level+'/replicate_'+str(reps)+'/'
    if create_directory_if_not_exists(output_transcriptomic_path):
        transcriptomic_data = data_filter(data_path,attribute_data_path,output_transcriptomic_path+'transcriptomic_data.csv',sacrifice_period=sacrific_period,dose_level=dose_level,compound_name=compound_name)
        transcriptomic_data_mapped = map_genes(transcriptomic_data,mapped_genes='data/microarray/gene_with_protein_product.txt',output=[output_transcriptomic_path+'transcriptomic_data_mapped.csv',output_transcriptomic_path+'transcriptomic_data_mapped_corrected.csv'])
        # transcriptomic_data_mapped_normalized = riptide.read_transcription_file(output_transcriptomic_path+'transcriptomic_data_mapped_corrected.csv',sep=',',header=True)

        # reaction_weights = dexom_python.apply_gpr(model,transcriptomic_data_mapped_reps,save=False)
        # print(reaction_weights)
        # dexom_solutions  = dexom_python.imat(model,reaction_weights)
        # print(dexom_solutions.fluxes["biomass_reaction"])
        # exit()

        for reps in range(replicates):
            transcriptomic_data_mapped_reps = {k:v[reps] for k,v in transcriptomic_data_mapped.items()}
            out = "results/riptide/recon2.2/maxfit/"+compound_name+'/'+str(samples)+'/'+sacrific_period.split(' ')[0]+'_'+dose_level+'/replicate_'+str(reps)+'/'

            max_fit = riptide.maxfit(model=model,transcriptome=transcriptomic_data_mapped_reps,objective=True,prune=True,gpr=True,samples=samples)
            if create_directory_if_not_exists(out):
                riptide.save_output(riptide_obj=max_fit, path=out,file_type='SBML')
        # reps = "NAN"
        # output_riptide_path = "results/riptide/recon2.2/"+compound_name+'/1000/'+sacrific_period.split(' ')[0]+'_'+dose_level+'/replicate_'+str(reps)+'/'

        # riptide_object = riptide.contextualize(model,transcriptome = transcriptomic_data_mapped,gpr=True,prune=True,objective=True,samples=1000,fraction=0.8)
        # if sampling_coverage:
        #     return riptide_object
        # if create_directory_if_not_exists(output_riptide_path):
        #     riptide.save_output(riptide_obj=riptide_object, path=output_riptide_path,file_type='SBML')
        # # print(riptide_object.flux_samples)
    else:
        print('not a valid path for ', output_transcriptomic_path)


def imat_solutions():
    from dexom_python import imat,apply_gpr
    model =  load_model("data/metabolic_networks/recon2v2_biomass_corrected.sbml")
    #define parameters
    eps = 1e-2  # threshold of activity for highly expressed reactions
    thr = 1e-5  # threshold of activity for all reactions
    obj_tol = 1e-5  # variance allowed for the objective_value
    tlim = 600  # time limit (in seconds) for the imat model.optimisation() call
    tol = 1e-6  # tolerance for the solver
    maxiter = 10
    model.solver = "cplex"

    #define reactions weights
    df = pd.read_csv('data/microarray/24_Control_amiodarone/transcriptomic_data_mapped_corrected.csv')
    transcriptomic_data_mapped = df_to_dict('data/microarray/24_Control_amiodarone/transcriptomic_data_mapped_corrected.csv','')
    transcriptomic_data_mapped = {k:v[0] for k,v in transcriptomic_data_mapped.items()}
    reaction_weights = apply_gpr(model,transcriptomic_data_mapped)

    #run imat
    dexom_solutions  = imat(model,reaction_weights,epsilon=eps)
    print(dexom_solutions.fluxes)


# imat_solutions() ## do not forget to transform INF bounds into 1000 and -1000
# for cpd in ['amiodarone','valproic acid']:
#     for dose in ['Control','Low','Middle','High']:
        # get_flux_samples('data/microarray/annotated_data_uniq_high_sd_flagged.tsv','data/microarray/open_tggates_cel_file_attribute.csv',sacrific_period='24 hr',dose_level=dose,compound_name=cpd)
