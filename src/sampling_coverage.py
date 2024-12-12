#
#
# Script to identify the sampling coverage by riptide 
# 19.11.2024
# Made by Maxime Lecomte
#
#

import riptide_sampling, riptide, cobra
import matplotlib.pyplot as plt
import arviz as az
import os
import pandas as pd
from statistics import mean
import csv

SAMPLING=10

def plot_samples(data,path):
    print(len(data.keys()))
    az.plot_trace(data,compact=False)
    plt.show()

def compute_R2(data_control,data_trtm,df_trmt_t,mol=str(),dose=str()):
    not_in_trmt = set()
    number_of_dar = set()
    with open("results/riptide/recon2.2/"+mol+"/DAR/ctrl_"+dose.lower()+".txt","w") as f:
        for k,v in data_control.items():
            if k in data_trtm:
                n_total_ctrl = len(v)
                n_total_trmt = len(data_trtm[k])
                n_active_ctrl = len(list(el for el in v if el != 0.))
                n_active_trmt = len(list(el for el in data_trtm[k] if el != 0.))
                
                r2 = ((n_active_ctrl/n_total_ctrl) - (n_active_trmt/n_total_trmt))**2
                if r2 > 0.2 : 
                    number_of_dar.add(k)
                    f.write(str(k)+'\t'+str(r2)+'\n')
                    # print(f"{k}\t{r2}")
            else:
                not_in_trmt.add(k)
        
        both,all_ctrl,all_trmt,diff_trmt_ctrl,react_count = key_points_of_comparaison(data_control,not_in_trmt,df_trmt_t,data_trtm)

        f.write(f"\nRatio of molecules find in both conditions: {round(both/len(set(data_control.keys())),3)}\n")
        f.write(f"Ratio of molecules not find in the trmt condition but present in the controle condition: {round(100-both/len(set(data_control.keys())),3)}\n")
        f.write(f"Ratio of molecules not find in the controle condition but present in the trmt condition: {round((len(diff_trmt_ctrl)/ len(all_trmt))*100,3)}\n")

        ### compute minimal number of activated reactions, maximal number of activated reaction, average number of activated reaction
        f.write(f"\nThe minimal number of activated reactions with a high dose is {min(list(react_count.values()))}\n")
        f.write(f"The maximal number of activated reactions with a high dose is {max(list(react_count.values()))}\n")
        f.write(f"The average number of activated reactions with a high dose is {round(mean(list(react_count.values())),0)}\n")
        dar_number = len(number_of_dar)
        f.write(f"\nThe number of DAR is {dar_number}\n")
        f.write(f"Total number of reactions: {len(list(data_trtm.keys()))}\n")



def key_points_of_comparaison(data_control,not_in_trmt,df_trmt_t,data_trmt):
    react_count = dict()

    # number of predicted active reactions 
    df_trmt_t.drop(index=df_trmt_t.index[0], axis=0, inplace=True)
    data_trmt_t = df_to_dict(df_trmt_t)
    for k,v in data_trmt_t.items():
        set_v = set(v)
        filtered_v = list(set_v)
        if 0.0 in filtered_v:
            index = filtered_v.index(0.0)
            filtered_v.pop(index)
        
        if filtered_v: react_count[k]=len(filtered_v)

    # present in both conditions
    both = (len(set(data_control.keys())) - len(not_in_trmt))*100
    # present in the trmt condition but not in the control condition, i.e. the added value of the trmt
    all_ctrl = set(data_control.keys())
    all_trmt = set(data_trmt.keys())
    diff_trmt_ctrl = all_trmt - all_ctrl

    return both,all_ctrl,all_trmt,diff_trmt_ctrl, react_count

def df_to_dict(df):
    data={}
    tmp = df.to_dict('list')
    for k,v in tmp.items():
        if k != 'Unnamed: 0':
            if k not in data:
                data[k] = v
            else:
                for el in v:
                    data[k].append(el)
    return data

def read_sample_file(dir_path):
    files = [f for f in os.listdir(dir_path) if 'samples' in f]
    for f in files:
        df=pd.DataFrame()
        df = pd.read_csv(dir_path+'/'+f,sep='\t')
        df_t = df.transpose()
        data = df_to_dict(df)

    return data,df_t

def save_file(riptide_object,run):
    output = "results/riptide/sampling_coverage/"+str(SAMPLING)
    if riptide_sampling.create_directory_if_not_exists(output):
        sampling_flux = riptide_object.flux_samples
        sampling_flux.to_csv(output+"/run_"+str(run), sep='\t')
    else:
        print(f"not a valid path: {output}")

def compute_dar_specificity_ratio():
    dars= read_dar_file()

    union_dar = dars["amiodarone"] | dars['valproic_acid']
    intersection_dar = dars["amiodarone"] & dars['valproic_acid']
    unique_dar_a = dars["amiodarone"] - dars['valproic_acid']
    unique_dar_v = dars['valproic_acid'] - dars["amiodarone"]

    dar_a = len(dars["amiodarone"])
    only_dar_a = len(unique_dar_a)

    dar_v = len(dars["valproic_acid"])
    only_dar_v = len(unique_dar_v)

    print(f"\nRatio of common dars between amiodarone and valproic acid: {round((len(intersection_dar)/len(union_dar))*100,2)}%")
    print(f"Ratio of unique dar of amiodarone: {round((only_dar_a/dar_a)*100,2)} %")
    print(f"Ratio of unique dar of valproic acid: {round((only_dar_v/dar_v)*100,2)}%")


def read_dar_file():
    dars = {}
    for mol in ["amiodarone", "valproic_acid"]:
        dars[mol] = set()
        path = "results/riptide/recon2.2/"+mol+"/DAR/ctrl_high.txt"
        with open(path) as f:
            dar_file = csv.reader(f, delimiter="\t")
            for line in dar_file:
                if len(line) > 1:
                    dars[mol].add(line[0])
            
    return dars
            
def r2_iteration():
    for mol in ["amiodarone","valproic_acid"]:
        data_ctrl,_ = read_sample_file("results/riptide/recon2.2/"+mol+"/24_Control")
        for dose in ["Low","Middle","High"]:
            data_trmt,df_trmt_t = read_sample_file('results/riptide/recon2.2/'+mol+'/24_'+dose)
            compute_R2(data_ctrl,data_trmt,df_trmt_t,mol=mol,dose=dose)


def assess_riptide_coverage():
    transcript_abundances_rpm = riptide.read_transcription_file('data/tests/transcriptome1.tsv')
    my_model = cobra.io.read_sbml_model('data/tests/model_riptide_complex.sbml')
    all_solutions = []
    for i in range(5):
        # riptide_object = riptide.contextualize(model=my_model, transcriptome=transcript_abundances_rpm,objective=True,fraction=0.8,gpr=True,samples=SAMPLING,direct=False) # if True, Sets previous objective function as a constraint with minimum flux equal to user input fraction. Briefly, if biomass function must be taking into account during the contextualisation step

        riptide_object = riptide_sampling.get_flux_samples('data/microarray/annotated_data_uniq_high_sd_flagged.tsv','data/microarray/open_tggates_cel_file_attribute.csv',sacrific_period='24 hr',dose_level='Control',compound_name='amiodarone',sampling_coverage=True)
        save_file(riptide_object,i)
        data = riptide_object.flux_samples.to_dict('list')
        all_solutions.append(data)
    plot_samples(all_solutions)

# assess_riptide_coverage()
# for p in ["results/riptide/recon2.2/amiodarone/24_Control/","results/riptide/recon2.2/amiodarone/24_High/","results/riptide/recon2.2/amiodarone/24_Low/","results/riptide/recon2.2/amiodarone/24_Middle/"]:
#     data=read_sample_file(p)
#     plot_samples(data,p)

# r2_iteration()
compute_dar_specificity_ratio()