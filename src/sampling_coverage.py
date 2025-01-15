#
#
# Script to identify the sampling coverage by riptide 
# 19.11.2024
# Made by Maxime Lecomte
#
#

import scipy.stats
import riptide_sampling, riptide, cobra
import matplotlib.pyplot as plt
import arviz as az
import os
import pandas as pd
from statistics import mean
import csv
import scipy
import numpy as np


# SAMPLING=10

def plot_samples(data,path):
    print(len(data.keys()))
    az.plot_trace(data,compact=False)
    plt.show()

def compute_R2(data_control,data_trmt,df_trmt_t,mol=str(),dose=str(),reps=str(),dar_path=str()):
    dar_dir = dar_path+"/DAR/"
    if riptide_sampling.create_directory_if_not_exists(dar_dir):
        dar_file = dar_dir+"ctrl_"+dose.lower()+"_"+reps+".txt"
        stats_dar = dar_dir+"KS_ctrl_"+dose.lower()+"_"+reps+".txt"
        print(dar_file)
        not_in_trmt = set()
        not_in_ctrl = set()
        number_of_dar = set()
        with open(dar_file,"w") as f:
            with open(stats_dar,"w") as stats_dar_file:
                for k,v in data_control.items():
                    if k in data_trmt:
                        n_total_ctrl = len(v)
                        n_total_trmt = len(data_trmt[k])
                        n_active_ctrl = len(list(el for el in v if el != 0.))
                        n_active_trmt = len(list(el for el in data_trmt[k] if el != 0.))
                    
                        compute_two_sample_KS(v,data_trmt[k],k,dose,stats_dar_file,tag='both')

                        r2 = ((n_active_ctrl/n_total_ctrl) - (n_active_trmt/n_total_trmt))**2
                        
                        if r2 > 0.2 : 
                            number_of_dar.add(k)
                            f.write(str(k)+'\t'+str(f"{dose} dose treatment: \u2705 \t control treatment: \u2705\t"+str(r2)+'\n'))

                    else:
                        compute_two_sample_KS(v,np.full(len(v),99999),k,dose,stats_dar_file,tag='control')
                        not_in_trmt.add(k)
                        number_of_dar.add(k)
                        f.write(str(k)+'\t'+str(f"{dose} dose treatment: \u274c \t control treatment: \u2705")+'\n')

                for k,v in data_trmt.items():
                    if k not in data_control:
                        compute_two_sample_KS(np.full(len(data_trmt[k]),99999),data_trmt[k],k,dose,stats_dar_file,tag='treatment')
                        not_in_ctrl.add(k)
                        number_of_dar.add(k)
                        f.write(str(k)+'\t'+str(f"{dose} dose treatment: \u2705 \t control treatment: \u274c")+'\n')


                both,all_ctrl,all_trmt,diff_trmt_ctrl,react_count = key_points_of_comparaison(data_control,not_in_trmt,df_trmt_t,data_trmt)

                f.write(f"\nRatio of molecules find in both conditions: {round(both/len(set(data_control.keys())),3)}\n")
                f.write(f"Ratio of molecules not find in the trmt condition but present in the controle condition: {round(100-both/len(set(data_control.keys())),3)}\n")
                f.write(f"Ratio of molecules not find in the controle condition but present in the trmt condition: {round((len(diff_trmt_ctrl)/ len(all_trmt))*100,3)}\n")

                ### compute minimal number of activated reactions, maximal number of activated reaction, average number of activated reaction
                f.write(f"\nThe minimal number of activated reactions with a {dose} dose is {min(list(react_count.values()))}\n")
                f.write(f"The maximal number of activated reactions with a {dose} dose is {max(list(react_count.values()))}\n")
                f.write(f"The average number of activated reactions with a {dose} dose is {round(mean(list(react_count.values())),0)}\n")
                dar_number = len(number_of_dar)
                f.write(f"\nThe number of DAR is {dar_number}\n")
                f.write(f"Total number of reactions: {len(list(data_trmt.keys()))}\n")

                #TODO
                # proportion of DAR reactions in both conditions

                # proporrtion of DAR reactions present onnnnnly in one condition

        return dar_file



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
            data[k] = v
    return data

def read_sample_file(dir_path):
    files = [f for f in os.listdir(dir_path) if 'samples' in f]
    for f in files:
        df=pd.DataFrame()
        df = pd.read_csv(dir_path+'/'+f,sep='\t')
        df_t = df.transpose()
        data = df_to_dict(df)

    return data,df_t,df

def save_file(riptide_object,run):
    output = "results/riptide/sampling_coverage/"+str(SAMPLING)
    if riptide_sampling.create_directory_if_not_exists(output):
        sampling_flux = riptide_object.flux_samples
        sampling_flux.to_csv(output+"/run_"+str(run), sep='\t')
    else:
        print(f"not a valid path: {output}")

def compute_dar_specificity_ratio(dars):
    union_dar_r1_l = dars["amiodarone"]["replicate_0"]['Low'] | dars['valproic acid']["replicate_0"]['Low']
    union_dar_r1_m = dars["amiodarone"]["replicate_0"]['Middle'] | dars['valproic acid']["replicate_0"]['Middle']
    union_dar_r1_h = dars["amiodarone"]["replicate_0"]['High'] | dars['valproic acid']["replicate_0"]['High']
    union_dar_r2_l = dars["amiodarone"]["replicate_1"]['Low'] | dars['valproic acid']["replicate_1"]['Low']
    union_dar_r2_m = dars["amiodarone"]["replicate_1"]['Middle'] | dars['valproic acid']["replicate_1"]['Middle']
    union_dar_r2_h = dars["amiodarone"]["replicate_1"]['High']  | dars['valproic acid']["replicate_1"]['High'] 
    
    intersection_dar_r1_l = dars["amiodarone"]["replicate_0"]['Low'] & dars['valproic acid']["replicate_0"]['Low']
    intersection_dar_r1_m = dars["amiodarone"]["replicate_0"]['Middle'] & dars['valproic acid']["replicate_0"]['Middle']
    intersection_dar_r1_h = dars["amiodarone"]["replicate_0"]['High'] & dars['valproic acid']["replicate_0"]['High']
    intersection_dar_r2_l = dars["amiodarone"]["replicate_1"]['Low'] & dars['valproic acid']["replicate_1"]['Low']
    intersection_dar_r2_m = dars["amiodarone"]["replicate_1"]['Middle'] & dars['valproic acid']["replicate_1"]['Middle']
    intersection_dar_r2_h = dars["amiodarone"]["replicate_1"]['High'] & dars['valproic acid']["replicate_1"]['High']
    
    unique_dar_a_r1_l = dars["amiodarone"]["replicate_0"]['Low'] - dars['valproic acid']["replicate_0"]['Low']
    unique_dar_a_r1_m = dars["amiodarone"]["replicate_0"]['Middle'] - dars['valproic acid']["replicate_0"]['Middle']
    unique_dar_a_r1_h = dars["amiodarone"]["replicate_0"]['High'] - dars['valproic acid']["replicate_0"]['High']
    unique_dar_a_r2_l = dars["amiodarone"]["replicate_1"]['Low'] - dars['valproic acid']["replicate_1"]['Low']
    unique_dar_a_r2_m = dars["amiodarone"]["replicate_1"]['Middle'] - dars['valproic acid']["replicate_1"]['Middle']
    unique_dar_a_r2_h = dars["amiodarone"]["replicate_1"]['High'] - dars['valproic acid']["replicate_1"]['High']
    
    unique_dar_v_r1_l = dars['valproic acid']["replicate_0"]['Low'] - dars["amiodarone"]["replicate_0"]['Low']
    unique_dar_v_r1_m = dars['valproic acid']["replicate_0"]['Middle'] - dars["amiodarone"]["replicate_0"]['Middle']
    unique_dar_v_r1_h = dars['valproic acid']["replicate_0"]['High'] - dars["amiodarone"]["replicate_0"]['High']
    unique_dar_v_r2_l = dars['valproic acid']["replicate_1"]['Low'] - dars["amiodarone"]["replicate_1"]['Low']
    unique_dar_v_r2_m = dars['valproic acid']["replicate_1"]['Middle'] - dars["amiodarone"]["replicate_1"]['Middle']
    unique_dar_v_r2_h = dars['valproic acid']["replicate_1"]['High'] - dars["amiodarone"]["replicate_1"]['High']

    dar_a_r1_l = len(dars["amiodarone"]["replicate_0"]['Low'])
    dar_a_r1_m = len(dars["amiodarone"]["replicate_0"]['Middle'])
    dar_a_r1_h = len(dars["amiodarone"]["replicate_0"]['High'])
    only_dar_a_r1_l = len(unique_dar_a_r1_l)
    only_dar_a_r1_m = len(unique_dar_a_r1_m)
    only_dar_a_r1_h = len(unique_dar_a_r1_h)

    dar_a_r2_l = len(dars["amiodarone"]["replicate_1"]['Low'])
    dar_a_r2_m = len(dars["amiodarone"]["replicate_1"]['Middle'])
    dar_a_r2_h = len(dars["amiodarone"]["replicate_1"]['High'])
    only_dar_a_r2_l = len(unique_dar_a_r2_l)
    only_dar_a_r2_m = len(unique_dar_a_r2_m)
    only_dar_a_r2_h = len(unique_dar_a_r2_h)

    dar_v_r1_l = len(dars["valproic acid"]["replicate_0"]['Low'])
    dar_v_r1_m = len(dars["valproic acid"]["replicate_0"]['Middle'])
    dar_v_r1_h = len(dars["valproic acid"]["replicate_0"])
    only_dar_v_r1_l = len(unique_dar_v_r1_l)
    only_dar_v_r1_m = len(unique_dar_v_r1_m)
    only_dar_v_r1_h = len(unique_dar_v_r1_h)

    dar_v_r2_l = len(dars["valproic acid"]["replicate_1"]['Low'])
    dar_v_r2_m = len(dars["valproic acid"]["replicate_1"]['Middle'])
    dar_v_r2_h = len(dars["valproic acid"]["replicate_1"]['High'])
    only_dar_v_r2_l = len(unique_dar_v_r2_l)
    only_dar_v_r2_m = len(unique_dar_v_r2_m)
    only_dar_v_r2_h = len(unique_dar_v_r2_h)

    print(f"\nRatio of common dars between amiodarone and valproic acid in replicates 1 for low dose comparison: {round((len(intersection_dar_r1_l)/len(union_dar_r1_l))*100,2)}%")
    print(f"Ratio of common dars between amiodarone and valproic acid in replicates 1 for middle dose comparison: {round((len(intersection_dar_r1_m)/len(union_dar_r1_m))*100,2)}%")
    print(f"Ratio of common dars between amiodarone and valproic acid in replicates 1 for high dose comparison: {round((len(intersection_dar_r1_h)/len(union_dar_r1_h))*100,2)}%")

    print(f"Ratio of common dars between amiodarone and valproic acid in replicates 2 for low dose comparison: {round((len(intersection_dar_r2_l)/len(union_dar_r2_l))*100,2)}%")
    print(f"Ratio of common dars between amiodarone and valproic acid in replicates 2 for middle dose comparison: {round((len(intersection_dar_r2_m)/len(union_dar_r2_m))*100,2)}%")
    print(f"Ratio of common dars between amiodarone and valproic acid in replicates 2 for high dose comparison: {round((len(intersection_dar_r2_h)/len(union_dar_r2_h))*100,2)}%")
    print(f"Ratio of unique dar of amiodarone in replicates 1 for low dose comparison: {round((only_dar_a_r1_l/dar_a_r1_l)*100,2)} %")
    print(f"Ratio of unique dar of amiodarone in replicates 1 for middle dose comparison: {round((only_dar_a_r1_m/dar_a_r1_m)*100,2)} %")
    print(f"Ratio of unique dar of amiodarone in replicates 1 for high dose comparison: {round((only_dar_a_r1_h/dar_a_r1_h)*100,2)} %")
    print(f"Ratio of unique dar of amiodarone in replicates 2 for low dose comparison: {round((only_dar_a_r2_l/dar_a_r2_l)*100,2)} %")
    print(f"Ratio of unique dar of amiodarone in replicates 2 for middle dose comparison: {round((only_dar_a_r2_m/dar_a_r2_m)*100,2)} %")
    print(f"Ratio of unique dar of amiodarone in replicates 2 for high dose comparison: {round((only_dar_a_r2_h/dar_a_r2_h)*100,2)} %")
    print(f"Ratio of unique dar of valproic acid in replicates 1 for low dose comparison: {round((only_dar_v_r1_l/dar_v_r1_l)*100,2)}%")
    print(f"Ratio of unique dar of valproic acid in replicates 1 for middle dose comparison: {round((only_dar_v_r1_m/dar_v_r1_m)*100,2)}%")
    print(f"Ratio of unique dar of valproic acid in replicates 1 for high dose comparison: {round((only_dar_v_r1_h/dar_v_r1_h)*100,2)}%")
    print(f"Ratio of unique dar of valproic acid in replicates 2 for low dose comparison: {round((only_dar_v_r2_l/dar_v_r2_l)*100,2)}%")
    print(f"Ratio of unique dar of valproic acid in replicates 2 for middle dose comparison: {round((only_dar_v_r2_m/dar_v_r2_m)*100,2)}%")
    print(f"Ratio of unique dar of valproic acid in replicates 2 for high dose comparison: {round((only_dar_v_r2_h/dar_v_r2_h)*100,2)}%")



def read_dar_file(mol,reps,dose,dars,dar_file):
    with open(dar_file) as f:
        dar_file = csv.reader(f, delimiter="\t")
        for line in dar_file:
            if len(line) > 1:
                dars[mol][reps][dose].add(line[0])
        
    return dars
            
def r2_iteration():
    dars={}
    for mol in ["amiodarone","valproic acid"]:
        dars[mol] = {}
        path_dar = "results/riptide/recon2.2/maxfit/"+mol+"/1000/"
        for reps in ["replicate_0","replicate_1"]:
            dars[mol][reps] = {}
            path_control = "results/riptide/recon2.2/maxfit/"+mol+"/1000/24_Control/"+reps
            data_ctrl,_,df_ctrl = read_sample_file(path_control)
            for dose in ["Low","Middle","High"]:
                dars[mol][reps][dose] = set()
                path_trmt = 'results/riptide/recon2.2/maxfit/'+mol+'/1000/24_'+dose+'/'+reps
                data_trmt,df_trmt_t,df_trmt = read_sample_file(path_trmt)
                # data_trmt,df_trmt_t = read_sample_file('results/riptide/recon2.2/maxfit/'+mol+'/1000/24_'+dose+'/replicate_0')
                # compute_R2(data_ctrl,data_trmt,df_trmt_t,mol=mol,dose=dose,reps="replicate_0")
                dar_file = compute_R2(data_ctrl,data_trmt,df_trmt_t,mol=mol,dose=dose,reps=reps,dar_path=path_dar)
                dars = read_dar_file(mol,reps,dose,dars,dar_file)
                # dars = read_dar_file(mol,"replicate_0",dose,dars)
    compute_dar_specificity_ratio(dars)

def compute_two_sample_KS(ctrl,trmt,k,dose, stat_dar_file,tag=str()):
    res = scipy.stats.ks_2samp(ctrl,trmt,mode = 'auto',alternative="two-sided") ## two_sided: the null hypothesis is that ctrl and trmt are identical (same continuous distribution)
    if res.pvalue < 0.05: # if over 0.05, the null hypothesis is valid, i.e. no difference
        if tag == "both": conditions = f"{dose} dose treatment: \u2705 \t control treatment: \u2705\t"
        elif tag == "control" : conditions = f"{dose} dose treatment: \u274c \t control treatment: \u2705\t"
        else: conditions = f"{dose} dose treatment: \u2705 \t control treatment: \u274c\t"

        stat_dar_file.write(str(k)+'\t'f"{conditions}\t{res.pvalue} \n")

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

r2_iteration()