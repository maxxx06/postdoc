#
#
# Script to identify the sampling coverage by riptide 
# 19.11.2024
# Made by Maxime Lecomte
#
#

import scipy.stats
import riptide_sampling
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import arviz as az
import os
import pandas as pd
from statistics import mean
import csv
import scipy
import numpy as np

MOLECULES = ["amiodarone","valproic acid"]
REPLICATES = ["replicate_0","replicate_1"]
DOSE_LEVEL = ["Low","Middle","High"]


def plot_samples(data,path):
    print(len(data.keys()))
    az.plot_trace(data,compact=False)
    plt.show()

def compute_R2(data_control,data_trmt,df_trmt_t,replicate_analysis,replicate_analysis_ks,mol=str(),dose=str(),reps=str(),dar_path=str()):
    dar_dir = dar_path+"/DAR/"
    if riptide_sampling.create_directory_if_not_exists(dar_dir):
        dar_file = dar_dir+"ctrl_"+dose.lower()+"_"+reps+".txt"
        stats_dar = dar_dir+"KS_ctrl_"+dose.lower()+"_"+reps+".txt"
        not_in_trmt = set()
        not_in_ctrl = set()
        number_of_dar = set()
        number_of_dar_KS = set()
        trmt_and_ctrl = set()
        trmt_and_ctrl_KS = set()
        with open(dar_file,"w") as f:
            with open(stats_dar,"w") as stats_dar_file:
                for k,v in data_control.items():
                    if k in data_trmt:
                        n_total_ctrl = len(v)
                        n_total_trmt = len(data_trmt[k])
                        n_active_ctrl = len(list(el for el in v if el != 0.))
                        n_active_trmt = len(list(el for el in data_trmt[k] if el != 0.))
                    
                        compute_two_sample_KS(v,data_trmt[k],k,dose,stats_dar_file,tag='both')
                        trmt_and_ctrl_KS.add(k)
                        number_of_dar_KS.add(k)
                        replicate_analysis_ks[mol][reps][dose].add(k)
                        
                        r2 = ((n_active_ctrl/n_total_ctrl) - (n_active_trmt/n_total_trmt))**2
                        
                        if r2 > 0.2 : 
                            number_of_dar.add(k)
                            trmt_and_ctrl.add(k)
                            replicate_analysis[mol][reps][dose].add(k)
                            f.write(str(k)+'\t'+str(f"{dose} dose treatment: \u2705 \t control treatment: \u2705\t"+str(r2)+'\n'))

                    else:
                        compute_two_sample_KS(v,np.full(len(v),99999),k,dose,stats_dar_file,tag='control')
                        not_in_trmt.add(k)
                        number_of_dar.add(k)
                        replicate_analysis[mol][reps][dose].add(k)
                        replicate_analysis_ks[mol][reps][dose].add(k)
                        f.write(str(k)+'\t'+str(f"{dose} dose treatment: \u274c \t control treatment: \u2705")+'\n')

                for k,v in data_trmt.items():
                    if k not in data_control:
                        compute_two_sample_KS(np.full(len(data_trmt[k]),99999),data_trmt[k],k,dose,stats_dar_file,tag='treatment')
                        not_in_ctrl.add(k)
                        number_of_dar.add(k)
                        replicate_analysis[mol][reps][dose].add(k)
                        replicate_analysis_ks[mol][reps][dose].add(k)
                        f.write(str(k)+'\t'+str(f"{dose} dose treatment: \u2705 \t control treatment: \u274c")+'\n')

                number_of_dar_KS = number_of_dar_KS | number_of_dar
                
                react_count= key_points_of_comparaison(data_control,not_in_trmt,df_trmt_t,data_trmt)
                tot_reactions = len(list(data_trmt.keys()))+len(data_control.keys())

                write_stats_file(trmt_and_ctrl_KS,not_in_ctrl,not_in_trmt,number_of_dar_KS,react_count,stats_dar_file,dose,tot_reactions)
                write_stats_file(trmt_and_ctrl,not_in_ctrl,not_in_trmt,number_of_dar,react_count,f,dose,tot_reactions)

        return dar_file,stats_dar,replicate_analysis,replicate_analysis_ks



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
    return react_count

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


def compute_dar_specificity_ratio_inter_molecules(dars,reps,dose,df,tag=""):
    union_dar = dars["amiodarone"][reps][dose] | dars['valproic acid'][reps][dose]
    intersection_dar = dars["amiodarone"][reps][dose] & dars['valproic acid'][reps][dose]
    unique_dar_a = dars["amiodarone"][reps][dose] - dars['valproic acid'][reps][dose]
    unique_dar_v = dars["valproic acid"][reps][dose] - dars['amiodarone'][reps][dose]

    plt.figure()
    venn2([dars["amiodarone"][reps][dose],dars['valproic acid'][reps][dose]], set_labels = ("Amiodarone",'Valproic acid'))
    plt.savefig("results/riptide/recon2.2/maxfit/inter_molecules/a_v_"+reps+"_"+dose+'_'+tag+'.png')

    dar_a = len(dars["amiodarone"][reps][dose])
    only_dar_a = len(unique_dar_a)
    dar_v = len(dars["valproic acid"][reps][dose])
    only_dar_v = len(unique_dar_v)

    df[f"union_a_v_{reps}_{dose}"] = pd.Series(list(union_dar))
    df[f"intersection_a_v_{reps}_{dose}"] = pd.Series(list(intersection_dar))
    df[f"unique_a_{reps}_{dose}"] = pd.Series(list(unique_dar_a))
    df[f"unique_v_{reps}_{dose}"] = pd.Series(list(unique_dar_v))

    print(f"\nProportion of common dars between amiodarone and valproic acid in {reps} for {dose} dose comparison: {round((len(intersection_dar)/len(union_dar))*100,2)}% ({len(intersection_dar)})")
    print(f"Proportion of unique dar of amiodarone in {reps} for {dose} dose comparison: {round((only_dar_a/dar_a)*100,2)} % ({only_dar_a})")
    print(f"Proportion of unique dar of valproic acid in {reps} for {dose} dose comparison: {round((only_dar_v/dar_v)*100,2)}% ({only_dar_v})")

    return df


def compute_dar_specificity_ratio_inter_molecules_iMAT(dars,df,tag=""):
    union_dar = dars[tag]["amiodarone"] | dars[tag]["valproic acid"]
    intersection_dar = dars[tag]["amiodarone"] & dars[tag]["valproic acid"]
    unique_dar_a = dars[tag]["amiodarone"] - dars[tag]["valproic acid"]
    unique_dar_v = dars[tag]["valproic acid"] - dars[tag]["amiodarone"]

    plt.figure()
    venn2([dars[tag]["amiodarone"],dars[tag]["valproic acid"]], set_labels = ("Amiodarone",'Valproic acid'))
    plt.savefig("results/iMAT/recon2.2/inter_molecules/a_v_High_"+tag+'.png')

    dar_a = len(dars[tag]["amiodarone"])
    only_dar_a = len(unique_dar_a)
    dar_v = len(dars[tag]["valproic acid"])
    only_dar_v = len(unique_dar_v)

    df[f"union_a_v_High"] = pd.Series(list(union_dar))
    df[f"intersection_a_v_High"] = pd.Series(list(intersection_dar))
    df[f"unique_a_High"] = pd.Series(list(unique_dar_a))
    df[f"unique_v_High"] = pd.Series(list(unique_dar_v))

    print(f"\nProportion of common dars between amiodarone and valproic acid for High dose comparison: {round((len(intersection_dar)/len(union_dar))*100,2)}% ({len(intersection_dar)})")
    print(f"Proportion of unique dar of amiodarone in for High dose comparison: {round((only_dar_a/dar_a)*100,2)} % ({only_dar_a})")
    print(f"Proportion of unique dar of valproic acid for High dose comparison: {round((only_dar_v/dar_v)*100,2)}% ({only_dar_v})")

    return df


def compute_dar_specificity_ratio_intra_molecules(dars,dars_ks,mol,reps,dose,df):
    union_dar = dars[mol][reps][dose] | dars_ks[mol][reps][dose]
    intersection_dar = dars[mol][reps][dose] & dars_ks[mol][reps][dose]
    unique_dar = dars[mol][reps][dose] - dars_ks[mol][reps][dose]
    unique_dar_ks = dars_ks[mol][reps][dose] - dars[mol][reps][dose]

    plt.rcParams['figure.max_open_warning'] = 50
    plt.figure()
    venn2([dars[mol][reps][dose],dars_ks[mol][reps][dose]],set_labels = ("R2","KS"))
    plt.savefig("results/riptide/recon2.2/maxfit/intra_molecules_inter_DAR/"+mol+'_'+reps+"_"+dose+'.png', bbox_inches='tight')

    df[f"union_{mol}_{reps}_{dose}"] = pd.Series(list(union_dar))
    df[f"intersection_{mol}_{reps}_{dose}"] = pd.Series(list(intersection_dar))
    df[f"unique_r2_{mol}_{reps}_{dose}"] = pd.Series(list(unique_dar))
    df[f"unique_ks_{mol}_{reps}_{dose}"] = pd.Series(list(unique_dar_ks))

    print(f"\nProportion of common dars for {mol} in {reps} for {dose} dose comparison: {round((len(intersection_dar)/len(union_dar))*100,2)}% ({len(intersection_dar)})")
    print(f"Proportion of unique dar of {mol} in {reps} for {dose} dose comparison for r2: {round((len(unique_dar)/len(union_dar))*100,2)} % ({len(unique_dar)})")
    print(f"Proportion of unique dar of {mol} in {reps} for {dose} dose comparison for ks: {round((len(unique_dar_ks)/len(union_dar))*100,2)} % ({len(unique_dar_ks)})")

    return df


def read_dar_file(mol,reps,dose,dars,dar_file):
    with open(dar_file) as f:
        dar_file = csv.reader(f, delimiter="\t")
        for line in dar_file:
            if len(line) > 1:
                dars[mol][reps][dose].add(line[0])
    return dars
            
def r2_iteration():
    dars={}
    dars_ks={}
    replicate_analysis = {}
    replicate_analysis_ks = {}
    imat_riptide = {}
    for mol in MOLECULES:
        dars[mol] = {}
        dars_ks[mol] = {}
        replicate_analysis[mol] = {}
        replicate_analysis_ks[mol] = {}
        imat_riptide[mol] = {}

        path_dar = "results/riptide/recon2.2/maxfit/"+mol+"/1000/"
        for reps in REPLICATES:
            dars[mol][reps] = {}
            dars_ks[mol][reps] = {}
            replicate_analysis[mol][reps] = {}
            replicate_analysis_ks[mol][reps] = {}
            imat_riptide[mol][reps] = {}

            path_control = "results/riptide/recon2.2/maxfit/"+mol+"/1000/24_Control/"+reps
            data_ctrl,_,df_ctrl = read_sample_file(path_control)

            for dose in DOSE_LEVEL:
                dars[mol][reps][dose] = set()
                dars_ks[mol][reps][dose] = set()
                replicate_analysis[mol][reps][dose] = set()
                replicate_analysis_ks[mol][reps][dose] = set()
                imat_riptide[mol][reps][dose] = {}


                path_trmt = 'results/riptide/recon2.2/maxfit/'+mol+'/1000/24_'+dose+'/'+reps
                data_trmt,df_trmt_t,df_trmt = read_sample_file(path_trmt)
                dar_file,stats_dar,replicate_analysis,replicate_analysis_ks = compute_R2(data_ctrl,data_trmt,df_trmt_t,replicate_analysis,replicate_analysis_ks,mol=mol,dose=dose,reps=reps,dar_path=path_dar)
                dars = read_dar_file(mol,reps,dose,dars,dar_file)
                dars_ks = read_dar_file(mol,reps,dose,dars_ks,stats_dar)
    
    ## quid replicats ?
    replicat_intra_molecule(replicate_analysis,tag='r2')
    replicat_intra_molecule(replicate_analysis_ks,tag='ks')

    # ## inter molecules & DAR identification methods on RIPTiDe
    inter_molecules_and_dar_with_riptide(dars,reps,dose,dars_ks)
    
    # ## intra molecules & inter DAR identification methods on RIPTiDe
    intra_molecules_and_inter_dar_with_riptide(dars,reps,dose,dars_ks)

    ## inter_molecules & DAR identification methods on iMAT
    inter_molecules_and_dar_with_imat()
    
    ## intra molecules & inter DAR identification methods on iMAT
    intra_molecules_and_inter_dar_with_imat()

    ## intra_molecules & inter DAR identification methods on iMAT and riptide
    
    ## inter_molecules & DAR identification methods on iMAT and riptide

    
def replicat_intra_molecule(replicate_analysis,tag=''):
    all_analysis = {}
    for mol,reps_dars in replicate_analysis.items():
        all_analysis[mol] = {}
        reps = list(reps_dars.keys())
        for rep in range(len(reps)):
            for doses in reps_dars.values():
                for dose in doses.items():
                    all_analysis[mol]["union"] = reps_dars[reps[rep]][dose[0]] | reps_dars[reps[rep+1]][dose[0]]
                    all_analysis[mol]["intersection"] = reps_dars[reps[rep]][dose[0]] & reps_dars[reps[rep+1]][dose[0]]
                    all_analysis[mol]["unique_rep_1"] = reps_dars[reps[rep]][dose[0]] - reps_dars[reps[rep+1]][dose[0]]
                    all_analysis[mol]["unique_rep_2"] = reps_dars[reps[rep+1]][dose[0]] - reps_dars[reps[rep]][dose[0]]
                    plt.figure()
                    venn2([reps_dars[reps[rep]][dose[0]],reps_dars[reps[rep+1]][dose[0]]], set_labels = (reps[rep],reps[rep +1]))
                    plt.savefig("results/riptide/recon2.2/maxfit/intra_molecules_inter_DAR/replicates_"+mol+"_"+dose[0]+"_"+tag+".png")
                    all_analysis_df = pd.DataFrame(all_analysis)
                    all_analysis_df.to_csv("results/riptide/recon2.2/maxfit/intra_molecules_inter_DAR/replicates_"+mol+"_"+dose[0]+"_"+tag+".tsv",sep='\t')
            break


def inter_molecules_and_dar_with_riptide(dars,reps,dose,dars_ks):
    df_r2 = pd.DataFrame()
    df_ks = pd.DataFrame()

    for reps in REPLICATES:
        for dose in DOSE_LEVEL:
            df_r2 = compute_dar_specificity_ratio_inter_molecules(dars,reps,dose,df_r2,tag='r2')
            df_ks = compute_dar_specificity_ratio_inter_molecules(dars_ks,reps,dose,df_ks,tag="ks")

    df_r2 = df_r2.apply(lambda col: col.sort_values().reset_index(drop=True))
    df_ks = df_ks.apply(lambda col: col.sort_values().reset_index(drop=True))
    
    result = []

    for i in range(len(df_ks.columns)):
        intersection = set(df_r2[df_ks.columns[i]]).intersection(set(df_ks[df_ks.columns[i]]))
        diff_r2_ks = set(df_r2[df_ks.columns[i]]).difference(set(df_ks[df_ks.columns[i]]))
        diff_ks_r2 = set(df_ks[df_ks.columns[i]]).difference(set(df_r2[df_ks.columns[i]]))
        result.append({'Colonnes': df_ks.columns[i],'Intersection': list(intersection), 'length intersection': len(intersection), 'difference_r2_ks': list(diff_r2_ks), 'length difference r2 ks': len(diff_r2_ks), 'difference_ks_r2': list(diff_ks_r2), 'length difference ks r2': len(diff_ks_r2)})

    df_comparison = pd.DataFrame(result)
    df_comparison.to_csv("results/riptide/recon2.2/maxfit/inter_molecules/compare_r2_Ks.tsv",sep='\t')

    df_r2.to_csv("results/riptide/recon2.2/maxfit/inter_molecules/r2_dar.tsv",sep='\t')
    df_ks.to_csv("results/riptide/recon2.2/maxfit/inter_molecules/ks_dar.tsv",sep='\t')



def inter_molecules_and_dar_with_imat():
    df_r2 = pd.DataFrame()
    df_ks = pd.DataFrame()
    dars = {}
    dars["r2"] = {}
    dars["ks"] = {}

    for mol in MOLECULES:
        tmp_df_r2 = pd.read_csv("results/iMAT/recon2.2/DAR/ctrl_High_"+mol+"_r2.txt",header=None)
        dars["r2"][mol] = set(tmp_df_r2[0].values)
        tmp_df_ks = pd.read_csv("results/iMAT/recon2.2/DAR/ctrl_High_"+mol+"_ks.txt",header=None)
        dars["ks"][mol] = set(tmp_df_ks[0].values)

    df_r2 = compute_dar_specificity_ratio_inter_molecules_iMAT(dars,df_r2,tag='r2')
    df_ks = compute_dar_specificity_ratio_inter_molecules_iMAT(dars,df_ks,tag="ks")

    df_r2 = df_r2.apply(lambda col: col.sort_values().reset_index(drop=True))
    df_ks = df_ks.apply(lambda col: col.sort_values().reset_index(drop=True))
    
    result = []

    for i in range(len(df_ks.columns)):
        intersection = set(df_r2[df_ks.columns[i]]).intersection(set(df_ks[df_ks.columns[i]]))
        diff_r2_ks = set(df_r2[df_ks.columns[i]]).difference(set(df_ks[df_ks.columns[i]]))
        diff_ks_r2 = set(df_ks[df_ks.columns[i]]).difference(set(df_r2[df_ks.columns[i]]))
        result.append({'Colonnes': df_ks.columns[i],'Intersection': list(intersection), 'length intersection': len(intersection), 'difference_r2_ks': list(diff_r2_ks), 'length difference r2 ks': len(diff_r2_ks), 'difference_ks_r2': list(diff_ks_r2), 'length difference ks r2': len(diff_ks_r2)})

    df_comparison = pd.DataFrame(result)
    df_comparison.to_csv("results/iMAT/recon2.2/inter_molecules/compare_r2_Ks.tsv",sep='\t')

    df_r2.to_csv("results/iMAT/recon2.2/inter_molecules/r2_dar.tsv",sep='\t')
    df_ks.to_csv("results/iMAT/recon2.2/inter_molecules/ks_dar.tsv",sep='\t')


def intra_molecules_and_inter_dar_with_riptide(dars,reps,dose,dars_ks):
    df_mol = pd.DataFrame()

    for mol in MOLECULES:
        for reps in REPLICATES:
            for dose in DOSE_LEVEL:
                df_mol = compute_dar_specificity_ratio_intra_molecules(dars,dars_ks,mol,reps,dose,df_mol)

    df_mol = df_mol.apply(lambda col: col.sort_values().reset_index(drop=True))
    df_mol.to_csv("results/riptide/recon2.2/maxfit/intra_molecules_inter_DAR/r2_ks_dar.tsv",sep='\t')


def intra_molecules_and_inter_dar_with_imat():

    ## no replicates here
    result = []
    for mol in MOLECULES:
        df_r2 = pd.read_csv("results/iMAT/recon2.2/DAR/ctrl_High_"+mol+"_r2.txt",header=None)
        df_ks = pd.read_csv("results/iMAT/recon2.2/DAR/ctrl_High_"+mol+"_ks.txt",header=None)
        df_r2.columns = ['dar']
        df_ks.columns = ['dar']

        intersection = set(df_r2['dar']).intersection(set(df_ks['dar']))
        diff_r2_ks = set(df_r2['dar']).difference(set(df_ks['dar']))
        diff_ks_r2 = set(df_ks['dar']).difference(set(df_r2['dar']))

        plt.figure()
        venn2([set(df_r2['dar']),set(df_ks['dar'])], set_labels = ("R2","KS"))
        plt.savefig("results/iMAT/recon2.2/intra_molecules_inter_DAR/"+mol+"_High.png", bbox_inches='tight')

        result.append({'Intersection': list(intersection), 'length intersection': len(intersection), 'difference_r2_ks': list(diff_r2_ks), 'length difference r2 ks': len(diff_r2_ks), 'difference_ks_r2': list(diff_ks_r2), 'length difference ks r2': len(diff_ks_r2)})

    df_comparison = pd.DataFrame(result)
    df_comparison.to_csv("results/iMAT/recon2.2/intra_molecules_inter_DAR/compare_r2_Ks.tsv",sep='\t')

    df_r2.to_csv("results/iMAT/recon2.2/intra_molecules_inter_DAR/r2_dar.tsv",sep='\t')
    df_ks.to_csv("results/iMAT/recon2.2/intra_molecules_inter_DAR/ks_dar.tsv",sep='\t')




def compute_two_sample_KS(ctrl,trmt,k,dose, stat_dar_file,tag=str()):
    res = scipy.stats.ks_2samp(ctrl,trmt,mode = 'auto',alternative="two-sided") ## two_sided: the null hypothesis is that ctrl and trmt are identical (same continuous distribution)
    if res.pvalue < 0.05: # if over 0.05, the null hypothesis is valid, i.e. no difference
        if tag == "both": conditions = f"{dose} dose treatment: \u2705 \t control treatment: \u2705\t"
        elif tag == "control" : conditions = f"{dose} dose treatment: \u274c \t control treatment: \u2705\t"
        else: conditions = f"{dose} dose treatment: \u2705 \t control treatment: \u274c\t"

        stat_dar_file.write(str(k)+'\t'f"{conditions}\t{res.pvalue} \n")


if __name__ == "__main__":
    r2_iteration()