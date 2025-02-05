#
#
# Script to compute DARs metrics
# 04.02.2025
# Made by Maxime Lecomte
#
#

import utils
import scipy.stats
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import numpy as np
import pandas as pd


def compute_R2(data_control,data_trmt,df_trmt_t,replicate_analysis,replicate_analysis_ks,mol=str(),dose=str(),reps=str(),dar_path=str()):
    dar_dir = dar_path+"/DAR/"
    if utils.create_directory_if_not_exists(dar_dir):
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

                utils.write_stats_file(trmt_and_ctrl_KS,not_in_ctrl,not_in_trmt,number_of_dar_KS,react_count,stats_dar_file,dose,tot_reactions)
                utils.write_stats_file(trmt_and_ctrl,not_in_ctrl,not_in_trmt,number_of_dar,react_count,f,dose,tot_reactions)

        return dar_file,stats_dar,replicate_analysis,replicate_analysis_ks

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

    union_dar = [x for x in list(union_dar) if str(x) != 'nan']
    intersection_dar = [x for x in list(intersection_dar) if str(x) != 'nan']
    unique_dar_a = [x for x in list(unique_dar_a) if str(x) != 'nan']
    unique_dar_v = [x for x in list(unique_dar_v) if str(x) != 'nan']

    df[f"union_a_v_{reps}_{dose}"] = pd.Series(union_dar)
    df[f"intersection_a_v_{reps}_{dose}"] = pd.Series(intersection_dar)
    df[f"unique_a_{reps}_{dose}"] = pd.Series(unique_dar_a)
    df[f"unique_v_{reps}_{dose}"] = pd.Series(unique_dar_v)


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

    union_dar = [x for x in list(union_dar) if str(x) != 'nan']
    intersection_dar = [x for x in list(intersection_dar) if str(x) != 'nan']
    unique_dar_a = [x for x in list(unique_dar_a) if str(x) != 'nan']
    unique_dar_v = [x for x in list(unique_dar_v) if str(x) != 'nan']

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

    plt.rcParams['figure.max_open_warning'] = 100
    plt.figure()
    venn2([dars[mol][reps][dose],dars_ks[mol][reps][dose]],set_labels = ("R2","KS"))
    plt.savefig("results/riptide/recon2.2/maxfit/intra_molecules_inter_DAR/"+mol+'_'+reps+"_"+dose+'.png', bbox_inches='tight')

    union_dar = [x for x in list(union_dar) if str(x) != 'nan']
    intersection_dar = [x for x in list(intersection_dar) if str(x) != 'nan']
    unique_dar = [x for x in list(unique_dar) if str(x) != 'nan']
    unique_dar_ks = [x for x in list(unique_dar_ks) if str(x) != 'nan']

    df[f"union_{mol}_{reps}_{dose}"] = pd.Series(list(union_dar))
    df[f"intersection_{mol}_{reps}_{dose}"] = pd.Series(list(intersection_dar))
    df[f"unique_r2_{mol}_{reps}_{dose}"] = pd.Series(list(unique_dar))
    df[f"unique_ks_{mol}_{reps}_{dose}"] = pd.Series(list(unique_dar_ks))

    print(f"\nProportion of common dars for {mol} in {reps} for {dose} dose comparison: {round((len(intersection_dar)/len(union_dar))*100,2)}% ({len(intersection_dar)})")
    print(f"Proportion of unique dar of {mol} in {reps} for {dose} dose comparison for r2: {round((len(unique_dar)/len(union_dar))*100,2)} % ({len(unique_dar)})")
    print(f"Proportion of unique dar of {mol} in {reps} for {dose} dose comparison for ks: {round((len(unique_dar_ks)/len(union_dar))*100,2)} % ({len(unique_dar_ks)})")

    return df


def key_points_of_comparaison(data_control,not_in_trmt,df_trmt_t,data_trmt):
    react_count = dict()

    # number of predicted active reactions 
    df_trmt_t.drop(index=df_trmt_t.index[0], axis=0, inplace=True)
    data_trmt_t = utils.df_to_dict(df_trmt_t)
    for k,v in data_trmt_t.items():
        set_v = set(v)
        filtered_v = list(set_v)
        if 0.0 in filtered_v:
            index = filtered_v.index(0.0)
            filtered_v.pop(index)
        if filtered_v: react_count[k]=len(filtered_v)
    return react_count


def compute_two_sample_KS(ctrl,trmt,k,dose, stat_dar_file,tag=str()):
    res = scipy.stats.ks_2samp(ctrl,trmt,mode = 'auto',alternative="two-sided") ## two_sided: the null hypothesis is that ctrl and trmt are identical (same continuous distribution)
    if res.pvalue < 0.05: # if over 0.05, the null hypothesis is valid, i.e. no difference
        if tag == "both": conditions = f"{dose} dose treatment: \u2705 \t control treatment: \u2705\t"
        elif tag == "control" : conditions = f"{dose} dose treatment: \u274c \t control treatment: \u2705\t"
        else: conditions = f"{dose} dose treatment: \u2705 \t control treatment: \u274c\t"

        stat_dar_file.write(str(k)+'\t'f"{conditions}\t{res.pvalue} \n")