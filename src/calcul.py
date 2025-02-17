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
import os
import sampling_coverage
from scipy.stats import chi2_contingency

def compute_freq_table_df(df,name,method=''):
    if method == 'mana':
        return pd.Series([sum(df.iloc[:, col])/df.shape[0] for col in range(df.shape[1])],name=name)
    elif method == 'riptide':
        vecteur = pd.Series()
        for col_i,col in enumerate(df.columns):
            temp_df = df[[col]].copy()
            temp_df = temp_df[temp_df[col] != 0.]
            vecteur = pd.concat([vecteur,pd.Series(temp_df.shape[0]/df.iloc[:, col_i].shape[0],name=name)])

        return vecteur


def compute_r2_df(r2,df,df_trmt):
    intersection_index = df.index.intersection(df_trmt.index)
    if len(intersection_index) == df.shape[0]: ## check if all reactions identified in controle condition are in the trmt condition
        for col in range(df.shape[1]):
            score = np.square(df.loc[intersection_index,col] - df_trmt.loc[intersection_index,col])
            r2 = pd.concat([r2,score])

        r2 = r2.rename(columns={0:'val'})

    elif len(intersection_index) != df.shape[0] or len(intersection_index) != df_trmt.shape[0]:
        df_intersection = pd.DataFrame()

        diff_df_df_trmt = df.index.difference(df_trmt.index)
        diff_df_trmt_df = df_trmt.index.difference(df.index)

        df_diff_df_df_trmt = pd.DataFrame(data=1, index=list(diff_df_df_trmt),columns=['val'])
        df_diff_df_trmt_df = pd.DataFrame(data=1, index=list(diff_df_trmt_df),columns=['val'])

        for col in range(len(intersection_index)):
            score = np.square(df.loc[intersection_index[col]] - df_trmt.loc[intersection_index[col]])
            df_intersection = pd.concat([df_intersection,score])

        df_intersection = df_intersection.set_axis(list(intersection_index))
        df_intersection = df_intersection.rename(columns={0:'val'})
        r2 = pd.concat([r2,df_diff_df_df_trmt,df_diff_df_trmt_df,df_intersection])

    r2 = r2.drop(r2[r2['val'] <= 0.2].index)
    r2 = r2.loc[(r2!=0).any(axis=1)]

    return r2

def chi2_independance_test(df1,df2):
    common_columns = df1.columns.intersection(df2.columns)
    for col in common_columns:
        contingency_table = pd.crosstab(df1[col], df2[col])
        if not contingency_table.empty:
            print("\nTableau de contingence:")
            print(contingency_table)
            # Effectuer le test du chi-2
            chi2, p, dof, expected = chi2_contingency(contingency_table)
            
            # Afficher les résultats
            print("\nRésultats du test du chi-2:")
            print(f"Chi-2: {chi2:.4f}")
            print(f"Degrés de liberté: {dof}")
            print(f"P-valeur: {p:.4f}")
            print("\nFréquences attendues:")
            print(pd.DataFrame(expected, index=contingency_table.index, columns=contingency_table.columns))
            
            # Interprétation du résultat
            alpha = 0.05
            if p < alpha:
                print("=> Rejet de l'hypothèse nulle : les distributions sont significativement différentes.")
            else:
                print("=> Échec du rejet de l'hypothèse nulle : aucune différence significative n'est détectée.")


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
        with open(dar_file,"w+") as f:
            with open(stats_dar,"w+") as stats_dar_file: 
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
                
                react_count = key_points_of_comparaison(data_control,not_in_trmt,df_trmt_t,data_trmt)
                tot_reactions = len(list(data_trmt.keys()))+len(data_control.keys())

                utils.write_stats_file(trmt_and_ctrl_KS,not_in_ctrl,not_in_trmt,number_of_dar_KS,react_count,stats_dar_file,dose,tot_reactions)
                utils.write_stats_file(trmt_and_ctrl,not_in_ctrl,not_in_trmt,number_of_dar,react_count,f,dose,tot_reactions)

        return dar_file,stats_dar,replicate_analysis,replicate_analysis_ks

def compute_dar_specificity_ratio_inter_molecules(path_dar,tool,reps,dose,df,MOLECULES,model,tag=""):
    if tool == 'iMAT' and tag == 'r2': 
        end_path = "/DAR/ctrl_High_"+reps+'.csv'
    elif tool == 'iMAT' and tag == 'ks':
        end_path="/DAR/"+tag.upper()+"_ctrl_High_"+reps+'.csv'
    elif tool == 'riptide' and tag == 'r2': 
        end_path="/10000/DAR/ctrl_"+dose+"_"+reps+'.csv'
    elif tool == 'riptide' and tag == 'ks':
        end_path="/10000/DAR/"+tag.upper()+"_ctrl_"+dose+"_"+reps+'.csv'

    if os.path.exists(path_dar+MOLECULES[0]+end_path):
        if tag == 'ks':
            mol1 = pd.read_csv(path_dar+MOLECULES[0]+end_path,sep='\t',index_col="reactions")
            mol2 = pd.read_csv(path_dar+MOLECULES[1]+end_path,sep='\t',index_col="reactions")
        else:
            mol1 = pd.read_csv(path_dar+MOLECULES[0]+end_path,sep='\t',index_col="Unnamed: 0")
            mol2 = pd.read_csv(path_dar+MOLECULES[1]+end_path,sep='\t',index_col="Unnamed: 0")


        result = []
        list_reaction_ks = list(mol1.index)
        # list_reaction_ks = list(dars_ks[mol][reps][dose])
        annot_df_ks = sampling_coverage.generate_annotation_table(list_reaction_ks,path_dar+"/inter_molecules/annotation/df_"+MOLECULES[0]+"_annotated_"+reps+"_"+dose+".tsv",model)

        list_reaction_r2 = list(mol2.index)
        # list_reaction_r2 = list(dars[mol][reps][dose])
        annot_df_r2 = sampling_coverage.generate_annotation_table(list_reaction_r2,path_dar+"/inter_molecules/annotation/df_"+MOLECULES[1]+"_annotated_"+reps+"_"+dose+".tsv",model)

        plt.figure()
        venn2([set(annot_df_ks["Pathway in model"].values),set(annot_df_r2["Pathway in model"].values)], set_labels = ("Amiodarone","acide valproic"))
        plt.savefig(path_dar+"/inter_molecules/annotation/df_a_v_annotated_"+reps+"_"+dose+"_"+tag+".png", bbox_inches='tight')
                
        intersection_0 = set(annot_df_r2["Pathway in model"]).intersection(set(annot_df_ks["Pathway in model"]))
        diff_r2_ks_0 = set(annot_df_r2["Pathway in model"]).difference(set(annot_df_ks["Pathway in model"]))
        diff_ks_r2_0 = set(annot_df_ks["Pathway in model"]).difference(set(annot_df_r2["Pathway in model"]))

        diff_ks_r2_0 = utils.remove_nan_from_list(diff_ks_r2_0)
        diff_r2_ks_0 = utils.remove_nan_from_list(diff_r2_ks_0)
        intersection_0 = utils.remove_nan_from_list(intersection_0)


        result.append({'Intersection_0': list(intersection_0), 'length intersection_0': len(intersection_0), 'difference_r2_ks_0': list(diff_r2_ks_0), 'length difference r2 ks': len(diff_r2_ks_0), 'difference_ks_r2_0': list(diff_ks_r2_0), 'length difference ks r2_0': len(diff_ks_r2_0)})

        df_comparison = pd.DataFrame(result)
        df_comparison.to_csv(path_dar+"/inter_molecules/annotation/df_dar_annotated_a_v_"+reps+"_"+dose+".tsv",sep='\t')

        union_dar_df = pd.concat([mol1,mol2])
        intersection_dar_df =  mol1.index.intersection(mol2.index)
        unique_dar_mol1 =  mol1.index.difference(mol2.index)
        unique_dar_mol2 =  mol2.index.difference(mol1.index)

        plt.figure()
        venn2([set(mol1.index),set(mol2.index)], set_labels = ("Amiodarone",'Valproic acid'))
        if tool == 'riptide':
            plt.savefig("results/riptide/recon2.2/maxfit/inter_molecules/df_a_v_"+reps+"_"+dose+'_'+tag+'.png')
        else:
            plt.savefig("results/iMAT/recon2.2/inter_molecules/df_a_v_"+reps+"_"+dose+'_'+tag+'.png')


        dar_a = mol1.shape[0]
        only_dar_a = len(unique_dar_mol1)
        dar_v = mol2.shape[0]
        only_dar_v = len(unique_dar_mol2)

        df[f"union_a_v_{reps}_{dose}"] = pd.Series(union_dar_df.index)
        df[f"intersection_a_v_{reps}_{dose}"] = pd.Series(intersection_dar_df)
        df[f"unique_a_{reps}_{dose}"] = pd.Series(unique_dar_mol1)
        df[f"unique_v_{reps}_{dose}"] = pd.Series(unique_dar_mol2)

        print(f"\nProportion of common dars between amiodarone and valproic acid in {reps} for {dose} dose comparison: {round((len(intersection_dar_df)/len(union_dar_df))*100,2)}% ({len(intersection_dar_df)})")
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


def compute_dar_specificity_ratio_intra_molecules(path_dar,tool,mol,reps,dose,df):
    if tool == 'iMAT': 
        path_dar_r2 = path_dar+mol+"/DAR/ctrl_High_"+reps+'.csv'
        path_dar_ks = path_dar+mol+"/DAR/KS_ctrl_High_"+reps+'.csv'

    elif tool == 'riptide': 
        path_dar_r2 = path_dar+mol+"/10000/DAR/ctrl_"+dose+"_"+reps+'.csv'
        path_dar_ks = path_dar+mol+"/10000/DAR/KS_ctrl_"+dose+"_"+reps+'.csv'

    
    if os.path.exists(path_dar_r2) and os.path.exists(path_dar_ks):
        output = path_dar+"/intra_molecules_inter_DAR/"
        mol_ks = pd.read_csv(path_dar_ks,sep='\t',index_col="reactions")
        mol_r2 = pd.read_csv(path_dar_r2,sep='\t',index_col="Unnamed: 0")


        union_dar_df = pd.concat([mol_r2,mol_ks])
        intersection_dar_df =  mol_r2.index.intersection(mol_ks.index)
        unique_dar_mol_ks =  mol_r2.index.difference(mol_ks.index)
        unique_dar_mol_r2 =  mol_ks.index.difference(mol_r2.index)

        plt.figure()
        venn2([set(mol_r2.index),set(mol_ks.index)],set_labels = ("R2","KS"))
        plt.savefig(output+"df_"+mol+'_'+reps+"_"+dose+'.png', bbox_inches='tight')

        df[f"union_{mol}_{reps}_{dose}"] = pd.Series(list(union_dar_df.index))
        df[f"intersection_{mol}_{reps}_{dose}"] = pd.Series(list(intersection_dar_df))
        df[f"unique_r2_{mol}_{reps}_{dose}"] = pd.Series(list(unique_dar_mol_r2))
        df[f"unique_ks_{mol}_{reps}_{dose}"] = pd.Series(list(unique_dar_mol_ks))

        print(f"\nProportion of common dars for {mol} in {reps} for {dose} dose comparison: {round((len(intersection_dar_df)/len(union_dar_df))*100,2)}% ({len(intersection_dar_df)})")
        print(f"Proportion of unique dar of {mol} in {reps} for {dose} dose comparison for r2: {round((len(unique_dar_mol_r2)/len(union_dar_df))*100,2)} % ({len(unique_dar_mol_r2)})")
        print(f"Proportion of unique dar of {mol} in {reps} for {dose} dose comparison for ks: {round((len(unique_dar_mol_ks)/len(union_dar_df))*100,2)} % ({len(unique_dar_mol_ks)})")

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

def compute_two_sample_KS_from_df(ctrl,trmt,dose,stat_dar_file):
    ks_df = pd.DataFrame([],columns=["reactions",'dose','control',"p-value"])
    with open(stat_dar_file,'w+') as f:
        for i in ctrl.columns.intersection(trmt.columns):
            res = scipy.stats.ks_2samp(list(ctrl[i]),list(trmt[i]),mode = 'auto',alternative="two-sided") ## two_sided: the null hypothesis is that ctrl and trmt are identical (same continuous distribution)

            if res.pvalue < 0.05: # if over 0.05, the null hypothesis is valid, i.e. no difference
                ks_df= pd.concat([ks_df,pd.DataFrame([[str(i),"\u2705","\u2705",res.pvalue]],columns=["reactions","dose","control","p-value"])])

        for i in ctrl.columns.difference(trmt.columns):
            res = scipy.stats.ks_2samp(list(ctrl[i]),np.full(len(list(ctrl[i])),99999),mode = 'auto',alternative="two-sided")

            if res.pvalue < 0.05:
                ks_df= pd.concat([ks_df,pd.DataFrame([[str(i),"\u274c","\u2705",res.pvalue]],columns=["reactions","dose","control","p-value"])])

        for i in trmt.columns.difference(ctrl.columns):
            res = scipy.stats.ks_2samp(np.full(len(list(trmt[i])),99999),list(trmt[i]),mode = 'auto',alternative="two-sided")

            if res.pvalue < 0.05:
                ks_df= pd.concat([ks_df,pd.DataFrame([[str(i),"\u2705","\u274c",res.pvalue]],columns=["reactions","dose","control","p-value"])])

    ks_df.to_csv(stat_dar_file, sep='\t')

def merge_replicates(rep1,rep2):
    df_rep1 = pd.read_csv(rep1)
    df_rep2 = pd.read_csv(rep2)
    return pd.concat([df_rep1,df_rep2])