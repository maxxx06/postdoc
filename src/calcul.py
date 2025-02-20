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
from statsmodels.stats import multitest

def compute_freq_table_df(df,name,method='',doses=str(),rep=str()):
    """
    Computes the frequency table for a given dataframe based on the selected method.
    
    Args:
        df (pd.DataFrame): Input dataframe.
        name (str): Name of the frequency table.
        method (str, optional): Method type ('mana' or 'riptide'). Defaults to ''.
        doses (str, optional): Doses information. Defaults to str().
        rep (str, optional): Replicate information. Defaults to str().
    
    Returns:
        pd.Series: Computed frequency table.
    """
    if method == 'mana':
        if not rep and not doses:
            return pd.Series([sum(df.iloc[:, col])/df.shape[0] for col in range(df.shape[1])],name=name)
        else:
            return pd.Series([sum(df.iloc[:, col])/df.shape[0] for col in range(df.shape[1]-2)],name=name)

    elif method == 'riptide':
        vecteur = pd.Series()
        for col_i,col in enumerate(df.columns[:-2]):
            temp_df = df[[col]].copy()
            temp_df = temp_df[temp_df[col] != 0.]
            vecteur = pd.concat([vecteur,pd.Series(temp_df.shape[0]/df.iloc[:, col_i].shape[0],name=name)])

        return vecteur

def compute_r2_df(r2,df,df_trmt):
    """
    Computes R-squared values between two dataframes.
    
    Args:
        r2 (pd.DataFrame): Dataframe to store R-squared values.
        df (pd.DataFrame): Control condition dataframe.
        df_trmt (pd.DataFrame): Treatment condition dataframe.
    
    Returns:
        pd.DataFrame: Updated R-squared dataframe.
    """
    intersection_index = df.index.intersection(df_trmt.index)
    if len(intersection_index) == df.shape[0]: ## check if all reactions identified in controle condition are in the trmt condition
        score = np.square(df.loc[intersection_index,'0'] - df_trmt.loc[intersection_index,'0'])
        r2 = pd.concat([r2,score])

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

    r2 = r2.rename(columns={0:'val'})
    r2 = r2.drop(r2[r2['val'] <= 0.2].index)
    r2 = r2.loc[(r2!=0).any(axis=1)]

    return r2

def chi2_independance_test(df1,df2,stat_dar_file):
    """
    Performs a chi-square independence test between two dataframes.
    
    Args:
        df1 (pd.DataFrame): First dataframe.
        df2 (pd.DataFrame): Second dataframe.
        stat_dar_file (str): Path to the output file.

    Saves:
        A TSV file containing reactions, p-values, and significance indicators based on the chi-square independence test.
    """
    if df1.shape[0] > df2.shape[0]:
        list_to_add = np.full(df1.shape[0] - df2.shape[0],'Nan')
        df2 = pd.concat([df2,pd.DataFrame(list_to_add)])
    if df1.shape[0] < df2.shape[0]:
        list_to_add = np.full(df2.shape[0] - df1.shape[0],'Nan')
        df1 = pd.concat([df1,pd.DataFrame(list_to_add)])

    with open(stat_dar_file,'w+') as f:
        chi_2_df = pd.DataFrame([],columns=["reactions",'dose','control',"p-value corrected","p-value"])
        for col in df1.columns.intersection(df2.columns):
            contingency_table = pd.crosstab(pd.Categorical(df1[col]), pd.Categorical(df2[col]),dropna=True)
            chi2, p, dof, expected = chi2_contingency(contingency_table)
            _,corrected_pvalue,_,_ = multitest.multipletests(p,alpha=0.01,method="bonferroni")
            
            if corrected_pvalue[0] < 0.01: # if over 0.05, the null hypothesis is valid, i.e. no difference
                chi_2_df= pd.concat([chi_2_df,pd.DataFrame([[str(col),"\u2705","\u2705",corrected_pvalue[0],p]],columns=["reactions",'dose','control',"p-value corrected","p-value"])])

        for col in df1.columns.difference(df2.columns):
            contingency_table = pd.crosstab(pd.Categorical(df1[col]), pd.Categorical(np.full(len(list(df1[col])),'Nan')),dropna=False)
            chi2, p, dof, expected = chi2_contingency(contingency_table)
            _,corrected_pvalue,_,_ = multitest.multipletests(p,alpha=0.01,method="bonferroni")

            if corrected_pvalue[0] < 0.01: # if over 0.05, the null hypothesis is valid, i.e. no difference
                chi_2_df= pd.concat([chi_2_df,pd.DataFrame([[str(col),"\u274c","\u2705",corrected_pvalue[0],p]],columns=["reactions",'dose','control',"p-value corrected","p-value"])])

        for col in df2.columns.difference(df1.columns):
            contingency_table = pd.crosstab(pd.Categorical(np.full(len(list(df2[col])),'Nan')), pd.Categorical(df2[col]),dropna=False)
            chi2, p, dof, expected = chi2_contingency(contingency_table)
            _,corrected_pvalue,_,_ = multitest.multipletests(p,alpha=0.01,method="bonferroni")

            if corrected_pvalue[0] < 0.01: # if over 0.05, the null hypothesis is valid, i.e. no difference
                chi_2_df= pd.concat([chi_2_df,pd.DataFrame([[str(col),"\u2705","\u274c",corrected_pvalue[0],p]],columns=["reactions",'dose','control',"p-value corrected","p-value"])])

        chi_2_df.to_csv(stat_dar_file, sep='\t')

def compute_dar_specificity_ratio_inter_molecules(path_dar,tool,reps,dose,df,MOLECULES,model,tag=""):
    """
    Compute the specificity number of DARs (Differentially Activated Reactions) between two molecules.
    It analyzes and compares the number of DARs of two molecules.
    
    Parameters:
        path_dar (str): Path to the DAR files.
        tool (str): The method used ('mana' or 'riptide').
        reps (str): Replication identifier.
        dose (str): Dose level.
        df (pd.DataFrame): DataFrame to store results.
        MOLECULES (list): List containing the two molecules to compare.
        model (object): Cobra Model used for annotation.
        tag (str, optional): Tag indicating the analysis type ('r2','chi2' or 'ks). Defaults to "".
    
    Returns:
        pd.DataFrame: Updated DataFrame with specificity number of DAR.

    Saves:
        A Venn diagram comparing DARs between Amiodarone and Valproic acid.
        A TSV files containing comparisons.
    """

    if tool == 'mana' and tag == 'r2': 
        end_path = "/DAR/files/ctrl_High_"+reps+'.tsv'
    elif tool == 'mana' and tag == 'chi2':
        end_path="/DAR/files/"+tag+"_ctrl_High_"+reps+'.tsv'
    elif tool == 'riptide' and tag == 'r2': 
        end_path="/10000/DAR/files/ctrl_"+dose+"_"+reps+'.tsv'
    elif tool == 'riptide' and tag == 'ks':
        end_path="/10000/DAR/files/"+tag+"_ctrl_"+dose+"_"+reps+'.tsv'

    if os.path.exists(path_dar+MOLECULES[0]+end_path):
        if tag in ['chi2','ks']:
            mol1 = pd.read_csv(path_dar+MOLECULES[0]+end_path,sep='\t',index_col="reactions")
            mol2 = pd.read_csv(path_dar+MOLECULES[1]+end_path,sep='\t',index_col="reactions")
        else:
            mol1 = pd.read_csv(path_dar+MOLECULES[0]+end_path,sep='\t',index_col="Unnamed: 0")
            mol2 = pd.read_csv(path_dar+MOLECULES[1]+end_path,sep='\t',index_col="Unnamed: 0")


        result = []
        list_reaction_ks = list(mol1.index)
        # list_reaction_ks = list(dars_ks[mol][reps][dose])
        annot_df_ks = sampling_coverage.generate_annotation_table(list_reaction_ks,path_dar+"/inter_molecules/annotation/files/df_"+MOLECULES[0]+"_annotated_"+reps+"_"+dose+".tsv",model)

        list_reaction_r2 = list(mol2.index)
        # list_reaction_r2 = list(dars[mol][reps][dose])
        annot_df_r2 = sampling_coverage.generate_annotation_table(list_reaction_r2,path_dar+"/inter_molecules/annotation/files/df_"+MOLECULES[1]+"_annotated_"+reps+"_"+dose+".tsv",model)

        plt.figure()
        venn2([set(annot_df_ks["Pathway in model"].values),set(annot_df_r2["Pathway in model"].values)], set_labels = ("Amiodarone","acide valproic"))
        plt.savefig(path_dar+"/inter_molecules/annotation/images/df_a_v_annotated_"+reps+"_"+dose+"_"+tag+".png", bbox_inches='tight')
                
        intersection_0 = set(annot_df_r2["Pathway in model"]).intersection(set(annot_df_ks["Pathway in model"]))
        diff_r2_ks_0 = set(annot_df_r2["Pathway in model"]).difference(set(annot_df_ks["Pathway in model"]))
        diff_ks_r2_0 = set(annot_df_ks["Pathway in model"]).difference(set(annot_df_r2["Pathway in model"]))

        diff_ks_r2_0 = utils.remove_nan_from_list(diff_ks_r2_0)
        diff_r2_ks_0 = utils.remove_nan_from_list(diff_r2_ks_0)
        intersection_0 = utils.remove_nan_from_list(intersection_0)


        result.append({'Intersection_0': list(intersection_0), 'length intersection_0': len(intersection_0), 'difference_r2_ks_0': list(diff_r2_ks_0), 'length difference r2 ks': len(diff_r2_ks_0), 'difference_ks_r2_0': list(diff_ks_r2_0), 'length difference ks r2_0': len(diff_ks_r2_0)})

        df_comparison = pd.DataFrame(result)
        df_comparison.to_csv(path_dar+"/inter_molecules/annotation/files/df_dar_annotated_a_v_"+reps+"_"+dose+".tsv",sep='\t')

        union_dar_df = pd.concat([mol1,mol2])
        intersection_dar_df =  mol1.index.intersection(mol2.index)
        unique_dar_mol1 =  mol1.index.difference(mol2.index)
        unique_dar_mol2 =  mol2.index.difference(mol1.index)

        plt.figure()
        venn2([set(mol1.index),set(mol2.index)], set_labels = ("Amiodarone",'Valproic acid'))
        if tool == 'riptide':
            plt.savefig("results/riptide/recon2.2/maxfit/inter_molecules/images/df_a_v_"+reps+"_"+dose+'_'+tag+'.png')
        else:
            plt.savefig("results/iMAT/recon2.2/inter_molecules/images/df_a_v_"+reps+"_"+dose+'_'+tag+'.png')


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

def compute_dar_specificity_ratio_intra_molecules(path_dar,tool,mol,reps,dose,df):
    """
    Computes the specificity ratio of DARs within a molecule for different methods (mana or Riptide).

    Args:
        path_dar (str): Path to the DAR files.
        tool (str): The tool used for DAR computation ('mana' or 'riptide').
        mol (str): The molecule name.
        reps (str): Replicate identifier.
        dose (str): Dose condition.
        df (pd.DataFrame): DataFrame to store the computed DAR statistics.

    Returns:
        pd.DataFrame: Updated DataFrame containing union, intersection, and unique DARs for different methods.

    Saves:
        A Venn diagram comparing DARs obtained from different statistical methods.
    """
    if tool == 'mana': 
        path_dar_r2 = path_dar+mol+"/DAR/files/ctrl_High_"+reps+'.tsv'
        path_dar_ks = path_dar+mol+"/DAR/files/chi2_ctrl_High_"+reps+'.tsv'

    elif tool == 'riptide': 
        path_dar_r2 = path_dar+mol+"/10000/DAR/files/ctrl_"+dose+"_"+reps+'.tsv'
        path_dar_ks = path_dar+mol+"/10000/DAR/files/ks_ctrl_"+dose+"_"+reps+'.tsv'

    
    if os.path.exists(path_dar_r2) and os.path.exists(path_dar_ks):
        output = path_dar+"/intra_molecules_inter_DAR/images/"
        mol_ks = pd.read_csv(path_dar_ks,sep='\t',index_col="reactions")
        mol_r2 = pd.read_csv(path_dar_r2,sep='\t',index_col="Unnamed: 0")


        union_dar_df = pd.concat([mol_r2,mol_ks])
        intersection_dar_df =  mol_r2.index.intersection(mol_ks.index)
        unique_dar_mol_ks =  mol_r2.index.difference(mol_ks.index)
        unique_dar_mol_r2 =  mol_ks.index.difference(mol_r2.index)


        if tool == 'mana': 
            labels = ('R2',"chi-2")
        else:
            labels = ('R2',"KS")
                        
        plt.figure()
        venn2([set(mol_r2.index),set(mol_ks.index)],set_labels = labels)
        plt.savefig(output+"df_"+mol+'_'+reps+"_"+dose+'.png', bbox_inches='tight')

        df[f"union_{mol}_{reps}_{dose}"] = pd.Series(list(union_dar_df.index))
        df[f"intersection_{mol}_{reps}_{dose}"] = pd.Series(list(intersection_dar_df))
        df[f"unique_r2_{mol}_{reps}_{dose}"] = pd.Series(list(unique_dar_mol_r2))
        df[f"unique_ks_{mol}_{reps}_{dose}"] = pd.Series(list(unique_dar_mol_ks))

        print(f"\nProportion of common dars for {mol} in {reps} for {dose} dose comparison: {round((len(intersection_dar_df)/len(union_dar_df))*100,2)}% ({len(intersection_dar_df)})")
        print(f"Proportion of unique dar of {mol} in {reps} for {dose} dose comparison for r2: {round((len(unique_dar_mol_r2)/len(union_dar_df))*100,2)} % ({len(unique_dar_mol_r2)})")
        print(f"Proportion of unique dar of {mol} in {reps} for {dose} dose comparison for ks: {round((len(unique_dar_mol_ks)/len(union_dar_df))*100,2)} % ({len(unique_dar_mol_ks)})")

    return df

def compute_two_sample_KS_from_df(ctrl,trmt,stat_dar_file):
    """
    Performs a two-sample Kolmogorov-Smirnov (KS) test on reaction activity between control and treatment groups.

    Args:
        ctrl (pd.DataFrame): DataFrame containing reaction activity data for the control group.
        trmt (pd.DataFrame): DataFrame containing reaction activity data for the treatment group.
        stat_dar_file (str): Path to save the results as a TSV file.

    Saves:
        A TSV file containing reactions, p-values, and significance indicators based on the KS test.
    """
    ks_df = pd.DataFrame([],columns=["reactions",'dose','control',"p-value corrected","p-value"])
    with open(stat_dar_file,'w+') as f:
        for i in ctrl.columns.intersection(trmt.columns):
            res = scipy.stats.ks_2samp(list(ctrl[i]),list(trmt[i]),mode = 'auto',alternative="two-sided") ## two_sided: the null hypothesis is that ctrl and trmt are identical (same continuous distribution)
            _,corrected_pvalue,_,_ = multitest.multipletests(res.pvalue,alpha=0.01,method="bonferroni")
            if corrected_pvalue[0] < 0.01: # if over 0.05, the null hypothesis is valid, i.e. no difference
                ks_df= pd.concat([ks_df,pd.DataFrame([[str(i),"\u2705","\u2705",corrected_pvalue[0],res.pvalue]],columns=["reactions",'dose','control',"p-value corrected","p-value"])])

        for i in ctrl.columns.difference(trmt.columns):
            res = scipy.stats.ks_2samp(list(ctrl[i]),np.full(len(list(ctrl[i])),99999),mode = 'auto',alternative="two-sided")
            _,corrected_pvalue,_,_ = multitest.multipletests(res.pvalue,alpha=0.01,method="bonferroni")

            if corrected_pvalue[0] < 0.01:
                ks_df= pd.concat([ks_df,pd.DataFrame([[str(i),"\u274c","\u2705",corrected_pvalue[0],res.pvalue]],columns=["reactions",'dose','control',"p-value corrected","p-value"])])

        for i in trmt.columns.difference(ctrl.columns):
            res = scipy.stats.ks_2samp(np.full(len(list(trmt[i])),99999),list(trmt[i]),mode = 'auto',alternative="two-sided")
            _,corrected_pvalue,_,_ = multitest.multipletests(res.pvalue,alpha=0.01,method="bonferroni")

            if corrected_pvalue[0] < 0.01:
                ks_df= pd.concat([ks_df,pd.DataFrame([[str(i),"\u2705","\u274c",corrected_pvalue[0],res.pvalue]],columns=["reactions",'dose','control',"p-value corrected","p-value"])])

    ks_df.to_csv(stat_dar_file, sep='\t')
