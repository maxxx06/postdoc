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

def compute_freq_table_df(df,name,method,doses=str(),rep=str()):
    """
    Computes the frequency table for a given dataframe based on the selected method.
    
    Args:
        df (pd.DataFrame): Input dataframe.
        name (str): Name of the frequency table.
        method (str): Method type ('mana' or 'riptide').
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
        vecteur = pd.Series(name=name)
        if not rep and not doses:
            temp_df = df.copy()
            for col_i,col in enumerate(df.columns):
                temp_df = temp_df[temp_df[col] != 0.]
                res = pd.Series([temp_df.shape[0]/df.iloc[:, col_i].shape[0]],index=[col_i],name=name)
                vecteur = pd.concat([vecteur,res])
        else:
            for col_i,col in enumerate(df.columns[:-2]):
                temp_df = df[[col]].copy()
                temp_df = temp_df[temp_df[col] != 0.]
                res = pd.Series([temp_df.shape[0]/df.iloc[:, col_i].shape[0]],index=[col_i],name=name)
                vecteur = pd.concat([vecteur,res])
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

def compute_dar_specificity_ratio_inter_molecules(path_dar,tool,reps,dose,df,molecules,model,path_dar_files,out_annot_files,out_annot_images,out_images,out_files,tag=""):
    """
    Compute the specificity number of DARs (Differentially Activated Reactions) between two molecules.
    It analyzes and compares the number of DARs of two molecules.
    
    Parameters:
        path_dar (str): Path to the DAR files.
        tool (str): The method used ('mana' or 'riptide').
        reps (str): Replication identifier.
        dose (str): Dose level.
        df (pd.DataFrame): DataFrame to store results.
        molecules (list): List containing the two molecules to compare.
        model (object): Cobra Model used for annotation.
        tag (str, optional): Tag indicating the analysis type ('r2','chi2' or 'ks). Defaults to "".
    
    Returns:
        pd.DataFrame: Updated DataFrame with specificity number of DAR.

    Saves:
        A Venn diagram comparing DARs between Amiodarone and Valproic acid.
        A TSV files containing comparisons.
    """
    if tool == 'mana' and tag == 'r2': 
        end_path = path_dar_files+"ctrl_High_"+reps+'.tsv'
        index="Unnamed: 0"
    elif tool == 'mana' and tag == 'chi2':
        end_path=path_dar_files+tag+"_ctrl_High_"+reps+'.tsv'
        index="reactions"
    elif tool == 'riptide' and tag == 'r2': 
        end_path="/10000/"+path_dar_files+"ctrl_"+dose+"_"+reps+'.tsv'
        index="Unnamed: 0"
    elif tool == 'riptide' and tag == 'ks':
        end_path="/10000/"+path_dar_files+tag+"_ctrl_"+dose+"_"+reps+'.tsv'
        index="reactions"
    else:
        raise Exception("invalid parameters")


    if os.path.exists(path_dar+molecules[0]+"/"+end_path):
        mol1 = pd.read_csv(path_dar+molecules[0]+"/"+end_path,sep='\t',index_col=index)
        mol2 = pd.read_csv(path_dar+molecules[1]+"/"+end_path,sep='\t',index_col=index)

        result_annot = []

        ## annotation part, both molecules 
        list_reaction_mol1 = list(mol1.index)
        annot_df_mol1 = sampling_coverage.generate_annotation_table(list_reaction_mol1,out_annot_files+"df_"+molecules[0]+"_annotated_"+reps+"_"+dose+".tsv",model)

        list_reaction_mol2 = list(mol2.index)
        annot_df_mol2 = sampling_coverage.generate_annotation_table(list_reaction_mol2,out_annot_files+"df_"+molecules[1]+"_annotated_"+reps+"_"+dose+".tsv",model)

        plt.figure()
        venn2([set(annot_df_mol1["Pathway in model"].values),set(annot_df_mol2["Pathway in model"].values)], set_labels = (molecules[0],molecules[1]))
        plt.savefig(out_annot_images+"df_a_v_annotated_"+reps+"_"+dose+"_"+tag+".png", bbox_inches='tight')
                
        intersection_0 = set(annot_df_mol2["Pathway in model"]).intersection(set(annot_df_mol1["Pathway in model"]))
        diff_r2_ks_0 = set(annot_df_mol2["Pathway in model"]).difference(set(annot_df_mol1["Pathway in model"]))
        diff_ks_r2_0 = set(annot_df_mol1["Pathway in model"]).difference(set(annot_df_mol2["Pathway in model"]))

        diff_ks_r2_0 = utils.remove_nan_from_list(diff_ks_r2_0)
        diff_r2_ks_0 = utils.remove_nan_from_list(diff_r2_ks_0)
        intersection_0 = utils.remove_nan_from_list(intersection_0)


        result_annot.append({'Intersection_0': list(intersection_0), 'length intersection_0': len(intersection_0), 'difference_r2_ks_0': list(diff_r2_ks_0), 'length difference r2 ks': len(diff_r2_ks_0), 'difference_ks_r2_0': list(diff_ks_r2_0), 'length difference ks r2_0': len(diff_ks_r2_0)})

        df_comparison = pd.DataFrame(result_annot)
        df_comparison.to_csv(out_annot_files+"df_dar_annotated_a_v_"+reps+"_"+dose+".tsv",sep='\t')

        ## Dar part, both molecules
        result_dar = []
        union_dar_df = pd.concat([mol1,mol2])
        intersection_dar_df =  mol1.index.intersection(mol2.index)
        unique_dar_mol1 =  mol1.index.difference(mol2.index)
        unique_dar_mol2 =  mol2.index.difference(mol1.index)

        result_dar.append({f'Union_{molecules[0]}_{molecules[1]}':union_dar_df,'intersection_dar_df': list(intersection_dar_df), 'length intersection_dar_df': len(intersection_dar_df), 'unique_dar_mol1': list(unique_dar_mol1), 'length unique_dar_mol1': len(unique_dar_mol1), 'unique_dar_mol2': list(unique_dar_mol2), 'length unique_dar_mol2': len(unique_dar_mol2)})

        plt.figure()
        venn2([set(mol1.index),set(mol2.index)], set_labels = (molecules[0],molecules[1]))
        plt.savefig(out_images+"df_a_v_"+reps+"_"+dose+'_'+tag+'.png')

        df_comparison = pd.DataFrame(result_dar)
        df_comparison.to_csv(out_files+"df_dar_a_v_"+reps+"_"+dose+".tsv",sep='\t')


        df[f"union_a_v_{reps}_{dose}"] = pd.Series(union_dar_df.index)
        df[f"intersection_a_v_{reps}_{dose}"] = pd.Series(intersection_dar_df)
        df[f"unique_a_{reps}_{dose}"] = pd.Series(unique_dar_mol1)
        df[f"unique_v_{reps}_{dose}"] = pd.Series(unique_dar_mol2)

        print(f"\nProportion of common dars between amiodarone and valproic acid in {reps} for {dose} dose comparison: {round((len(intersection_dar_df)/len(union_dar_df))*100,2)}% ({len(intersection_dar_df)})")
        print(f"Proportion of unique dar of amiodarone in {reps} for {dose} dose comparison: {round((len(unique_dar_mol1)/mol1.shape[0])*100,2)} % ({len(unique_dar_mol1)})")
        print(f"Proportion of unique dar of valproic acid in {reps} for {dose} dose comparison: {round((len(unique_dar_mol2)/mol2.shape[0])*100,2)}% ({len(unique_dar_mol2)})")

    else:
        print(path_dar+molecules[0]+"/"+end_path+" not exist.")

    return df

def compute_dar_specificity_ratio_intra_molecules(path_dar,tool,mol,reps,dose,df,out_images,tag,path_dar_files):
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
    if tool == "riptide":
        samples= "10000"
    else:
        samples = ""
    path_dar_r2 = path_dar+mol+"/"+samples+"/"+path_dar_files+"ctrl_"+dose+"_"+reps+'.tsv'
    path_dar_ks = path_dar+mol+"/"+samples+"/"+path_dar_files+tag+"_ctrl_"+dose+"_"+reps+'.tsv'

    if os.path.exists(path_dar_r2) and os.path.exists(path_dar_ks):
        mol_ks = pd.read_csv(path_dar_ks,sep='\t',index_col="reactions")
        mol_r2 = pd.read_csv(path_dar_r2,sep='\t',index_col="Unnamed: 0")

        union_dar_df = pd.concat([mol_r2,mol_ks])
        intersection_dar_df =  mol_r2.index.intersection(mol_ks.index)
        unique_dar_mol_ks =  mol_ks.index.difference(mol_r2.index)
        unique_dar_mol_r2 =  mol_r2.index.difference(mol_ks.index)

        plt.figure()
        venn2([set(mol_r2.index),set(mol_ks.index)],set_labels = ("R2",tag))
        plt.savefig(out_images+"df_"+mol+'_'+reps+"_"+dose+'.png', bbox_inches='tight')

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
