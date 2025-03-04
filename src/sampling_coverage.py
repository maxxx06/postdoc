#
#
# Script to identify the sampling coverage by riptide 
# 19.11.2024
# Made by Maxime Lecomte
#
#

import calcul
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import pandas as pd

import utils
import re
import time
import datetime
import os
import upsetplot

import logging

logging.getLogger("cobra.io.sbml").setLevel(logging.ERROR)
plt.rcParams['figure.max_open_warning'] = 200

def r2_iteration(path_samples,tool,tag,model,dose_level,replicates,molecules,path_dar_files):
    """
    Executes a multi-step workflow analysis for molecular data using the specified tool, model, and parameters.

    This function performs the following steps:
    1. Runs the analysis tool on the given data.
    2. Conducts intra-molecule replicate analysis.
    3. Performs inter-molecule and DAR identification.
    4. Conducts intra-molecule and inter-DAR analysis.

    Args:
        path_samples (str): Path to the directory containing sample data.
        tool (str): The name of the tool used for analysis.
        tag (str): A unique identifier for the analysis.
        model (str): The model used for analysis.
        dose_level (str): The dose level applied in the study.
        replicates (int): The number of replicates in the study.
        molecules (str): The type of molecules being analyzed.
        path_dar_files (str): Path to the directory containing DAR-related files.

    Returns:
        None: The function prints progress updates and analysis results.
    """
    print(f"Analysis of {tool} data, for {molecules} at {dose_level} dose and replicate {replicates} begin...")
    print(f"the Tag {tag} is choosen")
    start1 = time.time()
    run_tools(tool,dose_level,replicates,molecules,path_samples,path_dar_files,tag)
    print(f'{tag} computed: {str(datetime.timedelta(seconds=round(time.time() - start1)))}')

    ## quid replicats ?
    print("*******  replicat analysis **********")
    replicat_intra_molecule(path_samples,model,dose_level,replicates,molecules,path_dar_files,tag,tool)

    ## inter molecules & DAR identification methods
    print("\n******* inter molecules and dar **********")
    inter_molecules_and_dar(path_samples,model,dose_level,replicates,molecules,tool,tag,path_dar_files)

    # ## intra molecules & inter DAR identification methods
    print("\n******* intra molecules and inter dar **********")
    intra_molecules_and_inter_dar(path_samples,model,dose_level,replicates,molecules,tool,tag,path_dar_files)

    print(f'all workflow analysis performed in : {str(datetime.timedelta(seconds=round(time.time() - start1)))}')

def run_tools(tool,dose_level,replicates,molecules,path_samples_init,path_dar_files,tag):
    """
    This function processes the input files, computes frequency tables for each dose level and replicate,
    calculates statistical measures like R2 and KS (or chi-squared), executes the analysis of frequency tables and DAR computation for a given tool (either RIPTiDe or MANA). and stores the results in appropriate directories.
    Outputs the results of the computations as TSV files.

    Args:
        tool (str): The name of the tool used for analysis ('riptide' or 'mana').
        dose_level (list): List of dose levels.
        replicates (list): List of replicates.
        molecules (list): List of molecules to be analyzed.
        path_samples_init (str): Initial path to sample data.
        path_dar_files (str): Path for storing DAR-related files.
        tag (str): Unique identifier for the analysis.

    Returns:
        None: The function processes and saves output files based on computations.
    """
    for mol in molecules:
        if tool == 'riptide':
            path_samples = path_samples_init+"/"+mol+"/10000/"
        elif tool == 'mana':
            path_samples = path_samples_init+"/"+mol+"/"

        utils.create_directory_if_not_exists(path_samples+path_dar_files)
        for dose in dose_level:
            df_merged = pd.DataFrame()
            freq_merged = pd.DataFrame()
            for reps in replicates:
                freq_table = pd.DataFrame()
                path_dose = path_samples+"24_"+dose+"/"+reps+"/"
                if os.path.exists(path_dose):
                    r2 = pd.DataFrame()
                    _,df_dose = utils.read_sample_file(path_dose,tool=tool)
                    df_dose.to_csv(path_samples+path_dar_files+"/df_"+dose+"_"+reps+".tsv",sep='\t',mode='w+')
                    df_dose["dose"] = dose
                    df_dose["reps"] = reps
                    freq_table = pd.concat([freq_table,calcul.compute_freq_table_df(df_dose,dose,method=tool,doses=dose,rep=reps)])
                    freq_table = freq_table.set_axis(list(df_dose.columns)[:-2])
                    freq_table.to_csv(path_samples+path_dar_files+"freq_"+dose+"_"+reps+".tsv",sep='\t')
                    df_merged=pd.concat([df_merged,df_dose],ignore_index=True)

                    if dose != "Control":
                        df_ctrl = pd.read_csv(path_samples+path_dar_files+"df_Control_"+reps+".tsv",sep='\t')
                        df_freq_ctrl = pd.read_csv(path_samples+path_dar_files+"freq_Control_"+reps+".tsv",index_col='Unnamed: 0',sep='\t')
                        df_trmt = pd.read_csv(path_samples+path_dar_files+"df_"+dose+"_"+reps+".tsv",sep='\t')
                        df_freq_trmt = pd.read_csv(path_samples+path_dar_files+"freq_"+dose+"_"+reps+".tsv",index_col='Unnamed: 0',sep='\t')[:-2]

                        r2 = calcul.compute_r2_df(r2,df_freq_ctrl,df_freq_trmt)
                        r2.to_csv(path_samples+path_dar_files+"ctrl_"+dose+"_"+reps+".tsv",sep='\t')

                        stats_dar = path_samples+path_dar_files+tag+"_ctrl_"+dose+"_"+reps+".tsv"
                        if tool == 'mana': 
                            calcul.chi2_independance_test(df_ctrl,df_trmt,stats_dar)
                        elif tool == 'riptide':
                            calcul.compute_two_sample_KS_from_df(df_ctrl,df_trmt,stats_dar)

            if not df_merged.empty:
                if dose in df_merged["dose"].values:
                    df_merged_tmp = df_merged[df_merged['dose'] == dose]
                    freq_merged = pd.concat([freq_merged,calcul.compute_freq_table_df(df_merged_tmp,"dose_merged",method=tool,doses=dose)])
                    freq_merged = freq_merged.set_axis(list(df_merged.columns)[:-2])
                    freq_merged.to_csv(path_samples+path_dar_files+"freq_"+dose+".tsv",sep='\t')

def replicat_intra_molecule(dar_path,model,dose_level,replicates,molecules,path_dar_files,tag,tool):
    """
    Performs an intra-molecule replicat analysis for a given tool (RIPTiDe or MANA).
    Computes the intersection, unique reactions for each replicate, and outputs the results.
    
    Also generates and saves Venn diagrams representing the overlap between replicates for each dose level, at DAR and pathway levels.

    Args:
        dar_path (str): The path where DAR results are stored.
        model (cobra.Model): The metabolic network model used for annotation.
        tag (str): The type of analysis ('r2', 'ks', or 'chi2').
        tool (str): The tool used for analysis ('riptide' or 'mana').
    """
    analysis = pd.DataFrame(columns=['union','intersection','unique_rep_1','unique_rep_2'])

    for mol in molecules:
        if tool == 'mana':
            if tag == 'chi2':
                path_dar=dar_path+mol+"/"+path_dar_files+tag+'_'
                index = "reactions"
            else:
                path_dar=dar_path+mol+"/"+path_dar_files
                index = "Unnamed: 0"
                tag='r2'

        elif tool == 'riptide':
            if tag == 'ks':
                path_dar=dar_path+mol+"/10000/"+path_dar_files+tag+'_'
                index = "reactions"
            else:
                path_dar=dar_path+mol+"/10000/"+path_dar_files
                index = "Unnamed: 0"
                tag='r2'

        for dose in dose_level:
            try:
                if os.path.exists(path_dar+"/ctrl_"+dose+"_"+replicates[0]+".tsv"):
                    out_images = dar_path+"/intra_molecules_intra_DAR/images/"
                    out_files = dar_path+"/intra_molecules_intra_DAR/files/"
                    out_annot_images = dar_path+"/intra_molecules_intra_DAR/annotation/images/"
                    out_annot_files = dar_path+"/intra_molecules_intra_DAR/annotation/files/"

                    utils.create_directory_if_not_exists(out_images)
                    utils.create_directory_if_not_exists(out_files)
                    utils.create_directory_if_not_exists(out_annot_images)
                    utils.create_directory_if_not_exists(out_annot_files)

                    rep0 = pd.read_csv(path_dar+"ctrl_"+dose+"_"+replicates[0]+".tsv",sep='\t')[index]
                    rep1 = pd.read_csv(path_dar+"ctrl_"+dose+"_"+replicates[1]+".tsv",sep='\t')[index]

                    analysis["intersection"] = pd.Series(list(set(rep0).intersection(set(rep1))))
                    analysis["unique_rep_1"] = pd.Series(list(set(rep0).difference(set(rep1))))
                    analysis["unique_rep_2"] = pd.Series(list(set(rep0).difference(set(rep1))))
                    analysis["union"] = pd.concat([rep0,rep1],ignore_index=True)

                    plt.figure()
                    venn2([set(rep0),set(rep1)], set_labels = (replicates[0],replicates[1]))    
                    plt.savefig(out_images+"df_replicates_"+mol+"_"+dose+"_"+tag+".png") 

                    analysis.to_csv(out_files+"df_replicates_"+mol+"_"+dose+"_"+tag+".tsv",sep='\t')

                    list_reaction_rep0 = list(set(rep0.values))
                    list_reaction_rep1 = list(set(rep1.values))
                
                    annot_df_rep0 = generate_annotation_table(list_reaction_rep0,out_annot_files+"df_annot_replicates_"+mol+"_"+dose+"_"+replicates[0]+'_'+tag+".tsv",model)
                    annot_df_rep1 = generate_annotation_table(list_reaction_rep1,out_annot_files+"df_annot_replicates_"+mol+"_"+dose+"_"+replicates[1]+'_'+tag+".tsv",model)

                    plt.figure()
                    venn2([set(annot_df_rep0["Pathway in model"].values),set(annot_df_rep1["Pathway in model"].values)], set_labels = (replicates[0],replicates[1]))
                    plt.savefig(out_annot_images+"df_annot_replicates_"+mol+"_"+dose+"_"+tag+".png", bbox_inches='tight')
            
            except IOError:
                raise IOError('File does not exist: %s' % path_dar+"/ctrl_"+dose+"_"+replicates[0]+".tsv")

                

def inter_molecules_and_dar(path_dar,model,dose_level,replicates,molecules,tool,tag,path_dar_files):
    """
    Compares DAR specificity ratios (R2 and Chi2/ KS) across molecules and generates associated plots and annotations.

    Parameters:
    - path_dar (str): The directory path containing the DAR data files.
    - model (object): The model object used for generating annotations.
    - tool (str): The tool to use for computation ('mana' or 'riptide'). 

    This function processes the DAR specificity ratios for different molecules and doses across replicates. 
    It calculates the intersection and differences between the R2 and Chi2/ KS ratios for molecules, then 
    saves results to appropriate directories. It also generates and saves Venn diagrams for visualizing 
    the relationships between reactions and pathways in the model.
    """
    df_r2 = pd.DataFrame()
    df_ks = pd.DataFrame()

    out_files = path_dar+"/inter_molecules/files/"
    out_images = path_dar+"/inter_molecules/images/"
    out_annot_files = path_dar+"/inter_molecules/annotation/files/"
    out_annot_images = path_dar+"/inter_molecules/annotation/images/"

    utils.create_directory_if_not_exists(out_annot_files)
    utils.create_directory_if_not_exists(out_annot_images)
    utils.create_directory_if_not_exists(out_files)
    utils.create_directory_if_not_exists(out_images)


    for reps in replicates:
        if tool == 'mana':
            df_r2 = calcul.compute_dar_specificity_ratio_inter_molecules(path_dar,tool,reps,"High",df_r2,molecules,model,path_dar_files,out_annot_files,out_annot_images,out_images,out_files,tag='r2')
            df_ks = calcul.compute_dar_specificity_ratio_inter_molecules(path_dar,tool,reps,"High",df_ks,molecules,model,path_dar_files,out_annot_files,out_annot_images,out_images,out_files,tag=tag)
        else:
            for dose in dose_level:
                if dose != "Control":
                    df_r2 = calcul.compute_dar_specificity_ratio_inter_molecules(path_dar,tool,reps,dose,df_r2,molecules,model,path_dar_files,out_annot_files,out_annot_images,out_images,out_files,tag='r2')
                    df_ks = calcul.compute_dar_specificity_ratio_inter_molecules(path_dar,tool,reps,dose,df_ks,molecules,model,path_dar_files,out_annot_files,out_annot_images,out_images,out_files,tag=tag)

    df_r2 = df_r2.apply(lambda col: col.sort_values().reset_index(drop=True))
    df_ks = df_ks.apply(lambda col: col.sort_values().reset_index(drop=True))

    result = []
    for i in range(len(df_ks.columns)):

        intersection = set(df_r2[df_ks.columns[i]]).intersection(set(df_ks[df_ks.columns[i]]))
        diff_r2_ks = set(df_r2[df_ks.columns[i]]).difference(set(df_ks[df_ks.columns[i]]))
        diff_ks_r2 = set(df_ks[df_ks.columns[i]]).difference(set(df_r2[df_ks.columns[i]]))

        diff_ks_r2 = utils.remove_nan_from_list(diff_ks_r2)
        diff_r2_ks = utils.remove_nan_from_list(diff_r2_ks)
        intersection = utils.remove_nan_from_list(intersection)

        result.append({'Colonnes': df_ks.columns[i],'Intersection': intersection, 'length intersection': len(intersection), 'difference_r2_ks': diff_r2_ks , 'length difference r2 ks': len(diff_r2_ks), 'difference_ks_r2': diff_ks_r2, 'length difference ks r2': len(diff_ks_r2)})

    df_comparison = pd.DataFrame(result)
    df_comparison.to_csv(out_files+"df_compare_r2_"+tag+".tsv",sep='\t')
    df_r2.to_csv(out_files+"df_r2_dar.tsv",sep='\t')
    df_ks.to_csv(out_files+"df_"+tag+"_dar.tsv",sep='\t')


    for col in df_ks:
        list_reaction_ks = list(set(df_ks[col].values))
        list_reaction_r2 = list(set(df_r2[col].values))
        
        annot_df_ks = generate_annotation_table(list_reaction_ks,out_annot_files+"df_"+tag+"_dar_annotated_"+col+".tsv",model)
        annot_df_r2 = generate_annotation_table(list_reaction_r2,out_annot_files+"df_r2_dar_annotated_"+col+".tsv",model)
        plt.figure()
        venn2([set(annot_df_ks["Pathway in model"].values),set(annot_df_r2["Pathway in model"].values)], set_labels = (tag,"R2"))
        plt.savefig(out_annot_images+"df_dar_annotated_"+col+".png", bbox_inches='tight')
        

def intra_molecules_and_inter_dar(path_dar,model,dose_level,replicates,molecules,tool,tag,path_dar_files):
    """
    Computes intra-molecule DAR specificity ratios and performs annotation for inter-DAR data.

    Parameters:
    - path_dar (str): The directory path containing the DAR data files.
    - model (object): The model object used for generating annotations.
    - tool (str): The tool to use for computation ('mana' or 'riptide').

    This function calculates the DAR specificity ratio for intra-molecule interactions, considering different 
    replicates and doses. It processes data for either 'riptide' or 'mana' tools, applies sorting, and 
    generates results that are saved to the specified directory. The function also calls an annotation function 
    for further data interpretation.
    """
    out_files = path_dar+"intra_molecules_inter_DAR/files/"
    out_images = path_dar+"intra_molecules_inter_DAR/images/"
    utils.create_directory_if_not_exists(out_files)
    utils.create_directory_if_not_exists(out_images)

    df_mol = pd.DataFrame()
    for mol in molecules:
        for reps in replicates:
            for dose in dose_level:
                df_mol = calcul.compute_dar_specificity_ratio_intra_molecules(path_dar,tool,mol,reps,dose,df_mol,out_images,tag,path_dar_files)

    annotation_intra_molecules_and_inter_dar(path_dar,model,molecules,dose_level,replicates,tool,path_dar_files,tag)
            

    df_mol = df_mol.apply(lambda col: col.sort_values().reset_index(drop=True))
    if os.path.exists(out_files):
        df_mol.to_csv(out_files+"df_r2_"+tag+"_dar.tsv",sep='\t')
    else:
        raise Exception("Path not found.")
    return df_mol

def intra_molecules_inter_dar_and_context_methods(path_dar,path_dar_mana,model,molecules,replicates,path_dar_files,output):
    """
    Compares intra-molecule interactions across different context methods (e.g., 'mana' and 'riptide').

    Parameters:
    - path_dar (str): The directory path containing the DAR data files for 'riptide'.
    - path_dar_mana (str): The directory path containing the DAR data files for 'mana'.
    - model (object): The model object used for generating annotations.

    This function compares intra-molecule DAR specificity ratios across different context methods ('mana' and 
    'riptide') for each molecule and replicate. It generates and saves UpSet plots, Venn diagrams, and 
    annotations that highlight intersections and differences between the context methods. The results are saved 
    to respective files, and figures representing the comparisons are saved as images.
    """
    utils.create_directory_if_not_exists(output+"images/")
    utils.create_directory_if_not_exists(output+"files/")
    result = []
    result_annot = []
    for mol in molecules:
        for reps in replicates:
            imat_ks = pd.read_csv(path_dar_mana+mol+"/"+path_dar_files+"chi2_ctrl_High_"+reps+".tsv",sep='\t',index_col=["reactions"])
            imat_r2 = pd.read_csv(path_dar_mana+mol+"/"+path_dar_files+"ctrl_High_"+reps+".tsv",sep='\t',index_col=["Unnamed: 0"])
            riptide_ks = pd.read_csv(path_dar+mol+"/10000/"+path_dar_files+"ks_ctrl_High_"+reps+'.tsv',sep='\t',index_col=["reactions"])
            riptide_r2 = pd.read_csv(path_dar+mol+"/10000/"+path_dar_files+"ctrl_High_"+reps+'.tsv',sep='\t',index_col=["Unnamed: 0"])

            intersection_0 = set(riptide_r2.index).intersection(set(imat_r2.index))
            diff_riptide_imat_0 = set(riptide_r2.index).difference(set(imat_r2.index))
            diff_imat_riptide_0 = set(imat_r2.index).difference(set(riptide_r2.index))

            diff_imat_riptide_0 = utils.remove_nan_from_list(diff_imat_riptide_0)
            diff_riptide_imat_0 = utils.remove_nan_from_list(diff_riptide_imat_0)
            intersection_0 = utils.remove_nan_from_list(intersection_0)


            result.append({'Intersection_0': list(intersection_0), 'length intersection_0': len(intersection_0), 'difference_riptide_imat': list(diff_riptide_imat_0), 'length difference riptide imat': len(diff_riptide_imat_0), 'difference_imat_riptide': list(diff_imat_riptide_0), 'length difference imat riptide': len(diff_imat_riptide_0)})

            df_comparison = pd.DataFrame(result)
            df_comparison.to_csv(output+f"files/df_High_r2_{mol}_{reps}.tsv",sep='\t')

            set_names = [f"mana_r2_{reps.split('_')[-1]}",f"riptide_ks_{reps.split('_')[-1]}",f"riptide_r2_{reps.split('_')[-1]}"]
            all_elems = set(imat_r2.index).union(set(riptide_ks.index)).union(set(riptide_r2.index))
            df = pd.DataFrame([[e in set(imat_r2.index), e in set(riptide_ks.index), e in set(riptide_r2.index)] for e in all_elems], columns = set_names)
            df_up = df.groupby(set_names).size()
            upsetplot.plot(df_up, orientation='horizontal')
            current_figure = plt.gcf()
            current_figure.savefig(output+f"images/upset_High_{mol}_{reps}.png")

            # venny4py(sets=sets,out=f"results/comparison_between_context_method/df_High_{mol}_{reps}")

            plt.figure()
            venn2([set(imat_r2.index),set(riptide_r2.index)], set_labels = ("mana_r2","riptide_r2"))
            plt.savefig(output+f"images/df_High_r2_{mol}_{reps}.png", bbox_inches='tight')
            
            annot_imat_r2 = generate_annotation_table(list(imat_r2.index),output+f"files/df_annot_imat_r2_{reps}.tsv",model)
            annot_riptide_r2 = generate_annotation_table(list(riptide_r2.index),output+f"files/df_annot_riptide_r2_{reps}.tsv",model)

            annot_imat_ks = generate_annotation_table(list(imat_ks.index),output+f"files/annot_imat_chi2_{reps}.tsv",model)
            annot_riptide_ks = generate_annotation_table(list(riptide_ks.index),output+f"files/annot_riptide_ks_{reps}.tsv",model)

            set_names = [f"mana_chi2_{reps.split('_')[-1]}",f"mana_r2_{reps.split('_')[-1]}",f"riptide_ks_{reps.split('_')[-1]}",f"riptide_r2_{reps.split('_')[-1]}"]
            all_elems = set(annot_imat_ks["Pathway in model"]).union( set(annot_imat_r2["Pathway in model"])).union(set(annot_riptide_ks["Pathway in model"])).union(set(annot_riptide_r2["Pathway in model"]))
            df = pd.DataFrame([[e in set(annot_imat_ks["Pathway in model"]), e in set(annot_imat_r2["Pathway in model"]), e in set(annot_riptide_ks["Pathway in model"]), e in set(annot_riptide_r2["Pathway in model"])] for e in all_elems], columns = set_names)
            df_up = df.groupby(set_names).size()
            upsetplot.plot(df_up, orientation='horizontal')
            current_figure = plt.gcf()
            current_figure.savefig(output+f"images/upset_annot_{mol}_{reps}.png")

            counts_annot_riptide_r2=annot_riptide_r2["Pathway in model"].value_counts()
            counts_annot_riptide_ks=annot_riptide_ks["Pathway in model"].value_counts()
            counts_annot_imat_r2=annot_imat_r2["Pathway in model"].value_counts()
            counts_annot_imat_ks=annot_imat_ks["Pathway in model"].value_counts()

            counts_annot_riptide_r2.to_csv(output+f"files/df_count_annot_riptide_High_r2_{mol}_{reps}.tsv",sep='\t')
            counts_annot_riptide_ks.to_csv(output+f"files/df_count_annot_riptide_High_ks_{mol}_{reps}.tsv",sep='\t')
            counts_annot_imat_r2.to_csv(output+f"files/df_count_annot_mana_High_r2_{mol}_{reps}.tsv",sep='\t')
            counts_annot_imat_ks.to_csv(output+f"files/df_count_annot_mana_High_chi2_{mol}_{reps}.tsv",sep='\t')

            annot_riptide_r2 = annot_riptide_r2[~annot_riptide_r2["Pathway in model"].isin(counts_annot_riptide_r2[counts_annot_riptide_r2 < 5].index)]
            annot_riptide_ks = annot_riptide_ks[~annot_riptide_ks["Pathway in model"].isin(counts_annot_riptide_ks[counts_annot_riptide_ks < 5].index)]
            annot_imat_r2 = annot_imat_r2[~annot_imat_r2["Pathway in model"].isin(counts_annot_imat_r2[counts_annot_imat_r2 < 5].index)]
            annot_imat_ks = annot_imat_ks[~annot_imat_ks["Pathway in model"].isin(counts_annot_imat_ks[counts_annot_imat_ks < 5].index)]


            intersection_0 = set(annot_riptide_r2["Pathway in model"]).intersection(set(annot_imat_r2["Pathway in model"]))
            diff_riptide_imat_0 = set(annot_riptide_r2["Pathway in model"]).difference(set(annot_imat_r2["Pathway in model"]))
            diff_imat_riptide_0 = set(annot_imat_r2["Pathway in model"]).difference(set(annot_riptide_r2["Pathway in model"]))

            diff_imat_riptide_0 = utils.remove_nan_from_list(diff_imat_riptide_0)
            diff_riptide_imat_0 = utils.remove_nan_from_list(diff_riptide_imat_0)
            intersection_0 = utils.remove_nan_from_list(intersection_0)


            result_annot.append({'Intersection_0': list(intersection_0), 'length intersection_0': len(intersection_0), 'difference_riptide_imat': list(diff_riptide_imat_0), 'length difference riptide imat': len(diff_riptide_imat_0), 'difference_imat_riptide': list(diff_imat_riptide_0), 'length difference imat riptide': len(diff_imat_riptide_0)})

            df_comparison = pd.DataFrame(result_annot)
            df_comparison.to_csv(output+f"files/annot_df_High_r2_{mol}_{reps}.tsv",sep='\t')
            
            # venny4py(sets=sets,out=f"results/comparison_between_context_method/df_annot_{mol}_{reps}")
            
            plt.figure()
            venn2([set(annot_imat_r2["Pathway in model"]),set(annot_riptide_r2["Pathway in model"])], set_labels = ("annot_mana_r2","annot_riptide_r2"))
            plt.savefig(output+f"images/annot_df_High_r2_{mol}_{reps}.png", bbox_inches='tight')
            
def generate_annotation_table(reaction_list,output_path,model):
    """
    Generates an annotation table for a list of reactions, including pathway information, and saves it to a file.

    Parameters:
    - reaction_list (list): A list of reaction IDs to be annotated.
    - output_path (str): The path where the generated annotation table will be saved.
    - model (object): The model object that contains the reactions and their corresponding subsystems.

    Returns:
    - annot_df (DataFrame): A pandas DataFrame containing the reaction annotations, including the pathways in the model.

    This function processes a list of reactions and generates annotations based on their corresponding subsystems 
    in the model. Reactions that belong to exchange, demand, or sink types are labeled accordingly. For other reactions, 
    the pathway information is retrieved from the model. The resulting annotations are stored in a DataFrame and 
    saved to the specified output path as a tab-separated file.
    """
    annot = {}
    annot["Pathway in model"] = {}
    reaction_list = utils.remove_nan_from_list(list(reaction_list))
    for reaction in reaction_list:
        rid = re.sub(r'^R_','',str(reaction),count=1)
        if 'EX_' in rid or 'DM_' in rid or 'sink' in rid:
            annot["Pathway in model"][rid] = "Exchange/demand reaction"
        else:
            if "group_phosphotase_3" in rid: rid='R_'+rid
            if "ALA-DTDe" in rid: rid=rid.replace('-','_DASH_')

            r = model.reactions.get_by_id(rid)
            if "array([]," in r.subsystem: #if no subsystem, leave cell empty
                annot["Pathway in model"][rid] = ""
            else:
                annot["Pathway in model"][rid] = r.subsystem

    annot_df = pd.DataFrame(annot)
    annot_df.to_csv(output_path,sep='\t')
    return annot_df

def annotation_intra_molecules_and_inter_dar(path_dar,model,molecules,dose_level,replicates,tool,path_dar_files,tag):
    """
    Annotates intra-molecule DAR specificity ratios and generates pathway comparison reports for different context methods.

    Parameters:
    - path_dar (str): The directory path containing the DAR data files.
    - model (object): The model object used for generating annotations.
    - tool (str, optional): The tool to use for computation ('mana' or 'riptide'). Default is ''.

    This function processes DAR specificity ratios for different molecules, doses, and replicates across the specified 
    tool ('mana' or 'riptide'). It generates annotated tables for both R2 and Chi2/KS reactions, creates Venn diagrams 
    comparing the pathways associated with these reactions, and saves the results to files. The intersections and 
    differences in the pathways are also calculated and saved for further analysis.
    """

    if tool == "riptide":
        samples= "10000"
    else:
        samples = ""

    for mol in molecules:
        for dose in dose_level:
            for reps in replicates:
                path_dar_r2 = path_dar+mol+"/"+samples+"/"+path_dar_files+"ctrl_"+dose+"_"+reps+'.tsv'
                path_dar_ks = path_dar+mol+"/"+samples+"/"+path_dar_files+tag+"_ctrl_"+dose+"_"+reps+'.tsv'

                if os.path.exists(path_dar_r2) and os.path.exists(path_dar_ks):
                    output = path_dar+"/intra_molecules_inter_DAR/annotation/"
                    utils.create_directory_if_not_exists(output+'images/')
                    utils.create_directory_if_not_exists(output+'files/')

                    mol_ks = pd.read_csv(path_dar_ks,sep='\t',index_col="reactions")
                    mol_r2 = pd.read_csv(path_dar_r2,sep='\t',index_col="Unnamed: 0")

                    result = []
                    list_reaction_ks = list(mol_ks.index)
                    annot_df_ks = generate_annotation_table(list_reaction_ks,output+"files/df_"+tag+"_dar_annotated_"+mol+"_"+reps+"_"+dose+".tsv",model)

                    list_reaction_r2 = list(mol_r2.index)
                    annot_df_r2 = generate_annotation_table(list_reaction_r2,output+"files/df_r2_dar_annotated_"+mol+"_"+reps+"_"+dose+".tsv",model)

                    plt.figure()
                    venn2([set(annot_df_ks["Pathway in model"].values),set(annot_df_r2["Pathway in model"].values)], set_labels = (tag,'R2'))
                    plt.savefig(output+"images/df_dar_annotated_"+mol+"_"+reps+"_"+dose+".png", bbox_inches='tight')
                            
                    intersection_0 = set(annot_df_r2["Pathway in model"]).intersection(set(annot_df_ks["Pathway in model"]))
                    diff_r2_ks_0 = set(annot_df_r2["Pathway in model"]).difference(set(annot_df_ks["Pathway in model"]))
                    diff_ks_r2_0 = set(annot_df_ks["Pathway in model"]).difference(set(annot_df_r2["Pathway in model"]))

                    diff_ks_r2_0 = utils.remove_nan_from_list(diff_ks_r2_0)
                    diff_r2_ks_0 = utils.remove_nan_from_list(diff_r2_ks_0)
                    intersection_0 = utils.remove_nan_from_list(intersection_0)


                    result.append({'Intersection_0': list(intersection_0), 'length intersection_0': len(intersection_0), 'difference_r2_ks_0': list(diff_r2_ks_0), 'length difference r2 ks': len(diff_r2_ks_0), 'difference_ks_r2_0': list(diff_ks_r2_0), 'length difference ks r2_0': len(diff_ks_r2_0)})

                    df_comparison = pd.DataFrame(result)
                    df_comparison.to_csv(output+"files/df_dar_annotated_"+mol+"_"+reps+"_"+dose+".tsv",sep='\t')


if __name__ == "__main__":
    r2_iteration()