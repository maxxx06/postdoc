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
import cobra
import re
import time
import datetime
import os
import upsetplot

import logging

logging.getLogger("cobra.io.sbml").setLevel(logging.ERROR)

MOLECULES = ["amiodarone","valproic acid"]
REPLICATES = ["replicate_0","replicate_1"]
DOSE_LEVEL = ["Control","Low","Middle","High"]

plt.rcParams['figure.max_open_warning'] = 200

def r2_iteration():
    start1 = time.time()
    run_tools(tool="riptide")
    print(f'R2 et KS computed from RIPTiDe data: {str(datetime.timedelta(seconds=round(time.time() - start1)))}')
    start = time.time()
    run_tools(tool="mana")
    print(f'R2 et chi2 computed from MANA data: {str(datetime.timedelta(seconds=round(time.time() - start)))}')

    path_dar = "results/riptide/recon2.2/maxfit/"
    path_dar_mana = "results/iMAT/recon2.2/"

    ## quid replicats ?
    start = time.time()
    model = cobra.io.read_sbml_model("data/metabolic_networks/recon2.2.xml")

    print("*******  replicat analysis **********")
    replicat_intra_molecule(path_dar,model,tag='r2',tool='riptide')
    replicat_intra_molecule(path_dar,model,tag='ks',tool='riptide')

    replicat_intra_molecule(path_dar_mana,model,tag='r2',tool='mana')
    replicat_intra_molecule(path_dar_mana,model,tag='chi2',tool='mana')    
    ## inter molecules & DAR identification methods
    print("\n******* inter molecules and dar **********")
    inter_molecules_and_dar(path_dar,model,tool='riptide')
    inter_molecules_and_dar(path_dar_mana,model,tool='mana')
    
    # ## intra molecules & inter DAR identification methods
    print("\n******* intra molecules and inter dar **********")
    intra_molecules_and_inter_dar(path_dar,model,tool='riptide')
    intra_molecules_and_inter_dar(path_dar_mana,model,tool='mana')

    ## intra_molecules & inter DAR identification methods on mana and riptide
    print("\n******* between MANA and RIPTiDe **********")
    intra_molecules_inter_dar_and_context_methods(path_dar,path_dar_mana,model)
    print(f'analysis performed in : {str(datetime.timedelta(seconds=round(time.time() - start)))}')
    print(f'all workflow analysis performed in : {str(datetime.timedelta(seconds=round(time.time() - start1)))}')

def run_tools(tool=str()):
    for mol in MOLECULES:
        if tool == 'riptide':
            path_dar = "results/riptide/recon2.2/maxfit/"+mol+"/10000/"
        elif tool == 'mana':
            path_dar = "results/iMAT/recon2.2/"+mol+"/"
        for dose in DOSE_LEVEL:
            df_merged = pd.DataFrame()
            freq_merged = pd.DataFrame()
            for reps in REPLICATES:
                freq_table = pd.DataFrame()
                path_dose = path_dar+"24_"+dose+"/"+reps+"/"
                if os.path.exists(path_dose):
                    r2 = pd.DataFrame()
                    _,df_dose = utils.read_sample_file(path_dose,tool=tool)
                    df_dose.to_csv(path_dar+"DAR/files/df_"+dose+"_"+reps+".tsv",sep='\t')
                    df_dose["dose"] = dose
                    df_dose["reps"] = reps
                
                    freq_table = pd.concat([freq_table,calcul.compute_freq_table_df(df_dose,dose,method=tool,doses=dose,rep=reps)])
                    freq_table = freq_table.set_axis(list(df_dose.columns)[:-2])
                    freq_table.to_csv(path_dar+"DAR/files/freq_"+dose+"_"+reps+".tsv",sep='\t')
                    df_merged=pd.concat([df_merged,df_dose],ignore_index=True)

                    if dose != "Control":
                        df_ctrl = pd.read_csv(path_dar+"DAR/files/df_Control_"+reps+".tsv",sep='\t')
                        df_freq_ctrl = pd.read_csv(path_dar+"DAR/files/freq_Control_"+reps+".tsv",index_col='Unnamed: 0',sep='\t')
                        df_trmt = pd.read_csv(path_dar+"DAR/files/df_"+dose+"_"+reps+".tsv",sep='\t')
                        df_freq_trmt = pd.read_csv(path_dar+"DAR/files/freq_"+dose+"_"+reps+".tsv",index_col='Unnamed: 0',sep='\t')[:-2]

                        r2 = calcul.compute_r2_df(r2,df_freq_ctrl,df_freq_trmt)
                        r2.to_csv(path_dar+"DAR/files/ctrl_"+dose+"_"+reps+".tsv",sep='\t')

                        if tool == 'mana': 
                            stats_dar = path_dar+"DAR/files/chi2_ctrl_"+dose+"_"+reps+".tsv"
                            calcul.chi2_independance_test(df_ctrl,df_trmt,stats_dar)
                        elif tool == 'riptide':
                            stats_dar = path_dar+"DAR/files/ks_ctrl_"+dose+"_"+reps+".tsv"
                            calcul.compute_two_sample_KS_from_df(df_ctrl,df_trmt,stats_dar)

            if not df_merged.empty:
                if dose in df_merged["dose"].values:
                    df_merged_tmp = df_merged[df_merged['dose'] == dose]
                    freq_merged = pd.concat([freq_merged,calcul.compute_freq_table_df(df_merged_tmp,"dose_merged",method=tool,doses=dose)])
                    freq_merged = freq_merged.set_axis(list(df_merged.columns)[:-2])
                    freq_merged.to_csv(path_dar+"DAR/files/freq_"+dose+".tsv",sep='\t')

def replicat_intra_molecule(dar_path,model,tag='',tool=''):
    analysis = pd.DataFrame(columns=['union','intersection','unique_rep_1','unique_rep_2'])
    for mol in MOLECULES:
        if tool == 'mana' and tag == 'r2': 
            path_dar=dar_path+mol+"/DAR/files/"
        elif tool == 'mana' and tag == 'chi2':
            path_dar=dar_path+mol+"/DAR/files/"+tag+'_'

        elif tool == 'riptide' and tag == 'r2': 
            path_dar=dar_path+mol+"/10000/DAR/files/"

        elif tool == 'riptide' and tag == 'ks':
            path_dar=dar_path+mol+"/10000/DAR/files/"+tag+'_'

        for dose in DOSE_LEVEL:
            if os.path.exists(path_dar+"ctrl_"+dose+"_"+REPLICATES[0]+".tsv"):
                if tag == 'chi2':
                    rep0 = pd.read_csv(path_dar+"ctrl_"+dose+"_"+REPLICATES[0]+".tsv",sep='\t')["reactions"]
                    rep1 = pd.read_csv(path_dar+"ctrl_"+dose+"_"+REPLICATES[1]+".tsv",sep='\t')["reactions"]
                elif tag == 'ks':
                    rep0 = pd.read_csv(path_dar+"ctrl_"+dose+"_"+REPLICATES[0]+".tsv",sep='\t')["reactions"]
                    rep1 = pd.read_csv(path_dar+"ctrl_"+dose+"_"+REPLICATES[1]+".tsv",sep='\t')["reactions"]
                else:
                    rep0 = pd.read_csv(path_dar+"ctrl_"+dose+"_"+REPLICATES[0]+".tsv",sep='\t')["Unnamed: 0"]
                    rep1 = pd.read_csv(path_dar+"ctrl_"+dose+"_"+REPLICATES[1]+".tsv",sep='\t')["Unnamed: 0"]

                analysis["intersection"] = pd.Series(list(set(rep0).intersection(set(rep1))))
                analysis["unique_rep_1"] = pd.Series(list(set(rep0).difference(set(rep1))))
                analysis["unique_rep_2"] = pd.Series(list(set(rep0).difference(set(rep1))))
                analysis["union"] = pd.concat([rep0,rep1],ignore_index=True)

                plt.figure()
                venn2([set(rep0),set(rep1)], set_labels = (REPLICATES[0],REPLICATES[1]))    
                if tool == 'riptide':
                    plt.savefig("results/riptide/recon2.2/maxfit/intra_molecules_intra_DAR/images/df_replicates_"+mol+"_"+dose+"_"+tag+".png") 
                    analysis.to_csv("results/riptide/recon2.2/maxfit/intra_molecules_intra_DAR/files/df_replicates_"+mol+"_"+dose+"_"+tag+".tsv",sep='\t')
                else:
                    plt.savefig("results/iMAT/recon2.2/intra_molecules_intra_DAR/images/df_replicates_"+mol+"_"+dose+"_"+tag+".png") 
                    analysis.to_csv("results/iMAT/recon2.2/intra_molecules_intra_DAR/files/df_replicates_"+mol+"_"+dose+"_"+tag+".tsv",sep='\t')

                list_reaction_ks = list(set(rep0.values))
                list_reaction_r2 = list(set(rep1.values))
                
                if tool == 'mana':
                    annot_df_ks = generate_annotation_table(list_reaction_ks,"results/iMAT/recon2.2/intra_molecules_intra_DAR/annotation/files/df_annot_replicates_"+mol+"_"+dose+"_"+tag+".tsv",model)
                    annot_df_r2 = generate_annotation_table(list_reaction_r2,"results/iMAT/recon2.2/intra_molecules_intra_DAR/annotation/files/df_annot_replicates_"+mol+"_"+dose+"_"+tag+".tsv",model)

                    plt.figure()
                    venn2([set(annot_df_ks["Pathway in model"].values),set(annot_df_r2["Pathway in model"].values)], set_labels = (REPLICATES[0],REPLICATES[1]))
                    plt.savefig("results/iMAT/recon2.2/intra_molecules_intra_DAR/annotation/images/df_annot_replicates_"+mol+"_"+dose+"_"+tag+".png", bbox_inches='tight')
                    
                else:
                    annot_df_ks = generate_annotation_table(list_reaction_ks,"results/riptide/recon2.2/maxfit/intra_molecules_intra_DAR/annotation/files/df_annot_replicates_"+mol+"_"+dose+"_"+tag+".tsv",model)
                    annot_df_r2 = generate_annotation_table(list_reaction_r2,"results/riptide/recon2.2/maxfit/intra_molecules_intra_DAR/annotation/files/df_annot_replicates_"+mol+"_"+dose+"_"+tag+".tsv",model)

                    plt.figure()
                    venn2([set(annot_df_ks["Pathway in model"].values),set(annot_df_r2["Pathway in model"].values)], set_labels = (REPLICATES[0],REPLICATES[1]))
                    plt.savefig("results/riptide/recon2.2/maxfit/intra_molecules_intra_DAR/annotation/images/df_annot_replicates_"+mol+"_"+dose+"_"+tag+".png", bbox_inches='tight')

def inter_molecules_and_dar(path_dar,model,tool=''):
    df_r2 = pd.DataFrame()
    df_ks = pd.DataFrame()

    for reps in REPLICATES:
        if tool == 'mana':
            df_r2 = calcul.compute_dar_specificity_ratio_inter_molecules(path_dar,tool,reps,"High",df_r2,MOLECULES,model,tag='r2')
            df_ks = calcul.compute_dar_specificity_ratio_inter_molecules(path_dar,tool,reps,"High",df_ks,MOLECULES,model,tag='chi2')
        else:
            for dose in DOSE_LEVEL:
                if dose != "Control":
                    df_r2 = calcul.compute_dar_specificity_ratio_inter_molecules(path_dar,tool,reps,dose,df_r2,MOLECULES,model,tag='r2')
                    df_ks = calcul.compute_dar_specificity_ratio_inter_molecules(path_dar,tool,reps,dose,df_ks,MOLECULES,model,tag='chi2')


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
    if tool == 'mana':
        df_comparison.to_csv("results/iMAT/recon2.2/inter_molecules/files/df_compare_r2_chi2.tsv",sep='\t')
        df_r2.to_csv("results/iMAT/recon2.2/inter_molecules/files/df_r2_dar.tsv",sep='\t')
        df_ks.to_csv("results/iMAT/recon2.2/inter_molecules/files/df_chi2_dar.tsv",sep='\t')

    else:    
        df_comparison.to_csv("results/riptide/recon2.2/maxfit/inter_molecules/files/df_compare_r2_Ks.tsv",sep='\t')
        df_r2.to_csv("results/riptide/recon2.2/maxfit/inter_molecules/files/df_r2_dar.tsv",sep='\t')
        df_ks.to_csv("results/riptide/recon2.2/maxfit/inter_molecules/files/df_ks_dar.tsv",sep='\t')

    for col in df_ks:
        list_reaction_ks = list(set(df_ks[col].values))
        list_reaction_r2 = list(set(df_r2[col].values))
        
        if tool == 'mana':
            annot_df_ks = generate_annotation_table(list_reaction_ks,"results/iMAT/recon2.2/inter_molecules/annotation/files/df_chi2_dar_annotated_"+col+".tsv",model)
            annot_df_r2 = generate_annotation_table(list_reaction_r2,"results/iMAT/recon2.2/inter_molecules/annotation/files/df_r2_dar_annotated_"+col+".tsv",model)
            plt.figure()
            venn2([set(annot_df_ks["Pathway in model"].values),set(annot_df_r2["Pathway in model"].values)], set_labels = ("Chi-2","R2"))
            plt.savefig("results/iMAT/recon2.2/inter_molecules/annotation/images/df_dar_annotated_"+col+".png", bbox_inches='tight')
            
        else:
            annot_df_ks = generate_annotation_table(list_reaction_ks,"results/riptide/recon2.2/maxfit/inter_molecules/annotation/files/df_ks_dar_annotated_"+col+".tsv",model)
            annot_df_r2 = generate_annotation_table(list_reaction_r2,"results/riptide/recon2.2/maxfit/inter_molecules/annotation/files/df_r2_dar_annotated_"+col+".tsv",model)

            plt.figure()
            venn2([set(annot_df_ks["Pathway in model"].values),set(annot_df_r2["Pathway in model"].values)], set_labels = ("KS","R2"))
            plt.savefig("results/riptide/recon2.2/maxfit/inter_molecules/annotation/images/df_dar_annotated_"+col+".png", bbox_inches='tight')

def intra_molecules_and_inter_dar(path_dar,model,tool=''):
    df_mol = pd.DataFrame()
    for mol in MOLECULES:
        for reps in REPLICATES:
            if tool == 'riptide':
                tag='ks'
                for dose in DOSE_LEVEL:
                    if dose != 'Control':
                        df_mol = calcul.compute_dar_specificity_ratio_intra_molecules(path_dar,tool,mol,reps,dose,df_mol)
            elif tool == 'mana':
                tag = 'chi2'
                df_mol = calcul.compute_dar_specificity_ratio_intra_molecules(path_dar,tool,mol,reps,'High',df_mol)

    annotation_intra_molecules_and_inter_dar(path_dar,model,tool=tool)
            

    df_mol = df_mol.apply(lambda col: col.sort_values().reset_index(drop=True))
    df_mol.to_csv(path_dar+"/intra_molecules_inter_DAR/files/df_r2_"+tag+"_dar.tsv",sep='\t')

def intra_molecules_inter_dar_and_context_methods(path_dar,path_dar_mana,model):
    from venny4py.venny4py import venny4py

    result = []
    for mol in MOLECULES:
        for reps in REPLICATES:
            imat_ks = pd.read_csv(path_dar_mana+mol+"/DAR/files/chi2_ctrl_High_"+reps+".tsv",sep='\t',index_col=["reactions"])
            imat_r2 = pd.read_csv(path_dar_mana+mol+"/DAR/files/ctrl_High_"+reps+".tsv",sep='\t',index_col=["Unnamed: 0"])
            riptide_ks = pd.read_csv(path_dar+mol+"/10000/DAR/files/ks_ctrl_High_"+reps+'.tsv',sep='\t',index_col=["reactions"])
            riptide_r2 = pd.read_csv(path_dar+mol+"/10000/DAR/files/ctrl_High_"+reps+'.tsv',sep='\t',index_col=["Unnamed: 0"])
            
            sets={f"mana_chi2_{reps.split('_')[-1]}" : set(imat_ks.index),
                  f"mana_r2_{reps.split('_')[-1]}" : set(imat_r2.index),
                  f"riptide_ks_{reps.split('_')[-1]}" : set(riptide_ks.index),
                  f"riptide_r2_{reps.split('_')[-1]}" : set(riptide_r2.index)
                  }
            set_names = [f"mana_chi2_{reps.split('_')[-1]}",f"mana_r2_{reps.split('_')[-1]}",f"riptide_ks_{reps.split('_')[-1]}",f"riptide_r2_{reps.split('_')[-1]}"]
            all_elems = set(imat_ks.index).union( set(imat_r2.index)).union(set(riptide_ks.index)).union(set(riptide_r2.index))
            df = pd.DataFrame([[e in set(imat_ks.index), e in set(imat_r2.index), e in set(riptide_ks.index), e in set(riptide_r2.index)] for e in all_elems], columns = set_names)
            df_up = df.groupby(set_names).size()
            upsetplot.plot(df_up, orientation='horizontal')
            current_figure = plt.gcf()
            current_figure.savefig(f"results/comparison_between_context_method/images/upset_High_{mol}_{reps}.png")

            # venny4py(sets=sets,out=f"results/comparison_between_context_method/df_High_{mol}_{reps}")

            plt.figure()
            venn2([set(imat_r2.index),set(riptide_r2.index)], set_labels = ("mana_r2","riptide_r2"))
            plt.savefig(f"results/comparison_between_context_method/images/df_High_r2_{mol}_{reps}.png", bbox_inches='tight')
            
            annot_imat_r2 = generate_annotation_table(list(imat_r2.index),f"results/comparison_between_context_method/files/df_annot_imat_r2_{reps}.tsv",model)
            annot_riptide_r2 = generate_annotation_table(list(riptide_r2.index),f"results/comparison_between_context_method/files/df_annot_riptide_r2_{reps}.tsv",model)

            annot_imat_ks = generate_annotation_table(list(imat_ks.index),f"results/comparison_between_context_method/files/annot_imat_chi2_{reps}.tsv",model)
            annot_riptide_ks = generate_annotation_table(list(riptide_ks.index),f"results/comparison_between_context_method/files/annot_riptide_ks_{reps}.tsv",model)

            sets={f"mana_chi2_{reps.split('_')[-1]}" : set(annot_imat_ks["Pathway in model"]),
                  f"mana_r2_{reps.split('_')[-1]}" : set(annot_imat_r2["Pathway in model"]),
                  f"riptide_ks_{reps.split('_')[-1]}" : set(annot_riptide_ks["Pathway in model"]),
                  f"riptide_r2_{reps.split('_')[-1]}" : set(annot_riptide_r2["Pathway in model"])
                  }
            set_names = [f"mana_chi2_{reps.split('_')[-1]}",f"mana_r2_{reps.split('_')[-1]}",f"riptide_ks_{reps.split('_')[-1]}",f"riptide_r2_{reps.split('_')[-1]}"]
            all_elems = set(annot_imat_ks["Pathway in model"]).union( set(annot_imat_r2["Pathway in model"])).union(set(annot_riptide_ks["Pathway in model"])).union(set(annot_riptide_r2["Pathway in model"]))
            df = pd.DataFrame([[e in set(annot_imat_ks["Pathway in model"]), e in set(annot_imat_r2["Pathway in model"]), e in set(annot_riptide_ks["Pathway in model"]), e in set(annot_riptide_r2["Pathway in model"])] for e in all_elems], columns = set_names)
            df_up = df.groupby(set_names).size()
            upsetplot.plot(df_up, orientation='horizontal')
            current_figure = plt.gcf()
            current_figure.savefig(f"results/comparison_between_context_method/images/upset_annot_{mol}_{reps}.png")
            # df_sets = pd.DataFrame(sets)
            # upsetplot.plot(df_sets,orientation = 'horizontal')
            intersection_0 = set(annot_riptide_r2["Pathway in model"]).intersection(set(annot_imat_r2["Pathway in model"]))
            diff_riptide_imat_0 = set(annot_riptide_r2["Pathway in model"]).difference(set(annot_imat_r2["Pathway in model"]))
            diff_imat_riptide_0 = set(annot_imat_r2["Pathway in model"]).difference(set(annot_riptide_r2["Pathway in model"]))

            diff_imat_riptide_0 = utils.remove_nan_from_list(diff_imat_riptide_0)
            diff_riptide_imat_0 = utils.remove_nan_from_list(diff_riptide_imat_0)
            intersection_0 = utils.remove_nan_from_list(intersection_0)


            result.append({'Intersection_0': list(intersection_0), 'length intersection_0': len(intersection_0), 'difference_riptide_imat': list(diff_riptide_imat_0), 'length difference riptide imat': len(diff_riptide_imat_0), 'difference_imat_riptide': list(diff_imat_riptide_0), 'length difference imat riptide': len(diff_imat_riptide_0)})

            df_comparison = pd.DataFrame(result)
            df_comparison.to_csv(f"results/comparison_between_context_method/files/annot_df_High_r2_{mol}_{reps}.tsv",sep='\t')
            
            # venny4py(sets=sets,out=f"results/comparison_between_context_method/df_annot_{mol}_{reps}")
            
            plt.figure()
            venn2([set(annot_imat_r2["Pathway in model"]),set(annot_riptide_r2["Pathway in model"])], set_labels = ("annot_mana_r2","annot_riptide_r2"))
            plt.savefig(f"results/comparison_between_context_method/images/annot_df_High_r2_{mol}_{reps}.png", bbox_inches='tight')
            
def generate_annotation_table(reaction_list,output_path,model):
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

def annotation_intra_molecules_and_inter_dar(path_dar,model,tool=''):
    for mol in MOLECULES:
        for dose in DOSE_LEVEL:
            for reps in REPLICATES:
                if tool == 'mana': 
                    path_dar_r2 = path_dar+mol+"/DAR/files/ctrl_High_"+reps+'.tsv'
                    path_dar_ks = path_dar+mol+"/DAR/files/chi2_ctrl_High_"+reps+'.tsv'
                    
                elif tool == 'riptide': 
                    path_dar_r2 = path_dar+mol+"/10000/DAR/files/ctrl_"+dose+"_"+reps+'.tsv'
                    path_dar_ks = path_dar+mol+"/10000/DAR/files/ks_ctrl_"+dose+"_"+reps+'.tsv'

                if os.path.exists(path_dar_r2) and os.path.exists(path_dar_ks):
                    output = path_dar+"/intra_molecules_inter_DAR/annotation/"
                    if tool == "mana": tag = "chi2"
                    elif tool == "riptide": tag = "ks"
                    mol_ks = pd.read_csv(path_dar_ks,sep='\t',index_col="reactions")
                    mol_r2 = pd.read_csv(path_dar_r2,sep='\t',index_col="Unnamed: 0")

                    result = []
                    list_reaction_ks = list(mol_ks.index)
                    annot_df_ks = generate_annotation_table(list_reaction_ks,output+"files/df_"+tag+"_dar_annotated_"+mol+"_"+reps+"_"+dose+".tsv",model)

                    list_reaction_r2 = list(mol_r2.index)
                    annot_df_r2 = generate_annotation_table(list_reaction_r2,output+"files/df_r2_dar_annotated_"+mol+"_"+reps+"_"+dose+".tsv",model)

                    if tool == 'mana': 
                        labels = ("chi-2",'R2')
                    else:
                        labels = ("KS",'R2')

                    plt.figure()
                    venn2([set(annot_df_ks["Pathway in model"].values),set(annot_df_r2["Pathway in model"].values)], set_labels = labels)
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