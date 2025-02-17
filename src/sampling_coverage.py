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

import logging

logging.getLogger("cobra.io.sbml").setLevel(logging.ERROR)

MOLECULES = ["amiodarone","valproic acid"]
REPLICATES = ["replicate_0","replicate_1"]
DOSE_LEVEL = ["Low","Middle","High"]

plt.rcParams['figure.max_open_warning'] = 200

def r2_iteration():

    # start = time.time()
    # iteration_riptide()
    # print(f'R2 et ks computed from RIPTiDe data: {str(datetime.timedelta(seconds=round(time.time() - start)))}')
    start = time.time()
    iteration_mana()
    print(f'R2 et ks computed from MANA data: {str(datetime.timedelta(seconds=round(time.time() - start)))}')

    path_dar = "results/riptide/recon2.2/maxfit/"
    path_dar_mana = "results/iMAT/recon2.2/"

    ## quid replicats ?
    start = time.time()
    model = cobra.io.read_sbml_model("data/metabolic_networks/recon2.2.xml")

    replicat_intra_molecule(path_dar,model,tag='r2',tool='riptide')
    replicat_intra_molecule(path_dar,model,tag='ks',tool='riptide')

    replicat_intra_molecule(path_dar_mana,model,tag='r2',tool='iMAT')
    replicat_intra_molecule(path_dar_mana,model,tag='ks',tool='iMAT')    

    ## inter molecules & DAR identification methods
    inter_molecules_and_dar(path_dar,tool='riptide')
    inter_molecules_and_dar(path_dar_mana,tool='iMAT')
    
    # ## intra molecules & inter DAR identification methods
    intra_molecules_and_inter_dar(path_dar,tool='riptide')
    intra_molecules_and_inter_dar(path_dar_mana,tool='iMAT')

    ## intra_molecules & inter DAR identification methods on iMAT and riptide
    intra_molecules_inter_dar_and_context_methods(path_dar,path_dar_mana)
    print(f'All analysis performed in : {str(datetime.timedelta(seconds=round(time.time() - start)))}')


def iteration_riptide():
    for mol in MOLECULES:
        freq_ctrl_merged = pd.DataFrame()
        freq_trmt_merged = pd.DataFrame()
        df_ctrl_merged, df_trmt_merged= pd.DataFrame(),pd.DataFrame()

        path_dar = "results/riptide/recon2.2/maxfit/"+mol+"/10000/"
        for reps in REPLICATES:
            freq_ctrl = pd.DataFrame()

            path_control = path_dar+"24_Control/"+reps
            _,df_ctrl = utils.read_sample_file(path_control)

            freq_ctrl = pd.concat([freq_ctrl,calcul.compute_freq_table_df(df_ctrl,"control",method='riptide')])
            freq_ctrl = freq_ctrl.set_axis(list(df_ctrl.columns))
            freq_ctrl.to_csv(path_dar+"DAR/"+"freq_ctrl_"+reps+".tsv",sep='\t')

            for dose in DOSE_LEVEL:
                r2 = pd.DataFrame()
                freq_trmt = pd.DataFrame()
                
                path_trmt = path_dar+'24_'+dose+'/'+reps
                _,df_trmt = utils.read_sample_file(path_trmt)

                stats_dar = path_dar+"DAR/KS_ctrl_"+dose+"_"+reps+".csv"
                calcul.compute_two_sample_KS_from_df(df_ctrl,df_trmt,dose,stats_dar)

                freq_trmt = pd.concat([freq_trmt,calcul.compute_freq_table_df(df_trmt,dose,method='riptide')])
                freq_trmt = freq_trmt.set_axis(list(df_trmt.columns))
                freq_trmt.to_csv(path_dar+"DAR/"+"freq_"+dose+"_"+reps+".tsv",sep='\t')
                
                r2 = calcul.compute_r2_df(r2,freq_ctrl,freq_trmt)
                r2 = r2.rename(columns={'val':'r2_score'})
                r2.to_csv(path_dar+"DAR/ctrl_"+dose+"_"+reps+".csv",sep='\t')


                df_ctrl_merged=pd.concat([df_ctrl_merged,df_ctrl],ignore_index=True)
                df_trmt_merged=pd.concat([df_trmt_merged,df_trmt],ignore_index=True)

                # dar_file,stats_dar,replicate_analysis,replicate_analysis_ks = calcul.compute_R2(data_ctrl,data_trmt,df_trmt_t,replicate_analysis,replicate_analysis_ks,mol=mol,dose=dose,reps=reps,dar_path=path_dar)
                # dars = utils.read_dar_file(mol,reps,dose,dars,dar_file)
                # dars_ks = utils.read_dar_file(mol,reps,dose,dars_ks,stats_dar)

        freq_ctrl_merged = pd.concat([freq_ctrl_merged,calcul.compute_freq_table_df(df_ctrl_merged,"control_merged",method='mana')])
        freq_trmt_merged = pd.concat([freq_trmt_merged,calcul.compute_freq_table_df(df_trmt_merged,"trmt_merged",method='mana')])
        
        freq_ctrl_merged = freq_ctrl_merged.set_axis(list(df_ctrl_merged.columns))
        freq_trmt_merged = freq_trmt_merged.set_axis(list(df_trmt_merged.columns))

        freq_trmt_merged.to_csv(path_dar+"DAR/"+"freq_high.tsv",sep='\t')
        freq_ctrl_merged.to_csv(path_dar+"DAR/"+"freq_ctrl.tsv",sep='\t')
        
        r2_merged = pd.DataFrame()
        r2_merged = calcul.compute_r2_df(r2_merged,freq_ctrl_merged,freq_trmt_merged)
        r2_merged = r2_merged.rename(columns={0:'r2_score'})
        r2_merged.to_csv(path_dar+"DAR/ctrl_high.csv",sep='\t')

def iteration_mana():
    for mol in MOLECULES:
        freq_ctrl_merged = pd.DataFrame()
        freq_trmt_merged = pd.DataFrame()

        df_ctrl_merged, df_trmt_merged= pd.DataFrame(),pd.DataFrame()

        dar_dir_path = "results/iMAT/recon2.2/"+mol+"/"
        for reps in REPLICATES:
            freq_ctrl = pd.DataFrame()
            freq_trmt = pd.DataFrame()

            path_control = dar_dir_path+"24_Control/"+reps+"/"
            df_t_ctrl,df_ctrl = utils.read_sample_file_imat(path_control)
            dose = "High"
            path_trmt = dar_dir_path+"24_"+dose+"/"+reps+"/"
            df_trmt_t,df_trmt = utils.read_sample_file_imat(path_trmt)

            calcul.chi2_independance_test(df_ctrl,df_trmt)
            exit()
            stats_dar = dar_dir_path+"DAR/KS_ctrl_"+dose+"_"+reps+".csv"
            calcul.compute_two_sample_KS_from_df(df_ctrl,df_trmt,dose,stats_dar) ## super slow with mana


            freq_ctrl = pd.concat([freq_ctrl,calcul.compute_freq_table_df(df_ctrl,"control",method='mana')])
            freq_trmt = pd.concat([freq_trmt,calcul.compute_freq_table_df(df_trmt,dose,method='mana')])

            freq_ctrl = freq_ctrl.set_axis(list(df_ctrl.columns))
            freq_trmt = freq_trmt.set_axis(list(df_trmt.columns))

            freq_trmt.to_csv(dar_dir_path+"DAR/"+"freq_High_"+reps+".tsv",sep='\t')
            freq_ctrl.to_csv(dar_dir_path+"DAR/"+"freq_ctrl_"+reps+".tsv",sep='\t')
            
            r2 = pd.DataFrame()
            r2 = calcul.compute_r2_df(r2,freq_ctrl,freq_trmt)
            r2 = r2.rename(columns={0:'r2_score'})
            r2.to_csv(dar_dir_path+"DAR/ctrl_High_"+reps+".csv",sep='\t')


            # dar_file,stats_dar,replicate_analysis_imat,replicate_analysis_ks_imat = calcul.compute_R2(data_ctrl,data_trmt,df_trmt_t,replicate_analysis_imat,replicate_analysis_ks_imat,mol=mol,dose=dose,reps=reps,dar_path=dar_dir_path)
            # dars_imat = utils.read_dar_file(mol,reps,dose,dars_imat,dar_file)
            # dars_ks_imat = utils.read_dar_file(mol,reps,dose,dars_ks_imat,stats_dar)

            df_ctrl_merged=pd.concat([df_ctrl_merged,df_ctrl],ignore_index=True)
            df_trmt_merged=pd.concat([df_trmt_merged,df_trmt],ignore_index=True)

        freq_ctrl_merged = pd.concat([freq_ctrl_merged,calcul.compute_freq_table_df(df_ctrl_merged,"control_merged",method='mana')])
        freq_trmt_merged = pd.concat([freq_trmt_merged,calcul.compute_freq_table_df(df_trmt_merged,"trmt_merged",method='mana')])
        
        freq_ctrl_merged = freq_ctrl_merged.set_axis(list(df_ctrl_merged.columns))
        freq_trmt_merged = freq_trmt_merged.set_axis(list(df_trmt_merged.columns))

        freq_trmt_merged.to_csv(dar_dir_path+"DAR/"+"freq_High.tsv",sep='\t')
        freq_ctrl_merged.to_csv(dar_dir_path+"DAR/"+"freq_ctrl.tsv",sep='\t')
        
        r2_merged = pd.DataFrame()
        r2_merged = calcul.compute_r2_df(r2_merged,freq_ctrl_merged,freq_trmt_merged)
        r2_merged = r2_merged.rename(columns={0:'r2_score'})
        r2_merged.to_csv(dar_dir_path+"DAR/ctrl_High.csv",sep='\t')

def replicat_intra_molecule(dar_path,model,tag='',tool=''):
    analysis = pd.DataFrame(columns=['union','intersection','unique_rep_1','unique_rep_2'])

    for mol in MOLECULES:
        if tool == 'iMAT' and tag == 'r2': 
            path_dar=dar_path+mol+"/DAR/"
        elif tool == 'iMAT' and tag == 'ks':
            path_dar=dar_path+mol+"/DAR/"+tag.upper()+'_'

        elif tool == 'riptide' and tag == 'r2': 
            path_dar=dar_path+mol+"/10000/DAR/"

        elif tool == 'riptide' and tag == 'ks':
            path_dar=dar_path+mol+"/10000/DAR/"+tag.upper()+'_'

        for dose in DOSE_LEVEL:
            if os.path.exists(path_dar+"ctrl_"+dose+"_"+REPLICATES[0]+".csv"):
                if tag == 'ks':
                    rep0 = pd.read_csv(path_dar+"ctrl_"+dose+"_"+REPLICATES[0]+".csv",sep='\t')["reactions"]
                    rep1 = pd.read_csv(path_dar+"ctrl_"+dose+"_"+REPLICATES[1]+".csv",sep='\t')["reactions"]
                else:
                    rep0 = pd.read_csv(path_dar+"ctrl_"+dose+"_"+REPLICATES[0]+".csv",sep='\t')["Unnamed: 0"]
                    rep1 = pd.read_csv(path_dar+"ctrl_"+dose+"_"+REPLICATES[1]+".csv",sep='\t')["Unnamed: 0"]


                analysis["intersection"] = pd.Series(list(set(rep0).intersection(set(rep1))))
                analysis["unique_rep_1"] = pd.Series(list(set(rep0).difference(set(rep1))))
                analysis["unique_rep_2"] = pd.Series(list(set(rep0).difference(set(rep1))))
                analysis["union"] = pd.concat([rep0,rep1],ignore_index=True)

                plt.figure()
                venn2([set(rep0),set(rep1)], set_labels = (REPLICATES[0],REPLICATES[1]))

                if tool == 'riptide':
                    plt.savefig("results/riptide/recon2.2/maxfit/intra_molecules_intra_DAR/df_replicates_"+mol+"_"+dose+"_"+tag+".png") 
                    analysis.to_csv("results/riptide/recon2.2/maxfit/intra_molecules_intra_DAR/df_replicates_"+mol+"_"+dose+"_"+tag+".tsv",sep='\t')
                else:
                    plt.savefig("results/iMAT/recon2.2/intra_molecules_intra_DAR/df_replicates_"+mol+"_"+dose+"_"+tag+".png") 
                    analysis.to_csv("results/iMAT/recon2.2/intra_molecules_intra_DAR/df_replicates_"+mol+"_"+dose+"_"+tag+".tsv",sep='\t')

                list_reaction_ks = list(set(rep0.values))
                list_reaction_r2 = list(set(rep1.values))
                
                if tool == 'iMAT':
                    annot_df_ks = generate_annotation_table(list_reaction_ks,"results/iMAT/recon2.2/intra_molecules_intra_DAR/annotation/df_annot_replicates_"+mol+"_"+dose+"_"+tag+".tsv",model)
                    annot_df_r2 = generate_annotation_table(list_reaction_r2,"results/iMAT/recon2.2/intra_molecules_intra_DAR/annotation/df_annot_replicates_"+mol+"_"+dose+"_"+tag+".tsv",model)

                    plt.figure()
                    venn2([set(annot_df_ks["Pathway in model"].values),set(annot_df_r2["Pathway in model"].values)], set_labels = (REPLICATES[0],REPLICATES[1]))
                    plt.savefig("results/iMAT/recon2.2/intra_molecules_intra_DAR/annotation/df_annot_replicates_"+mol+"_"+dose+"_"+tag+".png", bbox_inches='tight')
                    
                else:
                    annot_df_ks = generate_annotation_table(list_reaction_ks,"results/riptide/recon2.2/maxfit/intra_molecules_intra_DAR/annotation/df_annot_replicates_"+mol+"_"+dose+"_"+tag+".tsv",model)
                    annot_df_r2 = generate_annotation_table(list_reaction_r2,"results/riptide/recon2.2/maxfit/intra_molecules_intra_DAR/annotation/df_annot_replicates_"+mol+"_"+dose+"_"+tag+".tsv",model)

                    plt.figure()
                    venn2([set(annot_df_ks["Pathway in model"].values),set(annot_df_r2["Pathway in model"].values)], set_labels = (REPLICATES[0],REPLICATES[1]))
                    plt.savefig("results/riptide/recon2.2/maxfit/intra_molecules_intra_DAR/annotation/df_annot_replicates_"+mol+"_"+dose+"_"+tag+".png", bbox_inches='tight')
       

def inter_molecules_and_dar(path_dar,tool=''):
    df_r2 = pd.DataFrame()
    df_ks = pd.DataFrame()

    model = cobra.io.read_sbml_model("data/metabolic_networks/recon2.2.xml")
    for reps in REPLICATES:
        if tool == 'iMAT':
            df_r2 = calcul.compute_dar_specificity_ratio_inter_molecules(path_dar,tool,reps,"High",df_r2,MOLECULES,model,tag='r2')
            df_ks = calcul.compute_dar_specificity_ratio_inter_molecules(path_dar,tool,reps,"High",df_ks,MOLECULES,model,tag='ks')
        else:
            for dose in DOSE_LEVEL:
                df_r2 = calcul.compute_dar_specificity_ratio_inter_molecules(path_dar,tool,reps,dose,df_r2,MOLECULES,model,tag='r2')
                df_ks = calcul.compute_dar_specificity_ratio_inter_molecules(path_dar,tool,reps,dose,df_ks,MOLECULES,model,tag='ks')


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
    if tool == 'iMAT':
        df_comparison.to_csv("results/iMAT/recon2.2/inter_molecules/df_compare_r2_Ks.tsv",sep='\t')
        df_r2.to_csv("results/iMAT/recon2.2/inter_molecules/df_r2_dar.tsv",sep='\t')
        df_ks.to_csv("results/iMAT/recon2.2/inter_molecules/df_ks_dar.tsv",sep='\t')

    else:    
        df_comparison.to_csv("results/riptide/recon2.2/maxfit/inter_molecules/df_compare_r2_Ks.tsv",sep='\t')
        df_r2.to_csv("results/riptide/recon2.2/maxfit/inter_molecules/df_r2_dar.tsv",sep='\t')
        df_ks.to_csv("results/riptide/recon2.2/maxfit/inter_molecules/df_ks_dar.tsv",sep='\t')

    for col in df_ks:
        list_reaction_ks = list(set(df_ks[col].values))
        list_reaction_r2 = list(set(df_r2[col].values))
        
        if tool == 'iMAT':
            annot_df_ks = generate_annotation_table(list_reaction_ks,"results/iMAT/recon2.2/inter_molecules/annotation/df_ks_dar_annotated_"+col+".tsv",model)
            annot_df_r2 = generate_annotation_table(list_reaction_r2,"results/iMAT/recon2.2/inter_molecules/annotation/df_r2_dar_annotated_"+col+".tsv",model)
            plt.figure()
            venn2([set(annot_df_ks["Pathway in model"].values),set(annot_df_r2["Pathway in model"].values)], set_labels = ("KS","R2"))
            plt.savefig("results/iMAT/recon2.2/inter_molecules/annotation/df_dar_annotated_"+col+".png", bbox_inches='tight')
            
        else:
            annot_df_ks = generate_annotation_table(list_reaction_ks,"results/riptide/recon2.2/maxfit/inter_molecules/annotation/df_ks_dar_annotated_"+col+".tsv",model)
            annot_df_r2 = generate_annotation_table(list_reaction_r2,"results/riptide/recon2.2/maxfit/inter_molecules/annotation/df_r2_dar_annotated_"+col+".tsv",model)

            plt.figure()
            venn2([set(annot_df_ks["Pathway in model"].values),set(annot_df_r2["Pathway in model"].values)], set_labels = ("KS","R2"))
            plt.savefig("results/riptide/recon2.2/maxfit/inter_molecules/annotation/df_dar_annotated_"+col+".png", bbox_inches='tight')

def intra_molecules_and_inter_dar(path_dar,tool=''):
    df_mol = pd.DataFrame()
    for mol in MOLECULES:
        for reps in REPLICATES:
            if tool == 'riptide':
                for dose in DOSE_LEVEL:
                    df_mol = calcul.compute_dar_specificity_ratio_intra_molecules(path_dar,tool,mol,reps,dose,df_mol)
            elif tool == 'iMAT':
                df_mol = calcul.compute_dar_specificity_ratio_intra_molecules(path_dar,tool,mol,reps,'High',df_mol)

    model = cobra.io.read_sbml_model("data/metabolic_networks/recon2.2.xml")
    annotation_intra_molecules_and_inter_dar(path_dar,model,tool=tool)
            

    df_mol = df_mol.apply(lambda col: col.sort_values().reset_index(drop=True))
    df_mol.to_csv(path_dar+"/intra_molecules_inter_DAR/df_r2_ks_dar.tsv",sep='\t')

def intra_molecules_inter_dar_and_context_methods(path_dar,path_dar_mana):
    from venny4py.venny4py import venny4py

    model = cobra.io.read_sbml_model("data/metabolic_networks/recon2.2.xml")
    result = []
    for mol in MOLECULES:
        for reps in REPLICATES:
            imat_ks = pd.read_csv(path_dar_mana+mol+"/DAR/KS_ctrl_High_"+reps+".csv",sep='\t',index_col=["reactions"])
            imat_r2 = pd.read_csv(path_dar_mana+mol+"/DAR/ctrl_High_"+reps+".csv",sep='\t',index_col=["Unnamed: 0"])
            riptide_ks = pd.read_csv(path_dar+mol+"/10000/DAR/KS_ctrl_High_"+reps+'.csv',sep='\t',index_col=["reactions"])
            riptide_r2 = pd.read_csv(path_dar+mol+"/10000/DAR/ctrl_High_"+reps+'.csv',sep='\t',index_col=["Unnamed: 0"])
            
            sets={f"mana_ks_{reps.split('_')[-1]}" : set(imat_ks.index),
                  f"mana_r2_{reps.split('_')[-1]}" : set(imat_r2.index),
                  f"riptide_ks_{reps.split('_')[-1]}" : set(riptide_ks.index),
                  f"riptide_r2_{reps.split('_')[-1]}" : set(riptide_r2.index)
                  }
            
            venny4py(sets=sets,out=f"results/comparison_between_context_method/df_High_{mol}_{reps}")

            plt.figure()
            venn2([set(imat_r2.index),set(riptide_r2.index)], set_labels = ("mana_r2","riptide_r2"))
            plt.savefig(f"results/comparison_between_context_method/df_High_r2_{mol}_{reps}.png", bbox_inches='tight')
            
            annot_imat_r2 = generate_annotation_table(list(imat_r2.index),f"results/comparison_between_context_method/df_annot_imat_r2_{reps}.tsv",model)
            annot_riptide_r2 = generate_annotation_table(list(riptide_r2.index),f"results/comparison_between_context_method/df_annot_riptide_r2_{reps}.tsv",model)

            annot_imat_ks = generate_annotation_table(list(imat_ks.index),f"results/comparison_between_context_method/annot_imat_ks_{reps}.tsv",model)
            annot_riptide_ks = generate_annotation_table(list(riptide_ks.index),f"results/comparison_between_context_method/annot_riptide_ks_{reps}.tsv",model)

            sets={f"mana_ks_{reps.split('_')[-1]}" : set(annot_imat_ks["Pathway in model"]),
                  f"mana_r2_{reps.split('_')[-1]}" : set(annot_imat_r2["Pathway in model"]),
                  f"riptide_ks_{reps.split('_')[-1]}" : set(annot_riptide_ks["Pathway in model"]),
                  f"riptide_r2_{reps.split('_')[-1]}" : set(annot_riptide_r2["Pathway in model"])
                  }
            
            intersection_0 = set(annot_riptide_r2["Pathway in model"]).intersection(set(annot_imat_r2["Pathway in model"]))
            diff_riptide_imat_0 = set(annot_riptide_r2["Pathway in model"]).difference(set(annot_imat_r2["Pathway in model"]))
            diff_imat_riptide_0 = set(annot_imat_r2["Pathway in model"]).difference(set(annot_riptide_r2["Pathway in model"]))

            diff_imat_riptide_0 = utils.remove_nan_from_list(diff_imat_riptide_0)
            diff_riptide_imat_0 = utils.remove_nan_from_list(diff_riptide_imat_0)
            intersection_0 = utils.remove_nan_from_list(intersection_0)


            result.append({'Intersection_0': list(intersection_0), 'length intersection_0': len(intersection_0), 'difference_riptide_imat': list(diff_riptide_imat_0), 'length difference riptide imat': len(diff_riptide_imat_0), 'difference_imat_riptide': list(diff_imat_riptide_0), 'length difference imat riptide': len(diff_imat_riptide_0)})

            df_comparison = pd.DataFrame(result)
            df_comparison.to_csv(f"results/comparison_between_context_method/annot_df_High_r2_{mol}_{reps}.tsv",sep='\t')
            
            venny4py(sets=sets,out=f"results/comparison_between_context_method/df_annot_{mol}_{reps}")
            
            plt.figure()
            venn2([set(annot_imat_r2["Pathway in model"]),set(annot_riptide_r2["Pathway in model"])], set_labels = ("annot_mana_r2","annot_riptide_r2"))
            plt.savefig(f"results/comparison_between_context_method/annot_df_High_r2_{mol}_{reps}.png", bbox_inches='tight')
            
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
        for reps in REPLICATES:
            for dose in DOSE_LEVEL:
                if tool == 'iMAT': 
                    path_dar_r2 = path_dar+mol+"/DAR/ctrl_High_"+reps+'.csv'
                    path_dar_ks = path_dar+mol+"/DAR/KS_ctrl_High_"+reps+'.csv'
                    
                elif tool == 'riptide': 
                    path_dar_r2 = path_dar+mol+"/10000/DAR/ctrl_"+dose+"_"+reps+'.csv'
                    path_dar_ks = path_dar+mol+"/10000/DAR/KS_ctrl_"+dose+"_"+reps+'.csv'

                if os.path.exists(path_dar_r2) and os.path.exists(path_dar_ks):
                    output = path_dar+"/intra_molecules_inter_DAR/annotation/"

                    mol_ks = pd.read_csv(path_dar_ks,sep='\t',index_col="reactions")
                    mol_r2 = pd.read_csv(path_dar_r2,sep='\t',index_col="Unnamed: 0")

                    result = []
                    list_reaction_ks = list(mol_ks.index)
                    # list_reaction_ks = list(dars_ks[mol][reps][dose])
                    annot_df_ks = generate_annotation_table(list_reaction_ks,output+"df_ks_dar_annotated_"+mol+"_"+reps+"_"+dose+".tsv",model)

                    list_reaction_r2 = list(mol_r2.index)
                    # list_reaction_r2 = list(dars[mol][reps][dose])
                    annot_df_r2 = generate_annotation_table(list_reaction_r2,output+"df_r2_dar_annotated_"+mol+"_"+reps+"_"+dose+".tsv",model)

                    plt.figure()
                    venn2([set(annot_df_ks["Pathway in model"].values),set(annot_df_r2["Pathway in model"].values)], set_labels = ("KS","R2"))
                    plt.savefig(output+"df_dar_annotated_"+mol+"_"+reps+"_"+dose+".png", bbox_inches='tight')
                            
                    intersection_0 = set(annot_df_r2["Pathway in model"]).intersection(set(annot_df_ks["Pathway in model"]))
                    diff_r2_ks_0 = set(annot_df_r2["Pathway in model"]).difference(set(annot_df_ks["Pathway in model"]))
                    diff_ks_r2_0 = set(annot_df_ks["Pathway in model"]).difference(set(annot_df_r2["Pathway in model"]))

                    diff_ks_r2_0 = utils.remove_nan_from_list(diff_ks_r2_0)
                    diff_r2_ks_0 = utils.remove_nan_from_list(diff_r2_ks_0)
                    intersection_0 = utils.remove_nan_from_list(intersection_0)


                    result.append({'Intersection_0': list(intersection_0), 'length intersection_0': len(intersection_0), 'difference_r2_ks_0': list(diff_r2_ks_0), 'length difference r2 ks': len(diff_r2_ks_0), 'difference_ks_r2_0': list(diff_ks_r2_0), 'length difference ks r2_0': len(diff_ks_r2_0)})

                    df_comparison = pd.DataFrame(result)
                    df_comparison.to_csv(output+"df_dar_annotated_"+mol+"_"+reps+"_"+dose+".tsv",sep='\t')


if __name__ == "__main__":
    r2_iteration()