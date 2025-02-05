#
#
# Script to identify the sampling coverage by riptide 
# 19.11.2024
# Made by Maxime Lecomte
#
#

import calcul
import matplotlib.pyplot as plt
from matplotlib_venn import venn2,venn3
import pandas as pd

import utils
import cobra
import re

MOLECULES = ["amiodarone","valproic acid"]
REPLICATES = ["replicate_0","replicate_1"]
DOSE_LEVEL = ["Low","Middle","High"]


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
            data_ctrl,_,df_ctrl = utils.read_sample_file(path_control)

            for dose in DOSE_LEVEL:
                dars[mol][reps][dose] = set()
                dars_ks[mol][reps][dose] = set()
                replicate_analysis[mol][reps][dose] = set()
                replicate_analysis_ks[mol][reps][dose] = set()
                imat_riptide[mol][reps][dose] = {}


                path_trmt = 'results/riptide/recon2.2/maxfit/'+mol+'/1000/24_'+dose+'/'+reps
                data_trmt,df_trmt_t,df_trmt = utils.read_sample_file(path_trmt)
                dar_file,stats_dar,replicate_analysis,replicate_analysis_ks = calcul.compute_R2(data_ctrl,data_trmt,df_trmt_t,replicate_analysis,replicate_analysis_ks,mol=mol,dose=dose,reps=reps,dar_path=path_dar)
                dars = utils.read_dar_file(mol,reps,dose,dars,dar_file)
                dars_ks = utils.read_dar_file(mol,reps,dose,dars_ks,stats_dar)
    
    # # quid replicats ?
    # replicat_intra_molecule(replicate_analysis,tag='r2')
    # replicat_intra_molecule(replicate_analysis_ks,tag='ks')

    # # ## inter molecules & DAR identification methods on RIPTiDe
    # inter_molecules_and_dar_with_riptide(dars,reps,dose,dars_ks)
    
    # # ## intra molecules & inter DAR identification methods on RIPTiDe
    # intra_molecules_and_inter_dar_with_riptide(dars,reps,dose,dars_ks)

    # # inter_molecules & DAR identification methods on iMAT
    # inter_molecules_and_dar_with_imat()
    
    # # intra molecules & inter DAR identification methods on iMAT
    # intra_molecules_and_inter_dar_with_imat()

    # ## intra_molecules & inter DAR identification methods on iMAT and riptide
    intra_molecules_inter_dar_and_context_methods(dars,dars_ks)

    # ## inter_molecules & DAR identification methods on iMAT and riptide
    inter_molecules_inter_dar_and_context_methods(dars,dars_ks)
    
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
                    plt.savefig("results/riptide/recon2.2/maxfit/intra_molecules_intra_DAR/replicates_"+mol+"_"+dose[0]+"_"+tag+".png")

                    all_analysis[mol]["union"] = utils.remove_nan_from_list(all_analysis[mol]["union"])
                    all_analysis[mol]["intersection"] = utils.remove_nan_from_list(all_analysis[mol]["intersection"])
                    all_analysis[mol]["unique_rep_1"] = utils.remove_nan_from_list(all_analysis[mol]["unique_rep_1"])
                    all_analysis[mol]["unique_rep_2"] = utils.remove_nan_from_list(all_analysis[mol]["unique_rep_2"])


                    all_analysis_df = pd.DataFrame(all_analysis)
                    all_analysis_df.to_csv("results/riptide/recon2.2/maxfit/intra_molecules_intra_DAR/replicates_"+mol+"_"+dose[0]+"_"+tag+".tsv",sep='\t')
            break


def inter_molecules_and_dar_with_riptide(dars,reps,dose,dars_ks):
    df_r2 = pd.DataFrame()
    df_ks = pd.DataFrame()

    for reps in REPLICATES:
        for dose in DOSE_LEVEL:
            df_r2 = calcul.compute_dar_specificity_ratio_inter_molecules(dars,reps,dose,df_r2,tag='r2')
            df_ks = calcul.compute_dar_specificity_ratio_inter_molecules(dars_ks,reps,dose,df_ks,tag="ks")

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
    df_comparison.to_csv("results/riptide/recon2.2/maxfit/inter_molecules/compare_r2_Ks.tsv",sep='\t')

    df_r2.to_csv("results/riptide/recon2.2/maxfit/inter_molecules/r2_dar.tsv",sep='\t')
    df_ks.to_csv("results/riptide/recon2.2/maxfit/inter_molecules/ks_dar.tsv",sep='\t')
    
    model = cobra.io.read_sbml_model("data/metabolic_networks/recon2.2.xml")

    for col in df_ks:
        if 'union_a_v_replicate_' in col:
            list_reaction_ks = list(df_ks[col].values)
            annot_df_ks = generate_annotation_table(list_reaction_ks,"results/riptide/recon2.2/maxfit/inter_molecules/ks_dar_annotated_"+col+".tsv",model)

            list_reaction_r2 = list(df_r2[col].values)
            annot_df_r2 = generate_annotation_table(list_reaction_r2,"results/riptide/recon2.2/maxfit/inter_molecules/r2_dar_annotated_"+col+".tsv",model)

            plt.figure()
            venn2([set(annot_df_ks["Pathway in model"].values),set(annot_df_r2["Pathway in model"].values)], set_labels = ("KS","R2"))
            plt.savefig("results/riptide/recon2.2/maxfit/inter_molecules/dar_annotated_"+col+".png", bbox_inches='tight')


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

    df_r2 = calcul.compute_dar_specificity_ratio_inter_molecules_iMAT(dars,df_r2,tag='r2')
    df_ks = calcul.compute_dar_specificity_ratio_inter_molecules_iMAT(dars,df_ks,tag="ks")

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

        result.append({'Colonnes': df_ks.columns[i],'Intersection': list(intersection), 'length intersection': len(intersection), 'difference_r2_ks': list(diff_r2_ks), 'length difference r2 ks': len(diff_r2_ks), 'difference_ks_r2': list(diff_ks_r2), 'length difference ks r2': len(diff_ks_r2)})

    df_comparison = pd.DataFrame(result)
    df_comparison.to_csv("results/iMAT/recon2.2/inter_molecules/compare_r2_Ks.tsv",sep='\t')
    
    model = cobra.io.read_sbml_model("data/metabolic_networks/recon2.2.xml")
    for col in df_ks:
        if 'union_a_v_' in col:
            list_reaction_ks = list(df_ks[col].values)
            annot_df_ks = generate_annotation_table(list_reaction_ks,"results/iMAT/recon2.2/inter_molecules/ks_dar_annotated_"+col+".tsv",model)

            list_reaction_r2 = list(df_r2[col].values)
            annot_df_r2 = generate_annotation_table(list_reaction_r2,"results/iMAT/recon2.2/inter_molecules/r2_dar_annotated_"+col+".tsv",model)

            plt.figure()
            venn2([set(annot_df_ks["Pathway in model"].values),set(annot_df_r2["Pathway in model"].values)], set_labels = ("KS","R2"))
            plt.savefig("results/iMAT/recon2.2/inter_molecules/dar_annotated_"+col+".png", bbox_inches='tight')

    df_r2.to_csv("results/iMAT/recon2.2/inter_molecules/r2_dar.tsv",sep='\t')
    df_ks.to_csv("results/iMAT/recon2.2/inter_molecules/ks_dar.tsv",sep='\t')


def intra_molecules_and_inter_dar_with_riptide(dars,reps,dose,dars_ks):
    df_mol = pd.DataFrame()

    for mol in MOLECULES:
        for reps in REPLICATES:
            for dose in DOSE_LEVEL:
                df_mol = calcul.compute_dar_specificity_ratio_intra_molecules(dars,dars_ks,mol,reps,dose,df_mol)

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

        diff_ks_r2 = utils.remove_nan_from_list(diff_ks_r2)
        diff_r2_ks = utils.remove_nan_from_list(diff_r2_ks)
        intersection = utils.remove_nan_from_list(intersection)

        plt.figure()
        venn2([set(df_r2['dar']),set(df_ks['dar'])], set_labels = ("R2","KS"))
        plt.savefig("results/iMAT/recon2.2/intra_molecules_inter_DAR/"+mol+"_High.png", bbox_inches='tight')

        result.append({'Intersection': list(intersection), 'length intersection': len(intersection), 'difference_r2_ks': list(diff_r2_ks), 'length difference r2 ks': len(diff_r2_ks), 'difference_ks_r2': list(diff_ks_r2), 'length difference ks r2': len(diff_ks_r2)})

    df_comparison = pd.DataFrame(result)
    df_comparison.to_csv("results/iMAT/recon2.2/intra_molecules_inter_DAR/compare_r2_Ks.tsv",sep='\t')

    df_r2.to_csv("results/iMAT/recon2.2/intra_molecules_inter_DAR/r2_dar.tsv",sep='\t')
    df_ks.to_csv("results/iMAT/recon2.2/intra_molecules_inter_DAR/ks_dar.tsv",sep='\t')

def intra_molecules_inter_dar_and_context_methods(dars,dars_ks):

    model = cobra.io.read_sbml_model("data/metabolic_networks/recon2.2.xml")

    for mol,reps in dars.items():
        imat_ks = pd.read_csv("results/iMAT/recon2.2/DAR/ctrl_High_"+mol+"_ks.txt",header=None)
        imat_r2 = pd.read_csv("results/iMAT/recon2.2/DAR/ctrl_High_"+mol+"_r2.txt",header=None)
        riptide_r2 = pd.DataFrame(reps)
        riptide_ks = pd.DataFrame(dars_ks[mol])

        imat_ks.columns = ['dar']
        imat_r2.columns = ['dar']

        plt.figure()
        venn3([riptide_r2.loc[["High"],["replicate_0"]].values[0][0],riptide_r2.loc[["High"],["replicate_1"]].values[0][0],set(imat_r2['dar'].values)], set_labels = ("RIPTiDe_r2_0","RIPTiDe_r2_1","imat_r2"))
        plt.savefig("results/comparison_between_context_method/High_r2_"+mol+'.png')

        plt.figure()
        venn3([riptide_ks.loc[["High"],["replicate_0"]].values[0][0],riptide_ks.loc[["High"],["replicate_1"]].values[0][0],set(imat_ks['dar'].values)], set_labels = ("RIPTiDe_ks_0","RIPTiDe_ks_1","imat_ks"))
        plt.savefig("results/comparison_between_context_method/High_ks_"+mol+'.png')

        annot_imat_r2 = generate_annotation_table(list(imat_r2['dar'].values),"results/comparison_between_context_method/annot_imat_r2.tsv",model)
        annot_riptide_r2_0 = generate_annotation_table(list(riptide_r2.loc[["High"],["replicate_0"]].values[0][0]),"results/comparison_between_context_method/annot_riptide_r2_0.tsv",model)
        annot_riptide_r2_1 = generate_annotation_table(list(riptide_r2.loc[["High"],["replicate_1"]].values[0][0]),"results/comparison_between_context_method/annot_riptide_r2_1.tsv",model)

        annot_imat_ks = generate_annotation_table(list(imat_ks['dar'].values),"results/comparison_between_context_method/annot_imat_ks.tsv",model)
        annot_riptide_ks_0 = generate_annotation_table(list(riptide_ks.loc[["High"],["replicate_0"]].values[0][0]),"results/comparison_between_context_method/annot_riptide_ks_0.tsv",model)
        annot_riptide_ks_1 = generate_annotation_table(list(riptide_ks.loc[["High"],["replicate_1"]].values[0][0]),"results/comparison_between_context_method/annot_riptide_ks_1.tsv",model)

        plt.figure()
        venn3([set(annot_riptide_r2_0["Pathway in model"].values),set(annot_riptide_r2_1["Pathway in model"].values),set(annot_imat_r2["Pathway in model"].values)], set_labels = ("RIPTiDe_r2_0","RIPTiDe_r2_1","imat_r2"))
        plt.savefig("results/comparison_between_context_method/pathway_High_r2_"+mol+'.png', bbox_inches='tight')

        plt.figure()
        venn3([set(annot_riptide_ks_0["Pathway in model"].values),set(annot_riptide_ks_1["Pathway in model"].values),set(annot_imat_ks["Pathway in model"].values)], set_labels = ("RIPTiDe_ks_0","RIPTiDe_ks_1","imat_ks"))
        plt.savefig("results/comparison_between_context_method/pathway_High_ks_"+mol+'.png', bbox_inches='tight')

def generate_annotation_table(reaction_list,output_path,model):
    annot = {}
    annot["Pathway in model"] = {}
    reaction_list = utils.remove_nan_from_list(list(reaction_list))
    for reaction in reaction_list:
        rid = re.sub(r'^R_','',str(reaction))
        if 'EX_' in rid or 'DM_' in rid or 'sink' in rid:
            annot["Pathway in model"][rid] = "Exchange/demand reaction"
        else:
            r = model.reactions.get_by_id(rid)
            if "array([]," in r.subsystem: #if no subsystem, leave cell empty
                annot["Pathway in model"][rid] = ""
            else:
                annot["Pathway in model"][rid] = r.subsystem

    annot_df = pd.DataFrame(annot)
    annot_df.to_csv(output_path,sep='\t')
    return annot_df

if __name__ == "__main__":
    r2_iteration()