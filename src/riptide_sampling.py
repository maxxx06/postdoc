#
#
# Script to use riptide on microarray data 
# 19.11.2024
# Made by Maxime Lecomte
#
#

import riptide
import data_management
import utils
import sampling_coverage
import yaml
import pandas as pd


def load_riptide(file_path):
    with open(file_path, 'r') as file:
        config = yaml.safe_load(file)

    data_path=config["DATA"]["path"]
    data_attributes=config["DATA"]["attributes"]
    data_mapped_gene_path=config["DATA"]["mapped_gene"]
    data_output_transcriptomic_folder=config["DATA"]["output_transcriptomic_folder"]
    data_sacrific_period=config["DATA"]["sacrific_period"]
    data_dose_level=[dose if '/' in config["DATA"]["dose_level"] else config["DATA"]["dose_level"][0] for dose in config["DATA"]["dose_level"].split('/')]
    data_compound=[dose if '/' in config["DATA"]["compound"] else config["DATA"]["compound"][0] for dose in config["DATA"]["compound"].split('/')]
    data_replicates=[dose if '/' in config["DATA"]["replicates"] else config["DATA"]["replicates"][0] for dose in config["DATA"]["replicates"].split('/')]
    tag = config["DATA"]["tag"] if 'tag' in config["DATA"] else ''

    output_riptide = config["RIPTiDe"]["out_riptide"]
    riptide_samples=config["RIPTiDe"]["samples"]
    riptide_maxfit=config["RIPTiDe"]["maxfit"]


    media = config["ANALYSIS"]["media"]
    model = utils.load_model(config["RIPTiDe"]['metabolic_model'])
    if "recon2v2" in config["RIPTiDe"]['metabolic_model']:
        model_version = "recon2.2"
        model.name = 'data/metabolic_networks/recon2v2_biomass_corrected_final.sbml'
    else:
        model_version = "human1"
        model.name = "data/metabolic_networks/Human-GEM.xml"
    
    # model = add_media_constrains(model,pd.read_csv(media,sep='\t'))

    for compound in data_compound:
        for dose_level in data_dose_level:
            output_transcriptomic_path = data_output_transcriptomic_folder+compound+'/'+data_sacrific_period.split(' ')[0]+'_'+dose_level+'/'
            utils.create_directory_if_not_exists(output_transcriptomic_path)
            transcriptomic_data = data_management.data_filter(data_path,data_attributes,output_transcriptomic_path+'transcriptomic_data.csv',sacrifice_period=data_sacrific_period,dose_level=dose_level,compound_name=compound,tag=tag)
            transcriptomic_data_mapped = data_management.map_genes(transcriptomic_data,data_mapped_gene_path,output_transcriptomic_path+'transcriptomic_data_mapped.csv',model_version,tag=tag)

            for replicates in data_replicates:
                transcriptomic_data_mapped_reps = {gene:float(exprs[int(replicates)]) for gene,exprs in transcriptomic_data_mapped.items()}
                out = output_riptide+compound+'/'+str(riptide_samples)+'/'+data_sacrific_period.split(' ')[0]+'_'+dose_level+'/replicate_'+str(replicates)+'/'

                utils.create_directory_if_not_exists(out)
                get_flux_samples(model,transcriptomic_data_mapped_reps,samples=riptide_samples,maxfit=riptide_maxfit,out_riptide=out)

def load_analysis(file_path):
    with open(file_path, 'r') as file:
        config = yaml.safe_load(file)

    tools = [el for el in config["ANALYSIS"]['tool'].split("/") ]
    model = utils.load_model(config["ANALYSIS"]['pathway_model_description'])
    path_samples = [el for el in config["ANALYSIS"]["samples"].split(',')]
    tag = [el for el in config["ANALYSIS"]["tag"].split("/")]
    dose_level = [str(el) for el in config["ANALYSIS"]["dose_level"].split("/") ]
    replicates = [str(el) for el in config["ANALYSIS"]["replicates"].split("/") ]
    molecules = [el for el in config["ANALYSIS"]["molecules"].split("/") ]
    path_dar_files = config["ANALYSIS"]["output_dar_files"]
    model = utils.load_model(config["ANALYSIS"]['pathway_model_description'])


    # for i_tool,tool in enumerate(tools):
    #     sampling_coverage.r2_iteration(path_samples[i_tool],tool,tag[i_tool],model,dose_level,replicates,molecules,path_dar_files)
    
    ## intra_molecules & inter DAR identification methods on mana and riptide
    print("\n******* between MANA and RIPTiDe **********")
    path_dar_riptide = config["COMPARISON"]["samples_riptide"]
    path_dar_iMAT = config["COMPARISON"]["samples_mana"]
    output = config["COMPARISON"]["output"]
    sampling_coverage.intra_molecules_inter_dar_and_context_methods(path_dar_riptide,path_dar_iMAT,model,molecules,replicates,path_dar_files,output)

def add_media_constrains(model,media):
    medium = model.medium
    new_medium = {}
    metabolite_name = media['reaction_id']
    for r_media in metabolite_name:
        if r_media in medium:
            new_medium[r_media] = 1000
    model.medium = new_medium


def get_flux_samples(model,transcriptome,samples=int(),maxfit=bool(),out_riptide=str(),frac_step=0.1):
    """ Function to run RIPTiDe from either TG-GATES MicroAray or RNA seq data. 

    Args:
        data_path (str): Path of the microArray data
        attribute_data_path (str): Path of the attribute microArray data
        sacrific_period (str): time of incubation of the molecules. Could be 3, 6, 9 or 24 hr . Defaults to str().
        dose_level (str): Dose level of the molecule. Could be Control, Low, Middle, High. Defaults to str().
        compound_name (str): Name of the molecule. Could be amiodarone or valproic acid. Defaults to str().
        sampling_coverage (bool, optional): Perform sampling metrics. True|False. Defaults to False.
        replicates (int): Number of replicates in the data. Defaults to int().
        samples (int): Number of desired samples to perform. If maxfit, 101 is the minimum. Defaults to int().
        maxfit (bool): Performs maxfit algorithm of RIPTiDe. Defaults to bool().
    """    
    
    print(frac_step)
    if maxfit:
        riptide_object = riptide.maxfit(model,transcriptome=transcriptome,objective=True,prune=True,gpr=True,samples=samples,frac_step=frac_step)
        if riptide_object.concordance['p'] > 0.05:
            frac_step-=0.02
            get_flux_samples(model,transcriptome,samples=samples,maxfit=maxfit,out_riptide=out_riptide,frac_step=frac_step)
    else:
        # fraction = utils.get_fraction("results/riptide/recon2.2/maxfit/"+compound_name+'/'+str(samples)+'/'+sacrific_period.split(' ')[0]+'_'+dose_level+'/replicate_'+str(reps)+'/parameters.txt')
        # out = "results/riptide/recon2.2/"+compound_name+'/'+str(samples)+'/'+sacrific_period.split(' ')[0]+'_'+dose_level+'/replicate_'+str(reps)+'/'
        riptide_object = riptide.contextualize(model,transcriptome=transcriptome,objective=True,prune=True,gpr=True,samples=samples,fraction=0.8)

    riptide.save_output(riptide_obj=riptide_object, path=out_riptide,file_type='SBML')
            


