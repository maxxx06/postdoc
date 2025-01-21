#
#
# Script to use riptide on microarray data 
# 19.11.2024
# Made by Maxime Lecomte
#
#

import cobra
import os
import riptide
from pathlib import Path
import data_management

def load_model(model_path):
    """ load and create cobra SBML model

    Args:
        model_path (str): Path of the SBML model

    Returns:
        cobra.core.Model
    """    
    return cobra.io.read_sbml_model(model_path)

def create_directory_if_not_exists(path):
    """ Create directory path from input path and check if the directory already exists

    Args:
        path (str): Path of the desired directory

    Returns:
        Bool
    """    

    if not os.path.exists(os.path.dirname(path)):
        os.makedirs(os.path.dirname(path))
        print(f"Le répertoire '{path}' a été créé.")
    else:
        print(f"Le répertoire '{path}' existe déjà.")
    return True


def get_fraction(path):
    """Get fraction of optimum objective value from previous maxfit optimization

    Args:
        path (str): Path of the desired directory

    Returns:
        int: Fraction of optimum objective value
    """    
    if os.path.exists(path):
        f = open(path)
        line = f.readline()
        return line.split(' ')[-1]
    else:
        print(f'No optimized fraction of optimum was find because the path {path} is invalid.\nFraction will be set to the defaut value of 0.8.')
        return 0.8
        

def get_flux_samples(data_path,attribute_data_path,sacrific_period=str(),dose_level=str(),compound_name=str(),sampling_coverage=False,replicates=int(),samples=int(),maxfit=bool()):
    """ Function to run RIPTiDe from TG-GATES MicroAray data. 

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
    
    output_transcriptomic_path = 'data/microarray/'+compound_name+'/'+sacrific_period.split(' ')[0]+'_'+dose_level+'/'
    model = load_model('data/metabolic_networks/recon2v2_biomass_corrected_final.sbml')   
    model.name = 'data/metabolic_networks/recon2v2_biomass_corrected_final.sbml'

    if create_directory_if_not_exists(output_transcriptomic_path):
        transcriptomic_data = data_management.data_filter(data_path,attribute_data_path,output_transcriptomic_path+'transcriptomic_data.csv',sacrifice_period=sacrific_period,dose_level=dose_level,compound_name=compound_name)
        transcriptomic_data_mapped = data_management.map_genes(transcriptomic_data,mapped_genes='data/microarray/gene_with_protein_product.txt',output=[output_transcriptomic_path+'transcriptomic_data_mapped.csv',output_transcriptomic_path+'transcriptomic_data_mapped_corrected.csv'])
        
        for reps in range(1,replicates,1):
            transcriptomic_data_mapped_reps = {k:v[reps] for k,v in transcriptomic_data_mapped.items()}
            if maxfit:
                out = "results/riptide/recon2.2/maxfit/"+compound_name+'/'+str(samples)+'/'+sacrific_period.split(' ')[0]+'_'+dose_level+'/replicate_'+str(reps)+'/'
                riptide_object = riptide.maxfit(model=model,transcriptome=transcriptomic_data_mapped_reps,objective=True,prune=True,gpr=True,samples=samples)

            elif maxfit == False:
                fraction = get_fraction("results/riptide/recon2.2/maxfit/"+compound_name+'/'+str(samples)+'/'+sacrific_period.split(' ')[0]+'_'+dose_level+'/replicate_'+str(reps)+'/parameters.txt')
                out = "results/riptide/recon2.2/"+compound_name+'/'+str(samples)+'/'+sacrific_period.split(' ')[0]+'_'+dose_level+'/replicate_'+str(reps)+'/'
                riptide_object = riptide.contextualize(model=model,transcriptome=transcriptomic_data_mapped_reps,objective=True,prune=True,gpr=True,samples=samples,fraction=fraction)

            if create_directory_if_not_exists(out):
                riptide.save_output(riptide_obj=riptide_object, path=out,file_type='SBML')

    else:
        print('not a valid path for ', output_transcriptomic_path)

    if sampling_coverage:
        import sampling_coverage
        sampling_coverage.r2_iteration()
