#
#
# Script to manipulate DEXOM_python method 
# 4.10.2024
# Made by Maxime Lecomte
#
#

import dexom_python
import cobra
import pandas as pd
import riptide
import tempfile

from dexom_python.default_parameter_values import DEFAULT_VALUES
from dexom_python.result_functions import read_solution
from dexom_python.enum_functions.rxn_enum_functions import rxn_enum
from dexom_python.enum_functions.icut_functions import icut
from dexom_python.enum_functions.diversity_enum_functions import diversity_enum

def load_model():
    model = cobra.io.read_sbml_model("data/tests/model_riptide_complex.sbml")
    pd_gene_w = pd.read_csv("data/tests/transcriptome1.tsv",sep='\t',header=None)
    gene_w = pd_gene_w.to_dict('index')
    gene_w_good = {list(v.values())[0]:list(v.values())[1] for k,v in gene_w.items()} # get the relationship gene:abundance

    return gene_w_good,model


def playGPR_rules():
    """Apply GPR rule according Dexom.

    and: min
    or : max

    Returns:
        tuple : reaction weights (CSV) and cobrapy model
    """
    gene_w_good,model = load_model()   
    return dexom_python.apply_gpr(model,gene_w_good,save=True,filename='data/tests/dexom_output/reaction_weigths')


def enumeration_methods(provenance="iMAT",enum_method = 'maxdist'):
    """perform reaction enumeration and diversity enumeration from iMAT or riptide solution.

    Args:
        min_max_reaction_weights (dict, optional): mandatory if tool option is riptide. Contains minimal and maximal reaction weights. 
        reaction_weights_iMAT (dict, optional): mandatory if tool option is iMAT. Contains reaction weights. 
        provenance (str, optional): indicate which methods is used to genrating solution. Defaults to "iMAT". Can be either 'riptide' or 'iMAT'.
    """    
    eps = DEFAULT_VALUES['epsilon']  # threshold of activity for highly expressed reactions
    thr = DEFAULT_VALUES['threshold']  # threshold of activity for all reactions
    obj_tol = DEFAULT_VALUES['obj_tol']  # variance allowed for the objective_value
    maxiter = 5  # maximum number of iterations
    dist_anneal = 0.9  # diversity-enumeration parameter

    _, model = load_model()

    if provenance == 'riptide':
        riptide_object = lunch_riptide('data/tests/transcriptome1.tsv',model)
        tempory_csv_one_solution = create_tempory_previous_solution_file(riptide_object.flux_samples)

        # read_solution allows the generation of a solution object from previous flux solution
        prev_solution = read_solution(tempory_csv_one_solution,solution_index=0)[0] # only the first element of this tuple is a solution class
        reaction_weights = riptide_object.reaction_weights
        print('prev_solution of RIPTIDE fluxes \n',prev_solution.fluxes)
        print('prev_solution of RIPTIDE objective value \n',prev_solution.objective_value)
        print('reaction weights RIPTIDE \n',reaction_weights)

    else:
        prev_solution=None
        reaction_weights = playGPR_rules()
        # reaction_weights = reaction_weights_iMAT
        print('reaction weights iMAT \n',reaction_weights)
        
    ### Use different enumeration methods of DEXOM
    
    if enum_method == 'rxn':
        #Reaction enum: Based on the idea of generating alternative solutions by single reaction changes.
        rxn_sol  = rxn_enum(model=model, prev_sol=prev_solution, reaction_weights=reaction_weights, eps=eps,thr=thr, obj_tol=obj_tol, tool='riptide', transcriptomic_file='data/tests/transcriptome1.tsv',save=True)

        uniques = pd.DataFrame(rxn_sol.unique_binary)
        uniques.columns = [r.id for r in model.reactions]
        uniques.to_csv('data/tests/dexom_output/'+provenance+'_rxnenum_solutions.csv')

    elif enum_method == 'icut':
        # Icut-enum: Using integer-cuts as constrains to discards previously solutions
        icut_sol = icut(model,reaction_weights, prev_sol=prev_solution, eps=eps,thr=thr, obj_tol=obj_tol, tool=provenance, transcriptomic_file='data/tests/transcriptome1.tsv',full=True)
        uniques = pd.DataFrame(icut_sol.binary)
        sol = pd.DataFrame(icut_sol.solutions)
        print(icut_sol.objective_value)
        print(icut_sol.solutions[0].fluxes)
        print(icut_sol.binary)

        uniques.columns = [r.id for r in model.reactions]
        uniques.to_csv('data/tests/dexom_output/'+provenance+'_icut_binary_sol.csv')

    elif enum_method == 'maxdist':
        # Most distant optimal solution with respect to the previous optimal solution, and using integer-cuts to avoid re-discovering the same distant solutions.
        #TODO
        pass

    elif enum_method == 'diversity':
        # takes the best of the other three techniques without their disadvantages.
        div_sol, div_res = diversity_enum(model=model, prev_sol=prev_solution, reaction_weights=reaction_weights, eps=eps,
                                        thr=thr, obj_tol=obj_tol, maxiter=maxiter, dist_anneal=dist_anneal)
        div_res.to_csv('data/tests/dexom_output/'+provenance+'_divenum_results.csv')
        sol = pd.DataFrame(div_sol.binary, columns=[r.id for r in model.reactions])
        sol.to_csv('data/tests/dexom_output/'+provenance+'_divenum_solutions.csv')
    



def create_tempory_previous_solution_file(df):
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as temp_file:
        df.to_csv(temp_file.name)
        temp_file_path = temp_file.name
        print(f"Fichier CSV temporaire créé : {temp_file_path}")

    return temp_file_path

def lunch_riptide(file_transcrit,model):
    transcript_abundances_rpm = riptide.read_transcription_file(file_transcrit)
    temp_sol_full = riptide.contextualize(model=model, transcriptome = transcript_abundances_rpm,objective=False,fraction=0.001,gpr=False,samples=100,prune=False) # objective = True --> for biomass objective function / fraction --> fraction of optimum objective value / gpr = True --> take into account GPR in the calculation of reaction weights / samples = 2 --> generate only 2 alternative solutions

    # riptide.save_output(riptide_obj=temp_sol_full, path='data/tests/example_riptide_output') 
    temp_sol_full.reaction_weights = {key: int(value[0]) for key, value in temp_sol_full.reaction_weights.items()}

    return temp_sol_full
    

enumeration_methods(provenance="iMAT",enum_method ='maxdist')