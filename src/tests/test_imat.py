
#
#
# Script to manipulate iMAT method implemented in dexom_python. Inspired from Louison project on git.
# 4.10.2024
# Made by Maxime Lecomte
#
#

from dexom_python import imat
from cobra import *

model = io.read_sbml_model('data/tests/model_riptide_complex.sbml')
#define parameters
eps = 1e-2  # threshold of activity for highly expressed reactions
thr = 1e-5  # threshold of activity for all reactions
obj_tol = 1e-5  # variance allowed for the objective_value
tlim = 600  # time limit (in seconds) for the imat model.optimisation() call
tol = 1e-6  # tolerance for the solver
maxiter = 10
model.solver = "cplex"

#define reactions weights
reaction_names = []
for reaction in model.reactions:
    reaction_names.append(reaction.id)
reaction_weights = {}
rh_reactions = ['rxn1']
rl_reactions = ['rxn2', 'rxn3']

for rname in reaction_names:
    if rname in rh_reactions:
        reaction_weights[rname] = 1.
    elif rname in rl_reactions:
        reaction_weights[rname] = -1.

#run imat
dexom_solutions  = imat(model,reaction_weights,epsilon=eps)
print(dexom_solutions.fluxes)
