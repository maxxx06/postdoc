
#
#
# Script to manipulate RIPTiDe method 
# 4.10.2024
# Made by Maxime Lecomte
#
#


import riptide
import cobra

def modelFromTutorial():
    my_model = cobra.io.read_sbml_model('data/tests/model_riptide.sbml')

    transcript_abundances_1 = riptide.read_transcription_file('data/tests/transcriptome1.tsv')
    # transcript_abundances_2 = riptide.read_transcription_file('data/tests/transcriptome1.tsv') # has replicates

    riptide_object_1_a = riptide.contextualize(model=my_model, transcriptome=transcript_abundances_1)
    riptide_object_1_b = riptide.contextualize(model=my_model, transcriptome=transcript_abundances_1, tasks=['rxn1'], exclude=['rxn2','rxn3'])

    riptide.save_output(riptide_obj=riptide_object_1_b, path='data/tests/example_riptide_output')    

def addComplexGPR():
    my_model = cobra.io.read_sbml_model('data/tests/model_riptide_complex.sbml')

    transcript_abundances_1 = riptide.read_transcription_file('data/tests/transcriptome1.tsv')
    # transcript_abundances_2 = riptide.read_transcription_file('data/tests/transcriptome1.tsv') # has replicates

    riptide_object_1_a = riptide.contextualize(model=my_model, transcriptome=transcript_abundances_1)
    riptide_object_1_b = riptide.contextualize(model=my_model, transcriptome=transcript_abundances_1, tasks=['rxn1'], exclude=['rxn2','rxn3'])

    riptide.save_output(riptide_obj=riptide_object_1_b, path='data/tests/example_riptide_output')  

def usingGPR():   
    my_model = cobra.io.read_sbml_model('data/tests/model_riptide_complex.sbml')

    transcript_abundances_1 = riptide.read_transcription_file('data/tests/transcriptome1.tsv')
    # transcript_abundances_2 = riptide.read_transcription_file('data/tests/transcriptome1.tsv') # has replicates
    max_fit = riptide.maxfit(model=my_model,transcriptome=transcript_abundances_1)
    riptide_object_1_a = riptide.contextualize(model=my_model, transcriptome=transcript_abundances_1)
    # riptide_object_1_b = riptide.contextualize(model=my_model, transcriptome=transcript_abundances_1, tasks=['rxn1'],gpr=True, additive=True)
    riptide.save_output(riptide_obj=riptide_object_1_a, path='data/tests/example_riptide_output')  
    fsamples_specific_model= riptide_object_1_a.flux_samples
    fraction_bounds_specific_model = max_fit.fraction_bounds
    print(fsamples_specific_model)
    print(fraction_bounds_specific_model)
