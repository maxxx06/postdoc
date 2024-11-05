
% 
% tutorial integration tools 
% Made by Maxime Lecomte
% Oct. 2024
%


%changeCobraSolver ('glpk', 'all');
modelFileName = 'ecoli_core_model.mat';
modelDirectory = getDistributedModelFolder(modelFileName); %Look up the folder for the distributed Models.
modelFileName= [modelDirectory filesep modelFileName]; % Get the full path. Necessary to be sure, that the right model is loaded
model = readCbModel(modelFileName);
[expressionRxns parsedGPR] = mapExpressionToReactions(model, expression);

load('dataEcoli');

solverName = 'ibm_cplex';
solverType = 'LP';
printLevel = 2;
validationLevel = 2;

[solverOK, solverInstalled] = changeCobraSolver(solverName,solverType, printLevel, validationLevel);

% iMat option
options = 'options_iMAT';
load(['options_methods' filesep options]);
iMAT_model = createTissueSpecificModel(model, options);

% GIMME option
options = 'options_GIMME';
load(['options_methods' filesep options]);
GIMME_model = createTissueSpecificModel(model, options);


% MBA option
options = 'options_MBA';
load(['options_methods' filesep options]);
options.tol = 1e-6;
MBA_model = createTissueSpecificModel(model, options);

% Fascore option
options = 'options_fastCore';
load(['options_methods' filesep options]);
FastCore_model = createTissueSpecificModel(model, options);

[minFlux, maxFlux]=fluxVariability(FastCore_model);
disp(FastCore_model.S);
disp(length(FastCore_model.rxns));
disp([minFlux,maxFlux]);