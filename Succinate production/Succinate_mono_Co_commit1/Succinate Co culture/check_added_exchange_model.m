clc
close 
clear
SolverOk = changeCobraSolver('gurobi','lp');
model = readCbModel('iJO1366.mat');

%% intermediate strain
model_up = model;
%% add exchange reaction for intermediate
lb = -100;
ub = 1000;
[model_up, EX_pep_e] = addExchangeRxn(model_up, 'pep_c', lb, ub);
%% test the flux through the exchange
model_up = changeRxnBounds(model_up,{'EX_glc__D_e','EX_o2_e'}...
    ,[-15 -20],{'l','l'});
model_up = changeObjective(model_up, 'EX_pep_c');

FBA1 = optimizeCbModel(model_up, 'max');
printFluxVector(model_up, FBA1.x, true, true,-1,[],[],true)
fprintf('\n')
FBA2 = optimizeCbModel(model_up, 'min');
printFluxVector(model_up, FBA2.x, true, true,-1,[],[],true)
%% Identify flux inconsistent reactions
[fluxConsistentMetBool, fluxConsistentRxnBool, fluxInConsistentMetBool, ...
fluxInConsistentRxnBool] = findFluxConsistentSubset(model_up);