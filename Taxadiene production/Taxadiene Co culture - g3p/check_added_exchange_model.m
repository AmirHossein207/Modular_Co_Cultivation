clc
close 
clear
SolverOk = changeCobraSolver('gurobi','lp');
model = readCbModel('iJO1366.mat');

%% intermediate strain
model_up = model;
%% add exchange reaction for intermediate
model_up  =  addReaction(model_up,'g3p_ex', {'g3p_c','g3p_e'}, [-1 1], true, -1000, 1000);
lb = 0;
ub = 1000;
[model_up, EX_g3p_e] = addExchangeRxn(model_up, 'g3p_e', lb, ub);
%% test the flux through the exchange

model_up = changeRxnBounds(model_up,{'EX_glc__D_e','EX_o2_e'}...
    ,[-15 -20],{'l','l'});
model_up = changeObjective(model_up, 'BIOMASS_Ec_iJO1366_core_53p95M');

rxnIDs = findRxnIDs(model_up,{'g3p_ex','DXPS','GAPD'});
S = model_up.S;
ss = zeros(1,length(model_up.rxns));
ss(1,rxnIDs(1)) = 1-5; 
ss(1,rxnIDs(2:end)) = 5;
S(end+1,:) = ss;
model_up.S = S;
model_up = addMetabolite(model_up, 'iodine_c');

FBA1 = optimizeCbModel(model_up, 'max');
printFluxVector(model_up, FBA1.x, true, true,-1,[],[],true)
% fprintf('\n')
% FBA2 = optimizeCbModel(model_up, 'min');
% printFluxVector(model_up, FBA2.x, true, true,-1,[],[],true)
%% Identify flux inconsistent reactions
% [fluxConsistentMetBool, fluxConsistentRxnBool, fluxInConsistentMetBool, ...
% fluxInConsistentRxnBool] = findFluxConsistentSubset(model_up);