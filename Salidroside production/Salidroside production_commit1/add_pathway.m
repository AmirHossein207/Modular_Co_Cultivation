clc
clear
close all

format short g
warning off all

%% Solver checking

solverOk = changeCobraSolver('gurobi','all');

%% Load model
load('iJO1366.mat');
model = iJO1366;
%% add pathway
model = changeRxnBounds(model, {'EX_glc__D_e','EX_o2_e'},[-10, -15],'l');

model = addReaction(model,'TYRCBOX', 'metaboliteList', {'h_c','tyr__L_c',...
    'co2_c','tym_c'}, 'stoichCoeffList', [-1;-1;1;1],...
    'lowerBound',0,'upperBound',1000);
model = addReaction(model,'TYROXDAc', 'metaboliteList', {'h2o_c','o2_c',...
    'tym_c','h2o2_c','nh4_c','4hoxpacd_c'}, 'stoichCoeffList', [-1;-1;-1;1;1;1],...
    'lowerBound',0,'upperBound',1000);
%%%%% add metabolite 4-tyrosol
model = addReaction(model,'TYROSOL', 'metaboliteList', {'4hoxpacd_c','h_c',...
    'nadh_c','nad_c','4_tyrosol_c'}, 'stoichCoeffList', [-1;-1;-1;1;1],...
    'lowerBound',0,'upperBound',1000);
%%%%% add metabolite salidroside
model = addReaction(model,'SALID', 'metaboliteList', {'udpg_c','4_tyrosol_c',...
    'udp_c','h_c','salid_c'}, 'stoichCoeffList', [-1;-1;1;1;1],...
    'lowerBound',0,'upperBound',1000);
model = addReaction(model,'SALIDpp', 'metaboliteList', {'salid_c','salid_e'},...
    'stoichCoeffList', [-1;1;], 'lowerBound',0,'upperBound',1000);
[model, EX_salid_e] = addExchangeRxn(model, 'salid_e',0,1000);

% model = deleteModelGenes(model,{'b1101','b1854','b1676','b2600','b2599'});
model = changeRxnBounds(model,{'EX_glc__D_e','EX_o2_e'}...
    ,[-10 -15],{'l','l'});
model = changeObjective(model, 'EX_salid_e');

FBA1 = optimizeCbModel(model, 'max');
printFluxVector(model, FBA1.x, true, true,-1,[],[],true)
fprintf('\n')
FBA2 = optimizeCbModel(model, 'min');
printFluxVector(model, FBA2.x, true, true,-1,[],[],true)
%% Identify flux inconsistent reactions
% [fluxConsistentMetBool, fluxConsistentRxnBool, fluxInConsistentMetBool, ...
%     fluxInConsistentRxnBool] = findFluxConsistentSubset(model);



