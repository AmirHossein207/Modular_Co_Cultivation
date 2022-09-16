clear 
clc

solverOK = changeCobraSolver('gurobi','lp');
load('iJO1366.mat');
model = iJO1366;

model = addReaction(model, 'G3P_C', 'metaboliteList', {'g3p_c',...
    'g3p_e'},'stoichCoeffList', [-1;1],'lowerBound',...
    0,'upperBound',1000);
model = addReaction(model,'EX_g3p_e', {'g3p_e'}, -1, true, 0,1000);

model = removeRxns(model,'TPI');
model = addReaction(model, 'TPI', 'metaboliteList', {'g3p_c','dhap_c'}...
    ,'stoichCoeffList', [-1;1],'lowerBound',...
    -1000,'upperBound',1000);

rxnIDs = findRxnIDs(model,{'G3P_C','GAPD','G3P_C','TALA','G3P_C',...
    'TPI','TALA','PGI','G6PDH2r'});
b1 = 0.99;

S = model.S;


ss = zeros(1,length(model.rxns));
ss(1,rxnIDs(1:2)) = [1-b1 -b1];
S(end+1,:) = ss;
% % ss = zeros(1,length(model.rxns));
% % ss(1,rxnIDs(3:4)) = [1-b1 -b1];
% % S(end+1,:) = ss;
ss = zeros(1,length(model.rxns));
ss(1,rxnIDs(5:6)) = [1-b1 -b1];
S(end+1,:) = ss;
model.S = S;

model = addMetabolite(model, 'iodine1_c');
model = addMetabolite(model, 'iodine2_c');
model = addMetabolite(model, 'iodine3_c');
save g3p_iJO1366.mat model

model = changeObjective(model,'BIOMASS_Ec_iJO1366_core_53p95M');
solution_1 = optimizeCbModel(model);
printFluxVector(model, solution_1.x, true, true,-1,[],[],false)
fprintf('\n')
[directionRxns_W, involvedMets_W, deadEnds_W] = draw_by_met(model, {'g3p_c'},...
    'true', 1, 'prod', {'iodine3_c', 'iodine2_c', 'iodine1_c'}, solution_1.x);
% [minFlux, maxFlux] = fluxVariability(model, 60, 'max', {'F6Pt6_2pp_reverse'})














