clear 
clc

solverOK = changeCobraSolver('gurobi','lp');
load('iJO1366.mat');
model = iJO1366;
% model = addReaction(model, 'PEPtex', 'metaboliteList', {'pep_p','pep_e'},...
%     'stoichCoeffList', [-1; 1],'lowerBound',-1000,'upperBound',1000);
% model = addReaction(model, 'PEPt6pp', 'metaboliteList', {'pep_c','pi_c','pi_p','pep_p'},...
%     'stoichCoeffList', [-1; 1; -1; 1],'lowerBound',-1000,'upperBound',1000);
% [model, EX_pep_e] = addExchangeRxn(model, 'pep_e',-10,0);


b1 = 0.99;
b2 = 0.99;
b3 = 0.99;
rxnIDs = findRxnIDs(model,{'ICL','ICDHyr',...
    'SUCCt3pp','SUCDi',...
    'FRD2','FUM'});

S = model.S;

ss = zeros(1,length(model.rxns));
ss(1,rxnIDs(1:2)) = [1-b1 -b1];
S(end+1,:) = ss;
ss = zeros(1,length(model.rxns));
ss(1,rxnIDs(3:4)) = [1-b2 -b2];
S(end+1,:) = ss;
ss = zeros(1,length(model.rxns));
ss(1,rxnIDs(5:6)) = [1-b3 -b3];
S(end+1,:) = ss;



model.S = S;

model = addMetabolite(model, 'iodine1_c');
model = addMetabolite(model, 'iodine2_c');
model = addMetabolite(model, 'iodine3_c');

model = changeObjective(model,'BIOMASS_Ec_iJO1366_core_53p95M');
solution_1 = optimizeCbModel(model);
printFluxVector(model, solution_1.x, true, true,-1,[],[],true)
fprintf('\n')
[directionRxns_D, involvedMets_D, deadEnds_D] = draw_by_met(model, {'fum_c'},...
    'true', 1, 'sub', {'iodine1_c','iodine2_c','iodine3_c'}, solution_1.x);
% [minFlux, maxFlux] = fluxVariability(model, 60, 'max', ...
%     {'EX_g6p_e','G6Pt6_2pp_reverse','UAGCVT','PSCVT','PPC','KDOPS','DDPA','DHAPT','PYK'})










