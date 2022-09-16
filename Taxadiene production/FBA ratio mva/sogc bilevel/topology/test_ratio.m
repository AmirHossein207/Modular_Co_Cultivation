clear
clc
close all
solverOK = changeCobraSolver('gurobi','lp');
load('iJO1366.mat');
model = iJO1366;
model = addReaction(model,'EX_mva_e', {'mva_e'}, -1, true, 0,1000);
% model = addReaction(model,'MEVRt','mva_c ⇌ mva_e');
model = addReaction(model,'MEVRt', {'mva_c','mva_e'}, [1 -1], true,-1000, 1000);
% model = addReaction(model,'MVNOR', 'coa_c + nad_c + mva_c → 4.0 h_c + nadh_c + hmgcoa_c');
% model_up = addReaction(model_up,'MVNOR', {'coa_c','nad_c','mva_c','h_c','nadh_c','hmgcoa_c'},[-1;-1;-1;4;1;1],false,0,1000);
% model = addReaction(model,'pksg', 'aacoa_c + accoa_c + h2o_c -> h_c + coa_c + hmgcoa_c');
model = addReaction(model,'pksg', {'aacoa_c','accoa_c','h2o_c','h_c','coa_c','hmgcoa_c'},[-1;-1;-1;1;1;1],false,0,1000);
% model = addReaction(model,'hmg1', 'hmgcoa_c  + 2 nadph_c  + 2 h_c  -> 2 nadp_c  + coa_c  + mva_c');
model = addReaction(model,'hmg1', {'hmgcoa_c','nadph_c','h_c','nadp_c', 'coa_c', 'mva_c'},[-1;-2;-2;2;1;1],false,0,1000);




S = model.S;
rxnIDs = findRxnIDs(model,{'pksg','ACCOAC','pksg','CS',...
    'pksg','IPPS','ggpps','IPDDI','DMATT','GRTT',...
    'OCTDPS','UDCPDPS','GAPD','DXPS','GAPD','TALA'});
b1 = 0.99;
b2 = 0.99;
b3 = 0.99;

ss = zeros(1,length(model.rxns));
ss(1,rxnIDs(1:2)) = [1-b1 -b1]; 
S(end+1,:) = ss;
ss = zeros(1,length(model.rxns));
ss(1,rxnIDs(3:4)) = [1-b2 -b2]; 
S(end+1,:) = ss;
ss = zeros(1,length(model.rxns));
ss(1,rxnIDs(5:6)) = [1-b3 -b3]; 
S(end+1,:) = ss;

% ss = zeros(1,length(model.rxns));
% ss(1,rxnIDs(7:12)) = [1-b4 -b4 -b4 -b4 -b4 -b4]; 
% S(end+1,:) = ss;

model.S = S;

model = addMetabolite(model, 'iodine1_c');
model = addMetabolite(model, 'iodine2_c');
model = addMetabolite(model, 'iodine3_c');
% model = addMetabolite(model, 'iodine4_c');
% model = addMetabolite(model, 'iodine5_c');

model = changeRxnBounds(model, {'EX_glc__D_e','EX_o2_e'},[-10, -15],'l');
model = changeObjective(model,'BIOMASS_Ec_iJO1366_core_53p95M');
solution_1 = optimizeCbModel(model);
printFluxVector(model, solution_1.x, true, true,-1,[],[],true)
fprintf('\n')
[directionRxns_W, involvedMets_W, deadEnds_W] = draw_by_met(model, {'accoa_c'},...
    'true', 1, 'prod', {'iodine1_c', 'iodine2_c', ...
    'iodine3_c', 'coa_c', 'pi_c','actp_c'}, solution_1.x);
% [minFlux, maxFlux] = fluxVariability(model, 60, 'max', {'g3p_ex','DXPS','GAPD'})














