clear
clc

solverOK = changeCobraSolver('gurobi','lp');
load('iJO1366.mat');
model = iJO1366;

% model = addReaction(model, 'G3P_C', 'metaboliteList', {'g3p_c',...
%     'g3p_e'},'stoichCoeffList', [-1;1],'lowerBound',...
%     -1000,'upperBound',1000);
% model = addReaction(model,'EX_g3p_e', {'g3p_e'}, -1, true, 0,1000);
model1 = model;
%% find candidate minimal flux for production




model = addReaction(model, 'G3P_C', 'metaboliteList', {'g3p_c',...
    'g3p_e'},'stoichCoeffList', [-1;1],'lowerBound',...
    0,'upperBound',1000);
model = addReaction(model,'EX_g3p_e', {'g3p_e'}, -1, true, 0,1000);

rxnIDs = findRxnIDs(model,{'G3P_C','GAPD','G3P_C','TALA','G3P_C',...
    'TPI','TALA','PGI','G6PDH2r'});
b1 = 0.99;

S = model.S;


ss = zeros(1,length(model.rxns));
ss(1,rxnIDs(1:2)) = [1-b1 -b1];
S(end+1,:) = ss;
% % ss = zeros(1,length(model_up.rxns));
% % ss(1,rxnIDs(3:4)) = [1-b1 -b1];
% % S(end+1,:) = ss;
ss = zeros(1,length(model.rxns));
ss(1,rxnIDs(5:6)) = [1-b1 -b1];
S(end+1,:) = ss;
model.S = S;

model = addMetabolite(model, 'iodine1_c');
model = addMetabolite(model, 'iodine2_c');
% model = addMetabolite(model, 'iodine3_c');
[non_zero_flux_D,FBA_solution_D,non_zero_flux_W,FBA_solution_W] =...
    pathway_min_norm(model1);
%% visualization of flux

[Involved_mets_D, Dead_ends_D] = draw_by_rxn(model1, non_zero_flux_D, 'true', ...
    'struc', {'g3p_c'}, {'nadp_c','nadph_c','pi_c','h_c','nad_c','nadh_c','h_e','h2o_c','h2o_e','h2o_p','h_p',...
    'atp_c','adp_c','co2_c','co2_e','co2_p'}, FBA_solution_D);
% [Involved_mets_W, Dead_ends_W] = draw_by_rxn(model1, non_zero_flux_W, 'true', ...
%     'struc', {'pep_p'}, {'h_c','nad_c','nadh_c','h_e','h2o_c','h2o_e','h2o_p','h_p',...
%     'atp_c','adp_c','co2_c','co2_e','co2_p'}, FBA_solution_W);

%% creat the graph of the flux
% [flux_model_D,mets_index_D,rxns_index_D] = creat_flux_graph(model1,Involved_mets_D,non_zero_flux_D);
% 
% [flux_model_W,mets_index_W,rxns_index_W] = creat_flux_graph(model1,Involved_mets_W,non_zero_flux_W);

% 
% [directionRxns_D, involvedMets_D, deadEnds_D] = draw_by_met(model1, {'f6p_c'},...
%     'true', 1, 'prod', {''}, FBA_solution_D);
[directionRxns_W, involvedMets_W, deadEnds_W] = draw_by_met(model1, {'g3p_c'},...
    'true', 1, 'prod', {''}, FBA_solution_W);
model1 = changeObjective(model1,'BIOMASS_Ec_iJO1366_core_53p95M');
solution_1 = optimizeCbModel(model1,'max');
printFluxVector(model1, solution_1.x, true, true,-1,[],[],true)
fprintf('\n')