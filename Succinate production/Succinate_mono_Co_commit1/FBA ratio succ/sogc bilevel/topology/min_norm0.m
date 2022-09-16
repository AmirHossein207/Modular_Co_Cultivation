clear
clc

solverOK = changeCobraSolver('gurobi','lp');
load('iJO1366.mat');
model = iJO1366;


model1 = model;
%% find candidate minimal flux for production
[non_zero_flux_D,FBA_solution_D,non_zero_flux_W,FBA_solution_W] =...
    pathway_min_norm(model1);
printFluxVector(model1, FBA_solution_W, true, true,-1,[],[],true)
%% visualization of flux
[Involved_mets_D, Dead_ends_D] = draw_by_rxn(model1, non_zero_flux_D, 'true', ...
    'struc', {'succ_c'}, {'h_c','nad_c','nadh_c','h_e','h2o_c','h2o_e','h2o_p','h_p',...
    'atp_c','adp_c','co2_c','co2_e','co2_p'}, FBA_solution_D);
% [Involved_mets_W, Dead_ends_W] = draw_by_rxn(model1, non_zero_flux_W, 'true', ...
%     'struc', {'pep_p'}, {'h_c','nad_c','nadh_c','h_e','h2o_c','h2o_e','h2o_p','h_p',...
%     'atp_c','adp_c','co2_c','co2_e','co2_p'}, FBA_solution_W);

%% creat the graph of the flux
% [flux_model_D,mets_index_D,rxns_index_D] = creat_flux_graph(model1,Involved_mets_D,non_zero_flux_D);
% 
% [flux_model_W,mets_index_W,rxns_index_W] = creat_flux_graph(model1,Involved_mets_W,non_zero_flux_W);

% 
% [directionRxns_D, involvedMets_D, deadEnds_D] = draw_by_met(model1, {'succ_c'},...
%     'true', 1, 'prod', {''}, FBA_solution_D);
[directionRxns_W, involvedMets_W, deadEnds_W] = draw_by_met(model1, {'succ_c'},...
    'true', 1, 'prod', {''}, FBA_solution_W);
