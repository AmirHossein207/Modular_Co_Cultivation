clc
close 
clear
SolverOk = changeCobraSolver('gurobi','lp');
model = readCbModel('iJO1366.mat');

%% add exchange reaction for intermediate


model = addReaction(model,'EX_mva_e', {'mva_e'}, -1, true, 0,1000);
% model = addReaction(model,'MEVRt','mva_c ⇌ mva_e');
model = addReaction(model,'MEVRt', {'mva_c','mva_e'}, [1 -1], true,-1000, 1000);
% model = addReaction(model,'MVNOR', 'coa_c + nad_c + mva_c → 4.0 h_c + nadh_c + hmgcoa_c');
% model_up = addReaction(model_up,'MVNOR', {'coa_c','nad_c','mva_c','h_c','nadh_c','hmgcoa_c'},[-1;-1;-1;4;1;1],false,0,1000);
% model = addReaction(model,'pksg', 'aacoa_c + accoa_c + h2o_c -> h_c + coa_c + hmgcoa_c');
model = addReaction(model,'pksg', {'aacoa_c','accoa_c','h2o_c','h_c','coa_c','hmgcoa_c'},[-1;-1;-1;1;1;1],false,0,1000);
% model = addReaction(model,'hmg1', 'hmgcoa_c  + 2 nadph_c  + 2 h_c  -> 2 nadp_c  + coa_c  + mva_c');
model = addReaction(model,'hmg1', {'hmgcoa_c','nadph_c','h_c','nadp_c', 'coa_c', 'mva_c'},[-1;-2;-2;2;1;1],false,0,1000);
% model = addRea

%% find candidate minimal flux for production

[non_zero_flux_D,FBA_solution_D,non_zero_flux_W,FBA_solution_W] =...
    pathway_min_norm(model);
%% visualization of flux
[Involved_mets_D, Dead_ends_D] = draw_by_rxn(model, non_zero_flux_D, 'true', ...
    'struc', {'g3p_c'}, {'h_c','nad_c','nadh_c','h_e','h2o_c','h2o_e','h2o_p','h_p',...
    'atp_c','adp_c','co2_c','co2_e','co2_p'}, FBA_solution_D);
% [Involved_mets_W, Dead_ends_W] = draw_by_rxn(model, non_zero_flux_W, 'true', ...
%     'struc', {'pep_p'}, {'h_c','nad_c','nadh_c','h_e','h2o_c','h2o_e','h2o_p','h_p',...
%     'atp_c','adp_c','co2_c','co2_e','co2_p'}, FBA_solution_W);

%% creat the graph of the flux
% [flux_model_D,mets_index_D,rxns_index_D] = creat_flux_graph(model,Involved_mets_D,non_zero_flux_D);
% 
% [flux_model_W,mets_index_W,rxns_index_W] = creat_flux_graph(model,Involved_mets_W,non_zero_flux_W);

% 
% [directionRxns_D, involvedMets_D, deadEnds_D] = draw_by_met(model, {'3pg_c'},...
%     'true', 1, 'prod', {''}, FBA_solution_D);
% [directionRxns_W, involvedMets_W, deadEnds_W] = draw_by_met(model, {'3pg_c'},...
%     'true', 1, 'prod', {''}, FBA_solution_W);
