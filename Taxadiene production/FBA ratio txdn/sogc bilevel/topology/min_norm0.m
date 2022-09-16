clear
clc
close all
solverOK = changeCobraSolver('gurobi','lp');
load('iJO1366.mat');
model = iJO1366;

model = addReaction(model,'pksg', 'aacoa_c + accoa_c + h2o_c -> h_c + coa_c + hmgcoa_c');
model = addReaction(model,'hmg1', 'hmgcoa_c  + 2 nadph_c  + 2 h_c  -> 2 nadp_c  + coa_c  + mva_c');
model = addReaction(model,'mvak1', 'mva_c  + atp_c -> adp_c + h_c mvap_c ');
model = addReaction(model,'mvak2', 'mvap_c  + atp_c -> adp_c  + mvapp_c ');
model = addReaction(model,'mvad', 'mvapp_c  + atp_c -> adp_c  + co2_c  + pi_c + ipdp_c ');
model = addReaction(model,'ggpps', 'frdp_c  + ipdp_c  -> ggpp_c  + ppi_c');
model = addReaction(model,'txs', 'ggpp_c  -> txdn_c  + ppi_c ');
model = addReaction(model,'txdn_ex', {'txdn_c','txdn_e'}, [-1 1], true, -1000, 1000);
model = addReaction(model,'Ex_txdn_e', {'txdn_e'}, -1, false, 0, 1000);


%% find candidate minimal flux for production
[non_zero_flux_D,FBA_solution_D,non_zero_flux_W,FBA_solution_W] =...
    pathway_min_norm(model);
printFluxVector(model, FBA_solution_D, true, true,-1,[],[],true)
fprintf('\n')
% printFluxVector(model, FBA_solution_W, true, true,-1,[],[],true)
 %% visualization of flux
% [Involved_mets_D, Dead_ends_D] = draw_by_rxn(model, non_zero_flux_D, 'true', ...
%     'struc', {'txdn_c'}, {'nadp_c','nadph_c','pi_c','h_c','nad_c','nadh_c','h_e','h2o_c','h2o_e','h2o_p','h_p',...
%     'atp_c','adp_c','co2_c','co2_e','co2_p'}, FBA_solution_D);
% [Involved_mets_W, Dead_ends_W] = draw_by_rxn(model, non_zero_flux_W, 'true', ...
%     'struc', {'txdn_c'}, {'h_c','nad_c','nadh_c','h_e','h2o_c','h2o_e','h2o_p','h_p',...
%     'atp_c','adp_c','co2_c','co2_e','co2_p'}, FBA_solution_W);

%% creat the graph of the flux
% [flux_model_D,mets_index_D,rxns_index_D] = creat_flux_graph(model,Involved_mets_D,non_zero_flux_D);
% 
% [flux_model_W,mets_index_W,rxns_index_W] = creat_flux_graph(model,Involved_mets_W,non_zero_flux_W);


% [directionRxns_D, involvedMets_D, deadEnds_D] = draw_by_met(model, {'hmgcoa_c'},...
%     'true', 1, 'sub', {'h_c','nad_c','nadph_c','nadh_c','h_e','h2o_c','h2o_e','h2o_p','h_p',...
% 'atp_c','adp_c','co2_c','co2_e','co2_p'}, FBA_solution_D);
[directionRxns_W, involvedMets_W, deadEnds_W] = draw_by_met(model, {'pyr_c'},...
    'true', 1, 'sub', {'h_c','nad_c','nadph_c','nadh_c','h_e','h2o_c','h2o_e','h2o_p','h_p',...
'atp_c','adp_c','co2_c','co2_e','co2_p'}, FBA_solution_W);

% [minFlux, maxFlux] = fluxVariability(model, 20, 'max', {'EX_txdn_c'})







