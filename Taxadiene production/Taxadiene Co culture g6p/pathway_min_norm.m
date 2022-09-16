function [non_zero_flux_D,FBA_solution_D,non_zero_flux_W,FBA_solution_W] =...
    pathway_min_norm(model_1)
%% wild type
model_1 = changeRxnBounds(model_1, {'EX_glc__D_e','EX_g6p_e','EX_o2_e'},[0, -10 -15],{'l','l','l'});
model_1 = changeObjective(model_1,'BIOMASS_Ec_iJO1366_core_53p95M');%
solution_2 = optimizeCbModel(model_1,'max');
% solution_2 = optimizeCbModel(model_1,'max','zero','0','log');
FBA_solution_W = solution_2.x;
index = find(FBA_solution_W~=0);
non_zero_flux_W = model_1.rxns(index);
% v1(abs(v1)<(1e-6)) = 0;
NO_of_non_zero_fluxes = sum(FBA_solution_W==0);
%% desired
% the bounds of the feed and the product are set using changeRxnBounds
model_1 = changeRxnBounds(model_1, {'EX_glc__D_e','EX_g6p_e','EX_o2_e'},[0, -10, -15],{'l','l','l'});
% model_1 = changeObjective(model_1, {'Ex_txdn_e',...
%     'BIOMASS_Ec_iJO1366_core_53p95M'},[0.48529,1.1121]);
model_1 = changeObjective(model_1, {'Ex_txdn_e'});
% % FBA1 = optimizeCbModel(model_1, 'max');
% solution_1 = optimizeCbModel(model_1,'max');
solution_1 = optimizeCbModel(model_1,'max','zero','0','log');
FBA_solution_D = solution_1.x;
index = find(FBA_solution_D~=0);
non_zero_flux_D = model_1.rxns(index);
% v1(abs(v1)<(1e-6)) = 0;
NO_of_non_zero_fluxes = sum(FBA_solution_D==0);


end