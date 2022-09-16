function [non_zero_flux_D,FBA_solution_D,non_zero_flux_W,FBA_solution_W] =...
    pathway_min_norm(model)
%% wild type
model_1 = model;
model_1 = changeRxnBounds(model_1, {'EX_glc__D_e','EX_xyl__D_e','EX_o2_e'...
    ,'EX_tyr__L_e','EX_4_tyrosol_e','EX_phe__L_e'},[-10,0,-15,-5,-5,5]...
    ,{'l','l','l','l','l','u'});
model_1 = changeObjective(model_1,'BIOMASS_Ec_iJO1366_core_53p95M');%
solution_2 = optimizeCbModel(model_1,'max');
% solution_2 = optimizeCbModel(model_1,'max','zero','0','log');
FBA_solution_W = solution_2.x;
index = find(FBA_solution_W~=0);
non_zero_flux_W = model_1.rxns(index);
% v1(abs(v1)<(1e-6)) = 0;
NO_of_non_zero_fluxes = sum(FBA_solution_W==0);
%% desired
model_1 = model;
% the bounds of the feed and the product are set using changeRxnBounds
model_1 = changeRxnBounds(model_1, {'EX_glc__D_e','EX_xyl__D_e','EX_o2_e'...
    ,'EX_tyr__L_e','EX_4_tyrosol_e','EX_phe__L_e'},[-10,0,-15,-5,-5,5]...
    ,{'l','l','l','l','l','u'});
% model_1 = changeObjective(model_1, {'EX_g3p_e',...
%     'BIOMASS_Ec_iJO1366_core_53p95M'},[0.013569, 0.26588]);
model_1 = changeObjective(model_1,'EX_salid_e');%
% FBA1 = optimizeCbModel(model_1, 'max');
% solution_1 = optimizeCbModel(model_1,'max');
solution_1 = optimizeCbModel(model_1,'max','zero','0','log');
FBA_solution_D = solution_1.x;
index = find(FBA_solution_D~=0);
non_zero_flux_D = model_1.rxns(index);
% v1(abs(v1)<(1e-6)) = 0;
NO_of_non_zero_fluxes = sum(FBA_solution_D==0);


end