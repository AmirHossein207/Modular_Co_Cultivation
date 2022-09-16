function [F_up,v_g_up,vz_up,v_int_up,v_o_up,v_phe_up,v_tyr_up] = fba_tyr(model_1)
rxnIDs = findRxnIDs(model_1,{'EX_glc__D_e','EX_o2_e','EX_4_tyrosol_e',...
    'BIOMASS_Ec_iJO1366_core_53p95M','EX_xyl__D_e','EX_phe__L_e','EX_tyr__L_e'});
solution = optimizeCbModel(model_1);

if ~isempty(solution.x)
    flux = solution.x;
    v_g_up = flux(rxnIDs(1));
    v_o_up = flux(rxnIDs(2));
    v_int_up = flux(rxnIDs(3));
    F_up = flux(rxnIDs(4));
    vz_up = flux(rxnIDs(5));
    v_phe_up = flux(rxnIDs(6));
    v_tyr_up = flux(rxnIDs(7));
else
    F_up = 0;
    v_g_up = 0;
    v_o_up = 0;
    v_int_up = 0;
    vz_up = 0;
    v_phe_up = 0;
    v_tyr_up = 0;
end
end