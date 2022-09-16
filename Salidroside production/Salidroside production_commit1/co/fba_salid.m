function [F_down,v_g_down,vz_down,v_int_down,v_salid,v_o_down,...
    v_phe_down,v_tyr_down] = fba_salid(model_1)

rxnIDs = findRxnIDs(model_1,{'EX_glc__D_e','EX_salid_e','EX_o2_e','EX_4_tyrosol_e',...
    'BIOMASS_Ec_iJO1366_core_53p95M','EX_xyl__D_e','EX_phe__L_e','EX_tyr__L_e'});
solution = optimizeCbModel(model_1);

if ~isempty(solution.x)
    flux = solution.x;
    v_g_down = flux(rxnIDs(1));
    v_salid = flux(rxnIDs(2));
    v_o_down = flux(rxnIDs(3));
    v_int_down = flux(rxnIDs(4));
    F_down = flux(rxnIDs(5));
    vz_down = flux(rxnIDs(6));
    v_phe_down = flux(rxnIDs(7));
    v_tyr_down = flux(rxnIDs(8));

else
    F_down = 0;
    v_g_down = 0;
    v_o_down = 0;
    v_salid = 0;
    v_int_down = 0;
    vz_down = 0;
    v_phe_down = 0;
    v_tyr_down = 0;
end
end
