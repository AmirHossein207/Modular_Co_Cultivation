function [F,v_g,v_int_m,v_int_g,v_txdn,v_o,v_z] = fba_txdn(model_1)
rxnIDs = findRxnIDs(model_1,{'EX_glc__D_e','Ex_txdn_e','EX_o2_e','EX_mva_e',...
    'EX_g6p_e','BIOMASS_Ec_iJO1366_core_53p95M','EX_xyl__D_e'});
solution = optimizeCbModel(model_1);

if ~isempty(solution.x)
    flux = solution.x;
    v_g = flux(rxnIDs(1));
    v_txdn = flux(rxnIDs(2));
    v_o = flux(rxnIDs(3));
    v_int_m = flux(rxnIDs(4));
    v_int_g = flux(rxnIDs(5));
    F = flux(rxnIDs(6));
    v_z = flux(rxnIDs(7));

else
    F = 0;
    v_g = 0;
    v_o = 0;
    v_txdn = 0;
    v_int_m = 0;
    v_int_g = 0;
    v_z = 0;
end
end
