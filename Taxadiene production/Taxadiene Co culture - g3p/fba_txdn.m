function [F,v_g,v_z,v_int,v_txdn,v_o] = fba_txdn(model_1)
rxnIDs = findRxnIDs(model_1,{'EX_glc__D_e','EX_txdn_c','EX_o2_e','EX_g3p_e',...
    'BIOMASS_Ec_iJO1366_core_53p95M','EX_xyl__D_e'});
solution = optimizeCbModel(model_1);

if ~isempty(solution.x)
    flux = solution.x;
    v_g = flux(rxnIDs(1));
    v_txdn = flux(rxnIDs(2));
    v_o = flux(rxnIDs(3));
    v_int = flux(rxnIDs(4));
    F = flux(rxnIDs(5));
    v_z = flux(rxnIDs(6));

else
    F = 0;
    v_g = 0;
    v_o = 0;
    v_txdn = 0;
    v_int = 0;
    v_z = 0;
end
end
