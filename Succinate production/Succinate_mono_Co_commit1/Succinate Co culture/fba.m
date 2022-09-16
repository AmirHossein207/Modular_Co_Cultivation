function [F,v_g,v_int,v_z,v_succ,v_o] = fba(model_1)
rxnIDs = findRxnIDs(model_1,{'EX_glc__D_e','EX_xyl__D_e',...
    'EX_succ_e','EX_o2_e','EX_g6p_e',...
    'BIOMASS_Ec_iJO1366_core_53p95M'});
solution = optimizeCbModel(model_1);

if ~isempty(solution.x)
    flux = solution.x;
    v_g = flux(rxnIDs(1));
    v_z = flux(rxnIDs(2));
    v_succ = flux(rxnIDs(3));
    v_o = flux(rxnIDs(4));
    v_int = flux(rxnIDs(5));
    F = flux(rxnIDs(6));


else
    F = 0;
    v_g = 0;
    v_z = 0;
    v_o = 0;
    v_succ = 0;
    v_int = 0;
end
end
