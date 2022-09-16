function [F,v_g,v_int,v_o,vz] = fba_mev(model_1)
rxnIDs = findRxnIDs(model_1,{'EX_glc__D_e','EX_o2_e','EX_mevR_e',...
    'BIOMASS_Ec_iJO1366_core_53p95M','EX_xyl__D_e'});
solution = optimizeCbModel(model_1);

if ~isempty(solution.x)
    flux = solution.x;
    v_g = flux(rxnIDs(1));
    v_o = flux(rxnIDs(2));
    v_int = flux(rxnIDs(3));
    F = flux(rxnIDs(4));
    vz = flux(rxnIDs(5));

else
    F = 0;
    v_g = 0;
    v_o = 0;
    v_int = 0;
    vz = 0;
end
end