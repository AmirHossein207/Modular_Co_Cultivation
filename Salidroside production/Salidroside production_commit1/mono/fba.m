function [F,v_g,v_txdn,v_o] = fba(model_1)
rxnIDs = findRxnIDs(model_1,{'EX_glc__D_e','EX_salid_e','EX_o2_e',...
    'BIOMASS_Ec_iJO1366_core_53p95M'});

model_1 = changeObjective(model_1, {'BIOMASS_Ec_iJO1366_core_53p95M'});
solution = optimizeCbModel(model_1);

if ~isempty(solution.x)
    flux = solution.x;
    v_g = flux(rxnIDs(1));
    v_txdn = flux(rxnIDs(2));
    v_o = flux(rxnIDs(3));
    F = flux(rxnIDs(4));

    
else
    F = 0;
    v_g = 0;
    v_o = 0;
    v_txdn = 0;
end
end
