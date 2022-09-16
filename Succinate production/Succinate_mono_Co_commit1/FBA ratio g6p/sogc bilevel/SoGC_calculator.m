function  [w,uu] = SoGC_calculator(model, targetRxn, biomassRxn)
% Strength of Growth Coupling
% -------------------------------------------------------------------------
model1 = model;
model1 = changeObjective(model1, biomassRxn, 1);
r1 = optimizeCbModel(model1); mu_max = r1.f; 

model2 = model;
model2 = changeObjective(model2, targetRxn, 1);
r2 = optimizeCbModel(model2); v_max = r2.f; 
% -------------------------------------------------------------------------

nPts = 1000;

[biomassValues, targetValues] = productionEnvelope_no_plot(model, [], 'b', targetRxn, biomassRxn, [], nPts);
targetLowerBound = targetValues(:,1);
biomass_vector = biomassValues;
clf
% -------------------------------------------------------------------------

first_point = targetLowerBound(1);

if (all(targetLowerBound == first_point))
    w = NaN;
else
    
    a = find(targetLowerBound > first_point);
    
    if ~isempty(a) 
        slope =  (targetLowerBound(end) - targetLowerBound(a(1)))/(biomass_vector(end) - biomass_vector(a(1)));
        w = (targetLowerBound(end))^2/(slope);
    else
        w = NaN;
    end
    
end

uu = w/(mu_max*v_max);
end
% -------------------------------------------------------------------------