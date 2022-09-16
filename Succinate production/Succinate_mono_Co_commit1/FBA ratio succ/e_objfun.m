function obj = e_objfun(x,~)

a = 1; %x(1);
b = x; %x(2);

global  model pep_max bio_max rxnID

model1 = model;
model1 = changeRxnBounds(model1,{'EX_glc__D_e','EX_o2_e'}...
    ,[-10 0],{'l','l'});
model1 = changeObjective(model1,{'BIOMASS_Ec_iJO1366_core_53p95M','EX_succ_e'}...
    ,[(a/bio_max) (b/pep_max)]);

r1 = optimizeCbModel(model1);

if ~isempty(r1.x)
    p = r1.x(rxnID(2));
    bio = r1.x(rxnID(1));   
else    
    p = 0;
    bio = 0;   
end

obj = -p*bio;

end