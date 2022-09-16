function obj = e_objfun(x,~)
global  model

b1 = x(1,1); 
b2 = x(1,2); 
b3 = x(1,3); 




model1 = model;

rxnIDs = findRxnIDs(model1,{'ICL','ICDHyr',...
    'SUCCt3pp','SUCDi',...
    'FRD2','FUM'});

S = model1.S;

ss = zeros(1,length(model1.rxns));
ss(1,rxnIDs(1:2)) = [1-b1 -b1];
S(end+1,:) = ss;
ss = zeros(1,length(model1.rxns));
ss(1,rxnIDs(3:4)) = [1-b2 -b2];
S(end+1,:) = ss;
ss = zeros(1,length(model1.rxns));
ss(1,rxnIDs(5:6)) = [1-b3 -b3];
S(end+1,:) = ss;




model1.S = S;

model1 = addMetabolite(model1, 'iodine1_c');
model1 = addMetabolite(model1, 'iodine2_c');
model1 = addMetabolite(model1, 'iodine3_c');


model1 = changeObjective(model1,'BIOMASS_Ec_iJO1366_core_53p95M');
model1 = changeRxnBounds(model1, {'EX_glc__D_e','EX_o2_e'},[-10, -15],{'l','l'});
[w2,~] = SoGC_calculator(model1,{'EX_g6p_e'},{'BIOMASS_Ec_iJO1366_core_53p95M'});
obj = -w2;

end