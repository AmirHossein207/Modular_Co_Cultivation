function obj = e_objfun(x,~)
global  model

b1 = x(1,1); 
b2 = x(1,2); 


model1 = model;


rxnIDs = findRxnIDs(model1,{'G3P_C','GAPD','G3P_C','TALA','G3P_C',...
    'TPI','TALA','PGI','G6PDH2r'});


S = model1.S;


ss = zeros(1,length(model1.rxns));
ss(1,rxnIDs(1:2)) = [1-b1 -b1];
S(end+1,:) = ss;

ss = zeros(1,length(model1.rxns));
ss(1,rxnIDs(5:6)) = [1-b2 -b2];
S(end+1,:) = ss;

model1.S = S;
model1 = changeObjective(model1,'BIOMASS_Ec_iJO1366_core_53p95M');
model1 = changeRxnBounds(model1, {'EX_glc__D_e','EX_o2_e'},[-10, -15],'l');
[w2,~] = SoGC_calculator(model1,{'EX_g3p_e'},{'BIOMASS_Ec_iJO1366_core_53p95M'});
obj = -w2;

end