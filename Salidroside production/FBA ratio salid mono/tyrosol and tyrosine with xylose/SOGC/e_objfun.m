function obj = e_objfun(x,~)
global  model

b1 = x(1,1); 
b2 = x(1,2); 




model1 = model;



rxnIDs = findRxnIDs(model1,{'PSCVT','DHAPT','PSCVT','KDOPS','PSCVT','PPC','PSCVT',...
    'PYK','TYRTA','ASPTA','PSERT','TYRCBOX','TYRt2rpp'});


S = model1.S;


ss = zeros(1,length(model1.rxns));
ss(1,rxnIDs(9:11)) = [1-b1 -b1 -b1];
S(end+1,:) = ss;
ss = zeros(1,length(model1.rxns));
ss(1,rxnIDs(12:13)) = [1-b2 -b2];
S(end+1,:) = ss;
model1.S = S;

model1 = addMetabolite(model1, 'iodine1_c');
model1 = addMetabolite(model1, 'iodine2_c');

model1 = changeRxnBounds(model1, {'EX_glc__D_e','EX_xyl__D_e','EX_o2_e'...
    ,'EX_phe__L_e','EX_tyr__L_e'},[0,-10,-15,-5,0],{'l','l','l','l','l'});
model1 = changeObjective(model1,'BIOMASS_Ec_iJO1366_core_53p95M');
[w2,~] = SoGC_calculator(model1,{'EX_4_tyrosol_e'},{'BIOMASS_Ec_iJO1366_core_53p95M'});
obj = -w2;

end