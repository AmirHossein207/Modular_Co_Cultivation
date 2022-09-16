function obj = e_objfun(x,~)
global  model

b1 = x(1,1); 
b2 = x(1,2); 
b3 = x(1,3); 
b4 = x(1,4); 
b5 = x(1,5); 



model1 = model;

rxnIDs = findRxnIDs(model1,{'G6Pt6_2pp_reverse','G6PDH2r',...
    'G6Pt6_2pp_reverse','PGI',...
    'G6Pt6_2pp_reverse','PGMT',...
    'G6Pt6_2pp_reverse','G6PP',...
    'GLCptspp','DHAPT','UAGCVT','PSCVT','PPC','KDOPS','DDPA','PYK'});

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

ss = zeros(1,length(model1.rxns));
ss(1,rxnIDs(7:8)) = [1-b4 -b4];
S(end+1,:) = ss;

ss = zeros(1,length(model1.rxns));
ss(1,rxnIDs(9:16)) = [1-b5 -b5 -b5 -b5 -b5 -b5 -b5 -b5];
S(end+1,:) = ss;

model1.S = S;

model1 = addMetabolite(model1, 'iodine1_c');
model1 = addMetabolite(model1, 'iodine2_c');
model1 = addMetabolite(model1, 'iodine3_c');
model1 = addMetabolite(model1, 'iodine4_c');
model1 = addMetabolite(model1, 'iodine5_c');

model1 = changeObjective(model1,'BIOMASS_Ec_iJO1366_core_53p95M');
model1 = changeRxnBounds(model1, {'EX_glc__D_e','EX_o2_e'},[-10, -15],{'l','l'});
[w2,~] = SoGC_calculator(model1,{'EX_g6p_e'},{'BIOMASS_Ec_iJO1366_core_53p95M'});
obj = -w2;

end