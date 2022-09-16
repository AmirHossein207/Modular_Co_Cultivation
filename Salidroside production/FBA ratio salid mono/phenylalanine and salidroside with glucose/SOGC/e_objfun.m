function obj = e_objfun(x,~)
global  model

b1 = x(1,1); 
b2 = x(1,2); 
b3 = x(1,3); 
b4 = x(1,4); 
b5 = x(1,5); 
b6 = x(1,6);


model1 = model;

rxnIDs = findRxnIDs(model1,{'PGMT','PGI','PGMT','G6PDH2r','GALUi',...
    'MLTP3','MLTP1','MLTP2','SALID','TRE6PS','HEX1','XYLI2','GLCDpp',...
    'PHEt2rpp','BIOMASS_Ec_iJO1366_core_53p95M','TYRL',...
    'BIOMASS_Ec_iJO1366_core_53p95M'});

S = model1.S;

ss = zeros(1,length(model1.rxns));
ss(1,rxnIDs(1:2)) = [1-b1 -b1];
S(end+1,:) = ss;
ss = zeros(1,length(model1.rxns));
ss(1,rxnIDs(3:4)) = [1-b2 -b2];
S(end+1,:) = ss;
ss = zeros(1,length(model1.rxns));
ss(1,rxnIDs(5:8)) = [1-b3 -b3 -b3 -b3];
S(end+1,:) = ss;
ss = zeros(1,length(model1.rxns));
ss(1,rxnIDs(9:10)) = [1-b4 -b4];
S(end+1,:) = ss;
ss = zeros(1,length(model1.rxns));
ss(1,rxnIDs(11:13)) = [1-b5 -b5 -b5];
S(end+1,:) = ss;
ss = zeros(1,length(model1.rxns));
ss(1,rxnIDs(14:15)) = [1-b6 -b6];
S(end+1,:) = ss;
% ss = zeros(1,length(model1.rxns));
% ss(1,rxnIDs(16:17)) = [1-b4 -b4];
% S(end+1,:) = ss;
model1.S = S;

model1 = addMetabolite(model1, 'iodine1_c');
model1 = addMetabolite(model1, 'iodine2_c');
model1 = addMetabolite(model1, 'iodine3_c');
model1 = addMetabolite(model1, 'iodine4_c');
model1 = addMetabolite(model1, 'iodine5_c');
model1 = addMetabolite(model1, 'iodine6_c');
% model1 = addMetabolite(model1, 'iodine7_c');

model1 = changeObjective(model1,'BIOMASS_Ec_iJO1366_core_53p95M');
model1 = changeRxnBounds(model1, {'EX_glc__D_e','EX_o2_e'},[-10, -15],{'l','l'});
[w2,~] = SoGC_calculator(model1,{'EX_salid_e'},{'BIOMASS_Ec_iJO1366_core_53p95M'});
obj = -w2;

end