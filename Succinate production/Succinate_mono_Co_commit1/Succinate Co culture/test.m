clear

clc

SolverOk = changeCobraSolver('glpk','lp');
load('iJO1366.mat');
model = iJO1366;
model = addReaction(model, 'PEPtex', 'metaboliteList', {'pep_p','pep_e'},...
    'stoichCoeffList', [-1; 1],'lowerBound',-1000,'upperBound',1000);
model = addReaction(model, 'PEPt6pp', 'metaboliteList', {'pep_c','pi_c','pi_p','pep_p'},...
    'stoichCoeffList', [-1; 1; -1; 1],'lowerBound',-1000,'upperBound',1000);
[model, EX_pep_e] = addExchangeRxn(model, 'pep_e',-100,1000);
model1 = model;
model_up = model;
rxnIDs = findRxnIDs(model_up,{'PEPt6pp','GLCptspp','UAGCVT','PSCVT','PPC',...
    'KDOPS','DDPA','PGM','PGCD','GAPD','TALA','DXPS',...
    'PFK','F6PA','GF6PTA','GLCt2pp','GLCptspp'});

b0 = 1.050000;
b1 = 3.940000;
b2 = 1.210000;
b3 = 2.150000;
b4 = 3.000000;

S = model_up.S;


ss = zeros(1,length(model_up.rxns));
ss(1,rxnIDs(1:7)) = [1-b0 ones(1,6)*(b0)];
S(end+1,:) = ss;
% ss = zeros(1,length(model_up.rxns));
% ss(1,rxnIDs(8:9)) = [1-b1 -b1];
% S(end+1,:) = ss;
% ss = zeros(1,length(model_up.rxns));
% ss(1,rxnIDs(10:12)) = [1-b2 -b2 -b2];
% S(end+1,:) = ss;
% ss = zeros(1,length(model_up.rxns));
% ss(1,rxnIDs(13:15)) = [1-b3 -b3 -b3];
% S(end+1,:) = ss;
% ss = zeros(1,length(model_up.rxns));
% ss(1,rxnIDs(16:17)) = [1-b4 -b4];
% S(end+1,:) = ss;
model_up.S = S;

model_up = addMetabolite(model_up, 'iodine1_c');
% model_up = addMetabolite(model_up, 'iodine2_c');
% model_up = addMetabolite(model_up, 'iodine3_c');
% model_up = addMetabolite(model_up, 'iodine4_c');
% model_up = addMetabolite(model_up, 'iodine5_c');

model1 = changeRxnBounds(model1,{'EX_glc__D_e','EX_o2_e','EX_pep_e'}...
    ,[-30 -20 0],{'l','l','b'});

model1 = changeObjective(model1,'EX_succ_e');
solution1 = optimizeCbModel(model1);
printFluxVector(model1, solution1.x, true, true,-1,[],[],true)
fprintf('\n')
model2 = model;
model2 = changeRxnBounds(model2,{'EX_glc__D_e','EX_o2_e','EX_pep_e'}...
    ,[-15 -20 -13.40],{'l','l','b'});
model2 = changeObjective(model2,'EX_succ_e');
solution2 = optimizeCbModel(model2);
printFluxVector(model2, solution2.x, true, true,-1,[],[],true)