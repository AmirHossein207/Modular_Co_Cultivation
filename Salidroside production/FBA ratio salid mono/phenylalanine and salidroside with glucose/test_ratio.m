clear
clc
close all
solverOK = changeCobraSolver('gurobi','lp');
load('iJO1366.mat');
model = iJO1366;


model = addReaction(model,'TYROSpp', 'metaboliteList', {'4_tyrosol_c','4_tyrosol_e'},...
    'stoichCoeffList', [-1;1], 'lowerBound',-1000,'upperBound',0);
[model, EX_4_tyrosol_e] = addExchangeRxn(model, '4_tyrosol_e',-1000,0);
model = addReaction(model,'SALID', 'metaboliteList', {'udpg_c','4_tyrosol_c',...
    'udp_c','h_c','salid_c'}, 'stoichCoeffList', [-1;-1;1;1;1],...
    'lowerBound',0,'upperBound',1000);
model = addReaction(model,'SALIDpp', 'metaboliteList', {'salid_c','salid_e'},...
    'stoichCoeffList', [-1;1;], 'lowerBound',0,'upperBound',1000);
[model, EX_salid_e] = addExchangeRxn(model, 'salid_e',0,1000);

model = removeRxns(model,'PGMT');
model = addReaction(model, 'PGMT', 'metaboliteList', {'g1p_c','g6p_c'}...
   ,'stoichCoeffList', [1;-1],'lowerBound',-1000,'upperBound',1000);
model = removeRxns(model,'MLTP2');
model = addReaction(model, 'MLTP2', 'metaboliteList', {'malthx_c','pi_c',...
    'g1p_c','maltpt_c'},'stoichCoeffList', [1;1;-1;-1],'lowerBound',...
    -1000,'upperBound',1000);
model = removeRxns(model,'MLTP1');
model = addReaction(model, 'MLTP1', 'metaboliteList', {'maltpt_c','pi_c',...
    'g1p_c','maltttr_c'},'stoichCoeffList', [1;1;-1;-1],'lowerBound',...
    -1000,'upperBound',1000);
model = removeRxns(model,'MLTP3');
model = addReaction(model, 'MLTP3', 'metaboliteList', {'malthp_c','pi_c',...
    'g1p_c','malthx_c'},'stoichCoeffList', [1;1;-1;-1],'lowerBound',...
    -1000,'upperBound',1000);
model = removeRxns(model,'PHEt2rpp');
model = addReaction(model, 'PHEt2rpp', 'metaboliteList', {'h_p','phe__L_p',...
    'h_c','phe__L_c'},'stoichCoeffList', [1;1;-1;-1],'lowerBound',...
    -1000,'upperBound',1000);
model = removeRxns(model,'TYRTA');
model = addReaction(model, 'TYRTA', 'metaboliteList', {'akg_c','tyr__L_c',...
    '34hpp_c','glu__L_c'},'stoichCoeffList', [-1;-1;1;1],'lowerBound',...
    0,'upperBound',1000);
rxnIDs = findRxnIDs(model,{'PGMT','PGI','PGMT','G6PDH2r','GALUi',...
    'MLTP3','MLTP1','MLTP2','SALID','TRE6PS','HEX1','XYLI2','GLCDpp',...
    'PHEt2rpp','BIOMASS_Ec_iJO1366_core_53p95M','TYRL',...
    'BIOMASS_Ec_iJO1366_core_53p95M'});
b1 = 0.80;
b2 = 0.99;
b3 = 0.01;
b4 = 0.001;
S = model.S;

ss = zeros(1,length(model.rxns));
ss(1,rxnIDs(1:2)) = [1-b1 -b1];
S(end+1,:) = ss;
ss = zeros(1,length(model.rxns));
ss(1,rxnIDs(3:4)) = [1-b1 -b1];
S(end+1,:) = ss;
ss = zeros(1,length(model.rxns));
ss(1,rxnIDs(5:8)) = [1-b2 -b2 -b2 -b2];
S(end+1,:) = ss;
ss = zeros(1,length(model.rxns));
ss(1,rxnIDs(9:10)) = [1-b2 -b2];
S(end+1,:) = ss;
ss = zeros(1,length(model.rxns));
ss(1,rxnIDs(11:13)) = [1-b2 -b2 -b2];
S(end+1,:) = ss;
ss = zeros(1,length(model.rxns));
ss(1,rxnIDs(14:15)) = [1-b3 -b3];
S(end+1,:) = ss;
% ss = zeros(1,length(model.rxns));
% ss(1,rxnIDs(16:17)) = [1-b4 -b4];
% S(end+1,:) = ss;
model.S = S;

model = addMetabolite(model, 'iodine1_c');
model = addMetabolite(model, 'iodine2_c');
model = addMetabolite(model, 'iodine3_c');
model = addMetabolite(model, 'iodine4_c');
model = addMetabolite(model, 'iodine5_c');
model = addMetabolite(model, 'iodine6_c');
% model = addMetabolite(model, 'iodine7_c');
model = changeRxnBounds(model,{'EX_glc__D_e','EX_xyl__D_e','EX_o2_e'...
    ,'EX_tyr__L_e','EX_4_tyrosol_e','EX_phe__L_e'},...
    [-10,0,-15,-1,-5,5]...
    ,{'l','l','l','l','l','u'});

% model = deleteModelGenes(model,{'b2600'});
% model = removeRxns(model,{'TYRTA'});
model = changeObjective(model,'BIOMASS_Ec_iJO1366_core_53p95M');
solution1 = optimizeCbModel(model);
printFluxVector(model, solution1.x, true, true,-1,[],[],false)
fprintf('\n')
[directionRxns_W, involvedMets_W, deadEnds_W] = draw_by_met(model,...
    {'tyr__L_c'},...
    'true', 1, 'prod', {''}, solution1.x);

% [minFlux, maxFlux] = fluxVariability(model,100, 'max', ...
%     {'EX_phe__L_e','EX_tyr__L_e'})










