clear
clc
close all
solverOK = changeCobraSolver('gurobi','lp');
load('iJO1366.mat');
model = iJO1366;

model = addReaction(model,'TYRCBOX', 'metaboliteList', {'h_c','tyr__L_c',...
    'co2_c','tym_c'}, 'stoichCoeffList', [-1;-1;1;1],...
    'lowerBound',0,'upperBound',1000);
model = addReaction(model,'TYROXDAc', 'metaboliteList', {'h2o_c','o2_c',...
    'tym_c','h2o2_c','nh4_c','4hoxpacd_c'}, 'stoichCoeffList', [-1;-1;-1;1;1;1],...
    'lowerBound',0,'upperBound',1000);
%%%%% add metabolite 4-tyrosol
model = addReaction(model,'TYROSOL', 'metaboliteList', {'4hoxpacd_c','h_c',...
    'nadh_c','nad_c','4_tyrosol_c'}, 'stoichCoeffList', [-1;-1;-1;1;1],...
    'lowerBound',0,'upperBound',1000);
%%%%% add metabolite salidroside
model = addReaction(model,'SALID', 'metaboliteList', {'udpg_c','4_tyrosol_c',...
    'udp_c','h_c','salid_c'}, 'stoichCoeffList', [-1;-1;1;1;1],...
    'lowerBound',0,'upperBound',1000);
model = addReaction(model,'SALIDpp', 'metaboliteList', {'salid_c','salid_e'},...
    'stoichCoeffList', [-1;1;], 'lowerBound',0,'upperBound',1000);
[model, EX_salid_e] = addExchangeRxn(model, 'salid_e',0,1000);


model = removeRxns(model,'TYRt2rpp');
model = addReaction(model, 'TYRt2rpp', 'metaboliteList', {'h_p','tyr__L_p', ...
    'tyr__L_c','h_c'},'stoichCoeffList', [1;1;-1;-1],'lowerBound',...
    -1000,'upperBound',1000);
model = removeRxns(model,'TYRTA');
model = addReaction(model, 'TYRTA', 'metaboliteList', {'akg_c','tyr__L_c',...
    '34hpp_c','glu__L_c'},'stoichCoeffList', [1;1;-1;-1],'lowerBound',...
    -1000,'upperBound',1000);
% model = removeRxns(model,'TALA');
% model = addReaction(model, 'TALA', 'metaboliteList', {'g3p_c','s7p_c',...
%     'e4p_c','f6p_c'},'stoichCoeffList', [1;1;-1;-1],'lowerBound',...
%     -1000,'upperBound',1000);
model = removeRxns(model,'ASPTA');
model = addReaction(model, 'ASPTA', 'metaboliteList', {'akg_c','asp__L_c',...
    'glu__L_c','oaa_c'},'stoichCoeffList', [1;1;-1;-1],'lowerBound',...
    -1000,'upperBound',1000);
model = removeRxns(model,'FBA');
model = addReaction(model, 'FBA', 'metaboliteList', {'fdp_c','dhap_c',...
    'g3p_c'},'stoichCoeffList', [1;-1;-1],'lowerBound',...
    -1000,'upperBound',1000);

rxnIDs = findRxnIDs(model,{'PSCVT','DHAPT','PSCVT','KDOPS','PSCVT','PPC','PSCVT',...
    'PYK','TYRTA','ASPTA','PSERT','TYRCBOX','TYRt2rpp'});
b1 = 0.99;
b2 = 0.60;
b3 = 0.99;

S = model.S;

ss = zeros(1,length(model.rxns));
ss(1,rxnIDs(1:2)) = [1-b2 -b2];
S(end+1,:) = ss;
ss = zeros(1,length(model.rxns));
ss(1,rxnIDs(3:4)) = [1-b1 -b1];
S(end+1,:) = ss;
ss = zeros(1,length(model.rxns));
ss(1,rxnIDs(5:6)) = [1-b1 -b1];
S(end+1,:) = ss;
ss = zeros(1,length(model.rxns));
ss(1,rxnIDs(7:8)) = [1-b1 -b1];
S(end+1,:) = ss;
ss = zeros(1,length(model.rxns));
ss(1,rxnIDs(9:11)) = [1-b1 -b1 -b1];
S(end+1,:) = ss;
ss = zeros(1,length(model.rxns));
ss(1,rxnIDs(12:13)) = [1-b3 -b3];
S(end+1,:) = ss;
model.S = S;

model = addMetabolite(model, 'iodine1_c');
model = addMetabolite(model, 'iodine2_c');
model = addMetabolite(model, 'iodine3_c');
model = addMetabolite(model, 'iodine4_c');
model = addMetabolite(model, 'iodine5_c');
model = addMetabolite(model, 'iodine6_c');
model = changeRxnBounds(model, {'EX_glc__D_e','EX_xyl__D_e','EX_o2_e'...
    },[-10,0,-15],{'l','l','l'});

% model = removeRxns(model,{'GLUDy','MOX'});
model = changeObjective(model,'BIOMASS_Ec_iJO1366_core_53p95M');
solution1 = optimizeCbModel(model);
printFluxVector(model, solution1.x, true, true,-1,[],[],true)
fprintf('\n')
[directionRxns_W, involvedMets_W, deadEnds_W] = draw_by_met(model, {'tyr__L_c'},...
    'true', 1, 'prod', {''}, solution1.x);












