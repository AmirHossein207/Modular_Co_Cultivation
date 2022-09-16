clear
clc

global model_up model_down Mw

SolverOk = changeCobraSolver('gurobi','lp');
load('iJO1366.mat');
model = iJO1366;

Mw = [0.18 0.3 0.138 0.165 0.15 0.138];


%% intermediate strain
model_up = model;
model_up = addReaction(model_up,'TYRCBOX', 'metaboliteList', {'h_c','tyr__L_c',...
    'co2_c','tym_c'}, 'stoichCoeffList', [-1;-1;1;1],...
    'lowerBound',0,'upperBound',1000);
model_up = addReaction(model_up,'TYROXDAc', 'metaboliteList', {'h2o_c','o2_c',...
    'tym_c','h2o2_c','nh4_c','4hoxpacd_c'}, 'stoichCoeffList', [-1;-1;-1;1;1;1],...
    'lowerBound',0,'upperBound',1000);
%%%%% add metabolite 4-tyrosol
model_up = addReaction(model_up,'TYROSOL', 'metaboliteList', {'4hoxpacd_c','h_c',...
    'nadh_c','nad_c','4_tyrosol_c'}, 'stoichCoeffList', [-1;-1;-1;1;1],...
    'lowerBound',0,'upperBound',1000);
model_up = addReaction(model_up,'TYROSpp', 'metaboliteList', {'4_tyrosol_c','4_tyrosol_e'},...
    'stoichCoeffList', [-1;1;], 'lowerBound',0,'upperBound',1000);
[model_up, EX_4_tyrosol_e] = addExchangeRxn(model_up, '4_tyrosol_e',0,1000);


model_up = removeRxns(model_up,'TYRt2rpp');
model_up = addReaction(model_up, 'TYRt2rpp', 'metaboliteList', {'h_p','tyr__L_p', ...
    'tyr__L_c','h_c'},'stoichCoeffList', [1;1;-1;-1],'lowerBound',...
    -1000,'upperBound',1000);
model_up = removeRxns(model_up,'TYRTA');
model_up = addReaction(model_up, 'TYRTA', 'metaboliteList', {'akg_c','tyr__L_c',...
    '34hpp_c','glu__L_c'},'stoichCoeffList', [1;1;-1;-1],'lowerBound',...
    -1000,'upperBound',1000);
% model = removeRxns(model,'TALA');
% model = addReaction(model, 'TALA', 'metaboliteList', {'g3p_c','s7p_c',...
%     'e4p_c','f6p_c'},'stoichCoeffList', [1;1;-1;-1],'lowerBound',...
%     -1000,'upperBound',1000);
model_up = removeRxns(model_up,'ASPTA');
model_up = addReaction(model_up, 'ASPTA', 'metaboliteList', {'akg_c','asp__L_c',...
    'glu__L_c','oaa_c'},'stoichCoeffList', [1;1;-1;-1],'lowerBound',...
    -1000,'upperBound',1000);
model_up = removeRxns(model_up,'FBA');
model_up = addReaction(model_up, 'FBA', 'metaboliteList', {'fdp_c','dhap_c',...
    'g3p_c'},'stoichCoeffList', [1;-1;-1],'lowerBound',...
    -1000,'upperBound',1000);

rxnIDs = findRxnIDs(model_up,{'PSCVT','DHAPT','PSCVT','KDOPS','PSCVT','PPC','PSCVT',...
    'PYK','TYRTA','ASPTA','PSERT','TYRCBOX','TYRt2rpp'});
b1 = 0.99;
b2 = 0.60;
b3 = 0.99;

S = model_up.S;

ss = zeros(1,length(model_up.rxns));
ss(1,rxnIDs(1:2)) = [1-b2 -b2];
S(end+1,:) = ss;
ss = zeros(1,length(model_up.rxns));
ss(1,rxnIDs(3:4)) = [1-b1 -b1];
S(end+1,:) = ss;
ss = zeros(1,length(model_up.rxns));
ss(1,rxnIDs(5:6)) = [1-b1 -b1];
S(end+1,:) = ss;
ss = zeros(1,length(model_up.rxns));
ss(1,rxnIDs(7:8)) = [1-b1 -b1];
S(end+1,:) = ss;
ss = zeros(1,length(model_up.rxns));
ss(1,rxnIDs(9:11)) = [1-b1 -b1 -b1];
S(end+1,:) = ss;
ss = zeros(1,length(model_up.rxns));
ss(1,rxnIDs(12:13)) = [1-b3 -b3];
S(end+1,:) = ss;
model_up.S = S;

model_up = addMetabolite(model_up, 'iodine1_c');
model_up = addMetabolite(model_up, 'iodine2_c');
model_up = addMetabolite(model_up, 'iodine3_c');
model_up = addMetabolite(model_up, 'iodine4_c');
model_up = addMetabolite(model_up, 'iodine5_c');
model_up = addMetabolite(model_up, 'iodine6_c');
model_up = changeRxnBounds(model_up, {'EX_glc__D_e','EX_xyl__D_e','EX_o2_e'...
    ,'EX_phe__L_e'},[0,-10,-15,-5],{'l','l','l','l'});

% model = removeRxns(model,{'GLUDy','MOX'});
model_up = changeObjective(model_up,'BIOMASS_Ec_iJO1366_core_53p95M');
solution1 = optimizeCbModel(model_up);
printFluxVector(model_up, solution1.x, true, true,-1,[],[],true)
fprintf('\n')
%% production strain
model_down = model;
model_down = addReaction(model_down,'TYROSpp', 'metaboliteList', {'4_tyrosol_c','4_tyrosol_e'},...
    'stoichCoeffList', [-1;1], 'lowerBound',-1000,'upperBound',0);
[model_down, EX_4_tyrosol_e] = addExchangeRxn(model_down, '4_tyrosol_e',-1000,0);
model_down = addReaction(model_down,'SALID', 'metaboliteList', {'udpg_c','4_tyrosol_c',...
    'udp_c','h_c','salid_c'}, 'stoichCoeffList', [-1;-1;1;1;1],...
    'lowerBound',0,'upperBound',1000);
model_down = addReaction(model_down,'SALIDpp', 'metaboliteList', {'salid_c','salid_e'},...
    'stoichCoeffList', [-1;1;], 'lowerBound',0,'upperBound',1000);
[model_down, EX_salid_e] = addExchangeRxn(model_down, 'salid_e',0,1000);

model_down = removeRxns(model_down,'PGMT');
model_down = addReaction(model_down, 'PGMT', 'metaboliteList', {'g1p_c','g6p_c'}...
   ,'stoichCoeffList', [1;-1],'lowerBound',-1000,'upperBound',1000);
model_down = removeRxns(model_down,'MLTP2');
model_down = addReaction(model_down, 'MLTP2', 'metaboliteList', {'malthx_c','pi_c',...
    'g1p_c','maltpt_c'},'stoichCoeffList', [1;1;-1;-1],'lowerBound',...
    -1000,'upperBound',1000);
model_down = removeRxns(model_down,'MLTP1');
model_down = addReaction(model_down, 'MLTP1', 'metaboliteList', {'maltpt_c','pi_c',...
    'g1p_c','maltttr_c'},'stoichCoeffList', [1;1;-1;-1],'lowerBound',...
    -1000,'upperBound',1000);
model_down = removeRxns(model_down,'MLTP3');
model_down = addReaction(model_down, 'MLTP3', 'metaboliteList', {'malthp_c','pi_c',...
    'g1p_c','malthx_c'},'stoichCoeffList', [1;1;-1;-1],'lowerBound',...
    -1000,'upperBound',1000);
model_down = removeRxns(model_down,'PHEt2rpp');
model_down = addReaction(model_down, 'PHEt2rpp', 'metaboliteList', {'h_p','phe__L_p',...
    'h_c','phe__L_c'},'stoichCoeffList', [1;1;-1;-1],'lowerBound',...
    -1000,'upperBound',1000);
% model = removeRxns(model,'ASPTA');
% model = addReaction(model, 'ASPTA', 'metaboliteList', {'akg_c','asp__L_c',...
%     'glu__L_c','oaa_c'},'stoichCoeffList', [1;1;-1;-1],'lowerBound',...
%     -1000,'upperBound',1000);
% 
rxnIDs = findRxnIDs(model_down,{'PGMT','PGI','PGMT','G6PDH2r','GALUi',...
    'MLTP3','MLTP1','MLTP2','SALID','TRE6PS','HEX1','XYLI2','GLCDpp',...
    'PHEt2rpp','BIOMASS_Ec_iJO1366_core_53p95M',''});
b1 = 0.80;
b2 = 0.99;
b3 = 0.1;

S = model_down.S;

ss = zeros(1,length(model_down.rxns));
ss(1,rxnIDs(1:2)) = [1-b1 -b1];
S(end+1,:) = ss;
ss = zeros(1,length(model_down.rxns));
ss(1,rxnIDs(3:4)) = [1-b1 -b1];
S(end+1,:) = ss;
ss = zeros(1,length(model_down.rxns));
ss(1,rxnIDs(5:8)) = [1-b2 -b2 -b2 -b2];
S(end+1,:) = ss;
ss = zeros(1,length(model_down.rxns));
ss(1,rxnIDs(9:10)) = [1-b2 -b2];
S(end+1,:) = ss;
ss = zeros(1,length(model_down.rxns));
ss(1,rxnIDs(11:13)) = [1-b2 -b2 -b2];
S(end+1,:) = ss;
ss = zeros(1,length(model_down.rxns));
ss(1,rxnIDs(14:15)) = [1-b3 -b3];
S(end+1,:) = ss;
model_down.S = S;

model_down = addMetabolite(model_down, 'iodine1_c');
model_down = addMetabolite(model_down, 'iodine2_c');
model_down = addMetabolite(model_down, 'iodine3_c');
model_down = addMetabolite(model_down, 'iodine4_c');
model_down = addMetabolite(model_down, 'iodine5_c');
model_down = addMetabolite(model_down, 'iodine6_c');
model_down = changeRxnBounds(model_down, {'EX_glc__D_e','EX_xyl__D_e','EX_o2_e'...
    ,'EX_tyr__L_e','EX_4_tyrosol_e','EX_phe__L_e','EX_phe__L_e'},[-10,0,-15,-5,-5,5,0]...
    ,{'l','l','l','l','l','u','l'});

% model_down = deleteModelGenes(model_down,{'b2600'});
model_down = changeObjective(model_down,'BIOMASS_Ec_iJO1366_core_53p95M');
solution1 = optimizeCbModel(model_down);
printFluxVector(model_down, solution1.x, true, true,-1,[],[],true)
fprintf('\n')
%% solve dynamic model
t0 = 0;
tEnd = 10;
%---------------- g/l mmol/l mmol/l
time_interval = [0 180];
% x0_up = 0.05:0.05:1;
% for i=1:length(x0_up)
    initial_values = [0.05;0.05;15;15;0;0;0;0];
    [t,y] = ode23s(@dynamic_ecoli,time_interval,initial_values);
%     x_up(i) = y(end,1);
%     x_do(i) = y(end,2);
%     txdn(i) = y(end,4);
%     inter(i) = y(end,5);
% end
figure
plot(t,y(:,1),'g')
hold on
plot(t,y(:,2),'r')
legend('biomass up','biomass down')
ylabel('g/l')
xlabel('h')
figure
plot(t,y(:,3),'b')
hold on
plot(t,y(:,4),'c')
hold on
plot(t,y(:,5),'g')
hold on
plot(t,y(:,6),'r')
legend on
ylabel('g/l')
xlabel('h')
legend('glucose','xylose','salidroside','tyrosol')
figure
plot(t,y(:,7),'g')
hold on
plot(t,y(:,8),'r')
legend('phenylalanin','tyrosine')
ylabel('g/l')
xlabel('h')
% titer = y(end,5);
% volume = 1;
% tss = 6;
% tend =71;
% Productivity = volume.*titer./(tend+tss);
% Yield = volume.*titer./volume./20;
%% write table
% Biomass_upstream = y(:,1);
% Biomass_Downstream = y(:,2);
% Glucose = y(:,3);
% Intermediate = y(:,5);
% Succinate = y(:,4);
% Oxygen = y(:,6);
% T = table(t,Biomass_upstream,Biomass_Downstream,Glucose,Succinate,...
%     Intermediate,Oxygen);
% writetable(T,'Succinate_Co_culture_results.xls')



