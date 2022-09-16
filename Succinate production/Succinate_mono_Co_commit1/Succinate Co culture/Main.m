clear

clc
global model_up model_down Mw
SolverOk = changeCobraSolver('gurobi','lp');
load('iJO1366.mat');
model = iJO1366;
Mw = [0.18 0.06 0.168 0.150];


%% intermediate strain
model_up = model;

% model_up = addReaction(model_up, 'PEPtex', 'metaboliteList', {'pep_p','pep_e'},...
%     'stoichCoeffList', [-1; 1],'lowerBound',-1000,'upperBound',1000);
% model_up = addReaction(model_up, 'PEPt6pp', 'metaboliteList', {'pep_c','pi_c','pi_p','pep_p'},...
%     'stoichCoeffList', [-1; 1; -1; 1],'lowerBound',-1000,'upperBound',1000);
% [model_up, EX_pep_e] = addExchangeRxn(model_up, 'pep_e',0,1000);
% 
% 
% rxnIDs = findRxnIDs(model_up,{'PEPt6pp','DDPA','KDOPS','PPC','PSCVT',...
%     'UAGCVT','DHAPT','GLCptspp'});
% 
% S = model_up.S;
% b1 = 0.99;
% ss = zeros(1,length(model_up.rxns));
% ss(1,rxnIDs(1:8)) = [1-b1 -b1 -b1 -b1 -b1 -b1 -b1 -b1];
% S(end+1,:) = ss;
% 
% model_up.S = S;
% 
% model_up = addMetabolite(model_up, 'iodine1_c');




model_up = removeRxns(model_up,'G6Pt6_2pp');
model_up = addReaction(model_up, 'G6Pt6_2pp_reverse', 'metaboliteList', {'pi_c',...
    'g6p_p','g6p_c','pi_p'},'stoichCoeffList', [2;1;-1;-2],'lowerBound',...
    0,'upperBound',1000);

model_up = removeRxns(model_up,'PGMT');
model_up = addReaction(model_up, 'PGMT', 'metaboliteList', {'g1p_c','g6p_c'}...
    ,'stoichCoeffList', [1;-1],'lowerBound',...
    -1000,'upperBound',1000);


b1 = 0.99;
b2 = 0.99;
b3 = 0.99;
b4 = 0.99;
b5 = 0.99;
b6 = 0.99;
% model_up = deleteModelGenes(model_up,...
%     {'b4151', 'b4014', 'b2976', 'b0529', 'b2913'});
% model_up = deleteModelGenes(model_up,...
%     {'b4015', 'b2913', 'b3952', 'b3114', 'b0902', 'b3951', 'b2579', 'b0903', 'b1897'});
rxnIDs = findRxnIDs(model_up,{'G6Pt6_2pp_reverse','G6PDH2r',...
    'G6Pt6_2pp_reverse','PGI',...
    'G6Pt6_2pp_reverse','PGMT',...
    'G6Pt6_2pp_reverse','G6PP',...
    'GLCptspp','DHAPT','UAGCVT','PSCVT','PPC','KDOPS','DDPA','PYK',...
    'GLCptspp','GLCt2pp'});


S = model_up.S;



ss = zeros(1,length(model_up.rxns));
ss(1,rxnIDs(1:2)) = [1-b1 -b1];
S(end+1,:) = ss;

ss = zeros(1,length(model_up.rxns));
ss(1,rxnIDs(3:4)) = [1-b2 -b2];
S(end+1,:) = ss;

ss = zeros(1,length(model_up.rxns));
ss(1,rxnIDs(5:6)) = [1-b3 -b3];
S(end+1,:) = ss;

ss = zeros(1,length(model_up.rxns));
ss(1,rxnIDs(7:8)) = [1-b4 -b4];
S(end+1,:) = ss;

ss = zeros(1,length(model_up.rxns));
ss(1,rxnIDs(9:16)) = [1-b5 -b5 -b5 -b5 -b5 -b5 -b5 -b5];
S(end+1,:) = ss;

model_up.S = S;

model_up = addMetabolite(model_up, 'iodine1_c');
model_up = addMetabolite(model_up, 'iodine2_c');
model_up = addMetabolite(model_up, 'iodine3_c');
model_up = addMetabolite(model_up, 'iodine4_c');
model_up = addMetabolite(model_up, 'iodine5_c');
model_up = changeObjective(model_up,'BIOMASS_Ec_iJO1366_core_53p95M');
model_up = changeRxnBounds(model_up, {'EX_glc__D_e','EX_o2_e',...
    'EX_xyl__D_e'},[-10 -15 0],{'l','l','l'});
solution_1 = optimizeCbModel(model_up);


printFluxVector(model_up, solution_1.x, true, true,-1,[],[],true)
fprintf('\n')

%% production strain
model_down = model;
% model_down = addReaction(model_down, 'PEPtex', 'metaboliteList', {'pep_p','pep_e'},...
%     'stoichCoeffList', [-1; 1],'lowerBound',-1000,'upperBound',1000);
% model_down = addReaction(model_down, 'PEPt6pp', 'metaboliteList', {'pep_c','pi_c','pi_p','pep_p'},...
%     'stoichCoeffList', [-1; 1; -1; 1],'lowerBound',-1000,'upperBound',1000);
% [model_down, EX_pep_e] = addExchangeRxn(model_down, 'pep_e',-20,0);
b1 = 0.99;
b2 = 0.99;
b3 = 0.99;
rxnIDs = findRxnIDs(model_down,{'ICL','ICDHyr',...
    'SUCCt3pp','SUCDi',...
    'FRD2','FUM'});

S = model_down.S;

ss = zeros(1,length(model_down.rxns));
ss(1,rxnIDs(1:2)) = [1-b1 -b1];
S(end+1,:) = ss;
ss = zeros(1,length(model_down.rxns));
ss(1,rxnIDs(3:4)) = [1-b2 -b2];
S(end+1,:) = ss;
ss = zeros(1,length(model_down.rxns));
ss(1,rxnIDs(5:6)) = [1-b3 -b3];
S(end+1,:) = ss;



model_down.S = S;

model_down = addMetabolite(model_down, 'iodine1_c');
model_down = addMetabolite(model_down, 'iodine2_c');
model_down = addMetabolite(model_down, 'iodine3_c');
% model_down = deleteModelGenes(model_down,{'b3115', 'b2296', 'b1849', 'b1612', 'b4122', 'b1611', 'b2029', 'b4388'});
model_down = changeRxnBounds(model_down,{'EX_glc__D_e','EX_o2_e','EX_pep_e','EX_xyl__D_e'}...
    ,[0 -20 0 -10],{'l','l','b','l'});
solution = optimizeCbModel(model_down);
printFluxVector(model_down, solution.x, true, true,-1,[],[],true)
fprintf('\n')
%% solve dynamic model
t0 = 0;
tEnd = 10;
%---------------- g/l mmol/l mmol/l
time_interval = [0 30];
% x0_up = 0.05:0.05:1;
% for i=1:length(x0_up)
initial_values = [0.05;0.05;15;15;0;0];
[t,y] = ode15s(@dynamic_ecoli,time_interval,initial_values);
%     x_up(i) = y(end,1);
%     x_do(i) = y(end,2);
%     suc(i) = y(end,4);
% end
figure
plot(t,y(:,1),'g')
hold on
plot(t,y(:,2),'r')
legend('biomass up','biomass down')
ylabel('g/l')
xlabel('h')
figure
plot(t,y(:,3),'c')
hold on
plot(t,y(:,4),'b')
hold on
plot(t,y(:,5),'r')
hold on
plot(t,y(:,6),'g')
legend on
ylabel('g/l')
xlabel('h')
legend('glucose','xylose','succ','intermediate','oxygen')

titer = y(end,5);
volume = 1;
tss = 6;
tend = 27.2053916296449;
Productivity = volume.*titer./(tend+tss);
Yield = volume.*titer./volume./20;
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



