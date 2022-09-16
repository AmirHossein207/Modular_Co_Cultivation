clear
clc

global model_up model_down Mw

SolverOk = changeCobraSolver('gurobi','lp');
load('iJO1366.mat');
model = iJO1366;

Mw = [0.18 0.06 0.272 0.170 0.15];


%% intermediate strain
model_up = model;

model_up = addReaction(model_up, 'G3P_C', 'metaboliteList', {'g3p_c',...
    'g3p_e'},'stoichCoeffList', [-1;1],'lowerBound',...
    0,'upperBound',1000);
model_up = addReaction(model_up,'EX_g3p_e', {'g3p_e'}, -1, true, 0,1000);

rxnIDs = findRxnIDs(model_up,{'G3P_C','GAPD','G3P_C','TALA','G3P_C',...
    'TPI','TALA','PGI','G6PDH2r'});
b1 = 0.99;

S = model_up.S;


ss = zeros(1,length(model_up.rxns));
ss(1,rxnIDs(1:2)) = [1-b1 -b1];
S(end+1,:) = ss;
% % ss = zeros(1,length(model_up.rxns));
% % ss(1,rxnIDs(3:4)) = [1-b1 -b1];
% % S(end+1,:) = ss;
ss = zeros(1,length(model_up.rxns));
ss(1,rxnIDs(5:6)) = [1-b1 -b1];
S(end+1,:) = ss;
model_up.S = S;

model_up = addMetabolite(model_up, 'iodine1_c');
model_up = addMetabolite(model_up, 'iodine2_c');
% model = addMetabolite(model, 'iodine3_c');

model_up = changeRxnBounds(model_up, {'EX_glc__D_e','EX_o2_e'}, [-10 -15],'l');
model_up = changeObjective(model_up,'BIOMASS_Ec_iJO1366_core_53p95M');
solution_1 = optimizeCbModel(model_up);
printFluxVector(model_up, solution_1.x, true, true,-1,[],[],false)
fprintf('\n')
%% production strain

model_down = model;
model_down = addReaction(model_down, 'G3P_C', 'metaboliteList', {'g3p_c',...
    'g3p_e'},'stoichCoeffList', [-1;1],'lowerBound',...
    -1000,'upperBound',1000);
model_down = addReaction(model_down,'EX_g3p_e', {'g3p_e'}, -1, true, -1000,0);
model_down = addReaction(model_down,'pksg', 'aacoa_c + accoa_c + h2o_c -> h_c + coa_c + hmgcoa_c');
model_down = addReaction(model_down,'hmg1', 'hmgcoa_c  + 2 nadph_c  + 2 h_c  -> 2 nadp_c  + coa_c  + mva_c');
model_down = addReaction(model_down,'mvak1', 'mva_c  + atp_c -> adp_c + h_c mvap_c ');
model_down = addReaction(model_down,'mvak2', 'mvap_c  + atp_c -> adp_c  + mvapp_c ');
model_down = addReaction(model_down,'mvad', 'mvapp_c  + atp_c -> adp_c  + co2_c  + pi_c + ipdp_c ');
model_down = addReaction(model_down,'ggpps', 'frdp_c  + ipdp_c  -> ggpp_c  + ppi_c');
model_down = addReaction(model_down,'txs', 'ggpp_c  -> txdn_c  + ppi_c ');
model_down = addReaction(model_down,'EX_txdn_c','txdn_c  -> ');

model_down = deleteModelGenes(model_down,...
    {'b3124', 'b3115', 'b2296', 'b1849', 'b0514', 'b2889'});



S = model_down.S;
rxnIDs = findRxnIDs(model_down,{'pksg','ACCOAC','pksg','CS',...
    'pksg','IPPS','ggpps','IPDDI','DMATT','GRTT',...
    'OCTDPS','UDCPDPS','GAPD','DXPS','GAPD','TALA'});
b1 = 0.99;
b2 = 0.99;
b3 = 0.99;
b4 = 0.99;
ss = zeros(1,length(model_down.rxns));
ss(1,rxnIDs(1:2)) = [1-b1 -b1];
S(end+1,:) = ss;
ss = zeros(1,length(model_down.rxns));
ss(1,rxnIDs(3:4)) = [1-b2 -b2];
S(end+1,:) = ss;
ss = zeros(1,length(model_down.rxns));
ss(1,rxnIDs(5:6)) = [1-b3 -b3];
S(end+1,:) = ss;

% ss = zeros(1,length(model_down.rxns));
% ss(1,rxnIDs(7:12)) = [1-b4 -b4 -b4 -b4 -b4 -b4];
% S(end+1,:) = ss;

model_down.S = S;

model_down = addMetabolite(model_down, 'iodine1_c');
model_down = addMetabolite(model_down, 'iodine2_c');
model_down = addMetabolite(model_down, 'iodine3_c');


% model_down = removeRxns(model_down,{'DDPGALA','EDA','F6PA',...
%     'TGBPA','TRPS3'});


model_down = changeRxnBounds(model_down,{'EX_glc__D_e','EX_o2_e','EX_g3p_e'...
    ,'EX_xyl__D_e'},[-10 -0.8*15 -10 -10],{'l','l','l','l'});
model_down = changeObjective(model_down,'BIOMASS_Ec_iJO1366_core_53p95M');
% save g3p_txdn_iJO1366.mat model_down
solution_1 = optimizeCbModel(model_down);
printFluxVector(model_down, solution_1.x, true, true,-1,[],[],false)
fprintf('\n')

%% solve dynamic model
t0 = 0;
tEnd = 10;
%---------------- g/l mmol/l mmol/l
time_interval = [0 50];
% x0_up = 0.05:0.05:1;
% for i=1:length(x0_up)
initial_values = [0.05;0.05;15;15;0;0];
[t,y] = ode15s(@dynamic_ecoli,time_interval,initial_values);
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
plot(t,y(:,6),'g')
legend on
ylabel('g/l')
xlabel('h')
legend('glucose','xylose','txdn','intermediate')
titer = y(end,4);
volume = 1;
tss = 6;
tend = 20.8391868387666;
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



