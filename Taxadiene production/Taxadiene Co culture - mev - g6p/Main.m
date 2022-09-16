clear
clc

global model_up model_down Mw

SolverOk = changeCobraSolver('gurobi','lp');
load('iJO1366.mat');
model = iJO1366;

Mw = [0.18 0.06 0.272 0.148 0.26 0.15] ;


%% intermediate strain
model_up = model;
model_up = removeRxns(model_up,'G6Pt6_2pp');
model_up = addReaction(model_up, 'G6Pt6_2pp_reverse', 'metaboliteList', {'pi_c',...
    'g6p_p','g6p_c','pi_p'},'stoichCoeffList', [2;1;-1;-2],'lowerBound',...
    -1000,'upperBound',1000);

model_up = removeRxns(model_up,'PGMT');
model_up = addReaction(model_up, 'PGMT', 'metaboliteList', {'g1p_c','g6p_c'}...
    ,'stoichCoeffList', [1;-1],'lowerBound',...
    -1000,'upperBound',1000);

% model_up = changeRxnBounds(model_up, {'EX_g6p_e','EX_g6p_e'...
%     },[-1000,0],{'l','u'});

% b1 = 0.91;
% b2 = 0.90;
% b3 = 0.91;
% b4 = 0.92;
% b5 = 0.94;
% b6 = 0.960;
b1 = 0.99;
b2 = 0.99;
b3 = 0.99;
b4 = 0.99;
b5 = 0.99;
model_up = deleteModelGenes(model_up,...
    {'b3115', 'b2296', 'b1849', 'b4015', 'b2913', 'b3519'});

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

% ss = zeros(1,length(model_up.rxns));
% ss(1,rxnIDs(17:18)) = [1-b6 -b6];
% S(end+1,:) = ss;

model_up.S = S;

model_up = addMetabolite(model_up, 'iodine1_c');
model_up = addMetabolite(model_up, 'iodine2_c');
model_up = addMetabolite(model_up, 'iodine3_c');
model_up = addMetabolite(model_up, 'iodine4_c');
model_up = addMetabolite(model_up, 'iodine5_c');
% model_up = addMetabolite(model_up, 'iodine6_c');
% model = addReaction(model,'EX_mva_e','mva_e ⇌');
model_up = addReaction(model_up,'EX_mva_e', {'mva_e'}, -1, true, 0,1000);
% model = addReaction(model,'MEVRt','mva_c ⇌ mva_e');
model_up = addReaction(model_up,'MEVRt', {'mva_c','mva_e'}, [1 -1], true,-1000, 1000);
% model = addReaction(model,'MVNOR', 'coa_c + nad_c + mva_c → 4.0 h_c + nadh_c + hmgcoa_c');
% model_up = addReaction(model_up,'MVNOR', {'coa_c','nad_c','mva_c','h_c','nadh_c','hmgcoa_c'},[-1;-1;-1;4;1;1],false,0,1000);
% model = addReaction(model,'pksg', 'aacoa_c + accoa_c + h2o_c -> h_c + coa_c + hmgcoa_c');
model_up = addReaction(model_up,'pksg', {'aacoa_c','accoa_c','h2o_c','h_c','coa_c','hmgcoa_c'},[-1;-1;-1;1;1;1],false,0,1000);
% model = addReaction(model,'hmg1', 'hmgcoa_c  + 2 nadph_c  + 2 h_c  -> 2 nadp_c  + coa_c  + mva_c');
model_up = addReaction(model_up,'hmg1', {'hmgcoa_c','nadph_c','h_c','nadp_c', 'coa_c', 'mva_c'},[-1;-2;-2;2;1;1],false,0,1000);
% model = addReaction(model,'mvak1', 'mva_c  + atp_c -> adp_c + h_c mvap_c ');
% model_up = addReaction(model_up,'mvak1',{'mva_c','atp_c','adp_c','h_c','mvap_c'},[-1;-1;1;1;1],false,0,1000);
% model = addReaction(model,'mvak2', 'mvap_c  + atp_c -> adp_c  + mvapp_c ');
% model_up = addReaction(model_up,'mvak2',{'mvap_c','atp_c','adp_c','mvapp_c'},[-1;-1;1;1],false,0,1000);
% model = addReaction(model,'mvad', 'mvapp_c  + atp_c -> adp_c  + co2_c  + pi_c + ipdp_c');
% model_up = addReaction(model_up,'mvad', {'mvapp_c','atp_c','adp_c','co2_c','pi_c', 'ipdp_c'},[-1;-1;1;1;1;1],false,0,1000);
% model = addReaction(model,'ggpps', 'frdp_c  + ipdp_c  -> ggpp_c  + ppi_c');
% model_up = addReaction(model_up,'ggpps', {'frdp_c','ipdp_c','ggpp_c','ppi_c'},[-1;-1;1;1],false,0,1000);
% % model = addReaction(model,'txs', 'ggpp_c  -> txdn_c  + ppi_c ');
% model = addReaction(model,'txs', {'ggpp_c','txdn_c','ppi_c'},[-1;1;1],false,0,1000);
% model = addReaction(model,'txdn_ex', {'txdn_c','txdn_e'}, [-1 1], true, -1000, 1000);
% model = addReaction(model,'Ex_txdn_e', {'txdn_e'}, -1, false, 0, 1000);

model_up = deleteModelGenes(model_up,...
    {'b0351', 'b1241', 'b4069', 'b0728'});

S = model_up.S;
rxnIDs = findRxnIDs(model_up,{'pksg','ACCOAC','pksg','CS',...
    'pksg','IPPS','ggpps','IPDDI','DMATT','GRTT',...
    'OCTDPS','UDCPDPS','GAPD','DXPS','GAPD','TALA'});
b1 = 0.99;
b2 = 0.99;
b3 = 0.99;

ss = zeros(1,length(model_up.rxns));
ss(1,rxnIDs(1:2)) = [1-b1 -b1]; 
S(end+1,:) = ss;
ss = zeros(1,length(model_up.rxns));
ss(1,rxnIDs(3:4)) = [1-b2 -b2]; 
S(end+1,:) = ss;
ss = zeros(1,length(model_up.rxns));
ss(1,rxnIDs(5:6)) = [1-b3 -b3]; 
S(end+1,:) = ss;


model_up.S = S;

model_up = addMetabolite(model_up, 'iodine6_c');
model_up = addMetabolite(model_up, 'iodine7_c');
model_up = addMetabolite(model_up, 'iodine8_c');

model_up = changeRxnBounds(model_up, {'EX_glc__D_e','EX_o2_e'}, [-10 -15],'l');
model_up = changeObjective(model_up,'BIOMASS_Ec_iJO1366_core_53p95M');
solution_1 = optimizeCbModel(model_up);
printFluxVector(model_up, solution_1.x, true, true,-1,[],[],false)
fprintf('\n')
%% production strain

model_down = model;

% model_down = removeRxns(model_down,'G6Pt6_2pp');
% model_down = addReaction(model_down, 'G6Pt6_2pp_reverse', 'metaboliteList', {'pi_c',...
%     'g6p_p','g6p_c','pi_p'},'stoichCoeffList', [2;1;-1;-2],'lowerBound',...
%     1000,'upperBound',1000);
% 
% model_down = removeRxns(model_down,'PGMT');
% model_down = addReaction(model_down, 'PGMT', 'metaboliteList', {'g1p_c','g6p_c'}...
%     ,'stoichCoeffList', [1;-1],'lowerBound',...
%     -1000,'upperBound',1000);


model_down = changeRxnBounds(model_down, {'EX_g6p_e','EX_g6p_e'...
    },[-1000,0],{'l','u'});

model_down = addReaction(model_down,'EX_mva_e', {'mva_e'}, -1, false, -1000,1000);
% % % % model = addReaction(model,'MEVRt','mva_c ⇌ mva_e');
model_down = addReaction(model_down,'MEVRt', {'mva_c','mva_e'}, [1 -1], true,0, 1000);

% model = addReaction(model,'pksg', 'aacoa_c + accoa_c + h2o_c -> h_c + coa_c + hmgcoa_c');
model_down = addReaction(model_down,'pksg', {'aacoa_c','accoa_c','h2o_c','h_c','coa_c','hmgcoa_c'},[-1;-1;-1;1;1;1],false,0,1000);
% % model = addReaction(model,'hmg1', 'hmgcoa_c  + 2 nadph_c  + 2 h_c  -> 2 nadp_c  + coa_c  + mva_c');
model_down = addReaction(model_down,'hmg1', {'hmgcoa_c','nadph_c','h_c','nadp_c', 'coa_c', 'mva_c'},[-1;-2;-2;2;1;1],false,0,1000);
% model = addReaction(model,'mvak1', 'mva_c  + atp_c -> adp_c + h_c mvap_c ');
model_down = addReaction(model_down,'mvak1',{'mva_c','atp_c','adp_c','h_c','mvap_c'},[-1;-1;1;1;1],false,0,1000);
% model = addReaction(model,'mvak2', 'mvap_c  + atp_c -> adp_c  + mvapp_c ');
model_down = addReaction(model_down,'mvak2',{'mvap_c','atp_c','adp_c','mvapp_c'},[-1;-1;1;1],false,0,1000);
% model = addReaction(model,'mvad', 'mvapp_c  + atp_c -> adp_c  + co2_c  + pi_c + ipdp_c');
model_down = addReaction(model_down,'mvad', {'mvapp_c','atp_c','adp_c','co2_c','pi_c', 'ipdp_c'},[-1;-1;1;1;1;1],false,0,1000);
% model = addReaction(model,'ggpps', 'frdp_c  + ipdp_c  -> ggpp_c  + ppi_c');
model_down = addReaction(model_down,'ggpps', {'frdp_c','ipdp_c','ggpp_c','ppi_c'},[-1;-1;1;1],false,0,1000);
% model = addReaction(model,'txs', 'ggpp_c  -> txdn_c  + ppi_c ');
model_down = addReaction(model_down,'txs', {'ggpp_c','txdn_c','ppi_c'},[-1;1;1],false,0,1000);
model_down = addReaction(model_down,'txdn_ex', {'txdn_c','txdn_e'}, [-1 1], true, 0, 1000);
model_down = addReaction(model_down,'Ex_txdn_e', {'txdn_e'}, -1, false, 0, 1000);

model_down = deleteModelGenes(model_down,...
    {'b3115', 'b2296', 'b1849', 'b3956', 'b2464', 'b0008', 'b0505'});
% 'b3744', 'b3956', 'b2407', 'b4384', 'b0505'
S = model_down.S;

rxnIDs = findRxnIDs(model_down,{'pksg','ACCOAC','pksg','CS',...
    'pksg','IPPS','ggpps','IPDDI','DMATT','GRTT',...
    'OCTDPS','UDCPDPS','GAPD','DXPS','GAPD','TALA','mvak1','mev_ex'});

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
model_down = addMetabolite(model_down, 'iodine4_c');
model_down = changeRxnBounds(model_down,{'EX_glc__D_e','EX_o2_e',...
    'EX_mva_e','EX_g6p_e','EX_xyl__D_e'},...
    [0 -0.8*15 -5 -5 -10],{'l','l','b','l','l'});
model_down = changeObjective(model_down,'BIOMASS_Ec_iJO1366_core_53p95M');
solution_1 = optimizeCbModel(model_down);
printFluxVector(model_down, solution_1.x, true, true,-1,[],[],false)
fprintf('\n')
% 
%% solve dynamic model
t0 = 0;
tEnd = 10;
%---------------- g/l mmol/l mmol/l
time_interval = [0 90];
% x0_up = 0.05:0.05:1;
% for i=1:length(x0_up)
initial_values = [0.05;0.05;15;15;0;0;0];
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
plot(t,y(:,6),'k')
hold on
plot(t,y(:,7),'r')
legend on
ylabel('g/l')
xlabel('h')
legend('glucose','xylose','txdn','mev','g6p')
titer = y(end,5)
volume = 1;
tss = 6;
tend = 84.65
Productivity = volume.*titer./(tend+tss)
Yield = volume.*titer./volume./20
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



