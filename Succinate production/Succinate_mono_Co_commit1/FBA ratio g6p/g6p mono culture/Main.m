clear

clc
global model_1 Mw
SolverOk = changeCobraSolver('gurobi','lp');
model = readCbModel('iJO1366.mat');

model = removeRxns(model,'G6Pt6_2pp');
model = addReaction(model, 'G6Pt6_2pp_reverse', 'metaboliteList', {'pi_c',...
    'g6p_p','g6p_c','pi_p'},'stoichCoeffList', [2;1;-1;-2],'lowerBound',...
    0,'upperBound',1000);

model = removeRxns(model,'PGMT');
model = addReaction(model, 'PGMT', 'metaboliteList', {'g1p_c','g6p_c'}...
    ,'stoichCoeffList', [1;-1],'lowerBound',...
    -1000,'upperBound',1000);


model = changeRxnBounds(model, {'EX_glc__D_e','EX_o2_e',...
    'EX_for_e'},[-10,-15,0],{'l','l','b'});

b1 = 0.91;
b2 = 0.90;
b3 = 0.91;
b4 = 0.92;
b5 = 0.94;
b6 = 0.960;
rxnIDs = findRxnIDs(model,{'G6Pt6_2pp_reverse','G6PDH2r',...
    'G6Pt6_2pp_reverse','PGI',...
    'G6Pt6_2pp_reverse','PGMT',...
    'G6Pt6_2pp_reverse','G6PP',...
    'GLCptspp','DHAPT','UAGCVT','PSCVT','PPC','KDOPS','DDPA','PYK',...
    'GLCptspp','GLCt2pp'});

S = model.S;



ss = zeros(1,length(model.rxns));
ss(1,rxnIDs(1:2)) = [1-b1 -b1];
S(end+1,:) = ss;

ss = zeros(1,length(model.rxns));
ss(1,rxnIDs(3:4)) = [1-b2 -b2];
S(end+1,:) = ss;

ss = zeros(1,length(model.rxns));
ss(1,rxnIDs(5:6)) = [1-b3 -b3];
S(end+1,:) = ss;

ss = zeros(1,length(model.rxns));
ss(1,rxnIDs(7:8)) = [1-b4 -b4];
S(end+1,:) = ss;

ss = zeros(1,length(model.rxns));
ss(1,rxnIDs(9:16)) = [1-b5 -b5 -b5 -b5 -b5 -b5 -b5 -b5];
S(end+1,:) = ss;

% ss = zeros(1,length(model.rxns));
% ss(1,rxnIDs(17:18)) = [1-b6 -b6];
% S(end+1,:) = ss;

model.S = S;

model = addMetabolite(model, 'iodine1_c');
model = addMetabolite(model, 'iodine2_c');
model = addMetabolite(model, 'iodine3_c');
model = addMetabolite(model, 'iodine4_c');
model = addMetabolite(model, 'iodine5_c');

model_1 = model;

Mw = [0.18 0.06 0.260];

%---------------- g/l mmol/l mmol/l
time_interval = [0 20];

initial_values = [0.05;15;0;40];
[t,y] = ode23s(@dynamic_ecoli,time_interval,initial_values);
figure
plot(t,y(:,1),'g')
hold on
plot(t,y(:,2),'r')
hold on
plot(t,y(:,3),'b')
hold on
plot(t,y(:,4),'c')
legend on
ylabel('g/l')
xlabel('h')
legend('biomass','glucose','succ','oxygen')
figure
plot(t,y(:,1))
title('biomass vs time')

%% calculating miu
% [vg,vo] = menten(y(:,2),y(:,4),y(:,3));
% for i=1:length(vg)
%     model_1 = changeRxnBounds(model_1,{'EX_glc__D_e','EX_o2_e'}...
%         ,[-vg(i) -vo(i)],{'l','l'});
%     solution = optimizeCbModel(model_1);
%     miu(i) = solution.f;
% end
% figure
% plot(t,miu)
% title('growth rate vs time')
% Biomass = y(:,1);
% Glucose = y(:,2);
% Succinate = y(:,3);
% Oxygen = y(:,4);
% % T = table(t,Biomass,Glucose,Succinate,Oxygen);
% writetable(T,'Succinate_Mono_culture_results.xls') 


