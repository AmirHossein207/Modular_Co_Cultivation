clear

clc
global model_1 Mw
SolverOk = changeCobraSolver('gurobi','lp');
model = readCbModel('iJO1366.mat');
model_1 = model;

%% engineering strain
model_1 = deleteModelGenes(model_1,{'b2029', 'b4388', 'b2297', 'b2458', 'b0721'});

b1 = 0.99;
b2 = 0.99;
b3 = 0.99;
rxnIDs = findRxnIDs(model_1,{'ICL','ICDHyr',...
    'SUCCt3pp','SUCDi',...
    'FRD2','FUM'});

S = model_1.S;

ss = zeros(1,length(model_1.rxns));
ss(1,rxnIDs(1:2)) = [1-b1 -b1];
S(end+1,:) = ss;
ss = zeros(1,length(model_1.rxns));
ss(1,rxnIDs(3:4)) = [1-b2 -b2];
S(end+1,:) = ss;
ss = zeros(1,length(model_1.rxns));
ss(1,rxnIDs(5:6)) = [1-b3 -b3];
S(end+1,:) = ss;



model_1.S = S;

model_1 = addMetabolite(model_1, 'iodine1_c');
model_1 = addMetabolite(model_1, 'iodine2_c');
model_1 = addMetabolite(model_1, 'iodine3_c');


Mw = [0.18 0.06];

%---------------- g/l mmol/l mmol/l
time_interval = [0 30];

initial_values = [0.1;30;0];
[t,y] = ode15s(@dynamic_ecoli,time_interval,initial_values);
figure

plot(t,y(:,2),'r')
hold on
plot(t,y(:,3),'b')

legend on
ylabel('g/l')
xlabel('h')
legend('glucose','succinate')
figure
plot(t,y(:,1))
title('biomass vs time')
titer = y(end,3);
xlabel('h')
volume = 1;
tss = 6;
tend = 23.5232942581512;
Productivity = volume.*titer./(tend+tss);
Yield = volume.*titer./volume./30;
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


