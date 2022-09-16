clear
clc
global model Mw
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
    },[-10,0,-15]...
    ,{'l','l','l'});

% model = deleteModelGenes(model,{'b2600'});
% model = changeObjective(model,'BIOMASS_Ec_iJO1366_core_53p95M');
% solution1 = optimizeCbModel(model);
% printFluxVector(model, solution1.x, true, true,-1,[],[],true)
% fprintf('\n')%% engineering strain


Mw = [0.18 0.06 0.3];
% time_interval = [0 12];t
%---------------- g/l mmol/l mmol/l
time_interval = [0 250];

initial_values = [0.1;30;0];
[t,y] = ode15s(@dynamic_ecoli,time_interval,initial_values);
figure
plot(t,y(:,2),'r')
hold on
plot(t,y(:,3),'b')
legend on
ylabel('g/l')
xlabel('h')
legend('glucose','salidroside')
figure
plot(t,y(:,1))
title('biomass vs time')
titer = y(end,3);
% volume = 1;
% tss = 6;
% tend = 21.2188609124807;
% Productivity = volume.*titer./(tend+tss);
% Yield = volume.*titer./volume./20;
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
% T = table(t,Biomass,Glucose,Succinate,Oxygen);
% writetable(T,'Taxadien_Mono_culture_results.xls') 


