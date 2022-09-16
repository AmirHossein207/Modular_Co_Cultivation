%% In the name of GOD
% -------------------------------------------------------------------------
clc
clear
close all

format short g
warning off all

%% Solver checking

solverOk = changeCobraSolver('gurobi','all');

global  model taxadiene_max bio_max rxnID

%% Load model
load('iJO1366.mat');
model = iJO1366;

model = changeRxnBounds(model, {'EX_glc__D_e','EX_o2_e'},[-10, -15],'l');
model = addReaction(model,'HMGCOAS', 'metaboliteList', {'aacoa_c','accoa_c',...
    'coa_c','h_c','h2o_c','hmgcoa_c'}, 'stoichCoeffList', [1;1;-1;-1;1;-1],...
    'lowerBound',-1000,'upperBound',1000);
model = addReaction(model,'HMGCOAR', 'metaboliteList', {'coa_c','h_c',...
    'nadp_c','nadph_c','hmgcoa_c','mev__R_c'}, 'stoichCoeffList', [-1;2;-2;2;1;-1],...
    'lowerBound',-1000,'upperBound',1000);
model = addReaction(model,'MEVK1', 'metaboliteList', {'adp_c','atp_c',...
    'h_c','5pmev_c','mev__R_c'}, 'stoichCoeffList', [1;-1;1;1;-1],...
    'lowerBound',-1000,'upperBound',1000);

model = addReaction(model,'PMEVK', 'metaboliteList', {'adp_c','atp_c',...
    '5dpmev_c','5pmev_c'}, 'stoichCoeffList', [1;-1;1;-1],...
    'lowerBound',0,'upperBound',1000);
model = addReaction(model,'DPMVD', 'metaboliteList', {'adp_c','atp_c',...
    'co2_c','ipdp_c','pi_c','5dpmev_c'}, 'stoichCoeffList', [1;-1;1;1;1;-1],...
    'lowerBound',0,'upperBound',1000);
model = addReaction(model,'FRTT', 'metaboliteList', {'frdp_c','ipdp_c',...
    'ppi_c','ggdp_c'}, 'stoichCoeffList', [-1;-1;1;1],...
    'lowerBound',0,'upperBound',1000);
model = addReaction(model,'TXDN', 'metaboliteList', {'ggdp_c','txdn_c',...
    'ppi_c'}, 'stoichCoeffList', [1;-1;-1],...
    'lowerBound',-1000,'upperBound',1000);
lb = -100;
ub = 1000;
[model, EX_txdn_c] = addExchangeRxn(model, 'txdn_c', lb, ub);
%% Optimization
model2 = model;
model3 = model;
%
model2 = changeObjective(model2,{'BIOMASS_Ec_iJO1366_core_53p95M'});
r2 = optimizeCbModel(model2); bio_max = r2.f;
model3 = changeObjective(model3,{'EX_txdn_c'});
r3 = optimizeCbModel(model3); taxadiene_max = r3.f;

rxnID = findRxnIDs(model, {'BIOMASS_Ec_iJO1366_core_53p95M','EX_txdn_c'});

% taxadiene_max_new = r1.x(rxnID(2)); biomass_max_new = r1.x(rxnID(1));

% VTR		"Value To Reach" (stop when ofunc < VTR)
VTR = 1.e-8;

% D		number of parameters of the objective function
D = 1;

% XVmin,XVmax   vector of lower and bounds of initial population
%    		the algorithm seems to work well only if [XVmin,XVmax]
%    		covers the region where the global minimum is expected
%               *** note: these are no bound constraints!! ***
XVmin = 0; % zeros(1,2);
XVmax = 5; % ones(1,2);

% y		problem data vector (remains fixed during optimization)
y = [];

% NP            number of population members
NP = 40*D; % 15

% itermax       maximum number of iterations (generations)
itermax = 100;
% F             DE-stepsize F ex [0, 2]
F = 0.8;


% CR            crossover probabililty constant ex [0, 1]
CR = 0.6;

% strategy       1 --> DE/best/1/exp           6 --> DE/best/1/bin
%                2 --> DE/rand/1/exp           7 --> DE/rand/1/bin
%                3 --> DE/rand-to-best/1/exp   8 --> DE/rand-to-best/1/bin
%                4 --> DE/best/2/exp           9 --> DE/best/2/bin
%                5 --> DE/rand/2/exp           else  DE/rand/2/bin

strategy = 6;

% refresh       intermediate output will be produced after "refresh"
%               iterations. No intermediate output will be produced
%               if refresh is < 1
refresh = 1;

resume = false;

[x,f,nf,iter] = f_devec5('e_objfun',VTR,D,XVmin,XVmax,y,NP,itermax,F,CR,...
    strategy,refresh,resume);

disp(x')
disp(f)


a = 1; b = x;

save('a','a'); save('b','b');

model_1 = changeObjective(model, {'EX_txdn_c',...
    'BIOMASS_Ec_iJO1366_core_53p95M'},[x/taxadiene_max,1/bio_max]);
model_1 = removeRxns(model_1,{'FADRx','FLVRx','FE3Ri','PPC',...
    'THD2pp','MGSA'});
FBA1 = optimizeCbModel(model_1, 'max');
printFluxVector(model_1, FBA1.x, true, true,-1,[],[],true)
fprintf('\n')