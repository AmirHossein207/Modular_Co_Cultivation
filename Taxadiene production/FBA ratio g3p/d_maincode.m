%% In the name of GOD
% -------------------------------------------------------------------------
clc
clear
close all

format short g


%% Solver checking

solverOk = changeCobraSolver('gurobi','all');

global  model pep_max bio_max rxnID

%% Load model
load('iJO1366.mat');
model = iJO1366;
model = addReaction(model, 'G3P_C', 'metaboliteList', {'g3p_c',...
    'g3p_e'},'stoichCoeffList', [-1;1],'lowerBound',...
    0,'upperBound',1000);
model = addReaction(model,'EX_g3p_e', {'g3p_e'}, -1, true, 0,1000);

%% Optimization
rxnID = findRxnIDs(model, {'BIOMASS_Ec_iJO1366_core_53p95M','EX_g3p_e'});


model2 = model;
model3 = model;

model2 = changeRxnBounds(model2,{'EX_glc__D_e','EX_o2_e'}...
    ,[-10 -15],{'l','l'});
model2 = changeObjective(model2,{'BIOMASS_Ec_iJO1366_core_53p95M'});
r2 = optimizeCbModel(model2); bio_max = r2.f;
model3 = changeRxnBounds(model3,{'EX_glc__D_e','EX_o2_e'}...
    ,[-10 -15],{'l','l'});
model3 = changeObjective(model3,{'EX_g3p_e'});
r3 = optimizeCbModel(model3); pep_max = r3.f;


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
NP = 60*D; % 15

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

[bestmem , bestval, nfeval , iter] = o_devec6 ( 'e_objfun' , VTR , D ,...
    XVmin , XVmax, y , NP , itermax , F , CR , strategy , refresh ,...
    resume);
% ---------------


a = 1; b = bestmem;

save('a','a'); save('b','b');
model_1 = changeRxnBounds(model,{'EX_glc__D_e','EX_o2_e'}...
    ,[-10 -15],{'l','l'});
model_1 = changeObjective(model_1,{'BIOMASS_Ec_iJO1366_core_53p95M','EX_g3p_e'}...
    ,[(a/bio_max) (b/pep_max)]);
FBA1 = optimizeCbModel(model_1);
printFluxVector(model_1, FBA1.x, true, true,-1,[],[],true)
fprintf('\n')