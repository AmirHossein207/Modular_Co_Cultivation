%% In the name of GOD
% -------------------------------------------------------------------------
clc
clear
close all

format short g
warning off all

%% Solver checking

solverOk = changeCobraSolver('gurobi','all');

global  model g3p_max bio_max rxnID

%% Load model
load('iJO1366.mat');
model = iJO1366;
% model = addReaction(model,'EX_mva_e','mva_e ⇌');
model = addReaction(model,'EX_mva_e', {'mva_e'}, -1, true, 0,1000);
% model = addReaction(model,'MEVRt','mva_c ⇌ mva_e');
model = addReaction(model,'MEVRt', {'mva_c','mva_e'}, [1 -1], true,-1000, 1000);
% model = addReaction(model,'MVNOR', 'coa_c + nad_c + mva_c → 4.0 h_c + nadh_c + hmgcoa_c');
% model_up = addReaction(model_up,'MVNOR', {'coa_c','nad_c','mva_c','h_c','nadh_c','hmgcoa_c'},[-1;-1;-1;4;1;1],false,0,1000);
% model = addReaction(model,'pksg', 'aacoa_c + accoa_c + h2o_c -> h_c + coa_c + hmgcoa_c');
model = addReaction(model,'pksg', {'aacoa_c','accoa_c','h2o_c','h_c','coa_c','hmgcoa_c'},[-1;-1;-1;1;1;1],false,0,1000);
% model = addReaction(model,'hmg1', 'hmgcoa_c  + 2 nadph_c  + 2 h_c  -> 2 nadp_c  + coa_c  + mva_c');
model = addReaction(model,'hmg1', {'hmgcoa_c','nadph_c','h_c','nadp_c', 'coa_c', 'mva_c'},[-1;-2;-2;2;1;1],false,0,1000);


%% Optimization
model2 = model;
model3 = model;
%
model2 = changeObjective(model2,{'BIOMASS_Ec_iJO1366_core_53p95M'});
r2 = optimizeCbModel(model2); bio_max = r2.f;
model3 = changeObjective(model3,{'EX_mva_e'});
r3 = optimizeCbModel(model3); g3p_max = r3.f;

rxnID = findRxnIDs(model, {'BIOMASS_Ec_iJO1366_core_53p95M','EX_mva_e'});

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

model_1 = changeObjective(model, {'EX_mva_e',...
    'BIOMASS_Ec_iJO1366_core_53p95M'},[x/g3p_max,1/bio_max]);
FBA1 = optimizeCbModel(model_1, 'max');
printFluxVector(model_1, FBA1.x, true, true,-1,[],[],true)
fprintf('\n')