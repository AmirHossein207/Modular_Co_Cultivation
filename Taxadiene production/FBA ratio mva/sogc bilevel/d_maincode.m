%% In the name of GOD
% -------------------------------------------------------------------------
clc
clear
global  model

format short g
warning off all

%% Solver checking

SolverOk = changeCobraSolver('gurobi','lp');
load('iJO1366.mat');
model = iJO1366;



%% Load model

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



% VTR		"Value To Reach" (stop when ofunc < VTR)
VTR = 1.e-8;

% D		number of parameters of the objective function
D = 3;

% XVmin,XVmax   vector of lower and bounds of initial population
%    		the algorithm seems to work well only if [XVmin,XVmax]
%    		covers the region where the global minimum is expected
%               *** note: these are no bound constraints!! ***
XVmin = zeros(1,3);
XVmax = ones(1,3);

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

% strategy       1 --> DE /best/1/exp           6 --> DE/best/1/bin
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
%
% [x,f,nf,iter] = f_devec5('e_objfun',VTR,D,XVmin,XVmax,y,NP,itermax,F,CR,...
%     strategy,refresh,resume);
[bestmem , bestval, nfeval , iter] = o_devec6 ( 'e_objfun' , VTR , D ,...
    XVmin , XVmax, y , NP , itermax , F , CR , strategy , refresh ,...
    resume);




% a = 1;
b = bestmem;

% save('a','a');
save('b','b');
