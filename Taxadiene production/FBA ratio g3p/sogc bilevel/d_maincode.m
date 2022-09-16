%% In the name of GOD
% -------------------------------------------------------------------------
clc
clear
close all

format short g
warning off all

%% Solver checking

solverOk = changeCobraSolver('gurobi','all');

global  model

%% Load model

model = addReaction(model, 'G3P_C', 'metaboliteList', {'g3p_c',...
    'g3p_e'},'stoichCoeffList', [-1;1],'lowerBound',...
    0,'upperBound',1000);
model = addReaction(model,'EX_g3p_e', {'g3p_e'}, -1, true, 0,1000);

%% Optimization



% VTR		"Value To Reach" (stop when ofunc < VTR)
VTR = 1.e-8;

% D		number of parameters of the objective function
D = 2;

% XVmin,XVmax   vector of lower and bounds of initial population
%    		the algorithm seems to work well only if [XVmin,XVmax]
%    		covers the region where the global minimum is expected
%               *** note: these are no bound constraints!! ***
XVmin = zeros(1,2);
XVmax = ones(1,2);

% y		problem data vector (remains fixed during optimization)
y = [];

% NP            number of population members
NP = 100*D; % 15

% itermax       maximum number of iterations (generations)
itermax = 60;
% F             DE-stepsize F ex [0, 2]
F = 0.8;


% CR            crossover probabililty constant ex [0, 1]
CR = 0.6;

% strategy       1 --> DE/best/1/exp           6 --> DE/best/1/bin
%                2 --> DE/rand/1/exp           7 --> DE/rand/1/bin
%                3 --> DE/rand-to-best/1/exp   8 --> DE/rand-to-best/1/bin
%                4 --> DE/best/2/exp           9 --> DE/best/2/bin
%                5 --> DE/rand/2/exp           else  DE/rand/2/bin

strategy = 1;

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
