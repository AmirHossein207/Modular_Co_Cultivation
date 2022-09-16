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

load('iJO1366.mat');
model = iJO1366;

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
model  =  addReaction(model,'taxadiene_ex', {'txdn_c','txdn_e'}, [-1 1], true, -1000, 1000);
lb = -100;
ub = 1000;
[model, EX_txdn_e] = addExchangeRxn(model, 'txdn_e', lb, ub);

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
NP = 50*D; % 15

% itermax       maximum number of iterations (generations)
itermax = 50;
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
