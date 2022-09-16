clear 
clc

solverOK = changeCobraSolver('gurobi','lp');
load('iJO1366.mat');
model = iJO1366;
model = removeRxns(model,'G6Pt6_2pp');
model = addReaction(model, 'G6Pt6_2pp', 'metaboliteList', {'pi_c',...
    'g6p_p','g6p_c','pi_p'},'stoichCoeffList', [2;1;-1;-2],'lowerBound',...
    -1000,'upperBound',1000);

% model = removeRxns(model,'PGMT');
% model = addReaction(model, 'PGMT', 'metaboliteList', {'g1p_c','g6p_c'}...
%     ,'stoichCoeffList', [1;-1],'lowerBound',...
%     -1000,'upperBound',1000);
model = removeRxns(model,'RPE');
model = addReaction(model, 'RPE', 'metaboliteList', {'ru5p__D_c','xu5p__D_c'}...
    ,'stoichCoeffList', [1;-1],'lowerBound',...
    -1000,'upperBound',1000);


% 
% model = deleteModelGenes(model,...
%     {'b0521', 'b0323', 'b2874', 'b2551', 'b0870', 'b1761', 'b4015'});

rxnIDs = findRxnIDs(model,{'G6Pt6_2pp','G6PDH2r',...
    'G6Pt6_2pp','PGI',...
    'G6Pt6_2pp','PGMT',...
    'G6Pt6_2pp','G6PP',...
    'GLCptspp','DHAPT','UAGCVT','PSCVT','PPC','KDOPS','DDPA','PYK',...
    'GLCptspp','GLCt2pp','TKT2','TKT1','RPE'});


S = model.S;
b1 = 0.60;
b2 = 0.99;
b3 = 0.99;
b4 = 0.50;
b5 = 0.99;

ss = zeros(1,length(model.rxns));
ss(1,rxnIDs(1:2)) = [1-b1 -b1];
S(end+1,:) = ss;

% ss = zeros(1,length(model.rxns));
% ss(1,rxnIDs(3:4)) = [1-b2 -b2];
% S(end+1,:) = ss;
% % 
% ss = zeros(1,length(model.rxns));
% ss(1,rxnIDs(5:6)) = [1-b3 -b3];
% S(end+1,:) = ss;
% % 
% ss = zeros(1,length(model.rxns));
% ss(1,rxnIDs(7:8)) = [1-b4 -b4];
% S(end+1,:) = ss;

% ss = zeros(1,length(model.rxns));
% ss(1,rxnIDs(9:16)) = [1-b5 -b5 -b5 -b5 -b5 -b5 -b5 -b5];
% S(end+1,:) = ss;
ss = zeros(1,length(model.rxns));
ss(1,rxnIDs(19:21)) = [1-b4 -b4 -b4];
S(end+1,:) = ss;
model.S = S;

model = addMetabolite(model, 'iodine1_c');
model = addMetabolite(model, 'iodine2_c');
model = addMetabolite(model, 'iodine3_c');
model = addMetabolite(model, 'iodine4_c');
model = addMetabolite(model, 'iodine5_c');
model = addMetabolite(model, 'iodine6_c');
model = changeRxnBounds(model, {'EX_glc__D_e','EX_o2_e',...
    'EX_xyl__D_e'},[0,-15,-10],{'l','l','l'});
model = changeObjective(model,'BIOMASS_Ec_iJO1366_core_53p95M');
solution_1 = optimizeCbModel(model);
printFluxVector(model, solution_1.x, true, true,-1,[],[],true)
fprintf('\n')
[directionRxns_D, involvedMets_D, deadEnds_D] = draw_by_met(model, {'g6p_c'},...
    'true', 1, 'prod', {''}, solution_1.x);
% [minFlux, maxFlux] = fluxVariability(model, 60, 'max', ...
%     {'EX_g6p_e','G6Pt6_2pp_reverse','UAGCVT','PSCVT','PPC','KDOPS','DDPA','DHAPT','PYK'})




% [w2,~] = SoGC_calculator(model,{'EX_g6p_e'},{'BIOMASS_Ec_iJO1366_core_53p95M'});









