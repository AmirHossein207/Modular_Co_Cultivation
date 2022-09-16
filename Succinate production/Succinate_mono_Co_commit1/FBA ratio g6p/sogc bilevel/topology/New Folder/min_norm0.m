clear
clc

solverOK = changeCobraSolver('gurobi','lp');
load('iJO1366.mat');
model = iJO1366;
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

model = deleteModelGenes(model,...
    {'b0521', 'b0323', 'b2874', 'b2551', 'b0870', 'b1761', 'b4015'});

rxnIDs = findRxnIDs(model,{'G6Pt6_2pp_reverse','G6PDH2r',...
    'G6Pt6_2pp_reverse','PGI',...
    'G6Pt6_2pp_reverse','PGMT',...
    'G6Pt6_2pp_reverse','G6PP',...
    'GLCptspp','DHAPT','UAGCVT','PSCVT','PPC','KDOPS','DDPA','PYK',...
    'GLCptspp','GLCt2pp'});


S = model.S;
b1 = 0.99;
b2 = 0.99;
b3 = 0.99;
b4 = 0.99;
b5 = 0.99;

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

model.S = S;

model = addMetabolite(model, 'iodine1_c');
model = addMetabolite(model, 'iodine2_c');
model = addMetabolite(model, 'iodine3_c');
model = addMetabolite(model, 'iodine4_c');
model = addMetabolite(model, 'iodine5_c');
%% find candidate minimal flux for production
[non_zero_flux_D,FBA_solution_D,non_zero_flux_W,FBA_solution_W] =...
    pathway_min_norm(model);
printFluxVector(model, FBA_solution_D, true, true,-1,[],[],true)
%% visualization of flux
[Involved_mets_D, Dead_ends_D] = draw_by_rxn(model, non_zero_flux_D, 'true', ...
    'struc', {'g6p_c'}, {'nadp_c','nadph_c','pi_c','h_c','nad_c','nadh_c','h_e','h2o_c','h2o_e','h2o_p','h_p',...
    'atp_c','adp_c','co2_c','co2_e','co2_p'}, FBA_solution_D);
% [Involved_mets_W, Dead_ends_W] = draw_by_rxn(model1, non_zero_flux_W, 'true', ...
%     'struc', {'pep_p'}, {'h_c','nad_c','nadh_c','h_e','h2o_c','h2o_e','h2o_p','h_p',...
%     'atp_c','adp_c','co2_c','co2_e','co2_p'}, FBA_solution_W);

%% creat the graph of the flux
% [flux_model_D,mets_index_D,rxns_index_D] = creat_flux_graph(model1,Involved_mets_D,non_zero_flux_D);
% 
% [flux_model_W,mets_index_W,rxns_index_W] = creat_flux_graph(model1,Involved_mets_W,non_zero_flux_W);

% 
% [directionRxns_D, involvedMets_D, deadEnds_D] = draw_by_met(model1, {'glc__D_p'},...
%     'true', 1, 'prod', {''}, FBA_solution_D);
% [directionRxns_W, involvedMets_W, deadEnds_W] = draw_by_met(model, {'g6p_c'},...
%     'true', 1, 'prod', {''}, FBA_solution_W);
