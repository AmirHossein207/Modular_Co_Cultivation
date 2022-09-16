clc
close 
clear
SolverOk = changeCobraSolver('gurobi','lp');
model = readCbModel('iJO1366.mat');

%% add exchange reaction for intermediate
model = addReaction(model,'TYROSpp', 'metaboliteList', {'4_tyrosol_c','4_tyrosol_e'},...
    'stoichCoeffList', [-1;1], 'lowerBound',-1000,'upperBound',0);
[model, EX_4_tyrosol_e] = addExchangeRxn(model, '4_tyrosol_e',-1000,0);
model = addReaction(model,'SALID', 'metaboliteList', {'udpg_c','4_tyrosol_c',...
    'udp_c','h_c','salid_c'}, 'stoichCoeffList', [-1;-1;1;1;1],...
    'lowerBound',0,'upperBound',1000);
model = addReaction(model,'SALIDpp', 'metaboliteList', {'salid_c','salid_e'},...
    'stoichCoeffList', [-1;1;], 'lowerBound',0,'upperBound',1000);
[model, EX_salid_e] = addExchangeRxn(model, 'salid_e',0,1000);

model = removeRxns(model,'PGMT');
model = addReaction(model, 'PGMT', 'metaboliteList', {'g1p_c','g6p_c'}...
   ,'stoichCoeffList', [1;-1],'lowerBound',-1000,'upperBound',1000);
model = removeRxns(model,'MLTP2');
model = addReaction(model, 'MLTP2', 'metaboliteList', {'malthx_c','pi_c',...
    'g1p_c','maltpt_c'},'stoichCoeffList', [1;1;-1;-1],'lowerBound',...
    -1000,'upperBound',1000);
model = removeRxns(model,'MLTP1');
model = addReaction(model, 'MLTP1', 'metaboliteList', {'maltpt_c','pi_c',...
    'g1p_c','maltttr_c'},'stoichCoeffList', [1;1;-1;-1],'lowerBound',...
    -1000,'upperBound',1000);
model = removeRxns(model,'MLTP3');
model = addReaction(model, 'MLTP3', 'metaboliteList', {'malthp_c','pi_c',...
    'g1p_c','malthx_c'},'stoichCoeffList', [1;1;-1;-1],'lowerBound',...
    -1000,'upperBound',1000);
model = removeRxns(model,'PHEt2rpp');
model = addReaction(model, 'PHEt2rpp', 'metaboliteList', {'h_p','phe__L_p',...
    'h_c','phe__L_c'},'stoichCoeffList', [1;1;-1;-1],'lowerBound',...
    -1000,'upperBound',1000);
% model = removeRxns(model,'ASPTA');
% model = addReaction(model, 'ASPTA', 'metaboliteList', {'akg_c','asp__L_c',...
%     'glu__L_c','oaa_c'},'stoichCoeffList', [1;1;-1;-1],'lowerBound',...
%     -1000,'upperBound',1000);
% 
rxnIDs = findRxnIDs(model,{'PGMT','PGI','PGMT','G6PDH2r','GALUi',...
    'MLTP3','MLTP1','MLTP2','SALID','TRE6PS','HEX1','XYLI2','GLCDpp',...
    'PHEt2rpp','BIOMASS_Ec_iJO1366_core_53p95M'});
b1 = 0.54;
b2 = 0.99;
b3 = 0.01;

S = model.S;

ss = zeros(1,length(model.rxns));
ss(1,rxnIDs(1:2)) = [1-b1 -b1];
S(end+1,:) = ss;
ss = zeros(1,length(model.rxns));
ss(1,rxnIDs(3:4)) = [1-b1 -b1];
S(end+1,:) = ss;
ss = zeros(1,length(model.rxns));
ss(1,rxnIDs(5:8)) = [1-b2 -b2 -b2 -b2];
S(end+1,:) = ss;
ss = zeros(1,length(model.rxns));
ss(1,rxnIDs(9:10)) = [1-b2 -b2];
S(end+1,:) = ss;
ss = zeros(1,length(model.rxns));
ss(1,rxnIDs(11:13)) = [1-b2 -b2 -b2];
S(end+1,:) = ss;
ss = zeros(1,length(model.rxns));
ss(1,rxnIDs(14:15)) = [1-b3 -b3];
S(end+1,:) = ss;
model.S = S;

model = addMetabolite(model, 'iodine1_c');
model = addMetabolite(model, 'iodine2_c');
model = addMetabolite(model, 'iodine3_c');
model = addMetabolite(model, 'iodine4_c');
model = addMetabolite(model, 'iodine5_c');
model = addMetabolite(model, 'iodine6_c');
% model = addMetabolite(model, 'iodine7_c');


%% find candidate minimal flux for production
% model = removeRxns(model,{'GLUDy','MOX'});
% model = deleteModelGenes(model,{'b2599'});
[non_zero_flux_D,FBA_solution_D,non_zero_flux_W,FBA_solution_W] =...
    pathway_min_norm(model);
printFluxVector(model, FBA_solution_W, true, true,-1,[],[],true)

%% visualization of flux
% [Involved_mets_D, Dead_ends_D] = draw_by_rxn(model, non_zero_flux_D, 'true', ...
%     'struc', {'tyr__L_c'}, {'pi_c','h_c','nadph_c','nadp_c','nad_c','nadh_c','h_e','h2o_c','h2o_e','h2o_p','h_p',...
%     'atp_c','adp_c','co2_c','co2_e','co2_p'}, FBA_solution_D);
% [Involved_mets_W, Dead_ends_W] = draw_by_rxn(model, non_zero_flux_W, 'true', ...
%     'struc', {'pep_p'}, {'h_c','nad_c','nadh_c','h_e','h2o_c','h2o_e','h2o_p','h_p',...
%     'atp_c','adp_c','co2_c','co2_e','co2_p'}, FBA_solution_W);

%% creat the graph of the flux
% [flux_model_D,mets_index_D,rxns_index_D] = creat_flux_graph(model,Involved_mets_D,non_zero_flux_D);
% 
% [flux_model_W,mets_index_W,rxns_index_W] = creat_flux_graph(model,Involved_mets_W,non_zero_flux_W);

% 
% [directionRxns_D, involvedMets_D, deadEnds_D] = draw_by_met(model, {'3pg_c'},...
%     'true', 1, 'prod', {''}, FBA_solution_D);
[directionRxns_W, involvedMets_W, deadEnds_W] = draw_by_met(model, {'phe__L_c'},...
    'true', 1, 'prod', {''}, FBA_solution_W);
