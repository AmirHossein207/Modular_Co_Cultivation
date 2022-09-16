clear
clc

solverOK = changeCobraSolver('gurobi','lp');
load('iJO1366.mat');
model = iJO1366;


model = addReaction(model,'pksg', {'aacoa_c','accoa_c','h2o_c','h_c','coa_c','hmgcoa_c'},[-1;-1;-1;1;1;1],false,0,1000);
model = addReaction(model,'hmg1', {'hmgcoa_c','nadph_c','h_c','nadp_c', 'coa_c', 'mva_c'},[-1;-2;-2;2;1;1],false,0,1000);
model = addReaction(model,'mvak1',{'mva_c','atp_c','adp_c','h_c','mvap_c'},[-1;-1;1;1;1],false,0,1000);
model = addReaction(model,'mvak2',{'mvap_c','atp_c','adp_c','mvapp_c'},[-1;-1;1;1],false,0,1000);
model = addReaction(model,'mvad', {'mvapp_c','atp_c','adp_c','co2_c','pi_c', 'ipdp_c'},[-1;-1;1;1;1;1],false,0,1000);
model = addReaction(model,'ggpps', {'frdp_c','ipdp_c','ggpp_c','ppi_c'},[-1;-1;1;1],false,0,1000);
model = addReaction(model,'txs', {'ggpp_c','txdn_c','ppi_c'},[-1;1;1],false,0,1000);
model = addReaction(model,'txdn_ex', {'txdn_c','txdn_e'}, [-1 1], true, 0, 1000);
model = addReaction(model,'Ex_txdn_e', {'txdn_e'}, -1, false, 0, 1000);


S = model.S;
rxnIDs = findRxnIDs(model,{'pksg','ACCOAC','pksg','CS',...
    'pksg','IPPS','ggpps','IPDDI','DMATT','GRTT',...
    'OCTDPS','UDCPDPS','GAPD','DXPS','GAPD','TALA','mvak1','mev_ex'});

b1 = 0.99;
b2 = 0.99;
b3 = 0.99;
b4 = 0.99;
ss = zeros(1,length(model.rxns));
ss(1,rxnIDs(1:2)) = [1-b1 -b1]; 
S(end+1,:) = ss;
ss = zeros(1,length(model.rxns));
ss(1,rxnIDs(3:4)) = [1-b2 -b2]; 
S(end+1,:) = ss;
ss = zeros(1,length(model.rxns));
ss(1,rxnIDs(5:6)) = [1-b3 -b3]; 
S(end+1,:) = ss;

model.S = S;

model = addMetabolite(model, 'iodine1_c');
model = addMetabolite(model, 'iodine2_c');
model = addMetabolite(model, 'iodine3_c');

model = changeRxnBounds(model, {'EX_glc__D_e','EX_o2_e','EX_xyl__D_e'},...
    [0, -15, -10],'l');


model = changeObjective(model,'BIOMASS_Ec_iJO1366_core_53p95M');
solution_1 = optimizeCbModel(model,'max');
printFluxVector(model, solution_1.x, true, true,-1,[],[],false)
fprintf('\n')
% [non_zero_flux_D,FBA_solution_D,non_zero_flux_W,FBA_solution_W] =...
%     pathway_min_norm(model);

%% visualization of flux
[non_zero_flux_D,FBA_solution_D,non_zero_flux_W,FBA_solution_W] =...
    pathway_min_norm(model);
% [Involved_mets_D, Dead_ends_D] = draw_by_rxn(model, non_zero_flux_D, 'true', ...
%     'struc', {'txdn_c'}, {'h_c','nad_c','nadh_c','h_e','h2o_c','h2o_e','h2o_p','h_p',...
%     'atp_c','adp_c','co2_c','co2_e','co2_p'}, FBA_solution_D);
[directionRxns_w, involvedMets_w, deadEnds_w] = draw_by_met(model, {'pyr_c'},...
    'true', 1, 'sub', {'iodine1_c'}, solution_1.x);
% max taxadien from this added reactions 2.06

% [minFlux, maxFlux] = fluxVariability(model, 60, 'max', {'IPDDI','DMATT','GRTT',...
%     'EX_mva_e','MVNOR','pksg','hmg1','mvak1','mvak2','mvad','ggpps','txs',...
%    'txdn_ex','Ex_txdn_e' })













