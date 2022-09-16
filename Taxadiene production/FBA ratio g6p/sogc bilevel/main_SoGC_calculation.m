% In the name of GOD

clc
clear
close all
format short g

solverOK = changeCobraSolver('gurobi','all');

%% read model
model = readCbModel('iJO1366.mat');
% load('e_coli_core','e_coli_core')
% % writeCbModel(e_coli_core, 'xls', 'e_coli_core')
% model = e_coli_core;

%% Gene Deletion
[num,txt,raw] = xlsread('ALL_triple_EX_succ_e_FVA.xlsx');

%% SoGC
for i=2:length(txt)
    A = txt{i,1};
    B = strsplit(A,',');
    DEL{i-1,1} = B;
end
environment = getEnvironment();
parfor (i = 1:length(txt)-1, 8)
    restoreEnvironment(environment);
    changeCobraSolver('gurobi','LP',0,-1);    
    model_1 = removeRxns(model,DEL{i});
    model_1 = changeRxnBounds(model_1, {'EX_glc__D_e','EX_o2_e'}, [-15, -18],'l');
    [w2,uu2] = SoGC_calculator(model_1,{'EX_succ_e'},{'BIOMASS_Ec_iJO1366_core_53p95M'});
    W2(i,1) = w2;
    UU2(i,1) = uu2;
end

T = table(deletions,W2,UU2);
writetable(T,'SoGC_3_gene_aerobic_succinate.xls')