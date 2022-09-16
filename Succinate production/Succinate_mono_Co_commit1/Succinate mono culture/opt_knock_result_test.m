clc
clear
close all 
changeCobraSolver('gurobi','all');
model = readCbModel('e_coli_core.xml');
modelx = model;
%% optKnock found a optKnock set of large 3 composed by ATPS4r, PYRt2 and CO2t
selectedRxnList = modelx.rxns;
changeObjective(modelx,'EX_succ_e');
modelx = changeRxnBounds(modelx, {'EX_glc__D_e','EX_o2_e'}, [-15,0], 'l');
fba = optimizeCbModel(modelx);
printFluxVector(modelx,fba.x,true,true,-1,[],[],true)
% %% engineered
% geneList = {'b0903','b2579','b0902','b0902','b3951','b2133','b1380',...
%     'b0351','b1241'};
modelx = singleRxnDeletion(modelx,{'ATPS4r','PYRt2','CO2t'});
% modelx = deleteModelGenes(modelx,geneList);
fba_optknock = optimizeCbModel(modelx);
printFluxVector(modelx,fba_optknock.x,true,true,-1,[],[],true)
