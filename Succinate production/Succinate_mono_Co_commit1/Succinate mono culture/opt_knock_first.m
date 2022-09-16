clc
clear
close all 
changeCobraSolver('gurobi','all');
model = readCbModel('e_coli_core.xml');
modelx = model;
biomass = 'BIOMASS_Ecoli_core_w_GAM';
threshold = 5;
selectedRxnList = modelx.rxns;
%% 
% Then, calculates the production of metabolites before running optKnock.

% determine succinate production and growth rate
lb = -100;
ub = 1000;
[modelx, EX_pep_e] = addExchangeRxn(modelx, 'pep_c', lb, ub);
modelx = changeRxnBounds(modelx, {'EX_glc__D_e','EX_o2_e'}, [-15,-20], 'l');
fbaWT = optimizeCbModel(modelx);
succFluxWT = fbaWT.x(strcmp(modelx.rxns, 'EX_pep_e'));
etohFluxWT = fbaWT.x(strcmp(modelx.rxns, 'EX_etoh_e'));
formFluxWT = fbaWT.x(strcmp(modelx.rxns, 'EX_for_e'));
lactFluxWT = fbaWT.x(strcmp(modelx.rxns, 'EX_lac__D_e'));
acetFluxWT = fbaWT.x(strcmp(modelx.rxns, 'EX_ac_e'));

growthRateWT = fbaWT.f;
fprintf('The production of succinate before optimization is %.1f \n', succFluxWT);
fprintf('The growth rate before optimization is %.1f \n', growthRateWT);
fprintf(['The production of other products such as ethanol, formate, lactate and'...
         'acetate are %.1f, %.1f, %.1f and %.1f, respectively. \n'], ...
        etohFluxWT, formFluxWT, lactFluxWT, acetFluxWT);
%% 
% *I) EXAMPLE 1 : SUCCINATE OVERPRODUCTION*
% 
% *Aim:* *finding optKnock reactions sets of size 2 for increasing production 
% of succinate*
%%
fprintf('\nFinding optKnock sets of size 3 or less...\n\n')
% Set optKnock options
% The exchange of succinate will be the objective of the outer problem
options = struct('targetRxn', 'EX_succ_e', 'numDel', 5);
% We will impose that biomass be at least 50% of the biomass of wild-type
constrOpt = struct('rxnList', {{biomass}},'values', 0.6*fbaWT.f, 'sense', 'G');
% We will try to find 10 optKnock sets of a maximun length of 2
previousSolutions = cell(10, 1);
contPreviousSolutions = 1;
nIter = 1;
while nIter < threshold
    fprintf('...Performing optKnock analysis...\n')
    if isempty(previousSolutions{1})
        optKnockSol = OptKnock(modelx, selectedRxnList, options, constrOpt);
    else
        optKnockSol = OptKnock(modelx, selectedRxnList, options, constrOpt, previousSolutions, 1);
    end
    
    % determine succinate production and growth rate after optimization
    succFluxM1 = optKnockSol.fluxes(strcmp(modelx.rxns, 'EX_pep_e'));
    growthRateM1 = optKnockSol.fluxes(strcmp(modelx.rxns, biomass));
    etohFluxM1 = optKnockSol.fluxes(strcmp(modelx.rxns, 'EX_etoh_e'));
    formFluxM1 = optKnockSol.fluxes(strcmp(modelx.rxns, 'EX_for_e'));
    lactFluxM1 = optKnockSol.fluxes(strcmp(modelx.rxns, 'EX_lac__D_e'));
    acetFluxM1 = optKnockSol.fluxes(strcmp(modelx.rxns, 'EX_ac_e'));
    setM1 = optKnockSol.rxnList;
    
    if ~isempty(setM1)
        previousSolutions{contPreviousSolutions} = setM1;
        contPreviousSolutions = contPreviousSolutions + 1;
        %printing results
        fprintf('optKnock found a optKnock set of large %d composed by ', length(setM1));
        for j = 1:length(setM1)
            if j == 1
                fprintf('%s', setM1{j});
            elseif j == length(setM1)
                fprintf(' and %s', setM1{j});
            else
                fprintf(', %s', setM1{j});
            end
        end
        fprintf('\n');
        fprintf('The production of succinate after optimization is %.2f \n', succFluxM1);
        fprintf('The growth rate after optimization is %.2f \n', growthRateM1);
        fprintf(['The production of other products such as ethanol, formate, lactate and acetate are' ...
                 '%.1f, %.1f, %.1f and %.1f, respectively. \n'], etohFluxM1, formFluxM1, lactFluxM1, acetFluxM1);
        fprintf('...Performing coupling analysis...\n');
%         [type, maxGrowth, maxProd, minProd] = analyzeOptKnock(model,setM1,'EX_etoh_e','BIOMASS_Ecoli_core_w_GAM');
%         fprintf('The solution is of type: %s\n', type);
%         fprintf('The maximun growth rate given the optKnock set is %.2f\n', maxGrowth);
%         fprintf(['The maximun and minimun production of succinate given the optKnock set is ' ...
%                  '%.2f and %.2f, respectively \n\n'], minProd, maxProd);
%         if strcmp(type, 'growth coupled')
%             singleProductionEnvelope(model, setM1, 'EX_etoh_e', biomass, 'savePlot', 1, 'showPlot', 1, ...
%                                      'fileName', ['etoh_ex1_' num2str(nIter)], 'outputFolder', 'OptKnockResults');
%         end
    else
        if nIter == 1
            fprintf('optKnock was not able to found an optKnock set\n');
        else
            fprintf('optKnock was not able to found additional optKnock sets\n');
        end
        break;
    end
    nIter = nIter + 1;
end

    

