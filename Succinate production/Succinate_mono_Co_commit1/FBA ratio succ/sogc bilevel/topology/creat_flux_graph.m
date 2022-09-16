function [flux_model,mets_index,rxns_index] = creat_flux_graph(model_1,Involved_mets,non_zero_flux)
%% creat the graph of the flux

S_reference = model_1.S;
lb_reference = model_1.lb;
ub_reference = model_1.ub;
mets_reference = model_1.mets;
rxns_reference = model_1.rxns;
rules_reference = model_1.rules;


for i=1:length(Involved_mets)
    mets_index(i,1) = find(strcmp(mets_reference,Involved_mets{i}));
end
for i=1:length(non_zero_flux)
    rxns_index(i,1) = find(strcmp(rxns_reference,non_zero_flux{i}));
    lb(i,1) = lb_reference(rxns_index(i,1));
    ub(i,1) = ub_reference(rxns_index(i,1));
end

% G = zeros(length(Involved_mets),length(non_zero1_flux));
% S = zeros(length(Involved_mets),length(non_zero1_flux));
S = S_reference(mets_index,rxns_index);

field1 = 'S';  value1 = {S};
field2 = 'mets';  value2 = {Involved_mets};
field3 = 'lb';  value3 = {lb};
field4 = 'ub';  value4 = {ub};
field5 = 'rxns';  value5 = {non_zero_flux};



flux_model = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5)
end