rxnIDs = findRxnIDs(model,{'PGMT','PGI','PGMT','G6PDH2r','GALUi',...
    'MLTP3','GALUi','MLTP1','GALUi','MLTP2','SALID','TRE6PS'});
b1 = 0.54;
b2 = 0.99;
b3 = 0.20;

S = model.S;

ss = zeros(1,length(model.rxns));
ss(1,rxnIDs(1:2)) = [1-b1 -b1];
S(end+1,:) = ss;
ss = zeros(1,length(model.rxns));
ss(1,rxnIDs(3:4)) = [1-b1 -b1];
S(end+1,:) = ss;
ss = zeros(1,length(model.rxns));
ss(1,rxnIDs(5:6)) = [1-b2 -b2];
S(end+1,:) = ss;
ss = zeros(1,length(model.rxns));
ss(1,rxnIDs(7:8)) = [1-b2 -b2];
S(end+1,:) = ss;
ss = zeros(1,length(model.rxns));
ss(1,rxnIDs(9:10)) = [1-b2 -b2];
S(end+1,:) = ss;
% ss = zeros(1,length(model.rxns));
% ss(1,rxnIDs(11:12)) = [1-b3 -b3];
% S(end+1,:) = ss;
model.S = S;
