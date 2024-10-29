% This code gets the count on number of reactions that are involved in TICs
% in each of the bigg models
clear
p=dir('./BiggModelsMat/'); % path to bigg models
p={p(3:end).name}';
tbl = table('Size',[0,3],'VariableTypes',{'string','double','double'},...
    'VariableNames',{'Model','No_of_reactions','No_of_reactions_in_TICs'});
for i=1:numel(p)
    load(['./BiggModelsMat/',p{i}])
    a=getTICrxnCount(model);
    tbl(end+1,:)={strrep(p{i},'.mat',''),numel(model.rxns),numel(a)};
    clear model
end
tbl.Fraction=tbl.No_ofreactions_in_TICs./tbl.No_of_reactions;
function a= getTICrxnCount(model)
tol=1;
[~,n] = size(model.S); % number of metabolites and reactions
IrR = model.ub<=0; % reactions irreversible in reverse direction
temp = model.ub(IrR);
model.S(:,IrR) = -model.S(:,IrR);
model.ub(IrR) = -1*model.lb(IrR);
model.lb(IrR) = -1*temp;
rev = ones(n,1);
rev(model.lb>=0) = 0;
model.rev=rev;
model.lb(model.lb>0)=0;
ind=findExcRxns(model);
% blocking all the exchange reactions
model.lb(ind) = 0; model.ub(ind) = 0;
a = sprintcc(model,tol); % reactions that are in TICs
end