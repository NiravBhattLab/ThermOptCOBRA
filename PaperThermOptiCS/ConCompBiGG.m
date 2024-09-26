% This code gets the result on number of connected components in the BiGG
% models and stores in the folder ConnectedCompsResults. 
clear
A = importdata('BiggModels.txt');
for i=1:numel(A)
    load(A{i}) % load the Bigg model
    [CC,sizeCC]=getCC(model);
    save(['./ConnectedCompsResults/',A{i}],'CC','sizeCC')
    clear model CC sizeCC
end
function [CC,sizeCC] =  getCC(model)
[~,n] = size(model.S); % number of metabolites and reactions
IrR = model.ub<=0; % reactions irreversible in reverse direction
temp = model.ub(IrR);
model.S(:,IrR) = -model.S(:,IrR);
model.ub(IrR) = -1*model.lb(IrR);
model.lb(IrR) = -1*temp;
rev = ones(n,1);
rev(model.lb>=0) = 0;
model.rev=rev;
% converting all the positive lower bounds to zero lower bounds
model.lb(model.lb>0)=0;

% normalising the bounds to max of 1000
temp1 = max(abs(model.lb));
temp2 = max(abs(model.ub));
model.lb(abs(model.lb)==temp1) = sign(model.lb(abs(model.lb)==temp1))*1000;
model.ub(abs(model.ub)==temp2) = sign(model.ub(abs(model.ub)==temp2))*1000;

tol=1;
ind=findExcRxns(model);
% blocking all the exchange reactions
model.lb(ind) = 0; model.ub(ind) = 0;
a = sprintcc(model,tol); % reactions that are in TICs
% model with only TICs
model = removeRxns(model,model.rxns(setdiff([1:numel(model.rxns)],a)));
S = model.S~=0; S=logical(S'*S);
g = graph(S);
[~,CC]=conncomp(g);
sizeCC = numel(CC);
end