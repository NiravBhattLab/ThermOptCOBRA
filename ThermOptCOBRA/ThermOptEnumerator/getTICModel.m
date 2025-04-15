function [allM1,blkdCore,allFlux,stat] = getTICModel(model,core,tol,dir,TICcons)
% Obtains the thermodynamically infeasible cycle that has the core reaction in it
% 
% USAGE: 
%   [m1,blkdCore,flux,stat] = getTICModel(model,core,tol,dir,TICcons)
%
% INPUTS:
%     model:     COBRA model structure
%     core:      Reaction ID for which the TIC has to be identified
%     tol:       Minimum positive number defined as non-zero
%     dir:       1/-1 to define the flux direction of the core reaction. 1: positive flux
%                -1: negative flux
%     TICcons:   Details on the previously identified TICs that involves the core reaction
%
% OUTPUTS:
%     allM1:     A cell of TICs obtained
%     blkdCore:  1/0 to define consistency of the core reaction
%     allFlux:   Flux though the reactions in the TICs in allM1
%     stat:      Status of the MILP optimization
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras

[~,n] = size(model.S); % number of metabolites and reactions

direction = zeros(n,1);
if dir ==1
    % for forward direction
    direction(core)=1;
    [reacInd,~,stat] = findConsistentReacIDTIC_feas(model,direction,tol,TICcons);
elseif dir==-1
    % for reverse direction
    direction(core)=-1;
    [reacInd,~,stat] = findConsistentReacIDTIC_feas(model,direction,tol,TICcons);
end

if stat~=1
    blkdCore=1;allFlux={};allM1={};
else
    blkdCore=0;
    model2 =removeRxns(model,setdiff(model.rxns,model.rxns(reacInd)));
    direction2 = zeros(numel(model2.rxns),1);
    direction2(ismember(model2.rxns,model.rxns(core)))=dir;
    TICcons2 = getNewTICcons(reacInd,TICcons,model,model2);
    [reacInd2,x,stat] = findConsistentReacIDTIC(model2,direction2,tol,TICcons2);
    allM1{1,1} = model2.rxns(reacInd2); allFlux{1,1} = x(reacInd2);
    TICcons2{end+1,1} = find(reacInd2);
    [allM1_,blkdCore_,allFlux_] = getTICModel(model2,ismember(model2.rxns,model.rxns(core)),tol,dir,TICcons2);
    if blkdCore_~=1
        allM1 = [allM1;allM1_]; allFlux = [allFlux; allFlux_];
    end
end
end
function TICcons2 = getNewTICcons(reacInd,TICcons,model,model2)
k=1;
reacInd = find(reacInd);
TICcons2={};
for i=1:numel(TICcons)
    if sum(ismember(TICcons{i,1},reacInd))==numel(TICcons{i,1})
        TICcons2{k,1} = find(ismember(model2.rxns,model.rxns(TICcons{i,1})));
        k=k+1;
    end
end
end

function [reacInd,x,stat] = findConsistentReacIDTIC(model,direction,tol,TICcons)

[m,n] = size(model.S);
dir0 = direction==0;
n_ = sum(dir0);
% objective
f = [zeros(n,1);ones(n_,1)];

% equalities
Aeq = [model.S, sparse(m,n_)];
beq = zeros(m,1);
csenseeq = repmat('E',m,1); % equality

% inequalities
temp1 = speye(n);
temp2 = -1*spdiag(model.lb(dir0));
Aineq1 = [temp1(dir0,:),temp2];
bineq1 = zeros(n_,1);
csenseineq1 = repmat('G',n_,1); % greater than

temp2 = -1*spdiag(model.ub(dir0));
Aineq2 = [temp1(dir0,:),temp2];
bineq2 = zeros(n_,1);
csenseineq2 = repmat('L',n_,1); % lesser than


Aineq3=[];bineq3=[];csenseineq3=[];
if ~isempty(TICcons)
    for i=1:size(TICcons,1)
        temp1 = sparse(1,n);
        Tids = TICcons{i,1};
        temp2 = sparse(1,n); temp2(Tids) = ones(numel(Tids),1);
        temp2 = temp2(dir0);
        Aineq3 = [Aineq3;temp1,temp2];
        bineq3 = [bineq3;sum(temp2)-1];
        csenseineq3 = [csenseineq3;'L']; % lesser than
    end
end

% bounds
lb = model.lb;
lb(direction==1)=max([tol*ones(sum(direction==1),1),lb(direction==1)],[],2);
lb = [lb;zeros(n_,1)];
ub = model.ub;
ub(direction==-1)=-tol*ones(sum(direction==-1),1);
ub = [ub;ones(n_,1)];

% Set up LP problem
MILPproblem.A=[Aeq;Aineq1;Aineq2;Aineq3];
MILPproblem.b=[beq;bineq1;bineq2;bineq3];
MILPproblem.lb=lb;
MILPproblem.ub=ub;
MILPproblem.c=f;
MILPproblem.osense=1;%minimise
MILPproblem.vartype = [repmat('C',n,1);repmat('B',n_,1)];
MILPproblem.csense = [csenseeq; csenseineq1; csenseineq2; csenseineq3];
solution = solveCobraMILP(MILPproblem);
stat=solution.stat;
if stat~=1
    x = [];
    reacInd = [];
else
    x=solution.full;
    reacInd = abs(x(1:n))>=tol*1e-3;
end
end


function [reacInd,x,stat] = findConsistentReacIDTIC_feas(model,direction,tol,TICcons)

[m,n] = size(model.S);
dir0 = direction==0;
n_ = sum(dir0);
% objective
f = [zeros(n+n_,1)];

% equalities
Aeq = [model.S, sparse(m,n_)];
beq = zeros(m,1);
csenseeq = repmat('E',m,1); % equality

% inequalities
temp1 = speye(n);
temp2 = -1*spdiag(model.lb(dir0));
Aineq1 = [temp1(dir0,:),temp2];
bineq1 = zeros(n_,1);
csenseineq1 = repmat('G',n_,1); % greater than

temp2 = -1*spdiag(model.ub(dir0));
Aineq2 = [temp1(dir0,:),temp2];
bineq2 = zeros(n_,1);
csenseineq2 = repmat('L',n_,1); % lesser than


Aineq3=[];bineq3=[];csenseineq3=[];
if ~isempty(TICcons)
    for i=1:size(TICcons,1)
        temp1 = sparse(1,n);
        Tids = TICcons{i,1};
        temp2 = sparse(1,n); temp2(Tids) = ones(numel(Tids),1);
        temp2 = temp2(dir0);
        Aineq3 = [Aineq3;temp1,temp2];
        bineq3 = [bineq3;sum(temp2)-1];
        csenseineq3 = [csenseineq3;'L']; % lesser than
    end
end

% bounds
lb = model.lb;
lb(direction==1)=max([tol*ones(sum(direction==1),1),lb(direction==1)],[],2);
lb = [lb;zeros(n_,1)];
ub = model.ub;
ub(direction==-1)=-tol*ones(sum(direction==-1),1);
ub = [ub;ones(n_,1)];

% Set up LP problem
MILPproblem.A=[Aeq;Aineq1;Aineq2;Aineq3];
MILPproblem.b=[beq;bineq1;bineq2;bineq3];
MILPproblem.lb=lb;
MILPproblem.ub=ub;
MILPproblem.c=f;
MILPproblem.osense=1;%minimise
MILPproblem.vartype = [repmat('C',n,1);repmat('B',n_,1)];
MILPproblem.csense = [csenseeq; csenseineq1; csenseineq2; csenseineq3];
solution = solveCobraMILP(MILPproblem);
stat=solution.stat;
if stat~=1
    x = [];
    reacInd = [];
else
    x=solution.full;
    reacInd = abs(x(1:n))>=tol*1e-3;
end
end