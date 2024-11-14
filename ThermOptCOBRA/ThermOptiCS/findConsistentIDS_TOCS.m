function [reacInd,x] = findConsistentIDS_TOCS(model,core,TICmat,rev2irrev,tol)
% Identifies the thermodynamically consistent IDs for the given input of irreversible reaction model
% 
%
% USAGE: 
%   reacInd = findConsistentIDS_TOCS(model,core,TICmat,rev2irrev,tol)
%
% INPUTS:
%     model:     COBRA model structure with irreversible reactions only
%     core:      Reaction IDs for which thermodynamic consistency has to be verified
%     TICmat:    Matlab matrix with information on TICs in the model
%     rev2irrev: Mapping between the reversible reactions to its corresponding irreversible reactions
%     tol:       Minimum positive number defined as non-zero
%
% OUTPUTS:
%     reacInd:   Indices of reactions that are identified to be thermodynamically consistent
%     x:         The  solution obtained after solving the MILP optimization problem
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras
[m,n] = size(model.S);
lbs= max([tol*ones(n,1),model.lb],[],2);
% objective
f = [zeros(n,1);double(ismember([1:n]',core))];

% equalities
Aeq = [model.S, sparse(m,n)];
beq = zeros(m,1);
csenseeq = repmat('E',m,1); % equality

% inequalities
temp1 = speye(n);
temp2 = -1*spdiag(lbs);
Aineq1 = [temp1,temp2];
bineq1 = zeros(n,1);
csenseineq1 = repmat('G',n,1); % greater than

temp2 = -1*spdiag(model.ub);
Aineq2 = [temp1,temp2];
bineq2 = zeros(n,1);
csenseineq2 = repmat('L',n,1); % lesser than
nTICs =size(TICmat,1);

% TIC constraints
Aineq3=[sparse(nTICs,n), TICmat];
bineq3=sum(TICmat,2)-1;
csenseineq3=repmat('L',nTICs,1);

Aineq4=[];bineq4=[];csenseineq4=[];
%%% type 1
for i=1:numel(rev2irrev)
    if numel(rev2irrev{i})>=2
        temp1 = sparse(1,n);
        temp2 = sparse(1,n); temp2(rev2irrev{i}) = ones(numel(rev2irrev{i}),1);

        Aineq4 = [Aineq4;temp1,temp2];
        bineq4 = [bineq4;1];
        csenseineq4 = [csenseineq4;'L']; % lesser than
    end
end


% bounds
lb = [model.lb;zeros(n,1)];
ub = [model.ub;ones(n,1)];

% Set up MILP problem
MILPproblem.A=[Aeq;Aineq1;Aineq2;Aineq3;Aineq4];
MILPproblem.b=[beq;bineq1;bineq2;bineq3;bineq4];
MILPproblem.lb=lb;
MILPproblem.ub=ub;
MILPproblem.c=f;
MILPproblem.osense=-1;%maximise
MILPproblem.vartype = [repmat('C',n,1);repmat('B',n,1)];
MILPproblem.csense = [csenseeq; csenseineq1; csenseineq2; csenseineq3; csenseineq4];
solution = solveCobraMILP(MILPproblem);
stat=solution.stat;
if stat~=1
    x = [];
    reacInd = [];
else
    x=solution.full;    
    reacInd = find(x(n+1:end)>0.5);
end
end

