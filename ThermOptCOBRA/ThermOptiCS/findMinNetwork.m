function IDS = findMinNetwork(model,core,P,TICmat,rev2irrev,tol,x0,timeLimit)
% Identifies the thermodynamically consistent minimal network
% 
% USAGE: 
%   IDS = findMinNetwork(model,core,P,TICmat,rev2irrev,tol,x0,timeLimit)
%
% INPUTS:
%     model:     COBRA model structure with irreversible reactions only
%     core:      Reaction IDs for which thermodynamic consistency has been verified
%     TICmat:    Matlab matrix with information on TICs in the model
%     rev2irrev: Mapping between the reversible reactions to its corresponding irreversible reactions
%     tol:       Minimum positive number defined as non-zero
%     x0:        Inital feasible point for the MILP problem
%     timeLimit: timeLimit for solving the minNetwork MILP problem
%
% OUTPUTS:
%     IDS:       Indices of reactions in the minimal network
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras


[m,n] = size(model.S);
% objective
temp = zeros(n,1);
temp(P) = 1;
lbs = max([tol*ones(numel(model.lb),1),model.lb],[],2);
f = [zeros(n,1);temp];

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

%%% type2
% rBlk=[];
% for i=1:numel(rev2irrev)
%     tp=rev2irrev{i};
%     if numel(tp)>=2
%         tp2 = tp(~ismember(tp,core))';
%         if numel(tp2)==1
%             rBlk=[rBlk;tp2];
%         end
%     end
% end

% bounds
lb = model.lb; lb(core)=tol;
lb = [lb;zeros(n,1)];
ub = [model.ub;ones(n,1)];

% Set up MILP problem
MILPproblem.A=[Aeq;Aineq1;Aineq2;Aineq3;Aineq4];
MILPproblem.b=[beq;bineq1;bineq2;bineq3;bineq4];
MILPproblem.lb=lb;
MILPproblem.ub=ub;
MILPproblem.c=f;
MILPproblem.osense=1; 
MILPproblem.vartype = [repmat('C',n,1);repmat('B',n,1)];
MILPproblem.csense = [csenseeq; csenseineq1; csenseineq2; csenseineq3;csenseineq4];
if ~exist('x0', 'var') || isempty('x0')
    MILPproblem.x0 = [];
else
    MILPproblem.x0  = x0;
end

solution = solveCobraMILP(MILPproblem,'timeLimit',timeLimit);
x=solution.full;
IDS = find(abs(x(1:n))>=1e-12 | x(n+1:end)>0.5);
end
