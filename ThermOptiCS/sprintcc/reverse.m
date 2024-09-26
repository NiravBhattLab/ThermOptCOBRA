function [flux,z] = reverse(model,core,tol)
% USAGE:
%   [flux,z] = reverse(model,core,tol)
%
% INPUTS:
%    model: COBRA model structure.
%    core:  A binary vector indicating whether a core or non-core reaction
%    tol:   tolerance level (minimum absolute flux that has to be carried
%           by a reaction for it to be defined as consistent)
%
% OUTPUTS:
%    flux:  Flux through the reactions in the COBRA model
%    z:     Solution for the decision variables
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras

[m,n] = size(model.S);
rev_core = model.rev & core; % reversible core reactions
n_rev_core = sum(rev_core);% number of reversible core reactions


% objective
f = [zeros(n,1);unifrnd(1,1.1,n_rev_core,1)];

% equalities
Aeq = [model.S, sparse(m,n_rev_core)];
beq = zeros(m,1);
csenseeq = repmat('E',m,1); % equality

% inequalities
temp1 = speye(n);
temp2 = speye(n_rev_core);
Aineq = [temp1(rev_core,:),-1*temp2 ];
bineq = zeros(n_rev_core,1);
csenseineq = repmat('L',n_rev_core,1); % lesser than

% bounds
lb = model.lb;
ub = model.ub;
lb = [lb;-tol*ones(n_rev_core,1)];
ub = [ub; Inf(n_rev_core,1)];

% Set up LP problem
LPproblem.A=[Aeq;Aineq];
LPproblem.b=[beq;bineq];
LPproblem.lb=lb;
LPproblem.ub=ub;
LPproblem.c=f;
LPproblem.osense=1;%minimise
LPproblem.csense = [csenseeq; csenseineq];

solution = solveCobraLP(LPproblem);
if solution.stat~=1
    fprintf('%s%s\n',num2str(solution.stat),' = solution.stat')
    fprintf('%s%s\n',num2str(solution.origStat),' = solution.origStat')
    warning('LP solution may not be optimal')
end

x=solution.full;
if ~isempty(x)
    flux = x(1:n);
    z = x(n+1:end);
else
    flux=nan(n,1);
    z=nan(n_rev_core,1);
end

