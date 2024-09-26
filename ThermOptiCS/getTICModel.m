function [m1,blkdCore,flux,stat] = getTICModel(model,core,tol,dir,TICcons)
[~,n] = size(model.S); % number of metabolites and reactions

direction = zeros(n,1);
if dir ==1
    % for forward direction
    direction(core)=1;
    [reacInd,x,stat] = findConsistentReacIDTIC(model,direction,tol,TICcons);
elseif dir==-1
    % for reverse direction
    direction(core)=-1;
    [reacInd,x,stat] = findConsistentReacIDTIC(model,direction,tol,TICcons);
end

if stat~=1
    blkdCore=1;flux={};m1={};
else
    blkdCore=0;
    flux = x(1:numel(model.rxns));
    m1 = model.rxns(reacInd);
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