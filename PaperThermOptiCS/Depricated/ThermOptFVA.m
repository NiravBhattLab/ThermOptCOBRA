function [minFlux, maxFlux] = ThermOptFVA(model,optPercentage,rxnNameList,TICs,Dir)

IrR = model.ub<=0; % reactions irreversible in reverse direction
temp = model.ub(IrR);
model.S(:,IrR) = -model.S(:,IrR);
model.ub(IrR) = -1*model.lb(IrR);
model.lb(IrR) = -1*temp;

% normalising the bounds to max of 1000
temp1 = max(abs(model.lb));
temp2 = max(abs(model.ub));
model.lb(abs(model.lb)==temp1) = sign(model.lb(abs(model.lb)==temp1))*1000;
model.ub(abs(model.ub)==temp2) = sign(model.ub(abs(model.ub)==temp2))*1000;
if ~exist('TICs', 'var') || isempty(TICs) || ~exist('Dir', 'var') || isempty(Dir)
    [TICs,Dir,~] = ThermoOptEnumMILP(model);
end
if ~exist('rxnNameList', 'var') || isempty(rxnNameList)
    rxnNameList=model.rxns;
end
% performing FBA on the given input model
sol = optimizeCbModel(model);
model.lb(find(model.c))=optPercentage*sol.x(find(model.c));
TICmat = zeros(numel(TICs),numel(model.rxns));
for i=1:numel(TICs)
TICmat(i,ismember(model.rxns,TICs{i})) = sign(Dir{i});
end
minFlux=[]; maxFlux=[];
for i=1:numel(rxnNameList)
    r=rxnNameList{i};
    id = find(ismember(model.rxns,r));
    Trxns = find(sum(TICmat(find(TICmat(:,id)),:),1));
    %minimizing with ll constraints
    sol = getLLFBA(model,id,1,Trxns);
    minFlux(i)=sol;
    %maximizing with ll constraints
    sol = getLLFBA(model,id,-1,Trxns);
    maxFlux(i)=sol;
end
end
function flux = getLLFBA(model,rxn,opt,Trxns)
[m,n] = size(model.S);
nSpace = null(full(model.S(:,Trxns)));
nT = numel(Trxns);

% objective
f = zeros(n+2*nT,1);
f(rxn)=1;

% equalities
Aeq = [model.S, sparse(m,2*nT)];
beq = zeros(m,1);
csenseeq = repmat('E',m,1); % equality

if nT==0
    Aineq1 = [];
    bineq1 = [];
    csenseineq1 = []; 
    Aineq2 = [];
    bineq2 = [];
    csenseineq2 = []; 
    Aineq3 = [];
    bineq3 = [];
    csenseineq3 = []; 
    Aineq4 = [];
    bineq4 = [];
    csenseineq4 = []; 
    Aineq5 = [];
    bineq5 = [];
    csenseineq5 = []; 
else
% inequalities
temp1 = speye(n);
temp1 = temp1(Trxns,:);
temp2 = sparse(nT);
temp3 = -1000*speye(nT);
Aineq1 = [temp1,temp2,temp3];
bineq1 = -1000*ones(nT,1);
csenseineq1 = repmat('G',n,1); 

% inequalities
temp1 = speye(n);
temp1 = temp1(Trxns,:);
temp2 = sparse(nT);
temp3 = -1000*speye(nT);
Aineq2 = [temp1,temp2,temp3];
bineq2 = zeros(nT,1);
csenseineq2 = repmat('L',n,1); 

% inequalities
temp1 = sparse(nT,n);
temp2 = speye(nT);
temp3 = (1000+1e-7)*speye(nT);
Aineq3 = [temp1,temp2,temp3];
bineq3 = 1e-7*ones(nT,1);
csenseineq3 = repmat('G',n,1); 

% inequalities
temp1 = sparse(nT,n);
temp2 = speye(nT);
temp3 = (1000+1e-7)*speye(nT);
Aineq4 = [temp1,temp2,temp3];
bineq4 = 1000*ones(nT,1);
csenseineq4 = repmat('L',n,1); 


% equalities
temp1 = sparse(size(nSpace',1),n);
temp2 = sparse(size(nSpace',1),n);
temp3 = nSpace';
Aineq5 = [temp1,temp2,temp3];
bineq5 = zeros(nT,1);
csenseineq5 = repmat('E',n,1); 
end

% bounds
lb = model.lb;
lb = [lb;zeros(2*nT,1)];
ub = [model.ub;1000*ones(nT,1);ones(nT,1)];

% Set up MILP problem
MILPproblem.A=[Aeq;Aineq1;Aineq2;Aineq3;Aineq4;Aineq5];
MILPproblem.b=[beq;bineq1;bineq2;bineq3;bineq4;bineq5];
MILPproblem.lb=lb;
MILPproblem.ub=ub;
MILPproblem.c=f;
MILPproblem.osense=opt;
MILPproblem.vartype = [repmat('C',n+nT,1);repmat('B',nT,1)];
MILPproblem.csense = [csenseeq; csenseineq1; csenseineq2; csenseineq3; csenseineq4;csenseineq5];
solution = solveCobraMILP(MILPproblem);
x=solution.full;
flux = x(rxn);
end






