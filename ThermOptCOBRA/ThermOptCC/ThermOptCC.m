function [a,TICs,Dir,modModel,TIC_Rxns] = ThermOptCC(model,tol,TICs,Dir)
% Identifies thermodynamically feasible flux directions for all the
% reactions in the input model
%
% USAGE: 
%   [a,TICs,Dir] = ThermOptiCC(model,tol)
%
% INPUTS:
%     model:     COBRA model structure for which thermodynamic feasibility
%                of the reactions has to be identified
%     tol:       Tolerance value (User defined non-zero value).
%
% OPTIONAL INPUTS:
%     TICs:       List of all the Thermodynamically infeasible cycles in
%                 the given input model
%     Dir:        The flux directions for reactions in the corresponding
%                 TICs
% OUTPUTS:
%     a:          A cell describing the thermodynamically feasible direction 
%                 of the reactions in the given input model                 
%     TICs:       List of all the Thermodynamically infeasible cycles in
%                 the given input model
%     Dir:        The flux directions for reactions in the corresponding
%                 TICs
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras


% checking for reactions irreversible in reverse direction
IrR = model.ub<=0; % reactions irreversible in reverse direction
temp = model.ub(IrR);
model.S(:,IrR) = -model.S(:,IrR);
model.ub(IrR) = -1*model.lb(IrR);
model.lb(IrR) = -1*temp;
modModel = model;
% normalising the bounds to max of 1000
temp1 = max(abs(model.lb));
temp2 = max(abs(model.ub));
model.lb(abs(model.lb)==temp1) = sign(model.lb(abs(model.lb)==temp1))*10000;
model.ub(abs(model.ub)==temp2) = sign(model.ub(abs(model.ub)==temp2))*10000;
if ~exist('TICs', 'var') || isempty(TICs) || ~exist('Dir', 'var') || isempty(Dir)
    [TICs,Dir,TIC_Rxns] = ThermOptEnumMILP(model);
end

% the model should not have any reaction with flux only in reverse
% direction
[modelIrrev, ~, rev2irrev] = convertToIrreversible(model);
[~,n] = size(modelIrrev.S);
% creating a TICmat that has info about which reactions (irreversible)
% participate in a TIC
TICmat = zeros(numel(TICs),n);
for i=1:numel(TICs)
    cTIC = TICs{i};cDir = Dir{i};
    cDir(cDir>0)=1;cDir(cDir<0)=-1;
    temp = rev2irrev(ismember(model.rxns,cTIC));
    ids=cellfun(@(x)temp{x}([1;-1]==cDir(x)),num2cell([1:numel(temp)]'));
    TICmat(i,ids)=1;
end

core = 1:n; newCore = [];
while numel(core)~=numel(newCore)
    newCore = core;
    IDS = findConsistentIDS(modelIrrev,core,TICmat,rev2irrev,tol);
    core = setdiff(core,IDS);
end
consIrr = double(~ismember([1:n],core)); % consistent irreversible reactions
a = cellfun(@(x)getConsDir(consIrr,x),rev2irrev,'UniformOutput',0);
end
function consDir = getConsDir(consIrr,x)
temp = consIrr(x);
if sum(temp)==2
    consDir = 'Reversible'; % consistent in both the direction
elseif sum(temp)==0
    consDir = 'Blocked'; % blocked reaction
elseif temp(1)==1
    consDir = 'Forward'; % consistent only in forward direction
elseif temp(2)==1
    consDir = 'Reverse'; % consistent only in reverse direction
end
end

function reacInd = findConsistentIDS(model,core,TICmat,rev2irrev,tol)
[m,n] = size(model.S);
lbs = zeros(n,1);
lbs(core) = max([tol*ones(numel(core),1),model.lb(core)],[],2);
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

for i=1:numel(rev2irrev)
    temp1 = sparse(1,n);
    temp2 = sparse(1,n); temp2(rev2irrev{i}) = ones(numel(rev2irrev{i}),1);
    
    Aineq4 = [Aineq4;temp1,temp2];
    bineq4 = [bineq4;1];
    csenseineq4 = [csenseineq4;'L']; % lesser than
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
    reacInd = intersect(find(x(n+1:end)>=0.5),core);
end
end
