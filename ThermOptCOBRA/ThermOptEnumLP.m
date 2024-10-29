function [TICs, Direction, TIC_Rxns,modModel] = ThermOptEnumLP(model)
% Enumerates all the Thermodynamically infeasible cycles in a given model
%
% USAGE: 
%   [TICs,Direction,TIC_Rxns,modModel] = ThermoOptEnumLP(model)
%
% INPUTS:
%     model:     COBRA model structure for which TICs has be found
%
% OUTPUTS:
%     TICs:       List of all the Thermodynamically infeasible cycles in
%                 the given input model
%     Direction:  Relative flux coefficients for reactions 
%                 in the corresponding TICs
%     TIC_Rxns:   Reaction list that participates in the TICs
%     modModel:   Modified model that has no irreversible reactions that
%                 carry flux in reverse direction
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras


% checking for reactions irreversible in reverse direction
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

modModel = model;
tol=1;
ind=findExcRxns(model);
% blocking all the exchange reactions
model.lb(ind) = 0; model.ub(ind) = 0;
if exist('sprintcc','file')
    a = sprintcc(model,tol); % reactions that are in TICs
else
    a = fastcc(model,tol);
end
if isempty(a)
    TICs={};Direction={};TIC_Rxns={};opt=1;
    return
end
TIC_Rxns = model.rxns(a);
% model with only TICs
model = removeRxns(model,model.rxns(setdiff([1:numel(model.rxns)],a)));

order = getOrderOfRxns(model);

TICs={};Direction={};
while ~isempty(model.rxns)
    % for flux in positive direction
    core=order(1);blkd=[];dir=1;
    [T,D]= getAllCurrTIC(model,core,blkd,dir,tol,{});
    TICs = [TICs;T]; Direction = [Direction;D];
    % for flux in negative direction
    core=order(1);blkd=[];dir=-1;
    [T,D]= getAllCurrTIC(model,core,blkd,dir,tol,{});
    TICs = [TICs;T]; Direction = [Direction;D];
    
    model.lb(order(1))=0; model.ub(order(1))=0;
    if exist('sprintcc','file')
        a = sprintcc(model,tol); 
    else
        a = fastcc(model,tol);
    end

    modelTemp = removeRxns(model,model.rxns(setdiff([1:numel(model.rxns)],a)));
    order = order(ismember(order,a));
    [~,order] =ismember(model.rxns(order),modelTemp.rxns);
    model = modelTemp;
end
end

function [currTIC,currDir,cache] = getAllCurrTIC(model,core,blkd,dir,tol,cache)
currTIC ={};currDir={};
blkdstr =regexprep(num2str(sort(blkd')),'\s+',',');
if ismember(blkdstr,cache)
    return
else
    cache =[cache;blkdstr];
end
modelTemp = model;
modelTemp.lb(blkd)=0;modelTemp.ub(blkd)=0;
if dir==1
    modelTemp.lb(core)=1;
elseif dir==-1
    modelTemp.ub(core)=-1;
end
[reacInd,stat,flux] = findConsistentReacIDTIC(modelTemp,tol);
if stat==1
    currTIC=[currTIC;reacInd];currDir = [currDir;flux(reacInd)];
    reacInd1 = setdiff(reacInd,core);
    for i=1:numel(reacInd1)
        blkd1=[blkd;reacInd1(i)];
        [currTIC1,currDir1,cache] = getAllCurrTIC(model,core,blkd1,dir,tol,cache);
        currTIC=[currTIC;currTIC1];currDir = [currDir;currDir1];
    end    
end
if isempty(blkd)
    [currTIC,ID] = getUniqueArrays(currTIC);
    currDir = currDir(ID);
    for i=1:numel(currTIC)
        currTIC{i}=model.rxns(currTIC{i});
    end
end
end

function [reacInd,stat,flux] = findConsistentReacIDTIC(model,tol)

[m,n] = size(model.S);

% objective
f = [zeros(n,1);ones(n,1)];

% equalities
Aeq = [model.S, sparse(m,n)];
beq = zeros(m,1);
csenseeq = repmat('E',m,1); % equality

% inequalities
temp1 = speye(n);
temp2 = speye(n);
Aineq1 = [temp1,temp2];
bineq1 = zeros(n,1);
csenseineq1 = repmat('G',n,1); % greater than

Aineq2 = [temp1,-1*temp2];
bineq2 = zeros(n,1);
csenseineq2 = repmat('L',n,1); % lesser than

% bounds
lb = model.lb;ub = model.ub;
lb = [lb;zeros(n,1)];
ub = [ub;Inf(n,1)];

% Set up LP problem
LPproblem.A=[Aeq;Aineq1;Aineq2];
LPproblem.b=[beq;bineq1;bineq2];
LPproblem.lb=lb;
LPproblem.ub=ub;
LPproblem.c=f;
LPproblem.osense=1;%minimise
LPproblem.csense = [csenseeq; csenseineq1; csenseineq2];


solution = solveCobraLP(LPproblem);

stat=solution.stat;
if stat~=1
    x = []; flux =[];
    reacInd = [];
else
    x=solution.full;flux = x(1:n);
    reacInd = find(abs(flux)>=tol*1e-3);
    
end
end

function [uniqueCellArray,ID] = getUniqueArrays(cellArray)
  % Initialize an empty cell array for unique elements
  uniqueCellArray = {};ID=[];
  % Iterate through the input cell array
  for i = 1:numel(cellArray)
      isUnique = ~any(cellfun(@(x)numel(union(x,cellArray{i}))==numel(cellArray{i}),uniqueCellArray));
      % Append if unique
      if isUnique
        uniqueCellArray{end + 1,1} = cellArray{i};
        ID(end+1,1) = i;
      end
  end
end

