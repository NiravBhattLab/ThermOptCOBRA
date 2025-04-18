function [TICs,Direction,TIC_Rxns,modModel,opt,nMILPS,nDVs] = ThermOptEnumMILP_count_milp_dvs(model,timeLimit)
% Enumerates all the Thermodynamically infeasible cycles in a given model
%
% USAGE: 
%   [TICs,Direction,TIC_Rxns,modModel,opt] = ThermoOptEnumMILP(model,timeLimit)
%
% INPUTS:
%     model:     COBRA model structure for which TICs has be found
%
% OPTIONAL INPUTS:
%     timeLimit: If the algorithm takes more than this time null set will
%                be returned for all the outputs
%
% OUTPUTS:
%     TICs:       List of all the Thermodynamically infeasible cycles in
%                 the given input model
%     Direction:  Relative flux coefficients for reactions 
%                 in the corresponding TICs
%     TIC_Rxns:   Reaction list that participates in the TICs
%     modModel:   Modified model that has no irreversible reactions that
%                 carry flux in reverse direction. The obtained TICs are
%                 for this modModel.
%     opt:        Says whether the provided solution is optimal or not
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
    a = sprintcc(model,tol); 
else
    a = fastcc(model,tol,0);
end

if isempty(a)
    TICs={};Direction={};TIC_Rxns={};opt=1;
    return
end
TIC_Rxns = model.rxns(a);
% model with only TICs
model = removeRxns(model,model.rxns(setdiff([1:numel(model.rxns)],a)));
rev = ones(numel(model.rxns),1);
rev(model.lb>=0) = 0;
model.rev=rev;
TICs={};Direction={};
i=1;
order = getOrderOfRxns(model);
RXNS = model.rxns(order);

if exist('timeLimit', 'var') && ~isempty(timeLimit)
     tic
else
    opt=1;
end
nMILPS=0;  % to count the number of MILPs used in the TOE
nDVs={}; % to count the number of decision vairables used in each of the MILP
while ~isempty(model.rxns)
    r = RXNS(1);
    core = ismember(model.rxns,r);
    nTICcons=0; % number of TIC constraints
    TICcons={}; % TIC constraints
    if exist('timeLimit', 'var') && ~isempty(timeLimit)
        if toc>timeLimit
            opt=0;
            TICs=[];Direction=[];TIC_Rxns=[];modModel=[];
            warning('TimeLimitReached')
            return 
        else
            opt=1;
        end
    end
    
    while 1
        if exist('timeLimit', 'var') && ~isempty(timeLimit)
            if toc>timeLimit
                opt=0;
                TICs=[];Direction=[];TIC_Rxns=[];modModel=[];
                warning('TimeLimitReached')
                return  
            else
                opt=1;
            end
        end
        
        dir=1;
        [allm1,blkdCore,allflux,~,nMILPS,nDVs] = getTICModel(model,core,tol,dir,TICcons,nMILPS,nDVs);
        if blkdCore
            break
        end
        for q =1:numel(allm1)
            m1 = allm1{q,1}; flux = allflux{q,1};
            TICs{i,1} = m1; 
            ids =find(ismember(model.rxns,m1));
            Direction{i,1} =flux/min(abs(flux));
            i=i+1;
            nTICcons =nTICcons+1;
            TICcons{nTICcons,1} = ids;
        end
    end
    [nTICcons,TICcons,TICs,Direction] = getReverseTIC(TICcons,TICs,Direction,model);
    i = numel(TICs)+1;
    while 1
        if exist('timeLimit', 'var') && ~isempty(timeLimit)
            if toc>timeLimit
                opt=0;
                TICs=[];Direction=[];TIC_Rxns=[];modModel=[];
                warning('TimeLimitReached')
                return 
            else
                opt=1;
            end
        end
        dir=-1;
        
        [allm1,blkdCore,allflux,~,nMILPS,nDVs] = getTICModel(model,core,tol,dir,TICcons,nMILPS,nDVs);
        if blkdCore
            break
        end

        for q =1:numel(allm1)
            m1 = allm1{q,1}; flux = allflux{q,1};
            TICs{i,1} = m1; 
            ids =find(ismember(model.rxns,m1));
            Direction{i,1} =flux/min(abs(flux));
            i=i+1;
            nTICcons =nTICcons+1;
            TICcons{nTICcons,1} = ids;
        end
    end
    model.lb(core) = 0;model.ub(core) = 0;
    if exist('sprintcc','file')
        a = sprintcc(model,tol); 
    else
        a = fastcc(model,tol,0);
    end

    rmRXNS = model.rxns(setdiff(1:numel(model.rxns),a));
    model = removeRxns(model,rmRXNS);
    rev = ones(numel(model.rxns),1);
    rev(model.lb>=0) = 0;
    model.rev=rev;
    RXNS = RXNS(~ismember(RXNS,rmRXNS));
end
end


function [nTICcons,TICcons,TICs,Direction] = getReverseTIC(pTICcons,TICs,Direction,model)
cache = cellfun(@(x)strjoin(sort(x)),TICs,'UniformOutput',false);
TICcons={};nTICcons=0;
k = numel(TICs)+1;
for t =1:numel(pTICcons)
    if sum(model.rev(pTICcons{t,1}))==numel(pTICcons{t,1})
        % all the reactions are reversible hence the reverse TIC has to be
        % there
        id = find(ismember(cache,strjoin(sort(model.rxns(pTICcons{t,1})))));
        TICs{k,1} = TICs{id,1};
        Direction{k,1} = -1*Direction{id,1};
        k=k+1;
        nTICcons = nTICcons+1;
        TICcons{nTICcons,1} = pTICcons{t,1};
    end   
end
end

function [order,bins,binsizes] = getOrderOfRxns(model)
% To get the order of reactions based on the connected components obtained
%
% USAGE: 
%   [order,bins,binsizes] = getOrderOfRxns(model)
%
% INPUTS:
%     model:     COBRA model structure for which the order of reaction has to be obtained
%
% OUTPUTS:
%     order:      Reaction IDs
%     bins:       Bins defining the IDs of connected components in which the reactions belong
%     binsizes:   Size of the connected components obtained
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras


S = double(model.S~=0); S=logical(S'*S);
g = graph(S);
try
    [bins,binsizes]=conncomp(g);
catch
    bins = conncomp(g);
    binsizes = 1:max(bins);
    binsizes = arrayfun(@(x)sum(bins==x),binsizes);
end
[binsizes,i]=sort(binsizes);
bins = arrayfun(@(x)find(i==x),bins);
order=[];
for i =1:numel(binsizes)
    order=[order;find(bins==i)'];
end
end

function [allM1,blkdCore,allFlux,stat,nMILPS,nDVs] = getTICModel(model,core,tol,dir,TICcons,nMILPS,nDVs)
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
    nMILPS=nMILPS+1;
    nDVs = [nDVs;[numel(model.rxns),numel(model.rxns)-1]];
    [reacInd,~,stat] = findConsistentReacIDTIC_feas(model,direction,tol,TICcons);
elseif dir==-1
    % for reverse direction
    direction(core)=-1;
    nMILPS=nMILPS+1;
    nDVs = [nDVs;[numel(model.rxns),numel(model.rxns)-1]];
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
    nMILPS=nMILPS+1;
    nDVs = [nDVs;[numel(model2.rxns),numel(model2.rxns)-1]];
    [reacInd2,x,stat] = findConsistentReacIDTIC(model2,direction2,tol,TICcons2);
    
    allM1{1,1} = model2.rxns(reacInd2); allFlux{1,1} = x(reacInd2);
    TICcons2{end+1,1} = find(reacInd2);
    [allM1_,blkdCore_,allFlux_,~,nMILPS,nDVs] = getTICModel(model2,ismember(model2.rxns,model.rxns(core)),tol,dir,TICcons2,nMILPS,nDVs);
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