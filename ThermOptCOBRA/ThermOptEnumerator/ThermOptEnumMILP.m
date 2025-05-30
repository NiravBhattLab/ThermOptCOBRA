function [TICs,Direction,TIC_Rxns,modModel,opt] = ThermOptEnumMILP(model,timeLimit)
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
        [allm1,blkdCore,allflux] = getTICModel(model,core,tol,dir,TICcons);
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
        
        [allm1,blkdCore,allflux] = getTICModel(model,core,tol,dir,TICcons);
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