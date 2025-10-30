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

core = 1:n; 
runtime=60;
while 1
    [IDS,stat] = findConsistentIDS(modelIrrev,core,TICmat,rev2irrev,tol,runtime);
    if isempty(IDS)
        if stat==3
            runtime = runtime+60;
            continue
        else
            break
        end
    else
        core = setdiff(core,IDS);        
    end
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


