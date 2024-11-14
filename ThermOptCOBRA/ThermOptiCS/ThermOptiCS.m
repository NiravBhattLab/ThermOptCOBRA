function [Model,bCoreRxns,TICs,Dir,CSM_TIC,CSM_Dir] = ThermOptiCS(model,core,tol,TICs,Dir,doTOCC,timeLimit)
% Thermodynamically feasible context-specific model building.
%
% USAGE: 
%   [Model,bCoreRxns,TICs,Dir] = ThermOptiCS(model,core,tol)
%
% INPUTS:
%     model:     COBRA model structure for which act as the GSMM
%     core:      Reaction IDs which are defined to be core (These
%                reactions will be present in the final model)
%     tol:       Tolerance value (User defined non-zero value).
%
% OPTIONAL INPUTS:
%     TICs:      TICs of the 'model' (output from ThermOptEnumMILP)
%     Dir:       Direction of the TICs (output from ThermOptEnumMILP)
%     doTOCC:    Boolean value.
%                  1: ThermOptCC will be run on the model and blocked
%                     reactions will be removed
%                  0: ThermOptCC will not be run on the model (This
%                     requires the input model to be thermodynamically-flux
%                     consistent)
%     timeLimit: Upper time limit in seconds to solve the minReac
%                optimization problem (Default: 120s)
% OUTPUTS:
%     Model:     Context-specific model
%     bCoreRxns: Core reactions that are thermodynamically blocked (or infeasible)
%     TICs:      List of all the Thermodynamically infeasible cycles in
%                the given input model
%     Dir:       The flux directions for reactions in the corresponding
%                TICs
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras


% checking for reactions irreversible in reverse direction
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

if ~exist('doTOCC','var')||isempty(doTOCC)
    doTOCC=1;
end

if ~exist('TICs', 'var') || isempty(TICs) || ~exist('Dir', 'var') || isempty(Dir)
    [a,TICs,Dir] = ThermOptCC(model,tol);
else
    if doTOCC
        [a,TICs,Dir] = ThermOptCC(model,tol,TICs,Dir);
    else
        a = repmat({''},numel(model.rxns),1);
    end
end
if ~exist('timeLimit','var')||isempty(timeLimit)
    timeLimit=120;
end

bRxns = find(ismember(a,'Blocked')); % thermodynamically blocked reactions
bCoreRxns = intersect(bRxns,core); % getting the blocked core reactions
if ~isempty(bCoreRxns)
    warning('Few reactions in the core reaction list are thermodynamically blocked')
    fprintf('\n Thermodynamically blocked reactions will be removed from core reaction list\n')
    core =  setdiff(core,bRxns);
end

[modelIrrev,matchRev, rev2irrev,irrev2rev] = convertToIrreversible(model);
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

minNetwork = [];
core = rev2irrev(core);core = horzcat(core{:}); % core reactions in terms of irreversible model
P = setdiff([1:numel(modelIrrev.rxns)],core); % penalty set
while ~isempty(core)
    [IDScc,x0] = findConsistentIDS_TOCS(modelIrrev,core,TICmat,rev2irrev,tol);
    IDS = findMinNetwork(modelIrrev,intersect(core,IDScc),P,TICmat,rev2irrev,tol,x0,timeLimit);
    minNetwork = [minNetwork;IDS];
    % getting all the reactions corresponding to those irreversible
    % reactions
    modIDS = matchRev(IDS);
    modIDS(modIDS==0)=[]; modIDS=union(IDS,modIDS);
    core = setdiff(core,modIDS);
    P = setdiff(P,modIDS);
end
minNetwork = unique(irrev2rev(minNetwork));
Model = removeRxns(model,model.rxns(setdiff([1:numel(model.rxns)],minNetwork)));

% getting the TICs for the built context-specific model
temp = [];
for i=1:numel(TICs)
    currTIC=TICs{i};
    if sum(ismember(Model.rxns,currTIC))==numel(currTIC)
        temp = [temp;i];
    end
end
CSM_TIC = TICs(temp);
CSM_Dir = Dir(temp);

end

