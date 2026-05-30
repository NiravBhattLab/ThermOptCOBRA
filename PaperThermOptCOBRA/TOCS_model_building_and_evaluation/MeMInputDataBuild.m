% This code builds the inputs required for different MeMs as preprocessed
% by localT2, Global and StanDep
clear
% initCobraToolbox(0);
load('iJO1366.mat')

% checking for reactions irreversible in reverse direction
IrR = model.ub<=0; % reactions irreversible in reverse direction
temp = model.ub(IrR);
model.S(:,IrR) = -model.S(:,IrR);
model.ub(IrR) = -1*model.lb(IrR);
model.lb(IrR) = -1*temp;

% changing all the negative lower bounds to its defaults
model.lb(model.lb<0)=-1000;
sol = optimizeCbModel(model);
model.lb(find(model.c))=0.8*sol.f;
DataSet = 'Kim_et_al';
mkdir(['./',DataSet,'/inputs'])
mkdir(['./',DataSet,'/inputs/FCC'])
mkdir(['./',DataSet,'/inputs/TCC'])

% loading the datset
tbl = readtable('64Samples.xlsx');
geneExp = struct();
geneExp.gene=tbl.Var1;
geneExp.valuebyTissue=table2array(tbl(:,2:end));
geneExp.valuebyTissue = 2.^(geneExp.valuebyTissue);
temp = max(geneExp.valuebyTissue(abs(geneExp.valuebyTissue)~=inf),[],'all');
geneExp.valuebyTissue(abs(geneExp.valuebyTissue)==inf) = temp*sign(geneExp.valuebyTissue(abs(geneExp.valuebyTissue)==inf));
geneExp.Tissue = arrayfun(@num2str,[1:size(geneExp.valuebyTissue,2)]','UniformOutput',false);
%% Using the consistent model from fastcc

% Removing all the blocked reactions
a = fastcc(model,1e-5);
model_fcc = removeRxns(model,model.rxns(setdiff(1:numel(model.rxns),a)));
save('consIJO1366.mat','model_fcc')
PreProcessExpressionData(model_fcc,'FCC',geneExp,DataSet)

%% Using the consistent model from ThermOptCC

% Removing all the T-blocked reactions
[a,TICs,Dir,model,TIC_Rxns] = ThermOptCC(model,1e-5);
% removing thermodynamically blocked reactions
model_tcc = removeRxns(model,model.rxns(ismember(a,'Blocked')));
save('TconsIJO1366.mat','model_tcc')
% load('TconsIJO1366.mat')
PreProcessExpressionData(model_tcc,'TCC',geneExp,DataSet)

function PreProcessExpressionData(model,preprocess_type,geneExp,DataSet)
    % for the case of standep
    modelData = getModelData(geneExp,model);
    spec = getSpecialistEnzymes(model);  
    prom = getPromEnzymes(model);
    enzymeData = comparePromiscuousSpecific(spec,prom,modelData);
    edgeX = [-3 -2 -1 0 1 2 3 4 5 6]; % bins  
    distMethod = 'euclidean'; % distance method  
    linkageMethod = 'complete'; % linkage metric for hierarchical clustering
    clustObj = geneExprDist_hierarchy(enzymeData,[],edgeX,7,distMethod,linkageMethod);
    core_rxn = models4mClusters1(clustObj,enzymeData.Tissue,model,edgeX,[],[],false,0,[1 1]);
    % fastcore and swiftcore
    save(['./',DataSet,'/inputs/',preprocess_type,'/SDinFASTCORE'],'core_rxn')
    % GIMME
    core_rxn = getGIMMEscores(clustObj,edgeX,model);
    core_rxn(core_rxn==-inf)=min(min(core_rxn(core_rxn~=-inf)));
    save(['./',DataSet,'/inputs/',preprocess_type,'/SDinGIMME'],'core_rxn')
    % mCADRE
    [core_rxn,~] = getUbiquityScore(clustObj,edgeX,model);
    save(['./',DataSet,'/inputs/',preprocess_type,'/SDinmCADRE'],'core_rxn')
    % MBA
    [H,M] = getMBAsets(clustObj,edgeX,model,0.1);
    save(['./',DataSet,'/inputs/',preprocess_type,'/SDinMBA'],'H','M')

    %  for the case of LocalT2
    lowerThs=10;
    upperThs=80;
    thrRxnData = getLocalT2_case(modelData,model,lowerThs,upperThs);
    % fastcore and swiftcore
    core_rxn = thrRxnData.value> 5*log(2);
    save(['./',DataSet,'/inputs/',preprocess_type,'/LTinFASTCORE'],'core_rxn')
    % GIMME and mCADRE
    core_rxn = thrRxnData.value;
    save(['./',DataSet,'/inputs/',preprocess_type,'/LTinGIMME'],'core_rxn')
    save(['./',DataSet,'/inputs/',preprocess_type,'/LTinmCADRE'],'core_rxn') 
    % MBA 
    core_rxn = thrRxnData.value;
    up_thr = prctile(core_rxn(core_rxn~=-2),75,'all');
    lw_thr = 5*log(2);
    H = core_rxn>=up_thr;
    M = core_rxn>=lw_thr & core_rxn<up_thr; 
    save(['./',DataSet,'/inputs/',preprocess_type,'/LTinMBA'],'H','M')
    
    %  for the case of Global80
    thr=80;
    thrRxnData = getGlobal_case(modelData,model,thr);
    % fastcore and swiftcore
    core_rxn = thrRxnData.value> 5*log(2);
    save(['./',DataSet,'/inputs/',preprocess_type,'/GLinFASTCORE'],'core_rxn')
    % GIMME and mCADRE
    core_rxn = thrRxnData.value;
    save(['./',DataSet,'/inputs/',preprocess_type,'/GLinGIMME'],'core_rxn')
    save(['./',DataSet,'/inputs/',preprocess_type,'/GLinmCADRE'],'core_rxn') 
    % MBA 
    core_rxn = thrRxnData.value;
    up_thr = prctile(core_rxn(core_rxn~=-2),75,'all');
    lw_thr = 5*log(2);
    H = core_rxn>=up_thr;
    M = core_rxn>=lw_thr & core_rxn<up_thr; 
    save(['./',DataSet,'/inputs/',preprocess_type,'/GLinMBA'],'H','M')
end