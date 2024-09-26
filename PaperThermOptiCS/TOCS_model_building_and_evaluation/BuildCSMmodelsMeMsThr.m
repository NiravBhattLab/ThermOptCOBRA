clear
% initCobraToolbox(0);
changeCobraSolverParams('MILP', 'feasTol', 1e-9);
changeCobraSolverParams('LP', 'feasTol', 1e-9);
changeCobraSolverParams('MILP', 'optTol', 1e-9);
changeCobraSolverParams('LP', 'optTol', 1e-9);

DataSet = 'Kim_et_al';
mkdir(['./',DataSet,'/models'])
manual_core = {'BIOMASS_Ec_iJO1366_core_53p95M','ATPM'};
%% For the case fastcc consistent model
load('consIJO1366')
model = model_fcc;
model=buildRxnGeneMat(model);
sol = optimizeCbModel(model);
pre_process = 'FCC';
mkdir(['./',DataSet,'/models/',pre_process])

% % for global threshold case
thr='GL';
BuildGlandLTModels(model,manual_core,DataSet,pre_process,thr)

% % for LocalT2 threshold case
thr='LT';
BuildGlandLTModels(model,manual_core,DataSet,pre_process,thr)

% % for StanDep threshold case
thr = 'SD';
BuildSDModels(model,manual_core,DataSet,pre_process,thr)

% % for LG threshold case
thr = 'LG';
BuildLGModels(model,manual_core,DataSet,pre_process,thr)

%% For the case TOCC consistent model
load('TconsIJO1366')
model = model_tcc;
model=buildRxnGeneMat(model);
sol = optimizeCbModel(model);

pre_process = 'TCC';
mkdir(['./',DataSet,'/models/',pre_process])

% % for global threshold case
thr='GL';
BuildGlandLTModels(model,manual_core,DataSet,pre_process,thr)

% % for LocalT2 threshold case
thr='LT';
BuildGlandLTModels(model,manual_core,DataSet,pre_process,thr)

% % for StanDep threshold case
thr = 'SD';
BuildSDModels(model,manual_core,DataSet,pre_process,thr)

% % for LG threshold case
thr = 'LG';
BuildLGModels(model,manual_core,DataSet,pre_process,thr)

%% Required functions
function BuildLGModels(model,manual_core,DataSet,pre_process,thr)

    % loading the datset
    tbl = readtable('64Samples.xlsx');
    geneExp = struct();
    geneExp.genes=tbl.Var1;
    geneExp.value=table2array(tbl(:,2:end));
    geneExp.value = 2.^(geneExp.value);
    temp = max(geneExp.value(abs(geneExp.value)~=inf),[],'all');
    geneExp.value(abs(geneExp.value)==inf) = temp*sign(geneExp.value(abs(geneExp.value)==inf));
    geneExp.context = arrayfun(@num2str,[1:size(geneExp.value,2)]','UniformOutput',false);

    filename_1 = ['./',DataSet,'/models/',pre_process,'/',thr,'/Fastcore/'];
    filename_2 = ['./',DataSet,'/models/',pre_process,'/',thr,'/Swiftcore/'];
    mkdir(filename_1)
    mkdir(filename_2)

    coreRxn = find(ismember(model.rxns,manual_core));
    [~,~] = buildContextmodels(geneExp,model,'FASTCORE',geneExp.context,80,10,1,coreRxn,filename_1,1:numel(model.rxns),1e-6);
    [~,~] = buildContextmodels(geneExp,model,'SWIFTCORE',geneExp.context,80,10,1,coreRxn,filename_2,1:numel(model.rxns),1e-6);
    if strcmp(pre_process,'TCC')
        filename = ['./',DataSet,'/models/',pre_process,'/',thr,'/TOCS/'];
        mkdir(filename)
        [~,~] = buildContextmodels(geneExp,model,'TOCS',geneExp.context,80,10,1,coreRxn,filename,1:numel(model.rxns),1e-6);
    end
    filename = ['./',DataSet,'/models/',pre_process,'/',thr,'/MBA/'];
    mkdir(filename)
    [~,~] = buildContextmodels(geneExp,model,'MBA',geneExp.context,80,10,1,coreRxn,filename,1:numel(model.rxns),1e-6);

    filename = ['./',DataSet,'/models/',pre_process,'/',thr,'/GIMME/'];
    mkdir(filename)
    [~,~] = buildContextmodels(geneExp,model,'GIMME',geneExp.context,80,10,1,coreRxn,filename,1:numel(model.rxns),1e-6);

    filename = ['./',DataSet,'/models/',pre_process,'/',thr,'/mCADRE/'];
    mkdir(filename)
    [~,~] = buildContextmodels(geneExp,model,'mCADRE',geneExp.context,80,10,1,coreRxn,filename,1:numel(model.rxns),1e-6);

end
function BuildGlandLTModels(model,manual_core,DataSet,pre_process,thr)
    % fastcore and swiftcore
    load(['./',DataSet,'/inputs/',pre_process,'/',thr,'inFASTCORE'])
    filename_1 = ['./',DataSet,'/models/',pre_process,'/',thr,'/Fastcore/'];
    filename_2 = ['./',DataSet,'/models/',pre_process,'/',thr,'/Swiftcore/'];
    mkdir(filename_1)
    mkdir(filename_2)
    if strcmp(pre_process,'TCC')
        filename_3 = ['./',DataSet,'/models/',pre_process,'/',thr,'/TOCS/'];
        mkdir(filename_3)
    end
    core_rxn(ismember(model.rxns,manual_core),:)=1;
    [TICs,Dir] = ThermOptEnumMILP(model);
    for i=1:size(core_rxn,2)
        core = find(core_rxn(:,i));

        m = fastcore(model,core,1e-6);
        fn = [filename_1,'m',num2str(i),'.mat'];
        save(fn,'m')

        m = swiftcore(model, core, ones(numel(model.rxns),1), 1e-6, 0);
        fn = [filename_2,'m',num2str(i),'.mat'];
        save(fn,'m')

        if strcmp(pre_process,'TCC')
            m = ThermOptiCS(model,core,1e-6,TICs,Dir,0,300);
            fn = [filename_3,'m',num2str(i),'.mat'];
            save(fn,'m')
        end
    end

    % GIMME
    load(['./',DataSet,'/inputs/',pre_process,'/',thr,'inGIMME'])
    filename = ['./',DataSet,'/models/',pre_process,'/',thr,'/GIMME/'];
    mkdir(filename)
    core_rxn(ismember(model.rxns,manual_core),:)=10*log(2);
    core_rxn(find(sum(model.rxnGeneMat,2)==0),:)=5*log(2);
    for i=1:size(core_rxn,2)
        core = core_rxn(:,i);
        m = GIMME(model,core,5*log(2));
        fn = [filename,'m',num2str(i),'.mat'];
        save(fn,'m')
    end

    % MBA
    load(['./',DataSet,'/inputs/',pre_process,'/',thr,'inMBA'])
    filename = ['./',DataSet,'/models/',pre_process,'/',thr,'/MBA/'];
    mkdir(filename)
    H(ismember(model.rxns,manual_core),:)=1;
    for i=1:size(core_rxn,2)
        curr_H = model.rxns(find(H(:,i)));
        curr_M = model.rxns(find(M(:,i)));
        m = MBA(model,curr_M,curr_H,1e-6);
        fn = [filename,'m',num2str(i),'.mat'];
        save(fn,'m')
    end
    
    % mCADRE
    load(['./',DataSet,'/inputs/',pre_process,'/',thr,'inmCADRE'])
    filename = ['./',DataSet,'/models/',pre_process,'/',thr,'/mCADRE/'];
    mkdir(filename)
    core_rxn(find(sum(model.rxnGeneMat,2)==0),:)=-1;
    core_rxn(ismember(model.rxns,manual_core),:)=10*log(2);
     
    for i=1:size(core_rxn,2)
        ubiq = core_rxn(:,i);
        m = mCADRE(model,ubiq,zeros(numel(model.rxns),1),{},0,1/3,1e-6);
        fn = [filename,'m',num2str(i),'.mat'];
        save(fn,'m')
    end
end

function BuildSDModels(model,manual_core,DataSet,pre_process,thr)
    % fastcore and swiftcore
    load(['./',DataSet,'/inputs/',pre_process,'/',thr,'inFASTCORE'])
    filename_1 = ['./',DataSet,'/models/',pre_process,'/',thr,'/Fastcore/'];
    filename_2 = ['./',DataSet,'/models/',pre_process,'/',thr,'/Swiftcore/'];
    mkdir(filename_1)
    mkdir(filename_2)
    if strcmp(pre_process,'TCC')
        filename_3 = ['./',DataSet,'/models/',pre_process,'/',thr,'/TOCS/'];
        mkdir(filename_3)
    end
    core_rxn(ismember(model.rxns,manual_core),:)=1;
    [TICs,Dir] = ThermOptEnumMILP(model);
    for i=1:size(core_rxn,2)
        core = find(core_rxn(:,i));

        m = fastcore(model,core,1e-6);
        fn = [filename_1,'m',num2str(i),'.mat'];
        save(fn,'m')

        m = swiftcore(model, core, ones(numel(model.rxns),1), 1e-6, 0);
        fn = [filename_2,'m',num2str(i),'.mat'];
        save(fn,'m')

        if strcmp(pre_process,'TCC')
            m = ThermOptiCS(model,core,1e-6,TICs,Dir,0,300);
            fn = [filename_3,'m',num2str(i),'.mat'];
            save(fn,'m')
        end
    end
    
    % GIMME
    load(['./',DataSet,'/inputs/',pre_process,'/',thr,'inGIMME'])
    filename = ['./',DataSet,'/models/',pre_process,'/',thr,'/GIMME/'];
    mkdir(filename)
    core_rxn(ismember(model.rxns,manual_core),:)=1;
    core_rxn(find(sum(model.rxnGeneMat,2)==0),:)=0;
    for i=1:size(core_rxn,2)
        core = core_rxn(:,i);
        m = GIMME(model,core,0);
        fn = [filename,'m',num2str(i),'.mat'];
        save(fn,'m')
    end
    
    % MBA
    load(['./',DataSet,'/inputs/',pre_process,'/',thr,'inMBA'])
    filename = ['./',DataSet,'/models/',pre_process,'/',thr,'/MBA/'];
    mkdir(filename)
    H(ismember(model.rxns,manual_core),:)=1;
    for i=1:size(core_rxn,2)
        curr_H = model.rxns(find(H(:,i)));
        curr_M = model.rxns(find(M(:,i)));
        m = MBA(model,curr_M,curr_H,1e-6);
        fn = [filename,'m',num2str(i),'.mat'];
        save(fn,'m')
    end
    
    % mCADRE
    load(['./',DataSet,'/inputs/',pre_process,'/',thr,'inmCADRE'])
    filename = ['./',DataSet,'/models/',pre_process,'/',thr,'/mCADRE/'];
    mkdir(filename)
    core_rxn(find(sum(model.rxnGeneMat,2)==0),:)=-1; 
    core_rxn(ismember(model.rxns,manual_core),:)=1;
    
    for i=1:size(core_rxn,2)
        ubiq = core_rxn(:,i);
        m = mCADRE(model,ubiq,zeros(numel(model.rxns),1),{},0,1/3,1e-6);
        fn = [filename,'m',num2str(i),'.mat'];
        save(fn,'m')
    end
end