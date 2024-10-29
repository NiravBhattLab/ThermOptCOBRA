clear
modelNames = dir('./models');
modelNames = {modelNames(3:end).name}';
for m = 1:numel(modelNames)
    modelName = modelNames{m};
    load(['.\models\',modelName])
    model.lb(model.lb<0)=-1000;
    model.ub(:)=1000;
    tic
    [TICs,Direction,TIC_Rxns,model] = ThermOptEnumMILP(model);
    ticComp = toc;
    
    model.TICs = TICs;
    model.Direction = Direction;
    model.lb(model.lb<0)=-1000;
    model.ub(:)=1000;
    % Set up sampling parameters
    options.numSamples    = 1e4;
    options.stepsPerPoint = 1e2;
    options.numDiscarded  = 0.1*options.numSamples;
    options.algorithm     = 'ADSB';
    options.loopless      = 1;
    options.warmUpFlag    = 0;
    options.parallelFlag  =0;
    % options.numCores = 4;
    options.diagnostics   = 1;
    
    
    tic
    options.loopcheck = 'NINT';
    sampleNINT = looplessFluxSampler(model,options);
    NintT=toc;
    tic
    options.loopcheck = 'TOF';
    sampleTOF = looplessFluxSampler(model,options);
    tofT = toc;
    tic
    options.loopcheck = 'LP';
    sampleLP = looplessFluxSampler(model,options);
    lpT = toc;
    
    clearvars -except modelNames
end
