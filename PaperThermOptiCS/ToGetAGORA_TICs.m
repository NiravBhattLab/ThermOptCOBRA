% This code enumerates on all the AGORA2 models and gets the TICs and
% their respective directions.
clear
% initCobraToolbox(0);
p = dir('./AGORA2/'); % path to agora models
agoraModels = {p(3:end).name}';
p2 = dir('./AGORA_TIC_Results/');
p2 = {p2(3:end).name}';
agoraModels = setdiff(agoraModels,p2);

for i =1:numel(agoraModels)
    load(['./AGORA2/',agoraModels{i}])
    
    [TICs,Direction,TIC_Rxns,modModel,opt] = ThermOptEnumMILP(model,1800);
    if opt
        save(['./AGORA_TIC_Results/',agoraModels{i}],'TICs','Direction','TIC_Rxns')
    end
    clear model TICs Direction TIC_Rxns modModel opt
end