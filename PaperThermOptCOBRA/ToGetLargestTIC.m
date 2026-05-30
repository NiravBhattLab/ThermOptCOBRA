% This code gets the details of the largest TIC present in BiGG and AGORA
% models
clear
path = ; % path to the AGORA TIC results (or) BiGG TIC results
p = dir(path);
p={p(3:end).name}';
largestTIC={''};
noTICModels={};nTICs=[];nTICRxns=[];

for z=1:numel(p)
    load([path,p{z}])
    if isempty(TICs)
        TICs={};
    end
    [~,I] = sort(cellfun(@length,TICs));
    nTICs(z) = numel(TICs);
    nTICRxns(z) = numel(TIC_Rxns);
    TICs = TICs(I);
    if isempty(TICs)
        noTICModels=[noTICModels;p{z}];
    else
        temp = TICs(end);
        if numel(temp{:})> numel(largestTIC{:})
            largestTIC=temp;
            CorrModel=p{z};
        end
    end
end