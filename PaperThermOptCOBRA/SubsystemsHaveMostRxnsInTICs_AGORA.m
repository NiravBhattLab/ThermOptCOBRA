clear
tbl = table('VariableNames',{'Model','Phylum','Subsystem','Frac_RXNS','Frac_TICs'},...
    'Size',[0,5],'VariableTypes',{'string','string','string','double','double'});
p = dir('./AGORA_TIC_Results/');
p = {p(3:end).name}';
p = strrep(p,'.mat','');
model_dir = 'D:\OneDrive_moved_files\Microbial_community_gap_filling\With_AGORA2\AGORA2\';

model2organism = readtable('./AGORAModelDetails.xlsx');

for i =1:numel(p)
    data = load(['./AGORA_TIC_Results/',p{i}]);
    if ~isfield(data,'TIC_rxns')
        data.TIC_rxns = data.TIC_Rxns;
    end
    load([model_dir,p{i}]);
    model.subSystems = cellfun(@cell2char,model.subSystems,'UniformOutput',false);
    TIC_SS = model.subSystems(ismember(model.rxns,data.TIC_rxns)); % subsystems of reactions involved in TICs
    emp = cellfun(@isempty,TIC_SS);
    TIC_SS = TIC_SS(~emp);
    if ~isempty(TIC_SS)
        [count,TIC_SS] = groupcounts(TIC_SS);
        [count2,TIC_SS2] = getTICfrac(data.TICs,model);
        for k=1:numel(TIC_SS)
            phylum = model2organism.Phylum(ismember(model2organism.MicrobeID,p{i}));
            id2 = find(ismember(TIC_SS2,TIC_SS{k}));
            new_row = {p{i},phylum,TIC_SS{k},count(k)/sum(count),count2(id2)/sum(count2)};
            tbl=[tbl;new_row];
        end
    end
    clear model data
end
writetable(tbl,'TICProp_SS_model_AGORA.csv');

function a = cell2char(a)
if ~ischar(a)
    a=a{1};
end
end
   
function [count,subsystem] = getTICfrac(TICs,model)
ss = unique(model.subSystems);
TIC_ss_mat = zeros(numel(TICs),numel(ss));
for t=1:numel(TICs)
    Css = unique(model.subSystems(ismember(model.rxns,TICs{t})));
    TIC_ss_mat(t,:) = ismember(ss,Css);
end
empty_cols = sum(TIC_ss_mat,1)==0;
ss = ss(~empty_cols);
TIC_ss_mat = TIC_ss_mat(:,~empty_cols);
count = sum(TIC_ss_mat);
subsystem=ss;
end
