clear
tbl = table('VariableNames',{'Model','Organism','Subsystem','Frac_RXNS','Frac_TICs'},...
    'Size',[0,5],'VariableTypes',{'string','string','string','double','double'});
p = dir('./BiggBlkdResults_tol_5/');
p = {p(3:end).name}';
p = strrep(p,'.mat','');
model_dir = '../Bigg models/';
model2organism = readtable('./BiggModelDetails.xlsx');
subsystem_tbl = readtable('./SubsystemStringMap.xlsx');
for i =1:numel(p)
    data = load(['./BiggBlkdResults_tol_5/',p{i}]);
    load([model_dir,p{i}]);
    model.subSystems = getCorrSubsystem(subsystem_tbl,model.subSystems);
    TIC_SS = model.subSystems(ismember(model.rxns,data.TIC_rxns)); % subsystems of reactions involved in TICs
    emp = cellfun(@isempty,TIC_SS);
    TIC_SS = TIC_SS(~emp);
    if ~isempty(TIC_SS)
        [count,TIC_SS] = groupcounts(TIC_SS);
        [count2,TIC_SS2] = getTICfrac(data.TICs,model,subsystem_tbl);
        for k=1:numel(TIC_SS)
            organism = model2organism.Organism(ismember(model2organism.BiGGID,p{i}));
            id2 = find(ismember(TIC_SS2,TIC_SS{k}));
            new_row = {p{i},organism,TIC_SS{k},count(k)/sum(count),count2(id2)/sum(count2)};
            tbl=[tbl;new_row];
        end
    end
    clear model data
end
writetable(tbl,'TICProp_SS_model_Bigg.csv');
function subsystem = getCorrSubsystem(tbl,subsystem)

for i=1:numel(subsystem)
   
    s = subsystem{i};
    if ~isempty(s)
        for j=1:width(tbl)
            [ind, ~] = find(strcmp(tbl.(tbl.Properties.VariableNames{j}), s));
            if ~isempty(ind)
                subsystem{i} = tbl{ind,1};
                break
            end
        end
        if ~strcmp(class(subsystem{i}),'char')
            subsystem{i}=subsystem{i}{1};
        end
    end
end

end
   
function [count,subsystem] = getTICfrac(TICs,model,tbl)
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
