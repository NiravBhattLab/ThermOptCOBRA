clear
thrs = {'GL','LG','LT','SD'};
mems = {'Fastcore','GIMME','MBA','mCADRE','Swiftcore','TOCS'};
tbl = table('VariableNames',{'MeM','Threshold','model','TIC_Count'},'Size',[numel(thrs)*numel(mems)*64,4],...
    'VariableTypes',{'string','string','string','double'});
id=1;

for t = 1:numel(thrs)
    thr = thrs{t};
    for m = 1:numel(mems)
        mem = mems{m};
        for i=1:64
            load(['./Kim_et_al/TCC_results/TCC/',thr,'/',mem,'/m',num2str(i)],'TICs');
            modelName = ['E',num2str(2198+i)];
            tbl(id,:) = {mem,thr,modelName,numel(TICs)};
            id=id+1;
            clear TICs
        end
    end
end
writetable(tbl,'TIC_count_CSMs.csv')