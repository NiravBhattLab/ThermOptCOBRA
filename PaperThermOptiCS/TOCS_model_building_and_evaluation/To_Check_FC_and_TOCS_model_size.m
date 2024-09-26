clear
Thrs = {'GL','LG','LT','SD'};
for t=1:numel(Thrs)
    t = Thrs{t};
    count=0;
    count2=0;
    for i=1:64
        m1 = load(['./Kim_et_al/models/TCC/',t,'/Fastcore/m',num2str(i)]);
        f = m1.m;
        m1 = load(['./Kim_et_al/models/TCC/',t,'/TOCS/m',num2str(i)]);
        to = m1.m;
        if numel(f.rxns)<numel(to.rxns)
            count = count+1;
            load(['./Kim_et_al/TCC_results/TCC/',t,'/Fastcore/m',num2str(i)],'a')
            if any(ismember(a,'Blocked'))
                count2 = count2+1;
            end
            clear a
        end
    end
    fprintf(['For the thresholding method ',t,'\n'])
    fprintf('No of models in which FC has lesser number of reactions is %d\n',count)
    fprintf('No of minimum sized fastcore models which has atleast 1 t-blkd reactions is %d\n',count2)
end