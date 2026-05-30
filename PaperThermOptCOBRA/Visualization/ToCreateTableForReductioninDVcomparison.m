clear
load('../Results_TestingMatlabOptfill.mat')
tbl = table('VariableNames',{'tp','count'},'VariableTypes',{'string','double'},'Size',[0,2]);
for i=1:4
    eval(['temp=nDVs_',num2str(i),';'])
    for j=1:numel(temp)
        tbl =[tbl;{['c',num2str(i)],temp{j}(1)}];
        tbl =[tbl;{['b',num2str(i)],temp{j}(2)}];
    end
end
writetable(tbl,'lp_milp2.csv')