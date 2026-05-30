clear
load('./FastSNP_vs_LLNS_vs_TOCC/iJN1463.mat');
tbl = table('VariableNames',{'model','TOE','TOCC','TOE_TOCC','Fast-SNP','LLL-FVA','Fastest_method'},...
    'VariableTypes',{'string','double','double','double','double','double','string'},...
    'size',[0,7]);
for i=1:numel(p)
    if ts_toe(i)+ts_tocc(i)<ts_fastsnp(i) && ts_toe(i)+ts_tocc(i)<ts_llns(i)
        temp2 = 'TOE_TOCC';
    elseif ts_fastsnp(i)<ts_toe(i)+ts_tocc(i) && ts_fastsnp(i)<ts_llns(i)
        temp2 = 'Fast-SNP';
    elseif ts_llns(i)<ts_toe(i)+ts_tocc(i) && ts_llns(i)<ts_fastsnp(i)
        temp2 = 'LLL-FVA';
    end
    temp = {p{i},ts_toe(i),ts_tocc(i),ts_toe(i)+ts_tocc(i),ts_fastsnp(i),...
        ts_llns(i),temp2};
    tbl = [tbl;temp];
end
