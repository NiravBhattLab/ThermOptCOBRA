clear
model = 'm2';
load(['..\Kim_et_al\models\TCC\LG\TOCS\',model,'.mat'])
[TICs,Direction,TIC_Rxns,m,opt] = ThermOptEnumMILP(m);
nSamples = 100000;
chrr_samples= chrrSampler(m,100,nSamples);

chrr_samples_no_tic=zeros(size(chrr_samples));
for i=1:nSamples
    f = chrr_samples(:,i);
    [flux] = ThermOptFlux(m,f,TICs,Direction);
    chrr_samples_no_tic(:,i)=flux;
    clear flux
end

tbl = table('Size',[nSamples,numel(TIC_Rxns)],'VariableNames',TIC_Rxns,'VariableTypes',repmat({'double'},1,numel(TIC_Rxns)));
tbl_chrr = tbl;
tbl_tof = tbl;
for i =1:numel(TIC_Rxns)
    id = find(ismember(m.rxns,TIC_Rxns{i}));
    tbl_chrr.(TIC_Rxns{i}) = chrr_samples(id,:)';
    tbl_tof.(TIC_Rxns{i}) = chrr_samples_no_tic(id,:)';
end
save(['TOF_samples_TIC_rxns_LG_TOCS2_',model])
writetable(tbl_tof,['TOF_samples_TIC_rxns_LG_TOCS2_',model,'.csv'])
writetable(tbl_chrr,['CHRR_samples_TIC_rxns_LG_TOCS2_',model,'.csv'])