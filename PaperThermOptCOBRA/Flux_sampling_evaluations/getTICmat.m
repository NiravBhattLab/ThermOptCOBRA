function sample = getTICmat(sample)
rxns = sample.rxns;
TICs = sample.TICs;
dir = sample.Direction;
TICmat = sparse(numel(rxns),numel(TICs));
for t = 1:numel(TICs)
    if all(ismember(TICs{t},rxns))
        [loca,locb]=ismember(rxns,TICs{t});
        d=sign(dir{t});
        TICmat(loca,t) = d(locb(locb~=0));
    end
end
emptyTICIDs = sum(abs(TICmat),1)==0;
TICmat(:,emptyTICIDs) =[];
sample.TICs(emptyTICIDs)=[];
sample.Direction(emptyTICIDs)=[];
sample.TICmat=TICmat;
end