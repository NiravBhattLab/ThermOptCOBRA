clear
Dataset = 'Kim_et_al';

mkdir(['./',Dataset,'/GeneEss'])
preProcess = {'FCC','TCC'};
Thrs = {'GL','LG','LT','SD'};
MeMs = {'Fastcore','GIMME','MBA','mCADRE','Swiftcore','TOCS'};
for p =1:numel(preProcess)
    p=preProcess{p};
    mkdir(['./',Dataset,'/GeneEss/',p])
    for t=1:numel(Thrs)
        t=Thrs{t};
        mkdir(['./',Dataset,'/GeneEss/',p,'/',t])
        for mem =1:numel(MeMs)
            mem=MeMs{mem};
            if strcmp(mem,'TOCS') && strcmp(p,'FCC')
                continue
            end
            mkdir(['./',Dataset,'/GeneEss/',p,'/',t,'/',mem])
            for i =1:64
                i=num2str(i);
                load(['./',Dataset,'/models/',p,'/',t,'/',mem,'/m',i,'.mat'])
                % setting the bounds for the biomass reaction
                m.lb(ismember(m.rxns,'BIOMASS_Ec_iJO1366_core_53p95M'))=0;
                % setting the bounds for glucose exchange
                m.lb(ismember(m.rxns,'EX_glc__D_e'))=-10;
                genes=m.genes(sum(m.rxnGeneMat,1)~=0);
	            [grRatio,grRateKO] = singleGeneDeletion(m,'FBA',genes);
                save(['./',Dataset,'/GeneEss/',p,'/',t,'/',mem,'/m',i,'.mat'],'grRatio','grRateKO')
                clear m grRatio grRateKO
            end
        end
    end
end