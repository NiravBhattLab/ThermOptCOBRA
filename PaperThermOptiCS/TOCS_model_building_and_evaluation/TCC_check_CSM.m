clear
Dataset = 'Kim_et_al';
changeCobraSolverParams('MILP', 'feasTol', 1e-9);
changeCobraSolverParams('LP', 'feasTol', 1e-9);
changeCobraSolverParams('MILP', 'optTol', 1e-9);
changeCobraSolverParams('LP', 'optTol', 1e-9);
tol=1e-5;
mkdir(['./',Dataset,'/TCC_results'])
preProcess = {'FCC','TCC'};
Thrs = {'GL','LG','LT','SD'};
MeMs = {'Fastcore','GIMME','MBA','mCADRE','Swiftcore','TOCS'};
for p =1:numel(preProcess)
    p=preProcess{p};
    mkdir(['./',Dataset,'/TCC_results/',p])
    for t=1:numel(Thrs)
        t=Thrs{t};
        mkdir(['./',Dataset,'/TCC_results/',p,'/',t])
        for mem =1:numel(MeMs)
            mem=MeMs{mem};
            if strcmp(mem,'TOCS') && strcmp(p,'FCC')
                continue
            end
            mkdir(['./',Dataset,'/TCC_results/',p,'/',t,'/',mem])
            for i =1:64
                i=num2str(i);
                load(['./',Dataset,'/models/',p,'/',t,'/',mem,'/m',i,'.mat'])
                m.lb(m.lb>0)=0;
                [a,TICs,Dir,~,TIC_Rxns] = ThermOptCC(m,tol);
                save(['./',Dataset,'/TCC_results/',p,'/',t,'/',mem,'/m',i,'.mat'],'a','TICs','Dir','TIC_Rxns')
                clear a TICs Dir TIC_Rxns m
            end
        end
    end
end