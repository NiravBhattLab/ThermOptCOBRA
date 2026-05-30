clear
Dataset = 'Kim_et_al';
preProcess = {'FCC','TCC'};
Thrs = {'GL','LG','LT','SD'};
MeMs = {'Fastcore','GIMME','MBA','mCADRE','Swiftcore','TOCS'};
tbl=readtable('./Ecoli_essentiality.csv');
thr=0.05;
for p =1:numel(preProcess)
    p=preProcess{p};
    for t=1:numel(Thrs)
        t=Thrs{t};
        for mem =1:numel(MeMs)
            mem=MeMs{mem};
            if strcmp(mem,'TOCS') && strcmp(p,'FCC')
                continue
            end
            for i =1:64
                i=num2str(i);
                load(['./',Dataset,'/GeneEss/',p,'/',t,'/',mem,'/m',i,'.mat'])
                load(['./',Dataset,'/models/',p,'/',t,'/',mem,'/m',i,'.mat'])
                genes=m.proteins(sum(m.rxnGeneMat,1)~=0);
                common_genes = intersect(genes,tbl.Gene);
                [~,ib] = ismember(common_genes,genes);
                grRatio_2 = grRatio(ib);
                grRateKO_2 = grRateKO(ib);
                [ia,ib] = ismember(common_genes,tbl.Gene);
                tbl_temp = tbl(ib,:);
                % true positives
                TP = sum(tbl_temp.Essential & grRatio_2<thr);
                % true negatives
                TN = sum(~tbl_temp.Essential & grRatio_2>thr);
                % false positives
                FP = sum(~tbl_temp.Essential & grRatio_2<thr);
                % false negatives
                FN = sum(tbl_temp.Essential & grRatio_2>thr);
                % calculate some metrics
                sensitivity = TP./(TP + FN);
                specificity = TN./(TN + FP);
                precision = TP./(TP+FP);
                accuracy = (TP + TN)./(TP + TN + FP + FN);
                F1 = 2*TP./(2*TP + FP + FN); %F1 coefficient
                MCC = ((TP.*TN) - (FP.*FN))./sqrt((TP+FP).*(TP+FN).*(TN+FP).*(TN+FN));  % Matthews correlation coefficient
                rocObj = rocmetrics(tbl_temp.Essential,double(grRatio_2<thr),[1]);
                roc_auc = rocObj.AUC;
                save(['./',Dataset,'/GeneEss/',p,'/',t,'/',mem,'/m',i,'.mat'],'grRatio','grRateKO','grRatio_2','grRateKO_2',...
                    'genes','common_genes','sensitivity','specificity','precision','accuracy','F1','MCC','roc_auc')
                clear m grRatio grRateKO grRatio_2 grRateKO_2 genes common_genes sensitivity specificity precision accuracy F1 MCC roc_auc
            end
        end
    end
end

