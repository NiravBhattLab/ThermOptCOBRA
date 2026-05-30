% This code generates the metabolic models used in OptFill paper in .mat format
clear
path = ; % path to the .txt file of the respective model (can be obtained from github page of OptFill paper)
model = importdata(path);
rxns={};rxnForm={};lb=[];ub=[];
for i =1:numel(model)
    temp = split(model{i});
    temp1 = temp{1};
    temp2 = strjoin(temp(2:end),' ');
    if contains(model{i},'<->')
        rxns{i,1} =temp1;
        rxnForm{i,1} = strrep(temp2,'<->','<=>');
        lb(i) = -1000; ub(i)=1000;
    elseif contains(model{i},'->')
        rxns{i,1} =temp1;
        rxnForm{i,1} = temp2;
        lb(i) = 0; ub(i)=1000;
    elseif contains(model{i},'<-')
        rxns{i,1} =temp1;
        rxnForm{i,1} = strrep(temp2,'<-','->');
        lb(i) = -1000; ub(i)=0;
    end
        
end
model = createModel(rxns,rxns,rxnForm,...
    'lowerBoundList',lb,'upperBoundList',ub);