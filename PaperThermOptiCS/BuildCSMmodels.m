% This code builds the context-specific models using FastCore and
% ThermOptiCS given the GEM and gene-expression data
clear
load('iML1515') % this model has to be present in the path
% initCobraToolbox(0);
changeCobraSolverParams('MILP', 'feasTol', 1e-9);
changeCobraSolverParams('LP', 'feasTol', 1e-9);
changeCobraSolverParams('MILP', 'optTol', 1e-9);
changeCobraSolverParams('LP', 'optTol', 1e-9);

[a,TICs,Dir,model,TIC_Rxns] = ThermOptCC(iML1515,1e-4);
% removing blocked reactions
Fid=ismember(a,'Forward');
model.lb(Fid)=max([zeros(sum(Fid),1),model.lb(Fid)],[],2);
Rid=ismember(a,'Reverse');
model.ub(Rid)=min([zeros(sum(Rid),1),model.ub(Rid)],[],2);
model = removeRxns(model,model.rxns(ismember(a,'Blocked')));

filename = ; % path to the gene-expression data file
tbl = readtable(filename);
geneExp = struct();
geneExp.genes=tbl.Var1;
geneExp.value=table2array(tbl(:,2:end));
geneExp.value = 2.^(geneExp.value);
geneExp.context = arrayfun(@num2str,[1:size(geneExp.value,2)]','UniformOutput',false);
[RxnImp] = GiniReactionImportance(geneExp,model,80,10,[find(model.c);find(contains(model.rxns,'ATPM'))]);


for i=1:size(RxnImp,2)
    c = find(RxnImp(:,i)>=1);
    m = fastcore(model,c,1e-4);
    save(['./FastCore_models/m',num2str(i)],'m')
    m = ThermOptiCS(model,c,1e-4);
    save(['./ThermOptiCS_models/m',num2str(i)],'m')
end