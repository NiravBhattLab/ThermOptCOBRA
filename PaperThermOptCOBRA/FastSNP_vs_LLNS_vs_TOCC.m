clear
changeCobraSolverParams('MILP', 'feasTol', 1e-9);
changeCobraSolverParams('LP', 'feasTol', 1e-9);
changeCobraSolverParams('MILP', 'optTol', 1e-9);
changeCobraSolverParams('LP', 'optTol', 1e-9);
p = dir('./BiggBlkdResults/');
p = {p(3:end).name}';
p = p(~startsWith(p,'Incomp'));
p2 = dir('./FastSNP_vs_LLNS_vs_TOCC/');
p = setdiff(p,{p2(3:end).name}');
ts_toe=[];ts_tocc=[];ts_fastsnp=[];ts_llns=[];
tol=1e-5;
p(ismember(p,'iJN1463.mat'))=[];
p{end+1} = 'iJN1463.mat';
for i=1:numel(p)
    i
    load(['/home/ibsc/pavan/TIC_experiments/github_upload/models/',p{i}])
    model.c(:)=0;
    tic
    [TICs,Dir,TIC_Rxns] = ThermOptEnumMILP(model);
    t_toe = toc;
    ts_toe(i)=t_toe;
    
    tic
    [Tcc] = ThermOptCC(model,tol,TICs,Dir);
    t_tocc = toc;
    ts_tocc(i)=t_tocc;
    
    tic
    [mi_fsnp,ma_fsnp]= fluxVariability(model,'allowLoops','fastSNP');
    t_fastsnp = toc;
    ts_fastsnp(i) = t_fastsnp;
    
    tic
    [mi_llns,ma_llns]= fluxVariability(model,'allowLoops','LLC-NS');
    t_llns = toc;
    ts_llns(i) = t_llns;
    
    save(['./FastSNP_vs_LLNS_vs_TOCC/',p{i}])
end
