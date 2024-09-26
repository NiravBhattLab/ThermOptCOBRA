% this code compares the efficiency of the TOCC and LLNS algorithms in terms of time
% this code requires all the Bigg Database models to be in the folder BiggModelsMat
clear
changeCobraSolverParams('MILP', 'feasTol', 1e-9);
changeCobraSolverParams('LP', 'feasTol', 1e-9);
changeCobraSolverParams('MILP', 'optTol', 1e-9);
changeCobraSolverParams('LP', 'optTol', 1e-9);
p = dir('./BiggBlkdResults/');
p = {p(3:end).name}';
p = p(~startsWith(p,'Incomp'));
p2 = dir('./LLNS_vs_TOCC/');
p = setdiff(p,{p2(3:end).name}');
t_toe=[];t_tocc=[];t_llns=[];
tol=1e-5;
for i=1:numel(p)
    load(['./BiggModelsMat/',p{i}])
    model.c(:)=0;
    tic
    [TICs,Dir,TIC_Rxns] = ThermOptEnumMILP(model);
    temp = toc;
    t_toe(i)=temp;
    
    tic
    [Tcc] = ThermOptCC(model,tol,TICs,Dir);
    temp = toc;
    t_tocc(i)=temp;
    
    tic
    [mi,ma]= fluxVariability(model,'allowLoops','LLC-NS');
    temp = toc;
    t_llns(i) = temp;
    
    save(['./LLNS_vs_TOCC/',p{i}])
end
