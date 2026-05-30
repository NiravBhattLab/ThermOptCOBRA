% This code is for the toy network that shows why TICs are there in a CSM
% built using ThermOptiCS
clear

% Building the toy model
R = {'R1','R2','R3','R4','B1','B2'};
Rform = {'A -> B','B -> C','B <=> F','F <=> C','-> A','C ->'};
model = createModel(R,R,Rform);

% defining the core reaction
core = find(ismember(model.rxns,{'R2','R3'}));
tol = 1e-4;

% Building the ThermOptiCS model
[TOModel,bCoreRxns,TICs,Dir,CSM_TIC,CSM_Dir] = ThermOptiCS(model,core,tol);
printRxnFormula(TOModel)