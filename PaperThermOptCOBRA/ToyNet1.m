% This code is for the toy network that gives overview about the
% ThermOptEnumerator
clear
R = {'R1','R2','R3','R4','R5','R6','R7','R8','R9','R10','B1','B2'};

Rform = {'A -> B','B -> C','I -> J','J -> K','D <=> E','H -> G',...
    'H -> I','F <=> E','K -> G','F <=> D','-> A','C ->'};

lbs = [0,0,0,0,-10,-10,0,-10,0,-10,0,0];
ubs = [10,10,10,10,10,0,10,10,10,10,10,10];

model = createModel(R,R,Rform,'lowerBoundList',lbs,'upperBoundList',ubs);

[TICs,Direction,TIC_Rxns,modModel,opt] = ThermOptEnumMILP(model);
