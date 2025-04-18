function [TICs,Direction] = OptFillTFP(model,database,mTFP)
% Enumerates all the Thermodynamically infeasible cycles in a given model
% using OptFill's algorithm. This code is adapted based on codes written in 
% GAMS and python in the paper "OptFill: A Tool for Infeasible Cycle-Free 
% Gapfilling of Stoichiometric Metabolic Models"
%
% USAGE: 
%   [TICs,Direction] = OptFillTFP(model,database,mTFP)
%
% INPUTS:
%     model:     COBRA model structure defining the actual model
%     database:  COBRA model structure defining the database
%     mTFP:      boolean value indicating whether to work on modified-TIC
%                finding problem (1) or TIC finding problem (0). In TIC
%                finding problem, atleast one reaction in a TIC will be
%                from the database.
%
% OUTPUTS:
%     TICs:       List of all the Thermodynamically infeasible cycles in
%                 the given input model-database pair
%     Direction:  Relative flux coefficients for reactions 
%                 in the corresponding TICs
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras
%       - Based on: Schroeder, W. L., & Saha, R. (2020). OptFill: a tool 
%                   for infeasible cycle-free gapfilling of stoichiometric 
%                   metabolic models. IScience, 23(1).

m_Db = mergeTwoModels(model,database); % the combined model
m_Db.origin = ~ismember(m_Db.rxns,database.rxns); % boolean vector indicating if a reaction is from model (1) or the database (0)

% parameters used were taken from the OptFill paper
epsilon=0.001;
M=1000;
tolerance=1e-6; 

changeCobraSolverParams('LP', 'feasTol', tolerance);

[~,max_phi] = size(m_Db.S); % number of reactions in the model-database pair
IrR = m_Db.ub<=0; % reactions irreversible in reverse direction
temp = m_Db.ub(IrR);
m_Db.S(:,IrR) = -m_Db.S(:,IrR);
m_Db.ub(IrR) = -1*m_Db.lb(IrR);
m_Db.lb(IrR) = -1*temp;
rev = ones(max_phi,1);
rev(m_Db.lb>=0) = 0;
m_Db.rev=rev;
% converting all the positive lower bounds to zero lower bounds
m_Db.lb(m_Db.lb>0)=0;

% normalising the bounds to max of M
temp1 = max(abs(m_Db.lb));
temp2 = max(abs(m_Db.ub));
m_Db.lb(abs(m_Db.lb)==temp1) = sign(m_Db.lb(abs(m_Db.lb)==temp1))*M;
m_Db.ub(abs(m_Db.ub)==temp2) = sign(m_Db.ub(abs(m_Db.ub)==temp2))*M;

ind = findExcRxns(m_Db); % indices of the exchange reactions

% blocking all the exchange reactions
m_Db.lb(ind) = 0;
m_Db.ub(ind) = 0;

TICs={}; Direction = {}; % empty cell to store the TICs and their respective direction
t = 1;
for phi = 2:max_phi % minimum size of a TIC has to be two
    found_all = false;
    while ~found_all
        [flux,stat] = findTIC(m_Db,phi,mTFP,epsilon,TICs,Direction);
        if stat==1
            TICs{t,1} = find(abs(flux)>epsilon*0.1);
            temp = flux/min(abs(flux));
            Direction{t,1} = temp(TICs{t,1});
            t=t+1;
        else
            found_all=1;
        end        
    end
end
end

function [flux,stat] = findTIC(model,phi,mTFP,epsilon,TICs,Direction)
    [m,n] = size(model.S);
    nTICs = numel(TICs);
    % objective function (the decision variables are in the order: v, 
    % eta, alpha, beta, gamma)
    f = [zeros(n,1);ones(n,1);zeros(n,1);zeros(n,1);zeros(nTICs,1)];
    
    % For equation (2) there are two inequality constraints
    Aineq1 = [-1*eye(n,n),epsilon*eye(n,n),zeros(n,n),diag(model.lb),zeros(n,nTICs)];
    bineq1 = zeros(n,1);
    csenseineq1 = repmat('L',n,1);
    
    Aineq2 = [eye(n,n),zeros(n,n),zeros(n,n),diag(epsilon+model.ub),zeros(n,nTICs)];
    bineq2 = model.ub;
    csenseineq2 = repmat('L',n,1);
    
    % For equation (3) there are two inequality constraints
    Aineq3 = [-1*eye(n,n),diag(model.lb),zeros(n,n),zeros(n,n),zeros(n,nTICs)];
    bineq3 = zeros(n,1);
    csenseineq3 = repmat('L',n,1);
    
    Aineq4 = [eye(n,n),-1*diag(model.ub),zeros(n,n),zeros(n,n),zeros(n,nTICs)];
    bineq4 = zeros(n,1);
    csenseineq4 = repmat('L',n,1);
    
    % For equation (4) there is one equality constraint
    Aeq1 = [model.S,zeros(m,n),zeros(m,n),zeros(m,n),zeros(m,nTICs)];
    beq1 = zeros(m,1);
    csenseeq1 = repmat('E',m,1);
    
    % For equation (5) there is one equality constraint
    Aeq2 = [zeros(1,n),ones(1,n),zeros(1,n),zeros(1,n),zeros(1,nTICs)];
    beq2 = phi;
    csenseeq2 = repmat('E',1,1);
    
    % For equation (6) there is one inequality constraint if mTFP==0
    if mTFP==0
        Aineq5 = [zeros(1,n),~model.origin',zeros(1,n),zeros(1,n),zeros(1,nTICs)];
        bineq5 = 1;
        csenseineq5 = repmat('G',1,1);
    else
        Aineq5 =[]; bineq5=[]; csenseineq5=[];
    end
    
    % For equation (7) there is one inequality constraint
    Aineq6 = [zeros(n,n),-1*eye(n,n),eye(n,n),zeros(n,n),zeros(n,nTICs)];
    bineq6 = zeros(n,1);
    csenseineq6 = repmat('L',n,1);
    
    % For equation (8) there is one inequality constraint
    Aineq7 = [zeros(n,n),zeros(n,n),eye(n,n),eye(n,n),zeros(n,nTICs)];
    bineq7 = ones(n,1);
    csenseineq7 = repmat('L',n,1);
    
    % For equation (9) there is one inequality constraint
    Aineq8 = [zeros(n,n),-1*eye(n,n),eye(n,n),eye(n,n),zeros(n,nTICs)];
    bineq8 = zeros(n,1);
    csenseineq8 = repmat('G',n,1);
    
    % For integer cut constraints (equation-10 and equation-11)
    Aineq9=[]; bineq9=[];csenseineq9=[];
    for t=1:nTICs
        cTIC = TICs{t}; cDir = Direction{t};
        pos_ids = cTIC(cDir>0);
        neg_ids = cTIC(cDir<0);
        % for postive flux reactions in the TIC
        temp1 = zeros(1,n);
        temp1(pos_ids)=1;
        temp2 = zeros(1,nTICs);
        temp2(t)=1;
        Aineq9 = [Aineq9;[zeros(1,n), zeros(1,n), temp1, zeros(1,n), temp2]];
        bineq9 = [bineq9;numel(pos_ids)];
        csenseineq9 = [csenseineq9;'L'];
        
        % for negative flux reactions in the TIC
        temp1 = zeros(1,n);
        temp1(neg_ids)=1;
        Aineq9 = [Aineq9;[zeros(1,n), zeros(1,n), zeros(1,n), temp1, -1*temp2]];
        bineq9 = [bineq9;numel(neg_ids)-1];
        csenseineq9 = [csenseineq9;'L'];
    end
    
% bounds
lb = model.lb;
lb = [lb;zeros(3*n+nTICs,1)];

ub = model.ub;
ub = [ub;ones(3*n+nTICs,1)];

% Set up LP problem
MILPproblem.A=[Aeq1;Aeq2;Aineq1;Aineq2;Aineq3;Aineq4;Aineq5;Aineq6;Aineq7;Aineq8;Aineq9];
MILPproblem.b=[beq1;beq2;bineq1;bineq2;bineq3;bineq4;bineq5;bineq6;bineq7;bineq8;bineq9];
MILPproblem.lb=lb;
MILPproblem.ub=ub;
MILPproblem.c=f;
MILPproblem.osense=1;%minimise
MILPproblem.vartype = [repmat('C',n,1);repmat('B',3*n+nTICs,1)];
MILPproblem.csense = [csenseeq1;csenseeq2;csenseineq1;csenseineq2;csenseineq3;csenseineq4;csenseineq5;...
    csenseineq6;csenseineq7;csenseineq8;csenseineq9];
solution = solveCobraMILP(MILPproblem);
stat=solution.stat;
if stat~=1
    flux = [];
else
    x=solution.full;
    flux = x(1:n);
end
end
