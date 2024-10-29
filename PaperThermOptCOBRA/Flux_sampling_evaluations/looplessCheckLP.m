function isFeasible = looplessCheckLP(points,sample)
    Nint = sample.loopMatrix;
    K = 1000;
    isFeasible = false(1,size(points,2));
    for i=1:numel(isFeasible)
        vint = points(:,i);
        vint = vint(sample.intRxns);
        % Set up LP problem
        LPproblem.A=Nint';
        LPproblem.b=zeros(size(Nint,2),1);
        LPproblem.lb=-K*(1+sign(vint))-sign(vint);
        LPproblem.ub=K*(1-sign(vint))-sign(vint);
        LPproblem.c=zeros(numel(vint),1);
        LPproblem.osense=1;%minimise
        LPproblem.csense = repmat('E',size(Nint,2),1);
        solution = solveCobraLP(LPproblem);
        if solution.stat==1
            isFeasible(i)=true;
        end
    end
end

