function isFeasible = looplessCheckTOF(points,sample)
% Checks the loopless condition using the TICmat
%
% Checks loopless condition of a series of points for a given sample structure
%
% USAGE:
%              isFeasible = looplessCheckTOF(points, sample)
%
% INPUTS:
%              points:  n x k Matrix of points
%              sample (structure):  (the following fields are required - others can be supplied)
%                                    * loopMatrix - `ni x p` Loop law matrix
%                                    * nnzEntries - `1 x p` Matrix of non-zero entries of the loop law matrix
%                                    * intRxns - `ni x 1` Matrix with indexes of internal rxns
%
% OUTPUT:
%              isFeasible:  1 x k boolean vector

nSamples = size(points,2);
intRxns = sample.intRxns;
TICmat = sample.TICmat(intRxns,:);
nnzTIC = sparse(repmat(sum(abs(TICmat),1),nSamples,1));
flux = sparse(sign(points(intRxns,:)));

TICfluxMat = sum((flux'*TICmat)==nnzTIC,2);
isFeasible = TICfluxMat'==0;