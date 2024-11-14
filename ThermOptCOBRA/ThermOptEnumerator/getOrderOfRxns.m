function [order,bins,binsizes] =getOrderOfRxns(model)
% To get the order of reactions based on the connected components obtained
%
% USAGE: 
%   [order,bins,binsizes] =getOrderOfRxns(model)
%
% INPUTS:
%     model:     COBRA model structure for which the order of reaction has to be obtained
%
% OUTPUTS:
%     order:      Reaction IDs
%     bins:       Bins defining the IDs of connected components in which the reactions belong
%     binsizes:   Size of the connected components obtained
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras



S = double(model.S~=0); S=logical(S'*S);
g = graph(S);
try
    [bins,binsizes]=conncomp(g);
catch
    bins = conncomp(g);
    binsizes = 1:max(bins);
    binsizes = arrayfun(@(x)sum(bins==x),binsizes);
end
[binsizes,i]=sort(binsizes);
bins = arrayfun(@(x)find(i==x),bins);
order=[];
for i =1:numel(binsizes)
    order=[order;find(bins==i)'];
end
end
