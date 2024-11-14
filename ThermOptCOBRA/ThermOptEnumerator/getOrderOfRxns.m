function [order,bins,binsizes] =getOrderOfRxns(model)
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