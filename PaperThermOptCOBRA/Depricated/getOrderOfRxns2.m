function [order,bins,binsizes] =getOrderOfRxns2(model,lb)
% this works based on the directed edges of the reaction graph. Results are
% similar for both type of edges. Hence we used the undirected one.
S = model.S;
rev=lb<0;
edgelist=[]; % first column is from edge and second column is to edge
for i=1:size(S,1)
    r = S(i,:);r=r';
    fr1 = find(r>0 & ~rev); to1 = find(r<0&~rev);
    fr2 = find(r~=0 & rev); to2 = find(r~=0 & rev);
    [t1,t2]=meshgrid([fr1;fr2],[to1;to2]);
    edges = [t1(:),t2(:)];
    edgelist=[edgelist;edges];
end

edgelist=unique(edgelist,'rows');
edgelist(find(edgelist(:,1)==edgelist(:,2)),:)=[];
g = digraph(edgelist(:,1),edgelist(:,2));
[bins,binsizes]=conncomp(g);
[binsizes,i]=sort(binsizes);
bins = arrayfun(@(x)find(i==x),bins);
order=[];
for i =1:numel(binsizes)
    order=[order;find(bins==i)'];
end
end