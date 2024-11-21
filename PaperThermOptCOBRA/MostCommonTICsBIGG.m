clear
%% to get the commonly occuring TICs in the BiGG models
unique_tics = {}; count = [];
p = dir('./BiggBlkdResults_tol_5/');
p = {p(3:end).name}';
p = strrep(p,'.mat','');

for i =1:numel(p)
    data = load(['./BiggBlkdResults_tol_5/',p{i}]);
    temp_tics={};
    for j = 1:numel(data.TICs)
        temp = strjoin(sort(data.TICs{j}));
        if ~ismember(temp,temp_tics)
            if ismember(temp,unique_tics)
                id = ismember(unique_tics,temp);
                count(id)=count(id)+1;
            else
                unique_tics = [unique_tics;temp];
                count = [count;1];
            end
            temp_tics = [temp_tics;temp];
        end
    end
end

% to sort the unique_tics based on the count
[count,i] = sort(count,'descend');
unique_tics = unique_tics(i);
save('Results_MostCommonTICsBIGG')