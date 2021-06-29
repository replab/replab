function res = uniqueStable(values)
    [~,ind,~] = unique(values, 'first');
    res = values(sort(ind));
end
