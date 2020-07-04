function t = sortByOrdering(s, ordering)
    [~, ind] = sort(ordering(s));
    t = s(ind);
end
