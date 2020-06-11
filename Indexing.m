function [inelements, table] = Indexing(index, table)
% Note: perm must be a permutation group element i.e. [1, 3, 2]

inelements = 0;
n = length(table);

% deal with first 2 elements
if n == 0
    table = [index];
    return
elseif n == 1
    if index < table{1}
        table = [index, table{1}];
    elseif index > table{1}
        table = [table{1}, index];
    else
        inelements = 1;
    end
    return
end
        

% Binary search existing hashtable
L = 1;
R = n;
while L <= R
    m = floor((L + R) / 2);
    if table(m) < index
        L = m + 1;
    elseif table(m) > index
        R = m - 1;
    else
        break
    end
end

% either hashtable(m) is hash or we add hash to the table
m = floor((L + R) / 2);
if table(m) == index
    inelements = 1;
else
    table = [table(1:m), index, table(m+1:n)];
end

end