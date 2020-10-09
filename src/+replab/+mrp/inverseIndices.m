function ind = inverseIndices(group, elements)
% Returns a vector of integer indexing element inverses
%
% Args:
%   group (`+replab.Group`): Group to consider
%   elements (cell(1,\*) of elements of ``group``): List of group elements
%
% Returns:
%   integer(1,\*): Returns ``ind`` such that, for each ``i``, ``elements(ind(i))`` is the inverse of ``elements(i)`` if such an element exists; otherwise ``ind(i) = 0``.
    n = length(elements);
    ind = zeros(1, n);
    for i = 1:n
        if ind(i) == 0
            ginv = group.inverse(elements{i});
            for j = 1:n
                if group.eqv(ginv, elements{j})
                    ind(i) = j;
                    ind(j) = i;
                    break
                end
            end
        end
    end
end
