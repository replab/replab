function grp = AlternatingGroup(n)
% Constructs the alternating group
%
% Args:
%   n (integer): Group degree
%
% Returns:
%   `+replab.PermutationGroup`: The alternating group of degree ``n``
    Sn = replab.S(n);
    if n <= 2 % special case: group is trivial
        grp = replab.PermutationGroup.trivial(n);
    else
        t = [2 3 1 4:n];
        if n == 3 % special case: it is a cyclic group, one generator only
            grp = Sn.subgroup({t}, Sn.order/2);
            return
        end
        % generators from page 2100 of
        % https://www.ams.org/journals/tran/2003-355-05/S0002-9947-03-03040-X/S0002-9947-03-03040-X.pdf
        if mod(n, 2) == 0
            s = [2 1 4:n 3];
        else
            s = [1 2 4:n 3];
        end
        grp = Sn.subgroup({s t}, Sn.order/2);
    end
end
