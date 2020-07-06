function grp = AlternatingGroup(n)
    Sn = replab.S(n);
    if n <= 2
        grp = replab.PermutationGroup.trivial(n);
    else
        c3 = [2 3 1 4:n];
        if mod(n, 2) == 0
            s = [1 3:n 2];
        else
            s = [2:n 1];
        end
        grp = Sn.subgroup({c3 s}, Sn.order/2);
    end
end
